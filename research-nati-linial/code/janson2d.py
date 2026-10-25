#!/usr/bin/env python3
"""
Janson (1999) minimal-path experiment, lifted from K_n (dim=1) to the full
2-dimensional simplicial complex on [n] (dim=2).

MODEL (dim=2)
-------------
* Every triangle t = {a,b,c} of [n] carries an i.i.d. weight w(t) ~ Exp(1).
* The "points" between which distances are measured are the C(n,2) EDGES
  (codimension-1 faces) -- they play the role of Janson's vertices.
* A path from edge e to edge f is a chain of triangles
        e = e_0 ⊂ t_1 ⊃ e_1 ⊂ t_2 ⊃ e_2 ⊂ ... ⊂ t_k ⊃ e_k = f,
  i.e. consecutive triangles share an edge. Its weight is  w(t_1)+...+w(t_k).
* TRANSITION RULE used by Dijkstra (see `_sssp`):
        an edge {a,b} settled at distance d  --(triangle {a,b,c})-->
        edges {a,c} and {b,c}, both relaxed to  d + w({a,b,c}).
  A triangle's weight is paid once, on entering it; leaving through either of
  its two other edges is free. Equivalently: Dijkstra on the bipartite graph
  edges <-> triangles, arcs edge->triangle cost w(t), triangle->edge cost 0.

MODEL (dim=1, sanity mode) -- Janson's original: vertices, edges of K_n.

METRICS (Janson's three)
  pair : d(e0, f0) for two fixed disjoint edges  (1D: two fixed vertices)
  ecc  : max_f d(e0, f)                          (fixed-to-furthest)
  diam : max_{e,f} d(e, f)                       (global maximum)

EFFICIENCY
  * Truncation: only simplices with w < T, T = c * log n / n, are generated
    (O(M) geometric skipping, never an n^3 loop). A computed distance < T is
    provably exact: any path using a dropped simplex weighs >= T. If a result
    is >= T the complex is extended to a larger T *on the same realisation*
    (memorylessness: an un-generated simplex has weight T + Exp(1)), so the
    retry introduces no selection bias.
  * CSR incidence: edge -> list of (third vertex c, w({a,b,c})).
  * Numba-compiled Dijkstra with a binary heap (lazy deletion).
  * Diameter: SSSP from the edges with the largest *escape time*
    (lightest incident triangle) -- the extremal mechanism in Janson's proof --
    run in parallel. Validated against exact all-pairs for small n
    (`--validate`).

Usage
  python janson2d.py --validate          # correctness checks (small n)
  python janson2d.py --dim 2             # main 2D experiment
  python janson2d.py --dim 1             # Janson 1D sanity run
  python janson2d.py --dim 2 --quick     # smoke test
  python janson2d.py --plot-only --dim 2 # re-plot from saved CSV
"""

import argparse
import csv
import json
import logging
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import numba as nb
from numba import njit, prange

PROJECT_NAME = "Janson2D_Simplicial_Experiment"

# Predicted constants under the normalisation  d * n / log n
# (1D: Janson's theorem; 2D: branching-process heuristic, see README in summary)
PREDICTED = {1: {"pair": 1.0, "ecc": 2.0, "diam": 3.0},
             2: {"pair": 1.0, "ecc": 3.0, "diam": 5.0}}

# (n, Monte-Carlo reps, compute diameter?)
SCHEDULE = {
    2: [(50, 200, True), (100, 100, True), (200, 60, True), (500, 30, True),
        (1000, 15, True), (2000, 8, True), (3000, 5, True), (5000, 4, False)],
    1: [(100, 400, True), (300, 300, True), (1000, 200, True),
        (3000, 100, True), (10000, 50, True), (30000, 30, True),
        (100000, 15, True)],
}
QUICK = {2: [(30, 20, True), (60, 10, True), (120, 5, True)],
         1: [(100, 20, True), (1000, 10, True)]}

C_TRUNC = {"pair_ecc": 3.8, "diam": 7.0}   # T = c * log n / n
MEM_BUDGET = 5.5e9                          # bytes
K_CAND = 24                                 # diameter candidate sources

METRICS = ["pair", "ecc", "diam"]


# ----------------------------------------------------------------------------
# Project folder on the Desktop (Windows / WSL / macOS / Linux)
# ----------------------------------------------------------------------------
def _is_wsl():
    try:
        return "microsoft" in Path("/proc/version").read_text().lower()
    except OSError:
        return False


def get_desktop() -> Path:
    cands = []
    if _is_wsl():
        try:
            prof = subprocess.run(["cmd.exe", "/c", "echo %USERPROFILE%"],
                                  capture_output=True, text=True, timeout=10,
                                  cwd="/mnt/c").stdout.strip()
            if prof and "%" not in prof:
                lin = subprocess.run(["wslpath", "-u", prof], capture_output=True,
                                     text=True, timeout=10).stdout.strip()
                cands += [Path(lin) / "Desktop", Path(lin) / "OneDrive" / "Desktop"]
        except (OSError, subprocess.SubprocessError):
            pass
        cands += sorted(Path("/mnt/c/Users").glob("*/Desktop"))
    else:
        home = Path.home()
        if platform.system() == "Windows":
            cands.append(home / "OneDrive" / "Desktop")
        cands.append(home / "Desktop")
    for c in cands:
        if c.is_dir() and c.parent.name not in ("Public", "Default", "All Users"):
            return c
    return Path.home()


def setup_project(out=None) -> Path:
    proj = Path(out) if out else get_desktop() / PROJECT_NAME
    for sub in ("code", "data", "plots", "logs"):
        (proj / sub).mkdir(parents=True, exist_ok=True)
    me = Path(__file__).resolve()
    dst = proj / "code" / me.name
    if me != dst.resolve():
        shutil.copy2(me, dst)
    (proj / "code" / "requirements.txt").write_text(
        "numpy\nnumba\nscipy\nmatplotlib\nnetworkx\n")
    return proj


def setup_logging(proj: Path, tag: str):
    log = logging.getLogger("janson")
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s  %(message)s", "%H:%M:%S")
    for h in (logging.StreamHandler(sys.stdout),
              logging.FileHandler(proj / "logs" / f"run_{tag}.log", mode="a")):
        h.setFormatter(fmt)
        log.addHandler(h)
    return log


# ----------------------------------------------------------------------------
# Indexing
# ----------------------------------------------------------------------------
@njit(cache=True, inline="always")
def _eidx(n, x, y):
    """Index of the edge {x,y} in 0..C(n,2)-1 (lexicographic, x != y)."""
    if x > y:
        x, y = y, x
    return x * n - (x * (x + 1)) // 2 + (y - x - 1)


@njit(cache=True)
def _edge_endpoints(n):
    N = n * (n - 1) // 2
    ea = np.empty(N, np.int32)
    eb = np.empty(N, np.int32)
    k = 0
    for a in range(n):
        for b in range(a + 1, n):
            ea[k] = a
            eb[k] = b
            k += 1
    return ea, eb


# ----------------------------------------------------------------------------
# Random truncated complex, generated in layers [lo, hi) of weight
# ----------------------------------------------------------------------------
# Each simplex independently lands in layer [lo,hi) w.p. p = 1 - e^{-(hi-lo)};
# its weight is then lo + (Exp(1) conditioned < hi-lo). Layers [0,T0),[T0,T1),...
# reproduce exactly "all simplices with w < T_last" with correct joint law.
# A simplex can be drawn in two layers; the extra (heavier) copy never affects
# a shortest path or an escape time, so duplicates are harmless.

@njit(cache=True)
def _gen2(n, lo, hi, seeds, fill, pos, nbr, wt):
    total = 0
    for L in range(lo.shape[0]):
        np.random.seed(seeds[L])
        width = hi[L] - lo[L]
        p = -math.expm1(-width)
        for a in range(n):
            for b in range(a + 1, n):
                c = b
                while True:
                    f = -math.log(1.0 - np.random.random()) / width  # geometric skip
                    if f >= n:
                        break
                    c += 1 + int(f)
                    if c >= n:
                        break
                    w = lo[L] - math.log1p(-np.random.random() * p)
                    total += 1
                    e_ab = _eidx(n, a, b)
                    e_ac = _eidx(n, a, c)
                    e_bc = _eidx(n, b, c)
                    if fill:
                        k = pos[e_ab]; nbr[k] = c; wt[k] = w; pos[e_ab] = k + 1
                        k = pos[e_ac]; nbr[k] = b; wt[k] = w; pos[e_ac] = k + 1
                        k = pos[e_bc]; nbr[k] = a; wt[k] = w; pos[e_bc] = k + 1
                    else:
                        pos[e_ab] += 1
                        pos[e_ac] += 1
                        pos[e_bc] += 1
    return total


@njit(cache=True)
def _gen1(n, lo, hi, seeds, fill, pos, nbr, wt):
    total = 0
    for L in range(lo.shape[0]):
        np.random.seed(seeds[L])
        width = hi[L] - lo[L]
        p = -math.expm1(-width)
        for a in range(n):
            v = a
            while True:
                f = -math.log(1.0 - np.random.random()) / width
                if f >= n:
                    break
                v += 1 + int(f)
                if v >= n:
                    break
                w = lo[L] - math.log1p(-np.random.random() * p)
                total += 1
                if fill:
                    k = pos[a]; nbr[k] = v; wt[k] = w; pos[a] = k + 1
                    k = pos[v]; nbr[k] = a; wt[k] = w; pos[v] = k + 1
                else:
                    pos[a] += 1
                    pos[v] += 1
    return total


class Complex:
    """CSR incidence: for state e (an edge in 2D, a vertex in 1D),
    nbr[off[e]:off[e+1]] are the third vertices c of its triangles {a,b,c}
    (1D: the neighbours v), with weights wt[...]."""

    def __init__(self, dim, n, layers):
        self.dim, self.n = dim, n
        lo = np.array([l[0] for l in layers], np.float64)
        hi = np.array([l[1] for l in layers], np.float64)
        sd = np.array([l[2] for l in layers], np.int64)
        self.T = float(hi[-1])
        N = n if dim == 1 else n * (n - 1) // 2
        self.N = N
        ntype = np.uint16 if (dim == 2 and n < 65536) else np.int32
        gen = _gen2 if dim == 2 else _gen1
        deg = np.zeros(N, np.int64)
        gen(n, lo, hi, sd, False, deg, np.empty(0, ntype), np.empty(0, np.float32))
        off = np.zeros(N + 1, np.int64)
        np.cumsum(deg, out=off[1:])
        del deg
        tot = int(off[-1])
        nbr = np.empty(tot, ntype)
        wt = np.empty(tot, np.float32)
        pos = off[:-1].copy()
        self.M = gen(n, lo, hi, sd, True, pos, nbr, wt)
        del pos
        self.off, self.nbr, self.wt = off, nbr, wt
        if dim == 2:
            self.ea, self.eb = _edge_endpoints(n)
        else:
            self.ea = self.eb = np.empty(0, np.int32)

    def sssp(self, src, target=-1):
        return _sssp(self.dim, self.n, self.off, self.nbr, self.wt,
                     self.ea, self.eb, src, target)

    def eccentricities(self, sources):
        return _ecc_many(self.dim, self.n, self.off, self.nbr, self.wt,
                         self.ea, self.eb, np.asarray(sources, np.int64))

    def escape(self):
        return _escape(self.off, self.wt)


def estimate_bytes(dim, n, T, threads):
    N = n if dim == 1 else n * (n - 1) // 2
    k = 3 if dim == 2 else 2
    M = (math.comb(n, dim + 1)) * (-math.expm1(-T))
    nb_ = 2 if (dim == 2 and n < 65536) else 4
    csr = k * M * (nb_ + 4) + 2 * 8 * N + (8 * N if dim == 2 else 0)
    per_thread = N * (8 + 1) + 3 * N * 16          # dist, done, heap (generous)
    return 1.4 * (csr + threads * per_thread), M     # 1.4: measured RSS overhead


# ----------------------------------------------------------------------------
# Dijkstra
# ----------------------------------------------------------------------------
@njit(cache=True, inline="always")
def _push(hk, hv, size, key, val):
    if size == hk.shape[0]:
        nk = np.empty(2 * size, np.float64); nk[:size] = hk; hk = nk
        nv = np.empty(2 * size, np.int64); nv[:size] = hv; hv = nv
    i = size
    while i > 0:
        par = (i - 1) >> 1
        if hk[par] <= key:
            break
        hk[i] = hk[par]; hv[i] = hv[par]
        i = par
    hk[i] = key; hv[i] = val
    return hk, hv, size + 1


@njit(cache=True)
def _sssp(dim, n, off, nbr, wt, ea, eb, src, target):
    N = off.shape[0] - 1
    dist = np.full(N, np.inf)
    done = np.zeros(N, np.bool_)
    hk = np.empty(1024, np.float64)
    hv = np.empty(1024, np.int64)
    size = 0
    dist[src] = 0.0
    hk, hv, size = _push(hk, hv, size, 0.0, src)
    while size > 0:
        d = hk[0]; e = hv[0]
        size -= 1
        if size > 0:                       # pop: sift last element down
            key = hk[size]; val = hv[size]; i = 0
            while True:
                l = 2 * i + 1
                if l >= size:
                    break
                m = l
                if l + 1 < size and hk[l + 1] < hk[l]:
                    m = l + 1
                if hk[m] < key:
                    hk[i] = hk[m]; hv[i] = hv[m]; i = m
                else:
                    break
            hk[i] = key; hv[i] = val
        if done[e]:
            continue                        # stale heap entry
        done[e] = True
        if e == target:
            break
        if dim == 1:
            for j in range(off[e], off[e + 1]):
                f = nbr[j]
                nd = d + wt[j]
                if nd < dist[f]:
                    dist[f] = nd
                    hk, hv, size = _push(hk, hv, size, nd, f)
        else:
            a = ea[e]; b = eb[e]
            for j in range(off[e], off[e + 1]):
                # ---- TRANSITION: enter triangle {a,b,c} (pay its weight once),
                # ---- exit through edge {a,c} or edge {b,c} (free).
                c = np.int64(nbr[j])
                nd = d + wt[j]
                f = _eidx(n, a, c)
                if nd < dist[f]:
                    dist[f] = nd
                    hk, hv, size = _push(hk, hv, size, nd, f)
                f = _eidx(n, b, c)
                if nd < dist[f]:
                    dist[f] = nd
                    hk, hv, size = _push(hk, hv, size, nd, f)
    return dist


@njit(cache=True, parallel=True)
def _ecc_many(dim, n, off, nbr, wt, ea, eb, sources):
    S = sources.shape[0]
    ecc = np.empty(S, np.float64)
    arg = np.empty(S, np.int64)
    for i in prange(S):
        dist = _sssp(dim, n, off, nbr, wt, ea, eb, sources[i], -1)
        k = np.argmax(dist)
        ecc[i] = dist[k]
        arg[i] = k
    return ecc, arg


@njit(cache=True, parallel=True)
def _escape(off, wt):
    """Escape time of each state: weight of its lightest incident simplex."""
    N = off.shape[0] - 1
    out = np.full(N, np.inf)
    for e in prange(N):
        m = np.inf
        for j in range(off[e], off[e + 1]):
            if wt[j] < m:
                m = wt[j]
        out[e] = m
    return out


# ----------------------------------------------------------------------------
# One Monte-Carlo replicate
# ----------------------------------------------------------------------------
def _seed(base, dim, n, rep, layer):
    return (base * 1000003 + dim * 7919 + n * 104729 + rep * 1299709
            + layer * 15485863) % (2 ** 31 - 1)


def run_one(dim, n, rep, do_diam, base_seed=12345, max_layers=6):
    t0 = time.time()
    c = C_TRUNC["diam" if do_diam else "pair_ecc"]
    T = c * math.log(n) / n
    layers = [(0.0, T, _seed(base_seed, dim, n, rep, 0))]
    src, tgt = (0, 1) if dim == 1 else (0, _eidx(n, 2, 3))
    while True:
        cx = Complex(dim, n, layers)
        dist = cx.sssp(src)
        pair, ecc = float(dist[tgt]), float(dist.max())
        ok = ecc < cx.T
        diam, cand_rank, diam_ok = np.nan, -1, True
        if do_diam and ok:
            esc = cx.escape()
            k = min(K_CAND, cx.N)
            top = np.argpartition(-esc, k - 1)[:k]
            top = top[np.argsort(-esc[top])]           # slowest escape first
            eccs, _ = cx.eccentricities(top)
            diam = float(max(eccs.max(), ecc))
            cand_rank = int(np.argmax(eccs))
            diam_ok = diam < cx.T
        if ok and diam_ok:
            break
        if len(layers) >= max_layers:
            break
        lo = layers[-1][1]
        layers.append((lo, lo * 1.5, _seed(base_seed, dim, n, rep, len(layers))))
    return dict(dim=dim, n=n, rep=rep, T=cx.T, layers=len(layers), M=cx.M,
                pair=pair, ecc=ecc, diam=diam, exact=bool(ok and diam_ok),
                diam_cand_rank=cand_rank, seconds=round(time.time() - t0, 3))


# ----------------------------------------------------------------------------
# Experiment driver
# ----------------------------------------------------------------------------
FIELDS = ["dim", "n", "rep", "T", "layers", "M", "pair", "ecc", "diam",
          "exact", "diam_cand_rank", "seconds"]


def run_experiment(proj, dim, schedule, log, threads, base_seed, tag, on_n_done=None):
    path = proj / "data" / f"runs_{tag}.csv"
    done = set()
    if path.exists():
        with open(path) as fh:
            for r in csv.DictReader(fh):
                done.add((int(r["n"]), int(r["rep"])))
    new = not path.exists()
    fh = open(path, "a", newline="")
    wr = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        wr.writeheader()
    for n, reps, do_diam in schedule:
        c = C_TRUNC["diam" if do_diam else "pair_ecc"]
        th = threads                     # fewer parallel SSSPs if memory is tight
        while th > 1 and estimate_bytes(dim, n, c * math.log(n) / n, th)[0] > MEM_BUDGET:
            th -= 1
        nb.set_num_threads(th)
        need, M = estimate_bytes(dim, n, c * math.log(n) / n, th)
        if need > MEM_BUDGET and do_diam:
            do_diam = False
            c = C_TRUNC["pair_ecc"]
            need, M = estimate_bytes(dim, n, c * math.log(n) / n, 1)
            log.info(f"n={n}: diameter dropped (memory)")
        if need > MEM_BUDGET:
            log.info(f"n={n}: SKIPPED, est. {need/1e9:.1f} GB > budget")
            continue
        log.info(f"--- dim={dim} n={n} reps={reps} diam={do_diam} "
                 f"est. M={M:.3g} simplices, {need/1e9:.2f} GB, {th} threads")
        for rep in range(reps):
            if (n, rep) in done:
                continue
            r = run_one(dim, n, rep, do_diam, base_seed)
            wr.writerow(r)
            fh.flush()
            s = n / math.log(n)
            log.info(f"n={n:6d} rep={rep:3d}  pair*n/ln n={r['pair']*s:6.3f}  "
                     f"ecc*n/ln n={r['ecc']*s:6.3f}  diam*n/ln n={r['diam']*s:6.3f}"
                     f"  M={r['M']:.3g} layers={r['layers']} exact={r['exact']}"
                     f" cand={r['diam_cand_rank']}  {r['seconds']:.1f}s")
        if on_n_done:
            on_n_done()                  # refresh graphs + results file
    fh.close()
    return path


def load_runs(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            rows.append({k: (float(v) if k not in ("exact",) else v == "True")
                         for k, v in r.items()})
    return rows


def summarize(rows, dim, out_csv):
    ns = sorted({int(r["n"]) for r in rows})
    summ = []
    for n in ns:
        rs = [r for r in rows if int(r["n"]) == n]
        s = n / math.log(n)
        row = dict(n=n, reps=len(rs), exact_frac=np.mean([r["exact"] for r in rs]),
                   mean_M=np.mean([r["M"] for r in rs]),
                   mean_seconds=np.mean([r["seconds"] for r in rs]))
        for m in METRICS:
            v = np.array([r[m] for r in rs if r["exact"] and np.isfinite(r[m])])
            if len(v) == 0:
                for k in ("mean", "sd", "ci95", "norm", "norm_ci95"):
                    row[f"{m}_{k}"] = np.nan
                continue
            sd = v.std(ddof=1) if len(v) > 1 else np.nan
            ci = 1.96 * sd / math.sqrt(len(v)) if len(v) > 1 else np.nan
            row.update({f"{m}_mean": v.mean(), f"{m}_sd": sd, f"{m}_ci95": ci,
                        f"{m}_norm": v.mean() * s, f"{m}_norm_ci95": ci * s})
        summ.append(row)
    with open(out_csv, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(summ[0].keys()))
        wr.writeheader()
        wr.writerows(summ)
    return summ


def extrapolate(summ, m):
    """Weighted LSQ fit  norm(n) = a + b / log n ; returns (a, b, a_se)."""
    pts = [(1 / math.log(r["n"]), r[f"{m}_norm"], r[f"{m}_norm_ci95"])
           for r in summ if np.isfinite(r[f"{m}_norm"])]
    if len(pts) < 3:
        return np.nan, np.nan, np.nan
    x, y, ci = map(np.array, zip(*pts))
    w = 1 / np.maximum(np.nan_to_num(ci / 1.96, nan=np.nanmax(ci)), 1e-9) ** 2
    A = np.vstack([np.ones_like(x), x]).T
    W = np.diag(w)
    cov = np.linalg.inv(A.T @ W @ A)
    a, b = cov @ A.T @ W @ y
    return float(a), float(b), float(math.sqrt(cov[0, 0]))


def loglog_slope(summ, m):
    pts = [(math.log(r["n"]), math.log(r[f"{m}_mean"])) for r in summ
           if np.isfinite(r[f"{m}_mean"])]
    if len(pts) < 2:
        return np.nan
    x, y = map(np.array, zip(*pts))
    return float(np.polyfit(x, y, 1)[0])


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------
COL = {"pair": "#2a78d6", "ecc": "#eb6834", "diam": "#1baf7a"}
LABEL = {"pair": "fixed pair", "ecc": "fixed → furthest", "diam": "diameter"}
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def _style(ax):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def make_plots(proj, dim, rows, summ, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK})
    pdir = proj / "plots"
    unit = "edges (triangle paths)" if dim == 2 else "vertices (Janson 1D)"
    ns = np.array([r["n"] for r in summ])

    # 1) normalised curves
    fig, ax = plt.subplots(figsize=(8, 5.2), facecolor=SURF)
    _style(ax)
    for m in METRICS:
        y = np.array([r[f"{m}_norm"] for r in summ])
        e = np.array([r[f"{m}_norm_ci95"] for r in summ])
        ok = np.isfinite(y)
        if not ok.any():
            continue
        ax.errorbar(ns[ok], y[ok], yerr=np.nan_to_num(e[ok]), color=COL[m], lw=2,
                    marker="o", ms=6, capsize=3, label=LABEL[m],
                    mec=SURF, mew=1.5)
        c = PREDICTED[dim][m]
        ax.axhline(c, color=COL[m], ls="--", lw=1, alpha=0.6)
        ax.annotate(f"{LABEL[m]}  (pred. {c:g})", (ns[ok][-1], y[ok][-1]),
                    xytext=(8, 0), textcoords="offset points", va="center",
                    fontsize=9, color=INK2)
    if dim == 2:
        for c in (2.0, 3.0):
            ax.axhline(c, color=INK2, ls=":", lw=0.8)
        ax.text(ns[0], 2.0, " Janson 1D: 1, 2, 3 (dotted)", fontsize=8,
                color=INK2, va="bottom")
    ax.set_xscale("log")
    ax.set_xlabel("n (number of vertices)")
    ax.set_ylabel(r"normalised distance   $d \cdot n / \log n$")
    ax.set_title(f"{dim}D: minimal path weights between {unit}", loc="left",
                 fontsize=11, color=INK)
    ax.legend(frameon=False, fontsize=9, loc="upper center",
              bbox_to_anchor=(0.5, -0.13), ncol=3)
    ax.set_xlim(ns.min() / 1.3, ns.max() * 3.5)
    fig.tight_layout()
    fig.savefig(pdir / f"{tag}_normalised.png", dpi=160)
    plt.close(fig)

    # 2) raw log-log with slope
    fig, ax = plt.subplots(figsize=(7, 5), facecolor=SURF)
    _style(ax)
    for m in METRICS:
        y = np.array([r[f"{m}_mean"] for r in summ])
        ok = np.isfinite(y)
        if not ok.any():
            continue
        sl = loglog_slope(summ, m)
        ax.plot(ns[ok], y[ok], color=COL[m], lw=2, marker="o", ms=6, mec=SURF,
                mew=1.5, label=f"{LABEL[m]}  (slope {sl:.3f})")
    ref = np.log(ns) / ns
    ax.plot(ns, ref, color=INK2, ls="--", lw=1, label=r"$\log n / n$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("n"); ax.set_ylabel("mean minimal path weight")
    ax.set_title(f"{dim}D: raw distances (log–log)", loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(pdir / f"{tag}_loglog.png", dpi=160)
    plt.close(fig)

    # 3) extrapolation in 1/log n
    fig, ax = plt.subplots(figsize=(7, 5), facecolor=SURF)
    _style(ax)
    xs = np.linspace(0, 1 / math.log(ns.min()), 50)
    for m in METRICS:
        y = np.array([r[f"{m}_norm"] for r in summ])
        e = np.array([r[f"{m}_norm_ci95"] for r in summ])
        ok = np.isfinite(y)
        if not ok.any():
            continue
        a, b, ase = extrapolate(summ, m)
        ax.errorbar(1 / np.log(ns[ok]), y[ok], yerr=np.nan_to_num(e[ok]),
                    color=COL[m], ls="none", marker="o", ms=6, capsize=3,
                    mec=SURF, mew=1.5,
                    label=f"{LABEL[m]}: a = {a:.2f} ± {1.96*ase:.2f}")
        if np.isfinite(a):
            ax.plot(xs, a + b * xs, color=COL[m], lw=1.2, alpha=0.7)
        ax.plot(0, PREDICTED[dim][m], marker="*", ms=12, color=COL[m], mec=INK)
    ax.set_xlim(left=0)
    ax.set_xlabel(r"$1/\log n$   (n → ∞ at the left edge)")
    ax.set_ylabel(r"$d \cdot n / \log n$")
    ax.set_title(f"{dim}D: extrapolation  y = a + b/log n   (★ = predicted)",
                 loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(pdir / f"{tag}_extrapolation.png", dpi=160)
    plt.close(fig)

    # 4) fluctuations at the largest n having each metric
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.8), facecolor=SURF)
    for ax, m in zip(axs, METRICS):
        _style(ax)
        good = [r for r in rows if r["exact"] and np.isfinite(r[m])]
        if not good:
            continue
        nmax = max(r["n"] for r in good)
        # pool the two largest n with many reps for a usable histogram
        nsel = sorted({r["n"] for r in good})
        nsel = [x for x in nsel if sum(1 for r in good if r["n"] == x) >= 30][-1:] or [nmax]
        v = np.array([r[m] * r["n"] - PREDICTED[dim][m] * math.log(r["n"])
                      for r in good if r["n"] in nsel])
        ax.hist(v, bins=25, color=COL[m], edgecolor=SURF, linewidth=1)
        ax.set_title(f"{LABEL[m]}  (n={int(nsel[0])}, {len(v)} reps)",
                     loc="left", fontsize=10)
        ax.set_xlabel(rf"$n\,d - {PREDICTED[dim][m]:g}\log n$")
    fig.tight_layout()
    fig.savefig(pdir / f"{tag}_fluctuations.png", dpi=160)
    plt.close(fig)

    # 5) second-order check: n*d - C*log n should level off if C is right
    fig, ax = plt.subplots(figsize=(7, 5), facecolor=SURF)
    _style(ax)
    for m in METRICS:
        pts = []
        for n in ns:
            v = np.array([r[m] * r["n"] - PREDICTED[dim][m] * math.log(r["n"])
                          for r in rows if r["n"] == n and r["exact"] and np.isfinite(r[m])])
            if len(v):
                ci = 1.96 * v.std(ddof=1) / math.sqrt(len(v)) if len(v) > 1 else 0.0
                pts.append((n, v.mean(), ci))
        if not pts:
            continue
        x, y, e = map(np.array, zip(*pts))
        ax.errorbar(x, y, yerr=e, color=COL[m], lw=2, marker="o", ms=6, capsize=3,
                    mec=SURF, mew=1.5,
                    label=rf"{LABEL[m]}:  $n d - {PREDICTED[dim][m]:g}\log n$")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("n"); ax.set_ylabel("mean of  n·d − C·log n")
    ax.set_title(f"{dim}D: second-order check (flat ⇒ predicted C is right)",
                 loc="left", fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(pdir / f"{tag}_centered.png", dpi=160)
    plt.close(fig)

    # 6) dashboard: everything on one page, with the numbers printed
    fig = plt.figure(figsize=(14, 9.5), facecolor=SURF)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.35, wspace=0.22)
    ax = fig.add_subplot(gs[0, 0]); _style(ax)
    for m in METRICS:
        y = np.array([r[f"{m}_norm"] for r in summ])
        e = np.array([r[f"{m}_norm_ci95"] for r in summ])
        ok = np.isfinite(y)
        if not ok.any():
            continue
        ax.errorbar(ns[ok], y[ok], yerr=np.nan_to_num(e[ok]), color=COL[m], lw=2,
                    marker="o", ms=5, capsize=3, mec=SURF, mew=1.2, label=LABEL[m])
        ax.axhline(PREDICTED[dim][m], color=COL[m], ls="--", lw=1, alpha=0.6)
        for xv, yv in zip(ns[ok], y[ok]):
            ax.annotate(f"{yv:.2f}", (xv, yv), xytext=(0, 7), textcoords="offset points",
                        ha="center", fontsize=7, color=COL[m])
    ax.set_xscale("log")
    ax.set_xlabel("n"); ax.set_ylabel(r"$d \cdot n / \log n$")
    ax.set_title("normalised distances (dashed = predicted)", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="lower right", bbox_to_anchor=(1, 0.2))

    ax = fig.add_subplot(gs[0, 1]); _style(ax)
    for m in METRICS:
        y = np.array([r[f"{m}_norm"] for r in summ])
        e = np.array([r[f"{m}_norm_ci95"] for r in summ])
        ok = np.isfinite(y)
        if not ok.any():
            continue
        a, b, ase = extrapolate(summ, m)
        ax.errorbar(1 / np.log(ns[ok]), y[ok], yerr=np.nan_to_num(e[ok]), color=COL[m],
                    ls="none", marker="o", ms=5, capsize=3, mec=SURF, mew=1.2,
                    label=f"{LABEL[m]}: limit {a:.2f} ± {1.96*ase:.2f}")
        if np.isfinite(a):
            ax.plot(xs, a + b * xs, color=COL[m], lw=1.2, alpha=0.7)
        ax.plot(0, PREDICTED[dim][m], marker="*", ms=12, color=COL[m], mec=INK)
    ax.set_xlim(left=0)
    ax.set_xlabel(r"$1/\log n$   (n = ∞ at 0)"); ax.set_ylabel(r"$d \cdot n / \log n$")
    ax.set_title("extrapolation to n → ∞  (★ = predicted)", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)

    ax = fig.add_subplot(gs[1, :]); ax.axis("off")
    head = ["n", "reps"] + [f"{LABEL[m]}" for m in METRICS]
    cells = [[f"{int(r['n'])}", f"{int(r['reps'])}"] +
             [("—" if not np.isfinite(r[f"{m}_norm"]) else
               f"{r[f'{m}_norm']:.3f} ± {np.nan_to_num(r[f'{m}_norm_ci95']):.3f}")
              for m in METRICS] for r in summ]
    fits = ["n → ∞ (fit)", ""]
    for m in METRICS:
        a, _, ase = extrapolate(summ, m)
        fits.append("—" if not np.isfinite(a) else f"{a:.3f} ± {1.96*ase:.3f}")
    cells += [fits, ["predicted", ""] + [f"{PREDICTED[dim][m]:g}" for m in METRICS]]
    tb = ax.table(cellText=cells, colLabels=head, loc="upper center", cellLoc="center")
    tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 1.35)
    for (i, j), cell in tb.get_celld().items():
        cell.set_edgecolor(GRID)
        if i == 0 or i >= len(cells) - 1:
            cell.set_text_props(weight="bold")
    fig.suptitle(f"{dim}D minimal paths between {unit}   —   updated "
                 f"{time.strftime('%Y-%m-%d %H:%M')}", x=0.06, ha="left", fontsize=12)
    fig.savefig(pdir / f"{tag}_dashboard.png", dpi=130, bbox_inches="tight")
    plt.close(fig)


def write_results(proj, dim, summ, tag):
    """Human-readable results file (+ JSON twin) in the project root."""
    fits = {m: dict(zip(("limit", "b", "limit_se"), extrapolate(summ, m)))
            for m in METRICS}
    lines = [f"# Results: {dim}D Janson experiment",
             f"Updated {time.strftime('%Y-%m-%d %H:%M')}. Values are d·n/log n, "
             "mean ± 95% CI. Raw runs: data/runs_" + tag + ".csv", "",
             "| n | reps | " + " | ".join(LABEL[m] for m in METRICS) + " |",
             "|---|---|" + "---|" * len(METRICS)]
    for r in summ:
        lines.append(f"| {int(r['n'])} | {int(r['reps'])} | " + " | ".join(
            "—" if not np.isfinite(r[f"{m}_norm"]) else
            f"{r[f'{m}_norm']:.3f} ± {np.nan_to_num(r[f'{m}_norm_ci95']):.3f}"
            for m in METRICS) + " |")
    lines += ["", "## Extrapolation  d·n/log n = a + b/log n", "",
              "| metric | a (n → ∞) | b | predicted | log-log slope |",
              "|---|---|---|---|---|"]
    for m in METRICS:
        f = fits[m]
        lines.append(f"| {LABEL[m]} | {f['limit']:.3f} ± {1.96*f['limit_se']:.3f} | "
                     f"{f['b']:.2f} | {PREDICTED[dim][m]:g} | "
                     f"{loglog_slope(summ, m):.3f} |")
    lines += ["", f"Plots: plots/{tag}_dashboard.png (all in one), "
              f"plots/{tag}_normalised.png, plots/{tag}_extrapolation.png, "
              f"plots/{tag}_loglog.png, plots/{tag}_fluctuations.png, "
              f"plots/{tag}_centered.png", ""]
    (proj / f"results_{tag}.md").write_text("\n".join(lines), encoding="utf-8")
    clean = lambda v: None if isinstance(v, float) and not math.isfinite(v) else v
    (proj / f"results_{tag}.json").write_text(json.dumps(
        dict(dim=dim, updated=time.strftime("%Y-%m-%d %H:%M"), predicted=PREDICTED[dim],
             by_n=[{k: clean(float(v)) for k, v in r.items()} for r in summ],
             fits={m: {k: clean(v) for k, v in f.items()} for m, f in fits.items()}),
        indent=2))


def report(proj, dim, tag, log=None):
    """Re-summarise everything saved so far, redraw plots, rewrite results."""
    path = proj / "data" / f"runs_{tag}.csv"
    if not path.exists():
        return None
    rows = load_runs(path)
    if not rows:
        return None
    summ = summarize(rows, dim, proj / "data" / f"summary_{tag}.csv")
    make_plots(proj, dim, rows, summ, tag)
    write_results(proj, dim, summ, tag)
    if log:
        log.info(f"[report] updated results_{tag}.md and plots/{tag}_dashboard.png")
    return summ


# ----------------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------------
def validate(proj, log):
    import networkx as nx
    report = {}

    # (a) Dijkstra vs networkx on the explicit bipartite edge/triangle graph,
    #     with NO truncation (T = 60: every simplex present).
    errs = []
    for dim in (1, 2):
        for n in (6, 8, 10, 12):
            cx = Complex(dim, n, [(0.0, 60.0, 1000 + n)])
            exp_M = math.comb(n, dim + 1)
            assert cx.M == exp_M, (dim, n, cx.M, exp_M)
            G = nx.DiGraph()
            if dim == 2:
                for e in range(cx.N):
                    a, b = int(cx.ea[e]), int(cx.eb[e])
                    for j in range(cx.off[e], cx.off[e + 1]):
                        c = int(cx.nbr[j])
                        t = ("t",) + tuple(sorted((a, b, c)))
                        G.add_edge(("e", e), t, weight=float(cx.wt[j]))
                        G.add_edge(t, ("e", e), weight=0.0)
            else:
                for u in range(n):
                    for j in range(cx.off[u], cx.off[u + 1]):
                        G.add_edge(("e", u), ("e", int(cx.nbr[j])),
                                   weight=float(cx.wt[j]))
            for src in range(min(cx.N, 5)):
                ref = nx.single_source_dijkstra_path_length(G, ("e", src))
                mine = cx.sssp(src)
                for e in range(cx.N):
                    errs.append(abs(ref[("e", e)] - mine[e]))
    report["dijkstra_vs_networkx_max_abs_err"] = max(errs)
    log.info(f"[validate] Dijkstra vs networkx: max |err| = {max(errs):.2e}")

    # (b) triangle count matches expectation (generator law)
    n, T = 200, 0.05
    Ms = [Complex(2, n, [(0.0, T, s)]).M for s in range(20)]
    exp = math.comb(n, 3) * (-math.expm1(-T))
    report["gen_M_mean_vs_expected"] = [float(np.mean(Ms)), exp]
    log.info(f"[validate] triangle count mean {np.mean(Ms):.0f} vs expected {exp:.0f}")

    # (c) weight law: weights in a layer should be Exp(1) | <T
    cx = Complex(2, 300, [(0.0, 0.03, 7)])
    w = cx.wt.astype(np.float64)
    report["gen_weight_mean_vs_expected"] = [
        float(w.mean()), 1 - 0.03 * math.exp(-0.03) / (-math.expm1(-0.03))]
    log.info(f"[validate] truncated weight mean {w.mean():.5f} vs "
             f"{report['gen_weight_mean_vs_expected'][1]:.5f}")

    # (d) diameter: candidate method vs exact all-pairs
    res = {}
    for dim, nlist in ((2, (15, 20, 30, 40, 60, 80)), (1, (50, 100, 200, 400))):
        for n in nlist:
            hits, gaps = 0, []
            for rep in range(20):
                cx = Complex(dim, n, [(0.0, 50.0, 5000 + 31 * n + rep)])
                eccs, _ = cx.eccentricities(np.arange(cx.N))
                exact = eccs.max()
                esc = cx.escape()
                k = min(K_CAND, cx.N)
                top = np.argsort(-esc)[:k]
                ecand, _ = cx.eccentricities(top)
                cand = max(ecand.max(), cx.sssp(0 if dim == 1 else 0).max())
                hits += abs(cand - exact) < 1e-12
                gaps.append((exact - cand) / exact)
            res[f"dim{dim}_n{n}"] = dict(match=hits / 20, max_rel_gap=max(gaps))
            log.info(f"[validate] diameter dim={dim} n={n}: candidate == exact "
                     f"in {hits}/20, max rel. gap {max(gaps):.2e}")
    report["diameter_candidate_vs_exact"] = res
    (proj / "logs" / "validation.json").write_text(json.dumps(report, indent=2))
    return report


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dim", type=int, default=2, choices=(1, 2))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--ns", type=str, default=None,
                    help="comma list n:reps[:nodiam], overrides schedule")
    ap.add_argument("--threads", type=int, default=nb.config.NUMBA_NUM_THREADS)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", type=str, default=None, help="project folder override")
    args = ap.parse_args()

    nb.set_num_threads(args.threads)
    proj = setup_project(args.out)
    log = setup_logging(proj, "validate" if args.validate else
                        f"dim{args.dim}" + ("_quick" if args.quick else ""))
    log.info(f"project folder: {proj}")

    if args.validate:
        validate(proj, log)
        return

    tag = f"dim{args.dim}" + ("_quick" if args.quick else "")
    if not args.plot_only:
        if args.ns:
            sched = []
            for tok in args.ns.split(","):
                p = tok.split(":")
                sched.append((int(p[0]), int(p[1]), not (len(p) > 2 and p[2] == "nodiam")))
        else:
            sched = (QUICK if args.quick else SCHEDULE)[args.dim]
        run_experiment(proj, args.dim, sched, log, args.threads, args.seed, tag,
                       on_n_done=lambda: report(proj, args.dim, tag, log))

    summ = report(proj, args.dim, tag, log)
    log.info("summary (normalised by log n / n):")
    for r in summ:
        log.info("  n={:6d} reps={:3d}  ".format(int(r["n"]), int(r["reps"])) +
                 "  ".join(f"{m}={r[f'{m}_norm']:.3f}±{r[f'{m}_norm_ci95']:.3f}"
                           for m in METRICS))
    for m in METRICS:
        a, b, ase = extrapolate(summ, m)
        log.info(f"  {m:5s}: extrapolated a={a:.3f}±{1.96*ase:.3f} (b={b:.2f}),"
                 f" predicted {PREDICTED[args.dim][m]:g}, log-log slope "
                 f"{loglog_slope(summ, m):.3f}")


if __name__ == "__main__":
    main()
