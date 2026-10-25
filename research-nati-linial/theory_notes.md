# Minimal triangle-paths in the random 2-complex: what to prove, and how

Setting. Every triangle $t\subset[n]$ gets an independent weight $w(t)\sim\mathrm{Exp}(1)$.
The points are the $N=\binom n2$ edges. A path $e=e_0\subset t_1\supset e_1\subset\dots\supset e_k=f$
costs $\sum w(t_i)$. Three quantities:

| quantity | definition | conjectured limit of $d\cdot n/\log n$ |
|---|---|---|
| pair | $d(e,f)$, $e,f$ fixed disjoint edges | **1** |
| ecc (fixed → furthest) | $\max_f d(e,f)$ | **3** |
| diam | $\max_{e,f} d(e,f)$ | **5** |

Janson (1999) proved the 1D analogue (vertices of $K_n$, Exp(1) edge weights): **1, 2, 3**.

---

## 0. A useful reformulation

From an edge $e=\{a,b\}$ the next edge can be any edge sharing exactly one vertex with $e$
(there are $2(n-2)$ of them), and the step $\{a,b\}\to\{a,c\}$ goes through the unique triangle $\{a,b,c\}$.
So **the model is first-passage percolation on the Johnson graph $J(n,2)$** (the line graph of $K_n$),
where the three $J(n,2)$-edges inside a "$K_n$-triangle" $\{ab,ac,bc\}$ share one weight.
$J(n,2)$ is $2(n-2)$-regular and its spectrum is explicit: $2n-4\ (\times1),\ n-4\ (\times n-1),\ -2\ (\times n(n-3)/2)$.
That is what makes the non-asymptotic bound in §3 work.

---

## 1. Why 1, 3, 5: the heuristic

Explore from $e$ (Dijkstra/first-passage). A triangle becomes *touched* when one of its edges is reached,
and by memorylessness every touched, not-yet-used triangle "rings" at rate 1; when it rings, its
unreached edges are reached.

1. **Growth.** Each reached edge has $n-2$ triangles, and each ring adds 2 new edges, so early on
   $\#\text{reached}(t)\approx W e^{2nt}$ (a Yule process with jumps of $+2$; $W\sim\Gamma(\tfrac12,2)$).
   The growth rate is $2n$ (in 1D it is $n$).
2. **Pair = 1.** Run the exploration from both ends. The two clouds meet when
   $|A|\cdot|B|\cdot\Pr[\text{two edges share a vertex}]\approx 1$, i.e. $e^{4nt}\cdot\frac4n\cdot\frac1n\approx1$,
   so each side needs $t\approx\frac{\log n}{2n}$ and the total is $\frac{\log n}{n}$.
   Equivalently, the time to reach a typical edge is $\frac{\log N}{2n}=\frac{2\log n}{2n}$.
   The doubling of the exponent ($N\approx n^2/2$) is exactly cancelled by the doubled growth rate.
3. **Escape time.** An edge leaves (or is entered) only through its $n-2$ triangles, so its escape time is
   $\min$ of $n-2$ Exp(1) values, i.e. $\mathrm{Exp}(n-2)$. The maximum over all $N\approx n^2/2$ edges is
   $\approx\frac{\log N}{n}\approx\frac{2\log n}{n}$, **twice** the 1D value ($\frac{\log n}{n}$ over $n$ vertices).
4. **ecc = 1 + 2 = 3.** The furthest edge from $e$ is the slowest edge. You pay its extra escape time $2\frac{\log n}{n}$
   on top of a typical distance $\frac{\log n}{n}$.
5. **diam = 2 + 1 + 2 = 5.** Both ends are slow edges.

The structure is identical to Janson's; only two numbers change. **The growth rate doubles (n → 2n), and the
number of points squares (n → n²/2).** The first change does not affect pair; the second doubles each escape term.

**General dimension (prediction).** In the $k$-dimensional complex the points are the $(k-1)$-faces, $N\approx n^k/k!$,
the growth rate is $kn$, and the escape rate is $n$. That gives
$$\text{pair}\to1,\qquad \text{ecc}\to1+k,\qquad \text{diam}\to1+2k .$$
$k=1$ recovers Janson (1, 2, 3), and $k=2$ gives (1, 3, 5). A 3D run (predicted 1, 4, 7) would be a strong independent test.

---

## 2. Proof plan for the limit theorem

**Target theorem.** As $n\to\infty$, in probability:
$$\frac{n\,d(e,f)}{\log n}\to1,\qquad\frac{n\max_f d(e,f)}{\log n}\to3,\qquad\frac{n\max_{e,f}d(e,f)}{\log n}\to5.$$

Each statement needs a lower bound and an upper bound. Janson's 1D proof is the template.

### Lower bounds (first moment; easiest)
* **pair ≥ 1−ε.** Count triangle-paths of length $k$ from $e$ to $f$ (walks in $J(n,2)$; see §3).
  A minimal path uses distinct triangles, so its weight is $\Gamma(k,1)$ and $\Pr[\le t]\le t^k/k!$.
  The expected number of paths with weight $\le t$ is about $e^{2nt}/N$. This is $\to0$ at $t=(1-\varepsilon)\frac{\log n}{n}$.
* **ecc ≥ 3−ε.** (a) Second moment: whp some edge $f$ has *all* $n-2$ triangles heavier than
  $s=(2-\varepsilon)\frac{\log n}n$. The expected count is $N e^{-(n-2)s}=n^{\varepsilon}/2\to\infty$. Two edges share at most
  one triangle, so the variance is small. (b) Condition on this event. Each of $f$'s triangles then has weight $s+\mathrm{Exp}(1)$
  (memorylessness). So $d(e,f)\ge s+(\text{time for } e\text{'s cloud to reach a constant fraction of } f\text{'s neighbourhood})$,
  which is $\ge(1-\varepsilon)\frac{\log n}n$ by the same first-moment count.
* **diam ≥ 5−ε.** Find two far-apart slow edges with the same second-moment argument, then repeat (b) at both ends.

### Upper bounds (exploration + concentration; the real work)
* **Growth phase.** Couple the exploration from $e$ with the Yule process with jumps $+2$ at rate $n-2$ per edge.
  This holds up to size $n^{2-\delta}$, because collisions (an exit edge already reached) are rare before that point.
  Concentration: after the first $O(1/n)$ time, the size is $W e^{2nt}(1+o(1))$.
* **Connection.** Run explorations from $e$ and from $f$ until each has about $n\sqrt{\log n}$ edges.
  Then there are about $|A||B|\cdot\frac4n\gg n\log n$ triangles with one edge in each cloud, each with a fresh Exp(1) weight (they are
  touched but not yet rung). The lightest has weight $o(\log n/n)$. Total: $(1+\varepsilon)\frac{\log n}n$ whp. **(pair upper bound)**
* **ecc/diam upper bounds.** Every edge escapes within $(2+\varepsilon)\frac{\log n}{n}$, by a union bound over $N$ edges:
  $N e^{-(n-2)(2+\varepsilon)\log n/n}\to0$. After escaping, an edge has 3 reached edges, and from a "typical start" the
  growth phase is fast. Glue: $d(e,f)\le \text{esc}(e)+\text{typical}+\text{esc}(f)$. This needs a uniform-in-$e$ version of
  the growth-phase concentration, which Janson does with a union bound over starts plus exponential tail estimates.

**Main technical differences from 1D.** (i) A triangle's weight is shared by its 3 edges, so after a ring two new
edges appear *together* (jumps of $+2$, correlated). (ii) The underlying graph is $J(n,2)$, not $K_n$. Two edges have
$n-2$ common $J$-neighbours if they share a vertex and 4 if they are disjoint, so the "collision" estimates have to
use this structure. Neither difference changes the leading constants; both change the $O(1)$ corrections.

---

## 3. A bound that holds for every $n$ (not just asymptotically)

Let $P_k$ be the number of walks of length $k$ in $J(n,2)$ from $e$ to $f$ (disjoint). The spectral decomposition of
$J(n,2)$ gives the **exact** count (checked numerically by matrix powers):
$$P_k=\frac{(2n-4)^k}{N}-\frac{4\,(n-4)^k}{n(n-2)}+\frac{2\,(-2)^k}{(n-1)(n-2)} .$$
**Proof of the walk count (proved).** Let $A$ be the adjacency matrix of $J(n,2)$, $J$ the all-ones matrix, and $B$ the
"disjoint" matrix ($B_{ef}=1$ iff $e\cap f=\emptyset$), so $I+A+B=J$. Two edges sharing a vertex have $n-2$ common
$J(n,2)$-neighbours; two disjoint edges have exactly $4$. Hence $A^2=(2n-4)I+(n-2)A+4B$, and substituting $B=J-I-A$,
$$(A-(n-4)I)(A+2I)=4J .$$
On $\mathbf 1^\perp$ the right side vanishes, so the only eigenvalues there are $n-4$ and $-2$; on $\mathbf 1$ the eigenvalue is $2n-4$.
The spectral projections are $E_0=J/N$, $E_1=\big(A+2I-\tfrac{2n-2}{N}J\big)/(n-2)$, $E_2=I-E_0-E_1$, and $A^k=(2n-4)^kE_0+(n-4)^kE_1+(-2)^kE_2$.
For disjoint $e,f$ we have $I_{ef}=A_{ef}=0$, so $(E_0)_{ef}=1/N$, $(E_1)_{ef}=-4/(n(n-2))$, $(E_2)_{ef}=2/((n-1)(n-2))$, which is the formula.
Check: $n=4$ ($J(4,2)$ is the octahedron) gives $P_2=4$, $P_3=8$; matrix powers for $n=4,5,6$, $k\le4$ agree.

Take a minimal path; it uses distinct triangles, so its weight is $\Gamma(k,1)$ and $\Pr[\le t]\le t^k/k!$.
A union bound over walks, after dropping the negative middle term, gives **for every $n\ge4$ and every $t>0$**:
$$\Pr\big[d(e,f)\le t\big]\;\le\;\frac{e^{2(n-2)t}-1}{N}+\frac{2\,(e^{2t}-1)}{(n-1)(n-2)} .$$
In centered form, with $t=(\log n-y)/n$:
$$\Pr\big[n\,d(e,f)\le\log n-y\big]\;\le\;\tfrac{2n}{n-1}\,e^{-2y}+O(\log n/n^3).$$
So "pair is at least $(1-\varepsilon)\log n/n$" fails with probability at most about $2n^{-2\varepsilon}$, **explicitly at each finite $n$**.
The simulated lower tails respect this at every $n$ with enough reps.

"For all $n$" cannot mean that $d\,n/\log n$ *equals* 1, 3, 5. At finite $n$ there is a deterministic correction of order $1/\log n$
plus random $O(1/\log n)$ fluctuations. The correct finite-$n$ statement has two parts: the limit theorem above, plus
$n\,d - C\log n = O_p(1)$, which is what the "centered" plot tests.

---

## 4. The picture as n → ∞ (second order)

The natural refinement, as in Janson, is a distributional limit:
$$n\,d_{\text{pair}}-\log n\;\Rightarrow\;X_{\text{pair}},\qquad n\,d_{\text{ecc}}-3\log n\Rightarrow X_{\text{ecc}},\qquad n\,d_{\text{diam}}-5\log n\Rightarrow X_{\text{diam}}$$
with non-degenerate limits. In 1D these are sums of Gumbel variables. A mean-field ODE for the fraction $p$ of reached edges,
$\dot p=n\,p(1-p)(2-p)$ (the 2D analogue of Janson's logistic equation $\dot p = np(1-p)$), gives approximate constants
$$\mathbb E\,n d_{\text{pair}}\approx\tfrac12\log\tfrac N2+0.635+\log2,\qquad
\mathbb E\,n d_{\text{ecc}}\approx\tfrac32\log N-\tfrac12\log2+0.635+\gamma .$$
These are heuristic. The ODE ignores the correlation between an edge and its triangle-neighbours, so the $O(1)$
constants may be slightly off. The leading constants are not affected.

**Why the diameter converges slowly from below.** $\log N=2\log n-\log2$, and the maximum of the escape times over $N$
edges sits at $\log N/n$, not $2\log n/n$. Both ends contribute a negative $O(1)$ term. So
$\text{diam}\cdot n/\log n\approx5-c/\log n$ with $c>0$. With $\log n\approx8$ at $n=3000$, a value near 4.8 is what
the theory predicts, and the gap shrinks only like $1/\log n$.

---

## 5. What the simulations say (run of 2026-10-25)

Full tables: `results_dim2.md`, `results_dim1.md`. Graphs: `plots/dim2_dashboard.png`, `plots/dim2_centered.png`.

**Control (1D, n up to 100 000).** The extrapolated limits are $1.006\pm0.08,\ 2.041\pm0.07,\ 3.027\pm0.07$, which
recover Janson's theorem (1, 2, 3). The pipeline and the extrapolation method are sound.

**2D (n up to 5000; diameter up to 3000).** Extrapolated limits of the normalised distance:

| fit | pair | ecc | diam |
|---|---|---|---|
| naive: $d\,n/\log n = a+b/\log n$ | 0.87 ± 0.14 | 2.80 ± 0.16 | 4.80 ± 0.18 |
| **natural 2D scale:** $(n-2)\,d\,/\,(\tfrac12\log N) = a + b/\log N$ | **0.93 ± 0.14** | **2.94 ± 0.16** | **5.02 ± 0.17** |
| predicted | 1 | 3 | 5 |

In 2D the true scales are $n-2$ (triangles per edge) and $\tfrac12\log N=\log n-\tfrac12\log2+O(1/n)$. With the naive
normalisation these produce $O(\log n/n)$ corrections that a pure $b/\log n$ fit cannot absorb, which pulls the
intercept down. With the natural normalisation all three limits are consistent with (1, 3, 5). (The same change
leaves 1D essentially unchanged.)

**Second-order check** (`plots/dim2_centered.png`): the mean of $n\,d-C\log n$ is flat in $n$ within error bars:
pair ≈ +0.8, ecc ≈ 0, diam ≈ −1.6. A wrong $C$ would make this curve grow like $\pm\log n$. So the data fit
$$\frac{n\,d}{\log n}=C+\frac{c}{\log n}+o\!\left(\frac1{\log n}\right),\qquad (C,c)\approx(1,+0.8),\ (3,0),\ (5,-1.6).$$
The mean-field formulas of §4 match pair and ecc at every $n$ within the CI (e.g. ecc at $n=50$: data 3.102, MF 3.071;
at $n=3000$: 2.964 vs 2.980). In particular they explain the slow drift of ecc from 3.10 down toward $3-0.17/\log n$.

**Bottom line as $n\to\infty$.**
$$\frac{n}{\log n}\,d_{\text{pair}}\to1,\qquad\frac{n}{\log n}\,d_{\text{ecc}}\to3,\qquad\frac{n}{\log n}\,d_{\text{diam}}\to5\quad\text{(in probability)},$$
with $O(1)$ fluctuations of $n\,d - C\log n$. Compared with 1D: the typical distance is unchanged, and each "slow end"
costs twice as much, because there are $\sim n^2$ edges instead of $n$ vertices and each still escapes at rate $n$.

**Suggested next steps.**
1. Prove the pair lower bound of §3 carefully (it is nearly complete), then the escape-time lower bounds for ecc/diam.
   These are the cleanest parts.
2. Upper bounds: adapt Janson's exploration argument to the jump-$+2$ Yule process on $J(n,2)$.
3. Experiment: 3D (points = triangles, weights on tetrahedra). Predicted **1, 4, 7**. This tests the general-$k$ formula
   $1,\ 1+k,\ 1+2k$.
4. More reps at $n=2000$–$5000$ to shrink the error bars on the large-$n$ points.
