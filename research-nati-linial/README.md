# Janson 2D: random shortest paths on triangles of K_n

Research with Nati Linial on a 2-dimensional version of Janson's 1999 theorem.
Each triangle of K_n gets an independent Exp(1) weight; the points are the
edges of K_n, and a path is a chain of triangles sharing edges.
Conjecture: normalised by log n / n, pair -> 1, ecc -> 3, diam -> 5.

- `code/janson2d.py`: Numba Dijkstra simulation (`pip install -r code/requirements.txt`)
- `results_dim1.*`, `results_dim2.*`: simulation summaries (1D reproduces Janson's 1, 2, 3)
- `theory_notes.md`: proof notes
- `data/`, `logs/`, `plots/`: raw runs, logs and figures
