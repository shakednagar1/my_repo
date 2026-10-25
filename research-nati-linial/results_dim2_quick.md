# Results: 2D Janson experiment
Updated 2026-10-25 12:53. Values are d·n/log n, mean ± 95% CI. Raw runs: data/runs_dim2_quick.csv

| n | reps | fixed pair | fixed → furthest | diameter |
|---|---|---|---|---|
| 30 | 20 | 1.522 ± 0.305 | 3.296 ± 0.293 | 4.931 ± 0.284 |
| 60 | 10 | 1.258 ± 0.216 | 3.000 ± 0.224 | 4.703 ± 0.224 |
| 120 | 5 | 1.262 ± 0.250 | 3.046 ± 0.068 | 4.700 ± 0.326 |

## Extrapolation  d·n/log n = a + b/log n

| metric | a (n → ∞) | b | predicted | log-log slope |
|---|---|---|---|---|
| fixed pair | 0.583 ± 1.136 | 3.02 | 1 | -0.888 |
| fixed → furthest | 2.608 ± 0.692 | 2.07 | 3 | -0.810 |
| diameter | 4.014 ± 1.255 | 3.01 | 5 | -0.788 |

Plots: plots/dim2_quick_dashboard.png (all in one), plots/dim2_quick_normalised.png, plots/dim2_quick_extrapolation.png, plots/dim2_quick_loglog.png, plots/dim2_quick_fluctuations.png, plots/dim2_quick_centered.png
