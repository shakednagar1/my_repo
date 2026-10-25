# Results: 2D Janson experiment
Updated 2026-10-25 13:39. Values are d·n/log n, mean ± 95% CI. Raw runs: data/runs_dim2.csv

| n | reps | fixed pair | fixed → furthest | diameter |
|---|---|---|---|---|
| 50 | 200 | 1.249 ± 0.057 | 3.102 ± 0.056 | 4.702 ± 0.062 |
| 100 | 100 | 1.212 ± 0.076 | 3.021 ± 0.068 | 4.639 ± 0.073 |
| 200 | 60 | 1.146 ± 0.079 | 2.983 ± 0.074 | 4.701 ± 0.090 |
| 500 | 30 | 1.138 ± 0.092 | 2.974 ± 0.086 | 4.737 ± 0.097 |
| 1000 | 15 | 1.090 ± 0.118 | 3.080 ± 0.166 | 4.756 ± 0.174 |
| 2000 | 8 | 1.176 ± 0.291 | 2.951 ± 0.155 | 4.725 ± 0.152 |
| 3000 | 5 | 1.001 ± 0.164 | 2.964 ± 0.200 | 4.808 ± 0.209 |
| 5000 | 4 | 1.033 ± 0.128 | 2.917 ± 0.283 | — |

## Extrapolation  d·n/log n = a + b/log n

| metric | a (n → ∞) | b | predicted | log-log slope |
|---|---|---|---|---|
| fixed pair | 0.873 ± 0.142 | 1.49 | 1 | -0.872 |
| fixed → furthest | 2.798 ± 0.159 | 1.12 | 3 | -0.842 |
| diameter | 4.797 ± 0.177 | -0.48 | 5 | -0.822 |

Plots: plots/dim2_dashboard.png (all in one), plots/dim2_normalised.png, plots/dim2_extrapolation.png, plots/dim2_loglog.png, plots/dim2_fluctuations.png, plots/dim2_centered.png
