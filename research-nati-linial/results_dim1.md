# Results: 1D Janson experiment
Updated 2026-10-25 12:59. Values are d·n/log n, mean ± 95% CI. Raw runs: data/runs_dim1.csv

| n | reps | fixed pair | fixed → furthest | diameter |
|---|---|---|---|---|
| 100 | 400 | 1.127 ± 0.045 | 2.220 ± 0.038 | 3.044 ± 0.037 |
| 300 | 300 | 1.123 ± 0.044 | 2.211 ± 0.041 | 3.070 ± 0.038 |
| 1000 | 200 | 1.065 ± 0.046 | 2.158 ± 0.036 | 3.056 ± 0.033 |
| 3000 | 100 | 1.110 ± 0.051 | 2.172 ± 0.046 | 3.051 ± 0.043 |
| 10000 | 50 | 1.058 ± 0.075 | 2.125 ± 0.061 | 3.016 ± 0.059 |
| 30000 | 30 | 1.068 ± 0.075 | 2.117 ± 0.053 | 3.059 ± 0.073 |
| 100000 | 15 | 1.009 ± 0.120 | 2.097 ± 0.092 | 3.008 ± 0.070 |

## Extrapolation  d·n/log n = a + b/log n

| metric | a (n → ∞) | b | predicted | log-log slope |
|---|---|---|---|---|
| fixed pair | 1.006 ± 0.082 | 0.58 | 1 | -0.883 |
| fixed → furthest | 2.041 ± 0.066 | 0.87 | 2 | -0.878 |
| diameter | 3.027 ± 0.065 | 0.15 | 3 | -0.871 |

Plots: plots/dim1_dashboard.png (all in one), plots/dim1_normalised.png, plots/dim1_extrapolation.png, plots/dim1_loglog.png, plots/dim1_fluctuations.png, plots/dim1_centered.png
