# Factorial: robustness mechanism × geometry (Lorenz tracking, m=27, k_e=10)

Seeds: 5; metric = steady-state ‖e‖ over the last 40% of a 30s run; mean ± 95% CI over seeds (max over window in brackets); 'div' = diverged runs.

| geometry | mechanism | mean SS ‖e‖ | median SS | max SS | max ‖W‖_F | min sep | div |
|---|---|---|---|---|---|---|---|
| fixed (W) | projection | 0.1584 ± 0.0072 | 0.1288 | 0.559 | 47.5 | 1.365 | 0 |
| fixed (W) | sigma-mod | 0.7327 ± 0.0223 | 0.5647 | 2.720 | 38.5 | 1.469 | 0 |
| fixed (W) | e-mod | 0.5619 ± 0.0306 | 0.4264 | 2.192 | 26.7 | 1.469 | 0 |
| adaptive (W,c) | projection | 0.1468 ± 0.0038 | 0.1020 | 0.751 | 40.9 | 0.039 | 0 |
| adaptive (W,c) | sigma-mod | 0.4004 ± 0.0050 | 0.3468 | 1.863 | 32.3 | 0.078 | 0 |
| adaptive (W,c) | e-mod | 0.2327 ± 0.0024 | 0.1993 | 0.718 | 28.1 | 0.099 | 0 |
| adaptive (W,c,σ) | projection | 0.1119 ± 0.0043 | 0.0852 | 0.512 | 38.5 | 0.172 | 0 |
| adaptive (W,c,σ) | sigma-mod | 0.3862 ± 0.0057 | 0.3422 | 1.655 | 30.9 | 0.215 | 0 |
| adaptive (W,c,σ) | e-mod | 0.2385 ± 0.0029 | 0.2103 | 0.820 | 27.7 | 0.239 | 0 |

## Paired (same seed) differences, projection minus baseline, mean SS ‖e‖

| geometry | baseline | mean diff | 95% CI | projection better in |
|---|---|---|---|---|
| fixed (W) | sigma-mod | -0.5743 | ±0.0232 | 5/5 seeds |
| fixed (W) | e-mod | -0.4035 | ±0.0310 | 5/5 seeds |
| adaptive (W,c) | sigma-mod | -0.2536 | ±0.0021 | 5/5 seeds |
| adaptive (W,c) | e-mod | -0.0858 | ±0.0025 | 5/5 seeds |
| adaptive (W,c,σ) | sigma-mod | -0.2743 | ±0.0043 | 5/5 seeds |
| adaptive (W,c,σ) | e-mod | -0.1266 | ±0.0031 | 5/5 seeds |

## Geometry effect within the projection mechanism (paired, mean SS ‖e‖)

| comparison | mean diff | 95% CI | adaptive better in |
|---|---|---|---|
| adaptive (W,c) − fixed (W) | -0.0116 | ±0.0056 | 5/5 |
| adaptive (W,c,σ) − adaptive (W,c) | -0.0349 | ±0.0026 | 5/5 |
