# Proposition D‴.36: the gauge under dither (`exp_gauge_dither.py`)

Proposition D‴.32's instance, one particle per class at ϑ^e(1/3) (m = 3), the replicator at rate γ, the positions dithered √(2τG), τ = 0.0001, Euler–Maruyama dt = 0.25; the masses W_k = β_k/3.

## A. The predicted diffusion of the masses along the manifold of exact fits (W = 1/3, per unit slow time)

| γ | D_11, D_22, D_33 | D_12, D_13, D_23 |
|---|---|---|
| 0.001 | 4.85e-10, 2.24e-10, 2.62e-10 | -2.23e-10, -2.62e-10, -7.17e-13 |
| 0.01 | 2.85e-8, 1.63e-8, 1.81e-8 | -1.33e-8, -1.51e-8, -2.95e-9 |
| 0.1 | 2.56e-7, 3.47e-7, 3.34e-7 | -1.35e-7, -1.21e-7, -2.13e-7 |
| 1 | 4.01e-7, 9.81e-7, 8.96e-7 | -2.43e-7, -1.58e-7, -7.38e-7 |
| 10 | 4.23e-7, 1.14e-6, 1.03e-6 | -2.64e-7, -1.59e-7, -8.73e-7 |

## B. The diffusion measured (64 seeds, γ = 0.01, 50000 slow units)

| s | mean-square displacement of W_1, W_2, W_3 | the displacement's variance |
|---|---|---|
| 5000 | 1.54e-4, 7.91e-5, 1.26e-4 | 1.53e-4, 8.03e-5, 1.25e-4 |
| 10000 | 3.07e-4, 1.69e-4, 2.80e-4 | 3.06e-4, 1.71e-4, 2.78e-4 |
| 20000 | 5.38e-4, 4.53e-4, 4.67e-4 | 5.34e-4, 4.56e-4, 4.73e-4 |
| 30000 | 7.96e-4, 6.70e-4, 6.01e-4 | 7.66e-4, 6.49e-4, 6.09e-4 |
| 50000 | 0.00151, 0.00116, 9.13e-4 | 0.00138, 0.00105, 9.27e-4 |

The variance's slope over s ∈ [10000, 50000]: 2.69e-8 ± 5.82e-9, 2.21e-8 ± 4.65e-9, 1.62e-8 ± 3.76e-9 (bootstrap standard errors); predicted (A): 2.85e-8, 1.63e-8, 1.81e-8; ratios 0.95, 1.36, 0.90.

The masses' mean at s = 50000 (from 1/3): 0.3212 ± 0.0046, 0.3445 ± 0.0041, 0.3342 ± 0.0038 (standard errors): the drift term of the limit, which (A) does not include.

## C–D. Long runs: no leak, the control γ = 0, and the leak on the masses (2·10⁵ slow units)

| γ | ε | seeds | the masses over the second half: mean | std | predicted std (the leak's Ornstein–Uhlenbeck limit at W = 1/3) | the amplitude coordinates' range (α_1, α_2, α_3) | the largest |a|/a_max | fit error (second-half mean) |
|---|---|---|---|---|---|---|---|---|
| 0.01 | 0 | 10 | 0.3141, 0.3673, 0.3186 | 0.0571, 0.0527, 0.0519 | — | [+2.32, +16.11], [-3.69, -1.52], [+1.93, +9.45] | 0.9994 | 0.0253 |
| 0.01 | 0.0001 | 4 | 0.3252, 0.3397, 0.3351 | 0.0152, 0.0103, 0.0117 | 0.0139, 0.0099, 0.0105 | [+3.36, +5.66], [-2.88, -2.00], [+2.35, +3.27] | 0.8887 | 0.0264 |
| 0.01 | 0.001 | 4 | 0.3329, 0.3338, 0.3333 | 0.0044, 0.0032, 0.0034 | 0.0044, 0.0031, 0.0033 | [+3.50, +4.57], [-2.89, -2.07], [+2.51, +3.13] | 0.8156 | 0.0265 |
| 0 | 0 | 4 | 0.3333, 0.3333, 0.3333 | 0.0000, 0.0000, 0.0000 | — | [+3.58, +4.43], [-2.83, -2.06], [+2.53, +3.08] | 0.8033 | 0.0269 |

Numbers are produced by `experiments/exp_gauge_dither.py`; nothing is edited by hand.
