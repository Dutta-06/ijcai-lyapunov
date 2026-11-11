# Proposition D‴.56: the transit under measurement noise (`exp_noisy_certificate.py`)

The sliding path on y = f* + ξ (21 Gaussian kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]) against its first-order response along the noiseless path (`path_response`); targets: three features at least 0.45 apart (D‴.55's recipe), certificate margin m > 0.8q. The noisy paths keep their centres in the samples' window ± 2σ (a relaxation that carries one out stops the path: runaway) and judge exits and stationarity on the data's scale max|y|; a kernel exits only if its correlation stays inward with its own contribution removed.

## A. The response is the first order

Per target (one noise direction ξ, unit variance per sample) and level h: the noisy path's entry times and its duals at the K-th entry against the first-order prediction — the largest error, and the first-order term's own size.

| target | h | same entries | entry times: error (first-order term) | duals: error (first-order term) |
|---|---|---|---|---|
| 0 | 0.001 | True | 1.26e-05 (1.24e-02) | 4.81e-07 (6.57e-04) |
| 0 | 0.0001 | True | 1.26e-07 (1.24e-03) | 4.82e-09 (6.57e-05) |
| 0 | 1e-05 | True | 1.26e-09 (1.24e-04) | 4.82e-11 (6.57e-06) |
| 1 | 0.001 | True | 2.06e-05 (1.84e-02) | 5.03e-07 (1.13e-03) |
| 1 | 0.0001 | True | 2.06e-07 (1.84e-03) | 5.02e-09 (1.13e-04) |
| 1 | 1e-05 | True | 2.06e-09 (1.84e-04) | 5.03e-11 (1.13e-05) |
| 2 | 0.001 | True | 1.43e-05 (1.53e-02) | 1.26e-06 (1.33e-03) |
| 2 | 0.0001 | True | 1.43e-07 (1.53e-03) | 1.26e-08 (1.33e-04) |
| 2 | 1e-05 | True | 1.43e-09 (1.53e-04) | 1.26e-10 (1.33e-05) |
| 3 | 0.001 | True | 3.06e-06 (7.89e-03) | 8.54e-07 (1.78e-03) |
| 3 | 0.0001 | True | 3.07e-08 (7.89e-04) | 8.54e-09 (1.78e-04) |
| 3 | 1e-05 | True | 3.07e-10 (7.89e-05) | 8.54e-11 (1.78e-05) |
| 4 | 0.001 | True | 7.64e-06 (1.06e-02) | 3.31e-07 (4.79e-04) |
| 4 | 0.0001 | True | 7.64e-08 (1.06e-03) | 3.30e-09 (4.79e-05) |
| 4 | 1e-05 | True | 7.64e-10 (1.06e-04) | 3.30e-11 (4.79e-06) |
| 5 | 0.001 | True | 3.55e-06 (6.23e-03) | 1.42e-07 (3.89e-04) |
| 5 | 0.0001 | True | 3.56e-08 (6.23e-04) | 1.42e-09 (3.89e-05) |
| 5 | 1e-05 | True | 3.56e-10 (6.23e-05) | 1.42e-11 (3.89e-06) |
| 6 | 0.001 | True | 3.87e-06 (6.60e-03) | 9.74e-08 (1.35e-03) |
| 6 | 0.0001 | True | 3.87e-08 (6.60e-04) | 9.74e-10 (1.35e-04) |
| 6 | 1e-05 | True | 3.87e-10 (6.60e-05) | 9.74e-12 (1.35e-05) |
| 7 | 0.001 | True | 3.16e-05 (2.55e-02) | 1.83e-06 (1.66e-03) |
| 7 | 0.0001 | True | 3.16e-07 (2.55e-03) | 1.83e-08 (1.66e-04) |
| 7 | 1e-05 | True | 3.16e-09 (2.55e-04) | 1.83e-10 (1.66e-05) |
| 8 | 0.001 | True | 2.84e-05 (2.08e-02) | 1.74e-06 (9.28e-04) |
| 8 | 0.0001 | True | 2.84e-07 (2.08e-03) | 1.74e-08 (9.28e-05) |
| 8 | 1e-05 | True | 2.84e-09 (2.08e-04) | 1.74e-10 (9.28e-06) |
| 9 | 0.001 | True | 1.42e-05 (1.33e-02) | 6.08e-07 (1.55e-03) |
| 9 | 0.0001 | True | 1.43e-07 (1.33e-03) | 6.08e-09 (1.55e-04) |
| 9 | 1e-05 | True | 1.43e-09 (1.33e-04) | 6.08e-11 (1.55e-05) |

The duals' error falls by 99.8–100.1 per decade of h (second order: 100).

## B. The order and the support

40 targets (40 with the noiseless path exactly E), 10 draws per level. The radius r₁ = min (1 − |u⁰|)/|U| over the order criterion's rows: 0.0085 … 0.2781 (median 0.1038); over the support criterion's alone: 0.0316 … 0.3644; the targets whose consecutive entries lie within 5% (τ_(k+1)/τ_k < 1.05): 10; draws whose K-th block ran away (not exact): 0.

Observed: the path enters exactly E first (as a set; in the noiseless order — the swaps the difference); predicted per draw: the order criterion max|u⁰_j + U_jξ| < 1 over the entries and the kernels not yet entered (E's later entrants included), and the support criterion over the kernels off E alone; the linear model's expected count of the order kept (Monte Carlo) ± its sd and the observed count's z; the union bound on the expected failures (Σ over targets of min(1, Σ_rows 2Φ̄)·draws) against the observed; disagreements of the order criterion, and of these the ones within the level's largest linear-dual error of 1.

| s | exact: set, order (swaps) | criterion < 1: order, support | agree per draw: order, support | expected order (linear model) ± sd, z | union bound: order [support] vs failures | disagreements (within the linear error of 1) |
|---|---|---|---|---|---|---|
| 0.003 | 400, 400 (0) | 400, 400 | 400, 400 | 400.0 ± 0.2, +0.2 | 0.1 [0.0] vs 0 [0] | 0 (0) |
| 0.01 | 400, 398 (2) | 398, 400 | 400, 400 | 396.1 ± 1.8, +1.0 | 7.7 [0.0] vs 2 [0] | 0 (0) |
| 0.03 | 399, 373 (26) | 375, 398 | 398, 399 | 380.5 ± 3.8, -2.0 | 40.0 [4.3] vs 27 [1] | 2 (2) |
| 0.1 | 361, 297 (64) | 301, 352 | 392, 381 | 302.0 ± 7.9, -0.6 | 219.6 [110.3] vs 103 [39] | 8 (8) |
| 0.3 | 208, 127 (81) | 141, 177 | 372, 335 | 115.3 ± 8.8, +1.3 | 400.0 [400.0] vs 273 [192] | 28 (28) |

The linear duals against the path's own (`path_duals`) at its first K entries, the kernels off E: their largest difference on the draws that keep the order (and over s²); on the swaps, the support criterion against the path's own largest dual off E.

| s | order kept: max error | / s² | swaps: max |crit − path's| |
|---|---|---|---|
| 0.003 | 4.59e-05 | 5.10 | — |
| 0.01 | 7.40e-04 | 7.40 | 0.031 |
| 0.03 | 6.73e-03 | 7.48 | 0.048 |
| 0.1 | 4.70e-02 | 4.70 | 0.193 |
| 0.3 | 2.54e-01 | 2.82 | 0.359 |

## B2. Small radius (r₁ in [0.002, 0.02])

20 targets (15 with the noiseless path exactly E), 20 draws per level. The radius r₁ = min (1 − |u⁰|)/|U| over the order criterion's rows: 0.0026 … 0.0192 (median 0.0128); over the support criterion's alone: 0.0026 … 0.2922; the targets whose consecutive entries lie within 5% (τ_(k+1)/τ_k < 1.05): 5; draws whose K-th block ran away (not exact): 0.

Observed: the path enters exactly E first (as a set; in the noiseless order — the swaps the difference); predicted per draw: the order criterion max|u⁰_j + U_jξ| < 1 over the entries and the kernels not yet entered (E's later entrants included), and the support criterion over the kernels off E alone; the linear model's expected count of the order kept (Monte Carlo) ± its sd and the observed count's z; the union bound on the expected failures (Σ over targets of min(1, Σ_rows 2Φ̄)·draws) against the observed; disagreements of the order criterion, and of these the ones within the level's largest linear-dual error of 1.

| s | exact: set, order (swaps) | criterion < 1: order, support | agree per draw: order, support | expected order (linear model) ± sd, z | union bound: order [support] vs failures | disagreements (within the linear error of 1) |
|---|---|---|---|---|---|---|
| 0.001 | 400, 400 (0) | 400, 400 | 400, 400 | 399.8 ± 0.4, +0.4 | 0.5 [0.5] vs 0 [0] | 0 (0) |
| 0.003 | 389, 389 (0) | 389, 389 | 400, 400 | 388.3 ± 3.1, +0.2 | 37.2 [37.2] vs 11 [11] | 0 (0) |
| 0.01 | 346, 344 (2) | 343, 345 | 399, 399 | 337.5 ± 6.7, +1.0 | 173.2 [164.4] vs 56 [54] | 1 (1) |

The linear duals against the path's own (`path_duals`) at its first K entries, the kernels off E: their largest difference on the draws that keep the order (and over s²); on the swaps, the support criterion against the path's own largest dual off E.

| s | order kept: max error | / s² | swaps: max |crit − path's| |
|---|---|---|---|
| 0.001 | 4.93e-06 | 4.93 | — |
| 0.003 | 6.89e-05 | 7.65 | — |
| 0.01 | 6.03e-04 | 6.03 | 0.000 |

## C. The window

On the draws that keep the order with the order criterion below 1 (the first order's regime): the first noise kernel's entry after the K-th, τ_{K+1} − τ_K, against the first-order W — the ratio's median, range and 5–95% quantiles, and their half-width over s (constant if the relative error is O(s)) —; the kernel as predicted, and on the draws where it is not, the relative gap between the two smallest predicted windows (W₍₂₎ − W₍₁₎)/W₍₁₎ (median; all draws' median in brackets); W·s (median, quartiles) against its limit law (`law_W`, 20000 samples per target; pooled over the targets: median 14, quartiles [10.1, 19.7]) — the Kolmogorov–Smirnov p-value of the draws' probability integral transforms under their targets' laws —; outside the regime: the draws that keep the order with the criterion ≥ 1, and those with the support exact but the order swapped (of each, those with W_pred < 0: the first-order dual already past its wall).

| s | draws | W_obs/W_pred: median [min, max] | 5–95% | half-width / s | next kernel as predicted | missed: (W₍₂₎ − W₍₁₎)/W₍₁₎ [all] | W·s: median [quartiles] | KS p | order kept, criterion ≥ 1 (W_pred < 0) | swaps (W_pred < 0) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.003 | 400 | 1.0001 [0.983, 1.030] | [0.9928, 1.0106] | 2.96 | 398/400 | 0.00783 [0.116] | 13.5 [9.73, 19.1] | 0.28 | 0 (0) | 0 (0) |
| 0.01 | 398 | 1.0002 [0.923, 1.165] | [0.9792, 1.0397] | 3.02 | 387/398 | 0.0143 [0.142] | 14.3 [9.99, 19.9] | 0.88 | 0 (0) | 2 (0) |
| 0.03 | 373 | 1.0010 [0.856, 1.361] | [0.9393, 1.1206] | 3.02 | 349/373 | 0.0386 [0.14] | 13 [9.28, 18.5] | 0.03 | 0 (0) | 26 (0) |
| 0.1 | 295 | 1.0078 [0.766, 3.801] | [0.9078, 1.5595] | 3.26 | 239/295 | 0.105 [0.161] | 12.7 [8.93, 17.3] | 5.5e-05 | 2 (2) | 64 (0) |
| 0.3 | 120 | 1.0801 [0.475, 87.782] | [0.7107, 5.7875] | 8.46 | 88/120 | 0.507 [0.335] | 8.96 [5.88, 13.1] | 4.5e-14 | 7 (3) | 81 (4) |

## D. Stopping and the oracle

The discrepancy rule: the first entry after which the relaxed residual's rms is ≤ 1.1s, on the draws with the support exact: at the K-th, not by the K-th (late), earlier; its first-order prediction — at the K-th iff rms(R_Kξ) ≤ 1.1s, R_Kξ the first-order residual (|R_Kξ|²/s² ~ χ²_(N−2K)) — and its agreement per draw; the path's late stops against the first order's binomial law, P(χ²_(N−2K) > 1.1²N) = 0.0461 per draw (expected ± sd, the binomial test's p); the relaxed residual after the K-th entry against rms(R_Kξ).

| s | draws | at K | late | early | first order: at K | agree per draw | late: expected ± sd, p | rms_K / rms(R_Kξ): median [5–95%] |
|---|---|---|---|---|---|---|---|---|
| 0.003 | 400 | 385 | 15 | 0 | 385 | 400 | 18.4 ± 4.2, 0.48 | 1.0000 [0.9999, 1.0001] |
| 0.01 | 400 | 388 | 12 | 0 | 388 | 400 | 18.4 ± 4.2, 0.15 | 1.0000 [0.9996, 1.0003] |
| 0.03 | 399 | 379 | 20 | 0 | 379 | 399 | 18.4 ± 4.2, 0.63 | 1.0000 [0.9989, 1.0009] |
| 0.1 | 361 | 343 | 18 | 0 | 343 | 361 | 16.6 ± 4.0, 0.71 | 1.0000 [0.9965, 1.0033] |
| 0.3 | 208 | 111 | 8 | 89 | 200 | 119 | 9.6 ± 3.0, 0.74 | 1.0004 [0.9946, 1.0090] |

The largest relative gap (τ_{k+1} − τ_k)/τ_k over the first K + 1 entries falling after the K-th; the K-th entry's fit (the K-kernel least-squares fit of y): its rms error against the noise's projection on the tangent space, rms(P_Tξ), and against s(2K/N)^{1/2}.

| s | draws | largest gap at K | error / rms(P_Tξ): median [range] | 5–95% | error / (s(2K/N)^{1/2}): median |
|---|---|---|---|---|---|
| 0.003 | 400 | 400/400 | 1.0001 [0.994, 1.004] | [0.9982, 1.0019] | 0.949 |
| 0.01 | 400 | 400/400 | 1.0001 [0.984, 1.017] | [0.9939, 1.0056] | 0.953 |
| 0.03 | 399 | 399/399 | 1.0005 [0.938, 1.049] | [0.9839, 1.0213] | 0.968 |
| 0.1 | 361 | 359/361 | 0.9986 [0.894, 1.195] | [0.9478, 1.0548] | 0.921 |
| 0.3 | 208 | 181/208 | 0.9865 [0.675, 1.402] | [0.8414, 1.1274] | 0.831 |

## E. Past the first noise kernel

Four families: 'base' (B's first targets), 'wide' (amplitudes 1, 0.3, 0.1 in random order and signs), 'K5' (five features), 'long' (B's first targets, the horizon K + 8); each path to its (K + more)-th entry under 3e+06 right-hand-side evaluations. Paths: finished (out of the budget), stopped by a runaway kernel before the horizon; the first K entries E as a set; the largest relative gap over the first K + 1, K + 4 (and K + 8) entries at the K-th (on the exact paths); its condition, W/τ_K above the features' own largest relative gap and above the later ones (on the exact paths that reach the horizon); the discrepancy rule over the relaxed fits (at the K-th fit, later, earlier, never); the fits the two rules stop at — the gap rule over the longest horizon — against the oracle (the K-kernel least-squares fit of y on the true support) on every finished path; the stopped fits unrelaxed or with a centre outside the samples' window. Later exits: paths with an exit after the (K + 1)-th entry.

| family | s | paths (budget) | runaway, later exits | exact | largest gap at K: K+1, K+4[, K+8] | W/τ_K > features' gaps, > later gaps | discrepancy: at K, later, early, never | gap-stopped error / oracle: median [max] | discrepancy-stopped: median [max] | stopped fits unrelaxed or outside: gap, discrepancy |
|---|---|---|---|---|---|---|---|---|---|---|
| base | 0.01 | 94 (6) | 19, 0 | 94 | 94/94, 75/75 | 75, 75 of 75 | 92, 1, 0, 1 | 1.000 [1] | 1.000 [1.69] | 0, 0 |
| base | 0.03 | 94 (6) | 26, 0 | 94 | 94/94, 68/68 | 68, 68 of 68 | 91, 1, 0, 2 | 1.000 [1] | 1.000 [1.22] | 0, 0 |
| base | 0.1 | 81 (19) | 20, 0 | 76 | 75/76, 52/59 | 58, 52 of 59 | 78, 2, 0, 1 | 1.000 [1.91] | 1.000 [1.44] | 1, 0 |
| base | 0.3 | 69 (31) | 22, 0 | 41 | 33/41, 16/30 | 23, 16 of 30 | 36, 3, 28, 2 | 1.000 [3] | 1.000 [5.06] | 8, 1 |
| wide | 0.01 | 19 (1) | 4, 0 | 18 | 18/18, 14/14 | 14, 14 of 14 | 19, 0, 0, 0 | 1.000 [1] | 1.000 [1] | 0, 0 |
| wide | 0.03 | 18 (2) | 3, 0 | 15 | 10/15, 8/12 | 8, 10 of 12 | 16, 1, 0, 1 | 1.000 [11] | 1.000 [1.18] | 0, 1 |
| wide | 0.1 | 13 (7) | 5, 0 | 1 | 0/1, 0/1 | 0, 0 of 1 | 0, 0, 13, 0 | 1.821 [5.38] | 1.415 [1.94] | 0, 0 |
| K5 | 0.01 | 11 (1) | 5, 0 | 11 | 11/11, 6/6 | 6, 6 of 6 | 11, 0, 0, 0 | 1.000 [1] | 1.000 [1] | 0, 0 |
| K5 | 0.03 | 10 (2) | 4, 0 | 9 | 9/9, 4/5 | 5, 4 of 5 | 10, 0, 0, 0 | 1.000 [1.03] | 1.000 [1] | 0, 0 |
| long | 0.03 | 9 (6) | 8, 0 | 9 | 9/9, 5/5, 1/1 | 1, 1 of 1 | 9, 0, 0, 0 | 1.000 [1] | 1.000 [1] | 0, 0 |
| long | 0.1 | 7 (8) | 6, 0 | 7 | 7/7, 6/6, 1/1 | 1, 1 of 1 | 7, 0, 0, 0 | 1.000 [1.69] | 1.000 [1] | 1, 0 |

The features' own relative gaps (τ_{k+1} − τ_k)/τ_k, k < K, and the later ones (k > K), on the exact paths that reach the horizon: their largest (median [max]), against W/τ_K (median [min]).

| family | s | features' largest | later largest | W/τ_K |
|---|---|---|---|---|
| base | 0.01 | 0.304 [0.912] | 1.29 [7.38] | 72.2 [34.2] |
| base | 0.03 | 0.334 [0.941] | 1.31 [3.59] | 23.9 [5.71] |
| base | 0.1 | 0.334 [0.816] | 1.76 [12.1] | 6.94 [0.0987] |
| base | 0.3 | 0.37 [1.49] | 1.4 [9.35] | 1.58 [0.00815] |
| wide | 0.01 | 2.38 [2.48] | 2.28 [4.14] | 11.9 [4.41] |
| wide | 0.03 | 2.44 [2.66] | 0.997 [2.06] | 4.87 [0.955] |
| wide | 0.1 | 3.3 [3.3] | 2.07 [2.07] | 0.502 [0.502] |
| K5 | 0.01 | 0.334 [0.639] | 1.9 [8.25] | 111 [10] |
| K5 | 0.03 | 0.401 [0.679] | 3.51 [19.6] | 38.9 [10.3] |
| long | 0.03 | 0.161 [0.161] | 0.942 [0.942] | 101 [101] |
| long | 0.1 | 0.354 [0.354] | 1.93 [1.93] | 28.3 [28.3] |

Numbers are produced by `experiments/exp_noisy_certificate.py`; nothing is edited by hand.
