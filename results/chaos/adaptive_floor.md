# Proposition D‴.24: the finite-size floor is the leak's bias; an adaptive prior removes it without a mass constraint (`exp_adaptive_floor.py`)

The three-bump instance (σ_crit = 0.04716); the warm clouds of Proposition D‴.12 (width 0.02), m = 15, 60, 240 against M = 510, 4 seeds; midpoint rule dt = 0.25, the plan every 4 steps.

## A. The floor as the weight-derivative of the leak's bias (matching prior, no dynamics)

F: Proposition D‴.9's function-level sensitivity along e₁ − e₃, e₂ − e₃ (central differences of the fixed-prior nearest-atom stationary point); F_bias = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰. The law per unit of 1/m − 1/M, tr(FC̃Fᵀ)/N, from F and as the sum of the mode terms h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j; the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ).

| σ₀/σ_crit | ‖F − F_bias‖/‖F‖ | law from F | sum of the mode terms | bound | law/bound | floor (1/m − 1/M)·law at m = 15, 60, 240 |
|---|---|---|---|---|---|---|
| 0.05 | 1.1e-06 | 1.3559e-02 | 1.3559e-02 | 2.7310e-02 | 0.496 | 8.77e-04, 1.99e-04, 2.99e-05 |
| 0.15 | 1.9e-07 | 6.1662e-02 | 6.1662e-02 | 8.1930e-02 | 0.753 | 3.99e-03, 9.07e-04, 1.36e-04 |
| 0.5 | 1.2e-08 | 2.2913e-01 | 2.2913e-01 | 2.7310e-01 | 0.839 | 1.48e-02, 3.37e-03, 5.05e-04 |

The parameter level: the fixed prior's weight sensitivity S (whitened weighted) against H_L(H_L + σ₀I)⁻¹S^e with S^e = ∂_wϑ^e (per mode the fraction h_j/(h_j + σ₀) of the gauge's displacement, the rest σ₀/(h_j + σ₀) being the bias), and Proposition D‴.5's coupled-distance law per unit of 1/m − 1/M, Σ tr(S̃C̃S̃ᵀ), for the fixed prior and for the adaptive one (S^e: the whole displacement in the parameters).

| σ₀/σ_crit | ‖S − H_L(H_L + σ₀I)⁻¹S^e‖/‖S‖ | coupled-distance law, fixed prior | adaptive prior | ratio |
|---|---|---|---|---|
| 0.05 | 1.4e-07 | 3.2274e+01 | 4.6325e+01 | 1.435 |
| 0.15 | 9.2e-08 | 2.1321e+01 | 4.6325e+01 | 2.173 |
| 0.5 | 2.9e-08 | 8.7922e+00 | 4.6325e+01 | 5.269 |

Mode by mode (whitened weighted Gauss–Newton eigenvalues h_j; the share of each term in the law):

| σ₀/σ_crit | mode 1 | mode 2 | mode 3 | mode 4 | mode 5 | mode 6 | mode 7 | mode 8 | mode 9 |
|---|---|---|---|---|---|---|---|---|---|
| 0.05: h_j/σ₀ | 0.13 | 0.83 | 1.32 | 6.42 | 8.95 | 13.77 | 37.58 | 56.91 | 82.78 |
| 0.05: share | 0.000 | 0.109 | 0.313 | 0.261 | 0.220 | 0.091 | 0.005 | 0.000 | 0.000 |
| 0.15: h_j/σ₀ | 0.04 | 0.28 | 0.44 | 2.14 | 2.98 | 4.59 | 12.53 | 18.97 | 27.59 |
| 0.15: share | 0.000 | 0.050 | 0.179 | 0.321 | 0.302 | 0.140 | 0.009 | 0.000 | 0.000 |
| 0.5: h_j/σ₀ | 0.01 | 0.08 | 0.13 | 0.64 | 0.89 | 1.38 | 3.76 | 5.69 | 8.28 |
| 0.5: share | 0.000 | 0.019 | 0.078 | 0.316 | 0.359 | 0.208 | 0.020 | 0.001 | 0.000 |

The weight h_jσ₀²/(h_j + σ₀)² is a band-pass centred at h = σ₀: it is at least half its peak σ₀/4 for h/σ₀ ∈ [3 − 2√2, 3 + 2√2] = [0.17, 5.83]. The shares of the modes with the weight and without it (z_jᵀC̃z_j alone: the bound's distribution, independent of σ₀):

| σ₀/σ_crit | modes 2–6: weighted share | unweighted | modes with h_j/σ₀ ∈ [1/5, 5]: weighted share | unweighted | the weight relative to its peak, modes 1–6 |
|---|---|---|---|---|---|
| 0.05 | 0.994 | 0.972 | 0.422 | 0.213 | 0.41, 0.99, 0.98, 0.47, 0.36, 0.25 |
| 0.15 | 0.990 | 0.972 | 0.990 | 0.972 | 0.16, 0.68, 0.85, 0.87, 0.75, 0.59 |
| 0.5 | 0.979 | 0.972 | 0.902 | 0.784 | 0.05, 0.28, 0.41, 0.95, 1.00, 0.97 |

Since S^e = −∂_wx⁰, the bound is σ₀/4 times the adaptive prior's coupled-distance law per unit of 1/m − 1/M, tr(∂_wx⁰C̃∂_wx⁰ᵀ) = 46.3249 (independent of σ₀): the fixed prior's function-level law is at most σ₀/4 times the adaptive prior's parameter-level law.

The fixed matching prior's replicated stationary point Θ*(w) at the sampled weights (σ₀ = 0.15σ_crit): Proposition D‴.5's splitting margin min_k(σ₀ − λ_max R^G_k) over the classes with at least two particles (the class attaining it), and the smallest eigenvalue of the whitened weighted K-atom Hessian. Θ*(w) is a strict local minimiser (w ∈ W_adm, "D‴.5-admissible") iff both are positive; otherwise it is a saddle and a cluster can split. This is a different set from W, where an exact representation exists ("representable").

| m | seed 0 | seed 1 | seed 2 | seed 3 |
|---|---|---|---|---|
| 15 | -3.48e-03 (class 1), 7.4e-03 **saddle** | -1.08e-03 (class 1), 7.3e-03 **saddle** | +3.95e-03 (class 3), 7.4e-03 | -6.14e-03 (class 1), 7.5e-03 **saddle** |
| 60 | +5.61e-03 (class 2), 7.4e-03 | -2.16e-03 (class 1), 7.3e-03 **saddle** | +4.17e-03 (class 3), 7.4e-03 | -3.98e-03 (class 1), 7.3e-03 **saddle** |
| 240 | +3.78e-03 (class 3), 7.4e-03 | +5.06e-03 (class 3), 7.4e-03 | +5.95e-03 (class 3), 7.4e-03 | +5.13e-03 (class 1), 7.4e-03 |
| 510 | +4.40e-03 (class 3), 7.4e-03 | +5.44e-03 (class 3), 7.4e-03 | +6.62e-03 (class 3), 7.4e-03 | +4.52e-03 (class 2), 7.4e-03 |

Representable weights (the set W): an exact representation at the weights w needs |A_k|/w_k < a_max = 4, i.e. w > (0.250, 0.175, 0.200). Per seed and size: the sampled counts, the margin min_k(w_k − |A_k|/a_max), and the exact representation's gradient residual and fit error.

| m | seed 0 | seed 1 | seed 2 | seed 3 |
|---|---|---|---|---|
| 15 | (6, 5, 4): +0.0667, 1e-16, 2e-16 | (6, 3, 6): +0.0250, 3e-17, 2e-16 | (5, 4, 6): +0.0833, 2e-16, 2e-16 | (8, 4, 3): +0.0000 (none) |
| 60 | (19, 21, 20): +0.0667, 3e-17, 2e-16 | (24, 15, 21): +0.0750, 7e-14, 4e-13 | (19, 18, 23): +0.0667, 2e-17, 9e-17 | (25, 17, 18): +0.1000, 4e-17, 1e-16 |
| 240 | (71, 76, 93): +0.0458, 7e-17, 2e-16 | (77, 76, 87): +0.0708, 2e-16, 3e-16 | (76, 81, 83): +0.0667, 2e-17, 1e-16 | (82, 79, 79): +0.0917, 5e-14, 6e-13 |
| 510 | (158, 161, 191): +0.0598, 1e-14, 4e-14 | (163, 166, 181): +0.0696, 3e-17, 1e-16 | (170, 167, 173): +0.0833, 5e-16, 3e-15 | (168, 184, 158): +0.0794, 3e-16, 3e-16 |

## B. The nearest-atom leak at the matching prior: fixed (η = 0) against adaptive priors (6000 slow units)

‖f̂_m − f̂_M‖² (the function-level floor) and the fit error at s = 600, 3000, 6000, per seed; W²(μ_m, μ_M) at the end against the W² between the exact representations of the two systems' sampled weights (the measure-level distance; its split, and the coupled distance D, below); the class means' whitened weighted distance from the exact representation of the sampled weights at the end; the counts (all conserved). —: no exact representation at those weights (outside W or on its boundary).

| η | m | ‖f̂_m − f̂_M‖² at s = 600, 3000, 6000 (per seed) | fit error at 6000 (per seed) | W²(μ_m, μ_M) at the end / exact-representation W² (per seed) | distance from ϑ^e(w) at the end | counts conserved |
|---|---|---|---|---|---|---|
| 0.00000 | 15 | 1.8e-03, 1.7e-03, 1.7e-03; 2.5e-03, 2.4e-03, 2.4e-03; 3.0e-04, 3.0e-04, 3.0e-04; 5.5e-03, 5.5e-03, 5.5e-03 | 4.2e-02; 7.4e-02; 2.9e-02; 1.2e-01 | 4.48e+00 / 5.92e+00; 5.51e+00 / 7.27e+00; 2.22e+00 / 2.32e+00; 7.99e+00 / — | 3.7e-01; 9.5e-01; 2.0e-01; — | yes |
| 0.00000 | 60 | 1.1e-04, 1.1e-04, 1.1e-04; 1.1e-03, 9.9e-04, 9.9e-04; 1.6e-04, 1.6e-04, 1.6e-04; 7.3e-04, 7.1e-04, 7.1e-04 | 1.3e-02; 4.3e-02; 2.1e-02; 3.8e-02 | 1.57e+00 / 1.59e+00; 3.70e+00 / 4.16e+00; 1.80e+00 / 1.86e+00; 3.63e+00 / 3.95e+00 | 1.1e-01; 3.8e-01; 1.4e-01; 3.6e-01 | yes |
| 0.00000 | 240 | 8.2e-05, 8.2e-05, 8.2e-05; 4.9e-06, 4.8e-06, 4.8e-06; 7.3e-05, 7.3e-05, 7.3e-05; 7.5e-05, 7.5e-05, 7.5e-05 | 3.6e-02; 1.4e-02; 1.3e-02; 5.8e-03 | 6.81e-01 / 7.59e-01; 3.52e-01 / 3.55e-01; 8.34e-01 / 8.67e-01; 1.37e+00 / 1.37e+00 | 3.3e-01; 9.0e-02; 1.1e-01; 4.1e-02 | yes |
| 0.00707 | 15 | 1.2e-04, 1.3e-06, 8.8e-09; 5.4e-04, 2.3e-05, 2.0e-06; 2.3e-05, 1.4e-07, 4.6e-09; 1.3e-03, 2.9e-04, 1.5e-04 | 8.7e-05; 2.2e-03; 9.5e-05; 1.9e-02 | 5.92e+00 / 5.92e+00; 7.15e+00 / 7.27e+00; 2.32e+00 / 2.32e+00; 1.34e+01 / — | 2.3e-03; 6.0e-02; 3.0e-03; — | yes |
| 0.00707 | 60 | 5.3e-06, 8.0e-08, 2.6e-09; 1.1e-04, 2.8e-07, 5.4e-10; 9.2e-06, 1.2e-07, 6.2e-09; 6.2e-05, 1.8e-07, 1.7e-09 | 8.0e-05; 5.7e-05; 1.1e-04; 4.4e-05 | 1.58e+00 / 1.59e+00; 4.16e+00 / 4.16e+00; 1.86e+00 / 1.86e+00; 3.95e+00 / 3.95e+00 | 3.0e-03; 1.8e-03; 3.7e-03; 8.3e-04 | yes |
| 0.00707 | 240 | 1.9e-05, 4.0e-07, 1.9e-08; 3.7e-07, 1.1e-09, 2.5e-12; 8.9e-06, 4.8e-08, 3.3e-09; 6.8e-06, 1.2e-08, 1.3e-09 | 3.7e-04; 7.9e-05; 7.7e-05; 1.9e-05 | 7.58e-01 / 7.59e-01; 3.55e-01 / 3.55e-01; 8.67e-01 / 8.67e-01; 1.37e+00 / 1.37e+00 | 1.3e-02; 2.8e-03; 2.8e-03; 6.8e-04 | yes |
| ∞ | 15 | 3.4e-05, 7.2e-09, 6.5e-10; 1.8e-04, 1.5e-06, 5.4e-08; 4.4e-06, 2.0e-09, 1.9e-10; 6.3e-04, 1.4e-04, 7.0e-05 | 1.1e-05; 3.9e-04; 1.1e-05; 1.3e-02 | 5.92e+00 / 5.92e+00; 7.26e+00 / 7.27e+00; 2.32e+00 / 2.32e+00; 1.48e+01 / — | 4.2e-04; 1.2e-02; 3.9e-04; — | yes |
| ∞ | 60 | 1.5e-06, 2.1e-09, 2.4e-10; 2.1e-05, 4.5e-09, 6.2e-10; 2.5e-06, 4.8e-09, 6.2e-10; 1.3e-05, 3.5e-09, 4.7e-10 | 2.7e-05; 7.4e-06; 2.9e-05; 1.8e-05 | 1.58e+00 / 1.59e+00; 4.16e+00 / 4.16e+00; 1.86e+00 / 1.86e+00; 3.95e+00 / 3.95e+00 | 1.0e-03; 2.5e-04; 1.0e-03; 6.2e-04 | yes |
| ∞ | 240 | 5.8e-06, 1.6e-08, 2.2e-09; 7.2e-08, 1.5e-12, 7.1e-14; 1.8e-06, 3.0e-09, 4.8e-10; 1.1e-06, 9.9e-10, 1.7e-10 | 1.2e-04; 3.1e-05; 2.4e-05; 4.5e-06 | 7.59e-01 / 7.59e-01; 3.55e-01 / 3.55e-01; 8.67e-01 / 8.67e-01; 1.37e+00 / 1.37e+00 | 4.5e-03; 1.1e-03; 8.9e-04; 1.7e-04 | yes |

### The fixed prior at the sampled weights: stability, splitting, and the floor as the exact bias

Per seed and size at s = 6000: the D‴.5 margin of the fixed prior's Θ*(w) (its class), the largest G⁻¹-diameter of a class at the end with the sub-cluster sizes of that class (single linkage at 0.05) for the fixed prior and the diameters for the adaptive priors; the observed ‖f̂_m − f̂_M‖² of the fixed prior against the exact bias at the sampled weights, ‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖² (Proposition D‴.24(i) without the linearisation), and against the delta method's ‖F(w^m − w^M)‖² (Proposition D‴.9's law per sample).

| seed | m | counts | D‴.5 margin of Θ*(w) | end diameter, fixed (sub-clusters) | adaptive η = σ₀, ∞ | ‖f̂_m − f̂_M‖², fixed | exact bias (observed/exact) | delta method (observed/delta) |
|---|---|---|---|---|---|---|---|---|
| 0 | 15 | (6, 5, 4) | -3.48e-03 (class 1) saddle | 1.11 (5, 1) | 0.0e+00, 0.0e+00 | 1.653e-03 | 1.896e-03 (0.8718) | 2.211e-03 (0.75) |
| 0 | 60 | (19, 21, 20) | +5.61e-03 (class 2) | 0.00 (19,) | 0.0e+00, 0.0e+00 | 1.062e-04 | 1.062e-04 (1.0000) | 1.231e-04 (0.86) |
| 0 | 240 | (71, 76, 93) | +3.78e-03 (class 3) | 0.00 (71,) | 0.0e+00, 0.0e+00 | 8.166e-05 | 8.166e-05 (1.0000) | 4.757e-05 (1.72) |
| 1 | 15 | (6, 3, 6) | -1.08e-03 (class 1) saddle | 0.64 (3, 1, 1, 1) | 0.0e+00, 0.0e+00 | 2.404e-03 | 2.471e-03 (0.9729) | 1.835e-03 (1.31) |
| 1 | 60 | (24, 15, 21) | -2.16e-03 (class 1) saddle | 0.73 (10, 6, 2, 2, 2, 1, 1) | 0.0e+00, 0.0e+00 | 9.913e-04 | 1.116e-03 (0.8880) | 1.463e-03 (0.68) |
| 1 | 240 | (77, 76, 87) | +5.06e-03 (class 3) | 0.00 (77,) | 0.0e+00, 0.0e+00 | 4.848e-06 | 4.848e-06 (1.0000) | 5.207e-06 (0.93) |
| 2 | 15 | (5, 4, 6) | +3.95e-03 (class 3) | 0.00 (5,) | 0.0e+00, 0.0e+00 | 3.000e-04 | 3.000e-04 (1.0000) | 2.785e-04 (1.08) |
| 2 | 60 | (19, 18, 23) | +4.17e-03 (class 3) | 0.00 (19,) | 0.0e+00, 0.0e+00 | 1.575e-04 | 1.575e-04 (1.0000) | 1.702e-04 (0.93) |
| 2 | 240 | (76, 81, 83) | +5.95e-03 (class 3) | 0.00 (76,) | 0.0e+00, 0.0e+00 | 7.275e-05 | 7.275e-05 (1.0000) | 6.100e-05 (1.19) |
| 3 | 15 | (8, 4, 3) | -6.14e-03 (class 1) saddle | 1.50 (2, 2, 2, 1, 1) | 0.0e+00, 0.0e+00 | 5.516e-03 | 5.460e-03 (1.0103) | 9.225e-03 (0.60) |
| 3 | 60 | (25, 17, 18) | -3.98e-03 (class 1) saddle | 0.95 (15, 8, 1, 1) | 0.0e+00, 0.0e+00 | 7.126e-04 | 8.786e-04 (0.8112) | 1.708e-03 (0.42) |
| 3 | 240 | (82, 79, 79) | +5.13e-03 (class 1) | 0.00 (82,) | 0.0e+00, 0.0e+00 | 7.502e-05 | 7.502e-05 (1.0000) | 7.574e-05 (0.99) |

5 of 12 fixed-prior m-systems sit at weights where Θ*(w) is a saddle; their observed floors are 0.811–1.010 of the replicated point's (they have split), against 1.0000–1.0000 at the stable ones. The delta method's prediction misses the observed floor by the factors 0.86–1.72 at the stable samples (0.42–1.31 at the saddles).

### The parameter level: the coupled distance D and W²

Per seed and size at s = 6000: Proposition D‴.5's coupled distance D of the shared particles, observed and at the replicated stationary points with the same classes (the fixed prior's Θ*(w); the adaptive prior's exact representations ϑ^e(w), —: one of the two weight vectors outside W); W²(μ_m, μ_M) by an optimal plan split into the transport between different classes (the cells of the exact atoms) and within classes, and ½‖w^m − w^M‖₁ (the cross-class mass of every optimal plan here).

| seed | m | ½‖w^m − w^M‖₁ | D fixed: observed / Θ*(w) | D adaptive η = σ₀, ∞: observed / ϑ^e(w) | D adaptive (∞) / fixed | W² fixed: cross + within | W² adaptive η = σ₀: cross + within | W² adaptive η = ∞: cross + within |
|---|---|---|---|---|---|---|---|---|
| 0 | 15 | 0.1078 | 0.709 / 0.749 | 1.563, 1.565 / 1.565 | 2.21 | 3.900 + 0.576 | 4.570 + 1.349 | 4.571 + 1.351 |
| 0 | 60 | 0.0412 | 0.055 / 0.055 | 0.094, 0.094 / 0.094 | 1.72 | 1.517 + 0.053 | 1.494 + 0.091 | 1.494 + 0.091 |
| 0 | 240 | 0.0140 | 0.015 / 0.015 | 0.072, 0.072 / 0.072 | 4.81 | 0.666 + 0.015 | 0.687 + 0.071 | 0.687 + 0.072 |
| 1 | 15 | 0.1255 | 0.624 / 0.624 | 2.389, 2.497 / 2.512 | 4.00 | 4.967 + 0.546 | 4.895 + 2.254 | 4.895 + 2.361 |
| 1 | 60 | 0.0804 | 0.428 / 0.433 | 0.931, 0.931 / 0.931 | 2.17 | 3.345 + 0.353 | 3.355 + 0.801 | 3.356 + 0.801 |
| 1 | 240 | 0.0088 | 0.002 / 0.002 | 0.004, 0.004 / 0.004 | 1.78 | 0.350 + 0.002 | 0.350 + 0.004 | 0.351 + 0.004 |
| 2 | 15 | 0.0608 | 0.133 / 0.133 | 0.258, 0.258 / 0.258 | 1.94 | 2.094 + 0.121 | 2.077 + 0.242 | 2.076 + 0.242 |
| 2 | 60 | 0.0441 | 0.070 / 0.070 | 0.131, 0.131 / 0.131 | 1.89 | 1.736 + 0.065 | 1.733 + 0.124 | 1.733 + 0.124 |
| 2 | 240 | 0.0167 | 0.018 / 0.018 | 0.052, 0.052 / 0.052 | 2.84 | 0.816 + 0.018 | 0.815 + 0.052 | 0.815 + 0.052 |
| 3 | 15 | 0.2039 | 1.974 / 2.091 | 6.383, 7.800 / — | 3.95 | 6.565 + 1.422 | 7.747 + 5.616 | 7.749 + 7.031 |
| 3 | 60 | 0.0873 | 0.448 / 0.452 | 0.775, 0.775 / 0.775 | 1.73 | 3.274 + 0.357 | 3.298 + 0.649 | 3.298 + 0.649 |
| 3 | 240 | 0.0316 | 0.032 / 0.032 | 0.062, 0.062 / 0.062 | 1.94 | 1.335 + 0.031 | 1.307 + 0.060 | 1.307 + 0.060 |

The fixed prior's D is 0.944–1.000 of the replicated Θ*(w)'s; the adaptive prior's is 0.951–1.000 of the exact representations' (both weight vectors in W); the adaptive/fixed ratio of D is 1.72–4.81. The cross-class transport is 0.82–0.99 of W² under the fixed prior and 0.67–0.99 under the adaptive ones where both weight vectors are in W (0.52–0.58 otherwise, the parameters running off); there the adaptive (η = ∞)/fixed ratio of W² is 1.00–1.32, of its within-class part 1.71–4.84.

## C. The design table at the misplaced prior (centres +0.15), 3000 slow units

The fit error (the bias; 0 at an exact representation) and ‖f̂_m − f̂_M‖² (the floor) at the end, means over the seeds (min–max); the OT leak toward the adaptive prior with the settled gate of Theorem D‴.23.

| leak | prior | m | fit error at the end | ‖f̂_m − f̂_M‖² at the end |
|---|---|---|---|---|
| nearest atom | fixed | 15 | 6.76e-02 (3.0e-02–1.2e-01) | 2.50e-03 (3.0e-04–5.6e-03) |
| nearest atom | fixed | 60 | 3.12e-02 (1.6e-02–4.5e-02) | 5.10e-04 (1.0e-04–1.0e-03) |
| nearest atom | fixed | 240 | 2.01e-02 (1.2e-02–3.7e-02) | 5.89e-05 (4.8e-06–8.2e-05) |
| nearest atom | adaptive (η = σ₀) | 15 | 8.93e-03 (5.8e-04–2.7e-02) | 7.96e-05 (1.4e-07–2.9e-04) |
| nearest atom | adaptive (η = σ₀) | 60 | 5.53e-04 (2.8e-04–7.3e-04) | 1.64e-07 (8.0e-08–2.8e-07) |
| nearest atom | adaptive (η = σ₀) | 240 | 5.97e-04 (1.4e-04–1.6e-03) | 1.14e-07 (1.1e-09–3.9e-07) |
| OT | fixed | 15 | 9.78e-03 (9.8e-03–9.8e-03) | 5.08e-22 (2.3e-22–8.0e-22) |
| OT | fixed | 60 | 9.78e-03 (9.8e-03–9.8e-03) | 2.47e-22 (2.4e-23–6.2e-22) |
| OT | fixed | 240 | 9.78e-03 (9.8e-03–9.8e-03) | 1.16e-21 (4.1e-23–2.4e-21) |
| OT | adaptive (η = σ₀, settled gate) | 15 | 1.59e-04 (1.5e-04–1.6e-04) | 1.05e-11 (2.6e-12–2.8e-11) |
| OT | adaptive (η = σ₀, settled gate) | 60 | 1.62e-04 (1.6e-04–1.6e-04) | 1.37e-11 (9.6e-15–3.2e-11) |
| OT | adaptive (η = σ₀, settled gate) | 240 | 1.65e-04 (1.6e-04–1.7e-04) | 1.44e-11 (2.3e-13–4.2e-11) |

The coupling level and the measure level at the end (means over the seeds, min–max): the coupled distance D of the shared particles and W²(μ_m, μ_M), with the cross-class part's share (classes: the cells of the exact atoms); the particles whose end cell differs from their initial class.

| leak, prior | m | D | W² | cross-class share of W² (—: W² below 10⁻¹²) | relabelled particles (all seeds) |
|---|---|---|---|---|---|
| nearest atom, fixed | 15 | 0.862 (0.133–1.988) | 4.98e+00 (2.2e+00–7.8e+00) | 0.82–0.95 | 0 |
| nearest atom, fixed | 60 | 0.251 (0.055–0.452) | 2.62e+00 (1.6e+00–3.6e+00) | 0.91–0.97 | 0 |
| nearest atom, fixed | 240 | 0.017 (0.002–0.032) | 8.10e-01 (3.5e-01–1.4e+00) | 0.98–0.99 | 0 |
| nearest atom, adaptive | 15 | 2.292 (0.256–5.259) | 6.83e+00 (2.3e+00–1.2e+01) | 0.63–0.90 | 0 |
| nearest atom, adaptive | 60 | 0.479 (0.093–0.921) | 2.88e+00 (1.6e+00–4.1e+00) | 0.81–0.94 | 0 |
| nearest atom, adaptive | 240 | 0.046 (0.004–0.066) | 8.35e-01 (3.5e-01–1.4e+00) | 0.91–0.99 | 0 |
| OT, fixed | 15 | 6.033 (2.619–9.448) | 8.92e-16 (0.0e+00–3.5e-15) | — | 7 |
| OT, fixed | 60 | 3.484 (2.057–5.402) | 0.00e+00 (0.0e+00) | — | 14 |
| OT, fixed | 240 | 0.997 (0.187–1.473) | 0.00e+00 (0.0e+00) | — | 26 |
| OT, adaptive | 15 | 6.019 (2.615–9.423) | 5.54e-09 (1.3e-09–1.5e-08) | 0.00 | 7 |
| OT, adaptive | 60 | 3.480 (2.053–5.394) | 6.94e-09 (8.8e-12–1.6e-08) | 0.00 | 14 |
| OT, adaptive | 240 | 0.994 (0.186–1.467) | 7.32e-09 (1.2e-10–2.1e-08) | 0.00 | 26 |

The fixed misplaced prior's Θ*(w) at the sampled weights (Part C's nearest-atom fixed design): the D‴.5 margin (its class) and the end diameter (sub-clusters).

| m | seed 0 | seed 1 | seed 2 | seed 3 |
|---|---|---|---|---|
| 15 | -2.87e-03 (class 1) saddle; 1.18 (5, 1) | -3.79e-04 (class 1) saddle; 0.25 (4, 1, 1) | +2.78e-03 (class 3); 0.00 (6,) | -5.76e-03 (class 1) saddle; 1.53 (5, 1, 1, 1) |
| 60 | +5.55e-03 (class 2); 0.00 (21,) | -1.52e-03 (class 1) saddle; 0.86 (18, 2, 1, 1, 1, 1) | +3.05e-03 (class 3); 0.00 (23,) | -3.41e-03 (class 1) saddle; 1.29 (20, 3, 1, 1) |
| 240 | +2.71e-03 (class 3); 0.00 (93,) | +3.90e-03 (class 3); 0.00 (87,) | +4.74e-03 (class 3); 0.00 (76,) | +5.74e-03 (class 1); 0.00 (82,) |

The ranking depends on the horizon: the OT leak toward the adaptive prior equals the fixed one until its settled gate opens (s = 1211–1304). The fit errors at s = 1250 and 2000 (means over the seeds, min–max):

| leak, prior | m | fit error at s = 1250 | at s = 2000 |
|---|---|---|---|
| nearest atom, fixed | 15 | 6.76e-02 (3.0e-02–1.2e-01) | 6.76e-02 (3.0e-02–1.2e-01) |
| nearest atom, fixed | 60 | 3.14e-02 (1.6e-02–4.6e-02) | 3.13e-02 (1.6e-02–4.5e-02) |
| nearest atom, fixed | 240 | 2.01e-02 (1.2e-02–3.7e-02) | 2.01e-02 (1.2e-02–3.7e-02) |
| nearest atom, adaptive | 15 | 1.75e-02 (3.3e-03–4.0e-02) | 1.24e-02 (1.5e-03–3.2e-02) |
| nearest atom, adaptive | 60 | 3.77e-03 (2.0e-03–6.0e-03) | 1.60e-03 (8.1e-04–2.3e-03) |
| nearest atom, adaptive | 240 | 2.77e-03 (8.3e-04–6.8e-03) | 1.39e-03 (3.6e-04–3.6e-03) |
| OT, fixed | 15 | 9.78e-03 (9.8e-03–9.8e-03) | 9.78e-03 (9.8e-03–9.8e-03) |
| OT, fixed | 60 | 9.78e-03 (9.8e-03–9.8e-03) | 9.78e-03 (9.8e-03–9.8e-03) |
| OT, fixed | 240 | 9.78e-03 (9.8e-03–9.8e-03) | 9.78e-03 (9.8e-03–9.8e-03) |
| OT, adaptive | 15 | 9.17e-03 (7.9e-03–9.8e-03) | 4.02e-04 (3.9e-04–4.2e-04) |
| OT, adaptive | 60 | 9.59e-03 (9.0e-03–9.8e-03) | 4.11e-04 (4.0e-04–4.2e-04) |
| OT, adaptive | 240 | 9.71e-03 (9.5e-03–9.8e-03) | 4.20e-04 (4.0e-04–4.3e-04) |

## D. The real closed loop (k_e = 25, dt = 0.001, 3000 time units), matching prior, nearest-atom leak, seed 0's sampled weights

The rms tracking error over the five periods ending at t = 31, 1000, 2000, 3000, against (1/k_e)‖f − f̂‖_ρ of the averaged loop from the same start at s = t/k_e (the slaved error of Theorem J at leading order); from the exact representation of w with the adaptive prior at its atoms the tracking error stays at the integrator's level.

| m | counts | start, prior | e_rms at those times | (1/k_e)‖f − f̂‖_ρ averaged at s = t/k_e | ratio at the end |
|---|---|---|---|---|---|
| 15 | (6, 5, 4) | Θ*(w), fixed prior | 1.173e-03, 1.180e-03, 1.180e-03, 1.182e-03 | 1.181e-03, 1.181e-03, 1.181e-03, 1.181e-03 | 1.001 |
| 15 | (6, 5, 4) | Θ*(w), adaptive prior (η = ∞) | 1.157e-03, 6.102e-04, 4.062e-04, 3.273e-04 | 1.151e-03, 6.056e-04, 4.093e-04, 3.317e-04 | 0.987 |
| 15 | (6, 5, 4) | ϑ^e(w), adaptive prior at ϑ^e(w) | 9.0e-12, 3.6e-12, 2.1e-12, 1.4e-12 | 6.9e-18, 6.8e-18, 6.8e-18, 6.8e-18 | — |
| 60 | (19, 21, 20) | Θ*(w), fixed prior | 3.427e-04, 3.429e-04, 3.433e-04, 3.441e-04 | 3.456e-04, 3.456e-04, 3.456e-04, 3.456e-04 | 0.996 |
| 60 | (19, 21, 20) | Θ*(w), adaptive prior (η = ∞) | 3.397e-04, 2.294e-04, 1.794e-04, 1.522e-04 | 3.401e-04, 2.312e-04, 1.812e-04, 1.534e-04 | 0.992 |
| 60 | (19, 21, 20) | ϑ^e(w), adaptive prior at ϑ^e(w) | 9.0e-12, 3.7e-12, 2.2e-12, 1.5e-12 | 8.7e-18, 8.5e-18, 8.5e-18, 8.6e-18 | — |

The fixed prior's Θ*(w) at seed 0's sampled weights: m = 15: D‴.5 margin -3.48e-03 (class 1), a saddle: the fixed-prior run sits on a saddle orbit (Proposition D‴.15), held there by the exact replication of its start; m = 60: D‴.5 margin +5.61e-03 (class 2), stable.

## E. The instance's probabilities: representability, the fixed prior's saddles, and the adaptive prior's expected limiting floor

Exact sums of the multinomial probabilities (q = 1/3) over the count vectors of probability above 0 (m = 15), 1e-12 (m = 60), 1e-09 (m = 240); the neglected mass is given. Representable: w ∈ W (an exact representation on the branch of ϑ^e(q)); boundary: some w_k = |A_k|/a_max; with a permutation: outside W, but a permutation of the classes over the bumps is representable (an exact representation exists off the branch). Hoeffding: Σ_kexp(−2m(q_k − |A_k|/a_max)²) bounds P(outside W ∪ boundary). The fixed matching prior's saddle: Θ*(w) not D‴.5-admissible (classes with at least two particles; a vector with an empty class has no K₀-cluster configuration). The adaptive prior's expected limiting floor: Σ p(w)·min_{|a| ≤ a_max} ⟨(f − f̂_w)²⟩_ρ over the vectors outside W (the constrained minimum on the branch of ϑ^e(q), by L-BFGS-B from the bumps, taken as the dynamics' limit — the dynamics from the warm clouds follows the branch; a local minimum, not the infimum, which is 0 where a permutation is representable; 0 on W and on its boundary, the reference taken exact), a lower bound on that limit since the sum is truncated; against the fixed prior's law (1/m − 1/M)·tr(FC̃Fᵀ)/N and its Monte Carlo expectation (Part F).

| m | vectors | neglected mass | P(outside W) | P(boundary) | Hoeffding | share of P(outside W) with a permutation | P(fixed-prior saddle) | P(empty class) | Newton failed | expected limiting floor, adaptive | fixed: law | fixed: Monte Carlo | adaptive/law | adaptive/Monte Carlo |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 15 | 136 | none (all vectors) | 0.3606 | 0.1179 | 1.870 | 0.31 | 0.3687 | 6.9e-03 | 1.5e-05 | 7.25e-03 | 3.99e-03 | 1.11e-02 | 1.818 | 0.655 |
| 60 | 1432 | 4.1e-11 | 0.0738 | 0.0534 | 0.602 | 0.87 | 0.2597 | 6.9e-11 | 0.0e+00 | 2.16e-04 | 9.07e-04 | 1.39e-03 | 0.238 | 0.155 |
| 240 | 4264 | 2.9e-07 | 0.0021 | 0.0012 | 0.036 | 1.00 | 0.1089 | 0.0e+00 | 0.0e+00 | 6.58e-07 | 1.36e-04 | 1.51e-04 | 0.005 | 0.004 |

## F. Monte Carlo: the exact floors and coupled distances against the laws (600 nested samples)

M = 510 labels i.i.d. uniform (seeds 10⁴ + i), the m-system the first m. At each sample the replicated stationary points: the fixed prior's Θ*(w) and the adaptive prior's ϑ^e(w). The ratios to the laws per unit of 1/m − 1/M (Part A, σ₀ = 0.15σ_crit), mean ± standard error: ‖f̂_m − f̂_M‖² of the fixed prior (all samples, and conditionally on both samples D‴.5-admissible), and the coupled distances (fixed: all samples; adaptive: both samples in W).

| m | E‖f̂_m − f̂_M‖²/law, fixed | conditional on D‴.5-admissible | samples | E D/law, fixed | E D/law, adaptive (in W) | samples | E D adaptive / E D fixed (same samples) |
|---|---|---|---|---|---|---|---|
| 15 | 2.776 ± 0.198 | 2.887 ± 0.262 | 595, 365 | 0.682 ± 0.023 | 0.529 ± 0.024 | 316 | 3.65 |
| 60 | 1.532 ± 0.106 | 1.564 ± 0.145 | 600, 414 | 0.925 ± 0.037 | 1.171 ± 0.061 | 533 | 3.30 |
| 240 | 1.114 ± 0.065 | 1.067 ± 0.071 | 600, 520 | 0.952 ± 0.041 | 1.259 ± 0.110 | 596 | 2.93 |

The tails (the largest sample's D/law; in parentheses the adaptive prior's mean without that sample): m = 15: fixed 3.36, adaptive 1.9 (0.525); m = 60: fixed 5.64, adaptive 8.5 (1.157); m = 240: fixed 7.47, adaptive 55.7 (1.167). The adaptive prior's gauge diverges at the boundary of W, so its coupled distance has a far heavier tail, and its mean over these samples is not a reliable estimate of its law.

Numbers are produced by `experiments/exp_adaptive_floor.py`; nothing is edited by hand.
