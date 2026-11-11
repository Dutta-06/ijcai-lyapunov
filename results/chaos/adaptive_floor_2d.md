# Proposition D‴.24 in two dimensions: the floor is the leak's bias, and an adaptive prior removes it (`exp_adaptive_floor_2d.py`)

The planar instance of `exp_ot_2d.py` (four bumps with vector amplitudes, the Lissajous reference, p = 5, G = diag(1.0, 1.0, 0.25, 0.25, 0.25), σ₀ = 0.05, the matching prior = the exact representation with masses 1/4) and the averaged slow loop of `src/core/planar_loop.py` (ρ: 400 times of one period).

## A. The floor as the weight-derivative of the leak's bias (matching prior, no dynamics)

The whitened weighted Gauss–Newton spectrum at the exact representation (h_j/σ₀): 0.004, 0.012, 0.017, 0.024, 0.052, 0.072, 0.165, 0.184, 0.218, 0.253, 0.281, 0.399, 0.414, 0.428, 0.661, 0.844, 1.092, 1.257, 1.312, 1.413.

The fixed prior's weighted stationary points at q ± 0.0001(e_j − e_4) by Newton (largest residual 2.1e-17; every atom in its own cell: True); at q itself the fit error is 1.3e-16 (residual 6.0e-18). F against the bias form Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰: ‖F − F_bias‖/‖F‖ = 6.0e-07. The law per unit of 1/m − 1/M: 3.0962e-01 from F, 3.0962e-01 as the sum of the mode terms, against the bound 1.4380e+00 (law/bound 0.215).

| mode j | h_j/σ₀ | share of the law |
|---|---|---|
| 1 | 0.004 | 0.033 |
| 2 | 0.012 | 0.015 |
| 3 | 0.017 | 0.020 |
| 4 | 0.024 | 0.015 |
| 5 | 0.052 | 0.065 |
| 6 | 0.072 | 0.161 |
| 7 | 0.165 | 0.047 |
| 8 | 0.184 | 0.014 |
| 9 | 0.218 | 0.019 |
| 10 | 0.253 | 0.050 |
| 11 | 0.281 | 0.038 |
| 12 | 0.399 | 0.116 |
| 13 | 0.414 | 0.090 |
| 14 | 0.428 | 0.054 |
| 15 | 0.661 | 0.184 |
| 16 | 0.844 | 0.005 |
| 17 | 1.092 | 0.036 |
| 18 | 1.257 | 0.010 |
| 19 | 1.312 | 0.017 |
| 20 | 1.413 | 0.009 |

Representable weights (the set W): an exact representation at w needs |A_k,c|/w_k < a_max, i.e. w > (0.200, 0.150, 0.175, 0.150).

## B. The nearest-atom leak from the warm clouds: fixed (η = 0) against adaptive priors (6000 slow units)

Midpoint rule dt = 0.25, the cells every 4 steps; η = ∞: the atoms reset to the cell means at every refresh of the cells (every 1 slow unit) and frozen in between — a lagged form of the η = ∞ prior.

‖f̂_m − f̂_M‖² (the function-level floor) and the fit error at s = 300, 1000, 3000, 6000; the sampled counts and their margin min_k(w_k − max_c|A_k,c|/a_max) (representable iff positive); at the end the class means' whitened weighted distance from the exact representation of the sampled weights, and W²(μ_m, μ_M) against the W² between the two systems' exact representations (—: no exact representation); the counts conserved.

| seed | m | counts (margin) | η | ‖f̂_m − f̂_M‖² at those s | fit error at those s | distance from ϑ^e(w) | W² / exact-representation W² | conserved |
|---|---|---|---|---|---|---|---|---|
| 0 | 16 | (4, 2, 6, 4) (-0.025) | 0 | 9.0e-03, 9.0e-03, 9.0e-03, 9.0e-03 | 2.8e-01, 2.8e-01, 2.8e-01, 2.8e-01 | — | 7.762 / — | yes |
| 0 | 16 | (4, 2, 6, 4) (-0.025) | 0.05 | 2.9e-03, 1.5e-03, 7.6e-04, 5.5e-04 | 1.4e-01, 9.3e-02, 6.2e-02, 5.2e-02 | — | 10.630 / — | yes |
| 0 | 16 | (4, 2, 6, 4) (-0.025) | ∞ | 2.0e-03, 9.6e-04, 5.5e-04, 4.2e-04 | 1.1e-01, 7.0e-02, 5.2e-02, 4.6e-02 | — | 12.173 / — | yes |
| 0 | 32 | (7, 4, 12, 9) (-0.025) | 0 | 8.7e-03, 8.7e-03, 8.7e-03, 8.7e-03 | 2.9e-01, 2.9e-01, 2.9e-01, 2.9e-01 | — | 7.285 / — | yes |
| 0 | 32 | (7, 4, 12, 9) (-0.025) | 0.05 | 2.9e-03, 1.5e-03, 6.9e-04, 4.9e-04 | 1.4e-01, 9.7e-02, 6.3e-02, 5.3e-02 | — | 10.058 / — | yes |
| 0 | 32 | (7, 4, 12, 9) (-0.025) | ∞ | 2.0e-03, 9.1e-04, 5.0e-04, 4.0e-04 | 1.2e-01, 7.3e-02, 5.3e-02, 4.7e-02 | — | 11.535 / — | yes |
| 0 | 64 | (14, 13, 21, 16) (+0.019) | 0 | 7.0e-04, 7.0e-04, 7.0e-04, 7.0e-04 | 1.5e-01, 1.5e-01, 1.5e-01, 1.5e-01 | 1.1e+00 | 2.236 / 2.636 | yes |
| 0 | 64 | (14, 13, 21, 16) (+0.019) | 0.05 | 1.2e-04, 4.1e-05, 1.1e-05, 4.3e-06 | 6.4e-02, 3.0e-02, 1.6e-02, 1.1e-02 | 5.2e-01 | 2.251 / 2.636 | yes |
| 0 | 64 | (14, 13, 21, 16) (+0.019) | ∞ | 6.7e-05, 1.7e-05, 4.4e-06, 2.1e-06 | 4.5e-02, 1.9e-02, 1.1e-02, 7.7e-03 | 3.9e-01 | 2.313 / 2.636 | yes |
| 1 | 16 | (3, 6, 1, 6) (-0.112) | 0 | 1.9e-02, 1.9e-02, 1.9e-02, 1.9e-02 | 4.2e-01, 4.2e-01, 4.2e-01, 4.2e-01 | — | 9.168 / — | yes |
| 1 | 16 | (3, 6, 1, 6) (-0.112) | 0.05 | 1.3e-02, 1.2e-02, 1.2e-02, 1.2e-02 | 3.0e-01, 2.8e-01, 2.6e-01, 2.6e-01 | — | 14.963 / — | yes |
| 1 | 16 | (3, 6, 1, 6) (-0.112) | ∞ | 1.2e-02, 1.3e-02, 1.2e-02, 1.2e-02 | 2.9e-01, 2.7e-01, 2.6e-01, 2.5e-01 | — | 17.551 / — | yes |
| 1 | 32 | (7, 10, 4, 11) (-0.050) | 0 | 5.4e-03, 5.4e-03, 5.4e-03, 5.4e-03 | 2.7e-01, 2.7e-01, 2.7e-01, 2.7e-01 | — | 4.741 / — | yes |
| 1 | 32 | (7, 10, 4, 11) (-0.050) | 0.05 | 3.0e-03, 2.5e-03, 2.3e-03, 2.0e-03 | 1.7e-01, 1.4e-01, 1.1e-01, 1.0e-01 | — | 8.521 / — | yes |
| 1 | 32 | (7, 10, 4, 11) (-0.050) | ∞ | 2.7e-03, 2.4e-03, 2.0e-03, 1.7e-03 | 1.5e-01, 1.2e-01, 1.0e-01, 9.6e-02 | — | 10.992 / — | yes |
| 1 | 64 | (15, 21, 10, 18) (-0.019) | 0 | 3.2e-03, 3.2e-03, 3.2e-03, 3.2e-03 | 2.0e-01, 2.0e-01, 2.0e-01, 2.0e-01 | — | 4.881 / — | yes |
| 1 | 64 | (15, 21, 10, 18) (-0.019) | 0.05 | 1.1e-03, 8.4e-04, 6.4e-04, 5.1e-04 | 1.1e-01, 8.5e-02, 6.1e-02, 5.1e-02 | — | 6.879 / — | yes |
| 1 | 64 | (15, 21, 10, 18) (-0.019) | ∞ | 9.4e-04, 7.2e-04, 5.1e-04, 4.1e-04 | 9.8e-02, 7.0e-02, 5.2e-02, 4.6e-02 | — | 7.946 / — | yes |
| 2 | 16 | (4, 5, 2, 5) (-0.050) | 0 | 1.4e-02, 1.4e-02, 1.4e-02, 1.4e-02 | 2.5e-01, 2.5e-01, 2.5e-01, 2.5e-01 | — | 7.172 / — | yes |
| 2 | 16 | (4, 5, 2, 5) (-0.050) | 0.05 | 6.3e-03, 3.8e-03, 2.5e-03, 2.1e-03 | 1.7e-01, 1.4e-01, 1.1e-01, 1.0e-01 | — | 12.076 / — | yes |
| 2 | 16 | (4, 5, 2, 5) (-0.050) | ∞ | 4.7e-03, 2.9e-03, 2.1e-03, 1.7e-03 | 1.5e-01, 1.2e-01, 1.0e-01, 9.6e-02 | — | 15.262 / — | yes |
| 2 | 32 | (7, 10, 9, 6) (+0.019) | 0 | 1.7e-03, 1.7e-03, 1.7e-03, 1.7e-03 | 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01 | 1.2e+00 | 2.980 / 6.394 | yes |
| 2 | 32 | (7, 10, 9, 6) (+0.019) | 0.05 | 2.6e-04, 9.6e-05, 3.9e-05, 1.8e-05 | 5.9e-02, 3.2e-02, 1.9e-02, 1.2e-02 | 5.2e-01 | 3.564 / 6.394 | yes |
| 2 | 32 | (7, 10, 9, 6) (+0.019) | ∞ | 1.4e-04, 5.9e-05, 1.8e-05, 1.0e-05 | 4.2e-02, 2.4e-02, 1.2e-02, 7.6e-03 | 3.8e-01 | 3.725 / 6.394 | yes |
| 2 | 64 | (13, 18, 16, 17) (+0.003) | 0 | 3.2e-04, 3.2e-04, 3.2e-04, 3.2e-04 | 9.0e-02, 9.0e-02, 9.0e-02, 9.0e-02 | 2.4e+00 | 0.988 / 0.981 | yes |
| 2 | 64 | (13, 18, 16, 17) (+0.003) | 0.05 | 5.3e-05, 9.3e-06, 1.3e-06, 1.0e-07 | 3.8e-02, 3.1e-02, 2.4e-02, 1.9e-02 | 1.9e+00 | 0.990 / 0.981 | yes |
| 2 | 64 | (13, 18, 16, 17) (+0.003) | ∞ | 2.2e-05, 3.2e-06, 1.1e-07, 1.2e-09 | 3.4e-02, 2.7e-02, 1.9e-02, 1.5e-02 | 1.7e+00 | 0.986 / 0.981 | yes |
| 3 | 16 | (7, 4, 2, 3) (-0.050) | 0 | 2.1e-02, 2.1e-02, 2.1e-02, 2.1e-02 | 3.2e-01, 3.2e-01, 3.2e-01, 3.2e-01 | — | 11.476 / — | yes |
| 3 | 16 | (7, 4, 2, 3) (-0.050) | 0.05 | 7.2e-03, 4.2e-03, 2.9e-03, 2.2e-03 | 1.8e-01, 1.5e-01, 1.2e-01, 1.1e-01 | — | 15.755 / — | yes |
| 3 | 16 | (7, 4, 2, 3) (-0.050) | ∞ | 5.3e-03, 3.3e-03, 2.2e-03, 1.7e-03 | 1.6e-01, 1.3e-01, 1.1e-01, 9.5e-02 | — | 18.634 / — | yes |
| 3 | 32 | (10, 8, 8, 6) (+0.038) | 0 | 2.2e-03, 2.2e-03, 2.2e-03, 2.2e-03 | 1.3e-01, 1.3e-01, 1.3e-01, 1.3e-01 | 1.0e+00 | 4.287 / 4.376 | yes |
| 3 | 32 | (10, 8, 8, 6) (+0.038) | 0.05 | 3.4e-04, 1.9e-04, 1.0e-04, 4.1e-05 | 4.7e-02, 3.5e-02, 2.4e-02, 1.5e-02 | 3.4e-01 | 4.327 / 4.376 | yes |
| 3 | 32 | (10, 8, 8, 6) (+0.038) | ∞ | 2.4e-04, 1.4e-04, 4.3e-05, 5.4e-06 | 4.0e-02, 2.9e-02, 1.5e-02, 5.5e-03 | 1.2e-01 | 4.337 / 4.376 | yes |
| 3 | 64 | (18, 18, 15, 13) (+0.053) | 0 | 1.3e-03, 1.3e-03, 1.3e-03, 1.3e-03 | 9.7e-02, 9.7e-02, 9.7e-02, 9.7e-02 | 6.9e-01 | 3.196 / 3.412 | yes |
| 3 | 64 | (18, 18, 15, 13) (+0.053) | 0.05 | 2.3e-04, 8.0e-05, 3.2e-05, 1.3e-05 | 3.4e-02, 2.3e-02, 1.5e-02, 8.8e-03 | 2.1e-01 | 3.341 / 3.412 | yes |
| 3 | 64 | (18, 18, 15, 13) (+0.053) | ∞ | 1.3e-04, 4.7e-05, 1.4e-05, 2.3e-06 | 2.7e-02, 1.8e-02, 9.0e-03, 3.6e-03 | 8.7e-02 | 3.379 / 3.412 | yes |
| 4 | 16 | (2, 3, 5, 6) (-0.075) | 0 | 1.5e-02, 1.5e-02, 1.5e-02, 1.5e-02 | 3.1e-01, 3.1e-01, 3.1e-01, 3.1e-01 | — | 10.080 / — | yes |
| 4 | 16 | (2, 3, 5, 6) (-0.075) | 0.05 | 4.9e-03, 2.8e-03, 2.0e-03, 1.9e-03 | 1.6e-01, 1.3e-01, 1.1e-01, 1.1e-01 | — | 14.383 / — | yes |
| 4 | 16 | (2, 3, 5, 6) (-0.075) | ∞ | 3.6e-03, 2.2e-03, 1.9e-03, 1.8e-03 | 1.4e-01, 1.2e-01, 1.1e-01, 1.0e-01 | — | 15.905 / — | yes |
| 4 | 32 | (5, 7, 7, 13) (-0.044) | 0 | 6.6e-03, 6.6e-03, 6.6e-03, 6.6e-03 | 2.6e-01, 2.6e-01, 2.6e-01, 2.6e-01 | — | 6.755 / — | yes |
| 4 | 32 | (5, 7, 7, 13) (-0.044) | 0.05 | 1.3e-03, 8.1e-04, 5.6e-04, 4.8e-04 | 1.0e-01, 8.2e-02, 6.5e-02, 5.8e-02 | — | 7.571 / — | yes |
| 4 | 32 | (5, 7, 7, 13) (-0.044) | ∞ | 9.8e-04, 6.5e-04, 4.8e-04, 4.5e-04 | 9.1e-02, 7.1e-02, 5.8e-02, 5.4e-02 | — | 8.326 / — | yes |
| 4 | 64 | (12, 15, 18, 19) (-0.013) | 0 | 3.0e-03, 3.0e-03, 3.0e-03, 3.0e-03 | 1.3e-01, 1.3e-01, 1.3e-01, 1.3e-01 | — | 4.185 / — | yes |
| 4 | 64 | (12, 15, 18, 19) (-0.013) | 0.05 | 8.3e-04, 2.5e-04, 1.1e-04, 8.3e-05 | 5.8e-02, 4.4e-02, 3.4e-02, 2.9e-02 | — | 4.855 / — | yes |
| 4 | 64 | (12, 15, 18, 19) (-0.013) | ∞ | 4.5e-04, 1.4e-04, 8.3e-05, 6.8e-05 | 5.0e-02, 3.8e-02, 2.9e-02, 2.4e-02 | — | 5.102 / — | yes |
| 5 | 16 | (4, 5, 4, 3) (+0.038) | 0 | 4.2e-03, 4.2e-03, 4.2e-03, 4.2e-03 | 1.2e-01, 1.2e-01, 1.2e-01, 1.2e-01 | 8.1e-01 | 4.660 / 5.404 | yes |
| 5 | 16 | (4, 5, 4, 3) (+0.038) | 0.05 | 7.1e-04, 2.1e-04, 4.4e-05, 9.4e-06 | 4.3e-02, 2.6e-02, 1.4e-02, 6.2e-03 | 1.3e-01 | 5.271 / 5.404 | yes |
| 5 | 16 | (4, 5, 4, 3) (+0.038) | ∞ | 3.8e-04, 8.7e-05, 1.0e-05, 8.6e-07 | 3.2e-02, 1.9e-02, 6.4e-03, 1.7e-03 | 3.9e-02 | 5.362 / 5.404 | yes |
| 5 | 32 | (10, 8, 5, 9) (-0.019) | 0 | 3.0e-03, 3.0e-03, 3.0e-03, 3.0e-03 | 1.8e-01, 1.8e-01, 1.8e-01, 1.8e-01 | — | 3.480 / — | yes |
| 5 | 32 | (10, 8, 5, 9) (-0.019) | 0.05 | 1.4e-03, 9.5e-04, 6.7e-04, 5.1e-04 | 1.2e-01, 8.8e-02, 6.4e-02, 5.3e-02 | — | 5.537 / — | yes |
| 5 | 32 | (10, 8, 5, 9) (-0.019) | ∞ | 1.1e-03, 7.7e-04, 5.2e-04, 4.0e-04 | 1.0e-01, 7.2e-02, 5.4e-02, 4.6e-02 | — | 6.654 / — | yes |
| 5 | 64 | (18, 15, 13, 18) (+0.028) | 0 | 3.7e-04, 3.7e-04, 3.7e-04, 3.7e-04 | 9.9e-02, 9.9e-02, 9.9e-02, 9.9e-02 | 9.1e-01 | 1.652 / 1.770 | yes |
| 5 | 64 | (18, 15, 13, 18) (+0.028) | 0.05 | 1.0e-04, 4.4e-05, 1.4e-05, 3.5e-06 | 5.5e-02, 3.2e-02, 1.5e-02, 7.9e-03 | 1.8e-01 | 1.758 / 1.770 | yes |
| 5 | 64 | (18, 15, 13, 18) (+0.028) | ∞ | 6.6e-05, 2.3e-05, 3.7e-06, 2.7e-07 | 4.3e-02, 2.1e-02, 8.1e-03, 3.0e-03 | 7.0e-02 | 1.775 / 1.770 | yes |

The reference systems (M = 128): seed 0: (27, 30, 39, 32) (+0.011); seed 1: (26, 36, 27, 39) (+0.003); seed 2: (26, 37, 34, 31) (+0.003); seed 3: (32, 33, 34, 29) (+0.050); seed 4: (29, 33, 28, 38) (+0.027); seed 5: (34, 32, 28, 34) (+0.044).

Reading (at s = 6000; a system counts as representable when its sample and its reference's are). The fixed prior's ‖f̂_m − f̂_M‖²: 3.2e-04–4.2e-03 at the representable systems, 3.0e-03–2.1e-02 at the others. The adaptive prior's: η = 0.05, representable (7 systems): 1.0e-07–4.1e-05, 0.000–0.019 of the fixed prior's, falling between the last two records in 7 of 7; η = 0.05, not representable (11 systems): 8.3e-05–1.2e-02, 0.027–0.641 of the fixed prior's, falling between the last two records in 11 of 11; η = ∞, representable (7 systems): 1.2e-09–1.0e-05, 0.000–0.006 of the fixed prior's, falling between the last two records in 7 of 7; η = ∞, not representable (11 systems): 6.8e-05–1.2e-02, 0.023–0.634 of the fixed prior's, falling between the last two records in 11 of 11.

### The fixed prior at the sampled weights: stability and the floor as the exact bias

Per seed and size: the fixed prior's replicated stationary point Θ*(w) at the sampled counts, Proposition D‴.5's splitting margin min_k(σ₀ − λ_max(G^{1/2}R_kG^{1/2})) over the classes with at least two particles (its class), the largest eigenvalue of the weighted K-atom drift's Jacobian (negative: the collective block is stable); Θ*(w) is a strict local minimiser (D‴.5-admissible) iff the margin is positive and the eigenvalue negative. For the m-systems: the fixed prior's observed ‖f̂_m − f̂_M‖² at the end against the exact bias at the sampled weights, ‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖² (ratio 1 for a system that stays at its replicated point).

| seed | m | counts | splitting margin (class) | collective: largest eigenvalue | observed / exact bias |
|---|---|---|---|---|---|
| 0 | 16 | (4, 2, 6, 4) | +1.55e-02 (class 3) | -2.59e-02 | 9.047e-03 / 9.047e-03 (1.0000) |
| 0 | 32 | (7, 4, 12, 9) | +1.55e-02 (class 3) | -2.59e-02 | 8.683e-03 / 8.683e-03 (1.0000) |
| 0 | 64 | (14, 13, 21, 16) | +2.65e-02 (class 3) | -3.37e-02 | 6.998e-04 / 6.998e-04 (1.0000) |
| 0 | 128 | (27, 30, 39, 32) | +3.31e-02 (class 3) | -3.92e-02 | — |
| 1 | 16 | (3, 6, 1, 6) | +2.77e-02 (class 2) | -4.12e-02 | 1.910e-02 / 1.910e-02 (1.0000) |
| 1 | 32 | (7, 10, 4, 11) | +3.63e-02 (class 4) | -4.37e-02 | 5.429e-03 / 5.429e-03 (1.0000) |
| 1 | 64 | (15, 21, 10, 18) | +3.39e-02 (class 2) | -4.42e-02 | 3.152e-03 / 3.152e-03 (1.0000) |
| 1 | 128 | (26, 36, 27, 39) | +4.17e-02 (class 4) | -4.66e-02 | — |
| 2 | 16 | (4, 5, 2, 5) | +3.66e-02 (class 2) | -4.53e-02 | 1.405e-02 / 1.405e-02 (1.0000) |
| 2 | 32 | (7, 10, 9, 6) | +3.62e-02 (class 2) | -4.50e-02 | 1.701e-03 / 1.701e-03 (1.0000) |
| 2 | 64 | (13, 18, 16, 17) | +4.27e-02 (class 2) | -4.83e-02 | 3.246e-04 / 3.246e-04 (1.0000) |
| 2 | 128 | (26, 37, 34, 31) | +4.09e-02 (class 2) | -4.77e-02 | — |
| 3 | 16 | (7, 4, 2, 3) | +3.18e-02 (class 1) | -4.56e-02 | 2.070e-02 / 2.070e-02 (1.0000) |
| 3 | 32 | (10, 8, 8, 6) | +4.15e-02 (class 1) | -4.81e-02 | 2.201e-03 / 2.201e-03 (1.0000) |
| 3 | 64 | (18, 18, 15, 13) | +4.26e-02 (class 2) | -4.81e-02 | 1.311e-03 / 1.311e-03 (1.0000) |
| 3 | 128 | (32, 33, 34, 29) | +4.50e-02 (class 3) | -4.86e-02 | — |
| 4 | 16 | (2, 3, 5, 6) | +3.09e-02 (class 3) | -3.72e-02 | 1.535e-02 / 1.535e-02 (1.0000) |
| 4 | 32 | (5, 7, 7, 13) | +3.06e-02 (class 4) | -3.97e-02 | 6.592e-03 / 6.592e-03 (1.0000) |
| 4 | 64 | (12, 15, 18, 19) | +4.02e-02 (class 3) | -4.49e-02 | 3.011e-03 / 3.011e-03 (1.0000) |
| 4 | 128 | (29, 33, 28, 38) | +4.26e-02 (class 4) | -4.73e-02 | — |
| 5 | 16 | (4, 5, 4, 3) | +3.63e-02 (class 2) | -4.55e-02 | 4.178e-03 / 4.178e-03 (1.0000) |
| 5 | 32 | (10, 8, 5, 9) | +3.91e-02 (class 1) | -4.81e-02 | 3.020e-03 / 3.020e-03 (1.0000) |
| 5 | 64 | (18, 15, 13, 18) | +4.33e-02 (class 1) | -4.87e-02 | 3.737e-04 / 3.737e-04 (1.0000) |
| 5 | 128 | (34, 32, 28, 34) | +4.65e-02 (class 1) | -4.97e-02 | — |

0 of 18 m-systems (and 0 of the 6 references) sit at weights where the fixed prior's replicated point is a saddle. At the stable pairs the observed floor is 1.0000–1.0000 of the exact bias at the sampled weights; at the others —.

## C. The real closed loop (k_e = 25, dt = 0.002, 2500 time units), nearest-atom leak, seed 0, m = 64, counts (14, 13, 21, 16) (representability margin +0.019)

The rms tracking error over the five periods ending at the times below, against (1/k_e)‖f − f̂‖_ρ of the averaged loop from the same start at s = t/k_e.

| start, prior | t | e_rms | (1/k_e)‖f − f̂‖_ρ averaged | ratio |
|---|---|---|---|---|
| warm cloud, fixed prior | 31 | 3.354e-03 | 3.326e-03 | 1.009 |
| warm cloud, fixed prior | 500 | 2.740e-03 | 2.748e-03 | 0.997 |
| warm cloud, fixed prior | 1000 | 2.598e-03 | 2.616e-03 | 0.993 |
| warm cloud, fixed prior | 2500 | 2.552e-03 | 2.570e-03 | 0.993 |
| warm cloud, adaptive prior (η = ∞) | 31 | 3.354e-03 | 3.325e-03 | 1.009 |
| warm cloud, adaptive prior (η = ∞) | 500 | 2.470e-03 | 2.479e-03 | 0.996 |
| warm cloud, adaptive prior (η = ∞) | 1000 | 1.974e-03 | 2.004e-03 | 0.985 |
| warm cloud, adaptive prior (η = ∞) | 2500 | 1.320e-03 | 1.348e-03 | 0.979 |
| ϑ^e(w), adaptive prior at ϑ^e(w) | 31 | 4.695e-10 | 3.510e-18 | — |
| ϑ^e(w), adaptive prior at ϑ^e(w) | 500 | 3.463e-10 | 3.510e-18 | — |
| ϑ^e(w), adaptive prior at ϑ^e(w) | 1000 | 2.835e-10 | 3.510e-18 | — |
| ϑ^e(w), adaptive prior at ϑ^e(w) | 2500 | 1.977e-10 | 3.510e-18 | — |

Numbers are produced by `experiments/exp_adaptive_floor_2d.py`; nothing is edited by hand.
