# Proposition D‴.43: the colour of the measurement noise selects the kernels (`exp_colour_selection.py`)

w = u⊙u − v⊙v on a fixed dictionary Φ (N × d), fitted by gradient descent with label noise; along the exact fits the reduced flow dP/ds = −ω⊙P, dw/ds = −Π_D(ω⊙w), ω_i = φ_iᵀf(K)φ_i, the gradient flow of R = tr Ψ(K) (`colour_selection`). Success: within 0.0001·max|w*| of the truth. The start: the noiseless gradient flow's interpolant from u = v = 0.5 (the hyperbolic entropy's minimiser on M, P = α²); the end interpolants likewise from each path's end state.

## A. The potential's identities

| N × d | τ | det(I + τK) against its subset expansion (relative) | ∂R/∂D against (4/N)ω (relative) | R's second differences along 20 random segments (max, min) | white: |second difference| |
|---|---|---|---|---|---|
| 4 × 7 | 0.3 | 0.0e+00 | 5.6e-10 | -7.12e-03, -2.54e-01 | 1.8e-15 |
| 4 × 7 | 10 | 7.8e-16 | 1.6e-09 | -3.25e-03, -7.93e-02 | 3.6e-15 |
| 4 × 7 | 1000 | 2.6e-15 | 2.8e-09 | -8.03e-05, -2.18e-03 | 8.9e-16 |
| 5 × 9 | 0.3 | 4.4e-16 | 7.4e-10 | -2.29e-02, -3.87e-01 | 0.0e+00 |
| 5 × 9 | 10 | 8.9e-16 | 2.5e-09 | -2.16e-02, -1.29e-01 | 3.6e-15 |
| 5 × 9 | 1000 | 3.3e-16 | 3.9e-09 | -1.66e-04, -1.62e-03 | 0.0e+00 |

On a k-sparse exact fit (P = 0, D = |w|): τ_cR − k log(4τ_c/N) against log(det(Φ_SᵀΦ_S)Π_{i∈S}|w_i|) —

| N × d | k | the limit | τ_c = 10², 10⁴, 10⁶ |
|---|---|---|---|
| 10 × 40 | 3 | -0.065900 | +0.012561, -0.065104, -0.065892 |
| 10 × 40 | 6 | -1.168243 | -0.941086, -1.165909, -1.168220 |
| 20 × 80 | 5 | -1.141188 | -0.788862, -1.137498, -1.141151 |

## B. Label-noise descent against the reduced flow (10 × 40, a 4-sparse truth, step 0.5, 16 paths, the first 8 seeds)

σ = 0.1: the paths' mean against the flow (largest deviation over w and the nine times s = 0, S/8, …, S), the flow's error at S, the fraction of paths whose end interpolant recovers the truth (within 0.0001·max|w*|), the end interpolants' errors (median over the paths) and supports' sizes (median; n.c. not converged), the raw ends' error (median), basis pursuit's error:

| seed | the noise | steps | s at the end | the paths' mean against the flow | the flow's error at S | paths recovering the truth | the end interpolants' error | their support | the raw ends' error | basis pursuit's error |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | fresh | 3000000 | 300 | 0.024 | 2.4e+00 | 0.00 | 2e+00 | 13 | 2.40 | 2.98 |
| 0 | AR(1), τ_c = 100 | 4000000 | 400 | 0.129 | 5.1e-12 | 1.00 | 4e-13 | 4 | 0.15 | 2.98 |
| 1 | fresh | 3000000 | 300 | 0.016 | 1.5e+00 | 0.00 | 2e+00 | 10 | 1.51 | 1.51 |
| 1 | AR(1), τ_c = 100 | 4000000 | 400 | 0.353 | 6.4e-12 | 1.00 | 7e-11 | 4 | 0.23 | 1.51 |
| 2 | fresh | 3000000 | 300 | 0.035 | 1.3e+00 | 0.00 | 1e+00 | 10 | 1.32 | 1.32 |
| 2 | AR(1), τ_c = 100 | 4000000 | 400 | 0.145 | 1.7e-11 | 1.00 | 1e-13 | 4 | 0.17 | 1.32 |
| 3 | fresh | 3000000 | 300 | 0.014 | 1.3e+00 | 0.00 | 1e+00 | 10 | 1.32 | 1.33 |
| 3 | AR(1), τ_c = 100 | 4000000 | 400 | 0.248 | 1.7e-11 | 1.00 | 2e-13 | 4 | 0.17 | 1.33 |
| 4 | fresh | 3000000 | 300 | 0.015 | 8.0e-01 | 0.00 | 8e-01 | 10 | 0.80 | 0.80 |
| 4 | AR(1), τ_c = 100 | 4000000 | 400 | 0.109 | 7.6e-01 | 0.00 | n.c. | 11 | 0.87 | 0.80 |
| 5 | fresh | 3000000 | 300 | 0.011 | 1.3e+00 | 0.00 | 1e+00 | 10 | 1.27 | 1.27 |
| 5 | AR(1), τ_c = 100 | 4000000 | 400 | 0.110 | 1.3e+00 | 0.00 | 1e+00 | 11 | 1.41 | 1.27 |
| 6 | fresh | 3000000 | 300 | 0.016 | 1.2e+00 | 0.00 | 1e+00 | 12 | 1.22 | 1.22 |
| 6 | AR(1), τ_c = 100 | 4000000 | 400 | 0.134 | 1.2e+00 | 0.00 | n.c. | 11 | 1.30 | 1.22 |
| 7 | fresh | 3000000 | 300 | 0.015 | 9.4e-01 | 0.00 | 9e-01 | 11 | 0.95 | 0.94 |
| 7 | AR(1), τ_c = 100 | 4000000 | 400 | 0.292 | 9.1e-01 | 0.12 | n.c. | 11 | 0.78 | 0.94 |

Per seed (AR(1), τ_c = 100): the fraction of paths recovering the truth at σ = 0.1, 0.3, 0.5 (the same rescaled time), against the reduced flow's own success and basis pursuit's:

| seed | the reduced flow recovers | basis pursuit recovers | σ = 0.1 | σ = 0.3 | σ = 0.5 |
|---|---|---|---|---|---|
| 0 | yes | no | 1.00 | 1.00 | 1.00 |
| 1 | yes | no | 1.00 | 0.00 | 0.00 |
| 2 | yes | no | 1.00 | 1.00 | 0.00 |
| 3 | yes | no | 1.00 | 0.81 | 0.00 |
| 4 | no | no | 0.00 | 0.44 | 0.00 |
| 5 | no | no | 0.00 | 0.00 | 0.00 |
| 6 | no | no | 0.00 | 0.00 | 0.00 |
| 7 | no | no | 0.12 | 0.06 | 0.00 |

## C. Recovery against sparsity (the reduced flows, continuum filters)

| N × d | k | instances | bp | white | slow 10 | slow 100 | slow 1000 | white, then slow 100 | the slow flow's (τ_c = 100) mean support size | basis pursuit's |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 × 40 | 1 | 40 | 0.97 | 0.97 | 1.00 | 1.00 | 1.00 | 0.97 | 1.0 | 1.2 |
| 10 × 40 | 2 | 40 | 0.90 | 0.90 | 0.95 | 0.97 | 1.00 | 0.90 | 2.2 | 2.8 |
| 10 × 40 | 3 | 40 | 0.38 | 0.38 | 0.62 | 0.75 | 0.78 | 0.38 | 4.8 | 7.3 |
| 10 × 40 | 4 | 40 | 0.17 | 0.17 | 0.42 | 0.60 | 0.62 | 0.17 | 6.4 | 8.9 |
| 10 × 40 | 5 | 40 | 0.00 | 0.00 | 0.07 | 0.15 | 0.17 | 0.00 | 9.2 | 10.0 |
| 10 × 40 | 6 | 40 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 10.1 | 10.0 |
| 20 × 80 | 2 | 20 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 2.0 | 2.0 |
| 20 × 80 | 4 | 20 | 0.90 | 0.90 | 0.95 | 1.00 | 1.00 | 0.95 | 4.0 | 5.6 |
| 20 × 80 | 6 | 20 | 0.35 | 0.35 | 0.65 | 0.75 | 0.75 | 0.40 | 9.4 | 15.1 |
| 20 × 80 | 8 | 20 | 0.10 | 0.10 | 0.20 | 0.35 | 0.30 | 0.10 | 15.6 | 18.6 |
| 20 × 80 | 10 | 20 | 0.05 | 0.05 | 0.00 | 0.10 | 0.10 | 0.05 | 18.9 | 19.4 |

The white flow run to convergence (max P < 10⁻⁸·max|w|, w still to 10⁻¹⁰, s ≤ 1e+05): converged in 336 of 340 instances, its largest distance to basis pursuit's solution 9.4e-09. Where a slow flow fails, the fraction of its ends with R below the truth's (the truth not the global minimiser) — τ_c = 10, 100:

| N × d | k | failures (10, 100) | R(end) < R(truth) |
|---|---|---|---|
| 10 × 40 | 1 | 0, 0 | —, — |
| 10 × 40 | 2 | 2, 1 | 0.00, 0.00 |
| 10 × 40 | 3 | 15, 10 | 0.00, 0.00 |
| 10 × 40 | 4 | 23, 16 | 0.39, 0.00 |
| 10 × 40 | 5 | 37, 34 | 0.68, 0.09 |
| 10 × 40 | 6 | 40, 40 | 0.85, 0.35 |
| 20 × 80 | 2 | 0, 0 | —, — |
| 20 × 80 | 4 | 1, 0 | 0.00, — |
| 20 × 80 | 6 | 7, 5 | 0.14, 0.00 |
| 20 × 80 | 8 | 16, 13 | 0.38, 0.00 |
| 20 × 80 | 10 | 20, 18 | 0.75, 0.11 |

## D. Gaussian kernel dictionaries on [0, 1] (40 centres, 20 equispaced samples), random k-atom truths (40 per row)

| width | the largest off-diagonal |φ_iᵀφ_j| | k | bp | white | slow 10 | slow 100 | mean support size: basis pursuit, slow 100 |
|---|---|---|---|---|---|---|---|
| 0.03 | 0.90 | 3 | 0.85 | 0.85 | 0.90 | 0.95 | 3.7, 3.2 |
| 0.03 | 0.90 | 4 | 0.80 | 0.80 | 0.80 | 0.95 | 5.0, 4.2 |
| 0.03 | 0.90 | 5 | 0.70 | 0.70 | 0.80 | 0.90 | 6.6, 5.3 |
| 0.03 | 0.90 | 6 | 0.45 | 0.45 | 0.55 | 0.68 | 8.8, 6.9 |
| 0.05 | 0.96 | 3 | 0.95 | 0.95 | 0.95 | 1.00 | 3.5, 3.0 |
| 0.05 | 0.96 | 4 | 0.85 | 0.85 | 0.93 | 1.00 | 5.3, 4.0 |
| 0.05 | 0.96 | 5 | 0.82 | 0.82 | 0.85 | 0.93 | 6.7, 5.6 |
| 0.05 | 0.96 | 6 | 0.65 | 0.65 | 0.60 | 0.78 | 9.2, 7.4 |

Label-noise descent on the first instances where basis pursuit fails (width 0.05, k = 4; 16 paths, σ = 0.1, the step min(0.5, 1/(2λ_max)) at the start):

| instance | step | the noise | steps | the end interpolants' errors (per path) | their supports' sizes | the ends' atoms alive (D > 10⁻⁸ max D) | the raw ends' distance to the truth, to basis pursuit's solution (max over paths) | the reduced flow at the same s: its error, support, distance to basis pursuit's solution | the reduced flow run to convergence: its error, support, the s reached |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.5 | fresh | 12000000 | n.c., n.c., n.c., n.c., n.c., 1e+00, 1e+00, n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | n.c., n.c., n.c., n.c., n.c., 14, 14, n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | 15, 15, 15, 14, 15, 15, 15, 15, 14, 14, 15, 15, 15, 15, 14, 14 | 1.20, 0.11 | 1.2e+00, 14, 2.2e-08 | 1.2e+00, 14, 600 |
| 0 | 0.5 | AR(1), τ_c = 100 | 16000000 | 2e-13, 3e-13, 6e-13, 3e-12, 5e-12, 2e-12, 1e-12, 5e-13, 6e-13, 3e-12, 9e-12, 4e-12, 4e-12, 6e-13, 2e-12, 4e-12 | 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 | 5, 6, 7, 7, 7, 10, 5, 6, 7, 9, 7, 5, 6, 6, 5, 5 | 0.70, 1.33 | 6.9e-12, 4, 1.2e+00 | 1.5e-11, 4, 800 |
| 7 | 0.5 | fresh | 12000000 | n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | 15, 15, 14, 14, 13, 14, 14, 14, 13, 14, 13, 14, 14, 14, 15, 14 | 0.73, 0.11 | 7.3e-01, 18, 1.3e-09 | 7.3e-01, 18, 600 |
| 7 | 0.5 | AR(1), τ_c = 100 | 16000000 | 6e-08, 5e-08, 5e-08, 4e-08, 3e-08, 4e-08, 4e-08, 4e-08, 4e-08, 4e-08, 4e-08, 4e-08, 7e-08, 4e-08, 4e-08, 5e-08 | 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 | 9, 8, 9, 8, 10, 7, 9, 8, 8, 8, 9, 9, 8, 8, 10, 8 | 0.48, 0.95 | 1.8e-08, 4, 7.3e-01 | 3.1e-13, 4, 1600 |
| 13 | 0.5 | fresh | 12000000 | n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c., n.c. | 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 | 1.09, 0.08 | 1.0e+00, 15, 3.8e-07 | 1.0e+00, 15, 2400 |
| 13 | 0.5 | AR(1), τ_c = 100 | 16000000 | 6e-01, 7e-01, 7e-01, 7e-01, 7e-01, 7e-01, 7e-01, 7e-01, 7e-01, 7e-01, 7e-01, 7e-01, 7e-01, 6e-01, 7e-01, 6e-01 | 12, 13, 12, 12, 12, 12, 12, 12, 12, 13, 12, 13, 13, 13, 12, 13 | 11, 12, 10, 12, 12, 12, 12, 13, 11, 11, 12, 14, 11, 12, 12, 10 | 1.06, 0.50 | 8.8e-03, 8, 1.0e+00 | 3.0e-13, 4, 1600 |

n.c.: the end's flow interpolant does not converge — the exact fit needs atoms whose masses have underflowed (fewer than N alive, a residual outside their span), its maximiser is at infinity; the raw end (residual O(σ)) is reported instead.

## E. Other colours at the step 0.5 (the reduced flow with each scheme's filter; 10 × 40, 20 instances)

| k | fresh | AR(1), τ_c = 10 | smooth (Gaussian), τ_c = 10 | kept 40 steps | smooth, τ_c = 100, 1× the slow horizon | smooth, τ_c = 100, 10× the slow horizon | basis pursuit |
|---|---|---|---|---|---|---|---|
| 3 | 0.30 | 0.65 | 0.70 | 0.70 | 0.00 | 0.85 | 0.45 |
| 4 | 0.20 | 0.40 | 0.55 | 0.40 | 0.00 | 0.65 | 0.25 |

## F. The start: the interpolant of u = v = α (10 × 40), the flows run to convergence (s ≤ 1e+05)

| k | α | the interpolant itself | the white flow | the slow flow (τ_c = 100) | converged (white, slow) | median s reached (white, slow) | basis pursuit |
|---|---|---|---|---|---|---|---|
| 3 | 0.1 | 0.00 | 0.45 | 0.70 | 1.00, 1.00 | 4e+03, 2.01e+03 | 0.45 |
| 3 | 0.5 | 0.00 | 0.45 | 0.75 | 1.00, 1.00 | 4e+03, 2.01e+03 | 0.45 |
| 3 | 2 | 0.00 | 0.45 | 0.35 | 1.00, 1.00 | 4e+03, 4.01e+03 | 0.45 |
| 4 | 0.1 | 0.00 | 0.25 | 0.40 | 1.00, 1.00 | 4e+03, 4.01e+03 | 0.25 |
| 4 | 0.5 | 0.00 | 0.25 | 0.50 | 1.00, 1.00 | 4e+03, 2.01e+03 | 0.25 |
| 4 | 2 | 0.00 | 0.25 | 0.30 | 1.00, 1.00 | 4e+03, 4.01e+03 | 0.25 |

## G. References on the same instances (10 × 40, 40 per k): reweighted ℓ¹ (Candès–Wakin–Boyd; 12 iterations, ε relative to max|w*|, from basis pursuit's solution or from the flows' start), the log-det's majorisation–minimisation (τ_c = 100), sparse Bayesian learning (Tipping's EM, λ = 10⁻¹⁰) — the fraction recovering the truth, and the slow flow's (τ_c = 100, Part C) beside them

| method | k = 2 | k = 3 | k = 4 | k = 5 | k = 6 | mean support (k = 4) |
|---|---|---|---|---|---|---|
| bp | 0.90 | 0.38 | 0.17 | 0.00 | 0.00 | 8.9 |
| reweighted ℓ¹ (ε = 0.01), from basis pursuit | 0.90 | 0.38 | 0.17 | 0.00 | 0.00 | 8.9 |
| reweighted ℓ¹ (ε = 0.01), from the start | 0.78 | 0.47 | 0.23 | 0.10 | 0.00 | 8.6 |
| reweighted ℓ¹ (ε = 0.1), from basis pursuit | 0.90 | 0.38 | 0.17 | 0.00 | 0.00 | 8.9 |
| reweighted ℓ¹ (ε = 0.1), from the start | 0.78 | 0.47 | 0.23 | 0.07 | 0.00 | 8.6 |
| reweighted ℓ¹ (ε = 1), from basis pursuit | 0.90 | 0.53 | 0.23 | 0.05 | 0.00 | 8.6 |
| reweighted ℓ¹ (ε = 1), from the start | 0.78 | 0.47 | 0.25 | 0.05 | 0.00 | 8.4 |
| log-det MM from basis pursuit | 0.93 | 0.53 | 0.33 | 0.05 | 0.00 | 8.0 |
| log-det MM from the start | 0.95 | 0.70 | 0.45 | 0.10 | 0.00 | 7.2 |
| sparse Bayesian learning (EM) | 0.97 | 0.88 | 0.55 | 0.12 | 0.05 | 6.6 |
| the slow flow (τ_c = 100) | 0.97 | 0.75 | 0.60 | 0.15 | 0.00 | 6.4 |

## H. Which exact fits the selection holds: the dual certificate max_{i∉S}|c_i|/ω_i (below 1: the fit is stable)

| k | instances | at the truth: the fraction certified, white; τ_c = 10, 100, 1000 | the median ratio at the truth, white; τ_c = 10, 100, 1000 | where basis pursuit fails: its solution's support size, and the fraction of its solutions certified for τ_c = 10, 100, 1000 | basis pursuit's success |
|---|---|---|---|---|---|
| 2 | 40 | 0.38; 0.97; 1.00; 1.00 | 1.2; 0.439; 0.0814; 0.0087 | 10.0; 0.00; 0.00; 0.00 | 0.90 |
| 3 | 40 | 0.03; 0.93; 1.00; 1.00 | 1.72; 0.585; 0.111; 0.0126 | 10.0; 0.24; 0.12; 0.24 | 0.38 |
| 4 | 40 | 0.03; 0.75; 1.00; 1.00 | 1.98; 0.788; 0.185; 0.0217 | 10.0; 0.21; 0.24; 0.30 | 0.17 |
| 5 | 40 | 0.00; 0.50; 1.00; 1.00 | 2.37; 0.997; 0.268; 0.0393 | 10.0; 0.25; 0.23; 0.23 | 0.00 |
| 6 | 40 | 0.00; 0.30; 0.93; 1.00 | 2.79; 1.19; 0.407; 0.0661 | 10.0; 0.12; 0.07; 0.10 | 0.00 |

The certificate is sufficient, not necessary (for white draws it is Fuchs's least-squares certificate): basis pursuit's own success is in the last column.

On the kernel dictionaries (40 centres, 20 samples; 20 instances per row): the fraction of truths certified, white and τ_c = 1e+02, 1e+03, 1e+04, 1e+05; the fraction to which the slow flow returns from P = 10⁻³ (τ_c = 100, 1000); basis pursuit's success:

| width | k | certified: white; 1e+02; 1e+03; 1e+04; 1e+05 | returns (100, 1000) | basis pursuit |
|---|---|---|---|---|
| 0.03 | 4 | 0.60; 1.00; 1.00; 1.00; 1.00 | 1.00, 1.00 | 0.75 |
| 0.03 | 6 | 0.15; 0.85; 0.90; 0.95; 1.00 | 0.90, 1.00 | 0.45 |
| 0.05 | 4 | 0.10; 0.90; 1.00; 1.00; 1.00 | 1.00, 1.00 | 0.80 |
| 0.05 | 6 | 0.00; 0.60; 0.85; 0.95; 1.00 | 1.00, 1.00 | 0.65 |
| 0.08 | 4 | 0.00; 0.70; 0.90; 1.00; 1.00 | 1.00, 1.00 | 0.85 |
| 0.08 | 6 | 0.00; 0.10; 0.60; 0.75; 0.85 | 0.90, 1.00 | 0.65 |

## I. A fit floor: y = Φw* + ε·e (|e| = 1; 10 × 40, the instances of Part C) — the fraction whose k largest atoms are the truth's support, and the support sizes (|w| > 10⁻³·max|w*|)

| k | ε | the slow flow (τ_c = 100) | basis pursuit | support sizes (slow, basis pursuit) |
|---|---|---|---|---|
| 3 | 0.001 | 0.75 | 0.45 | 4.9, 7.5 |
| 3 | 0.01 | 0.75 | 0.45 | 8.4, 9.3 |
| 3 | 0.05 | 0.70 | 0.40 | 9.7, 9.8 |
| 4 | 0.001 | 0.50 | 0.25 | 7.2, 8.8 |
| 4 | 0.01 | 0.50 | 0.25 | 8.9, 9.7 |
| 4 | 0.05 | 0.50 | 0.25 | 9.7, 9.9 |

## J. The selection is the explicit scheme's: label-noise descent at the same physical noise (AR(1) of correlation time 50 time units, per-step σ = 0.1) over the same physical time 200000 (16 paths; Part C's first 10 × 40, k = 4 instance where basis pursuit fails)

| scheme | step | steps | mean Δlog P off the truth's support (± s.e.) | the end interpolants recovering the truth |
|---|---|---|---|---|
| Euler | 0.5 | 400000 | -2.3907 ± 1.9e-03 | 0.00 |
| Heun | 0.5 | 400000 | -0.0967 ± 1.9e-05 | 0.00 |
| Euler | 0.125 | 1600000 | -0.2781 ± 8.8e-05 | 0.00 |

## K. Statistics and the competing fits (referee round 142)

K1. The recovery rates of Parts C and G with exact (Clopper–Pearson) 95% intervals — successes/instances, rate [interval]:

| N × d | k | bp | slow 10 | slow 100 | slow 1000 | sparse Bayesian learning (EM) |
|---|---|---|---|---|---|---|
| 10 × 40 | 1 | 39/40, 0.97 [0.87, 1.00] | 40/40, 1.00 [0.91, 1.00] | 40/40, 1.00 [0.91, 1.00] | 40/40, 1.00 [0.91, 1.00] | — |
| 10 × 40 | 2 | 36/40, 0.90 [0.76, 0.97] | 38/40, 0.95 [0.83, 0.99] | 39/40, 0.97 [0.87, 1.00] | 40/40, 1.00 [0.91, 1.00] | 39/40, 0.97 [0.87, 1.00] |
| 10 × 40 | 3 | 15/40, 0.38 [0.23, 0.54] | 25/40, 0.62 [0.46, 0.77] | 30/40, 0.75 [0.59, 0.87] | 31/40, 0.78 [0.62, 0.89] | 35/40, 0.88 [0.73, 0.96] |
| 10 × 40 | 4 | 7/40, 0.17 [0.07, 0.33] | 17/40, 0.42 [0.27, 0.59] | 24/40, 0.60 [0.43, 0.75] | 25/40, 0.62 [0.46, 0.77] | 22/40, 0.55 [0.38, 0.71] |
| 10 × 40 | 5 | 0/40, 0.00 [0.00, 0.09] | 3/40, 0.07 [0.02, 0.20] | 6/40, 0.15 [0.06, 0.30] | 7/40, 0.17 [0.07, 0.33] | 5/40, 0.12 [0.04, 0.27] |
| 10 × 40 | 6 | 0/40, 0.00 [0.00, 0.09] | 0/40, 0.00 [0.00, 0.09] | 0/40, 0.00 [0.00, 0.09] | 0/40, 0.00 [0.00, 0.09] | 2/40, 0.05 [0.01, 0.17] |
| 20 × 80 | 2 | 20/20, 1.00 [0.83, 1.00] | 20/20, 1.00 [0.83, 1.00] | 20/20, 1.00 [0.83, 1.00] | 20/20, 1.00 [0.83, 1.00] | — |
| 20 × 80 | 4 | 18/20, 0.90 [0.68, 0.99] | 19/20, 0.95 [0.75, 1.00] | 20/20, 1.00 [0.83, 1.00] | 20/20, 1.00 [0.83, 1.00] | — |
| 20 × 80 | 6 | 7/20, 0.35 [0.15, 0.59] | 13/20, 0.65 [0.41, 0.85] | 15/20, 0.75 [0.51, 0.91] | 15/20, 0.75 [0.51, 0.91] | — |
| 20 × 80 | 8 | 2/20, 0.10 [0.01, 0.32] | 4/20, 0.20 [0.06, 0.44] | 7/20, 0.35 [0.15, 0.59] | 6/20, 0.30 [0.12, 0.54] | — |
| 20 × 80 | 10 | 1/20, 0.05 [0.00, 0.25] | 0/20, 0.00 [0.00, 0.17] | 2/20, 0.10 [0.01, 0.32] | 2/20, 0.10 [0.01, 0.32] | — |

Paired on the same instances: the slow flow (τ_c = 100) against another method — the instances only the slow flow recovers, those only the other recovers, and the exact two-sided McNemar p-value (the binomial test on the discordant pairs):

| N × d | k | instances | against | slow only | the other only | p |
|---|---|---|---|---|---|---|
| 10 × 40 | 1 | 40 | bp | 1 | 0 | 1 |
| 10 × 40 | 1 | 40 | slow 1000 | 0 | 0 | 1 |
| 10 × 40 | 2 | 40 | bp | 3 | 0 | 0.25 |
| 10 × 40 | 2 | 40 | sparse Bayesian learning (EM) | 0 | 0 | 1 |
| 10 × 40 | 2 | 40 | slow 1000 | 0 | 1 | 1 |
| 10 × 40 | 3 | 40 | bp | 15 | 0 | 6.1e-05 |
| 10 × 40 | 3 | 40 | sparse Bayesian learning (EM) | 1 | 6 | 0.12 |
| 10 × 40 | 3 | 40 | slow 1000 | 0 | 1 | 1 |
| 10 × 40 | 4 | 40 | bp | 17 | 0 | 1.5e-05 |
| 10 × 40 | 4 | 40 | sparse Bayesian learning (EM) | 4 | 2 | 0.69 |
| 10 × 40 | 4 | 40 | slow 1000 | 0 | 1 | 1 |
| 10 × 40 | 5 | 40 | bp | 6 | 0 | 0.031 |
| 10 × 40 | 5 | 40 | sparse Bayesian learning (EM) | 3 | 2 | 1 |
| 10 × 40 | 5 | 40 | slow 1000 | 0 | 1 | 1 |
| 10 × 40 | 6 | 40 | bp | 0 | 0 | 1 |
| 10 × 40 | 6 | 40 | sparse Bayesian learning (EM) | 0 | 2 | 0.5 |
| 10 × 40 | 6 | 40 | slow 1000 | 0 | 0 | 1 |
| 20 × 80 | 2 | 20 | bp | 0 | 0 | 1 |
| 20 × 80 | 2 | 20 | slow 1000 | 0 | 0 | 1 |
| 20 × 80 | 4 | 20 | bp | 2 | 0 | 0.5 |
| 20 × 80 | 4 | 20 | slow 1000 | 0 | 0 | 1 |
| 20 × 80 | 6 | 20 | bp | 8 | 0 | 0.0078 |
| 20 × 80 | 6 | 20 | slow 1000 | 1 | 1 | 1 |
| 20 × 80 | 8 | 20 | bp | 5 | 0 | 0.062 |
| 20 × 80 | 8 | 20 | slow 1000 | 1 | 0 | 1 |
| 20 × 80 | 10 | 20 | bp | 1 | 0 | 1 |
| 20 × 80 | 10 | 20 | slow 1000 | 0 | 0 | 1 |

K2. Where the slow flow (τ_c = 100) fails at its fixed horizon (Part C, 10 × 40), the same flow run to convergence (s ≤ 1e+05): the failures, how many converged and still fail (agreeing with Part C), the end's support (median, range), the fraction of ends that are basic exact fits on N = 10 atoms, the fraction certified at the end (max_{i∉S}|c_i|/ω_i < 1) for τ_c = 100, 10³, 10⁴, 10⁵ and in the τ_c → ∞ limit, that limit's median:

| k | failures | converged, failing | the end's support | basic on N atoms | certified at the end: τ_c = 100; 10³; 10⁴; 10⁵; ∞ | the limit's median ratio |
|---|---|---|---|---|---|---|
| 3 | 10 | 10 | 10 (10–10) | 1.00 | 1.00; 0.70; 0.60; 0.30; 0.30 | 1.13 |
| 4 | 16 | 16 | 10 (10–10) | 1.00 | 1.00; 0.81; 0.62; 0.62; 0.56 | 0.852 |
| 5 | 34 | 34 | 10 (10–10) | 1.00 | 1.00; 0.68; 0.44; 0.38; 0.38 | 1.42 |
| 6 | 40 | 40 | 10 (10–10) | 1.00 | 1.00; 0.65; 0.57; 0.47; 0.47 | 1.03 |

R at the converged end against R at the truth, the same two points evaluated with the potential of each τ_c (the slow flow ran at τ_c = 100): the fraction of the failures whose end has the lower R, and the median of τ_c(R(end) − R(truth)):

| k | τ_c = 1 | τ_c = 10 | τ_c = 100 | τ_c = 1000 | τ_c = 10000 | τ_c = 100000 |
|---|---|---|---|---|---|---|
| 3 | 0.90; -0.0873 | 0.00; +0.459 | 0.00; +5.18 | 0.00; +15 | 0.00; +28.9 | 0.00; +44.4 |
| 4 | 0.94; -0.22 | 0.38; +0.179 | 0.00; +4.25 | 0.00; +13.4 | 0.00; +25.7 | 0.00; +39.3 |
| 5 | 0.91; -0.342 | 0.71; -0.383 | 0.09; +1.79 | 0.00; +8.11 | 0.00; +17.8 | 0.00; +28.8 |
| 6 | 0.93; -0.33 | 0.82; -0.824 | 0.38; +0.392 | 0.07; +5.46 | 0.00; +13.1 | 0.00; +22 |

Numbers are produced by `experiments/exp_colour_selection.py`; nothing is edited by hand.
