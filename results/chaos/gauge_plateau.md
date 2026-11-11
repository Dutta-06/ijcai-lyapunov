# Proposition D‴.39: dissipativity of the dithered gauge (`exp_gauge_plateau.py`)

Proposition D‴.36's instance (three bumps, one particle per class, m = 3, the cohesion leak); energies over τ = 0.0001; the box L = (5.0, 3.0, 4.0) in (α, γ, κ), σ_B = 1: |a| ≤ 3.393, the representable masses W_k ≥ 0.2947, 0.2063, 0.2358 (the saturation's 0.2500, 0.1750, 0.2000).

## A. The plateau: the lowest configuration with one amplitude fixed far out

| masses | class | α_k | ΔJ/τ above the exact fit | fit error | a | c | s |
|---|---|---|---|---|---|---|---|
| uniform | 1 | ±20 | 7.718 | 0.0350 | 4.00, -3.73, 3.15 | 1.09, 0.34, -1.41 | 0.72, 0.74, 0.58 |
| uniform | 2 | ±20 | 7.762 | 0.0351 | 3.11, -4.00, 4.00 | -1.42, 0.32, 1.04 | 0.58, 0.70, 0.75 |
| uniform | 3 | ±20 | 7.718 | 0.0350 | 3.15, -3.73, 4.00 | -1.41, 0.34, 1.09 | 0.58, 0.74, 0.72 |
| sharp | 1 | ±20 | 7.637 | 0.0348 | 4.00, -3.49, 4.00 | -1.32, -0.02, 1.28 | 0.66, 1.17, 0.67 |
| sharp | 2 | ±20 | 7.707 | 0.0350 | 2.78, -4.00, 3.89 | -1.41, 0.31, 1.10 | 0.58, 0.72, 0.73 |
| sharp | 3 | ±20 | 7.637 | 0.0348 | 4.00, -3.49, 4.00 | -1.32, -0.02, 1.28 | 0.66, 1.17, 0.67 |
| uniform | 1 | ±6 | 7.692 | 0.0349 | 3.62, -3.38, 3.14 | 1.13, 0.29, -1.41 | 0.72, 0.72, 0.58 |
| uniform | 2 | ±6 | 7.710 | 0.0350 | 3.14, -3.62, 3.83 | -1.41, 0.32, 1.10 | 0.58, 0.73, 0.73 |
| uniform | 3 | ±6 | 7.692 | 0.0349 | 3.14, -3.38, 3.62 | -1.41, 0.29, 1.13 | 0.58, 0.72, 0.72 |
| sharp | 1 | ±6 | 7.722 | 0.0350 | 3.62, 3.52, -3.92 | 1.08, -1.41, 0.35 | 0.73, 0.58, 0.74 |
| sharp | 2 | ±6 | 7.670 | 0.0349 | 2.77, -3.62, 3.51 | -1.42, 0.26, 1.15 | 0.58, 0.70, 0.72 |
| sharp | 3 | ±6 | 7.681 | 0.0349 | 2.78, -3.67, 3.62 | -1.41, 0.28, 1.14 | 0.58, 0.72, 0.72 |
| uniform | 1 | ±L_α (box) | 7.662 | 0.0349 | 3.39, -3.17, 3.14 | 1.16, 0.26, -1.41 | 0.71, 0.71, 0.58 |
| uniform | 2 | ±L_α (box) | 7.735 | 0.0350 | 3.09, -3.39, 3.39 | -1.43, 0.26, 1.11 | 0.57, 0.68, 0.74 |
| uniform | 3 | ±L_α (box) | 7.662 | 0.0349 | 3.14, -3.17, 3.39 | -1.41, 0.26, 1.16 | 0.58, 0.71, 0.71 |
| sharp | 1 | ±L_α (box) | 8.048 | 0.0357 | 3.39, -3.18, 3.39 | -1.32, 0.01, 1.27 | 0.64, 0.92, 0.66 |
| sharp | 2 | ±L_α (box) | 7.626 | 0.0348 | 2.75, -3.39, 3.27 | -1.42, 0.23, 1.18 | 0.58, 0.69, 0.71 |
| sharp | 3 | ±L_α (box) | 7.650 | 0.0348 | 2.79, -3.39, 3.39 | -1.41, 0.24, 1.18 | 0.58, 0.72, 0.70 |

## B. The minimum-energy path from the exact fit to the plateau (string method)

Masses uniform, toward class 3's saturation: the largest ΔJ/τ along the path 7.715 (image 39 of 39); at every fourth image the profile 0.00, 1.86, 5.33, 7.17, 7.57, 7.66, 7.70, 7.71, 7.71, 7.71, the fit error 0.0000, 0.0172, 0.0291, 0.0337, 0.0346, 0.0349, 0.0349, 0.0350, 0.0350, 0.0350, the second width 0.40, 0.46, 0.54, 0.63, 0.68, 0.71, 0.74, 0.75, 0.75, 0.74.

Masses sharp, toward class 3's saturation: the largest ΔJ/τ along the path 8.211 (image 6 of 39); at every fourth image the profile 0.00, 8.04, 8.19, 8.17, 8.12, 7.90, 7.79, 7.74, 7.70, 7.68, the fit error 0.0000, 0.0357, 0.0360, 0.0360, 0.0359, 0.0354, 0.0352, 0.0350, 0.0349, 0.0349, the second width 0.40, 0.77, 0.91, 0.93, 0.96, 1.07, 1.12, 1.14, 1.15, 1.16.

## C. The escapes at fixed masses 1/3 (γ = 0; 24 seeds, at most 400000 slow units; an escape: max|α_k| > 6)

| τ | escapes | the rate (per slow unit, ± s.e.) | the median time | the rate × e^{ΔJ/τ} | the fit error at the escape |
|---|---|---|---|---|---|
| 0.00015 | 10/24 | 1.33e-6 ± 4.22e-7 | 4e+05 | 2.25e-4 | 0.0456 |
| 0.0002 | 19/24 | 4.06e-6 ± 9.32e-7 | 1.6e+05 | 1.90e-4 | 0.0465 |
| 0.0003 | 24/24 | 1.09e-5 ± 2.23e-6 | 6.65e+04 | 1.42e-4 | 0.0538 |

The slope of log(rate) against 1/τ: 6.31e-4 (the lowest J on the escape set's boundary, ΔJ = 7.69e-4); weighted by the escape counts (log(rate)'s variance ≈ 1/n): 6.25e-4 ± 1.09e-4, 0.81 ± 0.14 of ΔJ (referee round 31).

## D. The dead zone at fixed masses 1/3 (γ = 0, τ = 0.0002, 16 seeds, 400000 slow units)

A record is on the plateau when its second width exceeds 0.6 (A and B: the plateaus' configurations and the path's).

| box | seeds that reach the plateau | fraction of the records on it (s_2 > 0.6), by quarters | fraction with max|α| > 6 | the largest |α|, by quarters | the fit error, by quarters |
|---|---|---|---|---|---|
| none | 16/16 | 0.179, 0.401, 0.535, 0.568 | 0.090, 0.358, 0.461, 0.559 | 15.0, 19.8, 21.0, 27.6 | 0.0395, 0.0425, 0.0439, 0.0448 |
| L = (5.0, 3.0, 4.0) | 16/16 | 0.123, 0.086, 0.185, 0.157 | 0.000, 0.000, 0.000, 0.000 | 5.0, 5.0, 5.0, 5.0 | 0.0391, 0.0392, 0.0403, 0.0403 |

## E. The full loop with the box (γ = 0.01, τ = 0.0001, 12 seeds, 1e+06 slow units)

The box's vertices: 0.5579, 0.2063, 0.2358; 0.2947, 0.4695, 0.2358; 0.2947, 0.2063, 0.4990; W♯ = 0.3765, 0.2988, 0.3247.

| κ_F | s | mean of W_1, W_2, W_3 (± s.e.) | std |
|---|---|---|---|
| 0 | 100000 | 0.3578 ± 0.0111, 0.3272 ± 0.0105, 0.3150 ± 0.0125 | 0.0385, 0.0364, 0.0434 |
| 0 | 300000 | 0.3672 ± 0.0158, 0.3172 ± 0.0147, 0.3156 ± 0.0113 | 0.0547, 0.0509, 0.0392 |
| 0 | 1e+06 | 0.3658 ± 0.0141, 0.3171 ± 0.0216, 0.3171 ± 0.0164 | 0.0490, 0.0747, 0.0569 |
| 2 | 100000 | 0.3692 ± 0.0077, 0.3070 ± 0.0074, 0.3237 ± 0.0077 | 0.0267, 0.0255, 0.0265 |
| 2 | 300000 | 0.3695 ± 0.0103, 0.3017 ± 0.0099, 0.3288 ± 0.0068 | 0.0358, 0.0342, 0.0236 |
| 2 | 1e+06 | 0.3803 ± 0.0093, 0.2961 ± 0.0122, 0.3236 ± 0.0102 | 0.0322, 0.0424, 0.0354 |

| κ_F | the masses at the end, by seed (the nearest vertex for κ_F = 0, its distance) | the largest |α| | the fit error (second half) |
|---|---|---|---|
| 0 | 0.342, 0.389, 0.269 (v2, 0.099); 0.371, 0.221, 0.408 (v3, 0.120); 0.296, 0.347, 0.356 (v2, 0.172); 0.470, 0.224, 0.306 (v1, 0.114); 0.354, 0.292, 0.355 (v3, 0.178); 0.316, 0.348, 0.336 (v2, 0.159); 0.316, 0.375, 0.309 (v2, 0.122); 0.388, 0.227, 0.384 (v3, 0.149); 0.366, 0.399, 0.235 (v2, 0.100); 0.364, 0.395, 0.242 (v2, 0.102); 0.431, 0.222, 0.347 (v1, 0.170); 0.377, 0.367, 0.257 (v2, 0.133) | 5.03 | 0.0257 |
| 2 | 0.373, 0.308, 0.319 (0.010); 0.395, 0.260, 0.345 (0.039); 0.300, 0.375, 0.325 (0.076); 0.404, 0.263, 0.333 (0.036); 0.391, 0.311, 0.299 (0.026); 0.350, 0.290, 0.360 (0.036); 0.361, 0.275, 0.364 (0.039); 0.410, 0.282, 0.308 (0.034); 0.417, 0.290, 0.294 (0.040); 0.376, 0.378, 0.246 (0.080); 0.386, 0.242, 0.373 (0.057); 0.402, 0.279, 0.318 (0.026) | 5.01 | 0.026 |

κ_F = 0: the records coincide with exp_gauge_drift.py's Part H (the same seeds without the box) up to the first record at which a coordinate has left B: 179000, 101000, 23000, 411000, 348000, 17000, 92000, 76000, 130000, 65000, 236000, 270000.

κ_F = 2: the records coincide with exp_gauge_fixman.py's Part E (the same seeds without the box) up to the first record at which a coordinate has left B: never, 372000, 24000, 691000, 484000, 70000, 316000, 257000, 207000, 354000, 658000, 279000.

## F. The two channels to the saturation in the stored runs

| runs | class | seeds with records at |α_k| > 10 | such records | of them with W_k − |A_k|/a_max < 0.01 | > 0.05 | their mean fit error |
|---|---|---|---|---|---|---|
| no leak (exp_gauge_drift.py H) | 1 | 8/12 | 2350 | 2349 | 0 | 0.0245 |
| no leak (exp_gauge_drift.py H) | 2 | 3/12 | 138 | 138 | 0 | 0.0244 |
| no leak (exp_gauge_drift.py H) | 3 | 7/12 | 1373 | 1373 | 0 | 0.0244 |
| the leak, no correction (exp_gauge_fixman.py F) | 1 | 2/16 | 516 | 0 | 516 | 0.0427 |
| the leak, no correction (exp_gauge_fixman.py F) | 2 | 0/16 | 0 | 0 | 0 | — |
| the leak, no correction (exp_gauge_fixman.py F) | 3 | 1/16 | 13 | 0 | 13 | 0.0408 |

## G. The escapes from a bounded domain with a smooth boundary (referee round 143; fixed masses 1/3, γ = 0, 48 seeds per τ)

The domain 𝒪 = {Σ_k,j ((ϑ_kj − ϑ*_kj)/L_j)^4 < 1}, ϑ* the exact fit, L = (2.5, 0.8, 1.2) in (α, γ, κ): bounded in the amplitudes', centres' and widths' coordinates, its boundary smooth. The lowest J on ∂𝒪 (L-BFGS from 200 random directions; 108 reach it to 10⁻³τ): ΔJ_𝒪 = 6.31e-4, ΔJ_𝒪/τ = 6.306 at τ = 0.0001, at the deviation from ϑ* 0.12, 0.09, 0.37; -0.19, 0.02, 1.19; 0.13, -0.14, 0.52 (a = 3.05, -2.23, 2.48, c = -1.43, 0.02, 1.40, s = 0.57, 0.59, 0.61; the fit error 0.0316).

The deterministic flow from 64 points of ∂𝒪 (the minimiser and random directions; (J − J(exact))/τ between 6.31 and 4057.9 there, τ = 0.0001) over 40000 slow units: the largest final G⁻¹-distance from the exact fit 5.37e-5, the largest final (J − J(exact))/τ 4.41e-9.

An exit: the first step (dt = 0.25) after which Σ((ϑ − ϑ*)/L)^p ≥ 1; the runs are integrated together by src/core/gauge_batch.py with each seed's noise as in run_dithered.

| τ | ΔJ_𝒪/τ | exits | the cap (slow units) | the total time | the rate (per slow unit, ± s.e.) | the mean exit time of the exits | the median time | the rate × e^{ΔJ_𝒪/τ} |
|---|---|---|---|---|---|---|---|---|
| 0.00025 | 2.52 | 41/48 | 60000 | 1.094e+06 | 3.75e-5 ± 5.85e-6 | 1.64e+04 | 1.87e+04 | 4.67e-4 |
| 0.0002 | 3.15 | 46/48 | 100000 | 1.805e+06 | 2.55e-5 ± 3.76e-6 | 3.49e+04 | 2.67e+04 | 5.97e-4 |
| 0.00016 | 3.94 | 37/48 | 150000 | 3.836e+06 | 9.64e-6 ± 1.59e-6 | 5.91e+04 | 7.38e+04 | 4.97e-4 |
| 0.00014 | 4.50 | 38/48 | 200000 | 5.328e+06 | 7.13e-6 ± 1.16e-6 | 8.76e+04 | 1.06e+05 | 6.45e-4 |
| 0.000125 | 5.05 | 34/48 | 250000 | 7.631e+06 | 4.46e-6 ± 7.64e-7 | 1.22e+05 | 1.77e+05 | 6.92e-4 |
| 0.00011 | 5.73 | 22/48 | 300000 | 1.149e+07 | 1.91e-6 ± 4.08e-7 | 1.68e+05 | 3e+05 | 5.91e-4 |

Weighted (by the exit counts; var log(rate) ≈ 1/n) least-squares fits of log(rate) = a − c·log τ − B/τ; the mean exit time extrapolated to τ = 10⁻⁴ against 1/(γτ) = 1e+06 at γ = 0.01:

| prefactor τ^(−c) | B (± s.e.) | B/ΔJ_𝒪 (± s.e.) | c | χ² (dof) | the extrapolated mean exit time at τ = 10⁻⁴ | over 1/(γτ) |
|---|---|---|---|---|---|---|
| 0 | 5.74e-4 ± 4.11e-5 | 0.91 ± 0.07 | 0 | 2.38 (4) | 7.74e+05 | 0.77 |
| 0.5 | 6.55e-4 ± 4.11e-5 | 1.04 ± 0.07 | 0.5 | 2.20 (4) | 8.17e+05 | 0.82 |
| 1 | 7.37e-4 ± 4.11e-5 | 1.17 ± 0.07 | 1 | 2.13 (4) | 8.63e+05 | 0.86 |
| free | 7.49e-4 ± 3.51e-4 | 1.19 ± 0.56 | 1.07 ± 2.13 (free) | 2.13 (3) | 8.69e+05 | 0.87 |

## H. The prefactor with the barrier fixed, the lag, the extrapolations' uncertainty, the minima on ∂𝒪 and the exit points (referee round 150; from Part G's stored runs)

With the barrier fixed at ΔJ_𝒪: the weighted (exit counts) least-squares fit of log(rate) + ΔJ_𝒪/τ = a − c·log τ (χ² is quadratic in c, so the profile's Δχ² = 1 and 4 intervals are ĉ ± s.e. and ĉ ± 2 s.e.):

| prefactor τ^(−c) | c (± s.e.) | Δχ² ≤ 1 | Δχ² ≤ 4 | χ² (dof) | p | z = (c − ĉ)/s.e. |
|---|---|---|---|---|---|---|
| fitted | 0.36 ± 0.25 | 0.11 to 0.61 | -0.14 to 0.86 | 2.24 (4) | 0.69 | — |
| 0 | 0 | — | — | 4.31 (5) | 0.506 | -1.44 |
| 0.5 | 0.5 | — | — | 2.56 (5) | 0.767 | +0.57 |
| 1 | 1 | — | — | 8.85 (5) | 0.115 | +2.57 |

The lag: the first t₀ time units of every run discarded (an exit counts after t₀; each run contributes max(t − t₀, 0)); t₀ = 1650, 3300, 6600 is ½, 1, 2 times the slowest relaxation time 1/(3.06·10⁻⁴):

| t₀ | exits by τ | B/ΔJ_𝒪 (± s.e.), c = 0 | c = ½ | c = 1 | ĉ with the barrier ΔJ_𝒪 (± s.e.) |
|---|---|---|---|---|---|
| 0 | 41, 46, 37, 38, 34, 22 | 0.909 ± 0.065 | 1.039 ± 0.065 | 1.169 ± 0.065 | 0.36 ± 0.25 |
| 1650 | 41, 46, 37, 38, 34, 22 | 0.931 ± 0.065 | 1.061 ± 0.065 | 1.191 ± 0.065 | 0.27 ± 0.25 |
| 3300 | 39, 45, 36, 37, 33, 22 | 0.946 ± 0.066 | 1.076 ± 0.066 | 1.205 ± 0.066 | 0.21 ± 0.25 |
| 6600 | 32, 42, 36, 37, 33, 22 | 0.937 ± 0.069 | 1.066 ± 0.069 | 1.194 ± 0.069 | 0.26 ± 0.27 |

The first exit by τ (and the number of exits before 3300): 0.00025: 2135.5 (2); 0.0002: 2204.25 (1); 0.00016: 2614 (1); 0.00014: 2782.75 (1); 0.000125: 3107.25 (1); 0.00011: 22553.2 (0).

The mean exit time extrapolated to τ = 10⁻⁴ (1/(γτ) = 1e+06 at γ = 0.01), with the standard error of the fitted log rate there (the fit's covariance, var log(rate) ≈ 1/n), the range e^{∓s.e.} and the barrier's part alone (10⁴ × its s.e.):

| fit | the mean exit time at τ = 10⁻⁴ | s.e. of its log | the range (one s.e.) | over 1/(γτ) | the barrier's part of the s.e. |
|---|---|---|---|---|---|
| G, c = 0 | 7.74e+05 | 0.17 | 6.55e+05–9.15e+05 | 0.77 (0.65–0.92) | 0.41 |
| G, c = 0.5 | 8.17e+05 | 0.17 | 6.91e+05–9.66e+05 | 0.82 (0.69–0.97) | 0.41 |
| G, c = 1 | 8.63e+05 | 0.17 | 7.3e+05–1.02e+06 | 0.86 (0.73–1.02) | 0.41 |
| G, the barrier ΔJ_𝒪, ĉ = 0.36 | 8e+05 | 0.14 | 6.94e+05–9.22e+05 | 0.80 (0.69–0.92) | — |
| C (the escape set max|α_k| < 6), c = ½ | 6.88e+06 | 0.61 | 3.74e+06–1.27e+07 | 6.88 (3.74–12.65) | 1.09 |

The local minima of J on ∂𝒪 (Part G's 200 L-BFGS runs from the same starts, clustered: (J − J(exact))/τ within 0.001 and (ϑ − ϑ*)/L within 0.05), and the deterministic flow from each over 40000 slow units (1 minus the largest Σ((ϑ − ϑ*)/L)^p after the first step: positive, the flow enters 𝒪 at once and stays in it):

| (J − J(exact))/τ | starts | (ϑ − ϑ*)/L by class (α, γ, κ) | 1 − the largest Σ((ϑ − ϑ*)/L)^p on the flow | the final G⁻¹-distance | the final (J − J(exact))/τ |
|---|---|---|---|---|---|
| 6.306 | 108 | +0.05, +0.11, +0.31; -0.08, +0.03, +0.99; +0.05, -0.17, +0.43 | 9.26e-5 | 3.83e-5 | 2.25e-9 |
| 6.397 | 92 | +0.03, -0.02, -0.14; -0.06, +0.00, -1.00; +0.02, +0.03, -0.18 | 1.49e-4 | 2.04e-5 | 6.39e-10 |

The exit points (the first state outside 𝒪, read from Part G's checkpoint and checked against the stored exit times): the term of Σ((ϑ − ϑ*)/L)^p that dominates there, by (class, coordinate), and the exits through the second width (its term dominant) broadened, on the side of the lowest minimum, and narrowed:

| τ | exits | the dominant term: (class, coordinate) counts | the second width broadened | narrowed | broadened, fraction |
|---|---|---|---|---|---|
| 0.00025 | 41 | (2, κ): 39; (3, κ): 2 | 28 | 11 | 0.68 |
| 0.0002 | 46 | (2, κ): 44; (3, κ): 2 | 30 | 14 | 0.65 |
| 0.00016 | 37 | (2, κ): 37 | 27 | 10 | 0.73 |
| 0.00014 | 38 | (2, κ): 38 | 29 | 9 | 0.76 |
| 0.000125 | 34 | (2, κ): 34 | 27 | 7 | 0.79 |
| 0.00011 | 22 | (2, κ): 22 | 19 | 3 | 0.86 |

Numbers are produced by `experiments/exp_gauge_plateau.py`; nothing is edited by hand.
