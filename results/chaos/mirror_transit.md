# Proposition D‴.48: the rich transit (`exp_mirror_transit.py`)

The free mass–amplitude law (Proposition D‴.45; m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025) from small masses, in the time rescaled by the transit's scale (`mirror_free.iss_path`, `rich_transit`).

## A. From α = 0: the inverse scale space flow on the scale Λ(ε) = (πa_max/(4G))(2/(cε))^{1/2}

Three 6 × 14 Gaussian designs, 3-sparse truths, masses εb_i (b_i log-uniform in [1/4, 4]); the path of Σ|w_i|/√b_i against the loop in τ = t/Λ, u = η/Λ, on 2001 times in [0, 1.5τ*] (τ* the last event): sup_τ|u − u_ISS|, the residual's difference |Φ(w − w_ISS)| over the samples — its rms, its median and the number of samples where it exceeds 10⁻³ (at an event the difference is O(1) over a τ-width ∝ ε^{1/2}, so its L² norm over τ is ∝ ε^{1/4}, while away from the events it is O(ε^{3/2}), the inactive atoms' coefficients; once the width falls below the grid's step the rms counts only the samples that fall within an entry), the largest difference of the entry times (the first sample where |w_i| > 10⁻³ max|w_bp|, on the loop and on the path; in brackets the atoms crossing on one of them only; the grid's step given), the end's distance to weighted basis pursuit.

Design 3: the path's events — enter 8 at 0.7057; enter 13 at 1.356; enter 7 at 2.028; its end against weighted basis pursuit 2.22e-16; the grid's step 0.00152.

| ε | Λ | sup|u − u_ISS| | residual's difference: rms, median, samples > 10⁻³ | entry times: largest difference | |w(1.5τ*) − w_bp|∞ |
|---|---|---|---|---|---|
| 0.01 | 70.25 | 0.536 | 1.17, 1.19, 1941 | 0.931 (4) | 0.368 |
| 0.001 | 222.1 | 0.203 | 0.728, 0.15, 1680 | 0.389 (3) | 0.0171 |
| 0.0001 | 702.5 | 0.0673 | 0.42, 4.59e-4, 871 | 0.14 (0) | 1.84e-4 |
| 1e-05 | 2221 | 0.0216 | 0.238, 3.92e-6, 273 | 0.0441 (0) | 2.82e-6 |
| 1e-06 | 7025 | 0.00685 | 0.133, 9.74e-8, 87 | 0.0137 (0) | 7.06e-8 |
| 1e-07 | 2.221e+04 | 0.00217 | 0.0726, 2.86e-9, 27 | 0.00456 (0) | 2.08e-9 |
| 1e-08 | 7.025e+04 | 6.87e-4 | 0.0402, 8.82e-11, 8 | 0.00152 (0) | 6.42e-11 |
| 1e-10 | 7.025e+05 | 6.87e-5 | 2.87e-4, 8.74e-14, 1 | 0.00152 (0) | 6.36e-14 |
| 1e-12 | 7.025e+06 | 6.70e-6 | 3.96e-6, 7.39e-16, 0 | 0 (0) | 1.11e-16 |

sup|u − u_ISS| in log ε: the slope 0.500 between ε = 10⁻⁶ and 10⁻⁸, 0.503 between 10⁻⁸ and 10⁻¹²; the residual's difference — rms: the slope 0.260 between 10⁻⁶ and 10⁻⁸ (the L² norm's ε^{1/4}: 0.25); median: 1.512 between 10⁻⁶ and 10⁻¹⁰ (ε^{3/2}: 1.5; at 10⁻¹² it nears the rounding floor) (1/Λ ∝ ε^{1/2}: 0.5). Each atom is integrated on its own curve (α alone, β = εb + 2c sinh²(α/a_max)), so Q is kept exactly.

Design 4: the path's events — enter 0 at 1.226; enter 4 at 1.453; enter 2 at 19.16; enter 11 at 38.43; enter 10 at 39.41; enter 3 at 55.48; exit 10 at 55.48; enter 7 at 78.31; its end against weighted basis pursuit 4.54e-16; the grid's step 0.0587.

| ε | Λ | sup|u − u_ISS| | residual's difference: rms, median, samples > 10⁻³ | entry times: largest difference | |w(1.5τ*) − w_bp|∞ |
|---|---|---|---|---|---|
| 0.01 | 70.25 | 0.61 | 0.172, 0.0217, 1996 | 57.6 (3) | 0.02 |
| 0.001 | 222.1 | 0.235 | 0.106, 0.00693, 1837 | 42.9 (0) | 0.00345 |
| 0.0001 | 702.5 | 0.0865 | 0.065, 0.00262, 1683 | 37.5 (0) | 5.04e-4 |
| 1e-05 | 2221 | 0.0284 | 0.0498, 5.72e-4, 718 | 11.6 (0) | 6.62e-5 |
| 1e-06 | 7025 | 0.00892 | 0.0154, 1.33e-5, 254 | 3.76 (0) | 4.92e-6 |
| 1e-07 | 2.221e+04 | 0.00282 | 0.00125, 2.46e-7, 88 | 1.17 (0) | 1.85e-7 |
| 1e-08 | 7.025e+04 | 8.91e-4 | 6.23e-4, 3.11e-9, 28 | 0.352 (0) | 6.15e-9 |
| 1e-10 | 7.025e+05 | 8.90e-5 | 2.19e-4, 2.23e-12, 3 | 0.0587 (0) | 6.28e-12 |
| 1e-12 | 7.025e+06 | 8.81e-6 | 1.33e-4, 9.24e-15, 1 | 0.0587 (0) | 6.09e-15 |

sup|u − u_ISS| in log ε: the slope 0.500 between ε = 10⁻⁶ and 10⁻⁸, 0.501 between 10⁻⁸ and 10⁻¹²; the residual's difference — rms: the slope 0.696 between 10⁻⁶ and 10⁻⁸ (the L² norm's ε^{1/4}: 0.25); median: 1.694 between 10⁻⁶ and 10⁻¹⁰ (ε^{3/2}: 1.5; at 10⁻¹² it nears the rounding floor) (1/Λ ∝ ε^{1/2}: 0.5). Each atom is integrated on its own curve (α alone, β = εb + 2c sinh²(α/a_max)), so Q is kept exactly.

Design 5: the path's events — enter 0 at 2.413; enter 9 at 2.556; enter 1 at 2.988; enter 2 at 3.18; exit 0 at 3.18; enter 12 at 4.348; enter 3 at 4.615; exit 9 at 4.615; enter 8 at 12.15; enter 10 at 58.54; its end against weighted basis pursuit 3.33e-16; the grid's step 0.0439.

| ε | Λ | sup|u − u_ISS| | residual's difference: rms, median, samples > 10⁻³ | entry times: largest difference | |w(1.5τ*) − w_bp|∞ |
|---|---|---|---|---|---|
| 0.01 | 70.25 | 1.15 | 0.235, 0.0203, 1331 | 55.6 (5) | 0.178 |
| 0.001 | 222.1 | 0.5 | 0.147, 0.0149, 1289 | 53 (0) | 0.00341 |
| 0.0001 | 702.5 | 0.166 | 0.0847, 5.17e-4, 841 | 26.7 (0) | 4.37e-6 |
| 1e-05 | 2221 | 0.0534 | 0.0453, 6.99e-6, 284 | 8.43 (0) | 2.17e-7 |
| 1e-06 | 7025 | 0.017 | 0.0169, 1.12e-7, 96 | 2.68 (0) | 8.00e-9 |
| 1e-07 | 2.221e+04 | 0.00537 | 0.0101, 2.65e-9, 30 | 0.834 (0) | 2.65e-10 |
| 1e-08 | 7.025e+04 | 0.0017 | 0.00484, 8.61e-11, 11 | 0.263 (0) | 8.48e-12 |
| 1e-10 | 7.025e+05 | 1.70e-4 | 1.67e-4, 1.66e-13, 1 | 0.0439 (0) | 8.22e-15 |
| 1e-12 | 7.025e+06 | 1.97e-5 | 2.20e-8, 2.63e-15, 0 | 0 (0) | 4.51e-16 |

sup|u − u_ISS| in log ε: the slope 0.500 between ε = 10⁻⁶ and 10⁻⁸, 0.484 between 10⁻⁸ and 10⁻¹²; the residual's difference — rms: the slope 0.271 between 10⁻⁶ and 10⁻⁸ (the L² norm's ε^{1/4}: 0.25); median: 1.458 between 10⁻⁶ and 10⁻¹⁰ (ε^{3/2}: 1.5; at 10⁻¹² it nears the rounding floor) (1/Λ ∝ ε^{1/2}: 0.5). Each atom is integrated on its own curve (α alone, β = εb + 2c sinh²(α/a_max)), so Q is kept exactly.

## B. On branches: the logarithmic scale

Design 3 of A; amplitudes s_i·U(0.2, 0.6) (the truth's signs on its atoms, random elsewhere), masses ε; the one-sided path (frozen signs, weights 1/(γ|a_i(0)|)) against the loop in τ = t/log(1/ε), u = (η − η(0))/log(1/ε). The path's events: enter 8 at 21.77; enter 13 at 35.26; enter 7 at 92.15; its end against the sign-constrained weighted basis pursuit 4.44e-16.

| ε | all on branches | sup|u − u_ISS| | residual: rms difference | entry times: largest difference | |w(1.5τ*) − w_lim|∞ |
|---|---|---|---|---|---|
| 0.0001 | yes | 15.7 | 0.652 | 55 (0) | 1.26e-5 |
| 1e-08 | yes | 7.97 | 0.465 | 27.9 (0) | 3.59e-9 |
| 1e-16 | yes | 3.98 | 0.328 | 14 (0) | 4.44e-15 |
| 1e-32 | yes | 1.99 | 0.232 | 6.98 (0) | 6.66e-16 |
| 1e-64 | yes | 0.996 | 0.164 | 3.52 (0) | 7.77e-16 |

A ± pair per atom (amplitudes ±0.4, equal masses, w(0) = u): the atom's common increment over log(1/u) against the two-sided path with the weights 1/(γa₀):

| u | sup|u − u_ISS| | residual: rms difference | |W(1.5τ*) − w_end|∞ | entries on the path |
|---|---|---|---|---|
| 0.0001 | 20.8 | 1.49 | 0.0396 | 3 |
| 1e-08 | 7.86 | 0.998 | 0.00712 | 3 |
| 1e-16 | 3.96 | 0.705 | 2.32e-4 | 3 |
| 1e-32 | 1.98 | 0.493 | 1.94e-6 | 3 |
| 1e-64 | 0.991 | 0.347 | 2.78e-10 | 3 |

## C. Moving centres

Gaussian kernels of width 0.12 on 61 samples of [−1, 1]; 21 kernels started on a grid of [−1, 1] (spacing 0.1), α = 0, masses ε; the target 3 bumps of the same width at [-0.52, 0.07, 0.61] with amplitudes [0.9, -0.6, 0.75]. A kernel has formed when |w| > 0.01·max|A*|.

(1) The paradigm's metric (ċ = Gβa⟨r∂_cφ⟩). The sliding inverse scale space: enter kernel 5 (start -0.5) at τ = 10.72; enter kernel 16 (start +0.6) at τ = 12.85; enter kernel 11 (start +0.1) at τ = 16.39; it ends with 3 kernels at [-0.52, 0.07, 0.61] with coefficients [0.9, -0.6, 0.75] (rms residual 1.36e-16). The loop (the entries at the first sample — step 0.0484 in τ — where |w| > 0.01·max|A*|):

| ε | Λ | entries (kernel: τ) | formed | their centres | their coefficients | rms residual | inactive: max |Δc|, max |w| |
|---|---|---|---|---|---|---|---|
| 0.001 | 222.1 | 5: 10.7, 4: 11.23, 6: 12.25, 16: 12.83, 17: 13.7, 15: 14.52, 11: 16.22, 10: 16.8 | 8 | [-0.52, -0.52, -0.5202, 0.07, 0.07, 0.6103, 0.6101, 0.6101] | [0.2103, 0.6566, 0.0323, -0.1823, -0.3884, 0.0382, 0.5812, 0.1299] | 0.00702 | 0.181, 0.00714 |
| 0.0001 | 702.5 | 5: 10.75, 16: 12.88, 11: 16.36, 10: 17.53 | 4 | [-0.52, 0.07, 0.07, 0.61] | [0.8969, -0.0127, -0.5872, 0.7493] | 1.21e-6 | 0.0877, 0.00312 |
| 1e-05 | 2221 | 5: 10.75, 16: 12.88, 11: 16.41 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 3.80e-7 | 0.0331, 2.91e-5 |
| 1e-06 | 7025 | 5: 10.75, 16: 12.88, 11: 16.41 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 1.55e-8 | 0.0026, 4.26e-7 |
| 1e-07 | 2.221e+04 | 5: 10.75, 16: 12.88, 11: 16.41 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 4.40e-10 | 2.34e-4, 1.12e-8 |
| 1e-08 | 7.025e+04 | 5: 10.75, 16: 12.88, 11: 16.41 | 3 | [-0.52, 0.07, 0.61] | [0.9, -0.6, 0.75] | 1.33e-11 | 2.26e-5, 3.37e-10 |

(2) The Wasserstein–Fisher–Rao metric (ċ = Ga⟨r∂_cφ⟩). The inactive kernels' limit (du/dτ = ⟨f*φ⟩, dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ⟩) gathers the 21 kernels at 3 positions by τ = 10 [-0.521, 0.068, 0.613]. The loop against that limit along the trajectory (the largest |c − c_lim| over τ < 0.9τ_first, τ_first the first time a kernel's |w| exceeds the threshold), and its end:

| ε | τ_first | max |c − c_lim| before | max |w| before | formed | their positions (clusters) | rms residual |
|---|---|---|---|---|---|---|
| 0.0001 | 10.67 | 1.67e-5 | 8.69e-5 | 15 | [-0.52, 0.07, 0.61] (3) | 5.85e-10 |
| 1e-05 | 10.65 | 3.43e-7 | 3.08e-6 | 12 | [-0.52, 0.07, 0.61] (3) | 1.25e-10 |
| 1e-06 | 10.65 | 2.59e-7 | 9.85e-8 | 8 | [-0.52, 0.07, 0.61] (3) | 1.05e-10 |
| 1e-08 | 10.65 | 1.96e-7 | 9.86e-11 | 3 | [-0.52, 0.07, 0.61] (3) | 3.37e-11 |
| 1e-10 | 10.65 | 2.02e-7 | 9.86e-14 | 3 | [-0.52, 0.07, 0.61] (3) | 1.77e-13 |

Unequal masses (εb_i, b_i log-uniform in [0.264, 3.82], ε = 1e-06): the limit with ρ_i = b_i^{−1/2} followed to 1.13e-7 before τ_first = 8.254 (the limit with the masses ignored, ρ = 1: 0.447); formed 3 kernels in 3 clusters.

(3) The centres held on the grid (ε = 1e-06, to τ = 63.9): 5 kernels formed, rms residual 0.0192 (the grid's least-squares fit: 8.38e-5).

(4) 16 random targets (three bumps of the kernels' width, centres in [−0.8, 0.8] at least 0.3 apart, amplitudes ±U(0.4, 1)): the sliding inverse scale space ends with exactly three kernels within 10⁻³ of the bumps' centres on 16 of 16 (0 exits in all); the loop with the paradigm's metric on 11, 13 and 16 of 16 at ε = 1e-05, 1e-06, 1e-08 (formed kernels [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5], [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4], [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]; the inactive kernels' largest displacement 0.12, 0.0699, 0.00553), the formed kernels those of the sliding path, in its order, on 15 of 16 at ε = 1e-08; with the Wasserstein–Fisher–Rao metric ε = 1e-06: [4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7] kernels, exactly three at the bumps on 0 of 16; ε = 1e-08: [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5] kernels, exactly three at the bumps on 10 of 16; ε = 1e-10: [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3] kernels, exactly three at the bumps on 16 of 16; with the centres held on the grid [4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6] kernels (rms residual 0.0114–0.0257).

(5) Ties, near ties, close bumps off the grid, wider bumps (the paradigm's metric): the sliding inverse scale space (its entries; a tie enters together), and the loop at each ε:

| target | sliding ISS: entries (τ, kernel), kernels, centres, coefficients, rms residual | ε | loop: entries (kernel: τ), formed, centres, coefficients (clusters), rms residual |
|---|---|---|---|
| one bump midway between two kernels (a tie) | [(12.48, 10), (12.48, 11)], 2, [0.05, 0.05], [0.4, 0.4], 6.72e-17 | 1e-06 | 10: 12.47, 11: 12.47; 2, [0.05, 0.05], [0.4, 0.4] (1), 1.12e-9 |
|  |  | 1e-08 | 10: 12.48, 11: 12.48; 2, [0.05, 0.05], [0.4, 0.4] (1), 1.09e-12 |
|  |  | 1e-10 | 10: 12.48, 11: 12.48; 2, [0.05, 0.05], [0.4, 0.4] (1), 1.09e-15 |
| one bump 0.001 off the midpoint (a near tie) | [(12.458, 11)], 1, [0.051], [0.8], 4.08e-17 | 1e-06 | 11: 12.45, 10: 12.5; 2, [0.051, 0.051], [0.062, 0.738] (1), 1.11e-9 |
|  |  | 1e-08 | 11: 12.46; 1, [0.051], [0.8] (1), 5.55e-8 |
|  |  | 1e-10 | 11: 12.46; 1, [0.051], [0.8] (1), 4.83e-11 |
| three bumps, one at a midpoint | [(10.734, 5), (12.822, 16), (16.822, 10)], 3, [-0.52, 0.05, 0.61], [0.9, -0.6, 0.75], 1.57e-16 | 1e-06 | 5: 10.74, 16: 12.84, 10: 16.81, 11: 16.84; 4, [-0.52, 0.05, 0.05, 0.61], [0.9, -0.446, -0.154, 0.75] (3), 7.66e-9 |
|  |  | 1e-08 | 5: 10.74, 16: 12.82, 10: 16.82; 3, [-0.52, 0.05, 0.61], [0.9, -0.6, 0.75] (3), 8.86e-8 |
|  |  | 1e-10 | 5: 10.74, 16: 12.82, 10: 16.82; 3, [-0.52, 0.05, 0.61], [0.9, -0.6, 0.75] (3), 6.53e-10 |
| two bumps 0.15 apart, off the grid | [(7.618, 10), (33.496, 11)], 2, [-0.03, 0.12], [0.8, 0.6], 8.06e-10 | 1e-06 | 10: 7.62, 11: 16.18; 2, [-0.03, 0.12], [0.8, 0.6] (2), 1.97e-9 |
|  |  | 1e-08 | 10: 7.62, 11: 29.59; 2, [-0.03, 0.12], [0.8, 0.6] (2), 3.07e-12 |
|  |  | 1e-10 | 10: 7.62, 11: 33.05; 2, [-0.03, 0.12], [0.8, 0.6] (2), 5.76e-15 |
| two bumps 0.25 apart, off the grid | [(8.879, 10), (17.065, 11)], 2, [-0.07, 0.18], [0.8, 0.6], 3.57e-14 | 1e-06 | 10: 8.886, 11: 12.06; 2, [-0.07, 0.18], [0.8, 0.6] (2), 1.11e-9 |
|  |  | 1e-08 | 10: 8.886, 11: 16.05; 2, [-0.07, 0.18], [0.8, 0.6] (2), 3.44e-12 |
|  |  | 1e-10 | 10: 8.886, 11: 16.95; 2, [-0.07, 0.18], [0.8, 0.6] (2), 5.55e-15 |
| three wider bumps (width 0.2) | 21 entries, at τ = 9.1, 11.2, 15.7, 18.8, 24.7, 86.5, 198, 514, 1.6e+03, 4.28e+03, 1.27e+05, 8.61e+05, 4.91e+06, 6.8e+06, 8.56e+06, 1.4e+08, 1.64e+08, 2.29e+08, 3.04e+08, 1.22e+10, 1.23e+10 (not terminated: 21 kernels active, 13 with |w| above the threshold, rms residual 3.32e-5) | 1e-05 | 5: 9.185, 4: 9.952, 16: 11.22, 17: 12.49, 11: 15.67, 10: 26.74, 1: 229.8, 20: 634.2, 19: 1142, 0: 1305, 12: 1.309e+05, 9: 1.327e+05, 2: 2.205e+05, 18: 6.19e+07, 14: 7.135e+08, 7: 1.468e+09; 16, (16), 2.37e-7 |

The target 'three wider bumps (width 0.2)' is not three kernels of width 0.12: the sliding path enters a kernel at a time, at times growing geometrically, its rms residual at the events 0.505, 0.37, 0.242, 0.132, 0.0831, 0.0356, 0.0238, 0.0121, 0.00259, 9.54e-4, 6.93e-4, 1.90e-4, 1.90e-4, 1.90e-4, 1.90e-4, 4.87e-5, 3.75e-5, 3.63e-5, 3.46e-5, 3.39e-5, 3.35e-5, 3.32e-5 (at τ = 0, 9.1, 11.2, 15.7, 18.8, 24.7, 86.5, 198, 514, 1.6e+03, 4.28e+03, 1.27e+05, 8.61e+05, 4.91e+06, 6.8e+06, 8.56e+06, 1.4e+08, 1.64e+08, 2.29e+08, 3.04e+08, 1.22e+10, 1.23e+10); 11 of its 21 relaxations became stationary within the budget. The loop at ε = 1e-05 (on a geometric τ grid to τ = 1.6e+10) enters 16 kernels, at τ = 9.19 to 1.47e+09: the first 1 of its entrants the path's in its order, the first 10 the path's first 10 (its entries at τ = 9.19, 9.95, 11.2, 12.5, 15.7, 26.7, 230, 634, 1.14e+03, 1.3e+03; the path's 9.1, 11.2, 15.7, 18.8, 24.7, 86.5, 198, 514, 1.6e+03, 4.28e+03); it ends with 16 kernels above the threshold at the rms residual 2.37e-7 (the path 3.32e-5).


(6) A lost near-tie (target 10 of (4)): the sliding path enters kernels [2, 14, 12] at τ = [11.492, 17.494, 21.386]. After the second entry an inactive kernel left near its wall lingers at a small margin μ and creeps; its instantaneous speed dc/dτ (the loop's own right-hand side) against the core's profile (πa_max²/(4c))εb^{3/2}cos(πμ/2)/sin³(πμ/2)·s⟨r∂_cφ⟩ and against the bound's μ^{−3} form (πa_max²/(4c))εμ^{−3}|⟨r∂_cφ⟩|:

| ε | kernel | τ | μ | ⟨r∂_cφ⟩ | dc/dτ | the profile's speed (ratio) | the μ^{−3} form (ratio) |
|---|---|---|---|---|---|---|---|
| 1e-08 | 3 | 17.72 | 0.00434 | -9.28e-6 | -9.14e-6 | -9.22e-6 (0.9915) | 3.57e-5 (0.2558) |
| 1e-08 | 3 | 18.58 | 0.00434 | -9.28e-6 | -9.14e-6 | -9.22e-6 (0.9915) | 3.57e-5 (0.2558) |
| 1e-08 | 3 | 19.46 | 0.00434 | -9.27e-6 | -9.14e-6 | -9.21e-6 (0.9915) | 3.57e-5 (0.2558) |
| 1e-08 | 3 | 20.35 | 0.00434 | -9.27e-6 | -9.13e-6 | -9.21e-6 (0.9915) | 3.57e-5 (0.2558) |
| 1e-08 | 15 | 17.72 | 0.00646 | -0.033 | 0.00988 | 0.00992 (0.9961) | 0.0385 (0.2570) |
| 1e-08 | 15 | 18.58 | 0.00549 | -0.0364 | 0.0178 | 0.0179 (0.9946) | 0.0693 (0.2566) |
| 1e-08 | 15 | 19.46 | 0.00392 | -0.0395 | 0.0525 | 0.0531 (0.9896) | 0.206 (0.2553) |
| 1e-08 | 15 | 20.35 | 4.49e-4 | 2.21e-5 | -0.0109 | -0.0197 (0.5515) | 0.0765 (0.1423) |
| 1e-10 | 3 | 17.72 | 0.00484 | -1.06e-5 | -7.59e-8 | -7.60e-8 (0.9999) | 2.94e-7 (0.2580) |
| 1e-10 | 3 | 18.89 | 0.00484 | -1.06e-5 | -7.59e-8 | -7.59e-8 (0.9999) | 2.94e-7 (0.2580) |
| 1e-10 | 3 | 20.09 | 0.00484 | -1.06e-5 | -7.59e-8 | -7.59e-8 (0.9999) | 2.94e-7 (0.2580) |
| 1e-10 | 3 | 21.29 | 0.00484 | -1.06e-5 | -7.59e-8 | -7.59e-8 (0.9999) | 2.94e-7 (0.2580) |
| 1e-10 | 15 | 17.72 | 0.00696 | -0.0326 | 7.84e-5 | 7.84e-5 (1.0000) | 3.04e-4 (0.2580) |
| 1e-10 | 15 | 18.89 | 0.00587 | -0.0327 | 1.31e-4 | 1.31e-4 (1.0000) | 5.08e-4 (0.2580) |
| 1e-10 | 15 | 20.09 | 0.00474 | -0.0327 | 2.49e-4 | 2.49e-4 (0.9999) | 9.66e-4 (0.2580) |
| 1e-10 | 15 | 21.29 | 0.0036 | -0.0329 | 5.72e-4 | 5.72e-4 (0.9999) | 0.00222 (0.2580) |

| ε | kernel: displacement over the window (since τ_2 + 0.2), since its start | the loop's entrants | the formed kernels' centres |
|---|---|---|---|
| 1e-08 | 3: -2.40e-5, -0.0038; 15: 0.161, 0.162 | [2, 14, 15] | [-0.7519, 0.1116, 0.4879] |
| 1e-10 | 3: -2.71e-7, -3.19e-5; 15: 7.94e-4, 8.02e-4 | [2, 14, 12] | [-0.7519, 0.1116, 0.4879] |

## D. Moving widths (σ̇ = Gβa⟨r∂_σφ⟩ beside ċ = Gβa⟨r∂_cφ⟩; kernels started at width 0.12)

The target: the three bumps of C at widths [0.2, 0.08]; τ ≤ 60.

| the bumps' width | ε | widths | entrants | formed | centres | widths | coefficients | rms residual |
|---|---|---|---|---|---|---|---|---|
| 0.2 | 1e-05 | moving | [5, 16, 4, 11] | 4 | [-0.658, -0.5119, 0.07, 0.6095] | [0.1837, 0.1972, 0.1991, 0.199] | [0.0586, 0.8562, -0.5988, 0.7481] | 6.64e-5 |
| 0.2 | 1e-05 | fixed | [5, 4, 16, 17, 11, 10] | 6 | [-0.6729, -0.4443, -0.0059, 0.1443, 0.5485, 0.7702] | [0.12, 0.12, 0.12, 0.12, 0.12, 0.12] | [0.5714, 0.7547, -0.366, -0.3404, 0.6221, 0.4448] | 0.0239 |
| 0.2 | 1e-06 | moving | [5, 16, 11] | 3 | [-0.52, 0.07, 0.61] | [0.2, 0.2, 0.2] | [0.9, -0.6, 0.75] | 2.33e-8 |
| 0.2 | 1e-06 | fixed | [5, 4, 16, 17, 11, 10] | 6 | [-0.6731, -0.4445, -0.0032, 0.1475, 0.5481, 0.7699] | [0.12, 0.12, 0.12, 0.12, 0.12, 0.12] | [0.5705, 0.7552, -0.3773, -0.3293, 0.6213, 0.446] | 0.0238 |
| 0.08 | 1e-05 | moving | [5, 16, 11] | 3 | [-0.52, 0.07, 0.61] | [0.08, 0.08, 0.08] | [0.9, -0.6, 0.75] | 5.74e-7 |
| 0.08 | 1e-05 | fixed | [5, 16, 11] | 3 | [-0.5194, 0.0713, 0.6083] | [0.12, 0.12, 0.12] | [0.707, -0.4756, 0.591] | 0.0943 |
| 0.08 | 1e-06 | moving | [5, 16, 11] | 3 | [-0.52, 0.07, 0.61] | [0.08, 0.08, 0.08] | [0.9, -0.6, 0.75] | 1.19e-8 |
| 0.08 | 1e-06 | fixed | [5, 16, 11] | 3 | [-0.5194, 0.0713, 0.6083] | [0.12, 0.12, 0.12] | [0.707, -0.4756, 0.591] | 0.0943 |

Numbers are produced by `experiments/exp_mirror_transit.py`; nothing is edited by hand.
