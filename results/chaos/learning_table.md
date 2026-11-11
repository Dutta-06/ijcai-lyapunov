# Theorem D‴: uniform-in-time propagation of chaos in the learning regime — the averaged slow loop (`exp_chaos_learning.py`)

n = 1, saturated features with a_max = 4.0, c_max = 3.0, s ∈ [0.3, 2.0], prior ϑ⁰ = 0, G = diag(1.0, 0.25, 0.25); plant f(x) = Σ_k A_k e^{−(x−C_k)²/(2S_k²)} with (A, C, S) = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)] (three bumps; the network can represent it exactly with three kernels of weights 3A_k); reference x_m = 2 sin t, ρ its occupation measure on a grid of 256 times; slow time s = t/k_e, (𝔄_m) integrated by the midpoint rule with step 0.02; the m-systems share their particles with an M = 2048 reference (the stand-in for the mean-field flow, as in Theorems D′ and D″); horizon 900 slow units.

Constants. Proposition K: the collapsed critical point (best single kernel, a = 0.331, c = -0.057, s = 1.218) is a strict local minimum of J_m iff σ₀ > σ_crit = 0.04716; the runs use σ₀ = 0.15σ_crit = 0.00707 (learning regime; the leak time 1/σ₀ = 141 slow units). Whitened Lipschitz constants by search (|x| ≤ 2.5, 6·10⁴ sampled ϑ, Nelder–Mead polish): L_ϑ^G = 4.044 (closed-form bound 4.658), L₂^G = 11.11 (closed-form bound λ_G a_max/s_min² = 44.44); K_Φ = 4.0, F_ρ = 0.999.

## 1. Hypothesis (H) at the three-atom stationary measure, and the window

σ₀ = 0.15σ_crit = 0.00707: the critical point of J₃ from the exact representation has atoms (a, c, s) = (2.4504, -1.4161, 0.8287), (-2.2879, 0.0120, 0.7009), (2.0478, 1.3427, 0.8810); sup-norm of ∇J₃ = 4.2e-17 (a stationary measure of the mean-field flow); fit error ‖f − f̂‖/‖f‖ = 0.2439; eigenvalues of G^{1/2}∇²J₃G^{1/2}: 0.0057, 0.0086, 0.0205, 0.0216, 0.0362, 0.0610, 0.0736, 0.1357, 0.2003 (a strict local minimum of J₃). κ_R at the atoms: 0.00004, -0.00227, 0.00411 (σ₀ − max = +0.00297); over the spheres of radius ϱ′ around the atoms (400 directions, radii ϱ′ and ϱ′/2): ϱ′ = 0.05: κ_R = 0.00461, margin σ₀ − κ_R = +0.00247; ϱ′ = 0.1: κ_R = 0.00512, margin σ₀ − κ_R = +0.00195; ϱ′ = 0.2: κ_R = 0.00619, margin σ₀ − κ_R = +0.00089; ϱ′ = 0.3: κ_R = 0.00730, margin σ₀ − κ_R = -0.00023.

For comparison, σ₀ = 0.3σ_crit = 0.01415: the critical point of J₃ from the exact representation has atoms (a, c, s) = (2.1082, -1.3689, 0.9251), (-2.1262, -0.0059, 0.7462), (1.7063, 1.1865, 1.0875); sup-norm of ∇J₃ = 9.7e-17 (a stationary measure of the mean-field flow); fit error ‖f − f̂‖/‖f‖ = 0.4057; eigenvalues of G^{1/2}∇²J₃G^{1/2}: 0.0032, 0.0099, 0.0274, 0.0375, 0.0498, 0.0597, 0.0753, 0.1768, 0.2241 (a strict local minimum of J₃). κ_R at the atoms: 0.00841, -0.00220, 0.01696 (σ₀ − max = -0.00282); over the spheres of radius ϱ′ around the atoms (400 directions, radii ϱ′ and ϱ′/2): ϱ′ = 0.05: κ_R = 0.01757, margin σ₀ − κ_R = -0.00342; ϱ′ = 0.1: κ_R = 0.01819, margin σ₀ − κ_R = -0.00404; ϱ′ = 0.2: κ_R = 0.01949, margin σ₀ − κ_R = -0.00534; ϱ′ = 0.3: κ_R = 0.02088, margin σ₀ − κ_R = -0.00673.

At σ₀ = 0.15σ_crit the window of Remark D‴.1 is realised: κ_R[μ*](N_0.1) = 0.00512 < σ₀ = 0.00707 < σ_crit = 0.04716: the collapsed configuration is unstable (the kernels separate) and the three-atom measure satisfies the strict second-order condition with margin κ(ϱ′ = 0.1) = 0.00195 (relaxation time 1/κ = 512 slow units). At σ₀ = 0.3σ_crit the three-atom measure fails the splitting condition at one atom (κ_R = 0.01696 > σ₀ = 0.01415): there the mean-field flow's limit has more atoms, and the hypothesis has to be checked at that configuration.

Theorem D‴'s constants for the warm start (Λ₁ = 0, s₁ = 0, κ := ½ min{κ(ϱ′), min_s a(s)} = 0.00070 with a(s) measured along the reference below): R₀ = 4.561, R_Θ = 2857.326, V_max = 40.428, ϱ = 1.064e-06, θ = 4.574e-11, Δ = 6.995e-14, C₁ = 4(L_ϑ^G)²K_Φ²/κ² = 2.163e+09, m₀ = 16K_Φ²/θ² = 1.224e+23. The bound C₁/m is uniform in time but far above the measured distances: C₁ carries the worst-case constants L_ϑ^G and K_Φ and the factor κ⁻², and the tail term needs m ≥ m₀.

## 2. The hypothesis along the reference flow

κ_R at the reference particles (max over the 2048 particles) and over the spheres of radius ϱ′ = 0.1 (G⁻¹ metric) around 256 of them (24 directions each; a lower estimate of the supremum over the neighbourhood), computed with the reference's own residual (which carries its sampling fluctuation, the term L₂^GY of (2.1)) and, for the warm start, with the residual r* of the stationary measure μ* (the proxy for the mean-field flow's residual: the hypothesis (H) is about μ_s); the margin a(s) = σ₀ − κ_R[μ*](spheres); ‖r_ref − r*‖_{L¹(ρ)}; the reference's fit error ‖f − f̂‖/‖f‖ on the orbit and its rms drift (G⁻¹ norm); means over the seeds.

| s | warm: κ_R particles (r_ref) | warm: κ_R spheres (r_ref) | warm: κ_R spheres (r*) | warm: a(s) | warm: ‖r_ref − r*‖_{L¹} | warm: fit | warm: rms drift | cold: κ_R particles | cold: κ_R spheres | cold: σ₀ − κ_R spheres | cold: fit | cold: rms drift |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.00486 | 0.00576 | 0.00566 | +0.00141 | 1.73e-02 | 0.2473 | 5.97e-03 | 0.21869 | 0.18640 | -0.17932 | 0.9478 | 1.47e-01 |
| 10 | 0.00476 | 0.00569 | 0.00568 | +0.00139 | 8.66e-03 | 0.2459 | 2.07e-03 | 0.10047 | 0.09501 | -0.08794 | 0.6416 | 7.58e-02 |
| 25 | 0.00473 | 0.00568 | 0.00566 | +0.00142 | 5.02e-03 | 0.2452 | 8.44e-04 | 0.06278 | 0.05967 | -0.05260 | 0.4321 | 3.78e-02 |
| 50 | 0.00468 | 0.00563 | 0.00561 | +0.00147 | 3.57e-03 | 0.2449 | 3.56e-04 | 0.04369 | 0.03862 | -0.03154 | 0.3093 | 1.79e-02 |
| 100 | 0.00463 | 0.00559 | 0.00556 | +0.00151 | 3.01e-03 | 0.2449 | 1.54e-04 | 0.03053 | 0.02472 | -0.01765 | 0.2592 | 7.44e-03 |
| 150 | 0.00461 | 0.00557 | 0.00554 | +0.00154 | 2.90e-03 | 0.2449 | 1.01e-04 | 0.02116 | 0.01795 | -0.01088 | 0.2511 | 4.30e-03 |
| 225 | 0.00459 | 0.00554 | 0.00552 | +0.00156 | 2.85e-03 | 0.2450 | 6.12e-05 | 0.01580 | 0.01352 | -0.00645 | 0.2478 | 2.57e-03 |
| 300 | 0.00456 | 0.00553 | 0.00550 | +0.00158 | 2.82e-03 | 0.2450 | 3.96e-05 | 0.01279 | 0.01143 | -0.00436 | 0.2464 | 1.82e-03 |
| 450 | 0.00452 | 0.00549 | 0.00546 | +0.00161 | 2.80e-03 | 0.2451 | 1.84e-05 | 0.01095 | 0.01069 | -0.00362 | 0.2450 | 1.12e-03 |
| 600 | 0.00449 | 0.00547 | 0.00543 | +0.00164 | 2.79e-03 | 0.2451 | 9.36e-06 | 0.01097 | 0.01025 | -0.00318 | 0.2444 | 8.06e-04 |
| 750 | 0.00447 | 0.00545 | 0.00541 | +0.00166 | 2.79e-03 | 0.2451 | 5.22e-06 | 0.01084 | 0.00958 | -0.00250 | 0.2440 | 5.99e-04 |
| 900 | 0.00445 | 0.00543 | 0.00540 | +0.00167 | 2.79e-03 | 0.2451 | 3.16e-06 | 0.01089 | 0.00959 | -0.00251 | 0.2438 | 4.71e-04 |

Warm start: min_s a(s) = +0.00139 > 0: (H) holds along the whole run with the stationary residual (s₁ = 0, Λ₁ = 0). With the reference's own residual the sphere statistic is larger by at most 0.00013 (its sampling fluctuation; the proof's bound (0.2) allows L₂^G‖r_ref − r*‖_{L¹} ≤ 0.2140).
Cold start: min_s a(s) = -0.19609 at s = 0, final a = -0.00251: (H) fails at the end of the run (the connected support reaches the ridges, Remark D‴.2).

## 3. Warm start near the stationary measure (μ₀ = the mixture of Gaussians of width 0.02 around the atoms; (H) from s = 0)

| m | max D_m on (0, 112.5] | max D_m on (112.5, 225] | max D_m on (225, 450] | max D_m on (450, 900] | time average on (450, 900] | bound C₁/m | mean ‖f̂_Θ − f̂_ref‖ on the last window | max_i ‖δ_i‖_{G⁻¹} at s = 900 (3ϱ′/4 = 0.075) | fraction of particles with ‖δ_i‖_{G⁻¹} > 1 at s = 900 |
|---|---|---|---|---|---|---|---|---|---|
| 16 | 3.957e-01 ± 2.4e-01 | 4.230e-01 ± 2.6e-01 | 4.410e-01 ± 2.8e-01 | 4.484e-01 ± 2.8e-01 | 4.459e-01 ± 2.8e-01 | 1.352e+08 | 3.073e-02 ± 9.7e-03 | 1.361e+00 | 0.188 ± 0.188 |
| 64 | 1.224e-01 ± 4.6e-02 | 1.308e-01 ± 4.9e-02 | 1.364e-01 ± 5.0e-02 | 1.383e-01 ± 5.0e-02 | 1.377e-01 ± 5.1e-02 | 3.379e+07 | 1.542e-02 ± 3.3e-03 | 6.197e-01 | 0.000 ± 0.000 |
| 256 | 1.193e-02 ± 5.0e-03 | 1.273e-02 ± 5.2e-03 | 1.337e-02 ± 5.5e-03 | 1.368e-02 ± 5.6e-03 | 1.358e-02 ± 5.5e-03 | 8.448e+06 | 4.674e-03 ± 1.2e-03 | 2.144e-01 | 0.000 ± 0.000 |

Slope of the last-window maximum of D_m against m (per-seed fits, mean ± s.e.): -1.28 ± 0.30 (Theorem D‴: −1 for the mean square; the function-level distance ‖f̂_Θ − f̂_ref‖ has slope -0.69 ± 0.16, m^{−1/2} being −0.5). Ratio of the last to the first window maxima: 1.13, 1.13, 1.15; last to the one before: 1.02, 1.01, 1.02 (uniform in time means no growth).

Proposition D‴.5 (the weight-sensitivity law). The sensitivity S_k = ∂ϑ*_k/∂w of the three-atom stationary point to the cluster weights (central differences of the weighted stationary point, whitened norms): ‖S_k‖ = 8.377, 4.657, 8.322. Predicted mean-square stationary distance E D_∞ = (1/m − 1/M)·Σ_k w̄_k tr(S̃_k C̃ S̃_kᵀ), measured D_m on the last window, and the per-seed prediction Σ_k w^m_k‖S̃_k(w^m − w^M)‖² from the realised cluster weights against the measured D_m at s = 900:

| m | predicted E D_∞ | measured time average on the last window | ratio measured/predicted | per seed: linearised prediction Σ_k w^m_k‖S̃_k(w^m − w^M)‖² → exact stationary distance Σ_k w^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖² → measured D_m at the end [splitting margin of Θ*(w^m)] |
|---|---|---|---|---|
| 16 | 8.233e-01 | 4.459e-01 ± 2.8e-01 | 0.54 | 1.229e-01 → 1.118e-01 → 1.116e-01 [+0.00484]; 1.868e-01 → 1.859e-01 → 1.856e-01 [+0.00266]; 2.064e-01 → 2.047e-01 → 2.047e-01 [+0.00266]; 2.255e+00 → 1.292e+00 → 1.292e+00 [-0.00219] |
| 64 | 2.010e-01 | 1.377e-01 ± 5.1e-02 | 0.69 | 1.638e-02 → 2.174e-02 → 2.124e-02 [+0.00064]; 2.921e-01 → 2.414e-01 → 2.410e-01 [+0.00497]; 1.036e-01 → 8.984e-02 → 8.963e-02 [+0.00065]; 2.736e-01 → 2.013e-01 → 2.013e-01 [+0.00383] |
| 256 | 4.538e-02 | 1.358e-02 ± 5.5e-03 | 0.30 | 3.032e-02 → 2.623e-02 → 2.610e-02 [+0.00082]; 8.723e-04 → 8.541e-04 → 8.540e-04 [+0.00215]; 8.147e-03 → 8.786e-03 → 8.745e-03 [+0.00201]; 1.874e-02 → 1.905e-02 → 1.903e-02 [+0.00279] |

Per seed, the measured D_m at s = 900 is within 25% of the linearised prediction from the realised weights in 9 of the 12 (seed, m) pairs; the largest deviation is the ratio 0.57 at m = 16, seed 3, the pair with the largest weight fluctuation |w^m − w^M| = 0.275 (the linearisation's remainder and the saturation a ≤ a_max). The exact stationary distance Σ_k w^m_k‖ϑ*_k(w^m) − ϑ*_k(w^M)‖² (the first equality of Proposition D‴.5(ii), with the weighted stationary points themselves) reproduces the measured D_m within 2.3% in all 12 pairs: the m-systems sit at the replicated weighted stationary configurations. The splitting margin of Θ*(w^m) (σ₀ − max_k κ_R at its atoms with the weighted residual; Proposition D‴.5(i)) is positive in 11 of the 12 pairs: the configuration is a strict local minimiser of J_m there. It is negative for m = 16, seed 3 (w^m = (0.5625, 0.2500, 0.1875), margin -0.00219): there the replicated configuration is a saddle of J_m (a cluster can split), the m-system's plateau at s = 900 is the saddle's, and the system leaves it later by migrating particles between basins (referee round 9 continued the run to s = 3600: D grew from 1.29 to 5.13 as two particles changed basins and the weights re-sorted). With weight fluctuations of root-mean-square (Σ_k w̄_k(1 − w̄_k)/m)^{1/2} = 0.204 (m = 16), 0.102 (m = 64), 0.051 (m = 256) the admissible region of Proposition D‴.5(iv) is left with appreciable probability at m = 16. The multinomial expectation E D_∞ is compared with a mean over 4 seeds, whose standard error is that of a heavy-tailed quantity (E|w^m − w^M|² is dominated by the rare large fluctuations; for m = 256 the four realised fluctuations are simply small); the per-seed comparison is the sharp test.

Per-cluster check on seed 0 (whitened mean deviation of the m-system's particles of cluster k from their reference copies at s = 900, against S_k(w^m − w^M)):
- m = 16, weights w^m = (0.3750, 0.3125, 0.3125), w^M = (0.3184, 0.3354, 0.3462): cluster 0: |measured| 0.434, |predicted| 0.455, cosine 1.000; cluster 1: |measured| 0.130, |predicted| 0.129, cosine 0.992; cluster 2: |measured| 0.338, |predicted| 0.358, cosine 0.999.
- m = 64, weights w^m = (0.2969, 0.3438, 0.3594), w^M = (0.3184, 0.3354, 0.3462): cluster 0: |measured| 0.189, |predicted| 0.173, cosine 0.999; cluster 1: |measured| 0.050, |predicted| 0.048, cosine 0.923; cluster 2: |measured| 0.165, |predicted| 0.137, cosine 0.987.
- m = 256, weights w^m = (0.3086, 0.3125, 0.3789), w^M = (0.3184, 0.3354, 0.3462): cluster 0: |measured| 0.119, |predicted| 0.122, cosine 0.997; cluster 1: |measured| 0.136, |predicted| 0.131, cosine 0.995; cluster 2: |measured| 0.205, |predicted| 0.232, cosine 0.999.

## 4. Cold start (a connected spread cloud: a ~ N(0, 0.8²), c ~ U(−2, 2), s ~ U(0.4, 1); outside (H))

| m | max D_m on (0, 112.5] | max D_m on (112.5, 225] | max D_m on (225, 450] | max D_m on (450, 900] | time average on (450, 900] | bound C₁/m | mean ‖f̂_Θ − f̂_ref‖ on the last window | max_i ‖δ_i‖_{G⁻¹} at s = 900 (3ϱ′/4 = 0.075) | fraction of particles with ‖δ_i‖_{G⁻¹} > 1 at s = 900 |
|---|---|---|---|---|---|---|---|---|---|
| 16 | 1.809e+00 ± 2.8e-01 | 3.130e+00 ± 5.0e-01 | 4.498e+00 ± 2.8e-01 | 5.332e+00 ± 1.6e-01 | 5.027e+00 ± 1.3e-01 | 1.352e+08 | 2.527e-02 ± 4.1e-03 | 6.659e+00 | 0.203 ± 0.047 |
| 64 | 5.048e-01 ± 2.8e-01 | 9.004e-01 ± 5.3e-01 | 1.418e+00 ± 7.8e-01 | 1.785e+00 ± 9.3e-01 | 1.661e+00 ± 8.8e-01 | 3.379e+07 | 8.620e-03 ± 2.5e-03 | 6.272e+00 | 0.082 ± 0.021 |
| 256 | 7.437e-02 ± 1.6e-02 | 1.616e-01 ± 4.6e-02 | 3.473e-01 ± 1.2e-01 | 5.828e-01 ± 1.9e-01 | 4.860e-01 ± 1.7e-01 | 8.448e+06 | 4.249e-03 ± 9.2e-04 | 6.202e+00 | 0.037 ± 0.010 |

Slope of the last-window maximum of D_m against m (per-seed fits, mean ± s.e.): -0.87 ± 0.14 (Theorem D‴: −1 for the mean square; the function-level distance ‖f̂_Θ − f̂_ref‖ has slope -0.66 ± 0.06, m^{−1/2} being −0.5). Ratio of the last to the first window maxima: 2.95, 3.54, 7.84; last to the one before: 1.19, 1.26, 1.68 (uniform in time means no growth).
Slope against m of the fraction of pairs more than one whitened unit apart: -0.64 ± 0.08; of the median per-particle deviation: -0.62 ± 0.07 (m^{−1/2}: −0.5).

Per-particle deviations ‖δ_i‖_{G⁻¹} at s = 900 (pooled over the seeds): median, 90th percentile, maximum.

| start | m | median | 90% | max |
|---|---|---|---|---|
| warm | 16 | 4.456e-01 | 1.233e+00 | 1.361e+00 |
| warm | 64 | 3.573e-01 | 5.921e-01 | 6.197e-01 |
| warm | 256 | 1.140e-01 | 1.837e-01 | 2.144e-01 |
| cold | 16 | 4.005e-01 | 5.536e+00 | 6.659e+00 |
| cold | 64 | 1.445e-01 | 4.377e-01 | 6.272e+00 |
| cold | 256 | 6.742e-02 | 2.136e-01 | 6.202e+00 |

## 5. The real closed loop at k_e = 25 (σ = σ₀/k_e = 0.00028, τ = 0), warm start near the stationary measure

The m-systems share their particles with the M = 2048 reference; dt = 0.001; D_m in the G⁻¹ metric as above; real time t (slow time s = t/25).

| m | max D_m on (0, 500] | max D_m on (500, 1000] | max D_m on (1000, 2000] | time average on (1000, 2000] |
|---|---|---|---|---|
| 16 | 3.277e-01 ± 1.1e-01 | 4.671e-01 ± 1.4e-01 | 5.897e-01 ± 1.7e-01 | 5.367e-01 ± 1.5e-01 |
| 64 | 6.054e-02 ± 2.9e-02 | 8.008e-02 ± 3.9e-02 | 9.087e-02 ± 4.5e-02 | 8.655e-02 ± 4.2e-02 |
| 256 | 1.882e-02 ± 9.1e-03 | 2.531e-02 ± 1.1e-02 | 2.936e-02 ± 1.1e-02 | 2.775e-02 ± 1.1e-02 |

Slope of the last-window maximum against m: -1.13 ± 0.15. Ratio of the last to the first window maxima: 1.80, 1.50, 1.56. Reference fit error on the orbit: 0.2462 at t = 0, 0.2451 at t = 1000, 0.2461 at t = 2000.

The real loop against the averaged warm runs (section 3) at equal slow times s = t/25, windows (0, 20], (20, 40], (40, 80] (the averaged D_m is recorded every slow unit; the two sets of runs use different initial clouds; means ± s.e. over the seeds); the horizon 80 slow units is 0.16 relaxation times 1/κ = 512, the early transient of the averaged runs, whose window maxima settle only by s ≈ 450.

| m | averaged: max D_m on (0, 20] | averaged: max D_m on (20, 40] | averaged: max D_m on (40, 80] | real loop: max D_m on (0, 20] | real loop: max D_m on (20, 40] | real loop: max D_m on (40, 80] | real/averaged on the last window |
|---|---|---|---|---|---|---|---|
| 16 | 2.494e-01 ± 1.5e-01 | 3.249e-01 ± 1.9e-01 | 3.768e-01 ± 2.3e-01 | 3.277e-01 ± 1.1e-01 | 4.671e-01 ± 1.4e-01 | 5.897e-01 ± 1.7e-01 | 1.57 |
| 64 | 7.362e-02 ± 2.8e-02 | 1.013e-01 ± 3.9e-02 | 1.172e-01 ± 4.5e-02 | 6.054e-02 ± 2.9e-02 | 8.008e-02 ± 3.9e-02 | 9.087e-02 ± 4.5e-02 | 0.78 |
| 256 | 7.445e-03 ± 3.3e-03 | 9.843e-03 ± 4.2e-03 | 1.139e-02 ± 4.8e-03 | 1.882e-02 ± 9.1e-03 | 2.531e-02 ± 1.1e-02 | 2.936e-02 ± 1.1e-02 | 2.58 |

On the SAME clouds (the averaged flow started from the real loop's own initial clouds, reproduced from their seeds), the ratio D_real(t)/D_averaged(t/25) at the window ends s = 20, 40, 80: m = 16: 0.993–1.004; m = 64: 0.996–1.002; m = 256: 0.996–1.000 (over the 3 seeds): the real loop at k_e = 25 follows its averaged slow dynamics to within 0.7%. Splitting margins of the replicated configurations at the real-loop clouds' initial weights: m = 16, seed 0: -0.00207; m = 64, seed 0: +0.00326; m = 256, seed 0: +0.00450; m = 16, seed 1: +0.00257; m = 64, seed 1: +0.00494; m = 256, seed 1: +0.00413; m = 16, seed 2: -0.00001; m = 64, seed 2: +0.00288; m = 256, seed 2: +0.00304.

Reading: on the warm start near the stationary measure the hypothesis holds along the whole run (min a(s) = +0.0014); the per-particle localisation of part (c) is not in force at these m (max_i‖δ_i‖_{G⁻¹} = 1.36, 0.62, 0.21 against 3ϱ′/4 = 0.075; Remark D‴.6: the particles of a cluster are displaced together by the weight fluctuation); the window maxima of D_m do not grow over the horizon run — for the inadmissible weights (m = 16, seed 3) this is a saddle's plateau, left later by particle migration (Proposition D‴.5(iv)) — (last/first 1.13, 1.13, 1.15), the m-law of the mean square has slope -1.28 (theorem: −1), and the measured values lie below the bound C₁/m (which carries no transient factor here, Λ₁ = 0). On the cold start the support threads the ridges and (H) fails along the run (min a = -0.1961); the coupled distance grows over the horizon (last/first 2.95, 3.54, 7.84), its m-law has slope -0.87, and at the end a fraction 0.203, 0.082, 0.037 (m = 16, 64, 256) of the pairs sit more than one whitened unit apart: pairs in different basins, the floor of Remark D‴.2, whose fraction decreases with m (exponent in the cold-start block). The real loop on the warm start follows its averaged slow dynamics: on the same initial clouds D_real(t)/D_averaged(t/k_e) stays within 0.7% of one at every window end for every seed and m (across different clouds the last-window ratios are 1.57, 0.78, 2.58, within the seed-to-seed spread): its growth over t ≤ 2000 (last/first 1.80, 1.50, 1.56; 80 slow units, 0.16 relaxation times) is the early transient of the averaged dynamics, which saturates by s ≈ 450 (Corollary D‴.7); Theorem D″'s cold start at the same k_e grew threefold over 1000 s and the cold start here grows and carries the floor.

Numbers are produced by `experiments/exp_chaos_learning.py`; nothing is edited by hand.
