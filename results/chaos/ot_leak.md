# Proposition D‴.12: the barycentric OT leak toward a prior measure — the mass constraint against the weight fluctuation (`exp_ot_leak.py`)

The three-bump instance of `exp_chaos_learning.py` (n = 1, a_max = 4.0, c_max = 3.0, s ∈ [0.3, 2.0], G = diag(1.0, 0.25, 0.25), (A, C, S) = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)], x_m = 2 sin t, ρ on 256 times) at σ₀ = 0.15σ_crit = 0.00707, with the matching prior ν = Σ_j q_jδ_{ϑ⁰_j} of `exp_prior_design.py` (three atoms at the exact representation (3A_k, C_k, S_k), masses q = (1/3, 1/3, 1/3), minimal inter-atom distance Δ_ν = 6.263 in the G⁻¹ metric). Two leaks on the same clouds: the nearest-atom leak (𝔄_m^ν) and the barycentric OT leak (𝔄_m^{OT}) with the plan recomputed every 10 steps of 0.02 slow units; warm start = the mixture of Gaussians of width 0.02 around the exact representation; 600 slow units (4.2/σ₀); 4 seeds; records every 1 slow units. The reference M stands in for the mean-field flow; the m-systems share their particles with it.

## 1. Commensurate sizes: m = 15, 60, 240 against M = 510 (mq ∈ ℕ³: every plan integral)

### Cluster masses: conserved by the nearest-atom leak, driven to q by the OT leak

| m | mean |w(0) − q|₁ over seeds (both leaks start from the same clouds) | ½|w^m(0) − w^M(0)|₁ (mean; the mass the two systems' clusters disagree on) | nearest-atom leak: max_s |w(s) − w(0)|_∞ | OT leak: |w(end) − q|_∞ (by nearest atom) | OT leak: equalisation time s_eq (first record from which w = q to within half a particle; per seed) | OT leak: end of the transit s_tr (first record from which no unsplit particle is assigned away from its nearest atom; per seed; recorded series — on the incommensurate sizes the recorded mask counted split rows, so only the commensurate column is exact) | OT leak: unsplit particles in transit at s = 0 → at the end (per seed; recomputed from the clouds) | split particles (OT, end) |
|---|---|---|---|---|---|---|---|---|
| 15 | 0.2333 | 0.1245 | 0.0e+00 | 0.0e+00 | 111, 94, 93, 164 | 111, 94, 93, 164 | 1, 2, 1, 3 → 0, 0, 0, 0 | 0, 0, 0, 0 |
| 60 | 0.1167 | 0.0632 | 0.0e+00 | 0.0e+00 | 96, 97, 117, 115 | 96, 97, 117, 115 | 1, 5, 3, 5 → 0, 0, 0, 0 | 0, 0, 0, 0 |
| 240 | 0.0542 | 0.0178 | 0.0e+00 | 0.0e+00 | 119, 108, 102, 98 | 119, 108, 102, 98 | 13, 7, 4, 2 → 0, 0, 0, 0 | 0, 0, 0, 0 |
| 510 (reference) | 0.0480 | 0.0000 | 0.0e+00 | 0.0e+00 | 113, 105, 97, 96 | 113, 105, 97, 96 | 21, 11, 3, 14 → 0, 0, 0, 0 | 0, 0, 0, 0 |

The pure leak carries a migrating particle across the Voronoi midpoint at s = ln 2/σ₀ = 98 slow units; the equalisation times above are where the last migrating particle's nearest atom flips.

### Fit errors ‖f − f̂‖/‖f‖ at the end of the run

| m | nearest-atom leak (mean over seeds; per seed) | OT leak (mean over seeds; per seed) | max whitened drift of the final configuration: nearest / OT (per seed) | migrated particles (nearest atom changed): count per seed; their distance to the new atom (G⁻¹; min–max over all) |
|---|---|---|---|---|
| 15 | 6.6735e-02 ± 2.0e-02; 0.0443, 0.0747, 0.0291, 0.1188 | 2.7757e-04 ± 1.1e-04; 3.0e-04, 1.5e-04, 8.9e-05, 5.8e-04 | 1.0e-03, 1.3e-04, 1.7e-05, 2.2e-04 / 7.3e-04, 5.3e-04, 5.0e-04, 8.1e-04 | 1, 2, 1, 3; 0.064–0.157 |
| 60 | 2.9377e-02 ± 7.3e-03; 0.0133, 0.0445, 0.0213, 0.0384 | 1.2374e-04 ± 3.5e-05; 3.6e-05, 1.1e-04, 1.4e-04, 2.1e-04 | 7.6e-06, 3.9e-04, 2.5e-05, 5.9e-04 / 6.2e-04, 5.7e-04, 7.0e-04, 7.4e-04 | 1, 5, 3, 5; 0.076–0.107 |
| 240 | 1.7155e-02 ± 6.5e-03; 0.0359, 0.0136, 0.0134, 0.0058 | 1.1824e-04 ± 4.7e-05; 2.5e-04, 9.9e-05, 1.0e-04, 2.4e-05 | 3.7e-05, 2.2e-05, 1.1e-05, 1.4e-05 / 6.8e-04, 6.8e-04, 6.6e-04, 6.7e-04 | 13, 7, 4, 2; 0.090–0.107 |
| 510 (reference) | 1.2205e-02 ± 4.1e-03; 0.0227, 0.0124, 0.0025, 0.0112 | 8.1276e-05 ± 3.6e-05; 1.7e-04, 1.0e-04, 9.4e-06, 4.1e-05 | 2.6e-05, 1.8e-05, 7.9e-06, 2.1e-05 / 6.8e-04, 6.7e-04, 5.2e-04, 6.2e-04 | 21, 11, 3, 14; 0.088–0.103 |

The pure leak's remaining approach of a migrated particle at s = 600: Δ_{kk'}e^{−σ₀·600} = 0.103, 0.090, 0.096 for the three atom pairs.

### The distances to the reference on the last window (450, 600] (means over the window and over seeds)

| m | ‖f̂_m − f̂_M‖²_{L²(ρ)}: nearest | OT | predicted (1/m − 1/M)tr(FC̃Fᵀ)/N (nearest, Proposition D‴.9(iii)) | W²_{G⁻¹}(μ_m, μ_M): nearest | OT | ½|w^m − w^M|₁·Δ_ν² (nearest: the cross-cluster transport of the weight mismatch, mean over seeds) | synchronous D: nearest | OT |
|---|---|---|---|---|---|---|---|---|
| 15 | 2.525e-03 ± 1.1e-03 | 4.927e-07 ± 4.4e-07 | 3.990e-03 | 5.146e+00 ± 1.2e+00 | 7.913e-03 ± 4.2e-03 | 4.884e+00 | 8.714e-01 ± 3.9e-01 | 5.765e+00 ± 1.3e+00 |
| 60 | 5.219e-04 ± 2.4e-04 | 4.890e-08 ± 1.1e-08 | 9.068e-04 | 2.677e+00 ± 5.8e-01 | 2.304e-03 ± 5.4e-04 | 2.480e+00 | 2.519e-01 ± 1.1e-01 | 3.316e+00 ± 7.5e-01 |
| 240 | 5.862e-05 ± 1.8e-05 | 7.466e-09 ± 4.0e-09 | 1.360e-04 | 8.071e-01 ± 2.1e-01 | 6.858e-04 ± 1.4e-04 | 6.970e-01 | 1.689e-02 ± 6.0e-03 | 9.431e-01 ± 2.7e-01 |

Nearest-atom leak, exact laws at the realised weights (last-window means over the predictions from the weighted stationary points Θ*(w^m), Θ*(w^M)): W² against the K × K transport between the two weighted stationary measures, measured/exact 0.997–0.999 in the 7 pairs whose realised weights are admissible and 0.846, 0.933, 0.950, 0.974, 0.993 at the 5 saddle weights (where the m-system is not at the replicated configuration); the formula ½|w^m − w^M|₁Δ_ν² is within -7% to +28% of the measured W² (measured/formula 0.928–1.275). Function level against ‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖²_{L²(ρ)}: measured/exact 1.000–1.001 (admissible pairs), 0.845, 0.969, 0.974, 0.999, 1.005 at the saddles; the w̄-law of the table is the leading term at the mean weights.

OT leak, the synchronous D: its plateau is the mass reassigned in one system and not in the other, Σ_{i: label_m(i) ≠ label_M(i)}Δ²_{kk'}/m from the end labels; measured last-window D/plateau = 0.947–0.964 over the 12 pairs with mismatched particles — D converges to the plateau from below (the migrated particles still approaching their atoms) and does not vanish.

Ratios OT/nearest of the function-level distance: 2.0e-04, 9.4e-05, 1.3e-04; of W²: 1.5e-03, 8.6e-04, 8.5e-04. Ratios of the last window to the one before ((300, 450]; 1 = stationary): function level nearest 1.00, 0.97, 1.00, OT 0.018, 0.035, 0.028; W² nearest 0.99, 1.00, 1.00, OT 0.118, 0.120, 0.121; synchronous D nearest 1.00, 1.00, 1.00, OT 1.079, 1.092, 1.104. Slopes in m (log–log, per seed, mean ± s.e.): function level nearest -1.36 ± 0.37, OT -1.31 ± 0.60; W² nearest -0.67 ± 0.13, OT -0.77 ± 0.17; synchronous D nearest -1.40 ± 0.26, OT -0.71 ± 0.20.

### Time course under the OT leak: window means of the function-level distance² and of W² (the transit, then the stationary state)

| m | ‖f̂_m − f̂_M‖² on (0, 75] | ‖f̂_m − f̂_M‖² on (75, 150] | ‖f̂_m − f̂_M‖² on (150, 300] | ‖f̂_m − f̂_M‖² on (300, 450] | ‖f̂_m − f̂_M‖² on (450, 600] | W² on (0, 75] | W² on (75, 150] | W² on (150, 300] | W² on (300, 450] | W² on (450, 600] |
|---|---|---|---|---|---|---|---|---|---|---|
| 15 | 9.699e-03 | 2.630e-03 | 9.228e-04 | 2.746e-05 | 4.927e-07 | 3.370e+00 | 1.721e+00 | 5.591e-01 | 6.709e-02 | 7.913e-03 |
| 60 | 2.337e-03 | 4.059e-04 | 5.972e-05 | 1.393e-06 | 4.890e-08 | 1.748e+00 | 6.859e-01 | 1.633e-01 | 1.925e-02 | 2.304e-03 |
| 240 | 2.113e-04 | 6.387e-05 | 9.737e-06 | 2.621e-07 | 7.466e-09 | 5.144e-01 | 2.029e-01 | 4.743e-02 | 5.654e-03 | 6.858e-04 |

## 2. Incommensurate sizes: m = 16, 64, 256 against M = 512 (mq ∉ ℕ³: one or two split particles)

### Cluster masses: conserved by the nearest-atom leak, driven to q by the OT leak

| m | mean |w(0) − q|₁ over seeds (both leaks start from the same clouds) | ½|w^m(0) − w^M(0)|₁ (mean; the mass the two systems' clusters disagree on) | nearest-atom leak: max_s |w(s) − w(0)|_∞ | OT leak: |w(end) − q|_∞ (by nearest atom) | OT leak: equalisation time s_eq (first record from which w = q to within half a particle; per seed) | OT leak: end of the transit s_tr (first record from which no unsplit particle is assigned away from its nearest atom; per seed; recorded series — on the incommensurate sizes the recorded mask counted split rows, so only the commensurate column is exact) | OT leak: unsplit particles in transit at s = 0 → at the end (per seed; recomputed from the clouds) | split particles (OT, end) |
|---|---|---|---|---|---|---|---|---|
| 16 | 0.2188 | 0.1133 | 0.0e+00 | 4.2e-02 | —, —, —, — | n/a (split rows) | 0, 0, 0, 2 → 0, 0, 0, 0 | 1, 2, 2, 1 |
| 64 | 0.1276 | 0.0645 | 0.0e+00 | 1.0e-02 | —, —, —, — | n/a (split rows) | 1, 4, 3, 4 → 0, 0, 0, 0 | 2, 2, 2, 2 |
| 256 | 0.0469 | 0.0161 | 0.0e+00 | 2.6e-03 | —, —, —, — | n/a (split rows) | 11, 4, 3, 2 → 0, 0, 0, 0 | 2, 2, 2, 2 |
| 512 (reference) | 0.0452 | 0.0000 | 0.0e+00 | 1.3e-03 | —, —, —, — | n/a (split rows) | 19, 9, 1, 12 → 0, 0, 0, 0 | 2, 2, 2, 2 |

The pure leak carries a migrating particle across the Voronoi midpoint at s = ln 2/σ₀ = 98 slow units; the equalisation times above are where the last migrating particle's nearest atom flips.

### Fit errors ‖f − f̂‖/‖f‖ at the end of the run

| m | nearest-atom leak (mean over seeds; per seed) | OT leak (mean over seeds; per seed) | max whitened drift of the final configuration: nearest / OT (per seed) | migrated particles (nearest atom changed): count per seed; their distance to the new atom (G⁻¹; min–max over all) |
|---|---|---|---|---|
| 16 | 6.1204e-02 ± 2.6e-02; 0.0245, 0.0408, 0.0408, 0.1387 | 2.8235e-02 ± 1.5e-03; 2.6e-02, 3.1e-02, 3.1e-02, 2.6e-02 | 8.6e-05, 2.4e-05, 1.6e-05, 1.6e-04 / 1.4e-03, 1.4e-03, 1.4e-03, 1.5e-03 | 1, 2, 2, 4; 0.083–3.651 |
| 64 | 3.4548e-02 ± 3.7e-03; 0.0323, 0.0390, 0.0251, 0.0418 | 9.3271e-03 ± 2.0e-04; 9.9e-03, 9.3e-03, 9.0e-03, 9.2e-03 | 1.3e-05, 5.6e-04, 2.0e-05, 5.2e-04 / 1.3e-03, 5.6e-04, 5.0e-04, 1.7e-03 | 3, 4, 3, 6; 0.046–2.235 |
| 256 | 1.3025e-02 ± 4.1e-03; 0.0243, 0.0091, 0.0135, 0.0052 | 2.3299e-03 ± 2.3e-04; 2.3e-03, 2.4e-03, 2.9e-03, 1.7e-03 | 2.7e-05, 9.7e-06, 7.9e-06, 8.3e-06 / 5.9e-04, 5.7e-04, 7.2e-04, 5.1e-04 | 11, 4, 5, 4; 0.041–2.464 |
| 512 (reference) | 1.1570e-02 ± 4.1e-03; 0.0219, 0.0117, 0.0020, 0.0107 | 1.0517e-03 ± 1.3e-04; 1.1e-03, 1.2e-03, 1.2e-03, 6.7e-04 | 2.5e-05, 1.4e-05, 6.6e-06, 2.0e-05 / 8.3e-04, 6.8e-04, 5.2e-04, 6.2e-04 | 21, 11, 3, 14; 0.071–2.456 |

The pure leak's remaining approach of a migrated particle at s = 600: Δ_{kk'}e^{−σ₀·600} = 0.103, 0.090, 0.096 for the three atom pairs.

### The distances to the reference on the last window (450, 600] (means over the window and over seeds)

| m | ‖f̂_m − f̂_M‖²_{L²(ρ)}: nearest | OT | predicted (1/m − 1/M)tr(FC̃Fᵀ)/N (nearest, Proposition D‴.9(iii)) | W²_{G⁻¹}(μ_m, μ_M): nearest | OT | ½|w^m − w^M|₁·Δ_ν² (nearest: the cross-cluster transport of the weight mismatch, mean over seeds) | synchronous D: nearest | OT |
|---|---|---|---|---|---|---|---|---|
| 16 | 2.476e-03 ± 1.7e-03 | 3.406e-04 ± 4.5e-05 | 3.733e-03 | 4.782e+00 ± 1.4e+00 | 1.157e+00 ± 1.1e-01 | 4.443e+00 | 8.146e-01 ± 5.4e-01 | 4.651e+00 ± 1.8e+00 |
| 64 | 4.970e-04 ± 2.0e-04 | 3.117e-05 ± 2.3e-06 | 8.430e-04 | 2.738e+00 ± 6.2e-01 | 2.854e-01 ± 1.3e-02 | 2.528e+00 | 2.642e-01 ± 1.2e-01 | 2.799e+00 ± 7.5e-01 |
| 256 | 4.426e-05 ± 2.5e-05 | 1.115e-06 ± 1.5e-07 | 1.204e-04 | 6.940e-01 ± 3.2e-01 | 4.114e-02 ± 5.0e-03 | 6.320e-01 | 1.776e-02 ± 1.2e-02 | 6.962e-01 ± 3.0e-01 |

Nearest-atom leak, exact laws at the realised weights (last-window means over the predictions from the weighted stationary points Θ*(w^m), Θ*(w^M)): W² against the K × K transport between the two weighted stationary measures, measured/exact 0.997–0.999 in the 8 pairs whose realised weights are admissible and 0.859, 0.911, 0.938, 0.960 at the 4 saddle weights (where the m-system is not at the replicated configuration); the formula ½|w^m − w^M|₁Δ_ν² is within -0% to +28% of the measured W² (measured/formula 0.996–1.279). Function level against ‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖²_{L²(ρ)}: measured/exact 1.000–1.001 (admissible pairs), 0.818, 0.906, 0.990, 1.009 at the saddles; the w̄-law of the table is the leading term at the mean weights.

OT leak, the synchronous D: its plateau is the mass reassigned in one system and not in the other, Σ_{i: label_m(i) ≠ label_M(i)}Δ²_{kk'}/m from the end labels; measured last-window D/plateau = 0.449–1.116 over the 12 pairs with mismatched particles — D does not vanish; on these sizes the split (bridge) particles of the two systems add to and subtract from the plateau, so the ratio is not one-sided.

Ratios OT/nearest of the function-level distance: 1.4e-01, 6.3e-02, 2.5e-02; of W²: 2.4e-01, 1.0e-01, 5.9e-02. Ratios of the last window to the one before ((300, 450]; 1 = stationary): function level nearest 1.00, 0.97, 1.00, OT 0.819, 0.898, 0.971; W² nearest 1.00, 0.99, 1.00, OT 0.973, 1.002, 0.977; synchronous D nearest 1.00, 1.00, 1.00, OT 1.110, 1.090, 1.101. Slopes in m (log–log, per seed, mean ± s.e.): function level nearest -1.59 ± 0.27, OT -2.07 ± 0.05; W² nearest -0.79 ± 0.12, OT -1.21 ± 0.04; synchronous D nearest -1.60 ± 0.26, OT -0.72 ± 0.19.

### Time course under the OT leak: window means of the function-level distance² and of W² (the transit, then the stationary state)

| m | ‖f̂_m − f̂_M‖² on (0, 75] | ‖f̂_m − f̂_M‖² on (75, 150] | ‖f̂_m − f̂_M‖² on (150, 300] | ‖f̂_m − f̂_M‖² on (300, 450] | ‖f̂_m − f̂_M‖² on (450, 600] | W² on (0, 75] | W² on (75, 150] | W² on (150, 300] | W² on (300, 450] | W² on (450, 600] |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 9.600e-03 | 2.938e-03 | 1.380e-03 | 4.159e-04 | 3.406e-04 | 3.420e+00 | 2.370e+00 | 1.582e+00 | 1.190e+00 | 1.157e+00 |
| 64 | 2.424e-03 | 3.574e-04 | 7.067e-05 | 3.471e-05 | 3.117e-05 | 1.817e+00 | 8.436e-01 | 3.848e-01 | 2.848e-01 | 2.854e-01 |
| 256 | 1.868e-04 | 2.847e-05 | 4.856e-06 | 1.148e-06 | 1.115e-06 | 4.440e-01 | 1.832e-01 | 6.858e-02 | 4.209e-02 | 4.114e-02 |

### The quantisation laws of (iv) against the mean-field stationary measure μ*(q) = ν (OT leak; the reference carries its own bridge particles)

| m | W²(μ_m, ν) on the last window (mean over seeds) | m·W²(μ_m, ν) | fit error ‖f̂_m − f‖/‖f‖ at the end (mean) | m·fit |
|---|---|---|---|---|
| 16 | 1.213e+00 ± 9.7e-02 | 19.4 | 2.8235e-02 ± 1.5e-03 | 0.45 |
| 64 | 3.280e-01 ± 1.3e-02 | 21.0 | 9.3271e-03 ± 2.0e-04 | 0.60 |
| 256 | 7.853e-02 ± 2.6e-03 | 20.1 | 2.3299e-03 ± 2.3e-04 | 0.60 |
| 512 (reference) | 3.808e-02 ± 8.0e-04 | 19.5 | 1.0517e-03 ± 1.3e-04 | 0.54 |

Slopes in m: W²(μ_m, ν) -0.98 ± 0.03 over m = 16, 64, 256 (-1.00 ± 0.03 with M = 512); fit error squared -1.81 ± 0.05 (-1.91 ± 0.07 with M): the order 1/m of (iv) in W² and 1/m² at function level, with the constants m·W² ≈ 20 = 0.51Δ_ν² and m·fit ≈ 0.55.

### The split particles at the end of the OT runs (plan rows m·γ_i; distance to the nearest prior atom in units of Δ_ν; distance to the barycentric target T̄_i)

- seed 0, m = 16: particle 6, weights (0.333, 0.333, 0.333), dist/Δ_ν = 0.580, ‖ϑ − T̄‖ = 0.578
- seed 0, m = 64: particle 17, weights (0.667, 0.333, 0.000), dist/Δ_ν = 0.357, ‖ϑ − T̄‖ = 0.293
- seed 0, m = 64: particle 50, weights (0.667, 0.000, 0.333), dist/Δ_ν = 0.320, ‖ϑ − T̄‖ = 0.258
- seed 0, m = 256: particle 180, weights (0.333, 0.000, 0.667), dist/Δ_ν = 0.338, ‖ϑ − T̄‖ = 0.114
- seed 0, m = 256: particle 182, weights (0.000, 0.333, 0.667), dist/Δ_ν = 0.323, ‖ϑ − T̄‖ = 0.073
- seed 0, m = 512: particle 375, weights (0.000, 0.667, 0.333), dist/Δ_ν = 0.342, ‖ϑ − T̄‖ = 0.056
- seed 0, m = 512: particle 476, weights (0.667, 0.000, 0.333), dist/Δ_ν = 0.367, ‖ϑ − T̄‖ = 0.072
- seed 1, m = 16: particle 6, weights (0.000, 0.667, 0.333), dist/Δ_ν = 0.277, ‖ϑ − T̄‖ = 0.413
- seed 1, m = 16: particle 12, weights (0.333, 0.667, 0.000), dist/Δ_ν = 0.386, ‖ϑ − T̄‖ = 0.025
- seed 1, m = 64: particle 39, weights (0.667, 0.000, 0.333), dist/Δ_ν = 0.276, ‖ϑ − T̄‖ = 0.509
- seed 1, m = 64: particle 55, weights (0.667, 0.333, 0.000), dist/Δ_ν = 0.344, ‖ϑ − T̄‖ = 0.336
- seed 1, m = 256: particle 137, weights (0.333, 0.000, 0.667), dist/Δ_ν = 0.337, ‖ϑ − T̄‖ = 0.116
- seed 1, m = 256: particle 199, weights (0.000, 0.333, 0.667), dist/Δ_ν = 0.323, ‖ϑ − T̄‖ = 0.073
- seed 1, m = 512: particle 3, weights (0.667, 0.000, 0.333), dist/Δ_ν = 0.364, ‖ϑ − T̄‖ = 0.051
- seed 1, m = 512: particle 361, weights (0.000, 0.667, 0.333), dist/Δ_ν = 0.342, ‖ϑ − T̄‖ = 0.052
- seed 2, m = 16: particle 3, weights (0.333, 0.667, 0.000), dist/Δ_ν = 0.386, ‖ϑ − T̄‖ = 0.026
- seed 2, m = 16: particle 11, weights (0.000, 0.667, 0.333), dist/Δ_ν = 0.277, ‖ϑ − T̄‖ = 0.414
- seed 2, m = 64: particle 27, weights (0.000, 0.333, 0.667), dist/Δ_ν = 0.309, ‖ϑ − T̄‖ = 0.198
- seed 2, m = 64: particle 34, weights (0.333, 0.000, 0.667), dist/Δ_ν = 0.299, ‖ϑ − T̄‖ = 0.357
- seed 2, m = 256: particle 67, weights (0.667, 0.000, 0.333), dist/Δ_ν = 0.354, ‖ϑ − T̄‖ = 0.033
- seed 2, m = 256: particle 125, weights (0.667, 0.333, 0.000), dist/Δ_ν = 0.385, ‖ϑ − T̄‖ = 0.070
- seed 2, m = 512: particle 11, weights (0.000, 0.667, 0.333), dist/Δ_ν = 0.341, ‖ϑ − T̄‖ = 0.049
- seed 2, m = 512: particle 337, weights (0.667, 0.000, 0.333), dist/Δ_ν = 0.360, ‖ϑ − T̄‖ = 0.029
- seed 3, m = 16: particle 4, weights (0.333, 0.333, 0.333), dist/Δ_ν = 0.583, ‖ϑ − T̄‖ = 0.595
- seed 3, m = 64: particle 54, weights (0.000, 0.333, 0.667), dist/Δ_ν = 0.325, ‖ϑ − T̄‖ = 0.073
- seed 3, m = 64: particle 58, weights (0.333, 0.000, 0.667), dist/Δ_ν = 0.352, ‖ϑ − T̄‖ = 0.087
- seed 3, m = 256: particle 16, weights (0.000, 0.667, 0.333), dist/Δ_ν = 0.341, ‖ϑ − T̄‖ = 0.046
- seed 3, m = 256: particle 192, weights (0.333, 0.667, 0.000), dist/Δ_ν = 0.393, ‖ϑ − T̄‖ = 0.063
- seed 3, m = 512: particle 278, weights (0.000, 0.333, 0.667), dist/Δ_ν = 0.342, ‖ϑ − T̄‖ = 0.057
- seed 3, m = 512: particle 382, weights (0.667, 0.333, 0.000), dist/Δ_ν = 0.392, ‖ϑ − T̄‖ = 0.057

Forest bound on the split mass: (K₀ − 1)/m = 0.1250, 0.0312, 0.0078 for m = 16, 64, 256; its worst-case W² contribution (K₀ − 1)Δ_ν²/m = 4.90, 1.23, 0.31.

## 3. The masses as a design input: the weighted stationary point Θ*(q_δ) for q_δ = (1/3 + δ, 1/3 − δ/2, 1/3 − δ/2) (no dynamics)

| δ | q_δ | fit error of Θ*(q_δ) (the bias of the mass constraint) | atoms (a, c, s) | κ_R at the atoms | κ_R over N_ϱ′ | margin κ = σ₀ − κ_R(N_ϱ′) | splitting margin σ₀ − max κ_R(atoms) (admissibility of q_δ) | ε = max‖ϑ*_k − ϑ⁰_k‖ against Δ_ν/2 (diagonal plan) |
|---|---|---|---|---|---|---|---|---|
| 0 | (0.3333, 0.3333, 0.3333) | 0.0000 | (3.000, -1.500, 0.500), (-2.100, 0.000, 0.400), (2.400, 1.500, 0.500) | +0.00000, +0.00000, -0.00000 | +0.00000 | +0.00707 | +0.00707 | 0.000 < 3.131: yes |
| 0.01 | (0.3433, 0.3283, 0.3283) | 0.0069 (0.69δ) | (2.945, -1.498, 0.497), (-2.121, -0.000, 0.400), (2.428, 1.501, 0.500) | +0.00230, +0.00005, -0.00002 | +0.00241 | +0.00466 | +0.00477 | 0.128 < 3.131: yes |
| 0.02 | (0.3533, 0.3233, 0.3233) | 0.0131 (0.66δ) | (2.890, -1.497, 0.494), (-2.143, -0.001, 0.401), (2.457, 1.502, 0.501) | +0.00428, +0.00009, -0.00005 | +0.00448 | +0.00259 | +0.00279 | 0.250 < 3.131: yes |
| 0.03 | (0.3633, 0.3183, 0.3183) | 0.0188 (0.63δ) | (2.835, -1.496, 0.491), (-2.166, -0.001, 0.401), (2.486, 1.502, 0.502) | +0.00596, +0.00013, -0.00008 | +0.00624 | +0.00084 | +0.00111 | 0.367 < 3.131: yes |
| 0.05 | (0.3833, 0.3083, 0.3083) | 0.0287 (0.57δ) | (2.725, -1.494, 0.488), (-2.212, -0.002, 0.401), (2.545, 1.504, 0.503) | +0.00856, +0.00019, -0.00014 | +0.00896 | -0.00188 | -0.00149 | 0.584 < 3.131: yes |
| 0.1 | (0.4333, 0.2833, 0.2833) | 0.0491 (0.49δ) | (2.467, -1.491, 0.482), (-2.334, -0.004, 0.403), (2.698, 1.510, 0.507) | +0.01193, +0.00029, -0.00036 | +0.01248 | -0.00541 | -0.00485 | 1.036 < 3.131: yes |

Reading: the stationary measure of the OT leak carries the prior's masses, so a mass mismatch is a bias (fit error 1.1e-16 at δ = 0, 0.049 at δ = 0.1) and, from a threshold between δ = 0.03 and 0.05, an instability: the over-weighted first cluster overshoots the plant, its splitting margin turns negative and the equal-mass clusters split (Remark D‴.11); the designed masses must be admissible — the same asymmetry as the matching prior's admissible set in `prior_design.md`. The diagonal plan is the unique optimal plan at every row (ε < Δ_ν/2).

Numbers are produced by `experiments/exp_ot_leak.py`; nothing is edited by hand.
