# Theorem D‴.21: a measure-level basin for the OT leak — the Polyak–Łojasiewicz inequality along the runs (`exp_ot_basin.py`)

The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit = 0.00707 (2σ₀ = 0.01415); sizes m = 15, 60, 240 and M = 510 from the warm clouds of Proposition D‴.12 (width 0.02 around the exact representation; the m-system the first m rows); 4 seeds; 600 slow units, midpoint rule dt = 0.02, leak target every 10 steps. E = 𝒥(μ_Θ) − 𝒥* with 𝒥 = L̄ + (σ₀/2)W²(·, ν) the OT functional, evaluated along both leaks' trajectories; D_a the class-respecting distance to μ*; the observed PL ratio |∇J_a|²_G/(mE), which along the OT flow is −d log E/ds; λ the theorem's rate from the plan-segment curvatures (9 points per segment, every 25 slow units). Entries: mean over the seeds (min–max).

## Prior: matching (centres shifted by 0)

μ*: Newton residual 2.3e-17, stationarity of the atoms 2.7e-17, 𝒥* = 5.0447e-33; diagonal plan: ε = 0.000 < Δ_ν/2 = 3.131. Identities (E), (P) along every run: largest relative residual 5.4e-15 (E) and 4.1e-15 (P).

| leak | m | E at s = 0, 25, 100, 300, 600 | D_a at s = 0, 100, 600 | W²(μ_Θ, μ*) at s = 0, 100, 600 | plan changes (particle relabellings) |
|---|---|---|---|---|---|
| ot | 15 | 3.88e-02, 1.51e-02, 8.00e-03, 6.29e-04, 8.08e-06 | 5.24e+00, 1.78e+00, 2.28e-03 | 5.24e+00, 1.78e+00, 2.28e-03 | 0–0 |
| ot | 60 | 1.44e-02, 7.11e-03, 2.84e-03, 1.59e-04, 2.25e-06 | 2.75e+00, 7.19e-01, 6.34e-04 | 2.75e+00, 7.19e-01, 6.34e-04 | 0–0 |
| ot | 240 | 5.34e-03, 3.17e-03, 1.41e-03, 7.59e-05, 1.06e-06 | 1.17e+00, 3.45e-01, 2.99e-04 | 1.17e+00, 3.45e-01, 2.99e-04 | 0–0 |
| ot | 510 | 4.26e-03, 2.60e-03, 1.05e-03, 5.85e-05, 8.29e-07 | 9.94e-01, 2.72e-01, 2.34e-04 | 9.94e-01, 2.72e-01, 2.34e-04 | 0–0 |
| nearest | 15 | 3.88e-02, 2.22e-02, 1.84e-02, 1.81e-02, 1.79e-02 | 5.24e+00, 4.80e+00, 4.71e+00 | 5.24e+00, 4.80e+00, 4.71e+00 | 0–2 |
| nearest | 60 | 1.44e-02, 1.06e-02, 9.39e-03, 9.20e-03, 9.04e-03 | 2.75e+00, 2.57e+00, 2.49e+00 | 2.75e+00, 2.57e+00, 2.49e+00 | 0–8 |
| nearest | 240 | 5.34e-03, 4.47e-03, 4.33e-03, 4.35e-03, 4.37e-03 | 1.17e+00, 1.19e+00, 1.21e+00 | 1.17e+00, 1.19e+00, 1.21e+00 | 0–10 |
| nearest | 510 | 4.26e-03, 3.69e-03, 3.59e-03, 3.61e-03, 3.62e-03 | 9.94e-01, 1.00e+00, 1.01e+00 | 9.94e-01, 1.00e+00, 1.01e+00 | 2–17 |

The PL inequality along the OT runs. s₊: the first segment record after which the margin α stays positive; the observed ratio and the bounds over the records with α > 0 (the chain ratio ≥ Q²/(DE) ≥ λ must hold); the late values: means over s ∈ [500, 600] (records with E > 10⁻²⁶).

| m | s₊ (seeds) | E(s₊)/E(0) | min over α > 0 of ratio/λ and Q²/(DE)/λ | late observed ratio | late λ | late α | late κ*, κ^Θ | a-posteriori bound E(0)e^{−∫λ₊} at s = 600 against E(600) |
|---|---|---|---|---|---|---|---|---|
| 15 | 150, 25, 25, 300 | 1.3e-01, 4.9e-01, 4.8e-01, 2.0e-02 | 1.00, 1.00 | 0.01417 | 0.01396 | 0.00703 | -3.5e-19, 9.4e-05 | 3.1e-04 ≥ 3.7e-06; 1.2e-04 ≥ 4.1e-06; 1.7e-05 ≥ 1.9e-06; 4.4e-03 ≥ 2.3e-05 |
| 60 | 0, 25, 125, 150 | 1.0e+00, 5.0e-01, 1.7e-01, 1.0e-01 | 1.00, 1.00 | 0.01417 | 0.01404 | 0.00705 | -3.5e-19, 5.2e-05 | 1.6e-06 ≥ 6.2e-07; 4.3e-05 ≥ 2.9e-06; 3.0e-05 ≥ 1.8e-06; 9.3e-05 ≥ 3.7e-06 |
| 240 | 150, 0, 0, 0 | 1.5e-01, 1.0e+00, 1.0e+00, 1.0e+00 | 1.00, 1.00 | 0.01417 | 0.01401 | 0.00704 | -3.5e-19, 7.0e-05 | 9.4e-05 ≥ 2.4e-06; 9.8e-06 ≥ 9.8e-07; 5.7e-06 ≥ 6.1e-07; 6.4e-07 ≥ 2.9e-07 |
| 510 | 125, 0, 0, 0 | 1.9e-01, 1.0e+00, 1.0e+00, 1.0e+00 | 1.00, 1.00 | 0.01416 | 0.01406 | 0.00705 | -3.5e-19, 4.6e-05 | 3.4e-05 ≥ 1.6e-06; 7.2e-06 ≥ 7.5e-07; 2.1e-07 ≥ 1.7e-07; 2.0e-06 ≥ 8.4e-07 |

The nearest-atom contrast: E at the end of the nearest-atom runs against the E of the replicated configuration of their sampled weights (a critical point of the nearest-atom energy; its drift under the nearest-atom leak and under the OT leak), and the admissibility of those weights: the smallest separation margin σ₀ − λ_max R^G_k(w) (the class), the smallest eigenvalue of the whitened weighted K-atom Hessian; inadmissible weights (a margin ≤ 0) have a saddle as their replicated configuration. E at s = 300, 450, 600; the end state's largest class spread (rms G⁻¹-distance to the class mean) and its largest distance from the atom of Θ*(w).

| m | seed | sampled counts | admissible | smallest separation margin (class) | K-atom Hessian min | E at s = 300, 450, 600 | E of the replicated configuration | E(600)/E_spur − 1 | class spread, largest distance | its drift: nearest / OT |
|---|---|---|---|---|---|---|---|---|---|---|
| 15 | 0 | (6, 5, 4) | no (saddle) | -0.00348 (0) | 0.00737 | 1.1428e-02, 1.1164e-02, 1.0760e-02 | 1.1802e-02 | -0.0883 | 0.173, 0.335 | 5e-17 / 2.3e-02 |
| 15 | 1 | (6, 3, 6) | no (saddle) | -0.00108 (0) | 0.00731 | 2.1047e-02, 2.1037e-02, 2.1025e-02 | 2.1148e-02 | -0.0058 | 0.067, 0.132 | 1e-17 / 4.4e-02 |
| 15 | 2 | (5, 4, 6) | yes | +0.00395 (2) | 0.00739 | 8.7195e-03, 8.7299e-03, 8.7354e-03 | 8.7421e-03 | -0.0008 | 0.003, 0.004 | 5e-17 / 3.6e-02 |
| 15 | 3 | (8, 4, 3) | no (saddle) | -0.00614 (0) | 0.00754 | 3.1122e-02, 3.0969e-02, 3.0918e-02 | 3.6728e-02 | -0.1582 | 0.619, 1.098 | 6e-14 / 4.4e-02 |
| 60 | 0 | (19, 21, 20) | yes | +0.00561 (1) | 0.00740 | 3.0928e-03, 3.0965e-03, 3.0980e-03 | 3.0989e-03 | -0.0003 | 0.001, 0.001 | 3e-15 / 4.4e-02 |
| 60 | 1 | (24, 15, 21) | no (saddle) | -0.00216 (0) | 0.00732 | 1.3577e-02, 1.3514e-02, 1.3452e-02 | 1.3772e-02 | -0.0233 | 0.161, 0.313 | 6e-17 / 4.4e-02 |
| 60 | 2 | (19, 18, 23) | yes | +0.00417 (2) | 0.00741 | 7.1898e-03, 7.2050e-03, 7.2128e-03 | 7.2213e-03 | -0.0012 | 0.003, 0.006 | 3e-17 / 3.6e-02 |
| 60 | 3 | (25, 17, 18) | no (saddle) | -0.00398 (0) | 0.00733 | 1.2942e-02, 1.2589e-02, 1.2395e-02 | 1.3925e-02 | -0.1098 | 0.354, 0.842 | 6e-17 / 4.4e-02 |
| 240 | 0 | (71, 76, 93) | yes | +0.00378 (2) | 0.00745 | 8.8488e-03, 8.8651e-03, 8.8741e-03 | 8.8855e-03 | -0.0013 | 0.003, 0.010 | 1e-17 / 3.6e-02 |
| 240 | 1 | (77, 76, 87) | yes | +0.00506 (2) | 0.00740 | 4.3000e-03, 4.3075e-03, 4.3108e-03 | 4.3136e-03 | -0.0006 | 0.001, 0.004 | 2e-14 / 3.6e-02 |
| 240 | 2 | (76, 81, 83) | yes | +0.00595 (2) | 0.00740 | 2.8479e-03, 2.8519e-03, 2.8535e-03 | 2.8546e-03 | -0.0004 | 0.001, 0.002 | 2e-15 / 4.4e-02 |
| 240 | 3 | (82, 79, 79) | yes | +0.00513 (0) | 0.00737 | 1.4181e-03, 1.4206e-03, 1.4217e-03 | 1.4226e-03 | -0.0006 | 0.001, 0.003 | 4e-17 / 4.4e-02 |
| 510 | 0 | (158, 161, 191) | yes | +0.00440 (2) | 0.00742 | 6.3525e-03, 6.3644e-03, 6.3703e-03 | 6.3765e-03 | -0.0010 | 0.002, 0.006 | 7e-17 / 3.6e-02 |
| 510 | 1 | (163, 166, 181) | yes | +0.00544 (2) | 0.00740 | 3.3250e-03, 3.3307e-03, 3.3332e-03 | 3.3350e-03 | -0.0006 | 0.001, 0.003 | 2e-15 / 3.6e-02 |
| 510 | 2 | (170, 167, 173) | yes | +0.00662 (2) | 0.00738 | 8.0901e-04, 8.1003e-04, 8.1041e-04 | 8.1062e-04 | -0.0003 | 0.001, 0.001 | 5e-17 / 3.6e-02 |
| 510 | 3 | (168, 184, 158) | yes | +0.00452 (1) | 0.00739 | 3.9423e-03, 3.9464e-03, 3.9484e-03 | 3.9502e-03 | -0.0005 | 0.001, 0.005 | 1e-14 / 4.4e-02 |

The a-priori basin of Theorem D‴.21(iv). Over the whole parameter space sup λ_max R^G[μ*] = 7.1420e-17 (approached at a = 4.000, c = 0.046, s = 0.300), so (H*) holds on the whole space with κ = 7.074e-03 and E₀ = κ²/(2(L₂^G)²) = 2.03e-07 (L₂^G = 11.111). Over the box of the end configurations, the priors and the atoms (not forward-invariant) sup = 5.8500e-17, κ = 0.00707, E₀ = 2.03e-07. A particle outside its cell is at G⁻¹-distance ≥ 3.1315 from its class atom, so {E ≤ E₀} holds a fraction p ≤ 5.8e-06 (whole space) or 5.8e-06 (box) of the particles in other cells: one such particle needs m ≥ 1.7e+05 or 1.7e+05. Local PL constant at μ*: 2min{λ_sep, λ_coll} = 0.01415 (λ_sep = 0.00707, λ_coll = 0.00738); the theorem's late λ (mean over the OT runs, s ≥ 500) is 0.991 of 2σ₀ = 0.991 (the slack of (ii)) × 1.000 (the curvature's cost, the local constant over 2σ₀).

## Prior: shifted (centres shifted by 0.15)

μ*: Newton residual 4.3e-17, stationarity of the atoms 5.5e-17, 𝒥* = 4.6447e-04; diagonal plan: ε = 0.389 < Δ_ν/2 = 3.162. Identities (E), (P) along every run: largest relative residual 5.6e-13 (E) and 4.2e-15 (P).

| leak | m | E at s = 0, 25, 100, 300, 600 | D_a at s = 0, 100, 600 | W²(μ_Θ, μ*) at s = 0, 100, 600 | plan changes (particle relabellings) |
|---|---|---|---|---|---|
| ot | 15 | 3.98e-02, 1.58e-02, 8.36e-03, 4.08e-04, 5.21e-06 | 5.25e+00, 1.71e+00, 1.51e-03 | 5.25e+00, 1.71e+00, 1.51e-03 | 0–0 |
| ot | 60 | 1.46e-02, 7.26e-03, 2.86e-03, 1.52e-04, 2.34e-06 | 2.76e+00, 7.04e-01, 6.84e-04 | 2.76e+00, 7.04e-01, 6.84e-04 | 0–0 |
| ot | 240 | 5.03e-03, 2.88e-03, 1.32e-03, 1.26e-04, 2.10e-06 | 1.16e+00, 3.76e-01, 5.96e-04 | 1.16e+00, 3.76e-01, 5.96e-04 | 0–0 |
| ot | 510 | 4.11e-03, 2.44e-03, 9.94e-04, 8.64e-05, 1.46e-06 | 9.93e-01, 2.89e-01, 4.17e-04 | 9.93e-01, 2.89e-01, 4.17e-04 | 0–0 |
| nearest | 15 | 3.98e-02, 2.31e-02, 1.93e-02, 1.89e-02, 1.86e-02 | 5.25e+00, 4.80e+00, 4.69e+00 | 5.25e+00, 4.80e+00, 4.69e+00 | 0–2 |
| nearest | 60 | 1.46e-02, 1.07e-02, 9.56e-03, 9.38e-03, 9.18e-03 | 2.76e+00, 2.57e+00, 2.49e+00 | 2.76e+00, 2.57e+00, 2.49e+00 | 0–8 |
| nearest | 240 | 5.03e-03, 4.18e-03, 4.04e-03, 4.06e-03, 4.07e-03 | 1.16e+00, 1.19e+00, 1.21e+00 | 1.16e+00, 1.19e+00, 1.21e+00 | 0–14 |
| nearest | 510 | 4.11e-03, 3.54e-03, 3.44e-03, 3.45e-03, 3.46e-03 | 9.93e-01, 1.00e+00, 1.01e+00 | 9.93e-01, 1.00e+00, 1.01e+00 | 2–23 |

The PL inequality along the OT runs. s₊: the first segment record after which the margin α stays positive; the observed ratio and the bounds over the records with α > 0 (the chain ratio ≥ Q²/(DE) ≥ λ must hold); the late values: means over s ∈ [500, 600] (records with E > 10⁻²⁶).

| m | s₊ (seeds) | E(s₊)/E(0) | min over α > 0 of ratio/λ and Q²/(DE)/λ | late observed ratio | late λ | late α | late κ*, κ^Θ | a-posteriori bound E(0)e^{−∫λ₊} at s = 600 against E(600) |
|---|---|---|---|---|---|---|---|---|
| 15 | 150, 125, 25, 275 | 9.5e-02, 1.5e-01, 4.6e-01, 1.6e-02 | 1.67, 1.64 | 0.01367 | 0.00766 | 0.00561 | 1.4e-03, 1.5e-03 | 1.7e-03 ≥ 2.3e-06; 1.6e-03 ≥ 4.7e-06; 3.8e-04 ≥ 3.4e-06; 1.3e-02 ≥ 1.0e-05 |
| 60 | 0, 25, 175, 125 | 1.0e+00, 5.1e-01, 1.0e-01, 1.3e-01 | 1.72, 1.69 | 0.01389 | 0.00778 | 0.00566 | 1.4e-03, 1.4e-03 | 5.0e-05 ≥ 8.5e-07; 6.6e-04 ≥ 2.6e-06; 7.6e-04 ≥ 3.6e-06; 9.2e-04 ≥ 2.3e-06 |
| 240 | 200, 25, 0, 0 | 1.0e-01, 5.8e-01, 1.0e+00, 1.0e+00 | 1.68, 1.64 | 0.01384 | 0.00778 | 0.00566 | 1.4e-03, 1.4e-03 | 1.2e-03 ≥ 4.9e-06; 3.4e-04 ≥ 2.1e-06; 2.0e-04 ≥ 1.2e-06; 3.1e-05 ≥ 2.2e-07 |
| 510 | 175, 0, 0, 0 | 1.2e-01, 1.0e+00, 1.0e+00, 1.0e+00 | 1.61, 1.58 | 0.01381 | 0.00779 | 0.00567 | 1.4e-03, 1.4e-03 | 6.7e-04 ≥ 3.3e-06; 2.5e-04 ≥ 1.6e-06; 8.4e-06 ≥ 3.3e-07; 6.8e-05 ≥ 6.4e-07 |

The nearest-atom contrast: E at the end of the nearest-atom runs against the E of the replicated configuration of their sampled weights (a critical point of the nearest-atom energy; its drift under the nearest-atom leak and under the OT leak), and the admissibility of those weights: the smallest separation margin σ₀ − λ_max R^G_k(w) (the class), the smallest eigenvalue of the whitened weighted K-atom Hessian; inadmissible weights (a margin ≤ 0) have a saddle as their replicated configuration. E at s = 300, 450, 600; the end state's largest class spread (rms G⁻¹-distance to the class mean) and its largest distance from the atom of Θ*(w).

| m | seed | sampled counts | admissible | smallest separation margin (class) | K-atom Hessian min | E at s = 300, 450, 600 | E of the replicated configuration | E(600)/E_spur − 1 | class spread, largest distance | its drift: nearest / OT |
|---|---|---|---|---|---|---|---|---|---|---|
| 15 | 0 | (6, 5, 4) | no (saddle) | -0.00287 (0) | 0.00737 | 1.2674e-02, 1.2476e-02, 1.2172e-02 | 1.3004e-02 | -0.0639 | 0.134, 0.263 | 4e-17 / 2.3e-02 |
| 15 | 1 | (6, 3, 6) | no (saddle) | -0.00038 (0) | 0.00733 | 2.1065e-02, 2.1062e-02, 2.1060e-02 | 2.1156e-02 | -0.0046 | 0.043, 0.083 | 2e-17 / 4.4e-02 |
| 15 | 2 | (5, 4, 6) | yes | +0.00278 (2) | 0.00743 | 8.2884e-03, 8.2983e-03, 8.3038e-03 | 8.3129e-03 | -0.0011 | 0.005, 0.009 | 6e-17 / 3.6e-02 |
| 15 | 3 | (8, 4, 3) | no (saddle) | -0.00576 (0) | 0.00749 | 3.3479e-02, 3.3162e-02, 3.3016e-02 | 3.9538e-02 | -0.1649 | 0.619, 1.190 | 5e-14 / 4.4e-02 |
| 60 | 0 | (19, 21, 20) | yes | +0.00555 (1) | 0.00742 | 2.9539e-03, 2.9572e-03, 2.9583e-03 | 2.9589e-03 | -0.0002 | 0.001, 0.002 | 2e-15 / 4.4e-02 |
| 60 | 1 | (24, 15, 21) | no (saddle) | -0.00152 (0) | 0.00732 | 1.3928e-02, 1.3878e-02, 1.3828e-02 | 1.4103e-02 | -0.0195 | 0.116, 0.217 | 3e-17 / 4.4e-02 |
| 60 | 2 | (19, 18, 23) | yes | +0.00305 (2) | 0.00744 | 6.7012e-03, 6.7170e-03, 6.7259e-03 | 6.7395e-03 | -0.0020 | 0.005, 0.012 | 3e-17 / 3.6e-02 |
| 60 | 3 | (25, 17, 18) | no (saddle) | -0.00341 (0) | 0.00733 | 1.3926e-02, 1.3514e-02, 1.3203e-02 | 1.4855e-02 | -0.1112 | 0.332, 0.897 | 4e-17 / 4.4e-02 |
| 240 | 0 | (71, 76, 93) | yes | +0.00271 (2) | 0.00748 | 8.1317e-03, 8.1484e-03, 8.1588e-03 | 8.1788e-03 | -0.0024 | 0.006, 0.018 | 4e-17 / 3.6e-02 |
| 240 | 1 | (77, 76, 87) | yes | +0.00390 (2) | 0.00742 | 3.9886e-03, 3.9967e-03, 4.0006e-03 | 4.0051e-03 | -0.0011 | 0.003, 0.009 | 2e-14 / 3.6e-02 |
| 240 | 2 | (76, 81, 83) | yes | +0.00474 (2) | 0.00742 | 2.6091e-03, 2.6138e-03, 2.6160e-03 | 2.6180e-03 | -0.0008 | 0.002, 0.004 | 2e-15 / 4.4e-02 |
| 240 | 3 | (82, 79, 79) | yes | +0.00574 (0) | 0.00739 | 1.5210e-03, 1.5234e-03, 1.5243e-03 | 1.5250e-03 | -0.0004 | 0.001, 0.002 | 7e-17 / 4.4e-02 |
| 510 | 0 | (158, 161, 191) | yes | +0.00329 (2) | 0.00745 | 5.8553e-03, 5.8678e-03, 5.8750e-03 | 5.8855e-03 | -0.0018 | 0.004, 0.011 | 1e-16 / 3.6e-02 |
| 510 | 1 | (163, 166, 181) | yes | +0.00426 (2) | 0.00742 | 3.0506e-03, 3.0569e-03, 3.0601e-03 | 3.0635e-03 | -0.0011 | 0.002, 0.007 | 3e-15 / 3.6e-02 |
| 510 | 2 | (170, 167, 173) | yes | +0.00529 (2) | 0.00740 | 7.7120e-04, 7.7217e-04, 7.7251e-04 | 7.7270e-04 | -0.0003 | 0.001, 0.003 | 1e-16 / 3.6e-02 |
| 510 | 3 | (168, 184, 158) | yes | +0.00449 (1) | 0.00741 | 4.1257e-03, 4.1302e-03, 4.1324e-03 | 4.1347e-03 | -0.0005 | 0.001, 0.005 | 1e-14 / 4.4e-02 |

The a-priori basin of Theorem D‴.21(iv). Over the whole parameter space sup λ_max R^G[μ*] = 6.6013e-03 (approached at a = 4.000, c = 1.824, s = 0.300), so (H*) holds on the whole space with κ = 4.731e-04 and E₀ = κ²/(2(L₂^G)²) = 9.06e-10 (L₂^G = 11.111). Over the box of the end configurations, the priors and the atoms (not forward-invariant) sup = 4.3902e-03, κ = 0.00268, E₀ = 2.92e-08. A particle outside its cell is at G⁻¹-distance ≥ 2.9471 from its class atom, so {E ≤ E₀} holds a fraction p ≤ 4.4e-07 (whole space) or 2.5e-06 (box) of the particles in other cells: one such particle needs m ≥ 2.3e+06 or 4.0e+05. Local PL constant at μ*: 2min{λ_sep, λ_coll} = 0.01135 (λ_sep = 0.00567, λ_coll = 0.00740); the theorem's late λ (mean over the OT runs, s ≥ 500) is 0.548 of 2σ₀ = 0.683 (the slack of (ii)) × 0.802 (the curvature's cost, the local constant over 2σ₀).

## The measure-level basin: particles moved into another cell (matching prior)

The lift of q with k particles of class 0 moved 0.6 of the way to atom 1 (into cell 1); the same fractions p = k/m at m = 30 and 300 (the same measures, so the same values). At s = 0: whether the plan keeps the moved particles in class 0, E, D_π, the current-residual curvature κ^Θ along the plan segments (129 points per segment), the margin α and the rate λ (κ* = 0 at the matching prior).

| m | k | p | kept in class 0 | E | D_π | κ^Θ | α | λ | observed ratio |
|---|---|---|---|---|---|---|---|---|---|
| 30 | 1 | 0.0333 | True | 3.346e-03 | 6.229e-01 | 9.572e-03 | +0.00229 | 0.00148 | 0.03141 |
| 30 | 2 | 0.0667 | True | 8.976e-03 | 1.246e+00 | 1.914e-02 | -0.00250 | 0.00000 | 0.03763 |
| 30 | 3 | 0.1000 | True | 1.689e-02 | 1.869e+00 | 2.872e-02 | -0.00728 | 0.00000 | 0.03952 |
| 30 | 5 | 0.1667 | True | 3.957e-02 | 3.114e+00 | 4.786e-02 | -0.01686 | 0.00000 | 0.03781 |
| 300 | 10 | 0.0333 | True | 3.346e-03 | 6.229e-01 | 9.572e-03 | +0.00229 | 0.00148 | 0.03141 |
| 300 | 20 | 0.0667 | True | 8.976e-03 | 1.246e+00 | 1.914e-02 | -0.00250 | 0.00000 | 0.03763 |
| 300 | 30 | 0.1000 | True | 1.689e-02 | 1.869e+00 | 2.872e-02 | -0.00728 | 0.00000 | 0.03952 |
| 300 | 50 | 0.1667 | True | 3.957e-02 | 3.114e+00 | 4.786e-02 | -0.01686 | 0.00000 | 0.03781 |

One particle moved (midpoint rule, ds = 0.5, plan every step): 

- m = 30, OT leak, 300 slow units: E 3.346e-03 → 2.546e-05 (E(0)exp(−∫λ₊) = 9.901e-05, λ evaluated at every step, smallest λ 0.00148); observed ratio at the end 0.01415; the moved particle labelled 0 by the OT plan, in cell 0 at the end.
- m = 30, nearest-atom leak, 3000 slow units (E is the OT functional's excess, not this leak's Lyapunov function): E 3.346e-03 → 5.7195e-03 (s = 300), 6.3042e-03 (s = 1000), 6.3126e-03 (s = 3000); the particle in cell 1; the counts (9, 11, 10), whose replicated configuration has E = 6.3126e-03 (admissible: smallest separation margin +0.00451, K-atom Hessian min 0.00743).
- m = 300, OT leak, 300 slow units: E 2.317e-04 → 3.093e-06 (E(0)exp(−∫λ₊) = 3.693e-06, λ evaluated at every step, smallest λ 0.01230); observed ratio at the end 0.01415; the moved particle labelled 0 by the OT plan, in cell 0 at the end.
- m = 300, nearest-atom leak, 3000 slow units (E is the OT functional's excess, not this leak's Lyapunov function): E 2.317e-04 → 5.5577e-04 (s = 300), 6.1284e-04 (s = 1000), 6.1328e-04 (s = 3000); the particle in cell 1; the counts (99, 101, 100), whose replicated configuration has E = 6.1328e-04 (admissible: smallest separation margin +0.00675, K-atom Hessian min 0.00738).

Numbers are produced by `experiments/exp_ot_basin.py`; nothing is edited by hand.
