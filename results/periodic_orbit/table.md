# Theorem J: periodic orbits of the real closed loop near the critical points of J_m (`exp_periodic_orbit.py`)

The period mean Θ̄_p of the exact periodic orbit (Newton shooting on the Poincaré map, residual < 1e-10) is compared with the critical point Θ* of J_m and with the first-order prediction Θ* − J⁻¹(b_frozen + b_osc) of Proposition H.7 and Remark H.7a. The parameter Floquet exponents of the orbit (from the monodromy matrix) are compared with the eigenvalues of the averaged linearisation. Distances and cosines are in the labelled G⁻¹ metric. The error-channel multipliers are below the resolution of the finite-difference monodromy (exponents of order k_e); their product is recovered from Liouville's identity Σ log|λ| = ∫tr ∂F/∂z dt along the orbit, and the mean error-channel exponent is reported against k_e (Theorem J(b): k_e − L with L = O(1)).

### Instance A (n = 2, m = 3, collapsed minimiser)

σ₀ = 0.5, G = diag(4.0, 4.0, 1.0, 1.0, 1.0); Θ* has J_m = 0.762775, λ_min(∇²J_m) = 0.0148, kernel spread 4.14e-14 (0 means a collapsed configuration); its kernels have amplitudes (-0.036, 0.871), (-0.036, 0.871), (-0.036, 0.871), centres (-0.006, -0.070), (-0.006, -0.070), (-0.006, -0.070) and widths 0.919, 0.919, 0.919.

| k_e | shooting residual | max‖e‖ on the orbit | ‖Θ̄_p − Θ*‖ | full prediction: size ratio, cosine, vector error | frozen only: size ratio, cosine, vector error | parameter Floquet exponents against (1/k_e)λ(G_m^{1/2}∇²J_mG_m^{1/2}): max rel. error, median rel. error | k_e × max rel. error | mean error-channel exponent from Liouville's identity, minus k_e |
|---|---|---|---|---|---|---|---|---|
| 10 | 1.0e-13 | 8.77e-02 | 4.088e-02 | 1.149, 0.9972, 0.169 | 1.282, 0.9693, 0.398 | 8.310e-02, 1.208e-02 | 0.831 | 9.840, -0.160 |
| 20 | 3.6e-15 | 4.31e-02 | 1.910e-02 | 1.074, 0.9993, 0.083 | 1.198, 0.9624, 0.360 | 7.335e-02, 7.286e-03 | 1.467 | 19.889, -0.111 |
| 40 | 1.2e-15 | 2.13e-02 | 9.218e-03 | 1.036, 0.9998, 0.041 | 1.157, 0.9584, 0.348 | 4.498e-02, 3.961e-03 | 1.799 | 39.914, -0.086 |
| 80 | 1.1e-15 | 1.06e-02 | 4.527e-03 | 1.018, 1.0000, 0.020 | 1.136, 0.9563, 0.343 | 2.461e-02, 2.064e-03 | 1.969 | 79.926, -0.074 |
| 160 | 8.6e-12 | 5.26e-03 | 2.244e-03 | 1.009, 1.0000, 0.010 | 1.126, 0.9552, 0.342 | 1.284e-02, 1.053e-03 | 2.054 | 159.932, -0.068 |

Fitted exponent of ‖Θ̄_p − Θ*‖ against k_e: -1.05 (Theorem J: −1). Vector error of the full prediction: 0.169, 0.083, 0.041, 0.020, 0.010 (Theorem J(c): O(1/k_e)).

### Instance B (n = 1, m = 2, two distinct kernels)

σ₀ = 0.05, G = diag(2.0, 1.0, 1.0); Θ* has J_m = 0.258249, λ_min(∇²J_m) = 0.0661, kernel spread 2.10e+00 (0 means a collapsed configuration); its kernels have amplitudes (1.947), (-1.667), centres (0.929), (-0.906) and widths 0.782, 0.765.

| k_e | shooting residual | max‖e‖ on the orbit | ‖Θ̄_p − Θ*‖ | full prediction: size ratio, cosine, vector error | frozen only: size ratio, cosine, vector error | parameter Floquet exponents against (1/k_e)λ(G_m^{1/2}∇²J_mG_m^{1/2}): max rel. error, median rel. error | k_e × max rel. error | mean error-channel exponent from Liouville's identity, minus k_e |
|---|---|---|---|---|---|---|---|---|
| 10 | 1.5e-14 | 7.26e-02 | 6.445e-02 | 2.818, 0.8797, 1.995 | 2.818, 0.8797, 1.995 | 9.475e-02, 4.013e-02 | 0.948 | 10.054, +0.054 |
| 20 | 3.9e-15 | 3.77e-02 | 2.221e-02 | 1.942, 0.9292, 1.078 | 1.942, 0.9292, 1.078 | 3.210e-02, 1.387e-02 | 0.642 | 20.074, +0.074 |
| 40 | 4.4e-15 | 1.91e-02 | 8.370e-03 | 1.464, 0.9696, 0.551 | 1.464, 0.9696, 0.551 | 1.156e-02, 5.440e-03 | 0.462 | 40.084, +0.084 |
| 80 | 3.4e-13 | 9.60e-03 | 3.507e-03 | 1.227, 0.9896, 0.277 | 1.227, 0.9896, 0.277 | 4.628e-03, 2.228e-03 | 0.370 | 80.089, +0.089 |
| 160 | 2.4e-14 | 4.81e-03 | 1.589e-03 | 1.111, 0.9969, 0.139 | 1.111, 0.9969, 0.139 | 2.026e-03, 9.503e-04 | 0.324 | 160.092, +0.092 |

Fitted exponent of ‖Θ̄_p − Θ*‖ against k_e: -1.33 (Theorem J: −1). Vector error of the full prediction: 1.995, 1.078, 0.551, 0.277, 0.139 (Theorem J(c): O(1/k_e)).

Numbers are produced by `experiments/exp_periodic_orbit.py`; nothing is edited by hand.
