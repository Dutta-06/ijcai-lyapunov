# Corollary D‴.14 in two dimensions: the real closed loop with the OT leak on a planar plant (`exp_ot_2d.py`)

Plant ẋ = f(x) + u on ℝ², f the sum of four Gaussian bumps with vector amplitudes (A_k, C_k, S_k) = ((0.8, 0.3), (-1.2, 0.8), 0.7), ((-0.5, 0.6), (1, 1), 0.6), ((0.4, -0.7), (0.8, -1.1), 0.6), ((-0.6, -0.4), (-1, -0.9), 0.7), reference x_m(t) = (2 sin t, 1.5 sin 2t), saturated network with n = 2 (a_max = 4.0, c_max = 3.0, s ∈ [0.3, 2.0]), G = diag(1.0, 1.0, 0.25, 0.25, 0.25), matching prior = the exact representation (four atoms with amplitudes 4A_k, masses q = 1/4). Real loop at k_e = 25, σ = σ₀/k_e with σ₀ = 0.05 (ln 2/σ₀ = 14 slow units), RK4 plant with the control and the reference evaluated at the stages and the parameters frozen over the step, Euler parameters, dt = 0.002, 2500 time units = 100 slow units = 5.0/σ₀ = 398 periods; warm start = the mixture of Gaussians of width 0.02 (G⁻¹ metric) around the exact representation with random class labels; the same clouds under the nearest-atom leak and under the OT leak (plan every 100 steps); sizes m = 8, 32 and M = 128; seeds 0, 1; statistics over the last 10 periods.

| seed | leak | m | cluster weights at the start → at the end | within-class spread at the end (rms, whitened) | rms tracking error (last 10 periods) | fit error along the orbit (end) | ‖f̂_m − f̂_M‖²_{L²(orbit)} (end) | W²_{G⁻¹}(μ_m, μ_M) (end) |
|---|---|---|---|---|---|---|---|---|
| 0 | ot | 8 | (0.375, 0.250, 0.250, 0.125) → (0.250, 0.250, 0.250, 0.250) | 1.5e-02 | 1.566e-05 | 8.386e-04 | 1.585e-07 | 4.177e-04 |
| 0 | ot | 32 | (0.219, 0.125, 0.375, 0.281) → (0.250, 0.250, 0.250, 0.250) | 1.8e-02 | 1.258e-05 | 6.763e-04 | 6.625e-08 | 2.684e-04 |
| 0 | ot | 128 | (0.211, 0.234, 0.305, 0.250) → (0.250, 0.250, 0.250, 0.250) | 1.2e-02 | 3.885e-06 | 2.094e-04 | 0.000e+00 | 1.166e-15 |
| 0 | nearest | 8 | (0.375, 0.250, 0.250, 0.125) → (0.375, 0.250, 0.250, 0.125) | 2.8e-04 | 4.466e-03 | 2.577e-01 | 1.850e-02 | 1.174e+01 |
| 0 | nearest | 32 | (0.219, 0.125, 0.375, 0.281) → (0.219, 0.125, 0.375, 0.281) | 3.1e-03 | 4.907e-03 | 2.855e-01 | 8.687e-03 | 7.305e+00 |
| 0 | nearest | 128 | (0.211, 0.234, 0.305, 0.250) → (0.211, 0.234, 0.305, 0.250) | 5.1e-04 | 1.863e-03 | 1.080e-01 | 0.000e+00 | 1.443e-15 |
| 1 | ot | 8 | (0.250, 0.125, 0.125, 0.500) → (0.250, 0.250, 0.250, 0.250) | 1.7e-02 | 1.601e-05 | 8.661e-04 | 1.105e-07 | 4.412e-04 |
| 1 | ot | 32 | (0.219, 0.312, 0.125, 0.344) → (0.250, 0.250, 0.250, 0.250) | 1.5e-02 | 1.275e-05 | 6.938e-04 | 4.977e-08 | 1.972e-04 |
| 1 | ot | 128 | (0.203, 0.281, 0.211, 0.305) → (0.250, 0.250, 0.250, 0.250) | 1.3e-02 | 5.016e-06 | 2.740e-04 | 0.000e+00 | 1.665e-15 |
| 1 | nearest | 8 | (0.250, 0.125, 0.125, 0.500) → (0.250, 0.125, 0.125, 0.500) | 1.3e-03 | 6.747e-03 | 3.924e-01 | 2.310e-02 | 1.179e+01 |
| 1 | nearest | 32 | (0.219, 0.312, 0.125, 0.344) → (0.219, 0.312, 0.125, 0.344) | 5.9e-04 | 4.653e-03 | 2.695e-01 | 5.461e-03 | 4.731e+00 |
| 1 | nearest | 128 | (0.203, 0.281, 0.211, 0.305) → (0.203, 0.281, 0.211, 0.305) | 3.6e-04 | 2.268e-03 | 1.311e-01 | 0.000e+00 | 1.721e-15 |

Equalisation of the OT leak's cluster weights (first weight record from which w = q to within half a particle, in slow units): seed 0, m = 8: 13; seed 0, m = 32: 14; seed 0, m = 128: 14; seed 1, m = 8: 13; seed 1, m = 32: 14; seed 1, m = 128: 14 (ln 2/σ₀ = 14).

Reading. In two dimensions as in one: the nearest-atom leak conserves the clouds' sampled cluster weights and the sizes learn different functions (‖f̂_m − f̂_M‖² = 5.5e-03–2.3e-02), the OT leak drives the weights to the prior's masses and the sizes agree (5.0e-08–1.6e-07, at the stopping time 5/σ₀ and still falling); the rms tracking error is 3.9e-06–1.6e-05 under the OT leak against 1.9e-03–6.7e-03 under the nearest-atom leak, 285–479 times larger per system.

Numbers are produced by `experiments/exp_ot_2d.py`; nothing is edited by hand.
