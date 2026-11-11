# Theorem D″: uniform-in-time propagation of chaos in the leak-dominated regime (`exp_chaos_uniform.py`)

n = 1, k_e = 25, G = diag(0.2, 0.05, 0.05), saturated features with a_max = 4.0, c_max = 3.0, s ∈ [0.3, 2.0], prior ϑ⁰ = 0; plant f(x) = 1.6e^{−(x−1)²/0.4} − 1.2e^{−(x+1)²/0.4} + 0.1x, reference x_m = 1.5 sin t + 0.5 sin 2t; the m-system shares its particles (and its Brownian increments) with an M = 2048 reference; dt = 0.001.

Constants. Lipschitz constants of the features in closed form (suprema over all ϑ and x): L_x = L_ϑ = a_max e^{−1/2}/s_min = 8.087, L₂ = L_xϑ = a_max/s_min² = 44.44; numerical search over 1.2·10⁵ sampled ϑ, |x| ≤ 2.5, Nelder–Mead polish: 8.087, 8.087, 44.44, 44.44. L_f = sup|f′| = 2.293 on |x| ≤ 2.5. A-priori tube (Theorem D(b), e(0) = 0): sup|x_m| = 1.7425, F = sup_tube|f| = 1.7006, E_X = (F + K_Φ)/k_e = 0.2280 (measured sup‖ē‖ on the τ = 0 reference runs: 0.0697; the runs stay within |x| ≤ 1.77 ≤ 2.5). κ_e = k_e − L_f − L_x = 14.620; L_b′ = λ_G(L_ϑ + E_X L_xϑ) = 3.644; the theorem applies (a priori) for σ > σ_cert = λ_G E_X L₂ + L_b′L_ϑ/κ_e = 4.043. Proposition K: the collapsed critical point (best single kernel, a = 0.136, c = 0.068, s = 1.150) is a strict local minimum of J_m iff σ > σ_K = σ₀_crit/k_e = 0.00144 (σ₀_crit = 0.03610); below σ_K the kernels separate: the learning regime. σ_cert/σ_K = 2800.

## σ = 4.8513, τ = 0.02: certified, dither τ = 0.02 (certified: det A = 11.821, κ_ϑ = 2.824); T = 80 s, 8 seeds

Reference network at T: particle spread in (a, c, s), max 1.07e-01 and rms 3.26e-02 (rms at t = 0: 1.53); relative fit error ‖f − f̂‖/‖f‖ on the orbit 1.000 (at t = 0: 1.001).

| m | max D_m on (10, 20] | max D_m on (20, 40] | max D_m on (40, 80] | time average on (40, 80] | bound K_Φ(κ_ϑ + L_b′)m^{−1/2}/det A |
|---|---|---|---|---|---|
| 16 | 6.670e-04 | 6.967e-04 | 7.855e-04 | 1.457e-04 | 5.472e-01 |
| 64 | 3.253e-04 | 3.577e-04 | 3.751e-04 | 7.247e-05 | 2.736e-01 |
| 256 | 1.645e-04 | 1.696e-04 | 1.917e-04 | 3.404e-05 | 1.368e-01 |

Slope of the last-window maximum against m (per-seed fits, mean ± s.e.): -0.51 ± 0.02 (Theorem D″: −0.5). Ratio of the last to the first window maxima: 1.18, 1.15, 1.17 (uniform in time means no growth). The bound column is sup_t E D_m; the table reports E max_t D_m over the window, which is larger.

## σ = 4.8513, τ = 0: certified, τ = 0 (certified: det A = 11.821, κ_ϑ = 2.824); T = 80 s, 8 seeds

Reference network at T: particle spread in (a, c, s), max 1.69e-14 and rms 1.69e-14 (rms at t = 0: 1.53); relative fit error ‖f − f̂‖/‖f‖ on the orbit 1.000 (at t = 0: 1.001).

| m | max D_m on (10, 20] | max D_m on (20, 40] | max D_m on (40, 80] | time average on (40, 80] | bound K_Φ(κ_ϑ + L_b′)m^{−1/2}/det A |
|---|---|---|---|---|---|
| 16 | 2.504e-25 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 5.472e-01 |
| 64 | 6.099e-26 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 2.736e-01 |
| 256 | 1.232e-26 | 0.000e+00 | 0.000e+00 | 0.000e+00 | 1.368e-01 |

Slope of the last-window maximum against m (per-seed fits, mean ± s.e.): — (zero to rounding) (Theorem D″: −0.5). Ratio of the last to the first window maxima: 0.00, 0.00, 0.00 (uniform in time means no growth). The bound column is sup_t E D_m; the table reports E max_t D_m over the window, which is larger.

## σ = 0.5053, τ = 0: σ_cert/8, τ = 0 (below the threshold); T = 80 s, 8 seeds

Reference network at T: particle spread in (a, c, s), max 3.18e-14 and rms 3.18e-14 (rms at t = 0: 1.54); relative fit error ‖f − f̂‖/‖f‖ on the orbit 1.000 (at t = 0: 1.001).

| m | max D_m on (10, 20] | max D_m on (20, 40] | max D_m on (40, 80] | time average on (40, 80] | bound K_Φ(κ_ϑ + L_b′)m^{−1/2}/det A |
|---|---|---|---|---|---|
| 16 | 7.246e-05 | 5.063e-07 | 4.306e-11 | 2.235e-12 | nan |
| 64 | 5.942e-05 | 4.142e-07 | 3.471e-11 | 1.800e-12 | nan |
| 256 | 3.077e-05 | 2.140e-07 | 1.786e-11 | 9.262e-13 | nan |

Slope of the last-window maximum against m (per-seed fits, mean ± s.e.): -0.33 ± 0.20 (Theorem D″: −0.5). Ratio of the last to the first window maxima: 0.00, 0.00, 0.00 (uniform in time means no growth). The bound column is sup_t E D_m; the table reports E max_t D_m over the window, which is larger.

## σ = 0.0005, τ = 0: σ_K/3, τ = 0 (learning regime) (below the threshold); T = 1000 s, 6 seeds

Reference network at T: particle spread in (a, c, s), max 3.39e+00 and rms 1.53e+00 (rms at t = 0: 1.54); relative fit error ‖f − f̂‖/‖f‖ on the orbit 0.632 (at t = 0: 1.000).

| m | max D_m on (100, 200] | max D_m on (200, 500] | max D_m on (500, 1000] | time average on (500, 1000] | bound K_Φ(κ_ϑ + L_b′)m^{−1/2}/det A |
|---|---|---|---|---|---|
| 16 | 6.054e-02 | 1.224e-01 | 1.991e-01 | 1.601e-01 | nan |
| 64 | 2.989e-02 | 5.858e-02 | 9.216e-02 | 7.498e-02 | nan |
| 256 | 1.325e-02 | 2.475e-02 | 3.858e-02 | 3.162e-02 | nan |

Slope of the last-window maximum against m (per-seed fits, mean ± s.e.): -0.62 ± 0.16 (Theorem D″: −0.5). Ratio of the last to the first window maxima: 3.29, 3.08, 2.91 (uniform in time means no growth). The bound column is sup_t E D_m; the table reports E max_t D_m over the window, which is larger.

Reading: at the certified leak the particle map is a contraction uniformly in time. Without dither the reference collapses to one kernel and the coupled distance vanishes to rounding; with dither the particles form the cloud of the prior and the coupled distance is of order m^{−1/2}, bounded in time and below the bound — the theorem's content at τ > 0. In both cases the network does not learn (its fit error stays at that of the empty network). At σ_cert/8 the particles still synchronise, more slowly, and the coupled distance decays with them. At σ_K/3 the collapsed configuration is unstable, the kernels stay separated (rms spread 1.53) and the network fits the plant (fit error 0.632): the learning regime, where the coupled distance follows the m^{−1/2} law at each time and grows over the horizon (last-to-first window ratios 3.29, 3.08, 2.91): on this horizon propagation of chaos in the learning regime is the finite-horizon statement of Theorem D′, not a uniform one; whether the growth saturates on longer horizons is open.

Numbers are produced by `experiments/exp_chaos_uniform.py`; nothing is edited by hand.
