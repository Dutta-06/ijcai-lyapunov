# Proposition D‴.45: the free mass–amplitude law is a mirror flow (`exp_mirror_free.py`)

Free masses (λ ≡ 0): dα_i = Gβ_ia′_ig_i, dβ_i = γβ_ia_ig_i, a = a_max tanh(α/a_max), w = β⊙a/m; m = 12, G = 1, γ = 0.1, a_max = 1 (c = γa_max²/(4G) = 0.025) unless said otherwise; Q_i = β_i − c cosh(2α_i/a_max); on a dictionary the mirror flow dη = Φᵀ(r + ξ)dt/N of Ψ_q, ψ′ = η = ∫dα/(Gβa′) (`mirror_free`).

## A. Identities

η's closed form against quadrature (q = 100c, 3, 0.5, c, c/2, 0, −c/2 and the branch q = −2c, 17 and 5 points): largest relative difference 1.51e-14; (ψ*)′ = w against finite differences 7.15e-11; dw/dη = mD 4.28e-10.

Q along the loop (5 × 12, 2-sparse truths, β(0) = 0.3, α(0) ~ N(0, 0.05²)), largest |Q(t) − Q(0)| over 60 times; the conserved law's shift: ΔQ_i along its integrated trajectory (t ≤ 200, the smooth noise) against the quadrature of −γβ_iλ, λ from the state (relative):

| instance | noiseless, t ≤ 10⁴ | a smooth noise (30 sinusoids, frequencies 10⁻³–3, mixed over the nodes, amplitude 0.3), t ≤ 3·10⁴ | the conserved law's shift |
|---|---|---|---|
| 0 | 2.04e-12 | 2.89e-10 | 4.44e-7 |
| 1 | 1.40e-12 | 2.50e-10 | 1.45e-7 |
| 2 | 1.41e-11 | 3.61e-10 | 3.05e-7 |

Q's drift along Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's update; at the lazy projection) against the Euler step's closed form (2|φ_i|²/(mN))γβ_i(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²)) per unit s: largest relative difference 5.27e-16, 9.11e-16, 7.88e-16.

## B. The implicit bias: the start's masses set the regime (5 × 12, 20 instances per k)

From α = 0 and equal masses β = s·m/d: the Bregman projection (Newton on the dual; a projection whose residual exceeds 1e-09(1 + |y|) is a failure, counted and left out of every aggregate) — ‖w‖₁ over basis pursuit's, the atoms used (|w| > 10⁻³), the fraction recovering the truth (10⁻³·max|w*|), the median and largest distances to basis pursuit, the median relative distance to the minimum-norm fit Φ⁺y, the largest |a|/a_max and β; for s ≥ 0.01 the noiseless loop's end at t = 1e+07 against the projection (its largest residual and distance).

| k | s | failures | ‖w‖₁/‖w_bp‖₁ | used | the truth | ‖w − w_bp‖∞: median, largest | median ‖w − Φ⁺y‖∞/‖Φ⁺y‖∞ | max |a|/a_max, max β | the loop's end: residual, distance to the projection |
|---|---|---|---|---|---|---|---|---|---|
| 2 | basis pursuit | | 1 | | 0.80 | | | | |
| 2 | 10 | 0 | 1.3715 | 11.8 | 0.00 | 0.14, 0.206 | 0.00499 | 0.241, 10 | 1.14e-16, 2.03e-12 |
| 2 | 1 | 0 | 1.5248 | 11.8 | 0.00 | 0.159, 0.212 | 0.197 | 0.979, 2.13 | 8.44e-17, 1.90e-12 |
| 2 | 0.1 | 0 | 1.1662 | 11.4 | 0.00 | 0.0973, 0.192 | 0.408 | 0.993, 3.76 | 5.69e-15, 8.51e-12 |
| 2 | 0.01 | 0 | 1.0274 | 8.1 | 0.00 | 0.0235, 0.182 | 0.814 | 0.994, 3.86 | 4.79e-15, 4.74e-11 |
| 2 | 0.001 | 0 | 1.0022 | 3.6 | 0.15 | 9.04e-4, 0.154 | 1.03 | 0.993, 3.58 | — |
| 2 | 0.0001 | 0 | 1.0004 | 3.0 | 0.75 | 1.38e-5, 0.0697 | 1.16 | 0.993, 3.59 | — |
| 3 | basis pursuit | | 1 | | 0.30 | | | | |
| 3 | 10 | 0 | 1.2770 | 11.9 | 0.00 | 0.111, 0.179 | 0.00548 | 0.331, 10 | 1.45e-16, 2.59e-12 |
| 3 | 1 | 0 | 1.4152 | 11.9 | 0.00 | 0.13, 0.19 | 0.187 | 0.989, 3.16 | 8.36e-17, 3.82e-12 |
| 3 | 0.1 | 0 | 1.1022 | 11.6 | 0.00 | 0.0659, 0.146 | 0.322 | 0.994, 4.18 | 8.78e-17, 5.68e-12 |
| 3 | 0.01 | 0 | 1.0175 | 8.0 | 0.00 | 0.0286, 0.109 | 0.587 | 0.994, 3.92 | 7.34e-17, 4.10e-11 |
| 3 | 0.001 | 0 | 1.0019 | 5.6 | 0.05 | 0.0041, 0.0629 | 0.727 | 0.993, 3.56 | — |
| 3 | 0.0001 | 0 | 1.0001 | 4.5 | 0.30 | 1.30e-4, 0.0217 | 0.747 | 0.993, 3.67 | — |

The loop at s = 0.001 (k = 2, instance 0): distance to its projection 0.277 at t = 1e+02, 0.277 at t = 1e+03, 0.00643 at t = 2e+04, 7.10e-4 at t = 2e+05, 1.50e-5 at t = 2e+06, 1.02e-12 at t = 3e+07; below 10⁻⁸ from t ≈ 8.5e+06.

Small random amplitudes (α(0) ~ N(0, 0.05²)) at tiny masses: the particles with Q < −c sit on branches (Proposition D‴.47) — the sign-constrained fits S, the coordinates vanishing on S (|Z|), the instances with a truth atom on a wrong-signed branch, whether the frozen signs still fit y exactly, and the distribution of the projection's distance to basis pursuit:

| k | s | branches (mean) | failures | exact fit kept | a truth atom on a wrong-signed branch | |Z| (mean) | ‖w − w_bp‖∞: ≤ 10⁻³, median, largest |
|---|---|---|---|---|---|---|---|
| 2 | 0.01 | 0.0 | 0 | 20/20 | 0/20 | 0.0 | 0/20, 0.0268, 0.132 |
| 2 | 0.0001 | 4.2 | 0 | 20/20 | 5/20 | 0.0 | 10/20, 6.58e-4, 0.37 |
| 2 | 1e-06 | 11.1 | 0 | 18/20 | 13/20 | 0.8 | 6/20, 0.252, 16.3 |
| 3 | 0.01 | 0.0 | 0 | 20/20 | 0/20 | 0.0 | 1/20, 0.0236, 0.215 |
| 3 | 0.0001 | 4.2 | 0 | 20/20 | 5/20 | 0.0 | 4/20, 0.147, 0.572 |
| 3 | 1e-06 | 11.4 | 0 | 19/20 | 18/20 | 0.4 | 1/20, 0.357, 0.727 |

Unequal masses β⁰_i = s·r_i (r_i log-uniform in [1/4, 4]), α = 0: the projection against weighted basis pursuit (min Σ|w_i|/√r_i) as s → 0 — and against the unweighted one — and against the weighted minimum norm (min Σw_i²/r_i²) as s → ∞ (relative); median and largest over the converged projections (failures in brackets):

| k | s = 0.001: weighted (median, largest), unweighted (median) | s = 1e-05: weighted (median, largest), unweighted (median) | s = 1e-07: weighted (median, largest), unweighted (median) | s = 100: weighted ℓ² (median, largest) | s = 10000: weighted ℓ² (median, largest) |
|---|---|---|---|---|---|
| 2 | 0.00622, 0.0892, 0.0203 [0] | 2.33e-6, 0.00374, 2.70e-5 [0] | 1.74e-9, 6.56e-6, 1.53e-8 [0] | 2.44e-5, 3.20e-4 [0] | 2.44e-9, 3.20e-8 [0] |
| 3 | 0.0128, 0.19, 0.0833 [0] | 3.04e-5, 0.17, 0.084 [0] | 2.06e-8, 0.169, 0.084 [0] | 2.52e-5, 2.90e-4 [0] | 2.52e-9, 2.90e-8 [0] |

## C. Noise

(1) White label noise at the intensity σ = 0.3 (dη_i = (Φᵀr)_idt/N + σ(ΦᵀdW)_i/N), 64 paths to t = 200 from random states (α ~ N(0, 0.3²), β ~ U(0.2, 2)): Q's mean change against the mean time integral of the Itô correction (γβ_i/2)(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²))σ²|φ_i|²/N² along the same paths, for three schemes as the step h → 0:

| instance | h | scheme | mean ΔQ | the Itô integral | ratio | max |ΔQ| |
|---|---|---|---|---|---|---|
| 0 | 0.1 | Euler (Itô) | -2.5763e-02 | -2.5493e-02 | 1.0106 | 0.137 |
| 0 | 0.1 | Heun (Stratonovich) | +1.3843e-05 | -2.7932e-02 | -0.0005 | 1.48e-4 |
| 0 | 0.1 | mirror | -1.2649e-18 | -2.7924e-02 | 0.0000 | 2.22e-16 |
| 0 | 0.03 | Euler (Itô) | -2.5681e-02 | -2.5586e-02 | 1.0037 | 0.134 |
| 0 | 0.03 | Heun (Stratonovich) | +4.1160e-06 | -2.8018e-02 | -0.0001 | 4.72e-5 |
| 0 | 0.03 | mirror | +7.2280e-19 | -2.8016e-02 | -0.0000 | 2.22e-16 |
| 0 | 0.01 | Euler (Itô) | -2.5686e-02 | -2.5658e-02 | 1.0011 | 0.132 |
| 0 | 0.01 | Heun (Stratonovich) | +1.3599e-06 | -2.8071e-02 | -0.0000 | 1.32e-5 |
| 0 | 0.01 | mirror | +1.8793e-18 | -2.8071e-02 | -0.0000 | 2.22e-16 |
| 1 | 0.1 | Euler (Itô) | -3.2788e-02 | -3.2635e-02 | 1.0047 | 0.0876 |
| 1 | 0.1 | Heun (Stratonovich) | +9.6836e-06 | -3.4582e-02 | -0.0003 | 1.31e-4 |
| 1 | 0.1 | mirror | +0.0000e+00 | -3.4577e-02 | -0.0000 | 2.22e-16 |
| 1 | 0.03 | Euler (Itô) | -3.2848e-02 | -3.2799e-02 | 1.0015 | 0.089 |
| 1 | 0.03 | Heun (Stratonovich) | +3.0278e-06 | -3.4780e-02 | -0.0001 | 4.13e-5 |
| 1 | 0.03 | mirror | +2.8912e-19 | -3.4779e-02 | -0.0000 | 2.22e-16 |
| 1 | 0.01 | Euler (Itô) | -3.2802e-02 | -3.2781e-02 | 1.0006 | 0.0919 |
| 1 | 0.01 | Heun (Stratonovich) | +9.7190e-07 | -3.4733e-02 | -0.0000 | 1.31e-5 |
| 1 | 0.01 | mirror | +5.0596e-19 | -3.4733e-02 | -0.0000 | 2.22e-16 |
| 2 | 0.1 | Euler (Itô) | -4.3296e-02 | -4.2928e-02 | 1.0086 | 0.123 |
| 2 | 0.1 | Heun (Stratonovich) | +1.8664e-05 | -4.6228e-02 | -0.0004 | 1.10e-4 |
| 2 | 0.1 | mirror | -3.9754e-19 | -4.6218e-02 | 0.0000 | 2.22e-16 |
| 2 | 0.03 | Euler (Itô) | -4.2979e-02 | -4.2810e-02 | 1.0039 | 0.117 |
| 2 | 0.03 | Heun (Stratonovich) | +5.6768e-06 | -4.6109e-02 | -0.0001 | 4.40e-5 |
| 2 | 0.03 | mirror | +0.0000e+00 | -4.6106e-02 | -0.0000 | 2.22e-16 |
| 2 | 0.01 | Euler (Itô) | -4.2771e-02 | -4.2683e-02 | 1.0021 | 0.117 |
| 2 | 0.01 | Heun (Stratonovich) | +1.8217e-06 | -4.5964e-02 | -0.0000 | 1.50e-5 |
| 2 | 0.01 | mirror | +1.9877e-18 | -4.5963e-02 | -0.0000 | 2.22e-16 |

(2) The Euler step on (α, log β) with label noise (σ = 0.3, at random states): its exact mean change of Q (Gauss–Hermite, 80 nodes, over each particle's Gaussian drive) against the closed form (h²E[g²]γβ/2)(γa² − Gβa′(1 + a²/a_max²)) — largest relative difference over the particles, three instances per h: h = 0.5: 0.0324, 0.0326, 0.036; h = 0.05: 0.00288, 0.00341, 0.00354; h = 0.005: 2.85e-4, 3.43e-4, 3.55e-4.

(3) The selection in the mirror coordinates: Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's update) from the lazy projection (s = 1) — ‖w‖₁, the distance to basis pursuit, the atoms used, the dead atoms (Q ≤ −c), Q's mean and range; at the end the support against basis pursuit's:

| instance (basis pursuit = truth) | s | ‖w‖₁ | ‖w − w_bp‖∞ | used | dead | Q: mean, min, max |
|---|---|---|---|---|---|---|
| 0 (yes) | 0e+00 | 0.7060 | 0.196 | 12 | 0 | 0.975, 0.975, 0.975 |
| 0 (yes) | 1e+02 | 0.7689 | 0.212 | 12 | 0 | 0.7922, 0.7497, 0.8296 |
| 0 (yes) | 1e+03 | 0.7813 | 0.215 | 12 | 0 | 0.3792, 0.3301, 0.4388 |
| 0 (yes) | 1e+04 | 0.7683 | 0.211 | 12 | 0 | 0.01402, 0.0121, 0.01648 |
| 0 (yes) | 1e+05 | 0.6146 | 0.154 | 12 | 0 | -0.001371, -0.003571, -0.000425 |
| 0 (yes) | 1e+06 | 0.4563 | 1.58e-6 | 2 | 0 | -0.02056, -0.02466, -0.0001863 |
| 0 (yes) | 1e+07 | 0.4563 | 4.80e-8 | 2 | 0 | -0.02085, -0.02497, -0.0001863 |
| 0 (yes) | 1e+08 | 0.4563 | 6.54e-9 | 2 | 0 | -0.02087, -0.025, -0.0001863 |
| 1 (yes) | 0e+00 | 0.5690 | 0.123 | 12 | 0 | 0.975, 0.975, 0.975 |
| 1 (yes) | 1e+02 | 0.6770 | 0.132 | 12 | 0 | 0.7806, 0.725, 0.8237 |
| 1 (yes) | 1e+03 | 0.7465 | 0.125 | 12 | 0 | 0.365, 0.2347, 0.4085 |
| 1 (yes) | 1e+04 | 0.7173 | 0.124 | 12 | 0 | 0.01347, 0.008037, 0.01562 |
| 1 (yes) | 1e+05 | 0.5053 | 0.0945 | 10 | 0 | -0.004291, -0.0195, -0.0004748 |
| 1 (yes) | 1e+06 | 0.3863 | 1.17e-6 | 2 | 0 | -0.02058, -0.02468, -0.0002576 |
| 1 (yes) | 1e+07 | 0.3863 | 4.27e-8 | 2 | 0 | -0.02085, -0.02497, -0.0002576 |
| 1 (yes) | 1e+08 | 0.3863 | 7.96e-9 | 2 | 0 | -0.02088, -0.025, -0.0002576 |
| 2 (yes) | 0e+00 | 0.7899 | 0.178 | 12 | 0 | 0.975, 0.975, 0.975 |
| 2 (yes) | 1e+02 | 0.8336 | 0.188 | 12 | 0 | 0.808, 0.7331, 0.874 |
| 2 (yes) | 1e+03 | 0.8541 | 0.192 | 12 | 0 | 0.3973, 0.307, 0.4913 |
| 2 (yes) | 1e+04 | 0.8382 | 0.189 | 12 | 0 | 0.01476, 0.01121, 0.0187 |
| 2 (yes) | 1e+05 | 0.6570 | 0.135 | 12 | 0 | -0.001433, -0.003571, -0.0003782 |
| 2 (yes) | 1e+06 | 0.4894 | 1.01e-6 | 2 | 0 | -0.02056, -0.02466, -0.0001892 |
| 2 (yes) | 1e+07 | 0.4894 | 3.55e-8 | 2 | 0 | -0.02084, -0.02497, -0.0001892 |
| 2 (yes) | 1e+08 | 0.4894 | 4.54e-9 | 2 | 0 | -0.02087, -0.025, -0.0001892 |
| 3 (yes) | 0e+00 | 0.7513 | 0.169 | 12 | 0 | 0.975, 0.975, 0.975 |
| 3 (yes) | 1e+02 | 0.8431 | 0.17 | 12 | 0 | 0.8125, 0.7251, 0.8806 |
| 3 (yes) | 1e+03 | 0.8939 | 0.174 | 12 | 0 | 0.407, 0.2752, 0.4858 |
| 3 (yes) | 1e+04 | 0.8651 | 0.173 | 12 | 0 | 0.01516, 0.00992, 0.01839 |
| 3 (yes) | 1e+05 | 0.6050 | 0.158 | 11 | 0 | -0.00311, -0.01572, -0.0003515 |
| 3 (yes) | 1e+06 | 0.4857 | 1.41e-6 | 2 | 0 | -0.02055, -0.02467, -0.0002009 |
| 3 (yes) | 1e+07 | 0.4857 | 4.65e-8 | 2 | 0 | -0.02084, -0.02497, -0.0002009 |
| 3 (yes) | 1e+08 | 0.4857 | 6.88e-9 | 2 | 0 | -0.02087, -0.025, -0.0002009 |
| 4 (yes) | 0e+00 | 0.7598 | 0.169 | 12 | 0 | 0.975, 0.975, 0.975 |
| 4 (yes) | 1e+02 | 0.8381 | 0.173 | 12 | 0 | 0.8126, 0.725, 0.8789 |
| 4 (yes) | 1e+03 | 0.8714 | 0.173 | 12 | 0 | 0.4032, 0.2126, 0.4692 |
| 4 (yes) | 1e+04 | 0.8430 | 0.17 | 12 | 0 | 0.01492, 0.006535, 0.01755 |
| 4 (yes) | 1e+05 | 0.5881 | 0.136 | 11 | 0 | -0.003355, -0.02053, -0.000337 |
| 4 (yes) | 1e+06 | 0.4273 | 9.99e-7 | 2 | 0 | -0.02058, -0.02469, -0.0001962 |
| 4 (yes) | 1e+07 | 0.4273 | 2.26e-8 | 2 | 0 | -0.02085, -0.02497, -0.0001962 |
| 4 (yes) | 1e+08 | 0.4273 | 2.67e-9 | 2 | 0 | -0.02087, -0.025, -0.0001962 |
| 5 (yes) | 0e+00 | 0.7133 | 0.187 | 12 | 0 | 0.975, 0.975, 0.975 |
| 5 (yes) | 1e+02 | 0.7984 | 0.201 | 12 | 0 | 0.7976, 0.7604, 0.843 |
| 5 (yes) | 1e+03 | 0.8123 | 0.205 | 12 | 0 | 0.3882, 0.3354, 0.4382 |
| 5 (yes) | 1e+04 | 0.7954 | 0.2 | 12 | 0 | 0.01438, 0.01229, 0.01636 |
| 5 (yes) | 1e+05 | 0.5995 | 0.123 | 12 | 0 | -0.001489, -0.004242, -0.0003432 |
| 5 (yes) | 1e+06 | 0.4299 | 1.17e-6 | 2 | 0 | -0.02057, -0.02466, -0.000187 |
| 5 (yes) | 1e+07 | 0.4299 | 5.98e-8 | 2 | 0 | -0.02085, -0.02497, -0.000187 |
| 5 (yes) | 1e+08 | 0.4299 | 1.16e-8 | 2 | 0 | -0.02087, -0.025, -0.000187 |
| 8 (no) | 0e+00 | 0.5493 | 0.0607 | 11 | 0 | 0.975, 0.975, 0.975 |
| 8 (no) | 1e+02 | 0.6840 | 0.0683 | 12 | 0 | 0.7799, 0.7252, 0.8183 |
| 8 (no) | 1e+03 | 0.7443 | 0.0757 | 12 | 0 | 0.3655, 0.2196, 0.4058 |
| 8 (no) | 1e+04 | 0.7131 | 0.0703 | 12 | 0 | 0.01346, 0.007105, 0.01507 |
| 8 (no) | 1e+05 | 0.4796 | 0.0405 | 11 | 0 | -0.003941, -0.02056, -0.000483 |
| 8 (no) | 1e+06 | 0.4228 | 0.0993 | 5 | 0 | -0.01507, -0.02463, -0.0003378 |
| 8 (no) | 1e+07 | 0.4228 | 0.0993 | 5 | 0 | -0.01518, -0.02497, -0.0003378 |
| 8 (no) | 1e+08 | 0.4228 | 0.0993 | 5 | 0 | -0.01519, -0.02499, -0.0003378 |
| 9 (no) | 0e+00 | 0.3881 | 0.0722 | 12 | 0 | 0.975, 0.975, 0.975 |
| 9 (no) | 1e+02 | 0.5189 | 0.0774 | 12 | 0 | 0.748, 0.725, 0.7771 |
| 9 (no) | 1e+03 | 0.5856 | 0.082 | 12 | 0 | 0.3176, 0.2489, 0.3837 |
| 9 (no) | 1e+04 | 0.5587 | 0.0809 | 12 | 0 | 0.01158, 0.008751, 0.01441 |
| 9 (no) | 1e+05 | 0.3742 | 0.0571 | 10 | 0 | -0.004401, -0.01804, -0.0007273 |
| 9 (no) | 1e+06 | 0.3277 | 2.76e-5 | 5 | 0 | -0.01479, -0.02467, -0.0003997 |
| 9 (no) | 1e+07 | 0.3277 | 2.59e-5 | 5 | 0 | -0.01491, -0.02497, -0.0003997 |
| 9 (no) | 1e+08 | 0.3277 | 2.58e-5 | 5 | 0 | -0.01492, -0.02498, -0.0003997 |
| 13 (no) | 0e+00 | 0.5207 | 0.142 | 12 | 0 | 0.975, 0.975, 0.975 |
| 13 (no) | 1e+02 | 0.6726 | 0.152 | 12 | 0 | 0.7767, 0.7251, 0.8387 |
| 13 (no) | 1e+03 | 0.7601 | 0.155 | 12 | 0 | 0.3686, 0.2712, 0.4174 |
| 13 (no) | 1e+04 | 0.7254 | 0.153 | 12 | 0 | 0.01364, 0.009735, 0.01577 |
| 13 (no) | 1e+05 | 0.4892 | 0.13 | 9 | 0 | -0.004545, -0.01746, -0.0005203 |
| 13 (no) | 1e+06 | 0.4714 | 0.275 | 2 | 0 | -0.01944, -0.02453, -0.0001877 |
| 13 (no) | 1e+07 | 0.4714 | 0.275 | 2 | 0 | -0.02084, -0.02497, -0.0001877 |
| 13 (no) | 1e+08 | 0.4714 | 0.275 | 2 | 0 | -0.02087, -0.025, -0.0001877 |
| 17 (no) | 0e+00 | 0.7098 | 0.167 | 12 | 0 | 0.975, 0.975, 0.975 |
| 17 (no) | 1e+02 | 0.7890 | 0.186 | 12 | 0 | 0.7978, 0.7301, 0.8413 |
| 17 (no) | 1e+03 | 0.8108 | 0.192 | 12 | 0 | 0.3868, 0.3109, 0.4466 |
| 17 (no) | 1e+04 | 0.7961 | 0.187 | 12 | 0 | 0.01433, 0.01136, 0.01675 |
| 17 (no) | 1e+05 | 0.6319 | 0.111 | 12 | 0 | -0.002056, -0.01115, -0.0003643 |
| 17 (no) | 1e+06 | 0.5075 | 6.58e-4 | 5 | 0 | -0.01495, -0.02467, -0.0002019 |
| 17 (no) | 1e+07 | 0.5075 | 6.50e-4 | 5 | 0 | -0.01508, -0.02497, -0.000202 |
| 17 (no) | 1e+08 | 0.5075 | 6.50e-4 | 5 | 0 | -0.01508, -0.025, -0.000202 |

The end supports: instance 0: [2, 3] against basis pursuit's [2, 3] (‖w‖₁ 0.4563 against 0.4563) — equal; instance 1: [1, 7] against basis pursuit's [1, 7] (‖w‖₁ 0.3863 against 0.3863) — equal; instance 2: [7, 10] against basis pursuit's [7, 10] (‖w‖₁ 0.4894 against 0.4894) — equal; instance 3: [7, 9] against basis pursuit's [7, 9] (‖w‖₁ 0.4857 against 0.4857) — equal; instance 4: [1, 7] against basis pursuit's [1, 7] (‖w‖₁ 0.4273 against 0.4273) — equal; instance 5: [0, 2] against basis pursuit's [0, 2] (‖w‖₁ 0.4299 against 0.4299) — equal; instance 8: [1, 5, 7, 8, 9] against basis pursuit's [2, 7, 8, 9, 11] (‖w‖₁ 0.4228 against 0.4156) — different; instance 9: [4, 6, 8, 9, 11] against basis pursuit's [4, 6, 8, 9, 11] (‖w‖₁ 0.3277 against 0.3277) — equal; instance 13: [2, 4] against basis pursuit's [0, 2, 6, 7, 8] (‖w‖₁ 0.4714 against 0.4516) — different; instance 17: [1, 3, 9, 10, 11] against basis pursuit's [1, 3, 9, 10, 11] (‖w‖₁ 0.5075 against 0.5070) — equal.

(4) The sampled-data loop at the step h = 0.5 with label noise drawn afresh (σ = 0.3, 16 paths, 20000 steps) from the lazy projection (α = 0, β = m/d): the paradigm's Euler update on (α, log β) against the mirror update η⁺ = η + hg on each particle's curve — Q's change (mean over the particles and paths, largest), the projection of the end state: its distance to the start's (mean over the paths), to basis pursuit's, its ‖w‖₁; the end's residual:

| instance | the start: ‖w‖₁, ‖w − w_bp‖∞ | the update | Q − Q(0): mean, max | ‖P − P(0)‖∞ | ‖P − w_bp‖∞ | ‖P‖₁ | |Φw − y| |
|---|---|---|---|---|---|---|---|
| 0 | 0.7060, 0.196 | Euler | -3.379e-01, 0.411 | 0.03 | 0.214 | 0.7782 | 0.021 |
| 0 | 0.7060, 0.196 | mirror | -2.313e-18, 1.11e-16 | 5.86e-12 | 0.196 | 0.7060 | 0.046 |
| 1 | 0.5690, 0.123 | Euler | -3.559e-01, 0.477 | 0.0429 | 0.132 | 0.7139 | 0.022 |
| 1 | 0.5690, 0.123 | mirror | -3.469e-18, 1.11e-16 | 8.98e-12 | 0.123 | 0.5690 | 0.048 |

Proposition D‴.44's reduced flow (the step's filter) to the same s = (m v/4)t: its ‖w‖₁ and mean change of Q, and the Euler paths' mean projection against its state: instance 0 (s = 270): 0.7793, -3.381e-01, 7.42e-4; instance 1 (s = 270): 0.7232, -3.552e-01, 0.00568.

## D. The real loop (k_e = 5, γ = 0.5, m = 6, six Gaussian atoms of width 0.4 on [−1.2, 1.2] and a second particle on the third, without cohesion between them; x_m = sin t + 0.4 sin 2.7t; 2000 time units)

The target's coefficient on the shared atom w*₃ = 0.05111, at the pair's Bregman split 0.03685 + 0.01426 (V's reference).

| the noise | V(0) | V(T) | V monotone | max |V̇ + k_ee²|/max k_ee² | the pair's mirror coordinates (max |Δη_i − Δη_j|, relative) | max |Q − Q(0)| | max |e| (last 10%) | the distinct atoms' |w − w*| | the pair's sum's error |
|---|---|---|---|---|---|---|---|---|---|
| noiseless | 0.175 | 0.022 | yes | 2.95e-9 | 2.39e-15 | 7.55e-15 | 0.00544 | 0.048 | 0.0272 |
| a smooth noise in the adaptation only | 0.175 | 0.136 | no | — | 6.26e-13 | 6.82e-12 | 0.148 | 0.247 | 0.0234 |
| a sensor noise (feedback, regressor, adaptation) | 0.175 | 0.377 | no | — | 9.20e-14 | 4.90e-13 | 0.868 | 0.255 | 0.139 |

## E. The conserved law (Σβ = m) from the same start (α = 0, β = m/d; t = 10⁶)

| instance | its residual | max |Q − Q(0)| | ‖w‖₁: conserved, free (the projection) | the ends' distance | Σβ |
|---|---|---|---|---|---|
| 0 | 4.16e-17 | 0.0817 | 0.7097, 0.7060 | 0.00158 | 12.000000 |
| 1 | 1.43e-17 | 0.0584 | 0.5721, 0.5690 | 0.00177 | 12.000000 |
| 2 | 6.80e-17 | 0.186 | 0.7886, 0.7899 | 0.00219 | 12.000000 |
| 3 | 4.13e-17 | 0.224 | 0.7454, 0.7513 | 0.00409 | 12.000000 |
| 4 | 4.86e-17 | 0.191 | 0.7529, 0.7598 | 0.00388 | 12.000000 |
| 5 | 7.47e-17 | 0.0984 | 0.7181, 0.7133 | 0.00209 | 12.000000 |

## F. The rich limit's leading order (unequal masses εr_i, the instances of B; Proposition D‴.45(iii))

Where weighted basis pursuit (weights r^{−1/2}) is a vertex with N atoms: its certificate's smallest relative slack μ_min off the support; the projection from α = 0 (exact fallback) — its distance |w − w_bp|∞ — against the leading order computed from the certificate (|w − w₁|∞, absolute: the rounding floor is ≈ 10⁻¹²) and against the closed form a_max(εr_j)^{3/2}cos(πμ_j/2)/(m(2c)^{1/2}sin³(πμ_j/2)) off the support (its largest relative error); medians and largest values over the instances with μ_min ≥ 0.01:

k = 2: 10 of 20 instances with an N-atom vertex, 9 with μ_min ≥ 0.01; μ_min median 0.0764, smallest 0.00491.

| ε | distance: median, largest | |w − w₁|∞: median, largest | the closed form: median, largest | projection: largest residual |
|---|---|---|---|---|
| 1e-05 | 3.37e-5, 7.22e-4 | 4.04e-9, 2.83e-6 | 1.09, 11.4 | 1.13e-12 |
| 1e-07 | 1.56e-8, 1.12e-7 | 4.55e-14, 1.81e-13 | 0.0758, 1.13 | 2.61e-13 |
| 1e-09 | 1.45e-11, 9.48e-11 | 7.61e-13, 1.08e-12 | 0.00731, 0.113 | 1.17e-12 |

The distance's ratio between ε = 10⁻⁷ and 10⁻⁹ (ε^{3/2}: 1000): median 1068, range 609–1252.

k = 3: 17 of 20 instances with an N-atom vertex, 15 with μ_min ≥ 0.01; μ_min median 0.0835, smallest 2.73e-4.

| ε | distance: median, largest | |w − w₁|∞: median, largest | the closed form: median, largest | projection: largest residual |
|---|---|---|---|---|
| 1e-05 | 1.17e-5, 0.0878 | 9.90e-10, 6.58 | 0.776, 138 | 1.38e-12 |
| 1e-07 | 1.43e-8, 1.50e-6 | 7.42e-14, 2.64e-11 | 0.0667, 1.38 | 7.04e-13 |
| 1e-09 | 1.46e-11, 6.84e-10 | 5.25e-13, 1.15e-12 | 0.00689, 0.0797 | 1.19e-12 |

The distance's ratio between ε = 10⁻⁷ and 10⁻⁹ (ε^{3/2}: 1000): median 984, range 185–2200.

The instances with μ_min < 0.01 (near-optimal second vertices; those whose distance at ε = 10⁻⁷ exceeds 10⁻³ followed to 10⁻¹¹):

| instance | μ_min | ε | distance | the leading order's distance | |w − w₁|∞ | support |
|---|---|---|---|---|---|---|
| k = 2, 1 | 0.00491 | 1e-05 | 0.00374 | 0.00424 | 4.98e-4 | [3, 6, 7, 9, 10, 11] |
| k = 2, 1 | 0.00491 | 1e-07 | 6.56e-6 | 6.56e-6 | 2.73e-10 | [3, 6, 7, 9, 10, 11] |
| k = 2, 1 | 0.00491 | 1e-09 | 5.67e-9 | 5.67e-9 | 5.39e-13 | [3, 6, 7, 9, 11] |
| k = 3, 1 | 0.00332 | 1e-05 | 0.025 | 0.16 | 0.135 | [0, 2, 3, 5, 8, 10, 11] |
| k = 3, 1 | 0.00332 | 1e-07 | 0.00331 | 0.00407 | 7.60e-4 | [2, 3, 5, 8, 10, 11] |
| k = 3, 1 | 0.00332 | 1e-09 | 1.03e-6 | 1.03e-6 | 6.20e-12 | [2, 5, 8, 10, 11] |
| k = 3, 1 | 0.00332 | 1e-11 | 8.62e-10 | 8.62e-10 | 1.28e-12 | [2, 5, 8, 10, 11] |
| k = 3, 13 | 2.73e-4 | 1e-05 | 0.17 | 292 | 292 | [0, 1, 2, 3, 5, 6, 8] |
| k = 3, 13 | 2.73e-4 | 1e-07 | 0.169 | 114 | 114 | [0, 1, 2, 3, 5, 6] |
| k = 3, 13 | 2.73e-4 | 1e-09 | 0.00427 | 0.00462 | 3.49e-4 | [0, 1, 2, 3, 5, 6] |
| k = 3, 13 | 2.73e-4 | 1e-11 | 3.20e-8 | 3.20e-8 | 4.37e-12 | [0, 2, 3, 5, 6] |

Numbers are produced by `experiments/exp_mirror_free.py`; nothing is edited by hand.
