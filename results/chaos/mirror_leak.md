# Proposition D‴.46: a leak in the mirror coordinates (`exp_mirror_leak.py`)

The free law with the leak −σ(η_i − η_p,i) on each particle's scalar drive (Q kept); η, B the mirror coordinate and Bregman divergence of Proposition D‴.45 (`mirror_free`).

## A. Identities

V̇ against −k_ee² + ed − n(f − f̂) − σΣ[B(w*, w) + B(w, w_p) − B(w*, w_p)] in the real loop (k_e = 5, a disturbance and a measurement noise in the adaptation, σ = 0 and 0.3, 12 random states): largest relative difference 9.54e-11.

The averaged equilibrium (k_e = 2) against its limits on random positive definite H (10 × 10, 4 draws), relative distances, the largest over the draws: the ridge with ρ = k_eσm/(Gβ⁰²) = 0.5 (its effect max|w_ridge − w*|/max|w*| given; at β⁰ = 3 the saturation's wall w_⊓ = β⁰a_max/m lies within the coefficients' range), the LASSO with penalty σλ₁(β⁰) = 0.02 (λ₁ exact; its leading term's relative error given):

| the limit | β⁰ | distance | the solver's residual | the ridge's effect, or the supports equal | w_⊓ against max|w*|, or λ₁'s leading term |
|---|---|---|---|---|---|
| ridge | 10000 | 1.31e-8 | 1.51e-17 | 0.307–0.554 | 833 against 0.295 |
| ridge | 100 | 1.31e-4 | 2.83e-17 | 0.307–0.554 | 8.33 against 0.295 |
| ridge | 3 | 0.147 | 8.61e-13 | 0.307–0.554 | 0.25 against 0.295 |
| LASSO | 0.0001 | 0.0152 | 4.54e-13 | 4 of 4 | 0.001 |
| LASSO | 1e-06 | 0.00151 | 1.57e-14 | 4 of 4 | 1.00e-5 |
| LASSO | 1e-08 | 1.52e-4 | 1.81e-13 | 4 of 4 | 1.00e-7 |

## B. (1) The real loop under an approximation error (0.1 sin 4x) and a disturbance (0.2 sin 3.1t + 0.1) (k_e = 5, m = 6, G = 1, γ = 0.5; six atoms of width 0.4; x_m's period 20π)

The orbit is persistently exciting: the filtered Gram H_L = k_e⟨φL[φ]ᵀ⟩ (L = 1/(s + k_e)) has eigenvalues 0.00269 to 0.678. Without the leak the averaged loop converges to the biased fit w* + H_L⁻¹d̄_L (d̄_L = k_e⟨φL[ε + d]⟩): max|w| = 0.337, ΣB(w*, w) = 0.465, its slowest time constant 1.52e+04. By quarters of each run: the largest |w|, the largest V (V = ½e² + ΣB(w*, w), w* the representable part), the rms error; the largest |Q − Q(0)|; for σ > 0 the late mean's distance to the averaged equilibrium w_σ^L, the comparison bound's largest ratio V/(V(0)e^{−λt} + (C/λ)(1 − e^{−λt})) (asserted ≤ 1) and the ultimate bound C/λ (d̄ = 0.4):

| σ | T | max |w| by quarters | max V by quarters | rms e by quarters | max |Q − Q(0)| | late mean − w_σ^L | V/bound (largest) | C/λ |
|---|---|---|---|---|---|---|---|---|
| 0 | 20000 | 0.309, 0.330, 0.337, 0.341 | 0.396, 0.43, 0.444, 0.45 | 0.0253, 0.0243, 0.0243, 0.0243 | 6.27e-11 | — | — | — |
| 0.01 | 3000 | 0.081, 0.081, 0.080, 0.081 | 0.182, 0.0847, 0.0847, 0.0847 | 0.0332, 0.0309, 0.0309, 0.0309 | 5.03e-12 | 0.0025 | 1.0000 | 1.74 |
| 0.1 | 3000 | 0.040, 0.040, 0.040, 0.040 | 0.182, 0.13, 0.13, 0.13 | 0.0382, 0.0366, 0.0366, 0.0366 | 1.62e-11 | 6.47e-4 | 1.0000 | 0.297 |

(2) A drift without persistent excitation: ė = −k_ee + (w* − w)φ(t) + d(t), φ(t) = (1 + t)^{−1/2}, w* = 0, one free particle (m = 1, G = 1, γ = 0.5, α(0) = 0.5, β(0) = 1), k_e = 1; d built so that without the leak η = η(0) + κ log(1 + t), κ = 0.5, and w(t) grows like (1 + t)^{1/4}; |d| ≤ 0.766. w at t = 1, 10², 10⁴, 10⁶ (and the construction's w(t) for σ = 0); for σ > 0 the comparison bound's largest ratio and C/λ:

| σ | w(1), w(10²), w(10⁴), w(10⁶) | V/bound (largest) | C/λ |
|---|---|---|---|
| 0 | 0.699, 2.31, 7.37, 23.3 | — | — |
| 0.01 | 0.698, 1.52, 0.515, 0.464 | 1.0000 | 29.5 |
| 0.1 | 0.687, 0.625, 0.467, 0.462 | 1.0000 | 3.07 |
| the construction (σ = 0) | 0.699, 2.31, 7.37, 23.3 | | |

## C. Slow adaptation (k_e = 20, G = 0.5, γ = 0.05, from α = 0 and equal masses β⁰, the leak toward the start)

The filtered Gram against the idealised one: max|H_L − H|/max|H| = 0.00888, asymmetry 3.22e-14. The real loop started at the averaged equilibrium w_σ^L (with H_L), its means over x_m's periods (20π): the last period's mean against w_σ^L and against the idealised w_σ (H/k_e), the first period's, the change between the last two periods' means, the last period's peak-to-peak oscillation; the lag |w_σ^L − w_σ| and the leak's bias |w_σ^L − w*|:

| start | β⁰ | σ | T (periods) | last mean − w_σ^L | last mean − w_σ | first mean − w_σ^L | the last periods' drift | oscillation | lag | bias |
|---|---|---|---|---|---|---|---|---|---|---|
| lazy | 3 | 0.0005 | 19981 (318) | 1.56e-4 | 6.94e-4 | 9.77e-5 | 1.62e-10 | 8.80e-4 | 8.17e-4 | 0.0382 |
| intermediate | 0.3 | 0.0001 | 59942 (954) | 1.74e-4 | 4.53e-4 | 1.32e-5 | 1.07e-9 | 3.64e-5 | 5.12e-4 | 0.123 |

## D. On an orbit's Gram (12 atoms of width 0.25 on [−1.4, 1.4]; eigenvalues 1.42e-4 to 0.553)

(1) A 3-sparse truth on atoms [2, 6, 9] (signs +, −, +; its irrepresentable quantity 0.771): the equilibrium against the ridge (lazy start, β⁰ = 30) and the LASSO of penalty μ = σλ₁(β⁰) (rich start, β⁰ = 10⁻⁶); its support (|w| > 10⁻³ max|w|), whether it keeps the truth's atoms, the atoms beyond them, max |w − w*| (absolute distances):

| start | μ (rich) or σ (lazy) | distance to the reference | support | the truth kept | extra atoms | max |w − w*| |
|---|---|---|---|---|---|---|
| lazy | 0.0001 | 3.86e-6 | [2, 6, 7, 8, 9, 10] | yes | 3 | 4.84e-4 |
| lazy | 0.0003 | 1.13e-5 | [2, 4, 5, 6, 7, 8, 9, 10, 11] | yes | 6 | 0.00143 |
| lazy | 0.001 | 3.50e-5 | [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11] | yes | 8 | 0.00458 |
| lazy | 0.003 | 8.70e-5 | [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] | yes | 9 | 0.0124 |
| lazy | 0.01 | 1.80e-4 | [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] | yes | 9 | 0.0312 |
| lazy | 0.03 | 2.49e-4 | [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11] | yes | 8 | 0.0571 |
| rich | 0.0001 | 3.71e-6 | [2, 6, 9] | yes | 0 | 7.52e-4 |
| rich | 0.0003 | 1.11e-5 | [2, 6, 9] | yes | 0 | 0.00226 |
| rich | 0.001 | 3.68e-5 | [2, 6, 9] | yes | 0 | 0.00752 |
| rich | 0.003 | 1.09e-4 | [2, 6, 9] | yes | 0 | 0.0226 |
| rich | 0.01 | 3.40e-4 | [2, 6, 9] | yes | 0 | 0.0752 |
| rich | 0.03 | 6.53e-4 | [2, 9] | no | 0 | 0.214 |

(2) Every 3-sparse support and sign pattern (1760, magnitudes 0.15–0.3), the rich start at μ = 0.001: the irrepresentable quantity ‖H_{S^cS}H_SS⁻¹s‖∞ is below 1 for 416; the LASSO's prediction (below 1 and no sign flip on S) holds for 416; the equilibrium's support is exactly the truth's for 415 — the prediction and the equilibrium agree on 1745 of 1760; where the support is not the truth's it carries 1.76 extra atoms on average (largest solver residual 9.98e-13).

(3) A biased disturbance correlated with the regressors (d̄ = Hv, v of size 0.05, a lazy start): the equilibrium's distance to w* and to w* + H⁻¹d̄ as σ → 0 — σ = 0.01: 0.0574, 0.0245; σ = 0.0001: 0.043, 3.70e-4; σ = 1e-06: 0.0428, 3.73e-6; σ = 1e-08: 0.0428, 3.73e-8.

(4) The averaged loop from the rich start (β⁰ = 10⁻⁶, μ = 10⁻³): the time at which max|w − w_σ| < 10⁻³ max|w_σ| — the truth above (irrepresentable quantity 0.771, the equilibrium's support [2, 6, 9]): t ≈ 3.2e+05; (0, 5, 11; +, −, +) (irrepresentable quantity 1.218, the equilibrium's support [0, 1, 5, 10, 11]): t ≈ 2e+06.

## E. The patterns of D(2) against the weighted LASSO itself (1760 patterns; the same Gram, magnitudes and rich start, β⁰ = 1e-06, μ = σλ₁ = 0.001, k_e = 1)

The LASSO min ½(w − w*)ᵀH(w − w*) + μ|w|₁ solved exactly (coordinate descent, then the active-set equations on its support, the KKT conditions checked); its exact support (w_j ≠ 0) and its support under D's threshold (|w| > 10⁻³ max|w|) against the truth's, the prediction of D(2) and the equilibrium's (threshold support, as in D):

| quantity | count |
|---|---|
| the LASSO's exact support is the truth's | 416 |
| — and agrees with the prediction | 1760 of 1760 |
| the LASSO's threshold support is the truth's | 437 |
| — and agrees with the prediction | 1739 of 1760 |
| the LASSO's threshold support equals its exact support | 1715 of 1760 |
| the equilibrium's support is the truth's (D(2)) | 415 |
| — and agrees with the prediction (D(2)) | 1745 of 1760 |
| the equilibrium's support equals the LASSO's exact support | 1717 of 1760 |
| the equilibrium's support equals the LASSO's threshold support | 1696 of 1760 |

The equilibrium against the LASSO, max|w_eq − w_LASSO| over the patterns: median 2.78e-4, largest 0.00197 (relative to max|w_LASSO|: median 0.00109, largest 0.00875); the LASSO's largest KKT residual 2.69e-14 (relative to μ), the equilibria's largest solver residual 9.98e-13.

Where the LASSO's threshold support differs from its exact support (45 patterns): the LASSO's smallest nonzero coefficient relative to its largest 5.59e-5–9.98e-4; irrepresentable quantity 1.0009–3.0503.

Where the equilibrium's support differs from the LASSO's exact support (43 patterns): atoms the equilibrium keeps beyond the LASSO's 31, LASSO atoms it drops 12; the LASSO's smallest relative slack 6.39e-4–0.355 (all patterns: median 0.128); the equilibrium's largest coefficient on the LASSO's inactive atoms, relative to its largest, 8.71e-9–0.00372, and the core's tangent profile there 8.71e-9–1.75; the LASSO's smallest nonzero coefficient, relative, 5.59e-5–0.831.

The equilibrium on the LASSO's inactive atoms against the core's tangent profile a_max β⁰^{3/2}/(m(2c)^{1/2})·cos(πμ_j/2)/sin³(πμ_j/2) (a_max β⁰^{3/2}/(m(2c)^{1/2}) = 3.73e-10) at the LASSO's relative slacks μ_j, the relative error per pattern: median of the medians 0.0259, largest 204.

Where the prediction and the equilibrium disagree (15 patterns), irrepresentable quantity 0.9930–1.0035:

- the prediction holds (8): the LASSO's exact support is the truth's for 8; the equilibrium keeps 8 atoms beyond it, its largest coefficient on the LASSO's inactive atoms relative to its largest 0.00108–0.00348 (above the threshold 10⁻³), at the LASSO's smallest relative slacks 0.00134–0.007;
- the prediction fails (7): the LASSO's exact support has 4 atoms, its smallest coefficient relative to its largest 5.59e-5–2.21e-4 (below the threshold), its threshold support is the truth's for 7, and the equilibrium's support equals the LASSO's threshold support for 7.

Numbers are produced by `experiments/exp_mirror_leak.py`; nothing is edited by hand.
