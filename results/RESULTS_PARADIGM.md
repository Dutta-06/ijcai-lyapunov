# Numerical results for PARADIGM.md (transport adaptive control)

Every number here is produced by the named script from raw files, and no table is edited by hand.

- `tests/test_transport.py`, `tests/test_langevin.py`, `tests/test_averaging.py`, `tests/test_excitation.py`,
  `tests/test_nonreciprocity.py`, `tests/test_periodic.py` and `tests/test_barycentric.py` check the inequalities at every
  sample of simulated trajectories, or on sampled states.
- Run them with `venv\Scripts\python -m pytest tests`. The current count is in `README.md`.

## Theorem B — label-free Lyapunov function (`exp_transport_permutation.py` → `transport_permutation/`)

Each network starts as a random relabelling of the ideal configuration plus noise 0.05, with mean-field gains and 3
seeds per m. "Excess" is |(1/m)Σ Taylor remainders|, the part that actually enters V̇.

| m | labelled error (t=0) | OT error W (t=0) | RMS ‖e‖ | actual excess | bound, matched error | same bound, labelled error |
|---|---|---|---|---|---|---|
| 8 | 2.675 | 0.107 | 0.011 | 2.0e-05 | 8.3e-02 | 2.9 |
| 16 | 2.532 | 0.104 | 0.007 | 9.6e-06 | 4.4e-02 | 1.9 |
| 32 | 2.684 | 0.112 | 0.008 | 2.7e-05 | 7.1e-02 | 2.1 |
| 64 | 2.764 | 0.114 | 0.009 | 1.2e-05 | 6.8e-02 | 2.4 |

- The bound is at least the actual excess at every sample of every run.
- The matched error is about 1/25 of the labelled one.
- The bound is 30–45× smaller when evaluated with the matched error.
- Nothing grows with m.
- The bounds remain far above the actual excess, because the Hessian supremum is global.

**Corollary B.5 along (L_σ)** (test `test_corollary_B5_fails_along_lsigma_and_its_leak_variant_holds`). The test
asserts every hypothesis of C(iv) and measures sup W ≤ η̄. There, the plain B.5 bound is exceeded by more than 3×,
while the leak variant with σ_m b²/4 holds.

## Theorem C — OT σ-modification (`exp_ot_sigma_mod.py`, `exp_sigma_mod_transient.py` → `ot_sigma_mod/`)

The two starts share plant, prior, geometry and x₀. They differ only in particle labels.

| start | mechanism | RMS ‖e‖ (20–40 s) | W_end | labelled error at end |
|---|---|---|---|---|
| aligned | projection | 0.0051 ± 0.0010 | 0.295 | 0.295 |
| aligned | labelled σ-mod | 0.0091 ± 0.0018 | 0.304 | 0.304 |
| aligned | OT σ-mod | 0.0091 ± 0.0018 | 0.304 | 0.304 |
| relabelled | projection | 0.0051 ± 0.0010 | 0.295 | 0.295 |
| relabelled | labelled σ-mod | 0.0091 ± 0.0018 | 0.304 | **2.059** |
| relabelled | OT σ-mod | 0.0091 ± 0.0018 | 0.304 | 0.304 |

Reshuffling transient from the relabelled start (5 seeds):

| leak | RMS ‖e‖ [0,5] s | [5,10] s | [10,20] s | [20,40] s | particle path length |
|---|---|---|---|---|---|
| labelled | 0.0771 | 0.0086 | 0.0093 | 0.0091 | 2.26 |
| OT | 0.0743 | 0.0072 | 0.0093 | 0.0091 | 0.39 |

- The labelled leak's longer path costs only about 4% of early RMS‖e‖. Theorem F does not price leak-driven transport
  (PARADIGM.md Proposition C.5).
- The C(iv) instance in the tests satisfies every hypothesis, asserted. V_OT stays below the envelope
  e^{−αt}V(0) + β/α at every sample and contracts by more than 5×.
- The C(ii) test asserts C ≤ b² + 4ηR, √C ≤ b + 2η, and C = b² under both gap conditions, on 300 configurations.

## Theorem D — convergence in m (`exp_mean_field_m.py` → `mean_field_m/`)

| m | sup_t ‖e^m − e^{2m}‖ | W(μ^m_0, μ^{2m}_0) | W(μ^m_T, μ^{2m}_T) | ratio |
|---|---|---|---|---|
| 16 vs 32 | 4.6e-3 | 0.218 | 0.366 | 1.68 |
| 32 vs 64 | 9e-4 | 0.109 | 0.179 | 1.64 |
| 64 vs 128 | 1e-4 | 0.055 | 0.092 | 1.69 |
| 128 vs 256 | <1e-4 | 0.027 | 0.046 | 1.68 |

This is a consistency check in m, not a check of the constant.

- The guaranteed factor e^{ΛT} is at least e^{496} on this box.
- The Van der Pol plant fails the tube condition of Theorem D(b), so the run is covered by the localised D(d).

## Theorem E — label-free convergence under excitation (`exp_pe_convergence.py` → `pe_convergence/`)

| reference | seed | λ_min, first 60 s window | min over sliding 60 s windows | W(T)/W(0) | fitted rate of W, t ≥ 100 s (1/s) | labelled(T)/labelled(0) |
|---|---|---|---|---|---|---|
| rich (6 tones) | 0 | 8.6e-3 | 4.1e-3 | 0.544 | 2.3e-3 | 0.992 |
| rich (6 tones) | 1 | 5.5e-3 | 1.0e-3 | 0.409 | 2.3e-3 | 0.962 |
| rich (6 tones) | 2 | 1.1e-3 | 5.9e-4 | 0.369 | 3.1e-3 | 1.007 |
| single sine | 0 | 2.5e-2 | 2.4e-2 | 0.274 | 4.1e-3 | 0.988 |
| single sine | 1 | 3.0e-2 | 2.7e-2 | 0.314 | 4.3e-3 | 0.967 |
| single sine | 2 | 6.7e-3 | 6.6e-3 | 0.226 | 2.9e-3 | 1.008 |

- The matched error decays in all six runs. The labelled error stays within 4% of its initial value.
- The rank correlation between the uniform excitation level and the fitted rate is 0.60 over the six runs. The PE
  level is a hypothesis of Theorem E, not a rate predictor (Remark E.1).
- The proof's linearisation is checked exactly: the Jacobian equals A(t), and QA + AᵀQ = diag(−k_eI, 0).

## Theorem F — transport cost (tests)

- Along (L), the particle path length (which is at least the Wasserstein length) stays below λ_G^{1/2}L₁∫‖e‖ at
  every sample.
- Along (L_σ) with the OT leak, W(μ_{t₁},μ_{t₀}) ≤ ∫(λ_G^{1/2}L₁‖e‖ + σ_mW(μ_r,μ⁰)) on a grid of pairs.
- Lemma P.1 holds on 2000 adversarial states, and a tight instance of F(i) reaches ratio ≥ 0.99.

## Proposition G — slow geometry (`exp_variable_projection.py`, `exp_slow_manifold_slope.py` → `variable_projection/`)

Ranges are along the run (5 samples). "a / b" means first-quarter / last-quarter medians. Velocities are centred
differences of period means.

| configuration | (H2a) ω_ref / joint fastest | (H2b) slowest weight / geometry | (H2c) S | G(a) direct check (cosine, G⁻¹ metric) | G(b) cosine | G(b) size ratio | L* | L̄ at period-mean weights |
|---|---|---|---|---|---|---|---|---|
| violating (m = 4) | 15–26 | 0.10–0.23 | 5.4–25 | 0.9994 / 0.9985 | −0.12 / 0.39 | 4.2 / 0.47 | **increases** (3.8e-3 → 5.9e-3) | non-increasing (0.219 → 0.008) |
| satisfying (m = 2) | 6.5–7.3 | 79–87 | 0.006 | 0.985 / 0.962 (geometry block 0.999 / 0.998) | 0.999 / 0.999 | 0.97 / 0.97 | decreases (1.20e-2 → 4.11e-3) | non-increasing |

Slope counterexample (`table_slope.md`). Here m = 1, k_e = 20, gains (g_a, δ, δ), and (H2b) holds in every row. The
closed loop and the joint flow start on the slow manifold of (b′).

| case | (H2b) ratio | S | closed / VarPro (b): size, cosine | closed / induced metric (b′): size, cosine | (a) / (b′) | closed / (a) |
|---|---|---|---|---|---|---|
| δ = 0.47, g_a = 20 | 81–178 | 2.0–2.8 | 0.42, 0.930 (centre moves the opposite way) | 1.12, 0.9986 | 1.015 | 1.10 |
| δ = 0.47/4, g_a = 5 (all gains ÷ 4) | 81–178 | 2.0–2.8 | 0.40, 0.916 | 1.065, 0.9999 | 1.015 | 1.05 |
| δ = 0.047, g_a = 20 | 654–1779 | 0.20–0.32 | 0.85, 0.9986 | 1.05, 1.0000 | 1.003 | 1.04 |
| δ = 0.0047, g_a = 20 | 13162–17790 | 0.020–0.023 | 1.02, 1.0000 | 1.04, 1.0000 | 1.000 | 1.04 |

- (a)/(b′) isolates the O(1/R) error of the induced-metric flow.
- closed/(a) isolates the averaging error, which shrinks with the gains, plus a 4.4% O(1/k_e) error-channel factor.

What the tables show:

- G(b) holds in the satisfying configuration, where (H2b) and S ≪ 1 both hold. The violating configuration fails both.
- The slope rows show that (H2b) alone is not sufficient.
- The induced-metric flow G(b′) tracks the closed loop at every S.
- In the satisfying configuration, the G(a) weight block is dominated by the O(1/k_e) shift of the weight equilibrium:
  cosine −0.21 / 0.50, carrying 92–95% of the residual. The geometry block matches.

## Theorem H — Langevin adaptive control (`exp_langevin_global.py`, `exp_langevin_tradeoff.py`)

The setting is the real closed loop with m = 32, a scalar plant and a periodic reference.

- The unregularised offline optimum of the occupation loss is about 7.6·10⁻⁶ (best of 5 L-BFGS restarts, an upper
  bound). One run gave 7.0e-6, because of floating-point summation order.
- The regularised optimum L̄(μ*_τ) is 0.11118 at (σ, τ) = (0.02, 0.004) and 0.00850 at (0.002, 5·10⁻⁴), grid-converged
  (`exp_gibbs_reference.py` → `langevin_global/gibbs_reference.md`).
- In the scaling σ = 0.2/k_e, τ = 0.04/k_e of Conjecture H.6, L̄(μ*_τ) = 0.11118 for every k_e. The exactly averaged
  stationary state gives L̄ = 0.11415, 0.11148 and 0.11101 at k_e = 10, 20 and 40. The L²(ρ) distance between the two
  learned functions is 6.8e-3, 1.4e-3 and 9.7e-4.

Speed–bias trade-off (σ = 0.002, 2 seeds, `langevin_tradeoff/table.md`):

| τ | escape of the collapsed start (s) | loss, collapsed | loss, spread |
|---|---|---|---|
| 0 | never | 0.166 | 0.0070 |
| 0.0005 | ∞, 1050 | 0.025 | 0.011 |
| 0.001 | 1480, 1005 | 0.022 | 0.015 |
| 0.002 | 1115, 800 | 0.029 | 0.023 |
| 0.004 | 970, 430 | 0.039 | 0.036 |
| 0.008 | 390, 295 | 0.060 | 0.065 |
| annealed (hot start) | 990, 530 | 0.032 | 0.032 |

- **Fit.** The seed-averaged escape time fits A + B·log(1/τ) with a 4.5% residual, against 25% for an Arrhenius fit.
  On the eight per-seed values both laws fit poorly (27% and 33%). The data do not identify the escape mechanism, and
  version 2's "saddle" reading is withdrawn (PARADIGM.md Remark H.2).
- **Caveats.** The threshold grows with τ, and the censored τ = 5·10⁻⁴ point is excluded from the fit.
- **Guaranteed rate.** The Holley–Stroock rate 2τρ_τ is at least e^{−8401} s⁻¹ at (0.002, 5·10⁻⁴) and at least
  e^{−1053} s⁻¹ at (0.02, 0.004).

Two fixed settings (3 seeds, `langevin_global/table_*.md`):

- **σ = 0.02, τ = 0.004.** Both starts equalise: the loss ratio is 1.03, and the final measures are 4.2× closer in W
  than without dither. The common loss is about 0.11, which matches the grid value L̄(μ*_τ) = 0.111. It is dominated by
  the regularisation bias.
- **σ = 0.002, τ = 5·10⁻⁴.** Dither frees the collapsed start (0.166 → 0.048), but the starts do not equalise within
  900 s (loss ratio 5.2).

## Remark H.0 and Proposition H.7 — the exactly averaged closed loop (tests)

Relative error of the frozen-μ averaged drift against (𝔄):

| instance | k_e = 10 | 20 | 40 | 80 | k_e × error |
|---|---|---|---|---|---|
| n = 1 (`test_langevin.py`) | 4.2% | 2.0% | 0.91% | 0.48% | ≈ 0.4 |
| n = 2, ellipse (`test_averaging.py`) | 43% | 21% | 10% | 5.1% | ≈ 4.2 |

Structure of the first correction (`test_averaging.py`):

| check | n = 1 | n = 2 (ellipse) |
|---|---|---|
| expansion error of k_e²(P_μ − P₀) against U_even + U_odd, k_e = 40 → 80 | 0.54 → 0.25 | 0.038 → 0.019 |
| circulation kernel K | ≡ 0 (asserted at 1e-12) | antisymmetric; equals the enclosed-area integral to 7e-15 (asserted at 1e-6) |
| orbit reversal (same ρ) | U_odd ≡ 0 | U_even unchanged, U_odd flipped, to 2e-15 (asserted at 1e-10) |
| k_e²·cyclic sum of the antisymmetric second variation, exact (k_e = 80) against predicted (asserted at 3%; the full matrix at 5%) | −0.0711 against −0.0704 | −0.0471 against −0.0464 |

## Remark H.7a — the orientation of the reference orbit (`exp_orbit_orientation.py` → `orbit_orientation/`)

The setup is deterministic: n = 2, m = 6, G = diag(4,4,1,1,1), σ = 2/k_e, no dither. The same elliptic reference
is traversed both ways, so the occupation measure is the same, and so is the initial configuration. The runs
converged to machine precision: the residual drift over the last 10 periods is about 1e-16. The six kernels converge
to a common configuration (kernel spread 1.05 at the start, at most 5.6e-17 at the end of every run), so the
learned network is effectively a single kernel and the test probes the effect at a collapsed equilibrium.
Predictions are linear responses at the (+) equilibrium.

| k_e | measured ‖Θ₊ − Θ₋‖ | frozen-measure prediction (H.7): size ratio, cosine | full prediction (H.7 + oscillation term): size ratio, cosine |
|---|---|---|---|
| 10 | 5.73e-3 | 2.67, 0.758 | 0.937, 0.931 |
| 20 | 2.84e-3 | 2.41, 0.830 | 0.965, 0.981 |
| 40 | 1.40e-3 | 2.29, 0.863 | 0.969, 0.994 |

- The vector error ‖measured − predicted‖/‖predicted‖ of the full prediction is 0.36, 0.20 and 0.11 at k_e = 10, 20
  and 40.
- The measured difference scales like k_e^{−1.01}, against a predicted −1.
- The full first-order theory accounts for it. The frozen-measure terms alone account for less than half.
- A first attempt with σ = 0.5/k_e and half the horizon had not converged (residual drift about a third of the
  difference). It was discarded and rerun.


## Lemma E.0 and Theorem E with explicit constants (`exp_explicit_rate.py` → `explicit_rate/`)

The explicit rate of Lemma E.0 against the exact Floquet rate of the linearised loop. It lies below the exact rate in
every row (a numerical ordering: Ψ_M, Ψ_D are grid maxima and ε a grid optimum; margins ≥ 17×). Full tables in
`results/explicit_rate/table.md`.

| instance | regime | exact / explicit |
|---|---|---|
| Remark E.1 (m = n = 1, G = I, k_e = 3), ω = 1 … 64 | Ψ_M²/k_e ≈ 1 | 954 → 17.6 (decreasing in ω) |
| the ideal of the Theorem E tests (m = 3, n = 2, k_e = 6), two-tone reference of period 20π, G = g·diag(2,2,1), g = 30 … 0.01 | fast → slow | 6.56·10¹⁰ → 4.22·10⁵ (optimised ε; 8.74·10¹⁰ → 5.63·10⁵ with the closed-form ε) |

- In Remark E.1's instance the exact rate falls like ω⁻²; the explicit rate approaches that law (log-log slope −1.56
  over ω ≥ 4, −1.92 over the last octave ω = 32 → 64). Lemma E.0(iii) derives the ω⁻² law from Ψ_D ∝ ω.
- In the slow-adaptation limit (TΨ_M² ≪ k_e) the ratio decreases towards a constant of order κ_Ψ(1 + Ψ_D/(k_eΨ_M))² for
  the optimised ε (3κ_Ψ(…)² for the closed-form bound (ii)), with κ_Ψ the condition number of the averaged Gramian
  (1.27·10⁵ for this ideal; formula values 4.9·10⁵ and 1.47·10⁶; measured 4.22·10⁵ with the optimised ε and 5.63·10⁵ with
  the closed-form ε at g = 0.01; on this instance the exact rate is 0.66 of λ_min⟨ΨᵀΨ⟩/k_e because k_e = 6 is not large
  against the reference's bandwidth, and the T²Ψ_M³ term is still 13% of c₁; the PE level is evaluated with 800
  Gauss–Legendre nodes, converged over the 20π window). The exact rate pays for ill-conditioning once, the bound twice
  (Remark E.0a).
- Lemma E.0′ (the slow regime; table 4 of the same file): on the g = 0.01 ideal (κ(M̄) = 7.8·10³; peak factors ρ = 69.9,
  ρ_D = 3.7·10³, ρ_DD = 3.7·10⁵; threshold k₀ = 226 for δ = 0.2, 321 for δ = 0.1) the guaranteed rate is within a factor
  1.26 (δ = 0.2) and 1.12 (δ = 0.1) of the exact Floquet rate at k = 400 and 600, where Lemma E.0's optimised rate is
  1.3·10⁵ times smaller; the certificate sup V̇/V ≤ −γ holds on the grid wherever the conditions hold; exact/(λ_min/k) =
  0.9885, 0.9971, 0.9993, 0.9997, 0.9998, 0.9999 at k = 50 … 600 (the k⁻² law of N₂ = 0).
- Corollary E.0′.2 (Theorem E in the slow regime; table 5): on the single-kernel instance at k_e = 40 the nonlinear rate from
  Lemma E.0′'s V is 4.53·10⁻⁴ (δ = 0.2) against 2.63·10⁻⁵ from Lemma E.0 and an exact linear rate of 1.35·10⁻³; the explicit
  ball radius is 1.53·10⁻⁶ against 6.15·10⁻⁶ (the anisotropy of the quadratic form), the overshoot 11.8 against 1.43.
- Tests (`tests/test_excitation.py`): the pointwise inequalities of Lemma E.0 on random states for four regressors, the
  envelope along trajectories, the ordering explicit ≤ exact, the ω⁻² law, the no-growth of the ratio as the gains
  shrink, and the explicit nonlinear envelope of Theorem E on the implemented closed loop (single-kernel instance with
  G = 0.3·I, where r₀ ≈ 10⁻⁵ in the w-coordinates and the explicit nonlinear rate is 3.7·10⁻⁵: the inequality holds at
  every sample over 8000 s, about 0.3 e-folds of the explicit envelope, while the state itself decays by far more).
  Lemma E.0′: the certificate and the sharpness on Ψ = (sin t, cos t) (k = 20, 50, 200; ratio → 1/(1 − δ)) and on the
  slow-adaptation ideal at k = 400 (ratio ≤ 1.5; more than 10⁴ above Lemma E.0's rate).
  Corollary E.0′.2: the explicit envelope on the implemented single-kernel loop at k_e = 40 over 1.5 e-folds (dt = 0.005,
  below the fixed-step discretisation floor).
- Corollary E.0′.3 (table 6): the explicit bound on the stationary parameter spread of the dithered (LTV), with Σ = (2·10⁻³)^{1/2}I,
  against the exact periodic stationary covariance: on Ψ = (sin t, cos t) the ratio bound/exact is 2.69, 1.58, 1.31 (δ = 0.2) and
  2.52, 1.41, 1.17 (δ = 0.1) at k = 20, 50, 200, falling toward 1/(1 − δ); on the Theorem E ideal (κ = 7.8·10³) it is 15.4–17.5,
  the price of q = 15 directions charged at the worst excitation; the exact value equals the averaged Langevin covariance k·tr X
  within 3·10⁻⁴ (ideal) and 1.25/k² (sinusoid; Σ isotropic, so ΣΣᵀ commutes with M̄ — otherwise the deviation is O(1/k)), and
  E‖e‖² = ‖Σ‖_F²/(2k) = tr(M̄X)/k to five digits (round-7 review).
- Lemma E.0″ (`exp_pe_general.py` → `explicit_rate/general_pe.md`): the exponential-average construction on two non-periodic
  regressors. Quasi-periodic Ψ = (sin t, cos √2 t): λ_T = 0.297, 0.391, 0.444 for windows 2π, 4π, 8π; conditions from k = 20;
  exact/guaranteed (δ = 0.2 / 0.1) 2.49/2.13 → 1.63/1.45 at T = 4π and 2.10/1.81 → 1.43/1.27 at T = 8π from k = 20 to 200; the
  exact rate is 0.5/k to within 1%. Chaotic Ψ = (x₁/20, (x₃ − 25)/20) along a Lorenz orbit traversed five times slower than standard
  time: λ_T = 0.015, 0.024, 0.051 for windows 5, 10, 20 against an effective excitation 0.15 (exact rates 0.00759/k·20 …);
  conditions from k = 100; exact/guaranteed 5.0/4.3 (k = 100) → 4.2/3.7 (k = 200) at T = 20, dominated by the window loss
  λ_∞/λ_T ≈ 3; the certificate holds at every admissible row.

## Theorem J — periodic orbits of the finite-m closed loop (`exp_periodic_orbit.py` → `periodic_orbit/`)

Newton shooting finds the exact periodic orbit near the critical point Θ* of J_m (residual < 10⁻¹¹); its period mean and
monodromy matrix are compared with the averaged system. Full tables in `results/periodic_orbit/table.md`.

| instance | k_e | ‖Θ̄_p − Θ*‖ | vector error of the first-order location (full / frozen only) | parameter Floquet exponents: max, median relative error |
|---|---|---|---|---|
| A (n = 2, m = 3, collapsed minimiser) | 10 / 20 / 40 / 80 / 160 | 4.09e-2 / 1.91e-2 / 9.22e-3 / 4.53e-3 / 2.24e-3 | 0.169 / 0.083 / 0.041 / 0.020 / 0.010 (frozen: 0.40 → 0.34) | 8.3e-2, 1.2e-2 → 1.3e-2, 1.1e-3 |
| B (n = 1, m = 2, two distinct kernels, fast sweep) | 10 / 20 / 40 / 80 / 160 | 6.44e-2 / 2.22e-2 / 8.37e-3 / 3.51e-3 / 1.59e-3 | 1.995 / 1.078 / 0.551 / 0.277 / 0.139 (b_osc = 0 by symmetry) | 9.5e-2, 4.0e-2 → 2.0e-3, 9.5e-4 |

- Fitted exponent of the offset against k_e: −1.05 (A), −1.33 (B; not yet asymptotic below k_e ≈ 80).
- Error channel: Liouville's trace identity along the orbit gives a mean error-channel exponent k_e − L with L = 0.16 → 0.07 (A) and
  L = −0.05 → −0.09 (B) for k_e = 10 → 160, as Theorem J(b) states (the finite-difference monodromy resolves these multipliers only as noise).
- The vector error of the first-order location halves with every doubling of k_e in both instances; in A the
  frozen-measure terms alone leave a 35% error at every k_e, so the oscillation term of Remark H.7a is an O(1) fraction.
- Tests (`tests/test_periodic.py`): ∇J_m against finite differences; shooting residual < 1e-10; the offset, the vector
  error and the Floquet error all fall by the required factors from k_e = 20 to 40 on instance B.

## Proposition K — collapsed configurations and the orientation of the orbit (`exp_orientation_stability.py` → `orientation_stability/`)

Instance A of Theorem J (n = 2, m = 3, elliptic reference traversed both ways). Floquet exponents of the collapsed
periodic orbit by Newton shooting, restricted to the separation subspace Σu_i = 0; the second-order prediction is the
linearisation of (1/k_e)⟨H⟩ + b_frozen + b_osc at its own equilibrium Θ̄₂. Full table in `results/orientation_stability/table.md`.

| | k_e = 10 | 20 | 40 | 80 |
|---|---|---|---|---|
| softest separation exponent, orientation (+) / (−), σ₀ = 0.5 | 3.563e-3 / 2.660e-3 | 1.800e-3 / 1.614e-3 | 9.277e-4 / 8.863e-4 | 4.737e-4 / 4.640e-4 |
| difference (+) − (−) measured / predicted at Θ̄₂ / at Θ* | 9.03e-4 / 7.77e-4 / 9.01e-4 | 1.86e-4 / 1.69e-4 / 2.08e-4 | 4.14e-5 / 3.93e-5 / 4.98e-5 | 9.71e-6 / 9.44e-6 / 1.22e-5 |
| k_e² × difference, measured / predicted at Θ̄₂ | 0.090 / 0.078 | 0.075 / 0.068 | 0.066 / 0.063 | 0.062 / 0.060 |

- Fitted exponent of the difference against k_e: −2.18 (Proposition K(iii): −2); the prediction at Θ̄₂ converges to the
  measured k_e²-scaled limit, the prediction at Θ* (third-derivative term omitted) converges to a value 30% too large.
- At σ_crit = 0.4606 both orientations are unstable (−4.1e-4 and −1.27e-3 at k_e = 10); the even part dominates.
- Designed σ₀ (pred₊ + pred₋ = 0): σ₀ − σ_crit = 0.0125 (k_e = 10), 0.0061 (k_e = 20); measured softest exponents
  +8.5e-4 / −2.5e-5 and +1.4e-4 / −3.7e-5: one orientation stable, the other unstable.
- Test (`tests/test_periodic.py::test_proposition_K_hessian_splitting_and_critical_leak`): the Hessian splitting at the
  collapsed point (1e-5), σ_crit's gap equation (1e-8), the two-fold zero separation mode at σ_crit.

## Theorem D′ — propagation of chaos for the dithered closed loop (`exp_chaos.py` → `chaos/`)

Synchronous coupling of the m-particle loop (n = 1, k_e = 10, σ = 0.02, T = 20) to an M = 4096 reference that shares its
first m particles and Brownian motions; sup_t of the coupled distance D_m = ‖x^m − x^M‖ + (1/m)Σ‖ϑ^m_i − ϑ^M_i‖, 24 seeds.
Full table in `results/chaos/table.md` (raw per-seed values in `raw_sup_D.npz`).

| | τ = 0 | τ = 2·10⁻³ |
|---|---|---|
| mean sup_t D_m, m = 16 → 1024 | 1.49e-1 → 2.17e-2 | 1.52e-1 → 2.32e-2 |
| slope of the mean against m (m ≤ 256), bootstrap s.e. | −0.47 ± 0.06 | −0.46 ± 0.05 |
| slope of the median | −0.46 | −0.45 |
| slope against (1/m − 1/M)^{1/2} (predicted 1) | 0.94 ± 0.08 | 0.92 ± 0.08 |

- The m^{−1/2} rate of Theorem D′(b) holds within one standard error at both temperatures; the shared-particle coupling
  predicts the mean ∝ (1/m − 1/M)^{1/2}, which the data follow up to m = M/4.
- Single seeds differ by more than an order of magnitude at fixed m; a three-seed version of the experiment gave slopes
  −0.25 and −0.21 and was withdrawn (review round 5, D2). The constant of the theorem (T K_Φ e^{ΛT}) is not checked.
- No test: the experiment is the check.
- Theorem D″ (`exp_chaos_uniform.py` → `chaos/uniform_table.md`; n = 1, k_e = 25, G = 0.05·diag(4,1,1), a two-bump plant that
  is not odd, M = 2048): exact Lipschitz constants L_x = L_ϑ = 8.087, L₂ = L_{xϑ} = 44.4 and the a-priori tube E_X = 0.228 give
  σ_cert = 4.04, against Proposition K's critical leak σ_K = 0.00144 of the collapsed critical point (ratio 2800). At 1.2σ_cert
  with dither τ = 0.02 (80 s, 8 seeds) the coupled distance is 6.7·10⁻⁴ → 7.9·10⁻⁴ (m = 16) and 1.6·10⁻⁴ → 1.9·10⁻⁴ (m = 256)
  between the first and the last window, slope −0.51 ± 0.02, bound 0.55 and 0.14; the network stays at its prior (fit error
  1.000, rms spread 0.033 = the dither cloud). Without dither the reference collapses (spread 10⁻¹⁴) and the distance is zero;
  at σ_cert/8 it decays to 10⁻¹¹. In the learning regime (σ_K/3, 1000 s, 6 seeds) the network stays spread (rms 1.53) and fits
  the plant (fit error 0.632), and the distance follows m^{−1/2} (slope −0.62 ± 0.16) but grows threefold between the windows
  (100, 200] and (500, 1000]: no empirical uniformity there (the first draft's "spread regime" run was a slow collapse; review
  round 6).

## Theorem D‴ — uniform-in-time propagation of chaos in the learning regime (`exp_chaos_learning.py` → `chaos/learning_table.md`)

- Instance: n = 1, G = diag(1, 0.25, 0.25), three-bump plant (1, −1.5, 0.5), (−0.7, 0, 0.4), (0.8, 1.5, 0.5), x_m = 2 sin t, σ_crit = 0.04716
  (Proposition K), σ₀ = 0.15σ_crit = 0.00707; L_ϑ^G = 4.044, L₂^G = 11.11, K_Φ = 4, F_ρ = 0.999. Averaged slow loop (𝔄_m), m = 16, 64, 256
  against M = 2048, 900 slow units, 4 seeds per start; real loop at k_e = 25 over 2000 time units, 3 seeds.
- The hypothesis at the three-atom stationary measure (fit error 0.244): κ_R at the atoms 0.00004, −0.00227, 0.00411; over spheres of
  radius 0.05, 0.1, 0.2, 0.3: 0.00461, 0.00512, 0.00619, 0.00730 — (H) holds with ϱ′ = 0.1 and margin 0.00195 (1/κ = 512), the window
  κ_R[μ*](N_{0.1}) = 0.00512 < σ₀ = 0.00707 < σ_crit = 0.04716 is realised; at 0.3σ_crit the three-atom measure fails splitting
  (κ_R = 0.01696 > 0.01415). Theorem constants: C₁ = 2.2·10⁹, m₀ = 1.2·10²³ (vacuous envelopes).
- Along the warm run (jittered mixture around the atoms) the margin a(s) stays in [+0.00139, +0.00167]; the own-residual statistic
  exceeds the r*-one by ≤ 0.00013 ((0.2) allows 0.214); rms drift 6·10⁻³ → 3·10⁻⁶. Along the cold run (connected cloud) the margin is
  negative throughout (−0.196 → −0.0025) while the fit converges to 0.244 as well.
- Warm start: window maxima of D_m 0.396 → 0.448 (m = 16), 0.122 → 0.138 (64), 0.0119 → 0.0137 (256) over (0, 112.5] → (450, 900]
  (last/first 1.13–1.15, last/previous 1.01–1.02: saturated); slope −1.28 ± 0.30 (mean square; theorem −1), learned function −0.69 ± 0.16;
  per-particle maxima 1.36, 0.62, 0.21 ≫ 3ϱ′/4 = 0.075 (Remark D‴.6).
- Proposition D‴.5: ‖S_k‖ = 8.38, 4.66, 8.32; E D_∞ predicted 0.823, 0.201, 0.0454 against means 0.446 ± 0.28, 0.138 ± 0.051,
  0.0136 ± 0.0055 (four seeds, heavy-tailed); per seed the linearised prediction is within 25% of the measured D in 9 of 12 pairs
  (worst ratio 0.57 at |w^m − w^M| = 0.275) and the exact stationary distance of the replicated weighted configuration within 2.3% in
  all 12; splitting margins positive in 11 of 12 (+0.0006 to +0.0050), negative (−0.0022) for m = 16, seed 3 (weights (9,4,3)/16: a
  saddle, left later by particle migration — review round 9); per-cluster displacements along S_k(w^m − w^M) with cosines 0.92–1.000
  and norms within 0.8–20%.
- Cold start: D grows 1.81 → 5.33 (16), 0.505 → 1.79 (64), 0.0744 → 0.583 (256) (last/first 2.95, 3.54, 7.84), slope −0.87 ± 0.14;
  bimodal deviations (medians 0.40, 0.14, 0.067; maxima 6.2–6.7); fraction of pairs > 1 whitened unit apart 0.203 ± 0.047,
  0.082 ± 0.021, 0.037 ± 0.010, decaying like m^{−0.64±0.08} (medians m^{−0.62±0.07}): the floor of Remark D‴.2.
- Real loop (k_e = 25, warm start): window maxima 0.328 → 0.590 (16), 0.0605 → 0.0909 (64), 0.0188 → 0.0294 (256) over (0, 500] →
  (1000, 2000] (slow time ≤ 80 = 0.16/κ), slope −1.13 ± 0.15; on the same initial clouds the averaged flow gives D_real/D_averaged
  between 0.993 and 1.004 at every window end (across different clouds the last-window ratios are 1.57, 0.78, 2.58): the real loop
  follows its averaged slow dynamics to within 0.7% (Corollary D‴.7); two of the three m = 16 clouds have inadmissible or marginal
  weights (margins −0.0021, −0.00001).

## Proposition D‴.9 — measure-valued priors and the design of the window (`exp_prior_design.py` → `chaos/prior_design.md`)

- Same instance and leak as Theorem D‴'s numerics (σ₀ = 0.15σ_crit = 0.00707), nearest-atom leak toward three prior atoms. Matching prior
  (exact representation): fit error 0, κ_R = 0, margin σ₀ = 0.00707 against the σ-modification's 0.00196 (the theorem's relaxation time 1/κ 141 against 511 slow units; observed late-time e-folding of the coupled distance 48–84 against 107–306); ‖S_k‖ = 13.01, 6.18, 8.91 against 8.38, 4.66, 8.32; E D_∞(m = 64) 0.323 against 0.201 (1.61×); function-level floor 9.33·10⁻⁴
  against 4.29·10⁻⁴ (2.18×).
- Misplaced priors: centres +0.15/+0.3/+0.45/+0.6 → margins +0.0055/+0.0037/+0.0017/−0.0006; amplitudes ×0.5/×0.7/×1.2 → +0.0074/+0.0073/−0.0317;
  widths ×0.7/×1.5/×2 → −0.0261/+0.0058/+0.0040; the shrunk priors (×0.5, ×0.7) have the largest margins and lower floors (E D_∞ 0.199, 0.230;
  function level 7.1·10⁻⁴, 7.7·10⁻⁴) at fit errors 0.09, 0.06.
- Admissible weights (splitting margin of Θ*(w)): (7,5,4)/16: +0.0026 / −0.0054 / +0.0073; (9,4,3)/16: −0.0022 / −0.0058 / +0.0057;
  (10,3,3)/16: −0.0038 / −0.0042 / +0.0054; (0.5,0.3,0.2): −0.0006 / −0.0065 / +0.0066 (σ-modification / matching / amplitudes ×0.7): the
  shrunk prior keeps all four admissible (above σ₀ only at (7,5,4)/16; its over-weighted atom overshoots by 1–3% at the other three), and
  the matching prior's set is smallest toward the largest bump (an excess 0.04 on atom 1 alone is inadmissible, 0.08 on atom 2 or 3 alone is not).
- Matching-prior coupled runs (m = 16, 64, 256 vs M = 2048, 300 slow units, 3 seeds): D saturates (last/previous 1.04–1.05), slope −1.23 ± 0.40;
  per-seed linearised-law ratios 0.81–1.11 (parameter level) and 0.61–1.69 (function level), exact law 0.987–0.998 in all nine pairs; two of
  nine realised weight vectors are saddles (margins −0.0005, −0.0030), left by a local split inside the cell with the coupled distance on its
  plateau (round 10); observed e-folding times 48–84 slow units (median 59) against 107–306 (median 195) for the σ-modification.
- Function-level law for the σ-modification against Theorem D‴'s warm runs: predicted RMS 0.0419, 0.0207, 0.0098 vs measured means
  0.0307 ± 0.0097, 0.0154 ± 0.0033, 0.0047 ± 0.0012 (ratios 0.73, 0.74, 0.48; RMS over seeds/predicted 0.84, 0.79, 0.52; per-seed law from the
  realised weights 0.60–1.57).

## Proposition D‴.12 — the OT leak's mass constraint against the weight fluctuation (`exp_ot_leak.py` → `chaos/ot_leak.md`)

- Same instance and matching prior as Proposition D‴.9's numerics (σ₀ = 0.15σ_crit, q = (1/3, 1/3, 1/3), Δ_ν = 6.26); the nearest-atom
  and the barycentric OT leak on the same warm clouds, 600 slow units (4.2/σ₀), 4 seeds; m = 15, 60, 240 vs M = 510 (commensurate) and
  m = 16, 64, 256 vs M = 512 (incommensurate).
- Cluster weights: conserved exactly under the nearest-atom leak (max deviation 0 in 32 systems); under the OT leak the mass in transit
  ½|w(0) − q|₁ = 0.117, 0.058, 0.027, 0.024 (seed means; 1–21 particles per system) reaches zero and the weights equal q from s_eq = 93–164
  slow units (ln 2/σ₀ = 98). Fit error at the end (seed means): 0.067, 0.029, 0.017, 0.012 (nearest) → 2.8·10⁻⁴, 1.2·10⁻⁴, 1.2·10⁻⁴,
  0.8·10⁻⁴ (OT; the migrated particles still 0.06–0.16 from their atoms, the OT finals not yet stationary).
- Commensurate sizes, last window: function-level distance² 2.5·10⁻³, 5.2·10⁻⁴, 5.9·10⁻⁵ (nearest; law predicts 4.0·10⁻³, 9.1·10⁻⁴,
  1.4·10⁻⁴) vs 4.9·10⁻⁷, 4.9·10⁻⁸, 7.5·10⁻⁹ (OT, still falling ×0.02–0.035 per 150 slow units); W² 5.1, 2.7, 0.81 (nearest; cross-cluster
  transport ½|w^m − w^M|₁Δ_ν² = 4.9, 2.5, 0.70; slope −0.67 ± 0.13) vs 7.9·10⁻³, 2.3·10⁻³, 6.9·10⁻⁴ (OT, ×0.12 per 150 slow units);
  synchronous D 0.87, 0.25, 0.017 (nearest) vs 5.8, 3.3, 0.94 (OT: the reassigned particles; at 0.95–0.96 of the mismatch plateau and
  still rising). Exact K × K transport between the weighted stationary measures at the realised weights: 0.997–0.999 of the measured W² at
  the seven admissible pairs.
- Incommensurate sizes: one or two split (bridge) particles per system at 0.28–0.58 Δ_ν from the nearest atom; against ν (the mean-field
  stationary measure): W²(μ_m, ν) = 1.21, 0.33, 0.079, 0.038, m·W² = 19.4–21.0 (slope −0.98 ± 0.03), fit error 0.028, 0.0093, 0.0023, 0.0011,
  m·fit = 0.45–0.60 (squared slope −1.81 ± 0.05); against the M = 512 reference (two bridge particles of its own) OT 3.4·10⁻⁴, 3.1·10⁻⁵,
  1.1·10⁻⁶ (function level; m = 16 still relaxing, plateau 2.8·10⁻⁴ in round 11's continued run) and 1.16, 0.29, 0.041 (W²) vs nearest
  2.5·10⁻³, 5.0·10⁻⁴, 4.4·10⁻⁵ and 4.8, 2.7, 0.69; nearest fit errors 0.061, 0.035, 0.013, 0.012.
- Masses as a design input: q_δ = (1/3 + δ, 1/3 − δ/2, 1/3 − δ/2) gives fit error 0.0069–0.049 for δ = 0.01–0.1, splitting margin +0.0048,
  +0.0028, +0.0011, −0.0015, −0.0049 (inadmissible between δ = 0.03 and 0.05), diagonal plan unique at every row (ε ≤ 1.04 < 3.13).

Tests (`tests/test_ot_leak.py`): the plan is a vertex with the prescribed masses (integral when commensurate, ≤ K₀ − 1 split rows
otherwise) and the single-atom prior is the σ-modification; the leak is the G-gradient flow of J_m^{OT} against finite differences; the
replicated weighted configuration is stationary for three commensurate sizes with the diagonal plan and the same learned function (W² = 0
between sizes), for matching and mismatched masses; conservation against constraint on a short run; the forest bound; the exact W² by
replication against the transportation LP.

## Corollary D‴.14 — the real loop with the OT leak (`exp_ot_real_loop.py` → `chaos/ot_real_loop.md`)

- Real loop at k_e = 25 (σ = σ₀/k_e, RK4 plant with the control and the reference evaluated at the stages, Euler parameters, dt = 10⁻³,
  100 time units = 16 periods) from the clouds at the end of Proposition D‴.12's averaged runs (contracted where the weights are admissible) (seeds 0, 1; m = 15, 60, M = 510) and from
  the collapsed lifts; last five periods on the reference orbit.
- OT leak: ‖f̂_m − f̂_M‖²_{L²(orbit)} = 7.5·10⁻⁹–2.7·10⁻⁸, against the 3-atom lift 5.2·10⁻¹⁰–3.6·10⁻⁸; fit error along the orbit
  3.5·10⁻⁵–2.9·10⁻⁴; rms tracking error 0.9·10⁻⁶–7.3·10⁻⁶ (the lift 8·10⁻¹²); W²(μ_m, μ_M) 3.3·10⁻⁴–1.4·10⁻³; within-class spreads 0.013–0.033.
- Nearest-atom leak: ‖f̂_m − f̂_M‖² = 1.1·10⁻⁴–2.5·10⁻³; fit errors 0.012–0.075; rms tracking errors 3.2·10⁻⁴–1.9·10⁻³; W²(μ_m, μ_M) 1.6–5.5.
  Three of six systems have inadmissible realised weights (margins −0.0035, −0.0011, −0.0022): their lifts are saddle orbits and the systems
  have split (spread 0.04–0.11; lift-to-lift ratios 0.960, 0.998, 0.961); at the one pair with both weights admissible the disagreement is
  the lifted orbits' distance to 1.0003, the admissible systems within 2.5·10⁻¹¹ of their lifts (spread 0.0006–0.0013).

Tests (`tests/test_ot_leak.py`): the collapsed 15-system and the 3-atom loop integrate identically under the OT leak at masses q (10⁻¹²
over 300 steps); the collapsed system with weights (6, 5, 4)/15 stays collapsed under the nearest-atom leak and loses a particle under the
OT leak; two commensurate lifts keep W² = 0 and the same learned function; the lifted orbit's class separation subspaces are invariant and
carry the ghost spectrum (Proposition D‴.15, m = 6).

## Corollary D‴.14 in two dimensions (`exp_ot_2d.py` → `chaos/ot_2d.md`)

- Planar plant (four Gaussian bumps with vector amplitudes), Lissajous reference, p = 5 parameters per kernel, matching prior (masses ¼),
  k_e = 25, σ₀ = 0.05, 2500 time units = 100 slow units, sizes 8, 32 vs 128 from the same jittered clouds, two seeds.
- OT leak: weights equalised at 13–14 slow units (ln 2/σ₀ = 14); ‖f̂_m − f̂_M‖² = 5.0·10⁻⁸–1.6·10⁻⁷; fit error along the orbit
  2.1·10⁻⁴–8.7·10⁻⁴; rms tracking error 3.9·10⁻⁶–1.6·10⁻⁵; W²(μ_m, μ_M) 2.0·10⁻⁴–4.4·10⁻⁴.
- Nearest-atom leak: weights conserved; ‖f̂_m − f̂_M‖² = 5.5·10⁻³–2.3·10⁻²; fit errors 0.11–0.39; rms tracking errors 1.9·10⁻³–6.7·10⁻³;
  W²(μ_m, μ_M) 4.7–11.8.

## Proposition D‴.15 — the transverse spectrum of a lifted orbit (`exp_lifted_orbit.py` → `chaos/lifted_orbit.md`)

- Real loop at k_e = 25; monodromies by central differences (DOP853, rtol 10⁻¹¹). OT leak at the matching masses: the 3-atom orbit (the exact
  representation, zero tracking error) lifted to m = 6, 9, 15 — return residuals ≤ 1.7·10⁻¹⁵; separation exponents σ₀/k_e = 2.830·10⁻⁴ in
  every class (= ghosts; invariance defects ≤ 1.8·10⁻⁸); collective exponents 2.935·10⁻⁴ … 7.934·10⁻³ identical across m to 10⁻⁴ relative,
  Theorem J(b)'s prediction 2.952·10⁻⁴ … 8.091·10⁻³ (0.6–3.2%).
- σ-modification at (7, 5, 4)/16 (margin +0.00257) and (9, 4, 3)/16 (−0.00219), lifts to m = 16 (residuals ≤ 4.9·10⁻¹⁴): separation exponents =
  ghosts to 2·10⁻⁵ relative, = (σ₀ − λ_j(R^G_k))/k_e within 0.2–4.1% (defects ≤ 5.7·10⁻⁸); the saddle's class 0 exponent −9.13·10⁻⁵ (predicted
  −8.77·10⁻⁵, multiplicity 8); period-start offsets from Θ*(w) 1.1·10⁻², 1.6·10⁻² (Theorem J(a)).

## Theorem D‴.18 — Theorem J uniformly in m (`exp_uniform_basin.py` → `chaos/uniform_basin.md`)

- Real loop at k_e = 25 (σ₀ = 0.15σ_crit, dt = 5·10⁻³), 8000 time units, many systems integrated together.
- The test of the rates (Part F, by shooting): the 3-atom orbit's parameter Floquet exponents within 0.56–3.22 % (matching prior) and
  0.57–3.15 % (prior centres +0.15) of λ_j(H_q)/k_e; the ghosts' exponents σ₀/k_e = 2.830·10⁻⁴ (matching) and (σ₀ − λ_j(R^G_k))/k_e to
  3·10⁻⁴ relative (shifted; slowest 2.2701·10⁻⁴ against 2.2696·10⁻⁴). The proof's radius is at most 0.019 (sampled C_N ≥ 0.092): the runs
  start 31–261 times outside the proved start radius.
- Matching prior at the matching masses (constant lift), sizes 3, 15, 60, 240, per-particle random displacements ρ: inside the cells
  (ρ = 0.3, 1.0, 2.5) no particle leaves its cell and both leaks are identical; fitted decay rates of the per-particle distance
  2.69–2.97·10⁻⁴ for m ≥ 15 (σ₀/k_e = 2.83·10⁻⁴; ratios 0.95–1.05), constants C(8000) = 1.03–1.29 (still growing at the horizon for
  m = 15 at ρ = 0.3, 1); at m = 3 (no separation modes) the runs are mixtures of collective modes (2.948·10⁻⁴ when the slowest carries the
  end, λ_coll/k_e = 2.952·10⁻⁴; 4.06·10⁻⁴ when the third does).
- Across a cell boundary (ρ = 4.5; 1, 3, 13, 55 particles start in another cell): the nearest-atom leak ends at counts (2, 0, 1), (6, 5, 4),
  (23, 19, 18), (87, 70, 83) — at m = 60, 240 through crossings during the transient — with W² 12.3, 2.81, 2.30, 1.92, already that of the
  lifts of the realised weights (3.21, 2.35, 1.93 for m = 15, 60, 240) though the clouds have not collapsed; (6, 5, 4)/15 and (23, 19, 18)/60
  are inadmissible (margins −0.00348, −0.00161). The OT leak ends at q for every m (last count change at t ≤ 2695), W² 0.010–0.23 falling at
  the leak's rate 2σ₀ per slow unit (the transit is over); the plan relabels 4, 15, 53 particles out of their initial class's cell.
- One particle across a boundary (Remark D‴.19; mean-square displacement 0.311, 0.078 at m = 60, 240): W² grows to 0.79, 0.20 under the
  nearest-atom leak (the lifts of (19, 21, 20)/60 and (79, 81, 80)/240: 0.8655, 0.2163) and falls to 3.0·10⁻³, 8.2·10⁻⁴ under the OT leak.
- Moving orbit (OT leak, prior centres +0.15): second-half fits 2.18–2.63·10⁻⁴ against 2.27·10⁻⁴ (ratios 0.96–1.16; mixtures of separation
  modes, consistent with the bound (1 − κ)μ), constants 1.04–1.23.

Tests (`tests/test_uniform_basin.py`, at the prior with centres shifted by +0.15, where the residual is not zero and λ_sep = 0.80σ₀): the
drift's second-order remainder ≤ C_N|||u|||² with C_N < 0.1 for random per-particle displacements at m = 6, 30, 150, 600, not growing with m;
the displaced particle's remainder exactly N_a + N_b/m with N_a the same for every m (non-zero here, zero at the matching prior), the others'
exactly O(1/m); a separation perturbation moves only its own particles at first order, with the response (GR_k − σ₀I)v at m = 6, 60, 300;
the averaged loop contracts |||·||| at μ = min{λ_coll, λ_sep} with one constant C ≤ 1.05 for m = 6, 30, 150.

## Proposition D‴.20 — the entropic OT leak (`exp_entropic_leak.py` → `chaos/entropic_leak.md`)

- Matching prior, σ₀ = 0.15σ_crit; the entropic plan by semi-dual Newton; every weighted stationary point continued from a converged
  neighbour and checked to 10⁻⁹.
- Critical excess at w_δ = q + δ(1.5e₁ − 0.5), by continuation in ε (residuals ≤ 8.7·10⁻¹³): ε_c = 5.11, 13.5, 23.1, 65.0 for δ = 0.02, 0.05,
  0.1, 0.2; at ε_c the outflow p = 0.057, 0.15, 0.33, 0.57 splits about equally between the two under-weighted atoms (Cov^G₁ of rank two);
  the exact threshold 2σ₀λ_max(A^{−1/2}Cov^G₁A^{−1/2}) equals ε_c; the rank-one formula 1.22–1.31 (atom 2) and 0.70–1.22 (atom 3) times ε_c.
  Excess from atom 3 alone: ε_c = 6.99, 14.6, 24.3 for δ = 0.02, 0.05, 0.1. At fixed ε: δ_c = 0.0042, 0.0083, 0.016, 0.030, 0.062, 0.14, 0.20
  at ε = 1, 2, 4, 8, 16, 32, 64, equal to w₁p_c up to ε = 8.
- Dynamics (m = 15, 60, 240 vs M = 510, 600 slow units, 4 seeds; ε = 4, 16, 64 continued to 6000): ε = 1 — all equalised, fit errors
  1.7–3.0·10⁻⁴; ε = 4 — five of the twelve larger systems at unstable weights at s = 600, every system at stable weights at s = 6000 with
  residual excesses up to 1/60, 3/240, 7/510 (δ_c(4) = 0.012–0.016); ε = 16 — the stable systems keep their sampled weights, three unstable
  ones move ((6, 5, 4)/15 → (5, 6, 4): the excess moved to an atom at q, not reduced); ε = 64 — only (8, 4, 3)/15 moves, fit errors
  0.12–0.14 (the barycentric bias).
- Floors: Proposition D‴.5's law at Θ*_ε(q) is 0.55 (ε = 16) and 0.74 (ε = 64) of the nearest-atom law; observed 0.23–0.84 and 0.47–0.91 of
  the nearest-atom leak's floors on the same clouds, per seed.

Tests (`tests/test_entropic_leak.py`): marginals and the agreement of the two solvers; the gradient of J^ε_m against finite differences; the
limits ε → 0 and ε → ∞ at their rates (the OT plan's uniqueness checked); the separation blocks of the finite-difference Hessian at a
collapsed configuration equal to I ⊗ (σ₀(I − (2/ε)Cov^G) − R^G) entrywise, the separation subspace invariant; the over-weighted class of
(6, 5, 4)/15, continued in ε with A ≻ 0 checked and every residual below 10⁻⁹, unstable at ε = 16, 17 and stable at 18, 20, 64.


## Theorem D‴.21 — a measure-level basin for the OT leak (`exp_ot_basin.py` → `chaos/ot_basin.md`)

- Averaged slow loop at σ₀ = 0.15σ_crit, the warm clouds of Proposition D‴.12, m = 15, 60, 240, 510, four seeds, 600 slow units, the matching
  prior and the prior with centres +0.15, both leaks. The identities (E), (P) hold along every run to 6·10⁻¹³ relative.
- OT leak: no particle relabelled; the margin α positive from s ≤ 300 on (from s ≤ 25 in 19 of 32 runs); the observed ratio at least λ at
  every record with α > 0 (smallest ratio/λ 1.00 at the matching prior, 1.61–1.72 at the shifted one); E from 4.1·10⁻³–4.0·10⁻² to
  8.3·10⁻⁷–8.1·10⁻⁶ (W²(μ_Θ, μ*) from 0.99–5.25 to 2.3·10⁻⁴–2.3·10⁻³), below the a-posteriori bound; at the end the observed rate is
  0.01416–0.01417 against 2σ₀ = 0.01415 and λ = 0.0140–0.0141 (matching prior: sharp); at the shifted prior λ = 0.0077–0.0078 against the
  observed 0.0137–0.0139, the local PL constant 2λ_sep = 0.01135 (curvature cost 0.802, the slack of (ii) 0.683).
- Nearest-atom leak: no particle changes cell. At the 22 admissible count vectors E(600) is 0.02–0.24 % below the energy of the replicated
  configuration of the sampled weights and rising toward it — a critical point of the nearest-atom energy (drift ≤ 6·10⁻¹⁴) but not of the OT
  energy (drift 2.3–4.4·10⁻²); at the 10 inadmissible ones (m ≤ 60; class-0 margins −0.0004 to −0.0061) class 0 splits and E ends
  0.5–16.5 % below it, still falling.
- The a-priori basin: E₀ = 2.0·10⁻⁷ at the matching prior (κ = σ₀); at the shifted prior (H*) holds on the whole parameter space
  (sup λ_max R^G[μ*] = 6.60·10⁻³ < σ₀, at a saturation corner; κ = 4.7·10⁻⁴, E₀ = 9.1·10⁻¹⁰). {E ≤ E₀} holds a fraction ≤ 5.8·10⁻⁶ (matching)
  and 4.4·10⁻⁷ (shifted) of the particles in other cells.
- The measure-level basin (matching prior): a fraction p of the particles moved 0.6 of the way into another cell stays in its class; κ^Θ =
  0.0096, 0.019, 0.029, 0.048 and α = +0.0023 at p = 1/30, negative from p = 1/15 (the same at m = 30 and m = 300: the same measures). One
  particle moved: E 3.3·10⁻³ → 2.5·10⁻⁵ (m = 30) and 2.3·10⁻⁴ → 3.1·10⁻⁶ (m = 300), at 2σ₀ at the end, below E(0)exp(−∫λ₊), the particle back
  in its cell; under the nearest-atom leak E rises to the spurious minima's 6.3126·10⁻³ and 6.1328·10⁻⁴ by s = 3000.

Tests (`tests/test_ot_basin.py`): (E) and (P) to 10⁻¹² for random assignments with the counts mq at both priors, and their exact defects at
other counts; the PL chain; the constant 2σ₀ attained by a separation to 10⁻⁴; the decay from a particle in another cell, with a positive
margin, below the a-posteriori bound, at 2σ₀ at the end, against the nearest-atom floor; the spurious minimum (11, 9, 10)/30 of the
nearest-atom energy, admissible (separation margins and K-atom Hessian positive), not critical for the OT leak.

## Theorem D‴.23 — adaptive priors (`exp_adaptive_prior.py` → `chaos/adaptive_prior.md`)

- Averaged slow loop at σ₀ = 0.15σ_crit, the OT leak toward the prior with centres +0.15 (masses 1/3); the prior's atoms move toward the
  class barycentres at the rate η.
- Linear theory at the exact representation: h_j = 3.06·10⁻⁴ … 0.195; a fixed prior keeps σ₀/(h_j + σ₀) = 0.959 … 0.035 of its error and
  relaxes at h_j + σ₀ ≥ 7.38·10⁻³; the adaptive prior's slowest rate λ_1^− = 2.68·10⁻⁵, 1.51·10⁻⁴, 2.78·10⁻⁴, 3.06·10⁻⁴ for η = σ₀/10,
  σ₀, 10σ₀, ∞. (v)'s radius: C_N ≥ 0.341 (sampled), so ρ ≤ 2.0·10⁻⁵, 1.1·10⁻⁴, 2.0·10⁻⁴, 2.2·10⁻⁴, while the K₀-atom runs start 94–6·10⁴
  times outside, the warm clouds 2.2·10⁴–3.6·10⁵ times and the separated starts 131–3344 times. The fixed prior's bias: fit error 9.78·10⁻³, within 4.8% of the first-order formula in norm, 0.85 and 0.73 of it in the weak
  modes 1 and 3; the error is first order in the shift (0.58% in norm, 0.98 and 0.97 in modes 1 and 3, at 1/8 of the shift).
- K₀-atom loop from Θ*(q), 6000 slow units: every mode above 10⁻⁶ (13 of 36) decays at λ_j^− (log-slopes 0.981–1.003 of it; ratios to the
  pair solutions 0.82–1.12, constant to 1.04% from s = 1000 on for η ≥ σ₀); fit error 9.78·10⁻³ → 3.6·10⁻⁴, 1.8·10⁻⁵, 8.0·10⁻⁶, 6.7·10⁻⁶.
- m = 15, 60, 240, 510 from the warm clouds, four seeds, 3000 slow units. Fixed prior: the lift of Θ*(q) (‖ȳ‖ = 2.10·10⁻², W² = 4.4·10⁻⁴).
  Ungated: the prior dragged by the transit (‖ȳ‖ 0.032–0.32 at η = σ₀, 0.048–0.71 at η = ∞, at s = 500); one m = 15 system at η = ∞
  on a plateau of the unregularised loss (Hessian eigenvalues −1.8·10⁻⁵, 7.2·10⁻⁵; fit error 3.46·10⁻² at s = 23000), which the
  unregularised gradient flow crosses to the exact representation (below 0.9 of the plateau's loss after 4.16·10⁴, half at 4.78·10⁴, below
  10⁻¹² at 7.97·10⁴ more slow units, the crossings located as events). Transit gate (opens at s = 88–149): ‖ȳ‖
  4.9·10⁻³–3.7·10⁻² at η = σ₀ against the K₀-atom run's 1.3·10⁻³. Settled gate (opens at s = 1211–1304, the systems replicated: rms
  separations ≤ 7.9·10⁻⁴): ‖ȳ‖ and the fit error equal to the K₀-atom run's at s − t_open for every m (per-seed ratios 0.995–1.000 at
  η = σ₀, 0.968–1.001 at η = ∞: the invariance of the replicated configurations); fit errors 1.6·10⁻⁴ and 3.4–3.5·10⁻⁵, W² 6.7–7.2·10⁻⁶
  and 9.3–9.9·10⁻⁷; in the seed means the settled gate is the lower except the fit error at η = σ₀ and m = 510 (1.64·10⁻⁴ against the
  transit gate's 1.55·10⁻⁴), and system by system the transit gate ends lower in all three of ‖ȳ‖, W² and the fit error in 6 of 32.
- Separated starts (the structure behind (v), beyond its ball; five draws): ‖ȳ‖ at s = 1750 departs from the K₀-atom run's by +0.13% and
  +3.1% (rms separation 0.02 and 0.1) at η = σ₀ and +0.77% and +18.5% at η = ∞ for m = 15 in the draw means (m = 510: +0.09%, +2.2%,
  +0.63%, +15.5%): second order in the separations draw by draw (factors 23.8–25.0), consistent with the same expected size at both m
  (the means, higher at m = 15, within 0.9–1.7 combined standard errors; five draws), the spread 3.3–6.0 times larger at m = 15.
- The limits of (i): a non-replicated stationary saddle on a four-bump target (one class split, held by a leak force 5.8·10⁻³); from the
  exact split fit the flow collapses the class (width 1.05, fit error 2.23·10⁻² at s = 4·10⁴), and the merged width coordinate then runs off
  like ln s, the width tending to s_max like 1/s: the solution is unbounded, its physical image approaching the minimiser of the loss on the
  face s = s_max (fit error 1.943·10⁻², a Karush–Kuhn–Tucker point with ∂L/∂s = −1.3·10⁻⁶); beyond the saturation (weighted amplitude
  4.8 > a_max) the fixed prior's solution converges while the adaptive prior's amplitude coordinate grows like 2 ln s; a prior permuted among
  atoms of equal mass is stationary at the exact representation.

Tests (`tests/test_adaptive_prior.py`): the joint drift is the gradient flow of 𝒥 in the stated metric (finite differences); the exact
representation stationary, the fixed prior's stationary point not; the linearisation's spectrum (σ₀ and the pair rates λ_j^±) and its
self-adjointness in the gradient-flow metric; the first-order bias; the transit gate and the settled gate, with 𝒥 non-increasing; the
non-replicated saddle, the permuted prior and the adaptive prior's run-off (the force on the amplitude tending to the positive multiplier of
the loss's face minimum, and the run-off law (a_max − a)s → a_max²(σ₀ + η)/(4ηG_αF_∞)).

## Proposition D‴.24 — the floor is the leak's bias; an adaptive prior without a mass constraint (`exp_adaptive_floor.py` → `chaos/adaptive_floor.md`)

- Averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit (0.05σ_crit and 0.5σ_crit for the law), the nearest-atom leak toward the
  matching prior unless stated; the warm clouds (four seeds), m = 15, 60, 240 against M = 510; the checks at the sampled weights from the
  stored final states.
- The floor as the bias (no dynamics): F equals Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ to the finite-difference error (1.1·10⁻⁶, 1.9·10⁻⁷, 1.2·10⁻⁸; exact
  at q); the law (B) per unit of 1/m − 1/M is 1.356·10⁻², 6.166·10⁻², 0.2291, the sum of its mode terms, at 0.496, 0.753, 0.839 of the
  bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ), itself σ₀/4 times the adaptive prior's coupled-distance law 46.32. Modes 2–6 carry 99.4%, 99.0%, 97.9% of
  the law and 97.2% unweighted; the weight (a band-pass centred at h = σ₀) moves the shares of the modes with h_j/σ₀ ∈ [1/5, 5] from 21%,
  97%, 78% to 42%, 99%, 90%. The law's value at m = 15, 60, 240: 3.99·10⁻³, 9.07·10⁻⁴, 1.36·10⁻⁴. The parameter split
  S = H_L(H_L + σ₀I)⁻¹S^e holds to ≤ 1.4·10⁻⁷; the adaptive prior's coupled-distance law is 1.44, 2.17, 5.27 times the fixed prior's.
- Two notions of admissibility. Representable (w > (0.25, 0.175, 0.2)): every sample but (8, 4, 3)/15, on the boundary. D‴.5-admissible
  for the fixed prior: 7 of the 12 m-systems; (6, 5, 4)/15, (6, 3, 6)/15, (24, 15, 21)/60, (8, 4, 3)/15, (25, 17, 18)/60 are saddles
  (margins −3.48, −1.08, −2.16, −6.14, −3.98·10⁻³, class 1, w₁ ≥ 0.4); the references and the others ≥ 3.78·10⁻³.
- Fixed against adaptive priors, 6000 slow units: the fixed prior's ‖f̂_m − f̂_M‖² stays (3.0·10⁻⁴–5.5·10⁻³, 1.1·10⁻⁴–9.9·10⁻⁴,
  4.8·10⁻⁶–8.2·10⁻⁵ at m = 15, 60, 240); at the seven D‴.5-admissible samples it is the exact bias at the sampled weights (ratio 1.0000),
  the delta method missing it by 0.86–1.72; the five saddles have split (class diameters 0.64–1.50) at 0.81–1.01 of the replicated
  floor. The adaptive prior's decreases without a floor (at s = 6000 at most 2.0·10⁻⁶ at η = σ₀ and 5.4·10⁻⁸ at η = ∞ over the
  representable samples); the boundary sample decays like 1/s (6000·‖f̂_m − f̂_M‖² = 0.88 and 0.42 against the predicted 0.8 and 0.4); the
  counts are conserved.
- The parameter level: the coupled distance D is 0.944–1.000 of the replicated stationary points' (fixed: Θ*(w); adaptive: ϑ^e(w)), and
  1.72–4.81 times larger under the adaptive prior; W² is dominated by the cross-class transport of ½‖w^m − w^M‖₁ (0.82–0.99 of it fixed,
  0.67–0.99 adaptive), so the two priors' W² differ by 1.00–1.32 and their within-class parts by 1.71–4.84.
- The design table (the misplaced prior, 3000 slow units, m = 15, 60, 240): nearest atom and fixed, fit error 6.8·10⁻²–2.0·10⁻² and floor
  2.5·10⁻³–5.9·10⁻⁵ (the same five samples saddles, margins −0.38 to −5.76·10⁻³); nearest atom and adaptive, 8.9·10⁻³–5.5·10⁻⁴ and
  8.0·10⁻⁵–1.1·10⁻⁷; OT and fixed, 9.78·10⁻³ and no floor; OT and adaptive (settled gate), 1.6·10⁻⁴ and below 5·10⁻¹¹. The OT designs'
  coupled distance is the largest (6.0, 3.5, 1.0: the transit's mismatch) with W² ≤ 2.1·10⁻⁸. At s = 1250, before the settled gate's
  effect, the nearest-atom adaptive prior leads at m = 60 and 240 (fit error 3.8·10⁻³ and 2.8·10⁻³ against 9.2–9.8·10⁻³), not at
  m = 15 (1.75·10⁻²); from s ≈ 2000 the OT design leads.
- The real loop (k_e = 25, 3000 time units, seed 0's weights): the fixed prior's rms tracking error 1.17–1.18·10⁻³ (m = 15, a saddle
  orbit) and 3.43–3.44·10⁻⁴ (m = 60), at the averaged prediction (ratios 1.001, 0.996 at the end); the adaptive prior's (η = ∞)
  1.16·10⁻³ → 3.27·10⁻⁴ and 3.40·10⁻⁴ → 1.52·10⁻⁴ along the averaged loop (0.987, 0.992); from the exact representation the five-period
  rms is at most 9.0·10⁻¹².
- The instance's probabilities (exact): outside W 0.361, 0.074, 0.0021 and on its boundary 0.118, 0.053, 0.0012 at m = 15, 60, 240, against
  Hoeffding's 1.87, 0.60, 0.036; a permutation of the classes representable for 31%, 87%, 100% of the mass outside W; the fixed prior's
  saddles 0.369, 0.260, 0.109 (comparable to the adaptive prior's failure set at m = 15, much larger at m ≥ 60). The adaptive prior's
  expected limiting floor (the constrained minima on the branch of ϑ^e(q) outside W; the global infimum is 0 where a permutation is
  representable) ≥ 7.3·10⁻³, 2.2·10⁻⁴,
  6.6·10⁻⁷: 1.8, 0.24, 0.005 times the fixed prior's law.
- Monte Carlo (600 nested samples, the exact replicated stationary points): the fixed prior's floor is 2.78 ± 0.20, 1.53 ± 0.11,
  1.11 ± 0.07 times the law (2.89, 1.56, 1.07 for D‴.5-admissible pairs); E D/law 0.68, 0.93, 0.95 (fixed) and 0.53, 1.17, 1.26
  (adaptive: a far heavier tail, one sample at 55.7 times the law at m = 240, 1.17 without it — the means not reliable), the adaptive
  prior's 2.9–3.7 times the fixed prior's.
- In two dimensions (`exp_adaptive_floor_2d.py` → `chaos/adaptive_floor_2d.md`; four bumps with vector amplitudes, p = 5, σ₀ = 0.05,
  m = 16, 32, 64 against M = 128, six seeds): F equals the bias form to 6.0·10⁻⁷, the law 0.310 per unit of 1/m − 1/M at 0.215 of the bound;
  every sample D‴.5-admissible (margins +0.016 to +0.047) and the fixed prior's floor equal to the exact bias at the sampled weights in all
  18 systems (ratio 1.0000); 7 of the 18 systems representable (w > (0.200, 0.150, 0.175, 0.150)), where the adaptive prior's
  function-level distance falls to
  1.0·10⁻⁷–4.1·10⁻⁵ (η = σ₀) and 1.2·10⁻⁹–1.0·10⁻⁵ (η = ∞: the atoms reset to the cell means every slow unit), at most 0.019 and 0.006 of
  the fixed prior's, elsewhere to 0.023–0.64 of it (still falling, toward the constrained minimum at the sampled weights);
  the real loop (k_e = 25, m = 64; the adaptive prior at η = ∞) at the averaged prediction (ratios 0.979–1.009) and invariant from the exact representation
  (2.0–4.7·10⁻¹⁰).

Tests (`tests/test_adaptive_floor.py`): the nearest-atom adaptive flow is the gradient flow of 𝒥_∘ with the cell masses (finite
differences); 𝒥_∘ is the minimum of the OT energy over the masses, strictly below the energies with masses moved by ±1/m; the floor law
as the bias derivative, below the bound, the weight moving the shares toward h ≈ σ₀ against the unweighted ones at two σ₀; at the
D‴.5-admissible sample (5, 4, 6)/15 the exact representation stationary under the adaptive prior and the fixed prior's stationary point a
strict local minimiser, biased; at the saddle (6, 5, 4)/15 the exact representation's margins σ₀; the parameter split; 𝒥_∘ decreasing
across a cell crossing; an empty cell staying empty with its atom frozen. `tests/test_planar_loop.py`: the planar averaged loop (the scalar
one at n = 1; the learning term against finite differences of L̄ and the orbit average of the real law), the exact representation of
representable weights stationary for the adaptive prior, the floor law as the bias derivative in two dimensions, and the finite-η branch
of its integrator second order against DOP853.

## Corollary D‴.26 — perfect tracking with an adaptive prior in the real loop (`exp_adaptive_floquet.py` → `chaos/adaptive_floquet.md`)

- The three-bump instance at σ₀ = 0.15σ_crit, x_m = 2 sin t, the sampled weights w = (6, 5, 4)/15 (amplitudes A_k/w_k = 2.5, −2.1, 3.0;
  representable, and a saddle of the fixed prior), the nearest-atom leak toward an adaptive prior in the real loop with the prior in the
  state; η = σ₀ and ∞, k_e = 25 and 50; the monodromy by central differences of the period map (step 10⁻⁵, DOP853 at rtol 10⁻¹³).
- Every lift returns after one period to ≤ 1.7·10⁻¹⁵. Separations: at m = 15 and 30 the 36 and 81 exponents equal σ₀/k_e to ≤ 1.7·10⁻⁶
  (relative). Collective exponents: the other 18 (η = σ₀) or 9 (η = ∞) equal the weighted K₀-atom loop's to 1.4·10⁻⁵–6.4·10⁻⁴.
- Part C, the variational equation on the collapsed subspace (its generator the real loop's linearisation to 5.8·10⁻¹¹, its exponents the
  differences' to 4.1·10⁻⁵), k_e = 12.5–400, on x_m = 2 sin t (where the first-order averaging correction vanishes identically) and on
  2 sin t + 0.8 sin 2t (where it does not: ‖Ā₁‖ = 5.6·10⁻³, skew): every exponent below its first-order prediction λ_j^±(w)/k_e here —
  the slowest mode by 13.6% and 3.84% at k_e = 25 and 50 (η = σ₀) — the deviations falling by → 4 per doubling of k_e (O(k_e⁻²) relative),
  the residuals of the second-order prediction by → 16 (O(k_e⁻⁴): the odd orders vanish; 4.9·10⁻⁷ at k_e = 400). The second-order
  coefficients are all negative here, c₁/λ₁ = −100.6 in the slowest mode, of which the error filter's part −102.1 = −0.989ω₁², ω₁ = 10.16
  (measured −100.5 at k_e = 400); with the adaptation gain scaled by 30 two are positive (+4.17, +12.65) and those exponents lie above
  their first-order prediction by c_j/(λ_jk_e²); at η = 10σ₀ Ā₀ is self-adjoint in diag(I, (σ₀/η)I) only. The error channel's exponent by Liouville's formula within O(k_e⁻³) of k_e − tr H_L(w)/k_e.

Tests (`tests/test_adaptive_floquet.py`): at the matching weights, η = ∞, m = 6: the return of one period; the nine separation exponents
at σ₀/k_e; the collective exponents equal to the weighted K₀-atom loop's, within the averaging correction of the prediction. At the
weights (2, 1, 1)/4 with the prior in the state, η = σ₀, m = 4: the real loop's Jacobian at the lift maps the separations to −σ₀/k_e times
themselves and the collapsed subspace into itself as the collapsed variational equation; the monodromy's 3 separation exponents at
σ₀/k_e and its 18 collective ones the variational ones, within the second-order prediction and below the first-order one (an instance
property). On an asymmetric reference: Ā₁ non-zero and skew, ℓ_jĀ₁r_j = 0, the second-order residual falling like k_e⁻⁴, the measured
second-order coefficient the explicit c₁, the error channel by Liouville's formula; at η = 10σ₀ Ā₀ self-adjoint in the metric and not in
I, Ā₁ skew in it.

## Lemma D‴.26a — the exact slow subspace: Corollary D‴.26(c) and Proposition D‴.30(ii) proved for the linearised loop (`exp_adaptive_floquet.py` Part D → `chaos/adaptive_floquet.md`)

- Corollary D‴.26's instance on the asymmetric reference x_m = 2 sin t + 0.8 sin 2t, η = 0 (K = Ĵ, Λ = σ₀I), σ₀ and ∞, k_e = 12.5–1600.
- Chang's slow subspace: the contraction in 3–9 iterations to an invariance residual ≤ 3.8·10⁻¹⁵ (at η = 0, k_e = 12.5 it stalls at
  rounding, 499 iterations); |P − P₂|/ε³ → 3.34·10³ and |P − P₍₃₎|/ε⁴ → 5.8·10⁴; the slow block's rates the full variational equation's
  to ≤ 1.0·10⁻⁸.
- The parity as an identity: the rates at ±ε equal to ≤ 5.9·10⁻¹² of the largest rate and ≤ 1.1·10⁻⁹ per mode, while their first-order
  deviation is 9.8·10⁻²–8.9·10⁻⁶ of the largest rate and per mode up to 29.6% (the slowest mode at η = ∞: −29.6%, −12.5%, −4.0%, −1.1%
  at k_e = 12.5–100); a mismatch in the error channel breaks the equality at O(1/k_e) (1.1·10⁻² → 7.8·10⁻⁵); the Magnus parameter
  ≤ 0.48 (0.004 at k_e = 1600).
- The slow noise: the decoupling residual ≤ 3·10⁻¹⁶; b against b₀ = −(K̇ + aK)ᵀ at O(1/k_e) (63% → 0.74%) and against b₀ + εb₁ at
  O(1/k_e²) (|b − b₀ − εb₁|/ε² → 3.34·10³); at t = 0 Cov(v)·2ε³ → 1, E[vY′ᵀ] → −bᵀ at O(1/k_e), Cov(Y′) → Σ̄ (0.35 → 9.4·10⁻⁵ at η = ∞;
  O(ε²) to k_e = 200–400, O(ε) beyond). The weakest mode's deviation (−34% → −6.6·10⁻⁵, O(ε²)) is the slow noise's attenuation,
  ⟨bbᵀ⟩ below Q̄ (−56% → −1.1·10⁻⁴), plus the rate's correction (+49% → +4.6·10⁻⁵): 2.4–2.5 to 1 and of opposite sign at η = ∞ (12–23 to 1
  at η = 0). Cov(Y) − Cov(Y′) is the jitter (ε/2)K(0)ᵀK(0) at O(ε) and with its next term at O(ε²) (0.30 → 0.001 at k_e = 100–1600).
- The floor: at t = 0 0.763 → 1.0000 of the slow-plus-fast value; over the period its average 0.766 → 0.99997 of Proposition D‴.30(ii)'s
  formula, pointwise in [0.9952, 1.0091] of it at k_e = 1600 — the predicted range of the fast jitter (ε/2)tr(HĴ(t)ᵀĴ(t)): the formula
  holds as the period average.

Tests (`tests/test_adaptive_floquet.py`): Chang's P with the exact slow multipliers and the ε³ truncation on the instance; on a synthetic
block with an O(1), non-symmetric Λ the truncations at ε³ and ε⁴, a flipped Λ and a dropped cubic term caught; the rates exactly even in
ε at η = σ₀ and ∞, a mismatch breaking it; the slow noise decoupled, b₀ against K̇ by the chain rule and the plant's own slope by
differences, b₁ at ε², the covariances at their orders.

## Proposition D‴.27 — a projected adaptive prior (`exp_adaptive_prior.py` Part H → `chaos/adaptive_prior.md`)

- The box 𝒵 (unconstrained coordinates): α ∈ [−5.889, 5.889], γ ∈ [−4.417, 4.417] (amplitudes and centres to 0.9 of a_max and c_max),
  κ ∈ [−6, 6]; the instance's exact atoms inside.
- Beyond the saturation (Part G(b)'s target, weighted amplitude 4.8 > a_max; the K₀-atom system from the exact atoms with the prior there):
  the projected prior converges at η = σ₀ (the projected midpoint rule, positions projected inside the step, dt = 0.25 and 0.125, to
  s = 10⁵; the two end states within 3.4·10⁻¹², their shared fixed point; the trajectory at s = 300 within 6.3·10⁻⁷, 1.9·10⁻⁷, 3.8·10⁻⁸ of
  a reference for dt = 0.25, 0.125, 0.0625) and at η = ∞ (LSODA) to the same bounded stationary point (drift ≤ 3.5·10⁻¹⁵): the prior's
  amplitude atom on the face (5.889), the particle's at 6.728, the fit error 0.1415.
- The particles' largest G⁻¹-distance from 𝒵 along the run: 0.92 (η = σ₀), 0.89 (η = ∞), against the a-priori bound B/σ₀ ≤ 3686
  (B ≤ 26.08 with the closed-form bound 4.658 of L_ϑ^G; 4.044 by search).
- The box's size: the amplitude bound 0.75, 0.8, 0.9, 0.95, 0.99 of a_max (boxes containing the exact atoms) gives the fit errors 0.1676,
  0.1611, 0.1415, 0.1268, 0.1101 — between the fixed prior's 0.1750 and the unprojected adaptive prior's 0.1050, reached only as its
  coordinate runs off; 0.7 (0.1736) excludes the exact atoms, and the tighter 0.65, 0.6, 0.5 are worse than the fixed prior (0.1792,
  0.1845, 0.1949): a box should contain the representation.
- At η = ∞ the plan is Lloyd-consistent, not jointly optimal: on six particles in a consistent, non-optimal pairing the selected energy
  falls (largest increment −3.4·10⁻⁶) while the reduced energy min_Z𝒥 rises from 1.251·10⁻³ to 1.574·10⁻³ before the plan becomes optimal.
- Part G(a)'s four-bump target: the box stops the merged atom's width coordinate at the face κ = 6 from s = 10⁹ (unprojected: κ = 9.75 at
  s = 10¹⁰ and growing like ln s), the fit error within 9.3·10⁻⁷ of the face minimiser's.
- Part C's K₀-atom loop at the misplaced prior (η = σ₀, 6000 slow units) is identical with and without 𝒵: the box is never approached
  (margins at least 16.9% of its width), so the runs coincide by construction.

Tests (`tests/test_adaptive_prior.py`): beyond the saturation the projected prior at η = ∞ keeps the particles within 2 of the box and
converges, while without the box the prior's coordinate passes 11 by s = 10⁴; at η = σ₀ and dt = 0.125 the prior stays in the box with its
amplitude atom on the face (the step that projects the velocity at the midpoint stalled it inside), the particles bounded and 𝒥
non-increasing; at the misplaced prior near the exact representation the projected and unprojected runs are identical.

## Proposition D‴.28 — the floor removal is a gauge property (`exp_gauge_floor.py` → `chaos/gauge_floor.md`)

- The three-bump instance plus a fourth bump ε·exp(−(x − 0.8)²/(2·0.3²)) that three kernels cannot represent, ε = 0–0.5, masses 1/3,
  σ₀ = 0.15σ_crit: the weighted approximant ϑ^b(q) continued in ε, a strict local minimum of the unregularised loss with the relative fit
  error 3.6·10⁻³–2.5·10⁻² (ε = 0.05–0.5).
- The gauge: at six weight vectors of W_g the critical point continued in w from ϑ^b(q) is the gauge transform (to ≤ 1.7·10⁻¹¹, f̂
  unchanged to ≤ 5.5·10⁻¹³); the adaptive prior's sensitivity 1.4·10⁻¹¹–2.7·10⁻¹⁰ of the fixed prior's; P(w^m ∈ W_g) = 0.52, 0.87,
  0.997 at m = 15, 60, 240 (ε ≤ 0.3); the branch reaches the saturation near ε = 0.68.
- The fixed prior's sensitivity equals σ₀Ĵ(Ã + σ₀I)⁻¹g̃ (the full Hessian) to ≤ 2.2·10⁻⁹, against 7.0·10⁻⁴–5.6·10⁻³ for the Gauss–Newton
  form; its floor law 6.17·10⁻² (ε = 0) to 8.52·10⁻² (ε = 0.5) per unit of 1/m − 1/M.
- The residual curvature: the amplitude row ≤ 7.4·10⁻¹⁴; the margins σ₀ − (q_k/w_k)λ⁺_k exact to ≤ 2.1·10⁻¹⁴ over every count vector of
  m = 15 and 60 in W_g, none negative; the dynamics' separation eigenvalues at a replicated configuration equal to them to ≤ 5.8·10⁻¹².
- The dynamics (m = 15, 60, 240 against M = 510, four seeds, 6000 slow units, ε = 0.1 and 0.3): the fixed prior's floor within 8% of its
  value at s = 600, its fit error above the approximant's (its bias; split classes at its saddles); the adaptive priors' floor falling (m = 60,
  ε = 0.1: to 1.0·10⁻⁹–1.4·10⁻⁸ at η = σ₀, 3.2·10⁻¹¹–3.8·10⁻⁹ at η = ∞) and the fit error tending to the approximant's (median relative
  excess 1.2·10⁻⁴, 1.2·10⁻⁴, 1.6·10⁻⁵ at m = 15, 60, 240 for the samples 0.02 inside W_g); at the sample outside W_g the systems approach the
  box minimum on the branch (1.107 and 1.586 times the approximant's fit error), at the sample 6·10⁻⁴ inside it the approach is slow.

Tests (`tests/test_gauge_floor.py`): at ε = 0.3, the gauge transform a critical point at two sampled weight vectors with the same
function, reached by continuation from ϑ^b(q); the adaptive prior's sensitivity negligible, the fixed prior's the full-Hessian formula and
not the Gauss–Newton form; the residual curvature's zero amplitude row, its scaling by q_k/w_k, the margins and the dynamics' separation
eigenvalues; adaptive-prior runs (η = ∞ and σ₀) converging to the gauge transform with the approximant's function.

## Proposition D‴.29 — second-order plants: the sliding-variable reduction; a cogging-torque servo (`exp_servo_cogging.py` → `servo/servo_cogging.md`)

- The servo ẍ = −τ_cog(x) + u, τ_cog(x) = 0.4 sin(3x + 0.3) + 0.1 sin 7x, x_m = sin t, s = ė + 5e; K₀ = 4 saturated Gaussian kernels on
  the position, masses 1/4, σ₀ = 0.01.
- ϑ^b(q), the nondegenerate local minimum BFGS reaches from a fixed start (no exact representation): fit error 1.46·10⁻², Ã spanning
  2.9·10⁻⁴σ₀–30σ₀ (four modes below 0.1σ₀; the adaptive prior's slowest rates 1.4·10⁻⁴σ₀ and 2.9·10⁻⁴σ₀ against the fixed prior's
  ≈ σ₀), λ⁺_k ≤ 0.0098σ₀, the floor law 1.29·10⁻². Globally the loss's infimum lies at the saturation: of 120 random starts 1 ends at
  ϑ^b(q) and 33 below its fit error, all saturated; the best 6.39·10⁻³, approached with the saturated amplitudes held at |α| = 6–14
  (7.14·10⁻³–6.40·10⁻³), not attained. A misplaced prior (phase 0.5, the engineer's model, fit error 0.200) whose replicated stationary
  point is a D‴.5 saddle (margin −3.97·10⁻³, fit error 0.0545).
- The averaged loop (m = 16, 64, 128 against M = 256, three seeds, 2000 slow units) from clouds around ϑ^b(q): the fixed prior's floor
  1.85·10⁻³, 2.24·10⁻⁴, 4.49·10⁻⁵ (fit error 0.174, 0.066, 0.038); the adaptive prior's (η = ∞) 38–780 times lower, the fit error 0.0298,
  0.0148, 0.0147 against the approximant's 0.0146; the OT leak at the approximant with no floor; the σ-modification toward 0 biased
  (0.328); no leak (σ₀ = 0) as good: 0.0180, 0.0147, 0.0146 (floors 7.7·10⁻⁶, 2.4·10⁻⁷, 6.6·10⁻⁸) — noise-free, the floors and biases
  the leaks remove are those a fixed leak creates. Misplaced: the fixed designs split their classes and end at 0.048 (OT) and 0.10–0.048
  (nearest atom); the adaptive prior is unaffected at η = ∞ (its prior never acts from this start), slower at η = σ₀ (0.035–0.019).
- From the engineer's model (clouds around it, the prior at it): the fixed designs 0.102, 0.060, 0.048 (OT 0.048); the adaptive prior
  0.027, 0.025, 0.022 (η = ∞) and 0.032, 0.028, 0.028 (η = σ₀) at s = 2000; no leak 0.0152, 0.0147, 0.0147. The collapsed K₀-atom flow
  from the engineer's model at η = ∞: 0.0297, 0.0213, 0.0169, 0.0155 at s = 600, 2000, 10⁴, 10⁵.
- Phase-plane kernels (five isotropic kernels): on the circle 5 exact zero Gauss–Newton eigenvalues and 8 below 10⁻⁶σ₀; ellipses ×1.1
  and ×2: 8 and 6 below 10⁻⁶σ₀, none exact; a two-harmonic orbit 3; centred kernels 16 exact of 20.
- The real servo (k_s = 25, 1000 time units; records and leak targets every 0.02 time units): position-error rms 1.7–2.1·10⁻⁵ (adaptive
  designs, matching OT, no leak), 7.8·10⁻⁵ (misplaced OT), 9.1·10⁻⁵–1.8·10⁻⁴ (fixed nearest atom), 5.4–5.5·10⁻⁴ (σ-modification), equal
  to the prediction (d/dt + λ)⁻¹(d/dt + k_s)⁻¹r(x_m) to 1.0001–1.0029 at t = 500 and 1000 (with the slaved r/k_s: 0.931–0.996, the
  filter's attenuation of the residual's broad fine scale).

Tests (`tests/test_sliding.py`): the phase-plane averaged loop is the gradient flow of L̄; on a circular orbit isotropic phase-plane
kernels in general position have one zero Gauss–Newton eigenvalue each; the servo's parameters follow the averaged loop over ten periods
(within 3% at k_s = 25, closer at 50), its position error within 1% of both predictions.

## Proposition D‴.30 — measurement noise at perfect tracking: coherent, high-passed, and the errors in variables (`exp_measurement_noise.py` → `chaos/measurement_noise.md`)

- The three-bump instance at the lift of the exact representation of q = 1/3 (σ₀ = 0.15σ_crit, x_m = 2 sin t, k_e = 25), the measurement
  y = x + ν in the controller and the adaptive law.
- The linear theory: tr Q = 8.784 (the derivative Gram 8.105, the slope term 0.679); the parameters' floor per unit N₀ 0.1597 (fixed
  prior), 0.1627 (η = σ₀), 0.1690 (10σ₀), 0.1757 = tr Q/(2k_e) (η = ∞); the weakest mode's parameter variance 0.0865 → 2.09 (×24.1;
  ×28.0 in the exact linearised loop at k_e = 25); the fast parts tr H² = 0.0657, ⟨|Ĵ|⁴⟩ = 0.299.
- The errors in variables (OU noise, B = 100): the averaged stationary function-level bias 491 (fixed) and 2397 (adaptive, every η > 0)
  per unit (Eν²)², above the floor of the parameters beyond the noise std 0.0026 and 0.0012; the parameters' fluctuation changes it by
  less than 10⁻⁵ (2397.68); the exact periodic second-order mean of the linearised loop 504 and 2653, the excess O(k_e⁻²) (10.7%, 2.65%,
  0.66%, 0.17% at k_e = 25–200, adaptive); the bias's growth in k_e at the local exponent 1.93–1.99 for B ≥ 20k_e, 1.69 and 1.47 at
  B = 8k_e and 4k_e.
- The exactly linearised loop: its periodic stationary covariance at the prediction to 0.976–0.980 (k_e = 25), 0.994–0.995 (50),
  0.9986–0.9988 (100), 0.9997–0.9998 (200), above 1 at 800 (+1.3·10⁻⁵): the relative error is O(1/k_e), its first-order term from the
  slope a in the noise input — with the noise input −k_eν the deviation is numerically even in 1/k_e (successive ratios 3.67 → 4.00 to
  k_e = 800; not proved).
- Invariance, not attraction: at the lift the separations' second-order drift λ_max(S_k) = 0.0148, 0.0101, 0.298 per Eν² in the three
  classes, the third's separations growing there above the noise std 0.031; along the noisy trajectory the mean error also carries the
  tracking error of the developing errors-in-variables deviation (rms 0.121 → 0.602 per Eν², λ_max(S₃) 0.307 → 0.349 over 300 time
  units). In the real loop (the lift m = 15, eight antithetic pairs of noise paths of 300 time units per std) the separations' rates are
  −2.50·10⁻⁴, +2.74·10⁻⁶ and +4.53·10⁻⁴ per time unit at std 0.01, 0.03, 0.05; λ_eff = (rate + σ₀/k_e)/Eν² is 0.3271 ± 0.0038,
  0.3175 ± 0.0037, 0.2943 ± 0.0034 against the lift's 0.2981 and the transient's 0.3264, 0.3270, 0.3281; the fit's small-noise limit
  0.3291 ± 0.0033, the next order −13.8 ± 2.0; the rate vanishes at std 0.0299 over the window (0.0308 by the lift's drift, 0.0294 by the
  transient's); the long-time threshold is open.
- The real loop (12 noise paths of 300 time units, common to the four cases): within 1.75 standard errors of the linear prediction plus
  the mean part at both noise levels, the mean part 42% (std 0.003) and 89% (std 0.01) at t = 300, the priors 6.0–6.4% apart (predicted
  5.3–6.1%); the lift m = 15 on a common noise path identical to the K₀-atom run to rounding (5.6·10⁻¹⁶ in float64).

Tests (`tests/test_measurement_noise.py`): the closed form against the Lyapunov solution, the trace identity at η = ∞ and the floor
increasing in η; the exact periodic floors of the linearised loop at k_e = 100 for η = 0, σ₀, ∞; a lift m = 6 under a noise path staying
replicated and following the K₀-atom loop; a separation of the weakest class decaying at std 0.01 and growing at 0.05 (a sign test);
the separations' drift against an independent quadrature of its two moments; along the second-order transient the mean parameters moved
by the averaged force and the mean error carrying the tracking error of the developing deviation; the OU noise's stationary variance and
correlation, and the law seeing the measurement; the responses (H + σ₀I)⁻¹ and H⁻¹ to a constant force, and the quasi-static force
against the full second-order mean's noise part; the exact periodic mean against the filtered-Gram response under both priors, the fixed
prior's branch and its leak in the bias; the mean function's curvature term against an independent contraction.

## Proposition D‴.31 — the adaptive prior is a washout leak; the desired regressor makes the noise's own persistent force a least-squares force (`exp_eiv_remedy.py` → `chaos/eiv_remedy.md`)

- The instance of Proposition D‴.30 with an Ornstein–Uhlenbeck noise of rate B; MR the measured regressor (Proposition D‴.30), DR the
  regressor at the reference, DR* the regressor and the compensation at the reference (the desired compensation adaptation law).
- The adaptive prior's function-level errors-in-variables bias per unit (Eν²)² at k_e = 25, B = 100: from the noise's own moments 2397
  (MR), 0.471 (DR), 0.0526 (DR*); with the parameters' fluctuation (the full second-order mean) 2398, 0.486, 0.0481 — the parameters'
  part ≤ 0.03% for MR, −8.7% to +3.3% for DR and DR* over the grid except +56% for DR at k_e = B = 100, where its noise part nearly
  cancels; its own response 3.9·10⁻⁴, 1.1·10⁻⁴, 2.9·10⁻⁵, 7.4·10⁻⁶ at k_e = 25–200 (∝ k_e⁻²: on the one-dimensional reference the law's
  leading term, a line integral along the orbit, vanishes), its share of DR's bias not falling monotonically. This part is the classical
  second-order bias of a constant-gain stochastic approximation (Benveniste–Métivier–Priouret, Dieuleveut–Durmus–Bach, Box's curvature
  bias); the desired regressor does not remove it. MR over DR 3.0·10³–2.1·10⁸ over k_e = 25–200 and B = 100, 1000;
  MR's growing like (k_eB/(k_e + B))², DR's and DR*'s at most 1.18, below ¼‖∂²_xf‖²_ρ = 1.31; DR's noise part at its least-squares bound
  to 0.01%; k_e²‖P E[e]‖²_ρ 2394 for MR (zero by construction for DR, DR*), the mean tracking error itself 3.83 → 0.66 (k_e = 25 → 200).
  The noise std at which the bias reaches the floor: 0.0012 (MR), 0.085 (DR), 0.27 (DR*).
- The prior family dz/ds = −η(z − T̂) − κ(z − z_a): RH + B = I to 1.1·10⁻¹⁴; the weakest mode's response 3269 → 136 as its bias factor
  goes 0 → 0.959, the trade binding on the 3 of 9 modes with h_j < σ₀; MR's bias 2397 → 491, DR's 0.471 → 0.211. Priors whose
  zero-frequency leak does not commute with H (fixed amplitudes or centres with η = σ₀ elsewhere, an active face at η = ∞, the anchor gain
  ηP_amp, a random anchor gain with K + Kᵀ ⪰ 0), R and B from their dynamic blocks: equal to (H + S₀)⁻¹ and (H + S₀)⁻¹S₀, S₀ by
  elimination, to 1.9·10⁻¹⁵; RH + B = I to 1.4·10⁻¹⁴, HR + B − I 0.47–0.79.
- The real loop on common noise paths (k_e = 25, std 0.01, η = ∞, 12 runs of 300 time units): E‖f̂ − f‖² at t = 300 2.43·10⁻⁶ (MR)
  against 2.80·10⁻⁷ (DR) and 2.75·10⁻⁷ (DR*), DR and DR* at their covariance part; the MR − DR run mean at its prediction to 0.5–5.3%,
  0.43–1.83 times the run mean's sampling error.

Tests (`tests/test_eiv_remedy.py`): the washout identity along a real-loop solution at finite η and across a plan switch of the averaged
loop; the collective block's transfer function against the washout formula, RH + B = I along the prior family and, R and B from their
dynamic blocks, for non-commuting leaks (HR + B = I failing there); the response bound under random leaks S ≥ 0, the smallest equivalent
target perturbation and the
projection; the lift exact under DR and DR*, on a common noise path DR's run within the second order of MR's (the difference ∝ the
noise's variance), the regressor's point in the law; DR's and DR*'s noise forces least-squares with their pointwise bounds on g, MR's F₁
three orders of magnitude above; the quasi-static forces against the full second-order mean, the parameters' part above 1% for DR and
DR*, below 10⁻⁵ for MR.

## Proposition D‴.32 — the cohesion leak: a leak on the unexcited subspace, without a prior (`exp_cohesion_leak.py` → `chaos/cohesion_leak.md`)

- The three-bump instance (σ₀ = 0.15σ_crit), the averaged slow loop. (A) At each of the 15 count vectors of m = 15 the lift of ϑ^e(w) is
  stationary (the largest drift 2.8·10⁻¹⁷ under the compact kernel, 6.2·10⁻¹⁰ under the Gaussian) with (iii)'s spectrum to 2.8·10⁻¹¹ (the
  finite differences' error) and 4.6·10⁻⁹ (the Gaussian kernel's cross-class terms): the separations at σ₀w_k, the collective modes the
  unregularised flow's. The Gaussian kernel's bias is 2.1·10⁻⁸, 1.3·10⁻³ and 0.068 (whitened) at R = 1, 1.5, 2, at its first order to
  6·10⁻⁸–1.0% (56% off at R = 3, where the classes overlap). The splitting margins σ₀w_k − (q_k/w_k)λ⁺_k are the Jacobian's to 2·10⁻¹¹; at
  σ₀/8 the class with the margin −1.22·10⁻⁴ splits (3 → 4 classes, the fit error 0.017 → 0.0148). Twins merge under the Gaussian kernel
  and persist under the compact one; a dead class joins a live one at R = 2 and 1.5 and does not move at R = 1, as (v)'s time scale says.
- (I) The twins' contraction on their gauge curve is (v)'s rate to 0.04%; twins 5.3 and 8.2 apart do not move in 4·10⁴ slow units;
  opposite-sign twins with a third class within reach do not merge in 4·10⁵ slow units and end at a biased four-class fit (the fit error
  3.2·10⁻⁴ → 4.5·10⁻³).
- (B) From warm clouds (m = 15, 60, 240) both kernels remove the floor as the adaptive prior does (‖f̂_m − f̂_M‖² at s = 6000
  1.9·10⁻¹⁰–5.4·10⁻⁸ over the representable samples at m = 15, against the fixed prior's floor 3.0·10⁻⁴–2.4·10⁻³); the one sample on
  the boundary of W (seed 3) ends at 7·10⁻⁵ under every leak but the fixed prior.
- (C) From broad starts the annealed Gaussian kernel (R = 2, then R = 1) ends at the target's three classes in 9 of 9 runs, the fit error
  1.8·10⁻⁷–2.3·10⁻⁶ — as the adaptive prior told K₀ = 3 (4.2·10⁻⁹–2.1·10⁻⁶); told K₀ = 2 it ends biased (0.13–0.19), told K₀ = 4 with a bump
  split (5.4·10⁻⁴–5.6·10⁻³); the σ-modification ends at 0.24, no leak at 40–48 classes, the compact kernel without the broad stage at 17–22.
  (J) At a fixed Gaussian scale it does not: 4–6 classes at R = 1, three biased classes at R = 2; three other annealed schedules end at
  three classes.
- (D) The real loop (k_e = 25): the lift returns after a period to 7.8·10⁻¹⁶, and its 45 slow Floquet exponents (36 separations at
  σ₀w_k/k_e) are the predicted ones to 5.5·10⁻⁶; from warm clouds the rms tracking error at t ≈ 2900 is 5.8·10⁻⁴ (m = 15) and 2.4·10⁻⁴
  (m = 60), against the fixed prior's 1.3·10⁻³ and 3.9·10⁻⁴.
- (E) Under the dither τ = 10⁻⁵ the within-class spread over the second half is 0.068 (compact) and 0.063 (Gaussian) against
  √((1 − 1/n_k)τ/(σ₀w_k)) = 0.0635 (the fixed prior's 0.036 against 0.0366; without a leak it diffuses, 0.44–0.56).
- (F) On the cogging servo the leak settles on 7–8 classes with the fit error 0.14%–0.35%, below the four-kernel designs of Proposition D‴.29
  (1.46%, 0.64%) and above the best equal-weight six- and eight-kernel fits.
- (G) σ halved stage by stage from σ₀ to σ₀/64 (the fourth bump, m = 30): no split down to σ₀/8 (the predicted threshold 0.121σ₀), the
  middle class splitting at σ₀/16 in every seed; no class with a positive margin split, all 8 with a margin below −5·10⁻⁵ (the detection
  threshold −ln 250/10⁵) did. A split sheds single particles (the kernel's concavity: shedding one particle costs 0.43 of the quadratic
  cohesion energy, halving 0.74): 8–10 classes at σ₀/64, five to seven of them of one particle, the fit error 0.0028–0.0031.
- (H) Two bumps: two classes, the fit error ≤ 1.6·10⁻⁸; the four-bump f_ε: three classes (0.0174), the fourth bump below the resolution;
  five bumps (Σ_k|A_k|/a_max = 0.875): 6–9 classes, biased (0.045–0.118) — the transit's masses are not representable. (K) Three equal
  bumps: three classes (9.4·10⁻⁵–1.8·10⁻⁴, the Gaussian kernel's bias, removed by the compact kernel: (L)); four alternating bumps: four
  (9.6·10⁻³–1.1·10⁻², a slow transient: (L)); three heavy bumps (0.95): 3–5, biased
  (0.058–0.20); two bumps closer than the broad scale: merged into one stable class (0.098–0.13). The count is not selected by the leak:
  it is the transit's and the schedule's.

Tests (`tests/test_cohesion_leak.py`): the drift equals −G∇J^coh_m (central differences) and J^coh_m decreases along the flow, for both
kernels; the forces sum to zero and equal σ₀ × the kernel mass × the mean-shift vector; at the exact representation for the counts
(3, 2, 2) and (4, 3, 3) the drift vanishes and the Jacobian's spectrum is the predicted one to 10⁻⁹; at the best one-atom fit the top
separation eigenvalue is λ⁺ − σ₀w, stable at σ₀ and split at λ⁺/2; twins contract at (v)'s rate to 1% under the Gaussian kernel and
persist beyond the compact range; under the compact kernel a pair at 0.99R attracts and at 1.01R does not, and the Gaussian pair force at
the distance R is e^{−1/2} of its linear part; single linkage recovers known clusters; the real loop's lift is T₀-periodic and a separation contracts by
e^{−σ₀w_kT₀/k_e} per period to 10⁻⁶.

## Corollary D‴.33 — the desired regressor keeps the replicated configuration attracting under measurement noise (`exp_eiv_remedy.py` Part D → `chaos/eiv_remedy.md`)

- Proposition D‴.30's instance (k_e = 25, OU noise B = 100), per unit Eν². At the lift the separations' second-order drift has the extreme
  eigenvalues [−3.7·10⁻¹, +1.5·10⁻²]; [−2.6·10⁻², +1.0·10⁻²]; [−3.5·10⁻², +3.0·10⁻¹] (MR), [−1.1·10⁻², −5.0·10⁻⁴]; [−1.5·10⁻², −1.5·10⁻⁴];
  [−7.2·10⁻³, −3.4·10⁻⁴] (DR) and [+1.6·10⁻⁴, +3.8·10⁻³]; [+5.1·10⁻⁵, +5.0·10⁻³]; [+1.1·10⁻⁴, +2.4·10⁻³] (DR*): DR's negative definite on the
  instance, DR*'s largest 60 times below MR's (the threshold std 0.238 against 0.031 under the adaptive prior, 0.138 against 0.018 under the
  cohesion leak's σ₀w_k).
- Against k_e: at B = 4k_e MR's threshold std 0.031 → 0.0054 (∝ k_e^{−1/2}, k_e = 25 → 800), k_e times DR's and DR*'s largest eigenvalue
  constant (−0.0038, 0.1245–0.1242); at the fixed B = 100 DR's drift positive definite from k_e = 200 (std 0.185 → 0.121), DR*'s threshold
  0.238 → 0.113, MR's 0.031 → 0.014. ‖S_k‖ ≤ 0.37 of (ii)'s bound: the threshold's floor std 0.069, uniform in k_e.
- Other targets at the lift: DR's centre entries negative for all six, the whole drift negative definite only for alternating-sign separated
  kernels (three same-sign bumps +2.9·10⁻⁴; wide, edge and close bumps up to +5.6·10⁻⁴).
- At the stationary point (the exact periodic second-order mean, adaptive prior): rms E[e] 6.70·10⁻⁴ (DR) and 6.44·10⁻⁴ (DR*), 2.4% and 7.0%
  of the lift's, its parameters' part 6.32·10⁻⁴ and 6.27·10⁻⁴ (the averaged noise-only residual 2.34·10⁻⁴, 2.60·10⁻⁴); the drift indefinite,
  DR's largest eigenvalue +3.68·10⁻⁵ (class 1); under the fixed prior DR's drift negative definite and DR*'s positive definite (std 0.435).
- Along DR's exact second-order transient (5.2·10⁵ time units) the largest eigenvalue turns positive at t = 4103, 2947, 2714 (classes 1–3)
  and reaches +1.05·10⁻⁴ (t ≈ 3.1·10⁴; formally std 1.64); the cross term −8% (DR) and −17% (DR*) of the class-2 drift over [0, 50), DR*'s
  threshold std there 0.277.
- The real loop on common noise paths (four antithetic pairs per family): MR's separations grow (+4.69·10⁻⁴, +1.25·10⁻³ at std 0.05, 0.1;
  every run), DR's decay faster than σ₀/k_e = 2.83·10⁻⁴ (along u −2.965·10⁻⁴ and −3.365·10⁻⁴, the noise part 0.98 and 0.97 of the
  transient's prediction, 0.82 and 0.81 of the lift's; along DR's own top eigenvectors faster by 0.14–0.36% at std 0.05), DR*'s slower (by
  3.2% and 12.6% along its top eigenvector, at the transient's prediction); no DR or DR* run grows.

Tests (`tests/test_eiv_remedy.py`): DR's separation drift equals MR's mean-error term at the lift and is negative definite in every class,
MR's 0.30 in class 3; DR*'s drift is O(1/k_e) (k_e times its largest eigenvalue the same at k_e = 25 and 50, B = 4k_e) and its independent
quadrature's to 10⁻⁴; at the adaptive prior's stationary point the exact periodic mean is the sum of its noise and parameters' parts, the
parameters' part dominant and largely in the tangent space, the drift indefinite; the period maps of `long_transient` continue
`second_order_transient` to 10⁻⁷, and DR's drift turns indefinite within 800 periods.

## Proposition D‴.34 — mass transfer: the replicator weight law with the cohesion leak (`exp_mass_transfer.py` → `chaos/mass_transfer.md`)

- Proposition D‴.32's instances, the weighted network f̂ = (1/m)Σβ_iΦ̃_i, the replicator law dβ_i/ds = γβ_i(u_i − λ) beside the transport
  law. (A) At fixed positions (the five bumps' atoms, two classes short of mass) the transport alone stalls at the fit error 0.388; with the
  replicator it reaches 2.9·10⁻⁴ (γ = 10⁻³) and 3.0·10⁻¹⁴ (γ = 10⁻²) at s = 4·10⁴; with the amplitudes scaled by 1.25 (Σ_k|A_k|/a_max = 1.094)
  it stalls at 0.071 with every class saturating, creeping toward the ℓ¹ ball's least error 0.0684 — the representability criterion (iv)
  for bounded trajectories.
- (B) The five-bump target from broad starts (m = 60, the annealed Gaussian schedule): 6–9 biased classes without the replicator
  (0.045–0.118), 5–7 at γ = 10⁻³ (0.0076–0.0139), five at the bumps in 3 of 3 runs at γ = 10⁻² (6.0·10⁻⁴–1.4·10⁻³; continued, seed 1 levels off at 2.8·10⁻⁴).
  (F) The three heavy bumps of D‴.32 (K) (Σ_k|A_k|/a_max = 0.95): three classes in 3 of 3, 7.8·10⁻⁸–1.6·10⁻⁷ (3–5 biased without); two
  bumps closer than the broad scale stay merged (0.068–0.13): the scale's obstruction is not the masses'. (D) The instance: three
  classes, 2.4·10⁻⁷–9.7·10⁻⁷, as without.
- (G) At the class level under the Gaussian kernel the replicator selects the masses and moves the bias (B's seed 1: 3.85·10⁻⁴ with the
  masses held, 2.63·10⁻⁴ with them free; three equal bumps: 9.39·10⁻⁵ and 9.17·10⁻⁵); under the compact kernel the stationary points are
  exact either way.
- (C) The finest resolution stage (σ₀/64): the same 8–10 classes with and without the replicator — the fitness vanishes at a critical
  point, so the shedding of D‴.32 (G) is untouched.
- (E) The real loop at the instance's atoms (k_e = 25, γ = 2): V = ½e² + (1/(γm))Σβ*log(β*/β) from 0.0418 to 2.1·10⁻⁷, monotone, V̇ = −k_ee²
  to 1.1·10⁻¹⁶ (relative), Σβ = m to 2.4·10⁻¹⁴. The averaged linearisation at the lift of ϑ^e(w), counts (6, 5, 4): 15 zero modes in
  (Θ, log β) — the 12 within-class redistributions and 3 gauges — the rest stable, at γ = 0 Proposition D‴.32's spectrum to 2.1·10⁻¹¹.

Tests (`tests/test_mass_transfer.py`): the weighted loop's drift equals the gradient of its energy in the Euclidean–Fisher–Rao metric (central
differences) and is Proposition D‴.32's at γ = 0; the weights' drift matches an independent quadrature of the fitness; in the real loop at fixed positions V̇ = −k_ee² pointwise from the vector field and Σβ is
conserved; at fixed positions Σβ = m is kept and the exact fit is reached with the replicator and not without it, and not below the ℓ¹ ball's
least error when Σ|A_k|/a_max > 1; at the lift of an exact representation the zero modes in (Θ, log β) are the predicted Σ(n_k − 1)
redistributions and K gauges, a separation keeps −σ₀w_k with γ > 0, and the rest of the spectrum is D‴.32's (γ = 0) and stable (γ > 0).

## Proposition D‴.35 — the split of a class is an S_n-symmetric transcritical bifurcation (`exp_split_bifurcation.py` → `chaos/split_bifurcation.md`)

- Proposition D‴.32's resolution instance (the lift of ϑ^b(q), m = 30, compact R = 1). (A) The middle class splits first (λ⁺ = 2.86·10⁻⁴,
  the threshold 0.1213σ₀); along its splitting direction the residual's third derivative is c₃ = −1.257·10⁻³.
- (B) Near the threshold Newton on the full drift finds the S_q × S_{n−q} split states (q = 1–4) where the reduced equation puts them (0.3%–3.0%
  at μ = 5·10⁻⁷), each a saddle, its eigenvalues converging to ∓μn/(n − 2q) and −μ as μ → 0.
- (C) The reduced quadratic equation names the first excursion's particle in 12 of 12 runs and its time to 1.1%; the linear equation names it
  in 6 and misses the time by up to 22.4%, the flipped one in 2. The first excursion is on the side sign(c₃)v in 10 of 12. At σ₀/16 a bound
  pair leaves in 3 of 6 runs (the top two within 1.2%–11%), a single particle in the others; the cascade stops within a stage once the main
  piece's margin is positive.
- (D) Under Part G's schedule: 8–10 classes with the cohesion leak, 6–7 with the replicator (γ = 0.1). (E) With the replicator at γ = 1:
  5, 6 and 8 classes (the first escapee's mass 0.134–0.136): the rate makes no difference in three seeds. (F) The mass-free leak continued to
  σ₀/256: 4 classes at σ₀/64 in 3 of 3 runs (the fit error 0.0149), then 9, 10, 9 at σ₀/128 and 9, 8, 9 at
  σ₀/256, the fit error 0.0028 — the cohesion leak's end state: the 4 classes were where the schedule stopped.

Tests (`tests/test_split_bifurcation.py`): the full drift along v against the reduced quadratic equation (the remainder cubic); the
single-particle split state a saddle with the predicted eigenvalues; the predicted particle leaving first, on the predicted side; the
reduced crossing against closed forms; the mass-free leak's blocks.

## Proposition D‴.36 — the gauge under dither (`exp_gauge_dither.py` → `chaos/gauge_dither.md`)

- Proposition D‴.32's three-bump instance, one particle per class at the exact representation ϑ^e(1/3), the replicator at γ = 10⁻², the
  positions dithered at τ = 10⁻⁴. (A–B) The masses diffuse along the manifold of exact fits at the predicted rate (the displacement's variance at 0.95, 1.36, 0.90 times the prediction over 64 seeds, within two standard errors); their means drift by about 2.6 standard errors in 5·10⁴ slow units (the limit's drift
  term, computed in Proposition D‴.37).
- (C) Without a leak the masses spread to 0.05 in 2·10⁵ slow units and the first amplitude runs into the saturation (|a|/a_max = 0.9994,
  against 0.80 with fixed weights), the fit error unchanged (the dither's).
- (D) The entropic leak on the masses pins them: at ε = 10⁻³ the spreads 0.0044, 0.0032, 0.0034 against the Ornstein–Uhlenbeck
  prediction 0.0044, 0.0031, 0.0033; at ε = 10⁻⁴ four seeds scatter too much to test it (0.0152, 0.0103, 0.0117 against 0.0139, 0.0099,
  0.0105), the 16-seed test being Proposition D‴.37 (G).

Tests (`tests/test_gauge_dither.py`): the leak a Fisher–Rao gradient; the two-class diffusion's closed form; the leak keeping the exact
representation stationary and stabilising the gauge; the leak's spread against the two-class closed form; no weight motion at γ = 0; D against the projector of the code's
drift Jacobian; run_dithered's noise 2τG dt.

## Proposition D‴.37 — the dither's drift along the gauge (`exp_gauge_drift.py` → `chaos/gauge_drift.md`)

- Proposition D‴.36's instance at τ = 10⁻⁴. (A) The drift from the limit map's second-order identities agrees with the integrated limit map
  to 10⁻³; (B) it is the amplitudes' noise's (each class's drains its own mass), the centres' and widths' 10³ times smaller.
- (C) As γ → 0 the drift over γτ becomes the replicator flow of the entropic free energy (τ/2) log det ∇²_ΘJ (0.3% at γ = 10⁻⁵),
  ℓ = 2Σ log(W_ka′_k) + const; (D) in a saturation layer of width ∝ √γ the drift of W turns outward (its zero at ≈ 0.56δ* as γ → 0), but
  it is the Itô term of the amplitude's driftless diffusion: the layer slows the descent, it does not stop it.
- (E) At γ = 10⁻² the mean drift has three attracting states, each with two classes near their saturation — not where the masses
  settle; (F) the 64-seed ensemble mean follows the mean path; (G) the leak's shift is approximately κ⁻¹b for ε ≳ 3τ, 1.8 times it at
  ε = τ; (H) without a leak the excursions into the saturation grow over 10⁶ slow units (the largest |α_1| by quarters 6.6 → 11.3,
  11–12% of the last quarter beyond the saturation).

Tests (`tests/test_gauge_drift.py`): the identities on flows with known limit maps and on an integrated one; the diffusion against
D‴.36's projection; the slow limit against the log-determinant law; the split-invariance of the log-determinant within a class.

## Proposition D‴.38 — a Fixman correction for the replicator under dither (`exp_gauge_fixman.py` → `chaos/gauge_fixman.md`)

- Proposition D‴.36's instance at τ = 10⁻⁴; the correction κ_F(τ/m)(1 + Ŵ_k²/W_k²)/(W_k − Ŵ_k²/W_k) on the replicator's fitness, Ŵ_k from the
  classes' contributions (the masses form) or from the current amplitudes (the plug-in form), equal on M. (A) At κ_F = 1 the corrected drift
  is 0.31, 0.07, 2·10⁻³ and 7·10⁻⁴ of the uncorrected at W = 1/3 for γ = 10⁻², 5·10⁻³, 10⁻³, 10⁻⁴ (O(γ²τ) against O(γτ)); from
  γ = 2·10⁻³ on it is reversed (the finite-γ projection weighs the correction's drift more than the dither's: 74% and 57% of their slow
  limits at γ = 10⁻²), at γ = 2·10⁻² by the factor 1.6.
- (B) At κ_F = 2 the corrected field's zero is the D-optimal fit W♯ = (0.3765, 0.2988, 0.3247) (the largest det ∇²_ΘJ; not the
  best-conditioned — the condition number is nearly flat along M) to 2·10⁻⁴ or better, attracting, the linearised spread there 0.040,
  0.014, 0.0046 at γ = 10⁻², 10⁻³, 10⁻⁴ (∝ √γ); at κ_F = 1 the remainder's zero is attracting at γ = 10⁻² but repelling at 10⁻³ and 10⁻⁴.
  (C) In Katzenberger's limit near the saturation the density's exponent is −0.99, −0.00, +0.98 for κ_F = 0, 1, 2 (predicted κ_F − 1);
  (H) across the implemented floor's cap (δ_cap = 0.00126) it falls back toward −1, as predicted to 10⁻³.
- (D, E) In the full loop (γ = 10⁻², 12 seeds, 10⁶ slow units; the plug-in form re-run with the implemented floor, the masses form with Â fixed) the correction keeps the amplitudes out of the saturation (the seeds' mean of the largest |α_1| 4–6 against 6.6 → 11.3 without, no
  record beyond a saturation; at κ_F = 1 excursions along M into a layer, 1% of the third class's last-quarter records with |α_3| > 10);
  at κ_F = 2 the second half's mean is within 1.07 s.e. of W♯ and of the mean path (plug-in) and within 1.0 (masses form); at κ_F = 1
  within 1.5 s.e. of both W♯ and the mean path, which the run cannot tell apart; equipartition with the rms fit error: p_eff = 8.96–9.04
  with the correction, 8.06 without (a saturated amplitude holds no fit energy), the dither's level read back to 0.4%.
- (F) With the leak (γ = 3·10⁻³, ε = 10⁻⁴, 16 seeds): two seeds' positions escape the exact fit's well in every form; the masses form keeps
  the masses where the leak holds them (0.31–0.35), the plug-in form moves them over 0.21–0.46, a running average of Â drifts with the
  escape; on the other 14 seeds the bias reverses as predicted (within 1.1 s.e. of the first order, 1.3 of the second).
- (G) Inside the leak's window (γ = 10⁻², ε = 3·10⁻⁴, 16 seeds): the shifts (−0.0013, +0.0011, +0.0002), (+0.0000, +0.0003, −0.0003),
  (+0.0012, −0.0005, −0.0007) at κ_F = 0, 1, 2 are within 1.03 s.e. of the second order — the drift balance plus the fluctuations'
  term through the leak's curvature, (−0.0004, +0.0002, +0.0002) at κ_F = 0 —, the paired differences within 0.93 s.e. of it (the first
  order 1.2–4.8 s.e. off): "the bias gone at κ_F = 1" is the finite-γ over-compensation cancelled by the fluctuation term. The same term
  accounts for Proposition D‴.37's Part G at ε = 10⁻⁴ and 10⁻³ (within 1.7 s.e.).

Tests (`tests/test_gauge_fixman.py`): the correction as (τ/2)∂_βℓ along M; no correction without dither; the slow-limit cancellation; the
gain's zero at W♯; the density exponents and the floor's cap; the general form reducing to the closed form; the masses form equal to the
plug-in form on M and blind to an amplitude's excursion at fixed masses; the tracked contributions fixed on M and relaxing at their time
constant; kappa() as the leak's projected stiffness; the analysis' statistics on synthetic records (`tests/test_runcache.py`: the cache
key carries the code's version).


## Proposition D‴.39 — dissipativity of the dithered gauge (`exp_gauge_plateau.py` → `chaos/gauge_plateau.md`)

- Proposition D‴.36's instance at τ = 10⁻⁴; the box L = (5, 3, 4) in the unconstrained coordinates (α, γ, κ), σ_B = 1. (A) The plateau of
  approximate fits with one amplitude saturated lies 7.6–7.8τ above the exact fit at W = 1/3 and at W♯ (fit error 0.035), 7.66τ at the
  box's edge: the box bounds the plateau, it does not raise it. (B) At W = 1/3 the minimum-energy path climbs to it monotonically.
- (C) At fixed masses the escapes occur at (1.33, 4.06, 10.9)·10⁻⁶ per slow unit for τ = 1.5, 2, 3·10⁻⁴ (10, 19, 24 of 24 seeds within
  4·10⁵), the Arrhenius slope (0.81 ± 0.14) of ΔJ = 7.69·10⁻⁴ at ΔJ/τ = 2.6–5.1 (weighted by the escape counts), within 1.3 standard errors of
  Freidlin–Wentzell's limit. (D) Without the box the plateau's occupation grows (0.18 → 0.57
  of the records by quarters, the largest |α| to 27.6); with it it stays at 0.09–0.19 and |α| ≤ 5. (E) In the full loop with the box
  (γ = 10⁻², 12 seeds, 10⁶ slow units) each run is the unboxed one until a coordinate leaves B; at κ_F = 0 the masses' spread at 10⁶ is
  0.05–0.075 against 0.08–0.15 without the box, on their way to P_B's vertices; at κ_F = 2 within 0.41 s.e. of W♯. (F) Without a leak
  all but one of the records with |α_k| > 10 have their mass within 0.01 of the saturation and the exact fit's error; with the leak
  every such record is off M (W_k more than 0.05 above the saturation, the fit error 0.041–0.043).

Tests (`tests/test_gauge_plateau.py`): the dead zone inactive inside the box, minus its G-gradient outside; J bounded on the unconstrained
coordinates and J + V coercive; the box's wall and vertices; a run inside the box bit for bit the unboxed one.


## Proposition D‴.40 — the gauge under general noise (`exp_gauge_noise.py` → `chaos/gauge_noise.md`)

- (A) The trace formula against the closed form at three points of M for eight noises (isotropic, the amplitudes only, the centres and widths only, levels (2, 1, ½)τ, three random correlated covariances, Part E's): the replicator flows agree to 1.3·10⁻⁸–1.0·10⁻⁶; with the centres and widths only the fitness is below 5.3·10⁻¹¹. (B) Katzenberger's limit drift at γ = 10⁻⁵ against the closed form's flow: within 1.1·10⁻⁴–8.8·10⁻⁴; the centres' and widths' noise alone: drift 8·10⁻¹⁶, diffusion 10⁻³⁶. (C) At finite γ the centres' and widths' dither alone moves the masses at O(γ²τ) (drift per γτ 3.2·10⁻⁴ … 8.0·10⁻⁷ at γ = 10⁻² … 10⁻⁵) with a diffusion below 10⁻²²: the gauge frozen.
- (D) The identity in simulation at fixed masses inside the box (2·10⁵ slow units; isotropic, amplitudes-only and correlated noise at τ = 10⁻³, 10⁻², the centres and widths alone at 10⁻²): the box-free form fails by up to 58 errors, its residual the box's term; the whole identity with the Euler chain's exact per-step bias (0.2–1.9% of the noise's term): χ² = 24.8 for 18 values (|z| ≤ 1.71; one seed for all runs, the diagonal-noise runs' twelve values three streams'), the boundary term below 0.03 errors; the centres and widths alone: the identity to 7·10⁻¹⁹; τ̂_k from the quadratic variation within 1.4% of τ_k.
- (E) The full loop inside the box (γ = 0.01, τ = 10⁻⁴, 12 seeds, 10⁶ slow units): the centres and widths dithered only — the masses within 2·10⁻⁴ of their start; with the barrier at κ_F = 1 — at W♯ within 8·10⁻⁴; the correlated noise at levels (2, 1, ½)τ without a correction — the masses drift (0.10 from W♯ over the second half); with the masses form at the classes' levels (κ_F = 2) — toward W♯_τ (the first mass 0.409, 0.428, 0.440 at 10⁵, 3·10⁵, 10⁶; the second half 0.022 from W♯_τ).

Tests (`tests/test_gauge_noise.py`): the noise mask; the correlated-noise path; the trace formula's reductions (isotropic, amplitudes
only, centres and widths only); only each amplitude's own variance at leading order under correlated noise; the levels and the per-class
correction; the pointwise gauge identity, the generator applied to ψ_k and the chain's closed-form expectation; the weighted analytic
centre.


## Proposition D‴.41 — measurement noise on the gauge (`exp_gauge_label.py` → `chaos/gauge_label.md`)

- Proposition D‴.36's instance without dither; the measured error r + ξ in both laws, ξ independent over the quadrature nodes, of intensity v.
  (A) At three points of M and γ = 10⁻², 10⁻⁵, Katzenberger's Itô drift (Σ = v⟨UUᵀ⟩_ρ) and the Wong–Zakai term (v/2)∂Φ⟨(U·∇)U⟩_ρ cancel to
  4.6·10⁻⁹–1.1·10⁻⁸ relative; the diffusion along M vanishes (≤ 2.1·10⁻²⁰); the closed form on the masses matches the Itô drift to
  2–4·10⁻⁴ at γ = 10⁻⁵. (A2) Other node covariances (an offset, an orbit-correlated noise, a variance profile): the Stratonovich drift over
  its Itô drift O(γ) from γ ≈ 10⁻³ on (from 10⁻⁴ for the 1 + 0.8 sin t profile) (2.4·10⁻⁵ at γ = 10⁻⁵ for an offset), 0.2–12% at γ = 10⁻¹; the coloured-noise identity to 2.5·10⁻⁸.
- (B) The slow limit's target of the per-update selection, W_ℓ = (0.3163, 0.4115, 0.2722), interior (|a|/a_max 0.79, 0.43, 0.73), far from
  the dither's D-optimal W♯ = (0.3765, 0.2988, 0.3247); at γ = 10⁻² the Itô drift vanishes at (0.3221, 0.3999, 0.2780); an Euler step in
  the amplitude a gives the factor (1 − q)(1 + 3q) (Katzenberger's drift in those coordinates within 3.4·10⁻⁴ of it) and the target
  (0.2967, 0.4513, 0.2520).
- (C) The full loop (γ = 10⁻², 4 seeds × 2·10⁵ slow units): the gauge's slope over the linear-response prediction averaged along the path
  0.996–1.008 for fresh draws and AR(1) draws at ρ = −1/2 … 0.6, 0.996–1.022 at ρ = 0.8, 0.9 (± up to 0.03), 1.004–1.027 at ρ = 0.95;
  the per-step law rejected at ρ = 0.9 (0.0462 ± 0.0006 of the Itô drift against 0.0526); Heun's scheme (2 ± 7)·10⁻⁴ of the Itô drift.
  (C2) At the step 8 (hλ_max ≈ 1.6): 0.979–1.049 of the prediction (± 0.004–0.053), where the step-free Itô drift is off by 42%.
- (D) f = uv: E[Q]'s decay rates within 0.3% of the linear-response prediction for every scheme (ρ up to 0.999, 0.19 of the per-step law
  there); Heun keeps Q to 1.2·10⁻⁵. (E) The diagonal network: the products' decay rates over the prediction 0.9988–0.9998 on average,
  the prediction differing across coordinates by a factor up to 2.17 at ρ = 0.999. (F) A tanh network: the sharpness's slope over the
  prediction 0.990–1.060 to ρ = 0.9 (slower noise below the measurement's floor).
- (G) The per-update selection's zero set at γ = 10⁻²: one interior target for fresh draws and τ_c = 10; three stable zeros, each with one
  class's mass large, from τ_c = 100; from τ_c = 10³ three saddles and an unstable node next to W♯ besides; W₀'s flow changing target
  between τ_c = 100 and 10³; at γ = 10⁻⁵ the targets approach the saturation (W/W_min − 1 ≈ 10⁻⁴ at τ_c = 10⁵). The slow limit τ_cΣ_f
  within 7·10⁻⁴ of its finite-step value; the slow noise's direction the dither's (cos 1.00000) at γ = 10⁻⁵; a smooth noise selects
  2–5% as much as a Markov one of the same correlation time.

Tests (`tests/test_gauge_label.py`): the noise along the fields; the fields orthogonal to the gauge in the flow's metric; the Wong–Zakai
cancellation at a point; W_ℓ interior and minimising Ψ_ℓ, and the amplitude coordinates' target; the schemes' reductions; the AR(1) draws
and Heun steps (correlated noise off M; dither) against hand-coded steps; the weighted and the covariance Wong–Zakai terms; the
linear-response factor against the exact stationary moments and its continuum limits; the filtered covariance against the pairwise
factors; Gaussian correlation sums; the finite-step slow limit.


## Proposition D‴.42 — several particles per class (`exp_gauge_classes.py` → `chaos/gauge_classes.md`)

- Proposition D‴.36's instance with 2 and 3 particles per class, together at the class's exact-fit position. (A) log det ∇²_ΘJ minus its
  one-per-class value at the same class masses is −46.131219 (m = 6) and −92.262438 (m = 9) at every class mass and split (to 3·10⁻⁷,
  4·10⁻⁷). (B) In Katzenberger's limit under the isotropic dither the class masses' drift and diffusion are 3/m times the one-per-class
  ones at γ = 10⁻⁵ (0.5000, 0.3333, n = 3's degenerate equal split included; the diffusion at γ = 10⁻² too); the fractions' diffusion is
  at most 10⁻¹⁶ per γ²τ, their drift 10⁻⁴ of the class drift at γ = 10⁻⁵ and at γ = 10⁻² O(γ) relative, toward the equal split.
- (C) The correction per class at κ_F = 1 cancels the class drift as for one particle per class, the fractions' drift unchanged.
  (D) Measurement noise at m = 6: at the equal split the Itô drift and the Wong–Zakai term cancel to 6·10⁻⁹ and the class Itô drift is m = 3's (1.0000); at unequal splits the residual is 2.2–2.6% of the Itô drift (up to 9.6% per component), at γ = 10⁻⁵ as at 10⁻²; the fractions' Itô drift is the full trace's flow less it.
- (E) The full loop at m = 6 (8 seeds × 4·10⁵ slow units): the class masses at (0.330, 0.349, 0.321) and (0.334, 0.315, 0.351) against m = 3 at half the time (0.340, 0.346, 0.314); the fractions still from the equal split, from (0.2, 0.7, 0.5) s₁ moving by −0.035 ± 0.012 — away from the equal split — as |α| reaches 26–28. (F) Anisotropic dither (Richardson-extrapolated): the per-class law at the equal split for every mask (1.0000–1.0004 at γ ≥ 10⁻⁵, 1.0075 at 10⁻⁶), the fractions still there; at the unequal split 1.4–3.5% off it for the amplitudes-only mask, 5.9·10⁴ times it for the centres-and-widths mask (which freezes one particle per class), the fractions drifting at 1–6·10⁻⁴ per γτ.

Tests (`tests/test_gauge_classes.py`): the determinant blind to the split; the class masses at the law of m and the fractions without
diffusion; the correction per class, off M and at the floor; the degenerate equal split.


## Proposition D‴.43 — the colour of the measurement noise selects the kernels (`exp_colour_selection.py` → `chaos/colour_selection.md`)

- (A) The identities: Cauchy–Binet to 2.6·10⁻¹⁵, the gradient (4/N)ω to 4·10⁻⁹, R's second differences negative along every random segment
  (−8·10⁻⁵ to −0.39, the continuum filters; zero for white draws), τ_cR − k log(4τ_c/N) → log(det(Φ_SᵀΦ_S)Π|w_i|) (−0.065892 against
  −0.065900 at τ_c = 10⁶).
- (B) Label-noise descent (10 × 40, k = 4 where basis pursuit fails, step 0.5, σ = 0.1, eight seeds × 16 paths) follows the reduced flow
  (the paths' mean within 0.011–0.035 of it for fresh draws, 0.11–0.35 for AR(1) draws at τ_c = 100); fresh draws end on dense fits
  (supports 10–13, errors 0.80–2.4, basis pursuit's 0.80–2.98); AR(1) draws recover the truth on all 16 paths of the four seeds whose
  reduced flow recovers it (the end interpolants' errors 10⁻¹³–7·10⁻¹¹) and on 0–2 of 16 on the other four; at σ = 0.3 the four recover on
  16, 0, 16 and 13 paths, at σ = 0.5 one seed still recovers — beyond the small-σ limit the descent prunes on its own.
- (C) Recovery (10 × 40, 40 instances per row): k = 3 basis pursuit 0.38, white 0.38, slow 0.62/0.75/0.78 (τ_c = 10/100/1000); k = 4 0.17,
  0.17, 0.42/0.60/0.62; k = 5 0, 0, 0.07/0.15/0.17; k = 6 none; at 20 × 80, k = 6 0.35 against 0.75; the white flow run to convergence ends
  at basis pursuit's solution (336 of 340 converged, to 9.4·10⁻⁹); from the white flow's end the slow flow does no better than basis pursuit
  at 10 × 40 (dropped atoms do not revive). Where the slow flow fails at τ_c = 100 its end has a larger R than the truth up to k = 4 (a
  local minimum; 0.09 and 0.35 of the failures end below the truth's R at k = 5, 6); at τ_c = 10, 0.39, 0.68 and 0.85 of them (k = 4, 5, 6)
  do — there the truth is not the finite-τ_c minimiser.
- (D) Gaussian kernel dictionaries (coherence 0.90, 0.96): the slow flow 0.95 against 0.80 (width 0.03, k = 4), 1.00 against 0.85 (width
  0.05, k = 4), 0.78 against 0.65 (k = 6), the white flow at basis pursuit's rate; label-noise descent with AR(1) draws recovers the truth
  on two of three instances where basis pursuit fails (all 16 paths, errors 2·10⁻¹³–7·10⁻⁸); on the third it is still on 12–13 atoms at the
  end, its reduced flow then at 8.8·10⁻³ and reaching the truth when run on (3·10⁻¹³ at s = 1600); fresh draws end within 0.08–0.11 of basis
  pursuit's solution.
- (E) At the step 0.5 a smooth noise of τ_c = 10 and a draw kept 40 steps select as well (0.70, 0.70 against fresh 0.30 and basis pursuit's
  0.45 at k = 3: the finite-step white comparator is not basis pursuit); a smooth noise of τ_c = 100 has not begun to select at the slow
  horizon (0.00) and recovers 0.85 and 0.65 at ten times it. (F) The flows run to convergence: from α = 0.1, 0.5 the slow flow recovers
  0.70, 0.75 (k = 3) and 0.40, 0.50 (k = 4), from α = 2 — toward the kernel regime — 0.35 and 0.30, basis pursuit's rate (0.45, 0.25), which
  the white flow matches from every α.
- (G) Sparse Bayesian learning (Tipping's EM) 0.97, 0.88, 0.55, 0.12, 0.05 for k = 2–6 against the slow flow's 0.97, 0.75, 0.60, 0.15, 0.00
  (mean supports at k = 4: 6.6, 6.4; basis pursuit 8.9); the log-det's MM 0.70, 0.45 (k = 3, 4) from the flows' start, 0.53, 0.33 from basis
  pursuit; reweighted ℓ¹ 0.38–0.53 and 0.17–0.25 at ε = 0.01, 0.1, 1 from either start.
- (H) The certificate (sufficient; for white draws Fuchs's, certifying 0.38 of the 2-sparse truths and 0–0.03 from k = 3 on) holds at every
  truth for τ_c = 100 at k = 2–5 (0.93 at k = 6) and for τ_c = 1000, the median ratio ∝ 1/τ_c (k = 2: 0.439, 0.0814, 0.0087 at τ_c = 10,
  100, 1000); on kernel dictionaries (widths 0.03–0.08, k = 4, 6) the certified fraction rises with τ_c to 0.85–1.00 at τ_c = 10⁵ (white:
  0–0.60), and the slow flow returns to the truth from P = 10⁻³ in 0.90–1.00 (τ_c = 100) and all (τ_c = 1000) of the instances.
- (I) Under a fit floor (y = Φw* + εe, |e| = 1, ε = 10⁻³–0.05) the slow flow's k largest atoms are the truth's support in 0.70–0.75 (k = 3)
  and 0.50 (k = 4) of the instances, basis pursuit's in 0.40–0.45 and 0.25. (J) At the same physical noise (AR(1) of correlation time 50,
  per-step σ = 0.1) Δlog P off the truth's support is −2.39 under Euler's step 0.5, −0.278 under Euler's step 0.125 and −0.097 under Heun's
  step 0.5 (25 times less): the selection is the explicit scheme's.

Tests (`tests/test_colour_selection.py`): the subset expansion and the gradient, also with a finite step's filter; the white and slow
potentials (finely sampled AR(1) draws tending to the slow one), R linear (white) and concave (slow) along segments; the support count; the
projection, singular cases included; the reduced flow on M, P's rate, finite weights at a sparse fit, an overshooting trial stage rejected,
white → basis pursuit, slow → the sparse truth; the flow as R's gradient flow in its metric (a doubled or halved w-drift fails); the exact
interpolant against small-step descent, the end interpolants (at infinity flagged, converging within the alive atoms' span); label-noise
descent against hand-coded AR(1) steps; the dictionary; the (u, v) round trip; the experiment's success threshold and dual certificate
against hand computations; sparse Bayesian learning's EM on a 2-sparse truth; the white flow run to convergence; Heun's noiseless step on M.

## Proposition D‴.45 — the free mass–amplitude law as a mirror flow (`exp_mirror_free.py` → `chaos/mirror_free.md`)

- (A) Unless stated: 5 × 12 Gaussian designs with unit columns, 2- or 3-sparse truths, m = 12, G = 1, γ = 0.1, a_max = 1 (c = 0.025); no
  projection failed. η's closed form agrees with the quadrature of dα/(Gβa′) to 1.5·10⁻¹⁴ (relative, free curves and a branch), (ψ*)′ = w
  and dw/dη = mD to 7.2·10⁻¹¹ and 4.3·10⁻¹⁰. Along the loop Q moves by at most 1.4·10⁻¹¹ noiseless and 3.6·10⁻¹⁰ under a smooth noise of 30
  sinusoids mixed over the nodes (the integrator's tolerance); the conserved law's shift equals the quadrature of −γβ_iλ to 1.5–4.4·10⁻⁷;
  along D‴.44's reduced flow dQ/ds is (iv)(b)'s closed form to 9.1·10⁻¹⁶.
- (B) From α = 0 and equal masses s·m/d the noiseless loop's end at t = 10⁷ is the Bregman projection to 4.7·10⁻¹¹ (20 instances per k). For
  k = 2 the projection is 0.5% from the minimum-norm fit at s = 10 and its ‖w‖₁ 1.0004 times basis pursuit's at s = 10⁻⁴, the truth
  recovered nearly as often as by basis pursuit (0.75 against 0.80); at s = 10⁻³ the loop is within 10⁻⁸ of its projection only from t ≈
  8.5·10⁶. Unequal masses s·r_i: weighted basis pursuit to 1.7·10⁻⁹ and 2.1·10⁻⁸ (k = 2, 3) at s = 10⁻⁷, the weighted minimum norm to
  2.4–2.5·10⁻⁹ at s = 10⁴; from small random amplitudes at s = 10⁻⁶ the particles sit on branches, the bias the sign-constrained one
  (D‴.47). (F) At N-atom vertices the leading order from the certificate agrees with the projection to the rounding floor (1.8·10⁻¹³,
  2.6·10⁻¹¹ at ε = 10⁻⁷) and the distance falls by 1068 and 984 from 10⁻⁷ to 10⁻⁹ (ε^{3/2}: 1000); a small slack (μ_min = 2.7·10⁻⁴) keeps a
  sixth atom, 0.17 away at 10⁻⁷.
- (C) White label noise (σ = 0.3, 64 paths): Euler–Maruyama changes Q on average by its Itô correction's integral (ratios 1.0106, 1.0037,
  1.0011 at h = 0.1, 0.03, 0.01), Heun's (Stratonovich) mean change is 1.4·10⁻⁵, 4.1·10⁻⁶, 1.4·10⁻⁶, the mirror update keeps Q to 2.2·10⁻¹⁶;
  the Euler step's exact mean change is the closed form's to 0.029–0.036% at h = 0.005. Along D‴.44's reduced flow from the lazy projection
  ‖w‖₁ first rises, then Q's mean falls to −0.021 by its time s = 10⁸ (no atom yet extinct): the state ends on basis pursuit where it is the
  truth (6 instances), elsewhere on its support (two) or on another sparse fit of larger ℓ¹ norm (two). At h = 0.5 the paradigm's Euler
  update lowers Q by 0.34–0.36 and moves the projection by 0.030–0.043, as the reduced flow does; the mirror update keeps it to 9.0·10⁻¹².
- (D) The real loop (k_e = 5, m = 6, six Gaussian atoms, a second particle on the third, 2000 time units): noiseless, V falls monotonically
  from 0.175 to 0.022 with V̇ = −k_ee² to 3.0·10⁻⁹, Q kept to 7.6·10⁻¹⁵ and the pair's mirror coordinates moving together to 2.4·10⁻¹⁵;
  under a smooth noise in the adaptation, or a sensor noise, V is no longer monotone (0.175 → 0.136, 0.377) while Q is kept to 6.8·10⁻¹² and
  4.9·10⁻¹³; at t = 2000 the noiseless coefficients are still 0.048 (the distinct atoms) and 0.027 (the pair's sum) from w*. (E) The
  conserved law (Σβ = m) from the same starts reaches exact fits 1.6–4.1·10⁻³ from the free law's projections, its Q moving by 0.058–0.22.

Tests (`tests/test_mirror_free.py`): η against quadrature in the three cases, on both sides of a branch and on the critical curve, and
against a 60-digit reference near κ = 1 and |t| = 1; the inverse's round trips and its range errors; the conjugate (ψ*)′ = w and the
metric dw/dη = mD; Q under a smooth node-correlated noise and the conserved law's shift dQ = −γβλ dt; the noiseless limit equal to the
Bregman projection; the rich projection converging to weighted basis pursuit; the Euler step's second-order change of Q; V̇ = −k_ee² in
the real loop (and a wrong mass law seen by V); the mirror update keeping Q and moving η by hg, on branches too.


## Proposition D‴.46 — a leak in the mirror coordinates (`exp_mirror_leak.py` → `chaos/mirror_leak.md`)

- (A) The identity of (i) — the leak with a disturbance and a measurement noise — holds to 9.5·10⁻¹¹ (relative) at twelve random states of
  the real loop (σ = 0 and 0.3). On random positive definite H (10 × 10, k_e = 2, four draws) the averaged equilibrium is the ridge with ρ =
  k_eσm/(Gβ⁰²) = 0.5, which moves w by 31–55%, to 1.3·10⁻⁸ (β⁰ = 10⁴) and 1.3·10⁻⁴ (10²), but departs from it by 0.15 at β⁰ = 3, where the
  saturation's wall w_⊓ = 0.25 lies inside the coefficients' range: a ridge only below the wall; at small masses it is the LASSO of penalty
  σλ₁ = 0.02 to 0.015, 1.5·10⁻³, 1.5·10⁻⁴ (β⁰ = 10⁻⁴, 10⁻⁶, 10⁻⁸; the remainder ∝ β⁰^{1/2}), with the LASSO's support in every draw.
- (B) On a persistently exciting orbit (the filtered Gram's eigenvalues 2.7·10⁻³ to 0.68) under an approximation error and a biased
  disturbance the free law without the leak creeps toward the biased fit w* + H_L⁻¹d̄_L (its largest coefficient 0.309, 0.330, 0.337, 0.341
  by quarters, V rising from 0.40 to 0.45); σ = 0.01 and 0.1 hold the coefficients at 0.081 and 0.040 and V at 0.085 and 0.13, below C/λ =
  1.74 and 0.30 and the comparison bound at every sample, the late means within 2.5·10⁻³ and 6.5·10⁻⁴ of w_σ^L, the rms error rising from
  0.024 to 0.031 and 0.037. Without persistent excitation the free law's coefficient grows as constructed, 0.70, 2.31, 7.37, 23.3 at t = 1,
  10², 10⁴, 10⁶, while σ = 0.01 and 0.1 bound it (0.46 at t = 10⁶).
- (C) Slow adaptation (k_e = 20): the error filter's Gram is 0.89% from the idealised ⟨φφᵀ⟩, symmetric to 3·10⁻¹⁴; started at the averaged
  equilibrium w_σ^L, the real loop's last period means lie within 1.6·10⁻⁴ (lazy start) and 1.7·10⁻⁴ (intermediate) of it, 4.4 and 2.6 times
  closer than to the idealised equilibrium (the leak's bias 0.038 and 0.12). (D) On an orbit's Gram the rich start's equilibrium is the
  LASSO to 3.7·10⁻⁶–6.5·10⁻⁴, keeping exactly the truth's three atoms for μ ≤ 10⁻² (one dropped at 3·10⁻²), the lazy start's the ridge to
  3.9·10⁻⁶–2.5·10⁻⁴, dense; over all 1760 3-sparse patterns the support is exactly the truth's for 415, the LASSO's prediction holds for
  416, the two agree on 1745; with a correlated disturbance bias the equilibrium tends to w* + H⁻¹d̄ as σ → 0, linearly; from the rich start
  it is reached at t ≈ 3.2·10⁵, at 2·10⁶ where the irrepresentable condition fails (five atoms kept).

Tests (`tests/test_mirror_free.py`): the Lyapunov identity with the leak, a disturbance and a measurement noise; Q kept by the leak; the
equilibrium's ridge limit, its biased limit w* + H⁻¹d̄ and its LASSO limit (same support) at k_e ≠ 1, the solver's residual asserted; the
comparison bound under a worst-case disturbance, coming within a factor 2 of its ultimate value; the filtered Gram's expansion
H − ⟨φ̇φ̇ᵀ⟩/k_e² + O(k_e⁻³) and its O(k_e⁻³) asymmetry.


## Proposition D‴.47 — the free law's branches (`exp_mirror_branch.py` → `chaos/mirror_branch.md`)

- (A) The slope dη/dlog|w| = a/(γa² + Gβa′²) holds to 2.3·10⁻⁹ (40 random states, free and on branches); at the branch point the temperature
  is 1/(γa₀) to ten digits (α₀ = 0.05, 0.3, 1, 2); B(w̄, w) → β/(γm) as w̄ → 0 (2·10⁻⁴ relative at δ̄ = 10⁻⁶, 5·10⁻¹⁰ at 10⁻¹²); with α(0)
  uniform on [−0.3, 0.3] the fraction of starts on branches is 0.52, 0.95, 0.995, 0.9996 at β(0) = 10⁻³, 10⁻⁵, 10⁻⁷, 10⁻⁹; a ± pair's
  coefficient follows 2u sinh(γa₀Δ) to 2·10⁻⁵ for Δ ≤ 10 (6·10⁻⁴ at Δ = 120, where the members leave their branch points); on the critical
  curve η·w^{1/3} → −K = −3.218 (−3.2182 at α = 10⁻⁴).
- (B) Random 6 × 12 Gaussian designs (the target on four atoms, masses 10⁻³, amplitudes 0.2–0.4, 8 seeds for each of three sign patterns):
  the flow at t = 3·10⁷ is the Bregman projection onto S, computed independently, to 5.0·10⁻¹⁰; where the frozen signs cannot fit y (5 of 8
  positive starts, 2 of 8 random ones) J reaches its sign-constrained minimum and the wrong-signed particles go extinct (masses below
  10⁻¹⁰³) at the rates γa₀|g(∞)| to 0.5–10%, B(w̄, w(t)) never increasing by more than 2.2·10⁻¹⁶; the mixed starts fit exactly. With
  duplicated atoms (S not a point; 8 seeds) the flow is the Bregman projection to 6.4·10⁻¹⁰, 1.6·10⁻⁴–0.20 from the equal split of the
  duplicates: the selection with extinction is the Bregman one, not the symmetric one.
- (C) On one 8 × 20 design: from branches (amplitudes 0.2–0.6, masses ε) the projection's distance to the sign-constrained basis pursuit
  weighted by 1/|a(0)| is 0.23, 0.030, 1.3·10⁻³, 5.1·10⁻⁶, 1.1·10⁻⁹ at ε = 10⁻⁴, 10⁻⁸, 10⁻¹⁶, 10⁻³², 10⁻⁶⁴, the leading order agreeing to
  4.1·10⁻⁴, 2.4·10⁻⁹, 3.4·10⁻¹⁴ from 10⁻¹⁶ on (μ_min = 0.104): the approach ε^{μ_j}; from α(0) = 0 the distance to basis pursuit is 0.29,
  0.14, 4.0·10⁻⁴, 5.2·10⁻⁸, 4.3·10⁻¹¹ at ε = 10⁻², …, 10⁻¹⁰: ε^{3/2}. In a mixed start the free particles carry ≈ ε² — or, where the branch
  atoms' frozen signs cannot fit, 2.27 on the power scale, the projection failing at ε = 10⁻¹² (residual 1.95). (D) A ± pair per atom tends
  to basis pursuit weighted by 1/a₀, slowly (slacks 0.016, 0.020): 0.29 at u = 10⁻⁴, 2.6·10⁻⁴ at 10⁻¹²⁸, 1.9·10⁻⁶ at 10⁻²⁵⁶; never are both
  members of an atom active.
- (E) Positive Gaussian dictionaries, positive branches: with 21 atoms of width 0.15 and 7 samples (margin μ = 0.015) the loop reaches the
  target from every start scale, β(0) = 10⁻⁶ to 10, within 6.0·10⁻⁵ at t = 10⁹, where from α(0) = 0 the rich start ends on a seven-atom fit
  1 away and the lazy one spreads over all 21 atoms; at μ = 0.003 it is still 1.0–1.3·10⁻³ away at 10⁹, at μ = 0.027 on the target's atoms;
  the decay is algebraic (t^{−0.9}–t^{−1}). Over all three-atom supports the median margins are 4.0·10⁻⁶, 3.8·10⁻⁹, 1.2·10⁻³ (25%, 45%, 4%
  at most 10⁻⁹); the loop reaches the supports at the quantiles whose margins exceed 10⁻³ by t = 10⁹, and on those below 10⁻⁴ is still
  0.37–0.83 away: the selection's speed is the margin. (F) The real loop on positive branches: V̇ = −k_ee² to 1.7·10⁻⁸, V falls
  monotonically from 3.64 while the unused atoms' masses rise to 0.78; after 4000 time units the coefficients are still 0.080 from w*.

Tests (`tests/test_mirror_free.py`): the branch criterion and δ(0) solved from β(0) (masses far below c's precision kept), η's inverse to
any depth, the temperature 1/(γa₀) at a branch point, B(w̄, w) → β/(γm) as w̄ → 0, a ± pair's W(Δ) = 2u sinh(γa₀Δ) at leading order;
the limit theorem: from positive branches that cannot fit the labels the loop converges to the projection onto the sign-constrained
fits (computed independently), the vanishing coordinates' masses below 10⁻³⁰.


## Proposition D‴.48 — the rich transit (`exp_mirror_transit.py` → `chaos/mirror_transit.md`)

- (A) From α = 0 on three 6 × 14 Gaussian designs (3-sparse truths, masses εb_i, b_i log-uniform in [1/4, 4]) the adaptive inverse scale
  space path of Σ|w_i|/√b_i ends at weighted basis pursuit after 3, 8 and 10 events, and the loop in τ = t/Λ(ε) follows it: sup_τ|u − u_ISS|
  falls from 0.54, 0.61, 1.15 at ε = 10⁻² to 6.7·10⁻⁶, 8.8·10⁻⁶, 2.0·10⁻⁵ at 10⁻¹² (log-log slopes 0.500 between 10⁻⁶ and 10⁻⁸: 1/Λ ∝
  ε^{1/2}); the entry times converge at that rate down to the grid's step; the residual's difference is O(1) at an event over a τ-width ∝
  ε^{1/2} and O(ε^{3/2}) away from the events (median slopes 1.51, 1.69, 1.46); the end at 1.5τ* is weighted basis pursuit to
  1.1·10⁻¹⁶–6.1·10⁻¹⁵ at 10⁻¹².
- (B) On branches (amplitudes 0.2–0.6) the loop follows the one-sided path in τ = t/log(1/ε): sup|u − u_ISS| = 15.7, 8.0, 4.0, 2.0, 1.0 and
  the entry times' largest difference 55, 28, 14, 7.0, 3.5 at ε = 10⁻⁴, 10⁻⁸, 10⁻¹⁶, 10⁻³², 10⁻⁶⁴, exactly ∝ 1/log(1/ε) (the path's weights
  1/(γ|a(0)|) between 17 and 50, the amplitudes' range), ending at the sign-constrained weighted basis pursuit (to 1.3·10⁻⁵ at ε = 10⁻⁴, below 10⁻¹⁴ from 10⁻¹⁶ on); ± pairs
  follow the two-sided path the same way (21, 7.9, 4.0, 2.0, 0.99 at u = 10⁻⁴ … 10⁻⁶⁴).
- (C) Moving centres: the sliding path enters three kernels that slide onto the three bumps and fit them exactly; the loop with the
  paradigm's metric forms 8 kernels at ε = 10⁻³, 4 at 10⁻⁴ and exactly these three for ε ≤ 10⁻⁵, the others displaced ∝ ε; with the
  Wasserstein–Fisher–Rao metric the centres follow (iv)'s limit and the loop forms 15, 12, 8, 3, 3 kernels at ε = 10⁻⁴, 10⁻⁵, 10⁻⁶, 10⁻⁸,
  10⁻¹⁰; with the centres held on the grid it forms 5 kernels, 0.019 from the target. On sixteen random targets the loop forms exactly three
  kernels at the bumps on 11, 13, 16 at ε = 10⁻⁵, 10⁻⁶, 10⁻⁸ (Wasserstein–Fisher–Rao: none at 10⁻⁶, all at 10⁻¹⁰); a tie is shared (0.4,
  0.4), a near tie takes two kernels at ε = 10⁻⁶; the creep's speed is the profile's to 0.990–0.996 (ε = 10⁻⁸) and 0.9999–1.0000 (10⁻¹⁰);
  wider bumps are fitted along a staircase. (D) Moving widths fit bumps of widths 0.2 and 0.08 with exactly three kernels at ε = 10⁻⁶, where
  fixed widths need six kernels or miss the amplitudes (0.71, −0.48, 0.59).

Tests (`tests/test_rich_transit.py`, `tests/test_mirror_free.py`): the critical curve's coefficient and its depth-3 form; the sliding
path recovers three bumps exactly; the loop at ε = 10⁻⁵ forms the same three kernels at the path's entry times; the inactive kernels'
coefficients below the core's bound before the first entry; the Wasserstein–Fisher–Rao metric's centres on their limit; the inverse
scale space path ends at weighted basis pursuit (one-sided: the sign-constrained one) and the loop at ε = 10⁻⁶ follows it.


## Proposition D‴.51 — the transit at a finite mass (`exp_transit_mass.py` → `chaos/transit_mass.md`)

- (A) One bump between the kernels at 0.0 and 0.1, δ from their midpoint (A = 0.3, 0.8, 2; δ = 0.0001, 0.001, 0.01; ε = 10⁻⁵ … 10⁻¹⁰;
  the loop to 1.32 τ_entry): the runner-up's share over the ε-free inner problem's at the same x = μ*/ε^{1/2} is 1.0000–1.0064 at
  δ = 0.0001, 1.0001–1.0668 at 0.001 and 1.0002–1.5448 at 0.01 (x from 0.11 to 3413; at δ = 0.01 1.3858–1.5448 at x = 10.79, the
  partners' approach creep at the actual geometry), the inner problem stationary on 51 of 54 and the other 3 (x ≈ 34.1–34.7) drifting
  by at most 9.8·10⁻⁴ of their share from s = 4.69·10⁶ to 1.88·10⁷. Whatever the amplitude the shares are 0.48, 0.437, 0.309, 0.073 at
  x = 0.11, 0.35, 1.1, 3.5 against the Gibbs split 1/(1 + e^{κx}) (κ = 0.7025) 0.481, 0.439, 0.316, 0.080; far behind, the runner-up
  keeps the frozen coefficient K'³ρ⁻³/(|A|(x − Δ)³) (K'³ = 0.0962), 3.04·10⁻¹² at x = 3413 (A = 0.8). It slides onto the bump with the
  entrant while x ≲ 11 and stays within 0.011 of its start from x ≈ 34.
- (B) The overshoot η_c(α_A)/(Λ₀ρ) is 3.6500, 4.6370, 5.2143, 5.9415, 6.6105, 7.1877, 7.9149, 8.9016 for A = 0.1, 0.2, 0.3, 0.5, 0.8, 1.2,
  2, 4, growing by 1.4240, 1.4237, 1.4236, 1.4235, 1.4235, 1.4235, 1.4235 per unit log A between consecutive amplitudes — the exponential
  branch's slope a/(4GcΛ₀ρ) = 1.4235 (Λ₀ = 7.0248); the inner profile against its core power law, 𝒲(ξ)/(K'³ξ⁻³) − 1 = −6.7%, −14.0%,
  −19.1%, −39.5%, −65.2% at ξ = 15, 10, 8.3, 5, 3. The matched pull on a frozen runner-up at the exact tie (x = 200), from the starts
  xi0 = 400, 800, 1600, 3200 behind the wall (the entrant's far-field creep: a/xi0 + b/xi0²) extrapolated: 2.2300, 3.1749, 3.7275,
  4.4238, 5.0644, 5.6171, 6.3134, 7.2582, so I∞ = 1.42–1.64, growing slowly with A; at δ = 0.001 2.2776 … 7.2970. The loop's (δ = 0.001;
  ε = 10⁻⁸, 10⁻⁹, 10⁻¹⁰, x = 34.66, 109.6, 346.6) with the runner-up's centre held fixed is the inner problem's to +0.0024 to +0.0028,
  +0.0008, −0.0007; free, it exceeds that by the runner-up's approach creep, 0.3753–0.4381, 0.1170–0.1228, 0.0366–0.0371, against
  (a²b/(π²c))(g′_r/g_r)²/x = 0.3667, 0.1159, 0.0367 (ratios 1.02–1.20, 1.01–1.06, 1.00–1.01: the excess falls tenfold as x triples, a
  1/x² term growing with A).
- (C) The creep on target 10 of D‴.48's random targets (centres −0.752, 0.112, 0.488, amplitudes 0.868, −0.412, −0.514; the path's
  kernels 2, 14, 12): predicted from the path, the entry threshold 6.12·10⁻⁷ (kernel 3 at the entry of 2, margin 0.00489) and the creep
  thresholds 4.88·10⁻⁹, 1.36·10⁻⁶, 5.56·10⁻⁵, 2.62·10⁻⁶ (kernels 15, 3, 11, 13): ε_c = 4.88·10⁻⁹, kernel 15's creep. The loop forms the
  path's kernels at ε/ε_c = 0.205, 0.41, 0.717 and enters kernel 15 in place of 12 at 1.43, 2.87, 6.14 (the inactive kernels' largest
  displacement 0.0884, 0.0884, 0.138, then 0.288).
- (D) 48 random three-bump targets (seed 51; the path exact — three kernels, rms residual < 10⁻⁸ — on all 48): the loop recovers 32, 38,
  42, 47 of them at ε = 10⁻⁵, 10⁻⁶, 10⁻⁷, 10⁻⁸, the estimate predicts 27, 36, 43, 47 and agrees on 43, 46, 47, 48 (5 and 2 recovered
  though not predicted at 10⁻⁵ and 10⁻⁶, 1 not recovered though predicted at 10⁻⁷) — a grid a decade apart, which does not measure the
  estimate's accuracy (G).
- (E) 400 random targets (seed 52): the fraction predicted not to recover is 0.7700, 0.4650, 0.2400, 0.1225, 0.0700, 0.0300, 0.0150,
  0.0075 at ε = 10⁻⁴, 10⁻⁵, …, 10⁻¹¹, by a near tie 0.6950, 0.2775, 0.1050, 0.0300, 0.0125, 0.0025 and by the creep 0.7025, 0.4250,
  0.2125, 0.1175, 0.0650, 0.0300 at ε = 10⁻⁴ … 10⁻⁹ (not exclusive); the least-squares slopes in log ε are 0.294 (all, ε = 10⁻⁵ … 10⁻⁹:
  ≈ ε^{1/4}), 0.462 (a near tie, 3.16·10⁻⁵ … 10⁻⁷: ≈ ε^{1/2}) and 0.285 (the creep, 10⁻⁵ … 10⁻⁹). These are the estimate's fractions;
  the loop's own were not measured.
- (F) Every target of E with 10⁻¹² ≤ ε_c ≤ 10⁻⁷ (48; 1 more below 10⁻¹², not run; 7 bound by an entry, 41 by the creep): at ε_c/3 all
  48 recover, as predicted; at 3ε_c 43 fail and five creep-bound targets recover (8, 115, 117, 216, 267) — their loop thresholds,
  bracketed and bisected as in G, at 9.11, 1.15·10⁵, 2.33·10³, 9.93 and 225 times ε_c: false alarms of the creep criterion.
- (G) The loop's own threshold on D's 48 targets, bracketed from ε_c by factors 4 and bisected 4 times in log ε (the bracket ×1.09):
  ε_loop/ε_c has the median 0.771 (quartiles 0.545, 0.917) and the range 0.131–94.5; it is below 1 on 41, within ×1.5 on 20 and within
  ×2 on 36. Where an entry binds (32 targets) 0.648 (0.545, 0.808), the range 0.131–43.3 — by the predicted mass 0.917 for ε_c ≤ 2·10⁻⁶
  (5 targets), 0.771 between 2·10⁻⁶ and 10⁻⁴ (17), 0.545 above (10): the bias grows with the binding margin —; where the creep binds
  (16) 0.844 (0.545, 8.24), the range 0.324–94.5. The outliers: target 2 at 0.131 (its first entry's runner-up migrating to an unfitted
  bump, H(3)); targets 15, 28, 45, 46 at 94.5, 25.8, 16.7, 19.9, the creep criterion's false alarms; target 27 at 43.3, a pull measured
  on a block relaxation slower than the time to the next entry.
- (H) Near-tie controls. (1) The single bump (A = 0.8) at δ = 0.001, 0.01, 0.03 from the midpoint (μ* = 0.0035, 0.0341, 0.0989): the
  runner-up's formation threshold x_θ = 6.392, 7.100, 8.189 free, 5.315, 5.819, 6.373 with the runner-up's centre held fixed until the
  entrant's arrival and 6.557, 6.546, 6.520 with both centres; the inner problem's own at the exact tie 6.404, 6.300, 6.272 for A = 0.3,
  0.8, 2, below the Gibbs value ln(|A|/θ − 1)/κ = 6.541; the free loop's rise over it, 0.091, 0.800, 1.889, is 0.597, 0.581, 0.528 times
  the partners' approach-creep difference (a²b/(π²c))((g′_r/g_r)² − (g′_e/g_e)²)/x_θ = 0.153, 1.376, 3.580; the estimate's own x for
  this runner-up is 6.252, 6.145, 5.853 (x_θκ_ij with κ_ij = 0.956, 0.939, 0.895, over the frozen pull's Δ + ξ_θ = 5.146, 5.443, 5.722).
  (2) Exact ties made by unequal masses: the farther kernel's share 0.8849 and 0.8854 at ε = 10⁻⁸ and 10⁻¹⁰ (b = 1.2000, the bump at
  0.07625; both kernels formed, as on the path) and 0.9972, 0.9973 (b = 1.6260, the bump at 0.08; kernel 12 alone, where the path forms
  11 and 12). (3) D's target 2 with its middle bump moved to c₂ = 0.043, 0.046, 0.048: the first entry's runner-up (kernel 11, μ* =
  0.01890, 0.00866, 0.00177) forms and ends on the unfitted bump at 0.632, the loop's threshold 1.32·10⁻⁶, 9.55·10⁻⁸, 3.27·10⁻¹⁰
  (ε_loop/μ*² = 0.00371, 0.00127, 0.000104; slope 3.52 in log–log); at c₂ = 0.03, 0.035, 0.04 the extra kernel is another (4, at
  −0.655). (4) E's target 117: kernel 14, between the third entrant (13) and its bump, reaches at most 10⁻⁷, 4.15·10⁻⁶, 6.53·10⁻⁴
  (θ = 0.00896) at x = 5.186, 1.640, 0.519 as the entrant slides through it (their closest approach 0.0053, 0.0716, 0.0000), where the
  Gibbs split would give it 0.026, 0.240, 0.410; the loop forms kernels 12, 13, 2, the path's.

Tests (`tests/test_transit_mass.py`): the inner profile is the small-mass limit of the mirror map (not the core's power law), with the
next orders 1 − 8(κξ)⁻² + 64(κξ)⁻⁴ and a dependence on ρξ alone; the overshoot's slope on the exponential branch and the Gibbs ratio of
two duals; the inner problem's tie (½), its Gibbs split, the frozen tail and its sign symmetry, and the loop's share against it; the
inner threshold below the Gibbs value; the partners' approach creep raising the loop's threshold (free against both centres fixed);
unequal masses breaking the tie; a near-side runner-up passed over; the path's duals reaching their walls at the entries; the predicted
mass on D‴.48's lost near tie (the creep's, kernel 15) and the creep's monotonicity in ε; `entry_pull` as the frozen runner-up's matched
pull, with its 1/xi0 tail; a rearranging entry's measured pull against the near-tie estimate (E's target 327); the estimate overstating
the loop's threshold at a large margin (D's target 7).


## Theorem D‴.53 — the sliding limit (`exp_transit_limit.py` → `chaos/transit_limit.md`)

- (B) D‴.51's 48 out-of-sample three-bump targets (seed 51), equal masses ε: 12 have an exact path with simple entries, no exits and
  margins ≥ 0.05. Over the 48 the approach speeds v_k lie in 0.00132–0.125, the slides in 0.000176–0.302, |s_θ + σ_k| in
  2.698–9.856·10⁴ and the largest pull per entry in 1.17–3704 — the extreme values on the slowest entrants (v_k ≤ 0.008, six excluded
  targets) —, λ_k in 0.0215–0.0883; s_θ and the pulls by second-order Richardson in xi0 from 1600, 3200, 6400 (`entry_pull`'s deviation
  form, Radau at rtol 10⁻¹⁰), their spread against the same from 800, 1600, 3200 at most 2.8·10⁻⁵ in s_θ and 1.3·10⁻⁶ in the pulls on
  the twelve (up to 13 and 0.52 on the slowest entrants); (H3) holds at every entry, the inner problem's end within 3.1·10⁻¹¹
  (coefficients) and 6.4·10⁻¹³ (centres) of the path's S_k, every inner problem stationary.
- (A) The loop at ε = 10⁻⁶ … 10⁻¹¹ (rtol 10⁻¹¹ for ε ≥ 10⁻⁸, 10⁻¹² at 10⁻⁹ and 10⁻¹³ for ε ≤ 10⁻¹⁰) on eight of the twelve (targets 39,
  24, 19, 6, 41, 10, 33, 4): it forms the path's three kernels throughout; the never-entering kernels' largest coefficient is
  28.2–391·ε^{3/2} from 10⁻¹⁰ on (31.5–527·ε^{3/2} at 10⁻⁶; slopes 1.501–1.504 in log ε) and their centres move by the creep's εC_j to
  within 0.08–0.20 at 10⁻⁶ and 5·10⁻⁴–1.5·10⁻³ at 10⁻¹¹ (slopes 1.001–1.002); the block is within 1.4·10⁻¹⁵–5.6·10⁻¹³ of S_k at the
  midpoints at 10⁻¹¹ (slopes 1.23–1.53) and within 4·10⁻¹⁰–3·10⁻⁷ at (b)'s window's end; the duals' deviation is the first order — the
  median obs/pred 0.998–1.0005 at 10⁻⁶ and 1.0000 from 10⁻⁹ on —, the largest |obs − pred| (in units of ε^{1/2}) falling by 1.6–5.5 per
  decade below 10⁻⁷, most near √10 (a remainder O(ε) in u), to 2.9·10⁻⁵–6.2·10⁻⁴ at 10⁻¹¹, while the pulls alone leave 0.0019–0.39:
  the entries' shifts are the rest of the first order; the coefficients cross θ at s_θ + σ_k to within 3.6·10⁻⁴–9.5·10⁻³ at 10⁻¹¹. At
  10⁻⁶ targets 6 and 4 lie outside the regime — their blocks 0.186 and 0.481 from S_k at a midpoint, an entry within ε^{1/2}|s_θ + σ|
  of the next (gaps 0.219 and 0.142 in τ) —, inside it from 10⁻⁷ on.

Tests (`tests/test_transit_limit.py`): an isolated entrant's inner problem descends the fit (E never increases, Lemma D‴.53a) and ends
on its feature; on target 39 at ε = 10⁻¹⁰ the loop forms exactly the path's kernels, the never-entering kernels' coefficients are the
core's profile at their path margins, their displacements the creep's εC_j, their duals off the path by ε^{1/2}·(the entries' pulls and
shifts) + εD_j — the shifts' term needed —, and the entries cross at τ_k + ε^{1/2}(s_θ + σ_k); on target 41 at 10⁻⁹ (an entrant that
slides, one that arrives late) the sign per entry and the shifts, with converged predictions, (H3) at every entry and the block at S_k.

## Theorem D‴.55 — exact recovery by the sliding transit (`exp_transit_certificate.py` → `chaos/transit_certificate.md`)

- The sliding inverse scale space `sliding_iss_nd` (in one dimension on the 1-D samples, with the exact Jacobian; a fast-flow chunk that
  LSODA has not finished within 2·10⁵ evaluations redone by Radau), its least-squares polish bounded; three features of the kernels'
  width; E the starting kernels of largest correlation; exact with E: the active set E, the residual below 10⁻⁸, no exits; the
  skeletons recomputed from the targets at report time.
- (A) One dimension (21 kernels of width 0.12 at spacing 0.1; 400 targets, minimal separations 0.25 … 0.6, amplitudes ±U(0.4, 1)):
  every path an exact fit with three kernels, exact with E on 262; exactness with E predicted by the fitted skeleton on 400, the linear
  skeleton on 397 (its misses at m_sk = −0.0036, −0.0095 and −0.188, where the fitted skeleton's margins are +0.006 … +0.008), the atoms
  subtracted on 390, the relative margin on 344, the certificate on 331 (its disagreements at |m| ≤ 0.195), the first entrant on 327
  (65 decided at the first entry), always E on 262; with the polish off the flow ends as the polished path on 400 of 400.
- (B) The amplitudes' range (40 configurations, |A_k| = 10^{−r·u_k}, r = 0 … 3, ranges up to 749): the outcome changes with the range on
  9 (all at |m|/q ≤ 0.67); both skeletons right on 160 of 160, the atoms subtracted on 140, the first entrant on 125, the certificate on
  119, always E on 113, the relative margin on 107.
- (C) The plane (17² kernels of width 0.15 at spacing 0.125; 120 targets): every path exact with three kernels; both skeletons right on
  120, the atoms subtracted on 117, the relative margin on 102, the first entrant on 98, the certificate on 89, always E on 83.
- (D) Designed families: (a) lattice pairs 0.3, 0.4, 0.5 apart at the ratios 10⁻³ … 10³ — at 0.3 the outcome changes with the ratio
  while the certificate's margins are −0.34 and +0.22 —, the fitted skeleton right on 54 of 54, the linear skeleton on 52; (b) 240 pairs at
  1.5σ … 4σ over a grid cell: exact with E on 84, both skeletons and the unpolished flow right on 240, the certificate on 151; (c) 30 of
  A's targets in all six orders of the amplitudes 1, 10^{−1/2}, 10^{−1}: no outcome changes with the order. Outside |m| ≤ 0.8q the
  certificate is wrong on 25 of 333 designed targets (|m|/q up to 1.62). Over all 1154 targets: the fitted skeleton 1154, the linear
  skeleton 1149, the atoms subtracted 1100, the relative margin 1014, the first entrant 977, the certificate 903, always E 763.
- (E) The two-feature map (4802 points): the fitted skeleton's threshold D*(ℓ) at most 4.67σ, without the second feature's own near
  ties (m₀ ≤ 0.05) at most 3.33σ (same signs) and 2.83σ (opposite); the outcome constant at both ends of the ratio on 48–49 of 49
  separations; the linear skeleton's sign the fitted skeleton's at 4719 points; |m_sk − m_fit| against the amplitude ratio of slope
  0.98–1.00 (first order in the tail); the path at 240 sampled points as the fitted skeleton at all 240.
- (F) The range 10⁶ (amplitudes 1, 10⁻³, 10⁻⁶ in all six orders at 5σ and 3σ, exits only below 10⁻¹² of the largest coefficient): at 5σ
  every path exact with E, at 3σ exact with E on 2 of 12 (one an exact fit with four kernels); both skeletons right on all 24.

Tests (`tests/test_transit_certificate.py`, 19): the certificate's and the skeletons' margins, the fitted skeleton on negative features
and on stored targets that break its mutants, the atoms baseline, the entrants near the window's edge, the ascent's end, the outcome's
counts and absolute tolerance, the bounded polish, the polish safeguard on stored paths in both implementations.

## Remark D‴.16 — dither and the mass constraint (`exp_dither_ot.py` → `chaos/dither_ot.md`)

- Dithered averaged slow loop (Langevin in slow time, Gibbs law exp(−J_m/τ₀)), m = 15 from the counts (6, 5, 4), matching prior,
  σ₀ = 0.15σ_crit; τ₀ = 0.005, 0.01, 0.02 over 50 000, 200 000, 50 000 slow units; two seeds.
- Nearest-atom leak: record-level label changes 1.0·10⁻⁵, 5.1·10⁻⁴, 5.6·10⁻³ per slow unit (about half undone at the next change;
  Arrhenius barrier 0.048, 0.043–0.052 at 95%, among the ridge minima 0.040–0.074); time at the balanced counts 0.00, 0.54, 0.25 against the
  energy–multinomial law 0.86, 0.58, 0.33 (standard errors ≈ 0.1 at τ₀ = 0.01; the law under-weights the saddle count vectors 2.7–3.2 times).
- OT leak: balanced counts in 100%, 99.4%, 88% of the records after the initial transit (99.6%, 99.2%, 88% of the whole runs); the excursions are returns
  of the crossing particle (121 of 121 at τ₀ = 0.01, 298 of 336 at 0.02, the rest exchanges); off-balance fraction 6.3·10⁻³, 0.117, effective
  barrier 0.058: the mass constraint makes a crossing reversible without raising its barrier.

## Theorem C′ — the barycentric OT leak (`exp_barycentric_leak.py` → `barycentric_leak/`)

| experiment | result |
|---|---|
| pure leak, m = 6 toward a 48-atom prior, 3 starts | W² non-increasing at every step; end point centroidal (V_bar² ≤ 2.2e-25, W² = Var_γ = 1.164, 1.213, 1.203); a Lloyd fixed point; best of 20 Lloyd runs 1.171 (between the end points: a local minimiser, as for Lloyd) |
| m₀ = m | exact decay of Theorem C(i) to 5e-16 |
| closed loop, representable plant, m = 4, k_e = 4, σ_m = 0.5: (A) OT leak to the ideal / (B) barycentric leak to a 32-atom cloud / (C) labelled leak to a random 4-subset | W(μ_T, μ*) = 0 / 0.167, 0.097, 0.157 (= b_bar = 0.168, 0.097, 0.159, against b = W(μ*, cloud) = 0.45, 0.33, 0.40) / 0.61, 1.02, 0.35 (≈ W(μ*, subset)); RMS‖e‖ on [20, 40]: 0 / 0.002–0.007 / 0.03–0.04 |

Tests (`tests/test_barycentric.py`): the plan reduces to the assignment for m₀ = m and W² = V_bar² + Var_γ; the pure leak
dissipates at every step with D⁺Z ≤ −σ_mV_bar² and ends at a Lloyd fixed point; √C_bar ≤ b + 2η on 90 configurations with
C_bar = b_bar² at small η; the Lyapunov inequality (iii) along a closed loop with a 24-atom prior.

## Theorem J′ — the stationary law of the dithered single-kernel loop (`exp_dither_stationary.py` → `dither/`)

Single-kernel instance (n = 1, m = 1, ϑ = (α, γ), frozen width, G = diag(2, 2), σ₀ = 0.4, τ₀ = 0.1), reference traversed both
ways; 2·10⁴ independent replicas per (k_e, orientation) started from the first-order law, k_e = 10, 20, 40, 80; the Gibbs
measure, the first-order correction g, the current j₁ and the circulation coefficient on a 201² grid. Full tables in
`results/dither/table.md` (raw in `raw.npz`).

| | k_e = 10 | 20 | 40 | 80 |
|---|---|---|---|---|
| (m_ε(α) − π(α))·k_e, prediction −0.255 | −0.28 | −0.36 | −0.36 | −0.37 |
| (m_ε(γ) − π(γ))·k_e, prediction −0.426 | −0.44 | −0.54 | −0.51 | −0.67 |
| residual m_ε(α) − first order (s.e. 0.0013) | −0.0021 | −0.0051 | −0.0025 | −0.0015 |
| residual m_ε(γ) − first order (s.e. 0.0018) | −0.0015 | −0.0056 | −0.0020 | −0.0030 |
| orientation-odd shift of γ, measured/predicted | 0.93 ± 0.04 | 1.10 ± 0.07 | 1.15 ± 0.15 | 1.38 ± 0.29 |
| k_e²E e²/(2π(L̄)), orientation +1 | 1.054 | 1.033 | 1.017 | 1.009 |
| circulation Λ_ε (+1), measured ± 0.0005 / predicted | +0.0037 / −0.0005 | +0.0006 / −0.0003 | −0.0004 / −0.0001 | +0.0003 / −0.0001 |

- The Gibbs discrepancy decays like k_e^{−0.87} (α), k_e^{−0.83} (γ), k_e^{−1.09} (αγ), and the first-order prediction removes
  93–97% of it at k_e = 10 and 70–85% at k_e = 20–40; the residual is at the Monte Carlo floor at k_e = 80 and its second-order
  decay is not resolved.
- The orientation-odd shift (the dither version of Remark H.7a) is predicted to 7–8% at k_e = 10 and within the 10–30%
  resolution at larger k_e.
- The tracking error is the Gibbs-averaged occupation loss over k_e²/2 with an O(1/k_e) correction of 0.5–0.7/k_e.
- The first-order circulation is below the Monte Carlo resolution for k_e ≥ 20 and dominated by higher orders at k_e = 10; the
  current machinery of Proposition J′.2 is checked on the first-order averaged SDE on the grid: ∫Θ∧(j_ε − j₀)/(ε∫Θ∧j₁) = 1.47,
  1.21, 1.09, 1.03 at ε = 0.2, 0.1, 0.05, 0.025 and ‖m_ε − m₀(1 + εg)‖₁/ε² = 0.31 constant.
- Referee round 6 (`review/REVIEW_round6_Jprime.md`): the first run's batch standard errors were biased low by 1.4–1.9 and its
  burn-in too short (JP-4–JP-6); the run reported here is the corrected one.
- Corollary J′.3 (`exp_gibbs_m.py` → `dither/gibbs_m.md`): μ*_τ on the grid (∫α dμ* = 0.4963, ∫γ dμ* = 0.1873, ∫α² dμ* = 0.9353,
  L̄(μ*) = 0.1394); the Gibbs measures π_m for m = 1–32 by MALA (4000 replicas): E[L̄(μ_Θ)] = 0.1960, 0.1681, 0.1540, 0.1470, 0.1432,
  0.1413 with (E[L̄] − L̄(μ*))·m = 0.057–0.061 (slope −0.99), E[γ₁] slope −1.01, E[α₁²] from 0.536 to 0.922 with rescaled deviation
  −0.41 ± 0.01; the m = 1 row matches the quadrature of π₁ within one standard error.

## Theorem I — reciprocity at order k_e⁻² (tests)

No new experiment: every statement is checked in `tests/test_nonreciprocity.py` against the generic second-variation
code of `src/core/averaging.py`.

| check | result |
|---|---|
| A(ϑ,ϑ') = aᵀB(ξ,ξ')a' against the antisymmetric part of S_even + K, n = 1 and 2 | agreement to 1e-10 |
| A(ϑ, ϑ_z) = 0 for the zero feature; cyclic sum over (ϑ₀, ϑ₁, ϑ_z) equals A(ϑ₀, ϑ₁) | 1e-14 |
| n = 1, f := f̂_μ (perfect fit on the orbit) | max |A| < 1e-13 over 20 random pairs |
| n = 1, non-zero residual | |A| > 1e-3 found; equals the closed form −aa'((c−c')/s²)∫rφ_cφ_c'dρ to 1e-12 |
| κ(ξ,ξ') against −(1/(2s²))e^{−|c−c'|²/(4s²)}(c−c')·F_s(m), n = 2 | 1e-10 |
| F₂((0.5, 0)) > 0 on the ellipse; κ ≡ 0 on a segment traversed back and forth | F₂ > 1e-3; |κ| < 1e-12 |
| n = 2 with f ≡ f̂_μ on ℝ² (residual and its Jacobian vanish on the orbit): B = 2κI | 1e-13; A ≠ 0 on the ellipse, A ≡ 0 (< 1e-13) on the segment |
| n = 2 with f = f̂_μ + he₁, h zero on the ellipse (residual zero on the orbit, Jacobian not symmetric): B − 2κI = −∫φφ'∂₂h dρ·(e₁e₂ᵀ − e₂e₁ᵀ) | 1e-12, and not small (the round-5 correction of Theorem I(iv)) |
| positivity witness f ≡ 0, μ = δ_ϑ̄ | integral > 0; generic code equals the closed form to 1e-12 |
