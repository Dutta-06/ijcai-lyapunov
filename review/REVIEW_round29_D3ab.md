# Referee report D3ab: Propositions D‴.41 and D‴.42

Object: Propositions D‴.41 (measurement noise on the gauge) and D‴.42 (several particles per class) in their drafts (scratchpad
`label/D41_draft.md`, `multi/D42_draft.md` and their pieces), `src/core/label_noise.py`, `src/core/mass_transfer.py` (`run_dithered`'s
`label_*` options, `fixman_classes`), `experiments/exp_gauge_label.py`, `experiments/exp_gauge_classes.py`, `tests/test_gauge_label.py`,
`tests/test_gauge_classes.py`, at 74e955e. The report as delivered (the referee's scratch paths are under the session's scratchpad).

Both propositions have proofs that hold only at their leading order and are then stated more broadly than they were derived. Five issues are major:

- **D‴.41(iv)** is false beyond short noise correlation times. The proof drops a cross term.
- **D‴.42(ii)/(iii)**: the frozen split holds only for the isotropic dither. The claim "D‴.40's statements hold per class likewise" fails for D‴.40's own masked noise.
- **D‴.42(iv)**: the cohesion's residual is misreported, and the pieces overstate the measurement-noise cancellation.
- **Tests**: 11 of 21 mutants survive, including a Heun scheme that silently becomes Itô-like.
- **`limit_drift`**: it silently returns garbage at degenerate eigenvalues.

## What I checked and how

- **Setup.** I worked only in a scratch copy, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review29/repo`.
  - I pulled the coordinator's commit 74e955e into it (the `wong_zakai(weights=…)` argument and Part A's heteroscedastic row) and re-read the updated D41_draft (ii).
  - The repository is untouched: `git status` is clean and HEAD is 74e955e.
  - All runs used OMP_NUM_THREADS=1, nice 10, at most 2 processes, no EXP_CACHE.
  - Scripts are in `…/scratchpad/review29/work/` and outputs in `…/scratchpad/review29/out/`.
- **Proofs.** I re-derived, line by line:
  - (i) and the antisymmetry of B.
  - The symmetry of K(−M∇²J)M for coloured noise.
  - The slow-limit projection ∂Φ and the components of ⟨(U·∇)U⟩; a′² − aa″ = 1 − q², ∂_β(βa′)² = 2β(1 − q²), h″ = 2 + 6c²/W⁴.
  - (iv)'s Euler expansion and (vi)'s identities.
  - D‴.42's Jacobian/block argument and the g-orthogonality of the fractions.
  - I checked numerically that g·A is symmetric to 10⁻¹¹ (m = 3, 6, 9). So ∂Φ is indeed the g-orthogonal projection.
- **Reproductions.**
  - D‴.41 Parts A, B, D, E, F (`GL_WORKERS=1 … --only A,B,D,E,F`, 68 s, `out/gl/gauge_label.md`).
  - D‴.42 Parts A–D (`out/gc/`): byte-identical to `results/chaos/gauge_classes.md`.
  - Reduced Part C: 2 seeds × 4·10⁴ slow units for fresh noise and AR(1) ρ = 1/3, 0.6, 0.8 (`c_reduced.py`).
- **New computations.**
  - (iv)'s dropped cross term: the scalar model's exact per-step identity up to ρ = 0.999 (`ar1_scalar.py`), and a linear-response evaluation of the Euler drift including h²∂²Φ[F, ⟨Uξ⟩] for the replicator (`ar1_replicator.py`, `ar1_tauc.py`).
  - The Stratonovich residual for heteroscedastic profiles, for orbit-correlated real-time noise (one-sided kernel, giving a symmetric part and a Lévy-area part) and for an offset noise, at γ = 10⁻¹ to 10⁻⁵ (`corr_noise.py`, `corr_gamma.py`).
  - The Itô selection with the Euler step in a instead of α, and the γ = 10⁻² zero of the Itô drift (`coords_target.py`).
  - For D‴.42 (`d42_*.py`):
    - an eigh-based (g-symmetrised) `limit_drift` at the degenerate points, including n = 3 equal splits;
    - D‴.40-type masked dither at γ = 10⁻³ to 10⁻⁶;
    - Part D per component;
    - the fractions' Itô drift against candidate potentials.
- **Tests.** The 10 new tests pass, as does the gauge set (38 tests). I ran 21 mutants (`mutate.py`, results in `out/mutants.txt`).
- **Literature.** Web search only: fetches of arxiv, neurips and pmlr fail on DNS. Contents come from search snippets or memory, as flagged.

## Verdicts

### D‴.41

**(i) No diffusion along M.** Proved and correct. Diffusion of 2·10⁻²⁰ reproduced.

**(ii) The Stratonovich cancellation.** The algebra is correct for homoscedastic noise independent over the nodes, and Part A reproduces it (4.6·10⁻⁹ to 1.1·10⁻⁸). The coloured, continuous-time case is correct by linear response but is tested nowhere. Problems:

- **The new heteroscedastic sentence is framed wrongly.** On this gauge, every non-white covariance I tried leaves a residual that is O(γ) relative to its own Itô drift.
  - At W₀, the ratio is about 2.4γ for an offset noise and 0.57γ for 1 + 0.8 cos 2t. For the coordinator's 1 + 0.8 sin t it is 1.5·10⁻⁵ at γ = 10⁻³ and 6·10⁻⁸ at γ = 10⁻⁵.
  - The ratio is 5–12% at γ = 0.1 and 2·10⁻⁵ or less at γ = 10⁻⁵.
  - So in the slow-replicator limit the cancellation does not "need" homoscedasticity. A second mechanism works there: Φ̃_k ∝ ∂_{α_k}Φ̃_k, so the masses' noise is the amplitudes' carried by the gauge, and ψ = O(γ²) near M.
  - The sin t profile is a weak test. 1 + 0.8 cos t gives exactly zero, because x = 2 sin t visits every x at t and π − t.
  - The quoted "0.13–0.23%" does not match the code. It gives 0.092–0.133% at γ = 10⁻² and 6–7·10⁻⁸ at γ = 10⁻⁵.
- **"Physical noise of any spectrum" (BULLET, ROW0, README, AUDIT) is overstated.**
  - Real-time measurement noise maps to node-independent, occupation-homoscedastic noise only when it is white on the orbit's resolution (τ_n ≪ T/N). Only then is it Stratonovich with a cancellation; this is the answer to the brief's question.
  - For τ_n ≳ T/N the averaged noise is correlated across nodes, C_ij = Σ_k R(t_i − t_j + kT). The residual is then 0.08%, 1.1% and 2.1% at τ_n/T = 0.016, 0.16 and ≥ 1.6 (γ = 10⁻²), again O(γ).
  - The Lévy-area part vanishes here only because x(π − t) = x(t).
  - If the feedback uses the measured error, the adaptation sees the noise high-passed by the loop. That noise is correlated over the orbit, so the hypothesis fails.

**(iii) The Itô selection and the closed form.** The derivation is correct, and Part B reproduces it (W_ℓ and the eigenvalues −0.214, −0.128). Caveats:

- **(a) The selection depends on the coordinates of the update.**
  - An Euler step in the physical amplitude a instead of α selects with (1 − q)(1 + 3q), not 1 − q². Katzenberger's drift confirms this to 2·10⁻⁴.
  - At W₀ (γ = 10⁻⁵) the drift is (−1.64, 2.30, −0.66)·10⁻⁷ instead of (−1.00, 1.55, −0.55)·10⁻⁷. At W = 1/3 it is 7.7× larger on W₁.
  - The slow-limit target becomes (0.297, 0.451, 0.252) instead of (0.316, 0.412, 0.272).
- **(b) W_ℓ is the slow-limit target.** At Part C's γ = 10⁻² the Itô drift vanishes at (0.3221, 0.3999, 0.2780).
- **(c) Near the saturation the label-noise fitness vanishes** like 1 − q², while the dither's diverges. Label noise cannot keep the masses off the saturation against the dither.
- **(d) "SGD" should read label-noise SGD.** Minibatch noise vanishes on M.

**(iv) The per-step-variance law.** False beyond λ_max τ_c ≪ 1.

- **The proof's error term is wrong.** It writes Φ(z_{n+1}) − Φ(z_n) = h⟨ψξ⟩ + (h²/2)∂²Φ[⟨Uξ⟩, ⟨Uξ⟩] + O(h²|F|²). It drops h²∂²Φ[F, ⟨Uξ⟩], which is O(h²|F||ξ|). The mean of that term vanishes for independent draws but not for correlated ones. Its relative size per normal mode is about hλρ/(1 − ρ + hλ), i.e. λτ_c/(1 + λτ_c).
- **The draft's own scalar model refutes it.** Its identity is exact: the decay rate is h·E[(ξ − r)²], not h·σ².
  - At ρ = 0.9, h = 0.01 (λτ_c = 0.38) the selection is 0.750 of (1 − ρ)/(1 + ρ).
  - At ρ = 0.99, h = 10⁻³ it is 0.704.
  - At τ_c ≈ 1 (λτ_c = 4.2, ρ = 0.99 at h = 0.01 or ρ = 0.999 at h = 10⁻³) it is 0.197 and 0.192.
  - The last case has hλ = 0.004, so the stated O(hλ_max) caveat is met.
  - The linear-response formula 1 − 2hλρ/(1 − aρ) + h²λ²(1 + aρ)/((1 − a²)(1 − aρ)), with a = 1 − hλ, matches each case to 0.3%.
- **The replicator, at Part C's settings, confirms it.**

  | noise | (iv)'s law | measured (2 seeds, on W₂) | my linear-response prediction |
  |---|---|---|---|
  | fresh | 1 | 1.011 | 1.0088 |
  | AR(1), ρ = 1/3 | 0.5 | 0.5007 | 0.5000 |
  | AR(1), ρ = 0.6 | 0.25 | 0.2461 | 0.2458 |
  | AR(1), ρ = 0.8 | 0.1111 | 0.1051 (−5.4%, about 5.7 s.e.) | 0.1048 |

  The queued full Part C will therefore contradict (iv) at ρ = 0.8.
- **The sampled-OU claim tanh(h/2τ_c) overstates the selection.** At Part C's h, my prediction gives 0.95, 0.75, 0.59 and 0.39 of tanh(h/2τ_c) at τ_c = 1, 5, 10 and 20.
- **The status "up to O(hλ_max)" is wrong.** The neglected term is O(λ_max · max(h, τ_c)).

**(v) Additivity with dither.** Correct at leading order, subject to (iii)(c).

**(vi) The examples.**
- The identities are exact.
- The decay-rate statements inherit (iv)'s restriction, including the "1/k for redraw every k steps".
- The scalar example is Ziyin–Li–Ueda's law of balance. Part E's model is Pillaud-Vivien–Reygner–Flammarion's. Neither is cited.
- The generic tanh network's ratio reproduces as 0.995.

**Numbers against the code.**
- A, B, D and E reproduce to the digits quoted.
- F's "0.996" is 0.995 in the code's output.
- E's "each within 1% over the 40 coordinates" fails for ρ = 0.6 (minimum 0.2471, −1.2%) and for k = 4 (0.2468, −1.3%).
- In E, ‖w‖₁ rises (1.977 → 1.998 for fresh noise, → 2.057 for ρ = −1/3). The selection descends the weighted ℓ¹ norm Σc_i|w_i|; it does not move toward the ℓ¹-minimal interpolant.
- F's AR(1) and block ratios (0.491, 0.237, 0.258) scatter by up to 5%.
- All of D, E and F sit at λτ_c ≤ 0.03, so they cannot test the "any correlation" claim.

### D‴.42

**(i) The determinant.** Proved.
- The constant also contains pK·log(m/3): the fit block is (m/3) times the one-per-class Hessian at W, not "at the masses B_k".
- c₉ = 2c₆ is an identity, not a coincidence.

**(ii) The 3/m law and the frozen fractions.**
- Proved, and reproduced, for the isotropic dither 2τG (the Gibbs case). The fractions' drift is ∝ γ from 10⁻³ to 10⁻⁶.
- The diffusion ≡ 0 holds for any position noise; the g-orthogonality argument is correct.
- The O(γτ) freezing fails for D‴.40's masked dither. At the split (0.2, 0.7, 0.5):
  - Mask (1,0,0): the fractions' drift per γτ tends to (−2.80, −5.35)·10⁻⁴, unchanged from γ = 10⁻⁴ to 10⁻⁶. That is O(γτ). The class drift misses the per-class law by +2.2% and +1.4%.
  - Mask (0,1,1): the class masses get (0.0041, −0.0017)·γτ. One particle per class gives only about 7·10⁻⁹·γτ (O(γ²τ)).
- The statement of (ii) has a sign slip: the fitness should be −(τ/2m)∂_{W_k}ℓ₁.

**(iii) The correction per class.** Correct (Part C reproduces). The sentence "D‴.40's statements hold per class likewise" is false; see (ii).

**(iv) Measurement noise with several particles per class.**
- The breaking at unequal splits is real, but it is misreported. "2.2–2.6%" is normalised by the fractions' Itô drift.
  - On the class masses the residual is 9.6% and 1.4% of their own Itô drift at γ = 10⁻⁵, split (0.2, 0.7, 0.5); 2.4% and 3.1% at split (0.3, 0.6, 0.8).
  - It is 8.9% at γ = 10⁻², so it does not vanish as γ → 0.
- The exchange-symmetry argument covers only the fractions' component. The cancellation on the class masses at the equal split follows instead from an invariant subspace: co-located, equal-mass particles get identical drift and identical noise, so the dynamics there is the m = 3 loop exactly. Consequently the equal-split Itô class drift is m-independent: (−6.421·10⁻⁹, 8.613·10⁻⁸), equal to m = 3's, not 3/m of it.
- The Itô potential "amplitudes' part ∝ Σβ_i²" is 4–13× too small.
  - The measured fractions' drift is (5.39, −3.81)·10⁻⁷. The full trace (v/4)⟨‖∇Φ̃_k‖²_G⟩Σβ_i² predicts (5.51, −3.92)·10⁻⁷; the amplitude part predicts (0.43, −0.50)·10⁻⁷.
  - The centres' and widths' parts are not constant along the fractions.

**Is `limit_drift` valid at m = 6 and 9?**
- At m = 6 yes, including the 12 degenerate pairs at the equal split: cond(V) = 182, and it agrees with the eigh-based basis to 1.4·10⁻⁸.
- At m = 9, at the splits the experiment uses, yes (4·10⁻⁷).
- At n = 3 equal splits it fails:
  - silently at γ = 10⁻², returning a class drift of about 3·10³⁰ with cond(V) = 3.7·10¹⁸ and no exception;
  - or with LinAlgError at W = (0.30, 0.35).
- The eigh-based basis gives the claimed results there: ratio 0.3333, fractions' drift ≤ 3·10⁻⁸ per γτ.

### Tests and mutants
- The 10 new tests pass.
- Killed: L1 to L4, L7, L8, R8, F3, F4.
- Survived (11), listed as mutant id: what the mutation does.
  - **L5:** `wong_zakai`'s weights ignored.
  - **L6:** `wong_zakai`'s weights squared.
  - **R1:** per-step variance without N.
  - **R2:** AR(1) factor inverted.
  - **R3:** AR(1) innovations non-stationary.
  - **R4:** Heun corrector with a fresh draw.
  - **R5:** Heun corrector without noise.
  - **R6:** Heun's final log-mass step reduced to Euler.
  - **R7:** Heun keeps `label_rho`.
  - **F1:** the classes branch ignores the masses form.
  - **F2:** the floor is not scaled by W_k.
  - **F5:** per-class τ_k averaged.
- The tests check only the degenerate reductions (ρ = 0, and Heun at v = 10⁻³⁰).
- Heun does converge to the Stratonovich solution in the flat coordinates. Its predictor evaluates the drift at un-normalised β, an O(γh) bias of about 0.25% here.

### Literature

Present in [S69], with details as given: Wong–Zakai 1965; van Kampen 1981; Horsthemke–Lefever 1984; Blanc et al., COLT 2020; Damian–Ma–Lee, NeurIPS 2021; HaoChen et al., COLT 2021; Li–Tai–E, JMLR 2019; Li–Malladi–Arora, NeurIPS 2021; Orvieto et al., ICML 2022. Shalova–Schlichting–Peletier needs pages: JMLR 27(146):1–56, 2026 (search-verified).

Missing:
- Ziyin–Li–Ueda, "Law of balance…", arXiv:2308.06671 (search-verified). Under SGD, d/dt(‖u‖² − ‖w‖²) is negative; this is (vi)'s scalar law.
- Ziyin–Wang–Li–Wu, "Parameter symmetry and noise equilibrium of SGD", NeurIPS 2024, arXiv:2402.07193 (search-verified). Its "Noether flow" covers the diagonal network's u_iv_i decay.
- Pillaud-Vivien–Reygner–Flammarion, COLT 2022, PMLR 178:2127–2159, arXiv:2206.09841 (search-verified). A continuous-time (Itô) label-noise flow on the quadratic parametrisation solves the Lasso. This is Part E's model, and it contradicts the README's "no selection in continuous time" unless the Itô/Stratonovich reading is made explicit.
- Kupferman–Pavliotis–Stuart, Phys. Rev. E 70:036120 (2004), and Pavliotis–Stuart (2005, arXiv math-ph/0504054) (search-verified). The limit is Itô, Stratonovich or neither depending on the ratio of the correlation time to the relaxation time. This is the crossover (iv) misses.
- Kühn–Rosenow, arXiv:2306.05300 (snippet). Correlated epoch-based SGD noise, with a curvature crossover.
- Woodworth et al., COLT 2020 (from memory). The effective initialisation scale in diagonal networks.

Novelty: that the regulariser of label-noise SGD scales with step × variance is standard (Blanc, Damian, Li–Wang–Arora). Stratonovich conservation of the examples' charges follows at once from Ziyin et al.'s Itô derivation. I did not find the general antisymmetry cancellation stated anywhere.

## Table

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ab-1 | major | (iv) is false beyond λ_max τ_c ≪ 1. The proof drops h²∂²Φ[F, ⟨Uξ⟩]; the status "O(hλ_max)" is wrong. The scalar model gives 0.75, 0.70, 0.20, 0.19 of the law; the replicator gives 0.1051 against 0.1111 at ρ = 0.8; tanh(h/2τ_c) is off by up to 2.5×. | D41_draft (iv) statement, proof and status: keep the cross term; state the condition hλ_max ≪ 1 − ρ and the linear-response factor. Add ρ = 0.9–0.99 and OU with τ_c ≳ 1/λ to Parts C and D. Fix the `ar1_factor` docstring and the ROW0, BULLET, README, AUDIT pieces. |
| D3ab-2 | major | D‴.42's frozen fractions and "D‴.40 per class" fail for non-isotropic position noise. Masks (1,0,0) and (0,1,1) give an O(γτ) fractions' drift and an O(γτ) class selection that is absent at m = 3. The proof relies on the Gibbs free energy. | D42_draft (ii), (iii) and the title: restrict to Σ = 2τG and state the counterexample. Add a mask part to `exp_gauge_classes.py`. |
| D3ab-3 | major | Part D's residual is normalised by the fractions' Itô drift. On the class masses it is 0.7–9.6% and does not vanish as γ → 0. BULLET, README and AUDIT ("the cancellation acts on the class masses as with one particle per class"; "continuous-time adaptation does not select") are false at unequal splits. | Report Part D per component (`exp_gauge_classes.report`, numerics). Correct pieces_D42, the D41 AUDIT and the README. |
| D3ab-4 | major | 11 of 21 mutants survive. All label-noise paths of `run_dithered` are untested beyond the reductions, including a Heun corrector with a fresh draw; also `wong_zakai`'s weights and `fixman_classes`' masses form and floor. | tests/test_gauge_label.py: per-step variance and AR(1) moments with a stub rng; one Heun step against a hand-coded step at non-zero noise; the weighted Wong–Zakai term against the formula. tests/test_gauge_classes.py: the masses form off M and the floor near saturation. |
| D3ab-5 | major | `gauge_drift.limit_drift` returns garbage silently at degenerate eigenvalues (n = 3 equal split, γ = 10⁻²: 3·10³⁰) or raises LinAlgError. | Use the eigh of g^{1/2}Ag^{−1/2}, or guard on cond(V). Add n = 3 equal-split rows to Part B, plus a test. |
| D3ab-6 | minor | The heteroscedastic caveat: the residual is O(γ) and the condition is not needed in the slow limit; sin t is the weakest profile; the quoted 0.13–0.23% does not match the code's 0.09–0.13%. | D41_draft (ii): state the γ-scaling and the second mechanism. Part A: add cos 2t and offset profiles and a γ-scan. |
| D3ab-7 | minor | Orbit-correlated (non-white) physical noise and the feedback path are not covered; "any spectrum" is overstated (2.1% at γ = 10⁻², about 7% at γ = 0.1); the Lévy-area term vanishes only by the orbit's symmetry. | D41_draft (ii): put "white on the orbit's scale" in the hypothesis. Correct the BULLET, ROW0 and README pieces. |
| D3ab-8 | minor | The Itô selection depends on the update's coordinates: an Euler step in a gives (1 − q)(1 + 3q) and target (0.297, 0.451, 0.252). | D41 (iii) and (iv): say "an Euler step in z = (Θ, log β)" and give the a-coordinate variant. |
| D3ab-9 | minor | D‴.42's Itô potential uses the amplitudes' part only. The full trace is right (to within 2–3%); the amplitude part is 4–13× too small. | D42_draft (iv) proof. |
| D3ab-10 | minor | D‴.42 (iv)'s symmetry argument covers the fractions only. The equal-split Itô class drift is m-independent, not "with m the number of particles". | Use the invariant-subspace argument in D42 (iv). Correct the README and ROW0 pieces. |
| D3ab-11 | minor | [S69] misses Ziyin et al. (two papers), Pillaud-Vivien–Reygner–Flammarion, Kupferman–Pavliotis–Stuart and Pavliotis–Stuart, Kühn–Rosenow, Woodworth et al. The per-update-variance novelty is overstated. | LIT, AUDIT and the D41 status: add these and narrow the novelty claim to the general antisymmetry identity. |
| D3ab-12 | minor | "The target of the per-update selection" W_ℓ is the slow-limit target; at γ = 10⁻² it is (0.322, 0.400, 0.278). | numerics_D41 (B) and the RESULTS piece. |
| D3ab-13 | minor | The numerics sit only where the claims hold trivially (λτ_c ≤ 0.03); F scatters by up to 5%; the coloured continuous-time claim is untested. | `exp_gauge_label.py` D, E, F: extend ρ; report F's deviations; add a Heun run with OU noise. |
| D3ab-14 | minor | Near the saturation the label-noise fitness → 0 while the dither's diverges; "against the dither's drive" needs this qualification. | D41 (iii) and (v). |
| D3ab-15 | cosmetic | F's 0.996 should be 0.995; E's "within 1%" is false; E's "ℓ¹-minimal" should be weighted ℓ¹ (‖w‖₁ rises). | numerics_D41; (vi). |
| D3ab-16 | cosmetic | D‴.42: the sign of the fitness in (ii); the constant c includes pK·log(m/3); "Hessian at the masses B_k" is imprecise; Ŵ_k is undefined. | D42_draft. |
| D3ab-17 | cosmetic | D‴.41: F is used for both the flow and the target; Φ_k is undefined; "equal h′_k" should be "equal P_k h′_k"; "SGD" should be "label-noise SGD"; [S69] lacks Shalova et al.'s pages; the pieces lack the homoscedasticity condition. | D41_draft; pieces_D41. |
| D3ab-18 | cosmetic | Heun's predictor evaluates the drift at un-normalised β (an O(γh) bias, about 0.25%). | `mass_transfer.run_dithered`: renormalise lp before the corrector. |

## Lesson

Each of the major failures comes from promoting a leading-order identity to a universal quantifier ("any stationary sequence", "likewise per class", "needs homoscedasticity") without testing it where its derivation's assumptions are weakest. In (iv) a dropped term was written as O(h²|F|²) when it is O(h²|F||ξ|), and its mean is zero only for independent draws. D‴.42's freezing borrowed the Gibbs free energy, which holds only for the isotropic dither. The heteroscedastic caveat rests on one profile at one γ. Each time, the numerics sat where the claim is trivially true: λτ_c ≤ 0.03, isotropic noise, sin t at γ = 10⁻². The tests pinned only the degenerate reductions of the very code paths that carry the novelty: Heun's shared draw, the AR(1) normalisation, the weights. The remedy is to probe every claim at the edge of its hypotheses and to write tests that would fail if the novel mechanism were silently replaced by the old one.
