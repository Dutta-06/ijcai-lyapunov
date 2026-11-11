# Referee report D3y: Proposition D‴.36

## Summary of what was checked and how

Object: Proposition D‴.36 (PARADIGM.md l. 4463 ff.), its rows in §0 (l. 146) and §7 (l. 6972), the §8 bullet (l. 7236), RESULTS_PARADIGM.md
(l. 932), [S64] and the audit item in LITERATURE_TRANSPORT.md (l. 84, l. 539), PAPERS.md paper 13, README l. 77–78; code
src/core/mass_transfer.py (`eps`, `run_dithered`), experiments/exp_gauge_dither.py, tests/test_gauge_dither.py, results/chaos/gauge_dither.md / _raw.npz.
Nothing in the repository was modified. All scripts are in
`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref26/`.

Independent work (my own code, not the repository's):
- `mymodel.py`: my own version of the instance (quadrature X = 2 sin(2πk/256), the three bumps, a = a_max·tanh(α/a_max) as in `net.physical`,
  c = c_max·tanh, s = s_min + (s_max − s_min)·sigmoid), J = m·½⟨r²⟩ + (ε/γ)Σβ log β, analytic gradient (checked against finite differences of J
  to 6·10⁻¹²), drift dϑ = −G∇_ϑJ, dβ = −γβ(∂_βJ − λ). It matches `ReplicatorCohesionLoop.drift` to 10⁻¹⁷ (`check1.py`).
- `check2.py`: the Jacobian of the drift on {Σδβ = 0}: its spectrum, kernel and spectral projector, compared with T(TᵀgT)⁻¹Tᵀg, and D.
- `check3.py`: the nonlinear limit map Φ by LSODA (rtol 10⁻¹², 1.5·10⁵ slow units) from perturbed points; its first derivative, and the
  Katzenberger drift τΣ_iG_i∂²_iΦ from second differences (h = 0.02 and 0.01, which agree to four digits).
- `check4.py`: the full 11-dimensional linear Ornstein–Uhlenbeck covariance with the leak (Lyapunov equation), against the reduced κ of (iii).
- `check5.py`: the minimiser with the leak where w is not representable (A_1 = 1.5), for several (γ, ε).
- `check6.py`, `check8*.py`: D and the normal spectrum near the saturation. `check7.py`: the predicted stationary mean under the leak.
  `check9.py`: two particles per collapsed class. `check10.py`: D averaged along the simulated paths.
- `mysim.py`: my own batched Euler–Maruyama, the same scheme (dt = 0.25), with more seeds than the repository: 256 seeds of B at τ = 10⁻⁴
  and 256 at τ = 2.5·10⁻⁵, and 64 seeds of D at each of ε = 10⁻⁴ and 10⁻³ (2·10⁵ slow units each).
- `repro.py`: the repository's `series` for B seeds 100, 137 and 163 and D seed 2 (ε = 10⁻³). All four match the raw records **bit for bit**
  (max |diff| = 0). `--from-raw` (with GD_OUT pointing at a copy) regenerates gauge_dither.md **byte-identically**.
- `mutate.py`, `mutate2.py`: mutants run against a `git archive` copy of the repository.
- Every number in the PARADIGM text, the §0 and §7 rows, the §8 bullet and RESULTS was checked against the md. All of them match, up to rounding.

## Verdicts

**(i) Correct for one particle per class; false as stated for collapsed classes.** At the instance M = {W_ka(α_k) = A_k, c and s fixed} has
dimension K − 1 = 2. f̂ is constant on it (residual 2·10⁻¹⁶). J is constant on it: the cohesion is inactive, because the pairwise
‖Δϑ‖²_{G⁻¹} = 39–52 ≫ R² = 1. The tangent δα = −AδW/(W²a′) with a′ = 1 − (a/a_max)² = sech²(α/a_max) is right in the code's coordinates:
‖J_F T‖ = 2·10⁻¹⁰. M is normally attracting: two zero eigenvalues, and nine normal ones in [−0.195, −3.06·10⁻⁴] at γ = 10⁻². The slowest of
these is γ-independent for γ ∈ [10⁻³, 10], and it does not degenerate toward the saturation (−3.0·10⁻⁴ still at W_1 − 0.25 = 10⁻³). With
collapsed classes the zero-mode space has dimension m − 1, not K − 1 (ID-1).

**(ii) Correct.**
- The metric: the code's drift is the gradient flow of J in g = G⁻¹ on ϑ and δβ²/(γβ) on β. With one particle per class δβ = mδW, so
  this is (m/γ)/W_k on the masses, as stated. There is no stray factor of m or 3.
- The projector: g·J_F is symmetric to 6·10⁻¹¹ on Σδβ = 0, and the spectral projector onto the kernel equals the g-orthogonal projector
  to 2·10⁻⁹.
- The limit map: the derivative of the nonlinear limit map from long integrations is (−0.0338, 0.0174, 0.0165) along α_1, the same as the
  projector. Along c and s it is ≈ 10⁻¹⁰.
- Table A is reproduced to all digits for γ = 10⁻³ to 10. The two-class closed form is right.
- The γ limits are right. Near the boundary D_11/δ² → 5.0·10⁻⁵ (quadratic), and in the amplitude coordinate D_α1 → 2τG_α exactly, the free
  dither's rate.
- The drift term, which the paper does not compute, is b_W = (−2.055, +1.580, +0.475)·10⁻⁷ per slow unit. It explains B's means (ID-4).

**(iii) Mostly correct; one wrong order.**
- The leak's drift is the Fisher–Rao gradient of (ε/γ)Σβ log β. The representation at β ≡ 1 is stationary with it. κ = ε(TᵀgT)⁻¹T_W is
  the right reduced stiffness, and κC + Cκᵀ = D is the right Lyapunov equation.
- The reduction agrees with the full 11-dimensional linear covariance to 0.5% at ε = 10⁻⁴ and to 4% at ε = 10⁻³.
- The O(τ/ε) mean shift has the right scaling. My computed shifts are 0.0039 (ε = τ) and 0.0004 (ε = 10τ), and they agree with 64-seed
  simulations.
- "Where w is not representable … a bias of order ε" is wrong: the bias is of order ε/γ (ID-2).
- The linearised spread misses by 14–23% at ε = τ (ID-3).

**Numerics.**
- Reproducible: bit-identical seeds and a byte-identical md. The text agrees with the md.
- Part A is exact. Part B's diffusion is confirmed to about ±10–20%: 256 seeds give variance-slope ratios 0.94 ± 0.08, 1.26 ± 0.12,
  0.99 ± 0.12 at τ = 10⁻⁴ and 1.21 ± 0.11, 0.94 ± 0.10, 1.18 ± 0.11 at τ = 2.5·10⁻⁵. D's variation along the path accounts for only
  −3%, +1%, +1%.
- The weak points: the 4-seed statistics of Part D and their interpretation (ID-3, ID-4, ID-11), the MSD including the drift (ID-5), and the
  thin time-scale separation (ID-6).

## Issues

**ID-1 (major). (i) "(or with its classes collapsed)": the dimension K − 1 and the two-class closed form are wrong for collapsed classes.**
- *Location:* PARADIGM.md l. 4464–4468, and the closed form for T_α at l. 4476.
- *What is wrong:* with n_k > 1 particles at a common position, the within-class weight redistributions (Σ_{i∈k}δβ_i = 0, positions fixed)
  leave f̂ and the cohesion energy unchanged. They are therefore further zero modes of J. The manifold of stationary points through the
  configuration has dimension m − 1, not K − 1. Proposition D‴.34 says so itself ("15 zero modes … 12 redistributions, 3 gauges", l. 144),
  and the §8 bullet lists "the gauge with several particles per class" as *open*.
- *Evidence* (`check9.py`; m = 6, two particles per class, compact cohesion, γ = 10⁻²): **5** zero modes on Σδβ = 0. The class masses' D from
  the spectral projector is (1.42, 0.814, 0.903)·10⁻⁸. The stated closed form, with T_α := Σ_k(A_k/(W_k²a′_k))²/G_α (no n_k), cannot
  reproduce this. It is reproduced exactly with T_α → Σ_k n_k(…)²/G_α and m = 6 in T_W.
- *Fix:* restrict (i)–(iii) to one particle per class. Alternatively, state that M has dimension m − 1, that the redistributions receive no
  noise at first order (they have no position component, so the projection gives them none; they move only through the drift term,
  consistent with D‴.34's Open), and put n_k into T_α.

**ID-2 (major). (iii) "where w is not representable the leak pulls the masses back toward the saturation trap … a bias of order ε": the
order is ε/γ, not ε.**
- *Location:* PARADIGM.md l. 4489–4490.
- *What is wrong:* the leak's energy is (ε/γ)Σβ log β. The replicator's metric scale γ drops out of the stationary set, so the minimiser of
  mL̄ + (ε/γ)Σβ log β depends on ε/γ alone.
- *Evidence* (`check5.py`; target A_1 = 1.5, so W_1 ≥ 0.375 > 1/3 and w is not representable):
  - (γ, ε) = (10⁻², 10⁻⁴), (10⁻¹, 10⁻³) and (1, 10⁻²) give identical minimisers: W = (0.37402, 0.31299, 0.31299), fit error 1.571·10⁻³,
    the first amplitude saturated (|a|/a_max = 1 to six digits, α_1 → ∞: the infimum is not attained at finite α).
  - ε/γ = 10⁻³ gives fit error 1.61·10⁻⁴; ε/γ = 10⁻¹ gives 1.30·10⁻². The fit error is linear in ε/γ.
  - At the experiment's γ = 10⁻², "order ε" understates the bias by a factor of 100.
- *Fix:* "a bias of order ε/γ in the fit, the infimum in the saturation (the class's amplitude at a_max)". If the claim is kept, add it to the
  Status as sketched.

**ID-3 (major). Part D at ε = τ: the 4-seed agreement with the Ornstein–Uhlenbeck spread, and the §8 bullet's "as predicted to 11%",
are chance. The linear prediction misses by 14–23%.**
- *Location:* PARADIGM.md Numerics (D); §0 row; §7 row; §8 bullet ("as predicted to 11%"); RESULTS l. 939–940; PAPERS paper 13 and README
  ("the Ornstein–Uhlenbeck spread ∝ ε^{−1/2}").
- *Evidence* (`mysim.py`, `anasim.py`; 64 seeds, 2·10⁵ slow units, the same scheme):
  - **ε = 10⁻⁴.** Pooled second-half std (0.0170 ± 0.0015, 0.0112 ± 0.0006, 0.0120 ± 0.0006) (bootstrap over seeds), against the predicted
    (0.0139, 0.0099, 0.0105). That is +23%, +14%, +14%, i.e. 2.1, 2.1 and 2.6 standard errors. The full 11-dimensional linear covariance
    gives 0.0138, 0.0098, 0.0105, so the reduction is not at fault: the linearisation is (the spread is 1/5 of the distance to the
    saturation at W_1 = 0.25).
  - **4-seed scatter.** Over 16 disjoint 4-seed subsets the std estimate scatters by 23%, 22% and 12% (relative SD). The repository's
    0.0152, 0.0103, 0.0117 therefore cannot support "to 11%".
  - **ε = 10⁻³.** (0.00429 ± 0.00005, 0.00306 ± 0.00003, 0.00326 ± 0.00003), against the reduced prediction (0.00439, 0.00312, 0.00333),
    2% high, and the full linear (0.00422, 0.00300, 0.00324). Here the 4-seed scatter is only 3–5%, and the repository's values are sound.
  - **The ε^{−1/2} law.** The measured ratio std(10⁻⁴)/std(10⁻³) is 3.96, 3.65, 3.68, against √10 = 3.16. The law holds only for ε ≫ τ.
- *Fix:* run at least 32–64 seeds for D (the batched integrator makes 64 seeds cost about 30 min), report bootstrap errors, and state that
  the linearised spread is accurate to about 2% at ε = 10τ but underestimates by 14–23% at ε = τ. Replace "as predicted to 11%" accordingly.

**ID-4 (minor). The drift term is computable, and the text's reading of the leak's mean shift is not supported by its own data.**
- *Location:* PARADIGM.md (ii), the Status "drift term not computed", Numerics (B) and (D) ("shifts the pinned state by up to 0.008 at
  ε = τ and by at most 0.0005 at ε = 10τ, as O(τ/ε) says"), and the §7 row.
- *Evidence:*
  - **The drift.** `check3.py` gives b_W = (−2.055, +1.580, +0.475)·10⁻⁷ per slow unit: the limit map by LSODA, second differences at
    h = 0.02 and 0.01 agreeing to four digits, about 30 min of CPU.
  - **Part B.** This predicts a mean displacement at 5·10⁴ of (−0.0103, +0.0079, +0.0024). The repository measures (−0.0121 ± 0.0046,
    +0.0112 ± 0.0041, +0.0009 ± 0.0038); my 256 seeds give (−0.0079 ± 0.0023, +0.0072 ± 0.0020, +0.0007 ± 0.0019). At τ = 2.5·10⁻⁵ the
    prediction is (−0.0026, +0.0020, +0.0006) and my 256 seeds give (−0.0018 ± 0.0013, +0.0023 ± 0.0009, −0.0005 ± 0.0010). The drift sends
    class 1, the one closest to its saturation, toward it: this is the "selection of the masses by the dither" in the §8 Open.
  - **Part D, the mean.** κ⁻¹(b + ½F″_leak[C]) predicts the means (0.3294, 0.3359, 0.3346) at ε = τ and (0.33294, 0.33359, 0.33346) at 10τ.
    The shifts are at most 0.0039 and 0.0004, and their ratio is exactly 10. My 64 seeds give (0.3281 ± 0.0013, 0.3370 ± 0.0007,
    0.3348 ± 0.0008) and (0.33298 ± 0.00009, 0.33361 ± 0.00006, 0.33342 ± 0.00006).
  - **The repository's means.** Its 4-seed means at ε = τ are 1.6 and 1.45 between-seed standard errors from 1/3 (SE 0.0050, 0.0044,
    0.0017). The 4-seed means scatter by 0.006, 0.0027, 0.0043, so the "0.008 vs 0.0005" ratio of 16 is noise, not evidence for O(τ/ε).
- *Fix:* compute and state the drift (one cheap computation), compare B's means with it, and replace the D sentence with the predicted shift
  and its error bars.

**ID-5 (minor). Part B's slope is the slope of the mean-square displacement, which contains the squared drift.**
- *Location:* gauge_dither.md B; the PARADIGM Numerics (B) ratios 1.05, 1.53, 0.88; RESULTS l. 935.
- *Evidence* (`rawstats.py`):
  - In the repository's own 64 seeds the drift² term is 10–12% of the slope for classes 1–2. The variance slope (mean removed) is
    (2.69, 2.21, 1.62)·10⁻⁸, so the ratios become **0.95, 1.36, 0.90**.
  - My 256 seeds give 0.94 ± 0.08, 1.26 ± 0.12, 0.99 ± 0.12 at τ = 10⁻⁴ and 1.21 ± 0.11, 0.94 ± 0.10, 1.18 ± 0.11 at τ = 2.5·10⁻⁵. These are
    consistent with 1 at the 2σ level, with no systematic class-2 excess.
- *Fix:* use the variance, not the MSD (or subtract the computed drift), and give the ratio's error, about ±0.1 with 256 seeds.

**ID-6 (minor). The time-scale separation behind "Katzenberger on the time scale 1/τ" is thin at the experiment's parameters; at
ε = 10⁻³ the leak is faster than the slowest normal mode.**
- *Location:* (ii) and (iii) as applied in the Numerics.
- *Evidence* (`check2.py`, `check4.py`):
  - The slowest normal eigenvalue is −3.06·10⁻⁴: mainly bump 2's width coordinate, near s_min (ds/dκ = 0.094). So λ_min/τ ≈ 3.1, and the
    normal relaxation time (3.3·10³) is a third of 1/τ = 10⁴. In the linear OU that mode fluctuates by 0.26 in κ_2 (whitened √(τ/λ) ≈ 0.57).
  - At ε = 10⁻³ the reduced stiffnesses are (7.4, 9.0)·10⁻⁴, above λ_min.
  - The reduction still works, because the soft mode barely couples to the masses: reduced and full linear spreads differ by 4%, and the
    simulation lies between them (ID-3).
- *Fix:* state λ_min/τ ≈ 3 and that ε ≪ λ_min fails at ε = 10⁻³, or give the full linear covariance (cheap) beside the reduced one.

**ID-7 (minor). The proof of (ii) writes Katzenberger's limit with an inconsistent factor.**
- *Location:* PARADIGM.md proof (ii): "dY = ∂Φ(Y)√(2τ)B dW + τ∂²Φ(Y)[BBᵀ]ds, BBᵀ = Σ".
- *What is wrong:* with BBᵀ = Σ = 2τG, the noise covariance would be 2τΣ = 4τ²G and the drift O(τ²).
- *Fix:* BBᵀ = G (that is, dY = ∂Φσ dW + ½∂²Φ[σσᵀ]ds with σσᵀ = Σ = 2τG). The drift is then τΣ_iG_i∂²_iΦ, as I used.

**ID-8 (minor). The boundary statements are partly loose and partly beyond the theorem.**
- *Location:* (ii) "at the boundary … D vanishes there quadratically …, the diffusion stays of order one: the class's excursion into the
  saturation is unbounded in α_k"; the Status "its boundary behaviour as stated".
- *Evidence* (`check6.py`):
  - For K = 3 only the saturating class's row and column vanish: at δ = W_1 − 0.25 = 10⁻³ and 10⁻⁴, D_11 = 5.02·10⁻¹¹ and 5.00·10⁻¹³
    (D_11/δ² = 5.0·10⁻⁵), while D_22 = D_33 → 7.15·10⁻⁹.
  - In α_1 the rate tends to exactly 2τG_α (1.85, 1.998, 2.000·10⁻⁴), which is sharper than "of order one".
  - "Unbounded" needs the drift near the boundary, which is not computed. Katzenberger's theorem holds up to exit from compact subsets of the
    basin, and M is non-compact in α.
- *Fix:* say which rows vanish, give D_α → 2τG_α, and mark the unboundedness as heuristic (consistent with α_1 reaching 16.1 in C) rather
  than proved.

**ID-9 (minor). The tests never tie the metric g to the code's drift, nor the dither's amplitude to τ.**
- *Location:* tests/test_gauge_dither.py.
- *Evidence* (`mutate.py`, `mutate2.py`):
  - **Caught:** dropping a′, m/γ → 1/γ, the leak's sign and factor 2, the leak Hessian's factor m, the factor 2 in the Lyapunov D, Σ = τG,
    and the position metric scaled by 4. The two-class closed forms and the Fisher–Rao test catch these.
  - **Survived:** the noise √(2τ) → √τ in `run_dithered`; noise without G; the position drift without G, without β, or ×10 (each changes the
    true metric while `gauge_diffusion` stays as it is); dropping the Σβ renormalisation; db instead of db/b in the log-β step.
  - The mutant log β + 1 also survives, correctly: it is equivalent.
- *Fix:*
  1. A test comparing `gauge_diffusion` with the spectral projector of the Jacobian of `loop.drift` on Σδβ = 0. In my check the diagonal of D
     agrees to all printed digits, and the projectors to 2·10⁻⁹.
  2. A one-step test of `run_dithered` at γ = 0: Var(ΔΘ) = 2τG·dt.
  3. Σβ = m after a dithered run with γ > 0.

**ID-10 (minor). Literature.**
- *Location:* [S64], LITERATURE_TRANSPORT.md l. 84 and l. 539.
- *Parsons–Rogers:* the first author is **Todd L. Parsons** ("T. L. Parsons"), not "E. T. Parsons" (the file admits "its authors from
  memory"). The paper is published in J. Phys. A: Math. Theor. 50(41):415601, 2017.
- *Correct:* Katzenberger (Ann. Probab. 19(4):1587–1628, 1991), Blanc–Gupta–Valiant–Valiant (COLT 2020, PMLR 125:483–513, checked) and
  Li–Wang–Arora (ICLR 2022, arXiv:2110.06914; the η⁻² time scale) are cited correctly and fairly.
- *Missing:*
  - Constable, McKane, Rogers, "Stochastic dynamics on slow manifolds", J. Phys. A 46 (2013) 295002 (the same projection, with
    noise-induced drift along the slow manifold).
  - Funaki (PTRF 102, 1995), on SPDE motion along a manifold of minimisers.
  - Kurtz, "Averaging for martingale problems and stochastic approximation" (1992).
  - Berglund–Gentz, *Noise-induced phenomena in slow–fast dynamical systems* (2006).
  - On the SGD side: Damian–Ma–Lee (NeurIPS 2021, label-noise SGD and flat minimisers) and HaoChen–Wei–Lee–Ma (COLT 2021, the shape of the
    noise covariance).
  - On the leak: the entropic regulariser on the replicator is the logit / perturbed-best-response (entropy-regularised) replicator, and the
    mutation term of replicator–mutator and Wright–Fisher-with-mutation models.
- *The Gibbs alternative:* dithering the weights too, with Fisher–Rao noise √(2τγβ) (Kimura / Wright–Fisher type), would make the dynamics
  reversible in g. The law on M would then be explicit (the Gibbs density times the volume factor), and the dither's selection of the masses
  computable. This is worth one sentence.
- *"New":* "the masses' diffusion coefficient in the flow's metric" is the obvious metric version of the Euclidean projection lemma, and
  should read as "computed" rather than "new". The application to the amplitude–mass gauge and the entropic pin is a defensible "new for
  this instance".

**ID-11 (minor). "The amplitudes back inside the range (|a|/a_max ≤ 0.89 and 0.82)" depends on the sample size.**
- *Location:* PARADIGM Numerics (D).
- *Evidence:* over my 64 seeds the largest |a|/a_max is 0.992 at ε = 10⁻⁴ and 0.957 at ε = 10⁻³. The leak at ε = τ does not keep every run
  out of the saturation.
- *Fix:* report a quantile, or the value over more seeds.

**ID-12 (cosmetic). Small wording and number points.**
- "|a|/a_max ≤ 0.80 at γ = 0", while the md has 0.8033.
- (C) compares a 10-seed std (relative uncertainty about 20–25%) with √(Ds) without an error. Over 2·10⁵ slow units the drift's displacement
  is 0.031 at s = 1.5·10⁵ for class 1, against √(Ds) = 0.065: not negligible, though the pooled std about the pooled mean removes it.
- The proof of (i) is circular in wording: "nondegenerate on the complement of T where H_L ≻ 0" uses the hypothesis as the conclusion.

**No contradiction with D‴.28 or D‴.34 beyond ID-1.** The gauge formula is D‴.28's closed form W_ka(α_k) = A_k (`exact_rep`), and the
metric is D‴.34's Euclidean–Fisher–Rao one, used consistently. The units are consistent: D and b are per slow time unit, both O(τ), and
κ = O(ε).
