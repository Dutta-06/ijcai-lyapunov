## Referee D3ah (round 35): report on Proposition D‴.46 (a leak in the mirror coordinates)

The algebra holds:
- the Lyapunov identity in (i) and its ultimate bound;
- B(w, 0) = ψ(w) and η(0) = 0;
- the ridge limit, and the LASSO limit in the scaling σλ₁ fixed.

What fails is what the numerics are said to show:
- Part B shows slow convergence, not drift.
- (ii)'s averaged loop leaves out the error filter's lag.
- Part C's gap is that lag in one case and an unconverged run in the other.
- Part D's pruning is one lucky instance.
- The mechanism has close precedents in learning in games that lit_S73 does not cite.

### What I checked and how

**Setup**
- Copied the working tree (`cp -r`) to `…/review35/repo` and worked only there. Every process ran with `nice -n 10` and `OMP_NUM_THREADS=1`, one heavy process at a time.
- Scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/review35/work/`; logs are in `…/review35/logs/` and mutants in `…/review35/mut/`.

**Tests and reproduction**
- `tests/test_mirror_free.py`: 8 passed in 61 s. The D‴.46 test takes 53 s and raises "overflow encountered in dot".
- The committed `mirror_leak.md` regenerates byte-identically from `_raw.npz` through `report()`.
- I reproduced Part C's w_σ (to 1.5e-7) and Part D's equilibria (identical to the module's).

**(i) Lyapunov identity and bound** (`t11_identity.py`)
- Rebuilt ψ by quadrature of η dw, independently of the module. The three-point identity (η(w) − η(w_p))(w* − w) = B(w*, w_p) − B(w*, w) − B(w, w_p) holds to 1.3e-15. The module's B agrees with it.
- B(w*, w(α)) grows without bound at both ends of each curve (≈ 5e15 at α = ±20, for q from −0.124 to 5).

**Part B** (`avglib.py`, `t4_partB.py`, `t4b_partB_pred.py`, `t7_partB_long.py`)
- Built the averaged loop with the error filter L = 1/(s + k_e), computed as the periodic steady state by FFT on the 20π-periodic orbit (H_L = k_e⟨φL[φ]ᵀ⟩, d̄_φ^L = k_e⟨φL[d]⟩).
- Integrated that averaged ODE to 2·10⁵.
- Continued the real σ = 0 loop to t = 18000 with the experiment's own DOP853 settings (221 s).

**Part C** (`t3_partC.py`, `t9_lag.py`, `t10_partC_osc.py`)
- Computed w_σ three ways: with the idealised H/k_e, with the filtered H_L, and from the exact frozen-parameter periodic error (nonlinear in e).
- Integrated the averaged ODE from the experiment's start over T_C = 4·10⁴.
- Measured the real loop's within-period oscillation around its equilibrium.
- Measured the order and the symmetry of the lag correction, also on an orbit that is not odd in t.

**(iii) shapes and Part D**
- `t1_shape.py`: η(0), B(w, 0) − ψ(w), exact versus asymptotic λ₁, and the small-β⁰ remainder of ψ.
- `t2_ridge_irrep.py`: the size of Part A's ridge effect; the irrepresentable condition (IC) of Part D's Gram over all 1760 three-sparse (support, sign) patterns; LASSO recovery along the path from the exact KKT conditions.
- `t5_partD.py`: the leak's equilibrium on IC-failing truths; coordinate-descent LASSO against the exact KKT solution; the σ → 0 limit with a nonzero d̄_φ.
- `t6*_richdyn.py`: rich-start averaged dynamics over 10⁸ time units.
- `t8_solver.py`: an instrumented copy of `leak_equilibrium`.

**Mutation testing**
- `mutate.py`: 25 mutants of `real_rhs_leak`, `leak_equilibrium` and `lasso`, each run against the D‴.46 test.
- `mut_report.py`: 8 mutants of the experiment's bound and report code, evaluated on the committed raw.

**Literature**
- Web searches; arxiv.org itself was not fetched. Items marked "memory" below were not search-verified.

### Verdicts

**(i) Holds as algebra; the "drift" reading of Part B does not.**
- V̇ = −k_ee² + ed − σΣ[B(w*, w) + B(w, w_p) − B(w*, w_p)] has the right sign: checked symbolically, by quadrature (1.3e-15) and in the module (1.2e-10).
- V ≥ 0, and V is proper on each curve for Q > −c. The draft does not state properness, but it is needed to get from bounded V to bounded (α, β).
- V ≤ max{V(0), (d̄²/(2k_e) + σB(w*, w_p))/min(k_e, σ)} is correct, with min(k_e, σ) the right rate.
- When f is not representable (Part B), w* is an arbitrary split and the residual 0.1 sin 4x goes into d. d̄ = 0.3 + 0.1 = 0.4 is valid for that split.
- "Without the leak the coefficients are not bounded a priori" is true as an a priori statement, but Part B does not exhibit it.
  - The orbit is persistently exciting.
  - The σ = 0 loop converges to the bounded, biased least-squares fit.
  - It would not continue: max|w| tends to about 0.34 with time constant 1.5·10⁴ (D3ah-1).

**(ii) Correct for the idealised averaged system, not for the real loop's average.**
- For the idealised system (slaved error) these hold: mirror flow, ∇²F_σ ≻ 0, compact sublevel sets from the entropy tail when Q > −c, uniqueness, and convergence by LaSalle.
- Slow adaptation alone does not give e ≈ (f − f̂ + d)/k_e. The averaged loop has H_L and d̄_φ^L instead (D3ah-2).
  - On a scalar periodic orbit H_L stays symmetric: its first-order lag term is an exact differential and vanishes. So the structure survives, but the equilibrium moves.
- "w_σ → w* as σ → 0" is false whenever d̄_φ ≠ 0 (D3ah-5).
- Part C's agreement is mostly the lag (lazy start) and non-convergence (intermediate start). T = 4·10⁴ is only 4 time constants for σ = 10⁻⁴ (D3ah-3).

**(iii) Limits mostly right; the selection claims are overstated.**
- B(w, 0) = ψ(w) and η(0) = 0 at α = 0 hold exactly, for β⁰ from 10⁻⁸ to 10⁴.
- The ridge limit is correct.
- The LASSO limit is correct when σλ₁ is held fixed (σ ∝ √β⁰), and the √β⁰ rate is confirmed. But the remainder this proof cites from D‴.45 has the wrong form (D3ah-7).
- λ₁ is given only to leading order (D3ah-8).
- "Soft thresholding" is a misnomer for a non-diagonal H (D3ah-9).
- "The leak prunes to the target's kernels" is not generic: on Part D's own Gram the irrepresentable condition holds for only 23.6% of patterns (D3ah-4).
- "A slow leak already regularises" is true of the equilibria. The time to reach them grows like k_eβ⁰^{−1/2} and is never simulated (D3ah-10).

**Novelty: overstated.**
- A leak in the dual coordinates of a continuous-time mirror flow, its regulariser-perturbed rest point, its Bregman or Fenchel-coupling Lyapunov function, and the vanishing-leak selection all exist in the literature (D3ah-6).
- What is new is narrower: the closed-form (α, β) implementation that keeps Q, the β⁰^{−1/2} ℓ¹ weight, and the bound for the adaptive loop.

**Tests and mutations.**
- The tests pass.
- 17 of 25 mutants are killed: the leak's sign; the leak only in α (killed by the Q-kept assertion); the leak only in β (killed by the identity); the leak toward 0; an extra leak in β; the Jacobian without σI, without the metric, with flipped sign, row-scaled, or replaced by the identity; no damping; it_max = 3; and five LASSO mutants.
- 8 survive: the noise sign (R6); noise removed from the adaptation (R7); a damping that accepts 10× residual increases (D2); σk_e on the leak (K1); H without /k_e (K2); the sign of d̄_φ (P1); tolerance 10⁻⁴ (T1); the LASSO threshold without k_e (L1).
- All 8 report-code mutants survive, since nothing checks the experiment's output.

### Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ah-1 | major | **Part B shows slow convergence, not drift.** The orbit is persistently exciting (eig H from 4.6e-3 to 0.69). The averaged σ = 0 loop with the filter reproduces the run's quarterly max\|w\| (0.201, 0.249, 0.275, 0.289 against 0.198, 0.246, 0.273, 0.288). It converges to w* + H_L⁻¹⟨φL[ε + d]⟩: max\|w\| 0.337, ΣB(w*, w) 0.465, slowest time constant 1.5·10⁴, far beyond T_B = 3000. The real loop continued to t = 18000 has block maxima 0.288, 0.315, 0.327, 0.334, 0.337, 0.339 (increments halving) and V going 0.37 → 0.447, approaching that limit. On a periodic orbit ⟨dφ⟩ lies in range H, so the averaged free law cannot drift at all. Part B shows only the leak's price, not its robustness. | Reword "drifts" as "converges slowly (τ ≈ 1.5·10⁴) to the biased fit, \|w\| → 0.34" in numerics (B), BULLET ("where the free masses alone drift"), ROW0, AUDIT and the commit text. Run to ≥ 5τ. To show robustness, add a case without persistent excitation (decaying or degenerate) where σ = 0 actually diverges (exp_mirror_leak.py Part B). |
| D3ah-2 | major | **(ii) is the k_e → ∞ limit, not the slow-adaptation limit.** "In the slow-adaptation limit … the loop is the mirror flow of F_σ" with H = ⟨φφᵀ⟩ and d̄_φ = ⟨dφ⟩ also needs k_e much larger than the orbit's frequencies. The averaged loop actually has H_L = k_e⟨φL[φ]ᵀ⟩ and d̄_φ^L = k_e⟨φL[d]⟩. On a scalar periodic orbit ⟨φφ̇ᵀ⟩ = 0 (exact differential), so H_L = H − ⟨φ̇φ̇ᵀ⟩/k_e² + …, which stays symmetric (asymmetry 1e-13, also on an orbit not odd in t). \|k_eM − H\|/\|H\| = 0.9% at k_e = 20. d̄_φ^L carries a first-order lag. Effects: in Part C (lazy) the equilibrium moves by 8.2e-4 (2.2% of the bias), which is more than the reported 6.9e-4 gap; using H_L the gap falls to 1.6e-4, and to 1.0e-4 with the exact frozen-parameter average. In Part B (σ = 0.01) max\|w_σ\| moves 0.0865 → 0.0754 (13%). ROW0's "(ii) proved" covers only the idealised system. | Restate (ii) with H_L and d̄_φ^L, or as the k_e → ∞ limit with explicit O((ω/k_e)²) and O(ω_d/k_e) errors. Prove the averaging step (an exponentially stable equilibrium gives O(G + γ) closeness). Compare Part C against w_σ(H_L) (D46_draft (ii) and its proof; Part C; ROW0). |
| D3ah-3 | major | **Part C: T is too short, and its "drift" is mislabelled.** For the intermediate start the slowest averaged rate is 1.01e-4 ≈ σ, so T/τ = 4.0. The averaged ODE from the same start is still 3.5e-4 from its own equilibrium over the last quarter (half of the reported 7.1e-4) and drifts 4.36e-4, exactly the real loop's reported 4.36e-4. For the lazy start the averaged ODE is converged (1e-12), and the reported drift 1.3e-4 is within-period oscillation (peak-to-peak 9.5e-4; differences between two instants up to 7.9e-4). The lazy start's fastest averaged time constant (41) is also shorter than the orbit period (63), so "slow adaptation" fails for the excited directions. | Run to ≥ 10/σ (10⁵), or start near w_σ. Measure drift on means over whole periods (20π). Report the gap split into lag, non-convergence and residual (Part C; numerics (C)). |
| D3ah-4 | major | **"The leak prunes to the target's kernels" is a lucky instance.** It holds only under the LASSO's irrepresentable condition ‖H_{S^cS}H_SS⁻¹sign(w*_S)‖_∞ < 1. On Part D's own Gram (coherence 0.86; neighbour correlations 0.76–0.86) that condition holds for 416 of 1760 three-sparse patterns (23.6%), and for 20% of supports with gaps ≥ 4. The chosen truth (2, 6, 9; +, −, +) has IC = 0.771. When the condition fails, the rich-start equilibrium keeps spurious atoms at every μ: (0, 5, 11; +, −, +), IC 1.22, gives support [0, 1, 5, 10, 11] for μ from 10⁻⁴ to 10⁻²; (2, 6, 9; +, +, +), IC 1.14, adds atom 8 for μ ≥ 10⁻³. | State the condition in (iii). Change the README, BULLET, ROW0 and AUDIT ("sparsity-selecting") to "selects the LASSO's support, which is the target's under the irrepresentable condition". Extend Part D over all or random supports and signs and report the recovery fraction. |
| D3ah-5 | major | **"w_σ → w* as σ → 0" is false when d̄_φ ≠ 0** (the case of (i) and Part B). The limit is w* + H⁺d̄_φ, Bregman-selected on null H. On Part D's Gram (β⁰ = 30): \|w_σ − w*\| → 0.026, while \|w_σ − (w* + H⁻¹d̄_φ)\| = 4.9e-7 at σ = 10⁻⁷. The d̄_φ path is never exercised: `dphi` is passed nowhere, and its sign mutation (P1) survives. | Correct (ii). Add a test, and an experiment row with a biased disturbance comparing the loop with w_σ(d̄_φ^L) (tests; Part B or C). |
| D3ah-6 | major | **Novelty and literature: the mechanism is established.** Continuous-time discounted mirror descent (a dual-coordinate leak) has as rest points the equilibria perturbed by the mirror map's regulariser: unique for monotone games, converging to the unperturbed equilibria as ε → 0 (Romano–Pavel, IEEE TAC 66(11):5451–5458, 2021; Gao–Pavel, CDC 2019; arXiv 2011.10682, exponential rates). For the replicator/exponentiated-gradient law it is the entropic, exponentially discounted payoff aggregate with quantal response (Coucheney–Gaujal–Mertikopoulos, MOR 40(3):611–633, 2015). A regulariser toward a reference policy with a strictly decreasing Lyapunov function is Pérolat et al. (ICML 2021). A divergence toward an anchor is APMD (Abe et al., ICML 2024). Explicit regularisation through mirror descent is RMD (Azizan–Lale–Hassibi, arXiv 2202.10788). The σ → 0 selection is viscosity/Tikhonov selection (Alvarez–Bolte–Brahic, SICON 2004; Attouch 1996 from memory). Leakage toward a prior θ₀ is classical. The replicator with a "memory-loss" entropy term (Sato–Akiyama–Farmer 2002; Tuyls et al.; Kianercy–Galstyan 2012) is from memory. lit_S73 and AUDIT cite none of these. "We found no derivation of a leak in a mirror-descent adaptation law's dual coordinates" is misleading. COMID (Duchi et al.) is a proximal method, not a dual leak. | Add these to [S73]. Narrow "new" to: the (α, β) implementation that keeps Q; λ₁ ∝ β⁰^{−1/2}; the adaptive-loop bound. Mark the memory items (lit_S73.md, ROW0, AUDIT). |
| D3ah-7 | minor | **The remainder (iii)'s proof cites from D‴.45 is wrong.** D‴.45(ii) has ψ = λ₁\|w\| + O(\|w\|(1 + \|log\|w\|\|)). In fact r = ψ − λ₁\|w\| ≈ −(3/2)(a_max²/(2Gc))(2c/(a_max²m))^{1/3}\|w\|^{2/3} (coefficient 4.83 here). Measured r/\|w\|^{2/3} = −4.82, −4.80, −4.69, −4.21 at \|w\| = 10⁻⁶ to 10⁻³, independent of β⁰ from 10⁻⁶ to 10⁻¹⁰, while the claimed ratio r/(\|w\|(1+\|log\|w\|\|)) = −32.5 at 10⁻⁶ and grows without bound. The LASSO limit survives (σr → 0). But the remainder is non-uniform (−6.9% of λ₁\|w\| at \|w\| = 10⁻⁶): small active coefficients are under-penalised, and pruning is soft just below the threshold. | Correct D‴.45(ii) (and its manuscript entry) and D‴.46(iii)'s proof. |
| D3ah-8 | minor | **λ₁ is only the leading term.** The exact weight, which the code uses, is πa_maxκ′/(4Gc) with κ′ = ((2c − β⁰)/β⁰)^{1/2}. The asymptotic form is off by 1.0% at β⁰ = 0.04c and by 10.6% at 0.4c. | Give the exact form or the O(β⁰/c) error in (iii). |
| D3ah-9 | minor | **"Soft thresholding" is the wrong name.** \|(H(w* − w))_i\|/k_e < σλ₁,i is the LASSO's KKT condition, implicit in the solution w. It is soft thresholding of w* only when H is diagonal. | Reword (iii). |
| D3ah-10 | minor | **The rich (LASSO) regime is never run dynamically.** Part C's starts have β⁰/c = 120 and 12, both on the lazy side; ROW7 calls the second a "small-mass start". In the averaged loop (Part D's Gram, k_e = 1, β⁰ = 10⁻⁶) the equilibrium is reached only at t ≈ 2.7–3.4·10⁵, and at 0.9–3·10⁶ for a truth that fails the condition. That time scales like k_eλ₁ ∝ k_eβ⁰^{−1/2}. "A slow leak already regularises" therefore describes equilibria only. | Add a rich-start run with its times; fix ROW7; state the escape time in (iii). |
| D3ah-11 | minor | **Part B's bound check has no power.** The bound is loose (1.74 against 0.085; 0.30 against 0.13) and nothing asserts it. Five single mutations of the bound formula all keep the printed bound above V: min → max, d̄ without the bias, d̄²/k_e, no σB(w*, w_p), and no V(0) (the last leaves the .md byte-identical). Only the double mutation "no V(0) and max" shows. | Assert V ≤ bound in the experiment and a test, and test it where it binds (large V(0), worst-case d) (exp_mirror_leak.py). |
| D3ah-12 | minor | **Test gaps: 8 of 25 mutants survive.** The noise path (R6, R7). The k_e scalings (K1, K2, L1), because every test uses k_e = 1. d̄_φ's sign (P1). Tolerance 10⁻¹³ → 10⁻⁴ (T1) and a damping that accepts 10× residual increases (D2): both pass in 2–5 s instead of 80 s. | Test with k_e ≠ 1, a nonzero d̄_φ and noise; assert the solver's residual (tests/test_mirror_free.py). |
| D3ah-13 | minor | **`leak_equilibrium` fails silently.** It returns no residual or iteration count and never warns. At β⁰ = 10⁻⁸ (Part A and the test) it hits it_max = 200 with 183 failed line searches and ends at residual 1.8e-13 > tol; this is where the test's 53 s go. On a failed line search it returns the rejected trial point while η is moved by the halved step. Trial steps overflow, and `alpha_of_eta` clamps to ±340 without notice. All other calls converge (3–20 iterations, residual ≤ 7e-15). | Return (w, α, it, ‖R‖) and assert in the experiment; keep the best iterate; cap trial steps in η; set tol relative to the attainable floor (mirror_free.py). |
| D3ah-14 | minor | **(i)'s hypotheses are incomplete.** It needs w*_i in each curve's range (D‴.45(v) had this) and V proper on each curve. "f = Σw*φ on the orbit" must hold along the actual trajectory (e(0) = 0.3 in Part B). For non-representable f, w* is a free choice and the bound should be minimised over it. | Add these to (i) and to Part B's text. |
| D3ah-15 | minor | **The simulated law is not the implementable one.** It computes η from the nominal q, whereas the law on (α, β) computes Q from the state. The sampled-data case, where Q drifts (D‴.45(iv)(b)) and η_p moves with it, is never run. | Compute q = Q(α, β) in `real_rhs_leak`; add an Euler run, or restrict "implementable" to continuous time. |
| D3ah-16 | minor | **Part A's ridge rows test only linear response.** The ridge parameter is 1.2e-9 (β⁰ = 10⁴) and 1.2e-5 (10²) against eig H of 0.15–2.4, so the ridge moves w by only 1e-9 to 4e-5 relative to w*. The 2.3e-13 agreement checks ψ″(0) only, not the ridge "up to the saturation's wall". Part D's lazy rows are the real ridge test. | Choose σ for an O(1) ridge effect, and add a case near w_⊓ (Part A). |
| D3ah-17 | cosmetic | Part D's ridge "support" (6–12 atoms) is a threshold artefact. Part D's distances are absolute and Part A's relative, which the paragraph does not say. "V rising from 0.18" quotes V(0). `np.where` in `eta` evaluates log1p on masked elements, producing RuntimeWarnings. | Edit numerics (B) and (D); evaluate the stable branch only where \|t\| < 0.9. |

### Lesson

In every one of these cases the experiment's own parameters put the effect it neglects at the size of the effect it reports:
- In Part C the error filter's lag (0.9% of H, 2.2% of the bias) and the slowest time constant (T/τ = 4) are as large as the percent-level agreement claimed for an idealised model.
- In Part B the slow excited mode (τ = 1.5·10⁴ against T = 3000) was read as divergence.
- In Part D one well-separated truth was offered as generic, while the same Gram fails the irrepresentable condition for three patterns in four.

Before a number is called confirmation, compute the neglected quantity for the same instance. Here that means:
- the filtered Gram;
- the run length in units of the slowest averaged time constant;
- the excitation's smallest eigenvalue before saying "drift";
- the irrepresentable condition over all supports before saying "prunes".

Each would have taken seconds, and each changes what the paragraph may claim.
