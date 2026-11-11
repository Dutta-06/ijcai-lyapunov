# Referee report D3o (round 20) — PARADIGM §5d: Proposition D‴.28 (statement (i)–(iii), proof, Status, *Numerics*), the §0 and §7 rows "| D‴.28 |", the §8 bullet "Done (Proposition D‴.28)", `results/RESULTS_PARADIGM.md` (Proposition D‴.28), the LITERATURE audit bullet, `experiments/exp_gauge_floor.py` → `results/chaos/gauge_floor.md` (+ raw), `tests/test_gauge_floor.py`; and Corollary D‴.26(c) as revised in round 19 (the statement, the proof, the Status, the *Numerics* Part C), `src/core/adaptive_floquet.CollapsedVariational` and `whitened_jacobian`, `experiments/exp_adaptive_floquet.py` Part C, the three new tests of `tests/test_adaptive_floquet.py`

Branch `claude/festive-bardeen-eopej2` at ee2ed27, HEAD throughout the review. The concurrent session added untracked files
(`experiments/exp_eiv_remedy.py`, `experiments/exp_measurement_noise.py`, `results/chaos/measurement_noise_raw.npz`,
`src/core/measurement_noise.py`, `tests/test_eiv_remedy.py`, `tests/test_measurement_noise.py`, and their bytecode). At the start no
tracked file differed from ee2ed27 (`git diff --stat ee2ed27` was empty). By the end it had uncommitted edits to five tracked files
(PARADIGM, RESULTS, LITERATURE_TRANSPORT, PAPERS, README; 154 insertions, 8 deletions). They add new material after Proposition D‴.29,
with its rows and §8 bullets, and extend the header's file lists. They delete or rewrite none of the reviewed passages, but they shift
the working tree's line numbers. The committed version was reviewed throughout. Every script that uses committed code imports the
ee2ed27 tree, extracted with `git archive` into the scratch directory. The two independent modules (`ginst.py`, `var26.py`) import no
repository code.

Line numbers are ee2ed27's:
- PARADIGM, D‴.28: l. 3242–3315 — statement l. 3242–3269 ((i) l. 3248, (ii) l. 3253, (iii) l. 3262), proof l. 3271–3281, Status
  l. 3283–3293, *Numerics* l. 3295–3315; §0 row l. 129, §7 row l. 5822, §8 bullet l. 6071–6077.
- PARADIGM, D‴.26: l. 3064–3159 — (c) l. 3076–3086, the stability consequence l. 3088–3091, proof l. 3093–3121, Status l. 3123–3139,
  *Numerics* l. 3141–3159; its §0 row l. 127, §7 row l. 5820, §8 bullet l. 6066–6070.
- RESULTS: D‴.26 l. 617–637, D‴.28 l. 662–682.
- LITERATURE_TRANSPORT: D‴.26 l. 387–399, D‴.28 l. 410–416.

No repository file was modified, and no git command touching the index or the tree was run. Every Python process ran on one thread
(OMP/OPENBLAS/MKL_NUM_THREADS = 1). pytest ran on the scratch tree with `-p no:cacheprovider` and PYTHONDONTWRITEBYTECODE=1. The scratch
scripts ran with `python -I`, which ignores PYTHONDONTWRITEBYTECODE, so their bytecode went to the scratch directories; none went to
the repository. The two `.pyc` files that appeared in the repository's `src/core/__pycache__` during the review (`adaptive_floquet`,
22:43; `measurement_noise`, 22:46) come from the concurrent session: its `measurement_noise.py` imports `adaptive_floquet`, and at that
time no script of mine imported it.

`--from-raw` on the scratch tree rewrites both md files byte-identically: `gauge_floor.md` md5 `446c9252…` and `adaptive_floquet.md` md5
`d689083f…`, both equal to `git show ee2ed27:…`.

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3o/scripts/`:

| script | what it does |
|---|---|
| `ginst.py` (module) | independent re-implementation of D‴.28's instance (three bumps + ε·fourth): features with analytic first and second parameter derivatives, L_w with analytic gradient and Hessian, R_k, Newton for critical points and for the fixed prior's weighted stationary point, the gauge, the continuation in ε, the averaged nearest-atom drift with a fixed or an adaptive prior |
| `check_A.py` | Part A recomputed independently against the raw record. The gauge is checked non-tautologically: Newton continued in w from ϑ^b(q), and plain Newton from ϑ^b(q). The sensitivities come from independently solved stationary points |
| `check_branch.py` | the branch ϑ^b(q) continued to the saturation (ε = 0.69); 200 multi-start minimisations at ε = 0.1, 0.3, 0.5 |
| `check_Wg_prob.py` | P(w^m ∈ W_g) for m = 15, 60, 240, 510 per ε; the count vectors in W_g |
| `check_vector_amp.py` | (iii) with vector amplitudes (a planar orbit, two outputs; own code) |
| `check_margins_dyn.py` (+ `check_margins_dyn_helpers.py`) | the separation eigenvalues of the dynamics against σ₀ − (q_k/w_k)λ⁺_k (independent and committed drifts, m = 15); the fixed prior's D‴.5 margins at every Part B sample against the end diameters; a split under the adaptive prior at σ₀ = 0.005σ_crit |
| `check_window_margin.py` | a leak inside the window at w = q with a negative adaptive margin at a sample of W_g |
| `recompute_B.py` | every Part B number from the raw records, independently of `report()`; the readings and time series |
| `check_B_outside.py` | five m = 15 Part B systems re-run bitwise with the committed code and continued to s = 6·10⁴; the minimum of L_w on the closed amplitude box, on the branch and over all assignments |
| `var26.py` (module) | independent implementation of D‴.26(c)'s collapsed variational equation: Ĵ(t), M(t) (optionally with a(t)), the monodromy by DOP853 and by a fourth-order Magnus integrator, Floquet exponents by Hill's method, the slaving and second-order averaging (own FFT code), c_j |
| `check_C1.py` | Part C's exponents against the raw record (all k_e, both references, both η); DOP853 at three tolerances; Magnus; Hill's method and the error channel |
| `check_C2.py` | the averaged generators' convergence in N and c_j against the raw record; ω₁ by Parseval; parity fits; a derivative mismatch a(t); the gain scaled (positive c_j); the generator against the committed real loop at k_e = 25, 400 |
| `check_C3.py` | Magnus at N = 4000–32000 (its rounding floor); G × 30 by the monodromy |
| `check_C4.py` | DOP853's noise floor at k_e = 200, 400 (five start times and a product of two half periods) |
| `check_A1_part.py` | Ā₁'s share of c_j on the asymmetric reference; the first-order shift a symmetric matrix of its size would cause |
| `recompute_C_md.py` | every Part C number of the md, the paragraph, the rows and RESULTS from the raw records |
| `check_tests_A.py`, `check_tests_B.py` | the asserted quantities of both test files, and mutants |

Also in `ref_D3o/`: `tests.out` (8 passed, 123.85 s), `from_raw_gauge.out`, `from_raw_floquet.out`, `md5_before.out`, `md5_after.out`, and
the ee2ed27 tree (`tree/`). The experiments were not rerun in full; single systems were reproduced bitwise (confirmation 7).

## Verdicts

| object | verdict |
|---|---|
| **D‴.28(i)** | **holds, exactly, as a statement about stationary points** (for scalar and for vector amplitudes). G_w is a diffeomorphism of {α : \|q_ka(α_k)/w_k\| < a_max} onto its image, and L_w∘G_w = L_q. At a critical point the Hessian's signature is kept (Sylvester), so nondegeneracy is kept. W_g excludes empty cells (w_k > q_k\|a^b_k\|/a_max > 0). The atoms are distinct, since two kernels with equal centre and width would make the Hessian singular. With the prior at its atoms the replicated configuration is stationary for the nearest-atom adaptive dynamics (drift 4.2·10⁻¹⁷, `check_margins_dyn.py`). Confirmed non-tautologically: Newton continued in w from ϑ^b(q) lands on G_wϑ^b(q) to ≤ 1.2·10⁻¹², f̂ to ≤ 5.8·10⁻¹⁴. The md's own check is tautological (D3o-1). "No function-level floor whatever the approximation error" needs two unstated hypotheses (a minimum, positive adaptive margins) and holds only on W_g, which shrinks with the error and ends at ε ≈ 0.69 (D3o-3). |
| **(ii)** | **holds**, re-derived: Φ = −Ju, B = −Au, (A + σ₀D)δΘ = Auδw, δf̂ = σ₀Ĵ(Ã + σ₀I)⁻¹g̃δw with g̃ = −D^{1/2}u, Φ = Ĵg̃. D's dependence on w drops at q since Θ − ϑ⁰ = 0 there. Ã ≻ 0 on the branch (smallest eigenvalue 4.6·10⁻⁵–3.1·10⁻⁴), so the fixed prior's point at q is a strict minimum. ϑ^b(q) is the lowest of 200 multi-start minima at ε = 0.1, 0.3, 0.5. Independent fixed-prior stationary points reproduce F to 1.9·10⁻⁹ (with the analytic gauge derivative u); the Gauss–Newton form misses it by 6.99·10⁻⁴–5.64·10⁻³, every md digit. |
| **(iii)** | **holds for a scalar amplitude per kernel.** The separation eigenvalues of the dynamics equal −σ₀ + (q_k/w_k)eig(G^{1/2}R_k(q)G^{1/2}) to 5·10⁻¹² (independent drift, m = 15, η = σ₀ and ∞). The separation subspace is invariant and the collective block stable. "Can split" is confirmed by runs (D3o-3). **Fails for vector amplitudes**, i.e. the saturated network with n ≥ 2 outputs: the amplitude–centre block survives and the scaling breaks (D3o-2). The admissible set {w_k > q_kλ⁺_k/σ₀} presumes ϑ^b a minimum (D3o-3). |
| Proof (D‴.28) | **correct** for a scalar amplitude: (i) by L_w∘G_w = L_q, (ii) by the implicit-function derivation, (iii) by ⟨rΦ̃⟩ = a(α)Ψ(γ, κ). The last presumes a scalar a, unstated (D3o-2). |
| *Status* (D‴.28) | **accurate** on what is proved. "For real plants" is overstated (D3o-2, D3o-3). The novelty is fair in substance but uncited (D3o-10). |
| *Numerics* (D‴.28), `gauge_floor.md` | **every number reproduces**, independently (confirmations 2, 7). The gauge check and the adaptive sensitivity test nothing beyond an identity (D3o-1). The outside-W_g reading is wrong (D3o-4). Rounding and wording slips (D3o-8). |
| margins table | **numbers hold**. It compares the code's curvature at G_wϑ^b(q) with the formula, not with the dynamics; the dynamics agree as well (confirmation 5). "None negative" is a property of σ₀ = 0.15σ_crit (D3o-3). |
| Part B | **holds.** Every number reproduces from the raw records, and five systems were re-run bitwise. The fixed prior splits exactly at its D‴.5 saddles. "Constant from s = 600" holds to within 8.3% (D3o-8). The near-boundary sample does converge to its gauge transform. Outside W_g the "kept error" quoted is a transient (D3o-4). |
| `exp_gauge_floor.py` | **holds**; regeneration is byte-identical. A docstring sign slip (D3o-8). |
| `tests/test_gauge_floor.py` | **passes (4/4)**. Test 4 discriminates a fixed-prior mutant. Two assertions are tautological, one tolerance is fragile, and only η = ∞ is tested (D3o-9). |
| §0, §7 rows, §8 bullet, RESULTS, audit (D‴.28) | **numbers hold**, except "7·10⁻⁴–3·10⁻³" in the rows (D3o-8). Wording: D3o-1, D3o-3, D3o-4, D3o-10. |
| **D‴.26(c) as revised** | **holds as a formal expansion, and its structure is right.** Re-derived: L₁ = −J̄, L₂ = J̄̇ − aJ̄, L₃ = −J̄̈ + J̄C − J̄J̄ᵀJ̄; B₀, B₁, B₂; Ā₀, Ā₁, Ā₂ with zero-mean U's; the skewness of Ā₁ (both parts live on the y block, where the metric is I); M(t)ᵀ = DM(t)D⁻¹ and the equal multipliers of M(−t); the parity count (#J̄ even, so d + i ≡ k mod 2); c_j's formula; Liouville. Numerically even to O(k_e⁻⁴) on both references. The sign of c_j is not determined (D3o-5). The filter's part at finite η is mislabelled (D3o-6). The parity's summary names an incomplete reason (D3o-7). |
| the a(t) claim | **holds.** With a(t) = 0.5 + 0.3 cos t the reversed system has the same multipliers (4.7·10⁻¹¹ relative), yet odd orders appear: b₁/λ = +0.44 in the slowest mode and +0.40 in the fastest. |
| Liouville's formula, the error channel | **holds.** Hill's method gives the error channel's exponent directly — 24.98083504 at k_e = 25 and 99.99509993 at 100 — equal to Liouville's. (err − (k_e − tr H_L/k_e))·k_e³ → 7.79 (sym) and 10.73 (asym): O(k_e⁻³). |
| *Status* (D‴.26) | **accurate**: "a formal asymptotic expansion whose structure is proved". The Checked list matches the tests, whose limits are D3o-9. |
| *Numerics* (D‴.26, Part C), `adaptive_floquet.md` Part C | **every number reproduces** (confirmations 9, 13). The resolution is adequate: DOP853's noise floor is 2.6·10⁻⁹ relative in the slowest rate at k_e = 400, so the 4.9·10⁻⁷ residual is resolved to 0.5%. "Every exponent below" is an instance property (D3o-5). The "−ω₁²" label is wrong at η = σ₀ (D3o-6). |
| `CollapsedVariational`, `whitened_jacobian` | **hold.** Ĵ(t) equals the independent implementation's to 1.1·10⁻¹⁶. The generator equals the committed real loop's central-difference linearisation to 9.0·10⁻¹¹ at k_e = 25 and 400. The averaged generators are converged at N = 256 (≤ 1.7·10⁻¹³ to N = 2048) and reproduced by independent code to 1.9·10⁻¹⁴. `second_order`'s filter part is computed correctly; only its label in the report is wrong (D3o-6). |
| `exp_adaptive_floquet.py` Part C | **holds**; regeneration is byte-identical. |
| the three new tests | **pass, not vacuous**: the mutants A₂ = ⟨B₂⟩ and "no filter term" fail. But an Ā₁ = 0 mutant passes test 4; the metric is never tested; `var < p1` asserts an instance property; the Liouville tolerance is 17–130 times loose (D3o-9). |
| novelty (both audits) | the gauge's mechanism is fair; "exact margin formula" needs "scalar amplitude"; the even expansion is the classical time-symmetry argument (D3o-10). |

## Confirmed issues

**D3o-1 (minor; D‴.28 *Numerics* l. 3298–3301 "At six weight vectors of W_g Newton's method started at G_wϑ^b(q) stays there (distance
0.0) with f̂ unchanged to ≤ 7.8·10⁻¹⁶, and the adaptive prior's function-level sensitivity is ≤ 2.3·10⁻¹¹ of the fixed prior's (the
finite-difference level)"; md Part A columns "gauge: distance, Δf̂" and "‖F_ad‖/‖F_fixed‖"; §0 row l. 129 "the gauge transform is
Newton's critical point to rounding with f̂ unchanged to 10⁻¹⁵, the adaptive prior's sensitivity at the finite-difference level (10⁻¹¹
relative)"; RESULTS l. 667–668; `exp_gauge_floor.py` l. 153–155 and l. 162; `test_gauge_floor.py` l. 59, l. 71 and l. 87).** Both checks
are tautological.
- *The distance.* `floor_bias._exact_weighted_atoms` compares the gradient with tol = 10⁻¹³ before its first step. Started at
  G_wϑ^b(q), whose gradient is ~10⁻¹⁶, it returns its starting point, so "distance 0.0" holds by construction. It is exactly 0.0 in all
  36 Part A checks, and test 2's Newton also moved by exactly 0.0 (`check_tests_A.py`).
- *The adaptive sensitivity.* It is the central difference of f̂_w(G_wϑ^b(q)), which equals f̂_q(ϑ^b(q)) for every ϑ by the gauge's
  definition. It therefore measures the rounding of tanh∘artanh divided by 2·10⁻⁵, and is exactly 0.0 at ε = 0 and 0.5. The same holds
  for Δf̂ and for test 1's f̂ assertion.

What these checks establish is only ‖∇L_w(G_wϑ^b)‖_∞ < 10⁻¹³, which is not printed. The content of (i) that the dynamics needs is that
the critical point on ϑ^b(q)'s branch at the weights w is G_wϑ^b(q). Checked non-tautologically (`check_A.py`):
- Newton continued in w from ϑ^b(q), in 20 steps, lands on G_wϑ^b(q) to ≤ 1.2·10⁻¹² with f̂ unchanged to ≤ 5.8·10⁻¹⁴, at six weight
  vectors 0.01 inside W_g for every ε.
- Critical points continued from ϑ^b(q) at q ± 10⁻⁵(e_j − e₃) give an adaptive sensitivity of 1.6·10⁻¹¹–4.0·10⁻¹¹ of the fixed prior's.
  The conclusion holds.
- Plain Newton from ϑ^b(q) at the same weights, without continuation, ends at another critical point of L_w in 1–2 of the 6 cases for
  every ε (distances 1.5–66; one singular Hessian at ε = 0.2). "Newton's critical point is the gauge transform" therefore depends on
  where Newton starts.

*Fix:*
- Print ‖∇L_w(G_wϑ^b(q))‖, and compute the adaptive sensitivity from critical points started at ϑ^b(q) (or continued from q).
- In test 2, start Newton at `tb`.
- Drop "stays there (distance 0.0)" and "the finite-difference level" as evidence, in the paragraph, the §0 row and RESULTS.

**D3o-2 (minor; (iii) l. 3262–3269; proof l. 3277–3281 "With Ψ(γ, κ) := ⟨rφ(·; c(γ), s(κ))⟩_ρ, ⟨rΦ̃(·;ϑ)⟩_ρ = a(α)Ψ"; Status l. 3286–3289
"the exact scaling of the residual curvature under the gauge"; §0 row l. 129; §8 l. 6074–6075; audit l. 414–415 "the exact scaling
R_k(G_wϑ) = (q_k/w_k)R_k(ϑ) of the residual curvature at a critical point, which turns the adaptive prior's splitting margins into a
closed formula").** (iii) needs a scalar amplitude.

The proposition covers every network in which "the masses enter only through the products w_ka(α_k), as they do for the saturated
network". For n ≥ 2 outputs the saturated network's amplitude is a vector a(α) ∈ ℝⁿ, saturated componentwise, and
⟨rᵀΦ̃⟩ = Σ_ra_r(α_r)Ψ_r(γ, κ) is not of the form a(α)Ψ. At a critical point:
- stationarity in α_r gives Ψ_r = 0 for every r;
- stationarity in γ gives only Σ_ra_r∂_γΨ_r = 0, not ∂_γΨ_r = 0;
- so the amplitude–centre block a_r′(α_r)∂_γΨ_r survives, and a_r′ = 1 − (a_r/a_max)² does not scale like q_k/w_k under the gauge.

Counterexample (`check_vector_amp.py`: orbit (2 sin t, 2 sin 2t), three vector bumps plus 0.3 times a fourth, K₀ = 3, a critical point
of L_q with gradient 2.2·10⁻¹⁷ and fit error 9.0%):
- The amplitude–centre entries of R_k reach 3.0·10⁻³ (max\|R\| = 4.7·10⁻³). The (α, α) block is 3.8·10⁻¹⁸.
- At w = (0.42, 0.30, 0.28) the gauge transform is a critical point (2.2·10⁻¹⁷) with the same function (3.3·10⁻¹⁶), so (i) holds.
- But \|R_k(G_wϑ) − (q_k/w_k)R_k(ϑ)\| = 6.2·10⁻⁴, 13% of max\|R\|. Only the centre–width block scales, to 1.5·10⁻¹³.
- λ_max(G^{1/2}R_kG^{1/2}) at G_wϑ misses (q_k/w_k)λ⁺_k by +0.1%, −11.6% and −16.8% in the three classes.
- No zero eigenvalue is forced: the smallest eigenvalues are −3.5·10⁻⁴, −2.9·10⁻³ and −1.3·10⁻³.

The instance has one output and D‴.29's servo one torque, so no number of the document changes.

*Fix:* add "with a scalar amplitude per kernel (one output)" to (iii), the Status, the rows, §8 and the audit. For vector amplitudes,
say that only the centre–width block scales by q_k/w_k, while in the physical amplitude coordinates the amplitude–centre block ∂_γΨ_r is
gauge-invariant, so no closed margin formula follows.

**D3o-3 (minor; (i) l. 3248–3252 "the function-level distance between two systems at their approximants vanishes whatever the
approximation error, and the fit error is ‖r_b‖_ρ at every sample of W_g"; the hypothesis "a nondegenerate critical point" l. 3242;
(iii)'s admissible set l. 3268–3269; §0 row l. 129 "no function-level floor whatever the approximation error"; §8 l. 6072–6074 "so the
function-level floor vanishes on the gauge set W_g whatever the approximation error"; Status l. 3286–3287 "for real plants, where no
exact representation exists").** The floor is the stationary E‖f̂_m − f̂_M‖² of two coupled systems (D‴.9, D‴.24). (i) proves that their
stationary points carry the same function. A zero floor needs three things more.
- *ϑ^b(q) a local minimum.* A nondegenerate saddle of L_q maps to a saddle of L_w, since the signature is kept, and then the adaptive
  prior's collective block is unstable. (iii)'s admissible set {w ∈ W_g : w_k > q_kλ⁺_k/σ₀} also presumes a minimum. On the instance
  ϑ^b(q) is one (confirmation 3).
- *Positive adaptive margins.* (iii) itself says that a class with w_k < q_kλ⁺_k/σ₀ can split, yet (i), §0 and §8 state the floor
  removal on all of W_g. This happens inside W_g (`check_window_margin.py`, `check_margins_dyn.py` (3)):
  - At ε = 0.5 with σ₀′ = 0.0165σ_crit, every margin at w = q is positive (5.7·10⁻⁴, 5.8·10⁻⁵, 7.7·10⁻⁴), so the leak is in the window
    there. Yet the sample (11, 9, 10)/30 ∈ W_g has the class-2 adaptive margin −2.17·10⁻⁵, and the dynamics' largest separation
    eigenvalue is +2.166·10⁻⁵, as predicted.
  - At σ₀′ = 0.005σ_crit, at w = q ∈ W_g (class-2 margin −4.84·10⁻⁴), a class of the adaptive prior's system (η = ∞) splits: its
    diameter grows from 5·10⁻³ to 1.1 by s = 2·10⁴. The fit error leaves the approximant's 2.5275·10⁻² and is 2.129·10⁻² at s = 4·10⁴.
    The function is no longer f̂_b.
- *A sample in W_g.* W_g shrinks with the approximation error, and at the sizes used its complement is not small (`check_Wg_prob.py`).
  - P(w^m ∈ W_g) = 0.52, 0.87, 0.997 at m = 15, 60, 240 for ε = 0.05–0.3, and 0.26, 0.71, 0.98 at ε = 0.5.
  - On the instance family the branch's class-2 amplitude reaches the saturation at ε ≈ 0.69 (`check_branch.py`): max\|a\|/a_max =
    0.977, 0.9997, 1.000 at ε = 0.65, 0.68, 0.69. There the smallest eigenvalue of Ã is 1.7·10⁻¹⁷, λ⁺/σ₀ = 0.27, and W_g's lower bound
    for class 2 reaches q₂ = 1/3: w = q itself leaves W_g.
  - "Whatever the approximation error" thus ends at a fit error of 3.7% on this family. The Status's exclusion of critical points at
    the boundary covers that endpoint, but not the shrinking before it.

On the instance at σ₀ = 0.15σ_crit and ε ≤ 0.5 the margin condition never binds. The thresholds q_kλ⁺_k/σ₀ are ≤ 0.034, against W_g's
lower bounds ≥ 0.175, which is why the margins table shows no negative margin.

*Fix:*
- In (i), §0 and §8, write "at the stationary points; the floor vanishes on {w ∈ W_g : w_k > q_kλ⁺_k/σ₀}, for a strict local minimum
  ϑ^b(q), when the systems converge (numerical)".
- Add "a minimum" to the hypothesis of (iii)'s admissible set.
- Quote P(W_g) at the sizes used, and the end of the instance's branch at ε ≈ 0.69.

**D3o-4 (minor; *Numerics* l. 3312–3315 "At the m = 15 sample (0.533, 0.267, 0.200), outside W_g (margins −0.0016 and −0.0097 at ε = 0.1
and 0.3), no gauge transform exists and the adaptive priors keep an error (fit errors 2.2·10⁻² and 1.6·10⁻² at ε = 0.1, against
6.84·10⁻³)"; §7 row l. 5822 "an error kept at the sample outside W_g"; RESULTS l. 676–677 "at the sample outside W_g an error remains";
the md's near-boundary table).** The values quoted are transients, not the kept error.

The systems were re-run bitwise with the committed code (max difference from the raw fit records 0.0) and continued to s = 6·10⁴
(`check_B_outside.py`):
- *ε = 0.1, seed 3 (margin −0.0016).* Between s = 6000 and 6·10⁴ the fit error falls from 2.17·10⁻² to 1.05·10⁻² (η = σ₀) and from
  1.65·10⁻² to 9.33·10⁻³ (η = ∞). It decreases at every record, while the class-3 amplitude coordinate runs off (α₃ = 8.5 → 11.2 at
  η = ∞). The limit is the minimum of L_w over the closed amplitude box on the branch: fit error 7.575·10⁻³, 1.107 times the
  approximant's. The kept error is therefore 11%, not the factors 3.2 and 2.4 quoted.
- *ε = 0.3, seed 3 (margin −0.0097).* At η = ∞ the fit error falls from 3.30·10⁻² to 2.78·10⁻², toward the box minimum 2.688·10⁻²
  (1.586 times the approximant's).
- *The sample 6·10⁻⁴ inside W_g* does converge to its gauge transform. Its excess is 0.35 and 0.18 at s = 6000 and 0.038 and 0.019 at
  s = 6·10⁴ (η = σ₀, ∞), and its box minimum equals the approximant's. "The approach is slow" is right. The two cases differ in their
  limits, not in how they look at s = 6000.
- *"No gauge transform exists" is true on the branch only.* At (0.533, 0.267, 0.200), with classes 2 and 3 exchanged, the global box
  minimum reaches the approximant's fit error exactly (amplitudes 1.88, 3.02, −3.56). A permuted gauge transform exists, which the
  warm-cloud dynamics does not reach. D3n-10 made the same correction for D‴.24.

*Fix:* "Outside W_g the systems approach, algebraically (the amplitude coordinate running off like ln s), the minimum on the branch at
the saturation: 1.11 and 1.59 times the approximant's fit error at ε = 0.1 and 0.3 (still 2.2·10⁻² at ε = 0.1, η = σ₀, and 3.3·10⁻² at ε = 0.3, η = ∞,
at s = 6000). A gauge
transform with two classes exchanged exists, but the dynamics does not reach it." Apply the same in the §7 row and RESULTS.

**D3o-5 (minor; (c) l. 3084–3086 "c_j contains the error filter's attenuation −v_jᵀ⟨Ĵ̇ᵀĴ̇⟩v_j ≤ 0, which slows every mode"; §0 row l. 127
"c_j, explicit, contains the error filter's attenuation, which slows every mode"; §8 l. 6069 "the error filter slowing every mode";
*Numerics* l. 3149–3150 "every exponent lies below its first-order prediction at every k_e"; RESULTS l. 626–629 "every exponent below its
first-order prediction … The second-order coefficients are all negative"; `test_adaptive_floquet.py` l. 108 `assert np.all(var < p1)  #
the correction slows every mode (the error filter)`).** (c) does not determine the sign of c_j. Its explicit form has terms of both
signs:
- *The slaving's term.* B₂'s −(J̄ᵀJ̄)², from L₃'s −J̄J̄ᵀJ̄, speeds every mode up. For a constant regressor it is the whole correction:
  the slow root of λ² + k_eλ + j² = 0 gives rate·k_e = j² + j⁴/k_e² + ….
- *The level repulsion.* −Σ_iK_jiK_ij/(λ_i − λ_j) = Σ_iK_ij²/(λ_i − λ_j) is ≥ 0 for the slowest mode.

The filter's part dominates on the instance because ω_j² ≫ \|Ĵ\|² there (ω₁² = 103 against tr H_L = 0.49), so "every c_j < 0" is a
property of the instance, not of (c). Scaling the adaptation gain G, a free design parameter that changes nothing else, reverses it
(`check_C2.py` (5), `check_C3.py`, η = ∞):
- *G × 30.* c₈/λ₈ = +4.17 and c₉/λ₉ = +12.65. The exact exponents (monodromy) lie *above* their first-order predictions by 4.17·10⁻⁴
  and 1.270·10⁻³ at k_e = 100, and by 1.04·10⁻⁴ and 3.17·10⁻⁴ at k_e = 200. c_j/(λ_jk_e²) predicts 4.17·10⁻⁴, 1.265·10⁻³ and
  1.04·10⁻⁴, 3.16·10⁻⁴.
- *G × 300.* c₉/λ₉ = +1785, and the exponent lies 4.7% above at k_e = 200 (predicted 4.5%).
- The slowest mode stays below in these runs. Whether the stability rate is over- or underestimated is therefore also an observation.

*Fix:*
- In (c): "c_j = (the filter's part, ≤ 0) + (the slaving's −(J̄ᵀJ̄)² part, ≥ 0) + (the averaging terms, of either sign); on the instance
  every c_j < 0".
- Add "here" to RESULTS and §8.
- Comment the test's `var < p1` as an instance property.

**D3o-6 (cosmetic; *Numerics* l. 3152–3154 "of which the error filter's part is −102.1 = −ω₁², ω₁ = 10.1 the slowest mode's
root-mean-square frequency in the regressor"; md Part C ("the filter's −102.14 = −ω₁², ω₁ = 10.11"; "−109.72 = −ω₁², ω₁ = 10.47");
`exp_adaptive_floquet.py` l. 234–235; (c) l. 3084–3086).** ω₁ does not depend on η. By Parseval ω₁² = Σn²\|Ĵ_nv₁\|²/Σ\|Ĵ_nv₁\|² = 103.276,
so ω₁ = 10.1625 (asymmetric reference: 111.467 and 10.558; `check_C2.py` (2)). At η = σ₀ the filter's part of c₁/λ₁ is
−ω₁²·h₁(dλ₁⁻/dh)/λ₁⁻ = −103.276 × 0.98900 = −102.140, the printed value. So "−102.1 = −ω₁², ω₁ = 10.1" (and 10.11, 10.47) is
mislabelled at finite η.

Likewise, (c)'s "c_j contains the error filter's attenuation −v_jᵀ⟨Ĵ̇ᵀĴ̇⟩v_j" is exact only at η = ∞. At finite η its contribution to c_j^±
is −(dλ_j^±/dh)·v_jᵀ⟨Ĵ̇ᵀĴ̇⟩v_j, with dλ^±/dh ∈ (0, 1); `second_order` computes this correctly.

*Fix:* print ω₁ from the η-independent Rayleigh quotient; at η = σ₀ write "−102.1 = −0.989ω₁²"; in (c) write "(dλ_j^±/dh)·v_jᵀ⟨Ĵ̇ᵀĴ̇⟩v_j".

**D3o-7 (cosmetic; (c) l. 3082–3084 "the odd orders vanish — the first because the first-order averaging correction is skew, all of them
because the variational equation is symmetric under time reversal"; proof l. 3111–3116; §0 row l. 127 "(the first-order correction is
skew, the variational equation time-reversal symmetric)"; §8 l. 6069; audit l. 395–396).** The parity argument is right as the proof
gives it. The statement and the summaries present it in three misleading ways.
- *The reversal symmetry alone is not the reason.* M(t)ᵀ = DM(t)D⁻¹ holds for any a(t), since a enters only the diagonal (δe, δe)
  entry, and odd orders then do appear.
  - With a(t) = 0.5 + 0.3 cos t, the slow rates of M(t) and M(−t) agree to 4.7·10⁻¹¹ (`check_C2.py` (4)).
  - Yet a fit of k_eμ_j(k_e) over k_e = 25–400 gives b₁/λ = +0.44 in the slowest mode and +0.40 in the fastest.
  - What makes the expansion even is the reversal symmetry *together with* the parity (−1)^k of every term of Ā_k under reversal. That
    parity follows from a dimension count (each term has an even number of J̄'s, hence d + i ≡ k mod 2), which a rate a(t) breaks.
  - The proof's last sentence says this. The statement and the summaries name only the symmetry.
- *D at η = ∞.* D = diag(−I_n, I, (σ₀/η)I) is singular at η = ∞. There the relation holds on (δe, y), with D = diag(−I_n, I).
- *Degenerate eigenvalues.* At a degenerate eigenvalue of Ā₀ the first-order splitting ±iεκ is imaginary, so the complex exponents
  have odd terms and only the rates are even. There c_j's formula needs degenerate perturbation theory. The instance has no degenerate
  eigenvalues and no complex multipliers.

*Fix:*
- "the odd orders vanish because M and its time reverse have the same multipliers and every term of Ā_k has the parity (−1)^k under
  reversal, which a ≢ 0 breaks";
- D = diag(−I_n, I) at η = ∞;
- "(the rates; c_j for simple eigenvalues)".

**D3o-8 (cosmetic; numbers and wording of D‴.28).**
- *λ⁺.* *Numerics* l. 3303 "λ⁺_k ≤ 7.2·10⁻⁴ (at most 0.11σ₀)": 7.199·10⁻⁴/7.0744·10⁻³ = 0.102, so "0.10σ₀". The §7 row's
  "λ⁺ ≤ 0.04σ₀ at ε ≤ 0.3" is 0.0404.
- *The Gauss–Newton gap.* The §0 and §7 rows say "the Gauss–Newton form off by 7·10⁻⁴–3·10⁻³". The paragraph and RESULTS say
  7.0·10⁻⁴–5.6·10⁻³, which includes ε = 0.5.
- *"Constant from s = 600."* *Numerics* l. 3306 "stays at its value from s = 600 on" (§7 row and RESULTS: "constant from s = 600"). Over
  [600, 6000], ‖f̂_m − f̂_M‖² changes by up to 8.3% at m = 15 and 7.2% at m = 60 (ε = 0.1; 6.1% and 6.3% at ε = 0.3), and by ≤ 0.08% at
  m = 240; h2(6000)/h2(600) = 0.917–1.003.
- *The margins table.* It compares the curvature computed by central differences at G_wϑ^b(q) with (q_k/w_k)λ⁺_k — the identity of
  (iii) — not with the dynamics' separation eigenvalues. Those agree too (confirmation 5), which the text could say.
- *A docstring.* `exp_gauge_floor.py` l. 116–117 writes "g̃ = D^{1/2}∂_w(G_wϑ^b)". The code (l. 132) and the text use
  g̃ = −D^{1/2}∂_w(G_wϑ^b).
- *An article.* §8 l. 6072 "the gauge transform is the critical point at the sampled weights" should read "a critical point": L_w has
  others, which D3o-1's plain Newton finds.

**D3o-9 (cosmetic; the tests).** `tests/test_gauge_floor.py` passes 4/4; its asserted quantities (`check_tests_A.py`):

| test: asserted | measured | tolerance |
|---|---|---|
| 1: ∇L_w(G_wϑ^b) | 3.4·10⁻¹⁶, 5.1·10⁻¹⁶ | 10⁻¹² |
| 1: f̂ unchanged | 1.1·10⁻¹⁶, 2.2·10⁻¹⁶ | 10⁻¹⁴ — an identity for every Θ |
| 2: ‖F_a‖/‖F_f‖ | 2.3·10⁻¹¹ (Newton moved 0.0) | 10⁻⁸ — tautological (D3o-1) |
| 2: the full-Hessian formula | 3.1·10⁻⁹ | 10⁻⁶ |
| 2: not the Gauss–Newton form | 3.00·10⁻³ | > 10⁻³ (a factor 3) |
| 3: amplitude row | 7.5·10⁻¹⁴ | 1.27·10⁻¹³ — a factor 1.7, the central-difference error of `residual_curvature` |
| 3: scaling | 1.5·10⁻¹¹ | 10⁻⁹ (relative) |
| 4: \|means − G_wϑ^b\| | 1.7·10⁻⁶ | 10⁻⁴ |
| 4: \|means − ϑ^b(q)\| | 1.228 | > 1.0 |
| 4: f̂ | 6.6·10⁻⁸ | 10⁻⁶ |

- A fixed-prior mutant fails test 4 (0.46, 0.77 and 5·10⁻²), which is good.
- Test 4 runs at η = ∞ only. At η = σ₀ the same run ends 3.6·10⁻⁴ from the gauge transform at s = 4·10⁴, above the tolerance.
- Test 3 checks (iii)'s identity in the code's own curvature, never the dynamics' margins.

`tests/test_adaptive_floquet.py`, the three new tests, all pass (`check_tests_B.py`):
- *Test 2:* separation residual 1.6·10⁻¹⁴ (tolerance 10⁻⁹); the generator to 4.8·10⁻¹¹ (10⁻⁸).
- *Test 3:* the values in its comments (0.133 and 0.0263).
- *Test 4:* ‖Ā₁‖ = 8.1·10⁻³, skew to 9.7·10⁻¹⁷; the residual ratio 15.55; res(200) = 1.94·10⁻⁵; meas/c₁ − 1 = −0.0071; Liouville
  1.2·10⁻⁵ and 1.5·10⁻⁶ against a tolerance of 2·10⁻⁴.
- *Mutants of the averaging.* A₂ = ⟨B₂⟩ fails (ratio 9.0, −0.0164), and removing the filter term fails (ratio 3.9), so the tests are not
  vacuous. Three weaknesses remain:
  - *An Ā₁ = 0 mutant passes test 4 unchanged* (ratio 15.44, −0.0072). Ā₁'s contribution to c_j is ≤ 2.8·10⁻⁴ of c_j on this
    reference (`check_A1_part.py`). The second-order checks therefore never involve the term the asymmetric reference was chosen for;
    only the direct skewness assertion does.
  - *The metric is never tested.* Skewness is checked at η = σ₀, where the metric diag(I, (σ₀/η)I) is I. At η = 10σ₀, Ā₁ is skew in
    every metric diag(I, cI), because it lives on the y block. Ā₀, in contrast, is self-adjoint only in the stated metric (8.7·10⁻¹⁸,
    against 6.4·10⁻² in I), and no test checks Ā₀.
  - *Test 3's `var < p1` asserts an instance property* (D3o-5).

*Fix:*
- Start test 2's Newton at `tb`.
- Assert the dynamics' separation eigenvalues against σ₀ − (q_k/w_k)λ⁺_k (a 15-particle Jacobian, under 1 s).
- Add η = σ₀ to test 4, with a longer horizon.
- In test 4 of the Floquet tests, check Ā₀'s self-adjointness at η ≠ σ₀.
- Label `var < p1` as an instance property.

**D3o-10 (cosmetic; novelty: D‴.28 Status l. 3285–3289, its §0 row l. 129 (last column), audit l. 410–416; D‴.26 §0 row l. 127 (last
column), audit l. 387–399).**
- *D‴.28.* The facts used are the rescaling-symmetry literature:
  - a symmetry of the loss maps critical points to critical points with the same function;
  - a rescaling changes the curvature at a minimum without changing the function: Dinh–Pascanu–Bengio–Bengio, *Sharp minima can
    generalize for deep nets*, ICML 2017 (rescalings of rectifier networks give equivalent models at arbitrarily sharper minima);
    Neyshabur–Salakhutdinov–Srebro, *Path-SGD*, NeurIPS 2015;
  - for permutation symmetries, the symmetry-induced critical points of Şimşek et al., ICML 2021.

  R_k → (q_k/w_k)R_k is the scalar-amplitude instance of that reparametrisation dependence, and the margin formula is its corollary
  (D3o-2). The audit's "a new identification of the mechanism" is fair. "An exact margin formula" should say "for a scalar amplitude".
- *D‴.26.* The evenness is the standard time-symmetry argument, which the audit presents as new:
  - an expansion that is even in a small parameter because the problem is invariant under a symmetry reversing that parameter's sign
    is the mechanism by which symmetric integrators have even error expansions (Hairer–Lubich–Wanner, *Geometric Numerical
    Integration*, 2nd ed., 2006, the chapter on symmetric integration and reversibility);
  - it is also why the time-symmetric Magnus exponent, Ω(t_f, t_0) = −Ω(t_0, t_f), has only odd terms about the midpoint
    (Blanes–Casas–Oteo–Ros 2009, the time-symmetry section — the reference the audit already cites for the commutator);
  - "M(t) and M(−t) have the same multipliers" for a symmetric DM follows from the adjoint system's monodromy P^{−T}, standard Floquet
    theory.

  The audit's "classical tools, a new structural statement (… and its even expansion)" therefore overstates the evenness. New are the
  zero-residual lift at every representable sample with its m-independent spectrum, and the explicit c_j of this loop.

*Fix:* cite these, and write "the even expansion, a consequence of the reversibility of the adjoint pair (as for symmetric integrators
and the Magnus expansion's time symmetry)".

Sources: Dinh et al. 2017 and Şimşek et al. 2021 were checked through a web search (bibliographic details and abstracts). Blanes et al.
2009's time-symmetry section was seen only through a search snippet. Hairer–Lubich–Wanner's even-expansion theorem was checked only
through a secondary source (arXiv:1507.01202, Thm 14); the book's numbering was not checked. Everything else is from memory.

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and tests.** `--from-raw` rewrites both md files byte-identically (md5 `446c9252…`, `d689083f…`). Both test files pass:
   8 tests in 123.85 s.
2. **Part A of `gauge_floor.md`, independently** (`check_A.py`). The independent code's analytic gradient and Hessian agree with
   finite differences to 7.0·10⁻¹¹ and 1.8·10⁻¹⁰. Every number of the two Part A tables reproduces to all printed digits:
   - the fit errors 3.5914·10⁻³ … 2.5275·10⁻² and the smallest eigenvalues of Ã, 3.0587·10⁻⁴ … 4.6022·10⁻⁵;
   - ‖R̃‖/‖H_L‖ = 5.290·10⁻⁴ … 2.826·10⁻³. At ε = 0 it is 6.4·10⁻¹⁶ analytically, against the md's finite-difference 1.95·10⁻¹¹;
   - W_g's thresholds, λ⁺_k (to 10⁻⁹ relative) and max\|R\| = 2.642·10⁻⁴ … 2.909·10⁻³;
   - the amplitude row: analytically ≤ 4.2·10⁻¹⁶;
   - the fixed prior's F against the full-Hessian formula: 1.9·10⁻⁹–2.0·10⁻⁹, with the analytic gauge derivative u;
   - against the Gauss–Newton form: 6.99·10⁻⁴, 1.32·10⁻³, 2.30·10⁻³, 3.00·10⁻³ and 5.64·10⁻³; the floor laws 6.16616·10⁻² … 8.52303·10⁻²
     are identical.

   The counts of the margins table (15, 231, 210, 6, 105 count vectors) also reproduce (`check_Wg_prob.py`).
3. **ϑ^b(q) is the best weighted approximant** (`check_branch.py`). Of 200 L-BFGS multi-starts per ε = 0.1, 0.3, 0.5, none ends below
   the branch, and 48, 59 and 48 reach it, all with positive-definite Hessians. Ã ≻ 0 along the branch up to ε = 0.65 (smallest
   eigenvalue 8.5·10⁻⁷ there).
4. **(ii) re-derived** (see the verdict). The derivation uses A = JᵀJ − blockdiag(q_kR_k), so Ã = H_L − blockdiag(G^{1/2}R_kG^{1/2}).
   At an exact representation g̃ = ∂_wx⁰ and Ã = H_L, which is D‴.24(i).
5. **(iii) at the level of the dynamics** (`check_margins_dyn.py` (1)). Independent averaged nearest-atom drift with the prior in the
   state, m = 15 particles at G_wϑ^b(q):
   - *ε = 0.3, counts (6, 4, 5).* The 36 separation eigenvalues equal −σ₀ + (q_k/w_k)eig(G^{1/2}R_k(q)G^{1/2}), each n_k − 1 times, to
     5.4·10⁻¹² at η = σ₀ and at η = ∞. The separation subspace is invariant to 1.7·10⁻⁹. The collective eigenvalues are all negative
     (largest −1.27·10⁻⁴ at η = σ₀, −2.55·10⁻⁴ at η = ∞).
   - *ε = 0.5, counts (5, 5, 5).* The same agreement, to 5.0·10⁻¹² (largest collective −2.3·10⁻⁵ and −4.6·10⁻⁵).
   - The committed `AdaptivePriorLoop.drift` gives the largest separation eigenvalue −6.716803·10⁻³, equal to minus the md's κ_R margin.
6. **The fixed prior's splits are its D‴.5 saddles** (`check_margins_dyn.py` (2); its stationary points solved independently). At
   ε = 0.1:
   - m = 15: seeds 0, 1 and 3 (margins −3.6·10⁻³, −1.3·10⁻³, −0.35) split, with diameters 1.1, 0.67 and 1.5; seed 2 (+4.4·10⁻³) does not.
   - m = 60: seeds 1 and 3 (−2.3·10⁻³, −4.1·10⁻³) split.
   - m = 240 and M = 510: positive margins in every seed.

   ε = 0.3 is the same.
7. **Part B from the raw records** (`recompute_B.py`). Every number of the md table, the near-boundary table, the excess sentence, the
   paragraph, the §7 row and RESULTS reproduces:
   - the medians 1.21·10⁻⁴, 1.21·10⁻⁴ and 1.63·10⁻⁵;
   - the maxima 0.113 (ε = 0.1, seed 1, margin 0.0218), 1.20·10⁻² and 1.29·10⁻³;
   - the counts conserved in all 24 jobs and 4 sizes.

   Five m = 15 systems were re-run with the committed code and reproduce the stored fit records bitwise (max difference 0.0): ε = 0.1
   seed 3 at η = σ₀ and ∞; ε = 0.3 seed 3 at η = ∞; ε = 0.3 seed 1 at η = σ₀ and ∞.
8. **The near-boundary sample converges to its gauge transform.** See D3o-4's bullet: excess 0.038 and 0.019 at s = 6·10⁴, and a box
   minimum equal to the approximant's. This confirms the paragraph's "the approach is slow".
9. **Part C's monodromies** (`check_C1.py`, `check_C4.py`).
   - The independent exponents equal the raw ones to ≤ 2.0·10⁻⁹ (k_e = 400) and ≤ 1.6·10⁻¹⁰ (k_e ≤ 200), with Liouville's error channel
     to ≤ 3·10⁻¹⁴.
   - DOP853 at rtol 10⁻¹¹, 10⁻¹² and 10⁻¹³ agrees on the slowest rate to 1.9·10⁻⁹ (k_e = 400).
   - Over six integrations (five start times and a product of two half periods) the slowest rate spreads by 1.0·10⁻¹⁵ absolute
     (2.6·10⁻⁹ relative) at k_e = 400, so the second-order residual (4.88–4.90·10⁻⁷) and the ratio 15.97 are resolved.
   - Hill's method (N_h = 40, 56) agrees to ≤ 8·10⁻⁹ at k_e = 25 and 100, and gives the error channel directly: 24.98083504 and
     99.99509993, Liouville's values.
   - A fourth-order Magnus integrator converges like h⁴ up to N = 8000 (successive differences of the slowest rate 3.2·10⁻¹², 1.8·10⁻¹³)
     and is rounding-limited beyond, at ~2·10⁻¹³ absolute (`check_C3.py`).
10. **The averaged generators** (`check_C2.py` (1)). The own FFT implementation of the slaving and the second-order averaging reproduces
    λ_j and c_j to ≤ 1.9·10⁻¹⁴, and N = 256 is already converged (≤ 1.7·10⁻¹³ to N = 2048). Every c_j is negative on the instance:
    c₁/λ₁ = −100.55, −101.67, −108.21 and −109.93 (sym/asym × η = σ₀/∞). The skewness: the symmetric parts of ⟨B₁⟩ and ⟨B̃₀U₁⟩ are
    ½⟨(J̄ᵀJ̄)˙⟩ and ½⟨(U₁²)′⟩ (U₁ is symmetric), both on the y block.
11. **The evenness, measured.** (k_eμ₁ − λ₁)k_e²/λ₁ converges to c₁/λ₁ with ratios per doubling of 2.65, 3.47, 3.84, 3.96, 3.99 (sym)
    and 1.99, 2.89, 3.59, 3.88, 3.97 (asym), i.e. like k_e⁻², so there is no k_e⁻³ term (`recompute_C_md.py`). A fit with odd terms on the
    asymmetric reference gives b₁/λ ≤ 2.6·10⁻³ (`check_C2.py` (3)). Ā₁'s first-order invisibility is a real test: a symmetric matrix of
    its size would move the slowest rate by up to 26λ₁/k_e, while the diagonal of Ā₁ in Ā₀'s eigenbasis is ≤ 9·10⁻¹⁷
    (`check_A1_part.py`).
12. **The generator against the committed real loop's central-difference linearisation**: 9.0·10⁻¹¹ (asymmetric reference, η = σ₀) at
    k_e = 25 and at 400. Ĵ(t) equals `whitened_jacobian`'s to 1.1·10⁻¹⁶.
13. **Every number of Part C, the D‴.26 *Numerics*, the rows and RESULTS** (`recompute_C_md.py`):
    - the slowest deviations 37.53, 13.56, 3.84, 0.994, 0.251 and 0.0628%;
    - the ratios 2.77, 3.53, 3.86, 3.96, 3.99 and 11.05, 14.00, 15.41, 15.84, 15.97;
    - the c_j/λ_j range −100.55 to −0.47 and the measured −99.37, −100.25, −100.48;
    - asym: −108.21, of which the filter's −109.72, and the ratios 15.54, 15.87;
    - (err − (k_e − tr H/k_e))k_e³ = 6.22 … 7.79 (sym) and 7.01 … 10.73 (asym): O(k_e⁻³).
14. **The tests' quantities** (D3o-9's tables).

## Not checked

- The experiments were not rerun in full: Part B's 24 jobs (five m = 15 systems reproduced bitwise), the margins job over every count
  vector (its identity recomputed; the dynamics' margins at one count vector per ε), and the first two tables of `exp_adaptive_floquet.py`
  (the finite-difference monodromies, checked in round 19; Part C agrees with them to 4.1·10⁻⁵).
- With vector amplitudes only the stationary curvature and the margins' formula were tested, not the dynamics.
- D‴.29's use of D‴.28 (the servo, a scalar amplitude, so D3o-2 does not apply there) was not reviewed.
- Not attempted:
  - the nonlinear convergence of the adaptive dynamics to the gauge transform (numerical by the authors' account), and the basin of
    perfect tracking uniformly in m;
  - a bound on the O(k_e⁻⁴) remainder, and c_j at degenerate eigenvalues.
- The concurrent session's untracked files were not reviewed.
- The literature statements of D3o-10 are from memory, except as noted there.
