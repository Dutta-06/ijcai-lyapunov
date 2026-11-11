# Referee report D3m (round 19) — PARADIGM §5d: Corollary D‴.26 (statement (a)–(c), the stability consequence, proof, Status, *Numerics*), the §0 and §7 rows "| D‴.26 |", the §8 bullet "Done (Corollary D‴.26)", `results/RESULTS_PARADIGM.md` (Corollary D‴.26), the LITERATURE audit bullet, `src/core/adaptive_floquet.py`, `tests/test_adaptive_floquet.py`, `experiments/exp_adaptive_floquet.py` → `results/chaos/adaptive_floquet.md` (+ raw); and the planar check of Proposition D‴.24 (the paragraph "*The planar check*", the clause of the D‴.24 §7 row, the RESULTS bullet "In two dimensions"), `src/core/planar_loop.py`, `tests/test_planar_loop.py`, `experiments/exp_adaptive_floor_2d.py` (Parts A–C and the checks stage) → `results/chaos/adaptive_floor_2d.md` (+ raw)

Branch `claude/festive-bardeen-eopej2` at 442cd6c, HEAD throughout the review. The concurrent session added only untracked files
(`experiments/exp_gauge_floor.py`, `results/chaos/gauge_floor_raw.npz`, `src/core/sliding.py`, `tests/test_gauge_floor.py`,
`tests/test_sliding.py`). No tracked file differs from 442cd6c (`git diff --stat 442cd6c` over the reviewed files is empty). Every script
imports the 442cd6c tree, extracted with `git archive` into the scratch directory. Line numbers are 442cd6c's:
- PARADIGM: D‴.26 l. 3047–3097 (the stability consequence l. 3061–3064, the proof l. 3065–3078, the Status l. 3079–3087, the *Numerics*
  l. 3089–3097); the planar check l. 3006–3019; D‴.24(vi)'s gloss l. 2899–2900; the §0 row l. 124; the §7 rows l. 5587 (D‴.24) and
  l. 5588 (D‴.26); the §8 bullet l. 5832–5835.
- RESULTS: l. 593–599 (2D bullet), l. 601–608 (tests), l. 610–621 (D‴.26).
- LITERATURE_TRANSPORT: l. 385–391.

No repository file was modified, and no git command touching the index or the tree was run. Python ran with `PYTHONDONTWRITEBYTECODE=1`
and one thread; pytest ran with `-p no:cacheprovider` on the scratch tree. `--from-raw` on the scratch tree rewrites both md files
byte-identically: `adaptive_floquet.md` md5 `e43e0633…`, `adaptive_floor_2d.md` md5 `664a380a…`, both equal to `git show 442cd6c:…`.

Scratch scripts and their `.out` files are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3m/scripts/`:

| script | what it does |
|---|---|
| `inst1d.py` (module) | independent re-implementation of the three-bump instance: feature, analytic parameter gradient, ϑ^e(w), H_L(w), pair rates, and the variational equation of (CL) along the lift (no repository code) |
| `check_instance.py` | f̂ ≡ f and f̂′ ≡ f′ at ϑ^e(w); ∇Φ̃ against finite differences; h_j(w) against the raw record |
| `recompute_floquet_md.py` | every number of `adaptive_floquet.md`, the paragraph, the rows and RESULTS from the raw records, independently of `report()`; per-mode deviations and their k_e-scaling |
| `var_monodromy.py` | the monodromy by the variational equation. K₀-atom loop at k_e = 12.5–400 (η = σ₀, ∞); lifts m = 15, 30 at k_e = 25, 50 against the raw finite-difference records; Liouville's formula for the error channel |
| `fd_resolution.py` | the experiment's monodromy (committed `AdaptiveRealLoop` + `monodromy`) at steps 10⁻⁴, 10⁻⁵, 10⁻⁶ and rtol 10⁻¹², 10⁻¹³ against the variational one (K₀-atom loop, η = σ₀, k_e = 50) |
| `var_eta_mech.py` | (A) η = 0.1, 1, 10; (B) the mechanism of the O(k_e⁻²) remainder (filtered Gram, the skew corrections on the instance); (C) a term a(t) in the linearised error equation |
| `var_asym_ref.py` | (c) with the reference 2 sin t + 0.8 sin 2t, whose first-order correction is non-zero |
| `var_boundary.py` | (c) at weights approaching ∂W |
| `check_proof_M.py` | the proof's M(t) against the linearisation |
| `check_saddle_1d.py` | the fixed prior's D‴.5 margins at (6, 5, 4)/15 and at ϑ^e(w) |
| `check_ot_plan.py` | the OT plan at the lift of ϑ^e(q) under perturbations (the OT extension) |
| `check_tests.py` | the asserted quantities of both test files; a fixed-prior mutant; the branches the tests never execute |
| `inst2d.py` (module) | independent re-implementation of the planar instance (force, Θ*(w) by Levenberg–Marquardt, margins, collective block) |
| `recompute_2d_md.py` | every number of `adaptive_floor_2d.md`, the paragraph, the §7 clause and RESULTS from the raw records |
| `check_planar.py` | `planar_loop` against the independent implementation, the real law, the 1-D nearest-atom loop and the replicated system; the checks stage for all 24 weight vectors; R_k by second differences |
| `check_partA_lag.py` | Part A's ‖F − F_bias‖/‖F‖ at three steps; the η = ∞ refresh lag of `planar_loop.run` |
| `check_finite_eta.py` | the finite-η branch of `planar_loop.run` against DOP853; its order |
| `check_partB_outside.py` | three non-representable Part B systems, re-run bitwise, continued to s = 2.4·10⁴–6·10⁴, against the box-constrained minimum |
| `check_partC.py` | Part C's exact start: bitwise to t = 160; the same at dt = 10⁻³ and 4·10⁻³ |
| `derived.py` | ratios computed from the printed outputs of the scripts above |

Also in `ref_D3m/`: `tests.out` (5 passed, 67.36 s), `from_raw_floquet.out`, `from_raw_2d.out`, `md5_before.out` and `md5_after.out`. The
experiments were not rerun in full. Single systems were reproduced bitwise (confirmations 3 and 12).

## Verdicts

| object | verdict |
|---|---|
| **D‴.26 (a)** | **holds, exactly**, for finite η and for η = ∞. Re-derived: at e ≡ 0 both G∇²_ϑ(eᵀΦ̃) and the x-derivative of the learning law vanish. A separation leaves Df̂[δΘ] = 0, hence δe ≡ 0. Since δϑ̄_j = 0 the prior's equation is untouched, so δz ≡ 0. By the variational equation the 36 and 81 separation exponents equal σ₀/k_e to 1.2·10⁻¹²–3.1·10⁻¹² relative. The md's 1.7·10⁻⁶ is the finite-difference resolution. |
| **(b)** | **holds, exactly.** The collapsed subspace is invariant with the prior in the state, and the loop on it is the weighted K₀-atom loop with the unweighted learning law. The lifts' collective exponents equal the K₀-atom loop's to ≤ 1.8·10⁻¹⁰ (variational), for both m and both η. "The same for every m" means every m with the same sampled weights. |
| **(c)** | **holds.** In slow time the remainder is O(k_e⁻²): the slowest mode's successive ratios for k_e = 12.5 … 400 are 2.77, 3.53, 3.86, 3.97, 3.99, and a fit in k_e⁻¹, k_e⁻², k_e⁻³ gives coefficients 0.078λ₁ and −115λ₁ for the first two (η = ∞). The skewness mechanism is confirmed where the first-order correction is non-zero: with an asymmetric reference the commutator term is 5.6·10⁻³ and skew, and the deviations still fall like k_e⁻². Three problems. On the instance the first-order correction vanishes identically, so the numerics cannot test the mechanism (D3m-1). The linearised error equation uses ∂_xf̂ = ∂_xf on the orbit, which D‴.24(vi)'s gloss does not supply (D3m-2). The proof's M(t), N(t) lack a factor W⁻¹ (D3m-4). |
| stability consequence | **holds.** For each representable w and every m (with these weights), all exponents are positive for k_e above a threshold. Local exponential stability for each m follows (the period map has a hyperbolic attracting fixed point, and the targets are locally constant). The rate formula holds, with a large negative constant in its O(k_e⁻²): −13.6 % at k_e = 25 (D3m-5). Near ∂W the relative correction stays bounded, while the rate itself tends to 0. Consistent with the Status's "basin uniform in m not proved"; the §0 row drops "uniform in m" (D3m-5). |
| the OT extension | **holds.** At the lift of ϑ^e(q) the optimal plan is the class assignment. It is unchanged under 200 of 200 random perturbations up to G⁻¹-radius 2 per particle and per atom (inter-atom distances 6.3–7.2). The plan is locally constant, so the linearisation is the nearest-atom one. |
| the fixed-prior claim | **holds.** With a fixed prior the lift solves (CL) iff every class's prior atom is ϑ^e_k(w). At (6, 5, 4)/15 the fixed matching prior's Θ*(w) has margins (−3.48·10⁻³, +6.20·10⁻³, +7.17·10⁻³), while ϑ^e(w) has σ₀ in every class. (a) holds where the fixed prior's orbit is a saddle, at a different configuration. |
| Proof | **correct in structure** (slaving, then second-order averaging of a self-adjoint periodic block). It has slips (D3m-4) and the unstated derivative condition (D3m-2). |
| *Status* | **accurate** on what is proved. The novelty is overstated (D3m-3). The Checked list omits that the test runs at η = ∞ and w = q only (D3m-6). |
| *Numerics* (D‴.26) | **every number reproduces** (confirmation 2). The reading "second order in 1/k_e, as (c) says" cannot fail on this instance (D3m-1). The "factor 3.80" is the fastest mode's; the sign of the remainder is unstated (D3m-5). |
| `src/core/adaptive_floquet.py` | **holds.** The variational equation of an independent implementation reproduces the experiment's exponents to the finite-difference resolution. The docstring's pointer to `lifted_orbit.monodromy_fd` is stale (D3m-5). |
| `experiments/exp_adaptive_floquet.py` | **holds**; regeneration is byte-identical. Measured resolution on the K₀-atom loop (η = σ₀, k_e = 50): the slowest exponent to 4.1·10⁻⁵ relative at the chosen step and tolerance; 1.0·10⁻⁵ at step 10⁻⁴, 6.6·10⁻⁴ at 10⁻⁶, 6.8·10⁻⁷ at 10⁻⁴ with rtol 10⁻¹². The docstring calls (6, 5, 4)/15 "admissible" (D3m-7). |
| `results/chaos/adaptive_floquet.md` | **numbers hold**; wording as in D3m-1. |
| `tests/test_adaptive_floquet.py` | **passes; not vacuous, but weaker than its name** (D3m-6). It runs at η = ∞ and w = q only. A fixed prior at the exact atoms passes three of its four assertions. The last tolerance is 0.2 against a measured 0.137. |
| §0, §7 rows, §8 bullet, RESULTS, audit bullet (D‴.26) | **numbers hold.** Wording: D3m-1 (§0, §7, RESULTS), D3m-3 (§0 last column, audit bullet), D3m-5 (§8 O(k_e⁻²), RESULTS 14 %/3.9 %, §0 "the nonlinear basin not proved"). |
| *The planar check* paragraph | **every number reproduces** (confirmation 10). The η = ∞ of Part B is a lagged prior, and Part C's adaptive prior is η = ∞ only, unsaid (D3m-8). "Still falling toward a constrained minimum" is **justified** (confirmation 13). |
| `src/core/planar_loop.py` | **holds.** The drift agrees with an independent implementation to 3.3·10⁻¹⁶ relative; the real planar law, the 1-D nearest-atom loop and the replicated system agree to ≤ 3·10⁻¹⁷. The finite-η branch is second order (1.6·10⁻⁷ against DOP853 at dt = 0.25). η = ∞ resets the atoms every `cell_every` steps (D3m-8). |
| Part A | **holds.** ‖F − F_bias‖/‖F‖ = 6.0·10⁻⁵, 6.0·10⁻⁷, 6.1·10⁻⁹ at δ_w = 10⁻³, 10⁻⁴, 10⁻⁵: the finite-difference error, as stated. The law, bound and shares hold. |
| Part B | **holds.** The margin's definition with vector amplitudes is right (max over components). "Falling" means between the last two records (7/7 and 11/11). Outside W the systems are not stalled: they stay replicated, their amplitude coordinates run off like ln s, and L̄ decreases toward the box-constrained minimum at the sampled weights (2.1 %, 14 % and 57 % above it at s = 6000 in the three systems continued). |
| checks stage | **holds.** Θ*(w) agrees with an independent Levenberg–Marquardt solution to ≤ 1.2·10⁻¹⁴ for all 24 weight vectors, and the margins and collective eigenvalues agree to all printed digits. R_k agrees with the Hessian of the scalar potential ⟨rᵀΦ̃(ϑ)⟩ by second differences (one class checked). Every sample is D‴.5-admissible. Observed/exact bias is 0.99999979–1.00000007 over the 18 systems. |
| Part C | **holds.** The windows, the averaged prediction and the ratios reproduce. The exact start is at the integrator's level: the five-period rms scales by 16.2–16.4 per halving of dt (RK4), and the stored e² series reproduces bitwise. |
| `tests/test_planar_loop.py` | **passes, not vacuous** (4/4). `run()`'s finite-η branch is never executed, and test 3 calls representable weights "admissible" (D3m-6, D3m-7). |
| D‴.24 §7 clause, RESULTS "In two dimensions" | **numbers hold.** Wording: "four kernels", "the adaptive prior's floor falls", "admissible weights" (D3m-7, D3m-8). |

## Confirmed issues

**D3m-1 (minor; *Numerics* l. 3094–3097 "a factor 3.80, i.e. k_e^{−1.93}: second order in 1/k_e, as (c) says — largest in relative terms in
the slowest mode"; the md's sentence "the averaging correction is of second order in 1/k_e in slow time, so the exponents are (1/k_e)λ_j^±(w)
+ O(k_e⁻³) (Corollary D‴.26(c))"; the §0 and §7 rows "(a factor 3.80, k_e^{−1.93})"; RESULTS l. 617–618).** The observed second order
cannot test (c)'s mechanism on this instance. Both first-order (in slow time) corrections vanish identically, for reasons that have nothing
to do with self-adjointness (`var_eta_mech.py` (B)):
- ⟨N⟩ = ⟨ĴᵀĴ̇⟩ = (1/T₀)∮ψ(x)dx = 0 for every one-dimensional reference (the integrand is an exact differential); its norm is 3.2·10⁻¹⁶.
- The commutator term ⟨F̃U⟩ vanishes by the symmetry t → π − t of x_m = 2 sin t; its norm is 2.8·10⁻¹⁸.

So the k_e⁻² scaling would be observed even if the self-adjointness argument were wrong. Two further gaps:
- The scaling rests on two values of k_e. The quoted factor 3.80 is the fastest mode's; per mode the factors are 3.53–3.95, and 3.54
  (k_e^{−1.82}) in the slowest mode, whose relative deviation is the one quoted.
- The text does not say what the O(k_e⁻²) remainder is.

The claim itself is right (`var_asym_ref.py`). With x_m = 2 sin t + 0.8 sin 2t:
- ⟨F̃U⟩ has norm 5.6·10⁻³ and is skew (symmetric part 7·10⁻¹⁸).
- Its first-order diagonal coefficients v_jᵀKv_j/h_j are ≤ 4.3·10⁻¹⁵.
- The exponents' deviations still fall like k_e⁻²: the slowest mode at −12.5 %, −3.96 %, −1.07 %, −0.273 %, −0.0686 % for k_e = 25 … 400
  (ratios 3.15, 3.70, 3.92, 3.98).
- A symmetric matrix of the same size would shift the slowest rate at first order (−5.9·10⁻³, −2.9·10⁻³, −1.5·10⁻³ at k_e = 25, 50, 100).

On the instance the remainder is the error filter's attenuation of the regressor's harmonics, which is not skew. The effective Gram is
Σ_nĴ_n^HĴ_n·k_e²/(k_e² + n²) = H_L − k_e⁻²⟨Ĵ̇ᵀĴ̇⟩ + …. Its relative coefficient is −v_jᵀ⟨Ĵ̇ᵀĴ̇⟩v_j/h_j = −103 in the slowest mode
(fitted c₂/λ₁ = −107 at η = ∞, −106 at η = σ₀). The exact filtered Gram reproduces the observed deviations at k_e = 25 to within 0.8–6.9 %
mode by mode (1.4 % in the slowest). This is why every rate, in every mode and at every k_e from 12.5 to 400, lies below the prediction.

*Fix:*
- Replace "as (c) says" by: "on this reference the first-order correction vanishes identically (⟨N⟩ = 0 in one dimension, the commutator
  term by the symmetry of 2 sin t); with x_m = 2 sin t + 0.8 sin 2t, where it does not vanish, the deviation still falls like k_e⁻²".
- Add k_e = 100 (and 12.5) to the experiment, or a non-symmetric reference.
- Quote the per-mode factors 3.53–3.95.
- Name the remainder: the error filter's attenuation −k_e⁻²⟨Ĵ̇ᵀĴ̇⟩, which always slows; −13.6 % in the slowest mode at k_e = 25.

**D3m-2 (minor; proof l. 3065 "δė = −k_eδe − Df̂[δΘ](x_m(t))"; statement (c) "the error channel (≈ k_e)"; D‴.24(vi) l. 2899–2900 "f̂ = f on
the reference orbit — the support of the occupation measure ρ, which is all the error equation sees").** The linearised error equation
is δė = (∂_xf − ∂_xf̂ − k_e)δe − Df̂[δΘ]. Proposition D‴.15's proof keeps this term. The corollary drops it, which needs ∂_xf̂ = ∂_xf
along the orbit.

That holds because the exact representation is an identity of functions on ℝⁿ (D‴.24's hypothesis). On the 1-D instance
max|f̂′ − f′| = 5.6·10⁻¹² on [−4, 4] (`check_instance.py`); on the planar instance the representation is the same sum of bumps by
construction. It does not follow from "f̂ = f on the reference orbit", which is what (vi)
says the error equation sees. For n ≥ 2, where supp ρ is a curve, f̂ = f on supp ρ allows a transverse derivative A(t) ≠ 0. Then:
- (a) and (b) survive, since separations have δe = 0.
- The slaved error becomes −(k_e − A)⁻¹Ĵy, which adds the non-skew term −(1/k_e)ĴᵀĀĴ to the slow-time block. The first-order correction
  no longer vanishes, and the error channel moves by O(1).

Demonstration with a synthetic constant A = 0.5 added to the linearised error equation (η = ∞, `var_eta_mech.py` (C), `derived.py`):
the slowest rate shifts by +0.0128, +0.0090, +0.0049, +0.0025 relative to the filter-only case at k_e = 25, 50, 100, 200. These tend to
a/k_e = 0.02, 0.01, 0.005, 0.0025: a first-order shift.

*Fix:*
- In (c) and the proof, write "with ∂_xf̂_{ϑ^e(w)} = ∂_xf along the orbit (f̂ = f in a neighbourhood of it; on both instances everywhere)".
- Amend (vi)'s gloss: e ≡ 0 needs f̂ = f on supp ρ; its linear stability also needs ∂_xf̂ = ∂_xf there.

**D3m-3 (minor; §0 row l. 124 last column "new is the exact, m-independent spectrum of perfect tracking with a leak, kept for its
finite-size and transport properties, and the vanishing of the first averaging correction by self-adjointness"; audit bullet l. 385–391
"What is specific here … a leak … coexists with exact tracking at every representable sample … Verdict: classical tools, a new structural
statement").** Most of the structure is inherited or classical.
- *(a) and (b) are the document's own earlier results at e ≡ 0.* Corollary D‴.14(ii) already states "at the matching masses, where the
  residual vanishes, the ghosts' exponent is exactly σ₀/k_e". Proposition D‴.15 gives the block structure and the m-independence for any
  lifted orbit. D‴.26 adds the prior's coordinates, which separations do not move, and the exact representation at every representable w.
- *A leak compatible with asymptotic tracking is classical.*
  - The e-modification (Narendra–Annaswamy, IEEE TAC 32(2):134–145, 1987; covered by [4], and listed in the document's own dictionary,
    PARADIGM l. 287) scales the leak by ‖e‖, so it vanishes at perfect tracking.
  - A leak toward a centre that follows the parameters is the proximal-point/EASGD structure of [S49], whose fixed points carry no bias.
- *"The vanishing of the first averaging correction by self-adjointness" is standard.* The second Magnus term is a commutator (Magnus,
  CPAM 7 (1954); Blanes–Casas–Oteo–Ros, Phys. Rep. 470 (2009)). It is skew for a self-adjoint generator, and a skew perturbation leaves
  the eigenvalues fixed at first order. That is textbook perturbation theory.
- *Rate via the PE level.* The averaged rate of an adaptive loop is set by the PE level of the averaged regressor Gram. This is the
  classical averaging analysis of adaptive systems (Kosut–Anderson–Mareels, IEEE TAC 32(1):26–34, 1987; Sastry–Bodson [S9]). Here that
  level is h_j(w).
- *Mis-attribution.* The audit bullet attributes σ-modification's robustness/tracking trade-off to [S49]. It is Ioannou–Kokotović's ([3],
  [4]); [S49] is the proximal-point method and its relatives.

*Fix:* cite the e-modification, the proximal-point fixed-point property, the Magnus commutator and the averaging/PE analysis as
classical, and cite D‴.14(ii) and D‴.15 for (a) and (b). Restrict "new" to the combination: with the adaptive prior the zero-residual lift
exists at every representable sampled w, not only at w = q under the OT leak, and its Floquet spectrum is exactly m-independent with
collective rates λ_j^±(w)/k_e. Replace "[S49]'s trade-off" by [3], [4]. The literature statements here are from memory; no full text was
fetched.

**D3m-4 (cosmetic; proof l. 3070–3073 "Φ′(t) := [w_k∇_ϑΦ̃(x_m(t); ϑ^e_k)]_k, … M(t) = −GΦ′(t)ᵀΦ′(t) − … and N(t) = GΦ′(t)ᵀΦ̇′(t)";
"(b) The collapsed subspace is invariant (Corollary D‴.14)"; "(c) Proposition D‴.15's averaging of the variational equation").**
- *The factor W⁻¹.* The class atoms' learning law has no weight, δϑ̇_k = G∇Φ̃_kᵀδe, while Φ′ carries w_k. The linearisation is therefore
  M = −GW⁻¹Φ′ᵀΦ′ (W = diag(w) ⊗ I_p), and likewise N = GW⁻¹Φ′ᵀΦ̇′. The literal formula is wrong (`check_proof_M.py`): its average has
  eigenvalues 0.275–0.389 times h_j, and its asymmetry in the metric WG⁻¹ is 3.1·10⁻² against 1.1·10⁻¹⁶ for the correct M.
- *The citations.* Corollary D‴.14(i) is the invariance under a fixed prior; with the prior in the state the invariance needs the
  (immediate) remark that the class atoms and their prior atoms follow the K₀-atom equations. Proposition D‴.15 averages the ghost equation
  at first order. The collective block here needs second-order averaging after the slaving of δe (the Chang/Lyapunov–Perron reduction of a
  linear periodic system with a fast stable block).

*Fix:* "M(t) = −GW⁻¹Φ′ᵀΦ′ − …, N(t) = GW⁻¹Φ′ᵀΦ̇′ (in the whitened weighted coordinates −ĴᵀĴ and ĴᵀĴ̇)". Cite D‴.14 "with the prior's
equation" and "second-order averaging of the slow block after the slaving of δe".

**D3m-5 (cosmetic; consistency of the summaries).**
- The §8 bullet l. 5834 says "(1/k_e)λ_j^±(w) + O(k_e⁻²)". (c) says O(k_e⁻³), i.e. slow-time rates λ_j^± + O(k_e⁻²).
- RESULTS l. 618 "(14% and 3.9%)" against PARADIGM l. 3097 "(13.6% and 3.8%)". The slowest mode is at 13.56 % and 3.84 % for η = σ₀ and at
  13.69 % and 3.88 % for η = ∞; the two texts quote different η.
- *"Linearly exponentially stable" and "the nonlinear basin not proved"* (§0 row, §8). Local exponential stability for each m follows from
  the spectrum, since the targets are locally constant near the lift. Only a basin uniform in m (in the per-particle norm) is open, as the
  Status correctly says.
- The sign and size of the remainder are unstated. Every rate is below the prediction: all 18 (η = σ₀) and 9 (η = ∞) modes, at k_e = 12.5,
  25, 50, 100, 200, 400 (`var_monodromy.py`). The stated rate (1/k_e)λ_1^−(w) therefore overestimates the decay, by 13.6 % at k_e = 25 and
  37.5 % at k_e = 12.5.
- The error channel's "≈ k_e" is not measurable by the differences, but Liouville's formula gives it from the trace: 24.9808 at k_e = 25 and
  49.9902 at k_e = 50, against k_e − tr H_L(w)/k_e = 24.9804 and 49.9902. It is fair to state it unmeasured; it can be printed.
- `adaptive_floquet.py`'s docstring says "`lifted_orbit.monodromy_fd` … compute[s] the actual ones"; the experiment uses
  `adaptive_floquet.monodromy`.

*Fix:* correct the §8 order. Quote one η, or both, consistently. Write "exponentially stable (locally, for each m)". Add "the correction
is negative in every mode (the error filter): −13.6 % in the slowest at k_e = 25". Print the Liouville value. Update the docstring.

**D3m-6 (cosmetic; `tests/test_adaptive_floquet.py`, Status l. 3085–3087, RESULTS l. 620–621; `tests/test_planar_loop.py`).**
`test_perfect_tracking_has_the_separation_exponent_and_the_K0_atom_spectrum_for_every_m` runs at η = ∞ and at w = q only.

| asserted quantity | measured | tolerance |
|---|---|---|
| return of one period | 1.5·10⁻¹⁵ | 10⁻¹² |
| separations at σ₀/k_e | 8.7·10⁻⁸ | 10⁻⁶ |
| collective vs K₀-atom loop | 1.08·10⁻⁵ | 10⁻⁴ |
| prediction (variational: 0.1367) | 0.1367 | 0.2 |

What it leaves untested, and why it is weaker than its name:
- The prior in the state, the corollary's object, is never executed: `rhs` returns before `Zdot` at η = ∞, and `lift` adds no Z block.
- At w = q the weighted f̂ coincides with the unweighted mean, so the weighted K₀-atom path is not tested either.
- A mutant with the prior frozen at the exact atoms (a fixed matching prior) passes the first three assertions (3.8·10⁻⁷, 1.2·10⁻⁶) and
  fails only the last (23.0). The adaptive prior is tested only through (c), at a tolerance 1.46 times the measured second-order
  correction; the comment calls that correction "O(1/k_e)".
- "For every m" is the single lift m = 6.

In `test_planar_loop.py` no test executes `run()`'s finite-η branch (independently checked: second order, `check_finite_eta.py`). The η = ∞
run from the lift is invariant for any integrator.

*Fix:* add η = σ₀ and the weights (6, 5, 4)/15 (the separations, the K₀-atom equality and the 18 exponents against λ_j^±/k_e). Use the
tolerance |rel| < 2·(k_e = 25 value) with the comment "O(k_e⁻²)". Name the test "… at m = 6". Add a finite-η step test to
`test_planar_loop.py`.

**D3m-7 (cosmetic; terminology "admissible" after round 18).**
- `exp_adaptive_floquet.py` l. 6 calls w = (6, 5, 4)/15 "admissible". It is representable and *not* D‴.5-admissible: a fixed-prior saddle,
  margin −3.48·10⁻³, as the PARADIGM text says.
- `tests/test_planar_loop.py` l. 3, l. 69 (the test's name) and l. 71 ("# admissible: w > (0.2, 0.15, 0.175, 0.15)") use "admissible" for
  the representability condition.
- RESULTS l. 608 says "the exact representation of admissible weights stationary for the adaptive prior".

*Fix:* "representable" in all four places, as D3k-1 required.

**D3m-8 (cosmetic; the planar check l. 3006–3019, the D‴.24 §7 clause l. 5587, RESULTS l. 593–599, md Part B, `planar_loop.run`).**
- *η = ∞ in Part B is a lagged prior.* `run` resets Z to the cell means only at the cell refresh (every 4 steps, one slow unit) and freezes it
  in between. `AdaptivePriorLoop` (1-D) instead sets the prior to the barycentres at every drift evaluation. On seed 5, m = 16 at s = 6000
  the fit error is 1.707·10⁻³ against 1.633·10⁻³ (refresh every step) and 1.608·10⁻³ (continuous): 6.1 % in the fit error and 12.7 % in
  its square (`check_partA_lag.py`, `derived.py`). The module's docstring says so; the md and the paragraph do not. The lag is equally
  present in Part C's averaged prediction (refresh every 0.2 slow units) and in its real loop (every 0.2 time units).
- *Part C's adaptive prior is η = ∞ only.* The paragraph says "the adaptive prior's follows the averaged decay (0.979–1.009)".
- *"Still falling toward a constrained minimum" is right, but these systems have a floor.* The floor is the error of the box-constrained
  minimum of the weighted loss at the sampled weights (confirmation 13). At s = 6000 the systems are 2.1 %, 14 % and 57 % above it in L̄.
  RESULTS's "the adaptive prior's floor falls to …" uses "floor" for the representable systems, which have none.
- *"Four kernels" in the §7 clause* means four bumps (atoms); the networks have 16–128 kernels.

*Fix:* write "η = ∞ (the atoms reset to the cell means every slow unit)" and "(η = ∞)" for Part C. Add "outside W the systems stay
replicated and approach the box-constrained minimum at the sampled weights algebraically (the amplitude coordinates running off like
ln s)". Write "four bumps"; "the function-level distance falls to".

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and tests.** `--from-raw` rewrites both md files byte-identically (md5 `e43e0633…`, `664a380a…`). Both test files pass:
   5 tests in 67.36 s.
2. **Every number of `adaptive_floquet.md`, the D‴.26 paragraph, the rows and RESULTS**, recomputed from the raw records
   (`recompute_floquet_md.py`). Return residuals ≤ 1.74·10⁻¹⁵ (K₀-atom loops) and ≤ 1.69·10⁻¹⁵ (lifts). Separations: 36 and 81, at
   ≤ 1.66·10⁻⁶. Collective against the K₀-atom loop: 1.38·10⁻⁵–6.39·10⁻⁴ per lift; no collective exponent is within 3.8 % of σ₀/k_e (the nearest
   at 3.88 %, k_e = 50), so the selection by proximity is unambiguous. Slow-time deviations: 4.437·10⁻³ and 1.166·10⁻³ (η = σ₀), 4.446·10⁻³ and 1.169·10⁻³ (η = ∞);
   ratio 3.8045 (k_e^{−1.9277}). The pair rates by `np.roots` equal the stored prediction to 10⁻¹².
3. **The instance.** The exact representation at (6, 5, 4)/15 has amplitudes (2.5, −2.1, 3.0), with max|f̂ − f| = 1.1·10⁻¹⁶ and
   max|f̂′ − f′| = 5.6·10⁻¹² on [−4, 4]. h_j(w) from the analytic Gram (4096 points) equals the raw record to 4.3·10⁻¹⁰ relative
   (`check_instance.py`).
4. **(a) and (b) exactly.** By the variational equation (`var_monodromy.py`): separations at σ₀/k_e to 1.2–3.1·10⁻¹², and the lifts'
   collective exponents equal the K₀-atom loop's to ≤ 1.8·10⁻¹⁰ (m = 15, 30; η = σ₀, ∞; k_e = 25, 50). The raw finite-difference exponents
   differ from the variational ones by ≤ 1.66·10⁻⁶ (separations) and ≤ 6.8·10⁻⁴ (collective, m = 30, η = σ₀, k_e = 50): "the resolution of
   the differences", as stated.
5. **(c) to k_e = 400.** The K₀-atom loop's slowest relative deviation is −37.5, −13.56, −3.84, −0.994, −0.251, −0.0628 % (η = σ₀) at
   k_e = 12.5, 25, 50, 100, 200, 400 (ratios 2.77, 3.53, 3.86, 3.97, 3.99), and −37.8, −13.69, −3.88, −1.005, −0.253, −0.0635 % at η = ∞.
   No complex multipliers.
   - *Uniform in η.* The λ^+ modes' deviation at k_e = 25 is 2.1·10⁻², 4.3·10⁻⁵ and 1.7·10⁻⁷ for η = 0.1, 1, 10 (ηT₀/k_e up to 2.5), while
     the λ^− modes keep 0.137.
   - *Uniform over W.* As w₃ → |A₃|/a_max = 0.2, h₁(w) falls from 3.1·10⁻⁴ to 8.9·10⁻⁶. The slowest relative deviation at k_e = 25 falls
     from −0.137 to −0.037, so the threshold does not degenerate (`var_boundary.py`).
6. **Resolution of the finite differences** (`fd_resolution.py`; K₀-atom loop, η = σ₀, k_e = 50). Against the variational monodromy:

   | step / rtol | max\|DP − M\| | slowest-exponent error (relative) |
   |---|---|---|
   | 10⁻⁵ / 10⁻¹³ (the experiment) | 3.0·10⁻⁸ | 4.1·10⁻⁵ |
   | 10⁻⁴ / 10⁻¹³ | 1.5·10⁻⁹ | 1.0·10⁻⁵ |
   | 10⁻⁶ / 10⁻¹³ | 3.2·10⁻⁷ | 6.6·10⁻⁴ |
   | 10⁻⁵ / 10⁻¹² | 1.3·10⁻⁸ | 5.6·10⁻⁵ |
   | 10⁻⁴ / 10⁻¹² | 7.8·10⁻¹⁰ | 6.8·10⁻⁷ |

   The integrator's error enters divided by the step, as the docstring says. A step of 10⁻⁴ would have been better still.
7. **The fixed prior at (6, 5, 4)/15** (`check_saddle_1d.py`). Θ*(w) has margins (−3.48275·10⁻³, +6.195·10⁻³, +7.169·10⁻³), fit error
   0.0455; the quoted −3.48·10⁻³ is class 1. At ϑ^e(w) every margin is σ₀ = 7.0744·10⁻³ (R = 0).
8. **The OT extension.** The plan is locally constant at the lift (verdict table; `check_ot_plan.py`).
9. **The proof's skew corrections** (`var_eta_mech.py` (B), `var_asym_ref.py`). Both terms of the second-order averaging are skew, as
   claimed, and on the asymmetric reference the second-order averaged matrix shifts the rates only at second order (9.8·10⁻⁶ and 2.45·10⁻⁶
   relative at k_e = 25 and 50).
10. **Every number of `adaptive_floor_2d.md`, the planar paragraph, the §7 clause and RESULTS** (`recompute_2d_md.py`).
    - Part A: h_j/σ₀ = 0.0042–1.4128; rel 6.007·10⁻⁷; law 0.30962 = Σ terms, against the bound 1.43800 (0.21531); largest share 0.1843.
      The lower bounds of W are (0.200, 0.150, 0.175, 0.150), sum 0.675.
    - Part B: 7 representable m-systems, all six references representable (margins 0.0031–0.050). Fixed floors: 3.246·10⁻⁴–4.178·10⁻³
      (representable) and 3.011·10⁻³–2.070·10⁻² (others).
    - Adaptive, representable: 1.030·10⁻⁷–4.107·10⁻⁵ (η = σ₀), 0.0003–0.0187 of the fixed prior's; 1.210·10⁻⁹–1.026·10⁻⁵ (η = ∞),
      ≤ 0.0060 of the fixed prior's. Adaptive, others: 0.0226–0.6414 of the fixed prior's. Falling between the last two records in 7/7 and
      11/11; the counts are conserved in all 72 rows.
    - Checks stage: margins 0.015519–0.043315 (m-systems) and up to 0.046533 with the references ("+0.016 to +0.047"); collective largest
      eigenvalues −4.97·10⁻² to −2.59·10⁻²; 24/24 stable.
    - Part C: fixed 1.0085, 0.9970, 0.9931, 0.9930; adaptive 1.0088, 0.9964, 0.9849, 0.9793; exact start 1.98–4.70·10⁻¹⁰, decreasing in
      every window.
11. **`planar_loop` against independent code** (`check_planar.py`).
    - Drift with weights and targets: 3.3·10⁻¹⁶ relative.
    - The averaged force equals k_e⟨G∇(eᵀΦ̃)⟩ along the orbit at e = r/k_e to 3·10⁻¹⁷ (the real planar law of `exp_ot_2d.py`); the orbit
      grid is the same 400 points.
    - n = 1 nearest-atom loop: equal to `MixturePriorLoop` to 8.7·10⁻¹⁸, with the same cells.
    - The weighted K-atom drift equals the replicated m-system's to 2.8·10⁻¹⁷.
    - The checks stage, with my Θ*(w): observed/exact bias 0.99999979–1.00000007.
12. **Reproducibility and Part C** (`check_partB_outside.py`, `check_partC.py`).
    - Three Part B systems are reproduced bitwise in their stored fit records (max difference 0.0): seed 1, m = 16, η = ∞; seed 4, m = 16,
      η = σ₀; seed 1, m = 64, η = ∞.
    - Part C's exact start: the stored e² series is reproduced bitwise on its first 1600 samples (t ≤ 160).
    - Its five-period rms is 2.90·10⁻¹¹, 4.70·10⁻¹⁰ and 7.70·10⁻⁹ at dt = 10⁻³, 2·10⁻³ and 4·10⁻³ (first window). The ratios 16.2 and 16.4
      per halving are RK4's dt⁴: the integrator's level. The parameters stay within 8·10⁻⁹ of the lift at dt = 2·10⁻³.
13. **Outside W: slow convergence toward the box-constrained minimum, not a stall** (`check_partB_outside.py`). In all three runs the counts
    are conserved, the class diameters are 0, the amplitude coordinates run off like ln s, and L̄ decreases at every record toward the
    box-constrained minimum of the weighted K₀-atom loss at the sampled weights.

    | system | margin | fit error at s = 6000 → end | max\|α\| | box-constrained minimum | L̄ above it at 6000 / end |
    |---|---|---|---|---|---|
    | seed 1, m = 16, η = ∞, (3, 6, 1, 6) | −0.112 | 0.25445 → 0.25211 (s = 6·10⁴) | 11.5 → 15.9 | 0.25180, three amplitude components at ±a_max | 2.1 % / 0.25 % |
    | seed 4, m = 16, η = σ₀ | −0.075 | 0.10628 → 0.10102 (s = 3·10⁴) | 8.9 → 11.7 | 0.09944, one bound active | 14 % / 3.2 % |
    | seed 1, m = 64, η = ∞ | −0.019 | 0.04584 → 0.03971 (s = 2.4·10⁴) | 8.4 → 10.6 | 0.03662, one bound active | 57 % / 18 % |
14. **The finite-η branch of `planar_loop.run`** (`check_finite_eta.py`). Against DOP853 over 300 slow units (seed 5, m = 16, η = σ₀) the
    errors are 1.60·10⁻⁷, 3.99·10⁻⁸ and 9.96·10⁻⁹ at dt = 0.25, 0.125, 0.0625 (ratios 4.02, 4.01: second order). The cells are unchanged.
15. **Tests' quantities** (`check_tests.py`). Floquet test as in D3m-6. Planar test 1: 5.6·10⁻¹⁷. Planar test 3: margin +0.050; lift fit
    1.9·10⁻¹⁶; drift 6.4·10⁻¹⁸; after 50 slow units 1.9·10⁻¹⁶ and |Z − e| = 2.2·10⁻¹⁶; the fixed prior's fit 0.101 > 10⁻². Planar test 4
    as Part A.

## Not checked

- The experiments were not rerun in full: the twelve finite-difference monodromies, the 36 Part B systems × 3 priors, and Part C's fixed and
  adaptive runs. Single systems were reproduced bitwise, the regenerations are byte-identical, and the monodromies were recomputed
  independently by the variational equation.
- The O(k_e⁻³) remainder's constant was not derived analytically. The filtered-Gram estimate and the fits are numerical.
- The nonlinear basin, uniform in m, is open by the authors' account and was not attempted.
- The Floquet spectrum of planar perfect tracking was not computed. The skewness was tested on a one-dimensional asymmetric reference, not
  on the Lissajous reference, where ⟨N⟩ ≠ 0 as well.
- The OT-leak extension was checked through the local constancy of the plan, not by its own monodromy.
- The literature statements of D3m-3 (Narendra–Annaswamy 1987; Magnus 1954; Blanes et al. 2009; Kosut–Anderson–Mareels 1987) are from
  memory; no full text was fetched.
- README's sentence on D‴.26 and the concurrent session's untracked files (`exp_gauge_floor.py`, `sliding.py` and their tests) are outside
  the material and were not reviewed.
