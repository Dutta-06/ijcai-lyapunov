# Referee report D3n (round 19) — PARADIGM §5d: Proposition D‴.27 ((i)–(iv), proof, Status, *Numerics*), the §0 and §7 rows "| D‴.27 |", the §8 bullet "Done (Proposition D‴.27)" with the edits in branch 7 and in Theorem D‴.23(i), the a-priori-bound column of Remark D‴.25, RESULTS (Proposition D‴.27), the LITERATURE audit bullet, the `box` option of `src/core/adaptive_prior.py`, `tests/test_adaptive_prior.py::test_a_projected_prior_keeps_the_solutions_bounded_and_is_inactive_inside_its_box`, `experiments/exp_adaptive_prior.py` Part H → `results/chaos/adaptive_prior.md`; and the application of referee round 18 (D3k-1 … D3k-11, D3l-1 … D3l-8)

**Setup.** Branch `claude/festive-bardeen-eopej2` at 442cd6c (HEAD throughout; `git log --oneline -5`: 442cd6c, d2092a5, 4e2e938, b76b3b9,
1aeaa56).
- At the start and at 20:58 UTC the working tree had no change to tracked files. A concurrent session had five untracked files
  (`experiments/exp_gauge_floor.py`, `results/chaos/gauge_floor_raw.npz`, `src/core/sliding.py`, `tests/test_gauge_floor.py`,
  `tests/test_sliding.py`). None of them was used or reviewed.
- The material was reviewed from 442cd6c extracted with `git archive` into `ref_D3n/tree/` (analysis; every script imports this tree,
  never the repository). A second copy, `ref_D3n/regen/`, was used for `--from-raw`.
- No repository file was modified, and no git command touching the index or the tree was run.
- Python ran with `PYTHONDONTWRITEBYTECODE=1` and one BLAS thread, one job at a time. pytest ran with `-p no:cacheprovider` on the
  scratch tree.

Line numbers are 442cd6c's:
- `PARADIGM.md`: Proposition D‴.27 l. 3099–3157 ((i) l. 3106, (ii) l. 3109, (iii) l. 3114, (iv) l. 3118, Proof l. 3123, *Status* l. 3135,
  *Numerics* l. 3148); §0 row l. 125; §7 row l. 5589; §8 branch 7 l. 5785–5790; §8 bullet l. 5836–5840; Theorem D‴.23(i) l. 2648;
  Proposition D‴.24 l. 2818–3019 (its *Numerics* l. 2946–3004); Remark D‴.25 l. 3021–3046; Remark D‴.10's pointer l. 1766–1767.
- `results/RESULTS_PARADIGM.md`: D‴.24 l. 554–608, D‴.27 l. 623–637. `LITERATURE_TRANSPORT.md`: l. 351–398.

Scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3n/scripts/`;
each `.out` is one level up, in `ref_D3n/`.

| script | what it checks |
|---|---|
| `from_raw.out`, `tests.out`, `collect.out` | `--from-raw` of both experiments (byte-identical); the 16 tests of the two files; the collection count |
| `common.py`, `inspect_raw.py` | shared setup (the archived tree); the raw files' structure |
| `partH_rerun.py` | the experiment's own `projected_job()` against the stored entry H |
| `check_LG.py` | L_ϑ^G: the closed form against the true supremum (own derivatives, 400-start search); max\|f\|, B, B/σ₀; the run's own largest learning force |
| `check_bound.py` | (i) on twelve m-particle systems (OT plan at m = 15 and 16, cells with and without an empty cell, η = σ₀ and ∞, two targets beyond the saturation, starts 18 outside the box): the comparison inequality at every step, the bound with B and with the run's own force, the leak targets and the prior in 𝒵 |
| `check_unbounded_Z.py` | (i) with a closed convex but unbounded 𝒵 |
| `check_energy.py` | (ii): 𝒥 at every step of the projected midpoint rule; the rule at three dt against an event-driven reference of the projected flow; LSODA at η = ∞ against DOP853 with an event at the kink; 𝒥 at every refresh at m = 15 |
| `check_stall.py` | the projected midpoint rule near a face at six dt (the code's own `run()`), and the mechanism |
| `check_kkt.py` | the end state's KKT condition and balance; the face minimum; (iv)'s identification; Part G(a)'s saddle with the box |
| `check_iv.py`, `check_iv_G4.py`, `check_between.py` | (iv): interior critical points of the K₀-atom loss; m = 15 from a spread start; the four-bump target with the box to 10¹⁰; the limit against the box's amplitude bound |
| `check_partHb.py`, `check_test10.py` | Part H(b) instrumented (clip and tangent activations, the margin to the faces); test 10's quantities and mutations |
| `mymodel.py`, `check_model.py` | an independent implementation of the instance (analytic first and second derivatives, Newton with continuation, D‴.5 margins, K-atom Hessian, exact representations, the three laws) against finite differences and the repository |
| `check_sampled.py` | Round 18's per-sample checks from the stored Part B finals: margins at all 16 vectors, exact bias, delta method, D, W² split (own OT, HiGHS) |
| `check_laws.py` | `bias_floor_law` with four bases; the weighted and unweighted shares; D‴.24(iv)'s boundary constant; test 9's run-off law at three η |
| `check_prob.py`, `check_prob240.py` | Part E: exact probabilities, Hoeffding, permutations, the fixed prior's saddles, the constrained minima on the branch and globally |
| `check_mc.py`, `check_mc_tail.py` | a Monte Carlo with my own seeds (200 nested samples); the heavy tail of the adaptive prior's coupled distance |
| `check_partC.py`, `check_partF.py`, `rerun_EFG.py` | Part C's fit errors at s = 1250, 2000, 3000 per seed; Part F's draw statistics; Parts G, E and two Part F jobs re-run against the raw file |

The long runs (exp_adaptive_floor.py Parts B–D and `--only P`, exp_adaptive_prior.py Parts C, D and 38 of the 40 Part F jobs) were not
rerun in full. Every number below comes from the scripts listed.

## Verdicts

### (A) Proposition D‴.27

| object | verdict |
|---|---|
| **(i)** the a-priori bound | **the distance bound holds; "every solution is bounded" needs 𝒵 bounded (D3n-1, minor).** Re-derived. ½d² is C¹ with G⁻¹-gradient ϑ − Πϑ. The projection inequality in the G⁻¹ inner product gives ⟨ϑ − Πϑ, ϑ − z⟩_{G⁻¹} ≥ d² for z ∈ 𝒵, and Cauchy–Schwarz gives d′ ≤ B − σ₀d where d > 0. The comparison needs only the region d > B/σ₀ > 0, so d = 0 causes no difficulty; the same comparison gives the sharper d(s) ≤ B/σ₀ + (d(0) − B/σ₀)₊e^{−σ₀s}. For a box and diagonal G the G⁻¹-projection is separable, i.e. the coordinate-wise clip, and the tangent-cone projection zeroes the outward components (the code's `project` and `tangent`). \|f̂\| ≤ a_max holds: f̂ = (1/m)Σa_iφ_i with \|a_i\| < a_max and 0 < φ_i ≤ 1. The closed form 4.658 is a valid upper bound on L_ϑ^G, not the supremum, which is 4.0435 (D3n-5). In twelve stress runs, at every step: the leak targets lie in 𝒵 (≤ 5.3·10⁻¹⁵ off it for the split rows at m = 16, the barycentric projection being a convex combination); the comparison inequality holds with the code's actual drift (largest residual 0.0); the bound with the run's own largest force holds with margins of 11 to 97. The runs cover the OT plan at m = 15 and m = 16, the cells with and without an empty cell (its atom frozen), η = σ₀ and ∞, and starts 18.0 outside the box. |
| **(ii)** the energy | **holds; the convergence is formal, as stated, and can be proved at η = ∞ (D3n-7).** The prior's part of d𝒥/ds is −‖Π_T(−grad)‖² by Moreau's decomposition in the G⁻¹ metric, which is the gradient flow's metric up to the factor σ₀q_j/η. At η = ∞ the clipped barycentre minimises Σ_iγ_ij‖ϑ_i − z‖²_{G⁻¹} = q_j‖z − T̂_j‖² + const over the box. 𝒥 never increases in the discrete runs: 12001 steps of the K₀-atom projected midpoint rule (largest increment −4.2·10⁻⁹), and 4001 refreshes at m = 15 for η = σ₀ and ∞ (−6.2·10⁻⁹ and −3.8·10⁻¹⁰). At Part H's end state the KKT condition holds: only (atom 1, α) is active, with T̂ − z = +0.8389 in the normal cone and the other coordinates at z = T̂. The particles balance to 2.4·10⁻¹⁷. The characterisation is that of a unique plan. |
| **(iii)** nothing changes inside | **holds, trivially.** For the local-convergence part, "proved" inherits Theorem D‴.23(v)'s "sketched" (D3n-7). Part H(b) tests only that the box is not reached (D3n-6). |
| **(iv)** the price outside | **the identification holds exactly**: the code's η = ∞ drift minus the learning term equals −σ₀(ϑ − Πϑ) to 0.0, and the gradient of (σ₀/2)dist²_{G⁻¹} to 6.3·10⁻¹¹. **"Stops on the boundary" is not implied, and "between" is an observation (D3n-3, minor).** Confirmed on a second target, which the paper does not show: on Part G(a)'s four-bump target the box stops the width coordinate's run-off at the face κ = 6. The novelty and the "continuous switching σ" wording are under D3n-4 (minor). |
| Proof | **correct where it is a proof.** Missing: compactness (D3n-1); the convexity step for split rows; Z(0) ∈ 𝒵. |
| *Status* | (i) "proved" needs compact 𝒵; (iii) partly sketched; the novelty is overstated (D3n-4); the test description claims an unasserted check (D3n-6). |
| *Numerics* (a), (b) | **every number reproduces bitwise** (Part H re-run in 82 s: all eleven rows, max_dist 0.9165 and 0.8880, B = 26.0786, B/σ₀ = 3686.3, diff_b 0.0). (a)'s η = σ₀ result depends on dt = 0.25 happening to reach the face (D3n-2, minor). (b) holds by construction (D3n-6). "Fit error 0.1415" is 0.141545; "the unprojected limit 0.1050" is the face minimum, 0.104951. |
| §0, §7 rows; §8 bullet and branch-7 edit; D‴.23(i) edit; Remark D‴.25's column; RESULTS; LITERATURE; README | **numbers hold.** Wording: "closed convex" (D3n-1), "stops on the boundary" / "between" (D3n-3), novelty (D3n-4), L_ϑ^G and B (D3n-5), "the projection is inactive" (D3n-6). The column is consistent with D‴.27, but the Reading and four other sentences were not updated (D3n-8). |
| `box` option of `adaptive_prior.py` | **`project`, `tangent` and the η = ∞ drift are correct**, and the `box=None` path is unchanged (`git diff d2092a5 442cd6c`: identity maps). **`run()`'s projected midpoint rule can stall the prior inside 𝒵 (D3n-2).** At η = ∞ the field is Lipschitz with a kink, not discontinuous; LSODA agrees with DOP853 using an event at the kink (s = 74.1308) to ≤ 1.0·10⁻⁹. |
| test 10 | **passes; its two boundedness parts discriminate.** Without the box: at η = ∞ the distance is 11.74 and ϑ_a 17.62; at η = σ₀ ϑ_a is 9.02 and the prior leaves the box. "The prior stays in its box" is the clip, and the third part holds by construction (D3n-6). Measured: distance 0.8874 (< 2), max ϑ_a 6.776 (< 8), drift 5.7·10⁻¹⁷; at η = σ₀ max ϑ_a 6.805 and no 𝒥 increase. |

### (B) The application of round 18

| object | verdict |
|---|---|
| Proposition D‴.24 as rewritten ((i)–(vi), proof, Status, *Numerics*) | **holds as revised.** Every per-sample number is reproduced by my independent model and OT: the 16 margins, the exact-bias ratios, the delta ratios, D, and the W² split. The laws are 6.166158·10⁻², 21.32109 and 46.32494, and the bound 8.193011·10⁻² = (σ₀/4)·46.32494. The boundary constant is confirmed: s·‖f̂ − f‖² = 0.4017 and 0.8054 at s = 10⁶. Readings: D3n-9, D3n-10, D3n-11 (cosmetic). |
| Remark D‴.25 (split columns, Reading) | **holds.** "The mass constraint removes the function- and measure-level floors" holds at commensurate m: OT/fixed W² ≤ 3.5·10⁻¹⁵ and floor ≤ 2.4·10⁻²¹. "The adaptive prior also removes the splitting instability" holds at the exact representations: margins σ₀, and every adaptive system ends replicated, including at the five saddle samples. The bound's identity with the coupled-distance law is exact. The Reading's a-priori sentence and "at every weight" need updating (D3n-8). |
| the checks stage of `exp_adaptive_floor.py` (post_checks, enumerate_job, mc_job; the report's sections, Parts E and F) | **reproduces.** Regeneration is byte-identical. The probabilities are exact, and the branch minima are reproduced to 0 relative difference. The "constrained minimum" is a branch minimum, not the infimum (D3n-10). The adaptive prior's coupled distance has a heavy tail (D3n-11). |
| `exp_adaptive_prior.py` Parts E, F, G as recomputed | **re-run bitwise.** Part G: saddle, flows, continuation κ = 2.159 … 9.751, face minimiser 1.943283·10⁻² with slope −1.3257·10⁻⁶, unbounded solutions, relabelling. Part E: records and crossings 41564.5186, 47796.9288, 79714.1458. Part F: two jobs. Part F's statistics recomputed: sd ratios 3.33–5.96, mean differences +0.86 to +1.66 SE, growth factors 23.77–25.03. |
| `src/core/floor_bias.py` | **fixed.** Default, explicit e_j − e_K, orthonormal and skewed bases all give 6.166158·10⁻². |
| `tests/test_adaptive_floor.py` (6), test 9 | **pass.** Test 9's law re-derived: (a_max − a)s/law = 0.9864 and F/F_∞ = 1.0032 at s = 10⁵, η = σ₀; at s = 10⁶ the ratio is 0.9969, 0.9982 and 0.9987 at η = σ₀/10, σ₀, ∞. A bounded limit would need F = 0, so the test discriminates. |
| §0/§7 rows of D‴.23 and D‴.24, §8 bullets, RESULTS, LITERATURE ([S49], [S51], audit bullets), README | **numbers hold** (README's 129 tests: 129 collected). Wording: D3n-8 to D3n-11. |

The round-18 fixes, item by item:

| id | fixed as logged? | correct and complete? |
|---|---|---|
| D3k-1 | yes: margins (matching and misplaced prior), diameters and sub-clusters, the exact bias against the observed floor, the delta method; the two admissibility notions throughout; (B) for D‴.5-admissible samples; the margins σ₀ at ϑ^e(w); Part D's saddle orbit; test 4 at (5, 4, 6)/15 | correct: margins and ratios recomputed independently, identical to the printed digits |
| D3k-2 | yes: two columns, the three faces, D and the W² split printed | correct (D and W² recomputed); the Reading's a-priori sentence was not updated for D‴.27 (D3n-8) |
| D3k-3 | yes: "the law's value", "at most", the Monte Carlo | correct (my Monte Carlo: 2.82 ± 0.33, 1.51 ± 0.18, 1.11 ± 0.11); the reading of the adaptive coupled distance overreaches (D3n-11) |
| D3k-4 | yes: Part E, the wordings | the probabilities are exact; the floors are branch minima, not "the infimum" (D3n-10) |
| D3k-5, D3k-6, D3k-7, D3k-8 | yes | correct (shares recomputed; test 6 checked) |
| D3k-9 | yes: the pointer, the fit errors at s = 1250 and 2000, the gate | the s = 1250 reading is wrong at m = 15 (D3n-9) |
| D3k-10 | yes, every item | correct (`bias_floor_law` basis-independent) |
| D3k-11 | yes | correct |
| D3l-1 | yes: the continuation to 10¹⁰, the face minimiser, the md, the numerics, RESULTS, the round-17 correction | correct (re-run bitwise) |
| D3l-2, D3l-4, D3l-5, D3l-6, D3l-7 | yes | correct |
| D3l-3 | yes: five draws, means, ranges, standard errors | numbers correct; "of the same expected size" should be "consistent with" (D3n-11) |
| D3l-8 | yes | correct (the law re-derived and confirmed at three η) |

## Confirmed issues

**D3n-1 (minor; D‴.27 l. 3099 "Let 𝒵 ⊂ ℝ^p be a closed convex set"; (i) l. 3106–3108 "every solution is bounded"; (ii) l. 3110; the
§0 row; the branch-7 edit l. 5789 "projecting the prior onto a convex set restores the bound"; the §8 bullet l. 5836 "Confining the
adaptive prior's atoms to a closed convex set keeps every particle within B/σ₀ of it"; LITERATURE l. 392–398).**

(i)'s distance bound holds for every closed convex 𝒵, but the conclusion "every solution is bounded" needs 𝒵 bounded. Otherwise
dist(·, 𝒵) ≤ const says nothing: 𝒵 = ℝ^p is Theorem D‴.23's unprojected prior, which runs off.

Evidence (`check_unbounded_Z.py`):
- Take 𝒵 = ℝ × [−4.417, 4.417] × [−6, 6], the instance's box without its amplitude bounds (closed, convex, unbounded), and Part G(b)'s
  K₀-atom system at η = ∞.
- The particles stay in 𝒵 (distance 0.0), while ϑ_a = 6.200, 9.128, 13.131, 17.62, 22.213 at s = 10², …, 10⁶.
- That is Part G(b)'s unprojected run to the printed digits.

(ii)'s convergence of "the bounded solutions of (i)" and the §8 bullet rest on the same hypothesis. Also implicit is Z(0) ∈ 𝒵^{K₀}: the
projected flow is undefined outside 𝒵, and the code clips at the first step.

*Fix:*
- Write "a compact convex set 𝒵 (a bounded box on the instance), with Z(0) ∈ 𝒵^{K₀}" in the statement, the branch-7 edit, §0, §8 and
  the audit bullet.
- Optionally state the sharper form of the same comparison: d_i(s) ≤ B/σ₀ + (d_i(0) − B/σ₀)₊e^{−σ₀s}.

**D3n-2 (minor; `src/core/adaptive_prior.py`: `run()` l. 159–170 with `drift` l. 112 and `tangent` l. 84–92 — the projected midpoint
rule; D‴.27 *Numerics* l. 3150–3152 "the projected midpoint rule, dt = 0.25, to s = 10⁵ … converge to the same bounded stationary point";
the md's (a); the η = σ₀ part of test 10).**

The midpoint stage evaluates the prior's velocity at Z_mid = Z + (dt/2)k1Z and projects it onto the tangent cone *at Z_mid*. Suppose the
base point Z lies less than (dt/2)\|k1Z\| below a face:
- Z_mid is on or past the face, so `tangent` zeroes k2Z's outward component;
- then Z_new = Π(Z + dt·k2Z) = Z, and the atom stays strictly inside 𝒵;
- the same happens at every later step, because the outward pull persists.

The scheme then has a spurious fixed point that is not a stationary point of the projected flow: its own drift stays at 5.9·10⁻³.

Evidence 1 — the code's own `run()` on Part H(a)'s system, η = σ₀, to s = 20000 (`check_stall.py`):
- dt = 0.25 and 0.05 reach the face (gap 0; drift 2.6·10⁻⁷ at 20000, still falling; at dt = 0.25 the md's 1.75·10⁻¹⁵ by 10⁵).
- dt = 0.2, 0.125, 0.1 and 0.0625 stall 2.405·10⁻⁴, 1.171·10⁻⁴, 2.399·10⁻⁴ and 1.169·10⁻⁴ below the face, with the code's drift
  5.94·10⁻³ at s = 20000.
- The mechanism at dt = 0.125, step 2381 (s = 297.625): z_a is 1.171·10⁻⁴ below the face, k1Z_a = 4.917·10⁻³, the midpoint is
  1.902·10⁻⁴ past the face, k2Z_a = 0, and z_a is unchanged.

Evidence 2 — against an event-driven reference of the projected flow (DOP853 at rtol 10⁻¹², free until s = 297.6488 and then held;
`check_energy.py`):
- dt = 0.25 agrees to 3.5·10⁻⁶, 6.3·10⁻⁷, 2.2·10⁻⁸ and 7.0·10⁻⁹ at s = 100, 300, 1000 and 3000;
- dt = 0.125 and 0.0625 are 1.2–1.8·10⁻⁴ off from s = 300 on, and stay off.

So Part H's η = σ₀ numbers are right only because, at dt = 0.25, the full step's clip happened to land on the face before a midpoint did.
That happened at two of the six steps tried. The scheme is not a consistent discretisation of the projected flow at a face. The discrete
𝒥 does stay non-increasing.

*Fix:*
- Project positions, not velocities, inside the step: Z_mid = Π(Z + (dt/2)v(Z)) and Z_new = Π(Z + dt·v(Z_mid)), with v the unprojected
  velocity. Alternatively, take the tangent cone at the base point Z.
- Assert in Part H and in test 10 that the end drift is ≤ 10⁻¹² (a stall gives 5.9·10⁻³), and run Part H(a) at a second dt.

**D3n-3 (minor; (iv) l. 3118–3120 "Where the target has no exact representation with atoms in 𝒵, the projected prior stops on the boundary
of 𝒵 … the limit keeps a bias, bounded by (i)"; the §0 row "outside the prior stops on the boundary and a bounded bias remains"; the §8
bullet; *Numerics* l. 3153–3154, the §0 and §7 rows and RESULTS "between the fixed prior's 0.1750 and the unprojected limit 0.1050"; the
md's (a) "the fit error lies between …").**

(a) *"Stops on the boundary" does not follow from (ii).* A limit whose prior is interior is a stationary point of Theorem D‴.23(i), not
necessarily an exact representation. Non-replicated stationary points with an interior prior exist where no exact K₀-atom representation
does (`check_kkt.py`):
- Part G(a)'s saddle on the four-bump target, with the instance's box: its class barycentres lie at least 1.986 inside every face.
- The plan is the class partition, and the projected joint drift is 4.4·10⁻¹⁷ at η = σ₀ and at η = ∞.
- It is a saddle, so only a non-generic set of starts converges to it. But (iv) is stated for every solution, and (ii) says nothing about
  which critical point is reached.

(b) *"Between" is a property of one run.*
- The K₀-atom limit does stay between 0.1750 and 0.104951 for amplitude bounds 0.7, 0.75, 0.8, 0.9, 0.95 and 0.99 of a_max: fit errors
  0.1736, 0.1676, 0.1611, 0.1415, 0.1268, 0.1101 (`check_between.py`). The box is a design trade-off.
- But an m = 15 system from a spread start ends at s = 20000 at non-replicated configurations (`check_iv.py`; every second particle's
  amplitude coordinate +20, some centre coordinates −9):
  - at η = σ₀, fit error 0.0534 — below the "unprojected limit" 0.1050. Class 3 is split into sub-clusters at c = −0.58 and 1.33, all
    three prior atoms are on the α-face, the drift is 2.2·10⁻⁸, and 𝒥 = 6.29·10⁻³ against 8.66·10⁻³ for the replicated configuration;
  - at η = ∞, fit error 0.1747 (two classes split), at the fixed prior's level.

(c) *Positive.* On Part G(a)'s target the box does stop the run-off (`check_iv_G4.py`):
- from the exact split fit at η = ∞, the merged atom's width coordinate passes 2.159, 3.843 and 5.529 at +10⁶, 10⁷ and 10⁸;
- it reaches the face κ = 6 by 10⁹ and stays: drift 4.5·10⁻¹⁷, fit error 1.943376·10⁻²;
- the unprojected run has κ = 9.751 at 10¹⁰.

This is a second confirmation of the restoration, not in the paper. (The m = 3 K₀-atom loss on Part G(b)'s target showed no interior
critical point inside the box in 300 BFGS starts.)

*Fix:*
- (iv): "If the limit has an atom on ∂𝒵 — necessarily so when no stationary point of Theorem D‴.23(i) has its prior in int 𝒵 — the
  particles of that class leak toward the face …".
- Call "between" an observation for the K₀-atom run from the exact atoms, with the box-size dependence.
- Add Part G(a) with the box to Part H.

**D3n-4 (minor; *Status* l. 3138–3142 "New is this combination: the adaptive prior keeps the robustness half of the classical trade-off";
the §0 row's last column; LITERATURE l. 392–398 "a new combination that restores the robustness half"; (iv) l. 3120–3121 "a continuous
form of the switching σ-modification [4], whose dead zone is 𝒵").**

(a) *The classical switching σ-modification already has both halves.* Its leak is off inside a ball \|θ\| ≤ M₀ and on outside. The
estimates are therefore bounded, and the adaptive law keeps its ideal (unbiased) properties when the true parameter lies in the ball. A web
search returned this from the abstract of a 1986 IEEE CDC paper of Ioannou's on robust adaptive laws: "in contrast to the dead-zone,
fixed-σ, and ε1-modifications, the switching-σ approach guarantees robustness without destroying the ideal properties of the adaptive
control schemes" ([record](https://gnosis.library.ucy.ac.cy/handle/7/56943?show=full)). The full text could not be reached.

D‴.27(iv) itself shows that, on collapsed classes at η = ∞, the construction is such a law: −σ₀(ϑ − Π_𝒵ϑ), i.e. σ_s(ϑ) = σ₀(1 − M₀/\|ϑ\|)₊
for a ball. So the combination claimed as new is the switching σ's defining property. What is specific here is its transfer to a learned
centre of a mass-constrained K₀-atom prior at finite η, and its compatibility with Theorem D‴.23's bias removal.

(b) *The exterior penalty.* There, the leak is the gradient of the exterior quadratic penalty (σ₀/2)dist²(ϑ, 𝒵), and "the bias outside" is
the penalty's constraint violation ‖force‖/σ₀. At Part H's end state the particle is 0.8389 outside 𝒵, and its leak force
σ₀·0.8389 = 5.935·10⁻³ balances the learning force (`check_kkt.py`).

(c) *"The learning law keeps its smooth gradient form."* The particles' law does. But at finite η the joint (Θ, Z) law is a projected
dynamical system whose prior velocity jumps at the faces; the code's docstring says so ("the projected flow is discontinuous at the
boundary"). That is the regularity cost of parameter projection, moved to the prior. The η = ∞ form is Lipschitz but not C¹: it has a kink,
crossed at s = 74.13 in Part H.

(d) *Wording.* The switching σ of [4] is already continuous (made so to keep the law Lipschitz), so "a continuous form of" misleads. And
"dead zone" is the name of a different classical modification, in which adaptation is switched off for small errors.

*Fix:*
- Restrict "new" to the finite-η coupling (a projected adaptive law for the prior's atoms, the mass constraint, K₀ classes) and to its
  compatibility with Theorem D‴.23's bias removal.
- Cite the switching σ's retention of the ideal properties ([4] and the 1986 paper) and the exterior quadratic penalty.
- Write "a switching σ-modification toward 𝒵, the leak off on 𝒵".

**D3n-5 (cosmetic; D‴.27 *Numerics* l. 3155 "(B = 26.08, L_ϑ^G = 4.658)"; the statement l. 3103–3105 "B := sup_ϑ …"; the §7 row;
RESULTS l. 632–633; the md's (a); the `projected_job` docstring).**

L_ϑ^G is defined as sup_{x,ϑ}‖G^{1/2}∇_ϑΦ̃‖. The code's closed form is
√(G_α·1 + G_γ(a_max e^{−1/2}/s_min)² + G_κ(a_max(s_max − s_min)/(2e·s_min))²) = √(1 + 0.25·8.0871² + 0.25·4.1693²) = 4.657893.

It takes the three suprema separately, but they are not attained together:
- the α-term needs a = 0;
- the κ-term needs s′ at its maximum, i.e. s = 1.15;
- the γ-term needs \|a\| → a_max and s = s_min.

The supremum is 4.0435 (400 Nelder–Mead starts over x ∈ [−2, 2] and ϑ; 4.0399 on the ρ-grid). It is approached by the γ-term alone, and it
is the "4.044 by search" of Theorem D‴'s numerics (l. 3168). So:
- 4.658 is an upper bound, and 26.08 bounds B from above;
- with the true constant, B ≤ 22.64 and B/σ₀ ≤ 3200;
- B is a supremum over the configurations Θ as well, since r depends on Θ;
- the run's own largest learning force is 0.0822 (η = ∞; 0.0859 in the m = 15 runs), 3.2·10⁻³ of the bound, i.e. B_run/σ₀ = 11.6–12.2
  against the observed 0.89–0.92 (`check_LG.py`, `check_bound.py`).

*Fix:* write "L_ϑ^G ≤ 4.658 (closed form; 4.044 by search), B ≤ 26.08", and "sup over ϑ and Θ".

**D3n-6 (cosmetic; *Numerics* (b) l. 3155–3157; the md's (b) "(the run stays inside 𝒵: True)"; the §0 and §7 rows and RESULTS "the
projection is inactive"; test 10's third part; the *Status*'s test description l. 3143–3146 and RESULTS "(the unprojected coordinate passes
11 by s = 10⁴)").**

(a) *Part H(b)'s identity is guaranteed by the implementation* whenever the prior never reaches a face: `project` is the identity inside,
and `tangent` acts only within 10⁻¹⁴ of a face. Instrumented over the 6000 slow units (`check_partHb.py`):
- the clip changed Z 0 times, and the tangent-cone projection 0 times;
- the prior stays at least 1.991, 2.562 and 3.226 from the α-, γ- and κ-faces (17%, 29% and 27% of the box widths), the particles at least
  1.989, 2.756 and 3.225.

So it is a code-path check, not a test of (iii), whose content is a tautology at that level. The md's "stays inside 𝒵: True" is evaluated
at the end state only.

(b) *In test 10* (`check_test10.py`):
- "the prior stays in its box" at η = σ₀ is enforced by the clip;
- the identical-runs part fails only if the box cuts the path (a γ-bound of 1.80, below the shifted prior's 1.855: max\|Δ\| 7.8·10⁻³);
- the boundedness assertions do discriminate (see the verdict table).

(c) *An unasserted check.* The *Status* and RESULTS list "the unprojected prior's coordinate passes 11 by s = 10⁴" among what the test
checks. It is a comment; no test asserts it. It is true: 13.13 at η = ∞ and 11.79 at η = σ₀ in Part G(b).

*Fix:*
- Print the margin to the faces.
- Describe (b) as "the box is never approached (margin ≥ 17% of its width), so the runs coincide by construction".
- Assert "passes 11" or drop it from the test descriptions.

**D3n-7 (cosmetic; *Status* l. 3135–3137; (ii) l. 3109–3113; Proof (i) l. 3123 "z_{π(i)} ∈ 𝒵"; Proof (ii) l. 3131–3134).**
- (a) *(iii)'s label.* Its second sentence — "the bias removal of (iii) and the local convergence of (v) hold unchanged …" — inherits
  Theorem D‴.23(v)'s "sketched". Only "the projected flow is Theorem D‴.23's while the prior is interior" is proved.
- (b) *The energy of the nearest-atom variant.* (ii) covers the nearest-atom prior of D‴.24(ii), but names only 𝒥 and the metric of
  D‴.23(i). There the energy is 𝒥_∘ and the metric uses the cell masses (degenerate at an empty cell).
- (c) *Ties.* The critical-point characterisation z_j = Π_𝒵(T̂_j) is that of a unique plan. At a tie the Clarke-critical set is the
  convex hull of the active pieces', as in D‴.17(iii).
- (d) *Split rows.* For incommensurate m a split row's target is mΣ_jγ_{ij}z_j, in 𝒵 by convexity — the step the proof of (i) skips. It is
  checked: ≤ 5.3·10⁻¹⁵ off the box at m = 16.
- (e) *η = ∞ needs no indicator.* The reduced energy is the minimum over the vertex plans of L̄ + (σ₀/2)Σ_j[within-class part +
  q_j dist²(T̂_j, 𝒵)]. Each piece is C^{1,1} (dist² of a convex set) and definable (the box is semialgebraic). So Proposition D‴.17's
  argument applies as it stands to the bounded solutions, and the "formal" label can be lifted at η = ∞.

*Fix:*
- Write "(iii) proved; the local convergence as Theorem D‴.23(v) (sketched)".
- Name 𝒥_∘ and the cell-mass metric in (ii), and add the tie case.
- Add the convexity remark to the proof of (i).
- Write "the convergence at η = ∞ as Proposition D‴.17, at finite η formal".

**D3n-8 (cosmetic; Remark D‴.25's Reading l. 3036–3042; D‴.24(v) l. 2895–2896; Remark D‴.10's pointer l. 1766–1767; LITERATURE l. 361 and
l. 379).**

The a-priori-bound column now says "no (yes when projected, D‴.27)". These sentences still say the adaptive prior gives the bound up:
- the Reading: "every adaptive prior pays with the rate … and with the leak's a-priori bound on the parameters (Theorem D‴.23(i))";
- D‴.24(v): "The price is the collective rate …, the leak's a-priori bound, and a larger coupled distance";
- the D‴.10 pointer: "at the price of … the leak's a-priori bound";
- the D‴.23 audit bullet: "whose robustness half (the leak's a-priori bound) the adaptive prior gives up";
- the D‴.24 audit bullet: "the leak's a-priori bound".

The Reading also says the adaptive prior removes the function-level floor "by removing the bias at every weight". That holds at every
representable weight; outside W the sample keeps an error, by (iii).

*Fix:* add "(restored by projecting the prior, D‴.27)" in these places, and write "at every representable weight".

**D3n-9 (cosmetic; D‴.24 *Numerics* l. 2987–2988; RESULTS l. 580–581 — introduced by the D3k-9 fix).** The text says "at s = 1250 the
ungated nearest-atom adaptive prior leads (1.75·10⁻², 3.8·10⁻³, 2.8·10⁻³ against 9.2–9.8·10⁻³)". From the raw records
(`check_partC.py`):
- at m = 15 its seed mean, 1.746·10⁻², is above both OT designs' (9.776·10⁻³ fixed, 9.167·10⁻³ adaptive);
- per seed at m = 15 it leads in seeds 0 and 2 (5.78·10⁻³, 3.35·10⁻³) and trails in seeds 1 and 3 (2.04·10⁻², 4.03·10⁻²);
- it leads at m = 60 and 240 only.

*Fix:* write "at s = 1250 the ungated nearest-atom adaptive prior leads at m = 60 and 240 (3.8·10⁻³ and 2.8·10⁻³ against 9.6–9.8·10⁻³), not
at m = 15 (1.75·10⁻² against 9.2–9.8·10⁻³; two of four seeds)".

**D3n-10 (cosmetic; `constrained_floor` docstring, `exp_adaptive_floor.py` l. 463–465 "at weights outside W the infimum of the
unregularised weighted loss"; the md's Part E l. 198 "min_{\|a\| ≤ a_max}⟨(f − f̂_w)²⟩_ρ"; PARADIGM l. 2998–3000 "Taking the constrained
minimum of the loss as the limit outside W, … at least 7.3·10⁻³, 2.2·10⁻⁴, 6.6·10⁻⁷"; RESULTS l. 588 "(constrained minima outside W) ≥").**

The values are the branch's local minima (one L-BFGS-B run from the bumps), not the infimum. Wherever a permutation of the classes over
the bumps is representable — 31%, 87% and 100% of the outside mass, as the same paragraph says — the constrained infimum is 0.

Evidence (`check_prob.py`):
- My recomputation reproduces the branch values exactly: 7.254·10⁻³, 2.156·10⁻⁴, and 5.169·10⁻⁷ on the vectors of probability > 10⁻⁶ at
  m = 240.
- With starts at all six assignments of the bumps plus random starts:
  - at (3, 6, 6)/15 the minimum is 1.3·10⁻³³, against the branch's 4.872·10⁻³;
  - at (6, 8, 1)/15, which has no representable permutation, it is 1.649·10⁻², against 5.340·10⁻²;
  - at m = 15, 34.6% of the outside mass has an exact fit, and the expected floor of the global minima is 9.89·10⁻⁴ — 7.3 times below
    7.25·10⁻³, and 0.25 rather than 1.8 times the fixed prior's law.

The text's reading — the dynamics from the warm clouds follows the branch — supports the branch values as the limit. So "at least" is a
lower bound conditional on that, not a bound on the infimum.

*Fix:*
- Docstring and md: "the constrained minimum on the branch of ϑ^e(q) (L-BFGS-B from the bumps; the global infimum is 0 where a permutation
  is representable)".
- Text: "taking the branch's constrained minimum as the limit (the dynamics from the warm clouds follows the branch), … at least …".

**D3n-11 (cosmetic; D‴.24 *Numerics* l. 3003–3004 "The coupled distances follow their laws more closely (E D/law 0.68, 0.93, 0.95 fixed;
0.53, 1.17, 1.26 adaptive)" and RESULTS l. 590–591; l. 2997–2998 "a failure set as large as the adaptive prior's"; Theorem D‴.23's
*Numerics* and §0 row "of the same expected size at both m").**

(a) *The adaptive prior's coupled distance has a heavy tail* at the boundary of W, where its gauge diverges (`check_mc.py`,
`check_mc_tail.py`):
- among the repository's own 600 samples, the largest D_e/law at m = 240 is 55.7 (mean 1.259; 1.167 without that one sample);
- among my 200 samples (seeds 990000 + i), one — counts (61, 99, 80), 0.0042 inside W — has 40.3 (mean 1.32 ± 0.22; 1.12 without it);
- so at m ≥ 60 the adaptive values (1.17, 1.26) are not "closer to their law" than the function-level 1.11 at m = 240, and their standard
  errors are not reliable;
- only the fixed prior's coupled distance follows its law closely: 0.68, 0.93, 0.95 (mine 0.76, 0.96, 0.94).

(b) *"A failure set as large as the adaptive prior's".* The fixed prior's saddles have probability 0.369, 0.260 and 0.109, against the
adaptive prior's outside-or-boundary mass 0.479, 0.127 and 0.0032. That is smaller at m = 15, and 2.0 and 34 times larger at m = 60 and 240.

(c) *D‴.23's Part F.* The draw means are higher at m = 15 in all four columns: by +0.86 to +1.66 combined standard errors, and per draw in
4/5, 4/5, 3/5 and 3/5 (`check_partF.py`). The columns share draws, as the md says.

*Fix:*
- "the fixed prior's coupled distance follows its law more closely; the adaptive prior's has a heavy tail at the boundary of W (one sample
  in 600 at 55.7 times the law at m = 240)";
- "a failure set comparable at m = 15 and much larger at m ≥ 60";
- "consistent with the same expected size at both m (five draws)".

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration, reproducibility, tests.**
   - `--from-raw` rewrites both `adaptive_prior.md` (md5 `44272ee9…`) and `adaptive_floor.md` (`2b15070d…`) byte-identically.
   - Part H re-runs bitwise.
   - Parts G, E and two Part F jobs of `exp_adaptive_prior.py` re-run bitwise (`rerun_EFG.py`).
   - 16 tests pass in 5.02 s; 129 are collected, as README says.
2. **D‴.27(i) mechanism** (twelve stress runs; see the verdict table).
   - The comparison inequality holds with the code's drift at every step.
   - The bound holds with the run's own force.
   - The leak targets and the prior stay in 𝒵, also at the split rows of m = 16.
   - An empty cell's atom stays frozen and the bound holds there.
3. **D‴.27(ii).**
   - 𝒥 is strictly decreasing along the projected midpoint rule (12001 steps) and at every refresh at m = 15 (η = σ₀, ∞).
   - The KKT condition and the balance hold at the end state.
   - The identity plan is optimal along the reference with gap ≥ 77.4.
   - LSODA at η = ∞ equals event-located DOP853 to ≤ 1.0·10⁻⁹.
4. **D‴.27(iv).**
   - The identification is exact.
   - Both η end at the same state (ϑ_a = 6.7277665395, z_a = 5.8888779583 on the face, fit error 0.141545).
   - The box stops the four-bump run-off at κ = 6.
5. **Round 18's per-sample checks, with an independent model** (`mymodel.py`: derivatives within 1.8·10⁻¹⁰, 7.1·10⁻¹⁰ and 1.75·10⁻⁹ of
   finite differences; Θ*(w) within 1.0·10⁻¹² of the repository's).
   - All 16 D‴.5 margins: −3.483, −1.085, +3.952, −6.142·10⁻³ (m = 15); +5.607, −2.156, +4.172, −3.977·10⁻³ (m = 60); +3.785, +5.059,
     +5.953, +5.132·10⁻³ (m = 240); +4.399, +5.441, +6.619, +4.516·10⁻³ (m = 510).
   - Exact-bias ratios: 1.0000 at the seven stable samples; 0.8718, 0.9729, 0.8880, 1.0103 and 0.8111 at the saddles.
   - Delta-method ratios: 0.86–1.72 (stable), 0.42–1.31 (saddles).
   - D: 0.944–1.000 of Θ*(w)'s (fixed) and 0.951–1.000 of ϑ^e(w)'s (adaptive); adaptive/fixed 1.7–4.8.
   - W² cross-class shares: 0.82–0.99 (fixed) and 0.67–0.99 (adaptive in W; 0.52 at the boundary sample); the priors' W² within
     1.00–1.32.
6. **The laws.**
   - (B) is 1.355928·10⁻², 6.166158·10⁻², 0.2291306 at 0.05, 0.15 and 0.5σ_crit, from finite differences of my own Θ*(w) and as the mode
     sum.
   - The coupled-distance laws are 21.32109 and 46.32494 (ratio 2.1727).
   - Law/bound is 0.7526, and the bound is exactly (σ₀/4)·46.32494.
   - Shares: modes 2–6 carry 0.994, 0.990 and 0.979 (0.972 unweighted); the band [1/5, 5] 0.422, 0.990 and 0.902 (0.213, 0.972, 0.784
     unweighted).
   - The law's values at m = 15, 60, 240 are 3.990·10⁻³, 9.068·10⁻⁴ and 1.360·10⁻⁴.
7. **Part E, exactly.**
   - m = 15: P(outside W) 0.3606, P(boundary) 0.1179, Hoeffding 1.870, permutation share 0.308, P(saddle) 0.3688.
   - m = 60: 0.0738, 0.0534, 0.602, 0.871, 0.2597.
   - m = 240: 0.002050, 0.001157 (over all vectors), 0.036, 1.000, 0.1089.
   - The branch minima are identical to the repository's.
8. **Monte Carlo with my own seeds** (200 nested samples; `check_mc.py`).
   - E‖f̂_m − f̂_M‖²/law: 2.82 ± 0.33, 1.51 ± 0.18, 1.11 ± 0.11 (repository 2.78, 1.53, 1.11).
   - The same, D‴.5-admissible pairs: 3.27 ± 0.47, 1.56 ± 0.24, 1.09 ± 0.12.
   - E D_f/law: 0.76, 0.96, 0.94. E D_e/law: 0.58, 1.12, 1.32.
9. **D‴.24(iv) and test 9.**
   - The boundary constant is confirmed: s·‖f̂ − f‖² = 0.4201, 0.4062, 0.4017 (η = ∞) and 0.8689, 0.8220, 0.8054 (η = σ₀) at s = 10⁴, 10⁵,
     10⁶, against 0.4 and 0.8; α₃(10⁶) = 13.487 = (a_max/4)ln s − 0.33.
   - The run-off law (a_max − a)s → a_max²(σ₀ + η)/(4ηG_αF_∞) is re-derived from the conservation law with α − z_α slaved.
   - It is confirmed at η = σ₀/10, σ₀ and ∞ (F_∞ = 0.033257, face fit error 0.104951).
10. **Wording fixes checked in place.**
    - D3k-6, D3k-7, D3k-8, D3k-10: present.
    - D3l-2, D3l-6: present.
    - The settled-gate sentences of D3l-5: present.
    - The D3l-1 correction to round 17's D3j-1 row: present.

## Not checked

- **Long runs.** Not rerun in full: `exp_adaptive_floor.py`'s Parts B–D and `--only P`, and `exp_adaptive_prior.py`'s Parts C, D and 38 of
  the 40 Part F jobs. The stored raw records were used, and single systems were re-integrated.
- **Material outside the brief.** Corollary D‴.26 and the planar check (`exp_adaptive_floor_2d.py`), except as quoted in the rows.
- **Theory not derived.**
  - The Kurdyka–Łojasiewicz argument with the indicator (D‴.27(ii) at finite η) was not proved.
  - The existence of solutions of the projected flow with plan switching was assumed, as in D‴.21.
  - Theorem D‴.23(v), Propositions D‴.5, D‴.9, D‴.12 and D‴.17 were used as given.
- **Literature.**
  - The Ioannou–Sun PDF (`flyingv.ucsd.edu`) was refused by the proxy (403), and arxiv.org and gnosis.library.ucy.ac.cy failed DNS.
  - The switching-σ statement in D3n-4 rests on the abstract text returned by a web search for a 1986 IEEE CDC paper of Ioannou's
    ([record](https://gnosis.library.ucy.ac.cy/handle/7/56943?show=full)) and on a survey's figure title "Switching σ-modification
    (continuous)" ([figure](https://www.oaepublish.com/articles/ir.2021.19/image/fig2)); it was not verified against the full texts.
- **Concurrent work.** The concurrent session's untracked files were not reviewed.
