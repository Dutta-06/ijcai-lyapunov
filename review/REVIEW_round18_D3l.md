# Referee report D3l (round 18, second referee) — Theorem D‴.23 as revised after round 17: PARADIGM §5d (statement (i)–(vi), Proof, Status, *Numerics*); the §0 and §7 rows "| D‴.23 |"; the §8 bullets "Done (Theorem D‴.23(ii))" (branch 7) and "Done (Theorem D‴.23)" (branch 8); `results/RESULTS_PARADIGM.md` (Theorem D‴.23); LITERATURE [S49] and the D‴.23 audit bullet; README's D‴.23 sentence; `src/core/adaptive_prior.py` (module docstring); `experiments/exp_adaptive_prior.py` → `results/chaos/adaptive_prior.md` (+ `adaptive_prior_raw.npz`); `tests/test_adaptive_prior.py`

**Setup.** Branch `claude/festive-bardeen-eopej2` at d2092a5 (HEAD throughout).

- At the start the working tree held five untracked files of the concurrent session: `exp_adaptive_floor_2d.py`, `exp_adaptive_floquet.py`, `src/core/adaptive_floquet.py`, `src/core/planar_loop.py` and `tests/test_planar_loop.py`. No tracked file was modified.
- During the review (18:52–19:12 UTC) that session left uncommitted edits to the five D‴.23 files: `adaptive_prior.py`, `exp_adaptive_prior.py`, the md, the raw file and the tests (186 insertions, 14 deletions), plus new untracked `adaptive_floquet` files. None of this was used or reviewed: the review is of d2092a5 as committed.
- Every script imports the d2092a5 tree, extracted with `git archive` into `ref_D3l/tree/`.
- No repository file was modified, and no git command touching the index or the tree was run.
- Python ran with `PYTHONDONTWRITEBYTECODE=1` and one BLAS thread; pytest ran with `-p no:cacheprovider` on the scratch tree.
- `--from-raw`, run on a second scratch copy (`ref_D3l/regen/tree/`), rewrites `adaptive_prior.md` byte-identically (md5 `733ead3c…`, the committed file's).

Line numbers are d2092a5's:
- PARADIGM: §0 row l. 119; Theorem D‴.23 l. 2631–2797 ((i) l. 2639, (v) l. 2683, Proof l. 2714, *Status* l. 2734, *Numerics* l. 2754–2797); §7 row l. 5367; §8 bullets l. 5564 and l. 5602.
- RESULTS l. 516–543; LITERATURE [S49] l. 72 and the audit bullet l. 351–366; README l. 50–52.
- `adaptive_prior.md`: Part A l. 5–32, Part D's opening state l. 143, Part E l. 145–163, Part F l. 165–174, Part G l. 176–218.

Scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref_D3l/scripts/`. Each `.out` file is one level up, in `ref_D3l/`.

| script | what it does |
|---|---|
| `common.py` | shared setup (imports the archived tree; the K₀-atom loss, gradient and drift) |
| `inspect_raw.py` | structure of the raw file |
| `recompute.py` | every number of the paragraph, the rows and RESULTS, recomputed from the raw records independently of `report()` (λ_1^− by `np.roots`, the pair solutions by `expm`); per-seed gate comparison |
| `g_a_endstate.py` | Part G(a): the gradient and the physical partials at the "refined" atoms; the minimiser of the loss on the face s = s_max; the unregularised flow from the refined atoms to s = 10⁸ |
| `g_a_flow.py` | Part G(a): the joint flow from the exact split fit re-run (bitwise against the raw), its state at s = 4·10⁴, the collapsed flow continued to 10¹⁰ |
| `g_a_saddle.py` | Part G(a): the saddle; the plan's optimality gap over all 90 assignments; the joint Jacobian at h = 10⁻⁴…10⁻⁷; the frozen-plan Hessian of 𝒥; the unstable eigenvector; the flow from the saddle perturbed along it |
| `g_b_asym.py` | Part G(b): the asymptotics against a = a_max tanh(α/a_max); the conservation law along the solution |
| `cn_scale.py` | (v): the sampled C_N reproduced, its scale dependence, the local Lipschitz constant at the exact representation |
| `part_f.py`, `part_f_ratio.py` | Part F re-run (five systems, bitwise): the plan at every refresh, the starts' norms, the separation decay, the dependence on the separation draw at m = 15 and 60; the second-order factor per draw |
| `part_e_dense.py` | Part E: the crossing times by event location |
| `check_tests.py` | the asserted quantities of the two new tests and the renamed one; mutations |
| `misc.py` | (i)'s relabelling with unequal masses; (ii) on a 1600-point grid; the norm of Part D's starts |
| `schedule.py` | 𝒥 under a random schedule η(s) ∈ {0, σ₀/10, σ₀, 10σ₀, ∞} |
| `open_state.py` | the settled gate's opening state, from the md's formula and directly |

Also in `ref_D3l/`: `tests.out` (9 passed, 2.91 s) and `from_raw.out`. `dl_rock/` is empty: the download of Rockafellar (1976) was refused by the proxy (403), and arxiv.org and the springeropen hosts failed DNS. The literature points below therefore rest on memory, and are marked so.

## Verdicts

| object | verdict |
|---|---|
| **(i)** | **holds for replicated configurations and for bounded solutions, as now stated.** Checked: the stationary equations G⟨r∇Φ̃(ϑ_i)⟩ = σ₀(ϑ_i − T̂_π(i)) with z = T̂; the replicated characterisation; the conservation law (residual ≤ 1.7·10⁻¹⁷ along an unbounded solution) and its proof; "𝒥 non-increasing under any schedule η(s) ≥ 0" (660 records under a random schedule mixing 0, σ₀/10, σ₀, 10σ₀ and ∞, no increase). Boundedness is stated as a hypothesis where it is used: (i), *Status*, §0. (ii)–(vi) do not need it, being local, formal or the monotonicity of 𝒥. Two defects remain. The relabelling parenthesis "a permuted prior at an exact representation is stationary" is false for unequal masses (D3l-2). The numerical illustration of the K₀-cluster bias ends at a point at infinity, presented as a critical point: it is a second unbounded solution, not recognised as one (D3l-1). |
| **(ii)** | **holds.** On 1600 (h, η) points, h ∈ [10⁻⁶, 10] and η ∈ [10⁻⁶, 10³], the strict bounds and dλ⁻/dh = (η − λ⁻)/(λ⁺ − λ⁻) > 0 hold with no violation, the derivative to 10⁻⁵ against central differences. "Tight as η → 0": ratio 0.9895 at η = σ₀ and 0.9996 at σ₀/100 for h_1. The limits h_jη = 0 and η = ∞ are right, and so is the proof (p(h) = −σ₀h, p(η) = −σ₀η, implicit differentiation). |
| **(iii)** | **holds to first order** (unchanged). The per-mode table and the shift table recomputed. |
| **(iv)** | **holds as stated (formal).** The new sentence "a target that the network fits exactly only with a class split is fitted with K₀ clusters" holds. Its example's end state is misdescribed (D3l-1). |
| **(v)** | **holds as a sketch; the bookkeeping is now right.** In |||·|||_η: ‖ζ_k‖ ≤ ρ(η/(σ₀q_k))^{1/2}; per-particle extent (1 + q_min^{−1/2})ρ. max_i‖u_i‖ + max_k‖ζ_k‖ < Δ*/2 suffices for the diagonal plan (triangle inequality: the shifts of both atoms enter with weight ½), so the bracket [1 + q_min^{−1/2} + (η/(σ₀q_min))^{1/2}] is a correct sufficient condition. The η = ∞ term margin/(C(1 + 2q_min^{−1/2})) is correct but not sharp: (1 + q_min^{−1/2}) suffices, since a particle's distance to its own atom is then ‖u_i^⊥‖. The degeneration like η^{−1/2} and the vanishing as η → 0 are right. C_N is η-independent. The sampled bound C_N ≥ 0.3406 holds at the proof's own scale: the quotient is the same from radius 10⁻⁵ to 5·10⁻², and the local Lipschitz constant at the exact representation is ≥ 0.394. So the quoted radius bounds are valid upper bounds, and slightly loose. "By these factors" undercounts Parts D–F (D3l-4). |
| **(vi)** | **holds in substance** (unchanged apart from D3j-9 and D3j-11). The comparison of the two gates holds in seed means only (D3l-5). |
| Proof | **correct where it is a proof**: (i), (ii), (iii), (iv) formal, (v) explicitly a sketch, (vi). |
| *Status* | the labels are accurate. The per-mode credits are right in substance. Two problems: "σ₀/(h_j + σ₀) is the Tikhonov filter factor" inverts the standard term, and the "new … limits of the bias removal" are, at least in part, the cited precedents' own properties (D3l-6). The test description says "unbounded solution" for a finite-horizon check (D3l-8). |
| *Numerics* paragraph | **every number reproduces** (confirmation 2). Readings that overreach: Part G(a)'s "critical point on the saturation boundary" (D3l-1); "the test of (v)" and "of the same size at both m" from one separation draw per m (D3l-3); "the runs below … 94 to 6·10⁴ times outside" (D3l-4); "the settled gate's ‖ȳ‖ and W² are the lower everywhere" (D3l-5). The Part E crossing times are grid upper bounds (D3l-7). |
| `src/core/adaptive_prior.py` (docstring) | **holds** (the D3j-1 and D3j-2 wording is in). The `plan='nearest'` path (Proposition D‴.24) was not reviewed; the OT path is unchanged in substance (`git diff 1aeaa56 d2092a5`). |
| `experiments/exp_adaptive_prior.py` | **holds as code**: regeneration byte-identical; Part F, the Part G(a) flows and a settled Part D system re-run bitwise. The Part G(a) "refinement by least squares" stops on a saturated slope, not at a zero (D3l-1). |
| `results/chaos/adaptive_prior.md` | **numbers hold**; wording as in D3l-1, D3l-3, D3l-4 and D3l-7. |
| `tests/test_adaptive_prior.py` | **9 passed (2.91 s); none vacuous.** Test 8 fails without the leak (σ₀ → 10⁻⁶σ₀: leak force 1.2·10⁻⁸, max Re 3.2·10⁻⁸). The renamed test 7 checks what its name says (drift 9.95·10⁻⁴ at the opening, 1.002·10⁻³ at the previous refresh, no earlier sub-tolerance record, transit over at s = 14). Test 9 is a finite-horizon proxy calibrated to η = σ₀ (D3l-8). |
| §0 row, §7 row, §8 bullets, RESULTS, LITERATURE, README | **numbers hold.** Wording: D3l-1 (RESULTS "fit error 1.94·10⁻² on the saturation boundary"; the round-17 log row), D3l-3 (§0, §7, RESULTS), D3l-4 (RESULTS), D3l-5 (README "provided"), D3l-6 (audit bullet, [S49]). |

### The round-17 fixes, item by item

| id | fixed as logged? | new problem introduced? |
|---|---|---|
| D3j-1 | yes: heading "replicated", the stationary equations, the plan's labelling, the non-replicated points with Part G(a), the tautology replaced; §0, §8 and the docstring | the parenthesis dropped the logged "(equal masses)" (D3l-2); Part G(a)'s end state is misdescribed and is itself an unbounded solution (D3l-1) |
| D3j-2 | yes: the hypothesis, the conservation law with its proof, Part G(b), test 9, §8 branch 7 | test 9 checks a finite horizon (D3l-8) |
| D3j-3 | yes: ρ(η) stated; the η = ∞ norm and the two limits; *Status* "sketched"; C_N ≥ 0.3406; the factors | the factors are Part C's start only (D3l-4) |
| D3j-4 | yes: the opening state printed (and right, confirmation 6); Part F added | "the test of (v)" from starts 133–3344 times outside its ball; one draw per m (D3l-3) |
| D3j-5 | yes: the quadratic analyses credited; Hanke–Groetsch added | "filter factor" inverted; the "limits" are claimed as new (D3l-6) |
| D3j-6 | yes | none |
| D3j-7 | yes (13 of 36 modes; 0.981–1.003; 1.04% / 5.96%; "(η = σ₀)" in §7) | none |
| D3j-8 | yes | none |
| D3j-9 | yes | none |
| D3j-10 | yes (renamed; it tests what the name says) | none |
| D3j-11 | yes: "should", the wait's cost, README "after the transit and the fixed-prior relaxation" | the comparison sentence holds in seed means only; README keeps "provided" (D3l-5) |
| D3j-12 | yes (LSODA flow, crossing times) | the times are grid upper bounds, by 0.3–1.8% (D3l-7) |

## Confirmed issues

**D3l-1 (minor; md Part G(a) l. 180–191 — the table's "refined" column, the two "refined class atoms" lines and "a critical point of the unregularised weighted K₀-atom loss on the boundary, flat in that direction"; `boundary_job` docstring l. 287–296, "refined at the end by least squares"; PARADIGM *Numerics* l. 2793–2794, "the merged atom then widens toward s_max, approaching a critical point of the unregularised K₀-atom loss on the saturation boundary (fit error 1.94·10⁻²)"; RESULTS l. 538, "fit error 1.94·10⁻² on the saturation boundary"; REVIEW_LOG round 17, D3j-1 row, "approaches a boundary critical point with the merged width at s_max").**

The end state of the joint flow from the exact split fit is a point at infinity of the unconstrained coordinates, and the solution is unbounded. (i) is stated in those coordinates, with its "bounded solutions" hypothesis and its Clarke-critical limits.

*The run itself does not reach the boundary.* Re-run bitwise (`g_a_flow.py`; fit, split and 𝒥 at all eleven records equal to the raw to 0.0):
- at s = 4·10⁴ the classes are replicated to 2.2·10⁻¹⁶ (η = ∞) and 4.4·10⁻¹⁶ (η = σ₀);
- the merged atom's width coordinate is κ = −0.245 (η = ∞) and −0.610 (η = σ₀), i.e. widths 1.05 and 0.90 against s_max = 2;
- the fit error is 2.227·10⁻² and 2.385·10⁻², not 1.94·10⁻²;
- the unregularised drift at the class atoms is still 1.4·10⁻⁵ and 2.7·10⁻⁵.

*Continued, the width coordinate runs off.* The collapsed dynamics (the unregularised K₀-atom flow, LSODA rtol 10⁻¹⁰) gives, at +10⁶, 10⁷, 10⁸, 10⁹, 10¹⁰ slow units:
- κ = 2.16, 3.84, 5.53, 7.53, 9.75;
- s_max − s = 0.176, 3.6·10⁻², 6.7·10⁻³, 9.1·10⁻⁴, 9.9·10⁻⁵;
- κ's increments per decade over ln 10: 0.70, 0.73, 0.73, 0.87, 0.97, tending to 1.

This is κ ~ ln s and s_max − s ~ 1/s: Part G(b)'s mechanism, with the sigmoid in place of tanh. The fit error tends to 1.943285·10⁻². That is the minimum of the loss on the face s₂ = s_max (BFGS on the other eight parameters: 1.943283·10⁻², gradient 2.7·10⁻¹²; `g_a_endstate.py`).

*The face minimiser is not a critical point, and the loss is not flat there.*
- The physical slope is ∂L/∂s₂ = −1.326·10⁻⁶ < 0 at the face minimiser (−1.336·10⁻⁶ at the "refined" atoms): the loss still decreases outward. It is a KKT point of the box-constrained problem.
- In the unconstrained coordinates the gradient is −3.7·10⁻¹⁰ in κ at κ = 8.72. It decays like e^{−κ} and never vanishes, so no critical point exists there.
- "Flat in that direction" (Hessian eigenvalue 2.8·10⁻¹⁰) is the saturation map's flatness.
- The least-squares "refinement" moves κ from −0.25 to 8.72 (η = ∞) or 8.86 (η = σ₀). It stops where the residual happens to fall to 5.6·10⁻¹⁰ and 4.8·10⁻¹⁰. The two "refined" points differ, as they must on a slope with no zero.
- From the η = ∞ "refined" atoms the flow moves κ by 0.027 in 10⁸ slow units: the run-off is slow, but real.

The example is therefore a second instance of the lost a-priori bound, and the more telling one. The target admits an exact m-particle representation, which is the starting point. Yet the adaptive prior's K₀-cluster bias drives the solution off to infinity, while 𝒥 decreases.

*Fix:*
- md, paragraph and RESULTS: "the merged atom then widens toward s_max — its width coordinate runs off like ln s and the width tends to s_max like 1/s, as in (b): the solution is unbounded — and its physical image approaches the minimiser of the loss on the face s = s_max (a KKT point, ∂L/∂s = −1.3·10⁻⁶, fit error 1.94·10⁻²), which is reached only as s → ∞ (width 1.05, fit error 2.23·10⁻² at s = 4·10⁴)".
- Replace the least-squares "refinement" by the face minimiser, or by the long flow. Drop "critical point … flat in that direction".
- In (i), cite Part G(a) next to (b) for "solutions can leave every compact set".

**D3l-2 (cosmetic; (i) l. 2649–2651, "labelled by the atoms of a given configuration this reads z_{σ(k)} = ϑ_k for the plan's class-to-atom bijection σ, and a permuted prior at an exact representation is stationary").** The theorem is stated for general masses q. A replicated configuration with W² = 0 needs q_σ(k) = q_k, so the statement holds only for permutations among atoms of equal mass. The logged fix (D3j-1) said "(equal masses)"; the text dropped it.

Counterexample (`misc.py`): q = (1/2, 1/3, 1/6), a target that the instance's atoms represent exactly with these weights (fit error 1.1·10⁻¹⁶), m = 6. For each of the three non-trivial permutations of the prior's atoms:

| permutation | W² | drift (Θ) | drift (Z) |
|---|---|---|---|
| cyclic | 16.10 | 4.40·10⁻² | 1.47·10⁻² |
| (0 1) | 8.65 | 4.40·10⁻² | 1.47·10⁻² |
| (1 2) | 6.54 | 3.61·10⁻² | 1.81·10⁻² |

None is stationary. *Fix:* "a prior permuted among atoms of equal mass".

**D3l-3 (cosmetic; md Part F l. 165–174; `sep_job` docstring l. 247; PARADIGM *Numerics* l. 2787–2790, "*The test of (v)* … of the same size at both m, as (v) and Theorem D‴.18's facts (1)–(2) predict"; §0 row and §7 row, "alike at m = 15 and 510"; RESULTS l. 535–536).** Part F is a good experiment, but it is not a test of (v), and its m-comparison rests on one separation draw per m.

(a) *The starts lie far outside (v)'s ball.* |||start|||_η (`part_f.py`):
- η = σ₀: 0.371 (the prior's offset), 3344 times the bound;
- η = ∞: 0.030 (m = 15, rms 0.02), 0.149 (m = 15, rms 0.1) and 0.246 (m = 510, rms 0.1), 133–1095 times the bound — not 94.

(v) asserts nothing out there. The second-order law is fact (1)'s, not (v)'s. Part A's own sentence, "they probe the dynamics beyond it", is the right reading.

(b) *At m = 15 the departure depends strongly on the draw.* Five other separation draws at η = ∞, rms 0.1:

| m | departures, draws 0–4 | mean | sd |
|---|---|---|---|
| 15 | +0.104, +0.278, +0.189, +0.050, +0.193 | 0.163 | 0.079 |
| 60 | +0.140, +0.127, +0.152, +0.222, +0.212 | 0.171 | 0.039 |

At m = 15 the departure varies by a factor 5.6 across draws, more than the factor 1.32 between the md's m = 15 and m = 510 values (0.109 and 0.145, from one draw). The m-independence holds in expectation: the means agree, and the spread shrinks with m, as the class covariances concentrate. The md's own explanation, "(the separations' class covariances, not their number, enter at second order)", says exactly this. The second-order law holds draw by draw: factors 24.17–24.91 (m = 15) and 24.13–24.97 (m = 60) for the factor 5 (`part_f_ratio.py`).

(c) *Other checks pass.*
- The five re-run systems are bitwise equal to the raw.
- The plan is the class partition at all 1751 refreshes.
- The class means are zero to 7·10⁻¹⁸, and the rms is in the G⁻¹ norm.
- The separations are gone by s = 1750 (max 10⁻⁶).
- The K₀-atom reference is the right comparison: a zero-separation m-system is the K₀-atom system by invariance, so the departure is precisely the separations' second-order kick, kept by the slow modes.

*Fix:*
- "A test of the structure behind (v) — facts (1)–(2) — beyond its ball".
- "Of the same size at both m in expectation", with several draws per m: their means and spread.

**D3l-4 (cosmetic; md Part A l. 32, "The runs of Parts C–F start outside the ball that the proof of (v) covers, by these factors"; PARADIGM l. 2760–2762, "the runs below start at |||·|||_η = 1.17, 0.37, 0.12, 0.021, 94 to 6·10⁴ times outside"; RESULTS l. 523, "while the runs start 94–6·10⁴ times outside").** The table's factors are those of Part C's start (Θ*(q), the shifted prior). Parts D and F start farther out.

- *Part D's warm clouds* (`misc.py`, seeds 0 and 3, m = 15 and 510) have largest separations 5.3–7.1 at s = 0, so |||start|||_η ≥ 5.3. That is 4.8·10⁴–6.4·10⁴ times the bound at η = σ₀ (the table says 3344), and 2.4·10⁴–3.2·10⁴ times at η = ∞ (the table says 94).
- *Part F* at η = ∞ starts 133–1095 times outside.

*Fix:* "by at least these factors (Part C's start; the warm clouds of Part D start 2.4·10⁴–6.4·10⁴ times outside, Part F 133–3344 times)".

**D3l-5 (cosmetic; PARADIGM *Numerics* l. 2785–2787, "the settled gate's ‖ȳ‖ and W² are the lower everywhere, and its fit error at η = ∞ and at m ≤ 240"; README l. 51–52, "provided the adaptation starts after the transit and the fixed-prior relaxation").**

The comparison holds in the seed means of the md's table, not system by system. At s = 3000 the transit gate ends lower than the settled gate in all three of ‖ȳ‖, W² and the fit error in 6 of the 32 systems (`recompute.py`):

| η | m | seed | ‖ȳ‖ settled / transit | fit error settled / transit |
|---|---|---|---|---|
| σ₀ | 60 | 0 | 2.58·10⁻³ / 1.23·10⁻³ | 1.58·10⁻⁴ / 6.6·10⁻⁵ |
| σ₀ | 240 | 2 | 2.59·10⁻³ / 1.90·10⁻³ | 1.59·10⁻⁴ / 1.03·10⁻⁴ |
| σ₀ | 240 | 3 | 2.72·10⁻³ / 1.53·10⁻³ | 1.68·10⁻⁴ / 4.6·10⁻⁵ |
| σ₀ | 510 | 2 | 2.67·10⁻³ / 1.55·10⁻³ | 1.65·10⁻⁴ / 6.1·10⁻⁵ |
| σ₀ | 510 | 3 | 2.59·10⁻³ / 2.07·10⁻³ | 1.58·10⁻⁴ / 5.9·10⁻⁵ |
| ∞ | 60 | 0 | 9.7·10⁻⁴ / 4.6·10⁻⁴ | 3.35·10⁻⁵ / 1.27·10⁻⁵ |

The README's "provided" also makes the gate a condition for the bias removal, which (vi) no longer says. In (vi) the gate is what one should use "when the parameters, not only the fit, are wanted at a given horizon". Ungated, the bias is still removed, more slowly.

*Fix:* "in the seed means; per system the transit gate ends lower in all three in 6 of the 32". README: "and is best started after the transit and the fixed-prior relaxation".

**D3l-6 (cosmetic; *Status* l. 2741 "σ₀/(h_j + σ₀) is the Tikhonov filter factor"; audit bullet l. 355 "the Tikhonov filter factor σ₀/(h_j + σ₀)"; [S49] l. 72 "contracts the error … by the filter factor σ/(h + σ)"; *Status* l. 2746 and audit bullet l. 363 "new are … the limits of the bias removal (the K₀-cluster bias, the lost a-priori bound)").**

(a) *Terminology.* In the standard terminology (Engl–Hanke–Neubauer; Hansen) the Tikhonov filter factor is h/(h + σ), the fraction of a component that the regularised solution keeps. σ/(h + σ) is its complement: the fraction of the prior's error that remains, i.e. the per-step error contraction of iterated Tikhonov. The round-17 report wrote "residual factor"; the fix renamed it to the inverse term.

(b) *The two "limits" are largely the precedents' own properties.*
- The lost a-priori bound is the particle-system form of the proximal point method's behaviour without a minimiser. For a maximal monotone operator the proximal sequence is bounded iff a zero exists — Rockafellar 1976, as it is usually cited. The paper could not be fetched (proxy 403), so this attribution rests on memory.
- The K₀-cluster bias is the quantisation error of soft weight-sharing's hard (K-means) limit, which the text itself cites.

The §0 row's novelty column omits both, so *Status* and the audit bullet are also inconsistent with it.

*Fix:* "the complement 1 − h_j/(h_j + σ₀) of the Tikhonov filter factor (the residual factor)" in all three places. Either drop the limits from "new", or say what about them is specific: the particle and mass-constrained form, and the run-off at a fit that is exact for the network (D3l-1).

The rest checks:
- the Moreau-envelope sentence: its curvature hσ₀/(h + σ₀), and the slow-prior rate ηh/(h + σ₀) with mobility η/σ₀;
- "empirical Bayes … is the same proximal iteration" (EM for the prior mean is iterated Tikhonov);
- EASGD's quadratic worker–centre analysis.

These rest on memory, as the file says.

**D3l-7 (cosmetic; md Part E l. 163 "the loss stays above 0.9 of its 2.5220e-04 until 4.18e+04, halves by 4.87e+04, falls below 10⁻¹² by 8.00e+04"; PARADIGM l. 2776–2779 "for 4.2·10⁴ more slow units, halves it by 4.9·10⁴"; §7 row and RESULTS "8.0·10⁴").** The times are the first grid points below each level, on a logarithmic grid of ratio 1.039. Event location from the same atoms gives (`part_e_dense.py`):

| level | event time | md (first grid time below) | md / event |
|---|---|---|---|
| 0.9 L₀ | 4.156·10⁴ | 4.178·10⁴ | 1.005 |
| 0.5 L₀ | 4.780·10⁴ | 4.867·10⁴ | 1.018 |
| 10⁻¹² | 7.971·10⁴ | 7.997·10⁴ | 1.003 |

"By" is correct as an upper bound. "Stays above … until 4.18·10⁴" is not: the loss is below 0.9 L₀ from 4.156·10⁴ on. The paragraph's "4.9·10⁴" rounds the upper bound; the crossing rounds to 4.8·10⁴.

*Fix:* locate the crossings by events, or write "between 4.02·10⁴ and 4.18·10⁴".

**D3l-8 (cosmetic; `tests/test_adaptive_prior.py::test_the_adaptive_prior_gives_up_the_leaks_a_priori_bound`, l. 216; the module docstring "the adaptive prior's amplitude coordinate runs off"; *Status* "the adaptive prior's unbounded solution").** The test asserts finite-time behaviour calibrated to η = σ₀:
- an increment of ϑ_a over [10³, 10⁴] above 3.0 (measured 3.80);
- a drift above 10⁻⁵ (measured 1.9·10⁻⁴);
- a(10⁴) > 3.97 (measured 3.978).

An adaptive prior with η = σ₀/10, whose solution runs off by the same asymptotics, fails two of the three (increment 2.37, a = 3.897). A bounded solution still converging at s = 10⁴ could pass. The test does not test unboundedness; it is a regression check of the run.

*Fix:* assert the mechanism, which a bounded solution cannot satisfy:
- the amplitude force ⟨rφ₀⟩ stays above a positive bound while a → a_max (0.0333 here);
- the conservation law holds;
- or (a_max − a)·s approaches 8/⟨rφ₀⟩ (η = σ₀).

Or rename the test "… the adaptive prior's amplitude coordinate is still growing at s = 10⁴, the fixed prior's converged".

## Independent confirmations (what was recomputed and agrees)

1. **Regeneration and tests.**
   - `--from-raw` rewrites `adaptive_prior.md` byte-identically (md5 `733ead3c…`).
   - 9 tests pass in 2.91 s on the archived tree.
2. **Every number of the paragraph, the §0/§7 rows and RESULTS**, recomputed from the raw records independently of `report()` (`recompute.py`).
   - Part A: h = 3.0587·10⁻⁴ … 0.19521; σ₀/(h + σ₀) = 0.959 … 0.035; h + σ₀ ≥ 7.3803·10⁻³; λ_1^− = 2.6844·10⁻⁵, 1.5128·10⁻⁴, 2.7796·10⁻⁴, 3.0587·10⁻⁴ (`np.roots`).
   - (v)'s radius: C_N ≥ 0.3406; ρ ≤ 1.970·10⁻⁵, 1.110·10⁻⁴, 2.040·10⁻⁴, 2.245·10⁻⁴; starts 1.1723, 0.3713, 0.1191, 0.0210; factors 59507, 3344, 584, 94; λ₀/σ₀ = 0.0038–0.0432.
   - Part B: fit error 9.7765·10⁻³; 4.805%; modes 1 and 3 at 0.849 and 0.727; the weakest mode keeps 0.814 against 0.959; shift table 4.805, 2.341, 1.157, 0.575%, with modes 1 and 3 at 0.977 and 0.969 at 0.01875.
   - Part C: 13 log-slopes, 0.9812–1.0031; ratios to the `expm` pair solutions 0.8201–1.1146; changes from s = 1000 on of 5.96, 1.04, 0.02, 0.01%; fit errors at 6000 3.612·10⁻⁴, 1.825·10⁻⁵, 8.035·10⁻⁶, 6.738·10⁻⁶; the weakest mode 1.450·10⁻³, 6.663·10⁻⁴, 2.979·10⁻⁴, 2.498·10⁻⁴ of 1.582·10⁻³.
   - Part D: every seed-mean range quoted (0.032–0.32 / 0.048–0.71; 0.0094–0.10 / 0.012–0.59; 4.9·10⁻³–3.7·10⁻²; 7.4·10⁻³, 4.3·10⁻³, 1.3·10⁻³); t_open 88–149 and 1211–1304; per-seed ratios 0.9953–1.0000 and 0.9676–1.0005 (‖ȳ‖), 0.9984–1.0001 and 0.9791–1.0002 (fit); fit errors 1.594–1.649·10⁻⁴ and 3.369–3.544·10⁻⁵; W² 6.748–7.161·10⁻⁶ and 9.262–9.889·10⁻⁷; at m = 510 and η = σ₀, 1.547·10⁻⁴ (transit) against 1.638·10⁻⁴ (settled).
   - Part E: ‖ȳ‖ = 2.124; fit error 3.4639·10⁻² and drift 1.38·10⁻⁵ at s = 23000; L = 2.5220·10⁻⁴; Hessian −1.841·10⁻⁵, 7.159·10⁻⁵, 9.461·10⁻⁴; 1/|ev| = 5.43·10⁴.
   - Part F: +0.084, +0.071, +2.093, +1.772% (η = σ₀) and +0.442, +0.587, +10.930, +14.460% (η = ∞); growth factors 24.88, 25.03, 24.72, 24.62.
   - Part G: all values of (a), (b), (c).
3. **Reproducibility.** Re-runs with the committed code are bitwise equal to the raw:
   - the five Part F systems (‖ȳ‖ and fit error, difference 0.0);
   - both Part G(a) flows (all eleven records, difference 0.0);
   - the settled seed 3, m = 15 system (t_open = 1211).

   Part G(b) to 1.3·10⁻⁸, with LSODA integrated to 10⁷ instead of 10⁶.
4. **The saddle of Part G(a) is genuine** (`g_a_saddle.py`).
   - The plan is optimal with a gap of 9.85 over all 90 assignments (the next best costs 10.07 against 0.224), and every particle is nearest its own atom. The frozen plan is therefore the plan in a neighbourhood.
   - The largest eigenvalue of the joint Jacobian, +3.661099·10⁻⁵ (η = σ₀) and +6.635507·10⁻⁵ (10σ₀), is real and identical to seven digits for h = 10⁻⁴…10⁻⁷.
   - M·J is symmetric to ≤ 1.1·10⁻⁹. The frozen-plan Hessian of 𝒥 has the eigenvalue −3.0508·10⁻⁵ at both η, as it must, since 𝒥 does not involve η.
   - It is stationary for every η > 0 (Z = T̂ and the η = ∞ equations; drift 4.4·10⁻¹⁷).
   - The unstable direction is the common width of the split pair and its atom (|v₄ + v₅| = 2.24 against |v₄ − v₅| = 0.09), not the merging. Perturbed by ±10⁻³ along it, the pair's distance drifts to 1.646 and 1.632 by s = 1.2·10⁵. The text claims nothing about the direction.
5. **Part G(b) re-derived** (`g_b_asym.py`).
   - At η = ∞ (m = 3, leak zero), dα/ds = ⟨rφ₀⟩sech²(α/a_max) with ⟨rφ₀⟩ → 0.03326 > 0. Hence α = (a_max/2)ln s + O(1) = 2 ln s + O(1), and (a_max − a)s → 4/⟨rφ₀⟩.
   - At η = σ₀ the conservation law halves the speed, giving 8/⟨rφ₀⟩.
   - Measured at s = 10⁷: (a_max − a)s = 120.26 against 120.27, and 240.50 against 240.54. α − 2 ln s → −5.420 and −6.807; the difference 1.387 is ln 4, the halved speed.
   - The other coordinates converge (drift 2.2·10⁻¹² and 4.5·10⁻¹² at 10⁷), and the conservation-law residual is ≤ 1.7·10⁻¹⁷.
6. **The opening state** (`open_state.py`). √(W²_ex − ‖ȳ‖²) is the rms separation: 7.9200·10⁻⁴ by the md's formula and directly, at s = 1200, seed 3, m = 15. At the opening itself (s = 1211) the rms is 7.37·10⁻⁴, the largest separation 1.52·10⁻³, and the class-mean offset 6.12·10⁻⁵. Over the 32 runs, 1.24·10⁻⁴–7.92·10⁻⁴ and 7.38·10⁻⁷–6.66·10⁻⁵, as printed.
7. **(v)'s C_N** (`cn_scale.py`).
   - 0.340644 is reproduced with the same points.
   - Over 40 random directions the quotient ‖J(P) − J(0)‖₂/‖y‖ is 0.2768–0.2830 from radius 10⁻⁵ to 5·10⁻² (mean 0.1505 at every radius): the Jacobian is linear in y over this range.
   - The local Lipschitz constant at the exact representation, maximised over directions, is ≥ 0.3944 (0.3944 at radius 10⁻⁴, 0.4053 at 5·10⁻²). So C_N ≥ 0.3406 holds on the proof's own ball (radius ≤ 2.2·10⁻⁴).
   - The collective-block quotient is a valid lower bound for the |||·|||_η constant: at replicated points the variation has no separation component, and it does not act on ζ.
   - With 0.3944 the radius bounds would be 1.70·10⁻⁵, 9.59·10⁻⁵, 1.76·10⁻⁴, 1.94·10⁻⁴.
8. **(i): 𝒥 under any schedule** (`schedule.py`). In an m = 15 system with mass in transit, η was redrawn every 10 slow units from {0, σ₀/10, σ₀, 10σ₀, ∞} (11, 15, 6, 15 and 13 segments). 𝒥 falls from 9.27·10⁻³ to 6.72·10⁻⁷ with no increase in 660 records.
9. **(ii)** on 1600 grid points: no violation of either strict bound or of the sign and value of dλ⁻/dh (`misc.py`).
10. **The tests' quantities and mutations** (`check_tests.py`).
    - Test 8: leak force 5.80·10⁻³ and unregularised gradient 5.80·10⁻³ (tolerances 10⁻³); max Re 3.66·10⁻⁵ (tolerance 10⁻⁶); the class stays split at 1.640. With σ₀ → 10⁻⁶σ₀ its leak, gradient and saddle assertions fail (1.2·10⁻⁸, 1.2·10⁻⁸, 3.2·10⁻⁸), so it discriminates.
    - Test 7: t_transit = 15, t_open = 441; maximum Δ𝒥 −3.0·10⁻¹⁰.

## Not checked

- The 30-minute experiment was not rerun in full. Single systems of Parts D, F and G were reproduced bitwise, and the regeneration is byte-identical.
- The constants κ, C and the margin of (v) were not derived, only the formula's internal consistency and the sampled C_N. The separations' block of C_N was not sampled; it can only raise C_N and so lower the radius bounds.
- Proposition D‴.17's KL argument and Theorem D‴.18's facts (1)–(3) were taken as given.
- The literature: Rockafellar (1976), Hanke–Groetsch (1998) and EASGD could not be fetched (proxy 403 for sites.math.washington.edu; DNS failures for arxiv.org and springeropen). D3l-6 and the confirmation of the precedents rest on memory.
- Proposition D‴.24, Remark D‴.25, the `plan='nearest'` code path, and the concurrent session's untracked files and uncommitted edits to the five D‴.23 files were not reviewed.
- Incommensurate sizes, the real loop with an adaptive prior, and learned masses are outside the material.
