# Referee report D3ch, round 87: Proposition D‴.59 as revised after round 83

## Scope

- **Export.** Commit 7be2cb0 (branch claude/festive-bardeen-eopej2), taken with `git archive`. W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref87/D3ch_work; the export is in W/repo.
- **Files read in full:**
  - `review/drafts/D59.md` (331 lines)
  - `experiments/exp_moving_width.py` (1887 lines)
  - `results/chaos/moving_width.md` (613 lines)
  - `tests/test_moving_width.py` (46 tests)
  - `review/REVIEW_round83_D3cd.md`, and the Round 83 entry of `review/REVIEW_LOG.md`
  - `git diff 3b6aac7 7be2cb0` for the script, the md, D59.md, the tests and the D59 row of `review/drafts/README.md`; saved to W/diffs/*.diff
  - `src/core/moving_width.py`, and in `src/core/rich_transit_nd.py`: `loop_nd`, `critical_rhs_nd`, `critical_jac_nd` and `_curve_terms`
  - The raw, loaded with `exp_cohesion_leak.load`
- **Settings.**
  - Every Python process ran with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, using the repo's venv.
  - pytest ran as `-q -p no:cacheprovider`, with `--basetemp` under W, and `-x` for mutants.
  - `MW_OUT` always pointed to W/mwout, or to W/mut/<id>/out, each holding a copy of the raw.
  - I ran at most two heavy processes at a time. A few raw-only scripts (seconds each) ran beside them.
- **I modified nothing in the repository.**
  - In the repository I ran only `git archive`, `log`, `show`, `diff` and `status`.
  - At the end `git status` was clean. D59.md, the script, the md, the raw, the tests, both src modules and the drafts README are byte-identical to the export.
  - During the review another session committed 066dc8f (D‴.60 round 85 report). It touches none of these files.
  - Everything I wrote is under W. Mutant trees in W/mut symlink to the export, never to the repository. Nothing of mine is still running.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my run or my computation; [snippet] a web-search result; [memory] not verified this round.

## Summary

**What checks out**

- **Reproducibility** [own]:
  - `--from-raw` regenerates the md byte for byte, and the 46 tests pass in 12.2 s.
  - Reruns are bit-identical to the raw in every stored field: `a2e_job` on both spans (rows and LSODA steps), all six `gt_job` rows, and `sq_job` on targets 12 and 27 (all 8 masses each, 16 cells).
- **A2E's saddle is real.**
  - At the symmetric state, central differences of the right-hand side give ∂ċ/∂c equal to the exact Jacobian to 6 digits (+1.0705·10⁴ … +1.0889·10⁴).
  - The c-row decouples (off-diagonals ≤ 10⁻⁶ relative). The eigenvalues are −1.48·10⁴ (σ), −1.2·10⁻³ (α) and +1.09·10⁴ (c).
- **A2E's route verdict for 5A at 0.25 holds.** Every independent integrator reaches the plateau at s = 0.0056 from the exact start on the span (0, 10):
  - Radau (rtol 10⁻¹²);
  - scipy BDF (rtol 10⁻¹²);
  - DOP853 (rtol 10⁻¹³);
  - RK45 (rtol 10⁻¹¹);
  - LSODA at 10⁻¹²/10⁻¹⁶.
  - So does LSODA with A2's tolerances from any start offset by ≥ 10⁻¹⁵σ* (at 0.0032 from 10⁻⁸σ*). Only the exact start under A2's LSODA on the short span stays [own].
- **The plateau is the two-sample interpolation.**
  - E₂ = 0.643650, 0.906294, 0.976497 and E₁ = 0.8119, 0.9530, 0.9882, from my own formulas.
  - min_c J_N at σ = 0.002 equals E₂ to 6 digits. The plateau E lies within 6.1·10⁻⁵ of E₂ [own].
- **Recomputed lines** [own]:
  - GL's one-energy remainder: −6.26·10⁻⁵ … +1.17·10⁻⁴, largest at τ = 574.3.
  - The own-crossing bracket: inside on 9, 18, 20, 27, 29; on 12 it lies +0.251 … +0.349 % above the path's κ.
  - The rate line and its exponents.
  - SQ: 34 of 48 cells on the path's sequence; τ₃ and τ₄ equal G's at 48 of 48. SQ's τ₃/τ₄ rule is `fine_phases`' by construction.
  - The HLL per-decade line.
  - The dissolution times: 144.977, 154.129, 61.443.
  - `fit_split` equals my brute-force best assignment on all 94 rows of H, HP, H3, HL and HLL.
  - Every number on D59.md's changed lines and in the README row traces to the md or the raw; I found no misquote.
- **The symmetric pair's no-cancellation is not a symmetry artefact** [own]. At ε = 10⁻⁸ it survives every perturbation of 10⁻⁹–10⁻⁶ that I tried, to τ = 10⁵:
  - in the amplitude, the centre, and the dictionary start;
  - but A₂ = 1 − 10⁻⁴ gives an asymmetric deficient state;
  - and A₂ = 0.999 a cancellation (+1.39 → +1.48).
- **Round 83's behaviour-changing survivors are now killed.** HJ1–HJ3 and GJe fail the tests; its ten "equivalent" ones still survive. This matches the log.

**What does not (majors first)**

1. **D3ch-1 (major): A2's 3A at σ*/σ₀ = 0.25 row is the A2E artefact itself.**
   - Its symmetric state on the sample c* is a saddle: ∂ċ/∂c = +629 … +665 at σ ≈ 0.235σ*, from s ≈ 0.003 to beyond 200.
   - LSODA's BDF steps hold it there. The flow goes onto the plateau, like 3A at the other two ratios.
   - So "(5A at every width ratio)", "off the plateau" and "σ to 0.0049–0.092σ*" are wrong. "A2's route is the flow's" fails within A2, and a test asserts the artefact.
2. **D3ch-2 (minor): A2E's damping mechanism is misdescribed.**
   - The first BDF steps are order 2 at hλ 2.4; they amplify the offset ×2.5 per step, to its largest value.
   - Only the order-1 steps at hλ 4–50 damp it.
3. **D3ch-3 (minor): the split criterion's sign is Σw's, which cannot tell the ends apart.**
   - HP's equal end and HLL's symmetric 10⁻⁸ cell fit none of the three classes.
   - Σw flips sign within one cell.
   - H's "exact, null cluster" fails the stated null test.
   - A single broad kernel is labelled one feature's "main kernel" by tie-break.
4. **D3ch-4 (minor): "approaching a null cluster" is asserted as fact.** By the md's own test it is not being approached: the others' sum over the residual goes 9.1·10³ → 4.2·10³, decelerating.
5. **D3ch-5 (minor): tests.**
   - Of 121 single-point mutants of this round's code, 35 pass the tests (24 caught only by the md byte comparison), and 3 behaviour-changing ones pass everything: SA3, SJ4, AE12.
   - GL1 reverts D3cd-1's fix and passes all 46 tests.
6. **Cosmetic.** D3ch-6: SQ's off-sequence categories misdescribe 8 of 14 cells. D3ch-7: the κ′ rate wording.

## What I checked and how

All paths are under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref87/D3ch_work.

1. **Reproducibility.**
   - `--from-raw` into W/mwout, then cmp with the committed md.
   - pytest: W/logs/pytest_full.log.
   - GT and SQ reruns: W/scripts/rerun_gt_sq.py → W/logs/rerun_gt.log, W/logs/rerun_sq.log.
   - Raw overview: W/scripts/explore.py → W/logs/explore.log.
2. **Independent recomputation from the raw** (own code): W/scripts/recompute.py → W/logs/recompute.log.
   - SQ sets rebuilt from the event lists; κ′, the bracket and the rate.
   - GL's remainder three ways; A2 E₂/E₁ and J_N minimised on a fine c-grid.
   - A2E's saddle and step lines; the dissolution times.
3. **A2E.** W/scripts/a2e_check.py → W/logs/a2e_check.log:
   - the rerun against the raw;
   - an FD Jacobian with eigenvalues at every symmetric output;
   - LSODA steps with method (iwork[18]) and order (iwork[13]);
   - the BDF-k amplification roots;
   - Radau, BDF, DOP853, RK45 and tight LSODA on both spans, with perturbed starts.
4. **A2 robustness.** W/scripts/a2_robust.py → W/logs/a2_robust.log covers all 15 cells:
   - the early phase on c = c*, with ∂ċ/∂c there;
   - LSODA from c* ± 10⁻¹²σ*, +10⁻⁸σ*, +10⁻⁴σ*;
   - Radau and BDF at A2's tolerances.
5. **3A at 0.25 in detail.** W/scripts/a2_3a.py → W/logs/a2_3a.log:
   - LSODA's steps on the saddle;
   - DOP853 and RK45 from the exact start;
   - LSODA from δ = 10⁻¹⁶ … 10⁻⁶σ*;
   - E along c at the symmetric state.
   - Solver failures at the plateau's escape: W/scripts/a2_fail.py → W/logs/a2_fail.log.
6. **(ii)(b).**
   - W/scripts/split_check.py → W/logs/split_check.log: brute-force assignment on all 94 rows; the criterion applied; the others' sum at the feature centres; the mass Σw|σ|√(2π).
   - W/scripts/hll_ratio.py → W/logs/hll_ratio.log: the others' sum against the residual.
   - W/scripts/hll_stale.py → W/logs/hll_stale.log: the stored HLL split fields are round 80's, so recomputing at render time is required.
   - Symmetric pair: W/scripts/sym_pair.py → W/logs/sym_pair_1e-8.log (perturbed targets and dictionary); W/scripts/sym_jac.py → W/logs/sym_jac.log (the loop's Jacobian with mirror parity).
7. **Mutants.**
   - W/scripts/mutants.py holds 121 single-point mutants of this round's new or changed code. W/scripts/run_mutants.py runs each in a slim tree W/mut/<id>, with `-x` on the test file and with `--from-raw` + byte comparison. Results: W/logs/mut_all.log, W/logs/mut_results_all.json; baseline: W/logs/mut_results_base.json.
   - Round 83's 14 survivors: W/scripts/mutants_r83.py and W/scripts/run_mutants_r83.py → W/logs/mut_r83.log.
   - Job mutants that `--from-raw` cannot see: W/scripts/sq_mutants.py → W/logs/sq_mut_SJ4.log, W/logs/sq_mut_SA3.log.
8. **Literature.** Two web searches; no page fetched.

## Verdicts per part

| claim (D59.md lines) | verdict | evidence |
|---|---|---|
| (i) the landscape; convergence from E(start) < ½‖f*‖² (8–15; proof 120–127) | holds | unchanged |
| (i) "from 2.5A it dips to 0.3σ*" (17; 222–223) | holds | the symmetric state is stable (∂ċ/∂c < 0 at every output); robust to every perturbation and solver [own] |
| (i) "from 3A on (5A at every width ratio) it collapses at once onto the two samples" (17–18); "from 3A at σ*/σ₀ = 0.25 to 0.235σ* (at the outputs), off the plateau … (σ to 0.0049–0.092σ*)" (223–224) | **wrong for 3A/0.25** (D3ch-1) | saddle λ ≈ 650. DOP853, RK45, Radau and perturbed LSODA all give the plateau at s = 1 … 3·10⁵ (E 0.644, σ 0.111σ*), escape by 10⁶ [own] |
| the plateau = lim_{σ→0} min_c J_N = E₂; E₁ (18–19; 224–225; Status 173) | holds | own formula; min_c J_N(σ = 0.002) = E₂ to 6 digits [own] |
| durations 3·10⁵, 3·10⁷, 10¹⁰/3·10¹⁰; escapes; 20A at 10¹³; ends 1.8–2.6·10⁻³² (20–21; 225–227) | holds | [raw]; reproduced from perturbed starts [own] |
| A2E: the symmetric state at s = 10⁻³ (σ/σ* 0.1846, E 3.480); ∂ċ/∂c +1.05 … +1.09·10⁴ and −383 … −279; ×e per 9.2·10⁻⁵ (227–229) | holds | FD = exact; eigenvalues −1.48·10⁴, −1.2·10⁻³, +1.09·10⁴ [own] |
| A2E: Adams steps hλ ≤ 0.55 grow at 1.092·10⁴; plateau from 0.0056 (229–230) | holds | 10917 per unit s [own] |
| A2E: "BDF at s = 0.0026 (hλ 2.4–50), whose steps damp the unstable offset (backward Euler's factor …)" (230–231) | **partly wrong** (D3ch-2) | order 2 at hλ 2.41: ×−1.67, −2.54, −2.43; order 1: ×−0.335, −0.152 [own] |
| "A2's route is the flow's, the span-10 run's a solver artefact" (231–232) | holds for 5A/0.25; **fails for A2's 3A/0.25** | five integrators and perturbed starts [own] |
| (ii)(b): the fit's split by the best assignment, signed w/A − 1, the others' Σw (41–46) | computation holds; **criterion incomplete** (D3ch-3) | brute force = `fit_split` on 94/94 [own] |
| HP and HL/HLL asymmetric numbers (46–53; 299–305; 309–322) | holds | traced [md] |
| symmetric pair: main kernels deficient, no cancellation in the loop (54–56; Status 178–179) | holds, robust | perturbations 10⁻⁹–10⁻⁶ keep it at 10⁵; A₂ 1 − 10⁻⁴ asymmetric, 0.999 cancellation [own] |
| "(Σw +0.036 → +0.033: a split)" (56) against "approaching a null cluster, not one at 10⁹" (61–62) | the contrast is unsupported (D3ch-3, -4) | same structure in both cells; others/residual 9.1·10³ → 4.2·10³ [own] |
| (2) dissolution times to five digits (87; 259) | holds | 144.977, 154.129, 61.443 [own] |
| (3) κ′ named; converges where the loop forms the path's kernels; O(ε^{1/3}) asymptotically; 12 pre-asymptotic, own crossing 0.25–0.35 % above (95–99; 262–276) | holds; "against 10^{1/3}" generous for 20 and 29 (D3ch-7) | [own] |
| SQ: 48/48 and 34/48; the categories of the other 14 (269–274) | counts hold; **categories misdescribed** (D3ch-6) | [own] |
| (4) parity "in the continuum" (94–95; 135–137) | holds | wording |
| GL: one-energy remainder −6.3·10⁻⁵ … +1.2·10⁻⁴; slopes 0.951–0.963 / 0.951–1.052 (279–287; Status 185) | holds | [own] |
| New, literature (169–171, 182) | fits | [snippet] |
| README D59 row | fits, except that A2E's summary omits 3A/0.25 and "as ε^{1/3} on four targets" is generous | — |

## Status of the round-83 issues

| D3cd | status | what remains |
|---|---|---|
| 1 (GL's remainder) | fixed | One energy throughout; I get −6.26·10⁻⁵ … +1.17·10⁻⁴. But mutant GL1, which restores the mixed slope, passes all 46 tests (D3ch-5). |
| 2 (`fit_split`) | partly fixed | The best assignment, the signs and Σw are in, and the HL A₂ 0.99 rows now pick the +1.58 kernel. The criterion's sign is Σw's, which fails on the symmetric pair; the null test mixes norms; a broad kernel becomes "main" by tie-break (D3ch-3). "Approaching" is still asserted (D3ch-4). |
| 3 (κ′, rate) | fixed | Named, bracketed, rated. Wording only (D3ch-7). |
| 4 (A2 plateau, saddle) | partly fixed | E₂, `a2_spike` by E and "(at the outputs)" are right. The damping mechanism is wrong (D3ch-2). The same saddle-hold sits in A2's 3A/0.25 row (D3ch-1). |
| 5 (surviving mutants) | fixed | HJ1–HJ3 and GJe now fail the tests; the ten equivalent ones survive, as expected. New gaps (D3ch-5). |
| 6 († a proxy) | fixed | SQ works by index, † comes from SQ, and '—†' is printed. Category wording (D3ch-6). |
| 7 (wording, rounding, literature) | fixed | Nothing. |

## Issues

| id | severity | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ch-1 | major | **A2's 3A at σ*/σ₀ = 0.25 is held on the saddle by LSODA (D59.md:17–19, :222–224, :231–232; md row "0.25 \| 3"; test AS6 in `test_a2_one_kernel_from_large_coefficients`).** <br>• **The saddle.** In A2's own run the kernel sits on c = c* from the start. ∂ċ/∂c at c* turns positive at s ≈ 0.0032 (σ/σ* 0.274) and stays positive at 12.8 … 665 to beyond s = 204 (+629 … +665 at the outputs, σ/σ* 0.2345–0.2354). Over LSODA's steps ∫λ ds ≈ 7.3·10⁴. <br>• **How LSODA holds it.** The Adams steps grow the offset to 4·10⁻¹⁴σ* by s = 0.0105. BDF from s = 0.0128 (hλ 0.7 → 7 → 2·10³) damps it back to roundoff by s = 0.088. <br>• **Energy.** E is maximal in c there: at s = 3, E(c* + 0.165σ*) − E(c*) = −0.0375, with E(c*) = 1.205 of ½⟨f*²⟩. <br>• **The flow.** From the exact start, DOP853 (rtol 10⁻¹²) and RK45 leave c* by s = 0.1 and sit on the plateau (offset −0.165σ*, σ 0.1111σ*, E 0.6436). Radau at A2's tolerances, and LSODA from c* + 10⁻¹², 10⁻⁹ or 10⁻⁶σ*, give the plateau at s = 1 … 3·10⁵, escape by 10⁶, an exact end — as 3A at ratios 1 and 4. Only scipy BDF and LSODA (from the exact start, or from 10⁻¹⁵ and 10⁻¹⁴σ*) stay. <br>• **Other cells.** Perturbed starts reproduce A2's other 14 rows. | **A2:** <br>• rerun from a start off the sample (c* + 10⁻⁹σ*), or check each cell with an explicit integrator; <br>• extend A2E to 3A/0.25. <br>**(i) and the Numerics:** <br>• "from 3A on, at every width ratio"; <br>• "σ to 0.0049–0.111σ*"; <br>• drop "off the plateau"; <br>• say which runs the saddle held. <br>**Tests and README:** <br>• invert AS6; <br>• README: "held by LSODA on a short span, and in A2 at 3A/0.25". |
| D3ch-2 | minor | **A2E's damping (D59.md:230–231; md lines 53, 55).** <br>• **The amplifying steps.** LSODA's first three BDF steps are order 2 at hλ = 2.41. They multiply the offset by −1.67, −2.54 and −2.43 (BDF2's parasitic root at z = 2.44 is −2.35). That takes it from 2.1·10⁻¹² to its largest value, 2.2·10⁻¹¹σ*, under BDF. <br>• **The damping steps.** The order-1 steps at hλ 3.99 and 7.56 damp it by ×−0.335 and ×−0.152 = 1/(1 − hλ); hλ 49.5 takes it to roundoff. <br>• **The general rule is false.** "An implicit step at hλ > 2 damps" holds for backward Euler only. BDF-k is unstable on 0 < hλ < 2, 4, 6.67, 10.67, 17.07 (k = 1…5) of the positive real axis [own]. <br>• **The growth rate's thin base.** On the span-10 run the offset sat at +6.9·10⁻¹⁶σ* for 45 Adams steps (s 1.01–1.60·10⁻³), the roundoff floor. The quoted growth of 1.101·10⁴ rests on two steps. | **`a2e_job`:** record the order (iwork[13]). <br>**D59 and md:** "its order-2 steps at hλ 2.4 amplify the offset ×2.5, its order-1 steps at hλ 4–50 damp it (1/(1 − hλ))". |
| D3ch-3 | minor | **The split criterion (D59.md:41–46, :56, :61; `fit_split`).** <br>• **(a) Ends that fit no class.** HP's equal end and HLL's symmetric 10⁻⁸ cell have deficient main kernels and others with Σw −1.05 … −0.96 (the other sign). Yet the others' sum at the feature centres is positive (+7.5·10⁻² → +2.0·10⁻⁵ of max f*) and their mass √(2π)Σw\|σ\| is +0.18 (one broad positive kernel, w 0.45, σ 1.6–1.95). By the stated criterion these ends are neither split, cancellation nor null cluster. <br>• **(b) Σw's sign flips within one structure.** HL equal at 10⁻⁶: Σw −0.043, +0.000, +0.036 at τ = 6400, 2.5·10⁴, 10⁵, for one structure (a positive kernel at c = 0, σ ≈ 0.18, carrying the deficit). So "a split" at 10⁻⁶ against "approaching a null cluster" at 10⁻⁸ rests on Σw's sign alone. <br>• **(c) The null test.** H's σ* 0.04 end is "exact, null cluster" by `outcome`, but fails the stated test: the others' sum is 1.0·10⁻¹¹ of max f* against a residual of 1.3·10⁻¹³ (a sup norm against an rms, with no tolerance). <br>• **(d) Broad kernels as main kernels.** In HL's transient rows a single broad kernel (c ≈ 0, σ 4–8σ*) becomes one feature's main kernel, with error 3–7. Which feature gets it is the tie-break: mutant FS10 (`<=`) moves it in 13 md lines and changes w/A − 1 on A₂ 0.99 from −0.026 to −0.016. | **Criterion:** <br>• judge the others' sign by their sum at the features, or by their mass; <br>• define the null cluster with one norm and a tolerance. <br>**`fit_split`:** <br>• leave a feature unassigned when the best error exceeds 1. <br>**Text:** <br>• restate (ii)(b) accordingly. |
| D3ch-4 | minor | **"Approaching a null cluster, not one at 10⁹" (D59.md:48, :61–62).** <br>• **The md's own test is not being approached.** HLL symmetric 10⁻⁸: the others' sum over the residual (both over max f*) goes 9.1·10³, 7.0·10³, 5.5·10³, 4.7·10³, 4.2·10³ (×1.30, 1.27, 1.18, 1.11 per decade). <br>• **Both decelerate together.** The sum's local exponents are −1.21, −0.93, −0.79, −0.67 per decade, and the residual's −1.10, −0.82, −0.72, −0.62. The others carry the main kernels' deficit (sum/error 1.23 → 1.59). <br>• **What the data support** is the Status's "undecided at τ = 10⁹". | **(ii)(b):** <br>• "falling, decelerating; undecided at 10⁹"; <br>• print the ratio to the residual; <br>• optionally locate the symmetric window: no cancellation to A₂ = 1 − 10⁻⁶, a cancellation at 0.999 [own]. |
| D3ch-5 | minor | **Tests (the round's code).** <br>• **Totals.** Of 121 mutants the tests kill 86 and the md alone 24 more. 11 survive both. <br>• **By area, among the round's new tests:** A2E's saddle and step lines (12 of 15 `a2_lines` mutants pass the tests); the HLL per-decade line (6 of 7); GL1 (D3cd-1's fix reverted; the test's 10⁻²³ floor against E ~ 10⁻²¹ cannot separate the two at .1e); GL3, GL5; FS12 (main_dw over As[0]); SC3 (the test data have main_dw = main_err). <br>• **Behaviour-changing, passing everything:** SA3 (`sq_args` passes the widths as centres; target 9's events change entirely); SJ4 (the formed threshold by max\|A\|; on a rerun SQ's τ₃/τ₄ equal G's at 32 of 48); AE12 (span order). <br>• **Equivalent:** AS1, AS2, FS4, HR4, SJ1, SK7, SK8, GL4. <br>• **No full-md test.** No test compares the full md with the raw. | **Tests:** <br>• a test rendering `report` from the committed raw and comparing bytes with the md (0.6 s; kills all 24 md-only mutants); <br>• `sq_args` field by field; <br>• `sq_job` with two amplitudes; <br>• GL's floor ≥ 1 % of E; <br>• `fit_split` with unequal amplitudes and a broad-kernel tie. |
| D3ch-6 | cosmetic | **SQ's other 14 (D59.md:271–274).** <br>• **Kernel formed before τ_d, not at τ₃.** At 18 at 10⁻⁵ and 10⁻⁶, and 27 at 10⁻⁵ and 10⁻⁶, kernel 9 forms before τ_d (23.5, 25.0, 25.3, 25.8 against 25.003 and 31.54). 12's kernel 10 also forms before τ_d (15.2, 23.2 < 24.35). There the loop's "τ₃" (25.1, 31.6, 24.4) is just the first output past τ_d. <br>• **9 at 10⁻⁵.** It has no extra kernel: it is the path's j4 (11) at 22.0, then j3 (3) at 26.2. <br>• **The slopes.** "4.40–4.47 and −25.2" covers only the four cells with a landing. <br>• **The proxy sentence.** Nine of the 14 cells lack a landing; the proxy misses three of them (20 at 10⁻⁵ and 10⁻⁶, 29 at 10⁻⁵). | Reword per cell type: before τ_d; at τ₃; reversed order; τ₄ entrant. |
| D3ch-7 | cosmetic | **Rate wording (D59.md:96–98, :274–276; README).** <br>• **(κ′ − κ)/ε^{1/3} by target:** 9 −8.6 → −5.8; 20 −1.8 → −2.7 (saturating); 27 −5.9 → −5.6; 29 −12.6 → −25.9 → −21.1; 12 +3.7 → +19.6, still rising. 18 changes sign (+0.0035, +0.0000, −0.0005). <br>• **The README.** "Falling as ε^{1/3} on four targets" fits 20 (ε^{0.27–0.31}) and 29 (ε^{0.07–0.36}) loosely. | Report (κ′ − κ)/ε^{1/3}. |

## Sources

- **[own]** The BDF-k instability intervals (0, Σ_{j≤k} 2^j/j) = (0, 2), (0, 4), (0, 6.67), (0, 10.67), (0, 17.07), from the boundary locus at ρ = −1, and the BDF2 roots. Background on BDF stability regions [snippet]:
  - [Lóczi, stability angles of BDF](https://www.math.u-szeged.hu/12QTDE/abspdf.php?tid=69)
  - [Akrivis–Katsoprinakis](https://www.cs.uoi.gr/~akrivis/AK2a.pdf)
  - [Auzinger–Herfort](https://repo.agh.edu.pl/items/32e2694a-7c8d-466c-91c4-89e67cce2588/full)
- **[snippet]** "Dendrogram of mixing measures: Hierarchical clustering and model selection for finite mixture models" (Do, Do, McKinley, Terhorst, Nguyen): [arXiv:2403.01684](https://arxiv.org/pdf/2403.01684), [ar5iv](https://ar5iv.labs.arxiv.org/html/2403.01684). The r̄(2) = 4 attribution is round 83's snippet, not re-seen this round.
- **[memory]** The draft's other references, unchanged this round.

## Lesson

This round's best idea — when two solver runs disagree, read the steps — was applied only to the run that disagreed with A2, never to A2 itself. The artefact's preconditions were a start on the sample c*, a saddle that forms as the width shrinks, and a stiff block that sends LSODA to large BDF steps. Every A2 cell has the first two, and in 3A at 0.25 a slow instability (λ ≈ 650) gives LSODA time to switch before the offset grows, so the solver keeps one of A2's own rows on the saddle. A diagnosed artefact should be checked on every run that shares its preconditions, using the cheapest discriminator available here: a start offset by 10⁻¹²σ*, or one explicit integration. Two related lessons:
- An explanation should be checked against the solver's actual state. The BDF order was one integer away in iwork, and recording it would have shown that the steps quoted as damping were order 2, which amplified the offset.
- A classification sign should be taken where it acts. The others' sum at the features tells a split from a cancellation; their summed coefficients are set by kernels far outside the window.
