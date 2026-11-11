# Referee report D3ce, round 84: Proposition D‴.58 as revised after round 81

## Scope

- **Commit.** 902fa4a ("D‴.58 round 81: B4 evaluated exactly; the anchored residual reaches the optimum to 10⁻¹²"), branch claude/festive-bardeen-eopej2.
  - I exported it with `git archive` into `W/repo`, where W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref84/D3ce_work`.
  - Everything I wrote is under W:
    - `scripts`: 30 scripts.
    - `logs`, `out`.
    - `mut`: one slim tree per mutant.
    - `wmout`: WM_OUT, holding a copy of the raw.
- **Files read in full.**
  - `review/drafts/D58.md`, `experiments/exp_width_mismatch.py`, `results/chaos/width_mismatch.md`.
  - `tests/test_width_mismatch.py` (44 tests) and `src/core/wide_feature.py`.
  - `review/REVIEW_round81_D3cb.md` and the Round 81 entry of `review/REVIEW_LOG.md` (:3183–3218).
  - `git diff beeefe5 902fa4a` of the script and the tests, and a word diff of D58.md.
- **Files read in part.**
  - `src/core/rich_transit.py`: `Bumps`, `polish_block`, `ascent_target`, and the polish in `sliding_iss` (:395–460).
  - The raw: all of B3, B4 and A5, and B's stored first polishes.
- **Run settings.**
  - venv Python 3.13.16, numpy 2.5.4, scipy 1.18.1. Long double is x87 (ε 1.08·10⁻¹⁹).
  - Every process ran with `nice -n 10` (15 for sub-minute evaluations), OMP_NUM_THREADS=1 and PYTHONDONTWRITEBYTECODE=1.
  - pytest ran as `-q -p no:cacheprovider`; mutants added `-x` and a `--basetemp` under W. My first full test run used pytest's default basetemp under /tmp/pytest-of-root; pytest has since rotated it away.
  - I never ran more than two heavy processes at once.
  - Exact arithmetic used `decimal` and `fractions` only.
- **Nothing in the repository was modified.**
  - I ran only `git archive/log/show/diff/status`.
  - These files have the same md5 in the repository as in my export: D58.md, the script, the md (1130e5fc…), the raw (efcc0825…), the tests, `wide_feature.py`, `rich_transit.py`, `exp_cohesion_leak.py` and the round-81 report.
  - The repository's `__pycache__` files for this experiment predate my session (02:31).
  - Meanwhile other sessions committed 6600db5 and 3773079, and left `REVIEW_LOG.md`, `D52.md` and the drafts README modified. None of this touches a D‴.58 file.
  - Every process I started has exited.
- **Marking.** [md] the committed md, [raw] the committed raw, [mine] my computations, [snippet] web-search snippets, [memory] recalled facts.

## Summary

**What checks out**
- **Reproduction [mine].**
  - `--from-raw` regenerates the md byte for byte.
  - 44/44 tests pass in 39.9 s.
  - All 12 B4 cells rerun bit-identical to the raw in every field, at 2–348 s each.
  - `part_A5` reruns bit-identical in 6.4 s.
- **The evaluator is exact [mine].**
  - I wrote an independent evaluator: each end's p and p_lo as exact Fractions (hi + lo), exact exponents, exp at 100 digits.
  - It reproduces all 72 stored B4 ratios bit for bit.
  - `_dec`'s split is exact: the lo part has at most 11 bits, so v − hi and its conversion to double are both exact.
- **B4's headline numbers hold exactly [mine].**
  - The anchored double ends within 1.28·10⁻¹² of the best end from all 12 starts.
  - A block's two starts end within 2.68·10⁻¹³ of each other.
  - 'nostall' ends up to 9.81·10⁻⁹ above.
  - From the blocks, direct subtraction ends 1.93·10⁻⁸ … 2.33·10⁻³ above in double and up to 2.20·10⁻⁸ above in long double.
  - `lm_rms` errs by −2.4·10⁻¹⁰ … +1.7·10⁻¹⁰; the extremes fall at the two blocks below 10⁻¹⁰.
- **The best end is a sound reference [mine].**
  - It is the 'anchored ld' end at all six blocks.
  - From the two starts those ends agree to 3.5·10⁻²¹ … 2.1·10⁻¹⁹, exactly.
  - A Gauss–Newton step from the exact residual predicts at most 1.4·10⁻¹⁷ of further decrease.
- **`lm_fit` is the solve described.**
  - My instrumented copy is bit-identical at 12/12.
  - It shows 19–1622 stall re-anchors per solve, each with μ reset, and 0–76 periodic re-anchors counted across them. The flat stop is live at 10 of 12.
  - The frozen-rounding mechanism is right: scaling the anchor's rounding by κ = ½ … 4 scales the end's excess as κ² (slopes 2.00 and 2.05 at the two clean rows).
- **B3 [raw, mine].**
  - Each row is judged against the best end over its own optimum.
  - The guard sits in `b3_restart`, right after the 'path' restart, and compares acceptance.
  - The 22 'path' restarts equal B's stored first polishes in status, evaluations and acceptance (22/22).
  - The labels and 1.0235/1.0484 hold.
- **A5's series is the rule's squared error [mine].**
  - I computed the exact L² error of `scaled_rule` independently: closed-form Gaussian inner products in Decimal, Christoffel-sum weights.
  - It is a constant multiple of S(K) to ≤ 7·10⁻¹⁵ (K ≤ 40, τ/σ ≥ ½).
  - It equals the 'far' ratios at K = 60–120 to every printed digit.
  - The quoted limits 17.165, 5.3094 and 2.33333 are right to their digits.
- **Numbers [md].** Every number D58.md changed this round traces to the md or the raw.

**What does not (majors first)**
1. **D3ce-1 (major): "linear in its rounding" (D58 :73, Status :162–163) fails a scaling test.**
   - The double's stall excess over the long double's is 1.1·10³ … 3.3·10⁷ (median 1.75·10⁵). A floor linear in the rounding would give 2048.
   - Scaling the residual's rounding by κ, the stall excess grows like κ^1.2–1.9.
   - Both floors grow superlinearly; what separates them is size (10³–10⁶), not order.
2. **D3ce-2 (major): 20 of 87 single-point mutants of this round's code survive both checks; 11 change published results on a rerun.**
   - `a5_far` is run but never checked.
   - The half-squeeze frame, `frame_series`' stop and the low-defect metric are unpinned.
   - LF02, which disables the flat stop's period, hides behind the stall-only pinned cells.
3. **Minors.**
   - **D3ce-3:** A5's "closed form open" is superseded. One formula, r(v) = max{h(τ²/v), h(vσ*²/(σ²τ²))}, fits all of (ii)'s rates and A5's limits. At τ/σ = ¼ the Mercer frame's rule also ends slower than t (65.056 < 65.985).
   - **D3ce-4:** `anchor_floor` takes e₀ at the end, not at the last anchor (×0.23–7.7 against the gap, versus ×0.40–2.0 at the last anchor). Part of e₀ is systematic (2σ² rounded). At the lowest optimum the double parameters alone cost 1.5·10⁻¹².
   - **D3ce-5:** the smallest B4 differences are the double ratios' rounding (2.2·10⁻¹⁶ quoted; 9.5·10⁻¹⁷ exactly).
   - **D3ce-6:** `frame_series` stops at a zero term at K = 1. This produces the md's "1.5e-11".
4. **D3ce-7:** cosmetics.

## What I checked and how

- **Reproduction** (`b4_rerun.py`, `a5_rerun.py`).
  - `--from-raw` into W/wmout, then `cmp` with the committed md.
  - The test file.
  - The real `b4_job` on all 12 cells (two workers), every field compared with `==` against the raw.
  - `part_A5` on the raw's A, compared with the raw's A5.
- **The evaluator** (`exact_ends.py`, `exact_gaps.py`).
  - Each end rebuilt from the raw's p and p_lo as Fractions; exponents −(x − c)²/(2σ²) exact; `Decimal.exp` at 100 digits.
  - The ratio formed with the script's own roundings.
  - Exact excesses, gains and two-start agreement.
  - A Gauss–Newton step from the exact residual at every end.
- **The anchored solve** (`lm_inst.py`, `anchor_last.py`, `floor_decomp.py`).
  - An instrumented copy of `lm_fit` recording every anchor (point, kind, iteration). Its end, iterations, accepted steps and stop are bit-identical to the raw at 12/12.
  - Floors computed from e₀ at the end and at the last anchor, against the actual gap.
  - e₀ split into the 2σ²-rounding term (an exact difference of two Gaussians) and the rest.
- **The floors' order** (`kappa_scan.py`, `kappa_robust.py`).
  - The residual's rounding scaled by κ: r_κ = r_ld + κ(r_d − r_ld), so κ = 1 is the double residual.
  - Direct runs: κ = 1/16 … 4 from the stored optima, and five starts per row.
  - Anchored runs: the anchor's r₀ = r_κ for κ = 0 … 4. At κ = 1 they reproduce B4's rows exactly.
  - The best end rounded to double, as the parameters' own floor.
- **B3** (`b3_refs.py`, `b3_guard_check.py`): per-row references and labels on the raw; the 22 'path' restarts against B's stored first polishes.
- **A5** (`scaled_rate.py`, `scaled_rate2.py`, `a5_validate*.py`, `a5_k1.py`, `fit_limits.py`, `gh_exp_asym.py`).
  - The scaled rule's exact L² error, up to 652 digits. Gauss–Hermite (GH) nodes by Newton from numpy's; weights by the Christoffel sum, independent of `gh_nodes`.
  - e²/S as a function of K.
  - The K = 1 sum carried out in full.
  - Limits at 10 (T, v) points.
  - The asymptotics of the GH error for exponentials.
- **Mutation testing** (`mutants_spec.py`, `mut_run.py`, `probe.py`, `probe_lines.py`).
  - 87 single-point mutants of this round's new or changed code, plus a no-op control.
  - Each runs in a tree under W/mut: the mutated script, `src` symlinked, the test file, and the raw in both its WM_OUT and `results/chaos`.
  - Each runs `pytest -x` (600 s timeout) and `--from-raw` with a byte comparison.
  - Survivors were probed:
    - their `part_A5` against the raw, with the md regenerated;
    - their `b4_job` on cells 9, 10 and 11.

## Verdicts per part

| claim (D58.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii), (iv), (v), the proof (:8–23, :98–147) | Unchanged this round; not re-derived | word diff |
| (ii) λ, λ/(2 − λ), q "numerical (A3)" (:15–23, :158–159) | Consistent with one closed form | D3ce-3 [mine] |
| (iii) 1.0235, 1.0484 over the best end; TRF from the path's blocks 1.5·10⁻², 5.3·10⁻³ on its step test (:55–59) | Correct (exactly 1.48·10⁻², 5.29·10⁻³; status 3) | [md], [mine] |
| (iii) "not the arithmetic"; lm_rms −2.4·10⁻¹⁰ … +1.7·10⁻¹⁰, low at 12/12 (:59–62) | Correct (−2.45·10⁻¹⁰ … +1.68·10⁻¹⁰ at the two blocks below 10⁻¹⁰ alone) | [raw] |
| (iii) direct subtraction 1.9·10⁻⁸ … 2.3·10⁻³; long double up to 2.2·10⁻⁸ (6.6·10⁻⁹, 1.4·10⁻⁸; 2.2·10⁻⁸) (:62–65) | Correct (1.93·10⁻⁸ … 2.33·10⁻³; 6.58·10⁻⁹, 1.39·10⁻⁸, 2.20·10⁻⁸) | exact_gaps [mine] |
| (iii) the anchored solve as described: stall re-anchors with μ reset, a 2000-step period across them, the flat stop (:65–68) | Correct: it is the solve that runs | instrumented copy, bit-identical |
| (iii) within 1.3·10⁻¹² from all 12; two starts within 2.7·10⁻¹³; 10 flat, 2 stall; 8011 … 1.6·10⁵ iterations (:68–69) | Correct (1.28·10⁻¹², 2.68·10⁻¹³). The best end is the 'anchored ld' end, so this and :71's gain are one measurement | [mine]; D3ce-7 |
| (iii) the floor "quadratic in the frozen rounding"; predicted 1.8·10⁻¹⁶ … 8.8·10⁻¹³; gain 2.2·10⁻¹⁶ … 1.3·10⁻¹² (:69–71) | Mechanism correct (κ² for κ ≥ ½). The prediction uses the end's rounding. 2.2·10⁻¹⁶ is rounding (exactly 9.5·10⁻¹⁷). At (¼, 6) the parameters' double rounding is as large | D3ce-4, D3ce-5 |
| (iii) 'nostall' up to 9.8·10⁻⁹; "the anchored residual, not the restarts" (:71–73) | Correct (9.81·10⁻⁹). Clean at the optimum starts (no re-anchor at all); from the blocks it still re-anchored periodically | [raw]; D3ce-7 |
| (iii) "linear in its rounding" (:73) | Contradicted | D3ce-1 |
| (iii) the non-monotone examples; from the stored optima ≤ 2.3·10⁻⁵ (:73–75) | Correct ("the same 1.6·10⁻⁹" is two blocks: 1.6114 and 1.6069·10⁻⁹) | [md] |
| (iii) stored optimum 1.6·10⁻³ above the best end at ¾, K = 11; the 0.1% test admits 2.6·10⁻³; ≤ 3.1·10⁻⁵ elsewhere (:75–78) | Correct (1.629·10⁻³, 2.631·10⁻³) | [md], [mine] |
| (iii) B3 against the best end: 6.3·10⁻⁹, 1.7·10⁻⁶, 2.7·10⁻⁵, 2.3·10⁻⁴, 2.3–4.8% (:78–80) | Correct, with per-row references | b3_refs [mine] |
| Status (iii) (:159–164) | "Linear in its rounding in double and in long double" contradicted. "≲ 10⁻¹²" holds but is also the parameters' double floor | D3ce-1, D3ce-4 |
| Status A5: the series; "reproduces the rule's actual errors to 6·10⁻¹⁴ … 3·10⁻⁵" (:187–189) | The series is the squared error (constant factor to ≤ 7·10⁻¹⁵). The quoted numbers are the trapezoid's limits | a5_validate [mine]; D3ce-7 |
| Status A5: first term qv′/2 = 0.375; tail above the first from K = 19, 8, 6; 0.38 at ¼ (:189–191) | Correct | [md] |
| Status A5: "converges, as K⁻² (to K = 120), to 17.165, 5.3094 and 2.33333" (:191–192) | Correct to the digits quoted (closed form: 17.16497, 5.309401, 7/3). K⁻² is clean at 2 and 1 only | D3ce-3, D3ce-7 |
| Status A5: "at ¼ it still rises at K = 40"; "slower … at ½, 1 and 2, its rate's closed form open" (:193–194) | True as stated but superseded: at ¼ the ratio peaks near K = 40–60 and falls to 65.0555 < 1/t | D3ce-3 |
| Numerics A5 (:224–231) | Trace to the md ("5.30943" from 5.3094251) | [md], [mine] |
| Numerics B3, B4 (:241–250) | Trace to the md; same caveats as (iii) | [md] |

## Status of the round-81 issues

| D3cb | status | note |
|---|---|---|
| 1 | fixed in substance | **Done:** the exact evaluator (bit-exact by my independent one at 72/72); every B4 end and mark exact; `lm_rms`'s error printed; the test pins 60 against 90 digits. **Residue:** excesses below ~4·10⁻¹⁶ are still differences of double ratios (D3ce-5). |
| 2 | partly | **Done:** the real `b4_job` on two cells, and fakes that return distinct values. My mutants of those lines are killed: ER01–04, BE01–03, AF02–05, BJ01–04, BJ06–11, LF12. **New:** 20 survivors, 11 behaviour-changing — A5's new code, and LF02 behind the stall-only cells (D3ce-2). |
| 3 | fixed | **Done:** the docstring is the code (instrumented copy bit-identical); the flat stop is live (10/12, after 1.5·10⁴–1.6·10⁵ iterations); `stall_reanchor`; (iii) reworded. |
| 4 | partly | **Done:** both floors are stated. **Residue:** the direct one as "linear" (D3ce-1); the anchored one predicted from the end's rounding (D3ce-4). |
| 5 | fixed | **Done:** the exact series to K = 40 and 120; onsets 19, 8, 6; the leading-order onset dropped. **Residue:** the K = 1 truncation (D3ce-6); a closed form is available (D3ce-3). |
| 6 | fixed | **Done:** per-row references. At (¾, 10) the second variant's reference is 0.999951765916 = best/6.807792·10⁻¹⁰ [mine]; a two-variant test was added. |
| 7 | fixed | **Done:** one reference below 2·10⁻⁹; the 0.1% test's leniency (2.6·10⁻³) is stated. |
| 8 | fixed | **Done:** all five items. G01–G04 (the guard) are killed. |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ce-1 | major | **"Linear in its rounding" fails a scaling test.** Where: D58 :73, Status :162–163; round 81's referee proposed "linear in ε/e" (D3cb-4) without measuring it.<br>**(a) Double against long double, same start [raw, exact].** The double's stall excess over the long double's is 1.1·10³, 2.8·10⁴, 9.6·10⁴, 1.5·10⁵, 1.7·10⁵, 1.8·10⁵, 2.9·10⁵, 5.0·10⁵, 8.5·10⁵, 3.3·10⁷ (10 rows; at the other 2 the long double ends at the best end). The median is 1.75·10⁵ against 2048 = ε_d/ε_ld; the exponent is 0.91–2.27, median 1.58.<br>**(b) Controlled [mine].** I scaled the rounding: r_κ = r_ld + κ(r_d − r_ld), with B4's LM and parameters in double. From the stored optima at (1, 12), (¼, 6), (¾, 11) the stall excess grows ×159, ×216, ×178 from κ = 1/16 to 1 (exponents 1.83, 1.94, 1.87). Geometric means over five starts give exponents 1.41, 1.17, 1.62, 1.48; the starts scatter by up to 10⁶ at one κ (1.7·10⁻¹² … 2.1·10⁻⁶).<br>**(c) The anchored solve, same scan.** Its floor goes as κ² (slopes 2.00, 2.05).<br>**So** both floors grow superlinearly. They differ in size, 10³–10⁶ at one rounding, not in order. | Withdraw "linear in its rounding" at :73 and :162–163. State what was measured: a floor of the residual's evaluation, growing like κ^1.2–1.9 and path-dependent. Optionally add the κ-scan to B4 (an `lm_fit` option; seconds per start). |
| D3ce-2 | major | **20 of 87 mutants survive both checks; 11 change results on a rerun.**<br>**Totals.** The tests kill 61; `--from-raw` kills 6 more (BL01, BL06, SL01, SL03, SL05, SL06). The control passes both.<br>**Survivors that change results (probed against the raw):**<br>• `a5_far` is untested. FA01 (S(K−2)/S(K)), FA02 (no square root) and FA05 (q from σ*/σ) rewrite md :123, e.g. 2.333302 … → 5.444296 …, and the limit 2.333334 → 5.444446. FA03 (q at 110 digits) and FA04 (S(K−1) at 110 digits) change its low defects: 3.7·10⁻²² → 1.2·10⁻²¹, and → 1.0.<br>• The half-squeeze frame is untested. AS03 (wrong t) and AS04 (v′ = 2.001) change 9 and 4 lines of A5's series table.<br>• FS01 (the 10-term guard dropped) turns md :117 1.9·10⁻¹⁰ into 2.7·10⁻⁴ and :119 1.5·10⁻¹¹ into 6.7·10⁻³. FS02 (absolute low defects) turns 1.5·10⁻⁵² into 1.5·10⁻⁵⁰ and 3.7·10⁻²² into 4.4·10⁹. GN04 (weights from the iterate before the last Newton step) gives 1.1·10⁻⁵² and 7.0·10⁻²².<br>• LF02 (the flat test's reference never updated): at (¾, 11) and (¼, 6) from the optima the anchored solve runs to the cap, 400000/71898 and 400000/72909, instead of stopping flat at 14936/4000 and 22046/4000. 13 fields differ.<br>**Raw-only changes.** BJ05 (marks by `lm_rms`: 1.0854221631258343 → …1011765, md unchanged), AS05, FA06, FS04.<br>**Harmless on the data.** ER05 (30 digits: 0/72 ratios change), LF01, AF01, BB02, SL02.<br>**Cause.** `test_part_a5_on_its_own_numbers` runs `a5_far` and both frames but asserts only Mercer-frame values. The pinned B4 cells stop at a stall. | Pin `part_A5` against the raw (6 s), every field of `series` and `far`. Add cell 9 or 11 (6–7 s, flat after one periodic re-anchor) to `test_b4_job_on_two_cells_against_the_raw`. Add a K = 1, ρ = t test of `frame_series`. |
| D3ce-3 | minor | **A5's "its rate's closed form open" (:194), and τ/σ = ¼ (:193).**<br>**Conjecture [mine].** `scaled_rule` with variance v falls at r(v) = max{h(τ²/v), h(vσ*²/(σ²τ²))} per kernel, h(y) = (y − 1)/(y + 1). In a Mehler frame (ρ, v′) this is r = ρv′/(2 − ρ\|v′ − 2\|). The two arguments multiply to σ*²/σ², so r ≥ t, with equality only at v_hs (a cusp).<br>**What it reproduces.** All four of (ii)'s rates: λ, t, λ/(2 − λ), q. For the Mercer frame's rule, 1/r = (β⁴ + β² + 2)/(β²(β² − 1)) = 65.0555, 17.16497, 5.309401 and 7/3.<br>**Exact L² errors to K = 120–160** (limit by L + c/K² + d/K³):<br>• Mercer frame: 17.164987, 5.30940135, 2.33333332.<br>• At ¼: 65.1135, 65.1118, 65.0926, 65.0788, 65.0705, 65.0628 for K = 40 … 160 — it peaks, then falls toward the formula's 65.0555 (fits 65.054–65.057), below 1/t = 65.985.<br>• New points: 2.428608 (17/7) and 1.999999997 (2) at τ/σ = 2; 3.0000045 (3) below the Fock scaling and 1.666666667 (5/3) beyond moment matching at τ/σ = 1.<br>• Controls: λ to 10⁻⁹; the half-squeeze by K⁻¹, 2.618035982 against 1/t = 2.618033989.<br>**Derivation sketch.**<br>• Σ_n D_n(v′)sⁿ/n! = e^{(v′−1)s²/2}ε̃_K(√v′ s), with ε̃_K(ξ) = 1 − e^{−ξ²/2}Q_K[e^{ξz}].<br>• S_K = (πρv′)⁻¹∫\|ε̃_K\|²e^{−Ax²−By²}dA, with A = (1 − ρ(v′−1))/(ρv′), B = (1 + ρ(v′−1))/(ρv′).<br>• Numerically ε̃_K(ξ)·(2K)!/(K!ξ^{2K}e^{−ξ²/4}) = 0.68–1.28 on both axes (K = 25–100, \|ξ\| ≤ 1.6√K; K-th roots → 1).<br>• Laplace's method then gives r = 1/(2 min(A + ½, B − ½)), which is the formula. It also explains the convergence: K⁻² off the cusp, K⁻¹ at v_hs. | State the formula as a conjecture with this evidence. At ¼, replace "still rises" with the turnaround and the limit. "Slower at ½, 1 and 2" holds at every τ/σ and every v ≠ v_hs. Add ¼ to `a5_far` (precision ≈ 4K + 40 digits). |
| D3ce-4 | minor | **The anchored floor's accounting (:69–71, :164, :249).**<br>**(a) The frozen rounding is the last anchor's.** At 12/12 the last anchor is a stall re-anchor 1–28 iterations before the end, 9·10⁻¹¹ … 4.5·10⁻⁹ away at the flat stops. Its ‖P_Je₀‖²/(2‖r‖²) matches the gap (an exact Gauss–Newton step at the end) within ×0.40–2.0. The end's value (`anchor_floor`) matches only within ×0.23–7.7. Example, (¾, 11) from the block: 4.7·10⁻¹⁴ predicted at the end, 3.0·10⁻¹³ at the last anchor, 3.6·10⁻¹³ actual gap.<br>**(b) Part of e₀ is systematic.** 2σ² rounded to double (relative 4.66·10⁻¹⁷) alone predicts 15–64% of the floor, and it is the same at every anchor.<br>**(c) The double parameters cost as much at the lowest optimum.** The best end rounded to double lies 1.54·10⁻¹² above it at (¼, 6), where the anchored ends are 1.22 and 1.28·10⁻¹² above, and 3.1·10⁻¹⁴ at (¾, 11). With an exact anchor (κ = 0) the anchored solve still stops 1.0·10⁻¹² and 1.8·10⁻¹³ above. | Return the last anchor from `lm_fit` and evaluate the floor there, per row beside the gap. Name the 2σ² term. Add the double-rounded best end to B4's summary. |
| D3ce-5 | minor | **B4's smallest differences are rounding (the residue of D3cb-1).** The excesses are differences of double ratios. Exactly [mine]:<br>• anchored: 9.5·10⁻¹⁷, 2.57·10⁻¹⁶, 1.56·10⁻¹⁶, 1.12·10⁻¹⁶ (md 2.2, 4.4, 2.2, 2.2·10⁻¹⁶);<br>• 'nostall': 2.40·10⁻¹⁶, 8.41·10⁻¹⁶ (md 4.4, 6.7·10⁻¹⁶);<br>• long double: 2.4·10⁻²², 1.10·10⁻¹⁶ (md 0.0, 2.2·10⁻¹⁶).<br>D58 :71 and :248 quote "2.2·10⁻¹⁶ … 1.3·10⁻¹²"; exactly it is 9.5·10⁻¹⁷ … 1.28·10⁻¹². | Keep the exact rms (a Decimal string, or hi + lo) and form differences exactly, or print "< 4·10⁻¹⁶". |
| D3ce-6 | minor | **`frame_series` stops on a zero term.** At K = 1 in the half-squeeze frame D_n = 0 for n ≡ 0 (mod 4). The test `terms[-1] < tol·S` therefore stops at n = 24, not at 10⁻⁴⁰ of the sum. S(1) comes out 2.9·10⁻¹¹ low at τ/σ = 2. The md's "1.5e-11 (K ≤ 10)" (:119) is this truncation; summed fully it is 2.2·10⁻¹⁶ at K = 2 and 1.4·10⁻¹⁴ for K ≥ 3. D58 does not quote it. | Stop only on a nonzero small term, and test K = 1 at ρ = t. |
| D3ce-7 | cosmetic | (a) The "6·10⁻¹⁴ … 3·10⁻⁵" of :188–189 and :228 are the trapezoid's limits on the ratios; the series itself is exact to ≤ 7·10⁻¹⁵.<br>(b) "As K⁻²": the successive-difference ratios are 0.450, 0.534 at τ/σ = 2 and 0.425, 0.516 at 1, against 0.463, 0.543; at ½ they are 0.406, 0.507, still pre-asymptotic. Md :123's "converges, as K⁻², to a rate below 1/t" is a fixed string, not a check.<br>(c) Say that the best end is the 'anchored ld' end. Its independent support: the two starts agree to ≤ 2.1·10⁻¹⁹, and an exact Gauss–Newton step gains ≤ 1.4·10⁻¹⁷.<br>(d) From the blocks, 'nostall' still re-anchored 10–75 times periodically.<br>(e) The low-defect metric bounds the cancellation for n < 2K, not the tail: with FA04 the far low is 1.0, yet the printed ratios are unchanged.<br>(f) `lm_fit`'s docstring omits that acceptance also needs pred > 0.<br>(g) "The same 1.6·10⁻⁹" (:74) refers to two blocks. | Reword as indicated. |

## Sources

- **Local and computed [memory].**
  - The GH weights by the Christoffel sum, w_k = 1/Σ_{j<K}He_j(z_k)²/j!. Verified: Σw = 1 to 10⁻¹²⁰.
  - The Gaussian inner products: (2π)^{1/2}ab(a² + b²)^{−1/2}e^{−(c₁−c₂)²/(2(a²+b²))}. Verified: they reproduce A5's series to 7·10⁻¹⁵.
  - Σ He_n(x)sⁿ/n! = e^{xs − s²/2}; the Fock norm Σn!|c_n|² = π⁻¹∫|F|²e^{−|z|²}dA.
  - Python decimal's exp is correctly rounded; Sterbenz's lemma for the hi + lo split.
- **Web [snippet], not fetched.**
  - [Karvonen–Oates–Girolami, arXiv:2004.12654](https://arxiv.org/pdf/2004.12654) and [Math. Comp. 90 (2021)](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/): scaled GH rules without the length-scale condition of Kuo–Woźniakowski 2012; numerical rates r ≈ 0.528 and 0.203.
  - [Karvonen–Särkkä, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532): exponential convergence at scaled GH nodes.
  - Neither snippet gives a closed form for a scaled rule's rate; D3ce-3's formula is my conjecture.
- **Not checked.** The papers' full texts (arxiv.org and ams.org are unreachable).

## Lesson

A scaling law needs a scaling experiment. "Linear in its rounding" entered the Status from a referee's estimate, while B4's own two arithmetics already disagreed: about 10⁵ between them where linearity predicts 2048. Scaling the rounding by κ costs seconds and separates the two floors by size, not by order. The same holds for references and predictions. A best end made by refining the solve it judges confirms that solve only once it is validated independently; here it is, to 10⁻¹⁷. And a floor computed at the end describes an anchor the solve never used. Finally, an "open" closed form can be one guess away: all five rates lie on one curve.
