# Referee report D3cb, round 81: Proposition D‴.58 as revised after round 77

## Scope

- **Commit.** dd88e96 ("D‴.58 round 77: the floor the residual's evaluation, not the arithmetic's"), branch claude/festive-bardeen-eopej2.
  - Exported with `git archive` into `W/repo`, where W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref81/D3cb_work`.
  - All my scripts, outputs, logs and mutants are under W (`scripts`, `out`, `logs`, `mut`).
- **Files read in full.**
  - `review/drafts/D58.md`, `experiments/exp_width_mismatch.py`, `results/chaos/width_mismatch.md`, `tests/test_width_mismatch.py` (41 tests).
  - `review/REVIEW_round77_D3bx.md` and the Round 77 entry of `review/REVIEW_LOG.md` (:3068–3107).
  - `git show dd88e96 --stat`, plus `git diff ed55ab1 dd88e96` of the script, D58.md (word diff) and the drafts README.
- **Files read in part.**
  - `src/core/rich_transit.py`: `Bumps`, `polish_block`, `ascent_target`, and the polish calls in `sliding_iss` (:400–445).
  - `exp_cohesion_leak.load/save`; `wide_feature` through its use.
- **Run settings.**
  - venv Python 3.13.16, numpy 2.5.4, scipy 1.18.1; long double is x87 (ε 1.08·10⁻¹⁹).
  - Every process ran under `nice -n 10`, OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1. pytest ran with `-q -p no:cacheprovider` (`-x` for mutants).
  - I never ran more than two heavy processes at once. WM_OUT pointed at `W/out`, which holds a copy of the raw.
  - Exact arithmetic used Python `decimal`: 60 digits for residuals, 110 for the A5 series. Long-double parameters were converted exactly as hi + lo doubles.
- **Nothing in the repository was modified.**
  - D58.md, the script, the md, the raw, the tests, `rich_transit.py`, `wide_feature.py`, `exp_cohesion_leak.py` and the round-77 report have the same md5 as the export (md 661e11ff…, raw 59efb734…). `git status` is clean.
  - During the review, other sessions' commits moved HEAD to 5201021: cd71612 (D‴.59 report) and 5201021 (D‴.52 round 78). They touch no D‴.58 file.
  - Every process I started has exited.
- **Marking.** (md) is the committed md, (raw) the committed raw, (mine) my runs.

## Summary

**What checks out**
- **Reproduction (mine).**
  - `--from-raw` regenerates the md byte for byte.
  - 41/41 tests pass in 10.9 s.
  - All 12 B4 cells, not a sample, rerun bit-identical to the raw in every field: e_start; double; long double; anchored (ratio, iterations, accepted steps, stop, marks); TRF (ratio, status, evaluations). Each took 1–489 s.
- **The anchored residual is an exact identity.**
  - φ(x; c) = φ(x; c₀)e^{Δ}, with Δ = (c − c₀)(2x − c − c₀)/(2σ²).
  - Hence wφ = φ₀[w₀ + (w − w₀) + w·expm1(Δ)] and f − Φ(c)w = r₀ − Σφ₀[(w − w₀) + w·expm1(Δ)].
  - `lm_fit` :498–506 implements it as written.
- **The central claim holds, and is stronger in exact arithmetic than published (mine, 60-digit Decimal at the stored ends).**
  - Every anchored end is at or below the long-double end from the block: −1.385·10⁻⁸ … +2.2·10⁻¹⁶.
  - The anchored ends from the block and from the optimum coincide to ≤ 2·10⁻¹³.
  - A long-double anchored solve started from them gains only 1.2·10⁻¹³ … 1.0·10⁻¹².
  - Re-anchoring every 500 instead of every 2000 accepted steps moves the end by ≤ 7.3·10⁻¹³. So re-anchoring does not bias the solve beyond ~10⁻¹².
- **B3.**
  - All 22 'path' restarts reproduce B's stored first polish in stop, evaluations and also acceptance (raw).
  - The 7/5/5/2/1 labels hold, as do 1.0235/1.0484 (md :209) and 1.5·10⁻²/5.3·10⁻³ (md :237).
  - Four of the six radius-bound blocks made 22, 23, 23 and 24 rejected polishes (raw).
- **B2.**
  - The bound σ_min(joint) ≤ σ_min((I − PP⁺)J_c) is stated correctly for Kaufman's form.
  - On the data it also holds for the Golub–Pereyra form actually computed: ratio ≥ 1.0000000110 (Kaufman) and ≥ 1.0000793 (GP) over 20 first blocks × 7 centre scalings (mine).
- **Numbers and wording.**
  - Every new or changed number in D58.md traces to the md.
  - The withdrawn phrasing is gone everywhere ("arithmetic resolves", "censored by the floor", "soft floor", the 0.1% crossing).
- **Round-77 survivors.** I re-created the ten behaviour-changing survivors (BR02, BR03, BR07, BR08, PP04, BS06, BS07, PB01, MN01, LM15) plus FN01. All are now killed by the tests (mine).

**What does not (majors first)**
1. **D3cb-1 (major): B4's numbers below ~10⁻¹⁰ measure the long-double rms, not the solves.**
   - `lm_rms` errs by −2.4·10⁻¹⁰ … +1.3·10⁻¹⁰ at the two blocks below 10⁻¹⁰.
   - It is biased low at all 12 long-double ends, because it is the very function that solve minimised.
   - In exact arithmetic the "+7·10⁻¹²" is −3.7·10⁻¹⁴, and "each end moving within 3.1·10⁻¹⁰" is ≤ 3.1·10⁻¹².
   - The md's mark values are noise, and 6 of the 7 capped runs had converged by 10⁵ iterations.
2. **D3cb-2 (major): 9 behaviour-changing mutants of B4's new code survive both checks.**
   - They sit in `b4_job` and in `lm_fit`'s anchored bookkeeping.
   - They are hidden by the test's fakes, and each silently changes published B4 rows on a rerun.
3. **Minors.**
   - D3cb-3: the anchored solve as described is not the one that ran. There were 7459 and 7657 re-anchors at stalls and none periodic. "To its own relative accuracy" is false along the slow modes, and the 'flat' stop is dead code.
   - D3cb-4: the accurate solves' own floors are not stated.
   - D3cb-5: A5's onset is off. The exact tail overtakes the first term from K = 6, 8, 19, so "pass only at τ/σ = 2" is false; the rate is computable.
   - D3cb-6: `b3_lines` keys its long-double reference by (τ/σ, K) across variants.
   - D3cb-7: two references for "the optimum" below 2·10⁻⁹.
4. **D3cb-8:** cosmetics.

## What I checked and how

- **Reproduction** (`scripts/b4_rerun.py`). `b4_jobs` on the raw, then `b4_job`'s calls with the end parameters kept, compared with `==` against the raw.
- **Exact evaluation** (`scripts/exact_ends.py`).
  - The rms of f − Φ(c)w over X_B, with every double taken exactly and exp in 60-digit Decimal (0.1 s per end).
  - Applied to all 48 B4 ends and 17 marks, and compared with `lm_rms`.
- **The anchored solve** (`scripts/lm_inst.py`, `inst_run.py`, `inst_analyse.py`).
  - An instrumented copy of `lm_fit`, bit-identical in end, iterations, accepted steps and stop to `WM.lm_fit` at cells 9 and 11.
  - It records every re-anchor (its kind, the period's computed decrease, the jump at re-anchoring), with sampled parameters evaluated exactly.
  - Variants: no stall re-anchoring ('nostall'); anchor_every = 500; the anchored residual in long double from the anchored double's end, as a reference.
  - The error of the change term against the exact change.
  - The predicted bias of a frozen anchor rounding e₀: ‖P_J e₀‖²/(2‖r‖²).
- **B3.** `b3_binds`/`b3_lines` on the raw. The reference mixing on the raw and on a synthetic two-variant case. The polish logs of the radius-bound blocks.
- **B2.** The Schur bound for Kaufman's and Golub–Pereyra's forms, under σ, x and 5 random diagonal centre scalings.
- **A5** (`scripts/mercer_series.py`).
  - The Mehler-frame series Σ_{n≥2K} ρⁿ/n!·D_n(v′)² in 110-digit Decimal.
  - Gauss–Hermite nodes by Newton from numpy's; weights K!/(K²He_{K−1}(z_k)²), with Σw = 1 to 10⁻⁶⁰.
  - Computed for the Mercer frame (ρ = q, v′ = 1 + q) and the half-squeeze frame (ρ = t, v′ = 2), K ≤ 40.
  - Validated against the md: it reproduces the md's A ratios and first-term shares, and A5's 5.3221, 5.3221, 5.3210 and 2.3304, 2.3314, 2.3320.
- **Mutation testing.**
  - 64 single-point mutants plus a no-op control (`scripts/mutants.py`), each in a slim tree under `W/mut`. Each ran against the test file (`-x`, 400 s timeout) and against `--from-raw`.
  - Ten re-created round-77 survivors.
  - Survivors were probed by running the mutant's `b4_job` on fast cells and comparing with the raw (`scripts/probe_b4.py`).

## Verdicts per part

| claim (D58.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii), (iv), (v), proof (:8–23, :89–138) | Unchanged by this commit; not re-derived | word diff |
| (iii) δ_t, nodes, C ≥ 0, flank peaks 2.0–2.6σ at K = 1 (:24–30) | Correct (sampled) | md :268–363: 2.02–2.57σ |
| (iii) entries, "its duals within 10⁻⁴" (:30–35) | Correct | md :240, :245 |
| (iii) ranges, variants, optima down to 6.8·10⁻¹⁰ (:35–41) | Unchanged; correct | md :126–165 |
| (iii) the restart reproduces the stored polish at all 22 (:41–43) | Correct | md :199; acceptance also agrees at 22/22 (raw) |
| (iii) what binds: 7/5/5/2/1; 1.10–1.64σ; 1.0017–1.0157; 1.012; 31769; 1.88σ (:43–50) | Correct | md :176–207; `b3_binds` on the raw (mine) |
| (iii) 24 polishes; 22–24 rejected at four of six (:50–52) | Correct | raw: 22, 23, 23, 24 rejected (mine) |
| (iii) the two ends below 10⁻¹⁰: 1.0235, 1.0484, step test or budget; TRF 1.5·10⁻² and 5.3·10⁻³ (:55–59) | Correct | md :209, :237; exact +1.477·10⁻², +5.287·10⁻³ (mine) |
| (iii) the double stalls 1.9·10⁻⁸ … 2.3·10⁻³; "1.0007 … passes the experiment's 0.1% test" (:59–62) | Numbers correct (exact agrees to 3 digits); the test's reference is 1.6·10⁻³ high | md :230, :233; D3cb-7 |
| (iii) anchored: "at the long-double end or below it … (−1.4·10⁻⁸ … +7·10⁻¹²…)"; "each end moving within 3.1·10⁻¹⁰" (:62–65) | Conclusion correct and stronger exactly; the +7·10⁻¹² and the 3.1·10⁻¹⁰ are evaluation noise; the solve is mis-described | D3cb-1, D3cb-3 |
| (iii) "the floor … the residual's … the solver's"; non-monotone; ≤ 2.3·10⁻⁵ (:65–67) | Correct; incomplete | md :233; D3cb-4 |
| (iii) optima converged only to the floor; B3 against the long-double end (:67–71) | Correct as quoted (the maxima); the md's lower ends mix references | md :211; D3cb-6 |
| (iii) "end at the polish's solver … not tested" (:71–72) | Fair | — |
| (iii) first blocks in the optimum's basin (:72–75) | Correct | md :167, :233 |
| (iii) B2: condition numbers; never faster in a common scaling (:75–82) | Correct; also holds for the GP form computed | ≥ 1.0000793 (mine) |
| (iii) 33/38, 1.001, 16.36–24.74, 116/116, 53/63 (:82–88) | Unchanged; correct | md :165, :167 |
| Status (iii) (:150–155) | Inherits D3cb-1, -3, -4 | — |
| Status A5: 5.3221…, 2.3304…, 63.4950, 17.1718; 0.375 (:176–180) | Correct at τ/σ = 1 and 2; at ¼ and ½ the last digits are the trapezoid floor's (exact 63.4971, 17.1721) | mine; D3cb-8 |
| Status A5: onset 70.0, 22.0, 10.2, 8.0; "pass only at τ/σ = 2"; "at the tail's onset"; "the rate … not shown" (:180–183) | Wrong | D3cb-5 |
| Numerics A5, B, B3, B4 (:212–233) | Numbers trace to the md; the B4 sentence inherits D3cb-1 and the A5 one D3cb-5 | md |

## Status of the round-77 issues

| D3bx | status | note |
|---|---|---|
| 1 | fixed in substance | **Done:** the anchored double passes the floor from all 12 starts (exact −1.385·10⁻⁸ … +2.2·10⁻¹⁶, mine); "not the arithmetic's" (D58 :59, :153); "no step … as double resolves it" deleted; optimum-start rows added (md :233). **New:** D3cb-1, -3, -4. |
| 2 | fixed | **Done:** `b3_binds` carries the stop (md :195, :197); B4 has a TRF column; "censored" and the crossing are dropped. **Residue:** the reference is keyed across variants (D3cb-6). |
| 3 | partly | **Done:** the ten re-created round-77 survivors and FN01 are killed (mine). **New:** 9 behaviour-changing survivors in this round's B4 code (D3cb-2). |
| 4 | fixed | **Done:** labels come from accepted reaching ends; "(and the radius)" only when every reaching end is rejected; "nothing" flagged as an inconsistency; the stop carried. **Residue:** D3cb-6. |
| 5 | fixed | **Done:** `b3_job` refuses a mismatching 'path' restart (BG01–03 killed); md :199 has "22 of 22"; n0 counts every call (FN02 killed). **Residue:** the md check line is untested (L302/L303 survive). |
| 6 | fixed | "first polish" throughout; D58 :51–52. |
| 7 | fixed | Verified; also holds for the GP form (mine). |
| 8 | partly | **Done:** an onset line is printed (md :108). **Not done:** it is a leading-order estimate, the inference drawn from it is wrong, and the requested exact series to K ≈ 40 was not computed (D3cb-5). |
| 9 | fixed | **Done:** "its duals within 10⁻⁴"; .4f; VP's condition number printed once; the summaries quoted against the long-double end. **Residue:** the .4f digits at ¼ and ½ are floor-limited (D3cb-8). |

## Issues

| id | sev. | issue | fix |
|---|---|---|---|
| D3cb-1 | major | **B4 is quoted below its instrument's resolution.** Where: D58:64–65, :231–232; md :233–235; `lm_rms` :557–561.<br>**The instrument (mine; 60-digit Decimal at all 48 ends from my bit-identical reruns).**<br>• `lm_rms` errs by −2.4·10⁻¹⁰ … +1.3·10⁻¹⁰ at 2.55·10⁻¹¹ and 7.75·10⁻¹¹, and by ≤ 3.1·10⁻¹¹ elsewhere.<br>• It is low at all 12 long-double ends (−2.3·10⁻¹² … −2.4·10⁻¹⁰). It is the same long-double function that solve minimised, so the reference is selection-biased.<br>**Exact results against the md.**<br>• The anchored range is −1.385·10⁻⁸ … +2.2·10⁻¹⁶ (md: … +7.2·10⁻¹²). The +7.2·10⁻¹² row (1, K = 12, from the optimum) is exactly −3.7·10⁻¹⁴: it is the difference of two evaluation errors, +2.4·10⁻¹² and −4.9·10⁻¹².<br>• The capped ends move ≤ 3.1·10⁻¹² over the last 2·10⁵ iterations (md: within 3.1·10⁻¹⁰).<br>• At 2.55·10⁻¹¹ from the block the md's marks −6.44, −6.19, −6.28, −6.50·10⁻⁹ are exactly −6.580, −6.581, −6.579, −6.582·10⁻⁹.<br>• The ends from the two starts differ in the md by 2.9·10⁻¹⁰ (2.55·10⁻¹¹) and 1.4·10⁻¹¹ (6.8·10⁻¹⁰); exactly they coincide to ≤ 2·10⁻¹³.<br>• 6 of the 7 capped runs had converged by 10⁵ iterations, the seventh by 2·10⁵. At cells 9 and 11 the exact objective is constant to 2·10⁻¹² from iteration ~1.4·10³ and ~3·10³ respectively.<br>**So** the conclusions survive and strengthen, but the published numbers measure the long-double evaluation. | • Evaluate every B4 ratio, including the marks, exactly (Decimal at 50–60 digits, 0.1 s per end) or with a double-double residual.<br>• Quote no difference below that evaluation's error.<br>• Restate D58:64–65: "−1.39·10⁻⁸ … 0; ends from both starts agree to 2·10⁻¹³; converged by 10⁵–2·10⁵ iterations".<br>• Add a test pinning the evaluator against an exact one at a 2.55·10⁻¹¹ point. |
| D3cb-2 | major | **The tests miss 9 behaviour-changing mutants of the round's B4 code.**<br>**Totals.** 64 mutants plus a control (the control survives, as it should). The tests kill 44; `--from-raw` kills 4 more (BL02, L304, L307, FN01). 16 survive both.<br>**Survivors that change B4 on a rerun (probed against the raw, mine):**<br>• LA17, the period counted by all accepted steps (:539): cell 0 goes from stall to flat, 10564/8120 → 7975/7966; cell 5 from stall to flat, 130766/12038 → 16663/2000, and its mark is lost.<br>• LA22, the first anchor 10⁻¹² off: cell 1 goes 8011/481 → 190/20.<br>• BJ07, anchor_every = 500 in `b4_job`: cell 0 goes 10564/8120 → 10542/8119.<br>• BJ08, the stop stored as 'cap': stall → cap.<br>• BJ09, iterations and accepted steps swapped.<br>• BJ05, the TRF given s = −1: the ratio goes 1.0000000104 → 6.2·10⁸.<br>• BJ06, the TRF ratio from the start's amplitudes: cell 0 goes 1.00000013 → 4.5·10⁷.<br>• BJ04 and MK01, the marks truncated (by construction).<br>**Harmless survivors.** LA07, LA14, LA21, BL05 and BB07 are equivalent on the data. L302/L303 weaken md :199's check.<br>**Cause.** `test_b4_jobs_job_and_lines` (:880–890) uses a fake `lm_fit` that returns iterations 7, accepted 5 for every call and 'cap' for the anchored one, and a `polish_block` stub that returns the start's amplitudes and ignores the signs. The real `b4_job` never runs. The period is tested only with anchor_every ≤ 20 from the start (:1004–1009). | • Run the real `b4_job` on cell 1 (1 s) and cell 0 (16 s), and pin iterations, accepted steps, stop, ratios and TRF status/evaluations against the raw.<br>• Make the fakes return distinct values per call and moved amplitudes.<br>• Test the period with a stall re-anchor inside it.<br>• Test md :199's line with one mismatching row. |
| D3cb-3 | minor | **The anchored solve as described is not the solve run, and its 'flat' stop is dead.**<br>**(a) The mechanism (mine; instrumented copy bit-identical to `lm_fit`).** D58:63 says "r₀ … once … re-anchored every 2000 accepted steps". At cells 9 and 11 the solve re-anchored 7459 and 7657 times, all at stalls with μ reset to 10⁻²⁰, and never periodically. An anchor period holds a median of 9 accepted steps (at most 25 over the last 200). The computed decrease per period is 4.2·10⁻¹² and 4.6·10⁻¹³ of F; the jump at re-anchoring is 4.4·10⁻⁷ and 1.3·10⁻⁷. The solve thus works within one frozen-rounding frame.<br>**(b) "The change … to its own relative accuracy"** (D58:63; docstring :484–485) is false along the slow modes. There the change is a cancelling sum, and its error is 4.8·10⁻⁹–2.3·10⁻⁸ of \|change\| (~10⁸ε), about 0.5ε·Σ_k\|terms\|. For a generic 10⁻⁹ displacement it is 1.2–2.3·10⁻¹⁶ (mine).<br>**(c) The 'flat' stop never runs.** Stall re-anchors reset the count, so 0 of 12 runs stopped flat (raw), and the capped runs ran ~3·10⁵ iterations after converging.<br>**(d) The control that separates residual and restarts.** Without stall re-anchoring, from the stored optima the solve stops after 95 and 1481 iterations, +1.6·10⁻¹² and +9.8·10⁻⁹ above the default end and still below the long-double end (mine). So the anchored residual, not the restarts, passes the floor. | • Describe the stall re-anchoring and the μ reset in D58.<br>• Replace "once" and "to its own relative accuracy" with "an error proportional to the parameters' displacement from the anchor".<br>• Run the flat test at stall re-anchors too.<br>• Add the 'nostall' row to B4. |
| D3cb-4 | minor | **The accurate solves' own floors are not stated.** Where: D58:62–67, :151–153; the Status calls the long-double residual "evaluated accurately".<br>**Long double (exact).** The long-double LM by direct subtraction stalls 6.58·10⁻⁹, 1.39·10⁻⁸ (from the blocks) and 2.20·10⁻⁸ (from the 7.75·10⁻¹¹ optimum) above the anchored end.<br>**Anchored double.** Its own floor is the frozen anchor rounding entering quadratically: ‖P_J e₀‖²/(2‖r‖²) predicts 7.0·10⁻¹³ (2.55·10⁻¹¹) and 1.1·10⁻¹³ (7.75·10⁻¹¹). A long-double anchored solve from its end gains 8.6·10⁻¹³–1.0·10⁻¹² and 1.2–2.3·10⁻¹³ (mine).<br>**So** every solve has a floor set by its residual's evaluation: linear in ε/e for direct subtraction, quadratic for the anchored one. | • State both floors.<br>• Use the anchored (or exact) end as the reference. |
| D3cb-5 | minor | **A5: the onset and the inference drawn from it.** Where: D58:180–183, :215–216; md :108; `a5_lines` :749–753.<br>**The series (mine; exact, and it reproduces md :105–106 to 4 decimals).**<br>**(a) The formula.** It drops R₁'s lower orders: R₁ = (v′−2)K² + (v′−3)K + (v′−1). The second term exceeds the first from K = 6.13, 8.74, 20.7, 68.7, not 8.0, 10.2, 22.0, 70.0. At τ/σ = 2, K = 8 the ratio is 1.55, not 1.00.<br>**(b) The tail.** The whole tail overtakes the first term from K = 6 (τ/σ = 2), 8 (1), 19 (½), and beyond 40 (¼).<br>**(c) The inference.** "Pass only at τ/σ = 2" is false: at τ/σ = 1 the tail is 1.13, 1.46 and 1.86 times the first term at K = 8–10. At τ/σ = 2 it is 2.5–4.4 times, which is not "the onset".<br>**(d) The rate.** "Not shown", but it is computable: the ratio per kernel is 2.3330/2.3332/2.3333 at K = 20/30/40 (τ/σ = 2, 1/t = 2.618); 5.3115/5.3100/5.3097 (1, 1/t = 5.828); 17.203/17.184/17.174 (½, 1/t = 17.94); and at ¼ still rising (65.11 at K = 40). The slower fall persists at τ/σ ≥ ½. | • Compute the frame series exactly in the script (`gh_defect_v` or a high-precision recurrence) to K = 40.<br>• Print tail/first and the ratio per kernel.<br>• Drop the leading-order onset, and state the slower fall at τ/σ = ½, 1, 2. |
| D3cb-6 | minor | **`b3_lines` keys its long-double reference by (τ/σ, K)** (:421–422), while the stored optima differ between variants at one (τ/σ, K).<br>• At ¾, K = 10 they are 6.807626·10⁻¹⁰ and 6.807792·10⁻¹⁰ (raw), 2.4·10⁻⁵ apart.<br>• So md :211 prints "+2.8e-06" where the row is +2.72·10⁻⁵ above the long-double end, and the 1.61·10⁻⁹ lower end is +3.88·10⁻⁸, not +3.1·10⁻⁸ (mine). D58 quotes only the maxima, which are right.<br>• Latent: in a synthetic case with optima 0.5% apart, an end 0.05% above the long-double end is printed "1.0055 times B4's long-double end" and labelled "neither (the budget spent)" (mine). | Store B4's absolute end (ratio × its e_opt) and divide by each row's own e_opt. |
| D3cb-7 | minor | **Two references for "the optimum" below 2·10⁻⁹.**<br>• D58:61–62 says the double end at ¾, K = 11, "1.0007 times the stored optimum, still passes the experiment's 0.1% test".<br>• Against B3's own reference below 2·10⁻⁹ that end is +2.33·10⁻³ (md +2.3·10⁻³), and it fails `b3_binds`' test. The stored optimum itself is +1.63·10⁻³ (md :233). | Use one reference below 2·10⁻⁹, or say the test's reference is 1.6·10⁻³ high there. |
| D3cb-8 | cosmetic | • D58:179, :214–215 quote 63.4950 and 17.1718: digits at the trapezoid floor (exact 63.4971, 17.1721; mine).<br>• "From the paths' last blocks" (D58:58) means the path's block at that K, not the paths' last blocks.<br>• md :211's "a negative value is the end below the stored optimum" is the long-double end below 2·10⁻⁹.<br>• `b3_job`'s guard runs after the four restarts (two of 5·10⁵ evaluations) and does not compare acceptance.<br>• The log and the docstring of `test_lm_fit_flat_stop` claim "each step's decrease" is tested; it is not asserted. | Fix each wording; move the guard before the large restarts and include acceptance. |

## Sources

- **Local and computed.**
  - The Gauss–Hermite weights K!/(K²He_{K−1}(z_k)²) for N(0, 1) [memory]. Verified: Σw = 1 to 10⁻⁶⁰, and the md's half-squeeze ratios and first-term shares are reproduced (mine).
  - The Mehler-frame series and D_{2K} = v′^K K! come from D58 and `gh_defect_v`; checked by reproducing md :103–106.
  - The projection bound σ_min([A, B]) ≤ σ_min((I − AA⁺)B) is elementary: take x = −A⁺By.
- **From memory [memory].**
  - Python `decimal`'s exp is correctly rounded.
  - Golub–Pereyra, SIAM J. Numer. Anal. 10 (1973); Kaufman, BIT 15 (1975).
  - Rounding of a sum ∝ Σ\|terms\| (Higham, *Accuracy and Stability of Numerical Algorithms*).
- **Not checked.** No web search this round (arxiv.org and ams.org are unreachable). The draft's literature claims are unchanged in this commit.

## Lesson

A claim about a floor has to be measured by an instrument finer than that floor. The round's evidence is real, and in exact arithmetic it is stronger than published: the anchored solve passes the double floor, and its ends from two starts agree to 10⁻¹³ and stop moving long before the cap. Yet every number B4 prints below ~10⁻¹⁰ — the range's +7·10⁻¹², the marks, "within 3.1·10⁻¹⁰" — measures the long-double rms. That rms is the very function the long-double solve optimised, so it is biased toward that solve. Evaluating the objective exactly costs 0.1 s an end. The same discipline applies to describing the solve: say which mechanism actually ran. Here it was thousands of stall re-anchors and never a periodic one.
