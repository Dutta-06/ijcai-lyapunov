# Referee report D3bx, round 77: Proposition D‴.58 as revised after round 70

## Scope

- **Commit.** I reviewed 39a8820 ("D‴.58 round 70: the range ends end by end; the floor in long double").
  - I exported it with `git archive` into `W/repo`, where W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref77/D3bx_work`.
  - All my work stayed in W (`repo`, `out`, `scripts`, `mut`, `logs`).
- **Files read in full.**
  - `review/drafts/D58.md`, `experiments/exp_width_mismatch.py`, `results/chaos/width_mismatch.md` and `_raw.npz`, `tests/test_width_mismatch.py` (35 tests).
  - `review/REVIEW_round70_D3bq.md`, the "Round 70" entry of `review/REVIEW_LOG.md`, and the commit's word diff of D58.md.
- **Files read in part.** In `src/core/rich_transit.py`: `Bumps`, `polish_block`, `ascent_target`, `sliding_iss`, `_alpha_of`. Also `exp_transit_mass.setup` and `exp_cohesion_leak.load/save`.
- **Run settings.**
  - venv Python: numpy 2.5.4, scipy 1.18.1; long double is x87, ε = 1.08·10⁻¹⁹.
  - Every run used `nice -n 10`, OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1; pytest ran with `-q -p no:cacheprovider`.
  - For `--from-raw`, WM_OUT pointed at `W/out`, which held a copy of the raw.
  - I never ran more than two heavy processes at once. Mutants ran in slim trees under `W/mut` (symlinks into the export).
- **Nothing in the repository was modified.** The repository's md, raw, script, draft and tests have the same md5 as the export. Every process I started has exited.
- **Marking.** Numbers are marked (md), (raw) or (mine).

## Summary

**What checks out**
- **Reproduction.**
  - `--from-raw` regenerates the md byte for byte (md5 4d85821c…; the raw is untouched, 3e22203d…).
  - 35/35 tests pass in 9.3 s (mine).
  - These reruns with the script's own functions are bit-identical to the raw (mine):
    - all 12 B4 cells (`lm_fit` ratios to 16 digits, iterations, stops);
    - `b3_job` with the stored block at (½, grid +0.05) and (¼, 5·10⁴ once), all four columns;
    - `b2_job` at three first non-optimal blocks, all 14 fields including the new `cond_joint_x`, `smin_ratio_x` and `smin_ratio_eq`.
- **B3 now restarts what the path did.**
  - At all 22 rows the restarted "path" polish reproduces B's stored first polish: status, evaluations and acceptance (raw). Examples: 1/3635/r at (1, grid), 1/30240/a at (¼, 5·10⁴).
  - So `first_nonoptimal`'s n0 points at the right `polish_block` call, repolishes included.
  - `seed_rule` is `sliding_iss`'s first-polish rule (rich_transit.py:415–419): the same 10⁻² seed set from w₀, the same residual and signs, and "start or ascent target within step·σ".
  - The counts 7/5/5/2/1 and every quoted B3 number match md :172–207.
- **The long-double residual is extended end to end (mine).**
  - `np.exp` on long double is within 6.9·10⁻²⁰ (relative) of a 60-digit exp, and the long-double matmul accumulates in long double.
  - At a B4-like point the residual's error is 4·10⁻²⁰, against 1.1·10⁻¹⁶ in double. The casts of f and X_B are exact.
  - The tests protect this: mutants that move the exp or the matmul to double (LM18, LM19, RM05) are killed.
- **Numbers.**
  - Every new number in (iii), the Status and the Numerics traces to the md: B3 :172–207, B4 :214–229, B2 :165, A5 :103–106, entries :232–241.
  - The three referee-only computations are now marked as not reproduced.

**What does not (majors first)**
1. **D3bx-1 (major): the floor is not "the arithmetic's".**
   - A pure-double solve passes it: the same LM with its residual evaluated about a fixed anchor reaches B4's long-double end from both below-floor blocks, and from all 12 B4 starts, within 3·10⁻⁸ (mine).
   - At the double stall points the objective computed in double is 0.39% and 0.46% lower at the long-double end, about 10⁴ times its rounding noise.
   - So "no step decreasing the objective as double resolves it" is false. The rounding destroys the step, not the objective.
2. **D3bx-2 (major): "censored by the floor", 'neither' and "soft, crossing 0.1% near 10⁻¹⁰" do not follow.**
   - At ¾, K = 11 the B3 restart spent its whole budget.
   - B4's own double LM from that block ends within OPT_TOL of the stored optimum.
   - The soft-floor curve is non-monotone and depends on the start and the solver; from the stored optima the double excess is at most 2.3·10⁻⁵.
3. **D3bx-3 (major): the tests still miss behaviour-changing mutants of the round's new code.** 11 surviving mutants silently change B3 or B4 rows on a rerun.
4. **Minors and cosmetics.**
   - D3bx-4: `b3_binds` labels.
   - D3bx-5: the rerun guard and n0 bookkeeping.
   - D3bx-6: B3 covers the first polish only.
   - D3bx-7: the B2 claim "whether … depends on the units" contradicts the Schur bound.
   - D3bx-8: A5's "falls more slowly" is shown at two widths only.
   - D3bx-9: cosmetics.

## What I checked and how

- **Reproduction.** As listed above. I also compared all 22 B3 "path" restarts with B's stored polish logs (raw).
- **The floor (B4).** I used all 12 of B4's starts: the two below-floor blocks, their stored optima, and the eight above 10⁻¹⁰. At each start:
  - `lm_fit` in both precisions (bit-identical to the raw);
  - the long-double end rounded to double, evaluated in long double (the true ratio) and in double (the objective F);
  - F's rounding noise, from 50 random 10⁻¹⁵ perturbations;
  - a 2401-value damping scan (μ = 10⁻⁴⁰ … 10²⁰) at the double stall;
  - `lm_fit` restarted from its own end;
  - the same LM in pure double with an anchored residual: r₀ = f − Φ(c₀)w₀ computed once, then r = r₀ − Σₖ φ₀ₖ[(wₖ − w₀ₖ) + wₖ·expm1(Δₖ)], with Δₖ = (cₖ − c₀ₖ)(2x − cₖ − c₀ₖ)/(2σ²). I ran it as is, re-anchored once, and re-anchored every 2000 steps.
- **The floor, below-floor blocks only.** Compensated summation of the residual (TwoSum and Dekker's TwoProduct), identity scaling, μ₀ = 10⁻¹⁰, and `polish_block` (5·10⁵ evaluations, no gradient test) from the final block.
- **`b3_binds`.** Nine synthetic edge cases.
- **B2.** The Schur bound σ_min(joint) ≤ σ_min((I − PP⁺)J_c), at the 20 first blocks under 22 centre scalings each.
- **A5.** The first-term rates qv′/2 and the tail's onset K, from q and v′.
- **Mutation testing.**
  - 99 single-point mutants plus a no-op control. They cover `lm_fit`, `lm_rms`, `b4_job`, `b4_jobs`, `b4_lines`, `b2_job` (new fields), `b2_line` (new sentence), `a5_lines`/`part_A5`, `b3_lines` (new summary), `b3_binds`, `first_nonoptimal`, `path_polishes`, `b3_job`/`b3_restart`, `b3_jobs`, `seed_rule`, `main`, `part_B4` and B4_MAX.
  - Each ran against the test file (`-x`) and against `--from-raw` (the md byte for byte).
  - Survivors in code that only a rerun reaches were probed with their own job functions:
    - `b4_job` on five B4 starts;
    - `b4_jobs`, `b3_jobs` and `first_nonoptimal` on the raw;
    - `b3_job` at (½, grid) and (½, 2σ), and on the sixth variant at ½;
    - `main` and `part_B4` with a fake pool.
- **Sources.** One web search (Nielsen's damping rule) and the installed scipy source (TRF's scaling).

## Verdicts per part

| claim (D58.md line) | verdict | evidence |
|---|---|---|
| (i), (ii), (iv), (v), proof (:8–23, :75–124) | Unchanged by this commit; not re-derived | word diff |
| (iii) δ_t, gap shares, nodes; C ≥ 0 and its flanks (:24–30) | Correct | md A, A3 |
| (iii) entries: 1 of 107, 316, 0.21, 0.963, 105 at 1.17σ, 2 at 0.989–0.992 (:30–35) | Numbers correct; the event definition drops the 10⁻⁴ dual condition | md :232–241; D3bx-9 |
| (iii) ranges K ≤ 4, 5, 7, 9, 12, 12, the variants, optima down to 6.8·10⁻¹⁰ (:35–40) | Correct | md :126–163 |
| (iii) B3 end by end: 7/5/5/2/1, 1.10–1.64σ, 1.0017–1.0157, 1.012, 31769, 1.88σ, 0.99998–1.00022, 3643–2.8·10⁵ (:40–50) | Correct for the first polish; the labelling logic is fragile | 22/22 restarts = stored polish (raw); D3bx-4, -6 |
| (iii) below the floor: 1.233, 1.0235, 1.0467 (:52–55) | Numbers correct; 1.0467 is a budget stop | md :195; D3bx-2 |
| (iii) "The floor is the arithmetic's (B4) … no step decreasing the objective as double resolves it" (:55–57) | Refuted | D3bx-1 |
| (iii) soft floor, 0.1% near 10⁻¹⁰, the restarts' excess (:57–60) | Numbers correct; the inference is unsupported | md :207, :229; D3bx-2 |
| (iii) stored optimum 1.6·10⁻³ / ≤ 3.1·10⁻⁵ above the long-double end (:60–62) | Correct | md :229 |
| (iii) "Those two ranges are censored by the floor" (:62) | Not shown | D3bx-2 |
| (iii) first blocks in the optimum's basin: VP 20/20, the polish given evaluations, long double from the block = from the optimum (:62–64) | Correct | md :165, :216–227 |
| (iii) condition numbers 1.3·10⁶–5·10⁸ and 5.4·10³–3.8·10⁵; σ_min ratios in three conventions (:64–68) | Numbers correct (rerun identical); "whether" is wrong | D3bx-7 |
| (iii) 33/38, 1.001, 16.36–24.74, 116/116, 53, 63 (:69–74) | Correct | md :163, :165 |
| Status (iii) (:136–140) | Inherits D3bx-1 and D3bx-2 | — |
| Status, Mercer frame: 5.322…, 2.33…, 63.49/17.17 against 63.58/17.49, qv′/2 = 0.375 < 0.382 (:159–166) | Numbers correct (0.375 mine); "falls more slowly" shown at τ/σ = 1, 2 only | md :83–86, :103–106; D3bx-8 |
| Status: referee-only numbers (:128–129, :134–135, :149–150) | Marked as not reproduced | — |
| Numerics A5, B, B3, B4 (:194–214) | Numbers correct; the B4 sentence inherits D3bx-1 | md |

## Status of the round-70 issues

| D3bq | status | note |
|---|---|---|
| 1 | partly | **Now killed:** the round-70 mutant classes for the gradient-test plumbing (my BR01), the seed rule's target branch (SR01, SR02), A5's Mercer-zero regression (A5e), B2's diagnostics (B2a–g, NL01–04), and the extended-precision core (LM18, LM19, RM05). **Survivors:** 33 of 99 survive both checks, against 45 of 91 in round 70. 11 of them change published B3/B4 rows, against 25 behaviour-changing in round 70 (D3bx-3). |
| 2 | fixed | Recorded starts, one setting at a time, every end judged by the path's own rule. 22/22 "path" restarts reproduce the stored first polish (raw). Residue: D3bx-4, -5, -6. |
| 3 | partly | **Fixed:** B4 is added and reproducible (12/12 bit-identical, mine); the ¾ reference is noted; "whatever the budget" is dropped. **Not fixed:** "arithmetic" fails a pure-double control, and the soft floor belongs to one solver from one kind of start (D3bx-1, -2). |
| 4 | mostly fixed | Three conventions are printed, with the condition numbers in the polish's units. "Whether" contradicts the Schur bound (D3bx-7). |
| 5 | fixed in the md | D58.md drops the 10⁻⁴ condition (D3bx-9). |
| 6 | fixed | Reworded; the levelling is shown at τ/σ = 1, 2 (D3bx-8). |
| 7 | fixed | B3 restarts the first polish only, which has no improvement test. Its "path" column is the path's settings from the recorded start. |
| 8 | fixed | The three computations only referees made are marked. |

## Issues

| id | sev. | issue | fix |
|---|---|---|---|
| D3bx-1 | major | **The floor is not the arithmetic's.** Where it is claimed: D58.md:55–57, :137–139, :210–212; md :212 ("stall: no step the arithmetic resolves decreases the objective"); `lm_fit` docstring :448; tests/test_width_mismatch.py:782–783.<br>**A pure-double solve passes the floor.** From both below-floor blocks I ran the same LM in pure double, with the residual evaluated about a fixed anchor: r = r₀ − Σₖ φ₀ₖ[(wₖ − w₀ₖ) + wₖ expm1(Δₖ)].<br>• It ends at 1.0000144 and 0.9984863.<br>• Re-anchored once, it ends at 0.99996913 and 0.99837354. The long-double ends are 0.99996913 and 0.99837352 (mine; md :224–227).<br>• At all 12 B4 starts, the re-anchored pure-double solve ends within 3·10⁻⁸ of the long-double end (mine).<br>**Other double variants (mine).** Compensated summation alone (TwoSum and TwoProduct, keeping exp's rounding) gives 1.00052 and 0.99881. Identity scaling and μ₀ = 10⁻¹⁰ do not help (1.0021, 1.0007).<br>**The double objective resolves the gap (mine).**<br>• At the two double stall points, F computed in double is 0.39% and 0.46% lower at the long-double end rounded to double. The rounding costs only 2·10⁻¹⁰ of the ratio.<br>• F's rounding noise there is 4.4·10⁻⁷ and 1.0·10⁻⁷.<br>• A 2401-value damping scan finds decreases of at most 1.5·10⁻⁶. So the step, computed from a gradient lost in rounding, stalls; the objective does not.<br>**B4's own table already shows it.** The double LM from the stored optima ends within 7·10⁻⁷ and 2.3·10⁻⁵ of the long-double end (md :225, :227). D58.md does not quote these rows. | • Replace "the arithmetic's" with "the floor of the residual as these solvers evaluate it (f − Φ(c)w by direct subtraction in double)".<br>• Delete "no step … as double resolves it" from D58.md:56, md :212 and both docstrings.<br>• Add a pure-double control to B4 (an anchored or compensated residual).<br>• Add the optimum-start rows to B4's summary. |
| D3bx-2 | major | **"Censored by the floor", 'neither' and "soft, crossing 0.1% near 10⁻¹⁰" do not follow.** Where: D58.md:52–62, :137–139; md :195, :205–207, :229; `b3_binds` :393.<br>**(a) The budget binds at ¾, K = 11.** The B3 restart with both settings changed spent its 5·10⁵ evaluations (0/500000, 1.0467: md :195). `b3_binds` ignores the stop's status and prints 'neither', which (iii) reads as the floor.<br>**(b) That block passes the experiment's own yardstick.** B4's double LM from it ends at 1.000700 times the stored optimum (md :226), which is within OPT_TOL. The stored optimum itself is 1.6·10⁻³ above the long-double end.<br>**(c) TRF does not measure the floor.** The polish's own TRF from the same final blocks stops on its step test at 1.0147 (3/256857) and 1.0037 (3/384829) (mine). That is 4–8 times LM's double excess, so the B3 ends (2.3%, 4.7%) are not B4's floor.<br>**(d) The curve has only six points.**<br>• It is not monotone: 1.9·10⁻⁸ and 1.2·10⁻⁶ at the same 1.6·10⁻⁹; 1.9·10⁻³ at 2.55·10⁻¹¹, below the 2.3·10⁻³ at 7.75·10⁻¹¹.<br>• It depends on the start: from the stored optima the double excess is at most 2.3·10⁻⁵ everywhere (md :217–227).<br>• The "crossing" interpolates a single gap, from 5.6·10⁻⁵ to 2.3·10⁻³.<br>• It vanishes with an anchored residual (D3bx-1).<br>**(e) "A single restart's excess … grows" (:59–60)** mixes in the budget-limited 4.7%. It also counts −2.1·10⁻⁵ at 6.8·10⁻¹⁰ (md :207), which is the stored optimum's own error. | • Describe the two ends as "not reached by the polish's TRF": its step test at ¼, its budget at ¾.<br>• Make `b3_binds` print the stop's status ("neither, budget spent").<br>• State the soft floor as one solver's, from one kind of start.<br>• Measure excesses against the long-double end, not the stored optimum. |
| D3bx-3 | major | **The tests miss behaviour-changing mutants of the new code.**<br>**Totals.** 99 mutants plus a control (the control survives, as it should). The tests kill 58. `--from-raw` kills 8 more (BL01, BL03, BL06, A5f, L3d, L3f, L3h, T04). 33 survive both.<br>**11 survivors silently change what a rerun publishes (mine):**<br>• BR03, acceptance step fixed at 1σ: (½, 2σ) flips from 1.51σ a to r, and "the budget" becomes "the budget (and the radius)". This changes the published 7/5 split.<br>• BR02, `seed_rule` given the end's amplitudes: (½, grid) flips from 1.10σ r to 0.00σ a.<br>• BR07, the gradient test forced to 10⁻¹⁵ in B3: on the sixth variant at ½, "budget" becomes 1/250115, 1.02256, labelled "the budget and the gradient test".<br>• PP04, `path_polishes`' real solve given 10⁻¹⁵: the rerun stops at a different polish with the same kernel count (start 2.677 instead of 445). The size guard lets it through.<br>• FN01, n0 += 1 per block: the test's synthetic blocks give n0 = 2 either way.<br>• BR08, the block's error left undivided: 14.19 becomes 8.9·10⁻⁷. The test's e_opt = 1 hides it.<br>• BS06: B4's start taken from the first path of that width.<br>• BS07: B4's start chosen by a 10⁻⁹ threshold (12 jobs become 8).<br>• PB01: `part_B4` runs half its jobs.<br>• MN01: `main` passes no B3, so a rerun's B4 is empty.<br>• LM15, damping that never relaxes: the (¼, K = 5) block start hits the 4·10⁵ cap at 1.087405 in both precisions.<br>**The other 22 survivors:**<br>• 7 change only B4's printed iteration counts (LM02, LM05, LM08, LM11–14).<br>• 5 relabel synthetic edge cases (T01, T05–07, T09).<br>• PP03 is caught at run time by the size guard.<br>• 9 are equivalent on these data. | • Test `b3_restart`'s acceptance with a move between 1σ and 2σ in a 2σ variant.<br>• Test it with a near-seed kernel whose end amplitude differs from its start.<br>• Run `b3_job` on a variant without the test, and assert the "path" and "budget" gtol.<br>• Give `first_nonoptimal` blocks with unequal polish counts; use e_opt ≠ 1.<br>• Test `b4_jobs` with two variants at one width and an optimum in [10⁻⁹, 2·10⁻⁹).<br>• Test `main --only B4` and `part_B4` with a fake pool.<br>• Run one long `lm_fit` descent from a far start. |
| D3bx-4 | minor | **`b3_binds` mislabels in general** (:385–395; mine, synthetic cases).<br>• "(and the radius)" is appended whenever any reaching end is rejected. With the budget-alone end accepted and the both-changed end rejected, the label is "the budget (and the radius)".<br>• "The budget or the gradient test (and the radius)" is printed when the test-alone end is accepted.<br>• "Nothing" is returned for a non-optimal block whose own end is accepted at the optimum. That is an inconsistency, not a label.<br>• "Reach" is one-sided, against a reference that is 1.6·10⁻³ too high at ¾, K = 11.<br>None of this is triggered above 10⁻¹⁰ on the data. | • Label from the accepted, reaching ends only.<br>• Flag "nothing" as an inconsistency.<br>• Below 2·10⁻⁹, use the long-double end as the reference. |
| D3bx-5 | minor | **The rerun guard and n0 bookkeeping are weak.**<br>• The rerun is checked only by its kernel count (:331). PP04 shows a rerun that stops at a different polish of the same size passes silently.<br>• The "path" restart is never compared with B's stored first polish (it agrees at all 22 rows: raw).<br>• `first_nonoptimal` counts polishes over the stored blocks, and `b_job` drops blocks with k = 0 (:187–188). | • Assert that the "path" restart's status and evaluations equal the stored first polish.<br>• Store R['polish_log'] per iteration, including empty blocks. |
| D3bx-6 | minor | **B3 restarts the first polish only.**<br>• The four radius-bound paths made 22–24 rejected polishes (½ base, ¾ base, ½ grid, 1 grid: md :132, :134, :138, :146). Their later polishes start from flowed states.<br>• So "the budget and the acceptance radius in 5" (:43–45) is shown for the first polish only.<br>• "All 24 polishes there end so" (:49) says more than the md shows: status 1 and rejected, with their end errors not logged. | Say "its first polish". Restart the last polish too, or log every polish's end error. |
| D3bx-7 | minor | **B2's "whether its slowest modes are slower … depends on the units" is wrong** (:66–68; md :165).<br>• In any common scaling of the centres, σ_min(joint) ≤ σ_min((I − PP⁺)J_c) (Schur complement).<br>• I find the ratio ≥ 1.000000 at all 20 first blocks, under σ, x and 20 random diagonal scalings, for both Kaufman's and Golub–Pereyra's Jacobian (mine).<br>• Only how much slower depends on the units. The unit-norm convention normalises the two problems' centre columns differently, so it is not a common scaling.<br>• "The polish's own units" omits TRF's Coleman–Li rescaling d = v^{1/2} of the bounded amplitudes. | • Write "never faster in any common scaling of the centres; slower by …".<br>• Mention TRF's rescaling. |
| D3bx-8 | minor | **A5's "falls more slowly" (:160–164) is visible at τ/σ = 1 and 2 only, on three K.**<br>• At τ/σ = 2 the ratio is still rising: 2.330, 2.331, 2.332 (md :106).<br>• The first terms nearly coincide at small τ/σ. 1/(qv′/2) = 65.998, 17.98, 5.887, 2.667 against 1/t = 65.985, 17.944, 5.828, 2.618 (mine).<br>• The tail overtakes the first term only from K ≈ 2/(q(2 − v′)) = 70, 22, 10, 8 (mine). That lies beyond the K the trapezoid allows at τ/σ ≤ ½, and only the onset is reached at 1 and 2.<br>• The first-term argument explains why the first term does not set the rate. It gives no rate. | Compute the Mercer frame's series exactly, or in extended precision, to K ≈ 40 in the script. State the widths where the slower rate is shown. |
| D3bx-9 | cosmetic | • D58.md:31–32 drops the md's condition that an event's duals agree within 10⁻⁴ (md :232).<br>• The A5 cell's .4g format prints "63.1" and "2.33" (:637).<br>• Variable projection's condition number is printed twice, though it is the same in both units by construction (md :165).<br>• B4's six-decimal ratios sit at the rounding level of e_opt itself: the optimum start's e_start is 1.0000006 at 2.55·10⁻¹¹ (raw). | • Add "with duals within 10⁻⁴".<br>• Use .4f.<br>• Drop the duplicate.<br>• Quote B4's ratios against the long-double end. |

## Sources

- **Nielsen's damping rule.** On acceptance: μ ← μ·max(1/3, 1 − (2ρ − 1)³) and ν ← 2. On rejection: μ ← μν and ν ← 2ν. This matches `lm_fit` :477–481. Confirmed from [GSL `multifit_nlinear/nielsen.c`](https://git.stg.centos.org/source-git/gsl/blob/2d5f50cfd2a1fdeaae06be7e1dbdbf009ee5897e/f/multifit_nlinear/nielsen.c) and [arXiv:1410.0793](https://arxiv.org/pdf/1410.0793) (web search). I did not read the 1999 DTU report itself.
- **TRF's working units.** Installed scipy 1.18.1, `scipy/optimize/_lsq/trf.py`: `v, dv = CL_scaling_vector(x, g, lb, ub)`, `g_norm = norm(g * v, ord=np.inf)` and `d = v**0.5 * scale`. These give the polish's bounds-scaled gradient test and its working units (local source).
- **From memory [memory]:**
  - Golub–Pereyra, SIAM J. Numer. Anal. 10 (1973): variable projection.
  - Kaufman, BIT 15 (1975): the projected Jacobian.
  - Moré, LNM 630 (1978): scaled Levenberg–Marquardt.
  - Dekker, Numer. Math. 18 (1971), and Knuth, TAOCP vol. 2: the error-free TwoProduct and TwoSum.
  - Ogita–Rump–Oishi, SIAM J. Sci. Comput. 26 (2005): accurate sums.
- **Not re-checked.** The draft's literature claims are unchanged in this commit.

## Lesson

Changing only the precision shows what causes a stall; it does not show that the arithmetic cannot do better. The long-double run proved that rounding stops this Levenberg–Marquardt solve. It could not prove that no double-precision method passes the same point, and B4's own optimum-start rows already hinted that one does. The decisive control for "the arithmetic's floor" is the best method in that arithmetic, not the same method in a wider one. Here a pure-double residual evaluated about an anchor removes the floor. Likewise, a verdict such as "neither" must carry the solver's stop: a run that spent its budget has not met a floor.
