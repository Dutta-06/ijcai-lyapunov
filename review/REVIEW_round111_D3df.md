# Referee report D3df, round 111: Proposition D‴.60 as revised after round 105

## Scope

- **Commit.** I reviewed 232ca83 on claude/vigilant-clarke-hf939b. I exported it with `git -C … archive 232ca83 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref111. HEAD moved while I worked (148026d, then 8139e12). `git diff 232ca83 8139e12` is empty for the draft, the script, the md, the raw, the tests and the mutant spec.
- **Read in full:**
  - review/drafts/D60.md (378 lines), especially (iii), the Status and the Numerics;
  - the C8R/C8D and C8A paragraphs of results/chaos/noise_resonance.md;
  - in experiments/exp_noise_resonance.py: `c8a_job`, `C8A_PROBE`, `C8A_SCAN`, `c8a_lines`, and the C8R clause of `c8_more_lines`. Also what they call: `state_start`, `offsets_start`, `win_run`, `edge_run`, `c_run`, `c8a_cases`, `c8d_alarms`, `left_periods`, `raised`, `DZ.fitted_state`;
  - in tests/test_noise_resonance.py: the round 98, 102 and 105 tests, and `test_round89_md_is_the_report…`;
  - tools/mutants/round105_D60_own.py and tools/mutants.py;
  - the Round 105 entry of REVIEW_LOG.md and review/REVIEW_round105_D3cz.md;
  - `git diff e68fcd3 232ca83` of these files, and the drafts README line.
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv (symlinked as W/repo/venv);
  - NR_OUT=W/out holding a copy of the raw, TMPDIR=W/tmp;
  - pytest `--basetemp` inside W (through `PYTEST_ADDOPTS` for tools/mutants.py);
  - at most two heavy processes of mine at a time, each under `timeout`.
- **Nothing in the repository was modified.** I ran only read-only git commands and killed no process.
- **Tags:** [md] the committed md or draft; [raw] the committed raw; [own] my runs or code; [snippet] a web-search result only; [memory] not checked this round.

## Summary

**Design.** The revision does four things:
- adds to C8A's 12 isolated-alarm windows, on the reduced base, kernel 1's centre alone at its period mean, at pT and at (p − 1)T, a six-point scan of that centre, and a "warm" start one period earlier. Two conclusions are printed from the counts: "the alarms are the cold restart's" and "kernel 1's centre at the restart is one coordinate against a per-window threshold";
- runs draw 5's probe over all 8 subsets and renames the noise clause;
- tightens the centre, 'needs' and C8R conditions;
- tests the rerun path.

**What holds**
- **Reproduction is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 9c2f943e…);
  - 80 tests pass (106 s);
  - all 7 C8A jobs, rerun with the script's own `c8a_cases` and `c8a_job`, equal raw['C8A'] in every field except `secs` (8.9–44.1 s each).
- **Every new number agrees with the raw:** 7, 10 and 11 removed; 11, 9, 8, 5, 0, 0 along the scan; 12 of 12 "monotone"; "warm" 0 of 12; the class ranges; the 8 probe rows.
- **The reduced base is right for every new start** [raw]:
  - the scan at offset 0 equals C8D's window at p, leave period and largest distance, to 12 digits at 12 of 12;
  - 'warm' equals C8D's window at p − 1, to 12 digits at 12 of 12.
- **At the 11 alarms at k_e = 20 both readings hold, on finer evidence than the six offsets** [own]:
  - On a 49-point grid (step 0.0025, −0.08 … +0.04) the outcome is monotone with one threshold in each window, 0.005–0.045 outward of the feature. The 22 instantaneous samples (pT and (p − 1)T) agree with it.
  - The reduced state started 1–5 periods early, with δ(p) or with its own δ, passes the alarm's period in 110 of 110 runs.
- **Draw 5:** kernel 1 leaves at 160 if and only if kernel 3's offset is applied (8 of 8 subsets).
- **Mutants:** the authors' 12 are all killed.

**Findings**
1. **D3df-1 (major).** Draw 7's alarm (k_e = 10, s = 0.02, from 107) breaks both new conclusions, and the code's conditions let them print:
   - no centre value from −0.08 to +0.04 removes it, so its "monotone" scan is constant;
   - 'warm' runs 10 periods from 106 and ends one period before the alarm's period 117;
   - the reduced state with the window's own δ(107), started 2–10 periods earlier, leaves at 117 in 9 of 9 runs.
2. **D3df-2 (minor).** 'Warm' is C8D's dense window at p − 1: a cold restart carrying δ(p − 1), already counted in C8D's "0 of 106 repeat it", and checked against the exact period only. At s = 0.06, draw 2 from 58 the warm pass is not kernel 1's centre.
3. **D3df-3 (minor).** The rerun-path test pins the order of (W, C) but not the period. Two mutants that move or flatten the period means on the rerun path survive the full suite and change the published C8A paragraph on regeneration.
4. **D3df-4 (minor).** Seven more survivors, equivalent on the data: the 'needs' set, the noise clause, the C8R boundary, ", not its own", the unsorted scan, the centre clause by leave period.
5. **D3df-5 (cosmetic).** Wording: "where it delayeds it", "+0: 0", "kernel 1's leave needs the coming noise" right after a noise-free leave of kernel 1, stale round citations, the base dependence unstated.
6. **D3df-6 (cosmetic).** Carry-overs: `kk = 8`, and draw 7's marginality, which D3df-1 now turns on.

## What I checked and how

All paths are under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref111.

- **Reproduction:**
  - W/logs/fromraw.log, then `cmp` against the committed md: identical.
  - W/logs/pytest_full.log: 80 passed, 1 old warning.
  - W/scripts/repro_c8a.py → W/logs/repro_a.log and W/logs/repro_b.log: 7 of 7 jobs identical (2 on the C3 path, 5 on the rerun path).
- **Recount (numpy only, `python -I`):** W/scripts/recount.py → W/logs/recount.log. Every new count; the scan at 0 against C8D(p); 'warm' against C8D(p − 1); δ(p) − δ(p − 1) = 0.00098–0.0022; the class ranges; the probe rows; the reruns.
- **Own probes.** Both use my own runner: `simulate`, the period means, and the first period beyond 0.3.
  - W/scripts/probe.py → W/logs/probe_a.jsonl and W/logs/probe_b.jsonl, analysed by W/scripts/analyse.py → W/logs/analyse.log. Per window:
    - the reduced base with kernel 1's centre on 49 offsets;
    - the base with the centre at the rerun's values at pT and (p − 1)T;
    - warm starts from p − 1 with δ(p − 1) and δ(p) (11 periods), and from p − 2 and p − 3 with δ(p), recording kernel 1's centre and coefficient at pT.
  - W/scripts/probe2.py → W/logs/probe2_a.jsonl (draw 7, k = 1…10) and W/logs/probe2_b.jsonl (k_e = 20, k = 1…5): the reduced state from p − k with δ(p) and with δ(p − k), lengthened to 10 + k periods.
- **Planted inputs:** W/scripts/synthetic.py → W/logs/synthetic.log.
- **Mutants:**
  - The authors' 12 plus 13 of my own (W/mut/spec111.py, combined in W/mut/allspec.py) on `-k "round105 or round102 or round98 or round89 or c8_more_lines or c8_windows"` → W/mut/w1/res_a.jsonl.
  - The 9 survivors on the full suite → W/mut/wfull/res_a.jsonl and W/mut/wfull/res_b.jsonl.
  - R01 and R02 regenerated (jobs 2–6, the rerun path) with W/scripts/regen_mut.py → W/mutdemo/R01/c8a_line.txt and W/mutdemo/R02/c8a_line.txt.
- **Literature:** two web searches.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii), EDGE, power law, C4, s*, C7, C7L, C7S, C8, C8N, C8R, C8D, C8K | unchanged since round 105, not reverified | [md] |
| Draw 5: kernel 1 leaves at 160 iff kernel 3's offset is applied, "all 8 subsets" (122–123) | verified | [raw], [own] |
| "Without it [the noise] … 1.2·10⁻⁶ … kernel 1's own ejects it in period 2"; kernel 3's direction "not recorded" (124–127) | numbers verified; the caveat accurate. The md's renamed clause reads ambiguously (D3df-5) | [raw], [md] |
| C8R clause "(implied by C8R's aggregate …)" (133–135) | condition now checks the aggregate with C8R's own "earlier" definition (lag < −1); true on the data | [raw], [own] |
| 'w' keeps 12; 'wc' 6/1/5; 'w+c1' = 'wc', 'w+c23' = 'w' by leave period; class ranges −0.0713 … −0.0325, −0.0250 … +0.0028, −0.0467 (136–142) | verified | [raw] |
| "That split is one coordinate's … 7 … 10 … 11 … 5, 8, 9, 11 … at 0 and +0.04 at none" (143–146) | numbers verified | [raw] |
| "The outcome monotone in the offset at every window, so kernel 1's centre … is one coordinate against a per-window threshold, and the period mean and the instantaneous value two samples of it" (146–148) | true at the 11 windows at k_e = 20 (fine grid; 22 of 22 samples agree). Vacuous at draw 7: no threshold, the coefficient's alarm | [own]; D3df-1 |
| "The reduced state started one period earlier leaves in the alarm's period at none of the 12: the alarms are the cold restart's" (148–149) | true and robust at the 11 (110 of 110 earlier starts pass). At draw 7 the warm window cannot reach the alarm's period, and with δ(107) earlier starts reproduce it (9 of 9): false there | [raw], [own]; D3df-1, D3df-2 |
| Numerics C8A (365–369); Status (iii); drafts README | as (iii); same overgeneralisation | [md] |
| Sources (Schultz et al. 2018; Feudel et al. 1998) | consistent with the snippets | [snippet] |

## Status of the round-105 issues

| D3cz id | status | what remains |
|---|---|---|
| D3cz-1 | fixed in part: the starts are in, the numbers right, the base right, and the k_e = 20 reading confirmed on a fine grid | the referee's text kept draw 7 apart ("at k_e = 20 all 11 …; draw 7's … a crossing of 0.310 that no centre move removes"). The revision generalises to 12 (D3df-1). The base dependence (draw 2 from 58) is unstated (D3df-5) |
| D3cz-2 | fixed in part: 8 subsets; the caveat on kernel 3's direction added | the new clause's wording (D3df-5) |
| D3cz-3 | fixed for the W/C swap | row-offset mutants on the rerun path survive and change the md on regeneration (D3df-3) |
| D3cz-4 | fixed: the C8R aggregate, the centre clause (leave periods, 'wc' ≠ 'w'), 'needs' on the four kernel-3-free subsets, planted inputs; the authors' 12 mutants killed | the round's own new conditions are vacuous or weak (D3df-1, D3df-2) |
| D3cz-5 | `c8a_job`'s docstring fixed | `kk = 8` (logged as not done); `c8a_lines`' docstring still says "rounds 98, 102" |
| D3cz-6 | the ranges by class; the Round 105 entry is in 232ca83 itself | draw 7's marginality still unstated (D3df-6); the generated "delayeds" (D3df-5) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3df-1 | major | **Draw 7's alarm is neither one coordinate's nor the cold restart's, and both conclusions print for all 12** (D60:143–149, 365–369; md C8A; drafts README; `c8a_lines`).<br>• **The scan is constant there.** [raw]: the scan leaves at 10 at all six offsets (0.3103–0.3104), as do 'r+c1m', 'r+c1i' and 'r+c1p'. [own]: 10 at all 49 offsets from −0.08 to +0.04 (0.31034–0.31038). `mono` counts a constant sequence as monotone, so "12 of 12" includes a window with no threshold.<br>• **The warm window cannot see the alarm.** [raw]: the alarm is in the window's period 10, absolute 117. 'Warm' runs NP_C7 = 10 periods from 106 and ends at 116, so `warm[0] == reduced + 1 = 11` cannot happen.<br>• **Earlier starts with the window's coefficient reproduce it.** [own, probe2_a]:<br>  – lengthened to 11 periods, the start from 106 passes (0.0204 at 117), with δ(106) and with δ(107);<br>  – the reduced state with δ(107) = 0.07667 started at 105 … 97 leaves exactly at 117 in 9 of 9 runs (0.3093–0.3101);<br>  – with its own δ(q) = 0.0683–0.0757 it passes in 10 of 10 (≤ 0.0205).<br>• So this alarm is δ(107) against a marginal crossing of the cut, the "its coefficient 1" of the split, not the restart.<br>• Planted [synthetic A]: an alarm in the window's last period with a constant scan prints both conclusions. Planted [synthetic C]: a scan removed at every offset prints "one coordinate". | **D60 (iii)**, replacing "That split is one coordinate's … the alarms are the cold restart's" with:<br>"At the 11 alarms at k_e = 20 the split is one coordinate's: from the reduced start each is removed once kernel 1's centre lies 0.005–0.045 outward of its feature (monotone on a grid of step 0.0025). The period mean is outward enough at 7, the run's value at pT at 10, at (p − 1)T at 11. The reduced state started 1–5 periods earlier, with δ(p) or its own δ, passes all 11 alarm periods: these alarms are the cold restart's. Draw 7's at k_e = 10, a crossing of 0.310, is neither: no centre from −0.08 to +0.04 removes it, and the reduced state with δ(107) started 2–10 periods earlier leaves at 117 too."<br>**Same fix** in the Numerics and the README line.<br>**`c8a_lines`:**<br>• count as "one coordinate" only monotone scans that contain both outcomes and whose offset 0 reproduces `reduced`; print that count;<br>• print "cold restart's" only for windows whose warm run reaches the alarm's period and does not leave within ±1 of it.<br>**`c8a_job`:**<br>• 'warm' with NP = NP_C7 + 1, plus a δ(p) twin (and from p − 2, p − 3), about 1 s each. |
| D3df-2 | minor | **'Warm' is C8D's window at p − 1, and its condition is weak** (`c8a_job`, `c8a_lines`; D60:148).<br>• [raw]: identical to C8D(p − 1) to 12 digits at 12 of 12. It is a cold restart with δ(p − 1), not a warm state, and at the 11 alarms with reduced = 2 it is among the 106 "other dense windows that cover an isolated alarm's leave period" that C8D already reports (0 repeat).<br>• The check is the exact period, while C8D's "repeat" uses ±1. [synthetic B]: a warm leave one period after the alarm's period still prints "the cold restart's".<br>• What holds [own]: at 10 of the 11, the warm loop's kernel 1 centre at pT is outward of the window's threshold (−0.040 … −0.098, within 0.001 of the run's value), so warm and cold differ by that coordinate.<br>• What does not: at s = 0.06, draw 2 from 58 it is −0.0257, inward of the threshold in (−0.045, −0.0425]. A cold restart with the run's pT centre (−0.0267) leaves there, yet the warm start passes, so there the cold/warm difference is not kernel 1's centre. | Call it "C8D's window at p − 1, lengthened" or cite C8D's repeat count. Use ±1 as in `c8d_alarms`. Add the δ(p) warm start. State the draw 2 from 58 exception in D60 (iii). |
| D3df-3 | minor | **The rerun-path test pins the order, not the period** (`test_round105_…`: Wp and Cp have identical rows, and only C[0] is compared).<br>Two mutants survive the full suite (80 passed):<br>• R01: `np.roll(Wp, −1, 0)`, `np.roll(Cp, −1, 0)` on the rerun path, i.e. the period means one period late;<br>• R02: every period's centre mean replaced by the first period's.<br>Regenerating jobs 2–6 under them changes the md's C8A paragraph [W/mutdemo]:<br>• R01: 'wc' "kept 7, removed 5" with no delayed alarm; 'r+c1m' "removed at 4 of 12"; "the 7 that 'wc' leaves"; ranges −0.0465 … −0.0122;<br>• R02: 'wc' "kept 11, removed 1"; 'r+c1m' "0 of 12". | Distinct rows in the mock (e.g. `Wp = arange`), asserting that the starts receive row p − 1. Better: rerun one rerun-path job against the raw (job 3, 15 s), as test_round102 does for draw 7. Add raw invariants: scan at 0 equals C8D(p) and 'warm' equals C8D(p − 1) at 12 of 12. |
| D3df-4 | minor | **Survivors equivalent on the data** (full suite):<br>• R03: the scan unsorted (C8A_SCAN happens to be ascending);<br>• R04, R05: the 'needs' set without (0, 1) or without ();<br>• R06: the noise clause without the noisy kernel-3 row;<br>• R07: the C8R aggregate's −1 → −2 (a loop two periods early would no longer break it);<br>• R11: ", not its own" unchecked;<br>• R13: the centre clause's 'wc' ≠ 'w' by leave period.<br>Among mine, R08, R09, R10 and R12 are killed. | Planted rows for R04 (subsets lacking (0, 1)), R07 (a loop at P − 2), R06 and R11. Record R03 and R13 as equivalent. |
| D3df-5 | cosmetic | **Wording.**<br>• md: "−0.0467 where it delayeds it" (`k_ + 's'`; "earliers" likewise); "+0: 0".<br>• "— kernel 1's leave needs the coming noise" follows directly "kernel 1's offset: kernel 1 leaves in period 2" (without noise) and reads as contradicting it. The quiet rows show the leave needs noise; that it needs this noise was the round-105 referee's (shifted noise, other draws), not recorded.<br>• The C8A heading and the probe sentence still cite "(referee rounds 98, 102)".<br>• D60:142–143 applies the scan, run on the reduced base, to a split made on the 'wc' base. Draw 2 from 58 ("kernels 2 and 3's state 1") is removed by kernel 1's mean centre on the reduced base and only delayed on the 'w' base. | Verb table ("delays", "keeps"); "0" for "+0"; "— kernel 1's leave with kernel 3's offset needs the noise (other noise not recorded)"; "(referee rounds 98, 102, 105)"; one clause on the base dependence. |
| D3df-6 | cosmetic | **Carry-overs.** `kk = 8` is still separate from `c_run`'s k (D3cz-5). Draw 7's marginality (0.310 against 0.3) is still unstated (D3cz-6), and D3df-1 turns on it. | Return k from `c_run`; state the margin in (iii). |

## Sources

- P. Schultz, F. Hellmann, K. N. Webster, J. Kurths, "Bounding the first exit from the basin: independence times and finite-time basin stability", Chaos 28, 043102 (2018) [snippet]. Finite-time basin stability is the probability that a perturbed trajectory returns within a given time. This agrees with D60's "the windowed counterpart of C7's sets".
  - [PIK record](https://publications.pik-potsdam.de/pubman/item/item_22275)
  - [Aberdeen record](https://abdn.elsevierpure.com/en/publications/bounding-the-first-exit-from-the-basin-independence-times-and-fin)
- U. Feudel, A. Witt, Y.-C. Lai, C. Grebogi, "Basin bifurcation in quasiperiodically forced systems", Phys. Rev. E 58, 3060–3066 (1998) [snippet]. The records describe fractal and Wada basin boundaries, and isolated "islands" of basins created by a basin-boundary bifurcation. This agrees with D60's "basins with islands and fractal or Wada boundaries".
  - [Aberdeen record](https://abdn.elsevierpure.com/en/publications/basin-bifurcation-in-quasiperiodically-forced-systems/)
  - [IREAP record](https://ireap.umd.edu/node/3175)
- Not rechecked this round [memory or earlier rounds' snippets]: Giesl–Rasmussen, Oljača–Ashwin–Rasmussen, Berger–Doan–Siegmund, McDonald–Grebogi–Ott–Yorke, Ashwin et al., Berglund–Gentz, Bode, Freudenberg–Looze, Hoppensteadt, Khalil.

## Lesson

- **A control must be able to see the event it rules out.** Here the warm window ended one period before draw 7's alarm, and the count read "not there" as "passed".
- **A control that changes two things must be split.** The warm start moved both the start and the coefficient, and with the coefficient held fixed it reproduced the alarm.
- **A monotonicity count is vacuous on a constant sequence.** Require both outcomes, and the base reproduced at offset 0, before calling something a threshold.
- **When a referee's proposed text carves out an exception, carry the exception.** Generalising "11 at k_e = 20" to "every window" turned a correct reading into a false one.
