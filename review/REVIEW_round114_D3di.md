# Referee report D3di, round 114: Proposition D‴.60 as revised after round 111

## Scope

- **Commit.** I reviewed 5be79c1 on claude/vigilant-clarke-hf939b. I exported it with `git -C … archive 5be79c1 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref114.
  - HEAD moved to d5c3c29 while I worked. For the reviewed files, `git diff 5be79c1 HEAD` touches only the script and the tests. Commit 0e4f9e6 ("D‴.60 round-111 leftovers") takes `kk` from `c_run` and plants rows for R06 and R07. I did not review it. It bears on the deferred items in the status table.
- **Read in full:**
  - review/drafts/D60.md (386 lines), especially (iii) at lines 136–155 and the Numerics at 370–376;
  - in experiments/exp_noise_resonance.py, line by line: `c8a_job` (the starts 'warm', 'warm_p', 'warm_pk') and `c8a_lines` (`coord`, `near`/`neark`, the cold/warm lists, the labels of the windows that are not one coordinate's, the verb table, the noise clause). Also what they call: `offsets_start`, `state_start`, `win_run`, `edge_run`, `c_run`, `left_periods`, `c8a_cases`;
  - results/chaos/noise_resonance.md (C8A paragraph at line 681);
  - in tests/test_noise_resonance.py: the round-89 md test, the round-98 and round-102 tests, and `test_round105_…` with the round-111 additions;
  - tools/mutants/round111_D60_own.py and tools/mutants.py;
  - review/REVIEW_round111_D3df.md and the Round 111 entry of REVIEW_LOG.md;
  - `git diff 232ca83 5be79c1` of these files (5 files, +134 −68), and the drafts README line.
  - I also read the round-111 referee's mutant spec (ref111/mut/spec111.py) to reuse the R06, R07 and R11 anchors.
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv (symlinked as W/repo/venv);
  - NR_OUT=W/out holding a copy of the raw, TMPDIR=W/tmp;
  - pytest `--basetemp` inside W (through `PYTEST_ADDOPTS` for tools/mutants.py);
  - at most two heavy processes of mine at a time, each under `timeout`.
- **Nothing in the repository was modified.** I ran only read-only git commands and killed no process.
- **Tags:** [md] the committed md or draft; [raw] the committed raw; [own] my runs or code; [snippet] a web-search result only; [memory] not checked this round.

## Summary

**Design.** The revision answers D3df-1 … D3df-6 in four places.
- **`c8a_job`:**
  - 'warm' (δ(p − 1) from p − 1) is now run NP_C7 + 1 periods;
  - a twin 'warm_p' starts from p − 1 with δ(p);
  - 'warm_pk' starts two and three periods earlier with δ(p).
- **`c8a_lines`:**
  - A window counts as one coordinate's only if its scan is monotone, contains both outcomes, and offset 0 reproduces the reduced window. The other windows are named with a reason.
  - An earlier start "reproduces" the alarm if it leaves within ±1 period of the alarm's period.
  - "The cold restart's" and "one coordinate" are printed only for the windows the counts support.
  - Verb table, "0", the noise clause and the round citations are fixed.
- **Tests:** a mock with distinct rows for the rerun path, and planted rows for a flat scan, an all-removed scan, offset 0 differing from the reduced window, a non-monotone scan, a warm start one period late, and a k = 2 start.
- **D60 (iii) and the Numerics:** rewritten. The 11 alarms at k_e = 20 are now separated from draw 7's.

**What holds**
- **Reproduction is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 5f8eea4a…);
  - 80 tests pass (112.7 s);
  - all 7 C8A jobs, rerun with the script's own `c8a_cases` and `c8a_job`, equal raw['C8A'] in every field except `secs`. This includes draw 7 at k_e = 10, whose 'warm_pk' rows are [2, 12, 0.31014] and [3, 13, 0.30998].
- **Every new number in the md and D60 agrees with the raw verbatim** [raw, own recount]:
  - 7, 10, 11; 11, 9, 8, 5, 0, 0;
  - "11 of 12", and draw 7 "no offset removes it, its largest distance 0.3104";
  - "1 of 12 (draw 7)", "the 11 others";
  - draw 2 from 58: 'r+c1m' −1 against 'w+c1' 5.
- **Both readings at the 11 alarms at k_e = 20 hold on finer evidence** [own]:
  - On a grid of step 0.0025 inside each scan bracket (61 points, plus the 66 raw points reproduced bit-identically by my runner), each window has one removal threshold, 0.005–0.045 outward of the feature.
  - The 33 samples (period mean, pT, (p − 1)T) all lie on the side of that threshold their outcome implies.
  - Starting the reduced state 1–5 periods earlier, with δ(p) or with its own δ (110 runs), never leaves within one period of the alarm; the largest distance at the alarm's period is 0.045.
  - The draft's "by kernel 1's centre at 10 of the 11" is confirmed:
    - the warm loops' kernel 1 centre at pT is −0.0395 … −0.0981, outward of the threshold, at 10 windows;
    - at s = 0.06, draw 2 from 58 it is −0.0257 … −0.0269 over all 10 earlier starts, inward of (−0.045, −0.0425].
- **Draw 7's exclusion is right, and "the coefficient's" is confirmed and sharpened** [own]:
  - lowering δ(107) = 0.07667 by 0.0004 removes the alarm; the threshold lies in (0.07627, 0.07637];
  - the run's own coefficient at 107T, 0.07613, lies below it.
- **No major issue.**

**Findings**
1. **D3di-1 (minor).** The headline "At 11 of the 12 that split is one coordinate's" holds on the reduced base, not on the split's own base. At draw 2 from 58 the split ('wc' delays, 'i' removes) is decided by kernels 2 and 3's instantaneous state, not by kernel 1's centre. The draft says so itself one line earlier ("kernels 2 and 3's state 1").
2. **D3di-2 (minor).** Draw 7's clause "started earlier with δ(p) leaves within one period" omits that from one period earlier with δ(p) the loop passes (raw 'warm_p' [−1, 0.0304]). What decides is kernel 1's coefficient at 107T against a threshold 0.0004 below δ(107): a second margin, besides 0.3104 against 0.3.
3. **D3di-3 (minor).** The δ(p) control also changes the coefficient at the alarm. An earlier start with δ(p) arrives at pT with δ(p) + 0.0008 … +0.0093 (k_e = 20). At k_e = 20, s = 0.08, draw 1 from 50 these starts leave at period 54, two periods after the alarm's 52 (0.30001, 0.30036; 0.30091 and 0.30126 from 4 and 5 periods back), while the own-δ starts pass. This is not reported.
4. **D3di-4 (minor).** The new conditions still cannot tell some readings apart:
   - the "cold restart's" condition counts an earlier start that leaves more than a period *before* the alarm's period as having passed it (planted P1, P5);
   - "one coordinate" does not check the samples against the scan (planted P7);
   - the warm loops' centres and the instantaneous samples that the draft cites are not in the raw.
5. **D3di-5 (minor).** Tests. The authors' 10 mutants are all killed. Of my 21:
   - 8 are killed;
   - 13 survive the full suite: L01–L06, L08 and L09 (each killed by one planted row of mine), L07 and J04 (equivalent), and R06, R07 and R11 (round 111's deferred survivors).
6. **D3di-6 (cosmetic).** Wording and robustness:
   - stale docstrings;
   - "with δ(p − 1) or δ(p), run to the alarm's period";
   - D60:124's "the coming noise";
   - the `p = 1` wrap in 'warm', and `offs.index(0.0)`;
   - a redundant `not all(seq)`;
   - the delayed class at the threshold (draw 2 from 58: −0.0425 → period 4).

## What I checked and how

All paths are under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref114.

- **Reproduction:**
  - W/logs/fromraw.log, then `cmp` against the committed md: identical.
  - W/logs/pytest_full.log: 80 passed, 1 old warning.
  - W/scripts/repro_c8a.py → W/logs/repro_0246.log and W/logs/repro_135.log: 7 of 7 jobs identical (9.7–35.4 s each).
- **Recount (numpy only, `python -I`):**
  - W/scripts/dump_c8a.py → W/logs/dump_c8a.log;
  - W/scripts/recount.py → W/logs/recount.log. It covers every count of the new sentence, the earlier starts' absolute leave periods, draw 2 from 58's base dependence, and C8D's dense windows around draw 7 from 107 and draw 1 from 50 (s = 0.08).
- **Own probes.** My runner uses the repository's `simulate`, `raised` and `c_noise`; the start states, period means (64 outputs a period) and first-period-beyond-0.3 are my own code.
  - W/scripts/probe.py, `grid` mode → W/logs/probe_grid.jsonl: the threshold grid, the run's instantaneous centre at pT and (p − 1)T from `c_run`, and the 66 raw scan points rechecked.
  - W/scripts/probe.py, `starts` mode → W/logs/probe_starts.jsonl: per window, k = 1…5 periods earlier with δ(p) and with δ(p − k), NP = 10 + k. Each run records its leave, its distance at the alarm's period, and kernel 1's centre and |w| offset at pT. The k = 1 and draw-7 rows reproduce the raw's 'warm', 'warm_p' and 'warm_pk' to all digits.
  - W/scripts/analyse_starts.py analyses both.
  - W/scripts/draw7_coef.py → W/logs/draw7_coef.jsonl: draw 7's cold restart at 107 with δ(107) − x, x = −0.001 … 0.001, and the run's coefficient at 106T and 107T.
- **Planted inputs:**
  - W/scripts/planted.py → W/logs/planted.log: P1–P8.
  - W/scripts/killplants.py → W/logs/killplants.log: one planted row per surviving `c8a_lines` mutant, each mutant applied to a private copy of the module.
- **Mutants:**
  - The authors' 11 (the control included) plus 21 of my own (W/mut/spec114.py, combined in W/mut/allspec.py) on `-k "round89_md or round98 or round102 or round105"` → W/mut/w1/res_a.jsonl and W/mut/w1/res_b.jsonl.
  - The 13 survivors on the full suite → W/mut/wfull/res_a.jsonl and W/mut/wfull/res_b.jsonl (80 passed each).
- **Literature:** three web searches (initialisation shock, warm and cold starts, finite-time basin stability).

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii); (iii) outside C8A; Proof; Status | unchanged since 232ca83 (the diff touches only 143–155 and 372–376); not reverified | [md] |
| md C8A paragraph (md:681), round-111 sentence | verified; regenerated byte for byte | [raw], [own] |
| "At 11 of the 12 that split is one coordinate's" (143–144) | true on the reduced base only. On 'wc''s base, at draw 2 from 58: 'wc' 5, 'wc+i1c' 2 (kernel 1's instantaneous centre restores the alarm at its own period), 'wc+i23' −1. 'r+c1i' is 2 on the reduced base too. D3di-1 | [raw] |
| 7, 10, 11; 5, 8, 9, 11; none at 0 and +0.04 (144–147) | verified | [raw] |
| "at the 11 at k_e = 20 … monotone … both outcomes, offset 0 the reduced window itself" (147–148) | verified. One threshold per window on the fine grid: (−0.03, −0.0275], (−0.045, −0.0425], (−0.0175, −0.015], (−0.01, −0.0075], (−0.0075, −0.005], (−0.02, −0.0175], (−0.045, −0.0425], (−0.01, −0.0075], (−0.0325, −0.03], (−0.03, −0.0275], (−0.0275, −0.025] | [raw], [own] |
| "the period mean and the instantaneous value two samples of it" (148) | true, 33 of 33 samples consistent. Not computed by the code; the instantaneous values are not in the raw. D3di-4 | [own] |
| "started one, two or three periods earlier … never leaves within one period of the alarm's: those alarms are the cold restart's" (148–150) | verified (raw: 44 of 44 starts; own: 110 of 110 with k = 1…5). At draw 1 from 50 (s = 0.08) the δ(p) starts leave two periods later, unreported. D3di-3 | [raw], [own] |
| "(by kernel 1's centre at 10 of the 11 … draw 2 from 58 … inward …)" (150–151) | verified, but the evidence is the referee's, not the raw's | [own] |
| "draw 2 from 58 is removed by kernel 1's mean centre on the reduced base and only delayed on the 'w' base" (152) | verified ('r+c1m' −1; 'w+c1' 5) | [raw] |
| Draw 7: "no offset of the centre removes it (0.3104 against 0.3: a marginal crossing)" (153–154) | verified (10 at all six offsets, 0.31035–0.31038) | [raw] |
| "the reduced state started earlier with δ(p) leaves within one period … the coefficient's, δ(107), not the restart's" (154–155) | the reading holds and is sharpened (threshold in (0.07627, 0.07637]). From one period earlier with δ(p) it passes (0.0304), which is omitted. D3di-2 | [raw], [own] |
| Numerics C8A (372–376) | verified; the same caveats apply | [raw] |
| drafts README round-111 line | consistent; "10 own mutants killed" matches the log's 9 plus T05 after its plant | [md] |
| Sources | unchanged; Schultz et al. consistent with the snippet | [snippet] |

The new conditions against the readings, with planted inputs [own, W/logs/planted.log]:

| plant | condition | the code's output | correct? |
|---|---|---|---|
| flat scan, all-removed scan, offset 0 ≠ reduced, non-monotone | coordinate | excluded and named | correct |
| P2: only 'warm_p' at the alarm's period | ±1 | counted | correct |
| P3: only k = 3 exactly at the alarm's period | ±1 | counted | correct |
| P4: warm two periods late | ±1 | not counted | correct |
| P1: every earlier start leaves in its own period 1, the alarm in period 5 | cold | prints "the alarms are the cold restart's" | **wrong** |
| P5: the k = 3 start leaves two periods before the alarm | cold | prints "the alarms are the cold restart's" | **wrong** |
| P7: the mean centre at −0.08 (outward of every removing offset) yet 'r+c1m' keeps the alarm | coordinate | still prints "one coordinate against a per-window threshold" | **wrong** |
| P6: the kept offsets delayed | coordinate | prints "one coordinate" (removal only) | acceptable |
| P8: no 0 in the scan | label | raises ValueError | crash |

Mutants:

| set | killed | survive the full suite |
|---|---|---|
| authors' 10 (T01–T08, R01, R02) | all; T04 and T07 also change the md | none |
| J-mutants of `c8a_job` | J01, J02, J05, J07, J08 (by test_round102's draw-7 rerun against the raw); J03, J06 (by the wiring test) | J04 (`p − k − 1 > 0`): equivalent on the data, smallest p = 6 |
| L-mutants of `c8a_lines` | L10 | L01–L09; L07 equivalent (the offset-0 check already excludes an all-removed scan) |
| round 111's deferred | — | R06, R07, R11 |

## Status of the round-111 issues

| D3df id | status at 5be79c1 | what remains |
|---|---|---|
| D3df-1 (major) | **fixed in substance.** Draw 7 excluded and named; the coordinate needs both outcomes and offset 0 the reduced window; the earlier starts now cover the alarm's period; D60 rewritten | the headline's base (D3di-1); draw 7's one-period pass and its coefficient margin (D3di-2); the δ(p) control's drift (D3di-3); early leaves counted as passes (D3di-4) |
| D3df-2 (minor) | fixed: "C8D's window at p − 1, lengthened"; ±1; the δ(p) twins; the draw 2 from 58 exception stated and confirmed here [own] | it rests on referee evidence not in the raw (D3di-4) |
| D3df-3 (minor) | fixed: R01 and R02 killed by `test_round105_…` [own] | none |
| D3df-4 (minor) | partly. R03 and R13 recorded as equivalent; the 'needs' set planted | R06, R07, R11 still survive the full suite at 5be79c1 [own]. 0e4f9e6 (after 5be79c1, not reviewed) plants R06 and R07; R11 has no plant at HEAD |
| D3df-5 (cosmetic) | fixed in the md: "delays", "0", the noise clause, "(referee rounds 98, 102, 105, 111)", the base-dependence clause | `c8a_lines`' docstring still says "(referee rounds 98, 102)" (line 2720); D60:124 still says "needing the coming noise" |
| D3df-6 (cosmetic) | draw 7's margin stated (0.3104 against 0.3) | `kk = 8` remains at 5be79c1 (line 1591); 0e4f9e6 takes it from `c_run` (not reviewed) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3di-1 | minor | **"At 11 of the 12 that split is one coordinate's" names the wrong base** (D60:143–144).<br>• The 11 are monotone on the *reduced* base.<br>• "That split" is 'wc' against 'i', on 'wc''s base. At draw 2 from 58 [raw]: 'wc' delays (5); 'wc+i1c' gives the alarm back at its own period (2); 'wc+i23' removes it. On the reduced base 'r+c1i' keeps it too (2).<br>• So there the instantaneous removal is kernels 2 and 3's, as D60:143 itself counts ("kernels 2 and 3's state 1").<br>• On its own base the split is kernel 1's centre's at 10 of 12: 6 via 'w+c1' and 4 via 'wc+i1c'. Draw 7's is the coefficient's; draw 2 from 58's is kernels 2 and 3's. | **D60 (iii):** "On the reduced base kernel 1's centre is one coordinate at the 11 alarms at k_e = 20 … On 'wc''s base the split is kernel 1's centre's at 10 of the 12: draw 7's is its coefficient's, draw 2 from 58's kernels 2 and 3's state's (there kernel 1's instantaneous centre keeps the alarm on either base)." Optionally add the fine thresholds, 0.005–0.045 outward (referee round 114). |
| D3di-2 | minor | **Draw 7: the δ(p) start one period earlier passes, and the coefficient margin is unstated** (D60:153–155).<br>• [raw] 'warm_p' (from 106 with δ(107)) is [−1, 0.0304]; only the starts two and three periods earlier leave (117).<br>• [own] The cold restart at 107 leaves with δ(107) − x for x ≤ 0.0003 and passes for x ≥ 0.0004, so the threshold is in (0.07627, 0.07637].<br>• [own] The run's own coefficient at 107T, 0.07613, is below it, which is why 'wc+i1w' and 'i' remove the alarm.<br>• [own] Earlier starts arrive at 107T with 0.07544–0.07627 (they pass: δ(p) from 106, and own δ from 102–106) or with 0.07724–0.07893 (they leave: δ(107) from 102–105).<br>• So the reading "the coefficient's" holds; the earlier δ(p) starts reproduce the alarm through their drift above δ(107). | **D60 (iii):** "… it is the coefficient's: the period mean δ(107) = 0.0767 lies 0.0005 above the run's coefficient at 107T (0.0761), across a threshold in (0.0763, 0.0764]. The reduced state started two and three periods earlier with δ(107) arrives above it and leaves at 117; from one period earlier it arrives below it and passes (0.0304)." Same in the Numerics. |
| D3di-3 | minor | **The δ(p) control moves the coefficient at the alarm; one consequence is unreported** (`c8a_job`, md:681, D60:148–150).<br>• [own] δ(p) starts arrive at pT with δ(p) + 0.0008 … +0.0093 at k_e = 20 (k = 1…5); own-δ starts with −0.0014 … +0.0017.<br>• This is conservative for "the cold restart's" (they pass anyway) and liberal for "the coefficient's" (D3di-2).<br>• [raw] At k_e = 20, s = 0.08, draw 1 from 50 'warm_pk' leaves at 54, two periods after the alarm's 52 (0.30001, 0.30036). [own] From 4 and 5 periods back: 54 (0.30091, 0.30126). The own-δ starts pass (0.2855–0.2864), as do C8D's windows 44–52 except 50.<br>• "Never leaves within one period" is literally true but hides a flagged control window. | **`c8a_job`:** record kernel 1's |w| and centre offsets at pT for each earlier start (a short `simulate` to pT, or a state return from `edge_run`).<br>**`c8a_lines`:** print the earlier starts that leave outside ±1 by name.<br>**D60:** "(at s = 0.08, draw 1 from 50 the δ(p) starts from 48 and 47 leave at 54, at 0.3000 and 0.3004, carried above δ(50) by their own drift)". |
| D3di-4 | minor | **The conditions cannot tell "passed" from "left before" or a contradicting sample** (`c8a_lines` `near`/`cold`, `coord`).<br>• [own, P1, P5] An earlier start that leaves more than one period before the alarm's period counts as having passed, and "the alarms are the cold restart's" prints. Round 111's "whose warm run reaches the alarm's period" was implemented as a length, not as a condition.<br>• [own, P7] A mean-centre sample outward of every removing offset that keeps the alarm still prints "one coordinate against a per-window threshold".<br>• D60's "the period mean and the instantaneous value two samples of it" and "by kernel 1's centre at 10 of the 11" rest on referee runs. The raw stores neither the instantaneous centres nor the warm loops' centres at pT.<br>• All of this is equivalent on the data: no earlier start leaves early, and 33 of 33 samples are consistent [own]. | **`c8a_lines`:** a window counts as cold only if every earlier start either stays or leaves after alarm + 1, i.e. `w[0] < 0 or w[0] > x['reduced'] + k + 1`. Otherwise name it "not tested (an earlier start left at …)".<br>**`c8a_job`:** store `c1i`, `c1p` (the offsets at pT and (p − 1)T) and the earlier starts' centre at pT.<br>**`c8a_lines`:** check each of 'r+c1m', 'r+c1i', 'r+c1p' against the scan's bracket before printing "one coordinate", and print the warm centres against the bracket. |
| D3di-5 | minor | **Survivors of the full suite** [own, W/mut/wfull]:<br>• L01: 'warm_p' ignored;<br>• L02: 'warm_pk' judged at reduced + 1;<br>• L03, L04: `near` one-sided;<br>• L05: `neark` with no lower bound;<br>• L06: `near` at reduced + 2;<br>• L08: the label's distance replaced by the scan's largest;<br>• L09: `neark` at reduced + k − 1;<br>• L07: `not all(seq)` dropped — equivalent;<br>• J04: `p − k − 1 > 0` — equivalent on the data;<br>• R06, R07, R11 (round 111's).<br>The planted 'near' (warm at R + 2) and 'deep' (k = 2 at R + 2) rows each pin one side of one condition only. | Add the rows of W/scripts/killplants.py, each verified to kill its mutant:<br>• 'warm_p' only at R + 1 (L01);<br>• k = 3 at R + 3 (L02);<br>• warm at R + 3 (L03), R − 1 (L04), R (L06);<br>• k = 3 at R + 1 (L05), R + 4 (L09);<br>• a flat scan whose offset-0 distance differs from its largest (L08).<br>Run the wiring test with p = 4 (J04). Plant R11 (R06 and R07 are planted at 0e4f9e6). Record L07 as equivalent, or drop `not all(seq)`. |
| D3di-6 | cosmetic | **Wording and robustness:**<br>• `c8a_lines`' docstring "(referee rounds 98, 102)" (line 2720);<br>• `c8a_job`'s docstring describes only 'warm' (not 'warm_p', 'warm_pk', the lengthening);<br>• `test_round105_…`'s docstring ("the cold restart's only if no warm start keeps the alarm; one coordinate only if every scan is monotone") is stale;<br>• D60:149 "(with δ(p − 1) or δ(p), run to the alarm's period)": δ(p − 1) is used one period earlier only, and the runs go to p + 10. The md's "(C8D's window at p − 1, lengthened)" fits the δ(p − 1) start only;<br>• D60:124 "its needing the coming noise" (the md now says "the noise (other noise not recorded)");<br>• 'warm' and 'warm_p' are unguarded at p = 1 (`W[p − 2]` wraps to the last period), while 'warm_pk''s guard excludes the valid p = k;<br>• `offs.index(0.0)` raises without 0 in the scan (P8);<br>• the reading is about removal only: at draw 2 from 58 the fine grid gives −0.0425 → period 4, as 'wc' delays it (2 → 5) [own];<br>• `kk = 8` (fixed after 5be79c1). | Update the docstrings; "one period earlier with δ(p − 1) or δ(p), two or three earlier with δ(p), each run to the window's end"; "the noise" in D60:124; guard p ≥ 2 and p − k ≥ 0; use `next(...)` for offset 0; one clause on the delayed class at the threshold. |

## Sources

- D. Mulholland, P. Laloyaux, K. Haines, M. Balmaseda, "Origin and impact of initialisation shocks in coupled atmosphere-ocean forecasts", Mon. Wea. Rev. 143(11), 4631–4644 (2015) [snippet].
  - A related search snippet separates drift, which "occurs regardless of the initialisation method", from shock, which "comes from an imperfect initialization method".
  - This is the analogy behind D60's "the cold restart's" (a shock of starting from the reduced, inconsistent state) against "the coefficient's" (the drift). An earlier start is the usual way to separate the two. D60 could cite it as an analogy.
  - Links: [CentAUR record](https://centaur.reading.ac.uk/40638), [ECMWF slides (Balmaseda)](https://events.ecmwf.int/event/335/contributions/3919/attachments/2414/4177/TCPR_Balmaseda_2023_Initialization.pdf).
- Cold against warm starts and spin-up in NWP [snippet]. A cold start from an analysis without a model-consistent fine-scale state spins up with a transient; a warm start reuses earlier model fields and shortens it.
  - Link: [UCAR record, reducing the spin-up of a regional NWP system](https://impacts.ucar.edu/en/publications/reducing-the-spin-up-of-a-regional-nwp-system-withoutdata-assimil/).
- P. Schultz, F. Hellmann, K. N. Webster, J. Kurths, "Bounding the first exit from the basin: independence times and finite-time basin stability", Chaos 28, 043102 (2018) [snippet].
  - The independence time is the time at which finite-time basin stability reaches its asymptotic value. The earlier starts here are a crude version of this: at the 11 alarms one period suffices for the restart's transient to stop deciding.
  - It is consistent with D60's existing citation.
  - Links: [PIK record](https://publications.pik-potsdam.de/pubman/item/item_22275), [arXiv 1711.03857](https://ar5iv.labs.arxiv.org/html/1711.03857).
- Not rechecked this round [memory or earlier rounds' snippets]: Feudel et al. 1998 (round 111's snippet), and the remaining citations of D60:230–245.

## Lesson

- **A control that moves the start moves what it carries to the event.** The earlier start with the window's δ(p) reached the alarm with a different coefficient. Draw 7's "reproduction" and draw 1's late leave come from that drift. Record the control's state at the event, not only its outcome.
- **A pass must be a pass through the event.** Lengthening a control lets it see the event's period, but a control that left earlier still never saw it. Put "reached the event" into the condition, not only into the run length.
- **A headline count must be on the base of the split it names.** "11 of 12" is true on the reduced base and 10 of 12 on the split's own base.
- **Plant both sides of every boundary.** One planted row per condition pinned ±1 on one side and k = 2 only, and eight single-point changes of the new comparisons survived the full suite.
