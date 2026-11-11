# Referee report D3cs, round 98: Proposition D‴.60 as revised after round 94

## Scope

- **Commit.** 99edc3d on claude/vigilant-clarke-hf939b, exported with `git -C … archive 99edc3d | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref98.
  - While I worked, other sessions moved HEAD to c72a571 and then fd9efda.
  - `git diff --stat 99edc3d HEAD` is empty for the draft, the script, the md, the raw, the tests and the round-94 mutant spec, so everything below holds at HEAD too.
- **Read in full:**
  - review/drafts/D60.md (351 lines);
  - results/chaos/noise_resonance.md (710 lines);
  - review/REVIEW_round94_D3co.md and the Round 94 entry of review/REVIEW_LOG.md (lines 3706–3739);
  - tools/mutants/round94_D60_own.py and tools/mutants.py;
  - `git diff 89e9e28 99edc3d` for these files.
- **tests/test_noise_resonance.py (75 tests):** the round-85, round-89 and round-94 tests line by line, the rest by name and docstring.
- **experiments/exp_noise_resonance.py:** all code changed since 89e9e28, line by line:
  - `c7l_lines`' header and superset sentence;
  - `c7s_lines`' least raise, overlap and far pairs;
  - the C8R miss list in `c8_more_lines`;
  - the C8D sentence;
  - `c8d_alarms`' repeat at a ± 1.

  I also read what these call: `gap_devs`, `c7_job`, `c7s_job`, `edge_run`, `win_run`, `c8_delta`, `c8r_job`, `c8d_job`, `offsets_start`, `c8k_job`, `c8k_cases`, `c_run`, `c_noise`, `left_periods`, `main`, and the stepping of the C integrator (`closed_loop.c`).
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with the repository's venv;
  - NR_OUT=W/out holding a copy of the raw, LCK_BUILD=W/build, and pytest `--basetemp` inside W;
  - at most two heavy processes of mine at a time.
- **Nothing in the repository was modified.** I ran only read-only git commands. The only things I added are under W: a symlink W/repo/venv to the venv, so that tools/mutants.py can run from my export. I killed nothing. The `tools/mutants.py --exp exp_moving_width` processes and the modified working-tree files belong to another session.
- **Tags:** [md] the committed md; [raw] the committed raw; [own] my runs or my own code; [snippet] a web-search result only; [memory] not checked this round.

## Summary

**Design.** This round reviews the revision after round 94, which made these changes:
- the C7L superset restated as true by construction;
- the start's effect on the set's low end: least-raise change ×4.08, median ×1.019; Jaccard overlap 0.975 against 0.831;
- "isolated" defined as a lead other than one period, with repeats tested at a ± 1 (0 of 106);
- the run's own distance cited as the evidence for "not the state's";
- draw 5's kernel 1 described as displaced;
- C8R's never-reproduced runs named;
- new sources.

My checks:
- **Reproduction:** regenerated the md from the raw, ran the tests, and reran 17 jobs with the script's own functions, comparing them to the raw.
- **Numbers:** checked every number new to the draft against the md and the raw, and recounted every new aggregate with my own code.
- **Logic:** three probes.
  1. Alternative pairings for the far comparison.
  2. The run's distance at each isolated alarm, and the alarm windows restarted from fuller states.
  3. Draw 5's window restarted with each formed kernel's offset separately.
- **Mutants:** 19 of my own, plus the authors' 7, against the full suite.
- **Literature:** searched the new sources.

**What holds**

- **Reproducibility is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 7b1bb64e5b717f164db1fa7f6a98705e).
  - The 75 tests pass in 105 s.
  - Every rerun equals the raw in every field except `secs`: two C7S scans, one C7L scan, six C8D windows (three of them alarms), two C8R jobs and one C8K job. The five C runs at k_e = 20 have bit-identical kernel-1 coefficients, leave periods and largest distances.
- **Every number new to the draft is printed by the md** [md], except "0.038–0.111" (D3cs-2).
- **My recounts agree** [own]:
  - ×4.0816 (s = 0.03, draw 7, 50 against 49); the median ×1.0193 over 26 neighbour pairs.
  - Jaccard median 0.9748 (least 0.333); the far median 0.8313 over 420 pairs.
  - 21 alarms = 9 with lead 1 + 12 with leads 11–86; 106 neighbours, 0 repeats.
  - C8R covers 9 alarms and passes all 9; C8R's never-reproduced runs are exactly the ones the md names.
- **"By construction" is right.** The C integrator clips its adaptive step to the next output time, and the final time enters only a diagnostic. The 10-period output grid is the 20-period grid's prefix, so the prefix runs are identical. My C7L rerun keeps every 10-period leave period [own].
- **"Displaced" is right, and more is true** [own]. Restarting draw 5's window from 150:
  - with kernel 3's offset alone, kernel 1 left at its fitted amplitude: kernel 3 leaves at 158 and kernel 1 at 160, as in the run;
  - with kernel 1's offset, or kernels 1 and 2's: nothing leaves.

**What does not hold (no major issue)**

1. **D3cs-1 (minor).** The comparison 0.975 against 0.831 is not like for like.
   - The neighbours are 26 pairs at the 13 gap windows, in 6 configuration-levels.
   - The 0.831 pools 420 pairs over all 15 configuration-levels. Nine of those hold no C7S window, among them near-noiseless levels with medians 0.984.
   - Anchored on the same base windows the far median is 0.706 (round 94's figure); over the same configuration-levels it is 0.756 [own]. The conclusion stands and is stronger like for like.
2. **D3cs-2 (minor).** The evidence for "not the state's":
   - "0.038–0.111" is produced neither by the md nor by the raw, and it is not the run's distance at the alarm's period: that is 0.0001–0.081 at k_e = 20 [own rerun] and 0.032 for draw 7 [raw C3].
   - Any value below the cut is implied by what an alarm is, so it says nothing about how close the run is to ejecting. Leaves are jumps: two periods before each of the 22 C3 band leaves, kernel 1 is within 0.000–0.048 of its feature [raw].
   - The discriminating restarts are not used: 'w' keeps all 11 alarms at k_e = 20; 'wc' removes 6, delays 1, and 4 remain; draw 7 survives both [own].
3. **D3cs-3 (minor).** Mutants: 9 of my 19 survive the full suite, all equivalent on the data or by construction, none changing the md or a rerun. Four are worth a planted edge; one (R04) points to a latent bug, which drops a least raise of 0 silently. Of the authors' 7, 5 are killed, and B02 and B05 survive as recorded [own].
4. **D3cs-4 to D3cs-7 (cosmetic):**
   - the md's C8D sentence attaches the early replay's explanation to "isolated"; `byc` still drops the metric;
   - C8R's 37 mixes in-sample runs into an out-of-sample sentence;
   - draw 5's exception to "a leave needs both" is not stated;
   - source details (year, the correction, forward against finite-time basins).

## What I checked and how

All paths are under W.

- **Reproduction:**
  - logs/fromraw.log, then `cmp` against the committed md: identical.
  - logs/pytest_full.log: 75 passed.
  - scripts/repro_jobs.py → logs/repro_jobs.log, comparing whole records with `secs` excluded:
    - `c7s_job` for s = 0.03 from 49 and for k_e = 20, s = 0.04 from 26;
    - `c7_job` over 20 periods for s = 0.03 from 175, plus both prefix checks;
    - `win_run` for k_e = 20, s = 0.08, draw 1 from 50, 58, 59; k_e = 20, s = 0.06, draw 3 from 166, 167; k_e = 10, s = 0.02, draw 7 from 107;
    - `c8r_job` for k_e = 20, s = 0.06, draw 7 and for k_e = 10, s = 0.02, draw 10;
    - `c8k_job` for s = 0.03, draw 5.
  - scripts/repro_runs.py → logs/repro_runs.log: `c_run` for k_e = 20 at s = 0.06 (draws 1, 2, 3) and s = 0.08 (draws 1, 2), against C8's `w`, `left` and `maxdev`.
- **Own recounts from the raw (numpy only, `python -I`):**
  - scripts/own_c7s.py → logs/own_c7s.log: neighbour Jaccard and least ratios; far pairs pooled, by configuration-level and anchored; Jaccard by raise length.
  - scripts/own_c8d.py → logs/own_c8d.log: alarms, leads, neighbours, repeats at ± 1 … ± 10, neighbours gone before a − 1, the windows at p ± 1, C8R coverage, C8R lags and the runs it never reproduces.
  - scripts/own_c8r_oos.py: C8R in and out of sample.
  - scripts/own_c3.py → logs/own_c3.log: period means around draw 5's and draw 11's leaves, and draw 7 at 116–118.
  - scripts/own_preleave.py → logs/own_preleave.log: kernel 1's distance before each C3 band leave.
- **Probes:**
  - scripts/repro_jobs.py `probe5` → logs/probe5.json: draw 5's window from 150 started from kernel 1 only, all three, 1 + 3, 1 + 2 and 3 only.
  - scripts/probe_wc.py → logs/probe_wc.log: the 12 isolated-alarm windows restarted with `offsets_start` as 'w' and 'wc'.
- **Mutants:**
  - mut/spec98.py: 19 and a control, run with tools/mutants.py from W/repo.
  - Stage 1 (`-k "round89 or round94 or c8_more_lines or c7l or c7_gap"`): mut/w_s1/res_s1a.jsonl.
  - Stage 2, the survivors on the full suite with `--from-raw` and `cmp`: mut/w_s2/res_s2a.jsonl.
  - The authors' tools/mutants/round94_D60_own.py on the full suite: mut/w_own/res_s1b.jsonl.
- **Literature:** web searches on Oljača–Ashwin–Rasmussen 2024 and the 2022 slides, and on Giesl–Rasmussen. Fetching ore.exeter.ac.uk and www.mpecdt.ac.uk failed (host outside this session's allowed network), so the full texts were not read.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii) and their proofs (9–49, 149–172) | unchanged since round 94, not reverified this round | [md] |
| The noise-free edge, power law, outside the window, C4, s* (51–83) | unchanged; numbers printed | [md] |
| Bands, cuts, 12 basin and 6 cut gaps (83–90) | unchanged | [md] |
| C7L: 14 of 23 change, 13 multi-band, 10 basin; the selection 18 + 4 + 1 (90–92) | verified | [md], [own] |
| "By construction … can only grow (23 of 23 a check of determinism)" (92–93) | correct; the integrator's steps do not depend on the horizon | [own] |
| 10 gaps survive, 3 close (93–94) | verified | [md] |
| C7S 6 / 7 / none; 99 → 0.56–1; 74 → 0.085–1 (94–96) | verified | [md] |
| "The start moving the set's low end", ×4.08, median ×1.019 (97–98) | verified: 10 of 26 pairs move more than ×1.1, 4 more than ×1.5, 6 have identical sets | [own] |
| "Jaccard median 0.975, against 0.831 for windows … 25 or more periods apart" (98–100) | numbers right; comparison not like for like | D3cs-1 |
| C8, C8N, the δ sweep, the noise-free edge, the six misses (101–116) | unchanged | [md] |
| C8R 33/37, 33/37, 29/33; "the other 4 never …" (105–108) | verified; 4 of the 37 are in sample | [own], D3cs-5 |
| C8D phases; the three robust misses (116–119) | unchanged | [md] |
| Draw 5: "reaches the first feature while kernel 1 is still at it … displaced" (119–120) | verified [raw C3]: at 159 kernel 3 at −0.393, kernel 1 at −0.613; at 160 kernel 3 at −0.487, kernel 1 at −0.971. Confirmed causally | [own] probe |
| Draw 11's contrast; C8K's recovery (120–122) | verified | [raw], [own] |
| 9 early + 12 isolated, leads 11–86, period 2 or last; "isolated meaning a lead other than one period" (122–127) | verified; classes separate cleanly (leads 1 or ≥ 11) | [own] |
| "None of the 106 … repeats it" (127–128) | verified; 0 within ± 5 periods; 1 of the 106 had left 8 periods before | [own] |
| "C8R … passes all 9 — implied … weaker than the run's own kernel 1 there, 0.038–0.111 …: not the state's" (128–130) | number not reproducible and wrong as stated; the evidence is implied by definition | D3cs-2 |
| C9, C9A, single tones (130–158) | unchanged | [md] |
| Status (iii) (185–196) | accurate | — |
| Classical: non-autonomous basins, Giesl–Rasmussen, "monotone by definition" (216–220) | apt; details missing | D3cs-7 [snippet] |
| New: "33 of 37 … from 20 and 40 periods before, and 29 of 33 from 80" (226–228) | verified | [md] |
| Numerics C7L / C7S (325–329) | printed; the comparison as D3cs-1 | [md] |
| Numerics C8D (339–342) | verified | [own] |

## Status of the round-94 issues

| D3co id | status | what remains |
|---|---|---|
| D3co-1 | fixed | `c7l_lines`, (iii) and the Numerics say "by construction". `test_round94_main_c7s_superset_repeat_and_overlap` asserts the prefix leave periods on all 23 windows [own: passes] |
| D3co-2 | fixed in substance | The least raise and the overlap are printed and in (iii). Remaining: the far comparison uses a pooled pairing, not the like-for-like one (D3cs-1). The sub-period scan was not added (it was optional) |
| D3co-3 | partly fixed | "Isolated" is defined; repeats are tested at a ± 1 and printed; the Numerics are fixed. Remaining: the "direct evidence" cited is the round-94 referee's number, not reproducible and not the distance at a, and like C8R's passes it is implied by definition. The 'wc' restart was not added (D3cs-2) |
| D3co-4 | fixed | My probe confirms it (D3cs-6 adds a sentence) |
| D3co-5 | fixed | Q30 is killed by the new test. Q05–Q08, Q14, Q20, Q23, Q26 and Q29 were recorded as equivalent in the log instead of getting the proposed planted edges; Q05's `rep_alarm` is no longer printed in the md. The authors' B01, B03, B04, B06 and B07 are killed; B02 and B05 survive, as recorded [own] |
| D3co-6 | mostly fixed | C7L header, "a gap of the basin survives" and the alarms' sentence are done. Remaining: `byc` still drops the metric, and the new clause misplaces a parenthesis (D3cs-4) |
| D3co-7 | fixed | — (in-sample detail: D3cs-5) |
| D3co-8 | mostly fixed | Remaining: Giesl–Rasmussen's year and pages (D3cs-7) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cs-1 | minor | **The far overlap is not like for like** (D60:98–100, 328–329; md C7S; `c7s_lines`' `far`).<br>• **The two sides differ:**<br>&nbsp;&nbsp;– the neighbours are 26 pairs around the 13 windows selected for a gap of the basin, in 6 configuration-levels;<br>&nbsp;&nbsp;– the 0.831 is the median of 420 pairs over all 15 configuration-levels of C7.<br>• **Nine of the 15 hold no C7S window**, and several are nearly noise-free: k_e = 10 at s = 0.0025 and 0.005 (medians 0.984), k_e = 20 at s = 0.02 (0.925), WFR at s = 0.08 (0.906). Per configuration-level the medians run from 0.234 to 0.984 [own].<br>• **Like-for-like values:**<br>&nbsp;&nbsp;– each of the 12 C7 base windows against C7's windows of its own configuration, level and draw 25 or more periods away: median 0.706 (84 pairs; round 94's figure, recorded only parenthetically in the log);<br>&nbsp;&nbsp;– all pairs in the same 6 configuration-levels: 0.756 (168).<br>• **The least raise on the same anchors:** far median ×1.361 against ×1.019 for neighbours. Weighted by raise length rather than grid points, the Jaccard is 0.993 against 0.879 [own].<br>• **Effect:** the conclusion holds and is stronger, but the printed baseline understates the contrast, and the reader takes it for the matched comparison. | • In `c7s_lines`, anchor the far pairs on `win`'s base windows (each against C7's windows with the same `key[:5]` and \|Δshift\| ≥ 25), or restrict to the base windows' configuration-levels.<br>• Print 0.706 and, optionally, the far least-raise median ×1.36.<br>• Update (iii) (D60:98–100), the Numerics (328–329) and the md.<br>• Add a test on the raw that asserts the anchored median. |
| D3cs-2 | minor | **"Not the state's": the cited evidence** (D60:128–130).<br>• **(a) The number.** "The run's own kernel 1 there, 0.038–0.111" comes from the round-94 referee [own] and is printed neither by the md nor by the raw: C keeps no period means for the k_e = 20 runs. From my reruns (bit-identical to C8's `w`, `left`, `maxdev`), kernel 1's distance at each isolated alarm's period a is:<br>&nbsp;&nbsp;– at k_e = 20: 0.0071 (s = 0.06, draw 1, a = 205), 0.0454 (draw 2, 60), 0.0388, 0.0119, 0.0001 (draw 3, a = 168, 195, 212), 0.0814, 0.0080, 0.0278 (s = 0.08, draw 1, 43, 52, 60) and 0.0061, 0.0097, 0.0381 (draw 2, 8, 20, 34);<br>&nbsp;&nbsp;– draw 7 at k_e = 10 (a = 117): 0.032 [raw C3].<br>So the range is 0.0001–0.081. The 0.111 is the period after a = 52.<br>• **(b) The logic.** An alarm is a window that leaves at a ≤ P − 11, where P is the run's first leave, so the run is under the cut at a by definition. The draft rightly calls C8R's passes "implied by its aggregate", and the same holds here. The distance also says nothing about how close the run is to ejecting: leaves are jumps, and two periods before each of the 22 C3 band leaves kernel 1 is within 0.000–0.048 of its feature [raw]. The full state is the run itself, which by definition does not leave at a.<br>• **(c) The informative evidence is unused:**<br>&nbsp;&nbsp;– the repeat test: 0 of 106, or 0 of 105 at risk (window 50 had left at 52, 8 periods before the alarm at 60); none within ± 5; the windows at p ± 1 never leave;<br>&nbsp;&nbsp;– the restarts from fuller states with the script's `offsets_start` [own, logs/probe_wc.log]: 'w' (every coefficient's offset) keeps all 11 alarms at k_e = 20; 'wc' (and the centres' period means) removes 6, moves 1 from window period 2 to 5, and keeps 4 (s = 0.06, draw 1 from 203; draw 3 from 166; s = 0.08, draw 1 from 50; draw 2 from 32); draw 7's alarm at k_e = 10 survives both.<br>&nbsp;&nbsp;– So about half the alarms come from resetting the centres, none from dropping kernels 2 and 3's coefficient offsets, and 5 of 12 come from what 'wc' still drops (the instantaneous state against period means, the unformed kernels, the error). | • In (iii), replace the parenthesis by: "by definition the run, the full state, does not leave at a; the windows started one period off never leave, and none of the 105 neighbours still at their feature repeats the alarm; restarted with every formed kernel's coefficient offset the 11 alarms at k_e = 20 remain, with the centres at their period means too 6 vanish, 1 is delayed and 4 remain".<br>• Make it the code's: either run C8K's `offsets_start` on the alarm windows (a small part, say C8A), or store C8's centre period means and print the distance with its definition.<br>• Drop "0.038–0.111" or recompute it. |
| D3cs-3 | minor | **Mutants of this round's code** (mut/spec98.py) [own].<br>• **Killed, 10:**<br>&nbsp;&nbsp;– R02, R03, R07, R09, R10, R11, R13, R19, by the md's byte comparison in `test_round89_md_is_the_report_and_planted_edges`;<br>&nbsp;&nbsp;– R16 (a repeat only at a itself) and R18 (only after a), by round 94's planted repeat.<br>• **Surviving the full suite, 9, none changing the md or a rerun (all report code):**<br>&nbsp;&nbsp;– equivalent by construction: R05 (the far pairs keep the noise-free windows: one per configuration), R08 (far pairs in both orders: the Jaccard is symmetric, so the median is unchanged), R17 (a never-leaving neighbour counted as a repeat: impossible, since p < a);<br>&nbsp;&nbsp;– equivalent on the data: R01 (Jaccard of two empty sets 0: no such pair), R04 (the least-raise guard `is not None`), R06 (far pairs across draws: C7 has only draw 7), R12 (the tail guard `or`), R14 (C8R misses beyond two periods: all lags are 0 or −1), R15 (the label keyed on 'ke').<br>• **R04 shows a latent bug:** `if least(b) and least(q)` silently drops a window whose least ejecting raise is 0. C7 has one such window (k_e = 20, s = 0.08 from 50: 0–0.11), though not in C7S.<br>• **The authors' 7:** B01, B03, B04, B06 and B07 are killed; B02 and B05 survive, equivalent as recorded (no neighbour leaves within 5 periods of an alarm; C7's starts lie 25 apart). | • Planted edges in a test: a far pair of empty sets (R01); a C7S window with least raise 0, with the guard written `is not None` and the ratio skipped explicitly (R04); a second draw in C7 (R06); a C8R row lagging 2 periods (R14).<br>• Record R05, R08 and R17 as equivalent by construction. |
| D3cs-4 | cosmetic | **md wording.**<br>• **A misplaced parenthesis.** The C8D sentence now reads "9 of them the replay one period early, the others isolated (a lead other than one period) (the window's loop leaving in its last period, the run in the next; …)". The second parenthesis explains the early replay but follows "isolated".<br>• **Leftovers.** `byc` still drops the metric (D3co-6), and `rep_alarm` is computed but no longer printed. | • Move "the others isolated (a lead other than one period)" to the "12 isolated starts" clause (`c8_more_lines`).<br>• Add the metric to `byc`; drop or print `rep_alarm`. |
| D3cs-5 | cosmetic | **C8R in and out of sample** (D60:105–108).<br>• **Mixed samples.** The 37 runs include draw 7's four leaving runs, two of them at 80 periods back, inside the out-of-sample sentence. Out of sample the counts are 30 of 33, 30 of 33 and 27 of 31; in sample 3 of 4, 3 of 4 and 2 of 2 [own].<br>• **Term order.** "C8's three robust misses" uses C8D's term ("missed at every phase", D60:118) before it is introduced.<br>• **A near-miss.** The never-reproduced k_e = 20, s = 0.06, draw 7 is the cut's: its largest distance is 0.293 at both 20 and 40 periods back [raw]. | • Split C8R's counts by sample, in (iii) and in the md line.<br>• Write "the three leaves C8D misses at every phase".<br>• Note the 0.293. |
| D3cs-6 | cosmetic | **Draw 5 is a collision, not kernel 1's drift** (D60:119–122). Restarting the window from 150 [own, logs/probe5.json]:<br>• kernel 1's offset alone (C8): nothing leaves;<br>• kernels 1 and 2: nothing leaves;<br>• kernels 1 and 3, or all three: kernel 3 leaves at 158 and kernel 1 at 160;<br>• kernel 3 alone, kernel 1 at its fitted amplitude: also 158 and 160.<br>Kernel 1's leave there needs kernel 3's offset, not its own drift. This confirms "displaced" and is an exception to "a leave needs both the drifted coefficient and the coming noise" (D60:110–111). | Add one clause in (iii): "kernel 3's offset alone, with kernel 1 at its fit, reproduces both leaves". |
| D3cs-7 | cosmetic | **Sources** (D60:216–220).<br>• **Giesl–Rasmussen:** J. Math. Anal. Appl. 390(1), 27–46, 2012, doi 10.1016/j.jmaa.2011.12.051 [snippet]. The draft gives no year or pages.<br>• **Oljača–Ashwin–Rasmussen:** J. Dyn. Diff. Equat. 36 (2024) 2375–2411, and a correction is listed [snippet]. The forward basin with fibres B⁺(A)(t₀) = {x : d(Φ(t, t₀)x, A(t)) → 0} is in the 2022 slides [snippet]; I did not see it in the paper (the full text was not reachable).<br>• **The asymptotic–finite-time distinction.** The forward basin is an asymptotic object, while the windowed set is a finite-time one, whose counterpart is the area of attraction. | • Add the year, pages and the correction.<br>• Attribute the fibre definition to the slides.<br>• Say "the finite-time counterpart (areas of attraction) of a basin indexed by the initial time". |

## Sources

- L. Oljača, P. Ashwin, M. Rasmussen, "Measure and statistical attractors for nonautonomous dynamical systems", J. Dyn. Diff. Equat. 36 (2024) 2375–2411, with a published correction [snippet]:
  - [Exeter record](https://ore.exeter.ac.uk/repository/handle/10871/130526?show=full)
  - [Imperial record](https://spiral.imperial.ac.uk/entities/publication/857d7eb1-182f-45e9-bccc-9d09bf818308)
  - [Correction](https://spiral.imperial.ac.uk/entities/publication/784321a8-d4b9-4aae-8c99-ab91cd847da9)
- L. Oljača, 2022 MPE CDT talk slides, with the forward basin's fibres B⁺(A)(t₀) [snippet; the fetch was blocked]: [slides](https://www.mpecdt.ac.uk/wp-content/uploads/2022/05/12-Early-Career-Speaker-OLJACA-Lea.pdf)
- P. Giesl, M. Rasmussen, "Areas of attraction for nonautonomous differential equations on finite time intervals", J. Math. Anal. Appl. 390(1) (2012) 27–46 [snippet]:
  - [Sussex record](https://sro.sussex.ac.uk/id/eprint/54585/)
  - [figshare](https://sussex.figshare.com/articles/journal_contribution/Areas_of_attraction_for_nonautonomous_differential_equations_on_finite_time_intervals/23417591)
  - [arXiv 1109.1162](https://arxiv.org/pdf/1109.1162), a search result whose summary gave the volume, pages, year and doi.
- Not rechecked this round [memory or earlier rounds' snippets]: Berger–Doan–Siegmund 2008, Rasmussen 2010, Feudel–Witt–Lai–Grebogi, McDonald–Grebogi–Ott–Yorke, Schultz–Hellmann–Webster–Kurths, Ashwin–Wieczorek–Vitolo–Cox, Berglund–Gentz, Bode, Freudenberg–Looze, Hoppensteadt and Khalil.

## Lesson

- **A baseline must be conditioned like the quantity it is compared with.** A median over all cells is not the control for a statistic taken at hand-picked cells. The matched comparison (0.706) was already on the record; the pooled one (0.831) understates the effect.
- **The test of "implied by definition" applies to every piece of evidence, not only to the one a referee flagged.** The draft dismissed C8R's passes as implied by C8R's aggregate, then put in their place the run's distance below the cut, which the definition of an alarm implies just as surely. For an abrupt event, a small distance just before it carries no information.
- **A number taken from a referee's own run belongs in the code, or out of the text.** Here, recomputed, it does not even match its stated definition.
- **Test a mechanism by removing its ingredients one at a time.** Restarting with kernel 3's offset alone settles "displaced" more firmly than the period means do.

Key files (W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref98):
- reproduction: W/logs/fromraw.log, W/logs/pytest_full.log, W/logs/repro_jobs.log, W/logs/repro_runs.log;
- recounts: W/logs/own_c7s.log, W/logs/own_c8d.log, W/logs/own_c3.log, W/logs/own_preleave.log;
- probes: W/logs/probe_wc.log, W/logs/probe5.json;
- mutants: W/mut/spec98.py, W/mut/w_s1/res_s1a.jsonl, W/mut/w_s2/res_s2a.jsonl, W/mut/w_own/res_s1b.jsonl.
