# Referee report D3cw, round 102: Proposition D‴.60 as revised after round 98

## Scope

- **Commit.** I reviewed 38fcf58 on claude/vigilant-clarke-hf939b, exported with `git -C … archive 38fcf58 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref102. While I worked, HEAD moved to ba0b3d1. `git diff --stat 38fcf58 HEAD` is empty for the draft, the script, the md, the raw, the tests and the mutant spec, so everything below also holds at HEAD.
- **Read in full:**
  - review/drafts/D60.md (365 lines);
  - results/chaos/noise_resonance.md (the C7S, C8R, C8D and C8A lines line by line, the rest against the draft);
  - review/REVIEW_round98_D3cs.md and the Round 98 entry of review/REVIEW_LOG.md (lines 3850–3884);
  - tools/mutants/round98_D60_own.py and tools/mutants.py;
  - `git diff 99edc3d 38fcf58` for all of these files.
- **experiments/exp_noise_resonance.py:** every changed line, and what it calls:
  - `c7s_overlap` and `c7s_lines`;
  - the by-sample counts in C8R;
  - the C8D sentence in `c8_more_lines`, and `c8d_alarms`;
  - `offsets_start`'s `only`;
  - `C8A_PROBE`, `c8a_job`, `c8a_cases`, `c8a_lines` and `part_C8A`;
  - the wiring in `report`, `KEYS` and `main`;
  - what these call: `c_run`, `win_run`, `edge_run`, `raised`, `c_noise`, `left_periods`, `c8_job`, `c8d_job`, `c8r_job`, and `closed_loop.simulate`'s state and `x_init`.
- **tests/test_noise_resonance.py (78 tests):** the three round-98 tests line by line, and the changed asserts in the older tests.
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with the repository's venv;
  - NR_OUT=W/out holding copies of the raw and the md, LCK_BUILD=W/build, TMPDIR=W/tmp;
  - pytest `--basetemp` inside W (for tools/mutants.py through `PYTEST_ADDOPTS`, one basetemp per slot);
  - at most two heavy processes of mine at a time, each with a time limit.
- **Nothing in the repository was modified.** I ran only read-only git commands. The only thing I added inside the export is a symlink, W/repo/venv, to the venv. I stopped only my own two mutant slots after their last results, each PID checked in /proc/PID/cmdline. The ref101 processes and the modified experiments/exp_moving_width.py belong to another session.
- **Tags:** [md] the committed md; [raw] the committed raw; [own] my runs or my own code; [snippet] a web-search result only; [memory] not checked this round.

## Summary

**Design.** This round reviews the revision after round 98, which made these changes:
- the far overlap is anchored on the base windows: 0.706, the pooled 0.831 kept beside it;
- part C8A restarts the 12 isolated alarms' windows from fuller states ('w', 'wc'), with draw 5's probe by kernel subsets;
- C8R's counts are split by sample;
- the "not the state's" evidence is replaced by the definition, the repeat test and C8A;
- wording fixes and new source details.

My checks:
- **Reproduction:** regenerated the md, ran the tests, and reran 5 C8A jobs with the script's own `c8a_cases` and `c8a_job`.
- **Numbers:** recounted every new number from the raw with my own code.
- **Probe (W/logs/probe_c8a\*.json):** restarted the 12 alarm windows dropping or restoring one ingredient at a time, 28 starts per window:
  - kernel 1's centre alone, kernels 2 and 3's centres alone, and the centres mirrored;
  - a 2×2×2 factorial: formed kernels at their period means or instantaneous at pT, unformed kernels at the fit or instantaneous, plant at x_m or at x(pT);
  - kernel 1's coefficient against its centre, kernel 1 against kernels 2 and 3;
  - δ ± 0.001 and ± 0.005;
  - the full state, as a sanity check.
- **Draw 5:** restarted with all 8 kernel subsets, under the noise and without it.
- **Mutants:** 19 of my own and the authors' 20.
- **Literature:** searched the new sources.

**What holds**

- **Reproducibility is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 5e088bc624117c1dc17c0ed6738fe02b).
  - 78 tests pass in 103 s.
  - Five C8A jobs equal the raw in every field except `secs`: draw 7 at k_e = 10; draw 5's probe; and k_e = 20 at s = 0.06 (draws 2 and 3) and at s = 0.08 (draw 2). Three of them rerun C's loop. Together they cover 8 of the 12 alarm windows and all 5 probe windows.
  - My probe reruns all 6 alarm runs' loops: kernel 1's period means are bit-identical to C8's `w`, and the leave periods match.
  - The reduced window equals C8D's (leave period and largest distance) in 12 of 12, and my 'w' and 'wc' starts equal C8A's raw in 12 of 12.
- **Every number new to the draft is printed by the md** [md] and agrees with my recounts [own]:
  - the anchored overlap 0.7056 over 84 pairs, the least raises' ratio median ×1.3609, the pooled 0.8313 over 420 pairs;
  - C8R out of sample 30/33, 30/33, 27/31, in sample 3/4, 3/4, 2/2;
  - C8A: 'w' keeps 12; 'wc' removes 6, delays 1 (draw 2 at s = 0.06 from 58: period 2 → 5, still 81 periods before its run's leave) and keeps 5;
  - draw 5's five subsets.
- **"Like for like" is right:**
  - The near pairs include C9's draw-3 window, which has no anchored partner (C7 holds only draw 7). The near median over the 12 C7 windows is 0.974.
  - The far pairs are not a decay with distance: by |Δshift| their medians are 0.698, 0.632, 0.768, 0.830, 0.572, 0.745, 0.596 at 25 … 175 [own].
- **"No alarm comes from dropping kernels 2 and 3's coefficient offsets" holds whichever way the ingredients are added** [own]. With the centres at their period means, removing those offsets ('c') changes one outcome (draw 2 from 58), and in the direction that removes the alarm.
- **Draw 5** [own]:
  - Under the noise kernel 1 leaves at 160 (and kernel 3 at 158) if and only if kernel 3's offset is applied, in all 8 subsets, with or without the centres.
  - Without noise, kernel 3's offset alone (+0.124, beyond its noise-free edge 0.075) ejects kernel 3 in period 1 but leaves kernel 1 in place (its largest distance below 10⁻⁴). The exception the draft states is precise: the leave does not need kernel 1's drifted coefficient, but it does need the coming noise.

**What does not hold (no major issue)**

1. **D3cw-1 (minor). C8A's remainder is left unresolved, and "the centres" is one centre.** My one-ingredient restarts [own]:
   - **The unformed kernels and the plant's error never change an outcome:** 12 of 12 windows, in all four combinations.
   - **The 6 that 'wc' removes are removed by kernel 1's centre alone:**
     - 'w' plus kernel 1's centre gives 'wc' in 12 of 12; 'w' plus kernels 2 and 3's centres gives 'w' in 12 of 12;
     - it is not a generic sensitivity: with the centres mirrored inward all 12 stay, and δ ± 0.005 keeps all 11 alarms at k_e = 20;
     - kernel 1's period-mean centre sits 0.0325–0.0713 outward in the removed ones and at most 0.025 in the kept ones.
   - **The rest go with the formed kernels' instantaneous state:**
     - at k_e = 20, kernel 1's instantaneous centre (0.040–0.098 outward) removes 4 of the 5, and from the reduced start it removes 10 of 11;
     - draw 7's alarm at k_e = 10 is removed by kernel 1's instantaneous coefficient, and is marginal: largest distance 0.3104 against the cut 0.3, in the window's last period, gone at δ − 0.001;
     - draw 2 from 58 at s = 0.06 is removed by kernels 2 and 3's instantaneous state.
   - **Two alarms do not involve the drifted coefficient:** at s = 0.08, draw 2 from 6 and from 18, δ is 0.0096 and 0.032.
   - **The draft's list of three candidates** is true by determinism but uninformative, and two of its three items never matter.
2. **D3cw-2 (minor). The suite is blind to C8A's computation.**
   - 7 of my mutants survive the full suite and change the published C8A results when the part is regenerated: T01–T05 (the rerun's draw and sign, the windows' draw, level and metric), T07 (W and C swapped) and T08 ('reduced' set to the absolute period).
   - Regenerated with `--only C8A`: under T03 'w' "removed 9"; under T08 'w' "earlier 12".
   - The md test only reads the raw, and the wiring test's mocks ignore s, opts, the metric and the start.
3. **D3cw-3 (minor). C8A's conclusions are fixed strings.**
   - "So no alarm comes from dropping kernels 2 and 3's coefficient offsets; those 'wc' removes come from resetting the centres …" and "— its leave there needs kernel 3's offset, not its own" are printed whatever the data say.
   - They are printed on contrary synthetic input, and under mutant T03 next to "'w' … removed 9" [own].
4. **D3cw-4 (cosmetic).** The md attaches the wrong reason to C8R's passes: "(implied: by an alarm's definition the run, the full state, does not leave there)". C8R's loop is the reduced loop; its passes follow from C8R's aggregate (none earlier), as D60 correctly says.
5. **D3cw-5 (cosmetic). Latent bugs, equivalent on the data:**
   - `c8a_job`'s rerun drops the metric: `c_run(m, ke, s, dict(opts))`, whereas `c8_job` passes `metric=mt`;
   - `ancr` silently drops a zero least raise (R04's pattern again).
6. **D3cw-6 to D3cw-8 (cosmetic):**
   - wording: the md's C8A parenthesis follows "removed 6" but describes the delayed alarm; the draft's broken dash pair and double space at lines 123–124; C8A is missing from Status (iii) and from the script's docstring and usage line (C7S too); "with kernel 2's" is ambiguous;
   - draw 5's noise dependence, not stated;
   - "(an area of attraction)" read as an identity rather than an analogy.

## What I checked and how

All paths are under W.

- **Reproduction:**
  - logs/fromraw.log, then `cmp` against the committed md: identical.
  - logs/pytest_full.log: 78 passed, 1 old RuntimeWarning (B3).
  - scripts/repro_c8a.py → logs/repro_c8a.log: `c8a_job` on 5 of `c8a_cases`' 7 jobs, compared with raw['C8A'] with `secs` excluded: all identical.
- **Own recounts from the raw (numpy only, `python -I`):** scripts/own_recount.py → logs/own_recount.log:
  - near Jaccard and least ratios; anchored far pairs (median, ratio, by distance); pooled pairs;
  - C8R by sample and the runs it never reproduces;
  - C8A's classes, leads and run leaves;
  - C3 against C8 for draws 7 and 5: kernel-1 period means bit-identical.
- **Probes:**
  - scripts/probe_c8a.py → logs/probe_c8a.json: the 12 alarm windows from 'red', 'w', 'wc', 'c', 'c1', 'w+c1', 'w+c23', 'wcm' (centres mirrored), the 8-cell factorial (formed: means or instantaneous; unformed: fit or instantaneous; plant: x_m or x(pT)) and δ ± 0.001, ± 0.005. The state at pT comes from `c_run`'s outputs (8 a period). The full-state restart never leaves, and its kernel-1 period-mean distance is within 0.014–0.043 of the run's (output grids of 64 against 8 a period).
  - scripts/probe_c8a2.py → logs/probe_c8a2.json: formed coefficients and centres instantaneous separately; kernel 1 against kernels 2 and 3; kernel 1's instantaneous centre or coefficient alone, on 'wc' and on the reduced start.
  - scripts/probe5.py → logs/probe5.json: draw 5 from 150, all 8 subsets, with and without the centres, under the noise and without it, all kernels' leave periods. Each window under the noise is cross-checked against `win_run`.
  - logs/c8a_lines_hardcoded.log: `c8a_lines` on contrary synthetic inputs.
- **Mutants:**
  - mut/spec102.py: 19 and a control.
  - Stage 1 (`-k "round98 or round89 or round94 or c8_more_lines or c8_windows or round85"`): mut/w_s1/res_s1a.jsonl.
  - Survivors on the full suite: mut/w_s2/res_s2a.jsonl and mut/w_s2r/res_s2r.jsonl, the latter from the reversed spec mut/spec102r.py.
  - The authors' tools/mutants/round98_D60_own.py: mut/w_own/res_own1.jsonl (stage 1) and mut/w_own2/res_own2.jsonl (S01 and S04 on the full suite).
  - Regeneration under a mutant: scripts/demo_mut.sh → mutdemo/T03/c8a_line.txt and mutdemo/T08/c8a_line.txt.
- **Literature:** web searches on Oljača–Ashwin–Rasmussen, Giesl–Rasmussen and the doi. Fetching link.springer.com failed because the host is not in this session's allowed network. It has to be allowed in the session's network settings before the full text can be read.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii) and their proofs (9–49, 162–187) | unchanged since round 98, not reverified | [md] |
| The noise-free edge, power law, C4, s*, bands, cuts, C7L (51–96) | unchanged | [md] |
| C7S: ×4.08, median ×1.019; near 0.975 against the anchored 0.706, ×1.36, "the like-for-like baseline"; pooled 0.831 (97–101, 338–339) | verified. Near over the 12 C7 windows 0.974. The far pairs show no trend with distance (0.57–0.83) | [md], [own] |
| C8R 33/37, 33/37, 29/33; by sample 30/33, 30/33, 27/31 and 3/4, 3/4, 2/2; "the three leaves C8D misses at every phase" (106–109, 348–349) | verified | [md], [own] |
| Draw 5: kernel 3's offset alone → kernel 1 at 160; kernel 1's alone or with kernel 2's → never; "an exception to a leave needing kernel 1's drifted coefficient" (120–124) | verified. Holds over all 8 subsets. The collision needs the coming noise | [raw], [own]; D3cw-7 |
| 21 alarms = 9 early + 12 isolated, leads 11–86, periods 2 or 10 (127–131) | verified | [md], [own] |
| "By an alarm's definition the run … does not leave there"; C8R "passes all 9 … (implied by its aggregate)" (131–132) | correct. The md gives the wrong reason | D3cw-4 |
| C8A: 'w' keeps 12; 'wc' removes 6, delays 1, keeps 5 (133–135, 353–354) | verified, bit-identically | [raw], [own] |
| "No alarm comes from dropping kernels 2 and 3's coefficient offsets" (135–136) | verified, whichever way the ingredients are added | [own] |
| "About half come from resetting the centres" (136) | true; precisely kernel 1's centre, and its outward offset | [own]; D3cw-1 |
| "The rest from what 'wc' still drops — the instantaneous state …, the unformed kernels, the error" (136–137) | true by determinism, unresolved. Only the formed kernels' instantaneous state matters | [own]; D3cw-1 |
| Status (iii) (198–206) | accurate, but C8A is missing | D3cw-6 |
| Sources (224–229) | year, pages and correction confirmed. "(an area of attraction)" loose | [snippet]; D3cw-8 |
| New: "(out of sample 30 of 33, 30 of 33 and 27 of 31)" (236) | verified | [md] |
| Numerics C8D, C8A (350–355) | verified | [md], [raw] |

## Status of the round-98 issues

| D3cs id | status | what remains |
|---|---|---|
| D3cs-1 | fixed | `c7s_overlap` prints 0.706, 84 pairs and ×1.36, with the pooled 0.831 beside them; D60 and the md say "like for like"; the test asserts the md string [own: passes]. Remaining: the near statistic includes C9's draw-3 window, which has no far partner (0.974 without it), a negligible effect |
| D3cs-2 | fixed in substance | "0.038–0.111" dropped; the definition is stated; C8A added and reproducible. Remaining: the remainder is unresolved (D3cw-1); the md's "implied" clause gives the wrong reason (D3cw-4); the conclusions are hard-coded (D3cw-3) |
| D3cs-3 | fixed | R01, R04, R06 and R14 are planted and kill their mutants; R05, R08 and R17 are recorded as equivalent. Of the authors' 20, 18 are killed; S01 (on the data) and S04 (by construction) survive the full suite, as recorded [own]. Remaining: `ancr` drops a zero least raise silently (D3cw-5) |
| D3cs-4 | fixed | The C8D sentence, `byc` and `rep_alarm` are done. The same misplacement recurs in C8A's line (D3cw-6) |
| D3cs-5 | fixed | Split by sample; the term order is fixed. The near-miss 0.293 was not added, with a reason recorded (acceptable) |
| D3cs-6 | fixed | `C8A_PROBE`; confirmed and extended (D3cw-7) |
| D3cs-7 | fixed | Year, pages, the correction, the slides' attribution and the finite-time sentence are in. Wording remains (D3cw-8) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cw-1 | minor | **C8A's remainder unresolved; "the centres" is kernel 1's centre** (D60:133–137; md C8A).<br>The 12 alarm windows restarted one ingredient at a time [own, logs/probe_c8a.json, logs/probe_c8a2.json]:<br>• **The unformed kernels and the plant's error never matter.** Formed kernels at their period means give 'wc''s outcome in 12 of 12, in all four unformed × error combinations; instantaneous formed kernels give no leave in 12 of 12, again in all four.<br>• **Kernel 1's centre alone does what 'wc' does.** 'w' plus kernel 1's centre equals 'wc' in 12 of 12; 'w' plus kernels 2 and 3's centres equals 'w' in 12 of 12.<br>• **Not a generic sensitivity.** The centres mirrored inward keep 12 of 12; δ ± 0.005 keeps all 11 alarms at k_e = 20.<br>• **The removed and the kept separate cleanly.** Kernel 1's period-mean offset is −0.0325 … −0.0713 (outward) in the 6 removed and −0.025 … +0.003 in the 5 kept.<br>• **The remainder goes with the instantaneous formed state:**<br>&nbsp;&nbsp;– kernel 1's instantaneous centre (−0.027 … −0.099) removes 4 of the 5 at k_e = 20, and 10 of 11 from the reduced start;<br>&nbsp;&nbsp;– draw 7's alarm (k_e = 10) goes with kernel 1's instantaneous coefficient. It is marginal: 0.3104 against the cut 0.3 in the window's last period, removed at δ − 0.001;<br>&nbsp;&nbsp;– draw 2 from 58 (s = 0.06) goes with kernels 2 and 3's instantaneous state.<br>• **Two alarms are ejections from nearly the exact fit** (s = 0.08, draw 2 from 6 and 18: δ = 0.0096 and 0.032). | • Add to C8A, per window (about 1 s each): 'w' + kernel 1's centre, 'w' + kernels 2 and 3's centres, and the factorial with the state at pT (rerun C's loop for the C3 cases too, since C3 stores only period means). Optionally δ ± 0.001.<br>• In D60:135–137 and the md: "those 6 by kernel 1's centre alone (its period mean 0.033–0.071 outward, the kept at most 0.025; mirrored inward none); the rest with the formed kernels' instantaneous state — 4 with kernel 1's instantaneous centre, draw 7's (a crossing of 0.310 in the window's last period, gone at δ − 0.001) with its coefficient, draw 2 from 58 with kernels 2 and 3's —, never with the unformed kernels or the error". |
| D3cw-2 | minor | **The tests do not protect C8A's computation** (tests/test_noise_resonance.py: `test_round98_c8a_job_wiring`, `…_the_fuller_restarts_and_draw_5`).<br>• **Survivors on the full suite, 14 of my 19:**<br>&nbsp;&nbsp;– **would change the published C8A results on regeneration, 7:** T01 (the rerun under the default draw), T02 (the rerun under −n), T03 (the windows under draw 7's noise), T04 (at s/2), T05 (under WFR), T07 (the rerun's W and C swapped), T08 ('reduced' set to the absolute period);<br>&nbsp;&nbsp;– **regenerated with `--only C8A` [own]:** T03 gives "'w' delayed 2, kept 1, removed 9"; T08 gives "'w' earlier 12";<br>&nbsp;&nbsp;– **equivalent on the data:** T09, T10, T12, T13, T21;<br>&nbsp;&nbsp;– **equivalent by construction:** T14, T20.<br>• **Killed:** T06, T11, T18, T19, T22.<br>• **Why they survive:** the md test reads the raw, and the mocks record only (p, d, start), never calling the start closure. The wiring test's docstring claims "under the run's noise from the window's start", which no assert checks. | • A real rerun test: `c8a_job` on draw 7's job (1 s) and the probe job (2 s), compared with raw['C8A'] with `secs` excluded. Plus one rerun job (k_e = 20, s = 0.08, draw 2: 12 s), or a `c_run` mock asserting (m, ke, s, seed, metric).<br>• In the wiring test, record and assert (m, ke, mt, s, opts) in the `win_run` mock.<br>• Assert `c8a_cases`' reduced periods against C8D's window leave periods. |
| D3cw-3 | minor | **Hard-coded conclusions** (`c8a_lines`).<br>• "So no alarm comes from dropping kernels 2 and 3's coefficient offsets; those 'wc' removes come from resetting the centres, and those it keeps from what it still drops (…)" and "— its leave there needs kernel 3's offset, not its own" are fixed strings.<br>• On synthetic input where 'w' removes both alarms, or where the kernel-3 window stays, the same sentences print [own, logs/c8a_lines_hardcoded.log]. Under T03 they print beside "'w' … removed 9" [own, mutdemo/T03].<br>• "Nothing is edited by hand" holds only for the numbers. | • Print each clause only when the counts imply it: 'w' keeps all; 'wc''s removals among 'w''s kept; every probe row with kernel 3 leaves and every row without it does not. Otherwise print the counts alone.<br>• Add contrary synthetic cases to `test_round98_c8a_the_fuller_restarts_and_draw_5`. |
| D3cw-4 | cosmetic | **The md's C8R clause gives the wrong reason** (md C8D sentence; `c8_more_lines`).<br>• "(implied: by an alarm's definition the run, the full state, does not leave there)" — C8R's loop is the reduced loop started earlier, not the run. Its passes follow from C8R's aggregate (0 earlier), as D60:132 says.<br>• The clause is printed whenever C8R covers an alarm, even if one failed. | • Print "(implied by C8R's aggregate: none of its loops leaves earlier than one period before its run)", only when all covered alarms pass.<br>• State the definition separately.<br>• Update the strings in `test_round89_alarms_by_kind_and_the_survivors` and `test_round98_…`. |
| D3cw-5 | cosmetic | **Latent bugs, equivalent on the data** (`c8a_job`, `c7s_overlap`).<br>• The rerun calls `c_run(m, ke, s, dict(opts))` without `metric=mt` (`c8_job` passes it). A WFR alarm would get the paradigm's period means. All 12 alarms are under the paradigm's metric; my reruns with the metric are bit-identical to C8.<br>• `ancr` drops a pair with a least raise of 0 without counting it, while the md prints "84 pairs; their least raises' ratio …". On the data 84 of 84. | • `c_run(m, ke, s, dict(opts, metric=mt))`.<br>• Count zero-raise anchored pairs apart (like `lr0`) and print `len(ancr)`.<br>• Add a planted test for each. |
| D3cw-6 | cosmetic | **Wording.**<br>• The md's C8A line: "'wc' delayed 1, kept 5, removed 6 (…draw 2 from 58: period 2 → 5)". `desc` appends the parenthesis after the sorted classes, so it follows "removed 6" (as in D3cs-4).<br>• D60:123–124: "— an exception to a leave needing kernel 1's drifted coefficient;  leaving first alone does not explain it: … —," has a broken dash pair and a double space.<br>• Status (iii) (D60:204–206) omits C8A. The script's docstring and usage line omit C8A and C7S.<br>• "With kernel 2's" (D60:123) is ambiguous: it was kernels 1 and 2 (kernel 2 alone also never leaves [own]). | • Attach the parenthesis to its class.<br>• "(C8A; referee round 98), an exception to a leave needing kernel 1's drifted coefficient — kernel 3 leaving first does not alone explain it: … —".<br>• Add C8A to Status (iii), the docstring and the usage line.<br>• "With kernels 1 and 2's". |
| D3cw-7 | cosmetic | **Draw 5's collision needs the noise** [own, logs/probe5.json].<br>• The offsets at 150 are +0.203, +0.057, +0.124.<br>• Without noise: kernel 3's offset ejects kernel 3 in period 1 and leaves kernel 1 in place (largest distance < 10⁻⁴); kernel 1's own +0.203 ejects kernel 1 in period 2.<br>• Under the noise: kernel 1 leaves at 160 if and only if kernel 3's offset is applied, in all 8 subsets, with or without the centres. | • Add to (iii): "without the noise kernel 3's offset (beyond its noise-free edge 0.075) ejects kernel 3 in period 1 and kernel 1 stays: the collision needs the coming noise".<br>• Extend C8A_PROBE to all subsets and a noise-free row. |
| D3cw-8 | cosmetic | **Sources** (D60:224–229).<br>• **Confirmed** [snippet]: Oljača–Ashwin–Rasmussen, JDDE 36 (2024) 2375–2411, with a correction; Giesl–Rasmussen, JMAA 390(1) (2012) 27–46 (via arXiv 1109.1162's reference list).<br>• **Giesl–Rasmussen's areas of attraction** are defined by finite-time (exponential) attractivity, with a time-varying Riemannian metric in the criterion [snippet]. C7's set is a threshold set of starting raises.<br>• **The fibre definition** in the paper is still unverified: link.springer.com is not reachable from this session. | • "The windowed set a threshold analogue of an area of attraction" rather than "(an area of attraction)". |

## Sources

- L. Oljača, P. Ashwin, M. Rasmussen, "Measure and statistical attractors for nonautonomous dynamical systems", J. Dyn. Diff. Equat. 36 (2024) 2375–2411, doi 10.1007/s10884-022-10196-5, with a correction [snippet]:
  - [Exeter record](https://ore.exeter.ac.uk/repository/handle/10871/130526?show=full)
  - [Imperial record](https://spiral.imperial.ac.uk/entities/publication/857d7eb1-182f-45e9-bccc-9d09bf818308)
  - [Correction](https://spiral.imperial.ac.uk/entities/publication/784321a8-d4b9-4aae-8c99-ab91cd847da9)
  - The Springer page could not be fetched: link.springer.com is not in this session's allowed network.
- L. Oljača, 2022 MPE CDT slides (seen only in the search listing) [snippet]: [slides](https://www.mpecdt.ac.uk/wp-content/uploads/2022/05/12-Early-Career-Speaker-OLJACA-Lea.pdf)
- P. Giesl, M. Rasmussen, "Areas of attraction for nonautonomous differential equations on finite time intervals", J. Math. Anal. Appl. 390(1) (2012) 27–46 [snippet]:
  - [Sussex record](https://sro.sussex.ac.uk/id/eprint/54585/)
  - [figshare](https://sussex.figshare.com/articles/journal_contribution/Areas_of_attraction_for_nonautonomous_differential_equations_on_finite_time_intervals/23417591)
  - [arXiv 1109.1162](https://arxiv.org/pdf/1109.1162): its reference list gives the volume, pages and year.
- M. Rasmussen, "Finite-time attractivity and bifurcation for nonautonomous differential equations", Differ. Equ. Dyn. Syst. 18 (2010) [snippet, via his publications page](https://www.ma.imperial.ac.uk/~mrasmuss/main.html)
- Not rechecked this round [memory or earlier rounds' snippets]: Berger–Doan–Siegmund 2008, Feudel–Witt–Lai–Grebogi, McDonald–Grebogi–Ott–Yorke, Schultz–Hellmann–Webster–Kurths, Ashwin–Wieczorek–Vitolo–Cox, Berglund–Gentz, Bode, Freudenberg–Looze, Hoppensteadt, Khalil.

## Lesson

- **Finish a decomposition by restoring each remaining ingredient alone.** "The rest from what the restart still drops" is true by determinism and says nothing. The one-at-a-time restarts cost a second each, and here they show that two of the three listed candidates never matter, and that "the centres" is one kernel's outward offset.
- **A conclusion the generator prints must be computed from the data it summarises.** Otherwise it is hand-written text that survives any change of the data, as mutant T03 shows.
- **A job tested only through mocks that ignore its arguments is untested.** When the real job takes a second, run it against the raw.

Key files (W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref102):
- reproduction: W/logs/fromraw.log, W/logs/pytest_full.log, W/logs/repro_c8a.log, W/logs/own_recount.log;
- probes: W/logs/probe_c8a.json, W/logs/probe_c8a2.json, W/logs/probe5.json, W/logs/c8a_lines_hardcoded.log;
- mutants: W/mut/spec102.py, W/mut/w_s1/res_s1a.jsonl, W/mut/w_s2/res_s2a.jsonl, W/mut/w_s2r/res_s2r.jsonl, W/mut/w_own/res_own1.jsonl, W/mut/w_own2/res_own2.jsonl, W/mutdemo/T03/c8a_line.txt, W/mutdemo/T08/c8a_line.txt.
