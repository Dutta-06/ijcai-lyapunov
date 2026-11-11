# Referee report D3co, round 94: Proposition D‴.60 as revised after round 89

## Scope

- **Commit:** 89e9e28 on branch claude/vigilant-clarke-hf939b, exported with `git -C … archive 89e9e28 | tar -x -C W/repo`.
  - W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref94. It holds repo, scripts, logs, out, build, pyt and mut.
  - While I worked, another session moved HEAD to 5a171c7. `git diff --stat 89e9e28 HEAD` is empty for the D‴.60 files: the draft, the script, the md, the raw and the tests. Everything below holds for both commits.
- **Files read in full:**
  - review/drafts/D60.md (342 lines);
  - results/chaos/noise_resonance.md (708 lines);
  - tests/test_noise_resonance.py (74 tests; the three round-89 tests line by line);
  - review/REVIEW_round89_D3cj.md and the Round 89 entry of review/REVIEW_LOG.md (lines 3566–3606);
  - tools/mutants/round89_D60_own.py, tools/mutants/round89_D60_referee.py and tools/mutants.py;
  - `git diff 5e42eee 89e9e28` for these files.
- **experiments/exp_noise_resonance.py (2772 lines):**
  - all code new or changed since 5e42eee, line by line: `C7S_CASES` and `C7S_D`, `c7s_job`, `part_C7S`, `c7s_lines`, `c8d_alarms`, the C8D block's split, `c7l_lines`' superset and persistence, C8K's offsets, `c9a_lines`' run count, `c7_gap_lines`' guard, the `edge_lines` wording, and `report`, `KEYS` and `main`;
  - the code these call: `c7_job`, `edge_run`, `c_noise`, `win_run`, `c8_delta`, `c8d_job`, `c8r_job`, `offsets_start`, `gap_devs`, `left_periods`, `c_run`, `simulate`.
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` with the repository's venv;
  - `NR_OUT=W/out`, holding a copy of the raw and the md; `LCK_BUILD=W/build`; pytest with `--basetemp` inside W;
  - at most two heavy processes of mine at a time, apart from a few seconds when two short reruns overlapped the start of the mutant slots.
- **Nothing in the repository was modified.** I ran only read-only git commands. I stopped only one process, my own Monitor tail, through TaskStop. Nothing of mine is left running.
- **Tags:** [md] the committed md; [raw] the committed raw; [own] my runs or my own code; [snippet] a web-search result only; [memory] not verified this round.

## Summary

**Design.** This round reviews the revision after round 89. That round added:
- C7S: the gap windows started one period earlier and later;
- the split of the dense early alarms into two kinds;
- the C7L superset and persistence counts;
- the cause of draw 5's miss and its offsets;
- the s* wording, and new sources.

My checks:
- **Reproduction:** regenerated the md from the raw, ran the tests, and reran 13 jobs or windows with the script's own job functions (six C7S scans among them).
- **Numbers:** checked every number in the draft against the md, and recounted every new aggregate from the raw with my own code.
- **Logic of the new classifications:** three probes.
  1. The run's own state at each isolated alarm, and the alarm windows restarted from fuller states.
  2. How close the sets of neighbouring starts are, compared with starts 25 or more periods apart.
  3. The noise slid by quarter periods between two neighbouring starts.
- **Mutants:** 33 new single-point mutants of round 89's code, plus the authors' 19 and 10 of the round-89 referee's.
- **Literature:** searched the new sources.

**What holds**

- **Reproducibility is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 2aa8882bcc3beaf788c0eb573842c5d6).
  - The 74 tests pass in 95 s.
  - Every rerun equals the raw in every field except `secs`: six C7S scans of 151 raises each, two C8R jobs, five C8D windows (among them three alarm windows), and the five C runs at k_e = 20 that I reran for a probe.
- **Every new number of (iii), the Status and the Numerics is printed by the md** [own, numcheck]. The exceptions are roundings (8.1%, 0.97–1.03, 1.25, 0.89, 1.22), the derived "from 99" and "from 74", and citation years.
- **My recounts agree** [own]:
  - C7: 18 windows with more than one band, 12 of them with a gap of the basin.
  - C7S: the per-window table, and the tally 6 / 7 / 0. The tally is also 6 / 7 / 0 when "survives" means that a neighbour's basin gap overlaps one of the base window's.
  - C7L: 14 changed, 13 with more than one band, 10 with a gap of the basin, 23 supersets; 10 gaps of the basin kept and the same 3 closed.
  - The dense alarms: 21 = 9 + 12, by configuration, leads 11–86, leaving in the window's period 2 or 10; 106 neighbours, 1 "leaves by then"; C8R covers 9 and passes 9.
  - C8K: 0.015–0.032 and 0.032–0.065 at s = 0.02, 0.040–0.078 and 0.089–0.152 over 10 of 11 draws at s = 0.03, draw 11's +0.144, draw 5's +0.057 and +0.124.
  - C9A: 1554 windows of 12 runs.
- **"The isolated alarms are the reduced start's, not the state's" survives a direct test** [own].
  - At each isolated alarm's leave period, the run's own kernel 1 is 0.038–0.111 from its feature, far from the cut of 0.3.
  - The windows one period before and after each isolated alarm never leave.
- **Round 89's issues are addressed.** D3cj-2, -6, -8 and -9 are fixed; D3cj-1, -3, -4, -5, -7 and -10 are fixed in substance, with what remains listed below. The seven round-89 survivors (R03, R05, R11, R16, R21, R22, R31) are now killed, and so are all 19 of the authors' own round-89 mutants [own].

**What does not hold (no major issue)**

1. **D3co-1 (minor).** "The 20-period set contains the 10-period one in all 23 — the set grows with the window" is true by construction, not measured.
   - The 20-period window is the 10-period window continued. Its first 10 periods are bit-identical (max |Δc| = 0 [own]), and all 23 windows keep every short leaver's leave period [raw].
   - The round-89 referee (D3cj-4) also offered this as a finding. That was a mistake.
2. **D3co-2 (minor).** The dependence on the start is stated only for gaps, but it moves the least ejecting raise too, and it varies on scales much shorter than the one-period step.
   - Between neighbouring starts the least ejecting raise changes ×4.1 (s = 0.03: 0.1225 from period 50, 0.5 from 49) and ×2.4 (0.2 from 175, 0.48 from 176) [raw].
   - With the noise slid by quarter periods between 49 and 50, the least raise runs 0.5, 0.1275, 0.58, 0.21, 0.1225 [own].
   - Neighbouring starts nonetheless share most of the set: median Jaccard 0.975, against 0.706 for starts 25 or more periods apart [own].
3. **D3co-3 (minor).** The alarms' classification and its evidence.
   - "Isolated" is defined by exclusion (any alarm whose lead is not 1) and is never tested as isolation.
   - The "1 leaves by then" is window 50's own alarm leaving at period 52, not a repeat of window 58's leave at 60. The Numerics' "but one alarm another" is therefore wrong.
   - C8R's 9 passes follow from its aggregate result, so they add no new evidence. The direct evidence (the run's distance; fuller starts) is not used.
4. **D3co-4 (minor).** Draw 5's cause: "jumps into the vacated first feature two periods before kernel 1 leaves" contradicts the period means [raw C3]. The feature is still occupied when kernel 3 arrives, and kernel 1 is displaced.
5. **D3co-5 (minor).** Of my 33 mutants, 10 survive the full suite and none changes the md.
   - Q30 changes a rerun and is seen by neither the tests nor the md: `main` would run C7L's part for `--only C7S`.
   - The other 9 are equivalent on the data, one of them (Q08) by construction.
6. **D3co-6 to D3co-8 (cosmetic).**
   - The md's C7L header still misdescribes the selection, plus some wording.
   - The four runs C8R never reproduces are not named (the open part of D3cj-10), and New compresses C8R's counts.
   - Reference details for the new sources.

## What I checked and how

All paths are under W.

- **Reproduction:**
  - logs/fromraw1.log with `cmp` against the committed md; logs/pytest_full.log.
  - scripts/repro.py → logs/repro_a.log and logs/repro_b.log. Field-by-field comparison with `secs` excluded, using the script's `c7s_job`, `c8r_job` and `win_run`:
    - C7S jobs 2, 4, 14, 17, 20, 25 (s = 0.03 from 49; s = 0.01 from 74; s = 0.03 from 174; k_e = 20, s = 0.04 from 26; m = 4, s = 0.08 from 99; draw 3 from 151);
    - C8R for k_e = 10, s = 0.02, draw 6 and for k_e = 20, s = 0.08, draw 2;
    - C8D windows (k_e = 20, s = 0.08, draw 2 from 6; k_e = 10, s = 0.02, draw 7 from 107; k_e = 20, s = 0.06, draw 3 from 210; m = 4, s = 0.08, draw 3 from 146; k_e = 20, s = 0.08, draw 1 from 58).
  - `repro.py PREFIX`: one raise integrated over 10 and over 20 periods, the first 10 compared bit for bit.
  - scripts/numcheck.py: every number on D60:51–147 and 176–342 looked up in the md.
- **Own recounts from the raw:**
  - scripts/own_c7.py → logs/own_c7.log: bands, gaps of the basin, the C7S tally, C7L changes, supersets, prefix leave periods, persistence, and the window period of each leave;
  - scripts/own_c8d.py → logs/own_c8d.log: each alarm's lead, its neighbours, the neighbours one period before and after, the C8R rows started before it, and for draw 7 the run's distance from C3;
  - scripts/own_c8k.py → logs/own_c8k.log: the offsets per draw;
  - scripts/own_draw5.py → logs/own_draw5.log: the period means around the leaves of draws 5, 8 and 11 at s = 0.03, with their C8 and C8K windows;
  - scripts/own_c8r.py: C8R's lags and the runs it never reproduces.
- **Probes:**
  - scripts/own_alarm_probe.py → logs/own_alarm_probe.log. It reran the five C runs at k_e = 20 that carry the 11 isolated alarms. For each alarm it gives the run's distance at the alarm's period, and the alarm window restarted from C8K's 'w' (every coefficient's offset) and 'wc' (offsets and centres).
  - scripts/own_c7s_overlap.py → logs/own_c7s_overlap.log: Jaccard similarity and the least raise, for neighbouring starts and for starts 25 or more periods apart; overlap of the gaps of the basin.
  - scripts/own_fracshift.py → logs/own_fracshift.log and fracshift_49.25_49.5_49.75.json: `c7_job` at shifts 49.25, 49.5 and 49.75 (s = 0.03, draw 7).
- **Mutants:**
  - mut/mutants94.py is tools/mutants.py with REPO set to W/repo, the venv's python and `--basetemp`. mut/spec94.py holds 33 mutants and a control.
  - Stage 1 ran the tests that read this code (`-k "round89 or c8_more_lines or c7l or c7_gap or c9a or case_lists or round85 or c8_windows"`, 11 tests, 3 s). Stage 2 ran the survivors against the full suite (each mutant also `--from-raw` and `cmp`).
  - Results: mut/w_spec94/res_*.jsonl, mut/w_round89_D60_own/res_*.jsonl, mut/w_full/res_*.jsonl, mut/w_ref/res_*.jsonl. The slot trees were deleted afterwards.
- **Literature:** web searches on Oljača–Ashwin–Rasmussen, Berger–Doan–Siegmund 2008, Rasmussen 2010, and finite-time areas of attraction.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii) and their proofs (9–49, 149–172) | unchanged; re-derived: |S|² − 1, |S(√K)|, the four bound values 1.41/1.38/1.24/1.11, the Floquet geometric sum, τ = −1/ln ρ | [own] |
| Noise-free edge, anatomy, power law, outside the window (51–72) | unchanged; numbers printed | [md] |
| C4 like for like; s* "cannot tell the draw's ⟨w₂⟩ from the configuration's", dash opened (72–83) | verified; D3cj-2 fixed | [md] |
| Bands, cuts, 12 basin / 6 cut gaps (83–90) | verified | [own] |
| C7L: 14 of 23, 13 multi-band, 10 basin; "C7's 18 … four controls and C9's draw 3" (90–92) | verified | [own] |
| "the 20-period set contains the 10-period one in all 23 — the set grows" (92) | true by construction | D3co-1 |
| 10 survive, 3 close (93–94) | verified, the same three | [own] |
| C7S: 6 / 7 / none; 99 → 0.56–1; 74 → 0.085–1; "depends on the start" (94–96) | verified; understated (start dependence also in the least ejecting raise; sub-period scale) | D3co-2 |
| Least ejecting raise along the draw (96–98) | numbers verified; 25-period sampling does not resolve its variation | D3co-2 |
| C8, C8R, C8N, the δ sweep, the noise-free edge, the six misses (98–111) | unchanged; C8R recounted (lags 0: 23/23/20, −1: 10/10/9) | [own] |
| C8D phases, the three robust misses (111–114) | unchanged | [md] |
| Draw 5: "the reduction's"; draw 11's contrast; C8K recovers it; offsets (114–118) | the reduction and the recovery verified; "vacated" wrong | D3co-4 |
| Dense alarms: 9 one period early (5, 2, 1, 1), 12 isolated (1, 5, 6), leads 11–86, period 2 or 10 (118–122) | verified | [own] |
| "106 …, 1 leaves by then, itself one of the isolated alarms" (122–123) | literally true; not a repeat | D3co-3 |
| "C8R covers 9 of the 12 and passes all 9: … the reduced start's, not the state's" (123–124) | conclusion supported by my direct test; the cited evidence is indirect | D3co-3 |
| C9, C9A (124–134) | unchanged; "12 runs" verified | [own] |
| Single tones, E only at 13.1 rad/s, 10 rad/s at 0.06 and 0.08 (134–147) | verified; D3cj-8 fixed | [md] |
| Status (iii) with C7S (185–193) | accurate | — |
| Classical: stochastic counterpart; non-autonomous basins (205–212) | apt; sources exist | D3co-8 [snippet] |
| New: "reproduces 33 of 37 leaves within one period, from up to 80 periods before" (218–219) | the denominators conflated | D3co-7 |
| Numerics C7L / C7S (317–320) | printed; "23 of 23" by construction | D3co-1 |
| Numerics C8D (329–332): "no neighbouring dense window repeats but one alarm another" | wrong: none repeats | D3co-3 |
| Numerics C8K, C9A (332–342) | verified | [own] |

## Status of the round-89 issues

| D3cj id | status | what remains |
|---|---|---|
| D3cj-1 | fixed | `c8d_alarms` and the split are printed; C7S has been run and is in (iii), the Status and the Numerics. Remaining: "isolated" is undefined, and the "repeat" statistic does not test repeats (D3co-3); the start dependence is not stated for the least raise (D3co-2) |
| D3cj-2 | fixed | — |
| D3cj-3 | fixed in substance | Draw 11's contrast, draw 5's own offsets and the range over draws still at the feature are in. Remaining: "vacated" (D3co-4) |
| D3cj-4 | fixed in the draft and the code | The selection is right in (iii); superset and persistence are printed. Remaining: the superset holds by construction (D3co-1, which also corrects my predecessor's reading); the md's C7L header (D3co-6) |
| D3cj-5 | fixed | The tests added in round 89 kill R03, R05, R11, R16, R21, R22 and R31 [own]. R01, R15 and R32 survive, as recorded [own]. All 19 of the authors' own mutants are killed [own]. New survivors: D3co-5 |
| D3cj-6 | fixed | — |
| D3cj-7 | fixed | Residual conflation of the counts (D3co-7) |
| D3cj-8 | fixed | — |
| D3cj-9 | fixed | Reference details (D3co-8) |
| D3cj-10 | mostly fixed | Five of its six points are done. Open: the four runs C8R never reproduces are not named (D3co-7) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3co-1 | minor | **The C7L superset is a tautology** (D60:92, 318; md C7L; `c7l_lines`).<br>• **Why it holds:** the 20-period window starts from the same state under the same noise. The integrator steps from t₀ to the same output times, so the first 10 periods are the 10-period window. I checked bit-identity on s = 0.02 from 150, Δ = 0.12: max |Δc| = 0 [own]. In all 23 windows every raise that leaves within 10 periods leaves in the same period over 20 [raw/own].<br>• **What follows:** "contains … in all 23" is a check of determinism. "Grows" and "gaps only filling or splitting" hold by definition. No gap actually splits in the data; all shrink [md].<br>• **What is measured:** that 14 of 23 sets change and that 3 of 13 gaps of the basin close. | • Write "by construction (the 20-period run's first 10 periods are the 10-period run) the set can only grow; over 20 periods it changes in 14 of 23, and 3 of the 13 gaps of the basin close" in (iii), the Numerics and `c7l_lines`.<br>• Turn the count into an assertion in a test. |
| D3co-2 | minor | **The start dependence is of the set's low end, at sub-period scale** (D60:94–98, 189, 319; md C7S).<br>• **Gaps are not the only casualty:** the least ejecting raise also changes between neighbouring starts:<br>&nbsp;&nbsp;– s = 0.03 from 50: 0.1225, from 49: 0.5 (×4.1);<br>&nbsp;&nbsp;– s = 0.03 from 175: 0.2, from 176: 0.48 (×2.4);<br>&nbsp;&nbsp;– k_e = 20, s = 0.04 from 25: 0.225, from 26: 0.37;<br>&nbsp;&nbsp;– m = 4, s = 0.08 from 100: 0.325, from 99: 0.56 [raw].<br>• **The draft's examples are the mild ones.** Its "least ejecting raise moves along the draw" (0.0675–0.325) is sampled every 25 periods, but neighbours one period apart differ as much.<br>• **Neighbours still share most of the set:** median Jaccard 0.975 (min 0.333), against 0.706 for starts 25 or more periods apart at the same configuration, level and draw. The least raise's median |log ratio| is 0.019 against 0.312 [own]. So the bulk of the set (the high raises) is robust, and the low end is not.<br>• **Why one period is coarse:** the sets are decided early. Of C7's 7435 band leaves, 58% fall in the window's period 2 and 71% in periods 1–3 [own].<br>• **The low end varies below one period:** with the noise slid by quarter periods between 49 and 50 (the orbit's phase kept), the least raise runs 0.5, 0.1275, 0.58, 0.21, 0.1225, and the gap at 50 is absent at 49.25, 49.5 and 49.75 [own].<br>• **The neighbourhood itself:** a one-period shift is the finest that keeps the start's orbit phase, so it is a natural neighbourhood, but it samples a set that varies faster. | • Add to (iii): "the start moves the least ejecting raise too (×4 between periods 49 and 50 at s = 0.03), while neighbouring starts share most of the set (median Jaccard 0.975 against 0.71 for starts 25 periods apart)".<br>• Optionally print the least raise and the Jaccard in `c7s_lines`.<br>• Optionally add a sub-period noise-phase scan as C7S's control. |
| D3co-3 | minor | **"Isolated", "repeat" and the evidence for "not the state's"** (D60:120–124, 331–332; `c8d_alarms`, C8D block).<br>• **"Isolated" is a residual class.** It is `kind = 'isolated'` for any alarm whose lead is not 1; isolation itself is never tested. On the data the classes separate cleanly: leads are 1 (9 alarms) or 11–86 (12), with none between 2 and 10 [own].<br>• **The one "leaves by then" is not a repeat.** It is window 50 of k_e = 20, s = 0.08, draw 1, whose loop leaves at period 52 while it covers window 58's alarm period, 60 [own]. No other dense window leaves at any alarm's period, and the windows one period before and after each of the 12 never leave (largest distance 0.03–0.29) [own]. So the Numerics' "no neighbouring dense window repeats but one alarm another" is wrong.<br>• **C8R's "passes all 9" adds nothing new.** It is implied by C8R's aggregate result: every row leaves within one period of P or never, and every alarm period is at least 11 periods before P.<br>• **The direct test supports the conclusion** [own]:<br>&nbsp;&nbsp;– the run's own kernel 1 at the alarm's period is 0.038–0.111 from its feature (draw 7, from 107: 0.038 [C3]);<br>&nbsp;&nbsp;– restarting the 11 alarm windows at k_e = 20 with every formed kernel's offset (C8K's 'w') keeps all 11;<br>&nbsp;&nbsp;– with the centres at their period means too ('wc'), 6 vanish, 1 leaves in period 5 instead of 2, and 4 remain (draw 1, s = 0.06 from 203; draw 3 from 166; draw 1, s = 0.08 from 50; draw 2 from 32).<br>&nbsp;&nbsp;– So the reset of the centres causes about half of them, and kernels 2 and 3's coefficient offsets none. | • Define "isolated" in (iii) and the code ("not the replay one period early").<br>• Test repeats as "another dense window leaving at a ± 1" (0 of 106) and print that.<br>• Fix the Numerics.<br>• Support "not the state's" by the run's own distance at a (0.04–0.11).<br>• Optionally add the 'wc' restart. |
| D3co-4 | minor | **Draw 5: "jumps into the vacated first feature two periods before kernel 1 leaves"** (D60:115). Period means at s = 0.03 [raw C3]:<br>• period 158: kernel 3 at +0.115 (its leave); kernel 1 at −0.533;<br>• period 159: kernel 3 at −0.393; kernel 1 at −0.613, still within 0.3 of −0.52;<br>• period 160: kernel 3 at −0.487; kernel 1 at −0.971 (its leave).<br>The feature is occupied when kernel 3 arrives, and kernel 1 leaves as it arrives: kernel 1 is displaced. This fits C8K: the window from 150 with every offset leaves in its period 10, i.e. at 160, while C8's reduced window from 150 stays (largest distance 0.049) [raw]. "Vacated" is right only for the end positions (D60:280). | Write "kernel 3 leaves at 158 and reaches the first feature while kernel 1 is still at it; kernel 1 leaves at 160, displaced (C3)". |
| D3co-5 | minor | **Mutants: 10 of 33 survive the full suite, none changes the md** (mut/spec94.py) [own].<br>• **Killed:** 23 (Q01–Q04, Q09–Q13, Q15–Q19, Q21, Q22, Q24, Q25, Q27, Q28, Q31–Q33). Most die in `test_round89_md_is_the_report_and_planted_edges`.<br>• **Changes a rerun, unseen (1):** Q30, `main` dispatching `C7S` to `part_C7L`. A rerun with `--only C7S` would store C7L's windows under C7S, and building the md would then fail (KeyError 'base').<br>• **Equivalent by construction (1):** Q08, early when lead ≤ 1 (lead ≥ 1 always).<br>• **Equivalent on the data (8):**<br>&nbsp;&nbsp;– Q05: `rep_alarm` by any (only one repeat);<br>&nbsp;&nbsp;– Q06: a C8R row leaving at a counted as passing (none leaves at a);<br>&nbsp;&nbsp;– Q07: C8R rows started at a (none starts there);<br>&nbsp;&nbsp;– Q14: the conclusion by any (all 9 pass);<br>&nbsp;&nbsp;– Q20: C7S's basin test by the largest gap (all neighbour gaps < 0.1);<br>&nbsp;&nbsp;– Q23: tones in C7S's base (no key collides);<br>&nbsp;&nbsp;– Q26: kernel 3 leaving exactly at the window's start (none does);<br>&nbsp;&nbsp;– Q29: `c7_gap_lines`' guard on multi (no window has more than one band without a gap of the basin).<br>• **The other specs:** all 19 of the authors' mutants are killed. Of the referee's, R03, R05, R11, R16, R21, R22 and R31 are killed; R01, R15 and R32 survive, as recorded. | • Add a test that `main --only C7S` calls `part_C7S` (as for C8D and E40).<br>• Add planted edges for Q05 (two repeats), Q06 and Q07 (a C8R row leaving at, or starting at, a), Q20 (a neighbour with gaps 0.05 and 0.2) and Q26 (kernel 3 leaving at pl).<br>• Record Q08 as equivalent. |
| D3co-6 | cosmetic | **md wording** (`c7l_lines`, the C8D block).<br>• C7L's header "the 23 windows of C7 and C9 whose set had a gap or that are controls": only one of C9's eight multi-band windows was rerun. This is D3cj-4's md half.<br>• "one survives over 20 in 10" should be "a gap of the basin survives …".<br>• "1 leaves by then (each itself one of …)".<br>• `byc` drops the metric from its labels (no WFR alarm on the data). | One clause each. |
| D3co-7 | cosmetic | **C8R in the text** (D60:102–103, 218–219).<br>• Not stated (D3cj-10's last point): the four runs C8R never reproduces at 20 and 40 periods back are C8's three robust misses plus k_e = 20, s = 0.06, draw 7. At 80 back, draw 10 at s = 0.02 replaces draw 7 [own].<br>• New's "33 of 37 … from up to 80 periods before" merges 33/37 (20 and 40 periods back) with 29/33 (80 back). | • Name the four in (iii).<br>• Write "33 of 37 from 20 and 40 periods before, 29 of 33 from 80" in New. |
| D3co-8 | cosmetic | **Sources** (D60:210–212).<br>• Berger–Doan–Siegmund (DCDS-B 9, 463–492, 2008) and Rasmussen (Differential Equations and Dynamical Systems 18, 57–78, 2010) exist as cited [snippet].<br>• "Forward basins: Oljača–Ashwin–Rasmussen" names no work. The paper is "Measure and statistical attractors for nonautonomous dynamical systems", J. Dyn. Diff. Equat. 36 (2024) 2375–2411 [snippet]. The definition of a basin whose fibre is indexed by the initial time I saw only in Oljača's 2022 slides [snippet].<br>• The closest counterpart of a windowed set is a finite-time area of attraction: Giesl–Rasmussen, J. Math. Anal. Appl. 390(1), 27–46 [snippet; year not seen].<br>• The length dependence is monotone by definition (D3co-1). | Give the paper (and the slides); add Giesl–Rasmussen beside Berger–Doan–Siegmund; say "by definition" for the length. |

## Sources

- L. Oljača, P. Ashwin, M. Rasmussen. The 2022 MPE CDT talk slides define the forward basin with fibres indexed by t₀, and note that a change of start time changes which initial conditions converge [snippet]: [slides](https://www.mpecdt.ac.uk/wp-content/uploads/2022/05/12-Early-Career-Speaker-OLJACA-Lea.pdf). The paper, "Measure and statistical attractors for nonautonomous dynamical systems", J. Dyn. Diff. Equat. 36 (2024) 2375–2411 [snippet]: [Imperial record](https://spiral.imperial.ac.uk/entities/publication/857d7eb1-182f-45e9-bccc-9d09bf818308), [Exeter record](https://ore.exeter.ac.uk/repository/handle/10871/130526?show=full), [Rasmussen's page](https://www.ma.imperial.ac.uk/~mrasmuss/).
- A. Berger, T. S. Doan, S. Siegmund, "Nonautonomous finite-time dynamics", DCDS-B 9(3–4), 463–492 (2008) [snippet]: [PDF](https://www.math.ualberta.ca/~aberger/papers/finite_time_dcdsb.pdf).
- M. Rasmussen, "Finite-time attractivity and bifurcation for nonautonomous differential equations", Differential Equations and Dynamical Systems 18, 57–78 (2010), DOI 10.1007/s12591-010-0009-7 [snippet, from a reference list; dx.doi.org is not reachable from this session]: [search result](https://dx.doi.org/10.1007/s12591-010-0009-7).
- P. Giesl, M. Rasmussen, "Areas of attraction for nonautonomous differential equations on finite time intervals", J. Math. Anal. Appl. 390(1), 27–46 [snippet]: [Sussex record](https://sro.sussex.ac.uk/id/eprint/54585/).
- Not rechecked this round: Feudel–Witt–Lai–Grebogi 1998, McDonald–Grebogi–Ott–Yorke 1985, Schultz–Hellmann–Webster–Kurths 2018, Ashwin–Wieczorek–Vitolo–Cox 2012 and Berglund–Gentz 2005 (checked in round 89 [snippet]); Bode, Freudenberg–Looze, Widrow et al., Glover, Yakubovich–Starzhinskii, Wereley, Hoppensteadt, Khalil, Egardt, Ioannou–Kokotović [memory].

## Lesson

- **A count of something that holds by construction is not evidence.** "23 of 23 supersets" checks only that the integrator is deterministic. What the window's length measures is how much the set changes and which gaps close.
- **Test a repeat by the event itself, not by "by then".** The only "neighbour leaving by then" was a different alarm, eight periods earlier.
- **Prefer direct evidence over indirect.** The run's own distance at the alarm (0.04–0.11), and the alarm window restarted from fuller states, say more than other reduced loops passing the same period.
- **Fit a neighbourhood to the scale on which the measured object varies.** The windowed set's low end changes within a quarter period of the noise's phase, so ±1 period samples the start dependence rather than resolving it. Neighbours still share 97.5% of the set: the start dependence sits in the low end of the set, the part a leave depends on.

Key files (W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref94):
- reproduction: W/logs/repro_a.log, W/logs/repro_b.log;
- recounts: W/logs/own_c7.log, W/logs/own_c8d.log, W/logs/own_draw5.log;
- probes: W/logs/own_alarm_probe.log, W/logs/own_c7s_overlap.log, W/logs/own_fracshift.log;
- mutants: W/mut/spec94.py, W/mut/w_spec94/res_s1a.jsonl, W/mut/w_spec94/res_s1b.jsonl, W/mut/w_full/res_s2a.jsonl, W/mut/w_full/res_s2b.jsonl.
