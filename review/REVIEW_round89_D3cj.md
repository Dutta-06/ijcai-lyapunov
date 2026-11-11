# Referee report D3cj, round 89: Proposition D‴.60 as revised after round 85

## Scope

- **Commit:** 5e42eee ("D‴.60 round 85 applied: log entry, mutation testing, tests"), branch claude/vigilant-clarke-hf939b.
  - Exported with `git -C … archive 5e42eee | tar -x -C W/repo`.
  - W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref89. It holds repo, scripts, logs, out, build and mut.
  - While I worked, another session moved the repository's HEAD to 1ac3d05. The four D‴.60 files are byte-identical between 5e42eee and that HEAD (`git diff --stat` is empty), so everything below holds for both.
- **Files read in full:**
  - review/drafts/D60.md (321 lines);
  - experiments/exp_noise_resonance.py (2646 lines; round 85's new code line by line);
  - results/chaos/noise_resonance.md (686 lines) and the raw;
  - tests/test_noise_resonance.py (71 tests);
  - review/REVIEW_round85_D3cf.md, and the Round 85 entry of review/REVIEW_LOG.md (lines 3438–3483);
  - `git diff 3773079 5e42eee` for these four files;
  - tools/mutants/round85_D60.py (109 mutants).
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 15` with the repository's venv;
  - `NR_OUT=W/out`, holding a copy of the raw; `LCK_BUILD=W/build`; pytest with `--basetemp` inside W.
- **I modified nothing in the repository** and ran no git command that writes.
  - At most two heavy processes of mine ran at a time.
  - The only processes I stopped were two of my own shell wait-loops, which matched their own `pgrep` pattern. I checked them through /proc/PID/cmdline and stopped them as their background tasks.
  - Nothing of mine is left running.
- **Tags:** [md] the committed md; [raw] the committed raw; [own] my runs or computations; [snippet] a web-search result only; [memory] not verified this round.

## Summary

**Design.** This round checks the revision after round 85. The checks:
- regenerated the md from the raw, ran the tests, and reran a sample of the round's new jobs bit for bit;
- recounted every new aggregate with my own code from the raw;
- probed the logic of the replay, the nulls and the bands with new computations:
  - the 21 dense early alarms, classified;
  - the δ-rule recalibrated per level;
  - a rank test of s*;
  - the gap windows rescanned from neighbouring starts;
- ran 35 single-point mutants of round 85's code that are not in the spec, and spot-checked 12 of the spec's;
- searched the sources added this round.

**What checks out**

- **Reproducibility is exact** [own]:
  - `--from-raw` regenerates the md byte for byte, twice (md5 f00dd9a4c16716fe32ce7a1d37c0f27f).
  - The 71 tests pass in 94 s.
  - 29 reruns with the script's own job functions equal the raw in every field: E40 ×2; C7L ×2; C8N ×2; C8R ×2; C8K ×1; C8D ×13 windows; C9A ×6 windows; C8 ×1 run.
  - `c9a_jobs` regenerates C9A's 1554 (seed, p, rise, δ) tuples exactly.
  - C8D's windows at p ≡ 0 mod 10 equal C8's 1526 windows exactly.
- **Every new number of (iii), the Status and the Numerics is printed by the md** [own, numcheck over D60:51–322]. The exceptions are roundings (8.1%, 0.97–1.03, 0.89, 1.22, 2.6·10⁻⁴, 1.085, 0.70, 0.49) and the years of citations.
- **Recounted independently and right** [own]:
  - the cut counts 13, 13, 14, 18, 19, 17, 17;
  - 12 basin gaps and 6 cut gaps, with no window mixing the two;
  - 100 one-band windows (98 reaching 1) and 2 empty;
  - C7L: 14 of 23 changed, 13 multi-band, 10 with a basin gap;
  - C8: 27/6/2/1182 out of sample, and none of the 20 runs that never left raises an alarm;
  - C8N: 0/33 and 1/1184; 4, 3, 3 of 33 with 19, 24, 15 of 1184; 2, 2, 1 false alarms;
  - the δ sweep: 16/33 at f = 1.96 (98 earlier) and 25/33 at f = 1.4 (247 earlier);
  - C8R: 33/37, 33/37, 29/33, with lags only −1 or 0;
  - C8D: per phase 27, 25, 26, 25, 27, 25, 25, 28, 29, 28 leaves flagged and 2, 2, 1, 3, 1, 1, 5, 0, 3, 1 earlier windows flagged; the three robust misses; 29–30 of 33 with a window every 5 periods;
  - C8K: only draw 5 at s = 0.03 changes (miss → hit), and no hit is lost;
  - the power law: α 0.5366 ± 0.0254 (HC3 0.0602, jackknife 0.0582); β 0.5552 ± 0.0277 (HC3 0.0525, jackknife 0.0507); F(2, 10) = 8.93, p = 0.0060;
  - the s* spread (max 30%, median 5%, 20 of 22 below 20%).
- **The window-length claim holds and is stronger than stated** [own].
  - In all 23 C7L windows the 20-period set contains the 10-period set: the set grows with the window.
  - 9 of C7's 12 basin gaps persist over 20 periods and 3 close.
- **"A leave needs both the drifted coefficient and the coming noise" survives a stronger competitor** [own]. A δ-only rule calibrated per configuration and level (leave-one-draw-out) reaches at best:
  - 23 of 33 leaves with 90 earlier windows flagged and 3 false-alarm runs, or
  - 27 of 33 with 192 earlier windows,
  - against the replay's 27 of 33 with 2.
- **The round-85 majors are fixed.** The title, (iii), the Status and New now describe a window-dependent edge with gaps counted against the cut, and C8 as a replay with its nulls.

**What does not hold (no major issue this round)**

1. **D3cj-1 (minor).** The "21 dense windows flagged before the leave" mix two kinds of flag, and the genuine early ones are artefacts of the window's start.
   - 9 are the replay one period early.
   - 12 are isolated single-start alarms, 11 of them at k_e = 20, 11–86 periods early. No neighbouring window covering the same period repeats any of them, and C8R's loop, started earlier, passes 9 of them.
   - The C7 gaps are partly start-sensitive too: 2 of 5 basin gaps vanish when the window starts one period earlier.
2. **D3cj-2 (minor).** "The prediction is the configuration's, not the draw's" is not supported by a null the draft itself calls powerless. The draws' s* order their first jumps in 40 of 58 untied pairs, but with four shared draws neither reading can be established.
3. **D3cj-3 (minor).** Draw 5's miss is attributed to "kernel 3 leaves first", but that also holds for draw 11 at s = 0.03 (30 periods first), and the reduced window hits that one. What sets draw 5 apart is that its kernel 3 jumps into the first feature. The quoted kernel-3 offset range is the range over all eleven draws, and it includes a kernel already ejected.
4. **D3cj-4 (minor).**
   - C7L's selection is misdescribed: one of C9's eight multi-band windows was rerun, not "C9's with more than one band".
   - The monotonicity, and the 9-of-12 persistence of the basin gaps, are not stated, although they are what makes "depends on the window" precise.
5. **D3cj-5 (minor).** Of my 35 mutants, 10 are killed by the tests. The other 25 survive the tests:
   - 4 change the published md and are caught only by `cmp`;
   - 3 change a rerun and are seen by neither the tests nor the md;
   - 18 are equivalent on the data.
6. **D3cj-6 to D3cj-10 (cosmetic).**
   - The md and the Numerics say C9A covers "20 runs"; it covers 12.
   - New's "to the period" should be "within one period".
   - Single-tone loose ends.
   - Sources: the cited tipping literature treats stochastic noise, while here the forcing is deterministic multi-tone; and the nonautonomous-basin literature is missing.
   - Wording.

## What I checked and how

All paths are under W.

- **Reproduction:**
  - logs/fromraw1.log and fromraw2.log, then `cmp` and `md5sum`;
  - logs/pytest_full.log;
  - scripts/repro.py → logs/repro.log: 29 jobs or windows compared field by field (exact, `secs` excluded) with the diff in scripts/common.py;
  - scripts/numcheck_full.py: every number on D60:51–322 looked up in the md.
- **Own recounts from the raw:**
  - scripts/own_c7.py → logs/own_c7.log: bands, cuts, gaps; C7L changes, supersets, persistence;
  - scripts/own_c8.py → logs/own_c8.log: C8, C8N, the δ sweep extended to f = 6, C8R lags, C8D phases and alarms, C8K verdicts, C3's kernel order;
  - scripts/own_misc.py: C8K offsets per draw; C9A's runs; s* spread;
  - scripts/own_fit.py: OLS, HC3, jackknife, curvature F and the cross term.
- **Logic probes:**
  - scripts/own_onset.py, own_onset2.py, own_isolated.py: the leave period inside flagged windows; isolation of each early alarm among the dense windows; C8R passing it;
  - scripts/own_deltarule.py: the per-level δ rule;
  - scripts/own_sstar.py → logs/own_sstar.log: concordance of s* with the first jumps, with within-configuration (5000) and global draw-label (24) permutations;
  - scripts/own_shift.py → logs/own_shift.jsonl: five basin-gap windows of C7/C9 rescanned with `c7_job` from one period earlier and one period later (10 scans).
- **Mutants:**
  - W/mut/spec89.py lists 35 mutants plus a control. W/mut/mutants89.py is tools/mutants.py with REPO set to W/repo, the venv's python, `nice -n 15` and `--basetemp`.
  - Each mutant runs the full test file (`-x`), then `--from-raw` and `cmp`. Results are in mut/res_a.jsonl and res_b.jsonl.
  - The spot check of 12 spec mutants is in mut/res_s1.jsonl and res_s2.jsonl.
  - The mutant trees were deleted afterwards.
- **Literature:** web searches on the five sources added this round, on nonautonomous basins and on finite-time attractivity.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| Title: "an edge that depends on the window" (1–3) | verified | C7L: 14 of 23 sets change; the 20-period set ⊇ the 10-period set in 23 of 23 [own] |
| Noise-free edge, anatomy, onset sign (51–59) | verified (unchanged) | [md]; the sign argument checked [own] |
| k_e = 2.5: returns once, detaches at 2.918, upward, +1.383 / +1.403 (59–61) | verified | E40 rerun bit-identical [own]; "comes back" inferred from a leave in period 2 against the detachment at 2.918 [md] |
| End positions under noise (61–63) | verified | C3 [md]; "with it" is loose for s = 0.03, draw 11, whose kernel 3 left 30 periods earlier (D3cj-10) |
| Power law with OLS, HC3 and jackknife errors, F, "all 9 at k_e = 10 negative", corners (63–68) | verified | refit [own] |
| Outside-window ratios; ω_ref^{0.71–0.77}; "an interpolation, not a law" (68–72) | verified | [md] |
| C4 like for like; "or, at 1.01, as it does" (72–75) | verified | [md] |
| s*: 33 of 77, null 33.56, low power (30%, 5%, 20 of 22) (76–82) | numbers verified; "not the draw's" unsupported | D3cj-2 |
| Bands: 100 of 120 one band (98 to 1), 2 empty, 18 multi-band, cuts 13 … 17, 12 basin / 6 cut, single tones (82–89) | verified | own recount [own]; partly start-sensitive (D3cj-1) |
| C7L: 14 of 23, 13 multi-band, 10 with a basin gap (89–91) | verified; selection misdescribed | D3cj-4 |
| Least ejecting raise along the draw; k_e = 20 from the exact fit (91–92) | verified | [md] |
| C8 a replay: 27 of 33, 2 of 1184, no false alarm, jumps 25 of 33 (92–96) | verified | [own] |
| C8R: 33/37, 33/37, 29/33, never earlier or later (96–97) | verified | lags only −1 (10, 10, 9) and 0 (23, 23, 20) [own] |
| C8N and the δ sweep; "a leave needs both" (97–102) | verified, and robust to a per-level δ rule | [own] |
| The noise-free edge as a rule that sees no noise (102–103) | verified | [md] |
| The six misses (103–105) | verified | [md] |
| C8D phases: 25–29 of 33, 0–5 earlier, step 5 29–30, three robust misses (105–108) | verified | [own] |
| Draw 5's miss "the reduction's", kernel 3 first, C8K (108–110) | recovery verified; the cause under-qualified | D3cj-3 |
| Dense alarms: 21 by configuration (110–112) | counts verified; their nature not stated | D3cj-1 |
| C9: bands decide; 48 pairs 20/0/3/25; 12 of 12; over 20 periods (112–120) | verified | [md]; C7L's draw-3 line [md] |
| C9A: 8 of 1434, 7 within one period, draw 7 from 107; 9–10 of 10 for 11 runs, 4 of 10 (120–122) | verified | [raw]; the 107 alarm is isolated (D3cj-1); "20 runs" in the md and Numerics wrong (D3cj-6) |
| Single tones: odd part's sign change, E per tone, E flat only at 13.1, pooled drifts (123–133) | numbers verified; two loose ends | D3cj-8 |
| Proof (iii) (158–160) | adequate as "numerical" | — |
| Status (iii) (171–179) | accurate | — |
| Classical (180–195) | the new sources exist as cited [snippet] | D3cj-9 |
| New (195–203) | mostly accurate; "to the period" overstated | D3cj-7 |
| Numerics C7L, C8, C8N, C8R, C8D, C8K, C9A (300–322) | printed by the md | "20 runs" (D3cj-6) |

## Status of the round-85 issues

| D3cf id | status | what remains |
|---|---|---|
| D3cf-1 (major) | fixed | Done: counts against the cut; basin against cut gaps; the 13.1 rad/s example marked as the cut's; C7L run; title, (iii), Status and New reworded. Remaining: D3cj-4 (C7L's selection and the monotonicity); D3cj-1 (some small gaps are start-sensitive) |
| D3cf-2 (major) | fixed | Done: "replays … rather than forecasting"; C8R and C8N; the sweep; New restricted. Remaining: New's "to the period" (D3cj-7) |
| D3cf-3 | fixed | Done: the phases; the three robust misses; C8K; the dense alarms by configuration. Remaining: draw 5's cause (D3cj-3); the dense alarms not sorted by kind (D3cj-1) |
| D3cf-4 | fixed | Done: C9A, with each alarm printed with its own lead; the 10/20-period qualifier. Remaining: "1554 windows of 20 runs" (D3cj-6) |
| D3cf-5 | fixed | Done: E40; the robust errors and the curvature test (refit identical [own]); the sign |
| D3cf-6 | mostly fixed | Done: the sign change; E flat only at 13.1; the pooled range named. Remaining: 10 rad/s at s = 0.06–0.08 still unlisted; "E·s² … only to leading order" (D3cj-8) |
| D3cf-7 | fixed in design | Spec of 109 plus tests. My spot check: CL4, ER2, EJ2, WR2, L75, MD6, MK2, RW8, ES2 and KC1 are killed by the tests; GD2 and SP4 survive with the md unchanged, as recorded. New survivors: D3cj-5 |
| D3cf-8 | fixed | Done: five sources added; New restricted; "one passage" limited to k_e ≥ 5. Remaining: D3cj-9 |
| D3cf-9 | fixed | Done: "1.01"; the s* power printed; one blank line; "1 false alarm"; the direction. Remaining: the s* conclusion (D3cj-2) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cj-1 | minor | **The dense windows' early alarms mix two kinds, and the genuine ones are the window's start, not the state** (D60:110–112, 312; `c8_more_lines` C8D block).<br>• **9 of the 21 are the replay one period early.** All are at k_e = 10 or m = 4, and the window's loop leaves in its last period, one period before the run. This is the same lag −1 that C8R shows in 10 of 33 runs [own].<br>• **12 are isolated single-start alarms**, 11 periods or more early: draw 7 at s = 0.02 from period 107 (16 periods), and 11 at k_e = 20 (s = 0.06: draws 1, 2, 3×3; s = 0.08: draws 1×3, 2×3), all leaving in the window's second period, 11–86 periods before the run [own].<br>• **None of them is repeated by a neighbour.** Of the nine other dense windows that also cover the alarm's leave period, none leaves by then (7 for the alarm from period 6), apart from one alarm window covering another [raw/own].<br>• **C8R passes them.** Its loop, started 20–80 periods earlier, passes 9 of the 12 alarm periods without leaving and leaves with the run [raw]. The same holds for C9A's draw 7 from period 107, and for C8's "2 of 1184" (k_e = 20, from periods 210 and 50).<br>• **The C7 sets are partly start-sensitive too.** I rescanned five basin-gap windows from one period earlier and later [own]:<br>&nbsp;&nbsp;– the gap persists at both neighbours for draw 7 at s = 0.02 from 125 and 175 (175's bands rearranged: from 174 they are 0.0725, 0.0925–0.095, 0.195–0.21, 0.38–1), and for draw 3 from 150;<br>&nbsp;&nbsp;– it vanishes at one neighbour for m = 4, s = 0.08 from 100 (from 99 the set is 0.56–1) and for s = 0.01 from 75 (from 74 it is 0.085–1). | • Split the 21: 9 the replay one period early, 12 isolated starts (11 at k_e = 20, leaving in the window's second period).<br>• State that no neighbouring start repeats them and that C8R passes them, so they are artefacts of the reduced start.<br>• Add to (iii) and the Status: "the windowed set depends on the window's start as well as its length; one-raise bands at the low end (s = 0.01 from 75; m = 4, s = 0.08 from 100) vanish when the window starts one period earlier".<br>• Print the split from `c8_more_lines`: the lead, and the neighbours' verdicts. |
| D3cj-2 | minor | **"The prediction is the configuration's, not the draw's" is not supported by the null** (D60:78–79, 282–284; `star_null`).<br>• The draft itself says the null has little power [md], and a powerless null cannot show "not the draw's".<br>• A rank statistic sees structure that the count does not. Within configurations the draws' s* order their first jumps (a no-jump counted as beyond the levels run) in 40 concordant against 18 discordant pairs; within-configuration permutation p = 0.019; restricted to pairs where both jumped, 22 against 10, p = 0.063 [own]. Example: at m = 12, k_e = 5, 10 and 20, draws 2 and 3 have the lower s* and jump at the lower level.<br>• The four draws are shared across configurations, so per-configuration permutations are anti-conservative. Under a global relabelling of the draws the observed concordance is reached by 8 of 24 labellings [own]. | Write "with four draws the test cannot tell the draw's ⟨w₂⟩ from the configuration's (the null has little power)", and drop "not the draw's" (D60:78; md EDGE paragraph). Optionally print the concordance with the global-permutation rank. |
| D3cj-3 | minor | **Draw 5's cause is under-qualified, and the quoted offsets are pooled** (D60:108–110, 313–314).<br>• **Kernel 3 leaving first is not unique to draw 5.** At s = 0.03, draw 11's kernel 3 leaves at 172, 30 periods before kernel 1 (202), yet the reduced window hits draw 11 (flags 200) [md].<br>• **What sets draw 5 apart:** its kernel 3 jumps downward into the vacated first feature (it ends at −0.513 [own]; "−0.50 … −0.54, in the vacated first feature" [md]) two periods before kernel 1 leaves.<br>• **The range is not draw 5's.** "+0.089 … +0.152" is the range over the eleven draws (draw 5's own +0.124 [own]).<br>• **It includes an ejected kernel.** For draw 11 kernel 3 has already left: its centre is at +1.359 in period 199 [own], so +0.144 is not an offset at its feature. All eleven kernel-3 offsets at s = 0.03 exceed kernel 3's own noise-free edge, 0.074 [own; md]. | • Write "there kernel 3 jumps into the first feature two periods before kernel 1 leaves (at 158 and 160)", and note that draw 11's kernel 3 leaves 30 periods first without the reduction missing.<br>• Give draw 5's own offset, and print the range over the draws whose kernel 3 is still at its feature (`c8_more_lines` C8K offsets). |
| D3cj-4 | minor | **C7L's reading** (D60:89–91, 175, 300–301; `c7l_lines`).<br>• **"C7's and C9's with more than one band" is wrong for C9.** C9 has eight multi-band windows (draw 3 from 150, 155, 165; draw 7 from 110, 115, 125, 130; draw 10 from 100) [md], and only one, draw 3 from 150, was rerun. The md's header has the same looseness.<br>• **Not stated: the sets grow with the window.** The 20-period set contains the 10-period set in all 23 windows [own]. As the window lengthens, gaps can only fill or split, and the least ejecting raise can only fall, reaching 0 once the window contains the un-raised run's leave.<br>• **Not stated: of C7's 12 basin gaps, 9 persist over 20 periods** and 3 close (draw 7 at s = 0.02 from 150; k_e = 20, s = 0.04 from 25; m = 4, s = 0.06 from 0) [own]. The "10 of them" includes C9's draw 3. | • "C7's 18 multi-band windows, four controls and C9's draw 3 from period 150".<br>• Add: "the 20-period set contains the 10-period one in all 23; of C7's 12 basin gaps 9 persist and 3 close" (print it from `c7l_lines`). |
| D3cj-5 | minor | **Mutants beyond the spec: 25 of 35 survive the tests** (spec89.py) [own].<br>• **Killed by the tests (10):** R07, R08, R20, R24, R25, R29, R30, R33, R34, R35.<br>• **Change the published md, caught only by `cmp` (4):**<br>&nbsp;&nbsp;– R03: `c7_gap_lines` counts a band ending at the grid's end as "ending below it";<br>&nbsp;&nbsp;– R05: C8R's lag −1 counted as earlier (0 → 10 earlier);<br>&nbsp;&nbsp;– R11: C9A's "within one period" made strict (7 → 0);<br>&nbsp;&nbsp;– R31: C8D's phases over the in-sample runs too.<br>• **Change a rerun, seen by neither the tests nor the md (3):**<br>&nbsp;&nbsp;– R16: `c8r_job` without the 80-period start; the tests use P = 60 and 30, where h = 80 is skipped anyway;<br>&nbsp;&nbsp;– R21, R22: `main` builds C8D and E40 from the stored C8 and EDGE (`old`) rather than this run's (`res`), so a full rerun with a changed C8 or EDGE would use stale inputs; RW7 and RW8 pin only C9A and C8K.<br>• **Equivalent on the data (18):** R01, R02 and R23 (no window mixes a basin gap with a cut gap); R04 and R27 (every C7L window has a short twin); R06; R09 (33 leaves at every phase); R10; R12 (every run has 10 containing starts); R13 (every C9A run left); R14 (no WFR run matches C3); R15 (exactly equivalent); R17; R18; R19; R26; R28 (98 < 100, 247 < 250); R32 (lags only −1 or 0). | Add tests for:<br>• `c7_gap_lines`' one-band-to-the-end count on a synthetic scan;<br>• C8R's "earlier" with a lag of −1;<br>• C9A's "within one period" with a lead of exactly 1;<br>• C8D's out-of-sample filter in the phase loop;<br>• `c8r_job` with P > 81 (all three h);<br>• `main --only C8,C8D,E40` passing this run's C8 and EDGE. |
| D3cj-6 | cosmetic | **C9A's run count** (D60:321; md last paragraph; `c9a_lines`). "1554 windows of 20 runs": `len(runs)` counts the lookup table of C3's and C6's 20 runs at s = 0.02, but C9A ran 12 (draws 3, 7, 10, each with Δ = 0, 0.03, 0.05, 0.065) [own]; (iii) itself says 12. | Count the runs present in C9A (e.g. `len({(x['seed'], x['rise']) for x in C9A})`) and fix the Numerics. |
| D3cj-7 | cosmetic | **New overstates the replay** (D60:201–202). "reproduces the leaves, to the period from up to 80 periods before": the replay reproduces 33 of 37 leaves (29 of 33 from 80 periods back), within one period, and is one period early in 10 of 33 runs (9 of 29 from 80 back) [own]. | "reproduces 33 of 37 leaves within one period, from up to 80 periods before". |
| D3cj-8 | cosmetic | **Single tones, the residue of D3cf-6** (D60:127–130).<br>• 10 rad/s also leaves the grid at s = 0.06 (+n: >1; −n: only 0.22–0.325, a band ending below 1), and its least raises at 0.08 are 0.84 and 0.0225 (odd part +5.109 per s) [md]. Neither is listed.<br>• "E·s² is a tone's edge only to leading order" is not shown except at 13.1 rad/s; at 15 rad/s E rises from s = 0.005 on [md]. | List 10 rad/s with the others; write "E·s² describes a tone's edge only at 13.1 rad/s over s = 0.01 … 0.03". |
| D3cj-9 | cosmetic | **Sources** (D60:192–195). The five added sources exist as cited [snippet]:<br>• Feudel–Witt–Lai–Grebogi: fractal and Wada basin boundaries, and islands;<br>• McDonald–Grebogi–Ott–Yorke;<br>• Schultz et al.: finite-time basin stability is the probability of returning within a given time;<br>• Ashwin et al. 2012;<br>• Berglund–Gentz: the hardcover is 2005, the eBook 2006.<br>Two gaps:<br>• Ashwin et al.'s N-tipping and Berglund–Gentz concern random (stochastic) noise, while here the forcing is a deterministic 16-tone signal. "Noise-induced tipping" is an analogy, not an instance.<br>• The start-dependence (D3cj-1) is the standard non-autonomous fact that a basin depends on the initial time: a forward basin is a family of sets indexed by the initial time (Oljača–Ashwin–Rasmussen); finite-time attractivity (Berger–Doan–Siegmund 2008; Rasmussen 2010) [snippet]. | Mark the tipping sources as the stochastic counterpart; cite the non-autonomous, finite-time framework beside Schultz et al. |
| D3cj-10 | cosmetic | **Wording** (D60):<br>• 79–80: a closing "—;" with no opening dash;<br>• 80: "m = 4–6 at k_e = 10" also takes in m = 4.83 and 5.5, which are not "later";<br>• 62–63: "the third kernel … with it" fails for s = 0.03, draw 11 (kernel 3 left 30 periods earlier);<br>• 94: "(on which C7 was built)": the replay rule fits nothing to C7, so in and out of sample matters only for the δ rules;<br>• 106: "0–5 earlier windows with them";<br>• it is not said that C8R never reproduces exactly the three robust misses, plus k_e = 20, s = 0.06, draw 7 [own]. | One clause each. |

## Sources

- U. Feudel, A. Witt, Y.-C. Lai, C. Grebogi, "Basin bifurcation in quasiperiodically forced systems", Phys. Rev. E 58, 3060–3066 (1998) — abstract: fractal and Wada boundaries, basin "islands" [snippet]: [UMD IREAP listing](https://ireap.umd.edu/node/3175), [Aberdeen record](https://abdn.elsevierpure.com/en/publications/basin-bifurcation-in-quasiperiodically-forced-systems/).
- S. W. McDonald, C. Grebogi, E. Ott, J. A. Yorke, "Fractal basin boundaries", Physica D 17, 125–153 (1985) [snippet]: [UMD IREAP listing](https://ireap.umd.edu/node/3586).
- P. Schultz, F. Hellmann, K. N. Webster, J. Kurths, "Bounding the first exit from the basin: independence times and finite-time basin stability", Chaos 28, 043102 (2018) [snippet]: [Aberdeen record](https://abdn.elsevierpure.com/en/publications/bounding-the-first-exit-from-the-basin-independence-times-and-fin), [PIK record](https://publications.pik-potsdam.de/pubman/item/item_22275).
- P. Ashwin, S. Wieczorek, R. Vitolo, P. Cox, "Tipping points in open systems: bifurcation, noise-induced and rate-dependent examples in the climate system", Phil. Trans. R. Soc. A 370, 1166–1184 (2012) [snippet]: [author PDF (search result)](https://empslocal.ex.ac.uk/people/staff/rv211/mypapers/Ashwinetal2012.pdf), [erratum record](https://research.ucc.ie/en/publications/erratum-tipping-points-in-open-systems-bifurcation-noise-included/).
- N. Berglund, B. Gentz, *Noise-Induced Phenomena in Slow-Fast Dynamical Systems: A Sample-Paths Approach*, Springer 2005/2006 [snippet]: [Springer page](https://link.springer.com/book/9781846280382), [author's book page](https://idpoisson.fr/berglund/book.html).
- L. Oljača, P. Ashwin, M. Rasmussen, slides on forward basins indexed by the initial time [snippet]: [MPE CDT slides](https://www.mpecdt.ac.uk/wp-content/uploads/2022/05/12-Early-Career-Speaker-OLJACA-Lea.pdf).
- A. Berger, T. S. Doan, S. Siegmund, "Nonautonomous finite-time dynamics", DCDS-B 9, 463–492 (2008) [snippet]: [PDF](https://www.math.ualberta.ca/~aberger/papers/finite_time_dcdsb.pdf); M. Rasmussen, "Finite-time attractivity and bifurcation for nonautonomous differential equations", DEDS 18, 57–78 (2010) [snippet]: [record](https://dx.doi.org/10.1007/s12591-010-0009-7).
- P. E. Kloeden, pullback attraction and Lyapunov functions (EJDE conf. proc. 05) [snippet]: [PDF](https://emis.univie.ac.at/journals/EJDE/conf-proc/05/k1/kloeden.pdf).
- Not rechecked this round: Widrow et al. 1975, Glover 1977, Bode 1945, Freudenberg–Looze 1985, Yakubovich–Starzhinskii, Wereley, Hoppensteadt 1966, Khalil §10.1, Egardt 1979, Ioannou–Kokotović 1984 [memory]; Riedle–Cyr–Kokotović, Fu–Sastry, Kokotović–Riedle–Praly, Anderson 1985 and Neishtadt 2014, checked in earlier rounds.

## Lesson

A window's verdict belongs to its start as well as its length. An alarm counts as the state's only if three things hold:
- a start one period earlier or later repeats it;
- the same loop run from further back passes through it;
- the leave lands where the replay would put it.

Here every one of the 12 genuine early alarms fails the first test, and the C8R loop passes 9 of them. The 9 other alarms are just the replay being one period early.

The same discipline applies to nulls. A permutation null that its authors call powerless cannot prove that an ingredient is absent. And when four draws recur across all configurations, permuting within each configuration overstates the evidence in both directions: only a relabelling of the draws themselves is honest.

Key files: W/logs/repro.log, W/logs/own_c8.log, W/logs/own_isolated.log, W/logs/own_shift.jsonl, W/logs/own_sstar.log, W/mut/spec89.py, W/mut/res_a.jsonl, W/mut/res_b.jsonl, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref89.
