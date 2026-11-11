# Referee report D3cf, round 85: Proposition D‴.60 as revised after round 76

## Scope

- **Commit:** 3773079 ("D‴.60 round 76: under noise the basin is bands, not an edge; the edge tested in time").
  - Exported with `git archive 3773079 | tar -x` into W/repo.
  - W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref85/D3cf_work (subdirectories repo, scripts, logs, out, mut, build).
  - HEAD was 3773079 and the working tree was clean.
- **Files read:**
  - review/drafts/D60.md (275 lines);
  - experiments/exp_noise_resonance.py (2100 lines), and what it uses from exp_closed_transit.py, exp_dead_zone.py (`fitted_state`, `noise_of`), exp_cohesion_leak.py (`load`, `save`) and src/core/closed_loop.py and closed_loop.c;
  - results/chaos/noise_resonance.md (642 lines) and noise_resonance_raw.npz;
  - tests/test_noise_resonance.py (58 tests);
  - review/REVIEW_round76_D3bw.md, the "Round 76" entry of review/REVIEW_LOG.md (lines 3220–3261), and the commit's diff.
- **Settings:**
  - `nice -n 10`, OMP_NUM_THREADS=1 and PYTHONDONTWRITEBYTECODE=1; pytest as `-q -p no:cacheprovider`;
  - the C extension was built into W/build (`LCK_BUILD`);
  - `NR_OUT` pointed at W/out/nr, which held a copy of the raw;
  - mutants ran in private trees under W/mut/trees, which were deleted afterwards.
- **I modified nothing in the repository** and ran no git command that writes; `git status` was clean at the end.
- At most two heavy processes of mine ran at a time. I killed only the explicit PIDs of my own mutant driver and its pytest, after checking /proc/PID/cmdline. Nothing is left running.
- **Marks:** [md] the committed md; [raw] the committed raw; [own] my runs or computations; [snippet] a web-search result; [memory] not verified this round.

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (md5 62f0396eab6b739c9b1f0d50f226fb2c) [own].
  - The 58 tests pass in 97 s [own].
  - 19 cells rerun with the script's own job functions equal the raw in every field [own]:
    - `edge_job` ×3: k_e = 2.5; ω_ref = 0.4; WFR at k_e = 40;
    - `c8_job` ×3: draw 10 at s = 0.02; k_e = 20, s = 0.08, draw 1; WFR at s = 0.32, draw 1;
    - `c9_own_job` ×7;
    - `c7_job` ×4: s = 0; 10 rad/s under −n at 0.06; draw 7 from period 125 at 0.02; k_e = 20 from period 50 at 0.08;
    - C9's scans of draw 3 from 150 and of draw 10 from 110.
  - My dense C8 windows equal the raw's 1563 windows wherever they coincide [own].
- **Every number in (iii), the Status, New and the Numerics is printed by the md** [own].
  - Four are not printed verbatim: 8.1% (= 0.0731/0.9); 0.97–1.03 and 0.89 (roundings of 0.972–1.030 and 0.888); and §10.1.
  - Recounted and right: 18 of 120, 27 of 33, 2 of 1184, 31 of 33, 465 of 1184, 20/0/3/25, 12 of 12, 33 of 77, 33.56, C3's end positions, and the anatomy's t/T and positions [md; own].
- **The mechanism of (iii) is supported, as a reduced model** [own].
  - The fit with kernel 1 raised by the run's δ(p), under the run's own coming noise, leaves when the run leaves.
  - At δ = 0 under the same noise it flags 0 of 33 leave windows.
  - With the run's δ but other noise of the same draw (from p + 75, 150 or 225) it flags 2, 6 and 2 of 33, and 13–21 of 1184 earlier windows.
  - No threshold on δ calibrated on C7 comes near it: at best 21 of 33 with 85 of 1184 earlier windows flagged, or 26 of 33 with 213.
  - So a leave needs both the drifted coefficient and the particular coming noise.
- **C8's aggregate numbers do not depend on the tiling** [own].
  - At all ten phases of the 10-period tiling, 25–29 of the 33 out-of-sample leaves are flagged and 0–5 of 1133–1184 earlier windows.
  - No run that never left raises an alarm at any phase.
  - With a window every 5 periods, 29–30 of 33 leaves fall in a flagged window.
- **Some gaps are real** [own].
  - In 12 of C7's windows the raises inside the gaps keep the kernel within 0.074 of its feature [raw].
  - These bands are unchanged under rtol 10⁻¹¹/atol 10⁻¹⁴ and at 128 outputs a period: not one raise changes.
  - Without noise the ejecting set stays one interval on finer grids too: 0.0732–1 at m = 12, k_e = 10 on a step of 0.0001 near the edge; 0.135–1, 0.105–1 and 0.325–1 at m = 4, k_e = 20 and WFR.
- **Further checks that hold** [own]:
  - the power-law fit α = 0.537 ± 0.025, β = 0.555 ± 0.028 and the ω_ref slopes 0.72, 0.77, 0.71 refit;
  - an independent recount of s*: 33/13/12/19, the null's mean 33.52 over 2000 permutations;
  - C9's 48 pairs equal the C9D records exactly;
  - below the edge the coefficient returns to the fit: +0.0722 → +0.0294 → +0.0062 → +0.00003 at periods 1, 100, 300 and 1000;
  - noise-free detachment does not return over 40 periods;
  - under noise the kernel ends on the side of the passage it detached in, in 21 of 22 runs;
  - "a drift, then a jump, as in bursting but without return" is accurate [md].

**What does not hold (majors first)**

1. **D3cf-1 (major).** "Under noise the basin is no longer an interval … bands with gaps in 18 of the 120 windows" overgeneralises.
   - Six of the 18 windows are artefacts of the 0.3 cut, including the quoted single tone 13.1 rad/s at s = 0.06.
   - The count depends on the cut: 13–19.
   - 102 of the 120 sets are one band.
   - The sets also depend on the window's length.
2. **D3cf-2 (major).** C8's window does not "anticipate" the leave; it re-simulates it with the run's own coming noise.
   - It leaves in the run's own leave period, and its lead equals the window's length.
   - Started 80 periods earlier, the same reduced loop still leaves within one period of the run in 29 of 33 runs.
   - The noise-free edge it is compared with sees no noise at all. The informative nulls (above) are not reported.
3. **Minors:**
   - the misses' causes, and their dependence on the tiling's phase (D3cf-3);
   - C9's windows were chosen around the leaves, and its gap holds over 10 periods only (D3cf-4);
   - at k_e = 2.5 the kernel ends outside the range, not "+0.214 inside" it; OLS standard errors on curved residuals; the onset's "≈" (D3cf-5);
   - the single-tone statements (D3cf-6);
   - 32 of 72 mutants survive the tests, 23 of them not equivalent (D3cf-7);
   - precedents and New (D3cf-8);
   - loose ends (D3cf-9).

## What I checked and how

All paths are under W.

- **Reproduction:**
  - logs/fromraw.log: the md regenerated, then `cmp`;
  - logs/pytest_full.log: the 58 tests;
  - scripts/repro.py → logs/repro.log: 19 cells compared field by field, using the diff in scripts/common.py;
  - scripts/numcheck.py → logs/numcheck.log: every number on D60 lines 52–107, 132–134, 145–151, 164–170 and 213–275 looked up in the md.
- **C8:**
  - scripts/win.py is a parametrised copy of `edge_run`; it is bit-identical to it at the default settings.
  - scripts/c8dense.py dense → out/c8_dense.jsonl: 15,436 windows, one starting at every period, up to the leave + 10 (all 300 for runs that never left). Analysed by scripts/c8dense_analysis.py and c8dense_fix.py (logs/c8dense_analysis.log, logs/c8dense_fix.log):
    - the ten tiling phases;
    - a window every 5 periods;
    - 5-period windows;
    - early alarms, with windows past the 300-period horizon excluded.
  - scripts/c8dense.py nulls → out/c8_nulls.jsonl (6,104 windows): δ = 0, and the run's δ under the draw's noise from p + 75, 150 and 225. Analysed by scripts/nulls_analysis.py → logs/nulls_analysis.log.
  - scripts/horizon.py → out/horizon.jsonl, logs/horizon_analysis.log: the reduced loop run from P − 10, 20, 40 and 80 to P + 2.
  - scripts/reduced.py → logs/reduced_analysis.log: on C3's 22 runs, windows that also carry kernels 2 and 3's coefficients, and the centres.
  - scripts/dropped.py: the size of what the reduction drops.
  - scripts/roc.py: thresholds on δ calibrated on C7.
  - scripts/c8_rawanalysis.py and misses.py: the misses' windows.
- **C7:**
  - scripts/c7_gaps.py, c7_cuts.py, c7_single.py: the largest distance inside each gap; the count of multi-band windows against the cut; sets that end below 1.
  - scripts/c7robust.py → out/c7_robust.jsonl, out/c7_noisefree.jsonl, logs/c7robust_analysis.log: five gap windows rescanned at the default settings, at tight tolerances, at 128 outputs a period and over 20 periods; noise-free fine scans of the four configurations.
- **C9:**
  - scripts/c9check.py: the pairs against C9D; C8's coefficient series against C3's;
  - scripts/c9dense.py → out/c9_dense.jsonl and logs/c9dense_analysis.log: `c9_own_job` at every start, for all 12 runs.
- **EDGE, s\* and the tones:**
  - scripts/edgefit.py: OLS, HC3 and jackknife errors, and a curvature F-test;
  - scripts/sstar.py: an independent recount, 2000 permutations, the configuration-mean ⟨w₂⟩, and s* scaled;
  - scripts/tones.py: C2's single-tone runs by tone and seed;
  - scripts/anatomy.py → logs/anatomy.log: 40-period traces on 1024 outputs a period, 1000-period runs below the edge, and the onset's error ratio;
  - scripts/direction.py → logs/direction.log: C3's 22 runs rerun with 1024 outputs a period over [P − 3, P + 1]; every leave period reproduced.
- **Mutants:**
  - W/mut/mutants.py lists 72 single-point mutants of this round's code;
  - W/mut/driver.py runs, for each: a private tree, the targeted tests, then `--from-raw` and `cmp` (results in W/mut/results.jsonl).
    - The first five survivors were also rerun against all 58 tests.
    - The full rerun was then dropped: by grep, no test outside a mutant's targeted list calls its functions.
  - W/mut/equiv.py → logs/equiv.log reruns representative jobs in each job-level survivor's tree and compares them with the raw.
- **Literature:** web searches on riddled and fractal basins, basins of quasiperiodically forced systems, basin stability and finite-time basin stability, noise-induced tipping, random attractors, and adaptive-control drift and bursting.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| Noise-free basin an interval; Δ_e brackets; one interval to Δ = 1 at four configurations (52–55, 236–238, 252–253) | verified | edge_job reruns bit-identical [own]; fine scans 0.0732–1, 0.135–1, 0.105–1, 0.325–1 [own]; the coefficient returns to the fit below the edge (+0.00003 at period 1000) [own] |
| Onset: repulsion, "ẽ ≈ −Δφ_i/k_e", of the sign of −w_i (55–57, 132–133) | sign verified, magnitude not | ẽ over −Δφ/k_e is 0.050 and 0.054 at the onset, median 0.27 and 0.18 over the passage [own] (D3cf-5) |
| Detachment: t/T 0.594, x_m −0.556, −1.334; kernel 3: 0.401, +0.581, −1.343 (57–60, 239–241) | verified | [md]; non-returning over 40 periods, at −1.404 and −1.412 at 40T [own] |
| k_e = 2.5 "comes back once and detaches for good only at t/T = 2.918, to +0.214 inside the range" (61, 241) | half right | it returns once and detaches for good, but ends at +1.403, outside the range [own] (D3cf-5) |
| "The direction is the passage's"; the noisy end positions −1.42 … −1.44 and +1.39 … +1.41 (60–63) | verified, one exception | end positions [md]; 21 of 22 runs follow the passage; the exception is s = 0.03, draw 5 [own] |
| Power law α, β ± errors, residuals, corners; Δ_e(v₁/k_e)^{1/2} = 0.0756 … 0.0982 (63–67, 245–246) | numbers verified; errors understated | refit [own]; HC3 0.060/0.053; curvature p = 0.006 [own] (D3cf-5) |
| Outside-window ratios; ω_ref^{0.71–0.77}; "an interpolation, not a law" (67–71, 242–247) | verified | [md]; slopes 0.724, 0.769, 0.706 [own] |
| C4: the rise is the even part; 1.02–1.59; 0.64–1.01 "before the drift saturates" (71–74, 232–235) | verified, one draw at 1.01 | [md] (D3cf-9) |
| s*: 33 of 77, null 33.56, widths 0.20–0.50, the lists; "the configuration's, not the draw's" (74–80, 247–249) | verified; the null has little power | recount 33/13/12/19 plus 11 untested; null 33.52; the configuration-mean ⟨w₂⟩ gives 33 [own] |
| "Under noise the basin is no longer an interval … bands with gaps in 18 of the 120 windows"; examples; Status's "not intervals" (80–83, 149, 262–263) | overstated | D3cf-1 |
| The least ejecting raise moves along the draw; k_e = 20, s = 0.08 ejects from the fit (83–85, 259–262) | verified | [md] |
| "A formed kernel leaves when the drift has carried its coefficient into the raises the coming noise ejects" (70–71) | supported as a reduced model | nulls and thresholds [own] |
| C8: 27/33, 2/1184, no false alarm; the null's 31/33, 465/1184, 10 (85–90, 263–269) | verified; robust to the tiling | 25–29 of 33 and 0–5 earlier at all ten phases [own] |
| "the window that ejects the kernel anticipates the leave" (86–87) | not established | D3cf-2 |
| C8's misses (90–92) | listed right; causes not | D3cf-3 |
| C9: draw 3's bands; 20/0/3/25; "0–1 periods before"; 12 of 12 (92–99, 269–275) | verified on the chosen windows; nearly true on all | D3cf-4 |
| Single tones: ±n parts, E per tone and level, the grid's resolution, "E·s² up to 0.02–0.03" (99–104, 253–258) | numbers verified; two readings overstated | D3cf-6 |
| "Single tones never jump"; drifts beyond the edges; "an edge from one start does not decide" (104–107) | verified, with a pooled range | D3cf-6 |
| Proof (iii) (132–134) | adequate as "numerical" | "≈", see the onset row |
| Status (iii) (145–151) | mostly accurate | "not intervals" (D3cf-1) |
| Classical (152–163) | accurate as far as checked; precedents missing | D3cf-8 |
| New (163–170) | overstated in three clauses | D3cf-1, D3cf-2, D3cf-8 |
| Numerics EDGE, C6, C7, C8, C9 (235–275) | printed by the md | numcheck [own] |

## Status of the round-76 issues

| D3bw id | status | what remains |
|---|---|---|
| D3bw-1 | fixed in design | C8 is as asked: time-resolved, on each run's own offset and noise, out of sample, with a null; "decides" is gone. Remaining: its reading as an anticipation (D3cf-2), and the misses, the phase and the dense alarms (D3cf-3) |
| D3bw-2 | fixed in substance | Fine scans without bisection, bands printed, C9 added, C6's conclusion rewritten; real gaps in 12 windows. Remaining: cut artefacts, overgeneralisation, window length (D3cf-1); the text still argues from the leaving (0.3) sets, although the jumping sets are printed |
| D3bw-3 | fixed | Window, errors, outside ratios, ω_ref, "interpolation", the law dropped from New. Remaining: the OLS errors (D3cf-5) |
| D3bw-4 | fixed | ±n, E per tone and level, '>1' kept, C2's levels. Remaining: D3cf-6 |
| D3bw-5 | partly | Onset-only wording, non-returning detachment (holds over 40 periods [own]), both directions (21 of 22 follow the passage [own]). Remaining: the k_e = 2.5 ending, "≈" (D3cf-5) |
| D3bw-6 | fixed | 77, the null, the widths, "the configuration's"; low power (D3cf-9) |
| D3bw-7 | not fixed as claimed | 32 of my 72 mutants survive, 23 of them not equivalent (D3cf-7) |
| D3bw-8 | fixed | — |
| D3bw-9 | partly | Bursting wording fixed; New overstated and precedents missing (D3cf-8) |
| D3bw-10 | fixed | — |

## Issues

| id | severity | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cf-1 | major | **The C7 band claim counts artefacts of the cut and overgeneralises** (D60:3, 80–83, 149, 167; `scan_sets`, `c7_lines`).<br>• **Real or cut artefact.** Of the 18 windows "not an interval" [md], in 12 the gap's raises stay within 0.074 of the feature [raw]. In 6 the largest distance inside the gap is 0.252–0.300, just under the 0.3 cut [raw]:<br>&nbsp;&nbsp;– m = 12, k_e = 10, s = 0.01 from period 25;<br>&nbsp;&nbsp;– k_e = 20: s = 0.06 from 75 and 125; s = 0.08 from 50 and 75;<br>&nbsp;&nbsp;– m = 4, s = 0.04 from 75.<br>• **The quoted single tone is an artefact.** At 13.1 rad/s, s = 0.06 ("0.36–0.48 and 0.96–1"), the raises in the gap move the kernel 0.238–0.289 [raw]. The 21.3 rad/s gap is real (0.009–0.013) [raw].<br>• **The count depends on the cut** [own]: 13 at cuts 0.1 and 0.2, 14 at 0.25, 18 at 0.3, 19 at 0.35, 17 at 0.4 and 0.6.<br>• **Most windows are one band.** 102 of the 120 windows have one band [md]:<br>&nbsp;&nbsp;– 98 run from the least raise to Δ = 1;<br>&nbsp;&nbsp;– 2 are empty;<br>&nbsp;&nbsp;– 2 are a single band ending below 1 (0.5–0.98; 0.62–0.7). Their basins are not intervals, yet they are not counted [raw].<br>What does vary is the least raise: 0.0675–0.325 along the draw at s = 0.02 [md].<br>• **The sets depend on the window's length** [own]. Over 20 periods instead of 10:<br>&nbsp;&nbsp;– draw 7's set from period 175 (s = 0.02) becomes 0.095 and 0.18–0.34;<br>&nbsp;&nbsp;– draw 3's gap from period 150 shrinks from (0.12, 0.1575) to (0.12, 0.1375). | • Count only gaps whose raises stay near the feature (largest distance below 0.1: 12 of 120), and print the count against the cut.<br>• Drop the 13.1 rad/s example, or mark it as the cut's.<br>• State that the sets depend on the window's length.<br>• Reword (iii), the Status, the title's "which is not a threshold" and New's "a set with gaps rather than an edge" to: "a window-dependent edge, whose set has gaps in about one window in ten". |
| D3cf-2 | major | **C8 re-simulates the leave rather than anticipating it; its null shares none of its information** (D60:85–90, 167–170; `c8_job`, `c8_lines`).<br>• **Same timing as the run.** The window carries the run's own coming noise. Its loop leaves in the run's leave period or one before it: lag −1 or 0 in nearly all windows, −4 to +9 in a few [own].<br>• **The lead is the window's length.** The earliest flagged window is the earliest one that contains the leave: it starts at P − 10 in 32 of the 33 runs flagged at all [own].<br>• **Longer reach.** Started 20, 40 or 80 periods before the leave and run to P + 2, the same reduced loop leaves within one period of the run in 33 of 37, 33 of 37 and 29 of 33 runs. It never leaves earlier or later; the others never leave [own]. So "2 of 1184 earlier windows" says that this re-simulation does not leave before the run.<br>• **The null is of another kind.** "Against the noise-free edge's 31 of 33 and 465 of 1184" compares this replay with a threshold that sees no noise.<br>• **The informative nulls** share all of the window's information but one ingredient [own]:<br>&nbsp;&nbsp;– δ = 0 under the coming noise: 0 of 33 leave windows, 1 of 1184 earlier;<br>&nbsp;&nbsp;– the run's δ under other noise of the draw: 2–6 of 33, 13–21 of 1184;<br>&nbsp;&nbsp;– thresholds on δ calibrated on C7: at best 21 of 33 with 85 of 1184.<br>These support (iii)'s mechanism, but the text does not report them. | • Describe C8 as what it is: "the fit with the first coefficient raised by the run's offset, under the run's own coming noise, leaves when the run leaves — within a period, from up to 80 periods before".<br>• Report the three nulls beside the noise-free edge, in C8's lines and in (iii).<br>• Drop "anticipates" from (iii), and "a criterion … against the noise-free edge's" from New. |
| D3cf-3 | minor | **The misses' causes, and their dependence on the tiling's phase** (D60:90–92; `c8_verdict`).<br>• **Robust misses.** Three of the six out-of-sample misses occur at all ten phases [own]: draws 5 and 8 at s = 0.03, and WFR at s = 0.32.<br>• **Phase-dependent misses** [own]:<br>&nbsp;&nbsp;– draw 10 at s = 0.02 is missed at 6 of 10 phases;<br>&nbsp;&nbsp;– k_e = 20, s = 0.08, draw 2 at 4 of 10;<br>&nbsp;&nbsp;– k_e = 20, s = 0.06, draw 1 only at the script's own phase: its leave at 221 falls in the first period of the window from 220, and the windows from 211–219 all reproduce it.<br>So "the two leaves without a jump" share no cause.<br>• **A kernel-3 cause.** Draw 5 at s = 0.03 is recovered once kernels 2 and 3's coefficient offsets join the window's state; the centres add nothing more [own]. In that run kernel 3 leaves first, at period 158, and kernel 1 at 160 [md]. The reduction drops kernel 3's offset: +0.032 … +0.065 at s = 0.02, +0.089 … +0.152 at 0.03 [own].<br>• **Dense windows raise alarms.** With a window at every period there are 12 genuine early alarms in 14,435 windows [own]. Eleven are at k_e = 20, s = 0.06 and 0.08 (e.g. draw 2 at s = 0.08 from periods 6, 18 and 32); one is draw 7 at s = 0.02 from period 107. | • Report the range over phases (25–29 of 33; 0–5 earlier windows) and name the three robust misses.<br>• Give the kernel-3 cause of draw 5's miss.<br>• Note that denser windows raise alarms at k_e = 20 (C8's lines; (iii)). |
| D3cf-4 | minor | **C9's windows were chosen around the leaves, and its gap holds over 10 periods only** (D60:92–99; C9_CASES).<br>• **All window starts.** Over every start for the 12 runs, 8 of the 1434 windows that end before a leave are flagged [own]. In 7 the loop leaves one period before the run. One is genuine: draw 7 un-raised, from period 107, where the loop leaves at 117 and the run at 133. The chosen windows (110, 115) skip it [md].<br>• **Windows that contain the leave** are flagged at 9–10 of 10 starts for 11 runs, and at 4 of 10 for draw 10 un-raised [own].<br>• **Over 20 periods** draw 3's set from period 150 contains the raised runs' offsets +0.1406 and +0.1501 [own; md], and those runs leave at 166 [md]. | • Run `c9_own_job` at every start (one 10-period run each) and print the alarms.<br>• Qualify "they stay" by the window's length (D60:93–94). |
| D3cf-5 | minor | **The anatomy and the scaling** (D60:61, 64–66, 132–133, 241).<br>• **The k_e = 2.5 ending.** +0.214 is the centre at the end of the 3-period trace, 0.08 periods after the detachment [md]. Traced for 40 periods [own]:<br>&nbsp;&nbsp;– the kernel leaves at t/T ≈ 1.59;<br>&nbsp;&nbsp;– it is back near its feature at 2.83–2.92;<br>&nbsp;&nbsp;– it detaches upward for good at 2.918;<br>&nbsp;&nbsp;– it is at +1.383 at 10T and +1.403 at 40T: outside the orbit's range.<br>• **The power law's errors.** OLS errors on systematic residuals [own]: all seven residuals at k_e = 10 are negative, and a curvature term is significant (F(2, 10) = 8.9, p = 0.006). HC3 and jackknife errors are 0.060/0.053 and 0.058/0.051, about 2.3 times larger. Adding one row (k_e = 2.5) moves α to 0.444 ± 0.045 and β to 0.497 ± 0.056 [own, mutant EF3].<br>• **"≈" at the onset.** When φ first reaches 0.5, ẽ is 0.050 (kernel 1) and 0.054 (kernel 3) of −Δφ/k_e [own]. Only the sign holds. | • Correct the k_e = 2.5 sentence: it returns once, then is carried out upward to +1.40.<br>• Give robust errors, or say the residuals are curved.<br>• Write "of the sign of −w_i" for "≈" in the proof. |
| D3cf-6 | minor | **Single tones** (D60:99–107).<br>• **The odd part changes sign.** At 10 rad/s the half-difference per s is +3.000, −3.000 and −3.750 at s = 0.0025, 0.005 and 0.01 [md]. "(3–3.75)s in magnitude" hides that it is not an O(s) part at the smallest level.<br>• **E·s² is flat only at 13.1 rad/s** (800–862.5) [md]. Within the claimed "s ≈ 0.02–0.03" [md]:<br>&nbsp;&nbsp;– at 10 rad/s, E falls from 725 to 417 (s = 0.01 → 0.03);<br>&nbsp;&nbsp;– at 15 rad/s it rises from 800 to 1138 (0.005 → 0.02);<br>&nbsp;&nbsp;– at 21.3 rad/s it rises from 212.5 to 312.5.<br>• **An unlisted case.** 10 rad/s also leaves the grid at s = 0.06 under +n and comes back to 0.43125 at 0.08 [md].<br>• **A pooled range.** "+0.494 … +0.964" [md] pools all six tones [md]. At s = 0.08, the 13.1 and 15 rad/s runs rise +0.637 … +0.732 and +0.567 … +0.747 [raw]. The 21.3 rad/s runs (+0.955 … +0.964 [raw]) lie below their own edge, which is beyond the grid [md]. | • Say that the 10 rad/s half-difference changes sign.<br>• Restrict "E·s² up to 0.02–0.03" to 13.1 rad/s, or give the ranges.<br>• Quote the 13.1 and 15 rad/s drifts in the comparison. |
| D3cf-7 | minor | **Mutants.** The log's "39 … all killed" does not carry over: of my 72 single-point mutants of this round's code, 32 survive every test [own]. 40 are killed.<br>• **Equivalent on the data (9):** ER1, ER3, SC2, SC7, C71, C77, C87, MN2; and MN1, unless C3 and C6 are rerun together with C9D.<br>• **Change the md, caught only by `cmp` (11):**<br>&nbsp;&nbsp;– EF3: k_e = 2.5 joins the fit;<br>&nbsp;&nbsp;– EL1: mean for median;<br>&nbsp;&nbsp;– EL4: the wrong corners;<br>&nbsp;&nbsp;– EL5: kernel 3 against the power law;<br>&nbsp;&nbsp;– SC1: the ω_ref rows enter s*;<br>&nbsp;&nbsp;– SC5: 88 testable for 77;<br>&nbsp;&nbsp;– C75: the gap list cut at 10;<br>&nbsp;&nbsp;– C76: mean for median;<br>&nbsp;&nbsp;– C89: draw 1 counted in sample;<br>&nbsp;&nbsp;– C8A: the null at lo;<br>&nbsp;&nbsp;– C8B: every leave marked a jump.<br>• **Change a rerun, seen by neither the tests nor the md (12):**<br>&nbsp;&nbsp;– the trace's criterion, sample, period and x_m (ER2, ER4–ER7): t/T 0.594 → 0.591; the centre after −1.334 → −1.348; x_m wrong at ω_ref = 0.4;<br>&nbsp;&nbsp;– TRACE_K (ER8);<br>&nbsp;&nbsp;– the anatomy at ω_ref computed without its frequency (EJ2);<br>&nbsp;&nbsp;– EDGE_OUT's kernel-3 row (EO3);<br>&nbsp;&nbsp;– SCAN_C7's step (S71), C7_MORE (S72), C7_CASES without −n (S73), and C8_CASES' seeds (C8C). | Add tests for:<br>• the trace's fields on a fake loop with a known crossing;<br>• `edge_job`'s anatomy at ω_ref ≠ 0.2;<br>• the case lists (EDGE_CASES, C7_CASES, C8_CASES, SCAN_C7) against their stated contents;<br>• `edge_fit` with rows below its window;<br>• `edge_lines`' median, corners and kernel-3 rows;<br>• `star_cases` with ω_ref rows and untested cases;<br>• `c7_lines`' median and gap list;<br>• `c8_lines`' split, the null's midpoint and the jump mark. |
| D3cf-8 | minor | **Sources and New** (D60:152–170).<br>• **Classical, and missing.** Basins with islands and fractal or Wada boundaries in quasiperiodically forced systems: Feudel–Witt–Lai–Grebogi 1998; fractal basin boundaries: McDonald–Grebogi–Ott–Yorke 1985 [snippet]. The noise here is a 16-tone quasiperiodic forcing, so C7's sets are an instance.<br>• **Finite-time basin stability** (Schultz–Hellmann–Webster–Kurths 2018) is the windowed version of C7's set [snippet].<br>• **Noise-induced tipping on a slow drift** (Ashwin–Wieczorek–Vitolo–Cox 2012; Berglund–Gentz 2005) is the generic "a drift, then the noise ejects" [snippet].<br>• **New overstates.** "A set with gaps rather than an edge" and "a criterion" (D3cf-1, D3cf-2); "carrying it out in one passage" fails at k_e = 2.5. | • Cite these in Classical.<br>• Restrict New to "in this loop: the window-dependent ejecting set, measured; and the reduced re-simulation that reproduces the leaves".<br>• Qualify "one passage" (k_e ≥ 5). |
| D3cf-9 | cosmetic | • "The leaves come before the drift saturates": draw 1's share is 1.01 [md].<br>• The direction follows the passage in 21 of 22 runs, not always. The exception, s = 0.03, draw 5, detaches in an upward passage and ends at −1.438 [own].<br>• s*'s null has little power: within a configuration s* differs by at most 31% (median 6%) [own], against level steps of 25–50%.<br>• C7's noise-free reference 0.075 sits 0.0018 above the 10-period edge 0.0732 [own], which shifts E by up to 290 at s = 0.0025 [own]. | One clause each in (iii), C4 and C7. |

## Sources

- U. Feudel, A. Witt, Y.-C. Lai, C. Grebogi, "Basin bifurcation in quasiperiodically forced systems", Phys. Rev. E 58, 3060 (1998) — [snippet]: [Aberdeen record](https://abdn.elsevierpure.com/en/publications/basin-bifurcation-in-quasiperiodically-forced-systems/), [UMD listing](https://ireap.umd.edu/node/3175).
- S. W. McDonald, C. Grebogi, E. Ott, J. A. Yorke, "Fractal basin boundaries", Physica D 17, 125 (1985) — [snippet]: [UMD listing](https://ireap.umd.edu/node/3586).
- A. N. Silchenko, D. G. Luchinsky, P. V. E. McClintock, noise-induced escape through a fractal basin boundary — [snippet]: [Lancaster eprint](https://eprints.lancs.ac.uk/id/eprint/10344), [arXiv nlin/0302014](https://arxiv.org/pdf/nlin/0302014) (search result only; arXiv is unreachable from here).
- P. J. Menck, J. Heitzig, N. Marwan, J. Kurths, "How basin stability complements the linear-stability paradigm", Nature Physics 9, 89–92 (2013) — [snippet]: [PIK record](https://publications.pik-potsdam.de/pubman/item/item_18411).
- P. Schultz, F. Hellmann, K. N. Webster, J. Kurths, "Bounding the first exit from the basin: independence times and finite-time basin stability", Chaos 28, 043102 (2018) — [snippet]: [Aberdeen record](https://abdn.elsevierpure.com/en/publications/bounding-the-first-exit-from-the-basin-independence-times-and-fin).
- P. Ashwin, S. Wieczorek, R. Vitolo, P. Cox, "Tipping points in open systems: bifurcation, noise-induced and rate-dependent examples in the climate system", Phil. Trans. R. Soc. A 370, 1166–1184 (2012) — [snippet]: [arXiv 1103.0169](https://arxiv.org/pdf/1103.0169) (search result), [erratum](https://research.ucc.ie/en/publications/erratum-tipping-points-in-open-systems-bifurcation-noise-included/).
- N. Berglund, B. Gentz, *Noise-Induced Phenomena in Slow-Fast Dynamical Systems*, Springer 2005 — [snippet]: [book page](https://idpoisson.fr/berglund/book.html).
- H. Crauel, F. Flandoli, "Attractors for random dynamical systems", PTRF 1994 — [snippet]: [overview](https://en.wikipedia.org/wiki/Pullback_attractor).
- Adaptive-control drift and bursting — [snippet]: [historical review, arXiv 2108.11336](https://arxiv.org/pdf/2108.11336); [Ioannou–Kokotović 1984 record](https://gnosis.library.ucy.ac.cy/handle/7/56932); [Rohrs et al., MIT LIDS-P-1240](https://dspace.mit.edu/bitstream/handle/1721.1/2820/P-1240-15681714.pdf).
- Not rechecked this round: Anderson 1985 (bursting) and Kokotović–Riedle–Praly 1985 were checked in round 76; Bode 1945, Freudenberg–Looze 1985, Hoppensteadt 1966 and Khalil §10.1 — [memory].

## Lesson

A replay is not a forecast. A window that carries the run's own coming noise, and leaves when the run leaves, measures how faithful a reduced state is. Here it is remarkably faithful: from 80 periods before, the fit with one raised coefficient reproduces 29 of 33 leaves to the period. But it anticipates nothing, and comparing it with a threshold that sees no noise mixes two kinds of rule. The informative nulls keep all of the window's information but one ingredient: the offset without its drift, or the drift under other noise. Likewise, a gap in a set defined by a cut has to be checked against the cut. A raise that moves the kernel 0.297 is not inside the basin in any sense that matters, and a third of the counted gaps were of that kind.
