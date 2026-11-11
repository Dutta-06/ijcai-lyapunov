# Referee report D3cz, round 105: Proposition D‴.60 as revised after round 102

## Scope

- **Commit.** I reviewed e68fcd3 on claude/vigilant-clarke-hf939b, exported with `git -C … archive e68fcd3 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref105. HEAD is now 5ee87f7, a D‴.59 commit. `git diff e68fcd3 5ee87f7` is empty for the draft, the script, the md, the raw, the tests and the mutant spec. Only REVIEW_LOG.md differs: it gains the round-102 and round-103 entries.
- **Read in full:**
  - review/drafts/D60.md (371 lines);
  - the C7S, C8R, C8D and C8A paragraphs of results/chaos/noise_resonance.md;
  - review/REVIEW_round102_D3cw.md;
  - the Round 102 entry of REVIEW_LOG.md, read at 5ee87f7 because it is absent at e68fcd3 (see D3cz-6);
  - tools/mutants/round102_D60_own.py and tools/mutants.py;
  - `git diff 38fcf58 e68fcd3` for all of these files.
- **experiments/exp_noise_resonance.py:** every changed line — `state_start`, `C8A_PROBE_QUIET`, `C8A_VARIANTS`, `c8a_job`, `c8a_lines`, `c7s_overlap`'s `ancr0`, the C8R clause in `c8_more_lines`, the docstring and the usage line — and what they call: `offsets_start`, `raised`, `c_run`, `edge_run`, `win_run`, `c_noise`, `c8a_cases`, `c8d_alarms`, `DZ.fitted_state` and `closed_loop.simulate`'s outputs.
- **tests/test_noise_resonance.py:** the round-102 tests line by line (`test_round98_c8a_job_wiring` as rewritten, `test_round102_state_start_the_rerun_and_computed_conclusions`, the C8R asserts in `test_round89_…`), and `test_round89_md_is_the_report_and_planted_edges`.
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv;
  - NR_OUT=W/out holding a copy of the raw, LCK_BUILD=W/build, TMPDIR=W/tmp;
  - pytest `--basetemp` inside W (for tools/mutants.py through `PYTEST_ADDOPTS`, one basetemp per slot);
  - at most two heavy processes of mine at a time, each under `timeout`. The machine's load was 6–7 from other sessions.
- **Nothing in the repository was modified.** I ran only read-only git commands. The only thing I added inside the export is a symlink, W/repo/venv, to the venv. I killed no process; all of mine ended on their own.
- **Tags:** [md] the committed md; [raw] the committed raw; [own] my runs or code; [snippet] a web-search result only; [memory] not checked this round.

## Summary

**Design.** The revision does four things:
- restarts C8A's 12 isolated-alarm windows one ingredient at a time (`state_start`, six variants);
- adds noise-free rows to draw 5's probe;
- computes `c8a_lines`' conclusions from the counts;
- makes small fixes (the rerun always with the metric, `ancr0`, the C8R clause, wording).

My checks:
1. **Reproduction:**
   - the md from the raw;
   - the tests;
   - all 7 C8A jobs with the script's own `c8a_cases` and `c8a_job`;
   - an independent restart routine of my own that reproduces every C8A row.
2. **Recounts:** every new number, from the raw with numpy only.
3. **Planted inputs** for every printed conclusion.
4. **The decomposition probed further:**
   - the 12 alarm windows from 23 starts each:
     - mean against instantaneous, kernel 1's instantaneous deviation scaled (λ = −1, 0.5, 1.5);
     - the instantaneous state from (p ± 1)T, i.e. the same orbit phase at another noise instant;
     - δ-perturbations;
     - the full state;
     - warm starts one period earlier;
   - from the reduced start, kernel 1's centre alone moved to its mean, to its instantaneous value at pT and (p − 1)T, and to fixed offsets −0.06 … +0.04;
   - other cuts (0.25, 0.35);
   - a specificity control: 36 windows that contain a real leave, from 12 runs.
5. **Draw 5:** all 8 subsets; without the noise, with each kernel's trajectory; under other noise (shifted segments, −n, other draws, other levels); kernel 3's offset scaled.
6. **Mutants:** the authors' 16 and 16 of my own, survivors on the full suite, and one survivor regenerated.
7. **Sources:** one search.

**What holds**

- **Reproducibility is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 d6bb5e5c85ffa2a85501265d6e708272);
  - 79 tests pass (121 s);
  - all 7 C8A jobs equal raw['C8A'] in every field except `secs` — 5 on the rerun path, 2 with C3's period means;
  - my own restart code gives the reduced window and all 8 variants at the 12 windows to 10⁻¹² (leave period and largest distance).
- **Every new number agrees with the raw** [raw, own]:
  - 'w+c1' = 'wc' and 'w+c23' = 'w' at every window, in class and in leave period (12 of 12);
  - kernel 1's period-mean offsets −0.0713 … −0.0325 (removed) and −0.0467 … +0.0028 (not removed);
  - 'i' removes 12;
  - the 6 that 'wc' leaves are removed by exactly one of the three variants each — 4 by 'wc+i1c', 1 by 'wc+i1w', 1 by 'wc+i23' — a partition;
  - noise-free draw-5 rows (−1, 1.178·10⁻⁶) and (2, 0.864);
  - 84 anchored pairs, 84 ratios, median ×1.3609, `ancr0` = 0, no empty set;
  - 7 of 7 reruns leave at C8's period.
- **The fuller starts are specific, not blanket suppressors** [own]. In 36 windows that contain a real leave (P − 8, P − 5, P − 2 of 12 runs, 7 at k_e = 20): 'wc', 'i' and the full state reproduce the run's leave within one period in 36 of 36; the reduced window and 'w' do so in 35.
- **Draw 5's "needs kernel 3's offset" holds over all 8 subsets** [own]: kernel 1 leaves at 160 if and only if kernel 3's offset is applied. Other noise leaves every kernel in place (D3cz-2), so the leave does need the coming noise.
- **The authors' 16 mutants are all killed by the subset of tests** [own], as the log says (U14 also changes the md).

**What does not hold (no major issue)**

1. **D3cz-1 (minor). The decomposition splits one coordinate in two, and its 6/4/1/1 attribution depends on the base and the cut** [own].
   - From the reduced start, the outcome at each of the 11 k_e = 20 alarm windows is monotone in kernel 1's centre offset at the restart, with a per-window threshold 0.0065–0.06 outward of its feature:
     - moved outward by 0.02, 0.03, 0.04 or 0.06 it removes 5, 8, 9 or 11; moved inward by 0.04 it removes none;
     - at its period mean it removes 7, at its instantaneous value at pT 10, at its value at (p − 1)T 11.
   - The period mean and the instantaneous centre are two samples of that one coordinate. At orbit phase 0 the run's kernel 1 sits 0.006–0.073 outward of its period mean at (p − 1)T in 12 of 12 windows, 10 of 12 at pT. Kernel 1's centre taken from (p − 1)T or (p + 1)T removes the same 10 as at pT, so what matters is the orbit phase, not the noise's instant.
   - Draw 2 from 58, attributed to "kernels 2 and 3's state", is removed by kernel 1's mean centre from the reduced start; only kernels 2 and 3's coefficient offsets turn it into a delay. The order of adding matters there.
   - "Its coefficient 1" is draw 7's crossing of 0.310 against the cut 0.3:
     - no centre displacement removes it;
     - a 5·10⁻⁴ lower kernel-1 coefficient removes it, and +0.005 on 'i' brings it back;
     - the instantaneous state from (p ± 1)T keeps it.
   - At the cut 0.25, 'i' and the full state no longer remove draw 1 at s = 0.08 from 50 (largest distance 0.287), and none of the three variants does. At 0.35 draw 7's window is not an alarm.
   - Warm-started one period earlier, the reduced state passes all 12 alarm periods: the alarms belong to the cold restart at pT.
2. **D3cz-2 (minor). Draw 5: what needs the coming noise is kernel 3's direction, not the collision's effect.** Without the noise, kernel 3's +0.124 ejects kernel 3 upward to +1.348, and it never meets kernel 1. At 0.7 of that offset (+0.087), kernel 3 is carried down to −1.387 and kernel 1 is ejected in period 2 (0.871), with no noise. The md's "— the collision needs the coming noise" misnames this. The single noise-free row cannot tell "no collision" from "a harmless collision".
3. **D3cz-3 (minor). The rerun path of `c8a_job` is untested.** Z01 (the rerun's W and C swapped) survives the full suite. Regenerating two of its five jobs changes the md's C8A line: 'w' "kept 8, removed 4", "earlier 4", "the 8 that 'wc' leaves", and the 'w' clause dropped [own].
4. **D3cz-4 (minor). Some printed conclusions are checked on their consequence or are vacuous** [own, planted inputs]:
   - the C8R clause states "none of its loops leaves more than one period before its run" while C8R's own line on the same input says "1 earlier";
   - "the centres' effect is kernel 1's centre alone" prints when the centres have no effect at all;
   - outcomes are compared by class only;
   - "needs kernel 3's offset" is decided on whatever subsets were run.

   Z07, Z08, Z09, Z13 and Z14 survive the full suite, all equivalent on the data.
5. **D3cz-5 (cosmetic).** `c8a_job`'s docstring is still round 98's. `kk = 8` duplicates `c_run`'s k.
6. **D3cz-6 (cosmetic).**
   - The two juxtaposed ranges overlap: −0.0467 is draw 2 from 58, delayed, against −0.0325 and −0.0365 removed.
   - Draw 7's marginality is not stated.
   - e68fcd3's message says "log entry", but the Round 102 entry first appears in 5ee87f7.

## What I checked and how

All paths are under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref105.

- **Reproduction:**
  - W/logs/fromraw.log, then `cmp` against the committed md: identical.
  - W/logs/pytest_full.log: 79 passed, 1 old RuntimeWarning (B3).
  - W/scripts/repro_c8a.py → W/logs/repro_c8a.log: `c8a_job` on all 7 of `c8a_cases`' jobs: identical in every field except `secs`.
- **Recounts (numpy only, `python -I`):**
  - W/scripts/recount.py → W/logs/recount.log: classes, leave periods, `c1` ranges, the partition, margins, the probe rows, reruns against C8.
  - W/scripts/c7s_recount.py: the 84 anchored pairs, J median 0.7056, ratio median ×1.3609, 0 zero, 0 empty.
  - W/scripts/c8d_view.py: C8D's windows containing each leave.
- **Probes** (own code in W/scripts/common.py: `edge_run`'s loop, returning per-period distances):
  - W/scripts/probe_alarms.py → W/logs/probe_alarms.json, analysed by W/scripts/analyse_alarms.py → W/logs/analyse_alarms.log. The 12 windows from 23 starts each:
    - red, w, wc and the 6 variants;
    - wc with kernel 1's centre at mean + λ(inst − mean) for λ = −1, 0.5, 1.5;
    - wc with kernel 1's centre from the instantaneous state at (p ± 1)T;
    - 'i' from (p ± 1)T;
    - 'i' with kernel 1's coefficient ± 0.005 and its centre ± 0.01;
    - the full state at pT (all α, c and x(pT));
    - red and wc warm-started at (p − 1)T;
    - the classes re-read at the cuts 0.25, 0.3 and 0.35.
  - W/scripts/probe_thresh.py → W/logs/probe_thresh.json, analysed by W/scripts/analyse_thresh.py → W/logs/analyse_thresh.log: the reduced start with kernel 1's centre at its mean, at its value at pT or (p − 1)T, or at −0.02, −0.03, −0.04, −0.06 or +0.04 from its feature; per-window monotonicity and threshold.
  - W/scripts/probe_hits.py → W/logs/probe_hits.json, analysed by W/scripts/analyse_hits.py → W/logs/analyse_hits.log: 36 windows containing real leaves (red, w, wc, w+c1, i, full); the reduced window equals C8D's in 36 of 36.
  - W/scripts/probe_draw5.py → W/logs/probe_draw5.log and W/logs/probe_draw5.json: draw 5 from 150 — the 8 subsets under the coming noise; 5 subsets without noise; kernel 3's offset under the noise shifted by +75, +150, −75 or −140 periods, under −n, under draws 1, 2, 3 and 7, and at s = 0.01, 0.02 and 0.04; kernel 3's offset × 0.5, 0.7 and 1.3, with and without the noise. All three kernels' period means and traces recorded.
- **Planted inputs:** W/scripts/synthetic.py → W/logs/synthetic.log (five cases, cited in D3cz-4).
- **Mutants:**
  - Specs: tools/mutants/round102_D60_own.py → W/mut/w_own/res_a.jsonl; my W/mut/spec105.py → W/mut/w105/res_b.jsonl.
  - Stage 1 used `-k "round102 or round98 or round89 or md_is_the_report or c8_more_lines or c8_windows"`.
  - Survivors on the full suite: W/mut/w105full/res_c.jsonl and W/mut/w105full/res_d.jsonl.
  - Z01 regenerated (jobs 3 and 6): W/scripts/regen_mut.py with the mutant copy W/mutdemo/Z01 → W/mutdemo/Z01/c8a_line.txt.
- **Literature:** one web search on Giesl–Rasmussen. arxiv.org is outside this session's allowed network; it has to be allowed in the session's network settings before that full text can be read.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii), their proofs, EDGE, power law, C4, s*, C7, C7L | unchanged since round 102, not reverified | [md] |
| C7S: anchored 0.706 over 84 pairs, ×1.36 "over 84", `ancr0` | verified | [raw], [own] |
| Draw 5: kernel 3's offset alone → kernel 1 at 160; never with kernel 1's or kernels 1 and 2's (120–123) | verified over all 8 subsets | [raw], [own] |
| "Though not to its needing the coming noise: without it … 1.2·10⁻⁶ … kernel 1's own ejects it in period 2" (124–125) | numbers verified. The leave needs the coming noise (other noise: no leave). The reason (kernel 3 goes up) is missing, and the md's "the collision needs the coming noise" misnames it | [raw], [own]; D3cz-2 |
| C8R "passes all 9 of the 12 it covers (implied by its aggregate)" (133–134) | true on the data. The md's clause states the aggregate without checking it | [md], [own]; D3cz-4 |
| 'w' keeps 12; 'wc' removes 6, delays 1, keeps 5; "no alarm from dropping kernels 2 and 3's coefficient offsets" (135–138) | verified | [raw], [own] |
| "The centres' effect is kernel 1's centre alone … at every window" (138–140) | true on the 'w' base, in class and leave period. On the reduced base kernel 1's mean centre removes 7, not 6 | [raw], [own]; D3cz-1 |
| "−0.0713 … −0.0325 where 'wc' removes … −0.0467 … +0.0028 where not" (140) | numbers verified. The ranges overlap; the kept alone are −0.0250 … +0.0028 | [raw]; D3cz-6 |
| "The 6 … go with the formed kernels' instantaneous state, which removes all 12 — kernel 1's instantaneous centre 4, its coefficient 1, kernels 2 and 3's state 1" (140–142) | true at the cut 0.3 and on the 'wc' base. One coordinate (kernel 1's centre at phase 0); base- and cut-dependent; the "coefficient 1" a marginal crossing | [raw], [own]; D3cz-1 |
| Status (iii) (210–211), Numerics C8A (359–363) | accurate | [md] |
| Sources: "a threshold analogue … (Giesl–Rasmussen's are defined by finite-time attractivity)" (233–234) | consistent with the snippets (finite-time (exponential) attractivity, a time-varying Riemannian metric) | [snippet] |

## Status of the round-102 issues

| D3cw id | status | what remains |
|---|---|---|
| D3cw-1 | fixed: the variants are in, the numbers are right, and the partition holds | the interpretation (D3cz-1) |
| D3cw-2 | fixed in part: a real rerun of draw 7's job (the C3 path), and the wiring asserts (m, k_e, s, draw, metric, output 8p, the quiet rows). The authors' 16 mutants are killed | the rerun path is untested: Z01 survives and changes the published C8A on regeneration (D3cz-3) |
| D3cw-3 | fixed: the clauses are computed, and contrary inputs are tested | some conditions are vacuous or weak (D3cz-4) |
| D3cw-4 | wording fixed | the condition checks the passes, not the aggregate it states (D3cz-4) |
| D3cw-5 | fixed: `metric=mt`; `ancr0` and "over 84" | none |
| D3cw-6 | fixed: parentheses, dash pair, "kernels 1 and 2's", Status (iii), docstring, usage line | `c8a_job`'s docstring (D3cz-5) |
| D3cw-7 | fixed in part: two noise-free rows | the mechanism and the noise-free collision (D3cz-2). The 8 subsets were not added (true on the data [own]) |
| D3cw-8 | fixed | none |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cz-1 | minor | **One coordinate, split by the restart's choice of base** (D60:138–142; md C8A).<br>From the reduced start, at each of the 11 k_e = 20 windows, the outcome is monotone in kernel 1's centre offset at the restart [own, logs/analyse_thresh.log]:<br>• thresholds per window: (−0.03, −0.02], (−0.0467, −0.04], (−0.02, −0.0065], (−0.02, 0] ×4, (−0.06, −0.04], (−0.04, −0.03], (−0.03, −0.02], (−0.03, −0.0236];<br>• fixed offsets −0.02, −0.03, −0.04, −0.06 remove 5, 8, 9, 11; +0.04 removes 0;<br>• the period mean removes 7, the value at pT 10, the value at (p − 1)T 11.<br>Phase, not instant: kernel 1's phase-0 value lies outward of its period mean by 0.006–0.073 at (p − 1)T (12 of 12). On 'wc', kernel 1's centre from (p ± 1)T removes the same 10 as from pT. λ = 0.5 or 1.5 removes the same 10; mirrored (λ = −1) keeps 7, i.e. 6 kept and 1 delayed [logs/analyse_alarms.log].<br>Base dependence: draw 2 from 58 ("kernels 2 and 3's state") is removed by kernel 1's mean centre without kernels 2 and 3's coefficient offsets, and only delayed with them.<br>Draw 7 ("its coefficient"): largest distance 0.310 against 0.3. No centre move removes it; 'i' with +0.005 on kernel 1's coefficient keeps it; 'i' from (p ± 1)T keeps it.<br>Cut and δ: at 0.25, 'i' and the full state do not remove draw 1 from 50 at s = 0.08 (0.287), and none of the three variants does; at 0.35 there are 11 alarms and the split is 4/0/1. 'i' with kernel 1's coefficient +0.005 removes 10 of 12.<br>Warm start: red and 'wc' started at (p − 1)T pass all 12. | In D60 (iii) and `c8a_lines`, replace the 4/1/1 reading by:<br>"the alarms are the cold restart's — started one period earlier none of the 12 leaves; at k_e = 20 all 11 go with kernel 1's centre at the window's start alone: from the reduced start each is removed once the centre lies 0.0065–0.06 outward of its feature (inward never); the period mean is outward enough in 7, the run's value at the window's orbit phase (0.006–0.073 outward of the period mean, the same one period earlier or later) in 10–11; draw 7's at k_e = 10 is a crossing of 0.310 that no centre move removes and a 5·10⁻⁴ lower coefficient does".<br>Add 'red+c1' (mean, at pT, at (p − 1)T) and the warm start to `C8A_VARIANTS` (about 1 s each). |
| D3cz-2 | minor | **Draw 5's noise dependence is kernel 3's direction** (D60:124–125; md "— the collision needs the coming noise"; `c8a_lines`) [own, logs/probe_draw5.log].<br>• Without noise, kernel 3's +0.124 leaves in period 1 upward (trace 0.610 … 1.349, end +1.348): no collision.<br>• Without noise at ×0.7 (+0.087), kernel 3 goes down to −1.387 and kernel 1 is ejected in period 2 (0.871).<br>• Under the coming noise, kernel 3 leaves at 158 downward (to −0.500, the first feature) and kernel 1 at 160.<br>• Under the draw's noise 75 or 150 periods later or 75 or 140 earlier, under −n, under draws 1, 2, 3 or 7, and at s = 0.04, nothing leaves. At s = 0.02 kernels 1 and 3 leave at 158; at s = 0.01 kernel 3 goes up.<br>• The quiet rows record only kernel 1, so the clause cannot tell "no collision" from "harmless collision". | D60: "— though not to its needing the coming noise, which sends kernel 3 down: without it kernel 3 is ejected upward (to +1.35) and never meets kernel 1, at 0.7 of its offset it is carried down and ejects kernel 1 in period 2 without any noise, and under other noise (the draw 75–150 periods on, −n, other draws) no kernel leaves".<br>`c8a_job`: record every kernel's leave and end in the quiet rows; add a scaled-offset row and a shifted-noise row.<br>`c8a_lines`: print "kernel 3 … upward/downward" and drop "the collision needs" unless kernel 3 reaches feature 1 in the quiet row. |
| D3cz-3 | minor | **The rerun path is untested** (`c8a_job` line 1584; tests).<br>The real rerun covers only draw 7's job, which takes C3's period means. Z01 (`W, C = Cp, Wp` on the rerun path) survives the full suite (W/mut/w105full/res_c.jsonl).<br>Regenerating jobs 3 and 6 under it changes the md's C8A line [W/mutdemo/Z01/c8a_line.txt]:<br>• "'w' … kept 8, removed 4";<br>• "earlier 4";<br>• "the 8 that 'wc' leaves";<br>• the "no alarm comes from …" clause disappears;<br>• "… where not -0.0250 … +1.4975". | Rerun one rerun-path job against the raw (job 3, 25 s under load), or assert in the mock test that the W and C passed to `offsets_start` and `state_start` are `Wp.tolist()` and `Cp.tolist()` in that order. |
| D3cz-4 | minor | **Weak or vacuous conditions** (`c8_more_lines` line 2645; `c8a_lines`) [logs/synthetic.log].<br>• A C8R loop leaving at 30 for a run leaving at 40, with alarms at 7: the clause "(implied by C8R's aggregate: none of its loops leaves more than one period before its run)" prints, while C8R's line says "1 earlier".<br>• With 'wc' = 'w' everywhere, "the centres' effect is kernel 1's centre alone" prints.<br>• 'w+c1' delaying to 9 where 'wc' delays to 5 still prints "gives 'wc''s outcome".<br>• "Needs kernel 3's offset" prints on the subsets {3} and {1} alone.<br>• Survivors on the full suite, equivalent on the data: Z07 (`all` → `any` in the C8R clause), Z08 (the 'i' check on the first window only), Z09 (an empty set counted as a zero raise), Z13 (", not its own" unchecked), Z14 (a C8R loop leaving at the alarm's period counted as a pass). Z05 and Z06 are equivalent by construction. | • Print the C8R clause only if C8R's own count "earlier" is 0 at every h, and say "implied by" only then.<br>• Require some window where 'wc' ≠ 'w' for the centre clause, and compare leave periods, not classes.<br>• Decide "needs" only when every subset without kernel 3 was run (C8A_PROBE with all 8).<br>• Planted tests for Z07, Z08, Z09, Z13, Z14. |
| D3cz-5 | cosmetic | **Stale code.** `c8a_job`'s docstring (lines 1573–1577) still says "else C's loop rerun" and lists only 'w' and 'wc'. `kk = 8` (line 1587) duplicates `c_run`'s `k = 8`; the test pins kk = 8 but not that the two agree. | Update the docstring (variants, the instantaneous state, the quiet rows); have `c_run` return k, or a module constant. |
| D3cz-6 | cosmetic | **Wording and log.**<br>• D60:140 juxtaposes overlapping ranges; the kept alone are −0.0250 … +0.0028, and −0.0467 is the delayed draw 2 from 58.<br>• Draw 7's marginality (0.310) is not stated.<br>• e68fcd3's message says "log entry", but REVIEW_LOG.md at e68fcd3 has no Round 102 entry; it arrived in 5ee87f7, a D‴.59 commit. | "−0.0713 … −0.0325 where removed, −0.0250 … +0.0028 where kept (the delayed one −0.0467)"; name draw 7's crossing; record in the log which commit carries the entry. |

## Sources

- P. Giesl, M. Rasmussen, "Areas of attraction for nonautonomous differential equations on finite time intervals" [snippet]. The records say the paper introduces new concepts of (exponential) attractivity on a finite time interval, and that its existence condition uses a time-varying Riemannian metric. This agrees with D60's "defined by finite-time attractivity".
  - [Sussex record](https://sro.sussex.ac.uk/id/eprint/54585/)
  - [Sussex figshare](https://sussex.figshare.com/articles/journal_contribution/Areas_of_attraction_for_nonautonomous_differential_equations_on_finite_time_intervals/23417591)
- M. Rasmussen's publications page [snippet]: lists "Finite-time attractivity and bifurcation for nonautonomous differential equations", Differ. Equ. Dyn. Syst. 18 (2010) 57–78. [page](https://www.ma.imperial.ac.uk/~mrasmuss/main.html)
- Bielefeld SFB 701 preprint on finite-time attraction and repulsion areas [snippet; listing only]. [listing](https://www.math.uni-bielefeld.de/~beyn/AG_Numerik/html/de/preprints/sfb_17_011.html)
- Not rechecked this round [memory or earlier rounds' snippets]: Oljača–Ashwin–Rasmussen (round 102 confirmed JDDE 36 (2024) 2375–2411 and the correction), Berger–Doan–Siegmund, Feudel–Witt–Lai–Grebogi, McDonald–Grebogi–Ott–Yorke, Schultz–Hellmann–Webster–Kurths, Ashwin–Wieczorek–Vitolo–Cox, Berglund–Gentz, Bode, Freudenberg–Looze, Hoppensteadt, Khalil.

## Lesson

- **Before splitting an effect into ingredients, scan the one continuous variable they may share.** Here "the centres at their period means" and "the instantaneous centre" are two samples of one coordinate, kernel 1's centre at the restart, against a per-window threshold. Moving it by fixed steps (a second per run) shows this at once. It also shows that the 6/4/1/1 split depends on the base the variants are added to and on the cut.
- **A noise-free control that "stays" must say why.** Here kernel 3 went the other way. A scaled offset turns the same control into a noise-free collision that ejects kernel 1.
- **A clause that says "implied by X" must check X, not only what X implies.** And a test of a job with two code paths must run both.
