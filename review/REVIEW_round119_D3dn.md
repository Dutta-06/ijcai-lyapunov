# Referee report D3dn, round 119: Proposition D‴.60 as revised after round 114

## Scope

- **Commit.** I reviewed a2f3b77 on claude/vigilant-clarke-hf939b. I exported it with `git -C … archive a2f3b77 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref119.
  - HEAD moved to b2088e8 while I worked (round 118, D‴.58). `git diff a2f3b77 HEAD` touches none of the reviewed files.
- **Read in full:**
  - review/drafts/D60.md (393 lines). The diff since 5be79c1 touches only (iii) at lines 124 and 143–163. The Numerics C8A paragraph (379–384) is unchanged.
  - experiments/exp_noise_resonance.py, every hunk since 5be79c1:
    - `c8a_job`: 'c1i'/'c1p', the p ≥ 2 guard, `kk` from `c_run`, the docstring;
    - `c8a_lines`: `coord`, `agree`, `mono`, `why`, `zd`, `starts`, `near`/`neark`, `early`, `late`, the warm/cold/untested lists, the late-leave clause, the docstring;
    - and what they call: `c_run`, `edge_run`, `win_run`, `offsets_start`, `state_start`, `raised`, `left_periods`, `c_noise`.
  - results/chaos/noise_resonance.md, the C8A paragraph at line 681.
  - tests/test_noise_resonance.py: `test_round105_…` (the j4 wiring and the round-114 plants, line by line), `test_round102_…` (R06's plant) and `test_round89_alarms_…` (R07's plants).
  - tools/mutants/round114_D60_own.py and tools/mutants.py.
  - review/REVIEW_round114_D3di.md, the Round 114 entry of review/REVIEW_LOG.md, and the drafts README line.
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv (symlinked as W/repo/venv);
  - NR_OUT=W/out holding a copy of the raw; TMPDIR=W/tmp;
  - pytest `--basetemp` inside W (through `PYTEST_ADDOPTS` for tools/mutants.py);
  - at most two heavy processes of mine at a time.
- **Nothing in the repository was modified.** I ran only read-only git commands and killed no process.
- **Tags:** [md] the committed md or draft; [raw] the committed raw; [own] my runs or code; [snippet] a web-search result only; [memory] not checked this round.

## Summary

**Design.** The revision answers D3di-1 … D3di-6 in four places.
- **`c8a_job`:**
  - stores the samples' offsets 'c1i' (at pT) and 'c1p' (at (p − 1)T);
  - guards 'warm'/'warm_p' at p ≥ 2;
  - takes `kk` from `c_run`.
- **`c8a_lines`:**
  - partitions each earlier start (k periods back) into early (< R + k − 1), near (±1) and late (> R + k + 1);
  - calls a window cold only if no start is near and none is early, and names early windows "not tested";
  - names late leaves, with the fixed label "the δ(p) start's own drift";
  - `agree` checks the three samples against the scan's removing and keeping offsets;
  - `why` names the reason a window is not one coordinate's; `zd` gives offset 0's distance.
- **Tests:** plants on both sides of the early/near/late boundaries, a sample against the scan on each side, a sample inside the bracket, a flat scan, and the wiring at p = 4.
- **D60:**
  - the headline is now on its base;
  - draw 7's one-period start and its coefficient margin;
  - the late leaves at draw 1 from 50;
  - "the noise (other noise not recorded)".

**What holds**
- **Reproduction is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 8ca842cf…);
  - 80 tests pass (105 s);
  - all 7 C8A jobs, rerun with the script's own `c8a_cases`/`c8a_job`, equal raw['C8A'] in every field but `secs`. This includes 'c1i' and 'c1p' at all 12 windows.
  - My own 64-outputs-a-period run of C's loop gives the stored 'c1i'/'c1p' to 1.3·10⁻¹⁰.
- **Every new number of D60 and the md agrees with the raw** [raw, own recount]: 7, 10, 11; 5, 8, 9, 11; 11 of 12; 0.3104; 0.3000 and 0.3004; "4 periods after the window's start, the alarm's at 2".
- **The three bases are right** [raw]:
  - On 'wc''s base the split is kernel 1's centre's at 10 of 12: 6 via 'w+c1' = 'wc', and 4 via 'wc+i1c' (draws 1 from 203 and 3 from 166 at s = 0.06; draws 1 from 50 and 2 from 32 at 0.08).
  - 'wc+i1w' removes only draw 7. 'wc+i23' removes only draw 2 from 58, where 'wc+i1c' and 'r+c1i' keep the alarm (2).
- **At the 11 alarms at k_e = 20 no earlier start leaves within one period of the alarm or before it** [raw]. Of the 44 starts, 42 pass and 2 are late (draw 1 from 50).
- **The numbers attributed to round 114 reproduce** [own]:
  - draw 7's cold threshold in δ is (0.07635, 0.07637], inside D60's (0.0763, 0.0764]; δ(107) = 0.07667; the run's coefficient at 107T is 0.07613;
  - at draw 2 from 58 the warm loops' centre at 58T is −0.0257 … −0.0268 over 8 starts, inward of (−0.045, −0.0425];
  - at four windows the fine thresholds equal round 114's;
  - own-δ starts at draw 1 from 50 pass (0.2855–0.2859, k = 1…4).
- **The early/near/late boundaries are right on both sides** for all four start types, at R + k − 2 … R + k + 3 (24 plants) [own].

**Findings**
1. **D3dn-1 (minor). Draw 7's "neither" is the δ(p) control's own drift too.**
   - With their own δ, all earlier starts at draw 7 pass (k = 1…4). They arrive at 107T with 0.0754–0.0762, below the cold threshold (0.07635, 0.07637].
   - Only the δ(107) starts from 105 and 104 reproduce the alarm. They arrive at 0.0772 and 0.0780, above it. The one from 106 arrives at 0.0763, below it, and passes.
   - This is the mechanism D60 names "the δ(p) starts' own drift" at draw 1 from 50. Measured with starts that track the run, all 12 alarms are the cold restart's: 11 by the centre (10 of them by kernel 1's centre), draw 7 by its coefficient (the period mean 0.0767 against the run's 0.0761, the threshold between them).
   - The md's "leaves within one period … at 1 of 12 (draw 7)" therefore counts the control, not the alarm. D60's "the earlier δ(107) starts arriving above it" holds for two of the three.
2. **D3dn-2 (minor). "The δ(p) start's own drift" is asserted by the repository, not shown, and its label is unconditional.**
   - The md appends the label to every late leave, including one by 'warm' (the δ(p − 1) start). The authors' own plant s7 prints it for such a start. Two "1 period earlier" entries cannot be told apart.
   - The raw holds no earlier start's state at pT (deferred).
   - My transplant shows the reading is right. At 50T, exchanging kernel 1's coefficient alone between a δ(p) start's state and an own-δ start's state exchanges the outcomes.
   - But the "leaves" are crossings of the cut by 1·10⁻⁵ and 3.6·10⁻⁴, not ejections. The cold restart with the drifted coefficient still produces the alarm at period 2 (1.1967). They do not carry anything "to the alarm".
3. **D3dn-3 (minor). "The three samples on their side of it" is checked only outside the scan's bracket.**
   - 29 of the 33 samples at k_e = 20 are tested by `agree`. 4 lie inside the bracket and are untested by the repository. They are on their side of the finer thresholds (round 114; mine at the 4 windows).
   - `agree` does not order samples against each other inside the bracket. Two contradicting samples there still print "one coordinate" (planted). This is equivalent on the data: at most one sample lies inside the bracket per window.
4. **D3dn-4 (minor). Twelve of my 21 mutants of this round's code survive the full suite:**
   - `agree`'s three boundaries (N01–N03);
   - its pT and (p − 1)T entries (N04, N05): no test varies 'c1i' or 'c1p';
   - two guards and `kk` (N15, N16), equivalent on the data;
   - the late/early handling of warm windows (N10–N12), equivalent on the data;
   - two equivalent by construction (N08, N09).
   - The authors' K09 also survives. Each non-equivalent survivor is killed by one planted row of mine. R06 and R07 are now killed. R11 is equivalent by construction.
5. **D3dn-5 (cosmetic).** Wording and robustness:
   - the Numerics C8A paragraph is not updated;
   - `test_round105_…`'s docstring is stale;
   - `near` is dead code;
   - the late clause's separators and "nan";
   - a row without 'warm' (p = 1) silently drops the whole sentence;
   - p = 0 still wraps;
   - the 'wc'-base identities are not in the md.
6. **D3dn-6 (cosmetic, deferred).** The delayed class at the threshold is confirmed: at draw 2 from 58, offset −0.0425 delays the alarm to period 4 (1.0959) [own]. "One coordinate against a per-window threshold" holds for removal only.
- **No major issue.**

## What I checked and how

All paths are under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref119.

- **Reproduction:**
  - W/logs/fromraw.log, then `cmp` against the committed md: identical.
  - W/logs/pytest_full.log: 80 passed, 1 old warning.
  - W/scripts/repro_c8a.py → W/logs/repro_0246.log and W/logs/repro_135.log: 7 of 7 jobs equal to the raw, field by field, 'c1i'/'c1p' included (9.8–33.6 s each).
- **Recount (numpy only, `python -I`):**
  - W/scripts/dump_c8a.py → W/logs/dump_c8a.log;
  - W/scripts/recount.py → W/logs/recount.log: per window, each earlier start's class, and each sample against the bracket.
- **Planted inputs:**
  - W/scripts/planted.py, W/scripts/planted2.py → W/logs/planted.log: 24 boundary plants, 15 sample and robustness plants, P8 (no offset 0) and two late-label plants.
- **Own probes** (W/scripts/probe.py → W/logs/probe_{d1,d7,d2,brk}.jsonl).
  - My runner uses the repository's `simulate`, `raised`, `c_noise` and `state_start`. The state capture at pT, the period means, the leave and the transplants are my own code.
  - The δ(p) runs reproduce `win_run` and the raw's 'warm_p'/'warm_pk' to all digits.
  - A first d1/d7 pass had a kernel-index bug in the pT coefficient only. It is kept as *_v1_buggy_dw.jsonl and superseded; its leaves and distances are unaffected.
  - **d1** (m = 12, k_e = 20, s = 0.08, draw 1 from 50):
    - δ(p) and own-δ starts, k = 1…4, with state at 50T;
    - transplants at 50T (full state; kernel 1's coefficient replaced);
    - the cold restart with the drifted coefficient;
    - a dose from 48 between δ(48) and 1.5·(δ(50) − δ(48)) above it.
  - **d7** (draw 7 from 107): the cold threshold in δ (7 levels), and δ(p) and own-δ starts, k = 1…4.
  - **d2** (s = 0.06, draw 2 from 58): the warm loops' centre at 58T.
  - **brk**: the fine grid around the 4 samples inside the coarse bracket, and an independent 64-output check of 'c1i'/'c1p'.
- **Mutants:**
  - The authors' 11 (K00 the control) and 21 of my own (W/mut/spec119.py; combined in W/mut/allspec.py) on `-k "round89_md or round98 or round102 or round105"` → W/mut/w1/res_a.jsonl.
  - The 13 survivors plus round 111's R06, R07 and R11 on the full suite → W/mut/wfull/res_a.jsonl and res_b.jsonl.
  - One planted row per survivor, each applied to a private copy of the module: W/scripts/killplants.py and W/scripts/kp_child.py → W/logs/killplants.log.
- **Literature:** three web searches (initialisation shock against drift; exit against escape; the threshold's effect on detected escapes).

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii); (iii) outside 124 and 143–163; Proof; Status | unchanged since 5be79c1; not reverified | [md] |
| md C8A paragraph (md:681) | verified; regenerated byte for byte | [raw], [own] |
| D60:124 "the noise (other noise not recorded)" | verified, matches the md | [md] |
| "On the reduced base kernel 1's centre is one coordinate at 11 of the 12"; 7, 10, 11; 5, 8, 9, 11; none at 0 and +0.04 (143–146) | verified | [raw] |
| "monotone … offset 0 the reduced window itself and the three samples … on their side of it" (147–148) | true. The repository tests 29 of the 33 samples. The 4 inside the bracket hold only on the finer grid: draw 2 from 58 c1 −0.0467 against (−0.045, −0.0425]; draw 3 from 166 c1 −0.0065 against (−0.0175, −0.015]; draw 1 from 50 c1i −0.0600 against (−0.045, −0.0425]; draw 2 from 32 c1 −0.0236 against (−0.0275, −0.025]. D3dn-3 | [raw], [own] |
| "(one threshold per window on a finer grid, between 0.005 and 0.045 … referee round 114 …)" | consistent; the 4 windows I reran give round 114's brackets. At draw 2 from 58, −0.0425 delays to period 4 (D3dn-6) | [own] |
| "on 'wc''s base … kernel 1's centre's at 10 of the 12 — draw 7's … coefficient's, draw 2 from 58's kernels 2 and 3's … (there kernel 1's instantaneous centre keeps the alarm on either base)" (150–151) | verified from the raw. The md names no window for 'wc+i1w'/'wc+i23' (D3dn-5) | [raw] |
| "There the reduced state … never leaves within one period of the alarm's nor before it: those alarms are the cold restart's" (151–153) | verified (44 starts: 42 pass, 2 late, 0 near, 0 early). My own k = 1…4 starts with both δ at two windows agree | [raw], [own] |
| "(by kernel 1's centre at 10 of the 11: … draw 2 from 58 … the warm loop's centre inward …)" (153–156) | verified: −0.0257 … −0.0268 at 58T, inward; at draw 1 from 50 −0.0586 … −0.0604, outward | [own] |
| "The δ(p) starts carry their own drift to the alarm: … leave two periods after … (0.3000 and 0.3004), while starts with their own δ pass" (156–158) | numbers verified. The drift is real; the repository does not show it, I do. The δ(p) starts arrive at 50T 0.0023, 0.0035 and 0.0060 above the run's 0.0930 (k = 2, 3, 4). Own-δ starts arrive 0.0012–0.0019 below it. "To the alarm" overreaches (crossings at period 4, not the alarm's jump). D3dn-2 | [raw], [own] |
| Draw 7: "no offset of the centre removes it (0.3104 …)" | verified | [raw] |
| "started two and three periods earlier with δ(107) … leaves …, from one period earlier it passes" | verified ('warm_p' [−1, 0.0304]; 'warm' passes too) | [raw] |
| "the cold restart's threshold … (0.0763, 0.0764], δ(107) = 0.0767, the run's coefficient 0.0761" | verified: (0.07635, 0.07637], 0.07667, 0.07613 | [own] |
| "the earlier δ(107) starts arriving above it" | true for those from 105 and 104 (0.0772, 0.0780); the one from 106 arrives below (0.0763) | [own] |
| "Draw 7's … is neither" | not supported once the control tracks the run: own-δ starts pass at k = 1…4, so draw 7 is the cold restart's by its coefficient. D3dn-1 | [own] |
| Numerics C8A (379–384) | true but stale: "(rounds 105, 111)", no samples, no late leaves, no one-period pass | [md] |

The new conditions against the readings, with planted inputs [own, W/logs/planted.log]:

| plant | condition | the code's output | correct? |
|---|---|---|---|
| each of 'warm', 'warm_p', k = 2, k = 3 at R + k − 2 / −1 / 0 / +1 / +2 / +3 | early / near / late | not tested / reproduced / reproduced / reproduced / late named / late named | correct, all 24 |
| a sample at max(rem) keeping; at min(keep) removing; between two removing offsets keeping | agree | "a sample against the scan" | correct |
| 'c1i' or 'c1p' outward keeping, or inward removing | agree | "a sample against the scan" | correct |
| a sample inside the bracket, either outcome | agree | one coordinate | correct (no information) |
| two samples inside the bracket in contradicting order (−0.015 keeps, −0.005 removes) | agree | "one coordinate" | **wrong** (equivalent on the data) |
| 'warm' (δ(p − 1)) late | late label | "the δ(p) start's own drift" | **wrong** |
| 'warm' and 'warm_p' both late | late text | two indistinguishable "1 period earlier" entries | ambiguous |
| a row without 'warm' (p = 1) | — | the whole round-105+ sentence disappears | silent |
| no offset 0 in the scan (P8) | label | "offset 0 not the reduced window, its largest distance nan" | no crash; cosmetic |

My own probes at draw 1 from 50 (m = 12, k_e = 20, s = 0.08) [own, W/logs/probe_d1.jsonl]:

| start (state at 50T) | kernel 1's coefficient at 50T | leave (window period) | largest distance |
|---|---|---|---|
| δ(p), k = 1, 2, 3, 4 | 0.0933, 0.0953, 0.0965, 0.0990 | −1, 4, 4, 4 | 0.2875, 0.30001, 0.30036, 0.30091 |
| own δ, k = 1, 2, 3, 4 | 0.0911, 0.0917, 0.0908, 0.0918 | −1 | 0.2853–0.2859 |
| δ(p) k = 2's state, unchanged (restart at 50T) | 0.0953 | 4 | 0.30001 (continuous run: same to 4·10⁻¹³) |
| δ(p) k = 2's state, coefficient := δ(50), own k = 2's, or the run's | 0.0914, 0.0917, 0.0930 | −1 | 0.2857, 0.2858, 0.2870 |
| own-δ k = 2's state, coefficient := δ(p) k = 2's | 0.0953 | 4 | 0.30001 |
| cold restart (centres at the fit), coefficient := 0.0953 | 0.0953 | 2 (the alarm) | 1.1967 |
| from 48, start coefficient δ(48) + f(δ(50) − δ(48)), f = 0 … 1.5 | 0.0917 … 0.0972 | −1 up to f = 0.75, then 4 | 0.2858, 0.2865, 0.2879, 0.2985, 0.30001, 0.3003, 0.3005 |

At draw 7 (m = 12, k_e = 10, s = 0.02) [own, W/logs/probe_d7.jsonl]:
- δ(107) starts arrive with 0.07627, 0.07724, 0.07797, 0.07881 (k = 1…4) and leave at −1, 10, 10, 10.
- Own-δ starts arrive with 0.07544, 0.07551, 0.07552, 0.07616 and all pass (0.0303–0.0304).
- The cold threshold is (0.07635, 0.07637].

Mutants:

| set | killed | survive the full suite |
|---|---|---|
| authors' K01–K10 | K01–K08, K10 (K03 also changes the md) | K09 (warm guard p ≥ 1): equivalent on the data; killed by a wiring plant at p = 1 [own] |
| mine N01–N21 | N06, N07, N13, N14, N17, N18, N19 (J04, now killed by the j4 plant), N20, N21 | N01, N02, N03 (`agree`'s boundaries); N04, N05 (the pT and (p − 1)T samples); N08 (late at alarm + 1) and N09 (a leave in period 0), both equivalent by construction; N10, N11, N12 (warm windows' late or early starts) and N15 (guard p ≥ 3), equivalent on the data; N16 (`kk = 8`), equivalent on the data and on the tests' fake `c_run` (2400 rows) |
| round 111's R06, R07, R11 | R06 (1 failed, 78 passed), R07 (1 failed, 71 passed) | R11 (", not its own" unchecked): equivalent by construction |

## Status of the round-114 issues

| D3di id | status at a2f3b77 | what remains |
|---|---|---|
| D3di-1 (minor) | **fixed.** The headline is on the reduced base; 'wc''s base split is stated (10 of 12, draw 7's coefficient, draw 2 from 58's kernels 2 and 3); the fine thresholds are attributed | the 'wc'-base identities are not printed by the md (D3dn-5) |
| D3di-2 (minor) | **fixed in D60.** The one-period start passes; the threshold, δ(107) and the run's coefficient are given (verified [own]) | the Numerics is unchanged; "the earlier δ(107) starts arriving above it" holds for two of the three; draw 7's "neither" (D3dn-1) |
| D3di-3 (minor) | **partly.** The late leaves are named in the md and D60. Recording each earlier start's state at pT is deferred | the label is unconditional and also covers δ(p − 1); the drift is unshown in the repository (shown here [own]); the leaves are cut crossings (D3dn-2) |
| D3di-4 (minor) | **fixed in substance.** Early leaves are "not tested"; 'c1i'/'c1p' are stored (C8A rerun, bit-identical [own]); the samples are checked against the scan | samples inside the bracket are not tested, nor ordered against each other (D3dn-3); the warm loops' centres are still not in the raw (deferred with D3di-3) |
| D3di-5 (minor) | **fixed for round 114's list.** Plants on both sides; K01–K08 and K10 killed; R06 and R07 killed; J04 killed; L07 removed | K09 survives (equivalent on the data); 12 of my mutants survive, 5 of them not equivalent (D3dn-4) |
| D3di-6 (cosmetic) | **mostly fixed.** `c8a_job`/`c8a_lines` docstrings; "one period earlier and one longer … each to the window's end"; "the noise"; the p ≥ 2 guard; offset 0 by `next` | `test_round105_…`'s docstring is still stale; the p = 0 wrap; **the delayed class at the threshold is deferred and still open**: −0.0425 → period 4 at draw 2 from 58 [own] (D3dn-6) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dn-1 | minor | **Draw 7's earlier-start "reproduction" is the δ(p) start's own drift** (D60:158–163; md:681 "at 1 of 12 (… draw 7 from 107)").<br>• [own] Own-δ starts from 106 … 103 pass, arriving at 107T with 0.07544, 0.07551, 0.07552, 0.07616. The cold threshold is (0.07635, 0.07637].<br>• [own] δ(107) starts arrive with 0.07627 (from 106, passes), 0.07724 (from 105, leaves at 10) and 0.07797 (from 104, leaves at 10).<br>• Same mechanism as D60's draw 1 from 50. With a control that tracks the run, all 12 alarms are the cold restart's; draw 7 by its coefficient (period mean 0.0767 against the run's 0.0761).<br>• "Is neither" contradicts D60's own definition ("the cold restart's" also covers draw 2 from 58, which is not by the centre). | **`c8a_job`:** add own-δ starts two and three periods back ('warm_k', `offsets_start(W, C, p − k, False, (0,))`).<br>**`c8a_lines`:** count a window as reproduced only if an own-δ start reproduces it; name δ(p)-only reproductions "the δ(p) start's".<br>**D60:** "Draw 7's is the cold restart's too, by its coefficient: with their own δ the earlier starts pass (k = 1…4, arriving at 0.0754–0.0762 below the threshold); the δ(107) starts from two and three periods earlier, arriving at 0.0772 and 0.0780, reproduce it by their own drift, as at draw 1 from 50 (referee round 119)". Same in the Numerics. |
| D3dn-2 | minor | **"The δ(p) start's own drift" is a constant label, and "to the alarm" overreaches** (`c8a_lines` late clause; md:681; D60:156–158).<br>• The label is appended to every late leave, including by 'warm' (δ(p − 1)). The authors' plant s7 and my plants print it for a δ(p − 1) start. Two "1 period earlier" entries are indistinguishable.<br>• The raw holds no state at pT.<br>• [own] Shown by transplant: exchanging kernel 1's coefficient alone at 50T between the δ(p) k = 2 and the own-δ k = 2 states exchanges the outcomes (0.30001 against 0.2857–0.2870).<br>• [own] Dose from 48: the largest distance rises monotonically with the start's coefficient and crosses 0.3 just below δ(50) (0.2985 at f = 0.75, 0.30001 at f = 1).<br>• [own] The late leaves are crossings by 6·10⁻⁶ and 3.6·10⁻⁴ without a jump. The cold restart with the drifted 0.0953 still produces the alarm at period 2 (1.1967). | **`c8a_lines`:** print the label only if every late start is a δ(p) start ('warm_p' or k ≥ 2); write "δ(p − 1)"/"δ(p)" in each entry.<br>**`c8a_job`:** record kernel 1's \|w\| at pT for each earlier start (the state at output k·64 in `edge_run`) and print its excess over the run's.<br>**D60:** "carry their own drift to pT (kernel 1's coefficient 0.0023 and 0.0035 above the run's there, referee round 119) and cross the cut by 1·10⁻⁵ and 4·10⁻⁴ two periods after the alarm, without a jump; exchanging kernel 1's coefficient alone at pT between a δ(p) and an own-δ start exchanges the outcomes". |
| D3dn-3 | minor | **The samples are not tested inside the scan's bracket** (`agree`; md:681 and D60:147–148 "the three samples … on their side of it").<br>• [raw] 29 of the 33 samples at k_e = 20 are outside the bracket, all on their side. 4 are inside (draw 2 from 58 c1; draw 3 from 166 c1; draw 1 from 50 c1i; draw 2 from 32 c1), untested by the code.<br>• [own] On the fine grid those 4 are on their side.<br>• [own, planted] Two samples inside one bracket in contradicting order still give "one coordinate". This is equivalent on the data (at most one sample inside the bracket per window). | **`agree`:** merge the 3 samples with the 6 scan points, sort by offset, and require the removing set to be a prefix; print "n of the 3N samples outside the bracket on their side, m inside it".<br>**D60:** "(29 of the 33 outside the scan's bracket on their side, the 4 inside it on their side of the finer thresholds: referee rounds 114, 119)". |
| D3dn-4 | minor | **Survivors of the full suite** [own, W/mut/wfull]:<br>• N01, N02 (`agree`'s ≤/≥ at the bracket's ends);<br>• N03 (max(rem) → min(rem));<br>• N04, N05 (the 'c1i'/'c1p' entries of `agree` dropped): the suite's plants vary only 'c1';<br>• K09, N15 (warm guard p ≥ 1, p ≥ 3);<br>• N16 (`kk = 8`);<br>• N10–N12 (late/early in warm windows);<br>• N08, N09 (equivalent by construction: a start at alarm + 1 makes the window warm; `left_periods` counts from 1).<br>R06 and R07 are killed; R11 is equivalent by construction. | Add the rows of W/scripts/kp_child.py, each verified to kill its mutant [own, W/logs/killplants.log]:<br>• c1 = −0.02 keeping (N01, N03); c1 = 0 removing (N02); c1 = −0.04 keeping (N03);<br>• 'c1i' outward keeping (N04); 'c1p' outward keeping (N05);<br>• the wiring at p = 1 and p = 2 (K09, N15);<br>• a fake `c_run` with 16 outputs a period (N16);<br>• P_N10–P_N12 to pin the intended behaviour of warm windows, or record them as equivalent on the data.<br>Record N08, N09 and R11 as equivalent. |
| D3dn-5 | cosmetic | **Wording and robustness:**<br>• Numerics C8A (D60:379–384) still says "(rounds 105, 111)", with no samples, late leaves or draw 7's one-period pass;<br>• `test_round105_…`'s docstring ("the cold restart's only if no warm start keeps the alarm; one coordinate only if every scan is monotone");<br>• `near` is now unused;<br>• the late clause separates its entries with ", " ("0.3000, 3 periods earlier");<br>• P8 prints "nan";<br>• a row without 'warm' (p = 1, now possible through the guard) drops the whole sentence silently (`all('warm' in x …)`);<br>• p = 0 still wraps (`C[p − 1]`, `W[p − 1]`, `r['c'][(p − 1)·kk]`);<br>• D60's 'wc'-base identities are not printed by the md. | Update the Numerics and the docstring; drop `near`; use "; " between entries; require 'warm' per row, or print "not tested (no earlier start)"; guard p ≥ 1 in `c8a_cases`; print the window each of 'wc+i1c'/'wc+i1w'/'wc+i23' removes. |
| D3dn-6 | cosmetic | **The delayed class at the threshold (deferred)** [own]: at draw 2 from 58, −0.045 removes, −0.0425 delays to period 4 (1.0959) and −0.04 keeps (2). The coordinate's outcome there has three values. | One clause in D60: "a threshold for removal; at draw 2 from 58 a delayed band between it and the kept offsets". |

## Sources

- D. Mulholland, P. Laloyaux, K. Haines, M. Balmaseda, "Origin and impact of initialisation shocks in coupled atmosphere-ocean forecasts", Mon. Wea. Rev. 143(11), 4631–4644 (2015) [snippet].
  - The paper separates shock (from an imperfect initialisation method) from drift, which appears whatever initialisation is used.
  - The δ(p) control is an imperfect initialisation (a period-p coefficient placed at p − k), and its drift is what decides draw 1's late leaves and draw 7's reproductions (D3dn-1, D3dn-2).
  - Links: [Reading eprint](https://reading-9.eprints-hosting.org/40638/1/mulhollandetal2015_initshock_revise2.pdf), [CentAUR record](https://centaur.reading.ac.uk/40638).
- L. Zonca, D. Holcman, "Exit versus escape in a stochastic dynamical system of neuronal networks explains heterogenous bursting intervals", arXiv:2010.06699 [snippet].
  - Reaching the basin's boundary is not escape. The durations do not match the time to reach the boundary, and the snippet reports an escape time 2–3 times the first exit time.
  - Analogy: the late "leaves" at 0.30001 and 0.30036 are exits across the cut, not the alarm's ejection (1.196) (D3dn-2).
  - Link: [arXiv 2010.06699](https://arxiv.org/abs/2010.06699v2).
- "Interplay between Detection Strategies and Stochastic Resonance Properties", arXiv:1505.07646 [snippet].
  - The threshold used to decide that a passage occurred changes the escape-time statistics. A margin ±θ is used so that crossings near the boundary are not counted.
  - It supports reporting the late leaves' margin over the cut.
  - Link: [arXiv 1505.07646](https://arxiv.org/pdf/1505.07646).
- P. Schultz, F. Hellmann, K. N. Webster, J. Kurths, Chaos 28, 043102 (2018) [snippet], seen again in the second search. It is consistent with D60's citation.
  - Link: [arXiv 1711.03857](https://arxiv.org/pdf/1711.03857).
- Not rechecked this round [memory or earlier rounds' snippets]: the remaining citations of D60:232–254.

## Lesson

- **A control's artefact is an artefact wherever it acts.** The δ(p) start's drift was named where it gave an unexpected outcome (draw 1's late leave), but read as evidence where it gave the expected one (draw 7's reproduction). Run the control that tracks the run (its own δ) beside it, and classify by that.
- **A conclusion printed by code must be computed by code.** "The δ(p) start's own drift" is a fixed string attached to any late leave, a δ(p − 1) one included. One transplant at pT (exchange one coordinate, keep the rest) would have shown it.
- **A check that cannot fail inside its bracket is not a check there.** Four of the 33 samples were "on their side" only on the referee's finer grid. Merge the samples into the scan.
- **A leave by 10⁻⁵ over the cut is a crossing, not an ejection.** Report the margin over the cut along with the period.
