# Referee report D3ds, round 124: Proposition D‴.60 as revised after round 119

## Scope

- **Commit.** I reviewed f346010 on claude/vigilant-clarke-hf939b.
  - Exported with `git -C … archive f346010 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref124. The repository's venv is symlinked as W/repo/venv.
  - HEAD moved to 5f56b63 while I worked (round 122, D‴.61). `git diff f346010 HEAD` touches none of the reviewed files.
- **Read in full:**
  - review/drafts/D60.md (398 lines; the diff since a2f3b77 touches (iii) 151–165 and the Numerics' C8A, 382–389).
  - experiments/exp_noise_resonance.py:
    - `c8a_job`: 'warm_k';
    - `c8a_lines`: `merged`, `agree`, `inside`, `starts`, `kname`, warm/dpo/cold/unt and the late clause;
    - what they call: `offsets_start`, `win_run`, `edge_run`, `c_run`, `state_start`, `raised`, `c8d_alarms`, `left_periods`.
  - results/chaos/noise_resonance.md, the C8A paragraph.
  - tests/test_noise_resonance.py: `test_round105_…` with the round-114 and round-119 plants, `test_round98_c8a_job_wiring`, `test_round102_…`, `test_round89_md_is_the_report…`.
  - tools/mutants/round119_D60_own.py and tools/mutants.py.
  - review/REVIEW_round119_D3dn.md, the Round 119 entry of review/REVIEW_LOG.md, and §5 of HANDOFF.md.
  - `git diff a2f3b77 f346010` for these files: commit 1ad8abb holds the code and md, cb786a1 the plants and D60's Numerics.
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv;
  - NR_OUT=W/out holding a copy of the raw; TMPDIR=W/tmp; pytest `--basetemp` inside W;
  - at most two heavy processes of mine at a time.
- **Nothing in the repository was modified.** I ran only read-only git there.
  - I killed no process. The only processes I stopped were my own two monitor tasks, after checking their cmdlines.
  - A modification of experiments/exp_width_mismatch.py appeared in the real working tree during my run. It is not mine.
- **Tags:** [md] the committed md or draft; [raw] the committed raw; [own] my runs or code; [snippet] a web-search result only; [memory] not checked this round.

## Summary

**Design.** The revision answers D3dn-1 … D3dn-6.
- **`c8a_job`** adds 'warm_k': the reduced state started two and three periods back with its own δ, `offsets_start(W, C, p − k, False, (0,))`.
- **`c8a_lines`:**
  - counts a window as reproduced only by an own-δ start ('warm', 'warm_k');
  - names δ(p)-only reproductions in `dpo` with the fixed text "the δ(p) start's own drift, not the alarm";
  - calls a window cold unless an own-δ start leaves near or before the alarm;
  - merges the three samples into the scan (`merged`; `agree` asks that the removing set be a prefix) and counts the samples strictly inside the bracket (`inside`);
  - labels each late leave by its δ (`kname`), and prints the drift label only if every late leave is a δ(p) start.
- **Tests:** plants for each change, and the wiring with 'warm_k'.
- **D60 (iii) and Numerics:** with its own δ, 0 of 12 reproduced; "all 12 alarms are the cold restart's — 11 by the centre, 10 of them by kernel 1's —, draw 7's by its coefficient"; round 114's "neither" withdrawn.

**What holds**
- **Reproduction is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 87803ec7…);
  - 80 tests pass (109.6 s);
  - all 7 C8A jobs, rerun with the script's own `c8a_cases`/`c8a_job`, equal raw['C8A'] in every field but `secs`, 'warm_k' included (21.6–75.9 s each).
  - My own two-segment runner reproduces 'warm', 'warm_p', 'warm_k' and 'warm_pk' at draw 7 and at draw 2 from 58: the same leave periods, largest distances to ≤ 1.1·10⁻¹².
- **Every new number of the md and D60 agrees with the raw** [raw]:
  - all 36 own-δ starts pass (0 of 12);
  - draw 7's δ(p) starts: [2, 12, 0.31014] and [3, 13, 0.30998], at R + k (1 more);
  - draw 1 from 50's δ(p) starts: [2, 6, 0.30001] and [3, 7, 0.30036], "4 periods after the window's start", late;
  - 4 of the 36 samples are inside their bracket: draw 2 from 58 c1 −0.0467; draw 3 from 166 c1 −0.0065; draw 1 from 50 c1i −0.0600; draw 2 from 32 c1 −0.0236;
  - "106 … repeats it": 106 neighbours, 0 repeats.
- **The δ(p) and own-δ arrivals that D60 attributes to round 119 reproduce** [own]:
  - at 107T the own-δ starts arrive with 0.07544, 0.07551, 0.07552 and the δ(p) starts with 0.07627, 0.07724, 0.07797;
  - the cold threshold in δ is (0.076350, 0.076351], inside (0.0763, 0.0764].
- **The authors' Y01–Y10 are all killed**, and so are 18 of my 23 mutants. No surviving mutant changes the published results.
- **No major issue.**

**Findings**
1. **D3ds-1 (minor; the leading finding). At draw 7 the own-δ control's pass is its arrival deficit, so "all 12 alarms are the cold restart's" and the withdrawal of "neither" do not hold there** [own].
   - An own-δ (or δ(p)) earlier start's state at 107T, with kernel 1's coefficient alone set to δ(107), leaves at the alarm's period 10 (0.31036, against the cold restart's 0.31037). This holds at k = 1, 2, 3 for both kinds of start.
   - Setting the warm state's kernel 1 centre, kernels 2 and 3, the unformed kernels or x to the cold values changes nothing.
   - The warm state's threshold, (0.076337, 0.076338], is 1.3·10⁻⁵ below the cold restart's, (0.076350, 0.076351]. The cold restart started with the own-δ arrivals passes.
   - So the own-δ starts pass because they arrive 0.0012 below δ(107). This is the mirror image of the δ(p) starts from 105 and 104, which reproduce the alarm by arriving above it. Neither control isolates the restart.
   - Draw 7's alarm is the coefficient's (δ(107), a period mean, placed at 107T), not the restart's. Round 114's "neither" was right.
   - The code prints "the δ(p) start's own drift, not the alarm" as a fixed string whenever `dpo` is non-empty (planted P5). It computes nothing.
2. **D3ds-2 (minor). The own-δ starts are C8D's dense windows at p − k, run k periods longer** [raw].
   - At all 36 starts the leave is identical. At 28 of them the largest distance agrees to < 10⁻¹²; 2 agree to 8 digits; 6 reach their maximum in the extra periods. C3's δ equals C8's (difference 0).
   - So at the 11 windows at k_e = 20 the new count restates D60:135's "none of the 106 other dense windows … repeats it".
   - It adds information only at draw 7, whose alarm period 117 no earlier dense window covers. D60 presents it as a separate control.
   - δ(p − k) is a period mean on the reduced state (kernel 1's centre, kernels 2 and 3 and x at the fit). It is not "the run's own state", which is what the round-119 log's lesson asks for.
3. **D3ds-3 (minor). The md no longer prints "one coordinate against a per-window threshold"** [md].
   - The clause is printed only if `all((r, x) in nco for r, x in cold)`. Draw 7, which is not in `nco`, has moved into `cold`, so the md at f346010 lost the sentence it printed at a2f3b77.
   - D60:144–149 still asserts it. The counts (11 of 12) still support it.
   - No test pins it, because the md test compares against the regenerated md.
4. **D3ds-4 (minor). "11 by the centre, 10 of them by kernel 1's": the 11th centre is kernel 3's, shown only by my transplant** [own].
   - At draw 2 from 58, setting kernel 3's centre alone to the fit (warm at +0.038 … +0.040) makes the warm state leave at period 2 (0.8538). Kernel 1's coefficient and centre, kernel 2's centre, kernels 2 and 3's coefficients, the unformed kernels and x do not.
   - At draw 1 from 50, draw 3 from 193 and draw 2 from 6, kernel 1's centre alone does it (1.196, 1.893, 1.124).
   - D60 names no centre for the 11th window, and the repository computes none.
5. **D3ds-5 (minor). Draw 7's "marginal crossing (0.3104 against 0.3)" is an ejection that the window's end cuts off** [own].
   - Run one period longer, the cold restart's loop reaches 0.822 at period 11.
   - At draw 7 (R = 10 = NP_C7) the earlier starts end at the alarm's period, so "within one period after it" is not run. I ran it: the own-δ starts still pass at period 11 (0.0303).
6. **D3ds-6 (minor). Five of my mutants survive the full suite** [own]:
   - Z08 (`dpo` without the warm exclusion), Z09 (the late label with `any`), Z12 (`inside` without 'c1p'): equivalent on the data, each killed by one planted row;
   - Z21, Z22: equivalent by construction;
   - none changes the md.
7. **D3ds-7 (cosmetic). Wording and robustness:**
   - "(and one period longer)" is wrong for k = 2, 3, which run k periods longer;
   - `kname` takes the first match, so two starts with identical outcomes share one label, and the drift label can print over an own-δ start (P1, P2);
   - "4 of the 36" counts draw 7's 3 samples, which have no bracket;
   - D60 (iii) is not updated with the 4 inside;
   - the docstrings of `c8a_job` and `c8a_lines` do not mention 'warm_k' or round 119;
   - "round 105's own count had generalised it too" is opaque;
   - "0.3000 and 0.3004 against 0.3" hides the margins (6·10⁻⁶, 3.6·10⁻⁴);
   - the scope of "(referee round 119, not reproduced …)" is ambiguous: the arrival "above the threshold" is also round 119's, and the coefficient exchange was done only at draw 1 from 50.

## What I checked and how

All paths are under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref124.

- **Reproduction:**
  - W/logs/fromraw.log, then `cmp` against the committed md: identical.
  - W/logs/pytest_full.log: 80 passed, 1 old warning.
  - W/scripts/repro_c8a.py → W/logs/repro_0.log and W/logs/repro_1.log: 7 of 7 jobs, 0 differences.
- **Recount and identity (numpy only, `python -I`):**
  - W/scripts/dump_c8a.py → W/logs/dump_c8a.log;
  - W/scripts/c8d_vs_warm.py → W/logs/c8d_vs_warm.log: own-δ starts against C8D's windows at p − k, C3's δ against C8's, and the neighbours covering each alarm.
- **Own probes.** W/scripts/probe.py, probe2.py and probe3.py use the repository's `simulate`, `raised`, `alpha_of`, `offsets_start` and `c_noise`. The two-segment runner, the state capture at pT, the transplants and the bisection are my code.
  - **probe_d7.jsonl, probe_d2.jsonl:** per start, kernel 1's |w| − w₀ and the centres at pT. Transplants: one component of the warm state set to the cold restart's, singly and cumulatively. The cold restart with each warm arrival coefficient. Each run one period beyond the window.
  - **probe2_{d1,d2,d3,d6}.jsonl:** the split by kernel and component on own-δ k = 1, 2 at draws 1 from 50, 2 from 58, 3 from 193 and 2 from 6.
  - **probe3_d7.jsonl:** 9 bisections of the coefficient threshold on the cold state and on the own-δ k = 2 state.
- **Plants:** W/scripts/plants.py → W/logs/plants.log, P1–P7.
- **Mutants.** The authors' Y00–Y10 plus my Z01–Z23 (W/mut/spec124.py; combined in W/mut/allspec.py):
  - on `-k "round89_md or round98 or round102 or round105"` → W/mut/k1/res_{a,b}.jsonl;
  - the survivors and the control on the full suite → W/mut/full/res_{a,b}.jsonl;
  - kill plants on private mutated copies: W/scripts/killcheck.py → W/logs/killcheck.log.
- **Literature:** four web searches. Fetching ntrs.nasa.gov and www.ecmwf.int failed (the hosts are not allowed in this session), so I cite snippets only.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii); (iii) outside 151–165; Proof; Status | unchanged since a2f3b77; not reverified | [md] |
| md C8A paragraph | regenerated byte for byte; every new number recounted. It no longer prints the one-coordinate conclusion (D3ds-3) | [raw], [own] |
| 135 "none of the 106 other dense windows … repeats it" | verified (106, 0). At draw 7 all 9 neighbours start after 107 | [raw] |
| 151–153 "With its own δ … never leaves within one period … nor before it, at any of the 12" | numbers verified. These starts are C8D's windows at p − k (D3ds-2); "(and one period longer)" is wrong for k = 2, 3 (D3ds-7); at draw 7 the period after the alarm is not run (my run: they pass, D3ds-5) | [raw], [own] |
| 153 "all 12 alarms are the cold restart's" | follows from the repository's counts by its own definition. Confirmed by transplant at 4 windows at k_e = 20. Refuted at draw 7: a warm state with δ(107) leaves at period 10 (D3ds-1) | [raw], [own] |
| 153 "11 by the centre, 10 of them by kernel 1's" | true where I tested: kernel 1's centre at 3 windows, kernel 3's at draw 2 from 58. Kernel 3 is not named; 7 of the 10 rest on round 114's centre positions (D3ds-4) | [own] |
| 153–156 "(draw 2 from 58 … warm loop's centre inward …; removed by kernel 1's mean centre on the reduced base, only delayed on the 'w' base)" | verified: −0.0257 … −0.0268; 'r+c1m' −1, 'w+c1' 5 | [own], [raw] |
| 156–157 "draw 7's … by its coefficient: no offset of the centre removes it" | verified and sharpened: the coefficient decides on the warm and the cold state alike | [raw], [own] |
| 157 "(0.3104 against the cut 0.3: a marginal crossing)" | 0.3104 verified; "marginal" is not: 0.822 one period later (D3ds-5) | [raw], [own] |
| 157–159 the threshold (0.0763, 0.0764], δ(107) = 0.0767, the run's coefficient 0.0761 | consistent: the cold threshold is (0.076350, 0.076351] | [own] |
| 159–161 "δ(p) starts carry their own drift to pT … reproduce its alarm (arriving above; from one period earlier below, passing)" | arrivals verified (0.07724, 0.07797; 0.07627). But a warm state carrying exactly δ(107) reproduces the alarm too, so the reproduction is not the control's artefact (D3ds-1). The arrivals are not in the repository | [own] |
| 161–162 draw 1 from 50 "cross the cut two periods after … (0.3000 and 0.3004 …, without a jump)" | verified; the margins are not printed (D3ds-7) | [raw] |
| 162–163 "exchanging kernel 1's coefficient alone … exchanges the outcomes" | round 119 did this at draw 1 from 50 only. Same direction at draw 7 [own]. At draw 2 from 58 the coefficient decides nothing | [own] |
| 163–165 "So a δ(p) start is the control's artefact wherever it acts — round 114's 'neither' counted the control" | not supported (D3ds-1) | [own] |
| 144–149 "kernel 1's centre is one coordinate at 11 of the 12 … one coordinate against a per-window threshold" | counts verified (`nco` = 11). The md no longer prints the conclusion (D3ds-3). The 4 samples inside the bracket are not mentioned in (iii) | [raw], [md] |
| Numerics C8A (382–389) | numbers match the md. "That start's own drift" inherits D3ds-1 | [md] |

## Status of the round-119 issues

| D3dn id | status at f346010 | what remains |
|---|---|---|
| D3dn-1 (minor) | **implemented as proposed:** 'warm_k', reproduction counted by own-δ starts only, `dpo`, "neither" withdrawn; reproduced bit-identically | the proposed reading is refuted at draw 7 by a matched-arrival transplant; the own-δ control is confounded too (D3ds-1) |
| D3dn-2 (minor) | **mostly fixed:** each entry names its δ; "; " separators; the label conditional; "to the alarm" removed; the transplant attributed | the margins are not printed; kernel 1's \|w\| at pT is still not recorded (deferred); a new fixed string in the `dpo` clause; `kname`'s first match (D3ds-7) |
| D3dn-3 (minor) | **fixed in the code and md:** `merged`, `agree` as a prefix; "4 of the 36 inside" verified | D60 (iii) is not updated; draw 7's 3 bracketless samples are counted in the 36 |
| D3dn-4 (minor) | **fixed:** the N01–N05 plants kill their mutants; Y01–Y10 killed; K09, N08–N12, N15, N16, R11 recorded | my Z08, Z09, Z12 survive, each with a kill plant; Z21, Z22 equivalent by construction (D3ds-6) |
| D3dn-5 (cosmetic) | **partly fixed:** `near` dropped; the test docstring and D60's Numerics updated | "nan", a row without 'warm', p = 0's wrap and the 'wc'-base identities deferred (still open); the `c8a_job`/`c8a_lines` docstrings are now stale |
| D3dn-6 (cosmetic) | deferred, open | — |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ds-1 | minor | **Draw 7's verdict is the own-δ control's arrival deficit** (D60:153, 156–165; md "the δ(p) start's own drift, not the alarm"; `c8a_lines`' `dpo` string).<br>• [own, probe_d7] Warm states at 107T (own-δ and δ(p), k = 1, 2, 3) with kernel 1's coefficient := δ(107) = 0.07667 leave at period 10 (0.31036). Kernel 1's centre, kernels 2 and 3, the unformed kernels or x := cold change nothing.<br>• [own, probe3] The warm threshold is (0.076337, 0.076338], the cold one (0.076350, 0.076351].<br>• [own] The cold restart started with the own-δ arrivals (0.07544–0.07552) passes; with 0.07724 and 0.07797 it leaves.<br>• So the own-δ pass is the arrival 0.0012 below δ(107), the mirror of the δ(p) starts' excess. The alarm is the period-mean coefficient's, not the restart's: round 114's "neither" stands.<br>• [planted P5] The `dpo` reading is a fixed string. | **`c8a_job`:** keep each earlier start's state at pT. Record kernel 1's \|w\| − w₀ and the formed centres. Run 'warm_m': that state with kernel 1's coefficient := w₀ + δ(p), NP_C7 periods under the noise from p.<br>**`c8a_lines`:** call a window the cold restart's only if 'warm_m' passes for every own-δ k, and "the coefficient's" if it reproduces the alarm. Compute the `dpo` sentence from the arrivals.<br>**D60:153–165 and Numerics:** "11 of the 12 are the cold restart's; draw 7's is the coefficient's, not the restart's: the earlier starts' states at 107T with kernel 1's coefficient set to δ(107) leave at the alarm's period (k = 1, 2, 3, either δ; the warm threshold (0.076337, 0.076338] against the cold (0.076350, 0.076351]); the own-δ starts pass by arriving 0.0012 below δ(107), as the δ(p) starts from 105 and 104 reproduce it by arriving above (referee round 124)". Delete "So a δ(p) start is the control's artefact wherever it acts …" and restore "neither".<br>**REVIEW_LOG:** amend round 119's lesson. |
| D3ds-2 | minor | **The own-δ starts are C8D's dense windows at p − k** [raw, c8d_vs_warm.log].<br>• Same leaves at 36 of 36; largest distances equal to < 10⁻¹² at 28, 8 digits at 2, maximum in the extra periods at 6. C3's δ equals C8's.<br>• At the 11 at k_e = 20 the count restates D60:135's 106/0. Only draw 7 (period 117, no earlier dense window) gains anything.<br>• δ(p − k) is a period mean on the reduced state, not the run's own state. | **D60:151:** "the earlier starts with their own δ are C8D's windows at p − 1 … p − 3 run k periods longer; at the 11 at k_e = 20 this restates 'none of the 106 …'; it adds draw 7, whose alarm in the window's last period no earlier dense window covers".<br>**Tests:** assert the identity of 'warm'/'warm_k' with C8D's `wins` at p − k.<br>**`c8a_job` docstring.**<br>**REVIEW_LOG:** "its own δ (C8's rule at p − k), not the run's state". |
| D3ds-3 | minor | **The one-coordinate conclusion vanished from the md** (md C8A at a2f3b77 against f346010; D60:144–149).<br>• The gate `all((r, x) in nco for r, x in cold)` fails once draw 7 is in `cold`.<br>• No test pins the sentence on the real raw. | **`c8a_lines`:** print "at the {len(nco)} … kernel 1's centre at the restart is one coordinate against a per-window threshold" from `nco`, independent of `cold`. Regenerate the md.<br>**Test:** assert the clause in `report(raw)`. |
| D3ds-4 | minor | **The 11th centre is kernel 3's** (D60:153) [own, probe2].<br>• At draw 2 from 58, c3 := fit alone → leaves at period 2 (0.8538). c2, w2, w3, the unformed kernels, kernel 1's w and c, and x → pass.<br>• At draws 1 from 50, 3 from 193 and 2 from 6, c1 := fit alone → the alarm (1.196, 1.893, 1.124).<br>• Not named in D60; not computed by the repository. | **D60:153:** "11 by a centre — kernel 1's at 10 (round 114's warm centres outward of the threshold; by transplant at draws 1 from 50, 3 from 193, 2 from 6), kernel 3's at draw 2 from 58 (warm at +0.038 … +0.040; referee round 124)".<br>**Optionally `c8a_job`:** a variant 'r+c3i' (the reduced start with kernel 3's centre at pT). |
| D3ds-5 | minor | **"A marginal crossing" is a truncated ejection** (D60:157) [own].<br>• The cold restart at draw 7 is at 0.822 one period after the window.<br>• At R = NP_C7 the earlier starts are not run past the alarm's period. Own run: they pass at period 11 (0.0303). | **D60:157:** "0.3104 in the window's last period, 0.822 one period later (referee round 124)".<br>**`c8a_job`:** run the earlier starts to p + R + 2, or print "the period after the alarm not run". |
| D3ds-6 | minor | **Survivors of the full suite** [own, W/mut/full]:<br>• Z08 (`dpo` without the warm exclusion), Z09 (`any` in the late label), Z12 (`inside` without 'c1p'): equivalent on the data;<br>• Z21 (`cold`'s early test with ≤), Z22 (`dpo` over every start): equivalent by construction;<br>• none changes the md. | **Tests:** add the rows that kill them [own, W/logs/killcheck.log]: P4 (an own-δ and a δ(p) start both near, against Z08), P3 (mixed late leaves, against Z09), P6 ('c1p' inside the bracket, against Z12).<br>**Mutant record:** Z21 and Z22 as equivalent by construction. |
| D3ds-7 | cosmetic | **Wording and robustness:**<br>• "(and one period longer)" (md, D60:152; k periods);<br>• `kname`'s first match (P1: two "δ(p − 1)" labels; P2: "own drift" printed over an own-δ start);<br>• "4 of the 36" includes draw 7's 3 bracketless samples;<br>• D60 (iii) does not mention the 4 inside;<br>• stale docstrings of `c8a_job`/`c8a_lines`;<br>• "round 105's own count had generalised it too";<br>• margins hidden by "0.3000 … against 0.3";<br>• the scope of round 119's attribution;<br>• the D3dn-5 deferrals. | **`kname`:** match by index, not by value. Print "of the 33 in a bracket". Write "k periods longer". Margins: "by 6·10⁻⁶ and 3.6·10⁻⁴". Attribute the arrivals and the exchange per window. Update the docstrings; clarify or drop the round-105 aside. |

## Sources

- **R. N. Hoffman and E. Kalnay (1983), lagged average forecasting** [snippet]. The search snippets say the lagged ensemble's initial states are the current analysis plus forecasts from earlier analyses verifying at the current time.
  - That is the own-δ design: each earlier start from its own state. The design compares starts that do not share an arrival state (D3ds-1, D3ds-2).
  - Seen through the ECMWF pages "Representing uncertainty in initial conditions: experimentation from lagged average forecasting to singular vectors" (https://www.ecmwf.int/node/28627) and the historical-perspective PDF (https://www.ecmwf.int/sites/default/files/elibrary/2017/17604-historical-perspective-earlier-ensembles-and-forecasting-forecast-skill-and-global-leading.pdf). Neither could be fetched.
- **"A Practical Probabilistic Benchmark for AI Weather Models"**, arXiv:2401.15305 [snippet]. Lagged members "are not interchangeable since the more recent members are more skillful". Relevant to D3ds-1: the earlier starts differ from the cold restart by their age-dependent arrival. https://arxiv.org/pdf/2401.15305
- **U. Stein and P. Alpert's factor separation** [snippets]:
  - Cambridge, *Factor Separation in the Atmosphere* (introduction page): one-at-a-time sensitivity tests "can lead to the incorrect answer";
  - S. O. Krichak and P. Alpert, JAS 2002: "four model simulations are sufficient to evaluate the role of two chosen factors and that of their synergistic interaction";
  - my transplants are this method applied at pT.
  - Links: https://www.cambridge.org/core/books/abs/factor-separation-in-the-atmosphere/introduction/8BE10EB67765B58AF3CDF24964A450BC and https://www.tau.ac.il/~pinhas/papers/2002/Krichak_Alpert_JAS_2002.pdf
- **D. Mulholland et al., Mon. Wea. Rev. 143 (2015)** [snippet, seen again]: initialisation shock against drift. https://reading-9.eprints-hosting.org/40638/1/mulhollandetal2015_initshock_revise2.pdf
- Not rechecked this round [memory or earlier rounds' snippets]: D60's own citations (232–254).

## Lesson

- **Match the arrival before blaming the restart.** Both earlier-start controls borrow a coefficient: the δ(p) start from a later period, the own-δ start from an earlier one. Each drifts its own way to pT. At draw 7 one arrives above the threshold and the other below, so neither isolates the restart.
- **The test that does is a transplant at pT** (factor separation): the warm state with the cold restart's value of one coordinate at a time. One such run would have shown round 119's reading and round 114's "neither" for what they are.
- **A published sentence gated on a set inclusion** (cold ⊆ nco) switches off silently when an unrelated reclassification changes the set. Pin each printed conclusion with a test on the real raw, not only on plants.
