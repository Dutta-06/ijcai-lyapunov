# Referee report D3ec, round 134: Proposition D‴.60 as revised after round 128

## Scope

- **Commit.** I reviewed 7d78bfa on claude/vigilant-clarke-hf939b. Between a15cc7e (round 128 applied) and 7d78bfa, none of D60.md, exp_noise_resonance.py, test_noise_resonance.py, noise_resonance.md, its raw or round128_D60_own.py changed (`git diff --stat`).
- **Export.** `git -C … archive 7d78bfa | tar -x -C W/repo`, with W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref134 and the venv symlinked as W/repo/venv.
- **How I ran things.**
  - Environment: `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`; NR_OUT=W/out; TMPDIR=W/tmp; LCK_BUILD=W/build.
  - Heavy processes ran one at a time, in sequence: tests, then reruns, then the factorial, then the mutants.
  - Nothing in the repository was modified; I used only read-only git there.
  - I killed no process. The pytest processes seen at the end (test_moving_width, `vver`) belong to other sessions.
- **Read in full:**
  - review/drafts/D60.md: (iii), the Proof, Status and the Numerics' C8A.
  - review/REVIEW_round128_D3dw.md and the Round 128 entry of review/REVIEW_LOG.md.
  - `git diff b361d6d a15cc7e` for the four files.
  - tools/mutants/round128_D60_own.py.
  - The code these changes call: `seg_run`, `c8a_pT`, `c8a_job`, `bracket_out`, `c8a_lines` (the pT block), `raised`, `edge_run`, `win_run`, `offsets_start`, `state_start`, `c_run` and `left_periods`.
  - The D60 row of review/drafts/README.md.
- **Tags:** [md] the committed md or draft; [raw] the committed raw; [own] my runs or code; [snippet] a web-search result only, no page fetched; [memory] not checked this round.

## Summary

**Design.** Round 128 added four things:
- **'cold_c'**: the cold restart with one own-δ start's formed centre at its value at pT, giving 9 runs per window;
- **'w1i'**: the run's own |w₁| − w₀ at pT;
- **`bracket_out`**: brackets rounded outward;
- **`c8a_lines` changes**: the matched arrival at every window that has one, "mixed over k", the margins, the 2² clause, "necessary and sufficient", and the header taken from `C8A_ROUNDS`.

D60 (iii) now reads:
- 11 of the 12 alarms are the cold restart's, each by one centre that is necessary and sufficient: kernel 1's at 10 windows, kernel 3's at draw 2 from 58.
- Draw 7's alarm is the coefficient's.

**Verdicts**
- **Reproduction is exact** [own]:
  - `--from-raw` gives the md byte for byte (md5 1a639156…);
  - the 87 tests pass (128.6 s);
  - all 7 C8A jobs, rerun with the script's own `c8a_cases`/`c8a_job`, equal raw['C8A'] in every field but `secs`.
- **'cold_c' changes exactly one centre and nothing else** [own, by reading the code]:
  - `cw = c0` with `cw[idx[j]] = c[idx[j]]`;
  - the cold restart's own α (`alc`) and the plant at x_m(0) (`x0=None`);
  - the same `c_noise(…, shift=p)` as the cold restart.
- **"Necessary and sufficient" follows from 'warm_c' and 'cold_c' together** at all 11 windows and all 3 own-δ k [raw]:
  - the named centre flips the outcome in both directions;
  - each other formed centre flips it in neither.
- **My full 2⁵ factor separation confirms the reading and closes its one gap** [own]. The gap: the repository moves only the three formed centres one at a time.
  - Setup: 36 starts × 32 corners = 1152 runs. The five factors are α₁, the named centre, the other two formed centres, the plant's x, and every other α together with the unformed kernels' centres.
  - Result: the outcome is a function of the named centre alone at all 33 starts at the 11 windows, and of α₁ alone at draw 7's 3 starts.
  - Stein–Alpert terms on the largest distance (by Möbius inversion): the named factor carries +1.000 … +1.001 of the total change, and every other main or interaction term is at most 0.001 in magnitude.
- **Draw 7's reading holds**:
  - all six 'warm_m' runs leave at the alarm's period 10;
  - the three own-δ 'cold_m' runs pass;
  - all 9 'cold_c' runs leave at period 10;
  - the outcome goes with α₁ at all 32 corners.
- **Every changed number of D60 and the md agrees with the raw** [raw]. No overreach found.
- **No major issue.**

**Findings**
1. **D3ec-1 (minor, latent). `bracket_out` rounds outward from the float's shortest repr, not from the double's exact value.**
   - Whenever repr(v) has at most 7 decimals, the printed end is inward of the true double in about half the cases. Of the 100,001 doubles nearest the 7-decimal grid in [0.07, 0.08], 50,000 would print lo inward and 49,999 would print hi inward [own].
   - On the data all four printed brackets contain the exact doubles [own].
   - But the tests pin the inward output in two plants: "(0.0763377, 0.0763379]" and "(0.0763300, 0.0763400]".
   - The authors' mutant H12, `Decimal(float(v))`, is the correct implementation, and it is "killed" only by those plants.
2. **D3ec-2 (minor, latent). Seven semantic survivors of my 22 mutants in round 128's printing logic.**
   - The survivors are R04–R10. All are equivalent on the data and none changes the md.
   - Six of them (R05–R10) are killed by my plants P05–P10.
3. **D3ec-3 (cosmetic). The scope of "necessary and sufficient".**
   - What the repository computes is two-factor: the centre against the rest of the state as a block, with only the formed centres moved one at a time.
   - My 2⁵ shows no hidden interaction. The text should state the scope, or the 2⁵ should be added (about 10 s per start).
4. **D3ec-4 (cosmetic). The α₁ contrast of the first 2² design is not printed.**
   - At the 11 windows |arrival − δ(p)| = 1.2·10⁻⁶ … 1.7·10⁻³; at some k that cell of the design is near-degenerate.
   - 'w1i' at the 11 windows (0.00055 … 0.00182 above δ(p)) is stored but not printed.
5. **D3ec-5 (cosmetic). 'cold_c' "passes" means "past the alarm's period" (`ps_`).**
   - At draw 1 from 50 the passing corners reach 0.2857 in the window's period 4, 0.014 below the cut. This is the same excursion that the δ(p) starts cross (0.3000, 0.3004).
6. **D3ec-6 (cosmetic). History and record residues:**
   - the README's round-124 clause still says "round 114's 'neither'";
   - "the fifth statement (105 … 124)" does not count round 128's confirmation;
   - "39 killed" counts H12;
   - at draw 7 'b' and the cold restart are integrated over 12 periods and the other cells over 10 (difference ≤ 5.6·10⁻¹³).

## What I checked and how

All paths below are under W.

- **Reproduction:**
  - logs/fromraw.log, then `cmp` against the committed md: identical.
  - logs/pytest_full.log: "87 passed, 1 warning in 128.58s".
  - scripts/repro_c8a.py → logs/repro.log: jobs 0–6 each "0 differences" (9–80 s); recursive comparison of every field but `secs`.
- **Recount:**
  - scripts/dump.py → logs/dump.log;
  - scripts/check.py: the per-window patterns of 'warm_m', 'cold_m', 'b', 'warm_c', 'cold_c', `cen`, `cenc` and w1i − δ(p).
- **Code reading:**
  - 'cold_c' (`c8a_pT`, after the bisections);
  - 'w1i' (`c8a_job`, from the rerun's `r['w'][p·kk, idx[0]]`; w₀ from `law(m)`, the same as `c8a_pT`'s);
  - `cenc`, `ncs`, `rcm`, `arv`, `inc`, `mgt`, `lows` and the header.
- **`bracket_out`** (scripts/bracket_adv.py):
  - the exact Decimal of each double against the printed ends;
  - the 100,001 grid doubles; 40,000 `nextafter` neighbours (0 inward);
  - negatives (correct); NaN prints "NaN"; inf raises InvalidOperation;
  - the four data brackets (all outward of their exact doubles);
  - H12's output against the repository's on the test plants.
- **Factorial** (scripts/fact5.py → logs/fact5.json, logs/fact5.log):
  - It uses the repository's `seg_run`, `raised`, `fitted_state` and `c8a_cases`; the 2⁵ design, the grouping and the outcome rule (a leave within one period of R) are mine.
  - Corners against the raw cells: 'cold_m' (10000), 'cold_c' (01000) and 'warm_c' (10111) are bit-identical, and all-cold leaves at R everywhere.
  - 'b' is bit-identical except at draw 7, where it differs by ≤ 5.6·10⁻¹³: the raw ran it 12 periods.
- **Margins:** scripts/margin50.py reran 20 periods of 'cold_c' and 'b' at draw 1 from 41 and from 50.
- **Mutants** (tools/mutants.py from W/repo):
  - my R00–R22 (scripts/spec134.py) on the full suite → mut/res_r.jsonl, logs/mut_r.log;
  - the authors' H00–H43 on `-k "round128 or round124 or round102 or round105"` → mut/res_h.jsonl, logs/mut_h.log;
  - kill plants for my survivors in-process on mutated copies (scripts/killplants.py; `c8a_lines` only).
- **Literature:** six web searches, snippets only. One fetch, of arxiv.org, failed: the host is not allowed in this session's network.

## Verdicts per part

| claim (D60) | verdict | evidence |
|---|---|---|
| (iii) 157–162: 'warm_m' passes at every k at the 11; 'cold_m' leaves; continued starts pass; "a 2² factor separation … the outcome goes with the rest of the state" | verified (33 of 33 each). The design is complete and correctly read; its α₁ contrast is small at some k (D3ec-4) | [raw], [own] |
| 163–168: by kernel 1's centre at 10, kernel 3's at draw 2 from 58; centres +0.0377 … +0.0396, −0.0264 … −0.0257 | verified | [raw] |
| 168–171: 'cold_c' passes at every k at all 11, no other centre alone; "necessary and sufficient" | verified as a two-factor statement. My 2⁵ extends it to x, the other α and the unformed kernels with no interaction (D3ec-3). "Passes" includes a later leave (D3ec-5) | [raw], [own] |
| 171–176: draw 7's 'warm_m' leave at the alarm's period (k = 1–3, either δ); deficits 0.00115 … 0.00122; margins 0.00082 … 0.00089 (exact 0.000822 … 0.000893); arrivals 0.07627 / 0.07724 / 0.07797; 'cold_m' rule; 'cold_c' leaves at the alarm's period | verified (all 9 'cold_c' runs leave at period 10) | [raw] |
| 176–178: thresholds (0.0763511, 0.0763513] and 0.0763370 … 0.0763380, "rounded outward"; inside (0.0763, 0.0764]; above 0.07613 | verified on the data, and outward also against the exact doubles. "Rounded outward" is not guaranteed in general (D3ec-1) | [raw], [own] |
| 180–183: round 111's reading confirmed in 114; "the fifth statement (105, 111, 114, 119, 124)" | consistent with the log; round 128's confirmation not counted (D3ec-6) | [md] |
| Status (iii): "two 2² factor separations" | right; the full 2ⁿ was not run by the repository (D3ec-3) | [md], [own] |
| Numerics C8A | numbers match the md and the raw | [md], [raw] |
| md:681 | regenerated byte for byte; the header from `C8A_ROUNDS` (98 … 128) | [own] |

## Status of the round-128 issues

| D3dw id | status at 7d78bfa | what remains |
|---|---|---|
| D3dw-1 | **fixed on the data**: outward brackets, a raw test, 200 random brackets | the rounding is from the repr, not the double; the plants assert inward output; H12 is a false kill (D3ec-1) |
| D3dw-2 | **fixed**: `cl`, `arv`, the header split; plant Q1 | survivors R05/R06 in the 2² clause (D3ec-2) |
| D3dw-3 | **fixed**: 'cold_c', `cenc`, `ncs`, the md and D60; reproduced bit-identically | scope wording (D3ec-3) |
| D3dw-4 | **fixed**: the kill plants; H40–H43 killed [own] | — |
| D3dw-5 | **fixed** in D60, the code and the log | README round-124 clause (D3ec-6) |
| D3dw-6 | **fixed**: header, margins, no empty clause, Status, 'w1i' 0.07613 | the margins printed rounded to nearest (0.000893 as 0.00089), as elsewhere in the md, acceptable |

## Overall assessment

- **Is the claim about the isolated alarms settled? Yes.** At each of the 12 windows the attribution is now overdetermined:
  - two complete 2² separations computed by the repository;
  - my 2⁵ over every state component, with zero interaction by outcome and ≤ 0.001 in the largest distance;
  - round 128's one-at-a-time runs.
- **What could still change a conclusion is definitional, not numerical:**
  - *(a) δ(p) is a period mean.* Draw 7's alarm exists because δ(107) = 0.07667 lies above the cold threshold (0.0763511, 0.0763513], while the run's own instantaneous coefficient (0.07613) lies below it. With δ taken at pT the window would pass, assuming monotonicity, which round 128 scanned. That is the content of "the coefficient's", not a threat to it.
  - *(b) The cut 0.3 and the "past the alarm's period" rule.* At draw 1 from 50 the passes are 0.014 below the cut; the verdict survives because a later leave counts as a pass.
- **Nothing found this round would change a published number or a verdict.** Every surviving mutant leaves the md unchanged.
- **I recommend closing C8A.** Apply the minor and cosmetic items without new runs, apart from optionally adding the 2⁵.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ec-1 | minor | **`bracket_out` is outward only relative to `repr`** (`Decimal(repr(float(v)))`).<br>• [own, bracket_adv.py] For a double whose shortest repr lies on the 7-decimal grid, the exact value is below the repr in 50,000 of the 100,001 grid doubles in [0.07, 0.08] (lo printed inward) and above it in 49,999 (hi printed inward). Example: 0.0763377 is exactly 0.07633769999…, yet it prints as lo 0.0763377.<br>• Off the grid, never inward: 0 of 40,000.<br>• Equivalent on the data: all four brackets contain their exact doubles.<br>• The plants in test_round128 (warm [0.0763377, 0.0763379] → "(0.0763377, 0.0763379]") and round 124's "(0.0763300, 0.0763400]" assert inward output.<br>• H12 (`Decimal(float(v))`, which gives "(0.0763376, 0.0763380]" and "(0.0763299, 0.0763401]") is killed only by them ("1 failed, 3 passed"; md unchanged). The 200-random test compares via `float()` and cannot see the defect.<br>• The docstring's "the printed interval contains the stored one" is false in exact arithmetic. | **Code:** `Decimal(float(v))`.<br>**Tests:** update the two plants' expectations (or move the plants off the grid); compare as Decimal of the exact doubles; add grid-point inputs.<br>**Mutant record:** H12 is the fix, not a killed mutant ("39 killed" → 38 and one false kill). |
| D3ec-2 | minor | **Seven survivors of my R00–R22 on the full suite (87 passed each), none changing the md** [own, mut/res_r.jsonl].<br>• R04: `cenc` takes a pass as no leave only.<br>• R05: `rcm` continued starts with `any`.<br>• R06: `rcm` 'cold_m' with `any`.<br>• R07: `ncs` takes the intersection, not the subset.<br>• R08: margins printed unless all negative.<br>• R09: 'w1i' compared with the cold bracket's upper end.<br>• R10: 'w1i' compared with the own brackets' upper ends.<br>Kill plants (killplants.py):<br>• P05: a continued start reproducing at k = 2 — the mutant prints the 2² clause;<br>• P06: 'cold_m' passing at k = 2 — likewise;<br>• P07: two centres in 'warm_c', one in 'cold_c' — the mutant prints "necessary and sufficient";<br>• P08: an own arrival above its threshold — the mutant prints "−0.00010 … 0.00094 below their own thresholds";<br>• P09 and P10: 'w1i' inside a bracket — the mutant prints "below every threshold".<br>Equivalent by construction: R01, R11, R16, R18, R20, R21. Killed: R02, R03, R14, R15, R19 (md changes); R12, R13, R17, R22 (tests only). Control R00 survives.<br>Authors' 40: all 39 killed and H00 survives, under the -k subset [own]. | **Tests:** add P05–P10.<br>**Code:** decide R04 as a convention: either keep `ps_` and say so (D3ec-5), or require no leave.<br>**Mutant record:** the R list. |
| D3ec-3 | cosmetic | **"Necessary and sufficient" is computed against the rest of the state as one block, moving the formed centres one at a time** (D60:168–171; md; Numerics).<br>• [own, fact5] The full 2⁵ (α₁ / named centre / other formed centres / x / other α and unformed centres; 1152 runs) gives the outcome as a function of the named centre alone at 33 of 33 starts (α₁ alone at draw 7's 3).<br>• Stein–Alpert (Möbius) on the largest distance: the named factor +1.000 … +1.001 of the change; every other term \|·\| ≤ 0.001. | **D60 (iii):** "necessary and sufficient between the cold restart and each own-δ start's state at pT (the centre against the rest of the state); a full 2⁵ separation over the state's components shows no interaction (referee round 134)".<br>Or add the 2⁵ to `c8a_pT`. |
| D3ec-4 | cosmetic | **The α₁ contrast is not printed** [raw].<br>• At the 11 windows \|arrival − δ(p)\| = 1.2·10⁻⁶ … 1.7·10⁻³ (signed −1.4 … +1.7·10⁻³), against draw 7's decisive 1.15–1.22·10⁻³.<br>• At draw 2 from 58, k = 2, the two levels differ by 1.2·10⁻⁶.<br>• 'w1i' at the 11 windows is 0.00055 … 0.00182 above δ(p), stored but not printed. | **`c8a_lines` and md:** print the contrast span with the 2² clause, and 'w1i' − δ(p) at the cold restart's windows. |
| D3ec-5 | cosmetic | **A "pass" in 'cold_c', 'warm_m' and `rest` is `ps_`** (no leave, or a leave after R + 1).<br>• [own, margin50.py] At draw 1 from 50 'cold_c' and 'b' reach 0.2857 in the window's period 4, 0.014 below the cut — the excursion the δ(p) starts cross (0.3000, 0.3004).<br>• At draw 1 from 41 they reach 0.2148 (period 5), and 0.286 at period 13, outside the window. | **D60 (iii):** "passes the alarm's period", or give the margin at draw 1 from 50. |
| D3ec-6 | cosmetic | **History and record residues:**<br>• README D60 row, round-124 clause: "round 114's 'neither' restored" (D3dw-5 attributes it to round 111);<br>• D60:181 "the fifth statement (105, 111, 114, 119, 124)" omits round 128's confirmation;<br>• REVIEW_LOG and README "39 killed" includes H12 (D3ec-1);<br>• at draw 7 the cold restart and 'b' run nb = 12 periods while 'cold_m', 'warm_c' and 'cold_c' run 10, so the 2² cells are not on one output grid ('b' ≤ 5.6·10⁻¹³ from a 10-period run). | **README:** "round 111's (confirmed in 114)".<br>**D60:** "confirmed in round 128".<br>**Log:** H12's status.<br>**Code (optional):** take every cell's window from one nb-period run. |

## Sources

- **Stein–Alpert factor separation** [snippet]: Lunt et al., "Multi-variate factorisation of numerical simulations", GMD 14 (2021) 4307. All on/off combinations need 2ⁿ runs; the non-linear remainder is grouped as "synergy"; linear-sum and shared-interaction factorisations. https://gmd.copernicus.org/articles/14/4307/2021/ , https://eprints.whiterose.ac.uk/176043
  - Relevance: my 2⁵ is this design. The Möbius-inversion terms I computed are its interaction terms.
- **arXiv:2609.11948, "Climate-model factor separation with Shapley values and efficient sampling"** [snippet; status and date not checked]: reports proving Lunt et al.'s conjecture via Möbius inversion and Shapley values. https://arxiv.org/pdf/2609.11948
- **Heimersheim & Nanda, "How to use and interpret activation patching", arXiv:2404.15255** [snippet]:
  - denoising (clean into corrupt) tests sufficiency; noising (corrupt into clean) tests necessity;
  - the two are not complements and neither implies the other;
  - aliases: interchange intervention (Geiger et al. 2020), causal tracing (Meng et al. 2022).
  - https://arxiv.org/pdf/2404.15255 ; https://www.lesswrong.com/posts/FhryNAFknqKAdDcYy/how-to-use-and-interpret-activation-patching
  - Relevance: 'cold_c' is the denoising cell (the warm centre into the cold run) and 'warm_c' the noising cell. Their agreement is the strong case.
- **Subspace activation patching and "interpretability illusions"**, ICLR 2024 [snippet]: patching a whole state can act through dormant pathways. https://proceedings.iclr.cc/paper_files/paper/2024/file/70b8505ac79e3e131756f793cd80eb8d-Paper-Conference.pdf
  - Relevance: the reason to partition the "rest" (D3ec-3).
- **Actual causation** [snippet]:
  - Halpern–Pearl's AC2(a) necessity and AC2(b) sufficiency; NESS ("a necessary element of a sufficient set");
  - Beckers, "Causal Sufficiency and Actual Causation", arXiv:2102.02311.
  - https://arxiv.org/pdf/2102.02311
- **One-at-a-time designs cannot detect interactions** [snippet]: Saltelli & Annoni, Environ. Model. Softw. 25 (2010) 1508, https://hero.epa.gov/reference/1060594 ; Ferretti et al. review, https://arxiv.org/pdf/1711.11359
  - Relevance: 'warm_c' and 'cold_c' are one-at-a-time; the 2⁵ closes that gap.
- Not rechecked this round [memory or earlier rounds]: D60's other citations (Status, 259–281).

## Lesson

- **Round outward from the number itself, not from its printed name.** The shortest repr of a double is a decimal on the other side of the double about half the time it lands on the grid, and a test written with grid-valued plants will then pin the wrong behaviour and "kill" the correct fix.
- **When two one-at-a-time transplants agree in both directions, run the full 2ⁿ once.** With five groups it costs 32 runs per start, and it turns "necessary and sufficient against the rest" into "no interaction anywhere in the state", which is what settles an attribution.
