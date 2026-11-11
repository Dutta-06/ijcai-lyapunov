# Referee report D3dw, round 128: Proposition D‴.60 as revised after round 124

## Scope

- **Commit.** I reviewed b361d6d on claude/vigilant-clarke-hf939b.
  - Exported with `git -C … archive b361d6d | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref128. The repository's venv is symlinked as W/repo/venv.
  - HEAD moved to eb6e7d3 (round 127, D‴.61) while I worked. `git diff b361d6d HEAD` touches none of the reviewed files.
- **Read in full:**
  - review/drafts/D60.md (416 lines). The revision touches (iii) at 150–179 and the Numerics' C8A at 396–407.
  - experiments/exp_noise_resonance.py: the new `NB_C8A`, `seg_run` and `c8a_pT` (the 'pT', 'warm_m', 'cold_m', 'warm_c', 'ext' and 'thr' fields); `c8a_job`'s new call; and `c8a_lines`' pT part:
    - the classes `rest` (the cold restart's), `coef` (the coefficient's), `mixed` ("neither") and `nt` (not run);
    - `fol`/`folc` (the arrival rules), `dpn`, `coefs`, `cen`, `side`, `extd`, `gone`, `brk`, and `kname` by index.
  - What they call: `offsets_start`, `win_run`, `edge_run`, `c_run`, `raised`, `c_noise`, and `closed_loop.simulate`. Its state is (x, α, c), n = 1 + 2d, with no gate.
  - results/chaos/noise_resonance.md, the C8A paragraph (md:681).
  - tests/test_noise_resonance.py: `test_round124_c8a_pT_wiring`, `test_round124_c8a_lines_at_pT_and_the_labels`, `test_round124_on_the_raw`; the round-105 plants as rewritten; `test_round98_c8a_job_wiring`; `test_round102_…`, which reruns draw 7's job end to end.
  - tools/mutants/round124_D60_own.py; review/REVIEW_round124_D3ds.md; the Round 111, 114, 119 and 124 entries of review/REVIEW_LOG.md; and the D60 row of review/drafts/README.md.
  - `git diff f346010 b361d6d` for these files.
- **How I ran things:**
  - `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv;
  - NR_OUT=W/out holding a copy of the raw; TMPDIR=W/tmp; LCK_BUILD=W/build;
  - at most two heavy processes of mine at a time.
- **Nothing in the repository was modified.** I ran only read-only git there.
  - I stopped only my own duplicate mutant slot (PIDs 29636 and 1077), after checking their /proc cmdlines.
- **Tags:** [md] the committed md or draft; [raw] the committed raw; [own] my runs or code; [snippet] a web-search result only; [memory] not checked this round.

## Summary

**Design.** The revision answers D3ds-1 … D3ds-7.
- **`c8a_pT`, called by `c8a_job` at every window with p ≥ 2.** For each earlier start (k = 1, 2, 3; own δ(p − k) or δ(p)) it:
  - runs the start to pT (`seg_run`) and keeps its state there: the arrival |w₁| − w₀ and the formed centres;
  - continues it ('b');
  - reruns it with kernel 1's α set to the cold restart's ('warm_m', the matched arrival);
  - runs the cold restart with the start's α₁ ('cold_m');
  - for the own-δ starts, runs the state with one formed centre at its feature ('warm_c');
  - runs to p + R + 2 where that passes the window's end ('ext');
  - bisects kernel 1's threshold where bracketed ('thr', 12 bisections).
- **`c8a_lines`** classes each `cold` window by 'warm_m' over the own-δ k: the cold restart's, the coefficient's, "neither", or not run. It prints the arrivals, the cold_m rule, the thresholds, the centre that explains each window and the periods past the window's end.
- **D60 (iii):** 11 of the 12 alarms are the cold restart's, each by a centre (kernel 1's at 10, kernel 3's at draw 2 from 58). Draw 7's is the coefficient's, not the restart's. Round 114's "neither" is restored.

**Verdicts**
- **Reproduction is exact** [own]:
  - `--from-raw` gives the md byte for byte;
  - the 83 tests pass;
  - all 7 C8A jobs, rerun with the script's own `c8a_cases`/`c8a_job`, equal raw['C8A'] in every field but `secs`, the 7 new fields included;
  - the raw differs from f346010's only by those fields; every older field is bit-identical.
- **'warm_m' is the right transplant** [own, snippet; checked by reading the code and by independent runs]:
  - only kernel 1's α is changed, to exactly the cold restart's;
  - the other 20 α, the 21 centres and the plant's x come from the start's own state at pT, which is the integrator's whole state;
  - the noise from p is the same `c_noise(…, shift=p)` call as the cold restart's, and phase-continuous with the start's first segment;
  - the four runs cold, 'cold_m', 'warm_m' and 'b' form a complete 2² factor separation (Stein–Alpert), with factors "kernel 1's coefficient" and "the rest of the state".
- **"The coefficient's" (draw 7) and "the cold restart's" (the 11) follow** [own]. My forward and reverse one-at-a-time transplants at 5 windows (450 runs) flip the outcome on exactly one coordinate at every start: α₁ at draw 7, kernel 3's centre at draw 2 from 58, kernel 1's centre at draws 1 from 50, 3 from 166 and 2 from 6.
- **Every new number of D60 and the md agrees with the raw** [raw]. One exception: a printed threshold bracket that excludes its threshold (D3dw-1).
- **No major issue.**

**Findings**
1. **D3dw-1 (minor). The printed threshold brackets are rounded to nearest, and two of them exclude the threshold** [own].
   - Cold restart: printed "(0.0763511, 0.0763512]" (md and D60:169). At δ = 0.0763512 the cold restart passes. So the threshold lies in (0.0763512, 0.07635121], outside the printed interval.
   - Own-δ k = 2: printed "(0.0763371, 0.0763373]". At 0.0763373 it passes.
2. **D3dw-2 (minor, latent). The matched-arrival test is applied only to windows already in `cold`** [own, plant].
   - An own-δ start that reproduces the alarm by arriving above δ(p) is the mirror of draw 7. In such a window 'warm_m' (δ(p) set at pT) passes, yet the window is printed as reproduced and gets no pT verdict.
   - Equivalent on the data: 0 of 12 windows are reproduced.
3. **D3dw-3 (minor). "By a centre" is computed in one direction only** ('warm_c': sufficiency on the warm base) [raw, own].
   - Necessity rests on the scan at 9 of the 10 kernel-1 windows. It is untested at draw 1 from 50 (the warm centres are inside the bracket) and at draw 2 from 58 (kernel 3; the scan moves only kernel 1's centre).
   - My reverse transplants close both gaps: both windows pass with the warm centre.
4. **D3dw-4 (minor). Four of my 26 mutants survive the full suite:** Z16, Z19, Z22, Z26.
   - All four are equivalent on the data, none changes the md, and each is killed by a plant.
   - The authors' 59 give 57 killed; G57 and G58 survive, equivalent by construction as recorded.
5. **D3dw-5 (cosmetic). The classification's vocabulary and history.**
   - "Neither" means the coefficient's in D60 and the log, but the code's "neither" class means mixed outcomes over k.
   - "Round 114's reading" was round 111's (D3df-1, major), which round 114 confirmed.
   - The log's "the third verdict on that window" is the fifth statement: 105 cold restart's, 111 neither, 114 neither, 119 cold restart's, 124 the coefficient's.
   - Round 119's lesson still cites "(here round 114's 'neither')" as an example, with no forward note.
6. **D3dw-6 (cosmetic). Wording:**
   - the md paragraph's header lists rounds "98, 102, 105, 111" only;
   - the arrival margin that decides is never printed: the own-δ starts arrive 0.00082–0.00089 below their own thresholds; D60 gives 0.0012 below δ(107);
   - a malformed "at 0 by a centre — ," appears when no centre explains a window;
   - the Status does not mention the pT transplants;
   - the run's coefficient at 107T (0.0761) is still not stored, though `c_run` reproduces it.

## What I checked and how

All paths are under W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref128.

- **Reproduction:**
  - W/logs/fromraw.log, then `cmp` against the committed md: identical (md5 585d5bef…).
  - W/logs/pytest_full.log: 83 passed, 1 old warning, 124.45 s.
  - W/scripts/repro_c8a.py → W/logs/repro_0.log and W/logs/repro_1.log: 7 of 7 jobs, 0 differences, 13.3–90.9 s each.
  - Raw against f346010's (`git show f346010:…raw.npz`, pickled per key): only C8A changed, by adding `dp, pT, warm_m, cold_m, warm_c, thr` (and `ext` at draw 7) to each of the 12 rows.
- **Recount:** W/scripts/dump_c8a.py → W/logs/dump_c8a.log. I recounted every new number of D60:150–179, the Numerics 396–407 and md:681 [raw]:
  - draw 7's arrivals 0.07544/0.07551/0.07552 and 0.07627/0.07724/0.07797; deficits 0.00115–0.00122;
  - draw 7's 'warm_m' all [10, 0.31036…]; 'cold_m' −1 below the threshold, 10 above it;
  - 'ext': cold 0.8218, 0.8341; own-δ 0.0073 … 0.0299;
  - the threshold brackets;
  - at the 11 windows at k_e = 20: 'warm_m' −1 at all 66 (both δ kinds); 'cold_m' leave at all 66;
  - 'warm_c': kernel 1's centre at 10 windows, kernel 3's at draw 2 from 58; centres +0.0377 … +0.0396 and −0.0264 … −0.0257;
  - draw 1 from 50: 0.30001 and 0.30036 (margins 6.1·10⁻⁶ and 3.6·10⁻⁴); arrivals 0.09533 and 0.09647 against 0.09138; 'warm_m' −1.
- **Own transplants.** W/scripts/probe.py → W/logs/probe_{0..4}.json, probe_a.log and probe_b.log. It uses the repository's `simulate`, `raised`, `c_noise` and `fitted_state`. The continuous run, the capture at pT, the transplants and the outcome rule are mine.
  - **Windows:** draw 7 from 107 (own δ and δ(p), k = 1–3), and at k_e = 20 draws 2 from 58, 1 from 50, 3 from 166 and 2 from 6 (own δ, k = 1–3).
  - **Capture at pT:** each start is run as one continuous integration from p − k and its state taken at output kT. Against the two-segment raw: arrivals agree to ≤ 1.3·10⁻¹³, centres to ≤ 3.7·10⁻¹², and 'b'/'warm_m' to ≤ 2·10⁻¹² in the largest distance. My cold restart equals C8D's window exactly.
  - **Per start, 25 transplants:**
    - on the matched base ('M'), then on the warm base ('W'), one more coordinate set to the cold value at a time: c₁, c₂, c₃, w₂, w₃, the 18 unformed kernels, x;
    - the reverse ('C'): the cold restart with one coordinate from the warm state at a time (w₁, c₁, c₂, c₃, w₂, w₃, unformed, x), with all warm centres, and with everything but w₁.
- **Thresholds.** W/scripts/probe2.py → W/logs/probe2.json:
  - 21-point δ scans (0.0755–0.0775) on the cold and the own-δ k = 2 states at 107T;
  - the printed and stored bracket ends;
  - `c_run`'s own coefficient at 107T.
- **Plants:**
  - W/scripts/plants.py: Q1 (an own-δ reproduction above δ(p) whose 'warm_m' passes) and Q2 (mixed outcomes over k);
  - W/scripts/killcheck.py: one plant per surviving mutant, run on private mutated copies in W/kc.
- **Mutants.** W/scripts/spec128.py holds my Z00–Z26; W/scripts/allspec.py combines them with the authors' G00–G59.
  - On `-k "round124 or round89_md or round98 or round105 or round102 or win_run or edge_run"` (11 tests) → W/mut/k1/res_a.jsonl and W/mut/k1/res_b.jsonl.
  - Survivors on the full suite → W/mut/full/res_c.jsonl.
- **Literature:** six web searches, snippets only; no page was fetched.

## Verdicts per part

| claim (D60.md lines) | verdict | evidence |
|---|---|---|
| (i), (ii); (iii) outside 150–179; Proof | unchanged since f346010; not reverified | [md] |
| md:681 C8A paragraph | regenerated byte for byte; every new number recounted; the header's round list is stale (D3dw-6) | [raw], [own] |
| 150–151 "4 of the 33 … inside the scan's bracket, untested there (draw 7's 3 have no bracket)" | verified (the round-124 test asserts it on the raw). At draw 1 from 50 my reverse transplant now tests the warm centre inside that bracket: it passes | [raw], [own] |
| 152–157 own-δ starts = C8D's windows at p − k (36 leaves; 28 / 2 / 6) | verified (test on the raw); "adds only draw 7" is right (R + k > NP_C7 only there) | [raw] |
| 157–159 "They do not isolate the restart … matched arrival … pass at every k at the 11" | verified: 66 of 66 'warm_m' pass at the 11 (both δ kinds); my runs reproduce them to 10⁻¹² | [raw], [own] |
| 159–165 "11 of the 12 … each by a centre — kernel 1's at 10 …, kernel 3's at draw 2 from 58" | verified as sufficiency ('warm_c'). Necessity is computed only through the scan (D3dw-3). My reverse transplants confirm necessity at 4 of the 11 windows | [raw], [own] |
| 165–169 draw 7 "the coefficient's, not the restart's", the arrivals, 'cold_m' | verified. My factor separation: α₁ alone flips the outcome, in every direction and at every k; nothing else does | [raw], [own] |
| 169–171 thresholds "(0.0763511, 0.0763512]" and "between 0.0763371 and 0.0763380" | the first is false by ≤ 1.5·10⁻⁹: the threshold is above 0.0763512 (D3dw-1). The second holds (all three warm thresholds are in (0.0763372, 0.0763380]). The bisection's monotonicity holds on my 21-point scans | [own] |
| 171 "above the run's own coefficient at 107T (0.0761 …, not reproduced by the repository's scripts)" | reproduced with the repository's `c_run`: 0.076126; δ(107) = 0.0766659 = C3's period mean | [own] |
| 172–173 "0.3104 … 0.8218 one period later — not a marginal crossing —, own-δ starts below the cut" | verified ('ext') | [raw] |
| 173–176 "round 114's reading, 'neither', stands …" | substance right. The attribution is round 111's (confirmed in 114), and the term clashes with the code's class (D3dw-5) | [md], log |
| 176–179 draw 1 from 50 (margins, arrivals, "gone with δ(50) set at pT") | verified | [raw] |
| Numerics 396–407 | numbers match the md; inherits D3dw-1 | [md], [raw] |
| Status (iii) | consistent; silent on the pT transplants (D3dw-6) | [md] |
| Classification across D60, Numerics, README, log, commit | consistent (11 cold restart's / draw 7 the coefficient's); vocabulary and history as D3dw-5 | [md] |

## Status of the round-124 issues

| D3ds id | status at b361d6d | what remains |
|---|---|---|
| D3ds-1 (minor) | **implemented as proposed**: `c8a_pT` ('warm_m', 'cold_m', 'thr'), the classes, the `dpo` and drift sentences computed from the arrivals, D60/Numerics/README/log updated; reproduced bit-identically and confirmed by my transplants | the brackets rounded inward (D3dw-1); the matched test gated on `cold` (D3dw-2); round 119's lesson not amended (D3dw-5) |
| D3ds-2 (minor) | **fixed**: D60:152–157, the docstring, the raw test (`n == [28, 2, 6]`) | — |
| D3ds-3 (minor) | **fixed**: the clause printed from `nco`; asserted on `report(raw)` | — |
| D3ds-4 (minor) | **fixed**: 'warm_c', printed by kernel and side; D60 names kernel 3 | the reverse direction not computed (D3dw-3) |
| D3ds-5 (minor) | **fixed**: 'ext'; "not a marginal crossing" | — |
| D3ds-6 (minor) | **fixed**: P3, P4, P6 kill G54, G53, G55 [own: killed]; G57, G58 equivalent by construction (agreed) | new survivors Z16, Z19, Z22, Z26 (D3dw-4) |
| D3ds-7 (cosmetic) | **mostly fixed**: "k periods longer", `kname` by index, "of the 33 in a bracket", the margins, the 4 inside in (iii), the docstrings, the round-105 aside dropped, the attribution per window | the md header's rounds; D3dn-5/6 deferrals (`zd`'s "nan", the row without 'warm', p = 0's wrap, the 'wc'-base identities, the delayed band) still open |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dw-1 | minor | **Brackets rounded to nearest exclude the threshold** (`coefs`' `:.7f`; md:681; D60:169–170, 404).<br>• [raw] cold [0.07635111803, 0.07635121487] printed (0.0763511, 0.0763512]; own-δ k = 2 [0.07633706, 0.07633735] printed (0.0763371, 0.0763373].<br>• [own, probe2] the cold restart at δ = 0.0763512 passes and at 0.07635121487 leaves; k = 2 at 0.0763373 passes. Both thresholds lie above the printed upper end.<br>• k = 1's printed lower end (0.0763377) is above its true one (0.07633767); k = 3's printed upper end (0.0763375) is below its true one (0.07633751). Not probed.<br>• Round 124's log called the referee's own printed bracket "a truncation"; the repository now prints the same defect. | **`c8a_lines`/`coefs`:** print lo rounded down and hi rounded up at 7 digits — (0.0763511, 0.0763513]; (0.0763376, 0.0763380], (0.0763370, 0.0763374], (0.0763372, 0.0763376] — or print 9 digits.<br>**D60:169 and 404** accordingly.<br>**Test** on the raw: the printed interval contains the stored one. |
| D3dw-2 | minor | **The matched-arrival verdict is applied only to `cold`** (`c8a_lines`: `rest`/`coef`/`mixed` ⊂ `cold`, which excludes `warm`).<br>• [own, plant Q1] An own-δ start that reproduces the alarm by arriving at 0.075 above δ(p) = 0.07 has 'warm_m' passing and 'cold_m' leaving. It is printed "leaves within one period of the alarm's at 1 of 2 (…)", with no pT verdict; the pT sentence covers "the 1 window where no own-δ start leaves near…".<br>• This is draw 7's confound mirrored: round 124's lesson applied to passes, not to reproductions. Equivalent on the data (`warm` is empty). | **`c8a_lines`:** class every window that has 'warm_m'. Print a reproduction whose matched arrival passes as "the arrival's, not the alarm".<br>**Test:** plant Q1 (W/scripts/plants.py). |
| D3dw-3 | minor | **"By a centre" is one-directional** (D60:159–164; md "by a centre — kernel 1's centre at 10 …").<br>• [raw] 'warm_c' shows sufficiency on the warm base. Necessity comes only through the scan (9 of 10), not at draw 1 from 50 (inside the bracket) and not at draw 2 from 58 (kernel 3, unscanned).<br>• [own, probe_{0..4}.json] At 5 windows × every k (450 runs) exactly one coordinate flips each outcome, forward on the matched base, forward on the warm base, and reverse on the cold base: α₁ (draw 7), c₃ (2 from 58), c₁ (1 from 50, 3 from 166, 2 from 6).<br>• The cold restart with the warm kernel-1 centre passes at draw 1 from 50 (k = 1, 2, 3); with the warm kernel-3 centre it passes at 2 from 58. c₂, w₂, w₃, the unformed kernels and x never matter.<br>• cold / 'cold_m' / 'warm_m' / 'b' is a 2² factor separation with no synergy by outcome. | **`c8a_pT`:** 'cold_c' (the cold restart with each own-δ start's formed centre, one at a time; 9 runs per window).<br>**`c8a_lines`:** print "necessary and sufficient" where both directions hold.<br>**D60 (iii)/Status:** say so; cite the 2² design (Stein–Alpert, via Lunt et al. 2021) [snippet], or quote these runs as the referee's (round 128). |
| D3dw-4 | minor | **Four survivors of the full suite (83 passed), all equivalent on the data, none changing the md** [own, W/mut/full/res_c.jsonl]:<br>• Z16: `rel` "above" with `any`;<br>• Z19: `late` with `>=`;<br>• Z22: `cen` with any leave;<br>• Z26: `folc` over the own-δ `cold_m` only.<br>Each is killed by one plant (W/scripts/killcheck.py):<br>• Z16: two δ(p) reproductions at 0.071 and 0.069, printed "against" vs "above";<br>• Z19: a δ(p) k = 3 start at R + 4, printed as late;<br>• Z22: 'warm_c' leaving at R + 5;<br>• Z26: an inconsistent δ(p) 'cold_m'.<br>The authors' 59 give 57 killed (W/mut/k1/res_a.jsonl). `seg_run`/`c8a_pT` numerics mutants are killed by `test_round102`'s end-to-end rerun of draw 7. | **Tests:** add the four plants.<br>**Mutant record:** Z16, Z19, Z22, Z26. |
| D3dw-5 | cosmetic | **Vocabulary and history of the classification** (D60:173–176; REVIEW_LOG Round 124 lesson; round-119 entry; `pts`).<br>• "Neither" is round 114's (neither the cold restart's nor the centre's: the coefficient's), but the code's "neither" is `mixed` over k (plant Q2 prints "neither at 1").<br>• The reading originated in round 111 (D3df-1, major), after round 105's applied text had called all 12 the cold restart's. The log's "the third verdict (114: …; 119: …; 124: …)" omits both.<br>• Round 119's lesson still offers "(here round 114's 'neither')" as an example of a borrowed state's effect. | **Code:** rename the `mixed` label ("mixed over k").<br>**D60:** "round 111's reading (confirmed in round 114), 'neither' — neither the cold restart's nor the centre's: the coefficient's".<br>**Log:** "the fifth statement on that window (105, 111, 114, 119, 124)"; a forward note on round 119's lesson. |
| D3dw-6 | cosmetic | **Wording and robustness:**<br>• the md header "(C8A, referee rounds 98, 102, 105, 111)";<br>• the arrival rule printed against δ(p): the own-δ starts pass 0.00082–0.00089 below their own thresholds, D60 says "0.00115 … 0.00122 below δ(107)";<br>• "leaves … at 0 by a centre — , by no single centre at 1" when no centre explains (Z22 plant);<br>• Status (iii) is silent on the pT transplants;<br>• the run's coefficient at 107T is not stored (mine: 0.076126);<br>• the D3dn-5/6 deferrals;<br>• HANDOFF §5 still lists round 124 "being applied". | **Code:** the header from the rounds list; print each arrival's margin to its threshold where 'thr' exists; drop the empty "by a centre" clause; store Wi[0] as 'w1i'.<br>**D60 Status:** "and at pT, the matched arrival, its reverse and one centre at a time".<br>**HANDOFF.** |

## Sources

- **Stein & Alpert (1993) factor separation** [snippets]:
  - Lunt et al., "Multi-variate factorisation of numerical simulations", GMD 14 (2021) 4307: the Stein–Alpert version needs all 2ⁿ runs, is unique and complete, and groups the non-linear remainder into an interaction ("synergy") term; the linear one-at-a-time alternative needs fewer runs. https://gmd.copernicus.org/articles/14/4307/2021/ and https://eprints.whiterose.ac.uk/176043/1/gmd-14-4307-2021.pdf
  - Krichak & Alpert, JAS 2002: "four model simulations are sufficient to evaluate the role of two chosen factors and that of their synergistic interaction". https://www.tau.ac.il/~pinhas/papers/2002/Krichak_Alpert_JAS_2002.pdf
  - Cambridge, *Factor Separation in the Atmosphere*, preface: 2ⁿ simulations for an n-factor problem. https://www.cambridge.org/core/books/abs/factor-separation-in-the-atmosphere/preface/EEAD86534A24B10032A296B8CDDDDFE1
  - Relevance: cold / 'cold_m' / 'warm_m' / 'b' is this 2² design (D3dw-3).
- **Saltelli & Annoni, "How to avoid a perfunctory sensitivity analysis", Environ. Model. Softw. 25 (2010) 1508–1517** [snippet]: one-at-a-time designs cannot detect interactions. https://hero.epa.gov/reference/1060594 and https://publications.jrc.ec.europa.eu/repository/handle/JRC60144
  - Relevance: 'warm_c' and my single-coordinate transplants are one-at-a-time. The full exchanges ('warm_m', the cold restart with all warm centres) agree with them, so no synergy shows at the outcome level.
- **Hoffman & Kalnay (1983), lagged average forecasting** [snippets]: members are forecasts from earlier analyses verifying at the same time; older members are less accurate and get less weight.
  - Relevance: the earlier starts are lagged members, and their age-dependent arrival is round 124's confound.
  - https://ntrs.nasa.gov/citations/19840014007 ; Kalnay's ECMWF history, https://www.ecmwf.int/sites/default/files/elibrary/2017/17604-historical-perspective-earlier-ensembles-and-forecasting-forecast-skill-and-global-leading.pdf
- **Mulholland, Laloyaux, Haines & Balmaseda, Mon. Wea. Rev. 143 (2015) 4631–4644** [snippet]: initialisation shock; the supplement's paired hindcasts differ only in the ocean's initial conditions (one component transplanted). https://centaur.reading.ac.uk/40638 and https://reading-clone.eprints-hosting.org/40638/2/mulhollandetal2015_supp_revise2.pdf
- Not rechecked this round [memory or earlier rounds]: D60's own citations (Status, 231–281).

## Lesson

- **Compute both directions of a transplant.** The four runs cold, cold with the warm coefficient, warm with the cold coefficient, and warm are a complete two-factor separation and settle "the coefficient or the restart" by outcome. Attributing the rest to one centre needs the reverse cell too (the cold restart with that warm centre); otherwise "by a centre" is sufficiency only.
- **Apply a new control to every class it could reclassify.** That includes the reproductions it was not written for, not only the passes that prompted it.
- **Round a printed bracket outward.** Printed to nearest at a fixed number of digits, a bracket can exclude the very quantity it brackets, as the cold threshold here does by 1.5·10⁻⁹.
