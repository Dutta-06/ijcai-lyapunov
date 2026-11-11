# Referee report D3eg, round 138: Proposition D‴.60 as a whole, after round 134

## Scope

- **Commit.** e53ec0d on claude/vigilant-clarke-hf939b. The D60 files last changed in 10e4244 (round 134 applied). I checked that `git diff 085250b 10e4244` covers exp_noise_resonance.py, test_noise_resonance.py, noise_resonance.md and its raw, D60.md, REVIEW_LOG.md and tools/mutants/round134_D60_own.py.
- **Export.** `git archive e53ec0d | tar -x -C W/repo`, with W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138 and the venv symlinked.
- **How I ran things.**
  - Environment: `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, NR_OUT=W/out.
  - One heavy process at a time. No process was killed, and the repository was not modified.
  - Another session's `pytest tests/test_moving_width.py` was running alongside; it is not mine.
- **Read in full:**
  - review/drafts/D60.md (all 440 lines);
  - review/REVIEW_round134_D3ec.md and the Round 134 entry of REVIEW_LOG.md;
  - the diff;
  - the D60 row of review/drafts/README.md;
  - PAPERS.md (every entry's "Rests on" and "Missing").
- **Code read:** `seg_run`, `c8a_pT`, `c8a_job`, `c8a_cases`, `c8f_job`, `c8f_cases`, `c8f_lines`, `bracket_out`, `ctr2`, `psm`, `cenc`, `left_periods`, `raised`, `c_noise`, `c_run`, `main` and `report`.
- **Not re-read line by line:** results/chaos/noise_resonance.md beyond C8A, C8F and spot checks.
- **Tags:** [md] the committed md or draft; [raw] the committed raw; [own] my runs or code; [snippet] a web-search result only, no page fetched; [memory] not checked this round.

## Summary

**Design.** D‴.60 has three parts, not four: there is no (iv).
- (i): the linear response, an oscillator with stiffness K(t), its sensitivity |S| and Bode's integral.
- (ii): the second order, the even and odd parts and the period means as Floquet sums.
- (iii): ejection from the fit's basin, numerical. It covers the noise-free edge, its scaling, the windowed ejecting sets (C7, C7L, C7S), the time-resolved replay and its nulls (C8, C8N, C8R, C8D, C8K), the attribution of the 12 isolated alarms (C8A plus round 134's C8F, the 2⁵ factor separation), C9 and C9A, and E per tone.

Round 134 added C8F: 36 starts × 32 corners over five groups of the state, with Möbius (Stein–Alpert) terms on the window's largest distance.

**Verdicts**
- **Reproduction is exact** [own]:
  - `--from-raw` gives the committed md byte for byte (md5 3097f825…);
  - the 91 tests pass (152 s);
  - 12 of the 36 C8F starts, recomputed with the script's own `c8f_cases`/`c8f_job`, are bit-identical to raw['C8F'] (every corner's leave period and largest distance). These are draw 7 from 107, draw 2 from 58, draw 1 from 50 at s = 0.08, and draw 3 from 210, each at k = 1, 2, 3.
- **C8F's design is right and its reading of the outcome holds** [own, raw]:
  - the five groups do partition the state (x, all 21 α, all 21 c);
  - the Möbius inversion is correct (recomputed independently: shares 0.99978 … 1.00123, other terms ≤ 1.22·10⁻³; from the 11111 baseline 0.99982 … 1.00143 and ≤ 1.43·10⁻³);
  - the 'cold_m', 'cold_c' and 'warm_c' corners are bit-identical to C8A at all 36 starts;
  - 'b' is bit-identical except at draw 7's 3 starts (5.6·10⁻¹³, the 12-period run).
- **What the "no interaction" claim adds is overstated** (D3eg-1). The Möbius result on the largest distance is the binary outcome restated, because that distance is bimodal. On a non-saturating response there are interactions.
- **The proposition's labels are mostly right.** (i) and (ii) are proved or derived; (iii) is numerical and says so; the scaling is "an interpolation, not a law".
- **No major issue.** Nothing found changes a published number or a verdict.

**Findings**

| id | severity | finding |
|---|---|---|
| D3eg-1 | minor | "Shows no interaction" holds for the outcome only. The Möbius terms on the saturating largest distance restate it; on the pre-alarm distance the named centre interacts with group C (+0.28–0.29 of the change at draw 2 from 58). |
| D3eg-2 | minor, latent | Six single-point mutants of `c8f_job`'s start segment survive the tests (J01–J05, J09). They are caught only after a C8F rerun, by the raw-identity test. |
| D3eg-3 | cosmetic | Two `c8f_lines` mutants (J06, J07) are equivalent on the data. |
| D3eg-4 | minor (journal) | (iii) is one 4,561-word paragraph with 26 "referee round" citations. About 1,500 of those words attribute 12 false alarms of a diagnostic. |
| D3eg-5 | minor (journal) | No PAPERS.md entry rests on D‴.60. |
| D3eg-6 | cosmetic | Stale script docstring and usage line; C8A_ROUNDS not updated for round 134. |
| D3eg-7 | minor | The headline mechanism is a coefficient-only section, yet C8A/C8F show the centre decides at 11 of the 12 isolated alarms. |
| D3eg-8 | cosmetic | (ii)'s averaged-loop identity is asserted in the Proof but only checked numerically; the Floquet bound assumes M is diagonalisable. |

## What I checked and how

All paths below are under W.

- **Reproduction:**
  - logs/fromraw.log, then `cmp` against the committed md: identical.
  - logs/pytest_full.log: "91 passed, 1 warning in 152.47s".
  - scripts/repro_c8f.py → logs/repro_c8f.log: "starts 12 differing 0" (21–48 s per window).
- **C8F on the raw** (scripts/check_c8f.py → logs/check_c8f.log):
  - Möbius terms from both baselines.
  - Corner-against-C8A differences per cell: only 'b' at (k_e 10, s 0.02, draw 7, p 107) is nonzero, max 5.611·10⁻¹³.
  - Largest distance per class: the leaving corners' spread is ≤ 1.86·10⁻³, the passing corners' ≤ 1.26·10⁻³, and the smallest gap between the classes is 0.280.
  - scripts/dump_c8f.py → logs/dump_c8f.log:
    - leaving corners 0.3104 (draw 7), 0.41–1.90 (k_e = 20);
    - passing corners 0.030–0.286.
- **A 2⁵ on a non-saturating response** (scripts/c8f_cont.py → logs/c8f_cont.log). Same corners as `c8f_job`; the response is kernel 1's period-averaged distance in each period before the alarm's.
  - Draw 2 from 58, period 1:
    - named centre (kernel 3) +1.749 … +1.796 of the change;
    - group C (kernel 1's and kernel 2's centres) −1.030 … −1.088;
    - interaction 01100 +0.283 … +0.292.
  - Draw 7 from 107, period 3: α₁ +1.198 … +1.238 and the named centre −0.198 … −0.223.
  - Draw 3 from 210 and draw 1 from 50 at s = 0.08: the named centre +1.000, other terms ≤ 0.002.
- **Mutants** (tools/mutants.py; `-k "round134 or round128 or round124 or round89"`, which contains every test that references `c8f_*`, C8F or KEYS, plus the md regeneration):
  - mine, scripts/spec138.py (J00–J10) → mut/res_k.jsonl;
  - the authors' K00–K28 → logs/mut_h.log: all 28 killed, control K00 survives;
  - reruns of C8F under J01 and J09 (copies mut_J01, mut_J09) → logs/mut_rerun.log.
- **Spot checks of D60 against the md** [md]: √K_i 15.76/10.51/13.13; |S(√K)| = 1.866 at m = 12; |S| = 1.41/1.38 at the bound (by hand); 0.99400; 1.495; α = 0.537; F(2, 10) = 8.9; 33.56; 0.0808; ×4.08; 0.706/0.975; 14805; C8R's 33/37, 33/37, 29/33; C8A and C8F verbatim.
- **Literature:** seven web searches, snippets only, no pages fetched.

## Verdicts per part

| claim (D60) | label in D60 | verdict | evidence |
|---|---|---|---|
| (i) the oscillator y″ + k_eẏ + Ky = K(ṅ − f̂′n), √K_i = (Gm)^{1/2}\|w_i\|/(σe^{1/2}) | derived; quasi-static | correct; the quasi-static law has no error bound, only measured (0.999–1.082) | [md], [own] arithmetic: 15.76, 10.51, 13.13 |
| (i) \|S\|, \|S\|²−1 = K(2ω²−K)/…, Bode's integral, \|S(√K)\| = (1+K/k_e²)^{1/2}, the bound's \|S\| | proved (frozen K) | correct. L has a pole at 0 on the axis (allowed: ln\|S\| ~ ln ω integrable); frozen-time object, as stated | [own] by hand; [md] quadrature 3·10⁻¹² |
| (i) the band's factor 1.495 → 1.225 | measured | consistent with the md | [md] |
| (ii) even/odd split, regular perturbation on a finite horizon | proved | correct | [md] |
| (ii) Floquet sums, bound, alternation, gτ/T offset | derived and checked (F) | correct for a diagonalisable M (real multipliers printed, distinctness not stated: D3eg-8) | [md] |
| (ii) "the constant part's response is the averaged loop's … exact along the orbit" | Proof asserts it; Status cites B5 numerically | the label is mixed (D3eg-8) | [md] |
| (iii) noise-free edge, anatomy, interpolation α, β, "not a law" | numerical | correctly hedged | [md] |
| (iii) "a formed kernel leaves when the drift has carried its coefficient into the raises the coming noise ejects" | numerical | supported as a replay (stated). The set is a section at centres-at-features; at 11 of the 12 isolated alarms the run's centre decides (D3eg-7) | [md], [raw] |
| (iii) C7/C7L/C7S, C8/C8N/C8R/C8D/C8K, C9/C9A | numerical | reproduced via md bytes; not re-derived this round | [md] |
| (iii) C8A: 11 of 12 alarms the cold restart's by a centre (necessary and sufficient); draw 7's the coefficient's | numerical | verified (round 134 and the C8F sample) | [raw], [own] |
| (iii) C8F: the outcome is one group's alone at 36 of 36 starts | numerical | verified, 12 starts recomputed bit-identically | [own] |
| (iii) C8F: "shows no interaction"; deciding term +1.000 … +1.001, others ≤ 1.2·10⁻³ | numerical | numbers right; the interpretation overreaches (D3eg-1) | [own] |
| (iii) "the corners that are C8A's runs give C8A's leave periods at all 36 (within 5.6·10⁻¹³)" | numerical | verified; bit-identical except 'b' at draw 7 | [raw] |
| Status, Numerics (C8F line) | — | consistent with the md and raw | [md] |
| Status, novelty "in the instance of the paradigm's kernel law" | literature | no counter-example found (D3eg sources); appropriately scoped | [snippet] |

## Status of the round-134 issues

| D3ec id | status at e53ec0d | remaining |
|---|---|---|
| D3ec-1 (`bracket_out` from the repr) | **fixed**: `Decimal(float(v))`; the plants corrected; a grid-double test (inw > 50); my J10 (lo ROUND_DOWN) killed | — |
| D3ec-2 (seven printing survivors) | **fixed**: P05–P10 as tests; R04 kept as a convention and pinned; the authors' K-list all killed [own] | new survivors in `c8f_job` (D3eg-2) |
| D3ec-3 (scope of "necessary and sufficient") | **fixed**: C8F added, the wording "(the centre against the rest of the state)", bit-identical to round 134's 2⁵ on my sample | the reading of the Möbius terms (D3eg-1) |
| D3ec-4 (α₁ contrast) | **fixed**: printed "1.2e-06 … 1.7e-03 apart", 'w1i' "0.00055 … 0.00182 above δ(p)" | — |
| D3ec-5 ("pass") | **fixed**: "passes the alarm's period"; margin 0.0143 printed | — |
| D3ec-6 (history) | **fixed** in D60, the README and the log | C8A_ROUNDS still ends at 128 although round 134 added clauses to that paragraph (D3eg-6) |

## Journal readiness

- **Where D‴.60 could be published.** No PAPERS.md entry rests on D‴.60.
  - Paper 12 (measurement noise) rests on D‴.26–D‴.33.
  - Paper 17 (rich transit) lists D‴.49, 50, 52, 54, 56, 57 and 58 under "Missing", but not D‴.60.
  - Only HANDOFF.md:33 groups it with D‴.48–D‴.60.
  - Its natural home would be a companion to Paper 17 once D‴.49 is final (the IEEE TAC venue named there), or a section of Paper 12.
- **What a referee would demand:**
  1. **A Proposition that contains only what is proved.** (i) and (ii) as a Proposition with explicit hypotheses: the frozen-K and quasi-static assumptions; M diagonalisable with distinct real multipliers; a derivation of the averaged loop as the zero-frequency part. (iii) becomes a numerical section ("Observations"), not part of a Proposition.
  2. **Tables and figures instead of prose** for EDGE, C7, C8 and C8N; (iii) is currently one 4,561-word paragraph.
  3. **The C8A/C8F attribution of 12 diagnostic false alarms** (about 1,500 words, rounds 98–134) compressed to two sentences in the text: "the coefficient-only replay misfires where the centre has drifted; a full factorial shows the centre decides at 11, the coefficient at 1". The detail goes to a supplement, and every "referee round N" and history clause ("the fifth statement on that window", "round 119's withdrawal was itself wrong") moves to the log.
  4. **Statistics.** Four draws per configuration and eleven at the threshold. A referee will ask for draw-level uncertainty on every count and on the s* test, whose null D60 already concedes has little power.
  5. **The framing.** The deterministic 16-tone forcing called "noise-induced" is labelled an analogy in D60; this must be kept prominent.
  6. **Positioning against bursting and parameter-drift theory:** the drift as an averaging result; the ejection as a basin-boundary crossing in a non-autonomous system.
- **Verdict.** Not yet submittable. The mathematics in (i) and (ii) is sound; the presentation of (iii) is a lab notebook.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3eg-1 | minor | **"Shows no interaction" (D60:176–179; Numerics C8F; md:683) holds for the outcome, not for the dynamics.**<br>• The largest distance is bimodal: passing corners 0.030–0.286, leaving 0.310–1.90. Within each class the spread is ≤ 1.9·10⁻³ and the gap between classes ≥ 0.280 [own, check_c8f]. So "+1.000 … +1.001, others ≤ 1.2·10⁻³" is the one-factor binary outcome restated, plus that spread; it is not a second test.<br>• On a non-saturating response (kernel 1's distance in period 1, before the alarm) at draw 2 from 58 [own, c8f_cont]: the named centre +1.75–1.80, group C −1.03–−1.09, and interaction 01100 +0.28–0.29 of the change.<br>• At draw 7, period 3: α₁ +1.20–1.24, named centre −0.20–−0.22.<br>• Cf. the interaction term in activation patching, which is negligible only where the response is locally affine [snippet, arXiv:2606.27510], and metric saturation [snippet, arXiv:2404.15255]. | **D60 (iii):** "shows no interaction in the outcome (the window's largest distance, bimodal, restates it)".<br>**Optional:** print one non-saturating response in `c8f_lines`, e.g. the pre-alarm distance or the named centre's threshold, and say that the dynamics before the decision do interact. |
| D3eg-2 | minor, latent | **Six single-point mutants of `c8f_job` survive all 14 relevant tests and leave the md unchanged** [own, mut/res_k.jsonl]:<br>• J01: the start with δ(p) instead of δ(p − k);<br>• J02: the start's noise from p;<br>• J03: the start one period longer;<br>• J04: δ from W[p − 2];<br>• J05: the largest distance without the last period;<br>• J09: group C over kernels 1 and 2 only.<br>Why: the stub in `test_round134_c8f_lines_cases_and_job` returns a fixed state for the start segment whatever its shift, length or δ, with a constant W, and only checks jb = 2.<br>Not equivalent on rerun [own, mut_rerun.log]: J01 changes draw 7's k = 2, 3 corners by 0.28 (outcome flips), J09 by ≤ 2.9·10⁻⁶. Both would be caught only after regenerating C8F, by `test_round134_c8f_on_the_raw`'s exact corner equalities. The authors' K11–K16 were in the corner loop, not the start. | **tests:** record the start-segment calls in the stub and assert shift = p − k, n = k and α₁ = raised(δ(p − k)); use a W that varies by period; a dev list whose last period carries the maximum; a case with jb = 0 asserting C sets both other formed centres.<br>**Mutant record:** add J01–J05, J09. |
| D3eg-3 | cosmetic | **J06 (`at` within 2 periods) and J07 (pass from R + 1, not R + 2) in `c8f_lines` survive** and are equivalent on the data: every corner leaves exactly at R or never [raw]. | Plant a corner leaving at R + 2 (must count as a pass) and at R ± 2 (must not count as at the alarm). |
| D3eg-4 | minor (journal) | **Readability of (iii)** [md]:<br>• one 4,561-word paragraph out of 10,879 words in D60;<br>• 26 "referee round" citations;<br>• history clauses in the statement (D60:188–193);<br>• the alarm attribution (D60:130–197) concerns a diagnostic's false alarms, not the dynamics. | Restructure as in Journal readiness, points 1–3. Keep the full text as an appendix or in the md. |
| D3eg-5 | minor (journal) | **D‴.60 has no paper.** No PAPERS.md "Rests on" names it, and Paper 17's "Missing" omits it [md]. | **PAPERS.md:** add it to Paper 17's "Missing" (or Paper 12's), or declare it internal. |
| D3eg-6 | cosmetic | **Stale records** [md]:<br>• exp_noise_resonance.py:2 says "Proposition D‴.60 (planned)";<br>• the Usage line (:47) and the parts list in the docstring stop at C8A, omitting C8F;<br>• C8A_ROUNDS (:1574) = (98 … 128), although round 134 added clauses to that paragraph (ctr2, psm), so md:681's header omits 134. | Update the docstring and Usage; add 134 to C8A_ROUNDS (md header only). |
| D3eg-7 | minor | **The headline mechanism is a section** (title, D60:72, Status 296–300).<br>• "the drift carries a coefficient into the set of raises that the coming noise ejects": the set is computed at the fit's centres (C7, C8).<br>• C8A/C8F show the run's own kernel 1 centre (period-mean offset −0.0713 … −0.0325, outward) or kernel 3's decides at 11 of the 12 isolated alarms [md, raw].<br>• D60 states this in (iii) but neither the title nor the Status qualifies the headline. | **Status (iii):** "the ejecting set a section of the basin at the features' centres; where a centre has drifted the section misfires (the 12 isolated alarms)". |
| D3eg-8 | cosmetic | **(ii)'s labels:**<br>• the Proof asserts "the constant part's response is the averaged loop's second order (… exact along the orbit)" without derivation, while Status cites B5 as a check;<br>• the bound Σ_v\|Pv\|\|g\|(1 + ρ^j)/\|ρ − e^{iΩT}\| assumes an eigenbasis (multipliers printed as real; distinctness not stated) [md]. | Mark the averaged-loop identity "derived in outline, checked (B5)", or add the derivation. State "M diagonalisable (distinct real multipliers, F)". |

## Sources

- **Factor separation:**
  - Lunt et al., "Multi-variate factorisation of numerical simulations", GMD 14 (2021) 4307 [snippet]: the Stein–Alpert notation, three-factor interaction terms. https://gmd.copernicus.org/articles/14/4307/2021/
  - arXiv:2609.11948, "Climate-model factor separation with Shapley values and efficient sampling" [snippet; status not checked]: identifies Stein–Alpert interactions with Boolean Möbius coefficients (Harsanyi dividends), the inversion C8F uses. https://arxiv.org/pdf/2609.11948
- **Interaction in patching designs:**
  - "The Curse of Multiple Mediators: Hidden Interaction Effects in Activation Patching", arXiv:2606.27510 [snippet]: the interaction term scales with the distance between clean and patched states, is negligible when the model is locally affine, and decomposes into pairwise and higher-order group terms. Relevant to D3eg-1: here a threshold makes the response non-affine, so the interaction shows up in the continuous dynamics and not in the outcome. https://arxiv.org/pdf/2606.27510
  - Heimersheim & Nanda, "How to use and interpret activation patching", arXiv:2404.15255 [snippet]: metrics can saturate and fail to show change past a threshold, which is the largest-distance issue. https://arxiv.org/pdf/2404.15255
- **Parameter drift and bursting** [snippet]:
  - "A Historical Perspective of Adaptive Control and Learning", arXiv:2108.11336. https://arxiv.org/pdf/2108.11336
  - Search summaries described drift along the manifold of explaining parameters under measurement noise until bursting, with deadzones, leakage and projection as remedies; I did not fetch the pages.
- **Noise-induced and rate-induced tipping** [snippet]:
  - escape times versus drift time scale: arXiv:2304.09160, https://arxiv.org/pdf/2304.09160
  - a forcing time scale with a resonance-like maximum of rate-induced tipping: arXiv:2603.15043, https://arxiv.org/pdf/2603.15043v1
  - enhancement of escape by a chaotic saddle in the basin: arXiv nlin/0309016, https://arxiv.org/abs/nlin/0309016v1
  - None combines a slow drift, quasiperiodic deterministic forcing and a finite-time basin as in C7/C8.
- **Moving-centre RBF adaptive control** [snippet]:
  - online centre selection with boundedness, CDC-ECC 2011: https://skoge.folk.ntnu.no/prost/proceedings/cdc-ecc-2011/data/papers/1046.pdf
  - dynamic centres for persistent excitation: https://researchconnect.buffalo.edu/en/publications/neural-network-based-control-of-large-structures/
  - No source analyses centre ejection under measurement noise, which is consistent with D60's novelty "in the instance".
- **LMS non-Wiener behaviour under narrowband interference** [snippet]: Virginia Tech thesis, https://vtechworks.lib.vt.edu/handle/10919/26772 . Consistent with D60's [snippet] citation; no Glover-specific source found.
- **Not rechecked this round** [memory]: Bode, Freudenberg–Looze, Yakubovich–Starzhinskii, Wereley, Hoppensteadt, Khalil, Oljača–Ashwin–Rasmussen, Giesl–Rasmussen and the other citations in the Status.

## Lesson

- **A Möbius decomposition of a saturating response restates the binary outcome.** To test "no interaction", decompose a response that has not yet crossed the threshold. Here it interacts (+0.29 at draw 2 from 58) while the outcome does not.
- **A stub that returns a fixed state hides every mutant upstream of it.** It must also record and assert its inputs, or the only guard left is a regeneration that nobody runs at test time.
- **Once a long attribution thread is closed, move it out of the statement.** A journal will read (iii) as a notebook until the history lives in the log and the claim fits in two sentences.

Files of this round, all under /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138:
- **scripts:** /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138/scripts/repro_c8f.py, /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138/scripts/check_c8f.py, /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138/scripts/c8f_cont.py, /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138/scripts/spec138.py
- **logs:** /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138/logs/ (fromraw.log, pytest_full.log, repro_c8f.log, check_c8f.log, dump_c8f.log, c8f_cont.log, mut_k.log, mut_h.log, mut_rerun.log)
- **mutant results:** /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138/mut/res_k.jsonl and /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref138/mut/res_h.jsonl
