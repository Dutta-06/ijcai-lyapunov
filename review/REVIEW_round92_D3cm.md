# Referee report D3cm, round 92: Proposition D‴.58 as revised after round 88

## Scope

- **Subject.** Commit 363855d on branch claude/vigilant-clarke-hf939b, which covers both 1ac3d05 (round 88 applied) and 363855d (D‴.61's round-91 results cited). I exported it with `git archive 363855d` into W/repo, where W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref92`.
  - WM_OUT was W/wmout, holding a copy of the raw.
  - Scripts are in W/scripts, logs in W/logs, mutant slots and results in W/mut and W/mut88.
- **Read in full.**
  - D58.md, `exp_width_mismatch.py`, `width_mismatch.md`, `test_width_mismatch.py` (55 tests);
  - REVIEW_round88_D3ci.md and the Round 88 entry of REVIEW_LOG.md;
  - `tools/mutants/round88_D58.py`, `round88_D58_referee.py` and `tools/mutants.py`;
  - `git diff 4bc7b07 363855d` for these files (the hunks, plus a word diff of D58).
- **Read for the citations.** D61.md's statement (i)–(v), its proof of (iii)–(iv), its Status and New, and the D58 and D61 rows of `review/drafts/README.md`.
- **Run settings.**
  - Python 3.13.18, numpy 2.5.4, scipy 1.18.1, from the repository's venv.
  - Every process ran with OMP_NUM_THREADS=1, PYTHONDONTWRITEBYTECODE=1 and `nice -n 10`, at most two of mine at a time.
  - pytest ran with `--basetemp` under W, also inside my copy of the mutant runner (W/scripts/mutants92.py, whose only changes are REPO, PY and `--basetemp`).
- **Nothing in the repository was modified.** I ran only `git archive/diff/log/show/status`.
  - Meanwhile another session moved HEAD to 89e9e28 (D‴.60). Neither of its commits touches a D‴.58 or D‴.61 file.
  - The working tree's modified moving_width files are not mine.
  - I stopped one of my own background scans (a slow first version of the domination-radius check) with TaskStop and replaced it with a vectorised one. All my processes have exited.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my computations; [snippet] web-search results only (arxiv.org does not resolve from this session); [memory] recalled facts.

## Summary

**Design.**
- Reproduced the md, the tests, and six jobs from the raw:
  - B5's censored start and its anchored solve;
  - A6 at the farthest rule and at the one near-cusp rule where the fits fail to bracket;
  - A6E at three (K, axis) cells.
- Recomputed with my own code:
  - the round's new quantities from the raw: A6's second fit, the bracket counts, the frame count, B5's censored ends;
  - D‴.61's domination radius, r*, the needed radii, and D‴.61 (iv)'s closed forms against A6's exact errors.
- Split A6's deviations by whether the rule has a frame, i.e. whether its limit is now proved.
- Ran 17 single-point mutants of this round's code plus a control, and re-ran round 88's P01, P02, R01, F24 and B02.

**Verdicts.**
- **Reproducibility is exact** [own]:
  - `--from-raw` regenerates the md byte for byte;
  - 55/55 tests pass in 53.3 s;
  - every rerun job is bit-identical to the raw in every field.
- **The round's numbers trace to the md and the raw** [own]:
  - the second fit: 6.1·10⁻⁵ and 1.4·10⁻⁴; +5.7·10⁻⁹ … +6.1·10⁻⁵ against −1.4·10⁻⁴ … −3.6·10⁻⁷; 13 of 13;
  - the frame count: 22 = 4 + 4 + 6 + 8, ρ = 1 at 4, |ρ| > 1 at 10, max t(v) 0.6345;
  - B5's marks: only (¼, 6, θ = 0) at κ = 2 and 4, with 2 accepted steps.
- **D58's citations of D‴.61 match D‴.61** in substance: the radius 1.8927, r* = 0.8053, σ*/σ < 9.27, the 22 rules, and "modulo a grid-checked domination condition". My own grid gives:
  - least radius 1.892733;
  - violation at |η| = 1.95 for arg ∈ [0.0864, 0.6174];
  - largest gaps −9.14·10⁻³, −3.29·10⁻³, −9.71·10⁻⁴ at |η| = 1.7, 1.78, 1.85;
  - needed radius 1.7773 at t(v) = 0.6345;
  - D‴.61's closed form over A6's e_K² at (2, ½·v_hs), K = 160: 1.19480, as D‴.61 prints.
- **The proofs are unchanged and right.** (i)'s sharp constant now follows from D‴.61 (iv): the half-squeeze frame ρ = t, v′ = 2, a = 2t is the cusp γx = γy = 1/(2t).
- **Round 88's issues are fixed or nearly so**: P01, P02, F24 and B02 are killed; R01 is md-equivalent. Residues are below.
- **No major issue.**

**Findings.**
1. **D3cm-1 (minor): A6's slow branch is summarised beyond the data, and since D‴.61 its worst deviation is at a proved rule.**
   - The Status's "the two fits bracketing it at every rule where t(v) is attained at τ²/v" is false at (½, Fock), v/v_hs = 0.9938, where both fits sit above (+3.0·10⁻⁷, +5.3·10⁻⁶). The md's 13 of 13 counts only |ln(v/v_hs)| ≥ 0.05.
   - "Within the fits' spread, 1.4·10⁻⁴" names the K⁻¹ fit's deviation; the spread there is 2.0·10⁻⁴.
   - The rule that sets both 6.1·10⁻⁵ and 1.4·10⁻⁴, (2, ½·v_hs), has a frame (ρ = −0.736, t(v) = 0.6345 < r*). Its limit is therefore proved, so the deviations measure the fits' bias.
   - The conjecture proper, the 14 rules without a frame, holds within 7.6·10⁻⁶ by the first fit and 4.5·10⁻⁵ by the second, bracketed at 5 of 5.
2. **D3cm-2 (minor): what D‴.61 gives D58 is stated unevenly.**
   - Statements (i) and (ii) say "proved by D‴.61 (iv)" without the grid proviso that the Status carries.
   - The Status's "(ii) … the rates themselves numerical (A3)" understates: the Fock and Mercer-zero rules always have a frame, so D‴.61 proves λ/(2 − λ) for σ*/σ < 4.19 and q for τ/σ < 6.52. Only moment matching's λ remains numerical.
   - "At the 22 of A6's 36 rules" understates a theorem that covers the whole frame interval where t(v) < r*.
   - New's "proved with D‴.61 for the rules with a frame" overstates: it omits t(v) < 0.8053, which fails near the interval's ends once σ*/σ > 3.045.
   - r* = 0.8053 and 9.27 are rounded up from 0.805257 and 9.2699.
   - Below t(v) = 0.0966 no grid check is needed at all: 7 of the 22 rules, and (i)'s constant for σ*/σ < 1.214.
3. **D3cm-3 (cosmetic): 5 of 17 new mutants survive, all equivalent on the data.** They are M03, C02, R11, R12 and R15; none changes a published result. The tests kill 12. Round 88's accounting holds.
4. **D3cm-4 (cosmetic): stale labels and unqualified phrases.**
   - The md's B4 line and the code still call the rounded best end "the double parameters' own floor" (D3ci-2's relabel was not done).
   - The md's headings still say "conjectural" for A6 and for D, though all four of D's t lie below r*.
   - The README row still says "conditional proof".
   - "Bounded above in K" (A6E) rests on three K, and on the imaginary axis the maximum is still rising, toward √2.
   - D58's two "to 2.2·10⁻¹⁶" lack the md's new qualification (D3ci-5f).
   - D‴.61 writes the frame interval as (τ²/σ*², τ²), with σ = 1 implicit; D58 has it right.
5. **D3cm-5 (minor): an uncited related literature.** Gauss–Hermite quadrature Fourier features for the Gaussian kernel come with Hermite-expansion error analyses (Mutný–Krause 2018; SLEIPNIR, arXiv:2003.02658; "Spectral Bounds for Kernel Quadrature", arXiv:2609.00553) [snippet]. They bear on (v)'s cosine sums and on New's exact-series claim, and are not checked.

## What I checked and how

**Reproduction.**
- `--from-raw` into W/wmout, then `cmp` with the committed md: identical (md5 1ce8eae5…; raw 959dc055…, unchanged since round 88).
- The test file: 55 passed in 53.31 s.
- `rerun.py`, comparing with `==` in every field except `secs`:
  - **B5, (¼, 6, θ = 0), direct** (the censored start): identical; accepted steps 51, 42, 40, 42, 35, 37, 2, 2.
  - **B5, (¼, 6), anchored**: the κ = ¼ end behind D58's "0.0972"; identical.
  - **A6, (2, ½·v_hs), to K = 160**: identical, 7 s.
  - **A6, (½, Fock), to K = 160**: identical, 55 s; precision up to 461 digits.
  - **A6E**: K = 40 on the real axis (all 20 s); K = 80 on the imaginary axis and K = 20 on the real axis, both at s = 2; identical.
- **B5's censoring [raw].**
  - The two starred ends have 2 accepted steps in 39 iterations. They lie 1.0000005904 times the stored optimum, i.e. 5.9·10⁻⁷ worse than their start; the excess +3.1·10⁻⁵ is the start's.
  - Every other direct end has ≥ 17 accepted steps, and every anchored end ≥ 64.
  - The star affects only the [¼, 4] slope 2.64. D58 quotes no number that includes it.

**A6 refit (`a6fits.py`, my own design matrices on the raw's ratios).**
- Both fits through K = 80, 120, 160, at all 36 rules: the md's numbers exactly (table in D3cm-1).
- **Bracketing.**
  - Over all 16 τ²/v rules: 15. The exception is (½, Fock), at +3.02·10⁻⁷ and +5.34·10⁻⁶.
  - Over the 3 near-cusp τ²/v rules: (¼, Fock) −2.9·10⁻⁵ against +4.2·10⁻⁶, and (¼, Mercer-zero) +3.0·10⁻¹⁰ against −2.1·10⁻⁸, both bracket.
- **Split by frame.**
  - Frame rules off the cusp (18): max |fit 1| 6.07·10⁻⁵, max |fit 2| 1.39·10⁻⁴, max spread 2.00·10⁻⁴.
  - No frame (14): 7.55·10⁻⁶, 4.45·10⁻⁵ and 5.20·10⁻⁵ (worst at τ/σ = 1, ½·v_hs).
- **Other ansätze at (2, ½·v_hs)**, whose limit D‴.61 proves:
  - L + aK^{−3/2} + bK^{−2} + cK^{−5/2} (4 K): −2.9·10⁻⁵;
  - L + aK^{−3/2} + bK^{−2}: −3.8·10⁻⁵;
  - free exponent (p = 1.368): −5.7·10⁻⁵.
  - No ansatz comes within 2.9·10⁻⁵.
- **The excess's local exponents** steepen from (80, 120) to (120, 160) at 13 of 13 far τ²/v rules, e.g. 1.248 → 1.525 at (2, ½·v_hs). At (¼, ½·v_hs) they are non-monotone: 2.69, 2.14, 1.86, 1.86, 1.89.
- **Frames.**
  - ρ = 1 + 2σ²/τ² − 2σ²/v gives |ρ| < 1 at 22 rules (4, 4, 6, 8), ρ = 1 at the 4 moment-matching rules, and |ρ| > 1 at 10.
  - The largest |ρ| among the frame rules is 0.9411; the largest t(v) is 0.63451.

**D‴.61's numbers as D58 cites them** (`domrad2.py`, my own vectorised `g_η + Λ⁺`, written from D‴.61's formulas).
- **Least radius.**
  - 361 angles × 4001 s: 1.892735 at arg 0.3447.
  - Refined over 81 angles: 1.892733 at arg 0.3432.
  - The same with 40001 s and a bounded local maximisation.
- **Sector and gaps.** Violated at |η| = 1.95 for arg ∈ [0.0864, 0.6174] (3601 angles). Largest gaps at |η| = 1.7, 1.78, 1.85: −9.137·10⁻³, −3.294·10⁻³, −9.708·10⁻⁴.
- **Thresholds.**
  - c_need(0.63451) = 1.77734.
  - c_need(0.8) = 1.88952, which is 0.0032 inside the radius.
  - r* = 0.805257 from 1.892733 (0.805203 from 1.8927 as printed).
  - σ*/σ < 9.2699.
- **Closed forms** (`d61pred.py`, (iv)'s forms with (v)'s constant, against A6's e₁).
  - (2, ½·v_hs): 1.32896, 1.24341, 1.19480 at K = 80, 120, 160 (D‴.61's R: 1.19480).
  - Cusps: 0.99991 (τ/σ = 2), 1.00062 (1), 1.00076 (½).
  - Imaginary-axis rules: 0.999–1.03.
  - So the real-axis branch's slow approach is D‴.61's slow real-axis convergence.
- **The unconditional range** (`uncond.py`). D‴.61's condition is void for |η| ≤ 1, and c_need < 1 for r < 0.09665. That covers 7 of A6's 22 frame rules (τ/σ = ¼: Fock, Mercer-zero, Mercer frame, v_hs; τ/σ = ½: Fock, Mercer frame, v_hs), and the half-squeeze rule for σ*/σ < 1.2140.
- **The Fock and Mercer-zero rules have frames at every τ/σ.**
  - Fock: ρ = 0.
  - Mercer-zero: ρ = 1 − 2(√(1 + 2x) − 1)/x ∈ (−1, 1), x = τ²/σ².
  - Their rates fall below r* for σ*/σ < 4.188 (Fock) and τ/σ < 6.517 (Mercer-zero).
  - The frame interval's ends, where t(v) → λ, exceed r* once σ*/σ > 3.045.

**A6E [raw].**
- Real axis: the minimum sits at s = 1.6, 1.7, 1.7; the maximum at s = 2 falls, 1.0831, 1.0516, 1.0207.
- Imaginary axis: the maximum at s = 2 rises, 1.4075, 1.4109, 1.4126, toward √2 = 1.41421.
- On the real axis at s = 1.9, over D‴.61's limit: 2.976, 2.716, 2.438.

**Mutation testing.**
- `W/mut/spec92.py`: 17 mutants plus a control (Y00) of the code changed in 4bc7b07..363855d, none in round 88's lists:
  - `b5_lines`' marks and footnote: M01–M04;
  - `a6_lines`' second fit, bracket count and farthest rule: S01–S06;
  - the cusp label: C02;
  - the frame count: R11–R16.
- Then round 88's spec at X00, P01, P02, R01, F24 and B02.

**Literature.** Six web searches (Sources). arxiv.org is unreachable by fetch (ENOTFOUND), so all are [snippet].

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| (i) "proved by D‴.61 (iv) for t < 0.8053 (σ*/σ < 9.27)" (:14) | Matches D‴.61 (iv). The statement omits "modulo the grid", the bounds are rounded up, and t < 0.0966 is unconditional | [own]; D3cm-2 |
| (ii) "(A6; proved by D‴.61 (iv) for the rules with a Mehler frame, v ∈ (τ²σ²/σ*², τ²), where t(v) < 0.8053, conjectural for the others)" (:23–24) | Correct; the interval is right (D‴.61's own (τ²/σ*², τ²) is σ = 1). Proviso missing; "where" is ambiguous | [own]; D3cm-2 |
| Status (i): constant for t < 0.8053, modulo grid; conjectural beyond (:169) | Correct (D‴.61: conditional on P beyond) | D61 Status |
| Status (ii): "the rates themselves numerical (A3)" (:180) | Understates: Fock and Mercer-zero are proved where < r* | D3cm-2 |
| Status (ii): 6.1·10⁻⁵ / 1.4·10⁻⁴ / 3.1·10⁻⁷; "bracketing it at every rule where t(v) is attained at τ²/v" (:181–182) | Numbers right; "every" false at (½, Fock) | [raw], [own]; D3cm-1 |
| Status (ii): "made a proof … at the 22 of A6's 36 rules …, modulo …; moment matching and the 10 with \|ρ\| > 1 … A6's numbers alone" (:183–184) | Correct; the restriction to 22 rules undersells the theorem | [own] |
| A6 paragraph: 13 rules, +5.7·10⁻⁹ … +6.1·10⁻⁵, −1.4·10⁻⁴ … −3.6·10⁻⁷, 13 of 13, K^{−1.52…−1.98} "steepening", +8.1·10⁻⁶ … +4.0·10⁻⁴; "within the fits' spread, 1.4·10⁻⁴" (:232–237) | Numbers right and steepening confirmed (13/13). "Spread" mislabelled (2.0·10⁻⁴); the 1.4·10⁻⁴ sits at a proved rule | [md], [own]; D3cm-1 |
| Sketch: A6E "bounded above" 0.677–1.021, 1.001–1.413; minimum 0.804, 0.741, 0.677 → (1 − s²/4)^{1/2} (:240–242) | Numbers right; "bounded above" observed at three K only | [raw]; D3cm-4 |
| "D‴.61 (iii) proves on both axes and on \|η\| < 1.8927 … so that (iv) proves … all 22 …; moment matching (4) and the 10" (:246–249) | Correct ("both axes" means the intervals \|η\| < 2) | [own] |
| (iii) B4/B5: rounded best end 1.1·10⁻¹² "a scale … 0.0972"; κ = 0 at 0.954 and 5.77 "which the parameters' rounding does not explain"; "(at κ = 0)" (:83–91) | Correct; the md label still says "floor" | [md], [raw]; D3cm-4 |
| Status (iii): "growing about as the square … (by least squares), ≲ 10⁻¹², below κ = ½ its solves stalling, unexplained" (:187–189) | Correct (flat stops 0/3, 0/3, 2/3, 3/3 … at κ = 0, ¼, ½, 1) | [md] |
| New: "every scaled rule's rate — proved with D‴.61 for the rules with a frame" (:265) | Overstates: omits t(v) < 0.8053 and the grid | D3cm-2 |
| Numerics A6, A6E, B5 (:300–325) | Trace to the md | [md] |
| Proofs (i), (v); (iii)–(iv) unchanged | Unchanged this round; (i)'s constant via D‴.61's cusp form checked | [own] |

## Status of the round-88 issues

| D3ci | status | note |
|---|---|---|
| 1 | fixed, with residues | The second fit and the bracket count are printed and "above it" is withdrawn [own: identical]. Residues: the Status's "every rule"; "spread" mislabelled; the worst deviation now sits at a proved rule (D3cm-1). The "same rule, two limits" item (τ/σ = ¼, Mercer frame: 65.054019 from A5, 65.052913 from A6) is not addressed. That rule's limit 65.055507 is now proved unconditionally (c_need 0.51), so both fits measure their own bias. |
| 2 | fixed in D58, residue in md/code | D58 :85 says "a scale, not the double parameters' floor" (0.0972 ✓). md :313 and the code (:751, :947) still say "the double parameters' own floor" (D3cm-4). |
| 3 | fixed | D‴.61 is cited (now its round-91 revision); the domain and the 22/4/10 count are confirmed [own]; A6E's minimum is falling toward (1 − s²/4)^{1/2}. Residues: D3cm-2, and A6E's "bounded above in K" from three K (D3cm-4). |
| 4 | fixed | Re-run: P01 and P02 are killed by `test_round88_precision_and_the_md`; F24 and B02 by `test_round88_the_second_fit_and_the_censored_ends`; R01 survives, the md unchanged (it alters only the stored `arg` at the three tie rules on a rerun), as the log says. |
| 5 | fixed, one residue | (a)–(e) and (g) done. For (f) the md is qualified, but D58 :215 and :295 still say "to 2.2·10⁻¹⁶" unqualified (D3cm-4). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cm-1 | minor | **"Every rule".** Status :182 says the fits bracket "at every rule where t(v) is attained at τ²/v". They do at 15 of 16: at (½, Fock), v/v_hs = 0.993808, both lie above (+3.02·10⁻⁷, +5.34·10⁻⁶) [raw, own]. The md's 13 of 13 is for \|ln(v/v_hs)\| ≥ 0.05.<br>**"Spread".** :237 "within the fits' spread, 1.4·10⁻⁴" (and the md's "the limits' spread over the two fits"): 1.4·10⁻⁴ is fit 2's deviation; \|L₁ − L₂\|·t(v) = 2.0·10⁻⁴ at (2, ½·v_hs) [own].<br>**What D‴.61 changes.** That rule has a frame (ρ = −0.736, t(v) = 0.6345), so its limit is 1/t(v) by D‴.61 (iv), modulo the grid. Its deviations (+6.1·10⁻⁵, −1.4·10⁻⁴; other ansätze −2.9 … −5.7·10⁻⁵) are the extrapolation's bias, consistent with e_K²/prediction = 1.329, 1.243, 1.195 at K = 80, 120, 160 [own].<br>**Where the conjecture still rests on numbers.** The 14 rules without a frame: within 7.6·10⁻⁶ by fit 1 and 4.5·10⁻⁵ by fit 2 (spread ≤ 5.2·10⁻⁵), bracketed at 5 of 5 [own]. | Status :182: "at every rule away from the cusp (\|ln(v/v_hs)\| ≥ 0.05) where …". :237: replace with "at the 18 frame rules off the cusp, whose limits D‴.61 proves, the two fits miss it by up to 6.1·10⁻⁵ and 1.4·10⁻⁴ (their bias); at the 14 without a frame, where the formula is still conjectural, within 7.6·10⁻⁶ and 4.5·10⁻⁵, bracketing it at 5 of 5". `a6_lines`: print this split by frame, and call 1.4·10⁻⁴ "the second fit's largest deviation". Optionally cross-refer D‴.61's account of the slow real-axis branch. |
| D3cm-2 | minor | **(a) Proviso.** Statements (i) :14 and (ii) :24 say "proved by D‴.61 (iv)" without "modulo the grid-checked domination condition"; the Status :169 and :183 carry it.<br>**(b) Status (ii) :180 understates.** "The rates themselves numerical (A3)" — but the Fock rule (ρ = 0) and the Mercer-zero rule (ρ ∈ (−1, 1) for every τ/σ) always have frames. D‴.61 (iv) proves λ/(2 − λ) for σ*/σ < 4.188 and q for τ/σ < 6.517 [own]; only moment matching (ρ = 1) stays numerical.<br>**(c) Scope.** Status :183's "at the 22 of A6's 36 rules" restricts a theorem that holds on the whole frame interval with t(v) < r*. New :265's "proved with D‴.61 for the rules with a frame" drops t(v) < 0.8053, which fails near the interval's ends once σ*/σ > 3.045 [own], and drops the proviso.<br>**(d) Rounding.** r* = 0.8053 and 9.27 are rounded up: the grid's radius 1.892733 gives r* = 0.805257, σ*/σ < 9.2699 [own; D‴.61's `r_star`].<br>**(e) No grid needed below 0.0966.** The condition is void for \|η\| ≤ 1, and c_need < 1 for r < 0.09665. So the grid proviso is not needed for 7 of the 22 rules, nor for (i)'s constant at σ*/σ < 1.214 [own]. | (a) Add the proviso in (i) and (ii), or "(D‴.61 (iv), modulo its grid-checked condition)". (b) Status (ii): "the rates λ/(2 − λ) and q proved by D‴.61 (iv) where they are below r*; moment matching's λ numerical". (c) Status: "every v ∈ (τ²σ²/σ*², τ²) with t(v) < r*, A6's 22 frame rules among them"; New: add "and t(v) < 0.8053, modulo …". (d) State t < 0.8052, σ*/σ < 9.26 (here and in D‴.61). (e) Optionally: "unconditionally for t(v) < 0.0966". |
| D3cm-3 | cosmetic | **Mutants.** 17 plus control: Y00 passes; the tests kill 12 (M01, M02, M04, S01, S03–S06, R13, R14, R16 also change the md; S02 only the tests, through a 10⁻¹⁶ residual on planted exact fits).<br>**5 survive both checks, all equivalent on the data; none changes a published result:**<br>• M03: the footnote over anchored rows (all anchored acc ≥ 64);<br>• C02: the cusp label by v_hs == 1.0 (exact in double);<br>• R11: \|ρ\| < 1 − 10⁻³ (largest frame \|ρ\| 0.9411);<br>• R12: ρ = 1 within 10⁻⁶ (exactly 1 in double);<br>• R15: any/all `'v'`. | Optional: a planted anchored end with acc ≤ 2 (M03); a planted rule with \|ρ\| = 0.9995 (R11). Record the rest as equivalent. |
| D3cm-4 | cosmetic | **Labels and phrasing.**<br>• md :313 label and the code (:751 docstring, :947) still say "the double parameters' own floor" (D3ci-2's relabel).<br>• md :128 "(referee round 84; conjectural)" and :510 "(the sharp constant, conjectural)" are stale: 22 of A6's rules, and all four of D's t (≤ 0.8 < r*), are now covered by D‴.61.<br>• `review/drafts/README.md`'s D58 row still says "its conditional proof".<br>• :240 "stays bounded above on both axes" and md "Bounded above in K" rest on K = 20, 40, 80. On the imaginary axis the maximum at s = 2 rises (1.4075, 1.4109, 1.4126) toward √2, and D‴.61 (iii)(b) holds only on compact subsets of \|η\| < 2.<br>• :246 "on both axes" means the intervals \|η\| < 2.<br>• :215 and :295 "to 2.2·10⁻¹⁶" lack the md's "the resolution of the stored doubles".<br>• D‴.61 (iv)–(v) write v ∈ (τ²/σ*², τ²) with σ = 1 implicit; D58's (τ²σ²/σ*², τ²) is the right form. | Relabel "(a scale; …)"; update the headings and the README row; "(observed at K = 20, 40, 80)"; "on both axes' intervals \|η\| < 2"; add the qualification; tell D‴.61 to restore σ. |
| D3cm-5 | minor | **Uncited related work.** Gauss–Hermite quadrature Fourier features for the squared-exponential kernel:<br>• Mutný–Krause (NeurIPS 2018);<br>• SLEIPNIR (arXiv:2003.02658);<br>• Li–Balakirsky–Mak (arXiv:2310.14544);<br>• "Spectral Bounds for Kernel Quadrature" (arXiv:2609.00553: Gauss–Hermite features with Hermite-function bases, a bound ≈ τⁿ/√(n!)).<br>One search summary attributes to an arXiv paper in this set an error expansion whose sum "starts at 2n" when the rule is exact for low-order Hermite functions [snippet; paper not identified]. This bears on (v)'s cosine sums at Hermite-zero frequencies and on New's "exact defect series". | Classical/New: cite as [snippet] and add to "to be checked against the papers"; compare (v) with QFF's cosine features. |

## Sources

- [Karvonen–Oates–Girolami, Integration in RKHS of Gaussian kernels, arXiv:2004.12654](https://arxiv.org/pdf/2004.12654); [Math. Comp. 90 (2021)](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/) [snippet]. Two-sided bounds exp(−c₁N^{1/d})N^{1/(4d)} ≲ · ≲ exp(−c₂N^{1/d})N^{−1/(4d)}, c₁ > c₂: no exact rate.
- [Karvonen–Särkkä, Gaussian kernel quadrature at scaled Gauss–Hermite nodes, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532); [Aalto record](https://aaltodoc.aalto.fi/items/0c4e072c-7344-4a4f-ac45-7a6e2cca9121) [snippet].
- [Scaling Optimized Hermite Approximation Methods, arXiv:2412.08044](https://arxiv.org/html/2412.08044v3) [snippet]. Optimal scaling by balancing spatial and frequency truncation; exp(−CN^{2/3}) for scaled Gauss–Hermite reported as claimed in unpublished work.
- [Derevianko–Plonka, Differential approximation of the Gaussian by short cosine sums with exponential error decay, arXiv:2307.13587](https://arxiv.org/pdf/2307.13587) [snippet]. Frequencies at zeros of a scaled Hermite polynomial; exponential decay in a weighted L²; published in Sampl. Theory Signal Process. Data Anal. (2026) per a citing reference.
- [Approximation by short exponential sums with geometric error decay based on Gauss quadrature, arXiv:2606.03855](https://arxiv.org/pdf/2606.03855) [snippet: title only; authors not seen].
- [Mutný–Krause, Efficient High Dimensional Bayesian Optimization with Additivity and Quadrature Fourier Features, NeurIPS 2018](https://papers.nips.cc/paper/2018/hash/4e5046fc8d6a97d18a5f54beaed54dea-Abstract.html); [SLEIPNIR, arXiv:2003.02658](https://arxiv.org/pdf/2003.02658); [Trigonometric QFF, arXiv:2310.14544](https://arxiv.org/pdf/2310.14544); [Spectral Bounds for Kernel Quadrature, arXiv:2609.00553](https://arxiv.org/pdf/2609.00553) [snippets] (D3cm-5).
- [Schrödinger picture analysis of the beam splitter: the Janszky representation, arXiv:1806.05748](https://arxiv.org/pdf/1806.05748) [snippet]. The squeezed vacuum as a Gaussian-weighted line integral of coherent states, as D58's squeezed-state citation says; no Gauss–Hermite discretisation was found.
- [Mehler kernel (Wikipedia)](https://en.wikipedia.org/wiki/Mehler_kernel) [snippet].
- **[memory].** Laplace's method; the bias of a fixed-exponent extrapolation when the true exponent lies between the assumed ones.

## Lesson

When a theorem arrives that covers part of a numerical test bed, re-partition the evidence. At the cases now proved, the deviations measure the instrument; here the two fits miss a proven limit by up to 1.4·10⁻⁴. The conjecture's evidence is what remains at the unproved cases, here 4.5·10⁻⁵ at the 14 rules without a frame. Quoting the worst deviation over both kinds mixes the two.

A citation must carry the cited result's hypotheses everywhere it is used: statement, Status and New alike. Its thresholds must be rounded toward safety: a strict bound of 0.805257 is 0.8052, not 0.8053. It also pays to read the hypotheses for where they are void, since part of the claimed range needs no grid check at all.
