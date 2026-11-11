# Referee report D3db, round 107: Proposition D‴.58 as revised after round 101

## Scope

- **Subject.** Commit ce7cfac on claude/vigilant-clarke-hf939b, which applies round 101 (D3cv), diffed against 0cb96bd.
  - Exported with `git archive ce7cfac | tar -x -C W/repo`, W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref107`.
  - WM_OUT = W/wmout and GP_OUT = W/gpout, each with a copy of its raw.
  - Scripts are in W/scripts, logs in W/logs, mutants in W/mut, W/mut2 and W/mut101.
- **Read in full.**
  - D58.md.
  - Every hunk of `exp_width_mismatch.py` changed since 0cb96bd: `out_e`, `D61_RSTAR_V`, `a6_lines`' new cells and sentences, the threshold line, a5's label, the headings, `a6_job`'s third argument, `A6X_RULES`/`A6X_K`, `part_A6X`, `a6x_lines`, and the report/main wiring.
  - width_mismatch.md; `tests/test_width_mismatch.py` (the three `test_round101_*` tests line by line); `tools/mutants/round101_D58_own.py`.
  - The Round 101 entry of REVIEW_LOG.md and REVIEW_round101_D3cv.md.
  - From D‴.61: D61.md (iii)(b), (iv), Status, New and Numerics; gh_peano.md parts D and V. For context I also read REVIEW_round108_D3dc.md, which archives D‴.61's round 108.
  - `git diff 0cb96bd ce7cfac` of these files, plus a word diff of D58.
- **Run settings.**
  - Every run used the repository's venv with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` and a `timeout`, at most two heavy processes of mine at a time.
  - pytest always ran with `--basetemp` inside W. My runner, W/scripts/mutants107.py, differs from `tools/mutants.py` only in REPO (= W/repo), PY and `--basetemp`.
  - Scripts run with `-I` wrote `__pycache__` inside W only.
- **The repository was not modified.** I ran only `git archive/log/diff/show/status/rev-parse`. I stopped no process, and all of mine have exited.
- **The repository moved during the review.**
  - The working tree holds other sessions' uncommitted edits (exp_gh_peano.py, D60 files), none of them D‴.58's.
  - HEAD moved from ce7cfac to 8139e12 (7 commits). One of them, 88ef79f, makes D‴.61's branch and bound rigorous in ball arithmetic (relevant to D3db-1).
  - None of the 7 touches D58.md, exp_width_mismatch.py, width_mismatch.md, its raw or test_width_mismatch.py.
- **Tags.**
  - [md]: the committed md or draft.
  - [raw]: the committed raw.
  - [own]: my own code, or the repository's functions on modified inputs.
  - [snippet]: web-search results only.
  - [memory]: recalled facts.

## Summary

**Design.**
- Regenerated both mds from the raw and ran the 61 tests.
- Reran five A6/A6X jobs at τ/σ = 2 against the raw, one of them with 40 more digits.
- Recomputed every number the revision added:
  - A6's cells, distinct counts, the outlier ratios, the slow-branch list and its neighbour reading, and the range ends, using exact-rational fits from the stored doubles, with ρ and t(v) in 50-digit Decimal from the scalings' definitions;
  - the tier-2 thresholds in closed form;
  - A6X's fits three ways (float solve, exact from the doubles, exact from the stored 25-digit errors), its intermediate windows and four-term fits.
- Planted limit offsets into the raw and ran the repository's `a6x_lines` on the result, to see what A6X's verdict tolerates.
- Cross-read D‴.61's status of part V at d5d8438 (the state ce7cfac cites) and at HEAD.
- Ran 24 own mutants plus a control, a second pass without the md regression, and round 101's spec as shipped.
- Did six web searches.

**Verdicts.**
- **Reproducibility is exact.** Both mds are byte-identical. 61/61 tests pass. All five reruns are identical to the raw, and the +40-digit run gives the same doubles.
- **Every new number is right.**
  - The cells and "(8, 7 distinct)".
  - The outlier ratios 41.4 and 2.0.
  - The 13-rule slow list, and "4 of the 5 between their proved neighbours'" (the fifth at 0.53–0.85 of its neighbours).
  - The outward range ends.
  - The tier-2 thresholds 19.924, 19.899, 140.712 and 14.107.
  - A6X's fits and shrink factors.
  - Three residues: one range end still rounded inward, one shrink factor printed past the data's resolution, and 198.99 for 199.
- **Question 3 (does A6X calibrate?).** A6X reproduces round 101's extension to every printed digit. But its computed verdict, "every deviation shrinks as the window moves out" (each shrink factor > 1 between two windows), is a necessary condition, not a calibration: with the outlier's limit planted 5·10⁻⁸ off 1/t(v) it still prints that verdict.
  - The raw already holds the calibration: the last window's deviations are 0.26–2.70 and 0.40–1.32 times the proved rule's bias; at the later windows the shrink factors match the proved rule's; and four-term fits through K = 160–320 put every limit within 6.5·10⁻¹³ of 1/t(v).
  - It covers 4 of the 9 fast-branch rules without a frame, and only 1 of the 5 with ρ > 1.
- **Question 4 (part V's weight).** D58 does not say it consistently, and twice says something false.
  - Since d5d8438, 32 minutes before ce7cfac, D‴.61 calls V rigorous in ball arithmetic, with its branch and bound the only floating-point step. Round 108 confirmed this.
  - D58's Status (i) still says both steps are "in floating point with a slack, not interval arithmetic", and the proof sketch says "a floating-point valley certificate".
  - Seven further places, all pinned by tests, carry "modulo also its valley certificate".
  - At ce7cfac the second tier therefore rests on exactly the first tier's proviso. At HEAD (88ef79f) D‴.61 calls both steps rigorous.
- **Tests and mutants.** Round 101's 19 mutants are all killed. Of mine, 18 are killed and 6 survive, all equivalent on the data; 4 of them weaken A6X's verdict. Nine of the 16 md-changing mutants are caught only by the md regression; two of those flip the published verdict, and one drops the A6X line from the report.
- **No major issue.**

**Findings.**
1. **D3db-1 (minor).** D58 and the md give D‴.61's part V the weight of an unverified floating-point step. Two statements are false. The tiers should merge, and 198.99 should read 199.
2. **D3db-2 (minor).** A6X's verdict does not calibrate the fast branch's unproved rules. "The four fast-branch rules without a frame" reads as all of them, but A6X covers 4 of 9. The calibration is already in the raw.
3. **D3db-3 (cosmetic).** Rounding residues: the exponent range "2.00 … 2.04" (true lower end 1.9961), A6X's ranges rounded to nearest, and "56280.4" printed past the stored doubles' resolution.
4. **D3db-4 (cosmetic).** Wording, counts and staleness: "20.7–20.7 … at the 1 proved"; counts that include Fock ≡ Mercer-zero twice; A6X missing from the Numerics; HANDOFF stale.
5. **D3db-5 (cosmetic).** Tests: A6X's verdict is pinned only by the md regression and an extreme planted case; the A6X test reads the committed md, not the report; `D61_RSTAR_V` is not pinned to D‴.61.

## What I checked and how

**1. Reproduction** [md, raw, own]
- `--from-raw` into W/wmout and W/gpout, then `cmp`: both identical.
  - width_mismatch.md md5 fb278cbc…; its raw f1079628….
  - gh_peano.md md5 0b478552….
- Full test file: 61 passed in 60.18 s.
- `W/scripts/rerun.py` runs the repository's `a6_job` and compares every field but `secs` with the raw by `==`:

| job | time (machine load 5–7) | result |
|---|---|---|
| A6 (2, moment matching) | 12.4 s | identical |
| A6 (2, 2·v_hs) | 6.8 s | identical |
| A6X (2, moment matching), K = 200, 240, 320 | 177.7 s | identical |
| A6X (2, 2·v_hs) | 182.2 s | identical |
| A6X (2, MM) with A6_MARGIN = 100 (+40 digits) | 133.2 s | ratios the same doubles at all three K; e₁ the same to all 25 stored digits |

**2. A6's new cells and sentences** (`W/scripts/a6check.py`) [own, raw]
- Method: three-point fits by Cramer's rule in Fractions from the stored doubles. ρ = 1 + 2/τ² − 2/v, t(v) and the attaining argument in 50-digit Decimal.
- Checks: the stored rates agree to < 10⁻¹⁵; every argument agrees except the cusp tie; 22 rules have a frame; ρ = 1 at 4 rules.
- **Cells** (each matches md :171):

| cell | n (distinct) | max \|fit 1\| | max \|fit 2\| | brackets |
|---|---|---|---|---|
| frame, τ²/v, far | 8 (7) | 6.0734·10⁻⁵ | 1.3900·10⁻⁴ | 8 |
| frame, τ²/v, near | 3 | 2.9114·10⁻⁵ | 5.3392·10⁻⁶ | 2 |
| frame, vσ*²/τ², far | 4 | 2.0981·10⁻⁹ | 1.7348·10⁻⁷ | 0 |
| frame, vσ*²/τ², near | 3 | 3.9877·10⁻⁵ | 3.6666·10⁻⁵ | 0 |
| no frame, τ²/v, far | 5 | 7.5472·10⁻⁶ | 4.4490·10⁻⁵ | 5 |
| no frame, vσ*²/τ², far | 8 | 8.7052·10⁻⁸ | 3.5976·10⁻⁷ | 4 |
| no frame, vσ*²/τ², near | 1 | 5.2612·10⁻¹¹ | 3.7751·10⁻⁸ | 0 |

  Fock ≡ Mercer-zero at τ/σ = 2 (v = 4/3) is the only coincidence.
- **The outlier.** Both maxima are at (¼, 2·v_hs). The ratios to the proved fast-far maxima (2.0981·10⁻⁹ and 1.7348·10⁻⁷, both at (2, Mercer frame)) are 41.4908 and 2.0738, printed rounded down as 41.4 and 2.0 ✓.
- **The slow-branch list** matches the md's 13 entries to every printed digit, in order. All 13 have fit 1 > 0 and fit 2 < 0.
- **The neighbour reading.**
  - The four rules (½, ½·v_hs), (¼, ½·v_hs), (½, 0.8·v_hs) and (¼, 0.8·v_hs) lie between their proved neighbours on both fits. The closest call is (¼, 0.8·v_hs) at 1.23 times its lower neighbour (½, Mercer-zero).
  - (1, ½·v_hs) lies at 0.684 and 0.533 times (2, 0.8·v_hs) and (2, Fock) on fit 1, and at 0.854 and 0.724 on fit 2: "below both, by less than a factor 2" ✓ (the largest factor is 1.88).
  - "No proved rule lies at τ/σ = 1/4" is structural, not a sampling accident. A frame needs v/v_hs ∈ (σ/σ*, σ*/σ). At τ/σ = ¼ that means |ln(v/v_hs)| < ln 1.0308 = 0.0303, below A6_NEAR = 0.05, so no far rule at τ/σ = ¼ can ever have a frame.
- **The "By the argument" ranges** (md :173). True ends against printed:

| branch (far) | quantity | true ends | printed | outward? |
|---|---|---|---|---|
| τ²/v | fit 1 | +5.652·10⁻⁹ … +6.0734·10⁻⁵ | +5.6e-09 … +6.1e-05 | ✓ |
| τ²/v | excess at K = 160 | +8.0717·10⁻⁶ … +4.0319·10⁻⁴ | +8.0e-06 … +4.1e-04 | ✓ |
| τ²/v | fit 2 | −1.3900·10⁻⁴ … −3.5905·10⁻⁷ | −1.4e-04 … −3.5e-07 | ✓ |
| vσ*²/τ² | fit 1 | −2.0981·10⁻⁹ … +8.7052·10⁻⁸ | −2.1e-09 … +8.8e-08 | ✓ |
| vσ*²/τ² | excess at K = 160 | −3.6234·10⁻⁶ … −1.7780·10⁻⁶ | −3.7e-06 … −1.7e-06 | ✓ |
| vσ*²/τ² | fit 2 | −1.7348·10⁻⁷ … +3.5976·10⁻⁷ | −1.8e-07 … +3.6e-07 | ✓ |
| τ²/v | exponent | 1.5246 … 1.9792 | 1.52 … 1.98 | outward only by luck |
| vσ*²/τ² | exponent | **1.9961** … 2.0399 | **2.00** … 2.04 | lower end inward (D3db-3) |

- Cusp: 3.1340·10⁻⁷, printed 3.2·10⁻⁷ ✓. The Mercer frame's constant is 2⁻⁵², printed 2.3·10⁻¹⁶ ✓. a5's t(v) is 3/7, printed 0.4286 by `up_f` and compared ✓.

**3. Thresholds** (`W/scripts/thresholds.py`, 50-digit Decimal, bisection) [own]
- r* from c² = 4 ln u/(u − 1), u = 1/r:
  - c = 1.995 gives 0.99004569. Its floor is 0.9900; D‴.61's `r_star(1.995)` gives 0.9900456947782 ✓.
  - c = 1.8927 gives 0.80520297; c = 1 gives 0.09664962.
- Closed forms at the two r*:

| r* | Fock σ*/σ < | (τ/σ <) | Mercer-zero τ/σ < | frame end σ*/σ > | half-squeeze σ*/σ < |
|---|---|---|---|---|---|
| 0.9900 | 19.92486 | 19.89975 | 140.71247 | 14.10674 | **199 exactly** |
| 0.8052 | 4.18735 | 4.06619 | 6.51445 | 3.04416 | 9.26694 |

  All are printed inward ✓.
- **198.99 is a float artefact.** In double, (1 + 0.99)/(1 − 0.99) = 198.99999999999983, whose floor is D58's "198.99" (:15, :171). D‴.61 prints "σ*/σ < 199.00 … in exact arithmetic from r's floor, not the float ratio's: referee round 104" (gh_peano.md V; D61.md :270).

**4. A6X** (`W/scripts/a6x_res.py`, `a6x_fit4.py`) [own, raw]
- **Fits by window**, both three-term fits per window (fit 1 = L + c/K² + d/K³, fit 2 = L + b/K + c/K²):

| rule | ρ | fit 1: 80–160, 120–200, 160–240, 200–320 | fit 2: same windows | excess·K², K = 160 → 320 |
|---|---|---|---|---|
| (¼, 2·v_hs) | 16.51 | +8.705·10⁻⁸, +3.517·10⁻¹⁰, +4.917·10⁻¹², +1.547·10⁻¹² | +3.598·10⁻⁷, −1.195·10⁻⁸, −6.664·10⁻⁹, −3.338·10⁻⁹ | −0.09276 → −0.09260 |
| (½, MM) | 1 | +1.929·10⁻¹⁰, +5.752·10⁻¹¹, +2.293·10⁻¹¹, +9.03·10⁻¹² | +1.912·10⁻⁸, +7.438·10⁻⁹, +3.661·10⁻⁹, +1.809·10⁻⁹ | −0.05183 → −0.05192 |
| (1, MM) | 1 | −7.517·10⁻¹¹, −2.247·10⁻¹¹, −8.969·10⁻¹², −3.535·10⁻¹² | −5.890·10⁻⁸, −2.348·10⁻⁸, −1.172·10⁻⁸, −5.850·10⁻⁹ | −0.08650 → −0.08622 |
| (2, MM) | 1 | −1.681·10⁻¹¹, −5.381·10⁻¹², −2.226·10⁻¹², −9.00·10⁻¹³ | −3.521·10⁻⁸, −1.407·10⁻⁸, −7.028·10⁻⁹, −3.512·10⁻⁹ | −0.04709 → −0.04692 |
| (2, 2·v_hs), proved | 0.941 | −6.938·10⁻¹¹, −2.101·10⁻¹¹, −8.443·10⁻¹², −3.345·10⁻¹² | −4.493·10⁻⁸, −1.790·10⁻⁸, −8.929·10⁻⁹, −4.456·10⁻⁹ | −0.05167 → −0.05146 |

- **Agreement with the md.** The md's cells and D58's 18.7–56280.4, 10.0–107.8, 20.7 and 10.1 match. My round-101 table is reproduced digit for digit.
- **The intermediate windows**, which the md does not print, are monotone at all five rules. From 160–240 to 200–320 the factors are 2.47–3.18 on fit 1 (proved 2.52) and 2.00–2.02 on fit 2 (proved 2.00).
- **The last window against the proved rule's bias:** 0.27, 1.06, 2.70 and 0.46 on fit 1, and 0.79, 1.31, 0.41 and 0.75 on fit 2 (MM 2, MM 1, MM ½, outlier).
- **Resolution of the stored doubles.** At 200–320, fit 1 at the outlier is +1.546763·10⁻¹² by float solve, +1.546704·10⁻¹² exact from the doubles, and +1.546441·10⁻¹² from the 25-digit errors. The shrink factor is therefore 56280.4, 56281.9 or 56292.3: three significant digits, not six.
- **Four-term fits through K = 160, 200, 240, 320**, from the 25-digit errors:
  - L + c/K² + d/K³ + e/K⁴: −6.5·10⁻¹³ (outlier), −3.9·10⁻¹⁴, +8.6·10⁻¹⁵, −3.6·10⁻¹⁴ (MM), −2.0·10⁻¹⁴ (proved);
  - with a K⁻¹ term: −1.1·10⁻¹¹, −4.3·10⁻¹¹, +1.7·10⁻¹¹, +4.1·10⁻¹², +1.6·10⁻¹¹.
- **Planted offsets.** I scaled the ratios by 1 + δ and ran the repository's `a6x_lines` on the result.
  - **At the outlier alone**, δ = 10⁻⁹, 10⁻⁸ and 5·10⁻⁸ all still print "every deviation shrinks as the window moves out". At 5·10⁻⁸ the factors are 2.7–21.4 and 8.8–10.6.
  - At every rule, δ = 10⁻¹⁰ prints "NOT …", because the moment-matching rules' small deviations are swamped.
- **Coverage.**
  - A6X has 4 of the 9 fast-branch rules without a frame: 3 moment-matching rules (ρ = 1) and 1 of the 5 with ρ > 1.
  - Not extended: (¼, 1.25·v_hs) ρ 6.61, (½, 1.25·v_hs) 1.85, (½, 2·v_hs) 4.53, (1, 2·v_hs) 1.59, and moment matching at ¼ (near).
  - Calibrators: 1 of the 4 proved fast-far rules. (1, 1.25·v_hs), (2, Mercer frame) and (2, 1.25·v_hs) are not extended.

**5. Part V's weight** [md]
- **D‴.61 at d5d8438** (16:42), which D58 at ce7cfac (17:14) cites:
  - (iii)(b): "certified in ball arithmetic (V) … computer-assisted, and rigorous";
  - (iv): "for r < 0.9900 with V's certificate … (modulo the branch and bound; V rigorous in ball arithmetic)";
  - Status: "a second computer-assisted step, and a rigorous one".
- **D58's tier-2 text was written at 04bca94** (15:59), when V was still floating point (b98e652, 15:18), and was not revised.
- **Round 108 (D3dc)** found "'Certified in ball arithmetic' and 'rigorous' are earned for V. 'Modulo the branch and bound alone' is accurate."
- **At HEAD (88ef79f)** D‴.61 reads "both computer-assisted steps in ball arithmetic".

**6. Mutants** (W/mut/spec107.py: R01–R24 plus control R00, aimed at A6X's code, which round 101's spec does not touch, and at part V's labels) [own]
- **Pass 1.** R00 passes. 18 are killed: the 16 md-changers, and R14 and R15 by unit tests without an md change.
  - 6 survive, all with the md unchanged:
    - R02, R03, R04 and R05 weaken the A6X verdict: and→or; fit 1 only; threshold 1 → 10⁻³, so growth up to 1000× would pass; the unproved rules only. All are equivalent on the data and unsound.
    - R19 (distinct by v at 3 decimals) and R21 (`out_e`'s guard 10⁻⁶) are equivalent on the data.
- **Pass 2**, the 16 md-changers with `-k 'not round88_precision_and_the_md'`:
  - the unit tests kill 7 (R01, R06, R11, R12, R16, R22, R23);
  - 9 are caught only by the md regression: R07, R08, R09, R10, R13, R17, R18, R20, R24.
  - Of those 9, R10 and R20 flip the published verdict to "NOT every deviation shrinks", R13 drops the A6X line from the report, and R24 shifts every deviation by 10⁻¹³ relative.
- **Round 101's spec as shipped:** X00 passes; 19/19 are killed, N02, N12 and N13 by the tests without an md change. This matches the log.

**7. Literature.** Six searches; the snippets are listed under Sources. D58's literature passages are unchanged since 0cb96bd.
- The snippets agree with D58's descriptions of Karvonen–Särkkä, Derevianko–Plonka and Hu–Yu.
- For Fasshauer–McCourt's Mercer/Hermite expansion the snippets show the method but not the β formula, so that formula stays unconfirmed.
- If D58 adopts "ball arithmetic" (D3db-1), the source for Arb is Johansson, arXiv:1611.02831 (ARITH 24, 2017) [snippet].

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| (i) :14–15: "proved by D‴.61 (iv) for t < 0.8052 (σ*/σ < 9.26; modulo its branch and bound, and without it for t < 0.0966) and for t < 0.9900 (σ*/σ < 198.99; modulo also its valley certificate, D‴.61 V)" | The numbers are right except that 198.99 should be 199. The weight is wrong: V is rigorous, so tier 2 carries tier 1's proviso alone | [own] thresholds; [md] D61 at d5d8438; D3db-1 |
| (ii) :24–26: "and with t(v) < 0.9900 modulo also its valley certificate" | Same weight error | D3db-1 |
| Status (i) :170–173: "modulo both of D‴.61's computer-assisted steps … both in floating point with a slack, not interval arithmetic" | **False** for V at ce7cfac, and for both at HEAD | [md] D61 (iii)(b), Status; gh_peano.md V |
| Status (ii) :182–185: tier-2 thresholds 19.924 and 140.712; "conjectural where a rule has no frame or t(v) ≥ r*" | Thresholds right and inward. "r*" names no tier | [own]; D3db-1 |
| Status, A6 paragraph :246–253: outward ranges; "8 (7 distinct)"; "9 rules without a frame (8 far, MM ¼ near)"; signs; "between their proved neighbours' at 4 of the 5"; "no proved rule at τ/σ = 1/4" | All right. The last is structural (no far rule at ¼ can have a frame) | [own] §2 |
| :253–257: outlier "41.4 and 2.0 times"; A6X "+1.55·10⁻¹² and −3.34·10⁻⁹, … −0.0926 … every deviation of the four fast-branch rules without a frame shrinks …, by 18.7–56280.4 and 10.0–107.8 times, as the proved … by 20.7 and 10.1" | Numbers right (56280.4 to ~3 digits). "The four" reads as all, though they are 4 of 9. The verdict is necessary, not a calibration. The outlier's 56280/108 is not "as the proved's" until the later windows | [own] §4; D3db-2, -3, -4 |
| :267–269: "on \|η\| ≤ 1.995 modulo also a floating-point valley certificate, its V" | **False**: V is certified in ball arithmetic | [md] D61; D3db-1 |
| :270–271, :294: "(t(v) < 0.9900 with the valley certificate)" | Neutral wording, consistent with D‴.61 | — |
| New :289–291: "… σ*/σ < 19.924 and τ/σ < 140.712 modulo also its valley certificate" | Weight wrong | D3db-1 |
| Numerics (A6) :320–326 | Ranges outward, except the exponent's "2.00". A6X not mentioned | D3db-3, -4 |
| md headings (A6 :128, A6E :185, D :514) and the threshold line :179: "modulo also its valley certificate" | Numbers right; weight wrong; pinned by tests :1780, :1789 | D3db-1, -5 |
| md A6X line :183 | Numbers right; "20.7–20.7 … at the 1 proved"; the verdict criterion weak | D3db-2, -4 |
| (iii)–(v), their proofs, and B, C, D | Unchanged since 0cb96bd; the tests pass | [md] |

What is proved, and how (at ce7cfac):
- **Proved:**
  - (i)'s identity, its bounds and the rate t;
  - (v);
  - the order of λ, λ/(2 − λ) and q against t.
- **Proved via D‴.61, modulo its floating-point branch and bound only:**
  - the sharp constant and every frame rule's rate for t(v) < 0.9900;
  - unconditionally for t < 0.0966.
  - Since 88ef79f D‴.61 calls that step rigorous too.
- **Numerical:**
  - moment matching's λ;
  - the 14 rules without a frame: slow branch by bracket and t(v) interleaving; fast branch by A6X's later windows and four-term fits, which the md does not print;
  - (iii) and (iv)'s numbers.
- **Conjectural:** t(v) at rules without a frame or with t(v) ≥ 0.99, and the sharp constant beyond.

## Status of the round-101 issues

| D3cv | status | note |
|---|---|---|
| 1 | fixed, residues | The md prints the outlier's ratios. Part A6X reproduces the extension exactly, and D58 cites it. Residues: the verdict criterion, the "four" count, one calibrator (D3db-2) |
| 2 | fixed | The slow list by t(v) with frame marks and both fits; "(8, 7 distinct)"; "9 (8 far, MM near)". All verified. The overall counts (22 of 36, 18, 32) still include the duplicate (D3db-4) |
| 3 | fixed where named | `out_e` on the line's three value ranges, 2.3·10⁻¹⁶, 0.4286 via `up_f`. The same line's exponent range was missed, and 198.99 was added (D3db-3) |
| 4 | fixed, new staleness | "or t(v) ≥ r*" (but which r*); the cusp marked proved; HANDOFF row 34 updated. That row now says "modulo two floating-point computer-assisted steps", stale since d5d8438 (D3db-1) |
| 5 | fixed | a5's "<" is computed and tested; the planted cells' maxima and brackets are asserted; 19/19 of round 101's mutants are killed (rerun). A6X's code was unmutated until this round (D3db-5) |
| 6 | not applied | Logged as needing the md to print D‴.61 (iv)'s bound first. That reason is thin: D61 (iv) states a/(2(1 + ρ)), and D58 cites D61 (iv)'s text elsewhere. It remains optional |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3db-1 | minor | **Part V's weight is understated, and two statements are false.**<br>• Since d5d8438, 32 minutes before ce7cfac, D‴.61 reads "V rigorous in ball arithmetic" and "a second computer-assisted step, and a rigorous one", with the branch and bound the only floating-point step. Round 108 found "'modulo the branch and bound alone' is accurate". At HEAD (88ef79f) both steps are rigorous [md].<br>• D58 Status (i) :171–173 says "both in floating point with a slack, not interval arithmetic"; :268 says "a floating-point valley certificate". Both are false.<br>• "Modulo also its valley certificate" appears at :15, :25 and :290, and in md :128, :179, :185 and :514. It is pinned by tests :1780 and :1789.<br>• So at ce7cfac tier 2 rests on exactly tier 1's proviso.<br>• "σ*/σ < 198.99" (:15, :171) is the double ratio's floor (198.99999999999983); D‴.61 gives 199.00 from r's floor in exact arithmetic [own].<br>• Status (ii)'s "t(v) ≥ r*" names no tier. HANDOFF row 34's "two floating-point … steps" is stale. The tier-2 text was written at 04bca94, before V's ball recheck. | **D58 (i), (ii), Status, sketch, New**, following D‴.61's current Status. At ce7cfac: "proved by D‴.61 (iv) for t < 0.9900 (σ*/σ < 199), without a computer-assisted step for t < 0.0966, beyond it modulo D‴.61's floating-point branch and bound on \|η\| ≤ 1.8927 (its valley certificate on 1.8927 ≤ \|η\| ≤ 1.995 rigorous in ball arithmetic)". At HEAD: "both computer-assisted steps rigorous in ball arithmetic (Arb)".<br>**Code:** derive the labels from one status read from gh_peano's raw (`V['ball']`, D's ball record), not fixed text.<br>**Tests:** pin `D61_RSTAR_V == floor(GP.r_star(1.995)·10⁴)/10⁴`; update the pinned strings; cite Arb |
| D3db-2 | minor | **A6X's verdict does not calibrate the fast branch.**<br>• The verdict is "each shrink factor > 1 between two windows". With the outlier's limit planted 10⁻⁹, 10⁻⁸ or 5·10⁻⁸ off 1/t(v), the repository's `a6x_lines` still prints "every deviation shrinks as the window moves out" [own].<br>• Mutants R02–R05 (and→or, fit 1 only, threshold 10⁻³, unproved only) survive.<br>• Coverage: 4 of the 9 fast-branch rules without a frame (3 MM, 1 of the 5 with ρ > 1). D58 :256's "the four fast-branch rules without a frame" reads as all of them.<br>• One proved calibrator, at τ/σ = 2. At τ/σ = ¼ no far rule can have a frame.<br>**What does calibrate, already in the raw** [own]:<br>• the last window against the proved rule's bias: 0.27–2.70 (fit 1) and 0.41–1.31 (fit 2);<br>• the later windows (160–240 → 200–320) shrinking as the proved rule's: 2.47–3.18 vs 2.52, and 2.00–2.02 vs 2.00;<br>• four-term fits through K = 160–320 from the 25-digit errors: within 6.5·10⁻¹³ (proved 2.0·10⁻¹⁴), and with a K⁻¹ term within 4.1·10⁻¹²–4.3·10⁻¹¹ (proved 1.6·10⁻¹¹).<br>• No K beyond 320 is needed for this, and the precision is adequate (+40 digits identical). | **`a6x_lines`:** print the intermediate windows, the last window over the proved rule's, and a four-term fit computed from `e0`/`e1` (not the doubles). Base the verdict on "within X times the proved rule's bias".<br>**D58 :253–257:** "the four of the nine extended (the other four, ρ = 1.59–6.61, not)"; "as the proved rule's in the later windows".<br>**Optionally,** extend A6X with (1, 1.25·v_hs), (½, 1.1·v_hs) (a far fast-branch frame rule at ½) and the four unextended far rules (≈45 min at τ/σ = ¼, ≈22 min each at ½, ≈6 min at 1).<br>**A test** with a planted one-rule offset |
| D3db-3 | cosmetic | **Rounding residues.**<br>• The exponent range "2.00 … 2.04" (md :173; D58 :243) has true lower end 1.9961: rounded inward. "1.52 … 1.98" is outward only by luck. `rng` with `:.2f` was missed by the sweep.<br>• A6X's ranges are rounded to nearest: 18.687 prints as "18.7", a lower end rounded up.<br>• "56280.4" has about 3 significant digits: 56281.9 exact from the doubles, 56292.3 from the 25-digit errors.<br>• 198.99 (see D3db-1). | Use `out_e` or a floor/ceil pair for the exponent and the shrink ranges. Print "≈5.6·10⁴", or fit from the 25-digit errors. D58 follows |
| D3db-4 | cosmetic | **Wording, counts, staleness.**<br>• md :183 has "20.7–20.7 and 10.1–10.1 at the 1 proved" and "at the 1 rules". "Proved" comes from the frame test alone, without checking t(v) < r*.<br>• D58's "22 of 36", "18 frame rules" and "32 rules" count Fock ≡ Mercer-zero at τ/σ = 2 twice (21 of 35; 17 distinct).<br>• The Numerics do not mention A6X.<br>• "As the proved (τ/σ = 2, 2·v_hs)'s does" is an awkward possessive.<br>• HANDOFF :183–185 at ce7cfac still says to run A6X "if its raw did not reach the repository", though it did.<br>• The log's Round 101 entry does not note that V became rigorous mid-round. | Singular/plural forms and "proved" requiring t(v) < r*; "(21 of 35 distinct)"; an A6X clause in the Numerics; update HANDOFF |
| D3db-5 | cosmetic | **Tests.**<br>• 9 of 16 md-changing mutants are caught only by the md regression, among them R10 and R20 (the verdict flipped to "NOT") and R13 (the A6X line dropped from the report).<br>• `test_round101_a6x_…` asserts strings of the committed md file, not of `report()`'s output.<br>• `D61_RSTAR_V` is unpinned, while `D61_RSTAR` is pinned to `GP.r_star(1.8927)`.<br>• R19 and R21 survive, equivalent on the data. | Assert on `WM.report(*raw parts)`. Add unit tests for the verdict (D3db-2), the fit forms and the windows. Pin `D61_RSTAR_V`. Record R19 and R21 as equivalent |

## Sources

- [Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic, arXiv:1611.02831](https://arxiv.org/pdf/1611.02831) [snippet]. Johansson; ball arithmetic; presented at ARITH 24 (2017); Arb merged into FLINT in 2023.
- [Derevianko–Plonka, arXiv:2307.13587](https://arxiv.org/pdf/2307.13587) [snippet]. Cosine sums with nodes at the zeros of a scaled Hermite polynomial; exponential decay in a weighted L², and in unweighted L² after truncation. A citing reference gives Sampl. Theory Signal Process. Data Anal. 24(6), 2026.
- [Hu–Yu, arXiv:2412.08044](https://arxiv.org/html/2412.08044v3) [snippet]. The optimal scaling balances the spatial and the frequency truncation errors; one source says "to appear in SIAM J. Numer. Anal. 2026".
- [Karvonen–Särkkä, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532), [BIT 59(4) 877–902, doi:10.1007/s10543-019-00758-3](https://link.springer.com/doi/10.1007/s10543-019-00758-3) [snippet]. Kernel quadrature at scaled Gauss–Hermite nodes through a truncated Mercer expansion.
- [Yurova–Kormann, arXiv:1709.02164](https://arxiv.org/pdf/1709.02164) [snippet]. Describes Fasshauer–McCourt's Mercer/Hermite (Gauss-QR) approach; the β formula is not shown in the snippets.
- Repository sources [md]:
  - D61.md at ce7cfac and at HEAD;
  - gh_peano.md, parts D and V;
  - REVIEW_round108_D3dc.md;
  - commits 04bca94, b98e652, d5d8438, ce7cfac and 88ef79f, for times and contents.
- [memory]: a fixed-exponent extrapolation's bias shows as drift with the window. A shrinking deviation bounds a limit offset only by the last window's size, not by the shrink factor.

## Lesson

**When a cited result's status changes during a round, re-derive the citing text's weight from the source's current Status, not from the draft written before the change.** Better still, read it from a flag in the source's raw. Here D‴.61 made V rigorous 32 minutes before D‴.58's commit, and D58 still calls V floating point. The tests pin the stale wording, so they would have failed on the correct wording.

**A shrinking deviation is necessary, not a calibration.** It passes for any limit offset smaller than the first window's deviation. Here an offset of 5·10⁻⁸ at the outlier still passes. What calibrates is the last window's size against a proved rule's bias in the same window, and a planted offset is the cheap test of whether a verdict line can tell the two apart.
