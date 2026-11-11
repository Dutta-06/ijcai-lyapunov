# Referee report D3cq, round 96: Proposition D‴.58 as revised after round 92

## Scope

- **Subject.** Commit 420a3ef on claude/vigilant-clarke-hf939b (round 92 applied), diffed against 363855d. I exported it with `git archive 420a3ef` into W/repo, with W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref96`.
  - WM_OUT was W/wmout and GP_OUT was W/gpout, each holding a copy of its raw.
  - Scripts are in W/scripts, logs in W/logs, mutant work in W/mut, W/mut2 and W/mut92.
- **Read in full.**
  - D58.md, `exp_width_mismatch.py` (every hunk changed since 363855d, line by line), `width_mismatch.md` and `test_width_mismatch.py` (56 tests).
  - REVIEW_round92_D3cm.md, the Round 92 entry of REVIEW_LOG.md (and Round 93's), `tools/mutants/round92_D58_own.py` and `tools/mutants.py`.
  - `git diff 363855d 420a3ef` for these files, plus a word diff of D58.
  - For the citations: D61.md as revised in round 93 (statement, proof of (iii)(b)–(iv), Status, New, Numerics), gh_peano.md part D and `exp_gh_peano.py`'s `c_need`, `r_star` and `d_lines`.
- **Run settings.**
  - Python 3.13.18, numpy 2.5.4, scipy 1.18.1, from the repository's venv.
  - Every process ran with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, at most two of mine at a time.
  - pytest always ran with `--basetemp` under W, including inside my copy of the mutant runner (W/scripts/mutants96.py; it differs from the repository's only in REPO, PY and `--basetemp`).
- **Nothing in the repository was modified.** I ran only `git archive/diff/log/show/status/rev-parse`. All my processes have exited; I stopped none.
  - Meanwhile another session moved HEAD to cace1fd (D‴.61 round 95 and others). Since 420a3ef no D58 file has changed. D61.md and gh_peano.md have changed: 49.92 → 49.89, and "|η| ≤ 1.98" → "the annulus 1.8927 < |η| ≤ 1.98". I review 420a3ef and flag the follow-up in D3cq-5.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my own computations (my own code; for the A6 extension, the repository's `a6_error` at larger K with my own fits); [snippet] web-search results only (proceedings.nips.cc and arxiv.org do not resolve from this session); [memory] recalled facts.

## Summary

**Design.**
- Reproduced:
  - both mds byte for byte;
  - the tests;
  - six jobs bit-identically from the raw.
- Recomputed independently:
  - A6's two fits, from the raw's ratios, with exact rational solves, and ρ, t(v) and the argument in 50-digit Decimal from τ and the scalings;
  - the split by frame, by branch and by nearness to the cusp;
  - D‴.61's thresholds r*(c), and the Fock, Mercer-zero and frame-end thresholds, by closed form and by root finding.
- Extended A6 beyond K = 160 at three rules, to test "their bias":
  - a proved rule, (2, ½·v_hs);
  - the unproved rule that sets 7.5·10⁻⁶ and 4.4·10⁻⁵, (1, ½·v_hs);
  - the unproved fast-branch outlier, (¼, 2·v_hs).
- Mutants:
  - 24 single-point mutants plus a control of this round's code;
  - a second pass of the md-changing ones without the md-regression test;
  - a rerun of round 92's own spec.
- Literature: six web searches.

**Verdicts.**
- **Reproducibility is exact** [own]:
  - `--from-raw` gives md5 56ac0748… = committed (raw 959dc055…, unchanged since round 88); gh_peano.md is identical too;
  - 56/56 tests pass in 54.9 s;
  - all six rerun jobs are identical in every field but `secs`.
- **Every number of the round traces to the md and the raw** [own]:
  - frame 18 / no frame 14;
  - 6.0734·10⁻⁵, 1.3900·10⁻⁴ and 1.9974·10⁻⁴ for the frame rules;
  - 7.547·10⁻⁶ and 4.449·10⁻⁵ for the rules without a frame;
  - brackets 10 of 18 and 9 of 14.
- **D58's claims about D‴.61 match D‴.61 at 420a3ef** [own]:
  - r* = 0.8052 (from 0.805203); σ*/σ < 9.26 (from 9.2671);
  - 0.0966 (from 0.096650);
  - 49.92 (from 49.922 at c = 1.98);
  - Fock σ*/σ < 4.187 (4.18735; τ/σ < 4.06619);
  - Mercer-zero τ/σ < 6.514 (6.51445);
  - "modulo its branch and bound" is D‴.61's proviso.
  - Only "> 3.044" is rounded the wrong way (true value 3.04417).
- **"Their bias" is right.** Extended to K = 320 at the proved rule (2, ½·v_hs), the two fits' deviations shrink toward the proved limit from both sides:
  - first fit: +6.07 → +3.22 → +1.93 → +1.10·10⁻⁵;
  - second fit: −1.39·10⁻⁴ → −8.76 → −6.00 → −4.00·10⁻⁵.
- **The proofs are unchanged and right**, and no claim overreaches in substance. No major issue.

**Findings.**
1. **D3cq-1 (minor): the re-partition stops at the A6 paragraph.** Status (ii) still cites the pooled maxima as the evidence for the conjectured curve: 6.1·10⁻⁵ and 1.4·10⁻⁴ over 32 rules, and A5's 2.3·10⁻⁵.
   - All three are set at rules whose limits D‴.61 now proves.
   - The "conjecture's values" (b² + b + 2)/(b(b − 1)) are now proved values.
2. **D3cq-2 (minor): the split pools two branches, so the new counts carry little.**
   - At the proved rules the fits bracket 8 of 8 on the slow branch (far from the cusp) and 0 of 4 on the fast one. Bracketing is therefore calibrated only on the slow branch, where the rules without a frame bracket 5 of 5: that is the evidence.
   - "9 of 14" adds 4 uninformative fast-branch brackets.
   - My extension makes the evidence much stronger:
     - at (1, ½·v_hs), K ≤ 320, the bracket narrows to [−8.3·10⁻⁶, +8.2·10⁻⁷];
     - at (¼, 2·v_hs) the +8.7·10⁻⁸ outlier is pre-asymptotic: +4.9·10⁻¹² through K = 160–240, bracketed.
3. **D3cq-3 (minor): hypotheses still uneven.**
   - Status (ii)'s Fock and Mercer-zero clause, :240 and the md's split label lack "modulo its branch and bound".
   - New's "the rate t, strictly faster than the Fock and Mercer-zero scalings'" is unqualified. Yet at the Fock rule D‴.61's unconditional bound a/2 never exceeds t: (a/2)/t = (s + 1)²/(2(s² + 1)) ≤ 1, s = σ*/σ. So the claim needs D‴.61's range.
4. **D3cq-4 (cosmetic): bounds rounded down.**
   - "within 7.5·10⁻⁶ and 4.4·10⁻⁵" (7.547, 4.449);
   - "within 3.1·10⁻⁷" at the cusp (3.134);
   - "once σ*/σ > 3.044" (3.04417);
   - HANDOFF has 7.6 and 4.5.
5. **D3cq-5 (cosmetic): stale labels.**
   - md :183 and code :1654–1655 still say "t(v) < 0.8053, modulo a grid-checked domination condition".
   - "the conjecture's" values at D58 :225–226 and :296 and md :123.
   - "checked on a grid, to hold on |η| ≤ 1.98" for what D‴.61 checked at 57 points on three circles.
   - D‴.61 has since moved to 49.89 and "the annulus".
6. **D3cq-6 (minor): the four QFF citations are described beyond what the sources show.**
   - TQFF (arXiv:2310.14544) uses a non-Gaussian trigonometric rule.
   - arXiv:2609.00553 (Cloninger–Le Gia–Mhaskar) is about eigenvalue bounds; Gauss–Hermite is not confirmed.
   - The authors are not named.
7. **D3cq-7 (cosmetic): mutants.**
   - 7 of 24 survive, all equivalent on the data; none changes a published result.
   - But 8 md-changing mutants are caught only by the md-as-report regression test. The split's maxima and bracket counts and the Mercer-zero threshold have no unit test of their semantics.

## What I checked and how

**Reproduction.**
- `--from-raw` into W/wmout, then `cmp` with the committed md: identical. The same for gh_peano (GP_OUT = W/gpout): identical (raw 199e6201…).
- Full test file: 56 passed in 54.90 s.
- `W/scripts/rerun.py` reruns a job and compares it with `==` in every field but `secs` [raw, own]. All six are identical:

| job | why this one | time |
|---|---|---|
| A6 (1, ½·v_hs), K ≤ 160 | the rule without a frame that sets 7.5·10⁻⁶ / 4.4·10⁻⁵ | 19.7 s |
| A6 (2, ½·v_hs) | the proved rule that sets 6.1·10⁻⁵ / 1.4·10⁻⁴ | 6.9 s |
| A6 (¼, Mercer frame) | a near-cusp frame rule, −4.0·10⁻⁵ | 139 s |
| A6 (½, moment matching) | a rule without a frame | 59.8 s |
| A6E K = 20, imaginary axis | — | — |
| B5 job 0 ((1, 12), θ = 0, direct) | — | 2.1 s |

**A6 split (`W/scripts/a6split.py`).**
- Method:
  - L + c/K² + d/K³ and L + b/K + c/K² through K = 80, 120, 160, solved exactly in Fractions from the stored doubles;
  - ρ = 1 + 2/τ² − 2/v, t(v) and the attaining argument in Decimal from τ and `a6_scalings`' definitions, not from the stored v.
- Checks:
  - stored rates agree to 1.6·10⁻¹⁵ and every stored argument agrees;
  - 22 rules have a frame; ρ = 1 exactly at the 4 moment-matching rules; the largest frame |ρ| is 0.94113; the largest frame t(v) is 0.634512.
- Results:

| set | n | max \|fit 1\| | max \|fit 2\| | max spread | bracket |
|---|---|---|---|---|---|
| frame, off cusp | 18 | 6.0734·10⁻⁵ | 1.3900·10⁻⁴ | 1.9974·10⁻⁴ | 10 |
| — slow (τ²/v), far | 8 | 6.07·10⁻⁵ | 1.39·10⁻⁴ | | 8 |
| — slow, near | 3 | 2.91·10⁻⁵ | 5.34·10⁻⁶ | | 2 ((½, Fock) not) |
| — fast (vσ*²/τ²), far | 4 | 2.10·10⁻⁹ | 1.73·10⁻⁷ | | 0 |
| — fast, near | 3 | 3.99·10⁻⁵ | 3.67·10⁻⁵ | | 0 |
| no frame, off cusp | 14 | 7.5472·10⁻⁶ | 4.4490·10⁻⁵ | 5.2037·10⁻⁵ | 9 |
| — slow, far | 5 | 7.55·10⁻⁶ | 4.45·10⁻⁵ | | 5 |
| — fast, far | 8 | 8.71·10⁻⁸ ((¼, 2·v_hs)) | 3.60·10⁻⁷ | | 4 |
| — fast, near | 1 (MM, ¼) | 5.3·10⁻¹¹ | 3.8·10⁻⁸ | | 0 |

- At the cusp the largest deviation is 3.134·10⁻⁷.
- The deviations on the slow branch grow with t(v), the same way with and without a frame:
  - (1, ½·v_hs), t(v) = 0.4776, no frame: +7.55·10⁻⁶ / −4.45·10⁻⁵;
  - (2, 0.8·v_hs), t(v) = 0.4730, frame: +1.10·10⁻⁵ / −5.21·10⁻⁵;
  - (2, Fock), t(v) = 0.5, frame: +1.42·10⁻⁵ / −6.15·10⁻⁵.

**A6 extended (`W/scripts/a6extend.py`)** [own]. The repository's `a6_error` in its own precision rule, with my fits through consecutive triples:

| rule | fit through | fit 1 | fit 2 |
|---|---|---|---|
| (2, ½·v_hs), proved | 80–160 | +6.073·10⁻⁵ | −1.390·10⁻⁴ |
| | 120–200 | +3.222·10⁻⁵ | −8.759·10⁻⁵ |
| | 160–240 | +1.928·10⁻⁵ | −6.005·10⁻⁵ |
| | 200–240–320 | +1.103·10⁻⁵ | −3.998·10⁻⁵ |
| (1, ½·v_hs), no frame | 80–160 | +7.547·10⁻⁶ | −4.449·10⁻⁵ |
| | 120–200 | +3.279·10⁻⁶ | −2.351·10⁻⁵ |
| | 160–240 | +1.678·10⁻⁶ | −1.412·10⁻⁵ |
| | 200–240–320 | +8.198·10⁻⁷ | −8.258·10⁻⁶ |
| (¼, 2·v_hs), no frame | 80–160 | +8.705·10⁻⁸ | +3.598·10⁻⁷ |
| | 120–200 | +3.517·10⁻¹⁰ | −1.195·10⁻⁸ |
| | 160–240 | +4.917·10⁻¹² | −6.664·10⁻⁹ |

- At (¼, 2·v_hs) the ratio's excess is −3.62, −2.32, −1.61·10⁻⁶ at K = 160, 200, 240.
- At K = 80, 120, 160 these reproduce the md's deviations exactly.

**Thresholds (`W/scripts/thresholds.py`; 50-digit Decimal, my own bisection and closed forms)** [own].
- r*(c), from c² = max(4 ln u/(u − 1), 4r), u = 1/r:

| c | r* | σ*/σ at the half-squeeze rule |
|---|---|---|
| 1.8927 | 0.8052030 | 9.26710 |
| 1.892733 | 0.8052572 | — |
| 1.98 | 0.9607245 | 49.9224 |
| 1 | 0.0966496 | 1.21398 |

- c_need(0.634512) = 1.777342; c_need(0.8052) = 1.892698 < 1.8927.
- Fock: λ* = 2r*/(1 + r*) = 0.8920895, σ*/σ < 4.187348, τ/σ < 4.066188. The same by bisection on t(v) at v = 2τ²/(2 + τ²).
- Mercer-zero: τ/σ < ((((1 + r*)/(1 − r*))² − 1)/2)^{1/2} = 6.514453. The same by t(v).
- Frame ends: λ > r* iff σ*/σ > ((1 + r*)/(1 − r*))^{1/2} = 3.044165 (3.044191 at r* = 0.805203).
- Fock and Mercer-zero lie inside the frame interval at τ = 0.1 … 100.
- Algebraically: ρ_Fock = 0; ρ_MZ = 1 − 2(√(1 + 2x) − 1)/x ∈ (−1, 1); Fock's t(v) = tanh((r + ln cosh r)/2) = λ/(2 − λ).

**D‴.61 cross-check.** D58's citations against D61.md at 420a3ef: statement (iii)(b) and (iv), the Status, New and Numerics, and gh_peano.md :135–143. They agree:
- r* = 0.8052 and "σ*/σ < 9.26" (d_lines: floor of (1 + 0.8052)/(1 − 0.8052));
- "for r < 0.0966 without it";
- "r < 0.9607 (σ*/σ < 49.92) if the valley condition holds on |η| ≤ 1.98";
- "floating point with slack 10⁻¹², not interval arithmetic";
- the valley condition at 57 points on |η| = 1.9, 1.95, 1.98.

**Unconditional lower bound at the Fock rule** [own]. ρ = 0 and a = τ²/(2 + τ²) (σ = 1). With s = σ*/σ, (a/2)/t = (s + 1)²/(2(s² + 1)) ≤ 1, with equality only at s = 1. D‴.61 (iv)'s unconditional a/2 therefore never shows the Fock rule slower than t; for example at τ/σ = 1, 0.1667 against 0.1716.

**Mutation testing.**
- **First pass**: `W/mut/spec96.py`, 24 mutants plus the control Y00, run with the full suite, `-x` and `--basetemp`. The mutants cover:
  - `fr_`: G01–G03, G24;
  - the split's lists, maxima and counts: G04–G08;
  - `cross` and the thresholds: G09–G15, G23;
  - D61_RSTAR and the literals: G16, G17, G20;
  - the spread: G18;
  - the farthest-rule key: G19;
  - the B4 and A6E labels: G21, G22.
- **Second pass**: the 15 md-changing mutants with `-k 'not round88_precision_and_the_md'`.
- **Round 92's spec**: `tools/mutants/round92_D58_own.py`, rerun as is.

**Literature.** Six searches (Sources). proceedings.nips.cc gave ENOTFOUND; everything is [snippet].

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| (i) "proved by D‴.61 (iv) for t < 0.8052 (σ*/σ < 9.26; modulo its branch and bound, and without it for t < 0.0966)" (:14–15) | Correct; rounded down safely | [own] 0.805203, 9.2671, 0.096650 |
| (ii) "(A6; proved by D‴.61 (iv), modulo its branch and bound, for every v ∈ (τ²σ²/σ*², τ²) … with t(v) < 0.8052, conjectural for the others)" (:24–25) | Correct (D‴.61 (iv) at 420a3ef) | [own] |
| Status (i): "conjectural beyond (for σ*/σ < 49.92 were D‴.61's valley condition, checked on a grid, to hold on \|η\| ≤ 1.98)" (:171–172) | Number right (49.922). "Checked on a grid" means 57 points on three circles; D‴.61 is now "the annulus", 49.89 | [own]; D3cq-5 |
| Status (ii): Fock and Mercer-zero "proved by D‴.61 (iv) while below its r* (… σ*/σ < 4.187 and τ/σ < 6.514)" (:181–182) | Numbers right. Proviso missing | [own]; D3cq-3 |
| Status (ii): conjectured curve's evidence "within 6.1·10⁻⁵ … 1.4·10⁻⁴ … A5 … 2.3·10⁻⁵" (:183–186) | Numbers right. All three are set at proved rules (not re-partitioned) | [own]; D3cq-1 |
| Status (ii): "made a proof … for every v with a frame and t(v) < 0.8052, A6's 22 frame rules among them, modulo …; … 10 rules with \|ρ\| > 1 … A6's numbers alone" (:186–188) | Correct | [own] |
| A6 paragraph: 13 / 12 rules, the ranges, 13 of 13, K^{−1.52…−1.98} (:235–240) | Unchanged; right | [md] |
| :240–242: "at the 18 frame rules … miss them by up to 6.1·10⁻⁵ and 1.4·10⁻⁴ — their bias; … at most 2.0·10⁻⁴ —, and at the 14 without a frame … within 7.5·10⁻⁶ and 4.4·10⁻⁵" | "Their bias" confirmed by extension. 7.5 and 4.4 rounded down. Proviso missing. Pooled across branches | [own]; D3cq-2, -3, -4 |
| md split line (bracketing at 10 of 18 and 9 of 14) | Counts right; uninformative without the branch | [own]; D3cq-2 |
| md :179: Fock / Mercer-zero thresholds; "> 3.044" | First three right; 3.044 should be 3.045 | [own]; D3cq-4 |
| :252–255 D‴.61 (iii) "on both axes' intervals \|η\| < 2 and on \|η\| ≤ 1.8927 (modulo a floating-point branch and bound …)" … "all 22 …" | Correct | D61 at 420a3ef |
| :225–226, :296 "the conjecture's values (b² + b + 2)/(b(b − 1))" | Now D‴.61-proved (Mercer-frame rule: ρ = q, t(v) ≤ 0.4286; unconditional at τ/σ = ¼, ½, t(v) = 0.0154 and 0.0583 < 0.0966) | [own]; D3cq-5 |
| New: "every scaled rule's rate — proved with D‴.61 for the rules with a frame and t(v) < 0.8052, modulo its branch and bound" (:272–273) | Fixed; correct | — |
| New: "the rate t, strictly faster than the Fock and Mercer-zero scalings'" (:269–270) | Unqualified; holds only in D‴.61's range (Fock's unconditional a/2 ≤ t) | [own]; D3cq-3 |
| Classical: QFF papers "[snippets, referee round 92]" (:257–259) | Cited and "unchecked" as asked; two of four misdescribed | [snippet]; D3cq-6 |
| (iii)–(v), their proofs; B, C, D numerics | Unchanged this round. B4's label fixed. Spot checks against the md pass (93.0, 107, 316, 0.963, 3643, 31769, 1.0235, 1.0484, 0.0972, 0.954, 5.77, 65.054019) | [md] |

## Status of the round-92 issues

| D3cm | status | note |
|---|---|---|
| 1 | fixed, with residues | Done: "away from the cusp (\|ln(v/v_hs)\| ≥ 0.05)"; "the second fit's largest deviation; the two fits differ by at most 2.0e-04" (1.9974·10⁻⁴ ✓); the split by frame in the md and at :240–242. Residues: the Status still quotes the pooled maxima (D3cq-1). The split pools branches and drops round 92's "bracketed 5 of 5" in favour of 9 of 14 (D3cq-2). 7.5 / 4.4 are rounded down from round 92's 7.55 / 4.45, where round 92 gave 7.6 / 4.5, as HANDOFF does (D3cq-4). |
| 2 | fixed, with residues | Done: (a) the proviso in (i), (ii), the Status's general clause and New; (b) Fock and Mercer-zero proved below r* (4.187, 6.514 ✓), moment matching numerical; (c) "every v with a frame and t(v) < 0.8052, A6's 22 among them"; (d) 0.8052 and 9.26, D61_RSTAR pinned by a test; (e) 0.0966 in (i), 49.92 added. Residues: the Fock/Mercer-zero clause, :240, the md's split label and New's "strictly faster" (D3cq-3); "> 3.044" (D3cq-4). |
| 3 | confirmed | Round 92's own spec rerun: X00 passes; F01, F02, F04, F05 and F06 are killed (each by the md-as-report test); F03 is equivalent, as the log says. Round 92's five equivalents were not rerun. |
| 4 | fixed, one residue | Done: B4's label (md :315, code :951 and `b4_job`'s docstring); A6 and D headings; README row; A6E "at these K" (md) and "at K = 20, 40, 80" (D58); "both axes' intervals \|η\| < 2"; "the resolution of the stored doubles" twice; D‴.61's σ restored. Residue: the A6E paragraph's last clause, md :183 and code :1654–1655, still reads "t(v) < 0.8053, modulo a grid-checked domination condition" (D3cq-5). |
| 5 | fixed, one residue | Cited as [snippets, referee round 92] and "unchecked … against (v)'s cosine sums". The descriptions do not fit two of the four (D3cq-6). |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cq-1 | minor | **The Status is not re-partitioned.** Status (ii) (:183–186) gives as the conjectured curve's evidence: the pooled 6.1·10⁻⁵ / 1.4·10⁻⁴ over 32 rules, both set at (2, ½·v_hs), ρ = −0.736, proved; and A5's 2.3·10⁻⁵, set at (¼, Mercer frame), ρ = q, t(v) = 0.0154 < 0.0966, proved even without the branch and bound [own]. These now measure the extrapolation, not the conjecture. The values called "the conjecture's" at :225–226 and :296 (md :123, code :1361) are D‴.61-proved at all four τ/σ. | Status (ii): "at the 14 rules without a frame (where the curve rests on A6 alone) within 7.6·10⁻⁶ and 4.5·10⁻⁵, the two fits bracketing 1/t(v) at all 5 on the slow branch; at the 18 frame rules the deviations, up to 6.1·10⁻⁵ and 1.4·10⁻⁴, are the fits' bias". Move A5's 2.3·10⁻⁵ to "the extrapolation's bias at a proved limit". "The conjecture's values" → "D‴.61's values". |
| D3cq-2 | minor | **The split pools branches.** At the proved rules the fits bracket 8/8 slow-far, 2/3 slow-near, 0/4 fast-far and 0/3 fast-near. Bracketing is a calibrated interval only on the slow branch. There the 5 rules without a frame bracket 5/5; "9 of 14" adds 4 of 8 fast-branch brackets, which mean nothing, and "10 of 18" hides the 8/8 [own].<br>Both no-frame maxima are set by one slow-branch rule, (1, ½·v_hs). It lies on the proved rules' bias-against-t(v) trend ((2, 0.8·v_hs): +1.1·10⁻⁵ / −5.2·10⁻⁵).<br>On the fast branch the no-frame outlier (¼, 2·v_hs), +8.7·10⁻⁸ / +3.6·10⁻⁷, exceeds every proved fast rule's deviation (≤ 2.1·10⁻⁹, ≤ 1.7·10⁻⁷).<br>**Extension** [own]:<br>• The proved rule's deviations shrink (+6.07 → +1.10·10⁻⁵; −1.39·10⁻⁴ → −4.00·10⁻⁵, fits ending at K = 320): "their bias" is confirmed.<br>• (1, ½·v_hs) shrinks alike, to +8.2·10⁻⁷ / −8.3·10⁻⁶, still bracketing.<br>• (¼, 2·v_hs) reaches +4.9·10⁻¹² / −6.7·10⁻⁹ through K = 160–240: pre-asymptotic at K = 80, not a deviation. | `a6_lines`: split by frame × branch (and near/far), with the bracket count per cell. D58 :240–242: "on the slow branch, where the fits bracket all 8 proved limits, they bracket 1/t(v) at all 5 rules without a frame, deviating as the proved rules of like t(v) do". Optionally add A6 at K = 200, 240, 320 for these rules (minutes each). |
| D3cq-3 | minor | **Hypotheses not everywhere.**<br>• Status (ii) :181–182 "the rates λ/(2 − λ) and q proved by D‴.61 (iv) while below its r*" has no "modulo its branch and bound" and no "conjectural beyond".<br>• :240 "whose limits D‴.61 proves" and the md split label "whose limits D‴.61 (iv) proves" have no proviso.<br>• New :269–270 "the rate t, strictly faster than the Fock and Mercer-zero scalings'" is unqualified. That those rules fall at λ/(2 − λ) and q is D‴.61's (σ*/σ < 4.187, τ/σ < 6.514, modulo the branch and bound). At the Fock rule the unconditional bound gives (a/2)/t = (s + 1)²/(2(s² + 1)) ≤ 1 [own], so outside that range the claim is conjectural. | Add "(modulo its branch and bound; conjectural beyond)" at :181–182. Add "modulo its branch and bound" at :240 and in the md label. New: "strictly faster than the Fock and Mercer-zero scalings' (proved for σ*/σ < 4.187 and τ/σ < 6.514 with D‴.61, conjectural beyond)". |
| D3cq-4 | cosmetic | **Bounds rounded the wrong way** [own]:<br>• "within 7.5·10⁻⁶ and 4.4·10⁻⁵" (7.547·10⁻⁶, 4.449·10⁻⁵): D58 :242, md :171, README;<br>• "within 3.1·10⁻⁷ of 1/t" (3.134·10⁻⁷): :185, :242, :302, md :171;<br>• "above r* once σ*/σ > 3.044" (3.044165): false on (3.044, 3.04417).<br>`:.1e` and `:.3f` round to nearest. 4.187, 4.066 and 6.514 are safe by luck. | Print upper bounds rounded up (ceil at the shown digit): 7.6·10⁻⁶, 4.5·10⁻⁵, 3.2·10⁻⁷, 3.045. Lower thresholds keep floor. Align HANDOFF. |
| D3cq-5 | cosmetic | **Stale labels.**<br>• md :183 / code :1654–1655: "t(v) < 0.8053, modulo a grid-checked domination condition".<br>• "The conjecture's" values (D3cq-1).<br>• Status (i) "checked on a grid, to hold on \|η\| ≤ 1.98": D‴.61 checked 57 points on \|η\| = 1.9, 1.95, 1.98.<br>• The A6 and D headings use the literal 0.8052, not D61_RSTAR (G17 and G20 are caught only by the md test).<br>• Since 420a3ef, D‴.61 (cace1fd) reads "σ*/σ < 49.89" and "on the annulus 1.8927 < \|η\| ≤ 1.98". | Update the A6E clause from D61_RSTAR. Write "checked at 57 points on three circles". Format the headings from D61_RSTAR. Next revision: follow D‴.61's 49.89 and "the annulus". |
| D3cq-6 | minor | **The QFF citations are misdescribed** [snippet]. Only Mutný–Krause (NeurIPS 2018) and SLEIPNIR (Angelis–Wenk–Schölkopf–Bauer–Krause, arXiv:2003.02658) are seen to be Gauss–Hermite QFF. Their bound, √π m^{−m}(e/(4ℓ²))^m for \|t_i − t_j\| ≤ 1, is uniform on a bounded set and superexponential: not (v)'s geometric \|t\|^K in L²(ℝ).<br>TQFF (Li–Balakirsky–Mak, AISTATS 2024, arXiv:2310.14544) uses "a novel non-Gaussian quadrature rule". arXiv:2609.00553 (Cloninger–Le Gia–Mhaskar) bounds eigenvalues of quadrature-discretised kernels (Gaussian kernel, NTKs on the sphere); Gauss–Hermite is not confirmed. | :257–259: name the authors and describe each by what was seen. Keep "unchecked" and add "to be checked against the papers". Optionally add Dao–De Sa–Ré (NeurIPS 2017, Gaussian quadrature for kernel features) [snippet]. |
| D3cq-7 | cosmetic | **Mutants.**<br>• Y00 passes. The tests kill 17 of 24 (16 also change the md; G19 only the tests).<br>• Survive both checks, all equivalent on the data, none changing a published result: G01 (no \|ρ\| in [1 − 10⁻¹², 1)); G03 (every row has 'v'); G09 (20 bisections move the thresholds < 5·10⁻⁵); G10 (return hi); G11 (≤ r*); G13 (hi = 8 > 6.514); G18 (the signed spread is positive at its maximum).<br>• Without the md-as-report test, 8 md-changing mutants survive the unit tests: G06 (no-frame "7.5e-06" → 6.1e-05), G07 (signed second fit: 1.4e-04 → 8.8e-06), G08 ("of the 32"), G12 and G15 (Mercer-zero threshold), G17, G20, G22. The round's new numbers are pinned only by the md regression; `test_round92_*` checks labels, Fock and the ends, but not Mercer-zero. | Add a unit test with planted A6 rows (known ρ, d1, d2) asserting the split's maxima and counts per cell. Assert the closed forms 4.18735, 6.51445 and 3.04417 independently of `rate_conj`. Record G01, G03, G09–G11, G13 and G18 as equivalent. |

## Sources

- [Mutný–Krause, NeurIPS 2018 slides](https://NIPS.cc/media/nips-2018/Slides/12676.pdf); [NeurIPS 2018 spotlight page](https://neurips.cc/virtual/2018/spotlight/12676); [QFF code](https://github.com/Mojusko/QFF) [snippet]. Deterministic quadrature Fourier features; the error decreases exponentially with the number of features. The proceedings PDF did not resolve.
- [SLEIPNIR, arXiv:2003.02658](https://arxiv.org/pdf/2003.02658) [snippet]. Angelis, Wenk, Schölkopf, Bauer, Krause. Restates the QFF bound E_m = √π m^{−m}(e/(4ℓ²))^m for \|t_i − t_j\| ≤ 1.
- [Trigonometric QFF, arXiv:2310.14544](https://arxiv.org/pdf/2310.14544); [PMLR v238](https://proceedings.mlr.press/v238/li24o.html) [snippet]. Li, Balakirsky, Mak. A non-Gaussian quadrature rule.
- [Spectral Bounds for Kernel Quadrature, arXiv:2609.00553](https://arxiv.org/abs/2609.00553) [snippet]. Cloninger, Le Gia, Mhaskar. Eigenvalue and spectral-norm bounds for quadrature-discretised kernels; Gaussian kernel and NTKs on S⁴.
- [Dao–De Sa–Ré, Gaussian Quadrature for Kernel Features, NeurIPS 2017](https://proceedings.neurips.cc/paper/2017/hash/62f91ce9b820a491ee78c108636db089-Abstract.html); [Belhadji et al., arXiv:2408.13231](https://arxiv.org/abs/2408.13231v1) [snippet]. Related, uncited.
- [Derevianko–Plonka, arXiv:2307.13587](https://arxiv.org/pdf/2307.13587) [snippet]. Exponential decay in a weighted L².
- [arXiv:2606.03855](https://arxiv.org/pdf/2606.03855) [snippet: title only]. D58's "Plonka–Riebe–Cuyt" authorship is still unseen.
- [Karvonen–Särkkä, arXiv:1803.09532](https://arxiv.org/html/1803.09532v4) and [arXiv:2412.05368](https://arxiv.org/pdf/2412.05368) [snippet: titles, as cited].
- **[memory].** The bias of a fixed-exponent extrapolation when the true correction exponent lies between the assumed ones.

## Lesson

**Partition by mechanism as well as by provenance.** Once a theorem covers part of a test bed, the proved cases calibrate the instrument, but only for the mechanism they share with the unproved ones. Here the fits bracket every proved limit on the slow branch and none on the fast one. The calibrated evidence is therefore 5 of 5 slow-branch brackets, not "9 of 14". Where an extrapolation's bias is the issue, more K is cheap and decisive: at K ≤ 320 the bias shrinks at the proved rule and the unproved one alike, and the one outlier turns out to be pre-asymptotic.

**Round printed bounds outward.** Thresholds you claim to stay below round down; bounds you claim to stay within round up. Printing with `%.1e` does neither.
