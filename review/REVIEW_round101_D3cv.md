# Referee report D3cv, round 101: Proposition D‴.58 as revised after round 96

## Scope

- **Subject.** Commit 0cb96bd on claude/vigilant-clarke-hf939b (round 96 applied), diffed against 420a3ef. I exported it with `git archive 0cb96bd` into W/repo, with W = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref101`.
  - WM_OUT = W/wmout and GP_OUT = W/gpout, each holding a copy of its raw.
  - Scripts are in W/scripts, logs in W/logs, mutants in W/mut.
- **Read in full.**
  - D58.md, `results/chaos/width_mismatch.md`, `tests/test_width_mismatch.py` (58 tests; the two new `test_round96_*` line by line).
  - Every hunk of `exp_width_mismatch.py` changed since 420a3ef: `up_e`, `dn_f`, `up_f`, the 18/14 lines, `a6_lines`' cells, the threshold line, `a5_series_lines`' D‴.61 label, and the A6, D and A6E labels from D61_RSTAR.
  - REVIEW_round96_D3cq.md, the Round 96 entry of REVIEW_LOG.md and `tools/mutants/round96_D58_own.py`.
  - `git diff 420a3ef 0cb96bd` of these files, plus a word diff of D58.
  - For the citations: D61.md and gh_peano.md at 0cb96bd.
- **Run settings.** The repository's venv. Every process ran with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10` and a `timeout`, at most two of mine at a time. pytest always ran with `--basetemp` under W, including inside my copy of the mutant runner (W/scripts/mutants101.py, which differs from `tools/mutants.py` only in REPO, PY and `--basetemp`).
- **Nothing in the repository was modified.** I ran only `git archive/diff/log/show/status/rev-parse`. All my processes have exited; I stopped none.
  - Meanwhile another session moved HEAD to f3e0340. No D‴.58 file has changed since 0cb96bd.
  - D61.md and gh_peano.md have changed (D‴.61 round 99, the chord check now by a Taylor ratio). D58's citations still hold at HEAD: 0.8052, 9.26, 0.0966, 49.89, the annulus 1.8927 < |η| ≤ 1.98, 57 points.
- **Tags.** [md] the committed md; [raw] the committed raw; [own] my own code, or the repository's `a6_error` at new K with my own fits; [snippet] web-search results only (arxiv.org and ar5iv.labs.arxiv.org give ENOTFOUND); [memory] recalled facts.

## Summary

**Design.**
- Reproduced:
  - both mds byte for byte;
  - the tests;
  - eight jobs bit-identically from the raw.
- Recomputed independently:
  - A6's split, in exact rationals from the raw's ratios, with ρ, t(v) and the attaining argument in 50-digit Decimal from τ and the scalings' definitions;
  - every threshold, by closed form;
  - (i)'s exact identity, in 80-digit arithmetic.
- Extended A6 beyond K = 160 at five fast-branch rules, to K = 200, 240, 320:
  - four without a frame: (¼, 2·v_hs), (½, MM), (1, MM), (2, MM);
  - one proved, for calibration: (2, 2·v_hs).
- Mutants: 30 single-point mutants plus a control of this round's code, a second pass without the md-as-report test, and round 96's own spec rerun.
- Literature: six web searches.

**Verdicts.**
- **Reproducibility is exact.** `--from-raw` gives md5 6d7c1538… = committed (raw 959dc055…, unchanged); gh_peano.md is identical too. 58/58 tests pass (56.1 s). All eight jobs are identical in every field but `secs`.
- **Every number of the round traces to the md and the raw** [own]. My cells match the md's seven cells digit for digit. The thresholds are inward and right.
- **The re-partitioned Status is substantially right.** The Mercer frame's values are D‴.61's: t(v) ≤ 3/7 < r*, and at τ/σ = ¼, ½ it lies below 0.0966, so those two need no branch and bound.
- **The calibration by branch is sound on the slow branch.**
  - Its claim that the rules without a frame "deviate as the proved rules of like t(v) do" is true [own], but the md does not show it (D3cv-2).
  - On the fast branch D58 rightly discards the brackets. But it does not say that its quoted maximum, 8.8·10⁻⁸, is 41 times the proved fast rules' largest deviation, 2.1·10⁻⁹.
  - My extension shows that this is pre-asymptotic, and that the fast-branch deviations shrink with K at the unproved rules just as at the proved one (D3cv-1).
- **Outward rounding is done where round 96 asked, not everywhere.** The same quantity is printed as 8.7·10⁻⁸ and as 8.8·10⁻⁸ (D3cv-3).
- **No major issue.**

**Findings.**
1. **D3cv-1 (minor).** The fast-branch evidence for the rules without a frame is uncalibrated, and its maximum, an outlier, goes unflagged.
2. **D3cv-2 (minor).** "Deviate as the proved rules of like t(v) do" is not shown by the md: no frame column and no per-rule second fit. The counts are also loose: the 8 proved slow-far limits are 7 distinct rules, and "among the 8" should read 8 far plus 1 near.
3. **D3cv-3 (cosmetic).** Outward rounding is incomplete:
   - the "By the argument" line's ranges;
   - "within 8.7·10⁻⁸" next to "8.8·10⁻⁸";
   - "2.2·10⁻¹⁶" for 2⁻⁵²;
   - 0.4286 is outward only by luck.
4. **D3cv-4 (cosmetic).** Status (ii) residues:
   - the head clause "conjectural where a rule has no frame" omits frame rules with t(v) ≥ r*;
   - the cusp's 3.2·10⁻⁷ sits among the conjecture's evidence, though D‴.61 proves the cusp rule;
   - HANDOFF row 34 is stale.
5. **D3cv-5 (cosmetic).** Code and test semantics:
   - a5's label prints "< r*" without comparing;
   - the planted-row test checks labels only, though its docstring says "bounds";
   - 14 of 22 md-changing mutants are caught only by the md regression;
   - 5 mutants survive, all equivalent on the data, none changing a published result.
6. **D3cv-6 (cosmetic, optional).** For τ/σ < 1.5039 the Mercer-zero rule is strictly slower than t unconditionally, by D‴.61 (iv)'s lower bound. The Fock rule never is.

## What I checked and how

**Reproduction.**
- `--from-raw` into W/wmout and W/gpout, then `cmp`: identical.
- Full test file: 58 passed in 56.14 s.
- `W/scripts/rerun.py` compares each rerun job with the raw by `==` in every field but `secs` [raw, own]. All eight are identical:

| job | time |
|---|---|
| A6 (2, 1.25·v_hs) | 6.2 s |
| A6 (1, moment matching) | 17.5 s |
| A6 (¼, 2·v_hs), the fast-branch outlier | 142.7 s |
| C (s = 10⁻³, draw 1) | 0.4 s |
| B3 job 9 ((½, grid +0.05), K = 6) | 6.9 s |
| B2 job 5 | 0.2 s |
| D (K = 80, t = 0.172) | — |
| A6E (K = 20, real axis) | — |

**A6 split (`W/scripts/split.py`)** [own].
- Method:
  - both fits through K = 80, 120, 160, solved by Cramer's rule in Fractions from the stored doubles;
  - ρ = 1 + 2/τ² − 2/v, t(v) and the attaining argument in Decimal from τ and the scalings' definitions, not from the stored v.
- Checks:
  - stored rates agree to 5.6·10⁻¹⁷; stored v exactly; every argument agrees;
  - 22 frames; ρ = 1 at 4 rules; the largest frame |ρ| is 0.94113 and the largest frame t(v) is 0.634512.
- Results, each matching the md's printed cell:

| cell | n | max \|fit 1\| | max \|fit 2\| | bracket | md prints |
|---|---|---|---|---|---|
| frame, τ²/v, far | 8 | 6.0734·10⁻⁵ | 1.3900·10⁻⁴ | 8 | 6.1e-05, 1.4e-04, 8 ✓ |
| frame, τ²/v, near | 3 | 2.9114·10⁻⁵ | 5.3392·10⁻⁶ | 2 | 3.0e-05, 5.4e-06, 2 ✓ |
| frame, vσ*²/τ², far | 4 | 2.0981·10⁻⁹ | 1.7348·10⁻⁷ | 0 | 2.1e-09, 1.8e-07, 0 ✓ |
| frame, vσ*²/τ², near | 3 | 3.9877·10⁻⁵ | 3.6666·10⁻⁵ | 0 | 4.0e-05, 3.7e-05, 0 ✓ |
| no frame, τ²/v, far | 5 | 7.5472·10⁻⁶ | 4.4490·10⁻⁵ | 5 | 7.6e-06, 4.5e-05, 5 ✓ |
| no frame, vσ*²/τ², far | 8 | 8.7052·10⁻⁸ (¼, 2·v_hs) | 3.5976·10⁻⁷ (same) | 4 | 8.8e-08, 3.6e-07, 4 ✓ |
| no frame, vσ*²/τ², near | 1 (MM, ¼) | 5.2612·10⁻¹¹ | 3.7751·10⁻⁸ | 0 | 5.3e-11, 3.8e-08, 0 ✓ |

- Cusp: 3.1340·10⁻⁷ (printed 3.2e-07). Largest spread: 1.9974·10⁻⁴ (printed 2.0e-04).
- The "By the argument" line's true ends:

| branch (far) | fit 1 | fit 2 | excess at K = 160 | exponent |
|---|---|---|---|---|
| τ²/v | 5.652·10⁻⁹ … 6.0734·10⁻⁵ | −1.3900·10⁻⁴ … −3.5905·10⁻⁷ | 8.0717·10⁻⁶ … 4.0319·10⁻⁴ | 1.5246 … 1.9792 |
| vσ*²/τ² | −2.098·10⁻⁹ … 8.705·10⁻⁸ | — | — | 1.9961 … 2.0399 |

- The Mercer frame's constant: max relative deviation exactly 2.220446·10⁻¹⁶ at all four τ/σ.

**Slow branch, far: unproved against proved rules of like t(v)** [own].

| t(v) | rule | frame | fit 1 | fit 2 | excess at K = 160 | exponent |
|---|---|---|---|---|---|---|
| 0.101 | (½, MZ) | yes | +5.65·10⁻⁹ | −3.59·10⁻⁷ | +8.07·10⁻⁶ | 1.979 |
| 0.126 | (¼, 0.8·v_hs) | no | +6.97·10⁻⁹ | −4.56·10⁻⁷ | +8.39·10⁻⁶ | 1.975 |
| 0.166 | (½, 0.8·v_hs) | no | +2.39·10⁻⁸ | −1.04·10⁻⁶ | +1.37·10⁻⁵ | 1.963 |
| 0.200 | (1, Fock) | yes | +2.13·10⁻⁷ | −3.55·10⁻⁶ | — | — |
| 0.277 | (1, 0.8·v_hs) | yes | +3.67·10⁻⁷ | −6.37·10⁻⁶ | — | — |
| 0.347 | (¼, ½·v_hs) | no | +9.29·10⁻⁷ | −1.07·10⁻⁵ | — | — |
| 0.382 | (½, ½·v_hs) | no | +1.59·10⁻⁶ | −1.68·10⁻⁵ | — | — |
| 0.473 | (2, 0.8·v_hs) | yes | +1.10·10⁻⁵ | −5.21·10⁻⁵ | +1.56·10⁻⁴ | 1.704 |
| 0.478 | (1, ½·v_hs) | no | +7.55·10⁻⁶ | −4.45·10⁻⁵ | +1.50·10⁻⁴ | 1.765 |
| 0.500 | (2, Fock) ≡ (2, MZ) | yes | +1.42·10⁻⁵ | −6.15·10⁻⁵ | — | — |
| 0.635 | (2, ½·v_hs) | yes | +6.07·10⁻⁵ | −1.39·10⁻⁴ | — | — |

- The signs are the same at all 13 rules.
- Four of the 5 unproved rules lie between their proved neighbours in t(v) on both fits. The fifth, (1, ½·v_hs), lies at 0.68 and 0.85 times its proved neighbour (2, 0.8·v_hs).
- The excess at K = 160 and its exponent interleave the same way.
- The unproved rules have ρ from −2.2 to −33. Two of them sit at τ/σ = ¼, where no proved slow-far rule exists.

**A6 extended (`W/scripts/a6extend.py`)** [own]. The repository's `a6_error` and `a6_scalings` in its own precision rule (2K·log₁₀(1/t) + 60 digits); my fits through consecutive triples.

| rule | ρ | 80–160 (fit 1 / fit 2) | 120–200 | 160–240 | 200–240–320 |
|---|---|---|---|---|---|
| (¼, 2·v_hs), no frame | +16.5 | +8.705·10⁻⁸ / +3.598·10⁻⁷ | +3.517·10⁻¹⁰ / −1.195·10⁻⁸ | +4.917·10⁻¹² / −6.664·10⁻⁹ | +1.547·10⁻¹² / −3.338·10⁻⁹ |
| (½, MM), no frame | 1 | +1.929·10⁻¹⁰ / +1.912·10⁻⁸ | +5.75·10⁻¹¹ / +7.44·10⁻⁹ | +2.29·10⁻¹¹ / +3.66·10⁻⁹ | +9.03·10⁻¹² / +1.81·10⁻⁹ |
| (1, MM), no frame | 1 | −7.52·10⁻¹¹ / −5.89·10⁻⁸ | −2.25·10⁻¹¹ / −2.35·10⁻⁸ | −8.97·10⁻¹² / −1.17·10⁻⁸ | −3.53·10⁻¹² / −5.85·10⁻⁹ |
| (2, MM), no frame | 1 | −1.68·10⁻¹¹ / −3.52·10⁻⁸ | −5.38·10⁻¹² / −1.41·10⁻⁸ | −2.23·10⁻¹² / −7.03·10⁻⁹ | −9.00·10⁻¹³ / −3.51·10⁻⁹ |
| (2, 2·v_hs), proved | +0.941 | −6.94·10⁻¹¹ / −4.49·10⁻⁸ | −2.10·10⁻¹¹ / −1.79·10⁻⁸ | −8.44·10⁻¹² / −8.93·10⁻⁹ | −3.35·10⁻¹² / −4.46·10⁻⁹ |

- My K ≤ 240 values at (¼, 2·v_hs) reproduce round 96's to every printed digit.
- At (¼, 2·v_hs) the excess times K² is −0.0928, −0.0927, −0.0926 at K = 160, 240, 320: a clean K⁻², with the limit 1/t(v) to about 10⁻¹².
- Times: τ/σ = ¼: 291, 704, 1850 s at K = 200, 240, 320; (½, MM): 107, 254, 951 s.

**Thresholds** [own; 40-digit Decimal and bisection].
- r* from c² = 4 ln u/(u − 1), u = 1/r:

| c | r* | σ*/σ at the half-squeeze rule |
|---|---|---|
| 1.8927 | 0.8052030 | 9.2671 |
| 1.98 | 0.9607245 | 49.922 (49.8906 at the rounded r = 0.9607, D‴.61's 49.89) |
| 1 | 0.0966496 | — |

- c_need(0.6345) = 1.77733.
- At r* = 0.8052:
  - Fock: σ*/σ < 4.187348 (τ/σ < 4.066188);
  - Mercer-zero: τ/σ < 6.514453;
  - frame end: σ*/σ > 3.044165.
- All are printed inward: 4.187, 4.066, 6.514, 3.045.

**(i)'s identity** [own; `W/scripts/series_check.py`].
- Compared, at τ/σ = ½, 1, 2 and K = 1, 2, 3, 5, 8:
  - the closed-form L² error (Gaussian inner products in 80-digit Decimal; the nodes from `gh_nodes`);
  - and (1 − t²)^{1/2}Σ_{n≥2K}(tⁿ/n!)D²_{n,K}.
- They agree to between 10⁻⁷⁸ and 2.2·10⁻⁶⁰ (the cancellation's limit). For example, at τ/σ = 1, K = 2 the error is 4.7745·10⁻² = md [md].

**Unconditional lower bound** [own]. D‴.61 (iv) gives a/(2(1 + ρ)) for ρ < 0 and a/2 otherwise, with a = v(1 − ρ²)/2. Against t:
- Mercer-zero: 1.98, 1.65 and 1.23 times t at τ/σ = 0.1, ½, 1; above t exactly for τ/σ < 1.50395.
- Fock: never above t (0.97 at τ/σ = 1).

**D‴.61 cross-check.** At 0cb96bd and at HEAD: (iii)(b) and (iv), the Status and New, and gh_peano.md lines 137–143.
- r* = 0.8052, 9.26, r < 0.0966.
- r < 0.9607 and 49.89 "if the valley condition holds on the annulus 1.8927 < |η| ≤ 1.98".
- 57 points at |η| = 1.9, 1.95, 1.98.
- "floating point with slack 10⁻¹², not interval arithmetic".
- D58 cites all of these correctly.

**Mutants.**
- **Pass 1**: W/mut/spec101.py, N00–N30, against the full suite.
  - The control N00 passes.
  - 25 are killed. 22 change the md, and the full suite kills all 22.
  - 5 survive, all with the md unchanged:

| mutant | what it changes | why it is equivalent on the data |
|---|---|---|
| N03 | `up_e`'s guard 10⁻¹² → 10⁻⁹ | no printed value lies within 10⁻⁹ of a printable one |
| N06 | a5's t(v) from the fit's L | 2.333333 either way |
| N08 | nearness `>=` → `>` | — |
| N10 | bracket `<= 0` → `< 0` | — |
| N20 | `<=` in the bisection | — |

- **Pass 2**: the 22 md-changers with `-k 'not round88_precision_and_the_md'`.
  - The unit tests kill 8: N09, N11, N12, N18, N19, N21, N22, N30.
  - The other 14 are caught only by the md regression:
    - N04: a5's deviation min instead of max;
    - N05: a5's t(v) min, which would print "at most 0.0154 < r*";
    - N07: 0.4286 printed to 3 decimals;
    - N13: the near cells relabelled with far data;
    - N14–N16: the 18/14 lines' maxima and bracket count;
    - N17: the spread's min;
    - N23–N28: the r* labels and the provisos.
  - N01, N02 and N29 are caught by the tests without changing the md.
- **Round 96's spec** (`round96_D58_own.py`), rerun as is: X00 passes; all 16 are killed, L02, L03 and L11 by the tests alone. This matches the log.

**Literature.** Six searches (Sources).
- Mutný–Krause (NeurIPS 2018, Gauss–Hermite QFF, from its reviews page) [snippet].
- SLEIPNIR's authors, and that it restates the QFF bounds, with bounds of the form (a/m)^m and a threshold involving e/(2ℓ²) [snippet].
- TQFF: Li, Balakirsky, Mak, AISTATS 2024, PMLR 238, "a novel non-Gaussian quadrature rule" [snippet].
- Dao–De Sa–Ré, NeurIPS 2017: deterministic Gaussian-quadrature feature maps, O(e^{e^γ} + ε^{−1/γ}) samples [snippet].
- Not reconfirmed this round:
  - the exact constant √π m^{−m}(e/(4ℓ²))^m;
  - the authors of arXiv:2609.00553 (its snippet shows a spectral-norm theorem and arccosine-kernel tables, no author list).
  D58 tags both "[snippets, referee rounds 92, 96; to be checked against the papers]".

## Verdicts per part

| claim (D58 lines) | verdict | evidence |
|---|---|---|
| (i), its proof and Status (i) "proved by D‴.61 (iv) for t < 0.8052 (σ*/σ < 9.26; … without it for t < 0.0966)"; "σ*/σ < 49.89 were D‴.61's valley condition, checked at 57 points on three circles, to hold on the annulus 1.8927 < \|η\| ≤ 1.98" (:14–15, :170–172) | Correct, rounded inward. The identity re-verified to 10⁻⁶⁰ | [own]; D61 at 0cb96bd and HEAD |
| (ii): "proved by D‴.61 (iv), modulo its branch and bound, for every v ∈ (τ²σ²/σ*², τ²) … t(v) < 0.8052, conjectural for the others" (:24–25) | Correct | [own] |
| Status (ii): Fock and Mercer-zero "proved … while below its r*, modulo its branch and bound, conjectural beyond (σ*/σ < 4.187 and τ/σ < 6.514)" (:181–183) | Correct and inward. Mercer-zero's order is unconditional for τ/σ < 1.5039 | [own]; D3cv-6 |
| Status (ii): the 14 rules without a frame (7.6·10⁻⁶, 4.5·10⁻⁵, 5 of 5 slow-far brackets, all 8 proved); the frame rules' 6.1·10⁻⁵ / 1.4·10⁻⁴ and A5's 2.3·10⁻⁵ "the fits' bias" (:183–188) | Numbers right; re-partitioned as asked. Head clause and the cusp's place are loose | [md, own]; D3cv-4 |
| A6 paragraph: "on the slow branch far from the cusp the two fits bracket all 8 proved limits and 1/t(v) at all 5 rules without a frame" (:246–247) | Right (8 = 7 distinct rules) | [own]; D3cv-2 |
| "… which deviate as the proved rules of like t(v) do" (:247–248) | True [own]: signs, magnitudes, excess and exponent interleave. Not shown by the md (no frame or ρ column, no per-rule second fit) | D3cv-2 |
| "on the fast branch … none of the 7 proved limits, so its 4 brackets among the 8 rules without a frame carry nothing, and the deviations there are at most 8.8·10⁻⁸ and 3.6·10⁻⁷" (:248–249) | Logic right; the bounds right. Two things missing: the 8.8·10⁻⁸ is 41 times the proved fast-far maximum, set at one pre-asymptotic rule; and the fast branch is calibrated by the shrinkage with K, which K ≤ 160 cannot show. "8" should read 8 far plus 1 near | [own] extension; D3cv-1 |
| A6 sentences unchanged except for the proviso (:239–246): 13/12 rules, the ranges, 13 of 13, "within 7.6·10⁻⁶ and 4.5·10⁻⁵ (… bounds rounded up, round 96)" | Right. "−2.1·10⁻⁹ … +8.7·10⁻⁸" is rounded inward at its upper end | [own]; D3cv-3 |
| A5: "the formula's values (b² + b + 2)/(b(b − 1)) — D‴.61 (iv)'s limits there, the frame's t(v) at most 0.4286 < r*, modulo its branch and bound" (:228–230, :307; md line 123) | Correct. t(v) = 0.015371, 0.058258, 0.188345, 3/7; ρ = q ∈ (0, 1); v ∈ frame interval. τ/σ = ¼, ½ need no B&B. 0.4286 outward by luck | [own]; D3cv-3, D3cv-5 |
| md threshold line (178) and New "strictly faster … (proved with D‴.61 for σ*/σ < 4.187 and τ/σ < 6.514, modulo its branch and bound; conjectural beyond)" (:280–281) | Correct, inward | [own] |
| Headings and A6E clause from D61_RSTAR, "modulo D‴.61's floating-point branch and bound" | Correct | [md] |
| QFF sources (:265–270) | As seen. The exact QFF constant and the authors of 2609.00553 not reconfirmed this round; tagged to be checked | [snippet] |
| (iii)–(v) and their proofs; B, C, D numerics | Unchanged since 420a3ef. Reruns identical; spot checks against the md pass (93.0, 0.963, 31769, 3643, 1.0484, 1.0235, 2.677, 8.35·10³) | [raw, md] |

What is proved, and how:
- **Proved:**
  - (i)'s identity, its bounds and the rate t;
  - (v);
  - the order of λ, λ/(2 − λ) and q against t.
- **Proved via D‴.61, modulo its floating-point branch and bound:**
  - the sharp constant for t < 0.8052 (unconditionally for t < 0.0966);
  - every frame rule's rate t(v) < 0.8052: Fock below σ*/σ = 4.187, Mercer-zero below τ/σ = 6.514, the Mercer frame's limits, the cusp, all 22 frame rules of A6.
- **Conditional on D‴.61's valley condition:** up to r = 0.9607 (σ*/σ < 49.89).
- **Numerical:** moment matching's λ (A3), the 14 rules without a frame (A6), (iii), and (iv)'s numbers.
- **Conjectural:** t(v) at rules without a frame, or with t(v) ≥ r*; the sharp constant beyond.

## Status of the round-96 issues

| D3cq | status | note |
|---|---|---|
| 1 | fixed, small residue | The Status cites the 14 rules without a frame; the frame rules' maxima and A5's 2.3·10⁻⁵ are called the fits' bias; "the formula's (D‴.61's)" in D58, the md and the code. Residue: the head clause, and the cusp's 3.2·10⁻⁷ among the evidence (D3cv-4) |
| 2 | fixed, residues | The cells are in the md (all seven verified) and D58 has the calibration sentence. Residues: "deviate as" not visible in the md; the counts; the fast-branch outlier unflagged; the extension in the log only (D3cv-1, D3cv-2) |
| 3 | fixed | The proviso in Status (ii), the A6 sentence, the md's frame label, the threshold line and New |
| 4 | fixed where named | 7.6, 4.5, 3.2 and 3.045 rounded outward; 4.187, 4.066 and 6.514 floored. Residues: the "By the argument" line, 2.2·10⁻¹⁶, 0.4286 (D3cv-3) |
| 5 | fixed | The A6E clause and the headings from D61_RSTAR; 49.89, "57 points on three circles", the annulus. HANDOFF row 34 is still stale (D3cv-4) |
| 6 | fixed | The authors named; each source described by what was seen; "to be checked against the papers" |
| 7 | fixed | Two new tests; the threshold closed forms (4.18735, 6.51445, 3.04417) asserted; round 96's 16 mutants all killed (rerun). The semantics of the 18/14 lines, the near cells and a5's label are still pinned only by the md regression (D3cv-5) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cv-1 | minor | **The fast branch is uncalibrated, and its maximum is an unflagged outlier.**<br>• D58 :248–249 offers "at most 8.8·10⁻⁸ and 3.6·10⁻⁷" as the fast-branch evidence for the rules without a frame.<br>• Both maxima are set by one rule, (¼, 2·v_hs), ρ = 16.5. They are 41 times and 2.1 times the proved fast-far maxima (2.1·10⁻⁹, 1.8·10⁻⁷) [own].<br>• The only resolution, round 96's extension, is in the log only.<br>**Extension** [own]:<br>• (¼, 2·v_hs) through K = 200–240–320: +1.5·10⁻¹² / −3.3·10⁻⁹; excess −0.0926/K² at K = 160, 240, 320.<br>• At (½, MM), (1, MM) and (2, MM) both fits' deviations shrink about 3–20 times from the 80–160 window to 200–320, as at the proved (2, 2·v_hs) (−6.9·10⁻¹¹ → −3.3·10⁻¹²; −4.5·10⁻⁸ → −4.5·10⁻⁹).<br>• On the fast branch the calibration is this shrinkage, not a bracket. | D58 :248–249 and the md's cells: "… at most 8.8·10⁻⁸ and 3.6·10⁻⁷, set at (¼, 2·v_hs), pre-asymptotic at K ≤ 160 (through K = 200–320 within 1.6·10⁻¹² and 3.4·10⁻⁹; referee rounds 96, 101, not reproduced by the repository's scripts); at the proved fast rules both fits lie below the limit and shrink with K, as at the rules without a frame". Better: add A6 at K = 200, 240, 320 for the fast-branch rules (minutes each, 31 min at τ/σ = ¼) and print the deviations by window |
| D3cv-2 | minor | **"Deviate as the proved rules of like t(v) do" (D58 :247–248) is not shown by the md.**<br>• The A6 table marks no rule as having a frame and prints only the first fit's deviation per rule.<br>• The claim is true [own]: same signs at all 13 slow-far rules; 4 of 5 between proved neighbours on both fits; (1, ½·v_hs) at 0.68 and 0.85 of (2, 0.8·v_hs); the excess and exponents interleave too.<br>• The slow-far calibration has no proved rule at τ/σ = ¼, where 2 of the 5 unproved rules lie.<br>• "All 8 proved limits" are 7 distinct rules: at τ/σ = 2 Fock ≡ Mercer-zero (v = 4/3).<br>• "Its 4 brackets among the 8 rules without a frame": the fast branch has 9 such rules (8 far, MM ¼ near), while "none of the 7 proved" counts near ones too. | `a6_lines`: add ρ (or a frame mark) and the second fit's deviation to the A6 table. D58 :246–249: "8 (7 distinct)", "among the 8 far ones (9 with MM at τ/σ = ¼)", and "deviate as … (A6's table: same signs, interleaved in t(v))" |
| D3cv-3 | cosmetic | **Outward rounding incomplete.** The "By the argument" line (md 173) still uses `:+.1e`:<br>• D58 :240 "−2.1·10⁻⁹ … +8.7·10⁻⁸" and :311 "within 8.7·10⁻⁸", for 8.705·10⁻⁸, printed 8.8·10⁻⁸ at :249;<br>• "+5.7·10⁻⁹" (5.652), "−3.6·10⁻⁷" (−3.591), and "+8.1·10⁻⁶ … +4.0·10⁻⁴" (8.072·10⁻⁶ … 4.032·10⁻⁴): range ends rounded inward.<br>Elsewhere:<br>• "to 2.2·10⁻¹⁶" (2⁻⁵² = 2.2204·10⁻¹⁶) at :221, :313 and md 180;<br>• "t(v) at most 0.4286" is `:.4f`, outward only because 3/7 = 0.428571. | Ranges outward (lower end down, upper end up) with an `up_e`/`dn_e` pair; `up_f(…, 4)` for t(v); `up_e` for the Mercer constant; D58 follows |
| D3cv-4 | cosmetic | **Status (ii) residues.**<br>• The head clause "conjectural where a rule has no frame" omits frame rules with t(v) ≥ r*: the Fock rule beyond σ*/σ = 4.187, which the same sentence calls conjectural, and the frame ends once σ*/σ > 3.045.<br>• "At the cusp within 3.2·10⁻⁷ of 1/t" stands among the conjecture's evidence, but the cusp rule is the frame ρ = t, which D‴.61 proves at all four τ/σ (t ≤ 0.382).<br>• HANDOFF row 34 still says "one conjectured formula t(v) for every scaled rule's rate". | "conjectural where a rule has no frame or t(v) ≥ r*"; "at the cusp, proved, within 3.2·10⁻⁷"; update HANDOFF |
| D3cv-5 | cosmetic | **Code and test semantics.**<br>• `a5_series_lines` prints "< r* = 0.8052" as fixed text, without comparing; the threshold line's "a frame at every τ/σ" is likewise text.<br>• `test_round96_…_thresholds`' planted rows assert labels and counts only, though the docstring says "counts and bounds".<br>• 14 of 22 md-changing mutants are caught only by the md regression (N04, N05, N07, N13–N17, N23–N28).<br>• Survivors N03, N06, N08, N10 and N20 are equivalent on the data; none changes a published result. | Compute "<"/"≥" from D61_RSTAR. Assert the planted cells' maxima and brackets, the 18/14 lines' numbers, the near cells and a5's t(v) maximum in unit tests. Record the five survivors as equivalent |
| D3cv-6 | cosmetic | **Optional strengthening.** D‴.61 (iv)'s unconditional bound a/(2(1 + ρ)) shows the Mercer-zero rule strictly slower than t for τ/σ < 1.50395 (1.23 times t at τ/σ = 1; 1.65 at ½) [own]. Fock's a/2 never does. | New / Status (ii): "(for the Mercer-zero rule unconditionally while τ/σ < 1.503)" |

## Sources

- [NeurIPS 2018 reviews of Mutný–Krause](https://papers.neurips.cc/paper_files/paper/2018/file/6e923226e43cd6fac7cfe1e13ad000ac-Reviews.html) [snippet]. Gauss–Hermite quadrature for Bochner's integral (QFF).
- [SLEIPNIR, arXiv:2003.02658](https://arxiv.org/abs/2003.02658) [snippet]. Angelis, Wenk, Schölkopf, Bauer, Krause; restates the QFF bound; bounds of the form (a/m)^m. The PDF did not resolve.
- [Trigonometric QFF, PMLR v238](https://proceedings.mlr.press/v238/li24o.html) and [arXiv:2310.14544](https://arxiv.org/abs/2310.14544v1) [snippet]. Li, Balakirsky, Mak: a non-Gaussian quadrature rule.
- [Spectral Bounds for Kernel Quadrature, arXiv:2609.00553](https://arxiv.org/pdf/2609.00553) [snippet]. A spectral-norm theorem and arccosine-kernel tables; no authors shown.
- [Dao–De Sa–Ré, NeurIPS 2017](https://proceedings.neurips.cc/paper/2017/hash/62f91ce9b820a491ee78c108636db089-Abstract.html) [snippet]. Deterministic Gaussian-quadrature feature maps.
- [Karvonen–Särkkä, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532) [snippet: title and node bound].
- **[memory].** A fixed-exponent extrapolation is biased when the true correction is pre-asymptotic, and the bias shows as drift with the fit window.

## Lesson

**A bracket that fails at every proved case is a symptom, not neutral evidence.** On the fast branch both fits sit on the same side of all 7 proved limits. That branch's calibration is therefore how the deviation shrinks as the fit window moves out in K. More K settles that cheaply: here the one outlier is pre-asymptotic, and the unproved rules converge as the proved one does. Print that, not a bound.

**Outward rounding is a property of every printed bound, not only of those a referee named.** One sweep for `:.1e` and `:.4f` in the lines that print bounds and ranges would have caught 8.7 against 8.8, 2.2·10⁻¹⁶ and 0.4286.
