# Referee report D3bv, round 75: Proposition D‴.52 as revised after round 73

## Scope

- **Commit reviewed:** 3dc9443 on branch claude/festive-bardeen-eopej2 ("D‴.52 round 73: the tie's zero counted; the part's errors; size by distance").
  - I exported it with `git archive 3dc9443 | tar -x` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/repo`.
  - All other work is under `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work` (scripts/, logs/, out/, mut/).
  - HEAD moved to eeef69c ("D‴.60 round 69") while I worked. `git diff 3dc9443 eeef69c` is empty for every file reviewed here [own], so I reviewed 3dc9443 only.
- **Files read:**
  - D52.md in full;
  - in `exp_close_features.py`: `e_prime_lines`, `own_knot`, `f_alloc_lines`, `u_job`, `knot_of`, `knot_of_kernel`, `e_job`, `s_job`, `f_job`, `rejections`, `equal_size_line`, `report` and `main`;
  - close_features.md (sections E, E′ and F);
  - `test_close_features.py` (47 tests);
  - REVIEW_round73_D3bt.md, and the Round 71 and Round 73 entries of REVIEW_LOG.md;
  - the diff 1f51d48..3dc9443, and `rich_transit.Bumps`.
- **Settings:**
  - the repository's venv, run with `nice -n 10`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`;
  - pytest as `python -m pytest -q -p no:cacheprovider`;
  - `CF_OUT` always pointed to a scratch directory holding a copy of the raw;
  - at most two heavy processes of mine at a time.
- **I modified nothing in the repository**, and nothing of mine is still running.
  - `git status` shows one modified file, `review/drafts/D58.md`. It is not mine: another session on the machine is editing it.
- **Tags:** [own] is my run or computation, [md] is close_features.md, [raw] is close_features_raw.npz.

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` regenerates the md byte for byte (`cmp`; md5 a4455d5d…) [own].
  - The 47 tests pass in 20.9 s [own].
  - The raw is byte-identical to round 71's (md5 670e6d17…) [own].
  - Rerun with the round-73 `u_job`, all 15 U jobs reproduce the raw's U in every field, mu0 included [own].
- **The numbers match the md.**
  - All 46 decimal numbers on the lines round 73 added to D52.md appear verbatim in the md [own]. ("0.12" is the magnitude of 0.1225.)
  - "3293 of 21000", "0 of 2400" and "6 of 6" are the md's [md].
  - Over the whole file: 521 numbers verbatim, 43 roundings and 3 derived (90.7%, 4.7%, 5.6·10⁻⁵) [own]. All 46 are on lines round 73 did not touch.
- **The inequality and the decomposition are right, and the code computes what the text says** [own].
  - Write F for the full knot and O for the pair's own knot.
  - F ≤ O draw by draw. Quantiles taken from the same samples therefore satisfy q_F ≤ q_O.
  - At a common threshold q: part = P(F < q) − P(O < q) = P(F < q ≤ O) − P(O < q ≤ F) = P(F < q ≤ O) ≥ 0.
  - Off the tie: part = P(F < q_F ≤ O) − P(q_F ≤ O < q_O). That is the minimum's gain at the lower threshold, less the cost of that lower threshold.
  - `betw` counts F < q₂ ≤ O on the H1 draws of the cells where q₂ = q_full. There q_full = q_law to 0.0. The count is 0 of 2400 [own = md].
  - The count is not vacuous: off the tie the same event occurs in 12–82 of the 400 draws per cell [own].
- **The excess statistic is round-off, as stated** [own].
  - In 7071 of the 21000 draws the own kernel sets the full knot, up to round-off. Among these, the path's knot is above U's knot in 3293, equal in 567, and below it by less than 10⁻¹⁰ relative in 3211.
  - That is a balanced sign, which is what round-off gives. The maximum is 6.5·10⁻¹².
  - Every knot is finite (21000 of 21000).
  - Across 15·10⁶ law samples, `knot_of` ≤ `knot_of_kernel` holds exactly.
- **The paired SE is the McNemar form, with ddof = 1.** Using ddof = 0 changes nothing at three decimals [own].
- **The tie's zero holds in the population, and the printed SEs are right** [own].
  - Population check: 10⁶ shared law samples and 10⁵ fresh H1 draws per cell.
  - The vectorised fits were validated on E's 21 000 draws: centres agree to 1.5·10⁻¹⁵, and knots to 5.3·10⁻⁷ relative.
  - At c̄ ≥ 0.045: q_full = q_own in 6 of 6 cells, and 0 of 600 000 draws fall between the two knots.
  - Off the tie (s = 0.003, 0.01, 0.03), SE ≈ 0.001:
    - c̄ = 0: −0.0833, −0.0711, −0.0816;
    - c̄ = 0.02: −0.0792, −0.0678, −0.0872;
    - c̄ = 0.032: −0.0344, −0.0340, −0.0383.
  - Standardised by the printed SEs, the md's nine parts give χ² = 8.0 on 9 df. The −0.1225 cell lies 2.0 SE from its population value.
- **Who sets the null's tail** (l. 70) [own].
  - At c̄ ≤ 0.032, other kernels set 34–82% of the law's 5% tail: the fit's other neighbour 9–19%, far kernels 25–63%.
  - At c̄ ≥ 0.045 the tie pair sets 100% of it.
- **The size by distance reproduces.**
  - F's line is 0.0512, 0.0456, 0.0398, 0.0379, 0.0329 [own = md].
  - Each mirror pair agrees within each null (Fisher p = 0.34–0.59) [own].
- **The round-73 code is well tested.** Of 44 single-point mutants of it, the tests kill 30 and the md byte diff catches 3 more. The 11 survivors are equivalent on this data or immaterial [own].
- **D3bt-1 and most of D3bt-5 are fixed. D3bt-4 is fixed as proposed**, though its new wording over-reads (D3bv-1).

**What is wrong (majors first)**

There is no major issue. There are four minor issues and three cosmetic ones.

1. **D3bv-1 (minor). The new size-by-distance wording pools two nulls that behave differently** (l. 54–58).
   - The fall with distance belongs to the weak null, |ŵ| ≈ 0.4: trend z = −3.83. The strong null, |ŵ| ≈ 1.4, shows no significant trend: z = −1.45, p = 0.15 [own].
   - At s = 0.03 the (0.8, −0.4) nulls sit 0.0097 from their c̄ label [raw], so for them "mirror" pooling does not pool mirrors.
2. **D3bv-2 (minor). (iii) still quotes the noisy extremes of the off-tie loss.** It says "0.025–0.12" [md]; the population range is 0.034–0.087 [own].
3. **D3bv-3 (minor). The literature fit of the new sources.**
   - Azaïs–Gassiat–Mercadier 2006 is known for local power close to the level, not for the N^{−1/4} order [snippet].
   - Dacunha-Castelle–Gassiat and Azaïs–Gassiat–Mercadier are results on population mixtures.
   - The tag [snippets] on Amar–Weiss's δ⁻⁴ is unverified.
4. **D3bv-4 (minor). Older mutants still escape the tests** [own].
   - 24 of the 27 mutants that round 73 caught only by the md diff still escape the tests.
   - A06, ES03, FR01 and FR03 still survive both checks.
5. **D3bv-5 to D3bv-7 (cosmetic):**
   - E′ prints zero standard errors at the tie;
   - the excess statistic silently drops infinite knots;
   - wording and log loose ends.

## What I checked and how

**Reproducibility**
- `--from-raw` with `CF_OUT=/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/out/fromraw`, followed by `cmp` against the commit's md and md5 against `git show 3dc9443:` and `1f51d48:`.
- The test file: output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/out/pytest_base.txt`.
- All 15 U jobs rerun with part_U's arguments, compared field by field. This is the identity mutant in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/mut/results_none.jsonl`.

**Numbers**
- Every decimal number in D52.md checked against the md, with the round-73 lines listed separately: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/scripts/numbers_check.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/scripts/numbers_round.py`.
- Output: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/out/numbers_check.txt`.

**E′'s new statistics**
- Script: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/scripts/eprime_check.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/out/eprime_check.txt`.
- It recomputes, per cell, from the raw:
  - the part, its SE (ddof 0 and 1) and the discordant counts;
  - `betw` at q₂ and at q_law;
  - the excess, split into positive, zero, negative below 10⁻¹⁰, and clearly negative;
  - any infinite or non-positive knots.
- The derivation of the inequality and the decomposition was done by hand (Summary).

**Population**
- Script: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/scripts/pop_part.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/out/pop_part.txt`.
- For each of the 15 cells:
  - 10⁶ law samples shared by `knot_of` and `knot_of_kernel`, with the inequality asserted on every sample;
  - 10⁵ fresh H1 draws, fitted with a vectorised `one_fit` and scored with the sliding-limit knots;
  - before use, the fits and knots were validated against E's stored fits, E's knots and U's knots on E's own regenerated draws.
- The run takes about 11 s per cell [own].

**Who sets the tail**
- Script: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/scripts/tail_who.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/out/tail_who.txt`.
- It uses E's own law samples (rng [seed, 9]).

**Size by distance**
- Scripts:
  - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/scripts/size_by_distance.py`;
  - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/scripts/trend.py`;
  - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/scripts/size_excl.py`.
- Outputs: `size_by_distance.txt`, `trend.txt` and `size_excl.txt` in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/out`.
- What they compute:
  - the size per H0 job, with the null's own centre c1 [raw];
  - the size grouped by distance, both from c̄ and from c1, per pair and pooled;
  - Cochran–Armitage trend tests;
  - Fisher tests between mirror offsets and between the two pairs;
  - the pooled sizes without the shifted nulls.

**Mutation**
- Definitions: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/mut/mutants.py`, 82 single-point mutants plus the identity:
  - **E01–E22** (`e_prime_lines`): the margin fallback; the '—' row; `own1` with ≤; the SE's threshold, ddof, n, or H0 draws; the excess misordered, filtered, inverted, counted with ≥, H1 only, or as |excess|; `betw` on H0, at q_law, over all cells, or with both rejecting; n1 over all cells; the SE printed only where nonzero.
  - **F01–F13** (`f_alloc_lines`): the spacing; no mirror; no rounding; H1 rows pooled; ≤ instead of <; the CI from k1; unsorted; first pair only; the labels; empty offsets; distance to the tie; rounding to two decimals.
  - **U01–U09** (`u_job`): j₁ instead of j10; signed argmax; j10 from f₁; argmin; max|g0|; no wall; signed denominator; g0 from f₁; the pair swapped in q₂'s law.
  - **38 ported** from earlier rounds: SM02, SM04, A06, FA01, FA03–FA07, SL02, SL03, SL05, EM02–EM07, EM09–EM11, FR01–FR05, ES03–ES05, ES09, ES10, ES12–ES18.
- Runner: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/mut/run_mut.py`, with `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/mut/behav.py`. Each mutant runs against:
  - the full test file;
  - `--from-raw`, with a byte diff of the md;
  - for `u_job` mutants, a rerun of the 15 U jobs with the md regenerated and diffed.
  - (s_job mutants that survived the tests would also have had a 15 E S job rerun. Neither survived.)
- The identity passes all three checks [own].
- Results: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/mut/results_a.jsonl` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/mut/results_b.jsonl`; summary in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref75/D3bv_work/out/mut_summary.txt`.
- Outcome [own]:

  | code | mutants | killed by the tests | caught only by the md diff or the U rerun | survive both |
  |---|---|---|---|---|
  | round 73: `e_prime_lines` | 22 | 16: E01, E02, E04, E06–E08, E11, E12, E15–E22 | — | 6: E03, E05, E09, E10, E13, E14 |
  | round 73: `f_alloc_lines` | 13 | 9: F01–F03, F06, F07, F10–F13 | 3: F04, F05, F09 | 1: F08 |
  | round 73: `u_job` | 9 | 5: U02, U04, U06, U07, U08 (U04 and U08 also change 15 md lines on rerun) | — | 4: U01, U03, U05, U09 |
  | earlier rounds', ported | 38 | 7: SM02, SM04, EM10, EM11, ES05, ES12, ES13 | 24: FA01, FA03–FA07, SL02, SL03, SL05, EM02, EM03, EM07, EM09, FR02, FR04, FR05, ES04, ES09, ES10, ES14–ES18 | 7: A06, EM04–EM06, FR01, FR03, ES03 |

- **Survivors that are equivalent in general or immaterial:**
  - E03: ≤ against an interpolated quantile;
  - E05: changes the SE by 0.1%;
  - E10: t₂ ≤ 0 cannot occur;
  - E14: q_law = q₂ in the summed cells;
  - U05: |g0[j10]| = max|g0| by construction;
  - EM04: q₂ ≥ q_full always.
- **Survivors that change behaviour only on inputs this data lacks:**
  - E09: an infinite path knot;
  - E13: H0 draws between the knots; here that count is also 0;
  - F08: offsets not listed in order of distance;
  - U01 and U03: a null whose first entrant is not f₁'s;
  - U09: a third kernel leading in a law sample;
  - EM05 and EM06: q_full ≠ q_law;
  - A06, ES03 and FR03: ŵ < 0;
  - FR01: an H0 extreme of s/|ŵ|.

**Literature**
- Web searches only. projecteuclid.org is not among the session's allowed hosts, so the abstract page could not be fetched.

## Verdicts per part

| claim (D52.md line) | verdict | evidence |
|---|---|---|
| Title, (i), (ii), (iv), (v) (l. 1–41, 83–102) | **unchanged; reproduces** | The raw is identical to round 71's, and the md regenerates byte for byte [own] |
| (iii) calibrated size by distance: "nominal at a kernel's centre (0.0512 …), 0.0456 … at 0.02, below it from 0.032 on (0.0398, 0.0379, 0.0329), lowest next to the tie"; "keeps its size across the offset" (l. 54–58) | **numbers right; the gradient belongs to one null; "keeps its size" is loose** | md = own. Strong null (\|ŵ\| ≈ 1.4): 0.0517, 0.0471, 0.0421, 0.0450, 0.0383, trend p = 0.15. Weak null (\|ŵ\| ≈ 0.4): 0.0508, 0.0442, 0.0375, 0.0308, 0.0275, p = 1·10⁻⁴. At distance ≥ 0.032 the two nulls differ (p = 0.006) [own]. The (0.8, −0.4) nulls at s = 0.03 sit at c̄ + 0.0097 [raw] (D3bv-1) |
| (iii) de facto one-sided, 75–100% (l. 58–60) | **holds** (unchanged) | [md] |
| (iii) the minimum's part off the tie: "a loss where other kernels' knots share the null's tail (c̄ ≤ 0.032: 0.025–0.12, standard errors 0.013–0.021 per cell)" (l. 69–71) | **holds; the range is noisy** | Other kernels set 34–82% of the tail [own]. Population 0.034–0.087 [own]. The SEs are right (χ² = 8.0 on 9 df) [own] (D3bv-2) |
| (iii) at the tie: the same quantile; never later; cannot lose; the part is the share between the knots; gains nothing (0 of 2400) (l. 71–75) | **holds; the logic is now right** | The derivation in the Summary. 0 of 2400 [own = md]. Population: equal quantiles with 10⁶ samples, and 0 of 600 000 draws between the knots [own] |
| Status: (iii)'s scope (l. 145–149) | **holds** | "The α allocation" now appears once [own] |
| Status: new sources (l. 156–161) | **partly holds** | Dacunha-Castelle–Gassiat and Azaïs–Gassiat–Mercadier treat population mixtures. Azaïs–Gassiat–Mercadier 2006's local-power result is "close to the level" [snippet]. Amar–Weiss's δ⁻⁴ is unverified [memory]. El Korso et al. use a first-order Taylor model and tie the limit to a CRB [snippet] (D3bv-3) |
| New (l. 173–180) | **holds** (unchanged) | "A loss away from the tie, none at it": population [own] |
| Numerics: F by distance (l. 259–260) | **matches the md** | [md = own]. The caveat of D3bv-1 applies |
| Numerics: E′'s part (l. 260–264) | **matches the md; one parenthesis is ambiguous** | "3293 of 21000" counts all 15 cells; at the tie cells the count is 2189 of 8400 [own] (D3bv-7) |

## Status of the round-73 issues

| D3bt id | status | what remains |
|---|---|---|
| D3bt-1 (minor) | **fixed** | E′ now counts the H1 draws between the knots (0 of 2400) and the excess (3293 of 21000, at most 6.5·10⁻¹²) [md, own]. (iii) and the Numerics now state the inequality's direction correctly, and the population confirms the zero (0 of 600 000) [own]. Loose ends: the "(0.000)" SEs (D3bv-5); the scope of "3293 of 21000" (D3bv-7); round 71's log entry is not corrected in place (D3bv-7) |
| D3bt-2 (minor) | **partly fixed** | Each cell's SE is now printed (0.013–0.021 [md]), and the SEs are right (χ² = 8.0 on 9 df [own]). (iii) still quotes the noisy extremes "0.025–0.12" [md] against the population's 0.034–0.087 [own] (D3bv-2) |
| D3bt-3 (minor) | **partly fixed** | Now killed by the tests: the targeted SM02, SM04, U71-04, U71-06, ES05, ES12 and ES13, and also EM10 and EM11. U71-05 is now equivalent by construction [own]. The `f_alloc_lines`, `smith_line` and `f_lines` tests that D3bt-3 proposed were not added: 24 mutants are still caught only by the md diff, and A06, ES03, FR01 and FR03 still survive both checks [own] (D3bv-4) |
| D3bt-4 (minor) | **fixed as proposed; the new statement over-reads** | F prints the size by distance with mirrors pooled (0.0512 … 0.0329 [md]), and the mirrors agree [own]. But the pooled gradient is the weak null's, and the shifted nulls are not at their mirrored positions (D3bv-1) |
| D3bt-5 (cosmetic) | **mostly fixed** | Done: the repeated "α allocation", round 67's numbers removed, j₁ taken from f₀ (the same kernel in 15 of 15 cells [own]), '—' printed without U, the four sources added. Remaining: the fit of the sources (D3bv-3); E′'s header and `s_job` still use f₁'s first entrant (D3bv-7) |

## Issues

| id | sev. | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bv-1 | minor | **The size-by-distance statement averages two nulls that behave differently, and some of its "mirrors" are not mirrors** (l. 54–58; Numerics l. 259–260; the md's F line; round 73's Lesson).<br>**The two nulls.** F's H0 nulls are the one-bump fits of the two pairs: \|ŵ\| = 1.385–1.399 and 0.403–0.431 [raw], so s/\|ŵ\| = 0.0021–0.0226 and 0.0074–0.0787 [own].<br>**Sizes at distances 0, 0.02, 0.032, 0.045, 0.049** [own]:<br>– strong null: 0.0517, 0.0471, 0.0421, 0.0450, 0.0383 (trend z = −1.45, p = 0.15);<br>– weak null: 0.0508, 0.0442, 0.0375, 0.0308, 0.0275 (z = −3.83, p = 1·10⁻⁴);<br>– at distance ≥ 0.032 they differ: 0.0425 [0.0375, 0.0479] against 0.0328 [0.0285, 0.0377], Fisher p = 0.006.<br>**E agrees for the strong null.** E's calibrated sizes for (0.8, 0.6) are 0.050, 0.053, 0.049, 0.042, 0.035 at c̄ = 0 … 0.049 [own, from md]: nominal through c̄ = 0.032, dipping only next to the tie.<br>**The shifted nulls.** At s = 0.03 the (0.8, −0.4) nulls sit at c̄ + 0.0097 [raw], as l. 254–255 itself says, yet F pools them by c̄.<br>– c̄ = 0.045 and 0.049 put the null past the tie, at 0.0547 and 0.0587 (distances 0.045 and 0.041); c̄ = 0.055 puts it at distance 0.035 [own].<br>– Binned by the null's own centre, the lowest bin is [0.035, 0.0425): 0.0175 (21 of 1200), made up of three of these jobs, not "next to the tie" [own].<br>– Without these nulls the pooled pattern persists: 0.0495, 0.0473, 0.0425, 0.0390, 0.0360 [own]. | **(iii) l. 54–58:** "holds its level"; then, by distance, "for the strong null (\|ŵ\| ≈ 1.4) 0.052 and 0.047, then 0.038–0.045 (0.0425 [0.0375, 0.0479] pooled); for the weak null (\|ŵ\| ≈ 0.4) falling to 0.028 next to the tie".<br>**F (`f_alloc_lines`):** print the by-distance line per pair, and group by the null's own centre `r['c1']`, or leave the shifted nulls out of the mirror pooling.<br>**Numerics l. 259–260:** the same.<br>**REVIEW_LOG round 73 Lesson:** qualify "falls smoothly with the distance". |
| D3bv-2 | minor | **(iii) still quotes the extremes of nine noisy estimates as the off-tie loss** (l. 71).<br>– The text says "0.025–0.12" [md].<br>– The population values are 0.0833, 0.0711, 0.0816 at c̄ = 0; 0.0792, 0.0678, 0.0872 at c̄ = 0.02; and 0.0344, 0.0340, 0.0383 at c̄ = 0.032, with SE ≈ 0.001 [own]. That is 0.034–0.087.<br>– The printed SEs are right (χ² = 8.0 on 9 df) [own], so this is about what (iii) quotes, not about the numbers. | **(iii) l. 71:** "≈ 0.08 at c̄ ≤ 0.02 and ≈ 0.035 at 0.032 (per-cell SE 0.013–0.021)".<br>**E′ (optional):** a population line from 10⁵ fresh H1 draws per cell with the vectorised fit. It costs about 11 s per cell [own]. |
| D3bv-3 | minor | **The new sources are cited for results they do not, or do not visibly, contain** (Status l. 156–161).<br>– **Azaïs–Gassiat–Mercadier 2006's** local-power result is that, in the contamination model, the asymptotic local power under contiguous alternatives "may be arbitrarily close to the asymptotic level" for large parameter sets, and equals it for an unbounded Gaussian mean [snippet]. It is not a source for "the mixture form of Δ₅₀ ∝ N^{−1/4}". That order is classical: Chen 1995; Chen–Chen–Kalbfleisch 2001 [memory].<br>– **Dacunha-Castelle–Gassiat (1997, 1999)** treat population mixtures, that is, i.i.d. draws from a mixture density [snippet]. Their weights are positive, so the LRT limit is the squared positive part of the supremum [memory]. The paper's model, a signed superposition in additive Gaussian noise, is the regression analogue (Davies 1977/1987; Liu–Shao 2003; Gassiat 2002 [memory]).<br>– **Amar–Weiss's "SNR ∝ δ⁻⁴"** is tagged [snippets], but neither D3bt nor I found a snippet stating it [own search]. It is [memory].<br>– **El Korso et al.** use a GLRT on a first-order Taylor model and tie the limit to a CRB [snippet]. That makes them kin to Liu–Nehorai and Smith, with the same amplitude caveat, not to the Δ² ∝ s detection law. They are also cited without a year. | **Status l. 156–161:**<br>– present Dacunha-Castelle–Gassiat as "the mixture analogue";<br>– take N^{−1/4} from Chen 1995 or Chen–Chen–Kalbfleisch 2001;<br>– cite Azaïs–Gassiat–Mercadier 2006 for the cost of searching a large parameter set, which is the analogue of the dictionary's minimum in (iii);<br>– retag Amar–Weiss [memory];<br>– move El Korso et al. next to Liu–Nehorai, with a year. |
| D3bv-4 | minor | **Round-73 code is well tested; older reporting code still is not.**<br>– For round 73's own code: 30 of 44 mutants are killed by the tests and 3 more (F04, F05, F09) by the md diff; the 11 survivors are equivalent on this data or immaterial [own].<br>– Still caught only by the md diff, as in round 73: FA01, FA03–FA07, SL02, SL03, SL05, EM02, EM03, EM07, EM09, FR02, FR04, FR05, ES04, ES09, ES10, ES14–ES18 (24 mutants) [own].<br>– Still surviving both checks: A06, ES03 and FR03 (ŵ < 0) and FR01 (an H0 extreme) [own].<br>– The tests D3bt-3 proposed for `f_alloc_lines`, `smith_line` and `f_lines` were not added. | **Tests:**<br>– `f_alloc_lines` on rows that differ by pair and offset, with H1 rows, qA ≠ qB, a `q_law` and a ŵ < 0 row, asserting the α = 0.01 figures;<br>– `smith_line` with a non-integer power and the unpatched `crb_sep`;<br>– `f_lines`' range and equal-size lines asserted;<br>– E′'s margin ratio with mu0 ≠ mu_true (EM02);<br>– the by-distance test with H1 rows and a p exactly at α (F04, F05, F09). |
| D3bv-5 | cosmetic | **E′ prints "(0.000)" as the standard error at the six tie cells** (md l. 210).<br>– This is a degenerate SE computed from zero discordant draws.<br>– The 95% one-sided bound is 0.0075 per cell (0 of 400) and 0.00125 pooled (0 of 2400) [own]. | **E′:** print the count and the bound instead. |
| D3bv-6 | cosmetic | **The excess statistic's filter silently drops non-finite knots** (exp l. 839).<br>– An infinite path knot next to a finite own knot is the extreme violation of "never later".<br>– The filter would drop such a draw and lower the denominator instead of counting it.<br>– E09 survives only because all 21000 knots are finite [own]. | **E′:** count dt = ∞ with a finite t₂ as an exceedance, and assert that the denominator equals the number of draws. |
| D3bv-7 | cosmetic | **Wording and log loose ends.**<br>– **Numerics l. 262–263:** "3293 of 21000" sits inside the c̄ ≥ 0.045 parenthesis but counts all 15 cells. At the six tie cells it is 2189 of 8400 [own]. The informative denominator is 7071, the draws where the two knots coincide [own].<br>– **(iii) l. 55:** "keeps its size" sits next to sizes of 0.0329–0.0512.<br>– **(iii) l. 74–75 and E′:** "no H1 draw falls between them". It is the threshold that falls between a draw's two knots.<br>– **E′'s header (exp l. 796–799)** still normalises by "the noiseless first entrant". `s_job`'s mu_u uses f₁'s entrant while `u_job` uses f₀'s. They are the same kernel in 15 of 15 cells [own], but the mismatch is latent.<br>– **REVIEW_LOG l. 2812 (round 71)** still says "as the inequality bounding it … requires". The correction lives only in round 73's D3bt-5 fix column.<br>– **Status l. 156–161** has three nested dash pairs. | Edit the text. Add an erratum in place at REVIEW_LOG l. 2812, as round 71 did for round 64's entry. |

## Sources

- **Read in the export:**
  - in full: D52.md, close_features.md (E, E′, F), `test_close_features.py` and REVIEW_round73_D3bt.md;
  - the Round 71 and Round 73 entries of REVIEW_LOG.md;
  - the functions of `exp_close_features.py` listed in Scope;
  - `rich_transit.py` (`Bumps`);
  - the diffs 1f51d48..3dc9443.
- **Web, snippets only** (projecteuclid.org is not among the session's allowed hosts):
  - [Dacunha-Castelle–Gassiat 1999, Ann. Statist. 27(4) 1178–1209](https://projecteuclid.org/journals/annals-of-statistics/volume-27/issue-4/Testing-the-order-of-a-model-using-locally-conic-parametrization/10.1214/aos/1017938921.short), and [their 1997 ESAIM P&S paper](https://numdam.org/item/PS_1997__1__285_0);
  - [Azaïs–Gassiat–Mercadier 2006, Bernoulli 12(5) 775–799](https://projecteuclid.org/euclid.bj/1161614946), with the [preprint listing](https://www.imo.universite-paris-saclay.fr/~biblio/ppo/2003/ppo2003-57.html): the abstract's local-power statements;
  - [Azaïs–Gassiat–Mercadier 2009, ESAIM P&S](https://www.numdam.org/item/PS_2009__13__301_0) (title only);
  - Amar–Weiss, IEEE TSP 56(11) 5309–5318 (2008): [journal listing](https://datalearner.com:8443/academic/journal-papers/1053-587X/volumes-and-issues/4652320), title only, and [dblp, AVSS 2007](https://dblp1.uni-trier.de/rec/conf/avss/AmarW07.html). The δ⁻⁴ law is not confirmed.
  - El Korso–Boyer–Renaux–Marcos: [hal-00771379](https://hal-supelec.archives-ouvertes.fr/hal-00771379) (GLRT, SRL tied to a CRB) and [hal-00771402](https://hal-centralesupelec.archives-ouvertes.fr/hal-00771402) (first-order Taylor model).
- **[memory]:** Chen 1995 (the n^{−1/4} rate); Chen–Chen–Kalbfleisch 2001; Davies 1977/1987; Liu–Shao 2003; Gassiat 2002; the positive-part form of the mixture LRT; Amar–Weiss's δ⁻⁴.

## Lesson

Symmetry tells you what may be pooled, not what must be. Round 73 pooled each offset with its mirror, which is right for the grid. But it also pooled two nulls of different amplitude, whose sizes differ (p = 0.006), and it labelled some nulls by an offset they do not sit at. The resulting "smooth fall with the distance from the nearest kernel" is the weak null's alone; the strong null holds its level until next to the tie. A pooled rate needs its strata before its gradient is read.

Two smaller points follow the same idea of matching a statement to its evidence. A count of zero deserves an upper bound, not a zero standard error. And a citation should carry the result its source is known for. Azaïs–Gassiat–Mercadier's local-power theorem is about the cost of searching a large parameter set, which is this paper's "dictionary's minimum", not about the N^{−1/4} order.
