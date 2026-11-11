# Referee report D3cc, round 82: Proposition D‴.52 as revised after round 78

## Scope

- **Commit reviewed:** 5201021 on claude/festive-bardeen-eopej2 ("D‴.52 round 78: the size falls toward the tie for both nulls").
  - I exported it with `git archive 5201021 | tar -x` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/repo`.
  - All other work is in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work`, in its scripts/, out/, logs/ and mut/ folders.
  - HEAD moved to 3b6aac7 while I worked, through other sessions' D‴.54, D‴.58 and D‴.59 commits. `git diff 5201021 3b6aac7` changes none of the files reviewed and no D‴.52 line of the README. The Round 73, 75 and 78 log entries are unchanged, apart from one blank line added after Round 78.
- **Files read:**
  - in full: D52.md; `exp_close_features.py`; close_features.md; `test_close_features.py` (53 tests); REVIEW_round78_D3by.md; the Round 73, 75 and 78 entries of REVIEW_LOG.md, with their errata;
  - `git diff 6764e41 5201021` for the code, the md, D52.md, the tests and the drafts README (the raw is unchanged since 6764e41);
  - in part: `rich_transit.Bumps` and the `sliding_iss` docstring; `exp_mirror_transit`'s NS_C, SIG_C, D_C; `exp_cohesion_leak.load/save`; D3by's mutant definitions, which I ported.
- **Settings:**
  - the repository's venv (Python 3.13.16, numpy 2.5.4, scipy 1.18.1);
  - every process run with `nice -n 10`, OMP_NUM_THREADS=1 and PYTHONDONTWRITEBYTECODE=1;
  - pytest as `-q -p no:cacheprovider`;
  - `CF_OUT` always pointed to a scratch copy of the raw;
  - at most two heavy processes of mine at a time.
- **I modified nothing in the repository** and ran no git command that writes. `git status --porcelain` was empty before and after.
- **All my processes have exited** (checked with `ps`). The Python and `tail` processes still running belong to another session.
- **Tags:** (md) is close_features.md, (raw) is close_features_raw.npz, (mine) is my own run or derivation. For the literature, [snippet] means a search result and [memory] my own knowledge. arxiv.org could not be reached (ENOTFOUND).

## Summary

**What checks out**

- **Reproduction:**
  - `--from-raw` regenerates the md byte for byte (`cmp`; md5 a476cc17…) (mine).
  - The 53 tests pass in 24.3 s (mine).
- **The round's new md lines recompute independently from the raw:**
  - the bins, the nulls' distances and the Clopper–Pearson intervals;
  - the Cochran–Armitage z values −1.45, −2.97, −3.28 and −3.24 (md). Mine are −1.4532, −2.9660, −3.2791 and −3.2374;
  - by direct maximum likelihood: the interaction z +1.418, the likelihood ratio 2.0016 with p 0.1571, and the common slope z −4.3705 (mine; md +1.42, 0.157, −4.37).
  - `trend_z` is the standard Cochran–Armitage statistic with binomial variance. `logit_fit` is Newton/IRLS with the right deviance. erfc(√(LR/2)) is the χ²₁ tail (mine).
- **E's and F's strong nulls are comparable:**
  - The null is the same: identical at s = 0.01 and 0.03, and at s = 0.003 within 6·10⁻⁴ in ŵ (E fits at Δ = 0.0175, F at 0.014) (raw).
  - The statistic is the same: `knot_walls`' minimum is `knot_of`.
  - The calibration is the same: `knot_law` with M_E = 2000 at the path's own fit, rejecting at p < 0.05. That rule equals the exact Monte Carlo rule (1 + #)/(M + 1) ≤ α [memory].
  - The seeds are independent, and the distances are those of the truth's one-bump fit (mine).
- **E′ holds:**
  - the pooled parts −0.1042 (0.0119), −0.0850 (0.0109) and −0.0292 (0.0081), and the contrast z −1.19 (md = mine; −1.14 to −1.17 with D3by's laws' error added);
  - the bounds 0.007461 / 0.009180 per cell and 0.001247 / 0.001536 for the mean share are the exact Clopper–Pearson ends (mine);
  - the AM–GM step Π(1 − p_i)⁴⁰⁰ ≤ (1 − p̄)²⁴⁰⁰ holds (on 10⁵ random vectors the largest log gap is −1.3·10⁻⁶) (mine);
  - "≈ 0.08, 0.08, 0.035" are D3by's per-cell means, 0.0787, 0.0777 and 0.0353, and indeed no repository script computes them (mine).
- **The numbers on the changed lines trace to the md.** All the decimal numbers on D52.md's 100 changed lines are in the md, verbatim or rounded. The exceptions are years, HAL ids, three derived values (90.7%, 100%, 6.5·10⁻¹²) and two |ŵ| ranges (D3cc-4) (mine).
- **Sources:** the squared positive part, Liu–Shao and Gassiat placed with the mixtures, and Davies alone as the regression analogue now fit [snippet/memory].
- **Fixed:** D3by-4, -5, -6 and -7.

**What does not (majors first)**

1. **D3cc-1 (major).** The round's central claim reads a trend in the distance into nulls that disagree at one distance. The claim is that the calibrated knot "grows conservative toward the tie for both nulls … the two slopes do not differ" (l. 54–60, 270–277).
   - **The spread is beyond binomial, and it replicates.** Nulls at the same distance disagree far beyond binomial noise: in the md's own bins, χ² = 39.6 on 8 df and 21.4 on 5 (mine). E and F reproduce each other's low cells (raw). The logistic model's Pearson dispersion is 2.1 (mine).
   - **With the null as the unit, the strong null's trend is not established:** permutation p = 0.054 and quasi-binomial z = −2.1 (mine).
   - **What orders the nulls is the runner-up's margin at the null over its noise, μ₀/sd(μ̂).**
     - The 12 nulls with 1.2 < μ₀/sd(μ̂) < 5 have size 0.022. The other 51 have 0.048, with no residual trend in the distance (mine).
     - The band was read off E's two low cells. It predicts F's strong nulls, the weak null out of sample, and two pre-registered new runs far from the tie: 0.027 and 0.0235, against the ≈ 0.046 the distance reading predicts (mine).
   - **It is the calibration's plug-in, not the knot.** On the same draws, the knot against its law at the truth stays at 0.045 and 0.044 (mine).
   - **"Do not differ" asserts equality from p = 0.16.** That is the converse of round 78's own Lesson.
2. **D3cc-2 (minor).** 13 of my 58 mutants of the round's code pass all 53 tests and change published lines, the central ones among them:
   - dropping E's nulls from the size lines;
   - any number of the logistic line;
   - the trend's weights (mine).

   Of the 28 earlier mutants named in D3bv-4, 11 are killed now (mine; D3by found 6).
3. **D3cc-3 (minor).** "Chen 1995's n^{−1/4} is the rate …, corrected by Heinrich–Kahn 2018" still asserts n^{−1/4}. Heinrich–Kahn give n^{−1/(4(m−m₀)+2)}, which is n^{−1/6} for one extra component [snippet].
4. **D3cc-4 to D3cc-7 (cosmetic):**
   - two |ŵ| ranges the md no longer prints;
   - a referee's estimate hard-coded in the md;
   - the README's two contradicting readings;
   - code and wording loose ends.

## What I checked and how

All paths below are written out in full. They are all under `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work`.

**Reproduction**
- `--from-raw` with `CF_OUT=/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/fromraw`, then `cmp` against the commit's md.
- Tests: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/pytest_base.txt`.

**Size by distance**
- Every H0 null of E and F, with its distance, ŵ, k and n: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/nulls.py`. Output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/nulls.json` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/nulls.txt`.
- Independent trend and maximum-likelihood logistic fits, within-bin heterogeneity, Pearson dispersion, quasi-binomial z, per-(s, bin) sizes: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/trend_check.py`, output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/trend_check.txt`.
- The trend with nulls as units (20 000 permutations of the nulls' distances): `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/perm_trend.py`, output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/perm_trend.txt`.

**The margin band**
- Each null's runner-up margin μ₀ at f₀, sd(μ̂) from 20 000 draws, and the logistic models in μ₀/sd(μ̂): `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/margin_size.py`. Output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/margin_size.txt` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/nulls_margin.json`.
  - My μ₀/sd values match E′'s, for example 2.11 against 2.08 and 0.74 against 0.73 (md).
- Out-of-sample checks, robustness to the band's edges, and a bootstrap over nulls: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/band_check.py`, output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/band_check.txt`.
- The knot against its law at the truth beside the calibrated knot, on the same draws: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/band_oracle.py`, output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/band_oracle.txt`.
- New H0 runs, calibrated exactly as E is:
  - Scripts: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/band_predict.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/band_predict2.py`. The second adds the law at the truth.
  - The predictions were written before the results to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/band_predict_PREREGISTERED.txt`.
  - Results: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/band_predict.jsonl` (its first line is a 20-draw timing run) and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/band_predict2.jsonl`.

**E′ and the numbers**
- E′'s parts, standard errors, contrast, population comparison, bounds and the AM–GM check: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/eprime_check.py`, output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/eprime_check.txt`.
- Every number on D52.md's changed lines against the md: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/numbers_check.py`, output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/numbers_check.txt`. I then checked each number on the round's new lines by hand against its md source.
- The first-entrant claim: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/scripts/first_entrants.py`, output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/first_entrants.txt`.

**Mutation**
- My mutants: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/mut/mutants.py`, 58 single-point mutants plus the identity:
  - T01–T08 (`trend_z`), L01–L07 (`logit_fit`), N01–N06 (`size_nulls`), S01–S16 (`size_by_distance`);
  - P01–P03 (`part_se`), C01–C05 (the pooled contrast and SEs), B01–B04 (the pooled bound);
  - F01–F06 (`f_alloc_lines`' per-offset and components lines), R01–R03 (passing E_ through).
- D3by's mutants that still apply, ported unchanged: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/mut/mutants2.py` (48; FA03 at its `f_alloc_lines` occurrence). It reads the copied definitions in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/mut/d3by_mutants.py`.
- Runner: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/mut/run_mut.py`. Each mutant runs against the test file with `-x`, then `--from-raw` with a line diff against the commit's md.
- Results: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/mut/results_a.jsonl`, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/mut/results_b.jsonl` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/mut/results_c.jsonl`.
- Summaries: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/mut_summary_a.txt` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/mut_summary_c.txt`.
- The identity passes every check.

## Verdicts per part

| claim (D52.md lines) | verdict | evidence |
|---|---|---|
| Title: the score at the bound with the sign known, ≈ 1.09 (1.04–1.18) without it (l. 1–4) | **holds** (unchanged) | \|z\| 1.085 [1.039, 1.135] and 1.094 [1.059, 1.179] (md) |
| (i), (ii) (l. 8–41) | **holds** (unchanged; sampled) | A's ratios 1.015 → 0.9995; B's 91.34 against 91.23 (md) |
| (iii): the knot's law, LARS, the null quantile across a cell (l. 43–54) | **holds** (unchanged) | sd(μ̂)/s 0.148–0.170 is E's (md). For the weak null it is 0.50–0.58 (mine) |
| (iii): the calibrated knot "holds its level … grows conservative toward the tie for both nulls … the two slopes do not differ" (l. 54–60) | **numbers right; the reading is not supported** (D3cc-1) | All numbers md = mine. Dispersion 2.1. Null-level p 0.054 for the strong null. The size is set by μ₀/sd(μ̂) and by the calibration's plug-in (mine) |
| (iii): one-sided by the margins (l. 60–62) | **holds** (unchanged) | 390/430 and 328/364 (md) |
| (iii): power comparisons (l. 62–69) | **holds** (unchanged; spot-checked) | 1.056/1.003 (0 of 24 ahead); 1.045 and 1.032 (4 of 24); 1.113/1.085; 1.064/1.019 (2 of 24); 0.987 (17 of 24); 2.093/1.257 (md) |
| (iii): the pooled parts, their contrast, the population values "not reproduced" (l. 70–75) | **holds** | md = mine; z −1.186. D3by's means are 0.0787, 0.0777, 0.0353, and no script computes them (mine) |
| (iii): none at the tie; the one- and two-sided bounds (l. 75–80) | **holds** | 0 of 2400 discordant (md = mine). The bounds are exact, and AM–GM holds (mine) |
| (iii): opposite-sign bias; the null moving with Δ (l. 80–86) | **holds** (unchanged) | Least power 0.000–0.025, below α in 24 of 24 (md) |
| (iv), (v), Proof (l. 88–146) | **not re-examined** (unchanged) | (v)'s 24Δ³ and 12Δ² as in D (md) |
| Status: the literature (l. 156–183) | **mostly holds** | Squared positive part, Liu–Shao/Gassiat, Davies and contiguous alternatives fit [snippet/memory]. The Chen sentence is wrong (D3cc-3) |
| Numerics A–D and Smith (l. 193–224, 283–284) | **holds** | The labels are now right: 474Δ^{2.99} unknown (0.8, 0.6), 109Δ^{1.98} known (0.8, −0.4) (md) |
| Numerics E (l. 224–243) | **holds** (unchanged) | Its 0.019 and 0.022 next to the tie are exactly the band's two E cells (md) |
| Numerics F, the size by distance (l. 266–277) | **numbers right except two orphaned ranges; reading as in D3cc-1** | md = mine. The \|ŵ\| ranges are not in the md (D3cc-4) |
| Numerics E′ (l. 277–283) | **holds** | md = mine; a hard-coded literal (D3cc-5) |

## Status of the round-78 issues

| D3by id | status | notes |
|---|---|---|
| D3by-1 (minor) | **done as asked, but superseded (D3cc-1)** | E's 15 000 strong-null draws are pooled with F's, with per-null trends and a logistic interaction (md = mine). The reading this produced rests on binomial statistics over nulls that are not exchangeable. The fall is a band in μ₀/sd(μ̂), and it comes from the calibration. "Do not differ" is overstated. |
| D3by-2 (minor) | **partly fixed** | **This round's code:** the tests kill 40 of my 58 mutants. 13 change published lines unseen (D3cc-2) (mine).<br>**D3by's round-75 mutants:** 17 of the 20 that still apply are killed now, among them X01–X06, Q01–Q05 and N01–N02. B02, D01 and D02 are equivalent on these data.<br>**The 28 named in D3bv-4:** 11 are killed (6 before).<br>– 14 still change only the md: FA03, FA05, EM03, EM07, EM09, FR04, ES04, ES09, ES10 and ES14–ES18.<br>– FR01, FR03 and ES03 survive both checks.<br>**Docstrings:** the f_alloc test still credits FA03 and FA05 (mine). |
| D3by-3 (minor) | **mostly fixed** | Done: the squared positive part; Liu–Shao and Gassiat with the mixtures; Davies alone; the contiguous alternatives. Left: the Chen sentence (D3cc-3). |
| D3by-4 (cosmetic) | **fixed** | The bins print as intervals with their nulls' distances (md). |
| D3by-5 (cosmetic) | **fixed** | The contrast z = −1.19 is printed. The population values are marked as not reproduced. The size of the laws' error is quoted, but hard-coded (D3cc-5). |
| D3by-6 (cosmetic) | **fixed** | One- and two-sided bounds, per cell and for the mean share, all exact (mine). |
| D3by-7 (cosmetic) | **fixed** | Smith's labels, "by the label", the first-entrant note (true for E's pair: D3cc-7), the dash, and El Korso's year marked. |

## Issues

| id | severity | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cc-1 | **major** | **The calibrated knot's size is said to fall with the distance, "for both nulls", with slopes that "do not differ".** Where: l. 54–60, l. 270–277, md l. 293–297, the errata at rounds 73 and 75, the Round-78 Lesson and the README. That reading is not what the data show.<br>**(a) Extra-binomial spread, replicated.**<br>– Strong null, E and F pooled: the nulls' sizes give χ² = 94.7 on 38 df (mine).<br>– Within the md's own bins, [0.0385, 0.047) gives 39.6 on 8 df (p < 10⁻⁵) and [0.047, 0.05] gives 21.4 on 5 (p = 0.0007). The three far bins are binomial: p = 0.30, 0.48, 0.27 (mine).<br>– E and F agree on the outliers. At c̄ 0.049, s 0.01: 19/1000 and 6/400. At c̄ 0.045, s 0.03: 22/1000 and 10/400, with the mirror 6/400 (raw).<br>– The logistic model's deviance is 133.6 on 59 df (Pearson dispersion 2.1). For the same reason the tie bin's [0.0305, 0.0420] is too narrow (mine).<br>**(b) The z values treat the draws as independent across nulls.** With the null as the unit (mine):<br>– strong null: permutation p = 0.054 (E alone 0.11, F 0.32); quasi-binomial z = −2.1 (E −1.7, F −1.0); a bootstrap over nulls z = −2.2;<br>– weak null: p = 0.014;<br>– common slope z = −3.0; interaction z = +0.98, F-test p = 0.33.<br>The md prints −3.28, −2.97, −1.45, −3.24 and −4.37 (md). All these groupings were also chosen after the data were seen.<br>**(c) The size is ordered by the runner-up's margin at the null over its noise, μ₀/sd(μ̂)** (E′ prints it for E):<br>– 1.2 < μ₀/sd(μ̂) < 5: 12 nulls at 0.022 [0.018, 0.026] (strong 0.020, weak 0.025). All the other nulls: 0.048 [0.045, 0.050] (mine).<br>– Outside the band there is no trend in the distance: z = −0.36 over 34 strong nulls, −1.06 over 17 weak (mine).<br>– A two-parameter band model fits (deviance 68.9 on 60 df, p = 0.20; the distance adds z = −0.96). The distance model does not (135.6 on 60) (mine).<br>– The band was read off E's two low cells (μ/sd 2.08 and 3.69 give 0.019 and 0.022 (md)). It predicts, out of sample (mine):<br>  – F's strong nulls: 0.018 against 0.049 (z = −4.7);<br>  – the weak null: 0.025 against 0.043 (z = −4.4);<br>  – two pre-registered new runs far from the tie: s = 0.1 at distance 0.032 (μ₀/sd 3.8) gives 0.027 [0.022, 0.033], and s = 0.06 at 0.040 (3.5) gives 0.0235 [0.017, 0.031]. The distance reading predicts ≈ 0.046 there;<br>  – controls at μ₀/sd 10.2 and 6.4: 0.0485 and 0.040.<br>**(d) The cause is the plug-in of the fitted null.** On the same in-band draws the knot against its law at the truth has 0.045 (strong) and 0.044 (weak). The md's per-cell columns show it, e.g. 0.068 against 0.015 (raw, mine). In the new runs it is 0.036 and 0.042, against 0.027 and 0.024 (mine).<br>**(e) "Toward the tie" is where the band falls when s ≤ 0.03, not a law in the distance.**<br>– At s = 0.003 the strong null has 0.043, 0.048, 0.042, 0.051 and 0.047 across the bins.<br>– At s = 0.03 it is 0.021 at 0.045 and back at 0.043 at 0.049, where μ₀/sd is 0.74 (raw).<br>**(f) "The two slopes do not differ" (p = 0.16) asserts equality from a test that failed to detect a difference.** That is the converse of the round-78 Lesson.<br>– The strong slope is half the weak's: −0.058 against −0.114 logit per 0.01.<br>– The difference's 95% interval is [−0.021, +0.134] (mine).<br>**(g) The Round-73 erratum's "F's strong null alone … lacked the power" is wrong:** F shows the band at z = −4.7 (mine). | **md:**<br>– print each H0 null's size with its s, its distance and μ₀/sd(μ̂). Take μ₀ at f₀ and sd(μ̂) by the delta method or a short Monte Carlo at f₀; no new paths are needed;<br>– print the knot against its law at the truth on the same draws beside it;<br>– test with the null as the unit (a permutation, or quasi-binomial with the dispersion printed).<br>**(iii) and the Numerics:** "holds its level (0.049, 0.043) except where the runner-up's margin at the null is a few times its noise (1.2 < μ₀/sd(μ̂) < 5). There the plug-in of the fitted centre makes it conservative (0.020, 0.025; against its law at the truth 0.045, 0.044). At these noises that band lies next to the tie for the strong null and 0.03–0.05 from a kernel for the weak one. Outside it there is no trend in the distance."<br>**Also:** drop "do not differ" (say "not resolved"), and make the same change in the errata, the Lesson and the README. |
| D3cc-2 | minor | **The tests miss the round's central lines.** 13 of my 58 mutants pass all 53 tests and change published lines (mine):<br>– **R01, R02:** `report` or `f_lines` drop E_. E's lines and the pooled line vanish, and the logistic model reverts to F (md l. 294–297).<br>– **S06–S10, S15:** each number of the logistic line can be wrong: +1.42 → +0.19 (the pair term printed), p → 0.045 (erfc(√LR)), p → 0.000 (the deviance itself), −1.42 under the same label (the coding flipped), −4.37 → −3.23 (the weak slope printed), and the wrong reduced model. The test asserts only the line's first words.<br>– **T01:** an unweighted x̄ in `trend_z` changes the pooled z from −3.28 to −3.82. The closed-form test uses equal group sizes.<br>– **L03, L06:** the deviance without its failures term or without its factor 2 changes p from 0.157 to 0.166 or 0.317. In the saturated test model the deviance is 0 either way.<br>– **F02, F06:** `<=` at 0.01 or at α in the per-offset line (e.g. 0.0512 → 0.0533). The raw has p-values exactly at α, but no test row does.<br>**The earlier mutants still caught only by the md diff:** FA03, FA05, EM03, EM07, EM09, FR04 and the ES group (D3by-2). | **`test_close_features.py`:**<br>– a report-level test with E_;<br>– the logistic line's three numbers against an independent maximum-likelihood fit (scipy) on a non-saturated design;<br>– `trend_z` with unequal n;<br>– per-offset rows with p = 0.05 and p = 0.01, and an H1 row at the second separation;<br>– the equal-size numbers through `f_lines` with both pairs.<br>**Docstrings:** credit only the mutants each test kills. |
| D3cc-3 | minor | **The Chen sentence (l. 170) still asserts n^{−1/4}.** "Chen 1995's n^{−1/4} is the rate of estimating the mixing distribution …, corrected by Heinrich–Kahn 2018" asserts the rate that was corrected. Heinrich–Kahn give the local minimax rate n^{−1/(4(m−m₀)+2)}, which is n^{−1/6} for one extra component, and call Chen's n^{−1/4} a mistake [snippet]. | **Status l. 170:** "Chen 1995 claimed n^{−1/4}; Heinrich–Kahn 2018 correct the local minimax rate to n^{−1/(4(m−m₀)+2)}, n^{−1/6} for one extra component — (v)'s s ∝ Δ³ analogue [snippets]". |
| D3cc-4 | cosmetic | **Two ranges no longer trace to the md** (l. 274, 276). "\|ŵ\| 1.385–1.399" and "0.403–0.431" were printed by round 75's size line, which `size_by_distance` replaced. They are right: 1.3852–1.3988 and 0.4027–0.4306 (mine). | **md (`size_by_distance`):** print \|ŵ\| per pair again.<br>**D52:** or drop the ranges. |
| D3cc-5 | cosmetic | **A referee's estimate is hard-coded in the md** (md l. 210; script l. 870–871). "0.003–0.006 per cell by referee round 78's estimate" is a literal in a file whose footer says its numbers are produced by the script. D52 l. 279 does not mark it as not reproduced, while l. 75 marks the population values so. | **D52 l. 279:** mark it as at l. 75.<br>**Or the script:** compute the laws' error from disjoint subsets of the shared law samples. |
| D3cc-6 | cosmetic | **The README contradicts itself.** Its D52 line keeps round 75's "the strong null holds its level …; the weak one falls" unmarked, beside round 78's "both nulls … slopes not different". Neither is the data's reading (D3cc-1). | **README:** mark the round-75 parenthesis as withdrawn, and state D3cc-1's reading. |
| D3cc-7 | cosmetic | **Loose ends:**<br>– `size_by_distance` builds `sets` and `lab` but never uses the label;<br>– `logit_fit` has no convergence flag, and a separated design would not raise `LinAlgError`;<br>– the md never names the Cochran–Armitage test or its binomial variance;<br>– E′'s "below half the spacing … they cannot differ" holds for E's pair (0 of 1600 centroids). Its reason is the fit's centre: for (0.8, −0.4), ĉ crosses the tie while c̄ < 0.05, and the entrants differ in 11 of 1600 (mine);<br>– "holds its level and grows conservative" needs "level" = size ≤ α. | **Code:** remove the dead code; add the convergence flag.<br>**md:** name the test.<br>**Text:** "for E's pair"; reword the level sentence. |

## Sources

- **Read in the export:** D52.md; `exp_close_features.py`; close_features.md; `test_close_features.py`; REVIEW_round78_D3by.md; the REVIEW_LOG entries for rounds 73, 75 and 78; the diffs 6764e41..5201021; `rich_transit.Bumps` and `sliding_iss`; `exp_mirror_transit`'s constants; `exp_cohesion_leak.load/save`; D3by's mutant file.
- **Web (snippets; arxiv.org unreachable):**
  - [Heinrich–Kahn 2018, Ann. Statist. 46(6A)](https://projecteuclid.org/euclid.aos/1536307235) and [its PDF record](https://projecteuclid.org/journals/annals-of-statistics/volume-46/issue-6A/Strong-identifiability-and-optimal-minimax-rates-for-finite-mixture-estimation/10.1214/17-AOS1641.pdf): the rate n^{−1/(4(m−m₀)+2)} and "corrects a previous paper by Chen". Search summaries add Chen's claimed n^{−1/4} and n^{−1/6} for one extra component [snippet].
  - [Liu–Shao 2003, Ann. Statist. 31(3)](https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-3/Asymptotics-for-likelihood-ratio-tests-under-loss-of-identifiability/10.1214/aos/1056562463.full), with [a 2022 review, arXiv 2206.15178](https://arxiv.org/pdf/2206.15178): the limit is the supremum of a squared left-truncated Gaussian process [snippet].
  - [Davies' publication list](https://robertnz.net/publicat.html): Biometrika 64 (1977) and 74 (1987); the maximum of the score over the nuisance parameter [snippet].
- **[memory]:**
  - overdispersion and inference with the cluster as the unit;
  - the Cochran–Armitage variance;
  - the exact Monte Carlo test rule;
  - the n^{−1/4} order of contiguous alternatives;
  - Chen–Chen–Kalbfleisch 2001.

## Lesson

Draws replicate only within a null. Round 78 pooled 24 600 draws of 39 nulls and tested a trend in the distance as if every draw were an independent trial at its null's distance. But nulls at one distance disagree far beyond binomial noise, and two independent experiments reproduce each other's outliers: the distance is not what sets the size. The variable that orders the nulls, the runner-up's margin over its noise, was already printed in E′ beside the two lowest sizes. Before fitting a trend across groups, check that groups at the same value agree, and make the group the unit of inference. A difference that a test fails to detect is unresolved, not absent: that is round 78's own rule for trends, applied to differences.
