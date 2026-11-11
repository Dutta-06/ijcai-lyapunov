# Referee report D3by, round 78: Proposition D‴.52 as revised after round 75

## Scope

- **Commit reviewed:** 6764e41 on branch claude/festive-bardeen-eopej2 ("D‴.52 round 75: size by each null's own distance; the off-tie loss pooled; the tie's bound").
  - I exported it with `git archive 6764e41 | tar -x` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/repo`.
  - All other work is in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work`, in its scripts/, out/, logs/ and mut/ folders.
  - HEAD moved while I worked: first to 7d2e76d, then to 74e57ea, through other sessions' D‴.54, D‴.58, D‴.59 and D‴.60 commits. `git diff 6764e41 74e57ea` is empty for every file reviewed here, and it touches no D‴.52 entry of REVIEW_LOG.md [own].
- **Files read:**
  - D52.md and `exp_close_features.py`, both in full;
  - close_features.md, sections A, B, the C/D summaries, E, E′ and F;
  - in `test_close_features.py`, the five round-75 tests and the older tests my mutants touch;
  - REVIEW_round75_D3bv.md;
  - the Round 71, 73 and 75 entries of REVIEW_LOG.md, with the errata added in place;
  - the diff 3dc9443..6764e41;
  - `rich_transit.Bumps` and the constants of `exp_mirror_transit`.
- **Settings:**
  - the repository's venv (Python 3.13.16, numpy 2.5.4, scipy 1.18.1), run with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`;
  - pytest as `python -m pytest -q -p no:cacheprovider`;
  - `CF_OUT` always pointed to a scratch copy of the raw;
  - at most two heavy processes of mine at a time.
- **I modified nothing in the repository** and ran no git command that writes.
  - During the run, `git status` showed seven modified D‴.59 files. They belonged to another session, which then committed them as 74e57ea. The working tree is now clean.
  - Nothing of mine is still running (checked with `ps`).
- **Tags:** [md] is close_features.md, [raw] is close_features_raw.npz, [own] is my run or derivation. For the literature, [snippet] means seen in a search result and [memory] means my own knowledge.

## Summary

**What checks out**

- **Reproducibility is exact** [own].
  - `--from-raw` regenerates the md byte for byte (`cmp`; md5 2067128b…).
  - The 51 tests pass in 18.2 s.
  - Against round 73's raw, parts A, B, C, E, S and F are identical. U differs only by the new fields j1 and j10.
  - All 15 U jobs rerun bit-identically in every field, with j1 = j10 = 10 in every cell.
- **Every number on the 32 lines round 75 changed in D52.md matches the md, and the md matches the raw** [own].
  - I recomputed from the raw:
    - the per-null size bins (0.0517 … 0.0383 and 0.0508 … 0.0325);
    - the per-cell parts, with McNemar standard errors;
    - the pooled parts per offset, −0.1042 (0.0119), −0.0850 (0.0109), −0.0292 (0.0081), using √(Σse²)/3;
    - the bounds 1 − 0.05^{1/400} = 0.00746 and 1 − 0.05^{1/2400} = 0.00125.
  - Over the whole file, all 598 decimal numbers are found in the md, as verbatim copies (548), roundings (46) or derived values (4: 90.7%, 4.7%, 5.6·10⁻⁵, 6.5·10⁻¹²).
- **The tie's zero holds in the population.** With 10⁶ shared law samples and 10⁵ fresh H1 draws per cell, the full and own α quantiles are equal in 6 of 6 tie cells, and 0 of 6·10⁵ draws fall between the two knots [own].
  - The vectorised fits were first checked against E's 21 000 stored fits: centres agree to 1.5·10⁻¹⁵, knots to 8.6·10⁻¹¹ relative, and the laws' quantiles exactly.
- **The off-tie parts are consistent with the population** [own].
  - Population parts per cell: 0.083, 0.071, 0.082 at c̄ = 0; 0.080, 0.071, 0.082 at 0.02; 0.034, 0.033, 0.039 at 0.032 (SE ≈ 0.001). They agree with D3bv's independent run to within 0.005.
  - The nine sample parts give χ² = 7.6 on 9 df.
  - The laws' Monte Carlo error, which the md does not include, is 0.003–0.006 per cell and raises the pooled SEs by at most 4%.
- **The pooled 2400-draw bound is valid** for the mean of the six cells' shares, by AM–GM: Π(1 − p_i)^400 ≤ (1 − p̄)^2400 [own].
- **The round-75 code changes are correct:**
  - the excess now counts an infinite path knot beside a finite own knot as an exceedance;
  - `null_distance`, `size_by_distance` and `part_se` compute what their docstrings say;
  - the strong null's mirror offsets agree.
- **The claimed laws hold:**
  - **Δ₅₀ ∝ N^{−1/4}** is confirmed numerically, though no run in the repository varies N [own]. I ran N = 61, 121, 241 and 481 samples, with the same 21 kernels, s = 0.01 and c̄ = 0.032. The slopes of log Δ₅₀ in log N are −0.247 for the law-calibrated knot, −0.249 for the score and −0.250 for the Neyman–Pearson bound. The ratios to N = 61 are 0.848, 0.713, 0.601 for the knot, against a prediction of 0.843, 0.709, 0.597.
  - **Smith's constants** reproduce with an independent finite-difference CRB: 500.5/504.9Δ³, 52.0Δ², 56.2Δ³, 119.2Δ², 49.6Δ², and the fitted 473.7Δ^{2.99} and 108.9Δ^{1.98} [own].
  - **(ii)'s continuum constants, (iv)'s cosine, the fit shift −T₃/(4Aσ²) and the LARS identity** all re-derive by hand [own].
  - **Who sets the null law:** far kernels set its median (84–92% of the median band). Toward the tie its 5% tail passes from far kernels (60–63% at c̄ = 0) to the neighbours (100% at c̄ ≥ 0.045) [own].
- **D3bv-2, -5, -6 and -7 are fixed. D3bv-1, -3 and -4 are partly fixed** (see the status table).

**What is wrong (majors first)**

There is no major issue. There are three minor issues and four cosmetic ones.

1. **D3by-1 (minor). "The size's fall with the distance is the weak null's" (l. 58–59) is contradicted by the repository's own data on the strong null.**
   - E has 15 000 H0 draws of the same null (0.8, 0.6) under the same calibration. Its sizes fall from 0.050 to 0.035 [0.029, 0.042] next to the tie (trend p = 0.003).
   - In F, the two nulls' slopes do not differ (p = 0.15).
2. **D3by-2 (minor). The new tests miss most of what they are credited for.**
   - The six FA mutants named in the new test's docstring all survive.
   - Undoing D3bv-1's fix (binning by the label c̄) survives the tests.
   - The pooled part and pooled bound that (iii) quotes are untested.
   - Of the 28 earlier mutants D3bv-4 named, the tests now kill 6.
3. **D3by-3 (minor). Three of the refitted sources still do not carry what they are cited for.**
   - Chen 1995 is an estimation rate, and Heinrich–Kahn 2018 (already cited) correct it.
   - Liu–Shao 2003 and Gassiat 2002 are density/mixture results, not "the regression analogue".
4. **D3by-4 to D3by-7 (cosmetic):**
   - the distance bins are intervals labelled by points;
   - the pooled part at c̄ = 0 sits 2.1 SE above the population;
   - the "at 95%" bounds are one-sided and per cell;
   - loose ends.

## What I checked and how

All paths below are under `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work`, and are written out in full.

**Reproducibility**
- `--from-raw` with `CF_OUT=/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/fromraw`, then `cmp` against the commit's md.
- Test output: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/pytest_base.txt`.
- The raw compared with round 73's, part by part: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/cmp_raw.py` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/cmp_u.py`.
- All 15 U jobs rerun: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/rerun_u.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/rerun_u.txt`.

**Numbers**
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/numbers_check.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/numbers_check.txt`.
- The script's matching is loose. I therefore checked every number on the round-75 lines by hand against its md source, and spot-checked the md tables for A, B, D and Smith.

**Size by distance**
- Per H0 job, with the null's own centre: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/size_jobs.py`.
- Cochran–Armitage trends, and logistic models with pair × distance and log(s/\|ŵ\|) terms: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/size_trend.py`.
- E's calibrated sizes: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/e_sizes.py`.
- Explicit distance intervals, with E and F pooled for the strong null: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/size_bins.py`.
- Outputs: the matching `.txt` files in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out`.

**E′**
- `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/eprime_check.py` recomputes the parts, the SEs (ddof 0 and 1, and the McNemar form), the gain and loss counts, the pooled values and the bounds.
- Output: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/eprime_check.txt`.

**Population**
- The vectorised `one_fit` and knots are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/vec.py`.
- They were validated against E's and U's stored values by `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/validate_vec.py`, with output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/validate_vec.txt`.
- The run itself is `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/pop_part.py`. For each of the 15 cells it uses:
  - 10⁶ law samples, shared by the full and the own knot, with full ≤ own asserted on every sample;
  - 10⁵ fresh H1 draws;
  - 50 disjoint 20 000-sample subsets of the law, to measure the laws' Monte Carlo error.
- Summary: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/pop_summary.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/pop_summary.txt`.

**Laws**
- N-scaling: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/n_scaling.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/n_scaling.txt`.
  - It uses the knot at the least-squares one-bump fit, which is the path's knot in the sliding limit.
  - Its vectorised score agrees with `score_stats` to 6.7·10⁻¹⁶.
- Smith's constants: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/smith_check.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/smith_check.txt`.
- Who sets the median and the tail: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/who_median.py`, output in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/who_median.txt`.
- By hand: (ii)'s Rayleigh quotients and continuum integrals; (iv)'s odd/even ratio ϑ/√6 and ∂_σφ = σ∂²_cφ + φ/σ; the shift −T₃/(4Aσ²), which is 104Δ³ for (0.8, −0.4); the ISS/lasso dual identity u_j(τ) = τc_j(1/τ).

**Mutation**
- Definitions: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/mut/mutants.py`, 63 single-point mutants plus the identity:
  - **35 of round 75's code:** X01–X06 (the excess), D01–D03 (nd), P01–P02 (`part_se`), Q01–Q05 (the pooled part), B01–B03 (the pooled bound), J01–J03 (the first entrants), N01–N02 (`null_distance`), S01–S11 (`size_by_distance`);
  - **the 28 earlier mutants named in D3bv-4,** ported unchanged.
- Runner: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/mut/run_mut.py`. Each mutant runs against the full test file and against `--from-raw` with a byte diff of the md. J02 and J03 also rerun the 15 U jobs (`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/mut/behav_u.py`).
- Results: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/mut/results_a.jsonl` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/mut/results_b.jsonl`, with a summary in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/out/mut_summary.txt`.
- The identity passes all checks, and so does its U rerun.

**Literature**
- Web searches only. projecteuclid.org, arxiv.org and academicnewsletter.sufe.edu.cn could not be fetched, so every literature statement below rests on search snippets or memory.

## Verdicts per part

| claim (D52.md lines) | verdict | evidence |
|---|---|---|
| Title: score at the bound with the sign known, ≈ 1.09 (1.04–1.18) without it (l. 1–4) | **holds** | \|z\| 1.085 [1.039, 1.135] and 1.094 [1.059, 1.179]; one-sided z 1.003 and 1.011 [md] |
| (i): the second entry's delay ∝ Δ⁻², the opposite-sign wall, the dipole (l. 8–21) | **holds** (unchanged) | A's ratios 1.015 → 0.9995 and 1.043 → 1.0002; with T₃ 1.229 → 1.003 [md]; ϑ = 25Δ [own] |
| (ii): λ_sep, λ_split, the mobilities, the two regimes (l. 23–41) | **holds** (unchanged) | Continuum integrals √π/(2σ³) and 3√π/(4σ⁵) [own]; λ_split/Δ⁶ = 91.2 against 91.23 [md]; λ_sep 58.1 for (0.8, 0.6) [own = md] |
| (iii): the knot's law, LARS = lasso for two knots, null quantile within a cell (l. 43–54) | **holds** | Dual identity re-derived [own]; median band 84–92% far kernels; 5% tail far 60–63% at c̄ = 0, 25% at 0.032, 0% at ≥ 0.045 [own] |
| (iii): the calibrated size by the null's own distance (l. 54–59) | **numbers right; attribution wrong; labels loose** | Numbers md = own. "The fall is the weak null's" is contradicted by E's strong null and by the interaction test (D3by-1). Bins are intervals (D3by-4) |
| (iii): one-sided by the margins, 75–100% (l. 59–60) | **holds** | 390/430, 328/364; per offset 75.4–100% [md] |
| (iii): power ½ at Δ² ∝ (\|A\|/\|A₁A₂\|)s/√N; 1.06/1.00; 1.045, 3.2% ≈ 13% in samples; 1.11/1.09; 1.06/1.02; 2.09/1.26 (l. 60–68) | **holds** | All match [md]. N^{−1/4} slopes −0.247, −0.249, −0.250 [own]; 1.032⁴ = 1.134 |
| (iii): the minimum's part off the tie, pooled 0.104 ± 0.012, 0.085 ± 0.011, 0.029 ± 0.008 (l. 68–72) | **numbers right; the c̄ = 0 value is high** | md = own. Population 0.079, 0.078, 0.035 (± 0.001); c̄ = 0 at z = −2.1; χ² = 7.6 on 9 df [own] (D3by-5) |
| (iii): none at the tie; 0 discordant; ≤ 0.0075 per cell, 0.0012 pooled; the straddle logic (l. 72–77) | **holds; bounds one-sided** | Population: equal quantiles, 0 of 6·10⁵ between the knots [own]. Bounds 0.00746 and 0.00125 one-sided [own] (D3by-6) |
| (iii): opposite-sign bias; the null moving across the tie (l. 77–83) | **holds** (unchanged) | Least power 0.000–0.025 in 24/24; walls 0.0034/0.0345 against 0.0223/0.0235 [md]. 104Δ³ gives 0.0016 at Δ = 0.025 against the md's +0.0015 [own] |
| (iv): the score's law, cosine, width and split readings (l. 85–97) | **holds** (unchanged) | cos = (1 + ϑ²/6)^{−1/2}; the ∂_σφ identity; the split's sign rule [own] |
| (v): s ≈ 24Δ³ and s ≈ 12Δ² (l. 99–104) | **holds** (unchanged) | 0.84–1.05 up to s/Δ³ = 24 [md]; s/Δ₅₀² = 12.8, 13.3, 11.5 from C [own from md] |
| Proof (l. 106–143) | **holds** | The parts listed above, re-derived by hand [own] |
| Status: what is proved and what is numerical (l. 145–152) | **holds** | Derived and measured claims are separated. N^{−1/4} is "derived" and is now also checked [own] |
| Status: literature (l. 152–178) | **partly holds** | AGM, DCG, Amar–Weiss and El Korso now fit. Chen 1995, Liu–Shao, Gassiat 2002 and "squared supremum" do not (D3by-3) |
| New (l. 179–186) | **holds** | Unchanged in substance |
| Numerics A–D, Smith (l. 188–219, 272–274) | **holds; one label missing** | Smith reproduces [own]. "474Δ^{2.99} and 109Δ^{1.98}" mixes unknown and known amplitudes (D3by-7) |
| Numerics E (l. 219–238) | **holds** (unchanged) | It records the strong null's 0.019 and 0.022 near the tie, which (iii) now overlooks (D3by-1) |
| Numerics F: size by the null's own distance (l. 263–267) | **numbers right; same over-read** | md = own (D3by-1, D3by-4) |
| Numerics E′: part, SEs, bounds, excess (l. 267–272) | **holds** | md = own; 3293 of 21000 with 21000 finite [md = own] |

## Status of the round-75 issues

| D3bv id | status | what remains |
|---|---|---|
| D3bv-1 (minor) | **partly fixed** | The code now bins each null by its own centre, per null [md = own]. The new wording over-reads the other way: E's strong null falls too (D3by-1). The bins are labelled by points (D3by-4). The mutant that undoes the fix (S01) passes the tests (D3by-2) |
| D3bv-2 (minor) | **fixed** | Pooled per offset with SEs [md = own]. Residual: the population values 0.079, 0.078, 0.035 would be the better statement (D3by-5). Q01–Q05 are untested (D3by-2) |
| D3bv-3 (minor) | **partly fixed** | Done: AGM 2006 for the search cost (local power near the level [snippet]); DCG as the mixture analogue; Amar–Weiss [memory]; El Korso beside Liu–Nehorai with HAL records (still no year). New misfits: Chen 1995; Liu–Shao and Gassiat as "the regression analogue"; the positive part (D3by-3) |
| D3bv-4 (minor) | **partly fixed** | Of the 28 mutants it named, 6 are now killed (A06, SL02, SL03, EM02, FR02, FR05). 19 are still caught only by the md diff, including the six FA mutants the new test names. FR01, FR03 and ES03 survive both checks (D3by-2) |
| D3bv-5 (cosmetic) | **fixed** | "(0 discordant of 400: ≤ 0.0075)" and "0.0012" are right [own]. They are one-sided and per cell (D3by-6) |
| D3bv-6 (cosmetic) | **fixed** | X01 and X02 are killed by the new test [own]. The "(N with a finite own knot)" denominator is never exercised (X04, X05: D3by-2) |
| D3bv-7 (cosmetic) | **fixed** | Done: the errata in place at rounds 71 (correct) and 73 (it repeats D3by-1's attribution); "over all cells"; the straddle wording; E′'s header; fewer nested dashes. The first-entrant line cannot fail on E's cells (D3by-7) |

## Issues

| id | severity | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3by-1 | minor | **The new attribution "the size's fall with the distance is the weak null's" is contradicted by the repository's own strong-null data** (l. 55–59; Numerics l. 265–267; README; the Round-75 Lesson and the erratum at Round 73).<br>**E** has 15 000 H0 draws of the same null (0.8, 0.6), with the same calibration (2000 law samples at each draw's fit) [raw].<br>– Its sizes are 0.0500, 0.0527, 0.0490, 0.0420, 0.0350 [0.0287, 0.0422] at c̄ = 0 … 0.049; trend z = −2.97, p = 0.003.<br>– At c̄ ≥ 0.045 the size is 0.0385 [0.0338, 0.0437], against 0.0513 at c̄ ≤ 0.02; Fisher p = 0.0008 [own].<br>– The Numerics themselves quote E's 0.019 and 0.022 next to the tie (l. 228).<br>**F's strong null** at distance ≥ 0.026 is 0.0425 [0.0375, 0.0479] [own].<br>– Its trend p = 0.15 is a failure to detect a trend, not evidence that there is none.<br>– The two nulls' slopes do not differ: in a logistic model the pair × distance term has z = −1.45 (LR p = 0.15), while the common slope has z = −3.22 [own].<br>**Pooling E and F**, the strong null gives 0.0505, 0.0502, 0.0459, 0.0433, 0.0360 [0.0305, 0.0420] (trend z = −3.44). The weak null gives 0.0508, 0.0429, 0.0379, 0.0296, 0.0325 [own].<br>– So both fall. The weak null falls from 0.032 on; next to the tie the two agree (0.036 and 0.033). | **(iii) l. 55–59:** "the knot holds its level across the offset and grows conservative toward the tie for both nulls: the strong null (E and F, 24 600 draws) 0.050 within 0.02 of a kernel, 0.046, 0.043 and 0.036 [0.031, 0.042] at 0.032, 0.045, 0.049; the weak null already 0.038 at 0.032 and 0.030–0.033 beyond".<br>**Numerics:** quote E's per-offset sizes beside F's.<br>**Also:** the same correction in the README line, the Round-75 Lesson and the Round-73 erratum. |
| D3by-2 | minor | **The new tests miss most of what they are credited for** (tests l. 1008–1133; Round-75 log, D3bv-4).<br>**Earlier mutants.** Of the 28 that D3bv-4 named, the tests now kill 6: A06, SL02, SL03, EM02, FR02, FR05.<br>– Still caught only by the md diff (19):<br>  – FA01, FA03–FA07: the docstring of `test_f_alloc_lines_on_rows_that_differ_by_pair_and_offset` names exactly these six, but the test never asserts the per-offset size line or the walls-components line. Its qA ≠ qB never reaches an assertion.<br>  – SL05: the new Smith test uses `crb_sep` as its own oracle, and the older test checks only slopes, so the known-amplitude constants 52Δ², 119Δ², 60.7Δ² are unpinned.<br>  – EM03, EM07, EM09 and FR04.<br>  – ES04, ES09, ES10, ES14–ES18: the `f_lines` test only counts the equal-size line.<br>– FR01, FR03 and ES03 survive both checks.<br>**Round 75's own code.** Of its 35 mutants, 17 are killed and 10 are caught only by the md diff:<br>– **S01** bins by the label c̄, which undoes D3bv-1. The test's shifted null (c̄ 0.045, c₁ 0.0547: distances 0.045 and 0.0453) lands in the same bin either way.<br>– S06 and S11.<br>– Q01–Q05: the pooled part and its SE, which (iii) quotes.<br>– B01 and B03: the pooled 0.0012.<br>The remaining 8 survive both checks and are equivalent on this data: X04 and X05 (every own knot is finite), D01, D02, B02, J02 and J03 (j₁ = j₁₀ in every cell), and S07 (ŵ > 0 everywhere) [own]. | **`test_close_features.py`:**<br>– assert the per-offset size line on rows that differ by pair and offset;<br>– assert the components line with qA ≠ qB and H1 rows at the two smallest separations;<br>– check `crb_sep` against a finite-difference CRB (`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref78/D3by_work/scripts/smith_check.py` reproduces the md's constants);<br>– assert `equal_size_line`'s numbers;<br>– add a `size_by_distance` row whose c₁ falls in another bin than its c̄ (e.g. c̄ = 0.032, c₁ = 0.0417);<br>– assert E′'s pooled line and pooled bound on the raw;<br>– add a row with an infinite own knot.<br>**Docstrings:** list only the mutants each test actually kills. |
| D3by-3 | minor | **Three of the refitted sources still do not carry what they are cited for** (Status l. 162–167).<br>– **Chen 1995** (Ann. Statist. 23, 221–233) gives n^{−1/4} as the best rate for *estimating the mixing distribution* when the number of components is unknown [snippet]. That is not a detection order.<br>– **Heinrich–Kahn 2018**, already cited at l. 178 for (v)'s orders, write "This corrects a previous paper by Chen [Ann. Statist. 23 (1995) 221–233]" [snippet; seen in the abstract as quoted by search results; one result said the published version may word the correction differently, unverified]. Their local minimax rate is n^{−1/(4(m−m₀)+2)}, which is n^{−1/6} for one extra component [snippet]. That is the analogue of (v)'s s ∝ Δ³, not of (iii)'s order.<br>– The detection order n^{−1/4} is that of the contiguous alternatives for testing homogeneity: weights × separation² ∝ n^{−1/2} [memory].<br>– **Liu–Shao 2003** is a quadratic approximation in a Hellinger neighbourhood of the true density, applied to mixtures [snippet]. **Gassiat 2002** is "Likelihood ratio inequalities with applications to various mixtures" [snippet]. Neither is "the regression analogue in additive noise". Davies 1977, 1987 is [memory].<br>– **"Squared supremum":** with positive weights the limit is the squared *positive part* of the supremum [memory]. D3bv-3 already said so. | **Status l. 162–167:**<br>– cite the n^{−1/4} detection order as the order of the contiguous alternatives for testing homogeneity, and drop Chen 1995 from it, or say that Chen 1995 is the (corrected) estimation rate;<br>– keep only Davies under "the regression analogue", and move Liu–Shao and Gassiat to the mixture analogue;<br>– write "the squared positive part of whose supremum". |
| D3by-4 | cosmetic | **F's distance bins are intervals labelled by points** (l. 56–58, 265–267; md F; `size_by_distance`).<br>– Each null goes to the nearest of {0, 0.02, 0.032, 0.045, 0.049}. In effect the bins are [0, 0.01), [0.01, 0.026), [0.026, 0.0385), [0.0385, 0.047) and [0.047, 0.05) [own].<br>– The weak null's "0" bin ("at a kernel's centre") holds its s = 0.03 null at 0.0097, a third of the bin. That null's mirror, at 0.0103, goes to "0.02".<br>– Its "0.045" bin spans 0.0413–0.0465 [own from raw]. | **md F (`size_by_distance`):** print each bin as a range, with its nulls' distances.<br>**(iii) and the Numerics:** "within 0.01 of a kernel" and "at 0.041–0.050". |
| D3by-5 | cosmetic | **The pooled off-tie part at c̄ = 0 sits 2.1 SE above the population** (l. 72, 270).<br>– The md's pooled values are 0.104, 0.085 and 0.029. The population values are 0.079, 0.078 and 0.035 (± 0.001) [own]; D3bv's independent run gave 0.079, 0.078 and 0.036.<br>– So the contrast between c̄ = 0 and 0.02 (0.104 against 0.085) is not real.<br>– The sample is still consistent overall (χ² = 7.6 on 9 df).<br>– The omitted laws' Monte Carlo error is 0.003–0.006 per cell; it raises the pooled SEs from 0.0119 to 0.0123 [own]. | **(iii) l. 72:** "≈ 0.08 at c̄ ≤ 0.02 and ≈ 0.035 at 0.032 (population, 10⁵ H1 draws per cell; the sample's pooled values agree within 2.1 SE)".<br>**Numerics l. 269:** replace "not included" with the size of the laws' error. |
| D3by-6 | cosmetic | **The tie bounds are one-sided and per cell, but are labelled only "at 95%"** (l. 73, 270–271; md E′).<br>– They are the exact one-sided bounds 1 − 0.05^{1/n}: 0.00746 and 0.00125 [own].<br>– Every other interval in the md is a two-sided Clopper–Pearson 95% interval. The two-sided upper ends would be 0.0092 and 0.0015 [own].<br>– The six per-cell bounds hold jointly only at about 0.95⁶ = 74%. | Write "one-sided 95%", or print the two-sided figures. |
| D3by-7 | cosmetic | **Loose ends.**<br>– **l. 272–274:** "474Δ^{2.99} and 109Δ^{1.98}" pair (0.8, 0.6) with *unknown* amplitudes with (0.8, −0.4) with *known* amplitudes [md, own]. Beside the Status's 501Δ³ and 56Δ³ (both unknown), this reads as the same case.<br>– **l. 264:** "at α = 0.05 0.0329–0.0512" still pools both nulls by the label c̄, which D3bv-1 retired.<br>– **E′ summary:** "f₁'s and f₀'s first entrants coincide in 15 of 15 cells" cannot fail at E's offsets, all of which lie in kernel 10's half-cell; J02 and J03 are equivalent on this data.<br>– **l. 162:** "its mixture analogue" follows Neyman–Scott's test of homogeneity, which is itself a mixture test.<br>– **l. 59:** the sentence ends in a dangling "—.".<br>– **El Korso et al.** still have no year. | Label the two fitted constants by pair and amplitude case. Mark the per-offset range as pooled over both nulls. Say the first-entrant check is trivial at c̄ < 0.05. Edit the remaining wording. |

## Sources

- **Read in the export:**
  - in full: D52.md and `exp_close_features.py`;
  - close_features.md, sections A, B, the C/D summaries, E, E′ and F;
  - in `test_close_features.py`, the round-75 tests and the older tests my mutants touch;
  - REVIEW_round75_D3bv.md;
  - the Round 71, 73 and 75 entries of REVIEW_LOG.md;
  - the diff 3dc9443..6764e41;
  - `rich_transit.Bumps`.
- **Web (snippets only):**
  - [Chen 1995, Ann. Statist. 23(1) 221–233](https://projecteuclid.org/euclid.aos/1176324464): n^{−1/4} is the best rate for estimating the mixing distribution when the order is unknown [snippet].
  - [Heinrich–Kahn 2018, Ann. Statist. 46(6A) 2844–2870](https://projecteuclid.org/euclid.aos/1536307235) and [arXiv 1504.03506](https://arxiv.org/pdf/1504.03506): the rate n^{−1/(4(m−m₀)+2)} and "This corrects a previous paper by Chen" [snippet].
  - [Azaïs–Gassiat–Mercadier 2006, Bernoulli 12(5) 775–799](https://projecteuclid.org/euclid.bj/1161614946): local power can be close to the level, and equals it for an unbounded normal mean [snippet].
  - [Dacunha-Castelle–Gassiat 1999, Ann. Statist. 27(4) 1178–1209](https://projecteuclid.org/journals/annals-of-statistics/volume-27/issue-4/Testing-the-order-of-a-model-using-locally-conic-parametrization/10.1214/aos/1017938921.short) and [the 1997 ESAIM P&S paper](https://numdam.org/item/PS_1997__1__285_0) [snippet].
  - [Liu–Shao 2003, Ann. Statist. 31(3) 807–832](https://projecteuclid.org/journals/annals-of-statistics/volume-31/issue-3/Asymptotics-for-likelihood-ratio-tests-under-loss-of-identifiability/10.1214/aos/1056562463.full): Hellinger neighbourhood of the true density; mixture application [snippet].
  - [Gassiat 2002, Ann. IHP 38(6)](https://www.numdam.org/item/AIHPB_2002__38_6_897_0): title [snippet].
  - [hal-00771379](https://hal-supelec.archives-ouvertes.fr/hal-00771379): resolvability of two point sources with a GLRT [snippet]. [hal-00771402](https://hal-centralesupelec.archives-ouvertes.fr/hal-00771402): first-order Taylor model with a GLRT [snippet].
  - [Milanfar, OSA-SRS-05](https://users.soe.ucsc.edu/%7emilanfar/publications/conf/OSA-SRS-05.pdf): SNR ≈ c/(N²d⁴) [snippet]. [Shahram–Milanfar, IEEE TIP 2004](https://users.soe.ucsc.edu/%7emilanfar/publications/journal/DiffLimitFinal.pdf).
  - [IEEE TSP 56(11) issue listing](https://www.datalearner.com/academic/journal-papers/1053-587X/volumes-and-issues/4652320): Amar–Weiss's title only; the δ⁻⁴ law is not confirmed.
  - [Chen–Chen–Kalbfleisch 2001, JRSS-B 63](https://ideas.repec.org/a/bla/jorssb/v63y2001i1p19-29.html): bibliographic record only; its local-alternative rate is not confirmed.
- **[memory]:**
  - the positive part in the limit of the mixture LRT;
  - Davies 1977, 1987;
  - the n^{−1/4} order of contiguous alternatives in testing homogeneity.

## Lesson

A non-significant trend is not a flat one, and a stratum is not exhausted by the cells that first suggested it. Round 75 split F's sizes by null and read the strong null's p = 0.15 as "holds its level", making the fall "the weak null's". But the two slopes do not differ (p = 0.15). E's 15 000 H0 draws of the same strong null, under the same calibration, show it falling to 0.035 next to the tie (p = 0.003). Before assigning a pattern to one stratum, test the interaction and use every sample of that stratum the repository already holds. Likewise, a test credited with killing a mutant must assert the line that mutant changes. Here the f_alloc test names six mutants and kills none of them, and the binning test cannot tell the label from the null's own centre.
