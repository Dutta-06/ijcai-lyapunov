# Referee report D3cg, round 86: Proposition D‴.52 as revised after round 82

## Scope

- **Commit reviewed:** dad96ae on claude/festive-bardeen-eopej2 ("D‴.52 round 82: the size is ordered by the runner-up's margin at the null, not by the distance").
  - I exported it with `git archive dad96ae | tar -x` to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/repo`.
  - All other work is in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work`, in its scripts/, out/, logs/, mut/, cfout/ and pytmp/ folders.
  - HEAD moved to 066dc8f while I worked, through three commits by other sessions (D‴.59 round 83; the round-84 and round-85 reports).
    - `git diff dad96ae 066dc8f` changes none of the reviewed files and not the README's D52 row.
    - REVIEW_LOG gains only a Round 83 entry, at its end.
- **Files read in full:**
  - D52.md (311 lines); `exp_close_features.py` (1805 lines); close_features.md (404 lines); `test_close_features.py` (62 tests); REVIEW_round82_D3cc.md;
  - the Round 82 entry of REVIEW_LOG.md, and the Round 73, 75 and 78 entries with their round-82 errata;
  - `git diff 5201021 dad96ae` for the code, the md, D52.md, the tests and the README's D52 row.
- **Read in part:** `rich_transit.Bumps`; `exp_cohesion_leak.load/save`; the constants of `exp_mirror_transit`. Also, read-only in its own scratch folder, round 82's referee's mutant list and its pre-registration note.
- **Settings:**
  - the repository's venv: Python 3.13.16, numpy 2.5.4, scipy 1.18.1;
  - every Python process run with `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`;
  - pytest as `-q -p no:cacheprovider` with `--basetemp` in my folder;
  - `CF_OUT` always a folder of mine holding a copy of the raw (a read-only copy for the mutant runs);
  - at most two heavy processes of mine at a time.
- **I modified nothing in the repository.** I ran only read-only git commands: archive, log, show, diff, status, rev-parse and check-ignore.
  - The working-tree edits present at the start (D59.md, the README, an untracked round-84 report) belonged to other sessions. Those sessions have since committed them, and `git status` is clean.
  - Files written in the repository during my session are other sessions' work: `exp_width_mismatch`, `exp_noise_resonance`, their tests and bytecode, and the round-85 report. I imported none of them. Every import of mine came from my export, with bytecode writing off.
  - The reviewed files are byte-identical at HEAD. All my processes have exited.
- **Tags:** (md) is close_features.md, (raw) is close_features_raw.npz, (mine) is my own computation, [snippet] is a web search result and [memory] is my own knowledge. I did not use arxiv.org or ams.org.

## Summary

**What checks out**

- **Reproduction:** `--from-raw` regenerates the md byte for byte (md5 8bab97a6336d7cb0b85d1c2a12e353c1). The 62 tests pass in 24.3 s (mine).
- **Every new md number recomputes from the raw with my own code.** I used my own kernels, margins, delta method, Clopper–Pearson, Cochran–Armitage, permutations and quasi-binomial maximum likelihood (scipy). That covers:
  - all 69 rows of F″ and its pooled lines;
  - the "where the calibration loses" line: 79 of 266 at 0.46, and 1258 of 1365 at 1.00;
  - the null-as-unit line. Its permutation p are within their 2000-permutation error: mine with 10⁵ are 0.0458, 0.0163, 0.7301 and 0.2841, against 0.042, 0.017, 0.736 and 0.277;
  - G's table and its pooled line.

  Every decimal number on D52.md's 96 changed lines, and on the README's D52 row, is in the md (mine).
- **The raw holds where I reran it (mine):**
  - I reran G's cells 0 and 2 (the strong band null at s = 0.1 and the weak one at s = 0.05), 3000 draws each, with the repository's `f_job`. All 12 stored fields of every draw match bit for bit, as do q_law, qA, qB and np.
  - The stored bootstrap p recompute exactly from the stored fits, for 4400 draws of 11 F nulls.
  - The delta-method sd(μ̂)/s matches a 20 000-draw Monte Carlo at each null's own s within 0.3–0.7%. No null changes side of the band.
- **The central claim holds: μ₀/sd(μ̂), not the distance, orders the sizes.**
  - At one distance the size is not monotone in s. At d = 0.049 the strong null gives, in F and E:
    - 0.055 and 0.044 at s = 0.003;
    - 0.015 and 0.019 at s = 0.01;
    - 0.045 and 0.042 at s = 0.03 (md).
    
    No monotone function of the distance, s, s/|ŵ|, μ₀ or sd(μ̂) produces this. μ₀/sd(μ̂), at 6.96, 2.10 and 0.74, does.
  - No single covariate comes near the band:
    - deviance 83.7–105.9 against 39.0 for the strong null, and 36.1–45.7 against 26.3 for the weak one;
    - none adds anything beside the band (|z| ≤ 0.52 and ≤ 1.15);
    - with distinct nulls as the units, the band keeps z = −6.47 and −3.81 (mine).
- **Out of sample, the same-draw comparison is decisive (mine).**
  - G's band nulls lose rejections against the law at the truth on their own draws: McNemar z = +6.62, +3.68 and +3.45 (exact p 4·10⁻¹², 3·10⁻⁴ and 7.5·10⁻⁴).
  - The other three do not: +0.83, +0.41 and +1.79.
- **The mechanism holds when tested directly (mine).** Take the bootstrap's α quantile at the law-rejected draws' own fitted nulls, over the law at the truth's:
  - for μ₀/sd(μ̂) = 1.01–3.98 it is 0.24–0.83, while the other draws give 0.98–1.16;
  - above the band it is 0.92–1.00.
- **Sources:** the Chen 1995 / Heinrich–Kahn 2018 sentence is now right, and Davidson–MacKinnon's description fits [snippet].
- **Fixed:** D3cc-1 in substance, and D3cc-3 to D3cc-7 (one overshoot). All 13 of D3cc's named surviving mutants are now killed by the tests (mine).

**What does not hold (no major issue)**

1. **D3cg-1 (minor). G is compared with the whole in-sample band, not like for like.**
   - G's band nulls sit at μ₀/sd(μ̂) = 3.5–3.9. At those values they replicate the in-sample band: size 0.0300 against 0.0253 (z = +1.32), and same-draw loss 0.0126 against 0.0170 (z = −1.22).
   - "The size does not replicate" comes from the in-sample nulls at 1.2–2.5 (0.0187), which G never sampled.
2. **D3cg-2 (minor). The null-as-unit tests count runs, not nulls.**
   - E and F run the same strong nulls, and the mirror offsets are nearly the same nulls. So 39 runs are 15 nulls, and the band's 5 are 2.
   - Permuting the distance over nulls gives p = 0.157, not 0.042.
   - The weak null's p = 0.017 shows a trend; it is not "weak".
   - The band's verdict survives.
3. **D3cg-3 (minor). The plug-in's loss begins below the band's lower edge.**
   - The null at 1.01 loses 11 of its 31 rejections and adds none (exact p = 0.001). That is as much net loss per draw as inside the band.
   - Its size reads 0.0500 only because its draws reject 0.0775 at the truth.
4. **D3cg-4 (minor). The printed mechanism evidence does not separate the band.**
   - Below the band, the law-rejected draws' fitted margins are also small (0.08–0.22 of all draws'), yet nothing is lost there.
   - What separates the band is the bootstrap's quantile at those draws' fitted nulls.
5. **D3cg-5 (minor). The tests miss the round's central lines.** Of my 122 mutants:
   - the tests kill 70;
   - 34 more change only the md, among them the whole mechanism line and the whole null-as-unit line;
   - 6 survive both checks and change a rerun of G: its cells, its draws, its seed, its job list, and whether a default run includes it.
6. **D3cg-6 (minor). The remedy "a double bootstrap … would recalibrate" is untested.** The distortion does not shrink with s at fixed μ₀/sd(μ̂): it is first order and local to the tie, not the higher-order distortion that iteration provably reduces.
7. **D3cg-7 (cosmetic). Wording and consistency:**
   - "at most α" beside an interval that reaches 0.0522;
   - two values of μ₀/sd(μ̂) given for the same nulls;
   - a mixed sd range;
   - a variable name printed in the md;
   - "79 of 266" without the 53 rejections the calibration adds;
   - a field that is never read.

## What I checked and how

**Reproduction**
- I ran `--from-raw` with CF_OUT = `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/cfout/fromraw` (a copy of the raw) and compared the output with `cmp`.
  - Log: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/logs/fromraw.log`.
- Tests: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/logs/pytest_base.txt`.

**F″ recomputed independently**
- Each null's j₂, μ₀ and sd(μ̂) from my own kernels; sizes, Clopper–Pearson intervals, the law at the truth, the kept rejections, the fitted-margin medians, the pooled lines and the mechanism line:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/recompute_fpp.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/recompute_fpp.txt` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/nulls_margin.json`.
- The structure of the raw: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/explore_raw.py`.

**Null-as-unit statistics**
- The Cochran–Armitage z, permutations over runs and over distinct nulls (10⁵ each), and quasi-binomial fits by direct maximum likelihood:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/units_check.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/units_check.txt`.
- The band model with distinct nulls as the units:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/cluster_band.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/cluster_band.txt`.

**What orders the sizes**
- Per null, the same-draw counts (lost and added against the law at the truth). Then, with and without G:
  - seven alternative covariates (linear and quadratic), each alone and beside the band;
  - the band's edges, moved anywhere in the gaps E's ratios allow;
  - the paired loss, pooled by μ₀/sd(μ̂).
  
  Script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/covariates.py`; output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/covariates.txt` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/nulls_paired.json`.
- The net loss per draw by μ₀/sd(μ̂):
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/loss_shape.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/loss_shape.txt`.
- G against the in-sample nulls at matched μ₀/sd(μ̂):
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/g_matched.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/g_matched.txt`.
- Exact McNemar tests: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/g_mcnemar.txt`.

**The mechanism**
- Each draw's bootstrap p recomputed from its stored fit with my own knot code, and the bootstrap's α quantile at the law-rejected draws' fitted nulls:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/mech_check.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/mech_check.txt`.
- The same quantile at every draw of 12 nulls:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/mech_all.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/mech_all.txt`.

**Reruns**
- G's cells 0 and 2, all 3000 draws each, with the repository's `f_job`, compared row by row:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/g_rerun.py`;
  - logs `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/logs/g_rerun_0.log` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/logs/g_rerun_2.log`.
- `null_margin`'s delta method against a Monte Carlo at each null's own s:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/sd_mc.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/sd_mc.txt`.

**Numbers**
- Every number on D52.md's added lines and on the README row, checked against the md:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/scripts/numbers_check.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/numbers_check.txt`;
  - inputs `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/d52_added.txt` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/readme_added.txt`.
- I then checked each number on the round's new lines by hand against its md line.

**Mutation testing**
- 122 single-point mutants of this round's new or changed code, plus the identity, in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/mut/mutants.py`. They cover:
  - `null_margin` (NM), `quasi_logit` (QL), `perm_trend` (PT), `dist_bins` (DB), `bin_index` (BI);
  - `margin_nulls` (MN), `in_band` and the band (IB), `pooled` (PO), `fmt_rate` (FR), `margin_lines` (ML);
  - `g_jobs`, `part_G`, `G_CELLS` and `R0_G` (GJ); `f_job`'s seed path (FJ); `logit_fit`'s convergence check (LF);
  - the |ŵ| range of `size_nulls` and `size_by_distance` (SN); E′'s discordances and pooled bound (EP); `report` and `main` (RM).
- Runner: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/mut/run_mut.py`. Each mutant runs in a slim tree (the three experiment modules, `src`, the test file and a read-only raw): first the test file with `-x`, then `--from-raw` with a byte comparison against the commit's md.
- Results: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/mut/results_w0.jsonl` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/mut/results_w1.jsonl`. Summary: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/mut_summary.txt`.
- Probes of the mutants that survive both checks:
  - script `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/mut/probe.py`;
  - output `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/probes.txt`.
- D3cc's 13 named survivors and its L07, ported unchanged:
  - definitions `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/mut/mutants_d3cc.py`, runner `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/mut/run_mut2.py`;
  - results `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/mut/results_d3cc.jsonl` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref86/D3cg_work/out/mut_summary_d3cc.txt`.
- The identity control passes both checks.

## Verdicts per part

| claim (D52.md lines) | verdict | evidence |
|---|---|---|
| Title, (i), (ii), (iv), (v), Proof (l. 1–41, 101–159) | not re-examined (unchanged) | A–D in the md match the Numerics, e.g. 1.0152 → 0.9995 and 91.34 against 91.23 (md) |
| (iii): the knot's law, the LARS identity, the null quantile across a cell (l. 43–54) | **holds** (unchanged) | sd 0.15–0.17s is the strong null's value (md) |
| (iii): calibrated size at the level outside 1.2 < μ₀/sd(μ̂) < 5 (0.0492, 0.0434) and conservative inside it (0.0197, 0.0246), against 0.0450 and 0.0436 at the truth on the same draws (l. 54–61) | **holds; numbers right.** The lower edge is D3cg-3, "at most α" is D3cg-7 | md = mine on all 69 rows. The band beats every single covariate (deviance 39.0 against ≥ 83.7; 26.3 against ≥ 36.1). With edges anywhere in E's gaps, the strong null's in-band runs stay the same 5 (0.0197); the weak one's in-band size moves within 0.0150–0.0285. F's own runs, independent of the selection, give 0.0183 [0.0115, 0.0276] (mine) |
| (iii): "the plug-in causes it", with 79 of 266 kept and fitted margins 0.46 against 1.00 (l. 61–64) | **holds; the printed evidence is not diagnostic** (D3cg-4) | q_fit/q_law is 0.24–0.83 at the law-rejected draws for μ₀/sd(μ̂) = 1.01–3.98, against 0.98–1.16 for the other draws (mine) |
| (iii): Davidson–MacKinnon; a double bootstrap "would recalibrate" (l. 64–66) | **the attribution fits; the remedy is unsupported** (D3cg-6) | D–M: the critical value function, and the independence of statistic and estimate [snippet]. The loss does not shrink with s at fixed μ₀/sd(μ̂) (mine) |
| (iii): with the nulls as the units, the band and not the distance orders the sizes (l. 66–69) | **holds; the p-values are misspecified** (D3cg-2) | Over distinct nulls: the band's z −6.47 and −3.81; the distance beside it −0.33 and −1.09; the distance alone p 0.157 (not 0.042) and 0.0165 (mine) |
| (iii): G out of sample, "by less than the band's in-sample sizes … not a constant" (l. 69–74) | **the direction holds strongly; the comparison is not like for like** (D3cg-1) | McNemar +6.62, +3.68, +3.45 against +0.83, +0.41, +1.79. At matched μ₀/sd(μ̂): 0.0300 against 0.0253 (z = +1.32). Cells 0 and 2 rerun bit for bit (mine) |
| (iii): one-sided by the margins; power; losses; opposite-sign bias; the null moving with Δ (l. 74–99) | **holds** (unchanged; spot-checked) | 390/430 and 328/364; 1.064 [0.987, 1.289]; 1.045 and 1.032 (4 of 24) (md l. 245–247, 285–286) |
| Status (l. 161–170) and New (l. 199–207) | **holds** | As above |
| Status: Chen 1995 / Heinrich–Kahn 2018 (l. 184–186) | **holds** | Chen: n^{−1/4}. Heinrich–Kahn: local minimax rate n^{−1/(4(m−m₀)+2)}, i.e. n^{−1/6} for one extra component [snippet] |
| Numerics E (l. 240–259) | **holds** (unchanged); μ₀/sd(μ̂) given two values (D3cg-7) | E′ gives 3.69 where F″ gives 3.52 (md) |
| Numerics F, size by distance (l. 286–293) | **holds** | The \|ŵ\| ranges are now printed (md l. 293–296). "Every z here with the draws as the units" |
| Numerics F″ (l. 294–299) | **numbers right**; the p-values as in D3cg-2 | md = mine |
| Numerics G (l. 299–302) | **numbers right** | md = mine; the two reruns |
| Numerics E′, Smith, s/\|ŵ\| (l. 302–311) | **holds** | The literal 0.003–0.006 is marked "not reproduced" (md l. 210) |
| README, D52 row | **numbers right**; "at most α" (D3cg-7) | md |

## Status of the round-82 issues

| D3cc id | status | notes |
|---|---|---|
| D3cc-1 (major) | **fixed in substance** | F″ prints every null with μ₀/sd(μ̂), its size and the same draws against the law at the truth, adds null-as-unit tests, and adds G. All its numbers recompute (mine). The band explains the spread in sample and replicates out of sample in the same-draw comparison. "Do not differ" has become "not resolved". The errata are in the Round 73, 75 and 78 entries, and the README is updated. Remaining: the G comparison (D3cg-1), the units (D3cg-2), the lower edge (D3cg-3) and the mechanism's evidence (D3cg-4). |
| D3cc-2 (minor) | **fixed for D3cc's survivors; new gaps** | 9 new tests. All 13 of D3cc's named survivors are now killed: T01, L03, L06, S06–S10, S15, F02, F06, R01, R02. Its L07 survives as equivalent on these data. Of my mutants of the new code, 34 change only the md and 6 change behaviour unseen (D3cg-5). The authors' list of 154 mutants is not in the repository, so its counts cannot be reproduced. |
| D3cc-3 (minor) | **fixed** | l. 184–186 read as proposed [snippet] |
| D3cc-4 (cosmetic) | **fixed** | The md prints 1.3852–1.3988 and 0.4027–0.4306 (md l. 293–296), and D52 l. 290–292 quotes them |
| D3cc-5 (cosmetic) | **fixed, by marking** | "not reproduced here" (md l. 210) and "not reproduced by the repository's scripts" (D52 l. 304–305). The value is still a literal in the script (l. 887–888). |
| D3cc-6 (cosmetic) | **fixed** | The README marks the round-75 and round-78 readings as withdrawn and states the band |
| D3cc-7 (cosmetic) | **fixed, with one overshoot** | The dead `sets`/`lab` code is gone. `logit_fit` now raises when it does not converge, and this is tested (my LF01, LF03 and LF04 are killed). The md names the Cochran–Armitage test and its binomial variance. "For E's pair" is in. The overshoot: the level sentence now says "at most α" against an interval that reaches 0.0522 (D3cg-7). |

## Issues

| id | severity | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3cg-1 | minor | **G is compared with the whole in-sample band, not like for like, so "the direction replicates, the size does not" misreads it.**<br>Where: D52 l. 72–74; md l. 393–402 (the column "the band predicts"); the commit message; the REVIEW_LOG Round 82 entry (the D3cc-1 row, and the Lesson's "out of sample the band loses half as much as in the cells it was read from").<br>– **What G sampled.** G's band nulls all sit at μ₀/sd(μ̂) = 3.52–3.91, and its other three at 6.32–9.85. G samples neither the band's lower part (1.2–2.5), where the in-sample loss is largest, nor anything below the band.<br>– **Where it sampled, G replicates (mine).** The in-sample runs with 2.5 < μ₀/sd(μ̂) < 5 give 76/3000 = 0.0253, against G's 270/9000 = 0.0300 (z = +1.32). Against the law at the truth, 0.0423 and 0.0426. Net same-draw loss per draw, 0.0170 and 0.0126 (z = −1.22). For the weak pair alone, 0.0317 against 0.0323.<br>– **Where the gap comes from (mine).** The in-sample runs at 1.2–2.5 give 0.0187 (net loss 0.0277; z = +3.31 against G). So the size does vary within the band — the draft's own conclusion — but G shows no failure to replicate.<br>– **G's weak band null is past the tested range.** It runs at s/\|ŵ\| = 0.110, beyond the 0.0021–0.0825 where the first-order law was tested, and its own law at the truth rejects only 0.0393 [0.0327, 0.0469] (md).<br>– **Provenance.** "The first two the round-82 referee's" use Δ = 0.1 and 0.07. The referee's pre-registered runs used Δ = 0.05 (its band_predict_PREREGISTERED.txt). G's nulls therefore sit at d = 0.0313 and 0.0398 with \|ŵ\| ≈ 1.34 and 1.37, not at 0.032 and 0.040 with 1.385. Nothing in the repository records the other four cells before their results; code, raw and md arrive in one commit. | **md (G table):** predict each G null from the in-sample nulls in its own range of μ₀/sd(μ̂) (e.g. 1.2–2.5, 2.5–5, 5–10), not from the band's pooled size. Print the same-draw counts lost and added.<br>**D52 l. 72–74:** "at the margins it sampled (3.5–3.9) G replicates the in-sample band (0.0300 against 0.0253, z = 1.3); in sample the band's lower part loses more (0.0187 at 1.2–2.5): the size varies within the band".<br>**Log:** correct the D3cc-1 row and the Lesson. Say that G's first two cells take the referee's (s, c̄) at another Δ. Record future cells and predictions in a committed file before running them. |
| D3cg-2 | minor | **"With the nulls as the units" counts runs, not nulls.**<br>Where: D52 l. 59, 66–69, 294–299; md l. 389; the round-82 erratum in the Round 73 entry. Round 73's own Lesson says it: "an offset and its mirror are one data point seen twice".<br>– **The units are duplicated (raw, mine).** E and F run the same strong nulls: at s = 0.01 and 0.03, E's five cells are F's, with identical ĉ, ŵ and μ₀/sd(μ̂); at s = 0.003 they agree within 6·10⁻⁴ in ŵ. F's offsets 0.055, 0.068 and 0.08 mirror 0.045, 0.032 and 0.02: distances within 2·10⁻⁴, and μ₀/sd(μ̂) within 3.2%. So the strong null's 39 "nulls" are 15 distinct (s, distance) nulls, and the band's 5 are 2.<br>– **Over distinct nulls the strong trend is not significant (mine).** Permuting the distance over the 15 nulls gives p = 0.157, against 0.046 over the 39 runs (10⁵ permutations each; md 0.042 with 2000). The quasi-binomial slope becomes z = −1.44 with φ = 5.14 on 13 df, not −2.12 with 2.40 on 37. The weak null's 24 runs are 20 nulls: p = 0.0165, z = −2.41.<br>– **"The distance's trend is weak" does not fit the weak null.** Its p = 0.017 and z = −2.54 show a trend. The band explains it: beside the band the distance's z is −1.14.<br>– **The conclusion stands over distinct nulls (mine):** the band alone has z = −6.47 and −3.81, and the distance beside it −0.33 and −1.09. | **md (`margin_lines`):** before the null-as-unit tests, pool the runs of one null (E and F at one (s, ĉ)), or permute within strata. Mark mirror offsets as near-replicates. Print the number of distinct nulls.<br>**D52 l. 59, 66–67, 294–297 and the Round-73 erratum:** "5 runs of 2 nulls". And: "the distance's trend (p = 0.16 over the strong null's 15 nulls, 0.017 over the weak null's 20) is the band's: beside it z = −0.33, −1.09". |
| D3cg-3 | minor | **The plug-in's loss does not stop at the band's lower edge, and the size against α hides it.**<br>Where: D52 l. 55–59 ("at most α … except where … 1.2 < μ₀/sd(μ̂) < 5") and l. 73 ("the band locates where the plug-in loses its level"); md l. 316.<br>– **The null at 1.01 loses as much as the band (mine).** The weak null at μ₀/sd(μ̂) = 1.01 ((0.8, −0.4), s = 0.03, c̄ = 0.045) keeps 20 of its 31 rejections under the law at the truth and adds none (McNemar z = +3.32, exact p = 0.001). Its net loss is 0.0275 per draw. The band's is 0.0296 at 1.7–2.5 and 0.0150 at 1.2–1.7.<br>– **Its size hides this (md, mine).** It reads 0.0500 only because its draws reject 0.0775 [0.053, 0.108] under the law at the truth.<br>– **Below it there is no loss (raw).** The three runs at 0.30 and 0.74 give 0.0000 (0.0018).<br>– **So the lower edge lies between 0.74 and 1.01, and it rests on four runs.** Above the band: 0.0006 at 5–7, then 0.0032 (0.0017) at 7–12, and −0.0003 beyond (mine). | **md:** print the same-draw counts lost and added against the law at the truth (n10, n01), per null and pooled by μ₀/sd(μ̂); my table is in loss_shape.txt.<br>**D52:** say "from about 1 to 5", or note that the effect at 1.01 is masked by that null's draws. |
| D3cg-4 | minor | **The printed evidence of the mechanism does not separate the band.**<br>Where: D52 l. 61–64; md l. 387.<br>– **Small fitted margins are not specific to the band (mine).** Below the band, the law-rejected draws' fitted margins are also a small share of all draws': 0.08–0.22 at μ₀/sd(μ̂) = 0.30–1.01, against a median of 0.46 in the band. Yet below the band nothing is lost except at 1.01: 16 of 16 and 18 of 18 kept, and E keeps 37 of 42 and adds 5. The "1.00 outside" is a median over 51 runs, 47 of them above the band.<br>– **What separates the band is the bootstrap's α quantile at those draws' own fitted nulls, over the law at the truth's (mine):**<br>  – 0.24–0.83 at 1.01–3.98, while the other draws give 0.98–1.16;<br>  – 0.99 at 4.43, where nothing is lost;<br>  – 0.92–1.00 above the band;<br>  – 1.04 and 0.77 below it, where the rejected knots lie deep in the tail (median 0.39–0.41 of q_law).<br>– The draft's mechanism holds, but this is the evidence that shows it, not "0.46 against 1.00". | **md (`margin_lines`):** per null and by μ₀/sd(μ̂), print the median q_fit/q_law of the law-rejected draws and of the others.<br>**D52 l. 63–64:** cite that. Also say that the early draws' fitted margins are just as low below the band, where the null law's quantile is flat in the margin. |
| D3cg-5 | minor | **The tests miss the round's central lines, and G's provenance (mine).** Of my 122 single-point mutants of this round's code, the 62 tests kill 70.<br>– **34 pass every test and change only the md** (caught by the md diff alone):<br>  – the whole null-as-unit line (md l. 389): PT04, PT07, ML11–ML19, ML28, ML29, ML35, ML36. This includes the trend in μ₀/sd(μ̂) instead of the distance, the band's and the distance's columns swapped, wrong standard errors and df, and the binomial p printed as the permutation p;<br>  – the whole mechanism line (l. 387): ML07, ML08, ML09, ML30, ML34;<br>  – the per-null table: MN03 (the distance by the label), MN05 (each draw's fitted margin replaced by the null's own), MN09, MN10, ML01, ML04, ML27;<br>  – PO01 (p ≤ α; 28 H0 draws have p exactly 0.05);<br>  – the band's edges (IB03, IB05);<br>  – G's distance bins (ML22, ML23) and G dropped from the md (RM01, RM03).<br>– **6 survive both checks and change behaviour** (probe.py):<br>  – FJ03: `f_job` seeds every G cell with args[5] = 3000, so rows and q_law differ from the raw;<br>  – GJ03: `part_G` runs F's 336 jobs;<br>  – GJ04, GJ06: a changed `G_CELLS`;<br>  – GJ05: `R0_G` = 2999;<br>  – RM05: a default run skips G.<br>  The tests compare `g_jobs` with `G_CELLS` itself, and nothing ties them to the raw.<br>– **12 survive both and are equivalent.** PO03, MN06 (a field never read), PT06 (an unreachable branch) and LF06 always. IB01, PO02, ML02, MN02, ML10, SN01, LF02 and LF05 on these data: no tie at an edge, no tie at q_law, no H1 row in G, no negative ŵ, and the fits agree to the printed precision. | **`test_close_features.py`:**<br>– tie `G_CELLS`, `R0_G`, the seeds and `part_G`'s job list to the raw's G (each row's (A, c̄, s, Δ, hyp) and its row count);<br>– rerun each G cell's first draw against the raw, about 0.3 s per cell;<br>– assert that `main`'s default parts include G;<br>– in `test_margin_lines_on_synthetic_nulls`: give draws fitted centres off the null's, a null whose centre is off its label, p-values exactly at α, and ratios near 1.2 and 5.5; then assert the "Where the calibration loses" and "With the nulls as the units" lines against `perm_trend`, `quasi_logit` and an independent fit. |
| D3cg-6 | minor | **The remedy is asserted for a regime where it is not assured.**<br>Where: D52 l. 64–66: "which a double bootstrap (Beran's prepivoting) would recalibrate [memory]". It is untested here.<br>– **The distortion does not shrink with s at fixed μ₀/sd(μ̂) (mine).** The net loss per draw is 0.015 at 1.44 with s = 0.003, 0.025–0.028 at 1.80–1.88 with s = 0.03, and 0.018 at 3.80 with s = 0.1.<br>– **So it is a first-order effect.** It comes from the local parameter μ₀/sd(μ̂) near the tie: the parametric bootstrap's non-uniformity near a non-smooth point of the null law, the setting of Andrews 2000 [snippet].<br>– **The cited results concern a different regime.** Davidson–MacKinnon's size distortion [snippet] is the error a bootstrap leaves in a regular model. Iteration removes an order of that error there (Beran 1988; Hall–Martin 1988 [snippet, memory]). Prepivoting can restore validity in some non-standard cases (Cavaliere et al. [snippet]), but nothing here shows it would in this one. | **D52 l. 64–66:** "which a double bootstrap might reduce (untested)", citing Andrews 2000 for the regime. Or run a double bootstrap on the band's nulls. |
| D3cg-7 | cosmetic | **Wording and consistency:**<br>– "the knot's size is at most α (0.0492 [0.0463, 0.0522] …)": the interval reaches 0.0522, so "at the level" fits, as the Round-73 erratum says. "At most α … except where … there it is conservative" also reads as a contradiction. The README has "at most α elsewhere".<br>– "the strong null's 34 nulls below": two of the 34 lie below the band, so "below" is ambiguous.<br>– Two values of μ₀/sd(μ̂) for the same nulls:<br>  – l. 257–258 quote E′'s empirical 6.95, 2.08, 0.73 and 3.69 (sd over 1000 draws), where F″ uses the delta method's 6.97, 2.10, 0.74 and 3.52; my 20 000-draw Monte Carlo gives 3.54;<br>  – likewise "sd(μ̂)/s 0.148–0.170" (l. 257) against "0.1533s–0.1679s" (l. 58).<br>– The weak range "0.4791s–0.5782s" (l. 58) includes G's null; the in-sample range is 0.4995–0.5782.<br>– md l. 391 prints "R0_G draws each" — a variable name.<br>– "79 of the law at the truth's 266 rejections" omits the 53 rejections the calibration adds in the band (90 outside).<br>– `margin_nulls`' w1 is never read. | D52 l. 55–59 and 257–259; the README; md l. 387 and 391; the code |

## Sources

- **Read in the export:**
  - D52.md, `exp_close_features.py`, close_features.md, `test_close_features.py`, REVIEW_round82_D3cc.md;
  - REVIEW_LOG.md: the Round 73, 75, 78 and 82 entries;
  - `git diff 5201021 dad96ae`;
  - `rich_transit.Bumps`, `exp_cohesion_leak.load/save`, `exp_mirror_transit`'s constants.
- **Read-only, in round 82's referee's folder:**
  - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/out/band_predict_PREREGISTERED.txt`;
  - `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref82/D3cc_work/mut/mutants.py`.
- **Web [snippet]:**
  - Davidson–MacKinnon, "The size distortion of bootstrap tests" (Econometric Theory 1999). The accuracy depends on the shape of the critical value function, and asymptotic independence of the statistic and the bootstrap DGP lowers the order of the distortion. Pages: [econpapers](https://econpapers.repec.org/RePEc:qed:wpaper:932), [QED wp 936](https://ageconsearch.umn.edu/record/273347/files/qed_wp_936.pdf), [QED wp 1035](https://ageconsearch.umn.edu/record/273505/files/qed_wp_1035.pdf).
  - Andrews 2000, Econometrica 68(2), 399–405: the bootstrap is inconsistent at the boundary ([econometricsociety.org](https://econometricsociety.org/publications/econometrica/2000/03/01/inconsistency-bootstrap-when-parameter-boundary-parameter-space)).
  - Beran 1988, JASA 83, 687–697: prepivoting, where each iteration lowers the level error by about n^{−1/2} (Hall–Martin 1988), as summarised in [arXiv 2512.10546](https://arxiv.org/pdf/2512.10546), [arXiv 2102.04423](https://arxiv.org/pdf/2102.04423.pdf) and [arXiv 1301.7026](https://arxiv.org/pdf/1301.7026). These are search snippets only; I did not fetch them.
  - Cavaliere et al., "Bootstrap inference in the presence of bias" ([pure.au.dk](https://pure.au.dk/ws/files/443099498/Bootstrap_Inference_in_the_Presence_of_Bias-1.pdf)).
  - Heinrich–Kahn 2018, Ann. Statist. 46(6A), 2844–2870: the rate n^{−1/(4(m−m₀)+2)}, and Chen 1995's n^{−1/4} ([projecteuclid](https://projecteuclid.org/euclid.aos/1536307235)).
- **[memory]:** pseudo-replication and the randomisation unit (Hurlbert 1984); McNemar's test for paired rejections; the non-uniformity of the bootstrap near non-smooth points (Andrews–Guggenberger 2009–10); iterated-bootstrap refinements in regular models (Hall–Martin 1988).

## Lesson

Compare two tests run on the same draws draw by draw. This round set the calibrated knot's size against α and against the law at the truth as two separate rates, and each rate carries its draws' own binomial noise. Counting, draw by draw, the rejections the calibration drops and those it adds removes the noise the two share. Counted that way, three things change:
- G's nulls replicate the in-sample loss at the margins where they were run. The comparison that seemed to fail set G's three nulls at μ₀/sd(μ̂) ≈ 3.5–3.9 against a band whose largest losses lie at 1.2–2.5.
- A null just below the lower edge loses as much as those inside the band. Its size reads at α only because its own draws rejected 7.75% under the truth.
- A test that takes the null as its unit needs distinct nulls. E and F run the same strong nulls, so 39 runs are 15 nulls, and permuting over runs gives p = 0.05 where permuting over nulls gives 0.16.
