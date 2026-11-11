# Referee report D3dy, round 130: Proposition D‴.59 as revised after round 125

## Scope

- **Export and work directory.** I exported commit 6465011 with `git -C … archive 6465011 | tar -x -C W/repo`, where W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref130. W/repo/venv is a symlink to the repository's venv (Python 3.13.18, numpy 2.5.4, scipy 1.18.1; `np.longdouble` is x87 extended, eps 1.08·10⁻¹⁹).
- **Read in full.**
  - D59.md (446 lines).
  - Every hunk of `git diff 783d55b 6465011` in exp_moving_width.py, and what those hunks call: `a2_exit`/`a2_exit_ld` (`perm`, `line`), `_LDBumps`, `BumpsCW`, `critical_rhs_nd`, `critical_jac_nd`, `_route_rates`, `a2_line_lam`, `a2p_job`, `part_A2P`, `a2p_stats`, `a2p_line`, `a2_tangent`, `a2_state`, `a2_seed`, `a2_sym_lag`, `a2_exits`, the exit and tangent blocks of `a2_lines`, and the mirror block of `split_summary`.
  - The md and its diff; the tests (`test_round125_…`, both round-106 tests, the round-112 test) and their diff; tools/mutants/round125_D59_own.py and tools/mutants.py.
  - review/REVIEW_round125_D3dt.md.
  - The Round 125 entry of REVIEW_LOG.md. It is not in 6465011; I read it from 967a0c9 with `git show`.
- **Settings.** Every heavy run used `OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 nice -n 10`, the repository's venv, and MW_OUT=W/out. pytest ran with `--basetemp` in W. At most two of my heavy processes ran at once.
  - My copy of the mutant tool (W/own/mutants_tool.py) sets REPO to W/repo and gives each slot its own `--basetemp`.
- **Process notes.** I killed nothing. Every process I started ran to completion. Other sessions' processes (PID 1782 `mut125p2`, PID 2099 `ref129`) were left alone. One Monitor expired at 30 min.
- **The repository is untouched.** I used only read-only git.
  - HEAD is 967a0c9.
  - `git diff 6465011` on the D59 files is empty.
  - The working tree holds another session's uncommitted D60 files.
- **Tags.**
  - [md]: the committed md or draft.
  - [raw]: the committed raw.
  - [own]: my computation.
  - [snippet]: a web-search result I saw (no page fetched).
  - [memory]: not verified this round.

## Summary

**Design.**
1. **Reproduction.**
   - md byte for byte from `--from-raw`.
   - The 56 tests.
   - A2P recomputed with `part_A2P` and compared row by row with the raw.
2. **A2P's soundness.**
   - Whether a permutation is arithmetic-equivalent (read from the code).
   - Each run's own effective seed S(s) = c(s)e^{−Λ(0,s)} along its trajectory, against the exit-based ln|S| (is Λ along the line the right integrand?).
   - Where along the route the seed forms.
3. **Fresh draws.** Two new permutation seeds (130 and 2026), 100 orders each, in all four cells: 800 paired draws, with the identity order dropped.
   - I compared round 125's own perm.py, which used `default_rng(125)` with the identity order first.
   - I tested the 3-sem rule on 80 blocks of ten.
   - I checked the +2.1 sem cell and the correlations across tolerances and arithmetics.
4. **The integrator's share of the seed.** 40 orders at 5A (DOP853's default steps, max_step h₀/4 and h₀/16) and at 3A (default and h₀/4), in both arithmetics. Also the step counts at both tolerances.
5. **Normal hyperbolicity.** My own long-double closed-form 3×3 Jacobian at A2E's 46 on-line states, on A's samples and on exactly symmetric ones. Eigenvalues, eigenvectors and the Hirsch–Pugh–Shub rate conditions.
6. **Mutants.**
   - The repository's F00–F29.
   - 18 of my own (G01–G18) and a control.
   - Each first against the three relevant tests. Survivors were then rerun against the whole suite with the md byte test deselected.
   - The tool's `--from-raw` comparison tells whether a published number changes.
7. **Web searches** (5).

**Verdicts.**
- **Reproducibility is exact.**
  - md identical (sha256 08d3c419…4e9263; 30.9 s).
  - 56 passed in 274.9 s (one LSODA warning).
  - `part_A2P` (31.0 s): all 80 exits bitwise identical, ln|S| and Λ differences 0.0, the line identical.
- **The A2P design is sound.**
  - The permutations are arithmetic-equivalent. ⟨·⟩ is a plain mean, f* is evaluated elementwise after the permutation, and on the symmetric samples the ±x terms are bitwise negatives. So the exact sum is 0, and only the association of the pairwise sums changes.
  - Λ along the line is the right integrand: each run's own seed agrees with the exit-based one to 0.997–0.9998.
- **"The seed scales with the arithmetic's roundoff" holds, and more sharply than claimed.** Over 800 fresh draws Λ = 7.652 ± 0.070, i.e. (1.004 ± 0.009)·11 ln 2. The 3A rtol 10⁻¹² cell's +2.1 sem is chance: with 200 fresh draws it is −0.35 sem, prefactors 0.75 against 0.80 uσ*.
- **But the stored draws are round 125's own draws.** Seed 125 is the referee's seed, so the rtol 10⁻¹² cells replicate round 125's table exactly. The 3-sem rule was proposed on those draws, so it is a post hoc check on them.
- **The seed's prefactor is the arithmetic's rounding as sampled by the integrator.** It moves by O(1) in ln|S| with DOP853's steps, while the paired runs share their steps.
- **`a2_tangent` is right.** It agrees with my Jacobian to 1.9·10⁻¹⁵. "Normally hyperbolic only in the repelling phase" is right by the HPS rate conditions.
- **Round 125's issues** are resolved or applied.
- **Mutants:**
  - The repository's 29 are all killed: 28 by tests, F24 by the md only, as logged.
  - Of my 18: 10 killed by tests, 2 by the md only, 6 survivors. Three survivors (G02, G15, G16) would change the published A2P on regeneration.
- **No major issue.**

**Findings.**
1. **D3dy-1 (minor): A2P's "distribution" is round 125's sample, and its test was formulated on that sample.**
   - The rtol 10⁻¹² cells reproduce round 125's numbers (8.95 ± 0.62 against 8.98 ± 0.63; 7.58 ± 0.57; the upward counts 5, 5 and 6, 3). The rtol 10⁻¹⁰ cells are new.
   - Fresh seeds confirm the conclusion: per cell t = +1.6, −0.4, +0.1, −0.8 (n = 200), pooled exponent 1.004 ± 0.009.
2. **D3dy-2 (minor): the seed's prefactor depends on the integrator's step sequence, not on u alone.**
   - At 5A, cutting the steps 16-fold lowers ln(S/uσ*) by 0.9 (double) and 1.3 (long double). At 3A, cutting them 4-fold raises it by 0.7 in both.
   - From rtol 10⁻¹⁰ to 10⁻¹² (steps ×0.72–0.77) ln|S| falls by 0.23 ± 0.07, pooled.
   - At max_step h₀/4 (5A) the two arithmetics' prefactors differ by e^{0.77} (4.5σ).
   - So "O(1) in units of uσ*" and Λ/(11 ln 2) = 1.00 are properties of DOP853's default steps at these tolerances. The pair's mean steps agree to within 7%, which is why Λ is fair.
3. **D3dy-3 (cosmetic): sample extremes stated as if they were properties.**
   - "Ranges within 5.32 … 12.51" and "12–38% later" are the extremes of ten draws.
   - Fresh draws give Λ from 0.79 to 14.02 (12% of them below 5.32) and lateness from 1.9% to 41.6%.
   - The test asserts 0.1 < lateness < 0.4 on the stored ten.
4. **D3dy-4 (cosmetic): the md's normal-hyperbolicity verdicts are hard-coded labels.** "Normally hyperbolic" and "not normally hyperbolic" are not computed from the rates printed beside them.
5. **D3dy-5 (cosmetic): "apart by 2.1·10⁻⁶ and 8.0·10⁻⁷" is relative.** The absolute differences in s are 1.0·10⁻⁷ and 2.8·10⁻⁹.
6. **D3dy-6 (minor): test gaps in A2P's regeneration path.**
   - G02 (the job's double exits at rtol 10⁻¹⁰ whatever its rtol) survives the whole suite.
   - G15 (RTOL_A2P = (1e-10, 1e-11)) and G16 (WMUL_A2P = (3, 4)) also survive the whole suite.
   - Each would silently change the published A2P on regeneration: the job is exercised live only at 5A, rtol 10⁻¹⁰, and nothing ties `part_A2P`'s configuration to the raw's.
   - G09 and G12 are guarded only by the md.

## What I checked and how

Scripts are in W/own/: rawp.py, a2p_re.py, fresh.py, stats.py, corr.py, hstep.py, steps.py, weight.py, selfseed.py, nh.py, apart.py, g04.py, spec130.py, mutants_tool.py.
- Data: fresh130.json, fresh2026.json, hstep5.json, hstep3.json, nh.json.
- Logs: W/logs/.
- Mutant records: W/mut/r125/res_a.jsonl, W/mut/own/res_b.jsonl, W/mut/own/res_c.jsonl.

1. **Reproducibility** [own]. `--from-raw` takes 30.9 s; `cmp` against the committed md: identical. pytest: 56 passed in 274.9 s. The slowest tests are the round-106 seed test (110 s), `test_a2e_…` (61 s) and the md byte test (31 s).

   `part_A2P` (a2p_re.py) takes 31.0 s:
   - every (s, ±1) exit equals the raw's;
   - max |Δ ln S|, |ΔΛ| = 0.0;
   - `a2p_line` is identical.
2. **Arithmetic-equivalence (code).**
   - `BumpsCW.corr` and `_LDBumps.corr` are `rP.sum(0)/N`: no endpoint weights, and no spacing is read from the permuted xs (`d` comes from XS).
   - f* is computed elementwise from the permuted xs in both functions.
   - On the exactly symmetric samples, u², P and r for ±x are bitwise equal and u/s² is bitwise negated. So Σ is exactly 0 in exact arithmetic, and the seed is the association's rounding: exact sums never leave, and reversing the order mirrors the run [md].
   - The only non-arithmetic difference within a pair is atol (10⁻¹⁴ for the double run, 10⁻¹⁴σ* for the long-double one). Steps (steps.py):

     | | rtol 10⁻¹⁰: steps (d / ld) | mean h (d / ld) | rtol 10⁻¹²: steps (d / ld) | mean h (d / ld) |
     |---|---|---|---|---|
     | 3A | 59 / 71 | 8.60 / 8.58·10⁻⁴ | 78 / 99 | 6.58 / 6.49·10⁻⁴ |
     | 5A | 56 / 66 | 7.09 / 6.84·10⁻⁵ | 74 / 94 | 5.22 / 4.86·10⁻⁵ |

     The mean steps agree to within 0–7%.
3. **Λ along the line** (selfseed.py; four fresh orders × two arithmetics × 3A, 5A) [own].
   - c is nonzero from the first step: unlike the code's order, which stays at exactly 0 to 4.8·10⁻⁴.
   - S(s) = c(s)e^{−Λ(0,s)} at the last two steps is 0.9973–0.9998 of the exit-based S. The 0.1–0.3% is the exit level's nonlinearity.
   - S(s) settles within 10% from s = 0.69–1.11·10⁻³ (5A) and 5.8–9.7·10⁻³ (3A). The turning points are 6.28·10⁻⁴ and 3.21·10⁻³, with Λ minima of −0.215 and −0.619 (weight.py). For white injection, the variance density e^{−2Λ} has its quartiles at 1.8, 2.8 and 3.7·10⁻³ (3A).
   - So the seed forms around and after the loss of stability, a slow passage. Λ along the line is the right integrand, and the code's order agrees with `a2_sym_lag` to 7.8·10⁻³ [md].
4. **Fresh draws** (fresh.py, the repository's `a2p_job` with seed 130 or 2026, k = 101, row 0 dropped; stats.py, corr.py) [own].

   | cell | stored (seed 125, n = 10) | seed 130 (n = 100) | seed 2026 (n = 100) | both (n = 200): t; gm d / ld (uσ*); lateness |
   |---|---|---|---|---|
   | 3A, 10⁻¹⁰ | 7.01 ± 0.40 | 7.63 ± 0.23 | 8.11 ± 0.21 | +1.6; 1.20 / 0.94; 1.9–38.5% |
   | 3A, 10⁻¹² | 8.95 ± 0.62 | 7.65 ± 0.24 | 7.48 ± 0.22 | −0.4; 0.75 / 0.80; 4.2–41.6% |
   | 5A, 10⁻¹⁰ | 7.42 ± 0.39 | 7.78 ± 0.17 | 7.50 ± 0.16 | +0.1; 0.68 / 0.67; 8.1–30.3% |
   | 5A, 10⁻¹² | 7.58 ± 0.57 | 7.72 ± 0.16 | 7.35 ± 0.18 | −0.8; 0.53 / 0.58; 8.1–30.2% |

   - Pooled over 800: Λ = 7.652 ± 0.070, against 11 ln 2 = 7.625: exponent 1.0036 ± 0.0091.
   - The largest |t| among the eight fresh cells is 2.35 (seed 2026, 3A 10⁻¹⁰), as chance predicts.
   - The difference between the arithmetics' log-prefactors is +0.25 ± 0.16, −0.06 ± 0.17, +0.02 ± 0.12 and −0.09 ± 0.12 by cell. So the 1.20 against 0.32 behind the +2.1 sem is a fluctuation.
   - The 3-sem rule fails in 0 of 80 blocks of ten (about 1.5% per cell is nominal for t₉). Its resolution with ten draws is ±1.2–1.9 in Λ, an exponent between about 0.75 and 1.25: it tells "scales" from "does not scale", not more.
   - Upward exits: 44–56%.
   - Per order, the double and long-double draws are uncorrelated (r between −0.21 and +0.24). The same order's draws at the two tolerances are also uncorrelated (r between −0.19 and +0.08; directions agree in 40/100 and 44/100). So the rtol 10⁻¹⁰ cells are independent new draws, and the rtol 10⁻¹² cells are round 125's draws.
   - Round 125's W/../ref125/own/perm.py uses `rng = np.random.default_rng(125)` with `perms = [np.arange(…)] + [rng.permutation(…) …]`, the same construction as `a2p_job`.
   - ln S has sd 1.1–1.85 and skew −0.3 … −1.4. This is not ln|N| (sd 1.111, skew −1.536 [own]; KS p ≤ 0.06): a scale mixture, broader at 3A.
5. **Step sequence** (hstep.py: the repository's exits with a max_step, ln|S| from `a2_line_lam`, 40 orders, seed 77) [own].

   | | double: mean h; ln(S/u_dσ*) (sd) | long double: mean h; ln(S/u_ldσ*) (sd) |
   |---|---|---|
   | 5A, default | 5.14·10⁻⁵; −0.59 ± 0.17 (1.09) | 4.79·10⁻⁵; −0.27 ± 0.16 (0.99) |
   | 5A, h₀/4 | 1.25·10⁻⁵; −1.02 ± 0.15 (0.98) | 1.25·10⁻⁵; −1.79 ± 0.08 (0.52) |
   | 5A, h₀/16 | 3.12·10⁻⁶; −1.47 ± 0.08 (0.49) | 3.12·10⁻⁶; −1.54 ± 0.07 (0.42) |
   | 3A, default | 6.85·10⁻⁴; −0.73 ± 0.34 (2.12) | 6.50·10⁻⁴; −0.08 ± 0.25 (1.61) |
   | 3A, h₀/4 | 1.85·10⁻⁴; −0.01 ± 0.26 (1.62) | 1.84·10⁻⁴; +0.61 ± 0.17 (1.07) |

   - The prefactor depends on the steps, in either direction: no simple h^{1/2} law. The directions stay a coin (15–25 of 40).
   - Across the arithmetics, at equal steps, the prefactors differ by −0.32, +0.77 (4.5σ), +0.07, −0.65 and −0.62. So "Λ ≈ 11 ln 2" holds to exponents of about 0.9–1.1 away from the default steps, and to 1.004 ± 0.009 at them.
6. **Normal hyperbolicity** (nh.py: my closed-form 3×3 Jacobian in long double at A2E's 46 on-line states) [own].
   - Against `critical_jac_nd`: max relative difference 1.9·10⁻¹⁵. λ agrees to 2.2·10⁻¹⁴, the fast tangential eigenvalue to 1.9·10⁻¹⁵, the slow one to 1.5·10⁻⁹.
   - The couplings ∂ċ/∂(α, σ) and ∂(α̇, σ̇)/∂c are ≤ 2.7·10⁻¹⁰ on A's samples and ≤ 1.6·10⁻¹⁶ on symmetric ones. The eigenvalues are the same on both.
   - Eigenvectors in the scaled coordinates (α/α, σ/σ): the fast eigenvalue is σ's (share 1.000). The slow one is α's rate (J_αα −0.00134 against −0.00124).
   - Rates, normal against tangential:

     | phase | 5A: normal | 5A: tangential | 3A: normal | 3A: tangential |
     |---|---|---|---|---|
     | repelling | 1.05–1.09·10⁴ | −1.446 … −1.477·10⁴ and −0.0012 | 629–665 | −2875 … −2961 and −0.00062 |
     | attracting | −279 … −383 | −644 … −886 and +0.013 … +0.022 | −24 … −230 | −338 … −1091 and +0.0003 … +0.0063 |

   - HPS conditions [memory]: an unstable normal bundle needs λ_u > k·max tangential, and a stable one needs λ_s < k·min tangential, for k = 0 … r.
     - In the repelling phase the max tangential rate is negative, so the condition holds for every r: normally hyperbolic.
     - In the attracting phase −279 < −886 fails at k = 1: not normally hyperbolic.
7. **D59 against the md.** Every number in the round-125 text is the md's and verbatim:
   - the Remark: 7.01 ± 0.40 … 7.58 ± 0.57; 5.32 … 12.51; 0.32–3.30; 5–6 and 3–5 of 10; 12–38%; 0.048297 and 0.0034820; 2.1·10⁻⁶ and 8.0·10⁻⁷;
   - A2E: −1.6 … +2.1 sem; the eight geometric means; 5, 5, 6, 6 and 5, 5, 4, 3; −886 … −644; −383 … −279;
   - the held states: 628.9, −0.00062, −2875; 1.05·10⁴, −0.0012, −1.446·10⁴.

   I recomputed the exits (apart.py): 0.04829717/0.04829727 and 0.0034820483/0.0034820511. The relative differences are 2.12·10⁻⁶ and 7.99·10⁻⁷; the absolute ones 1.02·10⁻⁷ and 2.78·10⁻⁹. "(ii)(a)" now starts its own paragraph.
8. **Mutants.** Tests first with `-k "round125 or round112 or round106_the_exit_s_conclusion"`, then the survivors against the whole suite with `-k "not test_round87_split_mass_bdf_and_the_md"`.
   - **The repository's spec:**
     - F00 is a control: it passes.
     - F01–F23 and F25–F29 are killed by tests.
     - F24 (`--from-raw` without A2P) is killed only by the md.
     - This matches the log.
   - **Mine:**
     - G00 (control) passes.
     - Killed by tests: G01 (seed 126), G03 (the job's long double at 10⁻¹²), G05 (sem over √(n − 1)), G06 (arithmetic mean for geometric), G07 (`any` for `all`), G08 (only the double's geometric mean checked), G10 (the eigenvalues over the repelling outputs), G13 (4 digits), G14 (u = eps), G18 (the paired Λ broken).
     - Killed only by the md: G09 (`min(dl)` printed), G12 ("apart by" absolute).
     - Survivors of the whole suite:
       - G04 (the line's Λ at 10⁻¹²): equivalent, ΔΛ ≤ 6·10⁻⁸ (g04.py).
       - G11 (the first exit for the mean): equivalent at five digits.
       - G17 (4 sem): equivalent. The +3.0 plant sits 5.2 sem out, so the threshold is pinned only within (2.14, 5.2).
       - G02, G15, G16: no published number changes now, but each changes A2P on regeneration (D3dy-6).
   - grep shows RTOL_A2P and WMUL_A2P are read only by `part_A2P`, which no test calls.
9. **Literature.** Five searches; no page fetched.

## Verdicts per part

| claim (D59 unless noted) | verdict | evidence |
|---|---|---|
| "over ten summation orders of the same data (the samples permuted)" | **holds**: plain means, f* elementwise, ±x terms bitwise antisymmetric. The same ODE up to atol (steps within 7%) | code, [own] |
| "Λ, the log-ratio of the two seeds" (Λ along the line, `a2_line_lam`) | **holds**: the runs' own seeds equal the exit-based ones to 0.997–0.9998; agreement with `a2_sym_lag` to 7.8·10⁻³ | [own], [md] |
| "Λ … 7.01 ± 0.40, 8.95 ± 0.62, 7.42 ± 0.39, 7.58 ± 0.57 … each mean within 3 sem of 11 ln 2" | **numbers verbatim; the conclusion holds and is sharper on fresh draws** (1.004 ± 0.009). The 10⁻¹² cells are round 125's own draws: post hoc (D3dy-1) | [raw], [own] |
| the 3A 10⁻¹² cell's +2.1 sem (1.20 against 0.32 uσ*) | **chance**: −0.35 sem and 0.75 against 0.80 on 200 fresh draws | [own] |
| "ranges within 5.32 … 12.51"; "12–38% later" | **true of the ten; sample extremes**: 0.79 … 14.02 and 1.9–41.6% fresh (D3dy-3) | [own] |
| "upward in 5–6 of 10 … 3–5 of 10" | **holds; a coin** (44–56% fresh) | [raw], [own] |
| "each arithmetic's seeds O(1) in units of … uσ* (0.32–3.30)" | **holds** (0.50–1.49 fresh; 0.21–1.8 at every step setting). **The prefactor is the integrator's sampling too** (D3dy-2) | [own] |
| "the seed scales with the arithmetic's roundoff" | **holds** at the paired, equal-step settings | [own] |
| "(the code's order one draw — …13–26% later, Λ 5.91 … 8.59)" | **holds** | [md] |
| A's long-double exits "0.048297 and 0.0034820, apart by 2.1·10⁻⁶ and 8.0·10⁻⁷" | **holds as relative differences** (D3dy-5) | [own] |
| "within 3·10⁻⁵ (a 10⁻³ change of S; the floor … no tolerance moves)" | **cited from round 125**, not rerun; the bound itself is live-tested | [md] |
| Fenichel: normally hyperbolic in the repelling phase; in the attracting phase the width contracts faster than c; a locally invariant manifold | **holds**: the Jacobian is exact to 10⁻¹⁵, and the HPS conditions hold for every r in the repelling phase and fail at k = 1 in the attracting one. The md's verdicts are hard-coded (D3dy-4) | [own], [snippet] |
| A2E's "σ's rate −886 … −644" | **holds**: that eigenvector is σ's (share 1.000) | [own] |
| (ii)(a) on its own paragraph; "the samples' and f*'s roundoff"; "(scipy, rtol 1e-10 unless stated)" | **holds** | [md] |
| What is proved and what numerical | **unchanged and correctly labelled** | re-read |
| tests | **good**: 56 pass and the repository's 29 mutants are killed. A2P's regeneration path is untested at 10⁻¹² and in its configuration (D3dy-6) | [own] |

## Status of the round-125 issues

| D3dt | status | what remains |
|---|---|---|
| 1 (single draws) | **applied**: A2P; mean ± sem, the seeds in uσ*, the counts, the code's order called one draw; tests on the distribution | The draws are round 125's own (D3dy-1). Ranges and lateness are sample extremes (D3dy-3). The prefactor's step dependence (D3dy-2) |
| 2 (the latent crash) | **resolved**: four return values; the plant passes | — |
| 3 (the resolution's cause, a noise digit) | **resolved**: five digits, the spread, the floor named | "apart by" is relative (D3dy-5) |
| 4 (the 40-digit citation) | **resolved** | — |
| 5 (Fenichel) | **resolved**: `a2_tangent` exact; the reading right | The verdicts are hard-coded (D3dy-4) |
| 6 (wording) | **resolved**: (a)–(e) all applied | — |
| 7 (md-only kills) | **resolved**: the semantic asserts kill F25–F29; the log records W09 and W12–W14 | G09 and G12 are md-only; G02, G15 and G16 survive (D3dy-6) |

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3dy-1 | minor | **A2P's draws replicate round 125's own; its test is post hoc on them.** <br>• `SEED_A2P = 125` with the identity order first is round 125's perm.py. The 10⁻¹² cells reproduce its table (8.95 ± 0.62 against 8.98 ± 0.63; 7.58 ± 0.57; up-counts 5, 5 and 6, 3) [raw], [own]. <br>• The 3-sem and (0.05, 20) rules were proposed on these draws. <br>• Fresh seeds 130 and 2026 (n = 200 per cell) give t = +1.6, −0.4, +0.1, −0.8; pooled Λ = 7.652 ± 0.070 = (1.004 ± 0.009)·11 ln 2 [own]. <br>• With ten draws the rule resolves the exponent only to about ±0.25. | **Code:** `SEED_A2P` set to a fresh seed and `K_A2P` to 40–100 (the cost is about 3 min for 100); a second seed as a replication. **md (`a2p_line`):** print the pooled exponent Λ/(11 ln 2) ± sem. **D59 Remark and A2E:** "(the referee's seed reproduces round 125's draws; on fresh orders …)". |
| D3dy-2 | minor | **The seed's prefactor is the arithmetic's rounding as sampled by DOP853's steps.** <br>• 5A: from the default steps to h₀/16, ln(S/uσ*) moves −0.59 → −1.47 (double) and −0.27 → −1.54 (long double). 3A: from default to h₀/4, −0.73 → −0.01 and −0.08 → +0.61. <br>• rtol 10⁻¹⁰ → 10⁻¹² (h ×0.72–0.77): −0.23 ± 0.07 pooled. <br>• At equal steps the arithmetics' prefactors differ by up to e^{0.77} (5A, h₀/4, 4.5σ) [own]. <br>• Λ is fair at A2P's settings because a pair's mean steps agree to within 7%. | **D59 Remark and A2E:** "the seeds O(1)uσ* at DOP853's steps (rtol 10⁻¹⁰–10⁻¹²); the prefactor moves by e^{±1} with the step sequence, the paired arithmetics sharing it". **md:** state the steps per arm. Optionally an A2P arm with a max_step. |
| D3dy-3 | cosmetic | **Sample extremes read as properties.** "Ranges within 5.32 … 12.51", "12–38% later" and the test's `0.1 < late < 0.4` are the extremes of ten draws. Fresh: Λ 0.79 … 14.02 (12% below 5.32, 1.3% above 12.51); lateness 1.9–41.6%, median 19–23% [own]. | **md and D59:** print the median and the interquartile range, or "(of ten)". **Test:** mark the lateness window as a regression on the stored draws, or drop it. |
| D3dy-4 | cosmetic | **Hard-coded verdicts.** "normally hyperbolic" and "not normally hyperbolic, … width contracting faster than c" are literals in `a2_lines`, with no comparison behind them. The tests assert the rates, not the verdict. | **`a2_lines`:** compute them. Repelling: λ > max tangential. Attracting: λ < min tangential, else "not normally hyperbolic". |
| D3dy-5 | cosmetic | **"apart by 2.1e-06 / 8.0e-07" is \|a/b − 1\|.** The absolute differences are 1.02·10⁻⁷ and 2.78·10⁻⁹ in s [own]. G12 (absolute) is caught only by the md. | **`lda` and D59:** "(apart by 2.1·10⁻⁶ relatively)". |
| D3dy-6 | minor | **A2P's regeneration path is untested.** G02 (the job's double exits at 10⁻¹⁰ whatever its rtol), G15 (RTOL_A2P = (1e-10, 1e-11)) and G16 (WMUL_A2P = (3, 4)) survive the whole suite. Each would change the published A2P on regeneration, because `a2p_job` is exercised live only at (5A, 10⁻¹⁰, k = 2), and no test ties `part_A2P`'s configuration to the raw's. G09 (`min(dl)`) and G12 are killed only by the md. G04, G11 and G17 are equivalent. | **Test:** `a2p_job(0.25, 3.0, 1e-12, k=2)` against the stored rows (about 5 s); `[(r['wmul'], r['rtol'], r['k'], r['seed']) for r in raw['A2P']] == [(w, q, K_A2P, SEED_A2P) for w in WMUL_A2P for q in RTOL_A2P]`; assert the printed "within {max(dl):.1e}" and "apart by" against the data. **Log:** record G09 and G12 as md-only. |

## Sources

- [snippet] Almada Monter and Bakhtin, "Normal forms approach to diffusion near hyperbolic equilibria" ([arXiv:1006.3000](https://arxiv.org/pdf/1006.3000)): under small noise near a hyperbolic point, the exit time and exit distribution have a joint scaling limit. Kifer (1981), as summarized in that search result, gives the exit time ~ λ₊⁻¹ ln ε⁻¹.
  - Related titles seen: "Long exit times near a repelling equilibrium" ([arXiv:1908.11840](https://arxiv.org/pdf/1908.11840)); Bakhtin and Pajor-Gyulai, "Tails of exit times from unstable equilibria on the line" ([J. Appl. Probab.](https://resolve.cambridge.org/core/journals/journal-of-applied-probability/article/abs/tails-of-exit-times-from-unstable-equilibria-on-the-line/D123EC7CC1FBD51CB04D464062F0648F)).
  - These are the theory behind A2P: Λ = ln(u_d/u_ld) plus the difference of two O(1) random terms.
- [snippet] Berglund and Gentz, "Pathwise description of dynamic pitchfork bifurcations with additive noise", Probab. Theory Relat. Fields 122 (2002) ([HAL](https://hal.archives-ouvertes.fr/hal-00130541)). Lythe ([adap-org/9707008](https://arxiv.org/abs/adap-org/9707008)): in a slow passage through a loss of stability, the dynamics are controlled by additive disturbances "even if they are extremely small". This matches item 3: the seed forms around the λ sign change.
- [snippet] Henrici's model of rounding errors as independent random variables, as described by Mosbach and Turner, "A quantitative probabilistic investigation into the accumulation of rounding errors in numerical ODE solution" ([arXiv:math/0512364](https://arxiv.org/pdf/math/0512364)), who single out saddle points. On Brouwer's law (roundoff as a random walk, sd ∝ √steps) and its limits: [arXiv:1702.03354](https://arxiv.org/pdf/1702.03354). These are the basis of D3dy-2. The step dependence I measured is not the simple h^{1/2}, consistent with the caveat that rounding errors are correlated.
- [snippet] Fenichel's theorem with boundary, where the perturbed manifold is locally invariant: [arXiv:0911.4990](https://arxiv.org/pdf/0911.4990), [arXiv:2012.06770](https://arxiv.org/pdf/2012.06770).
- [memory] The Hirsch–Pugh–Shub r-normal-hyperbolicity rate conditions. Golubitsky–Stewart–Schaeffer, as cited by D59.

## Lesson

A replicated draw is not a new draw.
- A test proposed on a sample confirms nothing on that sample. Run it again under a fresh seed before calling it a distribution.
- Here the fresh 800 draws turned a post hoc "within 3 sem" into an exponent of 1.004 ± 0.009.

A seed made by rounding belongs to the arithmetic and the integrator together.
- Compare arithmetics at equal steps, as A2P does, and say so.
- The prefactor in units of uσ* moved by e^{±1} when only the steps changed.
