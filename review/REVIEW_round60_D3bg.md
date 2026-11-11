# Referee report D3bg, round 60: Proposition D‴.58 as revised after round 54

I reviewed commit fef9c70. During the review, other sessions moved HEAD to ada0e1c, but `git diff --stat fef9c70 HEAD` is empty for every reviewed file (the draft, the three source files, the experiment, the md and raw, and the three test files). The review therefore applies unchanged at HEAD.

## Summary

**What checks out.**
- **Reproduction.**
  - `--from-raw` into my own `WM_OUT` regenerates `results/chaos/width_mismatch.md` byte for byte (md5 a50774c0…), and the raw is unchanged.
  - All 45 tests pass (11 + 12 + 22).
  - Three B paths, seven B2 jobs and four A rows rerun bit-identically: events, τ, errors, optima, polish logs and `entry_info`.
- **(i), (ii)'s orderings and (v).**
  - I re-derived the proofs, including the three gaps D3ba-9 asked to close:
    - the atoms' bound via Cramér;
    - the two-dimensional domination κ_ϑ ≤ (1 − ϑ²)^{−1/2}e^{a(z²+z′²)}, a = ϑ/(1 + ϑ) < ½, which is tight on the diagonal;
    - (v)'s measure on the imaginary axis.
  - q > t is correctly proved. The Fock and Mercer-zero variances are the closed forms the draft states.
  - Spot checks at untested ratios (τ/σ = 3 and 0.1; σ*/σ = 0.3, 0.65, 0.95) reach the double-precision floor.
- **(iv).**
  - `stop_law`'s expectation s²(b_K² + 1 − 2K/N) and scatter s²(2(N − 2K)/N² + 4b_K²/N)^{1/2} are the right first-order moments.
  - At the deciding K, the 10 draws per level sit at a mean z of −0.08 … +0.28 (sd 0.69–1.31).
  - The rms falls monotonically along all 30 noisy paths, as the law requires.
- **Quoted numbers.** Almost every number in (iii), (iv) and the numerics matches the md and raw; the exceptions are listed below.
  - "19 of 19" is not a tautology: VP from the block and VP from the rule agree to ≤ 1.03·10⁻⁴ in all 19 cases, and to ≤ 10⁻⁶ in 14.
  - C(c) ≥ 0 holds after all 258 optimal blocks of the five variants, to rounding.
  - The entrants' duals rebuilt from the raw match `entry_info` to 1·10⁻¹¹.
- **Earlier issues.** D3ba-3, -6, -8 and -9 are fixed.

**What is wrong, majors first.**
1. **D3bg-1 (major): the new explanation of where the optimal ranges end is refuted.** (iii) says the repeated polish "stopped on its gradient test 1.012, 2.677 and 1.121 times above the optimum: the joint least squares in (w, c) loses its descent where the optimum's error is 10⁻⁹–10⁻¹⁰, which variable projection, eliminating w, does not."
   - The "gradient test" is `polish_block`'s absolute `gtol = 1e-15` (rich_transit.py:218).
   - From the same three blocks, the same TRF with only that test disabled reaches the optimum.
   - Rerunning B's repeated variant with gtol off extends the ranges from 4 to 5 at τ/σ = ¼ and from 9 to 10 at ¾.
   - This is the same error class as D3ba-1: a solver setting presented as a property of the problem.
2. **D3bg-2 (major): "its later blocks leave the optimum's basin (variable projection reaches it from 13 of the 83)" is judged against optima below the experiment's own floor.** 65 of the 70 misses have reference optima below 10⁻¹⁰ (E_FLOOR), down to 9.7·10⁻¹⁷. Among the 18 later blocks with an optimum ≥ 10⁻¹⁰, VP reaches it from 13.
3. **D3bg-3 (major): 21 of 36 mutants of the new code survive.**
   - The survivors include the code behind "19 of 19", "68 of 390", "0.30", "339 of 360" and the P(stop) cells.
   - Every mutant of the repolish, budget, entry-log and count code in `sliding_iss` survives.

Minors:
- "95.7–100%" should be 93.0–100% (D3ba-5 not fixed).
- "not the stop" (D3ba-7 half fixed).
- The entry statistics pool first entries and non-optimal blocks.
- The "102 of 102 passed" test is vacuous for 39 of those blocks.
- "converged 1.2%".
- "within 10⁻⁶σ" holds for the base variant only.
- The "12" ranges are censored by the path's cap.
- Literature: the Gaussian-RKHS ≅ Hermite-space isomorphism literature is missing; the KW discrepancy is resolvable in one line; the Derevianko–Plonka comparison reduces to theirs directly.

## What I checked and how

- **Setup.**
  - I worked in a scratch copy (`git archive fef9c70`) under `scratchpad/ref60/D3bg_work`. I ran mutants in two separate copies of that copy and wrote nothing into the repository.
  - All processes used the venv with OMP_NUM_THREADS=1, nice 10 and PYTHONDONTWRITEBYTECODE=1, at most two heavy processes at a time.
  - The only process I stopped was my own monitor. Another referee's processes (ref61) were running and I left them alone.
  - mpmath was not needed. Extended-precision checks used `np.longdouble`.
- **Reproduction.**
  - `--from-raw` with WM_OUT in scratch; the three test files (45 tests).
  - `b_job` for (2, base), (1, polish 2σ) and (¼, 5·10⁴ repeated).
  - `b2_job` for B2 entries 0, 1, 3, 13, 17, 40 and 77.
  - A rows (½, 3), (1, 7), (2, 10) and (¼, 5): e_hs, e_fock, e_mercer, e_best, v_best and e_opt.
- **Numbers.** I checked every number in (iii), (iv), the status and the numerics against the md and the raw.
  - I recomputed δ_t(K), the gap shares, the node deficits and the B and B2 tallies.
  - Entry statistics were split by regime (first entries; entries from optimal blocks; entries from non-optimal blocks).
  - C(c) was computed on a 32001-point grid after all 258 optimal blocks.
  - The duals were rebuilt from the blocks and event times, and their charge attributed to "before the first entry", "within ±0.75σ of the current block's flank peak", and so on.
- **Where the ranges end.**
  - From the first non-optimal blocks of 9 paths I ran:
    - TRF exactly as polished;
    - TRF with gtol disabled (`gtol=None`), everything else equal;
    - TRF with `x_scale='jac'`;
    - LM on the joint (w, c), unbounded, with the analytic Jacobian;
    - VP-LM.
  - I recorded ‖g·v‖∞ (scipy's Coleman–Li-scaled gradient) and its rounding error (an 80-bit recomputation), and cond J for the joint and VP Jacobians.
  - I reran B's `b_job` for 'polish 2σ, 5·10⁴, repeated' at τ/σ = ¼, ½ and ¾ with only `polish_block`'s gtol disabled (monkeypatched).
  - I checked where the rejected polishes would land (τ/σ = 1 grid; ½ and ¾ base).
- **Later blocks.** I categorised B2's 70 later misses by reference optimum, distance and near-miss, and ran joint LM (6·10⁴ evaluations) from all 70.
- **Stationarity test.** I applied a Cauchy–Schwarz bound to `sliding_iss`'s test at all 102 non-optimal blocks, and recomputed the exact test at all 360 blocks (all pass).
- **(iv).** Per-draw rms²/s² against `stop_law`'s mean and sd for K = 1…8 at the three levels, with e_K both the rule's and B's optimum. I also checked the stops, the at-feature counts and the gap index per draw.
- **(i), (ii), (v).**
  - Step by step: Bargmann prefactors and exponents, the Fock Gaussian integral as Mehler's kernel at ρ = t, c_T² = ‖f*‖²(1 − t²)^{1/2}, D_{2K,K} = 2^K K!.
  - The atoms' bound E He_n(√2Z)² = n!²Σ_a 2^{n−2a}/((n − 2a)!4^a a!²) ≤ n!(3ⁿ + 1)/2.
  - The domination exponent a = ϑ/(1 + ϑ), and the continuation (smallest eigenvalue (1 − ρ)/(4(1 + ρ)) > 0; series radius ≥ 1).
  - (v) with v = iw, v′ = −iw′.
  - FM's Mercer variance, the Fock variance, λ/t = 1 + sech r, and q > t.
  - The classical rule (moment matching) at τ/σ = 2 and 3 to K = 30: ratios 1.4999 and 1.2222.
- **Mutation testing.** 36 single-point mutants:
  - 27 in `b2_job`, `part_B2`, `b2_line`, `polish_cell`, the B table and summary, `entry_lines`, `stop_law`, C's stop-law rows, `gap_share` and `gap_range`;
  - 9 in `sliding_iss` (repolish, `polish_nfev`, `entry_info`, counts, polish log).
  - Each mutant was run against the relevant test files and a `--from-raw` diff. For the 9 `sliding_iss` mutants I also ran three B paths, (¼ repeated), (2 repeated) and (¼ base), against the raw.
- **Literature.** Web-search snippets only.
  - WebFetch to ar5iv.labs.arxiv.org, ar5iv.arxiv.org, eprints.ncl.ac.uk and tskarvone.github.io failed (DNS: outside this session's network).
  - api.semanticscholar.org was refused by the proxy (403).
  - So nothing below is from full texts.

## Verdicts per part

| claim | verdict |
|---|---|
| (i) half-squeeze nodes and amplitudes; error = ‖f*‖²(1 − t²)^{1/2}Σ_{n≥2K}(tⁿ/n!)D²; D_{2K,K} = 2^K K!; lower bound; upper bound C_ϑ(t/ϑ)^K; K-th-root limit t | **Correct.** Re-derived; trapezoid = series to 10⁻¹⁶–10⁻¹⁵ at τ/σ = 3 (K ≤ 8) |
| (i) proof gaps (D3ba-9): atoms, two-dimensional domination, continuation | **Closed.** The two-dimensional Uspensky step is a sketch but standard (positive product rule, tails controlled by e^{bz²}, a < b < ½) |
| (ii) λ, λ/(2 − λ) = t/(1 − t + t²) and q exceed t; q > t via (x − 1)/(x + 1) at β² against σ*/σ; v_Fock = 2τ²σ²/(σ² + σ*²); v_Mercer = τ²(1 + 2τ²/σ²)^{−1/2}; λ/t = 1 + sech r | **Correct.** Rates numerical, as the status says |
| (iii) δ_t: 0.078→0.040%, 0.73→0.40%, 3.5→2.1%; nodes within 0.008–5.2% | **Correct** |
| (iii) "95.7–100% of the gap recovered" | **Wrong:** 93.0–100% (md A3) |
| (iii) C ≥ 0 with double zeros; flank peaks 2.0–2.6σ at K = 1 | **Correct.** After all 258 optimal blocks, to rounding |
| (iii) entrant = least time to the wall; 68/390; 0.30; 0.946; 339/360; 1.40σ; 21 between | **Numbers match**, but the first is the path's definition and the statistics pool regimes (D3bg-5) |
| (iii) "dual charged under the earlier blocks' flank peaks" | **Partly.** Median 53% of the charge within ±0.75σ of a flank peak; mean 30% before the first entry |
| (iii) optimal ranges: 5000 evaluations 4, 5, 7, 9, 12, 12 (one exception); 5·10⁴ 4, 7, 9, 11; repeated 4, 7, 9, 12 | **Numbers match.** The ranges remain solver-set: with gtol off, the repeated variant gives 5, 7, 10 at ¼, ½, ¾. "12" is the path's cap |
| (iii) VP from the first non-optimal block reaches the optimum (19/19) | **Verified** |
| (iii) "stopped on its gradient test 1.012, 2.677, 1.121" | **Numbers match.** At ½ and ¾ two polishes first spent 5·10⁴ evaluations each (status 0); the test is an absolute gtol |
| (iii) "the joint least squares in (w, c) loses its descent …, which VP does not" | **Refuted** (D3bg-1) |
| (iii) "later blocks leave the optimum's basin (13 of 83)" | **Unsupported** (D3bg-2) |
| (iii) "every non-optimal block passed the relaxation's test (102 of 102)" | **True but vacuous for 39** (D3bg-6) |
| (iv) expectation, scatter, P(stop at K), K(s) formula, bias–variance split | **Correct to first order;** consistent with C's draws |
| (iv) K(s) 5.26/3.91/2.55; P 0.993/0.789/0.986; 27.7 expected vs 26 observed; next P 0.001/0.210/0.012; 0.92/0.99/1.02; 0.57–1.35; ≤ 0.87; 1.10/2.66/1.93 | **All match** |
| (v) narrower features: cosine sum, series at \|t\| | **Correct.** Spot-checked at σ*/σ = 0.3, 0.65, 0.95 |
| Status | **Fair**, apart from (iii)'s mechanism sentence. "In the optimum's basin" means VP's basin |
| Numerics (A), (D) | **All numbers match** |
| Numerics (B) | **Match**, except "converged 1.2%" and "within 10⁻⁶σ of the optimum's" (base only) |
| Numerics (C) | **Match**, except "not the stop" |
| Literature / "New" | **Appropriately hedged on KOG**, but the Hermite-space isomorphism literature is missing (D3bg-8) |

### Status of the D3ba issues

| D3ba | status | note |
|---|---|---|
| 1 | partly | "Double precision/ε_mach" dropped; budgets varied; polish status recorded. But the replacement mechanism is again a solver setting (D3bg-1), and later blocks are judged against sub-floor optima (D3bg-2) |
| 2 | mostly | Least time to the wall stated and logged; statistics pooled across regimes (D3bg-5) |
| 3 | fixed | `if failed: continue` in `sliding_iss` (rich_transit.py:454) and `sliding_iss_nd`. Tested with a failing `relax_chunk` at the first and a later event, ±exit_kkt, 1-D and n-D |
| 4 | partly | Report cells now tested (15 of 36 mutants killed), but B2, the entry statistics, `stop_law`'s moments and all repolish/budget code are unprotected (D3bg-3) |
| 5 | not | "95.7–100%" against the md's 93.0%–100.0% |
| 6 | fixed | `stop_law`; formulas and probabilities verified |
| 7 | partly | "to 11 entries" fixed and the md prints "(the stop in 4)", but the draft still says "not the stop" |
| 8 | fixed | q > t proved; gain "up to a power of K" |
| 9 | fixed | Re-derived (above) |
| 10 | partly | Derevianko–Plonka (DP) citation and KW attribution fixed. KOG's scaling unchecked for a third round; Hermite-space literature missing |
| 11 | mostly | A's caption, "passed the relaxation's test" and the trapezoid attribution fixed. "The constant of (i) to within 1 − δ_t(K)" unchanged |

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3bg-1 | major | **The repeated polish's stops are polish_block's absolute gtol, not a loss of descent of the joint least squares (D58.md:40–42; numerics :146 "converged 1.2% above the optimum").** (a) The test is scipy TRF's ‖g·v‖∞ < gtol with `gtol=1e-15` (rich_transit.py:218). At the three blocks ‖g·v‖∞ = 5.2·10⁻¹⁶, 1.8·10⁻¹⁷ and 4.5·10⁻¹⁷; their rounding errors (80-bit recomputation) are only 2.3, 1.5 and 2.5·10⁻¹⁸, so the gradients are resolved. (b) From the same blocks, TRF with identical bounds, Jacobian and xtol/ftol and `gtol=None` reaches 1.000000 of the optimum at τ/σ = ¼, K = 5 (1530 evaluations), 1.00028 at ½, K = 8 (178024, xtol stop) and 0.99998 at ¾, K = 10 (52137). Joint LM reaches 1.000000 (515), 1.00004 (57969) and 0.99996 (19061). The base block at ¼, K = 5 (10.85×) is also stopped by gtol at 1.0016 after 25929 evaluations, and reaches 1.000000 with gtol off. (c) B's own `b_job`, repeated variant, with only gtol disabled: optimal to K = 5 (was 4) at ¼, with 1.6114·10⁻⁹ reached exactly; to K = 10 (was 9) at ¾; at ½ still 7, the K = 8 block at 1.233× instead of 2.677×, its four polishes now ending on the 5·10⁴ budget. (d) At ½ and ¾ the polishes logged "0/50000/a 0/50000/a 1/1/r": budget twice, then the test at the first evaluation. (e) What separates VP is conditioning, not descent: cond J(w, c) = 4.5·10⁷, 4.4·10⁸ and 2.0·10⁸ against cond J_VP = 1.2·10⁴, 2.3·10⁵ and 2.4·10⁵; VP-LM converges in 22, 75 and 36 evaluations. | (iii), (B), status: say the ranges end where the polish's budget, acceptance radius or absolute gradient tolerance binds. At e ≈ 10⁻⁹–10⁻¹⁰ the joint gradient is 10⁻¹⁶–10⁻¹⁷, so gtol = 10⁻¹⁵ trips above the optimum; the same solver without it, or VP, reaches the optimum. Code: make gtol relative (scale it by ‖r‖‖J‖, or rely on xtol/ftol), or polish by VP. Rerun B and quote the ranges as solver-limited. |
| D3bg-2 | major | **"Its later blocks leave the optimum's basin (variable projection reaches it from 13 of the 83)" (D58.md:42–43) uses optima below the experiment's own floor.** 65 of the 70 misses have reference optima below 10⁻¹⁰. That is E_FLOOR (exp_width_mismatch.py:53), below which the md's caption says VP "stops short of the optimum". 45 of the 70 are below 10⁻¹², down to 9.7·10⁻¹⁷. 6 misses end 1.0012–1.008 times the optimum, i.e. in its basin but unconverged. Of the 18 later blocks with an optimum ≥ 10⁻¹⁰, VP reaches 13, one more ends at 1.0014, and only 4 (τ/σ = ½, K = 8 and ¾, K = 10, base and grid) miss by 16–25×, each with a kernel 2.1–4.0σ from its optimal place. Joint LM from all 70 (6·10⁴ evaluations) ran out of budget in 52 and reached none, so it is inconclusive. B2 also only re-runs B's two VP candidates (bit-identical). | (iii), (B): restrict the count to optima ≥ E_FLOOR (or recompute them in extended precision). Say instead "VP reaches the optimum from 13 of 18 later blocks with an optimum ≥ 10⁻¹⁰; 4 hold a kernel ≥ 2σ from its place", or drop the sentence. |
| D3bg-3 | major | **21 of 36 single-point mutants of the new code survive the tests.** Survivors that change the committed md: M01 (B2's hit tolerance ×10, :205), M03 (B2's first block keyed by τ/σ only, :209), M15 (the "largest rate" tolerance 10⁻¹, :258; changes "68"), M18 (mean for median, :260; changes "0.30"), M19 (beyond/between at 0.5σ, :263; changes "339/21"), M20 (stop_law's scatter without the bias term, :375) and M21 (the mean with K/N, :374); the last two change the P(stop) cells. All 9 `sliding_iss` mutants survive. Seven change reruns of the B paths: M28 (a non-improving repolish accepted, :435), M30 (one extra polish, also in the base variant, :429), M31 (the first polish ignores `polish_nfev`, :410), M35 (no chunk reset after a repolish, :441), M32 and M36 (`entry_info`'s u0 and rate_max, :377–378) and M34 (the log's accepted flag, :415). M29 and M33 (failed chunks counted as redone, :399) do not change these paths, but no test returns how = 3 to the counter. M02, M05, M06, M07 and M25 are equivalent on the committed raw. Killed: M04, M08–M14, M16, M17, M22–M24, M26, M27. | Tests: repolish and `polish_nfev` on a small target (number of polishes, budget honoured, rejection of a non-improving polish, chunk reset); `entry_info` against duals rebuilt from the path (they agree to 10⁻¹¹ here); `b2_line` and `entry_lines` on asymmetric synthetic data; `stop_law` against closed forms; the counters with how = 3. |
| D3bg-4 | minor | **"95.7–100% of the gap recovered" (D58.md:26).** The md's A3 line prints 93.0%–100.0%: at τ/σ = 2 the share is 93.0–95.7% for K = 2–10, and at τ/σ = 1, K = 2 it is 97.7%. D3ba-5 is not fixed. | (iii): 93–100%. |
| D3bg-5 | minor | **The entry statistics pool regimes (D58.md:31–34).** "68 of 390" includes the 30 first entries, where all duals are zero and the largest rate is the least time to the wall by construction; it is 38 of the 360 later entries. After optimal blocks only 4 of 258 entrants had the largest rate (median 0.21 of it), against 34 of 102 after non-optimal blocks. "339 of 360 beyond, a median 1.40σ": after optimal blocks it is 256 of 258, median 1.19σ. 19 of the 21 between-node entrants follow non-optimal blocks (their duals 0.24–1.00); the 2 after optimal blocks are D3ba's (duals 0.992 and 0.989). "Least time to the wall" is `sliding_iss`'s entry rule (`dt.min()`), not a finding. | (iii): quote the optimal-block statistics; leave the pooled ones to (B). |
| D3bg-6 | minor | **"Every non-optimal block passed the relaxation's test (102 of 102)" (D58.md:43–44).** The test is absolute: \|⟨rφ_i⟩\|, \|w_i⟨r∂φ_i⟩\| ≤ 10⁻⁹·max\|⟨f*φ⟩\| (rich_transit.py:427). For 39 of the 102 blocks (e ≲ 4·10⁻¹⁰) it holds by Cauchy–Schwarz whatever the state, so it carries no information about stationarity there. | (iii): say so, or make the test relative to the residual. |
| D3bg-7 | minor | **Numerics wording.** (a) "Converged 1.2% above the optimum" (:146): it is the gtol stop (D3bg-1). (b) "Within 10⁻⁶σ of the optimum's" (:144) holds in the base variant (3.9·10⁻⁷σ); in the 5·10⁴ variants optimal blocks lie up to 1.05·10⁻³σ (τ/σ = 1, K = 12) and 4.3·10⁻⁵σ (½, K = 7) away. (c) "Not the stop" (:154): at s = 10⁻³ the largest gap is the stop in 4 of 10 draws. (d) The ranges "12" are K_B = 12, the path's cap, i.e. censored. (e) "C ≥ 0 to 10⁻⁶ relative" is base only; in the 5·10⁴ variants it is −8.1·10⁻⁶ relative, which is −5.6·10⁻¹⁶ absolute (rounding). | Numerics (B), (C); (iii). |
| D3bg-8 | minor | **Literature and novelty (D58.md:107–124), from snippets.** (a) Missing: the explicit isometric isomorphism, preserving point evaluations, between Gaussian-kernel RKHSs and Hermite spaces (Gnewuch–Hefter–Hinrichs–Ritter, J. Complexity 71 (2022)). Also missing: Gauss–Hermite (GH) integration in Hermite spaces of analytic functions with exponential rates (Irrgeher–Kritzer–Leobacher–Pillichshammer, J. Complexity 31 (2015)), and Gnewuch–Ritter–Rüßmann (arXiv:2412.05368). Any such isometry makes a GH rule's error an exact defect series, so "the frame … that makes its error an exact series" needs this comparison. By my own algebra (unverified against the papers), an isometry onto the Hermite space orthonormal under the integration measure itself gives the Mercer-zero variance and ratio q; at τ/σ = 1, ρ = 2 − √3 = q. The half-squeeze's other reference Gaussian may be the new element. (b) The KW/KS "discrepancy" resolves in one line. With weights ≥ 0 summing to 1 and k ≤ 1, the classical rule of the integration measure has e² ≤ 2, so KW's "exponential divergence for 2αγ² > 1" cannot be about it; it must be a GH rule of a fixed weight, reweighted. Here moment matching converges at λ (τ/σ = 2, 3: ratios 1.4999 and 1.2222 at K = 30). (c) (v)'s error is f* − f_K = e^{−x²/(2σ²)}(g − Σb_k cos(η_kx/σ²)), with g a Gaussian of width σσ*/(σ² − σ*²)^{1/2}. That is DP's weighted problem in L²(ℝ, e^{−x²/σ²}) (ρ = σ²/2), so the comparison is direct. "24(6)" should be article 6 of volume 24. (d) KOG's scaling is still unchecked. | Literature and novelty paragraph: add (a) and compare; state (b) and (c); keep "New" hedged. |
| D3bg-9 | cosmetic | K(s) denotes both the least integer and the root 5.26/3.91/2.55 (D58.md:50–54). "The constant of (i) to within 1 − δ_t(K)" (:25) is still ambiguous. `b2_line` prints "1–2.28e+06 times above it", where the 1 is 1.001 at `.3g` (exp_width_mismatch.py:214). B2 recomputes B's own two VP candidates. | Wording; the format string. |

Sources (search-result snippets only; full texts unreachable from this session):
- [Karvonen–Oates–Girolami, arXiv:2004.12654](https://arxiv.org/abs/2004.12654), [ResearchGate](https://www.researchgate.net/publication/340963491_Integration_in_reproducing_kernel_Hilbert_spaces_of_Gaussian_kernels), [Math. Comp. 90 (2021)](https://www.ams.org/journals/mcom/2021-90-331/S0025-5718-2021-03659-5/)
- [Karvonen–Särkkä, arXiv:1803.09532](https://arxiv.org/pdf/1803.09532), [BIT 59 (2019)](https://link.springer.com/article/10.1007/s10543-019-00758-3)
- [Kuo–Woźniakowski, BIT 52 (2012)](https://link.springer.com/article/10.1007/s10543-011-0358-9)
- [Derevianko–Plonka, STSIP 24 (2026)](https://link.springer.com/article/10.1007/s43670-026-00124-9), [arXiv:2307.13587](https://arxiv.org/abs/2307.13587)
- [Gnewuch–Hefter–Hinrichs–Ritter, arXiv:2110.05778](https://arxiv.org/pdf/2110.05778)
- [Irrgeher–Kritzer–Leobacher–Pillichshammer, arXiv:1403.5102](https://arxiv.org/pdf/1403.5102)
- [Gnewuch–Ritter–Rüßmann, arXiv:2412.05368](https://arxiv.org/pdf/2412.05368)
- [Szabó–Adam–Janszky–Domokos, PRA 53, 2698 (1996)](https://doi.org/10.1103/PhysRevA.53.2698) (citation only)

## Lesson

A stopping test is calibrated for a scale. An absolute gradient tolerance of 10⁻¹⁵ is met along a slow mode of any fit whose residual is 10⁻¹⁰, so there it certifies nothing. The revision replaced "double precision" with "the joint least squares loses its descent" without the check round 54 named: rerun the same solver from the same state with only the suspect setting changed. That check took two seconds per block and moved two of the three ranges. Likewise, a miss measured against an optimum below the experiment's own floor (E_FLOOR) is evidence about the solver, not about basins. And a statistic pooled over regimes (first entries, unfinished blocks) should not be quoted for the regime the sentence is about.
