# Referee report D3ar (round 45): Proposition D‴.56, the transit under measurement noise

## What I checked and how

**Setup.**
- I copied the repository, without `venv/` and `.git/`, to a scratch directory. Everything ran on the original `venv/bin/python` with
  `OMP_NUM_THREADS=1` and `nice -n 15`, one heavy job at a time, each under 15 minutes.
- **One rule slip, now undone:** a backgrounded `cd …&&` left my foreground shell in the repo, so one untracked scratch script
  (`anaJB.py`) landed in the repo root for under a minute. I moved it into the scratch directory. No tracked file was touched. The D‴.56
  objects are identical to HEAD, and no commit since f4c9a57 touches them.

**What I ran.**
1. **Numerics against the md and the raw file.** `--from-raw` regenerates `noisy_certificate.md` byte-identically. Every number in the
   numerics paragraph matches the md; the one wording slip is "52 of 100" at s = 0.3, which is 52 of 99 finished paths.
2. **Part B reproduced bit for bit** at s = 0.03, 0.1 and 0.3: 1200 paths, 0 mismatches in `exact` or `crit`. Per draw I also recorded
   the order of the first K entries, the noisy path's own largest off-E dual (via `path_duals`) and the `relaxed` flags.
3. **All 40 noiseless paths and their responses:**
   - entry gaps and margins;
   - trace R_K = 55 ± 10⁻¹⁴, projection error 10⁻¹⁵;
   - the Gaussian union bound evaluated per target;
   - a KS test of |R_Kξ|²/s² against χ²₅₅;
   - a PIT/KS test of the observed W·s against each target's limit law (Monte Carlo, 4·10⁴ samples).
4. **`path_response` against central finite differences** on six targets of my own, at h = 10⁻⁶ … 10⁻¹:
   - near ties of 0.86% and 1.1%;
   - a same-sign pair 0.30 apart, where ‖R₁ − R₁(no curvature)‖ = 2.8;
   - an opposite-sign pair 0.28 apart;
   - amplitudes 1 : 0.1;
   - the test file's own target.
5. **Mutation testing:** 11 mutants of `path_response`/`fit_hessian` and 8 of the experiment's classification and report. The survivors
   were also run against `test_dead_zone.py` and `test_rich_transit.py` (24 tests in all).
6. **Part E's raw file:** time-outs, the exact-only filter, short paths, Radau ("stiff") paths against rule failures, exits after the
   (K+1)-th entry, zero gaps. I reproduced three anomalous paths, with and without the polish.
7. **The two stopping rules on my own families**, each path run to K + 4 entries:
   - amplitudes (1, 0.3, 0.1): 5 targets × 3 levels × 4 draws;
   - K = 5 features: 4 targets × 2 levels × 3 draws;
   - Part E targets 0–2 at s = 0.1 to K + 8 entries.

   I also reran the one Part E path that went over the 180-s limit, without a limit.
8. **The oracle's bias:** 400 least-squares fits started at the truth, for 3 targets × 3 levels.
9. **Literature:** checked by web search. arxiv.org and ams.org could not be fetched, nor could the UCLA CAM report or deepai.org.
   ORXYY's stopping-time formula is from memory.

## Verdicts

**(i) The linear response — correct, but its radius is unstated.**
- The recursion is right and `path_response` implements it. The finite differences converge at second order everywhere, for example
  relative errors of 3.6·10⁻⁹, 10⁻⁷ and 10⁻⁵ at h = 10⁻⁴, 10⁻³ and 10⁻² on the 0.86% near tie. This includes blocks with large residual
  curvature.
- `path_response` correctly refuses paths with exits or simultaneous entries.
- What is missing is a quantitative "|ξ| small". The linear regime ends at h ≈ 0.01–0.03 (a unit-variance direction) on near ties, where
  the order of entries swaps, and on the close pair and the 1 : 10 amplitude target, where the support changes. On the test's own target
  it ends at 0.1.

**(ii) The support criterion — sound inside the regime, not "sharp" where the numerics use it.**
- It is exact at first order only while the noiseless order of entries is kept.
- 10 of the 40 targets have consecutive noiseless entries within 5% of each other (down to 0.86%). On exact draws the entries of E swap
  order on 23/396 at s = 0.03, 66/374 at 0.1 and 87/191 at 0.3.
- Swaps account for 1/1, 7/8 and 42/49 of the "exact but criterion ≥ 1" disagreements.
- On draws without a swap, the linear duals miss the path's own duals by at most 0.0015, 0.016 and 0.089. That is ∝ s², with coefficient
  ≈ 1.7, the same as Part A's. On swap draws the miss is up to 0.074, 0.19 and 0.39.
- "+O(s²)" in the probability bound is unproved and of the wrong form (see D3ar-2).
- The linear model's expected exact counts are rejected where failures actually happen:

  | s | observed exact | linear model | z |
  |---|---|---|---|
  | 0.1 | 374 | 359.9 ± 5.5 | 2.6 |
  | 0.3 | 191 | 159.9 ± 9.5 | 3.3 |

- "Sharp" can therefore only mean "the exact first-order boundary inside the linear regime". Per-draw agreement does not measure it: a
  wrong criterion (mutant M6) agrees better, on 395 and 367 draws against 388 and 338.

**(iii) The window — correct asymptotically, overstated at s ≥ 0.1.**
- W is exact given the K-th block, by construction of the path.
- The limit law of sW holds:
  - pooled limit law: median 14.0, quartiles [10.1, 19.7];
  - observed at s = 0.003: median 13.6, quartiles [10.1, 18.7];
  - PIT/KS p-values 0.27, 0.04, 0.20 at s = 0.003, 0.01, 0.03;
  - p < 10⁻⁴ at s = 0.1 and 0.3.
- The draft never makes this comparison.
- "Relative error O(s)" holds in the bulk (half-width/s ≈ 3). The tails are not O(s): the predicted window W_pred is negative on exact
  draws whose first-order dual is already past the wall.
- "W ∝ 1/s" is overstated at s = 0.1 and 0.3.

**(iv) Stopping and the oracle — oracle and discrepancy rule fine, largest-gap rule overstated.**
- **The oracle holds at first order.** The mean-square error over s²·2K/N is 0.99–1.06. Amplitude and centre biases are within 2
  standard errors at s ≤ 0.1, and +4.5 ± 2.7% at s = 0.3 on one target.
- **The discrepancy rule** is correct per draw (400, 400, 395/396, 374/374). The χ² evidence presented for its late stops is circular
  (D3ar-7).
- **The largest-gap rule is overstated:**
  - it is a post-hoc selection over a horizon chosen knowing K, not a stopping time;
  - it fails as the amplitude range grows and with more features;
  - its later gaps are computed on paths that are often numerically broken (D3ar-3, D3ar-4).
- **Part E's remaining selection:**
  - the one time-out does not matter;
  - Radau paths are not where the rules fail (0/24 and 0/32 at s ≤ 0.03), so the old 80-path selection did not bias the table;
  - what does matter is scoring only exact paths, against the path's own K-th fit rather than the true oracle;
  - and about 5.6% of the exact paths are numerically broken (D3ar-3).

**Novelty — overstated** (D3ar-8).

**Code — `path_response` is correct; `sliding_iss` breaks past the first noise kernel** (D3ar-3).

**Tests — the core recursion is guarded, the support criterion's rows and the classification are not** (D3ar-5).

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3ar-1 | major | **(ii) is linearised along the noiseless order, but the numerics apply it where that order breaks.** Of the 40 targets, 10 have consecutive entries within 5% (target 10: 0.86%; target 20: 1.1%). On exact draws, the entries of E swap order on 23/396 (s = 0.03), 66/374 (s = 0.1) and 87/191 (s = 0.3). Swaps cause 1/1, 7/8 and 42/49 of the "exact, criterion ≥ 1" disagreements. Example, target 20, s = 0.03, draw 7: noiseless order 16, 10, 3 (τ = 22.11, 22.36); noisy order 16, 3, 10. The criterion reads 1.0205 while the path's largest off-E dual is 0.946. Linear-dual error is up to 0.074 / 0.19 / 0.39 on swap draws against ≤ 0.0015 / 0.016 / 0.089 without swaps. | D56 (i)–(ii) and Status: state the radius of (i), in terms of min gap/‖T‖, margins/‖U‖ and λ_min(H_i). Add the order condition (T_{k+1} − T_k)ξ > −(τ_{k+1} − τ_k), or evaluate the criterion along the first-order-reordered skeleton. In `b_job`, record the order and the path's own largest dual (`path_duals`) and report swap draws separately. |
| D3ar-2 | major | **"1 − P(exact) ≤ Σ2Φ̄(…) + O(s²)" and "sharp at first order" are unproved.** The local expansion of (i) is integrated over the whole Gaussian law with no remainder bound. At fixed margins both sides are e^{−c/s²}, so "+O(s²)" is vacuous. Where failures are frequent (margins ≈ s‖U‖), an O(‖ξ‖²) band against an O(s) spread gives an O(s) error, plus P(swap). The linear model is rejected at z = 2.6 and 3.3 (s = 0.1, 0.3). The union bound is never computed in the experiment. Mine gives 4.3, 110 and 400 (capped; 1212 uncapped) expected failures against 4, 26 and 209 observed: valid but loose. | D56 (ii): restate two-sided — exact if crit < 1 − C‖ξ‖², not exact if crit > 1 + C‖ξ‖², for ‖ξ‖ ≤ r₀, with C and r₀ explicit — or as a large-deviation statement. Drop "+O(s²)". Part B: add a union-bound column, and small-margin targets so that failures occur at small s. |
| D3ar-3 | major | **`sliding_iss` past the first noise kernel returns non-fits flagged as relaxed.** An edge kernel escapes the sample window: the noisy fit has no minimiser. With the flow alone, target 16 (s = 0.1, draw 3) reaches c = 1.48, w = −333 and is not relaxed; with the polish, c = 1.70 and w = −3.46·10⁶. The `big`/exit thresholds are relative to max\|w\|, so the feature kernels (\|w\| = 0.39–0.54) count as small. Kernels 2, 4 and 10 exit and re-enter at the same τ (zero gaps), and the fit error is 5.8× the oracle's. Target 14 (s = 0.01, draw 1): c = −1.86, w = 2.65·10⁹, rms residual 411 against the data's 0.33, and `relaxed` is True. In Part E's raw data 19/342 exact paths carry this signature (later fit errors up to 1.2·10⁵ × the oracle's). In my runs a kernel sits outside [−1, 1] after entry K+1 on 8/20, 9/17 and 4/5 paths (wide amplitudes), 8/12 and 8/10 (K = 5) and 8/9 (K+8 horizon); in the K + 8 runs 8 of 9 paths have an unrelaxed block. The bounded polish (from D3ap-8) is defeated by the same relative threshold. The draft does not mention any of this, and `e_job` records neither `relaxed` nor positions. | `rich_transit.sliding_iss` / `polish_block`: bound the centres to the window ± a few σ, or stop with a flag. Use absolute thresholds (exit_tol × the data's scale) for `big`, for exits and for the polish's "seedy" set. Return a `runaway` flag. `e_job`: record `relaxed`, positions and max\|w\|; rerun Part E and report the affected paths. |
| D3ar-4 | major | **The largest relative gap is not a "stopping rule without s" as stated.** (a) It needs the path past the stop, over a horizon (K + 4) chosen knowing K. (b) It depends on the amplitude range, which the Status omits. With \|A\| = (1, 0.3, 0.1) it picks the K-th entry on 18/20, 12/17 and 0/5 draws at s = 0.01, 0.03, 0.1 (W/τ_K medians 15.8, 3.5, 1.4 against the features' own gaps ≈ 2.4); the discrepancy rule manages 17/20, 15/17 and 2/5. (c) With K = 5: 12/12 → 11/12 at s = 0.01 and 10/10 → 7/10 → 6/10 at s = 0.03 as the horizon grows from K + 1 to K + 4; the discrepancy rule gets 12/12 and 10/10. The failing blocks were relaxed and inside the window. (d) "The later gaps … each a relative gap O(1)": in Part E's own raw data the largest later gap reaches 14.4 (s = 0.01) and 9.9 (s = 0.03); with K = 5 a healthy path has gaps 11.3 then 19.5. | D56 (iv) and Status: state the rule as a selection over a fixed horizon, with the condition W/τ_K > max_k(τ_{k+1} − τ_k)/τ_k (≈ A_max/A_min − 1) and its dependence on K. Part E: add wide amplitude ranges, K = 5 and longer horizons. |
| D3ar-5 | major | **The tests do not protect the published numbers.** The criterion's intermediate-entry rows `U_all` are untested: mutants M6 (`U_all` aliased to the final U) and M11 (previous entry's rows) survive all 24 tests. M6 moves Part B's agreement from 388 → 395 and 338 → 367. No test calls `b_job`, `e_job` or `report()`: 7/7 classification mutants survive (exact without the exit condition, crit missing the K-th entry, the window's sign, the gap index, e_job's `exact`, the discrepancy threshold, E included in `off`). `test_window` re-implements the window formula instead of calling the experiment's code. Killed: 8/8 core recursion and Hessian mutants. The cross-curvature mutant is equivalent (Σ r ∂_cφ_i = 10⁻¹⁴ at stationary blocks). | `tests/test_noisy_certificate.py`: finite differences of the duals at every entry against `U_all[k]`; a seeded `b_job`/`e_job` draw with exact, crit, W_pred, k_disc and gaps asserted; `report()` on a synthetic raw file. |
| D3ar-6 | minor | **The limit law of sW (iii) is never tested.** Part C gives only medians and calls W·s = 12.3 and 8.94 "W ∝ 1/s". The PIT test passes for s ≤ 0.03 and rejects at s = 0.1 and 0.3 (p < 10⁻⁴). W_pred < 0 on some exact draws (ratios down to −2.3 and −16.3). | Part C: add the limit law's quantiles/PIT per level. Restrict "∝ 1/s" to s ≤ 0.03. Separate draws with crit ≥ 1. |
| D3ar-7 | minor | **The χ² evidence for the discrepancy rule's late stops is circular.** "17, 11, 18, 9, 18 on all draws" counts rms(R_Kξ) > 1.1s, the first-order statistic itself, which is χ²₅₅ by construction (KS p = 0.91). The path's own late stops are 17, 11, 17, 7, 9 against 18.4, 18.4, 18.2, 17.2, 8.8 expected; at s = 0.1, 7 against 17.2 has p = 0.004. | D56 (D) and the md: compare the path's late counts with the binomial expectation ± sd. |
| D3ar-8 | minor | **Novelty is overstated.** A no-false-positive time ∝ 1/σ and an unbiased oracle at the ISS stopping time are in Osher–Ruan–Xiong–Yao–Yin (ACHA 2016; τ̄ ∝ σ⁻¹√(n/log p), from memory). Early stopping of small-initialisation dynamics at the oracle rate without ℓ₁ shrinkage is Vaškevičius–Kanade–Rebeschini (NeurIPS 2019). Off-grid moving-centre greedy paths are sliding Frank–Wolfe and ADCG (cited under D‴.53 but not here). Knot-based stopping on paths is the covariance test (Lockhart–Taylor–Tibshirani–Tibshirani, AoS 2014); maximising a ratio of consecutive values is the eigenvalue-ratio estimator (Ahn–Horenstein, Econometrica 2013). "Unbiased" holds at first order only. | Novelty paragraph: restrict "new" to the off-grid first-order response and the limit law of sW; add these references; say "unbiased at first order". |
| D3ar-9 | minor | **Part E's remaining selection and scoring.** The rules are scored only on exact paths and against the path's own K-th fit, so 47 paths at s = 0.3 and 10 at 0.1 drop out. "52 of 100" should read 52 of 99. The path over the 180-s limit (target 8, s = 0.3, draw 4) still had not finished after 14 min, so Radau does not cure it. With K + 8 entries, 3 of 12 paths exceed 120 s. The limit is wall-clock and therefore depends on machine load. | `e_job`: score both rules on all finished paths against a known-support oracle (`polish_block` from the truth). Limit by right-hand-side evaluations, not seconds. |
| D3ar-10 | minor | **(D) overstates the gap rule and mislabels the error.** "The largest gap falls after the K-th on every exact draw for s ≤ 0.1" is computed over the first K + 1 entries only, so no noise gap competes. At s = 0.3 the "error/rms(P_Tξ)" is the K-th fit's error even where the rule stopped elsewhere. | D56 (D): say "over the first K + 1 entries" and "the K-th fit's error". |
| D3ar-11 | minor | **"Theorem D‴.53 applies to any data" is overstated.** D‴.53 needs (H1)–(H4). On noisy data its constants degrade with near ties (D3ar-1), and past τ_K + W its hypotheses fail (D3ar-3). | D56: restrict to [0, τ_K + W) and to draws whose path satisfies (H1)–(H4). |
| D3ar-12 | minor | **Index mismatches in the experiment's bookkeeping.** `report()` indexes fits by entry (`errs[argmax(gaps)]`) while `hist` holds one fit per iteration; after a re-entry burst there are more gaps than fits (7 vs 5 on target 16, s = 0.1, draw 3), a latent index error. `e_job`'s `exact` ignores exits after the (K+1)-th entry. | `exp_noisy_certificate.py`: map entries to iterations; flag paths with later exits. |
| D3ar-13 | cosmetic | The recursion's start (τ₀ = 0, δu(0) = 0, R₀ = I) and the noise model in (iv) are not stated. "Two rules stop there" is asserted in (iv), while the Status calls the rules' success numerical. | D56 (i), (iv). |

## Lesson

A linear-response criterion is only as good as its radius. Here that radius is set by the smallest gap between consecutive entries, not
by the noise against the margins. On a quarter of the targets, noise reorders the entries before it ever threatens the support. Per-draw
agreement at the noise levels where the support actually fails therefore mixes in-regime accuracy with out-of-regime luck — so much so
that a deliberately wrong criterion agrees better. Report the regime variable (order kept, distance of the criterion from 1 against Cs²)
next to every agreement count, and compute the bound you state. Numerics past the stop also need health checks: relaxed blocks, centres
inside the window, coefficients on the data's scale. Thresholds relative to max|w| silently accept a runaway kernel as a fit.

Sources: [ORXYY arXiv:1406.7728](https://arxiv.org/abs/1406.7728v5),
[Vaškevičius–Kanade–Rebeschini](https://proceedings.neurips.cc/paper/2019/hash/5cf21ce30208cfffaa832c6e44bb567d-Abstract.html),
[Sliding Frank–Wolfe](https://arxiv.org/pdf/1811.06416), [ADCG](https://ar5iv.arxiv.org/html/1507.01562),
[Lockhart et al.](https://arxiv.org/pdf/1301.7161),
[Ahn–Horenstein](https://www.econometricsociety.org/publications/econometrica/2013/05/01/eigenvalue-ratio-test-number-factors).
