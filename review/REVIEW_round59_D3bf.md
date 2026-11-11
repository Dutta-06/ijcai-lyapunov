# Referee report D3bf, round 59: Proposition D‴.52 as revised after round 55

**Scope.** I reviewed these files at 03d01bb:
- `review/drafts/D52.md`
- `experiments/exp_close_features.py`
- `results/chaos/close_features.md` and `close_features_raw.npz`
- `tests/test_close_features.py`
- `src/core/rich_transit.py` (`sliding_iss`, `polish_block`)
- `review/REVIEW_round55_D3bb.md` and the round-55 entry in `review/REVIEW_LOG.md`

None of these files changed between 03d01bb and the current HEAD a1f5d85 (`git diff --stat` is empty). I modified nothing in the repository. All work ran on a `git archive` copy in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref59/D3bf_work/`. Every run used `nice -n 10`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`, with at most two heavy processes at a time.

## Summary

**What checks out (verified):**
- **Reproducibility is exact.**
  - `--from-raw` into my own `CF_OUT` reproduces the committed md byte for byte (`cmp`).
  - The test file passes: 16 passed in 21.4 s.
  - Fresh reruns are bit-identical to the raw:
    - A: three rows.
    - B: three rows, all fields.
    - C/D: the first 41–60 draws of three cells, including 80 placements.
    - E: 50 H0 and 25 H1 draws in two cells, checking Δτ₂, the entrants, the first fit (ĉ, ŵ) and the per-draw bootstrap p, plus `q_law`, `np` and the H0 truth.
- **(i) and (ii) re-derived:**
  - V, T₃ and M₄ = (A₁A₂Δ⁴/A)(1 − 3A₁A₂/A²).
  - The fit's shift −T₃/(4Aσ²).
  - ⟨ψφ(c̄+d)⟩ = G(d)d²/(4σ⁴).
  - det ∂(M₂,M₃)/∂(δ,θ) = −W²θ(1−θ)δ⁴ = −w₁w₂δ⁴.
  - The split and separation images, the three continuum constants, r_split, and μ = 4Gc|w|/a.
  - Every number of A and B in the numerics matches the md, except "factor 18".
- **(i) holds on pairs it was not tested on:**
  - (0.8, 0.2): with T₃ the ratio is 1.0048, 1.0013, 1.0003.
  - (0.8, −0.6), where ϑ = 58Δ: with T₃ the ratio is 1.311, 1.0715, 1.0166. The remainder is ≈ 0.125ϑ², the same constant as (0.8, −0.4)'s.
- **The new null law is right.**
  - Under H0, `knot_of` is the path's second knot to first order. The relative error scales ∝ s at c̄ = 0, 0.032 and 0.049; at s = 10⁻³ its median is 0.9–1.1·10⁻⁵ and its maximum 2·10⁻³.
  - The law evaluated at the *true* fit calibrates the knot in all 15 E cells: sizes 0.039–0.063, pooled 0.0488 over 15 000 draws.
  - The fixed threshold's false alarms (0.21–0.45) and the 17× swing of q_α·s are confirmed.
- **(iv): placement at the Cramér–Rao bound is confirmed** in the resolved cells.

**What is wrong (majors first).**

1. **"Near-Neyman–Pearson power" is contradicted by Part E, and the matched-filter explanation does not account for the shortfall** (D3bf-1).
   - The bound exceeds the knot's power by 0.085–0.157 at c̄ = 0 and 0.02, and by 0.31 at (0.05, 0.03, 0.049).
   - The runner-up's filter is *best* aligned (cos 0.989) exactly where the power is worst.
   - A one-line score test computed from the same first fit has an analytic N(0,1) null. On E's own draws it sits at the bound in all 15 cells and beats the calibrated knot by 0.022–0.332 in every cell.
2. **For an opposite-sign pair the knot test is biased** (D3bf-2).
   - At (0.8, −0.4) its power is 0.010 and 0.005, below α, where the Neyman–Pearson bound is 0.47 and 0.91.
   - Its Δ₅₀ is 0.040, against 0.018 for the Neyman–Pearson bound.
   - (iii) is stated for general pairs, but only (0.8, 0.6) was tested.
3. **Novelty** (D3bf-3).
   - On a fixed dictionary the transit's second knot *is* the LARS second knot: τ₂ = 1/λ₂ to 5·10⁻¹⁵. So Δτ₂ is a spacing statistic, and spacing statistics have exact conditional nulls.
   - "Its calibration by that law at the first fit" is the standard parametric bootstrap at the fitted null.
   - The decision-theoretic resolution literature and the RBF overlap-singularity literature are missing.

**Minors:**
- The resolution limit depends on the offset: Δ₅₀ ≈ 0.076 at c̄ = 0.049 against 0.051 at c̄ = 0.032, for s = 0.03. The 1/s law needs a no-tie condition at the noise scale (D3bf-4).
- The bootstrap is valid but conservative near ties: 0.019 and 0.022. The cause is the plug-in (D3bf-5).
- D's position jump follows s/Δ² ≳ 12, not "s/Δ³ ≳ 240" (D3bf-6).
- In some cells the split's error lies below the Cramér–Rao bound (D3bf-7).
- C's 200-draw thresholds are liberal: realized size 0.080, and the power exceeds the Neyman–Pearson bound in 5 cells (D3bf-8).
- There are new misquotes and a floating-point count, "12 of 15" (D3bf-9).
- E's data-generating code is untested: 17 of 36 mutants survive, and 6 of them leave both tests and md unchanged (D3bf-10).

## What I checked and how

Scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref59/D3bf_work/scripts/`; `common.py` loads the copy and the raw.

- **Reproducibility.**
  - `CF_OUT=<work>/out_fromraw python experiments/exp_close_features.py --from-raw`, then `cmp` against the committed md.
  - `python -m pytest -q -p no:cacheprovider tests/test_close_features.py`.
  - Reruns with the experiment's own seeds:
    - `rerun_AB.py`: A rows (0.8,0.6,0.0125), (0.8,−0.8,0.025), (0.8,0.8,0.05); B rows (0.8,0.6,0.0125), (0.8,−0.4,0.2), (0.8,0.8,0.025).
    - `rerun_CD.py`: C/D (0.3,0.003,H1), (0.025,0.01,H0), (0.05,0.03,H1).
    - `rerun_E.py`: E (0.025,0.01,0.049) and (0.05,0.03,0.045), consuming the main rng exactly as `e_job` does.
- **Part E from the raw.**
  - `e_check.py`: Clopper–Pearson intervals of the sizes, the law-at-truth size, the exact power differences.
  - `e_summary.py`: the summary ranges, the sampling error of the draws' quantile, the fixed threshold without its in-sample cell.
  - `c_vs_e.py`: C's thresholds applied to E's H0 draws.
- **New analyses on E's draws** (regenerated from E's seeds):
  - `score_test.py`: z = ⟨r̂₁, ψ̂⟩/(s‖ψ̂‖), where ψ̂ is the second-derivative shape projected off span(φ, ŵ∂φ) at the path's own first fit; plus a t-version with ŝ.
  - `boot_variants.py`: the per-draw bootstrap three ways on the same bootstrap normals — at (ĉ, ŵ), which reproduces the raw p (asserted); at the true fit; and with ŝ.
  - `mf_cos.py`: the cosines between Rφ_j and f₁ − f₀.
  - `margin_noise.py`: the runner-up's margin and its noise.
  - `eff_crit.py`: the power-½ criterion's effective critical value.
- **Path-level checks.**
  - `knot_first_order.py`: `knot_of` against `sliding_iss` at 3 offsets × 4 values of s × 12 noise vectors.
  - `lars_iss.py`: ISS against LARS knots on a fixed dictionary.
  - `opp_expansion.py`: (i) on (0.8, 0.2) and (0.8, −0.6).
  - `detect.py`: new detection runs, 200 draws per hypothesis per Δ — (0.8, −0.4) at c̄ = 0.032, s = 0.01; and (0.8, 0.6) at c̄ = 0.049, s = 0.03.
  - `d_check.py`: D's errors against s/Δ², s/Δ³ and the Cramér–Rao bounds.
- **Mutation.**
  - 36 single-point mutants (`mutants.py`, `run_mut.py`, results in `mut/results_{1,2}.json`).
  - Each was run against the test file with `-x` and checked for a change in the md regenerated `--from-raw`.
- **Literature.**
  - WebFetch could not resolve ar5iv.labs.arxiv.org, hal-univ-tlse3.archives-ouvertes.fr or ideas.repec.org. The Semantic Scholar API returned 403 through the proxy.
  - All literature checks therefore rest on search-engine snippets, marked [snippet]; nothing was read in full.

## Verdicts per part

| claim | verdict |
|---|---|
| (i) expansion, moments, ψ-identity, ĉ − c̄ = −T₃/(4Aσ²) | **holds** (re-derived) |
| (i) leading order, error term O(ϑ) + O((Δ/σ)²(1+\|A₁A₂\|/A²)), T₃ remainder O(ϑ²), no-tie hypothesis | **holds**; also on (0.8, 0.2) and (0.8, −0.6), remainder ≈ 0.12–0.13ϑ² |
| (i) opposite-sign wall (margin 2 − μ); dipole ∝ Δ⁻¹ | **holds** |
| (ii) λ_split, λ_sep, det, confluent-Vandermonde orders | **holds** (orders classical) |
| (ii) flow −DH, r_split, the separation's two regimes | **holds**; "within a factor 18" is misquoted (D3bf-9) |
| (iii) projected rates; `knot_of`/`knot_law` first-order exact under H0 | **holds** (error ∝ s; law at the truth gives size 0.039–0.063) |
| (iii) null depends on the offset; a fixed threshold is miscalibrated | **holds** |
| (iii) "a delay ∝ √N/s" | **holds away from ties**; near a tie at the noise scale q_α·s falls 0.813 → 0.463 (D3bf-4) |
| (iii) power ½ where ½\|V\|⟨ψφ_{j₂}⟩ ≈ z₁₋α·sd | **approximate**: the effective critical value is 1.64–2.49 |
| (iii) runner-up ≈ matched filter ⇒ power close to the Neyman–Pearson bound | **fails** in 8 of 15 E cells, and for opposite-sign pairs the power is below α; the explanation is incomplete (D3bf-1, -2) |
| (iii) bootstrap "keeps its size … at nearly the oracle's power" | **valid but conservative near ties** (0.019, 0.022); 13 of 15 within 0.05 (D3bf-5) |
| (iv) orders s/Δ², s/Δ³; Cramér–Rao bound attained | **holds** while resolved (0.84–1.18); the split falls below the bound at the edge (D3bf-7) |
| (iv)/D: split lost at ≈ 15Δ³; positions jump at s/Δ³ ≳ 240 | 15Δ³ is a **chosen** threshold; the jump is at **s ≈ 12Δ²**, the detection limit (D3bf-6) |
| Status "New" | **partly**: the delay law with its offset dependence and the mobility-weighted rates are defensible; the statistic is a LARS spacing statistic, the calibration is standard, and near-Neyman–Pearson power is false (D3bf-3) |

**Status of the round-55 issues**

| D3bb | status | evidence |
|---|---|---|
| 1 calibration across the offset | **partly** | E added; the calibrated test is valid but conservative near ties (D3bf-5); "near-NP" not withdrawn (D3bf-1); Δ₅₀ still only at c̄ = 0.032 (D3bf-4) |
| 2 noise law, H0 mechanism | **fixed** | projected rates; `knot_of` first-order exact; 16.6% / 70.2% printed |
| 3 literature, novelty | **partly** | references added; "New" still claims the calibration and near-NP power; τ₂ = 1/λ₂; classical one-vs-two and RBF literature missing (D3bf-3) |
| 4 error term | **fixed** | T₃ column; holds on two further pairs |
| 5 no-tie hypothesis | **fixed** for the deterministic tie; the noise-scale analogue is missing (D3bf-4) |
| 6 r_sep regimes | **fixed** | `sep_regimes`; new misquote "factor 18" |
| 7 D saturation | **partly** | Cramér–Rao columns, constants and counts added; the jump's scaling is still wrong and "4–22" should read 5–22 (D3bf-6) |
| 8 misquotes | **fixed** | all eight verified; new ones in D3bf-9 |
| 9 weak tests | **partly** | analogues of the round-55 survivors are now killed (my M26, M29, M30, M36); the new E code is untested (D3bf-10) |
| 10 statistical precision | **partly** | C has bootstrap intervals and E uses 1000 H0 draws; C's q_α is still from 200 draws (D3bf-8); E prints no intervals |
| 11 cosmetics | **fixed** | the drafts README still says "Parts C and D queued" |

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3bf-1 | major | **The near-Neyman–Pearson claim fails across the grid cell, and the matched-filter explanation does not explain it** (D52.md:53–55, 107). E's own data: NP − max(oracle, calibrated) is 0.157 / 0.085 / 0.095 at c̄ = 0 and 0.122 / 0.120 / 0.090 at c̄ = 0.02 (all three cells), and 0.095 and 0.310 at c̄ = 0.049 for s = 0.01, 0.03. Only 7 of 15 cells are within 0.05. The filter argument predicts the opposite pattern: cos(Rφ_{j₂}, f₁ − f₀) runs from 0.94 at c̄ = 0 to 0.99 at 0.049, best where the power is worst. A one-sided test on n_{j₂} alone would give 0.65–0.69, 0.34–0.37 and 0.50–0.54 (NP 0.70, 0.375, 0.54). The losses come from what the knot adds to the filter. (a) The minimum over the dictionary: at c̄ = 0 and 0.02 the runner-up's margin/projected-noise ratio (0.73, 0.68) is close to the far kernels' ≈ 1, so they share the H0 tail; the effective critical value is 1.78–1.94 instead of 1.645. (b) The margin's own noise: sd ≈ 0.16s against μ_{j₂} = 0.0035 at c̄ = 0.049, giving cv 0.46 and 1.07 at s = 0.01, 0.03; the first entrant swaps in 23% of draws, and the effective critical value is 2.49. Benchmark on E's own draws: the score test z = ⟨r̂₁, ψ̂⟩/(s‖ψ̂‖) at the same first fit has size 0.036–0.064 (pooled 0.0493) and power within −0.057 … +0.052 of the bound in all 15 cells. It beats the calibrated knot by +0.022 … +0.332 in all 15 cells. Its t-version with ŝ has size 0.039–0.063. | (iii) l. 53–57 and New l. 107: drop "near-NP" or restrict it to c̄ ≈ 0.032–0.045, and replace the matched-filter sentence by the two loss mechanisms. Add the score test to `e_job` and the E table as the benchmark. |
| D3bf-2 | major | **Opposite-sign pairs: the knot test is biased.** (iii) states the power with \|V\| (l. 51–52), and the title says "counts two features", but C and E test only (0.8, 0.6). For (0.8, −0.4) at c̄ = 0.032, s = 0.01 (`detect.py`, 200 + 200 draws, oracle threshold), the knot's power is 0.010 [0.001, 0.036] at Δ = 0.018 and 0.005 [0.000, 0.028] at Δ = 0.025. That is below α, where the Neyman–Pearson bound is 0.471 and 0.908 and the one-sided score test 0.455 and 0.895 (size 0.045–0.06). At Δ = 0.035 it is 0.24 against NP 1.000, and 0.96 at Δ = 0.05. The knot's Δ₅₀ is ≈ 0.040 against the bound's 0.018: 2.2× in Δ, ≈ 5× in SNR. The mechanism is the drift: the pair pushes the A-sign neighbours' duals away from their wall. That removes H0's small-margin early-knot channel. The flank kernel then has to cover a margin ≈ 1 + u. | (iii): restrict the detection and power statements to same-sign pairs, and state the bias for A₁A₂ < 0. C/E: add an opposite-sign pair. For both signs, use a two-sided statistic (\|z\|, or each wall's knot against its own null). |
| D3bf-3 | major | **Novelty and literature** (l. 101–108). (a) On a fixed dictionary, after the first entry, u_j(τ) = c_j(1/τ)/(1/τ). So τ₂ = 1/λ₂, the LARS/lasso second knot: max\|τ₂λ₂ − 1\| = 4.9·10⁻¹⁵ over 200 random noisy pairs (`lars_iss.py`). Hence Δτ₂ = (λ₁ − λ₂)/(λ₁λ₂) is the spacing statistic, which has exact conditional nulls: Tibshirani–Taylor–Lockhart–Tibshirani, JASA 2016 (the spacing test is theirs, not Taylor–Loftus–Tibshirani's as at l. 103); its power, Azaïs–De Castro–Mourareau, Bernoulli 2018; grid-less, ADM 1706.00679. Sliding only re-fits the first atom. "Δτ₂ its analogue … [memory]" (l. 104) understates this. "Its first-order null law … its calibration by that law at the first fit" is a delta-method law plus the parametric bootstrap at the fitted null (McLachlan 1987, number of mixture components). (b) Testing one source against two with a likelihood-ratio or Neyman–Pearson benchmark is classical: Harris 1964; Helstrom 1973; Liu–Nehorai 2007 (GLRT resolution limit ∝ √CRB, the link credited to Smith at l. 95–96); Kulaitis–Munk–Werner 2021 (minimax testing of one vs two sources). (c) (ii)'s slow modes are the overlap singularity of Gaussian RBF networks under gradient learning (Amari–Park–Ozeki 2006; Wei–Amari 2008 [snippet]). EM slows the same way at it (Dwivedi et al. 2020), and (iv)'s orders have a mixture analogue (Heinrich–Kahn 2018). ISS under noise is not cited (Osher–Ruan–Xiong–Yao–Yin 2016). | Status: state τ₂ = 1/λ₂ on a fixed grid and cite the spacing tests as the statistic's origin. Limit "New" to the delay law with its offset dependence and its opposite-sign and dipole branches, and to the mobility-weighted rates, with the RBF references. Drop "calibration" and "near-NP" from New. |
| D3bf-4 | minor | **The resolution limit depends on the offset, and E does not measure that.** E has one Δ per s, and Δ₅₀ exists only at c̄ = 0.032 (C). At c̄ = 0.049, s = 0.03 (`detect.py`, 200 + 200 draws), the knot's power is 0.175, 0.405, 0.315, 0.73, 0.805 at Δ = 0.05, 0.06, 0.07, 0.085, 0.1. NP is 0.543, 0.807, 0.960, 0.9995, 1.000 and the score test 0.47, 0.79, 0.96, 1.00, 1.00. So Δ₅₀ ≈ 0.076, against 0.051 at c̄ = 0.032 and the bound's 0.048, and the power is 0.315 at a separation where the bound is 0.96. (iii)'s "∝ √N/s" (l. 48) and its power-½ criterion (l. 51) also need a no-tie condition at the noise scale, μ_{j₂} ≫ sd(μ̂) ≈ 0.16s here. At c̄ = 0.049, q_α·s = 0.813, 0.594, 0.463 for s = 0.003, 0.01, 0.03. | (iii): add the noise-scale condition. E: give the calibrated test's Δ₅₀ per offset over both halves of the cell, and say that the quoted limit is the favourable offset's. |
| D3bf-5 | minor | **The bootstrap is valid but conservative near the tie**, and the cause is the plug-in. Sizes: 0.019 [0.011, 0.030] at (0.025, 0.01, 0.049) and 0.022 [0.014, 0.033] at (0.05, 0.03, 0.045). It is never significantly liberal (maximum 0.062 [0.048, 0.079]), but the pooled size is 0.0457 [0.0424, 0.0492]. With the law at the true fit, on the same draws and the same bootstrap normals, these two cells give 0.045 and 0.039; across all cells 0.039–0.063, pooled 0.0488. Under H0, corr(ĉ − c, log Δτ₂) is −0.14, −0.17 and −0.20 in the three near-tie cells and \|·\| ≤ 0.06 elsewhere. The noise that moves ĉ toward the tie also shortens the knot, and the law at ĉ explains it away. The power loss, −0.09, falls in exactly these two cells. The bootstrap also uses the true s, which is not stated. With ŝ² = ‖r̂₁‖²/(N−2) it still holds: sizes 0.019–0.060, pooled 0.0435, power change −0.04 … 0. | (iii) l. 55–57: write "valid; conservative near a tie; s known (estimating it works)". E: print Clopper–Pearson intervals. Consider resampling only the rates n given the observed margins, since the duals at τ₁ are observed. |
| D3bf-6 | minor | **D's saturation is still mis-described** (l. 143–144; D3bb-7). "(s/Δ³ ≳ 240)" is the wrong variable. (0.0125, 0.001) at s/Δ³ = 512 and (0.0175, 0.003) at 560 have position rms 0.065 and 0.152, while (0.05, 0.03) at 240 has 0.293. The 9 jumped cells (0.29–0.41) are exactly the cells with s/Δ² ≥ 12; the other 31 have s/Δ² ≤ 9.8 and rms ≤ 0.152. This is the detection limit (C: s/Δ₅₀² = 12.8, 13.3, 11.5), and "both within σ" tracks C's power (r = 0.99). "4–22 of 40" should read 5–22. | (iv) l. 62–64 and numerics l. 143–144: the positions jump at s ≈ 12Δ², where the pair stops being detected, which is (iii)'s order. |
| D3bf-7 | minor | **"At the Cramér–Rao bound … 0.67–1.05 (split)" includes biased cells** (l. 141–142). (0.1, 0.03) is at 0.668 of the bound and (0.07, 0.01) at 0.793, about 3σ below with 40 draws. The split is bounded, so the estimator is no longer unbiased there. The split follows the bound (0.93–1.05) up to s/Δ³ ≈ 24 and drops below it from 29. "Lost at s ≈ 15Δ³ (the bound's std 0.1)" (l. 62) is a chosen threshold, not a measured transition. | Write "at the bound for s ≲ 25Δ³, below it (biased) beyond". |
| D3bf-8 | minor | **C's thresholds still come from 200 H0 draws.** C's q_α at (0.025, 0.01) has a realized size of 0.080 on E's 1000 independent H0 draws of the same cell. Hence C's power 0.42 exceeds the bound 0.375; E's oracle power in that cell is 0.325. Five C cells exceed the bound: 0.42 vs 0.375, 0.125 vs 0.077, 0.12 vs 0.062, 0.30 vs 0.295, 0.975 vs 0.96. Part of "Δ₅₀ 0.0274 … against the Neyman–Pearson bound's 0.0274" (l. 128–130) is this threshold noise. | Rerun C's H0 with 10³ draws, or calibrate it by the law at the truth, which E shows to be exact to about ±0.015. Then recompute Δ₅₀. |
| D3bf-9 | minor | **Misquotes and artefacts.** (1) l. 123, "within a factor 18": the md's amplitude-mode column (min μ_k⟨φ_k²⟩) gives 26.8, 27.8 and 22.7; 18 is (0.8, 0.6)'s second flow rate, and the three pairs give 12–21. (2) l. 139 and md:159, "12 of 15": the +0.050 cell is 0.34 − 0.29 = 0.05000000000000004 > 0.05 in floating point, so with the printed range it is 13. (3) l. 137, "0.027–0.050 at c̄ ≤ 0.032": the 0.050 is the reference cell, in-sample by construction; out of sample the range is 0.027–0.044. (4) l. 136, "0–22% … on the fit's two neighbours": one of the two is j₁ itself (0.000–0.057); the runner-up candidates are 6–22%. (5) l. 146, "the split's error 79": this is the cell's rms; the draw's own error is ≈ 500. | Correct the text, and compare with a tolerance in `report` (exp_close_features.py:656). |
| D3bf-10 | minor | **Tests.** Of 36 single-point mutants, 16 are killed and 3 are equivalent (M15 through float, M22 at the printed precision, M35 because `hist` has two entries). 17 non-equivalent mutants survive. **E's data path is untested**: in `e_job`/`part_E` (l. 298–337) the following pass the tests *and* leave the md from the raw unchanged — the p-value's tail flipped (M07); the law at the truth instead of the draw's fit, which silently turns the calibrated test into the oracle-law test (M08); H0 drawn from the pair (M09); the ratio without the square root (M10); one seed for all cells (M11); the law at s/2 (M33). The expected strings in `test_the_report_on_synthetic_results` (tests:245–249) are computed with `CF.np_power` and `CF.crb` themselves, so M12 (np_power's quantile) and M21 (s instead of s² in crb) survive. Also surviving: `d50_bootstrap` and `crossing` (M18–M20), report D (M23–M25), report C (M31) and two E columns (M16, M17). The one `knot_of` test (tests:151–167) rests on a single noise vector whose error is 10⁻⁸; over 36 vectors the median is 10⁻⁵ and the worst 2·10⁻³ at s = 10⁻³. The docstring and the round-55 log quote that 10⁻⁸. | Add a small `e_job` test with R0_E, R1_E and M_E monkeypatched, checking p against a direct `knot_law`, H0 ≠ H1, and distinct seeds. Use analytic expectations for `crb` (e.g. a diagonal J) and `np_power`. Test `d50_bootstrap` on a synthetic table. Test `knot_of` on several noise vectors with an error bound ∝ s. |
| D3bf-11 | cosmetic | The drafts README still says "Parts C and D queued". E's "std ratio" column does not name its kernels. The law/draws quantile ratio 0.891–1.080 lies inside the draws' own sampling error in 15 of 15 cells (bootstrap half-widths 3–21%); the realized size at the law's quantile (0.039–0.066, pooled 0.049) is the sharper check. | README; md caption; E summary. |

Sources (all from search-engine result snippets; no page was fetched in full):
- [Azaïs–De Castro–Mourareau, Testing Gaussian process with applications to super-resolution (arXiv 1706.00679)](https://arxiv.org/pdf/1706.00679)
- [Azaïs–De Castro–Mourareau, Power of the spacing test for LARS (Bernoulli 2018; arXiv 1503.05093)](https://arxiv.org/pdf/1503.05093)
- [Tibshirani–Taylor–Lockhart–Tibshirani, Exact post-selection inference for sequential regression procedures (JASA 2016; arXiv 1401.3889)](https://arxiv.org/pdf/1401.3889)
- [McLachlan 1987, bootstrapping the LRT for the number of mixture components (JRSS C 36:318)](https://academic.oup.com/jrsssc/article/36/3/318/6985676)
- [Kulaitis–Munk–Werner, What is resolution? (AoS 49(4), 2021)](https://projecteuclid.org/journals/annals-of-statistics/volume-49/issue-4/What-is-resolution-A-statistical-minimax-testing-perspective-on-superresolution/10.1214/20-AOS2037.full)
- [Liu–Nehorai, Statistical angular resolution limit for point sources (IEEE TSP 55(11), 2007)](https://profiles.wustl.edu/en/publications/statistical-angular-resolution-limit-for-point-sources/)
- [El Korso et al., GLRT-based resolvability of two point sources](https://hal-supelec.archives-ouvertes.fr/hal-00771379)
- [Harris 1964, Resolving power and decision theory (JOSA 54(5))](https://opg.optica.org/JOSA/issue.cfm?issue=5&volume=54)
- [On Helstrom 1973, one vs two sources (arXiv 1802.02300)](https://ar5iv.arxiv.org/html/1802.02300)
- [Tsang–Nair–Lu, Rayleigh's curse (PRX 2016; arXiv 1511.00552)](https://arxiv.org/abs/1511.00552v2)
- [Amari–Park–Ozeki 2006, Singularities affect dynamics of learning](https://research.knu.ac.kr/en/publications/singularities-affect-dynamics-of-learning-in-neuromanifolds/)
- [JMLR 2018/2020, citing Wei et al. 2008 and Wei–Amari (RBF overlap singularities)](https://jmlr2020.csail.mit.edu/papers/volume19/16-210/16-210.pdf)
- [Dwivedi et al., Singularity, misspecification and the convergence rate of EM (AoS 2020)](https://arxiv.org/pdf/1810.00828)
- [Heinrich–Kahn, Strong identifiability and optimal minimax rates for finite mixture estimation (AoS 2018)](https://projecteuclid.org/journals/annals-of-statistics/volume-46/issue-6A/Strong-identifiability-and-optimal-minimax-rates-for-finite-mixture-estimation/10.1214/17-AOS1641.pdf)
- [Osher–Ruan–Xiong–Yao–Yin, Sparse recovery via differential inclusions (arXiv 1406.7728)](https://arxiv.org/pdf/1406.7728)
- [Denoyelle–Duval–Peyré–Soubies, the sliding Frank–Wolfe algorithm (arXiv 1811.06416)](https://arxiv.org/pdf/1811.06416)

## Lesson

A benchmark beats an explanation. The round-55 story — the runner-up's projected kernel is nearly a matched filter, hence near-Neyman–Pearson power — was checked at one offset and one sign. Computing the matched filter itself on the same first fit takes one line, ⟨r̂₁, R̂e₂(ĉ)⟩/(s‖R̂e₂‖). It gives an analytically calibrated test at the Neyman–Pearson bound everywhere. It also shows what the knot adds to that filter: a noisy grid margin in the numerator and a minimum over the dictionary. Those cost 0.02–0.33 in power across one grid cell and make the statistic blind, with power below α, to opposite-sign pairs. When a fix adds a sweep over one nuisance parameter, re-test every earlier claim on that sweep, and treat the sign of the alternative as a nuisance too.
