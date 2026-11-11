# Referee report D3bk, round 64: Proposition D‴.52 as revised after round 59

## Scope

I reviewed these files at 8c826c6 (branch claude/festive-bardeen-eopej2):
- `review/drafts/D52.md` (196 lines)
- `experiments/exp_close_features.py` (parts A, B, C, C′, D, E, E′/S, F, and `--from-raw`)
- `results/chaos/close_features.md` and `close_features_raw.npz`
- `tests/test_close_features.py` (24 tests)
- `review/REVIEW_round59_D3bf.md` and the "Round 59" entry of `review/REVIEW_LOG.md`
- `src/core/rich_transit.py` (`sliding_iss`, `polish_block`) and `src/core/mirror_free.py` (`iss_path`)
- `review/drafts/D54.md` and `D59.md`, for the two readings in (iv)

While I worked, HEAD moved to 9dc2d24. Of the reviewed files only `rich_transit.py` changed: it gained a `gtol`/`polish_gtol` option whose default (10⁻¹⁵) keeps the reviewed behaviour. The working tree also carries uncommitted edits by others; I used none of them.

I modified nothing in the repository. All work ran on `git archive 8c826c6` in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref64/D3bk_work/repo`. I used the repo's venv with `nice -n 10`, `OMP_NUM_THREADS=1` and `PYTHONDONTWRITEBYTECODE=1`. At most two heavy processes of mine ran at a time, and none is left running.

## Summary

**What checks out**

- **Reproducibility is exact.**
  - `--from-raw` into a scratch `CF_OUT` reproduces the md byte for byte (`cmp`).
  - The 24 tests pass in 20.4 s.
  - Reruns with the experiment's own job functions and seeds are bit-identical to the raw:
    - all of S: 95 jobs and 37 000 draws, every row field, the laws' quantiles, j₂, μ and cos;
    - 18 F jobs: 2 H0 and 16 H1, 4 000 draws, every field including the bootstrap p-values;
    - B rows (0.8, 0.8, 0.0125) and (0.8, −0.4, 0.05);
    - A row (0.8, 0.6, 0.05).
  - The raw A, B, C and E parts are identical to round 59's, which D3bf reran. The A–E job functions are unchanged (checked by an `ast` diff).
  - One exception at round-off is listed under D3bk-11.
- **Numbers.** Every number in numerics (A)–(F) matches the md or raw, except two misquotes (D3bk-11).
- **Re-derived, then checked numerically:**
  - the second-knot identity u_j(τ) = τc_j(1/τ) (proof l. 99–101);
  - z ~ N(0, 1) under H0 to first order. Over all 42 200 H0 draws of C, E and F: mean −0.005, sd 0.995, KS p = 0.38;
  - z's H1 mean ‖f₁ − f₀‖cos/s. Here 1 − cos = ϑ²/12 + O(ϑ⁴), which sharpens "1 − O(ϑ²)": continuum 0.0832ϑ², measured 0.083ϑ² for four pairs;
  - reading (a): ∂_σφ = σ∂²_cφ + φ/σ, and the Gauss–Newton step σ̂ − σ = sz/(ŵσ‖R̂∂²_cφ‖). σ̂² − σ² → V/A holds: at Δ = 0.005 the ratio is 1.0002 for (±0.8, ±0.6), 0.994 for (0.8, −0.4) and 0.969 for (0.8, −0.6);
  - reading (b): the relative mode's rate at a merged state. I checked it by finite differences of ċ_k = Gmw_k⟨r∂_cφ(c_k)⟩, and it holds to 5·10⁻⁶ for both signs and both kinds of split.
- **The one-sided score sits at the Neyman–Pearson bound.**
  - F gives Δ₅₀/NP of 1.003 [0.952, 1.048] and 1.011 [0.965, 1.058].
  - At pairs the authors did not run: 1.036 and 1.076 for (0.8, −0.6) at s = 0.01 and 0.03; 1.045 for (0.8, 0.2).
- **The knot's opposite-sign bias is real, and worse at (0.8, −0.6):** power ≤ 0.015 at every Δ ≤ 0.028 (s = 0.01), where NP reaches 1.000.
- **The loss mechanisms in (iii) are right.** I confirmed them by a decomposition the repo does not make (D3bk-8).
- **Other checks that hold:**
  - the √N/s law: at N = 121 and 241 the null quantile grows by 1.403 and 1.990, against √(N/61) = 1.408 and 1.988;
  - D's placed blocks are at the least-squares optimum. Re-polishing without the gradient test moves the split by ≤ 5·10⁻⁷ of its Cramér–Rao bound (CRB);
  - (v) and the fixes of D3bf-4, -6, -7, -8 and -11.

**What is wrong (majors first)**

1. **D3bk-1 (major): the knot is a de facto one-sided test, and the draft's comparisons with the score ignore that.**
   - Under H0, 92% of the knot's rejections fall on the wall of the fit's sign (422/461 and 396/431; 75–100% per offset). The grid's margins put α there, not a choice.
   - This one allocation explains two things the draft treats as separate:
     - for same-sign pairs, the calibrated knot beats the sign-agnostic |z| in 17 of 24 F cells and in 10 of 15 E cells;
     - for opposite-sign pairs, the drift empties the wall where α was spent, which is the bias.
   - The score reaches the bound only when it is told sign(A₁A₂). Yet l. 60–61 compares the two-sided walls test with the one-sided z, and the title claims the NP bound without the sign condition.
2. **D3bk-2 (major): "testing each wall … at α/2 removes most of the bias" is false.**
   - Pooled over the 24 opposite-sign cells, its power at the two smallest separations is 0.036 [0.031, 0.042] and 0.039 [0.034, 0.045]. That is below α, where NP is 0.19–0.59.
   - Its least power is below α in 24 of 24 cells.

**Minor issues**

- D3bk-3: the Smith/Cramér–Rao link is mis-derived. With unknown, unequal amplitudes the separation's bound is ∝ s/Δ², so Smith's criterion sits at s ∝ Δ³.
- D3bk-4: "τ_k = 1/λ_k" holds for the first two knots only.
- D3bk-5: the mean of z has V's sign, not A₁A₂'s, and the notation clashes.
- D3bk-6: the non-monotone power comes from the null moving across the tie.
- D3bk-7: the calibrated knot is conservative beyond the tie.
- D3bk-8: "lower quantiles by the neighbours" depends on the offset, and the loss decomposition is not measured.
- D3bk-9: 33 of 69 mutants survive both the tests and `--from-raw`, 25 of them non-equivalent. They include all of S's data path and D3bf's M08.
- D3bk-10: reading (b) is the classical splitting criterion.

**Cosmetic:** D3bk-11.

## What I checked and how

Scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/ref64/D3bk_work/scripts/` (`common.py` loads the copy and the raw), and the mutation files in `…/D3bk_work/mut/`.

- **Reproducibility**
  - `cmp` of the md regenerated `--from-raw`, and `python -m pytest -q -p no:cacheprovider tests/test_close_features.py`.
  - `rerun_S.py`: `CF.part_S(C, E)`, all 95 jobs, compared recursively with the raw.
  - `rerun_F.py`: 18 `CF.f_job` jobs, comprising:
    - (0.8, −0.4) at c̄ = 0.032, s = 0.01, with H0 and all six separations;
    - (0.8, 0.6) next to the tie;
    - 8 random jobs.
  - `rerun_AB.py` and `rerun_A_diff.py`: two A rows and two B rows. `cmp_old.py` compares the raw with 03d01bb's.
- **Mathematics**
  - `lars_iss.py`: `iss_path` on the fixed dictionary against my LARS and lasso homotopies, over 300 noisy targets (one bump, same-sign pairs and opposite-sign pairs).
  - `crb_sep.py`: the CRB of c₂ − c₁, of the positions and of the split, for four pairs, with the amplitudes unknown and known.
  - `score_sign.py`, `width_check.py`, `cos_theta.py`, `theta12.py`, `relmode_check.py`, `sqrtN.py`, `s02_equiv.py`.
- **From the raw**
  - `f_analysis.py`: per-wall one-sided tests (`pA`, `pB` are in the raw) and the knot against |z|.
  - `f_detail.py`: the null fits ĉ, z's null law per cell, and q_law, H1 medians and entrants per separation.
  - `f_walls_bias.py`, `walls_split_h1.py`, `f_wall_split.py`, `f_alpha01.py`, `h0_tails.py`, `p_ties.py`, `cal_low.py`.
  - `e_vs_absz.py` and `e_median_adj.py`.
  - `e_decomp.py`: regenerates E's 21 000 labels from E's seeds and scores the runner-up's own knot against its own law at the truth.
  - `check_BD.py`, `mu_f0_f1.py`, `d_repolish.py`.
- **Probes at settings not run**
  - `probe_pairs.py` runs `f_job` itself, with new seeds, on (0.8, −0.6) at s = 0.01 and 0.03 and on (0.8, 0.2) at s = 0.01. Each has c̄ = 0.032, six separations, 200 H1 draws per separation and 400 H0 draws. `probe_analyse.py` and `probe_h0_fits.py` analyse them.
  - `z_null_snr.py` and `z_null_bias.py`: z's null against s/|ŵ|, with 4 000 and 20 000 draws per level.
  - The probe's s = 0.03 H0 cell had a one-sided size of 0.0975 (39/400). Fresh draws at the same null give 0.0535 (4 000) and 0.0512 (20 000), so it is a fluctuation.
- **Mutation**
  - 69 single-point mutants (`mut/mutants.py`). Each was run against the test file with `-x`, then `--from-raw` into a scratch `CF_OUT`, and the md compared (`run_mut.py`, `run_mut2.py`, `summary_all.py`).
- **Literature**
  - All literature checks rest on search snippets.
  - users.soe.ucsc.edu and apps.dtic.mil did not resolve, so the Shahram–Milanfar texts were not read.
  - The search budget ran out before the splitting-criterion query, so those items rest on memory.

## Verdicts per part

| claim (draft line) | verdict | evidence |
|---|---|---|
| (i) expansion, delay ∝ Δ⁻², T₃ remainder O(ϑ²), opposite-sign wall, dipole Δ⁻¹ (l. 7–20) | **holds** | Unchanged; raw A identical; rerun row identical; l. 140–147 match the md |
| (ii) λ_split, λ_sep, flow rates, the separation's two regimes (l. 22–40) | **holds** | Rows rerun bit-identical. "1.4–3.2% low" checks out: −1.44%, −2.36%, −3.23% |
| (iii) projected rates; the null law first-order exact (l. 42–48) | **holds** | D3bf; S and F reruns |
| (iii) "the transit's knots are the LARS knots, τ_k = 1/λ_k" (l. 45–46) | **holds for k ≤ 2 only** | D3bk-4 |
| (iii) null quantile moves ≈10× within a cell; near a tie it falls faster than 1/s (l. 48–51) | **holds** | q·s from 8.2 down to 0.46; 0.813, 0.594, 0.463 at c̄ = 0.049 |
| (iii) "its median set by far kernels and its lower quantiles by the fit's neighbours" (l. 48) | **median holds; lower quantiles offset-dependent** | Adjacent second entrants: middle 20% 0.09–0.34; earliest 5% 0.34–0.52 at c̄ ≤ 0.02, 0.76–0.80 at 0.032, 1.00 at 0.045–0.049 |
| (iii) calibrated knot keeps its size, "conservatively near a tie" (l. 51–52) | **valid; conservative more broadly** | D3bk-7 |
| (iii) same-sign power ½ at Δ² ∝ (\|A\|/\|A₁A₂\|)s/√N (l. 53–55) | **holds** | (0.8, 0.2) against (0.8, 0.6): knot Δ₅₀ ratio 1.56, predicted 1.46, NP ratio 1.49 |
| (iii) below NP "by what the knot adds to the runner-up's own filter" (l. 55–57) | **mechanisms right; the comparison partly credits the sign oracle** | D3bk-1, D3bk-8 |
| (iii) opposite-sign knot biased, least power 0.000–0.025 (l. 57–59) | **holds** | F; stronger at (0.8, −0.6) |
| (iii) "its power need not even grow with Δ" (l. 59) | **as measured; the cause is the moving null** | D3bk-6 |
| (iii) both walls at α/2 "removes most of the bias" (l. 59–60) | **fails** | D3bk-2 |
| (iii) walls' Δ₅₀ 1.707–1.840 "against the score's 0.965–1.058" (l. 60–61) | **numbers right; compares two-sided with one-sided** | \|z\| gives 1.059–1.179 |
| (iv) z ~ N(0, 1) under H0 (l. 64) | **holds** | 42 200 draws, as above. \|z\| size 0.0475 [0.0455, 0.0496]. At α = 0.01: 0.0090 (\|z\|) and 0.0089 (one-sided) |
| (iv) z's validity range | **remark** | The O(s) term is a mean shift: E[ŵz] = −0.034 ± 0.007 at s/\|ŵ\| = 0.14 and −0.057 ± 0.007 at 0.21. At 0.21 the one-sided sizes are 0.0575 toward the opposite sign and 0.0429 toward the same sign. F's range s/\|ŵ\| ≤ 0.07 is safe |
| (iv) H1 mean ‖f₁ − f₀‖/s·(1 − O(ϑ²)) (l. 64–65) | **holds** | 1 − cos = ϑ²/12 |
| (iv) "with the sign of A₁A₂" (l. 65); proof "V/A's sign" (l. 105) | **wrong for A < 0** | D3bk-5 |
| (iv) one-sided test at NP "for either sign"; \|z\| (l. 65–66) | **holds with the sign known** | Above; D3bk-1 |
| (iv)(a) width reading (l. 67–69, 106–109) | **holds** | Re-derived; five pairs checked |
| (iv)(b) relative-mode rate (l. 69–71, 109–110) | **holds as an identity** | D‴.54's forcing term vanishes at the fit. Notation per D3bk-5; prior art per D3bk-10 |
| (v) orders s/Δ², s/Δ³; CRB attained; split at its bound to ≈24Δ³; jump at ≈12Δ² (l. 73–78) | **holds** | md D; position CRB ∝ s/Δ² (0.0070); re-polish check |
| Status: Smith's criterion "∝ s/Δ … Δ² ∝ s" (l. 120–121) | **false for unequal amplitudes** | D3bk-3 |
| New (l. 133–137) | **partly** | D3bk-1, -8, -10 |
| Numerics (A)–(F) (l. 139–196) | **match the md/raw** | Two misquotes (D3bk-11) |

## Status of the round-59 issues

| D3bf | status | evidence |
|---|---|---|
| 1 (major): near-NP claim; matched filter | **fixed, with new overreach** | "Near-NP" is withdrawn and the two mechanisms are named. The score sits at the bound (E′ −0.057 … +0.052; F 1.003 and 1.011). The knot-vs-score comparisons ignore the knot's α allocation and the score's sign (D3bk-1). The decomposition is asserted, not measured (D3bk-8) |
| 2 (major): opposite-sign bias | **partly** | F measures the bias (least power 0.000–0.025 in 24 of 24 cells), and (iii) is restricted to same-sign pairs. The walls "remedy" is still biased (D3bk-2). The non-monotonicity is misattributed (D3bk-6) |
| 3 (major): novelty, literature | **mostly fixed** | τ₂ = 1/λ₂ is proved correctly. Spacing and covariance tests, McLachlan, Harris, Helstrom, Liu–Nehorai, KMW, Neyman–Scott, ORXYY, Amari et al., Dwivedi, Heinrich–Kahn are cited, and New is narrowed. New problems: τ_k for every k (D3bk-4); the Smith derivation (D3bk-3, already present and missed at round 59); the prior art of (iv)(b) (D3bk-10) |
| 4: offset dependence; noise-scale tie | **fixed** | F gives Δ₅₀ per offset; (iii) l. 50–51 |
| 5: bootstrap conservative near a tie; s known | **fixed in text** | The conservativeness is broader than stated (D3bk-7) |
| 6: D's jump variable | **fixed** | 9 cells with s/Δ² ≥ 12 at 0.293–0.409; 31 cells with ≤ 9.8 at ≤ 0.152 |
| 7: split below its bound | **fixed** | 0.84–1.05 up to s/Δ³ = 24, 0.79 at 29, 0.62–0.67 at 30–80 |
| 8: C's 200-draw thresholds | **fixed** | C′ (law at the truth): false alarms 0.0483 [0.0437, 0.0532]; Δ₅₀ with intervals |
| 9: misquotes | **fixed** | 17.8/21.4/11.9; 13 of 15; 0.027–0.044; column labels; 500.7. Two new misquotes (D3bk-11) |
| 10: tests | **partly** | D3bf's M07, M09–M12, M21, M33 (my X01, X03–X08) are now killed. M08 still survives (X02): with M_E = 64, the p at the draw's fit equals the p at the truth on all 4 test draws (laws differ ≤ 4%). See D3bk-9 |
| 11: cosmetics | **fixed** | |

## Issues

| id | severity | issue and evidence | proposed fix (where) |
|---|---|---|---|
| D3bk-1 | major | **The knot is a de facto one-sided test, and the draft's comparisons with the score ignore that.** <br>**Allocation.** Under H0, 92% of the knot's rejections in F fall on the wall of the fit's sign: 422/461 for (0.8, 0.6), 396/431 for (0.8, −0.4), 75–100% per offset. <br>**Same sign, no sign knowledge.** Against the sign-agnostic score \|z\|, the calibrated knot has the smaller Δ₅₀ in 17 of 24 F cells (median Δ₅₀/NP 1.064 vs 1.085; at α = 0.01, 1.037 vs 1.057). It has the larger power in 10 of 15 E cells (differences −0.175 … +0.103). Also at (0.8, 0.2): 1.069 vs 1.120. <br>**Same sign, one-sided.** Against the one-sided z (1.003), the one-sided own-wall knot (p_A < α, already in the raw) is at 1.056 [0.987, 1.289]. <br>**Opposite sign.** The same allocation is the bias: the drift empties the wall where α was spent. <br>**What the draft does.** l. 60–61 sets the two-sided walls test (1.707–1.840) against the one-sided z (0.965–1.058). Like for like, the comparisons are \|z\| at 1.059–1.179, or the one-sided other-wall knot (p_B) at 1.621–1.749 against z. The title's "detects the pair at the Neyman–Pearson bound" holds only with sign(A₁A₂) known. Without it, \|z\| costs 1.04–1.18 × NP. | **Title (l. 1–3) and (iv) l. 65–66:** add "with the pair's sign known". <br>**(iii) l. 53–61:** state the knot's α allocation. Compare one-sided with one-sided and two-sided with two-sided. <br>**F:** print the per-wall share of the null rejections and the one-sided p_A/p_B tests. <br>**New l. 135–136:** reword "the knot's losses against the efficient statistic". |
| D3bk-2 | major | **Both walls at α/2 stay biased.** <br>**Pooled power.** Over the 24 opposite-sign cells, at the two smallest separations (4 800 draws each): 0.036 [0.031, 0.042] and 0.039 [0.034, 0.045]; calibrated per draw, 0.034 and 0.040. NP there is 0.190–0.257 and 0.438–0.591, and \|z\| has 0.129 and 0.377. <br>**Per cell.** Least power 0.010–0.045, below α in 24 of 24 cells. <br>**At (0.8, −0.6), s = 0.01.** 0.035, 0.035, 0.015, 0.040 at Δ = 0.005–0.014, against NP 0.10, 0.18, 0.42, 0.86. <br>**Mechanism.** Under H1 the own-sign wall rejects 0.0054 and 0.0015 (0.0223 under H0); the other wall only 0.031 and 0.038 (0.0235 under H0). The draft's own 0.025–0.070 (l. 60, 193) shows it. | **(iii) l. 59–60 and F l. 193–194:** "both walls reach power ½ sooner (1.71–1.84 × NP) but remain biased (pooled 0.036 and 0.039 at the two smallest separations); only the score is unbiased for both signs". <br>**F:** print each test's least power with a pooled CP interval. |
| D3bk-3 | minor | **Smith's criterion is mis-derived** (Status l. 120–121). <br>**The error.** The bound of c₂ − c₁ is [(JᵀJ)⁻¹]_δδ, not 1/λ_sep. The Δ⁶ split moves the separation by (Δ/2)(1/w₂ − 1/w₁) per unit weight, since (ii) and (v) say it moves the centres by O(Δ). <br>**Unknown amplitudes.** sd(ĉ₂ − ĉ₁) = 1.98·10⁻³ s/Δ² for (0.8, 0.6) and 1.78·10⁻² s/Δ² for (0.8, −0.4); the slope is −1.997. Smith's criterion sd = Δ then gives s ≈ 506Δ³ and 56Δ³: the support-recovery order, not Δ² ∝ s. <br>**When it is ∝ s/Δ.** Only for equal amplitudes ((0.8, 0.8): 2.02·10⁻² s/Δ, s ≈ 49.6Δ²) or known amplitudes (1.92·10⁻² s/Δ). | **Status:** restate with the amplitudes' role, or keep only Liu–Nehorai's √CRB statement as cited. |
| D3bk-4 | minor | **"τ_k = 1/λ_k" holds for the first two knots only** (l. 45–46; Status l. 116 "the LARS identity proved"). The proof (l. 99–101) covers a one-atom block. <br>**300 noisy targets.** Knots 1–3 agree to 2·10⁻¹³. The 4th differs in 5/300 (up to 17%). The 5th differs from LARS in 65/300 and from the lasso in 68/300, by up to 95% and 280%. <br>**Why.** The ISS's sign-constrained block drops atoms at entry times; the lasso drops at its own knots; LARS never drops. | "The first two knots are the LARS (= lasso) knots; later ones agree while the active fits stay sign-consistent." |
| D3bk-5 | minor | **The sign of z, and the notation.** <br>**Sign.** The mean of z is ⟨½Vψ, R̂∂²φ⟩/(s‖·‖), with the sign of V, which is sign(A₁A₂)·sign(A). l. 65 ("with the sign of A₁A₂") and l. 105 ("V/A's sign") are wrong for A < 0: (−0.8, −0.6) gives z < 0 and ŵz > 0. `rejections` correctly multiplies z by sign ŵ, but every tested pair has A > 0, so mutant R01 (the sign dropped) survives. <br>**Notation.** ⟨r̂₁, ∂²φ⟩ is a Euclidean sum in z (l. 63) but the mean of l. 5 in (b) (l. 70, 110). In "dD/ds … s‖R̂∂²φ‖/N", s is both a time and the noise. | **(iv):** "ŵz has the sign of A₁A₂". **Proof:** "V's sign". Define both brackets; rename D‴.54's time. **Tests:** add a pair with A < 0. |
| D3bk-6 | minor | **The non-monotone power (l. 59, 191–193) is the null moving across the tie, not the drift.** <br>**Mechanism.** The closest single bump drifts by ĉ − c̄ ≈ −T₃/(4Aσ²) ∝ Δ³. <br>**s = 0.01, c̄ = 0.032.** At Δ = 0.0125 … 0.07, ĉ = 0.0322, 0.0325, 0.0335, 0.0358, 0.0417, 0.0537 runs up to and past the tie at 0.05. q_law·s goes 7.18, 7.25, 7.07, 6.59, 4.83, 1.97, while the H1 knot keeps falling (median 1281, 1393, 1223, 813, 408, 205). Hence the power 0.875 at Δ = 0.05 drops to 0.315 at 0.07. <br>**s = 0.03.** q·s = 5.86, 5.42, 2.98, 1.36 (power ≈ 0), and 7.32 at Δ = 0.1, past the tie (power 1). <br>**Labels.** F labels opposite-sign cells by c̄, but their nulls sit at c̄ + 0.0003, + 0.0015 and + 0.0097 at the H0 separations (s = 0.003, 0.01, 0.03). | **(iii) and F:** give this mechanism. Label opposite-sign cells by the null's ĉ, or choose c̄ so that ĉ sweeps the cell. |
| D3bk-7 | minor | **The calibrated knot is conservative beyond the tie** (l. 51–52). <br>**α = 0.05.** Pooled F size 0.0414 [0.0386, 0.0443] (0.0421 if p ≤ α). For the opposite-sign pair's null, 0.0379 [0.0342, 0.0419]. <br>**α = 0.01.** E 0.0071 [0.0058, 0.0085], F 0.0076 [0.0064, 0.0089]: 24–29% below nominal. <br>**Per cell.** 7 of 48 F nulls are at ≤ 0.020 (CP upper bounds 0.032–0.039), with ĉ up to 0.018 from the tie. <br>**For comparison.** The law at the truth gives 0.0465 [0.0435, 0.0495]. | **(iii):** "conservative (pooled 0.041 at α = 0.05, 0.008 at α = 0.01), most so near a tie"; report the α = 0.01 sizes. |
| D3bk-8 | minor | **Mechanism statements against measurements.** <br>**(a) Lower quantiles.** l. 48's "lower quantiles by the fit's neighbours" holds at c̄ ≥ 0.032 only. In the earliest 5% of H0 draws the second entrant is a neighbour in 0.34–0.52 of draws at c̄ = 0 and 0.02. l. 176 ("far kernels set the tail") is ambiguous. <br>**(b) "Measured" losses.** New calls the losses "measured", but the repo measures only the totals and z₂. My decomposition supports the mechanisms. Far from the tie, the runner-up's own knot (against its own law, size 0.039–0.064) equals z₂: at c̄ = 0, 0.647/0.640, 0.380/0.378, 0.552/0.542. The dictionary's minimum then costs 0.07–0.12 (full knot 0.530, 0.282, 0.438). At c̄ = 0.049 the runner-up's knot itself loses: 0.260 against z₂'s 0.400, and 0.223 against 0.530. | **(iii) l. 48:** make the wording offset-dependent. **E′:** add the runner-up's-knot column (`e_decomp.py`). |
| D3bk-9 | minor | **Tests: the new data paths are untested.** <br>**Counts.** Of 69 mutants, 13 are killed by the tests only, 17 by `--from-raw` only and 6 by both. 33 survive both; 8 of those are equivalent on this data (Z03, S02, E03, E04, G04, U01, U02, and R01, which is latent for A < 0). <br>**Non-equivalent survivors (25):** <br>– wall identity: walls swapped (K01, W01, D01); j₁ not excluded (W02; harmless only for the opposite-sign pair the tests use); <br>– score and width: Z04 (projector at c + 0.01), Z02 (t without √(N−2)), D02 (z_σ's sign), S03 and S05 (width escape); <br>– S's data path: all of `s_job`'s E and H0 paths (J01–J08: H0 drawn from the pair, z₂'s wall, swap, μ, qA, law at s/2, shifted fits, cos); <br>– `part_S`: seeds per hypothesis (P01), E's fit order (P02); <br>– F's data path: `f_job` H0 from the pair (F01), H0 with the H1 seed (F02), qA at α (F03), NP at s/2 (F04); `f_jobs` H0 separation (F07); <br>– D3bf's M08 (X02). <br>**Report sections.** C′/E′ are never called by a test (C01–C03, E01–E02 and E05 are caught only by `--from-raw`). `cp_int` is tested against itself. | **Add tests for:** <br>– the walls on a same-sign and an opposite-sign pair (tA = Δτ₂ < tB, and the reverse); <br>– `s_job` on a C H0 job and on a 2+2-draw E cell; <br>– `part_S` with a spy on its seeds and fits; <br>– `f_job` with hyp = 0; <br>– z against an independent formula to 10⁻¹⁰, and z_σ ≈ sign(ŵ)z; <br>– a spy on `knot_law`'s (c, w) in `e_job`; <br>– `report(C=…, S=…, E_=…)` on synthetic S rows. |
| D3bk-10 | minor | **Literature and novelty of (iv).** <br>– (b) is the splitting criterion at the overlap singularity: a duplicated unit's critical line is a minimum or a saddle according to the sign of the residual-weighted second derivative (Fukumizu–Amari 2000; Liu–Wu–Wang 2019's "splitting matrix") [memory]. <br>– The one-vs-two locally optimal test is Shahram–Milanfar's [snippet]. <br>– Neyman–Scott is dated 1965 in the snippet (Bull. ISI 41(1), 477–497). <br>– "The paradigm computes it twice" (l. 67) overstates (b): it needs a merged two-kernel block, and the transit here never forms one (the second entrant enters at its own grid kernel). | Cite these. In New, call (b) "the classical splitting criterion in the paradigm's metric". |
| D3bk-11 | cosmetic | **Misquotes and round-off.** <br>– l. 191 "ten times the bound's 0.01001": > 0.04 is 4.0×. <br>– l. 183–184 "below its noise (sd 0.137–0.167s)": true only at s = 0.03 (μ = 0.0036, sd 0.0041). At s = 0.01, μ = 0.0035 is 2.1 sd above its noise (sd 0.0016), consistent with 1.6% swaps. <br>– A row (0.8, −0.4, 0.0125) predates solver commits b3f434f/76eb48e and reproduces only to 2.4·10⁻¹⁴ in w; rms 3.45·10⁻¹⁶ against the md's 4.09·10⁻¹⁶. <br>– State z's tested range s/\|ŵ\| ≤ 0.07. | Correct the text; rerun A or name the solver commit. |

## Sources

All from search snippets; no page was read in full.
- [Gu–Koenker–Volgushev, Testing for homogeneity in mixture models (arXiv 1302.1805)](https://arxiv.org/pdf/1302.1805): C(α) tests against mixtures, and the score's singularity.
- [Springer chapter citing Neyman–Scott, "On the use of C(α) optimal tests…", Bull. ISI 41(1):477–497 (1965)](https://link.springer.com/chapter/10.1007/978-981-99-6141-2_9)
- [Liu–Nehorai, Statistical angular resolution limit for point sources, IEEE TSP 55(11), 2007](https://profiles.wustl.edu/en/publications/statistical-angular-resolution-limit-for-point-sources/): the limit is "proportional to the square root of the CRB".
- [Shahram–Milanfar, Imaging below the diffraction limit (IEEE TIP 2004)](https://users.soe.ucsc.edu/%7emilanfar/publications/journal/DiffLimitFinal.pdf) and [IEEE TIT 2006 (DTIC copy)](https://apps.dtic.mil/sti/pdfs/ADA460918.pdf): locally optimal tests for small separations. Both hosts were unreachable.
- [Kulaitis–Munk–Werner, What is resolution? (arXiv 2005.07450)](https://arxiv.org/pdf/2005.07450)
- [Osher–Ruan–Xiong–Yao–Yin, Sparse recovery via differential inclusions (arXiv 1406.7728)](https://arxiv.org/pdf/1406.7728): the ISS path can be unbiased, while the lasso path is biased.
- [Tibshirani–Taylor, The solution path of the generalized lasso (arXiv 1005.1971)](https://arxiv.org/pdf/1005.1971): LARS without modification differs from the lasso path.
- Memory only: Fukumizu & Amari, Neural Networks 13 (2000) 317–327; Liu, Wu & Wang, NeurIPS 2019 ("Splitting steepest descent"); Smith, IEEE TSP 53(5) (2005).

## Lesson

Where a test spends its α is part of what it is. The knot never chose its allocation: the grid's margins put 75–100% of its null rejections on the wall of the fit's sign. So it is a one-sided test toward same-sign pairs, as good there as the sign-agnostic score and blind to opposite-sign pairs. Round 59 found the bias. It then compared a two-sided remedy with the one-sided score, and judged that remedy by Δ₅₀ rather than by the property it was meant to restore: power at least α at small Δ. Compare like with like (one-sided with one-sided, sign known or not), and re-test every fix against its own defect.
