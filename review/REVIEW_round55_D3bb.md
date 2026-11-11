# Referee report — round 55, D3bb: Proposition D‴.52 (close features)

**Scope.** I reviewed these files at fa25a4f:
- `review/drafts/D52.md`
- `experiments/exp_close_features.py`
- `src/core/rich_transit.py` (`sliding_iss` with `max_entries`/`max_rhs`, `critical_rhs`, `polish_block`)
- `results/chaos/close_features.md` and `close_features_raw.npz`
- `tests/test_close_features.py`

Other sessions moved HEAD to 849f5fe during the review, but `git diff fa25a4f HEAD` on all of these files (and on the helpers they import) is empty. I modified nothing in the repository. All work ran in a scratch copy, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/d3bb/`, with `nice -n 10` and `OMP_NUM_THREADS=1`, never more than two processes at a time.

## Overall verdict

**What holds — the deterministic core, which is stronger than the draft shows:**
- **The path is the exact-residual law.** In all 18 monopole rows of A, the sliding path's Δτ₂ equals "the duals integrate the exact one-bump residual" to ≤ 3.3·10⁻⁷. So every deviation in table A is expansion error, not path error.
- **(i)'s building blocks check out.** The ψ-identity, the opposite-sign wall with its 2 − μ margin, and the dipole's Δ⁻¹ law all hold.
- **The O(ϑ) reading of (0.8, −0.4) is confirmed.** Adding the T₃ term turns the ratios 0.907 and 0.948 into 1.012 and 1.003.
- **(ii)'s algebra and spectra hold.** The determinant −W²ϑ(1−ϑ)δ⁴, the three continuum constants, λ_split, λ_sep and r_split all check out. The flow is exactly −DH.
- **The explanation of the separation excess is right.** Scaling the amplitude mobilities by 0.1 removes it.
- **The order claims of (iii) and (iv) are right, and more:**
  - At the tested offset the detection power equals the Neyman–Pearson oracle bound within Monte Carlo error.
  - Where the pair is resolved, the placement errors are at the Cramér–Rao bound.
- **Reproducibility is excellent.** The md regenerated from raw is byte-identical, and fresh reruns of cells from A, B, C and D are bit-identical.

**What does not hold:**
- **(iii) is false as written.**
  - The noise law is wrong: for the kernel that matters the std is 0.108 of the one stated.
  - The H0 mechanism is wrong: 83% of H0 second entrants are kernels far from the feature.
  - The test is calibrated at the true position. Its null quantile changes 13× within one grid cell, and a threshold set at c̄ = 0.032 gives 45% false alarms at c̄ = 0.049.
  - So "the transit counts the pair" is an oracle-calibrated statement, not a procedure.
- **Novelty is overstated.** Using a path knot as a spike-counting statistic has close prior art (Azaïs–De Castro–Mourareau; covariance and spacing tests). The resolution limits credited to "Liu–Ammari" are Liu–Zhang's.
- **Smaller defects.**
  - Two error terms are misstated: (1 + O(ϑ)) in (i) and r_sep's (1 + O(Δ)) in (ii).
  - (i) lacks a no-tie hypothesis.
  - The saturation in Part D is mis-described.
  - Six numbers are misquoted.
  - 19 non-equivalent mutants survive the tests.

## What I checked and how

Scripts are in the scratch directory above.

- **Reproducibility.**
  - `--from-raw` into scratch, compared with `cmp`: byte-identical.
  - Fresh reruns compared with raw:
    - A rows (0.8, −0.4, 0.00625), (0.8, 0.8, 0.00625) and (0.8, −0.8, 0.0125).
    - B rows (0.8, 0.6, 0.0125) and (0.8, −0.4, 0.2).
    - C cell (0.025, 0.01): both hypotheses, 400 draws.
    - D cells (0.05, 0.03) and (0.0125, 0.03), the second including the ±700 outlier.
    - All identical to raw (rms differs only at the 10⁻¹⁶ level).
  - pytest: 6 passed in 7.7 s.
- **(i)** (`checkA.py`).
  - Re-derived the moments: V, T₃, and M₄ = (A₁A₂Δ⁴/A)(1 − 3A₁A₂/A²).
  - Checked the identity ⟨ψφ(c̄+d)⟩ = G(d)d²/(4σ⁴). With ℓ = N·h it matches the sampled values to 4 digits, and ⟨ψφ_j⟩ ≥ 0 on all 21 kernels.
  - Built two predictions and compared them with the path:
    - an exact-residual prediction (least-squares one-bump fit, v_j = ⟨r₁φ_j⟩);
    - an O(ϑ) prediction (½Vψ + ⅙T₃ψ₃, with ψ₃ = ∂³φ projected off span(φ, ∂φ)).
- **(ii)** (`checkB.py`).
  - Re-derived the determinant symbolically (W = w₁ + w₂, ϑ = w₂/W, δ = c₂ − c₁, centroid fixed).
  - Re-derived the split and separation images, and the integrals √π/(2σ³), 3√π/(4σ⁵) and the 3/2 ratio.
  - Computed the eigenvalues of −D·JᵀJ analytically and compared them with b_job's finite-difference Jacobian.
  - Scaled the amplitude mobilities by κ = 1, 0.1, 0.01, and went down to Δ = 0.0005.
  - Read the split eigenvector's centre components directly.
- **(iii)**
  - Projected noise std per kernel (`checkC_noise.py`).
  - Which kernel enters second under H0, over all 8000 H0 draws (`inspect_C.py`).
  - Reran cell C(0.025, 0.01) with the experiment's seeds at c̄ = 0, 0.02, 0.032, 0.045, 0.049 (`cpos.py`, `ana_cpos.py`).
  - Neyman–Pearson oracle power Φ(‖f₁ − f₀‖/s − z₀.₉₅) for every cell (`np_bound.py`).
  - Bootstrap of Δ₅₀, its slope, and the realized size of the threshold (`boot.py`).
  - N = 61, 121, 241 (`cN.py`).
- **(iv)**
  - Cramér–Rao bound (s²/N)(JᵀJ/N)⁻¹ at the truth for every Δ, against the measured rms (`crb.py`).
  - What the placed blocks look like (`inspect_D.py`).
  - Unpolished vs polished placement on three noisy draws (`flow_vs_polish.py`): same end every time.
- **Mutation** (`mutants.py`, `run_mut.py`, results in `mut/results.json`).
  - 41 single-replacement mutants.
  - Each was run against `tests/test_close_features.py`, plus `test_noisy_certificate.py` for the solver mutants.
  - Each was also checked for whether the md regenerated from raw changes.
- **Literature.** Checked by web search where possible (arxiv.org itself is unreachable); items from memory are marked [memory].

## Verdicts per part

**(i) — holds, with a wrong error term and a missing hypothesis.**
- The path equals the exact-residual prediction: path/exact ∈ [0.99999971, 1.00000000].
- The O(ϑ) explanation is right. For (0.8, −0.4) the leading order gives 0.854, 0.907, 0.948 at Δ = 0.025, 0.0125, 0.00625; with the T₃ term these become 1.0526, 1.0123, 1.0029. The remainder falls 4.2× per halving, i.e. O(ϑ²).
- ϑ is not the only expansion parameter. For A₁ = A₂ (ϑ ≡ 0), path/leading − 1 = 0.184, 0.043, 0.0106, 0.0026, 0.00066, 0.00017, which is 0.061(Δ/σ)² (D3bb-4).
- The constant depends on the grid offset, and the law fails at ties (D3bb-5).
- The dipole law holds: path/shape = 1.0086, 1.0021, 1.0005 → 1, at O(Δ²).
- "The path's end is the exact fit" is numerical: true on the 24 polished paths. The unpolished flow reaches the same block on the three draws I tried.

**(ii) — holds, except r_sep's error term.**
- Holds:
  - Smallest Hessian eigenvalue over Δ⁶: 91.34 against 91.23.
  - Split rate against its prediction: 1.0054, 1.0013, 0.9993 at Δ = 0.025, 0.0125, 0.00625.
  - Flow = −DH: agreement to 4–5 digits.
  - The split eigenvector's centre components equal Δ/(2w_k) to 0.1%.
- The λ_sep shortfall of −2.4, −3.2 and −1.4% at Δ = 0.0125 is the amplitude-sum mode leaking into the Euclidean (w, c) metric: −(Δ²/8σ⁴)·w₁²w₂²/(w₁²+w₂²) predicts −2.2, −3.0, −1.2%. This is O(Δ²), consistent with the claim.
- The r_sep formula is valid only in a window around Δ ≈ 0.01. Below it the rate tends to 2/3 of the formula (D3bb-6).

**(iii) — order right, near-optimal at the tested offset, statement wrong, calibration missing.**
- **Order.** Δ² ∝ s/√N matches Shahram–Milanfar (SNR ≈ c/(N²d⁴)) and Liu–Zhang's exponent 1/(2n−2).
- **N-scaling.** At (0.025, 0.01):

  | N | H0 median | power |
  |---|---|---|
  | 61 | 1321 | 0.42 |
  | 121 | 1848 | 0.48 |
  | 241 | 2637 | 0.71 |

  The H0 median grows ×1.40 and ×2.00, i.e. ∝ √N, as predicted.
- **Near-optimal.** The Neyman–Pearson oracle (simple vs simple, the same f₀ and f₁) gives Δ₅₀ = 0.0152, 0.0279, 0.0484 against the measured 0.0153, 0.0274, 0.0510. Its cell-by-cell power matches the md's power table within Monte Carlo error. The reason: the runner-up's projected kernel Rφ₁₁ ≈ (d²/2)P⊥∂²φ is nearly parallel to the residual r₁ ∝ ψ, so its dual increment is close to a matched filter. This is worth stating in the draft.
- **Wrong as stated / uncalibrated:** see D3bb-1 and D3bb-2.
- **Is H0 a fair null?** As a simple null, the L²-closest single bump is the least favourable one, so it is fair for power. But it sits at the true offset (D3bb-1).
- **Are 200 draws enough?** See D3bb-10.

**(iv) — holds, with the Cramér–Rao bound actually attained; saturation mis-described.**
- Cramér–Rao constants:

  | Δ | 0.3 | 0.2 | 0.14 | 0.1 | 0.07 | 0.05 | ≤ 0.035 |
  |---|---|---|---|---|---|---|---|
  | position rms × Δ²/s | 0.0116 | 0.0088 | 0.0078 | 0.0074 | 0.0072 | 0.0071 | 0.0070 |

  The split rms × Δ³/s runs 0.0074 → 0.0068. The s/Δ² law is asymptotic, reached for Δ ≲ 0.1.
- In the 22 cells where the split's Cramér–Rao std is ≤ 0.2, the measured errors are 0.84–1.18 (positions) and 0.67–1.05 (split) of the bound. With 40 draws that is equal to it within ±11%.
- The split is lost at s ≈ 15–75Δ³: its Cramér–Rao std reaches 0.1 at s ≈ 15Δ³. This is the support-recovery order (Batenkov–Goldman–Yomdin SRF^{2ℓ−1} [memory]; Denoyelle–Duval–Peyré t^{2N−1}).
- The draft's description of the saturation is wrong (D3bb-7).

**Status and proof.**
- "(ii) the singular-value orders proved" refers to a classical fact: the confluent Vandermonde / Prony map.
- The determinant gives only the product of the singular values; the individual orders need a column-norm argument.
- The proof of (iii) at l. 63–64 is contradicted by the data (D3bb-2).

**Numerics paragraph.** Every number in (A), (B) and (C) matches the md, except the misquotes in D3bb-8. For (D), see D3bb-7.

**Tests.** All 6 pass, but 24 of 41 mutants survive:
- 5 are equivalent: M10, M11, M12, M14, M35. (M11 is equivalent even in b_job, because ⟨(t²−1)²⟩ = ⟨t⁴⟩ under e^{−t²}.)
- 19 are not (D3bb-9).

**Literature.** See D3bb-3.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3bb-1 | major | **(iii)/(C): the test is calibrated at the true position (a simple null at c̄ = 0.032), and Δτ₂'s null law depends on the feature's offset from the dictionary grid.** I reran cell (Δ, s) = (0.025, 0.01) with the experiment's seeds at c̄ = 0, 0.02, 0.032, 0.045, 0.049. q₀.₀₅·s = 8.16, 8.19, 8.02, 3.50, 0.63. The c̄ = 0.032 threshold gives false-alarm rates 0.04, 0.04, 0.05, 0.20, 0.45. Power at each offset's own threshold is 0.28, 0.34, 0.42, 0.38, 0.29, against a Neyman–Pearson oracle of 0.375 at every offset. Part C ran only at the tested offset with the highest power, with an oracle null. "The transit counts the pair" is therefore not a calibrated procedure. | (iii) l. 37–42 and Status: say the null is simple and oracle-calibrated. Calibrate on the first fit's (ĉ, ŵ) by parametric bootstrap, or derive the null law. In `c_job`/`part_C`: sweep c̄ across a grid cell and report size and power at each offset. Add the Neyman–Pearson column Φ(‖f₁ − f₀‖/s − z₁₋α). |
| D3bb-2 | major | **(iii)'s noise law and H0 mechanism are false.** Var⟨ξ̃φ_j⟩ = s²⟨(Rφ_j)²⟩/N, where R projects off the first fit's tangent span(φ, w∂φ); it is not s²⟨φ_j²⟩/N (l. 38). At c̄ = 0.032 the std ratio is 0.108 for j₂ = 11, 0.35 for 9, 0.51 for 12, and 1.00 for far kernels. The proof (l. 63–64) says "the second entry under H0 is the first of the kernels near the feature". In fact only 16.6% of the 8000 H0 second entrants are neighbours of j₁. Neighbours do make up 70% of the earliest 5%, which sets q_α, but the median H0 delay is set by far kernels. The criterion at l. 41 omits the margin μ_{j₂} and the projection. At the measured Δ₅₀, ½\|V\|⟨ψφ₁₁⟩ equals 1.56–1.80 × s(⟨(Rφ₁₁)²⟩/N)^{1/2}, close to z₀.₉₅ as for a matched filter. Against the draft's s(⟨φ²⟩/N)^{1/2} it is only 0.17×. So C depends on the grid offset and the dictionary, not on the pair alone. | Rewrite (iii): the rates are ½V⟨ψφ_j⟩ + N(0, s²⟨(Rφ_j)²⟩/N), and the H0 delay is min_j μ_j/\|n_j\| over the whole dictionary. Power is ½ at ½\|V\|⟨ψφ_{j₂}⟩ ≈ z₁₋α·s(⟨(Rφ_{j₂})²⟩/N)^{1/2}. Note that Rφ_{j₂} is nearly parallel to r₁, hence the near-optimal power. Fix the proof at l. 62–64. |
| D3bb-3 | major | **Literature and novelty (l. 42, 47, 69–75).** (a) Counting spikes from path knots is prior art. Azaïs–De Castro–Mourareau (arXiv 1706.00679) test the number of spikes in super-resolution from the first two LARS knots, with a calibrated null. See also the covariance test (Lockhart–Taylor–Tibshirani–Tibshirani, AoS 2014 [memory]) and the Kac–Rice/spacing tests (Taylor–Loftus–Tibshirani 2016 [memory]). The ISS gives τ₁ = 1/λ₁ exactly, and Δτ₂ is its analogue of the second knot. Round 45 (D3ar-8) already flagged knot-based tests. (b) The exponents (σ/m)^{1/(2n−2)} and (σ/m)^{1/(2n−1)} are Liu–Zhang's (arXiv 1912.05430, 2003.02917; multi-dimensional 2103.11632), not "Liu–Ammari 2022–24". Liu–He–Ammari (2211.13541) is the positive-source sequel. (c) Missing references: Smith 2005, the CRB-based resolution limit — with (ii)'s std(Δ̂) ∝ s/Δ it gives Δ² ∝ s and links (ii) to (iii). Demanet–Nguyen 2015 and Li–Liao 2021 [memory]. The Prony-map / confluent-Vandermonde conditioning behind (ii)'s orders: Batenkov–Yomdin, Akinshin–Batenkov–Yomdin 2015, Batenkov–Demanet–Goldman–Yomdin SIMAX 2020 [memory]. ADCG. Sliding Frank–Wolfe should be named as the algorithmic analogue of "enter, then polish". | Limit "New" to: the ISS delay law with its constant and its offset dependence; the near-Neyman–Pearson efficiency; and the mobility-weighted rates. Compare with the knot tests' calibrated nulls. Correct the attributions and add the references. |
| D3bb-4 | minor | **(i)'s error term (l. 7, 8, 12) is wrong.** ϑ is not the only expansion parameter. For A₁ = A₂ (ϑ ≡ 0), path/leading − 1 = 0.184, 0.043, 0.0106, 0.0026, 0.00066, 0.00017 for Δ = 0.2 … 0.00625, i.e. 0.061(Δ/σ)². The remainder O(\|A₁A₂\|Δ⁴/\|A\|) is not uniform in A: M₄ carries the factor 1 − 3A₁A₂/A², which is 7 for (0.8, −0.4). | Write (1 + O(ϑ) + O((Δ/σ)²(1 + \|A₁A₂\|/A²))). Optionally add the T₃ term ⅙T₃⟨ψ₃φ_j⟩; it takes (0.8, −0.4) to 1.012 and 1.003 at Δ = 0.0125 and 0.00625. |
| D3bb-5 | minor | **(i) has no non-degeneracy (no-tie) hypothesis.** The "explicit constant" is Δτ₂Δ² = μ_{j₂}/(½(V/Δ²)⟨ψφ_{j₂}⟩), so it depends on the grid offset. At Δ = 0.025 it is 0.873, 0.654, 0.276, 0.065 for c̄ = 0, 0.032, 0.045, 0.049. At an exact tie both kernels enter at τ₁, so Δτ₂ = 0 for every Δ. | (i) l. 8–13: assume a unique j₁ with μ_{j₂} bounded away from 0, and state how the constant depends on the offset. Novelty at l. 73. |
| D3bb-6 | minor | **(ii)'s r_sep (l. 31–33): "(1 + O(Δ))" is wrong.** The frozen-amplitude error is O(μ⟨φ²⟩/r_sep), which grows like Δ⁻². From exact −DH rates for (0.8, 0.6), rate/formula = 1.029, 1.150, 0.506, 0.663 at Δ = 0.0125, 0.00625, 0.003, 0.0005. Across the three pairs it is +12 to +26% at 0.00625 and −44 to −53% at 0.003. As Δ → 0 the rate tends to Gm·λ_sep (projected off e₀ too), which is 2/3 of the formula. The draft's explanation of the excess at 0.0125 is correct: multiplying μ_k by 0.1 takes +2.90/+2.35/+4.67% to +0.04/+0.02/+0.09%. | State both regimes and the crossover Δ ≈ 0.003–0.004, and drop "(1 + O(Δ))" for r_sep. |
| D3bb-7 | minor | **(iv) / numerics (D), l. 92–96: the saturation is mis-described.** (1) The claim that the positions saturate at 0.3–0.4 from s ≳ 10²Δ³ (and at Δ = 0.025 from s = 0.003) is wrong. The position rms is 0.052 at (0.025, 0.003) (s/Δ³ = 192), 0.065 at (0.0125, 0.001) (512) and 0.152 at (0.0175, 0.003) (560). The 0.3–0.4 plateau appears only at s/Δ³ ≈ 240–1900. (2) The plateau is not "the pair's scale". In these blocks one kernel carries ≈ A = 1.40 at c̄ and the other sits 0.1–0.9 away with \|w\| ≲ 0.04, fitting noise; the second entrant is often a far kernel. Only 9/40 blocks at (0.0125, 0.003) and 5/40 at (0.0125, 0.01) have both kernels within σ. (3) The Cramér–Rao constants (0.0070 and 0.0068) are never stated, though they are attained (Verdicts, (iv)). Where the sign constraints bind, the split rms falls far below the bound — 0.64× at (0.025, 0.001), 0.13× at (0.0125, 0.001) — so those rows do not confirm s/Δ³. | In report D, add Cramér–Rao columns, the number of blocks with both kernels within σ, and the second entrant. Rewrite the sentence: the split is lost at ≈ 15–75Δ³; the positions jump once the second kernel fits noise. |
| D3bb-8 | minor | **Misquoted numbers.** l. 86: "2.4–3.3% low" should be 1.4–3.2% ((0.8, −0.4) is −1.44%); "slope 1.89–1.93" should be 1.84–1.93 ((0.8, 0.8) is 1.837). l. 88: "2.9–4.7% above" should be 2.3–4.7% ((0.8, 0.8) is +2.35%). l. 92: "0.006–0.011" should be 0.006–0.012 ((0.3, 0.01) is 0.0116). l. 96: "w₁ + w₂ ≈ 0" is wrong: the outlier has w = (−698.6, +700.0), so w₁ + w₂ = 1.397 ≈ A, with centres 5.6·10⁻⁵ apart and entrant kernel 7 — a coalesced near-dipole fitting noise. l. 82–83: the md's ratio is path/shape, not the shape's over the path's. l. 79–81: at Δ = 0.2 the ratios 0.988 and 1.188 compare different kernels (predicted 9 vs path 11; predicted 8 vs path 15). | Correct the text. |
| D3bb-9 | minor | **The tests are weak: 24 of 41 single-replacement mutants survive.** 5 are equivalent. The 19 non-equivalent survivors: M4 (the opposite-sign wall branch dropped from `predicted_delay`, so the claimed "other wall" is untested); M17/M18 (the split-mode overlap ignores the centres, even reversed); M20–M22 (`detection`); M23–M27 and M41 (the report's ratios, slopes, C₁ and D metrics); M28 (H0 target := H1 target); M29 (no noise); M30 (`max_entries=None`, so delay runs relax without bound; the test passes 4× slower); M32 (`one_fit` takes one step); M36 (w not sorted with c); M37 (unrelaxed blocks counted as placed); M38 (τ₂ used for Δτ₂). M4, M17, M18, M28–M30, M32 and M36–M38 do not even change the md from raw. | In `tests/test_close_features.py`, add tests for: (0.8, −0.4) entrant 13 at wall −1 with ratio > 0.9 at Δ = 0.0125; H0 ≠ H1 targets and noise present; Δτ₂ = τ₂ − τ₁; w sorted with c; the relaxed flag; a small-Δ delay run with bounded nfev; `detection` on a synthetic table; an md-from-raw byte comparison; the split eigenvector's centre components. |
| D3bb-10 | minor | **Statistical precision is not reported.** From a bootstrap over 2000 resamples of both hypotheses per cell, the Δ₅₀ 95% intervals are [0.0139, 0.0176], [0.0253, 0.0302] and [0.0475, 0.0537]. The slope's is [0.46, 0.57] (quoted as "0.523 (½)"). A 200-draw 5% threshold applied to the other H0 cells at the same s gives false-alarm rates of 0.5–14% (median 5.3%). | In numerics (C) and the md: give the intervals, and use ≥ 10³ H0 draws for q_α or an analytic null. |
| D3bb-11 | cosmetic | (1) The proof (l. 52–53) reuses ϑ for w₂/W, and W and δ are never defined. (2) l. 56 says the images are orthogonal to span(e₀, e₁), but e₂ is not ⊥ e₀ (the P⊥ in the formula handles it). (3) l. 24 takes ℓ to be the interval's length. It must be N·h = 2.033: the constants are then exact to 10⁻⁵, whereas ℓ = 2 is 1.6% off (test l. 28–31 uses ℓ = 2). (4) "The singular-value orders proved" is a classical fact. (5) "In every case the path's end is the exact two-kernel fit" (l. 16) is numerical: 24 polished paths. (6) The md's "overlap 1.000" cannot see the split's centre components, though they do equal Δ/(2w_k). (7) In the tests the `__main__` block (l. 81–84) comes before the last test. | Define the notation; set ℓ = Nh; mark numerical claims as numerical; report the eigenvector components; move the main block to the end of the file. |

## Lesson

A statistic is only as good as its null. Here the null law is governed by a nuisance parameter the experiment held fixed: the feature's offset from the dictionary grid. Two cheap checks would have shown both sides at once:
- **A one-line Neyman–Pearson benchmark**, Φ(‖f₁ − f₀‖/s − z₁₋α), shows the strength: the power sits at the oracle bound.
- **A sweep of the offset** shows the weakness: q_α swings 13× within one grid cell.

Two corollaries:
- When an expansion has two small parameters, test the case where the advertised one vanishes (here A₁ = A₂).
- Entry-time and knot statistics have a calibrated prior literature, as round 45 (D3ar-8) already noted. Search it before claiming novelty.

**Sources (web-checked):**
- [Liu–Zhang, computational resolution limit, line spectra (arXiv 2003.02917)](https://arxiv.org/pdf/2003.02917)
- [Liu–Zhang, one dimension (arXiv 1912.05430)](https://arxiv.org/pdf/1912.05430)
- [Liu–Zhang, multi-dimensions (arXiv 2103.11632)](https://arxiv.org/pdf/2103.11632)
- [Liu–He–Ammari, positive sources (arXiv 2211.13541)](https://arxiv.org/abs/2211.13541v1)
- [Azaïs–De Castro–Mourareau (arXiv 1706.00679)](https://arxiv.org/pdf/1706.00679)
- [Batenkov–Goldman–Yomdin (arXiv 1904.09186)](https://arxiv.org/pdf/1904.09186)
- [Positive near-colliding sources (arXiv 2212.00536)](https://arxiv.org/pdf/2212.00536)
- [Milanfar group, SNR ≈ c/(N²d⁴) (OSA SRS 2005)](https://users.soe.ucsc.edu/%7emilanfar/publications/conf/OSA-SRS-05.pdf)
- [Shahram–Milanfar, IEEE TIT 2006 (DTIC ADA460918)](https://apps.dtic.mil/sti/pdfs/ADA460918.pdf)
- Smith 2005, IEEE TSP 53(5) (bibliographic listing only)
- Everything marked [memory] is unverified.
