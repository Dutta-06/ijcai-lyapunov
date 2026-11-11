# Referee report D3au, round 48: Proposition D‴.58 (the transit under width mismatch)

**Overall verdict.**
- **What holds.** The exact algebra of (i) is correct; I checked it independently to 30+ digits. The status's ₂F₁ conjecture is honestly
  stated and better supported than the draft says.
- **What does not.**
  - K* and the "non-optimal stationary point beyond K*" in (iii) are not properties of the dynamics. They come from a precision floor,
    the polish's acceptance rule, the starting grid, and silent LSODA failures in `relax_chunk`.
  - (iv)'s oracle reading is contradicted by its own data.
  - (v) is wrong as stated.
  - The novelty paragraph misplaces the closest precedents.
  - The tests do not protect the experiment's analysis layer.

## What I checked and how

- **Setup.** I worked in a scratch copy and did not modify the repository. Every process ran with OMP_NUM_THREADS=1 and nice 10, one heavy
  process at a time. The only process I killed was my own tail-ratio job, by explicit PID.
- **Reproduction.** `exp_width_mismatch.py --from-raw` with WM_OUT in scratch regenerates `results/chaos/width_mismatch.md`
  byte-identically. Both test files pass (9 tests).
- **(i), independently, in high precision.** mpmath at 80–160 digits, written from scratch:
  - Gauss–Hermite nodes by Newton's method;
  - the stated rule's L² error from closed-form Gaussian integrals;
  - the defects D_{n,K} computed from the nodes, not from the repository's integer code;
  - τ/σ = ¼, ½, 1, 2, 3 and K ≤ 20;
  - every constant of the proof by hand.
- **The family of scaled Gauss–Hermite rules** "GH of N(0,v) for the reweighted integrand", with closed-form errors up to K = 60, both
  with reweighted GH weights and with least-squares weights. I compared:
  - moment matching, v = τ²;
  - half-squeeze, v_hs = τ²σ/σ*;
  - the natural Fock/Taylor scaling, v_hs/cosh r;
  - Mercer-zero nodes τ²(1+2τ²/σ²)^{-1/2} (Fasshauer–McCourt parametrisation, from memory);
  - the best v for each K.
- **(v).** Analytic continuation to σ* < σ, checked by closed form and by mp.quad.
- **The ₂F₁ conjecture.**
  - P_j(K) exactly for j ≤ 10, K ≤ 40.
  - Part D's tail ratios recomputed independently in closed form for K ≤ 80.
  - Part D extended at t = 0.8 to K = 640 and 1280 with exact integers.
- **C(c) on ℝ** in closed form after the continuum optimal fits.
- **Part B, from the raw:**
  - relaxed flags and events;
  - analytic gradient and Hessian at every first non-optimal block;
  - the reduced objective along the segment from block to optimum;
  - variable-projection LM started from the block;
  - the experiment's own fast flow continued from the recorded states, with solve_ivp's status inspected.
- **Part B reruns**, changing one setting at a time: polish_step 1σ→2σ (all six widths), gtol 1e-9→1e-13, and the starting grid shifted
  by +0.05.
- **Part C, per draw from the raw:** stop index, bias/variance split, ratio to the oracle error.
- **Mutation testing:** 42 mutants (25 in `wide_feature.py`, 17 in `exp_width_mismatch.py`) run against the two test files.
- **Literature.** Web search only, so abstracts and snippets. ar5iv.labs.arxiv.org and deepai.org did not resolve, so I stopped trying
  other sites. I did not read Karvonen–Särkkä 2019 (KS2019), Karvonen–Oates–Girolami 2021 (KOG2021), Derevianko–Plonka (D–P) or
  Plonka–Riebe–Cuyt 2026 (PRC2026) in full. Anything marked [memory] is from memory.

## Verdicts

**(i) Correct, with small gaps in rigour.**
- What checks:
  - Nodes: y_k = (2σ² sinh r)^{1/2}z_k = τ(σ/σ*)^{1/2}z_k; the two forms agree to 70 digits.
  - Amplitudes a_k as stated.
  - The series ‖f*−f_K‖² = ‖f*‖²(1−t²)^{1/2}Σ_{n≥2K}(tⁿ/n!)D²_{n,K} matches the closed-form error to 31–46 digits wherever my working
    precision resolves it.
  - D_{n,K} = 0 for n < 2K and for odd n; D_{2K,K} = 2^K K!.
  - The proof's constants: Bargmann prefactor and exponents, ±t for widths e^{±r/2}, β, the Fock integral as Mehler at ρ = t, the
    undoing of steps (a) and (b), c_T² = ‖f*‖² sech(r/2), the eigenvalue (1−ϑ)/(2(1+ϑ)), Stirling.
  - The first-term shares match the md (e.g. 0.99946 at τ/σ=1, K=1).
- Gaps:
  - Step (c) integrates Mehler's series termwise against Δ⊗Δ. The N(0,2) part has a density relative to N(0,1) that is not in
    L²(N(0,1)), so this step needs a justification. The same interchange is used in E_ϑ(K) = Σ(ϑⁿ/n!)D².
  - Convergence of the product Gauss–Hermite rule for κ_ϑ is asserted without a cited theorem.
  - "Exactly the rate t" is proved only as a K-th-root limit. The per-kernel ratio e_{K−1}/e_K → 1/t needs the conjectured constant.

**Status: the ₂F₁ conjecture is fair, and the evidence is stronger than stated.**
- P_j(K) is integer-valued and of exact degree j, with leading coefficient −C_{j−1}/2^{j−1}, for every j ≤ 10 on K = 1…40 (exact
  arithmetic). The j-th differences are −(2j−3)!!.
- Part D is reproduced to every printed digit by my independent closed form: t = 0.8 at K = 10, 20, 40, 80 and t = 0.5 at K = 10, 40, 80.
  Doubling the number of terms at K = 320 changes nothing.
- Extension at t = 0.8:

| K | gap to ₂F₁ | K·gap |
|---|---|---|
| 320 | −3.31e-3 | −1.058 |
| 640 | −1.909e-3 | −1.222 |
| 1280 | −1.049e-3 | −1.343 |

  Successive gap ratios fall 0.67 → 0.61 → 0.58 → 0.55, so the data are consistent with ₂F₁.

**(ii) Numerically right, but stated as a fact.**
- Moment matching's rate per kernel at K = 60 is 2.99978 (against 3) and 1.49994 (against 1.5).
- "2^K for narrow mismatch" overstates the gain. The actual gain is about 0.6·K^{−1/4}·(λ/t)^K: 1388 at K = 12 for τ/σ = 1/20, against
  2^12 = 4096.

**(iii) Part right, part refuted.**
- **δ_t is misquoted and depends on K.** From the raw, for K = 1…10:

| τ/σ | δ_t in the raw | draft says |
|---|---|---|
| ½ | 0.078→0.040% | 0.05% |
| 1 | 0.728→0.395% | 0.4–0.5% |
| 2 | 3.48→2.12% | 2–3% |

- **Nodes** agree to 0.003–5.1%, not "0.1–5%".
- **"Rate t" for the optimal fit** is unproved, as the status says.
- **A constructive point.** Retuning only the variance, v ≈ v_hs(1 − 0.03/K), recovers almost all of δ_t: e_best/e_hs = 0.99631 at
  τ/σ=1 and 0.98218 at τ/σ=2 (K = 60). This is within about 0.1% of the series' first term alone. At K = 10, τ/σ = 2 the optimal fit is
  only 0.1% below this rescaled rule.
- **C ≥ 0 with double zeros** is confirmed on ℝ: |min C/max C| ≤ 2e-9 for τ/σ = ½…3, K ≤ 8. It is a property of the optimum only; after
  the half-squeeze rule itself, min C/max C runs from −0.6 to −4e-6.
- **"Moving out with K"** fails at wide mismatch: the peak goes 2.89σ→2.62σ at τ/σ = 3 and 2.57σ→2.52σ at τ/σ = 2.
- **"From alternating flanks"** is violated, e.g. at τ/σ = 1 the third kernel enters from the same side as the second (kernels 11 then
  13).
- **"Forms the optimal fits up to K*; beyond, relaxes to a non-optimal stationary point; K* = 5, 6, 9, 11; K* grows with τ/σ"** is
  refuted as dynamics. The transit does form the optimal fits by flank insertion; what stops it is numerical (D3au-1 to D3au-3).

**(iv) The prediction holds; its interpretation is false.**
- Stops and errors match the prediction: stops 6/4/3 against 6/4/3; median errors 8.7e-5, 1.71e-3, 1.08e-2 against predicted 9.5e-5,
  1.77e-3, 1.06e-2.
- But "the oracle's s(2K/N)^{1/2}; a wrong width costs a factor K(s)^{1/2}" does not hold, because bias dominates:

| s | predicted bias e_K | variance term | measured error / oracle(K) | cost vs right-width oracle (claimed √K) |
|---|---|---|---|---|
| 1e-4 | 5.3e-5 | 7.9e-5 | 1.0–1.33 | 2.7 (2.45) |
| 1e-3 | 1.65e-3 | 6.4e-4 | 2.6–2.8 (seven K=4 stops) | 5.3 (2.0) |
| 1e-2 | 9.0e-3 | 5.6e-3 | 1.7–2.3 (K=3 stops) | 3.4 (1.73) |

- The cause: with δ = 0.1 the discrepancy rule tolerates bias² up to 0.21s², while 2Ks²/N is only 0.02–0.04s².
- K(s) ≈ log(‖f*‖/s)/log(1/t) gives 4.4, 3.1, 1.8 against the actual 6, 4, 3.

**(v) Wrong.**
- For σ* < σ the continued rule has conjugate imaginary nodes, and the approximant is a real function: e^{−x²/(2σ²)} times a cosine sum
  with Hermite-zero frequencies.
- Its L² error is exactly the same series at |t| = (σ−σ*)/(σ+σ*). I checked this to 10 digits by closed form and by quadrature at
  σ*/σ = 0.5 and 0.8, K ≤ 12.
- This is Derevianko–Plonka's class, so the frame does not fail for narrower features.

**Literature and novelty: see D3au-6.**

**Numerics paragraph against the md.** Every quoted number matches the md, except:
- the δ_t and node ranges above;
- "10⁻¹⁵…10⁻⁵": the maximum is 3e-5;
- "the trapezoid's floor near 10⁻¹³": the trapezoid–series difference is about 1e-17 absolute, coming from the double-precision GH nodes
  and weights; 1e-13 is only the cut-off;
- "double-precision floor ≈ 10⁻⁸ at ¼": the K=6 block there reaches 2.4e-9, so 1e-8 is not a floor on the error (see D3au-1 for the
  actual floor).

**Tests: 42 mutants, 23 killed, 19 survived (W09 is equivalent).**
- All mutants of the rules, `params`, `l2_error`, `gh_defect`, `mehler_series`, `tail_ratio` and `first_term` are killed.
- Surviving in the library:
  - W12: `optimal_fit` drops the centre node at odd K.
  - W15: `ls_amplitudes` with rcond = 1e-3.
  - W18: `gh_defects` corrupted only for K > 21, which leaves all of Part D at K ≥ 40 unprotected.
  - W14: the unused symmetric=False path.
- Surviving in the experiment, 14 of 17:
  - E01 `predict_stop` with K/N instead of 2K/N;
  - E04 a coarser flank grid;
  - E06 the stop rule without δ;
  - E07 the optimal-block threshold that defines K*;
  - E08 `blocks()`;
  - E09 `sample_fit` sorting;
  - E10 Part D truncated to 10 terms at t = 0.8;
  - E11 a vacuous degree check in `lead_coeffs`;
  - E12 the sign of the A2 lead coefficient;
  - E13 the Catalan verdict forced true;
  - E14 σ*;
  - E15 `part_A`'s starting nodes;
  - E16 the error measured against noisy data;
  - E17 the at-feature radius.

## Issues

| id | severity | issue | proposed fix (where) |
|---|---|---|---|
| D3au-1 | major | **K* is not a property of the dynamics.** (a) The first non-optimal blocks at τ/σ = ½, ¾, 3/2 are relaxed=False in the raw, which the draft never mentions. (b) LM from every first non-optimal block's own centres reaches the optimal fit. (c) Changing polish_step from 1σ to 2σ, nothing else, moves K* from 6 to 7 at ¾ and from 11 to ≥12 at 3/2, and changes the "beyond" ratios (e.g. 144→4.7 at ½, K=7). (d) Shifting the grid by +0.05 with default settings moves K* from 6 to 7 at ¾ and from 11 to ≥12 at 3/2. (e) With polish_step = 2σ, every optimum with e_opt ≥ 1.0e-7 is reached and every one with e_opt ≤ 6.3e-8 is missed (1.6e-9, 6.3e-8, 5.2e-8, 5.2e-8). That is a floor near 7√ε_mach, so K* ≈ ln(1/e_floor)/ln(1/t). (f) gtol 1e-13 leaves the τ/σ=1 stall at 9.24e-8. | (iii), numerics (B), status, novelty: drop K*'s values and "law" and say the transit forms the optimum wherever double-precision relaxation resolves it (here e_opt ≳ 1e-7). Report per-block relaxed flags and polish rejections. Show invariance to grid, polish_step, gtol and integrator before naming any K*. |
| D3au-2 | major | **`rich_transit.relax_chunk` ignores solve_ivp failures.** LSODA fails ("Unexpected istate", t_end = 0) for any chunk ≳ 1–2·10¹¹ from relaxed states. I checked six states at three widths; there are no successful chunks at ≥ 2.1e11. `relax_chunk` returns the unchanged input state as if the chunk had run. So the last 3–4 of MAX_DOUBLINGS = 24 are no-ops (the md comment "≈ 1.7·10¹²" is false), and the τ/σ=3/2 stall is a frozen failure mid-descent (error 1.956e-6→1.777e-6 when the solver died). This is cross-cutting for every user of `sliding_iss`. | `relax_chunk`: check `sol.success`, fall back to Radau or split the chunk; count and report failures. Rerun Part B. |
| D3au-3 | major | **"Beyond K* the block relaxes to a non-optimal stationary point; the path wanders through exits and re-entries; 1.6–550×".** At τ/σ=3/2 the "frozen" state has \|∇F\| = 4.2e-11 with λ_min(H) = 2.9e-12 > 0, the error falls monotonically along the segment to the optimum, and LM reaches it. At ½ and ¾ the stalled blocks are flagged unrelaxed. The "exits and re-entries" are enter/exit cycles at a single τ: kernel 16 at ¼ (13 events at 2.499e10, error rising 1.3771e-8→1.4055e-8), kernel 3 at ½ (11 events), kernel 4 at ¾. They come from exiting without re-relaxing. This is also outside D‴.53's hypotheses (H1) no exits and (H2) nondegenerate blocks; Hessian condition numbers here are ~1e13. | (iii), (B): remove these sentences or describe them as an implementation artefact. Fix the exit rule in `sliding_iss` (re-relax after an exit). |
| D3au-4 | major | **(iv): "the oracle's s(2K/N)^{1/2} … a wrong width costs K(s)^{1/2}" is contradicted** (numbers in Verdicts). Also, ‖f*‖ is defined as the L²(ℝ) norm (0.439) but the code uses the sample rms (0.245). Taken literally, the predicted stop at s = 1e-3 becomes 5, not 4. | (iv): state error² ≈ e_K² + 2Ks²/N with e_K² ≤ ((1+δ)²−1+2K/N)s²/‖f*‖_N², i.e. error ≲ s((1+δ)²−1+4K/N)^{1/2}. The oracle claim needs δ = O(K/N). Use the sample norm ‖·‖_N. |
| D3au-5 | major | **(v) is wrong.** The frame continues exactly to σ* < σ: a Gaussian times a Hermite-zero cosine sum, with the error equal to the same series at \|t\|. Derevianko–Plonka's construction (Hermite-zero frequencies, coefficients from GH weights, exponential decay in the weighted L² norm, per the snippets) is this side's precedent. | Rewrite (v) as an exact statement. Compare D–P's rate with \|t\|. |
| D3au-6 | major | **Novelty and literature.** (a) D–P is a direct twin, not a side remark. PRC2026 (arXiv 2606.03855: exponential sums from GH quadrature of the Gaussian's integral representation, geometric decay) is missing. (b) "Their scaling possibly the half-squeeze one" is likely wrong. The natural Fock scaling has rate λ/(2−λ) (5.003 and 2.007 per kernel, against 1/t = 5.83 and 2.62), and its series' first term gives exactly a bound e^{−c₁N}N^{1/4} with e^{−c₁} = λ/2, which is the shape of KOG2021's abstract with c₁ > c₂. The Mercer-zero scaling (KS2019-type, [memory]) has rate q = (β²−1)/(β²+1) (3.74 against 1/q = 3.73). So the half-squeeze is likely a strictly faster scaling; confirm this from the papers. (c) Missing references: Kuo–Woźniakowski 2012 and Kuo–Sloan–Woźniakowski 2017 (relevant to (iii)'s open lower bound); Tang 1993 on Hermite scaling of Gaussians [memory]; Steinwart–Hush–Scovel 2006 and Fasshauer–McCourt 2012 (Fock/Mercer); Janszky–Vinogradov 1990 and Szabo et al. 1996 (discrete coherent-state superpositions) [memory]. (d) "Stalls beyond K*" is refuted (D3au-1–3); "the noise trade-off" is Morozov's discrepancy principle, not new. | Literature and novelty paragraph. |
| D3au-7 | major | **The experiment's analysis layer is untested.** 14 of 17 mutants survive, including the K* classifier, the stop rule, Part D's truncation and the Catalan verdict. In the library, W12, W15 and W18 survive. | Tests: `report()` on a synthetic raw; `optimal_fit` at odd K; `ls_amplitudes` against a closed-form Gram solve; `gh_defects` at K = 25–40 against `gh_defect` or high precision; a degree check that can fail; a stop-rule test sensitive to δ and to 2K/N. |
| D3au-8 | minor | δ_t and node ranges wrong at small K, and δ_t treated as constant (Verdicts, (iii)). | (iii) and numerics (A): quote δ_t(K); add the rescaled-rule observation as an upper bound. |
| D3au-9 | minor | Proof gaps: the termwise integration in (c) and in E_ϑ(K); convergence of the product GH rule; "exact rate" holds only as a K-th root. | Proof of (i): a Fock-space Parseval argument, or Cramér's bound with a slightly wider Gaussian; cite an Uspensky-type theorem. |
| D3au-10 | minor | (ii) is asserted unconditionally in the proposition; "2^K" should be ≈0.6K^{−1/4}(λ/t)^K. | Statement (ii). |
| D3au-11 | minor | (iii) wording: C ≥ 0 holds only for the optimum, not for the rule's residual; "moving out with K" and "alternating flanks" are not general. | (iii). |
| D3au-12 | minor | The "optimal fit" is the best of two LM starts in Part B and one symmetric start in Part A; it is not certified as global. | (iii) and (B): say "best local optimum found", or add multistarts. |
| D3au-13 | cosmetic | "Integer polynomials" should be integer-valued (P_2 = −(K²−25K−6)/2). "10⁻¹⁵…10⁻⁵" should be 3e-5 at most. The trapezoid "floor near 10⁻¹³" is really ~1e-17. "Double-precision floor ≈ 10⁻⁸" is misstated. The status understates its own P_j evidence (j ≤ 10, K ≤ 40, j-th differences −(2j−3)!!). | Numerics paragraph, md captions in `exp_width_mismatch.report`, status. |

## Lesson

The exact algebra was right and could be checked to 40 digits in an afternoon. Every serious problem sits in the interpretive sentences
built on numerics at the edge of double precision. A number reported as a property of a dynamics, here K*, must be shown invariant under
every numerical setting before it is interpreted: precision, tolerances, acceptance rules, solver success flags and the starting grid. In
this run a single floor at e_opt ≈ 1e-7, a polish acceptance radius, a grid offset and silent LSODA failures account for every K*.
Likewise, a novelty paragraph should first try the one-line continuations of the result. Taking σ* < σ would have exposed the
Derevianko–Plonka twin, and computing the competing scalings' rates would have settled the KOG2021 question.

Sources (search snippets only): [KOG2021 arXiv:2004.12654](https://arxiv.org/pdf/2004.12654),
[KS2019 arXiv:1803.09532](https://arxiv.org/pdf/1803.09532), [Derevianko–Plonka arXiv:2307.13587](https://arxiv.org/pdf/2307.13587),
[Plonka–Riebe–Cuyt arXiv:2606.03855](https://arxiv.org/pdf/2606.03855), [KSW2017 via unpaywall](https://unpaywall.org/10.1090%2FMCOM%2F3144),
[Gnewuch–Ritter–Rüßmann arXiv:2412.05368](https://arxiv.org/pdf/2412.05368),
[squeezed vacuum from Rabi interactions arXiv:2010.12411](https://arxiv.org/pdf/2010.12411).
