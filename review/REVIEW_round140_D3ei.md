# Referee report D3ei, round 140: the Paper 17 manuscript, first review

**Manuscript:** manuscript/paper17/paper17.tex, compiled to paper17.pdf (33 pp.). Title on the paper: "The small-mass transit of self-adaptive kernel networks: one kernel per feature, the sliding limit and exact recovery". Repository commit 9398242, branch claude/vigilant-clarke-hf939b. Sources: Propositions D‴.48 and D‴.51, Theorem D‴.53 with Lemma D‴.53a, Theorem D‴.55.
**Venue considered:** SIAM J. Math. Data Sci. or Found. Comput. Math.
**Workspace:** W = /tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref140. The repository was not modified.

## Summary of the submission

The paper studies Gaussian kernel networks. Each coefficient is w_i = β_i a_i/m, with a free mass β_i updated multiplicatively, a tanh-saturated amplitude a_i, and centres whose drift is proportional to the mass (eq. 2.2). The networks start from masses εb_i and zero amplitudes.

- **§2.** A conservation law Q_i, a mirror coordinate η_i with dη/dt = g_i, a dual potential, and the "wall" λ_{1,i} = Λ(ε)ρ_i(1+O(ε)), where Λ ∝ ε^{-1/2} (Lemma 2.1).
- **§3, fixed centres.** Theorem 3.2: the rescaled duals converge to the inverse scale space (ISS) flow of Σρ_i|w_i| on the time scale ε^{-1/2}. The proof uses Mosco convergence of an ε-dependent family of dual potentials, then Attouch and Brezis–Pazy.
- **§3, branches and creep.** Proposition 3.3 handles starts on branches, with a logarithmic scale and one-sided boxes. Proposition 3.4 gives the core coefficient and the creep speed of the centres.
- **§3, formal limit.** Formal statement 3.5 is the "sliding ISS": kernels that have not entered are frozen, and entered kernels slide onto the features.
- **§4, finite mass.** An ε-free inner profile at the wall (Proposition 4.1). An ε-free near-tie inner problem with a Gibbs split and a plateau tail (Proposition 4.2). The scales of the creep (Proposition 4.4). An estimate of the largest mass that gives exact recovery (Estimate 4.5).
- **§5, the sliding limit.** Under (H1)–(H4), Theorem 5.3 proves convergence to the sliding path at rate ε^{1/2}, with an explicit first-order correction: the pulls Δ_kj and the entry shifts σ_k. The remainder is O(ε^{3/4}). Lemma 5.2: the inner problem is a gradient flow of the fit.
- **§6, recovery.**
  - One feature: the path is exact (Proposition 6.2, continuum).
  - Two flow-free "skeletons" of the path (Proposition 6.3) and the two ends of the amplitude ratio (Proposition 6.4).
  - Amplitude-uniform exactness for separated features is left as Conjecture 6.5, with a proof sketch.
- **§7, numerics.** Taken verbatim from four md reports. The headline: the fitted skeleton is right on 1154/1154 targets, the amplitude-free certificate on 903.

The status of every result is given once, in Table 1. This is good practice, and on the whole the statuses are honest.

## Overall assessment (recommendation with reasons)

**Recommendation: major revision.**

**Strengths.**
- The mechanism is clean: the mass multiplies the drift, so kernels that have not entered are frozen. The two-scale picture is new as far as I know: power-law transit, sliding block, and ε^{1/2} inner layers whose fit is a Lyapunov function.
- Most of the analysis checks out. I re-derived these exactly:
  - Lemma 2.1(i)–(iii);
  - the proof of Theorem 3.2;
  - Proposition 4.1's expansions;
  - Lemma 5.2(i);
  - Proposition 6.2.
- Theorem 5.3's first-order prediction matches the loop to a median ratio of 1.0000 (md).
- The numerics are reproducible: all four md reports regenerate byte-identically from their raw files, and every number I sampled matches.

**Problems that need a revision.**
- **The entry shift is wrong for unequal masses.** σ_k in Theorem 5.3(c) lacks a factor 1/ρ_{j_k}. I resolved README inconsistency 1 by derivation: the sign in the paper is right and the ρ-factor is missing.
- **Theorem 5.3(b) is not proved at the left end of its window.** The proof's own Step 5 bound is O(ε^{1/4}) there, and Step 2 uses the wrong time origin.
- **(H3) is misdescribed.** The inner problem is called a "gradient flow in another metric". In fact it is the path's own relaxation, reparametrised by t = Λ₀s. So (H3) is the convergence of the relaxation that Definition 2.3 already presupposes, and Lemma 5.2(iii)'s last step is vacuous rather than wrong.
- **Proposition 3.4's error term is not sharp.** The stated O(ε^{1/2}/μ) should be the sharper −(4b/(π²ς))ε/μ². I confirmed this analytically and numerically, and the paper's own creep data fit it to four digits. §7.2 currently says the data are "consistent with O(ε^{1/2}/μ)", which the data contradict.
- **The recovery headline is largely circular.** "Right on all 1154 targets" mostly restates the basin hypothesis (Proposition 6.3(c) concedes that (c) is a tautology). The title ("one kernel per feature … exact recovery") promises more than is proved: one feature in the continuum, plus a conjecture.
- **The literature comparison misses the closest works on several fronts:**
  - Chizat's global result for moving particles;
  - continuous OMP exact recovery (Elvira–Gribonval–Soussen–Herzet);
  - positive spikes without separation (Schiebinger–Robeva–Recht);
  - Jacot et al. on saddle-to-saddle dynamics;
  - basins for non-convex spike estimation (Traonmilin–Aujol);
  - Berthier–Pillaud-Vivien's Theorem 4.2, which already gives Proposition 3.3's paired (two-sided) case.
- **The paper is not self-contained.** It depends on an unavailable companion [comp16] (Proposition 3.3; §1.1). It gives no motivation for the law, and at 33 pp. it is long.

**On venue.** For FoCM the rigorous core is narrow: Theorems 3.2 and 5.3, the latter conditional on (H3), which is proved only for small overlaps. SIMODS is the better fit, but it would want baselines (BLASSO with η_V, SFW/ADCG, continuous OMP) and a motivation.

## Major comments

**M1. The entry shift σ_k (Theorem 5.3(c), l. 486, p. 18; README inconsistency 1): resolved.**

*Derivation.* On the path the entrant's dual satisfies s_k u_{j_k}(τ_k) = ρ_{j_k}. Its approach rate in absolute units is s_k du/dτ = s_k⟨r_{k−1}φ⟩ = ρ_{j_k}v_k, because v_k (l. 437) is divided by ρ. The loop's dual is u + ε^{1/2}dev. Solving s_k(u + ε^{1/2}dev) = ρ to first order gives the arrival time:

  τ − τ_k = −ε^{1/2} s_k dev_{j_k}/(ρ_{j_k}v_k),

so

  **σ_k = −s_k dev_{j_k}/(ρ_{j_k} v_k) = −s_k dev_{j_k}/|a_{k−1,j_k}|.**

*Sign.* The paper's s_k is correct. The code (experiments/exp_transit_limit.py l. 126: `sk = -(np.sign(P['U'][k][j]) or 1.0) * dev[j] / v`) uses that sign. Its v is the absolute |⟨r_pre φ⟩|, so the code is right only because ρ = 1. The header of transit_limit.md ("σ_k = −dev_{j_k}/v_k") omits s_k; that is a typo in the report.

*Consistency.* (d) already carries ρ in the pulls (s_j ρ_j Δ_kj, "absolute units"), and round 41 (D3an-10) put "v_k with its ρ" into the definition of v_k. σ_k was not updated with it.

*Fix.*
- Correct (c).
- State in Theorem 5.3 whether general masses are covered. Equal masses are the only case run.
- Add one unequal-mass target to tests/test_transit_limit.py so the ρ-factor is exercised.
- Also write b_j for b in κ(μ) (l. 476).

**M2. Theorem 5.3(b) at the left end of its window (ll. 485, 500, 506; pp. 18–20).**

*What (b) claims.* The block is within O(ε^{3/2}μ_min^{-3} + ε^{3/4}v_{k+1}^{-3}) of S_k on [τ_k + Cε^{1/2}log(1/ε), τ_{k+1} − ε^{1/4}].

*What the proof gives.*
- Step 3 gives a perturbation of the inner solution of O(ε^{1/2}(X+S)) + O((1+g′²/v³)/X), which is O(ε^{1/4}) at X = ε^{-1/4}.
- Step 5 says the next outer interval starts "within O(e^{−λS}) + O(ε^{1/2}(X+S)) of S_k in the stable directions", which is again O(ε^{1/4}).
- Step 2's Tikhonov bound is written C e^{−λ_k(τ−τ_k)/ε^{1/2}}|y(τ_k+ε^{1/2}S) − S_k|. Its exponent measures time from τ_k, but the initial datum is taken at τ_k + ε^{1/2}S. With the correct origin, at τ = τ_k + ε^{1/2}S the bound is O(ε^{1/4}), not O(ε^{3/2}).

So the left end of (b) is not covered as written. The bound holds only because the ε^{3/4}v^{-3} term is uniform over the window, and that term physically belongs to the right end.

*Data (transit_limit.md Part A, column "window start").*
- At the window start the block's distance falls like ε^{0.80} on target 33 (8.42e-5 to 8.80e-9 over five decades) and like ε^{0.91} on target 41. That is far from ε^{3/2}.
- The md places the window start at C = 30. The theorem requires C ≥ 3/(2λ_k), which reaches 56.6 (target 33, λ₃ = 0.0265), 49.7 (target 41) and 45.2 (target 10). Seven of the eight Part-A targets have an entry with 3/(2λ_k) > 30.

*Fix.* Either
- (i) start the window at τ_k + (C + 5/(4λ_k))ε^{1/2}log(1/ε) and give a two-term bound, e^{−λ_k(τ−τ_k−ε^{1/2}S)/ε^{1/2}}·O(ε^{1/4}) + O(ε^{3/2}μ_min^{-3}) + O(ε^{3/2}(v_{k+1}(τ_{k+1}−τ))^{-3}); or
- (ii) prove that the error in the stable directions at s = S is O(ε^{1/2}). This looks plausible, because the translation error is multiplied by dy/ds = O(e^{−λS}).

Also report the window-start row in Table 6 (l. 709) and rerun Part A with C ≥ 3/(2λ_k).

**M3. (H3), Lemma 5.2(iii) and Definition 2.3 (ll. 205, 445–452, 464, 472; pp. 7, 16–17).**

*The identity.* In (α, c) coordinates the inner problem of (H3) is exactly the fast flow (2.1)–(2.2) on the critical curve, with t = Λ₀s:
- On the critical curve, η_c′(α) = coth²(α/a)/(2Gς) = 1/(Gβa′). From dξ/ds = −(s/ρ)g and η_c = −Λ₀ρξ one gets dα/ds = Λ₀Gβa′g = Λ₀α̇.
- dc/ds = Λ₀Gm w⟨r∇φ⟩ = Λ₀ċ.
- Equivalently, the w-mobility of Lemma 5.2(i) is |𝒲′|/ρ = Λ₀Gβa′w_c′.

*Consequences.*
- The closing sentence of Lemma 5.2(iii)'s proof, "the same gradient flow in another metric from the same start … ends at the same point", is wrong in its premise. It is the same flow in other coordinates.
- (H3) is therefore the convergence of Definition 2.3's relaxation (with the entrant leaving w = 0⁺ along the unique-modulo-translation trajectory). Definition 2.3 already presupposes this convergence.
- Lemma 5.2(iii)'s real content is convergence and localisation near (S_{k−1}, (A, c*)) for small overlaps.

*Gaps in the proof of (iii).*
- The comparison "at a fixed s₁ the coupled flow is within O(δ) of the decoupled one" needs control over (−∞, s₁]. The far-field argument of Step 3 is needed here too.
- The use of Lemma 5.2(ii) makes the continuum an implicit hypothesis; it should be stated.

*Fix.*
- State the identity.
- Rephrase (H3) as "the relaxation of Definition 2.3 converges (to a nondegenerate critical point S_k)".
- Remove "another metric".
- Make the continuum and the far-field control explicit in (iii).

**M4. Proposition 3.4's error term (ll. 285–299, p. 10; §7.2 l. 658; README inconsistency 4): resolved.**

*The proof.* It defines δ = α²/a² + ε + α/(Gςλ₁μ), with α/a = O(ε^{1/2}/μ) and α/(Gςλ₁) = O(ε/μ). Hence δ = O(ε/μ²) + O(ε), not "O(ε^{1/2}/μ)": the step on l. 299 confuses α/a with α²/a².

*The exact core coefficient [own].* I solved η(α) = (1−μ)λ₁ exactly (bisection, double precision) for ε = 10⁻⁶ to 10⁻¹⁰, μ = 0.004 to 0.1, and b = 1, 2. The relative error times μ²/(εb) tends to −16.21 = −4/(π²ς), while the relative error times μ/ε^{1/2} tends to 0.

*Interpretation.* The leading correction is exactly Proposition 4.1's −8/(κξ)², with ξ = μ/ε^{1/2} and κ = πρ√(2ς).

*Data [md].* With this constant the creep ratios of mirror_transit.md Part C are reproduced:

| μ (at ε = 10⁻⁸) | predicted ratio | md ratio |
|---|---|---|
| 0.00392 | 0.9895 | 0.9896 |
| 0.00434 | 0.9914 | 0.9915 |
| 0.00549 | 0.9946 | 0.9946 |
| 0.00646 | 0.9961 | 0.9961 |

The deficit scales like μ^{-2}, not μ^{-1}: the ratio of deficits is 2.67, against (0.00646/0.00392)² = 2.72.

*Fix.* State w_i = profile × (1 − (4b_i/(π²ς))ε/μ² + O(ε/μ) + O(ε²/μ⁴)), under the hypothesis ε^{1/2} ≪ μ (now missing from the statement). Change l. 658 accordingly.

**M5. Proposition 4.2 (ll. 370–388, p. 13) and its status in Table 1 (l. 95).**

- **(b).** "The critical points with both coefficients nonzero are the exact fits" is asserted without proof.
  - In a symmetric configuration one can show it: ⟨r∂_cφ⟩ = 0 forces Γ(2d) = 1.
  - A random search I ran (1-D continuum, 20 000 starts, w₁, w₂ > 0) found no other critical points.
  - Give a proof, or state it as a hypothesis or as numerical.
- **(c).** S = 1/(1+e^{κx}) is the exponential-branch asymptote. The exact split is 𝒲(ξ+x)/(𝒲(ξ)+𝒲(ξ+x)), with 𝒲(ξ)+𝒲(ξ+x) = |A|. Table 1's "the Gibbs split exactly" overstates this.
- **(d).** I_∞ = ∫(⟨rφ_e⟩ − ⟨rφ_r⟩)ds/ρ lacks the factor s_A. "I_∞ > 0" is not proved. While the entrant grows at its start, ⟨rφ_e⟩ drops faster than ⟨rφ_r⟩, so the integrand is negative early on. Positivity is numerical (md: I∞ = 1.42–1.60 for A = 0.1–2).

**M6. Exact recovery: the headline and the predictors (abstract l. 52, title l. 41, §6, §7.5 ll. 720–750).**

- **The headline restates the basin hypothesis.** The fitted skeleton is "the sign-constrained least-squares fit … started at the features". The computed paths are "the sliding ISS with a bounded least-squares polish of each block". By Proposition 6.3(c) the two coincide exactly when the basin hypothesis holds. So "predicts correctly on all 1154" mostly measures the basin hypothesis, and the abstract should say so. The linear skeleton (1149/1154) is the genuinely predictive object.
- **The polish-off check covers only part of the set.** It was done for 640 of the 1154 targets (A and D(b)), not B, C, D(a) or D(c).
- **D(c) is selected.** Its 30 targets were chosen with |ν_p| > 0.8q (transit_certificate.md l. 71), but the paper (l. 725) does not say so. D(c) adds 180 correct predictions to every column, the certificate's included.
- **The certificate baseline is weak.** The comparison is with an ad hoc amplitude-free p(c_j), not with the vanishing-derivative pre-certificate η_V of [DP15] or the BLASSO outcome. The linear skeleton's projection onto span{φ(c*), ∇φ(c*)} is exactly the Gram structure of η_V. Computing η_V and BLASSO on the same targets would show where the path and the convex program differ. That is the natural SIMODS question.
- **The title overclaims.** "One kernel per feature … exact recovery" is proved for one feature in the continuum; for several features it is conjectural or numerical. Retitle, for example "… the sliding limit and a predictor of recovery".

**M7. Novelty and literature (§1.3, ll. 110–127, pp. 4–5).** See also the Literature section below.

- **Theorem 3.2 is an instance of the Berthier–Pillaud-Vivien mechanism** [BPV26], which I read in full. Two points belong in the text:
  - BPV's Theorem 3.8 obtains only convergence weighted by s ds and warns that "unweighted convergence should not be expected". The paper's unweighted L² convergence of Φw^ε is legitimate because ζ(0) = 0 ∈ int P and D_ε(0) = D_∞(0) = 0. That is a genuine (small) advantage and should be spelled out, together with the energy identity of the limit flow (l. 265).
  - BPV's Theorem 4.2 (asymmetric diagonal networks: the ISS of a weighted ℓ¹ norm) already gives Proposition 3.3's last sentence.
- **Chizat 2022** [Chi22] is the closest precedent for particles with multiplicatively updated masses and moving positions, applied to sparse deconvolution. It is cited only for the WFR remark. The comparison should be result by result: global convergence to the BLASSO minimiser (with regularisation) there, against the unregularised sliding ISS path here.
- **AGF** [KMCKSDGM25]: the dormant neurons' utility ascent corresponds to the WFR variant of Remark 3.6, where massless kernels climb the correlation peaks. Under (2.2) they are frozen instead. This is worth stating.
- **Missing works:**
  - Elvira–Gribonval–Soussen–Herzet (continuous OMP recovers in exactly k steps for any amplitudes, for admissible kernels), the natural analogue of Conjecture 6.5;
  - Schiebinger–Robeva–Recht (same-sign Gaussian spikes recovered without separation by weighted BP), against Remark 6.7's same-sign midpoint failure of the path;
  - Jacot et al. 2021 (saddle-to-saddle dynamics in deep linear networks);
  - Traonmilin–Aujol (basins of attraction for non-convex spike estimation), relevant to Lemma 5.2(iii) and the basin hypothesis.
- **The GSD20 analogy (l. 111) is loose.** The depth-separation scale concerns the size of the initialisation needed for incremental learning, not a margin at an entry.

**M8. Self-containedness, motivation and form.**

- **The companion paper.** [comp16] is a placeholder, yet Proposition 3.3 and §1.1 rely on it. Either include the branch expansion and the mirror-flow interpretation, or make the companion available to the referees.
- **Motivation.** The law (2.1)–(2.2) mixes a gradient step in α with a multiplicative step in β and a mass-weighted drift. Nothing in the paper explains why it is natural (adaptive control? neural networks?) or what it is for.
- **Length and figures.** At 33 pp. the paper is long. I would move Remark 3.6, Remark 4.3, Estimate 4.5/Remark 4.6, §7.2's paragraphs on moving widths and on WFR, and Proposition 6.4 to a supplement. Alternatively, split the paper: (I) Theorems 3.2 and 5.3, rigorous; (II) finite mass and recovery predictors. The main phenomenon has no figure (centres and duals against τ, the walls, a layer). Figure 1 is a single convergence plot.

## Minor comments (numbered, with line or page numbers)

1. **l. 41 (title) vs. the task's title.** The manuscript's title differs from "The rich transit of self-adaptive kernels". See M6 on overclaiming.
2. **l. 164 (Lemma 2.1(iv), p. 6).** "η_i − α_i/(2Gς) → sgn(α_i)λ_{1,i}" is a limit whose right-hand side diverges. Correct: η − α/(2Gς) = sgn(α)λ₁ − a coth(α/a)/(2Gς) + o(1).
3. **l. 241/245 (Proposition 3.1, p. 8).**
   - "Gradient flow of u³" should read: the gradient flow of the fit as a function of u through w = u³, in time scaled by h.
   - The time-scale claims (logarithmic, ε^{-p}log(1/ε)) are heuristics, but Table 1 calls them "exact expansions".
4. **l. 263–265 (Theorem 3.2, p. 9).** Spell out the energy identity of the limit flow and the fact D_ε(0) = D_∞(0); see M7. The assumption that Φ has full row rank (l. 144) is needed only to reach basis pursuit, not for the convergence.
5. **l. 277 (Proposition 3.3, p. 10).** "The two dual potentials combine by infimal convolution". For a pair sharing one column, the *dual* potentials add, and the sum of two one-sided indicators is the two-sided box. The *primal* regularisers combine by infimal convolution. Cite BPV26 Theorem 4.2.
6. **l. 291 vs. l. 308.** "Frozen while ε ≪ μ³" (centre) against "creeps into it unless ε ≪ μ⁴" (dual reaching the wall, Proposition 4.4). Say that these are different criteria.
7. **l. 305 (Formal statement 3.5).** "At every ε" holds by symmetry only if the samples are also mirror-symmetric (cf. D3ao-7). State this.
8. **l. 385 (Proposition 4.2(b)).** "After the slide onto the feature the split is frozen" is asymptotic: the kernels coincide only as s → ∞.
9. **l. 418 (Estimate 4.5).** Mark the near-tie heuristic visibly, and justify the "four kernels left nearest their walls" or call it a parameter.
10. **l. 435 vs. l. 514, and Table 2 ll. 228–231.** A_k denotes both the active sets and the amplitudes, and Table 2 lists both. Further overloads:
    - κ: κ_i′ (l. 158), the Gibbs slope (l. 331), κ_ij (l. 418), κ(μ) (l. 476);
    - s: the inner time, the signs s_i, s_k, s_A, the sample index (l. 144), s_θ;
    - X: l. 295, Step 2, Lemma 5.2(iii);
    - ℓ: ll. 136, 269, 357, 583, 748;
    - t: time, tanh, tail;
    - R: the regulariser, the largest distance;
    - r: the residual, the runner-up, the exponent in §7.5 B.
    The README says the overloading was fixed; it largely was not.
11. **l. 476.** κ(μ) should carry b_j.
12. **l. 491.** "The O(ε^{3/4}) is what the matching proves": the term O(ε^{1/2}X) in Step 3 is not traced to a specific perturbation. Name its source.
13. **l. 512.** "In one dimension the loop follows it by Theorem 5.3". The statement and proof of Theorem 5.3 are dimension-free. Either restrict Theorem 5.3 or drop "in one dimension" (README inconsistency 3).
14. **l. 553/561 (Proposition 6.3(d)).** "With the projections dropped (Π_k = 0) … one feature fitted per phase" is contradictory: Π_k = 0 gives r_k = f*. What is meant is r_k = f* − Σ_{fitted}A_mφ(c*_m), with the cross-talk neglected.
15. **l. 593 (Remark 6.6).** ξ_θ = (k₃/θ)^{1/3}/ρ does not depend on |A_k|. "ξ_θ ∝ (θ/|A_k|)^{-1/3}" needs a normalisation of the amplitudes; state it.
16. **l. 597 (Remark 6.7).** The midpoint computation assumes a starting point at the midpoint; say so.
17. **l. 656 (§7.2, p. 26).**
    - The entry times 16.18, 29.59, 33.05 against 33.496 are those of the 0.15-apart pair only. For the 0.25 pair the md gives 12.06, 16.05, 16.95 against 17.065.
    - The "staircase" path for the wide bumps did not terminate: "not terminated: 21 kernels active, 13 above the threshold" (md). Say so.
18. **l. 658.** "Consistent with the correction O(ε^{1/2}/μ)": replace this per M4.
19. **l. 663.** Three of the 54 inner problems near x ≈ 34 were still drifting (md, Part A); mention it.
20. **l. 675 (§7.3).** The bisection to a factor 1.09 quantises the ratios: 0.917, 0.771, 0.648 and 0.545 are 1.09^{-1}, ^{-3}, ^{-5}, ^{-7}. Report the medians with this resolution.
21. **l. 692–700.** Note that the 8 Part-A targets have margins 0.063–0.151 and that the md's window start uses C = 30; see M2.
22. **l. 720.** The polish-rejection rule is a heuristic of the path's computation and affects what "the path" is. Discuss its effect, and the 14 rejected polishes on 3 paths (§7.5 map).
23. **l. 763–776.** The files named in the code and data statement all exist; I checked. The DOI and the [repo] entry are placeholders.
24. **Bibliography (ll. 779–862).**
    - From memory, the volumes and pages flagged % CHECK look right to me: BMBO13 Math. Comp. 82, 269–299; CFG14; DP15; DDPS20; BSR17; Chi22; PKP23; FT07; Yin10; YOGD08; KW97; KS01; Fen79; Mos69; BP72; RW98 Thm 7.17.
    - The AGF venue (NeurIPS 2025) is unverified.
    - Eckhaus/Kaplun are reasonable guesses.
    - All of these must still be checked against the sources before submission.

## Verification I did

**[own] derivations**
- Lemma 2.1(i)–(iv), including the conservation law, dη/dα = 1/(Gβa′) and λ₁ ≈ Λρ.
- Proposition 3.1's expansions and k³.
- Theorem 3.2's chain: D_ε convex; ∇D_ε = Φw − y; g_ε = O(ε^{3/2}δ^{-2}) inside; exponential growth beyond; Theorem 7.17 of Rockafellar–Wets; the minimal section is Definition 2.2; the energy argument is valid because D_ε(0) = 0.
- Proposition 4.1: both integrals; the inversion z = t(1 − 2t²/3 + 13t⁴/15); 𝒲 = …(1 − 2t² + 4t⁴), giving −8/(κξ)² and +64/(κξ)⁴; ℓ = 2/κ; (2ςa/m)ℓ³ = k₃ρ^{-3}; the exponential branch e^{2−κξ}; temperature a/(4Gς).
- Proposition 3.4's creep prefactor πa²εb^{3/2}/(4ς), and cos/sin³ ≤ (2/(πμ))³.
- Remark 3.6's prefactor.
- Lemma 5.2(i)–(ii).
- Proposition 6.2.
- Proposition 4.4's two thresholds.
- Remark 4.6: μ = 1 − e^{−h(h−2|δ|)/(4σ²)} and P(μ < y) ≈ 4σ²y/h².
- Remark 6.7's midpoint inequality.
- The constants: κ = 0.7025, k₃ = 0.0962, Λ₀ = 7.0248, 1/κ = 1.4235; the Gibbs value 6.541 at θ = 0.01·0.8.
- σ_k from first principles (M1).
- The identity "inner problem = fast flow with t = Λ₀s" (M3).
- The time-origin error in Step 2 and the O(ε^{1/4}) of Step 5 (M2).

**[own] numerics** (scripts in W/own)
- prop34.py: the exact core coefficient against Proposition 3.4. The relative error times μ²/(εb) → −16.21 = −4/(π²ς), for b = 1 and 2.
- twokernel.py: 20 000 Newton starts for the two-kernel critical points of one Gaussian (1-D continuum, w > 0). None besides the exact fits.
- Recompiled the paper: 33 pp., 0 warnings, text identical to the shipped PDF.
- Regenerated mirror_transit.md, transit_mass.md, transit_limit.md and transit_certificate.md with `--from-raw`: all byte-identical.

**[md]**
- Table 3: rows 10⁻², 10⁻⁴, 10⁻¹².
- Table 4: rows 10⁻³ to 10⁻⁵.
- Path entry times 10.72/12.85/16.39.
- Creep ratios and margins (and their fit to 4/(π²ς), M4).
- The 0.15/0.25 pairs and the staircase.
- Table 5's "entry" row; the medians 0.917/0.771/0.545; ε_c = 4.88e-9; the slopes 0.294/0.462/0.285.
- The thresholds 6.404/6.300/6.272/6.541; the shares 0.48/0.437/0.309 against 0.481/0.439/0.316; 1.5448; 0.8849.
- Table 6: all six rows over the eight targets, and the slope ranges 1.501–1.504, 1.001–1.002, 1.233–1.533.
- Part B: v_k, slides, |s_θ+σ|, pulls, λ_k, (H3) ranges.
- Table 7: rows A, D(a), D(b) and the totals.
- 25 of 333; 4719 and 3128; 4.67σ.

All sampled numbers match.

**[snippet]**
- experiments/exp_transit_limit.py ll. 105–133: σ_k computed with the entrant's sign; v absolute; dev updated as s_j·Δ + σ(a₀ − a₁).
- transit_certificate.md l. 71: D(c) is selected by |m| > 0.8q.

**[read]**
- README (all six inconsistencies; inconsistencies 1 and 4 resolved above).
- review/REVIEW_LOG.md round 41 (D3an-1…13, in particular D3an-10 "v_k with its ρ" and D3an-11/12 on (b)'s window).
- BPV26 (arXiv 2606.23198), full text of §§3–4: Theorem 3.8 with its weighted convergence; Theorems 4.1–4.2.
- arXiv abstracts of PF23, Ber23, AGF (2506.06489), Chizat (1907.10300), ADCG, CMR23, LMRVV24, Elvira et al. (1904.06311), Schiebinger–Robeva–Recht (1506.03144), Jacot et al. (2106.15933), Traonmilin–Aujol (1811.12000).
- Not read in full: PF23, Ber23, AGF, Chizat, BMBO13, DDPS20, DP15.

**[memory]**
- The bibliographic details in item 24.
- RW98 Theorem 7.17's content (equivalence of epi-convergence and uniform convergence on compacts that avoid the boundary of the domain).
- Elvira et al. in ACHA (2021); Schiebinger et al. in Inf. Inference (2018).

## Literature

- **ISS / linearised Bregman** [BGOX06, BMBO13, YOGD08, ORXYY16]. Correctly attributed. The contrast with the LBISS (no tails, no creep) is apt.
- **Saddle-to-saddle at small initialisation.**
  - PF23 (depth 2, log scale, LARS-like recursion) and Ber23 (incremental learning, restricted features) are described correctly per their abstracts.
  - BPV26 is the closest for Theorem 3.2. Its Theorem 4.2 already gives the weighted-ℓ¹ ISS for asymmetric diagonal networks, which is Proposition 3.3's pair case.
  - The paper's unweighted convergence improves on BPV's s ds weighting because the start is interior; say so.
  - Add Jacot et al. 2021.
- **AGF** [KMCKSDGM25]. "Dormant maximise a utility, active minimise a cost" is the right precedent for Theorem 5.3(b). Also link the dormant alignment to the WFR variant (Remark 3.6).
- **Sliding Frank–Wolfe / ADCG** [DDPS20, BSR17]. Fairly described as discrete-time analogues. A numerical comparison of outcomes (SFW from the same targets) would be informative.
- **Moving particles.** Chizat 2022 must be compared on its main result, not only for its geometry. The description "conic particle gradient descent in a WFR geometry" is consistent with its abstract (local analysis in Wasserstein space, perturbed mirror descent), so the % CHECK on it can be cleared.
- **Dual certificates** [CFG14, dCG12, DP15, PKP23]. The paper's p(c_j) is not DP15's η_V. Compute η_V and the BLASSO outcome as baselines (M6). Add:
  - Schiebinger–Robeva–Recht, against Remark 6.7;
  - Elvira–Gribonval–Soussen–Herzet, the greedy analogue of Conjecture 6.5;
  - Traonmilin–Aujol, for the basin hypothesis.
- **Novelty claims (i)–(iv), l. 122ff.**
  - (ii) and (iii) appear new.
  - (i) is new for this unit, but the mechanism is BPV's.
  - (iv): the fitted skeleton's success is largely its own definition (M6).

## Issues

| id | sev. | issue | proposed fix and where |
|---|---|---|---|
| D3ei-1 | major | σ_k lacks 1/ρ_{j_k} (v_k is relative, dev absolute); the md header lacks s_k; tested with equal masses only | Thm 5.3(c), l. 486: σ_k = −s_k dev_{j_k}/(ρ_{j_k}v_k) = −s_k dev/\|a_{k−1,j_k}\|; transit_limit.md header; one unequal-mass test |
| D3ei-2 | major | Thm 5.3(b)'s left end is not covered: Step 2's time origin, Step 5's O(ε^{1/4}); md window start at C = 30 < 3/(2λ_k) on 7 of 8 targets; observed ε^{0.8–0.9} there | ll. 485, 500, 506: correct origin; window start + 5/(4λ_k)·ε^{1/2}log(1/ε) or a stable-direction O(ε^{1/2}) lemma; two-term bound; Table 6 window-start row; rerun with C ≥ 3/(2λ_k) |
| D3ei-3 | major | The inner problem is the path's relaxation with t = Λ₀s; "another metric" is false; (H3) is the convergence presupposed by Def. 2.3; Lemma 5.2(iii) does not control (−∞, s₁] and assumes the continuum implicitly | ll. 205, 445–452, 464, 472: state the identity, rephrase (H3), fix the proof of (iii) |
| D3ei-4 | major | Prop. 3.4's stated O(ε^{1/2}/μ) is not sharp; the proof gives O(ε/μ²) with the constant −4b/(π²ς); §7.2 says the data are consistent with the wrong order | ll. 285, 299, 658: state the sharp term and ε^{1/2} ≪ μ; cite Prop. 4.1's −8/(κξ)² |
| D3ei-5 | major | "1154/1154" is the basin hypothesis restated; polish-off check on 640 targets only; D(c) selected by \|ν_p\| > 0.8q, not stated; no η_V/BLASSO baseline; title overclaims | abstract, title, §6, §7.5 ll. 720–727: reword; full polish-off check; disclose the selection; add η_V/BLASSO and continuous OMP |
| D3ei-6 | major | Literature gaps: Chizat's main result, Elvira et al., Schiebinger–Robeva–Recht, Jacot et al., Traonmilin–Aujol; BPV Thm 4.2 = Prop. 3.3's pair case; weighted vs unweighted convergence | §1.3, Prop. 3.3, Thm 3.2 proof |
| D3ei-7 | major | Not self-contained ([comp16]); no motivation for the law; 33 pp.; no figure of the transit | include the companion's needed parts or attach it; motivation paragraph; supplement or split; schematic figure |
| D3ei-8 | minor | Prop. 4.2: (b)'s classification of critical points unproved (numerically true in 1-D); (c)'s Gibbs split asymptotic, not exact; (d)'s I_∞ lacks s_A, its positivity numerical | ll. 372–387, Table 1 l. 95 |
| D3ei-9 | minor | Prop. 3.3: "infimal convolution" of dual potentials; they add | l. 277 |
| D3ei-10 | minor | Prop. 6.3(d): "Π_k = 0" contradicts the intended residual | ll. 553, 561 |
| D3ei-11 | minor | Lemma 2.1(iv) states a limit with a diverging right-hand side | l. 164 |
| D3ei-12 | minor | Notation overloads (A_k twice in Table 2; κ, s, X, ℓ, t, R, r); κ(μ) without b_j | Table 2, l. 476, passim |
| D3ei-13 | minor | Thm 5.3 is dimension-free, but §6 says "in one dimension" | l. 512 |
| D3ei-14 | minor | §7.2 conflates the 0.15 and 0.25 pairs; the staircase path did not terminate | l. 656 |
| D3ei-15 | minor | Remark 6.6: dimensions of ξ_θ; Remark 6.7: midpoint starting point assumed | ll. 593, 597 |
| D3ei-16 | minor | Bisection quantisation of the medians; 3 drifting inner problems not reported | ll. 663, 675 |
| D3ei-17 | minor | Thm 3.2: spell out D_ε(0) = D_∞(0) and the limit energy identity | l. 265 |
| D3ei-18 | minor | μ³ (frozen centre) vs μ⁴ (creep into the wall) not reconciled; FS 3.5 "at every ε" needs symmetric samples | ll. 291, 305, 308 |
| D3ei-19 | minor | Prop. 3.1: "gradient flow of u³"; its scale claims are heuristic, but Table 1 says "exact expansions" | ll. 88, 241, 245 |
| D3ei-20 | minor | Step 3's O(ε^{1/2}X) term not traced to a source; the remainder budget is not sharp (O(ε) observed) | l. 502 |
| D3ei-21 | minor | GSD20 analogy (margin ε^{1/2} as "the polynomial scale") is loose | l. 111 |
| D3ei-22 | cosmetic | 30 bibliographic % CHECK items; Chizat's description can be cleared; AGF venue unverified | bibliography |

## Lesson

- When a quantity is normalised (v_k divided by ρ), every formula that divides by it must carry the normalisation back. A test run only at ρ = 1 cannot catch the omission.
- Before calling two flows "the same gradient flow in another metric", write both in common coordinates. Here they are identical, and the hypothesis built on the difference, (H3), turns out to be the convergence that the definition already assumes.
- A creep ratio printed to four digits, compared with the stated μ-dependence, would have exposed Proposition 3.4's wrong error order at once. Fitting the residual's exponent is cheaper than arguing about it.
