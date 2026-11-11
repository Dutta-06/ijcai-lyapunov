# Referee report D3ej, round 141: the Paper 16 manuscript, first review

## Scope
I reviewed `manuscript/paper16/paper16.tex` (724 lines; PDF of 28 pages), `supplement.tex`, `README.md`, `make_figs.py` and `numcheck16.py` at commit a833d10. I exported the commit to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref141/repo` and did not touch the live repository.

I also read the sources: PARADIGM.md ll. 156–158 and 5499–5999 (D‴.45, D‴.46, D‴.47) and `review/drafts/D44.md` (for Cor. 5.3). I read the four md reports: `mirror_free`, `mirror_branch`, `mirror_leak` and `mass_dictionary`.

Line numbers below refer to `paper16.tex`. Page numbers refer to `paper16.pdf`.

## Summary of the submission
In the law under study, each kernel's coefficient is the product of two factors:
- a mass β, adapted multiplicatively by an exponentiated-gradient law;
- a saturating amplitude a = a_max·tanh(α/a_max), adapted by a gradient law.

The masses are free. One scalar drive moves both factors, so Q = β − ς·cosh(2α/a_max) is conserved, and w = βa/m follows a mirror flow. The potential is given in closed form:
- an ℓ² core weighted by the starting mass;
- a soft wall at the saturation and an entropy tail;
- an ℓ¹ limit with weight λ₁ ∝ (β⁰)^(−1/2), plus a −(3/2)K|w|^(2/3) correction;
- on "branches" (Q < −ς), an entropy with temperature 1/(γa₀) and the sign of w frozen.

The main results are:
- **Theorem 4.1.** From every start the dictionary loop converges to the Bregman projection of the start onto the sign-constrained least-squares fits. Wrong-signed kernels go extinct at the rate γa₀|g(∞)|.
- **Rich and lazy limits.** The lazy limit is the mass-weighted minimum-norm fit (Prop. 4.3(a)). The rich limit is weighted basis pursuit, approached like ε^(3/2) (Prop. 4.3(b), Derived st. 4.4). On branches the rich limit is sign-constrained basis pursuit, approached like ε^μ (Prop. 4.6).
- **Pairs and positive dictionaries.** A pair of particles per kernel gives the hyperbolic entropy (Prop. 4.8). On positive dictionaries the loop selects the target (Prop. 4.9).
- **Noise and sampled data.** Smooth (Stratonovich) noise does not move the limit map (Prop. 5.1). An Euler step changes Q by a closed-form O(h²) amount, which an exact step along the curve removes (Prop. 5.2). The reduced flow under label noise ends at the split Gβ(1 − t⁴) = γa² (Cor. 5.3).
- **Adaptive loop.** ½e² plus the Bregman divergences is a Lyapunov function (Prop. 6.1). A leak in the mirror coordinates keeps Q, gives ultimate boundedness (Prop. 6.2) and a unique attracting averaged equilibrium (Prop. 6.3). That equilibrium is a ridge from large starting masses and a weighted LASSO from small ones (Prop. 6.4).

## Overall assessment
**Recommendation: major revision.**

**What holds up.** The paper is carefully written and honest about status (Table 1, the "Derived statement" label). Its core identities are correct; I re-derived each of them and checked the main ones numerically. The best parts are:
- Theorem 4.1's argument for every start: the range argument identifying the limit, and the extinction rate. I confirmed the rate to 7 digits on my own example.
- The closed-form potential.
- The ε^μ versus ε^(3/2) contrast.
- The point that the choice of sampled-data integrator decides whether measurement noise selects a fit.

**What needs fixing.**
1. Proposition 4.9, as stated and proved, claims more than the hypotheses give. The certificate does not imply that the target is the only nonnegative fit.
2. The statement and discussion of the Euler change of Q (Prop. 5.2(a), l. 457, and the text at l. 470) misstate its sign. By the paper's own formula, the drift is strictly positive on the critical curve and on every branch. So, to leading order in h, a free particle cannot be pushed onto a branch by the Euler drift, which is the opposite of what the text says. The report's own reduced-flow runs agree with this (no dead atoms; min Q at −0.025 = −ς, its printed floor).
3. The 1745/1760 support check is interpreted wrongly. With the report's own support threshold, the LASSO itself matches the "prediction" on 1739 of 1760 patterns, and its support equals the truth's for 437 patterns, against 415 for the equilibrium. The support statistic therefore does not "check the limit of Proposition 6.4(b)".
4. The literature comparison misses several directly relevant works: label-noise SGD solves the Lasso for the u⊙v parametrisation; finite steps break conservation laws in closed form; the step size biases diagonal nets; NNLS gives sparse recovery on positive designs; and ℓ² penalties on a Hadamard product give the LASSO.

**Novelty and fit.** Much of Section 3 instantiates general theory (Li–Wang–Lee–Arora; Marcotte–Gribonval–Peyré) for one specific unit. The adaptive-control results (Prop. 6.1, the leak) are Boffi–Slotine and Gao–Pavel/σ-modification transplanted to that unit. The paper says so fairly. Even so, a SICON or TAC reader will ask why this particular law deserves a 25-page treatment, and the motivation (§1.1, l. 66) is thin on that point. At 28 pages, the paper is too long for a TAC regular paper. SICON is the better fit after shortening §1.3 and §7.

## Major comments

**M1. Proposition 4.9 (l. 425–434, p. 16): the certificate does not give uniqueness.**
- The proof (l. 433) shows that every nonnegative fit vanishes off supp w*. It then concludes "w* is the only nonnegative fit". That step also needs Φ restricted to supp w* to have full column rank (injective).
- Counterexample [own, analytic]: φ₁ = φ₂ = (1, 1), φ₃ = (0.5, 2), w* = (½, ½, 0), y = (1, 1).
  - Φᵀ(1, 1) = (2, 2, 2.5) > 0, so the row-space hypothesis holds.
  - v = (−1, 1) gives h = (0, 0, 1.5), so the certificate holds with μ = 1.5.
  - Yet every point (s, 1 − s, 0) with s ∈ [0, 1] is a nonnegative fit. So S is a segment, and the loop converges to the start-dependent Bregman projection on it, not to w*.
- The source D‴.47(v) (PARADIGM.md l. 5886–5888) assumed uniqueness ("if the target w* ≥ 0 is its only point — certified by …"); the paper dropped that assumption.
- The off-support bound and the convergence are unaffected.
- **Fix:** add "Φ_{supp w*} has full column rank" (it holds for Gaussian translates at distinct centres, the cited case), or state the result as convergence to the Bregman projection onto the face. Then correct the abstract (l. 51: "the frozen signs make the target the limit"), §1.2 (l. 70), §1.3 (l. 118) and Table 1 (l. 95).

**M2. The sign of the Euler change of Q (Prop. 5.2(a), l. 457; text l. 470; also l. 485).**
- Let b = γa² − Gβa′(1 + t²) be the bracket in Prop. 5.2(a).
  - On the critical curve, β = 2ς·sinh²(α/a_max), so b = 2Gς·t²(1 − t²) > 0.
  - On a branch (Q + ς = −b₀ < 0), b = 2Gς·t²(1 − t²) + Gb₀(1 − t⁴) > 0.
- Just above the critical curve (β⁰ ≪ ς), b < 0 only inside the core, |α| ≲ a_max(β⁰/(2ς))^(1/2). That is far below the saturation, so "negative below the saturation's layer" is wrong for small masses.
- Consequences:
  - (i) The Euler drift raises Q on every branch, which tends to unfreeze signs.
  - (ii) Near the corner, |ΔQ| ~ h²g²γGβ⁰², which is much smaller than the gap β⁰ to the critical curve. So Q → −ς only asymptotically.
  - (iii) The sentence at l. 470, "a particle whose Q_i reaches −ς moves onto a branch, where its sign freezes and its mass can die", is not supported, and to leading order in h it is false.
- Evidence:
  - [own] `own1.py`, check (7): b > 0 on Q = −ς, −1.5ς and −4ς over 400 points each; the exact Euler ΔQ is > 0 to rounding.
  - [own] `own2.py`: a confined particle under noisy Euler from β⁰ = 10⁻² and 10⁻³ at h up to 0.5 never crosses Q = −ς, and Q + ς only decreases toward 0. A particle started on a branch has its Q rise.
  - [md] `mirror_free.md` Part C(3), rows at τ = 10⁸: "dead atoms (Q ≤ −c)" is 0 on every instance, and the minimum of Q is printed as −0.025.
- Note: Itô read in (α, β) coordinates instead of (α, log β) would make the drift negative everywhere; see m4.
- **Fix:**
  - State the sign region exactly: b < 0 iff Gβ(1 − t⁴) > γa², which is never the case on the critical curve or on a branch.
  - Replace l. 470 by the correct mechanism: Q decreases only in the core, the masses approach the corner (Q, α) → (−ς, 0) from above, and branches drift back.
  - At l. 608, write "its minimum approaching −ς" instead of "reaching −ς".
  - Check D‴.45(iv)(b)'s "atoms die (Q → −c)" against this.

**M3. The 1745/1760 support study (l. 619; `mirror_leak.md` D(2)).**
- In the code (`exp_mirror_leak.py` ll. 258–276), the support is {|w| > 10⁻³·max|w|}. The equilibrium's support is compared with the irrepresentable-condition prediction only, never with the LASSO's own solution.
- [own] `own4.py`, using my coordinate-descent LASSO on the same orbit Gram, the same rng seed (23) and penalty 10⁻³:
  - The exact-zero LASSO support equals the truth for 416 patterns and agrees with the prediction on all 1760. The theory is exact here.
  - With the report's threshold, the LASSO's support equals the truth for 437 patterns, and agrees with the prediction on 1739. All 21 mismatches have IR quantity between 1.0009 and 1.0134; the LASSO's 4th coefficient falls below the threshold.
  - The equilibrium (β⁰ = 10⁻⁶) is "exact" on 415 patterns. So at the support level the equilibrium and the LASSO differ on at least 22 patterns.
- The sentence "This checks the limit of Proposition 6.4(b), not the LASSO's theory" therefore inverts what was measured. The statistic mixes a threshold artefact with finite-β⁰ effects; the −(3/2)K|w|^(2/3) under-penalisation is a plausible cause.
- **Fix:** report, per pattern, the distance between the equilibrium and the LASSO solution, and the exact-support agreement. Explain the 15 mismatches, or drop the count.

**M4. Literature and novelty: missing comparisons that bear on specific claims.**
- Pillaud-Vivien, Reygner, Flammarion, COLT 2022 (PMLR 178:2127–2159) [snippet]: label-noise SGD on the quadratic (u⊙v) parametrisation implicitly solves the Lasso, with support-recovery conditions. This is the closest prior result to Cor. 5.3 and to §7.3's "reweighting, not basis pursuit". It must be compared.
- Kunin, Sagastuy-Brena, Ganguli, Yamins, Tanaka, ICLR 2021 (arXiv:2012.04728) [snippet]: finite learning rates break symmetry-induced conservation laws, with exact modified-flow expressions. l. 120 attributes "finite steps break the conservation laws" to MGP23. The novelty "closed-form change of Q per step" (l. 120) should be positioned against Kunin et al. and Even et al.
- Nacson, Ravichandran, Srebro, Soudry, ICML 2022 (PMLR 162:16270–16295) [snippet]: the step size alone moves diagonal nets from the kernel to the rich regime. This is relevant to "which update the loop runs decides" (l. 470).
- Slawski–Hein, EJS 7 (2013) 3004–3056 [snippet], and Meinshausen, EJS 2013 [memory]: NNLS gives sparse recovery without regularisation on self-regularising (positive) designs. These are relevant to Prop. 4.9 and Theorem 4.1's sign-constrained fits.
- Hoff, CSDA 115 (2017) 186–198 [snippet]: an ℓ² penalty on a Hadamard (u∘v) parametrisation equals the LASSO. This is the explicit-regularisation analogue of Prop. 6.4(b)'s "leak gives a LASSO".
- Zou, JASA 101 (2006) 1418–1429 [snippet]: weighted ℓ¹ (the adaptive lasso) and its consistency conditions, for the weighted irrepresentable condition at l. 561.
- Bregman-divergence adaptive laws predating Boffi–Slotine: Lee et al. 2018 on natural adaptive robot control, and Wensing–Slotine 2018 [snippet, secondary]. Also Du–Hu–Lee 2018 on balancedness as the conservation law of u⊙v [memory].
- The comparisons that are present (LWLA22, MGP23, Azulay, Woodworth, Ghai–Hazan–Singer, Chou et al., Vaškevičius, Pesme, Even, Boffi–Slotine, Gao–Pavel, Azizan–Hassibi, σ-modification, Zhao–Yu/Wainwright) are fair and accurate as far as I checked.

**M5. Corollary 5.3 and the reduced flow's validity region (l. 472–485; Table 1 l. 98).**
- Cor. 5.3 rests on `review/drafts/D44.md`, whose re-review is still pending (README).
- Katzenberger's theorem applies on compact parts of an attracting manifold with nondegenerate normal rates. Here the normal rates are Φ·diag(mD)·Φᵀ with D ∝ β², and they degenerate as β → 0, which is exactly where §5's narrative sends the flow ("masses fall toward extinction").
- The corollary's necessary condition at rest points with β_i > 0 is fine. The paper should say that statements about extinction and about which sparse fit is reached lie outside the reduced-flow theorem's range, not only outside the proofs.

## Minor comments
1. **l. 457, l. 477:** write γa_max³|t_i|³, not t_i³. For t < 0 the right-hand side as written is negative.
2. **l. 442, l. 449 (Prop. 5.1):** "Read in Itô's sense" depends on the chart. Itô in (α, log β) gives Prop. 5.2(b)'s drift. Itô in (α, β) drops the +γ²βa²/2 term, which leaves the drift −(γGβ²/2)(1 − t⁴)·(…) < 0. Say which chart is meant (the Euler update's), and add one sentence on Euler in (α, β). The Stratonovich statement itself is correct: the chain rule holds and the noise is additive in η.
3. **l. 331 (Thm 4.1(c) proof):** the step "log|w_i| = s_iγa₀η_i + O(1)" uses Prop. 3.4(b), which holds for Q < −ς only. On the critical curve, |w| ~ t⁻³ and β ~ t⁻², so −t⁻¹·log β → 0. That agrees with the formula at a₀ = 0, but it needs its own line.
4. **l. 323:** boundedness follows from coercivity of B(w̄, ·), that is, from superlinearity of ψ* (exponential in η), not of ψ. Reword.
5. **l. 366, l. 372 (Derived st. 4.4):** the log-derivative of cos(πμ/2)/sin³(πμ/2) is −(π/2)tan(πμ/2) − (3π/2)cot(πμ/2). The relative error is therefore O(ε^(1/2)/(μ_j(1 − μ_j))); it blows up as μ_j → 1, where the leading term vanishes.
6. **l. 116 ("ε^(3/2) whatever the slack") and l. 399 ("reached far more slowly"):** the constant is μ_min⁻³ and the regime needs ε^(1/2) ≪ μ_min. Also, "slowly" here refers to ε, not time. Reword both.
7. **l. 284 (Prop. 3.4(e)):** the O(β(0)) in β⁰ = (1 − θ²)β(0)(1 + O(β(0))) is not uniform as θ ↑ 1; the relative error there is O(β(0)/(1 − θ²)). [own] At β(0) = 10⁻⁴ and θ = 0.95 the ratio is 0.9944, against 0.99999 at θ = 0.3. State that θ is fixed.
8. **l. 255 (Prop. 3.3(c)):** "the marginal cost rises by the factor β⁰/(4ς)" compares a value of η with a slope dη/dlog|w|. Say so.
9. **l. 560 (Prop. 6.4(a)):** a fixed ρ = k_eσm/(G(β⁰)²) needs σ ∝ (β⁰)². In `mirror_leak.md` A, σ ≈ 2·10⁶ at β⁰ = 10⁴. State the scaling, and note that averaging requires small adaptation gains, i.e. mD ≈ G(β⁰)²/m small.
10. **l. 561:** the limit object "σ·Σλ_{1,i}|w_i|" contains ε-dependent λ_{1,i} → ∞. Write it as σΛ(ε)·Σr_i^(−1/2)|w_i| with σΛ fixed.
11. **l. 528 (Prop. 6.2):** the drift example is stronger than stated. With d = ė + k_e·e + wφ, all three terms tend to 0, so d(t) → 0. Also, the ultimate bound is given only for n = 0. Because B is superlinear, the leak also absorbs a bounded measurement noise: n(w − w*)·φ ≤ (σ/2)B + C(n). Consider stating that.
12. **l. 552:** ψ = L[φ] clashes with the potential ψ. Rename.
13. **l. 608:** "its minimum reaching −ς" should read "approaching −ς (to the printed digits)". See M2.
14. **l. 608, Cor. 5.3 check:** the ratio Gβ(1 − t⁴)/(γa²) = 1.0000 confirms that the integrator reached a rest point. It is not independent evidence for the corollary; say so, as is done for the 2.22·10⁻¹⁶.
15. **l. 436 and bibliography l. 708:** SRR's journal version is Inf. Inference 7(1):1–30 (2018) [snippet], so that % CHECK can be closed. Their Gaussian theorem needs |S| > 2M samples and no separation [snippet], consistent with the paper.
16. **Figure 3 (left, p. 22):** the dotted reference line extends far beyond the data, to log₁₀(1/ε) ≈ 128. Clip it.
17. **Length:** 28 pages, with §1.3 about 3 pages. For TAC the paper must roughly halve. For SICON, shorten §1.3 and move §7.4's second paragraph to the supplement.
18. **§1.1 (l. 66):** the "measure on the kernels' parameters" motivation for the two-factor unit is asserted in one sentence. Give a reference or a short derivation of why a controls reader would use this unit rather than a plain gradient law with projection.

## Verification I did
- **[own] `own1.py`** (my code, not the repository's):
  - η's closed form (Prop. 3.2) against quadrature on 10 levels: Q = 100ς, 3, 0.5, ς, ς/2, 0, −ς/2, −ς (coth reading), −2ς (branch), −10ς. Largest relative error 4.95·10⁻¹⁵.
  - (ψ*)′ = w and dw/dη = mD to 1.8·10⁻¹⁰ and 2.5·10⁻¹⁰.
  - B(0, w) = β/(γm) on a branch to 12 digits.
  - The temperature 1/(γa₀) at the branch point (15.583885 vs 15.583874 at a distance 10⁻⁷ from the branch point).
  - K = 3.218298; η_c·w^(1/3) → −3.218298 at α = 10⁻⁴.
  - η − λ₁ → η_c (Prop. 3.3(d)).
  - The arcsin tilt (Prop. 3.4(e)) to 10⁻⁵ at β(0) = 10⁻⁶.
  - Q conserved under a random multi-sine drive to 2.5·10⁻¹³.
  - The exact Euler ΔQ against Prop. 5.2(a): relative error O(h), e.g. 1.67·10⁻³ and 1.67·10⁻⁴ at h = 10⁻² and 10⁻³.
  - The sign of the bracket (M2).
- **[own] `own2.py`:** noisy Euler simulations near and on the critical curve (M2).
- **[own] `own3.py`:** Theorem 4.1 on a 3 × 7 dictionary, with the loop integrated in (α, log β) (LSODA) and the projection computed from my own dual, using an accurate branch parametrisation.
  - Mixed start (3 free particles, 4 branches), exact fit: the end at T = 4·10⁴ is the Bregman projection to 4.9·10⁻¹⁰, with Q kept to 3.3·10⁻¹¹.
  - All-positive branches with an infeasible sign pattern: the limit is the sign-constrained least-squares point to 2.8·10⁻¹⁶. The five Z particles have extinction rates equal to γa₀|g(∞)| to 7 digits (e.g. 6.481146·10⁻³ against 6.481146·10⁻³).
- **[own] `own4.py`:** the 1760-pattern LASSO study (M3). I imported only the orbit Gram from the experiment module; the LASSO is my own coordinate descent.
- **[own] Re-derivations by hand:** Lemma 3.1; Prop. 3.2 (all readings); Prop. 3.3(a), (b), (d) and K; Prop. 3.4(a)–(e); Theorem 4.1(a)–(c) including existence; Prop. 4.3; Derived st. 4.4's tangent profile; Prop. 4.6(b)'s ε^{μ_j}; Prop. 4.8(b); Prop. 4.9's bound; Prop. 5.2(a), (b) (dt/dτ, E g²); the flat-update split 2Gm|w|(1 − t⁴) = γa_max³t³ as the minimiser of D; Remark 5.4; Props. 6.1–6.3 (three-point identity, comparison bound, H_L + H_Lᵀ, the 1/k_e expansion, the bounded roots as σ → 0); Prop. 6.4's ρ.
- **[raw→md] Report reproduction:** all four md reports regenerated with `--from-raw` (env vars MF_OUT, MB_OUT, ML_OUT, MD_OUT pointed into W). All four are byte-identical to the committed reports (`cmp`).
- **[own run] `numcheck16.py`:** 290 distinct decimal tokens; 1 not found (10.4, the ratio derived in the text), as the README states.

## Literature
- **Present and fair:** LWLA22 and MGP23 (general mirror-flow and conservation theory); Azulay et al.; Woodworth et al.; Ghai–Hazan–Singer; Chou et al. (NNLS, and CMS23); Vaškevičius et al.; Pesme et al.; Even et al.; Boffi–Slotine; Gao–Pavel; Azizan–Hassibi; Coucheney et al.; σ-modification (Ioannou–Kokotović, Ioannou–Sun); Zhao–Yu and Wainwright.
- **Novelty claims:** mostly accurate and modest (l. 110: "Our contribution here is the explicit form").
- **Missing or under-compared** (see M4): Pillaud-Vivien–Reygner–Flammarion 2022; Kunin et al. 2021; Nacson et al. 2022; Slawski–Hein 2013 (and Meinshausen 2013); Hoff 2017; Zou 2006; Bregman-Lyapunov adaptive laws before Boffi–Slotine (Lee et al. 2018; Wensing–Slotine 2018); Du–Hu–Lee 2018.
- **Weakest claim:** the "closed-form change of Q" novelty (l. 120) is the most exposed to Kunin et al. and Even et al.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ej-1 | major | Prop. 4.9 concludes "w* is the only nonnegative fit" from the certificate alone. Injectivity of Φ on supp w* is missing. Counterexample: φ₁ = φ₂, w* = (½, ½, 0) [own]. The source assumed uniqueness [read: PARADIGM l. 5886–5888]. | Add the rank hypothesis, or state convergence to the projection onto the face. Fix l. 51, 70, 95, 118, 425–434. |
| D3ej-2 | major | Prop. 5.2(a) "negative below the saturation's layer" and l. 470 "reaches −ς moves onto a branch". The bracket is > 0 on the critical curve and every branch, and < 0 only in the core [own1 (7), own2]. The md shows 0 dead atoms [md mirror_free C(3)]. | Restate the sign region exactly and rewrite l. 470, l. 485, l. 608; check D‴.45(iv)(b). |
| D3ej-3 | major | 1745/1760 "checks the limit of Prop. 6.4(b)". With the report's threshold the LASSO itself has 437 truth-supports against the equilibrium's 415, and 1739 agreements with the prediction; the exact LASSO agrees on 1760 [own4; read exp_mirror_leak.py ll. 258–276]. | Compare equilibrium with LASSO per pattern (distance, exact support); explain the mismatches; l. 619. |
| D3ej-4 | major | Missing literature directly relevant to Cor. 5.3, Prop. 5.2, Prop. 4.9 and Prop. 6.4 (Pillaud-Vivien et al.; Kunin et al.; Nacson et al.; Slawski–Hein; Hoff; Zou) [snippet]. | Add the comparisons in §1.3 (l. 112–122) and temper the novelty claim at l. 120. |
| D3ej-5 | major | The reduced-flow narrative (extinction, which fit is reached) lies where the normal rates degenerate (D ∝ β²), outside Katzenberger's range; Cor. 5.3 rests on an unreviewed draft [read README, D44.md]. | State the validity region; keep Cor. 5.3 to rest points with β > 0; l. 472–485, Table 1 l. 98. |
| D3ej-6 | minor | Itô reading is chart-dependent: (α, β) and (α, log β) give different drifts [own derivation]. | Name the chart at l. 449 and l. 470. |
| D3ej-7 | minor | Theorem 4.1(c)'s proof uses Prop. 3.4(b) on the critical curve, where it does not hold; there β ~ t⁻² [own derivation]. | Add a line for Q = −ς at l. 331. |
| D3ej-8 | minor | Derived st. 4.4's relative error degenerates as μ_j → 1 [own derivation]. | Use O(ε^(1/2)/(μ_j(1 − μ_j))) at l. 366 and l. 372. |
| D3ej-9 | minor | "ε^(3/2) whatever the slack" and "reached far more slowly" (l. 116, l. 399) overstate and confuse ε with time. | Reword both. |
| D3ej-10 | minor | Prop. 6.4(a) silently needs σ ∝ (β⁰)²; Prop. 6.4(b)'s limit object contains a diverging λ₁ (l. 560–561). | State the scalings. |
| D3ej-11 | minor | Prop. 6.2: the drift example's d actually tends to 0; the bound covers only n = 0 (l. 528). | Strengthen the example; add the bounded-n extension. |
| D3ej-12 | minor | Prop. 3.4(e)'s O(β(0)) is not uniform as θ → 1 [own: 0.9944 at θ = 0.95, β(0) = 10⁻⁴]. | State that θ is fixed (l. 284). |
| D3ej-13 | cosmetic | t³ should be \|t\|³ (l. 457, l. 477); ψ reused for L[φ] (l. 552); "superlinearity of ψ" should be ψ* (l. 323); Prop. 3.3(c) compares a value with a slope (l. 255). | Edit in place. |
| D3ej-14 | cosmetic | Fig. 3 left: the dotted reference line extends far beyond the data; the SRR15 % CHECK can be closed (Inf. Inference 7(1), 2018) [snippet]. | `make_figs.py`; bibliography l. 708. |
| D3ej-15 | minor | Length (28 pp.) and venue fit; the motivation for the unit is thin (l. 66). | Shorten §1.3 and §7; motivate the unit; target SICON. |

## Sources
- **[read]** `manuscript/paper16/paper16.tex`, `supplement.tex`, `README.md`, `numcheck16.py`; `fig1.pdf` and `fig3.pdf` (rendered); PARADIGM.md ll. 156–158, 5499–5660 and 5808–5945; `review/drafts/D44.md` ll. 1–70; `experiments/exp_mirror_leak.py` ll. 58–70, 238–280 and 395–420; `experiments/exp_mirror_free.py` ll. 40–62 and 645–670.
- **[md]** `results/chaos/mirror_leak.md` (all); `results/chaos/mirror_free.md` (Part A; Part C(3) rows).
- **[raw]** the four `*_raw.npz` files, regenerated to md and compared with `cmp`.
- **[own]** scripts and outputs in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref141/`: `own1.py`, `own2.py`, `own3.py`, `own4.py`, `own4.out`, `own4b.out`. Report reproductions are in `out_mirror_free/`, `out_mirror_branch/`, `out_mirror_leak/` and `out_mass_dictionary/` there.
- **[snippet]**
  - Pillaud-Vivien–Reygner–Flammarion: https://proceedings.mlr.press/v178/vivien22a.html , https://arxiv.org/pdf/2206.09841
  - Kunin et al.: https://arxiv.org/pdf/2012.04728
  - Slawski–Hein: https://www.doi.org/10.1214/13-EJS868
  - Hoff: https://arxiv.org/pdf/1611.00040
  - Nacson et al.: https://proceedings.mlr.press/v162/nacson22a.html
  - Schiebinger–Robeva–Recht: https://arxiv.org/abs/1506.03144
  - Zou: https://experts.umn.edu/en/publications/the-adaptive-lasso-and-its-oracle-properties/
  - Lee et al. / Wensing–Slotine, secondary only: https://arxiv.org/pdf/1806.06655
- **[memory]** Meinshausen 2013 (EJS); Du–Hu–Lee 2018 (NeurIPS).

## Lesson
When a paper gives a closed-form drift, check its sign in every regime the text later talks about, not only where the formula was derived. Here, one line of algebra on the critical curve and on a branch (b = 2Gς·t²(1 − t²) > 0) reverses the qualitative claim, and the report's own "0 dead atoms" column had already shown this. Likewise, for a headline count, recompute the reference object itself (here, the LASSO's own support under the same threshold) before reading agreement with a theory as a check of a limit.
