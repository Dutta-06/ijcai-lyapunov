# Referee report D3ek, round 142: the Paper 15 manuscript, first review

## Scope
This review covers `manuscript/paper15/paper15.tex` (25 pp.), `supplement.tex` (3 pp.), the compiled PDFs, `README.md`, `make_figs.py` and `numcheck15.py` at commit e82680b. I exported the commit to `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref142/repo` and symlinked the venv; the repository itself was not touched.

I also read:
- the sources: PARADIGM.md D‴.41 (l. 5016–5247) and D‴.43 (l. 5350–5440, statement, proof and status);
- the reports `results/chaos/gauge_label.md` and `colour_selection.md`, and the parts of `exp_colour_selection.py` and `exp_gauge_label.py` that write them.

Evidence tags: [read] read in the files; [md]/[raw] from the reports or the raw npz files; [own] my own computation; [snippet] web search snippets only; [memory] my background knowledge, not checked.

## Summary of the submission
The setting is a gradient flow driven by a measured error, ż = ⟨(r+ξ)U⟩. The noise enters along the descent's own fields, which are g-orthogonal to the manifold Γ of exact fits (Prop. 3.1).

- **Continuous time (Theorem 3.3).** Invariance of the limit map makes B(a,b) = ⟨∂_a r ∂_b ψ⟩ antisymmetric (Lemma 3.2). This cancels Katzenberger's Itô drift against the Wong–Zakai term, and a linear-response argument gives zero drift for any smooth stationary noise. Conclusion: no selection at order v.
- **Euler update (Theorem 4.2).** For an Euler update with stationary correlated draws, the selection is Katzenberger's drift of a filtered covariance with an explicit filter f_il, valid for every h < 2/λ_max.
- **Readings (Corollary 4.3).** Damian–Ma–Lee's regulariser for fresh draws; τ_c⁻¹ log det(I + τ_c H) for a finely sampled Markov noise; concavity for every colour in the continuum.
- **Finite step and slow noise.** Concavity fails at a finite step (Prop. 4.4). A slow Markov noise acts as a thermal noise in the flow's metric on the normal space (Prop. 4.5). Prop. 4.6 gives exact per-step identities.
- **Factorised coefficients (§5).** On w = u⊙u − v⊙v the reduced flow is the gradient flow of R = tr Ψ(K) (Theorem 5.1). White noise gives a weighted ℓ¹ norm; slow noise gives Fazel's log-det, which counts the support. A weighted dual certificate holds sparse fits for the reduced flow.
- **Kernel network (§6).** Closed forms on the masses of an adaptive kernel network, and a multistable coloured selection.
- **Numerics (§7).** The slow reduced flow recovers sparse truths where basis pursuit fails, at rates comparable to sparse Bayesian learning.

## Overall assessment
**Recommendation: major revision.** For SIAM J. Math. Data Sci. the core is a good fit after revision. For COLT/NeurIPS it needs a much shorter version built on §§3–5. For IEEE TAC/Automatica the adaptive-control reading is not yet supported (Major 2).

What holds up:
- The core mathematics is correct, and the antisymmetry mechanism is elegant and, as far as I can find, new.
- I checked the central closed forms independently: f_il, the fresh and continuum limits, the slow-noise limit, Cor. 4.3(b), the Euler/Heun identities, Cauchy–Binet, Theorem 5.1's reduced flow and the §6 closed forms.
- On a generic three-parameter model with no symmetry, a Monte Carlo of the Euler update matches Theorem 4.2 within 0.2% (fresh, AR(1) ρ = 0.6 and −0.5, and hλ_max ≈ 0.97). Itô and Wong–Zakai cancel there to 1·10⁻⁶ relative.
- The numbers are faithfully copied: both md reports regenerate byte-identically from their raw files, and the number check passes.
- The paper is unusually honest about where Katzenberger's range ends.

What needs revision:
1. Several results labelled "proved" are formal linear-response or two-scale derivations that no cited limit theorem covers.
2. The motivating coloured-noise examples (sensor drift, slowly changing operating conditions) violate the paper's own hypothesis of independence over the nodes.
3. Prop. 4.4 discards a correct statement from the source and misses an exact closed form that settles the finite-step question.
4. Prop. 5.6 generically covers only the truth; the competitors the flow actually ends at (N-atom fits) are outside it.
5. In §6 the multistable targets converge to the saturations, outside Hypothesis 2.1, and §7 hides magnitude disagreements behind cosines.
6. Several relevant literatures are missing, and the overlap with Shalova–Schlichting–Peletier is unresolved.

## Major comments

**M1. "Proved" is too strong for the coloured and finite-step results (Table 1, l. 75, 78, 81; Theorem 3.3(b), Theorem 4.2, Corollary 4.3, Prop. 4.5).**
- Katzenberger's theorem as stated (Theorem 2.2) is for Itô SDEs with white noise. Li–Wang–Arora cover SGD only as η → 0.
- Theorem 3.3(b) (smooth noise, time scale 1/v at fixed τ_c) and Theorem 4.2 (fixed h, v → 0) are obtained by expanding Φ to second order and averaging at linear response. The text itself says (l. 144): "for a discrete update we use the same name for the mean displacement … at linear response". These are formal derivations.
- My Monte Carlo [own] supports Theorem 4.2's formula, but evidence is not proof.
- Fix: relabel in Table 1 as "derived at linear response (formal)", or prove them.
  - For OU noise, and any noise with a rational spectrum, augment the state with the noise's Markov state (z, ξ). The augmented SDE then has the attracting manifold Γ × {0} and noise of size √v, so Katzenberger's theorem applies directly. This would turn Theorem 3.3(b) into a theorem for that class.
  - For the fixed-h Euler chain, a discrete Katzenberger-type limit (small noise, fixed map with a normally attracting manifold of fixed points) is needed. Otherwise state it as an assumption.

**M2. The hypothesis "homoscedastic and independent over the nodes" excludes the paper's own motivating coloured noises (l. 56, 121, 151–162; Remark 3.5).**
- In instance (ii) the nodes are quadrature points along the regressor's orbit, so one sensor's noise is sampled at successive phases.
- Independence over the nodes needs a correlation time τ_n ≪ T/N (the source says so explicitly, D‴.41 l. 5021). Such a noise is white on the adaptation's slow time scale. A sensor noise coloured on the adaptation time scale is necessarily correlated across the nodes.
- The motivating examples (l. 56: "sensor drifts, an operating condition that changes slowly") are common offsets, which is the case C ∝ 11ᵀ of Remark 3.5. That case is not covered: the Stratonovich drift is nonzero there, and only on the one-particle-per-class network is it O(γ).
- So Theorem 3.3(b) and Theorem 4.2 apply to per-sample label noise that is temporally correlated (e.g. augmentation reused across epochs, per-sample sensors), not to a single sensor in a control loop.
- Fix: say this in §1.1 and §2.3. Rewrite the adaptive-control claims (abstract, title reading, l. 58 "a loop adapting in continuous time does not select … whatever the noise's temporal spectrum") accordingly. Treat the coloured offset C = 11ᵀ explicitly.

**M3. Prop. 4.4 (l. 315–321) can be made exact, and README Inconsistency 4 discarded a correct source statement.**
- Substituting G = ρa/(1−ρa) and R = ρ/(1−ρa) gives, exactly, with x = hλ:

  f(x) = 1 + x(1 − 3ρ + ρx) / ((1 − ρ + ρx)(2 − x)) [own; checked against the paper's formula].

- Consequences for AR(1):
  - For ρ ≤ 1/3, f is increasing on (0, 2).
  - For ρ > 1/3, f < 1 exactly when hλ < 3 − 1/ρ.
  - f has its unique minimum at hλ* = (3ρ − 1)/(2ρ) < 1 [own: f′(x*) = 0 to 1e−10, sign change on both sides, for ρ = 0.4–0.999].
  - So Ψ is concave only on hλ < (3ρ−1)/(2ρ), and that threshold tends to 1 as ρ → 1.
- The source's "f rises again beyond hλ ≈ 0.7–1" is therefore correct for AR(1) with ρ ≥ 0.6, where the threshold lies in 0.67–1. The README replaced it as "unsupported" without checking the formula.
- The remainder O((hλ)²) has a constant of order (1−ρ)⁻² [own: c₂(1−ρ)² → 0.98 at ρ = 0.99]. The expansion is useful only for hλ ≪ 1−ρ, as Cor. 4.3(a) says.
- Fix: state the closed form, the minimiser and the concavity window in Prop. 4.4. Mark hλ* in Fig. 1(c).
- Practical relevance: the step-8 experiments (hλ_max ≈ 1.6) sit in the rising branch, where the log-det reading does not apply.

**M4. §5 and the "sparsity prior": what Prop. 5.6 actually covers (l. 413–429, 504–512).**
- Prop. 5.6 holds fits on k < N atoms in general position.
- For y = Φw* with w* k-sparse and Φ in general position, the truth is generically the only exact fit on fewer than N atoms. Every other basic fit, and so (by Prop. 5.5) every other local minimiser, has exactly N atoms, where Prop. 5.6 says "no conclusion follows".
- So the proposition says: for τ_c beyond a threshold the truth is a strict local minimum. Whether the flow reaches it is entirely a basin question, and the competing minima (where every failure ends; §7.3: R(end) > R(truth) in every failure up to k = 4) lie outside the theory.
- The source states this (D‴.43(v): "the truth is … the only one the slow selection holds robustly; where the flow ends is decided by the basins"). The paper dropped it.
- Fix: restore the sentence. Characterise the N-atom vertices (e.g. compare R at the truth with R at N-atom vertices as τ_c varies). Make clear in the abstract that recovery is a basin property of the reduced ODE, outside Katzenberger's range (Remark 5.2).
- The certificate numbers (Part H) show the gap: at k = 5 all truths are certified at τ_c = 100, yet only 0.15 are recovered.

**M5. §6: the multistable targets converge to the saturation, outside Hypothesis 2.1 (l. 459–476, 495; Supplement S1).**
- In `gauge_label.md` G1a/G1b [md], the stable zeros at γ = 10⁻⁵ are (0.2500, 0.1750, 0.5749) and permutations at τ_c = 10⁵, i.e. two classes at W_min to within about 1.3–1.6·10⁻⁴ relatively.
- The continued branch has W/W_min − 1 → 5.5·10⁻⁴ while λ_soft τ_c ≤ 0.03. The softest normal rate vanishes there, so Remark 2.3 applies: these are statements about an ODE outside the range.
- The main text describes the targets only as "each with one class's mass large" (l. 475, 495). That they hold two classes at their saturation appears only in the supplement.
- Prop. 6.2 predicts a drive "toward the flattest, saturated fits", which leaves its own domain of validity (λ_min τ_c ≫ 1).
- The source's caveat for Prop. 6.3, that "its constant grows like 1/a′ toward a saturation" (D‴.41(ii), l. 5040), was dropped. It matters precisely because the targets approach saturation.
- At γ = 10⁻⁵ the reported zero sets look incomplete by an index count: 3 stable + 1 saddle at τ_c = 100, and 3 stable + 1 unstable at τ_c = 10³ [md]. The γ = 0.01 sets satisfy the count. State the vector field's behaviour on the saturation boundary, or the search's completeness.
- Fix: say in the main text that the slow-noise targets sit at saturations; restore the 1/a′ caveat; consider moving §6 to a separate note (see M7).

**M6. Prop. 4.5 / 6.2 numerics: cosines hide magnitude disagreements (l. 495).**
- The paper reports only the cosine of the slow noise's drift against the dither's (1.00000; 0.99883–0.99989).
- Part G2 [md] also prints the magnitudes. At γ = 0.01:
  - the slow noise's drift is 1.83–3.27 times the dither's normal-part prediction, with −6.9 on one component at W = 1/3;
  - the dither's own drift is 2.8–4.1 times its normal part's, with −14 on that component.
- So at γ = 0.01 the tangent parts that Prop. 4.5 sets aside dominate the drift; only at γ = 10⁻⁵ are both ratios 1.00.
- Fix: print the magnitudes, and state that Prop. 4.5's identification is useful only as γ → 0.

**M7. Length, venue and self-containedness.**
- At 25 pp. with three loosely coupled themes — §3–4 (general mechanism), §5 (sparse recovery), §6 (kernel-network masses) — the paper is long for SIMODS and far over COLT/NeurIPS limits.
- §6 depends on a network whose bump parameters are not printed ("listed in the code release", l. 480), so its §7.3 numbers cannot be reproduced from the paper.
- Recommendation: §§2–5 plus the generic-network checks form a coherent SIMODS paper. §6 belongs in a separate note or the supplement. A TAC/Automatica version would need a closed-loop plant and a noise model that satisfies M2.

**M8. Literature (details below).**
- The overlap with Shalova–Schlichting–Peletier is left open (README "Open items"). That paper's results on how noise structure sets the time scale of the evolution along the zero-loss set (label noise, dropout, minibatch) must be compared result by result with Theorem 3.3.
- Ziyin–Wang–Li–Wu (NeurIPS 2024), cited in the source as the general frame for Prop. 4.6's identities, was dropped.
- The statement (l. 99) that adaptive-control robustness theory "treats the noise's effect on a unique equilibrium" is inaccurate. Parameter drift and bursting without persistent excitation are exactly noise-driven motion along a set of fits (Egardt; Ioannou–Kokotović; Anderson 1985 [snippet]); this literature must be cited and contrasted.

## Minor comments
1. l. 338, 346 (Prop. 4.6(a)). The proof writes "c_n = ξ_n − r_n" and "r_{n+1} = (1 − hλ)r_n + hλξ_n", which use r = f̂ − y. That contradicts §2.1 (r = f − f̂). With §2.1's convention c = r + ξ and r_{n+1} = (1−hλ)r_n − hλξ_n [own, numerical check]. The factor 1 − h²c² is unaffected; E[c²] is not.
2. l. 338. "the Itô equation's dQ = −σ_ξ²Q dt" should be dQ = −vQ dt (v = hσ_ξ² here; Table 2 defines σ_ξ² as the per-step variance). The md's σ² is the intensity.
3. l. 310 (proof of Cor. 4.3(b)). The factor m is dropped: Σ_f = m f(H)H and X = (m/2) f(H). The statement −(m v_h/4)P∇tr Ψ(H) is correct; I confirmed it independently against Katzenberger's ∂²Φ[Σ_f] to 1e−6 [own].
4. l. 382, 386 (Prop. 5.3). "Damian–Ma–Lee's strictly convex tr Ψ(K)": Ψ is strictly convex, but R is not strictly convex in D or (w, P) when D ↦ K has a kernel (d > N(N+1)/2). The proof shows only convexity, so the minimiser's uniqueness is unproved.
5. l. 232. Theorem 3.3's closing sentence needs "at order v" (cf. §8 item 1).
6. l. 316. Prop. 4.4's last sentence ("White noise is the continuum's only selection linear in H") belongs to Cor. 4.3(d).
7. l. 305, 406. "f ≈ (λτ_c)⁻²" holds only for λτ_c ≫ 1. For Gaussian correlations at λτ_c = 0.9, f = 0.377 against 1.23 [own].
8. l. 260 (Prop. 4.1(b)). The chart term is stated with Σ, which is correct for fresh draws as hλ → 0. Add "up to O(hλ_max)"; for coloured draws or large steps the term should use Σ_f.
9. l. 483. "Hypothesis 2.1 holds with room": λ_min = 3.06·10⁻⁴ against v = 5·10⁻⁵ in §7.2, so v/λ_min ≈ 0.16. Quantify the separation instead.
10. l. 488. §7.2 omits the ρ = 0.95 slopes (1.0043–1.0265 ± 0.02–0.07), which Fig. 3(b) shows. It also omits that the prediction at W₀ alone is off by up to 2.2% at step 0.25, and at h = 8, ρ = −0.5 gives 0.9141 and 1.4573 [md C2]. The path-averaged prediction uses the measured mean path; say so.
11. l. 504. "From a large start (α = 2) … the slow flow recovers 0.35 and 0.30 (k = 3, 4), basis pursuit's rate (0.45, 0.25)" — at k = 3 the flow is worse than basis pursuit. These are also the first 20 instances, not the 40 behind the headline. Part C's slow flows run to a fixed horizon S = 400 + 60√(1+τ_c), not to convergence; state it. On the 20 shared instances, horizon and converged results agree [raw, own].
12. Statistics (§7.3). Rates from 20–40 instances need binomial intervals (±0.07–0.08) and paired counts. From the raw file [own]:
    - slow (τ_c = 100) vs basis pursuit: 15/0 and 17/0 discordant instances at k = 3, 4 (exact McNemar p = 6·10⁻⁵, 1.5·10⁻⁵);
    - slow vs sparse Bayesian learning: 1/6 and 4/2 (p = 0.13, 0.69), i.e. indistinguishable.

    "Falls like 1/τ_c" (l. 512) holds only asymptotically: the median ratio falls by 5.4× between τ_c = 10 and 100.
13. l. 341, 510 (Prop. 4.6(d)). Heun's residual selection matters at practical steps. On my generic model, Heun with the same draw selects 0.44% of Euler's at hλ_max = 0.16 but 12.9% at hλ_max = 0.97 [own MC]. Qualify "Stratonovich-consistent updates do not select" (abstract, §1.1) as an h → 0 statement.
14. Remark 3.5 / Prop. 6.3. For noises cross-correlated over the nodes, the multidimensional Wong–Zakai limit carries a Lévy-area term. `gauge_label.md` A2 computes it (below 4·10⁻¹³ there); mention it.
15. Notation overload, which hampers reading:
    - R: correlation function, R(a) = G/a, and the potential;
    - P: ∂Φ|_Γ, u_i v_i, and ⟨Φ_k²⟩;
    - Q: I − P, Q_il, u² − v², and Q_k;
    - Φ: limit map, design, and bump;
    - K: matrix, number of particles, and 𝒦;
    - also X, s, h, G, Γ_k, and n vs N in Prop. 4.6(b).

    Table 2 does not resolve these. Rename.
16. Fig. 2. The sparse Bayesian learning markers at k = 2 and 5 are hidden under the flow curves; draw them on top or as a line. Fig. 1(c) should mark hλ* (M3).
17. l. 97. Intuition worth adding: for N = 1 the noisy Stratonovich flow follows the fixed integral curves of ∇f̂, so Φ is conserved pathwise. The theorem's content is the averaged cancellation for N ≥ 2, where span{∇f̂_x} is generally non-integrable. My test model is such a case, and the cancellation holds there [own].
18. Table 1, Prop. 5.6: "stability" holds for trajectories with every P_i ≠ 0, because the field is singular at the fit. Say so in the statement.

## Verification I did
- **Reproduction [own].** Copied both raw npz files to separate output directories, ran `exp_colour_selection.py --from-raw` (CS_OUT) and `exp_gauge_label.py --from-raw` (GL_OUT) under `python -I`, and compared with `cmp`: both byte-identical.
- **Number check [own].** `numcheck15.py`: 321 distinct decimal tokens, 0 not found. The match is by value only. By hand I traced §7.3's Parts B/C/G/H/J numbers and §7.2's C/C2/D/E/F numbers to their rows; all correct. The issues are selective reporting (Major 6, Minor 10–11), not mis-copying.
- **Closed forms [own]** (`scripts/check_filter.py`, exact rationals and brute-force series):
  - f_il against the defining series (AR(1) at ρ = −0.5, 0.3, 0.8; a generic MA correlation; hλ up to 1.9): agreement to 12 digits;
  - fresh draws f = 2/(2−hλ): exact;
  - AR(1) slope (1−3ρ)/(2(1−ρ)): exact, with curvature ∝ (1−ρ)⁻²;
  - 1 + 2G(−1) = (1−ρ)/(1+ρ);
  - continuum limit (1 + 2x_i x_l/(x_i + x_l))/((1+x_i)(1+x_l)): O(h) convergence;
  - slow limit 2/(λ_i + λ_l − hλ_iλ_l);
  - spectral formula for OU noise (exact) and Gaussian correlations; Ψ(∞) = √(π/2)/τ_c to 6 digits;
  - the closed form of M3 and its minimiser (3ρ−1)/(2ρ).
- **Identities [own]** (`scripts/check_identities.py`, `scripts/check_reduced.py`):
  - Euler factor 1 − h²c² and Heun's factor 1 − (h²/4)(c−c̃)² + (h⁴/4)c²c̃² (exact), with c − c̃ = hλc + h²c²uv exactly;
  - diagonal-network factor 1 − 4h²g² (2e−16);
  - Cauchy–Binet (1.4e−15); ∂R/∂D = (4/N)ω (2.4e−9); the sparse-fit asymptotics;
  - a′² − aa″ = 1 − q², a′² − 2aa″ = (1−q)(1+3q), ∂_W(W − c/W)² = 2W(1−q²), h″ = 2 + 6c²/W⁴, and the a-chart potential's derivative;
  - Theorem 5.1's reduced flow against −¼P_T∇R computed in (u, v): rel. 4e−10 (dP) and 1.5e−9 (dw).
- **Generic model, no symmetry [own]** (`scripts/mc_generic.py`; z ∈ ℝ³, two nodes, f̂_x = g_x(z) non-integrable, Γ a curve, gauge coordinate Φ(z)₃, Φ computed by RK4 to convergence):
  - Theorem 3.3(a): Itô ½∂²φ[Σ] = −0.0644821 and Wong–Zakai +0.0644820, sum −7.6e−8 (finite-difference floor).
  - Cor. 4.3(b) against Katzenberger's ∂²Φ[Σ_f]: agreement to 1e−6 for six (h, ρ) cases, including hλ_max = 0.97 with ρ = −0.5 and 0.9.
  - Monte Carlo of the Euler update (1000 paths, v_h = 1e−5, T = 1000), measured over predicted slope: 1.0002 ± 0.0014 (fresh, h = 0.2), 1.0002 ± 0.002 (ρ = 0.6), 1.0025 ± 0.0018 (ρ = −0.5), 1.0022 ± 0.0038 (fresh, h = 1.2), 0.9977 ± 0.005 (ρ = −0.5, h = 1.2).
  - Heun (same draw): 0.44% (h = 0.2) and 12.9% (h = 1.2) of Euler's, consistent with (hλ)²/4.
- **Raw analysis [own].** Part C's horizon flows against Part F's converged flows on the 20 shared instances: identical outcomes. Paired McNemar tests as in Minor 12.
- **Source comparison [read].** Dropped from the sources:
  - the 1/a′ caveat of Prop. 6.3;
  - D‴.43(v)'s "only the truth is held robustly";
  - G1's approach to saturation (main text);
  - G2's magnitude ratios;
  - C's at-W₀ deviations;
  - the Ziyin–Wang–Li–Wu citation.

  Correctly narrowed: D‴.43(v)'s attraction, and the "O(h³)" remainder. Wrongly discarded: D‴.43(iv)'s "0.7–1" (Major 3).
- **Rendered PDF and figures [own].** 25 pages. Fig. 1(c) agrees with the closed form (ρ = 0.5 dips to 0.889 at hλ = 0.5).

## Literature
- **Shalova–Schlichting–Peletier (JMLR 27, 2026) [snippet].** Their abstract states that the structure of the noise determines both the form and the time scale of the evolution along the zero-loss set (dropout, label noise, minibatch noise). According to a third-party summary (not checked against the paper), they recover Li–Wang–Arora for label noise. The paper cites them in one clause (l. 97). A comparison is required: does their small-step framework produce a Stratonovich/Itô distinction or a vanishing label-noise drift for any noise class? arxiv.org and jmlr.org were unreachable from this session, so I could not read the theorem; the hosts would have to be allowed in the session's network settings.
- **Missing, should be added:**
  - Ziyin–Wang–Li–Wu, Parameter symmetry and noise equilibrium of SGD (NeurIPS 2024) [snippet], directly relevant to Prop. 4.6 and cited in the source;
  - adaptive-control parameter drift and bursting without persistent excitation: Egardt (1979), Ioannou–Kokotović (1983/84), Anderson (1985, Automatica), Anderson et al. (1986, averaging), Sastry–Bodson (1989) [snippet/memory];
  - Markov-chain SGD: Sun–Sun–Yin (NeurIPS 2018), Even (ICML 2023) [memory];
  - finite-step and ordering regularisers: Smith–Dherin–Barrett–De (ICLR 2021) [memory];
  - step size on diagonal networks: Nacson et al. (ICML 2022) [memory];
  - large-step SGD and sparse features: Andriushchenko–Varre–Pillaud-Vivien–Flammarion (ICML 2023) [memory];
  - LMS and recursive estimation under correlated noise: Solo–Kong (1995) [memory].
- **Fair as stated:** Blanc et al., Damian–Ma–Lee (the white case of Theorem 4.2 is exactly their regulariser, credited), HaoChen et al., Li–Wang–Arora, Pillaud-Vivien–Reygner–Flammarion (Itô reading, credited), Even et al., Pesme et al., Wong–Zakai, Kupferman–Pavliotis–Stuart, Kunin et al., Jung–Hänggi and Maggi et al. for the 1/(1+λτ) filter and the determinant factor (credited).
- **Novelty.** The claims for the antisymmetry cancellation and the coloured filter appear justified, pending the SSP comparison. The novelty of "concavity for every colour" in the continuum is modest given the active-OU literature. Ljung, Benveniste–Métivier–Priouret and Kushner–Yin are cited, but the sentence characterising them (l. 99) needs the M8 correction.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3ek-1 | major | Theorem 3.3(b), Theorem 4.2, Cor. 4.3 and Prop. 4.5 are labelled "proved", but they are formal linear-response derivations; Katzenberger/LWA do not cover coloured noise at fixed τ_c or fixed h (l. 144 admits the naming) [read] | Table 1 l. 75, 78, 81: "derived (formal)". Prove the OU / rational-spectrum case by state augmentation plus Katzenberger; state the fixed-h limit as an assumption |
| D3ek-2 | major | Independence over the nodes is incompatible with a single sensor's noise coloured on the adaptation time scale; the motivating sensor drift is the C ∝ 11ᵀ case of Remark 3.5, not covered (source l. 5021: τ_n ≪ T/N) [read] | §1.1 l. 56, §2.3, abstract: restrict the claims to per-sample label noise; treat the coloured offset explicitly |
| D3ek-3 | major | Prop. 4.4 misses the exact form f = 1 + x(1−3ρ+ρx)/((1−ρ+ρx)(2−x)), with minimum at hλ = (3ρ−1)/(2ρ) and f < 1 iff hλ < 3 − 1/ρ; the source's "rises beyond 0.7–1" is correct for ρ ≥ 0.6 (README Inconsistency 4 is wrong); the O((hλ)²) constant is ∝ (1−ρ)⁻² [own] | Rewrite Prop. 4.4 (l. 315–321); correct the README; mark hλ* in Fig. 1(c) |
| D3ek-4 | major | Prop. 5.6 generically holds only the truth (the only exact fit on fewer than N atoms); all failures end at N-atom fits where it is silent; the source's sentence saying so was dropped (D‴.43(v)) [read][md] | Restore it after Prop. 5.6; analyse the N-atom vertices; abstract: recovery is a basin property of the reduced ODE |
| D3ek-5 | major | §6 slow-noise targets converge to the saturations (two classes within about 1.5·10⁻⁴ of W_min at γ = 10⁻⁵), outside Hyp. 2.1; the main text says only "one class's mass large"; Prop. 6.3's 1/a′ caveat was dropped [md G1][read l. 5040] | l. 475, 495: state the saturation; restore the caveat in Prop. 6.3; consider moving §6 out |
| D3ek-6 | moderate | G2 magnitudes are hidden: at γ = 0.01 the slow noise's drift is 1.8–3.3× the dither's normal part (−6.9 on one component), the dither's own drift 2.8–4.1× (−14); only cosines are printed [md] | l. 495: print the ratios; qualify Prop. 4.5/6.2 as γ → 0 statements |
| D3ek-7 | moderate | Literature: SSP overlap unresolved; Ziyin–Wang–Li–Wu dropped; adaptive-control parameter drift and bursting, Markov-chain SGD, Smith et al., Nacson et al., Andriushchenko et al. missing; l. 99's "unique equilibrium" characterisation is inaccurate [snippet][memory] | §1.3 |
| D3ek-8 | moderate | Recovery rates (20–40 instances) have no intervals; slow vs SBL is not significant (McNemar p = 0.13, 0.69) while slow vs BP is (p = 6e−5, 1.5e−5) [own, raw] | §7.3, Fig. 2: binomial CIs and paired counts |
| D3ek-9 | moderate | l. 504: the α = 2 result (0.35 < 0.45) is called "basis pursuit's rate", on a different 20-instance set; the Part C horizon S = 400 + 60√(1+τ_c) is unstated [md][read script] | Reword; state instance sets and horizons |
| D3ek-10 | moderate | Length 25 pp.; §6 not reproducible from the paper (bump parameters only in the code); three loosely coupled themes [read] | Split, or move §6 to the supplement; short version for COLT/NeurIPS |
| D3ek-11 | minor | Prop. 4.6(a) uses r = f̂ − y (c = ξ − r, r_{n+1} = … + hλξ), contradicting §2.1; "dQ = −σ_ξ²Q dt" should be −vQ dt [own] | l. 338, 346 |
| D3ek-12 | minor | Proof of Cor. 4.3(b) drops m (Σ_f = m f(H)H) [own] | l. 310 |
| D3ek-13 | minor | "strictly convex tr Ψ(K)": R is not strictly convex in D or (w, P) when D ↦ K has a kernel [own reasoning] | l. 382, 386 |
| D3ek-14 | minor | Theorem 3.3 closing sentence lacks "at order v"; Prop. 4.4's last sentence misplaced; "f ≈ (λτ_c)⁻²" only for λτ_c ≫ 1 (0.377 vs 1.23 at λτ_c = 0.9) [own] | l. 232, 316, 305, 406 |
| D3ek-15 | minor | Prop. 4.1(b)'s chart term is valid for fresh draws and small h only | l. 260: add O(hλ); use Σ_f in general |
| D3ek-16 | minor | "Hypothesis 2.1 holds with room": v/λ_min ≈ 0.16 at Part C's v [md] | l. 483: quantify |
| D3ek-17 | minor | §7.2 omits the ρ = 0.95 slopes (shown in Fig. 3b) and the at-W₀ deviations (0.9141, 1.4573 at h = 8, ρ = −0.5); the path-averaged prediction uses the measured path [md] | l. 488 |
| D3ek-18 | minor | Heun selects 12.9% of Euler at hλ_max ≈ 0.97 [own MC]; "Stratonovich-consistent updates do not select" holds only as h → 0 | Abstract, §1.1, Prop. 4.6(d) |
| D3ek-19 | minor | Lévy-area term for cross-correlated node noise not mentioned; the γ = 10⁻⁵ zero sets fail the index count (missing saddles or boundary zeros) [md] | Remark 3.5; Supplement S1 |
| D3ek-20 | minor | Notation overload: R ×3, P ×3, Q ×4, Φ ×3, K ×3, X, s, h, G, Γ_k, n/N [read] | Rename; complete Table 2 |
| D3ek-21 | minor | Fig. 2: SBL markers hidden at k = 2, 5 | `make_figs.py`: draw them on top |
| D3ek-22 | minor | Prop. 5.6's stability holds only for trajectories with all P_i ≠ 0 (singular field at the fit) | l. 418 |

## Sources
- [read] `manuscript/paper15/paper15.tex`, `supplement.tex`, `README.md`, `numcheck15.py`, `make_figs.py` (SBL part), the rendered `paper15.pdf` and fig1–3 PDFs.
- [read] `PARADIGM.md` l. 5016–5247 (D‴.41) and l. 5350–5440 (D‴.43); `experiments/exp_colour_selection.py` (c_job, f_job, main, constants); `experiments/exp_gauge_label.py` (main).
- [md]/[raw] `results/chaos/colour_selection.md` and `colour_selection_raw.npz`; `results/chaos/gauge_label.md` and `gauge_label_raw.npz`.
- [own] Scripts in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref142/scripts/` (`check_filter.py`, `check_identities.py`, `check_reduced.py`, `mc_generic.py`). Monte Carlo output: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref142/mc_out.txt`. Regenerated reports: `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref142/out_cs/colour_selection.md` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref142/out_gl/gauge_label.md`.
- [snippet] [SSP, JMLR](https://www.jmlr.org/beta/papers/v27/24-0571.html); [SSP, arXiv](https://arxiv.org/abs/2404.12293); [third-party summary of SSP](https://themoonlight.io/review/singular-limit-analysis-of-gradient-descent-with-noise-injection); [Ziyin–Wang–Li–Wu, NeurIPS 2024](https://proceedings.neurips.cc/paper_files/paper/2024/hash/aa933b5abc1be30baece1d230ec575a7-Abstract.html); [Anderson, lack of PE and bursting](https://openresearch-repository.anu.edu.au/entities/publication/b9a832b6-dae1-41f8-bf77-66c2fe5de1df); [historical perspective of adaptive control](https://arxiv.org/pdf/2108.11336); [Tsakalis, bursting](https://tsakalis.faculty.asu.edu/burst.htm); [Li–Wang–Arora](https://arxiv.org/abs/2110.06914). arxiv.org and jmlr.org could not be fetched from this session.

## Lesson
Before removing a source statement as "unsupported by the md", check it against the closed form. The discarded "f rises again beyond hλ ≈ 0.7–1" is exactly the minimiser (3ρ−1)/(2ρ) of the AR(1) filter.

Two related habits:
- A cosine near 1 is not agreement. Print the magnitudes the report already contains (G2's 1.8–3.3× ratios).
- When a theorem's hypothesis carries a modelling meaning (independence over the nodes means white on the orbit scale), test the paper's own motivating examples against it before writing the introduction.
