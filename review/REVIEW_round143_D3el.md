# Referee report D3el, round 143: the Paper 14 manuscript, first review

## Scope

I reviewed `manuscript/paper14/paper14.tex` (25 pp.), `supplement.tex` (7 pp.), the compiled PDFs, `README.md`, `make_figs.py` and `numcheck14.py` at commit 1de112e. I worked on an export at `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref143/repo` with the venv symlinked; the live repository was not touched.

Sources read:
- the PARADIGM.md summary rows ll. 148–154;
- the full entries of D‴.36–D‴.39 (ll. 4495–4926) and D‴.42 (ll. 5248–5340);
- the D‴.40 row, for the scope decision;
- the five md reports `gauge_dither`, `gauge_drift`, `gauge_fixman`, `gauge_plateau` and `gauge_classes`, plus the header of `cohesion_leak.md`.

Venues considered: SIAM J. Applied Dynamical Systems, Multiscale Model. Simul., and COLT/ALT. Line numbers below refer to `paper14.tex`; page numbers to the compiled `paper14.pdf`.

## Summary of the submission

The setting is a Gaussian-kernel network whose masses follow a replicator (Fisher–Rao) law and whose positions follow a gradient law. With one kernel per feature, the exact fits form a gauge manifold M of dimension K−1, parametrised by the masses. A white-noise dither of level τ on the positions, in unconstrained coordinates, then does the following:
- **Diffusion along M (Prop. 3.2).** In Katzenberger's limit the masses diffuse along M with covariance D = O(γ²τ), vanishing at the faces.
- **Drift and free energy.** The masses drift at the same order (Prop. 4.1). In the slow-replicator limit (formal, Derived st. 4.3) the drift is the replicator flow down (τ/2)ℓ, with ℓ = log det ∇²_Θ J. Lemma 4.2 gives the closed form ℓ = 2Σ log(W_k − W_min,k²/W_k) + const, which tends to −∞ at every saturation, so the dither selects saturated fits.
- **Boundary layer.** A layer of width O(γ^{1/2}) slows this but does not stop it (Derived st. 4.4).
- **Fixman-type correction.** Adding κ_F(τ/2)∂_β ℓ to the fitness gives the following:
  - a gradient flow with a log-barrier;
  - cancellation of the drift at κ_F = 1;
  - selection of the D-optimal fit W♯ at κ_F > 1;
  - a density δ^{κ_F−1} near the face (Prop. 5.1–5.3, Derived st. 5.4–5.5).
- **Confinement.** Without a prior there is no stationary law (Prop. 6.1, at fixed masses). The exit to plateaus of saturated approximate fits has Freidlin–Wentzell order (Prop. 6.2). A dead zone on a box restores a law and turns the saturation into a wall (Prop. 6.3).
- **Several particles per class.** With an isotropic dither the class masses carry the same laws, slowed by K/m (§7).
- **Numerics.** Everything numerical is on one three-bump instance (§8).

## Overall assessment

**Recommendation: major revision** for SIADS or MMS. **Not suitable for COLT/ALT** in this form: the proved content is Katzenberger's theorem applied to one architecture plus closed-form algebra, and the learning-theoretic claims are formal and rest on a single instance.

Strengths:
- The setting is clean.
- The closed forms are correct; I re-derived and recomputed every one I checked [own].
- The status labelling (Table 1) is better than average.
- The number provenance is exemplary: `numcheck14.py` reports 0 of 697 tokens unmatched, and the three reports I regenerated from their raw files are byte-identical.
- The Fixman/D-optimal angle is a nice idea.

Six problems need fixing before acceptance:
1. **A dropped caveat on the leak tests.** The instance violates the leak window that the paper itself states, and the source's explanation of why the tests still pass ("that mode barely couples to the masses") was dropped.
2. **Chart dependence of the headline.** The selection of saturated fits depends on which coordinates are dithered: the divergence of ℓ is exactly the Jacobian of the tanh map. The paper's chart discussion covers the masses' chart, not this one.
3. **Confinement versus the implemented floor.** The headline confinement exponent holds only without the floor. The floored correction that was actually run is not normalisable in the limit, and the full-loop runs never reach the regime where the exponent law applies.
4. **The deterministic reading fails at the instance's γ.** The Péclet number of the mass dynamics is about 1 at γ = 10⁻², so "deterministic flow down ℓ" does not describe the runs. My Monte Carlo agrees with Katzenberger's drift only after averaging the drift over the ensemble's spread.
5. **Where Katzenberger fails.** At fixed γ the normal rates do not degenerate at the faces; the obstruction is non-compactness in α. The text says the opposite.
6. **Literature gaps.** Missing: Shalova–Schlichting–Peletier, and the adaptive-control parameter-drift/bursting literature (Egardt, Ioannou–Kokotović, Anderson 1985), which is precisely the paper's motivation.

## Major comments

**M1. The leak's window is violated on the instance, and the source's explanation was dropped (Prop. 3.3, Derived st. 5.5, §8.1, S2.2).**
- **What the text requires.** Derived st. 5.5 (l. 402) assumes τ ≪ ε ≪ the slowest normal rate. Prop. 3.3(b) works in the regime "ε/τ fixed".
- **What the instance has [own, `lin.py`].**
  - The slowest normal rate is 3.059·10⁻⁴ at both γ = 10⁻³ and 10⁻². It is a width mode: mostly χ₂, with a mass share of the g-norm of 1.4·10⁻⁴ at γ = 10⁻² and 2.4·10⁻⁵ at γ = 10⁻³.
  - The leak's stiffness along M at γ = 10⁻² has eigenvalues 7.4·10⁻⁴ and 9.0·10⁻⁴ at ε = 10⁻³ (2.4–2.9 times the slowest rate), and 2.2·10⁻⁴ and 2.7·10⁻⁴ at ε = 3·10⁻⁴ (about equal to it).
- **Consequences for the paper.**
  - The "informative" spread test (§8.1, l. 488: 0.0044/0.0032/0.0034 against 0.0044/0.0031/0.0033) is at an ε outside the window.
  - Table S4's caption says "inside the window" (supplement l. 137) at ε = 3·10⁻⁴, which is the window's edge.
- **What the source says and the paper omits.** PARADIGM D‴.36's Status: "at ε = 10⁻³ the leak's stiffness exceeds it; the reduction holds because that mode barely couples to the masses". D‴.37(iv): "τ ≪ ε ≪ 3·10⁻⁴, which the instance barely meets". The paper keeps only "the time-scale separation … is therefore thin" (l. 185).
- **The window's lower end does not match the derivation.**
  - S1 (supplement l. 42) says the window keeps the spread O((τ/ε)^{1/2}) small.
  - The same derivation gives C = O(D/ε) = O(γ²τ/ε). That spread is O(γ(τ/ε)^{1/2}), confirmed by my OU computation: std 0.0044 at ε = 10⁻³, γ = 10⁻².
  - So the condition the derivation actually needs is spread ≪ δ, i.e. γ²τ/ε ≪ δ², not τ ≪ ε.
- **Fix.** Restore the source's caveat with the numbers. State the window in the derived form, or quantify it. Relabel S2.2 as "at the window's edge".

**M2. The selection of saturated fits comes from the dithered chart (Lemma 4.2, §1.3, abstract).**
- **Where the divergence comes from.** ℓ is the log-determinant in the unconstrained coordinates (α, ζ, χ), in which the noise is constant. Because ∂_α f̂ = W a′ φ, ℓ = ℓ_phys + 2Σ_k log a′_k, where ℓ_phys is the log-determinant in the physical coordinates (a, c, s). With the same noise in physical coordinates, ℓ_phys = 2Σ log W_k + const.
- **What that changes.** That law still pushes the masses away from 1/K, but it has no divergence at W_min,k. The −∞ at the faces, the "no stop" layer and the δ^{κ_F−1} exponent are all the Jacobian of a = a_max tanh(α/a_max), i.e. the push-forward of Lebesgue measure in α.
- **The gap in the text.** The paper's "Itô and Stratonovich" discussion (l. 156 and "Charts", l. 328) concerns the masses' chart. The decisive convention is the positions' noise metric, which the paper never discusses.
- **Novelty claim.** "New are the closed form, which identifies the flattest exact fits as the saturated ones" (l. 111) should say that flatness here is measured in the tanh-unconstrained chart. This is Fixman's own point: free energies on constraint manifolds depend on the reference measure.
- **Fix.** Add a remark that gives ℓ in a general dithered chart and the corresponding selection. Draw the design consequence: dithering physical coordinates removes the divergence, but needs a reflecting boundary.

**M3. The confinement exponent does not describe the implemented correction or the runs (Derived st. 5.4, abstract, §1.3, §8.3).**
- **The floored correction does not confine.** With the floor, the limit equation's density returns to δ^{−1}, which is not normalisable (l. 396). The source states the consequence explicitly: "Normalisability needs the unfloored correction, or a floor that shrinks with τ" (D‴.38(iii)). The paper drops that sentence and keeps an "as τ → 0" log-divergence phrase that is confusing for an equation already in the τ → 0 limit.
- **The exponent law is untested at the instance's τ.** The finite-τ threshold is (κ_Fτ/(m a_max²Γ₁₁))^{1/2} = 0.0030√κ_F [own]: 0.0030 at κ_F = 1 and 0.0042 at κ_F = 2. The tested δ = 0.002–0.004 and the floor's cap δ_cap = 0.00126 lie below it.
- **The full loop at κ_F = 2 never reaches the layer.** The largest |α| over all records is (6.2, 5.6, 8.1) (S2.1). For class 3, |α| = 8.1 means δ ≈ 0.0073; I computed this from the closed form. So the observed confinement (|α_1| ≈ 4.5) comes from the mean drift up ℓ outside the layer, not from δ^{κ_F−1}.
- **The headline should be qualified.** The abstract and §1.3 (l. 115) present "the confinement exponent δ^{κ_F−1}" as a consequence of the correction without saying this.
- **Fix.** State plainly that the implemented floored correction is not confining in Katzenberger's limit. Either use a τ-dependent floor or run unfloored. Separate the exponent law, which is untested at the instance's τ, from the full-loop evidence.

**M4. At γ = 10⁻² the slow law is not deterministic, and the ensemble mean is not the mean path (Derived st. 4.3(a), §8.1).**
- **The spread is not small.** The spread on the drift's time scale is (D·t)^{1/2}: at t = 1/(γτ) = 10⁶ that is √(2.85·10⁻⁸·10⁶) ≈ 0.17, against a distance of 0.083 from 1/3 to the first face.
- **The Péclet number is about 1.** |b|·L/D ≈ (2·10⁻⁷ × 0.08)/2.85·10⁻⁸ ≈ 0.6. The md confirms this: spreads 0.08–0.15 at 10⁶ (`gauge_drift.md` H).
- **My Monte Carlo [own, `mc.py`]** used an independent Euler–Maruyama implementation of (2.3): N = 800 seeds, γ = 10⁻², τ = 10⁻⁴, dt = 0.25, start on M at W = 1/3, window t ∈ [4000, 44000].
  - Measured mean drift over γτ: (−0.123 ± 0.029, +0.083 ± 0.025, +0.040 ± 0.026).
  - Katzenberger's drift at 1/3 (md): (−0.2055, +0.158, +0.0475).
  - The stored field (`gauge_drift_raw.npz` E) averaged over my ensemble, by the trapezoid rule between the window's end points: (−0.140, +0.125, +0.015). This agrees within 0.6, 1.7 and 1.0 standard errors.
  - The slow law (−0.362, +0.242, +0.120) is excluded by about 8 standard errors in the first component, as Derived st. 4.3(b) predicts.
- **Consequence for the paper.** Katzenberger's drift is confirmed, but only after averaging over a spread of 0.035 that already reaches the sign change near the layer. The paper's comparison of the 64-seed mean with the mean path (l. 486) is therefore not a sharp test.
- **Fix.** Say that (a)'s deterministic reading needs γ well below 10⁻², and quote the Péclet number. Compare ensembles with the limit SDE (the md already has 4000 paths in Part H) rather than with the mean path.

**M5. Where Katzenberger's theorem fails is misdescribed (l. 75, l. 182; Open problem 3).**
- **What the text says.** The text attributes the restriction to compact parts of M to normal rates that degenerate near the faces.
- **What I computed [own, `face.py`].** At fixed γ the normal rates do not degenerate approaching the face.
  - At γ = 10⁻² they are 3.05·10⁻⁴, 1.84·10⁻³ and 2.63·10⁻³ from δ = 0.08 down to δ = 5·10⁻⁴.
  - At γ = 10⁻⁴ the mass channel's rate falls to 3.4·10⁻⁵ and then stays there.
- **The actual obstruction** is that M is non-compact in α_k, and the masses' chart degenerates (|∂α/∂W| → ∞). Degeneracy of the rates appears only in the joint limit γ → 0 inside δ*.
- **Why it matters.** Near a face, M is a graph over α_k, so Katzenberger's theorem might be extended in an α-chart at fixed γ. That is exactly where Prop. 5.3, Derived st. 4.4 and Derived st. 5.4 live, and it would turn Open problem 3 from out of reach into tractable.
- **Fix.** Correct both sentences, and mention the α-chart route.
- **A related point (l. 185).** "The slowest normal rate is about 3τ" compares a rate with a noise level. The meaningful numbers are the normal rate over τ in Katzenberger's time (about 3), the Laplace parameter τ/λ_min(∇²_Θ J) = 10⁻⁴/0.0012 ≈ 0.08, and the near-zero mass coupling of that mode.

**M6. The escape statistics do not test Freidlin–Wentzell sharply (Prop. 6.2, §8.4, Fig. 3(b)).**
- **No prefactor in the fit.** The slope of log rate against 1/τ, 0.81 ± 0.14 ΔJ, comes from a fit without a prefactor. The exit has no saddle ("monotone climb"), so the generic boundary-exit prefactor is rate ∝ τ^{−1/2}e^{−ΔJ/τ}. The md's own prefactor data, rate·e^{ΔJ/τ} = 2.25, 1.90 and 1.42·10⁻⁴, scale as τ^{−0.67} [own].
- **Refitting with a prefactor [own, `arr.py`, weighted by 1/var(log rate)]:**
  - with τ^{−1/2}: slope (0.95 ± 0.14) ΔJ;
  - with τ^{−1}: (1.09 ± 0.14) ΔJ.
- **The barriers are not small-noise.** ΔJ/τ = 2.6–5.1 at the tested τ, so the three-point data cannot discriminate between these fits.
- **FW's hypotheses may fail for the domain used.** The escape set used, {max|α_k| < 6}, is unbounded in ζ and χ, which also saturate and form plateaus. FW needs a bounded domain whose closure flows to the exact fit.
- **The time scales are not well separated at the instance.** Extrapolating to τ = 10⁻⁴ (with the τ^{−1/2} prefactor) gives a mean exit time of about 7·10⁶, against 1/(γτ) = 10⁶. So "times short against the exit" in Derived st. 4.3 holds by a factor of about 7. Near W♯ the plateau is lower (7.637τ), so the margin is smaller there.
- **Fix.** Add the prefactor, add more τ values, use a bounded D₀, and state the separation of scales.

**M7. Statements about the long-run or stationary regime go beyond what is proved (Prop. 6.1, Prop. 6.3(b)–(d), abstract).**
- **(a) No stationary law.** Prop. 6.1 is proved at fixed masses. The abstract and §1.2 say "the dithered loop has no stationary law"; for the joint (Θ, β) process this is unproved.
- **(b) The box's stationary law.** Prop. 6.3(b) says that law gives the plateau "the weight of its Laplace order", and Table 1 marks (b) proved. That clause is not proved. The source adds that "there the box's force enters the masses' averaged fitness at the size of the noise's own selection". The paper dropped this with D‴.40, but it can be stated as a limitation without D‴.40. Without it, (c)–(d) read as long-run statements when they hold only before the first exit and at leading order.
- **(c) The vertex endpoint is not observed.** The runs end with each seed 0.099–0.178 from the nearest vertex of P_B.
- **Fix.** Restrict the abstract to fixed masses. Move the Laplace-weight clause to a remark. Add the limitation.

**M8. Labels and statements.**
- **Propositions that contain formal or numerical content.** Prop. 5.1(c)'s "together with Derived st. 4.3 … up to O(γ²τ)" has a numerical remainder. Prop. 5.2's "the masses of the slow loop follow it" rests on Derived st. 4.3. Either split these into Derived statements or put the condition in the statement's first line.
- **Prop. 7.2 omits the diffusion's scaling.** §8.5 cites the diffusion ratios 0.5000 and 0.3333 as evidence for Prop. 7.2, but the proposition states nothing about the diffusion. The ratio K/m = 1/n is an exact two-line identity: both TᵀgT and TᵀgΣgT scale by n. Add it.
- **The fractions' drift.** Remark 7.3 says the fractions' finite-γ drift goes "toward the equal split (numerical)". The only full-loop evidence (S2.4) moves s₁ away from it by −0.0346 ± 0.0119 (2.9 standard errors). The text explains this by saturation, but it should not be cited as support.

**M9. Literature (details under Literature).**
- **Missing:**
  - Shalova–Schlichting–Peletier, the closest modern analogue (the noise structure decides the limit and its time scale);
  - the adaptive-control parameter-drift and bursting literature: Egardt 1979, Ioannou–Kokotović 1983/84, Anderson 1985 (lack of persistent excitation leads to drift and bursting), which is the motivation of §1.1.
- **The direction of the Fixman analogy.** Fixman's correction makes the rigid system sample like the stiff one. Here the correction removes the stiff modes' free energy from the stiff, dithered system. That is the inverse direction, closer to exact ABF with a closed-form bias. One sentence should say so.

## Minor comments

1. l. 141: "the learning law of an adaptive loop averaged over a periodic reference orbit": give a reference or a derivation; it is not self-contained.
2. l. 185 (p. 8): "about 3τ": see M5; also give λ_min(∇²_Θ J) = 0.0012.
3. l. 269 against l. 402: Prop. 3.3's regime "ε/τ fixed" and Derived st. 5.5's "τ ≪ ε" are different regimes; reconcile them. Supplement l. 42 states the spread as O((τ/ε)^{1/2}) where the derivation gives O(γ(τ/ε)^{1/2}).
4. l. 275: the O(ε/γ) bias where w is not representable is "(sketched)" but is missing from Table 1.
5. l. 344 (p. 13): "Within about 0.1δ* … transverse fluctuations … as large as the distance": unsupported; derive it or drop it.
6. l. 396: "diverges logarithmically as τ → 0" is inside a limit equation; rephrase (M3).
7. l. 431: Prop. 6.2 should require D₀ bounded with a smooth boundary.
8. l. 470 and the abstract: "the class masses obey the same laws" should add "slowed by K/m".
9. l. 504 and the Fig. 2(b) caption: the O(γ²τ) remainder is supported only on γ ∈ [10⁻⁴, 3·10⁻⁴]. At γ = 5·10⁻⁴ and 10⁻³ the largest component is 8.2·10⁻⁴ and 7.2·10⁻⁴ (non-monotone, `gauge_fixman.md` A). Extend to γ = 10⁻⁵ or soften "falls like γ".
10. l. 510 and Table S4: the three gains share the same 16 seeds, so the "within 1.03 s.e." checks are not independent. At κ_F = 1 the check compares a near-zero prediction with a near-zero measurement.
11. l. 521: see M6. Quote ΔJ/τ = 2.6–5.1 in the main text.
12. Fig. 3(c) plots the plug-in form while the text argues with the masses form. Plot the masses form, which is the recommended one.
13. **Notation clashes** (the README says they were removed, but some remain):
    - c: centres c_k*, the correction c_k, chart coordinates c (proof of Prop. 3.2(c)), c = W_min² (proof of Lemma 4.2), and c_m;
    - s: width, and fractions s_i (§7);
    - t_k = a_k/a_max, and time t;
    - T: the tangent matrix, and the scalars T_α, T_W;
    - D: diffusion, and the domain D₀.
14. l. 109 and l. 111: the LWA22 isotropic log-pseudo-determinant corollary is correctly attributed (search-checked [snippet]). Check HaoChen et al. and Funaki against the full texts (README: abstracts only).
15. Code availability (l. 537): `exp_gauge_fixman.py --from-raw` silently omits one line of Part G unless `gauge_drift_raw.npz` sits in the same output directory [own]. Document this dependency or fail loudly.
16. l. 53: MSC 91A22 is fine; consider adding 93E15 (stochastic adaptive control) and 60J70.
17. Table 1, Prop. 5.3: "proved (given Thm 2.2)" is right for the formula's asymptotics, but its use inside the layer is outside the theorem; say so in the table.
18. §8.3, D-optimality: also say "maximal determinant of the positions' Fisher information", so ML readers recognise it.

## Verification I did

All [own] unless marked. Scripts are in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref143/own/`. I wrote my own implementation of the instance and imported no repository code, except the stored drift field used for the ensemble average.

- **ℓ in closed form** (`checks.py`):
  - log det of the finite-difference Hessian of J, minus the closed form, is −19.53724 to 3·10⁻⁹ at five points, including (0.26, 0.37, 0.37);
  - ℓ(1/3) = −29.3197 (md −29.320);
  - ℓ(0.3, 0.35, 0.35) − ℓ(1/3) = −0.561974 (md −0.562).
- **The slow law** at the md's four points: (−0.3619, +0.2416, +0.1203), (−0.9303, +0.5157, +0.4146), (+0.1958, −0.3761, +0.1803), (−1.799, +0.9452, +0.8538). These match `gauge_drift.md` C, and the factor m⁻¹ is confirmed (README Inconsistency 3).
- **W♯** = (0.37652, 0.29876, 0.32472), with equal derivatives 13.68802 and |a|/a_max = (0.664, 0.586, 0.616).
- **The h″ formula**: numeric over closed form = 1.0000000 at W = 0.26–0.9.
- **∂_α log(W a′) = −(1 + (a/a_max)²)/a**: ratio 1.000000 at α = 1, 3, 6, 12. The convexity of 1/a(α), and hence the sign of the inside-layer drift, is positive.
- **D** (`lin.py`) at W = 1/3:
  - γ = 10⁻²: D_kk = (2.846, 1.628, 1.807)·10⁻⁸ and off-diagonals (−1.334, −1.512, −0.295)·10⁻⁸, matching md A;
  - γ = 10⁻³: (4.85, 2.24, 2.62)·10⁻¹⁰.
- **The projected correction T_W(TᵀgT)⁻¹∇_WΨ / (γτ) at κ_F = 1**:
  - γ = 10⁻²: (0.26868, −0.1923, −0.07638), against md corrected minus uncorrected (0.26868, −0.19232, −0.07636);
  - γ = 10⁻³: (0.34971, −0.23534, −0.11438), against (0.34968, −0.2353, −0.11437).
- **The leak's OU spreads**: (0.0139, 0.0099, 0.0105) at ε = 10⁻⁴ and (0.0044, 0.0031, 0.0033) at ε = 10⁻³, matching md. The stiffness eigenvalues and the slowest-mode data are in M1.
- **Normal rates near the face**: see M5.
- **Layer and box numbers**:
  - δ* = 0.05774 at γ = 10⁻², giving zero/δ* = 0.664, and 0.566 at γ = 10⁻⁴;
  - δ_cap = 0.0012595;
  - finite-τ threshold 0.0030 and 0.0042;
  - W_B = (0.29471, 0.20630, 0.23577), |a| ≤ 3.3931;
  - Prop. 6.3(d) needs L_α > (3.20, 2.68, 2.87), which L_α = 5 meets.
- **Arrhenius refits** (`arr.py`): see M6.
- **Monte Carlo of the gauge SDE** (`mc.py`; N = 800, 887 s, nice 10, a single process): see M4. The variance slopes of the displacement were (2.77, 1.91, 2.21)·10⁻⁸ against D = (2.85, 1.63, 1.81)·10⁻⁸; the excess is consistent with D varying across the spread.
- **Regeneration from raw.** With the output env vars `GDR_OUT`, `GP_OUT` and `GF_OUT` pointing at scratch directories:
  - `gauge_drift.md` and `gauge_plateau.md` are byte-identical (cmp);
  - `gauge_fixman.md` is identical once `gauge_drift_raw.npz` is also present, and otherwise misses a line (minor 15).
- **`numcheck14.py`**: "697 distinct decimal tokens (per file); instance parameters: 0; not found in the md reports: 0".
- **Build**: paper14 built in scratch (three pdflatex passes): 25 pages, 0 warnings in the log, no undefined references.
- **Hand checks of §8 against the md**: §8.1 ratios and standard errors; the S2.2 "within 1.03" check; the D-row values 4.4 and 4.3; Table S2 against `gauge_drift.md` C and E [md].
- **Sources**: I compared the paper against PARADIGM D‴.36–D‴.39 and D‴.42 line by line for dropped caveats [raw]. Found: the leak window and the "barely couples" explanation (M1); "Normalisability needs the unfloored correction" (M3); the box's force in the stationary fitness (M7).

## Literature

- **Fair and adequate:**
  - Katzenberger [Kat91] together with LWA22. LWA22's isotropic-noise corollary (Brownian motion plus gradient flow of the log pseudo-determinant) is confirmed [snippet].
  - Blanc et al., Damian–Ma–Lee, HaoChen et al. ("Gaussian noise overfits to dense solutions" matches their abstract [memory]), Pesme et al., Pillaud-Vivien et al.
  - Funaki, Freidlin–Wentzell, Bovier et al.
  - Fixman, Helfand, Hinch, Ciccotti–Lelièvre–Vanden-Eijnden, ABF, metadynamics.
  - Foster–Young, Fudenberg–Harris, Imhof, Hofbauer–Imhof, Mertikopoulos–Viossat.
  - Boyd–Sastry, Krstić–Wang, Pomet–Praly, MYULA.
- **Missing, should be added:**
  - A. Shalova, A. Schlichting, M. Peletier, "Singular-limit analysis of gradient descent with noise injection", arXiv:2404.12293, JMLR 27 (2026) [snippet]. It extends LWA22 to general noise injection and shows that the noise structure decides both the limit and its time scale. This is directly relevant to §1.3, to the positions-only noise and to Open problem 1.
  - B. D. O. Anderson, "Adaptive systems, lack of persistency of excitation and bursting phenomena", Automatica 21(3) (1985) 247–258 [snippet; page range differs between indexes]. This is the classical account of the paper's motivating phenomenon.
  - B. Egardt, Stability of Adaptive Controllers, LNCIS 20, Springer 1979 [memory].
  - P. A. Ioannou and P. V. Kokotović, Adaptive Systems with Reduced Models (LNCIS 47, 1983) and Automatica 20 (1984) [memory]: parameter drift, and the original σ-modification, which the paper attributes only to IS96.
  - For the stiff/rigid free energy and effective dynamics: Bornemann–Schütte (1997, homogenization with a strong constraining potential), Legoll–Lelièvre (2010, effective dynamics) [memory].
- **Novelty claims.** "New are the closed form, which identifies the flattest exact fits as the saturated ones" needs the chart qualification (M2). "Confinement exponent" needs the floor qualification (M3). The other novelty claims (the gauge as the slow variable, the masses-versus-plug-in distinction, the dead zone's wall) are accurate as stated.

## Issues

| id | sev. | issue and evidence | proposed fix and where |
|---|---|---|---|
| D3el-1 | major | The leak window τ ≪ ε ≪ λ_slow is violated at ε = 10⁻³ (stiffness 7.4–9.0·10⁻⁴ against 3.06·10⁻⁴) and only met at the edge at ε = 3·10⁻⁴ (2.2–2.7·10⁻⁴). The source's "barely couples"/"barely meets" caveat is dropped. The lower end of the window does not match C = O(γ²τ/ε). [own lin.py; PARADIGM D‴.36 Status, D‴.37(iv)] | Restore the caveat with numbers in §2.5 and §8.1. Restate the window in Derived st. 5.5 and supplement S1. Fix the Table S4 caption. |
| D3el-2 | major | The divergence of ℓ at the faces is the Jacobian of the tanh chart: ℓ = ℓ_phys + 2Σ log a′, with ℓ_phys = 2Σ log W + const. The selection's form depends on the dithered chart, and this is not discussed. [own derivation] | Add a remark after Lemma 4.2. Qualify §1.3 l. 111 and the abstract. |
| D3el-3 | major | The floored correction is not normalisable in the limit; the source's sentence saying so was dropped. The tested δ and δ_cap lie below the finite-τ threshold 0.0030√κ_F, and the κ_F = 2 runs never enter the layer (max|α_3| = 8.1, δ ≈ 0.007). [own; md C, H, E; PARADIGM D‴.38(iii)] | Reword Derived st. 5.4, l. 396, the abstract and §1.3. Separate the full-loop evidence from the exponent law; use a τ-dependent floor or run unfloored. |
| D3el-4 | major | At γ = 10⁻² the Péclet number is about 0.6 and the spread at 1/(γτ) is about 0.17 against a distance of 0.08 to the face, so the reading "deterministic flow" fails. Own MC (N = 800): drift/(γτ) (−0.123 ± 0.029, +0.083 ± 0.025, +0.040 ± 0.026) matches the ensemble-averaged Katzenberger field (−0.140, +0.125, +0.015), not the value at the mean path (−0.2055, +0.158, +0.0475); the slow law is excluded. [own mc.py] | Derived st. 4.3(a): state the γ range. §8.1: compare with the limit SDE (md H) rather than the mean path. |
| D3el-5 | major | "Normal rates nondegenerate on compact parts" is wrong at fixed γ: the rates stay at ≥ 3.05·10⁻⁴ down to δ = 5·10⁻⁴ (γ = 10⁻²). The obstruction is non-compactness in α. [own face.py] | Fix ll. 75 and 182. Recast Open problem 3 as an α-chart extension. |
| D3el-6 | major | The Arrhenius test has no prefactor. With rate ∝ τ^{−1/2} the slope is 0.95 ± 0.14 ΔJ, with τ^{−1} it is 1.09; ΔJ/τ = 2.6–5.1. The escape set is unbounded in ζ and χ. The extrapolated exit time at τ = 10⁻⁴ (about 7·10⁶) is only about 7 times 1/(γτ). [own arr.py; md C] | §8.4, Fig. 3(b), Prop. 6.2: add the prefactor, more τ values and a bounded D₀; state the separation of scales. |
| D3el-7 | major | Overclaims about the stationary regime: the abstract's "no stationary law" for the loop (proved at fixed masses only); Prop. 6.3(b)'s Laplace-weight clause marked proved; the box force's effect on the stationary fitness omitted. | Abstract and §1.2; split Prop. 6.3(b); add the limitation to Prop. 6.3(c)–(d). |
| D3el-8 | major | Literature gaps: Shalova–Schlichting–Peletier; Egardt, Ioannou–Kokotović, Anderson 1985; the inverse direction of the Fixman analogy. [snippet, memory] | §1.1 and §1.3, bibliography. |
| D3el-9 | minor | Propositions mixing proved and formal content (5.1(c), 5.2). Prop. 7.2 omits the exact diffusion ratio K/m that §8.5 cites. Remark 7.3's "toward the equal split" is contradicted by the full loop (2.9 s.e. away). | Relabel; add the diffusion identity to Prop. 7.2; reword Remark 7.3. |
| D3el-10 | minor | The O(γ²τ) remainder is supported only on γ ∈ [10⁻⁴, 3·10⁻⁴]; the curve is non-monotone at 5·10⁻⁴–10⁻³. [md A] | l. 504 and the Fig. 2(b) caption: extend γ or soften. |
| D3el-11 | minor | "About 3τ" compares a rate with a noise level; the "0.1δ*" claim is unsupported; the "as τ → 0" phrasing. | ll. 185, 344, 396. |
| D3el-12 | minor | Notation clashes (c, s, t, T, D). | Table 2 and the text. |
| D3el-13 | minor | `exp_gauge_fixman.py --from-raw` silently depends on `gauge_drift_raw.npz` in the same output directory. [own] | Document it in the availability statement, or make the script fail loudly. |
| D3el-14 | minor | Fig. 3(c) shows the plug-in form; the recommended form is the masses form. | Replot from Part E. |
| D3el-15 | minor | Paired seeds across the three gains, and a null-against-null comparison at κ_F = 1, in the pinned-state tests. | Say so in S2.2 and §8.3. |
| D3el-16 | minor | The averaged-loop provenance at l. 141 is uncited; consider MSC 93E15. | §2.1 and front matter. |

## Sources

- [read] `manuscript/paper14/paper14.tex` (all 632 lines), `supplement.tex` (all), `README.md` (all), the compiled PDF (page mapping).
- [raw] PARADIGM.md ll. 148–154 and 4495–4926, D‴.42 ll. 5248–5340.
- [md] `results/chaos/gauge_dither.md`, `gauge_drift.md`, `gauge_fixman.md` (A–D), `gauge_plateau.md` (via the supplement and the from-raw cmp), `gauge_classes.md` (via S2.4).
- [own] Scripts `inst.py`, `checks.py`, `lin.py`, `face.py`, `mc.py` (with output `mc_out.npz`), `arr.py` in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref143/own/`.
- [own] Regenerated md reports in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref143/out_drift/`, `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref143/out_plat/` and `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref143/out_fix/`; scratch build in `/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/d38ad7aa-9a51-5611-a21f-499f621de00a/scratchpad/ref143/build/`.
- [snippet] Web searches: [Shalova–Schlichting–Peletier, JMLR](https://www.jmlr.org/papers/v27/24-0571.html), [arXiv:2404.12293](https://arxiv.org/abs/2404.12293?context=cs); [Anderson 1985, ANU repository](https://openresearch-repository.anu.edu.au/entities/publication/b9a832b6-dae1-41f8-bf77-66c2fe5de1df), [KB index with pp. 247–258](https://ibl.kb.nl/articles/0000000010030/21/3); [Li–Wang–Arora, arXiv:2110.06914](https://arxiv.org/pdf/2110.06914).
- [memory] Egardt 1979; Ioannou–Kokotović 1983/84; Bornemann–Schütte 1997; Legoll–Lelièvre 2010; the boundary-exit prefactor; HaoChen et al.'s abstract.

## Lesson

A closed form that checks to ten digits can still mislead in three ways. First, through the chart in which the noise is constant: here that chart is the whole reason ℓ diverges. Second, through the regime the runs actually visit: the floor is never engaged at κ_F = 2, and the Péclet number is about 1 at γ = 10⁻². Third, through caveats that were present in the source and lost in transcription. A Monte Carlo compared against the mean path rather than the ensemble-averaged field would have reported a spurious 3σ disagreement; compared properly, it confirms Katzenberger's drift and rejects the slow law. Compare like with like before declaring agreement or failure.
