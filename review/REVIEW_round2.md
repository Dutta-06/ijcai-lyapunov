# Adversarial review, round 2 (2 refuters + judge per target)

| target | judge verdict | fatal | major | minor | cosmetic |
|---|---|---|---|---|---|
| H | holds_with_fixes | 0 | 2 | 7 | 4 |
| G | holds_with_fixes | 0 | 2 | 7 | 3 |
| C | holds_with_fixes | 0 | 1 | 4 | 5 |
| D | holds_with_fixes | 0 | 0 | 4 | 5 |
| E | holds_with_fixes | 0 | 0 | 6 | 5 |
| AB | holds_with_fixes | 0 | 1 | 7 | 6 |
| F | holds_with_fixes | 0 | 0 | 3 | 6 |

## H — verdict: holds_with_fixes

Verdict: holds_with_fixes.

I judged <project root>\PARADIGM.md as saved at 22:37, which is after both reviews. I also checked <project root>\src\core\langevin.py, <project root>\tests\test_langevin.py (4/4 pass), the result tables and the extracted Chizat 2022 text (mfld.txt). I did not edit any project file.

My scripts are in <scratchpad>\H_judge\:
- bounds.py
- avg.py
- kernel.py
- gibbs.py
- exact_fp.py
- misc.py

**What holds (checked independently)**
- **Whitening.** ζ = G^{-1/2}ϑ turns (𝔄) into [S3] eq. (12): drift −∇_ζṼ, noise √(2τ)dB, confinement (σ/2)‖ζ−ζ⁰‖². W_{G⁻¹} becomes W2, and H shifts by ½log det G.
- **First variation.** It is exact, because L̄ is quadratic in μ.
- **Assumption 1 holds for every μ in P2.** For every probability measure, ‖r_μ‖ ≤ F_ρ+K_Φ on supp ρ. ∇_ϑΦ̃ and ∇²_ϑΦ̃ are bounded uniformly on ℝⁿ×ℝ^p: the sampled Hessian sup is about 160 and does not grow up to |ϑ| ≈ 500. Lipschitz continuity in μ follows from W1 ≤ W2.
- **Assumption 2 holds.**
- **Assumption 3 holds.**
  - |V_L| ≤ K_Φ(F_ρ+K_Φ); the largest value found is 10.6 against a bound of 21.0.
  - The Holley–Stroock oscillation is ≤ 2K_Φ(F_ρ+K_Φ)/(k_eτ).
  - The Bakry–Émery constant is σ/τ.
  - This is exactly [S3] Prop. 5.1 for the square loss, with λ→σ, K→K_Φ, K'→F_ρ and loss weight 1/k_e.
- **Citations.** Thm 3.2 and Cor 3.3 are transcribed correctly. Thm 4.1's hypothesis holds with C₀ = σ/τ₀ for τ ≤ τ₀, which is all its proof uses.
- **Constants.** F_ρ = 1.2469 and α* = 4.198. The log of the guaranteed rate is −8400.5 at (σ, τ) = (0.002, 5e-4) and −1052.6 at (0.02, 0.004).
- **Code.** The formulas in langevin.py match the theorem.
- **Current evidence is consistent.** For the current G = diag(4,1,1), a grid computation gives L̄(μ*_τ) = 0.111 (56³ and 72³ grids, fixed-point residual 1e-15). This equals the m = 32 common loss of 0.106–0.110.

**Decisions on each reviewer issue** (R1 = first reviewer, R2 = second; duplicates merged)
- **R1#1 + R2#3 (H(i) averaging): CONFIRMED, major.** I reject R2's 'minor' rating.
  - The word 'formal' appears only in REVIEW_LOG.md and LITERATURE_TRANSPORT.md, not in PARADIGM.md.
  - The cited test checks only an algebraic identity.
  - At frozen μ the drift error is 4.2% at k_e = 10 for n = 1 and 9.3% for n = 2. It scales as 1/k_e and is independent of the gains.
- **R1#2 + R2#1 (Conjecture H.6 ill-posed): CONFIRMED, major.**
  - With the current parameters at k_e = 10, the stationary state of the exactly averaged loop has L̄ = 0.1141, against L̄(μ*_τ) = 0.1112.
  - K is antisymmetric. Its cocycle residual is −3.5e-3 for n = 2, and K ≡ 0 for n = 1.
  - I did not reproduce R2's k_e = 3 figures, because my iteration did not converge there.
  - R2's first suggested fix (slow reference at fixed k_e) is REJECTED: the error from the shift x = x_m + e is O(1/k_e) regardless of ω_ref.
- **R1#3 + R2#9 (μ_0 in P2, the Dirac start, σ, τ > 0): CONFIRMED, minor.**
- **R1#4 + R2#9/#12 (f must be autonomous; the code uses f(x, 0.0)): CONFIRMED, minor.**
- **R1#5 + R2#9 (H(iii) mixes formal content and lacks the continuity step): CONFIRMED, minor.**
- **R1#6 + R2#10 (H(iv) gaps; the annealed runs are not instances): CONFIRMED, minor.** The remark that inf 𝒢 includes the leak is a clarification, not an error, since the text already writes inf 𝒢.
- **R1#7 + R2#5 (H(v) is a conditional inside the theorem): CONFIRMED, minor.**
- **R1#8 + R2#4 + part of R2#2 (Remark H.2 claims, including 'slows learning'): CONFIRMED, minor.**
- **R2#2, 'Evidence mixes the two settings / 5.5× / limit 0.194 = 62% of L̄(0)': REJECTED as stale.** The 22:37 text separates the two settings and quotes 4.2× and ≈ 0.11, which match the regenerated tables. The 0.194 figure was for the earlier G = diag(2, .5, .5).
- **R2#2, 'W₂ ≥ 1.08 at every τ shows the starts never equalise': UNCERTAIN.** No sampling floor for a W₂ between two 32-particle measures is given, and the losses do equalise at τ ≥ 0.004.
- **R1#10 + R2#7 (titles 'global optimality' and 'costs nothing'): CONFIRMED, minor.**
- **R1#9 + R2#11 (preamble 'stationary', the system the test runs, the Chizat–Bach separation hypothesis): CONFIRMED, cosmetic.**
- **R1#11 + R2#8 (tests): CONFIRMED, cosmetic.**
- **R1#12 + R2#6 (Remarks H.3 and H.4): CONFIRMED, cosmetic.**
- **R1#13 + R2#12 (bookkeeping and code):**
  - The [S6] reference to 'Remark H.5', the Assumption 1 bullet, the a_max notation, the LSI metric, the _sig overflow and the missing τ ≤ 0 guard: CONFIRMED, cosmetic.
  - 'RESULTS_PARADIGM.md has no Theorem H section': REJECTED as stale. The section was added at 22:37 with all six τ rows.
  - 'Stale results/langevin_global/table.md': REJECTED as stale. The file was removed, and the tables now compute their own reading.
- **New (mine): CONFIRMED, cosmetic.** The generated tables print '2τρ_τ ≥ 0.0e+00' and a hard-coded 'the observed convergence is far faster', even where their own computed reading says the starts did not equalise. I folded this into the Remark H.2 issue.

No statement about the averaged system (𝔄) in the H(i) identification, H(ii) or the H(iv) hypothesis is false. The defects are formal content presented as proved, an ill-posed conjecture, missing hypotheses, and overclaims in titles and remarks.

### H.1 [major] PARADIGM.md §6c Theorem H(i): sentence 'Moreover (𝔄) is the leading-order (Tikhonov + period-averaging) description ... checked numerically to 1e-7 by test_averaged_drift_is_the_mean_field_langevin_drift'; proof of (i), 'For the averaging statement'; §7 verification-map row 'H(i) identification'; §0 table row H; tests/test_langevin.py::test_averaged_drift_is_the_mean_field_langevin_drift

The step from the closed loop to (𝔄) sits inside a theorem whose status is 'proved', but it is only a formal substitution, and it is not labelled as such. The proof sets e = r_μ/k_e, evaluates the regressor at x_m instead of at x = x_m + e, and freezes μ over a period. It states no time-scale hypotheses (unlike Prop. G (H2a)) and gives no error order. The mean-field limit of the unprojected stochastic law (L_τ) is not established anywhere: Theorem D covers only the deterministic projected law (L). The word 'formal' appears only in review/REVIEW_LOG.md and LITERATURE_TRANSPORT.md, not in PARADIGM.md. The cited test imposes e = r/k_e at x = x_m and compares grad_g with finite differences of first_variation on the same quadrature. That is a check of an algebraic identity and of code consistency, and it cannot detect an averaging error.

**Argument.** Judge script avg.py: frozen μ with 8 random particles, the experiment's plant and network, G = diag(4,1,1). The comparison uses the true periodic error ė = r_μ(x_m+e) − k_e e with x = x_m + e.
- n = 1: the relative error of the period-averaged drift against (𝔄) is 4.2e-2, 1.3e-2, 5.7e-3, 2.7e-3 and 1.4e-3 at k_e = 10, 20, 40, 80 and 160; k_e·error tends to 0.22.
- n = 2 (ellipse): the error falls from 9.3e-2 to 5.8e-3; k_e·error = 0.93.
- e does not depend on G, so the error is the same for any gains and does not vanish as the gains go to 0.
- With the filter alone (x = x_m), the n = 1 error is second order (2.6e-2 down to 1.1e-4). So for n = 1 the first-order error comes from the shift x = x_m + e, not from ω_ref/k_e.

**Fix.** 1. In H(i), replace the 'Moreover ...' sentence with: 'Substituting e := r_μ(x_m)/k_e and x := x_m in the drift of (L_τ) and averaging over one reference period gives exactly −G∇V[μ] (an algebraic identity).'
2. Add 'Remark H.0 (formal derivation, not proved)'. Suggested content:
   - (𝔄) is formally the leading-order description of the mean-field closed loop under (L_τ), under two conditions.
   - Tikhonov condition: k_e ≫ max(ω_ref, Lip(f) + K_Φ/s_min, (F_ρ+K_Φ)/s_min), so that the shift x − x_m = O(1/k_e) is negligible on the kernel scale.
   - Period-averaging condition: T₀·sup‖G∇V[μ]‖ ≪ ℓ and τλ_G T₀ ≪ ℓ², where ℓ is the ϑ-scale of Φ̃.
   - The first corrections are O(1/k_e) relative to the drift, do not depend on the gains, and for n ≥ 2 contain a non-gradient interaction.
   - The exact formal limit is k_e → ∞ with G, k_eσ and k_eτ fixed, in time t/k_e.
   - The mean-field limit of (L_τ) itself is not proved here.
3. Change the verification-map row to: 'H(i) algebraic identity | with e := r/k_e, x := x_m: averaged drift = −G∇V (code gradients vs finite differences, 1e-7); does not test averaging'.
4. Change the §0 status of H to: 'proved for (𝔄); the reduction to (𝔄) is formal (Remark H.0); real loop: Conjecture H.6'.
5. Optionally, add a test that integrates the true periodic e at frozen μ and asserts that the drift error decreases like 1/k_e (k_e = 10, 20, 40, 80). Include an n = 2 orbit that encloses area.

### H.2 [major] PARADIGM.md §6c Conjecture H.6 ('with finite k_e and gains ... radius vanishes with the averaging parameters ... should follow from the exponential stability of (𝔄)'); §8 item 1 ('using H(ii) as the contraction'); Evidence paragraph

As stated the conjecture is ill-posed, and under its natural reading it is false.
- 'Averaging parameters' is never defined.
- At fixed k_e, let the gains go to 0, scaling G and σ together (which leaves μ*_τ unchanged). The limit is the exactly averaged system, whose drift uses the true periodic error e_μ and x = x_m + e_μ. That system is not (𝔄). Its stationary state is offset from μ*_τ by O(1/k_e). For n ≥ 2 its first correction is a non-reciprocal interaction, so it is not a Wasserstein gradient flow and [S3] does not apply.
- If instead 1/k_e → 0 at fixed σ and τ, μ*_τ itself tends to the prior N(ϑ⁰, (τ/σ)G), and the statement becomes degenerate.
- There is no growth or tube hypothesis on f and no bound on e(0), and the conjecture presupposes the mean-field limit of (L_τ).
- H(ii) gives decay of the F_τ-gap to the minimiser, not a contraction between solutions. Its only available ISS gain is 1/(τρ_τ) ≤ e^{α*/τ}/σ.
- All experiments use n = 1, where the non-gradient term vanishes identically, so the evidence cannot detect it.

**Argument.** (a) Stationary-state offset (judge script exact_fp.py). Current parameters: σ = 0.02, τ = 0.004, G = diag(4,1,1), k_e = 10, n = 1. The exactly averaged loop's stationary state was computed as a self-consistent Gibbs fixed point. This is exact for n = 1, because the frozen-μ drift is a ϑ-gradient; the residual is 1e-15. It has L̄ = 0.1141, against L̄(μ*_τ) = 0.1112 on the same grid (L1 distance 0.048). R2's particle runs, at R2's parameters, give 0.2016 against 0.1981 at k_e = 10.

(b) Non-gradient interaction (judge script kernel.py). K(ϑ,ϑ') = avg_t Φ̃(x_m;ϑ)ᵀ(d/dt)Φ̃(x_m;ϑ') is antisymmetric (|K+Kᵀ| ≤ 3e-13). For an n = 2 ellipse, K01+K12−K02 = −3.5e-3, with entries between 1.4e-3 and 7.3e-3. So the first Tikhonov correction (1/k_e²)G∇_ϑ∫K dμ is not of the gradient form a(ϑ)−a(ϑ'). For n = 1, K ≡ 0 (1e-13). Reversing the time direction of the reference flips K but leaves the x-shift terms unchanged. Hence for at least one orientation of the orbit the total first-order correction is non-gradient.

(c) Degenerate reading. Since k_eF_τ = L̄ + (k_eσ/2)M + k_eτH, μ*_τ tends to N(ϑ⁰, (τ/σ)G) as k_e → ∞ at fixed σ, τ.

(d) R2's suggested fix of a slow reference (T₀ → ∞) at fixed k_e fails. The x-shift error is O(1/k_e) independently of ω_ref (avg.py, n = 1).

**Fix.** Replace H.6 with:
'Conjecture H.6. Let f = f(x) satisfy the tube condition of Theorem D(b). Fix G ≻ 0 and σ₀, τ₀ > 0, and set σ := σ₀/k_e, τ := τ₀/k_e. Then μ*_τ, the minimiser of L̄ + (σ₀/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ + τ₀H, does not depend on k_e. In slow time s = t/k_e, (𝔄) is the formal limit of the closed loop. For the mean-field limit of the closed loop under (L_τ), every μ_0 in P2 and every e(0), there are S(μ_0) and r(k_e) → 0 such that W_{G⁻¹}(μ_t, μ*_τ) ≤ r(k_e) for all t ≥ k_e·S.'

Also add: 'At fixed k_e the analogous statement is false. The gains → 0 limit is the exactly averaged system (true periodic e_μ, x = x_m + e_μ). Its stationary state is offset from μ*_τ by O(1/k_e) (L̄ 0.1141 vs 0.1112 at k_e = 10), and the system is non-gradient for n ≥ 2.'

In §8 item 1, replace 'using H(ii) as the contraction' with: 'using the F_τ-gap as an ISS-Lyapunov function after a near-identity averaging transformation. For a drift perturbation δb, d/dt(F_τ − F*) ≤ −τρ_τ(F_τ − F*) + ½‖δb‖²_∞, hence limsup W ≤ ‖δb‖_∞/(τρ_τ). In the scaling above, ‖δb‖ = O(k_e⁻²) and τρ_τ = O(k_e⁻¹), so r(k_e) = O(1/k_e), with a constant of order e^{2K_Φ(F_ρ+K_Φ)/τ₀}/σ₀.'

In the Evidence, note that all runs are n = 1, where the non-gradient term vanishes, and add an n = 2 experiment.

### H.3 [minor] PARADIGM.md §6c Theorem H(ii) ('For every μ_0 with F_τ(μ_0) < ∞'); preamble 'Entropy removes the obstruction'

[S3] Theorem 3.2 assumes μ_0 in P2 and F_τ(μ_0) < ∞. The statement drops μ_0 in P2, and outside P2 the expression 𝒢 + τH can be ∞ − ∞. It also never states σ > 0 and τ > 0. More importantly, H(μ_0) = +∞ for every measure without a density. So the motivating collapsed start δ_ϑ̄, and every m-particle empirical start, lie outside (ii) as stated. That is exactly the case the preamble says entropy handles.

**Argument.** H(δ_ϑ̄) = +∞, so F_τ(δ_ϑ̄) = +∞: Cor. 3.3 becomes vacuous and the hypothesis of Thm 3.2 fails. With σ = 0, V_L is bounded, so e^{−V/τ} is not integrable.

**Fix.** 1. Restate (ii) as: 'Assume σ, τ > 0. For every μ_0 in P2(ℝ^p) with F_τ(μ_0) < ∞ ...'.
2. Write F_τ = (1/k_e)L̄ + τH(·|γ) + const with γ := N(ϑ⁰, (τ/σ)G). Then F_τ is well defined on all of P(ℝ^p) and finite only on P2^a.
3. Add the following, with its proof:
   - Statement: 'If μ_0 is only in P2 (for example a Dirac or empirical measure), then F_τ(μ_{t₀}) < ∞ for every t₀ > 0, so (ii) holds from t₀ on.'
   - Setup: in ζ-coordinates, compare the particle SDE with the OU process dY = −σ(Y−ζ⁰)dt + √(2τ)dB started from the same law. The two drifts differ by b = −(1/k_e)G^{1/2}∇V_L[μ_t], with ‖b‖ ≤ B := (1/k_e)λ_G^{1/2}(F_ρ+K_Φ)·sup‖∇_ϑΦ̃‖.
   - Entropy bound: Girsanov and data processing give H(μ_{t₀}|OU_{t₀}) ≤ B²t₀/(4τ). OU_{t₀} has density at most (2π(τ/σ)(1−e^{−2σt₀}))^{−p/2}, so H(μ_{t₀}) < ∞.
   - Moments: μ_{t₀} is in P2 because the drift grows linearly.

### H.4 [minor] PARADIGM.md §6c 'Averaged mean-field system' paragraph (ρ, F_ρ := sup_{supp ρ}‖f‖, r_μ := f − ∫Φ̃dμ); src/core/langevin.py occupation_loss and first_variation (plant.f(x, 0.0))

r_μ, F_ρ and the occupation measure of x_m all presuppose an autonomous plant, but Theorem H never says so. Standing hypothesis (S1) allows f(x,t), while Proposition G does state 'f = f(x) is time-invariant'. For a time-varying f, the period average is not a functional of ρ = law(x_m), and V_L as defined is meaningless. The code silently evaluates f(x, 0.0).

**Argument.** Take f(x,t) = sin(t)·x with x_m T₀-periodic. The averaged drift contains ∫f(x_m(t),t)Φ̃(x_m(t);ϑ)dt, which is not of the form ∫f(x)Φ̃(x;ϑ)dρ(x) for any f(x).

**Fix.** Add to the 'Averaged mean-field system' paragraph: 'Assume f = f(x) is continuous and time-invariant (as in Proposition G).'

Alternatively, to allow a T₀-periodic f(x,t), set:
- ρ := law of (x_m(t), t mod T₀);
- F_ρ := sup over supp ρ of ‖f(x,t)‖;
- r_μ(x,t) := f(x,t) − ∫Φ̃(x;·)dμ.

In langevin.py, either document 'autonomous plants only' in occupation_loss and first_variation, or pass the time samples to plant.f.

### H.5 [minor] PARADIGM.md §6c Theorem H(iii) and its proof '(iii) Combine (i) and (ii)'

Part (iii) mixes formal and rigorous content inside the theorem. (𝔄) is an equation for μ alone, so 'along (𝔄) the quasi-steady tracking error e ≈ r/k_e' and '⟨‖e‖²⟩_ρ ≈ 2L̄/k_e²' are formal. The rigorous claim L̄(μ_t) → L̄(μ*_τ) needs W-continuity of L̄, which the proof does not supply.

**Argument.** 'Combine (i) and (ii)' yields only that W_{G⁻¹}(μ_t, μ*_τ) → 0 and that the F_τ-gap → 0. Neither implies convergence of L̄ without a continuity estimate.

**Fix.** Replace (iii) with: 'Let ē_t := r_{μ_t}/k_e. Then ⟨‖ē_t‖²⟩_ρ = 2L̄(μ_t)/k_e², and |L̄(μ_t) − L̄(μ*_τ)| ≤ (F_ρ+K_Φ)L_Φ·W_{G⁻¹}(μ_t, μ*_τ) ≤ (F_ρ+K_Φ)L_Φ(2/(τρ_τ))^{1/2}e^{−τρ_τt}(F_τ(μ_0) − F_τ(μ*_τ))^{1/2}, where L_Φ := sup_{x,ϑ,‖u‖=1}‖G^{1/2}∇_ϑ(uᵀΦ̃)‖.'

Proof:
- L̄(μ)−L̄(ν) = ½∫(r_μ−r_ν)ᵀ(r_μ+r_ν)dρ.
- ‖r_μ+r_ν‖ ≤ 2(F_ρ+K_Φ).
- ‖r_μ−r_ν‖ ≤ L_Φ·W_{1,G⁻¹} ≤ L_Φ·W_{G⁻¹}.
- Conclude with Cor. 3.3.

Move the link 'e ≈ ē' into the formal Remark H.0.

### H.6 [minor] PARADIGM.md §6c Theorem H(iv); Remark H.2 ('A hot-start annealed schedule reaches 3.2e-2 from both starts')

Part (iv) has three gaps.
- The schedule τ_t = α/log t is undefined or ≤ 0 for t ≤ 1. [S3] Thm 4.1 needs τ_t smooth and decreasing, equal to α/log t only for large t, and μ_0 in P2 with F_{τ₀}(μ_0) < ∞. None of this is stated.
- 'Exactly the bound in (ii)' is inaccurate. (σ/τ)e^{−α*/τ} ≥ C₀e^{−α*/τ} holds with C₀ = σ/τ₀ only for τ ≤ τ₀. That is enough, because the proof in [S3] uses the LSI only along τ_t ≤ τ₀. The review log says this, but PARADIGM.md does not.
- The 'annealed' experiments are far from α > α*, so they are not instances of (iv), and the text does not say so.

The point that inf 𝒢 is the σ-regularised optimum is only a clarification, since the text already writes inf 𝒢.

**Argument.** With F_ρ = 1.2469, K_Φ = 4 and k_e = 10, α* = 4.198. The schedules actually used are:
- α = τ₁·log(e+20) = 1.56e-3 (exp_langevin_global, τ₁ = 5e-4);
- α = 1.25e-2 (exp_langevin_global, τ₁ = 0.004);
- α = 0.008·log(e+10) = 2.03e-2 (exp_langevin_tradeoff).

With α > α*, τ_t reaches 0.004 only when log t ≈ 1050.

**Fix.** Restate (iv): 'Let μ_0 be in P2 with F_{τ₀}(μ_0) < ∞, and let (τ_t) be smooth and decreasing with τ_t = α/log t for t ≥ t₁ > 1, where α > α* := 2K_Φ(F_ρ+K_Φ)/k_e. Then 𝒢(μ_t) − inf 𝒢 ≤ C' log log t/log t, where inf 𝒢 is the σ-regularised optimum. This is [S3] Thm 4.1. Its hypotheses hold: 𝒢 ≥ 0; Assumptions 1–3 hold for every τ > 0; and for τ in (0, τ₀], (ii) gives ρ_τ ≥ C₀e^{−α*/τ} with C₀ = σ/τ₀, which is all its proof uses.'

Add to Remark H.2: 'The annealed runs use α between about 2e-3 and 2e-2, far below α* = 4.2, so they are not instances of (iv).'

### H.7 [minor] PARADIGM.md §6c Theorem H(v) and proof '(iv) and (v) are the cited theorems'; Remark H.1 ('At finite m the fluctuation of f̂ is O(K_Φ/√m)')

Part (v) is a conditional whose hypotheses the text admits are unchecked ('must still be checked verbatim'). Yet it sits inside a theorem whose status is 'proved', and its proof is a citation. 'The m-particle averaged system' is never defined. Remark H.1's finite-m fluctuation bound holds uniformly in time only if (v) holds; for interacting particles, the finite-horizon bound grows with t.

**Argument.** As written, (v) says nothing about (𝔄) until the hypotheses of [S5] are verified. LITERATURE_TRANSPORT.md, open question 4, concedes this.

**Fix.** Move (v) out of the theorem into 'Remark H.5 (finite m, conditional)', as follows.
- Define the system: dϑ_i = −G∇V[μ^m_t](ϑ_i)dt + √(2τ)G^{1/2}dB_i, with μ^m_t := (1/m)Σδ_{ϑ_i}.
- State uniform-in-time propagation of chaos as conditional on the hypotheses of [S5].
- List the hypotheses already available: convexity of 𝒢; a uniform LSI, from (ii); and a bounded, smooth second variation (1/k_e)∫Φ̃(x;ϑ)ᵀΦ̃(x;ϑ')dρ.
- List those that remain to be checked.

In Remark H.1, write: 'O(K_Φ/√m) at each fixed t (finite-horizon propagation of chaos); uniformly in time under Remark H.5'.

### H.8 [minor] PARADIGM.md §6c Remark H.2 (k_e bullet 'it also slows learning'; 'underflows to 0 at τ = 5·10⁻⁴, while escape is observed within hundreds of seconds'; the τ table; the escape-time fit); experiments/exp_langevin_global.py table header

Several statements in Remark H.2 go beyond the data or the theorem.
(a) The vacuous bound at τ = 5e-4 is set against 'escape within hundreds of seconds'. But at τ = 5e-4 one seed escapes at 1050 s and the other does not escape within 1500 s, and at τ = 1e-3 escape takes 1005–1480 s.
(b) The displayed table omits the τ = 5e-4 (censored) and τ = 2e-3 rows, which RESULTS_PARADIGM.md does include.
(c) Two things confound the log(1/τ)-versus-Arrhenius comparison, which the text already calls 'indicative only':
- The escape threshold, twice the spread run's final loss at the same τ, grows with τ, and the escape time includes convergence after escape.
- The censored τ = 5e-4 point is dropped from the fit without comment.
(d) 'Higher feedback gain ... also slows learning' omits the main effect. Since k_eF_τ = L̄ + (k_eσ/2)M + k_eτH, both μ*_τ and ρ_τ depend only on (k_eσ, k_eτ). Raising k_e at fixed (σ, τ) therefore raises the effective leak and temperature. This improves the guaranteed rate 2σe^{−α*/τ} (α* ∝ 1/k_e) and pushes μ*_τ toward the prior.
(e) The generated tables print '2τρ_τ ≥ 0.0e+00' and a hard-coded 'the observed convergence is far faster'. This text appears even in table_sigma0.002_tau0.0005.md, whose own computed reading says the starts did not equalise.

**Argument.** results/langevin_tradeoff/table.md:
- τ = 0.0005: escape '∞, 1050'.
- τ = 0.001: escape '1480, 1005'.
- The spread run's final loss rises from 1.09e-2 (τ = 5e-4) to 6.53e-2 (τ = 8e-3), so the escape threshold moves from 0.022 to 0.13.

The k_e scaling identity is exact: multiply F_τ by k_e. R2's grid computation (G = diag(2, .5, .5)) gives L̄(μ*_τ) = 0.194, 0.283 and 0.303 at k_e = 10, 30 and 100, against L̄(0) = 0.311.

**Fix.** (a) Replace the sentence with: 'The guaranteed rate is 2τρ_τ ≥ e^{−8400} s⁻¹ at (σ, τ) = (0.002, 5·10⁻⁴) and e^{−1053} s⁻¹ at (0.02, 0.004). Observed escape of the collapsed start takes 295–1480 s for τ ≥ 10⁻³; at τ = 5·10⁻⁴, one of two seeds does not escape within 1500 s.'
(b) Show all six τ rows.
(c) Add: 'the escape threshold grows with τ, and the τ = 5·10⁻⁴ point is censored and excluded from the fit'. Alternatively, define escape by a τ-independent criterion, such as exit from a W-neighbourhood of the collapsed configuration.
(d) Replace 'but it also slows learning' with: 'Since k_eF_τ = L̄ + (k_eσ/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ + k_eτH, μ*_τ and ρ_τ depend only on (k_eσ, k_eτ). Raising k_e at fixed (σ, τ) improves ρ_τ and the guaranteed rate, but moves μ*_τ toward N(ϑ⁰, (τ/σ)G), which increases the regularisation bias. Meanwhile the loss-driven drift shrinks like 1/k_e.'
(e) In exp_langevin_global.py, print log(2τρ_τ) instead of the underflowed value, remove the hard-coded 'the observed convergence is far faster', and regenerate both tables.

### H.9 [minor] PARADIGM.md §6c title ('global optimality of kernel placement'); Remark H.1 title ('dither costs nothing in the mean-field limit')

The theorem proves convergence to the unique minimiser of the regularised free energy F_τ, not optimal kernel placement for L̄, and it gives no bound on L̄(μ*_τ) − inf L̄. Remark H.1's title contradicts the remark's own text ('the only effect of the dither is the entropic bias τH') and the project's data.

**Argument.** Judge grid computation (gibbs.py), current parameters σ = 0.02, τ = 0.004, G = diag(4,1,1): L̄(μ*_τ) = 0.111, against the offline unregularised optimum of 7.0e-6 and L̄(0) = 0.311. In results/langevin_tradeoff, the spread-start loss floor rises from 6.96e-3 (τ = 0) to 6.53e-2 (τ = 0.008), and the tracking RMS rises from 0.0116 to 0.0363.

**Fix.** 1. Retitle §6c 'Theorem H — Langevin adaptive control and global convergence to the regularised optimum'.
2. Retitle Remark H.1 'Dither adds no output noise in the mean-field limit; its cost is the entropic bias'.
3. Add the exact bias bound. With γ := N(ϑ⁰, (τ/σ)G), F_τ = (1/k_e)L̄ + τH(·|γ) + const, hence L̄(μ*_τ) ≤ inf over ν of {L̄(ν) + k_eτH(ν|γ)}.
4. Quote L̄(μ*_τ) = 0.111 against 7e-6 for the experiment.

### H.10 [cosmetic] PARADIGM.md §6c preamble ('the symmetric "collapsed" configuration is stationary and never separates without noise. This is checked by test_noiseless_...'); §7 row 'H, need for entropy'

The preamble misdescribes the collapsed configuration and the test that checks it.
- 'Stationary' is false. The collapsed configuration is invariant (it stays a Dirac mass), but it moves.
- The cited test runs the m = 6 real closed loop, not 'the averaged mean-field limit'. Invariance holds in both systems, by permutation symmetry and uniqueness of solutions.
- A Dirac start also violates the support-separation hypothesis of Chizat–Bach Thm 3.5. So it illustrates that hypothesis, not the failure of the boundary conditions 3.4(iii) discussed in the preceding sentence.

**Argument.** misc.py: with τ = 0 and the test's initial state, the rows stay exactly equal at every snapshot, while the common particle moves by (0.109, −0.412, 0.207) in 5 s. cb.txt, Thm 3.5: 'the support of µ0 is contained in [−r0, r0] × Θ and separates {−r0} × Θ from {r0} × Θ'.

**Fix.** Replace the sentences with: 'The symmetric collapsed configuration is invariant: it remains a (moving) Dirac mass and never separates without noise, by permutation symmetry and uniqueness of solutions. This is checked on the m = 6 real closed loop by test_noiseless_symmetric_initialisation_never_separates_but_dither_does. A Dirac start also violates the separation hypothesis on supp μ_0 of Chizat–Bach Thm 3.5, which is an obstruction distinct from the boundary conditions.' Update the §7 row to say 'finite-m real loop'.

### H.11 [cosmetic] tests/test_langevin.py::test_noiseless_symmetric_initialisation_never_separates_but_dither_does (lines 79–85) and ::test_lsi_lower_bound_formula; §7 rows 'H(ii) constant' and 'H(i) identification'

The tests check less than the verification map claims, and one is fragile.
- test_lsi_lower_bound_formula recomputes the same expression as the code. It is a transcription test, yet the verification map lists it as a check of H(ii).
- The symmetry test asserts np.std(...) == 0.0. That value can be nonzero for exactly identical rows, because the mean is rounded, so the test passes only because of the particular final values.
- Its τ = 0.01 half is trivially true, since noise always creates spread.
- No test checks the |V_L| bound, an LSI or Poincaré constant, or convergence of (𝔄).

**Argument.** np.std(np.full(6, 0.1)) = 1.39e-17. For 10,000 random values, 28% give a nonzero np.std over six identical copies. In the test's own τ = 0 run, np.std is nonzero at an intermediate snapshot although the rows are exactly equal (misc.py).

**Fix.** 1. For τ = 0, assert all(np.all(S == S[0]) for S in res['Th']).
2. For τ > 0, assert separation early (for example after 0.1 s), or drop that half.
3. Label the §7 row 'H(ii) constant: formula transcription'.
4. Add these tests:
   - (a) a random or adversarial maximisation asserting |V_L[μ](ϑ)| ≤ K_Φ(F_ρ+K_Φ);
   - (b) the frozen-μ averaging-error test described in the H(i) issue;
   - (c) optionally, a grid Gibbs fixed point μ*_τ compared with a particle simulation of (𝔄) from two starts at moderate τ.

### H.12 [cosmetic] PARADIGM.md §6c Remarks H.3 and H.4

Remark H.3: projecting the drift with the Pomet–Praly operator while keeping additive noise does not give reflected dynamics, because the noise still leaves Θ_K. The Neumann/Bakry–Émery route needs a Skorokhod reflection. On a bounded convex body, the uniform measure already satisfies an LSI of order 1/diam², so σ-confinement would not be needed there. That limits Remark H.4's claim that σ-modification is 'the reason the global result exists' to the saturated design on ℝ^p.

**Argument.** The SDE dϑ = Proj_K(ϑ, b)dt + √(2τ)G^{1/2}dB exits Θ_K with positive probability. With reflection and σ = 0, Holley–Stroock applied to the uniform measure of Θ_K gives LSI(c·diam⁻²·e^{−2K_Φ(F_ρ+K_Φ)/(k_eτ)}).

**Fix.** Remark H.3: 'Reflected (Skorokhod) dynamics on the convex body Θ_K, not the Pomet–Praly drift projection (which leaves the noise unconfined), satisfy an LSI with constant of order 1/diam², even for σ = 0 ...'.

Remark H.4: 'In the saturated design on ℝ^p, σ-modification supplies the confinement that makes e^{−V/τ} integrable and log-Sobolev. Under reflection on Θ_K, compactness would play that role.'

### H.13 [cosmetic] LITERATURE_TRANSPORT.md [S6] (line 23); PARADIGM.md §6c numbering, Assumption 1 bullet, design paragraph (a_max), H(ii) LSI statement; src/core/langevin.py _sig and lsi_constant_lower_bound

Bookkeeping and code nits:
(1) [S6] cites 'PARADIGM.md Remark H.5', which does not exist; the remarks run H.1–H.4, followed by Conjecture H.6.
(2) The Assumption 1 bullet omits the uniform bound ‖r_μ‖ ≤ F_ρ+K_Φ on supp ρ for every μ. That bound is what makes the ζ-Lipschitz constant independent of μ.
(3) In §1, a_max bounds ‖a‖; in §6c it bounds each coordinate, so K_Φ = a_max√n.
(4) The text never says that the LSI is for the whitened measures, i.e. with Fisher information in the G-metric in ϑ-coordinates.
(5) In Assumption 3 only, F_ρ can be replaced by ∫‖f‖dρ (0.69 instead of 1.25 here).
(6) _sig overflows with a RuntimeWarning for κ ≪ 0, although the value is correct.
(7) lsi_constant_lower_bound has no τ > 0 guard and returns nan at τ = 0.

**Argument.** Each point can be checked by reading the files. misc.py reproduces the RuntimeWarning and the nan.

**Fix.** 1. Point [S6] to the §6c preamble, or renumber.
2. In the Assumption 1 bullet, add 'since ‖r_μ‖ ≤ F_ρ+K_Φ on supp ρ for every μ in P(ℝ^p)'.
3. Use a separate symbol for the per-coordinate bound, or say explicitly that K_Φ = a_max√n uses a different convention from §1.
4. State 'LSI for the whitened measures (Fisher information in the G-metric)'.
5. Optionally use ∫‖f‖dρ in the Holley–Stroock oscillation.
6. Use scipy.special.expit for _sig.
7. Raise ValueError for τ ≤ 0 in lsi_constant_lower_bound.

Refuter verdicts: holds_with_fixes; holds_with_fixes

## G — verdict: holds_with_fixes

**Verdict: holds with fixes.** Part (a) holds as a formal statement, and its factors m/k_e (particle form) and 1/k_e (measure form) are right. Part (b) has two major problems: its hypothesis (H2b) is not sufficient, and its measure-form/[S2] sentence is false. Both are fixable by adding a hypothesis and rewording. No project files were edited. I reviewed the current files: PARADIGM.md and RESULTS_PARADIGM.md dated 22:37, table_satisfying.md dated 22:32.

**Targeted questions**
- **Factors m/k_e: correct.** Under law (L), ȧ_i = g_a·e·φ_i, with ∇_{a_i}L̄ = −(1/m)∫rφ_i dρ and e ≈ r/k_e. This gives ϑ̇_i = −(m/k_e)G∇_{ϑ_i}L̄ = −(1/k_e)G∇(δL̄/δμ)(ϑ_i). I checked this two ways:
  - On the stored violating run, finite-difference ∇L̄ against period-averaged joint velocity over 119 periods: cosine ≥ 0.9985 (G⁻¹ metric), magnitude ratio 0.93–1.00.
  - With parameters frozen, at the satisfying geometry: ratio 0.956 / 0.996 / 0.9996 at k_e = 20 / 80 / 320, cosine ≥ 0.9999.
- **Is (a) the noiseless case of H(i)?** Not literally. H(i) uses the saturated features Φ̃ in (α,γ,κ) coordinates, which gives a different metric. (a) is the same averaging computation with Φ in (a,c,s).
- **Consistency with the experiment.** The script uses the document's formulas: v_pred = −(m/k_e)δ∇L*, w_min = g_a·m·m₀/k_e, geo = (mδ/k_e)‖∇²L*‖. The numbers in table_violating.md match the text.

**Decisions on each reviewer issue (duplicates merged)**
- R1-1 + R2-4, [S2]/measure form: CONFIRMED, **major**.
- R2-1, slow-manifold slope: CONFIRMED, **major**. I reproduced it independently, both on the averaged flow and on the closed loop with the project code.
- R1-2, chain of theorems does not reach the O(1/δ) horizon: CONFIRMED, downgraded to **minor**. The proposition is explicitly called a formal reduction; what is wrong is the claim that the derivation chains standard theorems.
- R1-3 + R2-2, (H2a) mis-scaled: CONFIRMED, minor.
- R1-4 + R2-6, relation to H(i): CONFIRMED, minor.
- R1-5 + R2-5, missing hypotheses (stay in 𝒢, a* inside the weight box): CONFIRMED, minor.
- R1-10 + R2-3, Status/novelty wording: CONFIRMED, minor.
- R1-9, o(1) limit unspecified: CONFIRMED, minor.
- R1-7 / R2-8, Γ missing from (H2b): CONFIRMED, minor. R1-7's "O(δ)" point is merged into the slope issue, and its offset point into the o(1) issue.
- R1-6, "L̄ non-increasing" holds only for the averaged flow: CONFIRMED, cosmetic.
- R2-7, explanation of the 0.94–0.95 ratio: CONFIRMED, cosmetic, with a refined attribution.
- R1-11 + R2-8 rest, table and numerics labels: CONFIRMED, cosmetic.
- R1-8, "evidence for (b) does not exist": **REJECTED** as stale. table_satisfying.md and RESULTS §G now exist, their numbers match the text, and R2 reproduced the table byte for byte.
- No issue is UNCERTAIN.

**Key evidence for the slope issue (my own runs)**
- Setup: m = 1, x_m = 0.2 sin t, f = 2exp(−(x−1)²/0.5), k_e = g_a = 20, g_c = g_s = 0.47.
- Hypotheses hold at least as well as in the document's positive configuration:
  - (H2b) ratio 76–178 along the trajectory (document: 77.5).
  - ω_ref / fastest joint rate = 12.6 (document: 7.3).
  - But the slope parameter S = 2.0–2.9.
- Displacement Δ(c,s) over 232 s:

| Flow | Δc | Δs |
|---|---|---|
| Averaged flow (a) | −0.0036 | −0.0281 |
| Project closed loop | −0.0031 | −0.0297 |
| Prediction of (b) | +0.0170 | −0.0680 |
| Induced-metric flow | −0.0048 | −0.0259 |

- The closed loop moves 0.43× as far as (b) predicts, with cosine 0.94 and the centre moving the opposite way. The induced-metric flow matches the averaged flow (cosine 0.9985).
- The error tracks S, not the (H2b) ratio: at δ = 0.047 (S = 0.2–0.3) the closed loop moves 0.82× of (b); at δ = 0.0047 (S = 0.02) it moves 0.98×.
- The document's two configurations cannot tell (H2b) and S apart: S = 0.006 in the satisfying one and 27 → 5.6 in the violating one.

**Strongest correct form of (b)**
- Hypotheses: (H0), (H1), (H2a′), a Γ-weighted (H2b), and a new (H2c): S := (δ/g_a)·sup‖Da*·Γ^{1/2}‖² ≪ 1.
- Statement: for t ≤ τ₁k_e/(mδ), (c,s)(t) = (c̄,s̄)(δt) + o(1), as (ω_eff/k_e, ε_avg, S, 1/ratio) → 0. Here (c̄,s̄) is the stated gradient flow of L*, and the statement holds at fixed m.
- Without (H2c), the first-order slow flow is instead −(mδ/k_e)(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*.
- There is no measure-level form, and the flow is not [S2]'s.

**Scripts** (in <scratchpad>\G_judge\):
- slope_avg.py, slope_closed.py — the slope counterexample (averaged flows, project closed loop)
- doc_configs_S.py — S and rate ratios for all configurations
- frozen.py, h2a_check.py — frozen-parameter checks of the factors and of (H2a)
- check_a_violating.py — direct check of (a) on the stored violating run
- reduced_violating.py — reduced flow from the violating geometry (growth of a*)
- window_bias.py — forward-difference bias on the satisfying run

### G.1 [major] PARADIGM.md §6b: (H2b), Proposition G(b), Derivation step 3, and Numerics bullet 'Three scales need (H2b)'; results/RESULTS_PARADIGM.md §G, last sentence

CONFIRMED (R2 major; also absorbs R1's point that step 3's O(δ) misstates the small parameter).

(H2b) compares the slowest weight rate with the reduced-flow rate (mδ/k_e)‖∇²L*‖, which is a Schur complement. It does not control the slope of the slow manifold a*(c,s). Because L̄ = L* + ½(a−a*)ᵀM(a−a*):
- the weights lag behind, with a − a* ≈ (δ/g_a)M⁻¹Da*Γ∇L*;
- that lag feeds back through ∇_ξL̄ = ∇L* − Da*ᵀM(a−a*).

The first-order slow flow is therefore ξ̇ = −(mδ/k_e)(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*. This is the gradient flow of L* in the metric induced on the graph of a*, and it couples all kernels. It reduces to (b) only when S := (δ/g_a)‖Da*Γ^{1/2}‖² ≪ 1.

(H2b) does not imply S ≪ 1. The product (ratio)·S = m₀‖Da*Γ^{1/2}‖²/‖Γ^{1/2}∇²L*Γ^{1/2}‖ is fixed by the problem, and it is large whenever refitting the weights largely compensates a change of geometry. So the RESULTS sentence 'G(b) holds exactly when its hypothesis on the slowest weight mode holds' is false.

**Argument.** Judge's own runs (slope_avg.py, slope_closed.py).

Setup: m = 1, n = 1, x_m = 0.2 sin t, f = 2exp(−(x−1)²/0.5), k_e = 20, g_a = 20, g_c = g_s = 0.47, start at a = a*(1.75, 0.9). The projection is never active.

Hypotheses along the trajectory:
- (H2b) ratio 76–178 (the document's positive configuration has 77.5).
- ω_ref / fastest joint rate = 12.6 (document: 7.3).
- m₀ = 0.021–0.027.
- S = 2.0–2.9.

Displacement Δ(c,s) over 232 s:
- averaged flow (a): (−0.0036, −0.0281);
- project closed loop (TransportController + rk4_closed_loop, period means): (−0.0031, −0.0297);
- prediction of (b): (+0.0170, −0.0680). The closed loop moves 0.43× as far, cosine 0.94, and the centre moves the opposite way;
- induced-metric flow: (−0.0048, −0.0259), cosine 0.9985 with (a).

Scaling with δ:
- δ = 0.047: ratio 650–1780, S = 0.20–0.32, the averaged flow moves 0.82× of (b).
- δ = 0.0047: S = 0.02, 0.98×.

The document's configurations: satisfying S = 0.006; violating S = 27 → 5.6. The violating one fails both conditions, so it cannot isolate (H2b).

**Fix.** 1. Add (H2c), slow-manifold flatness: S := (δ/g_a)·sup_{𝒢′}‖D_{(c,s)}a*·Γ^{1/2}‖² ≪ 1. Because Da* = −M⁻¹∂²_{a,(c,s)}L̄, a checkable sufficient condition is δ·max(γ_c,γ_s)·‖∂²_{a,(c,s)}L̄‖² ≪ g_a m₀².
2. Add after (b): 'Without (H2c) the first-order slow flow is d(c,s)/dt = −(mδ/k_e)(Γ⁻¹ + (δ/g_a)Da*ᵀDa*)⁻¹∇L*, the gradient flow of L* in the metric induced on the graph of a*; the displayed flow is its S → 0 limit.'
3. In step 3, replace 'a = a*(c,s) + O(δ)' with 'a − a* = (δ/g_a)M⁻¹Da*Γ∇L* + h.o.t. = O(δ‖Da*‖‖Γ∇L*‖/(g_a m₀)) after the boundary layer'.
4. Make exp_variable_projection.py report S next to the (H2b) ratio.
5. Replace the RESULTS sentence with: 'G(b) holds in the satisfying configuration, where (H2b) and S ≪ 1 both hold; the violating configuration fails both (ratio 0.1, S = 27 → 5.6).' Rename the Numerics bullet to 'Three scales fail when (H2b) and (H2c) fail'.

### G.2 [major] PARADIGM.md §6b, Proposition G(b), last sentence ('In measure form, the geometry marginal follows the Wasserstein gradient flow of the variable-projection functional of [S2]'); LITERATURE_TRANSPORT.md §G; novelty table row G

CONFIRMED (R1 major and R2 minor, merged). Rated major because this is a false mathematical claim inside the proposition statement, and the 'bridge to [S2]' rests on it.

(i) [S2] (arXiv:2504.18208) studies the ridge-regularised functional ℒ^λ(μ) = min_u (1/2λ)‖Φ_μu − Y‖² + ‖u‖²_{L²(μ)}. Its λ→0 limit ℒ⁰ is a min-norm interpolation functional: +∞ off the range, χ²-type in the teacher–student case, and its flow is a weighted ultra-fast diffusion. The unregularised residual L* is not either of these.
(ii) L* has no Wasserstein gradient flow on P(geometries). See the argument below.
(iii) There is no mean-field version under (H2b), because m₀ ≤ 1/m².
(iv) When S is not small, the slow flow is not even particle-wise.

The fixed-m particle formula, including its factor m, is correct.

**Argument.** - Support dependence: as a functional of the geometry marginal ν, L*(ν) = inf_u ½‖f − ∫uφdν‖²_ρ depends only on supp ν, because the free weights absorb the masses.
- Discontinuity: take ν_ε = (1−ε)ν + εδ_ξ. Then W(ν_ε,ν) → 0, but L*(ν_ε) = L*(supp ν ∪ {ξ}) < L*(ν) whenever φ(·;ξ) is not ρ-orthogonal to the optimal residual. So L* is not W-continuous, has no first variation off the support, and has no Wasserstein gradient flow.
- Mean field: since 0 < φ ≤ 1, λ_min((1/m²)∫φφᵀdρ) ≤ 1/m². The slowest weight rate is therefore ≤ g_a/(k_e m), and the minimum eigenvalue of a Gaussian Gram matrix collapses quickly in m. Both reviewers give numbers: (H2b) ratio 4.7e3, 64, 6.2, 0.9 for m = 2, 4, 8, 12.
- The [S2] functional form and its λ→0 limit were checked against the arXiv abstract and HTML.

**Fix.** Replace the sentence with:

'For fixed m, on configurations of m distinct geometries in 𝒢, the reduced flow is the gradient flow of L*/k_e in the metric (1/m)Σ_i‖·‖²_{Γ⁻¹}, which is locally the Wasserstein metric restricted to uniform m-atom measures. L* is the unregularised profile loss. It is the λ→0 limit of λℒ^λ of [S2] (time rescaled by λ), not [S2]'s λ→0 functional ℒ⁰ (min-norm interpolation, whose flow is an ultra-fast diffusion). As a functional of the geometry marginal it depends only on the support, so the statement has no mean-field form, consistent with m₀ ≤ 1/m².'

Correct LITERATURE_TRANSPORT.md §G the same way. In the novelty table, change 'a bridge to [S1], [S2]' to 'a fixed-m, closed-loop analogue of the two-timescale/VarPro training of [S1], [S2]'.

### G.3 [minor] §6b Derivation (opening sentence and steps 1–3) and 'on time intervals of length O(1/δ)' in (b); LITERATURE_TRANSPORT.md §G 'Strongest objection'

CONFIRMED (R1, downgraded from major). The proposition is explicitly a formal reduction. What is wrong is the claim that the derivation chains three standard theorems: their horizons do not compose to O(1/δ).
- Khalil Thm 11.1 holds on a fixed interval [t₀,t₁].
- Thm 10.4(1) gives O(ε_avg) accuracy only on [0, b/ε_avg].
- Thm 10.4(2–3) and Thm 11.2 need an exponentially stable equilibrium. The joint averaged system has none: the geometry is slow and neutral.
- Step 3 is Tikhonov in slow time τ = δt on [0, τ₁]. It needs the averaged system to be valid up to t = τ₁/δ, which is much longer than b/ε_avg, and step 2 does not certify that.

Step 1's remainder O(ω_ref/k_e²) is also mis-scaled; see the (H2a) issue.

**Argument.** This follows from the standard statements of the theorems. Note that each step does hold on its own time scale, and the conclusion is plausibly provable: the averaged weight subsystem is uniformly exponentially stable under (H1), and the geometry Lipschitz constant is O(δ).

**Fix.** 1. Replace the opening sentence with: 'Each step is a standard theorem on its own time scale; their composition on the O(k_e/(mδ)) horizon is formal.'
2. Step 1: use the exact identity. For w := e − r(x,ϑ)/k_e, ẇ = −k_e w − ṙ/k_e, so ‖w(t)‖ ≤ e^{−k_e t}‖w(0)‖ + sup‖ṙ‖/k_e² for all t ≥ 0.
3. Steps 2–3: either give the long-horizon argument, or cite a mixed-time-scale averaging theorem (Sastry–Bodson §4.4, or averaging under attraction in Sanders–Verhulst–Murdock Ch. 5) and verify its hypotheses. The argument would combine a near-identity averaging transform of the joint system, uniform exponential stability of the averaged weights on 𝒢, Tikhonov in slow time with an O(ε_avg) perturbation, and Grönwall in (c,s).
4. Update 'Strongest objection' in LITERATURE_TRANSPORT.md to match.

### G.4 [minor] §6b (H2a) 'k_e ≫ ω_ref ≫ every parameter rate', and Derivation step 1 'e = r/k_e + O(ω_ref/k_e²)'

CONFIRMED (R1 and R2, merged). For e ≈ r(x_m)/k_e, k_e must dominate more than ω_ref:
- the bandwidth of t ↦ r(x_m(t))∇_ϑΦ(x_m(t)), which is about sup‖ẋ_m‖/s_min (the kernel-crossing rate);
- Lip_x(f − f̂);
- and ‖r‖/k_e must be ≪ s_min.
With k_e/ω_ref → ∞ and everything else fixed, (a) is still true, so the hypothesis is mis-scaled rather than the limit being false. The document's own configurations (A/s ≈ 3–8, small Lipschitz constant) are affected by only 1–4%.

**Argument.** Judge's frozen-parameter runs (h2a_check.py): exact periodic e, project grad_g, k_e/ω_ref = 20.

(i) m = 8, s = 0.3, x_m = A sin t:
- A/s = 2.9: ratio |avg|/|pred| = 1.009, cosine 1.000
- A/s = 8: 0.987, cosine 1.000
- A/s = 33: 0.859, cosine 0.999
- A/s = 100: 0.925, cosine 0.981
R1 reports up to 95% error in the weight Jacobian at A/s = 1000, and the slowest rate overestimated by 29×.

(ii) The document's reference, m = 2, f = Lx + sin 2x: ratio 0.97, 0.37, 0.18, 0.50, 0.16 for L = 0, 10, 15, −15, −60, with cosine as low as 0.42.

**Fix.** Replace (H2a) with (H2a′): k_e ≫ ω_eff := max{ω_ref, sup_t‖ẋ_m(t)‖/s_min, sup_X‖∂_x(f − f̂)‖} ≫ every rate of the joint linearisation, and sup‖r‖/k_e ≪ s_min. Write step 1 as e = r(x_m)/k_e + O(ω_eff·sup‖r‖/k_e²).

### G.5 [minor] §6b Proposition G(a): '…the noiseless and leak-free case of Theorem H(i)'

CONFIRMED (R1 minor and R2 cosmetic, merged). This answers the targeted question: (a) is not literally a special case of H(i).
- H(i) uses saturated features Φ̃ in unconstrained coordinates (α,γ,κ), with G acting on those coordinates.
- At σ = τ = 0, written in (a,c,s), the H flow is a gradient flow of L̄/k_e (built with Φ̃) under a different, state-dependent metric.
- The formula −(1/k_e)G∇(δL̄/δμ) is the same computation with Φ in place of Φ̃, and the factors m/k_e and 1/k_e are correct.
- Side note: test_averaged_drift_is_the_mean_field_langevin_drift substitutes e := r/k_e directly. It checks the first-variation identity, not the Tikhonov or averaging approximation.

**Argument.** By the chain rule, with a = a_max tanh(α/a_max) and s = s_min + (s_max−s_min)sig(κ):
- ȧ = −(g_a/k_e)sech⁴(α/a_max)∂_aV_L
- ṡ = −(g_s/k_e)[(s_max−s_min)sig′(κ)]²∂_sV_L
Neither extra factor is identically 1.

**Fix.** Reword: '…whose measure form −(1/k_e)G∇(δL̄/δμ) is the drift of (𝔄) at σ = τ = 0 with Φ in place of Φ̃ and (a,c,s) in place of (α,γ,κ), by the same averaging computation as in H(i); valid while the projection is inactive.' Also add a direct closed-loop check of (a) to the experiment. On the stored violating run it gives cosine ≥ 0.9985 and ratio 0.93–1.00.

### G.6 [minor] §6b Setting/Hypotheses; Proposition G(a) and G(b)

CONFIRMED (R1 and R2, merged). Part (b) also needs:
- the closed-loop and reduced trajectories to stay in a compact 𝒢′ ⋐ 𝒢, where (H1) holds (this is Thm 11.1's hypothesis (ii));
- the parameters to stay in int Θ_K^ε;
- a*(c,s) to lie strictly inside the weight box.
If a* reaches the box, the projection is active at the quasi-steady state: ∇_aL̄ ≠ 0 there, and step 4's envelope identity fails. Variable projection pushes Gaussians toward coalescence with growing |a*|, so these conditions are not automatic.

**Argument.** Reduced flow from the violating initial geometry (judge run, reduced_violating.py):
- max|a*| rises 7.8 → 9.8 (800 s) → 12.7 (4000 s);
- λ_min(M) falls 1.9e-3 → 7.0e-4.
R1 reports about 17 and 4.4e-4 later in the run, against the bound a_max − ε_K = 19.95. In the stored runs the projection is inactive (|a| ≤ 3.4).

**Fix.** Add (H0): on [0, τ₁/δ], (c,s)(t) and (c̄,s̄)(δt) lie in a compact 𝒢′ ⋐ 𝒢; ϑ_i(t) ∈ int Θ_K^ε; and sup_{𝒢′} max_i‖a*_i‖ < a_max − ε_K. State (a) under (H0) and (H2a′), and (b) under (H0), (H1), (H2a′), (H2b) and (H2c).

### G.7 [minor] §6b title and Status paragraph ('the time scales come from k_e and the reference period'); LITERATURE_TRANSPORT.md §G 'Distinction'

CONFIRMED (R1 and R2, merged). The weight rate and the geometry rate both carry the factor m/k_e. Their ratio R = g_a m₀/(δ‖Γ^{1/2}∇²L*Γ^{1/2}‖), and the slope S, therefore contain neither k_e nor ω_ref. The weight/geometry separation is the designer's gain ratio, exactly as in [S1, S2]. The loop supplies only the error-channel and averaging scales and the data distribution ρ. The title also asserts variable-projection descent unconditionally, but in the violating regime the geometry follows −∇L̄ at the current weights instead.

**Argument.** This follows from the algebra above. Changing k_e in the satisfying configuration rescales both rates by the same factor and leaves R = 77.5 unchanged.

**Fix.** Status: 'The error-channel and averaging scales, and the data distribution (the occupation measure of the reference), come from the closed loop. The weight/geometry separation is a gain choice, g_c, g_s ≪ g_a·min(m₀/‖Γ∇²L*‖, 1/‖Da*Γ^{1/2}‖²), as in [S1, S2].' Make the same change in LITERATURE 'Distinction'. Retitle: 'Proposition G — the slow geometry performs joint descent on L̄, and variable-projection descent when the geometry gains are small'.

### G.8 [minor] Proposition G(b): '(c,s)(t) = (c̄,s̄)(δt) + o(1)' on time intervals of length O(1/δ)

CONFIRMED (R1; also absorbs R1's δ-independent offset point and R2's O(1/δ) point).
- The o(1) does not vanish as δ → 0 alone. The finite-k_e error channel adds a δ-independent relative velocity error and shifts the weight equilibrium.
- The limit in which o(1) vanishes is not named.
- The initial condition (c̄,s̄)(0) = (c,s)(0) is not stated.
- The natural slow time is k_e/(mδ), not 1/δ.

**Argument.** With parameters frozen at the satisfying initial geometry and a = a*, the geometry velocity ratio is 0.956 at k_e = 20, 0.996 at 80 and 0.9996 at 320, whatever δ is. The weight velocity at a* is not zero with the exact e (4.2e-4), so the weight equilibrium is shifted. R1 reports ratios 0.970, 0.969 and 0.968 for δ = 0.004, 0.001 and 0.0005.

**Fix.** State: 'With (c̄,s̄)(0) = (c,s)(0), sup_{t ≤ τ₁k_e/(mδ)}‖(c,s)(t) − (c̄,s̄)(δt)‖ → 0 as (ω_eff/k_e, ε_avg, S, 1/R) → 0 with τ₁ fixed; formally O(ω_eff/k_e + ε_avg + S + 1/R).'

### G.9 [minor] (H2b) geometry rate '(mδ/k_e)‖∇²L*‖'; Numerics '−(m/k_e)Γ∇L*'; table header in exp_variable_projection.py

CONFIRMED (R1 and R2, merged).
- With gains (g_a, δγ_c, δγ_s), the reduced linearisation is −(mδ/k_e)Γ∇²L* with Γ = diag(γ_c I_n, γ_s). As written, (H2b) does not change when γ_c → ∞.
- Γ is never defined in §6b.
- The Numerics/table prediction '−(m/k_e)Γ∇L*' omits δ, unless Γ means diag(g_c, g_s) as in the script's convention. That reading conflicts with (b).
The experiments use γ = 1, so no reported number changes.

**Argument.** Linearise (b): d(c,s)/dt = −(mδ/k_e)Γ∇L*.

**Fix.** Define Γ := diag(γ_c I_n, γ_s) in the Setting. Restate (H2b) as: g_a m₀ m/k_e ≫ (mδ/k_e)‖Γ^{1/2}∇²L*Γ^{1/2}‖. In the Numerics and the table header, write '−(mδ/k_e)Γ∇L*'.

### G.10 [cosmetic] Proposition G(a), last sentence: 'In particular L̄ evaluated at the current weights is non-increasing'; Numerics 'decreases at every sample'

CONFIRMED (R1). This holds along the averaged flow, and for period-averaged parameters. It does not hold instantaneously along the closed loop.

**Argument.** On the stored violating run (check_a_violating.py), L̄ evaluated at the instantaneous parameters (20 samples per period) increases at 13.7% of steps, by at most +0.17%. Evaluated at period-mean parameters, it increases at 0 of 126 steps.

**Fix.** Write: 'L̄ is non-increasing along the averaged flow; along the closed loop its value at period-averaged parameters decreases, up to O(ε_avg) oscillation within a period.' In the Numerics, say 'non-increasing, sampled once per period' rather than 'decreases at every sample'.

### G.11 [cosmetic] §6b Numerics, third bullet: 'consistent with an averaging error of order (fastest weight)/ω_ref ≈ 0.14'

CONFIRMED (R2), with a refined attribution. The 5–6% shortfall does not come from averaging over parameter motion. About 2.5% is the bias of the script's 20-period forward difference. The rest is the O(1/k_e) error-channel correction, plus about 0.6% from the slope (S = 0.006).

**Argument.** - With parameters frozen, so no averaging is involved at all, the ratio is 0.956 at k_e = 20, 0.996 at 80 and 0.9996 at 320 (cosine ≥ 0.9999).
- On the stored satisfying run, comparing against the prediction at mid-window instead of at window start raises the median ratio from 0.944 to 0.969 (window_bias.py).

**Fix.** Write: 'Its magnitude is 0.94–0.95 of the prediction. About 2.5% is the bias of the 20-period forward difference (0.969 against the mid-window prediction); the rest is the O(1/k_e) error-channel correction (frozen-parameter ratio 0.956 at k_e = 20, 0.996 at k_e = 80).' Optionally, switch the script to a centred difference.

### G.12 [cosmetic] §6b Numerics, bullets 1–2; experiments/exp_variable_projection.py table code; table_violating.md and table_satisfying.md

CONFIRMED (R1 and R2, merged).
(i) 'Two scales hold' rests only on a necessary consequence (monotone L̄). The direct test passes and should be reported.
(ii) Time-scale ratios are reported only at the initial geometry.
(iii) '0.1×' is actually 0.14×.
(iv) The rows L*, L̄, ‖∇L*‖ and tracking RMS sit under 'first quarter | last quarter' but print the first and last samples.
(v) The line is labelled '(H2)' instead of (H2a)/(H2b).
(vi) The 'fastest weight mode' is not the fastest parameter rate that (H2a) refers to.

**Argument.** - Direct check of (a) on the violating run (119 periods): cosine ≥ 0.9985; magnitude ratio 0.93–1.00 (weights 0.945–1.007, geometry 0.922–0.984).
- Joint fastest rate (m/k_e)λ_max(G^{1/2}∇²L̄G^{1/2}) at the violating initial geometry, with a = a*: 0.068, against the reported 0.0163. R2 reports 0.019–0.033 along the run at the actual weights. So ω_ref/(fastest rate) is about 15–30, not 61.4, which still satisfies (H2a).
- In the satisfying configuration the joint rate is 0.137, about the reported value.

**Fix.** - Add the joint-velocity check of (a) and an S column to the script.
- Report the rate ratios along the run (minimum and maximum).
- Relabel the rows 'first sample | last sample'.
- Use the labels (H2a)/(H2b).
- Write 0.14×.
- Report the joint fastest rate.

Refuter verdicts: holds_with_fixes; holds_with_fixes

## C — verdict: holds_with_fixes

Verdict: holds_with_fixes. One major issue (the solution concept, which makes C(i) false for Filippov solutions), four minor and five cosmetic. No reviewer issue is rejected.

What holds, re-derived by hand:
- C(i) holds for selection solutions, and Lemma C.3's affine-in-λ argument holds for any G ≻ 0.
- C(ii): both C ≤ b(t)²+4ηR and the gap clause hold.
- Both displays of (iii) hold.
- (iv) holds: the three Young splits, the gain condition, α, β, forward invariance and comparison.
- The Scope arithmetic is right, and (v) is right for a fixed π.
- In Proposition C.5, D⁺V_OT = σ_m(δ−δ²), exact for δ < ½.
- Remark C.6 citations, checked against the source texts: AGS Thm 7.3.2 is the PC-space semiconcavity inequality, Ex. 9.1.5 is the two-atom ℝ² counterexample, Lemma 9.2.1 is 1-convexity along generalised geodesics. Santambrogio [71] §4.4 has W² = a² + min{(1−t)², (1+t)²}. Carlen–Craig [73] cite Ex. 9.1.5.

Checks run (project code, project venv, no project file edited; the three Theorem C tests pass). Scripts are in scratchpad\judgeC4: filippov_check.py, scope_check.py, c5_v_check.py.
- Tie equilibrium: ½v_id + ½v_swap = 0 exactly, for G = I and for an anisotropic diagonal G.
- Sliding with a common push: W(T) = 2.236 > max{W0, supΓ/σ} = 2.
- The selection solution from the same tie decays exactly: 0.006721 = (1−dt)^5000.
- 4000 random configurations: C ≤ b²+4ηR held, √C ≤ b+2η held (worst slack −1.0e−2), and 0 violations of Δ_gap ≤ 2R²/m in the finite case. A single-atom prior gives Δ_gap = ∞.
- Closed-loop labelled reshuffle: max|e| = 2.5e−13, V_OT rises from 0.005 to 0.125, path length 0.90 against 0.10 for the OT leak.
- Moving ideal in (v): b = 0.5, best fixed-π bias 1.0.

Decisions on Reviewer 1 (R1):
1. Filippov sliding breaks C(i): CONFIRMED, downgraded from fatal to major. (i) holds verbatim for selection solutions, and (ii)–(iv) hold for every Filippov solution.
2. Existence sentence mixes up Carathéodory and Filippov: CONFIRMED, minor, merged into issue 1. ẋ = −sgn x has no Carathéodory solution from 0.
3. (v) bias claim, the 'global strict' wording and the C.6 trade-off sentence: CONFIRMED, minor (reproduced).
4. Proposition C.5 and Theorem F: CONFIRMED, minor (reshuffle with e ≡ 0). The δ range and the 'every kernel' wording are CONFIRMED as cosmetic.
5. Scope paragraph:
   - Δ_gap = 0 whenever the ideal–prior matching switches: CONFIRMED, minor.
   - 'Generically': CONFIRMED as imprecise. The k-cycle worry is moot by R2's antisymmetry argument.
   - The d_min²/m refinement: correct but optional.
6. b(t) undefined, the last line of the (ii) proof, Δ overloading: CONFIRMED, cosmetic.
7. No 'Under (S)' line: CONFIRMED, cosmetic. C at ties is merged into issue 1.
8. Remark C.6 wording: CONFIRMED, cosmetic. The citations are correct.
9. Lemma C.3 strictness when A_ρ = 0, and the contraction identity being per factor: CONFIRMED, cosmetic.

Decisions on Reviewer 2 (R2):
1. Solution concept: CONFIRMED, major (duplicate of R1.1).
2. O(1/m) is an artefact of the R-dependent bound: CONFIRMED, minor. √C ≤ b(t)+2η is valid and never violated. The stated results remain true; whether m^{−1/2} is the 'true' scaling depends on the geometry (UNCERTAIN, not needed).
3. 'Generically' is unnecessary, and the Δ_gap = ∞ case is missed: CONFIRMED, minor. The antisymmetry argument was re-derived.
4. (v): CONFIRMED (duplicate). The Δ_gap > 0 ⇒ constant optimal set argument was re-derived.
5. Proposition C.5 and Theorem F, and 'steady states coincide': CONFIRMED (duplicate). The latter is empirical only.
6. δ ∈ (0,½), a_max−ε_K ≥ 1, ambiguous b_lab: CONFIRMED, cosmetic.
7. Lemma C.3 strictness and d_i: CONFIRMED (duplicate).
8. Notation b(t) and Δ_gap(t): CONFIRMED (duplicate).
9. Trade-off sentence: CONFIRMED, merged into the (v) issue.

UNCERTAIN, and not to be asserted in the text without proof: existence of selection solutions in the general closed loop. Two-way ties are never attracting because Proj is monotone, but multi-way ties, the degenerate case h = 1 and accumulating switch times are still open.

Confirmed issues by severity (exact fixes are in the issues list):
- MAJOR: the solution concept and C(i).
- MINOR:
  1. Scope paragraph: Δ_gap = ∞ case, 'generically', Δ_gap = 0 at switches.
  2. The R-dependent bound in (ii), and the O(1/m) claim that follows from it.
  3. The bias claim in (v), and the trade-off sentence in Remark C.6.
  4. 'As Theorem F prices it' in Proposition C.5.
- COSMETIC:
  1. Proposition C.5 example details.
  2. Notation.
  3. Hypotheses line.
  4. Lemma C.3 strictness and the per-factor contraction identity.
  5. Remark C.6 wording.

### C.1 [major] §4 'Solution concept' paragraph ('Carathéodory/Filippov ... any measurable selection ... so solutions exist ... every statement holds for every solution'), with Theorem C(i), the C of (ii)–(iii) at ties, and Lemma C.3's last sentence ('no chattering')

The paragraph admits Filippov solutions and claims that every statement holds for every solution. Its existence argument (a measurable selection gives a bounded field) proves existence of Filippov solutions only, not Carathéodory ones.

Tie surfaces of (L_σ) are crossing or repelling, never attracting, and Filippov solutions can slide on repelling ties. Along a sliding motion C(i) is false: the differential inequality, the ultimate bound and the exact pure-leak decay all fail. Linearity of (P) in y does not rescue (i), because its comparison point ϑ⁰_{τ(i)} differs between the tied fields.

Parts (ii)–(iv) do hold for every Filippov solution if C is read as the maximum over optimal τ. Lemma C.3 is correct for the fixed-τ ODE it states, but 'no chattering' does not exclude Filippov sliding. At ties, C in (iii) is ambiguous: with an arbitrary optimal pair it can be too small along a sliding motion.

I rule major rather than fatal. (i) holds verbatim for selection solutions: the class the code integrates, and the class Lemma C.3 constructs for the pure leak. The section's main results (ii)–(iv) survive for all solutions.

**Argument.** Setting: n = 1, m = 2, ideal = prior = {(a,c,s) = (0,0,1), (1,0,1)}, both particles at (½,0,1), e ≡ 0, projection inactive (a_max = 2, ε_K = 0.3).

1. Stationary Filippov solution. c_id = c_swap = ¼, and with project code ½v_id + ½v_swap = 0 exactly (also for G = diag(2.3, 0.7, 1.9)). So the constant state is a Filippov solution with W(μ_t,μ⁰) ≡ ½ and D⁺Z = 0, while C(i) and Lemma C.3 assert D⁺Z ≤ −2σ_mZ < 0 and W = e^{−σ_m t}/2.

2. Sliding with a push. Take coincident particles, prior a = ±1, G = I and a common push γ in a. This is realisable in closed loop with x(t) at the centre c and e constant, so Γ = |γ|. The sliding solution ẋ = γ − σ_m x from x = 0 gives W(T) = 2.236 for γ = 2, σ = 1, against the claimed bound max{W0, supΓ/σ} = 2.

3. Control. The selection solution produced by the code's tie-break from the same point decays exactly: W(5)/W(0) = 0.006721 = (1−dt)^5000.

Why two-way ties never attract: the normal rate of c_ρ − c_τ under the τ-field minus that under the ρ-field equals (2/(mσ_m))Σ⟨u^τ_i − u^ρ_i, G⁻¹(P_i(u^τ_i) − P_i(u^ρ_i))⟩ ≥ 0, by monotonicity of the per-factor Pomet–Praly and box maps.

Why (iii) survives for Filippov solutions: polarisation gives, for weights θ_k on tied optimal τ_k, D⁺Z ≤ −2σ_mZ + (σ_m/2)Σ_kθ_kA_{ττ_k} + (2Z)^{1/2}Γ. This is tight at the equilibrium (0 = −σ/4 + σ/4). (P) is termwise, and (ii) holds for every optimal τ.

**Fix.** 1. Replace the paragraph with:
'Solution concept. The right-hand side of (L_σ) is discontinuous at assignment ties. A solution is an absolutely continuous (x,ϑ) such that, for a.e. t, ϑ̇_i = Proj_K(ϑ_i, G∇g_t(ϑ_i) − σ_m(ϑ_i − ϑ⁰_{τ(t)(i)})) for some τ(t) optimal at ϑ(t) (a selection solution; this is what transport.py integrates). Filippov solutions exist by Filippov's theorem (measurable, locally bounded right-hand side), but they include sliding motions on repelling ties, along which (i) fails. Example: ideal = prior = {a=0, a=1}, both particles at a=½, e ≡ 0 is a stationary Filippov solution. Parts (ii)–(iv) hold for every Filippov solution, with C taken as the maximum over optimal τ. Part (i) and the exact decay hold for selection solutions. For a Filippov solution with weights θ_k on tied optimal τ_k, (i) weakens to D⁺Z ≤ −2σ_mZ + (σ_m/2)Σ_kθ_kA_{ττ_k} + (2Z)^{1/2}Γ, with A_{ττ'} := (1/m)Σ‖ϑ⁰_{τ(i)} − ϑ⁰_{τ'(i)}‖²_{G⁻¹}. V_OT and Z are absolutely continuous along solutions, and every inequality holds for a.e. t.'

2. Existence of selection solutions:
- State it for the pure leak (Lemma C.3 constructs one).
- State that two-way ties are never attracting (monotonicity of Proj).
- Assume existence in general, unless multi-way ties, h = 1 and accumulating switch times are handled.

3. In (iii), write 'C as in (ii), with τ the assignment used by the law at t (for Filippov solutions, the maximum over optimal τ)'.

4. In Lemma C.3, replace 'so there is no chattering' with 'so no selection solution of the pure leak switches; Filippov sliding on the tie is excluded by the solution concept'.

5. Delete 'every statement holds for every solution'. In the §0 table, write 'exact dissipation (selection solutions)'.

### C.2 [minor] §4 'Scope' paragraph ('Δ_gap ≤ 2R²/m generically ... Hence η̄ < R/(2m) ... small-m, near-prior local result') and the Notation definition of Δ_gap

The paragraph has three problems.
(a) The conclusion 'Hence η̄ < R/(2m) ... O(1/m) ... small-m' is false when Δ_gap = +∞. This happens for a single-atom prior μ⁰ = δ_{ϑ⁰}, which is the classical σ-modification. There the gap condition is void and (iv) is m-uniform, with η̄ fixed by 8K²λ_G²η̄² ≤ k_eσ_m alone.
(b) 'Generically' is unnecessary: a finite Δ_gap is always ≤ 2R²/m.
(c) The paragraph omits that Δ_gap = inf_t Δ_gap(t) is 0 whenever the optimal ideal–prior assignment switches, or nearly ties, along a moving ideal. Then the second clause of (ii), the second display of (iii) and all of (iv) are vacuous.

**Argument.** (b) Let ρ₀ be optimal and suppose every transposition of ρ₀ is cost-neutral. With B(i,j) = ⟨ϑ*_i, ϑ⁰_{ρ₀(j)}⟩_{G⁻¹}, the matrix S(i,j) := B(i,j) − (B(i,i)+B(j,j))/2 is antisymmetric. Then excess(σ) = −(2/m)Σ_iS(i,σ(i)) = −excess(σ⁻¹), and optimality forces every excess to be 0, so all assignments tie. Hence a finite Δ_gap(t) means some transposition has positive excess, which is ≤ 2R²/m. Numerically: 4000 finite-gap configurations, 0 violations. Single-atom prior: the spread of all assignment costs is exactly 0.
(c) If two assignments with different cost functions are both optimal at t₀, the non-optimal one exceeds b(t)² by an amount that tends to 0 as t → t₀, so Δ_gap = 0. Example: reviewer 1's rotating ideal, which ties at ωt = π/2.

**Fix.** Replace the first two sentences of Scope with:
'If Δ_gap < ∞ (not all ideal–prior assignments tie), then Δ_gap ≤ 2R²/m always: if every transposition of an optimal ρ₀ were cost-neutral, all assignments would tie. Hence η̄ < R/(2m) and non-vacuity forces b < Δ_gap/(4√2R) = O(1/m). If Δ_gap = +∞ (e.g. a single-atom prior, the classical σ-modification), the gap condition is void and (iv) is m-uniform. Δ_gap > 0 forces the set of optimal ideal–prior assignments to be constant in t; if it switches along a moving ideal, Δ_gap = 0 and the second clause of (ii), the second display of (iii) and (iv) are vacuous.'

Optionally add the sharper Δ_gap ≤ (2/m)·min over cost-changing pairs of ⟨ϑ*_j − ϑ*_k, ϑ⁰_{ρ₀(j)} − ϑ⁰_{ρ₀(k)}⟩_{G⁻¹}.

### C.3 [minor] Theorem C(ii) bound C ≤ b(t)² + 4ηR (and its proof), the gap conditions 4ηR < Δ_gap in (iii) and η̄ < Δ_gap/(4R) in (iv), and the Scope paragraph's O(1/m) conclusion

The statements are correct, but the O(1/m) basin and bias bound come from the R-dependent perturbation bound |c_ρ(ϑ) − c_ρ(ϑ*)| ≤ 2ηR, not from the problem itself. Yet the Scope paragraph presents them as intrinsic ('in contrast with Theorem B').

√c_ρ is a (1/m)-weighted G⁻¹ norm, so the reverse triangle inequality gives |√c_ρ(ϑ) − √c_ρ(ϑ*)| ≤ η. Hence √C ≤ b(t) + 2η, and C = b(t)² whenever 4η(b(t)+η) < Δ_gap. This removes R and gives η̄ < √Δ_gap/2 ≤ R/√(2m), with b < η̄/√2, i.e. O(m^{−1/2}).

**Argument.** Proof: √C = √c_τ(ϑ*) ≤ √c_τ(ϑ) + η ≤ √c_{ρ₀}(ϑ) + η ≤ √c_{ρ₀}(ϑ*) + 2η = b(t) + 2η. Then C − b(t)² ≤ 4η(b(t)+η) < Δ_gap ≤ Δ_gap(t) forces C = b(t)².

Numerically, on 4000 random configurations (m = 2..5, anisotropic G, all optimal π* and τ), max[√C − (b + 2η)] = −1.0e−2, never violated. The document's bound also holds there (max slack −0.12).

**Fix.** 1. In (ii) add: 'Moreover √C ≤ b(t) + 2η, so C = b(t)² whenever 4η(b(t)+η) < Δ_gap.'
2. In (iii), replace 'Under 4ηR < Δ_gap' with 'Under 4η(b+η) < Δ_gap'.
3. In (iv), replace 'η̄ < Δ_gap/(4R)' with '4η̄(b+η̄) < Δ_gap'.
4. In Scope, replace 'Hence η̄ < R/(2m) ... = O(1/m)' with 'Hence η̄ < √Δ_gap/2 ≤ R/√(2m) and b < η̄/√2, i.e. O(m^{−1/2}) in RMS, O(1/m) in V_OT'.
If this is not adopted, at least drop the implication that the O(1/m) is intrinsic.

### C.4 [minor] Theorem C(v) ('satisfies the global strict inequality of (iii) ... Choosing the minimising π gives bias b'), and the last sentence of Remark C.6 ('trade-off between a global guarantee with a labelled bias and a local guarantee with the optimal bias')

π is fixed in time, so the bias that enters (v) is sup_t b_lab(π)(t). Then min_π sup_t b_lab(π)(t) ≥ sup_t min_π b_lab(π)(t) = b, strictly whenever the optimal ideal–prior assignment changes over time. So 'the minimising π gives bias b' is false for such moving ideals. It is true when one assignment is optimal at every t: a static ideal, or whenever Δ_gap > 0 (costs minus b(t)² take values in {0} ∪ [Δ_gap, ∞) and are continuous in t, so the optimal set is constant).

In that case the labelled leak achieves the optimal bias b globally, for V_{π_b}. This contradicts Remark C.6's 'labelled bias vs optimal bias' framing. The real contrast is a global inequality for a label-dependent V_π, which is large at relabelled starts, against a local, gap-protected inequality for the label-free V_OT.

'The global strict inequality of (iii)' is also ambiguous, because (iii)'s second display is not global for the OT leak.

**Argument.** Moving ideal a*_1(t) = t, a*_2(t) = 1−t on [0,1], prior a ∈ {0,1}, G = I. Computed: b = sup_t min_π b_lab = 0.500, while min_π sup_t b_lab(π) = 1.000. Reviewer 1's rotating ideal gives 2r against b = √2·r. In both examples inf_t Δ_gap(t) = 0.

**Fix.** 1. Define b_lab(π)(t) := ((1/m)Σ‖ϑ*_{π(i)}(t) − ϑ⁰_i‖²_{G⁻¹})^{1/2} and b_lab(π) := sup_t b_lab(π)(t).
2. Replace 'satisfies the global strict inequality of (iii)' with 'satisfies the second display of (iii) at every state (no gap condition)'.
3. Replace 'Choosing the minimising π gives bias b.' with 'Choosing π to minimise b_lab(π) gives bias min_π sup_t b_lab(π)(t) ≥ b, with equality when one assignment is optimal for (μ*_t, μ⁰) at every t (e.g. a static ideal, or whenever Δ_gap > 0).'
4. End Remark C.6 with: 'So the trade-off between a global inequality for a labelled V_π (large at relabelled starts, hence the reshuffling of Proposition C.5) and a local, gap-protected inequality for the label-free V_OT reflects this geometry.'

### C.5 [minor] Proposition C.5, last sentences ('The steady states coincide, but the transient costs tracking error, as Theorem F prices it (results/ot_sigma_mod/transient.md: particle path 5.8× longer)')

'As Theorem F prices it' is a non sequitur. Theorem F bounds transport by integrated tracking error only for (L). Its (L_σ) clause adds σ_mW(μ_r,μ⁰). The labelled leak is neither: its motion is paid for by the labelled distance to the prior, not by e. So a long particle path implies no tracking-error cost. In the proposition's own example the reshuffle happens with zero tracking error.

'The steady states coincide' is an experimental observation, not a proved fact.

**Argument.** Closed loop with project code: plant f = f̂(μ*), ideal = prior = {(0,0.3,1), (1,0.3,1)}, x(0) = x_m(0), start (0.9, 0.1), σ_m = 0.8, RK4 for 8 s.
- Labelled leak: max|e| = 2.5e−13 (a_1 + a_2, hence f̂, is invariant), V_OT rises from 0.005 to 0.125 at t = 0.73, path length 0.899.
- OT leak: path length 0.100.
The cited experiment shows only RMS‖e‖ 0.0771 vs 0.0743 on [0,5] (+3.6%), for a path 5.8× longer.

**Fix.** Replace the sentence with: 'In the reshuffling experiment (results/ot_sigma_mod/transient.md) the steady states coincide, the labelled leak's particle path is 5.8× longer and the early RMS‖e‖ is about 4% larger (0.0771 vs 0.0743 on [0,5]). Theorem F does not price leak-driven transport: in the example above the reshuffle occurs with e ≡ 0.'

### C.6 [cosmetic] Proposition C.5 statement and Example

The example has four loose points:
- It needs 0 < δ < ½, because the swap is optimal only then (δ = ½ is a tie, and δ > ½ gives D⁺V_OT < 0).
- It needs a_max − ε_K ≥ 1, with the common (c, s) in the inner set, so that ideal = prior lies in Θ_K^ε.
- 'b_lab = 0' is ambiguous now that b_lab depends on π (b_lab(id) = 0, b_lab(swap) = 1).
- 'transports every kernel across the domain' overclaims, since fixed points of the relabelling do not move.

**Argument.** With σ_m = 0.8, project code gives D⁺V_OT = +0.038, +0.128, +0.198 at δ = 0.05, 0.2, 0.45, matching σ_m(δ−δ²) exactly. At δ = 0.55 and 0.8 it gives −0.162 and −0.032, where the formula would give +0.198 and +0.128.

**Fix.** Write: 'Take n = 1, G = I, a_max − ε_K ≥ 1, ideal = prior = {a = 0, a = 1} (other coordinates equal and in the inner set), and network weights (1−δ, δ) with 0 < δ < ½. Then b = b_lab(id) = 0 ...'. Replace 'transports every kernel across the domain' with 'moves every kernel in a non-trivial cycle of the relabelling to its prior position'.

### C.7 [cosmetic] Notation (b, b(t), Δ_gap); proof of (ii), last line; Lemma C.3 (Δ_i, d_i, μ_0 vs μ⁰)

There are four notation gaps:
- b(t) is used in Δ_gap, in (ii) and in the proof, but only b := sup_t W(μ*_t,μ⁰) is defined.
- The '+∞ if the set is empty' convention is per time but is attached to the infimum.
- The last line of the proof of (ii) ('No assignment cost lies strictly between b² and b² + Δ_gap') is false with b = sup_t when b(t) < b. It must read b(t)².
- Δ denotes three things (the diameter, Δ_gap and Δ_i), d_i is undefined, and μ_0 (the initial network) and μ⁰ (the prior) differ only in the index position.

**Argument.** Take b(t)² = 0 and Δ_gap = 1 at time t, and b² = 2 at another time. A cost of 2.5 at t is allowed, yet it lies in (b², b²+Δ_gap). Only (b(t)², b(t)²+Δ_gap) is guaranteed empty.

**Fix.** 1. Define b(t) := W(μ*_t, μ⁰) and b := sup_t b(t).
2. Define Δ_gap(t) := min{c_ρ(ϑ*(t)) − b(t)² : c_ρ(ϑ*(t)) > b(t)²} with min ∅ := +∞, and Δ_gap := inf_t Δ_gap(t).
3. In the proof of (ii), write '... ≤ b(t)² + 4ηR. No assignment cost at time t lies strictly between b(t)² and b(t)² + Δ_gap.'
4. In Lemma C.3, rename Δ_i to q_i and define d_i := ϑ_i(0) − ϑ⁰_{τ(i)}. Consider writing μ_init for the initial network.

### C.8 [cosmetic] Statement of Theorem C (no hypotheses line)

Theorem B opens with 'Under (S)'; Theorem C has no hypotheses line, yet it uses three standing hypotheses:
- (S4): (P) with ϑ' = ϑ*_{π*(i)}, and V_* < ∞.
- (S5): ε_max in (iii) and (iv).
- (S3): forward invariance of Θ_K, which bounds every distance by R.

**Argument.** Without (S4), step (P) with ϑ' = ϑ* in the proof of (iii) is unavailable. Without (S5), ‖e‖ε_max is undefined along the trajectory.

**Fix.** Open with 'Under (S), along every solution of (L_σ) in the sense above:'. Optionally note in (iv) that ‖e‖ < η̄ on the invariant set, so (S5) holds with X := {x : ‖x − x_m(t)‖ ≤ η̄}.

### C.9 [cosmetic] Lemma C.3, last sentence ('becomes strictly optimal immediately'); proof of (i) (contraction identity ‖Proj(ϑ,Gy)‖²_{G⁻¹} = ‖y‖²_G − h(2−h)(ϑᵀGy)²/(ϑᵀGϑ))

Strict optimality fails against any ρ with A_ρ := (1/m)Σ‖ϑ⁰_{τ(i)} − ϑ⁰_{ρ(i)}‖² = 0, which happens with repeated prior atoms. This is harmless, because such assignments induce the same leak.

The displayed identity is a per-factor identity for a norm-ball factor, valid only where that factor's projection is active. The s-box has the analogous form (1−h)². Only the inequality ‖Proj(ϑ,Gy)‖_{G⁻¹} ≤ ‖y‖_G is needed.

**Argument.** With ϑ⁰_1 = ϑ⁰_2, and τ, ρ differing by that transposition, c_ρ ≡ c_τ for all λ.

**Fix.** 1. In Lemma C.3, write: '... becomes strictly optimal immediately against every ρ with A_ρ > 0; assignments with A_ρ = 0 induce the same velocity.'
2. In the proof of (i), write: 'the projection is a G⁻¹-contraction, factor by factor: for an active norm-ball factor ‖Proj‖² = ‖y‖²_G − h(2−h)g(aᵀy)²/|a|², and for the s-box ‖Proj‖² = (1−h)²‖y‖²_G.'

### C.10 [cosmetic] Remark C.6, first and third bullets

The citations are correct. I checked them against the texts:
- AGS Thm 7.3.2 is the PC-space inequality.
- AGS Ex. 9.1.5 is the two-atom counterexample in ℝ² with μ¹ = ½(δ_(0,0) + δ_(0,−2)).
- AGS Lemma 9.2.1 is 1-convexity along generalised geodesics.
- Santambrogio [71] §4.4 has W² = a² + min{(1−t)², (1+t)²}.
- Carlen–Craig [73] cite Ex. 9.1.5.

The wording overreaches in three places:
- Bullet 1 reads as 'for every ν'. But for ν = δ_{ϑ⁰}, ½W²(·,ν) is a potential energy and is 1-convex.
- 'Semi-concave there' suggests a dependence on dimension, but Thm 7.3.2 holds in every dimension.
- Bullet 3 reads as a universal claim about label-free convex leaks, and as such it is false.

**Argument.** V(ϑ) := dist²_{G⁻¹}(ϑ, conv supp μ⁰) is convex. So ∫V dμ is a label-free, geodesically 0-convex leak whose velocity vanishes at μ⁰, and also at every measure supported in the hull. The correct version: a label-free functional that is λ-convex (λ > 0) on the lift ℝ^{mp} has a unique minimiser, which is permutation-invariant, so all its particles coincide. It therefore cannot vanish at a μ⁰ with two distinct atoms.

**Fix.** 1. Bullet 1: write 'For some ν (e.g. two-atom ν in dimension ≥ 2), ½W²(·,ν) is not λ-convex along geodesics for any λ (AGS Ex. 9.1.5; [71] §4.4), while it is 1-semiconcave in every dimension (AGS Thm 7.3.2).'
2. Bullet 3: write 'A label-free leak that is λ-convex (λ > 0) on the lift has a unique, permutation-invariant minimiser (all particles coincide), so it cannot vanish at a μ⁰ with distinct atoms. 0-convex leaks vanishing at μ⁰ exist (e.g. toward conv supp μ⁰), but they do not single out μ⁰.'

Refuter verdicts: holds_with_fixes; holds_with_fixes

## D — verdict: holds_with_fixes

**Verdict: holds_with_fixes.** Nothing is fatal or major, and no reviewer issue was refuted except one stale sentence. I re-read §5 (lines 373-431), §1, Golse [82] and the code, and re-ran the key checks with my own scripts in scratchpad\judgeD (j1_tube_vdp.py, j2_lip.py, j3_ball_lemma.py). No project file was modified.

Every part's substance holds:
- **(a)** is true, but its cited proof does not apply.
- **(b)** is correct once its hypothesis is read as an inequality.
- **(c)** is correct. I recomputed Λ = max(L_ϑ + L_2, L_x + L_1) from Ṡ ≤ (L_1 + L_x)A + (L_2 + L_ϑ)D.
- **(d)** is true as a localized statement on the existence interval of the limit, proved by a bootstrap.
- **(e)** is correct via the integrated form.
- **The L_ϑ = O(1/ε_K) claim holds.** The measured two-sided constant is ≈ 30/ε_K, i.e. Θ(1/ε_K).

**Decisions on reviewer issues (R1 = reviewer 1, R2 = reviewer 2):**
- **(b), circular F_max (R1-1 = R2-3): CONFIRMED, minor.** F = S(F) has two roots or none for f = x³ (recomputed 0.0081 and 28.57; none when k_e = 1). For the Van der Pol run, S(F) − F ≥ 51.6 for all F.
- **(d), hypothesis missing (R1-2 = R2-1; merged with R1-5 and the second half of R2-4 on (e)'s setting): CONFIRMED, minor.** f = x² blows up before t* = 0.4304. The condition S(F) ≤ F is not open in ‖e(0)‖ (threshold 2.40). I checked the bootstrap proof.
- **(a), Golse citation (R1-5): CONFIRMED, minor.** golse.txt shows (HK1)-(HK2) and a W₁ conclusion. The round-1 Picard fix (D.6(a)) was not applied.
- **L_1 versus L₁ (R1-3 = R2-6): CONFIRMED, rated minor.** The misread L_2 = 0.81 is below the true 2.0. The intended L_2 = 33.08 is valid (sampled sup 16.1).
- **L_ϑ (R1-4 = R2-2): CONFIRMED that the claim holds.** The one-sided constant is ε_K-free: measured 88-98 for ε_K from 0.2 to 0.005. I proved the ball-ramp lemma and checked it (worst ratio −0.75 ≥ −1). This is an optional strengthening, so I downgrade it to cosmetic.
- **(e) wording (R1-7 = R2-5, plus the first half of R2-4): CONFIRMED, cosmetic.** R2's "lossy sup_m" point is downgraded from minor, because applying (e) to tail subsequences already gives limsup.
- **(c) proof rigour and unstated standing assumptions (R1-8, R2-7): CONFIRMED, cosmetic.** PARADIGM.md never states k_e > 0 or the allowed range of ε_K.
- **§0 and numerics overstatement (R1-9): CONFIRMED, cosmetic.**
- **RESULTS_PARADIGM.md mislabels 1.68 as e^{Λt} (R2-8): REJECTED.** The current file holds only the table; its substance is merged into the numerics item.
- No UNCERTAIN items.

**Confirmed issues by severity** (exact fixes are in the issue list):
- Minor: (d) localization; (b) inequality form; (a) Picard proof; renaming L_1 to M_x.
- Cosmetic: L_ϑ one-sided strengthening and L_x definition; (e) wording; (c) proof via D² and D_δ; standing assumptions; §0 row and numerics caveat (e^{ΛT} ≥ e^{496}).

### D.1 [minor] PARADIGM.md §5, Theorem D(d) lines 404-405; proof line 425 '(d) By (b) and (c)'; consequential in (e) lines 407-410 and §1 (S5)

(d) is stated under (S1), (S2) and G ≻ 0 only, with T arbitrary. Its proof needs one compact X that holds every m-solution and the limit on [0,T].

'Part (b) supplies a common X' needs (b)'s extra hypothesis, which (d) does not state, and it fails in two further ways:
- even under that hypothesis, the tube radius depends on ‖e^m(0)‖, and S(F) ≤ F is not open in ‖e(0)‖;
- (b) does not cover the paper's own Theorem D experiment.

The conclusion is true on the existence interval of the limit, by a bootstrap that needs no tube. (e) inherits this setting and must use the same X. Merges R1-2, R1-5 (setting of (e)), R2-1 and the second half of R2-4.

**Argument.** All recomputed.

**Blow-up.** Take f = x², x_m = 0, k_e = a_max = 1, e(0) = 3. Then ė ≥ e² − e − 1, so every solution blows up before t* = ∫_3^∞ de/(e² − e − 1) = 0.4304, for every m and every μ_0. On [0,1] there is no limit object.

**Not open.** Take f = x², k_e = 10, a_max = 1. S(F) ≤ F is solvable iff 100 − 4(10‖e(0)‖ + 1) ≥ 0, i.e. ‖e(0)‖ ≤ 2.40. So x_0 = 2.4 with x^m_0 = 2.4 + 1/m gets no tube from (b) for any m. Yet even the worst case ė = −10e + e² + 1 decays.

**Paper's experiment.** For the Van der Pol run (μ = 1.5, A = 1.2, k_e = 4, a_max = 6, ‖e(0)‖ = 1.772), S(F) − F ≥ 51.6 for all F in [0, 1e6].

**Bootstrap.** Let τ_m be the first exit of x^m from X_1. On [0, τ_m], (c) gives ‖x^m − x‖ ≤ e^{Λ_1T}δ_m < 1, so x^m(τ_m) is interior to X_1, a contradiction unless τ_m = T.

**Fix.** Replace lines 404-405 by:

'(d) *Mean-field limit.* Let (x, μ) be the solution from (x_0, μ_0), and let T be smaller than its existence time. Put X_1 := {y : dist(y, x([0,T])) ≤ 1}, and let Λ_1 be the Λ of (c) computed on X_1 × [0,T]. Suppose x^m_0 → x_0 and μ^m_0 → μ_0 in W, and set δ_m := ‖x^m_0 − x_0‖ + W(μ^m_0, μ_0). Then for all large m the m-particle loop exists on [0,T] and stays in X_1, and sup_{t≤T}(‖x^m_t − x_t‖ + W(μ^m_t, μ_t)) ≤ e^{Λ_1T}δ_m → 0. Under the hypothesis of (b), e.g. for bounded f, the limit is global and T is arbitrary.'

Replace line 425 by:

'(d) Let τ_m ≤ T be the first time x^m leaves X_1, with τ_m = T if it never leaves. On [0, τ_m] both solutions lie in X_1, so (c) gives ‖x^m_t − x_t‖ ≤ e^{Λ_1T}δ_m < 1 once m is large. Then x^m(τ_m) is interior to X_1, so τ_m = T. By (a), the m-solution continues while it is bounded.'

Begin (e) with 'In the setting of (d), with X := X_1, ...'.

### D.2 [minor] PARADIGM.md §5, Theorem D(b) lines 387-391; proof lines 417-418

F_max is defined through a tube whose radius contains F_max. As written this is the equation F = S(F), where

S(F) := sup{‖f(x,t)‖ : t ≥ 0, ‖x − x_m(t)‖ ≤ ‖e(0)‖ + (F + a_max)/k_e}.

That equation can have two solutions or none, so it is not a definition.

The first-exit proof needs only a sub-fixed point, S(F) ≤ F, and under that reading the conclusion is correct. The variation-of-constants bound is a convex combination of ‖e(0)‖ and (F + a_max)/k_e. Since a_max > 0 and k_e > 0, it is strictly below the radius. Merges R1-1 and R2-3.

**Argument.** Recomputed.

For f = x³, x_m = 0, ‖e(0)‖ = 0.1, a_max = 1:
- with k_e = 10, F = S(F) has two roots, 0.0081 and 28.57;
- with k_e = 1 it has no root.

For the paper's Van der Pol run, S(0) = 51.6, S(10) = 195, S(50) = 2822 and S(1e6) = 9e15. No admissible F exists.

**Fix.** Replace lines 387-391 by:

'(b) *Adaptation-free tube (m-uniform, μ*-free).* For every μ ∈ P(Θ_K), ‖∫Φ(x;·)dμ‖ ≤ a_max. For F ≥ 0 let r(F) := ‖e(0)‖ + (F + a_max)/k_e and S(F) := sup{‖f(x,t)‖ : t ≥ 0, ‖x − x_m(t)‖ ≤ r(F)}. Suppose some F ≥ 0 satisfies S(F) ≤ F. Then every solution, for every m and every μ_0, exists for all t ≥ 0, satisfies ‖x(t) − x_m(t)‖ ≤ r(F), and

    ‖e(t)‖ ≤ e^{−k_e t}‖e(0)‖ + (1 − e^{−k_e t})(F + a_max)/k_e.

Examples:
- for bounded f, F = sup‖f‖ works for every e(0);
- if ‖f(x,t)‖ ≤ α + β‖x‖ with β < k_e, then F = [α + β(sup‖x_m‖ + ‖e(0)‖ + a_max/k_e)]/(1 − β/k_e) works;
- the condition can fail for superlinear f, e.g. the Van der Pol plant of the numerics.'

Replace the proof by:

'(b) While ‖e‖ ≤ r(F), ė = −k_e e + (f − f̂) with ‖f − f̂‖ ≤ F + a_max. Variation of constants gives ‖e(t)‖ ≤ max(‖e(0)‖, (F + a_max)/k_e) < r(F). So x never reaches the boundary of the tube (first exit), stays bounded, and (a) continues the solution.'

### D.3 [minor] PARADIGM.md §5, proof of Theorem D(a), lines 414-415 ('Dobrushin; Golse [82]')

The cited theorems do not cover this system. Golse [82] Thm 1.3.2 and Thm 1.4.3 assume:
- (HK1): a skew-symmetric pairwise kernel;
- (HK2): K ∈ C¹ with globally bounded first derivatives.
Their conclusion is in dist_MK,1 = W₁.

Theorem D differs on every count:
- there is no pairwise kernel; the coupling runs through x;
- the field is only Lipschitz, because of the (·)^+ kink and the ramps;
- its Lipschitz constant grows with ‖e‖ and 1/ε_K;
- it is defined only for s > 0;
- the metric is W₂ with the G⁻¹ cost.

Also, 'Lipschitz on bounded sets, by the proof of A(ii)' covers ϑ for a frozen trajectory only. Round-1 item D.6(a) asked for a Picard argument; the fix only changed the citation. The result itself is true. From R1-5; R2 did not raise it.

**Argument.** I checked the extracted text (scratchpad thmD/golse.txt):
- §1.3, around lines 975-1000, states (HK1)-(HK2) verbatim;
- Theorem 1.4.3, at lines 1586-1596, concludes dist_MK,1(μ₁(t), μ₂(t)) ≤ e^{2L|t|} dist_MK,1(μ₁in, μ₂in).

None of these hypotheses holds for law (L) coupled to the plant.

**Fix.** Replace lines 414-415 by:

'(a) Write the loop as an ODE for (x, Ξ) ∈ ℝⁿ × C(Θ_K; ℝ^p), with the sup norm (G⁻¹ on ℝ^p):

    ẋ = f(x,t) + ẋ_m − k_e(x − x_m) − ∫Φ(x; Ξ(ϑ))dμ_0(ϑ),   ∂_tΞ(ϑ) = v(Ξ(ϑ); x, t),   Ξ_0 = id,

so that μ_t = Ξ_t#μ_0. Extend v Lipschitz to a neighbourhood of Θ_K inside {s > s_min/2} (McShane). On {‖x‖ ≤ R}, with the constants of (c), the right-hand side is Lipschitz:

    ‖v(Ξ(·);x,t) − v(Ξ'(·);x',t)‖_∞ ≤ L_ϑ‖Ξ − Ξ'‖_∞ + L_x‖x − x'‖,

and the x-field changes by at most M_x‖x − x'‖ + L_2‖Ξ − Ξ'‖_∞, because W(Ξ#μ_0, Ξ'#μ_0) ≤ ‖Ξ − Ξ'‖_{L²(μ_0;G⁻¹)}. Picard–Lindelöf gives a unique local solution. Each characteristic stays in Θ_K (Nagumo, as in A(ii)), so the solution continues while x is bounded. Golse [82] Thms 1.3.2/1.4.3 and [77]-[79] are structural analogues (pairwise C¹ kernels, W₁), not the proof.'

### D.4 [minor] PARADIGM.md §5, Theorem D(c) lines 397, 401-402; proof line 422; (e) line 410

One theorem uses two different constants with nearly the same name:
- line 401 defines L_1 := Lip_X f + k_e + sup‖∇_x f̂‖;
- line 402 uses L₁, Lemma B.2's (1 + a_max²J²)^{1/2}, inside L_2.

If L₁ is read as the L_1 just defined, L_2 is not a valid W-Lipschitz constant, and the x-inequality d/dt‖x − x'‖ ≤ L_1‖x − x'‖ + L_2D fails.

Three smaller gaps in the same bullets:
- 'the W-Lipschitz constant' should be 'a';
- the domain of sup‖∇_x f̂‖ is unstated (it should be x ∈ ℝⁿ, μ ∈ P(Θ_K), which gives ≤ a_max J_c);
- Lip_X f must be uniform in t ≤ T.

Merges R1-3 (minor) and R2-6 (cosmetic). I rate it minor because one natural reading makes a constant in the statement false.

**Argument.** Take f = 0, k_e = 0.1, a_max = 0.5, s_min = 1, G = 4I.
- The misread constant is L_2 = 2(0.1 + 0.5e^{−1/2}) = 0.81.
- δ_ϑ against δ_{ϑ+(Δa,0,0)} at x = c gives ‖Δ∫Φ‖/W = |Δa|/(|Δa|/2) = 2.0 > 0.81.
- The intended L_2 = 2(1 + 0.25J²)^{1/2} = 2.22 is valid.

On the experiment's box I recomputed λ_G^{1/2}L₁ = 33.08 against a sampled sup of ‖∇_ϑΦ G^{1/2}‖ = 16.1. The blockwise value is 23.4.

**Fix.** Replace lines 401-402 by:

'- M_x := sup_{t≤T} Lip_X f(·,t) + k_e + a_max J_c, where a_max J_c bounds ‖∇_x∫Φ(x;·)dμ‖ over x ∈ ℝⁿ and μ ∈ P(Θ_K).
- L_2 := λ_G^{1/2}L₁, with L₁ = (1 + a_max²J²)^{1/2} of Lemma B.2, is a W-Lipschitz constant of μ ↦ ∫Φ(x;·)dμ. Indeed, for any coupling, ‖∫Φ(x;·)d(μ − μ')‖ ≤ L_2 times its L²(G⁻¹) cost.'

Also:
- on line 397 set Λ := max(L_ϑ + L_2, L_x + M_x);
- on line 422 write D⁺‖x − x'‖ ≤ M_x‖x − x'‖ + L_2D.

Optional sharpenings:
- L_2 may be replaced by (g_a + max(g_c,g_s)a_max²J²)^{1/2};
- with a one-sided Lipschitz constant of f, k_e can enter M_x with a minus sign.

### D.5 [cosmetic] PARADIGM.md §5, Theorem D(c) lines 399-400 (L_ϑ = O(1/ε_K); L_x 'in (x, e)'); proof lines 420-421

The claim L_ϑ = O(1/ε_K) is correct for the defined two-sided constant, and the 1/ε_K does come from the ramp (h' ≤ 2/ε_K). This is not an error.

The coupling step, however, uses v only through ⟨X − X', G⁻¹(v(X) − v(X'))⟩. So the one-sided constant L_ϑ^+ suffices, and it is ε_K-free. Using it makes Λ uniform as ε_K → 0. This is a recommended but optional strengthening.

Also unstated:
- the domain of L_ϑ and L_x (x ∈ X, t ≤ T);
- that both are proportional to Ē := sup_{X×[0,T]}‖x − x_m(t)‖;
- that L_x must be the Lipschitz constant, Euclidean → G⁻¹, of x ↦ v(ϑ; x, x − x_m(t)). Reading 'in (x, e)' jointly costs a factor √2 or 2.

Merges R1-4 and R2-2. Both reviewers rated it minor; I downgrade it because nothing stated is false.

**Argument.** Recomputed on TransportController.velocity (n = 2, the experiment's box and gains, ‖e‖ = 1), with the Jacobian whitened by G^{±1/2}:

| ε_K | two-sided sup σ_max | × ε_K | one-sided sup λ_max(sym) |
|---|---|---|---|
| 0.2 | 149 | 29.7 | 90 |
| 0.1 | 288 | 28.8 | 88 |
| 0.05 | 628 | 31.4 | 95 |
| 0.02 | 1696 | 33.9 | 96 |
| 0.01 | 3188 | 31.9 | 93 |
| 0.005 | 6771 | 33.9 | 98 |

The two-sided constant is Θ(1/ε_K), as claimed. The one-sided constant is flat in ε_K.

**Lemma** (proved, and checked on 4.8e5 pairs; worst ratio −0.75 ≥ −1). Write N := h(‖θ‖)(θ̂ᵀy)^+θ̂, r = ‖θ‖, c = θ̂ᵀθ̂', α = (θ̂ᵀy)^+. Then

⟨θ − θ', N(θ,y) − N(θ',y)⟩ = (r − r')(h(r)α − h(r')α') + (1 − c)(r'h(r)α + rh(r')α') ≥ −‖y‖‖θ − θ'‖²/(2(θ_max − ε_K)).

The proof uses three facts: h is monotone, |α − α'| ≤ ‖y‖√(2(1 − c)), and ‖θ − θ'‖² ≥ 2√(rr')|r − r'|√(2(1 − c)). The s-box correction is monotone.

**Fix.** **Minimal fix.** Replace lines 399-400 by:

'- L_ϑ is a Lipschitz constant of ϑ ↦ v in the G⁻¹ norm, uniform over x ∈ X and t ≤ T. It is proportional to Ē := sup_{X×[0,T]}‖x − x_m(t)‖ and is O(Ē/ε_K) because of the ramp (h' ≤ 2/ε_K).
- L_x is a Lipschitz constant, Euclidean → G⁻¹, of x ↦ v(ϑ; x, x − x_m(t)), uniform over ϑ ∈ Θ_K, t ≤ T and x ∈ X.'

**Recommended strengthening:**
- Define L_ϑ^+ := sup⟨ϑ − ϑ', G⁻¹(v(ϑ) − v(ϑ'))⟩/‖ϑ − ϑ'‖²_{G⁻¹}.
- Add the lemma above and state L_ϑ^+ ≤ Ē[λ_G H_K + max(g_a/(2(a_max − ε_K)), g_c a_max J_c/(2(c_max − ε_K)))]. The bound uses that P is 1-Lipschitz in y and positively homogeneous blockwise.
- Write the coupling step as d/dt D² ≤ 2L_ϑ^+D² + 2L_x‖x − x'‖D.
- Set Λ := max(L_ϑ^+ + L_2, L_x + M_x), and say that Λ is independent of m and bounded as ε_K → 0.

### D.6 [cosmetic] PARADIGM.md §5, Theorem D(e) lines 407-410; proof lines 427-428

Five wording and rigour points in (e):
- (i) The qualifier 'for every t where W(μ_t,μ*_t) is continuous' is vacuous; that function is continuous everywhere.
- (ii) The hypothesis 'μ*_t is a W-Lipschitz curve in P(Θ_K^{ε})' follows from the other hypotheses.
- (iii) 'The right-hand sides converge uniformly' is imprecise. R^m contains ε^m_max and V^m_*(t), which need not converge, so one must first majorise R^m by R̄^m, which uses ε_max and V_*.
- (iv) sup_m ε^m_max ≤ ε_max is lossy but correct. Applying (e) to tail subsequences already gives limsup_m, and with the actual sups one gets the limit's own sup_{X×[0,T]}‖f − ∫Φdμ*_t‖.
- (v) V_* depends on the approximating ideals.

Merges R1-7, R2-5 and the first half of R2-4. I downgrade R2-4 from minor to cosmetic because of the tail argument.

**Argument.** Each point checks out:
- (i) μ_t is a push-forward by the flow of a bounded field, and μ*_t is a uniform limit of continuous curves.
- (ii) The identity coupling gives W(μ*^m_t, μ*^m_s) ≤ V^m_*|t − s| ≤ V_*|t − s|. This passes to the uniform limit, and Θ_K^{ε} is closed.
- (iv) |ε^m − ε^∞| ≤ L_2 W(μ*^m_t, μ*_t) → 0, and (B) is monotone in ε_max.
- (v) m equally spaced atoms rotating at speed ω on a circle of radius R converge to a stationary measure, yet V^m_* = ωR.

**Fix.** Rewrite (e) as:

'(e) *Passage of (B) to the limit.* In the setting of (d), with X := X_1, let μ*^m_t = (1/m)Σ_jδ_{ϑ*^m_j(t)} satisfy (S4), and let ε^m_max ≥ sup_{X×[0,T]}‖f − ∫Φdμ*^m_t‖. Assume μ*^m_t → μ*_t uniformly in W on [0,T], sup_m ε^m_max ≤ ε_max and sup_m V^m_* ≤ V_*. Only tails matter, so sup_m may be replaced by limsup_m. Then μ*_t is V_*-Lipschitz in W with values in P(Θ_K^{ε}). The limit satisfies (B) in integrated form, and in Dini form for every t ∈ [0,T), with constants (k_e, ε_max, K, L₀, L₁, λ_G, V_*). V_* depends on the approximating ideals; the intrinsic version is Remark B.7.'

Proof:

'(e) Theorem B for each m gives V^m(t₂) − V^m(t₁) ≤ ∫R^m ≤ ∫R̄^m, where R̄^m is R^m with ε^m_max and V^m_*(t) replaced by ε_max and V_*. By (d) and the triangle inequality for W, e^m → e and W(μ^m_t, μ*^m_t) → W(μ_t, μ*_t) uniformly, so V^m → V and R̄^m → R̄ uniformly. Hence V(t₂) − V(t₁) ≤ ∫R̄. R̄ is continuous, which gives the Dini form.'

### D.7 [cosmetic] PARADIGM.md §5, proof of Theorem D(c), lines 420-422

Three rigour gaps in the proof of (c):
- D(t) is never defined.
- 'd/dt D' and 'd/dt‖x − x'‖' do not exist at their zeros.
- The proof does not say that the two solutions share (f, x_m, k_e, G, Θ_K, ε_K).

Also, X denotes both the compact set and the flow. Round-1 item D.6(c) raised these points and was not applied. Merges R1-8 and part of R2-7.

**Argument.** If μ_0 = μ'_0 and x_0 ≠ x'_0, then D(0) = 0 and D(t) ~ ct, so D has no derivative at 0.

The estimate itself is correct. With S = A + D, Ṡ ≤ (L_1 + L_x)A + (L_2 + L_ϑ)D ≤ ΛS.

**Fix.** Replace lines 420-422 by:

'(c) The two solutions share (f, x_m, k_e, G, Θ_K, ε_K). Let γ₀ be an optimal plan for (μ_0, μ'_0), let Ξ_t and Ξ'_t be the two flows, and set D(t)² := ∫‖Ξ_t(ϑ) − Ξ'_t(ϑ')‖²_{G⁻¹}dγ₀. Then (Ξ_t, Ξ'_t)#γ₀ couples μ_t and μ'_t, so W(μ_t,μ'_t) ≤ D(t). Split v_t(Ξ) − v'_t(Ξ') = [v_t(Ξ) − v_t(Ξ')] + [v_t(Ξ') − v'_t(Ξ')] and apply Cauchy–Schwarz:

    d/dt D² ≤ 2L_ϑD² + 2L_x‖x − x'‖D,   D⁺‖x − x'‖ ≤ M_x‖x − x'‖ + L_2D.

Apply Grönwall to D_δ + ‖x − x'‖, with D_δ := (D² + δ²)^{1/2}, and let δ ↓ 0.'

### D.8 [cosmetic] PARADIGM.md §1 (gain k_e, layer width ε_K, Admissible box paragraph) and Theorem D hypothesis line 382

PARADIGM.md never states k_e > 0, 0 < ε_K < min(a_max, c_max) or 2ε_K < s_max − s_min; grep finds none of them. They are needed:
- (b) divides by k_e;
- (a) and (c) need θ_max − ε_K > 0, so that θ̂ is only evaluated away from 0.

The code enforces these conditions in ProjectionOperator and BandwidthProjection. Merges part of R2-7 and part of R1-1.

**Argument.** If ε_K > a_max, then h(0+) = ((ε_K − a_max)/ε_K)² > 0. The correction h(θ̂ᵀy)^+θ̂ then jumps at a = 0, so L_ϑ = ∞ and the proof of (a) breaks. R2 observed that the flow still contracts, so this is a gap in the proof, not a counterexample.

**Fix.** Add to §1 (Admissible box): 'Throughout, k_e > 0, 0 < ε_K < min(a_max, c_max) and 2ε_K < s_max − s_min (as enforced in projection.py).'

### D.9 [cosmetic] PARADIGM.md §0 table row D (line 40); §5 numerics lines 430-431; §7 row D (line 746)

Row D of §0 says 'well-posed mean-field closed loop ... and a global tube bound', marked 'proved'. What is proved:
- (a) gives only local well-posedness;
- global existence and the tube need (b)'s self-consistency condition, which the Van der Pol plant of the D experiment fails.

The numerics paragraph reports the observed W_T/W_0 of 1.64-1.69 without noting that the theorem's factor e^{ΛT} is astronomically larger. The run checks consistency in m (linear dependence on W_0), not the constant.

Merges R1-9 and R2-8. R2's specific claim that results/RESULTS_PARADIGM.md calls 1.68 'the amplification factor e^{Λt}' is REJECTED: the current file holds only the table, and grep finds no e^{Λ...} text.

**Argument.** On the experiment's box, L_2 = λ_G^{1/2}L₁ = 33.08 alone gives e^{ΛT} ≥ e^{33.08·15} = e^{496}, before L_ϑ is counted. The observed W_T/W_0 is 1.64-1.69.

For (b), S(F) − F ≥ 51.6 for every F.

**Fix.** Change §0 row D to: 'local well-posedness of the mean-field closed loop; a Dobrushin estimate; a tube bound and global existence under a self-consistent growth condition (e.g. bounded f); the mean-field limit on the existence interval'.

After line 431 add: 'This is a consistency check only. Here the guaranteed factor e^{ΛT} is ≥ e^{496}. (b)'s condition fails for this plant (S(F) − F ≥ 51.6 for every F), so the run is covered by the localized (d) through the observed boundedness, not by (b).'

Refuter verdicts: holds_with_fixes; holds_with_fixes

## E — verdict: holds_with_fixes

Verdict: holds_with_fixes. Theorem E is true once its tube hypothesis is read as a uniform tube. The core argument checks out: the linearisation, the change of variables to ζ, the uniform exponential stability (UES) step, Khalil Thm 4.13, and the constants c₀ and C. I confirm 6 minor and 5 cosmetic issues and nothing fatal or major. No reviewer issue is rejected outright; I corrected three sub-claims.

**Independent checks.** Scripts are in <scratchpad>\judgeE\ (floquet_m1.py, chirp_and_tube.py, tie_and_norms.py, sliding_gramian.py). No project file was edited; pytest ran with no cache and no bytecode writes.
- **Algebra:** QA + AᵀQ = diag(−k_eI, 0) exactly. The ζ change of variables gives the stated skew-symmetric form, with |B + Bᵀ| = 2e-16. The PE level for m = 1 is 1.39, above αλ_min(G)/m = 0.566.
- **c₀, C and d²:** the norm conversions behind c₀ and C hold on 2000 random states each for m = 2, 3, 4. The closed form d² = (2/m)min‖ϑ*_i − ϑ*_j‖²_{G⁻¹} matches brute force. The collapse (tie) point has V_OT = d²/8 exactly, with two optimal assignments, for m = 2, 3, 4.
- **Reproduced reviewer numbers:** R2's Floquet rates (0.0597; the G-sweep from 3.84e-4 up to 0.210 and down to 1.49e-3; k_e = 0 gives |μ| = 1). R2's ρ-counterexample (|e| = 50 at t = 0.38, 0.51, 0.65). R1's chirp (ratios 0.9715, 0.9942, 0.9980).
- **Gramian windows:** for the rich references, the minimum over sliding 60 s windows is 2.7–5.9× below the t = 0 values that Remark E.1 quotes.
- **Literature:**
  - In JALG 2018 (Jenkins–Annaswamy–Lavretsky–Gibson, SIAM J. Control Optim. 56(4)), Anderson's '[1, proof of Theorem 1]' is used for θ̇ = −uuᵀθ, and they write that 'the rate of convergence is an existence one … No closed form expression is given'.
  - JALG apply [S8, Theorem 5] only to piecewise-smooth PE regressors.
  - The Morgan–Narendra abstract says Theorem 1 characterises uniform asymptotic stability (UAS) of the block system with A + Aᵀ negative definite: UAS iff B(t) is 'sufficiently rich'.
- **Test:** test_theorem_E_linearisation_and_lyapunov_identity passes.

**Adjudication** (reviewer item: decision, merged issue):
- R1-1 domain D: CONFIRMED. Khalil needs a full-state ball; D is a parameter-space set and r is undefined. → M1
- R1-2 tube: CONFIRMED. 'A neighbourhood' of a non-closed orbit, or of the graph, gives no uniform ρ; under the graph reading the statement itself fails. → M1
- R1-3 Corollary E.2: CONFIRMED. The right-hand side does not vanish as the two initial conditions coincide, and C suffices in place of 2C. → M5
- R1-4 dependencies: CONFIRMED for ρ, δ and the phrase 'determined by the linearisation'. Uniformity over the class of Ψ is CONFIRMED as unproved, though it is probably true. → M2
- R1-5 step-3 citations: CONFIRMED. The [S8] abstract supports Thm 1. → M4
- R1-6 8c₀ < d²: CONFIRMED, raised to minor because the step is invalid as written. → M3
- R1-7 notation for z: CONFIRMED. → C1
- R1-8 strictness and unargued conclusions: CONFIRMED. → C1
- R1-9 Remark E.1: CONFIRMED, raised to minor. Its 'decreases with k_e' holds only for large k_e; at m = 1 the rate rises from k_e = 0.3 to 3. → M6
- R1-10 Remark A.1: CONFIRMED. The Scope paragraph of C(iv) also contradicts it. → C3
- R1-11 PE wording: CONFIRMED. → C2
- R2-1 (= R1-6): CONFIRMED. → M3
- R2-2 (= R1-1, R1-2, R1-4): CONFIRMED. Its proposed hypothesis ('ε ≡ 0 for x within ρ of some x_m(t)') is stronger than needed; the moving-tube form suffices. → M1, M2
- R2-3 (= R1-3): CONFIRMED. → M5
- R2-4 (= R1-9): CONFIRMED. → M6
- R2-5 k_e > 0: CONFIRMED; it is stated nowhere in PARADIGM.md. → C4
- R2-6 d² at m = 1: CONFIRMED. → C5
- R2-7 (= R1-11): CONFIRMED. → C2
- R2-8 (= R1-7, R1-8 and the D part of R1-1): CONFIRMED, except the claim that 's ≥ s_min − r' 'contradicts' D ⊂ int Θ_K^ε. That claim is REJECTED: the bound is undefined and looser, not contradictory. → C1, M1
- R2-9 (= R1-5): CONFIRMED, raised to minor. → M4

**UNCERTAIN** (full texts are paywalled):
- whether Anderson (1977) has a theorem for the mixed error equation; its abstract says 'various' equations;
- the exact richness definition and regularity class in [S8];
- whether k and λ are truly uniform over the class of Ψ (probably true, but unproved).

**Statement the repaired proof supports.**
- Assume:
  - there is ρ > 0 with ε(x,t) = 0 whenever ‖x − x_m(t)‖ < ρ, for all t ≥ 0;
  - μ* is constant with ϑ*_j ∈ int Θ_K^ε; put δ := min_j dist(ϑ*_j, ∂Θ_K^ε);
  - k_e > 0;
  - x_m ∈ C¹, with x_m and ẋ_m bounded;
  - (PE) holds.
- Then there are k and λ, from the converse Lyapunov function of the linearisation, and r₀ ≤ min(ρ, δ). With c₀ := r₀²/(2(1 + mλ_G)), which is automatically < d²/8, and C as stated, all three bullets of Theorem E hold.
- The constants depend in addition on ρ and δ.
- Corollary E.2 holds with C in place of 2C, as an exponential-merging statement.

### E.1 [minor] PARADIGM.md §6. Assumptions, bullet 1 (line 456): 'ε ≡ 0 on a neighbourhood of the reference orbit'. Proof step 2 (lines 481–483): 'On a ball D ⊂ int Θ_K^{ε} around ϑ*' … 's ≥ s_min − r'. Step 4 (line 493).

[M1] Merges R1-1, R1-2 and the domain parts of R2-2 and R2-8. The proof needs a fixed Khalil ball, and neither the hypothesis nor D provides one.

(a) Hypothesis. 'A neighbourhood of the reference orbit' does not give a tube {‖x − x_m(t)‖ < ρ} that is uniform in t. ε depends on t, and the orbit of the document's own 6-tone reference is not closed.

(b) Domain.
- D is described as a set in parameter space, and 'D ⊂ int Θ_K^ε around ϑ*' mixes ℝ^{mp} with ℝ^p.
- D has no extent in e, and r is undefined.
- Khalil Thm 4.13 needs f to be C¹ on [0,∞) × {z = (e, ϑ̃) : ‖z‖ < r}, with a Jacobian that is bounded and Lipschitz uniformly in t.
- On that ball the field must equal the smooth closed loop. That requires ‖e‖ < ρ, so that f = f̂*, and ‖ϑ̃_i‖ < δ, so that Proj_K = id.

**Argument.** Under the graph reading, the statement itself fails.
- Let ρ(t) = ρ₀e^{−t} and ε(x,t) = 20(e − ρ(t)e/‖e‖)·1{‖e‖ > ρ(t)}, with e = x − x_m(t).
- This ε is continuous and Lipschitz in x. It vanishes on the open neighbourhood {‖x − x_m(t)‖ < ρ(t)} of the graph. It leaves Ψ, α and every listed constant unchanged.
- For any c₀ > 0, pick t₀ with 2ρ(t₀)² < c₀ and start at e = 2ρ(t₀), ϑ = ϑ*. Then ė(t₀) = (20 − 2k_e)ρ(t₀) > 0 for k_e < 10, and the error escapes.
- With ρ held constant (m = 1, k_e = 3, x_m = 1.5 sin t), |e| reaches 50 at t = 0.38, 0.51 and 0.65 for ρ = 0.1, 0.01 and 0.001. I reproduced these values.

Under the x-only reading, R1's bump construction applies. Bumps on B(p_k, w_k) with p_k → q ∈ cl(O)∖O satisfy the hypothesis, yet inf_t dist(x_m(t), supp ε) = 0, so no Khalil ball exists.

Once the tube and D are fixed, every Khalil hypothesis holds.
- F is C¹ in (t,z), because x_m ∈ C¹.
- The derivatives of Φ through order 3 are bounded for s ≥ s_min and ‖a‖ ≤ a_max.

**Fix.** Assumptions, bullet 1. Replace with: 'There is ρ > 0 such that ε(x,t) = 0 whenever ‖x − x_m(t)‖ < ρ, for all t ≥ 0. μ* is constant with ϑ*_j ∈ int Θ_K^{ε}; put δ := min_j dist(ϑ*_j, ∂Θ_K^{ε}) > 0.'

Step 2, first two sentences. Replace with: 'Let r := min(ρ, δ) and D := {z = (e, ϑ̃) ∈ ℝ^{n+mp} : ‖z‖₂ < r}. For z ∈ D, ‖e‖ < ρ, so f = f̂*; and ‖ϑ̃_i‖ < δ for every i, so ϑ_i ∈ int Θ_K^{ε} and Proj_K = id. On [0,∞) × D the closed loop is F(t,z) = (−k_e e + f̂*(x) − f̂(x), (G∇_ϑΦ(x;ϑ_i)ᵀe)_i), with x = x_m(t) + e. It is C¹ in (t,z), since x_m ∈ C¹, and C^∞ in z.'

Step 2, later sentence. Replace 'for s ≥ s_min − r' with 'for s ≥ s_min'.

Step 4. Append: 'Khalil's r₀ satisfies r₀ ≤ r, so the trajectories it covers stay in D, where the closed loop equals F.'

### E.2 [minor] PARADIGM.md §6, Theorem E: line 464 ('determined by the linearisation below') and line 473 ('The constants depend on m, G, k_e, T, α, sup‖Ψ‖, sup‖Ψ̇‖, a_max and s_min').

[M2] Merges R1-4 and the dependency part of R2-2. The dependency sentence is false, and 'determined by the linearisation' overstates.
- r₀, and hence c₀, also depends on the domain radius min(ρ, δ) and on the Lipschitz constant of the Jacobian on D.
- The list omits ρ and δ. It also omits d² while d² remains in c₀.
- The implied uniformity of k and λ over every Ψ with given (T, α, sup‖Ψ‖, sup‖Ψ̇‖) is not delivered by the cited existence-type theorems.

**Argument.** Dependence on ρ. Setup: m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = diag(2, 1, 0.5), k_e = 3, x_m = 1.5 sin t. I reproduced the following.
- f = f̂* + ε with ε = 20(e − ρ·sign e)·1{|e| > ρ} vanishes on the uniform ρ-tube. So it satisfies even the corrected hypothesis, and it leaves every listed quantity unchanged.
- From ϑ = ϑ* and e = 2ρ (so V_OT = 2ρ²), |e| escapes.
- With ε ≡ 0, the same starts give |e(40)| ≈ 4e-4, 4e-5 and 4e-6.
- So every valid c₀ satisfies c₀ < 2ρ², which tends to 0 as ρ → 0.

Dependence on δ. The proof forces r₀ ≤ δ. R2's runs suggest the true constants do not depend on δ. A quadratic-remainder argument would replace δ by ε_K.

Uniformity. JALG 2018 §3.1 say of Anderson's method that 'the rate of convergence is an existence one … No closed form expression is given relating ρ₂ to the original measures of PE'; I checked this in the text. The uniform claim is probably true, either from explicit uniform-complete-observability (UCO) and output-injection bounds or from compactness of {‖Ψ‖ ≤ b₀, ‖Ψ̇‖ ≤ b₁, PE(T, α)}. It is not proved here.

**Fix.** Line 464. Replace with: 'There are r₀, k, λ > 0 such that the following holds. Here k and λ come from the converse Lyapunov function of the linearisation (steps 3–4), and r₀ depends in addition on r = min(ρ, δ) and on the Lipschitz constant of the Jacobian on D.'

Line 473. Replace with: 'The constants depend on m, G, k_e, T, α, sup‖Ψ‖, sup‖Ψ̇‖, a_max, s_min, ρ and δ.'

To keep the reading that k and λ are uniform over the class, do one of the following:
- cite a quantitative UES estimate for this skew-symmetric class, for example Loría–Panteley, Syst. Control Lett. 47 (2002) 13–24, after checking its hypotheses;
- add the compactness argument.
Otherwise say only that k and λ depend on Ψ(·), k_e, m and G.

### E.3 [minor] PARADIGM.md §6, line 465 (c₀ := min{r₀²/(2(1 + mλ_G)), d²/8}) and proof step 5, line 499 ('≤ 8c₀ < d²').

[M3] Merges R1-6 (rated cosmetic there) and R2-1. The definition of c₀ gives only 8c₀ ≤ d², so the uniqueness step is invalid as written. At V_OT = d²/8, both the uniqueness bullet and the exponential bound fail. The theorem survives only because r₀²/(2(1 + mλ_G)) < d²/8 for every valid Khalil radius, and the proof never says so.

**Argument.** The tie point.
- Take the G⁻¹-closest pair (j,k). Set e = 0 and ϑ_j = ϑ_k = (ϑ*_j + ϑ*_k)/2, with all other particles ideal.
- Then V_id = V_(jk) = V_OT = d²/8, with two optimal assignments. I reproduced this for m = 2, 3, 4 on random configurations, with d² checked by brute force.
- By exchange symmetry and uniqueness of solutions, ϑ_j ≡ ϑ_k for all t, so W(μ_t, μ*) ≥ d/2 for all t.

Why the theorem survives.
- The tie point lies at Euclidean distance² ≤ 2mλ_G·d²/8 from both (0, ϑ*) and (0, ϑ*∘(jk)).
- If r₀²/(2(1 + mλ_G)) ≥ d²/8, that distance² is at most r₀²·mλ_G/(1 + mλ_G) < r₀². The solution would then converge to two distinct equilibria, which is impossible.
- So the d²/8 branch is never active, and 8c₀ < d² does hold.

**Fix.** Define c₀ := r₀²/(2(1 + mλ_G)).

Replace the last sentence of step 5 with: 'If V_π, V_π' ≤ c₀ with π ≠ π', then ‖z − z*_π‖² ≤ 2max(1, mλ_G)c₀ < r₀², and likewise for π'. Khalil's estimate would then send the solution to the distinct equilibria (0, ϑ*∘π) and (0, ϑ*∘π'). Hence π is unique.'

Add a remark: 'At the collapse point (e = 0, ϑ_j = ϑ_k = (ϑ*_j + ϑ*_k)/2 for the G⁻¹-closest pair), V_OT = d²/8, there are two optimal assignments, and ϑ_j ≡ ϑ_k for all t. The same argument shows c₀ < d²/8, and d²/8 is sharp along this direction.'

Minimal alternative: keep the min, but with d²/9 in place of d²/8.

### E.4 [minor] PARADIGM.md §6 proof step 3, lines 490–491: 'Uniform exponential stability then follows from Morgan–Narendra [S8] and Anderson (1977).'

[M4] Merges R1-5 and R2-9 (rated cosmetic there). The load-bearing LTV step rests on citations alone.
- No theorem numbers are given.
- The hypotheses are not matched: richness against (PE), the regularity class of B(t), and the passage from UAS to UES.
- 'Anderson (1977)' without a number risks repeating round 1's citation mismatch. JALG 2018 use Anderson's Theorem 1 for θ̇ = −uuᵀθ.

**Argument.** What the sources say.
- The [S8] abstract: 'we characterize (in Theorem 1) the uniform asymptotic stability' of [[A(t), −B(t)], [B(t)ᵀ, 0]] with A + Aᵀ negative definite, which holds iff B(t) is 'sufficiently rich'. That matches this system, with A = −k_eI.
- JALG (in the proof of their Theorem 4) apply [S8, Theorem 5] only to piecewise-smooth PE regressors, with u and u̇ bounded.

The bound on the derivative is essential.
- Ψ = sin(t²) is bounded and PE: its integral over any 10 s window is ≈ 5.0.
- Yet for ė = −e − Ψζ, ζ̇ = Ψe, the ratio ‖w(t₀+200)‖/‖w(t₀)‖ is 0.9715, 0.9942 and 0.9980 for t₀ = 5, 20 and 50 (reproduced). So that system is not UES.

Richness follows from the Gramian. For a unit vector w, ∫_t^{t+T}‖Ψ_G w‖ ≥ ∫‖Ψ_G w‖²/sup‖Ψ_G‖ ≥ αλ_min(G)/(m·sup‖Ψ_G‖).

The conclusion itself is true: every case I tested is UES numerically.

**Fix.** Replace the last sentence of step 3 with: 'By Morgan–Narendra [S8, Thm 1], which characterises uniform asymptotic stability of [[A, −B(t)], [B(t)ᵀ, 0]] with A + Aᵀ ≺ 0 by richness of B, applied with A = −k_eI and B = Ψ_G (bounded, with bounded derivative), the origin of the w = (e, ζ) system is uniformly asymptotically stable. (PE) supplies the richness, since ∫_t^{t+T}‖Ψ_G w‖ ≥ αλ_min(G)/(m sup‖Ψ_G‖) for every unit w. For linear systems this is UES (Khalil Thm 4.11). UES transfers to (e, ϑ̃) through the constant similarity diag(I, m^{1/2}(I⊗G^{1/2})). The bound on Ψ̇ cannot be dropped: Ψ = sin t² is bounded and PE, but the resulting system is not UES.'

Check the exact richness definition and regularity class against the paper. Either drop 'Anderson (1977)' or cite the specific theorem of that paper that covers the mixed error equation.

### E.5 [minor] PARADIGM.md §6 Corollary E.2, lines 506–508.

[M5] Merges R1-3 and R2-3.
- The title claims uniform-in-time continuity, but the right-hand side measures each initial state's distance to (x_m, μ*). It does not vanish as (x'_0, μ'_0) → (x_0, μ_0).
- What the inequality proves is exponential merging of two solutions in the basin.
- The factor 2 is unnecessary.
- It is left implicit that the two solutions belong to the same closed loop (same m, x_m, μ*).

**Argument.** Take (x'_0, μ'_0) = (x_0, μ_0) with 0 < V_OT(0) ≤ c₀. The left-hand side is identically 0, but the right-hand side is 4Ce^{−λt}(‖e(0)‖ + W(μ_0, μ*)) > 0.

In R2's runs, starts δ₀ = 1e-2, 1e-3 and 1e-4 apart give sup LHS = 4.4e-3, 4.4e-4 and 4.4e-5, while the right-hand side stays at 2C·0.090.

Why C suffices:
- ‖x_t − x'_t‖ = ‖e − e'‖ ≤ ‖e‖ + ‖e'‖;
- W(μ_t, μ'_t) ≤ W(μ_t, μ*) + W(μ'_t, μ*);
- each resulting bracket is at most Ce^{−λt}(…) by Theorem E.

**Fix.** Retitle it 'Corollary E.2 (uniform-in-time exponential merging under PE)' and state: 'For two solutions of the same closed loop with V_OT(0) ≤ c₀ and V'_OT(0) ≤ c₀, ‖x_t − x'_t‖ + W(μ_t, μ'_t) ≤ Ce^{−λt}R₀, where R₀ := ‖e(0)‖ + ‖e'(0)‖ + W(μ_0, μ*) + W(μ'_0, μ*).'

If a continuity statement is wanted, add the following.
- Let δ₀ := ‖x_0 − x'_0‖ + W(μ_0, μ'_0). Note that δ₀ ≤ R₀ ≤ CR₀.
- Apply Theorem D(c) on [0, T₀]; both solutions stay in the closed ρ-tube, which is compact.
- Apply the bound above on [T₀, ∞), with T₀ := (Λ + λ)^{−1} log(CR₀/δ₀).
- Together these give sup_{t≥0}(‖x_t − x'_t‖ + W(μ_t, μ'_t)) ≤ (CR₀)^{Λ/(Λ+λ)} δ₀^{λ/(Λ+λ)}.

### E.6 [minor] PARADIGM.md §6 Remark E.1, lines 502–503 ('The rate is set by αλ_min(G)/m … The 60 s Gramian has λ_min between 1e-3 and 3e-2'). Also experiments/exp_pe_convergence.py line 5 and tests/test_transport.py line 364.

[M6] Merges R1-9 (rated cosmetic there) and R2-4.
- The rate is not set by the PE level. It is non-monotone in G and in k_e, and it collapses for large G.
- α carries a factor of time, so αλ_min(G)/m is not even a rate.
- The quoted Gramian values are single windows starting at t = 0, whereas (PE) needs the infimum over t.

**Argument.** Exact Floquet rates, reproduced to 3 digits. Setup: m = n = 1, ϑ* = (1.5, 0.3, 0.8), x_m = 1.5 sin t, G = gI.
- With k_e = 3 and g = 0.01, 0.1, 1, 10, 100, 1000, the rate is 3.84e-4, 3.91e-3, 4.66e-2, 0.210, 2.32e-2 and 1.49e-3. Over the same range αg grows 10⁵-fold.
- With g = 0.01 and k_e = 0.3, 3, 30, the rate is 1.18e-4, 3.84e-4 and 5.97e-5.
- For large k_e the rate matches λ_min⟨Ψ_GᵀΨ_G⟩/k_e: at g = 1, k_e = 30 it is 5.99e-3 against 6.01e-3.

Gramian windows.
- exp_pe_convergence.gramian integrates over [0, 60) only.
- Over windows starting at t₀ = 0, 30, …, 570, the minimum of λ_min(∫ΨᵀΨ) for the rich seeds 0–2 is 3.2e-3, 9.4e-4 and 3.8e-4. The t₀ = 0 values are 8.6e-3, 5.5e-3 and 1.1e-3.

**Fix.** Replace the two sentences with: 'The rate λ is bounded below by a function of the PE level αλ_min(G)/m, k_e, T and sup‖Ψ_G‖, but it is not set by the PE level. It is not monotone in G or k_e, and it decreases for large G. When k_e is large relative to the reference bandwidth, it is approximately λ_min⟨Ψ_GᵀΨ_G⟩/k_e ≥ αλ_min(G)/(m k_e T). It is small in practice: over sliding 60 s windows, the minimum of λ_min(∫ΨᵀΨ) is about 4e-4 to 3e-2 (results/pe_convergence).'

Also:
- make exp_pe_convergence.py report the minimum over sliding windows;
- reword 'which sets the exponential rate' (exp_pe_convergence.py line 5) and 'The rate is set by the smallest PE eigenvalue' (tests/test_transport.py line 364) to match.

### E.7 [cosmetic] PARADIGM.md §6 proof steps 2–5 (lines 482, 488–499) and the statement bullets (lines 468–471).

[C1] Merges R1-7, R1-8 and the notation part of R2-8.
- (a) z means (e, ϑ̃) in steps 2 and 5 but (e, ζ) in step 3. The proof never says which coordinates Khalil's theorem is applied in, and the stated c₀ and C are correct only for (e, ϑ̃).
- (b) Step 5 writes '≤ r₀²', but Khalil's ball is open.
- (c) Three conclusions are never argued: 'optimal for all large t', 'ϑ_i(t) → ϑ*_π(i)', and the labelled-error limit. 'Labelled error' is not defined in §6.
- (d) 'Continuous in t' understates what Thm 4.13 needs, which is C¹ in (t, z); that does hold.

R2-8's claim that 's ≥ s_min − r' contradicts D ⊂ int Θ_K^ε is rejected. The bound is undefined and looser, not contradictory; it is fixed under M1.

**Argument.** I recomputed both conversions:
- ‖z‖² ≤ 2max(1, mλ_G)V_π;
- ‖z(t₀)‖ ≤ max(1, (mλ_G)^{1/2})(‖e‖ + W).
Both hold on 2000 random states each for m = 2, 3, 4, with worst ratios 0.97 and 0.98.

Strictness: ‖z(t₀)‖² ≤ 2max(1, mλ_G)c₀ < 2(1 + mλ_G)c₀ ≤ r₀².

In w = (e, ζ) coordinates, ‖w‖² = 2V_π exactly, ‖ζ(t₀)‖ = W(t₀) for the optimal π, and W ≤ ‖ζ‖. Applying Khalil there would give c₀ = r₀'²/2 and C = √2k'.

The missing lines are short. Khalil's estimate gives z → 0. For π' ≠ π, V_π'(t) tends to (1/2m)Σ‖ϑ*_π(i) − ϑ*_π'(i)‖² ≥ d²/2, while V_π(t) → 0.

**Fix.** - Name the step-3 state w := (e, ζ). In step 4, add: 'Khalil's theorem is applied in z = (e, ϑ̃).'
- In step 2, replace 'continuous in t' with 'C¹ in (t, z), since x_m ∈ C¹'.
- In step 5, write: '‖z(t₀)‖² ≤ 2max(1, mλ_G)c₀ < 2(1 + mλ_G)c₀ ≤ r₀²'.
- Append to step 5: 'Hence z(t) → 0, i.e. ϑ_i(t) → ϑ*_{π(i)}. For π' ≠ π, V_π'(t) → (1/2m)Σ‖ϑ*_{π(i)} − ϑ*_{π'(i)}‖²_{G⁻¹} ≥ d²/2 > 0, while V_π(t) → 0, so π is the unique optimal assignment for all large t. The labelled error E_lab(t) := ((1/m)Σ‖ϑ_i(t) − ϑ*_i‖²_{G⁻¹})^{1/2} therefore tends to ((1/m)Σ‖ϑ*_{π(i)} − ϑ*_i‖²_{G⁻¹})^{1/2}.'

### E.8 [cosmetic] PARADIGM.md §6, lines 461–462: '(PE) forces the ϑ*_j to be pairwise distinct, with a*_j ≠ 0. Otherwise two column blocks coincide or vanish'.

[C2] Merges R1-11 and R2-7. The claim is true but not sharp, and the stated reason is inaccurate.
- If a*_j = 0, only the n+1 columns for (c,s) in block j vanish; the a-columns φI_n remain.
- (PE) actually forces a stronger condition: the pairs (c*_j, s*_j) must be pairwise distinct.

**Argument.** ∇_ϑΦ = [φI_n, a∇_cφᵀ, a∂_sφ], which becomes [φI_n, 0, 0] when a = 0.

ϑ*_1 = (1, 0.3, 0.8) and ϑ*_2 = (−0.5, 0.3, 0.8) are distinct, but they share (c, s). Their a-columns coincide, so ∫ΨᵀΨ annihilates (v, 0, 0, −v, 0, 0); R2 found two eigenvalues ≈ 1e-16.

**Fix.** Replace with: '(PE) forces a*_j ≠ 0 and pairwise distinct pairs (c*_j, s*_j). If a*_j = 0, the (c,s)-columns of block j vanish. If (c*_j, s*_j) = (c*_k, s*_k), the vector with v in the a-slot of block j, −v in the a-slot of block k and zeros elsewhere lies in the kernel of ∫ΨᵀΨ. In particular, the ϑ*_j are pairwise distinct and d² > 0.'

### E.9 [cosmetic] PARADIGM.md Remark A.1, line 182 ('In this normalisation, all constants below are free of m'). It conflicts with §6 line 473 and with the Theorem C Scope paragraph (line 306).

[C3] R1-10. Remark A.1 says all later constants are free of m. Two parts of the document contradict this:
- Theorem E, whose constants depend on m, as §6 itself states;
- Theorem C(iv), which its own Scope paragraph calls a small-m result.

**Argument.** - Line 473 lists m among Theorem E's dependencies.
- c₀ ∝ 1/(1 + mλ_G).
- λ depends on m through the mp-dimensional Gramian.
- Line 306 gives η̄ < R/(2m) and b = O(1/m) for C(iv).

**Fix.** Replace the sentence with: 'In this normalisation, the constants of Theorems B, C(i)–(iii), D and F are free of m. C(iv) (see its Scope paragraph) and Theorem E are fixed-m statements.'

### E.10 [cosmetic] PARADIGM.md §1 Standing hypotheses (S) and §6 Assumptions.

[C4] R2-5. PARADIGM.md never states k_e > 0. It is used only implicitly, for example in the divisions by k_e in Theorems B and D(b). Theorem E is false for k_e = 0.

**Argument.** For k_e = 0 the (e, ζ) linearisation is purely skew-symmetric, so ‖w‖ is conserved. For m = 1 and G = I, the monodromy matrix has max|μ| = 1.0000000 (reproduced). R2 found the same for m = 3.

**Fix.** Add '(S6) k_e > 0' to the standing hypotheses, or add 'k_e > 0' to the assumptions of Theorem E.

### E.11 [cosmetic] PARADIGM.md §6, line 462 (definition of d²).

[C5] R2-6. For m = 1 the minimum is over the empty set, so d², and c₀ as currently defined, are undefined. A closed form would also make d² > 0 immediate.

**Argument.** Any permutation σ ≠ id moves at least two indices, and the transposition of the G⁻¹-closest pair attains the bound. Hence d² = (2/m)·min_{i≠j}‖ϑ*_i − ϑ*_j‖²_{G⁻¹}. This matched brute force for m = 2, 3, 4.

**Fix.** Write: 'd² := min_{π≠π'}(1/m)Σ‖ϑ*_{π(i)} − ϑ*_{π'(i)}‖²_{G⁻¹} = (2/m)·min_{i≠j}‖ϑ*_i − ϑ*_j‖²_{G⁻¹}, with d² := +∞ for m = 1.' With the M3 fix, d² leaves c₀, but it is still used in the 'all large t' line of C1.

Refuter verdicts: holds_with_fixes; holds_with_fixes

## AB — verdict: holds_with_fixes

Judge verdict: holds_with_fixes. Everything stated for law (L) holds: the A(i) gain map, the scaling identity (max error 9.4e-16 on 24,000 cases), Theorem B step (c), Theorem B, and Corollaries B.3, B.5 and the B.6 inequality. There is one false implication: B.5's pointer to Theorem C(iv).

What I checked:
- I read the target sections of PARADIGM.md, THEORY_REVISED.md Lemma 2, tests/test_transport.py, src/baselines/sigma_mod.py and src/core/projection.py.
- I ran two independent scripts in <scratchpad>\judge\: b5_judge.py and misc_judge.py.
- No project files were edited.

Decisions. R1 and R2 are the two reviewers; duplicates are merged.
1. **B.5 via C(iv)** [R1-1, R2-1]: CONFIRMED, major. I reproduced it independently along (L_σ), with every C(iv) hypothesis asserted and sup W ≤ η̄ measured. On [100, 300] the measured RMS‖e‖ exceeds the B.5 bound by:
   - 8.3× (m=1, b=0.002)
   - 7.2× (m=1, b=0.01)
   - 3.8× (m=1, b=0.05)
   - 6.0× (m=2, b=0.01)
   RMS scales linearly in b. The leak-corrected bound holds in all four runs.
2. **Remark B.4 numbers** [R1-2, R2-2]: CONFIRMED, minor. 533.25 and 382.03 are correct. The range 101–138 reproduces only with a joint centre radius c_max. On the per-kernel set, the classical constant and the per-particle block OT constant are both 214.37 for every m. R2's "true sup ≈ 175" was not checked and is not needed.
3. **A(i) "only if"** [R1-3, R2-9]: CONFIRMED, minor. It is true for m ≥ 2 but unproved. The cited test never activates a joint projection (classical radii are set to 1e6).
4. **Dictionary projection row** [R1-4]: CONFIRMED, minor. The corrected field is not a gradient: its Jacobian is asymmetric. The [46] sub-point is REJECTED, because the row cites [46] only for the exact W₂ projection.
5. **σ- and e-modification rows** [R1-5, R2-7 row 2]: CONFIRMED, minor. sigma_mod.py puts the leak inside the gain.
6. **Remark A.1, "all constants free of m"** [R1-6, R2-6]: CONFIRMED, minor. Theorem E itself says its constants depend on m.
7. **B.6, Zador/Graf–Luschgy for equal-weight atoms** [R1-7, R2-3]: CONFIRMED, minor. The O(m^{−1/p}) conclusion is true for p ≥ 3, but the cited theorem has free weights.
8. **B.6 rate and the curse-of-dimensionality sentence** [R1-8, R2-4]: CONFIRMED, minor. Centroid quantisers give a second-order bound, O(m^{−2/p}); in a 1-D check ε_max/W² ≈ 0.31 is constant. Maurey sampling gives an O(m^{−1/2}) bound whose exponent is free of the dimension.
9. **B.6 for time-varying f** [R1-9, R2-5]: REJECTED. B.6 posits one fixed measure μ*_∞, with any time dependence in ε_∞. A constant μ*_m has V_* = 0, and moving ideals are handled by the explicit speed hypothesis of D(e).
10. **(S) bookkeeping and §1 "C¹ at the inner edge"** [R1-10, R2-8]: CONFIRMED, cosmetic. The s-box ramp is a linear clip. The conditions on ε_K and k_e are not stated in the text, though the code enforces the ε_K condition.
11. **B.3 "=" should be "≤"**: CONFIRMED, cosmetic.
12. **B.5 "sharp" and "measured along the run"**: CONFIRMED, cosmetic.
13. **B.6 Kantorovich–Rubinstein attribution**: CONFIRMED, cosmetic.
14. **OT σ-modification row: missing factor σ and the superdifferential at ties**: CONFIRMED, cosmetic.
15. **Remark B.7**: CONFIRMED, cosmetic. A coupling argument replaces the two-curve formula whose reference is still open.

Also correct the "3.9–5.3×" row (B.6) in review/REVIEW_LOG.md.

### AB.1 [major] PARADIGM.md §3, Corollary B.5, hypothesis 'a quantity measured along the run, or controlled a priori by Theorem C(iv) or Theorem E'

B.5 integrates (B), which holds only along (L). Theorem C(iv) controls W only along (L_σ); its gain condition 8K²λ_G²η̄² ≤ k_eσ_m forces σ_m > 0. The Lyapunov inequality of (L_σ) is C(iii), which carries the extra leak term −σ_mη² + σ_mηC^{1/2}, positive when b > η. Applied along (L_σ), the B.5 bound is false: the leak causes a steady tracking bias of order b, while B.5 predicts order η̄² = O(b²). B.5 is correct for (L), including with η̄ taken from Theorem E.

**Argument.** Judge's script scratchpad/judge/b5_judge.py.

Setup: a_max=1, c_max=2, s_min=1, s_max=2, ε_K=0.15, G=I, k_e=4, σ_m=1. The plant is representable (ε_max=0) and the ideal is constant (V_*=0). The run starts at the ideal with e(0)=0. The prior is the ideal with its weights shifted by b. η̄ = 1.02·(2β/α)^{1/2}.

All three C(iv) hypotheses were asserted, and sup_t W ≤ η̄ was measured. Integration: RK4, dt = 2e-3, window [100, 300]. Measured RMS‖e‖ against the B.5 bound:
- m=1, b=0.002: 3.65e-4 vs 4.4e-5 (8.3×)
- m=1, b=0.01: 1.82e-3 vs 2.54e-4 (7.2×)
- m=1, b=0.05: 9.12e-3 vs 2.41e-3 (3.8×)
- m=2, Δ_gap=5.3 with the gap condition met, b=0.01: 1.60e-3 vs 2.66e-4 (6.0×)

RMS scales linearly in b. Both reviewers' independent runs agree: R1 on the C(iv) test instance, R2 with ratios of 4.9–11.6×.

Mechanism: at the leaky steady state, G∇g ≈ σ_m(ϑ−ϑ⁰) forces e = O(b).

**Fix.** Replace the clause with: 'along a solution of (L); η̄ can be measured in oracle simulations or controlled a priori by Theorem E'.

Add the following variant.

*Variant for (L_σ).* Along a solution of (L_σ) with sup_{[T₁,T₂]} W(μ_t,μ*_t) ≤ η̄ and 4η̄R < Δ_gap (for example, η̄ from Theorem C(iv)):
RMS_{[T₁,T₂]}(‖e‖) ≤ D̄/(2k_e) + ( D̄²/(4k_e²) + (V_*η̄ + σ_m b²/4)/k_e + V_OT(T₁)/(k_e(T₂−T₁)) )^{1/2}.

Proof: apply C(iii) with C = b(t)² ≤ b², and use −σ_mη² + σ_mηb ≤ σ_m b²/4.

Without the gap condition, replace b² by b² + 4η̄R. The limit of the bound is ε_max/(2k_e) + (ε_max²/(4k_e²) + σ_m b²/(4k_e))^{1/2}, so the leak bias is first order in b. When b = 0 and 4η̄R < Δ_gap, the leak term is ≤ 0 and B.5 holds unchanged.

The corrected bound holds in all four runs above: 5.0e-4, 2.55e-3, 1.36e-2 and 2.55e-3.

### AB.2 [minor] PARADIGM.md §3, Remark B.4, first bullet ('533, against 101–138 … ≈ 382'); also review/REVIEW_LOG.md row B.6 ('3.9–5.3×')

533 and 382 are correct arithmetic. The classical range 101–138, however, is Lemma 2 with W_max = a_max/√m, which circumscribes the per-kernel product, combined with a joint centre radius c_max instead of √m·c_max. That centre set is √m times smaller than the per-kernel product that A(i) and Remark A.1 prescribe. With consistent bookkeeping the two global constants coincide. Also, 382 is not the best block-structured constant.

**Argument.** Judge's script misc_judge.py. Test box: a_max=8, c_max=4, s_min=0.3, s_max=3, G=(2,1,0.5), F_max=a_max. Then J=3.1784, L₁=25.447, λ_G=2, Δ=14.373.

OT constants:
- 533.25: crude form L₀ + L₁λ_G^{1/2}Δ
- 382.03: merged J, through W
- 322.95: J_c and J_s separate, through W
- 246.69: per particle, merged J
- 214.37: per particle with J_c and J_s separate, 4a_max + a_max(2c_max J_c + (s_max−s_min)J_s)

Classical Lemma 2:
- With joint centre radius c_max: 214.37 (m=1), 137.80 (m=6), 101.15 (m=64).
- On the per-kernel product (centre radius √m·c_max): 214.37 for every m.

The round-1 log already contains the m-uniform version (review/REVIEW_round1.md line 351) next to the inconsistent one (lines 133 and 191).

**Fix.** Replace the bullet with the following text:

'Globally the OT constant is not smaller than the classical one. For the test box (a_max = 8, c_max = 4, s_min = 0.3, s_max = 3, G = diag(2, 1, 0.5), F_max := a_max), on the same per-kernel admissible set, the classical Lemma 2 constant is F_max + 3a_max + 2a_max c_max J_c + a_max(s_max − s_min)J_σ = 214 for every m. The OT linear branch is 533 in the form L₀ + L₁λ_G^{1/2}Δ and 382 with the block split through W. It is 214, identical to the classical constant, when |R_i| ≤ ‖e‖(2a_max + ‖ã_i‖ + a_max(J_c‖c̃_i‖ + J_s|s̃_i|)) is bounded per particle. OT gains only locally (Corollary B.5).'

Also correct the '3.9–5.3×' row in review/REVIEW_LOG.md.

### AB.3 [minor] PARADIGM.md §2, Theorem A(i) 'if and only if'; last two sentences of the proof of (i); §7 verification-map row 'A(i) projected rates, and sharpness'

The 'only if' direction has only a heuristic sentence as support. test_identification_projected_rates_per_kernel sets the classical radii to 1e6, so the classical projection is inactive, and applies per-kernel projections by hand. No joint projection is ever compared; the test's only negative check is the unrescaled layer. The direction is true for m ≥ 2. For m = 1 the two operators coincide when W_max = a_max and the layer is ε_K. The scaling-identity proof cites only the ramp identity.

**Argument.** Proof (R1's two-case argument, verified). Let r_in = (a_max−ε_K)/m. Let the joint ball have radius W_max and layer ε_W, and set ρ = W_max−ε_W.
- Case ρ < √2·r_in: take two rows of equal norm r with √2r ∈ (ρ, min(√2r_in, W_max)], parallel to e. The per-kernel field is inactive and the joint field is active.
- Case ρ ≥ √2·r_in: take one row with norm in (r_in, min(a_max/m, ρ)). The per-kernel field is active and the joint field is inactive.

Explicit example (misc_judge.py): m=2, n=1, ϑ₁=(a_max,0,1), ϑ₂=(0,0,1), x=0.3, e=1, g_a=2. Law (L) gives ȧ = (0, 1.912). Joint balls with (W_max, ε_W) = (0.707, 0.05), (0.707, 0.3), (1.0, 0.6) and (0.9, 0.8) give m·Ẇ₁ = 1.912, 1.729, 1.859 and 1.434 respectively.

The scaling identity itself was confirmed to 9.4e-16.

**Fix.** Proof and statement:
- Write 'for m ≥ 2' and add the two-case argument to the proof; the same argument applies to vec(c).
- In the scaling-identity proof add: 'the direction a(aᵀy)/‖a‖² is homogeneous of degree 0 in a and 1 in y, and the activation aᵀy > 0 is scale-invariant'.

Test:
- Add an assertion to the test that compares against ProjectionOperator.project_matrix on an active joint ball, using the m = 2 example above.
- Until that assertion exists, change the proof sentence and the verification-map row to 'checks the identity and that the layer must be rescaled'.

### AB.4 [minor] PARADIGM.md §2, A(iii) dictionary row 'projection onto Θ_K | ε_K-regularised support constraint', and the proof of (iii)

The proof of (iii), 'the rows follow from ∇_{G⁻¹}V = G∇V', covers only the potential-energy rows. For ε_K > 0 the Pomet–Praly correction depends on the driving velocity: it acts only when ϑᵀGy > 0. Its G⁻¹-rescaling is not a gradient field in ϑ. So the projected law is not −G∇(δR/δμ) for any functional R, even a time-dependent one, and 'regularised support constraint' has no proved variational meaning. (R1's sub-point on [46] is rejected: the row cites [46] only for the exact W₂ projection.)

**Argument.** Judge's script misc_judge.py: m=1, n=1, G=I, x=0.7, e=1, ϑ=(0.93, 0.2, 1.5), with a inside the layer and c, s interior. In the Jacobian of G⁻¹(Proj_K(ϑ,G∇g) − G∇g), the c-derivative of the a-component is −0.060 while the a-derivative of the c-component is 0. The asymmetry of 0.060 means the field is not a gradient.

Analytically, the a-correction depends on c through φ(x;c,s), while the c-correction vanishes when h(‖c‖) = 0.

**Fix.** Relabel the row: 'velocity-level (Pomet–Praly) regularisation of the support constraint μ_t ∈ P(Θ_K); for ε_K > 0 not the gradient flow of any functional of μ ([46] Prop. 3.1 for the exact W₂ projection)'.

In the proof of (iii), write 'The potential-energy rows follow from …'.

Optionally add, marked as formal: as ε_K → 0 the field tends to the G⁻¹-projection of G∇g_t onto the tangent cone of Θ_K, which is the velocity of the W_{G⁻¹}-gradient flow of F_t + ι_{P(Θ_K)}.

### AB.5 [minor] PARADIGM.md §2, A(iii) dictionary rows 'σ-modification toward ϑ⁰' and 'e-modification'

These rows give the leak −σ(ϑ−ϑ⁰) placed outside the gain, which is the leak used in (L_σ) and Theorem H. The classical device as implemented in src/baselines/sigma_mod.py is Ẇ = Γ_W(φeᵀ − σW), with the leak inside the gain (the Ioannou form). Under the A(i) dictionary it maps to a different, block- and m-dependent object, so the table pairs the classical device with the wrong measure-level object.

**Argument.** Apply W_i = a_i/m, γ_W = g_a/m, γ_c = m·g_c and γ_σ = m·g_s to the classical leaks. They become ȧ_i ⊃ −σ(g_a/m)a_i, ċ_i ⊃ −σ·m·g_c·c_i and ṡ_i ⊃ −σ·m·g_s·s_i.

That is the velocity −σGD(ϑ−ϑ⁰) with D = diag(I_n/m, mI_n, m), the G⁻¹-gradient of (σ/2)∫(ϑ−ϑ⁰)ᵀD(ϑ−ϑ⁰)dμ. In the G⁻¹ metric it is λ-convex with λ = σ·min(g_a/m, m·g_c, m·g_s). It equals −σ(ϑ−ϑ⁰) only when GD = I.

**Fix.** Rename the row 'σ-modification in the G⁻¹ metric (leak outside the gain, ϑ̇ = G∇g − σ(ϑ−ϑ⁰), as in (L_σ) and H)'.

Add one line: 'The classical Γ-scaled leak Γ(y − σθ) corresponds under A(i) to (σ/2)∫(ϑ−ϑ⁰)ᵀdiag(I_n/m, mI_n, m)(ϑ−ϑ⁰)dμ, whose convexity modulus σ·min(g_a/m, m·g_c, m·g_s) depends on m.'

Make the same change for the e-modification row, with coefficient σ‖e‖.

### AB.6 [minor] PARADIGM.md §2, Remark A.1 ('all constants below are free of m'; 'fixed γ's corresponds to g_c, g_s ∝ 1/m')

The claim that all later constants are free of m is contradicted elsewhere in the document:
- Theorem E states that its constants depend on m: c₀ and C contain mλ_G, and Remark E.1 gives a rate of order αλ_min(G)/m.
- The Scope paragraph of C(iv) gives Δ_gap ≤ 2R²/m and η̄ < R/(2m).
- ε_max and V_* are data of the m-atom ideal; B.6 writes ε_max(m).
The last sentence of the remark also omits that g_a = mγ_W ∝ m.

**Argument.** The contradicting statements are quoted from Theorem E, Remark E.1 and the Scope paragraph of Theorem C.

With fixed classical γ's, G = diag(mγ_W, γ_c/m, γ_σ/m). Then λ_G ∝ m, and Δ² = (2a_max)²/(mγ_W) + (2c_max)²m/γ_c + (s_max−s_min)²m/γ_σ ∝ m.

**Fix.** Replace the claim with: 'In this normalisation the structural constants k_e, K, L₀, L₁, λ_G and Δ, and hence the constants of Theorems B, D and F, are free of m for fixed G and a_max. ε_max and V_* are data of the chosen m-atom ideal, and the local constants of C(ii), C(iv) and E depend on m.'

Replace the last sentence with: 'The classical convention of fixed γ's corresponds to g_a = mγ_W ∝ m and g_c, g_s ∝ 1/m. Weights speed up, centres and widths freeze, and λ_G ∝ m while Δ ∝ √m.'

### AB.7 [minor] PARADIGM.md §3, Corollary B.6, 'With μ*_m an optimal m-point quantiser, W(μ*_m, μ*_∞) = O(m^{−1/p}) for absolutely continuous μ*_∞ (Zador; Graf–Luschgy [S9])'

The m-atom ideals of §1 have equal weights 1/m. Optimal quantisers in the sense of Zador and Graf–Luschgy carry Voronoi-cell weights, so the cited theorem is used outside its hypotheses. The conclusion is still true: for compactly supported μ*_∞ and p ≥ 3 > 2, equal-weight approximations achieve O(m^{−1/p}) in W₂. Absolute continuity is what makes the order unimprovable (the lower bound); the upper bound does not need it. Atoms chosen as cell barycentres lie in the convex set Θ_K^ε, so admissibility is automatic.

**Argument.** Graf–Luschgy Lemma 3.4 identifies e_{m,2} with inf W₂ over measures with at most m atoms and free weights.

Weights matter (R2's example): for p=3, μ = q·U(box₁) + (1−q)·U(box₂) with q = 1/√2 and m = 2. The free-weight optimum has W₂ = 0.172, while every equal-weight measure has W₂ ≥ 0.670.

Equal-weight upper bound: an equal-mass dyadic decomposition gives W₂² ≲ Σ_j 2^{pj}4^{−j}/m + 4^{−J} = O(m^{−2/p}) for p > 2.

**Fix.** Replace the sentence with: 'For compactly supported μ*_∞ and p ≥ 3 there are equal-weight μ*_m, with atoms (cell barycentres) in Θ_K^ε, such that W(μ*_m, μ*_∞) = O(m^{−1/p}) (equal-mass dyadic decomposition; J. Chevallier, J. Appl. Probab. 55 (2018) 1037–1045; for p > 4 also Fournier–Guillin 2015 for i.i.d. atoms, in expectation). By Zador's theorem (Graf–Luschgy Thm 6.2 with Lemma 3.4) the order cannot be improved, even with free weights, when μ*_∞ has a non-trivial absolutely continuous part.'

### AB.8 [minor] PARADIGM.md §3, Corollary B.6, last sentence ('the curse of dimensionality sits in the number of atoms … and in no Lyapunov constant'); §8 item 7

The m^{−1/p} rate belongs to the first-order Lipschitz/W bound, not to ε_max(m). For a centroid μ*_m the first-order term cancels, so ε_max is second order in W. This covers in particular a best equal-weight approximation, whose atoms are barycentres of their transport cells. Separately, i.i.d. sampling (Maurey) gives an existence bound whose exponent does not depend on the dimension. The sentence therefore presents an artefact of the proof route as intrinsic.

**Argument.** Let γ = Σγ_j be an optimal plan with ϑ_j the barycentre of γ_j. Then ∫Φdμ*_∞ − (1/m)ΣΦ(ϑ_j) = Σ_j∫[Φ(ϑ) − Φ(ϑ_j) − ∇Φ(ϑ_j)(ϑ−ϑ_j)]dγ_j. By (1) and convexity of Θ_K^ε, ε_max(m) ≤ ε̄_∞ + min{Kλ_G W², L₁λ_G^{1/2}W} = ε̄_∞ + O(m^{−2/p}).

1-D check (misc_judge.py; uniform centres on [−1,1], midpoint atoms): ε_max/W₂² = 0.31, 0.31, 0.30, 0.31 and 0.32 for m = 4 to 64, while ε_max/W₂ halves with each doubling of m.

Maurey: draw i.i.d. atoms from μ*_∞, which lie in Θ_K^ε. Hoeffding per coordinate, a net of X at scale m^{−1/2} and Lipschitz continuity in x give, with positive probability, ε_max(m) ≤ ε̄_∞ + O(a_max(n log(m·diam X)/m)^{1/2} + a_max/(s_min√m)).

R2's simulations at n=3 (p=7) show the W-route bound overshooting by about 300×.

**Fix.** After the W bound, add: 'If μ*_m is a centroid quantiser for an optimal plan (in particular a best equal-weight approximation), then ε_max(m) ≤ ε̄_∞ + min{Kλ_G W², L₁λ_G^{1/2}W} = ε̄_∞ + O(m^{−2/p}). Independently of p, i.i.d. atoms from μ*_∞ give some m-atom ideal with ε_max(m) ≤ ε̄_∞ + O((n log m/m)^{1/2}) (Maurey).'

Replace the last sentence with: 'The exponent 1/p is a property of controlling ε_max through W (as Theorem D(e) requires), not of the approximation error itself.'

In §8 item 7, target the second-order quantisation error and rates that do not go through W.

### AB.9 [cosmetic] PARADIGM.md §1: projection paragraph ('C¹ at the inner edge of the layer'), standing hypotheses (S1)–(S5) and the sentence after them; Corollary B.3 limsup bullet

Five bookkeeping points:
(a) 'Under (S1)–(S3) … each V_π is C¹' also needs (S4), because V_π contains ϑ*_{π(i)}(t).
(b) Discharging (S5) by Theorem D(b) means taking X ⊇ ∪_t𝒯(t). Then ε_max, and every constant of B.3 and B.5, depends on ‖e(0)‖.
(c) The limsup in B.3 needs (S5) on [0, ∞).
(d) The conditions k_e > 0 and 0 < ε_K < min{a_max, c_max, (s_max−s_min)/2} are never stated. Without the second, Θ_K^ε can be empty.
(e) The s-box ramp clip((s−s_min)/ε_K, 0, 1) is piecewise linear, so Proj_K is not C¹ at the inner edge for s. This is harmless, since the proofs use only Lipschitz continuity.

**Argument.** Direct reading of the text. In src/core/projection.py, BandwidthProjection (lines 169–185) uses a linear clip ramp and raises a ValueError when 2ε ≥ σ_max − σ_min. The tube 𝒯 of Theorem D(b) has a radius that contains ‖e(0)‖.

**Fix.** Make these edits:
- Write '(S1)–(S4)' in the sentence after the hypotheses.
- Add (S0): k_e > 0, G ≻ 0 block-isotropic, s_min > 0, 0 < ε_K < min{a_max, c_max, (s_max−s_min)/2}.
- In (S5), add: 'For bounded f, take X := closure of ∪_t𝒯(t) from Theorem D(b); ε_max is the sup over this X and depends on ‖e(0)‖.'
- In the limsup bullet of B.3, add 'under (S5) on [0, ∞)'.
- In §1, restrict 'C¹ at the inner edge' to the ball operators.

### AB.10 [cosmetic] PARADIGM.md §3, Corollary B.3, 'using V_max^{OT} = Δ²/2'

Only W ≤ Δ is known, so the parameter part of V_OT is bounded by Δ²/2; it is not equal to it.

**Argument.** Because μ*_t ⊂ Θ_K^ε, sup ½W² can be strictly smaller than Δ²/2. The limsup argument needs only an upper bound.

**Fix.** Write 'using V_max^{OT} ≤ Δ²/2'.

### AB.11 [cosmetic] PARADIGM.md §3, Corollary B.5 title ('a-posteriori sharp bound') and 'a quantity measured along the run'

Nothing shows that the adaptive loop attains ε_max/k_e. That value is attained by the frozen-parameter channel ė = −k_e e + ε. Also, W(μ_t, μ*_t) requires the unknown μ*_t, so it can be measured only in oracle simulations.

**Argument.** With e ≠ 0 the adaptive law keeps moving the parameters. So η̄ → 0 and RMS → ε_max/k_e are not shown to occur together.

**Fix.** Change the title to 'Corollary B.5 (conditional bound)'. Write 'a quantity available in oracle simulations (it requires μ*_t)'. Write 'which tends to ε_max/k_e, the value of the frozen-parameter error channel, as η̄ → 0 and T₂ → ∞'. Optionally add verification-map rows for B.5 (a window check along (L)) and B.6.

### AB.12 [cosmetic] PARADIGM.md §3, Corollary B.6, 'Kantorovich–Rubinstein gives W₁ ≤ W₂ for the same ground metric'

The citation is attached to the wrong step. W₁ ≤ W₂ is Jensen's inequality. Kantorovich–Rubinstein duality concerns scalar Lipschitz functions, and Φ is ℝⁿ-valued. The bound itself is correct.

**Argument.** ‖∇_ϑΦ‖_op ≤ (φ² + a_max²‖∇_{(c,s)}φ‖²)^{1/2} ≤ L₁, and ‖v‖ ≤ λ_G^{1/2}‖v‖_{G⁻¹}. A direct coupling argument gives the inequality. The block constant (g_a + a_max²(J_c²g_c + J_s²g_s))^{1/2} is 21.4 for the test box, against L₁λ_G^{1/2} = 36.0.

**Fix.** Replace the phrase with: 'For any coupling γ, ‖∫Φ d(μ*_m − μ*_∞)‖ ≤ ∫‖Φ(ϑ) − Φ(ϑ')‖dγ ≤ L₁λ_G^{1/2}∫‖ϑ − ϑ'‖_{G⁻¹}dγ ≤ L₁λ_G^{1/2}(∫‖ϑ − ϑ'‖²_{G⁻¹}dγ)^{1/2}; minimise over γ.' Optionally use the block constant.

### AB.13 [cosmetic] PARADIGM.md §2, A(iii) dictionary row 'OT σ-modification (Theorem C)'

The W-gradient of (σ/2)W(·,μ⁰)² at atom i is σ(ϑ_i − ϑ⁰_{τ(i)}); the row drops the factor σ. At ties, ½W² is a minimum of smooth functions with distinct gradients. Its Fréchet subdifferential is then empty, so the measurable selection comes from the superdifferential.

**Argument.** Let f = min(f₁, f₂) at a tie with ∇f₁ ≠ ∇f₂. Testing a subgradient p along both d and −d forces ⟨∇f₁ − ∇f₂, d⟩ = 0 for every d, which is a contradiction.

**Fix.** Rewrite the row: '(σ/2)W(μ,μ⁰)²; for atomic μ with a unique optimal assignment its W-gradient is σ(ϑ_i − ϑ⁰_{τ(i)}); at ties (L_σ) uses a measurable selection from the superdifferential {σ(ϑ_i − ϑ⁰_{τ(i)}) : τ optimal}'.

### AB.14 [cosmetic] PARADIGM.md §3, Remark B.7 (continuum form)

Only a one-sided bound is needed, and it follows from a coupling argument, which closes the open 'two-curve formula' reference. Three smaller points:
- 'sup_t' should be 'ess sup_t', because w_t is defined only for a.e. t.
- The smallest admissible V_* is the metric derivative of μ*_t.
- μ_0 ∈ P(Θ_K) should be stated, since (P) needs it.

**Argument.** Push an optimal plan γ_t forward by (X_{t,t+h}, Y_{t,t+h}). Here X is the flow of v, and Y is the flow of w, or the superposition measure of AGS Thm 8.2.1 when w is not Lipschitz. This gives a coupling of (μ_{t+h}, μ*_{t+h}), so D⁺(½W²) ≤ ∫⟨ϑ−ϑ', v_t(ϑ) − w_t(ϑ')⟩_{G⁻¹}dγ_t. The bound holds at every t when both fields are continuous in t, and at a.e. t in general. Integrating (P) and Lemma B.2 against γ_t then gives (B).

**Fix.** Replace the appeal to the two-curve formula with this coupling argument. Define V_* := ess sup_t ‖w_t‖_{L²(μ*_t;G⁻¹)}; its smallest value is the metric derivative |μ̇*_t|_W (AGS Thm 8.3.1). Add μ_0 ∈ P(Θ_K).

Refuter verdicts: holds_with_fixes; holds_with_fixes

## F — verdict: holds_with_fixes

VERDICT: holds_with_fixes. Theorem F (i)–(iii) is true along (L). Its (L_σ) clause is true for the OT leak, and tight: on a pure-leak run the ratio is 1.0000. The proof is valid once the cited identity is replaced by a block-wise version with a positive part. After merging duplicates there are no major issues: 3 minor and 6 cosmetic.

RECOMPUTATION for the norm-ball operator as implemented. P = Gy − hGθ(θᵀGy)/(θᵀGθ) when h > 0 and θᵀGy > 0; otherwise P = Gy.
- P = Gq with q = y − hθ(θᵀGy)/(θᵀGθ), so ‖P‖²_{G⁻¹} = qᵀGq = ‖y‖²_G − h(2−h)(θᵀGy)²/(θᵀGθ). The constant h(2−h) is right: the sympy residual is 0 for general symmetric G, and the numerical error is 2e-15 on 14899 active states.
- The identity is valid everywhere only with a positive part, ‖y‖²_G − h(2−h)[(θᵀGy)₊]²/(θᵀGθ) (error 2.6e-15). Without the positive part it is off by up to 100% of ‖y‖²_G on 15011 states with h > 0 and θᵀGy < 0.
- The isotropic path that TransportController uses equals the general-G path for G = gI (error 3.6e-15).
- The box operator multiplies g_s·y_s by a linear ramp r and acts only when y_s points outward. That is the same form with n_s = ∓1 and h_s = 1 − r (error 5.8e-16).
- For Proj_K (factor-wise, with G block-diagonal) the block-sum identity holds to 6.4e-16 on 30480 active particle-states, and the largest ratio ‖Proj_K‖²_{G⁻¹}/‖y‖²_G is 1.000000.
- The velocity bound (i) has worst ratio 0.98 over 6000 adversarial states. The other pointwise bounds hold too: ‖∇g‖ ≤ L₁‖e‖ (worst 0.9987), ‖∇g‖_G ≤ L_G‖e‖ (0.9991), and the (L_σ)-OT pointwise clause (0.9992).
- Along (L), (ii) holds on all sampled pairs and (iii) on all windows, with worst ratio 0.034 (seed 0) and 0.032 (seed 1).

ADJUDICATION, reviewer 1:
- R1-1, identity wrong (no positive part, factor-wise operator, box form): CONFIRMED, minor. Both literal counterexamples reproduce. Merged with R2-1.
- R1-2, the Reading holds only along (L): CONFIRMED, minor. Reproduced exactly: sup‖e‖ = 7.1e-14 and ∫‖e‖ = 4.1e-13, yet W(μ_10,μ_0) = 1.0535, whereas the Reading requires ∫‖e‖ ≥ 0.0293.
- R1-3, C.5 applies F outside its hypotheses: CONFIRMED, minor. In the δ=0 labelled example e ≡ 0 (1.5e-13) while the path length is 1.000 and W = 0.498; F's (L_σ) clause fails by 39.7× on sampled pairs. On the experiment's real seed-3 trajectory it fails by 1.6× (integrated) and 1.75× (pointwise), and holds with d_lab (0.77).
- R1-3's figure '127× at the relabelled start of exp_sigma_mod_transient.py' is REJECTED as stated. It used a hand-set ‖e‖ = 1.4e-4 (x = x_ref + 1e-4·1). The experiment's real start has ‖e(0)‖ = 0.353, where the t=0 ratio is 0.05.
- R1-4, hypotheses not stated: CONFIRMED, cosmetic. The standing hypotheses (S) partly cover this.
- R1-5, gradient bound cited from inside a proof: CONFIRMED, cosmetic.
- R1-6, loose constant and W_∞ strengthening: CONFIRMED, cosmetic, optional. For the test box L_G = 21.356 against 35.988. The constant is not uniformly loose: with g_a = λ_G and small a_max it is reached to 99.9%.
- R1-7, test slack: CONFIRMED, cosmetic. Ratios are 0.034 and 0.029, and in seed 0 the projection is active at 0 of 300 samples.
- R1-8, box not C¹ at the inner edge: CONFIRMED, cosmetic. The one-sided derivatives are −6.667 and 0; the ball operator is C¹ there.
- R1-9, (L_σ) bookkeeping, Filippov case, AGS metric: CONFIRMED, cosmetic.
- R1-10, notation: CONFIRMED, cosmetic. A grep finds no definition of RMS in PARADIGM.md.
- R1's sub-fix 'edit REVIEW_round1.md': REJECTED. That file holds the raw reports, so log the erratum in review/REVIEW_LOG.md instead.

ADJUDICATION, reviewer 2:
- R2-1: CONFIRMED, duplicate of R1-1. Counterexample (b) would need h(2−h) = 1.25 > 1.
- R2-2: CONFIRMED, duplicate of R1-2 and R1-3.
- R2-3, test slack: CONFIRMED, severity lowered from minor to cosmetic, because nothing stated is false and F is proved.
- R2-4: CONFIRMED, cosmetic, duplicate of R1-8.
- R2-5: CONFIRMED, cosmetic, duplicate of R1-4, R1-5, R1-9 and R1-10.
- R2's side figures ('20×' in a random (L_σ) run, 'worst ratio 0.81'): UNCERTAIN as exact numbers but immaterial, since the exact e ≡ 0 examples cover them.

The judge scripts are j1_identity.py, j2_velocity.py and j3_dynamics.py in <scratchpad>\judge_F\. No project source or document was edited.

### F.1 [minor] PARADIGM.md l.320–321 (proof of Theorem C(i)), cited by the proof of Theorem F at l.444. Operators: src/core/projection.py l.85–114 (ball) and l.175–187 (box), applied factor-wise by TransportController.project in src/core/transport.py l.204–211.

Merges R1-1 and R2-1; CONFIRMED. F(i) relies on the display ‖Proj(ϑ,Gy)‖²_{G⁻¹} = ‖y‖²_G − h(2−h)(ϑᵀGy)²/(ϑᵀGϑ). That display is exact only for one norm-ball factor, and only on its active set {h>0, θᵀGy>0}. Read as written, for Proj_K and the lifted ϑ, it is false in three ways:
(a) It has no positive part. When θᵀGy ≤ 0 the operator returns Gy unchanged.
(b) Proj_K works factor by factor: an a-ball, a c-ball and an s-box, each with its own ramp. No single h with the whole-ϑ quotient reproduces it.
(c) The box factor has a linear ramp and acts only on outward y_s, with normal −1 at the lower face. It is not the ball formula with θ = s.
The inequality F actually needs, ‖Proj_K(ϑ,Gy)‖_{G⁻¹} ≤ ‖y‖_G, is true. It holds only because G = diag(g_aI, g_cI, g_s) is block-diagonal in (a,c,s) and every ramp lies in [0,1]. The proof does not state either fact.

**Argument.** Recomputation: P = Gq with q = y − hθ(θᵀGy)/(θᵀGθ). Then ‖P‖²_{G⁻¹} = qᵀGq = ‖y‖²_G − 2h(θᵀGy)²/(θᵀGθ) + h²(θᵀGy)²/(θᵀGθ). So the constant h(2−h) is correct on the active set: sympy residual 0, numerical error 2e-15 on 14899 states.
(a) Take θ = (2,0,0) on the boundary, so h = 1, with G = I and y = (−1,0,0). Then Proj = (−1,0,0) and ‖Proj‖² = 1, but the display gives 0. On 15011 random states with h > 0 and θᵀGy < 0 the error reaches 100% of ‖y‖²_G.
(b) Take n = 1, G = I, a_max = 2, ϑ = (2,0,1), y = (1,0,0). TransportController.project returns 0. The whole-ϑ form would need h(2−h) = ‖y‖²·ϑᵀϑ/(ϑᵀy)² = 1.25, but h(2−h) ≤ 1 for every h.
(c) For s in (s_min, s_min+ε_K) the box returns y·(s−s_min)/ε_K, and only when y < 0.
The block-sum identity in the fix matches TransportController.project to 6.4e-16 on 30480 active particle-states, and the largest ratio ‖Proj_K‖²_{G⁻¹}/‖y‖²_G is 1.000000.

**Fix.** Insert before Theorem F:

'Lemma P.1 (Proj_K is non-expansive). Let G = diag(g_aI_n, g_cI_n, g_s). For every ϑ = (a,c,s) and y = (y_a,y_c,y_s),
‖Proj_K(ϑ,Gy)‖²_{G⁻¹} = ‖y‖²_G − Σ_{b∈{a,c,s}} h_b(2−h_b)[(n_bᵀG_by_b)₊]²/(n_bᵀG_bn_b) ≤ ‖y‖²_G.
Here n_a = a and n_c = c, with h_a = h(‖a‖; a_max, ε_K) and h_c = h(‖c‖; c_max, ε_K) the quadratic ramps.
For s < s_min+ε_K, n_s = −1 and h_s = 1 − clip((s−s_min)/ε_K, 0, 1).
For s > s_max−ε_K, n_s = +1 and h_s = 1 − clip((s_max−s)/ε_K, 0, 1).
Otherwise h_s = 0.
Proof. G is block-diagonal, so both squared norms split over the blocks. On an active ball block, Proj_b = G_bq_b with q_b = y_b − h_bn_b(n_bᵀG_by_b)/(n_bᵀG_bn_b), and expanding q_bᵀG_bq_b gives the subtracted term. On an inactive block, Proj_b = G_by_b. The s-block is the one-dimensional case. Since h_b ∈ [0,1], h_b(2−h_b) ∈ [0,1]. ∎'

Then:
- In the proof of C(i), replace 'because the projection is a G⁻¹-contraction: ‖Proj(ϑ,Gy)‖²_{G⁻¹} = …' with 'by Lemma P.1'.
- In the proof of F, replace '(i) uses the projection's G⁻¹-contraction property (Theorem C(i)) and Lemma B.2' with '(i) By Lemma P.1, ‖ϑ̇_i‖_{G⁻¹} ≤ ‖∇g_t(ϑ_i)‖_G ≤ λ_G^{1/2}‖∇g_t(ϑ_i)‖ ≤ λ_G^{1/2}L₁‖e(t)‖ (the last step is Lemma B.2(b)); square and average over i.'
- Log the erratum in review/REVIEW_LOG.md. REVIEW_round1.md l.525 uses the joint form, but it is a raw report; leave it verbatim.

### F.2 [minor] PARADIGM.md l.447–448 (*Reading* of Theorem F), placed directly after the (L_σ) clause at l.443

Merges R1-2 and the first half of R2-2; CONFIRMED. The Reading says that moving the kernel measure by W costs at least W/(λ_G^{1/2}L₁) of integrated tracking error. That follows from (ii) only along (L). The Reading carries no qualifier and sits right after the (L_σ) clause. Under any leak, the leak term can pay for the whole motion, so F gives no lower bound on ∫‖e‖.

**Argument.** Run with the project's code (TransportController(leak='ot'), RK4, dt = 1e-3):
- n = 1, m = 2, test box and gains, σ_m = 0.5.
- Ideal: two copies of (1,0,1). Prior: a⁰ = (1.5, 0.5). Start: a = (3, −1).
- All particles have (c,s) = (0,1), and x(0) = x_m(0).
The leak preserves a₁ + a₂, so f̂ ≡ f and e ≡ 0.
Result: sup‖e‖ = 7.1e-14 and ∫₀¹⁰‖e‖ = 4.1e-13, while W(μ_10,μ_0) = 1.0535. The Reading would require ∫‖e‖ ≥ 0.0293.
The integrated (L_σ) clause still holds, with equality (ratio 1.0000), because the motion is two straight lines of equal length.

**Fix.** Replace l.447–448 with:
'*Reading.* Along (L), moving the kernel measure by W costs at least W/(λ_G^{1/2}L₁) of integrated tracking error. This is an exchange rate between learning and tracking. Along a leaky law the proved rate is ∫_{t₀}^{t₁}(λ_G^{1/2}L₁‖e‖ + σ_m d_r)dr ≥ W(μ_{t₁},μ_{t₀}), where d_r = W(μ_r,μ⁰) for (L_σ) and d_r = d_lab(r) for the labelled leak. This gives no lower bound on ∫‖e‖: the leak can pay for the whole motion. Under (L_σ) there are runs with e ≡ 0 and W(μ_10,μ_0) = 1.05.'
Move 'Checked at every sample' to the verification map.

### F.3 [minor] PARADIGM.md l.356–358 (Proposition C.5, sentence 'the transient costs tracking error, as Theorem F prices it'); F's leak clause at l.443

Merges R1-3 and the second half of R2-2; CONFIRMED with a numeric correction. C.5 applies F to the labelled leak, ϑ̇_i = Proj_K(ϑ_i, G∇g − σ_m(ϑ_i − ϑ⁰_i)). That law is neither (L) nor (L_σ). F's (L_σ) clause, with W(μ_r,μ⁰), is false for it; the correct added term is σ_m·d_lab(r), where d_lab := ((1/m)Σ‖ϑ_i−ϑ⁰_i‖²_{G⁻¹})^{1/2} ≥ W(μ_r,μ⁰). Even with that correction, F under a leak only bounds motion from above by tracking error plus the leak term. It cannot 'price' a transient in tracking error.

**Argument.** (1) The δ = 0 instance of C.5:
- n = 1, G = I, ideal = prior = {(0,0,1),(1,0,1)}, start = the swap, σ_m = 1, leak = 'labelled'.
- max|e| = 1.5e-13, but the particles travel a path of length 1.000 and max_t W(μ_t,μ_0) = 0.498.
- F's (L_σ) clause applied to this law fails by 39.7× on sampled pairs, and the ratio is unbounded as t₁ → t₀ = 0, because W(μ_0,μ⁰) = 0.

(2) The experiment itself (exp_sigma_mod_transient.py, seed 3, its real initial state, [0,15]):
- The integrated clause with W(μ_r,μ⁰) is violated by up to 1.60×, and the pointwise clause by up to 1.75×.
- Both hold with d_lab (0.77).
- The OT run satisfies the clause (0.11).

(3) Correction to R1: the '127×' figure comes from a hand-set ‖e‖ = 1.4e-4 (x = x_ref + 1e-4·1). The experiment actually starts at ‖e(0)‖ = 0.353, where the t = 0 ratio is about 0.05.

(4) transient.md:
- The labelled path is 5.8× longer (2.26 vs 0.39), but RMS‖e‖ on [0,5] is only 4% higher (0.0771 vs 0.0743), and equal from t = 10 on.
- For seed 3 on [0,15], λ_G^{1/2}L₁∫‖e‖ = 9.1, far above the path length 2.16, so the (L) bound is not binding either.

**Fix.** (a) In F, after l.443 add: 'For the labelled leak ϑ̇_i = Proj_K(ϑ_i, G∇g − σ_m(ϑ_i − ϑ⁰_i)), add instead σ_m d_lab(r), where d_lab(r) := ((1/m)Σ‖ϑ_i(r) − ϑ⁰_i‖²_{G⁻¹})^{1/2} ≥ W(μ_r,μ⁰). The term W(μ_r,μ⁰) does not suffice there.'

(b) In C.5, replace 'The steady states coincide, but the transient costs tracking error, as Theorem F prices it (`results/ot_sigma_mod/transient.md`: particle path 5.8× longer).' with: 'The steady states coincide, but the labelled leak forces a reshuffle. In `results/ot_sigma_mod/transient.md` the particle path is 5.8× longer (2.26 vs 0.39), while RMS‖e‖ on [0,5] is only 4% higher (0.0771 vs 0.0743). Theorem F does not price this motion: with a leak it bounds transport from above by tracking error plus the leak term, and the leak can pay for all of it. With δ = 0 above, e ≡ 0 while every particle crosses the domain.'

### F.4 [cosmetic] tests/test_transport.py l.307–321 (test_theorem_F_wasserstein_length_bounded_by_integrated_error); PARADIGM.md l.448 ('Checked at every sample'), l.742 (verification-map row F)

Merges R1-7 and R2-3; CONFIRMED. R2's severity 'minor' is lowered to cosmetic, because nothing stated is false and F is proved.
The only check is that particle path length ≤ λ_G^{1/2}L₁∫‖e‖ along (L), and it passes with about a 29× margin, so it cannot detect a wrong constant. The projection step is never tested on its own, and in seed 0 no factor is ever active. The (L_σ) clause is never checked. The verification map calls the checked quantity a 'Wasserstein path length'.

**Argument.** Replicated:
- Max path-length/bound is 0.0342 (seed 0) and 0.0290 (seed 1).
- Max pointwise ratio is 0.034 and 0.032.
- The projection is active at 0/300 samples in seed 0 and 77/300 in seed 1.
- The existing test passes (24.9 s).
The pointwise bound is nearly attainable elsewhere: 0.98 over random states, and 0.9989 analytically for m = n = 1, a_max = 0.05, s_min = 1, g_a = λ_G, x = c.

**Fix.** At every sample, assert:
- N(ctrl.velocity) ≤ N(G∇g_t)(1+1e-12), with equality when no factor is active. This tests Lemma P.1.
- N(G∇g_t) ≤ λ_G^{1/2}L₁‖e‖.
Add a tight instance (m = n = 1, a_max = 0.05, s_min = 1, g_a = λ_G, x at the centre) asserting ratio ≥ 0.99.
Add an (L_σ, OT) run asserting W(μ_{t₁},μ_{t₀}) ≤ ∫(λ_G^{1/2}L₁‖e‖ + σ_mW(μ_r,μ⁰))dr on a grid of pairs.
Rename the map entry 'particle path length (≥ Wasserstein length)'.

### F.5 [cosmetic] PARADIGM.md l.71–73 (§1: 'It is Lipschitz, C¹ at the inner edge of the layer'); src/core/projection.py l.33 and l.156 ('smooth' box)

Merges R1-8 and R2-4; CONFIRMED. BandwidthProjection multiplies by the linear ramp clip((s−s_min)/ε_K, 0, 1), and symmetrically at the upper face. So the box has a kink at s = s_min+ε_K. Only the norm-ball factors, with their quadratic ramp, are C¹ at the inner edge. The proofs use only Lipschitz continuity and ramp values in [0,1], so nothing downstream breaks.

**Argument.** BandwidthProjection(0.3, 3.0, 0.15), y = −1, at s = 0.45: the left derivative is −6.6667 (= −1/ε_K) and the right derivative is 0. The ball operator (θ_max = 2, ε = 0.3, y = +1) at r = 1.7 has one-sided derivatives 0 and −1.1e-5, so it is C¹ there.

**Fix.** Replace the §1 sentence with: 'It is Lipschitz. The norm-ball factors are C¹ at the inner edge of their layers; the box factor (linear ramp), and every factor at the outer edge, is only Lipschitz.'
Alternatively, use r(s) = 1 − (1 − (s−s_min)/ε_K)² on the lower layer, and symmetrically at the upper. It stays in [0,1], so (P) and Lemma P.1 still hold.
Drop 'smooth' from the box docstrings.

### F.6 [cosmetic] PARADIGM.md l.437 (statement of F), l.444 (proof of F), l.200–201 (proof of Lemma B.2)

Merges R1-4, R1-5 and R2-5 item 4; CONFIRMED.
- F states no hypotheses and no time interval. Its proof needs ϑ_i(t) ∈ Θ_K, from (S3) plus forward invariance (Nagumo), so that ‖a_i‖ ≤ a_max and s_i ≥ s_min feed L₁. It also needs existence on [t₀,t₁], from (S1)–(S2) on the maximal interval. Calling (S) 'standing' covers this only implicitly.
- The bound ‖∇g_t‖ ≤ L₁‖e‖ appears only as an unproved clause inside the proof of Lemma B.2.
- The step ‖G∇g‖_{G⁻¹} = ‖∇g‖_G ≤ λ_G^{1/2}‖∇g‖ is not written.
- The projection step is cited as 'Theorem C(i)', but it appears only in C(i)'s proof.

**Argument.** The two-line derivation holds: ∇_ag = φe and ∇_{(c,s)}g = (aᵀe)∇φ, with φ ≤ 1 and ‖∇φ‖ ≤ J on s ≥ s_min. Numerically, max ‖∇g_i‖/(L₁‖e‖) = 0.9987 over 6000 adversarial states. Without (S3), s_i(0) < s_min allows |∂_sφ| > J_s, and (i) can fail at t = 0.

**Fix.** Begin F with: 'Assume (S1)–(S3). On the maximal interval of existence (all t ≥ 0 under the hypotheses of Theorem D(b)), for 0 ≤ t₀ ≤ t₁:'. Add that μ*, (S4), (S5) and X are not needed.
Add 'Lemma B.2(b): on Θ_K, ‖∇_ϑg_t(ϑ)‖ ≤ L₁‖e(t)‖', with the two-line proof above.
Cite Lemma P.1, not Theorem C(i).

### F.7 [cosmetic] PARADIGM.md l.443 ('Along (L_σ), add σ_mW(μ_r,μ⁰) inside the integral') and the proof of (ii), l.444–446

Merges R1-9 and R2-5 item 3; CONFIRMED. The (L_σ) clause has three gaps, each true and easy to close:
- It is meaningless for (i) and (iii), which contain no integral.
- It does not address Filippov solutions at assignment ties.
- AGS Thm 8.3.1 is Euclidean, while W uses the G⁻¹ cost.

**Argument.** The pointwise (i′) for the OT leak holds, with max ratio 0.9992 over 6000 random states. It is tight: the pure-leak run in issue 2 gives ratio 1.0000. The Filippov argument (convexity of N) and the isometry ϑ ↦ G^{−1/2}ϑ are standard.

**Fix.** State, for a.e. t:
- (i′) N(ϑ̇) ≤ λ_G^{1/2}L₁‖e‖ + σ_mW(μ_t,μ⁰);
- (ii′) W(μ_{t₁},μ_{t₀}) ≤ ∫_{t₀}^{t₁}(λ_G^{1/2}L₁‖e‖ + σ_mW(μ_r,μ⁰))dr;
- (iii′) W(μ_{T₂},μ_{T₁}) ≤ λ_G^{1/2}L₁(T₂−T₁)RMS‖e‖ + σ_m∫_{T₁}^{T₂}W(μ_r,μ⁰)dr.
Add one sentence: 'At ties, Filippov velocities are convex combinations of fields with τ optimal at the current point, each satisfying (i′), and N is convex.'
Add another: 'For general measures, use the Lagrangian coupling (X_{t₀},X_{t₁})#μ_0 with Minkowski's inequality, or AGS 8.3.1 after the isometry ϑ ↦ G^{−1/2}ϑ.'

### F.8 [cosmetic] PARADIGM.md l.437–442 (time variable s, RMS), l.42 (§0 'Wasserstein length'), l.320 and l.444 ('G⁻¹-contraction')

Merges R1-10 and R2-5 items 1, 2 and 5; CONFIRMED.
- The letter s denotes both time ('0 ≤ s ≤ t', '∫_s^t') and the width coordinate of ϑ.
- RMS_{[T₁,T₂]} is used at l.227, 247 and 441 but never defined in PARADIGM.md, and (iii) needs T₂ > T₁.
- §0 claims a bound on the 'Wasserstein length', but F states only the endpoint displacement.
- 'G⁻¹-contraction' names what is a norm bound (non-expansiveness).

**Argument.** A grep finds no definition of RMS in PARADIGM.md. Summing (ii) over any partition gives Σ_kW(μ_{t_k},μ_{t_{k+1}}) ≤ λ_G^{1/2}L₁∫‖e‖. That is the missing one-line length statement.

**Fix.** - Use t₀ ≤ t₁ for times and r as the integration variable.
- Define RMS_{[T₁,T₂]}‖e‖ := ((T₂−T₁)⁻¹∫_{T₁}^{T₂}‖e‖²)^{1/2} for T₂ > T₁ at its first use (Corollary B.3).
- After (ii) add: 'Hence the W-length of r ↦ μ_r on [t₀,t₁] is ≤ λ_G^{1/2}L₁∫_{t₀}^{t₁}‖e‖.'
- Replace 'G⁻¹-contraction' with 'non-expansive (‖Proj_K(ϑ,Gy)‖_{G⁻¹} ≤ ‖y‖_G)'.

### F.9 [cosmetic] PARADIGM.md l.437–442 (constant λ_G^{1/2}L₁ in (i)–(iii))

Merges R1-6 and R2's optional remark; CONFIRMED, optional. The constant is valid but not sharp for the test box. The same proof gives the smaller block constant L_G and a stronger W_∞ bound.

**Argument.** For the test box:
- λ_G^{1/2}L₁ = 35.988 and L_G = 21.356.
- The sharp supremum of ‖∇g‖_G/‖e‖ is 20.44. The Euclidean supremum of ‖∇g‖/‖e‖ is 24.44, against L₁ = 25.45.
- Over random states, max ‖∇g_i‖_G/(L_G‖e‖) = 0.9991.
L_G ≤ λ_G^{1/2}L₁ always. The two nearly coincide when g_a = λ_G and a_max is small (0.9989), so the constant is not uniformly loose.

**Fix.** Optionally add: '(i)–(iii) also hold with L_G := (g_a + a_max²(g_cJ_c² + g_sJ_s²))^{1/2} ≤ λ_G^{1/2}L₁. Since max_i‖ϑ̇_i‖_{G⁻¹} ≤ L_G‖e‖, also W_∞(μ_{t₁},μ_{t₀}) ≤ L_G∫_{t₀}^{t₁}‖e‖.'

Refuter verdicts: holds_with_fixes; holds_with_fixes
