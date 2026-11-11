<!-- Archived verbatim from the referee's final message (round IS2, 2026-11-07). The coordinator added only this comment. -->

# REVIEW round IS2: adversarial referee report on Theorem "separated regime" (supp_sep.tex) and its proof outline

**Scope and provenance.** I read the current text of `supp_sep.tex` (edited 21:01:03; I re-read it after the edit), `supp_theory.tex`, `supp_theory2.tex`, `sec_theory.tex` (21:01:22), `abstract.tex`, `sec_intro.tex`, `sec_conclusion.tex` and `tab_sep.tex`. From the repository I read `src/core/rich_transit_nd.py`, `experiments/exp_separated_regime.py` and `results/chaos/separated_regime.md` with its `_raw.npz` (regenerated 20:59). Nothing in the repository or the paper directory was modified. My scripts and logs are in `/tmp/claude-0/-home-user/aa355708-d4bf-50e9-ad81-220b64c70ccb/scratchpad/referee_IS2/` (list at the end).

Notation: Γ(d)=exp(−d²/4σ²), h = lattice spacing, d_1 = distance from a feature to its best kernel, m = min(μ_0, g).

## 0. Bottom line

1. The architecture of the proof is right: perturb the decoupled path P_0, with only μ_0, g and block nondegeneracy mattering, and with r_0 independent of Δ. I found no counterexample to the supplement's theorem in about 1100 runs. They cover 1-D and 2-D, amplitude ratio R up to 10, lattice spacing 0.83σ and 2.5σ, and large and small margins.
2. One clause of the main-text theorem, "T_k(1+O(Γ(Δ)))", is false as worded (finding 1). The repair is a wording change to the supplement's own δ_Δ.
3. The theorem as proved covers none of the experiment. On the proof's own displayed constants, Δ_0 is at least 8.2σ (lower estimate) and about 12–14σ (with δ_Δ). The tested separations are at most 8σ for K=2 and at most 5.8σ for K=3, in a window 11.7σ wide. 0 of 720 targets are covered, and no K=3 target in that window could ever be covered (finding 2).
4. Three pieces are not established as written: the entry-layer convergence (H3, finding 3), (H4) for all levels (finding 4) and the mobility-weighted nondegeneracy in (H2) (finding 5). "(H1–H4) hold" and "proved" are therefore over-claims today (finding 8 and §4).

## 1. Verified (re-derived or recomputed; these agree with the text)

- **V1. Hessian.** At the one-feature minimum the Hessian is diag(Z, A²Z/(2σ²)·I_n).
  - For general n it follows from ⟨∂_aφ ∂_bφ⟩ = δ_ab Z/(2σ²) and ⟨φ ∂_cφ⟩ = 0.
  - Finite differences in 1-D (A=0.8, σ=0.12) give diag(0.026587, 0.5908) against the prediction (Z, A²Z/2σ²) = (0.026587, 0.5908), off-diagonal 1.6e−15.
  - Cross-feature Gram entries are O((1+Δ²/σ²)Γ(Δ)).
- **V2. Step 1 algebra.** F_S = ½⟨(f_S−Σw_iφ_i)²⟩ − Σw_i⟨φ_i f_U⟩ + const is correct.
- **V3. Lemma "monotone entrant".** I re-derived ḋ = −κ w A d Γ/(2σ²), dΓ(d(s))/ds = κ w A d²Γ²/(4σ⁴) ≥ 0 and ġ = Λ_0 G m Z² A² w d² Γ²/(4σ⁴).
  - The first-zero argument is valid: g>0 before s_0 forces ġ(s_0) ≤ 0, but ġ(s_0) > 0 unless d=0.
  - d=0 is invariant, so d(0) ≠ 0 implies d(s) ≠ 0 for all finite s.
  - It also gives w < AΓ(d) ≤ A.
  - Numerics: I integrated the ε-free inner flow of block plus entrant (all of `inner.py` and `inner2.py`) in 24 configurations: Δ = 2.5–5.4σ, A_2/A_1 ∈ {0.6, 0.1}, both signs, entrant started up to 0.035–0.045 off its feature.
    - The flow converged to the exact target (rms 1e−16 to 1e−9).
    - The entrant coefficient w_2(s) was monotone (0 decreasing steps above 1e−9 relative).
- **V4. Step 2 classes.**
  - The three classes exhaust the lattice, and r_0 is independent of Δ.
  - r_0σ = 2σ√ln(4KR/(μ_0 Γ(h√n/2))) is explicit.
  - The covering radius h√n/2 is right for hZ^n.
  - The bounds are per kernel, so no count of kernels within r_0σ is needed. For information, the count is ≈ (2r_0σ/h+1)^n: about 12 in the 1-D benchmark with r_0 ≈ 4.7σ, and about 100 for the 2-D benchmark.
- **V5. Step 3 competitors.** Only these enter: the runner-up of the same feature (μ_0), the next feature's best kernel (g), kernels frozen near fitted features (they need time ≥ μ_0/(Cδ) ≫ T_max), and the far class. Nothing else competes.
- **V6. Step 5.** Coefficients stay near A_l, so no unit exits.
- **V7. Step 6 functional form.** Δ_0 ~ √ln(1/m) is the right form. The data below give exactly this form with offset h/2 instead of (r_0+2)σ.
- **V8. Experiment recomputed from the raw file.** All of the following agree with the md: exact 720/720, 302/302 above 5Γ, per-row counts, and slope 1.01 for K=2.
  - The K=3 slope between 4σ and 5σ is 0.73. The script no longer prints it.
- **V9. First-entry time error.** For fixed local data (strong feature A=1, best kernel displaced d_1=0.04 toward its neighbour, neighbour A_2=0.6) the first-entry relative error equals (A_2/A_1)Γ(Δ−d_1)/Γ(d_1).
  - Observed/predicted = 0.998, 1.000, 1.000, 1.000, 1.000, 1.000 at Δ = 5, 5.8, 6.7, 8.3, 9.2, 10σ (`f1.py`).
- **V10. Dormant-kernel margins.** The path's actual margin of never-entering kernels is ≥ μ_0/2 on all targets at ≥ 6σ.
  - The ratio min(actual margin)/μ_0 is 0.74 at 6σ and 1.000 at 8σ.
  - At 5σ it is ≥ μ_0/2 on 29/30 targets, and at 4σ on 24/30 (`margin.py`).
- **V11. 2-D behaves like 1-D.** I ran 120 targets: K=2, σ=0.15, h=0.125 (0.83σ), 289 kernels, separations 3–6σ.
  - Exact 120/120.
  - Same entry set: 10, 25, 29, 30 of 30.
  - All failures have m ≤ 1.2Γ.
  - Entry-time error ≈ 1.1–1.5Γ.

## 2. Findings

### F1. [BLOCKER as worded; one-line repair] The main-text time accuracy O(Γ(Δ)) is false for fixed local data
- **Where.** `sec_theory.tex` Thm "separated regime" ("at the times T_k(1+O(Γ(Δ)))"). Also the md/CLAIMS statements "entry-time error ~ Γ(Δ)", and `supp_sep.tex` l. 38 "(theory: O(Γ))". The supplement's theorem (i) says O(δ_Δ), so the two disagree.
- **Evidence.** With the local data fixed and Δ moved by lattice translations, the first entrant's dual rate is a_{jk}+a_{jl}.
  - Its relative time error is (A_l/A_k)·Γ(Δ−d_1)/Γ(d_1), where d_1 is the best kernel's offset toward the neighbour.
  - That is Γ(Δ)·exp((2Δd_1−d_1²)/4σ²)·const, which is not O(Γ(Δ)) uniformly in Δ for any d_1>0.
  - The deterministic run (V9) gives error/Γ(Δ) = 1.38, 1.59, 1.82, 2.41, 2.77, 3.18 at Δ = 5, 5.8, 6.7, 8.3, 9.2, 10σ. It grows without bound while matching the formula to 3–4 digits.
  - In the report's own K=2 table, the maximum error divided by Γ(Δ−h/2) is 0.96, 0.99, 0.99, 0.89 at Δ = 5, 6, 7, 8σ. Divided by Γ(Δ) it is 2.60, 3.32, 4.08, 4.51.
  - The "slope 1.01" is the slope of medians over random offsets on 4–8σ.
- **Fix.** Replace by O(δ_Δ), or prove and state O(Γ(Δ−h√n/2)·poly). Report the maximum error against Γ(Δ−h/2), not the median against Γ(Δ).

### F2. [MAJOR] The proven regime lies above every experiment, and the experiment tests a different, sharper statement
- **Where.** `supp_sep.tex` Step 6, the Remark and l. 38; `exp_separated_regime.py` (MIN_C=5, "the theorem needs both to exceed CΓ"); `sec_theory.tex` record sentence.
- **Evidence.** Minimal explicit constants read off the displayed conditions (`d0est.py`), for the benchmark (σ=0.12, h=0.1, n=1, so Γ(h/2)=0.9575):
  - (a) Lower estimate: the cross-talk of a near kernel, K R Γ(Δ−r_0σ)/Γ(h/2), must stay below m/4, the same constant as in the definition of r_0. This gives Δ_0 ≥ r_0σ + 2σ√ln(4KR/(m Γ(h/2))).
  - (b) Face value: K R (1+Δ²/σ²) Γ(Δ−(r_0+2)σ)/Γ(h/2) ≤ m/4. Note that δ_Δ ≥ 1 for all Δ ≤ 11σ (δ_Δ ≈ 48 at 8σ, ≈ 0.2 at 12σ).
  - Over the 720 targets (R median 1.3, max 2.4; μ_0 median 0.047; g median 0.15; r_0 = 4.0–6.8σ, median 4.7σ):

    | Estimate | Median | Min | Max |
    |---|---|---|---|
    | (a) | 9.6σ | 8.2σ | 13.6σ |
    | (b) | 13.4σ | 12.2σ | 17.1σ |

  - The margin cannot be large on this lattice: μ_k ≤ 1−Γ(h)=0.159 always. Best case (K=2, R=1, μ_0=g=0.159) gives (a) 8.0σ and (b) 12.0σ. Typical (R=2, μ_0=g=0.05) gives 9.6σ and 13.5σ.
  - The window for positions is [−0.7, 0.7] = 11.7σ. K=2 separations tested are ≤ 8σ; K=3 minimal separations are ≤ 5.8σ.
  - Result: 0 of 720 targets satisfy Δ ≥ Δ_0 under either estimate, and no K=3 target in the window can, whatever its margins.
  - Observed: exactness from 2σ; same entry set 60/60 from 6σ (K=2), 59/60 at 5σ, 60/60 for K=3 at 5σ.
  - The 5Γ gate is not the theorem's condition (δ_Δ with the (r_0+2)σ shift and polynomial factor).
  - It is also lattice-specific. On a sparse lattice (h=0.3=2.5σ, 7 kernels, K=2, 40 targets per row at 3, 4, 5, 6, 8σ; `sparse.py`) margins are large (μ_0 median 0.29–0.54, max 0.79). The rows are exact 200/200, but 3 of the 153 targets above 5Γ depart from the decoupled set: m = 0.13 and 0.10 at 4σ, and m ≈ 0.01 at 5σ, all with d_1 ≈ 1.2σ.
  - 245 of the 302 targets passing the 5Γ gate sit in rows at Δ ≥ 5σ where almost everything agrees anyway. The informative rows are 3.5–4σ.
  - A parameter-free gate does better. Take min(μ_0,g) > Γ(Δ−h√n/2), the first-order cross-talk of a kernel half a cell from its feature:
    - 720 targets: 356/356. κ=0.75 gives 374/377.
    - Independent random targets with R up to 8.9 (K=2 at 4σ and 5σ, K=3 at ≥ 4σ): 88/88.
    - Sparse lattice: 139/139.
  - All same-set failures have m ≤ 2.0Γ(Δ) (720 targets) or ≤ 2.3Γ(Δ) (my extra 140 random targets).
  - The observed threshold is Δ*(m) ≈ h/2 + 2σ√ln(1/m). That is the functional form of Step 6 without the ≈ 4–7σ offset (r_0+2)σ.
- **Fix.**
  1. State in the text that the experiment is not a test of the theorem's Δ_0. The 5Γ statement "no target with margins above 5Γ departs" is true only on this lattice.
  2. Replace the gate by Γ(Δ−h√n/2) and present it as a conjectured sharp form.
  3. In the proof, use two radii. Units that enter are within h√n/2+2σ of their feature, which gives a small shift in δ_Δ. The margin conditions need the sup over kernels of the own rate plus cross-talk of the dual at T_max, not r_0.
     - That sup is attained at the runner-up (R Γ(Δ−d_2)/Γ(h/2) ≲ μ_0) or at the midpoint kernel ((1+R)Γ(Δ/2)/Γ(h/2) ≲ ½, i.e. Δ ≳ 4σ√ln(2(1+R)/Γ(h/2))).
     - Heuristic value for R=2.5, μ_0=0.05: about 5.6σ, which agrees with the observed onset at 5–6σ.
     - This is my suggested route, not a result.

### F3. [MAJOR] The coupled entry (H3), Step 4 via Lemma 3(iii), is a sketch whose hypotheses do not match what Step 4 supplies
- **Where.** `supp_theory.tex` l. 36 and its proof (iii); `supp_sep.tex` Step 4.
- **Evidence.**
  - (iii) assumes "the block sits at a nondegenerate minimum of its own fit". Step 4 supplies the critical point of F_S under cross-talk, which is O(δ) away. The proof's F_× omits the term ⟨(f_B*−f̂_B) f_⋆*⟩ (the block-dependent overlap with the entrant's feature). It is not multiplied by w_⋆, so it shifts the block by O(δ) for all s.
  - "for small δ" is never quantified. The closing step ("F is a Lyapunov function bounded away from {w_⋆=0}, IFT closes") is not carried out.
    - A finite-horizon Gronwall continuation would make δ_0 ~ exp(−L·S_1), with S_1 the slide length. The paper's own Remark rem:regime says the analogous near-field estimate forces ε_0 exponentially small.
    - The slide length depends on d_0/σ ≤ h√n/(2σ), on the amplitudes and on the mobilities. None of these is in the list of parameters of Δ_0. The relaxation rate drops with amplitude: in my inner flows the final rms after the same horizon is 7.6e−11 at A_2=0.1 against 1e−16 to 1e−15 at A_2=0.6.
    - A compactness route avoids the exponential but gives a non-explicit δ_0. The margin-independent threshold then makes Δ_0 = max(margin term, unknown constant).
  - "converges within O(δ_Δ) of the decoupled slide" is true in state space, not in inner time. Crossing times s^θ shift by O(δ s^θ).
  - Numerics support (H3): the 24 configurations in V3 converge, and the 720 paths end at the exact target.
- **Fix.**
  - State Lemma 3(iii) as: there is δ_0(K, n, h/σ, A_min, A_max, mobilities) such that for δ ≤ δ_0 …; add these to Δ_0's dependencies.
  - Carry out the Lyapunov classification of critical points in the sublevel set, including the manifold {w_⋆=0}.
  - Correct the hypothesis to "critical point of the block's fit to the whole target".
  - Call (H3) proved only for the pair "small δ_0, non-explicit".

### F4. [MAJOR] (H4) "once and transversally for every level θ" is not what Lemma (mono) gives
- **Where.** `supp_sep.tex` Step 4, last sentence ("…at fixed θ (H4)"), and Thm (iv).
- **Evidence.**
  - Monotonicity is proved for the decoupled problem only. The cross-talk breaks the first-zero argument near the end, where the decoupled term κ w A d²Γ² is smaller than the O(δ) perturbation.
  - C¹-closeness on compacts gives the crossing property only for fixed θ, with a threshold depending on θ. Levels near 0 (far field) and near ŵ_k are not controlled uniformly. The theorem asserts a θ-independent Δ_0.
  - Theorem 1(c) needs fixed θ only, so the repair is easy. Numerically no non-monotonicity was seen (V3).
- **Fix.** Either state "(H4) holds at each fixed level θ∈(0,|ŵ_k|), with Δ_0(θ)". Or prove monotonicity of the perturbed entrant: an exponentially weighted estimate of ẇ^δ−ẇ^0 using the non-vanishing leading mode. Overshoot beyond ŵ_k does not violate (H4).

### F5. [MAJOR] (H2) is the mobility-weighted Hessian; only the plain Hessian is bounded, and λ_0 and the smallness conditions are not dimensionally consistent
- **Where.** `supp_theory.tex` (H2) (λ_k is the smallest eigenvalue of M_k^{1/2}H_kM_k^{1/2}); `supp_sep.tex` Thm statement, Steps 1 and 6.
- **Evidence.**
  - The new proof bounds H_k ≥ λ_0/2 and never M_k. The coefficient mobility |𝒲'(ξ)|/ρ depends on the coefficient value; the centre mobility is Λ_0 G m.
  - The list of dependencies of Δ_0 omits G, m, v_s, a_max, b and the amplitude scale.
  - λ_0 = Z·min(1, A²/(2σ²)) mixes units: the (w,w) and (c,c) Hessian entries have different units, so the smallest eigenvalue of a mixed-unit matrix is not meaningful.
  - "Cδ_Δ ≤ cλ_0²σ" and "Cδ_Δ ≤ c|A|_min" compare a dimensionless δ_Δ with dimensional quantities. A_min appears as `A_min` in λ_0 and as `|A|_min` elsewhere.
  - In the benchmark's units A²/(2σ²) ≈ 22–35 ≫ 1, so λ_0 = Z for every amplitude, while the layer's relaxation keeps slowing with A (F3).
  - The Gauss–Newton Hessian is also not the full Hessian of F_S: the residual-curvature term is O(δ) and is to be listed.
- **Fix.** Scale the coordinates (w/A_l, c/σ) and normalise (Z=1, A_max=1). Define λ_k as in (H2) and give λ_k ≥ m_min λ_0/2 with m_min = min(mobility of the coefficient at |w|=A_min/2, Λ_0Gm). List the model constants in Δ_0's dependencies.

### F6. [MINOR, but false as written] Step 1: uniqueness "in the ball |θ−θ*_S| ≤ σ"
- **Where.** `supp_sep.tex` Step 1.
- **Evidence.**
  - "Hessian ≥ λ_0/2" is only at θ*. Newton–Kantorovich gives uniqueness in a ball of radius ~ λ_0/L_H, not σ.
  - On the ball of radius σ the Hessian of the block's own fit fails the λ_0/2 bound. For A=0.8, σ=0.12, Hessian of the one-feature fit (`hess.py`), λ_0/2 = 0.0133:

    | Displacement d | w | Min eigenvalue |
    |---|---|---|
    | 0.7σ | A−0.4 | 0.0105 |
    | σ | A | 0.0054 (0.2λ_0) |
    | σ | A−0.4 | −0.0112 |
    | 1.25σ | A | −0.0239 |

    It holds up to about 0.5σ.
  - The ball also mixes units (w and c in one norm).
  - NK gives existence and local uniqueness. Reaching the point by the flow comes only from Lemma 3(iii) (F3).
  - The C² bound by Cδ_Δ itself is fine. In fact |⟨φ f_U⟩| ≤ Z K|A|_max Γ(Δ−σ) there.
- **Fix.** State uniqueness in B(θ*, r_*) with r_* ≈ 0.5(σ, A) or the NK radius, and say that the dynamics come from Lemma 3(iii).

### F7. [MAJOR] Hidden uniformities and quantifiers
- **(a) Local data and random positions.**
  - The theorem fixes local data and sends Δ to ∞ along lattice translations. That is coherent (the lattice must be infinite, or a finite one containing the windows), but the experiment draws positions at random, so every target has its own Δ_0.
  - For random positions the theorem implies only a statement on {min(μ_0,g) ≥ m(Δ)}, with Δ_0(m(Δ)) = Δ. This is a probability bound ≈ 1 − O(m(Δ)/0.1) only after explicit constants. With the proof's constants m(Δ) = 1 for Δ ≲ 12σ, so nothing follows at Δ ≤ 12σ.
  - Typical margins are small because μ_k ≈ h·x/(2σ²) ≤ h²/(4σ²). As the dictionary refines, Δ_0 grows like σ√ln(σ/h).
- **(b) h/σ enters exponentially.**
  - C ∝ 1/Γ(h√n/2) = exp(h²n/16σ²), so the constants inside Δ_0 grow like h√n.
  - Genericity (μ_0>0, g>0) does not exclude two features sharing the same best kernel when Δ < h. The proof's own conditions force Δ_0 ≥ h√n, but the statement never says Δ>h√n (needed for "one kernel per feature", distinct j_k).
- **(c) Equal walls "change the notation only" is false.**
  - With unequal walls the best kernel minimises ρ_j/|a_jk|, not distance, and the constants carry ρ_max/ρ_min.
  - Main-text Prop "one feature" says "the unit of largest correlation enters", which is right only for equal walls.
- **(d) Finite dictionary.**
  - Theorem 1 is for a finite dictionary; the corollary needs a finite sub-lattice containing the K best kernels.
  - ε_0 then depends on the number of dormant kernels (N_ker·ε^{3/2}μ^{-3} in the residual) and is non-uniform in the layer length, per the paper's own Remark.
  - The corollary holds for fixed target as ε→0 only.
- **(e) Scale non-invariance.** The limit object's entry structure is invariant under f*→κf* (dual times scale as 1/κ), but (H3) and the stated Δ_0 are not, because 𝒲 has an absolute amplitude scale.
- **(f) Amplitude ratio.**
  - R enters only through r_0 = O(√ln(KR/μ_0)) and C, so Δ_0 grows like σ√ln R. The Remark's "not uniform in the amplitude ratio" is true but mild.
  - Empirically no R-dependence shows in the gate. Lattice-aligned R=10 (K=2 and K=3, `r10.py`, 48 runs, min separation 3.4–6.6σ): exact 48/48, same set 48/48, same order except 3 near-ties with g/Γ ≤ 0.04. Later entries see only the earlier blocks' misfit, so their error stays ≈ Γ (first entry ≈ 0.1Γ, weakest ≈ 1.1–1.6Γ).
  - Random R up to 8.9 (140 targets, `rdep.py`): exact 140/140; failures only at m ≤ 2.3Γ.
- **Fix.** Add "Δ>h√n" to genericity. Add "finite dictionary" and "for fixed target" to the corollary. Make the equal-wall assumption explicit. State the dependence on h/σ.

### F8. [MAJOR as wording] The corollary and Theorem 1's status
- **Where.** `supp_sep.tex` Thm (iv); `sec_theory.tex` "(H1–H4) hold, so training from ε→0 forms exactly these K units".
- **Evidence.**
  - What Steps 1–5 establish:
    - (H1) with μ_min ≥ μ_0/2.
    - (H2) in the plain Hessian sense only.
    - (H3) as a sketch for small δ (F3).
    - (H4) at fixed levels (F4).
  - Theorem 1 is itself a five-step outline whose remainder O(ε^{3/4}) comes from matching. The paper says the ε^{1/2} rate for fixed centres is "observed, not proved". The corollary inherits the outline status.
- **Fix.** See §4 (W1, W2, W6).

### F9. [MINOR] Step 2/3 details
- (i) Step 3's inductive hypothesis "duals within Cδ_ΔT_max of P_0's" is false for far-class kernels. P_0 gives them 0, while the true dual is up to μ_0/4. Restrict it to the near classes.
- (ii) "Rates within Cδ_Δ of P_0's" fails on the interval between T_k and the true τ_k. Compare duals, not rates.
- (iii) The displayed relative error "Cδ_ΔT_max/(Z|A|_minΓ(h√n/2))" has dimension of time (T_max appears twice). It should be Cδ_Δ·Z|A|_max·T_max ≤ CδR/Γ(h√n/2).
- (iv) The μ_0 in r_0 is an artifact. Far kernels only need to stay below 1−μ_0/2, not μ_0/4. Also, r_0 forces Δ_0 ≥ 2r_0σ.
- (v) Step 6's "(1+o(1))" hides ln(1+Δ²/σ²) ≈ 4.5–5.5, which exceeds ln(1/m) ≈ 1.8–3.9 at the benchmark's margins. The asymptotic formula is therefore off by about 3–4σ there.
- (vi) The Remark "a feature midway has μ_k=0 and two kernels enter together" should read "up to 2^n at a cell corner".

### F10. [MINOR] Notation
- g is both the gap (setting) and the drive g(s)=⟨rφ⟩ (Lemma mono), and g_k' appears in the Hypothesis.
- σ_k (entry shifts) clashes with the width σ.
- A_min / |A|_min are used inconsistently.
- Λ_0 is never defined in the supplement (it is Λ(1), as in `lam0 = transit_scale(L,1.0)`).
- R is the amplitude ratio in the theory but a radius in the `certificate` code.
- `Theorem 2 of the paper in full` depends on counter order.
- Two `\tbd` placeholders remain: the "Numerical test … \tbd{numbers}" is now filled, but the abstract and intro still carry `\tbd{headline numbers}`.

### F11. [MINOR] Experiment presentation
- It is 1-D only, with 61 samples on [−1,1]. It does not check intermediate block states (iii) or the margins (iv), only the set, the order, the times and the final fit. (I checked (iv): V10.)
- `terr` compares the k-th path entry with the k-th smallest T even when the sets differ, so rows ≤ 4σ mix different events.
- The md has no K=3 slope line (the slope would be 0.73 on two points).
- The 2-D record is only in the 1154-target record, not in the new experiment. My 2-D check is consistent (V11).
- "amplitude ratio ≤ 2.4" leaves the R-dependence untested in the new experiment; my R ≤ 10 runs show none.

## 3. Numerics vs theory, summary

Nothing contradicts the theorem. Exactness holds in every run: 720 + 140 + 48 + 200 (sparse lattice) + 120 (2-D) + 24 inner flows. No exits occur anywhere.

| Question | Result |
|---|---|
| Large margins but entries differ from the decoupled prediction? | Not in the regime Δ ≥ Δ_0, which the experiment cannot reach. On the sparse lattice there are targets with μ_0 = 0.13–0.32 that differ at 3–4σ. All have d_1 ≈ 1.0–1.2σ and first-order error Γ(Δ−d_1)/Γ(d_1) ≥ μ_0 (e.g. 0.18 against 0.13 at 4σ), so they are consistent with the first-order law. |
| K=3 and R=10? | Nothing breaks (F7f). |
| Margin law | Same set iff m ≳ Γ(Δ−h√n/2) (F2). |
| Entry-time error | Maximum ≈ Γ(Δ−h/2) with prefactor 0.9–1.0 (F1). |

## 4. Replacement wording (exact LaTeX)

- **W1. `sec_theory.tex`, Thm "separated regime".** Replace by:

```
\begin{theorem}[separated regime]\label{thm:sep}
Let the starting kernels form a lattice of spacing $h$ with equal walls, and let the target be $K$ Gaussian features of the kernels' width, in the continuum, with generic local data: each feature has a unique best starting kernel (margin $\mu_0>0$ over the runner-up), these kernels are distinct, and the $K$ decoupled hitting times $T_k$ are distinct (relative gap $g>0$). There is $\Delta_0=O\bigl(\sigma\sqrt{\ln(1/\min(\mu_0,g))}\bigr)$, with constants depending on $K$, $n$, $h/\sigma$ (through $e^{h^2n/16\sigma^2}$), the amplitude ratio and the unit's constants, such that for minimal separation $\Delta\ge\Delta_0$ the path has exactly $K$ entries, one starting kernel per feature, in the order of the $T_k$ and at the times $T_k(1+O(\delta_\Delta))$, where $\delta_\Delta$ is the overlap of kernels of different features within a few widths of their features, no kernel leaves, and the final block is the target with zero residual. (H1)--(H3) hold, and (H4) holds at each fixed level, so for a finite dictionary Thm.~\ref{thm:slide} gives that training from $\eps\to0$ forms exactly these $K$ units.
\end{theorem}
```
  - Add: "The proof (outline in the supplement) perturbs the decoupled path; its entry-layer step is a sketch for small overlaps, and $\Delta_0$ is not computed. With the proof's own constants it exceeds $8\sigma$ for the benchmark's lattice, above every separation tested below."
- **W2. `abstract.tex`.**
  - Replace "which we prove under stated hypotheses" by "for which we give a proof outline under stated hypotheses".
  - Replace "(proved for generic features beyond a separation; observed on 1154 benchmark targets down to $2\sigma$)" by "(a perturbative argument, outlined in the supplement, covers generic features beyond a separation that exceeds our benchmarks'; observed on 1154 benchmark targets down to $2\sigma$)".
- **W3. `sec_intro.tex`, contribution 2.** Replace "For generic separated features the path ends with exactly one unit per feature, at the feature's centre and amplitude (Thm.~\ref{thm:sep});" by "For generic features beyond a separation $\Delta_0$ the path ends with exactly one unit per feature, at the feature's centre and amplitude (Thm.~\ref{thm:sep}, proof outline in the supplement; $\Delta_0$ is not computed and, on the proof's own constants, exceeds our benchmarks' separations);".
- **W4. Scope.** Replace "a proof of exact recovery below the separated regime of Thm.~\ref{thm:sep} (there the record is numerical)" by "a proof of exact recovery at the separations we test ($2$--$8\sigma$): the separation $\Delta_0$ of Thm.~\ref{thm:sep} lies above them on its own constants, so the benchmark record is numerical throughout".
- **W5. `sec_conclusion.tex` (ii).** Replace "exact recovery is proved for generic features beyond a separation $\Delta_0$ that is not sharp and not uniform in the amplitude ratio (below it the record is numerical)" by "exact recovery has a proof outline for generic features beyond a separation $\Delta_0$ that is not computed, lies above our benchmarks' separations, and depends on the features' position against the starting grid and, logarithmically, on the amplitude ratio (at the tested separations the record is numerical)".
- **W6. `sec_theory.tex` record sentence.**
  - "always" → "on all 1154".
  - "on all 302 targets whose decoupled margins exceed $5\Gamma(\Delta)$" → "on all 356 targets whose decoupled margins exceed $\Gamma(\Delta-h/2)$, the first-order cross-talk of a kernel half a cell from its feature".
  - "the entry times follow $\Gamma(\Delta)$ (log--log slope $1.01$" → "the largest relative error of the entry times is $0.9$--$1.0\,\Gamma(\Delta-h/2)$ at $5$--$8\sigma$".
- **W7. `supp_sep.tex` l. 38.**
  - "(theory: $O(\Gamma)$)" → "(the first-order cross-talk; the theorem proves $O(\delta_\Delta)$)".
  - "The margin condition is quantitative: no target with margins above $5\Gamma$ departs from the decoupled prediction." → "On this lattice no target with margins above $\Gamma(\Delta-h/2)$ departs from the decoupled prediction (356 of 356); the constant $5$ is specific to this lattice."
- **W8. `supp_sep.tex` Thm (iv).** "(H1)--(H4) hold, with $\mu_{\min}\ge\mu_0/2$" → "(H1)--(H3) hold, and (H4) at every fixed level $\theta$ with a threshold depending on $\theta$, with $\mu_{\min}\ge\mu_0/2$; so by Thm.~\ref{thm:slide-full} (whose proof is an outline), for a finite dictionary, the networks trained from $\eps$ form exactly these $K$ units as $\eps\to0$".
- **W9. `supp_sep.tex` Setting.** "(other walls change the notation only)" → "(for unequal walls the best kernel minimises $\rho_j/|a_{jk}|$, not the distance, and the constants carry $\rho_{\max}/\rho_{\min}$)". Main-text Prop "one feature": "the unit of largest correlation" → "the unit of smallest hitting time (largest correlation for equal walls)".
- **W10. CLAIMS.md D15 and `supp_sep.tex` heading.** "upgraded to a THEOREM" → "theorem with a perturbative proof outline; $\Delta_0$ not computed".

## 5. Not checked

1. The proof of Lemma 3(iii) beyond the sketch: Gronwall or Lyapunov constants, and the far-field ξ⁻³ estimate in the variational equation.
2. Theorem 1's own outline (Mosco, Attouch, Brezis–Pazy, layer matching). The corollary is conditional on it.
3. Networks trained at finite ε. I ran only the limit object and the ε-free inner flow.
4. The 1154-target record. I did not re-run it.
5. n ≥ 3, non-cubic lattices, unequal walls, K ≥ 4, R > 10.
6. The polish-acceptance guard in `sliding_iss_nd`. It is a heuristic that identifies the stationary point, and I did not audit it. My inner-flow check is independent of it.
7. LaTeX compilation and cross-references.
8. Constants. The proof gives none, so every Δ_0 figure above is a lower estimate under my reading of the displayed conditions (C = 4KR/Γ(h√n/2)); larger constants only worsen it.
9. Edge truncation of ⟨φ²⟩ by the finite window (features ≥ 0.3 inside the [−1,1] window). Not quantified.

**Private artifacts** (in the `referee_IS2/` directory named at the top): `gate.py`, `gate2.py`, `d0est.py`, `errscale.py`, `hess.py`, `r10.py`/`r10.log`/`r10_out.json`, `rdep.py`/`rdep_out.json`, `sparse.py`/`sparse_0.3.json`, `twod.py`/`twod.log`, `inner.py`, `inner2.py`, `margin.py`, `f1.py`.
