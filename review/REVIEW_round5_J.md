# Adversarial review of Theorem J (PARADIGM.md §6d)

Reviewed: PARADIGM.md §6d (statement (a)–(d), Steps 1–6, Numerics, Remark J.1; lines ≈ 2154–2370, shifted by +3 during the
review because the file was being edited elsewhere), Proposition H.7(a)(b), Remarks H.7a and H.0; `src/core/periodic.py`,
`src/core/averaging.py`, `src/core/langevin.py`; `tests/test_periodic.py`; `experiments/exp_periodic_orbit.py`,
`results/periodic_orbit/{table.md,raw.npz}`. No repository file was modified. Scratch scripts and outputs:
`/tmp/claude-0/-home-user-Lyapunov-Convergence-Kernels-Showcase/bab27c21-f0c6-5e06-ae9a-246f40b4fb64/scratchpad/s1_raw.py`
(raw.npz re-analysis), `s2_instB.py` → `s2_out.txt` (slaving, oscillation, variational monodromy, Liouville, Riccati, b_osc),
`s3_saddle.py` → `s3_out.txt` (second collapsed minimum), `s4_misc.py` → `s4_out.txt` (transient, σ₀ = 0 degeneracy, basin),
`s6_saddle3.py` → `s6_out.txt` (hyperbolic index-1 saddle), all in the same directory.

## Verdict

Theorem J holds as stated, in all four parts, for the stated hypotheses (f ∈ C², x_m ∈ C², saturated features, σ = σ₀/k_e,
hyperbolic Θ*). I found no counterexample; every precise consequence I could test numerically on instance B (n = 1, m = 2,
state dimension 7) came out with the predicted order and, where a constant was predicted, the predicted constant:

- slaving (Step 2): e(t) − e_{Θ_p(t)}(t) = 0.30·ε³ at k_e = 20, 40, 80 (ratios 7.9, 8.1 per doubling);
- oscillation (a): Θ_p(t) − Θ̄_p − εu(Θ*,t) = 1.1·ε² (ratios 3.9, 4.0);
- location (d): k_e²·‖Θ̄_p − Θ* − prediction‖ → 0.57 (A), 5.08 (B): the remainder is O(k_e⁻²);
- Floquet (b): k_e²·(exponent − ν_j/k_e) converges for every j; the fast channel, inaccessible to the FD monodromy, is
  recovered through Liouville's identity: log|λ_fast| = −126.13 = −(k_e + 0.074)T₀ at k_e = 20 (bound e^{−(k_e−L)T₀} with any L ≥ 0);
- Riccati (Step 5): the slow invariant subspace of the monodromy is the graph δe = S(0)δΘ with k_e‖S(0)‖ = 0.92–0.94 and
  ‖S(0) − ε∂_Θr‖/‖ε∂_Θr‖ = 0.067, 0.033, 0.016 (halving: S = ε∂_Θr + O(ε²));
- unstable case (b),(d), never tested by the experiment: at an index-1 hyperbolic saddle (σ₀ = 0.003) the orbit exists, has
  exactly one multiplier > 1 (1.00168 vs predicted 1.00185 at k_e = 20) and k_e²‖Θ̄_p − Θ* − pred‖ = 3.80 → 3.49.

The proof, however, is a sketch with one genuine gap: Steps 4 and 6 use *additive* O(ε³) remainder bounds where *Lipschitz
(C¹) control* of the remainders is needed, both for the implicit-function argument on an ε-independent ball (existence and
uniqueness) and for the contraction (basin). With additive bounds alone, Step 6 only gives convergence into an O(ε²) ball
around the fixed point, not to it. The gap is routine to close with f ∈ C² (sketch below), but as written the proof does not
establish (a)'s uniqueness nor (c). Besides this there are one wrong intermediate estimate (Step 1's transient), one false
number in the Numerics paragraph, one imprecise claim in Remark J.1, and cosmetic items.

## Issues

| id | severity | issue | evidence | proposed fix |
|---|---|---|---|---|
| J-1 | major (fixable) | Steps 4, 6: remainders are bounded only in sup norm. "Ã_ε = ⟨H⟩ + O(ε)" does not give DÃ_ε = D⟨H⟩ + O(ε), which the IFT/uniqueness on a ball of ε-independent radius needs; "nonlinear remainder εT₀·O(‖Y₂−Y₂*‖²) + O(ε³)" with an additive O(ε³) is not a contraction estimate (crude Lipschitz bound of the remainder is O(ε), same order as the contraction gain). | Step 4 (line ≈2263), Step 6 (≈2298). The fixed-point structure in (Y₂, w) (first solve w = P_w(Y₂, w), then Ã_ε with w*(Y₂)) is also not spelled out. | Prove the remainders of (S)/(N) are O(ε³) *in the Lipschitz norm* in (Θ, ζ), ζ := e − e_Θ(t): with f ∈ C², the Taylor remainder Q(h) = r(x_m+h) − r(x_m) − ∂_x r(x_m)h has Lip_h Q = O(ε) on ‖h‖ ≤ Cε, e_Θ is O(ε)-Lipschitz in Θ, and the 1/k_e gain of the fast equation gives Lip_Θ(z_Θ) = O(ε³) for the ε³-remainder of e_Θ; then Y ↦ Y − D⟨H⟩(Θ*)⁻¹Ã_ε(Y) is a contraction on B(Θ*, ρ₀) (Lipschitz IFT) and P is a contraction in the product norm ‖·‖_{G⁻¹} + ‖ζ‖ (couplings: Lip_ζ(P_Y) = O(ε), Lip_Y(P_ζ) = O(ε³)). Alternatively assume f ∈ C³ and cite the C¹-remainder version of the averaging theorem (SVM Thm 2.9.x with one more derivative). |
| J-2 | minor | Step 1: "absorbs every solution of the tube within t₁ = O(ε log(1/ε)), during which Θ moves by O(ε² log(1/ε))" is false when ‖e(0)‖ = O(1): ∫₀^{t₁}‖w‖dt ≈ (2/k_e)‖w(0)‖, so Θ moves by O(ε)‖e(0)‖. Step 6 itself uses the correct O(ε)‖e(0) − e_p(0)‖. | Measured (s4, instance B, start at (x_m(0)+e₀, Θ_p(0)), t₁ = 0.5): ‖ΔΘ‖ = 1.88, 1.74 × e₀/k_e at k_e = 20 (e₀ = 0.3, 0.6) and 1.93, 1.78 × e₀/k_e at k_e = 40; ε² log(1/ε) = 7.5e-3, 2.3e-3 is 4–7× smaller and does not scale with e₀. | Replace by "Θ moves by O(ε)(‖e(0)‖ + ε log(1/ε))", and note that this is why 𝒦 is taken with radius 2ρ₀. Harmless for the argument. |
| J-3 | minor | Step 2/3: (S) and (N) are claimed "uniformly O(ε³) on W_ε × 𝒦", but on W_ε only w = O(ε²) is known; the quasi-static slaving w = w_qs + O(ε³) needs a further O(ε log 1/ε) transient. From W_ε the Poincaré-map remainder is O(ε³ log(1/ε)). | Step 2 (≈2234–2240), Step 3 (≈2252), Step 4 (≈2258). The claim itself is true on the thinner set: on the orbit e − e_{Θ_p(t)}(t) = 0.295, 0.300, 0.294 × ε³ (k_e = 20, 40, 80). | Use the fast variable ζ := e − e_Θ(t) (or w − w_qs(Θ,t)) and the forward-invariant set {‖ζ‖ ≤ Cε³}; or remark that a periodic solution has no transient so the location argument is unaffected, and absorb the log into Step 6's κ. |
| J-4 | minor | Step 5/(b): "Its eigenvalues are 1 − εT₀ν_j + O(ε²)" for a matrix I + εA + O(ε²) needs A diagonalisable (Bauer–Fike); A = −T₀G_m∇²J_m is similar to a symmetric matrix, so it holds, but at collapsed Θ* (instance A) the ν_j are (m−1)-fold and the O(ε²) term may split a multiplet into a complex pair. | Instance A: ν₁ = ν₂ = 0.0389 (2-fold), and the text's exponents are −log|λ|/T₀. | Add one sentence: "G_m∇²J_m(Θ*) is diagonalisable (similar to G_m^{1/2}∇²J_mG_m^{1/2}), so by Bauer–Fike the eigenvalues of I + εT₀D⟨H⟩ + O(ε²) are 1 − εT₀ν_j + O(ε²), multiplicities counted; exponents are understood as −log|λ|/T₀." |
| J-5 | minor | Numerics: "residual ≤ 3·10⁻¹³ in every case" is false. | `results/periodic_orbit/table.md`: instance A, k_e = 160: 8.6e-12 (also B, k_e = 80: 3.4e-13). The table's own header says < 1e-10. | Write "≤ 10⁻¹¹ in every case (8.6·10⁻¹² at k_e = 160)". |
| J-6 | minor | Remark J.1, bullet 3: "with σ₀ = 0 the critical points of mL̄_m are degenerate along the symmetric directions of the permutation group only when kernels coincide" names the wrong directions and over-generalises. At a collapsed point the separation Hessian is σ₀G⁻¹ − R₂ (as `periodic.first_variation_hessian` states), (m−1)-fold. For n = 1 and σ₀ = 0 the first-order conditions give ⟨rφ⟩ = ∂_γ⟨rφ⟩ = ∂_κ⟨rφ⟩ = 0, so the whole α-row of R₂ vanishes and the (m−1) *amplitude-redistribution* directions are exactly flat (the other separation eigenvalues are not: −0.031, +0.040). For n ≥ 2 only the αα-block vanishes and the collapsed point is generically non-degenerate even at σ₀ = 0. | s4 (H): σ₀ = 0, collapsed critical point of J_2 (a = 1.691, c = 1.049, s = 0.472): eig ∇²J_2 = (−0.0308, 2.6e-13, 0.0380, 0.0398, 0.135, 0.829); R₂[α,α] = 5.6e-14; same at the negative bump (eigenvalue 1.1e-13). | State: "for n = 1 the collapsed critical points of mL̄_m are degenerate in the m − 1 amplitude-redistribution directions (⟨rφ⟩ = 0 at a critical point without leak); the leak lifts this, σ_crit of `periodic.py` being the marginal value; for n ≥ 2 the degeneracy is generically absent." |
| J-7 | minor | (b), error channel: the experiment's "below numerical resolution, consistent with (b)" is vacuous — nothing about the n fast multipliers is checked. (c) is not checked at all (the text says so). | FD monodromy's smallest multiplier is noise: 2.7e-9, −9.9e-9, −7.5e-9 at k_e = 20, 40, 80 (max‖DP_fd − M_variational‖ = 3e-8). Liouville: log|λ_fast| = ∫₀^{T₀}tr DF dt − Σ_j log λ_j = −126.13, −251.86, −503.21 = −(k_e − L_eff)T₀ with L_eff = −0.074, −0.084, −0.089 (⟨∂_x r⟩_orbit = −0.12, −0.11, −0.10). Basin run (k_e = 20, ‖Θ(0) − Θ*‖_{G⁻¹} = 0.5 = 25× the offset, e(0) = 0.3): distance 0.5 → 0.104 (20 periods) → 0.0144 (60); per-period factors 0.87 → 0.951 → 0.953 inside the predicted band [e^{−T₀ν₆/k_e}, e^{−T₀ν₁/k_e}] = [0.795, 0.971]. | Add the trace identity to `exp_periodic_orbit.py`/`test_periodic.py` as the check of the error channel (it is exact and cheap: tr DF along the orbit = −nk_e + tr ∂_x r + tr ∂_ΘF), and report L_eff; optionally add the basin run. |
| J-8 | cosmetic | (i) "∂_ΘF = O(ε) (F is O(‖e‖ + ε) and e = O(ε))": right conclusion, wrong reason (smallness of F does not bound its derivative; F is linear in e and ∂_{ϑ_j}F_i = 0 for j ≠ i). (ii) ρ₀ measured in ‖·‖ in (a) and in ‖·‖_{G⁻¹} in (c). (iii) Step 1/6 work on {‖w‖ ≤ 1}; if R < 1 this leaves the tube — intersect with {‖e‖ ≤ R}. (iv) Θ ↦ Y₂ is a diffeomorphism onto an O(ε)-neighbourhood of 𝒦, not "of 𝒦". (v) "sup‖ẋ_m‖/s ≈ 3–10": with the minimiser's s ≈ 0.77–0.78 it is 2.6; 10 is 2/s_min. (vi) Remark J.1 bullet 1 "the stationary states of the real closed loop are O(1/k_e)-perturbations of the critical points" states a converse the theorem does not prove (only orbits near hyperbolic critical points are covered). (vii) The 2-fold "separation pair with the smallest exponent (0.0018 at k_e = 20)" is verified (pred 0.00194, measured 0.00180). | lines ≈2278, 2190/2203, 2230/2299, 2252, 2335, 2357 | Edit the phrases. |

## Details

### (a) Proof gaps, step by step

**Step 1.** Forward invariance of {‖e‖ ≤ R} for k_e > M_r/R: correct (eᵀr − k_e‖e‖² ≤ RM_r − k_eR² < 0), with the usual
circularity (M_r is a sup over Θ ∈ 𝒦, so Θ must stay in 𝒦 — true by continuation since Θ moves O(ε) in the transient, see J-2).
The w-equation is derived correctly: with w := e − εr_Θ(x), r_Θ(x) − k_e e = −k_e w and ẋ = ẋ_m + ė = ẋ_m − k_e w give
ẇ = −k_e w + ∂_x r·w − ε∂_x r·ẋ_m − ε∂_Θ r·Θ̇ (using εk_e = 1). The estimate d/dt‖w‖ ≤ −(k_e − L_r − O(ε))‖w‖ + ε(L_r sup‖ẋ_m‖ + O(ε))
≤ −(k_e/2)‖w‖ + C₁ε for k_e ≥ 2L_r + C is right, hence ‖w(t)‖ ≤ e^{−k_e t/2}‖w(0)‖ + 2C₁ε² and W_ε forward invariant. Wrong
only: the transient displacement of Θ (J-2).

**Step 2.** The quasi-static solution w_qs = −ε²∂_x r(x_m)ẋ_m + O(ε³) (the Θ̇-term is O(ε²) forcing, hence O(ε³) in w), the lag
O(ε·ε²) (forcing varies at O(1) rate since ẋ = ẋ_m + O(ε), needs f ∈ C²), and e = εr(x_m) + ε²∂_x r r + w + O(ε³) = e_Θ(t) + O(ε³)
with e_Θ the frozen-Θ periodic error of H.7(a): all correct; verified numerically to the constant (0.30ε³, see (b)(iii)).
(S) with ⟨H₂⟩_i = G∇_ϑ(U_even + U_odd)[μ_Θ](ϑ_i): consistent with H.7(b) — the ϑ-gradient acts on the explicit ϑ_i of Φ̃ with
r_Θ frozen, exactly as `offset_prediction` differentiates (`grad_at_particles` perturbs the evaluation rows with `fm` fixed).
Caveat J-3 about where (S) holds uniformly.

**Step 3.** The normal-form computation checks: (I + εD_Yu)Ẏ = ε⟨H⟩ + ε²[H₂ + D_ΘH·u] + O(ε³) ⇒ Ẏ = ε⟨H⟩ + ε²[H₂ + D_ΘH·u − D_Yu·⟨H⟩]
+ O(ε³); ⟨D_Yu·⟨H⟩(Y)⟩_t = D_Y⟨u⟩·⟨H⟩ = 0 because ⟨H⟩(Y) is t-independent; the second transformation removes the fluctuating
part. H, H₂ are C^∞ in Θ (Gaussian features; f enters only through its values along x_m), so the regularity needed here is
free. The near-identity map is a diffeomorphism for small ε (J-8(iv)). The remainder's *Lipschitz* control is what is
missing (J-1).

**Step 4.** The time-T₀ expansion Y₂(T₀) = Y₂(0) + εT₀⟨H⟩ + ε²T₀Φ₂ + ε²(T₀²/2)D⟨H⟩⟨H⟩ + O(ε³) is correct (Picard to second
order). The IFT needs DÃ_ε = D⟨H⟩ + O(ε) uniformly on the ball: not established (J-1); given it, ρ₀ is set by
‖D⟨H⟩(Θ*)⁻¹‖ and Lip(D⟨H⟩), hence ε-independent, and the argument is sign-independent (unstable case fine). The derivation
Y₂* = Θ* − εD⟨H⟩⁻¹Φ₂(Θ*) + O(ε²) is right ((T₀/2)D⟨H⟩⟨H⟩(Y₂*) = O(ε)), and the factors of k_e match: J = εD⟨H⟩(Θ*),
b_frozen + b_osc = ε²Φ₂(Θ*), so J⁻¹b = εD⟨H⟩⁻¹Φ₂; the code computes J as the FD Jacobian of (1/k_e)(−G∇J_m) and b with the
explicit 1/k_e², consistent. Period mean: Y₂(t) = Y₂* + O(ε²) on the orbit (Ẏ₂ = O(ε²)), and ε⟨u(Y₂(t),t)⟩ = ε⟨u(Y₂*,t)⟩ + O(ε³) = O(ε³),
so Θ̄_p = Y₂* + O(ε²): correct.

**Step 5.** Riccati: substituting δe = η + SδΘ gives the δΘ-coefficient (∂_x r − k_e)S + ∂_Θ r − Ṡ − S∂_eF S − S∂_ΘF, whose
vanishing is exactly the stated equation, and η̇ = [(∂_x r − k_e) − S∂_eF]η, δΘ̇ = ∂_eF η + [∂_eF S + ∂_ΘF]δΘ: brackets recomputed,
correct. Existence of the small periodic S for large k_e: S = 𝒯(g(S)) with ‖𝒯‖ ≤ 1/k_e and Lip g = L_r + O(1/k_e) on the ball
‖S‖ ≤ 2‖∂_Θr‖_∞/k_e: a contraction, fine. Block-triangular time-varying systems have block-triangular fundamental matrices
(X₁₂ ≡ 0 solves Ẍ₁₂ = AX₁₂, X₁₂(0) = 0), the coordinate change is T₀-periodic, so the multipliers are the union: correct.
The δΘ-block: S = ε∂_Θ r + O(ε²), ∂_eF S = εG∇_ϑ[(∂_{ϑ_j}r)ᵀΦ̃_i] + O(ε²), ∂_ΘF = block-diagonal ε∂_{ϑ_i}G∇_ϑ[rᵀΦ̃_i] − εσ₀I + O(ε²),
whose sum is εD_ΘH(Θ_p(t),t) + O(ε²) = εD_ΘH(Θ*,t) + O(ε²) because Θ_p = Θ* + O(ε): correct; monodromy I + εT₀D⟨H⟩ + O(ε²),
eigenvalues via Bauer–Fike (J-4). Numerically S(0) is a real graph with S(0) = ε∂_Θr(x_m(0)) + O(ε²) ((b)(iv)).

**Step 6.** Norm of I − εT₀G_m^{1/2}∇²J_mG_m^{1/2} is 1 − εT₀λ_* once εT₀ν_{mp} < 1: correct. Coupling to w is O(ε)‖w‖ (∫e^{−k_et/2}dt = 2ε)
and w's dependence on Y₂(0) is O(ε²)-Lipschitz: a product-norm contraction exists *given* J-1. The rate μ = (λ_*/k_e)(1 − κ)
follows from the per-period factor 1 − εT₀λ_*(1−κ) ≈ e^{−εT₀λ_*(1−κ)}: correct. The transient handling (O(ε)‖e(0) − e_p(0)‖)
is correct and contradicts Step 1's O(ε² log 1/ε) (J-2).

**Statement (b)** "n multipliers of modulus ≤ e^{−(k_e − L)T₀}" is consistent with the experiment's "below resolution" (at k_e = 10 the
bound is already e^{−63−}); the Liouville computation makes it quantitative (J-7). **Statement (c)**: the derivation of μ is
right; ρ₀ independent of k_e follows from the uniform contraction (modulo J-1).

### (b) Numerical checks (instance B unless stated; scripts s1, s2, s3, s4, s6)

(i) Scaling. From raw.npz: instance A, k_e·‖Θ̄_p − Θ*‖ = 0.409, 0.382, 0.369, 0.362, 0.359 and k_e²·‖Θ̄_p − Θ* − full‖ = 0.603,
0.592, 0.580, 0.573, 0.569 (k_e = 10…160): offset O(1/k_e), remainder O(k_e⁻²). Instance B: k_e·‖off‖ = 0.645, 0.444, 0.335,
0.281, 0.254 (→ ≈0.23 = k_e‖full‖) and k_e²·‖off − full‖ = 4.56, 4.93, 5.05, 5.07, 5.08: same orders, but remainder/leading ≈ 22/k_e,
so the first-order location is only quantitative from k_e ≈ 80 — exactly as the text says; the table's fitted exponent −1.33 for
B (not quoted in the text) is this pre-asymptotic effect. b_osc is 30.5% of the full prediction in A at every k_e (text: "O(1)
fraction", 35% vector error of the frozen part: 0.398 → 0.342). All numbers of the Numerics paragraph match table.md except J-5.

(ii) Floquet exponents. k_e²·(exponent − ν_j/k_e) converges to a constant for every j (A, smallest 2-fold pair: −0.032, −0.057,
−0.070, −0.077, −0.080; B, largest: −0.69, −0.47, −0.34, −0.27, −0.24, still drifting, i.e. a large k_e⁻³ term). Variational
monodromy (56-dim DOP853, rtol 1e-11) agrees with the FD monodromy to 3e-8 entrywise and gives the same parameter exponents to
5 digits, so the table's Floquet errors are genuine O(k_e⁻²) effects, not FD noise.

(iii) Slaving. max_t |e(t) − e_{Θ_p(t)}(t)| with e_Θ the exact frozen periodic error (RK4, dt = T₀/8000, 3 periods, Θ frozen at
Θ_p(t), same t): 3.69e-5, 4.68e-6, 5.75e-7 at k_e = 20, 40, 80 = 0.295, 0.300, 0.294 × ε³ (ratios 7.88, 8.14). Against the
second-order expansion εr + ε²(∂_x r r − ṙ): 9.96, 10.15, 10.02 × ε³; against εr alone: 2.15, 2.38, 2.50 × ε². Step 2's claim is
confirmed with its constant.

(iv) Monodromy structure. FD monodromy's 7th multiplier: 2.7e-9, −9.9e-9, −7.5e-9 (noise; variational: 0, −1e-18, 5e-19).
Liouville: ∫tr DF dt = −126.71, −252.15, −503.36; log|λ_fast| = −126.13, −251.86, −503.21 = −(k_e − L_eff)T₀, L_eff = −0.074,
−0.084, −0.089 (so |λ_fast| < e^{−k_eT₀}; ⟨∂_x r⟩ along the orbit = −0.120, −0.107, −0.101, sup|∂_x r| = 1.49). The six slow
multipliers are real and ≈ exp(−T₀ν_j/k_e) with relative errors 0.07–3.2% (k_e = 20), 0.04–1.2% (40), 0.04–0.5% (80). Riccati:
the 6-dimensional slow invariant subspace is the graph δe = S(0)δΘ (imaginary part 0), k_e‖S(0)‖ = 0.915, 0.933, 0.940,
‖S(0) − ε∂_Θr(x_m(0))‖/‖ε∂_Θr‖ = 0.067, 0.033, 0.016.

(v) b_osc = 0 on instance B: `oscillation_correction` gives max|k_e²b_osc| = 7.1e-15 (the frozen coefficients are 0.004–0.06);
H(Θ*, π − t) = H(Θ*, t) to 3.4e-15 on the grid; U_odd = 5e-17 (n = 1). The reflection argument is correct (u(π − t) = −u(t) follows
from H even about π/2 and ⟨u⟩ = 0). Note it needs x_m(π − t) = x_m(t), i.e. a sweep symmetric about its turning points, not just n = 1.

(vi) Oscillation (a): max_t |Θ_p(t) − Θ̄_p − εu(Θ*,t)| = 2.71e-3, 6.96e-4, 1.75e-4 = 1.08, 1.11, 1.12 × ε² against amplitudes
1.1e-2, 5.6e-3, 2.8e-3.

(vii) Unstable case (not in the experiment). Collapsed critical point of J_2 at σ₀ = 0.003 (a = 1.6345, c = 1.0459, s = 0.4889),
ν = (−0.00587, 0.00491, 0.0585, 0.0690, 0.216, 0.804), index 1, hyperbolic. Newton shooting converges (residual 5e-15 at k_e = 20,
7e-12 at 40), the orbit stays collapsed (4.7e-13), has exactly one multiplier > 1: 1.00168 vs exp(T₀·0.00587/20) = 1.00185 (k_e = 20),
1.00089 vs 1.00092 (40); k_e²(exponent − ν₁/k_e) = 0.0106, 0.0094. (d): k_e²‖off − full‖ = 3.80, 3.49, vector error 1.23 → 0.57,
cosine 0.006 → 0.84 — the first-order location is orthogonal to the true offset at k_e = 20 because the leading offset is small
(k_e‖off‖ ≈ 0.11) relative to the O(k_e⁻²) constant; consistent with the theorem, but a useful caution for Remark J.1's second
bullet. A second collapsed *minimum* of instance B (σ₀ = 0.05, a = 1.2022, c = 0.9875, s = 0.6540, ν_min = 0.043) behaves as
the main instance: max rel Floquet error 3.5% → 1.2%, k_e²‖off − full‖ = 5.49 → 5.53.

### (c) Hypotheses and scope

- **Scaling σ = σ₀/k_e** is essential for the statement, not for a periodic orbit: with fixed σ the parameter equation is
  Θ̇ = −σ(Θ − Θ⁰) + O(ε) — the leak is not slow, the orbit is a regular perturbation of Θ⁰ (always stable), and the critical points
  of J_m play no role; with fixed σ ≪ 1/k_e one is in the σ₀ = 0 case of the theorem (critical points of mL̄_m, see J-6 for the
  n = 1 degeneracy). The scaling is what puts mL̄_m and the leak at the same order.
- **f ∈ C² vs C³.** All sup-norm remainders need only f ∈ C² and x_m ∈ C² (the O(ε³) term of e_Θ uses ∂²_x f and ẍ_m through ė₂).
  The Lipschitz control of J-1 also holds with C² (the second-order Taylor remainder of r is O(ε)-Lipschitz on an O(ε) ball). The
  classical C¹ implicit-function theorem would need f ∈ C³. H and H₂ are smooth in Θ regardless. So the stated hypothesis is
  sufficient, but the text must say which route it takes.
- **Hyperbolicity in the unstable case**: "∇²J_m(Θ*) non-singular" suffices; the IFT is sign-blind and the Poincaré multipliers
  1 − εT₀ν_j + O(ε²) avoid the unit circle for ν_j ≠ 0 and k_e ≥ k₀(Θ*). Confirmed numerically ((b)(vii)).
- **Instance A's collapsed minimiser** is non-degenerate: eig ∇²J_m = 0.01479 (×2), 0.1429 (×2), 0.285, 0.388, 0.523 (×2), …,
  ν₁ = ν₂ = 0.0389; the (m−1) = 2-fold multiplets are the separation modes, as claimed. The permutation symmetry causes
  multiplicity, not singularity; the proof is unaffected except for the eigenvalue-perturbation remark (J-4). The collapsed
  subspace is invariant under (CL), so the orbit is itself collapsed (verified: 1e-13), which is implicit in the text.
- **Overclaims.** "Makes Remark H.7a a theorem for this system" is justified for the *period-mean location of the periodic orbit*
  (the orientation experiment does measure period means, `Th[-k:].mean(0)` over one period), modulo J-1; Remark H.7a's dynamical
  statement (finite-horizon second-order averaging of Y) is not what is proved and the text should say "its stationary
  prediction". Remark J.1 bullet 1 overstates (J-8(vi)). "residual ≤ 3·10⁻¹³" is false (J-5). The Checked paragraph and the
  verification-map row match the test's assertions (<0.65, <0.6, <0.02 factors) and the table.

### Summary of what would make the section pass as a proof

1. Add the Lipschitz (or C¹) remainder estimate and the two-step fixed-point argument in (Y₂, ζ) (J-1); define ζ := e − e_Θ(t) with
   its O(ε³) invariant set (J-3). 2. Fix the transient estimate (J-2) and the Bauer–Fike sentence (J-4). 3. Correct the residual
   number (J-5) and rewrite Remark J.1 bullet 3 (J-6). 4. Optionally add the Liouville check of the error channel (J-7).
