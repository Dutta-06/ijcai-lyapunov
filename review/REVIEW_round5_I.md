# Adversarial referee report on Theorem I (PARADIGM.md §6c, "reciprocity at order k_e⁻²")

Scope: PARADIGM.md Proposition H.7(e), Remark H.7a, Theorem I with proof, the *Checked* paragraph and Remark I.1;
`src/core/averaging.py`, `src/core/nonreciprocity.py`, `tests/test_nonreciprocity.py` (6/6 pass, 0.7 s).
Line numbers refer to the working copy at review time (2233 lines; the Theorem I block is byte-identical to HEAD
`ca69ced` shifted by +2 lines, because §6d/Theorem J was being inserted concurrently). Nothing in the repository was
modified. Scratch scripts and raw outputs: `scriptA_cyclic.py` … `scriptH_n1_perfectfit_triangle.py` in this directory.

## Verdict

Theorem I holds in parts (i), (ii), (iii) and (v), and in the second half of (iv) (K ≡ 0 ⟺ T_𝒪 = 0 ⟺ w_𝒪 = 0 a.e.
for n = 2). Every formula was re-derived by hand and checked against the *exact* second variation of the averaged
potential (finite differences on the periodic error, Richardson k_e = 80/160), including at adversarial parameters
(saturated weights, s = s_min, centres outside the orbit): the bilinear form aᵀB(ξ,ξ')a' agrees with the generic code to
6e-15 and with the exact antisymmetric second variation to 2e-4–9e-4 relative; the pairwise (zero-feature) criterion
is confirmed at the exact P by the triangle-of-measures circulation (ratio to the centroid formula 0.986 at k_e = 80,
0.993 at k_e = 160); the quoted cyclic sums −0.0704 and −0.0464 (and −0.0711/−0.0471 exact, −0.089/+0.043 split) are
those of `test_averaging.py`. **One sub-claim is false as stated:** (iv) "If r_μ = 0 ρ-a.e., then A = 2K" (line 1757;
proof line 1791 "then A_even ≡ 0"), together with its echoes (line 1763 "at perfect fit the semi-gradient defect has
vanished", summary row I line 89, H.7(e) bullet line 1606, Remark I.1 lines 1823–1827, line 2157). For n ≥ 2 the curl
term −∫Φ̃ᵀ(Dr_μ − Dr_μᵀ)Φ̃' dρ of (ii) survives when the residual vanishes on the orbit but its Jacobian does not
(vanishing on a curve constrains only the tangential derivative). An explicit instance (f := f̂_μ + h·e₁ with h = 0
on the ellipse) has sup_𝒪|r_μ| = 4e-16 yet A_even ≠ 0 at 38 % of the kernel, and the exact second variation confirms
A = A_even + 2K (1.6e-4) rather than 2K. The test `test_perfect_fit_n2_leaves_only_the_circulation_kernel` cannot see
this because it takes f := f̂_μ, i.e. r_μ ≡ 0 on all of ℝ². The fix is local and keeps the message: add the hypothesis
Dr_μ = Dr_μᵀ ρ-a.e. (e.g. r_μ ≡ 0 near the orbit), or, better, state the full perfect-fit characterisation
A ≡ 0 ⟺ (curl r_μ = 0 ρ-a.e. and T_𝒪 = 0), which the same Gaussian-convolution argument proves. A second, smaller
problem is one of definition: "reciprocal at order k_e⁻² at a frozen μ" is never defined; what is proved is the
pointwise statement A_μ ≡ 0 (closedness of the 1-form P at μ), which is necessary but not sufficient for P to be a
first variation near μ. At a perfect-fit μ in n = 1 the H.7(e) triangle test still fails on every neighbourhood
(measured ratio 1.0002 to the prediction −(ε²/2)·εA_ν), so the 1-d averaged loop is never a gradient flow up to
O(k_e⁻³) on any open set of measures; the sentences "reciprocal … exactly when μ interpolates f" (line 1753) and
"reciprocity at order k_e⁻² is equivalent to A_even + 2K ≡ 0" (line 2156) need that definition to be true.

## Issues

| id | severity | issue | evidence | proposed fix |
|---|---|---|---|---|
| I-1 | **major** | (iv) first claim false for n ≥ 2: r_μ = 0 ρ-a.e. does not give A_even ≡ 0; the curl term −∫Φ̃ᵀ(Dr_μ−Dr_μᵀ)Φ̃'dρ survives (lines 1757, 1763, 1791; echoed at 89, 1606, 1823–1827, 2157). The *Checked* claim "at perfect fit for n = 2, B = 2κI (1e-13)" (line 1815) holds only for r_μ ≡ 0 on ℝⁿ. | `scriptB_perfectfit.py`: f := f̂_μ + h e₁, h = (x₁/1.2)²+(x₂/0.8)²−1 on the test ellipse; sup_𝒪|r_μ| = 4.4e-16, sup_𝒪‖Dr−Drᵀ‖ = 2.5; B − 2κI has entries 0.02–0.12 vs |2κ| = 0.05–0.13 and equals −∫φφ'∂₂h dρ·[[0,1],[−1,0]] to all digits; exact Richardson k_e²(S−Sᵀ) matches A_even+2K to 1.6e-4 relative and differs from 2K by 37.7 %. | Replace (iv) by: "If r_μ = 0 ρ-a.e., then A_even(ϑ,ϑ') = −∫Φ̃_ϑᵀ(Dr_μ − Dr_μᵀ)Φ̃_ϑ'dρ (only the Lie-bracket term dies), so B = −∫(Dr_μ−Dr_μᵀ)φ_ξφ_ξ'dρ + 2κI. An antisymmetric matrix and I are linearly independent (take the trace), hence A ≡ 0 iff (a) Dr_μ = Dr_μᵀ ρ-a.e. [same Weierstrass argument applied to the matrix measure (Dr_μ−Dr_μᵀ)ρ] and (b) T_𝒪 = 0. In particular A = 2K when r_μ ≡ 0 on a neighbourhood of 𝒪." Propagate to lines 89, 1606, 1763, 1815, 1823–1827, 2157; add a test with r_μ = 0 on the orbit but curl r_μ ≠ 0 (the instance above). |
| I-2 | minor | "Reciprocal at order k_e⁻², at a frozen μ" (1753), "the set of reciprocal measures" (1768), "reciprocity … is equivalent to A_even + 2K ≡ 0" (2156) are undefined/overstated: (i) proves only first variation ⇒ A ≡ 0, and H.7(e)'s hypothesis is on a neighbourhood (the triangle), so the necessary condition is A_ν ≡ 0 for all ν near μ. The "if" direction at a single μ is not a reciprocity statement about the dynamics. | `scriptH_n1_perfectfit_triangle.py`: n = 1, f := f̂_μ (A_μ(ϑ₀,ϑ₁) = −1e-16); centroid measure μ_c = (1−ε)μ + (ε/3)(δ₀+δ₁+δ_z), ε = 0.3: A_{μ_c}(ϑ₀,ϑ₁) = −0.0612 = ε·A_ν; exact k_e²·circulation = 2.754e-3 at k_e = 80 and 160 vs predicted −(ε²/2)A_{μ_c} = 2.753e-3 (ratios 1.0002, 1.0003). So H.7(e)'s contradiction persists around every perfect-fit μ. | Define: "μ is *reciprocal at order k_e⁻²* if A_μ ≡ 0 (the 1-form P is closed at μ at that order)". In (i) say "only if A_ν ≡ 0 for every ν in a neighbourhood of μ, in particular A_μ ≡ 0". Add to (iii): "The reciprocal set {r_μ = 0 on 𝒪} has empty interior (εδ_ϑ changes r_μ by −εaφ_ξ), so for n = 1 the averaged loop is not a first variation up to O(k_e⁻³) on any open set of measures, perfect fit or not; A_μ ≡ 0 at perfect fit removes only the pointwise obstruction." Rephrase 2156 as "is a first variation up to O(k_e⁻³) only if A_even + 2K ≡ 0". |
| I-3 | minor | Remark I.1 bullet 2 (1826): "the performative defect is proportional to the residual, so it is small at a well-fitting equilibrium" — true for n = 1 only; for n ≥ 2 the curl term is not controlled by the residual on the orbit. | Same instance as I-1: ∫|r_μ|dρ = 8e-17, ∫‖Dr_μ−Dr_μᵀ‖dρ = 1.59, A_even entries up to 0.07. | State the general bound (checked, `scriptG_bound.py`, max ratio 0.083): \|A_even(ϑ,ϑ')\| ≤ \|a\|\|a'\|[∫‖Dr_μ−Dr_μᵀ‖dρ + (2e^{−1/2}/s_min)∫\|r_μ\|dρ]; "small at a well-fitting equilibrium in C¹ on the orbit". |
| I-4 | minor | Remark I.1 bullet 3 (1830): "the correction is O(k_e⁻²) relative to the leading drift" contradicts Remark H.0 ("O(1/k_e) relative to the drift"): the leading drift −(1/k_e)V_L is itself proportional to the residual, so the ratio is O(1/k_e). | P = (1/k_e)∫r_μᵀΦ̃dρ + (1/k_e²)(U_even+U_odd); both terms ∝ residual at μ*_τ. | "O(k_e⁻²) in absolute terms, O(1/k_e) relative to the leading drift". |
| I-5 | minor | (v) statement (1766–1769) omits "connected": the either/or dichotomy for a real-analytic zero set needs a connected parameter manifold (the proof, 1809, says it; the statement does not). "Real-analytic subvariety" for an intersection of uncountably many zero sets is fine locally (Noetherian) but "closed null set" suffices. | — | Insert "connected" in the statement; optionally replace "real-analytic subvariety" by "closed set of Lebesgue measure zero". |
| I-6 | minor | Line 99 "Theorem I settles when it is non-reciprocal" overclaims: for general n ≥ 2 only the dichotomy (v) and the (corrected) perfect-fit case are proved. | — | "settles it for n = 1, at perfect fit for n ≥ 2, and generically otherwise". |
| I-7 | cosmetic | (iii) and row I (89): "non-reciprocity is proportional to the residual" is an upper bound only; A can be arbitrarily small relative to ∫\|r_μ\|dρ. | `scriptE_n1.py`: residual 0.5 sin(ωx), s ∈ [0.3,0.6]: sup\|B\|/∫\|r\|dρ = 0.57, 0.12, 0.058, 0.022 for ω = 2, 8, 20, 40; random μ: ≤ 0.74 vs bound 4.04. | "bounded by the residual". |
| I-8 | cosmetic | H.7(e) step 5 (1653–1655) concludes cyc = 0 at the *centroid* μ_c(ε); passing to μ needs "let ε → 0 and use that A is affine in μ" (one clause). The proof of (i) (1775) inherits this. | — | Add the clause. |
| I-9 | cosmetic | *Checked* (1816) "F₂((0.5,0)) > 0": the subscript reads as width 2, but the test uses s = 0.7 and asserts the second *component*. | `tests/test_nonreciprocity.py:120–121`. | "the y-component of F_{0.7}((0.5,0)) is positive". |
| I-10 | cosmetic | (i) relies on ϑ_z = (α = 0, ξ) being in the parameter space; true by design (α ∈ ℝ unconstrained) but worth one clause, since with a lower bound on \|a\| the pairwise criterion would not follow. | — | "(the zero feature is admissible because α ranges over ℝ)". |

Everything else attacked holds; details below.

## Details

### (a)(i) Pairwise criterion — holds

- Algebra: A_even(ϑ,ϑ_z): both terms carry Φ̃_{ϑ_z} ≡ 0 (first as the direction (Φ̃_{ϑ_z}·∇), second inside
  r_μᵀΦ̃_{ϑ_z}); K(ϑ,ϑ_z) = ⟨Φ̃_ϑᵀ(d/dt)0⟩ = 0; so A(·,ϑ_z) = A(ϑ_z,·) = 0, and a coboundary a(ϑ) − a(ϑ_z) = 0 forces a ≡ const.
  The converse (A ≡ 0 is the coboundary of a ≡ 0) is indeed trivial.
- The triangle argument of H.7(e) applies verbatim with ϑ₂ := ϑ_z: μ_{a,b} = (1−ε)μ + (ε−a−b)δ_{ϑ₀} + aδ_{ϑ₁} + bδ_{ϑ_z}
  are probability measures; P_{μ_{a,b}} depends on b only through the mass removed from ϑ₀. Re-deriving dω for the
  1-form ω = ⟨P, dμ⟩ gives dω = −[A(ϑ₀,ϑ₁)+A(ϑ₁,ϑ₂)+A(ϑ₂,ϑ₀)] in general, and with ϑ₂ = ϑ_z the exact second variation
  has S(ϑ_z,·) = 0 (P(ϑ_z) ≡ 0) and S(·,ϑ_z) = 0 (mass at ϑ_z leaves r_μ unchanged), so dω = −A(ϑ₀,ϑ₁).
- Numerical confirmation at the *exact* P (`scriptF_triangle.py`, n = 2 ellipse, ε = 0.05, 6-point Gauss on each
  edge): predicted k_e²·circulation = −(ε²/2)A_{μ_c}(ϑ₀,ϑ₁) = −9.071e-5; exact k_e²·circulation −8.942e-5 (k_e = 80)
  and −9.005e-5 (k_e = 160), ratios 0.986 and 0.993, converging like 1/k_e. Exact S(·,ϑ_z) and S(ϑ_z,·) are 0.0 to the
  last bit (`scriptD_adversarial.py`).
- Gap I-8 (cosmetic) and definitional issue I-2 (minor) above.

### (a)(ii) Structural identity and B — hold

Differentiating the definition: (Φ̃'·∇)(rᵀΦ̃) = Σ_{k,l}Φ̃'_l(∂_l r_k)Φ̃_k + Σ_{k,l}Φ̃'_l r_k ∂_lΦ̃_k = Φ̃ᵀ(Dr)Φ̃' + rᵀ(Φ̃'·∇)Φ̃
with (Dr)_{kl} = ∂_l r_k; swapping gives Φ̃ᵀ(Dr)ᵀΦ̃' + rᵀ(Φ̃·∇)Φ̃'; the difference is Φ̃ᵀ(Dr−Drᵀ)Φ̃' + rᵀ[Φ̃',Φ̃] with
[X,Y] = (X·∇)Y − (Y·∇)X, so A_even = −∫{…}dρ with the stated sign. For Φ̃ = aφ: (Φ̃'·∇)Φ̃ = aφ'(a'ᵀ∇φ), giving
−∫r(∇φ)ᵀφ'dρ and +∫(∇φ')rᵀφdρ; the antisymmetric part of the k_e⁻² second variation of U_odd is K(ϑ,ϑ') − K(ϑ',ϑ) = 2K
= 2(aᵀa')κ, hence +2κI. The code (`B_matrix`, lines 78–88) implements exactly these four terms; `second_variation_formula`
and `expansion_terms` implement S_even, K, U_even, U_odd correctly (index check of every einsum). Agreement with the generic
code at adversarial points (a = ±2.985, c = ±1.9995 outside the orbit, s = 0.300 and 1.5, weights up to 0.6 and −0.3):
6.2e-15 (n = 1), 2.2e-15 (n = 2). Agreement with the exact second variation (Richardson 80/160): 8.8e-4 (n = 1
adversarial), 1.6e-4 (n = 2, I-1 instance); the n = 2 adversarial set has features ~1e-6 on the orbit (c = (2,±2),
s = 0.3), so the 1e-2 relative figure there is finite-difference noise on a vanishing quantity, not a discrepancy.
Also "a ranges over the open cube (−a_max,a_max)ⁿ ⇒ A ≡ 0 iff B ≡ 0" is correct.

### (a)(iii) n = 1 characterisation — holds

- φ_cφ_c' = e^{−(c−c')²/(4s²)}e^{−(x−m)²/s²} (checked: (x−c)²+(x−c')² = 2(x−m)² + (c−c')²/2) and
  φ_c'Dφ_c − φ_cDφ_c' = φ_cφ_c'(c−c')/s² are right; K ≡ 0 for n = 1.
- The Weierstrass-transform step is complete: ν = r_μρ is a finite signed measure (compact support is not even needed
  for m ↦ ∫e^{−(x−m)²/s²}dν to be entire: |e^{−(x−m)²/s²}| ≤ e^{(Im m)²/s²}); vanishing for m ∈ (−c_max,c_max) (c = m∓δ
  admissible for small δ) gives vanishing on ℝ; G_s∗ν ∈ L¹ with transform Ĝ_s·ν̂, Ĝ_s has no zeros, so ν = 0; r_μ
  continuous and ρ-a.e. zero means zero on supp ρ = the orbit image. The width restriction is harmless (one s suffices);
  the converse covers unequal widths since it is stated for ν = 0.
- Bound: sup|Dφ_ξ| = e^{−1/2}/s ≤ e^{−1/2}/s_min, so |A| ≤ 2e^{−1/2}/s_min·|a||a'|∫|r_μ|dρ (constant 4.04 at
  s_min = 0.3). Measured max of |B|/∫|r|dρ over 6 random μ × 400 pairs: 0.74 (`scriptE_n1.py`). No μ with
  ∫|r|dρ > 0 and B ≡ 0 was found (sup|B| ≥ 0.34 in every instance), as the proof requires. Only the wording "proportional"
  (I-7) and the definitional gloss (I-2) need attention.

### (a)(iv) n ≥ 2 perfect fit — first claim false (I-1), second claim holds

- Closed form re-derived: with equal widths ∇φ_ξ' = −((x−c')/s²)φ_ξ', so κ = −(1/s²)e^{−|c−c'|²/(4s²)}
  ⟨e^{−|x_m−m|²/s²}(x_m−c')·ẋ_m⟩; x_m − c' = (x_m−m) + (c−c')/2; the first part is (d/dt)[−(s²/2)e^{−|x_m−m|²/s²}]
  and integrates to zero over a period; hence κ = −(1/(2s²))e^{−|c−c'|²/(4s²)}(c−c')·F_s(m). Confirmed to 1e-10 by the
  test and on a figure-eight (`scriptC_fig8.py`).
- Injectivity: F_s = (1/T₀)G_s∗T_𝒪 with T_𝒪 a compactly supported finite ℝⁿ-valued measure; κ ≡ 0 on admissible pairs
  gives v·F_s(m) = 0 for all directions v and all ‖m‖ < c_max (c = m∓δv in the cube), F_s entire ⇒ F_s ≡ 0 ⇒
  Ĝ_s·T̂_𝒪 = 0 ⇒ T_𝒪 = 0. Converse immediate. For n = 2: ∮v·dx = ∫w_𝒪 curl v dA holds for every closed C¹ curve
  (∂(w_𝒪[ℝ²]) = T_𝒪), every g ∈ C_c^∞ is curl of the smooth v = (0, ∫_{−∞}^{x₁}g(t,x₂)dt) — v need not be compactly
  supported since T_𝒪 is — so T_𝒪 = 0 ⇔ w_𝒪 = 0 a.e. (w_𝒪 ∈ L¹, orbit image Lebesgue-null). Legitimate.
- Numerical attacks found nothing against the second claim: 300 random equal-width pairs on the ellipse give
  max|κ| = 0.097; a segment traversed back and forth and a circle traversed forwards then backwards (T_𝒪 = 0) give
  |K| ≤ 1e-16; a figure-eight (sin t, ½sin 2t) has signed area exactly 0 but max|K| = 0.049 and F_{0.5}(±0.5,0) =
  (0, 0.050) — consistent with the text's careful use of "circulation"/"winding function" rather than "area".
- The first claim: see I-1. The error is the carry-over of the n = 1 reasoning (where Dr − Drᵀ ≡ 0) to n ≥ 2.

### (a)(v) Real-analyticity and dichotomy — hold (I-5 wording)

tanh(α/a_max), tanh(γ/c_max) and sig(κ) are holomorphic on the strips |Im α| < πa_max/2, |Im γ| < πc_max/2, |Im κ| < π;
in |Im κ| < π, s = s_min + (s_max−s_min)sig(κ) never vanishes (zeros would need e^{−κ} real negative, i.e. Im κ = ±π);
for complex c, |exp(−Σ(x_k−c_k)²/(2s²))| ≤ exp(C|Im c|²) uniformly for x on the compact orbit. So a uniform thin
complex neighbourhood of ℝ^p exists on which every integrand of B is holomorphic with uniform bounds; integration
against ρ (any Borel probability measure on the orbit) and against ẋ_m dt preserves holomorphy. The zero set of a
real-analytic function on a connected manifold is the whole manifold or Lebesgue-null; the reciprocal set is the
intersection over (ξ,ξ') of such sets, hence closed and either everything or null. Connectedness must be in the statement.

### (b) Counterexample search — summary

(i) n = 1: none (impossible by (iii); confirmed on 6 random measures × 400 geometry pairs). (ii) n = 2: no orbit with
T_𝒪 ≠ 0 and K ≡ 0, no zero-current orbit with K ≠ 0; the figure-eight shows zero signed area ≠ zero circulation.
(iii) B vs generic formula: 6e-15 at adversarial points; B vs exact second variation: 1.6e-4–8.8e-4 relative where A is
O(1); exact S(·,ϑ_z) = S(ϑ_z,·) = 0 exactly. **The counterexample that does exist is to (iv)'s perfect-fit lemma (I-1).**

### (c) Hypotheses and scope

- Regularity: (i) needs the O(k_e⁻³) remainder of H.7(b) uniformly in the triangle measures, hence H.7's f ∈ C²,
  x_m ∈ C² (and k₀ independent of μ, which H.7 states); (ii)–(v) use only f ∈ C¹, x_m ∈ C¹. Inheriting H.7's hypotheses
  is legitimate; a one-line remark would help.
- ρ: no density is assumed anywhere; ν = r_μρ and T_𝒪 are treated as measures, and the Fourier argument is for finite
  measures. Correct.
- "First variation up to O(k_e⁻³) in the sense of H.7(e)": used consistently, but that sense is a hypothesis on a
  neighbourhood of μ with constants uniform in k_e (ℱ and c may depend on k_e); the theorem should say so (I-2).
- Overclaims: I-1 (perfect fit, n ≥ 2), I-2 (equivalence at 2156), I-3, I-4, I-6. The text does *not* claim
  non-reciprocity "for all μ" in n ≥ 2: line 1825 ("fails only on a null set of any analytic family … that contains one
  non-reciprocal member") and H.7(e)'s "whenever these conditions fail" are correctly hedged.
- Numbers (`scriptA_cyclic.py`, the exact `_setup`/`_plant_*` of `test_averaging.py`, N = 2000, sub = 2, TE[:3]):
  predicted cyclic sums −0.07035 (n = 1; A_even only, K part 0) and −0.04637 (n = 2; A_even −0.08935, 2K +0.04298);
  exact k_e²·cyc at k_e = 80: −0.07107 and −0.04706. All quoted figures (−0.0704, −0.0464, −0.0711, −0.0471, −0.089,
  +0.043) are reproduced.

### Code notes (no bugs found)

`nonreciprocity.py` and `averaging.py` are consistent with the formulas; `kappa`, `current_F`, `kappa_closed_form`,
`A_even_1d_closed_form`, `positivity_witness_1d` are correct. Test coverage gap: the perfect-fit tests
(`test_n1_reciprocal_iff_residual_vanishes_on_the_orbit`, `test_perfect_fit_n2_leaves_only_the_circulation_kernel`)
use f := f̂_μ (r_μ ≡ 0 on ℝⁿ), which is strictly stronger than the theorem's hypothesis r_μ = 0 ρ-a.e. and hides I-1;
add the f := f̂_μ + h e₁ instance as a test of the corrected (iv).
