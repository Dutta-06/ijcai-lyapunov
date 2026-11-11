# Adversarial review of Lemma E.0 and Theorem E (with explicit constants)

Reviewed: `PARADIGM.md` §6 (working copy; §6 is byte-identical to HEAD `ca69ced` apart from a line shift — current
anchors: §6 l.1086, Lemma E.0 l.1110, proof l.1139, Remark E.0a l.1173, Theorem E l.1198, proof steps 3/4/5 at
l.1244/1259/1273, Remark E.1 l.1284, Remark E.3 l.1324), `src/core/excitation.py`, `tests/test_excitation.py`,
`experiments/exp_explicit_rate.py`, `results/explicit_rate/table.md`. Baseline `pytest tests/test_excitation.py`:
5 passed (222 s). While this review ran, uncommitted edits appeared in the working copy of `excitation.py` (module
docstring), `tests/test_excitation.py` (docstrings, `n_nodes=800`) and `exp_explicit_rate.py` (`n_nodes=800`, a
closed-form column, reworded saturation text); where an issue is already addressed by those edits I say so. No file in
the repository was modified by this review. Scratch scripts and raw outputs (same directory as this file):
`lemma_attack.py/.out`, `constants_attack.py/.out`, `bandwidth_check.py/.out`, `cancellation_check.py`,
`fd_and_degenerate.py`.

## Verdict

**Lemma E.0 holds.** Every inequality in the proof is correct as written: Q̇ = ∫_t^{t+T}M − TM(t) (Leibniz; the
upper-limit term vanishes because the weight t+T−r is zero there), the three derivatives in (c), the exact
cancellation of the ±εT‖Ψζ‖² terms, the four cross-term bounds summing to εc₁ (2Q_MΨ_M = T²Ψ_M³), the Young step,
the 2×2 form (valid for all real (x, y); for the closed-form ε the matrix is ⪰ diag(k/4, εα/2), so λ_min ≥ γ_cf as the
test asserts), p_lo/p_hi, the envelope exponent γ/(2p_hi), and √3·e^{−2γt/3}. Numerically I checked the lemma in its
strongest pointwise form: V(t,·) = ½wᵀH(t)w and V̇(t,·) = wᵀN(t)w are quadratic forms, so I verified λ_max(N(t)) ≤ −γ
and p_lo ≤ ½λ(H(t)) ≤ p_hi over **all** w, with my own composite quadrature and Hessian (independent of the module), at
40 random times per case, for ten regressor families (q > n, q < n, n = q = 1; k from 10⁻³ to 10⁴; Ψ_D/Ψ_M = 100 with a
window of 50 oscillations; quasi-periodic non-periodic PE on [0, 300]; amplitude 30; near-degenerate PE with
α ≈ 5·10⁻⁶; slow adaptation with T = 20π) and eight ε per case (closed form, optimised, 10⁻⁶…0.999 of the admissible
maximum). No violation: for every ε with γ > 0 the worst scale-free margin (λ_max(N)+γ)/γ was −4.5·10⁻⁶ (at
ε = 10⁻⁶ε_max; closed form ≤ −0.30, optimised ≤ −0.81), and the V bounds never failed, even at p_lo = 5·10⁻⁴. The
module's V, V̇ agree with the independent forms to 10⁻¹⁵–10⁻¹² (except when its fixed 200-node quadrature is
under-resolved, E0-4) and its chain-rule V̇ agrees with finite differences along DOP853 trajectories to 10⁻⁵ (the FD
truncation level). **Theorem E with explicit constants holds**, given a true Lipschitz constant L_J of the Jacobian
on D (as the statement says): A(t), the similarity S A S⁻¹ = [[−k_eI, −Ψ_G],[Ψ_Gᵀ, 0]], the PE-data transfer, the
remainder bound C̃_R, the ball r_w, the forward invariance of {V ≤ p_lo r_w²} ⊃ {‖w‖ ≤ r₀}, ‖w‖² = 2V_π,
C = (2p_hi/p_lo)^{1/2}, λ = γ/(4p_hi) ≥ γ/3 ≥ (1/3)min{k_e/4, εα_G/2}, and the uniqueness argument are all correct; the
regressor derivative, the Floquet code and the ω⁻² asymptotics of (iii) check out numerically, and c₀ < d²/8 comes
out automatically on the m = 3 ideal (by 27 orders of magnitude). The genuine problems are not in the mathematics but
in the numerics behind the tables and in Remark E.0a's quantitative narrative: (1) the 2×2 λ_min in
`strict_lyapunov_constants` is computed with catastrophic cancellation, so the tabulated explicit rates in the
slow-adaptation regime carry 1 %–130 % round-off error (table 2's g = 30 row and its "8.4·10¹⁰" are pure round-off:
the module returns exactly 2⁻⁵⁰), and γ underflows to 0 for ill-conditioned PE; (2) `optimise_eps` returns `None` in
that regime; (3) the saturation accounting of Remark E.0a attributes to "the optimised ε" a factor that is in fact the
product of four unrelated effects, and rests on a PE level that is 21 % under-resolved and inconsistent with the κ
column. Nothing is fatal; two items are major because quoted numbers/explanations are wrong.

## Issues

| id | severity | issue | evidence | proposed fix |
|---|---|---|---|---|
| E0-1 | **major** (code → table 2, Remark E.0a numbers) | `gamma = 0.5(a+b) − sqrt(0.25(a−b)² + 0.25c²)` (`excitation.py` l.70) is numerically unstable when b, c ≪ a (every slow-adaptation row: a ≈ k_e = 6, γ ~ 10⁻¹⁴–10⁻¹⁶, floor a·ε_mach = 1.3·10⁻¹⁵). The returned γ/`rate` is garbage near the floor: table 2's g = 30 row reports 8.882e-16 = exactly 2⁻⁵⁰ while the exact λ_min is 3.81·10⁻¹⁶ (error +133 %), so the quoted "8.4·10¹⁰ at g = 30" (PARADIGM l.1181, table.md) should be 1.96·10¹¹; errors are 0.14 % (g = 1), 0.7 % (g = 0.01), 2.5 % (0.003), 20 % (0.001), 160 % (10⁻⁴). For ill-conditioned PE γ underflows to exactly 0 and `rate = 0` although γ_cf = εα/2 > 0. Since the computed "explicit rate" can exceed the one the lemma guarantees, "rigorous lower bound in every row" is not literally true (the ordering explicit ≤ exact survives only because the margin is ≥ 10⁵) | `cancellation_check.py` (60-digit Decimal reference) and `fd_and_degenerate.py`: near-degenerate case γ_cf = 3.2·10⁻¹⁹, exact λ_min = 4.8·10⁻¹⁹, module γ = 0, rate = 0 | compute λ_min as det/λ_max: `lmax = 0.5(a+b) + sqrt(0.25(a−b)² + 0.25c²); gamma = (a·b − 0.25c²)/lmax` (relative error 10⁻¹⁶ in all rows tested); regenerate table 2 and the numbers quoted from it; also return `max(gamma, gamma_cf)` for the closed-form ε |
| E0-2 | **major** (Remark E.0a text + table 2; no theorem depends on it) | the slow-adaptation accounting mis-attributes the discrepancy: (a) 3κ(1+Ψ_D/(kΨ_M))² is the leading-order value of exact/(2γ_cf/3), but table 2 tabulates the optimised-ε rate, asymptotically kα²/c₁² = 3× the closed-form value (the module's "closed-form rate" γ_2×2/(2p_hi) is 2.25×, and the optimised one 4/3 more: 3 in all — verified: opt-rate = kα²/c₁² to 0.3 %); (b) "Here the exact rate is λ_min⟨ΨᵀΨ⟩/k" (l.1178) is off by a g-independent factor 0.658 on this instance (the two-tone reference is not slow relative to k_e = 6): filtered averaging λ_min(sym⟨Ψᵀ(h∗Ψ)⟩), h(s) = e^{−k_e s}, reproduces the Floquet rate to 5 %; (c) table 2's α (200 GL nodes over 20π, min over 80 t's) is 21 % below the converged value (1.07e-3 vs 1.3596e-3 at g = 1; 800/1600 nodes and a 200×20 composite rule agree), whereas κ was computed with 400 nodes: the inconsistency enters the ratio as 1.27² = 1.61; (d) at g = 0.01 the T²Ψ_M³ term is still 13 % of c₁ (c₁/T = 2.235 vs kΨ_M+Ψ_D = 1.949; TΨ_M² = 1.73 vs k = 6 — the regime condition is TΨ_M² ≪ k, not Ψ_M² ≪ k), (c₁/c₁^∞)² = 1.31. Product (1/3)·0.658·1.61·1.31 = 1/2.15: this *is* the "remaining factor of about 2", and it is not "the optimised ε". The ratio has not saturated at g = 0.01 (already below the quoted 1.47·10⁶ and still falling); with the table's α convention it tends to ≈ 5.2·10⁵ (rows: 6.77e5, 5.73e5 at g = 0.003, consistent with 5.2e5·(c₁/c₁^∞)²), with a converged α to ≈ 0.66·κ(1+Ψ_D/(kΨ_M))² ≈ 3.2·10⁵ | `constants_attack.out` [5] (α convergence, exact/avg = 0.660/0.659/0.658 at g = 1/0.1/0.01, c₁/T vs kΨ_M+Ψ_D per row), `bandwidth_check.out` (exact/naive 0.658, exact/filtered 0.953) | state the formula for the rate actually tabulated (κ(1+Ψ_D/(kΨ_M))² for the optimised ε, or tabulate 2γ_cf/3 — the working copy adds a closed-form column); replace "the exact rate is λ_min/k" by the filtered-averaging value or quote the bandwidth factor; use a converged α (the working copy's `n_nodes=800` does this; the measured ratios then drop by ≈ 1.6 and every number quoted in l.1181–1183 must be regenerated); write the regime as TΨ_M² ≪ k; add rows g = 10⁻³, 10⁻⁴ (after E0-1 is fixed, otherwise they are round-off); delete "the optimised ε accounts for the remaining factor of about 2" |
| E0-3 | minor (code) | `optimise_eps` returns `None` when every ε with γ > 0 lies below its grid floor `eps_V·1e-9`; all callers index the result and would crash | near-degenerate regressor (n = q = 2, two almost parallel columns), k = 1, T = 2π, Ψ_M = 6.16, Ψ_D = 10.0, α = 5.27·10⁻⁶: closed-form ε = 6.0·10⁻¹⁴ < grid floor 6.5·10⁻¹³, all grid points have γ ≤ 0 (`lemma_attack.out` l.136–137; reproduced in `fd_and_degenerate.py`) | add the closed-form ε to the candidates, start the grid at ≤ 10⁻³·ε_cf, fall back to the closed-form constants; never return `None` |
| E0-4 | minor (numerics) | the module's fixed 200-node Gauss–Legendre rule on [t, t+T] (`gramian`, `Q_matrix`, `pe_level`, `StrictLyapunov`) is under-resolved when the window holds many feature widths; errors of unknown sign of O(10 %) then enter α, V and V̇ | table 2: 21 % in α (E0-2c); case "fast (2,3), frequencies ×50, T = 2π": module V̇ vs exact quadratic form differs by 4·10⁻⁶ (closed-form ε) up to 52 % (ε = 0.9ε_max) (`lemma_attack.out` l.75–89) | composite rule with panels ∝ T·(max frequency), or a 2×-refinement convergence check in `pe_level`/`Q_matrix` that raises on failure; document the assumption |
| E0-5 | minor (overclaim) | "The explicit rate is a rigorous lower bound of the exact rate in every row" (`table.md` l.3, `exp_explicit_rate.py` l.104), "the maximum is a rigorous bound" (`excitation.py` l.93, PARADIGM l.1169): the inputs are grid maxima for Ψ_M, Ψ_D (80–100 samples/period), a grid minimum of an under-resolved quadrature for α, a 300-point ε grid, and a γ with round-off (E0-1); the ordering is safe only because the margins are ≥ 17× | `constants_attack.out` [6]: 100-point vs 20000-point sup norms: Ψ_D 3.689 vs 3.690 (ω = 1), 236.1 vs 236.2 (ω = 64) — harmless here, but not certified | write "numerical lower bound (grid-estimated data)", or certify Ψ_M, Ψ_D by Lipschitz-corrected grid maxima and α by a converged quadrature with an error bound |
| E0-6 | minor (text) | Lemma E.0(iii) and step 3 assert that (LTV) with Ψ = sin t² "is not uniformly exponentially stable" on the strength of three simulated ratios (PARADIGM l.1133–1134, l.1254–1256) | the quoted 0.9715/0.9942/0.9980 equal 1 − (1/8)(1/t₀ − 1/(t₀+200)) to three digits, i.e. the local ω⁻² law | add the two-line proof: with F(s) := ∫_s^∞ sin r² dr = O(1/s) and e(t₀) = 0, integrating by parts in e(t) = −∫_{t₀}^t e^{−k(t−s)} sin(s²) ζ(s) ds (using |ζ̇| ≤ |e|) gives |e(t)| ≤ (C/t₀)‖w(t₀)‖ for all t ≥ t₀, so E(t₀+Δ)/E(t₀) ≥ 1 − 2kC²Δ/t₀² → 1; or write "numerically" |
| E0-7 | cosmetic | g_V := 1 + 2ε(Q_M + TΨ_M) is valid but loose: ∇_wV = H(t)w with ‖H‖ ≤ 1 + ε(TΨ_M + 2Q_M); the statement of Theorem E writes Q_M, Ψ_M where Ψ_{G,M} is meant | PARADIGM l.1201–1202; `excitation.py` l.80; `constants_attack.out` [8]: sup‖H‖ = 1.0017, code g_V = 1.505, sharper bound 1.5 (both far from sharp; the loss in r_w is at most the factor 1.5/1.0017) | use 1 + ε(TΨ_{G,M} + 2Q_M), or the exact sup‖H(t)‖; write Ψ_{G,M} |
| E0-8 | cosmetic | step 4: r_w = r/‖S⁻¹‖ allows ‖w‖ = r_w, i.e. ‖z‖ = r on ∂D while D is the open ball | PARADIGM l.1259–1262 | take D closed (the Lipschitz bound extends to the closure) or use a strict inequality |
| E0-9 | cosmetic (docs) | the module docstring at HEAD states a different closed form, ε := min{1/(4(Q_M+TΨ_M/2)), k/(2(TΨ_M²+c₁²/(2α)))}, γ = min{k/2, εα/2} (also valid, but neither the lemma's nor the code's) | `excitation.py` l.31–35 at HEAD | already fixed in the uncommitted working copy |
| E0-10 | cosmetic (tests) | "the bound is sharp up to a factor" (the factor is 10⁵–10⁶) and `ratios[-1] < 1.5*ratios[-2]` (observed 0.21) do not test saturation | `tests/test_excitation.py` l.114, l.135 | already relabelled a "no-growth check" in the working copy; a saturation test needs g = 10⁻³, 10⁻⁴ and the stable γ of E0-1 |
| E0-11 | cosmetic | (iii) implicitly needs Ψ_D > 0; "c₁ → TΨ_D·(1 + o(1))" should read c₁ = TΨ_D(1 + O(1/ω)); the o(1) in the rate is O(1/ω) too | PARADIGM l.1135–1136; `constants_attack.out` [4]: (2γ_cf/3)/(kα²/(3ω²T²Ψ_D²)) = 0.650, 0.970, 0.997, 0.9997 at ω = 10, 10², 10³, 10⁴ | reword |

## Details

### (a) Proof of Lemma E.0 — checked line by line

- *(a) Gramian weight.* Q(t) = ∫_0^T (T−s)M(t+s)ds; by Leibniz on ∫_t^{t+T}(t+T−r)M(r)dr the upper limit contributes
  (t+T−(t+T))M = 0, the lower −T·M(t), the t-derivative of the integrand ∫M. Correct, and only continuity of M is used.
  0 ⪯ Q ⪯ Ψ_M²∫_0^T s ds = Q_M. Correct.
- *(b)* Correct; p_hi is slightly conservative since −εζᵀQζ ≤ 0.
- *(c)* d/dt ½‖w‖² = −k‖e‖² (the eᵀΨζ terms cancel exactly); d/dt(−ζᵀQζ) = −ζᵀ(∫M)ζ + T‖Ψζ‖² − 2ζᵀQΨᵀe;
  d/dt(eᵀΨζ) = −k eᵀΨζ − ‖Ψζ‖² + eᵀΨ̇ζ + ‖Ψᵀe‖². With weights 1, ε, εT the ‖Ψζ‖² terms cancel exactly (this is
  what fixes the factor T in the cross term). Cross terms: 2εQ_MΨ_M + εTkΨ_M + εTΨ_D = εc₁. Correct.
- *2×2 form.* −ax²−by²+cxy ≤ −λ_min[[a,−c/2],[−c/2,b]](x²+y²) for all real (x, y). For the closed-form ε Young gives
  the matrix ⪰ diag(k/4, εα/2), hence γ_2×2 ≥ γ_cf, which `test_lemma_E0_pointwise…` asserts. Correct (the *formula*
  used by the code for λ_min is the unstable one, E0-1).
- *(d)* ε ≤ kα/c₁² ⇒ εc₁²/(2α) ≤ k/2; ε ≤ k/(4TΨ_M²) ⇒ a ≥ 3k/4; V ≤ ¾‖w‖² ⇒ V̇ ≤ −(4γ/3)V; ‖w‖² ≤ 4V ⇒
  ‖w(t)‖² ≤ 3e^{−4γ(t−t₀)/3}‖w(t₀)‖². Correct.
- *(e)* ∫_t^{t+T/ω}Ψ_ωᵀΨ_ω = ω⁻¹∫_{ωt}^{ωt+T}ΨᵀΨ ⪰ (α/ω)I; c₁(ω) = T²Ψ_M³/ω² + TkΨ_M/ω + TΨ_D; for large ω the
  active term of ε is kα_ω/c₁² and γ_cf = εα_ω/2, giving 2γ_cf/3 = kα²/(3ω²T²Ψ_D²)(1+O(1/ω)). Verified to 3·10⁻⁴ at
  ω = 10⁴ ([4]); the module's `rate` is asymptotically 9/4 of this (5.37e-11 vs 2.39e-11), as derived under E0-2(a).
- *λ ≥ (1/3)min{k_e/4, εα_G/2}.* λ = γ/(4p_hi) with p_hi ≤ 3/4 and γ ≥ γ_cf: consistent.

### (a) Theorem E, steps 3–5

- Step 2/3: finite-difference Jacobian of the implemented field at z = 0 vs A(t) = [[−k_eI, −Ψ/m],[(I⊗G)Ψᵀ, 0]]:
  rel. err 4.8·10⁻¹²; S A S⁻¹ vs `ltv_matrix(Ψ_G, k_e)`: 1.7·10⁻¹⁷ ([2]). PE transfer:
  (I⊗G^{1/2})(∫ΨᵀΨ)(I⊗G^{1/2})/m ⪰ (α/m)(I⊗G) ⪰ (αλ_min(G)/m)I; the norm bounds by submultiplicativity. Correct.
- Step 4: Taylor with a Lipschitz Jacobian on the convex D gives ‖R‖ ≤ (L_J/2)‖z‖²; ‖R̃(t,w)‖ = ‖S R(t, S⁻¹w)‖ ≤
  ‖S‖(L_J/2)‖S⁻¹‖²‖w‖²; the *pointwise* linear inequality V̇_lin(t,w) ≤ −γ‖w‖² plus ∇_wV·R̃ ≤ g_V C̃_R‖w‖³ gives
  V̇ ≤ −(γ/2)‖w‖² for ‖w‖ ≤ γ/(2g_VC̃_R); {V ≤ p_lo r_w²} ⊂ {‖w‖ ≤ r_w} by the lower bound and ⊃ {‖w‖ ≤ r₀} by the
  upper bound (p_hi r₀² = p_lo r_w²); V̇ ≤ −(γ/(2p_hi))V gives λ = γ/(4p_hi). Correct (E0-8 is a boundary
  technicality). `theorem_E_constants` implements exactly this (‖S‖ = max(1, (mλ_min(G))^{−1/2}),
  ‖S⁻¹‖ = max(1, (mλ_G)^{1/2}), r0_z = r0_w/‖S‖ is a correct sufficient radius in z).
- Step 5: ‖ζ‖² = (1/m)Σ‖G^{−1/2}ϑ̃_i‖² = W_π², so ‖w‖² = 2V_π; optimality of π at t₀ gives V_π(t₀) = V_OT(t₀) ≤ c₀ =
  r₀²/2 and W(t₀) = W_π(t₀); ‖e‖+W ≤ ‖e‖+‖ζ‖ ≤ √2‖w‖ and ‖w(t₀)‖ ≤ ‖e(t₀)‖+W(t₀); hence C = √2·(p_hi/p_lo)^{1/2}.
  Uniqueness: V_π, V_π' ≤ c₀ put the state in both invariant balls (step 4 needs only ‖w‖ ≤ r₀, not optimality),
  forcing convergence to two distinct equilibria (the ϑ*_j are distinct by PE). Correct. Consistency with Remark E.3
  (c₀ < d²/8 must come out automatically): on the m = 3 ideal with sampled L_J, c₀ = 1.4·10⁻²⁸ vs d²/8 = 0.165
  (g = 0.3) and 2.2·10⁻³⁷ vs 4.9·10⁻³ (g = 10) ([7]); r_w is set by γ/(2g_VC̃_R), not by r/‖S⁻¹‖, by 13 orders of
  magnitude.
- The constants are explicit modulo L_J, a well-defined supremum that the tests only estimate from 150 sampled pairs
  (stated in the text and the docstring). The Theorem E test is therefore a consistency check of the envelope, not a
  certification of r₀; its horizon (8000 s ≈ 0.3 e-folds of the envelope, as the working-copy docstring now says) is
  honest.

### (b) Counterexample search — none found

`lemma_attack.py`: exact quadratic-form check over all w at 40 random times per (case, ε); own composite
Gauss–Legendre quadrature for ∫M and Q; rigorous sup-norm bounds for the sinusoidal regressors (Σ_j‖[A_j B_j]‖,
Σ_j ω_j‖[A_j B_j]‖); α = min over a 60-point grid of λ_min(∫_t^{t+T}M). Cases: periodic (2,3) k = 3; (1,5) k = 0.5
(q > n); (4,2) k = 100 (q < n); (1,1) k = 10⁻³; fast (2,3) k = 1 with frequencies ×50, T = period/50 (Ψ_D/Ψ_M = 100)
and T = 2π (50 oscillations per window); quasi-periodic (2,3) k = 2 (tones 1, √2, 0.35π; window 10; t ∈ [0, 300]);
(3,3) k = 10⁴, amplitude 0.05; (2,2) k = 1, amplitude 30 (α = 3963); near-degenerate PE (α = 5.3·10⁻⁶); slow (2,5)
k = 6, amplitude 0.02, T = 20π. ε ∈ {closed form, optimised, 10⁻⁶, 10⁻³, 0.03, 0.3, 0.9, 0.999 × ε_max}. For every ε
with γ > 0, max_t (λ_max(N(t))+γ)/γ ≤ −4.5·10⁻⁶ (closed form ≤ −0.30, optimised ≤ −0.81); for every ε, including
0.999·ε_max where p_lo = 5·10⁻⁴, p_lo ≤ ½λ_min(H) and ½λ_max(H) ≤ p_hi. For ε with γ ≤ 0 the lemma claims nothing
and V̇ is indeed not negative definite. Finite differences of V along DOP853 trajectories agree with
`StrictLyapunov.Vdot` to 4·10⁻⁶–2·10⁻⁵ (closed-form, optimised and 0.3·ε_max). The optimised ε never produced a
failing Lyapunov function; the only failures found are the code/numerics items E0-1, E0-3, E0-4.

### (c) Hypotheses, constants, numbers

- PE data of Ψ_G: correct (above). `TheoremERegressor`: Ψ_G equals m^{−1/2}[∇_ϑΦ(x_m;ϑ*_j)]_j(I⊗G^{1/2}) built
  from `transport.grad_Phi` to 4·10⁻¹⁶; `Psidot` vs central differences: 1.6·10⁻⁹ ([1]).
- Floquet code: RK4 monodromy, rate = −log max|λ_i|/period; equals a DOP853 (rtol 10⁻¹²) monodromy to all printed
  digits (4.660231e-02) and the fitted decay of a trajectory over 240 periods (4.6619e-02) ([3]).
- ω-rescaling: the data (α/ω, T/ω, Ψ_M, ωΨ_D) are reproduced on the rescaled regressor; the asymptotic formula holds
  with relative error 0.35, 0.03, 0.003, 0.0003 at ω = 10…10⁴ ([4]).
- Remark E.0a derivation: with T = period, α = Tλ_min, c₁ ≈ T(kΨ_M+Ψ_D), active term ε = kα/c₁², γ_cf = εα/2:
  2γ_cf/3 = kα²/(3c₁²) and (λ_min/k)/(kα²/(3c₁²)) = 3(kΨ_M+Ψ_D)²/(k²λ_min) = 3κ(1+Ψ_D/(kΨ_M))². The algebra is right;
  what the four approximations cost on the table's own instance is itemised under E0-2. For the optimised ε
  (maximise εα − ε²c₁²/(4k); optimum ε = 2kα/c₁²) the rate is ≈ kα²/c₁², so the limit of the tabulated ratio is
  κ(1+Ψ_D/(kΨ_M))² times the bandwidth factor (0.66 here).
- Text vs table: every number quoted in §6 (954 → 17.6; slopes −1.97/−1.56; −1.92 over the last octave;
  8.4·10¹⁰ → 6.8·10⁵; 1.47·10⁶; κ = 1.27·10⁵; 4.7e-2 → 2.9e-5; best window = half a period) agrees with `table.md`
  as committed; 8.4·10¹⁰ is round-off (E0-1) and the table-2 numbers change with a converged α (E0-2), so l.1181–1183
  must be regenerated after the fixes.
- Overclaims: E0-5, E0-6, E0-10; otherwise the text is careful ("valid but conservative", "estimate from samples").
