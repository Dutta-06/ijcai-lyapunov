# Transport adaptive control: adaptive approximators as evolving measures

**Version 2** (after adversarial review round 1 and the literature audit).

**Companion files.**

- **Code.** `src/core/transport.py` covers laws (L) and (L_σ) and the oracle Lyapunov report. `src/core/langevin.py` covers Theorem H.
- **Tests.** `tests/test_transport.py` and `tests/test_langevin.py`; run `venv\Scripts\python -m pytest tests`.
- **Experiments.** `experiments/exp_transport_*.py`, `exp_ot_sigma_mod.py`, `exp_sigma_mod_transient.py`, `exp_mean_field_m.py`, `exp_pe_convergence.py`, `exp_variable_projection.py` and `exp_langevin_global.py`. Their numbers are collected in `results/RESULTS_PARADIGM.md`.
- **Literature audit.** `LITERATURE_TRANSPORT.md` gives claim-by-claim novelty verdicts. Reference numbers [1]–[89] point to `review/LITERATURE_refs.md` and [S1]–[S9] to the audit file.
- **Review log.** `review/REVIEW_LOG.md` lists every confirmed issue and its fix. The raw reports are in `review/REVIEW_round1.md`.

---

## 0. What this is, what is new, and what is not

**The idea.**  An adaptive approximator with nonlinearly parametrised features is an
*empirical measure* μ = (1/m)Σδ_{ϑ_i} on a lifted parameter space. Here that means Gaussian RBF kernels with
adaptive weights, centres and widths. The adaptive law transports this measure. The control output
f̂ = ∫Φ dμ is linear in μ. Stability, robustness modifications, excitation, the many-kernel limit,
the slow evolution of kernel geometry and global optimality of kernel placement all acquire
formulations on Wasserstein space. These formulations are independent of how the kernels are labelled.

**Why the classical theory needs it.**  The audit of `lyapunov.pdf` (`ASSESSMENT.md`) found two structural
defects. First, the classical Lyapunov analyses are label dependent: a network that is a relabelling of
the ideal has zero function error and large parameter error. Second, there is no object describing where
the kernels go or whether their placement is good. A third apparent defect, that constants grow with the
number of kernels, turned out to be a normalisation artefact. The mean-field normalisation removes it for
labelled and label-free analyses alike (Remark A.1).

**Contributions and their honest status.**  "Proved" means proved here. "Checked" means the inequality
is evaluated at every sample of simulated trajectories by a test. Each novelty column follows the
verdict in `LITERATURE_TRANSPORT.md`.

| result | statement | status | novelty (audit) |
|---|---|---|---|
| A | the adaptive law is the particle form of the Wasserstein gradient flow of a μ-linear potential | proved and checked; raw and projected identification | interpretive; ingredients known [29], [41], [46], [25] |
| B | V_OT = ½‖e‖² + ½W² satisfies a Dini inequality with a two-regime matched remainder | proved and checked, including a moving ideal, ε ≠ 0, assignment switches and an active projection | packaging new; ingredients known [62], [64], [20], [49] |
| C | OT σ-modification: exact dissipation, gap-protected strict inequality, a computable local exponential estimate, and reshuffling by the labelled leak | proved and checked; the C(iv) instance satisfies every hypothesis | novel as a control construction; C(i) and the obstruction are known [71], [73] |
| D | well-posed mean-field closed loop, a Dobrushin estimate and a global tube bound | proved; consistency in m checked numerically | technique anticipated [77]–[79] |
| E | label-free local exponential convergence under tangential PE | proved; linearisation checked; convergence observed qualitatively | partially anticipated [49], [54] |
| F | Wasserstein length of learning ≤ λ_G^{1/2} L₁ ∫‖e‖ | proved and checked | elementary |
| G | the slow geometry performs variable-projection descent on the occupation-measure loss | formal three-scale reduction; numerical check | a bridge to [S1], [S2] |
| H | the averaged mean-field loop with saturated features, σ-modification and dither is mean-field Langevin dynamics; it converges globally and exponentially, with the log-Sobolev constant in control quantities | proved for the averaged system; the real closed loop is Conjecture H.6; numerics | a bridge to [S3], [S4], new in adaptive control |

The genuinely new mathematics required by this viewpoint lies in the open problems of §8. They are
infinite-horizon averaging for McKean–Vlasov closed loops, uniform-in-time propagation of chaos with a
plant in the loop, and global convergence without dither.

---

## 1. Setting

**Plant, reference, control, error.**

    ẋ = f(x,t) + u,    e = x − x_m(t),    u = ẋ_m − k_e e − f̂(x).

**Lifted parameter.**  ϑ = (a, c, s) ∈ ℝⁿ × ℝⁿ × ℝ_{>0}, with p := 2n+1. The feature is
Φ(x; ϑ) := a φ(x; c, s) ∈ ℝⁿ, where φ(x; c, s) = exp(−‖x−c‖²/2s²).

**Network as a measure.**  f̂(x) = (1/m)Σ_i Φ(x; ϑ_i) = ∫Φ(x;ϑ) dμ(ϑ), with μ := (1/m)Σ_i δ_{ϑ_i}.
The map μ ↦ f̂ is linear.

**Metric.**  G = diag(g_a I_n, g_c I_n, g_s) ≻ 0 is the gain matrix. Distances are measured in
‖v‖²_{G⁻¹} = vᵀG⁻¹v, and λ_G := λ_max(G).

**Admissible box and projection.**  Θ_K := {‖a‖ ≤ a_max} × {‖c‖ ≤ c_max} × [s_min, s_max], with
s_min > 0. The boundary-layer width is ε_K. It is renamed from "ε" to avoid a clash with the approximation
error. The inner set is Θ_K^{ε} (each bound shrunk by ε_K), and int Θ_K^{ε} is its interior.

Proj_K(ϑ, y) is the per-particle, factor-wise projection: Pomet–Praly norm-ball operators for a and c and
the box operator for s (`projection.py`). It is Lipschitz, C¹ at the inner edge of the layer and only
Lipschitz at the outer edge. The proofs below need only Lipschitz continuity.

**Property (P).**  For ϑ ∈ Θ_K, ϑ' ∈ Θ_K^{ε} and all y,

    ⟨ϑ − ϑ', G⁻¹ Proj_K(ϑ, G y)⟩ ≤ ⟨ϑ − ϑ', y⟩.

On ∂Θ_K the projected velocity has no outward component, so Θ_K is forward invariant (Nagumo; Pomet–Praly 1992).

**Potential.**  g_t(ϑ) := e(t)ᵀ Φ(x(t); ϑ) and F_t(μ) := −∫ g_t dμ = −e(t)ᵀ f̂(x(t)).

**Adaptive law (the algorithm).**

    ϑ̇_i = Proj_K(ϑ_i, G ∇_ϑ g_t(ϑ_i)),   i = 1..m.                                                  (L)

**Ideal measure.**  For j = 1..m, ϑ*_j(t) ∈ Θ_K^{ε} is C¹ in t, and μ*_t := (1/m)Σδ_{ϑ*_j(t)}. Write
ε(x,t) := f(x,t) − ∫Φ(x;·)dμ*_t, with ‖ε‖ ≤ ε_max on the operating set X, and

    V_*(t) := ((1/m)Σ_j ‖ϑ̇*_j‖²_{G⁻¹})^{1/2},   V_* := sup_t V_*(t).

The proofs never use optimality of μ*. It is an analytic device, and any admissible C¹ curve may be used.

**Optimal-transport distance.**  For uniform atomic μ = (1/m)Σδ_{ϑ_i} and ν = (1/m)Σδ_{ϑ'_j},

    W(μ, ν)² := min_{π ∈ S_m} (1/m) Σ_i ‖ϑ_i − ϑ'_{π(i)}‖²_{G⁻¹}.

This is the Kantorovich 2-Wasserstein distance with ground cost ‖·‖²_{G⁻¹}: for equal masses the
optimum is attained at a permutation, by Birkhoff's theorem. It is computed by linear assignment in
O(m³). For general measures on Θ_K, W denotes the same Kantorovich distance.

**Label-free Lyapunov function.**

    V_OT(t) := ½‖e(t)‖² + ½ W(μ_t, μ*_t)² = min_π V_π(t),
    V_π := ½‖e‖² + (1/2m) Σ_i ‖ϑ_i − ϑ*_{π(i)}‖²_{G⁻¹}.

**Standing hypotheses (S).**

- **(S1)** f is continuous in (x,t) and locally Lipschitz in x, uniformly in t on bounded time intervals.
- **(S2)** x_m ∈ C¹, and x_m and ẋ_m are bounded.
- **(S3)** ϑ_i(0) ∈ Θ_K.
- **(S4)** ϑ*_j are C¹ with values in Θ_K^{ε}, and V_* < ∞.
- **(S5)** x(t) ∈ X on the interval considered. For bounded f this is discharged by Theorem D(b).

Under (S1)–(S3) the right-hand side of (L) is continuous in t and locally Lipschitz in (x, ϑ). Solutions
are classical (C¹) and unique, and each V_π is C¹.

**Constants of the Gaussian** (closed forms; `transport.gaussian_constants`).  With u := ‖x−c‖²/s² and
φ = e^{−u/2}, on s ≥ s_min:

    ‖∇_c φ‖ ≤ J_c := e^{−1/2}/s_min,   |∂_s φ| ≤ J_s := 2e^{−1}/s_min,   J := (J_c² + J_s²)^{1/2},
    ‖∇²_{(c,s)} φ‖ ≤ H_φ := (max(h_cc, h_ss) + h_cs)/s_min²,
    h_cc = 1,  h_ss = 2e^{−1/2},  h_cs = √u e^{−u/2}|u−2| evaluated at u = (5−√17)/2.

J and H_φ are upper bounds. The true sup of ‖∇_{(c,s)}φ‖ is about 0.92/s_min, and that of
‖∇²_{(c,s)}φ‖ about 1.57/s_min². The Hessian of the lifted potential is the block matrix
[[0, u∇φᵀ],[∇φuᵀ, (uᵀa)∇²φ]]. Hence

    sup_{x ∈ ℝⁿ, ϑ ∈ Θ_K, ‖u‖=1} ‖∇²_ϑ (uᵀΦ(x;ϑ))‖ ≤ H_K := J + a_max H_φ,     K := ½H_K.          (1)

A brute-force check over the full lifted Hessian gives sup ≈ 0.76·H_K (review round 1).

---

## 2. Theorem A — identification

**Theorem A.**

(i) Assume scalar classical gains Γ_W = γ_W I, Γ_c = γ_c I and Γ_σ = γ_σ I. Under the dictionary
W_i = a_i/m, γ_W = g_a/m, γ_c = m g_c and γ_σ = m g_s, the *raw* (unprojected) rates of (L) and of the
classical law of `THEORY_REVISED.md` §2 coincide. The *projected* rates coincide if and only if the
classical joint projections are replaced by per-kernel projections. The classical law projects W onto
one Frobenius ball and vec(c) onto one Euclidean ball. The per-kernel version projects W_i onto
‖W_i‖ ≤ a_max/m with layer ε_K/m, c_i onto ‖c_i‖ ≤ c_max with layer ε_K, and σ_i onto [s_min, s_max]
with layer ε_K.

(ii) Fix a closed-loop trajectory (x(·), e(·)). The field v_t(ϑ) := Proj_K(ϑ, G∇g_t(ϑ)) is
Lipschitz in ϑ on Θ_K, and Θ_K is forward invariant. For every μ_0 ∈ P(Θ_K), μ_t := (X_t)#μ_0 is a
distributional solution of ∂_tμ + ∇·(μ v_t) = 0. On int Θ_K^{ε}, v_t = −G∇(δF_t/δμ). So, at frozen time,
μ_t is the Wasserstein gradient flow (metric G⁻¹) of the linear functional F_t, and (L) is its particle
form. Because F_t is generated by the closed loop, the coupled system is the one studied in Theorem D.

(iii) **Dictionary**, all in the G⁻¹ metric:

| classical device | measure-level object |
|---|---|
| projection onto Θ_K | ε_K-regularised support constraint μ_t ∈ P(Θ_K) ([46] Prop. 3.1 for the exact W₂ projection) |
| σ-modification toward ϑ⁰ | potential energy (σ/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ (geodesically σ-convex), velocity −σ(ϑ−ϑ⁰) |
| σ-modification toward a labelled (ϑ⁰_i) | not a functional of μ; it is a functional on the labelled lift ℝ^{mp} |
| OT σ-modification (Theorem C) | (σ/2) W(μ, μ⁰)²; for atomic μ with a unique optimal assignment its gradient is ϑ_i − ϑ⁰_{τ(i)}, and at ties a measurable selection is used |
| e-modification | the same functionals with coefficient σ‖e‖ |
| parameter dither (Theorem H) | entropy τ∫μ log μ |

*Proof.* (i) We have ∇_a g = φe, ∇_c g = (aᵀe)∇_cφ and ∂_s g = (aᵀe)∂_sφ. Matching these with the
classical rates gives the gain map. For the projections, the scaling identity
(1/m)Proj_{r,ε}(a, y) = Proj_{r/m, ε/m}(a/m, y/m) holds because h_{r/m,ε/m}(‖a‖/m) = h_{r,ε}(‖a‖). With the
layer not rescaled the ramps differ. A joint projection couples different kernels, and a per-kernel one
cannot reproduce that. `test_identification_projected_rates_per_kernel` checks this with active projections,
and also checks that the statement is sharp.

(ii) Proj_K is Lipschitz in ϑ. The ramp vanishes for ‖ϑ‖ ≤ θ_max − ε_K, so the direction ϑ/‖ϑ‖ is only
evaluated away from 0. Also ∇g_t is C¹ on s ≥ s_min. Forward invariance follows from Nagumo's theorem,
since on the boundary ϑᵀProj ≤ 0. The push-forward of μ_0 by a Lipschitz flow solves the continuity
equation (Golse [82]; AGS §8.1). On int Θ_K^{ε} the projection is inactive, so v_t = G∇g_t = −G∇(δF_t/δμ).

(iii) The rows follow from the computation in the metric G⁻¹: for any V, ∇_{G⁻¹}V = G∇V. For example,
V = (σ/2)‖ϑ−ϑ⁰‖²_{G⁻¹} gives G∇V = σ(ϑ−ϑ⁰). The fixed-grid reading of the earlier version has been
withdrawn as unproved.  ∎

*Remark A.1 (normalisation).*  The algorithm's raw rates are unchanged. Two things do change: the
admissible set becomes a product of per-kernel balls, and a_max (not ‖W‖_F) is held fixed in m. In this
normalisation, all constants below are free of m. The classical convention of fixed γ's corresponds to
g_c, g_s ∝ 1/m, so centres and widths freeze as m grows.

---

## 3. Theorem B — the optimal-assignment Lyapunov function

**Lemma B.1 (Dini derivative of a minimum).**  Suppose each V_π is C¹ and V_OT = min_π V_π. Then V_OT is
locally Lipschitz, and D⁺V_OT(t) ≤ V̇_π(t) for every π active at t.
*Proof.*  For h > 0, V_OT(t+h) ≤ V_π(t+h) and V_OT(t) = V_π(t).  ∎

**Lemma B.2 (matched remainder, two regimes).**  Let L₀ := 2a_max and L₁ := (1 + a_max²J²)^{1/2}. For
ϑ, ϑ' ∈ Θ_K, the remainder R := g_t(ϑ') − g_t(ϑ) − ⟨∇g_t(ϑ), ϑ' − ϑ⟩ satisfies

    |R| ≤ ‖e(t)‖ · min{ K‖ϑ − ϑ'‖², L₀ + L₁‖ϑ − ϑ'‖ }.

*Proof.*  If e(t) = 0 both sides vanish. Otherwise set h := ϑ' − ϑ. Θ_K is convex, so the Lagrange
remainder of the scalar function τ ↦ g_t(ϑ + τh) gives R = ½hᵀ∇²g_t(ϑ + τ₀h)h for some τ₀ ∈ (0,1). Since
g_t = ‖e‖·(uᵀΦ) with u = e/‖e‖, bound (1) gives the quadratic branch. For the linear branch,
|g_t| ≤ ‖e‖a_max and ‖∇g_t‖ ≤ ‖e‖L₁.  ∎

**Theorem B.**  Under (S), along every solution of (L), for every t, with W := W(μ_t, μ*_t):

    D⁺V_OT ≤ −k_e‖e‖² + ‖e‖ ( ε_max + min{ Kλ_G W², L₀ + L₁ λ_G^{1/2} W } ) + V_*(t) W.        (B)

*Proof.*  Fix t and an optimal π, and write ϑ̃_i := ϑ_i − ϑ*_{π(i)}. Then

    V̇_π = eᵀė + (1/m)Σ⟨ϑ̃_i, G⁻¹ϑ̇_i⟩ − (1/m)Σ⟨ϑ̃_i, G⁻¹ϑ̇*_{π(i)}⟩.

(a) eᵀė = −k_e‖e‖² + eᵀε + (1/m)Σ_i [g_t(ϑ*_{π(i)}) − g_t(ϑ_i)].

(b) By (P) with ϑ' = ϑ*_{π(i)}, (1/m)Σ⟨ϑ̃_i, G⁻¹ϑ̇_i⟩ ≤ (1/m)Σ⟨ϑ̃_i, ∇g_t(ϑ_i)⟩.

(c) −(1/m)Σ⟨ϑ̃_i, G⁻¹ϑ̇*_{π(i)}⟩ ≤ ((1/m)Σ‖ϑ̃_i‖²_{G⁻¹})^{1/2}((1/m)Σ_i‖ϑ̇*_{π(i)}‖²_{G⁻¹})^{1/2} = W·V_*(t). Here π is a bijection and is optimal.

Summing, V̇_π ≤ −k_e‖e‖² + ‖e‖ε_max + (1/m)Σ_i R_i + W V_*(t), where each R_i is a matched remainder.
Average Lemma B.2, using (1/m)Σmin{A_i, B_i} ≤ min{avg A, avg B},
(1/m)Σ‖ϑ̃_i‖² ≤ λ_G W² and (1/m)Σ‖ϑ̃_i‖ ≤ λ_G^{1/2}W. Lemma B.1 then completes the proof.  ∎

**Corollary B.3 (ultimate bounds).**  Let Δ := diam_{G⁻¹}(Θ_K), so W ≤ Δ. Define
D^{OT} := ε_max + min{Kλ_GΔ², L₀ + L₁λ_G^{1/2}Δ}, C^{OT} := V_*Δ and
R^{OT} := D^{OT}/(2k_e) + (D^{OT 2}/(4k_e²) + C^{OT}/k_e)^{1/2}. Then:

- D⁺V_OT < 0 whenever ‖e‖ > R^{OT}.
- limsup ‖e‖ ≤ (R^{OT 2} + Δ²)^{1/2}, using V_max^{OT} = Δ²/2.
- RMS_{[T₁,T₂]}(‖e‖) ≤ D^{OT}/(2k_e) + (D^{OT 2}/(4k_e²) + C^{OT}/k_e + V_OT(T₁)/(k_e(T₂−T₁)))^{1/2}.

*Proof.*  As in `THEORY_REVISED.md` Theorem 1. The Dini comparison lemma applied with the continuous,
V-independent right-hand side gives the integrated form.  ∎

*Remark B.4 (what OT buys and what it does not).*

- Globally the OT constant is comparable to the classical one, and not smaller. For the test box the
  global constant L₀ + L₁λ_G^{1/2}Δ is 533, against 101–138 for the classical Lemma 2 constant under the
  dictionary. A block-structured linear branch, |R_i| ≤ ‖e‖(2a_max + ‖ã_i‖ + a_maxJ‖(c̃_i,s̃_i)‖), gives
  L₀ + (g_a + a_max²J²max(g_c,g_s))^{1/2}Δ ≈ 382.
- Locally (Corollary B.5) the excess disturbance is quadratic in the *matched* error.
- A labelled theory can make the same local statement only for a labelling fixed in advance, and only
  while the optimal assignment does not switch. (B) holds uniformly, with the minimum over π taken at each t.
- m-uniformity is a property of the normalisation (Remark A.1), not of OT.

**Corollary B.5 (conditional, a-posteriori sharp bound).**  Suppose sup_{[T₁,T₂]} W(μ_t, μ*_t) ≤ η̄, a
quantity measured along the run, or controlled a priori by Theorem C(iv) or Theorem E. With
D̄ := ε_max + min{Kλ_Gη̄², L₀ + L₁λ_G^{1/2}η̄},

    RMS_{[T₁,T₂]}(‖e‖) ≤ D̄/(2k_e) + ( D̄²/(4k_e²) + V_*η̄/k_e + V_OT(T₁)/(k_e(T₂−T₁)) )^{1/2},

which tends to ε_max/k_e as η̄ → 0 and T₂ → ∞.

**Corollary B.6 (approximation error as quantisation error).**  Suppose
f = ∫Φ dμ*_∞ + ε_∞ for a probability measure μ*_∞ on Θ_K^{ε}, with ‖ε_∞‖ ≤ ε̄_∞. Then for every m-atom
μ*_m supported in Θ_K^{ε}, ε_max(m) ≤ ε̄_∞ + L₁λ_G^{1/2} W(μ*_m, μ*_∞). Here ϑ ↦ Φ(x;ϑ) is L₁-Lipschitz,
and Kantorovich–Rubinstein gives W₁ ≤ W₂ for the same ground metric. With μ*_m an optimal m-point
quantiser, W(μ*_m, μ*_∞) = O(m^{−1/p}) for absolutely continuous μ*_∞ (Zador; Graf–Luschgy [S9]). In this
formulation the curse of dimensionality sits in the number of atoms needed to quantise the ideal measure
on the lifted space, and in no Lyapunov constant.

*Remark B.7 (continuum form, sketch).*  For μ_t = (X_t)#μ_0 and an absolutely continuous ideal curve μ*_t
in (P(Θ_K^{ε}), W) with velocity field w_t, set V_* := sup_t ‖w_t‖_{L²(μ*_t;G⁻¹)}. Then (B) holds for
a.e. t. This uses the two-curve derivative formula
d/dt ½W² = ∫⟨ϑ − ϑ', v_t(ϑ) − w_t(ϑ')⟩_{G⁻¹}dγ_t for optimal plans γ_t, with (P) and Lemma B.2 integrated
against γ_t. It is a sketch; the exact reference for the two-curve formula still needs to be fixed.

---

## 4. Theorem C — optimal-transport σ-modification

Let μ⁰ = (1/m)Σδ_{ϑ⁰_k} be a prior configuration with ϑ⁰_k ∈ Θ_K^{ε}. The **OT σ-modified law** is

    ϑ̇_i = Proj_K( ϑ_i, G∇g_t(ϑ_i) − σ_m(ϑ_i − ϑ⁰_{τ(i)}) ),   τ(t) ∈ argmin_ρ (1/m)Σ_i‖ϑ_i − ϑ⁰_{ρ(i)}‖²_{G⁻¹}.        (L_σ)

**Solution concept.**  The right-hand side is discontinuous at assignment ties. Solutions are taken in
the Carathéodory/Filippov sense: any measurable selection of optimal assignments gives a bounded field, so
solutions exist. V_OT and Z below are locally Lipschitz, and every inequality is understood for a.e. t
(or for D⁺ at times where the optimal assignment is locally unique). Because (P) is linear in y, it holds
for Filippov convex combinations. Uniqueness may fail at ties, and every statement holds for every solution.

**Notation.**  Write c_ρ(ϑ) := (1/m)Σ‖ϑ_i − ϑ⁰_{ρ(i)}‖²_{G⁻¹}. Let η := W(μ_t, μ*_t), R := Δ = diam_{G⁻¹}(Θ_K),
b := sup_t W(μ*_t, μ⁰) and Δ_gap := inf_t min{ c_ρ(ϑ*(t)) − b(t)² : c_ρ(ϑ*(t)) > b(t)² }, with
Δ_gap := +∞ if the set is empty.

**Theorem C.**

(i) *Dissipation toward the prior.*  Let Z := ½W(μ_t, μ⁰)². Then D⁺Z ≤ −2σ_mZ + (2Z)^{1/2}Γ(t), with
Γ := ((1/m)Σ‖G∇g_t(ϑ_i)‖²_{G⁻¹})^{1/2} ≤ λ_G^{1/2}L₁‖e‖. Hence
W(μ_t,μ⁰) ≤ max{W(μ_0,μ⁰), sup Γ/σ_m}. For the pure-leak flow (Γ ≡ 0), the decay
W(μ_t,μ⁰) = e^{−σ_m t}W(μ_0,μ⁰) is exact for every G ≻ 0 (Lemma C.3).

(ii) *Induced coupling.*  For optimal π* (μ to μ*) and τ (μ to μ⁰),
C := (1/m)Σ‖ϑ*_{π*(i)} − ϑ⁰_{τ(i)}‖²_{G⁻¹} ≤ b(t)² + 4ηR, and C = b(t)² whenever 4ηR < Δ_gap.

(iii) *Lyapunov inequality.*

    D⁺V_OT ≤ −k_e‖e‖² + ‖e‖(ε_max + Kλ_Gη²) + V_*η − σ_mη² + σ_mηC^{1/2}.

Under 4ηR < Δ_gap this becomes

    D⁺V_OT ≤ −k_e‖e‖² + ‖e‖(ε_max + Kλ_Gη²) + V_*η − (σ_m/2)η² + (σ_m/2)b².

(iv) *Local exponential estimate.*  Let η̄ < Δ_gap/(4R) and 8K²λ_G²η̄² ≤ k_eσ_m. Set
α := min(k_e, σ_m/2) and β := ε_max²/k_e + (σ_m/2)b² + 2V_*²/σ_m. If max{V_OT(0), β/α} < η̄²/2, then
{V_OT ≤ c} is forward invariant for every c ∈ [max{V_OT(0), β/α}, η̄²/2), and
V_OT(t) ≤ e^{−αt}V_OT(0) + β/α.

**Scope.**  Δ_gap ≤ 2R²/m generically: swapping two matched particles changes the averaged cost by at most
2R²/m. Hence η̄ < R/(2m), and non-vacuity forces b < Δ_gap/(4√2R) = O(1/m). Part (iv) is therefore a
small-m, near-prior local result, in contrast with Theorem B.

(v) *The labelled leak.*  The labelled σ-modification (leak toward ϑ⁰_i) satisfies the global strict
inequality of (iii) with bias b_lab(π) := ((1/m)Σ‖ϑ*_{π(i)} − ϑ⁰_i‖²_{G⁻¹})^{1/2}, but for the labelled
function V_π (η replaced by the π-labelled error), for any fixed π. Choosing the minimising π gives bias b.
It does **not** satisfy (iii) for V_OT (Proposition C.5): the labelled leak can increase the
optimal-transport Lyapunov function.

*Proof.*

(i) Apply (P) with ϑ' = ϑ⁰_{τ(i)} and Cauchy–Schwarz in the G⁻¹ inner product:
Ż_τ ≤ (2Z)^{1/2}Γ − 2σ_mZ. Write y := (2Z)^{1/2}. Where Z > 0, D⁺y ≤ −σ_m y + Γ. At Z = 0,
D⁺y ≤ ((1/m)Σ‖ϑ̇_i‖²_{G⁻¹})^{1/2} ≤ Γ, because the projection is a G⁻¹-contraction:
‖Proj(ϑ,Gy)‖²_{G⁻¹} = ‖y‖²_G − h(2−h)(ϑᵀGy)²/(ϑᵀGϑ). The bound on Γ uses Lemma B.2. The pure-leak case is
Lemma C.3, and the projection is inactive there: in the layer aᵀ(a⁰ − a) ≤ ‖a‖(‖a⁰‖ − ‖a‖) ≤ 0, and
similarly for c, while for s the leak points away from the nearby face.

(ii) Relabel so that (1/m)Σ‖ϑ_i − ϑ*_i‖²_{G⁻¹} = η². For any ρ, |c_ρ(ϑ) − c_ρ(ϑ*)| ≤ 2ηR. With ρ₀ optimal
for (ϑ*, ϑ⁰), C = c_τ(ϑ*) ≤ c_τ(ϑ) + 2ηR ≤ c_{ρ₀}(ϑ) + 2ηR ≤ b² + 4ηR. No assignment cost lies strictly
between b² and b² + Δ_gap.

(iii) Theorem B with y = ∇g_t − σ_mG⁻¹(ϑ_i − ϑ⁰_{τ(i)}) in (P). The extra term is
−σ_m(1/m)Σ⟨ϑ̃_i, ϑ̃_i + ϑ*_{π*(i)} − ϑ⁰_{τ(i)}⟩_{G⁻¹} ≤ −σ_mη² + σ_mηC^{1/2}. Then use σ_mηb ≤ (σ_m/2)(η² + b²).

(iv) Apply Young's inequality three times:
‖e‖ε_max ≤ (k_e/4)‖e‖² + ε_max²/k_e;
‖e‖Kλ_Gη² ≤ (k_e/4)‖e‖² + (σ_m/8)η² on {η ≤ η̄};
V_*η ≤ (σ_m/8)η² + 2V_*²/σ_m.
This gives D⁺V_OT ≤ −αV_OT + β. Each set {V_OT ≤ c} with β/α < c < η̄²/2 lies inside {η ≤ η̄} and has a
negative boundary derivative. Intersect over c > β/α for the endpoint, and apply the comparison lemma.

(v) Use the classical completion of squares for V_π, since steps (a)–(c) of Theorem B hold for any fixed π.
The negative statement is Proposition C.5.  ∎

**Lemma C.3 (no switching along the pure leak; exact decay, any G ≻ 0).**  Let τ be optimal for (μ_0, μ⁰)
and ϑ̇_i = −σ_m(ϑ_i − ϑ⁰_{τ(i)}). Then τ stays optimal and W(μ_t,μ⁰) = e^{−σ_m t}W(μ_0,μ⁰).
*Proof.*  ϑ_i(t) − ϑ⁰_{τ(i)} = λd_i with λ = e^{−σ_m t}. For any ρ,
c_ρ(t) − c_τ(t) = (2λ/m)Σ⟨d_i, Δ_i⟩_{G⁻¹} + (1/m)Σ‖Δ_i‖²_{G⁻¹}, where Δ_i := ϑ⁰_{τ(i)} − ϑ⁰_{ρ(i)}. This
is affine in λ, non-negative at λ = 1 by optimality and at λ = 0 trivially, hence non-negative on [0,1].
At an exact tie at t = 0 the chosen assignment becomes strictly optimal immediately, so there is no
chattering.  ∎
Checked: zero switches, and exact decay to integrator accuracy.

**Proposition C.5 (the labelled leak reshuffles).**  There are configurations with b = b_lab = 0 and e = 0
at which the labelled leak gives D⁺V_OT > 0.
*Example.*  Take n = 1, G = I, ideal = prior = {a = 0, a = 1} (other coordinates equal), and network
weights (1−δ, δ). The optimal assignment to the ideal is the swap. The labelled leak moves each particle
toward the *other* ideal atom, so D⁺V_OT = σ_m(δ − δ²) > 0. In general, from a relabelled start the
labelled leak transports every kernel across the domain to "its" prior position. The steady states
coincide, but the transient costs tracking error, as Theorem F prices it
(`results/ot_sigma_mod/transient.md`: particle path 5.8× longer).

*Remark C.6 (the obstruction is known).*  Several standard facts limit what a label-free leak can do:

- ½W²(·, ν) is not λ-convex along geodesics for any λ in dimension ≥ 2. W² is semi-concave there
  (AGS Thm 7.3.2; Example 9.1.5; Santambrogio [71] §4.4, two-atom switching example).
- It is 1-convex along generalised geodesics with base ν (AGS Lemma 9.2.1; Carlen–Craig [73]).
- Label-free geodesically convex leaks exist, for example toward the barycentre of μ⁰, but they do not
  vanish on the S_m-orbit of μ⁰.

So the trade-off between a global guarantee with a labelled bias and a local guarantee with the optimal
bias reflects this geometry.

---

## 5. Theorem D — the mean-field closed loop

The state is (x, μ) ∈ ℝⁿ × P(Θ_K), with

    ẋ = f(x,t) + ẋ_m − k_e(x − x_m) − ∫Φ(x;ϑ)dμ_t,    ∂_tμ_t + ∇·(μ_t v_t) = 0,    v_t = Proj_K(ϑ, G∇_ϑ[(x−x_m)ᵀΦ(x;ϑ)]),

and solutions are Lagrangian, μ_t = (X_t)#μ_0. This is law (L). The OT leak is multivalued at ties and is
excluded here.

**Theorem D.**  Assume (S1), (S2) and G ≻ 0.

(a) *Local well-posedness.*  For every (x_0, μ_0) there is a unique Lagrangian solution, which continues
while x stays bounded.

(b) *Adaptation-free tube (m-uniform, μ*-free).*  For every μ ∈ P(Θ_K), ‖∫Φ(x;·)dμ‖ ≤ a_max. Suppose
F_max := sup ‖f(x,t)‖ over the tube 𝒯 := {x : ‖x − x_m(t)‖ ≤ ‖e(0)‖ + (F_max + a_max)/k_e} is finite.
Then every solution, for every m and every μ_0, stays in 𝒯 for all t ≥ 0, and

    ‖e(t)‖ ≤ e^{−k_e t}‖e(0)‖ + (F_max + a_max)/k_e.

(c) *Dobrushin estimate.*  Let two solutions remain in a compact X on [0,T]. Then

    ‖x_t − x'_t‖ + W(μ_t,μ'_t) ≤ e^{Λt}(‖x_0 − x'_0‖ + W(μ_0,μ'_0)),

with Λ := max(L_ϑ + L_2, L_x + L_1) independent of m:

- L_ϑ is the Lipschitz constant of v in ϑ in the G⁻¹ norm. It is O(1/ε_K) because of the projection ramp.
- L_x is the Lipschitz constant of v in (x, e).
- L_1 := Lip_X f + k_e + sup‖∇_x f̂‖.
- L_2 := λ_G^{1/2}L₁ is the W-Lipschitz constant of μ ↦ ∫Φ(x;·)dμ.

(d) *Mean-field limit.*  If x^m_0 → x_0 and μ^m_0 → μ_0 in W, the m-particle loop converges uniformly
on [0,T] to the measure-valued loop. Part (b) supplies a common X.

(e) *Passage of (B) to the limit.*  Assume m-atom ideals μ*^m_t → μ*_t uniformly in W on [0,T], where μ*_t
is a W-Lipschitz curve in P(Θ_K^{ε}), together with sup_m ε^m_max ≤ ε_max and sup_m V^m_* ≤ V_*. Then the
limit satisfies (B) in integrated form, and for every t where W(μ_t,μ*_t) is continuous it satisfies (B)
in Dini form with constants (k_e, ε_max, K, L₀, L₁, λ_G, V_*).

*Proof.*

(a) The vector field is Lipschitz on bounded sets, by the proof of A(ii) and (S1). This gives a fixed
point in C([0,T]; ℝⁿ × P₂) (Dobrushin; Golse [82]).

(b) ė = −k_e e + (f − f̂) with ‖f − f̂‖ ≤ F_max + a_max as long as x ∈ 𝒯. Variation of constants and a
first-exit argument finish the proof.

(c) Couple with γ₀ optimal and push forward by (X_t, X'_t), so W(μ_t,μ'_t) ≤ D(t). Then
d/dt D ≤ L_ϑD + L_x‖x − x'‖, using Cauchy–Schwarz on the coupling, and
d/dt‖x − x'‖ ≤ L_1‖x − x'‖ + L_2D. Apply Grönwall to the sum. The stated constants assume Lipschitz
continuity in x uniform in t. Without it the estimate fails (review round 1: an explicit counterexample).

(d) By (b) and (c).

(e) Each V^m is locally Lipschitz, so V^m(t₂) − V^m(t₁) ≤ ∫R^m. The right-hand sides converge uniformly
by the triangle inequality for W, so V(t₂) − V(t₁) ≤ ∫R, and R is continuous wherever t ↦ W(μ_t,μ*_t) is.  ∎

Numerics (`results/mean_field_m`): halving W(μ^m_0, μ^{2m}_0) halves W(μ^m_T, μ^{2m}_T), with a constant
factor of 1.68. sup_t‖e^m − e^{2m}‖ falls from 4.6e-3 (16 vs 32) to below 1e-4 (128 vs 256).

---

## 5b. Theorem F — the transport cost of adaptation

**Theorem F.**  Along (L), for 0 ≤ s ≤ t,

- (i) ((1/m)Σ‖ϑ̇_i‖²_{G⁻¹})^{1/2} ≤ λ_G^{1/2}L₁‖e(t)‖;
- (ii) W(μ_t,μ_s) ≤ λ_G^{1/2}L₁∫_s^t‖e‖;
- (iii) W(μ_{T₂},μ_{T₁}) ≤ λ_G^{1/2}L₁(T₂ − T₁)·RMS_{[T₁,T₂]}‖e‖.

Along (L_σ), add σ_mW(μ_r,μ⁰) inside the integral.
*Proof.*  (i) uses the projection's G⁻¹-contraction property (Theorem C(i)) and Lemma B.2. (ii) follows
from the identity coupling and Minkowski's inequality, or AGS Thm 8.3.1 for general measures. (iii) is
Cauchy–Schwarz.  ∎
*Reading.*  Moving the kernel measure by W costs at least W/(λ_G^{1/2}L₁) of integrated tracking error.
This is an exchange rate between learning and tracking. Checked at every sample.

---

## 6. Theorem E — label-free local convergence under tangential excitation

**Assumptions.**

- ε ≡ 0 on a neighbourhood of the reference orbit, and μ* is constant with ϑ*_j ∈ **int** Θ_K^{ε}.
- x_m ∈ C¹, with x_m and ẋ_m bounded.
- Define Ψ(t) := [∇_ϑΦ(x_m(t);ϑ*_1), …, ∇_ϑΦ(x_m(t);ϑ*_m)] ∈ ℝ^{n×mp}.
- **(PE)** There are T, α > 0 with ∫_t^{t+T}ΨᵀΨ ⪰ αI_{mp} for all t ≥ 0.

(PE) forces the ϑ*_j to be pairwise distinct, with a*_j ≠ 0. Otherwise two column blocks coincide or
vanish, and the Gramian is singular. Let d² := min_{π≠π'}(1/m)Σ‖ϑ*_{π(i)} − ϑ*_{π'(i)}‖²_{G⁻¹} > 0.

**Theorem E.**  There are r₀, k, λ > 0, determined by the linearisation below, such that the following
holds. Put c₀ := min{ r₀²/(2(1 + mλ_G)), d²/8 } and C := √2·k·max{1,(mλ_G)^{1/2}}·max{1,(mλ_min(G))^{−1/2}}.
If V_OT(t₀) ≤ c₀, then:

- the optimal assignment π at t₀ is unique, and it is the optimal assignment for all large t;
- ‖e(t)‖ + W(μ_t,μ*) ≤ Ce^{−λ(t−t₀)}(‖e(t₀)‖ + W(μ_{t₀},μ*));
- ϑ_i(t) → ϑ*_{π(i)}, and the labelled error tends to ((1/m)Σ‖ϑ*_{π(i)} − ϑ*_i‖²_{G⁻¹})^{1/2}, which is
  non-zero for π ≠ id.

The constants depend on m, G, k_e, T, α, sup‖Ψ‖, sup‖Ψ̇‖, a_max and s_min.

*Proof.*

1. *Symmetry.*  Relabelling the network, (e,ϑ) ↦ (e, ϑ∘σ), maps solutions to solutions. It fixes e and is
   an isometry of the matched metric, and it conjugates the dynamics near (0,ϑ*∘π) to those near (0,ϑ*).
   It therefore suffices to treat π = id.

2. *Linearisation.*  On a ball D ⊂ int Θ_K^{ε} around ϑ* the projection is the identity. There the closed
   loop is C^∞ in (e,ϑ̃) and continuous in t. Its second derivatives are bounded uniformly in t, because the
   Gaussian's derivatives through order three are bounded uniformly in x for s ≥ s_min − r. Its Jacobian at
   0 is A(t) = [[−k_eI, −(1/m)Ψ], [(I_m⊗G)Ψᵀ, 0]], and the remainder is O(‖z‖²). The Jacobian and the
   identity QA + AᵀQ = diag(−k_eI, 0) with Q = diag(I/2, (I⊗G⁻¹)/(2m)) are checked numerically by
   `test_theorem_E_linearisation_and_lyapunov_identity`.

3. *LTV stability.*  Put ζ := m^{−1/2}(I_m⊗G^{−1/2})ϑ̃ and Ψ_G := m^{−1/2}Ψ(I⊗G^{1/2}). Then ė = −k_ee − Ψ_Gζ
   and ζ̇ = Ψ_Gᵀe. This is ż = [A₀ + B(t)]z with A₀ = diag(−k_eI, 0) and B skew-symmetric, bounded and with
   bounded derivative. Moreover ∫Ψ_GᵀΨ_G ⪰ (αλ_min(G)/m)I. Uniform exponential stability then follows from
   Morgan–Narendra [S8] and Anderson (1977).

4. *Nonlinear system.*  Khalil Thm 4.13 applies, since the Jacobian is bounded and Lipschitz on D uniformly
   in t. It gives local uniform exponential stability with (r₀, k, λ).

5. *Transfer to (‖e‖, W).*  V_π(t₀) ≤ c₀ gives ‖z(t₀)‖² ≤ 2c₀(1 + mλ_G) ≤ r₀². The bound
   W ≤ (mλ_min(G))^{−1/2}‖ϑ̃‖ and the optimality of π at t₀ give C. If two assignments satisfied
   V_π, V_π' ≤ c₀, the matched triangle inequality would give
   (1/m)Σ‖ϑ*_{π(i)} − ϑ*_{π'(i)}‖²_{G⁻¹} ≤ 8c₀ < d².  ∎

*Remark E.1.*  "The parameters converge" is false for the labelled error and true for the matched one.
The rate is set by αλ_min(G)/m and is small in practice. The 60 s Gramian has λ_min between 1e-3 and 3e-2
(`results/pe_convergence`). The convergence test is therefore qualitative; the proof's linearisation is
the part that is machine-checked.

**Corollary E.2 (uniform-in-time continuity under PE).**  If V_OT(0) ≤ c₀ and V'_OT(0) ≤ c₀, then
‖x_t − x'_t‖ + W(μ_t,μ'_t) ≤ 2Ce^{−λt}(‖e(0)‖ + ‖e'(0)‖ + W(μ_0,μ*) + W(μ'_0,μ*)), by the triangle
inequality through (x_m, μ*).

---

## 6b. Proposition G — the slow geometry performs variable-projection descent

**Setting.**

- f = f(x) is time-invariant and C².
- x_m is T₀-periodic with occupation measure ρ := (1/T₀)∫₀^{T₀}δ_{x_m(t)}dt.
- Law (L) acts in the interior, with gains (g_a, δγ_c, δγ_s).
- L̄(a,c,s) := ½∫‖f − (1/m)Σa_iφ(·;c_i,s_i)‖²dρ.
- M(c,s) := (1/m²)∫φφᵀdρ ⊗ I_n.
- a*(c,s) := argmin_a L̄, and L*(c,s) := L̄(a*(c,s),c,s) is the profile (variable-projection) loss.

**Hypotheses.**

- (H1) M ⪰ m₀I on a compact set 𝒢 of geometries.
- (H2a) Two scales: k_e ≫ ω_ref := 2π/T₀ ≫ every parameter rate.
- (H2b) Three scales: additionally, the *slowest* weight mode g_a m₀ m/k_e is much faster than the
  geometry rate (mδ/k_e)‖∇²L*‖.

(H2b) is a condition on λ_min of the feature Gram matrix, not λ_max. Overlapping kernels make it fail
(see the negative case below).

**Proposition G (formal reduction, two parts).**

(a) *Two scales, under (H2a).*  The period-averaged closed loop is the joint gradient flow
ϑ̇_i = −(m/k_e)G∇_{ϑ_i}L̄(ϑ_1..ϑ_m), whose measure form is −(1/k_e)G∇(δL̄/δμ), the noiseless and
leak-free case of Theorem H(i). In particular L̄ evaluated at the current weights is non-increasing.

(b) *Three scales, under (H2b).*  On time intervals of length O(1/δ), (c,s)(t) = (c̄,s̄)(δt) + o(1), where

    c̄_i' = −(mγ_c/k_e)∇_{c_i}L*(c̄,s̄),    s̄_i' = −(mγ_s/k_e)∂_{s_i}L*(c̄,s̄),

so L* is non-increasing along the reduced flow and the weights track a*(c̄,s̄). In measure form, the
geometry marginal follows the Wasserstein gradient flow of the variable-projection functional of [S2].

*Numerics* (`results/variable_projection`).

- **Two scales hold.**  The configuration with m = 4 overlapping kernels satisfies (H2a) but violates
  (H2b): its slowest weight mode is 0.1× the geometry rate. There, L̄ at the current weights decreases
  at **every** sample (0.219 → 0.008).
- **Three scales need (H2b).**  In the same run L* *increases* (3.8e-3 → 5.9e-3), and the geometry
  velocity is not aligned with −∇L* (median cosine −0.17, then 0.40).
- **Three scales hold under (H2b).**  The configuration has m = 2 well-separated kernels, and its
  measured rate ratios are k_e/ω_ref = 20, ω_ref/(fastest weight) = 7.3 and
  (slowest weight)/(geometry) = 77.5 (`table_satisfying.md`, 3000 s).
  - The period-averaged geometry velocity is aligned with −(m/k_e)Γ∇L* at every sample (median cosine 0.999).
  - Its magnitude is 0.94–0.95 of the prediction, consistent with an averaging error of order
    (fastest weight)/ω_ref ≈ 0.14.
  - The weights track a*(c,s) to 0.1%, and L* decreases monotonically from 1.20e-2 to 4.11e-3.

*Derivation.*  The reduction chains three standard theorems; error terms are not made explicit.

1. Tikhonov for the error channel (Khalil Thm 11.1): e = r/k_e + O(ω_ref/k_e²).
2. Period averaging for the weights (Khalil Thm 10.4; Sastry–Bodson Ch. 4):
   ȧ = −(g_am/k_e)∇_aL̄ + o(·). This is linear in a and exponentially stable under (H1).
3. Two-time-scale Tikhonov for weights versus geometry (Kokotović–Khalil–O'Reilly) gives a = a*(c,s) + O(δ).
   Substituting, the slow equation is ċ_i = −(δγ_cm/k_e)∇_{c_i}L̄(a*,c,s) after averaging.
4. The envelope theorem gives ∇_{c_i}L̄(a*(c,s),c,s) = ∇_{c_i}L*.

*Status.*  This is a bridge. In machine learning the same reduction is a *design* choice of training
algorithms on a given data distribution [S1, S2]. Here it is produced by the closed loop: the time scales
come from k_e and the reference period, and the data distribution is the occupation measure of the
reference. The adaptive controller is thus a variable-projection learner on its own excitation. Numerical
check: `results/variable_projection`.

---

## 6c. Theorem H — Langevin adaptive control and global optimality of kernel placement

**Why something new is needed.**  Theorems B–G are local, or describe the dynamics. A global statement
about kernel placement needs either the Chizat–Bach theory or entropy.

The Chizat–Bach theory for noiseless Wasserstein flows (partially 1-homogeneous case, [S6]) requires
boundary conditions: Neumann conditions on a bounded convex domain, or C¹-uniform behaviour at infinity.
Gaussian kernels with free centres and widths satisfy neither. In the averaged mean-field limit, the
symmetric "collapsed" configuration is stationary and never separates without noise. This is checked by
`test_noiseless_symmetric_initialisation_never_separates_but_dither_does`. Entropy removes the obstruction.

**Design.**  The particles are unconstrained, ϑ = (α, γ, κ) ∈ ℝ^p, and the features are smoothly
saturated:

    a = a_max tanh(α/a_max),    c = c_max tanh(γ/c_max),    s = s_min + (s_max − s_min)·sig(κ),

so Φ̃(x;ϑ) = aφ(x;c,s), with ‖Φ̃‖ ≤ K_Φ := a_max√n and all ϑ-derivatives bounded. The law is

    dϑ_i = [ G∇_ϑg_t(ϑ_i) − σ(ϑ_i − ϑ⁰) ] dt + √(2τ) G^{1/2} dB_i,                                    (L_τ)

that is, the classical gradient law, σ-modification written in the G⁻¹ metric, and parameter dither of
temperature τ.

**Averaged mean-field system.**  Let ρ be the occupation measure of a T₀-periodic reference, and
F_ρ := sup_{supp ρ}‖f‖. Define r_μ := f − ∫Φ̃(·;ϑ')dμ(ϑ'),

    L̄(μ) := ½∫‖r_μ‖²dρ,    V_L[μ](ϑ) := −∫r_μ(x)ᵀΦ̃(x;ϑ)dρ(x),
    V[μ] := (1/k_e)V_L[μ] + (σ/2)‖ϑ − ϑ⁰‖²_{G⁻¹},

and the averaged system

    ∂_tμ = ∇·(μG∇V[μ]) + τ∇·(G∇μ).                                                                  (𝔄)

**Theorem H.**

(i) *Identification.*  V[μ] is the first variation of 𝒢(μ) := (1/k_e)L̄(μ) + (σ/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ. Hence
(𝔄) is the mean-field Langevin dynamics, in the G⁻¹ metric, of the free energy F_τ := 𝒢 + τH, where
H(μ) = ∫μ log μ. Moreover (𝔄) is the leading-order (Tikhonov + period-averaging) description of the
mean-field limit of the closed loop under (L_τ): replacing e by r/k_e and averaging over one reference
period turns the drift of (L_τ) into −G∇V[μ]. This is checked numerically to 1e-7 by
`test_averaged_drift_is_the_mean_field_langevin_drift`.

(ii) *Global exponential convergence.*  F_τ has a unique minimiser, μ*_τ ∝ exp(−V[μ*_τ]/τ). For every
μ_0 with F_τ(μ_0) < ∞,

    F_τ(μ_t) − F_τ(μ*_τ) ≤ e^{−2τρ_τ t}(F_τ(μ_0) − F_τ(μ*_τ)),
    W_{G⁻¹}(μ_t, μ*_τ)² ≤ (2/(τρ_τ))e^{−2τρ_τ t}(F_τ(μ_0) − F_τ(μ*_τ)),

with log-Sobolev constant

    ρ_τ ≥ (σ/τ) · exp( −2K_Φ(F_ρ + K_Φ)/(k_e τ) ).

(iii) *Tracking.*  Along (𝔄) the quasi-steady tracking error e ≈ r_{μ_t}/k_e satisfies
⟨‖e‖²⟩_ρ ≈ 2L̄(μ_t)/k_e². Here L̄(μ_t) converges to L̄(μ*_τ), and μ*_τ minimises the σ- and
entropy-regularised occupation-weighted approximation problem
(1/k_e)L̄(μ) + (σ/2)∫‖ϑ−ϑ⁰‖²_{G⁻¹}dμ + τH(μ). The entropy is taken in ϑ-coordinates; whitening changes it
only by the constant ½log det G.

(iv) *Annealing.*  If τ_t = α/log t with α > α* := 2K_Φ(F_ρ + K_Φ)/k_e, then
𝒢(μ_t) − inf 𝒢 ≤ C'·log log t/log t. This is [S3] Theorem 4.1, whose hypothesis
ρ_τ ≥ C₀e^{−α*/τ} is exactly the bound in (ii).

(v) *Finite m.*  Under the hypotheses of [S5] (functional convexity, uniform log-Sobolev inequality,
regularity), the m-particle averaged system satisfies uniform-in-time propagation of chaos. The
hypotheses must still be checked verbatim for (𝔄) (`LITERATURE_TRANSPORT.md` open question 4).

*Proof.*

(i) Compute d/dε L̄((1−ε)μ + εν)|₀ = −∫r_μᵀ∫Φ̃d(ν−μ)dρ = ∫V_Ld(ν−μ). The confinement is linear in μ.
In the G⁻¹ metric, ∇·(μG∇(V + τlog μ)) = ∇·(μG∇V) + τ∇·(G∇μ). For the averaging statement, the drift of
(L_τ) with e = r/k_e is (1/k_e)G∇_ϑ(rᵀΦ̃) − σ(ϑ − ϑ⁰). Its period average is −(1/k_e)G∇V_L − σ(ϑ−ϑ⁰),
because r_μ does not depend on ϑ for fixed μ.

(ii) Whiten with ζ := G^{−1/2}ϑ, so the metric is Euclidean and the confinement is (σ/2)‖ζ − ζ⁰‖². Then
verify Assumptions 1–3 of [S3]:

- *Assumption 1.* ∇V[μ](ζ) is Lipschitz in (ζ, μ). In ζ this holds because ∇²Φ̃ is bounded (saturation,
  s ≥ s_min). In μ it holds because ∫Φ̃dμ is W₁-Lipschitz, W₁ ≤ W₂, and ∇Φ̃ is bounded.
- *Assumption 2.* 𝒢 is convex: a linear map followed by a convex quadratic, plus a linear term. F_τ has a
  minimiser by confinement and lower semicontinuity, as in [S3] Prop. 5.1.
- *Assumption 3.* e^{−V[μ]/τ} = e^{−(σ/2τ)‖ζ−ζ⁰‖²}·e^{−V_L[μ]/(k_eτ)}. The Gaussian factor satisfies
  LSI(σ/τ) (Bakry–Émery). Also |V_L[μ]| ≤ ∫‖r_μ‖‖Φ̃‖dρ ≤ (F_ρ + K_Φ)K_Φ, so the perturbation has
  oscillation at most 2K_Φ(F_ρ + K_Φ)/(k_eτ). Holley–Stroock then gives the stated ρ_τ.

[S3] Theorem 3.2 and Corollary 3.3 give the two inequalities. Uniqueness follows from strict convexity of H.

(iii) Combine (i) and (ii).

(iv) and (v) are the cited theorems.  ∎

*Remark H.1 (dither costs nothing in the mean-field limit).*  f̂ = ∫Φ̃dμ is deterministic in the limit,
by the law of large numbers. The only effect of the dither is the entropic bias τH. At finite m the
fluctuation of f̂ is O(K_Φ/√m).

*Remark H.2 (the control knobs in the theorem).*

- σ (σ-modification) is the strong-convexity constant of the confinement.
- τ (dither) is the temperature.
- k_e enters as 1/k_e on the loss. Higher feedback gain flattens the effective landscape, which improves
  ρ_τ, but it also slows learning.
- a_max sets K_Φ.

The guaranteed rate is exponentially small in 1/τ, a known feature of Holley–Stroock bounds. At
practical temperatures it is vacuous: `results/langevin_global` reports a guaranteed rate that underflows
to 0 at τ = 5·10⁻⁴, while escape is observed within hundreds of seconds.

The temperature sweep (`results/langevin_tradeoff`, m = 32, σ = 0.002, 2 seeds) shows the trade-off
directly:

| τ | escape of the collapsed start (s) | loss floor of the spread start |
|---|---|---|
| 0 | never | 7.0e-3 |
| 0.001 | 1005–1480 | 1.5e-2 |
| 0.004 | 430–970 | 3.6e-2 |
| 0.008 | 295–390 | 6.5e-2 |

A hot-start annealed schedule reaches 3.2e-2 from both starts.

The escape time fits t ≈ A + B·log(1/τ) with a 4.5% relative residual, against 25% for an Arrhenius law
exp(b/τ). With four temperatures and two seeds this is indicative only. It suggests that the symmetric
configuration is a *saddle* of the averaged landscape, not a barrier-protected minimum. That would
explain why practice is far faster than the worst-case bound.

*Remark H.3 (why saturation and not projection).*  The mean-field Langevin results are on ℝ^d.
Projection would give reflected dynamics on a convex body, where log-Sobolev inequalities also hold
(Bakry–Émery with a Neumann boundary condition). That route connects to [25] and is left open.

*Remark H.4 (the σ-modification is now essential).*  In the classical theory σ-modification is a
robustness device with a bias. Here it is the confinement that makes the log-Sobolev inequality hold, and
so it is *the reason* the global result exists. The bias is the price, made explicit as a regularisation
term in the variational problem that μ*_τ solves.

**Conjecture H.6 (the real closed loop).**  Consider the non-averaged closed loop under (L_τ), with finite
k_e and gains, in the mean-field limit. For every μ_0 it enters, after a finite time, and stays in a W-ball
around μ*_τ whose radius vanishes with the averaging parameters. Infinite-horizon averaging should follow
from the exponential stability of (𝔄).

**Evidence** (the real closed loop, m = 32, not the averaged system):

- *Necessity of dither.*  Without dither the collapsed start is trapped forever. Its loss stays at 0.166,
  against 0.007 from the spread start (`results/langevin_tradeoff`).
- *Equalisation at moderate temperature.*  At σ = 0.02, τ = 0.004 (3 seeds,
  `results/langevin_global/table_sigma0.02_tau0.004.md`), collapsed and spread starts reach the same loss
  and tracking error (loss ratio 1.03). Their final measures are 4.2× closer in W than without dither
  (0.40 against 1.69). The common loss (≈ 0.11) is, however, dominated by the regularisation bias: the
  unregularised optimum is 7e-6.
- *Metastability at low temperature.*  At σ = 0.002, τ = 5·10⁻⁴ (`table_sigma0.002_tau0.0005.md`), dither
  frees the collapsed start (0.166 → 0.048), but the starts have not equalised within 900 s.
- The temperature sweep above quantifies the trade-off. Conjecture H.6 is therefore supported
  qualitatively, but the equalisation time can exceed practical horizons at low τ.

---

## 7. Verification map

| statement | what is checked | where |
|---|---|---|
| constants (1), closed forms | ∇²_{(c,s)}φ ≤ H_φ by brute force | `test_transport.py::test_gaussian_constants_are_suprema` |
| A(i) raw rates | lifted raw rates = classical raw rates under the gain map (1e-10) | `::test_identification_with_classical_law` |
| A(i) projected rates, and sharpness | per-kernel projection with rescaled layer matches to 1e-10 with active projections; an unrescaled layer does not | `::test_identification_projected_rates_per_kernel` |
| B with constant ideal and ε = 0 | identity to fp precision, cross ≤ 0, remainder ≤ Lemma B.2, V̇_π ≤ (B) at every sample | `::test_theorem_B_inequality_along_trajectories` |
| B with moving ideal, ε ≠ 0, switches and active projection | the same plus V_*W, at least one switch, projection active at ≥ 10 samples, Dini quotient ≤ (B) | `::test_theorem_B_moving_ideal_nonrepresentable_switches_and_active_projection` |
| B, Dini form | forward-difference quotient of V_OT ≤ (B) | `::test_dini_derivative_of_V_OT_bounded_by_active_branch` |
| B versus labelling | matched error ≪ labelled error, and the bound is 30–45× smaller | `::test_permuted_ideal_labelled_vs_OT`, `exp_transport_permutation.py` |
| C(i), Lemma C.3 | zero switches and exact e^{−σt} decay (explicit Euler, pure leak) | `::test_theorem_C_i_exact_dissipation_pure_leak` |
| C(ii) | induced-coupling bound and the gap case on random configurations | `::test_lemma_C_ii_induced_coupling_bound` |
| C(iv) | an instance satisfying **every** hypothesis, asserted, with the envelope at every sample and actual contraction | `::test_theorem_C_iv_exponential_estimate` |
| C(v), Prop. C.5 | label invariance of the OT leak, and the reshuffling transient of the labelled leak | `exp_ot_sigma_mod.py`, `exp_sigma_mod_transient.py` |
| D | consistency in m (constant amplification 1.68) | `exp_mean_field_m.py` |
| E linearisation | Jacobian of the implemented loop = A(t); QA + AᵀQ = diag(−k_eI, 0) | `::test_theorem_E_linearisation_and_lyapunov_identity` |
| E convergence | matched error decays, labelled error does not (qualitative) | `::test_theorem_E_matched_error_converges_labelled_does_not`, `exp_pe_convergence.py` |
| F | Wasserstein path length ≤ λ_G^{1/2}L₁∫‖e‖ at every sample | `::test_theorem_F_wasserstein_length_bounded_by_integrated_error` |
| G(a), two scales | L̄ at the current weights non-increasing at every sample, even when (H2b) fails | `exp_variable_projection.py` → `table_violating.md` |
| G(b), three scales | velocity ∥ −∇L* and L* decreasing, in a configuration whose time-scale ratios are measured and satisfy (H2b) | `exp_variable_projection.py` → `table_satisfying.md` |
| H(i) identification | averaged drift = mean-field Langevin drift (1e-7) | `test_langevin.py::test_averaged_drift_is_the_mean_field_langevin_drift` |
| H, need for entropy | noiseless symmetric start never separates; dither separates it | `test_langevin.py::test_noiseless_symmetric_initialisation_never_separates_but_dither_does` |
| H(ii) constant | Holley–Stroock ∘ Bakry–Émery formula | `test_langevin.py::test_lsi_lower_bound_formula` |
| H, global behaviour | collapsed versus spread starts with and without dither, two (σ, τ) settings | `exp_langevin_global.py` |
| H, speed–bias trade-off | escape time and loss floor against τ; saddle versus Arrhenius fit | `exp_langevin_tradeoff.py` |

---

## 8. Research program

Each branch is stated with what is proved, what is open, and the first concrete step.

1. **Infinite-horizon averaging for McKean–Vlasov closed loops (Conjecture H.6).**
   - Proved: the averaged system converges globally and exponentially (H(ii)).
   - Open: whether the real loop inherits this.
   - First step: Khasminskii averaging for the mean-field SDE with periodic fast forcing, on finite
     horizons, then extension to infinite horizons using H(ii) as the contraction.

2. **Uniform-in-time propagation of chaos with a plant in the loop.**
   - Theorem D is finite-horizon. Uniform-in-time results exist for mean-field Langevin dynamics [S5] and
     for shallow networks [36].
   - Open: the closed loop, where the plant couples all particles through e.

3. **Noiseless global convergence.**
   - Open: identify a closed-loop separation condition under which the noiseless flow avoids spurious
     stationary points despite the failure of the Chizat–Bach boundary conditions. Candidates are a
     compactified width coordinate and reflecting boundaries.

4. **Unbalanced transport, meaning online structure adaptation.**
   - Treat the weights as masses (Hellinger–Kantorovich / Wasserstein–Fisher–Rao) to get birth and death
     of kernels with Lyapunov guarantees. The conic particle descent of [49] and [55] is the static analogue.

5. **Excitation theory in measure space.**
   - Open: the continuum form of (PE) (E), the trajectory kernel K_T(ϑ,ϑ') = ∫Φ(x(s);ϑ)Φ(x(s);ϑ')ds,
     identifiability of measures up to symmetry, and the relation to [16].

6. **Transport-regularised adaptive laws.**
   - Wasserstein proximal (JKO) discrete-time adaptive laws.
   - The OT leak toward priors with a different number of atoms, using optimal plans with barycentric
     projection.
   - Transfer: transporting a learned measure to a new plant.

7. **Quantisation and the curse of dimensionality (Corollary B.6).**
   - Open: rates in m and n for optimal quantisers of the ideal measure on the lifted space, and adaptive
     laws that track optimal quantisers (the centroidal Voronoi tessellation connection [14]).

8. **Symmetric systems beyond kernels.**
   - The OT quotient applies verbatim to swarms of identical agents, mixtures and multi-head models
     ([88] and [89] are the control-side precedents).

9. **Deep approximators.**
   - Hierarchical (multi-layer) mean-field measures, with a layer-wise OT Lyapunov function.

---

## 9. References

See `LITERATURE_TRANSPORT.md`: [1]–[89] in `review/LITERATURE_refs.md`, and [S1]–[S9]. Key mathematical
references used in proofs: Ambrosio–Gigli–Savaré (2008); Khalil (2002); Kokotović–Khalil–O'Reilly (1986);
Morgan–Narendra (1977); Anderson (1977); Chizat (2022, TMLR); Nitanda–Wu–Suzuki (2022); Chizat–Bach (2018);
Holley–Stroock (1987); Bakry–Émery (1985); Graf–Luschgy (2000); Golse (2013); Pomet–Praly (1992).
