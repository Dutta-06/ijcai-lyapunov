"""
Reciprocity of the exactly averaged closed loop at order 1/k_e² (PARADIGM.md Theorem I)
=======================================================================================

Proposition H.7(e) shows that the k_e⁻² coefficient U_even + U_odd of the averaged potential can be a first
variation only if the antisymmetric part of its second variation,

    A(ϑ, ϑ') := A_even(ϑ, ϑ') + 2K(ϑ, ϑ'),
    A_even(ϑ, ϑ') = −∫ [ (Φ̃_ϑ'·∇)(r_μᵀΦ̃_ϑ) − (Φ̃_ϑ·∇)(r_μᵀΦ̃_ϑ') ] dρ,
    K(ϑ, ϑ')      = ⟨ Φ̃(x_m; ϑ)ᵀ (d/dt) Φ̃(x_m; ϑ') ⟩_t ,

is a coboundary a(ϑ) − a(ϑ').  Theorem I sharpens and characterises this.

Pairwise criterion.  The feature family contains the zero feature ϑ_z = (a = 0, ξ) (α = 0), and every term of
A(ϑ, ϑ_z) carries Φ̃_{ϑ_z} ≡ 0 or one of its derivatives, so A(ϑ, ϑ_z) = 0 for every ϑ.  A coboundary with this
property is identically zero.  Hence: reciprocal at order k_e⁻² ⟺ A ≡ 0.

Structure.  Writing Dr for the Jacobian of r_μ and [X, Y] := (X·∇)Y − (Y·∇)X for the Lie bracket of vector fields,

    A_even(ϑ, ϑ') = −∫ { Φ̃_ϑᵀ (Dr_μ − Dr_μᵀ) Φ̃_ϑ' + r_μᵀ [Φ̃_ϑ', Φ̃_ϑ] } dρ ,

so the performative (semi-gradient) non-reciprocity is the curl of the residual paired with the features, plus the
residual paired with the Lie bracket of the feature fields.  For Gaussian features Φ̃_ϑ = a φ_ξ(x), ξ := (c, s),
A is a bilinear form in the weights:

    A(ϑ, ϑ') = aᵀ B(ξ, ξ') a',
    B(ξ, ξ') = −∫ (Dr − Drᵀ) φ_ξ φ_ξ' dρ − ∫ r (∇φ_ξ)ᵀ φ_ξ' dρ + ∫ (∇φ_ξ') rᵀ φ_ξ dρ + 2 κ(ξ, ξ') I_n,
    κ(ξ, ξ') := ⟨ φ_ξ(x_m) (d/dt) φ_ξ'(x_m) ⟩_t ,

and A ≡ 0 ⟺ B ≡ 0 on admissible pairs of geometries.

Two characterisations (proofs in PARADIGM.md):
  * n = 1: A ≡ 0 ⟺ r_μ = 0 ρ-a.e. on the orbit.  (B(ξ, ξ') = −((c − c')/s²) ∫ r_μ φ_c φ_c' dρ for equal widths; the
    product of two equal-width Gaussians is a Gaussian centred at (c + c')/2, and the Gaussian convolution of the finite
    signed measure r_μ ρ is entire and injective.)
  * n ≥ 2 at perfect fit (r_μ = 0 on the orbit): B = −∫(Dr_μ − Dr_μᵀ)φφ' dρ + 2κ I (the Lie-bracket terms vanish, the
    curl term survives because r_μ = 0 on the orbit constrains only its tangential derivative), so A ≡ 0 ⟺ K ≡ 0 and
    Dr_μ = Dr_μᵀ ρ-a.e.; A = 2K when r_μ vanishes near the orbit.  K ≡ 0 ⟺ the 1-current of the oriented orbit vanishes,
    T_𝒪 := ∫₀^{T₀} δ_{x_m(t)} ẋ_m(t) dt = 0.  (For equal widths κ(ξ, ξ') = −(1/(2s²)) e^{−|c−c'|²/(4s²)} (c − c')·F_s(m),
    m := (c + c')/2, F_s(m) := ⟨ e^{−|x_m − m|²/s²} ẋ_m ⟩_t = (G_s ∗ T_𝒪)(m) up to a constant, and F_s ≡ 0 on an open
    set forces T_𝒪 = 0.)  For n = 2, T_𝒪 = 0 ⟺ the winding-number function of the orbit vanishes a.e.
  * Nowhere a gradient flow: A_ν is affine in ν, and along (1 − δ)μ + δ δ_{ϑ'} (ϑ' the zero feature, or a feature with
    a' ≠ 0) A_ν ≢ 0 for all small δ > 0, so no open set of measures is reciprocal, for any n.
"""
import numpy as np
from .averaging import features, circulation_kernel, second_variation_formula


# ---------------------------------------------------------------- Gaussian geometry on the orbit
def gauss_phi_grad(X, c, s):
    """φ_ξ(X) and ∇_xφ_ξ(X) for ξ = (c, s): returns (N,), (N, n)."""
    D = X - c[None, :]
    phi = np.exp(-np.sum(D ** 2, axis=1) / (2.0 * s ** 2))
    return phi, -(phi[:, None] * D) / s ** 2


def residual(f_vec, jac_vec, fm, X):
    """r_μ = f − f̂_μ and its Jacobian Dr (N, n, n) with Dr[t, k, l] = ∂r_k/∂x_l on the sample points X."""
    return f_vec(X) - fm.f_hat(X), jac_vec(X) - fm.jac(X)


def kappa(orbit, c, s, c_p, s_p):
    """κ(ξ, ξ') = ⟨φ_ξ(x_m) (d/dt)φ_ξ'(x_m)⟩_t  (antisymmetric; zero for n = 1)."""
    phi, _ = gauss_phi_grad(orbit.xm, c, s)
    _, gp = gauss_phi_grad(orbit.xm, c_p, s_p)
    return float(np.mean(phi * np.einsum('tn,tn->t', gp, orbit.xmd)))


def current_F(orbit, m, s):
    """F_s(m) := ⟨exp(−|x_m(t) − m|²/s²) ẋ_m(t)⟩_t, the Gaussian-smoothed 1-current of the oriented orbit at m."""
    w = np.exp(-np.sum((orbit.xm - m[None, :]) ** 2, axis=1) / s ** 2)
    return np.einsum('t,tn->n', w, orbit.xmd) / orbit.N


def kappa_closed_form(orbit, c, s, c_p):
    """Equal widths: κ(ξ, ξ') = −(1/(2s²)) e^{−|c−c'|²/(4s²)} (c − c')·F_s((c + c')/2)."""
    m = 0.5 * (c + c_p)
    return float(-(1.0 / (2.0 * s ** 2)) * np.exp(-np.sum((c - c_p) ** 2) / (4.0 * s ** 2)) * (c - c_p) @ current_F(orbit, m, s))


# ---------------------------------------------------------------- the bilinear form B and the kernel A
def B_matrix(f_vec, jac_vec, fm, orbit, c, s, c_p, s_p):
    """B(ξ, ξ') ∈ ℝ^{n×n}, with A(ϑ, ϑ') = aᵀ B(ξ, ξ') a'."""
    X = orbit.xm
    r, Jr = residual(f_vec, jac_vec, fm, X)
    phi, g = gauss_phi_grad(X, c, s)
    phi_p, g_p = gauss_phi_grad(X, c_p, s_p)
    n = X.shape[1]
    curl = np.einsum('t,tkl->kl', phi * phi_p, Jr - np.transpose(Jr, (0, 2, 1))) / orbit.N
    lie1 = np.einsum('tk,tl->kl', r * phi_p[:, None], g) / orbit.N          # ∫ r (∇φ)ᵀ φ' dρ
    lie2 = np.einsum('tk,tl->kl', g_p * phi[:, None], r) / orbit.N          # ∫ (∇φ') rᵀ φ dρ
    return -curl - lie1 + lie2 + 2.0 * kappa(orbit, c, s, c_p, s_p) * np.eye(n)


def A_kernel(f_vec, jac_vec, fm, orbit, theta, theta_p):
    """A(ϑ, ϑ') for two unconstrained lifted particles (through the saturation maps of fm.net)."""
    a, c, s = fm.net.physical(np.atleast_2d(theta)); a_p, c_p, s_p = fm.net.physical(np.atleast_2d(theta_p))
    B = B_matrix(f_vec, jac_vec, fm, orbit, c[0], s[0], c_p[0], s_p[0])
    return float(a[0] @ B @ a_p[0])


def A_from_second_variation(f_vec, jac_vec, fm, orbit, Th):
    """The antisymmetric part S − Sᵀ of the k_e⁻² second-variation coefficient from averaging.second_variation_formula,
    on the rows of Th (reference implementation for the tests)."""
    S_e, K = second_variation_formula(f_vec, jac_vec, fm, orbit, np.atleast_2d(Th))
    S = S_e + K
    return S - S.T


def cyclic_sum(A):
    """A(ϑ₀,ϑ₁) + A(ϑ₁,ϑ₂) + A(ϑ₂,ϑ₀) for a 3×3 antisymmetric matrix of pairwise values."""
    return float(A[0, 1] + A[1, 2] + A[2, 0])


def zero_feature(net, xi_unconstrained):
    """The zero-feature particle ϑ_z = (α = 0, geometry) in unconstrained coordinates."""
    n = net.n
    th = np.zeros(2 * n + 1)
    th[n:] = xi_unconstrained
    return th


# ---------------------------------------------------------------- the one-dimensional closed form
def A_even_1d_closed_form(f_vec, fm, orbit, a, c, s, a_p, c_p):
    """n = 1, equal widths: A(ϑ, ϑ') = −a a' ((c − c')/s²) ∫ r_μ φ_c φ_c' dρ."""
    X = orbit.xm
    r = (f_vec(X) - fm.f_hat(X))[:, 0]
    phi = np.exp(-(X[:, 0] - c) ** 2 / (2 * s ** 2)); phi_p = np.exp(-(X[:, 0] - c_p) ** 2 / (2 * s ** 2))
    return float(-a * a_p * (c - c_p) / s ** 2 * np.mean(r * phi * phi_p))


def positivity_witness_1d(orbit, a0, c0, a1, c1, abar, cbar, s):
    """
    f ≡ 0, μ = δ_{ϑ̄} with ϑ̄ = (ā, c̄, s), equal widths: r_μ = −ā φ_c̄ and
        A(ϑ₀, ϑ₁) = a₀ a₁ ā ((c₀ − c₁)/s²) ∫ φ_c̄ φ_c₀ φ_c₁ dρ = a₀ a₁ ā ((c₀ − c₁)/s²) · I,   I > 0,
    because the integrand is strictly positive.  Returns (A, I).
    """
    X = orbit.xm[:, 0]
    I = float(np.mean(np.exp(-((X - cbar) ** 2 + (X - c0) ** 2 + (X - c1) ** 2) / (2 * s ** 2))))
    return a0 * a1 * abar * (c0 - c1) / s ** 2 * I, I
