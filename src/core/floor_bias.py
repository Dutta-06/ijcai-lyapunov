"""
Proposition D‴.24 (PARADIGM.md §5d): the finite-size floor is the leak's bias at the sampled weights
====================================================================================================

Let the target admit, for the weights w near q, an exact weighted K₀-atom representation ϑ^e(w) — f = Σ_kw_kΦ̃(·;ϑ^e_k(w)), a critical point
of the unregularised weighted loss L_w(Θ) := ½‖f − Σ_kw_kΦ̃(·;ϑ_k)‖² with an invertible Gauss–Newton block — and let the prior be matching,
ϑ⁰ = ϑ^e(q). Under the nearest-atom leak toward that fixed prior the weighted stationary point Θ*(w) of Proposition D‴.5 is biased away from
ϑ^e(w): in the whitened weighted coordinates y_k = √w_kG^{−1/2}(ϑ*_k(w) − ϑ^e_k(w)), to first order in the prior's offset
x⁰_k(w) := √w_kG^{−1/2}(ϑ⁰_k − ϑ^e_k(w)),

    y(w) = σ₀(H_L(w) + σ₀I)⁻¹x⁰(w),    f̂_{Θ*(w)} − f = Ĵ(w)y(w),    Ĵ(w) := [√w_k∇_ϑΦ̃(·;ϑ^e_k(w))G^{1/2}]_k,  ĴᵀĴ = H_L,

and x⁰(q) = 0, so the function-level weight-sensitivity of Propositions D‴.5 and D‴.9 is the weight-derivative of the bias,
F = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ at q, and its floor decomposes on the eigenvectors v_j of H_L (eigenvalues h_j):

    E‖f̂_m − f̂_M‖² = (1/m − 1/M)·Σ_j h_jσ₀²/(h_j + σ₀)²·z_jᵀC̃z_j,    z := Vᵀ∂_wx⁰,    h_jσ₀²/(h_j + σ₀)² ≤ σ₀/4,

each mode weighted by its share of the bias times its excitation; the weight peaks at h_j = σ₀, so the floor is at most
(1/m − 1/M)(σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ) whatever the excitation spectrum. This module computes ϑ^e(w) (Newton on the unregularised weighted
gradient), ∂_wϑ^e along zero-sum directions (central differences), Ĵ, H_L and the decomposition.
"""
import numpy as np

from .ot_basin import feature_tables


def unregularised_grad(loop, w, Th):
    """∇_ϑ𝒱_L[μ_w](ϑ_k) = −⟨r_w∇_ϑΦ̃(·;ϑ_k)⟩ with r_w = f − Σ_kw_kΦ̃(·;ϑ_k): the weighted loss's gradient with the weights divided out, (K, 3)."""
    a, c, s, diff, phi = loop.features(Th)
    return loop.grad_V_L(loop.F - phi @ (np.asarray(w, float) * a), Th)


def exact_weighted_atoms(loop, w, Th0, iters=60, h=1e-6, tol=1e-13, n_cont=4, w0=None, fit_tol=1e-10):
    """ϑ^e(w): Newton (finite-difference Jacobian, backtracking) on the unregularised weighted gradient from Th0, the weights moved from w0
    (default: uniform) to w in n_cont steps; if the critical point reached is not exact (relative fit error above fit_tol — near the
    saturation a coarse continuation can jump to another critical point) the continuation is refined to 4n_cont and 16n_cont steps.
    Returns (ϑ^e(w), sup-norm residual, relative fit error) of the last attempt."""
    for nc in (n_cont, 4 * n_cont, 16 * n_cont):
        out = _exact_weighted_atoms(loop, w, Th0, iters, h, tol, nc, w0)
        if out[2] < fit_tol:
            break
    return out


def _exact_weighted_atoms(loop, w, Th0, iters, h, tol, n_cont, w0):
    K = Th0.shape[0]; d = 3 * K
    w = np.asarray(w, float)
    w0 = np.ones(K) / K if w0 is None else np.asarray(w0, float)
    v = np.asarray(Th0, float).ravel().copy()
    for j in range(n_cont):
        wj = w0 + (w - w0) * (j + 1) / n_cont
        Fn = lambda v, wj=wj: unregularised_grad(loop, wj, v.reshape(K, 3)).ravel()
        for _ in range(iters):
            g = Fn(v); gn = np.abs(g).max()
            if gn < tol:
                break
            J = np.zeros((d, d))
            for i in range(d):
                e = np.zeros(d); e[i] = h
                J[:, i] = (Fn(v + e) - Fn(v - e)) / (2 * h)
            step = np.linalg.solve(J, g); lam = 1.0
            while lam > 1e-8:
                vn = v - lam * step
                if np.abs(Fn(vn)).max() < gn:
                    break
                lam *= 0.5
            v = vn
    Th = v.reshape(K, 3)
    a, c, s, diff, phi = loop.features(Th)
    r = loop.F - phi @ (w * a)
    return Th, float(np.abs(unregularised_grad(loop, w, Th)).max()), float(np.sqrt(np.mean(r ** 2) / np.mean(loop.F ** 2)))


def whitened_jacobian(loop, w, atoms):
    """Ĵ (N, 3K): the columns √w_k∇_ϑΦ̃(·;ϑ_k)G^{1/2}, k-major (the ordering of adaptive_prior.whitened); ĴᵀĴ/N = H_L at an exact representation."""
    Phi, dPhi = feature_tables(loop, np.asarray(atoms, float))
    return (np.sqrt(np.asarray(w, float))[None, :, None] * dPhi * np.sqrt(loop.G)[None, None, :]).reshape(loop.N, -1)


def bias_floor_law(loop, q, atoms, sigma0, h=1e-5, directions=None):
    """At a matching prior (atoms = ϑ^e(q) = the prior): ∂_wϑ^e along zero-sum directions spanning the zero-sum plane (default e_j − e_K; the
    covariance of the coordinates is Vd⁺C̃Vd⁺ᵀ, so the law does not depend on the basis), ∂_wx⁰, H_L = ĴᵀĴ/N with its
    eigenpairs, the bias formula F_bias = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ (N, K − 1), the floor's per-mode terms h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j (the law per
    unit of 1/m − 1/M, with C̃ the multinomial covariance restricted to the directions) and the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ)."""
    q = np.asarray(q, float); K = q.size; G = loop.G
    atoms = np.asarray(atoms, float)
    Vd = np.array([np.eye(K)[j] - np.eye(K)[K - 1] for j in range(K - 1)]).T if directions is None else np.asarray(directions, float)
    dthe = []
    for col in range(Vd.shape[1]):
        tp, rp, _ = exact_weighted_atoms(loop, q + h * Vd[:, col], atoms, n_cont=1, w0=q)
        tm, rm, _ = exact_weighted_atoms(loop, q - h * Vd[:, col], atoms, n_cont=1, w0=q)
        assert rp < 1e-10 and rm < 1e-10
        dthe.append((tp - tm) / (2 * h))
    dx0 = np.column_stack([(-np.sqrt(q)[:, None] * d / np.sqrt(G)[None, :]).ravel() for d in dthe])       # ∂_v x⁰ at q
    Jh = whitened_jacobian(loop, q, atoms)
    HL = Jh.T @ Jh / loop.N
    hs, V = np.linalg.eigh(0.5 * (HL + HL.T))
    F_bias = Jh @ (sigma0 * np.linalg.solve(HL + sigma0 * np.eye(HL.shape[0]), dx0))
    Ct = np.diag(q) - np.outer(q, q)
    Vp = np.linalg.pinv(Vd)                                    # δw = Vd·c on the zero-sum plane: c = Vd⁺δw, Cov(c) = Vd⁺C̃Vd⁺ᵀ (any spanning basis;
    C2 = Vp @ Ct @ Vp.T                                        # for the default e_j − e_K it is C̃ restricted to (w₁, …, w_{K−1}))
    z = V.T @ dx0
    terms = hs * (sigma0 / (hs + sigma0)) ** 2 * np.einsum('ja,ab,jb->j', z, C2, z)
    bound = 0.25 * sigma0 * float(np.trace(dx0 @ C2 @ dx0.T))
    return dict(dtheta=np.array(dthe), dx0=dx0, HL=HL, h=hs, V=V, F_bias=F_bias, C2=C2, z=z, terms=terms, law=float(terms.sum()), bound=bound)
