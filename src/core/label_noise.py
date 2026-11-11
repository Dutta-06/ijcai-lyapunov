"""
Proposition D‴.41 (PARADIGM.md §5d): measurement noise on the gauge — the selection is set by the noise's per-update variance
==========================================================================================================================

The measured error r + ξ drives both laws of (𝔉), so the noise enters along the descent's own vector fields: the drift is
⟨rU⟩_ρ = (1/N)Σ_x r(x)U_x and the noise (1/N)Σ_x ξ(x)U_x, with U_x = mM∇f̂(x), M = diag(G, γ/β) on (Θ, log β) projected on Σβ = m:

    U_x = (β_iG∇_ϑΦ̃(x; ϑ_i),  γ(Φ̃(x; ϑ_i) − f̂(x)))_i.

On the manifold M of exact fits (r ≡ 0) the fields are orthogonal to M in the flow's metric, and the limit map's invariance makes the
Wong–Zakai term cancel the Itô drift: with Stratonovich (or any smooth, coloured) noise the gauge does not move at order v, with noise
drawn afresh at each update it moves by −(v/2)∂Φ⟨(U·∇)U⟩_ρ, which on the masses in the slow-replicator limit is the replicator flow of

    ū_k = −(v/2)β_kG_α⟨Φ_k²⟩_ρ(1 − (a_k/a_max)⁴) = −∂_{β_k}Ψ_ℓ,   Ψ_ℓ = (v/4)Σ_kG_α⟨Φ_k²⟩_ρ(β_ka′_k)²,

the amplitudes' part of (v/4)tr(G∇²J) (the rest of the trace is constant on M); its minimiser on the simplex is interior.

A correlated noise sequence (step correlations ρ_d) selects through the per-step second moment of the whole increment h(F + ⟨Uξ⟩_ρ),
not of the noise alone: the normal deviation n responds to the past draws and, at linear response, the pair (i, l) of normal modes (decay
rates λ_i, λ_l, a = 1 − hλ) is weighted by

    f_il = 1 − h(λ_iR(a_i) + λ_lR(a_l)) + h²λ_iλ_l(1 + G(a_i) + G(a_l))/(1 − a_ia_l),   G(a) = Σ_{d≥1}ρ_da^d,  R(a) = G(a)/a,

while the cross term ⟨∇ψ·n ξ⟩ vanishes by the same antisymmetry as the Wong–Zakai cancellation (tr(R(A)M·B) = 0, R(A)M symmetric). The
per-step-variance law (1 − ρ)/(1 + ρ) is the case f ≡ 1 (λ_maxτ_c ≪ 1); for AR(1) noise of correlation time τ_c sampled finely,
f_ii → 1/(1 + λ_iτ_c), and on a flat-metric gradient flow with label noise (the covariance ∝ the Hessian H on M) the selected potential
becomes tr Ψ(H), Ψ′ = f: the trace for white noise, (1/τ_c) log det(I + τ_cH) in the continuum limit. The log-determinant is the
AR(1) (Markov, white-innovation) noise's: a smooth noise of the same correlation time (Gaussian step correlations, `corr_sums`'
('gauss', k)) has f_ii ≈ (λτ_c)⁻² for λτ_c ≫ 1 and a bounded Ψ (referee round 30).

This module: the fields, the Wong–Zakai term (for any node covariance), the closed-form selection, Ψ_ℓ and its minimiser, the per-step
factor of AR(1) noise, the correlation sums, the linear-response factor and the filtered covariance, its potential.
"""
import numpy as np


def _sig(x):
    return 1.0 / (1.0 + np.exp(-x))


def fields(loop, Th, beta):
    """U_x for every quadrature node in the coordinates (Θ, log β): (UT (N, m, p), Ul (N, m)); `loop` a `ReplicatorCohesionLoop`."""
    base = loop.base
    net = base.net
    a, c, s, diff, phi = base.features(Th)
    da = 1.0 - (a / net.a_max) ** 2
    dc = 1.0 - (c / net.c_max) ** 2
    sk = _sig(Th[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
    gT = np.empty((base.N,) + Th.shape)
    gT[:, :, 0] = phi * da
    gT[:, :, 1] = a * dc / s ** 2 * phi * diff
    gT[:, :, 2] = a * ds / s ** 3 * phi * diff ** 2
    UT = beta[None, :, None] * loop.G[None, None, :] * gT
    Phit = phi * a[None, :]
    fh = Phit @ beta / Th.shape[0]
    Ul = loop.gamma * (Phit - fh[:, None])
    return UT, Ul


def selection(loop, Th, beta, v=1.0, coord='alpha'):
    """The closed form on M in the slow-replicator limit: ū_k = −(v/2)β_kG_α⟨Φ_k²⟩_ρ(a′_k² − a_ka″_k) = −(v/2)β_kG_α⟨Φ_k²⟩_ρ(1 − q_k²),
    q_k = (a_k/a_max)², for an Euler step in the network's amplitude parameter α (a = a_max tanh(α/a_max)); with coord='a' (an Euler step in
    the amplitude a itself) the update's own Itô term adds −aa″: (a′_k² − 2a_ka″_k) = (1 − q_k)(1 + 3q_k). (m,)"""
    a, c, s, diff, phi = loop.base.features(Th)
    q = (a / loop.base.net.a_max) ** 2
    f = (1 - q ** 2) if coord == 'alpha' else (1 - q) * (1 + 3 * q)
    return -0.5 * v * beta * loop.G[0] * np.mean(phi ** 2, axis=0) * f


def psi(W, P, W_min, v=1.0, m=3):
    """Ψ_ℓ at the masses W (Σ W = 1) on M: (v/4)Σ_kP_k(β_ka′_k)², β_ka′_k = m(W_k − W_min,k²/W_k), P_k = G_α⟨Φ_k²⟩_ρ."""
    W = np.asarray(W, float)
    return 0.25 * v * float(np.sum(np.asarray(P, float) * (m * (W - np.asarray(W_min, float) ** 2 / W)) ** 2))


def w_label(P, W_min, coord='alpha'):
    """The minimiser W_ℓ of Ψ_ℓ on the simplex: P_k(W_k − W_min,k²/W_k)(1 + W_min,k²/W_k²) equal (increasing in W_k from 0 at W_min,k),
    so W_ℓ,k > W_min,k for every k when Σ_kW_min,k < 1; coord='a' (`selection`): P_k(W_k − W_min,k²/W_k)(1 + 3W_min,k²/W_k²) equal, the
    minimiser of the strictly convex (v/2)Σ_kG_αP_k(β_k²/2 + 2c_k log β_k + 3c_k²/(2β_k²)), c_k = (mA_k/a_max)²."""
    from scipy.optimize import brentq
    P = np.asarray(P, float); Wm = np.asarray(W_min, float)
    t = 1.0 if coord == 'alpha' else 3.0
    g = lambda W, k: P[k] * (W - Wm[k] ** 2 / W) * (1 + t * Wm[k] ** 2 / W ** 2)
    W_of = lambda mu, k: brentq(lambda W: g(W, k) - mu, Wm[k], 1.0 + 10 * Wm[k])
    hi = min(g(1.0, k) for k in range(len(P)))
    mu = brentq(lambda mu: sum(W_of(mu, k) for k in range(len(P))) - 1.0, 1e-12 * hi, hi)
    return np.array([W_of(mu, k) for k in range(len(P))])


def ar1_factor(rho):
    """The per-step variance of AR(1) noise (successive draws correlated by ρ) over that of noise drawn afresh, at the same intensity:
    (1 − ρ)/(1 + ρ) — the ratio of their selections when hλ_max ≪ 1 − ρ (λ_maxτ_c ≪ 1); beyond, times the linear-response factors
    (`lr_factor`, `filtered_covariance`)."""
    return (1.0 - rho) / (1.0 + rho)


def corr_sums(a, corr):
    """(G(a), R(a)) = (Σ_{d≥1}ρ_da^d, Σ_{d≥1}ρ_da^{d−1}) for the noise's step correlations ρ_d: corr = ('fresh',), ('ar1', ρ) (ρ_d = ρ^d),
    ('block', k) (ρ_d = (1 − d/k)₊, a draw kept k steps, averaged over the phase) or ('gauss', k) (ρ_d = exp(−½(d/k)²), a smooth noise of
    correlation time kh, the sums truncated at d = 9k); `a` a number, an array (elementwise) or a square matrix (the matrix functions;
    for 'gauss' through its eigendecomposition)."""
    a = np.asarray(a, float)
    I = np.eye(a.shape[0]) if a.ndim == 2 else 1.0
    mul = (lambda x, y: x @ y) if a.ndim == 2 else (lambda x, y: x * y)
    if corr[0] == 'fresh':
        R = 0.0 * a
    elif corr[0] == 'ar1':
        rho = float(corr[1])
        R = rho * (np.linalg.inv(I - rho * a) if a.ndim == 2 else 1.0 / (1.0 - rho * a))
    elif corr[0] == 'block':
        k = int(corr[1]); R = 0.0 * a; pw = I + 0.0 * a
        for d in range(1, k):
            R = R + (1.0 - d / k) * pw; pw = mul(pw, a)
    elif corr[0] == 'gauss':
        k = float(corr[1]); d = np.arange(1, int(np.ceil(9 * k)) + 1)
        rd = np.exp(-0.5 * (d / k) ** 2)
        sums = lambda z: (rd[None, :] * np.power(z.reshape(-1, 1).astype(complex), d[None, :] - 1)).sum(axis=1).reshape(z.shape)
        if a.ndim == 2:
            ev, V = np.linalg.eig(a)
            R = ((V * sums(ev)[None, :]) @ np.linalg.inv(V)).real
        else:
            R = sums(np.atleast_1d(a)).real.reshape(a.shape)
    else:
        raise ValueError(corr)
    return mul(a, R), R


def lr_factor(lam, h, corr):
    """f_ii: the per-step second moment of the increment along a normal mode of decay rate λ (an Euler step h) over that of the noise, at
    linear response: 1 − 2hλR(a) + h²λ²(1 + 2G(a))/(1 − a²), a = 1 − hλ (`corr_sums`); 1 + hλ/(2 − hλ) for fresh noise, → 1/(1 + λτ_c) for
    AR(1) noise of correlation time τ_c = −h/log ρ as h → 0."""
    lam = np.asarray(lam, float)
    a = 1.0 - h * lam
    G, R = corr_sums(a, corr)
    return 1.0 - 2 * h * lam * R + h * lam * (1.0 + 2 * G) / (2.0 - h * lam)      # (hλ)²/(1 − a²), stable as hλ → 0


def filtered_covariance(A, P, Sig, h, corr):
    """Σ_f, the per-step second moment of the increment h(F + ⟨Uξ⟩_ρ) near M over h²(the noise's per-step variance), at linear response,
    for the covariance Sig of the noise's fields (normal to M): in an orthonormal basis of the normal space range(I − P) (A the drift's
    Jacobian at the point of M, P the projection onto TM along range A, B = I + hA),
        Σ_f = S + h(A R(B)S + S R(B)ᵀAᵀ) + h²A X Aᵀ,   X − BXBᵀ = S + G(B)S + SG(B)ᵀ,   S = (I − P)Sig(I − P)ᵀ,
    the pair (i, l) of normal modes weighted by f_il (the module's docstring). Katzenberger's drift of Σ_f, times the per-step variance's
    ratio to the fresh noise's, is the scheme's selection at order v."""
    from scipy.linalg import solve_discrete_lyapunov
    n = A.shape[0]
    k = int(round(np.trace(P)))
    Q = np.eye(n) - P
    Nb = np.linalg.svd(Q)[0][:, :n - k]
    AN = Nb.T @ A @ Nb
    S = Nb.T @ Q @ Sig @ Q.T @ Nb
    B = np.eye(n - k) + h * AN
    G, R = corr_sums(B, corr)
    X = solve_discrete_lyapunov(B, S + G @ S + S @ G.T)
    Sf = S + h * (AN @ R @ S + S @ R.T @ AN.T) + h * h * AN @ X @ AN.T
    return Nb @ Sf @ Nb.T


def lr_potential(lam, h, corr, nodes=40):
    """Ψ(λ) = ∫_0^λ f_ii (Gauss–Legendre): on a flat-metric gradient flow the scheme's selection descends tr Ψ(H) along M (λ for white
    noise in the limit h → 0, log(1 + λτ_c)/τ_c for finely sampled AR(1) noise)."""
    x, wq = np.polynomial.legendre.leggauss(nodes)
    lam = np.atleast_1d(np.asarray(lam, float))
    out = np.array([0.5 * l * np.sum(wq * lr_factor(0.5 * l * (x + 1), h, corr)) for l in lam])
    return out


def flat_fields(loop, x, beta_of, m=3):
    """U_x in the flat coordinates (Θ (m·p), log(β_j/β_m), j < m) of `exp_gauge_drift.instance`: (N, m·p + m − 1)."""
    Th = x[:m * 3].reshape(m, 3)
    UT, Ul = fields(loop, Th, beta_of(x[m * 3:]))
    return np.concatenate([UT.reshape(UT.shape[0], -1), Ul[:, :m - 1] - Ul[:, m - 1:]], axis=1)


def field_derivatives(loop, x, beta_of, h=1e-4, m=3):
    """D[j, i] = ∂_{U_j}U_i at the flat point x (central differences), (N, N, m·p + m − 1)."""
    U0 = flat_fields(loop, x, beta_of, m)
    return np.array([(flat_fields(loop, x + h * U0[j], beta_of, m) - flat_fields(loop, x - h * U0[j], beta_of, m)) / (2 * h)
                     for j in range(U0.shape[0])])


def wong_zakai(loop, x, beta_of, h=1e-4, m=3, weights=None, cov=None):
    """⟨w(U·∇)U⟩_ρ = (1/N)Σ_x w_x ∂_{U_x}U_x at the flat point x (central differences), (m·p + m − 1,); w ≡ 1 unless `weights` (the noise's
    variance profile over the nodes, heteroscedastic noise) is given; with `cov` (N × N, the noise's covariance over the nodes, symmetric)
    (1/N)Σ_{ij}C_ij∂_{U_j}U_i."""
    if cov is not None:
        return np.einsum('ij,jik->k', np.asarray(cov, float), field_derivatives(loop, x, beta_of, h, m)) / len(cov)
    U0 = flat_fields(loop, x, beta_of, m)
    wt = np.ones(U0.shape[0]) if weights is None else np.asarray(weights, float)
    out = np.zeros(U0.shape[1])
    for i in range(U0.shape[0]):
        out += wt[i] * (flat_fields(loop, x + h * U0[i], beta_of, m)[i] - flat_fields(loop, x - h * U0[i], beta_of, m)[i]) / (2 * h)
    return out / U0.shape[0]
