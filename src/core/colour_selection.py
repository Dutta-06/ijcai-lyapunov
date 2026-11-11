"""
Proposition D‴.43 (PARADIGM.md §5d): the colour of the measurement noise selects the kernels — from the trace (a weighted ℓ¹ norm) to the
log-determinant (a support count)
=========================================================================================================================================

A kernel dictionary φ_1, …, φ_d (fixed centres and widths) sampled at N < d points (the design Φ, N × d), the coefficient of atom i the
difference of two adaptive masses, w_i = u_i² − v_i²: gradient descent on (u, v) is the unbalanced Fisher–Rao (Hellinger) law of the
signed kernel measure's two parts, the masses' half of the paradigm's laws on a fixed dictionary (the replicator's mass–amplitude gauge in
flat coordinates). Fitted to y by L = (1/2N)|Φw − y|² with the measured error r + ξ (label noise: ξ independent over the samples, per-step
variance σ², stationary in time with step correlations ρ_d), the parameters reach the manifold M = {Φw = y} and then drift along it at
order v_h = hσ²/N (Proposition D‴.41(iv)(b)): in s = 4v_ht/N, with D_i = u_i² + v_i² = √(w_i² + 4P_i²), P_i = u_iv_i,

    dP/ds = −ω⊙P,   dw/ds = −Π_D(ω⊙w),   ω_i = φ_iᵀ f(K) φ_i,   K = (4/N)Φ diag(D) Φᵀ,   Π_D = I − DΦᵀ(ΦDΦᵀ)⁻¹Φ,

the gradient flow on M of R = tr Ψ(K) (Ψ′ = f, the filter of the draws' correlation, `label_noise.lr_factor`) in the flat metric, which in
(w, P) is diag(1/(4D), 1/D). White draws: f = 1, R = (4/N)Σ_i|φ_i|²D_i, on M as P → 0 the weighted ℓ¹ norm (4/N)Σ_i|φ_i|²|w_i|; AR(1)
draws of correlation time τ_c sampled finely: f = (I + τ_cK)⁻¹, R = (1/τ_c)log det(I + τ_cK) = (1/τ_c)log Σ_{|S|≤N}(4τ_c/N)^{|S|}
det(Φ_SᵀΦ_S)Π_{i∈S}D_i (Cauchy–Binet), concave in D, ≈ min(|w|₀, N)·log(τ_c)/τ_c for large τ_c — the support's size.

This module: the filter weights, the potential and its subset expansion, the reduced flow, the gradient-descent interpolant, weighted basis
pursuit, label-noise gradient descent (paths in parallel), and Gaussian kernel dictionaries.
"""
import itertools

import numpy as np

from . import label_noise as LN


def gram(Phi, D):
    """K = (4/N)Φ diag(D)Φᵀ, the Gauss–Newton matrix in sample space (N × N)."""
    return (4.0 / Phi.shape[0]) * (Phi * D[None, :]) @ Phi.T


def filter_of(lam, h, corr):
    """The filter f(λ) on K's eigenvalues: the Euler scheme's linear-response factor at the step h for the draws `corr`
    (`label_noise.lr_factor`), or, with h = 0, its continuum limit — 1 for fresh draws, 1/(1 + λτ_c) for ('ou', τ_c)."""
    lam = np.clip(np.asarray(lam, float), 0.0, None)
    if h == 0:
        return np.ones_like(lam) if corr[0] == 'fresh' else 1.0 / (1.0 + lam * float(corr[1]))
    return LN.lr_factor(lam, h, corr)


def weights(Phi, D, h, corr):
    """ω_i = φ_iᵀ f(K) φ_i (d,): the coordinate's selection rate; |φ_i|² for white draws."""
    kap, Q = np.linalg.eigh(gram(Phi, D))
    return ((Q.T @ Phi) ** 2 * filter_of(kap, h, corr)[:, None]).sum(axis=0)


def potential(Phi, D, h, corr):
    """R = tr Ψ(K), Ψ′ = f (`label_noise.lr_potential`; in the continuum (1/τ_c)log det(I + τ_cK) for ('ou', τ_c), tr K for fresh)."""
    kap = np.clip(np.linalg.eigvalsh(gram(Phi, D)), 0.0, None)
    if h == 0:
        return float(kap.sum()) if corr[0] == 'fresh' else float(np.log1p(kap * float(corr[1])).sum() / float(corr[1]))
    return float(LN.lr_potential(kap, h, corr).sum())


def det_subsets(Phi, D, tau):
    """det(I + τK) by its subset expansion Σ_{|S|≤N}(4τ/N)^{|S|}det(Φ_SᵀΦ_S)Π_{i∈S}D_i (Cauchy–Binet; small d only)."""
    N, d = Phi.shape
    tot = 1.0
    for k in range(1, min(N, d) + 1):
        for S in itertools.combinations(range(d), k):
            S = list(S)
            tot += (4 * tau / N) ** k * np.linalg.det(Phi[:, S].T @ Phi[:, S]) * np.prod(D[S])
    return tot


def project(Phi, D, q, rcond=1e-13):
    """Π_D q = q − DΦᵀ(ΦDΦᵀ)⁻¹Φq, stably: with B = ΦD^{1/2}, Π_D = D^{1/2}(I − P_B)D^{−1/2}, P_B the orthogonal projector onto range(Bᵀ)
    (a rank-truncated SVD: near a support smaller than N, ΦDΦᵀ is singular to working precision)."""
    sq = np.sqrt(D)
    Vt = np.linalg.svd(Phi * sq[None, :], full_matrices=False)[1:]
    sv, Vt = Vt[0], Vt[1]
    Vr = Vt[sv > rcond * sv.max()]
    z = q / sq
    return sq * (z - Vr.T @ (Vr @ z))


def rhs(s, z, Phi, h, corr):
    """The reduced flow in s = 4v_ht/N on z = (w, log P)."""
    d = Phi.shape[1]
    w, P = z[:d], np.exp(z[d:])
    D = np.maximum(np.sqrt(w * w + 4 * P * P), 1e-300)                         # an atom underflowing to D = 0 stays where it is
    if not np.all(np.isfinite(D)):
        return np.full_like(z, np.nan)                                          # an overshooting trial stage: the integrator rejects it
    om = weights(Phi, D, h, corr)
    return np.r_[-project(Phi, D, om * w), -om]


def reduced_flow(Phi, w0, P0, S, h, corr, t_eval=None, rtol=1e-8, atol=1e-11, method='DOP853', chunks=20, stiff_nfev=20000):
    """The reduced flow from (w0, P0) over [0, S] → (s, w (d × len(s)), P), at t_eval or at S: an explicit integrator (the flow is mostly
    not stiff in (w, log P), and LSODA's finite-difference Jacobians cost 2d evaluations each), over `chunks` pieces, switching to LSODA
    for the rest once a piece needs more than `stiff_nfev` evaluations (a coherent dictionary makes ΦDΦᵀ ill-conditioned, the flow stiff)."""
    from scipy.integrate import solve_ivp
    d = Phi.shape[1]
    t_eval = None if t_eval is None else np.asarray(t_eval, float)
    z, meth, ts, ys = np.r_[w0, np.log(P0)], method, [], []
    edges = np.linspace(0.0, S, chunks + 1)
    for j, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        last = j == chunks - 1
        want = None if t_eval is None else t_eval[(t_eval >= a) & ((t_eval <= b) if last else (t_eval < b))]
        te = None if t_eval is None else np.unique(np.r_[want, b])
        sol = solve_ivp(rhs, (a, b), z, args=(Phi, h, corr), method=meth, rtol=rtol, atol=atol, t_eval=te)
        z = sol.y[:, -1]
        if t_eval is not None:
            keep = np.isin(sol.t, want)
            ts.append(sol.t[keep]); ys.append(sol.y[:, keep])
        if sol.nfev > stiff_nfev:
            meth = 'LSODA'
    if t_eval is None:
        return np.array([S]), z[:d, None], np.exp(z[d:, None])
    y = np.concatenate(ys, axis=1)
    return np.concatenate(ts), y[:d], np.exp(y[d:])


def flow_to_convergence(Phi, w0, P0, S0, h, corr, S_max=1e5, tol_P=1e-8, tol_w=1e-10):
    """The reduced flow run until it has converged: first to S0, then over doubling horizons until max P < tol_P·max|w| and w moved by
    less than tol_w·max|w| over the last horizon, or s reaches S_max → (w, P, s reached, converged)."""
    w, P = np.asarray(w0, float), np.asarray(P0, float)
    s, dS = 0.0, float(S0)
    while True:
        _, W, Pp = reduced_flow(Phi, w, np.maximum(P, 1e-300), dS, h, corr)
        w1, P1 = W[:, -1], Pp[:, -1]
        s += dS
        scale = max(np.abs(w1).max(), 1e-300)
        done = P1.max() < tol_P * scale and np.abs(w1 - w).max() < tol_w * scale
        w, P = w1, P1
        if done or s >= S_max:
            return w, P, s, bool(done)
        dS = min(s, S_max - s)


def uv_of(w, P):
    """(u, v) with u² − v² = w, uv = P (u, v > 0 for P > 0)."""
    u = np.sqrt((np.sqrt(w * w + 4 * P * P) + w) / 2)
    return u, P / u


def gd_interpolant(Phi, y, alpha=None, u=None, v=None, h=None, tol=1e-11, max_steps=2_000_000):
    """The noiseless gradient descent's interpolant from u = v = α (or from (u, v)) → (w, P); the step h, or 0.2/λ_max(K) renewed every
    500 steps (from a small start the coefficients grow and K with them), and never above 0.1/max|g| (so no u_i or v_i changes sign)."""
    N, d = Phi.shape
    u = np.full(d, alpha, float) if u is None else np.array(u, float)
    v = np.full(d, alpha, float) if v is None else np.array(v, float)
    hk = h
    for k in range(max_steps):
        r = Phi @ (u * u - v * v) - y
        if k % 500 == 0:
            if np.abs(r).max() < tol:
                break
            if h is None:
                hk = 0.2 / np.linalg.eigvalsh(gram(Phi, u * u + v * v)).max()
        g = Phi.T @ r / N
        he = min(hk, 0.1 / max(np.abs(g).max(), 1e-300))
        u, v = u - 2 * he * g * u, v + 2 * he * g * v
    return u * u - v * v, u * v


def flow_interpolant(Phi, y, u, v, tol=1e-13, max_iter=200):
    """The noiseless gradient flow's interpolant from (u, v), exactly: along the flow u_i = u_i(0)e^{−s_i/2}, v_i = v_i(0)e^{s_i/2} with
    s = Φᵀμ (μ = 4∫r dt/N), so w_i(μ) = u_i²e^{−s_i} − v_i²e^{s_i}, P_i = u_iv_i conserved, and μ solves Φw(μ) = y — the maximiser of
    the strictly concave −Σ_i(u_i²e^{−s_i} + v_i²e^{s_i}) − yᵀμ (Newton with backtracking) → (w, P). From u = v = α: w = −2α²sinh(s), the
    minimiser of the hyperbolic entropy on {Φw = y} (Woodworth et al.). When y needs atoms whose masses have underflowed (a noisy end
    with fewer than N atoms alive and a residual outside their span) the maximiser is at infinity: the iteration stops at max_iter and
    the caller checks the residual (`end_interpolants`)."""
    u2, v2 = np.asarray(u, float) ** 2, np.asarray(v, float) ** 2
    N = Phi.shape[0]
    mu = np.zeros(N)

    def obj(m_):
        with np.errstate(over='ignore', invalid='ignore'):
            f = -float(np.sum(u2 * np.exp(-(Phi.T @ m_)) + v2 * np.exp(Phi.T @ m_))) - float(y @ m_)
        return f if np.isfinite(f) else -np.inf
    for it in range(max_iter):
        sv = Phi.T @ mu
        w = u2 * np.exp(-sv) - v2 * np.exp(sv)
        F = Phi @ w - y
        if np.abs(F).max() < tol * max(1.0, np.abs(y).max()):
            break
        Hm = (Phi * (u2 * np.exp(-sv) + v2 * np.exp(sv))[None, :]) @ Phi.T       # −∇²(objective), positive semi-definite
        step = np.linalg.lstsq(Hm, F, rcond=None)[0]                               # ascent μ ← μ + step (minimal norm: with fewer
        f0, t = obj(mu), 1.0                                                       # than N atoms alive Hm is singular)
        while obj(mu + t * step) < f0 + 0.25 * t * float(F @ step) and t > 1e-12:
            t *= 0.5
        mu = mu + t * step
    sv = Phi.T @ mu
    with np.errstate(over='ignore', invalid='ignore'):
        return u2 * np.exp(-sv) - v2 * np.exp(sv), np.asarray(u, float) * np.asarray(v, float)


def end_interpolants(Phi, y, U, V, tol=1e-9):
    """The flow interpolants of the rows of (U, V) (paths' ends) → (W (paths × d), converged (paths,)): a row whose interpolant does not
    reach |Φw − y| ≤ tol·max(1, |y|) (it would revive underflowed atoms) is NaN."""
    W, ok = [], []
    for u, v in zip(U, V):
        w = flow_interpolant(Phi, y, u, v)[0]
        c = bool(np.all(np.isfinite(w)) and np.abs(Phi @ w - y).max() <= tol * max(1.0, np.abs(y).max()))
        W.append(w if c else np.full(len(w), np.nan)); ok.append(c)
    return np.array(W), np.array(ok)


def basis_pursuit(Phi, y, c):
    """min Σ_ic_i|w_i| subject to Φw = y (a linear programme in w = w⁺ − w⁻)."""
    from scipy.optimize import linprog
    d = Phi.shape[1]
    res = linprog(np.r_[c, c], A_eq=np.c_[Phi, -Phi], b_eq=y, bounds=(0, None), method='highs')
    return res.x[:d] - res.x[d:]


def label_noise_gd(Phi, y, u0, v0, h, sigma, rho, n_steps, rng, paths, rec_every=None):
    """Gradient descent with label noise on (u, v), `paths` paths in parallel: the labels y + ξ_n, ξ AR(1) with per-step variance σ² and
    step correlation ρ (ρ = 0: drawn afresh), started stationary → (u, v) at the end (paths × d), and the paths' mean of w every
    `rec_every` steps."""
    N = Phi.shape[0]
    u = np.tile(np.asarray(u0, float), (paths, 1)); v = np.tile(np.asarray(v0, float), (paths, 1))
    xi = sigma * rng.standard_normal((paths, N)); a = np.sqrt(1 - rho * rho) * sigma
    rec = []
    for k in range(n_steps + 1):
        if rec_every and k % rec_every == 0:
            rec.append((u * u - v * v).mean(axis=0))
        if k == n_steps:
            break
        g = ((u * u - v * v) @ Phi.T - y[None, :] - xi) @ Phi / N
        u, v = u - 2 * h * g * u, v + 2 * h * g * v
        xi = rho * xi + a * rng.standard_normal((paths, N))
    return u, v, np.array(rec)


def label_noise_heun(Phi, y, u0, v0, h, sigma, rho, n_steps, rng, paths):
    """As `label_noise_gd`, with Heun's scheme, the same draw in both stages (the ODE driven by the draws' piecewise-constant path, whose
    selection is O((hλ)²) of the Euler scheme's, Proposition D‴.41(iv)(d)) → (u, v) at the end (paths × d)."""
    N = Phi.shape[0]
    u = np.tile(np.asarray(u0, float), (paths, 1)); v = np.tile(np.asarray(v0, float), (paths, 1))
    xi = sigma * rng.standard_normal((paths, N)); a = np.sqrt(1 - rho * rho) * sigma
    grad = lambda u, v: ((u * u - v * v) @ Phi.T - y[None, :] - xi) @ Phi / N
    for _ in range(n_steps):
        g1 = grad(u, v)
        u1, v1 = u - 2 * h * g1 * u, v + 2 * h * g1 * v
        g2 = grad(u1, v1)
        u, v = u - h * (g1 * u + g2 * u1), v + h * (g1 * v + g2 * v1)
        xi = rho * xi + a * rng.standard_normal((paths, N))
    return u, v


def sbl_em(Phi, y, lam=1e-10, iters=3000, tol=1e-12):
    """Sparse Bayesian learning (Tipping's EM for the hyperparameters γ, noise variance λ fixed): μ = ΓΦᵀΣ_y⁻¹y, Σ_ii = γ_i − γ_i²φ_iᵀΣ_y⁻¹φ_i,
    γ ← μ² + diag Σ, Σ_y = λI + ΦΓΦᵀ → the posterior mean μ (its log-determinant term log|Σ_y| is R's at γ ∝ |w|, λ = 1/τ_c)."""
    N, d = Phi.shape
    g = np.ones(d)
    mu = np.zeros(d)
    for _ in range(iters):
        Sy = lam * np.eye(N) + (Phi * g[None, :]) @ Phi.T
        A = np.linalg.solve(Sy, np.c_[y, Phi])
        mu_new = g * (Phi.T @ A[:, 0])
        sig = g - g * g * np.einsum('ij,ij->j', Phi, A[:, 1:])
        g = np.maximum(mu_new ** 2 + sig, 1e-300)
        if np.abs(mu_new - mu).max() < tol * max(np.abs(mu_new).max(), 1e-300):
            mu = mu_new
            break
        mu = mu_new
    return mu


def rbf_dictionary(x, centres, width):
    """Gaussian kernels φ_i(x) = exp(−(x − c_i)²/(2 width²)) at the sample points x (N × d), each column scaled to unit norm."""
    Phi = np.exp(-(np.asarray(x, float)[:, None] - np.asarray(centres, float)[None, :]) ** 2 / (2 * width ** 2))
    return Phi / np.linalg.norm(Phi, axis=0)[None, :]
