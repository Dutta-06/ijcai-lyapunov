"""
Proposition D‴.20 (PARADIGM.md §5d): the entropic OT leak — a soft mass constraint with a critical excess
=========================================================================================================

The barycentric OT leak of Theorem C′ / Proposition D‴.12 uses a vertex optimal plan, which switches discontinuously.  Its entropic
regularisation replaces the plan by the Sinkhorn plan π^ε from μ_Θ = (1/m)Σδ_{ϑ_i} to ν = Σ_jq_jδ_{ϑ⁰_j} for the cost C_{ij} = ‖ϑ_i − ϑ⁰_j‖²_{G⁻¹},

    π^ε := argmin{⟨C, π⟩ + ε KL(π | μ⊗ν) : π1 = 1/m, πᵀ1 = q},        T^ε_i := mΣ_jπ^ε_{ij}ϑ⁰_j,

and the leak −σ₀(ϑ_i − T^ε_i) is the G-gradient flow of J^ε_m := mL̄(μ_Θ) + (σ₀m/2)OT_ε(μ_Θ, ν) (the envelope theorem: ∂OT_ε/∂ϑ_i =
Σ_jπ^ε_{ij}∂C_{ij}/∂ϑ_i = (2/m)G⁻¹(ϑ_i − T^ε_i)), smooth in Θ.  It interpolates between the OT leak (ε → 0) and the σ-modification toward
the barycentre Σ_jq_jϑ⁰_j of ν (ε → ∞, π^ε → μ⊗ν).  This module implements the entropic plan (Newton's method on the semi-dual, gauge-fixed,
for a few columns — the prior's atoms — and log-domain Sinkhorn otherwise), the entropic leak loop (a subclass of
`prior_design.MixturePriorLoop`), the separation block of a collapsed class, and the critical regularisation of Proposition D‴.20.
"""
import warnings

import numpy as np

from .chaos_learning import _sig
from .prior_design import MixturePriorLoop
from .ot_leak import _cost_matrix


def _logsumexp(A, axis):
    mx = np.max(A, axis=axis, keepdims=True)
    return np.squeeze(mx, axis=axis) + np.log(np.sum(np.exp(A - mx), axis=axis))


def _semidual_newton(C, a, b, eps, g0=None, tol=1e-13, max_iter=200):
    """The entropic plan by Newton's method on the concave semi-dual h(g) = Σ_ia_if_i(g) + Σ_jb_jg_j, f_i(g) = −ε log Σ_jb_j e^{(g_j − C_ij)/ε}
    (the row constraint exact for every g); the gauge g_1 = 0 is fixed. Quadratically convergent for any ε: used when the number of columns is
    small (the prior's atoms)."""
    n = C.shape[1]
    g = np.zeros(n) if g0 is None else np.array(g0, float) - g0[0]
    lb = np.log(b)
    def plan(g):
        A = (g[None, :] - C) / eps + lb[None, :]
        mx = A.max(axis=1, keepdims=True)
        P = np.exp(A - mx); s = P.sum(axis=1, keepdims=True)
        cond = P / s                                                      # π(j|i)
        f = -eps * (mx[:, 0] + np.log(s[:, 0]))
        return a[:, None] * cond, cond, f
    for it in range(max_iter):
        pi, cond, f = plan(g)
        grad = b - pi.sum(axis=0)                                          # ∂h/∂g_j
        if np.max(np.abs(grad)) < tol:
            break
        H = (pi.T @ cond - np.diag(pi.sum(axis=0))) / eps                  # ∂²h/∂g∂g (negative semidefinite, kernel = constants)
        Hr, gr = H[1:, 1:], grad[1:]
        step = np.zeros(n)
        try:
            step[1:] = np.linalg.solve(Hr, -gr)                            # Newton ascent direction
            if not np.all(np.isfinite(step)):
                raise np.linalg.LinAlgError
        except np.linalg.LinAlgError:                                      # a nearly deterministic plan: minimum-norm step in the curved directions
            step[1:] = np.linalg.lstsq(Hr, -gr, rcond=1e-12)[0]
        hval = lambda gg: float(np.sum(a * plan(gg)[2]) + np.sum(b * gg))
        h0, lam = hval(g), 1.0
        while lam > 1e-10 and hval(g + lam * step) < h0 - 1e-15 * max(1.0, abs(h0)):
            lam *= 0.5
        g = g + lam * step
    pi, cond, f = plan(g)
    return pi, f, g, it + 1


def sinkhorn_plan(C, a, b, eps, tol=1e-13, max_iter=100000, f0=None, g0=None):
    """The entropic plan π = diag(a)exp((f ⊕ g − C)/ε)diag(b) with marginals (a, b); returns (π, f, g, iterations), f, g the dual potentials
    (warm-startable). For at most eight columns (the prior's atoms) Newton's method on the semi-dual, with ε-scaling from the cost range down
    when no warm start is given (robust for small ε); otherwise log-domain Sinkhorn iterations."""
    if C.shape[1] <= 8:
        a = np.asarray(a, float); b = np.asarray(b, float)
        if g0 is None:
            e_hi = max(eps, float(C.max() - C.min()))
            g = None; e_cur = e_hi; its = 0
            while e_cur > eps * (1 + 1e-12):
                pi, f, g, it = _semidual_newton(C, a, b, e_cur, g0=g, tol=1e-10); its += it
                e_cur = max(eps, e_cur / 2.0)
            pi, f, g, it = _semidual_newton(C, a, b, eps, g0=g, tol=tol)
            return pi, f, g, its + it
        return _semidual_newton(C, a, b, eps, g0=g0, tol=tol)
    la, lb = np.log(a), np.log(b)
    f = np.zeros(C.shape[0]) if f0 is None else np.array(f0, float)
    g = np.zeros(C.shape[1]) if g0 is None else np.array(g0, float)
    for it in range(max_iter):
        g = -eps * _logsumexp((f[:, None] - C) / eps + la[:, None], axis=0)
        f_new = -eps * _logsumexp((g[None, :] - C) / eps + lb[None, :], axis=1)
        if np.max(np.abs(f_new - f)) < tol * max(1.0, np.max(np.abs(f_new))):
            f = f_new
            break
        f = f_new
    logpi = (f[:, None] + g[None, :] - C) / eps + la[:, None] + lb[None, :]
    return np.exp(logpi), f, g, it + 1


def ot_eps(C, a, b, eps, **kw):
    """OT_ε(μ, ν) = ⟨C, π⟩ + ε KL(π | a⊗b) at the Sinkhorn plan, and the plan."""
    pi, f, g, it = sinkhorn_plan(C, a, b, eps, **kw)
    with np.errstate(divide='ignore', invalid='ignore'):
        kl = np.sum(np.where(pi > 0, pi * (np.log(pi) - np.log(np.outer(a, b))), 0.0))
    return float(np.sum(pi * C) + eps * kl), pi


class EntropicPriorLoop(MixturePriorLoop):
    """(𝔄_m^ε): the averaged slow loop with the entropic OT leak toward ν = Σ_jq_jδ_{prior_j} at regularisation eps."""

    def __init__(self, net, f_vec, x_ref, period, G, sigma0, prior, eps, q=None, n_quad=256):
        super().__init__(net, f_vec, x_ref, period, G, sigma0, prior, n_quad=n_quad)
        self.K0 = self.prior.shape[0]
        self.q = np.ones(self.K0) / self.K0 if q is None else np.asarray(q, float) / np.sum(q)
        self.eps = float(eps)
        self._dual = None

    def plan(self, Th, warm=True):
        m = Th.shape[0]
        C = _cost_matrix(np.asarray(Th, float), self.prior, self.G)
        f0 = g0 = None
        if warm and self._dual is not None and self._dual[0].shape[0] == m:
            f0, g0 = self._dual
        pi, f, g, it = sinkhorn_plan(C, np.full(m, 1.0 / m), self.q, self.eps, f0=f0, g0=g0)
        self._dual = (f, g)
        return pi, C

    def leak_target(self, Th):
        pi, C = self.plan(Th)
        return Th.shape[0] * (pi @ self.prior)

    def J(self, Th):
        """J^ε_m = mL̄ + (σ₀m/2)OT_ε(μ_Θ, ν)."""
        m = Th.shape[0]
        C = _cost_matrix(np.asarray(Th, float), self.prior, self.G)
        val, _ = ot_eps(C, np.full(m, 1.0 / m), self.q, self.eps)
        return m * self.L_bar(Th) + 0.5 * self.sigma0 * m * val

    def drift(self, Th, T=None):
        T = self.leak_target(Th) if T is None else T
        a, c, s, diff, phi = self.features(Th)
        m = Th.shape[0]
        r = self.F - phi @ a / m
        net = self.net
        da = 1.0 - (a / net.a_max) ** 2
        dc = 1.0 - (c / net.c_max) ** 2
        sk = _sig(Th[:, 2]); ds = (net.s_max - net.s_min) * sk * (1 - sk)
        g = np.empty_like(Th)
        g[:, 0] = (phi.T @ r) * da
        g[:, 1] = a * dc / s ** 2 * ((phi * diff).T @ r)
        g[:, 2] = a * ds / s ** 3 * ((phi * diff ** 2).T @ r)
        return self.G[None, :] * g / self.N - self.sigma0 * (Th - T)

    def grad_J(self, Th):
        return -self.drift(Th) / self.G[None, :]


# ---------------------------------------------------------------------- collapsed configurations and the critical excess
def weighted_targets(prior, w, atoms, q, eps, G):
    """For the K-atom measure Σ_kw_kδ_{atoms_k}: the K × K₀ entropic plan to ν and the class targets T_k = Σ_j(π_kj/w_k)ϑ⁰_j, and the
    whitened conditional covariances Cov^G_k = G^{−1/2}Cov_{π(·|k)}(ϑ⁰)G^{−1/2}."""
    C = _cost_matrix(atoms, prior, G)
    pi, f, g, it = sinkhorn_plan(C, np.asarray(w, float), np.asarray(q, float), eps)
    cond = pi / pi.sum(axis=1, keepdims=True)
    T = cond @ prior
    sq = 1.0 / np.sqrt(np.asarray(G, float))
    cov = []
    for k in range(atoms.shape[0]):
        d = (prior - T[k][None, :]) * sq[None, :]
        cov.append(np.einsum('j,ja,jb->ab', cond[k], d, d))
    return pi, T, np.array(cov)


def weighted_energy_grad_entropic(loop, w, Th, eps, q):
    """∇_{ϑ_k}𝒢^ε_w / w_k for the K-atom measure with weights w: the fit's first variation plus σ₀G⁻¹(ϑ_k − T_k), T_k the class target of the
    K × K₀ entropic plan (the collapsed m-particle configuration's targets, by the class symmetry of the Sinkhorn plan)."""
    a, c, s, diff, phi = loop.features(Th)
    r = loop.F - phi @ (np.asarray(w, float) * a)
    _, T, _ = weighted_targets(loop.prior, w, Th, q, eps, loop.G)
    return loop.grad_V_L(r, Th) + loop.sigma0 * (Th - T) / loop.G[None, :]


def weighted_stationary_atoms_entropic(loop, w, Th0, eps, q, iters=100, h=1e-6, tol=1e-12, warn_tol=1e-9):
    """Newton with backtracking on the K-atom stationarity condition of the entropic leak; returns (atoms, residual). The caller must check
    the residual: a start far from the branch can stall (warned above warn_tol); continuation in ε or in the weights from a converged point is
    the reliable route (referee round 15, D3h-1)."""
    K = Th0.shape[0]; d = 3 * K
    v = Th0.ravel().copy()
    F = lambda v: weighted_energy_grad_entropic(loop, w, v.reshape(K, 3), eps, q).ravel()
    for _ in range(iters):
        gv = F(v); gn = np.abs(gv).max()
        if gn < tol:
            break
        J = np.zeros((d, d))
        for j in range(d):
            e = np.zeros(d); e[j] = h
            J[:, j] = (F(v + e) - F(v - e)) / (2 * h)
        step = np.linalg.solve(J, gv)
        lam = 1.0
        while lam > 1e-8:
            vn = v - lam * step
            if np.abs(F(vn)).max() < gn:
                break
            lam *= 0.5
        v = vn
    res = float(np.abs(F(v)).max())
    if res > warn_tol:
        warnings.warn(f"weighted_stationary_atoms_entropic: Newton did not converge (residual {res:.1e}, ε = {eps})", RuntimeWarning)
    return v.reshape(K, 3), res


def separation_blocks_entropic(loop, w, atoms, eps, q):
    """The separation blocks of the collapsed configuration under the entropic leak: σ₀(I − (2/ε)Cov^G_k) − R^G_k for each class k (the
    potentials g are stationary to first order under a zero-sum perturbation of a class); returns (blocks (K,3,3), their smallest
    eigenvalues, Cov^G)."""
    a, c, s, diff, phi = loop.features(atoms)
    r = loop.F - phi @ (np.asarray(w, float) * a)
    R = loop.residual_curvature(r, atoms)
    Rw = loop.sqG[None, :, None] * R * loop.sqG[None, None, :]
    _, T, cov = weighted_targets(loop.prior, w, atoms, q, eps, loop.G)
    blocks = loop.sigma0 * (np.eye(3)[None] - (2.0 / eps) * cov) - Rw
    return blocks, np.array([np.linalg.eigvalsh(B).min() for B in blocks]), cov
