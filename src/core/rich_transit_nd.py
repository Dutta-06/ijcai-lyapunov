"""
The rich transit in several dimensions (Proposition D‴.54): Gaussian kernels exp(−|x − c|²/(2σ²)) on samples in ℝⁿ, their centres moving
by the paradigm's position law, the masses small
====================================================================================================================================

The one-dimensional objects of `rich_transit` with the centres c_i ∈ ℝⁿ. A dictionary supplies the synthesis Φ(c)w and the correlations
⟨rφ_i⟩, ⟨r∇_cφ_i⟩ (and, for the Jacobian, ⟨r∇²_cφ_i⟩ and the Gram matrices ⟨φ_iφ_j⟩, ⟨φ_i∇_cφ_j⟩, ⟨∇_cφ_i∇_cφ_jᵀ⟩): `BumpsND` on any
samples, `BumpsGrid2` on a tensor grid in the plane (separable: every product through the one-dimensional factors). On it: the free law
on its curves with its exact Jacobian (`critical_rhs_nd`, `critical_jac_nd`: β = b0 + 2c sinh²(α/a_max), α̇ = curve mobility·⟨rφ_i⟩,
ċ_i = Gmw_i⟨r∇_cφ_i⟩ — b0 = 0 the active block's fast flow, b0 = εb the whole loop from α = 0, `loop_nd`), the block's least-squares
polish (`polish_block_nd`), the sliding inverse scale space (`sliding_iss_nd`), and the ε-free inner problem of k kernels entering
one feature together (`inner_k`): the k-fold near tie, whose shares for coinciding kernels are the profile split w_i = 𝒲(ξ + x_i),
Σ_iw_i = A — on the exponential branch the softmax e^{−κx_i}/Σ_je^{−κx_j} (`profile_split`).
"""
import time

import numpy as np
from scipy.integrate import solve_ivp

from src.core.rich_transit import (critical_w, curve_mobility, coef_crit, eta_crit, inner_coef, _alpha_of, transit_scale)
from src.core import rich_transit as _RT


class BumpsND:
    """Gaussian kernels of one width on the samples xs (N × n); ⟨f g⟩ the mean over the samples. Centres c: (d, n)."""

    def __init__(self, xs, sig):
        self.xs, self.sig = np.atleast_2d(np.asarray(xs, float)), float(sig)
        self.N, self.n = self.xs.shape

    def phi(self, c):
        """φ_j at the samples: (N, d)."""
        c = np.atleast_2d(np.asarray(c, float))
        d2 = ((self.xs[:, None, :] - c[None, :, :]) ** 2).sum(-1)
        return np.exp(-d2 / (2 * self.sig ** 2))

    def dphi(self, c):
        """∇_cφ_j at the samples: (N, d, n)."""
        c = np.atleast_2d(np.asarray(c, float))
        dx = self.xs[:, None, :] - c[None, :, :]
        return self.phi(c)[:, :, None] * dx / self.sig ** 2

    def synth(self, c, w):
        return self.phi(c) @ w

    def corr(self, c, r, second=False):
        """⟨rφ_j⟩ (d), ⟨r∇_cφ_j⟩ (d, n) and, with second, ⟨r∇²_cφ_j⟩ (d, n, n)."""
        c = np.atleast_2d(np.asarray(c, float))
        P = self.phi(c)
        dx = (self.xs[:, None, :] - c[None, :, :]) / self.sig ** 2
        rP = r[:, None] * P
        g = rP.sum(0) / self.N
        h = np.einsum('sj,sjk->jk', rP, dx) / self.N
        if not second:
            return g, h
        Q = np.einsum('sj,sjk,sjl->jkl', rP, dx, dx) / self.N - g[:, None, None] * np.eye(self.n) / self.sig ** 2
        return g, h, Q

    def grams(self, c):
        """P = ⟨φ_iφ_j⟩ (d, d), D[i, j] = ⟨φ_i∇_cφ_j⟩ (d, d, n), H[i, j] = ⟨∇_cφ_i∇_cφ_jᵀ⟩ (d, d, n, n)."""
        P = self.phi(c)
        dP = self.dphi(c)
        return (P.T @ P / self.N, np.einsum('si,sjk->ijk', P, dP) / self.N, np.einsum('sik,sjl->ijkl', dP, dP) / self.N)


class BumpsGrid2(BumpsND):
    """Gaussian kernels on the tensor grid x1 × x2 of samples in the plane (sample (a, b) at index a·len(x2) + b): φ_j = g_j ⊗ k_j, every
    synthesis, correlation and Gram matrix through the one-dimensional factors g_j(x1), k_j(x2)."""

    def __init__(self, x1, x2, sig):
        self.x1, self.x2 = np.asarray(x1, float), np.asarray(x2, float)
        super().__init__(np.stack(np.meshgrid(self.x1, self.x2, indexing='ij'), -1).reshape(-1, 2), sig)

    def _fac(self, c):
        c = np.atleast_2d(np.asarray(c, float))
        s2 = self.sig ** 2
        d1, d2 = self.x1[:, None] - c[None, :, 0], self.x2[:, None] - c[None, :, 1]
        g, k = np.exp(-d1 ** 2 / (2 * s2)), np.exp(-d2 ** 2 / (2 * s2))
        return g, k, g * d1 / s2, k * d2 / s2, d1, d2

    def phi(self, c):
        g, k = self._fac(c)[:2]
        return (g[:, None, :] * k[None, :, :]).reshape(self.N, -1)

    def dphi(self, c):
        g, k, gp, kp = self._fac(c)[:4]
        return np.stack([(gp[:, None, :] * k[None, :, :]).reshape(self.N, -1), (g[:, None, :] * kp[None, :, :]).reshape(self.N, -1)], -1)

    def synth(self, c, w):
        g, k = self._fac(c)[:2]
        return ((g * w) @ k.T).ravel()

    def corr(self, c, r, second=False):
        g, k, gp, kp, d1, d2 = self._fac(c)
        R = r.reshape(self.x1.size, self.x2.size)
        T, Tp = g.T @ R, gp.T @ R
        gg = (T * k.T).sum(1) / self.N
        h = np.stack([(Tp * k.T).sum(1), (T * kp.T).sum(1)], 1) / self.N
        if not second:
            return gg, h
        s2 = self.sig ** 2
        gpp, kpp = g * (d1 ** 2 / s2 - 1) / s2, k * (d2 ** 2 / s2 - 1) / s2
        q12 = (Tp * kp.T).sum(1)
        Q = np.stack([np.stack([((gpp.T @ R) * k.T).sum(1), q12], 1), np.stack([q12, (T * kpp.T).sum(1)], 1)], 1) / self.N
        return gg, h, Q

    def grams(self, c):
        g, k, gp, kp = self._fac(c)[:4]
        Gg, Gk = g.T @ g, k.T @ k
        Ggp, Gkp = g.T @ gp, k.T @ kp                                  # [i, j] = Σ g_i g′_j
        P = Gg * Gk / self.N
        D = np.stack([Ggp * Gk, Gg * Gkp], -1) / self.N
        H = np.empty(P.shape + (2, 2))
        H[..., 0, 0] = (gp.T @ gp) * Gk
        H[..., 0, 1] = Ggp.T * Gkp
        H[..., 1, 0] = Ggp * Gkp.T
        H[..., 1, 1] = Gg * (kp.T @ kp)
        return P, D, H / self.N


def _curve_terms(L, al, b0):
    """On the curve from α = 0 with the mass b0: w, dw/dα, the mobility M = Gβa′ and dM/dα (each without cancellation)."""
    x = np.clip(np.asarray(al, float) / L.A, -350.0, 350.0)
    sh, ch, th = np.sinh(x), np.cosh(x), np.tanh(x)
    beta = b0 + 2 * L.c * sh ** 2
    w = beta * L.A * th / L.m
    dw = (4 * L.c * sh ** 2 + beta / ch ** 2) / L.m
    M = L.G * (b0 / ch ** 2 + 2 * L.c * th ** 2)
    dM = 2 * L.G * th / (L.A * ch ** 2) * (2 * L.c - b0)
    return w, dw, M, dM


def critical_rhs_nd(L, K, fstar, b0=0.0):
    """The free law with moving centres on its curves in ℝⁿ, z = (α (d), c (d·n)): α̇ = G(b0 sech² + 2c tanh²)⟨rφ⟩, ċ = Gmw⟨r∇_cφ⟩,
    β = b0 + 2c sinh²(α/a_max) (each kernel kept on its curve; `rich_transit.critical_rhs` in ℝⁿ)."""
    n = K.n

    def f(t, z):
        d = z.size // (1 + n)
        al, c = z[:d], z[d:].reshape(d, n)
        w = critical_w(L, al, b0)
        g, h = K.corr(c, fstar - K.synth(c, w))
        return np.r_[curve_mobility(L, al, b0) * g, (L.G * L.m * w[:, None] * h).ravel()]
    return f


def critical_jac_nd(L, K, fstar, b0=0.0):
    """The exact Jacobian of `critical_rhs_nd`: with r = f* − Σ_jw_jφ_j, g_i = ⟨rφ_i⟩, h_i = ⟨r∇_cφ_i⟩,
    ∂α̇_i/∂α_j = δ_ijM′_ig_i − M_iw′_jP_ij, ∂α̇_i/∂c_j = M_i(δ_ijh_i − w_jD_ij), ∂ċ_i/∂α_j = Gm(δ_ijw′_ih_i − w_iw′_jD_ji),
    ∂ċ_i/∂c_j = Gmw_i(δ_ij⟨r∇²_cφ_i⟩ − w_jH_ij)."""
    n = K.n

    def jac(t, z):
        d = z.size // (1 + n)
        al, c = z[:d], z[d:].reshape(d, n)
        w, dw, M, dM = _curve_terms(L, al, b0)
        g, h, Q = K.corr(c, fstar - K.synth(c, w), second=True)
        P, D, H = K.grams(c)
        Gm = L.G * L.m
        J = np.zeros((d * (1 + n), d * (1 + n)))
        J[:d, :d] = -M[:, None] * P * dw[None, :] + np.diag(dM * g)
        Jac = -(M[:, None, None] * w[None, :, None]) * D
        Jac[np.arange(d), np.arange(d)] += M[:, None] * h
        J[:d, d:] = Jac.reshape(d, d * n)
        Jca = -Gm * (w[:, None, None] * dw[None, :, None]) * D.transpose(1, 0, 2)    # [i, j, k] = −Gm w_i w′_j D[j, i, k]
        Jca[np.arange(d), np.arange(d)] += Gm * dw[:, None] * h
        J[d:, :d] = Jca.transpose(0, 2, 1).reshape(d * n, d)
        Jcc = -Gm * (w[:, None] * w[None, :])[:, :, None, None] * H                  # [i, j, k, l]
        Jcc[np.arange(d), np.arange(d)] += Gm * w[:, None, None] * Q
        J[d:, d:] = Jcc.transpose(0, 2, 1, 3).reshape(d * n, d * n)
        return J
    return jac


def loop_nd(L, K, fstar, c0, eps, taus, b=None, rtol=1e-9, atol=1e-12, method='LSODA'):
    """The free law with moving centres from α = 0, masses εb (b = 1 by default), in τ = t/Λ(ε) on the grid taus (`exp_mirror_transit.
    run_moving` in ℝⁿ, with the exact Jacobian). Returns (taus, W (T, d), C (T, d, n), Λ, α (T, d))."""
    c0 = np.atleast_2d(np.asarray(c0, float))
    d, n = c0.shape
    b0 = eps * (np.ones(d) if b is None else np.asarray(b, float))
    Lam = transit_scale(L, eps)
    taus = np.asarray(taus, float)
    sol = solve_ivp(critical_rhs_nd(L, K, fstar, b0), (0, Lam * taus[-1]), np.r_[np.zeros(d), c0.ravel()], method=method,
                    jac=critical_jac_nd(L, K, fstar, b0), rtol=rtol, atol=atol, t_eval=Lam * taus)
    al = sol.y[:d].T
    return taus, critical_w(L, al, b0), sol.y[d:].T.reshape(-1, d, n), Lam, al


def polish_block_nd(K, fstar, w, c, s, info=False):
    """The block's fit from (w, c (k, n)) by sign-constrained least squares (trust-region reflective): min ½⟨(f* − Φ(c)w)²⟩ with
    s_iw_i ≥ 0 (`rich_transit.polish_block` in ℝⁿ). Returns (w, c), with info also the solver's status and evaluations."""
    from scipy.optimize import least_squares
    c = np.atleast_2d(np.asarray(c, float))
    k, n, N = w.size, K.n, K.N
    lo = np.r_[np.where(s > 0, 0.0, -np.inf), np.full(k * n, -np.inf)]
    hi = np.r_[np.where(s > 0, np.inf, 0.0), np.full(k * n, np.inf)]
    z0 = np.r_[np.where(s > 0, np.maximum(w, 1e-300), np.minimum(w, -1e-300)), c.ravel()]

    def res(z):
        return (fstar - K.synth(z[k:].reshape(k, n), z[:k])) / np.sqrt(N)

    def jac(z):
        cc = z[k:].reshape(k, n)
        dP = K.dphi(cc) * z[:k][None, :, None]
        return -np.c_[K.phi(cc), dP.reshape(N, k * n)] / np.sqrt(N)
    sol = least_squares(res, z0, jac=jac, bounds=(lo, hi), method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
    if info:
        return sol.x[:k], sol.x[k:].reshape(k, n), dict(status=int(sol.status), nfev=int(sol.nfev))
    return sol.x[:k], sol.x[k:].reshape(k, n)


def sliding_iss_nd(L, K, fstar, c0, rho, seed=1e-7, t_fast=1e5, tol=1e-9, exit_tol=1e-6, max_events=60, tie_tol=1e-9, gtol=1e-9,
                   max_doublings=14, polish=True, max_entries=None, polish_step=1.0, polish_wfrac=1e-2, exit_kkt=False, max_rhs=None):
    """The sliding inverse scale space in ℝⁿ (`rich_transit.sliding_iss` with the centres c (d, n)): the inactive kernels stay where they
    are, their u_j growing at ⟨rφ(c_j)⟩; the entrants (every kernel whose |u_j| reaches ρ_j within tie_tol of the first: ties enter
    together) join the block with the coefficient `seed`; the block relaxes by the fast flow (polished by least squares after the first
    t_fast, then checked) until stationary; a kernel whose coefficient falls below exit_tol·max leaves if its correlation points into its
    box. Returns dict(events = [(τ, j, 'enter'/'exit', c_j)], hist = [(τ, c, w, rms)], c, w, act, relaxed, polish_rejected,
    stiff_redone — the chunks of the flow redone by Radau, `rich_transit.relax_chunk`). The polish
    is accepted only if every kernel whose coefficient is below polish_wfrac·max|w| lands within polish_step kernel widths of its start
    or of the end of its correlation's ascent (`rich_transit.ascent_target`), else tried again after each doubling (`sliding_iss`).
    max_entries: stop as soon as that many kernels have entered, before the block relaxes. exit_kkt and the counts 'solver_redone',
    'solver_failed': as `rich_transit.sliding_iss` (referee round 48; a block whose chunk failed is left unrelaxed, without exits: round
    54). max_rhs: a budget of right-hand-side evaluations for the whole
    path, as `rich_transit.sliding_iss` ('out_of_budget', 'nfev'): a stiff block — fixed widths on a narrow feature, which LSODA hands to
    Radau chunk after chunk — stops the path before the block that spent it is relaxed, instead of running for hours; that block's state
    when the budget ran out is returned as 'stalled' (its kernels' indices, centres and coefficients; referee round 57), and each
    block's polishes as 'polish_log' (the solver's status and evaluations, accepted or not, the largest coefficient over the block's)."""
    c0 = np.atleast_2d(np.asarray(c0, float))
    d, n = c0.shape
    rho = np.broadcast_to(np.asarray(rho, float), (d,))
    u = np.zeros(d); c = np.array(c0, float); al = np.zeros(d)
    act = np.zeros(d, bool)
    tau, ev, hist, relaxed, n_rej, n_stiff, n_redo, n_fail = 0.0, [], [], [], 0, 0, 0, 0
    polish_log, stalled = [], None
    scale = float(np.abs(K.corr(c, fstar)[0]).max())
    rhs, jac = critical_rhs_nd(L, K, fstar), critical_jac_nd(L, K, fstar)
    nfev, out_of_budget = [0], False
    if max_rhs is not None:
        rhs_ = rhs

        def rhs(s_, y_):
            nfev[0] += 1
            if nfev[0] > max_rhs:
                raise _RT._OutOfBudget
            return rhs_(s_, y_)
    for _ in range(max_events):
        w = np.where(act, critical_w(L, al), 0.0)
        r = fstar - K.synth(c, w)
        rms = float(np.sqrt(np.mean(r ** 2)))
        hist.append((tau, c.copy(), w.copy(), rms))
        if rms < tol:
            break
        v = K.corr(c, r)[0]
        dt = np.full(d, np.inf)
        ina = ~act & (v != 0)
        dt[ina] = (rho[ina] - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])
        dmin = dt.min()
        if not np.isfinite(dmin):
            break
        ent = np.flatnonzero(dt <= dmin * (1 + tie_tol) + tie_tol * np.finfo(float).tiny)
        u[~act] += dmin * v[~act]
        tau += dmin
        for j in ent:
            u[j] = np.sign(v[j]) * rho[j]
            act[j] = True
            al[j] = np.sign(v[j]) * (seed * L.m * L.A ** 2 / (2 * L.c)) ** (1 / 3)
            ev.append((tau, int(j), 'enter', c[j].copy()))
        if max_entries is not None and sum(e[2] == 'enter' for e in ev) >= max_entries:
            break
        idx = np.flatnonzero(act)
        k = idx.size
        z = np.r_[al[idx], c[idx].ravel()]
        elapsed, chunk, ok, polished, failed = 0.0, float(t_fast), False, not polish, False
        plog = []
        polish_log.append(plog)
        try:
            for _k in range(max_doublings + 1):
                z, how = _RT.relax_chunk(rhs, z, chunk, jac=jac)
                n_stiff += how == 1
                n_redo += how == 2
                elapsed += chunk
                if how == 3:
                    n_fail += 1
                    failed = True
                    break
                if not polished:
                    w0 = critical_w(L, z[:k])
                    cz = z[k:].reshape(k, n)
                    wp, cp, pinfo = polish_block_nd(K, fstar, w0, cz, np.sign(u[idx]), info=True)
                    seedy = np.flatnonzero(np.abs(w0) <= polish_wfrac * np.abs(w0).max())
                    r0 = fstar - K.synth(cz, w0)
                    good = all(np.linalg.norm(cp[i] - cz[i]) <= polish_step * K.sig or
                               np.linalg.norm(cp[i] - _RT.ascent_target(K, r0, cz[i], np.sign(u[idx[i]]))) <= polish_step * K.sig
                               for i in seedy)
                    plog.append(dict(pinfo, accepted=bool(good), wmax_ratio=float(np.abs(wp).max() / np.abs(w0).max())))
                    if good:
                        z = np.r_[np.sign(wp) * _alpha_of(L, np.abs(wp), coef_crit), cp.ravel()]
                        polished = True
                    else:
                        n_rej += 1
                wi = critical_w(L, z[:k])
                ci = z[k:].reshape(k, n)
                g, h = K.corr(ci, fstar - K.synth(ci, wi))
                big = np.abs(wi) > exit_tol * np.abs(wi).max()
                si = np.sign(u[idx])
                if (np.all(np.abs(g[big]) <= gtol * scale) and np.all(np.abs(wi[big, None] * h[big]) <= gtol * scale)
                        and np.all(si[~big] * g[~big] <= gtol * scale)):
                    ok = True
                    break
                chunk = elapsed
        except _RT._OutOfBudget:
            out_of_budget = True
            relaxed.append(False)
            stalled = dict(idx=idx.copy(), c=z[k:].reshape(k, n).copy(), w=critical_w(L, z[:k]))   # the block when the budget ran out
            break
        relaxed.append(ok)
        al[idx] = z[:k]; c[idx] = z[k:].reshape(k, n)
        if failed:                    # a chunk failed: the block is left unrelaxed and no kernel exits on it (referee round 54)
            continue
        own = None
        if exit_kkt:
            Pc = K.phi(ci)
            own = (Pc ** 2).sum(0) / Pc.shape[0]
        out = _RT.leaves(si, g, wi, big, own)
        for q, i in enumerate(idx):
            if out[q]:
                act[i] = False; al[i] = 0.0
                ev.append((tau, int(i), 'exit', c[i].copy()))
    return dict(events=ev, hist=hist, c=c, w=np.where(act, critical_w(L, al), 0.0), act=act, relaxed=relaxed, polish_rejected=n_rej,
                stiff_redone=n_stiff, solver_redone=n_redo, solver_failed=n_fail, out_of_budget=out_of_budget, nfev=nfev[0],
                polish_log=polish_log, stalled=stalled)


def path_duals_nd(K, fstar, R):
    """`rich_transit.path_duals` in ℝⁿ: the sliding path's duals (`sliding_iss_nd`) reconstructed at its event times — on each interval the
    block's residual r_k is fixed, the inactive kernels' u grows at ⟨r_kφ(c_j)⟩, an entrant is pinned at its wall (±1: equal masses).
    Returns dict(taus, res, centres, U, ever)."""
    hist = R['hist']
    taus = [h[0] for h in hist]
    res = [fstar - K.synth(h[1], h[2]) for h in hist]
    cs = [h[1] for h in hist]
    d = len(cs[0])
    ever = sorted({j for t, j, k, c in R['events'] if k == 'enter'})
    u = np.zeros(d); act = np.zeros(d, bool); U = [u.copy()]
    ev = R['events']; e = 0
    for k in range(len(taus) - 1):
        v = K.corr(cs[k], res[k])[0]
        u = u + np.where(act, 0.0, (taus[k + 1] - taus[k]) * v)
        while e < len(ev) and ev[e][0] <= taus[k + 1] + 1e-12:
            t, j, kind, c = ev[e]
            if kind == 'enter':
                act[j] = True; u[j] = np.sign(u[j])
            else:
                act[j] = False
            e += 1
        U.append(u.copy())
    return dict(taus=taus, res=res, centres=cs, U=np.array(U), ever=ever)


def creep_constants_nd(L, K, P, j, c0j, taus, b=1.0):
    """`rich_transit.creep_constants` in ℝⁿ: the leading-order displacement of a kernel that does not enter, a vector, c_j = c0j + εC_j(τ)
    with dC_j/dτ = κ(μ_j)s_j⟨r_k∇_cφ(c0j)⟩ (κ(μ) = (πa_max²/(4c))b^{3/2}cos(πμ/2)/sin³(πμ/2), μ_j = 1 − |u_j| the path's margin), and its
    dual's gain εD_j, dD_j/dτ = C_j·⟨r_k∇_cφ(c0j)⟩, on each of the path's intervals (P from `path_duals_nd`). Returns (C (len(taus), n),
    D (len(taus)))."""
    pref = np.pi * L.A ** 2 / (4 * L.c) * b ** 1.5
    c0j = np.atleast_2d(np.asarray(c0j, float))
    n = c0j.shape[1]
    T = list(P['taus'])
    taus = np.asarray(taus, float)
    C, D = np.full((taus.size, n), np.nan), np.full(taus.size, np.nan)
    z = np.zeros(n + 1)
    for k in range(len(T)):
        t0 = T[k]
        t1 = T[k + 1] if k + 1 < len(T) else max(taus.max(), t0) + 1e-12
        g, h = K.corr(c0j, P['res'][k])
        v, gp = float(g[0]), h[0]
        u0 = P['U'][k][j]

        def f(t, z, u0=u0, v=v, gp=gp, t0=t0):
            u = u0 + (t - t0) * v
            mu = 1 - abs(u)
            return np.r_[pref * np.cos(np.pi * mu / 2) / np.sin(np.pi * mu / 2) ** 3 * np.sign(u) * gp, z[:n] @ gp]
        sel = (taus >= t0) & (taus <= t1)
        if t1 > t0:
            sol = solve_ivp(f, (t0, t1), z, method='LSODA', rtol=1e-10, atol=1e-14, dense_output=True)
            if sel.any():
                y = sol.sol(taus[sel]).T
                C[sel], D[sel] = y[:, :n], y[:, n]
            z = sol.y[:, -1]
    return C, D


class _OverBudget(Exception):
    """inner_k's wall-clock budget ran out inside an integration."""


def inner_k(L, K, fstar, cs, x, rho=1.0, xi0=400.0, s_end=2000.0, tol=1e-10, max_doublings=12, share_tol=None, max_secs=None):
    """The ε-free inner problem of k kernels entering one feature together (D‴.51(ii) with k kernels; Lemma D‴.53a: a gradient flow of the
    fit): the kernels at the centres cs (k, n), the first (the entrant) from the inner coordinate xi0 behind its wall and the j-th x_j
    behind the entrant (x_0 = 0), in s = (τ − τ_e)/ε^{1/2}: the fast flow on the critical curve times Λ₀ (dα/ds = Λ₀·2Gc tanh²(α/a)⟨rφ⟩,
    dc/ds = Λ₀Gmw⟨r∇_cφ⟩), the residual the feature fstar (nothing else fitted). Integrated over xi0/(the entrant's correlation) + s_end,
    then doubling until every rate is below tol·max|f*| — and, with share_tol, until the shares moved by less than share_tol over the last
    doubling (the absolute rates pass early while a tail kernel's tiny coefficient still slides: referee round 42). max_secs: a
    wall-clock budget — checked inside the integration too, so one long chunk cannot overrun it; on running out the state is the last
    completed chunk's and stationary is False. Returns dict(w, c, shares |w_j|/Σ|w|, stationary, s, drift: the shares' largest change
    over the last doubling, drifts: those changes per doubling, stop: why it stopped — 'converged' (the rates and the drift within their
    tolerances), 'cap' (max_doublings spent) or 'clock' (max_secs) —, remaining: the shares' change still to come if the drifts fall
    geometrically, drift²/(drift_prev − drift) (inf if they do not fall); a drift within share_tol over one doubling does not by itself
    bound the remaining change: referee round 47)."""
    lam0 = transit_scale(L, 1.0)
    cs = np.atleast_2d(np.asarray(cs, float))
    k, n = cs.shape
    sgn = np.sign(float(np.mean(fstar)) or 1.0)
    al0 = sgn * _alpha_of(L, -lam0 * rho * (xi0 + np.asarray(x, float)), eta_crit)
    z = np.r_[al0, cs.ravel()]
    f0, j0 = critical_rhs_nd(L, K, fstar), critical_jac_nd(L, K, fstar)
    t_end = np.inf if max_secs is None else time.time() + max_secs

    def rhs(s, z):
        if time.time() > t_end:
            raise _OverBudget
        return lam0 * f0(s, z)

    jac = lambda s, z: lam0 * j0(s, z)
    v_e = abs(float(K.corr(cs[:1], fstar)[0][0]))
    s_, chunk, ok = 0.0, xi0 * rho / v_e + float(s_end), False
    scale = tol * float(np.abs(fstar).max())
    sh_prev, drift, drifts, stop = None, np.inf, [], 'cap'
    for _ in range(max_doublings + 1):
        try:
            sol = solve_ivp(rhs, (s_, s_ + chunk), z, method='LSODA', jac=jac, rtol=1e-10, atol=1e-13)
        except _OverBudget:
            stop = 'clock'
            break
        z, s_ = sol.y[:, -1], s_ + chunk
        w_ = np.abs(critical_w(L, z[:k]))
        sh = w_ / w_.sum()
        drift = float(np.abs(sh - sh_prev).max()) if sh_prev is not None else np.inf
        if sh_prev is not None:
            drifts.append(drift)
        sh_prev = sh
        if np.all(np.abs(lam0 * f0(s_, z)) <= scale) and (share_tol is None or drift <= share_tol):
            ok, stop = True, 'converged'
            break
        if time.time() > t_end:
            stop = 'clock'
            break
        chunk = s_
    w = critical_w(L, z[:k])
    rem = (drifts[-1] ** 2 / (drifts[-2] - drifts[-1]) if len(drifts) > 1 and drifts[-2] > drifts[-1] else
           (0.0 if drifts and drifts[-1] == 0 else np.inf))
    return dict(w=w, c=z[k:].reshape(k, n), shares=np.abs(w) / np.abs(w).sum(), stationary=ok, s=s_, drift=drift, drifts=drifts,
                stop=stop, remaining=float(rem))


def profile_split(L, A, x, rho=1.0):
    """The k-fold split of coinciding kernels (D‴.54(iii)): they receive the same correlation, so their inner coordinates keep their gaps
    x_j (x_0 = 0) and the shares are w_j = 𝒲(ξ + x_j) with Σ_jw_j = |A| (ξ by bisection) — on the exponential branch the softmax
    e^{−κx_j}/Σ_ie^{−κx_i}, κ = πρ(2c)^{1/2}. Returns the shares w_j/|A|."""
    x = np.asarray(x, float)
    lo, hi = -1e3, 1e3
    for _ in range(200):
        m_ = 0.5 * (lo + hi)
        if inner_coef(L, m_ + x, rho).sum() > abs(A):
            lo = m_
        else:
            hi = m_
    return inner_coef(L, 0.5 * (lo + hi) + x, rho) / abs(A)
