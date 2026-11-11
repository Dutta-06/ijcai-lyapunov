"""
A faster evaluation of the same sliding-inverse-scale-space fit: the same algorithm, the same equations, a cheaper arithmetic
====================================================================================================================================

`imaging.sliding_path` (and its guarded copy `imaging_guard.sliding_path_guarded`) spends nearly all its time inside `scipy.optimize.least_squares`
(method 'trf', the exact trust-region solver, x_scale 'jac'): per iteration one SVD of the (N + n) x n augmented Jacobian (N pixels of the fitted window, n = 3k
or 4k variables) and one Jacobian assembled from broadcast products. This module computes the SAME path with four changes, none of which alters an equation:

  1. `trf_gram`: the trust-region reflective iteration of scipy (`scipy.optimize._lsq.trf.trf_bounds`, the same Coleman-Li scaling, trust-radius update,
     reflective step selection, termination tests; its helper functions are scipy's own) with the SVD of the augmented Jacobian replaced by the
     symmetric eigendecomposition of its n x n Gram matrix, A_h = D (J^T J) D + diag(diag_h): the trust-region sub-problem of Moré (`solve_lsq_trust_region`)
     only reads s, V and s*U^T f = V^T g_h, all functions of A_h and of the gradient g_h = D J^T f, which is computed from J directly. Cost per iteration: one
     n x N matrix product (N n^2 flops) instead of a Householder QR + bidiagonal SVD (about 4 N n^2). The Gram matrix squares the conditioning of the Jacobian
     (the columns are scaled to unit norm first, `x_scale='jac'`), so the singular values below about 1e-8 of the largest are not resolved: the rank test of
     the Gauss-Newton step therefore uses `RANK_RATIO` instead of scipy's eps * m.
  2. The Jacobian, the model and the correlations are assembled by batched matrix products (`Kern`): a kernel is a sum of T separable terms, so its image is
     a (ny x T)(T x nx) product; the erf and exp of the pixel edges are evaluated once per edge (n + 1 values for n pixels), the pixel's two edges sharing them.
     These are the same numbers as `PixelPSF.comps` to the last digit, the sums over the terms are re-ordered.
  3. `sliding_path_fast`: the guarded path (`cycle_guard=False` is the unguarded path) with these two replacements and, optionally, an EXACT early exit:
     `stop_kappa`: when the consumer reads the path with `stop_snr(…, kappa)`, the first entrant whose matched-filter significance is below kappa ends
     the path before it is fitted. `stop_snr` returns the block just before that entrant, which does not depend on anything that follows it, so the chosen block
     is the one of the full path; only `taus`, `snr` and `hist` end at that entrant. (Ignored if kappa < snr_min: the path's own stop is earlier.)
  4. The remaining Python overhead (repeated array creation) is removed where it is the same operation.

What is NOT identical: the iterates of the trust-region iteration differ in the rounding of the step (a Gram matrix instead of an SVD), so a block fit
agrees with scipy's to the tolerance of the iteration (the fit is stopped by the same tests on the same quantities), not bit for bit; the floating-point
level is studied in `tests/test_imaging_fast.py` and `results/synthetic/speed2.md`. `trf_gram` is also exposed alone (`fit_block_fast`) so that an A/B test
of a single block fit is possible.
"""
import numpy as np
from scipy.optimize import nnls
from scipy.optimize._lsq.common import (CL_scaling_vector, check_termination, make_strictly_feasible, solve_lsq_trust_region, update_tr_radius)
from scipy.optimize._lsq.trf import select_step
from scipy.special import erf

from . import imaging as IM

_SQ2 = np.sqrt(2.0)
_SQ2PI = np.sqrt(2.0 * np.pi)
_SQPI = np.sqrt(np.pi)
RANK_RATIO = 1e-6        # the Gauss-Newton step is tried if the smallest singular value of the (scaled, augmented) Jacobian exceeds this fraction of the largest
_EPS = np.finfo(float).eps


# ----------------------------------------------------------------------------------------------------------------------------------
# the kernels: separable terms, evaluated by edges
# ----------------------------------------------------------------------------------------------------------------------------------
def _axis(p, c, sg, deriv_s, deriv=True):
    """Ψ(p; c), ∂Ψ/∂c (and ∂Ψ/∂σ) for the pixel centres p (consecutive integers) and the kernels c with widths sg: `PixelPSF.axis`, from the n + 1 pixel edges.
    deriv=False: the value only (a 1-tuple), without the exponentials."""
    edges = np.empty(len(p) + 1)
    edges[:-1] = p - 0.5
    edges[-1] = p[-1] + 0.5
    u = (edges[:, None] - c[None, :]) / (_SQ2 * sg[None, :])
    E = erf(u)
    if not deriv:
        return (0.5 * (E[1:] - E[:-1]),)
    e = np.exp(-u ** 2)
    out = (0.5 * (E[1:] - E[:-1]), -(e[1:] - e[:-1]) / (_SQ2PI * sg[None, :]))
    if deriv_s:
        ue = u * e
        out += (-(ue[1:] - ue[:-1]) / (_SQPI * sg[None, :]),)
    return out


class Kern:
    """The separable terms of the kernels of a `PixelPSF` window: φ_i = Σ_t a_t g_y,t,i ⊗ g_x,t,i with t over the parts (and their high-pass partners)."""

    def __init__(self, psf):
        self.psf = psf
        root = psf._parent or psf
        self.root = root
        self.hp = psf.hp
        self.py, self.px = (root.py, root.px) if psf.hp else (psf.py, psf.px)
        self.crop = psf._crop if (psf.hp and psf._parent is not None) else None
        self.ny, self.nx, self.N = psf.ny, psf.nx, psf.N
        self.coef = np.array([a for a, _ in psf.parts for _ in range(2 if psf.hp else 1)]) * (np.array([1.0, -1.0] * len(psf.parts)) if psf.hp else 1.0)

    def terms(self, c, deriv_s=False, deriv=True):
        """(GY, GX, DGY, DGX[, SGY, SGX]) each (T, n, k): the factors of every term (the partner terms smoothed by the high-pass and cropped to the window);
        deriv=False: (GY, GX) only."""
        c = np.atleast_2d(np.asarray(c, float))
        c = c.reshape(-1, c.shape[-1])
        psf, root = self.psf, self.root
        k = len(c)
        s0 = c[:, 2] if c.shape[1] > 2 else np.full(k, psf.sigma)
        nf = (3 if deriv_s else 2) if deriv else 1
        P, step = len(psf.parts), 2 if self.hp else 1
        T = P * step
        if self.crop is not None:
            ya, yb, xa, xb = self.crop
        else:
            ya, yb, xa, xb = 0, self.ny, 0, self.nx
        FY = np.empty((nf, T, yb - ya, k))
        FX = np.empty((nf, T, xb - xa, k))
        for m, (a, kap) in enumerate(psf.parts):
            ax = _axis(self.px, c[:, 0], kap * s0, deriv_s, deriv)
            ay = _axis(self.py, c[:, 1], kap * s0, deriv_s, deriv)
            if self.hp:
                # one product for the value and the derivatives: S (n, n) @ [f_0 | f_1 | f_2] (n, nf k)
                sx = (root.Sx @ np.concatenate(ax, axis=1))[xa:xb].reshape(xb - xa, nf, k)
                sy = (root.Sy @ np.concatenate(ay, axis=1))[ya:yb].reshape(yb - ya, nf, k)
                for j in range(nf):
                    f = kap if j == 2 else 1.0
                    FX[j, m * 2] = ax[j][xa:xb] * f
                    FY[j, m * 2] = ay[j][ya:yb] * f
                    FX[j, m * 2 + 1] = sx[:, j] * f
                    FY[j, m * 2 + 1] = sy[:, j] * f
            else:
                for j in range(nf):
                    f = kap if j == 2 else 1.0
                    FX[j, m] = ax[j] * f
                    FY[j, m] = ay[j] * f
        if not deriv:
            return [FY[0], FX[0]]
        out = [FY[0], FX[0], FY[1], FX[1]]
        if deriv_s:
            out += [FY[2], FX[2]]
        return out

    def _mm(self, FY, FX):
        """Σ_t a_t FY_t ⊗ FX_t for every kernel: (k, N)."""
        A = (self.coef[:, None, None] * FY).transpose(2, 1, 0)           # (k, ny, T)
        B = FX.transpose(2, 0, 1)                                         # (k, T, nx)
        return np.matmul(A, B).reshape(A.shape[0], self.N)

    def phiT(self, c):
        """φ_i on the window, (k, N)."""
        GY, GX = self.terms(c, deriv=False)[:2]
        return self._mm(GY, GX)

    def synth(self, c, w, f=None):
        """Σ_i w_i φ(c_i) on the window, (ny, nx); f: the kernel factors `terms(c, …)` if already evaluated."""
        if np.size(w) == 0:
            return np.zeros((self.ny, self.nx))
        GY, GX = (self.terms(c, deriv=False) if f is None else f)[:2]
        T, ny, k = GY.shape
        Y = (self.coef[:, None, None] * GY * np.asarray(w, float)[None, None, :]).transpose(1, 0, 2).reshape(ny, T * k)
        X = GX.transpose(1, 0, 2).reshape(self.nx, T * k)
        return Y @ X.T

    def jacT(self, c, w, sw, f=None):
        """The Jacobian of the residual sw∘(y − Σ w_iφ_i) in z = (w, x, y[, σ]), transposed: (3k or 4k, N); f: the kernel factors `terms(c, deriv_s=<c has a width>)`."""
        c = np.atleast_2d(np.asarray(c, float))
        wid = c.shape[1] > 2
        k = len(w)
        w = np.asarray(w, float)
        f = self.terms(c, deriv_s=wid) if f is None else f
        GY, GX, DGY, DGX = f[:4]
        JT = np.empty(((4 if wid else 3) * k, self.N))
        JT[:k] = self._mm(GY, GX)
        JT[k:2 * k] = self._mm(GY, DGX)
        JT[2 * k:3 * k] = self._mm(DGY, GX)
        if wid:
            JT[3 * k:] = self._mm(f[4], GX) + self._mm(GY, f[5])
        JT[k:] *= np.tile(w, 3 if wid else 2)[:, None]
        JT *= -np.asarray(sw, float).reshape(1, -1)
        return JT

    def corr(self, c, Rw):
        """v_i = Σ_p Rw_p φ_i(p) and ∂v/∂(x_i, y_i): `PixelPSF.corr`."""
        GY, GX, DGY, DGX = self.terms(c)[:4]
        v = vx = vy = 0.0
        for t in range(GY.shape[0]):
            a = self.coef[t]
            A, B = GY[t].T @ Rw, DGY[t].T @ Rw
            v, vx, vy = v + a * (A * GX[t].T).sum(1), vx + a * (A * DGX[t].T).sum(1), vy + a * (B * GX[t].T).sum(1)
        return v, vx, vy


class _Op:
    """J_h = J D as an operator (only `dot`, which `select_step` calls)."""
    __slots__ = ('JT', 'd')

    def __init__(self, JT, d):
        self.JT, self.d = JT, d

    def dot(self, v):
        return self.JT.T @ (self.d * v)


# ----------------------------------------------------------------------------------------------------------------------------------
# the trust-region reflective iteration on the Gram matrix
# ----------------------------------------------------------------------------------------------------------------------------------
def trf_gram(fun, jacT, x0, lb, ub, ftol, xtol, gtol, max_nfev):
    """`scipy.optimize.least_squares(fun, x0, jac, bounds=(lb, ub), method='trf', x_scale='jac', ftol, xtol, gtol, max_nfev)` (the exact trust-region solver),
    with the jacobian given transposed (n, m) and the SVD of the augmented Jacobian replaced by the eigendecomposition of its Gram matrix. Returns (x, nfev, status)."""
    x = make_strictly_feasible(np.array(x0, float), lb, ub)
    f = fun(x)
    nfev = 1
    JT = jacT(x)
    n, m = JT.shape
    cost = 0.5 * np.dot(f, f)
    g = JT @ f
    A = JT @ JT.T
    scale_inv = np.sqrt(np.diag(A))
    scale_inv[scale_inv == 0] = 1
    scale = 1 / scale_inv
    v, dv = CL_scaling_vector(x, g, lb, ub)
    v[dv != 0] *= scale_inv[dv != 0]
    Delta = np.linalg.norm(x * scale_inv / v ** 0.5)
    if Delta == 0:
        Delta = 1.0
    if max_nfev is None:
        max_nfev = x.size * 100
    alpha = 0.0
    status = None
    m_rank = RANK_RATIO / _EPS
    while True:
        v, dv = CL_scaling_vector(x, g, lb, ub)
        g_norm = np.linalg.norm(g * v, ord=np.inf)
        if g_norm < gtol:
            status = 1
        if status is not None or nfev == max_nfev:
            break
        v[dv != 0] *= scale_inv[dv != 0]
        d = v ** 0.5 * scale
        diag_h = g * dv * scale
        g_h = d * g
        Ah = A * d[:, None] * d[None, :]
        Ah[np.diag_indices(n)] += diag_h
        lam, V = np.linalg.eigh(Ah)
        lam, V = lam[::-1], V[:, ::-1]
        s = np.sqrt(np.maximum(lam, 0.0))
        suf = V.T @ g_h                                            # = s * (U^T f_aug)
        with np.errstate(divide='ignore', invalid='ignore'):
            uf = np.where(s > 0, suf / np.where(s > 0, s, 1.0), 0.0)
        J_h = _Op(JT, d)
        theta = max(0.995, 1 - g_norm)
        actual_reduction = -1
        while actual_reduction <= 0 and nfev < max_nfev:
            p_h, alpha, _ = solve_lsq_trust_region(n, m_rank, uf, s, V, Delta, initial_alpha=alpha)
            p = d * p_h
            step, step_h, predicted_reduction = select_step(x, J_h, diag_h, g_h, p, p_h, d, Delta, lb, ub, theta)
            x_new = make_strictly_feasible(x + step, lb, ub, rstep=0)
            f_new = fun(x_new)
            nfev += 1
            step_h_norm = np.linalg.norm(step_h)
            if not np.all(np.isfinite(f_new)):
                Delta = 0.25 * step_h_norm
                continue
            cost_new = 0.5 * np.dot(f_new, f_new)
            actual_reduction = cost - cost_new
            Delta_new, ratio = update_tr_radius(Delta, actual_reduction, predicted_reduction, step_h_norm, step_h_norm > 0.95 * Delta)
            step_norm = np.linalg.norm(step)
            status = check_termination(actual_reduction, cost, step_norm, np.linalg.norm(x), ratio, ftol, xtol)
            if status is not None:
                break
            alpha *= Delta / Delta_new
            Delta = Delta_new
        if actual_reduction > 0:
            x, f, cost = x_new, f_new, cost_new
            JT = jacT(x)
            g = JT @ f
            A = JT @ JT.T
            scale_inv = np.maximum(np.sqrt(np.diag(A)), scale_inv)
            scale = 1 / scale_inv
    return x, nfev, (0 if status is None else status)


def fit_block_fast(psf, y, sw, c, w, c_lo, c_hi, joint=True, max_nfev=200, tol=1e-12, kern=None):
    """`imaging.block_fit` with `trf_gram`. Returns (c, w, info)."""
    k, nc = len(w), c.shape[1]
    if k == 0:
        return c, w, dict(nfev=0, status=0)
    K = kern or Kern(psf)
    swv = np.asarray(sw, float) * np.ones((psf.ny, psf.nx))
    yv = (swv * y).ravel()
    sv = swv.reshape(1, -1)
    w1, _ = nnls((K.phiT(c) * sv).T, yv, maxiter=200 * k)
    if not joint:
        return c, w1, dict(nfev=0, status=0)
    lo = np.r_[np.zeros(k), c_lo.T.ravel()]
    hi = np.r_[np.full(k, np.inf), c_hi.T.ravel()]
    z0 = np.r_[w1, np.clip(c, c_lo, c_hi).T.ravel()]

    memo = {}

    def fun(z):
        cc = z[k:].reshape(nc, k).T
        f = K.terms(cc, deriv_s=nc > 2)
        memo['z'], memo['f'] = z.copy(), f
        return (swv * (y - K.synth(cc, z[:k], f))).ravel()

    def jacT(z):
        cc = z[k:].reshape(nc, k).T
        f = memo['f'] if ('z' in memo and np.array_equal(memo['z'], z)) else None
        return K.jacT(cc, z[:k], swv, f)
    try:
        x, nfev, status = trf_gram(fun, jacT, z0, lo, hi, tol, tol, tol, max_nfev)
    except np.linalg.LinAlgError:                       # the eigensolver did not converge: scipy's own solver (SVD) for this block
        return IM.block_fit(psf, y, sw, c, w, c_lo, c_hi, joint=joint, max_nfev=max_nfev, tol=tol)
    return x[k:].reshape(nc, k).T.copy(), x[:k], dict(nfev=int(nfev), status=int(status))


def refine_fast(psf, y, sw, c, w, sig_bounds=None, box=1.0, sw_fn=None, rounds=2, tol=1e-10, max_nfev=100, sw_on='model'):
    """`imaging.refine` with `fit_block_fast`."""
    swv = np.asarray(sw, float) * np.ones((psf.ny, psf.nx))
    c, w = np.array(c, float), np.array(w, float)
    if len(w) == 0:
        return c, w, swv
    K = Kern(psf)
    for _ in range(rounds):
        lo, hi = c - box, c + box
        if c.shape[1] == 3:
            lo[:, 2], hi[:, 2] = sig_bounds
        c, w, _ = fit_block_fast(psf, y, swv, c, w, lo, hi, max_nfev=max_nfev, tol=tol, kern=K)
        if sw_fn is not None:
            swv = np.asarray(sw_fn(K.synth(c, w)) if sw_on == 'model' else sw_fn(c, w), float) * np.ones((psf.ny, psf.nx))
    return c, w, swv


# ----------------------------------------------------------------------------------------------------------------------------------
# the lattice correlations
# ----------------------------------------------------------------------------------------------------------------------------------
class FastLattice(IM.Lattice):
    """`imaging.Lattice` with a cheaper `corr`. With a high-pass the kernel is φ − S φ_smoothed term by term, so Σ_t a_t g_y,tᵀ R g_x,t over the parts and their partners equals
    Σ_m a_m g_y,mᵀ (R − S_yᵀ R S_x) g_x,m: the high-pass is applied (transposed) to the residual once instead of to every lattice factor, and the parts are one stacked product.
    Same numbers as `Lattice.corr` up to the order of the sums."""

    def __init__(self, psf, h=0.5, margin=1.0):
        super().__init__(psf, h, margin)
        root = psf._parent or psf
        self._root = root
        P = len(psf.parts)
        step = 2 if psf.hp else 1
        self._gy = [self.terms[m * step][1] for m in range(P)]          # (ny, nly) of the part m
        self._gx = [self.terms[m * step][2] for m in range(P)]
        self._a = [self.terms[m * step][0] for m in range(P)]
        self._Y = np.ascontiguousarray(np.concatenate(self._gy, axis=1).T)  # (P nly, ny)
        self._nly = self._gy[0].shape[1]
        self._ok = (psf._parent is None)

    def corr(self, Rw):
        if not self._ok:
            return super().corr(Rw)
        root = self._root
        Rh = Rw - root.Sy.T @ Rw @ root.Sx if self.psf.hp else Rw
        B = self._Y @ Rh
        out = 0.0
        for m, (a, gx) in enumerate(zip(self._a, self._gx)):
            out = out + a * (B[m * self._nly:(m + 1) * self._nly] @ gx)
        return out


# ----------------------------------------------------------------------------------------------------------------------------------
# the path
# ----------------------------------------------------------------------------------------------------------------------------------
def sliding_path_fast(psf, y, sw, lat, rho=1.0, max_entries=40, box=2.0, tie_tol=1e-9, vtol=1e-12, chi_stop=0.0, joint=True,
                      sw_fn=None, max_nfev=200, rtrunc=5.5, local=True, jump=None, extra=3, sig_bounds=None, tol=1e-12, passes=2, r_ex=0.0, sw_on='model',
                      snr_min=None, snr_extra=3, order='dual', fit_radius=None, tau_max=None, cycle_guard=True, stop_kappa=None, incremental=True, refresh=64):
    """`imaging_guard.sliding_path_guarded` (all its options, the same result) with `fit_block_fast` for the block fits and, if `stop_kappa` is given and not
    below `snr_min`, the exact early exit described in the module docstring: the path ends at the first entrant of matched-filter significance below `stop_kappa`,
    before that entrant is fitted (the result's `taus` and `snr` then end with that entrant, `hist` with the block before it).
    incremental: the model frame is carried from step to step (the refitted units' old image subtracted, the new one added; the others' image is `model − old image` for the
    window of the next fit; a unit's exit subtracts its image; every `refresh` steps the model is summed afresh from all units) instead of being summed over all units at every step and at
    every fit (the cost per step then does not grow with the number of units): the same numbers to the rounding of a few additions (about 1e-13 of the model's peak); False: summed afresh each time."""
    ny, nx = psf.ny, psf.nx
    swv = np.asarray(sw, float) * np.ones((ny, nx))
    U = np.zeros(lat.shape)
    taken = np.zeros(lat.shape, bool)
    barred = np.zeros(lat.shape, bool)
    last_exit, bars = {}, []
    nc = 2 if sig_bounds is None else 3
    ids, c, w, c0 = [], np.zeros((0, nc)), np.zeros(0), np.zeros((0, nc))
    tau, n_ent, events, hist, taus, snrs = 0.0, 0, [], [], [], []
    kern_full = Kern(psf)
    early = stop_kappa is not None and snr_min is not None and stop_kappa >= snr_min

    def bounds(c0_):
        lo, hi = c0_ - box, c0_ + box
        if nc == 3:
            lo[:, 2], hi[:, 2] = sig_bounds
        return lo, hi

    def start(b, a):
        return np.r_[lat.pos(b, a), psf.sigma] if nc == 3 else lat.pos(b, a)
    reach = psf.reach(psf.sigma if sig_bounds is None else sig_bounds[1])
    dint = 2.0 * rtrunc * reach
    model = np.zeros((ny, nx))

    def chi2():
        return float(np.mean((swv * (y - model)) ** 2))
    hist.append(dict(tau=0.0, n_ent=0, c=c.copy(), w=w.copy(), chi2=chi2(), ids=[], nfev=0))
    vscale, stop_at, best, stop_snr_at = None, None, (0.0, 0), None
    done = False
    for _ in range(11 * max_entries):
        Rw = swv ** 2 * (y - model)
        V = lat.corr(Rw)
        if vscale is None:
            vscale = float(np.abs(V).max())
            if vscale <= 0.0:
                break
        blocked = taken | barred if cycle_guard else taken
        if r_ex > 0.0 and len(w):
            blocked = blocked | lat.near(c[:, :2], r_ex)
        v = np.where(blocked, 0.0, V)
        grow = v > vtol * vscale
        if not grow.any():
            break
        if order == 'greedy':
            z = np.where(grow, v / np.sqrt(np.maximum(lat.norm2(swv ** 2), 1e-300)), -np.inf)
            ent = np.argwhere(z == z.max())[:1]
            tau += 1.0
        else:
            dt = np.full(lat.shape, np.inf)
            dt[grow] = (rho - U[grow]) / v[grow]
            dmin = float(dt.min())
            if tau_max is not None and tau + dmin > tau_max:
                break
            ent = np.argwhere(dt <= dmin * (1.0 + tie_tol) + tie_tol * 1e-300)
            U[~blocked] += dmin * v[~blocked]
            tau += dmin
        seeds = []
        for b, a in ent:
            U[b, a] = rho
            taken[b, a] = True
            seeds.append(len(ids))
            ids.append((int(b), int(a)))
            c0 = np.vstack([c0, start(b, a)[None, :]])
            c = np.vstack([c, start(b, a)[None, :]])
            w = np.r_[w, 0.0]
            events.append((tau, (int(b), int(a)), 'enter'))
            ph_ = psf.phi(start(b, a)[None, :])[:, 0]
            snrs.append(float(v[b, a] / np.sqrt(np.sum(swv.ravel() ** 2 * ph_ ** 2))))
            if early and snrs[-1] < stop_kappa:
                done = True
            if snr_min is not None and snrs[-1] < snr_min and stop_snr_at is None:
                stop_snr_at = n_ent + 1 + snr_extra
            if taus and tau / taus[-1] > best[0]:
                best = (tau / taus[-1], n_ent)
                stop_at = None if (jump is None or best[0] < jump) else n_ent + 1 + extra
            taus.append(tau)
            n_ent += 1
        if done:
            break
        if not local:
            S = np.arange(len(w))
        elif fit_radius is not None:
            S = np.flatnonzero(np.linalg.norm(c[:, None, :2] - c[seeds][None, :, :2], axis=2).min(1) < fit_radius)
        else:
            S = IM._cluster(c, seeds, dint)
        nfev = 0
        inc = incremental and local and len(S) < len(w)
        pre = kern_full.synth(c[S], w[S]) if inc else None
        for _pass in range(1 if sw_fn is None else passes):
            if local:
                r = rtrunc * reach
                ya, yb = max(int(np.floor(c[S, 1].min() - r)), 0), min(int(np.ceil(c[S, 1].max() + r)) + 1, ny)
                xa, xb = max(int(np.floor(c[S, 0].min() - r)), 0), min(int(np.ceil(c[S, 0].max() + r)) + 1, nx)
                if inc:
                    m_oth = (model - pre)[ya:yb, xa:xb]
                else:
                    oth = np.setdiff1d(np.arange(len(w)), S)
                    m_oth = kern_full.synth(c[oth], w[oth])[ya:yb, xa:xb]
                pw = psf.window(ya, yb, xa, xb)
                lo, hi = bounds(c0[S])
                cS, wS, info = fit_block_fast(pw, y[ya:yb, xa:xb] - m_oth, swv[ya:yb, xa:xb], c[S], w[S], lo, hi,
                                              joint=joint, max_nfev=max_nfev, tol=tol)
                c[S], w[S] = cS, wS
                if inc:
                    post = kern_full.synth(cS, wS)
                    model = model - pre + post
                    pre = post
                else:
                    model = kern_full.synth(c, w)
            else:
                lo, hi = bounds(c0)
                c, w, info = fit_block_fast(psf, y, swv, c, w, lo, hi, joint=joint, max_nfev=max_nfev, tol=tol, kern=kern_full)
                model = kern_full.synth(c, w)
            nfev += info['nfev']
            if sw_fn is not None:
                swv = np.asarray(sw_fn(model) if sw_on == 'model' else sw_fn(c, w), float) * np.ones((ny, nx))
        Rw = swv ** 2 * (y - model)
        cand = w <= 1e-12 * max(w.max(), 1e-300)
        gone = np.zeros(len(w), bool)
        if cand.any():
            ic = np.flatnonzero(cand)
            gone[ic] = kern_full.corr(c[ic], Rw)[0] <= vtol * vscale
        gone_img = kern_full.synth(c[gone], w[gone]) if (incremental and gone.any()) else None
        for i in np.flatnonzero(gone)[::-1]:
            events.append((tau, ids[i], 'exit'))
            taken[ids[i]] = False
            U[ids[i]] = rho
            if cycle_guard:
                if abs(tau - last_exit.get(ids[i], -np.inf)) <= 1e-12 * max(1.0, abs(tau)):
                    barred[ids[i]] = True
                    bars.append((tau, ids[i]))
                last_exit[ids[i]] = tau
            del ids[i]
            c, w, c0 = np.delete(c, i, 0), np.delete(w, i), np.delete(c0, i, 0)
        if gone.any():
            if incremental:
                model = model - gone_img
            else:
                model = kern_full.synth(c, w)
        if incremental and refresh and len(hist) % refresh == 0:
            model = kern_full.synth(c, w)
        hist.append(dict(tau=tau, n_ent=n_ent, c=c.copy(), w=w.copy(), chi2=chi2(), ids=list(ids), nfev=nfev))
        if (hist[-1]['chi2'] <= chi_stop or n_ent >= max_entries or (stop_at is not None and n_ent >= stop_at)
                or (stop_snr_at is not None and n_ent >= stop_snr_at)):
            break
    return dict(events=events, hist=hist, taus=np.array(taus), snr=np.array(snrs), barred=bars)
