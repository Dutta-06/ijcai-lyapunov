"""
The numerics of the real-window mechanism test (`exp_mechreal.py`): the network with a choice of readout, its gradient flow, and two discrete optimizers
=======================================================================================================================================================

Self-contained on purpose (numpy, scipy and `src/core/imaging.PixelPSF` / `Lattice` only), so that the same file runs on a remote CPU worker; `tests/test_mechreal.py` pins it to
`experiments/exp_network_imaging.py` (right-hand side, time unit and a whole trajectory of the paper's unit on a simulated frame).

The model (as `exp_network_imaging`): units i = 1..d on the start lattice (spacing 1 px, margin 1 px) of an ny x nx frame, kernel phi(c_i) = the pixel-integrated PSF (a sum of
Gaussians, high-passed with the border-consistent filter H if hp is set), loss F = 1/2 sum_p W_p (y - sum_i w_i phi_i)_p^2 with fixed weights W. The readouts w_i = W(p_i):
    P0  the paper's unit      p = (alpha, nu), w = nu^2 a_max tanh(alpha / a_max) / m, rates (G m, gamma m / 4), start (0, sqrt(eps));     clock Lambda(eps) = (pi a_max / (4G)) (2 / (v_s eps))^(1/2)
    D2  depth 2               p = (u, v),      w = u^2 - v^2,                           rates (1, 1),            start (sqrt(eps), sqrt(eps)); clock ln(1 / eps) / 4
    D1  depth 1               p = (w),         w = p,                                   rate  G m,               start 0;                      clock Lambda(eps) (the paper unit's physical time)
(the constants of `exp_parametrization`: m = 12, G = 1, gamma = 0.1, a_max = 1, the centres' rate G m in every readout, `factorized_rbf.Product2(1.0)`, `Linear(LR_C)`).
The flow: p_l' = lr_l dW/dp_l g_i, c_i' = lr_c w_i grad_c g_i, g_i = sum_p W_p r_p phi_i(p); the time of a snapshot is t = clock(eps) tau.

Solver: Radau (BDF from the last accepted state if Radau fails), rtol 1e-8, atol 1e-10 sqrt(eps) on the readout parameters (1e-10 for D1, whose amplitudes start at 0) and 1e-12 on the centres, with an ANALYTIC Jacobian
(`jac`): for the active units (those whose amplitude matters) the exact columns from the Gram matrices of the kernels and their centre derivatives, for the others the own block
with the drive frozen (the structure of the harness's approximate Jacobian, there by finite differences); all units active gives the exact dense Jacobian (used for D1, where every
unit moves). A wall-clock budget ends a run with status 'timeout' (the last accepted state is kept).

Discrete optimizers (the conventions of `factorized_rbf.train_adam` and `theory3f_escape.Opt`): 'adam' steps on the raw negative gradient of F (unit rates, one step size for every
parameter and centre, per-coordinate second moment, bias-corrected, eps_adam 1e-8); 'adam_g' the same with ONE scalar second moment, the mean of the squared gradient after the
gradient is multiplied by the flow's group rates (it keeps the relative size of the coordinates' steps).
"""
import os
import sys
import time

import numpy as np
import scipy.sparse as sp
from scipy.integrate import BDF, Radau
from scipy.special import erf

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)
try:
    from src.core import imaging as IM            # noqa: E402
except ImportError:                               # the remote worker ships imaging.py next to this file
    import imaging as IM                          # noqa: E402

LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)   # exp_parametrization.LAW
RTOL, ATOL_P, ATOL_C = 1e-8, 1e-10, 1e-12
ACT_FACTOR = 1e3                                  # P0: a unit is active when nu^2 > ACT_FACTOR eps (the harness's rule)
ACT_W = 1e-4                                      # D2: a unit is active when |w| > ACT_W (data in units of the flux limit)
FD_C = 1e-5                                       # px: the central difference of the own second derivative of the drive in the centre
LAT_H, LAT_MARGIN = 1.0, 1.0
STATUS = {'ok': 0.0, 'timeout': 1.0, 'failed': 2.0, 'error': 3.0, 'diverged': 4.0}


def lam_of(eps, law=LAW):
    """Lambda(eps) = (pi a_max / (4G)) (2 / (v_s eps))^(1/2), v_s = gamma a_max^2 / (4G): the paper unit's time unit (`rich_transit.transit_scale`)."""
    vs = law['gamma'] * law['a_max'] ** 2 / (4.0 * law['G'])
    return float(np.pi * law['a_max'] / (4.0 * law['G']) * np.sqrt(2.0 / (vs * eps)))


def sech2(x):
    x = np.minimum(np.abs(np.asarray(x, float)), 300.0)
    return 1.0 / np.cosh(x) ** 2


class Paper:
    """P0: w = nu^2 a_max tanh(alpha / a_max) / m, p = (alpha, nu)."""
    name, q = 'P0', 2

    def __init__(self, law=LAW):
        self.m, self.G, self.gam, self.A = law['m'], law['G'], law['gamma'], law['a_max']
        self.lr = (self.G * self.m, self.gam * self.m / 4.0)
        self.lr_c = self.G * self.m

    def init(self, d, eps):
        return [np.zeros(d), np.full(d, np.sqrt(eps))]

    def w(self, p):
        return p[1] ** 2 * self.A * np.tanh(p[0] / self.A) / self.m

    def dw(self, p):
        a, ap = self.A * np.tanh(p[0] / self.A), sech2(p[0] / self.A)
        return [p[1] ** 2 * ap / self.m, 2.0 * p[1] * a / self.m]

    def ddw(self, p):
        a, ap = self.A * np.tanh(p[0] / self.A), sech2(p[0] / self.A)
        app = -2.0 * a * ap / self.A ** 2
        off = 2.0 * p[1] * ap / self.m
        return [[p[1] ** 2 * app / self.m, off], [off, 2.0 * a / self.m]]

    def clock(self, eps):
        return lam_of(eps, dict(m=self.m, G=self.G, gamma=self.gam, a_max=self.A))

    def active(self, p, eps):
        return p[1] ** 2 > ACT_FACTOR * eps

    def atol(self, eps, atol_p):
        return atol_p * np.sqrt(eps)                                    # the harness's rule: (alpha, nu) are resolved down to the start nu = sqrt(eps)


class Depth2:
    """D2: w = u^2 - v^2, p = (u, v), started at (s, s), s^2 = eps (`factorized_rbf.Product2(1.0)`)."""
    name, q = 'D2', 2

    def __init__(self, law=LAW):
        self.lr = (1.0, 1.0)
        self.lr_c = law['G'] * law['m']

    def init(self, d, eps):
        return [np.full(d, np.sqrt(eps)), np.full(d, np.sqrt(eps))]

    def w(self, p):
        return p[0] ** 2 - p[1] ** 2

    def dw(self, p):
        return [2.0 * p[0], -2.0 * p[1]]

    def ddw(self, p):
        z = np.zeros_like(p[0])
        return [[2.0 + z, z], [z, -2.0 + z]]

    def clock(self, eps):
        return float(np.log(1.0 / eps) / (4.0 * self.lr[0]))

    def active(self, p, eps):
        return np.abs(self.w(p)) > ACT_W

    def atol(self, eps, atol_p):
        return atol_p * np.sqrt(eps)


class Depth1:
    """D1: w = p, started at 0, rate G m (`factorized_rbf.Linear(LR_C)`); its snapshots are taken at the paper unit's physical time."""
    name, q = 'D1', 1

    def __init__(self, law=LAW):
        self.lr = (law['G'] * law['m'],)
        self.lr_c = law['G'] * law['m']
        self._law = law

    def init(self, d, eps):
        return [np.zeros(d)]

    def w(self, p):
        return p[0].copy()

    def dw(self, p):
        return [np.ones_like(p[0])]

    def ddw(self, p):
        return [[np.zeros_like(p[0])]]

    def clock(self, eps):
        return lam_of(eps, self._law)

    def active(self, p, eps):
        return np.ones(len(p[0]), bool)

    def atol(self, eps, atol_p):
        return atol_p                                                   # the amplitudes are O(1) (units of the flux limit) and start at 0: no sqrt(eps) scale


READOUTS = {'P0': Paper, 'D2': Depth2, 'D1': Depth1}


def make_problem(y, sw, sigma, parts, hp):
    """The imaging problem of a window (y: the high-passed data, sw: sqrt weights; both ny x nx): dict(P, lat, c0 (d, 2), y, W)."""
    y = np.asarray(y, float)
    ny, nx = y.shape
    P = IM.PixelPSF(ny, nx, float(sigma), parts=tuple((float(a), float(k)) for a, k in np.asarray(parts, float).reshape(-1, 2)), hp=(float(hp) if hp else None))
    lat = IM.Lattice(P, h=LAT_H, margin=LAT_MARGIN)
    X, Y = np.meshgrid(lat.xs, lat.ys)
    swv = np.asarray(sw, float) * np.ones((ny, nx))
    return dict(P=P, lat=lat, c0=np.c_[X.ravel(), Y.ravel()], y=y, sw=swv, W=swv ** 2)


_SQ2, _SQ2PI = np.sqrt(2.0), np.sqrt(2.0 * np.pi)


def _axis_edges(p, c, s):
    """`PixelPSF.axis(p, c, s)` (the pixel-integrated Gaussian of width s and its centre derivative, both (len(p), len(c))) for unit pixels p = p0, p0 + 1, ...: the error
    function and the exponential are evaluated once at the len(p) + 1 pixel edges instead of twice per pixel."""
    e = np.r_[p - 0.5, p[-1] + 0.5]
    u = (e[:, None] - c[None, :]) / (_SQ2 * s)
    er, ex = erf(u), np.exp(-u ** 2)
    return 0.5 * (er[1:] - er[:-1]), -(ex[1:] - ex[:-1]) / (_SQ2PI * s)


def fast_comps(P, c):
    """`P.comps(c)` for a full-frame PixelPSF and centres c (d, 2): the separable terms (coefficient, gy, gx, dgy, dgx) of the kernels, one per Gaussian part and, with a
    high-pass, its smoothed partner. Same numbers as `imaging.PixelPSF.comps` (a test pins them), about half the special-function calls."""
    c = np.asarray(c, float)
    out = []
    for a, kap in P.parts:
        s = kap * P.sigma
        fx, fy = _axis_edges(P.px, c[:, 0], s), _axis_edges(P.py, c[:, 1], s)
        out.append((a, fy[0], fx[0], fy[1], fx[1]))
        if P.hp:
            out.append((-a, P.Sy @ fy[0], P.Sx @ fx[0], P.Sy @ fy[1], P.Sx @ fx[1]))
    return out


class _Budget(Exception):
    pass


class Net:
    """The network of a readout on an imaging problem; state z = (p_1 (d), ..., p_q (d), x (d), y (d))."""

    def __init__(self, prob, readout, eps):
        self.P, self.y, self.W = prob['P'], prob['y'], prob['W']
        self.c0 = np.asarray(prob['c0'], float)
        self.d, self.eps, self.ro, self.q = len(self.c0), float(eps), readout, readout.q
        self.nrhs = self.njac = 0
        self.deadline = None

    def init(self):
        return np.concatenate(self.ro.init(self.d, self.eps) + [self.c0[:, 0], self.c0[:, 1]])

    def split(self, z):
        d, q = self.d, self.q
        return [z[k * d:(k + 1) * d] for k in range(q)], z[q * d:(q + 1) * d], z[(q + 1) * d:]

    def w_of(self, z):
        return self.ro.w(self.split(z)[0])

    def centres(self, z):
        _, cx, cy = self.split(z)
        return np.c_[cx, cy]

    def pieces(self, z):
        """(comps, p, w, weighted residual W r, v, vx, vy): v_i = sum_p W r phi_i and its centre derivatives."""
        p, cx, cy = self.split(z)
        w = self.ro.w(p)
        comps = fast_comps(self.P, np.c_[cx, cy])
        model = 0.0
        for co, gy, gx, dgy, dgx in comps:
            model = model + co * ((gy * w) @ gx.T)
        Rw = self.W * (self.y - model)
        v = vx = vy = 0.0
        for co, gy, gx, dgy, dgx in comps:
            A_, B_ = gy.T @ Rw, dgy.T @ Rw
            v, vx, vy = v + co * (A_ * gx.T).sum(1), vx + co * (A_ * dgx.T).sum(1), vy + co * (B_ * gx.T).sum(1)
        return comps, p, w, Rw, v, vx, vy

    def neg_grad(self, z):
        """-dF/dz with unit rates: (dW/dp_l g, w d_x g, w d_y g)."""
        _, p, w, _, v, vx, vy = self.pieces(z)
        return np.concatenate([dwl * v for dwl in self.ro.dw(p)] + [w * vx, w * vy])

    def rates(self):
        """The flow's group rates as a vector of the state's shape."""
        d = self.d
        return np.concatenate([np.full(d, r) for r in self.ro.lr] + [np.full(2 * d, self.ro.lr_c)])

    def rhs(self, t, z):
        self.nrhs += 1
        if self.deadline is not None and (self.nrhs & 31) == 0 and time.time() > self.deadline:
            raise _Budget()
        return self.rates() * self.neg_grad(z)

    def loss(self, z):
        _, _, _, Rw, _, _, _ = self.pieces(z)
        return 0.5 * float(np.sum(Rw ** 2 / np.where(self.W > 0, self.W, 1.0)))

    def jac(self, t, z, dense=None):
        """d rhs / d z. For the active units (`readout.active`) the exact columns; for the others the own block of the readout parameters with the drive frozen. With every unit
        active the matrix is the exact Jacobian and is returned dense; otherwise sparse (csc). The own second derivative of the drive in the centre is a central difference."""
        self.njac += 1
        d, q, ro = self.d, self.q, self.ro
        comps, p, w, Rw, v, vx, vy = self.pieces(z)
        dw, ddw = ro.dw(p), ro.ddw(p)
        lr, lc = ro.lr, ro.lr_c
        act = np.flatnonzero(ro.active(p, self.eps))
        full = len(act) == d if dense is None else bool(dense)
        if full:
            act = np.arange(d)
        idx = np.arange(d)
        n = (q + 2) * d
        rows, cols, vals = [], [], []
        # the own blocks (drive frozen): parameters and centres against the parameters
        for k in range(q):
            for l in range(q):
                rows.append(l * d + idx)
                cols.append(k * d + idx)
                vals.append(lr[l] * ddw[l][k] * v)
            rows += [q * d + idx, (q + 1) * d + idx]
            cols += [k * d + idx, k * d + idx]
            vals += [lc * dw[k] * vx, lc * dw[k] * vy]
        R, C, V = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
        if len(act) == 0:
            return sp.coo_matrix((V, (R, C)), shape=(n, n)).tocsc()
        # the exact columns of the active units
        N = self.P.ny * self.P.nx
        Wf = self.W.ravel()
        Ph = sum(co * (gy[:, None, :] * gx[None, :, :]).reshape(N, d) for co, gy, gx, dgy, dgx in comps)
        Px = sum(co * (gy[:, None, :] * dgx[None, :, :]).reshape(N, d) for co, gy, gx, dgy, dgx in comps)
        Py = sum(co * (dgy[:, None, :] * gx[None, :, :]).reshape(N, d) for co, gy, gx, dgy, dgx in comps)
        WPh, WPx, WPy = Wf[:, None] * Ph[:, act], Wf[:, None] * Px[:, act], Wf[:, None] * Py[:, act]
        G00, Gx0, Gy0 = Ph.T @ WPh, Px.T @ WPh, Py.T @ WPh            # <W phi_i phi_j>, <W d_x phi_i phi_j>, <W d_y phi_i phi_j>     (i all, j active)
        G0x, Gxx, Gyx = Ph.T @ WPx, Px.T @ WPx, Py.T @ WPx            # ... against d_x phi_j
        G0y, Gxy, Gyy = Ph.T @ WPy, Px.T @ WPy, Py.T @ WPy            # ... against d_y phi_j
        c = np.c_[z[q * d:(q + 1) * d], z[(q + 1) * d:]]
        ex, ey = np.array([FD_C, 0.0]), np.array([0.0, FD_C])
        _, vx_px, vy_px = self.P.corr(c + ex, Rw)
        _, vx_mx, vy_mx = self.P.corr(c - ex, Rw)
        _, vx_py, vy_py = self.P.corr(c + ey, Rw)
        _, vx_my, vy_my = self.P.corr(c - ey, Rw)
        vxx, vxy, vyy = (vx_px - vx_mx) / (2 * FD_C), (vx_py - vx_my) / (2 * FD_C), (vy_py - vy_my) / (2 * FD_C)
        J = np.zeros((n, len(act) * (q + 2)))
        na = len(act)
        ar = np.arange(na)
        wa = w[act]
        for k in range(q):                                              # columns: parameter k of the active units
            cs = slice(k * na, (k + 1) * na)
            dwk = dw[k][act]
            for l in range(q):
                blk = -lr[l] * dw[l][:, None] * G00 * dwk[None, :]
                blk[act, ar] += lr[l] * ddw[l][k][act] * v[act]
                J[l * d:(l + 1) * d, cs] = blk
            bx = -lc * w[:, None] * Gx0 * dwk[None, :]
            bx[act, ar] += lc * dwk * vx[act]
            by = -lc * w[:, None] * Gy0 * dwk[None, :]
            by[act, ar] += lc * dwk * vy[act]
            J[q * d:(q + 1) * d, cs], J[(q + 1) * d:, cs] = bx, by
        for kk, (G0c, Gxc, Gyc, own_v, own_x, own_y) in enumerate(((G0x, Gxx, Gyx, vx, vxx, vxy), (G0y, Gxy, Gyy, vy, vxy, vyy))):   # columns: x, then y, of the active units
            cs = slice((q + kk) * na, (q + kk + 1) * na)
            for l in range(q):
                blk = -lr[l] * dw[l][:, None] * G0c * wa[None, :]
                blk[act, ar] += lr[l] * dw[l][act] * own_v[act]
                J[l * d:(l + 1) * d, cs] = blk
            bx = -lc * w[:, None] * Gxc * wa[None, :]
            bx[act, ar] += lc * wa * own_x[act]
            by = -lc * w[:, None] * Gyc * wa[None, :]
            by[act, ar] += lc * wa * own_y[act]
            J[q * d:(q + 1) * d, cs], J[(q + 1) * d:, cs] = bx, by
        colidx = np.concatenate([k * d + act for k in range(q + 2)])
        if full:
            out = np.zeros((n, n))
            out[:, colidx] = J
            return out
        colmask = np.zeros(n, bool)
        colmask[colidx] = True
        keep = ~colmask[C]
        nzr, nzc = np.nonzero(J)
        return sp.coo_matrix((np.concatenate([V[keep], J[nzr, nzc]]), (np.concatenate([R[keep], nzr]), np.concatenate([C[keep], colidx[nzc]]))), shape=(n, n)).tocsc()


def advance(net, z0, t0, t1, rtol=RTOL, atol_p=ATOL_P, atol_c=ATOL_C):
    """Integrate the flow of `net` from (t0, z0) to t1: Radau, then BDF from the last accepted state if Radau fails (the harness's `advance`). Returns (z, t, status, counts) with
    status 'ok', 'timeout' or 'failed: <message>'."""
    d, q = net.d, net.q
    atol = np.concatenate([np.full(q * d, net.ro.atol(net.eps, atol_p)), np.full(2 * d, atol_c)])
    cnt = dict(nfev=0, njev=0, nlu=0, steps=0, fallback=0)
    z, t = np.array(z0, float), float(t0)
    if t1 <= t0:
        return z, t, 'ok', cnt
    status = 'ok'
    for meth in (Radau, BDF):
        if meth is BDF:
            cnt['fallback'] = 1
        try:
            solver = meth(net.rhs, t, z, t1, rtol=rtol, atol=atol, jac=net.jac)
            while solver.status == 'running':
                solver.step()
                if solver.status == 'failed':
                    break
                z, t = solver.y.copy(), float(solver.t)
                cnt['steps'] += 1
        except _Budget:
            return z, t, 'timeout', cnt
        cnt['nfev'] += solver.nfev
        cnt['njev'] += solver.njev
        cnt['nlu'] += solver.nlu
        if solver.status == 'finished':
            return z, t1, 'ok', cnt
        status = 'failed: ' + str(solver.message)[:60]
    return z, t, status, cnt


RUN_FIELDS = ('status', 'secs', 'nfev', 'njev', 'nlu', 'steps', 'fallback', 't_reached', 'n_snap_done', 'clock')


def run_flow(prob, readout_name, eps, taus, budget, rtol=RTOL):
    """The gradient flow of a readout from eps to the snapshots tau (ascending; t = clock(eps) tau). Returns dict(run = RUN_FIELDS vector, W (n_snap, d), C (n_snap, d, 2),
    done (n_snap,) whether the snapshot was reached, loss (n_snap,)): a snapshot that was not reached (timeout, solver failure) holds the last accepted state."""
    t_start = time.time()
    ro = READOUTS[readout_name]()
    net = Net(prob, ro, eps)
    net.deadline = t_start + budget
    clock = ro.clock(eps)
    taus = np.asarray(taus, float)
    z, t = net.init(), 0.0
    cnt = dict(nfev=0, njev=0, nlu=0, steps=0, fallback=0)
    status, kind = 'ok', 'ok'
    Ws, Cs, done, loss = [], [], [], []
    for tau in taus:
        if status == 'ok':
            try:
                z, t, status, c_ = advance(net, z, t, clock * tau, rtol=rtol)
                for k in ('nfev', 'njev', 'nlu', 'steps'):
                    cnt[k] += c_[k]
                cnt['fallback'] = max(cnt['fallback'], c_['fallback'])
                kind = status.split(':')[0]
            except Exception as e:                                       # noqa: BLE001 - a recorded failure of the run
                status, kind = 'error', 'error'
                print('flow failed: %s: %s' % (type(e).__name__, e), flush=True)
        ok = status == 'ok' and bool(np.all(np.isfinite(z)))
        Ws.append(net.w_of(z))
        Cs.append(net.centres(z))
        done.append(float(ok))
        loss.append(net.loss(z) if np.all(np.isfinite(z)) else np.nan)
    run = dict(status=STATUS.get(kind, 2.0), secs=time.time() - t_start, nfev=cnt['nfev'], njev=cnt['njev'], nlu=cnt['nlu'], steps=cnt['steps'], fallback=cnt['fallback'],
               t_reached=t, n_snap_done=float(sum(done)), clock=clock)
    return dict(run=np.array([run[f] for f in RUN_FIELDS], float), W=np.array(Ws), C=np.array(Cs), done=np.array(done), loss=np.array(loss))


class Opt:
    """'adam': p += lr mhat / (sqrt(vhat) + eps_adam), per-coordinate second moment; 'adam_g': the gradient is first multiplied by `pre`, and the second moment is ONE scalar,
    the mean of its square (`theory3f_escape.Opt`). The caller passes the NEGATIVE gradient."""

    def __init__(self, kind, lr, beta1=0.9, beta2=0.999, eps_adam=1e-8, pre=None):
        assert kind in ('adam', 'adam_g')
        self.kind, self.lr, self.b1, self.b2, self.ea, self.pre = kind, lr, float(beta1), float(beta2), float(eps_adam), pre
        self.m = self.v = None
        self.t = 0

    def step(self, p, gneg):
        self.t += 1
        if self.pre is not None and self.kind == 'adam_g':
            gneg = self.pre * gneg
        if self.m is None:
            self.m, self.v = np.zeros_like(p), np.zeros_like(p)
        self.m = self.b1 * self.m + (1 - self.b1) * gneg
        gg = gneg ** 2 if self.kind == 'adam' else np.mean(gneg ** 2) * np.ones_like(gneg)
        self.v = self.b2 * self.v + (1 - self.b2) * gg
        mh, vh = self.m / (1 - self.b1 ** self.t), self.v / (1 - self.b2 ** self.t)
        den = np.sqrt(vh) + self.ea
        return p + self.lr * np.divide(mh, den, out=np.zeros_like(den), where=den > 0)


def run_discrete(prob, kind, eps, steps, lr, budget, blow=1e9):
    """The paper's unit trained by a discrete optimizer ('adam' or 'adam_g') from eps, snapshots after the step counts `steps` (ascending integers). Returns the dict of
    `run_flow` (t_reached = the steps done, clock = 1)."""
    t_start = time.time()
    net = Net(prob, Paper(), eps)
    opt = Opt(kind, lr, pre=net.rates() if kind == 'adam_g' else None)
    z = net.init()
    steps = [int(s) for s in steps]
    Ws, Cs, done, loss = [], [], [], []
    status, k = 'ok', 0
    for target in steps:
        while status == 'ok' and k < target:
            z = opt.step(z, net.neg_grad(z))
            k += 1
            if (k & 255) == 0:
                if not np.all(np.isfinite(z)) or np.abs(z).max() > blow:
                    status = 'diverged'
                elif time.time() - t_start > budget:
                    status = 'timeout'
        fin = bool(np.all(np.isfinite(z)))
        if status == 'ok' and not fin:
            status = 'diverged'
        Ws.append(net.w_of(z) if fin else np.full(net.d, np.nan))
        Cs.append(net.centres(z) if fin else np.full((net.d, 2), np.nan))
        done.append(float(status == 'ok'))
        loss.append(net.loss(z) if fin else np.nan)
    run = dict(status=STATUS[status], secs=time.time() - t_start, nfev=k, njev=0, nlu=0, steps=k, fallback=0, t_reached=float(k), n_snap_done=float(sum(done)), clock=1.0)
    return dict(run=np.array([run[f] for f in RUN_FIELDS], float), W=np.array(Ws), C=np.array(Cs), done=np.array(done), loss=np.array(loss))


ARMS = ('A', 'B', 'C', 'D', 'E')
ARM_READOUT = {'A': 'P0', 'B': 'D2', 'C': 'D1'}
ARM_OPT = {'D': 'adam', 'E': 'adam_g'}


def run_arm(prob, arm, eps, grid, budget, lr=1e-3):
    """One arm on one problem: A, B, C the flow of P0, D2, D1 to the snapshots `grid` (tau); D, E the discrete optimizers to the snapshots `grid` (step counts)."""
    if arm in ARM_READOUT:
        return run_flow(prob, ARM_READOUT[arm], eps, grid, budget)
    return run_discrete(prob, ARM_OPT[arm], eps, grid, lr, budget)


def job(spec):
    """One remote job: spec = dict(key, arm, eps, grid, budget, lr, y, sw, sigma, parts, hp). Returns (key, result dict of `run_arm`)."""
    prob = make_problem(spec['y'], spec['sw'], spec['sigma'], spec['parts'], spec['hp'])
    try:
        res = run_arm(prob, spec['arm'], spec['eps'], spec['grid'], spec['budget'], spec.get('lr', 1e-3))
    except Exception as e:                                               # noqa: BLE001
        print('job %s failed: %s: %s' % (spec['key'], type(e).__name__, e), flush=True)
        n, d = len(spec['grid']), len(prob['c0'])
        run = np.zeros(len(RUN_FIELDS))
        run[0] = STATUS['error']
        res = dict(run=run, W=np.full((n, d), np.nan), C=np.full((n, d, 2), np.nan), done=np.zeros(n), loss=np.full(n, np.nan))
    return spec['key'], res


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# job files (small derived arrays only: the normalised window, its weights, the kernel's parameters) and their results
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
W_STORE = {'A': 1e-7, 'B': 0.05, 'C': 0.05, 'D': 0.05, 'E': 0.05}     # a snapshot keeps the units with |w| above this (arm A in full precision down to the dormant units)
JOB_ARRAYS = ('y', 'sw', 'parts', 'grid')
JOB_SCALARS = ('arm', 'eps', 'budget', 'lr', 'sigma', 'hp')


def save_jobs(path, specs):
    """Write job specs (dicts of `job`) to an npz (no pickle): 'keys' and, per job i, 'i|<name>'."""
    arrays = {'keys': np.array([s['key'] for s in specs])}
    for i, s in enumerate(specs):
        for k in JOB_ARRAYS:
            arrays['%d|%s' % (i, k)] = np.asarray(s[k], float)
        arrays['%d|arm' % i] = np.array([s['arm']])
        for k in ('eps', 'budget', 'lr', 'sigma', 'hp'):
            arrays['%d|%s' % (i, k)] = np.array([float(s[k])])
    np.savez_compressed(path, **arrays)


def load_jobs(path):
    with np.load(path, allow_pickle=False) as z:
        keys = [str(k) for k in z['keys']]
        out = []
        for i, key in enumerate(keys):
            s = dict(key=key, arm=str(z['%d|arm' % i][0]))
            for k in JOB_ARRAYS:
                s[k] = z['%d|%s' % (i, k)]
            for k in ('eps', 'budget', 'lr', 'sigma', 'hp'):
                s[k] = float(z['%d|%s' % (i, k)][0])
            out.append(s)
    return out


def pack_result(key, arm, res):
    """The stored form of a run: 'R|<key>|run', 'done', 'loss', 'grid' is the caller's; the snapshots as concatenated sparse lists 'n' (units kept per snapshot), 'idx', 'w', 'c'."""
    thr = W_STORE[arm]
    dt = np.float64 if arm == 'A' else np.float32
    n, idx, w, c = [], [], [], []
    for Wk, Ck in zip(res['W'], res['C']):
        f = np.flatnonzero(np.isfinite(Wk) & (np.abs(Wk) > thr))
        n.append(len(f))
        idx.append(f)
        w.append(Wk[f])
        c.append(Ck[f])
    pre = 'R|%s|' % key
    return {pre + 'run': np.asarray(res['run'], float), pre + 'done': np.asarray(res['done'], float), pre + 'loss': np.asarray(res['loss'], float), pre + 'n': np.array(n, np.int32),
            pre + 'idx': np.concatenate(idx).astype(np.int16) if idx else np.zeros(0, np.int16), pre + 'w': np.concatenate(w).astype(dt) if w else np.zeros(0, dt),
            pre + 'c': np.concatenate(c).astype(dt).reshape(-1, 2) if c else np.zeros((0, 2), dt)}


def unpack_result(raw, key):
    """dict(run, done, loss, snaps = [(idx, w, c), ...]) of a stored run (None if absent)."""
    pre = 'R|%s|' % key
    if pre + 'run' not in raw:
        return None
    n = np.asarray(raw[pre + 'n'], int)
    ends = np.cumsum(n)
    snaps = [(np.asarray(raw[pre + 'idx'][e - k:e], int), np.asarray(raw[pre + 'w'][e - k:e], float), np.asarray(raw[pre + 'c'][e - k:e], float)) for k, e in zip(n, ends)]
    return dict(run=np.asarray(raw[pre + 'run'], float), done=np.asarray(raw[pre + 'done'], float), loss=np.asarray(raw[pre + 'loss'], float), snaps=snaps)


def _job_packed(spec):
    key, res = job(spec)
    return pack_result(key, spec['arm'], res)


def work(jobfile, outfile, procs=1):
    """Run every job of a job file and write the packed results (after every job, so that an interrupted worker keeps what is done)."""
    specs = load_jobs(jobfile)
    done = {}
    if os.path.exists(outfile):
        with np.load(outfile, allow_pickle=False) as z:
            done = {k: z[k] for k in z.files}
    todo = [s for s in specs if 'R|%s|run' % s['key'] not in done]
    print('%d jobs, %d to run, %d processes' % (len(specs), len(todo), procs), flush=True)

    def flush():
        tmp = outfile[:-4] + '.tmp.npz'
        np.savez_compressed(tmp, **done)
        os.replace(tmp, outfile)
    t0 = time.time()
    if procs > 1:
        import multiprocessing as mp
        with mp.Pool(procs) as pool:
            for out in pool.imap_unordered(_job_packed, todo):
                done.update(out)
                flush()
                k = [q for q in out if q.endswith('|run')][0]
                print('%7.0f s  %s  status %g  %.0f s  steps %g' % (time.time() - t0, k, out[k][0], out[k][1], out[k][5]), flush=True)
    else:
        for s in todo:
            out = _job_packed(s)
            done.update(out)
            flush()
            k = 'R|%s|run' % s['key']
            print('%7.0f s  %s  status %g  %.0f s  steps %g' % (time.time() - t0, k, out[k][0], out[k][1], out[k][5]), flush=True)
    flush()
    return len(todo)


if __name__ == '__main__':
    work(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 1)
