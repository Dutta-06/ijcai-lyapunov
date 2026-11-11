"""
Round theory3 / direction E: the limit path (Theorem K recursion) and the network flow for a GENERAL Gram, in the continuum, centres in R^n.

Generic-kernel version of the Gaussian code of exp_theory_t1.py (limit_path, crit_block, simulate): the kernel enters only through Gamma(|c - c'|) and Gamma' (theory3e_kernels.Gram:
fG, fdG in sigma-equivalent units, -Gamma''(0) = 1/2).  The unit (the repository's P0 law, Lemma S1: closed forms in the dual coordinate x) is taken from exp_theory_t1.Unit unchanged.

  limit_path        the recursion of Theorem K: piecewise linear duals (wall 1 in tau = t / lambda), entrant = argmin (1 - sgn(p) u)/|p|, assigned to the nearest unfitted feature, block = the
                    critical point of the fit over the entrants (Newton from the features), with the finite checks (a) unique minimiser with the feature's sign, (b) |u_j| < 1 at every entry
  simulate_network  the continuum network  dx_j/dt = (G/a_max)(eps + (2 v_s - eps) tanh^2 x_j) g_j,  dc_j/dt = mu_c w_j h_j  (LSODA), centres in R^n
  compare           formed set, entry times against lambda tau_k, final fit error, dormant coefficients
  inner_run         the inner problem of one entry layer (block at its critical point, entrant started at v_0 at offset d_0) integrated for a given positive mobility: box invariance and
                    convergence (the independent numerical test of Proposition G)
"""
import os
import sys

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import exp_theory_t1 as T1            # noqa: E402   (the unit: Unit, LAW, THETA_REL, H_LAT; unchanged)

H_LAT = T1.H_LAT                       # lattice spacing in widths (0.1 / 0.12)


# ----------------------------------------------------------------------------------------------------------------------------------
# kernels on point sets
# ----------------------------------------------------------------------------------------------------------------------------------
def _dist(a, b):
    """Pairwise distances and unit vectors from b to a (a: (p, n), b: (q, n)): |a_i - b_j| and (a_i - b_j)/|.| (0 at coincidence)."""
    d = a[:, None, :] - b[None, :, :]
    r = np.sqrt((d ** 2).sum(-1))
    u = d / np.where(r > 0, r, 1.0)[..., None]
    return r, u


def gam_mat(kern, a, b):
    r, _ = _dist(a, b)
    return kern.fG(r)


def grad_mat(kern, a, b):
    """(p, q, n): gradient with respect to a_i of Gamma(|a_i - b_j|)."""
    r, u = _dist(a, b)
    return kern.fdG(r)[..., None] * u


def as_pts(x, n=None):
    x = np.asarray(x, float)
    return x[:, None] if x.ndim == 1 else x


# ----------------------------------------------------------------------------------------------------------------------------------
# the limit object
# ----------------------------------------------------------------------------------------------------------------------------------
def crit_block(kern, cs, As, feats, init=None):
    """Critical point (w, c) of the fit over len(feats) units started at (A_f, c*_f): F_w = 0, F_c = 0."""
    k, n = len(feats), cs.shape[1]
    f = np.asarray(feats)

    def Fg(z):
        w, c = z[:k], z[k:].reshape(k, n)
        g = gam_mat(kern, c, cs) @ As - gam_mat(kern, c, c) @ w
        Gc = grad_mat(kern, c, c)                                     # (k, k, n); i = l entries vanish (Gamma'(0) = 0 or the unit vector 0)
        h = (grad_mat(kern, c, cs) * As[None, :, None]).sum(1) - (Gc * w[None, :, None]).sum(1)
        return np.r_[g, (w[:, None] * h).ravel()]
    z0 = np.r_[As[f], cs[f].ravel()] if init is None else init
    s = root(Fg, z0, method='hybr', tol=1e-14)
    return s.x[:k], s.x[k:].reshape(k, n), float(np.abs(Fg(s.x)).max()) < 1e-10


def limit_path(kern, cs, As, c0):
    """The recursion of Theorem K for the continuum Gram with equal walls.  cs: (K, n) feature centres, As: (K,) amplitudes, c0: (J, n) starting centres.
    Returns events (kernel, feature, tau, relative gap of the entry times, offset d0, the entrant's rate p_e, sign_ok, block w, c, block ok), breakpoints tau_br, duals U_br, rates P."""
    cs, c0, As = as_pts(cs), as_pts(c0), np.asarray(As, float)
    K, J = len(cs), len(c0)
    u = np.zeros(J)
    tau = 0.0
    active, assigned = [], []
    w, c = np.zeros(0), np.zeros((0, cs.shape[1]))
    ev, unfit = [], list(range(K))
    tau_br, U_br, P = [0.0], [u.copy()], []
    G0 = gam_mat(kern, c0, cs)
    for k in range(K):
        p = G0 @ As - (gam_mat(kern, c0, c) @ w if len(w) else 0.0)
        p = np.where(np.isin(np.arange(J), active), 0.0, p)
        dt = np.full(J, np.inf)
        for j in range(J):
            if j in active or abs(p[j]) < 1e-14:
                continue
            dt[j] = (1 - np.sign(p[j]) * u[j]) / abs(p[j])
        order = np.argsort(dt)
        j1 = int(order[0])
        gap = (dt[order[1]] - dt[order[0]]) / dt[order[0]] if np.isfinite(dt[order[1]]) else np.inf
        P.append(p.copy())
        u = u + dt[j1] * p
        tau += dt[j1]
        f = min(unfit, key=lambda m: np.linalg.norm(c0[j1] - cs[m]))
        unfit.remove(f)
        active.append(j1)
        assigned.append(f)
        w, c, okb = crit_block(kern, cs, As, assigned)
        u[j1] = np.sign(p[j1])
        tau_br.append(tau)
        U_br.append(u.copy())
        ev.append(dict(k=k + 1, kernel=j1, feature=f, tau=tau, gap=float(gap), d0=float(np.linalg.norm(c0[j1] - cs[f])), p_e=float(abs(p[j1])),
                       sign_ok=bool(np.sign(p[j1]) == np.sign(As[f])), w=w.copy(), c=c.copy(), block_ok=okb))
    pend = G0 @ As - gam_mat(kern, c0, c) @ w
    P.append(np.where(np.isin(np.arange(J), active), 0.0, pend))
    ent = set(active)
    non = [j for j in range(J) if j not in ent]
    mu = float(min(1 - np.abs(np.array(U_br)[:, non]).max(axis=1))) if non else np.inf        # margin (b): 1 - max_j |u_j| over the whole path
    return dict(events=ev, tau_br=np.array(tau_br), U_br=np.array(U_br), P=np.array(P), active=active, K=K, J=J, margin=mu,
                final_w=w, final_c=c, resid=float(np.abs(pend[non]).max()) if non else 0.0)


def make_target(n, K, sep, R, seed, margin=12.0, h=H_LAT):
    """Features on a line (direction random in R^n), separations about `sep` (offsets +-h/2 from the lattice), amplitudes (1, R) or (1, sqrt R, R) x 0.7, random order and signs; the
    dictionary: the cubic lattice of spacing h over the features' bounding box +- margin at a random phase (n = 1) / restricted to within `margin` of the features (n >= 2)."""
    rng = np.random.default_rng(seed)
    e = np.array([1.0]) if n == 1 else rng.normal(size=n)
    e = e / np.linalg.norm(e) if n > 1 else e
    t = np.arange(K) * sep + rng.uniform(-h / 2, h / 2, K)
    t = t - t.mean()
    cs = t[:, None] * e[None, :]
    a = np.array([1.0]) if K == 1 else (np.array([1.0, R]) if K == 2 else np.array([1.0, np.sqrt(R), R]))
    a = a[rng.permutation(K)] * 0.7
    As = a * rng.choice([-1.0, 1.0], K)
    lo_, hi_ = cs.min(0) - margin, cs.max(0) + margin
    ax = [np.arange(lo_[i], hi_[i], h) + rng.uniform(0, h) for i in range(n)]
    grid = np.stack(np.meshgrid(*ax, indexing='ij'), -1).reshape(-1, n)
    dmin = np.sqrt(((grid[:, None, :] - cs[None, :, :]) ** 2).sum(-1)).min(1)
    return cs, As, grid[dmin <= margin]


# ----------------------------------------------------------------------------------------------------------------------------------
# the network
# ----------------------------------------------------------------------------------------------------------------------------------
def simulate_network(kern, cs, As, c0, eps, path, tau_end=None, rtol=1e-10, nglob=3000, nloc=1500):
    """Integrate the continuum network from x = 0 (alpha = 0), nu^2 = eps, c = c0.  Returns (sol, U); y = (x_1..x_J, c_1..c_J flattened (J, n))."""
    cs, c0, As = as_pts(cs), as_pts(c0), np.asarray(As, float)
    J, n = c0.shape
    U = T1.Unit(eps)
    lam, muc = U.lam, U.mu_c
    ev = path['events']
    tau_end = (1.05 * ev[-1]['tau'] + 1.0) if tau_end is None else tau_end
    Tfin = lam * tau_end
    tk = lam * path['tau_br'][1:]

    def f(t, z):
        x, c = z[:J], z[J:].reshape(J, n)
        w = (U.amax / U.m) * np.tanh(x) * (eps + 2 * U.vs * np.sinh(x) ** 2)
        g = gam_mat(kern, c, cs) @ As - gam_mat(kern, c, c) @ w
        h = (grad_mat(kern, c, cs) * As[None, :, None]).sum(1) - (grad_mat(kern, c, c) * w[None, :, None]).sum(1)
        return np.r_[(U.G / U.amax) * (eps + (2 * U.vs - eps) * np.tanh(x) ** 2) * g, (muc * w[:, None] * h).ravel()]
    loc = [np.linspace(max(0, t - 600), min(Tfin, t + 900), nloc) for t in tk]
    ts = np.unique(np.concatenate([np.linspace(0, Tfin, nglob)] + loc))
    sol = solve_ivp(f, (0, Tfin), np.r_[np.zeros(J), c0.ravel()], method='LSODA', rtol=rtol, atol=np.r_[np.full(J, 1e-14), np.full(J * n, 1e-12)], t_eval=ts)
    return sol, U


def wfun(U, x, eps):
    return (U.amax / U.m) * np.tanh(x) * (eps + 2 * U.vs * np.sinh(x) ** 2)


def compare(cs, As, c0, eps, path, sol, U):
    """The network against the limit path: does it form exactly the entrants (|w| > theta at the end), exact fit, entry times against lambda tau_k."""
    cs, c0, As = as_pts(cs), as_pts(c0), np.asarray(As, float)
    J, n = c0.shape
    ev = path['events']
    x, c = sol.y[:J, :], sol.y[J:, :].reshape(J, n, -1)
    w = wfun(U, x, eps)
    theta = min(T1.THETA_REL * np.abs(As).max(), 0.1 * np.abs(As).min())
    ent = [e['kernel'] for e in ev]
    formed = sorted(np.flatnonzero(np.abs(w[:, -1]) > theta).tolist())
    ok = (formed == sorted(ent)) and sol.status == 0
    err_w = max(abs(abs(w[e['kernel'], -1]) - abs(As[e['feature']])) / abs(As[e['feature']]) for e in ev)
    err_c = max(float(np.linalg.norm(c[e['kernel'], :, -1] - cs[e['feature']])) for e in ev)
    # entry times: first time |w_e| > theta, against lambda tau_k
    shifts = []
    for e in ev:
        idx = np.flatnonzero(np.abs(w[e['kernel']]) > theta)
        t_net = sol.t[idx[0]] if len(idx) else np.nan
        shifts.append(float(t_net / (U.lam * e['tau']) - 1.0))
    non = [j for j in range(J) if j not in ent]
    dorm = float(np.abs(w[non, :]).max()) if non else 0.0
    drift = float(np.linalg.norm(c[non, :, -1] - c0[non], axis=1).max()) if non else 0.0
    return dict(formed_ok=bool(ok), exact=bool(ok and err_w < 2e-2 and err_c < 0.05), err_w=float(err_w), err_c=float(err_c), shifts=shifts, dormant_max=dorm, dormant_drift=drift,
                lam_sqrt_eps=float(U.lam * np.sqrt(eps)), order_ok=bool(np.all(np.argsort([e['tau'] for e in ev]) == np.argsort(np.array([np.nanmax([sol.t[np.flatnonzero(np.abs(w[e['kernel']]) > theta)[0]] if len(np.flatnonzero(np.abs(w[e['kernel']]) > theta)) else np.nan]) for e in ev])))))


# ----------------------------------------------------------------------------------------------------------------------------------
# the inner problem of one entry layer (independent numerical test of Proposition G)
# ----------------------------------------------------------------------------------------------------------------------------------
def inner_run(kern, K, Delta, d0, R, mob='v43', muc=1.0, v0=1e-2, seed=0, tmax=4000.0, n=1):
    """Features at 0, Delta, 2 Delta ... (amplitudes 1, R, 1, ... positive) along a line in R^n; units 1..K-1 start at their critical points (the block of the previous layers: computed by crit_block from the
    features), the entrant (unit K) at the offset d0 from the last feature (towards the previous one: the worst for the coupling) with coefficient v0; the flow
        v_i' = -M(v_i) dF/dv_i,   c_i' = -muc dF/dc_i        (M(v) = 1, v, v^{4/3}, 1 + v^2: all positive)
    Returns the trajectory's max normalised offset max|d_i|, the final |v - A|, |d|, monotonicity of the entrant's coefficient and the time series (for the box test)."""
    rng = np.random.default_rng(seed)
    a = np.ones(K)
    a[1::2] = R
    cs = (np.arange(K) * Delta)[:, None] * np.eye(n)[0][None, :]
    As = a.copy()
    if K > 1:
        w0, c0_, _ = crit_block(kern, cs, As, list(range(K - 1)))
    else:
        w0, c0_ = np.zeros(0), np.zeros((0, n))
    ent_c = cs[K - 1] - d0 * np.eye(n)[0]
    z0 = np.r_[w0, v0, c0_.ravel(), ent_c]
    Mf = {'1': lambda v: np.ones_like(v), 'v': lambda v: v, 'v43': lambda v: v ** (4 / 3), 'q': lambda v: 1 + v * v}[mob]

    def f(t, z):
        v, c = z[:K], z[K:].reshape(K, n)
        g = gam_mat(kern, c, cs) @ As - gam_mat(kern, c, c) @ v                    # -dF/dv
        h = (grad_mat(kern, c, cs) * As[None, :, None]).sum(1) - (grad_mat(kern, c, c) * v[None, :, None]).sum(1)      # grad_c of sum_k A_k Gamma - sum_l v_l Gamma  (centre force / v)
        return np.r_[Mf(v) * g, (muc * v[:, None] * h).ravel()]
    ts = np.linspace(0, tmax, 4001)
    sol = solve_ivp(f, (0, tmax), z0, method='LSODA', rtol=1e-10, atol=1e-13, t_eval=ts)
    v = sol.y[:K]
    c = sol.y[K:].reshape(K, n, -1)
    d = np.linalg.norm(c - cs[:, :, None], axis=1)
    return dict(t=sol.t, v=v, d=d, max_d=float(d.max()), final_v_err=float(np.abs(v[:, -1] - As).max()), final_d=float(d[:, -1].max()),
                ent_monotone=bool(np.all(np.diff(v[K - 1]) >= -1e-12)), vmin_after=float(v[:, 1:].min()) if K > 1 else float(v[:, 1:].min()),
                v_max=float(v.max()), converged=bool(sol.status == 0))


def box_run(kern, K, R, Delta, r, vlo, vbar, seed, mob='v43', muc=1.0, n=1, tmax=2500.0, entrant_low=True, adv=False):
    """Proposition G(i), (ii), (iii) tested on a random initial point of the box B = {|d_i| <= r, |w_i| in [v_lo, v_bar]} (the last unit's coefficient from (0, v_lo] when entrant_low): K features on a line at
    spacing Delta, amplitudes (1, R, 1, ...) with random signs, one unit per feature, mobility law `mob` (all positive), centre mobility muc.  The trajectory has to stay in B (|d_i| <= r, |w_i| <= v_bar,
    |w_i| >= min(start, v_lo)), converge to the target (w = A, d = 0), and the entrant's coefficient has to increase while it is below v_lo."""
    rng = np.random.default_rng(seed)
    a = np.ones(K)
    a[1::2] = R
    sg = rng.choice([-1.0, 1.0], K)
    As = a * sg
    e1 = np.eye(n)[0]
    cs = (np.arange(K) * Delta)[:, None] * e1[None, :]
    dirs = rng.normal(size=(K, n))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    rad = r * rng.uniform(0, 1, K) ** (1.0 / n) * 0.999
    c0 = cs + dirs * rad[:, None]
    v0 = rng.uniform(vlo, vbar, K)
    if adv:                                                      # adversarial corner: offsets on the boundary along +-e1, coefficients at the faces
        rad = np.full(K, r * 0.999)
        dirs = np.zeros((K, n))
        dirs[:, 0] = rng.choice([-1.0, 1.0], K)
        c0 = cs + dirs * rad[:, None]
        v0 = rng.choice([vlo, vbar], K)
    if entrant_low:
        v0[-1] = rng.uniform(1e-3, vlo)
    w0 = sg * v0
    Mf = {'1': lambda v: np.ones_like(v), 'v': lambda v: v, 'v43': lambda v: v ** (4 / 3), 'q': lambda v: 1 + v * v}[mob]

    def f(t, z):
        w, c = z[:K], z[K:].reshape(K, n)
        g = gam_mat(kern, c, cs) @ As - gam_mat(kern, c, c) @ w
        h = (grad_mat(kern, c, cs) * As[None, :, None]).sum(1) - (grad_mat(kern, c, c) * w[None, :, None]).sum(1)
        return np.r_[Mf(np.abs(w)) * g, (muc * w[:, None] * h).ravel()]
    ts = np.linspace(0, tmax, 3001)
    sol = solve_ivp(f, (0, tmax), np.r_[w0, c0.ravel()], method='LSODA', rtol=1e-10, atol=1e-13, t_eval=ts)
    w = sol.y[:K]
    c = sol.y[K:].reshape(K, n, -1)
    d = np.linalg.norm(c - cs[:, :, None], axis=1)
    v = np.abs(w)
    lowbound = np.minimum(v0, vlo)
    ent = v[K - 1]
    below = ent < vlo
    # entrant monotone while below v_lo: the increments of the part of the series before the first time it reaches v_lo
    first = np.argmax(~below) if (~below).any() else len(ent)
    mono = bool(np.all(np.diff(ent[:max(first, 1)]) >= -1e-12))
    return dict(max_d=float(d.max()), r=float(r), in_d=bool(d.max() <= r * (1 + 1e-9)), in_hi=bool(v.max() <= vbar * (1 + 1e-9)), in_lo=bool(np.all(v >= lowbound[:, None] * (1 - 1e-9))),
                final_w_err=float(np.abs(w[:, -1] - As).max()), final_d=float(d[:, -1].max()), mono=mono, converged=bool(sol.status == 0 and np.abs(w[:, -1] - As).max() < 1e-6 and d[:, -1].max() < 1e-6))
