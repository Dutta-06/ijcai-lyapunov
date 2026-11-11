"""
The entering orbit of the coupled entry layer (round theory3, direction Y): the layer flow from the dormant state, the box constants, the explicit
hitting-time and convergence bounds, and the checks against the integrated flows.

Setting (THEORY_NOTE.md Prop. A, THEORY_T1.md sections 2-4; the statements and proofs of this round are in its NOTE.md, Theorems E1-E4).  Continuum Gram
<phi(c) phi(c')> = Z Gamma(c - c') (Z = 1, time t' = Z t), Gamma(d) = exp(-d^2/4) (lengths in widths), n = 1, K Gaussian features at separation >= Delta.
The *layer flow* of the k-th entry is the flow of the k-1 formed units of the block and the entrant e (every other kernel dormant: Omega = 0), in the
coordinates x_i = alpha_i/a_max and c_i, on the critical curve of eps >= 0:
        dx_i/dt = (G/a_max) (eps + (2 v_s - eps) tanh^2 x_i) g_i,      dc_i/dt = mu_c w_i h_i,      w_i = W_eps(x_i) = (a_max/m) tanh x (eps + 2 v_s sinh^2 x),
        g_i = sum_m A_m Gamma(c_i - c*_m) - sum_j w_j Gamma(c_i - c_j),   h_i = sum_m A_m Gamma'(c_i - c*_m) - sum_{j != i} w_j Gamma'(c_i - c_j).
eps > 0 is the network's own flow of those units (the entrant at x = 0, i.e. w = 0 exactly, the dormant state); eps = 0 is the limit object's layer, whose
dormant state {x_e = 0} consists of equilibria and is left only asymptotically (w_e ~ |t|^-3 as t -> -infinity).  The block starts at the critical point
of the previous layer, the entrant at its lattice kernel at the offset d0 from its feature.

What is computed (everything is a finite computation; nothing is fitted):
  * the box data of Prop. A with the margins of the relaxation set (exp_theory_t1.box_t1, ball arithmetic) for each layer, the radius r > d0 chosen as the
    smallest grid value for which (L), (C), (N-hat) hold ((N-hat): the least eigenvalue of the own block at d-hat exceeds the cross-Hessian bound);
  * the explicit constants of Theorem E3: the level w_# = a_min v_lo / 2, the times T_2, T_3 (balls), the rate rho = m_lo nu, the prefactor sqrt(H_up/nu);
  * the integrated layer flows: eps = 1e-6 ... 1e-12 from the dormant state x_e = 0, and eps = 0 from w_e = v_0 a_min, v_0 = 1e-3 ... 1e-12, time-aligned at the
    first crossing of w_#; the checks of Theorems E1-E4 on every flow (the box, the monotone entrant, the hitting time of the relaxation set against T_2 + T_3, the
    exponential decay against rho, the time to enter D*) and the aligned distances between the flows (Cauchy in v_0; convergence in eps).
"""
import os
import sys

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T          # noqa: E402
import theory_t1b_wide as W        # noqa: E402

try:
    from flint import arb
    HAVE_ARB = True
except Exception:                      # pragma: no cover
    HAVE_ARB = False

LAW = T.LAW
SIG = T.SIG
MC = LAW['G'] * LAW['m'] / SIG ** 2                       # centre mobility (widths), t' units
VS = LAW['gamma'] * LAW['a_max'] ** 2 / (4 * LAW['G'])
VERIFIED = (0, 1, 2, 3, 4, 5, 6, 8)                       # the 8 targets of theory3n (indices into W.Q_IDX)
D_TARGETS = [tuple(T.C_TARGETS[W.Q_IDX[n]]) for n in VERIFIED]
R_SEED = 20260111
EPS_LIST = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
V0_LIST = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
V0_REF = 1e-14
S2 = 4000.0                                               # duration of the second stage (after the level crossing)
DT = 0.25                                                 # sampling step of the second stage
R_GRID = np.arange(0.2, 1.31, 0.05)
DHAT_FRAC, VSIG = 0.05, 0.01                              # the relaxation set's margins (as in Theory T1)
CROSS_RTOL = 1e-12


def random_targets(n=8):
    """n further targets (K in {2, 3}, separation 5.5-9.5 widths, amplitude ratio 1-2.4, seeds 7001...), drawn from a fixed generator."""
    rng = np.random.default_rng(R_SEED)
    out = []
    for i in range(n):
        K = int(rng.choice([2, 3]))
        out.append((K, round(float(rng.uniform(5.5, 9.5)), 2), round(float(rng.uniform(1.0, 2.4)), 2), 7001 + i))
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# the unit on the critical curve of eps >= 0 (Lemma S1; THEORY_T1.md Lemma 2.1)
# ----------------------------------------------------------------------------------------------------------------------------------
def W_of_x(x, eps):
    x = np.asarray(x, float)
    return (LAW['a_max'] / LAW['m']) * np.tanh(x) * (eps + 2 * VS * np.sinh(x) ** 2)


def x_of_w(w, eps):
    """The inverse of W_eps (odd, increasing), by bisection."""
    s = 1.0 if w >= 0 else -1.0
    w = abs(float(w))
    lo, hi = 0.0, 30.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if float(W_of_x(mid, eps)) < w:
            lo = mid
        else:
            hi = mid
    return s * 0.5 * (lo + hi)


def M_of_x(x, eps):
    """Mobility dw/d eta = (G/m) q (q + 4 v_s sinh^2 x), q = eps + (2 v_s - eps) tanh^2 x."""
    q = eps + (2 * VS - eps) * np.tanh(x) ** 2
    return (LAW['G'] / LAW['m']) * q * (q + 4 * VS * np.sinh(x) ** 2)


def M_of_w(w, eps):
    return float(M_of_x(x_of_w(w, eps), eps))


# ----------------------------------------------------------------------------------------------------------------------------------
# layers of a target
# ----------------------------------------------------------------------------------------------------------------------------------
def layers_of(target):
    """The entry layers of the limit path of a target (K, sep, R, seed): for the k-th entry the block (critical point of the previous layer), the entrant's kernel and
    offset, the layer's critical point, the units' features and signs."""
    K, sep, R, seed = target
    cs, As, c0 = T.make_target(K, sep, R, seed)
    prep = T.prepare(cs, As, c0)
    ev = prep['L']['events']
    out = []
    for k in range(K):
        feats = [e['feature'] for e in ev[:k + 1]]
        out.append(dict(target=tuple(target), k=k, K=K, cs=cs, As=As, feats=feats, sgn=np.sign(As[feats]), a=np.abs(As[feats]),
                        cb=np.asarray(ev[k - 1]['c']) if k > 0 else np.zeros(0), wb=np.asarray(ev[k - 1]['w']) if k > 0 else np.zeros(0),
                        ce=float(c0[ev[k]['kernel']]), d0=float(ev[k]['d0']), what=np.asarray(ev[k]['w']), chat=np.asarray(ev[k]['c']),
                        R=prep['R'], Delta=prep['Delta'], amin=prep['amin'], sign_ok=bool(ev[k]['sign_ok']), mu_min=prep['mu_min']))
    return out


STRESS_DELTAS = (4.9, 5.2, 5.6)
STRESS_RS = (1.0, 2.4)
STRESS_D0 = 0.40


def stress_layers():
    """Second layers of K = 2 targets built directly at separations just above Prop. A's threshold (Delta = 4.9, 5.2, 5.6 widths, R = 1, 2.4, the weak or the strong feature first, the entrant
    at the offset 0.40 TOWARD the formed unit, equal and opposite signs at Delta = 5.2): the block is the critical point of the first unit's fit to the whole target, the layer's critical point
    is the target.  Target label (2, Delta, R, code), code = 100 (weak first) + 10 (opposite signs) + 1."""
    out = []
    for Delta in STRESS_DELTAS:
        for R in STRESS_RS:
            for weak in (True, False):
                for opp in ((False, True) if Delta == 5.2 else (False,)):
                    a = 0.7 * np.array([1.0, R] if weak else [R, 1.0])
                    As = a * np.array([1.0, -1.0 if opp else 1.0])
                    cs = np.array([-Delta / 2, Delta / 2])
                    wb, cb = T.crit_block(cs, As, [0])
                    out.append(dict(target=(2, Delta, R, 100 * int(weak) + 10 * int(opp) + 1), k=1, K=2, cs=cs, As=As, feats=[0, 1], sgn=np.sign(As), a=np.abs(As), cb=np.asarray(cb),
                                    wb=np.asarray(wb), ce=float(cs[1] - STRESS_D0), d0=STRESS_D0, what=As.copy(), chat=cs.copy(), R=float(R), Delta=float(Delta), amin=0.7, sign_ok=True,
                                    mu_min=float('nan')))
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# the box constants and the explicit bounds (Theorem E3)
# ----------------------------------------------------------------------------------------------------------------------------------
def layer_box(L):
    """The smallest grid radius r > d0 for which the box conditions (L), (C), (N-hat) hold with Omega = 0 (ball arithmetic in box_t1); the dict of its constants, or
    None.  (N-hat) is nu = lambda_own(d-hat) - H_c > 0 and the positivity of the lower faces of the relaxation set; the Gershgorin margins of Lemma 4.1 are not needed."""
    K = L['K']
    Delta = L['Delta'] if np.isfinite(L['Delta']) else 0.0
    for r in R_GRID:
        r = float(r)
        if r <= L['d0'] + 1e-9:
            continue
        b = T.box_t1(K, L['R'], Delta, r, 0.0, DHAT_FRAC, VSIG)
        if 'nu' in b and b['nu'] > 0 and b['vD_lo'] > 0:
            b = dict(b)
            b['r'] = r
            return b
    return None


def _arb_W(x, eps):
    x = arb(x)
    return arb(LAW['a_max']) / arb(LAW['m']) * x.tanh() * (arb(eps) + 2 * arb(VS) * x.sinh() ** 2)


def _arb_M(x, eps):
    x = arb(x)
    e = arb(eps)
    q = e + (2 * arb(VS) - e) * x.tanh() ** 2
    return arb(LAW['G']) / arb(LAW['m']) * q * (q + 4 * arb(VS) * x.sinh() ** 2)


def M_lower(w, eps):
    """A certain lower bound of M_eps(w): M increases with x, so M at x_lo with W(x_lo) < w (checked in ball arithmetic) is below M(W^-1(w))."""
    x = x_of_w(w, eps)
    f = 1.0 - 1e-9
    for _ in range(80):
        xl = x * f
        if bool(_arb_W(xl, eps) < arb(w)):
            return float(_arb_M(xl, eps).lower())
        f -= 1e-6 * (80 - _ + 1) / 80 + 1e-9
    raise ArithmeticError('no certain lower bound of M')


def bounds(L, box, eps):
    """The explicit constants of Theorem E3 for the layer L (with box = layer_box(L)) and the unit of eps >= 0: the crossing level w_# = a_min v_lo/2, the times
    T_2, T_3 (so that the relaxation set D-hat is reached at most T_2 + T_3 after w_e = w_#), the decay rate rho and the prefactor of the Lojasiewicz estimate.
    Upper bounds of times and lower bounds of rates are computed in ball arithmetic from the midpoints of box_t1's constants."""
    amin = L['amin']
    G2 = lambda d: float(np.exp(-d * d / 4.0))
    r, vlo, Ev, vbar, R = box['r'], box['vlo'], box['Ev'], box['vbar'], L['R']
    dt, dh, nu, Hup = box['dtilde'], box['dhat'], box['nu'], box['Hup']
    w_hash = 0.5 * amin * vlo
    out = dict(w_hash=w_hash)
    kd = arb(MC) * arb(w_hash) * arb(amin) * (arb(-r * r / 4)).exp() / 2
    u = max(L['d0'], dt)
    t2 = max(0.0, float(arb(np.log((u - dt) / (dh - dt))))) / float(kd.lower()) if u > dt else 0.0
    Vt = G2(dh) - Ev
    v_up = R + Ev
    mV = M_lower(w_hash, eps)
    mVp = M_lower(amin * v_up, eps)
    t3a = np.log(max((amin * Vt - w_hash) / (amin * VSIG), 1.0)) / mV if amin * Vt > w_hash else 0.0
    t3b = np.log(max((amin * (vbar - v_up)) / (amin * VSIG), 1.0)) / mVp
    t3 = max(t3a, t3b)
    m_lo = min(M_lower(amin * box['vD_lo'], eps), MC * amin ** 2)
    rho = m_lo * nu
    dstar = ((np.exp(-dt * dt / 4) - Ev), (R + Ev))
    out.update(T2=float(t2), T3=float(t3), That=float(t2 + t3), rho=float(rho), pref=float(np.sqrt(Hup / nu)), m_lo=float(m_lo), mV=float(mV), mVp=float(mVp),
               kd=float(kd.lower()), dstar_lo=float(dstar[0]), dstar_hi=float(dstar[1]))
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# the integrated layer flows
# ----------------------------------------------------------------------------------------------------------------------------------
def _rhs(L, eps):
    cs, As = L['cs'], L['As']
    G, amax = LAW['G'], LAW['a_max']
    n = len(L['wb']) + 1
    Gam, dGam = T.Gam, T.dGam

    def f(t, z):
        x, c = z[:n], z[n:]
        w = W_of_x(x, eps)
        gg = Gam(c[:, None] - cs[None, :]) @ As - Gam(c[:, None] - c[None, :]) @ w
        M = dGam(c[:, None] - c[None, :])
        np.fill_diagonal(M, 0.0)
        h = dGam(c[:, None] - cs[None, :]) @ As - M @ w
        return np.r_[(G / amax) * (eps + (2 * VS - eps) * np.tanh(x) ** 2) * gg, MC * w * h]
    return f


def initial_state(L, eps, v0=None):
    """Block at the previous layer's critical point, the entrant at its kernel: x_e = 0 (eps > 0: the dormant state) or W_0(x_e) = v0 a_min sgn (eps = 0)."""
    xb = [x_of_w(w, eps) for w in L['wb']]
    if eps > 0:
        xe = 0.0
    else:
        xe = x_of_w(L['sgn'][-1] * v0 * L['amin'], 0.0)
    return np.r_[xb, xe, L['cb'], L['ce']]


def integrate_layer(L, eps, v0=None, w_hash=None, s2=S2, dt=DT, rtol=1e-11):
    """Stage 1 (Radau, rtol 1e-11): from the start to the first crossing of |w_e| = w_# (terminal event); stage 2: s2 further.  Returns the time-aligned samples s = t - t_# (the stage-1
    samples are log-spaced before t_#), the states x (units x samples), c, and t_#."""
    n = len(L['wb']) + 1
    f = _rhs(L, eps)
    z0 = initial_state(L, eps, v0)
    atol = np.r_[np.full(n, 1e-15), np.full(n, 1e-14)]

    def lev(t, z):
        return abs(float(W_of_x(z[n - 1], eps))) - w_hash
    lev.terminal = True
    lev.direction = 1
    s1 = solve_ivp(f, (0.0, 1e13), z0, method='Radau', rtol=rtol, atol=atol, events=lev, dense_output=True)
    if s1.status != 1:
        raise RuntimeError('the entrant never reached the level')
    th = float(s1.t_events[0][0])
    zh = s1.y_events[0][0]
    pre = th - np.geomspace(1e-3, th, 4000)[::-1] if th > 2e-3 else np.array([0.0])
    pre = np.unique(np.r_[0.0, pre[pre < th]])
    ts2 = th + np.arange(0.0, s2 + 1e-9, dt)
    s2s = solve_ivp(f, (th, th + s2), zh, method='Radau', rtol=rtol, atol=atol, t_eval=ts2)
    Z1 = s1.sol(pre)
    Z = np.concatenate([Z1, s2s.y], axis=1)
    s = np.concatenate([pre - th, s2s.t - th])
    return dict(s=s, x=Z[:n], c=Z[n:], t_hash=th, n=n)


def derived(L, tr, eps, box):
    """Normalised coordinates along a trajectory: v'_i = sgn_i w_i / a_min, d_i = c_i - c*_f(i)."""
    n = tr['n']
    w = W_of_x(tr['x'], eps)
    v = (L['sgn'][:, None] * w) / L['amin']
    d = tr['c'] - L['cs'][L['feats']][:, None]
    return v, d


def check_layer(L, box, bd, tr, eps):
    """The checks of Theorems E1-E4 on one integrated flow (the layer flow of eps from the dormant state, or of eps = 0 from v0)."""
    v, d = derived(L, tr, eps, box)
    s = tr['s']
    n = tr['n']
    r, vlo, vbar, Ev, R = box['r'], box['vlo'], box['vbar'], box['Ev'], L['R']
    dt, dh, vsig = box['dtilde'], box['dhat'], box['varsigma']
    ve = v[-1]
    ok = {}
    # E1: box, monotone entrant, lower faces
    ok['box_d'] = float(np.abs(d).max()) < r
    ok['box_vhi'] = float(v.max()) <= vbar * (1 + 1e-9)               # v_bar is R(1+kappa)/(1-kappa) (1 + 0^+): the entrant converges from below
    ok['box_vlo_block'] = bool(n == 1 or v[:-1].min() >= vlo)
    i_lo = int(np.argmax(ve >= vlo)) if (ve >= vlo).any() else len(ve)
    ok['mono_entrant'] = bool(np.all(np.diff(ve[:i_lo + 1]) >= -1e-13))
    ok['entrant_ge_vlo_after'] = bool(i_lo == len(ve) or ve[i_lo:].min() >= vlo - 1e-12)
    # E3: hitting time of the relaxation set (every unit; the block is in it from the start)
    in_hat = (np.abs(d) <= dh) & (v >= np.exp(-dh * dh / 4) - Ev - vsig) & (v <= R + Ev + vsig)
    allin = in_hat.all(axis=0)
    bad = np.flatnonzero(~allin)
    i_hat = int(bad[-1]) + 1 if len(bad) else 0
    t_hat = float(s[i_hat]) if i_hat < len(s) else float('inf')
    ok['block_in_hat'] = bool(in_hat[:-1].all()) if n > 1 else True
    # E0: the tail of the entrant before the crossing, w_e(s)^(-1/3) >= w_#^(-1/3) + (c_M w_#/3) |s|
    w_e = np.abs(W_of_x(tr['x'][-1], eps))
    wh = bd['w_hash']
    pre_ = (s <= 0) & (w_e > 0)
    ok['tail'] = bool(np.all(w_e[pre_] ** (-1.0 / 3.0) >= (wh ** (-1.0 / 3.0) + c_tail(wh) * wh / 3.0 * np.abs(s[pre_])) * (1 - 1e-6)))
    # E4: Lojasiewicz decay after t_hat
    yh = np.array(L['what'])
    vh = np.abs(yh) / L['amin']
    dhat_c = L['chat'] - L['cs'][L['feats']]
    err = np.sqrt(((v - vh[:, None]) ** 2).sum(axis=0) + ((d - dhat_c[:, None]) ** 2).sum(axis=0))
    if i_hat < len(s):
        env = bd['pref'] * err[i_hat] * np.exp(-bd['rho'] * (s[i_hat:] - s[i_hat]))
        viol = err[i_hat:] - env - 1e-9
        ok['decay'] = bool(viol.max() <= 0)
        ratio = float(np.max(err[i_hat:] / np.maximum(env, 1e-9)))
    else:
        ok['decay'] = False
        ratio = float('nan')
    # time from which every unit stays in D*
    in_star = (np.abs(d) <= dt) & (v >= np.exp(-dt * dt / 4) - Ev) & (v <= R + Ev)
    bad2 = np.flatnonzero(~in_star.all(axis=0))
    i_star = int(bad2[-1]) + 1 if len(bad2) else 0
    t_star = float(s[i_star]) if i_star < len(s) else float('inf')
    # the a priori bound of Theorem E3(iii): T_2 + T_3 plus the Lojasiewicz time from the diameter of D-hat to the ball of radius `marg` (the distance of the end point to the boundary of D*)
    marg = float(min(np.min(dt - np.abs(dhat_c)), np.min(vh - (np.exp(-dt * dt / 4) - Ev)), np.min(R + Ev - vh)))
    wid = (R + Ev + vsig) - (np.exp(-dh * dh / 4) - Ev - vsig)
    diam = float(np.sqrt(n * (wid ** 2 + (2 * dh) ** 2)))
    t_star_bd = bd['That'] + max(0.0, np.log(bd['pref'] * diam / marg)) / bd['rho'] if marg > 0 else float('inf')
    vec = dict(t_hash=tr['t_hash'], t_hat=t_hat, That=bd['That'], ratio_hat=t_hat / bd['That'] if bd['That'] > 0 else float('nan'), ratio_decay=ratio,
               t_star=t_star, t_star_bd=float(t_star_bd), err_end=float(err[-1]), err_hat=float(err[i_hat]) if i_hat < len(s) else float('nan'),
               dmax_over_r=float(np.abs(d).max() / r), vmax_over_vbar=float(v.max() / vbar), margin_star=marg)
    return ok, vec


# ----------------------------------------------------------------------------------------------------------------------------------
# aligned distances between flows
# ----------------------------------------------------------------------------------------------------------------------------------
def aligned_distance(trA, trB, LA, epsA, epsB, s_lo, s_hi, n=4001):
    """sup over s in [s_lo, s_hi] of the distance (v', d) between two time-aligned flows of the same layer (linear interpolation of the samples)."""
    grid = np.concatenate([np.linspace(s_lo, 0.0, n // 2), np.linspace(0.0, s_hi, n // 2)[1:]])
    n_u = trA['n']
    VA, DA = derived(LA, trA, epsA, None)
    VB, DB = derived(LA, trB, epsB, None)
    A = np.vstack([np.vstack([np.interp(grid, trA['s'], VA[i]) for i in range(n_u)]), np.vstack([np.interp(grid, trA['s'], DA[i]) for i in range(n_u)])])
    B = np.vstack([np.vstack([np.interp(grid, trB['s'], VB[i]) for i in range(n_u)]), np.vstack([np.interp(grid, trB['s'], DB[i]) for i in range(n_u)])])
    return float(np.sqrt(((A - B) ** 2).sum(axis=0)).max())


def orbit_distance(trA, trB, LA, epsA, epsB, v_start, frac=0.9, n=4001):
    """Distance of two flows as ORBITS (invariant under a time shift): both are graphs over the entrant's coefficient v'_e (strictly increasing until it reaches its end value);
    sup over the levels v'_e in [v_start, frac * v'_hat_e] (geometric grid) of the distance (block v', d, entrant d) of the points at the same level."""
    n_u = trA['n']
    VA, DA = derived(LA, trA, epsA, None)
    VB, DB = derived(LA, trB, epsB, None)
    top = frac * float(np.abs(LA['what'][-1]) / LA['amin'])
    lv = np.geomspace(max(v_start, 1e-300), top, n)
    def at(V, D_):
        ve = np.maximum.accumulate(V[-1])
        keep = np.r_[True, np.diff(ve) > 0]
        res = []
        for i in range(n_u):
            if i < n_u - 1:
                res.append(np.interp(lv, ve[keep], V[i][keep]))
            res.append(np.interp(lv, ve[keep], D_[i][keep]))
        return np.vstack(res)
    return float(np.sqrt(((at(VA, DA) - at(VB, DB)) ** 2).sum(axis=0)).max())


# ----------------------------------------------------------------------------------------------------------------------------------
# the tail of the entrant (Lemma E0): M_eps(w) >= c_M w^(4/3) for w <= w_#, uniformly in eps
# ----------------------------------------------------------------------------------------------------------------------------------
def c_tail(w_hash):
    """c_M = (G m / a_max^2) (2 v_s a_max / (m cosh^4 x_#))^(2/3), x_# = x_0(w_#): the explicit uniform constant of the tail lemma (NOTE.md Lemma E0)."""
    xh = x_of_w(w_hash, 0.0)
    return (LAW['G'] * LAW['m'] / LAW['a_max'] ** 2) * (2 * VS * LAW['a_max'] / (LAW['m'] * np.cosh(xh) ** 4)) ** (2.0 / 3.0)
