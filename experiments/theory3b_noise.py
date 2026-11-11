"""
A deterministic, non-asymptotic noise theorem for the limit object (the sliding inverse scale space): entry order, entry times, the stopping window
and the oracle fit with an explicit second-order remainder, and its high-probability form
====================================================================================================================================

Setting (the continuum separated regime of Theorem S-K; every length in units of the kernel width sigma, every coefficient in units of the smallest
feature amplitude A_min, the inner product normalised so that <phi_c phi_c'> = Gamma(c - c'), Gamma(d) = exp(-d^2/4), equal walls rho = 1, time in
units of 1/(Z A_min), Z = <phi^2>): K Gaussian features of the kernels' width, separation >= Delta, amplitude ratio R, a finite dictionary on a lattice.
The data are y = f* + xi, xi ANY element of the inner product space (no statistics) with the dual-norm gauges, over the box radius r about the K features'
centres N_r and the dictionary D,
        e0 = max( sup_{c in N_r u D} |<xi phi_c>| ,  e1 = sigma sup_{N_r} |<xi d_c phi_c>| ,  e2 = sigma^2 sup_{N_r} |<xi d_c^2 phi_c>| )  / (Z A_min).
S(e0, e1, e2; r) is a FINITE explicit system of inequalities on the noiseless path's data (the stages' blocks, rates, duals, hitting times) and on the box
constants of Proposition S-A with the noise added to its couplings. If it holds, the noisy path (with the box relaxation of Prop. S-A, which converges for
the noisy fit as well) has (i) the same K entries in the same order and no other entry, entry times within eta_k of the noiseless ones, blocks within
Theta_k of the noiseless blocks; (ii) after the K-th entry a window of length >= W_lo in which no other unit enters and the block is the K-unit least-squares
fit; (iii) in it ||f_hat - f* - P_T xi|| <= Rem, P_T the projection onto the tangent space of the K-unit model at the truth (the first-order oracle). Every bound is
monotone increasing in (e0, e1, e2). The proof is in the author's working note; this file evaluates the constants (python-flint ball arithmetic for the inequalities, the
limit-object data of the noiseless path in double precision, widened by 1e-9), validates the claims on the repository's path computation (`baselines.iss_path`) in
a wide window whose sampled Gram is the continuum one, and gives the Gaussian tail version (Lemma G).

Usage:  python experiments/theory3b_noise.py --run [--parts V,G,W] | --from-raw      -> results/synthetic/theory3b.md (+ _raw.npz)
Environment: T3B_OUT (output directory).
"""
import argparse
import math
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.optimize import root

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

try:
    from flint import arb
    HAVE_ARB = True
except Exception:                                    # python-flint missing: floats only (the report says which)
    HAVE_ARB = False

OUT = os.environ.get('T3B_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory3b_raw.npz'), os.path.join(OUT, 'theory3b.md')

SIG = 0.12                                           # the repository's kernel width
HALF, HX = 3.0, 1.0 / 30.0                           # wide window [-3, 3], sample spacing 1/30: the sampled Gram of every kernel with centre in [-1.6, 1.6] is the continuum one
DICT_HALF, DICT_H = 1.6, 0.1                         # dictionary: lattice of spacing 0.1 on [-1.6, 1.6] (33 kernels), the repository's lattice
WIDEN = 1e-9                                         # relative widening of the limit-object data (double precision) in the ball-arithmetic run
S2, S6 = math.sqrt(2.0), math.sqrt(6.0)
M_J = (1.0, 0.5, 0.75, 15.0 / 8.0, 105.0 / 16.0)     # m_j = <(sigma^j d^j phi)^2>/<phi^2>: 1, 1/2, 3/4, 15/8, 105/16


# ------------------------------------------------------------------------------------------------------------------------------------
# arithmetic backends: floats, or python-flint balls (a comparison a > b is then True only if it is certain)
# ------------------------------------------------------------------------------------------------------------------------------------
class Fl:
    name = 'float'

    @staticmethod
    def c(x, rad=0.0):
        return float(x)

    exp = staticmethod(math.exp)
    sqrt = staticmethod(math.sqrt)

    @staticmethod
    def mx(a, b):
        return a if a >= b else b

    @staticmethod
    def up(a):
        return float(a)

    @staticmethod
    def lo(a):
        return float(a)


class Ar:
    name = 'arb'

    @staticmethod
    def c(x, rad=0.0):
        return arb(float(x), float(rad)) if rad else arb(float(x))

    @staticmethod
    def exp(a):
        return a.exp()

    @staticmethod
    def sqrt(a):
        return a.sqrt()

    @staticmethod
    def mx(a, b):
        return a.max(b)

    @staticmethod
    def up(a):
        return float(a.upper())

    @staticmethod
    def lo(a):
        return float(a.lower())


def gam(x):
    return np.exp(-np.asarray(x, float) ** 2 / 4.0)


def psi(x):
    return 0.5 * np.asarray(x, float) * np.exp(-np.asarray(x, float) ** 2 / 4.0)


# ------------------------------------------------------------------------------------------------------------------------------------
# the box of Proposition S-A with the noise added to its couplings
# ------------------------------------------------------------------------------------------------------------------------------------
def box_noise(B, K, R, Delta, r, e0, e1, e2):
    """Proposition S-A's conditions (L), (C), (N) when the residual carries the noise: the correlation <xi phi_i> (<= e0) adds to the coefficient coupling P_i, <xi d_c phi_i>
    (<= e1) to the centre coupling C_i, and the Hessian of the fit gains, in each unit's own 2x2 block, the entries -<xi d_c phi_i> (<= e1) off the diagonal and
    -w_i <xi d_c^2 phi_i> (<= vbar e2) on the centre's diagonal. Arguments are backend numbers. Returns dict(ok, lam = the Hessian floor on the convex set D*, dstar, vlo, vbar, ...)."""
    out = dict(ok=False)
    c = B.c
    G = lambda x: B.exp(-(x * x) / 4)
    ps = lambda x: x / 2 * G(x)
    G2 = lambda x: (x * x / 4 - c(0.5)) * G(x)
    rr, Rr, Kk = r, R, c(K - 1)
    if not (rr < c(S2)):
        return out
    ell = c(Delta) - 2 * rr
    if not (ell >= c(S6)):
        return out
    kap = Kk * G(ell)
    if not (kap < c(1)):
        return out
    vbar = (Rr * (1 + kap) + e0) / (1 - kap) * c(1 + 1e-9)
    Ev = kap * (Rr + vbar) + e0
    Ec = Kk * (Rr + vbar) * ps(ell) + e1
    vlo = G(rr) - Ev
    Cm = ps(rr) - Ec
    out.update(vbar=vbar, Ev=Ev, Ec=Ec, vlo=vlo, Cm=Cm)
    if not (vlo > c(0) and Cm > c(0)):
        return out
    dstar = 2 * Ec / G(rr)
    if not (dstar < rr):
        return out
    vlod = G(dstar) - Ev
    if not (vlod > c(0)):
        return out
    s_max = Rr * ps(dstar) + e1
    t_min = vlod * (c(0.5) - dstar * dstar / 4) * G(dstar) - vbar * e2
    lam_own = (1 + t_min - B.sqrt((1 - t_min) * (1 - t_min) + 4 * s_max * s_max)) / 2
    rowv = Kk * (G(ell) + (Rr + vbar) * ps(ell) + vbar * ps(ell))
    rowd = Kk * ((Rr + vbar) * ps(ell) + vbar * (Rr + vbar) * G2(ell) + vbar * ps(ell) + vbar * vbar * G2(ell))
    Hc = B.mx(rowv, rowd)
    lam = lam_own - Hc
    out.update(dstar=dstar, vlod=vlod, lam_own=lam_own, Hc=Hc, lam=lam, ok=bool(lam > c(0)))
    return out


# ------------------------------------------------------------------------------------------------------------------------------------
# the noiseless path of the limit object in the continuum (normalised units)
# ------------------------------------------------------------------------------------------------------------------------------------
def _crit(cs, a, feats, c_init, v_init):
    """Critical point of the fit of the units started at features `feats` (indices) to sum_l a_l phi(cs_l): unknowns (v_i, c_i); the stationarity conditions
    <r phi_i> = 0, <r d_c phi_i> = 0 (v_i != 0). Newton by scipy.optimize.root with the exact Jacobian, to ~1e-15."""
    k = len(feats)

    def F(x):
        v, c = x[:k], x[k:]
        g = np.zeros(k)
        h = np.zeros(k)
        for i in range(k):
            d = c[i] - cs                               # c_i - c*_l
            g[i] = (a * gam(d)).sum() - (v * gam(c[i] - c)).sum()
            # <r d_c phi_i> = sum_l a_l d_{c_i} Gamma(c_i - c*_l) - sum_m v_m d_{c_i} Gamma(c_i - c_m),  d Gamma(x)/dx = -psi(x) (psi odd continuation)
            h[i] = -(a * psi(d)).sum() + (v * psi(c[i] - c)).sum()
        return np.r_[g, h]
    s = root(F, np.r_[v_init, c_init], method='hybr', tol=1e-15)
    return s.x[:k], s.x[k:], bool(np.abs(F(s.x)).max() < 1e-13)


def noiseless_path(cs_sig, As, c0_sig):
    """The limit object's path on noiseless data in the continuum (normalised units: lengths in sigma, A_min = 1, Z = 1, walls 1). cs_sig: the features' centres / sigma (sorted),
    As the amplitudes (any scale), c0_sig the dictionary centres / sigma. Returns dict with per-stage data (stage k = after k entries, k = 0..K):
      a (signed amplitudes / A_min), R, Delta (min separation), K, c0, cs
      act[k] (dictionary indices), feat[k] (assigned features), v[k], c[k] (block coefficients, centres), p[k] (rates of all dictionary kernels; entries of active kernels unused),
      u[k] (duals at tau_k), tau[k], ent[k] = entrant of stage k (k < K), T[k] = tau_{k+1} - tau_k, d0[k] (entrant's offset from its feature, widths), mu[k] (runner-up margin),
      gap[k] (relative hitting-time gap), ok (signs agree, Newton converges, unique entrant)."""
    cs, As, c0 = np.asarray(cs_sig, float), np.asarray(As, float), np.asarray(c0_sig, float)
    K, d = len(cs), len(c0)
    amin = np.abs(As).min()
    a = As / amin
    Delta = float(np.diff(np.sort(cs)).min()) if K > 1 else 1e3                   # K = 1: no other feature (any large separation)
    R = float(np.abs(a).max())
    P = dict(a=a, R=R, Delta=Delta, K=K, c0=c0, cs=cs, act=[[]], feat=[[]], v=[np.zeros(0)], c=[np.zeros(0)], p=[], u=[np.zeros(d)], tau=[0.0], ent=[], T=[], d0=[], mu=[], gap=[], ok=True, why='')
    u = np.zeros(d)
    tau = 0.0
    act, feat = [], []
    v, c = np.zeros(0), np.zeros(0)
    for k in range(K + 1):
        p = (a[None, :] * gam(c0[:, None] - cs[None, :])).sum(1)
        if len(v):
            p = p - (v[None, :] * gam(c0[:, None] - c[None, :])).sum(1)
        P['p'].append(p)
        if k == K:
            break
        ina = np.array([j not in act for j in range(d)])
        dt = np.full(d, np.inf)
        ok = ina & (np.abs(p) > 1e-14)
        dt[ok] = (1 - np.sign(p[ok]) * u[ok]) / np.abs(p[ok])
        o = np.argsort(dt)
        j = int(o[0])
        T = float(dt[j])
        P['gap'].append(float((dt[o[1]] - dt[o[0]]) / dt[o[0]]))
        u = np.where(ina, u + T * p, u)
        u[j] = np.sign(p[j])
        oth = [q for q in range(d) if ina[q] and q != j]
        P['mu'].append(float(1 - np.abs(u[oth]).max()))
        free = [i for i in range(K) if i not in feat]
        f_as = free[int(np.argmin(np.abs(cs[free] - c0[j])))]
        P['d0'].append(float(abs(c0[j] - cs[f_as])))
        P['ent'].append(j)
        P['T'].append(T)
        if np.sign(p[j]) != np.sign(a[f_as]):
            P['ok'], P['why'] = False, 'sign'
        tau += T
        act.append(j)
        feat.append(f_as)
        v, c, okc = _crit(cs, a, feat, cs[feat], a[feat])
        if not okc:
            P['ok'], P['why'] = False, 'newton'
        P['act'].append(list(act))
        P['feat'].append(list(feat))
        P['v'].append(v.copy())
        P['c'].append(c.copy())
        P['u'].append(u.copy())
        P['tau'].append(tau)
    return P


# ------------------------------------------------------------------------------------------------------------------------------------
# the Gram of the K-unit model at the truth and the constants of (iii)
# ------------------------------------------------------------------------------------------------------------------------------------
def gram_star(cs, a):
    """H* = J*^T J* of the model sum_k v_k phi(c*_k + d_k), parameters (v_k, d_k): <phi_k phi_l> = Gamma, <phi_k d phi_l> = ((c_k - c_l)/2) Gamma, <d phi_k d phi_l> = (1/2 - (c_k - c_l)^2/4) Gamma;
    the d-columns carry the factor a_k."""
    K = len(cs)
    H = np.zeros((2 * K, 2 * K))
    for k in range(K):
        for l in range(K):
            x = cs[k] - cs[l]
            g = float(gam(x))
            H[k, l] = g                                    # <phi_k phi_l>
            H[k, K + l] = a[l] * (x / 2) * g               # <phi_k, a_l d_c phi_l>
            H[K + k, l] = a[k] * (-x / 2) * g              # <a_k d phi_k, phi_l> = a_k d/dc_k Gamma(c_k - c_l) = -a_k (x/2) g
            H[K + k, K + l] = a[k] * a[l] * (0.5 - x * x / 4) * g
    return 0.5 * (H + H.T)


# ------------------------------------------------------------------------------------------------------------------------------------
# the finite system S(e; r)
# ------------------------------------------------------------------------------------------------------------------------------------
def _cap(B, x, y):
    """Gamma-bar and psi-bar of |x| +- y: max of Gamma and psi over [max(0, |x| - y), |x| + y] (backend numbers), psi's maximum sqrt(2) included if the interval can contain it."""
    c = B.c
    lo = B.mx(c(0), x - y)
    hi = x + y
    G = lambda t: B.exp(-(t * t) / 4)
    ps = lambda t: t / 2 * G(t)
    gb = G(lo)
    pb = B.mx(ps(lo), ps(hi))
    if B.lo(hi) >= S2 and B.up(lo) <= S2:
        pb = B.mx(pb, c(math.exp(-0.5) / S2 * (1 + 1e-12)))
    return gb, pb


def system(P, e0, e1, e2, r, B=Fl, ez=None):
    """Evaluate S(e0, e1, e2; r) on the noiseless data P (from `noiseless_path`) with the backend B. Returns dict with ok (bool) and, when computed, eta[k] (entry-time shifts),
    Theta[k], W_lo, rem (the remainder of (iii)), lam (Hessian floor), the first failing condition `why`, and detail."""
    c = B.c
    K, R, Delta = P['K'], P['R'], P['Delta']
    rad = lambda x: WIDEN * max(1.0, abs(float(x))) if B is Ar else 0.0
    cc = lambda x: c(x, rad(x))
    E0, E1, E2 = c(e0), c(e1), c(e2)
    out = dict(ok=False, why='')
    if not P['ok']:
        out['why'] = 'noiseless path not valid'
        return out
    rr = c(r)
    if not (min(P['d0']) >= 0 and rr > c(max(P['d0']))):
        out['why'] = 'r <= d0'
        return out
    bx = box_noise(B, K, c(R), Delta, rr, E0, E1, E2)
    if not bx['ok']:
        out['why'] = 'box (L)(C)(N)'
        out['box'] = bx
        return out
    lam, vbar = bx['lam'], bx['vbar']
    out.update(lam=B.up(lam) if B is Ar else lam, vbar=B.up(vbar) if B is Ar else vbar)
    dstar = bx['dstar']
    # the noiseless blocks lie in D*: |d_i| <= dstar, |v_i| in [Gamma(dstar) - Ev, R + Ev]
    for k in range(1, K + 1):
        for i, f in enumerate(P['feat'][k]):
            dd = abs(P['c'][k][i] - P['cs'][f])
            vv = abs(P['v'][k][i])
            if not (cc(dd) < dstar and cc(vv) > bx['vlod'] and cc(vv) < c(R) + bx['Ev']):
                out['why'] = 'noiseless block outside D*'
                return out
    d = len(P['c0'])
    c0 = P['c0']
    # Theta_k and rate errors eps_j^{(k)}
    Theta, eps = [], []
    for k in range(K + 1):
        vk = np.abs(P['v'][k])
        if k == 0:
            Theta.append(c(0))
            eps.append([E0] * d)
            continue
        th = B.sqrt(sum((E0 * E0 + cc(vi) * cc(vi) * E1 * E1 for vi in vk), c(0))) / lam
        Theta.append(th)
        ek = []
        for j in range(d):
            s = c(0)
            for i, ci in enumerate(P['c'][k]):
                gb, pb = _cap(B, cc(abs(c0[j] - ci)), th)
                s = s + gb * gb + cc(vk[i]) * cc(vk[i]) * pb * pb
            ek.append(E0 + th * B.sqrt(s))
        eps.append(ek)
    out['Theta'] = [B.up(t) for t in Theta]
    # recursion over the stages: eta_k (entry-time shift), ebar_{k,j} (dual state error at the entry), the competitors stay inside the walls
    eta = [c(0)]
    ebar = [c(0)] * d
    etas = [0.0]
    for k in range(K):
        j0 = P['ent'][k]
        act_next = set(P['act'][k + 1])
        pk = P['p'][k]
        uk, uk1 = P['u'][k], P['u'][k + 1]
        Tk = cc(P['T'][k])
        ps0 = cc(abs(pk[j0]))
        ej = eps[k]
        if not (ps0 > ej[j0]):
            out['why'] = 'entrant rate <= rate error (stage %d)' % k
            return out
        e_new = (ebar[j0] + ej[j0] * Tk + (ps0 + ej[j0]) * eta[k]) / (ps0 - ej[j0])
        # competitors j in the inactive set (not j0, not already active)
        for j in range(d):
            if j == j0 or j in P['act'][k]:
                continue
            pj = cc(abs(pk[j]))
            lo_ends = B.mx(cc(abs(uk[j])) + pj * eta[k], cc(abs(uk1[j])) + pj * e_new)
            span = B.mx(eta[k], Tk + e_new)
            Mk = lo_ends + ebar[j] + ej[j] * span + (pj + ej[j]) * eta[k]
            if not (Mk < c(1)):
                out['why'] = 'order: kernel %d not inside its wall at stage %d' % (j, k)
                out['eta'] = etas
                return out
        eb_new = list(ebar)
        for j in range(d):
            if j == j0 or j in P['act'][k]:
                continue
            pj = cc(abs(pk[j]))
            eb_new[j] = ebar[j] + ej[j] * Tk + (pj + ej[j]) * (eta[k] + e_new)
        eta.append(e_new)
        etas.append(B.up(e_new))
        ebar = eb_new
    out['eta'] = etas
    # (iii): the remainder of the fit error against the first-order oracle P_T xi
    a = P['a']
    H = gram_star(P['cs'], a)
    ev = np.linalg.eigvalsh(H)
    lmin, lmax = float(ev[0]), float(ev[-1])
    if lmin <= 0:
        out['why'] = 'Gram of the tangent space singular'
        return out
    lmin_c, lmax_c = c(lmin * (1 - 1e-9)), c(lmax * (1 + 1e-9))
    ThK = Theta[K]
    M2 = c(1 / S2) + c(math.sqrt(3) / 2) * vbar
    LJ = B.sqrt(c(0.5) + (c(1 / S2) + c(math.sqrt(3) / 2) * c(R)) * (c(1 / S2) + c(math.sqrt(3) / 2) * c(R)))
    q = M2 * ThK * ThK / 2
    rem = q + (1 / B.sqrt(lmin_c)) * ((c(S2) * E1 + c(R) * E2) * ThK + LJ * ThK * (B.sqrt(lmax_c) * ThK + M2 * ThK * ThK / 2))
    out['rem'] = B.up(rem)
    out['lmin'], out['lmax'] = lmin, lmax
    # the window after the K-th entry: crude (every rate error by the triangle inequality) and, if the orthogonalised gauge ez is given, refined
    # (the inactive rates are <xi, (I - P_T) phi_j> - <rho, phi_j>, |<rho, phi_j>| <= rem)
    pK = P['p'][K]
    uK = P['u'][K]
    Wlo = Wz = None
    ezv = None if ez is None else (np.full(d, float(ez)) if np.ndim(ez) == 0 else np.asarray(ez, float))
    for j in range(d):
        if j in P['act'][K]:
            continue
        num = c(1) - cc(abs(uK[j])) - ebar[j]
        if not (num > c(0)):
            out['why'] = 'window: kernel %d has no room after the K-th entry' % j
            return out
        Wj = num / (eps[K][j] + cc(abs(pK[j])))
        Wlo = Wj if Wlo is None else (Wj if Wj < Wlo else Wlo)
        if ezv is not None:
            Wjz = num / (c(ezv[j]) + rem + cc(abs(pK[j])))
            Wz = Wjz if Wz is None else (Wjz if Wjz < Wz else Wz)
    out['W_lo'] = B.lo(Wlo) if Wlo is not None else float('inf')
    if Wz is not None:
        out['W_lo_z'] = max(B.lo(Wz), out['W_lo'])                         # both are valid lower bounds on the window
    out['ok'] = True
    return out


def e_crit(P, r, ray=(1.0, 1.0, 1.0), B=Fl, hi=1.0, iters=44):
    """The largest scale t such that S(t * ray; r) holds, by bisection (S is monotone in t: every bound increases with the gauges); 0 if S fails at t = 1e-9."""
    f = lambda t: system(P, t * ray[0], t * ray[1], t * ray[2], r, B)['ok']
    if not f(1e-9):
        return 0.0
    lo, h = 1e-9, hi
    while f(h) and h < 1e3:
        lo, h = h, h * 2
    for _ in range(iters):
        m = math.sqrt(lo * h)
        lo, h = (m, h) if f(m) else (lo, m)
    return lo


# ------------------------------------------------------------------------------------------------------------------------------------
# the repository's window and path, and the gauges of an actual noise vector
# ------------------------------------------------------------------------------------------------------------------------------------
_ST = {}


def setup():
    """The wide window [-3, 3] at the repository's sample spacing, the kernel object, the dictionary (lattice 0.1 on [-1.6, 1.6]), Z = <phi^2> (sampled, equal to the continuum value to 1e-12)."""
    if not _ST:
        from src.core import rich_transit_nd as RN
        xs = np.round(np.arange(-HALF, HALF + 1e-9, HX), 12)
        Kobj = RN.BumpsND(xs[:, None], SIG)
        c0 = np.round(np.arange(-DICT_HALF, DICT_HALF + 1e-9, DICT_H), 10)
        Z = float(np.mean(Kobj.phi(np.array([[0.123]]))[:, 0] ** 2))
        _ST.update(K=Kobj, c0=c0, Z=Z, xs=xs, N=xs.size)
    return _ST['K'], _ST['c0'], _ST['Z']


def _he(j, t):
    return (np.ones_like(t), t, t * t - 1, t ** 3 - 3 * t, t ** 4 - 6 * t * t + 3)[j]


def gauges(xi, cs_abs, c0_abs, r, amin, mesh=0.01):
    """RIGOROUS upper bounds on the gauges of the noise vector xi (samples) over the box N_r = U_k [c*_k - r sigma, c*_k + r sigma] and the dictionary: e0, e1, e2 =
    (sup |<xi sigma^j d_c^j phi_c>|)/(Z A_min), j = 0, 1, 2 (e0 also over the dictionary kernels). The suprema over the continuum of centres: the maximum on a grid of spacing `mesh` (in
    widths) plus (mesh/2) times a bound on the next derivative's supremum, the last one (j = 4) by Cauchy-Schwarz, <(sigma^4 d^4 phi)^2>^(1/2) ||xi||. Returns (e0, e1, e2)."""
    Kobj, _, Z = setup()
    xs = Kobj.xs[:, 0]
    N = xs.size
    xi = np.asarray(xi, float)
    nrm = float(np.sqrt(np.mean(xi * xi)))
    grid = np.unique(np.concatenate([cs_abs[:, None] + SIG * np.arange(-math.ceil(r / mesh), math.ceil(r / mesh) + 1)[None, :] * mesh]).ravel())
    grid = grid[(np.abs(grid[:, None] - cs_abs[None, :]) <= r * SIG + 1e-12).any(1)]
    t = (xs[:, None] - grid[None, :]) / SIG
    ph = np.exp(-t * t / 2)
    gj = [np.abs(((xi[:, None] * _he(j, t) * ph).sum(0) / N)) / (Z * amin) for j in range(4)]
    m = [g.max() for g in gj]
    M4 = nrm * math.sqrt(Z * M_J[4]) / (Z * amin)
    M3 = m[3] + mesh / 2 * M4
    M2 = m[2] + mesh / 2 * M3
    M1 = m[1] + mesh / 2 * M2
    M0 = m[0] + mesh / 2 * M1
    t0 = (xs[:, None] - c0_abs[None, :]) / SIG
    dct = float(np.abs((xi[:, None] * np.exp(-t0 * t0 / 2)).sum(0) / N).max() / (Z * amin))
    return max(M0, dct), M1, M2


def run_path(y, K_features, max_extra=2):
    """The repository's sliding inverse scale space (baselines.iss_path) on the data y in the wide window, K_features + max_extra entries. Returns dict(order, tau (normalised: raw tau Z A_min given later),
    blocks)."""
    from src.core import baselines as BL
    Kobj, c0, _ = setup()
    R = BL.iss_path(Kobj, y, c0, max_entries=K_features + max_extra)
    return R


def target_f(cs_abs, As):
    Kobj, _, _ = setup()
    return Kobj.synth(np.asarray(cs_abs, float)[:, None], np.asarray(As, float))


def make_target(K, sep, R, mode, seed):
    """K features in the window: 'worst': every feature about 0.375 widths (0.045 +- 10 %) from a lattice point, the gaps m h (m = round(sep sigma / h)) plus the offset differences (distinct offsets: no
    exact tie of two features' hitting times); 'random': the first feature at a uniform position, the others at exactly `sep` widths. Amplitudes: magnitudes (1, R) or (1, sqrt R, R) in a random
    order, scaled 0.7, random signs; the configuration at random shift, at most 1.05 from the origin. Returns (centres (sorted, absolute), amplitudes)."""
    rng = np.random.default_rng(seed)
    h = DICT_H
    if K == 1:
        cs = np.array([rng.uniform(-0.5, 0.5)])
    else:
        if mode == 'worst':
            offs = rng.choice([-1.0, 1.0], K) * 0.045 * rng.uniform(0.9, 1.1, K)
            m = int(round(sep * SIG / h))
            span = (K - 1) * m * h
        else:
            span = (K - 1) * sep * SIG
        lim = 1.05 - 0.1
        shift_max = max(0.0, min(0.3, lim - span / 2))
        shift = rng.uniform(-shift_max, shift_max)
        if mode == 'worst':
            cs = np.round((np.arange(K) * m * h - span / 2 + shift) / h) * h + offs
        else:
            cs = np.arange(K) * sep * SIG - span / 2 + shift + rng.uniform(-h / 2, h / 2)
    mags = ([1.0], [1.0, R], [1.0, math.sqrt(R), R])[K - 1]
    mags = np.array(mags)[rng.permutation(K)] * 0.7
    return np.sort(cs), mags * rng.choice([-1.0, 1.0], K)


R_GRID = (0.5, 0.6, 0.8, 1.0)                       # box radii tried (each gives a valid instance of the theorem; the tightest bound of those that hold is used)


def noise_direction(P, cs_abs, fam, rng, amin):
    """A noise vector (raw samples) of a given family, unit-free (scaled later): 'white' i.i.d. N(0,1); 'comp' = the sign-aligned kernel of the stage-0 runner-up
    (the kernel that comes next after the first entrant in the hitting-time order); 'mid' = a kernel between the features (a unit that would fit noise)."""
    Kobj, c0, Z = setup()
    N = Kobj.N
    if fam == 'white':
        return rng.standard_normal(N)
    p0 = P['p'][0]
    if fam == 'comp':
        dt = np.where(np.abs(p0) > 1e-14, (1 - np.sign(p0) * 0) / np.maximum(np.abs(p0), 1e-300), np.inf)
        dt[P['ent'][0]] = np.inf
        j = int(np.argmin(dt))
        return np.sign(p0[j]) * Kobj.phi(np.array([[c0[j]]]))[:, 0]
    if fam == 'mid':
        cm = 0.5 * (cs_abs[0] + cs_abs[-1]) if len(cs_abs) > 1 else cs_abs[0] + 4 * SIG
        return Kobj.phi(np.array([[round(cm / DICT_H) * DICT_H]]))[:, 0] * (1 if rng.random() < 0.5 else -1)
    raise ValueError(fam)


def best_instance(P, cs_abs, xi0, amin, rgrid=R_GRID):
    """For the unit-scale noise vector xi0: for each box radius r the gauges g0(r) (linear in the scale), the largest scale t(r) with S(t g0(r); r) true (bisection, floats) and the best r.
    Returns (t_best, r_best, g0_best, table of (r, t, g0))."""
    c0 = setup()[1]
    tab = []
    for r in rgrid:
        g0 = gauges(xi0, cs_abs, c0, r, amin)
        t = e_crit(P, r, ray=g0, hi=1e-3)
        tab.append((r, t, g0))
    b = max(tab, key=lambda x: x[1])
    return b[1], b[0], b[2], tab


def actual(P, cs_abs, As, f, xi, K, amin):
    """Run the repository's path on f + xi and extract the quantities the theorem bounds (normalised units)."""
    Kobj, c0, Z = setup()
    R = run_path(f + xi, K, 2)
    order = list(R['order'])
    tau = np.array(R['taus']) * Z * amin
    out = dict(order_ok=bool(order[:K] == list(P['ent']) and len(order) >= K), n=len(order))
    if len(tau) >= K and out['order_ok']:
        out['t_err'] = np.abs(tau[:K] - np.array(P['tau'][1:K + 1]))
        out['win'] = float(tau[K] - tau[K - 1]) if len(tau) > K else np.inf
        th = []
        for k in range(1, K + 1):
            b = R['blocks'][k]
            idx = sorted(P['act'][k])
            # the block of the path is sorted by dictionary index
            wv = np.asarray(b['w']) / amin
            cv = np.asarray(b['c']).ravel() / SIG
            ref_v = np.zeros(len(idx))
            ref_c = np.zeros(len(idx))
            for i, jj in enumerate(P['act'][k]):
                q = idx.index(jj)
                ref_v[q], ref_c[q] = P['v'][k][i], P['c'][k][i]
            th.append(float(np.sqrt(((wv - ref_v) ** 2).sum() + ((cv - ref_c) ** 2).sum())))
        out['theta'] = np.array(th)
        # the fit error against the first-order oracle
        b = R['blocks'][K]
        fhat = Kobj.synth(np.asarray(b['c']).reshape(-1, 1), np.asarray(b['w']))
        t_ = (Kobj.xs[:, 0][:, None] - cs_abs[None, :]) / SIG
        ph = np.exp(-t_ * t_ / 2)
        a_ = As / amin
        Jc = amin * np.concatenate([ph, ph * (t_ * a_[None, :])], 1)
        coef = np.linalg.lstsq(Jc, xi, rcond=None)[0]
        PT = Jc @ coef
        out['rem'] = float(np.sqrt(np.mean((fhat - f - PT) ** 2)) / (amin * math.sqrt(Z)))
        out['lin'] = float(np.sqrt(np.mean(PT ** 2)) / (amin * math.sqrt(Z)))
        out['err'] = float(np.sqrt(np.mean((fhat - f) ** 2)) / (amin * math.sqrt(Z)))
    return out


def tangent_cols(cs_abs, As, amin):
    """Raw columns of J* (sample space): d/dv_k = A_min phi(c*_k), d/dd_k = A_min a_k He_1(t) phi (d in widths): shape (N, 2K)."""
    Kobj = setup()[0]
    t_ = (Kobj.xs[:, 0][:, None] - cs_abs[None, :]) / SIG
    ph = np.exp(-t_ * t_ / 2)
    return amin * np.concatenate([ph, ph * (t_ * (As / amin)[None, :])], 1)


def zeta_family(cs_abs, As, amin):
    """zeta_j = (I - P_T) phi_j for every dictionary kernel j (sample space, T = the K-unit model's tangent space at the truth): array (N, d), and the normalised norms ||zeta_j||/sqrt(Z)."""
    Kobj, c0, Z = setup()
    Jc = tangent_cols(cs_abs, As, amin)
    Phi = Kobj.phi(c0[:, None])
    Zeta = Phi - Jc @ np.linalg.lstsq(Jc, Phi, rcond=None)[0]
    return Zeta, np.sqrt(np.mean(Zeta ** 2, 0) / Z)


def zeta_gauge(xi, cs_abs, As, amin):
    """The orthogonalised dictionary gauges e_zeta,j = |<xi, zeta_j>|/(Z A_min), one per dictionary kernel (array of length d)."""
    Kobj, c0, Z = setup()
    Zeta, _ = zeta_family(cs_abs, As, amin)
    return np.abs((np.asarray(xi)[:, None] * Zeta).mean(0)) / (Z * amin)


# ------------------------------------------------------------------------------------------------------------------------------------
# Lemma G: sub-Gaussian noise
# ------------------------------------------------------------------------------------------------------------------------------------
def gauss_gauges(st, N, K, d, r, p, zeta_norms=None):
    """High-probability gauge bounds for i.i.d. sub-Gaussian noise of variance proxy s^2 (e.g. N(0, s^2)), st = s/(A_min sqrt(N Z)) the noise per kernel width: with probability >= 1 - p,
    e0 <= E0, e1 <= E1, e2 <= E2 (the gauges over the box of radius r about the K features and the d dictionary kernels) and e_zeta,j <= Ez[j] (if the norms ||zeta_j||/sqrt(Z) are given).
    Proof (note, Lemma G): each of nev = 3 K (2M+1) + d [+ d] Gaussian events <xi, sigma^j d^j phi_c> at the net points c* + i r/M (|i| <= M) is sub-Gaussian with proxy st^2 m_j (m_j = 1, 1/2, 3/4;
    for zeta_j st^2 ||zeta_j||^2): |g| <= st sqrt(2 m_j ln(4 nev/p)) fails with probability <= p/(2 nev); the continuum between net points by
    sup|g_j| <= max_net|g_j| + (m/2) sup|g_{j+1}| (m = r/M), the chain closed by Cauchy-Schwarz sup|g_3| <= st sqrt(N m_3 (1 + 2 sqrt(tau/N) + 2 tau/N)), tau = ln(2/p),
    which holds with probability >= 1 - p/2 (Hsu-Kakade-Zhang's tail of ||xi||^2 <= s^2 (N + 2 sqrt(N tau) + 2 tau))."""
    M = max(10, int(math.ceil(4.0 * r * math.sqrt(N))))
    m = r / M
    nev = 3 * K * (2 * M + 1) + d + (d if zeta_norms is not None else 0)
    L = math.sqrt(2.0 * math.log(4.0 * nev / p))
    t = [st * math.sqrt(M_J[j]) * L for j in range(3)]
    tau = math.log(2.0 / p)
    M3 = st * math.sqrt(N * M_J[3] * (1 + 2 * math.sqrt(tau / N) + 2 * tau / N))
    E2 = t[2] + m / 2 * M3
    E1 = t[1] + m / 2 * E2
    E0 = t[0] + m / 2 * E1
    Ez = None if zeta_norms is None else st * np.asarray(zeta_norms) * L
    return E0, E1, E2, Ez


def st_crit(P, r, p, N, d, zeta_norms=None, B=Fl, hi=1e-3, iters=40):
    """The largest st with S(E(st, p); r) true (bisection; S is monotone in the gauges and the gauges are increasing in st)."""
    def ok(st):
        E0, E1, E2, Ez = gauss_gauges(st, N, P['K'], d, r, p, zeta_norms)
        return system(P, E0, E1, E2, r, B, ez=Ez)['ok']
    if not ok(1e-9):
        return 0.0
    lo, h = 1e-9, hi
    while ok(h) and h < 1e2:
        lo, h = h, h * 2
    for _ in range(iters):
        mm = math.sqrt(lo * h)
        lo, h = (mm, h) if ok(mm) else (lo, mm)
    return lo


# ------------------------------------------------------------------------------------------------------------------------------------
# the experiments
# ------------------------------------------------------------------------------------------------------------------------------------
FAMS = ('white', 'comp', 'mid')
MULTS = (0.1, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 13.0, 20.0, 32.0, 50.0, 80.0)
CELLS = ([(1, 0.0, 1.0, 'random', s) for s in (11, 12, 13, 14)]
         + [(2, sep, R, mode, 100 * i + s) for i, (sep, R, mode) in enumerate([(6.0, 1.0, 'worst'), (6.0, 1.0, 'random'), (6.0, 2.4, 'worst'), (6.0, 2.4, 'random'),
                                                                                  (8.0, 1.0, 'worst'), (8.0, 1.0, 'random'), (8.0, 2.4, 'worst'), (8.0, 2.4, 'random')]) for s in (21, 22)]
         + [(3, sep, R, mode, 1000 + 100 * i) for i, (sep, R, mode) in enumerate([(7.0, 1.0, 'worst'), (7.0, 2.4, 'random'), (8.0, 1.0, 'random'), (8.0, 2.4, 'worst')])])
G_CELLS = (0, 4, 5, 9, 12, 18, 20)                   # indices into CELLS used by the Gaussian part
G_LEVELS = (0.05, 0.2, 0.5, 1.0)
G_SMALL = (1e-2, 1e-3, 1e-4, 1e-5)
X_TIDS, X_DRAWS = tuple(range(0, 24, 2)), 25
G_P, G_NMC = 0.05, 60
TOL = 1e-6


def _prep(cell):
    K, sep, R, mode, seed = cell
    Kobj, c0, Z = setup()
    cs, As = make_target(K, sep, R, mode, seed)
    P = noiseless_path(cs / SIG, As, c0 / SIG)
    return cs, As, P, float(np.abs(As).min()), target_f(cs, As)


def _ratios(P, S, a, K):
    """The ratios actual/bound of the claims for one noise vector (nan where not defined)."""
    out = dict(order=float(a['order_ok']), terr=np.nan, theta=np.nan, winz=np.nan, win=np.nan, rem=np.nan)
    if not (S['ok'] and a['order_ok'] and 'theta' in a):
        return out
    eta = np.array(S['eta'][1:])
    out['terr'] = float((a['t_err'] / np.maximum(eta, 1e-300)).max())
    out['theta'] = float((a['theta'] / np.maximum(np.array(S['Theta'][1:]), 1e-300)).max())
    w = a['win']
    out['win'] = float(S['W_lo'] / w) if np.isfinite(w) else 0.0
    out['winz'] = float(S.get('W_lo_z', S['W_lo']) / w) if np.isfinite(w) else 0.0
    out['rem'] = float(a['rem'] / max(S['rem'], 1e-300))
    return out


def part_V():
    warnings.filterwarnings('ignore')
    Kobj, c0, Z = setup()
    rows, tgt = [], []
    for tid, cell in enumerate(CELLS):
        cs, As, P, amin, f = _prep(cell)
        K, sep, R, mode, seed = cell
        tb_f = []
        for fi, fam in enumerate(FAMS):
            rng = np.random.default_rng(seed * 7 + fi + 5)
            xi0 = noise_direction(P, cs, fam, rng, amin)
            tb, rb, g0, tab = best_instance(P, cs, xi0, amin)
            tb_f.append(tb)
            ez0 = zeta_gauge(xi0, cs, As, amin)
            for mult in MULTS:
                if tb <= 0:
                    break
                sc = tb * mult
                xi = xi0 * sc
                a = actual(P, cs, As, f, xi, K, amin)
                S = system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, ez=ez0 * sc)
                q = _ratios(P, S, a, K)
                viol = float(S['ok'] and (not a['order_ok'] or max([x for x in (q['terr'], q['theta'], q['winz'], q['win'], q['rem']) if np.isfinite(x)] + [0.0]) > 1 + TOL))
                rows.append([tid, K, sep, R, float(mode == 'worst'), fi, mult, tb, rb, g0[0] * sc, g0[1] * sc, g0[2] * sc, float(ez0.max() * sc), float(S['ok']), float(a['order_ok']), a['n'],
                             q['terr'], q['theta'], q['winz'], q['win'], q['rem'], viol, S.get('W_lo_z', np.nan), a.get('win', np.nan)])
        tgt.append([tid, K, sep, R, float(mode == 'worst'), P['Delta'], max(P['d0']), min(P['mu']), min(P['gap']), float(P['ok']), P['tau'][K],
                    min(abs(P['p'][k][P['ent'][k]]) for k in range(K))] + tb_f)
    return dict(V_rows=np.array(rows, float), V_tgt=np.array(tgt, float))


def part_A():
    """Ball arithmetic: S evaluated with python-flint balls (limit-object data widened by 1e-9 relative) at 0.9 of the certified scale, against the float evaluation."""
    if not HAVE_ARB:
        return dict(A_rows=np.zeros((0, 8)))
    warnings.filterwarnings('ignore')
    rows = []
    for tid, cell in enumerate(CELLS):
        cs, As, P, amin, f = _prep(cell)
        K, sep, R, mode, seed = cell
        for fi, fam in enumerate(FAMS):
            rng = np.random.default_rng(seed * 7 + fi + 5)
            xi0 = noise_direction(P, cs, fam, rng, amin)
            tb, rb, g0, tab = best_instance(P, cs, xi0, amin)
            if tb <= 0:
                continue
            ez0 = zeta_gauge(xi0, cs, As, amin)
            sc = 0.9 * tb
            Sf = system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, ez=ez0 * sc)
            Sa = system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, B=Ar, ez=ez0 * sc)
            rows.append([tid, fi, float(Sf['ok']), float(Sa['ok']), Sa['eta'][-1] / max(Sf['eta'][-1], 1e-300) if Sa['ok'] else np.nan,
                         Sa['W_lo_z'] / Sf['W_lo_z'] if Sa['ok'] else np.nan, Sa['rem'] / Sf['rem'] if Sa['ok'] else np.nan, Sa['lam'] / Sf['lam'] if Sa['ok'] else np.nan])
    return dict(A_rows=np.array(rows, float))


def part_G():
    warnings.filterwarnings('ignore')
    Kobj, c0, Z = setup()
    N, d = Kobj.N, len(c0)
    rows, summ, wsm = [], [], []
    for tid in G_CELLS:
        cell = CELLS[tid]
        cs, As, P, amin, f = _prep(cell)
        K = cell[0]
        _, zn = zeta_family(cs, As, amin)
        best = (0.0, None)
        for r in R_GRID:
            sc = st_crit(P, r, G_P, N, d, zn)
            if sc > best[0]:
                best = (sc, r)
        sc, r = best
        summ.append([tid, K, cell[1], cell[2], sc, r, sc * amin * math.sqrt(N * Z)])
        if sc <= 0:
            continue
        for li, lev in enumerate(G_SMALL):                       # the window at small noise: W^z st tends to min_j (1 - |u_j|)/(||zeta_j|| L)
            st = lev * sc
            E0, E1, E2, Ez = gauss_gauges(st, N, K, d, r, G_P, zn)
            S = system(P, E0, E1, E2, r, ez=Ez)
            Mm = max(10, int(math.ceil(4.0 * r * math.sqrt(N))))
            Lg = math.sqrt(2.0 * math.log(4.0 * (3 * K * (2 * Mm + 1) + 2 * d) / G_P))
            inact = [j for j in range(d) if j not in P['act'][K]]
            lim = min((1 - abs(P['u'][K][j])) / max(zn[j] * Lg, 1e-300) for j in inact)
            wsm.append([tid, lev, st, S.get('W_lo_z', np.nan), S.get('W_lo_z', np.nan) * st, lim])
        for li, lev in enumerate(G_LEVELS):
            st = lev * sc
            E0, E1, E2, Ez = gauss_gauges(st, N, K, d, r, G_P, zn)
            S = system(P, E0, E1, E2, r, ez=Ez)
            rng = np.random.default_rng(900000 + 100 * tid + li)
            for draw in range(G_NMC):
                xi = st * amin * math.sqrt(N * Z) * rng.standard_normal(N)
                g = gauges(xi, cs, c0, r, amin)
                ez = zeta_gauge(xi, cs, As, amin)
                within = float(g[0] <= E0 and g[1] <= E1 and g[2] <= E2 and bool((ez <= Ez).all()))
                a = actual(P, cs, As, f, xi, K, amin)
                q = _ratios(P, S, a, K)
                viol = float(S['ok'] and (not a['order_ok'] or max([x for x in (q['terr'], q['theta'], q['winz'], q['win'], q['rem']) if np.isfinite(x)] + [0.0]) > 1 + TOL))
                rows.append([tid, K, lev, st, E0, E1, E2, float(S['ok']), S.get('W_lo', np.nan), S.get('W_lo_z', np.nan), S.get('eta', [np.nan])[-1], S.get('rem', np.nan), within,
                             g[0] / E0, g[1] / E1, g[2] / E2, float(a['order_ok']), q['terr'], q['theta'], q['winz'], q['rem'], viol, a.get('win', np.nan)])
    return dict(G_rows=np.array(rows, float), G_summ=np.array(summ, float), G_w=np.array(wsm, float))


def part_X():
    """A counterexample search: random adversarial noise (a sum of 1-4 kernels with random off-lattice centres and signs, half of the draws with one kernel within 1.5 widths of a feature)
    scaled to 0.5 and 1.0 of its own certified scale; the claims are checked on every draw."""
    warnings.filterwarnings('ignore')
    Kobj, c0, Z = setup()
    rows = []
    for tid in X_TIDS:
        cell = CELLS[tid]
        cs, As, P, amin, f = _prep(cell)
        K = cell[0]
        rng = np.random.default_rng(555000 + tid)
        for draw in range(X_DRAWS):
            m = int(rng.integers(1, 5))
            cm = rng.uniform(-1.5, 1.5, m)
            if draw % 2 == 0:
                cm[0] = cs[int(rng.integers(0, K))] + rng.uniform(-1.5, 1.5) * SIG
            co = rng.choice([-1.0, 1.0], m) * rng.uniform(0.2, 1.0, m)
            xi0 = Kobj.phi(cm[:, None]) @ co
            tb, rb, g0, tab = best_instance(P, cs, xi0, amin)
            if tb <= 0:
                continue
            ez0 = zeta_gauge(xi0, cs, As, amin)
            for mult in (0.5, 1.0):
                sc = tb * mult
                a = actual(P, cs, As, f, xi0 * sc, K, amin)
                S = system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb, ez=ez0 * sc)
                q = _ratios(P, S, a, K)
                viol = float(S['ok'] and (not a['order_ok'] or max([x for x in (q['terr'], q['theta'], q['winz'], q['win'], q['rem']) if np.isfinite(x)] + [0.0]) > 1 + TOL))
                rows.append([tid, draw, mult, float(S['ok']), float(a['order_ok']), q['terr'], q['theta'], q['winz'], q['win'], q['rem'], viol])
    return dict(X_rows=np.array(rows, float))


B_MUC = (1.0, 10.0)                                   # centre mobilities of the flows of part B (Proposition N1 holds for every positive mobility)


def part_B():
    """Proposition N1 on the noisy entry layer: K = 2 (cells of CELLS with K = 2, the 'random' family), the first unit at its NOISY stage-1 critical point, the second entering from v0 = 1e-4 at its
    lattice kernel (offset d0), the gradient flow of the noisy fit  v' = <r phi>/(Z A_min), d' = mu_c sigma <r d_c phi>/(Z A_min)  (the gradient flow of the fit with constant positive mobilities, w' <- v' and c' = mu_c w <r d_c phi>), at 0.9 of the certified scale: does the
    trajectory stay in the box B (|d| <= r, v in [vlo, vbar]), does it end at the noisy critical point, within Theta_2 of the noiseless block?"""
    from scipy.integrate import solve_ivp
    from scipy.optimize import least_squares
    warnings.filterwarnings('ignore')
    Kobj, c0, Z = setup()
    xs = Kobj.xs[:, 0]
    N = xs.size
    rows = []
    for tid, cell in enumerate(CELLS):
        if cell[0] != 2:
            continue
        cs, As, P, amin, f = _prep(cell)
        rng = np.random.default_rng(cell[4] * 7 + 5)
        xi0 = noise_direction(P, cs, 'white', rng, amin)
        tb, rb, g0, tab = best_instance(P, cs, xi0, amin)
        if tb <= 0:
            continue
        sc = 0.9 * tb
        xi = xi0 * sc
        S = system(P, g0[0] * sc, g0[1] * sc, g0[2] * sc, rb)
        bx = box_noise(Fl, 2, P['R'], P['Delta'], rb, g0[0] * sc, g0[1] * sc, g0[2] * sc)
        y = f + xi
        phi = lambda c: np.exp(-(xs[:, None] - np.asarray(c, float)[None, :]) ** 2 / (2 * SIG ** 2))
        j1, j2 = P['ent'][0], P['ent'][1]

        def grad_res(z, k):
            w, c = z[:k], z[k:]
            return (y - phi(c) @ w) / (math.sqrt(N))

        def jacf(z, k):
            w, c = z[:k], z[k:]
            ph = phi(c)
            dph = ph * (xs[:, None] - c[None, :]) / SIG ** 2
            return -np.concatenate([ph, dph * w[None, :]], 1) / math.sqrt(N)
        z1 = least_squares(lambda z: grad_res(z, 1), np.r_[P['v'][1] * amin, P['c'][1] * SIG], jac=lambda z: jacf(z, 1), xtol=1e-15, ftol=1e-15, gtol=1e-15).x
        z2 = least_squares(lambda z: grad_res(z, 2), np.r_[P['v'][2] * amin, P['c'][2] * SIG], jac=lambda z: jacf(z, 2), xtol=1e-15, ftol=1e-15, gtol=1e-15).x
        feats = P['feat'][2]
        for muc in B_MUC:
            def rhs(t, z):
                w, c = z[:2], z[2:]
                ph = phi(c)
                r = y - ph @ w
                dph = ph * (xs[:, None] - c[None, :]) / SIG ** 2
                gv = (r[:, None] * ph).mean(0) / (Z * amin)
                gd = (r[:, None] * dph).mean(0) * SIG / (Z * amin)
                return np.r_[gv, muc * (w / amin) * gd * SIG]               # c' = mu_c w <r d_c phi> (the network's law, w = A_min v): the gradient flow of the fit with a constant mobility
            zi = np.r_[z1[0], 1e-4 * amin * np.sign(P['a'][feats[1]]), z1[1], c0[j2]]
            sol = solve_ivp(rhs, (0, 4000.0), zi, method='LSODA', rtol=1e-10, atol=1e-13, dense_output=True)
            T_ = np.linspace(0, sol.t[-1], 4000)
            Y = sol.sol(T_)
            dd = np.abs(Y[2:] - np.array([cs[feats[0]], cs[feats[1]]])[:, None]) / SIG
            vv = np.abs(Y[:2]) / amin
            fin = Y[:, -1]
            vlo_f, vbar_f = float(bx['vlo']), float(bx['vbar'])
            above = [int(np.argmax(vv[i] >= vlo_f)) if (vv[i] >= vlo_f).any() else vv.shape[1] for i in range(2)]
            vmin_after = min(vv[i, above[i]:].min() if above[i] < vv.shape[1] else 0.0 for i in range(2))
            ok_in = float((dd <= rb + 1e-9).all() and vmin_after >= vlo_f - 1e-9 and (vv <= vbar_f + 1e-9).all())
            dist = float(np.sqrt((((fin[:2] - z2[:2]) / amin) ** 2).sum() + (((fin[2:] - z2[2:]) / SIG) ** 2).sum()))
            th2 = float(np.sqrt((((fin[:2] / amin) - P['v'][2]) ** 2).sum() + (((fin[2:] - P['c'][2] * SIG) / SIG) ** 2).sum()))
            rows.append([tid, muc, float(bx['ok']), ok_in, dd.max(), vmin_after, vv.max(), vlo_f, vbar_f, rb, dist, th2, S['Theta'][2] if S['ok'] else np.nan, sol.t[-1]])
    return dict(B_rows=np.array(rows, float))


# ------------------------------------------------------------------------------------------------------------------------------------
# report (a pure function of the raw)
# ------------------------------------------------------------------------------------------------------------------------------------
def _q(x, qs=(0, 0.5, 1)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [float(np.quantile(x, q)) for q in qs] if x.size else [float('nan')] * len(qs)


def _f(x, f='%.3g'):
    return 'nan' if not np.isfinite(x) else (f % x)


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return '\n'.join(out)


def report(raw):
    V, Tg = raw['V_rows'], raw['V_tgt']
    A, B, G, Gs = raw['A_rows'], raw['B_rows'], raw['G_rows'], raw['G_summ']
    ok = V[:, 13] == 1
    L = []
    L.append('# Theory 3B: a deterministic, non-asymptotic noise theorem for the limit object, and its Gaussian form\n')
    L.append('Generated by `experiments/theory3b_noise.py` (regenerate: `--from-raw`). The statements, proofs and the status of each claim are in the author\'s working note (outside the repository); '
             'this report contains the **numbers**: the finite system `S` of the theorem evaluated on %d noiseless limit paths (continuum Gram, `n = 1`, `K = 1,2,3`), '
             'validated against the repository\'s path computation (`baselines.iss_path`, the sliding ISS in its limit form) on noisy data in a wide window whose sampled Gram is the continuum one '
             '(window [-3, 3], sample spacing 1/30, dictionary lattice 0.1 on [-1.6, 1.6], width 0.12, N = %d).\n' % (len(Tg), raw['meta'][0]))
    L.append('## The statements (hypotheses exactly)\n')
    L.append('**Theorem N (deterministic; PROVED given the finite system `S`, which is COMPUTER-ASSISTED).** Continuum Gram, `n = 1`, equal walls, `K` Gaussian features of the kernels\' width at separation `Delta`, amplitude ratio `R`, '
             'a finite dictionary, the noiseless limit path simple (Thm S-K (a),(b): unique entrants with the right signs, dormant duals below their walls) and `r` larger than every entrant offset. '
             '`xi` is **any** noise element with gauges (in units of `Z A_min`, lengths in widths) `e0 = sup_{N_r u D}|<xi phi_c>|`, `e1 = sup_{N_r} |<xi d_c phi_c>|`, `e2 = sup_{N_r}|<xi d_c^2 phi_c>|` '
             '(`N_r` = the `r`-neighbourhoods of the features, `D` = the dictionary centres) and `ez_j = |<xi, (I - P_T) phi_j>|` (`P_T` the projection onto the tangent space of the K-unit model at the truth). '
             'If `S(e0,e1,e2,ez; r)` holds (a finite list of inequalities in the noiseless path\'s stages, increasing in the gauges), the limit object on `y = f* + xi` has: '
             '**(i)** the same `K` entries in the same order, no other entry before the K-th, no exit; the k-th entry time within `eta_k = O(e)` of the noiseless one, the k-th block within `Theta_k = O(e)`; '
             '**(ii)** after the K-th entry a window of length `>= W_lo^z = O(1/e)` in which no other unit enters and the block is the K-unit least-squares fit; '
             '**(iii)** in the window `||f_hat - f* - P_T xi|| <= Rem`, `Rem = O(e^2)` explicit (`P_T xi` is the first-order oracle: mean square `2K s^2/N` for white noise). '
             '**Theorem N-G (PROVED given `S`).** For independent sub-Gaussian noise of variance proxy `s^2` (`st = s/(A_min sqrt(N Z))`), with probability `>= 1 - p` the gauges are below explicit `E(st, p) = O(st sqrt(ln(nev/p)))` '
             '(`nev = 3K(2M+1) + 2d` Gaussian events on a net, a chain of Taylor remainders closed by Cauchy-Schwarz and the tail of `||xi||`), so if `S(E)` holds (i)-(iii) hold with probability `>= 1 - p`; the window is `>= min_j (1 - |u_j| - ebar_j)/(st ||zeta_j|| sqrt(2 ln(4 nev/p)) + Rem)`: `~ 1/(s sqrt(ln(nev/p)))`. '
             'All constants are evaluated for floating-point limit-object data (critical points to 1e-13, as in Thm S-T1b); the inequalities of `S` were also evaluated in python-flint ball arithmetic with the data widened by 1e-9 (Part A).\n')
    L.append('## Part V: the claims against the repository\'s path\n')
    nt = len(Tg)
    pairs = len({(int(r[0]), int(r[5])) for r in V})
    cert = {(int(t[0]), fi) for t in Tg for fi in range(3) if t[12 + fi] > 0}
    L.append('%d targets (K = 1: 4; K = 2: 16 at separations 6, 8 widths, amplitude ratios 1 and 2.4, "worst" (feature 0.375 widths off the lattice) and random offsets; K = 3: 4 at 7, 8 widths), three noise families '
             '(white Gaussian; `comp` = the kernel of the runner-up, sign-aligned, the adversarial push on the order; `mid` = a kernel between the features, the push toward a noise unit). '
             'For each (target, family) the certified scale `t_cert` is the largest multiple of the noise vector with `S` true (best of the radii r in %s); the path is then run at %d multiples of `t_cert` '
             '(from 0.1 to 80). **Certified (target, family) pairs: %d of %d** (a pair is not certified when a box condition (L), (C), (N) fails at the target\'s separation, or the margins are below the 1e-9 floor). '
             'Rows with `S` true: **%d**; **violations (S true and a claimed bound exceeded, or the order/entry set differs): %d**.\n'
             % (nt, str(R_GRID), len(MULTS), len(cert), nt * 3, int(ok.sum()), int(V[:, 21].sum())))
    rows = []
    for name, c in (('entry time shift / eta', 16), ('block distance / Theta', 17), ('window: W_lo^z / actual window', 18), ('window: W_lo (crude) / actual', 19), ('oracle remainder / Rem', 20)):
        mn, md, mx = _q(V[ok, c])
        rows.append((name, int((ok & np.isfinite(V[:, c])).sum()), _f(mn), _f(md), _f(mx)))
    L.append('Ratio (actual / bound) over the rows with `S` true (a ratio above 1 would refute the bound; the maximum is the sharpest case):\n')
    L.append(table(('claim', 'rows', 'min', 'median', 'max'), rows) + '\n')
    L.append('How far beyond the certified scale the path keeps the noiseless entries (first multiple of `t_cert` at which the first K entries differ; `inf` = never up to 80):\n')
    rows = []
    for fi, fam in enumerate(FAMS):
        ff = []
        for t in Tg:
            m = (V[:, 0] == t[0]) & (V[:, 5] == fi)
            if not m.any():
                continue
            v = V[m]
            bad = v[v[:, 14] == 0]
            ff.append(bad[0, 6] if len(bad) else np.inf)
        fin = np.array([x for x in ff if np.isfinite(x)])
        rows.append((fam, len(ff), len(fin), _f(float(fin.min())) if fin.size else 'nan', _f(float(np.median(fin))) if fin.size else 'nan', _f(float(fin.max())) if fin.size else 'nan'))
    L.append(table(('family', 'certified pairs', 'order lost within 80x', 'min multiple', 'median', 'max'), rows) + '\n')
    L.append('The certificate is therefore conservative for the entry order by a factor of **%s to %s** (adversarial `comp` family: down to %s).\n'
             % (_f(min(float(r[3]) for r in rows if r[3] != 'nan'), '%.1f'), _f(max(float(r[5]) for r in rows if r[5] != 'nan'), '%.0f'), rows[1][3]))
    L.append('### The certified scale against the margins\n')
    L.append('`e0_cert` = the gauge `e0` of the white-noise vector at `t_cert`; `mu` the smallest runner-up margin `1 - max|u_j|` at an entry; `gap` the smallest relative hitting-time gap; `tau_K` the last entry time; `v_min` the smallest entrant approach speed. '
             '`c = e0_cert tau_K / mu` is the explicit multiple of the margin (the theorem needs `e0 tau_K` below `c mu`).\n')
    rows = []
    for t in Tg:
        m = (V[:, 0] == t[0]) & (V[:, 5] == 0) & (V[:, 6] == 1.0)
        e0c = V[m, 9][0] if m.any() else np.nan
        rows.append((int(t[0]), int(t[1]), _f(t[2], '%.0f'), _f(t[3], '%.1f'), 'worst' if t[4] else 'random', _f(t[5], '%.2f'), _f(t[6], '%.2f'), _f(t[7]), _f(t[8]), _f(t[10]), _f(t[11]),
                     _f(e0c) if e0c > 0 else '0 (box fails)', _f(e0c * t[10] / t[7]) if e0c > 0 else '-'))
    L.append(table(('target', 'K', 'sep', 'R', 'offsets', 'Delta (min sep)', 'max d0', 'mu', 'gap', 'tau_K', 'v_min', 'e0_cert (white)', 'c = e0 tau_K / mu'), rows) + '\n')
    cs_ = np.array([float(r[-1]) for r in rows if r[-1] != '-'])
    L.append('Over the %d certified targets `c = e0_cert tau_K / mu` ranges %s to %s (median %s).\n' % (cs_.size, _f(cs_.min()), _f(cs_.max()), _f(float(np.median(cs_)))))
    L.append('### The window against the noise level\n')
    rows = []
    for mult in (0.1, 0.3, 0.6, 1.0):
        m = ok & (V[:, 6] == mult) & np.isfinite(V[:, 18]) & (V[:, 18] > 0)
        rows.append((mult, int(m.sum()), _f(float(np.median(V[m, 18]))) if m.any() else 'nan', _f(float(np.max(V[m, 18]))) if m.any() else 'nan',
                     _f(float(np.median(V[m, 19]))) if m.any() else 'nan'))
    L.append(table(('multiple of t_cert', 'rows', 'median W_lo^z / actual', 'max', 'median W_lo (crude) / actual'), rows) + '\n')
    L.append('(At the multiple 1.0 the margins are exhausted by construction, so the window bound vanishes there.) As the noise falls the refined window `W_lo^z` approaches the actual window (ratio 1; the bound is sharp to first order: `1 - |u_j|` over the orthogonalised rate), while the crude window, which uses the triangle inequality on the block correction, loses the cancellation.\n')
    L.append('## Part A: ball arithmetic\n')
    if A.shape[0]:
        L.append('`S` evaluated with python-flint balls (limit-object data widened by 1e-9 relative) at 0.9 `t_cert`: **%d of %d** (target, family) pairs are certainly true in ball arithmetic (float evaluation: %d); '
                 'ratios ball/float of `eta_K`, `W_lo^z`, `Rem`, `lam`: within [%s, %s].\n' % (int(A[:, 3].sum()), A.shape[0], int(A[:, 2].sum()), _f(float(np.nanmin(A[:, 4:])), '%.6f'), _f(float(np.nanmax(A[:, 4:])), '%.6f')))
    else:
        L.append('python-flint not available: floats only.\n')
    L.append('## Part B: the noisy entry layer (Proposition N1)\n')
    L.append('K = 2 targets at 0.9 `t_cert` (white noise), the first unit at its noisy critical point, the entrant started at `v0 = 1e-4` at its lattice kernel, the network\'s flow `w\' = <r phi>`, `c\' = mu_c w <r d_c phi>` '
             'with `mu_c = 1, 10`: **%d of %d** runs stay in the box `B` (`|d| <= r`, `vlo <= |v| <= vbar` after the entrant passes `vlo`) and **all** end at the noisy critical point (largest distance %s), '
             'within `Theta_2` of the noiseless block (largest ratio actual/Theta_2: %s).\n' % (int(B[:, 3].sum()), B.shape[0], _f(float(B[:, 10].max())), _f(float(np.nanmax(B[:, 11] / B[:, 12])))))
    L.append('## Part G: sub-Gaussian noise (Theorem N-G), p = %.2f\n' % G_P)
    L.append('For 7 targets: `st_crit` = the largest `st = s/(A_min sqrt(N Z))` with `S(E(st, p))` true, `s_crit/A_min` the corresponding noise standard deviation relative to the smallest amplitude, '
             'then %d white-noise draws at each fraction of `st_crit`. "within" = fraction of draws whose actual gauges (`gauges`, a rigorous supremum) are below `E`; "viol" = draws with `S` true and a bound exceeded.\n' % G_NMC)
    rows = []
    for g in Gs:
        tid = int(g[0])
        for lev in G_LEVELS:
            m = (G[:, 0] == tid) & (G[:, 2] == lev)
            if not m.any():
                continue
            gg = G[m]
            rows.append((tid, int(g[1]), _f(g[2], '%.0f'), _f(g[3], '%.1f'), _f(g[4]), _f(g[4] * math.sqrt(raw['meta'][0] * raw['meta'][1])), lev, _f(gg[0, 3]), _f(gg[0, 9]), _f(gg[0, 9] * gg[0, 3]),
                         _f(float(np.median(gg[:, 22]) * gg[0, 3])), _f(float(gg[:, 12].mean()), '%.2f'), int(gg[:, 21].sum()), _f(float(np.nanmax(gg[:, 17]))), _f(float(np.nanmax(gg[:, 19]))), _f(float(np.nanmax(gg[:, 20])))))
    L.append(table(('target', 'K', 'sep', 'R', 'st_crit', 's_crit/A_min', 'level', 'st', 'W_lo^z(E)', 'W_lo^z st', 'median actual window x st', 'within', 'viol', 'max time ratio', 'max window ratio', 'max rem ratio'), rows) + '\n')
    Gw = raw['G_w']
    L.append('The window at small noise (`W_lo^z(E) st` against its limit `c_W = min_j (1 - |u_j^0|)/(||zeta_j|| sqrt(2 ln(4 nev/p)))`, both from the deterministic system, no path run): the window is `~ c_W/st`, i.e. `~ 1/s`.\n')
    rows = []
    for g in Gs:
        m = Gw[:, 0] == g[0]
        if not m.any():
            continue
        w = Gw[m]
        rows.append([int(g[0]), int(g[1]), _f(g[2], '%.0f'), _f(g[3], '%.1f')] + [_f(x) for x in w[:, 4]] + [_f(w[0, 5])])
    L.append(table(['target', 'K', 'sep', 'R'] + ['W st at st_crit x %g' % x for x in G_SMALL] + ['limit c_W'], rows) + '\n')
    L.append('Total over Part G: %d draws, gauges within `E`: %d, `S` true on all of them, claim violations: %d, entry order different from the noiseless one: %d.\n'
             % (G.shape[0], int(G[:, 12].sum()), int(G[:, 21].sum()), int((G[:, 16] == 0).sum())))
    X = raw['X_rows']
    L.append('## Part X: a counterexample search\n')
    L.append('Random adversarial noise (a sum of 1-4 kernels with random off-lattice centres and random signs; every other draw puts one kernel within 1.5 widths of a feature), %d draws on %d targets, each at 0.5 and 1.0 of its own certified scale '
             '(the certificate is evaluated for each draw with the draw\'s own gauges): **%d runs, `S` true on %d, violations %d, entry order different from the noiseless one %d**. Largest ratios actual/bound: '
             'entry times %s, blocks %s, refined window %s, remainder %s.\n'
             % (X_DRAWS, len(X_TIDS), X.shape[0], int(X[:, 3].sum()), int(X[:, 10].sum()), int((X[:, 4] == 0).sum()), _f(float(np.nanmax(X[:, 5]))), _f(float(np.nanmax(X[:, 6]))), _f(float(np.nanmax(X[:, 7]))), _f(float(np.nanmax(X[:, 9])))))
    return '\n'.join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--parts', default='V,A,B,X,G')
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if a.run:
        Kobj, c0, Z = setup()
        raw = {}
        for pt in a.parts.split(','):
            t0 = time.time()
            raw.update({'V': part_V, 'A': part_A, 'B': part_B, 'G': part_G, 'X': part_X}[pt]())
            print('part', pt, 'done in %.0f s' % (time.time() - t0), flush=True)
        raw['meta'] = np.array([Kobj.N, Z])
        np.savez(RAW, **raw)
    if a.run or a.from_raw:
        d = np.load(RAW)
        open(MD, 'w').write(report({k: d[k] for k in d.files}) + '\n')


if __name__ == '__main__':
    main()
