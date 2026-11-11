"""
Fully rigorous interval-arithmetic verification of the certificate S(eps) of Theorem T1 (K <= 3, n = 1, equal masses), round theory3, direction C.

Differences from theory_t1b_rigorous.py (referee round T1b, findings 2 and 5):
  * the limit object (block critical points, entry times, duals at the breakpoints, rates, margins, entrant/feature assignment, offsets d_0, sign test) comes from
    `theory3c_limit.limit_rig` as BALLS (Krawczyk enclosure of every block critical point, strict-argmin entry decisions); nothing is read from floating-point code except the target,
    which is taken as the exact binary numbers of its doubles;
  * the hypotheses  r > d_0,  'EVERY stage's block critical point lies in the relaxation set D-hat (offset <= d-hat, |coefficient| in [v_D^-, v_D^+] a_min, sign = sign of the feature)',
    theta < v_#  are checked on the enclosures (certain comparisons);
  * EVERY comparison is a certain arb comparison (True only if it holds for all points of the balls); no float() of a ball enters a decision or an edge; the quadrature edges
    (np.linspace) are outward-rounded (fdown/fup) so that the subintervals cover the integration range; the sum of the sub-intervals' weights is exact;
  * no undecided ball is replaced by 0: every guard uses the upper end of the ball (pos_up, logpos_up, bmax_up); the bisection inversions use certain comparisons;
  * the window-start gating of the drift terms is made conservative: a window that MAY start inside a piece is split off with a clamped fraction, and a sub-piece counts every window
    that may have started before its end; overlapping window-start balls are refused (verification fails) instead of being ordered by their midpoints;
  * the interval chain of eps-balls has outward-rounded endpoints (a segment [lo, hi] is the ball between the double below 10^lo and the double above 10^hi).
Soundness of the verification is that of python-flint/arb (ball arithmetic) and of the paper's hand proof that S(eps) implies the theorem.
"""
import os
import sys
import time

import numpy as np
from flint import arb, ctx

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T      # noqa: E402
import theory3c_limit as LM    # noqa: E402

ctx.prec = 100
SIGMA2_INV = arb(625) / 9                            # sigma = 0.12 = 3/25
THETA_REL_B = arb(0.01).union(arb(1) / 100)          # the formed threshold 0.01 (as the double and as the rational)
TENTH_B = arb(0.1).union(arb(1) / 10)


def A(x):
    return x if isinstance(x, arb) else arb(x)


def up(x):
    return A(x).upper()


def lo(x):
    return A(x).lower()


def fup(x):
    """the smallest double >= the ball's upper end"""
    u = up(x)
    f = float(u)
    if arb(f) < u:
        f = float(np.nextafter(f, np.inf))
    return f


def fdown(x):
    """the largest double <= the ball's lower end"""
    l_ = lo(x)
    f = float(l_)
    if arb(f) > l_:
        f = float(np.nextafter(f, -np.inf))
    return f


def sq(z):
    """z^2 as a nonnegative ball (arb's z*z of a ball containing 0 can have a negative lower end, which breaks divisions by q(x))."""
    z = A(z)
    l_, u_ = lo(z), up(z)
    if l_ <= 0 and u_ >= 0:
        a_, b_ = up(l_ * l_), up(u_ * u_)
        return arb(0).union(a_ if a_ >= b_ else b_).nonnegative_part()
    return z * z


def clamp01(x):
    """a ball containing clamp(x, 0, 1) for every x in the ball"""
    l_, u_ = lo(x), up(x)
    l_ = l_ if l_ > 0 else arb(0)
    l_ = l_ if l_ < 1 else arb(1)
    u_ = u_ if u_ < 1 else arb(1)
    u_ = u_ if u_ > 0 else arb(0)
    return l_.union(u_)


def pos_up(x):
    """an upper bound of max(x, 0) as a point ball (an undecided x gives its upper end, never 0)"""
    u_ = up(x)
    return u_ if u_ > 0 else arb(0)


def logpos_up(x):
    """an upper bound of log(max(x, 1)) as a nonnegative ball (an undecided x gives the log of its upper end, never 0)"""
    u_ = up(x)
    if u_ <= 1:
        return arb(0)
    return arb(0).union(u_.log())


def bmax_up(a, b):
    """the ball among a, b with the larger upper end (use only for upper bounds)"""
    return a if up(a) >= up(b) else b


class RU:
    """The unit (closed forms) in ball arithmetic; eps may be a ball."""

    def __init__(self, eps):
        self.eps = A(eps)
        self.vs = arb(1) / 40
        self.G, self.m, self.amax = arb(1), arb(12), arb(1)
        self.kap = ((2 * self.vs - self.eps) / self.eps).sqrt()
        self.lam = arb.pi() * self.amax * self.kap / (4 * self.G * self.vs)
        self.muc = self.G * self.m * SIGMA2_INV
        self.c1 = self.amax * self.kap / (2 * self.G * self.vs)

    def q(self, x):
        return self.eps + (2 * self.vs - self.eps) * sq(x.tanh())

    def w(self, x):
        return (self.amax / self.m) * x.tanh() * (self.eps + 2 * self.vs * sq(x.sinh()))

    def s(self, x):
        """s(x) = (a_max/(2 G v_s)) (g(y)/tanh x - x), y = 1/(kappa tanh x), g(y) = atan(y)/y  (series 1 - y^2/3 + [0, y^4/5] for y < 1/2, alternating bounds)"""
        t = x.tanh()
        y = 1 / (self.kap * t)
        if up(y) < arb(1) / 2:
            g = (1 - y * y / 3) + (y ** 4 / 5).union(arb(0))
        else:
            g = y.atan() / y
        return self.amax / (2 * self.G * self.vs) * (g / t - x)

    def deta_dx(self, x):
        return self.amax / (self.G * self.q(x))

    def M(self, x):
        q = self.q(x)
        return (self.G / self.m) * q * (q + 4 * self.vs * sq(x.sinh()))

    def Psi(self, x):
        return self.amax ** 2 * sq(x.sinh()) / (2 * self.G * self.m)

    def _bisect(self, f, target, decreasing):
        """x with f(x) = target for every target in the ball, f monotone, as a ball [lo, hi] of doubles; certain comparisons only."""
        lo_, hi_ = 0.0, 1.0
        tl, tu = lo(target), up(target)
        while True:
            fv = f(arb(hi_))
            below = (up(fv) < tl) if decreasing else (lo(fv) > tu)
            if below:
                break
            hi_ *= 2
            if hi_ > 64:
                raise ValueError('no bracket: target [%s, %s], f(64) in [%s, %s]' % (tl, tu, lo(fv), up(fv)))
        for _ in range(90):
            mid = 0.5 * (lo_ + hi_)
            fv = f(arb(mid))
            if decreasing:
                gt, lt = lo(fv) > tu, up(fv) < tl
            else:
                gt, lt = up(fv) < tl, lo(fv) > tu           # f(mid) < target  <=> x > mid
            if gt:
                lo_ = mid
            elif lt:
                hi_ = mid
            else:
                break
        return arb(lo_).union(arb(hi_))

    def x_of_s(self, target):
        target = A(target)
        if lo(target) >= lo(self.lam):
            raise ValueError('target above the wall')
        return self._bisect(self.s, target, True)

    def x_of_w(self, target):
        return self._bisect(self.w, A(target), False)

    def B(self, s):
        return self.w(self.x_of_s(s))

    def Psi_s(self, s):
        return self.Psi(self.x_of_s(s))


def Gam_r(x):
    return (-sq(x) / 4).exp()


def psi_r(x):
    return A(x) / 2 * Gam_r(x)


PSI_MAX_R = arb(2).sqrt() / 2 * (arb(-1) / 2).exp()
SQ2, SQ6 = arb(2).sqrt(), arb(6).sqrt()
E15 = (arb(-3) / 2).exp()


def psi_sup_r(X):
    """sup of psi over [X, inf) for X >= 0, X a ball (its lower end is used); psi increases on [0, sqrt 2], decreases beyond"""
    Xl = lo(X)
    if not (Xl > up(SQ2)):
        return PSI_MAX_R
    return psi_r(Xl)


def g2_sup_r(X):
    """sup of |Gamma''| over [X, inf): |Gamma''(x)| = |x^2/4 - 1/2| exp(-x^2/4), decreasing beyond sqrt 6 where it is <= 0.2231 = e^(-3/2)"""
    Xl = lo(X)
    v = abs(sq(Xl) / 4 - arb(1) / 2) * Gam_r(Xl)
    if not (Xl >= up(SQ6)):
        return v if v.upper() > E15.upper() else E15
    return v


def rig_box(K, R, Delta, r, Om, dfrac, vsig):
    """Balls of the box constants and Gershgorin margins (as T.box_t1); None if a condition is not certain."""
    R, r, Om = A(R), A(r), A(Om)
    out = {}
    if not (up(r) < lo(SQ2)):                                   # r < sqrt 2 (THEORY_T1.md 3.0)
        return None
    if K > 1:
        ell = A(Delta) - 2 * r
        if not (lo(ell) >= up(SQ6)):
            return None
        kap = (K - 1) * Gam_r(ell)
        psl, gl = psi_r(ell), Gam_r(ell)
        G2l = abs(sq(ell) / 4 - arb(1) / 2) * Gam_r(ell)
    else:
        ell, kap, psl, gl, G2l = None, arb(0), arb(0), arb(0), arb(0)
    if not (up(kap) < 1):
        return None
    vbar = (R * (1 + kap) + Om) / (1 - kap)
    Ev = kap * (R + vbar) + Om
    Ec = (K - 1) * (R + vbar) * psl + arb(43) / 100 * Om
    vlo = Gam_r(r) - Ev
    dt = 2 * Ec / Gam_r(r)
    dh = dt + A(dfrac) * (r - dt)
    vlo_d = Gam_r(dh) - Ev - A(vsig)
    s_max = R * psi_r(dh)
    t_min = vlo_d * (arb(1) / 2 - dh * dh / 4) * Gam_r(dh)
    rowv = (K - 1) * (gl + (R + vbar) * psl + vbar * psl) + arb(43) / 100 * Om
    rowd = (K - 1) * ((R + vbar) * psl + vbar * (R + vbar) * G2l + vbar * psl + vbar * vbar * G2l) + arb(43) / 100 * Om + vbar * Om / 2
    bv, bd = 1 - s_max - rowv, t_min - s_max - rowd
    conds = [lo(vlo) > 0, lo(psi_r(r) - Ec) > 0, up(dt) < lo(r), lo(vlo_d) > 0, lo(bv) > 0, lo(bd) > 0]
    if not all(conds):
        return None
    out.update(vbar=vbar, Ev=Ev, Ec=Ec, vlo=vlo, dt=dt, dh=dh, vD_lo=vlo_d, vD_hi=R + Ev + A(vsig), bv=bv, bd=bd, ell=ell, psl=psl, gl=gl, vsig=A(vsig), R=R, r=r, Om=Om,
               rowv=rowv, rowd=rowd, s_max=s_max, t_min=t_min)
    return out


def rig_rgs(cs, As, LR, stage, c_lo, c_hi, nsub=200):
    """sup |ghat'| and the oscillation of ghat^(stage) on [c_lo, c_hi] (balls) by balls of subintervals.  The edges are doubles that cover the range (outward rounded);
    the limit's block (w, c) at the stage is the ball enclosure of the critical point."""
    if stage == 0:
        w, cc = [], []
    else:
        ev = LR['events'][stage - 1]
        w, cc = ev['w'], ev['c']
    a, b = fdown(c_lo), fup(c_hi)
    edges = np.linspace(a, b, nsub + 1)
    edges[0], edges[-1] = a, b
    supg, gmax, gmin = arb(0), None, None
    for x0, x1 in zip(edges[:-1], edges[1:]):
        c = arb(float(x0)).union(arb(float(x1)))
        g, gp = arb(0), arb(0)
        for m_ in range(len(cs)):
            x = c - cs[m_]
            e = Gam_r(x)
            g = g + As[m_] * e
            gp = gp - As[m_] * (x / 2) * e
        for i_ in range(len(w)):
            x = c - cc[i_]
            e = Gam_r(x)
            g = g - w[i_] * e
            gp = gp + w[i_] * (x / 2) * e
        supg = bmax_up(supg, up(abs(gp)))
        gmax = up(g) if gmax is None else bmax_up(gmax, up(g))
        gmin = lo(g) if gmin is None else (gmin if gmin <= lo(g) else lo(g))
    return supg, up(gmax) - gmin


def dhat_unique(LR, Bx, amin, As_f, k):
    """Krawczyk operator on the whole relaxation set D-hat of the k-block (a box containing it: offsets |d_i| <= upper(d-hat), coefficients between amin v_D^- and amin v_D^+ with the sign of the
    feature): K(D-hat) in the interior of D-hat proves that the critical equations of the k-block have exactly ONE root in D-hat (the block critical point used by the theorem, which is the
    enclosed one when the enclosure lies in D-hat)."""
    ev = LR['events']
    feats = [e['feature'] for e in ev[:k]]
    wr = lo(amin * Bx['vD_lo']).union(up(amin * Bx['vD_hi']))
    dr = (-up(Bx['dh'])).union(up(Bx['dh']))
    X = [wr if As_f[f] > 0 else -wr for f in feats] + [LR['cs'][f] + dr for f in feats]
    try:
        _, inside = LM.krawczyk(X, LR['cs'], LR['As'])
    except Exception:
        return False
    return bool(inside)


def hypotheses_rig(LR, cs_f, As_f, Bx, r, dh, vsharp, amin, theta):
    """The hypotheses of the theorem, on the ENCLOSURES (certain comparisons):
      'r>d0'        r exceeds every entrant offset d_0;
      'sign'        sgn(p_e) = sgn(A_f) for every entrant (from the limit enclosure);
      'crit in D'   at EVERY stage k the critical point of the k-block lies in the relaxation set D-hat: offsets |c_i - c*_f(i)| <= d-hat, |w_i| in [v_D^-, v_D^+] a_min,
                    sgn w_i = sgn A_f(i);
      'unique in D' the Krawczyk operator maps D-hat (a box containing it) into its interior for every stage's block: ONE critical point of the k-block lies in D-hat;
      'theta<v#'    the formed threshold theta is below v_# = v_-/2."""
    ev = LR['events']
    h = {}
    h['r>d0'] = bool(all(r > e['d0'] for e in ev))
    h['sign'] = bool(all(e['sign_ok'] for e in ev))
    inD = True
    for st in ev:
        feats = [e['feature'] for e in ev[:st['k']]]
        for i, f in enumerate(feats):
            off = abs(st['c'][i] - LR['cs'][f])
            wv_ = abs(st['w'][i])
            sg_ok = (st['w'][i] > 0) if As_f[f] > 0 else (st['w'][i] < 0)
            inD = inD and bool(off <= dh) and bool(wv_ >= amin * Bx['vD_lo']) and bool(wv_ <= amin * Bx['vD_hi']) and bool(sg_ok)
    h['crit in D'] = bool(inD)
    h['unique in D'] = bool(all(dhat_unique(LR, Bx, amin, As_f, st['k']) for st in ev))
    h['theta<v#'] = bool(theta < vsharp)
    return h


def verify_rig(prep, LR, eps, par, Dv, nq=1500, Tfin_tau=None):
    """Interval-arithmetic check of S(eps) with the candidate bounds Dv (list over kernels, doubles) and the BALL limit object LR (theory3c_limit.limit_rig).  eps: float or ball.
    Returns (ok, info)."""
    K = prep['K']
    cs_f, c0_f, As_f = prep['cs'], prep['c0'], prep['As']
    cs = [arb(float(x)) for x in cs_f]
    c0 = [arb(float(x)) for x in c0_f]
    As = [arb(float(x)) for x in As_f]
    J = len(c0)
    if K > 3:
        raise ValueError('K <= 3 only')
    if not LR['ok']:
        return False, dict(why='limit object not certified: ' + LR['why'])
    ev = LR['events']
    info = {}
    r, Om, d1, dc = A(par['r']), A(par['Om']), A(par['delta1']), A(par['dc'])
    dfrac, vsig = par['dfrac'], par['vsig']
    absA = np.abs(np.asarray(As_f, float))
    amin = arb(float(absA.min()))
    amax = arb(float(absA.max()))
    R = amax / amin
    if K > 1:
        Delta = None
        for i in range(K):
            for j in range(i + 1, K):
                d_ = abs(cs[i] - cs[j])
                Delta = d_ if Delta is None else Delta.min(d_)
    else:
        Delta = arb(0)
    Bx = rig_box(K, R, Delta, r, Om, dfrac, vsig)
    if Bx is None:
        return False, dict(why='box')
    U = RU(eps)
    lam, muc = U.lam, U.muc
    taus = LR['tau_br']
    tau_end = arb(fup(arb(105) / 100 * ev[-1]['tau'] + 1)) if Tfin_tau is None else A(Tfin_tau)
    Tfin = lam * tau_end
    tk = [lam * t_ for t_ in taus[1:]]
    vlo, vbar = amin * Bx['vlo'], amin * Bx['vbar']
    Ev, Ec, dt_, dh = Bx['Ev'], Bx['Ec'], Bx['dt'], Bx['dh']
    vsharp = vlo / 2
    # --- layer constants
    x1 = U.x_of_s(d1)
    w1 = U.w(x1)
    xs = U.x_of_w(vsharp)
    gminus = vlo - w1
    if not (up(w1) < lo(vsharp) and lo(gminus) > 0):
        return False, dict(why='w1')
    edges = np.linspace(fdown(x1), fup(xs), nq + 1)
    edges[0], edges[-1] = fdown(x1), fup(xs)
    T1, Atr = arb(0), arb(0)
    for a_, b_ in zip(edges[:-1], edges[1:]):
        xi = arb(float(a_)).union(arb(float(b_)))
        f = U.deta_dx(xi) / (vlo - U.w(xi))
        wd_ = arb(float(b_)) - arb(float(a_))
        T1 = T1 + up(f) * wd_
        Atr = Atr + up(U.w(xi) * f) * wd_
    d0max = ev[0]['d0']
    for e in ev[1:]:
        d0max = bmax_up(d0max, e['d0'])
    u0 = bmax_up(up(d0max), up(dt_))
    kd = muc * vsharp * amin * Gam_r(r) / 2
    ratio = (u0 - dt_) / (dh - dt_)
    T2 = logpos_up(ratio) / kd
    Vt = amin * (Gam_r(dh) - Ev)
    vup_ = amin * (R + Ev)
    mV = U.M(lo(xs)).lower()
    Y0 = pos_up(Vt - vsharp)
    q1 = Y0 / (amin * Bx['vsig'])
    T3lo = logpos_up(q1) / mV
    if not (up(vbar - vup_) <= 0):
        xu = U.x_of_w(vup_)
        mV2 = U.M(lo(xu)).lower()
        Y0p = pos_up(vbar - vup_)
        q2 = Y0p / (amin * Bx['vsig'])
        T3hi = logpos_up(q2) / mV2
        envhi = Y0p / mV2
    else:
        T3hi, envhi = arb(0), arb(0)
    T3 = bmax_up(T3lo, T3hi)
    Wv0 = amin * (R + 2 * Ev - Gam_r(dh))
    Env = Y0 / mV + envhi + Wv0 * T3
    Tlay = T1 + T2 + T3
    # --- stability constants
    xd = U.x_of_w(amin * Bx['vD_lo'])
    Mmin = U.M(lo(xd)).lower()
    k1 = Mmin * lo(Bx['bv'])
    k2 = muc * amin * amin * lo(Bx['bd'])
    kappa = k1 if up(k1) < up(k2) else k2
    Vdiam = bmax_up(2 * dh, Bx['vD_hi'] - Bx['vD_lo'])
    WvD = amin * (Bx['vD_hi'] - Bx['vD_lo'])
    vbp = Bx['vbar']
    ell = Bx['ell']
    pre = U.Psi(x1) / gminus
    cF = (Bx['gl'] / Bx['bv'] + vbp * Bx['psl'] / Bx['bd']) if K > 1 else arb(0)
    cN = 1 / Bx['bv'] + arb(43) / 100 * vbp / Bx['bd']
    ent = [e['kernel'] for e in ev]
    non = [j for j in range(J) if j not in ent]
    feat = [e['feature'] for e in ev]
    vhat = [abs(e['w'][-1]) for e in ev]
    vhat_fin = [abs(x_) for x_ in ev[-1]['w']]
    # weights
    wv = [[None] * K for _ in range(J)]
    wd = [[None] * K for _ in range(J)]
    for j in range(J):
        for i in range(K):
            if j in ent:
                k_ = ent.index(j)
                X = arb(0) if i == k_ else abs(cs[feat[k_]] - cs[feat[i]]) - 2 * r
            else:
                X = abs(c0[j] - cs[feat[i]]) - r - dc
            Xl = lo(X) if lo(X) > 0 else arb(0)
            wv[j][i], wd[j][i] = Gam_r(Xl), psi_sup_r(Xl)
    # --- entrants' pre-layer drift (Lemma 5.3)
    pre_drift = []
    x1h = fup(x1)
    for k in range(K):
        f_ = feat[k]
        Hs, osc = rig_rgs(cs, As, LR, k, cs[f_] - u0, cs[f_] + u0)
        if K > 1:
            hdev = (K - 1) * (psi_r(ell) * vbar + g2_sup_r(ell) * vbar * 2 * r) + arb(43) / 100 * amin * Om
        else:
            hdev = arb(43) / 100 * amin * Om
        He = Hs + hdev
        edg = np.concatenate([[0.0], np.geomspace(1e-9, max(x1h, 2e-9), 3000)])
        edg[-1] = max(x1h, 2e-9)
        acc = arb(0)
        for a_, b_ in zip(edg[:-1], edg[1:]):
            xi = arb(float(a_)).union(arb(float(b_)))
            dcx = muc * He * U.Psi(xi) / gminus
            dcx = dcx if up(dcx) < up(2 * r) else 2 * r
            t_ = up(Hs * dcx)
            t_ = t_ if up(t_) < up(osc) else up(osc)
            acc = acc + up(t_ * U.deta_dx(xi) / gminus) * (arb(float(b_)) - arb(float(a_)))
        pre_drift.append(up(acc))

    # --- kernels: pieces, mass, peaks
    def inv_of(x, y, dt_tau):
        d_ = abs(x) - abs(y)
        if (lo(d_) > 0 or up(d_) < 0) and abs(float(d_.mid())) > 1e-6 * max(abs(float(x.mid())), abs(float(y.mid())), 1e-300):
            return dt_tau / d_
        return None

    def pieces(j):
        """linear pieces (stage, t_a, t_b, dist_a, dist_b, dtau, inv, back) with inv = dtau/(|u_a| - |u_b|) = dt/ds (independent of lambda) or None for a flat piece;
        back = lambda x (tau_stage_end - tau_b): the time from the end of the piece to the end of its stage (0 for the last piece of the stage), as an expression free of the eps-dependency of lambda tau"""
        out = []
        for i in range(K + 1):
            tau_a = taus[i]
            tau_b = taus[i + 1] if i < K else tau_end
            ta = lam * tau_a
            tb_ = lam * tau_b if i < K else Tfin
            ua = LR['U_br'][i][j]
            ub = ua + LR['P'][i][j] * (tau_b - tau_a)
            sa_pos, sb_pos = lo(ua) > 0, lo(ub) > 0
            sa_neg, sb_neg = up(ua) < 0, up(ub) < 0
            dtau = tau_b - tau_a
            if (sa_pos and sb_neg) or (sa_neg and sb_pos):
                f = abs(ua) / (abs(ua) + abs(ub))
                tc = ta + (tb_ - ta) * f
                out.append((i, ta, tc, lam * (1 - abs(ua)), lam, dtau * f, inv_of(ua, arb(0), dtau * f), lam * dtau * (1 - f)))
                out.append((i, tc, tb_, lam, lam * (1 - abs(ub)), dtau * (1 - f), inv_of(arb(0), ub, dtau * (1 - f)), arb(0)))
            else:
                out.append((i, ta, tb_, lam * (1 - abs(ua)), lam * (1 - abs(ub)), dtau, inv_of(ua, ub, dtau), arb(0)))
        return out
    mass_tot, peak_sum = arb(0), arb(0)
    per = {}
    for j in non:
        Dj = A(Dv[j])
        mj, pk = arb(0), arb(0)
        pcs = []
        for (st, ta, tb_, da, db, dtau, inv, back) in pieces(j):
            sa, sb = da - Dj, db - Dj
            if not (lo(sa) > 1 and lo(sb) > 1):
                return False, dict(why='dist<=D', j=j)
            Ba, Bb = U.B(sa), U.B(sb)
            pk = bmax_up(pk, bmax_up(Ba, Bb))
            if inv is not None:
                mj = mj + abs(U.Psi_s(sb) - U.Psi_s(sa)) * abs(inv)
            else:
                mj = mj + lam * dtau * bmax_up(Ba, Bb)
            pcs.append((st, ta, tb_, sa, sb, dtau, inv, back))
        mass_tot = mass_tot + up(mj)
        peak_sum = peak_sum + up(pk)
        per[j] = (pcs, up(pk))
    mass_tot = up(mass_tot)
    omega_pk = up(peak_sum / amin)
    Vtail = Vdiam / kappa + cF * K * pre / amin + cN * omega_pk * Tfin
    # --- windows and unit L1 bounds
    cum = {j: arb(0) for j in range(J)}
    De, Wk, Sp, Sm, Tvl, Tdl = [], [], [], [], [], []
    for k in range(K):
        j = ent[k]
        De_k = cum[j] + sum((wv[j][i] * pre for i in range(k, K)), arb(0)) + mass_tot + pre_drift[k]
        De.append(up(De_k))
        pk_ = ev[k]['p_e']
        Spk = bmax_up((De[k] - d1) / pk_, arb(0))
        Smk = bmax_up((De[k] + d1) / pk_, arb(0))
        Sp.append(Spk)
        Sm.append(Smk)
        Wlen = (abs(d1) + De[k]) / pk_ + Tlay
        Wk.append(Wlen)
        Tv = pre + Atr + vhat[k] * T1 + vbar * bmax_up(Spk, Smk) + vbar * T2 + Env
        Td = (r + dt_) * Atr + vhat[k] * (r + dt_) * T2 + vhat[k] * 2 * dh * T3 + (r + dt_) * w1 * Spk
        Tvl.append(Tv + amin * Vtail)
        Tdl.append(Td + vhat[k] * Vtail)
        for jj in range(J):
            c_ = wv[jj][k] * Tv + wd[jj][k] * Td + (wv[jj][k] * amin + wd[jj][k] * vhat[k]) * Vtail
            for i in range(k):
                c_ = c_ + (wv[jj][i] * WvD + wd[jj][i] * vhat_fin[i] * 2 * dh) * Wlen + (wv[jj][i] * amin + wd[jj][i] * vhat_fin[i]) * Vtail
            cum[jj] = cum[jj] + c_
    Ttv = [Tvl[i] + sum((WvD * Wk[k] + amin * Vtail for k in range(i + 1, K)), arb(0)) for i in range(K)]
    Ttd = [Tdl[i] + sum((vhat_fin[i] * 2 * dh * Wk[k] + vhat_fin[i] * Vtail for k in range(i + 1, K)), arb(0)) for i in range(K)]
    # --- drift of the non-entrants
    gpts = np.unique(np.r_[np.linspace(0, 1, 25), np.geomspace(1e-6, 0.5, 10), 1 - np.geomspace(1e-6, 0.5, 10)])
    drift = {}
    dcmax = arb(0)
    for j in non:
        pcs, pk = per[j]
        HO = [rig_rgs(cs, As, LR, k, c0[j] - dc, c0[j] + dc, 60) for k in range(K + 1)]
        Hst = [up(h_[0] + arb(43) / 100 * amin * omega_pk) for h_ in HO]
        Ost = [up(h_[1] + arb(43) / 100 * amin * omega_pk * 2 * dc) for h_ in HO]
        hd = []
        for i in range(K):
            Xj = abs(c0[j] - cs[feat[i]]) - r - dc
            Xl = lo(Xj) if lo(Xj) > 0 else arb(0)
            hd.append(psi_sup_r(Xl) * Ttv[i] + g2_sup_r(Xl) * Ttd[i])
        Icum, dsum = arb(0), arb(0)
        for (st, ta, tb_, sa, sb, dtau, inv, back) in pcs:
            # sub-pieces by the start of the one window that can start inside this piece (window st, since the windows are ordered: condition 'order').  Window i < st has
            # started before the piece, window i > st starts after it.  back = time from the piece end to the stage end; Sm[st] = time from the stage end back to the window's start.
            seg = lam * dtau
            subs = []                                             # (f0, f1, number of windows counted: i < nw)
            if st >= K:
                subs.append((arb(0), arb(1), K))
            else:
                dback = Sm[st] - back
                if lo(dback) > up(seg):                           # the window started before the piece
                    subs.append((arb(0), arb(1), st + 1))
                elif up(dback) < 0:                               # the window starts after the piece
                    subs.append((arb(0), arb(1), st))
                else:                                             # it may start inside: split at the (clamped) fraction
                    fs = clamp01(1 - dback / seg)
                    subs.append((arb(0), fs, st))
                    subs.append((fs, arb(1), st + 1))
            for (f0, f1, nw) in subs:
                dtt_tau = dtau * (f1 - f0)
                if dtt_tau.is_zero():
                    continue
                hinc = pk * (sum((hd[i] for i in range(nw)), arb(0)) + arb(43) / 100 * mass_tot)
                sA, sB = sa * (1 - f0) + sb * f0, sa * (1 - f1) + sb * f1
                prevg = 0.0
                PsA = U.Psi_s(sA)
                for g in gpts[1:]:
                    ga = arb(float(g))
                    if inv is not None:
                        Ig = Icum + abs(U.Psi_s(sA * (1 - ga) + sB * ga) - PsA) * abs(inv)
                    else:
                        Ig = Icum + bmax_up(U.B(sA), U.B(sB)) * lam * dtt_tau * ga
                    Ig = up(Ig)
                    dcb = muc * (Hst[st] * Ig + hinc)
                    t_ = up(Hst[st] * dcb)
                    t_ = t_ if up(t_) < up(Ost[st]) else Ost[st]
                    dsum = dsum + up(t_ * lam * dtt_tau * (ga - arb(float(prevg))))
                    dcmax = bmax_up(dcmax, dcb)
                    prevg = g
                    last = Ig
                Icum = up(last)
        drift[j] = up(dsum)
    # --- totals and conditions
    res = {}
    ok = True
    worst = 0.0
    for j in non:
        Etot = up(cum[j] + mass_tot + drift[j])
        res[j] = float(Etot)
        if not (Etot * (arb(10001) / 10000) < A(Dv[j])):         # strict, as the bootstrap of THEORY_T1.md section 6 needs
            ok = False
        worst = max(worst, float(Etot) / float(Dv[j]))
    cond = {}
    cond['fixed point E<D/1.0001'] = bool(ok)
    cond['Omega'] = bool(omega_pk <= Om)
    cond['dc'] = bool(up(dcmax) <= dc)
    p = [e['p_e'] for e in ev]
    ordok = bool(tk[0] >= bmax_up((d1 + De[0]) / p[0], arb(0)))
    for k in range(K - 1):
        lo_ = bmax_up((d1 + De[k + 1]) / p[k + 1], arb(0))
        hi_ = bmax_up((De[k] - d1) / p[k], arb(0)) + Tlay
        ordok = ordok and bool(lam * (taus[k + 2] - taus[k + 1]) >= up(lo_ + hi_))
    cond['order'] = bool(ordok)
    theta = (THETA_REL_B * amax).min(TENTH_B * amin)
    pkmax = None
    for j in non:
        pkmax = per[j][1] if pkmax is None else bmax_up(pkmax, per[j][1])
    cond['theta'] = bool(pkmax < theta and amin * Bx['vD_lo'] > theta) if pkmax is not None else bool(amin * Bx['vD_lo'] > theta)
    cond.update(hypotheses_rig(LR, cs_f, np.asarray(As_f, float), Bx, r, dh, vsharp, amin, theta))
    # (N) of THEORY_T1.md 3.0 is not used by the final proof (Lemma 4.2 is not used); reported for completeness
    lam_own = (1 + Bx['t_min'] - ((1 - Bx['t_min']) ** 2 + 4 * Bx['s_max'] * Bx['s_max']).sqrt()) / 2
    Hc = Bx['rowv'] if up(Bx['rowv']) >= up(Bx['rowd']) else Bx['rowd']
    info.update(cond=cond, Etot=res, worst_ratio=worst, T1=float(up(T1)), T2=float(up(T2)), T3=float(up(T3)), Tlay=float(up(Tlay)), De=[float(x) for x in De], mass=float(mass_tot),
                omega=float(omega_pk), dcmax=float(up(dcmax)), kappa=float(lo(kappa)), Atr=float(up(Atr)), pre=float(up(pre)), prerad=float(pre.rad()), lam=float(lam.mid()),
                N_condition=bool(lam_own - Hc > 0), theta=float(up(theta)), tau_end=float(tau_end.mid()),
                drift=[float(drift[j]) for j in non], cum=[float(up(cum[j])) for j in non], wk=[float(up(w_)) for w_ in Wk], vtail=float(up(Vtail)), env=float(up(Env)))
    return all(cond.values()), info


# ----------------------------------------------------------------------------------------------------------------------------------
# the cover of [10^bot, 10^top] by eps-balls
# ----------------------------------------------------------------------------------------------------------------------------------
def pow10(e):
    return arb(10) ** arb(float(e))


def eps_ball(e_lo, e_hi):
    """the ball [double below 10^e_lo, double above 10^e_hi] (outward rounded)"""
    return arb(fdown(pow10(e_lo))).union(arb(fup(pow10(e_hi))))


def cover(prep, LR, e_top, e_bot, wmin=1.0 / 512, budget=None, facs=(1.1, 1.5, 2.0, 4.0)):
    """Interval-arithmetic cover of eps in [10^e_bot, 10^e_top] by balls on which the ball-valid certificate holds (candidates D_j from the floating-point `cert` at the upper end of each
    ball, enlarged by the factors in facs; the ball width adapts; the balls overlap at their outward-rounded ends).  Returns (list of (lo, hi), complete, lowest exponent reached, number of
    verify calls)."""
    cur, segs, calls = e_top, [], 0
    t_end = None if budget is None else time.time() + budget
    wlast = 0.5
    ent = [e['kernel'] for e in LR['events']]
    while cur > e_bot + 1e-9:
        if t_end is not None and time.time() > t_end:
            return segs, False, cur, calls
        o, pi = T.best_cert(prep, 10.0 ** cur)
        if pi < 0:
            return segs, False, cur, calls
        par = T.PARS[pi]
        w, done = min(1.0, 2.0 * wlast), False
        while w >= wmin and not done:
            if t_end is not None and time.time() > t_end:
                return segs, False, cur, calls
            lo_ = max(cur - w, e_bot)
            ball = eps_ball(lo_, cur)
            for fac in facs:
                Dv = [1.0 if j in ent else max(fac * o['Etot'][j], 1.0) for j in range(len(prep['c0']))]
                calls += 1
                try:
                    ok, info = verify_rig(prep, LR, ball, par, Dv)
                except (ValueError, ZeroDivisionError):
                    ok = False
                if ok:
                    segs.append((lo_, cur))
                    cur, done = lo_, True
                    wlast = w
                    break
            w /= 2
        if not done:
            return segs, False, cur, calls
    return segs, True, cur, calls
