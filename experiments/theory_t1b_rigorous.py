"""
Interval-arithmetic verification of the certificate S(eps) of Theorem T1 (K <= 3, n = 1, equal masses): python-flint balls, eps may be a BALL (then S is verified for every
real eps in it).

The verification is the standard "find in floating point, verify in balls": `exp_theory_t1.cert` supplies the candidate per-kernel bounds D_j; `verify` re-evaluates
every ingredient of S in ball arithmetic with those D_j and checks the inequalities (6.1), (S1)-(S5) with the upper/lower endpoints.  Replaced by rigorous enclosures:
  * the unit's functions (s(x), w(x), M(x), Psi(x)) and their inversions (bisection on the monotone functions, with three-valued comparisons);
  * the quadratures T_1, A_tr, the pre-layer drift integral (upper Riemann sums of the positive integrand on balls of subintervals);
  * sup |g-hat'| and the oscillation of g-hat (balls of subintervals; no Lipschitz margin needed);
  * the dormant sum's peak (sum of the per-kernel peaks: slightly more conservative than the sampled sum of the float version);
  * the drift double integral (upper right-endpoint sum of the monotone cumulative mass);
  * the weights, the box constants and the Gershgorin margins (balls).
NOT rigorous (inputs taken as exact numbers): the limit object's data (starting centres, the block's critical points, the entry times and duals), which are floating-point
solutions of the limit recursion, as in Theorem K's own certificate, and every final comparison (`float(ball.upper())`).  The hypotheses r > d_0, "the block critical point
lies in the relaxation set D-hat" and theta < v_# are asserted by `hypotheses` (floating-point limit data against ball constants); the bootstrap test is strict,
E_j < D_j / 1.0001.  Guards that used to replace an undecided ball by 0 (T_2, T_3, Y_0) now use the upper end of the ball.
"""
import os
import sys

import numpy as np
from flint import arb, ctx

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T      # noqa: E402

ctx.prec = 100
LAWR = dict(m=12, G=1, gam=(1, 10), amax=1)
SIGMA2_INV = arb(625) / 9                            # sigma = 0.12 = 3/25


def A(x):
    return x if isinstance(x, arb) else arb(x)


def up(x):
    return A(x).upper()


def lo(x):
    return A(x).lower()


def sq(z):
    """z^2 as a nonnegative ball (arb's z**2 / z*z of a ball containing 0 can have a negative lower end, which breaks divisions by q(x))."""
    z = A(z)
    l_, u_ = lo(z), up(z)
    if float(l_) <= 0.0 <= float(u_):
        a_, b_ = up(l_ * l_), up(u_ * u_)
        return arb(0).union(a_ if float(a_) > float(b_) else b_)
    return z * z


def clamp01(x):
    return arb(max(float(lo(x)), 0.0)).union(arb(min(float(up(x)), 1.0)))


def pos_up(x):
    """an upper bound of max(x, 0) as a point ball (x may be undecided: then its upper end is used, never 0)"""
    return arb(max(float(up(x)), 0.0))


def logpos_up(x):
    """an upper bound of log(max(x, 1)) as a nonnegative ball (x undecided around 1: the log of its upper end, never 0)"""
    u_ = float(up(x))
    if u_ <= 1.0:
        return arb(0)
    return arb(0).union(up(x).log())


def bmax_up(a, b):
    """the ball among a, b with the larger upper end (use only for upper bounds)"""
    return a if float(a.upper()) >= float(b.upper()) else b


class RU:
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
        """s(x) = (a_max/(2 G v_s)) (g(y)/tanh x - x), y = 1/(kappa tanh x), g(y) = atan(y)/y  (= kappa-form of Lemma 2.1; written so that a ball for eps does not blow up)"""
        t = x.tanh()
        y = 1 / (self.kap * t)
        if float(up(y)) < 0.5:
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
        """x with f(x) = target for a monotone f, as a ball [lo, hi]; three-valued comparisons."""
        lo_, hi_ = 0.0, 1.0
        tl, tu = float(lo(target)), float(up(target))
        while True:
            fv = f(arb(hi_))
            below = (float(up(fv)) < tl) if decreasing else (float(lo(fv)) > tu)
            if below:
                break
            hi_ *= 2
            if hi_ > 64:
                raise ValueError('no bracket: target [%s, %s], f(64) in [%s, %s]' % (tl, tu, float(lo(fv)), float(up(fv))))
        for _ in range(90):
            mid = 0.5 * (lo_ + hi_)
            fv = f(arb(mid))
            if decreasing:
                gt, lt = float(lo(fv)) > tu, float(up(fv)) < tl
            else:
                gt, lt = float(up(fv)) < tl, float(lo(fv)) > tu           # f(mid) < target  <=> x > mid
            if gt:
                lo_ = mid
            elif lt:
                hi_ = mid
            else:
                break
        return arb(lo_).union(arb(hi_))

    def x_of_s(self, target):
        target = A(target)
        if float(lo(target)) >= float(lo(self.lam)):
            raise ValueError('target above the wall')
        return self._bisect(self.s, target, True)

    def x_of_w(self, target):
        return self._bisect(self.w, A(target), False)

    def B(self, s):
        return self.w(self.x_of_s(s))

    def Psi_s(self, s):
        return self.Psi(self.x_of_s(s))


def Gam_r(x):
    return (-(A(x) * A(x)) / 4).exp()


def psi_r(x):
    return A(x) / 2 * Gam_r(x)


PSI_MAX_R = arb(2).sqrt() / 2 * (arb(-1) / 2).exp()
SQ2, SQ6 = arb(2).sqrt(), arb(6).sqrt()
E15 = (arb(-3) / 2).exp()


def psi_sup_r(X):
    """X: ball >= 0 (a lower bound is used)"""
    Xl = A(lo(X))
    if float(Xl) <= float(SQ2.lower()):
        return PSI_MAX_R
    return psi_r(Xl)


def g2_sup_r(X):
    Xl = A(lo(X))
    if float(Xl) < float(SQ6.lower()):
        v = ((Xl * Xl / 4 - arb(1) / 2).__abs__()) * Gam_r(Xl)
        return v if float(up(v)) > float(up(E15)) else E15
    return ((Xl * Xl / 4 - arb(1) / 2).__abs__()) * Gam_r(Xl)


def rig_box(K, R, Delta, r, Om, dfrac, vsig):
    """Balls of the box constants and Gershgorin margins (as T.box_t1, returning balls); None if a condition is not certain."""
    R, r, Om = A(R), A(r), A(Om)
    out = {}
    if K > 1:
        ell = A(Delta) - 2 * r
        if not (float(lo(ell)) >= float(SQ6.upper())):
            return None
        kap = (K - 1) * Gam_r(ell)
        psl, gl = psi_r(ell), Gam_r(ell)
        G2l = ((ell * ell / 4 - arb(1) / 2).__abs__()) * Gam_r(ell)
    else:
        ell, kap, psl, gl, G2l = None, arb(0), arb(0), arb(0), arb(0)
    if not float(up(kap)) < 1:
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
    conds = [float(lo(vlo)) > 0, float(lo(psi_r(r) - Ec)) > 0, float(up(dt)) < float(lo(r)), float(lo(vlo_d)) > 0, float(lo(bv)) > 0, float(lo(bd)) > 0]
    if not all(conds):
        return None
    out.update(vbar=vbar, Ev=Ev, Ec=Ec, vlo=vlo, dt=dt, dh=dh, vD_lo=vlo_d, vD_hi=R + Ev + A(vsig), bv=bv, bd=bd, ell=ell, psl=psl, gl=gl, vsig=A(vsig), R=R, r=r, Om=Om)
    return out


def rig_rgs(prep, stage, c_lo, c_hi, nsub=200):
    """sup |ghat'| and oscillation of ghat^(stage) on [c_lo, c_hi] by balls of subintervals (rigorous, no Lipschitz margin)."""
    cs, As, L = prep['cs'], prep['As'], prep['L']
    if stage == 0:
        w, cc = [], []
    else:
        ev = L['events'][stage - 1]
        w, cc = ev['w'], ev['c']
    edges = np.linspace(float(c_lo), float(c_hi), nsub + 1)
    supg, gmax, gmin = 0.0, -1e300, 1e300
    for a, b in zip(edges[:-1], edges[1:]):
        c = arb(a).union(arb(b))
        g, gp = arb(0), arb(0)
        for m_ in range(len(cs)):
            x = c - arb(float(cs[m_]))
            e = Gam_r(x)
            g = g + arb(float(As[m_])) * e
            gp = gp - arb(float(As[m_])) * (x / 2) * e
        for i_ in range(len(w)):
            x = c - arb(float(cc[i_]))
            e = Gam_r(x)
            g = g - arb(float(w[i_])) * e
            gp = gp + arb(float(w[i_])) * (x / 2) * e
        supg = max(supg, float(up(abs(gp))))
        gmax, gmin = max(gmax, float(up(g))), min(gmin, float(lo(g)))
    return arb(supg), arb(gmax) - arb(gmin)


def hypotheses(prep, par, Bx, r, dh, vsharp):
    """The hypotheses of the theorem that the float code used to take for granted, now asserted inside the verifier (floating-point limit data against ball constants):
      'r>d0'      the radius exceeds every entrant offset d_0;
      'crit in D' the critical point of the final K-block lies in the relaxation set D-hat (offsets <= d-hat, coefficients in [v_D^-, v_D^+] a_min);
      'theta<v#'  the formed threshold theta is below v_# = v_-/2 (so a unit above v_# is formed)."""
    amin_f = prep['amin']
    ev = prep['L']['events']
    h = {}
    h['r>d0'] = bool(float(lo(r)) > max(prep['d0']))
    fin = ev[-1]
    inD = True
    for i, e in enumerate(ev):
        off = abs(float(fin['c'][i]) - float(prep['cs'][e['feature']]))
        wv_ = abs(float(fin['w'][i]))
        inD = inD and off <= float(lo(dh)) and wv_ >= float(up(amin_f * Bx['vD_lo'])) and wv_ <= float(lo(amin_f * Bx['vD_hi']))
    h['crit in D'] = bool(inD)
    theta = min(T.THETA_REL * float(prep['a'].max()), 0.1 * amin_f)
    h['theta<v#'] = bool(theta < float(lo(vsharp)))
    return h


def verify(prep, eps, par, Dv, nq=1500, ngrid=1, Tfin_tau=None, detail=False):
    """Interval-arithmetic check of S(eps) with the candidate bounds Dv (list over kernels, floats).  eps: float or arb ball.  Returns (ok, info)."""
    K, amin_f, L = prep['K'], prep['amin'], prep['L']
    cs, c0 = prep['cs'], prep['c0']
    J = len(c0)
    ev = L['events']
    info = {}
    if K > 3:
        raise ValueError('K <= 3 only')
    r_f, Om_f, delta1_f, dfrac, vsig, dc_f = par['r'], par['Om'], par['delta1'], par['dfrac'], par['vsig'], par['dc']
    r, Om, d1, dc = A(r_f), A(Om_f), A(delta1_f), A(dc_f)
    amin = A(amin_f)
    a = prep['a']
    R = A(prep['R'])
    Delta = prep['Delta']
    Bx = rig_box(K, R, Delta if np.isfinite(Delta) else 0.0, r, Om, dfrac, vsig)
    if Bx is None:
        return False, dict(why='box')
    U = RU(eps)
    lam, muc = U.lam, U.muc
    tau_end = (1.05 * ev[-1]['tau'] + 1.0) if Tfin_tau is None else Tfin_tau
    taus = L['tau_br']
    Tfin = lam * arb(tau_end)
    tk = [lam * arb(float(t_)) for t_ in taus[1:]]
    vlo, vbar = amin * Bx['vlo'], amin * Bx['vbar']
    Ev, Ec, dt_, dh = Bx['Ev'], Bx['Ec'], Bx['dt'], Bx['dh']
    vsharp = vlo / 2
    # --- layer constants
    x1 = U.x_of_s(d1)
    w1 = U.w(x1)
    xs = U.x_of_w(vsharp)
    gminus = vlo - w1
    if not (float(up(w1)) < float(lo(vsharp)) and float(lo(gminus)) > 0):
        return False, dict(why='w1')
    xa, xb = float(lo(x1)), float(up(xs))
    edges = np.linspace(xa, xb, nq + 1)
    T1, Atr = arb(0), arb(0)
    for a_, b_ in zip(edges[:-1], edges[1:]):
        xi = arb(a_).union(arb(b_))
        f = U.deta_dx(xi) / (vlo - U.w(xi))
        wd_ = arb(b_) - arb(a_)
        T1 = T1 + up(f) * wd_
        Atr = Atr + up(U.w(xi) * f) * wd_
    u0 = A(max(max(prep['d0']), float(up(dt_))))
    kd = muc * vsharp * amin * Gam_r(r) / 2
    ratio = (u0 - dt_) / (dh - dt_)
    T2 = logpos_up(ratio) / kd
    Vt = amin * (Gam_r(dh) - Ev)
    vup_ = amin * (R + Ev)
    mV = U.M(arb(float(lo(xs))))
    mV = mV.lower()
    Y0 = Vt - vsharp
    Y0 = pos_up(Y0)                                  # upper bound of max(Y0, 0) (an undecided ball is not replaced by 0)
    q1 = Y0 / (amin * Bx['vsig'])
    T3lo = logpos_up(q1) / mV
    if float(up(vbar - vup_)) > 0:
        xu = U.x_of_w(vup_)
        mV2 = U.M(arb(float(lo(xu)))).lower()
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
    Mmin = U.M(arb(float(lo(xd)))).lower()
    k1 = Mmin * Bx['bv'].lower()
    k2 = muc * amin * amin * Bx['bd'].lower()
    kappa = k1 if float(up(k1)) < float(up(k2)) else k2
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
    vhat = [A(abs(float(e['w'][-1]))) for e in ev]
    vhat_fin = [A(abs(float(x_))) for x_ in ev[-1]['w']]
    # weights
    wv = [[None] * K for _ in range(J)]
    wd = [[None] * K for _ in range(J)]
    for j in range(J):
        for i in range(K):
            if j in ent:
                k_ = ent.index(j)
                if i == k_:
                    X = arb(0)
                else:
                    X = abs(arb(float(cs[feat[k_]])) - arb(float(cs[feat[i]]))) - 2 * r
            else:
                X = abs(arb(float(c0[j])) - arb(float(cs[feat[i]]))) - r - dc
            Xl = A(max(float(lo(X)), 0.0))
            wv[j][i], wd[j][i] = Gam_r(Xl), psi_sup_r(Xl)
    # --- entrants' pre-layer drift (Lemma 5.3)
    ell_f = float(lo(ell)) if ell is not None else None
    pre_drift = []
    for k in range(K):
        f_ = feat[k]
        uu = float(up(u0))
        Hs, osc = rig_rgs(prep, k, cs[f_] - uu, cs[f_] + uu)
        if K > 1:
            hdev = (K - 1) * (psi_r(ell) * vbar + g2_sup_r(ell) * vbar * 2 * r) + arb(43) / 100 * amin * Om
        else:
            hdev = arb(43) / 100 * amin * Om
        He = Hs + hdev
        x1h = float(up(x1))
        edges = np.concatenate([[0.0], np.geomspace(1e-9, max(x1h, 2e-9), 3000)])
        acc = arb(0)
        for a_, b_ in zip(edges[:-1], edges[1:]):
            xi = arb(a_).union(arb(b_))
            dcx = muc * He * U.Psi(xi) / gminus
            dcx = dcx if float(up(dcx)) < float(up(2 * r)) else 2 * r
            t_ = up(Hs * dcx)
            t_ = t_ if float(up(t_)) < float(up(osc)) else up(osc)
            acc = acc + up(t_ * U.deta_dx(xi) / gminus) * (arb(b_) - arb(a_))
        pre_drift.append(up(acc))
    # --- kernels: pieces, mass, peaks
    def pieces(j):
        """linear pieces (stage, t_a, t_b, dist_a, dist_b, dtau, inv) with inv = dtau/(|u_a| - |u_b|) = dt/ds (independent of lambda) or None for a flat piece"""
        out = []
        for i in range(K + 1):
            tau_a = arb(float(taus[i]))
            tau_b = arb(float(taus[i + 1])) if i < K else arb(tau_end)
            ta = lam * tau_a
            tb_ = lam * tau_b if i < K else Tfin
            ua = arb(float(L['U_br'][i][j]))
            ub = ua + arb(float(L['P'][i][j])) * (tau_b - tau_a)
            sa_pos, sb_pos = float(lo(ua)) > 0, float(lo(ub)) > 0
            sa_neg, sb_neg = float(up(ua)) < 0, float(up(ub)) < 0
            dtau = tau_b - tau_a

            def inv_of(x, y, dt_tau):
                d_ = abs(x) - abs(y)
                if (float(lo(d_)) > 0 or float(up(d_)) < 0) and abs(float(d_.mid())) > 1e-6 * max(abs(float(x.mid())), abs(float(y.mid())), 1e-300):
                    return dt_tau / d_
                return None
            if (sa_pos and sb_neg) or (sa_neg and sb_pos):
                f = abs(ua) / (abs(ua) + abs(ub))
                tc = ta + (tb_ - ta) * f
                out.append((i, ta, tc, lam * (1 - abs(ua)), lam, dtau * f, inv_of(ua, arb(0), dtau * f)))
                out.append((i, tc, tb_, lam, lam * (1 - abs(ub)), dtau * (1 - f), inv_of(arb(0), ub, dtau * (1 - f))))
            else:
                out.append((i, ta, tb_, lam * (1 - abs(ua)), lam * (1 - abs(ub)), dtau, inv_of(ua, ub, dtau)))
        return out
    mass_tot, peak_sum = arb(0), arb(0)
    per = {}
    for j in non:
        Dj = A(Dv[j])
        mj, pk = arb(0), arb(0)
        pcs = []
        for (st, ta, tb_, da, db, dtau, inv) in pieces(j):
            sa, sb = da - Dj, db - Dj
            if not (float(lo(sa)) > 1 and float(lo(sb)) > 1):
                return False, dict(why='dist<=D', j=j)
            Ba, Bb = U.B(sa), U.B(sb)
            pk = bmax_up(pk, bmax_up(Ba, Bb))
            if inv is not None:
                mj = mj + abs(U.Psi_s(sb) - U.Psi_s(sa)) * abs(inv)
            else:
                mj = mj + lam * dtau * bmax_up(Ba, Bb)
            pcs.append((st, ta, tb_, sa, sb, dtau, inv))
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
        pk_ = arb(float(ev[k]['p_e']))
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
    wstart = [tk[i] - Sm[i] for i in range(K)]
    gpts = np.unique(np.r_[np.linspace(0, 1, 25), np.geomspace(1e-6, 0.5, 10), 1 - np.geomspace(1e-6, 0.5, 10)])
    drift = {}
    dcmax = arb(0)
    for j in non:
        pcs, pk = per[j]
        HO = [rig_rgs(prep, k, c0[j] - dc_f, c0[j] + dc_f, 60) for k in range(K + 1)]
        Hst = [up(h_[0] + arb(43) / 100 * amin * omega_pk) for h_ in HO]
        Ost = [up(h_[1] + arb(43) / 100 * amin * omega_pk * 2 * dc) for h_ in HO]
        hd = []
        for i in range(K):
            Xj = abs(arb(float(c0[j])) - arb(float(cs[feat[i]]))) - r - dc
            Xl = A(max(float(lo(Xj)), 0.0))
            hd.append(psi_sup_r(Xl) * Ttv[i] + g2_sup_r(Xl) * Ttd[i])
        Icum, dsum = arb(0), arb(0)
        for (st, ta, tb_, sa, sb, dtau, inv) in pcs:
            ds = sb - sa
            dtp = tb_ - ta
            pts = [(ta, arb(0))]
            for wsi in wstart:
                if float(up(ta)) < float(lo(wsi)) < float(lo(tb_)):
                    pts.append((wsi, clamp01((wsi - ta) / dtp)))
            pts.append((tb_, arb(1)))
            pts = sorted(pts, key=lambda z: float(z[0].mid()))
            for (q0, f0), (q1, f1) in zip(pts[:-1], pts[1:]):
                hinc = pk * (sum((hd[i] for i in range(K) if float(lo(wstart[i])) <= float(up(q0))), arb(0)) + arb(43) / 100 * mass_tot)
                sA, sB = sa * (1 - f0) + sb * f0, sa * (1 - f1) + sb * f1
                dtt_tau = dtau * (f1 - f0)
                prevg = 0.0
                PsA = U.Psi_s(sA)
                dsb = sB - sA
                for g in gpts[1:]:
                    if inv is not None:
                        Ig = Icum + abs(U.Psi_s(sA * (1 - arb(float(g))) + sB * arb(float(g))) - PsA) * abs(inv)
                    else:
                        Ig = Icum + bmax_up(U.B(sA), U.B(sB)) * lam * dtt_tau * arb(float(g))
                    Ig = up(Ig)
                    dcb = muc * (Hst[st] * Ig + hinc)
                    t_ = up(Hst[st] * dcb)
                    t_ = t_ if float(up(t_)) < float(up(Ost[st])) else Ost[st]
                    dsum = dsum + up(t_ * lam * dtt_tau * arb(float(g - prevg)))
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
        if not (float(Etot) < float(Dv[j]) / 1.0001):          # strict, as the bootstrap of THEORY_T1.md section 6 needs
            ok = False
        worst = max(worst, float(Etot) / float(Dv[j]))
    cond = {}
    cond['fixed point E<=D'] = ok
    cond['Omega'] = float(omega_pk) <= float(Om_f)
    cond['dc'] = float(up(dcmax)) <= float(dc_f)
    p = [arb(float(e['p_e'])) for e in ev]
    ordok = float(lo(tk[0])) >= float(up(bmax_up((d1 + De[0]) / p[0], arb(0))))
    for k in range(K - 1):
        lo_ = bmax_up((d1 + De[k + 1]) / p[k + 1], arb(0))
        hi_ = bmax_up((De[k] - d1) / p[k], arb(0)) + Tlay
        ordok = ordok and float(lo(lam * (arb(float(taus[k + 2])) - arb(float(taus[k + 1]))))) >= float(up(lo_ + hi_))
    cond['order'] = bool(ordok)
    theta = min(T.THETA_REL * float(a.max()), 0.1 * amin_f)
    cond['theta'] = bool(max(float(up(per[j][1])) for j in non) < theta and float(lo(amin * Bx['vD_lo'])) > theta)
    cond.update(hypotheses(prep, par, Bx, r, dh, vsharp))
    info.update(cond=cond, Etot=res, worst_ratio=worst, T1=float(up(T1)), T2=float(up(T2)), T3=float(up(T3)), Tlay=float(up(Tlay)), De=[float(x) for x in De], mass=float(mass_tot), omega=float(omega_pk),
                dcmax=float(up(dcmax)), kappa=float(lo(kappa)), Atr=float(up(Atr)), pre=float(up(pre)), prerad=float(pre.rad()), lam=float(lam.mid()))
    return all(cond.values()), info
