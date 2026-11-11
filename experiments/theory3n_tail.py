"""
Theorem T1, the small-eps tail on the ENCLOSED limit object (round theory3, direction N).

This is the tail verifier of round theory3 A (`theory3a_tail.py`: one ball h = 1/kappa in [0, h_hi] covers every real eps in (0, eps_hi]; the dormant side is bounded through the eps-free envelope of Lemma
12.4, B_eps(s) <= C/(s (s^2 - sigma_0^2)), s >= 25) ported to the ball limit object of round theory3 C (`theory3c_limit.limit_rig`) and to the certain comparisons of `theory3c_rigorous`:

  * pieces, margins c_min, the rates p_e, the coefficients vhat, the entry times and the dual breakpoints are read from the ENCLOSURES of the limit object (balls), not from floats;
  * the hypotheses (r > d_0, every stage's critical point in D-hat and unique there, sgn p = sgn A, theta < v_#) are checked on the enclosures (`theory3c_rigorous.hypotheses_rig`);
  * every decision is a certain ball comparison (no float() of a ball);  g(y) = atan(y)/y is enclosed with certain case distinction and outward-rounded evaluation points;
  * the corrected accounting (f), (g), (6), kappa by the smaller LOWER end, as in `theory3n_system` (the displacement is stage-resolved; the pre-layer coefficient of a window not yet started
    is kept in hinc;  gates are the EARLIEST possible start of each window, tau_{i+1} - Sm_i L_hi, rounded down to a double and used as exact cut points;  a window counts as not started on a
    sub-interval only if its gate is >= the sub-interval's upper end, certainly).

Usage: verify_tail_rig(prep, LR, par, Dv, h_hi) -> (ok, info) : S(eps) with the bounds D_j for EVERY real eps in (0, eps_of_h(h_hi)];  solve_D_rig(...) iterates the candidate bounds.
"""
import os
import sys

import numpy as np
from flint import arb, ctx

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T          # noqa: E402
import theory3c_rigorous as R3     # noqa: E402

ctx.prec = 100
A, up, lo, fup, fdown, sq = R3.A, R3.up, R3.lo, R3.fup, R3.fdown, R3.sq
clamp01, pos_up, logpos_up, bmax_up = R3.clamp01, R3.pos_up, R3.logpos_up, R3.bmax_up
Gam_r, psi_r, psi_sup_r, g2_sup_r = R3.Gam_r, R3.psi_r, R3.psi_sup_r, R3.g2_sup_r

VS = arb(1) / 40
SIG0 = arb(20)                       # a_max/(2 G v_s)
C_ENV = arb(100) / 3                 # 4/(G gamma^2 m)
S_MIN = arb(25)                      # the envelope is used for distances s >= 25 only


def sqn(z):
    """z^2 as a nonnegative ball (tight for a ball containing 0)"""
    return sq(z)


H_MARGIN = arb(1) + arb(2) ** -40                  # relative margin that makes eps_of_h(h_of_eps(eps)) >= eps decidable (the ball radii are about 2^-100 relative)


def min_lower(k1, k2):
    """the smaller of two lower bounds, selected by the LOWER ends (the ball with the smaller lower end is a lower bound of the minimum of the two quantities)"""
    return k1 if lo(k1) < lo(k2) else k2


def h_of_eps(eps):
    """an upper bound (point ball) of h = sqrt(eps/(2 v_s - eps)) for a float eps, enlarged by the relative margin 2^-40"""
    e = A(eps)
    return arb(up((e / (2 * VS - e)).sqrt() * H_MARGIN))


def eps_of_h(h):
    h2 = sqn(A(h))
    return 2 * VS * h2 / (1 + h2)


def g_atan(y, nterm=7):
    """atan(y)/y for a ball y >= 0 (may contain 0): 1 - y^2/3 + y^4/5 - ... with the first omitted term as the error bar if y < 1/2 (certain), else the monotone enclosure
    [atan(b)/b, atan(a)/a] with a = fdown(y) (1 if a <= 0), b = fup(y)"""
    if up(y) < arb(1) / 2:
        y2 = sqn(y)
        s_, pw = arb(0), arb(1)
        for n in range(nterm):
            s_ = s_ + (pw / (2 * n + 1) if n % 2 == 0 else -pw / (2 * n + 1))
            pw = pw * y2
        err = pw / (2 * nterm + 1)
        return s_ + (-err).union(err)
    a_f, b_f = fdown(y), fup(y)
    hi_ = arb(1) if a_f <= 0.0 else arb(a_f).atan() / arb(a_f)
    b_ = arb(b_f)
    return (b_.atan() / b_).union(hi_)


class UH:
    """The unit in the variable h = 1/kappa (a ball in [0, h_hi]):  q = 2 v_s (h^2 + tanh^2 x)/(1 + h^2),  w = (a_max/m) tanh x (eps + 2 v_s sinh^2 x),
    s(x) = sigma_0 (g(h/tanh x)/tanh x - x),  M = (G/m) q (q + 4 v_s sinh^2 x),  Psi = a_max^2 sinh^2 x/(2 G m).  L = 1/lambda."""

    def __init__(self, h):
        self.h = A(h).nonnegative_part()
        self.h2 = sqn(self.h)
        self.vs = VS
        self.G, self.m, self.amax = arb(1), arb(12), arb(1)
        self.eps = 2 * self.vs * self.h2 / (1 + self.h2)
        self.L = (4 * self.G * self.vs * self.h / (arb.pi() * self.amax)).nonnegative_part()
        self.muc = self.G * self.m * R3.SIGMA2_INV

    def q(self, x):
        return 2 * self.vs * (self.h2 + sqn(x.tanh())) / (1 + self.h2)

    def w(self, x):
        return (self.amax / self.m) * x.tanh() * (self.eps + 2 * self.vs * sqn(x.sinh()))

    def s(self, x):
        t = x.tanh()
        y = self.h / t
        return SIG0 * (g_atan(y) / t - x)

    def deta_dx(self, x):
        """an UPPER bound of d eta/dx = a_max/(G q): q >= 2 v_s tanh^2 x for every h in the ball (the eps = 0 value)"""
        return self.amax / (self.G * 2 * self.vs * sqn(x.tanh()))

    def M(self, x):
        q = self.q(x)
        return (self.G / self.m) * q * (q + 4 * self.vs * sqn(x.sinh()))

    def Psi(self, x):
        return self.amax ** 2 * sqn(x.sinh()) / (2 * self.G * self.m)

    _bisect = R3.RU._bisect

    def x_of_s(self, target):
        return self._bisect(self.s, A(target), True)

    def x_of_w(self, target):
        return self._bisect(self.w, A(target), False)


_RGS = {}


def rgs(LR, stage, c_lo, c_hi, nsub):
    """memoised R3.rig_rgs on the limit object LR (it does not depend on eps)"""
    key = (id(LR), stage, repr(c_lo), repr(c_hi), nsub)
    if key not in _RGS:
        _RGS[key] = (LR, R3.rig_rgs(LR['cs'], LR['As'], LR, stage, c_lo, c_hi, nsub))
    return _RGS[key][1]


def kernel_margin(LR, j, tau_end):
    """c_j = 1 - max over the dual breakpoints and the horizon of |u_j| (|u| is convex on each linear piece, the maximum is at an endpoint), as a lower bound (point ball)"""
    K = LR['K']
    taus = LR['tau_br']
    m = abs(LR['U_br'][0][j])
    for i in range(1, K + 1):
        m = bmax_up(m, abs(LR['U_br'][i][j]))
    uend = LR['U_br'][K][j] + LR['P'][K][j] * (tau_end - taus[K])
    m = bmax_up(m, abs(uend))
    return lo(arb(1) - up(m))


def dormant_pieces(LR, j, tau_end):
    """c(tau) = 1 - |u_j(tau)| of the non-entrant j as linear pieces (stage, tau_a, tau_b, c_a, c_b) on [0, tau_end] (balls), split where the dual CERTAINLY changes sign.  An unsplit piece
    whose dual might change sign is still a valid (lower) chord: c = 1 - |u| is concave on the piece."""
    K = LR['K']
    taus = LR['tau_br']
    out = []
    for i in range(K + 1):
        ta = taus[i]
        tb_ = taus[i + 1] if i < K else tau_end
        ua = LR['U_br'][i][j]
        ub = ua + LR['P'][i][j] * (tb_ - ta)
        if (lo(ua) > 0 and up(ub) < 0) or (up(ua) < 0 and lo(ub) > 0):
            f = abs(ua) / (abs(ua) + abs(ub))
            tc = ta + (tb_ - ta) * f
            out.append((i, ta, tc, 1 - abs(ua), arb(1)))
            out.append((i, tc, tb_, arb(1), 1 - abs(ub)))
        else:
            out.append((i, ta, tb_, 1 - abs(ua), 1 - abs(ub)))
    return out


def layer(prep, LR, par, hb, nq=1500):
    """Everything in S(eps) that does not involve the dormant kernels' mass or drift, for every h in the ball hb: the box, the layer constants, the stability constants, the weights and the
    entrants' pre-layer drift, from the ENCLOSED limit object LR.  Returns a dict, or dict(why=...) when a certain condition fails."""
    K = prep['K']
    cs, As, c0 = LR['cs'], LR['As'], LR['c0']
    J = len(c0)
    ev = LR['events']
    r, Om, d1, dc = A(par['r']), A(par['Om']), A(par['delta1']), A(par['dc'])
    dfrac, vsig = par['dfrac'], par['vsig']
    absA = np.abs(np.asarray(prep['As'], float))
    amin, amax = arb(float(absA.min())), arb(float(absA.max()))
    R = amax / amin
    if K > 1:
        Delta = None
        for i in range(K):
            for j in range(i + 1, K):
                d_ = abs(cs[i] - cs[j])
                Delta = d_ if Delta is None else Delta.min(d_)
    else:
        Delta = arb(0)
    Bx = R3.rig_box(K, R, Delta, r, Om, dfrac, vsig)
    if Bx is None:
        return dict(why='box')
    U = UH(hb)
    muc = U.muc
    vlo, vbar = amin * Bx['vlo'], amin * Bx['vbar']
    Ev, Ec, dt_, dh = Bx['Ev'], Bx['Ec'], Bx['dt'], Bx['dh']
    vsharp = vlo / 2
    x1 = U.x_of_s(d1)
    w1 = U.w(x1)
    xs = U.x_of_w(vsharp)
    gminus = vlo - w1
    if not (up(w1) < lo(vsharp) and lo(gminus) > 0):
        return dict(why='w1')
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
    xd = U.x_of_w(amin * Bx['vD_lo'])
    Mmin = U.M(lo(xd)).lower()
    k1 = Mmin * lo(Bx['bv'])
    k2 = muc * amin * amin * lo(Bx['bd'])
    kappa = min_lower(k1, k2)
    Vdiam = bmax_up(2 * dh, Bx['vD_hi'] - Bx['vD_lo'])
    WvD = amin * (Bx['vD_hi'] - Bx['vD_lo'])
    vbp = Bx['vbar']
    ell = Bx['ell']
    pre = U.Psi(x1) / gminus
    cF = (Bx['gl'] / Bx['bv'] + vbp * Bx['psl'] / Bx['bd']) if K > 1 else arb(0)
    cN = 1 / Bx['bv'] + arb(43) / 100 * vbp / Bx['bd']
    ent = [e['kernel'] for e in ev]
    feat = [e['feature'] for e in ev]
    vhat = [abs(e['w'][-1]) for e in ev]
    vhat_fin = []                                                # (6): the stage maximum of the old units' coefficients
    for i in range(K):
        m_ = None
        for e in ev:
            if len(e['w']) > i:
                v_ = abs(e['w'][i])
                m_ = v_ if m_ is None or up(v_) >= up(m_) else m_
        vhat_fin.append(m_)
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
    # the entrants' pre-layer drift (Lemma 5.3); the part x < 1e-9 through Psi * deta/dx <= a_max^3 cosh^2 x/(4 G^2 m v_s)
    x1h = max(fup(x1), 2e-9)
    edges2 = np.geomspace(1e-9, x1h, 3000)
    edges2[0], edges2[-1] = 1e-9, x1h
    pre_drift = []
    for k in range(K):
        f_ = feat[k]
        Hs, osc = rgs(LR, k, cs[f_] - u0, cs[f_] + u0, 200)
        if K > 1:
            hdev = (K - 1) * (psi_r(ell) * vbar + g2_sup_r(ell) * vbar * 2 * r) + arb(43) / 100 * amin * Om
        else:
            hdev = arb(43) / 100 * amin * Om
        He = Hs + hdev
        c_in = U.amax ** 3 * sqn(arb(1e-9).cosh()) / (4 * U.G ** 2 * U.m * U.vs)
        acc = up(Hs * muc * He / gminus ** 2 * c_in) * arb(1e-9)
        for a_, b_ in zip(edges2[:-1], edges2[1:]):
            xi = arb(float(a_)).union(arb(float(b_)))
            dcx = muc * He * U.Psi(xi) / gminus
            dcx = dcx if up(dcx) < up(2 * r) else 2 * r
            t_ = up(Hs * dcx)
            t_ = t_ if up(t_) < up(osc) else up(osc)
            acc = acc + up(t_ * U.deta_dx(xi) / gminus) * (arb(float(b_)) - arb(float(a_)))
        pre_drift.append(up(acc))
    return dict(why='', Bx=Bx, U=U, r=r, Om=Om, d1=d1, dc=dc, amin=amin, amax=amax, R=R, vlo=vlo, vbar=vbar, Ev=Ev, Ec=Ec, dt=dt_, dh=dh, vsharp=vsharp, x1=x1, w1=w1, xs=xs, gminus=gminus,
                T1=T1, Atr=Atr, T2=T2, T3=T3, Env=Env, Tlay=Tlay, kappa=kappa, Vdiam=Vdiam, WvD=WvD, pre=pre, cF=cF, cN=cN, ent=ent, feat=feat, vhat=vhat, vhat_fin=vhat_fin,
                wv=wv, wd=wd, pre_drift=pre_drift, ell=ell, par=par)


def verify_tail_rig(prep, LR, par, Dv, h_hi, nq=1500, tau_end=None):
    try:
        if not LR['ok']:
            return False, dict(why='limit object not certified: ' + LR['why'])
        return _verify_tail(prep, LR, par, Dv, h_hi, nq, tau_end)
    except (ValueError, ZeroDivisionError, ArithmeticError) as e_:
        return False, dict(why='%s: %s' % (type(e_).__name__, e_))


def _verify_tail(prep, LR, par, Dv, h_hi, nq, tau_end):
    """S(eps) for every real eps in (0, eps_hi], eps_hi = eps_of_h(h_hi), with the candidate bounds Dv (list over kernels, floats), on the enclosed limit object.  Returns (ok, info)."""
    K = prep['K']
    cs, As, c0 = LR['cs'], LR['As'], LR['c0']
    J = len(c0)
    ev = LR['events']
    info = {}
    hb = arb(0).union(A(h_hi))
    ly = layer(prep, LR, par, hb, nq)
    if ly['why']:
        return False, dict(why=ly['why'])
    U, Bx = ly['U'], ly['Bx']
    r, Om, d1, dc, amin, amax = ly['r'], ly['Om'], ly['d1'], ly['dc'], ly['amin'], ly['amax']
    vbar, dt_, dh, vsharp = ly['vbar'], ly['dt'], ly['dh'], ly['vsharp']
    T1, Atr, T2, T3, Env, Tlay, kappa, Vdiam, WvD, pre, cF, cN = (ly[k_] for k_ in ('T1', 'Atr', 'T2', 'T3', 'Env', 'Tlay', 'kappa', 'Vdiam', 'WvD', 'pre', 'cF', 'cN'))
    ent, feat, vhat, vhat_fin, wv, wd, pre_drift, w1 = (ly[k_] for k_ in ('ent', 'feat', 'vhat', 'vhat_fin', 'wv', 'wd', 'pre_drift', 'w1'))
    Lb, muc = U.L, U.muc                                       # Lb = 1/lambda: a ball [0, L_hi]
    taus = LR['tau_br']
    TE = arb(fup(arb(105) / 100 * ev[-1]['tau'] + 1)) if tau_end is None else A(tau_end)
    non = [j for j in range(J) if j not in ent]
    # dormant kernels through the eps-free envelope:  c(tau) = 1 - |u_j(tau)| >= cmin_j,  a_j = cmin_j - D_j L,  B(dist - D) <= C L^3 kap_j c(tau)^(-3),  kap_j = cmin^3/(a_j (a_j^2 - sigma_0^2 L^2)),
    # int_{linear piece} c^(-3) dtau = dtau (c_a + c_b)/(2 c_a^2 c_b^2) exactly;  mass_j = lambda int B dtau <= C L^2 kap_j Q_j.
    Dj = {j: A(Dv[j]) for j in non}
    cmin = {j: kernel_margin(LR, j, TE) for j in non}
    aj, kj, ej = {}, {}, {}
    for j in non:
        aj[j] = cmin[j] - Dj[j] * Lb
        if not (lo(aj[j]) > up((S_MIN + 1) * Lb)):
            return False, dict(why='dist<=D+25', j=j)
        a_lo, L_hi = lo(aj[j]), up(Lb)
        den = a_lo * (a_lo * a_lo - SIG0 * SIG0 * L_hi * L_hi)
        if not (lo(den) > 0):
            return False, dict(why='den<=0', j=j)
        kj[j] = up(cmin[j] ** 3 / den)
        ej[j] = up(C_ENV / den)                                  # peak = L^3 e_j
    pk = {j: up(Lb ** 3 * ej[j]) for j in non}
    spl = {j: dormant_pieces(LR, j, TE) for j in non}
    for j in non:
        for (_, ta, tb_, ca, cb) in spl[j]:
            if not (lo(ca) > 0 and lo(cb) > 0):
                return False, dict(why='margin<=0', j=j)
    Qj = {j: up(sum(((tb_ - ta) * (ca + cb) / (2 * ca * ca * cb * cb) for (_, ta, tb_, ca, cb) in spl[j]), arb(0))) for j in non}
    massL = {j: up(C_ENV * sqn(Lb) * kj[j] * Qj[j]) for j in non}
    mass_tot = up(sum((massL[j] for j in non), arb(0)))
    omega_pk = up(sum((pk[j] for j in non), arb(0)) / amin)
    omega_T = up(sum((TE * sqn(Lb) * ej[j] for j in non), arb(0)) / amin)       # omega_pk * Tfin = omega_pk tau_end / L
    Vtail = Vdiam / kappa + cF * K * pre / amin + cN * omega_T
    for nm_, v_ in [('mass_tot', mass_tot), ('omega_pk', omega_pk), ('omega_T', omega_T), ('Vtail', Vtail), ('Tlay', Tlay), ('pre', pre)] + [('pre_drift', x) for x in pre_drift]:
        if not v_.is_finite():
            return False, dict(why='not finite: %s' % nm_)
    # windows and unit L1 bounds
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
    # drift of the non-entrants:  |c_j - c_j^0| <= dcb_j(t) = mu_c ( W_j(t) + p_j (sum_{windows started} hd_i + sum_{not started} hpre_i + 0.43 mass_tot) ),  W = int H^(st(t')) B dt' (stage-resolved),
    # dual deviation <= int H^(st) dcb dt.  Gates: the EARLIEST possible start of window i is tau_{i+1} - Sm_i L_hi (a double, rounded down).
    drift = {}
    dcmax = arb(0)
    gate = [fdown(taus[i + 1] - up(Sm[i] * Lb)) for i in range(K)]
    for j in non:
        HO = [rgs(LR, k, c0[j] - dc, c0[j] + dc, 60) for k in range(K + 1)]
        Hst = [up(h_[0] + arb(43) / 100 * amin * omega_pk) for h_ in HO]
        hd, hpre = [], []
        for i in range(K):
            Xj = abs(c0[j] - cs[feat[i]]) - r - dc
            Xl = lo(Xj) if lo(Xj) > 0 else arb(0)
            hd.append(psi_sup_r(Xl) * Ttv[i] + g2_sup_r(Xl) * Ttd[i])
            hpre.append(psi_sup_r(Xl) * pre)
        Wq, dsum = arb(0), arb(0)
        for (st, ta, tb_, ca, cb) in spl[j]:
            dtp = tb_ - ta
            cuts = sorted(g_ for g_ in set(gate) if lo(ta) < g_ < up(tb_))
            bnds = [(arb(0), None)] + [(clamp01((arb(g_) - ta) / dtp), arb(g_)) for g_ in cuts] + [(arb(1), tb_)]
            for (f0, _), (f1, end_) in zip(bnds[:-1], bnds[1:]):
                c0_, c1_ = ca + (cb - ca) * f0, ca + (cb - ca) * f1
                Qs = (f1 - f0) * dtp * (c0_ + c1_) / (2 * c0_ * c0_ * c1_ * c1_)
                DM = up(Wq + Hst[st] * Qs)                                  # displacement integral at the end of the sub-interval (each earlier stage with its own H)
                eu = up(end_)
                hs = up(sum((hpre[i] if arb(gate[i]) >= eu else hd[i] for i in range(K)), arb(0)) + arb(43) / 100 * mass_tot)
                part1 = C_ENV * Lb * kj[j] * DM                                   # lambda * C L^2 kap DM
                part2 = sqn(Lb) * ej[j] * hs                                     # lambda * p_j * (...)
                dsum = dsum + up((f1 - f0) * dtp * Hst[st] * muc * (part1 + part2))
                dcb = up(muc * (C_ENV * sqn(Lb) * kj[j] * DM + Lb ** 3 * ej[j] * hs))
                dcmax = bmax_up(dcmax, dcb)
                Wq = DM
        drift[j] = up(dsum)
    for nm_, v_ in [('T1', T1), ('Atr', Atr), ('Tlay', Tlay), ('pre', pre), ('Env', Env), ('kappa', kappa), ('Vtail', Vtail), ('mass', mass_tot), ('omega', omega_pk)] + [('De', x) for x in De] + [('pre_drift', x) for x in pre_drift] + [('drift', drift[j]) for j in non] + [('cum', cum[j]) for j in non]:
        if not v_.is_finite():
            return False, dict(why='not finite: %s' % nm_)
    res = {}
    ok = True
    worst = 0.0
    for j in non:
        Etot = up(cum[j] + mass_tot + drift[j])
        res[j] = float(Etot)
        if not (Etot * (arb(10001) / 10000) < A(Dv[j])):            # strict, as the bootstrap of THEORY_T1.md section 6 needs
            ok = False
        worst = max(worst, float(Etot) / float(Dv[j]))
    cond = {}
    cond['fixed point E<D/1.0001'] = bool(ok)
    cond['Omega'] = bool(omega_pk <= Om)
    cond['dc'] = bool(up(dcmax) <= dc)
    p = [e['p_e'] for e in ev]
    # lambda (tau_{k+1} - tau_k) >= X  <=>  tau_{k+1} - tau_k >= X L   (L is a ball [0, L_hi])
    ordok = bool(lo(taus[1]) >= up(bmax_up((d1 + De[0]) / p[0], arb(0)) * Lb))
    for k in range(K - 1):
        lo_ = bmax_up((d1 + De[k + 1]) / p[k + 1], arb(0))
        hi_ = bmax_up((De[k] - d1) / p[k], arb(0)) + Tlay
        ordok = ordok and bool(lo(taus[k + 2] - taus[k + 1]) >= up((lo_ + hi_) * Lb))
    cond['order'] = bool(ordok)
    theta = (R3.THETA_REL_B * amax).min(R3.TENTH_B * amin)
    pkmax = None
    for j in non:
        pkmax = pk[j] if pkmax is None else bmax_up(pkmax, pk[j])
    cond['theta'] = bool(pkmax < theta and amin * Bx['vD_lo'] > theta) if pkmax is not None else bool(amin * Bx['vD_lo'] > theta)
    cond.update(R3.hypotheses_rig(LR, prep['cs'], np.asarray(prep['As'], float), Bx, r, dh, vsharp, amin, theta))
    info.update(cond=cond, Etot=res, worst_ratio=worst, T1=float(up(T1)), T2=float(up(T2)), T3=float(up(T3)), Tlay=float(up(Tlay)), De=[float(x) for x in De], mass=float(mass_tot),
                omega=float(omega_pk), omega_T=float(omega_T), Vtail=float(up(Vtail)), tau_end=float(TE.mid()), dcmax=float(up(dcmax)), kappa=float(lo(kappa)), Atr=float(up(Atr)), pre=float(up(pre)),
                prerad=float(pre.rad()), Lhi=float(up(Lb)), eps_hi=float(up(eps_of_h(hb))), pre_drift=[float(x) for x in pre_drift], drift_max=max([float(drift[j]) for j in non]),
                cum_max=max([float(up(cum[j])) for j in non]), Dv=[float(x) for x in Dv])
    return all(cond.values()), info


def solve_D_rig(prep, LR, par, h_hi, Efloat, fac=1.0003, ladder=(1.05, 1.15, 1.3, 1.6, 2.0), iters=6, **kw):
    """Candidate bounds D_j for the tail ball: start from ladder[i] x Efloat_j (the corrected float certificate's E_j at the ball's upper end), iterate D <- fac E(D) and stop at the first start for which
    `verify_tail_rig` accepts.  Only the final call decides: returns (ok, info, Dv), `ok` being the verdict of `verify_tail_rig` for the returned Dv (strict: E_j < D_j/1.0001)."""
    ent = {e['kernel'] for e in LR['events']}
    J = len(prep['c0'])
    last = (False, dict(why='no start'), [1.0] * J)
    for lf in ladder:
        Dv = [1.0 if j in ent else max(lf * Efloat[j], 1.0) for j in range(J)]
        for _ in range(iters):
            ok, info = verify_tail_rig(prep, LR, par, Dv, h_hi, **kw)
            last = (ok, info, Dv)
            if ok:
                return ok, info, Dv
            if 'Etot' not in info:
                break
            new = [Dv[j] if j in ent else float(fac * info['Etot'][j]) for j in range(J)]
            if all(abs(new[j] - Dv[j]) <= 1e-9 * Dv[j] for j in range(J)):
                break
            Dv = new
    return last
