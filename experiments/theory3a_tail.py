"""
Theorem T1, the small-eps tail: a ball-arithmetic proof that the certificate S(eps) holds for EVERY real eps in (0, eps_hi] with one set of per-kernel bounds D_j.

Why a tail argument is needed.  `theory_t1b_rigorous.verify` takes eps as a ball [a, b] with a > 0; the chain of such balls ends at 1e-24 and the unit's functions blow up
(kappa = sqrt((2 v_s - eps)/eps), wall lambda ~ eps^(-1/2)) for a ball containing 0.  Here the unit is parametrised by

    h = 1/kappa = sqrt(eps/(2 v_s - eps))        (eps = 2 v_s h^2/(1 + h^2)),        L = 1/lambda = 4 G v_s h/(pi a_max),

and h is a ball [0, h_hi].  Every function that enters the layer constants is then a smooth function of h^2 on the x-range where it is used (x >= x_1 > 0, where
y = h/tanh x is small), so the ball is a rigorous "sup over eps <= eps_hi" majorant; the dormant side (coefficient, mass, drift of the inactive kernels) is written in
L = 1/lambda with explicit powers (so that nothing is divided by a ball containing 0) and bounded by the eps-free envelope of Lemma 12.4,

    B_eps(s) <= env(s) := C/(s (s^2 - sigma_0^2)),   C = 4/(G gamma^2 m) = 100/3,   sigma_0 = a_max/(2 G v_s) = 20,   s >= 25,

which is Lemma 2.2 with the Adamovic-Mitrinovic inequality (re-proved in THEORY3A.md).  The dormant terms are bounded crudely and monotonically: the mass of kernel j by
(time horizon) x env(smallest distance), the peak by env(smallest distance), the drift by (horizon) x H x (centre-displacement bound at the largest mass/peak); the window-start
gating of `theory_t1b_rigorous.verify` is dropped (an upper bound).  All of these are O(L^2) or smaller, i.e. O(eps) or O(eps^(1/2)) in the sum E_j, against a layer term of order 10^2.

What is and is not rigorous: exactly as in `theory_t1b_rigorous` (limit-object data are floating-point inputs taken as exact, comparisons through `float(ball.upper())`).
"""
import os
import sys

import numpy as np
from flint import arb, ctx

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T          # noqa: E402
import theory_t1b_rigorous as RG   # noqa: E402

ctx.prec = 100
A, up, lo, sq = RG.A, RG.up, RG.lo, RG.sq
Gam_r, psi_r, psi_sup_r, g2_sup_r = RG.Gam_r, RG.psi_r, RG.psi_sup_r, RG.g2_sup_r
pos_up, logpos_up = RG.pos_up, RG.logpos_up


def sqn(z):
    """z^2 as a ball with a nonnegative lower end (arb balls built by `union` have a lower end a few ulps below 0, which turns a denominator 2 v_s (h^2 + tanh^2 x) into a ball containing 0)"""
    z = A(z)
    return (z * z).nonnegative_part()


def bmax_up(a, b):
    """the ball among a, b with the larger upper end (for upper bounds only); a NaN is an error, never silently replaced"""
    if a.is_nan() or b.is_nan():
        raise ValueError('nan in bmax_up')
    return a if float(a.upper()) >= float(b.upper()) else b


VS = arb(1) / 40
SIG0 = arb(20)                       # a_max/(2 G v_s)
C_ENV = arb(100) / 3                 # 4/(G gamma^2 m) = 2 v_s a_max sigma_0^3/m
S_MIN = arb(25)                      # the envelope is used for distances s >= 25 only


def h_of_eps(eps):
    """an upper bound (point ball) of h = sqrt(eps/(2 v_s - eps)) for a float eps"""
    e = A(eps)
    return arb((e / (2 * VS - e)).sqrt().upper())


def eps_of_h(h):
    h2 = sqn(A(h))
    return 2 * VS * h2 / (1 + h2)


def g_atan(y, nterm=7):
    """atan(y)/y for a ball y >= 0 (may contain 0): alternating series 1 - y^2/3 + y^4/5 - ... with the first omitted term as the error bar for y <= 1/2, atan(y)/y (decreasing in y) otherwise"""
    yl, yu = max(float(lo(y)), 0.0), float(up(y))
    if yu < 0.5:
        y2 = sqn(y)
        s_, pw = arb(0), arb(1)
        for n in range(nterm):
            s_ = s_ + (pw / (2 * n + 1) if n % 2 == 0 else -pw / (2 * n + 1))
            pw = pw * y2
        err = pw / (2 * nterm + 1)                      # first omitted term (nonneg ball): the truncation error is within +-err
        return s_ + (-err).union(err)
    a_ = arb(yl)
    hi_ = arb(1) if yl == 0.0 else a_.atan() / a_
    b_ = arb(yu)
    lo_ = b_.atan() / b_
    return lo_.union(hi_)


class UH:
    """The unit in the variable h = 1/kappa (a ball in [0, h_hi]):  q = 2 v_s (h^2 + tanh^2 x)/(1 + h^2),  w = (a_max/m) tanh x (eps + 2 v_s sinh^2 x),
    s(x) = sigma_0 (g(h/tanh x)/tanh x - x),  g(y) = atan(y)/y,  M = (G/m) q (q + 4 v_s sinh^2 x),  Psi = a_max^2 sinh^2 x/(2 G m).  L = 1/lambda."""

    def __init__(self, h):
        self.h = A(h).nonnegative_part()
        self.h2 = sqn(self.h)
        self.vs = VS
        self.G, self.m, self.amax = arb(1), arb(12), arb(1)
        self.eps = 2 * self.vs * self.h2 / (1 + self.h2)
        self.L = 4 * self.G * self.vs * self.h / (arb.pi() * self.amax)
        self.muc = self.G * self.m * RG.SIGMA2_INV

    def q(self, x):
        return 2 * self.vs * (self.h2 + sqn(x.tanh())) / (1 + self.h2)

    def w(self, x):
        return (self.amax / self.m) * x.tanh() * (self.eps + 2 * self.vs * sqn(x.sinh()))

    def s(self, x):
        t = x.tanh()
        y = self.h / t
        return SIG0 * (g_atan(y) / t - x)

    def deta_dx(self, x):
        """an UPPER bound of d eta/dx = a_max/(G q): q = 2 v_s (h^2 + tanh^2 x)/(1 + h^2) >= 2 v_s tanh^2 x (because (h^2 + t^2)/(1 + h^2) >= t^2 for t <= 1), i.e. the eps = 0 value; valid for every h in the ball
        and for x arbitrarily close to 0 (where the ball of q itself would contain 0)"""
        return self.amax / (self.G * 2 * self.vs * sqn(x.tanh()))

    def M(self, x):
        q = self.q(x)                                   # used only at x >= x_# > 0.1 for a LOWER bound of M (increasing in q)
        return (self.G / self.m) * q * (q + 4 * self.vs * sqn(x.sinh()))

    def Psi(self, x):
        return self.amax ** 2 * sqn(x.sinh()) / (2 * self.G * self.m)

    _bisect = RG.RU._bisect

    def x_of_s(self, target):
        return self._bisect(self.s, A(target), True)

    def x_of_w(self, target):
        return self._bisect(self.w, A(target), False)


def dormant_pieces(prep, j, tau_end):
    """the limit's c(tau) = 1 - |u_j(tau)| of the non-entrant j as linear pieces (stage, tau_a, tau_b, c_a, c_b) on [0, tau_end], split where the dual changes sign (floats taken as exact)"""
    Lm = prep['L']
    K = prep['K']
    taus = Lm['tau_br']
    out = []
    for i in range(K + 1):
        ta_f = float(taus[i])
        tb_f = float(taus[i + 1]) if i < K else float(tau_end)
        ta, tb_ = arb(ta_f), arb(tb_f)
        ua = arb(float(Lm['U_br'][i][j]))
        ub = ua + arb(float(Lm['P'][i][j])) * (tb_ - ta)
        if (float(lo(ua)) > 0 and float(up(ub)) < 0) or (float(up(ua)) < 0 and float(lo(ub)) > 0):
            f = abs(ua) / (abs(ua) + abs(ub))
            tc = ta + (tb_ - ta) * f
            out.append((i, ta, tc, 1 - abs(ua), arb(1)))
            out.append((i, tc, tb_, arb(1), 1 - abs(ub)))
        else:
            out.append((i, ta, tb_, 1 - abs(ua), 1 - abs(ub)))
    return out


_RGS = {}


def rgs(prep, stage, c_lo, c_hi, nsub):
    """memoised RG.rig_rgs (it does not depend on eps)"""
    k = (id(prep), stage, float(c_lo), float(c_hi), nsub)
    if k not in _RGS:
        _RGS[k] = RG.rig_rgs(prep, stage, c_lo, c_hi, nsub)
    return _RGS[k]


def kernel_margin(prep, j, tau_end):
    """c_j = min over the limit's dual breakpoints and the horizon of 1 - |u_j|  (|u| is convex on each linear piece, so its maximum is at an endpoint); floats taken as exact"""
    Lm = prep['L']
    taus = Lm['tau_br']
    K = prep['K']
    vals = [abs(float(Lm['U_br'][i][j])) for i in range(K + 1)]
    uend = float(Lm['U_br'][K][j]) + float(Lm['P'][K][j]) * (tau_end - float(taus[K]))
    vals.append(abs(uend))
    return arb(1) - arb(max(vals))


def layer(prep, par, hb, nq=1500):
    """Everything in S(eps) that does not involve the dormant kernels' mass or drift, for every h in the ball hb (a ball in [0, h_hi]): the box, the layer constants of Lemma 3.3, the stability
    constants of Lemma 4.1, the weights of Lemma 5.1 and the entrants' pre-layer drift of Lemma 5.3.  Returns a dict, or None with the reason under the key 'why'."""
    K, amin_f, L_ = prep['K'], prep['amin'], prep['L']
    cs, c0 = prep['cs'], prep['c0']
    J = len(c0)
    ev = L_['events']
    r_f, Om_f, delta1_f, dfrac, vsig, dc_f = par['r'], par['Om'], par['delta1'], par['dfrac'], par['vsig'], par['dc']
    r, Om, d1, dc = A(r_f), A(Om_f), A(delta1_f), A(dc_f)
    amin = A(amin_f)
    R = A(prep['R'])
    Delta = prep['Delta']
    Bx = RG.rig_box(K, R, Delta if np.isfinite(Delta) else 0.0, r, Om, dfrac, vsig)
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
    if not (float(up(w1)) < float(lo(vsharp)) and float(lo(gminus)) > 0):
        return dict(why='w1')
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
    mV = U.M(arb(float(lo(xs)))).lower()
    Y0 = pos_up(Vt - vsharp)
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
    feat = [e['feature'] for e in ev]
    vhat = [A(abs(float(e['w'][-1]))) for e in ev]
    vhat_fin = [A(abs(float(x_))) for x_ in ev[-1]['w']]
    wv = [[None] * K for _ in range(J)]
    wd = [[None] * K for _ in range(J)]
    for j in range(J):
        for i in range(K):
            if j in ent:
                k_ = ent.index(j)
                X = arb(0) if i == k_ else abs(arb(float(cs[feat[k_]])) - arb(float(cs[feat[i]]))) - 2 * r
            else:
                X = abs(arb(float(c0[j])) - arb(float(cs[feat[i]]))) - r - dc
            Xl = A(max(float(lo(X)), 0.0))
            wv[j][i], wd[j][i] = Gam_r(Xl), psi_sup_r(Xl)
    # the entrants' pre-layer drift (Lemma 5.3); the part x < 1e-9 is bounded through Psi * deta/dx <= a_max^3 cosh^2 x/(4 G^2 m v_s)
    x1h = float(up(x1))
    edges2 = np.geomspace(1e-9, max(x1h, 2e-9), 3000)
    pre_drift = []
    for k in range(K):
        f_ = feat[k]
        uu = float(up(u0))
        Hs, osc = rgs(prep, k, cs[f_] - uu, cs[f_] + uu, 200)
        if K > 1:
            hdev = (K - 1) * (psi_r(ell) * vbar + g2_sup_r(ell) * vbar * 2 * r) + arb(43) / 100 * amin * Om
        else:
            hdev = arb(43) / 100 * amin * Om
        He = Hs + hdev
        c_in = U.amax ** 3 * sqn(arb(1e-9).cosh()) / (4 * U.G ** 2 * U.m * U.vs)
        acc = up(Hs * muc * He / gminus ** 2 * c_in) * arb(1e-9)
        for a_, b_ in zip(edges2[:-1], edges2[1:]):
            xi = arb(a_).union(arb(b_))
            dcx = muc * He * U.Psi(xi) / gminus
            dcx = dcx if float(up(dcx)) < float(up(2 * r)) else 2 * r
            t_ = up(Hs * dcx)
            t_ = t_ if float(up(t_)) < float(up(osc)) else up(osc)
            acc = acc + up(t_ * U.deta_dx(xi) / gminus) * (arb(b_) - arb(a_))
        pre_drift.append(up(acc))
    return dict(why='', Bx=Bx, U=U, r=r, Om=Om, d1=d1, dc=dc, amin=amin, R=R, vlo=vlo, vbar=vbar, Ev=Ev, Ec=Ec, dt=dt_, dh=dh, vsharp=vsharp, x1=x1, w1=w1, xs=xs, gminus=gminus,
                T1=T1, Atr=Atr, T2=T2, T3=T3, Env=Env, Tlay=Tlay, kappa=kappa, Vdiam=Vdiam, WvD=WvD, pre=pre, cF=cF, cN=cN, ent=ent, feat=feat, vhat=vhat, vhat_fin=vhat_fin,
                wv=wv, wd=wd, pre_drift=pre_drift, ell=ell, par=par)


def verify_tail(prep, par, Dv, h_hi, nq=1500, tau_end=None):
    try:
        return _verify_tail(prep, par, Dv, h_hi, nq, tau_end)
    except ValueError as e_:
        return False, dict(why='ValueError: %s' % e_)


def _verify_tail(prep, par, Dv, h_hi, nq, tau_end):
    """S(eps) for every real eps in (0, eps_hi], eps_hi = eps_of_h(h_hi), with the candidate bounds Dv (list over kernels, floats).  Returns (ok, info)."""
    K, amin_f, L_ = prep['K'], prep['amin'], prep['L']
    cs, c0 = prep['cs'], prep['c0']
    J = len(c0)
    ev = L_['events']
    info = {}
    a = prep['a']
    hb = arb(0).union(A(h_hi))
    ly = layer(prep, par, hb, nq)
    if ly['why']:
        return False, dict(why=ly['why'])
    U, Bx = ly['U'], ly['Bx']
    r, Om, d1, dc, amin, R = ly['r'], ly['Om'], ly['d1'], ly['dc'], ly['amin'], ly['R']
    vbar, dt_, dh, vsharp = ly['vbar'], ly['dt'], ly['dh'], ly['vsharp']
    T1, Atr, T2, T3, Env, Tlay, kappa, Vdiam, WvD, pre, cF, cN = (ly[k_] for k_ in ('T1', 'Atr', 'T2', 'T3', 'Env', 'Tlay', 'kappa', 'Vdiam', 'WvD', 'pre', 'cF', 'cN'))
    ent, feat, vhat, vhat_fin, wv, wd, pre_drift, w1 = (ly[k_] for k_ in ('ent', 'feat', 'vhat', 'vhat_fin', 'wv', 'wd', 'pre_drift', 'w1'))
    dc_f, Om_f = par['dc'], par['Om']
    Lb, muc = U.L, U.muc                                       # Lb = 1/lambda: a ball [0, L_hi]
    tau_end = (1.05 * ev[-1]['tau'] + 1.0) if tau_end is None else tau_end
    taus = L_['tau_br']
    TE = arb(tau_end)
    non = [j for j in range(J) if j not in ent]
    # --- dormant kernels through the eps-free envelope.  With c(tau) = 1 - |u_j(tau)| (piecewise linear, >= cmin_j) and a_j = cmin_j - D_j L:
    #     B(dist - D) <= env(s) <= C L^3 kap_j c(tau)^(-3),  s = (c - D L)/L,  kap_j = cmin^3/(a_j (a_j^2 - sigma_0^2 L^2)) >= 1,
    # and  int_{linear piece} c^(-3) dtau = dtau (c_a + c_b)/(2 c_a^2 c_b^2)  exactly.  mass_j = lambda int B dtau <= C L^2 kap_j Q_j.
    Dj = {j: A(Dv[j]) for j in non}
    cmin = {j: kernel_margin(prep, j, tau_end) for j in non}
    aj, kj, ej = {}, {}, {}
    for j in non:
        aj[j] = cmin[j] - Dj[j] * Lb                             # = L (dist - D) at the closest approach
        if not float(lo(aj[j])) > float(up((S_MIN + 1) * Lb)):          # needs dist - D >= 26 > 25 for the envelope
            return False, dict(why='dist<=D+25', j=j)
        # den(a, L) = a (a^2 - sigma_0^2 L^2) increases with a and decreases with L: its smallest value on the ball is at (a_lo, L_hi) (endpoint evaluation, not midpoint-radius products)
        a_lo, L_hi = lo(aj[j]), up(Lb)
        den = a_lo * (a_lo * a_lo - SIG0 * SIG0 * L_hi * L_hi)
        if not float(lo(den)) > 0:
            return False, dict(why='den<=0', j=j)
        kj[j] = up(cmin[j] ** 3 / den)
        ej[j] = up(C_ENV / den)                                  # peak = L^3 e_j
    pk = {j: up(Lb ** 3 * ej[j]) for j in non}
    spl = {j: dormant_pieces(prep, j, tau_end) for j in non}      # list of (stage, tau_a, tau_b, c_a, c_b)
    Qj = {j: up(sum((dtau * (ca + cb) / (2 * ca * ca * cb * cb) for (_, ta, tb_, ca, cb) in spl[j] for dtau in [tb_ - ta]), arb(0))) for j in non}
    massL = {j: up(C_ENV * sqn(Lb) * kj[j] * Qj[j]) for j in non}
    mass_tot = up(sum((massL[j] for j in non), arb(0)))
    omega_pk = up(sum((pk[j] for j in non), arb(0)) / amin)
    omega_T = up(sum((TE * sqn(Lb) * ej[j] for j in non), arb(0)) / amin)       # omega_pk * Tfin = omega_pk tau_end / L
    Vtail = Vdiam / kappa + cF * K * pre / amin + cN * omega_T
    for nm_, v_ in [('mass_tot', mass_tot), ('omega_pk', omega_pk), ('omega_T', omega_T), ('Vtail', Vtail), ('Tlay', Tlay), ('pre', pre)] + [('pre_drift', x) for x in pre_drift] + [('cum', cum_) for cum_ in []]:
        if not v_.is_finite():
            return False, dict(why='not finite: %s' % nm_)
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
    # --- drift of the non-entrants:  |c_j - c_j^0| <= dcb_j(t) = mu_c (H^(st) I_j(t) + p_j (sum_{windows started} hd_i + 0.43 mass_tot)),  dual deviation <= int H^(st) dcb_j dt
    #     (stage-resolved H, I_j(t) <= C L^2 kap_j (cumulative int c^-3), no cap by the oscillation, windows gated by their earliest possible start tau_i - Sm_i L_hi)
    drift = {}
    dcmax = arb(0)
    gate = [float(taus[i + 1]) - float(up(Sm[i] * Lb)) for i in range(K)]
    for j in non:
        HO = [rgs(prep, k, c0[j] - dc_f, c0[j] + dc_f, 60) for k in range(K + 1)]
        Hst = [up(h_[0] + arb(43) / 100 * amin * omega_pk) for h_ in HO]
        hd = []
        for i in range(K):
            Xj = abs(arb(float(c0[j])) - arb(float(cs[feat[i]]))) - r - dc
            Xl = A(max(float(lo(Xj)), 0.0))
            hd.append(psi_sup_r(Xl) * Ttv[i] + g2_sup_r(Xl) * Ttd[i])
        Wq, dsum = arb(0), arb(0)                               # Wq = sum over the earlier sub-intervals of H^(stage) x (their cumulative-mass integral): the displacement integral int H^(st(t')) B dt' / (C L^2 kap)
        for (st, ta, tb_, ca, cb) in spl[j]:
            dtp = tb_ - ta
            fr = sorted([(float(((arb(g_) - ta) / dtp).mid()), g_) for g_ in gate if float(lo(ta)) < g_ < float(up(tb_))])
            frs = [arb(0)] + [(arb(g_) - ta) / dtp for _, g_ in fr] + [arb(1)]
            for f0, f1 in zip(frs[:-1], frs[1:]):
                q0 = float(up(ta + f0 * dtp))
                c0_, c1_ = ca + (cb - ca) * f0, ca + (cb - ca) * f1
                Qs = (f1 - f0) * dtp * (c0_ + c1_) / (2 * c0_ * c0_ * c1_ * c1_)
                DM = up(Wq + Hst[st] * Qs)                                  # displacement integral at the end of the sub-interval (each earlier stage with its own H)
                hs = up(sum((hd[i] for i in range(K) if gate[i] <= q0), arb(0)) + arb(43) / 100 * mass_tot)
                part1 = C_ENV * Lb * kj[j] * DM                                   # lambda * C L^2 kap DM   (= C L kap DM)
                part2 = sqn(Lb) * ej[j] * hs                                     # lambda * p_j * (...)       (= L^2 e_j (...))
                dsum = dsum + up((f1 - f0) * dtp * Hst[st] * muc * (part1 + part2))
                dcb = up(muc * (C_ENV * sqn(Lb) * kj[j] * DM + Lb ** 3 * ej[j] * hs))
                dcmax = bmax_up(dcmax, dcb)
                Wq = DM
        drift[j] = up(dsum)
    # --- totals and conditions
    for nm_, v_ in [('T1', T1), ('Atr', Atr), ('Tlay', Tlay), ('pre', pre), ('Env', Env), ('kappa', kappa), ('Vtail', Vtail), ('mass', mass_tot), ('omega', omega_pk)] + [('De', x) for x in De] + [('pre_drift', x) for x in pre_drift] + [('drift', drift[j]) for j in non] + [('cum', cum[j]) for j in non]:
        if not v_.is_finite():
            return False, dict(why='not finite: %s' % nm_)
    res = {}
    ok = True
    worst = 0.0
    for j in non:
        Etot = up(cum[j] + mass_tot + drift[j])
        res[j] = float(Etot)
        if not (float(Etot) < float(Dv[j]) / 1.0001):
            ok = False
        worst = max(worst, float(Etot) / float(Dv[j]))
    cond = {}
    cond['fixed point E<=D'] = ok
    cond['Omega'] = float(omega_pk) <= float(Om_f)
    cond['dc'] = float(up(dcmax)) <= float(dc_f)
    p = [arb(float(e['p_e'])) for e in ev]
    # lambda (tau_{k+1} - tau_k) >= X  <=>  tau_{k+1} - tau_k >= X L   (L is a ball [0, L_hi])
    ordok = float(lo(arb(float(taus[1])))) >= float(up(bmax_up((d1 + De[0]) / p[0], arb(0)) * Lb))
    for k in range(K - 1):
        lo_ = bmax_up((d1 + De[k + 1]) / p[k + 1], arb(0))
        hi_ = bmax_up((De[k] - d1) / p[k], arb(0)) + Tlay
        ordok = ordok and float(lo(arb(float(taus[k + 2])) - arb(float(taus[k + 1])))) >= float(up((lo_ + hi_) * Lb))
    cond['order'] = bool(ordok)
    theta = min(T.THETA_REL * float(a.max()), 0.1 * amin_f)
    cond['theta'] = bool(max(float(up(pk[j])) for j in non) < theta and float(lo(amin * Bx['vD_lo'])) > theta)
    cond.update(RG.hypotheses(prep, par, Bx, r, dh, vsharp))
    info.update(cond=cond, Etot=res, worst_ratio=worst, T1=float(up(T1)), T2=float(up(T2)), T3=float(up(T3)), Tlay=float(up(Tlay)), De=[float(x) for x in De], mass=float(mass_tot),
                omega=float(omega_pk), omega_T=float(omega_T), Vtail=float(up(Vtail)), tau_end=float(tau_end), dcmax=float(up(dcmax)), kappa=float(lo(kappa)), Atr=float(up(Atr)), pre=float(up(pre)), prerad=float(pre.rad()), Lhi=float(up(Lb)),
                eps_hi=float(up(eps_of_h(hb))), T1lo=float(lo(T1)), prelo=float(lo(pre)), Atrlo=float(lo(Atr)), pre_drift=[float(x) for x in pre_drift],
                drift_max=max([float(drift[j]) for j in non]), cum_max=max([float(up(cum[j])) for j in non]))
    return all(cond.values()), info



def solve_D(prep, par, h_hi, Efloat, fac=1.0003, ladder=(1.05, 1.15, 1.3, 1.6, 2.0), iters=6, **kw):
    """Candidate bounds D_j for the tail ball.  Start from ladder[i] * Efloat_j (the floating-point certificate's E_j at the ball's upper end), iterate D <- fac * E(D) (E depends on D
    only through the dormant terms, with a tiny slope) and stop at the first start for which `verify_tail` accepts.  Only the final call decides: returns (ok, info, Dv), `ok` being the
    verdict of `verify_tail` for the returned Dv (strict: E_j < D_j/1.0001 for every non-entrant kernel j)."""
    ent = {e['kernel'] for e in prep['L']['events']}
    J = len(prep['c0'])
    last = (False, dict(why='no start'), [1.0] * J)
    for lf in ladder:
        Dv = [1.0 if j in ent else max(lf * Efloat[j], 1.0) for j in range(J)]
        for _ in range(iters):
            ok, info = verify_tail(prep, par, Dv, h_hi, **kw)
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
