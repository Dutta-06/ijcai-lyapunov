"""
Rigorous (ball-arithmetic) enclosure of the LIMIT OBJECT of Theorem T1 (round theory3, direction C).

The verifier of theory_t1b_rigorous.py took the limit object's data (block critical points, entry times, duals at the breakpoints, rates, margins, the entrant/feature assignment)
from floating-point code (`exp_theory_t1.limit_path`, root finder `hybr`).  This module re-derives all of it in python-flint arb arithmetic.  The target (centres c*_m,
amplitudes A_m, dictionary centres c0_j) is taken as the exact binary numbers that the doubles are.

  (1) Each block critical point y = (w, c) (k units, all K features in the residual) is enclosed by the Krawczyk operator
          K(X) = z0 - C F(z0) + (I - C J(X)) (X - z0),     C = mid(J(z0))^{-1},
      applied to the reduced system  F = (g, h) = 0  (g_i = <phi_i, residual>, h_i = <d phi_i, residual>; the critical equations w_i h_i = 0 of the code, with w_i != 0), whose Jacobian
      J is written out in closed form (checked against finite differences in the tests).  K(X) inside the interior of X proves: exactly one root of F lies in X (existence AND
      uniqueness in X).  The box is then contracted by X <- K(X) cap X (the root stays inside).  Floating-point (hybr + a few Newton steps in arb) is only used to PLACE the box.
  (2) The recursion of Theorem K is run on the enclosed blocks: rates p_j, hitting times dt_j = (1 - sgn(p_j) u_j)/|p_j|, the entrant = the unique minimiser (strict: the upper end of its
      time is below the lower end of every other time; a rate whose enclosure contains 0 only gets the lower bound (1 - |u_j|)/sup|p_j|), duals u_j at the breakpoints (balls),
      the entry times (balls), the margin |u_j| < 1 (strict) of Theorem K's (b), the nearest-unfitted-feature assignment (strict), the sign test sgn p_e = sgn A_f, the offsets d_0.
  (3) The data are returned as balls for `theory3c_rigorous.verify_rig`, which checks the hypotheses r > d_0, 'every block critical point lies in the relaxation set D-hat', theta < v_#.

Not part of this module: Theorem K's global-minimiser statement (a).  What T1 uses is the critical point of each block in D-hat (Lemma 4.1: the Gershgorin margins make it the unique one there),
and that is what is enclosed.
"""
import os
import sys

import numpy as np
from flint import arb, arb_mat, ctx

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))

LIMIT_PREC = 220                  # bits used inside the limit-object enclosure


class LimitFail(Exception):
    """A decision of the limit object could not be certified (reason in args[0])."""


def A(x):
    return x if isinstance(x, arb) else arb(x)


def sqb(z):
    """z^2 as a nonnegative ball (arb's z*z of a ball containing 0 can have a negative lower end)."""
    z = A(z)
    l_, u_ = z.lower(), z.upper()
    if l_ <= 0 and u_ >= 0:
        a_, b_ = (l_ * l_).upper(), (u_ * u_).upper()
        return arb(0).union(a_ if a_ >= b_ else b_).nonnegative_part()
    return z * z


def Gam(x):
    return (-sqb(x) / 4).exp()


def dGam(x):
    x = A(x)
    return -(x / 2) * Gam(x)


def d2Gam(x):
    x = A(x)
    return (sqb(x) / 4 - arb(1) / 2) * Gam(x)


def FJ(w, c, cs, As, jac=True):
    """Reduced critical system of the fit  F(y) = 1/2 || sum_m A_m phi*_m - sum_i w_i phi_i ||^2  over the units (w_i, c_i):
         g_i = sum_m A_m G(c_i - c*_m) - sum_j w_j G(c_i - c_j),    h_i = sum_m A_m G'(c_i - c*_m) - sum_{j != i} w_j G'(c_i - c_j)       (G = Gaussian Gram, G' its derivative)
    (the code's  Fg = (g, w h)).  Returns F = (g_1..g_k, h_1..h_k) and, if jac, the Jacobian in the variables (w_1..w_k, c_1..c_k) as an arb_mat.  w, c: lists of arb (balls allowed)."""
    k, M = len(w), len(cs)
    g, h = [], []
    for i in range(k):
        gi = arb(0)
        hi = arb(0)
        for m in range(M):
            x = c[i] - cs[m]
            gi = gi + As[m] * Gam(x)
            hi = hi + As[m] * dGam(x)
        for j in range(k):
            if j == i:
                gi = gi - w[i]                     # G(0) = 1
            else:
                x = c[i] - c[j]
                gi = gi - w[j] * Gam(x)
                hi = hi - w[j] * dGam(x)
        g.append(gi)
        h.append(hi)
    F = g + h
    if not jac:
        return F
    n = 2 * k
    ent = [[arb(0)] * n for _ in range(n)]
    for i in range(k):
        h2 = arb(0)
        for m in range(M):
            h2 = h2 + As[m] * d2Gam(c[i] - cs[m])
        for j in range(k):
            if j == i:
                ent[i][i] = arb(-1)                                     # dg_i/dw_i
                ent[i][k + i] = h[i]                                    # dg_i/dc_i = h_i
                ent[k + i][i] = arb(0)                                  # dh_i/dw_i
            else:
                x = c[i] - c[j]
                ent[i][j] = -Gam(x)                                     # dg_i/dw_j
                ent[i][k + j] = w[j] * dGam(x)                          # dg_i/dc_j
                ent[k + i][j] = -dGam(x)                                # dh_i/dw_j
                ent[k + i][k + j] = w[j] * d2Gam(x)                     # dh_i/dc_j
                h2 = h2 - w[j] * d2Gam(x)
        ent[k + i][k + i] = h2                                          # dh_i/dc_i
    Jm = arb_mat(n, n, [e for row in ent for e in row])
    return F, Jm


def _mat_of_col(v):
    return arb_mat(len(v), 1, list(v))


def krawczyk(X, cs, As):
    """One Krawczyk step on the box X (list of 2k balls, variables (w, c)).  Returns (Kvec, inside) with inside = K(X) lies in the INTERIOR of X (then F has exactly one root in X,
    and it lies in K(X))."""
    n = len(X)
    k = n // 2
    z0 = [x.mid() for x in X]
    F0 = FJ(z0[:k], z0[k:], cs, As, jac=False)
    _, J0 = FJ(z0[:k], z0[k:], cs, As)
    Cm = J0.mid().inv().mid()
    _, JX = FJ(X[:k], X[k:], cs, As)
    Id = arb_mat(n, n, [arb(1) if i == j else arb(0) for i in range(n) for j in range(n)])
    Mx = Id - Cm * JX
    d = _mat_of_col([X[i] - z0[i] for i in range(n)])
    Kv = _mat_of_col(z0) - Cm * _mat_of_col(F0) + Mx * d
    Kl = [Kv[i, 0] for i in range(n)]
    inside = all(Kl[i].lower() > X[i].lower() and Kl[i].upper() < X[i].upper() for i in range(n))
    return Kl, inside


def newton_polish(z, cs, As, iters=6):
    """Point Newton iteration in arb arithmetic (only places the box; it proves nothing)."""
    n = len(z)
    k = n // 2
    z = [x.mid() for x in z]
    for _ in range(iters):
        F, Jm = FJ(z[:k], z[k:], cs, As)
        dz = Jm.mid().solve(_mat_of_col([f.mid() for f in F]))
        z = [(z[i] - dz[i, 0]).mid() for i in range(n)]
    return z


def enclose_root(z_guess, cs, As, radii=(1e-14, 1e-11, 1e-9, 1e-7, 1e-5), contract=10, polish=True):
    """Certify and tighten a root of the reduced block system near z_guess (floats or arb).  Returns the list of 2k balls; raises LimitFail if no box passes the Krawczyk test."""
    z0 = [A(float(x)) if not isinstance(x, arb) else x for x in z_guess]
    if polish:
        z0 = newton_polish(z0, cs, As)
    last = None
    for rho in radii:
        X = [x.mid().union(x.mid() + arb(rho)).union(x.mid() - arb(rho)) for x in z0]
        Kl, ok = krawczyk(X, cs, As)
        last = (rho, ok)
        if ok:
            Xc = Kl
            for _ in range(contract):
                Kl2, _ = krawczyk(Xc, cs, As)
                Xn = []
                for a_, b_ in zip(Kl2, Xc):
                    try:
                        inter = a_.intersection(b_)
                    except Exception:
                        raise LimitFail('empty intersection in the contraction (cannot happen if the root is inside the box)')
                    Xn.append(inter)
                Xc = Xn
            return Xc
    raise LimitFail('Krawczyk test failed on all boxes (last radius %g, inside %s)' % last)


def exact_data(cs, As, c0):
    return [arb(float(x)) for x in cs], [arb(float(x)) for x in As], [arb(float(x)) for x in c0]


def sign_of(x):
    """+1 / -1 if the ball is certainly positive / negative, else 0"""
    if x.lower() > 0:
        return 1
    if x.upper() < 0:
        return -1
    return 0


def limit_rig(cs, As, c0, Lf=None, crit_block=None):
    """The limit object of Theorem K for the target (cs, As) and the dictionary c0, in ball arithmetic.  Lf: the floating-point `exp_theory_t1.limit_path` output (only used to place the boxes
    of the blocks and for the final cross-check; no decision is taken from it).  Returns a dict with balls:
        K, J, events [dict(k, kernel, feature, tau, p_e, w, c, d0, sign_ok)], tau_br [K+1], U_br [K+1][J], P [K+1][J], mu (min margin 1 - |u_j|, non-entrants), ok (bool), checks (dict), why.
    A decision that cannot be certified sets ok = False and records the reason (the structure is returned as far as it was built)."""
    import exp_theory_t1 as T
    with ctx.workprec(LIMIT_PREC):
        return _limit_rig(cs, As, c0, T, crit_block)


def _limit_rig(cs_f, As_f, c0_f, T, crit_block):
    cs, As, c0 = exact_data(cs_f, As_f, c0_f)
    K, J = len(cs), len(c0)
    u = [arb(0)] * J
    tau = arb(0)
    active, unfit = [], list(range(K))
    w, c = [], []
    tau_br, U_br, P, events = [tau], [list(u)], [], []
    chk = dict(krawczyk=True, strict_entrant=True, strict_feature=True, wall=True, sign=True)
    why = ''
    mu = arb(1)
    crit = crit_block or T.crit_block
    assigned = []

    def rates(w, c, active):
        p = []
        for j in range(J):
            if j in active:
                p.append(arb(0))
                continue
            pj = arb(0)
            for m in range(K):
                pj = pj + As[m] * Gam(c0[j] - cs[m])
            for i in range(len(w)):
                pj = pj - w[i] * Gam(c0[j] - c[i])
            p.append(pj)
        return p

    for k in range(K):
        p = rates(w, c, active)
        times, lows = {}, {}
        for j in range(J):
            if j in active:
                continue
            s = sign_of(p[j])
            if s != 0:
                times[j] = ((1 - s * u[j]) / (s * p[j]), s)
            else:
                pu = abs(p[j]).upper()
                den = (1 - abs(u[j])).lower()
                lows[j] = den / pu if pu > 0 else arb.pos_inf()
        if not times:
            why = 'stage %d: no kernel with a decided rate' % (k + 1)
            chk['strict_entrant'] = False
            break
        j1 = min(times, key=lambda j: float(times[j][0].mid()))
        dt1, s1 = times[j1]
        strict = dt1.lower() > 0
        for j, (dt, _) in times.items():
            if j != j1 and not (dt1.upper() < dt.lower()):
                strict = False
        for j, lb in lows.items():
            if not (dt1.upper() < lb):
                strict = False
        if not strict:
            why = 'stage %d: the entrant is not the strict argmin (kernel %d)' % (k + 1, j1)
            chk['strict_entrant'] = False
            break
        oth = [dt.lower() for j, (dt, _) in times.items() if j != j1] + [lb.lower() for lb in lows.values()]
        gap = float(((min(oth, key=lambda z: float(z)) if oth else arb(1e30)) - dt1.upper()) / dt1.upper()) if oth else float('inf')
        un = [u[j] + dt1 * p[j] for j in range(J)]
        un[j1] = arb(s1)
        for j in range(J):
            if j in active or j == j1:
                continue
            if not (abs(un[j]).upper() < 1):
                chk['wall'] = False
                why = why or 'stage %d: |u_j| < 1 not certified for kernel %d' % (k + 1, j)
            mu = mu.min((1 - abs(un[j])).lower())
        # nearest unfitted feature
        dist = {m: abs(c0[j1] - cs[m]) for m in unfit}
        mbest = min(dist, key=lambda m: float(dist[m].mid()))
        for m in unfit:
            if m != mbest and not (dist[mbest].upper() < dist[m].lower()):
                chk['strict_feature'] = False
                why = why or 'stage %d: nearest unfitted feature is not strict' % (k + 1)
        f = mbest
        sgn_ok = (s1 == (1 if float(As_f[f]) > 0 else -1))
        if not sgn_ok:
            chk['sign'] = False
            why = why or 'stage %d: sgn(p_e) != sgn(A_f)' % (k + 1)
        unfit.remove(f)
        active.append(j1)
        assigned.append(f)
        # block critical point (float guess, then Krawczyk)
        try:
            wf, cf = crit(np.asarray(cs_f, float), np.asarray(As_f, float), assigned)
            z = enclose_root(list(wf) + list(cf), cs, As)
            w, c = z[:len(assigned)], z[len(assigned):]
        except LimitFail as e:
            chk['krawczyk'] = False
            why = why or 'stage %d: %s' % (k + 1, e.args[0])
            break
        tau = tau + dt1
        events.append(dict(k=k + 1, kernel=j1, feature=f, tau=tau, p_e=abs(p[j1]), w=list(w), c=list(c), d0=abs(c0[j1] - cs[f]), sign_ok=bool(sgn_ok), gap=gap,
                           rad=max(float(x.rad()) for x in list(w) + list(c))))
        tau_br.append(tau)
        U_br.append(list(un))
        P.append(p)
        u = un
    ok = (len(events) == K) and all(chk.values())
    if len(events) == K:
        P.append(rates(w, c, active))
    return dict(K=K, J=J, events=events, tau_br=tau_br, U_br=U_br, P=P, mu=mu, ok=bool(ok), checks=chk, why=why, cs=cs, As=As, c0=c0)


def agrees_with_floats(LR, Lf):
    """largest discrepancy between the enclosure and the hybr/float limit object (a sanity check; it decides nothing)"""
    d = 0.0
    for e, ef in zip(LR['events'], Lf['events']):
        if e['kernel'] != ef['kernel'] or e['feature'] != ef['feature']:
            return float('inf')
        d = max(d, abs(float(e['tau'].mid()) - ef['tau']) / max(1.0, ef['tau']))
        for a_, b_ in zip(list(e['w']) + list(e['c']), list(ef['w']) + list(ef['c'])):
            d = max(d, abs(float(a_.mid()) - float(b_)))
    return d
