"""
The corrected certificate S(eps) of Theorem T1 (round theory3, direction N: merged, corrected verifications; supersedes the verifiers of theory_t1b_rigorous, theory3a and theory3c as history).

Findings of review/REVIEW_roundT3AC.md that are built in here (no existing file is edited: each corrected function is the SOURCE of the existing one with anchored replacements, each of which
must occur exactly once, compiled at import; `legacy` flags switch each repair off so that the tests can show that the old and the repaired quantities differ):

  (f) displacement accounting (Lemma 5.2).  For a dormant kernel l the centre drifts by  |c_l(t) - c_l^0| <= mu_c ( int_0^t H^(st(t')) B_l(t') dt' + peak_l * (window terms) ),  the weight H belonging to the
      stage at which the mass was accumulated.  The old code used  H^(st(t)) * int_0^t B  (the CURRENT stage's H on the whole past mass).  Now  W = sum over earlier pieces of H^(stage) x (their mass)
      + H^(st) x (mass so far in the current piece), dcb = mu_c (W + hinc).
  (g) window gating.  hinc counted the unit i only after its window has started (wstart_i).  Before that the entrant's pre-layer coefficient still enters d_c R_l with int v_i <= pre, so for every window not
      yet started the term  psi_sup(X_li) * pre  is added to the sum in hinc.
  (6) the old units' limit coefficient in a later window is the one of the CURRENT stage (Lemma 5.4), not the final block's: vhat_fin[i] = max over the stages s >= i of |w_i^(s)|.
  (5/8) kappa = min(kappa_1, kappa_2) is selected by the smaller LOWER end (the old code compared upper ends).
  (8) `inv_of`: the float gate is replaced by a certain ball comparison (lo(|d|) > 1e-6 (|x| + |y|)); both branches are valid bounds, so the gate never decided a verdict; now it uses no float of a ball.

Exports:  cert_n / cert2d_n / best_cert_n / best_cert2d_n  (floating-point certificate, n = 1 and planar),  verify_rig_n  (ball arithmetic, enclosed limit object),  build_* (with flags),
          cover_n  (chain of eps-balls with verify_rig_n).  The ball tail (0, eps_hi] is in theory3n_tail.py.
"""
import inspect
import os
import sys
import time

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T          # noqa: E402
import theory_t1b_wide as W        # noqa: E402,F401  (installs the cache of the residual-gradient sup in exp_theory_t1; speed only)
import theory_t1b_planar as PL     # noqa: E402
import theory3c_rigorous as R3     # noqa: E402
from flint import arb              # noqa: E402
from theory3n_tail import min_lower as _min_lower   # noqa: E402

# ----------------------------------------------------------------------------------------------------------------------------------
# anchored source replacements
# ----------------------------------------------------------------------------------------------------------------------------------
F_FLOAT = [("            Icum, dsum = 0.0, 0.0\n", "            Icum, dsum, Wprev = 0.0, 0.0, 0.0\n"),
           ("                    dcb = muc * (Hst[st] * Ipts + hinc)\n", "                    dcb = muc * (Wprev + Hst[st] * (Ipts - Icum) + hinc)\n"),
           ("                    Icum = float(Ipts[-1])\n", "                    Wprev = Wprev + Hst[st] * (float(Ipts[-1]) - Icum)\n                    Icum = float(Ipts[-1])\n")]
G_FLOAT = [("                    hinc = peak_j * (float(sum(hd[i] for i in range(K) if wstart[i] <= q0 + 1e-9)) + 0.43 * mass_ne)\n",
            "                    hinc = peak_j * (float(sum(hd[i] for i in range(K) if wstart[i] <= q0 + 1e-9)) + float(sum(psi_sup(Xj)[i] * pre_mass for i in range(K) if not (wstart[i] <= q0 + 1e-9))) + 0.43 * mass_ne)\n")]
V_FLOAT = [("    vhat_fin = np.abs(ev[-1]['w'])\n", "    vhat_fin = _stage_max_vhat(ev, K)\n")]

F_RIG = [("        Icum, dsum = arb(0), arb(0)\n", "        Icum, dsum, Wprev = arb(0), arb(0), arb(0)\n"),
         ("                    dcb = muc * (Hst[st] * Ig + hinc)\n", "                    dcb = muc * (Wprev + Hst[st] * (Ig - Icum) + hinc)\n"),
         ("                Icum = up(last)\n", "                Wprev = up(Wprev + Hst[st] * (up(last) - Icum))\n                Icum = up(last)\n")]
G_RIG = [("        hd = []\n", "        hd = []\n        hpre = []\n"),
         ("            hd.append(psi_sup_r(Xl) * Ttv[i] + g2_sup_r(Xl) * Ttd[i])\n", "            hd.append(psi_sup_r(Xl) * Ttv[i] + g2_sup_r(Xl) * Ttd[i])\n            hpre.append(psi_sup_r(Xl) * pre)\n"),
         ("                hinc = pk * (sum((hd[i] for i in range(nw)), arb(0)) + arb(43) / 100 * mass_tot)\n",
          "                hinc = pk * (sum((hd[i] for i in range(nw)), arb(0)) + sum((hpre[i] for i in range(nw, K)), arb(0)) + arb(43) / 100 * mass_tot)\n")]
V_RIG = [("    vhat_fin = [abs(x_) for x_ in ev[-1]['w']]\n", "    vhat_fin = _stage_max_vhat_rig(ev, K)\n")]
KAPPA_RIG = [("    kappa = k1 if up(k1) < up(k2) else k2\n", "    kappa = _min_lower(k1, k2)\n")]
INV_RIG = [("    def inv_of(x, y, dt_tau):\n        d_ = abs(x) - abs(y)\n        if (lo(d_) > 0 or up(d_) < 0) and abs(float(d_.mid())) > 1e-6 * max(abs(float(x.mid())), abs(float(y.mid())), 1e-300):\n            return dt_tau / d_\n        return None\n",
            "    inv_of = _inv_of\n")]

# the planar replacements of theory_t1b_planar._build_cert2 (its source is not recoverable: it is exec'd from a string)
PLANAR_REP = [("np.abs(c0[:, None] - cs[featarr][None, :])", "_nrm(c0[:, None, :] - cs[featarr][None, :, :])"),
              ("abs(cs[feat[k]] - cs[feat[i]])", "_nrm(cs[feat[k]] - cs[feat[i]])"),
              ("resid_grad_sup(prep, k, cs[f_] - u0, cs[f_] + u0)", "resid_grad_sup(prep, k, cs[f_], u0)"),
              ("resid_grad_sup(prep, k, c0[j] - dc, c0[j] + dc)", "resid_grad_sup(prep, k, c0[j], dc)"),
              ("np.abs(c0[j] - cs[featarr])", "_nrm(c0[j][None, :] - cs[featarr])")]


def _apply(src, edits):
    for a, b in edits:
        if src.count(a) != 1:
            raise RuntimeError('anchor found %d times: %r' % (src.count(a), a))
        src = src.replace(a, b)
    return src


def _stage_max_vhat(ev, K):
    """vhat_fin[i] = max over the stages s >= i of |w_i^(s)| (floats)"""
    return np.array([max(abs(float(e['w'][i])) for e in ev if len(e['w']) > i) for i in range(K)])


def _stage_max_vhat_rig(ev, K):
    """the same in ball arithmetic: a list of balls whose upper ends are the maxima of the upper ends"""
    out = []
    for i in range(K):
        m = None
        for e in ev:
            if len(e['w']) > i:
                v = abs(e['w'][i])
                m = v if m is None or R3.up(v) >= R3.up(m) else m
        out.append(m)
    return out


def _inv_of(x, y, dt_tau):
    """dt/ds = dtau/(|x| - |y|) of a linear piece of a dormant dual from |x| to |y|, or None (flat piece: the caller then uses peak x duration).  Both are valid bounds; the choice is by a CERTAIN comparison
    (the difference of the two absolute values is certainly nonzero by more than 1e-6 (|x| + |y|)); no float of a ball."""
    d_ = abs(x) - abs(y)
    tol = R3.up(abs(x) + abs(y)) / 1000000
    if R3.lo(d_) > tol or R3.up(d_) < -tol:
        return dt_tau / d_
    return None


def build_cert(legacy_f=False, legacy_g=False, legacy_v=False, planar=False):
    """The floating-point certificate `exp_theory_t1.cert` (or its planar version) with the repairs (f), (g), (6); legacy_* = True switches that repair off (all three True = the old function)."""
    src = inspect.getsource(T.cert)
    if not legacy_f:
        src = _apply(src, F_FLOAT)
    if not legacy_g:
        src = _apply(src, G_FLOAT)
    if not legacy_v:
        src = _apply(src, V_FLOAT)
    ns = dict(T.__dict__)
    ns['_stage_max_vhat'] = _stage_max_vhat
    if planar:
        src = _apply(src, PLANAR_REP)
        ns.update(_nrm=PL._nrm, resid_grad_sup=PL.rgs2_cached, G2_sup=PL.G2_sup2)
    exec(compile(src, 'theory3n_cert_planar' if planar else 'theory3n_cert', 'exec'), ns)
    return ns['cert']


def build_rig(legacy_f=False, legacy_g=False, legacy_v=False, legacy_kappa=False, legacy_inv=False):
    """`theory3c_rigorous.verify_rig` with the repairs; legacy_* = True switches one off (all True = the old function)."""
    src = inspect.getsource(R3.verify_rig)
    for flag, ed in ((legacy_f, F_RIG), (legacy_g, G_RIG), (legacy_v, V_RIG), (legacy_kappa, KAPPA_RIG), (legacy_inv, INV_RIG)):
        if not flag:
            src = _apply(src, ed)
    ns = dict(R3.__dict__)
    ns['_stage_max_vhat_rig'] = _stage_max_vhat_rig
    ns['_min_lower'] = _min_lower
    ns['_inv_of'] = _inv_of
    exec(compile(src, 'theory3n_verify_rig', 'exec'), ns)
    return ns['verify_rig']


cert_n = build_cert()
cert2d_n = build_cert(planar=True)
verify_rig_n = build_rig()


def best_cert_n(prep, eps, pars=None):
    """first parameter set of PARS for which the corrected float certificate holds (or the last failure): (out, index or -1)"""
    return _best(cert_n, prep, eps, pars)


def best_cert2d_n(prep, eps, pars=None):
    return _best(cert2d_n, prep, eps, pars)


def _best(cf, prep, eps, pars):
    pars = T.PARS if pars is None else pars
    o = None
    for i, par in enumerate(pars):
        o = cf(prep, eps, par)
        if o['ok']:
            return o, i
    return o, -1


# ----------------------------------------------------------------------------------------------------------------------------------
# the cover of [10^bot, 10^top] by eps-balls with the corrected verifier
# ----------------------------------------------------------------------------------------------------------------------------------
def cover_n(prep, LR, e_top, e_bot, wmin=1.0 / 512, budget=None, facs=(1.1, 1.5, 2.0, 4.0), verify=None):
    """Chain of eps-balls with outward-rounded ends on which `verify_rig_n` accepts the candidate bounds (D_j = fac x the corrected float certificate's E_j at the ball's upper end).
    Returns (list of (lo, hi), complete, lowest exponent reached, number of verify calls).  Same procedure as theory3c_rigorous.cover; only the verifier and the float certificate differ."""
    verify = verify_rig_n if verify is None else verify
    cur, segs, calls = e_top, [], 0
    t_end = None if budget is None else time.time() + budget
    wlast = 0.5
    ent = [e['kernel'] for e in LR['events']]
    while cur > e_bot + 1e-9:
        if t_end is not None and time.time() > t_end:
            return segs, False, cur, calls
        o, pi = best_cert_n(prep, 10.0 ** cur)
        if pi < 0:
            return segs, False, cur, calls
        par = T.PARS[pi]
        w, done = min(1.0, 2.0 * wlast), False
        while w >= wmin and not done:
            if t_end is not None and time.time() > t_end:
                return segs, False, cur, calls
            lo_ = max(cur - w, e_bot)
            ball = R3.eps_ball(lo_, cur)
            for fac in facs:
                Dv = [1.0 if j in ent else max(fac * o['Etot'][j], 1.0) for j in range(len(prep['c0']))]
                calls += 1
                try:
                    ok, info = verify(prep, LR, ball, par, Dv)
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
