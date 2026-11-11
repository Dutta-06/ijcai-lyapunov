"""
Round theory3 / direction E: the certificate of Theorem T1 (exp_theory_t1.cert: the trained network forms exactly the limit object's units) recomputed for a GENERAL Gram, n = 1.

exp_theory_t1.py is not modified.  This module builds, for a given Gram, a private copy of it by exec of its source with the Gaussian-specific pieces replaced:
  * the elementary functions Gam, dGam, psi, psi_sup, G2_sup  -> the Gram's Gamma, Gamma', -Gamma', m_1, m_2 (tail envelopes, theory3e_kernels);
  * the constants 0.43 (sup psi) and 1/2 (sup |Gamma''|, the Lipschitz constant of the residual's gradient) -> M_1, M_2 of the Gram;
  * Gamma(r)/2 (the centre contraction rate in T_2) -> lambda_psi(r);   Gamma(ell), psi(ell), Gamma(X) (couplings and weights at distance >= ell, X) -> the tail envelopes m_0, m_1;
  * box_t1 (the box of Prop. A and the two Gershgorin margins of Lemma 4.1) -> `box_t1_gen`, the same formulas with the Gram's own numbers (row-sum cross-Hessian, as Lemma 4.1 needs).
Every other line (the unit of Lemma S1, the limit path, the layer durations T_1, T_2, T_3, the dual-deviation fixed point, the order conditions) is the repository's own, run unchanged.
The proof of Theorem T1 uses the kernel only through the replaced quantities (THEORY_T1.md sections 3-6; see NOTE.md of the round, section 9) -- but that this is so is a reading of the proof,
the replaced script is a *computation of the certificate for the Gram*, and its guarantees are those of exp_theory_t1.cert (double precision except the box constants, here ball arithmetic).
"""
import os
import sys
import types

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from flint import arb

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..')))
import theory3e_kernels as TK     # noqa: E402

SRC = os.path.join(HERE, 'exp_theory_t1.py')

_REPL = [
    ("sec = 0.5 * (np.abs(As).sum() + np.abs(w).sum())", "sec = M2G * (np.abs(As).sum() + np.abs(w).sum())"),
    ("kd = muc * vsharp * amin * float(Gam(r)) / 2", "kd = muc * vsharp * amin * LAMPSI_R(r)"),
    ("cF = (float(Gam(ell)) / Bx['beta_v'] + vbp * float(psi(ell)) / Bx['beta_d'])", "cF = (float(TAIL0(ell)) / Bx['beta_v'] + vbp * float(TAIL1(ell)) / Bx['beta_d'])"),
    ("wv, wd = Gam(X), psi_sup(X)", "wv, wd = TAIL0(X), psi_sup(X)"),
    ("hdev = ((K - 1) * (float(psi(ell)) * vbar", "hdev = ((K - 1) * (float(TAIL1(ell)) * vbar"),
    ("0.43", "M1G"),
]


def _vec(f):
    return np.vectorize(f, otypes=[float])


def box_t1_gen(kern, K, R, Delta, r, Omega, dhat_frac=0.05, varsigma=0.01):
    """The Gram version of exp_theory_t1.box_t1 (same keys): (L), (C), (N), and the Gershgorin margins beta_v = 1 - s_max - rowv, beta_d = t_min - s_max - rowd."""
    out = dict(ok=False)
    ell = Delta - 2 * r
    if K > 1 and ell <= TK.L_MIN:
        return out
    one = arb(1)
    Rr, Kk, Om = arb(R), arb(K - 1), arb(Omega)
    M1, M2 = arb(kern.M1()), arb(kern.M2(1))
    if K > 1:
        gl, psl, G2l = arb(TK.up(arb(kern.tail(0, ell)))), arb(TK.up(arb(kern.tail(1, ell)))), arb(TK.up(arb(kern.tail(2, ell))))
    else:
        gl = psl = G2l = arb(0)
    kap = Kk * gl
    if not (kap < 1):
        return out
    vbar = (Rr * (one + kap) + Om) / (one - kap)
    Ev = kap * (Rr + vbar) + Om
    Ec = Kk * (Rr + vbar) * psl + M1 * Om
    g_r, psi_r, lps_r = arb(kern.g_min(r)), arb(kern.psi_at(r)), arb(kern.lampsi(r))
    vlo = g_r - Ev
    Cm = psi_r - Ec
    mid = lambda x: float(x.mid())
    pos = lambda x: bool(x > 0)
    out.update(vbar=mid(vbar), Ev=mid(Ev), Ec=mid(Ec), vlo=mid(vlo), Cmargin=mid(Cm), kappa=mid(kap), psl=mid(psl), gl=mid(gl))
    if not (pos(vlo) and pos(Cm) and pos(lps_r)):
        return out
    dt = Ec / lps_r
    if not (dt < r):
        return out
    dh = dt + arb(dhat_frac) * (arb(r) - dt)
    dhu = TK.up(dh)
    vlo_d = arb(kern.g_min(dhu)) - Ev - arb(varsigma)
    if not pos(vlo_d):
        return out
    s_max = Rr * arb(kern.psihat(dhu))
    t_min = vlo_d * arb(kern.rho_min(dhu))
    lam = (one + t_min - TK._sqrt((one - t_min) * (one - t_min) + 4 * s_max * s_max)) / 2
    rowv = Kk * (gl + (Rr + vbar) * psl + vbar * psl) + M1 * Om
    rowd = Kk * ((Rr + vbar) * psl + vbar * (Rr + vbar) * G2l + vbar * psl + vbar * vbar * G2l) + M1 * Om + vbar * M2 * Om
    Hc = rowv if mid(rowv) > mid(rowd) else rowd
    nu = lam - Hc
    bv = one - s_max - rowv
    bd = t_min - s_max - rowd
    tmax = Rr * vbar * M2
    Hup = (one if mid(tmax) < 1 else tmax) + s_max + Hc
    out.update(dtilde=mid(dt), dhat=mid(dh), lam_own=mid(lam), Hc=mid(Hc), rowv=mid(rowv), rowd=mid(rowd), nu=mid(nu), Hup=mid(Hup), smax=mid(s_max), tmin=mid(t_min),
               beta_v=mid(bv), beta_d=mid(bd), vD_lo=mid(vlo_d), vD_hi=mid(Rr + Ev + arb(varsigma)), varsigma=varsigma, ok=bool(pos(nu) and pos(bv) and pos(bd)))
    return out


def build(kern):
    """A private copy of exp_theory_t1 for the Gram `kern` (n = 1).  Returns the module; use module.cert / best_cert / simulate / follows / prepare / make_target / PARS as in the original."""
    assert kern.n == 1
    src = open(SRC).read()
    for a, b in _REPL:
        assert a in src, a
        src = src.replace(a, b)
    mod = types.ModuleType('exp_theory_t1_gen_' + kern.name)
    mod.__file__ = SRC
    sys.modules[mod.__name__] = mod
    exec(compile(src, SRC, 'exec'), mod.__dict__)
    mod.Gam = lambda d: kern.fG(np.abs(np.asarray(d, float)))
    mod.dGam = lambda d: np.sign(np.asarray(d, float)) * kern.fdG(np.abs(np.asarray(d, float)))
    mod.psi = lambda x: -kern.fdG(np.asarray(x, float))
    mod.psi_sup = _vec(lambda x: kern.tail(1, float(x)))
    mod.G2_sup = _vec(lambda x: kern.tail(2, float(x)))
    mod.TAIL0 = _vec(lambda x: kern.tail(0, float(x)))
    mod.TAIL1 = _vec(lambda x: kern.tail(1, float(x)))
    mod.M1G = kern.M1()
    mod.M2G = kern.M2(1)
    mod.LAMPSI_R = lambda r: kern.lampsi(float(r))
    mod.box_t1 = lambda K, R, Delta, r, Omega, dhat_frac=0.05, varsigma=0.01: box_t1_gen(kern, K, R, Delta, r, Omega, dhat_frac, varsigma)
    mod._BOXC.clear()
    return mod
