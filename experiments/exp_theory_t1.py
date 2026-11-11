"""
Theorem 1 (the network follows the limit object): an explicit epsilon_0 in the controlled regime, and the numbers
====================================================================================================================================

Companion of THEORY_T1.md (the statements and proofs; lemma numbers below refer to it).  The setting is the one in which the entry layers
are controlled (THEORY_NOTE.md, Prop. A / Thm K): the continuum Gram <phi(c) phi(c')> = Z Gamma(c - c'), Gamma(d) = exp(-d^2/4) (lengths in
widths sigma), K Gaussian features of the kernels' width at separation >= Delta, noiseless data, equal masses eps b_j = eps, n = 1; the
trained network is the repository's P0 unit (Lemma S1 of the supplement) from alpha = 0, nu^2 = eps.  Time is t' = Z t (the factor Z only rescales time).

What the script computes (everything is a finite computation on the target; nothing here is fitted):

  A. The unit.  Lemma 2.1-2.3 numerically: the wall lambda(eps), the distance-to-wall s(x), the exact dormant coefficient B(s) against its explicit
     bound (the remainder of the core formula, Lemma 2.2), B(s) s^3 -> 4/(G gamma^2 m), the mass identity  int W d eta = a_max^2 sinh^2 x/(2 G m).
  B. The regime.  The separation Delta_T1(K, R, d_0) from which the box conditions (L), (C), (N) of Prop. A *and* the two Gershgorin margins of
     Lemma 4.1 hold (ball arithmetic when python-flint is present), next to Delta_entry of THEORY_NOTE.
  C. The certificate.  For each target (K = 1, 2, 3, continuum, random signs and offsets) and each eps of a half-decade grid, the sufficient
     condition `cert` of Theorem T1 (explicit constants: box, unit, limit path, layer durations T_1 T_2 T_3, window lengths, dual deviations D_j,
     dormant masses and centre drift, order of the entries) -> eps_0 = the largest grid value from which `cert` holds at every smaller one, next to
     the measured eps_max (the continuum network integrated in the eta-formulation: the formed set is the path's from eps_max downwards) and the
     number of violations (eps with cert true and the network not following the path: must be 0).
  D. The ingredients against the integrated network at eps_0 and eps_0/10: every inequality of the proof checked on the trajectory (dual deviations
     against D_j, entry shifts against S^-, S^+, the transit T_1, the relaxation T_2 + T_3, monotone entrants, g_e >= v_lo - v_e, centre offsets <= r,
     dormant mass <= Omega a_min, the final fit).
  E. The repository's network.  `factorized_rbf.train_flow` with the P0 readout (sampled inner product, 91 samples on [-1.5, 1.5], 31 kernels at
     spacing 0.1, the targets of Part I of exp_theory_h3 in the interior of the window) on the same eps grid: formed set against the continuum
     path, and `cert`'s eps_0 against the repository's own eps_max.
  F. The first-order formula  eps_* = (c_u min(q_1, q_2, q_3))^2  (q_1 = min_j mu_j/(D_j^inf + s_Omega), q_2, q_3 the order conditions),  c_u = lambda sqrt(eps) = pi a_max/(2 sqrt(2) G sqrt(v_s)),
     against eps_0 (an upper bound for it).
  G. eps_0 against the separation for two features with the local data fixed (Delta = 4.8 ... 11.5 widths).
  H. What sets the measured threshold: the actual dual deviation E_act of the dormant kernels against the margin lambda mu at eps_max and above it.

Usage: python experiments/exp_theory_t1.py --run [--quick] [--parts A,B,C,D,H,G,E,F] | --from-raw
       -> results/synthetic/theory_t1.md (+ theory_t1_raw.npz).  Environment: TT1_OUT (output directory).  Deterministic: every random draw is seeded;
       the report is regenerated from the raw alone (byte for byte).  --quick: a reduced run (2 targets, coarse grid), used by the tests.
Numerical caveat (stated in the report): the constants of Prop. A are evaluated in ball arithmetic (python-flint) when available; the rest is double
precision with relative margins of 0.1 % on the quadratures; the continuum network is integrated by LSODA (rtol 1e-10).
"""
import argparse
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

OUT = os.environ.get('TT1_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory_t1_raw.npz'), os.path.join(OUT, 'theory_t1.md')

try:
    from flint import arb
    HAVE_ARB = True
except Exception:                                    # python-flint missing: floats with a safety margin (the report says which)
    HAVE_ARB = False

# ----------------------------------------------------------------------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------------------------------------------------------------------
LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)       # the repository's law (exp_parametrization.LAW)
SIG = 0.12                                            # the kernels' width of the benchmark
H_LAT = 0.1 / SIG                                     # lattice spacing in widths
PSI_MAX = 0.5 * np.sqrt(2.0) * np.exp(-0.5)           # sup psi = psi(sqrt 2) = 0.4289 (< 0.43 of Prop. A)
THETA_REL = 0.01                                      # the formed threshold of the repository (factorized_rbf.THETA)
EX_GRID = np.arange(-3.0, -14.01, -0.5)               # log10 eps of the grid
EX_QUICK = np.arange(-4.0, -10.01, -1.0)
PARS = [dict(r=r, Om=Om, delta1=d1, dfrac=0.05, vsig=vs, dc=0.3)                       # the parameters of `cert`, tried in this order
        for r in (0.5, 0.6, 0.7, 0.8) for Om in (0.02, 0.05) for d1 in (-15.0, 0.0) for vs in (0.05,)]
C_TARGETS = [(1, 0.0, 1.0, 1), (1, 0.0, 1.0, 2), (1, 0.0, 1.0, 3), (1, 0.0, 1.0, 4),
             (2, 6.0, 1.0, 11), (2, 6.0, 1.0, 21), (2, 6.0, 2.4, 12), (2, 8.0, 1.0, 13), (2, 8.0, 1.0, 22), (2, 8.0, 2.4, 14),
             (2, 10.0, 1.0, 23), (2, 10.0, 2.4, 24), (3, 8.0, 1.0, 31), (3, 8.0, 1.5, 15), (3, 10.0, 2.0, 32)]
C_QUICK = [(1, 0.0, 1.0, 2), (2, 8.0, 1.0, 13)]
E_CELLS = [(2, 6.0, 2.4, 'random', 0), (2, 6.0, 2.4, 'random', 1), (2, 8.0, 2.4, 'random', 3), (2, 8.0, 2.4, 'random', 2),
           (2, 6.0, 1.0, 'random', 2), (2, 8.0, 1.0, 'random', 1), (2, 8.0, 2.4, 'worst', 0), (3, 8.0, 2.4, 'random', 0)]      # the cell (K, sep, R, mode, idx) of Part I's generator
E_QUICK = [(2, 8.0, 2.4, 'random', 3)]
E_SEED = 11000
B_RS, B_D0S = (1.0, 2.4, 10.0), (0.2, 0.4167, 0.6)


# ----------------------------------------------------------------------------------------------------------------------------------
# elementary functions
# ----------------------------------------------------------------------------------------------------------------------------------
def Gam(d):
    return np.exp(-np.asarray(d, float) ** 2 / 4.0)


def dGam(d):
    return -np.asarray(d, float) / 2.0 * np.exp(-np.asarray(d, float) ** 2 / 4.0)


def psi(x):
    return 0.5 * np.asarray(x, float) * np.exp(-np.asarray(x, float) ** 2 / 4.0)


def psi_sup(x0):
    """sup_{x >= x0} psi(x), x0 >= 0 (psi increases on [0, sqrt 2] and decreases beyond)."""
    x0 = np.maximum(np.asarray(x0, float), 0.0)
    return np.where(x0 <= np.sqrt(2.0), PSI_MAX, psi(x0))


def G2_sup(x0):
    """sup_{x >= x0} |Gamma''(x)|, |Gamma''(x)| = |x^2/4 - 1/2| exp(-x^2/4) (0.5 at 0, zero at sqrt 2, local maximum 0.2231 at sqrt 6)."""
    x0 = np.maximum(np.asarray(x0, float), 0.0)
    g = np.abs((x0 ** 2 / 4 - 0.5) * np.exp(-x0 ** 2 / 4))
    return np.where(x0 < np.sqrt(6.0), np.maximum(g, 0.2231), g)


# ----------------------------------------------------------------------------------------------------------------------------------
# the unit (Lemma 2.1-2.4): closed forms in t' units, equal masses eps (b = 1)
# ----------------------------------------------------------------------------------------------------------------------------------
class Unit:
    """w = (a_max/m) tanh x (eps + 2 v_s sinh^2 x), x = alpha/a_max; eta = (a_max/(2 G v_s)) (x + kappa arctan(kappa tanh x)), kappa^2 = (2 v_s - eps)/eps;
    wall lambda = pi a_max kappa/(4 G v_s); distance to the wall s(x) = lambda - eta(x) = (a_max kappa/(2 G v_s)) (arctan(1/(kappa tanh x)) - x/kappa);
    mobility M = dw/d eta; mass identity  int_0^eta W = a_max^2 sinh^2 x/(2 G m).  mu_c = G m/sigma^2 (centres in widths, time t' = Z t)."""

    def __init__(self, eps, law=LAW):
        self.m, self.G, self.gam, self.amax = law['m'], law['G'], law['gamma'], law['a_max']
        self.vs = self.gam * self.amax ** 2 / (4 * self.G)
        self.mu_c = self.G * self.m / SIG ** 2
        self.eps = float(eps)
        vs = self.vs
        self.kap = np.sqrt((2 * vs - eps) / eps)
        self.lam = np.pi * self.amax * self.kap / (4 * self.G * vs)
        xs = np.concatenate([[0.0], np.geomspace(1e-9, 16.0, 80000)])
        ss = self.s_of_x(xs)
        self._xs, self._ss_asc, self._xs_desc = xs, ss[::-1], xs[::-1]

    def w_of_x(self, x):
        x = np.asarray(x, float)
        return (self.amax / self.m) * np.tanh(x) * (self.eps + 2 * self.vs * np.sinh(x) ** 2)

    def s_of_x(self, x):
        x = np.asarray(x, float)
        k = self.kap
        th = np.tanh(x)
        th = np.where(th > 0, th, 1e-300)
        return (self.amax * k / (2 * self.G * self.vs)) * (np.arctan(1.0 / (k * th)) - x / k)

    def x_of_s(self, s):
        s = np.asarray(s, float)
        if np.any(s < self._ss_asc[0]):
            raise ValueError('distance below the table (x > 16)')
        return np.where(s >= self.lam, 0.0, np.interp(s, self._ss_asc, self._xs_desc))

    def B(self, s):
        """|w| of a unit at distance s from its wall (the sup of |w| over eta <= lambda - s)."""
        return self.w_of_x(self.x_of_s(s))

    def Psi_x(self, x):
        return self.amax ** 2 * np.sinh(np.asarray(x, float)) ** 2 / (2 * self.G * self.m)

    def Psi_s(self, s):
        return self.Psi_x(self.x_of_s(s))

    def deta_dx(self, x):
        e = self.eps
        return self.amax / (self.G * (e + (2 * self.vs - e) * np.tanh(x) ** 2))

    def M_of_x(self, x):
        e = self.eps
        t = np.tanh(x)
        return (self.G / self.m) * (e + (2 * self.vs - e) * t ** 2) * ((e + 2 * self.vs * np.sinh(x) ** 2) / np.cosh(x) ** 2 + 4 * self.vs * np.sinh(x) ** 2)

    def x_of_w(self, w):
        xs = self._xs
        return np.interp(w, self.w_of_x(xs), xs)

    def M_range(self, wlo, whi):
        """min and max of M(w) on [wlo, whi] (grid of 400 points; M increases with w)."""
        x = np.linspace(float(self.x_of_w(wlo)), float(self.x_of_w(whi)), 400)
        M = self.M_of_x(x)
        return float(M.min()), float(M.max())

    def core_bound(self, mu):
        """Lemma 2.2: for eta <= (1 - mu) lambda and s_e = cot(pi mu/2)/kappa < 1,
        |w| <= (a_max/m) s_e (eps + 2 v_s s_e^2/(1 - s_e^2)) <= rho^{3/2}/(1 - s_e^2) (a_max eps^{3/2}/(m sqrt(2 v_s))) cos(th)/sin^3(th), th = pi mu/2."""
        th = np.pi * np.asarray(mu, float) / 2
        se = 1.0 / (np.tan(th) * self.kap)
        ok = se < 1
        se = np.where(ok, se, 0.5)
        ex = (self.amax / self.m) * se * (self.eps + 2 * self.vs * se ** 2 / (1 - se ** 2))
        rho = 2 * self.vs / (2 * self.vs - self.eps)
        cl = rho ** 1.5 / (1 - se ** 2) * self.amax * self.eps ** 1.5 / (self.m * np.sqrt(2 * self.vs)) * np.cos(th) / np.sin(th) ** 3
        return np.where(ok, ex, np.inf), np.where(ok, cl, np.inf)


# ----------------------------------------------------------------------------------------------------------------------------------
# the box constants of Proposition A (normalisation a_min = 1) and the quantities of Lemmas 3.3, 4.1
# ----------------------------------------------------------------------------------------------------------------------------------
if HAVE_ARB:
    def _mk(x): return arb(x)
    def _ex(x): return x.exp()
    def _sq(x): return x.sqrt()
    def _mid(x): return float(x.mid())
    def _pos(x): return bool(x > 0)
else:
    def _mk(x): return float(x)
    def _ex(x): return float(np.exp(x))
    def _sq(x): return float(np.sqrt(x))
    def _mid(x): return float(x)
    def _pos(x): return bool(x > 1e-9)


def _e(x): return _ex(-(_mk(x) * _mk(x)) / 4)
def _ps(x): return _mk(x) / 2 * _e(x)
def _G2(x): return (_mk(x) * _mk(x) / 4 - _mk(1) / 2) * _e(x)


def box_t1(K, R, Delta, r, Omega, dhat_frac=0.05, varsigma=0.01):
    """(L), (C), (N) of Prop. A and the Gershgorin margins of Lemma 4.1 for K features of amplitudes in [1, R], separation >= Delta, box radius r,
    dormant budget Omega (units of a_min): ell = Delta - 2r >= sqrt 6, kappa_c = (K-1) Gamma(ell) < 1,
      vbar = (R(1+kappa_c) + Omega)/(1-kappa_c), E_v = kappa_c (R + vbar) + Omega, E_c = (K-1)(R + vbar) psi(ell) + 0.43 Omega,
      v_lo = Gamma(r) - E_v, dtilde = 2 E_c/Gamma(r), dhat = dtilde + dhat_frac (r - dtilde)     [D-hat = {|d| <= dhat, v in [Gamma(dhat) - E_v - varsigma, R + E_v + varsigma]}],
    rowv, rowd = the row sums of the cross-Hessian (the two terms whose maximum is H_c of THEORY_NOTE), the own block [[1, s], [s, t]] at dhat,
      beta_v = 1 - s_max - rowv,  beta_d = t_min - s_max - rowd   (> 0: Gershgorin margins of Lemma 4.1),  nu = lambda_own - H_c (> 0: (N)).
    `ok` = all inequalities certainly true (ball arithmetic) / true with a 1e-9 margin (floats)."""
    out = dict(ok=False)
    if K > 1:
        if r >= np.sqrt(2.0) or Delta - 2 * r < np.sqrt(6.0):
            return out
        ell = _mk(Delta - 2 * r)
    else:
        if r >= np.sqrt(2.0):
            return out
        ell = None
    rr, Rr, Kk, Om = _mk(r), _mk(R), _mk(K - 1), _mk(Omega)
    kap = Kk * _e(ell) if K > 1 else _mk(0)
    if not (kap < 1):
        return out
    vbar = (Rr * (1 + kap) + Om) / (1 - kap)
    Ev = kap * (Rr + vbar) + Om
    psl = _ps(ell) if K > 1 else _mk(0)
    G2l = _G2(ell) if K > 1 else _mk(0)
    gl = _e(ell) if K > 1 else _mk(0)
    Ec = Kk * (Rr + vbar) * psl + _mk(0.43) * Om
    vlo = _e(rr) - Ev
    Cm = _ps(rr) - Ec
    out.update(vbar=_mid(vbar), Ev=_mid(Ev), Ec=_mid(Ec), vlo=_mid(vlo), Cmargin=_mid(Cm), kappa=_mid(kap), psl=_mid(psl), gl=_mid(gl))
    if not (_pos(vlo) and _pos(Cm)):
        return out
    dt = 2 * Ec / _e(rr)
    if not (dt < rr):
        return out
    dh = dt + _mk(dhat_frac) * (rr - dt)
    vlo_d = _e(dh) - Ev - _mk(varsigma)
    if not _pos(vlo_d):
        return out
    s_max = Rr * _ps(dh)
    t_min = vlo_d * (_mk(1) / 2 - dh * dh / 4) * _e(dh)
    lam = (1 + t_min - _sq((1 - t_min) ** 2 + 4 * s_max * s_max)) / 2
    rowv = Kk * (gl + (Rr + vbar) * psl + vbar * psl) + _mk(0.43) * Om
    rowd = Kk * ((Rr + vbar) * psl + vbar * (Rr + vbar) * G2l + vbar * psl + vbar * vbar * G2l) + _mk(0.43) * Om + vbar * Om / 2
    Hc = rowv if _mid(rowv) > _mid(rowd) else rowd
    nu = lam - Hc
    bv = 1 - s_max - rowv
    bd = t_min - s_max - rowd
    tmax = Rr * vbar / 2
    Hup = (_mk(1) if _mid(tmax) < 1 else tmax) + s_max + Hc
    out.update(dtilde=_mid(dt), dhat=_mid(dh), lam_own=_mid(lam), Hc=_mid(Hc), rowv=_mid(rowv), rowd=_mid(rowd), nu=_mid(nu), Hup=_mid(Hup),
               smax=_mid(s_max), tmin=_mid(t_min), beta_v=_mid(bv), beta_d=_mid(bd), vD_lo=_mid(_e(dh) - Ev - _mk(varsigma)), vD_hi=_mid(Rr + Ev + _mk(varsigma)),
               varsigma=varsigma, ok=bool(_pos(nu) and _pos(bv) and _pos(bd)))
    return out


_BOXC = {}


def get_box(K, R, Delta, r, Om, dfrac, vsig):
    key = (K, round(R, 9), round(Delta, 9) if np.isfinite(Delta) else 0.0, round(r, 9), round(Om, 9), dfrac, vsig)
    if key not in _BOXC:
        _BOXC[key] = box_t1(K, R, Delta if np.isfinite(Delta) else 0.0, r, Om, dfrac, vsig)
    return _BOXC[key]


def delta_t1(K, R, d0, Om=0.02, step=0.1, lo=3.0, hi=24.0):
    """Smallest Delta (grid `step`) for which some r in (d0, 1.3] (grid 0.05) has box_t1 ok; nan if none."""
    rs = [x for x in np.arange(0.2, 1.31, 0.05) if x > d0 + 1e-9]
    for D in np.arange(lo, hi, step):
        for r in rs:
            if get_box(K, R, float(D), float(r), Om, 0.05, 0.01).get('ok'):
                return float(D)
    return float('nan')


# ----------------------------------------------------------------------------------------------------------------------------------
# the limit object (Theorem K recursion, continuum, Z = 1, walls 1) and the target's data
# ----------------------------------------------------------------------------------------------------------------------------------
def crit_block(cs, As, feats):
    """Critical point (w, c) of the fit over len(feats) units started at the features."""
    k = len(feats)

    def Fg(z):
        w, c = z[:k], z[k:]
        g = Gam(c[:, None] - cs[None, :]) @ As - Gam(c[:, None] - c[None, :]) @ w
        M = dGam(c[:, None] - c[None, :])
        np.fill_diagonal(M, 0.0)
        h = dGam(c[:, None] - cs[None, :]) @ As - M @ w
        return np.r_[g, w * h]
    s = root(Fg, np.r_[As[feats], cs[feats]], method='hybr', tol=1e-14)
    return s.x[:k], s.x[k:]


def limit_path(cs, As, c0):
    """The recursion of Theorem K: duals piecewise linear in tau (wall 1), entrant = argmin (1 - sgn(p) u)/|p|, block = critical point near the features.
    Returns the events (kernel, feature, tau, relative gap, offset d0, entrant's rate p_e on its last pre-entry interval, block w, c), the breakpoints
    tau_br, the duals U_br at the breakpoints and the rates P of each interval (active kernels have rate 0)."""
    K, J = len(cs), len(c0)
    u = np.zeros(J)
    tau = 0.0
    active, assigned = [], []
    w, c = np.zeros(0), np.zeros(0)
    ev, unfit = [], list(range(K))
    tau_br, U_br, P = [0.0], [u.copy()], []
    for k in range(K):
        p = Gam(c0[:, None] - cs[None, :]) @ As - (Gam(c0[:, None] - c[None, :]) @ w if len(w) else 0.0)
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
        f = min(unfit, key=lambda m: abs(c0[j1] - cs[m]))
        unfit.remove(f)
        active.append(j1)
        assigned.append(f)
        w, c = crit_block(cs, As, assigned)
        u[j1] = np.sign(p[j1])
        tau_br.append(tau)
        U_br.append(u.copy())
        ev.append(dict(k=k + 1, kernel=j1, feature=f, tau=tau, gap=float(gap), d0=float(abs(c0[j1] - cs[f])), p_e=float(abs(p[j1])),
                       sign_ok=bool(np.sign(p[j1]) == np.sign(As[f])), w=w.copy(), c=c.copy()))
    pend = Gam(c0[:, None] - cs[None, :]) @ As - Gam(c0[:, None] - c[None, :]) @ w
    P.append(np.where(np.isin(np.arange(J), active), 0.0, pend))
    return dict(events=ev, tau_br=np.array(tau_br), U_br=np.array(U_br), P=np.array(P), active=active, K=K, J=J)


def prepare(cs, As, c0):
    cs, As, c0 = np.asarray(cs, float), np.asarray(As, float), np.asarray(c0, float)
    K = len(cs)
    a = np.abs(As)
    amin = float(a.min())
    Delta = float(np.min(np.abs(cs[:, None] - cs[None, :]) + np.eye(K) * 1e9)) if K > 1 else np.inf
    L = limit_path(cs, As, c0)
    ent = set(L['active'])
    non = [j for j in range(L['J']) if j not in ent]
    mu = float(min(1 - np.abs(L['U_br'][:, non]).max(axis=1))) if non else np.inf
    return dict(cs=cs, As=As, c0=c0, K=K, a=a, amin=amin, R=float(a.max() / amin), Delta=Delta, L=L, mu_min=mu,
                gap_min=float(min(e['gap'] for e in L['events'])), d0=[e['d0'] for e in L['events']], sign_ok=all(e['sign_ok'] for e in L['events']))


def resid_grad_sup(prep, stage, c_lo, c_hi, n=401):
    """sup |ghat'| and the oscillation of ghat^(stage)(c) = sum_m A_m Gam(c - c*_m) - sum_{i<=stage} w_i Gam(c - chat_i) on [c_lo, c_hi] (grid + Lipschitz margin)."""
    cs, As, L = prep['cs'], prep['As'], prep['L']
    if stage == 0:
        w, cc = np.zeros(0), np.zeros(0)
    else:
        ev = L['events'][stage - 1]
        w, cc = ev['w'], ev['c']
    c = np.linspace(c_lo, c_hi, n)
    g = Gam(c[:, None] - cs[None, :]) @ As - (Gam(c[:, None] - cc[None, :]) @ w if len(w) else 0.0)
    gp = dGam(c[:, None] - cs[None, :]) @ As - (dGam(c[:, None] - cc[None, :]) @ w if len(w) else 0.0)
    step = (c_hi - c_lo) / (n - 1)
    sec = 0.5 * (np.abs(As).sum() + np.abs(w).sum())
    return float(np.abs(gp).max() + 0.5 * step * sec), float(g.max() - g.min() + step * float(np.abs(gp).max()))


def kernel_pieces(prep, j, lam, Tfin):
    """The limit's distance-to-wall of the non-entrant kernel j as pieces (t_a, t_b, dist_a, dist_b), linear, split where the dual changes sign."""
    L = prep['L']
    tb = lam * L['tau_br']
    ts = np.r_[tb, Tfin]
    out = []
    for i in range(L['K'] + 1):
        ta_, tb_ = ts[i], ts[i + 1]
        ua = L['U_br'][i][j] * lam
        ub = ua + L['P'][i][j] * (tb_ - ta_)
        if ua * ub < 0:
            tc = ta_ + (tb_ - ta_) * abs(ua) / (abs(ua) + abs(ub))
            out.append((ta_, tc, lam - abs(ua), lam))
            out.append((tc, tb_, lam, lam - abs(ub)))
        else:
            out.append((ta_, tb_, lam - abs(ua), lam - abs(ub)))
    return out


_GRID = np.unique(np.r_[np.linspace(0, 1, 120), np.geomspace(1e-10, 0.5, 120), 1 - np.geomspace(1e-10, 0.5, 120)])


def kernel_mass(U, pieces, D):
    """Lemma 5.2: mass = int B(dist - D) dt over the pieces (exact through the mass identity), peak = max B.  None if dist - D <= 1 somewhere."""
    mass, peak = 0.0, 0.0
    for (ta, tb, da, db) in pieces:
        sa, sb = da - D, db - D
        if min(sa, sb) <= 1.0:
            return None
        dt = tb - ta
        if dt <= 0:
            continue
        if abs(sb - sa) <= 1e-9 * max(sa, sb):
            Bv = float(U.B(sa))
            peak = max(peak, Bv)
            mass += Bv * dt
        else:
            peak = max(peak, float(U.B(sa)), float(U.B(sb)))
            mass += abs(float(U.Psi_s(sb)) - float(U.Psi_s(sa))) * dt / abs(sb - sa)
    return mass, peak


# ----------------------------------------------------------------------------------------------------------------------------------
# the certificate of Theorem T1
# ----------------------------------------------------------------------------------------------------------------------------------
def cert(prep, eps, par, Tfin_tau=None):
    """The sufficient condition of Theorem T1 at eps for the parameters par = (r, Om, delta1, dfrac, vsig, dc).  Returns a dict: ok, why, D (largest dual
    deviation bound over the non-entrant kernels, in eta units), Dv (per kernel), De (the entrants' pre-layer deviation bounds), the window lengths,
    Sp/Sm (shift bounds), T1 T2 T3, kappa, and the elementary conditions.  Notation as in THEORY_T1.md section 6."""
    K, amin, a, R, Delta, L = prep['K'], prep['amin'], prep['a'], prep['R'], prep['Delta'], prep['L']
    cs, c0 = prep['cs'], prep['c0']
    J = len(c0)
    ev = L['events']
    r, Om, delta1, dfrac, vsig, dc = par['r'], par['Om'], par['delta1'], par['dfrac'], par['vsig'], par['dc']
    out = dict(ok=False, why='')
    if not prep['sign_ok']:
        out['why'] = 'sign'
        return out
    if r <= max(prep['d0']) + 1e-9:
        out['why'] = 'r<=d0'
        return out
    Bx = get_box(K, R, Delta, r, Om, dfrac, vsig)
    if not Bx.get('ok'):
        out['why'] = 'box'
        return out
    U = Unit(eps)
    lam, muc = U.lam, U.mu_c
    tau_end = (1.05 * ev[-1]['tau'] + 1.0) if Tfin_tau is None else Tfin_tau
    Tfin = lam * tau_end
    tk = lam * L['tau_br'][1:]
    vlo, vbar = amin * Bx['vlo'], amin * Bx['vbar']
    Ev, Ec, dt_, dh = Bx['Ev'], Bx['Ec'], Bx['dtilde'], Bx['dhat']
    vsharp = vlo / 2
    # ---- Lemma 3.3: the layer (transit T1, centre T2, coefficient T3)
    x1 = float(U.x_of_s(delta1))
    w1 = float(U.w_of_x(x1))
    xs = float(U.x_of_w(vsharp))
    if w1 >= vsharp:
        out['why'] = 'w1>=vsharp'
        return out
    gminus = vlo - w1
    xg = np.linspace(x1, xs, 8001)
    f = U.deta_dx(xg) / (vlo - U.w_of_x(xg))
    T1 = float(np.trapezoid(f, xg)) * 1.001
    Atr = float(np.trapezoid(U.w_of_x(xg) * f, xg)) * 1.001
    u0 = max(max(prep['d0']), dt_)
    kd = muc * vsharp * amin * float(Gam(r)) / 2
    T2 = max(0.0, np.log(max((u0 - dt_) / (dh - dt_), 1.0))) / kd
    Vt = amin * (float(Gam(dh)) - Ev)
    vup = amin * (R + Ev)
    mV, _ = U.M_range(vsharp, max(Vt, vsharp * 1.0001))
    Y0 = max(Vt - vsharp, 0.0)
    T3lo = np.log(max(Y0 / (amin * vsig), 1.0)) / mV
    if vbar > vup * 1.0001:
        mV2, _ = U.M_range(vup, vbar)
        Y0p = vbar - vup
        T3hi, envhi = np.log(max(Y0p / (amin * vsig), 1.0)) / mV2, Y0p / mV2
    else:
        T3hi = envhi = 0.0
    T3 = max(T3lo, T3hi)
    Wv0 = amin * (R + 2 * Ev - float(Gam(dh)))
    Env = Y0 / mV + envhi + Wv0 * T3
    Tlay = T1 + T2 + T3
    # ---- Lemma 4.1: the weighted-max-norm stability constants
    vD_lo, vD_hi = amin * Bx['vD_lo'], amin * Bx['vD_hi']
    Mmin, Mmax = U.M_range(vD_lo, vD_hi)
    kappa = min(Mmin * Bx['beta_v'], muc * amin ** 2 * Bx['beta_d'])
    Vdiam = max(2 * dh, Bx['vD_hi'] - Bx['vD_lo'])
    WvD = vD_hi - vD_lo
    vbp = Bx['vbar']
    ell = Delta - 2 * r if K > 1 else np.inf
    pre_mass = float(U.Psi_x(x1)) / gminus * 1.001
    cF = (float(Gam(ell)) / Bx['beta_v'] + vbp * float(psi(ell)) / Bx['beta_d']) if K > 1 else 0.0
    cN = 1.0 / Bx['beta_v'] + 0.43 * vbp / Bx['beta_d']
    out.update(Bx=Bx, lam=lam, tk=tk, Tfin=Tfin, T1=T1, T2=T2, T3=T3, Tlay=Tlay, Atr=Atr, Env=Env, x1=x1, w1=w1, gminus=gminus, kappa=kappa, Vdiam=Vdiam,
               pre_mass=pre_mass, mV=mV, Mmin=Mmin, Mmax=Mmax, vlo=vlo, vbar=vbar, vsharp=vsharp)
    # ---- kernels and weights (Lemma 5.1)
    ent = [e['kernel'] for e in ev]
    non = [j for j in range(J) if j not in ent]
    feat = [e['feature'] for e in ev]
    featarr = np.array(feat)
    X = np.maximum(0.0, np.abs(c0[:, None] - cs[featarr][None, :]) - r - dc)
    for k in range(K):
        for i in range(K):
            X[ent[k], i] = 0.0 if i == k else max(0.0, abs(cs[feat[k]] - cs[feat[i]]) - 2 * r)
    wv, wd = Gam(X), psi_sup(X)
    vhat = np.array([abs(e['w'][-1]) for e in ev])
    vhat_fin = np.abs(ev[-1]['w'])
    # ---- Lemma 5.3: the entrants' pre-layer centre drift
    pre_drift = []
    for k in range(K):
        f_ = ev[k]['feature']
        Hs, osc = resid_grad_sup(prep, k, cs[f_] - u0, cs[f_] + u0)
        hdev = ((K - 1) * (float(psi(ell)) * vbar + float(G2_sup(ell)) * vbar * 2 * r) if K > 1 else 0.0) + 0.43 * amin * Om
        He = Hs + hdev
        xg2 = np.concatenate([[0.0], np.geomspace(1e-9, max(x1, 2e-9), 20000)])
        dcx = np.minimum(2 * r, muc * He * U.Psi_x(xg2) / gminus)
        pre_drift.append(float(np.trapezoid(np.minimum(osc, Hs * dcx) * U.deta_dx(xg2) / gminus, xg2)) * 1.001)
    # ---- the fixed point for the per-kernel deviation bounds D_j
    Dv = np.full(J, 1.0)
    tt = np.r_[0.0, tk, Tfin]
    for it in range(120):
        mass_ne, fail, Bpk = 0.0, False, 0.0
        sumB = np.zeros(len(tt))
        per = []
        for j in non:
            pcs = kernel_pieces(prep, j, lam, Tfin)
            res = kernel_mass(U, pcs, Dv[j])
            if res is None:
                fail = True
                break
            mass_j, peak_j = res
            mass_ne += mass_j
            Bpk = max(Bpk, peak_j)
            per.append((j, pcs, peak_j))
            for ti, t in enumerate(tt):
                for (ta, tb_, da, db) in pcs:
                    if ta <= t <= tb_ * (1 + 1e-12) + 1e-9:
                        sumB[ti] += float(U.B((da + (db - da) * ((t - ta) / (tb_ - ta) if tb_ > ta else 0.0)) - Dv[j]))
                        break
        if fail:
            out.update(why='dist<=D', D=float(Dv[non].max()))
            return out
        omega_peak = float(sumB.max()) / amin
        Vtail = Vdiam / kappa + cF * K * pre_mass / amin + cN * omega_peak * Tfin
        cum = np.zeros(J)
        c_own, c_owntail, c_oldw, c_oldtail = np.zeros(J), np.zeros(J), np.zeros(J), np.zeros(J)
        De, Wk, Sp, Sm, Tv_l, Td_l = [], [], [], [], [], []
        for k in range(K):
            j = ent[k]
            De_k = cum[j] + sum(wv[j, i] * pre_mass for i in range(k, K)) + mass_ne + pre_drift[k]
            De.append(De_k)
            pk = ev[k]['p_e']
            Sp_k, Sm_k = max(0.0, (De_k - delta1) / pk), max(0.0, (De_k + delta1) / pk)
            Sp.append(Sp_k)
            Sm.append(Sm_k)
            Wlen = (abs(delta1) + De_k) / pk + Tlay
            Wk.append(Wlen)
            Tv = pre_mass + Atr + vhat[k] * T1 + vbar * max(Sp_k, Sm_k) + vbar * T2 + Env
            Td = (r + dt_) * Atr + vhat[k] * (r + dt_) * T2 + vhat[k] * 2 * dh * T3 + (r + dt_) * w1 * Sp_k
            Tv_l.append(Tv + amin * Vtail)
            Td_l.append(Td + vhat[k] * Vtail)
            c_own = c_own + wv[:, k] * Tv + wd[:, k] * Td
            c_owntail = c_owntail + (wv[:, k] * amin + wd[:, k] * vhat[k]) * Vtail
            for i in range(k):
                c_oldw = c_oldw + (wv[:, i] * WvD + wd[:, i] * vhat_fin[i] * 2 * dh) * Wlen
                c_oldtail = c_oldtail + (wv[:, i] * amin + wd[:, i] * vhat_fin[i]) * Vtail
            cum = c_own + c_owntail + c_oldw + c_oldtail
        # total L1 deviations of unit i (own window and tails, later windows): Ttv[i] >= int |v_i - vhat_i| dt, Ttd[i] >= int min(v_i, vhat_i)|d_i - dhat_i| dt
        Ttv = [Tv_l[i] + sum(WvD * Wk[k] + amin * Vtail for k in range(i + 1, K)) for i in range(K)]
        Ttd = [Td_l[i] + sum(vhat_fin[i] * 2 * dh * Wk[k] + vhat_fin[i] * Vtail for k in range(i + 1, K)) for i in range(K)]
        # centre drift of the non-entrants: |c_j - c_j^0| <= dcb_j(t) = mu_c (int H B dt + peak (sum_{windows started} psi_sup Ttv + G2_sup Ttd + 0.43 mass)); dual deviation <= int min(osc, H dcb) dt
        drift = np.zeros(J)
        dcmax = 0.0
        tstage = np.r_[0.0, tk, Tfin]
        wstart = [tk[i] - Sm[i] for i in range(K)]
        for (j, pcs, peak_j) in per:
            HO = [resid_grad_sup(prep, k, c0[j] - dc, c0[j] + dc) for k in range(K + 1)]
            Hst = [h_[0] + 0.43 * amin * omega_peak for h_ in HO]
            Ost = [h_[1] + 0.43 * amin * omega_peak * 2 * dc for h_ in HO]
            Xj = np.maximum(0.0, np.abs(c0[j] - cs[featarr]) - r - dc)
            hd = psi_sup(Xj) * np.array(Ttv) + G2_sup(Xj) * np.array(Ttd)
            Icum, dsum = 0.0, 0.0
            for (ta, tb_, da, db) in pcs:
                st = min(max(int(np.searchsorted(tstage, ta * (1 + 1e-12) + 1e-12, side='right') - 1), 0), K)
                sa, sb = da - Dv[j], db - Dv[j]
                dt = tb_ - ta
                if dt <= 0:
                    continue
                cuts = sorted(set([ta, tb_] + [w_ for w_ in wstart if ta < w_ < tb_]))
                for q0, q1 in zip(cuts[:-1], cuts[1:]):
                    f0, f1 = (q0 - ta) / dt, (q1 - ta) / dt
                    hinc = peak_j * (float(sum(hd[i] for i in range(K) if wstart[i] <= q0 + 1e-9)) + 0.43 * mass_ne)
                    sA, sB = sa + (sb - sa) * f0, sa + (sb - sa) * f1
                    dtt = q1 - q0
                    if abs(sB - sA) <= 1e-9 * max(sA, sB):
                        Ipts = Icum + float(U.B(sA)) * dtt * _GRID
                    else:
                        Ipts = Icum + np.abs(U.Psi_s(sA + (sB - sA) * _GRID) - float(U.Psi_s(sA))) * dtt / abs(sB - sA)
                    dcb = muc * (Hst[st] * Ipts + hinc)
                    dsum += float(np.sum(np.minimum(Ost[st], Hst[st] * dcb[1:]) * (_GRID[1:] - _GRID[:-1]) * dtt))
                    Icum = float(Ipts[-1])
                    dcmax = max(dcmax, float(dcb[-1]))
            drift[j] = dsum
        Etot = cum + mass_ne + drift
        D = float(Dv[non].max()) if non else 0.0
        jst = max(non, key=lambda j: Etot[j]) if non else 0
        out.update(comp=np.array([c_own[jst], c_owntail[jst], c_oldw[jst], c_oldtail[jst], mass_ne, drift[jst]]), D=D, Dnew=float(max(Etot[j] for j in non)) if non else 0.0, Dv=Dv.copy(), De=De, Wk=Wk, Sp=Sp, Sm=Sm, mass_ne=mass_ne, omega_peak=omega_peak,
                   dcmax=dcmax, Vtail=Vtail, Etot=Etot, drift=drift, cum=cum, Bpk=Bpk)
        if all(Etot[j] <= Dv[j] for j in non):
            break
        Dv = np.where(np.isin(np.arange(J), non), np.maximum(Etot, Dv) * 1.0001, Dv)
    else:
        out['why'] = 'D diverges'
        return out
    # ---- the elementary conditions
    p = [e['p_e'] for e in ev]
    theta = min(THETA_REL * float(a.max()), 0.1 * amin)
    ordok = bool(tk[0] >= max(0.0, (delta1 + De[0]) / p[0]))
    for k in range(K - 1):
        lo = max(0.0, (delta1 + De[k + 1]) / p[k + 1])
        hi = max(0.0, (De[k] - delta1) / p[k]) + Tlay
        ordok &= bool(tk[k + 1] - tk[k] >= lo + hi)
    conds = dict(Omega=bool(out['omega_peak'] <= Om), dc=bool(out['dcmax'] <= dc), order=ordok, theta=bool(out['Bpk'] < theta and vD_lo > theta))
    out['conds'] = conds
    out['ok'] = all(conds.values())
    out['why'] = '' if out['ok'] else 'conds'
    return out


def best_cert(prep, eps):
    """The first parameter set of PARS for which `cert` holds (or the last failure)."""
    o = None
    for i, par in enumerate(PARS):
        o = cert(prep, eps, par)
        if o['ok']:
            return o, i
    return o, -1


# ----------------------------------------------------------------------------------------------------------------------------------
# targets and the integrated networks
# ----------------------------------------------------------------------------------------------------------------------------------
def make_target(K, sep, R, seed):
    """Continuum target: K features, separations about sep widths (random offsets +- h/2 from the lattice), amplitudes (1, R) or (1, sqrt R, R) times 0.7, random order and
    signs; the dictionary is the lattice of spacing 0.1/sigma over the features' hull +- 8 widths, at a random phase."""
    rng = np.random.default_rng(seed)
    cs = np.arange(K) * sep + rng.uniform(-H_LAT / 2, H_LAT / 2, K)
    cs = cs - cs.mean()
    a = np.array([1.0]) if K == 1 else (np.array([1.0, R]) if K == 2 else np.array([1.0, np.sqrt(R), R]))
    a = a[rng.permutation(K)] * 0.7
    As = a * rng.choice([-1.0, 1.0], K)
    c0 = np.arange(cs.min() - 8, cs.max() + 8, H_LAT) + rng.uniform(0, H_LAT)
    return cs, As, c0


def _wfun(x, eps, U):
    return (U.amax / U.m) * np.tanh(x) * (eps + 2 * U.vs * np.sinh(x) ** 2)


def simulate(prep, eps, Tfin_tau=None, rtol=1e-10, nlocal=2500, nglob=4000):
    """The continuum network in the x-formulation (dx/dt = (G/a_max)(eps + (2 v_s - eps) tanh^2 x) g, dc/dt = mu_c w h) from x = 0, c = c0, to T_fin, recorded on a
    grid refined around the limit's entry times."""
    cs, As, c0 = prep['cs'], prep['As'], prep['c0']
    J = len(c0)
    L = prep['L']
    ev = L['events']
    U = Unit(eps)
    lam, muc = U.lam, U.mu_c
    tau_end = (1.05 * ev[-1]['tau'] + 1.0) if Tfin_tau is None else Tfin_tau
    Tfin = lam * tau_end
    tk = lam * L['tau_br'][1:]

    def f(t, z):
        x, c = z[:J], z[J:]
        w = _wfun(x, eps, U)
        g = Gam(c[:, None] - cs[None, :]) @ As - Gam(c[:, None] - c[None, :]) @ w
        M = dGam(c[:, None] - c[None, :])
        np.fill_diagonal(M, 0.0)
        h = dGam(c[:, None] - cs[None, :]) @ As - M @ w
        return np.r_[(U.G / U.amax) * (eps + (2 * U.vs - eps) * np.tanh(x) ** 2) * g, muc * w * h]
    loc = [np.linspace(max(0, t - 600), min(Tfin, t + 900), nlocal) for t in tk]
    ts = np.unique(np.concatenate([np.linspace(0, Tfin, nglob)] + loc))
    sol = solve_ivp(f, (0, Tfin), np.r_[np.zeros(J), c0], method='LSODA', rtol=rtol, atol=np.r_[np.full(J, 1e-14), np.full(J, 1e-12)], t_eval=ts)
    return sol, U


def follows(prep, eps, sol, U):
    """Does the integrated network form exactly the limit's entrants (|w| > theta at the end) with the exact fit?  Returns (formed_ok, exact)."""
    As, cs, c0 = prep['As'], prep['cs'], prep['c0']
    J = len(c0)
    ev = prep['L']['events']
    ent = sorted(e['kernel'] for e in ev)
    x, c = sol.y[:J, -1], sol.y[J:, -1]
    w = _wfun(x, eps, U)
    thr = min(THETA_REL * np.abs(As).max(), 0.1 * np.abs(As).min())
    formed = sorted(np.flatnonzero(np.abs(w) > thr).tolist())
    ok = (formed == ent) and sol.status == 0
    ex = bool(ok and all(abs(abs(w[e['kernel']]) - abs(As[e['feature']])) < 1e-3 * abs(As[e['feature']]) and abs(c[e['kernel']] - cs[e['feature']]) < 0.05 for e in ev))
    return bool(ok), ex


# ----------------------------------------------------------------------------------------------------------------------------------
# A. the unit
# ----------------------------------------------------------------------------------------------------------------------------------
def part_A():
    rows = []
    for eps in (1e-4, 1e-8, 1e-12):
        U = Unit(eps)
        # mass identity by quadrature in x: int W d eta vs a^2 sinh^2 x/(2 G m)
        x = np.linspace(0, 3.0, 200001)
        eta = U.lam - U.s_of_x(np.maximum(x, 1e-12))
        eta[0] = 0.0
        q = np.trapezoid(U.w_of_x(x) * np.gradient(eta, x), x)
        ident = abs(q - U.Psi_x(3.0)) / U.Psi_x(3.0)
        # the core bound: B(s) vs the explicit bound, at mu = 0.02 .. 0.5 (s = mu lambda)
        for mu in (0.02, 0.05, 0.1, 0.3, 0.5):
            s = mu * U.lam
            b = float(U.B(s))
            ex, cl = U.core_bound(mu)
            rows.append([eps, U.lam, U.lam * np.sqrt(eps), mu, b, float(ex), float(cl), b * s ** 3, ident])
    return dict(A_rows=np.array(rows))


# ----------------------------------------------------------------------------------------------------------------------------------
# B. the regime: Delta_T1
# ----------------------------------------------------------------------------------------------------------------------------------
def part_B():
    rows = []
    for K in (2, 3):
        for R in B_RS:
            for d0 in B_D0S:
                rows.append([K, R, d0, delta_t1(K, R, d0, 0.0), delta_t1(K, R, d0, 0.02)])
    return dict(B_rows=np.array(rows))


# ----------------------------------------------------------------------------------------------------------------------------------
# C. the certificate against the integrated network
# ----------------------------------------------------------------------------------------------------------------------------------
def thr_first(ok, ex):
    """The largest grid value (first index of the descending grid) from which every smaller one is True; nan if none."""
    for i in range(len(ex)):
        if all(ok[i:]):
            return float(ex[i])
    return float('nan')


def scan_target(prep, ex_grid, runner=None):
    """cert and the measured following for every eps of the grid (and D / parameter index when cert holds)."""
    okc, okm, okx, Dd, Pi = [], [], [], [], []
    for ex in ex_grid:
        eps = 10.0 ** ex
        o, pi = best_cert(prep, eps)
        okc.append(bool(o['ok']))
        Dd.append(o['D'] if o['ok'] else np.nan)
        Pi.append(pi)
        sol, U = simulate(prep, eps, nlocal=300)
        m, e = follows(prep, eps, sol, U)
        okm.append(m)
        okx.append(e)
    return np.array(okc, float), np.array(okm, float), np.array(okx, float), np.array(Dd), np.array(Pi, float)


def part_C(quick=False):
    warnings.filterwarnings('ignore')
    tg = C_QUICK if quick else C_TARGETS
    ex = EX_QUICK if quick else EX_GRID
    meta, OKC, OKM, OKX, DD, PI, COMP, FINE = [], [], [], [], [], [], [], []
    for (K, sep, R, seed) in tg:
        cs, As, c0 = make_target(K, sep, R, seed)
        prep = prepare(cs, As, c0)
        okc, okm, okx, Dd, Pi = scan_target(prep, ex)
        e0 = thr_first(okc, ex)
        comp = np.full(10, np.nan)
        if np.isfinite(e0):
            o0, pi0 = best_cert(prep, 10.0 ** e0)
            lam0 = o0['lam']
            comp = np.r_[o0['comp'], o0['D'], max(max(o0['Sp']), max(o0['Sm'])) / lam0, o0['T1'] + o0['T2'] + o0['T3'], float(pi0)]
        COMP.append(comp)
        fine = [float(best_cert(prep, 10.0 ** (e0 - dd))[1] >= 0) if np.isfinite(e0) else np.nan for dd in (0.125, 0.25, 0.375)]
        FINE.append(fine)
        meta.append([K, sep, R, seed, prep['Delta'] if K > 1 else np.nan, prep['R'], prep['mu_min'], prep['gap_min'], max(prep['d0']), len(c0),
                     min(e['p_e'] for e in prep['L']['events']), e0, thr_first(okm, ex)])
        OKC.append(okc)
        OKM.append(okm)
        OKX.append(okx)
        DD.append(Dd)
        PI.append(Pi)
    return dict(C_ex=ex, C_meta=np.array(meta, float), C_okc=np.array(OKC), C_okm=np.array(OKM), C_okx=np.array(OKX), C_D=np.array(DD), C_pi=np.array(PI), C_comp=np.array(COMP), C_fine=np.array(FINE))


# ----------------------------------------------------------------------------------------------------------------------------------
# D. the ingredients against the trajectory
# ----------------------------------------------------------------------------------------------------------------------------------
D_COLS = ('K', 'sep', 'ex', 'cert_ok', 'D', 'E_over_Dv', 'E_max', 'shift_ok', 'shift_min_slack', 'T1_ok', 'T1_ratio', 'rel_ok', 'rel_ratio', 'Ee_ok', 'Ee_ratio',
          'mono', 'gslack_min', 'd_ok', 'w_ok', 'omega_act', 'omega_bound', 'formed_ok', 'final_err', 'c_err', 'S_max_over_lambda')


def diagnostics(prep, eps, par):
    """Every ingredient inequality of the proof on the integrated trajectory (returns a row of D_COLS[3:])."""
    o = cert(prep, eps, par)
    if not o['ok']:
        return None
    cs, As, c0 = prep['cs'], prep['As'], prep['c0']
    J, K, L = len(c0), prep['K'], prep['L']
    ev = L['events']
    amin = prep['amin']
    sol, U = simulate(prep, eps, nlocal=1500)
    lam = U.lam
    t = sol.t
    x, c = sol.y[:J], sol.y[J:]
    w = _wfun(x, eps, U)
    eta = np.where(x >= 0, lam - U.s_of_x(np.abs(x)), -(lam - U.s_of_x(np.abs(x))))
    ent = [e['kernel'] for e in ev]
    non = [j for j in range(J) if j not in ent]
    tk = o['tk']
    tb = np.r_[0.0, tk, t[-1]]

    def ubar(j):
        u = np.zeros_like(t)
        for i in range(K + 1):
            m = (t >= tb[i]) & (t <= tb[i + 1])
            u[m] = L['U_br'][i][j] * lam + L['P'][i][j] * (t[m] - tb[i])
        return u
    ratio = max(np.abs(eta[j] - ubar(j)).max() / o['Dv'][j] for j in non)
    Emax = max(np.abs(eta[j] - ubar(j)).max() for j in non)
    gg = (Gam(c[:, None, :] - cs[None, :, None]) * As[None, :, None]).sum(1) - np.einsum('jit,it->jt', Gam(c[:, None, :] - c[None, :, :]), w)
    Bx = o['Bx']
    r, vlo, vbar = par['r'], o['vlo'], o['vbar']
    sh_ok, sh_slack, T1_ok, T1_ratio, rel_ok, rel_ratio, Ee_ok, Ee_ratio, mono, gsl, d_ok, w_ok = True, np.inf, True, 0.0, True, 0.0, True, 0.0, True, np.inf, True, True
    for k, e in enumerate(ev):
        j, f = e['kernel'], e['feature']
        wj = np.abs(w[j])
        ta = t[np.argmax(wj >= o['w1'])]
        t1 = t[np.argmax(wj >= o['vsharp'])]
        d = np.abs(c[j] - cs[f])
        inD = (d <= Bx['dhat'] + 1e-12) & (wj >= Bx['vD_lo'] * amin) & (wj <= Bx['vD_hi'] * amin)
        bad = np.flatnonzero(~inD & (t >= t1))
        trel = t1 if not len(bad) else (t[bad[-1] + 1] if bad[-1] + 1 < len(t) else np.inf)
        shift = ta - tk[k]
        sh_ok &= bool(-o['Sm'][k] - 1e-6 <= shift <= o['Sp'][k] + 1e-6)
        sh_slack = min(sh_slack, min(shift + o['Sm'][k], o['Sp'][k] - shift))
        T1_ok &= bool(t1 - ta <= o['T1'] + 1e-6)
        T1_ratio = max(T1_ratio, (t1 - ta) / o['T1'] if o['T1'] > 0 else 0.0)
        rel_ok &= bool(trel - t1 <= o['T2'] + o['T3'] + 1e-6)
        rel_ratio = max(rel_ratio, (trel - t1) / (o['T2'] + o['T3']))
        Ee = abs(-par['delta1'] - e['p_e'] * (ta - tk[k]))
        Ee_ok &= bool(Ee <= o['De'][k] + 1e-9)
        Ee_ratio = max(Ee_ratio, Ee / o['De'][k])
        mono &= bool(np.all(np.diff(wj[t <= trel]) >= -1e-13 * vbar))
        below = wj < vlo
        if below.any():
            gsl = min(gsl, float((gg[j] * np.sign(As[f]) - (vlo - wj))[below].min()))
        d_ok &= bool(d.max() <= r + 1e-12)
        w_ok &= bool(wj.max() <= vbar + 1e-12)
    om = float(np.abs(w[non]).sum(0).max() / amin)
    thr = min(THETA_REL * np.abs(As).max(), 0.1 * np.abs(As).min())
    formed = sorted(np.flatnonzero(np.abs(w[:, -1]) > thr).tolist())
    ferr = max(abs(abs(w[ent[k], -1]) - abs(As[ev[k]['feature']])) for k in range(K))
    cerr = max(abs(c[ent[k], -1] - cs[ev[k]['feature']]) for k in range(K))
    Smax = max(max(o['Sp']), max(o['Sm'])) / lam
    return [float(o['ok']), o['D'], ratio, Emax, float(sh_ok), sh_slack, float(T1_ok), T1_ratio, float(rel_ok), rel_ratio, float(Ee_ok), Ee_ratio, float(mono), gsl, float(d_ok),
            float(w_ok), om, o['omega_peak'], float(formed == sorted(ent)), ferr, cerr, Smax]


def part_D(quick=False, rawC=None):
    warnings.filterwarnings('ignore')
    tg = C_QUICK if quick else C_TARGETS
    rows = []
    meta = rawC['C_meta']
    for ti, (K, sep, R, seed) in enumerate(tg):
        e0 = meta[ti, 11]
        if not np.isfinite(e0):
            continue
        cs, As, c0 = make_target(K, sep, R, seed)
        prep = prepare(cs, As, c0)
        for dexp in (0.0, -1.0):
            exv = e0 + dexp
            o, pi = best_cert(prep, 10.0 ** exv)
            if pi < 0:
                continue
            row = diagnostics(prep, 10.0 ** exv, PARS[pi])
            if row is not None:
                rows.append([K, sep, exv] + row)
    return dict(D_rows=np.array(rows, float))


# ----------------------------------------------------------------------------------------------------------------------------------
# G. eps_0 against the separation (local data fixed)
# ----------------------------------------------------------------------------------------------------------------------------------
G_MS = (6, 7, 8, 9, 10, 12, 14)                        # separations m h, h = 0.1/sigma (5.0 ... 11.7 widths), fixed offsets and amplitudes
G_QUICK = (8,)


def part_G(quick=False):
    warnings.filterwarnings('ignore')
    ms = G_QUICK if quick else G_MS
    ex = EX_QUICK if quick else EX_GRID
    rows = []
    for m in ms:
        cs = np.array([-0.5 * m * H_LAT + 0.06, 0.5 * m * H_LAT - 0.11])
        As = np.array([0.7, -0.7 * 1.5])
        c0 = np.arange(cs.min() - 8, cs.max() + 8, H_LAT) + 0.3 * H_LAT
        prep = prepare(cs, As, c0)
        okc, okm, okx, Dd, Pi = scan_target(prep, ex)
        e0, em = thr_first(okc, ex), thr_first(okm, ex)
        D0 = float(Dd[np.flatnonzero(ex == e0)[0]]) if np.isfinite(e0) else np.nan
        rows.append([m, prep['Delta'], prep['mu_min'], prep['gap_min'], max(prep['d0']), e0, em, D0])
    return dict(G_rows=np.array(rows, float))


# ----------------------------------------------------------------------------------------------------------------------------------
# H. what sets the measured threshold: the actual dual deviation against the margin lambda mu
# ----------------------------------------------------------------------------------------------------------------------------------
def dual_dev(prep, eps, sol, U):
    """max over the non-entrant kernels and over time of |eta_j - ubar_j| (eta units) in an integrated run."""
    cs, c0 = prep['cs'], prep['c0']
    J, K, L = len(c0), prep['K'], prep['L']
    lam = U.lam
    x = sol.y[:J]
    t = sol.t
    eta = np.where(x >= 0, lam - U.s_of_x(np.abs(x)), -(lam - U.s_of_x(np.abs(x))))
    tb = np.r_[0.0, lam * L['tau_br'][1:], t[-1]]
    ent = [e['kernel'] for e in L['events']]
    best = 0.0
    for j in range(J):
        if j in ent:
            continue
        u = np.zeros_like(t)
        for i in range(K + 1):
            m = (t >= tb[i]) & (t <= tb[i + 1])
            u[m] = L['U_br'][i][j] * lam + L['P'][i][j] * (t[m] - tb[i])
        best = max(best, float(np.abs(eta[j] - u).max()))
    return best


def part_H(quick=False, rawC=None):
    warnings.filterwarnings('ignore')
    tg = C_QUICK if quick else C_TARGETS
    rows = []
    meta = rawC['C_meta']
    for ti, (K, sep, R, seed) in enumerate(tg):
        em = meta[ti, 12]
        if not np.isfinite(em):
            continue
        cs, As, c0 = make_target(K, sep, R, seed)
        prep = prepare(cs, As, c0)
        for dexp in (0.0, 0.5, 1.0):
            eps = 10.0 ** (em + dexp)
            sol, U = simulate(prep, eps, nlocal=1500)
            ok, ex = follows(prep, eps, sol, U)
            rows.append([K, sep, R, ti, em + dexp, U.lam * prep['mu_min'], dual_dev(prep, eps, sol, U), float(ok)])
    return dict(H_rows=np.array(rows, float))


# ----------------------------------------------------------------------------------------------------------------------------------
# E. the repository's network (P0, sampled inner product, wide window)
# ----------------------------------------------------------------------------------------------------------------------------------
def part_E(quick=False):
    warnings.filterwarnings('ignore')
    import exp_theory_h3 as H3
    import exp_parametrization as P
    from src.core import factorized_rbf as FR
    from src.core import rich_transit_nd as RN
    half = 1.5
    xs = np.linspace(-half, half, int(round(2 * half * 30)) + 1)
    Kw = RN.BumpsND(xs[:, None], SIG)
    c0w = np.round(np.arange(-half, half + 1e-9, 0.1), 10)
    Zs = float((Kw.phi(np.array([0.0])) ** 2).mean())
    ro = P.make_readout('P0')
    cells = E_QUICK if quick else E_CELLS
    ex = EX_QUICK if quick else EX_GRID
    meta, OKC, OKR, OKX = [], [], [], []
    for (K, sep, R, mode, idx) in cells:
        cs_phys, As = H3.i_target(K, sep, R, mode, E_SEED + 100 * idx + 5)       # the seed of the cell: E_SEED + 100 idx + 5
        cs, c0 = cs_phys / SIG, c0w / SIG
        prep = prepare(cs, As, c0)
        f = Kw.synth(cs_phys[:, None], As)
        tau_end = (1.05 * prep['L']['events'][-1]['tau'] + 1.0) / Zs
        ent = sorted(e['kernel'] for e in prep['L']['events'])
        thr = min(THETA_REL * np.abs(As).max(), 0.1 * np.abs(As).min())
        okc, okr, okx = [], [], []
        for x_ in ex:
            eps = 10.0 ** x_
            o, pi = best_cert(prep, eps)
            okc.append(bool(o['ok']))
            try:
                res = FR.train_flow(ro, Kw, f, c0w[:, None], eps, tau_end, P.LR_C)
                W, C = res['W'], res['C'][:, :, 0]
                formed = sorted(np.flatnonzero(np.abs(W[-1]) > thr).tolist())
                okr.append(bool(formed == ent and res['status'] == 'ok'))
                okx.append(bool(okr[-1] and all(abs(abs(W[-1][e['kernel']]) - abs(As[e['feature']])) < 1e-3 * abs(As[e['feature']]) and abs(C[-1][e['kernel']] - cs_phys[e['feature']]) < 0.05 * SIG
                                                for e in prep['L']['events'])))
            except Exception:
                okr.append(False)
                okx.append(False)
        meta.append([K, sep, R, float(mode == 'worst'), prep['Delta'], prep['mu_min'], prep['gap_min'], max(prep['d0']), thr_first(okc, ex), thr_first(okr, ex)])
        OKC.append(okc)
        OKR.append(okr)
        OKX.append(okx)
    return dict(E_ex=ex, E_meta=np.array(meta, float), E_okc=np.array(OKC, float), E_okr=np.array(OKR, float), E_okx=np.array(OKX, float), E_Zs=np.array([Zs]))


# ----------------------------------------------------------------------------------------------------------------------------------
# F. the first-order formula
# ----------------------------------------------------------------------------------------------------------------------------------
def part_F(quick=False, rawC=None):
    warnings.filterwarnings('ignore')
    tg = C_QUICK if quick else C_TARGETS
    rows = []
    meta = rawC['C_meta']
    cu = np.pi * LAW['a_max'] / (2 * np.sqrt(2.0) * LAW['G'] * np.sqrt(LAW['gamma'] * LAW['a_max'] ** 2 / (4 * LAW['G'])))
    eps_small = 1e-22
    for ti, (K, sep, R, seed) in enumerate(tg):
        cs, As, c0 = make_target(K, sep, R, seed)
        prep = prepare(cs, As, c0)
        o, pi = best_cert(prep, eps_small)
        if pi < 0:
            continue
        par = PARS[pi]
        L = prep['L']
        ev = L['events']
        ent = [e['kernel'] for e in ev]
        non = [j for j in range(len(c0)) if j not in ent]
        # s_Omega: distance at which the non-entrants' dormant mass n B(s), B(s) ~ 4/(G gamma^2 m s^3), would reach Omega a_min
        sOm = float((len(non) * 4.0 / (LAW['G'] * LAW['gamma'] ** 2 * LAW['m']) / (par['Om'] * prep['amin'])) ** (1 / 3))
        mu_j = np.array([min(1 - abs(L['U_br'][i][j]) for i in range(len(L['U_br']))) for j in non])
        q1 = float(min(mu_j / (o['Dv'][non] + sOm)))
        p = [e['p_e'] for e in ev]
        q2 = np.inf
        for k in range(K - 1):
            lo = max(0.0, (par['delta1'] + o['De'][k + 1]) / p[k + 1])
            hi = max(0.0, (o['De'][k] - par['delta1']) / p[k]) + o['Tlay']
            q2 = min(q2, (L['tau_br'][k + 2] - L['tau_br'][k + 1]) / (lo + hi))
        w0 = (par['delta1'] + o['De'][0]) / p[0]
        q3 = L['tau_br'][1] / w0 if w0 > 0 else np.inf
        eps_star = float((cu * min(q1, q2, q3)) ** 2)
        rows.append([K, sep, R, ti, q1, q2, q3, eps_star, meta[ti, 11], float(o['Dv'][non].max())])
    return dict(F_rows=np.array(rows, float), F_cu=np.array([cu]))


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(str(c) for c in r) + ' |')
    return '\n'.join(out)


def fe(ex):
    return 'none' if not np.isfinite(ex) else '1e%+.1f' % ex


BANNER = 'SUPERSEDED by results/synthetic/theory3n.md: the displacement accounting of the certificate was corrected (see review/REVIEW_roundT3AC.md, REVIEW_roundT3N.md); numbers below use the legacy accounting and are kept as history'


def report(raw):
    L = []
    P = L.append
    P(BANNER + '\n')
    P('# Theorem 1: an explicit epsilon_0 in the controlled regime (generated by experiments/exp_theory_t1.py)\n')
    P('The statements and proofs are in THEORY_T1.md.  Setting: continuum Gram, Gaussian features of the kernels\' width, noiseless, equal masses, n = 1; the network is the '
      'paper\'s P0 unit.  `cert` is the sufficient condition of Theorem T1 (all constants explicit); eps_0 is the largest grid value from which `cert` holds at every smaller '
      'one of the grid (half decades); eps_max is the measured analogue (the network forms exactly the limit object\'s units from there down).  Constants of Prop. A: %s.  '
      'The other constants are double precision with 0.1 %% margins on the quadratures.\n' % ('ball arithmetic (python-flint)' if bool(raw.get('HAVE_ARB', [0])[0]) else 'floats with a 1e-9 margin'))
    if 'A_rows' in raw:
        P('## A. The unit (Lemma 2.1-2.3)\n')
        P('lambda sqrt(eps) -> c_u = pi a_max/(2 sqrt 2 G sqrt(v_s)) = %.4f; B(s) s^3 -> 4/(G gamma^2 m) = %.3f (the dormant coefficient at distance s from the wall, independent of eps).  '
          '`bound` is the explicit bound of Lemma 2.2 (first form), `closed` its closed form rho^{3/2}/(1 - s_e^2) (a_max eps^{3/2}/(m sqrt(2 v_s))) cos th/sin^3 th; the mass identity '
          'int W d eta = a_max^2 sinh^2 x/(2 G m) holds to the relative error in the last column.\n' % (np.pi * LAW['a_max'] / (2 * np.sqrt(2.0) * LAW['G'] * np.sqrt(LAW['gamma'] * LAW['a_max'] ** 2 / (4 * LAW['G']))),
                                                                                     4 / (LAW['G'] * LAW['gamma'] ** 2 * LAW['m'])))
        rows = [['%.0e' % r[0], '%.5g' % r[1], '%.4f' % r[2], '%.2f' % r[3], '%.4e' % r[4], '%.4e' % r[5], '%.4e' % r[6], '%.3f' % r[7], '%.1e' % r[8]] for r in raw['A_rows']]
        P(table(['eps', 'lambda', 'lambda sqrt(eps)', 'mu = s/lambda', 'B(s) exact', 'bound', 'closed', 'B s^3', 'mass-identity err'], rows) + '\n')
    if 'B_rows' in raw:
        P('## B. The regime: Delta_T1 (box conditions (L), (C), (N) of Prop. A and the Gershgorin margins of Lemma 4.1)\n')
        P('Smallest separation (grid 0.1 sigma) with some radius r > d_0 for which every inequality holds in ball arithmetic, at dormant budget Omega = 0 and 0.02 a_min.  '
          'For comparison Delta_entry (THEORY_NOTE, Prop. A alone; d_0 = 0.4167): K=2: 4.66 (R=1), 5.50 (2.4), 6.70 (10); K=3: 5.08 (1), 5.84 (2.4).\n')
        rows = [['%d' % r[0], '%g' % r[1], '%.4g' % r[2], '%.1f' % r[3] if np.isfinite(r[3]) else 'none', '%.1f' % r[4] if np.isfinite(r[4]) else 'none'] for r in raw['B_rows']]
        P(table(['K', 'R', 'd_0', 'Delta_T1 (Omega=0)', 'Delta_T1 (Omega=0.02)'], rows) + '\n')
    if 'C_meta' in raw:
        M, ex = raw['C_meta'], raw['C_ex']
        P('## C. The certificate against the integrated network (%d targets, eps = 10^%.1f ... 10^%.1f)\n' % (len(M), ex[0], ex[-1]))
        P('mu = the limit\'s smallest dormant margin 1 - |u_j|, gap = the smallest relative gap of the entry times, d_0 = the largest entrant offset (widths).  eps_0 = `cert`; eps_max = '
          'measured (continuum network); ratio = eps_max/eps_0.  *violations* = grid values with `cert` true and the network not forming exactly the limit\'s units (the theorem says: none).\n')
        rows = []
        viol = 0
        for i, m in enumerate(M):
            v = int(((raw['C_okc'][i] > 0) & (raw['C_okm'][i] == 0)).sum())
            viol += v
            ratio = 10 ** (m[12] - m[11]) if np.isfinite(m[11]) and np.isfinite(m[12]) else np.nan
            Dd = raw['C_D'][i]
            D0 = Dd[np.flatnonzero(ex == m[11])[0]] if np.isfinite(m[11]) else np.nan
            rows.append(['%d' % m[0], '%g' % m[1] if m[0] > 1 else '-', '%.2f' % m[5], '%d' % m[3], '%.4f' % m[6], '%.3f' % m[7], '%.2f' % m[8], '%d' % m[9],
                         fe(m[11]), fe(m[12]), '%.0f' % ratio if np.isfinite(ratio) else '-', '%.0f' % D0 if np.isfinite(D0) else '-', '%d' % v])
        P(table(['K', 'sep', 'R', 'seed', 'mu', 'gap', 'd_0', 'kernels', 'eps_0 (cert)', 'eps_max (measured)', 'ratio', 'D at eps_0', 'violations'], rows) + '\n')
        C_ = raw['C_comp']
        P('Where the bound D comes from, at eps_0 (eta units, the kernel with the largest bound; *own*: the entrants\' own windows (pre-layer mass, transit, shift, settling); *tails*: the '
          'weighted-max-norm relaxation after each window (Lemma 4.1); *old*: the previously formed units while a later entrant forms; *mass*, *drift*: dormant coefficients and their '
          'centre drift (Lemma 5.2)); S/lambda = the proved bound on the entry-time shift in tau units; T_layer = T_1 + T_2 + T_3 (fast time); params = (r, Omega, delta_1) of PARS.\n')
        rows = []
        for i, m in enumerate(M):
            c = C_[i]
            if not np.isfinite(c[0]):
                continue
            pr = PARS[int(c[9])]
            rows.append(['%d' % m[0], '%g' % m[1] if m[0] > 1 else '-', '%.2f' % m[5], '%d' % m[3], '%.1f' % c[0], '%.1f' % (c[1]), '%.1f' % c[2], '%.1f' % c[3], '%.2g' % c[4], '%.2g' % c[5],
                         '%.1f' % c[6], '%.2e' % c[7], '%.0f' % c[8], '(%g, %g, %g)' % (pr['r'], pr['Om'], pr['delta1'])])
        P(table(['K', 'sep', 'R', 'seed', 'own', 'own tails', 'old', 'old tails', 'mass', 'drift', 'D', 'S/lambda', 'T_layer', 'params'], rows) + '\n')
        Fi = raw['C_fine']
        P('Refinement below eps_0: `cert` also holds at eps_0/10^(1/8), eps_0/10^(1/4), eps_0/10^(3/8) on %d of %d (target, point) pairs (targets without an eps_0 are not counted).\n' % (int(np.nansum(Fi)), int(np.isfinite(Fi).sum())))
        okr = np.array([10 ** (m[12] - m[11]) for m in M if np.isfinite(m[11]) and np.isfinite(m[12])])
        if len(okr):
            P('Total violations over %d (target, eps) pairs: **%d**.  Ratio eps_max/eps_0: median %.0f, range %.0f-%.0f (a factor %.1f-%.1f in eps^(1/2)).  Measured eps_max over the targets: '
              'median 10^%.1f, range 10^%.1f to 10^%.1f (Part E of theory_h3.md: median 1e-5 at 2-8 sigma).\n'
              % (raw['C_okc'].size, viol, np.median(okr), okr.min(), okr.max(), np.sqrt(okr.min()), np.sqrt(okr.max()), np.nanmedian(M[:, 12]), np.nanmin(M[:, 12]), np.nanmax(M[:, 12])))
    if 'D_rows' in raw and len(raw['D_rows']):
        R_ = raw['D_rows']
        P('## D. The ingredients against the trajectory (at eps_0 and eps_0/10)\n')
        P('Each column is the worst case over the rows; *ratio* = measured/bound (must be <= 1).\n')
        cols = D_COLS
        idx = {c: i for i, c in enumerate(cols)}
        n = len(R_)
        P('rows: %d.  `cert` true on all: %d.' % (n, int(R_[:, idx['cert_ok']].sum())))
        P('- dual deviation of the non-entrant kernels / its bound D_j: max ratio **%.3f** (max measured E = %.1f, bound D up to %.0f);' % (R_[:, idx['E_over_Dv']].max(), R_[:, idx['E_max']].max(), R_[:, idx['D']].max()))
        P('- entry shift within [-S^-, S^+]: %d of %d (smallest slack %.1f);' % (int(R_[:, idx['shift_ok']].sum()), n, R_[:, idx['shift_min_slack']].min()))
        P('- transit time to v_lo/2 <= T_1: %d of %d (max ratio %.2f);' % (int(R_[:, idx['T1_ok']].sum()), n, R_[:, idx['T1_ratio']].max()))
        P('- relaxation into D-hat after the transit <= T_2 + T_3: %d of %d (max ratio %.2f);' % (int(R_[:, idx['rel_ok']].sum()), n, R_[:, idx['rel_ratio']].max()))
        P('- entrant dual deviation at the layer start <= D^e_k: %d of %d (max ratio %.2f);' % (int(R_[:, idx['Ee_ok']].sum()), n, R_[:, idx['Ee_ratio']].max()))
        P('- entrants nondecreasing until formed (Lemma 3.1): %d of %d; smallest slack of g_e >= v_lo - v_e: %.3g;' % (int(R_[:, idx['mono']].sum()), n, R_[:, idx['gslack_min']].min()))
        P('- centre offsets <= r: %d of %d; coefficients <= v_bar: %d of %d;' % (int(R_[:, idx['d_ok']].sum()), n, int(R_[:, idx['w_ok']].sum()), n))
        P('- dormant mass: measured max %.2e, bound %.2e (units of a_min), ratio max %.2f;' % (R_[:, idx['omega_act']].max(), R_[:, idx['omega_bound']].max(), (R_[:, idx['omega_act']] / np.maximum(R_[:, idx['omega_bound']], 1e-300)).max()))
        P('- final state: formed set = limit\'s on %d of %d, coefficient error <= %.1e, centre error <= %.1e;' % (int(R_[:, idx['formed_ok']].sum()), n, R_[:, idx['final_err']].max(), R_[:, idx['c_err']].max()))
        P('- the proved bound on the entry-time shift |t - t_k|/lambda: at most %.2e (in tau units, = O(eps^(1/2))).\n' % R_[:, idx['S_max_over_lambda']].max())
    if 'H_rows' in raw and len(raw['H_rows']):
        P('## H. What sets the measured threshold: the actual dual deviation against the margin\n')
        P('At the measured eps_max (and 0.5, 1 decade above) the integrated network\'s largest dual deviation of a non-entrant, E_act = max_{j,t}|eta_j - ubar_j|, against the margin '
          'lambda mu_min (eta units).  *follows* = the network forms exactly the limit\'s units.\n')
        R_ = raw['H_rows']
        rows = [['%d' % r[0], '%g' % r[1] if r[0] > 1 else '-', '%g' % r[2], '%d' % r[3], '%+.1f' % r[4], '%.0f' % r[5], '%.1f' % r[6], '%.2f' % (r[5] / r[6]), 'yes' if r[7] > 0 else 'no'] for r in R_]
        P(table(['K', 'sep', 'R', 'target', 'log10 eps', 'lambda mu', 'E_act', 'lambda mu/E_act', 'follows'], rows) + '\n')
        off = np.array([r[4] - R_[R_[:, 3] == r[3], 4].min() for r in R_])
        at, up = R_[np.isclose(off, 0.0)], R_[np.isclose(off, 1.0)]
        P('At eps_max: lambda mu/E_act has median %.2f (range %.2f-%.2f); one decade above eps_max the network no longer follows on %d of %d targets, and there lambda mu/E_act has median %.2f.\n'
          % (np.median(at[:, 5] / at[:, 6]), (at[:, 5] / at[:, 6]).min(), (at[:, 5] / at[:, 6]).max(), int((up[:, 7] == 0).sum()), len(up), np.median(up[:, 5] / up[:, 6])))
    if 'G_rows' in raw:
        P('## G. eps_0 against the separation (two features, amplitudes 0.7 and -1.05, offsets 0.06 and 0.11 widths from the lattice, separation m h with h = 0.1/sigma)\n')
        rows = [['%d' % r[0], '%.2f' % r[1], '%.4f' % r[2], '%.3f' % r[3], '%.2f' % r[4], fe(r[5]), fe(r[6]), ('%.0f' % (10 ** (r[6] - r[5]))) if np.isfinite(r[5]) and np.isfinite(r[6]) else '-',
                 '%.0f' % r[7] if np.isfinite(r[7]) else '-'] for r in raw['G_rows']]
        P(table(['m', 'Delta', 'mu', 'gap', 'd_0', 'eps_0 (cert)', 'eps_max (measured)', 'ratio', 'D at eps_0'], rows) + '\n')
    if 'E_meta' in raw:
        M, ex = raw['E_meta'], raw['E_ex']
        P('## E. The repository\'s network (P0, 91 samples on [-1.5, 1.5], 31 kernels), interior targets of Part I of theory_h3\n')
        P('Same grid; the `cert` (continuum constants) against the repository\'s own formed set.  Z (sample <phi^2>) = %.5f.\n' % raw['E_Zs'][0])
        rows = []
        viol = 0
        for i, m in enumerate(M):
            v = int(((raw['E_okc'][i] > 0) & (raw['E_okr'][i] == 0)).sum())
            viol += v
            ratio = 10 ** (m[9] - m[8]) if np.isfinite(m[8]) and np.isfinite(m[9]) else np.nan
            rows.append(['%d' % m[0], '%g' % m[1], '%g' % m[2], 'worst' if m[3] else 'random', '%.2f' % m[4], '%.4f' % m[5], '%.3f' % m[6], '%.2f' % m[7], fe(m[8]), fe(m[9]),
                         '%.0f' % ratio if np.isfinite(ratio) else '-', '%d' % v])
        P(table(['K', 'sep', 'R', 'offsets', 'Delta', 'mu', 'gap', 'd_0', 'eps_0 (cert)', 'eps_max (P0, repository)', 'ratio', 'violations'], rows) + '\n')
        P('Violations over %d (target, eps) pairs: **%d**.\n' % (raw['E_okc'].size, viol))
    if 'F_rows' in raw and len(raw['F_rows']):
        P('## F. The first-order formula  eps_* = (c_u min(q_1, q_2, q_3))^2\n')
        P('c_u = lambda sqrt(eps) = %.4f.  q_1 = min_j mu_j/(D_j + s_Omega) (the dormant margins against the deviation bound and the budget distance), q_2 = min_k (tau_{k+1} - tau_k)/(window_k + lead_{k+1}) '
          '(the order condition), q_3 = tau_1/(first lead); all at eps = 1e-22, where the O(eps) terms (dormant mass, centre drift) vanish.  The first-order threshold is eps_* = (c_u min(q_1, q_2, q_3))^2 '
          '(lambda mu >= D + s_Omega  <=>  sqrt(eps) <= c_u mu/(D + s_Omega)); eps_0 is the grid value.\n' % raw['F_cu'][0])
        rows = [['%d' % r[0], '%g' % r[1] if r[0] > 1 else '-', '%g' % r[2], '%.3g' % r[4], '%.3g' % r[5], '%.3g' % r[6], '%.1e' % r[7], fe(r[8]), '%.0f' % r[9]] for r in raw['F_rows']]
        P(table(['K', 'sep', 'R', 'q_1', 'q_2', 'q_3', 'eps_*', 'eps_0 (grid)', 'D_inf'], rows) + '\n')
    return '\n'.join(L)


# ----------------------------------------------------------------------------------------------------------------------------------
def load_raw():
    z = np.load(RAW, allow_pickle=False)
    return {k: z[k] for k in z.files}


def save_raw(raw):
    os.makedirs(OUT, exist_ok=True)
    np.savez(RAW, **raw)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--parts', default='A,B,C,D,H,G,E,F')
    a = ap.parse_args(argv)
    if a.run:
        raw = {'HAVE_ARB': np.array([float(HAVE_ARB)])}
        if a.parts != ap.get_default('parts') and os.path.exists(RAW):                 # a partial run updates the stored raw
            raw.update(load_raw())
        parts = a.parts.split(',')
        t0 = time.time()
        if 'A' in parts:
            raw.update(part_A())
        if 'B' in parts:
            raw.update(part_B())
        if 'C' in parts:
            raw.update(part_C(a.quick))
        if 'D' in parts and 'C_meta' in raw:
            raw.update(part_D(a.quick, raw))
        if 'H' in parts and 'C_meta' in raw:
            raw.update(part_H(a.quick, raw))
        if 'G' in parts:
            raw.update(part_G(a.quick))
        if 'E' in parts:
            raw.update(part_E(a.quick))
        if 'F' in parts and 'C_meta' in raw:
            raw.update(part_F(a.quick, raw))
        save_raw(raw)
        print('raw written (%.0f s)' % (time.time() - t0))
    else:
        raw = load_raw()
    with open(MD, 'w') as fh:
        fh.write(report(raw))
    print('report written:', MD)


if __name__ == '__main__':
    main()
