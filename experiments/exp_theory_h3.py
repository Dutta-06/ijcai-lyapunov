"""
The coupled entry layer (H3) and the separation threshold of Theorem 2, made explicit: a proposition with an explicit delta_0, a certificate
for the K-feature separated regime, and the numbers
====================================================================================================================================

The supplement's Lemma "a coupled entry" (the least complete step of Theorem 2) was a sketch: "a finite-horizon continuation leaves a
threshold delta_0 that may be exponentially small", "the explicit delta_0 ... not carried out". This script checks, and computes the
constants of, the replacement argument (THEORY_NOTE.md in the author's working directory): the inner problem of the entry layer is a
gradient flow of the fit F with a positive diagonal mobility, F = S_1 + ... + S_k + (cross terms of the units' error functions), so a box in
(coefficients, centre offsets) is positively invariant as soon as three sign conditions hold on its faces; none of them involves the
mobility, and none is a finite-horizon estimate.

  A. The constants of the proposition (explicit, mobility-free): for K features at separation >= Delta, amplitude ratio R and the entrant's
     starting offset d_0 (in widths), the smallest Delta at which the sufficient conditions (L), (C), (N) hold for a box of radius r >= d_0
     (the conditions are evaluated in ball arithmetic when python-flint is installed). The table is the proved threshold of the *entry
     layer* alone, for the units' error functions bounded by the triangle inequality.
  B. The inner problem itself (the free law's critical-curve flow with moving centres, ratio of the repository's FreeLaw), two units, the
     first at its critical point and the second entering at a lattice offset d_0 from its feature: whether it ends at the target, the
     largest excursion of each unit, and whether the box of A holds the whole trajectory, at separations on both sides of A's threshold.
  C. A certificate for the limit object ("sliding inverse scale space", `rich_transit_nd.sliding_iss_nd`): from the target's features and the
     starting kernels it computes the entries by the piecewise-linear duals with the blocks at their critical points near the features,
     and checks (a) the next entrant is unique with a relative gap, (b) every other dual is below its wall, (c) A's conditions hold at each
     entry. Theorem K (the note) says: if all hold, the limit object has exactly these entries. MIXED TEST: (a), (b) and the blocks are evaluated in the
     repository's 61-sample inner product, but (c) is the CONTINUUM condition (the proposition assumes <phi_c phi_c'> = Z Gamma(c - c')); "continuum-certified" means that,
     and "interior-certified" adds that every feature is >= 3.5 widths from the window's edge (Part W: the sampled Gram departs from the continuum within ~2 widths of
     the edge) with (c) required under a budget Omega = 0.01 a_min. Compared with the repository's path on the regenerated targets of `exp_separated_regime`
     (K = 2 at 2..14 sigma, K = 3 at 3..7 sigma). The agreement of the certificate's entries with the path's is NOT an independent test of Proposition A (the path
     replaces its flow by a least-squares polish; both choose "the critical point near the features"): see G.
  D. The threshold of the certificate over the local data (the offsets of the two features from the lattice, the amplitudes): the smallest
     separation from which every larger separation of the family (features moved apart by multiples of h) is certified, in sigma, against
     the decoupled margins mu_0, gap; and the same for the dual conditions alone (a)-(b).
  E. The trained network (the paper's unit P0, `exp_parametrization.cell`) on K = 2 targets at 2..10 sigma, eps = 1e-3 ... 1e-10: the largest
     eps below which the formed set is the path's at every smaller eps, per target, against the certificate and the margin.
  F. The mismatch margin of Part C of `exp_parametrization` (200 targets): the per-target eps at which the formed set stops being the
     path's, against the target's runner-up margin, from the stored runs (eps 1e-6, 1e-8) and new P0 runs at eps 1e-5, 1e-7, 1e-9, 1e-10.
  W. The window: the departure of the sampled Gram from the continuum, against the distance of a feature to the window's edge.
  G. The path WITHOUT the polish (`sliding_iss_nd(polish=False)`: the fast flow alone) on 40 targets, against the certificate: a test of the entry flow's convergence itself.
  H. The trained network (P0) against the box of (L), (C): K = 2, eps = 1e-4, 1e-6, 1e-8, the units' excursions and the dormant mass over the run.
  (The skeleton of `exp_transit_certificate` -- fitted variant -- is compared with the certificate's entries in C, column 'skeleton = path'.)

Usage: python experiments/exp_theory_h3.py --run [--parts A,B,C,D,E,F,G,H,W,I] | --from-raw      -> results/synthetic/theory_h3.md (+ _raw.npz)
Environment: TH3_OUT (output directory). Deterministic: every random draw is seeded; the report is regenerated from the raw alone.
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

OUT = os.environ.get('TH3_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory_h3_raw.npz'), os.path.join(OUT, 'theory_h3.md')
PARAM_RAW = os.path.join(REPO, 'results', 'chaos', 'parametrization_raw.npz')

SIG = 0.12                                            # the repository's kernel width; Z below is its sample inner product <phi phi>
RMS_TOL = 1e-6
A_KS, A_RS, A_D0 = (2, 3), (1.0, 2.4, 10.0, 100.0, 1000.0), (0.0, 0.2, 0.4167, 0.6, 0.8)
D0_BENCH = 0.4167                                      # h / (2 sigma) of the benchmark lattice: 0.05 / 0.12
WIN_EDGE = 3.5                                         # 'interior': every feature at least this many widths from the window edge (|c| <= 1 - 3.5 sigma)
OMEGA_WIN = 0.01                                       # the budget (units of a_min) the interior targets' window deviations (<= 3e-5) are required to fit in
A_OM_R, A_OM_D = (1.0, 2.4, 10.0), (5.0, 6.0, 7.0, 8.0, 10.0)
B_DELTAS, B_D0S, B_AMPS, B_V0S = (2.0, 3.0, 4.0, 5.0, 6.0, 8.0), (0.0, 0.4167, 0.8), ((0.7, 0.9), (0.4, 1.0), (0.9, -0.6)), (1e-2, 1e-4)
C_SEEDS = dict(A=1000, B=2000, C=3000)
C_SEPS_A, C_SEPS_B, C_SEPS_C = (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0), (3.0, 4.0, 5.0), (10.0, 12.0, 14.0)
C_SEPS_B3 = (6.0, 7.0)                                 # added here: K = 3 at the separations of the entry layer's own threshold
C_N, C_WIN = 60, 0.95
D_N, D_SEED, D_AMPS = 150, 7001, (0.4, 1.0)
E_SEPS, E_N, E_SEED, E_EPS = (2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0), 16, 8001, tuple(10.0 ** -k for k in range(3, 11))
F_EXTRA_EPS = (1e-5, 1e-7, 1e-9, 1e-10)
F_MAXMARGIN = 0.07                                     # F: targets whose runner-up margin is below this get the extra runs (the rest never mismatch)
_ST = {}


# ------------------------------------------------------------------------------------------------------------------------------------
# A. the entry layer's box: explicit, mobility-free sufficient conditions
# ------------------------------------------------------------------------------------------------------------------------------------
try:
    from flint import arb
    HAVE_ARB = True
except Exception:                                      # python-flint missing: floats with a safety margin (the report says which)
    HAVE_ARB = False


def gam(x):
    return np.exp(-np.asarray(x, float) ** 2 / 4.0)


def psi(x):
    return 0.5 * np.asarray(x, float) * np.exp(-np.asarray(x, float) ** 2 / 4.0)


def g2(x):
    return (np.asarray(x, float) ** 2 / 4.0 - 0.5) * np.exp(-np.asarray(x, float) ** 2 / 4.0)


def _lam_floor(s_max, t_min):
    """Smallest eigenvalue of [[1, s], [s, t]] is increasing in t and decreasing in |s|: its value at (s_max, t_min) bounds it below."""
    return 0.5 * ((1.0 + t_min) - np.sqrt((1.0 - t_min) ** 2 + 4.0 * s_max ** 2))


def box(K, R, Delta, r, use_arb=None, Omega=0.0, need_N=True):
    """The sufficient conditions of the entry-layer proposition for K features of amplitudes in [1, R] (widths = 1), separation >= Delta, a box of
    radius r (|d_i| <= r for every unit): ell = Delta - 2r >= sqrt 6 and kappa = (K - 1) Gamma(ell) < 1; the coefficients are capped at
    vbar = R (1 + kappa)/(1 - kappa); the couplings are bounded by E_v = kappa (R + vbar) (on d/dv) and E_c = (K - 1)(R + vbar) psi(ell) (on d/dc / v), with
    Gamma(x) = e^{-x^2/4}, psi(x) = (x/2) Gamma(x);
      (L)  Gamma(r) > E_v                          the entrant's and the block's coefficients stay above v_lo = Gamma(r) - E_v > 0,
      (C)  psi(r) > E_c                            the centres' offsets point inward on the faces |d| = r,
      (N)  lambda_own - ||H_cross|| > 0            every critical point of the box is a nondegenerate minimum (so there is exactly one).
    need_N=False checks (L), (C) only (enough for the invariance of the box, not for the uniqueness of the critical point). Omega >= 0: the total |coefficient| (in units of a_min) of further units that are not in the box (the network's dormant units): they add Omega to E_v, 0.43 Omega
    (0.43 > sup|Gamma'|) to E_c and the corresponding terms to the cross Hessian. Returns dict(ok, ...) with the margins; with arb the inequalities are evaluated in ball
    arithmetic and `ok` means certainly true."""
    use_arb = HAVE_ARB if use_arb is None else use_arb
    out = dict(ok=False, ell=float(Delta - 2 * r))
    if r >= np.sqrt(2.0) or Delta - 2 * r < np.sqrt(6.0):
        return out
    if use_arb:
        e = lambda x: (-(arb(x) * arb(x)) / 4).exp()
        ps = lambda x: arb(x) / 2 * e(x)
        G2 = lambda x: (arb(x) * arb(x) / 4 - arb(1) / 2) * e(x)
        ell, rr, Rr, Kk = arb(Delta - 2 * r), arb(r), arb(R), arb(K - 1)
        kap = Kk * e(ell)
        if not (kap < 1):
            return out
        Om = arb(Omega)
        vbar = (Rr * (1 + kap) + Om) / (1 - kap)
        Ev = kap * (Rr + vbar) + Om
        Ec = Kk * (Rr + vbar) * ps(ell) + arb(0.43) * Om
        vlo_r = e(rr) - Ev
        Cm = ps(rr) - Ec
        out.update(Ev=float(Ev.mid()), Ec=float(Ec.mid()), vbar=float(vbar.mid()), vlo=float(vlo_r.mid()), Cmargin=float(Cm.mid()))
        if not (vlo_r > 0 and Cm > 0):
            return out
        if not need_N:
            out['ok'] = True
            return out
        dstar = 2 * Ec / e(rr)                          # psi(d) >= (d/2) Gamma(r) on [0, r]: every critical point has |d| <= dstar
        if not (dstar < rr):
            return out
        vlo_d = e(dstar) - Ev
        if not (vlo_d > 0):
            return out
        s_max = Rr * ps(dstar)
        t_min = vlo_d * (arb(1) / 2 - dstar * dstar / 4) * e(dstar)
        lam = (1 + t_min - ((1 - t_min) ** 2 + 4 * s_max * s_max).sqrt()) / 2
        rowv = Kk * (e(ell) + (Rr + vbar) * ps(ell) + vbar * ps(ell)) + arb(0.43) * Om
        rowd = Kk * ((Rr + vbar) * ps(ell) + vbar * (Rr + vbar) * G2(ell) + vbar * ps(ell) + vbar * vbar * G2(ell)) + arb(0.43) * Om + vbar * Om / 2
        Hc = rowv if float(rowv.mid()) > float(rowd.mid()) else rowd
        Nm = lam - Hc
        out.update(dstar=float(dstar.mid()), lam=float(lam.mid()), Hc=float(Hc.mid()), Nmargin=float(Nm.mid()), ok=bool(Nm > 0))
        return out
    ell = Delta - 2 * r
    kap = (K - 1) * float(gam(ell))
    if kap >= 1:
        return out
    vbar = (R * (1 + kap) + Omega) / (1 - kap)
    Ev = kap * (R + vbar) + Omega
    Ec = (K - 1) * (R + vbar) * float(psi(ell)) + 0.43 * Omega
    vlo_r, Cm = float(gam(r)) - Ev, float(psi(r)) - Ec
    out.update(Ev=Ev, Ec=Ec, vbar=vbar, vlo=vlo_r, Cmargin=Cm)
    if vlo_r <= 1e-12 or Cm <= 1e-12:
        return out
    if not need_N:
        out['ok'] = True
        return out
    dstar = 2 * Ec / float(gam(r))
    if dstar >= r:
        return out
    vlo_d = float(gam(dstar)) - Ev
    if vlo_d <= 1e-12:
        return out
    lam = _lam_floor(R * float(psi(dstar)), vlo_d * (0.5 - dstar ** 2 / 4) * float(gam(dstar)))
    rowv = (K - 1) * (float(gam(ell)) + (R + vbar) * float(psi(ell)) + vbar * float(psi(ell))) + 0.43 * Omega
    rowd = (K - 1) * ((R + vbar) * float(psi(ell)) + vbar * (R + vbar) * float(g2(ell)) + vbar * float(psi(ell)) + vbar * vbar * float(g2(ell))) + 0.43 * Omega + vbar * Omega / 2
    Hc = max(rowv, rowd)
    out.update(dstar=dstar, lam=lam, Hc=Hc, Nmargin=lam - Hc, ok=bool(lam - Hc > 1e-12))
    return out


def box_exists(K, R, Delta, d0, n_r=36, use_arb=None, Omega=0.0, need_N=True):
    """The smallest r >= d0 (on a grid up to 1.35) for which `box` holds (a float prefilter, then the ball arithmetic if available), else None."""
    use_arb = HAVE_ARB if use_arb is None else use_arb
    for r in np.linspace(max(d0, 1e-3) * 1.000001 + 1e-9, 1.35, n_r):
        if box(K, R, Delta, float(r), False, Omega, need_N)['ok'] and (not use_arb or box(K, R, Delta, float(r), True, Omega, need_N)['ok']):
            return float(r)
    return None


def entry_threshold(K, R, d0, step=0.02, lo=3.0, hi=24.0):
    """The smallest Delta (on a grid of `step`, from `lo`) at which a box of radius >= d0 exists; nan if none up to `hi`. (The conditions are
    monotone in Delta beyond it: ell grows, every coupling falls, the margins grow.)"""
    for D in np.arange(lo, hi, step):
        if box_exists(K, R, float(D), d0) is not None:
            return float(D)
    return float('nan')


def omega_max(K, R, Delta, d0, need_N, n_r=40):
    """The largest Omega (bisection, floats; confirmed in ball arithmetic) at which (L), (C) [and (N)] hold for some radius r >= d0 (a grid up to 1.35); nan if none at Omega = 0."""
    best_r, best = None, -1.0
    for r in np.linspace(d0 * 1.000001 + 1e-9, 1.35, n_r):
        if not box(K, R, Delta, float(r), False, 0.0, need_N)['ok']:
            continue
        lo, hi = 0.0, 3.0
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if box(K, R, Delta, float(r), False, mid, need_N)['ok'] else (lo, mid)
        if lo > best:
            best_r, best = float(r), lo
    if best_r is None:
        return float('nan')
    lo, hi = best * 0.999, best * 1.001 + 1e-9
    for _ in range(25):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if box(K, R, Delta, best_r, True, mid, need_N)['ok'] else (lo, mid)
    return lo


def part_A():
    res = np.full((len(A_KS), len(A_RS), len(A_D0)), np.nan)
    for a, K in enumerate(A_KS):
        for b, R in enumerate(A_RS):
            for c, d0 in enumerate(A_D0):
                res[a, b, c] = entry_threshold(K, R, d0)
    # the dormant mass the box absorbs: the largest Omega (units of a_min, total |coefficient| of units outside the box) at which the conditions still hold at (K = 2, R, Delta),
    # for the best radius r >= d_0 (a larger r raises psi(r), a smaller one raises Gamma(r)); with (N) (uniqueness) and without it (invariance of the box alone)
    om = np.zeros((len(A_OM_R), len(A_OM_D)))
    omlc = np.zeros((len(A_OM_R), len(A_OM_D)))
    for tab, need_N in ((om, True), (omlc, False)):
        for b, R in enumerate(A_OM_R):
            for c, D in enumerate(A_OM_D):
                tab[b, c] = omega_max(2, R, D, D0_BENCH, need_N)
    return dict(A_thr=res, A_om=om, A_omlc=omlc)


# ------------------------------------------------------------------------------------------------------------------------------------
# B. the inner problem itself: two units, the first at its critical point, the second entering at an offset
# ------------------------------------------------------------------------------------------------------------------------------------
def _law():
    from src.core.mirror_free import FreeLaw
    return FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)


Z_CONT = SIG * np.sqrt(np.pi) / 2                      # the continuum <phi phi> per unit length of the window of length 2: sigma sqrt(pi) / 2


def _Gam(x):
    return np.exp(-np.asarray(x, float) ** 2 / (4 * SIG ** 2))


def _dGam(x):
    return -np.asarray(x, float) / (2 * SIG ** 2) * _Gam(x)


def _gh(c, w, cs, As):
    """g_i = <r phi_i>, h_i = <r d_c phi_i> in the continuum (Z = Z_CONT), r = f* - sum_j w_j phi_j."""
    g = Z_CONT * (_Gam(c[:, None] - cs[None, :]) @ As) - Z_CONT * (_Gam(c[:, None] - c[None, :]) @ w)
    M = _dGam(c[:, None] - c[None, :])
    np.fill_diagonal(M, 0.0)
    h = Z_CONT * (_dGam(c[:, None] - cs[None, :]) @ As) - Z_CONT * (M @ w)
    return g, h


def _crit(cs, As, feats, c_init, w_init):
    k = len(feats)

    def Fg(x):
        g, h = _gh(x[k:], x[:k], cs, As)
        return np.r_[g, x[:k] * h] / Z_CONT
    s = root(Fg, np.r_[w_init, c_init], method='hybr', tol=1e-14)
    return s.x[:k], s.x[k:], bool(np.abs(Fg(s.x)).max() < 1e-12)


def rip_run(Delta, d0, A, v0, t_chunk=3000.0, max_chunks=60):
    """Two units, features at 0 and Delta sigma with amplitudes A; unit 1 at its critical point (unit 2 absent), unit 2 starting at the offset d0
    sigma from feature 2 with the coefficient v0 and flowing by the critical-curve law (`rich_transit.critical_rhs` in the continuum). Returns the final
    (w, c) in sigmas, the largest offsets |d_i| and the coefficients' range over the run, and whether the run reached a stationary point."""
    from src.core.rich_transit import critical_w, curve_mobility, coef_crit, _alpha_of
    L = _law()
    cs, As = np.array([0.0, Delta * SIG]), np.asarray(A, float)
    w1, c1, ok1 = _crit(cs, As, [0], [0.0], [As[0]])
    sg = np.sign(As[1])
    a2 = sg * float(_alpha_of(L, np.array([v0]), coef_crit)[0])
    a1 = np.sign(w1[0]) * float(_alpha_of(L, np.array([abs(w1[0])]), coef_crit)[0])
    z = np.r_[a1, a2, c1[0], cs[1] + d0 * SIG]

    def f(t, z):
        al, c = z[:2], z[2:]
        w = critical_w(L, al)
        g, h = _gh(c, w, cs, As)
        return np.r_[curve_mobility(L, al) * g, L.G * L.m * w * h]
    t, dmax, wmin, wmax, still = 0.0, np.zeros(2), np.full(2, np.inf), np.zeros(2), False
    for _ in range(max_chunks):
        sol = solve_ivp(f, (t, t + t_chunk), z, method='LSODA', rtol=1e-10, atol=1e-14, t_eval=np.linspace(t, t + t_chunk, 61))
        W = critical_w(L, sol.y[:2])
        D = np.abs(sol.y[2:] - cs[:, None]) / SIG
        dmax, wmin, wmax = np.maximum(dmax, D.max(1)), np.minimum(wmin, np.abs(W).min(1)), np.maximum(wmax, np.abs(W).max(1))
        z, t = sol.y[:, -1], t + t_chunk
        w = critical_w(L, z[:2])
        g, h = _gh(z[2:], w, cs, As)
        if np.abs(g).max() < 1e-9 * Z_CONT * np.abs(As).max() and np.abs(w * h).max() < 1e-9 * Z_CONT * np.abs(As).max():
            still = True
            break
    w = critical_w(L, z[:2])
    return dict(w=w, c=z[2:] / SIG, dmax=dmax, wmin=wmin, wmax=wmax, still=still, d_first=abs(c1[0]) / SIG, w_first=abs(w1[0]))


def part_B():
    rows = []
    for D in B_DELTAS:
        for d0 in B_D0S:
            for A in B_AMPS:
                for v0 in B_V0S:
                    a = np.abs(A)
                    R = float(a.max() / a.min())
                    # the box in the normalisation a_min = 1: radius r >= d0 (in widths) and, for the first unit, its own offset
                    r = box_exists(2, R, D, d0)
                    run = rip_run(D, d0, A, v0)
                    amin = float(a.min())
                    exact = bool(run['still'] and np.allclose(np.sort(run['w']), np.sort(np.asarray(A, float)), atol=1e-6)
                                 and np.allclose(run['c'], [0.0, D], atol=1e-6))
                    # inside the box (when it exists): offsets within r, coefficients within [v_lo, vbar] a_min-normalised (the entrant's lower face is its start)
                    inbox = np.nan
                    if r is not None:
                        b = box(2, R, D, r)
                        # units of a_min; unit 1 (the block) stays above v_lo, the entrant above min(its start v_0, v_lo)
                        inbox = float(run['dmax'].max() <= r + 1e-9 and (run['wmax'] / amin).max() <= b['vbar'] + 1e-9 and run['wmin'][0] / amin >= b['vlo'] - 1e-9
                                      and run['wmin'][1] / amin >= min(v0 / amin, b['vlo']) * (1 - 1e-6) - 1e-12)
                    rows.append([D, d0, A[0], A[1], v0, R, float(r is not None), (r if r is not None else np.nan), float(exact), run['dmax'][0], run['dmax'][1],
                                 run['d_first'], inbox, float(run['still'])])
    return dict(B_rows=np.array(rows))


# ------------------------------------------------------------------------------------------------------------------------------------
# C, D. the certificate for the limit object
# ------------------------------------------------------------------------------------------------------------------------------------
def setup():
    if not _ST:
        import exp_separated_regime as S
        L, K, c0, Phi0 = S.setup()
        _ST.update(S=S, L=L, K=K, c0=np.asarray(c0, float), Phi0=Phi0, h=float(c0[1] - c0[0]))
    return _ST['S'], _ST['L'], _ST['K'], _ST['c0'], _ST['Phi0']


class Inner:
    """An inner product with the Gaussian Gram structure: `kind` 'samples' (the repository's 61 samples of [-1, 1], with the window's truncation) or
    'continuum' (<phi(c) phi(c')> = Z Gamma(c - c'), Z = Z_CONT)."""

    def __init__(self, kind, xs=None):
        self.kind = kind
        self.half = 1.0
        if kind == 'samples':
            if xs is None:
                _, _, K, c0, _ = setup()
                xs = K.xs
            self.xs = np.ravel(xs)
            self.N = self.xs.size
            self.half = float(np.abs(self.xs).max())              # the window's half-width (interior: every feature >= WIN_EDGE widths inside it)

    def rates(self, cj, cs, As, w, c):
        """<r phi(c_j)> for the residual r = sum_k A_k phi(c*_k) - sum_i w_i phi(c_i) (arrays of kernels c_j)."""
        if self.kind == 'continuum':
            p = Z_CONT * (_Gam(cj[:, None] - cs[None, :]) @ As)
            return p - Z_CONT * (_Gam(cj[:, None] - c[None, :]) @ w) if len(w) else p
        phi = lambda cc: np.exp(-(self.xs[:, None] - np.asarray(cc, float)[None, :]) ** 2 / (2 * SIG ** 2))
        r = phi(cs) @ As - (phi(c) @ w if len(w) else 0.0)
        return phi(cj).T @ r / self.N

    def block(self, cs, As, feats):
        """The critical point of the fit over len(feats) units started at the features (nan, nan, False if Newton fails)."""
        k = len(feats)
        if self.kind == 'continuum':
            w, c, ok = _crit(cs, As, feats, cs[feats], As[feats])
            return w, c, ok
        phi = lambda cc: np.exp(-(self.xs[:, None] - np.asarray(cc, float)[None, :]) ** 2 / (2 * SIG ** 2))

        def Fg(x):
            w, c = x[:k], x[k:]
            P = phi(c)
            r = phi(cs) @ As - P @ w
            return np.r_[P.T @ r, w * ((P * (self.xs[:, None] - c[None, :]) / SIG ** 2).T @ r)] / self.N
        s = root(Fg, np.r_[As[feats], cs[feats]], method='hybr', tol=1e-14)
        return s.x[:k], s.x[k:], bool(np.abs(Fg(s.x)).max() < 1e-12)


def certify(cs, As, c0, inner, tol=1e-12):
    """The certificate of Theorem K for the limit object on the starting kernels c0 (equal walls). MIXED TEST: (a), (b) (duals, signs, margins) and the blocks'
    critical points are evaluated in the inner product `inner` (the repository's 61 samples, or the continuum), but (c), the box of Proposition A, is the
    CONTINUUM condition (the proposition assumes <phi_c phi_c'> = Z Gamma(c - c')). So `certified` means "continuum (c) holds and (a),(b) hold in `inner`"; it is a
    theorem about the sampled path only where the window's departure from the continuum Gram is negligible: `certified_int` = (a),(b) in the samples, (c) with the
    budget Omega = OMEGA_WIN (units of a_min) to absorb that departure, and every feature at least WIN_EDGE widths from the window's edge (departures <= 1e-4 relative).
    Returns dict with
    entries (kernel indices in order), times, per entry the relative gap to the runner-up (gap), the offset d0 of the entrant from its feature (widths),
    mu = the smallest 1 - |u_j| over the inactive kernels at an entry (the runner-up margin of the paper), dual_ok = (a) and (b) with the signs, box_ok = (c),
    certified = dual_ok and box_ok."""
    cs, As, c0 = np.asarray(cs, float), np.asarray(As, float), np.asarray(c0, float)
    K = len(cs)
    sep = np.abs(np.subtract.outer(cs, cs))[np.triu_indices(K, 1)].min() / SIG
    a = np.abs(As)
    R = float(a.max() / a.min())
    u, act, units = np.zeros(len(c0)), [], []
    w, c = np.zeros(0), np.zeros(0)
    tau = 0.0
    out = dict(entries=[], times=[], gap=[], d0=[], mus=[], sep=float(sep), R=R, dual_ok=True, why='')
    for step in range(K):
        p = inner.rates(c0, cs, As, w, c)
        ina = np.array([j not in act for j in range(len(c0))])
        dt = np.full(len(c0), np.inf)
        ok = ina & (np.abs(p) > 0)
        dt[ok] = (1 - np.sign(p[ok]) * u[ok]) / np.abs(p[ok])
        o = np.argsort(dt)
        j = int(o[0])
        gap = float((dt[o[1]] - dt[o[0]]) / dt[o[0]])
        u = np.where(ina, u + dt[j] * p, u)
        tau += dt[j]
        u[j] = np.sign(p[j])
        oth = [q for q in range(len(c0)) if ina[q] and q != j]
        mu = float(1 - np.abs(u[oth]).max()) if oth else 1.0
        fitted = [f for _, f in units]
        free = [i for i in range(K) if i not in fitted]
        f_as = free[int(np.argmin(np.abs(cs[free] - c0[j])))]
        d0 = float(abs(c0[j] - cs[f_as]) / SIG)
        out['entries'].append(j)
        out['times'].append(float(tau))
        out['gap'].append(gap)
        out['d0'].append(d0)
        out['mus'].append(mu)
        if np.sign(p[j]) != np.sign(As[f_as]):
            out.update(dual_ok=False, why='sign')
            break
        units.append((j, f_as))
        act.append(j)
        feats = [f for _, f in units]
        w, c, okc = inner.block(cs, As, feats)
        if not okc:
            out.update(dual_ok=False, why='newton')
            break
    out['mu'] = float(min(out['mus'])) if out['mus'] else float('nan')
    out['dual_ok'] = bool(out['dual_ok'] and len(out['entries']) == K and out['mu'] > tol and all(g > tol for g in out['gap']))
    rs = [box_exists(K, R, sep, d0) for d0 in out['d0']]
    out['box_ok'] = bool(all(r is not None for r in rs) and len(rs) == K)
    rr = [box_exists(K, R, sep, d0, Omega=OMEGA_WIN) for d0 in out['d0']]
    out['box_ok_rob'] = bool(all(r is not None for r in rr) and len(rr) == K)
    out['interior'] = bool(inner.kind == 'continuum' or np.abs(cs).max() <= inner.half - WIN_EDGE * SIG + 1e-9)
    out['certified'] = bool(out['dual_ok'] and out['box_ok'])                   # continuum conditions (c) + (a),(b) in the inner product `inner`
    out['certified_int'] = bool(out['dual_ok'] and out['box_ok_rob'] and out['interior'])
    return out


def path_of(cs, As):
    """The repository's path on a target: entry list [(tau, kernel)], number of exits, final rms, final active set."""
    import exp_transit_certificate as XC
    import exp_mirror_transit as X
    _, L, K, c0, _ = setup()
    f = X.bumps_target(K, cs, As)
    R = XC.path_1d(L, K, f, c0)
    ev = [(float(t), int(j)) for t, j, k, c in R['events'] if k == 'enter']
    return ev, int(sum(1 for t, j, k, c in R['events'] if k == 'exit')), float(R['hist'][-1][3]), sorted(int(j) for j in np.flatnonzero(R['act']))


def skeleton_entrants(cs, As):
    """The fitted skeleton of `exp_transit_certificate` (its entry list)."""
    import exp_transit_certificate as XC
    _, L, K, c0, _ = setup()
    sk = XC.skeleton(K, c0, np.asarray(cs, float), np.asarray(As, float), fit=True)
    return [int(j) for j in sk['entrants']], bool(sk['ok'])


def c_targets():
    S, L, K, c0, _ = setup()
    out = []
    for sep in C_SEPS_A:
        out += [('A', sep, i, *S.target_a(C_SEEDS['A'] + i, sep, K.sig)) for i in range(C_N)]
    for sep in C_SEPS_C:
        out += [('C', sep, i, *S.target_c(C_SEEDS['C'] + i, sep, K.sig)) for i in range(C_N)]
    for sep in C_SEPS_B:
        out += [('B', sep, i, *S.target_b(C_SEEDS['B'] + i, sep, K.sig)) for i in range(C_N)]
    for sep in C_SEPS_B3:
        out += [('B3', sep, i, *_target_k3(C_SEEDS['B'] + 500 + i, sep, K.sig)) for i in range(C_N)]
    return out


def _target_k3(seed, dmin, sig):
    rng = np.random.default_rng(seed)
    while True:
        cs = rng.uniform(-C_WIN, C_WIN, 3)
        if min(abs(cs[i] - cs[j]) for i in range(3) for j in range(i)) >= dmin * sig:
            return cs, rng.choice([-1.0, 1.0], 3) * rng.uniform(0.4, 1.0, 3)


C_COLS = ('part', 'sep_req', 'sep', 'R', 'dual_ok', 'box_ok', 'certified', 'match_entries', 'time_err', 'path_exact', 'exits', 'skel_match', 'mu', 'gap_min', 'd0_max', 'dec_mu0', 'dec_gap',
          'interior', 'certified_int')


def c_job(item):
    """One target: the certificate (samples), the path, the skeleton; returns a row of C_COLS."""
    part, sep, i, cs, As = item
    S, L, K, c0, Phi0 = setup()
    warnings.filterwarnings('ignore')
    order = np.argsort(cs)
    cs, As = np.asarray(cs)[order], np.asarray(As)[order]
    ce = certify(cs, As, c0, Inner('samples'))
    ev, nex, rms, act = path_of(cs, As)
    pe = [j for t, j in ev]
    match = pe == ce['entries']
    terr = max(abs(t - tt) / tt for (t, j), tt in zip(ev, ce['times'])) if match else np.nan
    sk, skok = skeleton_entrants(cs, As)
    import exp_mirror_transit as X
    dec = S.decoupled(Phi0, lambda c: X.bumps_target(K, [c], [1.0]), cs, As)
    _, _, dmu0, dgap = S.decoupled_summary(dec)
    return [{'A': 0, 'B': 1, 'C': 2, 'B3': 3}[part], sep, ce['sep'], ce['R'], float(ce['dual_ok']), float(ce['box_ok']), float(ce['certified']), float(match), terr,
            float(rms < RMS_TOL and len(act) == len(cs) and nex == 0), float(nex), float(sk == pe), ce['mu'], float(min(ce['gap'])) if ce['gap'] else np.nan,
            float(max(ce['d0'])) if ce['d0'] else np.nan, dmu0, dgap, float(ce['interior']), float(ce['certified_int'])]


def part_C(workers=1):
    items = c_targets()
    rows = [c_job(it) for it in items]
    return dict(C_rows=np.array(rows, float))


def decoupled_entries(cs, As, c0):
    """The decoupled prediction (K one-feature paths in parallel, continuum): per feature its best starting kernel and hitting time T = 1/(Z a Gamma); the kernels in the order of T;
    also mu_0 = the smallest margin 1 - second/best over the features and gap = the smallest relative gap between two features' T."""
    best, T, mu = [], [], []
    for c, A in zip(cs, As):
        a = np.abs(A) * _Gam(np.asarray(c0) - c)
        o = np.argsort(-a)
        best.append(int(o[0]))
        T.append(1.0 / (Z_CONT * a[o[0]]))
        mu.append(1.0 - a[o[1]] / a[o[0]])
    order = np.argsort(T)
    Ts = np.array(T)[order]
    return [best[i] for i in order], float(min(mu)), (float(np.min(np.diff(Ts)) / Ts[-1]) if len(Ts) > 1 else 1.0)


def part_D():
    """The threshold of the certificate over the local data: the continuum inner product, the lattice of the repository (h = 0.1, sigma = 0.12), feature 1 at
    the lattice point -0.8 plus an offset, feature 2 at m h later plus an offset; K = 2, amplitudes +-U(0.4, 1); the smallest separation from which every larger
    separation of the family (m up to the window) is certified (a)-(c), satisfies (a)-(b) alone, or has the decoupled prediction's entries."""
    S, L, K, c0, Phi0 = setup()
    inner = Inner('continuum')
    h = _ST['h']
    rng = np.random.default_rng(D_SEED)
    ms = np.arange(3, 17)
    rows, fam = [], []
    for q in range(D_N):
        o1, o2 = rng.uniform(-h / 2, h / 2, 2)
        A = rng.choice([-1.0, 1.0], 2) * rng.uniform(*D_AMPS, 2)
        cert, dual, dec, seps = [], [], [], []
        for m in ms:
            cs = np.array([-0.8 + o1, -0.8 + m * h + o2])
            ce = certify(cs, A, c0, inner)
            cert.append(ce['certified'])
            dual.append(ce['dual_ok'])
            seps.append(ce['sep'])
            order, mu0, gap = decoupled_entries(cs, A, c0)
            dec.append(ce['entries'] == order)
            if m == 8:
                loc_mu, loc_gap = mu0, gap
        a = np.abs(A)
        rows.append([o1 / h, o2 / h, A[0], A[1], float(a.max() / a.min()), loc_mu, loc_gap])
        fam.append(np.r_[np.array(cert, float), np.array(dual, float), np.array(dec, float), np.array(seps)])
    return dict(D_local=np.array(rows), D_fam=np.array(fam), D_ms=ms)


# ------------------------------------------------------------------------------------------------------------------------------------
# W, G, H. the window, the unpolished flow, the trained network inside the box
# ------------------------------------------------------------------------------------------------------------------------------------
W_EDGE = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)


def part_W():
    """How far the repository's sampled Gram (61 samples of [-1, 1]) departs from the continuum Z Gamma(c - c') for a feature e widths from the window's edge: the largest
    relative departure of <phi_c^2>, of sigma <phi_c d_c phi_c> (0 in the continuum), of sigma^2 <(d_c phi_c)^2> (Z/2) and of the overlap with a kernel 6 widths inward."""
    S, L, K, c0, Phi0 = setup()
    xs = np.ravel(K.xs)
    phi = lambda c: np.exp(-(xs - c) ** 2 / (2 * SIG ** 2))
    dphi = lambda c: phi(c) * (xs - c) / SIG ** 2
    Zi = float(np.mean(phi(0.0) ** 2))
    rows = []
    for e in W_EDGE:
        dev = np.zeros(4)
        for c in np.linspace(1.0 - e * SIG - 0.06, 1.0 - e * SIG + 0.06, 25):      # the offsets are those of a lattice kernel near the feature
            dev[0] = max(dev[0], abs(np.mean(phi(c) ** 2) / Zi - 1))
            dev[1] = max(dev[1], abs(np.mean(phi(c) * dphi(c)) * SIG / Zi))
            dev[2] = max(dev[2], abs(np.mean(dphi(c) ** 2) * SIG ** 2 / Zi / 0.5 - 1))
            c2 = c - 6 * SIG
            dev[3] = max(dev[3], abs(np.mean(phi(c) * phi(c2)) / (Zi * float(_Gam(6 * SIG))) - 1))
        rows.append(dev)
    return dict(W_dev=np.array(rows))


G_GROUPS = (('A', (5.0, 6.0, 7.0, 8.0), 4), ('C', (10.0, 12.0, 14.0), 4), ('B3', (6.0, 7.0), 6))


def part_G():
    """The path's flow WITHOUT the least-squares polish (`sliding_iss_nd(polish=False)`: the fast flow alone, to its stationarity test) on 40 targets spread over Parts A (5-8 sigma),
    C (10-14 sigma) and the K = 3 targets at 6, 7 sigma, against the certificate's entries and the target."""
    import exp_transit_certificate as XC
    import exp_mirror_transit as X
    S, L, K, c0, _ = setup()
    items = c_targets()
    rows = []
    for part, seps, n in G_GROUPS:
        for sep in seps:
            for it in [x for x in items if x[0] == part and x[1] == sep][:n]:
                _, _, _, cs, As = it
                o = np.argsort(cs)
                cs, As = np.asarray(cs)[o], np.asarray(As)[o]
                f = X.bumps_target(K, cs, As)
                t0 = time.time()
                R = XC.path_1d(L, K, f, c0, polish=False)
                secs = time.time() - t0
                ev = [(float(t), int(j)) for t, j, k, c in R['events'] if k == 'enter']
                nex = sum(1 for t, j, k, c in R['events'] if k == 'exit')
                ce = certify(cs, As, c0, Inner('samples'))
                act = np.flatnonzero(R['act'])
                kk = len(cs)
                if len(act) == kk:
                    cf, wf = np.asarray(R['c'], float).reshape(len(c0), -1)[act, 0], np.asarray(R['w'], float).reshape(-1)[act]
                    oo, os_ = np.argsort(cf), np.argsort(cs)
                    cerr, aerr = float(np.abs(cf[oo] - cs[os_]).max()), float(np.abs(wf[oo] - As[os_]).max())
                else:
                    cerr = aerr = np.nan
                t_ev = max(abs(t - tt) / tt for (t, j), tt in zip(ev, ce['times'])) if [j for t, j in ev] == ce['entries'] else np.nan
                rows.append([{'A': 0, 'B3': 3, 'C': 2}[part], sep, kk, float([j for t, j in ev] == ce['entries']), t_ev, float(R['hist'][-1][3]), float(len(act) == kk and nex == 0),
                             float(all(R['relaxed'])), cerr, aerr, float(ce['certified_int']), secs])
    return dict(G_rows=np.array(rows))


H_SEPS, H_N, H_EPS, H_SEED = (6.0, 8.0), 4, (1e-4, 1e-6, 1e-8), 9001


def part_H():
    """The trained network (P0) against the box: K = 2 targets, eps = 1e-4, 1e-6, 1e-8; for the two units that end on the features, over the whole run, |c - c*| <= r, |w| <= v_bar a_min
    and, once |w| exceeds v_lo a_min, |w| >= v_lo a_min, with the box of (L), (C) alone (need_N = False) and Omega = the largest total |w| of the other units over the run (units of a_min)."""
    import exp_parametrization as P
    from src.core import factorized_rbf as FR
    warnings.filterwarnings('ignore')
    S, L, K, c0, Phi0 = setup()
    Lp, Kp, c0p = P.setup()
    ro = P.make_readout('P0')
    rows = []
    for a, sep in enumerate(H_SEPS):
        for i in range(H_N):
            cs, As = S.target_a(H_SEED + 100 * a + i, sep, K.sig)
            o = np.argsort(cs)
            cs, As = cs[o], As[o]
            f = Kp.synth(cs[:, None], As)
            ref = P.reference(Lp, Kp, c0p, f)
            tau_end = 1.05 * max(ref['tau']) + 1.0
            amin, R = float(np.abs(As).min()), float(np.abs(As).max() / np.abs(As).min())
            for eps in H_EPS:
                res = FR.train_flow(ro, Kp, f, c0p, eps, tau_end, P.LR_C)
                W, C = res['W'], res['C'][:, :, 0]
                units = [int(j) for j in ref['set']]
                feat = [int(np.argmin(np.abs(cs - c0p[j]))) for j in units]
                if sorted(feat) != [0, 1]:
                    rows.append([sep, eps, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
                    continue
                dorm = [j for j in range(W.shape[1]) if j not in units]
                om = float((np.abs(W[:, dorm]).sum(1) / amin).max())
                d0 = max(abs(c0p[j] - cs[f]) / SIG for j, f in zip(units, feat))
                r = box_exists(2, R, sep, d0, Omega=om, need_N=False)
                if r is None:
                    rows.append([sep, eps, om, np.nan, np.nan, np.nan, np.nan, np.nan])
                    continue
                b = box(2, R, sep, r, True, om, False)
                dmax, wmax, below = 0.0, 0.0, 0
                for j, f_ in zip(units, feat):
                    dmax = max(dmax, float((np.abs(C[:, j] - cs[f_]) / SIG).max()))
                    v = np.abs(W[:, j]) / amin
                    wmax = max(wmax, float(v.max()))
                    up = np.flatnonzero(v >= b['vlo'])
                    below += int(up.size > 0 and (v[up[0]:] < b['vlo'] * (1 - 1e-9)).any())
                rows.append([sep, eps, om, r, dmax, wmax / b['vbar'], float(below), float(dmax <= r and wmax <= b['vbar'] and below == 0)])
    return dict(H_rows=np.array(rows, float))


# ------------------------------------------------------------------------------------------------------------------------------------
# I. a wider window: the theorem's regime on the real path code, with interior targets for K = 2 and K = 3
# ------------------------------------------------------------------------------------------------------------------------------------
I_HALF, I_SEED, I_N = 1.5, 11000, 12
I_K2_SEPS, I_K3_SEPS, I_RS, I_MODES = (4.5, 5.0, 5.5, 6.0, 8.0), (5.0, 5.5, 6.0, 7.0, 8.0), (1.0, 2.4), ('worst', 'random')
I_COLS = ('K', 'sep_nom', 'R', 'mode', 'sep', 'interior', 'dual_ok', 'box_ok', 'certified_int', 'ent_T', 'ent_F', 'exact_T', 'exact_F', 'cerr_F', 'aerr_F', 'diff_TF', 'time_F', 'd0max', 'gap_min', 'mu')


def i_cells():
    return [(2, d, R, m) for d in I_K2_SEPS for R in I_RS for m in I_MODES] + [(3, d, R, m) for d in I_K3_SEPS for R in I_RS for m in I_MODES]


def i_target(K, sep, R, mode, seed):
    """Features in the window [-1.5, 1.5] (sample spacing 1/30, lattice spacing 0.1, as in the repository). 'worst': every feature about 0.375 widths (0.045 +-10%) from a lattice point -- the nearest-kernel
    offset of the benchmark, 0.4167 widths, would be a tie of two kernels -- the gaps being m h with m = round(sep sigma / h) plus the offset differences; 'random': the first feature at a uniform
    position, the others at exactly `sep` widths. Amplitudes: magnitudes (1, R) or (1, sqrt R, R) in a random order, scaled 0.7, random signs. The configuration is shifted at random
    but stays at least WIN_EDGE widths (plus 0.1) from the edge: every target of Part I is interior."""
    rng = np.random.default_rng(seed)
    h = 0.1
    if mode == 'worst':
        offs = rng.choice([-1.0, 1.0], K) * 0.045 * rng.uniform(0.9, 1.1, K)      # distinct magnitudes: no exact tie of two features' hitting times
        m = int(round(sep * SIG / h))
        span = (K - 1) * m * h
    else:
        span = (K - 1) * sep * SIG
    lim = I_HALF - WIN_EDGE * SIG - 0.1
    shift_max = max(0.0, min(0.3, lim - span / 2))
    shift = rng.uniform(-shift_max, shift_max)
    if mode == 'worst':
        cs = np.round((np.arange(K) * m * h - span / 2 + shift) / h) * h + offs
    else:
        cs = np.arange(K) * sep * SIG - span / 2 + shift + rng.uniform(-h / 2, h / 2)
    mags = np.array([1.0, R] if K == 2 else [1.0, np.sqrt(R), R])[rng.permutation(K)] * 0.7
    return cs, mags * rng.choice([-1.0, 1.0], K)


def part_I():
    import exp_transit_certificate as XC
    import exp_mirror_transit as X
    from src.core import rich_transit as RT
    warnings.filterwarnings('ignore')
    S, L, K0, c0_, _ = setup()
    xs = np.linspace(-I_HALF, I_HALF, int(round(2 * I_HALF * 30)) + 1)
    Kw = RT.Bumps(xs, SIG)
    c0w = np.round(np.arange(-I_HALF, I_HALF + 1e-9, 0.1), 10)
    inner = Inner('samples', xs)
    rows = []
    for ci_, (K, sep, R, mode) in enumerate(i_cells()):
        for i in range(I_N):
            cs, As = i_target(K, sep, R, mode, I_SEED + 100 * ci_ + i)
            ce = certify(cs, As, c0w, inner)
            f = X.bumps_target(Kw, cs, As)
            out = {}
            for pol in (True, False):
                t0 = time.time()
                Rr = XC.path_1d(L, Kw, f, c0w, polish=pol)
                secs = time.time() - t0
                ev = [int(j) for t, j, k, c in Rr['events'] if k == 'enter']
                nex = sum(1 for t, j, k, c in Rr['events'] if k == 'exit')
                act = np.flatnonzero(Rr['act'])
                if len(act) == K:
                    cf = np.asarray(Rr['c'], float).reshape(len(c0w), -1)[act, 0]
                    wf = np.asarray(Rr['w'], float).reshape(-1)[act]
                    oo = np.argsort(cf)
                    cerr, aerr = float(np.abs(cf[oo] - cs[np.argsort(cs)]).max()), float(np.abs(wf[oo] - As[np.argsort(cs)]).max())
                else:
                    cf = wf = None
                    cerr = aerr = np.nan
                exact = float(Rr['hist'][-1][3] < RMS_TOL and len(act) == K and nex == 0)
                out[pol] = (ev, exact, cerr, aerr, cf, wf, secs)
            dif = (max(np.abs(np.sort(out[True][4]) - np.sort(out[False][4])).max(), np.abs(out[True][5][np.argsort(out[True][4])] - out[False][5][np.argsort(out[False][4])]).max())
                   if out[True][4] is not None and out[False][4] is not None else np.nan)
            rows.append([K, sep, R, float(mode == 'worst'), ce['sep'], float(ce['interior']), float(ce['dual_ok']), float(ce['box_ok']), float(ce['certified_int']),
                         float(out[True][0] == ce['entries']), float(out[False][0] == ce['entries']), out[True][1], out[False][1], out[False][2], out[False][3], dif, out[False][6],
                         float(max(ce['d0'])) if ce['d0'] else np.nan, float(min(ce['gap'])) if ce['gap'] else np.nan, ce['mu']])
    return dict(I_rows=np.array(rows, float))


# ------------------------------------------------------------------------------------------------------------------------------------
# E, F. the trained network
# ------------------------------------------------------------------------------------------------------------------------------------
def part_E():
    import exp_parametrization as P
    warnings.filterwarnings('ignore')
    S, L, K, c0, Phi0 = setup()
    Lp, Kp, c0p = P.setup()
    ro = P.make_readout('P0')
    rng_seed = E_SEED
    res = np.full((len(E_SEPS), E_N, len(E_EPS), 5), np.nan)       # set_equal, exact, order_equal, tau_rel, n_formed
    meta = np.full((len(E_SEPS), E_N, 8), np.nan)                  # certified, dual_ok, mu (path), margin_never, R, sep, d0max, cert predicted entries == path
    for a, sep in enumerate(E_SEPS):
        for i in range(E_N):
            cs, As = S.target_a(rng_seed + 100 * a + i, sep, K.sig)
            order = np.argsort(cs)
            cs, As = cs[order], As[order]
            f = Kp.synth(cs[:, None], As)
            ref = P.reference(Lp, Kp, c0p, f)
            ce = certify(cs, As, c0, Inner('samples'))
            meta[a, i] = [float(ce['certified']), float(ce['dual_ok']), ref['margin'], ref['margin_never'], ce['R'], ce['sep'], max(ce['d0']) if ce['d0'] else np.nan,
                          float([j for t, j in sorted(zip(ref['tau'], ref['order']))] == ce['entries'])]
            tau_end = 1.05 * max(ref['tau']) + 1.0
            for b, eps in enumerate(E_EPS):
                c = P.cell(ro, Kp, f, c0p, cs, As, eps, ref, tau_end)
                res[a, i, b] = [float(c['set_equal']), float(c['exact']), float(c['order_equal']), c['tau_rel'], float(c['n_formed'])]
    return dict(E_res=res, E_meta=meta)


def part_F():
    """Per-target eps-thresholds on the 200 targets of Part C of `exp_parametrization`: the stored P0/D3 runs at 1e-6, 1e-8 and new P0 runs."""
    import exp_cohesion_leak as EC
    import exp_parametrization as P
    warnings.filterwarnings('ignore')
    raw = EC.load(PARAM_RAW)
    Cc = list(raw['C'])
    Lp, Kp, c0p = P.setup()
    ro = P.make_readout('P0')
    n = len(Cc)
    eps_all = np.array(sorted(set([1e-6, 1e-8] + list(F_EXTRA_EPS)), reverse=True))
    P0 = np.full((n, len(eps_all)), np.nan)               # set_equal for P0
    D3 = np.full((n, 2), np.nan)
    mar = np.zeros((n, 2))
    for t, tg in enumerate(Cc):
        mar[t] = [tg['ref']['margin'], tg['ref']['margin_never']]
        for cell in tg['cells']:
            if cell['name'] == 'P0' and cell['eps'] in (1e-6, 1e-8):
                P0[t, int(np.flatnonzero(eps_all == cell['eps'])[0])] = float(cell['set_equal'])
            if cell['name'] == 'D3' and cell['eps'] in (1e-6, 1e-8):
                D3[t, 0 if cell['eps'] == 1e-6 else 1] = float(cell['set_equal'])
        if not (tg['ref']['margin'] < F_MAXMARGIN or np.isnan(tg['ref']['margin'])):
            P0[t, :] = np.where(np.isnan(P0[t]), 1.0, P0[t])            # margin above the cut: not run at the extra eps (the stored runs never mismatch there)
            continue
        cs, As = np.asarray(tg['cs'], float), np.asarray(tg['As'], float)
        f = Kp.synth(cs[:, None], As)
        ref = tg['ref']
        tau_end = 1.05 * max(ref['tau']) + 1.0
        for b, eps in enumerate(eps_all):
            if np.isnan(P0[t, b]):
                c = P.cell(ro, Kp, f, c0p, cs, As, float(eps), ref, tau_end)
                P0[t, b] = float(c['set_equal'])
    return dict(F_eps=eps_all, F_P0=P0, F_D3=D3, F_margin=mar)


# ------------------------------------------------------------------------------------------------------------------------------------
# the report
# ------------------------------------------------------------------------------------------------------------------------------------
def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    return '\n'.join(out + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows])


def fmt(x, f='%.2f'):
    return '-' if x is None or (isinstance(x, float) and np.isnan(x)) else f % x


def report(raw):
    L = ['# The coupled entry layer (H3) and the separation threshold: explicit constants, a certificate, and the numbers', '',
         'Generated by `experiments/exp_theory_h3.py` from `theory_h3_raw.npz`; every number below is printed by the script (protocol: its docstring). '
         'Lengths are in kernel widths sigma; Gamma(d) = exp(-d^2 / 4) is the overlap of two kernels at the distance d sigma. Interval arithmetic (python-flint): %s.' %
         ('used' if bool(raw['have_arb'][0]) else 'not available, floats with a safety margin'), '']
    if 'A_thr' in raw:
        t = raw['A_thr']
        L += ['## A. The entry layer: the box of the proposition and its threshold', '',
              'Smallest separation Delta (in widths) at which the sufficient conditions (L), (C), (N) hold for a box of radius r >= d_0 (the entrant starts d_0 widths from its '
              'feature; K features at separation >= Delta, amplitudes within a ratio R). The conditions do not involve the units\' mobilities, the masses, the walls or the '
              'length of the slide; they bound each error function by the triangle inequality.', '']
        for a, K in enumerate(A_KS):
            L += ['K = %d' % K, '', table(['R \\ d_0'] + ['%g' % d for d in A_D0], [['%g' % R] + [fmt(t[a, b, c]) for c in range(len(A_D0))] for b, R in enumerate(A_RS)]), '']
        if 'A_om' in raw:
            om = raw['A_om']
            L += ['The dormant mass the box absorbs (K = 2, d_0 = 0.4167; the largest Omega, in units of a_min, of the total |coefficient| of units outside the box at which (L), (C), (N) still hold; '
                  'it is what makes the box invariant for the trained network, whose dormant units carry coefficients):', '',
                  table(['R \\ Delta'] + ['%g' % d for d in A_OM_D], [['%g' % R] + [fmt(om[b, c], '%.4f') for c in range(len(A_OM_D))] for b, R in enumerate(A_OM_R)]), '']
        if 'A_omlc' in raw:
            oml = raw['A_omlc']
            L += ['The same for the invariance of the box alone, (L) and (C) without (N) (what the network needs; uniqueness of the critical point is not needed there):', '',
                  table(['R \\ Delta'] + ['%g' % d for d in A_OM_D], [['%g' % R] + [fmt(oml[b, c], '%.4f') for c in range(len(A_OM_D))] for b, R in enumerate(A_OM_R)]), '']
    if 'B_rows' in raw:
        B = raw['B_rows']
        L += ['## B. The inner problem, two units, against the box', '',
              'The free law (m = 12, G = 1, gamma = 0.1, a_max = 1) on its critical curve with moving centres, continuum inner product; features at 0 and Delta sigma; unit 1 at its critical '
              'point (unit 2 absent); unit 2 enters at the offset d_0 with the coefficient v_0. "exact": the run stops at the target (coefficients and centres to 1e-6). '
              '"box": the proposition\'s box exists at (Delta, d_0, R); "in box": the whole trajectory stays in it (offsets <= r, coefficients within [v_lo, v_bar] a_min).', '']
        rows = []
        for D in B_DELTAS:
            m = B[:, 0] == D
            bx = m & (B[:, 6] > 0)
            rows.append(['%g' % D, '%d / %d' % (int(B[m, 8].sum()), int(m.sum())), '%d / %d' % (int(bx.sum()), int(m.sum())),
                         '%d / %d' % (int(B[bx, 8].sum()), int(bx.sum())), '%d / %d' % (int(np.nansum(B[bx, 12])), int(bx.sum())),
                         fmt(float(np.nanmax(B[m, 9])), '%.3f'), fmt(float(np.nanmax(B[m, 10])), '%.3f'), fmt(float(np.nanmax(B[m, 11])), '%.3f')])
        L += [table(['Delta', 'ends at the target', 'box exists', 'exact where the box exists', 'in box where it exists', 'max |d_1|', 'max |d_2|', 'unit 1 offset at start'], rows), '']
        bad = B[(B[:, 8] == 0)]
        if len(bad):
            L += ['Runs that do not end at the target (%d of %d): Delta in %s, amplitudes %s; in them the first unit\'s critical point (unit 2 absent) lies %.1f-%.1f widths from its feature '
                  '(the one-unit fit has no critical point near the weak feature next to a feature 2.5 times larger, and Newton lands on the strong one), so the premise of the entry layer '
                  '(a block at a critical point near its feature) fails; at Delta >= 4 every run ends at the target, and wherever the box exists the whole trajectory stays in it.' % (
                      len(bad), len(B), [float(x) for x in sorted(set(bad[:, 0].tolist()))], [(float(x), float(y)) for x, y in sorted(set((bad[i, 2], bad[i, 3]) for i in range(len(bad))))],
                      float(bad[:, 11].min()), float(bad[:, 11].max())), '']
    if 'C_rows' in raw:
        C = raw['C_rows']
        ci = {c: i for i, c in enumerate(C_COLS)}
        L += ['## C. The certificate against the repository\'s path', '',
              'MIXED TEST, to be read as such. The certificate computes the entries from the duals with the blocks at their critical points near the features and checks the dual conditions '
              '(a)-(b) -- all in the repository\'s inner product (61 samples of [-1, 1], the 21 starting kernels, equal walls) -- and the box (c) of Proposition A, which is a CONTINUUM condition '
              '(the proposition assumes <phi_c phi_c\'> = Z Gamma(c - c\')). "continuum-certified" = (c) holds in the continuum and (a)-(b) hold in the samples. The sampled Gram departs from the continuum '
              'near the window\'s edge (Part W): "interior" = every feature at least %g widths from the edge (|c| <= %.2f), where the departures are <= 1e-4 (Part W) and are taken as absorbed by requiring (c) with the budget Omega = %g a_min (a perturbation of the Gram entries by delta moves the gradient and the Hessian by O(delta (R + v_bar)), dominated by Omega when delta <= 1e-3; this absorption is argued, not proved here); '
              '"interior-certified" = interior, (a)-(b), and (c) with that budget: only these are covered by the theorem for the sampled path. '
              '"entries = path": the certificate\'s kernels, in order, are the path\'s. THIS AGREEMENT IS NOT AN INDEPENDENT TEST OF PROPOSITION A: `sliding_iss_nd` integrates the fast flow and then '
              'REPLACES it by a least-squares polish, and both select "the critical point near the features"; Part G runs the path without the polish. '
              'Times: the largest relative difference between the certificate\'s entry times and the path\'s over targets with equal entries. '
              '"skeleton": the fitted skeleton of `exp_transit_certificate` lists the same entries.' % (WIN_EDGE, 1.0 - WIN_EDGE * SIG, OMEGA_WIN), '']
        for part, title in ((0, 'K = 2, separation exactly Delta (part A of `exp_separated_regime`; centres in [-0.7, 0.7])'), (2, 'K = 2, wider window (part C of `exp_separated_regime`, added after the full-draft review; centres in [-0.95, 0.95])'),
                            (1, 'K = 3, minimal separation at least Delta (part B of `exp_separated_regime`; the rows are the draw\'s lower bound)'), (3, 'K = 3 at 6 and 7 sigma, centres in [-0.95, 0.95] (added here)')):
            rows = []
            for s_ in sorted(set(C[C[:, 0] == part, 1].tolist())):
                m = (C[:, 0] == part) & (C[:, 1] == s_)
                n = int(m.sum())
                te = C[m & (C[:, ci['match_entries']] > 0), ci['time_err']]
                rows.append(['%g' % s_, n, '%d' % int(C[m, ci['dual_ok']].sum()), '%d' % int(C[m, ci['match_entries']].sum()), '%d' % int(C[m, ci['skel_match']].sum()),
                             '%d' % int(C[m, ci['box_ok']].sum()), '%d' % int(C[m, ci['certified']].sum()), '%d' % int(C[m, ci['interior']].sum()),
                             '%d' % int(C[m, ci['certified_int']].sum()), '%d' % int((C[m, ci['certified']] * (1 - C[m, ci['interior']])).sum()),
                             '%d' % int((C[m, ci['certified']] * C[m, ci['match_entries']]).sum()), '%d' % int(C[m, ci['path_exact']].sum()),
                             ('%.1e' % te.max()) if te.size else '-'])
            L += ['**%s**' % title, '', table(['Delta (requested)', 'targets', '(a)-(b) hold', 'entries = path', 'skeleton = path', '(c) continuum box', 'continuum-certified', 'interior targets',
                                                'interior-certified', 'edge targets continuum-certified', 'continuum-certified and = path', 'path exact, K units', 'max time difference'], rows), '']
        tot, ti = int(C[:, ci['certified']].sum()), int(C[:, ci['certified_int']].sum())
        edge = int((C[:, ci['certified']] * (1 - C[:, ci['interior']])).sum())
        L += ['Over the %d targets: continuum-certified %d (interior %d, edge %d), interior-certified %d; among the continuum-certified the path\'s entries differ from the certificate\'s on %d; '
              'among the interior-certified on %d; the certificate\'s entries are the path\'s on %d of %d.' % (
                  len(C), tot, int((C[:, ci['certified']] * C[:, ci['interior']]).sum()), edge, ti, int((C[:, ci['certified']] * (1 - C[:, ci['match_entries']])).sum()),
                  int((C[:, ci['certified_int']] * (1 - C[:, ci['match_entries']])).sum()), int(C[:, ci['match_entries']].sum()), len(C)), '']
    if 'W_dev' in raw:
        Wd = raw['W_dev']
        L += ['## W. The window: how far the sampled Gram is from the continuum', '',
              'Largest relative departure, over feature positions e widths from the edge of the 61-sample window [-1, 1], of: <phi_c^2> from its interior value; sigma<phi_c, d_c phi_c> (zero in the continuum; '
              'units of Z); sigma^2<(d_c phi_c)^2> from Z/2; the overlap with a kernel 6 widths inward from Z Gamma(6). Proposition A\'s budget for such departures is Omega (Part A: (L),(C) tolerate '
              'Omega up to about 0.6-0.7 a_min, (N) up to 0.04-0.26).', '',
              table(['e (widths from the edge)', '<phi^2>', 'sigma<phi, d phi>', 'sigma^2<(d phi)^2>', 'overlap at 6 widths'], [['%g' % e] + ['%.1e' % x for x in Wd[i]] for i, e in enumerate(W_EDGE)]), '',
              'A feature closer than about 2 widths to the edge has departures of 0.1-0.9, larger than the budget: the continuum conditions are NOT shown to cover those targets (a corrected '
              'derivation would have to carry the window-dependent self-terms, whose centre force sigma<phi, d phi> is of order 1 there, through (C) and (N); not done). '
              'At %g widths the departures are <= 1e-4, far inside the budget Omega = %g a_min (the rows for e <= 3 show up to 1e-3 in the curvature).' % (WIN_EDGE, OMEGA_WIN), '']
    if 'G_rows' in raw:
        Gr = raw['G_rows']
        L += ['## G. The path without the polish', '',
              '`sliding_iss_nd(polish=False)` (the fast flow alone to its stationarity test; no least-squares step replaces it) on %d targets spread over Part A at 5-8 sigma (4 per separation), Part C at 10-14 sigma (4) '
              'and the K = 3 targets at 6, 7 sigma (6): entries equal to the certificate\'s, final rms, exact fit with K units and no exit, every block relaxed (the flow\'s own test), '
              'largest centre / coefficient error against the target.' % len(Gr), '',
              table(['targets', 'K', 'entries = certificate', 'max time difference', 'max final rms', 'exact, K units, no exit', 'all blocks relaxed', 'max centre error', 'max coefficient error', 'interior-certified', 'max seconds'],
                    [['%d' % len(Gr[Gr[:, 2] == k]), '%d' % k, '%d' % int(Gr[Gr[:, 2] == k, 3].sum()), '%.1e' % np.nanmax(Gr[Gr[:, 2] == k, 4]), '%.1e' % Gr[Gr[:, 2] == k, 5].max(),
                      '%d' % int(Gr[Gr[:, 2] == k, 6].sum()), '%d' % int(Gr[Gr[:, 2] == k, 7].sum()), '%.1e' % np.nanmax(Gr[Gr[:, 2] == k, 8]), '%.1e' % np.nanmax(Gr[Gr[:, 2] == k, 9]),
                      '%d' % int(Gr[Gr[:, 2] == k, 10].sum()), '%.1f' % Gr[Gr[:, 2] == k, 11].max()] for k in (2, 3)]), '']
    if 'H_rows' in raw:
        Hr = raw['H_rows']
        L += ['## H. The trained network (P0) inside the box of (L), (C)', '',
              'K = 2 targets (%d per separation, seeded), eps = %s. For the two units that end on the features: over the whole run, |c - c*| <= r, |w| <= v_bar a_min and, once |w| >= v_lo a_min, |w| stays >= v_lo a_min; '
              'the box has (L), (C) only (invariance does not need (N)) and Omega = the largest total |w| of the other 19 units over the run, in units of a_min. "box exists": it exists for that Omega and the starting offsets. '
              '"max |d|/r": the largest centre excursion over r.' % (H_N, ', '.join('%g' % e for e in H_EPS)), '']
        rows = []
        for sep in H_SEPS:
            for e in H_EPS:
                m = (Hr[:, 0] == sep) & (Hr[:, 1] == e)
                bx = m & ~np.isnan(Hr[:, 3])
                rows.append(['%g' % sep, '%g' % e, '%d' % int(m.sum()), '%d' % int(bx.sum()), '%d' % int(np.nansum(Hr[bx, 7])), '%.1e' % np.nanmax(Hr[m, 2]),
                             ('%.2f' % np.nanmax(Hr[bx, 4] / Hr[bx, 3])) if bx.any() else '-', ('%.2f' % np.nanmax(Hr[bx, 5])) if bx.any() else '-'])
        L += [table(['Delta', 'eps', 'runs', 'box exists', 'trajectory inside the box', 'max Omega observed', 'max |d|/r', 'max |w|/(v_bar a_min)'], rows), '']
    if 'I_rows' in raw:
        Ir = raw['I_rows']
        ic = {c: i for i, c in enumerate(I_COLS)}
        beh = (Ir[:, ic['ent_F']] * Ir[:, ic['exact_F']] > 0)
        L += ['## I. A wider window: the theorem\'s regime on the real path code, every target interior', '',
              'Window [-%g, %g] with the repository\'s sample spacing (1/30) and kernel lattice (spacing 0.1, width 0.12), so that every feature is at least %g widths (plus 0.1) from the edge: all %d targets are interior '
              '(`certify`\'s `interior` is %d of %d); (c) is the continuum condition with the budget Omega = %g a_min as before, (a)-(b) are evaluated in the samples. %d targets per cell, deterministic seeds. '
              'K = 2 at nominal separation %s widths, K = 3 at %s; amplitude ratio R = %s (magnitudes (1, R) or (1, sqrt R, R), scaled 0.7, random order and signs); "worst": every feature 0.375 widths '
              '(0.045, +-10%%) from a lattice point, the gaps m h plus the offset differences (the largest nearest-kernel offset, 0.4167 widths, would tie two kernels); "random": the first feature at a uniform position, '
              'the others exactly at the nominal separation. `sep` is the actual minimal separation (it differs from the nominal one by up to ~0.4 widths in "worst"). For each target the real path code `sliding_iss_nd` is run '
              'WITH the least-squares polish and WITHOUT it (the flow alone). "behaves as the theorem says": the unpolished path\'s entries are the certificate\'s, in order, and it ends at an exact fit with K units, '
              'no exit (final rms < %g). "cert. but not behaving" counts interior-certified targets where it does not (the theorem\'s prediction: 0); "behaves, not certified" shows the conservativeness.' % (
                  I_HALF, I_HALF, WIN_EDGE, len(Ir), int(Ir[:, ic['interior']].sum()), len(Ir), OMEGA_WIN, I_N, ', '.join('%g' % x for x in I_K2_SEPS), ', '.join('%g' % x for x in I_K3_SEPS),
                  ', '.join('%g' % x for x in I_RS), RMS_TOL), '']
        rows = []
        for K in (2, 3):
            for d in sorted(set(Ir[Ir[:, 0] == K, 1].tolist())):
                m = (Ir[:, 0] == K) & (Ir[:, 1] == d)
                rows.append(['%d' % K, '%g' % d, int(m.sum()), '%.2f-%.2f' % (Ir[m, ic['sep']].min(), Ir[m, ic['sep']].max()), int(Ir[m, ic['dual_ok']].sum()), int(Ir[m, ic['box_ok']].sum()),
                             int(Ir[m, ic['certified_int']].sum()), int(Ir[m, ic['ent_T']].sum()), int(Ir[m, ic['ent_F']].sum()), int(Ir[m, ic['exact_T']].sum()), int(Ir[m, ic['exact_F']].sum()),
                             int(beh[m].sum()), int((Ir[m, ic['certified_int']] * (~beh[m])).sum()), int(((Ir[m, ic['certified_int']] == 0) & beh[m]).sum()),
                             '%.1e' % np.nanmax(Ir[m, ic['diff_TF']]), '%.1e' % np.nanmax(Ir[m, ic['cerr_F']])])
        L += [table(['K', 'nominal Delta', 'targets', 'actual sep', '(a)-(b) hold', '(c) box', 'interior-certified', 'entries = cert. (polish)', 'entries = cert. (no polish)', 'exact (polish)',
                     'exact (no polish)', 'behaves as the theorem says', 'cert. but not behaving', 'behaves, not certified', 'max |polished - unpolished| (c, w)', 'max centre error (no polish)'], rows), '']
        tot = int(Ir[:, ic['certified_int']].sum())
        L += ['Over the %d targets: interior-certified %d (K = 2: %d of %d; K = 3: %d of %d); the path (no polish) behaves as the theorem says on %d, with the polish on %d; interior-certified targets where it does not: %d; '
              'targets where it behaves although not certified: %d (the conservativeness of Delta_entry). The polished and unpolished paths have the same entries on %d of %d.' % (
                  len(Ir), tot, int(Ir[Ir[:, 0] == 2, ic['certified_int']].sum()), int((Ir[:, 0] == 2).sum()), int(Ir[Ir[:, 0] == 3, ic['certified_int']].sum()), int((Ir[:, 0] == 3).sum()),
                  int(beh.sum()), int((Ir[:, ic['ent_T']] * Ir[:, ic['exact_T']] > 0).sum()), int((Ir[:, ic['certified_int']] * (~beh)).sum()), int(((Ir[:, ic['certified_int']] == 0) & beh).sum()),
                  int((Ir[:, ic['ent_T']] == Ir[:, ic['ent_F']]).sum()), len(Ir)), '']
        rows = []
        for K in (2, 3):
            for R in I_RS:
                for mode in (1.0, 0.0):
                    for d in sorted(set(Ir[Ir[:, 0] == K, 1].tolist())):
                        m = (Ir[:, 0] == K) & (Ir[:, 2] == R) & (Ir[:, 3] == mode) & (Ir[:, 1] == d)
                        rows.append(['%d' % K, '%g' % R, 'worst' if mode else 'random', '%g' % d, int(m.sum()), '%d' % int(Ir[m, ic['certified_int']].sum()), '%d' % int(beh[m].sum()),
                                     '%.2f' % Ir[m, ic['sep']].mean(), '%.1e' % np.nanmin(Ir[m, ic['gap_min']]), '%.1e' % np.nanmin(Ir[m, ic['mu']])])
        L += ['By cell (certified versus behaving; `sep` is the mean actual minimal separation; the last columns are the smallest relative gap of the entry times and the smallest runner-up margin in the cell):', '',
              table(['K', 'R', 'offsets', 'nominal Delta', 'targets', 'interior-certified', 'behaves', 'mean sep', 'min gap', 'min mu'], rows), '']
    if 'D_fam' in raw:
        F = raw['D_fam']
        ms = raw['D_ms']
        nm = len(ms)
        cert, dual, dec, seps = F[:, :nm], F[:, nm:2 * nm], F[:, 2 * nm:3 * nm], F[:, 3 * nm:4 * nm]
        loc = raw['D_local']

        def thr(flag):
            out = np.full(len(F), np.nan)
            for q in range(len(F)):
                bad = np.flatnonzero(flag[q] == 0)
                k = (bad.max() + 1) if bad.size else 0
                if k < nm:
                    out[q] = seps[q, k]
            return out
        tc, td, tdec = thr(cert), thr(dual), thr(dec)
        mm = np.minimum(loc[:, 5], loc[:, 6])
        L += ['## D. The threshold over the local data (K = 2, continuum inner product)', '',
              '%d draws of the local data (the two features\' offsets from the lattice in [-h/2, h/2], amplitudes +-U(0.4, 1), ratio R <= 2.5); for each, the features are moved apart by multiples of h '
              '(m = %d..%d, separations %.1f..%.1f sigma) and the threshold is the smallest separation from which every larger one of the family (a) is certified (a)-(c), (b) satisfies '
              '(a)-(b) alone (the recursion is generic: unique minimisers, signs, margins), or (c) has the *decoupled prediction\'s* entries (the best starting kernel of each feature, in the order of the '
              'hitting times). nan: not within the window.' % (len(F), ms[0], ms[-1], np.nanmin(seps), np.nanmax(seps)), '']
        q = lambda x: [fmt(float(np.nanquantile(x, p))) for p in (0.0, 0.25, 0.5, 0.75, 0.9, 1.0)]
        L += [table(['threshold (sigma)', 'min', '25%', 'median', '75%', '90%', 'max', 'not reached'],
                    [['certified (a)-(c): the path IS the recursion'] + q(tc) + ['%d' % int(np.isnan(tc).sum())],
                     ['dual conditions (a)-(b) alone'] + q(td) + ['%d' % int(np.isnan(td).sum())],
                     ['entries are the decoupled prediction\'s'] + q(tdec) + ['%d' % int(np.isnan(tdec).sum())]]), '']
        bins = [(0.0, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 1.01)]
        rows = []
        for lo, hi in bins:
            m = (mm >= lo) & (mm < hi)
            ok = m.any() and not np.all(np.isnan(tc[m]))
            rows.append(['%g-%g' % (lo, hi), int(m.sum())] + ([fmt(float(np.nanmedian(tc[m]))), fmt(float(np.nanmax(tc[m]))), fmt(float(np.nanmedian(tdec[m]))), fmt(float(np.nanmax(tdec[m])))] if ok else ['-'] * 4))
        L += ['By the decoupled margin min(mu_0, gap) of the local data (at the lattice offsets of the draw):', '',
              table(['min(mu_0, gap)', 'draws', 'certified: median', 'max', 'decoupled entries: median', 'max'], rows), '']
        okd = ~np.isnan(tdec) & (mm > 0)
        gam = lambda d: np.exp(-np.maximum(d - 0.5 * 0.1 / SIG, 0.0) ** 2 / 4.0)                      # Gamma(Delta - h/2) with h = 0.1 = 0.83 widths
        qd = mm[okd] / gam(tdec[okd])
        L += ['For the decoupled-prediction threshold, min(mu_0, gap) / Gamma(Delta_dec - h/2) has quantiles 10%% / 50%% / 90%% = %.2f / %.2f / %.2f over the %d draws that reach it '
              '(the repository\'s empirical gate of `exp_separated_regime` is 1).' % (*[float(np.quantile(qd, p)) for p in (0.1, 0.5, 0.9)], int(okd.sum())), '']
    if 'E_res' in raw:
        R_, M_ = raw['E_res'], raw['E_meta']
        L += ['## E. The trained network (P0) on K = 2 targets: the largest eps below which the formed set is the path\'s', '',
              '%d targets per separation (seeded), eps = %s. For a target, eps_max is the largest grid eps such that the formed set equals the path\'s (and K units sit on the '
              'features) at it and at every smaller grid eps; "none" if it fails at 1e-10. mu is the path\'s runner-up margin (all later entrants included).' % (
                  R_.shape[1], ', '.join('%g' % e for e in E_EPS)), '']
        rows = []
        ok = (R_[..., 0] > 0) & (R_[..., 1] > 0)
        for a, sep in enumerate(E_SEPS):
            emax = np.full(R_.shape[1], np.nan)
            for i in range(R_.shape[1]):
                k = np.flatnonzero(~ok[a, i])
                kmax = k.max() if k.size else -1
                emax[i] = E_EPS[kmax + 1] if kmax + 1 < len(E_EPS) else np.nan
            cer = M_[a, :, 0] > 0
            r1 = lambda m: (fmt(float(np.nanmedian(emax[m])), '%.1e') if m.any() and not np.all(np.isnan(emax[m])) else '-')
            r2 = lambda m: (fmt(float(np.nanmin(emax[m])), '%.1e') if m.any() and not np.all(np.isnan(emax[m])) else '-')
            dd = M_[a, :, 2] / np.sqrt(emax)
            rows.append(['%g' % sep, R_.shape[1], '%d' % int(cer.sum()), '%d' % int(np.isnan(emax).sum()), r1(np.ones(len(emax), bool)), r2(np.ones(len(emax), bool)),
                         r1(cer), r2(cer), fmt(float(np.nanmedian(dd)), '%.1f') if np.any(~np.isnan(dd)) else '-', fmt(float(np.nanmax(dd)), '%.1f') if np.any(~np.isnan(dd)) else '-'])
        L += [table(['Delta', 'targets', 'continuum-certified', 'none', 'eps_max median', 'eps_max min', 'certified: median', 'certified: min', 'mu / eps_max^(1/2): median', 'max'], rows), '']
        # mu / sqrt(eps_max): the "mismatch margin in units of eps^(1/2)" per target
    if 'F_P0' in raw:
        eps = raw['F_eps']
        P0, D3, mar = raw['F_P0'], raw['F_D3'], raw['F_margin']
        n = len(P0)
        L += ['## F. The mismatch margin, per target (200 targets of Part C of `exp_parametrization`)', '',
              'P0 runs at eps = %s (stored at 1e-6 and 1e-8, new at the others for the %d targets with margin below %g; the others are taken as set = path, as in the stored runs). '
              'For a target, eps_c is the largest eps such that the formed set is the path\'s at it and at every smaller eps of the grid.' % (
                  ', '.join('%g' % e for e in eps), int(((mar[:, 0] < F_MAXMARGIN) | np.isnan(mar[:, 0])).sum()), F_MAXMARGIN), '']
        ok = P0 > 0
        mu = mar[:, 0]
        # largest mismatch margin at each eps (the paper's quantity)
        rows = []
        for b, e in enumerate(eps):
            mis = (P0[:, b] == 0)
            mx = float(np.nanmax(mu[mis])) if mis.any() else np.nan
            rows.append(['%g' % e, '%d' % int(mis.sum()), fmt(mx, '%.4f'), fmt(mx / np.sqrt(e), '%.1f') if mis.any() else '-'])
        L += [table(['eps', 'targets with set != path', 'largest margin of a mismatch', 'in units of eps^(1/2)'], rows), '',
              'The largest mismatch margin at eps = %g (%g) is close to the cut %g, below which the new runs were made: it is a lower bound.' % (eps[0], float(np.nanmax(mu[P0[:, 0] == 0])), F_MAXMARGIN), '']
        mis_eps = [(e, float(np.nanmax(mu[P0[:, b] == 0]))) for b, e in enumerate(eps) if (P0[:, b] == 0).any() and e <= 1e-6]
        if len(mis_eps) >= 2:
            x, y = np.log(np.array([m[0] for m in mis_eps])), np.log(np.array([m[1] for m in mis_eps]))
            L += ['Log-log slope of the largest mismatch margin against eps over the %d values eps <= 1e-6 with a mismatch (%s): %.3f; between the stored values 1e-6 and 1e-8 alone: %.3f '
                  '(a margin proportional to eps^(1/2) would give 0.5).' % (len(mis_eps), ', '.join('%g' % m[0] for m in mis_eps), float(np.polyfit(x, y, 1)[0]),
                                                                           np.log(dict(mis_eps)[1e-6] / dict(mis_eps)[1e-8]) / np.log(100.0)), '']
        # the collapse test: the fraction of mismatches as a function of x = mu / eps^(1/2)
        edges = [0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, np.inf]
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            row = ['%g-%g' % (lo, hi) if np.isfinite(hi) else '>%g' % lo]
            for b, e in enumerate(eps):
                x = mu / np.sqrt(e)
                m = (x >= lo) & (x < hi)
                if e == eps[0] and np.isfinite(hi) and hi * np.sqrt(e) > F_MAXMARGIN or (e == eps[0] and not np.isfinite(hi)):
                    row.append('-')                                              # the bin contains targets above the cut: not run at this eps
                else:
                    row.append('%d / %d' % (int((P0[m, b] == 0).sum()), int(m.sum())) if m.any() else '-')
            rows.append(row)
        L += ['If the mismatch margin were proportional to eps^(1/2), the fraction of mismatches would depend on x = mu / eps^(1/2) only, the same in every column. Mismatches / targets by x:', '',
              table(['x = mu / eps^(1/2)'] + ['%g' % e for e in eps], rows), '']
        # per target: eps_c = the largest eps of the grid below which the set is always the path's; the true threshold lies between eps_c and the next larger grid value
        epsc = np.full(n, np.nan)
        for t in range(n):
            k = np.flatnonzero(~ok[t])
            if not k.size:
                continue                                                          # never a mismatch on the grid: censored above (eps_c >= the grid's top)
            kmax = k.max()
            epsc[t] = eps[kmax + 1] if kmax + 1 < len(eps) else 0.0
        fin = np.isfinite(epsc) & (epsc > 0) & (mu > 0)
        L += ['Per target: %d targets have a mismatch somewhere on the grid and a threshold inside it (eps_c: the largest grid eps from which every smaller one has the path\'s set; the true '
              'threshold lies between eps_c and the next larger grid value, so eps_t = eps_c * sqrt(10) is used); %d targets never mismatch on the grid (censored above), %d mismatch at %g (censored below).' % (
                  int(fin.sum()), int(np.isnan(epsc).sum()), int((epsc == 0).sum()), eps.min()), '']
        if fin.sum() >= 3:
            et = epsc[fin] * np.sqrt(10.0)
            x_, y_ = np.log(et), np.log(mu[fin])
            sl, ic = np.polyfit(x_, y_, 1)
            res = y_ - (sl * x_ + ic)
            d_t = mu[fin] / np.sqrt(et)
            L += ['- log-log slope of the runner-up margin mu on eps_t: %.3f (eps^(1/2): 0.5), standard deviation of the log10 residual %.2f decades; '
                  'D_t = mu / eps_t^(1/2) has quantiles 10%% / 50%% / 90%% = %.1f / %.1f / %.1f and range %.2f - %.0f.' % (
                      sl, float(np.std(res / np.log(10.0))), *[float(np.quantile(d_t, p)) for p in (0.1, 0.5, 0.9)], float(d_t.min()), float(d_t.max())), '']
            ex = np.abs(res) <= 1e9
            hi_ = mu[fin] >= np.median(mu[fin])
            sl_hi = np.polyfit(x_[hi_], y_[hi_], 1)[0]
            sl_lo = np.polyfit(x_[~hi_], y_[~hi_], 1)[0]
            L += ['- the same slope for the half of these targets with the larger margins: %.3f, with the smaller margins: %.3f.' % (sl_hi, sl_lo), '']
    return '\n'.join(L) + '\n'


# ------------------------------------------------------------------------------------------------------------------------------------
def load_raw():
    with np.load(RAW) as z:
        return {k: z[k] for k in z.files}


def save_raw(raw):
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(RAW, **raw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--parts', default='A,B,C,D,E,F,G,H,W,I')
    a = ap.parse_args()
    if a.run:
        warnings.filterwarnings('ignore')
        raw = load_raw() if os.path.exists(RAW) else {}
        raw['have_arb'] = np.array([float(HAVE_ARB)])
        for p in a.parts.split(','):
            t0 = time.time()
            raw.update({'A': part_A, 'B': part_B, 'C': part_C, 'D': part_D, 'E': part_E, 'F': part_F, 'W': part_W, 'G': part_G, 'H': part_H, 'I': part_I}[p]())
            save_raw(raw)
            print('part', p, 'done in %.0f s' % (time.time() - t0), flush=True)
    if a.run or a.from_raw:
        with open(MD, 'w') as f:
            f.write(report(load_raw()))
        print('wrote', MD)
    else:
        ap.error('give --run or --from-raw')


if __name__ == '__main__':
    main()
