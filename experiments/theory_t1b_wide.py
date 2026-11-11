"""
Theorem T1, widened: general K, the planar case, interval arithmetic, monotonicity in eps  (companion of THEORY_T1.md, section 12)
====================================================================================================================================

Reuses `exp_theory_t1` (the certificate `cert`, the limit recursion, the integrated network, the box constants) without modifying it.

  BK. Delta_T1(K, R, d_0) for K = 2 ... 64: the regime (box conditions (L), (C), (N) and the two Gershgorin margins) in ball arithmetic.
  K.  The protocol of Part C of exp_theory_t1 for K = 4, 5, 6 (15 targets, 1-D, continuum): cert against the integrated network on the half-decade grid, the ingredient
      inequalities of the proof on the trajectories at eps_0 and eps_0/10 (Part D protocol), the measured threshold.
  (P, Q, M: added below as they are completed: planar certificate, ball-arithmetic evaluation, monotonicity.)

Usage: python experiments/theory_t1b_wide.py --run [--quick] [--parts BK,K,...] | --from-raw   -> results/synthetic/theory_t1b.md (+ theory_t1b_raw.npz); environment TT1B_OUT.
A partial run (--parts) updates the stored raw.  Deterministic; the report is regenerated from the raw byte for byte.
"""
import argparse
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import exp_theory_t1 as T      # noqa: E402
import theory_t1b_planar as PL  # noqa: E402
import theory_t1b_rigorous as RG  # noqa: E402
from flint import arb  # noqa: E402

OUT = os.environ.get('TT1B_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory_t1b_raw.npz'), os.path.join(OUT, 'theory_t1b.md')

# the heavy sup/oscillation of the limit's rate function does not depend on eps: memoised
_orig_rgs = T.resid_grad_sup
_RG = {}


def _cached_rgs(prep, stage, c_lo, c_hi, n=401):
    k = (id(prep), stage, round(float(c_lo), 12), round(float(c_hi), 12), n)
    if k not in _RG:
        _RG[k] = _orig_rgs(prep, stage, c_lo, c_hi, n)
    return _RG[k]


T.resid_grad_sup = _cached_rgs

K_TARGETS = [(4, 8.0, 1.5, 1), (4, 8.0, 1.5, 2), (4, 10.0, 1.0, 4), (4, 8.0, 2.4, 5), (4, 8.0, 1.5, 3),
             (5, 8.0, 1.5, 1), (5, 8.0, 1.5, 2), (5, 10.0, 2.0, 4), (5, 8.0, 1.0, 5), (5, 8.0, 1.5, 6),
             (6, 8.0, 1.5, 2), (6, 10.0, 1.0, 4), (6, 8.0, 2.4, 5), (6, 8.0, 1.5, 6), (6, 8.0, 1.5, 3)]
K_QUICK = [(4, 8.0, 1.5, 1)]
BK_KS, BK_RS = (2, 3, 4, 5, 6, 8, 10, 16, 32, 64), (1.0, 2.4, 10.0)


def make_target_k(K, sep, R, seed):
    """As exp_theory_t1.make_target for K <= 3 (identical draws); for K >= 4 amplitudes R^(i/(K-1)) times 0.7 in random order, random signs."""
    if K <= 3:
        return T.make_target(K, sep, R, seed)
    rng = np.random.default_rng(seed)
    cs = np.arange(K) * sep + rng.uniform(-T.H_LAT / 2, T.H_LAT / 2, K)
    cs = cs - cs.mean()
    a = R ** (np.arange(K) / (K - 1.0))
    a = a[rng.permutation(K)] * 0.7
    As = a * rng.choice([-1.0, 1.0], K)
    c0 = np.arange(cs.min() - 8, cs.max() + 8, T.H_LAT) + rng.uniform(0, T.H_LAT)
    return cs, As, c0


def part_BK():
    rows = []
    for K in BK_KS:
        for R in BK_RS:
            rows.append([K, R, T.delta_t1(K, R, 0.4167, 0.02), T.delta_t1(K, R, 0.2, 0.02)])
    return dict(BK_rows=np.array(rows, float))


def part_K(quick=False):
    warnings.filterwarnings('ignore')
    tg = K_QUICK if quick else K_TARGETS
    ex = T.EX_QUICK if quick else T.EX_GRID
    meta, OKC, OKM, DD, COMP, FINE, DROWS = [], [], [], [], [], [], []
    for (K, sep, R, seed) in tg:
        cs, As, c0 = make_target_k(K, sep, R, seed)
        prep = T.prepare(cs, As, c0)
        okc, okm, okx, Dd, Pi = T.scan_target(prep, ex)
        e0 = T.thr_first(okc, ex)
        comp = np.full(10, np.nan)
        drow = [np.nan] * 22
        fine = [np.nan] * 3
        if np.isfinite(e0):
            o0, pi0 = T.best_cert(prep, 10.0 ** e0)
            comp = np.r_[o0['comp'], o0['D'], max(max(o0['Sp']), max(o0['Sm'])) / o0['lam'], o0['Tlay'], float(pi0)]
            fine = [float(T.best_cert(prep, 10.0 ** (e0 - dd))[1] >= 0) for dd in (0.125, 0.25, 0.375)]
            for dexp in (0.0, -1.0):
                o1, p1 = T.best_cert(prep, 10.0 ** (e0 + dexp))
                if p1 >= 0:
                    row = T.diagnostics(prep, 10.0 ** (e0 + dexp), T.PARS[p1])
                    if row is not None:
                        DROWS.append([K, sep, e0 + dexp] + row)
        COMP.append(comp)
        FINE.append(fine)
        # the actual dual deviation at the measured threshold (Part H protocol)
        em = T.thr_first(okm, ex)
        eact = np.nan
        if np.isfinite(em):
            sol, U = T.simulate(prep, 10.0 ** em, nlocal=1500)
            eact = T.dual_dev(prep, 10.0 ** em, sol, U)
        meta.append([K, sep, R, seed, prep['Delta'], prep['R'], prep['mu_min'], prep['gap_min'], max(prep['d0']), len(c0),
                     min(e['p_e'] for e in prep['L']['events']), e0, em, eact])
        OKC.append(okc)
        OKM.append(okm)
        DD.append(Dd)
    return dict(K_ex=ex, K_meta=np.array(meta, float), K_okc=np.array(OKC), K_okm=np.array(OKM), K_D=np.array(DD), K_comp=np.array(COMP), K_fine=np.array(FINE),
                K_drows=np.array(DROWS, float) if DROWS else np.zeros((0, 25)))


P_TARGETS = [(1, 0.0, 1.0, 1), (1, 0.0, 1.0, 2), (1, 0.0, 1.0, 3),
             (2, 8.0, 1.0, 1), (2, 8.0, 1.5, 2), (2, 10.0, 2.4, 3),
             (3, 8.0, 1.0, 1), (3, 8.0, 1.5, 2), (3, 10.0, 2.4, 3)]
P_QUICK = [(1, 0.0, 1.0, 1)]
EX_P = np.arange(-3.0, -11.01, -1.0)
EX_PQ = np.array([-4.0, -7.0])


def planar_hessian_check(K, R, Delta, r, n=200, seed=0):
    """Random points of the planar relaxation set D-hat (where the margins are claimed): the exact Hessian of the fit F'(v, d_1, d_2) by finite differences of the analytic gradient; the cross operator-norm row sums
    against row_v, row_d and the block-Gershgorin margins (lambda_min of the diagonal block minus the sum of the operator norms of the off-diagonal blocks) against beta_v, beta_d.
    Returns (max cross row sum v / row_v, max cross row sum d / row_d, min margin v - beta_v, min margin d - beta_d)."""
    b = T.box_t1(K, R, Delta, r, 0.0)
    rng = np.random.default_rng(seed)
    cs = np.array([[Delta * i, 0.37 * i * (-1) ** i] for i in range(K)])
    cs[:, 1] = rng.uniform(-1, 1, K) * Delta * 0.3
    while True:
        D_ = np.linalg.norm(cs[:, None] - cs[None], axis=2) + np.eye(K) * 1e9
        if D_.min() >= Delta:
            break
        cs[:, 1] *= 0.9
    a = np.r_[1.0, rng.uniform(1, R, K - 1)] if K > 1 else np.array([1.0])
    As = a * np.where(rng.random(K) < 0.5, -1, 1)
    Dl = min(D_.min(), 1e9)

    def grad(y):
        v, d = y[:K], y[K:].reshape(K, 2)
        c = cs + d
        w = v * np.sign(As)
        Dd = c[:, None] - cs[None]
        E = c[:, None] - c[None]
        gD, gE = PL.gam2(PL._nrm(Dd)), PL.gam2(PL._nrm(E))
        g = gD @ As - gE @ w
        gradE = -(E / 2) * gE[..., None]
        gradE[np.arange(K), np.arange(K)] = 0.0
        h = np.einsum('ikd,k->id', -(Dd / 2) * gD[..., None], As) - np.einsum('ijd,j->id', gradE, w)
        return np.r_[-np.sign(As) * g, (-w[:, None] * h).ravel()]
    ra = rb = 0.0
    mv = md = np.inf
    for _ in range(n):
        ang = rng.uniform(0, 2 * np.pi, K)
        rad = b['dhat'] * np.sqrt(rng.uniform(0, 1, K))
        y = np.r_[rng.uniform(b['vD_lo'], b['vD_hi'], K), (rad[:, None] * np.c_[np.cos(ang), np.sin(ang)]).ravel()]
        N = len(y)
        H = np.zeros((N, N))
        for j in range(N):
            e = np.zeros(N)
            e[j] = 1e-6
            H[:, j] = (grad(y + e) - grad(y - e)) / 2e-6
        H = 0.5 * (H + H.T)
        blocks = [[i] for i in range(K)] + [[K + 2 * i, K + 2 * i + 1] for i in range(K)]
        own = np.zeros_like(H)
        for i in range(K):
            ai, di, vi = abs(As[i]), y[K + 2 * i:K + 2 * i + 2], y[i]
            dn = np.linalg.norm(di)
            s_ = ai * di / 2 * np.exp(-dn ** 2 / 4)
            own[i, i] = 1
            own[i, K + 2 * i:K + 2 * i + 2] = s_
            own[K + 2 * i:K + 2 * i + 2, i] = s_
            own[K + 2 * i:K + 2 * i + 2, K + 2 * i:K + 2 * i + 2] = vi * ai * (np.eye(2) / 2 - np.outer(di, di) / 4) * np.exp(-dn ** 2 / 4)
        C = H - own
        for bi, blk in enumerate(blocks):
            cross = sum(np.linalg.norm(C[np.ix_(blk, o)], 2) for oi, o in enumerate(blocks) if oi != bi)
            if bi < K:
                rv = cross
                ra = max(ra, rv / max(b['rowv'], 1e-300))
            else:
                rb = max(rb, cross / max(b['rowd'], 1e-300))
            lam = np.linalg.eigvalsh(H[np.ix_(blk, blk)]).min()
            off = sum(np.linalg.norm(H[np.ix_(blk, o)], 2) for oi, o in enumerate(blocks) if oi != bi)
            if bi < K:
                mv = min(mv, lam - off)
            else:
                md = min(md, lam - off)
    return ra, rb, mv - b['beta_v'], md - b['beta_d'], Dl


def part_P(quick=False):
    warnings.filterwarnings('ignore')
    tg = P_QUICK if quick else P_TARGETS
    ex = EX_PQ if quick else EX_P
    meta, OKC, OKM, DD, DROWS = [], [], [], [], []
    for (K, sep, R, seed) in tg:
        cs, As, c0 = PL.make_target2(K, sep, R, seed)
        prep = PL.prepare2(cs, As, c0)
        okc, okm, Dd = [], [], []
        for x_ in ex:
            eps = 10.0 ** x_
            o, pi = PL.best_cert2(prep, eps)
            okc.append(float(o['ok']))
            Dd.append(o['D'] if o['ok'] else np.nan)
            sol, U = PL.simulate2(prep, eps, nlocal=300)
            okm.append(float(PL.follows2(prep, eps, sol, U)[0]))
        e0, em = T.thr_first(okc, ex), T.thr_first(okm, ex)
        if np.isfinite(e0):
            for dexp in (0.0, -1.0):
                o1, p1 = PL.best_cert2(prep, 10.0 ** (e0 + dexp))
                if p1 >= 0:
                    row = PL.diagnostics2(prep, 10.0 ** (e0 + dexp), T.PARS[p1])
                    if row is not None:
                        DROWS.append([K, sep, e0 + dexp] + row)
        meta.append([K, sep, R, seed, prep['Delta'], prep['R'], prep['mu_min'], prep['gap_min'], max(prep['d0']), len(c0), e0, em])
        OKC.append(okc)
        OKM.append(okm)
        DD.append(Dd)
    hs = [planar_hessian_check(K, R, D, r, 60 if quick else 200) for (K, R, D, r) in ((2, 1.4, 6.2, 0.6), (3, 1.5, 6.5, 0.6), (2, 2.4, 6.5, 0.5))]
    return dict(P_ex=ex, P_meta=np.array(meta, float), P_okc=np.array(OKC), P_okm=np.array(OKM), P_D=np.array(DD), P_drows=np.array(DROWS, float) if DROWS else np.zeros((0, 17)),
                P_hess=np.array(hs, float), P_delta=np.array([[R, d0, T.delta_t1(2, R, d0, 0.02), T.delta_t1(3, R, d0, 0.02)] for R in (1.0, 2.4, 10.0) for d0 in (0.2, 0.4167, 0.59)], float))


def part_U():
    """Lemma 12.4 (eps-uniform envelope of the dormant coefficient): B_eps(s) <= 4/(G gamma^2 m) (1 - (sigma_0/s)^2)^(-1) s^(-3), sigma_0 = a_max/(2 G v_s) = 20, for every eps and s >= sigma_0
    (Lemma 2.2 with the Adamovic-Mitrinovic inequality (sin t / t)^3 > cos t).  Numerically: max over a grid of s of B s^3/envelope for several eps."""
    rows = []
    for eps in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-16):
        U = T.Unit(eps)
        s = np.geomspace(30.0, 0.999 * U.lam, 4000)
        s = s[s > 25]
        B = U.B(s)
        env = 4 / (T.LAW['G'] * T.LAW['gamma'] ** 2 * T.LAW['m']) / (1 - (20.0 / s) ** 2) / s ** 3
        rows.append([eps, float((B / env).max()), float((B * s ** 3).max()), float(s[np.argmax(B * s ** 3)])])
    t = np.linspace(1e-4, np.pi / 2 - 1e-9, 200001)
    am = float((np.sin(t) ** 3 / (t ** 3 * np.cos(t))).min())       # >= 1: (sin t / t)^3 >= cos t
    return dict(U_rows=np.array(rows, float), U_am=np.array([am]))


Q_IDX = (0, 1, 3, 4, 7, 8, 12, 9, 6, 13)       # indices into exp_theory_t1.C_TARGETS (K <= 3) for the interval-arithmetic cover
Q_BOT, Q_BOT_QUICK = -24.0, -10.0


def cover(prep, e_top, e_bot, wmin=1.0 / 512, budget=None):
    """Interval-arithmetic cover of eps in [10^e_bot, 10^e_top] by balls on which the ball-valid certificate holds (floating-point candidates D_j from `cert` at the upper end of each
    ball, enlarged by the factors 1.1, 1.5, 2, 4; the ball width adapts).  Returns (list of (lo, hi), complete, lowest exponent reached, number of verify calls)."""
    cur, segs, calls = e_top, [], 0
    t_end = None if budget is None else time.time() + budget
    wlast = 0.5
    ent = [e['kernel'] for e in prep['L']['events']]
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
            ball = arb(10.0 ** lo_).union(arb(10.0 ** cur))
            for fac in (1.1, 1.5, 2.0, 4.0):
                Dv = [1.0 if j in ent else max(fac * o['Etot'][j], 1.0) for j in range(len(prep['c0']))]
                calls += 1
                try:
                    ok, info = RG.verify(prep, ball, par, Dv)
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


Q_BUDGET = 240.0        # seconds per attempt (per target and per top exponent)


def part_Q(quick=False):
    warnings.filterwarnings('ignore')
    craw = T.load_raw()
    meta = craw['C_meta']
    idx = Q_IDX[:1] if quick else Q_IDX
    e_bot = Q_BOT_QUICK if quick else Q_BOT
    rows, segs_all = [], []
    for n_, ti in enumerate(idx):
        K, sep, R, seed = T.C_TARGETS[ti]
        cs, As, c0 = T.make_target(K, sep, R, seed)
        prep = T.prepare(cs, As, c0)
        e0 = meta[ti, 11]
        best = None
        for dtop in (0.0, -0.5, -1.0):
            segs, complete, cur, calls = cover(prep, e0 + dtop, e_bot, budget=None if quick else Q_BUDGET)
            if complete:
                best = (e0 + dtop, len(segs), calls, segs, cur, 1.0)
                break
            if best is None or cur < best[4]:
                best = (e0 + dtop, len(segs), calls, segs, cur, 0.0)
        top, nb, calls, segs, cur, comp = best
        wmax = max([w_[1] - w_[0] for w_ in segs]) if segs else 0.0
        rows.append([K, sep, R, seed, e0, top, nb, calls, wmax, cur, comp])
        segs_all.append(np.array([[n_, w_[0], w_[1]] for w_ in segs], float).reshape(-1, 3))
        print('Q target', n_, (K, sep, R, seed), 'top', top, 'lowest', cur, 'complete', comp, flush=True)
    return dict(Q_rows=np.array(rows, float), Q_bot=np.array([e_bot]), Q_segs=np.concatenate(segs_all) if segs_all else np.zeros((0, 3)))


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(str(c) for c in r) + ' |')
    return '\n'.join(out)


def informative(raw, key, i_eps0, i_max):
    """(pairs with S true (cert), pairs scanned, violations among them, whether every informative pair lies at eps <= eps_0 of its target, smallest log10(eps_max/eps_0) over the targets)"""
    okc, okm, M, ex = raw[key + '_okc'], raw[key + '_okm'], raw[key + '_meta'], raw[key + '_ex']
    inf = okc > 0
    below = all((not inf[i].any()) or ex[inf[i]].max() <= M[i, i_eps0] + 1e-9 for i in range(len(M)) if np.isfinite(M[i, i_eps0]))
    gaps = M[:, i_max] - M[:, i_eps0]
    return int(inf.sum()), int(okc.size), int((inf & (okm == 0)).sum()), bool(below), float(np.nanmin(gaps))


def fe(ex):
    return 'none' if not np.isfinite(ex) else '1e%+.1f' % ex


PARTS_ORDER = ('BK', 'K', 'P', 'Q', 'U')


BANNER = 'SUPERSEDED by results/synthetic/theory3n.md: the displacement accounting of the certificate was corrected (see review/REVIEW_roundT3AC.md, REVIEW_roundT3N.md); numbers below use the legacy accounting and are kept as history'


def report(raw):
    L = []
    P = L.append
    P(BANNER + '\n')
    P('# Theorem T1, widened (generated by experiments/theory_t1b_wide.py)\n')
    P('Companion of THEORY_T1.md section 12.  Same certificate `cert` as `exp_theory_t1.py` (not modified); %s.\n' % ('ball arithmetic (python-flint) for the box constants' if bool(raw.get('HAVE_ARB', [0])[0]) else 'floats'))
    if 'BK_rows' in raw:
        P('## BK. The regime against the number of features\n')
        P('Smallest separation (grid 0.1 sigma, dormant budget Omega = 0.02 a_min, ball arithmetic) with a radius r > d_0 for which (L), (C), (N) and the two Gershgorin margins hold; '
          'K enters through kappa_c = (K-1) Gamma(Delta - 2r) and the row sums (proportional to K-1).\n')
        rows = [['%d' % r[0], '%g' % r[1], '%.1f' % r[2], '%.1f' % r[3]] for r in raw['BK_rows']]
        P(table(['K', 'R', 'Delta_T1 (d_0 = 0.4167)', 'Delta_T1 (d_0 = 0.2)'], rows) + '\n')
    if 'K_meta' in raw:
        M, ex = raw['K_meta'], raw['K_ex']
        P('## K. The certificate against the integrated network for K = 4, 5, 6 (%d targets, eps = 10^%.1f ... 10^%.1f)\n' % (len(M), ex[0], ex[-1]))
        P('Columns as in Part C of theory_t1.md; E_act = the largest dual deviation of a dormant kernel of the integrated network at the measured eps_max.\n')
        rows, viol = [], 0
        for i, m in enumerate(M):
            v = int(((raw['K_okc'][i] > 0) & (raw['K_okm'][i] == 0)).sum())
            viol += v
            ratio = 10 ** (m[12] - m[11]) if np.isfinite(m[11]) and np.isfinite(m[12]) else np.nan
            Dd = raw['K_D'][i]
            D0 = Dd[np.flatnonzero(ex == m[11])[0]] if np.isfinite(m[11]) else np.nan
            rows.append(['%d' % m[0], '%g' % m[1], '%.2f' % m[5], '%d' % m[3], '%.2f' % m[4], '%.4f' % m[6], '%.3f' % m[7], '%.2f' % m[8], '%d' % m[9], fe(m[11]), fe(m[12]),
                         '%.0f' % ratio if np.isfinite(ratio) else '-', '%.0f' % D0 if np.isfinite(D0) else '-', '%.1f' % m[13] if np.isfinite(m[13]) else '-', '%d' % v])
        P(table(['K', 'sep', 'R', 'seed', 'Delta', 'mu', 'gap', 'd_0', 'kernels', 'eps_0 (cert)', 'eps_max (measured)', 'ratio', 'D at eps_0', 'E_act', 'violations'], rows) + '\n')
        C_ = raw['K_comp']
        rows = []
        for i, m in enumerate(M):
            c = C_[i]
            if not np.isfinite(c[0]):
                continue
            pr = T.PARS[int(c[9])]
            rows.append(['%d' % m[0], '%d' % m[3], '%.1f' % c[0], '%.1f' % c[1], '%.1f' % c[2], '%.1f' % c[3], '%.2g' % c[4], '%.2g' % c[5], '%.1f' % c[6], '%.2e' % c[7], '%.0f' % c[8],
                         '(%g, %g, %g)' % (pr['r'], pr['Om'], pr['delta1'])])
        P('Where D comes from at eps_0 (largest kernel): own windows, own tails, old units, old tails, mass, drift; S/lambda; T_layer; parameters (r, Omega, delta_1).\n')
        P(table(['K', 'seed', 'own', 'own tails', 'old', 'old tails', 'mass', 'drift', 'D', 'S/lambda', 'T_layer', 'params'], rows) + '\n')
        okr = np.array([10 ** (m[12] - m[11]) for m in M if np.isfinite(m[11]) and np.isfinite(m[12])])
        Fi = raw['K_fine']
        ni, nt, nv, below, gmin = informative(raw, 'K', 11, 12)
        P('Violations: **%d** among the %d informative (target, eps) pairs (S true, i.e. `cert` true; a pair with `cert` false cannot be a violation) out of %d pairs scanned.  %s'
          '  The check therefore has no power near the boundary: it compares the formed set at the final time only (the ingredient inequalities below cover the rest, on few trajectories).  '
          'Refinement below eps_0 (eps_0/10^(1/8), /10^(1/4), /10^(3/8)): %d of %d pairs hold.'
          % (nv, ni, nt, ('Every informative pair lies at eps <= eps_0 of its target, i.e. at least %.1f decades below the measured threshold (smallest log10 eps_max/eps_0 over the targets, from the table above).' % gmin) if below
             else 'NOT every informative pair lies at eps <= eps_0 of its target (smallest log10 eps_max/eps_0 over the targets: %.1f).' % gmin, int(np.nansum(Fi)), int(np.isfinite(Fi).sum())))
        if len(okr):
            P('Ratio eps_max/eps_0: median %.0f, range %.0f-%.0f.  eps_0 over the targets: 10^%.1f to 10^%.1f; measured eps_max: 10^%.1f to 10^%.1f.\n'
              % (np.median(okr), okr.min(), okr.max(), np.nanmin(M[:, 11]), np.nanmax(M[:, 11]), np.nanmin(M[:, 12]), np.nanmax(M[:, 12])))
        R_ = raw['K_drows']
        if len(R_):
            idx = {c: i for i, c in enumerate(T.D_COLS)}
            n = len(R_)
            P('Ingredients on %d trajectories (eps_0 and eps_0/10): dual deviation / bound max ratio **%.3f** (max measured E = %.1f, bound D up to %.0f); entry shift in [-S^-, S^+]: %d of %d; transit <= T_1: %d of %d (max ratio %.2f); '
              'relaxation <= T_2 + T_3: %d of %d (max ratio %.2f); entrant lead <= D^e_k: %d of %d (max ratio %.2f); monotone entrants: %d of %d; offsets <= r: %d of %d; coefficients <= v_bar: %d of %d; '
              'dormant mass ratio max %.4f; formed set = limit\'s: %d of %d; final coefficient error <= %.1e.\n'
              % (n, R_[:, idx['E_over_Dv']].max(), R_[:, idx['E_max']].max(), R_[:, idx['D']].max(), int(R_[:, idx['shift_ok']].sum()), n, int(R_[:, idx['T1_ok']].sum()), n, R_[:, idx['T1_ratio']].max(),
                 int(R_[:, idx['rel_ok']].sum()), n, R_[:, idx['rel_ratio']].max(), int(R_[:, idx['Ee_ok']].sum()), n, R_[:, idx['Ee_ratio']].max(), int(R_[:, idx['mono']].sum()), n,
                 int(R_[:, idx['d_ok']].sum()), n, int(R_[:, idx['w_ok']].sum()), n, (R_[:, idx['omega_act']] / np.maximum(R_[:, idx['omega_bound']], 1e-300)).max(),
                 int(R_[:, idx['formed_ok']].sum()), n, R_[:, idx['final_err']].max()))
    if 'P_meta' in raw:
        M, ex = raw['P_meta'], raw['P_ex']
        P('## P. The planar (n = 2) certificate against the integrated planar network (%d targets, eps = 10^%.1f ... 10^%.1f, decades)\n' % (len(M), ex[0], ex[-1]))
        P('Planar box and regime: the constants of `box_t1` are n-independent (THEORY_T1.md 12.2); Delta_T1 for the planar offsets d_0 up to h/sqrt2 = 0.59 (lattice spacing 0.83 sigma):\n')
        rows = [['%g' % r[0], '%.4g' % r[1], '%.1f' % r[2], '%.1f' % r[3]] for r in raw['P_delta']]
        P(table(['R', 'd_0', 'Delta_T1 (K=2)', 'Delta_T1 (K=3)'], rows) + '\n')
        pd_ = raw['P_delta']
        P('At the planar offset d_0 = %.2f the separation Delta_T1 runs from %.1f sigma to %.1f sigma over the K = 2, 3 rows above (R = 1 ... 10), against %.1f to %.1f sigma at d_0 = %.4g.\n'
          % (pd_[:, 1].max(), pd_[pd_[:, 1] == pd_[:, 1].max()][:, 2:].min(), pd_[pd_[:, 1] == pd_[:, 1].max()][:, 2:].max(), pd_[pd_[:, 1] == 0.4167][:, 2:].min(), pd_[pd_[:, 1] == 0.4167][:, 2:].max(), 0.4167))
        h = raw['P_hess']
        P('Planar Hessian check (finite differences of the analytic gradient, random points of the box, 3K x 3K): max (cross operator-norm row sum)/(row_v) = %.2f, /(row_d) = %.2f; '
          'smallest block-Gershgorin margin minus beta_v = %.3f, minus beta_d = %.3f (>= 0 means the bounds hold) over %d cases.\n' % (h[:, 0].max(), h[:, 1].max(), h[:, 2].min(), h[:, 3].min(), len(h)))
        rows, viol = [], 0
        for i, m in enumerate(M):
            v = int(((raw['P_okc'][i] > 0) & (raw['P_okm'][i] == 0)).sum())
            viol += v
            ratio = 10 ** (m[11] - m[10]) if np.isfinite(m[10]) and np.isfinite(m[11]) else np.nan
            Dd = raw['P_D'][i]
            D0 = Dd[np.flatnonzero(ex == m[10])[0]] if np.isfinite(m[10]) else np.nan
            rows.append(['%d' % m[0], '%g' % m[1] if m[0] > 1 else '-', '%.2f' % m[5], '%d' % m[3], '%.2f' % m[4] if m[0] > 1 else '-', '%.4f' % m[6], '%.3f' % m[7], '%.2f' % m[8], '%d' % m[9],
                         fe(m[10]), fe(m[11]), '%.0f' % ratio if np.isfinite(ratio) else '-', '%.0f' % D0 if np.isfinite(D0) else '-', '%d' % v])
        P(table(['K', 'sep', 'R', 'seed', 'Delta', 'mu', 'gap', 'd_0', 'kernels', 'eps_0 (cert)', 'eps_max (measured)', 'ratio', 'D at eps_0', 'violations'], rows) + '\n')
        R_ = raw['P_drows']
        ni, nt, nv, below, gmin = informative(raw, 'P', 10, 11)
        P('Violations: **%d** among the %d informative (target, eps) pairs (S true) out of %d pairs scanned.  %s'
          '  Dictionary: the square-lattice kernels within 3.5 sigma of a feature (%d to %d kernels, not a full planar window), double precision.\n'
          % (nv, ni, nt, ('Every informative pair lies at eps <= eps_0 of its target, i.e. at least %.1f decades below the measured threshold (smallest log10 eps_max/eps_0 over the targets, from the table above).' % gmin) if below
             else 'NOT every informative pair lies at eps <= eps_0 of its target (smallest log10 eps_max/eps_0 over the targets: %.1f).' % gmin, int(M[:, 9].min()), int(M[:, 9].max())))
        if len(R_):
            n = len(R_)
            P('Ingredients on %d planar trajectories (eps_0 and eps_0/10): dual deviation / bound max ratio **%.3f**; entry shift in [-S^-, S^+]: %d of %d; transit <= T_1: %d of %d (max ratio %.2f); '
              'entrant lead <= D^e_k: %d of %d (max ratio %.2f); monotone entrants: %d of %d; offsets <= r: %d of %d; coefficients <= v_bar: %d of %d; dormant mass / bound max %.3f; formed set = limit\'s: %d of %d; final coefficient error <= %.1e.\n'
              % (n, R_[:, 5].max(), int(R_[:, 6].sum()), n, int(R_[:, 7].sum()), n, R_[:, 8].max(), int(R_[:, 9].sum()), n, R_[:, 10].max(), int(R_[:, 11].sum()), n, int(R_[:, 12].sum()), n, int(R_[:, 13].sum()), n,
                 R_[:, 14].max(), int(R_[:, 15].sum()), n, R_[:, 16].max()))
    if 'Q_rows' in raw:
        P('## Q. Interval-arithmetic cover of eps (python-flint balls; K <= 3, n = 1)\n')
        P('For each target the certificate S(eps) was verified in ball arithmetic on a chain of balls covering **every real eps in [10^%.0f, 10^top]** (top = the highest exponent from which the chain closes; '
          '"float eps_0" is the half-decade grid value of theory_t1.md).  Candidate D_j from the floating-point certificate (x1.1 ... x4), then every ingredient re-evaluated with outward-rounded balls '
          '(see theory_t1b_rigorous.py).  Inputs taken as exact floating-point numbers: the limit object\'s data (critical points, entry times, duals at the breakpoints, rates, margins, the entrant/feature assignment); '
          'the final comparisons use the float upper end of a ball; the hypotheses r > d_0, the block critical point in the relaxation set D-hat and theta < v_# are checked against ball constants with floating-point limit data '
          '(asserted in the current `verify`; see theory_t1b_reverify.md for the rerun of these covers with the amended verifier).  Statement: S(eps) evaluated in interval arithmetic, with floating-point limit-object data, for every real eps in the ball.\n' % raw['Q_bot'][0])
        rows = [['%d' % r[0], '%g' % r[1] if r[0] > 1 else '-', '%g' % r[2], '%d' % r[3], '%+.1f' % r[4], '%+.2f' % r[5], '%d' % r[6], '%d' % r[7], '%.3f' % r[8],
                 '%+.2f' % r[9], 'complete' if r[10] > 0 else 'incomplete'] for r in raw['Q_rows']]
        P(table(['K', 'sep', 'R', 'seed', 'float log10 eps_0', 'top (log10)', 'balls', 'verify calls', 'widest ball (decades)', 'lowest exponent reached', 'status'], rows) + '\n')
        nc = int(sum(r[10] > 0 for r in raw['Q_rows']))
        P('%d of %d targets: the certificate holds, in ball arithmetic, for every real eps in [10^%.0f, 10^top] (the chain of balls above; top < float eps_0 when the interval tests fail close to the float threshold).  '
          'For the others the chain stopped (interval test failed at the smallest ball width tried, 1/512 decade, or the time budget of 240 s per attempt ran out): it covers every real eps in [10^lowest, 10^top] only.  '
          'These stops are outcomes of the search budget and of the ball width (the window-start gating of the drift terms becomes undecided on wide balls), not evidence that S fails there (referee round T1b, finding 4).\n' % (nc, len(raw['Q_rows']), raw['Q_bot'][0]))
    if 'U_rows' in raw:
        P('## U. The eps-uniform envelope of the dormant coefficient (Lemma 12.4)\n')
        P('This is the envelope version of the dormant-coefficient bound: it bounds B_eps(s) by an eps-uniform envelope; it does not by itself make the exact mass and drift monotone in eps, and the layer constants remain unsigned.  Valid for lambda >= 25, i.e. eps <~ 0.03 (and eps < 2 v_s).\n')
        P('B_eps(s) <= 4/(G gamma^2 m) (1 - (20/s)^2)^-1 s^-3 for all eps and s >= 25 (Lemma 2.2 + Adamovic-Mitrinovic: min over a 2e5-point grid of (sin t)^3/(t^3 cos t) on (0, pi/2) = %.6f >= 1).\n' % raw['U_am'][0])
        rows = [['%.0e' % r[0], '%.5f' % r[1], '%.4f' % r[2], '%.0f' % r[3]] for r in raw['U_rows']]
        P(table(['eps', 'max B/envelope (s from 25 to lambda)', 'max B s^3', 'at s ='], rows) + '\n')
    return '\n'.join(L)


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
    ap.add_argument('--parts', default=','.join(PARTS_ORDER))
    a = ap.parse_args(argv)
    if a.run:
        raw = {'HAVE_ARB': np.array([float(T.HAVE_ARB)])}
        if a.parts != ap.get_default('parts') and os.path.exists(RAW):
            raw.update(load_raw())
        parts = a.parts.split(',')
        t0 = time.time()
        if 'BK' in parts:
            raw.update(part_BK())
        if 'K' in parts:
            raw.update(part_K(a.quick))
        if 'P' in parts:
            raw.update(part_P(a.quick))
        if 'Q' in parts:
            raw.update(part_Q(a.quick))
        if 'U' in parts:
            raw.update(part_U())
        save_raw(raw)
        print('raw written (%.0f s)' % (time.time() - t0))
    else:
        raw = load_raw()
    with open(MD, 'w') as fh:
        fh.write(report(raw))
    print('report written:', MD)


if __name__ == '__main__':
    main()
