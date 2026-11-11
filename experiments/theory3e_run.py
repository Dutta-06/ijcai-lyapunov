"""
Round theory3 / direction E: the coupled entry layer and the limit path for a GENERAL Gram -- tables, numerical tests, deterministic report.

Companion of NOTE.md of the round (Proposition G, Corollary G1, Proposition G'' for a cusped Gram, Theorem G = Theorem K for a general Gram).  Everything below is a finite computation:

  A. thresholds  Delta_entry(K, R, d_0; Gamma) of the box conditions (L), (C), (N) (ball arithmetic, theory3e_kernels.box) for Gaussian, Laplace-basis (exponential), Matern-3/2- and
                 Matern-5/2-basis Grams in R^n, n = 1, 2, 3 (Gram = Matern-nu_G, nu_G = 2 nu_basis + n/2), the rational-quadratic Gram (algebraic decay) and the Wendland bases (compact Gram);
                 the Gaussian row in the Gaussian note's row-sum form reproduces Table `tab_entry` of the paper (4.66, 5.50, 6.70, 5.08, 5.84);  clean (coherence-form) thresholds of Corollary G1
  B. starting offset d_0 (0 .. 0.8, and the worst offset sqrt(n) h / 2 of a cubic lattice of spacing h)
  C. decay moduli and the clean threshold from them
  D. the cusped Gram exp(-d): conditions of Proposition G'', thresholds; the smoothed cusp integrated (locking of the centres)
  E. anisotropic Gaussian with a common covariance (condition number of the centre mobility)
  F. Proposition G(i)-(iii) on random and adversarial initial points of the box, for four mobility laws (independent of the proof: the inner flow is integrated)
  G. the limit path (recursion) against the network flow at eps = 1e-4, 1e-6, 1e-8 on certified and uncertified targets

Usage: python experiments/theory3e_run.py --run [--quick] [--parts A,B,...] | --from-raw    -> results/synthetic/theory3e.md (+ theory3e_raw.npz).  Environment: T3E_OUT.
Deterministic: every random draw is seeded; the report is regenerated from the raw alone, byte for byte.
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
import theory3e_kernels as TK      # noqa: E402
import theory3e_wendland as TW     # noqa: E402
import theory3e_path as TP         # noqa: E402
import theory3e_t1gen as TG        # noqa: E402

OUT = os.environ.get('T3E_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory3e_raw.npz'), os.path.join(OUT, 'theory3e.md')

D0_BENCH = 0.4167                                   # the benchmark's worst starting offset h / (2 sigma) = 0.05 / 0.12
COMBOS = [(2, 1.0), (2, 2.4), (2, 10.0), (3, 1.0), (3, 2.4)]
ROWS = [  # label, factory
    ('Gaussian (n=1)', lambda: TK.family('gauss', 1)),
    ('Gaussian (n=2)', lambda: TK.family('gauss', 2)),
    ('Gaussian (n=3)', lambda: TK.family('gauss', 3)),
    ('Laplace basis, Gram Matern-3/2 (n=1)', lambda: TK.family('exp', 1)),
    ('Laplace basis, Gram Matern-2 (n=2)', lambda: TK.family('exp', 2)),
    ('Laplace basis, Gram Matern-5/2 (n=3)', lambda: TK.family('exp', 3)),
    ('Matern-3/2 basis, Gram Matern-7/2 (n=1)', lambda: TK.family('mat32', 1)),
    ('Matern-3/2 basis, Gram Matern-4 (n=2)', lambda: TK.family('mat32', 2)),
    ('Matern-3/2 basis, Gram Matern-9/2 (n=3)', lambda: TK.family('mat32', 3)),
    ('Matern-5/2 basis, Gram Matern-11/2 (n=1)', lambda: TK.family('mat52', 1)),
    ('Matern-5/2 basis, Gram Matern-6 (n=2)', lambda: TK.family('mat52', 2)),
    ('Matern-5/2 basis, Gram Matern-13/2 (n=3)', lambda: TK.family('mat52', 3)),
    ('Gram rational-quadratic beta=2 (n=1)', lambda: TK.RationalQuadratic(2, 1)),
    ('Gram rational-quadratic beta=3 (n=1)', lambda: TK.RationalQuadratic(3, 1)),
    ('Gram rational-quadratic beta=5 (n=1)', lambda: TK.RationalQuadratic(5, 1)),
    ('Wendland psi_{1,1} (n=1)', lambda: TW.Wendland(1)),
    ('Wendland psi_{3,1} (n=3)', lambda: TW.Wendland(3)),
]
ROWS_QUICK = [0, 3, 15]
B_ROWS = [0, 3, 6]                                   # Gaussian, Laplace basis n=1, Matern-3/2 basis n=1
B_D0 = [0.0, 0.2, 0.4167, 0.6, 0.8]
E_CONDS = [1.0, 2.0, 4.0, 10.0, 100.0]
F_CASES = [(0, 2, 1.0), (0, 3, 2.4), (3, 2, 2.4), (6, 2, 1.0), (9, 2, 1.0), (13, 2, 2.4), (15, 2, 1.0), (4, 2, 1.0)]       # (row index, K, R)
F_QUICK = [(0, 2, 1.0), (15, 2, 1.0)]
F_OFF = [0.0, -1.0, -2.0]
F_RUNS = 8
G_ROWS = [0, 3, 6, 13, 15]                        # n = 1 kernels of the path/network comparison
G_PLANAR = [4]                                    # Gaussian n=2, Laplace basis n=2
G_KR = [(2, 1.0), (2, 2.4), (3, 1.5)]
G_REL = [0.6, -1.5]
G_EPS = [1e-4, 1e-6, 1e-8]
G_SEED0 = 4000
ONLY = None                                          # restrict the row-indexed parts to these Gram rows (--only); the results are merged into the stored raw
H_ROWS = [0, 3, 6, 15]                               # Gaussian, Laplace basis, Matern-3/2 basis, Wendland psi_{1,1}
H_EX = [-5.0, -6.0, -7.0, -8.0, -10.0]   # log10 eps of the certificate grid
H_NPARS = 4                                          # parameter sets of exp_theory_t1.PARS tried (the repository tries 16: fewer = more conservative)
H_SEED = 9000


def rowdef(i):
    return ROWS[i][1]()


_KC = {}


def kern(i):
    if i not in _KC:
        _KC[i] = rowdef(i)
    return _KC[i]


# ------------------------------------------------------------------------------------------------------------------------------------
def part_A(quick=False):
    rows = ROWS_QUICK if quick else list(range(len(ROWS)))
    if ONLY is not None:
        rows = [i for i in rows if i in ONLY]
    thr = np.full((len(ROWS), len(COMBOS)), np.nan)
    rst = np.full((len(ROWS), len(COMBOS)), np.nan)
    cln = np.full((len(ROWS), len(COMBOS)), np.nan)
    for i in rows:
        k = kern(i)
        for c, (K, R) in enumerate(COMBOS):
            t, r = TK.entry_threshold(k, K, R, D0_BENCH, lo_=1.0, hi=60.0, return_r=True)
            thr[i, c], rst[i, c] = t, (r if r is not None else np.nan)
            cln[i, c] = TK.clean_threshold(k, K, R, D0_BENCH, hi=80.0)
    g = TK.Gaussian(1)
    gauss_rowsum = np.array([TK.entry_threshold(g, K, R, D0_BENCH, lo_=3.0, hi=24.0, hc='rowsum') for (K, R) in COMBOS])
    paper = np.array([4.66, 5.50, 6.70, 5.08, 5.84])
    return dict(A_thr=thr, A_r=rst, A_clean=cln, A_gauss_rowsum=gauss_rowsum, A_paper=paper)


def part_B(quick=False):
    thr = np.full((len(B_ROWS), len(B_D0)), np.nan)
    for a, i in enumerate(B_ROWS[:1] if quick else B_ROWS):
        k = kern(i)
        for b, d0 in enumerate(B_D0[:3] if quick else B_D0):
            thr[a, b] = TK.entry_threshold(k, 2, 1.0, d0, lo_=1.0, hi=60.0)
    # the worst offset of a cubic lattice of spacing h = 0.8333 sigma: sqrt(n) h / 2
    fam = [(1, 0), (2, 1), (3, 2)]
    lat = np.full((len(ROWS), 2), np.nan)
    for i in ([0, 3] if quick else [0, 1, 2, 3, 4, 5, 6, 7, 8]):
        n = kern(i).n
        d0 = np.sqrt(n) * D0_BENCH
        lat[i, 0] = d0
        lat[i, 1] = TK.entry_threshold(kern(i), 2, 1.0, d0, lo_=1.0, hi=60.0)
    return dict(B_thr=thr, B_lat=lat)


def part_C(quick=False):
    out = np.full((len(ROWS), 5), np.nan)             # q, a, C(a,q), c_*(r=0.45, R=1) , Delta_clean-from-modulus (K=2, R=1)
    rows = ROWS_QUICK if quick else list(range(len(ROWS)))
    for i in rows:
        k = kern(i)
        name = ROWS[i][0]
        if name.startswith('Gaussian'):
            q, a = 2.0, 0.24
        elif name.startswith('Gram rational'):
            q, a = 0.0, float(k.beta) * 2          # algebraic: m(l) <= C l^{-2 beta}
        elif name.startswith('Wendland'):
            q, a = 1.0, 0.0
        else:
            q, a = 1.0, 0.9 * k.s
        if name.startswith('Gram rational'):
            b = k.beta
            k.tail(0, 3.0)
            ls = k._tab[0]
            C = 0.0
            for j in range(len(ls) - 1):
                if ls[j] < 3.0:
                    continue
                m = max(k._tab[1][j], k._tab[2][j], k._tab[3][j] if k.n == 1 else k._tab[4][j])
                C = max(C, m * ls[j + 1] ** (2 * b))
        elif name.startswith('Wendland'):
            C = 0.0
        else:
            C = TK.decay_modulus(k, q, a, l0=3.0)
        cst = TK.clean_cstar(k, 1.0, 0.45)
        # Delta from the modulus: 2 r + ell, ell solving (K-1) R^2 C exp(-a ell^q) = c_*  (algebraic: (K-1) R^2 C ell^{-2 beta} = c_*)
        if cst > 0 and C > 0:
            if q == 0.0:
                ell = (C / cst) ** (1.0 / a)
            else:
                ell = (np.log(C / cst) / a) ** (1.0 / q)
            dl = 2 * 0.45 + ell
        elif C == 0.0:
            dl = np.nan
        else:
            dl = np.nan
        out[i] = [q, a, C, cst, dl]
    return dict(C_mod=out)


def part_D(quick=False):
    e = TK.Exp1(1)
    thr = np.array([TK.entry_threshold_cusp(e, K, R, D0_BENCH / 2.402) for (K, R) in COMBOS])      # raw correlation lengths; 1 unit ~ 2.402 sigma (half-height width 1.665 sigma / 0.693)
    # the smoothed cusp  Gamma_delta(d) = exp(-(sqrt(d^2 + delta^2) - delta)): two features, entrant at offset d0, K = 2, R = 1: locking
    res = []
    for delta in ([0.02] if quick else [0.05, 0.02, 0.005]):
        class Sm(TK.Gram):
            name = 'smooth-cusp'
            n = 1
            lam0 = 0.5

            def __init__(self):
                self.s = 1.0
                self._tab = None

            def raw(self, order, x):
                return None

            def rawf(self, order, x):
                x = np.asarray(x, float)
                rr = np.sqrt(x * x + delta * delta)
                e = np.exp(-(rr - delta))
                if order == 0:
                    return e
                if order == 1:
                    return -x / rr * e
                return -(delta * delta) / rr ** 3 * e + (x / rr) ** 2 * e
        sm = Sm()
        D, d0 = 6.0, 0.35
        o = TP.inner_run(sm, 2, D, d0, 1.0, mob='v43', muc=1.0, v0=1e-2, tmax=400.0)
        dent = o['d'][1]                                              # the entrant's offset against time
        idx = np.flatnonzero(dent < 0.02)
        t_lock = float(o['t'][idx[0]]) if len(idx) else np.nan
        cm = TK.box_cusp(TK.Exp1(1), 2, 1.0, D, d0)['Cmargin']           # psi_min(r) - E_c of Proposition G'' for the exact cusp, Delta = 6, r = d_0
        res.append([delta, o['max_d'], o['final_d'], o['final_v_err'], t_lock, d0 / (1.0 * 1e-2 * cm)])
    return dict(D_thr=thr, D_smooth=np.array(res))


def part_E(quick=False):
    g = TK.Gaussian(1)
    out = np.full((len(E_CONDS), 2), np.nan)
    for a, c in enumerate(E_CONDS[:2] if quick else E_CONDS):
        out[a, 0] = TK.entry_threshold(g, 2, 1.0, D0_BENCH, lo_=3.0, hi=60.0, cond_c=c)
        out[a, 1] = TK.entry_threshold(g, 2, 2.4, D0_BENCH, lo_=3.0, hi=60.0, cond_c=c)
    return dict(E_thr=out)


def part_F(quick=False):
    cases = F_QUICK if quick else F_CASES
    if ONLY is not None:
        cases = [c for c in cases if c[0] in ONLY]
    res = np.full((len(F_CASES), len(F_OFF), 9), np.nan)
    # per case/offset: Delta, r, v_lo, v_bar, runs, n in_d, n in_hi, n in_lo, n mono + converged (packed)
    for a, (i, K, R) in enumerate(F_CASES):
        if (i, K, R) not in cases:
            continue
        k = kern(i)
        thr, r = TK.entry_threshold(k, K, R, D0_BENCH, lo_=1.0, hi=60.0, return_r=True)
        for b, off in enumerate(F_OFF):
            D = thr + off
            bx = TK.box(k, K, R, D, r)
            vlo = bx['vlo'] if bx.get('vlo', 0) > 0 else 0.5
            vbar = bx.get('vbar', R * 1.2)
            runs = 4 if quick else F_RUNS
            cnt = np.zeros(5)
            for sd in range(runs):
                mob = ['1', 'v', 'v43', 'q'][sd % 4]
                muc = [0.3, 1.0, 3.0][sd % 3]
                o = TP.box_run(k, K, R, D, r, vlo, vbar, 100 * a + sd, mob=mob, muc=muc, n=k.n, tmax=1200.0, adv=(sd % 3 == 2))
                cnt += np.array([o['in_d'], o['in_hi'], o['in_lo'], o['mono'], o['converged']], float)
            res[a, b] = [D, r, vlo, vbar, runs, cnt[0], cnt[1] + cnt[2], cnt[3], cnt[4]]
    return dict(F_res=res)


def g_targets():
    items = []
    for i in G_ROWS + G_PLANAR:
        for (K, R) in (G_KR if i in G_ROWS else [(2, 1.0)]):
            for rel in G_REL if i in G_ROWS else [0.6]:
                items.append((i, K, R, rel))
    return items


def part_G(quick=False):
    items = g_targets()
    ne = len(G_EPS)
    cols = 14
    res = np.full((len(items), cols + 4 * ne), np.nan)
    for a, (i, K, R, rel) in enumerate(items):
        if quick and not (i in (0, 15) and K == 2 and R == 1.0 and rel > 0):
            continue
        if ONLY is not None and i not in ONLY:
            continue
        k = kern(i)
        n = k.n
        thr0 = TK.entry_threshold(k, K, R, 0.5, lo_=1.0, hi=60.0)
        sep = float(np.ceil((thr0 + rel) * 10) / 10)
        seed = G_SEED0 + a
        margin = 12.0 if n == 1 else 5.0
        cs, As, c0 = TP.make_target(n, K, sep, R, seed, margin=margin)
        L = TP.limit_path(k, cs, As, c0)
        ev = L['events']
        dmax = max(e['d0'] for e in ev)
        Dact = float(min(np.linalg.norm(cs[p] - cs[q]) for p in range(K) for q in range(p)))
        Ract = float(np.abs(As).max() / np.abs(As).min())
        thr_t = TK.entry_threshold(k, K, Ract, dmax * 1.0000001, lo_=1.0, hi=60.0)
        certified = bool(Dact >= thr_t) and all(e['sign_ok'] for e in ev) and L['margin'] > 0
        fw = L['final_w']
        feats = [e['feature'] for e in ev]
        err_fin = float(max(np.abs(fw[j] - As[feats[j]]) for j in range(K)) + max(np.linalg.norm(L['final_c'][j] - cs[feats[j]]) for j in range(K)))
        distinct = len(set(feats)) == K
        res[a, :cols] = [i, K, R, rel, sep, Dact, Ract, dmax, thr_t, float(certified), float(L['margin']), err_fin, float(distinct), float(all(e['sign_ok'] for e in ev))]
        for b, eps in enumerate(G_EPS[:2] if quick else G_EPS):
            if n > 1 and eps < 1e-6 - 1e-12:
                continue
            sol, U = TP.simulate_network(k, cs, As, c0, eps, L)
            c = TP.compare(cs, As, c0, eps, L, sol, U)
            res[a, cols + 4 * b: cols + 4 * b + 4] = [float(c['formed_ok']), float(c['exact']), c['err_w'], max(abs(x) for x in c['shifts'])]
    return dict(G_res=res)


def part_H(quick=False):
    """The certificate of Theorem T1 (exp_theory_t1.cert, run for the Gram by theory3e_t1gen) on one target per Gram: the largest decade eps_0 from which it holds at every smaller one, the measured
    eps_max from which the integrated network forms exactly the path's units at every smaller one, and the number of violations (cert true, network not following)."""
    rows = [0, 15] if quick else H_ROWS
    if ONLY is not None:
        rows = [i for i in rows if i in ONLY]
    exs = H_EX[:4] if quick else H_EX
    res = np.full((len(H_ROWS), 8 + 2 * len(H_EX)), np.nan)
    for a, i in enumerate(H_ROWS):
        if i not in rows:
            continue
        k = kern(i)
        M = TG.build(k)
        M.PARS = M.PARS[:H_NPARS]
        d_t1 = M.delta_t1(2, 1.0, 0.5, 0.02, step=0.1, lo=3.0, hi=40.0)
        sep = float(np.ceil(d_t1 + 1.0)) if d_t1 == d_t1 else 12.0
        cs, As, c0 = M.make_target(2, sep, 1.0, H_SEED + a)
        prep = M.prepare(cs, As, c0)
        okc, okm = [], []
        for ex in exs:
            eps = 10.0 ** ex
            o, _ = M.best_cert(prep, eps)
            sol, U = M.simulate(prep, eps, nlocal=300)
            m, e = M.follows(prep, eps, sol, U)
            okc.append(bool(o['ok']))
            okm.append(bool(m))
        def first_all(ok):
            for j in range(len(exs)):
                if all(ok[j:]):
                    return exs[j]
            return np.nan
        viol = sum(1 for c_, m_ in zip(okc, okm) if c_ and not m_)
        res[a, :8] = [i, d_t1, sep, prep['Delta'], prep['mu_min'], first_all(okc), first_all(okm), viol]
        for j in range(len(exs)):
            res[a, 8 + 2 * j], res[a, 9 + 2 * j] = float(okc[j]), float(okm[j])
    return dict(H_res=res)


# ------------------------------------------------------------------------------------------------------------------------------------
def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    out += ['| ' + ' | '.join(r) + ' |' for r in rows]
    return '\n'.join(out)


def f2(x, f='%.2f'):
    return '-' if x != x else f % x


STATUS = [
    'Scope: every statement for a non-Gaussian Gram is for the CONTINUUM inner product (L2(R^n)) with the data spacing small against the kernel\'s correlation length. On sampled data the Gram of the Laplace '
    'basis has a cusp in c at every data point, so Proposition G does not apply to a sampled Laplace network (it behaves like a cusped Gram near data points). No non-Gaussian network exists in the '
    'repository; the network evidence for these Grams is the continuum simulator `theory3e_path.simulate_network`.',
    'Wording: "Laplace BASIS function exp(-|x|)", not "Laplace kernel". If "kernel" were read as the Gram exp(-d), the case would be the cusp Proposition G\'\', not Proposition G. The rational-quadratic and '
    'Wendland rows are Gram classes (the rational-quadratic Gram has no natural basis function; Wendland is the basis function, its Gram is the support-2 autocorrelation).',
    'Tail envelopes: |Gamma|, |Gamma\'|, |Gamma\'\'| beyond L = 40 sigma-equivalents are ASSUMED non-increasing (checked numerically; hand-proved for the half-integer Matern Grams only). '
    'The tabulated rows are computer-assisted sufficient thresholds conditional on this assumption.',
    'T1 beyond the Gaussian: only the computed instances of Part H (one target each for four n = 1 Grams), not Theorem T1 for general Grams; no slowly decaying (algebraic) Gram was run, and for '
    'the rational-quadratic lattice sums 2 beta > n is needed.',
    'The network-following evidence of Part G does NOT test certification: uncertified targets follow the recursion as well as certified ones (counts in Part G); the threshold Delta_entry is a '
    'sufficient, conservative bound, not a transition separation. It is not used as validation.',
]


def report(raw):
    L = []
    L.append('# Theory round 3, direction E: the entry layer and the limit path for a general Gram\n')
    L.append('Generated by `experiments/theory3e_run.py` (kernels: `theory3e_kernels.py`, `theory3e_wendland.py`; simulator: `theory3e_path.py`); the statements and proofs are in the round\'s NOTE.md '
             '(Proposition G, Corollary G1, Proposition G\'\', Theorem G). Lengths are in sigma-equivalents, defined once here: the curvature-matched unit in which -Gamma\'\'(0) = 1/2, so that the Gaussian of width sigma has '
             'Gamma(d) = exp(-d^2/4) (for a non-Gaussian basis function it is NOT the width of anything); the starting offset is d_0 = 0.4167, half the dictionary spacing 0.833 in these units '
             '(the benchmark\'s h/(2 sigma)); every threshold is conditional on that spacing, unless stated. Ball arithmetic: python-flint; every kernel number is a directed '
             'float bound turned into an exact ball (tails: upper; own minima: lower), so `ok` means the strict inequalities are certainly true (beyond L = 40 sigma the three envelopes '
             '|Gamma|, |Gamma\'|, |Gamma\'\'| are assumed non-increasing; hand-proved for the half-integer Matern Grams only).\n')
    L.append('**Status and scope** (referee round T3E, `review/REVIEW_roundT3E.md`).\n')
    for s_ in STATUS:
        L.append('* ' + s_)
    L.append('')

    if 'A_thr' in raw:
        thr, rs, cl = raw['A_thr'], raw['A_r'], raw['A_clean']
        L.append('## A. Entry thresholds Delta_entry(K, R, d_0 = 0.4167) for the Gram classes (n = dimension of the data space)\n')
        gr, pp = raw['A_gauss_rowsum'], raw['A_paper']
        L.append('Regression: the Gaussian Gram with the Gaussian note\'s row-sum cross-Hessian bound (n = 1), computed by the general code, against the paper\'s table `tab_entry`: ' +
                 ', '.join('(%d, %.1f): %.2f vs %.2f' % (K, R, g_, p_) for (K, R), g_, p_ in zip(COMBOS, gr, pp)) + (' -- agree on all five.' if np.allclose(gr, pp, atol=1e-9) else ' -- DO NOT AGREE.') + '\n')
        hdr = ['Gram (basis function)'] + ['K=%d, R=%g' % c for c in COMBOS]
        rows = []
        for i, (nm, _) in enumerate(ROWS):
            if thr[i].min() == thr[i].min():
                rows.append([nm] + [f2(x) for x in thr[i]])
        L.append('Sharp thresholds (Proposition G, blockwise spectral-norm cross-Hessian, any n), first Delta on a grid of 0.02 for which some radius r > d_0 satisfies (L), (C), (N):\n')
        L.append(table(hdr, rows) + '\n')
        rows = []
        for i, (nm, _) in enumerate(ROWS):
            if cl[i].min() == cl[i].min():
                rows.append([nm] + [f2(x) for x in cl[i]])
        L.append('Clean thresholds (Corollary G1: (K-1) R^2 max(m0, m1, m2)(Delta - 2r) <= c_*(r)), same grid; "-" = no r with c_* > 0 or no separation up to 80:\n')
        L.append(table(hdr, rows) + '\n')
        L.append('(`Gram` = the Gram <phi_c phi_c\'>/Z of the basis function: the Laplace basis exp(-|x|) in R^n has the Gram Matern-(1 + n/2), the Matern-nu basis the Gram Matern-(2 nu + n/2); '
                 'verified against numerical autocorrelations in the tests. The Matern-3/2 and Matern-5/2 *as Grams* are the rows `Laplace basis (n=1)` and `Laplace basis (n=3)`.)\n')
    if 'B_thr' in raw:
        L.append('## B. The starting offset\n')
        t = raw['B_thr']
        rows = [[ROWS[i][0]] + [f2(x) for x in t[a]] for a, i in enumerate(B_ROWS)]
        L.append('Delta_entry for K = 2, R = 1 against d_0:\n')
        L.append(table(['Gram'] + ['d_0=%g' % d for d in B_D0], rows) + '\n')
        lat = raw['B_lat']
        rows = [[ROWS[i][0], f2(lat[i, 0], '%.3f'), f2(lat[i, 1])] for i in range(len(ROWS)) if lat[i, 0] == lat[i, 0]]
        L.append('Worst offset sqrt(n) h / 2 of a cubic lattice of spacing h = 0.8333 sigma (the benchmark\'s 0.1/0.12) in R^n, K = 2, R = 1:\n')
        L.append(table(['Gram', 'd_0 = sqrt(n) 0.4167', 'Delta_entry'], rows) + '\n')
    if 'C_mod' in raw:
        L.append('## C. Decay moduli\n')
        C = raw['C_mod']
        rows = []
        for i, (nm, _) in enumerate(ROWS):
            if C[i, 0] == C[i, 0]:
                q, a, c, cs_, dl = C[i]
                form = ('C l^-%g' % a if q == 0.0 else ('compact support' if a == 0.0 else 'C exp(-%.3f l^%g)' % (a, q)))
                rows.append([nm, form, f2(c, '%.3g'), f2(cs_, '%.4f'), f2(dl)])
        L.append('m(l) = max(m0, m1, m2)(l) <= C exp(-a l^q) (or C l^{-2 beta}) for l in [3, 40] (from the table of tail envelopes); c_*(0.45) of Corollary G1 for R = 1; and the explicit separation 2 r + l with '
                 'K = 2, R = 1, r = 0.45 and (K-1) R^2 C exp(-a l^q) = c_*:\n')
        L.append(table(['Gram', 'modulus', 'C', 'c_*(r=0.45)', 'Delta from modulus'], rows) + '\n')
    if 'D_thr' in raw:
        L.append('## D. The cusped Gram exp(-d) (Proposition G\'\')\n')
        t = raw['D_thr']
        L.append('Conditions (L), (C\'), kappa < 1 (no Hessian condition; the centres lock onto the features). Delta_entry in correlation lengths (1 unit ~ 2.402 sigma by half-height width) and in sigma-equivalents: ')
        L.append(table(['K, R'] + ['%d, %g' % c for c in COMBOS], [['Delta_entry (corr. lengths)'] + [f2(x) for x in t], ['in half-height sigma units (x 2.402)'] + [f2(2.402 * x) for x in t]]) + '\n')
        d = raw['D_smooth']
        L.append('The cusp smoothed, Gamma_delta(d) = exp(-(sqrt(d^2 + delta^2) - delta)), inner problem (K = 2, Delta = 6, d_0 = 0.35, mobility v^{4/3}) integrated: the entrant and the block end at the target with the centre offsets '
                 'of order delta (locking), never leaving radius d_0 + O(delta):\n')
        L.append(table(['delta', 'max offset over the run', 'final offset', 'final coefficient error', 'time to offset < 0.02', 'bound r / (mu_c v_0 (psi_min - E_c))'],
                       [[f2(x[0], '%g'), f2(x[1], '%.3f'), f2(x[2], '%.1e'), f2(x[3], '%.1e'), f2(x[4], '%.1f'), f2(x[5], '%.0f')] for x in d]) + '\n')
    if 'E_thr' in raw:
        L.append('## E. Gaussian with a common anisotropic covariance S\n')
        t = raw['E_thr']
        L.append('(C) becomes psi(r) > cond(S) E_c for the centre mobility mu S^{-1} (the Euclidean-gradient flow in the whitened coordinates); the separation is the Mahalanobis one. Delta_entry (whitened units), K = 2:\n')
        L.append(table(['cond(S)', 'R = 1', 'R = 2.4'], [[f2(c, '%g'), f2(t[a, 0]), f2(t[a, 1])] for a, c in enumerate(E_CONDS) if t[a, 0] == t[a, 0]]) + '\n')
    if 'F_res' in raw:
        L.append('## F. Proposition G on the inner problem: random and adversarial initial points of the box\n')
        F = raw['F_res']
        rows = []
        for a, (i, K, R) in enumerate(F_CASES):
            for b, off in enumerate(F_OFF):
                x = F[a, b]
                if x[0] != x[0]:
                    continue
                rows.append([ROWS[i][0], '%d, %g' % (K, R), '%+.1f' % off, f2(x[0]), f2(x[1], '%.3f'), f2(x[2], '%.3f'), f2(x[3], '%.3f'), '%d/%d' % (x[5], x[4]), '%d/%d' % (x[6], 2 * x[4]), '%d/%d' % (x[7], x[4]), '%d/%d' % (x[8], x[4])])
        L.append('Start: every unit at a random point of B (|d_i| <= r, |w_i| in [v_lo, v_bar]; the last unit\'s coefficient in (0, v_lo]; every third run on the boundary, offsets +-0.999 r, coefficients at the faces), random signs, mobility law in '
                 '{1, v, v^{4/3}, 1 + v^2}, centre mobility in {0.3, 1, 3}; the offset column is relative to Delta_entry (0 = the certified threshold). Counts of runs: offsets stay <= r; coefficients stay in [min(start, v_lo), v_bar] (two tests per run); '
                 'entrant increasing while below v_lo; convergence to the target (error < 1e-6). For negative offsets the box does not exist (the conditions fail), and the v_lo, v_bar columns are then the failed evaluation\'s numbers or the defaults 0.5 and 1.2 R; '
                 'the runs are the same experiment and show how conservative Delta_entry is.\n')
        L.append(table(['Gram', 'K, R', 'Delta - Delta_entry', 'Delta', 'r', 'v_lo', 'v_bar', 'offsets in box', 'coeff. in box', 'entrant monotone', 'converged'], rows) + '\n')
    if 'G_res' in raw:
        L.append('## G. The limit path against the network flow\n')
        G = raw['G_res']
        items = g_targets()
        cols = 14
        rows = []
        for a, (i, K, R, rel) in enumerate(items):
            x = G[a]
            if x[0] != x[0]:
                continue
            e = []
            for b, eps in enumerate(G_EPS):
                y = x[cols + 4 * b: cols + 4 * b + 4]
                e.append('-' if y[0] != y[0] else ('%s/%s %.0e' % ('Y' if y[0] > 0.5 else 'n', 'Y' if y[1] > 0.5 else 'n', y[3])))
            rows.append([ROWS[i][0], '%d, %.2f' % (K, x[6]), f2(x[5]), f2(x[8]), 'yes' if x[9] > 0.5 else 'no', f2(x[10], '%.4f'), f2(x[11], '%.0e')] + e)
        L.append('One target per row (`theory3e_path.make_target`, fixed seeds): K features on a line (random direction in R^n) at separation Delta, the dictionary a cubic lattice of spacing 0.833 sigma. '
                 '"threshold" is Delta_entry for the target\'s own K, R and its largest entrant offset; "certified" = Delta >= threshold, every entrant has its feature\'s sign (a) and all inactive duals stay below the wall (b); '
                 '"margin" = 1 - max |u_j|; "final err" = |block - target| of the recursion. Network columns: formed exactly the recursion\'s units / exact fit (error < 1e-3 in the coefficients, < 0.05 in the centres), '
                 'and the largest relative entry-time shift |t_net/(lambda tau_k) - 1|.\n')
        L.append(table(['Gram', 'K, R', 'Delta', 'threshold', 'certified', 'margin', 'final err'] + ['eps=%g' % e for e in G_EPS], rows) + '\n')
        cert = [G[a] for a in range(len(items)) if G[a, 0] == G[a, 0] and G[a, 9] > 0.5]
        unc = [G[a] for a in range(len(items)) if G[a, 0] == G[a, 0] and G[a, 9] < 0.5]
        L.append('Recursion on certified targets: %d targets, exactly K distinct entrants with the right signs and a final block equal to the target (error < 1e-9) on %d. On uncertified targets: %d, of which %d have the same properties. ' % (
            len(cert), sum(1 for x in cert if x[11] < 1e-9 and x[12] > 0.5), len(unc), sum(1 for x in unc if x[11] < 1e-9 and x[12] > 0.5)))
        for b, eps in enumerate(G_EPS):
            ok = [x for x in cert if x[cols + 4 * b] == x[cols + 4 * b]]
            L.append('At eps = %g the network formed exactly the recursion\'s units on %d of %d certified targets (%d of %d uncertified).' % (
                eps, sum(1 for x in ok if x[cols + 4 * b] > 0.5), len(ok), sum(1 for x in unc if x[cols + 4 * b] == x[cols + 4 * b] and x[cols + 4 * b] > 0.5), sum(1 for x in unc if x[cols + 4 * b] == x[cols + 4 * b])))
        L.append('')
    if 'H_res' in raw:
        L.append('## H. Theorem T1\'s certificate recomputed for the Gram (n = 1; `theory3e_t1gen`: the repository\'s `exp_theory_t1.cert` run unchanged with the Gaussian pieces replaced)\n')
        Hh = raw['H_res']
        rows = []
        for a, i in enumerate(H_ROWS):
            x = Hh[a]
            if x[0] != x[0]:
                continue
            exs = H_EX
            cert_s = ''.join('Y' if x[8 + 2 * j] > 0.5 else ('.' if x[8 + 2 * j] == x[8 + 2 * j] else ' ') for j in range(len(exs)))
            net_s = ''.join('Y' if x[9 + 2 * j] > 0.5 else ('.' if x[9 + 2 * j] == x[9 + 2 * j] else ' ') for j in range(len(exs)))
            rows.append([ROWS[i][0], f2(x[1], '%.1f'), f2(x[2], '%.0f'), f2(x[4], '%.4f'), '1e%+d' % x[5] if x[5] == x[5] else 'none', '1e%+d' % x[6] if x[6] == x[6] else 'none', '%d' % x[7], cert_s, net_s])
        L.append('One K = 2, R = 1 target per Gram at separation ceil(Delta_T1 + 1) (Delta_T1: the separation from which the box and both Gershgorin margins hold at d_0 = 0.5, dormant budget 0.02); eps on the decade grid '
                 + ', '.join('1e%d' % e for e in H_EX) + ' (strings: Y = holds, . = fails); `cert` with the first %d parameter sets of the repository\'s list. eps_0 = the largest decade from which the certificate holds at all smaller '
                 'ones; eps_max = the same for the integrated network forming exactly the path\'s units; violations = eps with the certificate true and the network not following.\n' % H_NPARS)
        L.append(table(['Gram', 'Delta_T1', 'sep', 'margin', 'eps_0 (cert)', 'eps_max (net)', 'violations', 'cert', 'network'], rows) + '\n')
    return '\n'.join(L) + '\n'


def load_raw():
    with np.load(RAW) as z:
        return {k: z[k] for k in z.files}


def save_raw(raw):
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(RAW, **raw)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--parts', default='A,B,C,D,E,F,G,H')
    ap.add_argument('--only', default=None, help='comma-separated Gram row indices: recompute only these rows of the row-indexed parts and merge into the stored raw')
    a = ap.parse_args(argv)
    if a.run:
        global ONLY
        warnings.filterwarnings('ignore')
        raw = load_raw() if os.path.exists(RAW) else {}
        ONLY = None if a.only is None else {int(x) for x in a.only.split(',')}
        for p in a.parts.split(','):
            t0 = time.time()
            new = {'A': part_A, 'B': part_B, 'C': part_C, 'D': part_D, 'E': part_E, 'F': part_F, 'G': part_G, 'H': part_H}[p](a.quick)
            if ONLY is not None:                                    # merge: keep the stored values where this run left NaN
                new = {k: (np.where(np.isnan(v), raw[k], v) if k in raw and raw[k].shape == v.shape else v) for k, v in new.items()}
            raw.update(new)
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
