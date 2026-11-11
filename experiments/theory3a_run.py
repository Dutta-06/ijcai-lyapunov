"""
Theory round 3, direction A: close "S(eps) holds for all eps <= eps_0" of Theorem T1 (the tail eps -> 0).

Parts (raw data in results/synthetic/theory3a_raw.npz, report results/synthetic/theory3a.md, `--from-raw` regenerates the report byte for byte):
  U  the eps-free envelope of the dormant coefficient (Lemma 12.4) evaluated in BALL arithmetic on a grid of (eps, x): w(x) <= C/(s (s^2 - sigma_0^2)) for every grid point with s >= 25;
  S  the direction in which each layer ingredient of S(eps) moves with eps (point balls from 1e-2 to 1e-24 against the ball enclosing all eps <= 1e-24): which are monotone, which are not;
  T  the TAIL: for each of the 8 targets of Part Q of theory_t1b.md the largest eps_hi on a quarter-decade grid for which `theory3a_tail.verify_tail` proves S(eps) for EVERY real eps in (0, eps_hi]
     with one set of bounds D_j (h = 1/kappa a ball [0, h_hi]; the dormant side by the eps-free envelope, in powers of L = 1/lambda);
  C  the chain [eps_tail, eps_top]: the ball cover of `theory_t1b_wide.cover` re-run with `theory3a_chain.verify_fixed` (the sound displacement accounting of Lemma 5.2).
Usage: python experiments/theory3a_run.py --run [--quick] [--parts U,S,T,C] | --from-raw      (environment TH3A_OUT for the output directory)
"""
import argparse
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from flint import arb

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_theory_t1 as T          # noqa: E402
import theory_t1b_wide as W        # noqa: E402
import theory3a_tail as TL         # noqa: E402
import theory3a_chain as CH        # noqa: E402

OUT = os.environ.get('TH3A_OUT', os.path.join(HERE, '..', 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory3a_raw.npz'), os.path.join(OUT, 'theory3a.md')
QRAW = os.path.join(os.environ.get('TH3A_Q', os.path.join(HERE, '..', 'results', 'synthetic')), 'theory_t1b_raw.npz')
Q_TARGETS = (0, 1, 2, 3, 4, 5, 6, 8)          # rows of Part Q of theory_t1b_raw.npz whose stored chain is complete
S_TARGETS = (0, 3)                            # one K = 1 and one K = 2 target for the direction table
S_EX = np.array([-2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0, -9.0, -10.0, -12.0, -14.0, -16.0, -20.0, -24.0])
U_EX = (-2.0, -4.0, -6.0, -8.0, -12.0, -16.0, -20.0, -24.0)
S_NAMES = ('x1', 'w1', 'xs', 'T1', 'Atr', 'T2', 'T3', 'Env', 'Tlay', 'kappa', 'pre', 'pre_drift', 'mV')


def part_U(quick=False):
    """w(x) <= env(s(x)) in ball arithmetic: for each eps and each x of a geometric grid with lo(s) >= 25, the upper end of w against the lower end of C/(s(s^2 - sigma_0^2))."""
    rows = []
    for e in (U_EX[:2] if quick else U_EX):
        U = TL.UH(TL.h_of_eps(10.0 ** e))
        n_pts, worst = 0, 0.0
        for x in np.geomspace(1e-14, 1.6, 120 if quick else 600):
            xa = arb(float(x))
            s = U.s(xa)
            if not float(TL.lo(s)) >= 25.0:
                continue
            env_lo = TL.C_ENV / (TL.up(s) * (TL.up(s) * TL.up(s) - TL.SIG0 * TL.SIG0))             # env decreases in s: its lower end is at the upper end of s
            ratio = float(TL.up(U.w(xa) / env_lo))
            n_pts += 1
            worst = max(worst, ratio)
        rows.append([e, n_pts, worst])
    return dict(U_rows=np.array(rows, float))


def layer_vec(prep, par, hb):
    ly = TL.layer(prep, par, hb)
    if ly['why']:
        return None
    ex = lambda b: [float(b.lower()), float(b.upper())]
    return ex(ly['x1']) + ex(ly['w1']) + ex(ly['xs']) + ex(ly['T1']) + ex(ly['Atr']) + ex(ly['T2']) + ex(ly['T3']) + ex(ly['Env']) + ex(ly['Tlay']) + ex(ly['kappa']) + ex(ly['pre']) + ex(ly['pre_drift'][0]) + ex(ly['U'].M(arb(float(TL.lo(ly['xs'])))))


def part_S(quick=False):
    """For two targets: the layer ingredients (enclosures [lo, up]) at point eps = 10^e and for the ball of all eps <= 1e-24."""
    craw = np.load(QRAW, allow_pickle=False)['Q_rows']
    out = []
    for n in (S_TARGETS[:1] if quick else S_TARGETS):
        K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
        prep = T.prepare(*T.make_target(K, sep, R, seed))
        o, pi = T.best_cert(prep, 10.0 ** craw[n, 5])
        par = T.PARS[pi]
        lim = layer_vec(prep, par, arb(0).union(TL.h_of_eps(1e-24)))
        out.append([n, -99.0] + lim)
        for e in (S_EX[[0, 3, 6, 9, 13]] if quick else S_EX):
            hb = TL.h_of_eps(10.0 ** e)
            try:
                v = layer_vec(prep, par, TL.A(hb))
            except ValueError:
                v = None
            if v is None:
                continue
            out.append([n, e] + v)
        print('S target', n, flush=True)
    return dict(S_rows=np.array(out, float))


def tail_one(n, quick=False):
    """Largest quarter-decade eps_hi (from e0 + 0.5 down) for which verify_tail holds; then verification at 1e-12 and 1e-24 as well (regression: the verdict must not depend on the grid point)."""
    craw = np.load(QRAW, allow_pickle=False)['Q_rows']
    K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
    prep = T.prepare(*T.make_target(K, sep, R, seed))
    ent = {e['kernel'] for e in prep['L']['events']}
    J = len(prep['c0'])
    e0, top = craw[n, 4], craw[n, 5]
    start = np.floor((e0 + 0.5) * 4) / 4
    res = None
    for e in np.arange(start, -24.01, -0.25):
        o, pi = T.best_cert(prep, 10.0 ** e)
        if pi < 0:
            continue
        ok, info, Dv = TL.solve_D(prep, T.PARS[pi], TL.h_of_eps(10.0 ** e), o['Etot'], ladder=(1.1, 1.3, 1.6, 2.0), iters=3)
        if ok:
            res = (e, pi, info, Dv)
            break
    if res is None:
        return None
    e, pi, info, Dv = res
    extra = []
    for e2 in ((-12.0,) if quick else (-12.0, -24.0)):
        o2, pi2 = T.best_cert(prep, 10.0 ** e2)
        ok2, info2, Dv2 = TL.solve_D(prep, T.PARS[pi2], TL.h_of_eps(10.0 ** e2), o2['Etot'], ladder=(1.1, 1.3, 1.6, 2.0), iters=3)
        extra.append([e2, float(ok2), info2.get('T1', np.nan), info2.get('pre', np.nan) - info2.get('prelo', np.nan), max(info2['Etot'].values()) if 'Etot' in info2 else np.nan])
    Dn = np.array([Dv[j] for j in range(J) if j not in ent])
    En = np.array([info['Etot'][j] for j in range(J) if j not in ent])
    row = [n, K, sep, R, seed, e0, top, e, pi, float(TL.h_of_eps(10.0 ** e).upper()), info['Lhi'], En.max(), Dn.max(), (En / Dn).max(), info['drift_max'], info['mass'], info['dcmax'], info['omega'],
           info['kappa'], info['T1lo'], info['T1'], info['prelo'], info['pre'], info['Tlay'], info['cum_max']]
    return dict(row=row, extra=extra, D=np.array(Dv, float), E=np.array([info['Etot'].get(j, 0.0) for j in range(J)], float))


T_COLS = ('n', 'K', 'sep', 'R', 'seed', 'e0', 'top', 'etail', 'par', 'h_hi', 'L_hi', 'Emax', 'Dmax', 'E/D', 'drift', 'mass', 'dcmax', 'omega', 'kappa', 'T1lo', 'T1up', 'prelo', 'preup', 'Tlay', 'cum')


def part_T(quick=False):
    warnings.filterwarnings('ignore')
    rows, extra, Ds, Es = [], [], {}, {}
    for n in (Q_TARGETS[:1] if quick else Q_TARGETS):
        t0 = time.time()
        r = tail_one(n, quick)
        if r is None:
            print('T target', n, 'no tail found', flush=True)
            continue
        rows.append(r['row'])
        extra.append([n] + r['extra'][0] + (r['extra'][1] if len(r['extra']) > 1 else [np.nan] * 5))
        Ds['T_D%d' % n] = r['D']
        Es['T_E%d' % n] = r['E']
        print('T target', n, 'etail', r['row'][7], 'top', r['row'][6], '%.0fs' % (time.time() - t0), flush=True)
    out = dict(T_rows=np.array(rows, float), T_extra=np.array(extra, float))
    out.update(Ds)
    out.update(Es)
    return out


def part_C(raw, quick=False):
    """The chain from eps_top down to eps_tail with the sound verifier (monkeypatched into `theory_t1b_wide.cover`, which calls RG.verify)."""
    warnings.filterwarnings('ignore')
    W.RG.verify = CH.verify_fixed
    rows, segs_all = [], []
    for r in raw['T_rows']:
        n = int(r[0])
        K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
        top, etail = r[6], r[7]
        if etail >= top:
            rows.append([n, top, etail, 1.0, 0.0, 0.0, 0.0, etail, 0.0])
            segs_all.append(np.zeros((0, 3)))
            continue
        prep = T.prepare(*T.make_target(K, sep, R, seed))
        t0 = time.time()
        segs, complete, cur, calls = W.cover(prep, top, etail, budget=None)
        wmin = min([b - a for a, b in segs]) if segs else 0.0
        rows.append([n, top, etail, float(complete), float(len(segs)), float(calls), wmin, cur, time.time() - t0])
        segs_all.append(np.array([[n, a, b] for a, b in segs], float).reshape(-1, 3))
        print('C target', n, 'top', top, 'etail', etail, 'complete', complete, 'balls', len(segs), 'calls', calls, '%.0fs' % (time.time() - t0), flush=True)
    return dict(C_rows=np.array(rows, float), C_segs=np.concatenate(segs_all) if segs_all else np.zeros((0, 3)))


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(r) + ' |')
    return '\n'.join(out)


def sgn_dir(vals, rtol=1e-9):
    """classify a sequence (ordered by decreasing eps): 'non-increasing' / 'non-decreasing' / 'constant' / 'non-monotone' (up to a relative 1e-13)"""
    v = np.asarray(vals, float)
    d = np.diff(v)
    tol = rtol * np.maximum(np.abs(v[:-1]), 1e-300)
    if np.all(np.abs(d) <= tol):
        return 'constant'
    if np.all(d <= tol):
        return 'non-increasing'
    if np.all(d >= -tol):
        return 'non-decreasing'
    return 'non-monotone'


BANNER = 'SUPERSEDED by results/synthetic/theory3n.md: the displacement accounting of the certificate was corrected (see review/REVIEW_roundT3AC.md, REVIEW_roundT3N.md); numbers below use the legacy accounting and are kept as history'


def report(raw):
    L = [BANNER + '\n', '# Theorem T1, the small-eps tail: ball-arithmetic proof of S(eps) for every real eps in (0, eps_hi]  (generated by experiments/theory3a_run.py)\n']
    L.append('Round theory3 A.  Statements and their status are in the report text below; the proofs by hand are in the working note THEORY3A.md of the author (outside the repository).  Verifier: '
             '`theory3a_tail.verify_tail` (python-flint, 100 bits); chain: `theory_t1b_wide.cover` with `theory3a_chain.verify_fixed`.  The limit object\'s data (critical points, entry times, duals) are floating-point '
             'inputs taken as exact, exactly as in `theory_t1b_rigorous` (not enclosed here).\n')
    if 'U_rows' in raw:
        L.append('## U. The eps-free envelope of the dormant coefficient in ball arithmetic\n')
        L.append('For each eps (the unit with h = 1/kappa a point ball) and a geometric grid of x, the ball enclosure of w(x) is compared with the lower end of C/(s(s^2 - sigma_0^2)) at the ball enclosure of the distance s(x) '
                 '= lambda - eta(x), for all grid points with s >= 25 (C = 100/3 = 4/(G gamma^2 m), sigma_0 = 20).  "worst" = the largest upper end of w/env (must be <= 1).\n')
        L.append(table(['eps', 'grid points with s>=25', 'worst w/env (upper end)'], [['1e%+.0f' % r[0], '%d' % r[1], '%.9f' % r[2]] for r in raw['U_rows']]))
        L.append('\nAll ratios are <= 1: ' + ('yes' if bool(np.all(raw['U_rows'][:, 2] <= 1.0)) else 'NO') + '.\n')
    if 'S_rows' in raw:
        L.append('## S. How the layer ingredients of S(eps) move with eps\n')
        L.append('The ingredients of Lemma 3.3 / 4.1 / 5.3 (enclosures at the point eps = 10^e, for the same parameters pi) against the ball of ALL eps <= 1e-24 ("limit").  For an ingredient that enters S as an upper bound '
                 '(T1, Atr, T2, T3, Env, Tlay, pre, pre_drift) the table uses the upper end of the enclosure, for one that enters as a divisor (kappa, mV) the lower end, for the closed-form positions (x1, w1, xs) the midpoint.  '
                 '"adverse as eps decreases" = the ingredient gets WORSE for S (an upper-type ingredient larger, a lower-type smaller) as eps decreases, i.e. S(eps_1) does not by itself bound the ingredient at smaller eps; '
                 '"excess" = the relative amount by which the limit is worse than the value at eps (positive = adverse).  T1, Atr, T2, T3, Env and pre_drift are upper Riemann sums / guards on a fixed number of cells '
                 '(1500, 3000): their relative discretisation error (about 1e-7 .. 1e-6) limits what the last digits of their direction mean; pre, x1, w1, xs, kappa, mV are closed-form evaluations.\n')
        S = raw['S_rows']
        kind = {'x1': 'mid', 'w1': 'mid', 'xs': 'mid', 'kappa': 'lower', 'mV': 'lower'}
        tol = {'T1': 1e-5, 'Atr': 1e-5, 'Tlay': 1e-5, 'pre_drift': 1e-3}
        for n in sorted(set(S[:, 0])):
            Sn = S[S[:, 0] == n]
            lim = Sn[Sn[:, 1] == -99.0][0]
            pts = Sn[Sn[:, 1] != -99.0]
            K, sep, R, seed = T.C_TARGETS[W.Q_IDX[int(n)]]
            L.append('\n### target (K, sep, R, seed) = (%d, %g, %g, %d)\n' % (K, sep, R, seed))
            hdr = ['ingredient', 'type', 'limit (eps -> 0)', 'value at 1e-4', 'value at 1e-8', 'excess at 1e-4', 'excess at 1e-8', 'direction as eps decreases', 'adverse as eps decreases']
            rows = []
            for i, nm in enumerate(S_NAMES):
                k_ = kind.get(nm, 'upper')
                col = 3 + 2 * i if k_ == 'upper' else 2 + 2 * i
                val = (0.5 * (pts[:, 2 + 2 * i] + pts[:, 3 + 2 * i])) if k_ == 'mid' else pts[:, col]
                lv = 0.5 * (lim[2 + 2 * i] + lim[3 + 2 * i]) if k_ == 'mid' else lim[col]
                v4 = val[np.argmin(np.abs(pts[:, 1] + 4.0))]
                v8 = val[np.argmin(np.abs(pts[:, 1] + 8.0))]
                sg = -1.0 if k_ == 'lower' else 1.0
                ex = lambda v: sg * (lv - v) / max(abs(v), 1e-300)
                d = sgn_dir(val, tol.get(nm, 1e-9))
                if k_ == 'mid':
                    adv = 'n/a (position)'
                elif d == 'constant':
                    adv = 'no change'
                else:
                    adv = 'yes' if (d == 'non-decreasing') == (k_ == 'upper') and d != 'non-monotone' else 'no' if d != 'non-monotone' else 'mixed (mesh-limited)'
                rows.append([nm, k_, '%.10g' % lv, '%.10g' % v4, '%.10g' % v8, '%.2e' % ex(v4), '%.2e' % ex(v8), d, adv])
            L.append(table(hdr, rows))
        L.append('\nReading: "adverse" ingredients (pre, T3, Env, kappa, mV) are monotone in the wrong direction by a relative amount proportional to eps (about 1.6 eps for pre), so the value at eps_1 does not bound them at smaller eps; '
                 'the supremum over eps <= eps_1 is the limit value (eps = 0), which is what the ball [0, h_hi] of Part T encloses.  T1, Atr, Tlay and pre_drift move in the favourable direction (larger at larger eps).\n')
    if 'T_rows' in raw:
        L.append('## T. The tail: S(eps) for every real eps in (0, eps_hi]\n')
        L.append('`verify_tail` with the unit in h = 1/kappa a ball [0, h_hi], one set of bounds D_j, strict test E_j < D_j/1.0001.  eps_hi is the largest value on the quarter-decade grid (from e0 + 0.5 down) '
                 'for which the verifier accepts; e0 = the floating-point certified threshold of theory_t1.md, top = the top of the stored interval cover (Part Q of theory_t1b.md).  E = the proved upper bound '
                 'of the dual deviation, D = the candidate bound.\n')
        hdr = ['K', 'sep', 'R', 'seed', 'log10 e0', 'log10 top', 'log10 eps_hi', 'max E_j', 'max D_j', 'max E_j/D_j', 'drift', 'mass', 'dc max', 'omega', 'T1 (upper sum)', 'pre enclosure width', 'L_hi=1/lambda']
        rows = []
        for r in raw['T_rows']:
            rows.append(['%d' % r[1], '%g' % r[2] if r[1] > 1 else '-', '%g' % r[3], '%d' % r[4], '%+.2f' % r[5], '%+.2f' % r[6], '%+.2f' % r[7], '%.3f' % r[11], '%.3f' % r[12], '%.4f' % r[13],
                         '%.3g' % r[14], '%.3g' % r[15], '%.3g' % r[16], '%.3g' % r[17], '%.9f' % r[20], '%.2e' % (r[22] - r[21]), '%.3g' % r[10]])
        L.append(table(hdr, rows))
        n_ge = int(sum(r[7] >= r[6] for r in raw['T_rows']))
        L.append('\n%d of %d targets: eps_hi >= eps_top (the tail alone covers the whole range of the stored cover).\n' % (n_ge, len(raw['T_rows'])))
        L.append('Regression and the explicit remainder: the verifier is re-run at eps_hi = 1e-12 and 1e-24 (its own D_j).  The upper bound of T_1 (an upper Riemann sum) and the enclosure width of pre (closed form; '
                 'the ingredient with the largest gap to its limit) shrink with eps_hi: the remainder of the layer constants is O(eps_hi).  max E converges to the eps -> 0 value of the deviation bound:\n')
        hdr2 = ['target', 'ok at 1e-12', 'T1 bound at 1e-12', 'pre width at 1e-12', 'max E at 1e-12', 'ok at 1e-24', 'T1 bound at 1e-24', 'pre width at 1e-24', 'max E at 1e-24']
        rows2 = []
        for r in raw['T_extra']:
            rows2.append(['%d' % r[0], 'yes' if r[2] > 0 else 'NO', '%.9f' % r[3], '%.2e' % r[4], '%.4f' % r[5], 'yes' if r[7] > 0 else 'NO' if not np.isnan(r[7]) else '-', '%.9f' % r[8], '%.2e' % r[9], '%.4f' % r[10]])
        L.append(table(hdr2, rows2))
    if 'C_rows' in raw:
        L.append('\n## C. The chain between eps_tail and eps_top (sound displacement accounting)\n')
        L.append('`theory_t1b_wide.cover(prep, top, eps_tail)` with `theory3a_chain.verify_fixed` (the displacement integral of Lemma 5.2 with the stage-resolved H).  "balls" = number of eps-balls, "calls" = verify calls.\n')
        hdr3 = ['target', 'log10 top', 'log10 eps_tail', 'chain needed', 'complete', 'balls', 'calls', 'narrowest ball (decades)']
        rows3 = []
        for r in raw['C_rows']:
            need = r[2] < r[1]
            rows3.append(['%d' % r[0], '%+.2f' % r[1], '%+.2f' % r[2], 'yes' if need else 'no', 'yes' if r[3] > 0 else 'NO', '%d' % r[4], '%d' % r[5], '%.4f' % r[6] if need else '-'])
        L.append(table(hdr3, rows3))
        ncomp = int(sum(r[3] > 0 for r in raw['C_rows']))
        L.append('\n**Result.** For %d of %d targets S(eps) is verified, in ball arithmetic with floating-point limit-object data, for EVERY real eps in (0, eps_top]: the tail of Part T covers (0, eps_tail] '
                 'and the chain of Part C covers [eps_tail, eps_top].\n' % (ncomp, len(raw['C_rows'])))
    return '\n'.join(L) + '\n'


PARTS_ORDER = ('U', 'S', 'T', 'C')


def load_raw():
    z = np.load(RAW, allow_pickle=False)
    return {k: z[k] for k in z.files}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--parts', default=','.join(PARTS_ORDER))
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    if a.run:
        raw = {}
        if a.parts != ap.get_default('parts') and os.path.exists(RAW):
            raw.update(load_raw())
        parts = a.parts.split(',')
        t0 = time.time()
        if 'U' in parts:
            raw.update(part_U(a.quick))
        if 'S' in parts:
            raw.update(part_S(a.quick))
        if 'T' in parts:
            raw.update(part_T(a.quick))
        if 'C' in parts:
            raw.update(part_C(raw, a.quick))
        np.savez(RAW, **raw)
        print('raw written (%.0f s)' % (time.time() - t0))
    else:
        raw = load_raw()
    with open(MD, 'w') as fh:
        fh.write(report(raw))
    print('report written:', MD)


if __name__ == '__main__':
    main()
