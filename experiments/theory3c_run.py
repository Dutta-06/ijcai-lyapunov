"""
Fully rigorous re-verification of the interval covers of Theorem T1 (round theory3, direction C).

For each of the 8 interval-verified targets of Part Q of theory_t1b.md (K <= 3, n = 1): (1) the limit object is enclosed in ball arithmetic (Krawczyk enclosure of every block critical point,
strict entry decisions; experiments/theory3c_limit.py), (2) the certificate S(eps) is verified in ball arithmetic with those enclosures, certain comparisons and outward-rounded quadrature edges
(experiments/theory3c_rigorous.py) on a chain of eps-balls from the stored top exponent down to 10^-24; if the chain does not close, the top is lowered by half decades.

Usage: python experiments/theory3c_run.py --run [--targets 0,1,..] [--budget S] | --from-raw      -> results/synthetic/theory3c.md (+ theory3c_raw.npz); environment T3C_OUT, T3C_Q.
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
import theory_t1b_wide as W        # noqa: E402  (only its cache of the float certificate's rate-function sup; Q_IDX)
import theory3c_limit as LM        # noqa: E402
import theory3c_rigorous as R3     # noqa: E402

OUT = os.environ.get('T3C_OUT', os.path.join(HERE, '..', 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory3c_raw.npz'), os.path.join(OUT, 'theory3c.md')
QRAW = os.path.join(os.environ.get('T3C_Q', os.path.join(HERE, '..', 'results', 'synthetic')), 'theory_t1b_raw.npz')
VERIFIED = (0, 1, 2, 3, 4, 5, 6, 8)           # the rows of Part Q whose chain closed (the 8 interval-verified targets)
DTOPS = (0.0, -0.5, -1.0, -1.5, -2.0, -2.5)
BUDGET = 1800.0
EV_COLS = ('n', 'k', 'kernel', 'feature', 'tau', 'tau_rad', 'p_e', 'p_e_rad', 'd0', 'gap_rel_lower', 'block_rad', 'sign_ok', 'ok', 'mu_lower', 'wall_ok', 'strict_ok')
HY_COLS = ('n', 'r', 'r_minus_d0', 'dhat_minus_off', 'w_minus_vDlo', 'vDhi_minus_w', 'vsharp_minus_theta', 'dhat', 'vD_lo', 'vD_hi', 'unique_in_Dhat')


def hyp_slacks(prep, LR, par):
    """lower bounds of the slacks of the hypotheses r > d_0, 'every stage's critical point in D-hat', theta < v_# at the parameters par (positive = holds)"""
    K = prep['K']
    As = np.asarray(prep['As'], float)
    absA = np.abs(As)
    amin, amax = arb(float(absA.min())), arb(float(absA.max()))
    cs = LR['cs']
    Delta = arb(0)
    if K > 1:
        Delta = None
        for i in range(K):
            for j in range(i + 1, K):
                d_ = abs(cs[i] - cs[j])
                Delta = d_ if Delta is None else Delta.min(d_)
    r = arb(par['r'])
    Bx = R3.rig_box(K, amax / amin, Delta, r, arb(par['Om']), par['dfrac'], par['vsig'])
    if Bx is None:
        return None
    ev = LR['events']
    s_r = min(float((r - e['d0']).lower()) for e in ev)
    s_off, s_lo, s_hi = 1e9, 1e9, 1e9
    for st in ev:
        feats = [e['feature'] for e in ev[:st['k']]]
        for i, f in enumerate(feats):
            s_off = min(s_off, float((Bx['dh'] - abs(st['c'][i] - cs[f])).lower()))
            s_lo = min(s_lo, float((abs(st['w'][i]) - amin * Bx['vD_lo']).lower()))
            s_hi = min(s_hi, float((amin * Bx['vD_hi'] - abs(st['w'][i])).lower()))
    theta = (R3.THETA_REL_B * amax).min(R3.TENTH_B * amin)
    s_th = float(((Bx['vlo'] * amin / 2) - theta).lower())
    uniq = float(all(R3.dhat_unique(LR, Bx, amin, As, st['k']) for st in ev))
    return [float(par['r']), s_r, s_off, s_lo, s_hi, s_th, float(Bx['dh'].mid()), float(Bx['vD_lo'].mid()), float(Bx['vD_hi'].mid()), uniq]


def run(targets, budget):
    warnings.filterwarnings('ignore')
    Q = np.load(QRAW, allow_pickle=False)['Q_rows']
    trows, segs_all, ev_rows, hy_rows = [], [], [], []
    for n_ in targets:
        K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n_]]
        cs, As, c0 = T.make_target(K, sep, R, seed)
        prep = T.prepare(cs, As, c0)
        t0 = time.time()
        LR = LM.limit_rig(cs, As, c0)
        mu_lo = float(LR['mu'].lower())
        for e in LR['events']:
            ev_rows.append([n_, e['k'], e['kernel'], e['feature'], float(e['tau'].mid()), float(e['tau'].rad()), float(e['p_e'].mid()), float(e['p_e'].rad()), float(e['d0'].mid()), e['gap'], e['rad'],
                            float(e['sign_ok']), float(LR['ok']), mu_lo, float(LR['checks']['wall']), float(LR['checks']['strict_entrant'] and LR['checks']['strict_feature'])])
        print('limit object', n_, (K, sep, R, seed), 'ok', LR['ok'], LR['why'], flush=True)
        old = Q[n_]
        best = None
        if LR['ok']:
            for dtop in DTOPS:
                top = old[5] + dtop
                segs, complete, cur, calls = R3.cover(prep, LR, top, W.Q_BOT, budget=budget)
                if complete:
                    best = (top, len(segs), calls, segs, cur, 1.0)
                    break
                if best is None or cur < best[4]:
                    best = (top, len(segs), calls, segs, cur, 0.0)
                print('  attempt top', top, 'stopped at', cur, 'calls', calls, flush=True)
        if best is None:
            best = (old[5], 0, 0, [], old[5], 0.0)
        top, nb, calls, segs, cur, comp = best
        wmax = max([w_[1] - w_[0] for w_ in segs]) if segs else 0.0
        trows.append([n_, K, sep, R, seed, old[4], old[5], top, nb, calls, wmax, cur, comp, time.time() - t0])
        segs_all.append(np.array([[n_, w_[0], w_[1]] for w_ in segs], float).reshape(-1, 3))
        if segs:
            _, pi = T.best_cert(prep, 10.0 ** top)
            hs = hyp_slacks(prep, LR, T.PARS[pi])
            if hs is not None:
                hy_rows.append([n_] + hs)
        print('target', n_, (K, sep, R, seed), 'old top', old[5], 'new top', top, 'lowest', cur, 'complete', comp, 'balls', nb, 'calls', calls, flush=True)
    return dict(T_rows=np.array(trows, float), T_segs=np.concatenate(segs_all) if segs_all else np.zeros((0, 3)), EV_rows=np.array(ev_rows, float), HY_rows=np.array(hy_rows, float).reshape(-1, len(HY_COLS)),
                budget=np.array([budget]), bot=np.array([W.Q_BOT]))


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(r) + ' |')
    return '\n'.join(out)


BANNER = 'SUPERSEDED by results/synthetic/theory3n.md: the displacement accounting of the certificate was corrected (see review/REVIEW_roundT3AC.md, REVIEW_roundT3N.md); numbers below use the legacy accounting and are kept as history'


def report(raw):
    Tr, Ev, Hy = raw['T_rows'], raw['EV_rows'], raw['HY_rows']
    L = [BANNER + '\n',
         '# Interval verification of Theorem T1 for the 8 verified targets, with the legacy displacement accounting; superseded by theory3n.md (generated by experiments/theory3c_run.py)\n',
         'Round theory3, direction C.  Everything that `theory_t1b_rigorous.py` took from floating point is now enclosed: the **limit object** (each block critical point by a Krawczyk operator on the '
         'reduced critical system, with existence and uniqueness in the enclosing box; the entry times, the duals at the breakpoints, the rates, the margins, the entrant/feature assignment by a '
         'strict-argmin decision on the enclosures) and every **final comparison** (certain arb comparisons, no float() of a ball, outward-rounded quadrature edges and eps-ball ends, '
         'no undecided ball replaced by 0).  The hypotheses r > d_0, sgn p_e = sgn A_f, every stage\'s critical point in the relaxation set D-hat, and theta < v_# are checked on the enclosures inside '
         '`verify_rig`; in addition the Krawczyk operator is applied to the relaxation set D-hat itself (a box containing it) for every stage, which proves that the block has exactly ONE critical point '
         'in all of D-hat (existence and uniqueness), not only in the small enclosing box.  The target (centres, amplitudes, dictionary) is the exact binary value of its doubles.\n',
         'Not enclosed: Theorem K\'s statement that the block critical point is the GLOBAL minimiser of the fit under the sign constraint (what T1 uses, and what is proved here, is that the block has '
         'exactly one critical point in the relaxation set D-hat and that it is the enclosed one; the limit object is the recursion with these critical points), the monotonicity of S in eps beyond '
         'the chains below, and the soundness of python-flint/arb itself.\n',
         '## Limit object (ball arithmetic, 220 bits; entry events)\n']
    rows = []
    for e in Ev:
        rows.append(['%d' % e[0], '%d' % e[1], '%d' % e[2], '%d' % e[3], '%.15f' % e[4], '%.1e' % e[5], '%.12f' % e[6], '%.3f' % e[8], '%.4f' % e[9], '%.0e' % e[10], 'yes' if e[11] > 0 else 'NO',
                     'yes' if (e[12] > 0 and e[14] > 0 and e[15] > 0) else 'NO'])
    L.append(table(['target', 'k', 'entrant kernel', 'feature', 'tau_k', 'tau_k rad', 'p_k', 'd_0 (sigma)', 'relative gap to the 2nd hitting time (lower bound)', 'block box radius', 'sgn p = sgn A', 'strict & wall <1'], rows))
    L.append('\n"relative gap" = (lower end of the second-smallest hitting time - upper end of the entrant\'s) / entrant\'s: strictly positive means the entrant, hence the whole recursion, is decided in '
             'ball arithmetic.  "block box radius" = the largest radius of the enclosed block critical point (w, c) after contraction.\n')
    L.append('## Hypotheses on the enclosures (slacks = lower bounds of the margin; positive = holds)\n')
    rows = []
    for h in Hy:
        rows.append(['%d' % h[0], '%.2f' % h[1], '%.3f' % h[2], '%.3e' % h[3], '%.3e' % h[4], '%.3e' % h[5], '%.3e' % h[6], '%.4f' % h[7], '%.4f' % h[8], '%.3f' % h[9], 'yes' if h[10] > 0 else 'NO'])
    L.append(table(['target', 'r', 'r - max d_0', 'd-hat - max offset (all stages)', 'min (|w| - a_min v_D^-)', 'min (a_min v_D^+ - |w|)', 'v_# - theta', 'd-hat', 'v_D^- / a_min', 'v_D^+ / a_min',
                    'one critical point in all of D-hat (Krawczyk on D-hat, every stage)'], rows))
    L.append('\n## Interval covers of eps (every real eps in [10^-24, 10^top])\n')
    L.append('Candidate bounds D_j from the floating-point certificate (x1.1 ... x4), then every ingredient re-evaluated in ball arithmetic with the enclosed limit object.  The chain of balls has outward-rounded '
             'ends, consecutive balls overlap.  "old top" = the top of Part Q of theory_t1b.md (floating-point limit-object data); when the chain does not close from the old top, the top is lowered by half '
             'decades (attempt budget %g s per top).\n' % raw['budget'][0])
    rows = []
    for t in Tr:
        rows.append(['%d' % t[0], '%d' % t[1], '%s' % ('%g' % t[2] if t[1] > 1 else '-'), '%g' % t[3], '%d' % t[4], '%+.2f' % t[6], '%+.2f' % t[7], '%d' % t[8], '%d' % t[9], '%.3f' % t[10], '%+.2f' % t[11],
                     'complete' if t[12] > 0 else 'INCOMPLETE'])
    L.append(table(['target', 'K', 'sep', 'R', 'seed', 'old top (log10 eps)', 'NEW top', 'balls', 'verify calls', 'widest ball (decades)', 'lowest exponent', 'status'], rows))
    ncomp = int(sum(1 for t in Tr if t[12] > 0))
    nsame = int(sum(1 for t in Tr if t[12] > 0 and abs(t[6] - t[7]) < 1e-9))
    L.append('\n%d of %d targets: S(eps) is verified, in ball arithmetic with the enclosed limit object and certain comparisons, for every real eps in [10^-24, 10^NEW top]; the top is unchanged for %d of them.\n' % (ncomp, len(Tr), nsame))
    return '\n'.join(L)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--targets', default=','.join(str(v) for v in VERIFIED))
    ap.add_argument('--budget', type=float, default=BUDGET)
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    if a.run:
        raw = run([int(x) for x in a.targets.split(',')], a.budget)
        np.savez(RAW, **raw)
    else:
        z = np.load(RAW, allow_pickle=False)
        raw = {k: z[k] for k in z.files}
    with open(MD, 'w') as fh:
        fh.write(report(raw))
    print('report written:', MD)


if __name__ == '__main__':
    main()
