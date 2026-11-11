"""
The entering orbit of the coupled entry layer: runner and report (round theory3, direction Y; library theory3y_entry.py, the separated-regime part theory3y_sep.py).

Parts:
  A  the hypotheses of Theorems E1-E4 (Prop. A's box conditions (L), (C), (N-hat) with Omega = 0 at a radius r > d0) on every entry layer of the 8 targets of theory3n and of
     8 further random targets, and the constants (box, level w_#, T_2, T_3, decay rate rho);
  B  the layer flows from the DORMANT state of the network (x_e = 0, w_e = 0 exactly) at eps = 1e-6 ... 1e-12: the box, the monotone entrant, the hitting time of the relaxation set
     against T_2 + T_3, the exponential decay against rho, the time to enter D*;
  C  the limit object's layer (eps = 0): flows from w_e = v_0 a_min, v_0 = 1e-3 ... 1e-12, time-aligned at the first crossing of w_#, against the reference flow of v_0 = 1e-14: the aligned
     distance is the Cauchy property of the entering orbit (the same checks as in B on every flow);
  D  the network's layer flows of part B against the reference flow of eps = 0: the aligned distance as eps -> 0;
  S  the separated regime (K = 2, n = 1, equal amplitudes): the closed-form sufficient condition of Theorem S6' of the note against the exact recursion of Theorem K (theory3y_sep.py).

Usage: python experiments/theory3y_run.py --run [--quick] [--parts A,B,C,D,S] | --from-raw      -> results/synthetic/theory3y.md (+ theory3y_raw.npz); environment T3Y_OUT.
Deterministic (no random draw except the fixed generator of the random targets); the report is regenerated from the raw alone, byte for byte.
"""
import argparse
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import theory3y_entry as Y          # noqa: E402

OUT = os.environ.get('T3Y_OUT', os.path.join(HERE, '..', 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory3y_raw.npz'), os.path.join(OUT, 'theory3y.md')
A_COLS = ['K', 'k', 'Delta', 'R', 'd0', 'amin', 'r', 'vlo', 'vbar', 'Ev', 'Ec', 'dtilde', 'dhat', 'nu', 'Hup', 'vD_lo', 'regime', 'sign_ok']
CHK = ['box_d', 'box_vhi', 'box_vlo_block', 'mono_entrant', 'entrant_ge_vlo_after', 'block_in_hat', 'tail', 'decay']
VEC = ['t_hash', 't_hat', 'That', 'ratio_hat', 'ratio_decay', 't_star', 't_star_bd', 'err_end', 'err_hat', 'dmax_over_r', 'vmax_over_vbar', 'margin_star']
BCOLS = ['w_hash', 'T2', 'T3', 'That', 'rho', 'pref', 'm_lo', 'mV', 'mVp', 'kd']


def label(L):
    t = L['target']
    return '(%d,%g,%g,%d)/%d' % (t[0], t[1], t[2], t[3], L['k'] + 1)


def all_layers(quick):
    """All entry layers: set D (the 8 targets of theory3n), set R (8 random targets), set S (stress layers).  quick: the first layer of D and the first stress layer only (the tests)."""
    if quick:
        return [(Y.layers_of(Y.D_TARGETS[0])[0], 'D'), (Y.stress_layers()[0], 'S')]
    out = []
    for t in Y.D_TARGETS:
        out += [(L, 'D') for L in Y.layers_of(t)]
    for t in Y.random_targets(8):
        out += [(L, 'R') for L in Y.layers_of(t)]
    out += [(L, 'S') for L in Y.stress_layers()]
    return out


def run(quick=False, parts='A,B,C,D'):
    parts = set(parts.split(','))
    layers = all_layers(quick)
    nL = len(layers)
    eps_l = [1e-6, 1e-10] if quick else Y.EPS_LIST
    v0_l = [1e-3, 1e-6] if quick else Y.V0_LIST
    raw = {}
    boxes = [Y.layer_box(L) for (L, _) in layers]
    A = np.full((nL, len(A_COLS)), np.nan)
    for i, ((L, s), b) in enumerate(zip(layers, boxes)):
        row = dict(K=L['K'], k=L['k'] + 1, Delta=L['Delta'] if np.isfinite(L['Delta']) else 0.0, R=L['R'], d0=L['d0'], amin=L['amin'], regime=float(b is not None), sign_ok=float(L['sign_ok']))
        if b is not None:
            row.update({k: b[k] for k in ('r', 'vlo', 'vbar', 'Ev', 'Ec', 'dtilde', 'dhat', 'nu', 'Hup', 'vD_lo')})
        A[i] = [row.get(c, np.nan) for c in A_COLS]
    raw['A_lab'] = np.array([label(L) for (L, _) in layers])
    raw['A_set'] = np.array([s for (_, s) in layers])
    raw['A_num'] = A
    raw['eps_list'] = np.array(eps_l)
    raw['v0_list'] = np.array(v0_l)
    raw['quick'] = np.array(int(quick))
    ne, nv = len(eps_l), len(v0_l)
    B_ok = np.zeros((nL, ne, len(CHK)), bool)
    B_vec = np.full((nL, ne, len(VEC)), np.nan)
    B_bd = np.full((nL, ne, len(BCOLS)), np.nan)
    C_ok = np.zeros((nL, nv, len(CHK)), bool)
    C_vec = np.full((nL, nv, len(VEC)), np.nan)
    C_dist = np.full((nL, nv), np.nan)
    C_orb = np.full((nL, nv), np.nan)
    D_orb = np.full((nL, ne), np.nan)
    D_dist = np.full((nL, ne), np.nan)
    D_win = np.full(nL, np.nan)
    C_win = np.full((nL, nv), np.nan)
    t0 = time.time()
    for i, (L, s) in enumerate(layers):
        b = boxes[i]
        if b is None:
            continue
        bd0 = Y.bounds(L, b, 0.0)
        ref = Y.integrate_layer(L, 0.0, v0=Y.V0_REF, w_hash=bd0['w_hash'])
        trB = {}
        if 'B' in parts or 'D' in parts:
            for j, eps in enumerate(eps_l):
                bd = Y.bounds(L, b, eps)
                tr = Y.integrate_layer(L, eps, w_hash=bd['w_hash'])
                trB[j] = tr
                ok, vec = Y.check_layer(L, b, bd, tr, eps)
                B_ok[i, j] = [ok[c] for c in CHK]
                B_vec[i, j] = [vec[c] for c in VEC]
                B_bd[i, j] = [bd[c] for c in BCOLS]
        if 'C' in parts:
            for j, v0 in enumerate(v0_l):
                tr = Y.integrate_layer(L, 0.0, v0=v0, w_hash=bd0['w_hash'])
                ok, vec = Y.check_layer(L, b, bd0, tr, 0.0)
                C_ok[i, j] = [ok[c] for c in CHK]
                C_vec[i, j] = [vec[c] for c in VEC]
                lo = 0.9 * tr['t_hash']
                C_win[i, j] = lo
                C_dist[i, j] = Y.aligned_distance(tr, ref, L, 0.0, 0.0, -lo, Y.S2)
                C_orb[i, j] = Y.orbit_distance(tr, ref, L, 0.0, 0.0, v0)
        if 'D' in parts:
            lo = 0.9 * min(trB[j]['t_hash'] for j in trB)
            D_win[i] = lo
            for j, eps in enumerate(eps_l):
                D_dist[i, j] = Y.aligned_distance(trB[j], ref, L, eps, 0.0, -lo, Y.S2)
                vs_ = float(np.interp(-0.9 * trB[j]['t_hash'], trB[j]['s'], Y.derived(L, trB[j], eps, None)[0][-1]))
                D_orb[i, j] = Y.orbit_distance(trB[j], ref, L, eps, 0.0, max(vs_, 1e-12))
        print('layer %d/%d %s  %.0fs' % (i + 1, nL, label(L), time.time() - t0), flush=True)
    raw.update(B_ok=B_ok, B_vec=B_vec, B_bd=B_bd, C_ok=C_ok, C_vec=C_vec, C_dist=C_dist, C_orb=C_orb, D_orb=D_orb, D_dist=D_dist, D_win=D_win, C_win=C_win)
    if 'S' in parts:
        import theory3y_sep as SP
        raw.update(SP.run_part(quick))
    return raw


def load_raw():
    z = np.load(RAW, allow_pickle=False)
    return {k: z[k] for k in z.files}


def save_raw(raw):
    os.makedirs(OUT, exist_ok=True)
    np.savez_compressed(RAW, **raw)


def fe(x):
    return 'nan' if not np.isfinite(x) else '%.3g' % x


def fmax(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(a.max()) if len(a) else float('nan')


def fmin(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(a.min()) if len(a) else float('nan')


def fmed(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if len(a) else float('nan')


def _slope(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(y) & (y > 0)
    return float(np.polyfit(np.log10(x[ok]), np.log10(y[ok]), 1)[0]) if ok.sum() >= 2 else float("nan")


def table(head, rows):
    out = ['| ' + ' | '.join(head) + ' |', '|' + '|'.join(['---'] * len(head)) + '|']
    out += ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows]
    return out


def report(raw):
    L = []
    A = raw['A_num']
    ac = {c: A[:, i] for i, c in enumerate(A_COLS)}
    lab, st = raw['A_lab'], raw['A_set']
    eps_l, v0_l = raw['eps_list'], raw['v0_list']
    reg = ac['regime'] == 1
    nL = len(lab)
    L += ['# The entering orbit of the coupled entry layer (generated by experiments/theory3y_run.py)', '',
          'Round theory3, direction Y.  Statements and proofs: the round note (outside the repository), Theorems E0-E4; this report is the numerical side: the layer flow of the block and the entrant '
          'integrated from the dormant state (library `theory3y_entry.py`), checked against the constants of the theorems.  Continuum Gram (Z = 1), n = 1, Gaussian features, Omega = 0 '
          '(no dormant kernel in the layer), the repository law (m, G, gamma, a_max) = (12, 1, 0.1, 1).', '',
          '## A. The hypotheses and the constants', '',
          'Every entry layer (k-th entry) of the 8 targets of theory3n (set D: `(K,sep,R,seed)/k`) and of 8 further targets drawn from a fixed generator (set R: K in {2, 3}, separation 5.5-9.5 widths, '
          'amplitude ratio 1-2.4) and 16 stress layers (set S: second layers of K = 2 targets built at Delta = 4.9, 5.2, 5.6 widths, R = 1, 2.4, the entrant 0.40 widths from its feature toward the formed unit; label `(2,Delta,R,code)`, code = 100 weak feature first + 10 opposite signs + 1).  `regime` = the box conditions (L) v_lo > 0, (C) psi(r) > E_c, (N-hat) nu = lambda_own(d-hat) - H_c > 0 with positive lower faces of the relaxation set, in ball arithmetic '
          '(exp_theory_t1.box_t1, midpoints printed), at the smallest grid radius r > d0 (step 0.05).  Layers outside the regime are not covered by the theorems and are not integrated.', '']
    rows = []
    for i in range(nL):
        rows.append([lab[i], st[i], int(ac['K'][i]), fe(ac['Delta'][i]), fe(ac['R'][i]), fe(ac['d0'][i]), int(reg[i]),
                     fe(ac['r'][i]), fe(ac['vlo'][i]), fe(ac['vbar'][i]), fe(ac['Ev'][i]), fe(ac['Ec'][i]), fe(ac['dtilde'][i]), fe(ac['dhat'][i]), fe(ac['nu'][i])])
    L += table(['layer', 'set', 'K', 'Delta', 'R', 'd0', 'in regime', 'r', 'v_lo', 'v_bar', 'E_v', 'E_c', 'd-tilde', 'd-hat', 'nu'], rows)
    L += ['', 'Layers in the regime: **%d of %d** (set D: %d of %d, set R: %d of %d, set S: %d of %d); `sign_ok` (the entrant has the sign of its feature, Thm K (a)): %d of %d.' % (
        int(reg.sum()), nL, int(reg[st == 'D'].sum()), int((st == 'D').sum()), int(reg[st == 'R'].sum()), int((st == 'R').sum()), int(reg[st == 'S'].sum()), int((st == 'S').sum()),
        int((ac['sign_ok'] == 1).sum()), nL), '']
    ii = np.flatnonzero(reg)
    if len(ii) == 0:
        return '\n'.join(L) + '\n'
    # ---- B
    B_ok, B_vec, B_bd = raw['B_ok'], raw['B_vec'], raw['B_bd']
    vi = {c: i for i, c in enumerate(VEC)}
    bi = {c: i for i, c in enumerate(BCOLS)}
    L += ['## B. Layer flows from the dormant state (eps > 0)', '',
          'The network\'s own flow of the block and the entrant: the entrant at x = 0 (w = 0 exactly), the block at the previous critical point, integrated (Radau, rtol 1e-11) until the entrant reaches the level '
          'w_# = a_min v_lo / 2 (the first crossing, time t_#) and 4000 time units further.  Times below are measured from t_# (a fast time; t_# itself is the pre-layer, of order eps^(-1/2)).  '
          '`checks` = all of: the box |d| < r and v <= v_bar, the block above v_lo, the entrant nondecreasing until it reaches v_lo and above v_lo afterwards, the block inside the relaxation set D-hat '
          'throughout, the tail of Lemma E0 (w_e(s)^(-1/3) >= w_#^(-1/3) + c_M w_# |s|/3 for s <= 0), and the decay of Theorem E4 (|y - y_hat| <= sqrt(H_up/nu) |y(t_hat) - y_hat| exp(-rho (s - t_hat)) for s >= t_hat, slack 1e-9).  `t_hat` = the first time from which every unit stays in D-hat; '
          '`bound` = T_2 + T_3 (balls); `t_D*` = the time from which every unit stays in D*, its bound t_hat + ln(sqrt(H_up/nu) diam(D-hat) / margin)/rho (Theorem E3(iii), a priori); nan for K = 1 (d-tilde = 0: D* is the face d = 0, '
          'entered only asymptotically).', '']
    rows = []
    allok = True
    for j, e in enumerate(eps_l):
        okj = B_ok[ii, j].all(axis=1)
        vv = B_vec[ii, j]
        rh = vv[:, vi['ratio_hat']]
        rd = vv[:, vi['ratio_decay']]
        with np.errstate(invalid='ignore', divide='ignore'):
            rs = vv[:, vi['t_star']] / vv[:, vi['t_star_bd']]
        th = vv[:, vi['t_hash']]
        rows.append(['%g' % e, len(ii), int(okj.sum()), fe(fmax(rh)), fe(fmax(rd)), fe(fmax(rs)), fe(fmax(vv[:, vi['err_end']])), fe(th.min()), fe(th.max()),
                     fe(fmax(vv[:, vi['dmax_over_r']]))])
        allok = allok and int(okj.sum()) == len(ii)
    L += table(['eps', 'layers', 'all checks', 'max t_hat/bound', 'max decay ratio', 'max t_D*/bound', 'max final error', 'min t_#', 'max t_#', 'max |d|/r'], rows)
    L += ['', 'Violations of any check or bound over %d layers x %d values of eps: **%d**.' % (len(ii), len(eps_l), int(sum(len(ii) - int(B_ok[ii, j].all(axis=1).sum()) for j in range(len(eps_l))))), '']
    # per-layer table at the smallest eps
    j = len(eps_l) - 1
    rows = []
    for i in ii:
        vv = B_vec[i, j]
        rows.append([lab[i], st[i], fe(vv[vi['t_hash']]), fe(vv[vi['t_hat']]), fe(B_bd[i, j, bi['That']]), fe(vv[vi['ratio_hat']]), fe(vv[vi['ratio_decay']]),
                     fe(vv[vi['t_star']]), fe(vv[vi['t_star_bd']]), fe(B_bd[i, j, bi['rho']]), fe(vv[vi['err_hat']])])
    L += ['Per layer at eps = %g:' % eps_l[j], '']
    L += table(['layer', 'set', 't_#', 't_hat', 'T_2 + T_3', 'ratio', 'decay ratio', 't_D*', 'bound', 'rho', '|y(t_hat) - y_hat|'], rows)
    L += ['']
    # ---- C
    C_ok, C_vec, C_dist = raw['C_ok'], raw['C_vec'], raw['C_dist']
    L += ['## C. The limit object\'s layer: flows from v_0 -> 0 (eps = 0), time-aligned', '',
          'The flow of the critical curve (eps = 0) started at w_e = v_0 a_min with the block at the previous critical point, aligned at the first crossing of w_#.  `distance` = sup over the aligned window '
          's in [-0.9 t_#(v_0), 4000] of the distance in (v\', d) to the reference flow v_0 = 1e-14 (a Cauchy sequence if the entering orbit is unique up to a time shift).  The checks of part B are '
          'applied to every flow (the bound is that of eps = 0).', '']
    rows = []
    for j, v0 in enumerate(v0_l):
        okj = C_ok[ii, j].all(axis=1)
        vv = C_vec[ii, j]
        dd = C_dist[ii, j]
        oo = raw['C_orb'][ii, j]
        rows.append(['%g' % v0, len(ii), int(okj.sum()), fe(fmax(vv[:, vi['ratio_hat']])), fe(fmax(vv[:, vi['ratio_decay']])), fe(fmax(dd)), fe(fmed(dd)), fe(fmax(oo)), fe(fmed(oo)),
                     fe(fmin(vv[:, vi['t_hash']])), fe(fmax(vv[:, vi['t_hash']]))])
    L += table(['v_0', 'layers', 'all checks', 'max t_hat/bound', 'max decay ratio', 'max aligned distance', 'median aligned', 'max orbit distance', 'median orbit', 'min t_#', 'max t_#'], rows)
    sa = _slope(v0_l[3:-1], np.array([fmax(C_dist[ii, j]) for j in range(3, len(v0_l) - 1)]))
    so = _slope(v0_l[3:-1], np.array([fmax(raw['C_orb'][ii, j]) for j in range(3, len(v0_l) - 1)]))
    th_ = raw['C_vec'][ii][:, :, vi['t_hash']]
    j6, j9 = (list(v0_l).index(1e-6), list(v0_l).index(1e-9)) if (1e-6 in list(v0_l) and 1e-9 in list(v0_l)) else (0, len(v0_l) - 1)
    L += ['', 'Log-log slopes against v_0 over v_0 = %g ... %g (the last row, closest to the reference, is left out): aligned distance %s, orbit distance %s (both about 2/3: the finite-v_0 flow starts with '
          'the entrant\'s centre at its lattice position, so it lacks the centre drift int_(-inf)^(-t_#(v_0)) of order mu_c v_0 t_#(v_0) = O(v_0^(2/3)); the distance at v_0 = 1e-12 is the sampling floor of '
          'the linear interpolation, about 3e-6).  The ratio of the median t_# between v_0 = %g and %g is %s (v_0^(-1/3) gives 10: Lemma E0 gives the lower bound t_# >= (v_0^(-1/3) - w_#^(-1/3))/kappa_0).' % (
              v0_l[min(3, len(v0_l) - 1)], v0_l[max(len(v0_l) - 2, 0)], fe(sa), fe(so), v0_l[j9], v0_l[j6], fe(float(np.median(th_[:, j9] / th_[:, j6])))), '']
    # ---- D
    D_dist = raw['D_dist']
    L += ['## D. The network\'s layer flows against the limit object\'s (eps -> 0)', '',
          'Distance (v\', d), aligned at t_#, sup over s in [-0.9 min t_#, 4000], between the layer flow of eps (part B) and the reference flow of eps = 0 (v_0 = 1e-14).', '']
    rows = []
    for j, e in enumerate(eps_l):
        dd = D_dist[ii, j]
        oo = raw['D_orb'][ii, j]
        rows.append(['%g' % e, fe(fmax(dd)), fe(fmed(dd)), fe(fmax(oo)), fe(fmed(oo))])
    L += table(['eps', 'max aligned distance', 'median aligned', 'max orbit distance', 'median orbit'], rows)
    sl = _slope(eps_l[:2], np.array([fmax(D_dist[ii, j]) for j in range(min(2, len(eps_l)))])) if len(eps_l) >= 2 else float('nan')
    L += ['', 'Log-log slope of the maximal aligned distance against eps over the first two values: %s (the distance is proportional to eps until it reaches the sampling floor of about 3e-6; '
          'the relative difference of w_eps and w_0 at equal x is eps/(2 v_s sinh^2 x)).' % fe(sl), '']
    if 'S_rows' in raw:
        import theory3y_sep as SP
        L += SP.report_part(raw)
    L += ['## What this report does and does not show', '',
          '* It integrates the layer flow itself (not a least-squares polish) from the dormant state at eps = 1e-6 ... 1e-12 and from v_0 = 1e-3 ... 1e-12 at eps = 0, on the layers in the regime, '
          'and checks the box, the monotone entry, the hitting time of the relaxation set against T_2 + T_3 and the decay against rho.  A check is not a proof: the proofs are in the note; '
          'a violation here would be a bug of a theorem.', '* The constants are those of Prop. A / THEORY_T1.md evaluated in ball arithmetic (midpoints printed), T_2, T_3, rho in balls from those midpoints; '
          'the flows are floating point.',
          '* Layers outside the regime (Delta below Prop. A\'s threshold or d0 >= the admissible radius) and K > 3, n > 1, sampled Grams, unequal masses, noise are not covered.', '']
    return '\n'.join(L) + '\n'


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--parts', default='A,B,C,D,S')
    a = ap.parse_args(argv)
    if a.run:
        raw = run(a.quick, a.parts)
        save_raw(raw)
    elif a.from_raw:
        raw = load_raw()
    else:
        ap.error('give --run or --from-raw')
    os.makedirs(OUT, exist_ok=True)
    with open(MD, 'w') as f:
        f.write(report(raw))
    print('wrote', MD)


if __name__ == '__main__':
    main()
