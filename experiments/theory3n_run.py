"""
Theorem T1, corrected and merged verifications (round theory3, direction N; supersedes theory_t1/theory_t1b scans of `cert`, theory3a and theory3c as history).

The referee of round T3 (review/REVIEW_roundT3AC.md) found that the displacement accounting of the certificate S(eps) was wrong in every verifier ((f): the term H^(st(t)) int_0^t B must be the
stage-resolved integral of H^(st(t')) B; (g): the pre-layer coefficient of an entrant whose window has not started enters d_c R_l), plus minor points (stage maximum of the old units' coefficients,
kappa by the smaller lower end, a float gate in `inv_of`).  `theory3n_system.py` is the one corrected system; `theory3n_tail.py` ports the eps-free small-eps tail of round theory3 A to the
ENCLOSED limit object of round theory3 C; `theory3n_kernels.py` provides the Gram-tail correction below the table start (E-2 of review/REVIEW_roundT3E.md; not used by any result of the report).  Parts:

  F  floating-point certificate: `cert` / `cert2d` (n = 1, K <= 3; K = 4..6; planar) re-run on ALL grid pairs of theory_t1.md / theory_t1b.md with the corrected system; eps_0, counts, violations,
     the effect of each repair on dcmax and E at the pairs that were informative before; the network is re-run on every pair whose status changed;
  B  the 364 stored eps-balls of theory3c re-verified with the legacy and the corrected ball verifier;
  E  the 8 targets, enclosed limit object: chain of eps-balls from the top down to 10^EHI with `verify_rig_n` and ONE tail ball (0, 10^EHI] with `verify_tail_rig`;
  T  how far up one tail ball reaches on the enclosed data (comparison with theory3a.md);
  N  the bounds dcmax and E against the integrated network at eps_0 and one decade below;
  G  the full chain [10^-24, top] of theory3c.md with the corrected verifier.

Usage: python experiments/theory3n_run.py --run [--quick] [--parts F,B,E] [--targets 0,1,..] | --from-raw      -> results/synthetic/theory3n.md (+ theory3n_raw.npz); environment T3N_OUT.
Part files are merged into the raw (a --run of one part keeps the others).
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
import theory_t1b_planar as PL     # noqa: E402
import theory3c_limit as LM        # noqa: E402
import theory3c_rigorous as R3     # noqa: E402
import theory3c_run as C3          # noqa: E402
import theory3n_system as S        # noqa: E402
import theory3n_tail as TL         # noqa: E402

OUT = os.environ.get('T3N_OUT', os.path.join(HERE, '..', 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'theory3n_raw.npz'), os.path.join(OUT, 'theory3n.md')
SRC = os.environ.get('T3N_SRC', os.path.join(HERE, '..', 'results', 'synthetic'))           # where the earlier raw files are read from
VERIFIED = (0, 1, 2, 3, 4, 5, 6, 8)
EHI = -12.0                                   # the chain ends and the tail ball starts at 10^EHI
DTOPS = (0.0, -0.5, -1.0, -1.5, -2.0, -2.5)
FACS = (1.1, 1.5, 2.0, 4.0)
BUDGET = 3000.0


def _load(name):
    z = np.load(os.path.join(SRC, name), allow_pickle=False)
    return {k: z[k] for k in z.files}


def thr_first(ok, ex):
    return T.thr_first(list(ok), ex)


# ----------------------------------------------------------------------------------------------------------------------------------
# F. floating-point certificate on all grid pairs
# ----------------------------------------------------------------------------------------------------------------------------------
def _sets(quick):
    r1, r2 = _load('theory_t1_raw.npz'), _load('theory_t1b_raw.npz')
    sets = {'C': dict(targets=T.C_TARGETS, ex=r1['C_ex'], okc=r1['C_okc'], okm=r1['C_okm'],
                      mk=lambda t: T.prepare(*T.make_target(*t)), leg=T.cert, new=S.cert_n, sim=lambda p, e: T.simulate(p, e, nlocal=300), fol=T.follows),
            'K': dict(targets=W.K_TARGETS, ex=r2['K_ex'], okc=r2['K_okc'], okm=r2['K_okm'],
                      mk=lambda t: T.prepare(*W.make_target_k(*t)), leg=T.cert, new=S.cert_n, sim=lambda p, e: T.simulate(p, e, nlocal=300), fol=T.follows),
            'P': dict(targets=W.P_TARGETS, ex=r2['P_ex'], okc=r2['P_okc'], okm=r2['P_okm'],
                      mk=lambda t: PL.prepare2(*PL.make_target2(*t)), leg=PL.cert2d, new=S.cert2d_n, sim=lambda p, e: PL.simulate2(p, e, nlocal=300), fol=lambda p, e, s, u: PL.follows2(p, e, s, u))}
    if quick:
        sets = {'C': dict(sets['C'], targets=[T.C_TARGETS[7]], okc=sets['C']['okc'][7:8], okm=sets['C']['okm'][7:8]),
                'P': dict(sets['P'], targets=[W.P_TARGETS[0]], okc=sets['P']['okc'][0:1], okm=sets['P']['okm'][0:1])}
    return sets


def _best(cf, prep, eps):
    o = None
    for i, par in enumerate(T.PARS):
        o = cf(prep, eps, par)
        if o['ok']:
            return o, i
    return o, -1


def run_F(quick=False, log=print):
    warnings.filterwarnings('ignore')
    raw = {}
    nrows = []
    for key, st in _sets(quick).items():
        ex = st['ex']
        nt, ne = len(st['targets']), len(ex)
        A_ = {nm: np.full((nt, ne), np.nan) for nm in ('okc_new', 'pi_old', 'pi_new', 'dc_old', 'dc_same', 'E_old', 'E_same', 'D_new', 'dc_new')}
        for ti, t in enumerate(st['targets']):
            t0 = time.time()
            prep = st['mk'](t)
            for ei, e in enumerate(ex):
                eps = 10.0 ** e
                on, pn = _best(st['new'], prep, eps)
                A_['okc_new'][ti, ei] = float(on['ok'])
                A_['pi_new'][ti, ei] = pn
                if on['ok']:
                    A_['D_new'][ti, ei] = on['D']
                if on.get('dcmax') is not None:
                    A_['dc_new'][ti, ei] = on['dcmax']
                old_ok = bool(st['okc'][ti, ei] > 0.5)
                if old_ok:
                    ol, pl = _best(st['leg'], prep, eps)
                    A_['pi_old'][ti, ei] = pl
                    if ol['ok']:
                        o3 = st['new'](prep, eps, T.PARS[pl])
                        A_['dc_old'][ti, ei], A_['dc_same'][ti, ei] = ol['dcmax'], o3.get('dcmax', np.nan)
                        A_['E_old'][ti, ei], A_['E_same'][ti, ei] = ol['Dnew'], o3.get('Dnew', np.nan)
                if old_ok != bool(on['ok']):
                    sol, U = st['sim'](prep, eps)
                    m = st['fol'](prep, eps, sol, U)
                    m = m[0] if isinstance(m, tuple) else m
                    nrows.append([{'C': 0, 'K': 1, 'P': 2}[key], ti, e, float(old_ok), float(on['ok']), float(st['okm'][ti, ei]), float(bool(m))])
            log('F', key, ti, t, 'informative old %d new %d (%.0fs)' % (int((st['okc'][ti] > 0.5).sum()), int(A_['okc_new'][ti].sum()), time.time() - t0), flush=True)
        raw['F%s_targets' % key] = np.array(st['targets'], float)
        raw['F%s_ex' % key] = np.array(ex, float)
        raw['F%s_okc_old' % key] = np.array(st['okc'], float)
        raw['F%s_okm' % key] = np.array(st['okm'], float)
        for nm, v in A_.items():
            raw['F%s_%s' % (key, nm)] = v
    raw['F_net'] = np.array(nrows, float).reshape(-1, 7)
    return raw


# ----------------------------------------------------------------------------------------------------------------------------------
# B. the stored balls of theory3c, legacy against corrected
# ----------------------------------------------------------------------------------------------------------------------------------
def run_B(targets, quick=False, log=print):
    warnings.filterwarnings('ignore')
    z = _load('theory3c_raw.npz')
    segs = z['T_segs']
    rows = []
    for n in targets:
        K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
        cs, As, c0 = T.make_target(K, sep, R, seed)
        prep = T.prepare(cs, As, c0)
        LR = LM.limit_rig(cs, As, c0)
        ent = [e['kernel'] for e in LR['events']]
        sg = sorted([(b, a) for (m, a, b) in segs if int(m) == n], reverse=True)
        if quick:
            sg = sg[:2]
        t0 = time.time()
        passL = passN = same = failN = 0
        rdc, rdr = [], []
        wL = wN = 0.0
        for (hi, lo_) in sg:
            o, pi = S.best_cert_n(prep, 10.0 ** hi)
            par = T.PARS[pi if pi >= 0 else 0]
            ball = R3.eps_ball(lo_, hi)
            fL = fN = None
            for fac in FACS:
                Dv = [1.0 if j in ent else max(fac * o['Etot'][j], 1.0) for j in range(len(c0))]
                res = {}
                for key, vf in (('L', R3.verify_rig), ('N', S.verify_rig_n)):
                    try:
                        res[key] = vf(prep, LR, ball, par, Dv)
                    except (ValueError, ZeroDivisionError):
                        res[key] = (False, {})
                if res['L'][0] and fL is None:
                    fL = fac
                if res['N'][0] and fN is None:
                    fN = fac
                    if 'dcmax' in res['L'][1] and 'dcmax' in res['N'][1] and res['L'][1]['dcmax'] > 0:           # same candidate bounds D_j for both verifiers
                        rdc.append(res['N'][1]['dcmax'] / res['L'][1]['dcmax'])
                        rdr.append(max(res['N'][1]['drift']) / max(res['L'][1]['drift']))
                    wN = max(wN, res['N'][1]['worst_ratio'])
                    if res['L'][0]:
                        wL = max(wL, res['L'][1]['worst_ratio'])
                if fL is not None and fN is not None:
                    break
            passL += fL is not None
            passN += fN is not None
            same += (fL is not None and fN == fL)
            failN += fN is None
        rows.append([n, len(sg), passL, passN, same, failN, min(rdc) if rdc else np.nan, float(np.median(rdc)) if rdc else np.nan, max(rdc) if rdc else np.nan, max(rdr) if rdr else np.nan, wL, wN, time.time() - t0])
        log('B', n, rows[-1], flush=True)
    return dict(B_rows=np.array(rows, float).reshape(-1, 13))


# ----------------------------------------------------------------------------------------------------------------------------------
# E. enclosed chain + tail ball
# ----------------------------------------------------------------------------------------------------------------------------------
def run_E(targets, budget=BUDGET, quick=False, log=print):
    warnings.filterwarnings('ignore')
    z = _load('theory3c_raw.npz')
    old = {int(r[0]): r for r in z['T_rows']}
    trows, segs_all, ev_rows, hy_rows, extra = [], [], [], [], {}
    for n in targets:
        K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
        cs, As, c0 = T.make_target(K, sep, R, seed)
        prep = T.prepare(cs, As, c0)
        t0 = time.time()
        LR = LM.limit_rig(cs, As, c0)
        mu_lo = float(LR['mu'].lower())
        for e in LR['events']:
            ev_rows.append([n, e['k'], e['kernel'], e['feature'], float(e['tau'].mid()), float(e['tau'].rad()), float(e['p_e'].mid()), float(e['p_e'].rad()), float(e['d0'].mid()), e['gap'], e['rad'],
                            float(e['sign_ok']), float(LR['ok']), mu_lo, float(LR['checks']['wall']), float(LR['checks']['strict_entrant'] and LR['checks']['strict_feature'])])
        ehi = EHI if not quick else -9.0
        top0 = old[n][7]
        best = None
        if LR['ok']:
            for dtop in DTOPS:
                top = top0 + dtop
                sg, complete, cur, calls = S.cover_n(prep, LR, top, ehi, budget=budget)
                if complete:
                    best = (top, len(sg), calls, sg, cur, 1.0)
                    break
                if best is None or cur < best[4]:
                    best = (top, len(sg), calls, sg, cur, 0.0)
                log('  attempt top', top, 'stopped at', cur, 'calls', calls, flush=True)
        if best is None:
            best = (top0, 0, 0, [], top0, 0.0)
        top, nb, calls, sg, cur, comp = best
        wmax = max([w_[1] - w_[0] for w_ in sg]) if sg else 0.0
        # the tail ball (0, eps_meet], eps_meet = the lower end of the chain's last ball
        eps_meet = R3.fdown(R3.pow10(ehi))
        hhi = TL.h_of_eps(eps_meet)
        o, pi = S.best_cert_n(prep, eps_meet)
        tail = [np.nan] * 10
        Dtail, Etail = np.zeros(len(c0)), np.zeros(len(c0))
        if LR['ok'] and pi >= 0:
            ok, inf, Dv = TL.solve_D_rig(prep, LR, T.PARS[pi], hhi, o['Etot'])
            join = bool(R3.lo(TL.eps_of_h(hhi)) >= arb(eps_meet))
            Dtail = np.array(Dv, float)
            if 'Etot' in inf:
                for j, v in inf['Etot'].items():
                    Etail[j] = v
            tail = [float(ok), float(pi), float(hhi.upper()), float(R3.up(TL.eps_of_h(hhi))), inf.get('worst_ratio', np.nan), inf.get('dcmax', np.nan), inf.get('omega', np.nan),
                    inf.get('drift_max', np.nan), float(join), float(all(inf.get('cond', {}).values())) if 'cond' in inf else 0.0]
        trows.append([n, K, sep, R, seed, old[n][7], top, nb, calls, wmax, cur, comp, ehi, eps_meet] + tail + [time.time() - t0])
        segs_all.append(np.array([[n, w_[0], w_[1]] for w_ in sg], float).reshape(-1, 3))
        extra['E_D%d' % n] = Dtail
        extra['E_E%d' % n] = Etail
        if sg:
            _, pi2 = S.best_cert_n(prep, 10.0 ** top)
            hs = C3.hyp_slacks(prep, LR, T.PARS[pi2])
            if hs is not None:
                hy_rows.append([n] + hs)
        log('E target', n, (K, sep, R, seed), 'top', top, 'lowest', cur, 'complete', comp, 'balls', nb, 'tail', tail, flush=True)
    out = dict(E_rows=np.array(trows, float), E_segs=np.concatenate(segs_all) if segs_all else np.zeros((0, 3)), EV_rows=np.array(ev_rows, float).reshape(-1, len(C3.EV_COLS)),
               HY_rows=np.array(hy_rows, float).reshape(-1, len(C3.HY_COLS)), E_budget=np.array([budget]), E_ehi=np.array([EHI if not quick else -9.0]))
    out.update(extra)
    return out


def run_N(raw, quick=False, log=print):
    """the bounds against the integrated network (n = 1, K <= 3 set) at the corrected eps_0 of each target and one decade below: actual maximal centre displacement of the dormant kernels and actual
    dual deviation, against dcmax and E of the legacy and of the corrected accounting (same parameters)"""
    warnings.filterwarnings('ignore')
    ex, tg = raw['FC_ex'], raw['FC_targets']
    nn = raw['FC_okc_new'] > 0.5
    rows = []
    for i, t in enumerate(tg):
        e0 = thr_first(nn[i], ex)
        if not np.isfinite(e0):
            continue
        prep = T.prepare(*T.make_target(*[int(t[0]), float(t[1]), float(t[2]), int(t[3])]))
        ent = [e['kernel'] for e in prep['L']['events']]
        non = [j for j in range(len(prep['c0'])) if j not in ent]
        for de in (0.0, -1.0):
            e = e0 + de
            on, pn = S.best_cert_n(prep, 10.0 ** e)
            ol, pl = T.best_cert(prep, 10.0 ** e)
            if pn < 0:
                continue
            par = T.PARS[pn]
            sol, U = T.simulate(prep, 10.0 ** e)
            J = len(prep['c0'])
            act = float(np.abs(sol.y[J:][non] - prep['c0'][non][:, None]).max())
            dd = T.dual_dev(prep, 10.0 ** e, sol, U)
            oL = T.cert(prep, 10.0 ** e, par)
            rows.append([i, e, pn, act, oL['dcmax'], on['dcmax'], dd, oL['Dnew'], on['Dnew']])
            log('N', rows[-1], flush=True)
            if quick:
                break
        if quick:
            break
    return dict(N_rows=np.array(rows, float).reshape(-1, 9))


def run_G(raw, targets, budget=BUDGET, log=print):
    """the full chain of theory3c.md, [10^-24, 10^top], with the corrected verifier (top = the chain top of Part E)"""
    warnings.filterwarnings('ignore')
    tops = {int(r[0]): r[6] for r in raw['E_rows']}
    rows, segs_all = [], []
    for n in targets:
        K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
        cs, As, c0 = T.make_target(K, sep, R, seed)
        prep = T.prepare(cs, As, c0)
        LR = LM.limit_rig(cs, As, c0)
        t0 = time.time()
        sg, complete, cur, calls = S.cover_n(prep, LR, tops[n], W.Q_BOT, budget=budget)
        wmin = min([b - a for a, b in sg]) if sg else 0.0
        rows.append([n, tops[n], len(sg), calls, wmin, cur, float(complete), time.time() - t0])
        segs_all.append(np.array([[n, a, b] for a, b in sg], float).reshape(-1, 3))
        log('G', rows[-1], flush=True)
    return dict(G_rows=np.array(rows, float).reshape(-1, 8), G_segs=np.concatenate(segs_all) if segs_all else np.zeros((0, 3)))


def run_T(targets, quick=False, log=print):
    """how far up the single tail ball reaches on the enclosed limit data: the largest eps_hi on a quarter-decade grid (from the chain top down) at which `solve_D_rig` accepts; A's stored value for comparison"""
    warnings.filterwarnings('ignore')
    z3c, z3a = _load('theory3c_raw.npz'), _load('theory3a_raw.npz')
    old3c = {int(r[0]): r for r in z3c['T_rows']}
    old3a = {int(r[0]): r for r in z3a['T_rows']}
    rows = []
    for n in targets:
        K, sep, R, seed = T.C_TARGETS[W.Q_IDX[n]]
        cs, As, c0 = T.make_target(K, sep, R, seed)
        prep = T.prepare(cs, As, c0)
        LR = LM.limit_rig(cs, As, c0)
        top = old3c[n][7]
        row = [n, top, old3a[n][7], np.nan, np.nan, np.nan, np.nan, 0.0]
        for e in np.arange(top, -12.01, -0.25):
            o, pi = S.best_cert_n(prep, 10.0 ** e)
            if pi < 0:
                continue
            row[7] += 1
            ok, inf, Dv = TL.solve_D_rig(prep, LR, T.PARS[pi], TL.h_of_eps(10.0 ** e), o['Etot'])
            if ok:
                row[3:7] = [e, inf['worst_ratio'], float(pi), max(Dv)]
                break
        rows.append(row)
        log('T', row, flush=True)
    return dict(T_rows=np.array(rows, float).reshape(-1, 8))


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(r) + ' |')
    return '\n'.join(out)


def fe(x):
    return 'none' if not np.isfinite(x) else '1e%+.1f' % x


SETNAME = {'C': 'K <= 3, n = 1 (theory_t1.md)', 'K': 'K = 4..6, n = 1 (theory_t1b.md, Part K)', 'P': 'planar (theory_t1b.md, Part P)'}


def report_F(raw, L):
    P = L.append
    P('## F. The floating-point certificate with the corrected system on all grid pairs\n')
    P('`cert` / `cert2d` of `theory3n_system` (repairs (f), (g), (6)) on EVERY (target, eps) pair of the grids of theory_t1.md and theory_t1b.md, not only the informative ones; "old" = the stored '
      'status of those reports (legacy accounting).  eps_0 = the largest grid value from which the certificate holds at every smaller one; eps_max = the same for the measured network (stored); '
      '*violations* = grid pairs with the corrected certificate true and the network not forming exactly the limit\'s units.\n')
    tot = dict(old=0, new=0, scanned=0, flip_off=0, flip_on=0, viol=0)
    for key in ('C', 'K', 'P'):
        if 'F%s_ex' % key not in raw:
            continue
        ex, tg = raw['F%s_ex' % key], raw['F%s_targets' % key]
        oo, nn, om = raw['F%s_okc_old' % key] > 0.5, raw['F%s_okc_new' % key] > 0.5, raw['F%s_okm' % key] > 0.5
        P('### %s\n' % SETNAME[key])
        rows = []
        for i, t in enumerate(tg):
            e_old, e_new, e_max = thr_first(oo[i], ex), thr_first(nn[i], ex), thr_first(om[i], ex)
            dco, dcs = raw['F%s_dc_old' % key][i], raw['F%s_dc_same' % key][i]
            mk = np.isfinite(dco) & np.isfinite(dcs) & (dco > 0)
            fac = (dcs[mk] / dco[mk]) if mk.any() else np.array([np.nan])
            i0 = int(np.flatnonzero(ex == e_old)[0]) if np.isfinite(e_old) else -1
            fac0 = ('%.2f' % (dcs[i0] / dco[i0])) if (i0 >= 0 and np.isfinite(dcs[i0]) and dco[i0] > 0) else '-'
            Dn = raw['F%s_D_new' % key][i]
            D0 = ('%.1f' % Dn[int(np.flatnonzero(ex == e_new)[0])]) if np.isfinite(e_new) else '-'
            rows.append(['%d' % t[0], '%g' % t[1] if t[0] > 1 else '-', '%g' % t[2], '%d' % t[3], fe(e_old), fe(e_new), fe(e_max), ('%.1f' % (e_max - e_new)) if np.isfinite(e_new) else '-',
                         '%d / %d' % (oo[i].sum(), nn[i].sum()), '%d' % int((nn[i] & ~om[i]).sum()), fac0, '%.2f-%.2f' % (np.nanmin(fac), np.nanmax(fac)) if mk.any() else '-', D0])
        P(table(['K', 'sep', 'R', 'seed', 'eps_0 old', 'eps_0 corrected', 'eps_max (measured)', 'log10 ratio eps_max/eps_0', 'informative pairs old / corrected', 'violations',
                 'dcmax factor at old eps_0', 'dcmax factor over all old-informative pairs (same parameters)', 'D at corrected eps_0'], rows))
        n_old, n_new, n_all = int(oo.sum()), int(nn.sum()), int(oo.size)
        f_off, f_on = int((oo & ~nn).sum()), int((~oo & nn).sum())
        viol = int((nn & ~om).sum())
        Eold, Esame = raw['F%s_E_old' % key], raw['F%s_E_same' % key]
        mk = np.isfinite(Eold) & np.isfinite(Esame)
        dco, dcs = raw['F%s_dc_old' % key], raw['F%s_dc_same' % key]
        m2 = np.isfinite(dco) & np.isfinite(dcs) & (dco > 0)
        P('\nInformative pairs (certificate true): old %d, corrected %d, of %d scanned; true -> false: %d; false -> true: %d; violations among the corrected informative pairs: **%d**.  '
          'Over the old-informative pairs at unchanged parameters: dcmax factor (corrected / old) median %.2f, range %.2f-%.2f; max E relative change %.2e.\n' % (
              n_old, n_new, n_all, f_off, f_on, viol, np.median(dcs[m2] / dco[m2]) if m2.any() else np.nan, np.min(dcs[m2] / dco[m2]) if m2.any() else np.nan,
              np.max(dcs[m2] / dco[m2]) if m2.any() else np.nan, np.max(Esame[mk] / Eold[mk] - 1) if mk.any() else np.nan))
        for nm, v in zip(('old', 'new', 'scanned', 'flip_off', 'flip_on', 'viol'), (n_old, n_new, n_all, f_off, f_on, viol)):
            tot[nm] += v
    P('**Total over the three sets: informative pairs old %d, corrected %d, of %d scanned; true -> false %d, false -> true %d; violations %d.**\n' % (tot['old'], tot['new'], tot['scanned'], tot['flip_off'], tot['flip_on'], tot['viol']))
    net = raw['F_net']
    P('Network re-run on every pair whose status changed (the network is the same simulator as in the earlier reports):\n')
    if len(net):
        rows = []
        for r in net:
            key = 'CKP'[int(r[0])]
            t = raw['F%s_targets' % key][int(r[1])]
            rows.append([key, '(%d, %g, %g, %d)' % tuple(t), '%+.1f' % r[2], 'yes' if r[3] > 0 else 'no', 'yes' if r[4] > 0 else 'no', 'yes' if r[5] > 0 else 'no', 'yes' if r[6] > 0 else 'no'])
        P(table(['set', 'target (K, sep, R, seed)', 'log10 eps', 'old certificate', 'corrected certificate', 'network follows (stored)', 'network follows (re-run)'], rows))
        P('\nNo changed pair is a counterexample: %s.\n' % ('YES' if all((r[4] < 0.5) or (r[6] > 0.5) for r in net) else 'NO - SEE ROWS'))
    else:
        P('(none)\n')


def report_B(raw, L):
    P = L.append
    P('## B. The 364 stored eps-balls of theory3c, legacy against corrected ball verifier\n')
    P('Same limit enclosure, same ball, same candidate bounds D_j (corrected float certificate x 1.1 ... 4); "pass" = the verifier accepts some factor; "same factor" = the corrected verifier accepts the '
      'factor the legacy verifier accepted first.  dcmax and drift ratios (corrected / legacy) are taken at the factor at which the corrected verifier first accepts, i.e. for identical D_j; '
      'worst E/D = the largest ratio over the balls at the first accepted factor.\n')
    rows = []
    for r in raw['B_rows']:
        rows.append(['%d' % r[0], '%d' % r[1], '%d' % r[2], '%d' % r[3], '%d' % r[4], '%d' % r[5], '%.2f / %.2f / %.2f' % (r[6], r[7], r[8]), '%.2f' % r[9], '%.4f / %.4f' % (r[10], r[11])])
    P(table(['target', 'stored balls', 'pass legacy', 'pass corrected', 'corrected passes with the legacy factor', 'corrected fails', 'dcmax ratio min / median / max', 'max drift ratio', 'worst E/D legacy / corrected'], rows))
    P('\nTotals: %d stored balls, %d accepted by the legacy verifier, %d by the corrected one.  A ball on which the corrected verifier fails is re-covered by narrower balls in Parts E and G (the chains below are fresh covers, not these balls).\n' % (
        raw['B_rows'][:, 1].sum(), raw['B_rows'][:, 2].sum(), raw['B_rows'][:, 3].sum()))


def report_E(raw, L):
    P = L.append
    Tr, Ev, Hy = raw['E_rows'], raw['EV_rows'], raw['HY_rows']
    ehi = raw['E_ehi'][0]
    P('## E. Enclosed verification: chain of eps-balls plus one tail ball\n')
    P('Limit object enclosed as in round theory3 C (Krawczyk enclosure of every block critical point, strict entry decisions); the hypotheses r > d_0, sgn p = sgn A, every stage\'s critical point '
      'in D-hat (and unique there), theta < v_# are checked on the enclosures.  Chain: `theory3n_system.verify_rig_n` on balls of eps from the top down to 10^%g (outward-rounded ends, consecutive balls '
      'overlap).  Tail: `theory3n_tail.verify_tail_rig` on the single ball h = 1/kappa in [0, h_hi], h_hi an upper bound of sqrt(eps/(2 v_s - eps)) at the double eps_meet = the lower end of the '
      'chain\'s last ball: it covers EVERY real eps in (0, eps_of_h(h_hi)], and the join test is eps_of_h(h_hi) >= eps_meet (certain).\n' % ehi)
    P('### Limit object (ball arithmetic, 220 bits)\n')
    rows = []
    for e in Ev:
        rows.append(['%d' % e[0], '%d' % e[1], '%d' % e[2], '%d' % e[3], '%.15f' % e[4], '%.1e' % e[5], '%.12f' % e[6], '%.3f' % e[8], '%.4f' % e[9], '%.0e' % e[10], 'yes' if e[11] > 0 else 'NO',
                     'yes' if (e[12] > 0 and e[14] > 0 and e[15] > 0) else 'NO'])
    P(table(['target', 'k', 'entrant kernel', 'feature', 'tau_k', 'tau_k rad', 'p_k', 'd_0 (sigma)', 'relative gap to the 2nd hitting time (lower bound)', 'block box radius', 'sgn p = sgn A', 'strict & wall < 1'], rows))
    P('\n### Hypotheses on the enclosures at the top parameter set (slacks: positive = holds)\n')
    rows = []
    for h in Hy:
        rows.append(['%d' % h[0], '%.2f' % h[1], '%.3f' % h[2], '%.3e' % h[3], '%.3e' % h[4], '%.3e' % h[5], '%.3e' % h[6], '%.4f' % h[7], '%.4f' % h[8], '%.3f' % h[9], 'yes' if h[10] > 0 else 'NO'])
    P(table(['target', 'r', 'r - max d_0', 'd-hat - max offset (all stages)', 'min (|w| - a_min v_D^-)', 'min (a_min v_D^+ - |w|)', 'v_# - theta', 'd-hat', 'v_D^- / a_min', 'v_D^+ / a_min',
             'one critical point in all of D-hat (Krawczyk, every stage)'], rows))
    P('\n### Chain\n')
    rows = []
    for t in Tr:
        rows.append(['%d' % t[0], '%d' % t[1], '%s' % ('%g' % t[2] if t[1] > 1 else '-'), '%g' % t[3], '%d' % t[4], '%+.2f' % t[5], '%+.2f' % t[6], '%d' % t[7], '%d' % t[8], '%.3f' % t[9], '%+.2f' % t[10],
                     'complete' if t[11] > 0 else 'INCOMPLETE'])
    P(table(['target', 'K', 'sep', 'R', 'seed', 'top in theory3c.md', 'TOP (this report)', 'balls', 'verify calls', 'widest ball (decades)', 'lowest exponent', 'status'], rows))
    P('\n### Tail ball\n')
    rows = []
    for t in Tr:
        D = raw['E_D%d' % int(t[0])]
        Ee = raw['E_E%d' % int(t[0])]
        Dm = D[D > 1.0].max() if (D > 1.0).any() else float('nan')
        rows.append(['%d' % t[0], '%+.2f' % t[12], '%.17g' % t[13], '%.6e' % t[16], '%.6e' % t[17], 'yes' if t[14] > 0 else 'NO', '%d' % t[15], '%.1f' % Dm, '%.1f' % Ee.max(), '%.4f' % t[18], '%.2e' % t[19], '%.2e' % t[20], '%.3f' % t[21],
                     'yes' if t[22] > 0 else 'NO', 'yes' if t[23] > 0 else 'NO'])
    P(table(['target', 'log10 eps_hi', 'eps_meet (double)', 'h_hi', 'eps_of_h(h_hi)', 'strict E < D/1.0001 and all conditions', 'parameter set', 'max D_j', 'max E_j', 'worst E/D', 'dcmax', 'omega', 'drift (max)', 'join eps_of_h >= eps_meet', 'all conditions'], rows))
    P('\n(E_j, D_j are the bounds in eta units of the dual deviation of the dormant kernels; the tail ball has ONE set of D_j for all eps <= eps_hi.)\n')
    ncomp = int(sum(1 for t in Tr if t[11] > 0 and t[14] > 0 and t[22] > 0 and t[23] > 0))
    P('**%d of %d targets: S(eps) is verified, in ball arithmetic with Krawczyk-enclosed limit data and certain comparisons, for every real eps in (0, 10^TOP] (chain [10^%g, 10^TOP] plus the tail ball (0, 10^%g]).**\n' % (ncomp, len(Tr), ehi, ehi))
    for t in Tr:
        P('* target %d (K=%d, sep=%g, R=%g, seed=%d): top 10^%+.2f (theory3c.md: 10^%+.2f); %s' % (t[0], t[1], t[2], t[3], t[4], t[6], t[5], 'EARNED' if (t[11] > 0 and t[14] > 0 and t[22] > 0 and t[23] > 0) else 'NOT EARNED'))
    P('')


def report_N(raw, L):
    P = L.append
    P('## N. The bounds against the integrated network (n = 1, K <= 3 set)\n')
    P('At the corrected eps_0 of each target and one decade below (same certificate parameters for both accountings): the largest centre displacement max_t |c_l(t) - c_l^0| of the non-entrant '
      'kernels in the integrated continuum network, and the largest dual deviation max |eta_l - ubar_l| (eta units), against the bounds dcmax and E (max over non-entrants) of the legacy and of the '
      'corrected accounting.  A bound must lie above the measured value; the ratio shows how crude it is.\n')
    rows = []
    for r in raw['N_rows']:
        t = raw['FC_targets'][int(r[0])]
        rows.append(['(%d, %g, %g, %d)' % tuple(t), '%+.1f' % r[1], '%d' % r[2], '%.3e' % r[3], '%.3e' % r[4], '%.3e' % r[5], 'yes' if r[5] >= r[3] else 'NO', '%.3f' % r[6], '%.2f' % r[7], '%.2f' % r[8],
                     'yes' if r[8] >= r[6] else 'NO'])
    P(table(['target (K, sep, R, seed)', 'log10 eps', 'parameter set', 'measured max displacement', 'dcmax legacy', 'dcmax corrected', 'corrected >= measured', 'measured dual deviation', 'E legacy', 'E corrected', 'corrected >= measured'], rows))
    P('\nRows where the corrected dcmax is below the measured displacement: %d; where the legacy dcmax is below it: %d (a trajectory can only fail to contradict a bound: the legacy form is not a bound of '
      'Lemma 5.2\'s displacement integral, but the measured displacements here do not exceed it); rows where the corrected E is below the measured dual deviation: %d.\n' % (
          int(sum(1 for r in raw['N_rows'] if r[5] < r[3])), int(sum(1 for r in raw['N_rows'] if r[4] < r[3])), int(sum(1 for r in raw['N_rows'] if r[8] < r[6]))))


def report_G(raw, L):
    P = L.append
    P('## G. The full chain of theory3c.md with the corrected verifier\n')
    P('The same cover as `theory3c.md` (enclosed limit object, balls from the chain top down to 10^%g), re-run with `verify_rig_n`: the tops are those of Part E.\n' % W.Q_BOT)
    rows = []
    for r in raw['G_rows']:
        rows.append(['%d' % r[0], '%+.2f' % r[1], '%d' % r[2], '%d' % r[3], '%.5f' % r[4], '%+.2f' % r[5], 'complete' if r[6] > 0 else 'INCOMPLETE'])
    P(table(['target', 'top', 'balls', 'verify calls', 'narrowest ball (decades)', 'lowest exponent', 'status'], rows))
    P('\n%d of %d full chains complete down to 10^%g.\n' % (int(sum(1 for r in raw['G_rows'] if r[6] > 0)), len(raw['G_rows']), W.Q_BOT))


def report_T(raw, L):
    P = L.append
    ehi = raw['E_ehi'][0] if 'E_ehi' in raw else EHI
    P('## T. How far up one tail ball reaches on the enclosed limit data\n')
    P('The largest eps_hi on a quarter-decade grid from the chain top downwards at which `solve_D_rig` proves S(eps) for EVERY real eps in (0, eps_hi] with one set of bounds D_j (enclosed limit data, '
      'corrected accounting).  "theory3a.md" = the value of round theory3 A (floating-point limit data, accounting (f) repaired, (g) missing).  The merged statement does not need this reach: the chain '
      'covers [eps_hi, top] and the tail ball (0, eps_hi] with eps_hi = 10^%g; the table says how much of the chain could be replaced by the tail.\n' % ehi)
    rows = []
    for r in raw['T_rows']:
        rows.append(['%d' % r[0], '%+.2f' % r[1], '%+.2f' % r[2], ('%+.2f' % r[3]) if np.isfinite(r[3]) else 'none', ('%.4f' % r[4]) if np.isfinite(r[4]) else '-', ('%d' % r[5]) if np.isfinite(r[5]) else '-',
                     ('%.1f' % r[6]) if np.isfinite(r[6]) else '-', '%d' % r[7]])
    P(table(['target', 'chain top', 'eps_tail of theory3a.md', 'eps_tail here (enclosed data, corrected accounting)', 'worst E/D', 'parameter set', 'max D_j', 'grid values tried'], rows))
    P('')


def report_W(raw, L):
    P = L.append
    P('## W. What changed relative to the earlier reports\n')
    rows = []
    for key in ('C', 'K', 'P'):
        if 'F%s_ex' % key not in raw:
            continue
        ex, tg = raw['F%s_ex' % key], raw['F%s_targets' % key]
        oo, nn, om = raw['F%s_okc_old' % key] > 0.5, raw['F%s_okc_new' % key] > 0.5, raw['F%s_okm' % key] > 0.5
        rows.append([SETNAME[key], '%d' % oo.sum(), '%d' % nn.sum(), '%d' % oo.size, '%d' % (oo & ~nn).sum(), '%d' % (~oo & nn).sum(), '%d' % (nn & ~om).sum()])
        for i, t in enumerate(tg):
            e_old, e_new, e_max = thr_first(oo[i], ex), thr_first(nn[i], ex), thr_first(om[i], ex)
            if not (e_old == e_new or (not np.isfinite(e_old) and not np.isfinite(e_new))):
                P('* eps_0 of %s target (K=%d, sep=%g, R=%g, seed=%d): %s -> %s; log10(eps_max/eps_0) %.1f -> %.1f.' % (key, t[0], t[1], t[2], t[3], fe(e_old), fe(e_new), e_max - e_old, e_max - e_new))
    P('')
    P(table(['set', 'informative pairs old', 'informative pairs corrected', 'pairs scanned', 'true -> false', 'false -> true', 'violations'], rows + [[
        'total', '%d' % sum(int(r[1]) for r in rows), '%d' % sum(int(r[2]) for r in rows), '%d' % sum(int(r[3]) for r in rows), '%d' % sum(int(r[4]) for r in rows), '%d' % sum(int(r[5]) for r in rows),
        '%d' % sum(int(r[6]) for r in rows)]]))
    P('')
    if 'B_rows' in raw:
        B = raw['B_rows']
        P('* the stored balls of theory3c: %d of %d accepted by the legacy verifier, %d of %d by the corrected one (%d with the legacy factor).' % (B[:, 2].sum(), B[:, 1].sum(), B[:, 3].sum(), B[:, 1].sum(), B[:, 4].sum()))
    if 'E_rows' in raw:
        Tr = raw['E_rows']
        P('* tops of theory3c.md: unchanged for %d of %d targets (%s).' % (sum(1 for t in Tr if abs(t[5] - t[6]) < 1e-9), len(Tr), ', '.join('%d: %+.2f -> %+.2f' % (t[0], t[5], t[6]) for t in Tr)))
        P('* the statement of theory3c.md "S(eps) verified for every real eps in [10^-24, 10^top]" becomes "for every real eps in (0, 10^top]": the interval (0, 10^%g] is now one tail ball on the enclosed data.' % raw['E_ehi'][0])
    if 'T_rows' in raw:
        Tt = raw['T_rows']
        P('* tail of theory3a.md (floating-point limit data, accounting (g) missing): eps_tail %s; here, enclosed data, corrected accounting: %s.' % (
            ', '.join('%d: %+.2f' % (t[0], t[2]) for t in Tt), ', '.join('%d: %s' % (t[0], ('%+.2f' % t[3]) if np.isfinite(t[3]) else 'none') for t in Tt)))
    P('')


GAPS = """## Status and gaps

* COMPUTER-ASSISTED (this report).  Hypotheses: continuum Gram of the Gaussian kernel of width sigma = 0.12 (lengths in widths), n = 1, K <= 3 features at the stored targets, equal masses, noiseless data, the
  repository's law (m, G, gamma, a_max) = (12, 1, 0.1, 1), the stored lattice dictionary, the P0 unit from alpha = 0, horizon [0, lambda (1.05 tau_K + 1)].  Statement: for each of the 8 targets, the finite
  system S(eps) of Theorem S-T1 (with the repairs (f), (g), (6)) holds, in ball arithmetic (python-flint) with Krawczyk-enclosed limit data (existence and uniqueness of each stage block's critical point in the
  relaxation set, strict entry decisions, hypotheses on the enclosures, certain comparisons), for EVERY real eps in (0, 10^top]: a chain of eps-balls on [10^-12, 10^top] and one tail ball on (0, 10^-12].
* NOT verified: (i) that the enclosed critical points are the GLOBAL minimisers of Theorem K(a) (what T1 uses is the critical point in D-hat, which is unique there); (ii) the hand proof that S(eps) implies the
  conclusions of Theorem T1 (the referee re-derived only the displacement term of Lemma 5.2 and Lemmas E and M; the other steps of sections 2-6 of THEORY_T1.md were not independently re-derived, and the two
  omissions (f), (g) show that transcription gaps are possible there); (iii) any monotonicity of S in eps (not needed: the cover is explicit); (iv) the soundness of arb itself; (v) K >= 4 targets and the
  planar targets (floating point on a half-decade grid only: no enclosure, no tail); (vi) sampled (non-continuum) Gram, other kernels (the Gram tail fix of `theory3n_kernels.py` is tested but the
  non-Gaussian certificates are not re-run here), noise, other optimisers.
* NOT regenerated: the ingredient tables (Parts A, B, D, E, F, G, H of theory_t1.md; the refinement below eps_0 and the ingredient diagnostics of theory_t1b.md; Part U of theory_t1b.md) were computed with the
  legacy accounting; only the grids of Part C of theory_t1.md, Parts K and P of theory_t1b.md (floating-point certificate on the grid) and the interval covers are redone here.
* The candidate bounds D_j and the parameter sets come from the floating-point certificate; they only propose, the ball verifier decides (every accepted ball satisfies E_j < D_j/1.0001 strictly).
"""


def report(raw):
    L = ['# Theorem T1: corrected and merged verifications (generated by experiments/theory3n_run.py)\n',
         'Round theory3, direction N.  This report supersedes the verifications of `theory_t1.md`, `theory_t1b.md`, `theory3a.md` and `theory3c.md`, which stay as history.  The referee of round T3 '
         '(review/REVIEW_roundT3AC.md) found that the displacement accounting of the certificate S(eps) was wrong in every verifier: (f) the term H^(st(t)) int_0^t B must be the stage-resolved '
         'integral of H^(st(t\')) B; (g) the pre-layer coefficient psi_sup(X_li) pre of every entrant whose window has not yet started enters d_c R_l and was dropped by the window gating.  Also built '
         'in: (6) the old units\' coefficient in a later window is the one of the CURRENT stage (stage maximum, not the final block\'s), kappa is selected by the smaller LOWER end, and the float gate of '
         '`inv_of` is replaced by a certain comparison.  The correction of the Gram tails below the table start 0.4 (E-2 of review/REVIEW_roundT3E.md) is provided and tested in `theory3n_kernels.py`; '
         'not used by any result of this report.\n']
    if any(k.startswith('F') and k.endswith('_ex') for k in raw):
        report_F(raw, L)
    if 'B_rows' in raw:
        report_B(raw, L)
    if 'E_rows' in raw:
        report_E(raw, L)
    if 'T_rows' in raw:
        report_T(raw, L)
    if 'G_rows' in raw:
        report_G(raw, L)
    if 'N_rows' in raw:
        report_N(raw, L)
    report_W(raw, L)
    L.append(GAPS)
    return '\n'.join(L)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--parts', default='F,B,E,T')
    ap.add_argument('--targets', default=','.join(str(v) for v in VERIFIED))
    ap.add_argument('--budget', type=float, default=BUDGET)
    a = ap.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    raw = {}
    if os.path.exists(RAW):
        z = np.load(RAW, allow_pickle=False)
        raw = {k: z[k] for k in z.files}
    if a.run:
        tg = [int(x) for x in a.targets.split(',')]
        parts = a.parts.split(',')
        if 'F' in parts:
            raw.update(run_F(a.quick))
            np.savez(RAW, **raw)
        if 'B' in parts:
            raw.update(run_B(tg, a.quick))
            np.savez(RAW, **raw)
        if 'E' in parts:
            raw.update(run_E(tg, a.budget, a.quick))
            np.savez(RAW, **raw)
        if 'T' in parts:
            raw.update(run_T(tg, a.quick))
            np.savez(RAW, **raw)
        if 'N' in parts:
            raw.update(run_N(raw, a.quick))
            np.savez(RAW, **raw)
        if 'G' in parts:
            raw.update(run_G(raw, tg, a.budget))
            np.savez(RAW, **raw)
    with open(MD, 'w') as fh:
        fh.write(report(raw))
    print('report written:', MD)


if __name__ == '__main__':
    main()
