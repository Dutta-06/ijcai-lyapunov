"""
POST HOC merge analysis of the KAN spike: how much of "72% against 28%" is unmerged duplicate units of the L1 + prune control?
==============================================================================================================================

The referee's report on the four ML spikes (review/REVIEW_roundML.md, section 1, findings 1.1 - 1.4) found that the control's exact-edge share (0.283) is the
empty edges (K = 0), that its surplus units sit on the non-empty edges and are mostly coincident duplicates, and that merging units whose centres lie within a
small radius closes most or all of the gap.  This script produces those numbers from the stored raw of the KAN spike (nothing is re-trained; no torch).

READS only results/ml/kan_raw.npz (the records of experiments/exp_kan_spike.py, stage test).  A record of the 1-layer synthetic runs stores the centres of its
active units (`centres_0.02`, one list per edge, activity threshold tau = 0.02), the unit counts (`counts_0.02`) and the truth (`truth`: number of bumps per edge;
`truth_c`: their centres).  Every ours / control row is selected by the spike module's OWN rules (ours: the single record of the run; control: one-standard-error
rule on the validation MSE (`one_se_pick`); oracle control: the penalty with the highest exact share, ties by fewer active parameters), so that the delta = 0 rows
reproduce results/ml/kan.md exactly.  No new raw file is needed: the report regenerates directly from kan_raw.npz.

POST HOC.  The analysis was designed after the test set was scored and after reading the referee's report; the radius delta is a free choice, so the sweep is
shown over the whole grid and no single delta is selected.  Only the stored centres are used (no widths, no weights, no amplitudes): a merge here means "count
single-linkage clusters of the active centres of an edge", not "refit the merged unit".

    python experiments/exp_kan_merge.py --from-raw [RAW.npz] [--md OUT.md]     (numpy only; default KAN_OUT/kan_raw.npz -> KAN_OUT/kan_merge.md; default KAN_OUT: results/ml)
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exp_kan_spike as K  # noqa: E402

DELTAS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)    # merge radii of the sweep
GAP_THRESH = (0.1, 0.3)                    # adjacent-gap thresholds of the duplicate statistics
MATCH_TOLS = (0.2, 0.5)                    # tolerances of the centre precision / recall
N_BOOT, BOOT_SEED = 20000, 0               # cluster bootstrap over targets
CKEY = 'centres_%g' % K.TAU
NKEY = 'counts_%g' % K.TAU
METHODS = (('ours', 'ours (factorised, small start)'), ('ctrl', 'control: L1 + prune, penalty by validation'), ('ctrl_or', 'control: penalty by the truth (oracle, unmerged)'))


def fmt(x, n=3):
    return K.fmt(x, n)


# ------------------------------------------------------------------------------------------------------------------------------------
# selection of the rows (the spike module's own rules)
# ------------------------------------------------------------------------------------------------------------------------------------
def select_runs(recs, tau=K.TAU):
    """One dict per complete synthetic run, in the order of the spike module's synth_summary: key, truth, truth_c, the selected records ours / ctrl / ctrl_or and
    the list of all control records `ctrl_all`."""
    runs = {}
    for r in recs:
        if r['kind'] == 'synth':
            runs.setdefault(K.run_key(r), []).append(r)
    out = []
    for k in sorted(runs):
        rs = runs[k]
        o = [r for r in rs if r['method'] == 'ours']
        ctrl = [r for r in rs if r['method'] == 'ctrl']
        if not o or not ctrl or not any(r['method'] == 'fixed' for r in rs):
            continue                                            # an incomplete run is not scored
        o = o[0]
        cands = [dict(val=r['val_mse'], P=r['P_%g' % tau], r=r) for r in ctrl]
        cv = K.one_se_pick(cands, ctrl[0]['n_val'])['r']
        co = max(ctrl, key=lambda r: (K.exact_share(r, tau), -r['P_%g' % tau]))
        out.append(dict(key=k, truth=list(o['truth']), truth_c=o['truth_c'], ours=o, ctrl=cv, ctrl_or=co, ctrl_all=ctrl))
    return out


# ------------------------------------------------------------------------------------------------------------------------------------
# merging
# ------------------------------------------------------------------------------------------------------------------------------------
def sorted_centres(cs):
    return sorted(float(c) for c in cs)


def gaps(cs):
    """Adjacent gaps of the sorted centres of one edge (rounded to 6 decimals: the stored centres have 4)."""
    s = sorted_centres(cs)
    return [round(b - a, 6) for a, b in zip(s[:-1], s[1:])]


def merged_count(cs, delta):
    """Number of single-linkage clusters of the centres of one edge at radius delta: adjacent sorted centres closer than delta (gap < delta) are in one cluster.
    delta = 0 merges nothing."""
    n = len(cs)
    if n == 0:
        return 0
    return 1 + sum(1 for g in gaps(cs) if g >= delta)


def counts_matrix(runs, which, delta):
    """(runs, edges) matrix of merged unit counts of the selected record `which` ('ours', 'ctrl', 'ctrl_or')."""
    return np.array([[merged_count(cs, delta) for cs in r[which][CKEY]] for r in runs])


def counts_matrix_rec(recs, delta):
    return np.array([merged_count(cs, delta) for cs in recs[CKEY]])


def truth_matrix(runs):
    return np.array([r['truth'] for r in runs])


def check_consistency(runs):
    """The stored count of every record equals the number of stored centres of the edge (so that delta = 0 is the published count)."""
    n = 0
    for r in runs:
        for rec in [r['ours']] + r['ctrl_all']:
            assert [len(c) for c in rec[CKEY]] == list(rec[NKEY]), (r['key'], rec['conf'])
            n += 1
    return n


def oracle_merged(runs, delta):
    """Control with the penalty chosen per run by the truth AFTER merging at `delta` (the highest exact share of the merged counts, ties by fewer active
    parameters, then the first): (runs, edges) matrix of merged counts."""
    out = []
    for r in runs:
        t = np.array(r['truth'])
        best = max(r['ctrl_all'], key=lambda rec: (float(np.mean(counts_matrix_rec(rec, delta) == t)), -rec['P_%g' % K.TAU]))
        out.append(counts_matrix_rec(best, delta))
    return np.array(out)


def gain_ci(runs, share_a, share_b):
    """Mean paired gain a - b of per-run exact shares, cluster bootstrap over the targets (each target's gain averaged over its noise levels and seeds)."""
    tids = np.array([r['key'][0] for r in runs])
    d = np.asarray(share_a, float) - np.asarray(share_b, float)
    by = [float(np.mean(d[tids == t])) for t in sorted(set(tids.tolist()))]
    return K.boot_ci(by, n_boot=N_BOOT, seed=BOOT_SEED)


# ------------------------------------------------------------------------------------------------------------------------------------
# statistics
# ------------------------------------------------------------------------------------------------------------------------------------
def strata(T, C):
    """Exact share, number of edges, surplus (over-counted) and missing (under-counted) units per run on the empty (K = 0) and non-empty (K > 0) edges, and per
    edge of the stratum."""
    res = {}
    for nm, mask in (('empty', T == 0), ('nonempty', T > 0)):
        n = int(mask.sum())
        ex = (C == T)
        over = np.maximum(C - T, 0) * mask
        under = np.maximum(T - C, 0) * mask
        res[nm] = dict(n=n, exact_n=int(ex[mask].sum()), share=float(ex[mask].mean()) if n else float('nan'), over_run=float(over.sum(1).mean()),
                       under_run=float(under.sum(1).mean()), over_edge=float(over.sum() / n) if n else float('nan'))
    return res


def dup_stats(runs, which):
    """Adjacent gaps of the sorted active centres within an edge, over all runs and edges: (number of gaps, share below each GAP_THRESH)."""
    g = [x for r in runs for cs in r[which][CKEY] for x in gaps(cs)]
    g = np.array(g)
    return len(g), [float(np.mean(g < t)) for t in GAP_THRESH]


def match_stats(runs, which, tol):
    """Mean over runs of the spike module's per-run (recall, precision) of the active centres against the true bump centres."""
    m = np.array([K.centre_match(dict(r[which], truth_c=r['truth_c']), tol) for r in runs])
    return float(m[:, 0].mean()), float(m[:, 1].mean())


def seed_sd(runs, which, share):
    cells = {}
    for r, s in zip(runs, share):
        cells.setdefault(r['key'][:2], []).append(s)
    return float(np.mean([np.std(v) for v in cells.values()]))


def compute(recs):
    """All numbers of the report."""
    runs = select_runs(recs)
    T = truth_matrix(runs)
    n_checked = check_consistency(runs)
    res = dict(runs=runs, T=T, n_checked=n_checked)
    res['sweep'] = {}
    for d in DELTAS:
        C = {nm: counts_matrix(runs, nm, d) for nm, _ in METHODS}
        C['ctrl_best'] = oracle_merged(runs, d)
        sh = {nm: (C[nm] == T).mean(1) for nm in C}
        res['sweep'][d] = dict(C=C, share=sh, all5={nm: float(np.mean((C[nm] == T).all(1))) for nm in C},
                               units={nm: float(C[nm].sum(1).mean()) for nm in C}, nonempty={nm: float(((C[nm] == T)[T > 0]).mean()) for nm in C},
                               gain={nm: gain_ci(runs, sh['ours'], sh[nm]) for nm in ('ctrl', 'ctrl_or', 'ctrl_best')})
    return res


# ------------------------------------------------------------------------------------------------------------------------------------
# report
# ------------------------------------------------------------------------------------------------------------------------------------
def first_delta(sw, pred):
    for d in DELTAS:
        if pred(sw[d]):
            return d
    return None


def report(meta, recs):
    R = compute(recs)
    runs, T, sw = R['runs'], R['T'], R['sweep']
    nr = len(runs)
    ne = T.size
    L = ['# KAN spike: POST HOC merge analysis of the "72% against 28%" exact-edge result', '',
         'Generated by `experiments/exp_kan_merge.py --from-raw` from `results/ml/kan_raw.npz` only (stage: %s, %d records, %d complete synthetic runs, %d edges). No new raw file is '
         'needed: the stored centres of every record (`centres_%g`) are enough, and the report regenerates byte-identically from kan_raw.npz.' % (meta.get('stage'), len(recs), nr, ne, K.TAU), '',
         'This analysis is **POST HOC**: it was designed after the test set was scored, in response to the referee report (`review/REVIEW_roundML.md`, section 1, findings 1.1 - 1.4). '
         'It uses only the stored centres of the active units (no widths, no weights, no amplitudes), so "merging" means counting the single-linkage clusters of the active '
         'centres of an edge (adjacent sorted centres closer than delta are one cluster; delta = 0 merges nothing); nothing is refitted. The rows are the spike module\'s own '
         'selections (ours; control with the penalty chosen by the one-standard-error rule on the validation MSE; control with the penalty chosen by the truth), so the delta = 0 '
         'rows are the published ones. Consistency check: in all %d inspected records (ours and every control penalty of every run) the stored count equals the number of stored centres.' % R['n_checked'], '']
    # reproduction of the headline
    shs = {nm: float(np.mean(sw[0.0]['share'][nm])) for nm, _ in METHODS}
    rows0 = K.synth_summary(recs)
    pub = {nm: float(np.mean([r[nm]['share'] for r in rows0])) for nm in ('ours', 'ctrl', 'ctrl_or')}
    pm, plo, phi = K.paired_gain(rows0, 'ours', 'ctrl')
    L += ['## 0. Reproduction of the published headline (kan.md)', '',
          'Exact-edge share from the stored centres at delta = 0: ours %s, control %s, control (oracle) %s. The spike module\'s own summary gives ours %s, control %s, control (oracle) %s '
          '(match: %s). Its paired gain ours - control is %s [%s, %s] (10000 resamples, seed 12345).' % (
              fmt(shs['ours']), fmt(shs['ctrl']), fmt(shs['ctrl_or']), fmt(pub['ours']), fmt(pub['ctrl']), fmt(pub['ctrl_or']),
              'yes' if all(abs(shs[nm] - pub[nm]) < 1e-12 for nm in pub) else 'NO', fmt(pm), fmt(plo), fmt(phi)), '']
    # ---- 1. strata
    L += ['## 1. The exact-edge share split by empty and non-empty edges (unmerged, delta = 0)', '',
          'An edge is exact when the number of active units equals the number of true bumps K. %d of the %d edges have K = 0 (irrelevant inputs); %d have K > 0.' % (
              int((T == 0).sum()), ne, int((T > 0).sum())), '',
          '| method | share, all edges | share, K = 0 edges (n=%d) | share, K > 0 edges (n=%d) | surplus units / run on K = 0 edges | surplus units / run on K > 0 edges | surplus units per K > 0 edge | missing units / run on K > 0 edges | runs exact on all 5 edges |' % (
              int((T == 0).sum()), int((T > 0).sum())),
          '|---|---|---|---|---|---|---|---|---|']
    S0 = {}
    for nm, lab in METHODS:
        st = strata(T, sw[0.0]['C'][nm])
        S0[nm] = st
        L.append('| %s | %s | %s (%d/%d) | %s (%d/%d) | %s | %s | %s | %s | %d of %d (%s) |' % (
            lab, fmt(float(np.mean(sw[0.0]['share'][nm]))), fmt(st['empty']['share']), st['empty']['exact_n'], st['empty']['n'], fmt(st['nonempty']['share']),
            st['nonempty']['exact_n'], st['nonempty']['n'], fmt(st['empty']['over_run'], 2), fmt(st['nonempty']['over_run'], 2), fmt(st['nonempty']['over_edge'], 2),
            fmt(st['nonempty']['under_run'], 2), int(round(sw[0.0]['all5'][nm] * nr)), nr, fmt(sw[0.0]['all5'][nm])))
    ex_emp = {nm: S0[nm]['empty']['exact_n'] for nm, _ in METHODS}
    L += ['', 'Contribution of the empty edges to the all-edges exact share (exact empty edges / all edges): %s.' % '; '.join(
        '%s %s' % (lab, fmt(ex_emp[nm] / ne)) for nm, lab in METHODS), '',
          '### Exact-edge share by the number of true bumps K per edge (n = number of edges)', '',
          '| K | n edges | ours | control | control (oracle) | ours, delta=0.3 | control, delta=0.3 |', '|---|---|---|---|---|---|---|']
    for k in range(int(T.max()) + 1):
        m = T == k
        if not m.any():
            continue
        L.append('| %d | %d | %s | %s | %s | %s | %s |' % (k, int(m.sum()), *(fmt(float((sw[0.0]['C'][nm] == T)[m].mean())) for nm, _ in METHODS),
                                                           *(fmt(float((sw[0.3]['C'][nm] == T)[m].mean())) for nm in ('ours', 'ctrl'))))
    L.append('')
    # ---- 2. duplicates
    L += ['## 2. Duplicate units and the position of the active centres', '',
          'Adjacent gaps of the sorted active centres within an edge (all runs, all edges with at least two active units). Recall: share of true bump centres with an active unit within the '
          'tolerance; precision: share of active units within the tolerance of a true bump centre (mean over runs, as in kan.md).', '',
          '| method | adjacent gaps (n) | share of gaps below %g | share of gaps below %g | recall %g | precision %g | recall %g | precision %g |' % (*GAP_THRESH, MATCH_TOLS[0], MATCH_TOLS[0], MATCH_TOLS[1], MATCH_TOLS[1]),
          '|---|---|---|---|---|---|---|---|']
    D2 = {}
    for nm, lab in METHODS:
        n, sh = dup_stats(runs, nm)
        ms = [match_stats(runs, nm, t) for t in MATCH_TOLS]
        D2[nm] = (n, sh, ms)
        L.append('| %s | %d | %s | %s | %s | %s | %s | %s |' % (lab, n, fmt(sh[0]), fmt(sh[1]), fmt(ms[0][0]), fmt(ms[0][1]), fmt(ms[1][0]), fmt(ms[1][1])))
    L += ['', 'The generator keeps the true centres of one edge at least 1.0 apart (`make_target`), so two active units less than 0.3 apart that both lie within 0.2 of a true bump lie within 0.2 of the same bump (different bumps would put them at least 0.6 apart).', '']
    # ---- 3. sweep
    L += ['## 3. The merge sweep', '',
          'Active units of each edge merged by single linkage at radius delta (same delta for both methods; merged count = number of clusters). The paired gain is ours minus the control '
          '(validation penalty), per-run exact-edge share, cluster bootstrap over the %d targets (%d resamples, seed %d). The columns "K > 0" give the share on the non-empty edges.' % (
              len(set(r['key'][0] for r in runs)), N_BOOT, BOOT_SEED), '',
          '| delta | ours exact | control exact | gain ours - control [95% interval] | ours exact, K > 0 | control exact, K > 0 | ours units / run | control units / run | ours all 5 edges | control all 5 edges |',
          '|---|---|---|---|---|---|---|---|---|---|']
    for d in DELTAS:
        s = sw[d]
        g = s['gain']['ctrl']
        L.append('| %g | %s | %s | %s [%s, %s] | %s | %s | %s | %s | %s | %s |' % (
            d, fmt(float(np.mean(s['share']['ours']))), fmt(float(np.mean(s['share']['ctrl']))), fmt(g[0]), fmt(g[1]), fmt(g[2]), fmt(s['nonempty']['ours']), fmt(s['nonempty']['ctrl']),
            fmt(s['units']['ours'], 2), fmt(s['units']['ctrl'], 2), fmt(s['all5']['ours']), fmt(s['all5']['ctrl'])))
    allg = [x for r in runs for nm, _ in METHODS[:2] for cs in r[nm][CKEY] for x in gaps(cs)]
    ties = sum(1 for x in allg for d in DELTAS[1:] if abs(x - d) < 1e-9)
    L += ['', 'Ties: %d of the %d adjacent gaps of ours and the control equal a radius of the grid exactly, so "closer than delta" and "at most delta" give the same sweep.' % (ties, len(allg))]
    L += ['', 'Same sweep against the two controls with a penalty chosen by the truth. "control, oracle chosen at delta 0" is the penalty that is best by the truth at delta = 0 (the row of kan.md), '
          'merged afterwards; "best penalty per run at the same delta" chooses, in every run and at every delta, the one of the %d stored penalties whose merged counts have the highest exact-edge share '
          '(a selection by the truth that no practitioner has; it uses only the stored per-penalty rows, which exist for all %d runs).' % (len(runs[0]['ctrl_all']), nr), '',
          '| delta | control, oracle chosen at delta 0 | gain ours - it [95% interval] | control, best penalty per run at this delta | gain ours - it [95% interval] |', '|---|---|---|---|---|']
    for d in DELTAS:
        s = sw[d]
        go, gb = s['gain']['ctrl_or'], s['gain']['ctrl_best']
        L.append('| %g | %s | %s [%s, %s] | %s | %s [%s, %s] |' % (d, fmt(float(np.mean(s['share']['ctrl_or']))), fmt(go[0]), fmt(go[1]), fmt(go[2]), fmt(float(np.mean(s['share']['ctrl_best']))),
                                                                  fmt(gb[0]), fmt(gb[1]), fmt(gb[2])))
    L.append('')
    # ---- 4. reading
    d_ci = first_delta(sw, lambda s: s['gain']['ctrl'][1] <= 0.0)
    d_pt = first_delta(sw, lambda s: s['gain']['ctrl'][0] <= 0.0)
    d_s = first_delta(sw, lambda s: not (s['gain']['ctrl'][0] >= K.GO_SHARE_GAIN and s['gain']['ctrl'][1] > 0.0))
    g0 = sw[0.0]['gain']['ctrl']
    gtxt = '; '.join('%g: %s' % (d, fmt(sw[d]['gain']['ctrl'][0])) for d in DELTAS[1:])
    rm = {nm: float(np.mean([r[nm]['test_rmse'] for r in runs])) for nm, _ in METHODS}
    worse = float(np.mean([r['ours']['test_rmse'] > r['ctrl']['test_rmse'] for r in runs]))
    cap = int(sum(r['ours']['steps'] >= K.STEPS for r in runs))
    sd_o = seed_sd(runs, 'ours', sw[0.0]['share']['ours'])
    sd_c = seed_sd(runs, 'ctrl', sw[0.0]['share']['ctrl'])
    rec_o, rec_c = D2['ours'][2][0][0], D2['ctrl'][2][0][0]
    b5 = sw[0.5]['gain']['ctrl_best']
    ver = K.verdict(rows0, K.real_summary(recs, K.PRIMARY))
    L += ['## 4. Reading', '',
          '* The published gain in exact-edge share, ours minus the validation-penalised control, is %s (%s -> %s). After merging the active units of each edge, the gain at delta = %s is '
          'as follows (the point estimate at each radius): %s. So the published gain of %s shrinks to %s at delta = 0.1 and to %s at delta = 0.2.' % (
              fmt(g0[0]), fmt(float(np.mean(sw[0.0]['share']['ours']))), fmt(float(np.mean(sw[0.0]['share']['ctrl']))), ', '.join('%g' % d for d in DELTAS[1:]), gtxt, fmt(g0[0]),
              fmt(sw[0.1]['gain']['ctrl'][0]), fmt(sw[0.2]['gain']['ctrl'][0])),
          '* The 95%% interval of the gain first includes 0 at delta = %s, and the point estimate first is not positive at delta = %s. %s' % (
              '%g' % d_ci if d_ci is not None else 'none in the grid', '%g' % d_pt if d_pt is not None else 'none in the grid',
              'The plan\'s success criterion for the synthetic part (gain at least %g with the interval above 0) first fails at delta = %s.' % (K.GO_SHARE_GAIN, '%g' % d_s if d_s is not None else 'none')),
          '* The published "%s against %s" is largely the empty edges and unmerged duplicates: the control is exact on %d of the %d K = 0 edges (as is ours) and on %d of the %d non-empty edges '
          '(share %s), against ours on %d (share %s); the control puts %s surplus units per run on K > 0 edges (%s per edge) and %s on K = 0 edges; %s of its adjacent active-unit gaps are below %g and %s '
          'below %g (ours %s and %s).' % (
              fmt(float(np.mean(sw[0.0]['share']['ours'])), 2), fmt(float(np.mean(sw[0.0]['share']['ctrl'])), 2), S0['ctrl']['empty']['exact_n'], S0['ctrl']['empty']['n'], S0['ctrl']['nonempty']['exact_n'],
              S0['ctrl']['nonempty']['n'], fmt(S0['ctrl']['nonempty']['share']), S0['ours']['nonempty']['exact_n'], fmt(S0['ours']['nonempty']['share']), fmt(S0['ctrl']['nonempty']['over_run'], 2),
              fmt(S0['ctrl']['nonempty']['over_edge'], 2), fmt(S0['ctrl']['empty']['over_run'], 2), fmt(D2['ctrl'][1][0]), GAP_THRESH[0], fmt(D2['ctrl'][1][1]), GAP_THRESH[1],
              fmt(D2['ours'][1][0]), fmt(D2['ours'][1][1])),
          '* With the control\'s penalty chosen per run by the truth after merging at delta = 0.5 the control reaches %s (ours %s); the gain ours minus it is %s [%s, %s]. This selection uses the truth '
          'and the radius is a free choice, so it is an upper bound on what a merge step could give the control, not a protocol.' % (
              fmt(float(np.mean(sw[0.5]['share']['ctrl_best']))), fmt(float(np.mean(sw[0.5]['share']['ours']))), fmt(b5[0]), fmt(b5[1]), fmt(b5[2])),
          '* The count is not bought for free. Clean-test RMSE: ours %s, control (validation) %s, control (oracle) %s; ours is worse than the validation control in %s of the runs. Recall of the true '
          'bump centres within 0.2: ours %s, control %s. This is the opposite direction of the count result and must not be described as better model selection.' % (
              fmt(rm['ours']), fmt(rm['ctrl']), fmt(rm['ctrl_or']), fmt(worse, 2), fmt(rec_o), fmt(rec_c)),
          '* Stability and budget: ours is exact on all 5 edges in %d of %d runs at delta = 0 (control %d of %d); the sd of ours\' exact-edge share over seeds is %s against %s for the control '
          '(mean over target x noise cells); ours reached the %d-step cap in %d of %d runs, so its count partly reflects the step budget.' % (
              int(round(sw[0.0]['all5']['ours'] * nr)), nr, int(round(sw[0.0]['all5']['ctrl'] * nr)), nr, fmt(sd_o), fmt(sd_c), K.STEPS, cap, nr),
          '* Limits of this analysis: it is POST HOC; it uses only the stored centres; the radius delta is a choice (the true centres are at least 1.0 apart, the initial width and grid spacing are '
          '5/9 = 0.556, so a radius between about 0.3 and 0.55 is a natural "same bump" scale); the control has no merge step in its protocol, and a merge rule applied to ours is not part of its '
          'protocol either; only one rule (single linkage on the centres) and one activity threshold (tau = %g) are used; the real-data part of the spike (verdict %s, non-inferior on %d of %d sets) '
          'is untouched by it.' % (K.TAU, 'GO' if ver['go'] else 'KILL', ver['n_ni'], ver['n_sets']), '']
    return '\n'.join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--from-raw', nargs='?', const='', required=True)
    ap.add_argument('--md')
    a = ap.parse_args()
    outd = os.environ.get('KAN_OUT', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'ml'))
    meta, recs = K.load_raw(a.from_raw or os.path.join(outd, 'kan_raw.npz'))
    txt = report(meta, recs)
    md = a.md or os.path.join(outd, 'kan_merge.md')
    with open(md, 'w') as f:
        f.write(txt + '\n')


if __name__ == '__main__':
    main()
