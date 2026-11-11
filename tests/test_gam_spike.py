"""Tests of the pure parts of experiments/exp_gam_spike.py (the torch / interpret / xgboost parts run on a private Kaggle notebook;
the torch-dependent tests below are skipped where torch is absent)."""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
import exp_gam_spike as G  # noqa: E402

OUTD = os.environ.get('GAM_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'ml'))
REPO_ML = os.path.join(os.path.dirname(__file__), '..', 'results', 'ml')


def _res(name):
    """A results file: next to the raw (GAM_OUT) if it is there, else in the repository's results/ml (the mutation harness copies it there)."""
    for base in (OUTD, REPO_ML):
        if os.path.exists(os.path.join(base, name)):
            return os.path.join(base, name)
    return None


def test_choose_datasets_rule():
    info = {}
    for name, _, kind in G.CANDIDATES:
        info[name] = dict(n=5000, d_enc=10, kind=kind)
    chosen = G.choose_datasets(info)
    regs = [n for n, _, k in G.CANDIDATES if k == 'reg'][:8]
    clfs = [n for n, _, k in G.CANDIDATES if k == 'clf'][:8]
    assert chosen == regs + clfs
    info[regs[0]]['n'] = G.N_MIN - 1                      # too small: skipped, the next candidate takes its place
    info[regs[1]]['d_enc'] = G.D_MAX + 1                  # too wide
    info[regs[2]]['d_enc'] = 2                            # too narrow
    assert G.choose_datasets(dict(info, **{regs[1]: dict(n=5000, d_enc=G.D_MAX, kind='reg')}), {'reg': 8, 'clf': 0})[0] == regs[1]   # D_MAX is inclusive
    del info[clfs[0]]                                     # did not load
    chosen = G.choose_datasets(info)
    allreg = [n for n, _, k in G.CANDIDATES if k == 'reg']
    allclf = [n for n, _, k in G.CANDIDATES if k == 'clf']
    assert chosen[:8] == allreg[3:11] and chosen[8:] == allclf[1:9]
    assert G.choose_datasets(info, {'reg': 1, 'clf': 0}) == [allreg[3]]
    info[allreg[3]]['n'] = G.N_MIN                        # the boundary is inclusive
    assert G.choose_datasets(info, {'reg': 1, 'clf': 0}) == [allreg[3]]
    info[allreg[3]]['kind'] = 'clf'                       # a mislabelled kind never matches
    assert G.choose_datasets(info, {'reg': 1, 'clf': 0}) == [allreg[4]]


def test_assign_roles():
    names = [f'r{i}' for i in range(8)] + [f'c{i}' for i in range(8)]
    kinds = {n: ('reg' if n[0] == 'r' else 'clf') for n in names}
    r = G.assign_roles(names, kinds)
    assert sum(v == 'dev' for v in r.values()) == 4 and sum(v == 'test' for v in r.values()) == 12
    assert sum(r[n] == 'dev' for n in names if kinds[n] == 'reg') == 2
    assert r == G.assign_roles(list(reversed(names)), kinds)          # independent of the order given
    assert r != G.assign_roles(names, kinds, seed=1) or True
    assert G.assign_roles(names, kinds, n_dev={'reg': 3, 'clf': 0}).values().__len__() == 16
    r2 = G.assign_roles(names, kinds, n_dev={'reg': 3, 'clf': 0})
    assert sum(v == 'dev' for v in r2.values()) == 3


def test_subsample_and_split():
    i = G.subsample_idx(25000)
    assert len(i) == G.CAP and np.all(np.diff(i) > 0) and np.array_equal(i, G.subsample_idx(25000))
    assert np.array_equal(G.subsample_idx(500), np.arange(500))
    tr, va, te = G.split_idx(1000, 3)
    assert (len(tr), len(va), len(te)) == (600, 200, 200)
    assert len(set(tr) | set(va) | set(te)) == 1000
    assert np.array_equal(G.split_idx(1000, 3)[2], te) and not np.array_equal(G.split_idx(1000, 4)[2], te)


def test_encode_columns():
    nan = float('nan')
    cols = [('a', np.array([1.0, nan, 3.0, 5.0]), False), ('const', np.ones(4), False),
            ('two', np.array(['x', 'y', 'y', None], dtype=object), True),
            ('three', np.array(['p', 'q', 'r', 'p'], dtype=object), True),
            ('many', np.array([str(i) for i in range(4)], dtype=object), True),
            ('one', np.array(['z'] * 4, dtype=object), True)]
    X, names = G.encode_columns(cols, max_levels=3)
    assert names == ['a', 'two=y', 'three=p', 'three=q', 'three=r']
    assert np.allclose(X[:, 0], [1, 3, 3, 5])                              # missing -> median
    assert np.allclose(X[:, 1], [0, 1, 1, 0]) and np.allclose(X[:, 2:].sum(1), 1)
    X0, n0 = G.encode_columns([('c', np.ones(3), False)])
    assert X0.shape == (3, 0) and n0 == []


def test_quantizer():
    rng = np.random.RandomState(0)
    X = np.column_stack([rng.exponential(size=2000), rng.randint(0, 2, 2000).astype(float), rng.randint(0, 4, 2000).astype(float)])
    Q = G.Quantizer().fit(X)
    Z = Q.transform(X)
    assert Z.min() >= 0 and Z.max() <= 1
    assert set(np.unique(Z[:, 1])) == {0.0, 1.0}
    assert abs(np.mean(Z[:, 0]) - 0.5) < 0.02 and np.all(np.diff(Z[np.argsort(X[:, 0]), 0]) >= -1e-12)     # uniform, monotone
    assert abs(np.sort(Q.transform(np.array([[1e9, 0, 0]]))[:, 0])[0] - 1.0) < 1e-12                     # beyond the data -> 1
    assert len(np.unique(Z[:, 2])) == 4
    Zc = G.Quantizer().fit(np.ones((5, 1))).transform(np.ones((5, 1)))
    assert np.all(Zc == 0)


def test_metrics():
    y = np.array([0, 0, 1, 1, 1.0])
    assert G.auc_score(y, np.array([0.1, 0.4, 0.35, 0.8, 0.9])) == pytest.approx(5 / 6)
    assert G.auc_score(y, np.zeros(5)) == pytest.approx(0.5)
    assert G.auc_score(y, -np.arange(5.0)) == pytest.approx(0.0 + (0.0))   # reversed ordering
    assert np.isnan(G.auc_score(np.ones(3), np.arange(3.0)))
    assert G.rel_rmse([1, 2, 3], [1, 2, 5]) == pytest.approx(np.sqrt(4 / 3))
    assert G.score_of('reg', 0.3) == 0.3 and G.score_of('clf', 0.8) == pytest.approx(0.2)


def test_counts():
    w = np.array([[0.0, 0.02, -0.5], [0.01, -0.011, 0.0]])
    m = G.active_mask(w)
    assert m.sum() == 3 and not m[1, 0] and m[1, 1]                        # |w| > THETA strictly
    assert G.ours_params(m) == 6
    assert G.segments([1, 1, 1]) == 1 and G.segments([1, 2, 2, 3]) == 3 and G.segments([]) == 0 and G.segments([1, 1 + 1e-12]) == 1
    ts = [np.array([0, 1, 1, 2, 2, 0.0]), np.array([0, 5, 5, 5, 5, 0.0])]       # the end bins (missing / unknown) are dropped
    assert G.ebm_params(ts) == (2 * 2 - 1) + (2 * 1 - 1)
    assert G.nam_params(3, 2) == 3 * 13 + 1
    assert G.nam1_params(3, 2) == 3 * 7 + 1 and G.nam1_params(1, 1) == 5
    assert G.nam_matched_h(3, 22) == 2 and G.nam_matched_h(3, 23) == 3 and G.nam_matched_h(3, 10 ** 9) == G.NAM1_LADDER[-1]
    assert G.nam_matched_h(3, G.nam1_params(3, 4)) == 4                        # >= : an exact match is accepted
    assert G.nam_matched_h(1, 1) == 1


def test_n_extrema():
    g = np.linspace(0, 1, 100)
    assert G.n_extrema(np.zeros(100)) == 0 and G.n_extrema(g) == 0
    assert G.n_extrema(np.sin(2 * np.pi * g)) == 2                            # a maximum and a minimum inside, counted at reversals
    assert G.n_extrema(np.sin(4 * np.pi * g)) == 4
    noisy = g + 0.001 * np.sin(200 * g)
    assert G.n_extrema(noisy) == 0                                            # small wiggles are ignored
    assert G.n_extrema(-np.sin(2 * np.pi * g)) == 2


def test_shape_stability():
    g = np.linspace(0, 1, 50)
    a = np.sin(2 * np.pi * g)
    same = np.stack([np.stack([a, np.zeros(50)]) for _ in range(3)])
    s, used = G.shape_stability(same)
    assert s == pytest.approx(1.0) and used == 1                              # the flat feature is skipped
    flip = same.copy()
    flip[1, 0] = -a
    s2, _ = G.shape_stability(flip)
    assert s2 < 0.5
    one_flat = same.copy()
    one_flat[2, 0] = 0
    s3, _ = G.shape_stability(one_flat)
    assert s3 == pytest.approx((1 + 0 + 0) / 3)                               # pairs (0,1)=1, (0,2)=0, (1,2)=0
    assert np.isnan(G.shape_stability(np.zeros((3, 2, 10)))[0])


def test_paired_p_and_gap():
    assert G.paired_one_sided_p([1, 1, 1, 1, 1]) == 0.0 and G.paired_one_sided_p([-1, -1, -1]) == 1.0
    assert G.paired_one_sided_p([1.0]) == 1.0
    p = G.paired_one_sided_p([0.1, 0.2, 0.15, 0.12, 0.18])
    assert 0 < p < 0.01
    assert G.paired_one_sided_p([0.1, -0.2, 0.15, -0.12, 0.05]) > 0.2
    assert G.rel_gap('reg', 0.5, 0.51) == pytest.approx(0.02) and G.rel_gap('clf', 0.1, 0.103) == pytest.approx(0.003)


def _sc(kind, ours, nam, matched, ebm, c1, c2):
    f = (lambda v: np.array(v, float))
    return {'ours_gd': f(ours), 'nam': f(nam), 'nam_matched': f(matched), 'ebm': f(ebm), 'ctrl_l1': f(c1), 'ctrl_l2': f(c2)}


def test_criteria_boundaries():
    par = dict(ours_gd=100, nam=1000, nam_matched=120, ebm=1000, ctrl_l1=300, ctrl_l2=800)
    five = [0.5] * 5
    # A: within 2 % of the better of NAM / EBM and <= 0.7 x its parameters
    c = G.criteria_for_dataset('reg', _sc('reg', [0.5099] * 5, [0.5] * 5, [0.6] * 5, [0.52] * 5, five, five), par)
    assert c['best'] == 'nam' and c['noninferior'] and c['fewer'] and c['A']
    c = G.criteria_for_dataset('reg', _sc('reg', [0.5101] * 5, [0.5] * 5, [0.6] * 5, [0.52] * 5, five, five), par)
    assert not c['noninferior'] and not c['A']
    c = G.criteria_for_dataset('reg', _sc('reg', five, [0.5] * 5, five, five, five, five), dict(par, ours_gd=700, nam=1000))
    assert c['fewer'] and c['A']                                                  # 700 = 0.7 x 1000: inclusive
    c = G.criteria_for_dataset('reg', _sc('reg', five, [0.5] * 5, five, [0.5] * 5, five, five), dict(par, ours_gd=701))
    assert not c['fewer'] and not c['A']
    # classification margins are absolute on 1 - AUC
    c = G.criteria_for_dataset('clf', _sc('clf', [0.105] * 5, [0.1] * 5, [0.2] * 5, [0.1] * 5, [0.1] * 5, [0.1] * 5), par)
    assert c['noninferior'] and c['best'] == 'nam'
    c = G.criteria_for_dataset('clf', _sc('clf', [0.1051] * 5, [0.1] * 5, [0.2] * 5, [0.1] * 5, [0.1] * 5, [0.1] * 5), par)
    assert not c['noninferior']
    # B: the matched NAM, one-sided paired t, and it must have at least ours' parameters
    ours = [0.50, 0.51, 0.49, 0.50, 0.52]
    c = G.criteria_for_dataset('reg', _sc('reg', ours, ours, [o + 0.03 + 0.001 * i for i, o in enumerate(ours)], ours, ours, ours), par)
    assert c['B'] and c['p_B'] < 0.05
    c = G.criteria_for_dataset('reg', _sc('reg', ours, ours, [o + 0.03 for o in ours], ours, ours, ours), dict(par, nam_matched=99))
    assert not c['B']                                                            # fewer parameters than ours: not a matched comparison
    c = G.criteria_for_dataset('reg', _sc('reg', ours, ours, [o + 0.03 * (-1) ** i for i, o in enumerate(ours)], ours, ours, ours), par)
    assert not c['B']
    # C: accuracy against the better control, tie band 0.5 % (regression) / 0.002 (AUC), parsimony 0.7 x the L1 control
    c = G.criteria_for_dataset('reg', _sc('reg', [0.5] * 5, five, five, five, [0.5026] * 5, [0.6] * 5), par)
    assert c['acc'] == 'win' and c['C'] == 'win'
    c = G.criteria_for_dataset('reg', _sc('reg', [0.5] * 5, five, five, five, [0.5024] * 5, [0.6] * 5), par)
    assert c['acc'] == 'tie' and c['par_win'] and c['C'] == 'win'                # tie in accuracy and 100 <= 0.7 * 300
    c = G.criteria_for_dataset('reg', _sc('reg', [0.5] * 5, five, five, five, [0.5024] * 5, [0.6] * 5), dict(par, ctrl_l1=100))
    assert c['acc'] == 'tie' and not c['par_win'] and c['C'] == 'tie'
    c = G.criteria_for_dataset('reg', _sc('reg', [0.5026] * 5, five, five, five, [0.5] * 5, [0.6] * 5), par)
    assert c['acc'] == 'loss' and c['C'] == 'loss'
    c = G.criteria_for_dataset('reg', _sc('reg', [0.5] * 5, five, five, five, [0.9] * 5, [0.5] * 5), dict(par, ctrl_l1=100))   # the better control is the L2 one
    assert c['ctrl_best'] == 'ctrl_l2' and c['acc'] == 'tie' and c['C'] == 'tie'                      # parsimony is judged against the L1 control
    c = G.criteria_for_dataset('clf', _sc('clf', [0.1] * 5, five, five, five, [0.1021] * 5, [0.2] * 5), par)
    assert c['acc'] == 'win'
    c = G.criteria_for_dataset('clf', _sc('clf', [0.1] * 5, five, five, five, [0.1019] * 5, [0.2] * 5), par)
    assert c['acc'] == 'tie'


def test_verdict_thresholds():
    def mk(nA, nB, nC, n=12):
        return {f'd{i}': dict(A=i < nA, B=i < nB, C='win' if i < nC else 'tie') for i in range(n)}
    assert G.verdict_from(mk(8, 0, 7))['go'] and G.verdict_from(mk(0, 8, 7))['go']
    assert not G.verdict_from(mk(7, 7, 12))['go'] and not G.verdict_from(mk(12, 12, 6))['go']
    assert G.verdict_from(mk(8, 0, 7)) == dict(nA=8, nB=0, nC=7, n=12, go=True)


def test_pick_configs():
    recs = []
    for ds, base in (('a', 1.0), ('b', 10.0)):
        for cid, mult in ((0, 1.2), (1, 1.0), (2, 1.5)):
            for s in (0, 1):
                recs.append(dict(method='m', cfg_id=cid, dataset=ds, seed=s, score=base * mult * (1 + 0.1 * s)))
    recs += [dict(method='n', cfg_id=0, dataset='a', seed=0, score=1.0), dict(method='n', cfg_id=1, dataset='a', seed=0, score=1.0)]
    p = G.pick_configs(recs)
    assert p == {'m': 1, 'n': 0}                                                  # scale-free over datasets; ties go to the lowest id
    recs2 = [r for r in recs if not (r['method'] == 'm' and r['dataset'] == 'b' and r['cfg_id'] == 1)]
    recs2 += [dict(method='m', cfg_id=1, dataset='b', seed=0, score=100.0)]
    assert G.pick_configs(recs2)['m'] == 0


def _synthetic(seed=0, ds=('r1', 'r2', 'c1'), seeds=(0, 1, 2, 3, 4), ours_better=True):
    rng = np.random.RandomState(seed)
    recs, shapes = [], {}
    for d in ds:
        kind = 'clf' if d.startswith('c') else 'reg'
        for m in G.METHODS:
            sh = []
            for s in seeds:
                base = 0.5 if kind == 'reg' else 0.8
                met = base + (0.01 if (m == 'ours_gd') else 0.0) * (-1 if kind == 'reg' else 1) * (1 if ours_better else -1) + 0.002 * rng.randn()
                recs.append(dict(dataset=d, kind=kind, role='test', seed=s, method=m, cfg_id=0, metric=met, score=G.score_of(kind, met),
                                 params=(20 if m.startswith('ours') else 500), units=[2, 0, 3] if m.startswith(('ours', 'ctrl')) else None, steps=10,
                                 time=1.0, d=3, n_used_feats=2, extrema=1.5, extra={'attempts': 0}))
                if m in ('ours_gd', 'nam', 'ebm'):
                    g = np.linspace(0, 1, G.GRID)
                    sh.append(np.stack([np.sin(2 * np.pi * g) + 0.1 * rng.randn(G.GRID), np.zeros(G.GRID), g - 0.5]).astype(np.float32))
            if sh:
                shapes[f'shp__{d}__{m}'] = np.stack(sh)
    return recs, shapes


def test_report_roundtrip_and_verdict(tmp_path, monkeypatch):
    recs, shapes = _synthetic()
    meta = dict(note='synthetic')
    p = tmp_path / 'gam_raw.npz'
    G.save_raw(str(p), recs, shapes, meta)
    recs2, shapes2, meta2 = G.load_records(str(p))
    assert recs2 == [json.loads(json.dumps(r, sort_keys=True)) for r in recs] and meta2 == meta
    assert set(shapes2) == set(shapes)
    md1 = G.render_report(recs2, shapes2, meta2)
    assert md1 == G.render_report(recs2, shapes2, meta2)
    assert '**GO**' in md1 or '**KILL**' in md1
    monkeypatch.setattr(G, 'OUT', str(tmp_path))
    G.main(['--from-raw'])
    assert (tmp_path / 'gam.md').read_text() == md1
    # with ours slightly worse than everything the verdict has to be KILL
    recs3, shapes3 = _synthetic(ours_better=False)
    assert '**KILL**' in G.render_report(recs3, shapes3, meta)
    # a record that is an error is ignored, not counted
    recs4 = recs + [dict(dataset='r1', kind='reg', seed=9, method='ebm', cfg_id=0, error='boom')]
    assert 'Errors: 1.' in G.render_report(recs4, shapes, meta)


def test_registry_matches_rule():
    path = _res('gam_probe.json')
    if path is None or not G.DATASETS:
        pytest.skip('probe / registry not frozen yet')
    info = json.load(open(path))
    chosen = G.choose_datasets({k: v for k, v in info.items() if v.get('ok')})
    assert [d[0] for d in G.DATASETS] == chosen
    roles = G.assign_roles(chosen, {n: info[n]['kind'] for n in chosen})
    assert {d[0]: d[3] for d in G.DATASETS} == roles
    assert sum(d[3] == 'test' for d in G.DATASETS) == 12 and sum(d[3] == 'dev' for d in G.DATASETS) == 4


def test_frozen_matches_dev(tmp_path):
    fz, dev = _res('gam_frozen.json'), _res('gam_dev_all.jsonl')
    if fz is None or dev is None:
        pytest.skip('not frozen yet')
    recs = [json.loads(l) for l in open(dev) if l.strip()]
    pick = G.pick_configs([r for r in recs if 'error' not in r and r['method'] in G.GRIDS])      # the first-round gradient-descent controls are not candidates
    frozen = json.load(open(fz))
    for m, cid in pick.items():
        assert frozen[m] == G.GRIDS[m][cid] and frozen[m + '_id'] == cid


# ---- torch parts (skipped locally) -----------------------------------------------------------------------------------------------
def _toy(n=600, d=4, seed=0):
    rng = np.random.RandomState(seed)
    X = rng.rand(n, d).astype(np.float32)
    y = 2 * np.exp(-0.5 * ((X[:, 0] - 0.4) / 0.07) ** 2) + 0.5 * (X[:, 1] - 0.5)
    y = (y - y.mean()) / y.std() + 0.1 * rng.randn(n)
    y = y.astype(np.float32)
    a, b = 360, 480
    return dict(kind='reg', d=d, n=n, Xtr=X[:a], ytr=y[:a], Xva=X[a:b], yva=y[a:b], Xte=X[b:], yte=y[b:], b0=0.0)


def test_pins_and_extra_kills():
    # the split is a pinned function of (n, seed): the offset 1000 and the order train / validation / test matter
    tr, va, te = G.split_idx(20, 0)
    assert list(te) == [1, 0, 7, 19] and list(va) == [8, 9, 17, 13]
    assert (len(G.split_idx(1003, 0)[1]), len(G.split_idx(1003, 0)[2])) == (200, 201)
    # missing numeric values go to the median (not the mean)
    X, _ = G.encode_columns([('a', np.array([1.0, np.nan, 2.0, 9.0]), False)])
    assert X[1, 0] == 2.0
    # a three-valued column is not treated as binary
    Q = G.Quantizer().fit(np.array([[0.0], [1.0], [2.0], [2.0], [1.0], [0.0]]))
    assert len(np.unique(Q.transform(np.array([[0.0], [1.0], [2.0]])))) == 3
    # a wiggle below the tolerance is ignored, one above it counts
    g = np.linspace(0, 1, 400)
    assert G.n_extrema(g + 0.02 * np.sin(200 * g)) == 0 and G.n_extrema(g + 0.2 * np.sin(200 * g)) > 5
    # the stability is weighted by the size of the shapes
    t = np.linspace(0, 1, 50)
    big = np.sin(2 * np.pi * t) * 10
    small = np.stack([np.sin(2 * np.pi * t), np.sin(2 * np.pi * t), -np.sin(2 * np.pi * t)]) * 0.1
    sh = np.stack([np.stack([big, small[k]]) for k in range(3)])
    s, used = G.shape_stability(sh)
    w_big, w_small = 10 * np.sin(2 * np.pi * t).std(), 0.1 * np.sin(2 * np.pi * t).std()
    assert used == 2 and s == pytest.approx((w_big * 1.0 + w_small * (1 - 1 - 1) / 3) / (w_big + w_small))
    # pick_configs is scale-free over datasets and takes the mean over seeds
    r = []
    for ds, sc in (('a', (1.0, 1.5)), ('b', (100.0, 80.0))):
        for cid, v in enumerate(sc):
            for sd in (0, 1):
                r.append(dict(method='m', cfg_id=cid, dataset=ds, seed=sd, score=v))
    assert G.pick_configs(r) == {'m': 0}
    r2 = [dict(method='m', cfg_id=0, dataset='a', seed=0, score=1.0), dict(method='m', cfg_id=0, dataset='a', seed=1, score=1.0),
          dict(method='m', cfg_id=1, dataset='a', seed=0, score=0.5), dict(method='m', cfg_id=1, dataset='a', seed=1, score=3.0)]
    assert G.pick_configs(r2) == {'m': 0}
    # the paired p-value is the one-sided paired t-test
    from scipy.stats import ttest_1samp
    for d in ([0.1, 0.05, 0.2, -0.05, 0.12], [0.3, 0.1, 0.25, 0.0, 0.2], [-0.1, 0.0, 0.3]):
        assert G.paired_one_sided_p(d) == pytest.approx(ttest_1samp(d, 0, alternative='greater').pvalue, rel=1e-9)
    # B needs p < 0.05: p = 0.0559 does not pass
    ours = np.array([0.5, 0.51, 0.49, 0.5, 0.52])
    par = dict(ours_gd=100, nam=1000, nam_matched=120, ebm=1000, ctrl_l1=300, ctrl_l2=800)
    c = G.criteria_for_dataset('reg', _sc('reg', ours, ours, ours + np.array([0.1, 0.05, 0.2, -0.05, 0.12]), ours, ours, ours), par)
    assert 0.05 < c['p_B'] < 0.06 and not c['B']
    # A compares with the better of NAM and EBM
    c = G.criteria_for_dataset('reg', _sc('reg', [0.52] * 5, [0.6] * 5, [0.6] * 5, [0.5] * 5, [0.5] * 5, [0.5] * 5), par)
    assert c['best'] == 'ebm' and not c['noninferior']


def test_report_with_incomplete_dataset(tmp_path):
    recs, shapes = _synthetic()
    recs = [r for r in recs if not (r['dataset'] == 'c1' and r['method'] == 'nam_matched')]
    md = G.render_report(recs, shapes, {})
    assert 'on 0 of 2' in md or 'of 2 (needs' in md                                  # the criteria are evaluated on the datasets with every method
    assert '| c1 | nam' not in md.split('## Pre-registered criteria')[1]


def test_posthoc_keywords():
    par = dict(ours_gd=100, nam=1000, nam_matched=120, ebm=1000, ctrl_l1=300, ctrl_l2=800)
    five = [0.5] * 5
    # tie band scaling: a 0.6 % gap is a loss at the pre-registered band (0.5 %) and a tie at twice the band
    sc = _sc('reg', [0.5] * 5, five, five, five, [0.497] * 5, [0.6] * 5)
    assert G.criteria_for_dataset('reg', sc, par)['acc'] == 'loss'
    assert G.criteria_for_dataset('reg', sc, par, tie_scale=2.0)['acc'] == 'tie'
    # the control set: with only the weaker control the same ours wins
    sc = _sc('reg', [0.5] * 5, five, five, five, [0.52] * 5, [0.499] * 5)
    assert G.criteria_for_dataset('reg', sc, par)['ctrl_best'] == 'ctrl_l2'
    assert G.criteria_for_dataset('reg', sc, par, ctrls=('ctrl_l1',))['ctrl_best'] == 'ctrl_l1'
    assert G.criteria_for_dataset('reg', sc, par, ctrls=('ctrl_l1',))['acc'] == 'win'
    # margin scaling and the parameter condition
    sc = _sc('reg', [0.52] * 5, [0.5] * 5, five, [0.55] * 5, five, five)
    assert not G.criteria_for_dataset('reg', sc, par)['noninferior'] and G.criteria_for_dataset('reg', sc, par, margin_scale=2.0)['noninferior']
    sc = _sc('reg', [0.5] * 5, [0.5] * 5, five, [0.5] * 5, five, five)
    assert not G.criteria_for_dataset('reg', sc, dict(par, ours_gd=900))['A']
    assert G.criteria_for_dataset('reg', sc, dict(par, ours_gd=900), a_params=False)['A']
    # B at p < 0.10
    ours = np.array([0.5, 0.51, 0.49, 0.5, 0.52])
    c0 = G.criteria_for_dataset('reg', _sc('reg', ours, ours, ours + np.array([0.1, 0.05, 0.2, -0.05, 0.12]), ours, ours, ours), par)
    c1 = G.criteria_for_dataset('reg', _sc('reg', ours, ours, ours + np.array([0.1, 0.05, 0.2, -0.05, 0.12]), ours, ours, ours), par, p_alpha=0.10)
    assert not c0['B'] and c1['B']


def test_report_posthoc_section():
    recs, shapes = _synthetic()
    md = G.render_report(recs, shapes, {})
    sec = md.split('## POST HOC sensitivity')[1].split('## Supplementary')[0]
    assert 'pre-registered rule' in sec and 'C: tie band x2' in sec and 'against ctrl_l1 only' in sec and 'B: p < 0.10' in sec
    assert 'Ratio of mean stored parameters ours_gd / ctrl_l1: 0.04' in sec        # 20 vs 500 in the synthetic records
    assert 'ctrl_l1 3 of 3' in sec


def _crafted_records():
    """Three regression datasets with hand-set per-seed metrics and parameters (lower = better), to pin the post hoc table."""
    seeds = range(5)
    off = np.array([0.0, 0.01, -0.01, 0.0, 0.02])
    d = np.array([0.1, 0.05, 0.2, -0.05, 0.12])
    spec = {
        'r1': {'ours_gd': (0.5 + off, 100), 'ctrl_l1': (0.497 + off, 200), 'ctrl_l2': (0.6 + off, 300), 'nam': (0.6 + off, 1000),
               'nam_matched': (0.5 + off + d, 120), 'ebm': (0.4875 + off, 1000)},
        'r2': {'ours_gd': (0.5 + off, 100), 'ctrl_l1': (0.497 + off, 200), 'ctrl_l2': (0.6 + off, 300), 'nam': (0.6 + off, 1000),
               'nam_matched': (0.5 + off, 120), 'ebm': (0.6 + off, 1000)},
        'r3': {'ours_gd': (0.5 + off, 100), 'ctrl_l1': (0.55 + off, 200), 'ctrl_l2': (0.497 + off, 300), 'nam': (0.6 + off, 1000),
               'nam_matched': (0.5 + off, 120), 'ebm': (0.6 + off, 1000)},
    }
    recs = []
    for ds, ms in spec.items():
        for m in G.METHODS:
            met, par = ms.get(m, (0.7 + off, 50))
            for i, s in enumerate(seeds):
                recs.append(dict(dataset=ds, kind='reg', role='test', seed=s, method=m, cfg_id=0, metric=float(met[i]), score=float(met[i]), params=par,
                                 units=[1, 1] if m.startswith(('ours', 'ctrl')) else None, steps=1, time=1.0, d=2, n_used_feats=1, extrema=1.0, extra={}))
    return recs


def test_posthoc_table_values():
    md = G.render_report(_crafted_records(), {}, {})
    sec = md.split('## POST HOC sensitivity')[1].split('## Supplementary')[0]
    rows = {l.split('|')[1].strip(): [c.strip() for c in l.split('|')[2:-1]] for l in sec.splitlines() if l.startswith('| ') and '---' not in l}
    # columns: A, B, C wins, C losses, would be GO
    # A: r1 is 2.56 % behind the better of NAM / EBM (outside 2 %), r2 and r3 are ahead with 10 % of the parameters; B: p = 0.056 on r1 only;
    # C: ours is 0.6 % behind the better control on every dataset (a loss at the 0.5 % band; the better one is ctrl_l1 on r1, r2 and ctrl_l2 on r3)
    assert rows['pre-registered rule'] == ['2', '0', '0', '3', 'no']
    assert rows['C: tie band x2 (1 % / 0.004)'][2:4] == ['3', '0']
    assert rows['C: tie band x4 (2 % / 0.008)'][2:4] == ['3', '0']
    assert rows['C: against ctrl_l1 only'][2:4] == ['1', '2']
    assert rows['C: against ctrl_l2 only'][2:4] == ['2', '1']
    assert rows['A: margin x1.5, no parameter condition'][0] == '3'
    assert rows['B: p < 0.10'][1] == '1' and rows['pre-registered rule'][1] == '0'
    assert 'ctrl_l1 1 of 3' in sec and 'ctrl_l2 2 of 3' in sec
    assert 'Ratio of mean stored parameters ours_gd / ctrl_l1: 0.50' in sec


def test_dev_summary_in_report():
    dev = []
    for ds, base in (('a', 1.0), ('b', 10.0)):
        for cid, mult, par in ((0, 1.2, 10), (1, 1.0, 20)):
            for s in (0, 1):
                dev.append(dict(method='ours_gd', cfg_id=cid, dataset=ds, seed=s, score=base * mult, params=par))
    dev.append(dict(method='ctrlgd_l1', cfg_id=0, dataset='a', seed=0, score=1.0, params=77))
    recs, shapes = _synthetic()
    md = G.render_report(recs, shapes, dict(dev=dev))
    sec = md.split('## Development (tuning) summary')[1].split('## POST HOC')[0]
    assert '| ours_gd | 1* |' in sec and '| ours_gd | 0 |' in sec
    assert '1.2000 | 10 |' in sec and '1.0000 | 20 |' in sec
    assert '| ctrlgd_l1 | 0 | first round, plain gradient descent | - | 77 |' in sec
    assert 'Development (tuning) summary' not in G.render_report(recs, shapes, {})


# ---- torch parts (skipped locally; run on the notebook) ---------------------------------------------------------------------------

def test_iss_selects_few_units_torch():
    pytest.importorskip('torch')
    D = _toy()
    r = G.run_method('ours_gd', D, dict(eps=1e-3, sigma=0.07, lr=0.5), 0)
    assert r['units'][2] == 0 and r['units'][3] == 0                 # features that carry nothing stay dormant
    assert 1 <= r['units'][0] <= 12 and 1 <= r['units'][1] <= 12      # the two informative features switch some bumps on, of 16
    assert r['params'] == 2 * sum(r['units'])
    assert G.rel_rmse(D['yte'], r['pred']) < 0.2
    # a smaller start switches fewer bumps on (the scale of the entry times), a larger one more
    small = G.run_method('ours_gd', D, dict(eps=1e-4, sigma=0.07, lr=1.0), 0)
    large = G.run_method('ours_gd', D, dict(eps=1e-2, sigma=0.05, lr=1.0), 0)
    assert sum(small['units']) < sum(r['units']) < sum(large['units']) + 5
    # Adam on the same parametrization switches on units everywhere
    ad = G.run_method('ours_adam', D, dict(eps=1e-4, sigma=0.07, lr=1e-2), 0)
    assert sum(ad['units']) > 2 * sum(r['units'])


def test_hand_gradients_match_autograd_torch():
    pytest.importorskip('torch')
    assert G.check_grads() < 1e-12
