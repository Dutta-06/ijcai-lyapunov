"""Tests of experiments/exp_kan_spike.py: the pure-numpy parts (synthetic targets, the activity rule and counts, the statistics, the real-set rule, the splits, the
verdict, the report and its byte-identical regeneration) with planted records; the torch parts are tested only where torch is installed (importorskip)."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments')
sys.path.insert(0, HERE)
import exp_kan_spike as E  # noqa: E402

SCRIPT = E.__file__


# ---------------------------------------------------------------------------------------------------------------------------------
# targets
# ---------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('tid', list(E.TEST_TARGETS) + list(E.DEV_TARGETS))
def test_target_constraints(tid):
    t = E.make_target(tid)
    K = t['K']
    assert len(K) == E.D_SYN and K.min() == 0 and K.max() >= 3 and 6 <= K.sum() <= 11 and K.max() <= 4
    for k, b in zip(K, t['bumps']):
        assert b.shape == (k, 3)
        if k:
            assert np.all((b[:, 1] >= -1.8) & (b[:, 1] <= 1.8)) and np.all((b[:, 2] >= 0.25) & (b[:, 2] <= 0.45))
            assert np.all((np.abs(b[:, 0]) >= 0.8) & (np.abs(b[:, 0]) <= 1.5))
        if k >= 2:
            c = b[:, 1]
            assert np.min(np.abs(c[:, None] - c[None]) + 9 * np.eye(k)) >= 1.0


def test_target_deterministic_and_distinct():
    a, b, c = E.make_target(3), E.make_target(3), E.make_target(4)
    assert a['K'].tolist() == b['K'].tolist() and a['scale'] == b['scale']
    assert not np.array_equal(a['K'], c['K']) or a['scale'] != c['scale']


def test_target_unit_variance_and_value():
    t = E.make_target(2)
    X = np.random.default_rng(5).uniform(E.X_LO, E.X_HI, size=(200000, E.D_SYN))
    assert abs(np.std(E.target_clean(t, X)) - 1.0) < 0.02
    # value of one planted bump
    t0 = dict(tid=0, K=np.array([1, 0]), bumps=[np.array([[2.0, 0.5, 0.3]]), np.zeros((0, 3))], scale=2.0)
    X = np.array([[0.5, 9.0], [0.8, 0.0]])
    assert np.allclose(E.target_clean(t0, X), [1.0, np.exp(-0.5)])


def test_sample_synth_noise_and_split():
    t = E.make_target(1)
    d = E.sample_synth(t, 0.4, 1)
    assert d['Xtr'].shape == (E.N_TRAIN, 5) and d['Xva'].shape == (E.N_VAL, 5) and d['Xte'].shape == (E.N_TEST_SYN, 5)
    assert np.all((d['Xtr'] >= E.X_LO) & (d['Xtr'] <= E.X_HI))
    assert abs(np.std(d['ytr'] - E.target_clean(t, d['Xtr'])) - 0.4) < 0.03
    assert abs(np.std(d['yva'] - E.target_clean(t, d['Xva'])) - 0.4) < 0.06
    assert np.allclose(d['yte'], E.target_clean(t, d['Xte']))            # test targets are clean
    d2 = E.sample_synth(t, 0.4, 1)
    assert np.array_equal(d['ytr'], d2['ytr'])
    assert not np.array_equal(d['Xtr'], E.sample_synth(t, 0.4, 2)['Xtr'])
    assert not np.array_equal(d['Xtr'], E.sample_synth(t, 0.1, 1)['Xtr'])
    assert not np.array_equal(d['Xtr'], E.sample_synth(E.make_target(2), 0.4, 1)['Xtr'])


# ---------------------------------------------------------------------------------------------------------------------------------
# activity rule
# ---------------------------------------------------------------------------------------------------------------------------------
def test_unit_contrib_and_rms():
    x = np.array([[0.0], [1.0]])
    w, c, s = np.array([[[2.0]]]), np.array([[[0.0]]]), np.array([[[1.0]]])
    U = E.unit_contrib(w, c, s, x)
    assert U.shape == (2, 1, 1, 1)
    assert np.allclose(U[:, 0, 0, 0], [2.0, 2.0 * np.exp(-0.5)])
    rms, sd, act = E.layer_activity(w, c, s, x)
    assert np.allclose(rms[0, 0, 0], np.sqrt(np.mean(U ** 2)))
    assert np.allclose(sd[0], np.std(U[:, 0, 0, 0]))


def test_activity_threshold_and_dormant_units():
    rng = np.random.default_rng(0)
    x = rng.uniform(-2.5, 2.5, (4000, 2))
    w = np.zeros((1, 2, 3))
    c = np.zeros((1, 2, 3))
    s = np.full((1, 2, 3), 0.4)
    w[0, 0, 0], c[0, 0, 0] = 1.0, -1.0           # strong unit on edge 0
    w[0, 0, 1], c[0, 0, 1] = 1.0, 1.0            # second strong unit on edge 0
    w[0, 1, 0], c[0, 1, 0] = 1.0, 0.0            # strong unit on edge 1
    w[0, 1, 1], c[0, 1, 1] = 1e-4, 1.0           # dormant
    rms, sd, act = E.layer_activity(w, c, s, x)
    assert act[0, 0, :].tolist() == [True, True, False] and act[0, 1, :].tolist() == [True, False, False]
    assert E.edge_counts(act).tolist() == [2, 1]
    assert rms[0, 1, 1] < 1e-3 * rms[0, 1, 0] * 2
    # a unit just above / below the threshold relative to the node's sd
    for wt, expect in ((0.5, True), (0.001, False)):
        w2, c2 = w.copy(), c.copy()
        w2[0, 1, 2], c2[0, 1, 2] = wt, 0.0
        r2, sd2, a2 = E.layer_activity(w2, c2, s, x)
        assert bool(a2[0, 1, 2]) == (r2[0, 1, 2] >= E.TAU * sd2[0]) and bool(a2[0, 1, 2]) == expect
    # exactly zero weight is never active, whatever tau
    _, _, a0 = E.layer_activity(np.zeros((1, 1, 1)), c[:, :1, :1], s[:, :1, :1], x[:, :1], tau=0.0)
    assert not a0.any()
    # the threshold relative to the node's own spread: scaling the whole layer does not change the activity pattern
    _, _, a_s = E.layer_activity(7.0 * w, c, s, x)
    assert np.array_equal(a_s, act)


def test_tau_monotone():
    rng = np.random.default_rng(1)
    x = rng.uniform(-2.5, 2.5, (1000, 3))
    w = rng.standard_normal((2, 3, 5)) * np.array([1, 0.3, 0.05, 0.01, 0.0])
    c = rng.uniform(-2, 2, (2, 3, 5))
    s = np.full((2, 3, 5), 0.5)
    n = [E.layer_activity(w, c, s, x, tau)[2].sum() for tau in (0.0, 0.01, 0.02, 0.05, 0.2)]
    assert n == sorted(n, reverse=True) and n[0] > n[-1]


def test_param_counts():
    assert E.active_params(7, 2) == 23 and E.active_params(7, 2, per_unit=1) == 9 and E.active_params(0, 1) == 1
    assert E.mlp_params(5, 3) == 5 * 3 + 3 + 9 + 3 + 3 + 1


# ---------------------------------------------------------------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------------------------------------------------------------
def test_boot_ci():
    m, lo, hi = E.boot_ci([0.2] * 10)
    assert m == pytest.approx(0.2) and lo == pytest.approx(0.2) and hi == pytest.approx(0.2)
    v = [0.1, 0.3, 0.2, 0.4, -0.1, 0.25, 0.15, 0.05, 0.3, 0.2]
    m, lo, hi = E.boot_ci(v)
    assert m == pytest.approx(np.mean(v)) and lo < m < hi
    assert (lo, hi) == E.boot_ci(v)[1:]                       # deterministic
    assert E.boot_ci(v, seed=1)[1:] != (lo, hi)
    mm, l2, h2 = E.boot_ci(v, n_boot=2000)
    assert abs(l2 - lo) < 0.03
    m0, lo0, hi0 = E.boot_ci([-0.2, -0.3, -0.1, -0.2])
    assert hi0 < 0
    # wide interval with spread data
    assert E.boot_ci([-1, 1, -1, 1, -1, 1])[1] < 0 < E.boot_ci([-1, 1, -1, 1, -1, 1])[2]


def test_best_at_most():
    pts = [(10, 0.5), (20, 0.3), (40, 0.35), (80, 0.2)]
    assert E.best_at_most(pts, 40) == (0.3, False)
    assert E.best_at_most(pts, 20) == (0.3, False)
    assert E.best_at_most(pts, 19) == (0.5, False)
    assert E.best_at_most(pts, 1000) == (0.2, False)
    assert E.best_at_most(pts, 5) == (0.5, True)               # none that small: the smallest member, flagged
    assert E.best_at_most([(7, 0.4), (7, 0.2)], 3) == (0.2, True)


def test_one_se_pick():
    cands = [dict(val=0.100, P=50), dict(val=0.104, P=30), dict(val=0.2, P=5), dict(val=0.109, P=20), dict(val=0.1, P=60)]
    # n_val = 500: tolerance sqrt(2/500) = 6.3 percent: allowed val <= 0.1063
    assert E.one_se_pick(cands, 500)['P'] == 30
    assert E.one_se_pick(cands, 500000)['P'] == 50
    assert E.one_se_pick(cands, 20)['P'] == 20                  # a wide tolerance: 0.1 * (1 + 0.316) = 0.1316
    assert E.one_se_pick([dict(val=0.1, P=10), dict(val=0.1, P=10)], 500)['P'] == 10


def test_split_real():
    tr, va, te = E.split_real(10000, 0)
    assert len(tr) == 2400 and len(va) == 800 and len(te) == 800 and len(set(tr) | set(va) | set(te)) == 4000
    tr2, va2, te2 = E.split_real(10000, 0)
    assert np.array_equal(tr, tr2) and not np.array_equal(tr, E.split_real(10000, 1)[0])
    tr, va, te = E.split_real(1000, 2)
    assert (len(tr), len(va), len(te)) == (600, 200, 200)


# ---------------------------------------------------------------------------------------------------------------------------------
# the real-set rule
# ---------------------------------------------------------------------------------------------------------------------------------
def fake_listing():
    names = ['a', 'b', 'fri_c0', 'KIN8NM', 'c', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']
    out = []
    for i, nm in enumerate(names):
        out.append((100 + 3 * i, nm, 2000, 8, 8, 0))
    out.append((5, 'sym', 2000, 8, 7, 1))                       # symbolic feature: dropped
    out.append((6, 'sym2', 2000, 8, 7, 0))                      # nn != nf: dropped
    return out


def test_select_real():
    def loaded(did):
        rng = np.random.default_rng(did)
        if did == 100 + 3 * 1:                                  # 'b': a classification target
            return rng.standard_normal((50, 3)), rng.integers(0, 3, 50).astype(float)
        if did == 100 + 3 * 6:                                  # 'd': load failure
            return None
        return rng.standard_normal((200, 3)), rng.standard_normal(200)
    dev, test = E.select_real(fake_listing(), loaded)
    names_all = [n for _, n in dev + test]
    assert 'fri_c0' not in names_all and 'KIN8NM' not in names_all and 'sym' not in names_all and 'b' not in names_all and 'd' not in names_all
    assert len(dev) == 3 and len(test) == 8
    qual = ['a', 'c', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
    assert [n for _, n in dev] == [qual[2], qual[5], qual[8]]
    assert [n for _, n in test] == [q for i, q in enumerate(qual) if (i + 1) not in E.DEV_RANKS]
    # ascending id even if the listing is shuffled
    lst = fake_listing()[::-1]
    assert E.select_real(lst, loaded) == (dev, test)
    # constant or nonfinite target is dropped
    def loaded2(did):
        X = np.random.default_rng(did).standard_normal((200, 3))
        return X, (np.ones(200) if did == 100 else np.random.default_rng(did).standard_normal(200))
    d2, t2 = E.select_real(fake_listing(), loaded2)
    assert 'a' not in [n for _, n in d2 + t2]
    def loaded3(did):
        X = np.random.default_rng(did).standard_normal((200, 3))
        if did == 100:
            X[0, 0] = np.nan
        return X, np.random.default_rng(did).standard_normal(200)
    assert 'a' not in [n for _, n in sum(E.select_real(fake_listing(), loaded3), [])]
    # few distinct target values: dropped
    def loaded4(did):
        return np.random.default_rng(did).standard_normal((200, 3)), (np.arange(200) % 49).astype(float) if did == 100 else np.random.default_rng(did).standard_normal(200)
    assert 'a' not in [n for _, n in sum(E.select_real(fake_listing(), loaded4), [])]
    def loaded5(did):
        return np.random.default_rng(did).standard_normal((200, 3)), (np.arange(200) % 50).astype(float) if did == 100 else np.random.default_rng(did).standard_normal(200)
    assert 'a' in [n for _, n in sum(E.select_real(fake_listing(), loaded5), [])]


# ---------------------------------------------------------------------------------------------------------------------------------
# planted records, verdict and report
# ---------------------------------------------------------------------------------------------------------------------------------
def planted(share_ours=1.0, share_ctrl=0.4, rmse_ratio=1.0, n_ni=8, seed=0, centres=False):
    """Synthetic records: 10 targets x 2 noise x 3 seeds; 8 real sets x 3 seeds (arch B) with ours at `rmse_ratio` x the best baseline on the first n_ni sets."""
    rng = np.random.default_rng(seed)
    recs = []
    for tid in E.TEST_TARGETS:
        K = E.make_target(tid)['K'].tolist()
        for nz in E.NOISES:
            for s in E.SEEDS:
                def counts(share):
                    c = list(K)
                    for i in range(len(c)):
                        if rng.random() > share:
                            c[i] += int(rng.integers(1, 3))
                    return c
                base = dict(kind='synth', tid=tid, noise=nz, seed=s, arch='A', truth=K, n_val=500, steps=100, secs=1.0, test_rmse=0.1, val_mse=0.01)
                def add(method, conf, share, P, val=0.01, rm=0.1):
                    r = dict(base, method=method, conf=json.dumps(conf, sort_keys=True), val_mse=val, test_rmse=rm)
                    for tau in E.TAUS:
                        c = counts(share)
                        r['counts_%g' % tau] = c
                        r['P_%g' % tau] = P + sum(c)
                        r['n_%g' % tau] = sum(c)
                    if centres:
                        tc = [b[:, 1].tolist() for b in E.make_target(tid)['bumps']]
                        r['truth_c'] = tc
                        r['centres_0.02'] = [(tc[i] + [9.0] * 5)[:r['counts_0.02'][i]] for i in range(len(K))]
                    recs.append(r)
                add('ours', dict(kind='gd'), share_ours, 10)
                for lam, sh in zip(E.LAMS, (share_ctrl,) * len(E.LAMS)):
                    add('ctrl', dict(lam=lam, lr=0.01), sh, 10, val=0.01 + 0.0001 * lam)
                for G in E.GS:
                    add('fixed', dict(G=G, lr=0.01), 0.1, 10, val=0.02 - 0.001 * G)
    for k in range(8):
        for s in E.SEEDS:
            ratio = rmse_ratio if k < n_ni else 1.5
            base = dict(kind='real', did=1000 + k, name='set%d' % k, seed=s, arch='B', d=6 + k, n_train=3600, n_val=1200, steps=50, secs=2.0)
            recs.append(dict(base, arch='H', method='ours', conf='{}', val_mse=0.1, test_rmse=0.4 * ratio, **{'P_%g' % t: 100 for t in E.TAUS}))
            for G, rm in zip(E.GS, (0.9, 0.7, 0.5, 0.4, 0.4, 0.4)):
                recs.append(dict(base, method='fixed', conf=json.dumps(dict(G=G, lr=0.01)), val_mse=rm, test_rmse=rm * 0.8, **{'P_%g' % t: 20 * G for t in E.TAUS}))
                recs.append(dict(base, method='fixed', conf=json.dumps(dict(G=G, lr=0.03)), val_mse=rm * 2, test_rmse=rm * 0.1, **{'P_%g' % t: 20 * G for t in E.TAUS}))
            for h, rm in ((2, 0.8), (8, 0.6)):
                recs.append(dict(base, method='mlp', conf=json.dumps(dict(h=h, lr=0.01)), val_mse=rm, test_rmse=rm, **{'P_%g' % t: 50 * h for t in E.TAUS}))
    return recs


def test_exact_share():
    r = {'counts_0.02': [1, 2, 0, 3], 'truth': [1, 2, 1, 3]}
    assert E.exact_share(r) == 0.75
    r = {'counts_0.05': [1, 1, 1], 'truth': [1, 1, 1]}
    assert E.exact_share(r, 0.05) == 1.0


def test_synth_summary_and_gain():
    recs = planted(1.0, 0.0)
    rows = E.synth_summary(recs)
    assert len(rows) == 60
    assert all(r['ours']['share'] == 1.0 for r in rows)
    assert np.mean([r['ctrl']['share'] for r in rows]) < 0.3
    m, lo, hi = E.paired_gain(rows, 'ours', 'ctrl')
    assert m > 0.7 and lo > 0.5
    # the oracle control is at least as good as the validation-picked one
    assert all(r['ctrl_or']['share'] >= r['ctrl']['share'] for r in rows)
    # equal shares: interval straddles zero or gain ~ 0
    rows2 = E.synth_summary(planted(0.5, 0.5, seed=3))
    m, lo, hi = E.paired_gain(rows2, 'ours', 'ctrl')
    assert m < 0.05 and lo < 0.0                                  # the control's best-of-five pick is no worse: no gain


def test_ctrl_pick_by_validation_with_one_se():
    # the control's candidates: the sparse one has slightly higher validation error but within one SE: it is picked; a much worse one is not
    recs = []
    K = [1, 0, 2, 0, 0]
    base = dict(kind='synth', tid=1, noise=0.1, seed=0, arch='A', truth=K, n_val=500, steps=1, secs=1.0, test_rmse=0.1)
    for lam, val, P, c in ((0.0, 0.100, 60, [3, 2, 3, 2, 2]), (1e-3, 0.103, 30, [1, 0, 2, 0, 0]), (1e-2, 0.300, 10, [0, 0, 1, 0, 0])):
        r = dict(base, method='ctrl', conf=json.dumps(dict(lam=lam)), val_mse=val)
        for tau in E.TAUS:
            r['counts_%g' % tau], r['P_%g' % tau] = c, P
        recs.append(r)
    r = dict(base, method='ours', conf='{}', val_mse=0.1)
    f = dict(base, method='fixed', conf=json.dumps(dict(G=4)), val_mse=0.2)
    for x in (r, f):
        for tau in E.TAUS:
            x['counts_%g' % tau], x['P_%g' % tau] = K, 20
    rows = E.synth_summary(recs + [r, f])
    assert rows[0]['ctrl']['conf'] == json.dumps(dict(lam=1e-3)) and rows[0]['ctrl']['share'] == 1.0
    assert rows[0]['ctrl_or']['share'] == 1.0 and rows[0]['fixed']['conf'] == json.dumps(dict(G=4))


def test_real_summary_matched_params():
    recs = planted(n_ni=8, rmse_ratio=1.0)
    real = E.real_summary(recs)
    assert len(real) == 8
    r = real[0]
    # ours has 100 active parameters; fixed grids with <= 100: G = 2, 3, 4 (P = 40, 60, 80); per capacity the lr with lower validation MSE (0.01) -> test 0.8*rm; mlp h=2 (P=100): 0.8
    assert r['who'] == ['fixed'] and r['base_mean'] == pytest.approx(0.4)           # G = 4: 0.5 * 0.8
    assert r['ratio'] == pytest.approx(1.0) and r['ni'] and r['P_cv'] == 0.0 and r['clamped'] == 0
    # worse than the margin
    real = E.real_summary(planted(n_ni=0, rmse_ratio=1.0))
    assert not any(x['ni'] for x in real)
    real = E.real_summary(planted(n_ni=8, rmse_ratio=1.049))
    assert all(x['ni'] for x in real)
    real = E.real_summary(planted(n_ni=8, rmse_ratio=1.06))
    assert not any(x['ni'] for x in real)
    assert E.real_summary(planted(), arch='A') == []


def test_verdict():
    recs = planted(1.0, 0.3)
    v = E.verdict(E.synth_summary(recs), E.real_summary(recs))
    assert v['go'] and v['go_syn'] and v['go_real'] and v['n_ni'] == 8
    recs = planted(1.0, 0.3, n_ni=5)
    v = E.verdict(E.synth_summary(recs), E.real_summary(recs))
    assert not v['go'] and v['go_syn'] and not v['go_real'] and v['n_ni'] == 5
    recs = planted(1.0, 0.3, n_ni=6)
    assert E.verdict(E.synth_summary(recs), E.real_summary(recs))['go']
    recs = planted(0.6, 0.55)
    v = E.verdict(E.synth_summary(recs), E.real_summary(recs))
    assert not v['go'] and not v['go_syn']
    # gain real but below the 0.10 floor
    recs = planted(0.8, 0.72, seed=4)
    v = E.verdict(E.synth_summary(recs), E.real_summary(recs))
    assert not v['go_syn']


def test_report_text_and_byte_identity(tmp_path):
    recs = planted(0.9, 0.4)
    meta = dict(stage='test', jobs=1)
    txt = E.report(meta, recs)
    assert txt == E.report(meta, recs)
    assert '**GO**' in txt and 'exact-edge share' in txt and 'non-inferior on 8 of 8' in txt
    assert '**KILL**' in E.report(meta, planted(0.5, 0.5, n_ni=0))
    raw = tmp_path / 'raw.npz'
    E.save_raw(str(raw), meta, recs)
    m2, r2 = E.load_raw(str(raw))
    assert m2 == meta and r2 == json.loads(json.dumps(recs))
    out1, out2 = tmp_path / 'a.md', tmp_path / 'b.md'
    for o in (out1, out2):
        subprocess.run([sys.executable, '-I', SCRIPT, '--report', '--from-raw', str(raw), '--md', str(o)], check=True)
    assert out1.read_bytes() == out2.read_bytes() == (txt + '\n').encode()
    # the harness form: no path, the output directory from the environment
    d = tmp_path / 'o'
    d.mkdir()
    E.save_raw(str(d / 'kan_raw.npz'), meta, recs)
    subprocess.run([sys.executable, '-I', SCRIPT, '--from-raw'], check=True, env=dict(os.environ, KAN_OUT=str(d)))
    assert (d / 'kan.md').read_bytes() == (txt + '\n').encode()


def ours_synth_rec(conf, counts, K=(1, 0, 2, 0, 0)):
    r = dict(kind='synth', tid=101, noise=0.1, seed=0, method='ours', conf=json.dumps(conf, sort_keys=True), truth=list(K), val_mse=0.1, test_rmse=0.1)
    r['counts_0.02'] = counts
    return r


def real_rec(method, conf, rm, arch='B', val=0.1, did=1):
    return dict(kind='real', did=did, name='x%d' % did, seed=0, arch=arch, method=method, conf=json.dumps(conf, sort_keys=True), val_mse=val, test_rmse=rm)


def test_choose_freeze():
    A, B, C = dict(kind='gd', lr=3.0), dict(kind='adam', lr=0.01), dict(kind='gd', lr=1.0, l2_ratio=0.1)
    K = [1, 0, 2, 0, 0]
    recs = [ours_synth_rec(A, K), ours_synth_rec(B, [1, 0, 2, 1, 0])]
    recs += [real_rec('ours', A, 0.5, 'H'), real_rec('ours', C, 0.4, 'H'), real_rec('ours', B, 0.01, 'B')]       # B (pure) is not eligible: wrong architecture
    for lr, v in ((0.003, 0.3), (0.01, 0.1), (0.03, 0.2)):
        recs.append(real_rec('ctrl', dict(lr=lr, lam=0.0), 0.3, 'B', val=v))
        recs.append(real_rec('fixed', dict(lr=lr, G=4), 0.35, 'B', val=0.1 / (1 + lr * 10)))
        recs.append(real_rec('mlp', dict(lr=lr, h=2), 0.32, 'B', val=0.1))
        for fam, conf in (('ctrl', dict(lr=lr, lam=0.0)), ('fixed', dict(lr=lr, G=4))):
            x = dict(kind='synth', tid=101, noise=0.1, seed=0, method=fam, conf=json.dumps(conf, sort_keys=True), val_mse=0.01 * (1 + (lr == 0.003)), test_rmse=0.1)
            recs.append(x)
        recs.append(dict(kind='synth', tid=101, noise=0.1, seed=0, method='mlp', conf=json.dumps(dict(lr=lr, h=2)), val_mse=0.01 * (1 + (lr == 0.01)), test_rmse=0.1))
    fz = E.choose_freeze(recs)
    assert fz['ours_synth'] == A                                       # the exact one
    assert fz['ours'] == C                                             # the lowest log ratio to the best baseline (0.4 / 0.3 against 0.5 / 0.3)
    assert list(fz['ours_scores']['real']) == sorted([json.dumps(A, sort_keys=True), json.dumps(C, sort_keys=True)])
    assert fz['ours_scores']['real'][json.dumps(A, sort_keys=True)] == pytest.approx(np.log(0.5 / 0.3))      # against the baselines only, never against ours
    assert fz['lr_real'] == {'ctrl': 0.01, 'fixed': 0.03, 'mlp': 0.003}
    assert fz['lr_synth'] == {'ctrl': 0.01, 'fixed': 0.01, 'mlp': 0.003}
    # without any real run of ours the real choice falls back to the synthetic one
    fz2 = E.choose_freeze([r for r in recs if not (r['method'] == 'ours' and r['kind'] == 'real')])
    assert fz2['ours'] == fz2['ours_synth'] == A


def test_choose_freeze_logs_validation():
    # the learning rate with the lower mean LOG validation MSE wins even if another has a lower arithmetic mean
    recs = []
    for did, (v1, v2) in enumerate(((1e-4, 1e-2), (1.0, 0.5)), start=1):
        recs.append(real_rec('ctrl', dict(lr=0.01, lam=0.0), 0.3, 'B', val=v1, did=did))
        recs.append(real_rec('ctrl', dict(lr=0.03, lam=0.0), 0.3, 'B', val=v2, did=did))
    recs.append(ours_synth_rec(dict(kind='gd', lr=1.0), [1, 0, 2, 0, 0]))
    fz = E.choose_freeze(recs)
    assert fz['lr_real']['ctrl'] == 0.01                               # geometric mean 0.01 against 0.0707; arithmetic mean 0.5 against 0.255


def test_report_dev_section():
    A, C = dict(kind='gd', lr=1.0), dict(kind='gd', lr=1.0, l2_ratio=0.1)
    recs = [ours_synth_rec(A, [1, 0, 2, 0, 0]), ours_synth_rec(dict(kind='adam', lr=0.01), [2, 1, 2, 0, 0]), real_rec('ours', C, 0.3, 'H')]
    for fam in ('ctrl', 'fixed', 'mlp'):
        recs.append(real_rec(fam, dict(lr=0.01, h=2), 0.3, 'B'))
        recs.append(dict(kind='synth', tid=101, noise=0.1, seed=0, method=fam, conf=json.dumps(dict(lr=0.01, h=2)), val_mse=0.01, test_rmse=0.1))
    txt = '\n'.join(E.dev_section(recs))
    assert '## Development stage' in txt and txt.count('chosen') == 2 and 'chosen (synthetic and [d, 1])' in txt
    assert txt.index('"gd"') < txt.index('chosen')
    assert 'synthetic {"ctrl": 0.01, "fixed": 0.01, "mlp": 0.01}; real {"ctrl": 0.01, "fixed": 0.01, "mlp": 0.01}' in txt
    assert 'Development stage' not in E.report(dict(stage='test'), planted(0.9, 0.4))
    assert 'Development stage' in E.report(dict(stage='dev'), planted(0.9, 0.4))


def test_registries():
    assert E.TAU == 0.02 and E.TAU in E.TAUS and len(E.TEST_TARGETS) == 10 and len(E.DEV_TARGETS) == 3 and not set(E.TEST_TARGETS) & set(E.DEV_TARGETS)
    assert E.NOISES == (0.1, 0.4) and E.SEEDS == (0, 1, 2) and E.NI_MARGIN == 1.05 and E.GO_REAL_SETS == 6 and E.GO_SHARE_GAIN == 0.10
    assert E.LAMS[0] == 0.0 and len(E.LAMS) == 6 and len(E.GS) == 6 and len(E.HS) == 8 and E.DEV_RANKS == (3, 6, 9)


# ---------------------------------------------------------------------------------------------------------------------------------
# torch parts (skipped without torch)
# ---------------------------------------------------------------------------------------------------------------------------------
def test_torch_models():
    pytest.importorskip('torch')
    torch, nn, KAN, MLP = E.build_torch()
    m = KAN([3, 1], 6, 'fact', eps=1e-4, seed=0)
    x = torch.randn(50, 3)
    assert float(m(x).abs().max().item()) < 1e-6                                    # small start: the readout is dormant
    assert float(m.layers[0].weights().abs().max()) == 0.0                  # alpha(0) = 0: exactly zero
    m.layers[0].alpha.data.fill_(1.0)
    w = m.layers[0].weights()
    assert torch.allclose(w, torch.full_like(w, 1e-4 * E.AMP_CAP * float(np.tanh(1.0 / E.AMP_CAP))), rtol=1e-4)
    f = KAN([3, 1], 5, 'fixed', seed=0)
    assert 'layers.0.c' not in dict(f.named_parameters()) and 'layers.0.ls' not in dict(f.named_parameters()) and 'layers.0.w' in dict(f.named_parameters())
    assert 'layers.0.c' in dict(KAN([3, 1], 5, 'free', seed=0).named_parameters())
    assert torch.allclose(f.layers[0].c[0, 0], torch.linspace(E.X_LO, E.X_HI, 5))
    fr = KAN([3, 4, 1], 4, 'free', seed=0)
    assert fr(x).shape == (50,) and len(fr.layers) == 2
    # active parameters of a planted fixed network
    xt = x
    f.layers[0].w.data.zero_(); f.layers[0].w.data[0, 0, 2] = 1.0
    s = E.kan_summary(f, xt, 1)
    assert s['n_0.02'] == 1 and s['P_0.02'] == 2 and s['counts_0.02'] == [1, 0, 0]
    # a trained 'free' network recovers one planted bump and the pruning masks the rest
    rng = np.random.default_rng(0)
    X = rng.uniform(-2.5, 2.5, (600, 3)).astype(np.float32)
    y = (np.exp(-0.5 * ((X[:, 0] - 1.0) / 0.4) ** 2)).astype(np.float32)
    data = {k: torch.tensor(v) for k, v in dict(Xtr=X[:400], ytr=y[:400], Xva=X[400:], yva=y[400:], Xte=X[400:], yte=y[400:]).items()}
    mm = KAN([3, 1], 10, 'free', seed=0)
    E.fit(mm, data, 'adam', 0.03, steps=800, l1=3e-3)
    before = E.kan_summary(mm, data['Xtr'], 3)['n_0.02']
    E.prune_and_finetune(mm, data, 0.03, data['Xtr'])
    assert E.rmse(mm, data['Xte'], data['yte']) < 0.15
    after = E.kan_summary(mm, data['Xtr'], 3)
    assert after['n_0.02'] <= before and float((mm.layers[0].mask == 0).sum()) >= 1       # pruning removed units, and the removed ones have zero weight
    assert float(mm.layers[0].weights()[mm.layers[0].mask == 0].abs().max()) == 0.0
    assert after['counts_0.02'][0] >= 1
    hy = KAN([3, 4, 1], 5, 'fact', eps=1e-3, hybrid=True)
    assert hy.layers[0].mode == 'fact' and hy.layers[1].mode == 'free' and hy(x).shape == (50,)
    assert KAN([3, 4, 1], 5, 'fact', eps=1e-3).layers[1].mode == 'fact'
    mlp = MLP(3, 4)
    assert sum(p.numel() for p in mlp.parameters()) == E.mlp_params(3, 4)


def test_centre_match():
    r = {'centres_0.02': [[0.1, 0.9], [], [5.0]], 'truth_c': [[0.0, 1.0], [], [1.0]]}
    assert E.centre_match(r) == (2 / 3, 2 / 3)                       # 2 of 3 true bumps found; 2 of 3 active units near a bump
    r = {'centres_0.02': [[0.25]], 'truth_c': [[0.0]]}
    assert E.centre_match(r) == (0.0, 0.0) and E.centre_match(r, tol=0.3) == (1.0, 1.0)
    assert E.centre_match({'centres_0.02': [[]], 'truth_c': [[]]}) == (1.0, 1.0)
    assert E.centre_match({'counts_0.02': [1]}) is None


def test_report_centre_table():
    txt = E.report(dict(stage='test'), planted(1.0, 0.4, centres=True))
    assert 'Do the active units sit on the true bumps' in txt
    assert 'Do the active units sit' not in E.report(dict(stage='test'), planted(1.0, 0.4))
    line = [l for l in txt.splitlines() if l.startswith('| ours | ') and l.count('|') == 4][0]
    assert line == '| ours | 1.000 | 1.000 |'


def test_real_summary_architectures():
    # ours on H is compared with the baselines of B; ours on B (the pure small start) never enters the primary summary; A is compared within A
    recs = planted(n_ni=8, rmse_ratio=1.0)
    extra = []
    for r in recs:
        if r['kind'] == 'real' and r['method'] == 'ours':
            extra.append(dict(r, arch='B', test_rmse=1.0, **{'P_%g' % t: 5 for t in E.TAUS}))                 # pure two-layer small start: stuck, worse
            extra.append(dict(r, arch='A', test_rmse=0.9, **{'P_%g' % t: 100 for t in E.TAUS}))
            extra.append(dict(r, arch='A', method='mlp', conf=json.dumps(dict(h=2, lr=0.01)), val_mse=0.1, test_rmse=0.5, **{'P_%g' % t: 100 for t in E.TAUS}))
    real = E.real_summary(recs + extra)
    assert all(x['ni'] for x in real) and [x['ours_mean'] for x in real] == pytest.approx([0.4] * 8)
    assert [x['ratio'] for x in E.real_summary(recs + extra, 'A')] == pytest.approx([0.9 / 0.5] * 8)
    assert not any(x['ni'] for x in E.real_summary(recs + extra, 'A'))
    assert 'Pure small start' in E.report(dict(stage='test'), recs + extra) and 'Pure small start' not in E.report(dict(stage='test'), recs)


# ---------------------------------------------------------------------------------------------------------------------------------
# tests added to kill mutants (tools/mutants/wp_kanspike_own.py)
# ---------------------------------------------------------------------------------------------------------------------------------
def test_registry_literals():
    assert (E.N_TRAIN, E.N_VAL, E.N_TEST_SYN) == (1500, 500, 2000) and E.N_BOOT == 10000 and E.REAL_CAP == 4000 and E.D_SYN == 5
    assert 'titanic' in E.EXCLUDE_REAL and 'kin8nm' in E.EXCLUDE_REAL and E.TAUS == (0.01, 0.02, 0.05) and E.M_UNITS == 10
    assert E.DEV_SEEDS == (0, 1) and E.TEST_TARGETS == tuple(range(1, 11))


def test_target_counts_pinned():
    expected = {1: [2, 0, 0, 4, 0], 2: [2, 0, 0, 3, 1], 3: [0, 1, 2, 2, 3], 4: [3, 0, 2, 1, 3], 5: [0, 4, 1, 2, 3], 6: [0, 1, 4, 2, 4], 7: [0, 3, 2, 0, 2],
                8: [3, 1, 0, 1, 1], 9: [4, 0, 3, 1, 1], 10: [3, 1, 2, 0, 2], 101: [0, 2, 0, 4, 1], 102: [2, 0, 3, 1, 1], 103: [1, 3, 2, 0, 3]}
    for tid, K in expected.items():
        assert E.make_target(tid)['K'].tolist() == K
    # amplitudes have both signs (45 of the 80 test bumps are negative)
    amps = np.concatenate([b[:, 0] for t in E.TEST_TARGETS for b in E.make_target(t)['bumps']])
    assert (amps < 0).sum() == 45 and (amps > 0).sum() == 35


def test_unit_contrib_values_off_centre():
    x = np.array([[0.0], [2.0]])
    w, c, s = np.array([[[3.0]]]), np.array([[[0.5]]]), np.array([[[2.0]]])
    U = E.unit_contrib(w, c, s, x)[:, 0, 0, 0]
    assert np.allclose(U, [3.0 * np.exp(-0.5 * (0.5 / 2.0) ** 2), 3.0 * np.exp(-0.5 * (1.5 / 2.0) ** 2)])
    x2 = np.array([[0.0], [4.0]])
    U2 = E.unit_contrib(np.array([[[1.0]]]), np.array([[[0.0]]]), np.array([[[1.0]]]), x2)[:, 0, 0, 0]
    assert np.allclose(U2, [1.0, np.exp(-8.0)])                            # u = 4: the square matters


def test_layer_activity_sd_over_the_node():
    rng = np.random.default_rng(3)
    x = rng.uniform(-2.5, 2.5, (500, 2))
    w = rng.standard_normal((2, 2, 2))
    c = rng.uniform(-2, 2, (2, 2, 2))
    s = np.full((2, 2, 2), 0.7)
    rms, sd, act = E.layer_activity(w, c, s, x)
    U = E.unit_contrib(w, c, s, x)
    for j in range(2):
        node = U[:, j].sum(axis=(1, 2))                                    # the node's pre-activation: sum over inputs and units
        assert sd[j] == pytest.approx(np.std(node))
        for i in range(2):
            for k in range(2):
                assert rms[j, i, k] == pytest.approx(np.sqrt(np.mean(U[:, j, i, k] ** 2)))
                assert bool(act[j, i, k]) == bool(rms[j, i, k] >= E.TAU * sd[j])


def test_boot_ci_width_matches_normal_theory():
    v = np.random.default_rng(11).standard_normal(400)
    m, lo, hi = E.boot_ci(v)
    se = np.std(v) / np.sqrt(len(v))
    assert (hi - lo) == pytest.approx(2 * 1.96 * se, rel=0.06)             # a 95% interval, not 90% or 99%
    assert m == pytest.approx(np.mean(v)) and lo < m < hi


def test_boot_ci_is_clustered_by_target():
    # the gain is constant within a target: the interval must be that of ten clusters (not of sixty runs)
    gains = np.array([0.5, 0.1, 0.4, 0.0, 0.6, 0.2, 0.3, 0.7, -0.1, 0.4])
    rows = []
    for t, g in enumerate(gains):
        for rep in range(6):
            rows.append({'key': (t, 0.1, rep), 'a': {'share': 0.5 + g}, 'b': {'share': 0.5}})
    m, lo, hi = E.paired_gain(rows, 'a', 'b')
    assert m == pytest.approx(gains.mean())
    assert (hi - lo) == pytest.approx(2 * 1.96 * gains.std() / np.sqrt(10), rel=0.15)


def test_one_se_pick_ties_prefer_lower_validation():
    cands = [dict(val=0.104, P=30), dict(val=0.100, P=30), dict(val=0.102, P=30)]
    assert E.one_se_pick(cands, 500)['val'] == 0.100


def test_select_real_symbolic_with_consistent_counts():
    lst = [(100, 'a', 2000, 8, 8, 1)]                                       # nn == nf but one symbolic feature
    assert E.select_real(lst, lambda d: (np.random.default_rng(d).standard_normal((200, 3)), np.random.default_rng(d).standard_normal(200))) == ([], [])


def synth_row(tid, noise, seed, truth, counts_ours, ctrl_variants, fixed_variants, ours_conf='{}'):
    """Records of one synthetic run: ours; ctrl (list of (lam, val, P, counts)); fixed (list of (G, val, test, counts))."""
    base = dict(kind='synth', tid=tid, noise=noise, seed=seed, arch='A', truth=truth, n_val=500, steps=1, secs=1.0, test_rmse=0.1)
    recs = []

    def add(method, conf, counts, val, P, test=0.1):
        r = dict(base, method=method, conf=json.dumps(conf, sort_keys=True), val_mse=val, test_rmse=test)
        for tau in E.TAUS:
            r['counts_%g' % tau], r['P_%g' % tau] = counts, P
        recs.append(r)
    add('ours', {}, counts_ours, 0.1, 20)
    for lam, val, P, counts in ctrl_variants:
        add('ctrl', dict(lam=lam), counts, val, P)
    for G, val, test, counts in fixed_variants:
        add('fixed', dict(G=G), counts, val, 30, test)
    return recs


def test_synth_summary_selection_and_counts():
    K = [1, 0, 2, 0, 0]
    recs = synth_row(1, 0.1, 0, K, [2, 0, 1, 0, 1],
                     ctrl_variants=[(0.0, 0.100, 50, [3, 1, 2, 0, 0]), (1e-3, 0.101, 40, [1, 0, 2, 1, 0]), (1e-2, 0.101, 30, [1, 0, 2, 0, 0]), (3e-2, 0.5, 5, [0, 0, 1, 0, 0])],
                     fixed_variants=[(4, 0.2, 0.5, [4, 4, 4, 4, 4]), (8, 0.1, 0.9, [8, 8, 8, 8, 8])])
    rows = E.synth_summary(recs)
    assert len(rows) == 1
    r = rows[0]
    assert r['ours']['over'] == 2 and r['ours']['under'] == 1 and r['ours']['share'] == pytest.approx(0.4)      # counts [2,0,1,0,1] against [1,0,2,0,0]
    assert r['ctrl']['conf'] == json.dumps(dict(lam=1e-2)) and r['ctrl']['share'] == 1.0                          # fewest parameters within one SE
    assert r['ctrl_or']['conf'] == json.dumps(dict(lam=1e-2)) and r['ctrl_or']['share'] == 1.0
    assert r['fixed']['conf'] == json.dumps(dict(G=8)) and r['fixed']['rmse'] == 0.9                               # chosen by validation, not by test
    # oracle ties go to the smaller network
    recs2 = synth_row(1, 0.1, 0, K, K, ctrl_variants=[(0.0, 0.1, 50, K), (1e-3, 0.1, 20, K), (1e-2, 0.9, 80, K)], fixed_variants=[(4, 0.2, 0.9, K)])
    assert E.synth_summary(recs2)[0]['ctrl_or']['conf'] == json.dumps(dict(lam=1e-3))
    # an incomplete run (no control or no fixed grid) is not scored
    assert E.synth_summary(recs + synth_row(2, 0.1, 0, K, K, [], [])) == rows
    only_ours = [r for r in synth_row(3, 0.1, 0, K, K, [(0.0, 0.1, 50, K)], []) if r['method'] != 'fixed']
    assert len(E.synth_summary(recs + only_ours)) == 1


def test_real_summary_geometric_mean_and_tau():
    # three splits with ratios 1.0, 1.0, 1.1 against baselines of different size: the mean LOG ratio decides (1.032 <= 1.05 although the ratio of means is 1.057)
    recs = []
    for s, (base_rm, ours_rm) in enumerate(((0.2, 0.2), (0.4, 0.4), (0.8, 0.88))):
        b = dict(kind='real', did=1, name='x', seed=s, arch='B', d=5, n_train=100, steps=1, secs=1.0)
        recs.append(dict(b, method='mlp', conf=json.dumps(dict(h=2, lr=0.01)), val_mse=0.1, test_rmse=base_rm, **{'P_%g' % t: 50 for t in E.TAUS}))
        recs.append(dict(b, arch='H', method='ours', conf='{}', val_mse=0.1, test_rmse=ours_rm, **{'P_0.01': 10, 'P_0.02': 100, 'P_0.05': 1000}))
    r = E.real_summary(recs)[0]
    assert r['ni'] and r['clamped'] == 0 and r['ratio'] == pytest.approx(np.exp(np.log(1.1) / 3))
    assert np.mean([x['ours'] for x in r['rows']]) / np.mean([x['base'] for x in r['rows']]) > E.NI_MARGIN
    # the parameters of ours are those at tau = 0.02: with 1000 (tau 0.05) or 10 (tau 0.01) the mlp (50) could not be matched
    assert r['P_mean'] == 100 and r['clamped'] == 0


def test_real_summary_exact_margin(monkeypatch):
    monkeypatch.setattr(E, 'NI_MARGIN', 2.0)
    recs = planted(n_ni=8, rmse_ratio=2.0)
    assert all(x['ni'] for x in E.real_summary(recs))


def fake_rows(gains_by_target, n=6, ctrl_share=0.3, oracle_extra=0.0):
    rows = []
    for t, g in enumerate(gains_by_target):
        for rep in range(n):
            rows.append({'key': (t, 0.1, rep), 'ours': {'share': ctrl_share + g}, 'ctrl': {'share': ctrl_share}, 'ctrl_or': {'share': ctrl_share + g - oracle_extra}})
    return rows


def ni_rows(k):
    return [dict(ni=i < k) for i in range(8)]


def test_verdict_synthetic_conditions():
    # lower bound above 0 but gain below the floor
    v = E.verdict(fake_rows([0.05] * 10), ni_rows(8))
    assert v['gain'][1] > 0 and v['gain'][0] < E.GO_SHARE_GAIN and not v['go_syn'] and not v['go']
    # gain above the floor but the interval contains 0
    v = E.verdict(fake_rows([0.6, -0.5, 0.6, -0.5, 0.6, -0.5, 0.6, 0.6, -0.5, 0.6]), ni_rows(8))
    assert v['gain'][0] >= E.GO_SHARE_GAIN and v['gain'][1] <= 0 and not v['go_syn']
    # both: go
    v = E.verdict(fake_rows([0.3] * 10), ni_rows(8))
    assert v['go_syn'] and v['go']
    # the verdict is against the validation-penalised control, not the oracle one
    v = E.verdict(fake_rows([0.3] * 10, oracle_extra=0.3), ni_rows(8))
    assert v['gain'][0] == pytest.approx(0.3) and v['go']
    assert not E.verdict(fake_rows([0.3] * 10), ni_rows(5))['go']


def test_report_counts_tau_and_stability():
    recs = planted(0.9, 0.4, n_ni=6)
    txt = E.report(dict(stage='test'), recs)
    assert 'non-inferior on 6 of 8 sets' in txt and '(tau = 0.02)' in txt
    rows = E.synth_summary(recs)
    cells, cc = {}, {}
    for r in rows:
        cells.setdefault(r['key'][:2], []).append(r['ours']['share'])
        cc.setdefault(r['key'][:2], []).append(r['ctrl']['share'])
    line = 'ours %s, control %s' % (E.fmt(np.mean([np.std(v) for v in cells.values()])), E.fmt(np.mean([np.std(v) for v in cc.values()])))
    assert line in txt.splitlines()


def test_synth_over_under_signs():
    K = [3, 0, 1, 0, 0]
    rows = E.synth_summary(synth_row(1, 0.1, 0, K, [1, 2, 1, 0, 0], [(0.0, 0.1, 5, K)], [(4, 0.1, 0.1, K)]))
    assert rows[0]['ours']['over'] == 2 and rows[0]['ours']['under'] == 2
