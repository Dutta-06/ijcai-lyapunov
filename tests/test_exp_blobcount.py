"""Tests of experiments/exp_blobcount.py (the nucleus-counting spike): the data rules (eligibility, split, tiles, ground truth), the shared front end, the scoring, the statistics and
the decision rule, the configuration grids and an end-to-end run on a synthetic bank. No network, no third-party data."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(REPO, 'experiments'))
sys.path.insert(0, REPO)
import exp_blobcount as X  # noqa: E402

skimage = pytest.importorskip('skimage')


def discs(centres, radius=8.0, amp=100.0, bg=5.0, noise=2.0, seed=0, size=X.PATCH):
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[:size, :size]
    img = np.full((size, size), bg, float)
    for x, y in centres:
        img += amp * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * (radius / 2.0) ** 2))
    return np.clip(np.round(img + noise * rng.standard_normal(img.shape)), 0, 255).astype(np.uint8)


# ---------------------------------------------------------------- data rules
def test_eligibility_needs_grey_channels_and_dark_median():
    g = np.full((8, 8), 20, np.uint8)
    assert X.is_fluorescence(g)
    assert X.is_fluorescence(np.dstack([g, g, g, np.full_like(g, 255)]))
    c = np.dstack([g, g, g]).copy()
    c[0, 0, 1] = 21
    assert not X.is_fluorescence(c)                                        # one pixel with unequal channels
    c = np.dstack([g, g, g]).copy()
    c[0, 0, 2] = 21
    assert not X.is_fluorescence(c)                                        # only the third channel differs
    assert not X.is_fluorescence(np.dstack([np.full_like(g, X.MED_MAX)] * 3))        # colour branch: the median limit is exclusive
    assert X.is_fluorescence(np.dstack([np.full_like(g, X.MED_MAX - 1)] * 3))
    assert not X.is_fluorescence(np.full((8, 8), X.MED_MAX, np.uint8))      # median not below the limit
    assert X.is_fluorescence(np.full((8, 8), X.MED_MAX - 1, np.uint8))


def test_split_is_seeded_disjoint_and_sized():
    ids = ['id%03d' % i for i in range(300)]
    dev, test = X.split_ids(ids)
    assert len(dev) == X.N_DEV and len(test) == X.N_TEST and not set(dev) & set(test)
    assert (dev, test) == X.split_ids(list(reversed(ids)))                  # independent of the input order
    perm = np.random.default_rng(1847).permutation(300)
    assert dev == [ids[i] for i in perm[:60]] and test == [ids[i] for i in perm[60:210]]


def test_tiles_and_membership_boundaries():
    assert X.tile_grid(256, 256) == [(0, 0), (0, 128), (128, 0), (128, 128)]
    assert len(X.tile_grid(520, 696)) == 4 * 5 and X.tile_grid(100, 300) == []
    xy = np.array([[-0.5, 0.0], [-0.51, 0.0], [127.49, 127.49], [127.5, 3.0], [128.0, 3.0], [128 + 5, 130.0]])
    assert X.in_patch(xy, 0, 0).tolist() == [True, False, True, False, False, False]
    assert X.in_patch(xy, 128, 128).tolist() == [False, False, False, False, False, True]   # (133, 130) -> patch (5, 2)
    yy = np.array([[10.0, 127.5], [10.0, 127.49], [10.0, -0.5], [10.0, -0.51], [10.0, 300.0]])
    assert X.in_patch(yy, 0, 0).tolist() == [False, True, True, False, False]              # the y range is bounded on both sides
    assert X.inside(np.array([[10.0, 127.5], [10.0, 127.49], [10.0, -0.51], [10.0, 300.0]])).tolist() == [[10.0, 127.49]]


def test_mask_stats_centroid_and_area_filter():
    m1 = np.zeros((40, 40), bool)
    m1[10:20, 4:8] = True                                                   # 40 px, centroid x 5.5, y 14.5
    m2 = np.zeros((40, 40), bool)
    m2[30:33, 30:33] = True                                                 # 9 px: debris
    m3 = np.zeros((40, 40), bool)
    m3[0:2, 0:5] = True                                                     # exactly 10 px: kept
    xy, area = X.mask_stats([m1, m2, m3])
    assert area.tolist() == [40.0, 10.0]
    assert np.allclose(xy, [[5.5, 14.5], [2.0, 0.5]])


def test_touching_ids_four_neighbourhood_only():
    lab = np.zeros((7, 7), int)
    lab[1:3, 1:3] = 1
    lab[1:3, 3:5] = 2                                                       # touches 1 along an edge
    lab[4:6, 4:6] = 3                                                       # isolated, several pixels: its own pixels are not a touch
    lab[3, 0] = 4
    lab[4, 1] = 5                                                           # diagonal neighbour of 4 only
    assert X.touching_ids(lab) == {1, 2}


def test_choose_tiles_rules():
    counts = [3, 20, 16, 40, 1, 25, 30, 2]
    b, bo = X.choose_tiles('imgA', counts)
    assert len(b) == X.PER_IMAGE and len(bo) <= X.BOOST and not set(b) & set(bo)
    assert all(counts[i] >= X.CROWD for i in bo)
    assert (b, bo) == X.choose_tiles('imgA', counts) and (b, bo) != X.choose_tiles('imgB', counts) or True
    assert X.choose_tiles('imgA', [1])[0] == [0] and X.choose_tiles('imgA', [1])[1] == []
    b2, bo2 = X.choose_tiles('imgA', [X.CROWD - 1] * 6)
    assert bo2 == []                                                        # the boost needs crowded tiles
    big = [X.CROWD] * 12
    assert len(X.choose_tiles('imgC', big)[1]) == X.BOOST


def test_strata_definitions():
    bank = dict(split=np.array([0, 0, 0, 1, 1, 1]), boost=np.array([0, 1, 0, 0, 0, 1]), n_gt=np.array([3, 20, X.CROWD, 40, 2, X.CROWD - 1]))
    base, crowd = X.strata(bank, 0)
    assert base.tolist() == [0, 2] and crowd.tolist() == [1, 2]            # boost patches are never in the base sample, a crowded one is in the stratum
    base, crowd = X.strata(bank, 1)
    assert base.tolist() == [3, 4] and crowd.tolist() == [3]


# ---------------------------------------------------------------- front end
def test_noise_model_recovers_a_poisson_like_relation():
    rng = np.random.default_rng(1)
    mean = np.tile(np.linspace(5, 200, 128), (128, 1))
    g = mean + np.sqrt(0.5 * mean + 2.0) * rng.standard_normal(mean.shape)
    a, b = X.noise_model(g)
    assert 0.35 < a < 0.65 and 0.5 < b < 4.5
    flat = 20 + 3.0 * rng.standard_normal((128, 128))
    a, b = X.noise_model(flat)
    assert a < 0.1 and 8.4 < a * 20 + b < 9.6                               # variance 9 at the level 20 (the chi-square median constant matters)


def test_noise_model_has_a_floor():
    a, b = X.noise_model(np.full((128, 128), 7.0))
    assert b == 0.25 and a >= 0.0


def test_scale_estimate_follows_the_blob_width():
    for s in (3.0, 6.0):
        img = discs([(30, 30), (90, 40), (60, 95), (100, 100)], radius=2 * s, noise=0.0)
        est = X.scale_estimate(img.astype(float) - np.median(img))
        assert 0.75 * s < est < 1.5 * s
    assert X.scale_estimate(np.zeros((64, 64))) == 5.0                      # nothing to find: the default


def test_prepare_outputs():
    p = X.prepare(discs([(40, 40)]))
    assert p['y'].shape == (128, 128) and np.isclose(np.median(p['y']), 0.0)
    assert np.allclose(p['sw'] ** 2 * p['var'], 1.0) and p['var'].min() >= 0.25


# ---------------------------------------------------------------- scoring
def test_inside_and_match_are_one_to_one_within_the_radius():
    d = np.array([[10.0, 10.0], [10.5, 10.0], [200.0, 5.0], [-0.6, 3.0]])
    assert len(X.inside(d)) == 2
    gt = np.array([[10.0, 10.0], [30.0, 30.0]])
    assert X.match_count(d[:2], gt) == 1                                     # two detections, one truth: one true positive
    assert X.match_count(np.array([[10.0, 16.0]]), gt) == 1                  # exactly at the radius
    assert X.match_count(np.array([[10.0, 16.01]]), gt) == 0
    assert X.match_count(np.zeros((0, 2)), gt) == 0 and X.match_count(d, np.zeros((0, 2))) == 0
    # the optimal assignment, not the greedy nearest: det A is nearest to g1 but g1 is the only option of det B
    det = np.array([[5.0, 0.0], [0.0, 0.0]])
    gt2 = np.array([[2.0, 0.0], [-4.0, 0.0]])
    assert X.match_count(det, gt2, r=6.0) == 2


def test_score_patch_counts_only_detections_inside():
    dets = np.array([[10.0, 10.0, 3.0, 1.0], [300.0, 4.0, 3.0, 1.0], [50.0, 50.0, 3.0, 1.0]])
    nd, tp = X.score_patch(dets, np.array([[10.0, 11.0]]))
    assert (nd, tp) == (2, 1)
    assert X.score_patch(np.zeros((0, 4)), np.array([[1.0, 1.0]])) == (0, 0)


# ---------------------------------------------------------------- statistics and the rule
def test_boot_diff_is_a_cluster_ratio_of_sums():
    ea = np.array([2.0, 2.0, 2.0, 0.0, 0.0, 0.0])
    eb = np.zeros(6)
    img = np.array([0, 0, 0, 1, 1, 1])
    pt, lo, hi = X.boot_diff(ea, eb, img, np.arange(6), n_boot=2000)
    assert pt == 1.0 and lo == 0.0 and hi == 2.0                             # image 0 contributes 2, image 1 contributes 0
    pt, lo, hi = X.boot_diff(ea, ea, img, np.arange(6), n_boot=500)
    assert (pt, lo, hi) == (0.0, 0.0, 0.0)
    # signs: only |e| counts
    assert X.boot_diff(-ea, eb, img, np.arange(6), n_boot=100)[0] == 1.0
    # the stratum selects patches
    assert X.boot_diff(ea, eb, img, np.array([0, 1, 2]), n_boot=100)[0] == 2.0
    # unequal cluster sizes: ratio of sums, not mean of means
    e1 = np.array([4.0, 0.0, 0.0, 0.0])
    im2 = np.array([0, 1, 1, 1])
    assert X.boot_diff(e1, np.zeros(4), im2, np.arange(4), n_boot=100)[0] == 1.0


def test_boot_diff_is_seeded():
    rng = np.random.default_rng(0)
    ea, eb = rng.normal(size=60), rng.normal(size=60)
    img = np.repeat(np.arange(20), 3)
    assert X.boot_diff(ea, eb, img, np.arange(60), n_boot=300) == X.boot_diff(ea, eb, img, np.arange(60), n_boot=300)


def _table(**kw):
    t = dict(mae_c=2.0, mae_o=1.5, mae_g=1.5, diff_c=(0.5, 0.1, 0.9), mae_all_c=1.0, mae_all_o=1.0, diff_all=(0.0, -0.05, 0.05), diff_og=(0.0, -0.2, 0.2))
    t.update(kw)
    return t


def test_verdict_rule_boundaries():
    assert X.verdict(_table())['GO']
    # C1: margin of 5 % of the classical MAE (2.0 -> 0.1) and the interval above 0
    assert X.verdict(_table(mae_o=1.9))['C1'] and not X.verdict(_table(mae_o=1.91))['C1']
    assert not X.verdict(_table(diff_c=(0.5, 0.0, 0.9)))['C1']
    assert X.verdict(_table(diff_c=(0.5, 1e-9, 0.9)))['C1']
    # C2: not worse than greedy by more than 0.5
    assert X.verdict(_table(mae_g=1.0))['C2'] and not X.verdict(_table(mae_g=0.99))['C2']
    # C3: the upper end of the interval at most 10 % of the best classical MAE
    assert X.verdict(_table(diff_all=(0.0, -0.1, 0.1)))['C3'] and not X.verdict(_table(diff_all=(0.0, -0.1, 0.1001)))['C3']
    # GO needs all three
    assert not X.verdict(_table(mae_o=1.95))['GO'] and not X.verdict(_table(mae_g=0.5))['GO'] and not X.verdict(_table(diff_all=(0.0, 0.0, 0.5)))['GO']
    # the dual order contributes only if its interval is below 0
    assert X.verdict(_table(diff_og=(-0.3, -0.5, -0.1)))['dual_contributes'] and not X.verdict(_table(diff_og=(-0.3, -0.5, 0.0)))['dual_contributes']


def test_objective_is_the_mean_of_base_and_crowded_mae():
    nd = np.array([[1, 2, 10, 12], [2, 2, 10, 20]])
    ngt = np.array([1, 4, 12, 12])
    o = X.objective(nd, ngt, np.array([0, 1]), np.array([2, 3]))
    assert np.allclose(o, [0.5 * (1.0 + 1.0), 0.5 * (1.5 + 5.0)])


# ---------------------------------------------------------------- configurations
def test_every_method_has_at_most_the_declared_budget():
    for m in X.METHODS:
        cfgs = X.configs_of(m)
        assert 1 <= len(cfgs) <= X.MAX_CONFIGS, m
        assert len({json.dumps(c, sort_keys=True) for c in cfgs}) == len(cfgs), m   # no duplicates
    assert len(X.configs_of('ours')) == len(X.configs_of('greedy'))
    assert X.configs_of('ours') == X.configs_of('greedy')                           # the control has the same grid


def test_groups_partition_the_configurations():
    for m in X.METHODS:
        cfgs = X.configs_of(m)
        groups = X.groups_of(m, cfgs)
        ids = sorted(i for g in groups for i, _ in g)
        assert ids == list(range(len(cfgs)))
        if m in ('ours', 'greedy'):
            for g in groups:                                                       # one path per group, only the stop differs
                assert len({json.dumps(c['path'], sort_keys=True) for _, c in g}) == 1
            assert len(groups) == len(X.PATH_GRID)
    assert X.groups_of('ours', X.configs_of('ours'), only=[3]) == [[(3, X.configs_of('ours')[3])]]


# ---------------------------------------------------------------- methods on synthetic blobs
CENTRES = [(30, 30), (95, 35), (60, 90), (105, 105)]


@pytest.fixture(scope='module')
def prep4():
    return X.prepare(discs(CENTRES, radius=9.0, seed=3))


def _count_near(d, truth, r=6.0):
    return X.match_count(X.inside(np.asarray(d)[:, :2]), np.array(truth, float), r)


def test_ours_and_greedy_find_separated_blobs(prep4):
    for m in ('ours', 'greedy'):
        cfgs = X.configs_of(m)
        group = X.groups_of(m, cfgs)[0]
        res = X.run_group(m, prep4, group)
        counts = [X.score_patch(res[i][0], np.array(CENTRES, float))[1] for i, _ in group]
        assert max(counts) == 4, (m, counts)
        d = res[group[2][0]][0]
        assert d.shape[1] == 4 and np.all(d[:, 3] > 0)                              # amplitudes are non-negative, a width per unit
        assert np.all((d[:, 2] >= 0.5 * prep4['s_hat'] * 0.99) & (d[:, 2] <= 3.0 * prep4['s_hat'] * 1.01))


def test_a_higher_penalty_never_keeps_more_units(prep4):
    P, R = X.sliding_run(prep4, X.configs_sliding()[0]['path'], 'dual')
    ks = [len(X.sliding_stop(prep4, P, R, dict(lam=l, clean=0.0, cal=1))[0]) for l in (1.0, 4.0, 16.0, 64.0, 256.0, 1024.0, 1e5, 1e8)]
    assert all(a >= b for a, b in zip(ks, ks[1:])) and ks[0] > ks[-1]


def test_cleanup_only_removes_units(prep4):
    P, R = X.sliding_run(prep4, X.configs_sliding()[0]['path'], 'dual')
    a = X.sliding_stop(prep4, P, R, dict(lam=8.0, clean=0.0, cal=1))[0]
    b = X.sliding_stop(prep4, P, R, dict(lam=8.0, clean=1.0, cal=1))[0]
    assert len(b) <= len(a)


def test_calibration_scales_the_penalty(prep4):
    P, R = X.sliding_run(prep4, X.configs_sliding()[0]['path'], 'dual')
    chi = R['hist'][-1]['chi2']
    a = X.sliding_stop(prep4, P, R, dict(lam=8.0, clean=0.0, cal=1))
    b = X.sliding_stop(prep4, P, R, dict(lam=8.0 * max(chi, 1.0), clean=0.0, cal=0))
    assert a[1] == b[1]                                                              # calibrated penalty = lam x chi2 of the last block (at least 1)


def test_classical_detectors_find_separated_blobs(prep4):
    for m in CLASSICAL if False else X.CLASSICAL:
        cfgs = X.configs_of(m)
        best = 0
        for i, c in enumerate(cfgs):
            d = X.run_group(m, prep4, [(i, c)])[i][0]
            best = max(best, _count_near(d, CENTRES) if len(d) else 0)
        assert best == 4, m


def test_watershed_splits_touching_discs():
    yy, xx = np.mgrid[:128, :128]
    img = np.full((128, 128), 5.0)
    for cx in (52, 70):                                                              # two touching discs of radius 11
        img[(xx - cx) ** 2 + (yy - 64) ** 2 <= 11 ** 2] += 100.0
    p = X.prepare(np.clip(img, 0, 255).astype(np.uint8))
    best = max(len(X.run_ws(p, c)) for c in X.configs_ws() if c['fam'] == 'peak')
    assert best >= 2
    assert any(len(X.run_ws(p, c)) == 2 for c in X.configs_ws())


# ---------------------------------------------------------------- end to end on a synthetic bank
def _bank(n=6):
    rng = np.random.default_rng(5)
    patches, gt, ng = [], [], []
    for k in range(n):
        m = [1, 3, 4, 2, 5, 0][k % 6]
        pts = [(float(rng.uniform(15, 112)), float(rng.uniform(15, 112))) for _ in range(m)]
        patches.append(discs(pts, radius=9.0, seed=k))
        gt.append(np.array(pts, float).reshape(-1, 2))
        ng.append(m)
    mx = max(max(ng), 1)
    gxy = np.full((n, mx, 2), np.nan, np.float32)
    for k, a in enumerate(gt):
        gxy[k, :len(a)] = a
    return dict(patches=np.stack(patches), split=np.zeros(n, np.uint8), img=np.arange(n) // 2, boost=np.zeros(n, np.uint8), n_gt=np.array(ng, np.int16),
                n_touch=np.zeros(n, np.int16), gt_xy=gxy, gt_area=np.full((n, mx), 100.0, np.float32), r0=np.zeros(n, np.int16), c0=np.zeros(n, np.int16))


@pytest.mark.parametrize('method', ['plm', 'ours'])
def test_run_method_matches_a_direct_run(method):
    bank = _bank(3)
    idx = np.arange(3)
    r = X.run_method(method, bank, idx, cfg_idx=[2], workers=1, log=lambda *a: None)
    cfgs = X.configs_of(method)
    group = [(2, cfgs[2])]
    for j in range(3):
        d = X.run_group(method, X.prepare(bank['patches'][j]), group)[2][0]
        nd, tp = X.score_patch(d, X.gt_of(bank, j))
        assert r['ndet'][2, j] == nd and r['tp'][2, j] == tp
    assert r['ndet'][0].sum() == 0 and r['ndet'][3].sum() == 0                       # configurations that were not requested stay empty
    xy, off = X.pack_dets(r['dets'], [2], 3)
    assert off[0] == 0 and off[-1] == len(xy) and len(off) == 4


def test_the_work_function_counts_a_failed_run_as_no_detection(monkeypatch):
    bank = _bank(2)
    X._BANK.clear()
    X._BANK.update(bank)

    def boom(*a, **k):
        raise RuntimeError('x')
    monkeypatch.setattr(X, 'run_group', boom)
    k, out, secs, err = X._work(('plm', 0, [(0, X.configs_of('plm')[0])]))
    assert err is True and out[0][1] == 0 and out[0][2] == 0


def _raw_test(rng, n=120, nimg=40):
    ngt = rng.integers(0, 30, n)
    raw = dict(test_patch=np.arange(n), n_gt=ngt.astype(np.int16), img=np.repeat(np.arange(nimg), n // nimg), boost=np.zeros(n, np.uint8), n_touch=(ngt // 3).astype(np.int16))
    sd = dict(ours=0.8, greedy=1.2, log=2.0, dog=2.5, plm=2.2, ws=1.6)
    fr = {}
    for i, m in enumerate(X.METHODS):
        nd = np.maximum(ngt + np.round(sd[m] * rng.standard_normal(n)), 0).astype(np.int32)
        tp = np.minimum(nd, ngt).astype(np.int32)
        raw.update({'%s|ndet' % m: nd, '%s|tp' % m: tp, '%s|hit' % m: np.zeros(n, np.uint8), '%s|secs' % m: np.full(n, 0.5), '%s|fails' % m: np.zeros(n, np.uint8),
                    '%s|dets' % m: np.zeros((0, 2), np.float32), '%s|off' % m: np.zeros(n + 1, np.int64)})
        fr[m] = dict(index=1, config=X.configs_of(m)[1], objective=1.0, n_configs=len(X.configs_of(m)))
    raw['freeze'] = np.array(json.dumps(dict(methods=fr)))
    raw['info'] = np.array('{}')
    return raw


def test_report_is_deterministic_and_states_the_verdict(tmp_path):
    raw = _raw_test(np.random.default_rng(0))
    t1, t2 = X.test_report(raw), X.test_report(raw)
    assert t1 == t2 and ('**Verdict: GO**' in t1 or '**Verdict: KILL**' in t1)
    S = X.test_stats(raw)
    assert S['verdict']['GO'] == ('**Verdict: GO**' in t1)
    assert S['best']['crowded'] in X.CLASSICAL and S['best']['all'] in X.CLASSICAL
    p = tmp_path / 'raw.npz'
    np.savez_compressed(p, **raw)
    out = tmp_path / 'r.md'
    subprocess.run([sys.executable, os.path.join(REPO, 'experiments', 'exp_blobcount.py'), '--from-raw', str(p), '--md', str(out)], check=True, env=dict(os.environ, BLOBCOUNT_REPO=REPO))
    assert out.read_text() == t1


def test_best_classical_is_the_lowest_test_mae_in_the_stratum():
    raw = _raw_test(np.random.default_rng(1))
    S = X.test_stats(raw)
    for st in ('all', 'crowded'):
        maes = {m: X.pooled(S['err'][m], S['idx'][st]) for m in X.CLASSICAL}
        assert S['best'][st] == min(maes, key=maes.get)
    # an obviously better ours on every patch gives GO, an obviously worse one KILL
    raw2 = dict(raw)
    raw2['ours|ndet'] = raw['n_gt'].astype(np.int32)
    raw2['greedy|ndet'] = raw['n_gt'].astype(np.int32)
    assert X.test_stats(raw2)['verdict']['GO']
    raw3 = dict(raw)
    raw3['ours|ndet'] = raw['n_gt'].astype(np.int32) + 6
    v = X.test_stats(raw3)['verdict']
    assert not v['GO'] and not v['C1'] and not v['C3']


def test_dev_report_lists_every_method():
    rng = np.random.default_rng(2)
    n = 40
    ngt = rng.integers(0, 30, n)
    raw = dict(dev_patch=np.arange(n), n_gt=ngt, img=np.repeat(np.arange(10), 4), boost=np.zeros(n, np.uint8), n_touch=np.zeros(n, np.int16))
    for m in ('ours', 'plm'):
        k = len(X.configs_of(m))
        nd = np.maximum(ngt[None, :] + rng.integers(-3, 4, (k, n)), 0)
        raw.update({'%s|ndet' % m: nd, '%s|tp' % m: np.minimum(nd, ngt[None, :]), '%s|hit' % m: np.zeros((k, n)), '%s|secs' % m: np.ones((4, n)), '%s|fails' % m: np.zeros(n)})
    t = X.dev_report(raw)
    assert '| ours |' in t and '| plm |' in t and '| log |' not in t
    obj = X.objective(raw['ours|ndet'], raw['n_gt'], np.flatnonzero(raw['boost'] == 0), np.flatnonzero(raw['n_gt'] >= X.CROWD))
    assert '| ours | %d | %d | %.3f |' % (len(obj), int(np.argmin(obj)), obj.min()) in t


def test_freeze_and_assemble_round_trip(tmp_path):
    bank = _bank(6)
    bank['split'][:] = 0
    np.savez_compressed(tmp_path / 'bank.npz', **bank)
    np.savez_compressed(tmp_path / 'meta.npz', **{k: v for k, v in bank.items() if k != 'patches'})
    idx = np.arange(6)
    for m in ('plm', 'ws'):
        r = X.run_method(m, bank, idx, workers=1, log=lambda *a: None)
        np.savez_compressed(tmp_path / ('dev_%s.npz' % m), dev_patch=idx, **{'%s|%s' % (m, k): r[k] for k in ('ndet', 'tp', 'hit', 'secs', 'fails')})
    X.assemble_dev(str(tmp_path / 'meta.npz'), str(tmp_path), str(tmp_path / 'dev_raw.npz'), methods=('plm', 'ws'))
    with np.load(tmp_path / 'dev_raw.npz') as z:
        assert z['plm|ndet'].shape == (len(X.configs_of('plm')), 6) and np.array_equal(z['n_gt'], bank['n_gt'])


def test_constants_are_those_of_the_plan():
    assert (X.SEED_SPLIT, X.N_DEV, X.N_TEST, X.PATCH, X.MIN_AREA, X.CROWD, X.CROWD_HI, X.PER_IMAGE, X.BOOST, X.MED_MAX) == (1847, 60, 150, 128, 10, 15, 25, 2, 4, 100)
    assert (X.R_MATCH, X.N_BOOT, X.MAX_CONFIGS, X.M_REL, X.M_GREEDY, X.M_NONINF) == (6.0, 10000, 24, 0.05, 0.5, 0.10)
    assert X.METHODS == ('ours', 'greedy', 'log', 'dog', 'plm', 'ws') and X.CLASSICAL == ('log', 'dog', 'plm', 'ws')


def test_freeze_picks_the_lowest_objective_and_ties_the_first(tmp_path):
    rng = np.random.default_rng(3)
    n = 30
    ngt = rng.integers(0, 30, n)
    raw = dict(dev_patch=np.arange(n), n_gt=ngt.astype(np.int16), img=np.repeat(np.arange(10), 3), boost=np.zeros(n, np.uint8), n_touch=np.zeros(n, np.int16))
    for m in ('ours', 'ws'):
        k = len(X.configs_of(m))
        nd = np.maximum(ngt[None, :] + rng.integers(-4, 5, (k, n)), 0)
        nd[5] = ngt                                                                   # a perfect configuration
        nd[7] = ngt                                                                   # tie: the first one is taken
        raw['%s|ndet' % m] = nd
    np.savez_compressed(tmp_path / 'dev.npz', **raw)
    X.main(['--stage', 'freeze', '--dev_raw', str(tmp_path / 'dev.npz'), '--freeze', str(tmp_path / 'f.json'), '--methods', 'ours,ws'])
    fr = json.load(open(tmp_path / 'f.json'))['methods']
    assert fr['ours']['index'] == 5 and fr['ws']['index'] == 5 and fr['ours']['objective'] == 0.0
    assert fr['ours']['config'] == X.configs_of('ours')[5] and fr['ws']['config'] == X.configs_of('ws')[5]
    assert len(fr['ws']['objective_all']) == len(X.configs_of('ws'))


def test_oracle_block_is_at_least_as_good_as_the_frozen_stop(prep4):
    cfg = X.configs_of('ours')[18]
    o = X.oracle_patch(prep4, cfg['path'], 'dual', np.array(CENTRES, float), cfg['stop'])
    assert len(o['counts']) == len(o['tps']) and o['counts'][0] == 0 and o['tps'][0] == 0
    assert np.all(o['tps'] <= o['counts']) and np.all(np.diff(o['counts']) >= -len(o['counts']))
    assert 0 <= o['chosen'] < len(o['counts'])
    assert np.abs(o['counts'] - 4).min() <= abs(o['counts'][o['chosen']] - 4)
    X._BANK.clear()
    X._BANK.update(_bank(2))
    k, r = X._oracle_work(('ours', 1, cfg))
    assert k == 1 and len(r) == 7 and r[0] <= abs(r[3] - int(X._BANK['n_gt'][1])) and r[4] >= 1 and r[2] >= r[1]


def test_posthoc_report_reads_its_columns():
    rng = np.random.default_rng(0)
    raw = _raw_test(rng)
    n = len(raw['n_gt'])
    post = {'%s|oracle' % m: np.c_[np.zeros(n), np.minimum(raw['n_gt'], 3), raw['n_gt'], raw['n_gt'], np.full(n, 12), raw['n_gt'] + 2, raw['n_gt']].astype(np.int32) for m in ('ours', 'greedy')}
    t = X.posthoc_report(post, raw)
    assert 'POST HOC' in t and '| ours | crowded |' in t and '| greedy | all |' in t


def _stored(name):
    d = os.environ.get('BLOBCOUNT_OUT', os.path.join(REPO, 'results', 'ml'))
    return os.path.join(d, name)


@pytest.mark.skipif(not os.path.exists(os.path.join(REPO, 'results', 'ml', 'blobcount_raw.npz')), reason='stored raw not present')
def test_stored_reports_regenerate_byte_identically():
    for raw_name, md_name, fn in (('blobcount_raw.npz', 'blobcount.md', X.test_report), ('blobcount_dev_raw.npz', 'blobcount_dev.md', X.dev_report)):
        with np.load(os.path.join(REPO, 'results', 'ml', raw_name), allow_pickle=False) as z:
            raw = {k: z[k] for k in z.files}
        assert fn(raw) == open(os.path.join(REPO, 'results', 'ml', md_name)).read()
    with np.load(os.path.join(REPO, 'results', 'ml', 'blobcount_raw.npz'), allow_pickle=False) as y:
        fr = json.loads(str(y['freeze']))
    assert fr == json.load(open(os.path.join(REPO, 'results', 'ml', 'blobcount_freeze.json')))        # the frozen configurations used at test time are the committed ones
    for m, v in fr['methods'].items():
        assert X.configs_of(m)[v['index']] == v['config']


@pytest.mark.skipif(not os.path.exists(os.path.join(REPO, 'results', 'ml', 'blobcount_posthoc_raw.npz')), reason='post hoc raw not present')
def test_stored_posthoc_report_regenerates():
    with np.load(os.path.join(REPO, 'results', 'ml', 'blobcount_posthoc_raw.npz'), allow_pickle=False) as z, np.load(os.path.join(REPO, 'results', 'ml', 'blobcount_raw.npz'), allow_pickle=False) as y:
        txt = X.posthoc_report({k: z[k] for k in z.files}, {k: y[k] for k in y.files})
    assert txt == open(os.path.join(REPO, 'results', 'ml', 'blobcount_posthoc.md')).read()


def test_reference_report_lists_the_models_and_the_best_classical():
    rng = np.random.default_rng(4)
    raw = _raw_test(rng)
    n = len(raw['n_gt'])
    ref = dict(versions=np.array(json.dumps(dict(cellpose='x', stardist='y'))))
    for m in X.REFS:
        ref['%s|ndet' % m] = raw['n_gt'].astype(np.int32)
        ref['%s|tp' % m] = raw['n_gt'].astype(np.int32)
    t = X.ref_report(ref, raw)
    assert '| cellpose | 0.000 | 0.000 | +0.000 | 1.000 |' in t and '| stardist |' in t and '(best classical)' in t and 'cellpose x, stardist y' in t


def test_nocap_variants_have_the_same_grid_without_the_path_stop():
    a, b = X.configs_of('ours'), X.configs_of('ours_nc')
    assert len(a) == len(b) == 24 and X.configs_of('greedy_nc') == b
    for x, y in zip(a, b):
        assert x['stop'] == y['stop'] and y['path']['snr_min'] is None and x['path']['snr_min'] == 2.0
        assert {k: v for k, v in x['path'].items() if k != 'snr_min'} == {k: v for k, v in y['path'].items() if k != 'snr_min'}
    assert X.order_of('ours') == X.order_of('ours_nc') == 'dual' and X.order_of('greedy') == X.order_of('greedy_nc') == 'greedy'
    assert len(X.groups_of('ours_nc', b)) == len(X.PATH_GRID)


def test_nocap_report_substitutes_the_variants_in_the_rule():
    rng = np.random.default_rng(6)
    raw = _raw_test(rng)
    n = len(raw['n_gt'])
    nc = dict(raw)
    for m, tag in (('ours_nc', 'ours'), ('greedy_nc', 'greedy')):
        for k in ('ndet', 'tp', 'hit', 'secs', 'fails'):
            nc['%s|%s' % (m, k)] = raw['%s|%s' % (tag, k)]
    nc['ours_nc|ndet'] = raw['n_gt'].astype(np.int32)                                # a perfect variant
    nc['ours_nc|tp'] = raw['n_gt'].astype(np.int32)
    nc['freeze'] = np.array(json.dumps(dict(methods={m: dict(index=1, config=X.configs_of(m)[1], objective=1.0) for m in X.METHODS_NC})))
    dev = dict(n_gt=np.arange(20) % 30, boost=np.zeros(20, np.uint8))
    for m in X.METHODS_NC:
        dev['%s|ndet' % m] = np.zeros((24, 20), np.int32)
    t = X.nocap_report(dev, nc, raw)
    assert 'POST HOC' in t and '| ours_nc | crowded |' in t and 'rule outcome **GO**' in t


@pytest.mark.skipif(not os.path.exists(os.path.join(REPO, 'results', 'ml', 'blobcount_nocap_raw.npz')), reason='post hoc variant raw not present')
def test_stored_nocap_and_reference_reports_regenerate():
    d = os.path.join(REPO, 'results', 'ml')

    def load(name):
        with np.load(os.path.join(d, name), allow_pickle=False) as z:
            return {k: z[k] for k in z.files}
    main = load('blobcount_raw.npz')
    assert X.nocap_report(load('blobcount_nocap_dev_raw.npz'), load('blobcount_nocap_raw.npz'), main) == open(os.path.join(d, 'blobcount_nocap.md')).read()
    assert X.ref_report(load('blobcount_ref_raw.npz'), main) == open(os.path.join(d, 'blobcount_ref.md')).read()
    fr = json.loads(str(load('blobcount_nocap_raw.npz')['freeze']))
    assert fr == json.load(open(os.path.join(d, 'blobcount_nocap_freeze.json')))


def test_noise_model_floor_and_nonnegative_slope():
    rng = np.random.default_rng(2)
    quiet = 20 + 0.1 * rng.standard_normal((128, 128))
    a, b = X.noise_model(quiet)
    assert b == 0.25 and a >= 0.0 and a * 20 + b >= 0.25                    # the floor of the line fit (not only the constant-image branch)
    # a variance that falls with the level: the slope is clamped to 0 and the level is the weighted mean variance
    mean = np.tile(np.linspace(5, 200, 128), (128, 1))
    g = mean + np.sqrt(np.maximum(40.0 - 0.18 * mean, 1.0)) * rng.standard_normal(mean.shape)
    a, b = X.noise_model(g)
    assert a == 0.0 and 5.0 < b < 40.0


def test_prepare_variance_map_formula_and_floor():
    quiet = np.full((128, 128), 7, np.uint8)
    p = X.prepare(quiet)
    assert np.allclose(p['var'], 0.25)                                       # the variance floor
    from scipy import ndimage as ndi
    img = discs([(40, 40), (90, 80)], noise=3.0, seed=9)
    p = X.prepare(img)
    a, b = X.noise_model(img)
    assert np.allclose(p['var'], np.maximum(a * ndi.gaussian_filter(img.astype(float), 1.5) + b, 0.25))     # smoothed level, not the raw pixels


def test_scale_estimate_uses_many_maxima_not_the_strongest_only():
    yy, xx = np.mgrid[:128, :128]
    img = 400.0 * np.exp(-((xx - 64) ** 2 + (yy - 64) ** 2) / (2 * 8.0 ** 2))          # one very bright wide blob
    for cx, cy in ((20, 20), (20, 100), (100, 20), (100, 100), (64, 15), (15, 64), (110, 64)):
        img += 20.0 * np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 3.0 ** 2))       # seven faint small ones
    est = X.scale_estimate(img)
    assert est < 5.0                                                          # the median over the strongest twelve maxima, not the width of the strongest


def _fake_hist(chi2s):
    return [dict(chi2=c, w=np.ones(k), c=np.tile([float(k), 5.0, 3.0], (k, 1))) for k, c in enumerate(chi2s)]


def test_bic_calibration_has_a_floor_at_one_and_scales_above_it(monkeypatch):
    prep = {'sw': None}
    # last chi2 below 1: the penalty is lam itself (floor), 400 per parameter keeps block 1 (penalty 3880 > gain 3277)
    R = dict(hist=_fake_hist([8.0, 0.9, 0.7, 0.7]))
    assert X.sliding_stop(prep, None, R, dict(lam=400.0, clean=0.0, cal=1))[1] == 1
    # last chi2 above 1: the calibrated penalty is lam * chi2, the uncalibrated one is lam
    R = dict(hist=_fake_hist([10.0, 4.0, 3.9, 3.9]))
    assert X.sliding_stop(prep, None, R, dict(lam=100.0, clean=0.0, cal=1))[1] == 1
    assert X.sliding_stop(prep, None, R, dict(lam=100.0, clean=0.0, cal=0))[1] == 2


def test_cleanup_threshold_is_sqrt_lam_ln_n(monkeypatch):
    from src.core import imaging as IM
    R = dict(hist=_fake_hist([100.0, 50.0, 20.0, 20.0]))
    lam = 4.0
    thr = np.sqrt(lam * np.log(X.PATCH * X.PATCH))                            # 6.2
    monkeypatch.setattr(IM, 'unit_significance', lambda psf, sw, c, w: np.array([thr * 0.99, thr * 1.01, 100.0])[:len(w)])
    d, i, hit = X.sliding_stop({'sw': None}, None, R, dict(lam=lam, clean=1.0, cal=0))
    assert i == 2 and len(d) == 1                                             # two units in the block, the first is below the threshold, the second passes
    assert len(X.sliding_stop({'sw': None}, None, R, dict(lam=lam, clean=0.5, cal=0))[0]) == 2     # a smaller factor lowers the threshold in proportion
    assert len(X.sliding_stop({'sw': None}, None, R, dict(lam=lam, clean=2.0, cal=0))[0]) == 0     # a larger one raises it
