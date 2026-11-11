"""Tests of experiments/exp_noise_order.py: scoring of an empty result, a failed trial recorded not fatal, a tiny end-to-end run and its report."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_noise_order as X      # noqa: E402


def test_score_of_an_empty_result_is_order_zero_nothing_matched():
    r = X._score(0, np.zeros(0), np.zeros((0, 1)), np.array([0.1, 0.5]), np.array([1.0, -1.0]))
    assert r[0] == 0 and r[1] == 0 and np.isnan(r[2]) and np.isnan(r[3])
    r = X._score(2, np.array([1.0, -0.9]), np.array([[0.1], [0.52]]), np.array([0.1, 0.5]), np.array([1.0, -1.0]))
    assert r[0] == 2 and r[1] == 2 and r[2] < 0.03 and r[3] < 0.15


def test_trial_failure_is_recorded_as_nan_rows(monkeypatch):
    def boom(args):
        raise RuntimeError('planted')
    monkeypatch.setattr(X, '_trial', boom)
    out, secs, meta = X.trial((3, 1, 3001, 0))
    assert out.shape == (len(X.COLS), 4) and np.all(np.isnan(out)) and meta[0] == 3 and meta[2] == 3001


def test_tiny_run_report_excludes_failed_trials_and_is_deterministic(tmp_path, monkeypatch):
    monkeypatch.setattr(X, 'KS', (2,))
    monkeypatch.setattr(X, 'NOISE', (0.03,))
    path = str(tmp_path / 'raw.npz')
    X.run(1, 1, 1, out_path=path)
    z = dict(np.load(path))
    assert z['meta'][0, 0] == 2 and z['meta'][0, 2] == 2000 and z['res'].shape == (1, len(X.COLS), 4)      # K = 2, target seed 1000 K + 0
    res = z['res'].copy()
    res = np.concatenate([res, np.full((1, len(X.COLS), 4), np.nan)])           # one failed trial appended
    z.update(res=res, secs=np.concatenate([z['secs'], np.full((1, z['secs'].shape[1]), np.nan)]), meta=np.concatenate([z['meta'], np.array([[2, 0, 1, 0, np.nan, np.nan]])]))
    np.savez_compressed(path, **z)
    md = X.report(path)
    assert md == X.report(path) and 'Failed trials (an exception): 1' in md and '| sliding ISS, BIC |' in md
    assert '## Where the two orders differ' in md and '| identical entry sequences |' in md or '| same kernels, another order |' in md or '| different kernels |' in md


def test_order_strata_are_identical_for_one_feature_and_in_range_otherwise():
    pairs = [(1, 1001), (1, 1002), (3, 3001), (3, 3002), (4, 4001), (4, 4002)]
    st = X.order_strata(pairs)
    assert st[(1, 1001)] == 0 and st[(1, 1002)] == 0                         # one feature: one entry, the same kernel
    assert set(st) == set(pairs) and set(st.values()) <= {0, 1, 2}
    assert st == X.order_strata(pairs + pairs)                               # deterministic, duplicates ignored
    ct, at = X.make_target(np.random.default_rng(3001), 3)
    K, grid = X.setup()
    f = K.synth(ct[:, None], at)
    a = [int(j) for j in X.B.iss_path(K, f, grid, max_entries=3)['order'][:3]]
    b = [int(j) for j in X.B.omp_slide(K, f, grid, max_units=3)['order'][:3]]
    assert st[(3, 3001)] == (0 if a == b else (1 if sorted(a) == sorted(b) else 2))


def test_wilson_interval():
    lo, hi = X.wilson(5, 10)
    assert abs(lo - 0.2366) < 1e-3 and abs(hi - 0.7634) < 1e-3
    assert X.wilson(10, 10)[1] == 1.0 or abs(X.wilson(10, 10)[1] - 1.0) < 1e-12
    assert np.isnan(X.wilson(0, 0)[0])


def test_report_rates_biases_and_strata_on_a_known_raw(monkeypatch):
    """Eight trials (K = 1, 2; two noise levels; two draws): the exact-order rates, biases and strata are read off by hand."""
    monkeypatch.setattr(X, 'KS', (1, 2))
    monkeypatch.setattr(X, 'NOISE', (0.03, 0.1))
    cols = list(X.COLS)
    nc = len(cols)
    meta = np.array([[1, 0, 1001, 0, np.nan, 1.0], [1, 0, 1001, 1, np.nan, 1.0], [1, 1, 1001, 0, np.nan, 1.0], [1, 1, 1001, 1, np.nan, 1.0],
                     [2, 0, 2001, 0, 0.5, 0.5], [2, 0, 2001, 1, 0.5, 0.5], [2, 1, 2001, 0, 0.5, 0.5], [2, 1, 2001, 1, 0.5, 0.5]], float)
    res = np.zeros((8, nc, 4))
    res[:, :, 0] = meta[:, 0][:, None]                               # every method exact ...
    res[:, :, 1] = meta[:, 0][:, None]
    res[:, :, 2] = 0.01
    res[:, :, 3] = 0.02
    j = cols.index('iss|disc')
    res[[2, 3], j, 0] = 2                                            # ... except 'iss|disc' over-estimates the two K = 1, s = 0.1 trials
    res[[6], j, 0] = 1                                               # and under-estimates one K = 2, s = 0.1 trial
    secs = np.ones((8, 8)) * 0.5
    secs[:, 0] = [0.1, 0.2, 0.3, 0.4, 10.0, 0.5, 0.5, 0.5]            # median 0.45, mean 1.5
    z = dict(res=res, meta=meta, secs=secs, cols=np.array(cols), noise=np.array(X.NOISE), n_targets=1, n_draws=2)
    import tempfile
    d = tempfile.mkdtemp()
    path = os.path.join(d, 'raw.npz')
    np.savez_compressed(path, **z)
    md = X.report(path)
    lab = X.LABEL.get('iss|disc', 'iss|disc')
    ok = X.LABEL.get('omp|disc', 'omp|disc')
    lo, hi = X.wilson(5, 8)
    # by noise level: s = 0.03 -> 4 of 4 (100 %), s = 0.1 -> 1 of 4 (25 %), all -> 5 of 8
    lo4, hi4 = X.wilson(1, 4)
    assert '| %s | 100 %% [%.0f, 100] | %.0f %% [%.0f, %.0f] | %.0f %% [%.0f, %.0f] |' % (lab, 100 * X.wilson(4, 4)[0], 25, 100 * lo4, 100 * hi4, 100 * 5 / 8, 100 * lo, 100 * hi) in md
    assert '| %s | 100 %% [%.0f, 100] | 100 %% [%.0f, 100] | 100 %% [%.0f, 100] |' % (ok, 100 * X.wilson(4, 4)[0], 100 * X.wilson(4, 4)[0], 100 * X.wilson(8, 8)[0]) in md
    assert '| sliding ISS (limit path, 10 entries) | 0.450 |' in md and '| matching pursuit + sliding | 0.500 |' in md
    # by K: K = 1 -> 2 of 4 (50 %), K = 2 -> 3 of 4 (75 %)
    assert '| %s | 50 %% | 75 %% |' % lab in md
    # bias: K-hat - K = +1 on two trials, -1 on one: mean +1/8 = +0.12 (0.125 -> '+0.12'), 25 % over, 12 % under (12.5 -> 12)
    assert '| %s | +0.12 | 25 | 12 |' % lab in md
    assert '| %s | +0.00 | 0 | 0 |' % ok in md
    # clustered by target: 2 targets; 'iss|bic' is exact everywhere, 'iss|disc' on 2 of 4 trials of target 1 and 3 of 4 of target 2: per-target differences 0.5 and 0.25, mean 37.5 points, the interval of
    # the bootstrap over two targets 25 to 50; a column exact everywhere differs by 0; 'iss|bic' itself is not compared with itself
    assert '## Exact order: differences clustered by target' in md and 'over the 2 targets' in md
    assert '| %s | +37.5 [+25.0, +50.0] |' % lab in md and '| %s | +0.0 [+0.0, +0.0] |' % ok in md
    assert '| %s | +' % X.LABEL.get('iss|bic', 'iss|bic') not in md.split('## Exact order: differences clustered by target')[1].split('## Exact order, by the true number')[0]
    assert md.split('## Exact order: differences clustered by target')[1].split('## Exact order, by the true number')[0].count('\n| ') == len(cols) - 1 + 1       # the header row and one row per other column
    # strata: one row per occupied stratum with its number of targets and trials
    assert '## Where the two orders differ' in md
    st = X.order_strata([(1, 1001), (2, 2001)])
    for g in sorted(set(st.values())):
        n_t = sum(1 for v in st.values() if v == g)
        n_tr = 4 * n_t
        assert '| %s | %d | %d |' % (X.STRATA[g], n_t, n_tr) in md
    assert sum(md.count('| %s |' % name) for name in X.STRATA.values()) == len(set(st.values()))


def test_make_target_properties():
    for K in (1, 2, 4, 6):
        for seed in (1000 * K, 1000 * K + 1, 1000 * K + 2):
            c, a = X.make_target(np.random.default_rng(seed), K)
            assert len(c) == len(a) == K and (np.diff(c) > 0).all() and c.min() >= -0.85 and c.max() <= 0.85
            if K > 1:
                assert np.diff(c).min() >= X.MIN_SEP
            assert abs(np.abs(a).max() - 1.0) < 1e-12 and np.abs(a).min() >= 0.3 - 1e-12             # the largest has magnitude 1, the others in [0.3, 1]
    signs = np.concatenate([np.sign(X.make_target(np.random.default_rng(s), 4)[1]) for s in range(40)])
    assert (signs > 0).any() and (signs < 0).any()


def test_target_cluster_diffs_resamples_targets_not_trials():
    # 4 targets of 3 trials: the base column is exact on targets 0 and 1 (all trials), the other column on target 0 only; the per-target differences are 0, 1, 0, 0
    tid = np.repeat(np.arange(4), 3)
    exact = np.zeros((12, 3), bool)
    exact[:6, 0] = True
    exact[:3, 1] = True
    exact[:, 2] = True
    (d, lo, hi), (d2, lo2, hi2) = X.target_cluster_diffs(exact, tid, 0, [1, 2], B=4000, seed=1)
    assert abs(d - 0.25) < 1e-12 and lo == 0.0 and abs(hi - 0.75) < 1e-12                       # the mean of the target differences 0, 1, 0, 0; resamples of 4 targets: 0, .25, .5, .75, 1
    assert abs(d2 - (-0.5)) < 1e-12 and hi2 <= 0.0 + 1e-12 and lo2 >= -1.0 - 1e-12             # exact on 6 of 12 trials against all 12
    # deterministic, and the interval is that of the targets: trial-level resampling of the same data would give a narrower one
    assert X.target_cluster_diffs(exact, tid, 0, [1], B=4000, seed=1) == [(d, lo, hi)]
    rng = np.random.default_rng(0)
    trial_boot = np.array([(exact[ix, 0].mean() - exact[ix, 1].mean()) for ix in (rng.integers(0, 12, 12) for _ in range(4000))])      # paired resampling of the 12 trials
    assert np.percentile(trial_boot, 97.5) - np.percentile(trial_boot, 2.5) < hi - lo
    # a column compared with itself: no difference at all
    assert X.target_cluster_diffs(exact, tid, 0, [0], B=200, seed=0) == [(0.0, 0.0, 0.0)]


def test_target_ids_separate_targets_of_the_same_order_and_join_their_trials():
    ids = X.target_ids(np.array([2, 2, 2, 1]), np.array([5.0, 6.0, 5.0, 5.0]))
    assert ids.tolist() == [2000005, 2000006, 2000005, 1000005] and len(set(ids.tolist())) == 3


def test_target_cluster_diffs_uses_the_documented_resamples():
    rng = np.random.default_rng(3)
    tid = np.repeat(np.arange(12), 4)
    exact = rng.uniform(size=(48, 2)) < np.array([0.8, 0.6])
    (d, lo, hi), = X.target_cluster_diffs(exact, tid, 0, [1], B=500, seed=0)
    per = np.array([exact[tid == t, 0].mean() - exact[tid == t, 1].mean() for t in range(12)])
    boot = per[np.random.default_rng(0).integers(0, 12, size=(500, 12))].mean(axis=1)
    assert abs(d - per.mean()) < 1e-15 and abs(lo - np.percentile(boot, 2.5)) < 1e-15 and abs(hi - np.percentile(boot, 97.5)) < 1e-15
