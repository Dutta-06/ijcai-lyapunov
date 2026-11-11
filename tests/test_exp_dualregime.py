"""Tests of experiments/exp_dualregime.py: the one-dimensional engine (kernels, block fit, the two orders of the sliding path, the exit rule, the resolution limit), the stopping rules,
the backward step, the scenes (geometry, significance scaling), the sliding Frank-Wolfe comparator, the statistics (paired bootstrap, verdicts, Clopper-Pearson), the declared grids
and seeds, part T, the two-dimensional arm's bookkeeping, the raw files, and the report (planted raws, and the stored raws when present).  Fast: single scenes and tiny grids."""
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_dualregime as X      # noqa: E402

NOISE0 = 0.188279       # the noise of the planted pair scene below


# ---- kernels and the block fit -----------------------------------------------------------------------------------------------------------------------------
def test_constants_are_the_declared_ones():
    assert (X.DX, X.NS, X.H, X.MARGIN) == (0.5, 96, 0.25, 2.0)
    assert np.isclose(X.LEN, 47.5, atol=0) and len(X.LAT) == 207 and np.isclose(X.LAT[0], -2.0, atol=0) and np.isclose(X.LAT[-1], 49.5, atol=0)
    assert np.isclose(X.XS[1] - X.XS[0], 0.5, atol=0) and len(X.XS) == 96
    assert X.TOL == 0.7 and X.SAME_TOL == 0.05
    assert np.isclose(X.NORM2, 3.5449077018110318, atol=1e-12) and np.isclose(X.NORM2, np.sum(np.exp(-(X.XS - X.XS[48]) ** 2)), atol=1e-12)
    assert np.isclose(X.PEN_BIC, 2 * np.log(96), atol=0) and X.PEN_HI == 16.0 and X.KAPPA_NAT == 4.5 and (X.PRUNE_EXTRA, X.EXTRA_ENTRIES) == (1, 3)
    assert X.METHODS == ('dual', 'greedy', 'sfw') and X.RULES == ('bic', 'pen', 'disc', 'nat', 'bicbw', 'penbw', 'oracleK', 'path')
    assert X.WIDTHS == (0.4, 3.0) and (X.BOOT, X.BOOT_SEED) == (2000, 0)


def test_gauss_values_shape_and_broadcast():
    G = X.gauss([10.0, 20.0], 1.0)
    assert G.shape == (96, 2)
    assert np.isclose(G[20, 0], 1.0, atol=0) and np.isclose(G[21, 0], np.exp(-0.125), atol=1e-15) and np.isclose(G[22, 0], np.exp(-0.5), atol=1e-15)
    G2 = X.gauss([10.0, 20.0], [1.0, 2.0])
    assert np.isclose(G2[22, 1], np.exp(-0.5 * (2.0 / 2.0) ** 2) * 0 + np.exp(-0.5 * ((11.0 - 20.0) / 2.0) ** 2), atol=1e-15)
    assert np.array_equal(X.gauss(5.0)[:, 0], X.gauss([5.0], 1.0)[:, 0])
    assert np.array_equal(X._PHI_LAT, X.gauss(X.LAT)) and X._PHI_LAT.shape == (96, 207)


def test_feat_norm_is_the_kernel_norm():
    assert np.isclose(X.feat_norm(1.0) ** 2, X.NORM2, atol=1e-12)
    assert np.isclose(X.feat_norm(2.0) ** 2, np.sum(np.exp(-(X.XS - X.XS[48]) ** 2 / 4.0)), atol=1e-12)
    assert X.feat_norm(2.0) > X.feat_norm(1.0)


def test_block_fit_recovers_a_noise_free_atom_and_a_width():
    y = (X.gauss([20.1]) * 2.5)[:, 0]
    c, s, w = X.block_fit(y, np.array([20.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), None)
    assert np.isclose(c[0], 20.1, atol=1e-7) and np.isclose(w[0], 2.5, atol=1e-7) and np.isclose(s[0], 1.0, atol=0)
    y2 = (X.gauss([20.1], [1.4]) * 2.0)[:, 0]
    c, s, w = X.block_fit(y2, np.array([20.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), (0.4, 3.0))
    assert np.allclose([c[0], s[0], w[0]], [20.1, 1.4, 2.0], atol=1e-6)
    c, s, w = X.block_fit(y2, np.array([20.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), (0.4, 1.2))
    assert s[0] <= 1.2 + 1e-12 and np.isclose(s[0], 1.2, atol=1e-6)                      # the width bound is active


def test_block_fit_box_and_nonnegativity():
    y = (X.gauss([20.0]) * 2.0)[:, 0]
    c, s, w = X.block_fit(y, np.array([20.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), None, box=1.5)
    assert np.isclose(c[0], 20.0, atol=1e-7)
    far = (X.gauss([24.0]) * 2.0)[:, 0]
    c, s, w = X.block_fit(far, np.array([20.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), None, box=1.5)
    assert c[0] <= 21.5 + 1e-9                                                          # cannot leave the box around the start
    c, s, w = X.block_fit(-y, np.array([20.0]), np.array([1.0]), np.array([0.5]), np.array([20.0]), None)
    assert (w >= 0).all() and np.isclose(w[0], 0.0, atol=1e-9)                          # a negative signal gets weight 0
    c, s, w = X.block_fit(y, np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), None)
    assert len(c) == len(s) == len(w) == 0


# ---- the path -----------------------------------------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope='module')
def pair_scene():
    rng = np.random.default_rng(100002)
    c, a, s, noise = X.scene(rng, dict(kind='pair', d=1.6, ratio=1, snr=10))
    y = X.gauss(c, s) @ a + noise * rng.standard_normal(X.NS)
    return c, a, noise, y


def test_single_atom_both_orders_enter_the_same_lattice_point():
    y = (X.gauss([20.1]) * 2.5)[:, 0]
    for o in ('dual', 'greedy'):
        R = X.path(y, 1e-3, o, max_entries=2)
        assert R['entered'][0] == int(np.argmin(np.abs(X.LAT - 20.0))) and len(R['blocks']) >= 2
        assert np.isclose(R['blocks'][1]['c'][0], 20.1, atol=1e-6) and np.isclose(R['blocks'][1]['w'][0], 2.5, atol=1e-6)
        assert R['blocks'][0]['c'].size == 0 and np.isclose(R['blocks'][0]['rss'], y @ y, atol=0)
    assert np.isclose(X.path(y, 1e-3, 'dual', max_entries=1)['taus'][0], 0.11312036, atol=1e-7)
    assert np.array_equal(X.path(y, 1e-3, 'greedy', max_entries=1)['taus'], [1.0])


def test_greedy_taus_are_the_step_counts_and_dual_taus_do_not_decrease(pair_scene):
    c, a, noise, y = pair_scene
    G, D = X.path(y, noise, 'greedy', max_entries=5), X.path(y, noise, 'dual', max_entries=5)
    assert np.array_equal(G['taus'], np.arange(1.0, len(G['taus']) + 1))
    assert (np.diff(D['taus']) >= 0).all() and D['taus'][0] > 0
    assert len(G['snr']) == len(G['taus']) and len(D['snr']) == len(D['taus'])


def test_planted_case_the_orders_differ_and_the_dual_keeps_the_second_feature(pair_scene):
    """Seed 100002 of the pair d = 1.6, ratio 1, significance 10: after the first atom the greedy order enters a noise peak (lattice point 12.25) before the second feature, the dual
    order enters the second feature first: the memory of the early correlation."""
    c, a, noise, y = pair_scene
    assert np.allclose(c, [38.53777883, 36.93777883], atol=1e-8) and np.isclose(noise, NOISE0, atol=1e-6)
    D, G = X.path(y, noise, 'dual', max_entries=5), X.path(y, noise, 'greedy', max_entries=5)
    pos = lambda R: [float(X.LAT[j]) for j in R['entered']]      # noqa: E731
    assert pos(D) == [37.75, 39.25, 12.25, 42.5, 8.75] and pos(G) == [37.75, 12.25, 40.0, 42.5, 8.75]
    assert np.allclose(D['snr'][:3], [17.64, 1.92, 2.94], atol=0.01) and np.allclose(G['snr'][:3], [17.64, 2.94, 2.46], atol=0.01)
    bd, bg = D['blocks'][2], G['blocks'][2]
    assert np.allclose(bd['c'], [37.179, 38.762], atol=1e-3) and np.allclose(bd['w'], [1.274, 0.762], atol=1e-3)
    assert np.allclose(bg['c'], [37.694, 12.2], atol=0.01)
    assert X.score(bd, c)[0] == 1 and X.score(bg, c)[0] == 0
    assert X.select_pen(D['blocks'], noise, X.PEN_BIC) == 2 and X.select_pen(G['blocks'], noise, X.PEN_BIC) == 3 and X.select_pen(D['blocks'], noise, 16.0) == 1


def test_resolution_limit_blocks_lattice_points_near_an_active_centre():
    """A wide Gaussian (width 2.5) is not representable by one kernel of width 1: the residual has wings.  The dual order without a resolution limit enters the ghosts next to the
    first atom (their duals have accumulated to almost the wall); the limit r_ex keeps every entrant at least r_ex from the active centre at its entry (lattice start); the greedy
    order does not enter the ghosts even without the limit."""
    y = (X.gauss([20.0], [2.5]) * 3.0)[:, 0]
    pos = lambda R: [float(X.LAT[j]) for j in R['entered']]      # noqa: E731
    assert pos(X.path(y, 1e-3, 'dual', r_ex=0.0, max_entries=3)) == [20.0, 19.75, 20.25]
    assert pos(X.path(y, 1e-3, 'dual', r_ex=1.0, max_entries=3)) == [20.0, 19.0, 19.75]
    assert pos(X.path(y, 1e-3, 'dual', r_ex=2.0, max_entries=3)) == [20.0, 18.0, 23.5]
    assert pos(X.path(y, 1e-3, 'greedy', r_ex=0.0, max_entries=3)) == [20.0, 17.25, 24.0]
    for r_ex in (1.0, 2.0, 3.0):
        R = X.path(y, 1e-3, 'dual', r_ex=r_ex, max_entries=2)
        assert abs(X.LAT[R['entered'][1]] - R['blocks'][1]['c'][0]) >= r_ex - 1e-12


def test_a_ghost_leaves_when_its_weight_falls_to_zero():
    """Noise-free two atoms entered through a ghost: the exit rule removes a unit of weight 0 with non-positive correlation, the blocks never contain a unit with w = 0 and a
    positive gradient."""
    c = np.array([20.0, 22.0])
    y = X.gauss(c) @ np.array([1.0, 0.6])
    R = X.path(y, 1e-3, 'dual', r_ex=0.0, max_entries=5)
    for b in R['blocks']:
        assert (b['w'] >= 0).all()
    last = R['blocks'][-1]
    assert np.isclose(last['rss'], 0.0, atol=1e-8) or len(last['w']) >= 2


def test_path_stops_at_max_entries_and_on_an_empty_signal():
    y = (X.gauss([20.0, 30.0, 40.0]) @ np.array([1.0, 1.0, 1.0]))
    R = X.path(y, 1e-3, 'dual', max_entries=2)
    assert len(R['taus']) == 2 and len(R['blocks']) == 3
    R0 = X.path(np.zeros(X.NS), 1.0, 'dual')
    assert len(R0['blocks']) == 1 and len(R0['taus']) == 0


# ---- the stopping rules -------------------------------------------------------------------------------------------------------------------------------------
def planted_blocks(rss):
    return [dict(c=np.arange(k, dtype=float), s=np.ones(k), w=np.ones(k), rss=float(r)) for k, r in enumerate(rss)]


def test_select_pen_first_minimum_and_penalty():
    B = planted_blocks([100.0, 40.0, 20.0, 10.0])
    noise = 1.0
    assert X.select_pen(B, noise, 5.0) == 3                     # 100, 45, 30, 25
    assert X.select_pen(B, noise, 10.0) == 2                    # 100, 50, 40, 40 -> the first minimum
    assert X.select_pen(B, noise, 1000.0) == 0
    assert X.select_pen(B, 2.0, 5.0) == 1 and X.select_pen(B, 20.0, 5.0) == 0 and X.select_pen(B, 0.5, 5.0) == 3     # chi^2 = rss / noise^2: 25, 10, 5, 2.5 plus 5 k, first minimum at 1


def test_select_sig_and_select_time():
    B = planted_blocks([9, 5, 3, 2, 1])
    assert X.select_sig(B, [6.0, 5.0, 3.0, 6.0], 4.5) == 2 and X.select_sig(B, [6.0, 5.0, 4.6, 4.5], 4.5) == 4 and X.select_sig(B, [2.0, 5.0], 4.5) == 0
    noise = 0.2
    t = 1.0 / (4.5 * noise * np.sqrt(X.NORM2))
    assert X.select_time(B, [0.5 * t, 0.9 * t, 1.1 * t], 4.5, noise) == 2 and X.select_time(B, [1.5 * t], 4.5, noise) == 0 and X.select_time(B, [0.5 * t, 0.6 * t, 0.7 * t, 0.8 * t], 4.5, noise) == 4


def test_select_disc_first_block_within_the_noise():
    noise = 0.5
    thr = noise ** 2 * (1 + 2 * np.sqrt(2 / X.NS)) * X.NS
    B = planted_blocks([10 * thr, 3 * thr, 0.99 * thr, 0.5 * thr])
    assert X.select_disc(B, noise) == 2
    assert X.select_disc(planted_blocks([10 * thr, 5 * thr]), noise) == 1       # none within the noise: the last block


def test_match_and_score():
    assert X.match(np.array([10.0, 12.0]), np.array([10.2, 11.5])) == (2, pytest.approx(np.sqrt((0.04 + 0.25) / 2)))
    n, rms = X.match(np.array([10.0, 12.0]), np.array([10.2, 13.0]))
    assert n == 1 and np.isclose(rms, 0.2, atol=1e-12)
    assert X.match(np.array([10.0]), np.array([10.7001]))[0] == 0 and X.match(np.array([10.0]), np.array([10.7]))[0] == 1
    n, rms = X.match(np.array([]), np.array([1.0]))
    assert n == 0 and np.isnan(rms)
    n, rms = X.match(np.array([10.0, 10.3]), np.array([10.1]))                  # one to one: the second estimate is not matched
    assert n == 1
    b = dict(c=np.array([10.0, 12.0]), s=np.ones(2), w=np.ones(2), rss=0.0)
    assert X.score(b, np.array([10.2, 11.5]))[0] == 1 and np.isclose(X.score(b, np.array([10.2, 11.5]))[1], np.sqrt(0.145), atol=1e-12)
    assert X.score(b, np.array([10.2, 13.0]))[0] == 0 and np.isnan(X.score(b, np.array([10.2, 13.0]))[1])
    assert X.score(dict(c=np.array([10.0]), s=np.ones(1), w=np.ones(1), rss=0.0), np.array([10.0, 12.0]))[0] == 0     # K-hat != K
    assert X.score(dict(c=np.array([10.0, 12.0, 14.0]), s=np.ones(3), w=np.ones(3), rss=0.0), np.array([10.0, 12.0]))[0] == 0


def test_same_block():
    a = dict(c=np.array([10.0, 12.0]), w=np.ones(2))
    assert X.same_block(a, dict(c=np.array([12.04, 10.0]), w=np.ones(2))) == 1
    assert X.same_block(a, dict(c=np.array([12.06, 10.0]), w=np.ones(2))) == 0
    assert X.same_block(a, dict(c=np.array([10.0]), w=np.ones(1))) == 0
    e = dict(c=np.zeros(0), w=np.zeros(0))
    assert X.same_block(e, e) == 1


# ---- the backward step --------------------------------------------------------------------------------------------------------------------------------------
def test_prune_removes_a_noise_unit_and_keeps_a_real_one():
    rng = np.random.default_rng(7)
    noise = 0.2
    y = (X.gauss([20.0, 30.0]) @ np.array([3.0, 3.0])) + noise * rng.standard_normal(X.NS)
    c = np.array([20.0, 30.0, 40.0])
    w = np.array([3.0, 3.0, 0.2])
    r = y - X.gauss(c) @ w
    b = dict(c=c, s=np.ones(3), w=w, rss=float(r @ r))
    p = X.prune_block(y, b, noise, 9.0)
    assert len(p['w']) == 2 and np.allclose(np.sort(p['c']), [20.0, 30.0], atol=0.2) and p['rss'] < 4.0
    p0 = X.prune_block(y, b, noise, -1e9)                                        # a penalty no removal can beat: nothing is removed
    assert len(p0['w']) == 3 and np.array_equal(p0['c'], c)
    e = dict(c=np.zeros(0), s=np.zeros(0), w=np.zeros(0), rss=float(y @ y))
    assert len(X.prune_block(y, e, noise, 9.0)['w']) == 0


def test_prune_keeps_strong_units():
    noise = 0.2
    y = X.gauss([20.0, 30.0]) @ np.array([3.0, 3.0])
    b = dict(c=np.array([20.0, 30.0]), s=np.ones(2), w=np.array([3.0, 3.0]), rss=0.0)
    assert len(X.prune_block(y, b, noise, 16.0)['w']) == 2


# ---- the scenes ---------------------------------------------------------------------------------------------------------------------------------------------
def test_pair_scene_geometry_and_significance():
    for seed in range(8):
        rng = np.random.default_rng(seed)
        c, a, s, noise = X.scene(rng, dict(kind='pair', d=1.6, ratio=3, snr=10))
        assert np.isclose(abs(c[0] - c[1]), 1.6, atol=1e-12) and np.array_equal(a, [3.0, 1.0]) and np.array_equal(s, [1.0, 1.0])
        assert 8.0 <= (c[0] + c[1]) / 2 <= X.LEN - 8.0
        assert np.isclose(np.min(a * np.sqrt(X.NORM2)) / noise, 10.0, atol=1e-12)
    c, a, s, noise = X.scene(np.random.default_rng(1), dict(kind='pair', d=2.0, ratio=1, snr=5, wr=2.0))
    assert np.array_equal(s, [1.0, 2.0]) and np.isclose(min(a[0] * X.feat_norm(1.0), a[1] * X.feat_norm(2.0)) / noise, 5.0, atol=1e-12)


def test_cluster_scene_min_separation_ratio_widths():
    for seed in range(10):
        c, a, s, noise = X.scene(np.random.default_rng(seed), dict(kind='cluster', K=4, dmin=2.2, ratio=5, snr=8, wr=1.5))
        assert np.diff(np.sort(c)).min() >= 2.2 and (c >= 6.0).all() and (c <= X.LEN - 6.0).all()
        assert np.isclose(a.min(), 1.0, atol=0) and a.max() <= 5.0 + 1e-12
        assert (s >= 1 / 1.5 - 1e-12).all() and (s <= 1.5 + 1e-12).all()
        assert np.isclose(np.min(a * np.array([X.feat_norm(x) for x in s])) / noise, 8.0, atol=1e-12)
    c, a, s, noise = X.scene(np.random.default_rng(0), dict(kind='cluster', K=3, dmin=2.0, ratio=1, snr=8))
    assert np.array_equal(a, np.ones(3)) and np.array_equal(s, np.ones(3))


def test_crowd_scene_window():
    for seed in range(10):
        c, a, s, noise = X.scene(np.random.default_rng(seed), dict(kind='crowd', K=5, gap=2.2, dmin=1.0, ratio=6, snr=8))
        assert c.max() - c.min() <= 5 * 2.2 + 1e-12 and np.diff(np.sort(c)).min() >= 1.0 and c.min() >= 6.0 and c.max() <= X.LEN - 6.0 + 1e-9
    with pytest.raises(RuntimeError):
        X.scene(np.random.default_rng(0), dict(kind='cluster', K=40, dmin=5.0, ratio=1, snr=8))


# ---- the comparator, the rules together, the trial -----------------------------------------------------------------------------------------------------------
def test_sfw_blocks_finds_a_strong_atom():
    y = (X.gauss([20.0]) * 3.0)[:, 0]
    B = X.sfw_blocks(y, 0.1)
    assert len(B) == 1 + len(X.LAM_REL) and len(B[0]['w']) == 0 and np.isclose(B[0]['rss'], y @ y, atol=0)
    assert any(len(b['w']) == 1 and abs(b['c'][0] - 20.0) < 0.1 and abs(b['w'][0] - 3.0) < 0.3 for b in B[1:])
    assert all(b['rss'] >= -1e-12 for b in B)
    B2 = X.sfw_blocks(0.01 * np.ones(X.NS) * np.cos(X.XS), 0.1)                    # weak, below min_sig: all spikes dropped
    assert all(len(b['w']) == 0 for b in B2)


def test_apply_rules_on_a_planted_path():
    c = np.array([20.0, 30.0])
    noise = 0.2
    rng = np.random.default_rng(3)
    y = X.gauss(c) @ np.array([3.0, 3.0]) + noise * rng.standard_normal(X.NS)
    R = X.path(y, noise, 'greedy', max_entries=5)
    out = X.apply_rules(y, R['blocks'], c, noise, None, nat=2)
    assert out.shape == (8, 2) and out[X.RULES.index('bic'), 0] == 1 and out[X.RULES.index('oracleK'), 0] == 1 and out[X.RULES.index('path'), 0] == 1
    assert out[X.RULES.index('nat'), 0] == 1 and out[X.RULES.index('bicbw'), 0] == 1 and out[X.RULES.index('penbw'), 0] == 1
    assert np.isfinite(out[:, 1]).all() and (out[:, 1] < 0.7).all()
    nobw = X.apply_rules(y, R['blocks'], c, noise, None, bw=False)
    assert np.isnan(nobw[X.RULES.index('bicbw')]).all() and np.isnan(nobw[X.RULES.index('nat')]).all() and nobw[X.RULES.index('bic'), 0] == 1


def test_trial_is_deterministic_and_shaped():
    fam = dict(kind='pair', d=1.6, ratio=1, snr=10)
    o1, a1 = X.trial_full(fam, 5)
    o2, a2 = X.trial_full(fam, 5)
    assert o1.shape == (3, 8, 2) and np.array_equal(o1, o2, equal_nan=True) and np.array_equal(a1, a2) and a1.shape == (2,)
    assert np.array_equal(X.trial(fam, 5), o1, equal_nan=True)
    assert np.isnan(o1[2, X.RULES.index('nat')]).all() and np.isnan(o1[2, X.RULES.index('bicbw')]).all()          # sFW has no natural rule and no backward step
    assert np.array_equal(X.trial_full(fam, 5, methods=('dual', 'greedy'))[0], o1[:2], equal_nan=True)
    assert np.array_equal(X.trial_full(fam, 6)[0][:, :, 0], X.trial_full(fam, 6)[0][:, :, 0], equal_nan=True)
    assert np.isclose(o1[0, X.RULES.index('bic'), 0], 1.0) or o1[0, X.RULES.index('bic'), 0] == 0.0


def test_trial_planted_values_at_seed_100002():
    fam = dict(kind='pair', d=1.6, ratio=1, snr=10)
    o, ag = X.trial_full(fam, 100002)
    b = X.RULES.index('bic')
    assert o[0, b, 0] == 1 and o[1, b, 0] == 0                                  # the dual recovers the pair under BIC, the bare greedy order does not
    assert o[1, X.RULES.index('bicbw'), 0] == 1 and o[0, X.RULES.index('oracleK'), 0] == 1 and o[1, X.RULES.index('oracleK'), 0] == 0
    assert ag[0] == 0


def test_trial_with_widths_and_crowd():
    o, ag = X.trial_full(dict(kind='crowd', K=3, gap=3.0, dmin=1.0, ratio=1, snr=14), 11)
    assert o.shape == (3, 8, 2)
    o2, _ = X.trial_full(dict(kind='cluster', K=2, dmin=2.5, ratio=1, snr=12, wr=1.5), 12, widths=X.WIDTHS, methods=('dual', 'greedy'))
    assert o2.shape == (2, 8, 2)


# ---- statistics ---------------------------------------------------------------------------------------------------------------------------------------------
def test_paired_boot_and_verdict():
    d = np.r_[np.ones(60), np.zeros(40)] * 100.0
    m, lo, hi = X.paired_boot(d)
    assert np.isclose(m, 60.0, atol=0) and lo < 60.0 < hi and 49.0 < lo < 52.0 and 68.0 < hi < 71.0
    assert X.paired_boot(d) == X.paired_boot(d)
    assert X.paired_boot(np.zeros(10)) == (0.0, 0.0, 0.0)
    assert all(np.isnan(v) for v in X.paired_boot(np.array([np.nan])))
    m2, _, _ = X.paired_boot(np.array([1.0, np.nan, 3.0]))
    assert np.isclose(m2, 2.0, atol=0)
    assert X.verdict(0.1, 2.0) == 'win' and X.verdict(-2.0, -0.1) == 'loss' and X.verdict(-1.0, 1.0) == 'tie' and X.verdict(0.0, 1.0) == 'tie' and X.verdict(-1.0, 0.0) == 'tie'


def test_clopper_pearson():
    lo, hi = X.clopper_pearson(5, 10)
    assert np.isclose(lo, 0.18709, atol=1e-4) and np.isclose(hi, 0.81291, atol=1e-4)
    lo, hi = X.clopper_pearson(0, 10)
    assert lo == 0.0 and np.isclose(hi, 0.30850, atol=1e-4)
    lo, hi = X.clopper_pearson(10, 10)
    assert hi == 1.0 and np.isclose(lo, 0.69150, atol=1e-4)
    lo, hi = X.clopper_pearson(800, 800)
    assert np.isclose(lo, 0.99540, atol=2e-4)


def test_diff_stats_and_pooling():
    E = np.zeros((10, 3, 8, 2))
    E[:5, 0, 0, 0] = 1
    E[:2, 1, 0, 0] = 1
    m, se = X.diff_stats(E, 'dual', 'greedy', 'bic')
    assert np.isclose(m, 30.0, atol=1e-12) and np.isclose(se, 100 * np.std([1, 1, 0, 0, 0, 0, 0, 0, 0, 0] if False else [0, 0, 1, 1, 1, 0, 0, 0, 0, 0], ddof=1) / np.sqrt(10), atol=1e-12)
    raw = {'G|0|E': E, 'G|1|E': E.copy()}
    assert X.cell_keys(raw, 'G') == [0, 1] and X.cell_keys(raw, 'H') == []
    assert X.pooled(raw, 'G', 'dual', 'bic').shape == (20,) and X.pooled(raw, 'G', 'dual', 'bic', [1]).shape == (10,)
    assert np.array_equal(X.exact_of(E, 'dual', 'bic'), E[:, 0, 0, 0])
    rms = np.full((4, 3, 8, 2), np.nan)
    rms[:, 1, 0, 1] = [1, 2, 3, 4]
    assert np.array_equal(X.pooled_rms({'G|0|E': rms}, 'G', 'greedy', 'bic'), [1, 2, 3, 4])


def test_fam_label_and_fam_of():
    assert X.fam_label(dict(kind='pair', d=1.6, ratio=3, snr=10)) == 'd=1.6 ratio=3 snr=10'
    assert X.fam_label(dict(kind='crowd', K=5, gap=2.2, dmin=1.0, ratio=6, snr=8)) == 'K=5 gap=2.2 dmin=1 ratio=6 snr=8'
    assert X.fam_of({'G|3|fam': np.array(json.dumps(dict(kind='pair', d=1.0)))}, 'G', 3) == dict(kind='pair', d=1.0)
    assert X.pct(0.456) == '46'


# ---- the declared grids and seeds ---------------------------------------------------------------------------------------------------------------------------
def test_declared_grids():
    assert {g: (n, len(c)) for g, (n, c) in X.DEV_GRIDS.items()} == {'A': (60, 96), 'B': (40, 54), 'C': (40, 32), 'D': (50, 24)}
    assert {g: (n, len(c)) for g, (n, c) in X.TEST_GRIDS.items()} == {'P': (250, 8), 'E': (100, 8), 'Q': (100, 24), 'W': (120, 16)}
    P = [f for f, w in X.TEST_GRIDS['P'][1]]
    assert sorted(set(f['d'] for f in P)) == [1.2, 1.6] and sorted(set(f['ratio'] for f in P)) == [1, 3] and sorted(set(f['snr'] for f in P)) == [8, 12]
    E = [f for f, w in X.TEST_GRIDS['E'][1]]
    assert sorted(set(f['d'] for f in E)) == [3.5, 5.0] and sorted(set(f['ratio'] for f in E)) == [1, 10] and sorted(set(f['snr'] for f in E)) == [8, 12]
    assert X.TEST_GRIDS['Q'][1] == X.DEV_GRIDS['D'][1]
    W = X.TEST_GRIDS['W'][1]
    assert all(f['K'] == 4 and w == X.WIDTHS for f, w in W) and W == [c for c in X.DEV_GRIDS['C'][1] if c[0]['K'] == 4] and sorted(set(f['wr'] for f, w in W)) == [1.5, 2.0]
    assert all(w == X.WIDTHS for f, w in X.DEV_GRIDS['C'][1]) and all(w is None for g in 'ABD' for f, w in X.DEV_GRIDS[g][1]) and all(f['wr'] in (1.5, 2.0) for f, w in X.DEV_GRIDS['C'][1])
    assert sorted(set(f['d'] for f, w in X.DEV_GRIDS['A'][1])) == [1.2, 1.6, 2.0, 2.5, 3.5, 5.0]
    assert sorted(set(f['K'] for f, w in X.DEV_GRIDS['B'][1])) == [4, 7] and sorted(set(f['K'] for f, w in X.DEV_GRIDS['D'][1])) == [3, 5]
    assert all(f['kind'] == 'crowd' and f['dmin'] == 1.0 for f, w in X.DEV_GRIDS['D'][1])


def test_seed_blocks_do_not_overlap_and_development_is_below_test():
    dev = X.grid_seeds(X.DEV_GRIDS, X.DEV_SEED0, 1000)
    test = X.grid_seeds(X.TEST_GRIDS, X.TEST_SEED0, 10000)
    assert dev['A'][0] == 1000000 and dev['A'][1] == 1001000 and dev['B'][0] == 1096000 and dev['D'][-1] == 1000000 + 1000 * (96 + 54 + 32 + 23)
    assert test['P'][0] == 7000000 and test['P'][1] == 7010000 and test['E'][0] == 7080000 and test['Q'][0] == 7160000 and test['W'][0] == 7400000 and test['W'][-1] == 7550000
    top_dev = max(s + X.DEV_GRIDS[g][0] for g in dev for s in dev[g])
    assert top_dev < X.TEST_SEED0 and X.TEST_SEED0 + 10000 * 56 < X.TWOD_SEED0
    all_test = sorted(s for g in test for s in test[g])
    assert all(b - a >= 10000 for a, b in zip(all_test, all_test[1:])) and max(X.TEST_GRIDS[g][0] for g in X.TEST_GRIDS) <= 10000
    all_dev = sorted(s for g in dev for s in dev[g])
    assert all(b - a >= 1000 for a, b in zip(all_dev, all_dev[1:])) and max(X.DEV_GRIDS[g][0] for g in X.DEV_GRIDS) <= 1000


def test_comparisons_are_the_declared_ones():
    ids = [c[0] for c in X.COMPARISONS]
    assert ids == [str(k) for k in range(1, 11)]
    d = {c[0]: (c[2], c[3]) for c in X.COMPARISONS}
    assert d['1'] == (('dual', 'bicbw'), ('greedy', 'bicbw')) and d['2'] == (('dual', 'bic'), ('greedy', 'bic')) and d['3'] == (('dual', 'bic'), ('greedy', 'bicbw'))
    assert d['4'] == (('dual', 'bic'), ('sfw', 'bic')) and d['9'] == (('dual', 'oracleK'), ('greedy', 'oracleK')) and d['10'] == (('dual', 'nat'), ('greedy', 'nat'))
    assert d['6'] == (('dual', 'pen'), ('greedy', 'pen')) and d['7'] == (('dual', 'penbw'), ('greedy', 'penbw')) and d['8'] == (('dual', 'pen'), ('sfw', 'pen')) and d['5'] == (('dual', 'bicbw'), ('sfw', 'bic'))


# ---- part T and the two-dimensional arm's bookkeeping ---------------------------------------------------------------------------------------------------------
def test_part_t_declared_and_shaped(monkeypatch):
    assert X.T_REX == (0.0, 1.5) and X.T_DS == (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0) and X.T_RATIOS == (1.0, 2.0, 5.0, 20.0) and (X.T_N, X.T_SEED) == (12, 31)
    monkeypatch.setattr(X, 'T_REX', (1.5,))
    monkeypatch.setattr(X, 'T_DS', (1.0, 6.0))
    monkeypatch.setattr(X, 'T_RATIOS', (2.0,))
    monkeypatch.setattr(X, 'T_N', 3)
    T = X.part_t()
    assert T.shape == (1, 2, 1, 3, 4)
    assert (T[0, 0, 0, :, 0] == 0).all() and (T[0, 1, 0, :, 0] == 1).all()        # separation 1: the entry sequences differ; separation 6: the same
    assert (T[0, :, 0, :, 2:] == 1).all()                                        # both blocks exact


def test_twod_constants():
    assert X.TWOD_STOP == ('bic', 4.0) and X.TWOD_SEPS == (1.2, 1.6, 2.5, 3.0) and X.TWOD_PHOTONS == (0, 1) and X.TWOD_SCENES == 100 and X.TWOD_SEED0 == 9000000
    assert X.TWOD_CLOSE == (1.2, 1.6) and X.TWOD_MID == (2.5, 3.0)


def test_twod_lines_on_planted_classes():
    tw = np.ones((4, 2, 100, 2), int) * 2
    tw[0:2, :, :, 0] = 1                                                          # close separations: the dual resolves, the greedy order merges
    tw[2:4, :, :, 1] = 1                                                          # mid separations: the greedy order resolves
    L, out = X.twod_lines(tw)
    assert out == {'close': 'win', 'mid': 'loss'}
    assert sum(1 for l in L if l.startswith('| 1.2')) == 2 and any('+100.0 [+100.0, +100.0]' in l for l in L) and any('-100.0 [-100.0, -100.0]' in l for l in L)


# ---- the SMLM post hoc table ---------------------------------------------------------------------------------------------------------------------------------
def planted_smlm():
    rng = np.random.default_rng(0)
    n = 400
    nn = np.tile([100.0, 250.0, 350.0, 600.0], 100)                             # 100 frames, one emitter of each band per frame
    fr = np.repeat(np.arange(100), 4)
    hi = rng.random(n) < 0.5
    hg = hi.copy()
    hg[1::4] = rng.random(100) < 0.95                                            # greedy better in the 200-300 band
    return {'S|nn': nn, 'S|frame': fr, 'S|hit_iss': hi, 'S|hit_gr': hg}


def test_smlm_lines_planted():
    L = X.smlm_lines({k.replace('S|', 'MT0.N1.HD|'): v for k, v in planted_smlm().items()})
    rows = [l for l in L if l.startswith('| MT0.N1.HD')]
    assert len(rows) == 4 and 'tie' in rows[0] and 'loss' in rows[1] and rows[1].count('[200, 300)') == 1
    assert any('Tally' in l and '[200, 300): 0 win, 0 tie, 1 loss' in l and '[0, 200): 0 win, 1 tie, 0 loss' in l for l in L)


# ---- raw files and the report -------------------------------------------------------------------------------------------------------------------------------
def test_save_and_load_raw_merge_atomically(tmp_path):
    p = str(tmp_path / 'sub' / 'r.npz')
    assert X.load_raw(p) == {}
    X.save_raw(p, {'a': np.arange(3)})
    X.save_raw(p, {'b': np.array('x')})
    d = X.load_raw(p)
    assert set(d) == {'a', 'b'} and np.array_equal(d['a'], np.arange(3)) and str(d['b']) == 'x'
    X.save_raw(p, {'a': np.arange(2)})
    assert np.array_equal(X.load_raw(p)['a'], np.arange(2)) and not os.path.exists(p + '.tmp.npz')


def tiny_grids(monkeypatch):
    fam1 = dict(kind='pair', d=1.6, ratio=1, snr=12)
    fam2 = dict(kind='pair', d=5.0, ratio=1, snr=12)
    monkeypatch.setattr(X, 'DEV_GRIDS', {'A': (3, [(fam1, None), (fam2, None)])})
    monkeypatch.setattr(X, 'TEST_GRIDS', {'P': (4, [(fam1, None), (fam2, None)]), 'E': (3, [(fam2, None)]), 'Q': (3, [(dict(kind='crowd', K=3, gap=3.0, dmin=1.0, ratio=1, snr=14), None)])})


def test_run_grid_layout_and_reproducibility(monkeypatch):
    tiny_grids(monkeypatch)
    r1 = X.run_grid(X.TEST_GRIDS, X.TEST_SEED0, 10000, ['P'])
    r2 = X.run_grid(X.TEST_GRIDS, X.TEST_SEED0, 10000, ['P'])
    assert set(r1) == {'P|0|fam', 'P|0|widths', 'P|0|E', 'P|0|A', 'P|1|fam', 'P|1|widths', 'P|1|E', 'P|1|A'}
    assert r1['P|0|E'].shape == (4, 3, 8, 2) and r1['P|0|A'].shape == (4, 2) and json.loads(str(r1['P|0|widths'])) is None
    assert all((np.array_equal(r1[k], r2[k], equal_nan=True) if r1[k].dtype.kind == 'f' else str(r1[k]) == str(r2[k])) for k in r1)
    o, a = X.trial_full(dict(kind='pair', d=1.6, ratio=1, snr=12), X.TEST_SEED0 + 2)
    assert np.array_equal(r1['P|0|E'][2], o, equal_nan=True) and np.array_equal(r1['P|0|A'][2], a)
    o, a = X.trial_full(dict(kind='pair', d=5.0, ratio=1, snr=12), X.TEST_SEED0 + 10000 + 3)
    assert np.array_equal(r1['P|1|E'][3], o, equal_nan=True)
    r3 = X.run_grid(X.TEST_GRIDS, X.TEST_SEED0, 10000, ['E'])
    o, a = X.trial_full(dict(kind='pair', d=5.0, ratio=1, snr=12), X.TEST_SEED0 + 20000)       # the E cell follows the two P cells in the seed order
    assert np.array_equal(r3['E|0|E'][0], o, equal_nan=True)


def test_report_on_a_tiny_run(monkeypatch, tmp_path):
    tiny_grids(monkeypatch)
    raw_p, dev_p = str(tmp_path / 'raw.npz'), str(tmp_path / 'dev.npz')
    X.save_raw(dev_p, X.run_grid(X.DEV_GRIDS, X.DEV_SEED0, 1000, ['A']))
    X.save_raw(raw_p, X.run_grid(X.TEST_GRIDS, X.TEST_SEED0, 10000, ['P', 'E', 'Q']))
    monkeypatch.setattr(X, 'T_REX', (1.5,))
    monkeypatch.setattr(X, 'T_DS', (1.0, 6.0))
    monkeypatch.setattr(X, 'T_RATIOS', (2.0,))
    monkeypatch.setattr(X, 'T_N', 2)
    X.save_raw(raw_p, {'T': X.part_t()})
    X.save_raw(raw_p, {'smlm|' + k: v for k, v in planted_smlm().items()} | {'smlm|S|x': np.zeros(1)} if False else {'smlm|' + k.replace('S|', 'MT0.N1.HD|'): v for k, v in planted_smlm().items()})
    t1 = X.report(raw_p, dev_p)
    t2 = X.report(raw_p, dev_p)
    assert t1 == t2 and t1.startswith('# Is there a regime') and t1.endswith('\n')
    for h in ('## Part T', '## Held-out run', '### Regime P (close pairs)', '### Regime Q (crowded windows)', '### Regime E', '## POST HOC: the SMLM', '## Development scan', '### Grid A (3 scenes per cell'):
        assert h in t1
    assert 'Plan rule: the dual order matters in P iff P1, P3, P4 are all wins' in t1 and 'Plan rule: equivalent iff the lower end is at least 98 %' in t1
    assert t1.count('| P1 |') == 1 and t1.count('| Q10 |') == 1


def test_test_lines_verdict_rule_on_planted_cells():
    """A regime P in which the dual beats everything is 'yes'; if the comparator sFW is better it is 'no'."""
    n = 200
    fam = np.array(json.dumps(dict(kind='pair', d=1.6, ratio=1, snr=10)))
    E = np.zeros((n, 3, 8, 2))
    ib = X.RULES.index('bic')
    for r in ('bic', 'bicbw', 'pen', 'penbw', 'oracleK', 'nat'):
        E[:, 0, X.RULES.index(r), 0] = 1
    E[:100, 1, ib, 0] = 1
    E[:100, 1, X.RULES.index('bicbw'), 0] = 1
    E[:100, 2, ib, 0] = 1
    E[:, :, :, 1] = 0.1
    raw = {'P|0|E': E, 'P|0|fam': fam, 'P|0|widths': np.array('null')}
    L, res = X.test_lines(raw)
    assert res['P']['1'] == res['P']['3'] == res['P']['4'] == 'win' and res['P']['regime'] == 'yes'
    E2 = E.copy()
    E2[:, 2, ib, 0] = 1
    L, res = X.test_lines({'P|0|E': E2, 'P|0|fam': fam, 'P|0|widths': np.array('null')})
    assert res['P']['4'] == 'tie' and res['P']['regime'] == 'no'
    E3 = E.copy()
    E3[:, 1, X.RULES.index('bicbw'), 0] = 1
    L, res = X.test_lines({'P|0|E': E3, 'P|0|fam': fam, 'P|0|widths': np.array('null')})
    assert res['P']['3'] == 'tie' and res['P']['regime'] == 'no'


def test_equivalence_rule_on_planted_cells():
    n = 400
    fam = np.array(json.dumps(dict(kind='pair', d=5.0, ratio=1, snr=10)))
    E = np.zeros((n, 3, 8, 2))
    A = np.ones((n, 2), int)
    L, res = X.test_lines({'E|0|E': E, 'E|0|A': A, 'E|0|fam': fam, 'E|0|widths': np.array('null'), 'E|1|E': E, 'E|1|A': A, 'E|1|fam': fam, 'E|1|widths': np.array('null')})
    assert res['equiv'] == 'yes' and any('800 of 800' in l for l in L)
    A2 = A.copy()
    A2[:30, 0] = 0
    L, res = X.test_lines({'E|0|E': E, 'E|0|A': A2, 'E|0|fam': fam, 'E|0|widths': np.array('null'), 'E|1|E': E, 'E|1|A': A2, 'E|1|fam': fam, 'E|1|widths': np.array('null')})
    assert res['equiv'] == 'no'


def test_stored_raw_regenerates_the_stored_md(tmp_path):
    raw = os.path.join(HERE, '..', 'results', 'synthetic', 'dualregime_raw.npz')
    dev = os.path.join(HERE, '..', 'results', 'synthetic', 'dualregime_dev_raw.npz')
    md = os.path.join(HERE, '..', 'results', 'synthetic', 'dualregime.md')
    if not (os.path.exists(raw) and os.path.exists(md)):
        pytest.skip('no stored report yet')
    assert X.report(raw, dev if os.path.exists(dev) else str(tmp_path / 'none.npz')) == open(md).read()


def test_stored_raw_has_the_declared_shapes():
    raw = os.path.join(HERE, '..', 'results', 'synthetic', 'dualregime_raw.npz')
    if not os.path.exists(raw):
        pytest.skip('no stored raw yet')
    d = X.load_raw(raw)
    for g in ('P', 'E', 'Q', 'W'):
        if '%s|0|E' % g in d:
            n, cells = X.TEST_GRIDS[g]
            assert X.cell_keys(d, g) == list(range(len(cells)))
            for j, (fam, wd) in enumerate(cells):
                assert d['%s|%d|E' % (g, j)].shape == (n, 3, 8, 2) and d['%s|%d|A' % (g, j)].shape == (n, 2)
                assert X.fam_of(d, g, j) == fam
    if 'twod' in d:
        assert d['twod'].shape == (4, 2, X.TWOD_SCENES, 2)
    if 'T' in d:
        assert d['T'].shape == (2, 9, 4, 12, 4)


# ---- the generated summary and the development tallies -------------------------------------------------------------------------------------------------------
def planted_dev(monkeypatch):
    monkeypatch.setattr(X, 'DEV_GRIDS', {'A': (100, [(dict(kind='pair', d=1.6, ratio=1, snr=10), None)] * 4)})
    fam = np.array(json.dumps(dict(kind='pair', d=1.6, ratio=1, snr=10)))
    dev = {}
    ib, iw = X.RULES.index('bic'), X.RULES.index('bicbw')
    for j in range(4):
        E = np.zeros((100, 3, 8, 2))
        if j == 0:                                   # dual ahead by 30 under BIC
            E[:60, 0, ib, 0] = 1
            E[:30, 1, ib, 0] = 1
        elif j == 1:                                 # greedy ahead
            E[:30, 0, ib, 0] = 1
            E[:60, 1, ib, 0] = 1
        elif j == 2:                                 # identical
            E[:50, 0, ib, 0] = 1
            E[:50, 1, ib, 0] = 1
        else:                                        # a difference below two standard errors
            E[:50, 0, ib, 0] = 1
            E[:46, 1, ib, 0] = 1
            E[50:54, 1, ib, 0] = 1
        E[:50, 2, ib, 0] = 1                         # sFW at 50 % everywhere
        dev['A|%d|E' % j], dev['A|%d|fam' % j], dev['A|%d|widths' % j] = E, fam, np.array('null')
    return dev


def test_dev_tally_counts_cells(monkeypatch):
    dev = planted_dev(monkeypatch)
    t = X.dev_tally(dev)
    up, dn, nl, tot = t['dual - greedy, BIC']
    assert (up, dn, nl, tot) == (1, 1, 1, 4)                        # cell 3: 4 points apart with a standard error of about 6.6: neither
    assert t['dual + backward step - greedy + backward step'] == (0, 0, 4, 4)
    assert t['dual - sFW, BIC'][3] == 4 and t['dual - sFW, BIC'][0] == 1


def test_dev_lines_rows(monkeypatch):
    dev = planted_dev(monkeypatch)
    L, flags = X.dev_lines(dev)
    rows = [l for l in L if l.startswith('| d=1.6')]
    assert len(rows) == 4 and '+30 (' in rows[0] and '-30 (' in rows[1] and '+0 (0)' in rows[2] and len(flags) == 4
    assert '| 60 | 30 | 50 |' in rows[0]


def test_summary_lines_from_planted_results(monkeypatch):
    dev = planted_dev(monkeypatch)
    T = np.zeros((2, 9, 4, 12, 4))
    T[:, 7:, :, :, 0] = 1                                            # the sequences agree from 5 widths (indices 7, 8)
    raw = {'T': T}
    res = {'P': dict(regime='no', **{str(k): 'tie' for k in range(1, 11)}), 'equiv': 'yes', 'twod': dict(close='win', mid='loss')}
    res['P']['2'] = 'win'
    L = X.summary_lines(raw, res, dev)
    txt = '\n'.join(L)
    assert 'Regime P (close pairs): (i) the order effect with everything else equal' in txt and 'as the pipelines stand (P2) a win' in txt
    assert 'the same backward step on both, P1): a tie; as the pipelines stand (P2) a win; against greedy + backward step (P3) a tie; first block with K units (P9, oracle K) a tie;' in txt
    assert '(ii) the comparison with sliding Frank-Wolfe, a different algorithm and not an order comparison (P4): a tie.' in txt and 'Registered plan rule (P1, P3, P4 all wins): **no**, neither (i) nor (ii) is met.' in txt
    assert 'P10 (each order\'s own stopping rule) is not an order comparison and is not to be quoted.' in txt
    assert 'equivalent (plan rule) **yes**.  Scope: one-dimensional pairs with two features, the final BIC blocks (sorted centres within 0.05 widths), not the entry sequences and not images.' in txt
    assert 'equivalent (plan rule) **yes**' in txt and 'close separations (1.2, 1.6 widths) dual vs greedy **win**, mid separations (2.5, 3.0 widths) **loss**' in txt
    assert 'r_ex 0): the entry sequences of the two orders are the same in 96 of 96 phases at separations of 5 widths and more, and in 0 of 288 phases at separations of 3.5 widths and less' in txt
    assert 'Development scan, dual - greedy, BIC: of 4 cells, 1 with the difference at least two standard errors above 0, 1 at least two below, 1 never differ' in txt
    assert L[-1] == ''


def test_conclusion_lines_follow_the_verdicts():
    base = {str(k): 'tie' for k in range(1, 11)}
    res = {'P': dict(base, regime='no', **{'2': 'win', '4': 'loss'}), 'Q': dict(base, regime='yes'), 'W': dict(base, regime='no'), 'equiv': 'yes', 'twod': dict(close='win', mid='loss')}
    L = X.conclusion_lines(res)
    txt = '\n'.join(L)
    assert L[-1] == '' and L[0].startswith('1. **Equivalence.**') and 'plan rule met' in txt
    assert 'is a win (P2)' in txt and 'comparison with a different algorithm, not an order comparison, and the dual is a loss there (P4)' in txt
    assert 'matters in P is NOT met' in txt and 'matters in Q is  met' in txt
    assert 'Close separations: dual against greedy a win; mid separations: a loss' in txt and '3. **Regime W' in txt and 'The plan rule (W1 and W3) is NOT met' in txt
    assert '5. **Mechanism (a hypothesis, not a finding).**' in txt and '6. **What the paper may claim.**' in txt
    assert X.conclusion_lines({'equiv': 'no'})[0].endswith('(see its table).')


# ---- tests added after the first mutation run ---------------------------------------------------------------------------------------------------------------
def test_an_exit_is_recorded_and_the_blocks_follow_it():
    """Seed 14 of a crowd (3 features, gap 1.6, ratio 6, significance 14), dual order, no resolution limit: the first entrant (lattice point 34.75) leaves when its neighbours have entered
    (its weight falls to 0 with a non-positive correlation); the sizes of the blocks, the active ids and the chi^2 after the exit are pinned."""
    rng = np.random.default_rng(14)
    c, a, s, noise = X.scene(rng, dict(kind='crowd', K=3, gap=1.6, dmin=1.0, ratio=6, snr=14))
    y = X.gauss(c, s) @ a + noise * rng.standard_normal(X.NS)
    R = X.path(y, noise, 'dual', r_ex=0.0, max_entries=6)
    assert [float(X.LAT[j]) for j in R['entered']] == [34.75, 35.0, 34.5, 32.5, 40.5, 12.25]
    assert [float(X.LAT[j]) for j in R['ids']] == [35.0, 34.5, 32.5, 40.5, 12.25]
    assert [len(b['w']) for b in R['blocks']] == [0, 1, 2, 3, 3, 4, 5]
    assert np.allclose([b['rss'] for b in R['blocks']], [39.3687, 6.3203, 1.4138, 1.2448, 1.1506, 1.1167, 1.0851], atol=2e-4)
    assert np.allclose(R['blocks'][4]['c'], [36.5, 35.01, 32.951], atol=2e-3) and np.allclose(R['blocks'][4]['w'], [0.7433, 2.1678, 1.289], atol=2e-3)
    for b in R['blocks']:
        assert len(b['c']) == len(b['s']) == len(b['w'])
    # the stored chi^2 of every block is the residual of that block
    for b in R['blocks']:
        r = y - (X.gauss(b['c'], b['s']) @ b['w'] if len(b['w']) else 0.0)
        assert np.isclose(float(r @ r), b['rss'], rtol=1e-9)


def test_block_fit_matches_a_numerical_jacobian_fit_on_noisy_data():
    """The analytic Jacobian of the centre and width columns: with noise the fit must reach the stationary point of the true cost (a wrong column stops elsewhere)."""
    from scipy.optimize import least_squares
    rng = np.random.default_rng(5)
    y = (X.gauss([20.0, 22.5], [1.2, 0.8]) @ np.array([2.0, 1.5])) + 0.3 * rng.standard_normal(X.NS)
    c0 = np.array([20.0, 22.5])
    c, s, w = X.block_fit(y, c0, np.ones(2), np.zeros(2), c0, (0.4, 3.0), tol=1e-12, max_nfev=400)
    sol = least_squares(lambda z: y - X.gauss(z[2:4], z[4:6]) @ z[:2], np.r_[w, c, s], xtol=1e-14, ftol=1e-14, gtol=1e-14)
    assert np.allclose(np.r_[w, c, s], sol.x, atol=1e-4)
    z = np.r_[w, c, s]
    cost = lambda q: 0.5 * np.sum((y - X.gauss(q[2:4], q[4:6]) @ q[:2]) ** 2)      # noqa: E731
    grad = np.array([(cost(z + 1e-6 * e) - cost(z - 1e-6 * e)) / 2e-6 for e in np.eye(6)])
    assert np.abs(grad).max() < 3e-6                                                  # a stationary point of the true cost (a wrong Jacobian column stops at about 7e-6)
    c1, s1, w1 = X.block_fit(y, c0, np.ones(2), np.zeros(2), c0, None, tol=1e-12, max_nfev=400)
    sol1 = least_squares(lambda z: y - X.gauss(z[2:4], 1.0) @ z[:2], np.r_[w1, c1], xtol=1e-14, ftol=1e-14, gtol=1e-14)
    assert np.allclose(np.r_[w1, c1], sol1.x, atol=1e-4) and np.array_equal(s1, np.ones(2))
    z1 = np.r_[w1, c1]
    cost1 = lambda q: 0.5 * np.sum((y - X.gauss(q[2:4], 1.0) @ q[:2]) ** 2)         # noqa: E731
    grad1 = np.array([(cost1(z1 + 1e-6 * e) - cost1(z1 - 1e-6 * e)) / 2e-6 for e in np.eye(4)])
    assert np.abs(grad1).max() < 3e-6


def test_block_fit_bounds_on_both_sides_and_infeasible_start():
    y = (X.gauss([20.0], [0.3]) * 2.0)[:, 0]
    c, s, w = X.block_fit(y, np.array([20.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), (0.4, 3.0))
    assert np.isclose(s[0], 0.4, atol=1e-9)                                        # the lower width bound is active
    c, s, w = X.block_fit(y, np.array([20.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), (1.2, 3.0))
    assert np.isclose(s[0], 1.2, atol=1e-9)                                        # a start outside the bounds is clipped into them
    far = (X.gauss([16.0]) * 2.0)[:, 0]
    c, s, w = X.block_fit(far, np.array([20.0]), np.array([1.0]), np.array([0.0]), np.array([20.0]), None, box=1.5)
    assert c[0] >= 18.5 - 1e-9 and np.isclose(c[0], 18.5, atol=1e-6)               # the lower side of the box


def test_rules_at_equality_and_with_several_bad_entrants():
    B = planted_blocks([9, 5, 3, 2, 1])
    assert X.select_sig(B, [6.0, 3.0, 6.0, 3.0], 4.5) == 1                          # the FIRST entrant below kappa
    assert X.select_sig(B, [6.0, 4.5, 4.6], 4.5) == 3                               # an entrant exactly at kappa counts
    noise = 0.2
    t = 1.0 / (4.5 * noise * np.sqrt(X.NORM2))
    assert X.select_time(B, [t], 4.5, noise) == 1                                   # an entry exactly at the stopping time counts
    assert X.match(np.array([10.0]), np.array([10.5]), tol=0.5)[0] == 1 and X.match(np.array([10.0]), np.array([10.5]), tol=0.4999)[0] == 0
    a = dict(c=np.array([10.0]), w=np.ones(1))
    assert X.same_block(a, dict(c=np.array([10.5]), w=np.ones(1)), tol=0.5) == 1 and X.same_block(a, dict(c=np.array([10.5]), w=np.ones(1)), tol=0.4999) == 0


def test_prune_can_empty_a_block():
    noise = 0.2
    y = 0.05 * np.cos(X.XS)
    r = y - X.gauss([20.0]) @ np.array([0.05])
    b = dict(c=np.array([20.0]), s=np.ones(1), w=np.array([0.05]), rss=float(r @ r))
    assert len(X.prune_block(y, b, noise, 9.0)['w']) == 0


def planted_path_blocks(gain):
    """Blocks with 0, 1, 2 units; the second unit lowers chi^2 (noise 1) by `gain`; the first by 100."""
    return [dict(c=np.arange(k, dtype=float) + 10.0, s=np.ones(k), w=np.ones(k), rss=r) for k, r in enumerate([200.0, 100.0, 100.0 - gain])]


def test_apply_rules_penalties_and_missing_blocks():
    c_true = np.array([10.0, 11.0])
    B = planted_path_blocks(11.0)                                                   # a gain of 11: kept by the BIC penalty 9.1, dropped by 13.7 (three parameters) and by 16
    out = X.apply_rules(np.zeros(X.NS), B, c_true, 1.0, None, bw=False)
    assert out[X.RULES.index('bic'), 0] == 1 and out[X.RULES.index('pen'), 0] == 0
    out3 = X.apply_rules(np.zeros(X.NS), B, c_true, 1.0, X.WIDTHS, bw=False)
    assert out3[X.RULES.index('bic'), 0] == 0                                        # three parameters per unit: penalty 13.7 > 11
    B2 = planted_path_blocks(15.0)                                                  # a gain of 15 between 13.7 and 16
    assert X.apply_rules(np.zeros(X.NS), B2, c_true, 1.0, X.WIDTHS, bw=False)[X.RULES.index('bic'), 0] == 1
    assert X.apply_rules(np.zeros(X.NS), B2, c_true, 1.0, None, bw=False)[X.RULES.index('pen'), 0] == 0
    only1 = B[:2]                                                                    # no block with K = 2 units: both oracle rules fail
    o = X.apply_rules(np.zeros(X.NS), only1, c_true, 1.0, None, bw=False)
    assert o[X.RULES.index('oracleK'), 0] == 0 and o[X.RULES.index('path'), 0] == 0 and np.isnan(o[X.RULES.index('oracleK'), 1])


def test_trial_runs_paths_to_k_plus_three_entries(monkeypatch):
    seen = []
    real = X.path

    def spy(*a, **k):
        seen.append(k.get('max_entries'))
        return real(*a, **k)
    monkeypatch.setattr(X, 'path', spy)
    X.trial_full(dict(kind='pair', d=3.0, ratio=1, snr=12), 3, methods=('dual', 'greedy'))
    assert seen == [5, 5]
    seen.clear()
    X.trial_full(dict(kind='crowd', K=3, gap=3.0, dmin=1.0, ratio=1, snr=14), 3, methods=('dual', 'greedy'))
    assert seen == [6, 6]


def test_scene_pins_cluster_and_shuffle():
    c, a, s, noise = X.scene(np.random.default_rng(3), dict(kind='cluster', K=4, dmin=2.2, ratio=5, snr=8, wr=1.5))
    assert np.allclose(c, [34.4452, 14.4068, 9.0405, 26.6668], atol=1e-4) and np.allclose(a, [1.8005, 1.077, 2.7164, 1.0], atol=1e-4)
    assert np.allclose(s, [0.9156, 1.0137, 0.9453, 1.0729], atol=1e-4) and np.isclose(noise, 0.24378, atol=1e-5)
    assert any(not np.all(np.diff(X.scene(np.random.default_rng(k), dict(kind='cluster', K=4, dmin=2.2, ratio=5, snr=8))[0]) > 0) for k in range(10))   # the centres are shuffled, not sorted


def test_trial_pins_the_natural_rules_and_the_agreement_flags():
    o, ag = X.trial_full(dict(kind='pair', d=2.0, ratio=1, snr=8), 0, methods=('dual', 'greedy'))
    n = X.RULES.index('nat')
    assert o[0, n, 0] == 1 and o[1, n, 0] == 0                                    # the dual's stopping time keeps the pair, the entrant significance of the greedy order does not
    o, ag = X.trial_full(dict(kind='pair', d=1.2, ratio=1, snr=12), 9, methods=('dual', 'greedy'))
    assert list(ag) == [0, 1]                                                      # the BIC blocks differ, the penalty-16 blocks agree


def test_apply_rules_disc_oracle_and_path_rms():
    c_true = np.array([10.0, 11.0])
    B = planted_path_blocks(11.0)
    out = X.apply_rules(np.zeros(X.NS), B, c_true, 1.0, None, bw=False)
    assert out[X.RULES.index('disc'), 0] == 0 and out[X.RULES.index('bic'), 0] == 1          # the discrepancy rule stops at one unit (chi^2 / N within the noise)
    good = dict(c=np.array([10.0, 11.0]), s=np.ones(2), w=np.ones(2), rss=50.0)
    good2 = dict(c=np.array([10.3, 11.0]), s=np.ones(2), w=np.ones(2), rss=40.0)
    bad = dict(c=np.array([10.0, 14.0]), s=np.ones(2), w=np.ones(2), rss=60.0)
    B2 = [dict(c=np.zeros(0), s=np.zeros(0), w=np.zeros(0), rss=200.0), dict(c=np.array([10.0]), s=np.ones(1), w=np.ones(1), rss=100.0), bad, good2, good]
    o = X.apply_rules(np.zeros(X.NS), B2, c_true, 1.0, None, bw=False)
    assert o[X.RULES.index('oracleK'), 0] == 0                                     # the FIRST block with K units is the bad one
    assert o[X.RULES.index('path'), 0] == 1
    assert np.isclose(o[X.RULES.index('path'), 1], 0.0, atol=1e-12)               # the best rms over the exact blocks (the exact one, rms 0), not the worst (0.3 / sqrt 2)


def test_sfw_ladder_and_significance_threshold():
    assert np.allclose(X.LAM_REL, [0.5, 0.2407, 0.1159, 0.0558, 0.0269, 0.0129, 0.0062, 0.003], atol=1e-4) and len(X.LAM_REL) == 8
    rng = np.random.default_rng(2)
    y = X.gauss([20.0, 30.0]) @ np.array([3.0, 2.0]) + 0.2 * rng.standard_normal(X.NS)
    B = X.sfw_blocks(y, 0.2)
    assert [len(b['w']) for b in B] == [0, 2, 2, 2, 2, 2, 2, 2, 3]
    assert np.allclose([b['rss'] for b in B][:4], [49.84, 19.61, 7.26, 4.4], atol=0.02)
    for b in B:
        r = y - (X.gauss(b['c']) @ b['w'] if len(b['w']) else 0.0)
        assert np.isclose(float(r @ r), b['rss'], rtol=1e-9)
    weak = (X.gauss([20.0]) * 0.5)[:, 0]
    assert [len(b['w']) for b in X.sfw_blocks(weak, 0.1)] == [0, 1, 1, 1, 1, 1, 1, 1, 1]          # significance 0.5 * 1.88 / 0.1 = 9.4 >= 3: kept
    assert all(len(b['w']) == 0 for b in X.sfw_blocks(weak, 0.5))                                  # 1.9 < 3: dropped


def test_part_2d_jobs_and_seeds(monkeypatch):
    import exp_order_pairs as OP
    seen = []
    monkeypatch.setattr(OP, '_init', lambda: None)
    monkeypatch.setattr(OP, 'trial', lambda a: (seen.append(a), [(1, 2), (2, 2)])[1])
    monkeypatch.setattr(X, 'TWOD_SCENES', 3)
    out = X.part_2d()
    assert out.shape == (4, 2, 3, 2) and (out[..., 0] == 1).all() and (out[..., 1] == 2).all()
    assert len(seen) == 24 and len(set(a[3] for a in seen)) == 24
    assert seen[0] == (1.2, 1.0, 0, 9000000, ('bic', 4.0), None) and seen[3] == (1.2, 1.0, 1, 9001000, ('bic', 4.0), None) and seen[-1] == (3.0, 1.0, 1, 9007002, ('bic', 4.0), None)


def test_stored_part_t_recomputes():
    raw = os.path.join(HERE, '..', 'results', 'synthetic', 'dualregime_raw.npz')
    if not os.path.exists(raw):
        pytest.skip('no stored raw yet')
    d = X.load_raw(raw)
    if 'T' in d:
        assert np.array_equal(X.part_t(), d['T'])


# ---- the wording round after the referee report C+D ---------------------------------------------------------------------------------------------------------------
def test_comparison_status_marks_what_may_not_be_quoted():
    assert X.comparison_status('P', '1') == 'order effect, everything else equal; backward step capped at K_true + 1 units'
    assert X.comparison_status('P', '2') == ''
    assert X.comparison_status('P', '3') == 'backward step capped at K_true + 1 units'
    assert X.comparison_status('P', '4') == 'different algorithm (sliding Frank-Wolfe), not an order comparison'
    assert X.comparison_status('P', '5') == 'different algorithm (sliding Frank-Wolfe), not an order comparison; backward step capped at K_true + 1 units'
    assert X.comparison_status('P', '9') == 'oracle K (order only, uses the true number of units)'
    assert X.comparison_status('P', '10') == X.comparison_status('Q', '10') == X.comparison_status('W', '10') == 'not an order comparison (different stopping rules), not to be quoted'
    assert X.comparison_status('W', '4').endswith('misspecified here (fixed width), not to be quoted') and X.comparison_status('W', '8').endswith('misspecified here (fixed width), not to be quoted')
    assert 'misspecified' not in X.comparison_status('P', '8') and 'misspecified' not in X.comparison_status('Q', '4')
    assert X.comparison_status('W', '6') == 'penalty-16 sensitivity in the regime W rows W4-W8, not to be quoted' and 'not to be quoted' not in X.comparison_status('P', '6')
    assert X.comparison_status('W', '7').startswith('penalty-16 sensitivity') and X.comparison_status('W', '7').endswith('backward step capped at K_true + 1 units')


def planted_pair_raw(n=200):
    fam = np.array(json.dumps(dict(kind='pair', d=1.6, ratio=1, snr=10)))
    E = np.zeros((n, 3, 8, 2))
    for r in ('bic', 'bicbw', 'pen', 'penbw', 'oracleK', 'nat'):
        E[:, 0, X.RULES.index(r), 0] = 1
    E[:100, 1, X.RULES.index('bic'), 0] = 1
    E[:100, 1, X.RULES.index('bicbw'), 0] = 1
    E[:100, 2, X.RULES.index('bic'), 0] = 1
    E[:, :, :, 1] = 0.1
    return {'P|0|E': E, 'P|0|fam': fam, 'P|0|widths': np.array('null')}


def test_the_pooled_table_has_a_status_column_and_the_numbers_are_returned():
    L, res = X.test_lines(planted_pair_raw())
    head = [l for l in L if l.startswith('| id |')][0]
    assert head.endswith('| verdict | status |')
    assert [l for l in L if l.startswith('| P1 |')][0].endswith('| win | order effect, everything else equal; backward step capped at K_true + 1 units |')
    assert [l for l in L if l.startswith('| P10 |')][0].endswith('not an order comparison (different stopping rules), not to be quoted |')
    num = res['P']['num']
    assert sorted(num) == sorted(str(k) for k in range(1, 11)) and num['1'][0] == pytest.approx(50.0) and num['1'][1] <= num['1'][0] <= num['1'][2]


def test_the_summary_names_the_failing_part_of_the_registered_rule():
    base = {str(k): 'win' for k in range(1, 11)}
    raw = {}
    def txt(**over):
        res = {'P': dict(base, regime='no', **over)}
        return '\n'.join(X.summary_lines(raw, res, None))
    assert 'the failing part is (ii) alone' in txt(**{'4': 'loss'})
    assert 'the failing part is (i) alone' in txt(**{'3': 'tie'})
    assert 'neither (i) nor (ii) is met' in txt(**{'3': 'tie', '4': 'loss'})
    assert '**yes**, met' in '\n'.join(X.summary_lines(raw, {'P': dict(base, regime='yes')}, None))
    w = '\n'.join(X.summary_lines(raw, {'W': dict(base, regime='no', **{'3': 'loss'})}, None))
    assert 'Regime W (unequal widths): (i)' in w and '(W4-W8) is not an order comparison and, with its fixed width, misspecified here: not to be quoted.' in w
    assert 'Registered plan rule (W1, W3 all wins): **no**, the failing part is (i) alone' in w and 'first block with K units (W9, oracle K) a win' in w
    w2 = '\n'.join(X.summary_lines(raw, {'W': dict(base, regime='no', **{'3': 'loss', '4': 'loss'})}, None))
    assert 'the failing part is (i) alone' in w2                                      # the fixed-width sFW comparator is never counted in W


def test_the_equivalence_line_states_its_scope_and_the_two_dimensional_loss_at_three_widths():
    tw = np.ones((4, 2, 100, 2), int) * 2
    tw[3, :, :, 1] = 1                                                            # 3 widths: the greedy order resolves; 2.5 widths: neither does
    t3 = X.twod_sep(tw, 3.0)
    assert t3 == pytest.approx((-100.0, -100.0, -100.0)) and X.twod_sep(tw, 3.5) is None
    assert X.twod_sep(tw, 1.2) == pytest.approx((0.0, 0.0, 0.0)) and X.twod_sep(tw, 2.5) == pytest.approx((0.0, 0.0, 0.0)) and X.twod_sep(tw, 1.6) == pytest.approx((0.0, 0.0, 0.0))
    txt = '\n'.join(X.summary_lines({'twod': tw}, {'equiv': 'yes', 'equiv_n': (797, 800)}, None))
    assert 'equivalent (plan rule) **yes** (797 of 800 scenes).  Scope: one-dimensional pairs with two features' in txt
    assert 'in two dimensions at 3 widths the dual minus greedy difference is -100.0 [-100.0, -100.0] points (a loss) and 3.5 widths and more were not run' in txt
    L, out = X.twod_lines(tw)
    assert any(l.startswith('Separation 3 widths alone (both photon levels pooled): dual minus greedy -100.0 [-100.0, -100.0] points: **loss**.') for l in L)


def test_the_conclusions_state_the_scope_the_hypothesis_status_the_reduction_and_the_caveats():
    base = {str(k): 'tie' for k in range(1, 11)}
    num = lambda a, b: {'1': (a, a - 1, a + 1), '2': (b, b - 1, b + 1)}                       # noqa: E731
    res = {'P': dict(base, regime='no', **{'1': 'win', '2': 'win', 'num': num(2.1, 7.8)}), 'Q': dict(base, regime='no', **{'1': 'win', '2': 'win', 'num': num(3.4, 8.0)}), 'equiv': 'yes', 'equiv_n': (797, 800), 'twod': dict(close='win', mid='loss')}
    txt = '\n'.join(X.conclusion_lines(res))
    assert '(797 of 800 scenes returned the same block) and for the FINAL BIC blocks, not for the entry sequences' in txt and 'not for images' in txt and 'is a loss at the mid separations (2.5, 3 widths)' in txt
    assert '5. **Mechanism (a hypothesis, not a finding).**' in txt and 'CONSISTENT WITH' in txt and 'this was not tested' in txt and 'Where the order matters it is memory' not in txt
    assert 'The backward step reduces the lead over the bare greedy order but does not remove it: from +7.8 (P2) and +8.0 (Q2) points to +2.1 (P1) and +3.4 (Q1) points, that is to about 2-3 points (P1 win, Q1 win).' in txt
    assert 'is removed' not in txt
    assert '8. **Design caveats' in txt and 'K_true + 1 units' in txt and '(comparisons 1, 3, 5, 7)' in txt and 'inflated by about 0 to 3 points' in txt and 'not recomputed here' in txt
    assert 'a comparison with a different algorithm, not an order comparison, and the dual is a tie there (P4)' in txt
    assert 'The backward step reduces' not in '\n'.join(X.conclusion_lines({'P': dict(base, regime='no'), 'equiv': 'no'}))
    only_p = {'P': dict(base, regime='no', num=num(2.1, 7.8)), 'Q': dict(base, regime='no'), 'equiv': 'no'}              # the reduction is stated only when P and Q both have their numbers
    assert 'The backward step reduces' not in '\n'.join(X.conclusion_lines(only_p))
