"""Tests of experiments/exp_baselines_sanity.py: the fair count (merging), the floored BIC, the three baselines' unit extraction on hand-made raws and on planted cases with a known answer,
the sanity targets, the job lists, the development tuning, the resumable runner, failed trials, a tiny end-to-end run and the byte-identical report."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_baselines_sanity as E      # noqa: E402

X, B = E.X, E.B
N = X.NS


# ---------------------------------------------------------------------------------------------------------------------------
# merging and the information criterion
# ---------------------------------------------------------------------------------------------------------------------------
def test_merge_units_clusters_masses_positions_and_the_mass_threshold():
    W, C = E.merge_units([0.04, 0.5, -0.2, 0.5], [0.9, 0.31, 0.60, 0.30], 0.12, 0.1)          # unsorted on purpose
    assert np.allclose(W, [1.0, -0.2], atol=0, rtol=1e-12) and np.allclose(C, [0.305, 0.60], atol=1e-12, rtol=0)     # 0.04 < 0.1 x 1.0 dropped
    W, C = E.merge_units([0.04, 0.5, -0.2, 0.5], [0.9, 0.31, 0.60, 0.30], 0.12, 0.0)
    assert len(W) == 3 and np.isclose(W[2], 0.04, atol=0, rtol=1e-12)                       # no mass threshold: kept
    W, C = E.merge_units([3.0, 1.0], [0.0, 1.0], 1.0, 0.0)                                  # the gap equals the link: one cluster; position |w|-weighted
    assert len(W) == 1 and np.isclose(W[0], 4.0, atol=0, rtol=1e-14) and np.isclose(C[0], 0.25, atol=1e-14, rtol=0)
    W, C = E.merge_units([3.0, 1.0], [0.0, 1.0001], 1.0, 0.0)                               # just beyond the link: two
    assert len(W) == 2 and np.allclose(W, [3.0, 1.0], atol=0, rtol=0)


def test_merge_units_mass_threshold_is_inclusive_and_zero_units_do_not_bridge():
    W, C = E.merge_units([1.0, 0.1], [0.0, 1.0], 0.12, 0.1)                                 # a mass exactly at thm x the largest is kept
    assert len(W) == 2 and np.allclose(W, [1.0, 0.1], atol=0, rtol=0)
    W, C = E.merge_units([1.0, 0.0, 1.0], [0.0, 0.1, 0.2], 0.12, 0.0)                       # a zero-coefficient unit in the middle does not join two spikes 0.2 apart
    assert len(W) == 2 and np.allclose(C, [0.0, 0.2], atol=0, rtol=0)


def test_merge_units_mixed_signs_cancel_zeros_are_ignored_and_empty_is_empty():
    W, C = E.merge_units([1.0, -1.0], [0.3, 0.31], 0.12, 0.0)
    assert len(W) == 0 and len(C) == 0                                                      # zero total mass: no spike
    W, C = E.merge_units([0.0, 2.0, 0.0], [0.0, 0.5, 0.02], 0.12, 0.0)                     # zero coefficients neither join nor move a cluster
    assert len(W) == 1 and np.isclose(C[0], 0.5, atol=0, rtol=1e-14)
    W, C = E.merge_units(np.zeros(3), np.array([0.0, 0.1, 0.2]), 0.12, 0.1)
    assert len(W) == 0
    W, C = E.merge_units([2.0, -2.0], [0.3, 0.6], 0.12, 0.1)                                # opposite signs far apart: two spikes
    assert len(W) == 2 and np.allclose(W, [2.0, -2.0], atol=0, rtol=0)


def test_ic_value_floor_and_select_multiplier():
    assert np.isclose(E.ic_value(0.0, N, 3), N * np.log(E.FLOOR ** 2) + 3 * np.log(N), atol=0, rtol=1e-14)       # the floor binds at a vanishing residual
    assert np.isclose(E.ic_value(5.0, N, 3), N * np.log(5.0 / N) + 3 * np.log(N), atol=0, rtol=1e-14)           # and not otherwise
    assert E.ic_select([1.0, 0.0, 0.0, 0.0], 1, N) == 1                                     # noiseless: stops at the exact block, no extra unit fits "zero"
    rss = [10.0, 1.0, 0.9]
    assert E.ic_select(rss, 1, N, 1.0) == 2 and E.ic_select(rss, 1, N, 2.0) == 1             # the multiplier scales the penalty (gain of the 2nd unit 6.4 vs 4.1 and 8.2)
    assert E.ic_select(rss, 2, N, 1.0) == 1                                                  # two parameters per unit: 8.2 > 6.4


# ---------------------------------------------------------------------------------------------------------------------------
# the baselines on planted cases with a known answer
# ---------------------------------------------------------------------------------------------------------------------------
def _step_raw(ct, at, s=0.0, seed=0):
    Kd, grid = X.setup()
    f = Kd.synth(np.array(ct)[:, None], np.array(at))
    y = f + (s * np.random.default_rng(seed).standard_normal(Kd.N) if s > 0 else 0.0)
    S = B.stepwise_bic(Kd, y, grid, max_units=E.MAXU)
    sw = np.full((E.MAXU, E.MAXU), np.nan)
    for k in range(1, E.MAXU + 1):
        sw[k - 1, :k] = S['blocks'][k]['w']
    return Kd, y, np.array(S['order'], int), sw, np.array([b['rss'] for b in S['blocks']])


def test_forward_selection_is_exact_on_the_grid_and_needs_several_units_off_it_but_merges_to_one():
    Kd, y, idx, sw, rss = _step_raw([-0.6, 0.6], [1.0, -0.8])                               # on the grid, noiseless
    kh, W, C = E.units_step(idx, sw, rss, {'c': 1}, False, N)
    assert kh == 2 and np.allclose(sorted(C), [-0.6, 0.6], atol=1e-12) and np.allclose(sorted(np.abs(W)), [0.8, 1.0], atol=1e-6)
    kh1, W1, C1 = E.units_step(idx, sw, rss, {'c': 1, 'link': 0.12, 'thm': 0.1}, True, N)
    assert kh1 == 2
    Kd, y, idx, sw, rss = _step_raw([0.35], [1.0])                                          # half a grid step off, noiseless
    kh, W, C = E.units_step(idx, sw, rss, {'c': 1}, False, N)
    assert kh >= 3                                                                          # model class: an off-grid bump needs several grid kernels
    kh1, W1, C1 = E.units_step(idx, sw, rss, {'c': 1, 'link': 0.12, 'thm': 0.1}, True, N)
    assert kh1 == 1 and abs(C1[0] - 0.35) < 0.02 and abs(W1[0] - 1.0) < 0.05               # merged: one feature, in the right place, with the right mass
    kh2, W2, C2 = E.units_step(idx, sw, rss, {'c': 1, 'link': 0.06, 'thm': 0.1}, True, N)     # a link below the grid step merges nothing: V1 = V0 up to the mass threshold
    assert kh2 >= 3


def test_slide_merges_units_that_end_up_within_a_width_but_not_further():
    Kd, grid = X.setup()
    y = Kd.synth(np.array([[0.30], [0.39]]), np.array([1.0, 1.0]))                           # two features 0.09 apart
    kh, W, C = E.slide_units(Kd, y, np.array([1.0, 1.0]), np.array([0.30, 0.39]))
    assert kh == 1 and np.isclose(W[0], 2.0, atol=1e-6, rtol=0)                              # closer than the link 0.12: one spike
    y = Kd.synth(np.array([[0.30], [0.50]]), np.array([1.0, 1.0]))
    kh, W, C = E.slide_units(Kd, y, np.array([1.0, 1.0]), np.array([0.30, 0.50]))
    assert kh == 2 and np.allclose(sorted(C), [0.30, 0.50], atol=1e-6, rtol=0)


def test_slide_refines_the_merged_position():
    Kd, y, idx, sw, rss = _step_raw([0.35], [1.0])
    kh1, W1, C1 = E.units_step(idx, sw, rss, {'c': 1, 'link': 0.12, 'thm': 0.1}, True, N)
    kh2, W2, C2 = E.slide_units(Kd, y, W1, C1)
    assert kh2 == 1 and abs(C2[0] - 0.35) < 1e-6 and abs(W2[0] - 1.0) < 1e-6              # noiseless: the slide lands on the truth
    assert E.slide_units(Kd, y, np.zeros(0), np.zeros(0))[0] == 0
    W3, C3 = np.array([0.5, 0.5]), np.array([0.34, 0.36])                                  # two clusters that came together after the slide are merged again
    kh3, W3, C3 = E.slide_units(Kd, y, W3, C3)
    assert kh3 == 1 and abs(C3[0] - 0.35) < 1e-3


def test_units_pen_rule_theta_and_merge_on_a_hand_made_net():
    w = np.zeros((E.NPEN, E.ND))
    c = np.tile(E.GRID, (E.NPEN, 1)).copy()
    w[1, [3, 4, 5, 12]] = [0.6, 0.35, 0.05, 1e-4]                                           # penalty 1: three units close together and a speck
    w[2, [3]] = 0.9                                                                         # penalty 2: one unit
    c[1, 3], c[1, 4], c[1, 5] = 0.30, 0.30, 0.33
    rss = np.array([2.0, 0.05, 0.06, 1.0, 1.0, 1.0, 1.0])
    kh, W, C = E.units_pen(w, c, rss, {'rule': ('bic', 1), 'theta': 1e-6}, False, N)         # BIC: penalty 1 (rss 0.05, 4 nonzero) beats penalty 2 (0.06, 1 nonzero)? check by hand below
    vals = [E.ic_value(rss[q], N, 2 * int((np.abs(w[q]) > 1e-6 * max(np.abs(w[q]).max(), 1e-12)).sum())) for q in range(E.NPEN)]
    j = int(np.argmin(vals))
    assert kh == int((np.abs(w[j]) > 1e-6 * max(np.abs(w[j]).max(), 1e-12)).sum())
    rss2 = np.array([2.0, 0.01, 0.06, 1.0, 1.0, 1.0, 1.0])
    assert E.units_pen(w, c, rss2, {'rule': ('bic', 1), 'theta': 1e-6}, False, N)[0] == 4        # BIC x1: the four-unit fit (rss 0.01) beats the one-unit fit (0.06) ...
    assert E.units_pen(w, c, rss2, {'rule': ('bic', 8), 'theta': 1e-6}, False, N)[0] == 1        # ... BIC x8: the other way round
    assert E.units_pen(w, c, rss2, {'rule': ('bic', 8), 'link': 0.12, 'thm': 0.1}, True, N)[0] == 1
    kh, W, C = E.units_pen(w, c, rss, {'rule': ('lam', 1), 'theta': 1e-6}, False, N)         # fixed penalty 1: four nonzero coefficients
    assert kh == 4 and np.allclose(sorted(W), [1e-4, 0.05, 0.35, 0.6], atol=0, rtol=0)
    kh, W, C = E.units_pen(w, c, rss, {'rule': ('lam', 1), 'theta': 1e-3}, False, N)         # the count threshold removes the speck (1e-4 < 1e-3 x 0.6)
    assert kh == 3
    kh, W, C = E.units_pen(w, c, rss, {'rule': ('lam', 1), 'theta': 0.1}, False, N)          # and 0.05 < 0.1 x 0.6
    assert kh == 2
    kh, W, C = E.units_pen(w, c, rss, {'rule': ('lam', 1), 'link': 0.12, 'thm': 0.1}, True, N)    # merged: 0.3, 0.3, 0.33 and the speck at 0.2 (0.1 away) -> one cluster of mass 1.0001
    assert kh == 1 and np.isclose(W[0], 1.0001, atol=0, rtol=1e-12) and abs(C[0] - (0.6 * 0.30 + 0.35 * 0.30 + 0.05 * 0.33 + 1e-4 * 0.2) / 1.0001) < 1e-12
    kh, W, C = E.units_pen(w, c, rss, {'rule': ('lam', 0), 'theta': 1e-6}, False, N)         # an all-zero net: nothing
    assert kh == 0 and len(W) == 0


def test_units_prune_scale_threshold_and_merge_on_a_hand_made_net():
    w = np.zeros((2, E.ND))
    c = np.tile(E.GRID, (2, 1)).copy()
    w[0, :6] = [1.0, 0.5, 0.2, -0.3, 0.09, 0.0]
    w[1, :2] = [0.4, 0.4]
    kh, W, C = E.units_prune(w, c, {'scale': 0, 'thr': 0.1}, False)                           # |w| >= 0.1: 1.0, 0.5, 0.2, -0.3 (0.09 and 0 removed)
    assert kh == 4 and np.allclose(sorted(np.abs(W)), [0.2, 0.3, 0.5, 1.0], atol=0, rtol=0)
    kh, W, C = E.units_prune(w, c, {'scale': 0, 'thr': 0.3}, False)                           # 1.0, 0.5, -0.3 (>= 0.3 is inclusive)
    assert kh == 3
    kh, W, C = E.units_prune(w, c, {'scale': 1, 'thr': 0.1}, False)                           # the other scale's net
    assert kh == 2 and np.allclose(W, [0.4, 0.4], atol=0, rtol=0)
    kh, W, C = E.units_prune(w, c, {'scale': 0, 'link': 0.12, 'thr': 0.1}, True)              # grid 0.1: units 0..4 chain into one cluster of mass 1.0+0.5+0.2-0.3+0.09 = 1.49
    assert kh == 1 and np.isclose(W[0], 1.49, atol=0, rtol=1e-12)
    kh, W, C = E.units_prune(w, c, {'scale': 0, 'link': 0.06, 'thr': 0.1}, True)              # no merging below the grid step; masses 1.0, .5, .2, -.3 above 0.1; .09 below
    assert kh == 4


def test_config_grids_start_at_the_comparisons_setting_and_have_the_declared_sizes():
    sizes = {('step', False): len(E.C_STEP), ('step', True): len(E.C_STEP) * len(E.LINKS) * len(E.THMS), ('pen', False): len(E.PEN_RULES) * len(E.THETAS),
             ('pen', True): len(E.PEN_RULES) * len(E.LINKS) * len(E.THMS), ('prune', False): 2 * len(E.THRS), ('prune', True): 2 * len(E.LINKS) * len(E.THRS)}
    assert sizes == {('step', False): 8, ('step', True): 160, ('pen', False): 66, ('pen', True): 220, ('prune', False): 18, ('prune', True): 90}
    for (m, merged), n in sizes.items():
        g = E.config_grid(m, merged)
        assert len(g) == n and g[0] == E.DEFAULT[m][1 if merged else 0] and g.count(g[0]) == 1
    assert E.DEFAULT['step'][0] == {'c': 1} and E.DEFAULT['pen'][0] == {'rule': ('bic', 1), 'theta': 1e-6} and E.DEFAULT['prune'][0] == {'scale': 0, 'thr': 0.1}
    assert E.DEFAULT['step'][1] == {'c': 1, 'link': 0.12, 'thm': 0.1} and E.DEFAULT['pen'][1] == {'rule': ('bic', 1), 'link': 0.12, 'thm': 0.1}
    assert E.DEFAULT['prune'][1] == {'scale': 0, 'link': 0.12, 'thr': 0.1}
    assert E.SCALES == (0.3, 0.1) and E.STEPS_PEN == 2500 and E.STEPS_PRUNE == 5000 and len(E.PEN_LAMS) == 7 and E.PEN_LAMS == X.PEN_LAMS
    assert E.cfg_str('step', {'c': 2}) == 'BIC x2' and 'fixed lam' in E.cfg_str('pen', {'rule': ('lam', 2), 'theta': 0.1}) and E.cfg_str('prune', {'scale': 1, 'thr': 0.5}).startswith('scale 0.1')


class _FakeStage:
    """A stand-in for a development Stage: the exact-order rate is a given function of the configuration."""
    bad = np.zeros(1, bool)
    K = np.ones(1, int)

    def __init__(self, f):
        self.f = f

    def score(self, method, merged, cfg):
        return np.array([[1.0 if self.f(method, merged, cfg) else 0.0, 1.0, 0.0, 0.0]])


def test_exact_rate_leaves_failed_trials_out():
    class St:
        bad = np.array([False, True, False])
        K = np.array([1, 1, 2])

        def score(self, method, merged, cfg):
            return np.array([[1.0, 1, 0, 0], [np.nan] * 4, [3.0, 2, 0, 0]])
    assert E.exact_rate(St(), 'step', False, {'c': 1}) == 0.5


def test_tune_takes_the_best_and_the_first_on_ties():
    st = _FakeStage(lambda m, mg, cfg: m == 'step' and not mg and cfg['c'] == 4)
    t = E.tune(st)
    assert t[('step', False)][0] == {'c': 4} and t[('step', False)][1] == 1.0 and t[('step', False)][2] == 0.0
    assert t[('pen', False)][0] == E.DEFAULT['pen'][0] and t[('pen', False)][1] == 0.0                 # no setting better: the comparison's own (the first)
    assert t[('step', True)][0] == E.DEFAULT['step'][1] and t[('prune', True)][0] == E.DEFAULT['prune'][1]
    st = _FakeStage(lambda m, mg, cfg: m == 'prune' and not mg and cfg['thr'] >= 0.5)
    t = E.tune(st)
    assert t[('prune', False)][0] == {'scale': 0, 'thr': 0.5}                                          # the first (listed) of the maximisers


# ---------------------------------------------------------------------------------------------------------------------------
# the sanity targets and the job lists
# ---------------------------------------------------------------------------------------------------------------------------
def test_sanity_targets_have_the_declared_placement_separation_and_amplitudes():
    for seed in range(40):
        c, a = E.sanity_target(1, 0, seed)
        assert len(c) == 1 and abs(abs(a[0]) - 1.0) < 1e-12 and abs(c[0] / 0.1 - round(c[0] / 0.1)) < 1e-9 and abs(c[0]) <= 0.8 + 1e-9
        c, a = E.sanity_target(1, 1, seed)
        assert abs((c[0] - 0.05) / 0.1 - round((c[0] - 0.05) / 0.1)) < 1e-9 and abs(c[0]) <= 0.85 + 1e-9
        c, a = E.sanity_target(1, 2, seed)
        assert abs(c[0]) <= 0.85
        for plc in range(3):
            c, a = E.sanity_target(2, plc, seed)
            assert len(c) == 2 and c[1] - c[0] >= E.SAN_SEP - 1e-9 and c.min() >= -0.85 - 1e-9 and c.max() <= 0.85 + 1e-9
            assert abs(np.abs(a).max() - 1.0) < 1e-12 and 0.5 <= np.abs(a).min() <= 1.0 and np.abs(a).min() < 1.0
            if plc == 0:
                assert np.allclose(c / 0.1, np.round(c / 0.1), atol=1e-9)
            if plc == 1:
                assert np.allclose((c - 0.05) / 0.1, np.round((c - 0.05) / 0.1), atol=1e-9)
    assert E.SAN_SEP == pytest.approx(1.2) and E.sanity_target(2, 2, 5)[0].tolist() == E.sanity_target(2, 2, 5)[0].tolist()
    signs = {tuple(np.sign(E.sanity_target(2, 0, s)[1])) for s in range(60)}
    assert len(signs) == 4                                                                     # all four sign patterns occur
    big = {int(np.argmax(np.abs(E.sanity_target(2, 0, s)[1]))) for s in range(60)}
    assert big == {0, 1}                                                                       # the largest amplitude is on either feature
    pos = {int(np.sign(E.sanity_target(1, 2, s)[1][0])) for s in range(40)}
    assert pos == {-1, 1}


def test_job_lists_match_the_designs():
    sj, dj, mj = E.sanity_jobs(), E.dev_jobs(), E.main_jobs()
    assert len(sj) == 6 * E.SAN_T * (1 + 2 + 2) == 600 and len(dj) == 6 * 8 * 4 == 192 and len(mj) == 432
    assert mj == [(2, K, ni, 1000 * K + t, d, -1, X.NOISE[ni]) for K in X.KS for t in range(6) for ni in range(4) for d in range(3)]
    assert set(j[3] for j in dj).isdisjoint(set(j[3] for j in mj)) and set(j[3] for j in dj) == {1000 * K + 100 + t for K in X.KS for t in range(8)}
    assert all(j[4] == 0 for j in dj) and all(j[4] == 0 for j in sj if j[6] == 0) and {j[4] for j in sj if j[6] > 0} == {0, 1}
    assert {j[6] for j in sj} == {0.0, 0.001, 0.01} and {j[1] for j in sj} == {1, 2} and {j[5] for j in sj} == {0, 1, 2}
    assert len({j[3] for j in sj}) == 6 * E.SAN_T                                              # one target seed per (cell, t), shared by the noise levels
    assert E.main_jobs(1, 1)[0] == (2, 1, 0, 1000, 0, -1, 0.01) and len(E.main_jobs(1, 1)) == 24
    assert sj[0] == (0, 1, 0, 500000, 0, 0, 0.0) and sj[1] == (0, 1, 1, 500000, 0, 0, 0.001) and sj[2] == (0, 1, 1, 500000, 1, 0, 0.001) and sj[4] == (0, 1, 2, 500000, 1, 0, 0.01)
    assert sj[5][3] == 500001 and sj[5 * E.SAN_T][3] == 501000 and sj[-1] == (0, 2, 2, 505000 + E.SAN_T - 1, 1, 2, 0.01) and sj[5 * E.SAN_T * 3][3] == 503000
    assert dj[0] == (1, 1, 0, 1100, 0, -1, 0.01) and dj[5] == (1, 1, 1, 1101, 0, -1, 0.03)


def test_cluster_rate_ci_is_over_targets():
    x = np.array([1, 1, 1, 0, 0, 0, 1, 0], float)
    tid = np.array([7, 7, 7, 8, 8, 8, 9, 9])
    mu, lo, hi = E.cluster_rate_ci(x, tid)
    assert np.isclose(mu, (1.0 + 0.0 + 0.5) / 3, atol=0, rtol=1e-14) and 0.0 <= lo <= mu <= hi <= 1.0 and (lo, hi) == E.cluster_rate_ci(x, tid)[1:]
    means = np.array([1.0, 0.0, 0.5])
    boot = means[np.random.default_rng(0).integers(0, 3, size=(2000, 3))].mean(axis=1)             # 2000 resamples of the three targets, seed 0
    assert np.isclose(lo, np.percentile(boot, 2.5), atol=0, rtol=1e-14) and np.isclose(hi, np.percentile(boot, 97.5), atol=0, rtol=1e-14)
    assert E.cluster_rate_ci(np.ones(6), np.arange(6) // 2) == (1.0, 1.0, 1.0)


# ---------------------------------------------------------------------------------------------------------------------------
# the implementation checks: the gradient of the training the baselines use
# ---------------------------------------------------------------------------------------------------------------------------
def test_checks_gradient_matches_finite_differences_and_the_planted_solution_is_a_fixed_point(monkeypatch):
    monkeypatch.setattr(E, 'STEPS_PEN', 10)
    monkeypatch.setattr(E, 'STEPS_PRUNE', 10)
    out = E.checks()
    assert out['chk_grad'].shape == (2,) and np.all(out['chk_grad'] < 1e-5)                    # the update is the gradient of the loss, in w and in c
    assert np.all(out['chk_fixed'] < 1e-12)
    assert out['chk_conv'].shape == (3, 4) and out['chk_prune'].shape == (3, 4)
    assert np.all(out['chk_conv'][:, 1] >= 1) and np.all(out['chk_conv'][:, 0] > 0)


# ---------------------------------------------------------------------------------------------------------------------------
# runner: failed trials, resume, end-to-end, report
# ---------------------------------------------------------------------------------------------------------------------------
def test_trial_data_truth_and_networks_are_those_of_the_stated_design(monkeypatch):
    monkeypatch.setattr(E, 'STEPS_PEN', 40)
    monkeypatch.setattr(E, 'STEPS_PRUNE', 40)
    Kd, grid = X.setup()
    # a sanity trial: the noise is s x the largest amplitude times the standard normals of seed 1000 seed + draw + 1
    r = E._trial((0, 1, 2, 500002, 1, 2, 0.01))
    ct, at = E.sanity_target(1, 2, 500002)
    f = Kd.synth(ct[:, None], at)
    y = f + 0.01 * np.abs(at).max() * np.random.default_rng(500002 * 1000 + 1 + 1).standard_normal(Kd.N)
    assert np.isclose(r['step_rss'][0], y @ y, atol=0, rtol=1e-12) and r['ct'][0] == ct[0] and np.isnan(r['ct'][1]) and r['at'][0] == at[0]
    assert list(r['meta'][:6]) == [1, 0.01, 500002, 1, 2, 2] and np.isclose(r['meta'][6], 0.01 * abs(at[0]), atol=0, rtol=1e-14)
    net = B.penalised_net(Kd, y, grid, E.PEN_LAMS, steps=40, seed=500002)
    assert np.array_equal(r['pen_w'][3], net[3]['w']) and np.array_equal(r['pen_c'][3], net[3]['c'].ravel()) and np.isclose(r['pen_rss'][3], net[3]['rss'], atol=0, rtol=1e-12)
    for i, sc in enumerate(E.SCALES):
        Z = B.prune_net(Kd, y, grid, steps=40, seed=500002, scale=sc)
        assert np.array_equal(r['pr_w'][i], Z['w']) and np.array_equal(r['pr_c'][i], Z['c'].ravel())
    S = B.stepwise_bic(Kd, y, grid, max_units=E.MAXU)
    assert list(r['step_idx']) == list(S['order']) and np.array_equal(r['step_w'][2, :3], S['blocks'][3]['w']) and np.isnan(r['step_w'][2, 3])
    # noiseless: the data is the signal itself; controls found the planted feature; the secs are positive
    r = E._trial((0, 1, 0, 500000, 0, 0, 0.0))
    ct, at = E.sanity_target(1, 0, 500000)
    f = Kd.synth(ct[:, None], at)
    assert np.isclose(r['step_rss'][0], f @ f, atol=0, rtol=1e-12) and r['ctrl'][0, 0] == 1 and r['ctrl'][1, 0] == 1 and r['ctrl'][0, 1] == 1
    assert np.all(r['secs'] > 0) and r['secs'].shape == (5,)
    # a dev / main trial: the target is the comparison's, the noise draw its own
    r2 = E._trial((0, 1, 2, 500000, 0, 0, 0.01))                                               # noisy, on the grid: the positive controls (BIC with two parameters per unit) find one feature
    assert r2['ctrl'][0, 0] == 1 and r2['ctrl'][1, 0] == 1
    r = E._trial((2, 3, 1, 3002, 2, -1, 0.03))
    ct, at = X.make_target(np.random.default_rng(3002), 3)
    f = Kd.synth(ct[:, None], at)
    y = f + 0.03 * np.abs(at).max() * np.random.default_rng(3002 * 1000 + 2 + 1).standard_normal(Kd.N)
    assert np.isclose(r['step_rss'][0], y @ y, atol=0, rtol=1e-12) and np.array_equal(r['ct'][:3], ct) and np.isnan(r['ctrl']).all()
    for i, m in enumerate(('step', 'pen', 'prune')):                                           # V2 = the default merged units slid, each baseline in its own column
        if m == 'step':
            un = E.units_step(r['step_idx'], r['step_w'], r['step_rss'], E.DEFAULT['step'][1], True, Kd.N)
        elif m == 'pen':
            un = E.units_pen(r['pen_w'], r['pen_c'], r['pen_rss'], E.DEFAULT['pen'][1], True, Kd.N)
        else:
            un = E.units_prune(r['pr_w'], r['pr_c'], E.DEFAULT['prune'][1], True)
        kh, W, C = E.slide_units(Kd, y, un[1], un[2])
        assert np.allclose(r['v2'][i], X._score(kh, W, C[:, None], ct, at), atol=1e-12, rtol=0, equal_nan=True)


def test_main_stage_refuses_a_job_list_that_is_not_the_stored_comparisons(monkeypatch, tmp_path):
    jobs = E.main_jobs(1, 1)[:2]
    meta = np.array([[j[1], j[2], j[3] + (1 if i else 0), j[4], np.nan, 1.0] for i, j in enumerate(jobs)], float)      # the second seed differs
    np.savez_compressed(str(tmp_path / 'ref.npz'), res=np.zeros((2, 3, 4)), meta=meta, cols=np.array(['a', 'b', 'c']))
    monkeypatch.setattr(E, 'REF_RAW', str(tmp_path / 'ref.npz'))

    def fake(args):
        out = {k: np.full(v, 0.0) for k, v in E._shapes().items()}
        out['meta'] = np.array([args[1], args[6], args[3], args[4], args[5], args[2], 0.0])
        return out
    monkeypatch.setattr(E, 'trial', fake)
    with pytest.raises(AssertionError):
        E.run_stage('main', out_path=str(tmp_path / 'raw.npz'), jobs=jobs)


def test_failed_trial_is_recorded_as_nan_with_its_meta(monkeypatch):
    def boom(args):
        raise RuntimeError('planted')
    monkeypatch.setattr(E, '_trial', boom)
    out = E.trial((0, 2, 1, 500123, 1, 2, 0.001))
    assert set(out) == set(E._shapes()) and np.isnan(out['meta'][6]) and out['meta'][0] == 2 and out['meta'][2] == 500123 and out['meta'][4] == 2
    assert np.isnan(out['pen_w']).all() and out['pen_w'].shape == (E.NPEN, E.ND) and out['ctrl'].shape == (2, 4)


def test_resume_continues_a_partial_run_and_runs_each_job_once(monkeypatch, tmp_path):
    calls = []

    def fake(args):
        calls.append(args)
        out = {k: np.full(v, float(args[3])) for k, v in E._shapes().items()}
        out['meta'] = np.array([args[1], args[6], args[3], args[4], args[5], args[2], 0.0])
        return out
    monkeypatch.setattr(E, 'trial', fake)
    jobs = [(0, 1, 0, 500000 + i, 0, 0, 0.0) for i in range(5)]
    path = str(tmp_path / 'raw.npz')
    state = {'n': 0}
    real = np.savez_compressed

    def stop_after_first_chunk(p, **kw):
        real(p, **kw)
        state['n'] += 1
        if state['n'] == 1:
            raise KeyboardInterrupt
    monkeypatch.setattr(np, 'savez_compressed', stop_after_first_chunk)
    with pytest.raises(KeyboardInterrupt):
        E.run_stage('sanity', out_path=path, chunk=2, jobs=jobs)
    monkeypatch.setattr(np, 'savez_compressed', real)
    assert len(calls) == 2 and os.path.exists(path + '.sanity.part.npz')
    E.run_stage('sanity', out_path=path, chunk=2, jobs=jobs)
    assert [c[3] for c in calls] == [500000 + i for i in range(5)]                         # the first two not repeated
    z = np.load(path)
    assert z['sanity_meta'].shape == (5, 7) and z['sanity_pen_w'].shape == (5, E.NPEN, E.ND) and list(z['sanity_meta'][:, 2]) == [500000 + i for i in range(5)]
    assert not os.path.exists(path + '.sanity.part.npz') and float(z['sanity_wall']) >= 0


@pytest.fixture(scope='module')
def mini_raw(tmp_path_factory):
    """A raw of the four stages with a handful of real trials each (the networks trained for a few steps only)."""
    mp = pytest.MonkeyPatch()
    mp.setattr(E, 'STEPS_PEN', 30)
    mp.setattr(E, 'STEPS_PRUNE', 30)
    d = tmp_path_factory.mktemp('bs')
    path = str(d / 'raw.npz')
    sj = [(0, 1, 0, 500000, 0, 0, 0.0), (0, 1, 1, 500000, 0, 0, 0.001), (0, 1, 1, 500000, 1, 0, 0.001), (0, 1, 2, 500001, 0, 1, 0.01), (0, 2, 0, 503000, 0, 0, 0.0),
          (0, 2, 2, 504000, 1, 1, 0.01)]
    dj = [(1, K, ni, 1000 * K + 100, 0, -1, X.NOISE[ni]) for K in (1, 2) for ni in range(4)]
    mj = E.main_jobs(1, 1)[:8] + E.main_jobs(1, 1)[8:10]
    # a stored-comparison stand-in with the columns the report reads
    cols = [c for c, _ in E.REF_ROWS]
    res = np.zeros((len(E.main_jobs(1, 1)[:8]) + 2, len(cols), 4))
    meta = np.array([[j[1], j[2], j[3], j[4], np.nan, 1.0] for j in mj], float)
    res[:, :, 0] = meta[:, 0][:, None]
    res[:, :, 1] = meta[:, 0][:, None]
    res[0, 0, 0] = 3
    np.savez_compressed(str(d / 'ref.npz'), res=res, meta=meta, cols=np.array(cols))
    mp.setattr(E, 'REF_RAW', str(d / 'ref.npz'))
    E.run_checks(path)
    E.run_stage('sanity', out_path=path, jobs=sj)
    E.run_stage('dev', out_path=path, jobs=dj)
    E.run_stage('main', out_path=path, jobs=mj)
    mp.undo()
    return path


def test_stage_arrays_truth_and_default_scores(mini_raw):
    z = np.load(mini_raw)
    assert z['sanity_meta'].shape == (6, 7) and z['dev_meta'].shape == (8, 7) and z['main_meta'].shape == (10, 7)
    assert z['ref_cols'].shape[0] == len(E.REF_ROWS) and z['ref_res'].shape[0] == 10
    st = E.Stage(z, 'sanity')
    assert not st.bad.any() and list(st.K) == [1, 1, 1, 1, 2, 2] and list(st.plc) == [0, 0, 0, 1, 0, 1]
    ct, at = st.truth(4)
    c, a = E.sanity_target(2, 0, 503000)
    assert np.array_equal(ct, c) and np.array_equal(at, a)
    r = st.score('step', False, E.DEFAULT['step'][0])                                         # forward selection, on the grid, noiseless: exact (K = 1 trial 0 and K = 2 trial 4)
    assert r[0, 0] == 1 and r[0, 1] == 1 and r[4, 0] == 2 and r[4, 1] == 2
    assert st.score('step', False, E.DEFAULT['step'][0]) is st.score('step', False, E.DEFAULT['step'][0])      # cached
    r1 = st.score('step', True, E.DEFAULT['step'][1])
    assert r1[3, 0] == 1 and r1[3, 1] == 1                                                    # half a step off the grid: merged to one, in the right place
    assert st.score('step', False, E.DEFAULT['step'][0])[3, 0] > 1
    e, pl = st.exact(r)
    assert e[0] and pl[0]
    assert abs(E.exact_rate(st, 'step', False, E.DEFAULT['step'][0]) - float((r[:, 0] == st.K).mean())) < 1e-15
    assert np.isfinite(z['sanity_ctrl']).all() and z['sanity_ctrl'].shape == (6, 2, 4) and np.isnan(z['dev_ctrl']).all()
    # V2 columns are scores of slid clusters: finite, K-hat integer
    assert np.isfinite(z['sanity_v2'][:, :, 0]).all()
    assert np.all(z['sanity_v2'][0, 0] == z['sanity_v2'][0, 0])
    assert z['sanity_wall'] > 0 and z['chk_wall'] > 0


def test_report_is_byte_identical_and_has_every_section(mini_raw):
    a = E.report(mini_raw)
    assert a == E.report(mini_raw)
    for h in ('## What each baseline does', '## A0. Implementation checks', '## A1. Sanity cases', '### Exact order by case', '### Failure modes, on-grid features',
              '### Failure modes, off-grid features', '### Sensitivity of the merged count to the link', '### Optimality of the penalised network',
              '## A3. Hyper-parameters retuned', '## A4. The main comparison', '### Difference from the sliding ISS (BIC), clustered by target', '## Cost'):
        assert h in a, h
    z = np.load(mini_raw)
    mn = E.Stage(z, 'main')
    cols = [str(c) for c in z['ref_cols']]
    exp = [int((mn.score(m, False, E.DEFAULT[m][0])[:, 0] == z['ref_res'][:, cols.index(cc), 0]).sum()) for m, cc in (('step', 'step|bic'), ('pen', 'pen|bic'), ('prune', 'prune'))]
    assert 'K-hat identical to the stored raw in: forward selection %d / 10, penalised network %d / 10, random start + prune %d / 10 trials' % tuple(exp) in a
    assert '| sliding ISS, BIC | ' in a and 'Failed trials' not in a
    assert '240 trials' not in a and '(600 trials' not in a and '(6 trials, 0 failed)' in a


def test_report_counts_a_failed_trial_and_leaves_it_out(mini_raw, tmp_path):
    z = dict(np.load(mini_raw))
    for k in E._shapes():
        z['sanity_' + k] = np.concatenate([z['sanity_' + k], np.full((1,) + z['sanity_' + k].shape[1:], np.nan)])
    z['sanity_meta'][-1, :6] = [1, 0.01, 500009, 0, 0, 2]
    p = str(tmp_path / 'raw.npz')
    np.savez_compressed(p, **z)
    a = E.report(p)
    assert '(7 trials, 1 failed)' in a and a == E.report(p)
    b = E.report(mini_raw)
    sec = lambda t: t[t.index('### Exact order by case'):t.index('### Optimality of the penalised network')]      # noqa: E731
    assert sec(a) == sec(b)                                                                    # the failed trial is left out of every sanity table


def test_report_known_values_on_a_handmade_sanity_raw(tmp_path):
    """Four sanity trials (K = 1, noiseless, on the grid) with forward selection, penalised and prune outputs written by hand: the table cells are read off."""
    n = 5
    z = {}
    shapes = E._shapes()
    for k, v in shapes.items():
        z['sanity_' + k] = np.full((n,) + v, np.nan)
    z['sanity_meta'][:] = [1, 0.0, 500000, 0, 0, 0, 0.0]
    z['sanity_ct'][:, 0] = 0.3
    z['sanity_at'][:, 0] = 1.0
    z['sanity_step_idx'][:] = np.arange(E.MAXU)
    z['sanity_step_w'][:] = np.nan
    for k in range(1, E.MAXU + 1):
        z['sanity_step_w'][:, k - 1, :k] = 1.0 / k
    z['sanity_step_rss'][:] = np.r_[1.0, np.zeros(E.MAXU)]
    z['sanity_pen_w'][:] = 0.0
    z['sanity_pen_c'][:] = E.GRID
    z['sanity_pen_rss'][:] = 1.0
    z['sanity_pen_rss'][:, 1] = 1e-9
    z['sanity_pen_w'][:, 1, 13] = 0.6                                                          # penalty 1: two units at 0.3 (grid index 13 is 0.3)
    z['sanity_pen_w'][:, 1, 14] = 0.4
    z['sanity_pen_c'][:, 1, 13] = 0.3
    z['sanity_pen_c'][:, 1, 14] = 0.3
    z['sanity_pen_w'][4, 1, 13] = 1.0                                                          # the fifth trial: a single unit already
    z['sanity_pen_w'][4, 1, 14] = 0.0
    z['sanity_pr_w'][:] = 0.0
    z['sanity_pr_c'][:] = E.GRID
    z['sanity_pr_w'][:, :, 13] = 1.0
    z['sanity_ctrl'][:] = [[1, 1, 0.001, 0.001], [1, 1, 0.001, 0.001]]
    z['sanity_v2'][:] = [[1, 1, 0.0, 0.0]] * 3
    z['sanity_secs'][:] = 1.0
    p = str(tmp_path / 'raw.npz')
    np.savez_compressed(p, **z)
    a = E.report(p)
    # forward selection: rss 0 from one unit on, the floored BIC stops at 1 unit (weight 1.0 at grid index 0 = -1.0, not at 0.3!): count right, place wrong
    lines = {ln.split(' | ')[0]: ln for ln in a.splitlines() if ln.startswith('| ')}
    key = '| ' + MLABEL('step') + ', V0 (the comparison\'s count)'
    row = [ln for ln in a.splitlines() if ln.startswith(key)]
    assert len(row) >= 3                                                                       # the case table, the noise table and the failure-mode tables
    fm = [ln for ln in row if '| 5 |' in ln][0]                                                # failure-mode table: trials = 5
    cells = [c.strip() for c in fm.strip('|').split('|')]
    assert cells[2].startswith('100 %') and cells[3].startswith('0 %') and cells[6].startswith('100 %')      # exact order, not placed, right count wrong place
    pen = [ln for ln in a.splitlines() if ln.startswith('| ' + MLABEL('pen') + ', V0') and '| 5 |' in ln][0]
    cells = [c.strip() for c in pen.strip('|').split('|')]
    assert cells[2].startswith('20 %') and cells[4].startswith('80 %') and cells[7] == '+0.80'   # two nonzero units for one feature: over-estimates in four of five trials
    pen1 = [ln for ln in a.splitlines() if ln.startswith('| ' + MLABEL('pen') + ', V1') and '| 5 |' in ln][0]
    cells = [c.strip() for c in pen1.strip('|').split('|')]
    assert cells[2].startswith('100 %') and cells[3].startswith('100 %')                       # merged: one spike at 0.3, exact and placed
    assert 'Penalised network, over-estimating trials (V0 K-hat > K): 4 of 5 sanity trials; 4 of them are exact once the units are merged' in a
    assert 'J(network) <= J(truth)' in a


def MLABEL(m):
    return E.MLABEL[m]


def test_committed_report_is_regenerated_byte_identically():
    """The committed md equals the report regenerated from the committed raw (in a mutant slot both live in the output directory)."""
    if not os.path.exists(E.RAW) or not os.path.exists(E.MD):
        pytest.skip('no committed raw / md')
    assert open(E.MD).read() == E.report(E.RAW)
