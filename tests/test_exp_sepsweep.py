"""Tests of experiments/exp_sepsweep.py (the separation sweep on the trained network): the declared design (cells, drop rules, seeds), the targets (geometry, window, paired noise), the
outcome of a snapshot (formed threshold, merge, location, exactness), the stopping rules (the K-th entry, the window, the BIC), the path and one network on a real target, the statistics
(Clopper-Pearson, the transition and its bootstrap interval), the development rule, the registered predictions on planted raws, the checkpoints and the raw, and the report (planted raws,
the stored raws when present, byte-identical regeneration)."""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_sepsweep as X      # noqa: E402

SCRIPT = X.__file__
RES = os.path.join(HERE, '..', 'results', 'synthetic')


# ---- the declared design ---------------------------------------------------------------------------------------------------------------------------------
def test_constants_are_the_declared_ones():
    assert X.SIG == 0.12 and X.H_LAT == 0.1 and np.isclose(X.HX, 1.0 / 30.0, atol=0) and X.MARGIN == 0.5
    assert X.LAW == dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0) and X.A_MIN == 0.7
    assert X.KS == (2, 3, 4) and X.DELTAS == (1.5, 2.5, 3.5, 4.5, 5.5, 7.0) and X.RS == (1.0, 2.4, 10.0, 100.0)
    assert X.NOISES == (0.0, 0.01, 0.05, 0.2) and X.EPS == (1e-4, 1e-6, 1e-8)
    assert X.RS_K4 == (1.0, 10.0) and X.RS_MID == (1.0, 10.0) and (X.N_TEST, X.N_DEV) == (40, 6)
    assert (X.DEV_SEED0, X.TEST_SEED0) == (7_100_000, 8_200_000)
    assert (X.N_T, X.OMEGA_MAX, X.MAX_SECS, X.MAX_RHS, X.EXTRA_ENTRIES) == (601, 0.5, 300, 600_000, 2)
    assert (X.LOC_TOL, X.AMP_TOL, X.EXACT_POS, X.EXACT_AMP, X.LEVEL, X.CONF) == (0.5, 0.5, 0.05, 0.05, 0.9, 0.95)
    assert (X.BOOT, X.BOOT_SEED, X.BIC_FLOOR, X.D_THETA) == (2000, 0, 1e-10, 1e-3)
    assert X.STOPS == ('A05', 'A25', 'A50', 'BIC') and X.THRS == ('P', 'D') and X.RHOS == (0.0, 0.25, 0.5, 0.75)
    assert X.OMEGAS == {'A05': 0.05, 'A25': 0.25, 'A50': 0.5}
    assert X.N_CFG == 32 and X.CFGS[0] == ('A05', 'P', 0.0) and X.CFGS[-1] == ('BIC', 'D', 0.75) and X.CFGS[X.CFGS.index(('A25', 'P', 0.5))] == ('A25', 'P', 0.5)


def test_cells_and_drop_rules():
    C = X.CELLS
    assert len(C) == 192 and len(set(C)) == 192
    assert sum(c[0] == 2 for c in C) == 72 and sum(c[0] == 3 for c in C) == 72 and sum(c[0] == 4 for c in C) == 48
    assert all(X.RS[c[2]] in (1.0, 10.0) for c in C if c[0] == 4)
    assert all(X.RS[c[2]] in (1.0, 10.0) for c in C if X.NOISES[c[3]] in (0.01, 0.05))
    for K in X.KS:
        for iR, R in enumerate(X.RS):
            if K == 4 and R not in (1.0, 10.0):
                continue
            for iN in (0, 3):                                       # noise 0 and 0.2 in every (K, R)
                assert all((K, iD, iR, iN) in C for iD in range(6))
    assert len(X.targets('test', 40)) == 7680 and len(X.targets('dev', 6)) == 1152
    t = X.targets('test', 2)
    assert t[0] == (2, 0, 0, 0, 0) and t[1] == (2, 0, 0, 0, 1) and t[2] == (2, 1, 0, 0, 0)


def test_thresholds_match_theory3e_md_and_the_code():
    assert set(X.THRESH) == {(K, R) for K in X.KS for R in X.RS}
    md = open(os.path.join(RES, 'theory3e.md')).read()
    row = [ln for ln in md.splitlines() if ln.startswith('| Gaussian (n=1) | 4.56 |')][0]
    vals = [float(x) for x in row.split('|')[2:7]]
    assert vals == [X.THRESH[k] for k in ((2, 1.0), (2, 2.4), (2, 10.0), (3, 1.0), (3, 2.4))]
    assert X.THRESH_TAB == {(2, 1.0), (2, 2.4), (2, 10.0), (3, 1.0), (3, 2.4)}
    pytest.importorskip('flint')
    import theory3e_kernels as TK
    g = TK.Gaussian(1)
    for (K, R) in ((2, 100.0), (3, 10.0), (4, 1.0), (4, 100.0)):
        assert TK.entry_threshold(g, K, R, 0.4167, lo_=1.0, hi=80.0) == X.THRESH[(K, R)]


# ---- the targets -----------------------------------------------------------------------------------------------------------------------------------------
def test_make_target_geometry_window_and_amplitudes():
    for (K, iD, iR) in ((2, 5, 3), (3, 3, 1), (4, 0, 2)):
        tg = X.make_target(K, iD, iR, 0, 5, X.DEV_SEED0)
        cs, As = tg['cs'], tg['As']
        assert np.allclose(np.diff(cs), X.DELTAS[iD] * X.SIG, atol=1e-12) and 0.5 <= cs[0] < 0.6
        assert np.allclose(np.sort(As), 0.7 * X.RS[iR] ** (np.arange(K) / (K - 1)), atol=1e-12) and np.all(As > 0)
        win = tg['win']
        assert win * 0.1 >= cs[-1] + 0.5 - 1e-9 and (win - 1) * 0.1 < cs[-1] + 0.5
        _, Kobj, c0 = X.setup(win)
        assert Kobj.N == 3 * win + 1 and len(c0) == win + 1 and np.isclose(c0[1] - c0[0], 0.1) and np.isclose(Kobj.xs[1, 0] - Kobj.xs[0, 0], 1 / 30)
        assert np.array_equal(tg['y'], tg['f']) and np.allclose(tg['f'], Kobj.synth(cs[:, None], As))
    a = X.make_target(3, 3, 1, 0, 5, 1)
    b = X.make_target(3, 3, 1, 0, 5, 1)
    c = X.make_target(3, 3, 1, 0, 6, 1)
    d = X.make_target(3, 3, 1, 0, 5, 2)
    assert np.array_equal(a['cs'], b['cs']) and np.array_equal(a['As'], b['As']) and not np.array_equal(a['cs'], c['cs']) and not np.array_equal(a['cs'], d['cs'])


def test_noise_is_scaled_and_paired_across_levels():
    t1 = X.make_target(3, 2, 1, 1, 3, 77)
    t2 = X.make_target(3, 2, 1, 2, 3, 77)
    t3 = X.make_target(3, 2, 1, 3, 3, 77)
    assert np.array_equal(t1['cs'], t3['cs']) and np.array_equal(t1['As'], t3['As']) and np.array_equal(t1['f'], t3['f'])
    z = (t3['y'] - t3['f']) / (0.7 * 0.2)
    assert np.allclose((t1['y'] - t1['f']) / (0.7 * 0.01), z, atol=1e-10) and np.allclose((t2['y'] - t2['f']) / (0.7 * 0.05), z, atol=1e-10)
    assert 0.5 < z.std() < 1.6 and abs(z.mean()) < 0.5
    tr = X.make_target(2, 1, 3, 0, 0, 77)                              # R = 100: the noise is relative to the SMALLEST amplitude
    tn = X.make_target(2, 1, 3, 3, 0, 77)
    assert np.allclose((tn['y'] - tn['f']) / 0.14, (X.make_target(2, 1, 3, 3, 0, 77)['y'] - tr['f']) / 0.14) and 0.5 < ((tn['y'] - tn['f']) / 0.14).std() < 1.6


# ---- one snapshot ---------------------------------------------------------------------------------------------------------------------------------------
def test_formed_threshold():
    A = np.array([0.7, 7.0])
    assert np.isclose(X.formed_threshold('P', None, A), min(0.07, 0.07)) and np.isclose(X.formed_threshold('P', None, np.array([0.7, 70.0])), 0.07)
    assert np.isclose(X.formed_threshold('P', None, np.array([1.0, 1.0])), 0.01)
    assert np.isclose(X.formed_threshold('P', None, np.array([1.0, 100.0])), 0.1)               # 0.1 min|A*| = 0.1 < 0.01 max = 1
    assert np.isclose(X.formed_threshold('D', np.array([0.5, -2.0, 1.0]), A), 2e-3)


def test_merge_units_clusters_weights_and_radius():
    w = np.array([1.0, 3.0, 0.0, 2.0, 1e-9, 5.0])
    c = np.array([0.30, 0.34, 0.50, 0.90, 0.91, 1.60])
    cen, amp, idx = X.merge_units(w, c, 1e-6, 0.05)
    assert list(idx) == [0, 1, 3, 5] and cen.size == 3 and np.allclose(amp, [4.0, 2.0, 5.0]) and np.isclose(cen[0], (1 * 0.30 + 3 * 0.34) / 4)
    cen0, amp0, _ = X.merge_units(w, c, 1e-6, 0.0)
    assert cen0.size == 4 and np.allclose(cen0, [0.30, 0.34, 0.90, 1.60])
    cen1, amp1, _ = X.merge_units(np.ones(3), np.array([0.0, 0.25, 1.0]), 0.0, 0.25)       # a gap of exactly the radius joins (<= radius)
    assert cen1.size == 2 and np.allclose(cen1, [0.125, 1.0])
    assert X.merge_units(np.ones(3), np.array([0.0, 0.25, 1.0]), 0.0, 0.2499)[0].size == 3
    cen2, _, _ = X.merge_units(w, c, 1e-6, 0.6)                      # single linkage: 0.30, 0.34, 0.90 chain (gaps 0.04, 0.56), 1.60 is 0.70 away
    assert cen2.size == 2 and X.n_merged(w, c, 1e-6, 0.6) == 2 and X.merge_units(w, c, 1e-6, 1.0)[0].size == 1
    e, a, i = X.merge_units(np.zeros(3), np.arange(3.0), 0.0, 1.0)
    assert e.size == 0 and a.size == 0 and i.size == 0
    cen3, amp3, _ = X.merge_units(np.array([2.0, -1.0]), np.array([0.0, 0.01]), 0.0, 0.1)     # signed sum, |w|-weighted centre
    assert np.isclose(amp3[0], 1.0) and np.isclose(cen3[0], 0.01 / 3)
    o_w, o_c = np.array([1.0, 2.0]), np.array([0.5, 0.1])                                        # unsorted input is sorted by centre
    cen4, amp4, _ = X.merge_units(o_w, o_c, 0.0, 0.0)
    assert np.allclose(cen4, [0.1, 0.5]) and np.allclose(amp4, [2.0, 1.0])


def test_located_tolerances():
    cs, As = np.array([0.5, 1.0]), np.array([1.0, 2.0])
    s = X.SIG
    assert X.located(np.array([0.5, 1.0]), np.array([1.0, 2.0]), cs, As, 0.5, 0.5)
    assert X.located(np.array([0.5 + 0.49 * s, 1.0 - 0.49 * s]), np.array([1.49, 1.01]), cs, As, 0.5, 0.5)
    assert not X.located(np.array([0.5 + 0.51 * s, 1.0]), np.array([1.0, 2.0]), cs, As, 0.5, 0.5)
    assert not X.located(np.array([0.5, 1.0]), np.array([1.51, 2.0]), cs, As, 0.5, 0.5)
    assert not X.located(np.array([0.5]), np.array([1.0]), cs, As, 0.5, 0.5) and not X.located(np.array([0.5, 0.7, 1.0]), np.array([1.0, 1.0, 2.0]), cs, As, 0.5, 0.5)
    assert X.located(np.array([0.5, 1.0]), np.array([1.0, 2.0]), np.array([1.0, 0.5]), np.array([2.0, 1.0]), 0.5, 0.5)      # features are sorted by position first
    assert not X.located(np.array([0.5, 1.0]), np.array([2.0, 1.0]), np.array([1.0, 0.5]), np.array([2.0, 1.0]), 0.5, 0.5)
    assert not X.located(np.array([0.5, 1.0]), np.array([-1.0, 2.0]), cs, As, 0.5, 0.5)           # wrong sign
    assert X.located(np.array([0.5 + 0.04 * s]), np.array([1.04]), np.array([0.5]), np.array([1.0]), X.EXACT_POS, X.EXACT_AMP)
    assert not X.located(np.array([0.5 + 0.06 * s]), np.array([1.0]), np.array([0.5]), np.array([1.0]), X.EXACT_POS, X.EXACT_AMP)


def test_snapshot_outcome_fields():
    cs, As = np.array([0.5, 1.0]), np.array([1.0, 2.0])
    w = np.array([0.6, 0.4, 2.0, 1e-7, 0.0])
    c = np.array([0.50, 0.505, 1.0, 0.0, 0.2])
    ok, idx, nf, nm, ex = X.snapshot_outcome(w, c, cs, As, 'P', 0.5)
    assert ok and idx == (0, 1, 2) and nf == 3 and nm == 2 and ex
    ok0, idx0, nf0, nm0, ex0 = X.snapshot_outcome(w, c, cs, As, 'P', 0.0)
    assert not ok0 and idx0 == (0, 1, 2) and nm0 == 3 and not ex0
    ok2, idx2, nf2, nm2, ex2 = X.snapshot_outcome(w, c, cs, As, 'D', 0.5)          # D: 1e-3 * 2.0 = 2e-3: the 1e-7 unit stays out
    assert ok2 and idx2 == (0, 1, 2)
    w3 = np.array([0.6, 0.4, 2.0, 5e-3, 0.0])
    assert X.snapshot_outcome(w3, c, cs, As, 'D', 0.0)[2] == 4 and X.snapshot_outcome(w3, c, cs, As, 'P', 0.0)[2] == 3     # P: 0.01 * 2 = 0.02 > 5e-3
    assert isinstance(idx, tuple) and all(isinstance(j, int) for j in idx)


def test_bic_value():
    assert np.isclose(X.bic_value(2.0, 100, 3, 1e-3), 100 * math.log(0.02 + 1e-3) + 6 * math.log(100))
    assert X.bic_value(2.0, 100, 4, 0.0) > X.bic_value(2.0, 100, 3, 0.0)


# ---- the stopping rules ----------------------------------------------------------------------------------------------------------------------------------
def test_index_of_kth_entry():
    st = [(0.0, None, None, 1.0), (1.0, None, None, 0.5), (2.0, None, None, 0.2), (2.5, None, None, 0.1)]
    assert X.index_of_kth_entry(st, [1.0, 2.0, 2.5], 2) == (2, 2.0, True)
    assert X.index_of_kth_entry(st, [1.0, 2.0, 2.5], 1) == (1, 1.0, True)
    assert X.index_of_kth_entry(st, [1.0, 2.0, 2.5], 3) == (3, 2.5, True)
    assert X.index_of_kth_entry(st, [1.0, 1.0, 2.5], 2) == (1, 1.0, True)        # ties enter together: both at the state of tau = 1
    i, t, found = X.index_of_kth_entry(st, [1.0], 2)
    assert (i, t, found) == (3, 2.5, False)


def _planted_res(win, cs, As, extra=True):
    """A flow record with 5 snapshots: nothing, the K true units, the K true units plus one tiny noise unit, ..."""
    _, Kobj, c0 = X.setup(win)
    d = len(c0)
    T = 5
    W = np.zeros((T, d))
    C = np.tile(c0[None, :, None], (T, 1, 1)).astype(float)
    idx = [int(np.argmin(np.abs(c0 - c))) for c in cs]
    for t in range(1, T):
        for j, a in zip(idx, As):
            W[t, j] = a
            C[t, j, 0] = cs[idx.index(j)]
    if extra:
        k = (idx[0] + 5) % d
        W[3, k] = 0.3
        W[4, k] = 0.3
    return dict(taus=np.linspace(0.0, 4.0, T), W=W, C=C, status='ok')


def test_evaluate_network_rule_a_picks_the_window_snapshot_and_bic_the_best():
    tg = X.make_target(2, 5, 0, 0, 0, 5)
    cs, As, y, win = tg['cs'], tg['As'], tg['y'], tg['win']
    res = _planted_res(win, cs, As)
    psets = [(0, 1)] * X.N_CFG
    out = X.evaluate_network(win, cs, As, y, res, psets, 2.2)               # tau_K = 2.2: A05 -> tau 2.31 (snapshot 2), A25 -> 2.75 (snapshot 3: the extra unit), A50 -> 3.3 (snapshot 3)
    ic = X.CFGS.index(('A05', 'P', 0.5))
    assert out[ic][0] == 1 and out[ic][2] == 2                                # snapshot 2: exactly the K true units
    for stop in ('A25', 'A50'):
        ic = X.CFGS.index((stop, 'P', 0.5))
        assert out[ic][0] == 0 and out[ic][2] == 3 and out[ic][3] == 3         # snapshot 3: a third (noise) unit, far from the features
    ic = X.CFGS.index(('BIC', 'P', 0.5))
    assert out[ic][0] == 1 and out[ic][2] == 2 and out[ic][4] == 1             # BIC: the K-unit snapshot (the extra unit only adds a penalty, the fit is exact)
    c_idx = X.CFGS.index(('A05', 'P', 0.5))
    psets2 = [tuple(int(np.argmin(np.abs(X.setup(win)[2] - c))) for c in cs)] * X.N_CFG
    out2 = X.evaluate_network(win, cs, As, y, res, psets2, 2.2)
    assert out2[c_idx][1] == 1 and out[c_idx][1] == 0                        # agreement: the formed index set equals the path's


def test_evaluate_path_rules():
    tg = X.make_target(2, 5, 0, 0, 0, 5)
    cs, As, y, win = tg['cs'], tg['As'], tg['y'], tg['win']
    _, Kobj, c0 = X.setup(win)
    d = len(c0)
    j = [int(np.argmin(np.abs(c0 - c))) for c in cs]

    def state(tau, units):
        w, c = np.zeros(d), c0.copy()
        for (k, a, pos) in units:
            w[k], c[k] = a, pos
        r = y - Kobj.synth(c[:, None], w)
        return (tau, c, w, float(np.sqrt(np.mean(r ** 2))))
    states = [state(0.0, []), state(1.0, [(j[0], As[0], cs[0])]), state(2.0, [(j[0], As[0], cs[0]), (j[1], As[1], cs[1])]),
              state(3.0, [(j[0], As[0], cs[0]), (j[1], As[1], cs[1]), ((j[0] + 6) % d, 0.2, c0[(j[0] + 6) % d])])]
    out, sets, tauK, found, iK = X.evaluate_path(win, 2, cs, As, y, states, [1.0, 2.0, 3.0])
    assert iK == 2 and tauK == 2.0 and found
    ia = X.CFGS.index(('A25', 'P', 0.5))
    ib = X.CFGS.index(('BIC', 'P', 0.5))
    assert out[ia][0] == 1 and out[ia][2] == 2 and sets[ia] == tuple(sorted(j)) and out[ia][3] == 1
    assert out[ib][0] == 1 and out[ib][2] == 2                                # BIC: the second state (the third adds a unit and no fit)


# ---- the path and the network on a real target -----------------------------------------------------------------------------------------------------------
@pytest.fixture(scope='module')
def real():
    tg = X.make_target(2, 5, 0, 0, 0, 5)                                      # K = 2, Delta = 7, no noise
    states, ent, info = X.run_path(tg['win'], 2, tg['y'])
    return tg, states, ent, info


def test_path_on_a_separated_pair(real):
    tg, states, ent, info = real
    assert len(ent) == 2 and not info['oob'] and info['failed'] == 0 and info['exits'] == 0
    t, c, w, r = states[-1]
    assert r < 1e-8 and np.sum(w != 0) == 2
    assert np.allclose(np.sort(c[w != 0]), tg['cs'], atol=1e-6) and np.allclose(np.sort(w[w != 0]), np.sort(tg['As']), atol=1e-6)
    assert states[0][2].max() == 0.0 and len(states) == 3 and states[-1][0] == ent[-1]
    out, sets, tauK, found, iK = X.evaluate_path(tg['win'], 2, tg['cs'], tg['As'], tg['y'], states, ent)
    assert found and iK == 2 and tauK == ent[1] and out[:, 0].min() == 1 and out[:, 3].min() == 1
    assert X.index_of_kth_entry(states, ent, 2) == (2, ent[1], True)


def test_run_target_on_a_separated_pair_forms_the_pair_at_every_eps():
    row = (2, 5, 0, 0, 0)
    o = X.run_target(row, 5)
    assert o['path'].shape == (X.N_CFG, 4) and o['net'].shape == (3, X.N_CFG, 5) and o['diag'].shape == (17,)
    ia = X.CFGS.index(('A25', 'P', 0.5))
    assert o['path'][ia][0] == 1
    assert [int(o['net'][ie][ia][0]) for ie in range(3)] == [1, 1, 1]          # success at every eps (merge radius 0.5)
    assert [int(o['net'][ie][ia][1]) for ie in range(3)] == [0, 1, 1]          # but at eps = 1e-4 the second feature is carried by two co-located units: not the path's set
    assert o['net'][0][ia][2] == 3 and o['net'][0][ia][3] == 2 and o['net'][1][ia][2] == 2
    assert list(o['diag'][8::3]) == [0.0, 0.0, 0.0] and o['diag'][1] == 1.0 and o['diag'][7] == 2.0 and o['diag'][0] > 0
    assert o['net'][2, :, 0][[X.CFGS.index(c) for c in X.CFGS if c[2] > 0 and c[0] != 'BIC']].min() == 1      # eps = 1e-8: the merging configurations with a window succeed


def test_status_code():
    assert X.status_code('ok') == 0 and X.status_code('timeout') == 1 and X.status_code('failed: x') == 2


# ---- statistics -----------------------------------------------------------------------------------------------------------------------------------------
def test_clopper_pearson():
    assert X.cp_interval(0, 10) == (0.0, pytest.approx(0.3084971, abs=1e-6))
    assert X.cp_interval(10, 10)[1] == 1.0 and X.cp_interval(10, 10)[0] == pytest.approx(0.6915029, abs=1e-6)
    lo, hi = X.cp_interval(5, 10)
    assert lo == pytest.approx(0.1870860, abs=1e-6) and hi == pytest.approx(0.8129140, abs=1e-6)
    assert X.cp_interval(40, 40)[0] == pytest.approx(0.9119, abs=1e-3) and X.cp_interval(39, 40)[0] < 0.9


def test_crossing_interpolation_flags_and_last_run():
    g = X.DELTAS
    assert X.crossing([0, 0, 1, 1, 1, 1]) == (pytest.approx(2.5 + 0.9), '')
    assert X.crossing([0, 0.5, 1, 1, 1, 1])[0] == pytest.approx(2.5 + (0.9 - 0.5) / 0.5 * 1.0)
    assert X.crossing([1, 1, 1, 1, 1, 1]) == (1.5, '<') and X.crossing([0, 0, 0, 0, 0, 0]) == (float('inf'), '>') and X.crossing([1, 1, 1, 1, 1, 0.9]) == (float('inf'), '>')
    assert X.crossing([1, 0.5, 1, 1, 1, 1])[0] == pytest.approx(2.5 + 0.4 / 0.5)     # the smallest grid point from which EVERY larger one exceeds the level
    assert X.crossing([0.95, 0.95, 0.5, 0.95, 0.95, 0.95])[0] == pytest.approx(3.5 + (0.9 - 0.5) / 0.45 * 1.0)
    assert X.crossing([0.9, 1, 1, 1, 1, 1])[0] == 1.5 + 0.0 and X.crossing([0.9, 1, 1, 1, 1, 1])[1] == ''      # not above at 1.5 (0.9 is not > 0.9): interpolated from the left end
    assert g == (1.5, 2.5, 3.5, 4.5, 5.5, 7.0)


def test_crossing_many_equals_crossing_row_by_row():
    rng = np.random.default_rng(3)
    R = rng.random((600, 6))
    R[:100] = np.sort(R[:100], 1)
    R[100:150] = 1.0
    R[150:200] = 0.0
    R[200:300] = np.round(R[200:300], 1)
    R[300:350, 3:] = 0.95
    a = X.crossing_many(R)
    b = np.array([X.crossing(r)[0] for r in R])
    assert np.array_equal(a, b)
    assert np.isinf(a[150:200]).all() and (a[100:150] == 1.5).all()


def test_sensitivity_rows_cover_every_configuration():
    raw = planted('test', 10, net_from=(5.5, 4.5, 3.5), path_from=3.5)
    rows = X.sensitivity_rows(raw)
    assert len(rows) == 32 and [r[0] for r in rows] == list(range(32)) and rows[0][1:4] == ['A05', 'P', '0']
    assert all(r[4] == 'held' and r[5] == 'held' and r[6] == 'failed' and r[7] == 'held' and r[8] == 10 for r in rows)
    assert rows[0][9] == '3.40' and rows[0][10] == '3.40'


def test_transition_and_its_bootstrap_interval():
    n = [40] * 6
    t = X.transition([0, 0, 40, 40, 40, 40], n)
    assert t[0] == pytest.approx(3.4) and t[1] == '' and t[2] <= t[0] <= t[4] and t[3] == '' and t[5] == ''
    t2 = X.transition([40] * 6, n)
    assert t2[:2] == (1.5, '<') and t2[2] == 1.5 and t2[4] == 1.5
    t3 = X.transition([0] * 6, n)
    assert t3[:2] == (float('inf'), '>') and t3[5] == '>' and t3[3] == '>'
    t4 = X.transition([2, 10, 30, 36, 38, 39], n)
    assert t4 == X.transition([2, 10, 30, 36, 38, 39], n)                       # seeded: reproducible
    assert t4[2] < t4[4] and t4[1] == '' and 3.5 <= t4[0] <= 7.0
    assert t4[:2] == (pytest.approx(4.5), '') and t4[2] == pytest.approx(4.0) and t4[4] == pytest.approx(6.1) and t4[3] == '' and t4[5] == ''        # pins the percentiles and the seed
    t5 = X.transition([0, 5, 10, 10, 10, 10], [10] * 6)
    assert t5[0] == pytest.approx(3.3) and t5[2] == pytest.approx(3.0) and t5[4] == pytest.approx(3.375)
    assert X.fmt_t(float('inf'), '>') == '>7' and X.fmt_t(1.5, '<') == '<=1.5' and X.fmt_t(3.456, '') == '3.46'
    assert X.fmt_tr((3.4, '', 2.5, '', 4.5, '')) == '3.40 [2.50, 4.50]'
    assert X.flag_of(float('inf')) == '>' and X.flag_of(1.5) == '<' and X.flag_of(2.0) == ''
    assert X.diff(float('inf'), float('inf')) == 0.0 and X.diff(5.0, 3.0) == 2.0 and X.diff(3.0, float('inf')) == -float('inf')
    assert X.fmt_diff(float('inf')) == 'censored (>7)' and X.fmt_diff(0.25) == '+0.25'


# ---- planted raws: aggregate, development rule, predictions, report ------------------------------------------------------------------------------------------
IC = X.CFGS.index(('A25', 'P', 0.5))


def planted(split='test', n=10, net_from=(3.5, 3.5, 3.5), path_from=2.5, noisy_shift=0.0, only_cells=None, seed=0):
    """A raw whose outcomes are deterministic: the network (eps index e) succeeds iff Delta >= net_from[e] + noisy_shift * noise_index (all configurations), the path iff
    Delta >= path_from; agreement = network success; exact = network success at noise 0."""
    rows = [r for r in X.targets(split, n) if only_cells is None or r[:4] in only_cells]
    T = len(rows)
    rows = np.array(rows, np.int16)
    path = np.zeros((T, X.N_CFG, 4), np.int16)
    net = np.zeros((T, 3, X.N_CFG, 5), np.int16)
    diag = np.zeros((T, 17))
    diag[:, 1] = 1.0
    for t, (K, iD, iR, iN, idx) in enumerate(rows):
        D = X.DELTAS[iD]
        path[t, :, 0] = D >= path_from
        path[t, :, 3] = (D >= path_from) and iN == 0
        for e in range(3):
            ok = D >= net_from[e] + noisy_shift * iN
            net[t, e, :, 0] = ok
            net[t, e, :, 1] = ok
            net[t, e, :, 4] = ok and iN == 0
    return {split + '_rows': rows, split + '_path': path, split + '_net': net, split + '_diag': diag}


def test_aggregate_counts_and_discordance():
    raw = planted('test', 10, net_from=(5.5, 4.5, 3.5), path_from=4.5)
    agg = X.aggregate(raw, 'test', ('A25', 'P', 0.5))
    assert len(agg) == 192 and all(d['n'] == 10 for d in agg.values())
    d = agg[(2, 0, 3, 5)]                                                         # K = 2, R = 1.. (iR = 0 -> R = 1; index tuple is (K, iR, iN, iD)): see below
    d = agg[(2, 0, 0, 3)]                                                         # Delta = 4.5: path yes; network eps 1e-4 no (5.5), 1e-6 yes, 1e-8 yes
    assert d['path'] == 10 and d['net'] == [0, 10, 10] and d['only_path'] == [10, 0, 0] and d['only_net'] == [0, 0, 0] and d['agree'] == [0, 10, 10]
    assert d['exact'] == [0, 10, 10] and d['path_exact'] == 10 and d['bad'] == [0, 0, 0] and d['nofound'] == 0
    d = agg[(2, 0, 0, 2)]                                                         # Delta = 3.5: path no; network eps 1e-8 yes
    assert d['path'] == 0 and d['net'] == [0, 0, 10] and d['only_net'] == [0, 0, 10] and d['only_path'] == [0, 0, 0]
    d = agg[(2, 0, 3, 3)]                                                         # noise 0.2: exact never
    assert d['exact'] == [0, 0, 0] and d['path_exact'] == 0
    ks, ns = X.curve(agg, 2, 0, 0, 'net', 2)
    assert ks == [0, 0, 10, 10, 10, 10] and ns == [10] * 6
    ks, ns = X.curve(agg, 2, 0, 0, 'path')
    assert ks == [0, 0, 0, 10, 10, 10]
    assert X.group_keys(agg)[0] == (2, 0, 0) and len(X.group_keys(agg)) == 32


def test_aggregate_reads_the_chosen_configuration_and_diagnostics():
    raw = planted('test', 2, only_cells=[(2, iD, 0, 0) for iD in range(6)])
    jc = X.CFGS.index(('BIC', 'D', 0.25))
    raw['test_net'][:, :, jc, 0] = 1 - raw['test_net'][:, :, jc, 0]               # flip one configuration only
    raw['test_diag'][0, 4] = 1.0
    raw['test_diag'][1, 5] = 2.0
    raw['test_diag'][2, 1] = 0.0
    raw['test_diag'][3, 6] = 1.0
    raw['test_diag'][4, 8] = 1.0
    raw['test_diag'][5, 14] = 3.0
    a1 = X.aggregate(raw, 'test', ('A25', 'P', 0.5))
    a2 = X.aggregate(raw, 'test', ('BIC', 'D', 0.25))
    assert a1[(2, 0, 0, 0)]['net'][2] == 0 and a2[(2, 0, 0, 0)]['net'][2] == 2
    assert a1[(2, 0, 0, 0)]['oob'] == 2 and a1[(2, 0, 0, 1)]['nofound'] == 1 and a1[(2, 0, 0, 1)]['exits'] == 1 and a1[(2, 0, 0, 2)]['bad'][0] == 1 and a1[(2, 0, 0, 2)]['bad'][2] == 1
    with pytest.raises(ValueError):
        X.aggregate(raw, 'test', ('A25', 'P', 0.6))


def test_dev_objective_and_selection_rule():
    raw = planted('dev', 6, net_from=(3.5, 3.5, 3.5), path_from=2.5)
    T = raw['dev_rows'].shape[0]
    assert T == 1152
    obj0 = X.dev_objective(raw)
    assert len(obj0) == 32 and all(o == (pytest.approx(0.5 * (4 / 6 + 5 / 6)), pytest.approx(4 / 6), pytest.approx(5 / 6)) for o in obj0)     # ALL separations count
    assert X.select_config(raw) == 0                                                  # a tie of all: the first in order
    raw['dev_net'][:, 2, 5, 0] = 1
    raw['dev_path'][:, 5, 0] = 1
    obj = X.dev_objective(raw)
    assert obj[5] == (1.0, 1.0, 1.0) and X.select_config(raw) == 5
    raw['dev_net'][:, 2, 9, 0] = 1
    raw['dev_path'][:, 9, 0] = 1                                                      # a tie at the maximum: the first in order wins
    assert X.select_config(raw) == 5
    raw['dev_net'][:, 2, 2, 0] = 1
    raw['dev_path'][:, 2, 0] = 1
    assert X.select_config(raw) == 2
    raw['dev_net'][:, 2, 1, 0] = 1
    raw['dev_path'][:, 1, 0] = 1
    raw['dev_path'][:5, 1, 0] = 0                                                     # configuration 1 misses 5 targets (objective 0.9978): within 0.005 of the maximum, earlier in the order
    assert X.dev_objective(raw)[1][0] == pytest.approx(1.0 - 0.5 * 5 / T) and X.select_config(raw) == 1
    raw['dev_path'][:20, 1, 0] = 0                                                    # misses 20 targets (0.9913): outside the tolerance
    assert X.select_config(raw) == 2
    # the objective reads only eps = 1e-8, and the network's success only through field 0
    raw2 = planted('dev', 6, net_from=(3.5, 3.5, 3.5), path_from=2.5)
    raw2['dev_net'][:, 0, :, 0] = 0
    raw2['dev_net'][:, 1, :, 0] = 0
    raw2['dev_net'][:, 2, :, 1:] = 0
    assert X.dev_objective(raw2) == obj0
    low = raw2['dev_rows'][:, 1] == 0
    raw2['dev_net'][low, 2, 3, 0] = 1                                                 # success at Delta = 1.5 raises the objective of that configuration only
    o2 = X.dev_objective(raw2)
    assert o2[3][1] == pytest.approx(5 / 6) and o2[4] == obj0[4]


def test_predictions_on_planted_raws():
    # network (eps 1e-8 / 1e-6 / 1e-4) transitions 3.4 (from 3.5), path 2.4 (from 2.5): P1 holds (3.4 <= every threshold), P2 holds only if the difference is <= 1: it is exactly 1.0 (floating point)
    raw = planted('test', 10, net_from=(5.5, 4.5, 3.5), path_from=3.5)
    P = X.predictions(X.aggregate(raw, 'test', ('A25', 'P', 0.5)))
    r = P['rows'][0]
    assert (r['K'], r['R'], r['thr']) == (2, 1.0, 4.56) and len(P['rows']) == 10
    assert r['net8'][0] == pytest.approx(3.4) and r['net4'][0] == pytest.approx(5.4) and r['net6'][0] == pytest.approx(4.4) and r['path'][0] == pytest.approx(3.4)
    assert P['P1'] and r['p1'] and P['P2'] and r['p2_diff'] == pytest.approx(0.0, abs=1e-12)
    assert P['P4'] and P['P4_counts'] == (10, 10, 10) and r['p4_diff'] == pytest.approx(2.0)
    assert not P['P3']                                                              # noise does not move the transition: 3.4 < thresholds
    assert P['P3_n'] == 7 and sum(1 for x in P['rows'] if x['p3_testable']) == 7
    assert [x['p3'] for x in P['rows'] if not x['p3_testable']] == [None, None, None] and not any(x['thr'] <= 7 for x in P['rows'] if not x['p3_testable'])
    # shifting the network to 7 for noise 0.2 (noisy_shift 1.25 per noise index: 3.5 + 3 * 1.25 = 7.25 > 7: never) moves it above the thresholds
    raw = planted('test', 10, net_from=(3.5, 3.5, 3.5), path_from=1.5, noisy_shift=1.25)
    P = X.predictions(X.aggregate(raw, 'test', ('A25', 'P', 0.5)))
    assert P['P3'] and P['P1'] and not P['P2'] and all(x['net8_n02'][1] == '>' for x in P['rows'])
    assert not P['P4'] and P['P4_counts'] == (10, 0, 10) and all(x['p2_diff'] == pytest.approx(3.4 - 1.5) for x in P['rows'])
    # P1 fails when the network transitions above the threshold: network from 7 (threshold 4.56 .. 7.10 < 7 except the R = 100 rows)
    raw = planted('test', 10, net_from=(7.0, 7.0, 7.0), path_from=1.5)
    P = X.predictions(X.aggregate(raw, 'test', ('A25', 'P', 0.5)))
    assert not P['P1'] and [x['p1'] for x in P['rows']] == [x['net8'][0] <= x['thr'] for x in P['rows']]
    assert all(x['net8'][0] == pytest.approx(5.5 + 0.9 * 1.5) for x in P['rows']) and sum(x['p1'] for x in P['rows']) == 4      # 6.85: the thresholds of (2, 100), (3, 10), (3, 100), (4, 10) hold


def test_predictions_p4_rule_boundaries():
    groups = [(K, iR) for K in X.KS for iR in range(4) if not (K == 4 and X.RS[iR] not in (1.0, 10.0))]
    assert len(groups) == 10

    def p4(strict_groups, lower_groups=()):
        raw = planted('test', 10, net_from=(3.5, 3.5, 3.5), path_from=1.5)
        for t, (K, iD, iR, iN, idx) in enumerate(raw['test_rows'].tolist()):
            if iN == 0 and (K, iR) in strict_groups:
                raw['test_net'][t, 0, :, 0] = X.DELTAS[iD] >= 5.5                    # eps 1e-4 from 5.5: strictly higher than 1e-8 (3.5)
            if iN == 0 and (K, iR) in lower_groups:
                raw['test_net'][t, 0, :, 0] = X.DELTAS[iD] >= 2.5                    # lower
        return X.predictions(X.aggregate(raw, 'test', ('A25', 'P', 0.5)))
    P = p4(set(groups[:7]))
    assert P['P4_counts'] == (10, 7, 10) and P['P4']                                 # 7 of 10 is two thirds (rounded up)
    P = p4(set(groups[:6]))
    assert P['P4_counts'] == (10, 6, 10) and not P['P4']
    P = p4(set(groups[:9]), lower_groups={groups[9]})
    assert P['P4_counts'] == (9, 9, 10) and not P['P4']                              # one lower: fails whatever the rest
    P = p4(set(groups))
    assert P['P4_counts'] == (10, 10, 10) and P['P4']


def test_predictions_p4_small_differences_are_not_strict():
    raw = planted('test', 10, net_from=(3.5, 3.5, 3.5), path_from=1.5)
    for t, (K, iD, iR, iN, idx) in enumerate(raw['test_rows'].tolist()):
        if iN == 0 and iD == 1:
            raw['test_net'][t, 2, :, 0] = idx < 8                                    # eps 1e-8: 8/10 at Delta = 2.5 -> transition 3.0
            raw['test_net'][t, 0, :, 0] = idx < 7                                    # eps 1e-4: 7/10 at Delta = 2.5 and perfect from 3.5 -> 3.1667 (higher by 0.17)
        elif iN == 0 and iD >= 2:
            raw['test_net'][t, 0, :, 0] = 1
    P = X.predictions(X.aggregate(raw, 'test', ('A25', 'P', 0.5)))
    r = P['rows'][0]
    assert r['net8'][0] == pytest.approx(3.0) and r['net4'][0] == pytest.approx(2.5 + 0.2 / 0.3) and 0.1 < r['p4_diff'] < 0.25
    assert P['P4_counts'] == (10, 0, 10) and not P['P4']


def test_run_path_with_noise_runs_to_k_plus_extra_entries_and_run_target_sets_the_flow_horizon(monkeypatch):
    tg = X.make_target(2, 5, 0, 3, 0, 5)                                            # K = 2, Delta = 7, noise 0.2
    states, ent, info = X.run_path(tg['win'], 2, tg['y'])
    assert len(ent) == 2 + X.EXTRA_ENTRIES and len(states) == 2 + X.EXTRA_ENTRIES + 1 and states[-1][0] == ent[-1] and states[-1][3] > 1e-3
    seen = []
    orig = X.run_network

    def spy(win, y, eps, tau_end):
        seen.append((eps, tau_end))
        return orig(win, y, eps, tau_end)
    monkeypatch.setattr(X, 'run_network', spy)
    o = X.run_target((2, 5, 0, 0, 0), 5)
    tauK = o['diag'][0]
    assert [e for e, t in seen] == [1e-4, 1e-6, 1e-8] and all(np.isclose(t, 1.5 * tauK + 1.0) for e, t in seen)


def test_make_target_assigns_the_amplitudes_at_random():
    pos = {int(np.argmax(X.make_target(3, 3, 2, 0, i, 11)['As'])) for i in range(30)}
    assert pos == {0, 1, 2}


def test_report_planted_test_raw_and_regeneration(tmp_path):
    old = X.FROZEN
    X.FROZEN = ('A25', 'P', 0.5)
    try:
        raw = planted('test', 10, net_from=(5.5, 4.5, 3.5), path_from=3.5)
        txt = X.report({'test': raw})
        assert txt == X.report({'test': raw}) and txt.startswith('# Separation sweep on the trained network')
        assert '192 cells' in txt and '**A25, threshold P, merge radius 0.5**' in txt
        assert '| 2 | 1 | 4.56* | 3.40 [' in txt and '| 4 | 10 | 7.10 |' in txt and '| 3 | 100 | 8.42 |' in txt
        assert '**P1**' in txt and 'HELD' in txt and '**P3**' in txt and 'FAILED' in txt
        assert '## 4. Transitions at every noise level' in txt and '## 7. Where the network and the path part' in txt
        assert '0/10 [0.00, 0.31]' in txt and '10/10 [0.69, 1.00]' in txt
        assert txt.count('\n| 2 | 1 | 0 |') >= 3
        # incomplete cells are not analysed
        part = planted('test', 10, only_cells=[c for c in X.CELLS if not (c[0] == 2 and c[2] == 0 and c[3] == 3 and c[1] == 4)])
        t2 = X.report({'test': part})
        assert '| 2 | 1 | 0.2 |' not in t2 and 'Cells missing or incomplete (no claim is made for them): (K=2, R=1, noise=0.2).' in t2
        assert 'Cells missing' not in txt and '## 8. Sensitivity (registered)' in txt and txt.count('\n| ') > 350
    finally:
        X.FROZEN = old
    # the CLI regenerates the md from the raws byte for byte
    out = tmp_path / 'o'
    out.mkdir()
    np.savez_compressed(out / 'sepsweep_raw.npz', **raw)
    env = dict(os.environ, SS_OUT=str(out))
    if old is None:                                                                      # before the freeze the report cannot analyse a test raw
        pytest.skip('the shipped script has no frozen configuration yet')
    for _ in range(2):
        subprocess.run([sys.executable, SCRIPT, '--from-raw'], env=env, check=True, capture_output=True)
        first = (out / 'sepsweep.md').read_bytes()
    assert first == (out / 'sepsweep.md').read_bytes()


def test_report_dev_section_names_the_selected_configuration():
    raw = planted('dev', 6)
    txt = X.report({'dev': raw})
    assert '## 0. Development' in txt and 'Selected by the rule:' in txt and 'not yet frozen' in txt or 'Frozen in the script' in txt
    assert '| 1 | A05 | P | 0 |' in txt
    assert '## 1. Design and size' not in txt


# ---- checkpoints, shards, the raw ------------------------------------------------------------------------------------------------------------------------
def test_pending_rows_and_shards():
    allr = X.pending_rows('test')
    assert len(allr) == 7680 and allr[0] == (2, 0, 0, 0, 0)
    s = [X.pending_rows('dev', None, (i, 3)) for i in range(3)]
    assert sum(len(x) for x in s) == 1152 and len(set(r for x in s for r in x)) == 1152
    assert all(X.CELLS.index(r[:4]) % 3 == i for i, x in enumerate(s) for r in x)
    one = X.pending_rows('dev', [0, 5])
    assert len(one) == 12 and {r[:4] for r in one} == {X.CELLS[0], X.CELLS[5]}


def test_run_priority_tiers_and_incomplete_cells():
    allr = X.pending_rows('test')
    tiers = [X.tier(r[:4]) for r in allr]
    assert tiers == sorted(tiers) and set(tiers) == {0, 1, 2}
    assert X.tier((2, 0, 0, 0)) == 0 and X.tier((2, 0, 0, 3)) == 0 and X.tier((2, 0, 0, 2)) == 1 and X.tier((2, 0, 0, 1)) == 2
    first = allr[0]
    assert first[:4] == (2, 0, 0, 0) and [r[:4] for r in allr if X.tier(r[:4]) == 0][:2] == [(2, 0, 0, 0)] * 2
    s = X.pending_rows('test', None, (1, 3))
    assert all(X.CELLS.index(r[:4]) % 3 == 1 for r in s) and [X.tier(r[:4]) for r in s] == sorted(X.tier(r[:4]) for r in s)
    # a cell with fewer targets than the fullest makes its group incomplete
    raw = planted('test', 4)
    keep = np.ones(len(raw['test_rows']), bool)
    t0 = np.flatnonzero((raw['test_rows'][:, :4] == np.array([2, 4, 0, 0])).all(1))
    keep[t0[:2]] = False
    raw = {k: v[keep] for k, v in raw.items()}
    agg = X.aggregate(raw, 'test', ('A25', 'P', 0.5))
    assert (2, 0, 0) not in X.group_keys(agg) and (2, 0, 3) in X.group_keys(agg) and len(X.group_keys(agg)) == 31
    assert X.group_keys({}) == []


def test_assemble_orders_dedups_and_skips_absent_cells(tmp_path):
    rows = X.targets('dev', 2)
    sel = [rows[5], rows[1], rows[0], rows[200]]

    def rec(r, v):
        return dict(row=list(r), path=np.full((X.N_CFG, 4), v).tolist(), net=np.full((3, X.N_CFG, 5), v).tolist(), diag=[float(v)] * 17)
    f1, f2 = tmp_path / 'ckpt_dev_s0.jsonl', tmp_path / 'ckpt_dev_s1.jsonl'
    f1.write_text('\n'.join(json.dumps(rec(r, 1)) for r in sel[:3]) + '\n')
    f2.write_text('\n'.join(json.dumps(rec(r, 2)) for r in [sel[0], sel[3]]) + '\n')
    raw = X.assemble('dev', [str(f1), str(f2)])
    assert raw['dev_rows'].shape == (4, 5) and [tuple(r) for r in raw['dev_rows']] == [tuple(x) for x in sorted(sel, key=lambda r: rows.index(r))]
    assert raw['dev_path'].shape == (4, 32, 4) and raw['dev_net'].shape == (4, 3, 32, 5) and raw['dev_diag'].shape == (4, 17)
    assert raw['dev_path'][0, 0, 0] == 1 and raw['dev_path'][-1, 0, 0] == 2           # the first record of a duplicated target is kept
    assert raw['dev_path'].dtype == np.int16 and raw['dev_diag'].dtype == float


def test_save_and_load_raw_roundtrip(tmp_path):
    p = str(tmp_path / 'x.npz')
    X.save_raw(p, dict(a=np.arange(3)))
    assert np.array_equal(X.load_raw(p)['a'], np.arange(3)) and X.load_raw(str(tmp_path / 'none.npz')) == {}
    assert not os.path.exists(p + '.tmp.npz')


# ---- the stored raws (when present) ----------------------------------------------------------------------------------------------------------------------
STORED = os.path.exists(X.RAW) and os.path.exists(X.DEV_RAW)


@pytest.mark.skipif(not STORED, reason='raws not stored yet')
def test_stored_report_regenerates_byte_for_byte():
    assert X.FROZEN is not None
    raws = X.load_raws()
    assert X.report(raws) == open(X.MD).read()
    sel = X.select_config(raws['dev'])
    assert X.CFGS[sel] == X.FROZEN                                                   # the frozen configuration is the development rule's choice
    t = raws['test']
    assert t['test_rows'].shape[0] == t['test_path'].shape[0] == t['test_net'].shape[0] == t['test_diag'].shape[0]
    assert t['test_path'].shape[1:] == (32, 4) and t['test_net'].shape[1:] == (3, 32, 5)
    assert set(tuple(r[:4]) for r in t['test_rows'].tolist()) <= set(X.CELLS)
    d = raws['dev']
    assert d['dev_rows'].shape[0] <= 1152 and set(tuple(r[:4]) for r in d['dev_rows'].tolist()) <= set(X.CELLS)
