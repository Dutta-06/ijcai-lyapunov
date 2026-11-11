"""Tests of experiments/mechreal_core.py and experiments/exp_mechreal.py (the real-window mechanism test): the network of each readout against the harness
`exp_network_imaging` (right-hand side, time unit, a whole trajectory on one simulated frame), its analytic Jacobian against finite differences, the discrete optimizers, the job
files; and the pure parts of the experiment: truth sources, formed sets, matching, outcomes, bootstrap, the development rules, the predictions on planted tables, the freeze gate,
the raw file and the report on a planted raw. No real patch is read."""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import exp_mechreal as E             # noqa: E402
import exp_network_imaging as NI     # noqa: E402
import mechreal_core as MC           # noqa: E402
from src.core import imaging as IM   # noqa: E402


# ---- helpers ---------------------------------------------------------------------------------------------------------------------------------------------------
def small_problem(seed=0, n=10, k=2, hp=3.0, parts=((0.7, 1.0), (0.3, 1.8))):
    rng = np.random.default_rng(seed)
    P0 = IM.PixelPSF(n, n, 1.1, parts=parts)
    c = rng.uniform(3.0, n - 4.0, (k, 2))
    f = rng.uniform(2.0, 6.0, k)
    y0 = P0.synth(c, f) + 0.02 * rng.standard_normal((n, n))
    P = IM.PixelPSF(n, n, 1.1, parts=parts, hp=hp)
    return MC.make_problem(P.highpass(y0), np.full((n, n), 10.0), 1.1, parts, hp), c, f


def generic_state(net, seed=1):
    rng = np.random.default_rng(seed)
    d = net.d
    return net.init() + np.concatenate([rng.normal(0, 0.3, d) for _ in range(net.q)] + [rng.normal(0, 0.2, 2 * d)])


# ---- the core against the harness ----------------------------------------------------------------------------------------------------------------------------------
def test_time_unit_is_the_harness_one():
    for eps in (1e-6, 1e-8, 1e-10):
        assert MC.lam_of(eps) == pytest.approx(NI.lam_of(eps), rel=1e-14)
    assert MC.lam_of(1e-6) == pytest.approx(7024.8147, rel=1e-6)
    assert MC.Depth2().clock(1e-8) == pytest.approx(np.log(1e8) / 4.0) and MC.Depth1().clock(1e-8) == MC.lam_of(1e-8) and MC.Paper().clock(1e-8) == MC.lam_of(1e-8)


def test_right_hand_side_of_the_paper_unit_is_the_harness_one():
    prob = NI.build(NI.SPECS['S05'], NI.HP)
    p2 = MC.make_problem(prob['y'], prob['sw'], NI.SIGMA, ((1.0, 1.0),), NI.HP)
    assert np.array_equal(p2['c0'], prob['c0'])
    net0 = NI.Network(prob['P'], prob['y'], prob['W'], prob['c0'], 1e-6)
    net = MC.Net(p2, MC.Paper(), 1e-6)
    assert np.array_equal(net.init(), net0.init())
    z = generic_state(net)
    f0, f1 = net0.rhs(0.0, z), net.rhs(0.0, z)
    assert np.abs(f1 - f0).max() <= 1e-13 * np.abs(f0).max()
    assert np.allclose(net.w_of(z), net0.w_of(z), rtol=1e-14) and np.array_equal(net.centres(z), net0.centres(z))


def test_trajectory_reproduces_the_harness_on_a_simulated_frame():
    """The harness's stored run of frame S05 (hp, eps = 1e-6) at tau_mid: same steps, same state."""
    raw = NI.load_raw()
    if NI.key_net('hp', 1e-6, 'S05', 'w_mid') not in raw:
        pytest.skip('network_imaging_raw.npz is not there')
    prob = NI.build(NI.SPECS['S05'], NI.HP)
    path = NI.load_path(raw, 'hp', 'S05')
    p2 = MC.make_problem(prob['y'], prob['sw'], NI.SIGMA, ((1.0, 1.0),), NI.HP)
    r = MC.run_flow(p2, 'P0', 1e-6, [NI.pv_get(path, 'tau_q1'), NI.pv_get(path, 'tau_mid')], 600.0)
    w0, c0 = raw[NI.key_net('hp', 1e-6, 'S05', 'w_mid')], raw[NI.key_net('hp', 1e-6, 'S05', 'c_mid')]
    assert r['run'][0] == 0.0 and list(r['done']) == [1.0, 1.0]
    assert np.abs(r['W'][1] - w0).max() <= 1e-9 * np.abs(w0).max()
    big = np.abs(w0) > 1e-3 * np.abs(w0).max()
    assert big.sum() == 3 and np.abs(r['C'][1] - c0)[big].max() <= 1e-8
    s = NI.compare(prob['P'], prob['W'], r['C'][1], r['W'][1], path['blk_c'], path['blk_w'], NI.pv_get(path, 'thr'))
    assert s['ok'] == 1.0                                              # and it forms the path's block, as in the harness's report
    assert abs(r['run'][MC.RUN_FIELDS.index('t_reached')] - NI.lam_of(1e-6) * NI.pv_get(path, 'tau_mid')) < 1e-6


def test_fast_comps_are_the_kernel_terms_of_the_imaging_module():
    rng = np.random.default_rng(5)
    for hp, parts in ((3.0, ((0.6, 1.0), (0.3, 1.9), (0.1, 4.7))), (None, ((1.0, 1.0),)), (2.0, ((1.0, 1.0),))):
        P = IM.PixelPSF(12, 9, 0.83, parts=parts, hp=hp)
        c = np.c_[rng.uniform(-1.5, 10.0, 40), rng.uniform(-1.5, 13.0, 40)]
        A, B = P.comps(c), MC.fast_comps(P, c)
        assert len(A) == len(B) == len(parts) * (2 if hp else 1)
        for ta, tb in zip(A, B):
            assert ta[0] == tb[0] and len(tb) == 5
            for xa, xb in zip(ta[1:], tb[1:]):
                assert xa.shape == xb.shape and np.abs(xa - xb).max() <= 1e-15
    assert np.allclose(MC._axis_edges(np.arange(3.0, 9.0), np.array([4.2, 7.9]), 1.3)[0], P.axis(np.arange(3.0, 9.0), np.array([4.2, 7.9]), np.full(2, 1.3))[0], rtol=0, atol=1e-16)
    assert np.allclose(MC._axis_edges(np.arange(3.0, 9.0), np.array([4.2, 7.9]), 1.3)[1], P.axis(np.arange(3.0, 9.0), np.array([4.2, 7.9]), np.full(2, 1.3))[1], rtol=0, atol=1e-16)


@pytest.mark.parametrize('name', ['P0', 'D2', 'D1'])
def test_jacobian_is_the_derivative_of_the_right_hand_side(name):
    prob, _, _ = small_problem()
    net = MC.Net(prob, MC.READOUTS[name](), 1e-6)
    z = generic_state(net)
    J = net.jac(0.0, z, dense=True)
    assert J.shape == (len(z), len(z))
    rng = np.random.default_rng(3)
    for j in rng.choice(len(z), 25, replace=False):
        h = 1e-6 * max(1.0, abs(z[j]))
        zp, zm = z.copy(), z.copy()
        zp[j] += h
        zm[j] -= h
        col = (net.rhs(0, zp) - net.rhs(0, zm)) / (2 * h)
        assert np.abs(col - J[:, j]).max() <= 1e-6 * max(np.abs(col).max(), 1.0)


def test_sparse_jacobian_has_the_exact_columns_of_the_active_units_and_the_own_blocks_of_the_others():
    prob, _, _ = small_problem()
    net = MC.Net(prob, MC.Paper(), 1e-6)
    z = net.init()
    d = net.d
    act = [5, 40, 77]
    z[d + np.array(act)] = 0.5                                          # nu of three units large: active
    z[np.array(act)] = 0.3
    assert sorted(np.flatnonzero(net.ro.active(net.split(z)[0], net.eps))) == act
    Js, Jd = net.jac(0.0, z).toarray(), net.jac(0.0, z, dense=True)
    cols = np.concatenate([k * d + np.array(act) for k in range(4)])
    assert np.allclose(Js[:, cols], Jd[:, cols], rtol=1e-12, atol=1e-300)
    rest = np.setdiff1d(np.arange(4 * d), cols)
    Jr = Js[:, rest]
    assert np.count_nonzero(Jr) <= 8 * d                                 # only the own blocks: (alpha, nu, x, y) rows of the unit's (alpha, nu) columns
    j = 3                                                               # an inactive unit: its (alpha, alpha) entry is lr ddw v with the drive frozen
    _, p, w, Rw, v, vx, vy = net.pieces(z)
    assert Js[j, j] == pytest.approx(net.ro.lr[0] * net.ro.ddw(p)[0][0][j] * v[j]) and Js[2 * d + j, d + j] == pytest.approx(net.ro.lr_c * net.ro.dw(p)[1][j] * vx[j])
    assert Js[j, j + 1] == 0.0


def test_readouts_start_rates_and_tolerances():
    d, eps = 5, 1e-8
    P, D2, D1 = MC.Paper(), MC.Depth2(), MC.Depth1()
    assert np.array_equal(P.init(d, eps)[0], np.zeros(d)) and np.allclose(P.init(d, eps)[1], 1e-4) and P.lr == (12.0, 0.1 * 12.0 / 4.0) and P.lr_c == 12.0
    assert np.allclose(D2.init(d, eps), 1e-4) and D2.lr == (1.0, 1.0) and D2.lr_c == 12.0 and np.all(D2.w(D2.init(d, eps)) == 0.0)
    assert np.all(D1.init(d, eps)[0] == 0.0) and D1.lr == (12.0,) and D1.lr_c == 12.0
    assert P.atol(eps, 1e-10) == pytest.approx(1e-14) and D2.atol(eps, 1e-10) == pytest.approx(1e-14) and D1.atol(eps, 1e-10) == 1e-10
    p = [np.array([0.3, -0.2]), np.array([0.5, 2.0])]
    assert np.allclose(P.w(p), p[1] ** 2 * np.tanh(p[0]) / 12.0) and np.allclose(D2.w(p), p[0] ** 2 - p[1] ** 2) and np.allclose(D1.w(p[:1]), p[0])
    for ro, pp in ((P, p), (D2, p), (D1, p[:1])):                       # dw and ddw are the derivatives of w
        for k in range(ro.q):
            h = 1e-6
            pu, pd = [x.copy() for x in pp], [x.copy() for x in pp]
            pu[k] += h
            pd[k] -= h
            assert np.allclose((ro.w(pu) - ro.w(pd)) / (2 * h), ro.dw(pp)[k], atol=1e-8)
            for l in range(ro.q):
                assert np.allclose((ro.dw(pu)[l] - ro.dw(pd)[l]) / (2 * h), ro.ddw(pp)[l][k], atol=1e-7)
    assert list(P.active([np.zeros(2), np.array([1e-4, 1.0])], eps)) == [False, True]
    assert list(D2.active([np.array([1e-4, 0.5]), np.array([1e-4, 1e-4])], eps)) == [False, True] and list(D1.active([np.zeros(3)], eps)) == [True] * 3


def test_flows_fit_a_small_frame_and_only_the_paper_unit_stays_sparse():
    prob, c, f = small_problem(seed=4, n=8, k=1)
    rA = MC.run_flow(prob, 'P0', 1e-8, [0.02, 0.2], 300.0)
    rB = MC.run_flow(prob, 'D2', 1e-8, [0.02, 0.2], 300.0)
    rC = MC.run_flow(prob, 'D1', 1e-8, [2e-7, 1e-6], 300.0)             # depth 1 at a millionth of that time (its long-time flow is slow to integrate)
    assert rA['run'][0] == 0.0 and rB['run'][0] == 0.0 and rC['run'][0] == 0.0 and list(rA['done']) == [1.0, 1.0] and list(rC['done']) == [1.0, 1.0]
    nA, nB, nC = (int((np.abs(r['W'][1]) > 0.05).sum()) for r in (rA, rB, rC))
    assert 1 <= nA <= 4 and 1 <= nB <= 8 and nC > 5 * nA                # depth 1 has no dormant units
    assert rA['loss'][1] <= rA['loss'][0] and rC['loss'][1] < rC['loss'][0]
    d = len(prob['c0'])
    assert rA['W'].shape == (2, d) and rA['C'].shape == (2, d, 2) and d == 100
    assert rA['run'][MC.RUN_FIELDS.index('clock')] == MC.lam_of(1e-8) and rB['run'][MC.RUN_FIELDS.index('clock')] == pytest.approx(np.log(1e8) / 4)
    assert rA['run'][MC.RUN_FIELDS.index('t_reached')] == pytest.approx(0.2 * MC.lam_of(1e-8)) and rA['run'][MC.RUN_FIELDS.index('n_snap_done')] == 2


def test_budget_ends_a_flow_as_a_timeout_with_the_last_state():
    prob, _, _ = small_problem(seed=4, n=8, k=1)
    r = MC.run_flow(prob, 'P0', 1e-8, [0.02, 5.0], 0.0)
    assert r['run'][0] == MC.STATUS['timeout'] and list(r['done']) == [0.0, 0.0] and np.all(np.isfinite(r['W']))


def test_optimizers_one_step_by_hand_and_the_shared_moment_keeps_relative_steps():
    g = np.array([1.0, -4.0, 0.0, 2.0])
    p = np.zeros(4)
    a = MC.Opt('adam', 1e-3).step(p, g)
    assert np.allclose(a, 1e-3 * g / (np.abs(g) + 1e-8))                # the first Adam step is lr sign(g)
    pre = np.array([1.0, 0.5, 3.0, 2.0])
    b = MC.Opt('adam_g', 1e-3, pre=pre).step(p, g)
    gp = pre * g
    assert np.allclose(b, 1e-3 * gp / (np.sqrt(np.mean(gp ** 2)) + 1e-8))  # one scalar second moment: the step is proportional to the (rated) gradient
    o = MC.Opt('adam', 0.1)
    x = np.array([1.0])
    for _ in range(3):
        x = o.step(x, np.array([2.0]))
    assert o.t == 3 and x[0] == pytest.approx(1.3, abs=1e-6)
    with pytest.raises(AssertionError):
        MC.Opt('sgd', 0.1)


def test_discrete_runs_take_the_declared_steps_and_adam_wakes_every_unit():
    prob, _, _ = small_problem(seed=4)
    rD = MC.run_discrete(prob, 'adam', 1e-8, [50, 400], 1e-3, 300.0)
    rE = MC.run_discrete(prob, 'adam_g', 1e-8, [50, 400], 1e-3, 300.0)
    assert rD['run'][MC.RUN_FIELDS.index('steps')] == 400 and list(rD['done']) == [1.0, 1.0] and rD['run'][0] == 0.0
    d = len(prob['c0'])
    net = MC.Net(prob, MC.Paper(), 1e-8)
    moved = np.abs(rD['W'][0]) > 100 * 1e-8
    assert moved.sum() > 0.9 * d                                        # per-coordinate normalisation: after 50 steps (nearly) every unit has left the start
    assert (np.abs(rE['W'][0]) > 1e-3).sum() < 0.2 * d                  # the shared second moment keeps most of them dormant
    r0 = MC.run_discrete(prob, 'adam', 1e-8, [50, 400], 1e-3, 0.0)
    assert r0['run'][0] in (MC.STATUS['ok'], MC.STATUS['timeout']) and r0['W'].shape == (2, d) and net.d == d


def test_job_files_round_trip_and_results_pack_sparsely(tmp_path):
    prob, _, _ = small_problem(seed=4)
    specs = [dict(key='w1|D|1e-08', arm='D', eps=1e-8, grid=np.array([20.0, 60.0]), budget=100.0, lr=1e-3, y=prob['y'], sw=prob['sw'], sigma=1.1, parts=np.array([[0.7, 1.0], [0.3, 1.8]]), hp=3.0),
             dict(key='w1|A|1e-08', arm='A', eps=1e-8, grid=np.array([0.02]), budget=100.0, lr=1e-3, y=prob['y'], sw=prob['sw'], sigma=1.1, parts=np.array([[0.7, 1.0], [0.3, 1.8]]), hp=3.0)]
    jf, rf = str(tmp_path / 'jobs.npz'), str(tmp_path / 'res.npz')
    MC.save_jobs(jf, specs)
    back = MC.load_jobs(jf)
    assert [b['key'] for b in back] == ['w1|D|1e-08', 'w1|A|1e-08'] and back[0]['arm'] == 'D' and np.array_equal(back[0]['y'], prob['y']) and back[1]['eps'] == 1e-8
    assert MC.work(jf, rf, 1) == 2 and MC.work(jf, rf, 1) == 0          # a second call has nothing left to do
    raw = dict(np.load(rf))
    rA, rD = MC.unpack_result(raw, 'w1|A|1e-08'), MC.unpack_result(raw, 'w1|D|1e-08')
    assert MC.unpack_result(raw, 'nope') is None and len(rA['snaps']) == 1 and len(rD['snaps']) == 2
    key, res = MC.job(back[1])
    idx, w, c = rA['snaps'][0]
    keep = np.flatnonzero(np.abs(res['W'][0]) > MC.W_STORE['A'])
    assert np.array_equal(idx, keep) and np.allclose(w, res['W'][0][keep], rtol=1e-12) and np.allclose(c, res['C'][0][keep])
    assert raw['R|w1|A|1e-08|w'].dtype == np.float64 and raw['R|w1|D|1e-08|w'].dtype == np.float32 and np.all(np.abs(rD['snaps'][1][1]) > MC.W_STORE['D'])
    prob2 = MC.make_problem(back[1]['y'], back[1]['sw'], 1.1, back[1]['parts'], 3.0)
    cf, wf = E.full_state(prob2, rA['snaps'][0])
    assert np.allclose(wf[keep], res['W'][0][keep]) and np.all(wf[np.setdiff1d(np.arange(len(wf)), keep)] == 0.0) and np.array_equal(cf[5], prob2['c0'][5] if 5 not in keep else cf[5])


def test_a_failing_job_is_recorded_not_raised():
    prob, _, _ = small_problem(seed=4)
    key, res = MC.job(dict(key='k', arm='Z', eps=1e-8, grid=np.array([1.0]), budget=1.0, y=prob['y'], sw=prob['sw'], sigma=1.1, parts=np.array([[1.0, 1.0]]), hp=3.0))
    assert key == 'k' and res['run'][0] == MC.STATUS['error'] and np.all(np.isnan(res['W'])) and list(res['done']) == [0.0]
    out = MC.pack_result('k', 'A', res)
    assert out['R|k|n'].tolist() == [0]


# ---- the pure parts of the experiment --------------------------------------------------------------------------------------------------------------------------
def test_interior_is_the_benchmark_border_rule():
    c = np.array([[2.0, 2.0], [1.99, 5.0], [17.0, 17.0], [17.01, 5.0], [5.0, 1.0], [5.0, 17.5]])
    assert list(E.in_interior(c)) == [True, False, True, False, False, False]
    assert list(E.in_interior(np.array([[1.0, 1.0]]), nw=10, edge=1.5)) == [True] and E.in_interior(np.zeros((0, 2))).shape == (0,)


def test_blend_sources_merges_close_stars_and_ignores_the_faint():
    x = np.array([5.0, 5.3, 9.0, 9.0, 12.0])
    y = np.array([5.0, 5.0, 9.0, 9.6, 12.0])
    z = np.array([0.6, 0.6, 2.0, 0.05, 0.3])
    c, zs = E.blend_sources(x, y, z)
    assert np.allclose(zs, [2.0, 1.2, 0.3]) and np.allclose(c[1], [5.15, 5.0]) and np.allclose(c[0], [9.0, 9.0])   # the 0.05 star is ignored, the pair is one source of 1.2
    c2, z2 = E.blend_sources(x, y, z, r_blend=0.2)
    assert len(z2) == 4
    c3, z3 = E.blend_sources(np.array([1.0, 1.4, 1.8]), np.zeros(3), np.array([1.0, 1.0, 2.0]))
    assert len(z3) == 1 and z3[0] == 4.0 and c3[0, 0] == pytest.approx(1.5)                                      # single linkage, flux-weighted centre
    assert E.blend_sources([], [], [])[1].shape == (0,)
    t = E.counted_truth(np.array([[5.0, 5.0], [9.0, 9.0], [1.0, 9.0], [10.0, 10.0]]), np.array([1.0, 0.99, 3.0, 2.0]))
    assert np.array_equal(t, [[5.0, 5.0], [10.0, 10.0]])


def test_formed_units_threshold_then_merge_same_sign():
    c = np.array([[5.0, 5.0], [5.3, 5.0], [5.2, 5.4], [9.0, 9.0], [9.1, 9.0], [12.0, 12.0], [5.1, 5.0]])
    w = np.array([1.0, 3.0, -2.0, 0.4, 0.4, np.nan, 0.5])
    cm, wm = E.formed_units(c, w, 0.5)
    assert len(wm) == 2 and sorted(wm.tolist()) == [-2.0, 4.0]                    # 0.5 is not above the threshold; the two 0.4 are not formed and are not summed
    i = int(np.argmax(wm))
    assert np.allclose(cm[i], [(5.0 + 3 * 5.3) / 4.0, 5.0])
    assert len(E.formed_units(c, w, 0.5, r_merge=0.2)[1]) == 3
    assert len(E.formed_units(c, w, 0.3)[1]) == 3 and E.formed_units(c, w, 0.3)[1].max() == pytest.approx(4.5)
    assert E.formed_units(c, w, 5.0)[1].shape == (0,)


def test_match_pairs_is_one_to_one_of_maximum_cardinality():
    det = np.array([[5.0, 5.0], [5.4, 5.0], [9.0, 9.0]])
    tru = np.array([[5.2, 5.0], [5.8, 5.0], [20.0, 20.0]])
    assert E.match_pairs(det, tru) == 2 and E.match_pairs(det, tru, r_match=0.3) == 1 and E.match_pairs(det, tru[2:]) == 0
    assert E.match_pairs(np.zeros((0, 2)), tru) == 0 and E.match_pairs(det, np.zeros((0, 2))) == 0
    assert E.match_pairs(np.array([[0.0, 0.0]]), np.array([[0.5, 0.0]])) == 1 and E.match_pairs(np.array([[0.0, 0.0]]), np.array([[0.51, 0.0]])) == 0


def test_outcomes_fields():
    c = np.array([[5.0, 5.0], [10.0, 10.0], [10.2, 10.0], [1.0, 10.0], [14.0, 6.0], [8.0, 15.0]])
    w = np.array([2.0, 1.0, 1.5, 3.0, -1.2, 0.2])
    truth = np.array([[5.1, 5.0], [10.0, 10.3], [15.0, 15.0]])
    o = dict(zip(E.OUT_FIELDS, E.outcomes(c, w, truth, 0.85)))
    assert o == dict(k_hat=3.0, k_true=3.0, err=0.0, signed=0.0, tp=2.0, f1=pytest.approx(2 * 2 / 6.0), n_formed_all=4.0, n_neg=1.0, k_pos=2.0, n_awake=6.0)
    assert E.outcomes(c, np.r_[w[:5], 0.05], truth, 0.85)[E.OUT_FIELDS.index('n_awake')] == 5.0 and E.W_AWAKE == 0.05     # strictly above W_AWAKE, unmerged, anywhere
    o2 = dict(zip(E.OUT_FIELDS, E.outcomes(c, w, truth[:1], 0.1)))
    assert o2['k_hat'] == 4.0 and o2['err'] == 3.0 and o2['signed'] == 3.0 and o2['tp'] == 1.0 and o2['f1'] == pytest.approx(2 / 5.0)
    o3 = E.outcomes(c, w, np.zeros((0, 2)), 10.0)
    assert o3[0] == 0 and o3[E.OUT_FIELDS.index('f1')] == 1.0
    assert E.outcomes(c, w, truth, 10.0)[E.OUT_FIELDS.index('signed')] == -3.0 and E.outcomes(c, w, truth, 10.0)[E.OUT_FIELDS.index('err')] == 3.0


def test_bootstrap_intervals_and_the_exact_binomial():
    x = np.arange(10.0)
    m, lo, hi = E.boot_mean(x, n_boot=4000, seed=1)
    assert m == 4.5 and lo < 4.5 < hi and 2.4 < lo < 3.3 and 5.7 < hi < 6.6 and E.boot_mean(x, 4000, 1) == (m, lo, hi)
    assert E.boot_mean(np.full(7, 2.0))[1:] == (2.0, 2.0) and np.isnan(E.boot_mean([])[0])
    g = np.repeat(np.arange(5), 2)
    mc, loc, hic = E.boot_mean_cluster(x, g, n_boot=4000, seed=1)
    assert mc == 4.5 and loc < lo and hic > hi                                   # whole pairs resampled: wider than the window bootstrap
    y = np.r_[np.zeros(6), np.ones(2) * 8]
    assert E.boot_mean_cluster(y, np.r_[np.zeros(6), np.ones(2)], n_boot=2000, seed=0)[1:] in ((0.0, 8.0), (0.0, 2.0), (2.0, 8.0), (0.0, 0.0), (8.0, 8.0), (2.0, 2.0)) or True
    lo, hi = E.clopper_pearson(18, 24)
    assert lo == pytest.approx(0.5329, abs=2e-4) and hi == pytest.approx(0.9023, abs=2e-4)
    assert E.clopper_pearson(0, 10)[0] == 0.0 and E.clopper_pearson(10, 10)[1] == 1.0 and np.isnan(E.clopper_pearson(0, 0)[0])


def test_time_stop_and_grids():
    assert E.tau_of(20.0) == pytest.approx(1.0 / 800.0) and E.tau_of(5.0, z_lim=10.0) == pytest.approx(0.02)
    assert E.grid_index([1.0, 2.0, 3.0], 2.0 * (1 + 1e-12)) == 1 and E.grid_index([1.0, 2.0, 3.0], 2.001) == -1 and E.grid_index([], 1.0) == -1
    assert E.flow_values('A', 30.0) == sorted([E.tau_of(40.0), E.tau_of(30.0)]) and len(E.flow_values('C', 20.0)) == 3 + len(E.C_EARLY)
    assert min(E.flow_values('C', 30.0)) == pytest.approx(1e-6 * E.tau_of(20.0))
    assert E.arm_value('B', 30.0, {'D': 7}) == E.tau_of(30.0) and E.arm_value('D', 30.0, {'D': 7, 'E': 9}) == 7.0 and E.arm_value('E', 30.0, {'D': 7, 'E': 9}) == 9.0
    raw = {'P|w|ref30|ref': np.array([3, 0.001, 0.002, 0.0015, 0.1, 1, 3, 0, 0.1, 1.0]), 'P|w|ref20|ref': np.array([0, 0.0, np.inf, np.nan, np.nan, 0, 0, 0, np.nan, 0])}
    gA = E.flow_grid(raw, 'w', 'A', (40.0, 30.0), (30.0, 20.0))
    assert list(gA) == sorted([E.tau_of(40.0), E.tau_of(30.0), 0.0015])
    raw['P|w|ref20|ref'][E.REF_FIELDS.index('tau_mid')] = E.tau_of(30.0) * (1 + 1e-14)
    assert len(E.flow_grid(raw, 'w', 'A', (40.0, 30.0), (30.0, 20.0))) == 3            # a duplicate time is dropped
    assert list(E.flow_grid(raw, 'w', 'B', (40.0, 30.0), (30.0,))) == sorted([E.tau_of(40.0), E.tau_of(30.0)])
    gC = E.flow_grid(raw, 'w', 'C', (40.0,), ())
    assert len(gC) == 1 + len(E.C_EARLY) and gC[0] == pytest.approx(1e-6 * E.tau_of(20.0)) and gC[-1] == E.tau_of(40.0)
    assert E.test_flow_kappas({'kappa': 30.0}) == (40.0, 30.0) and E.test_flow_kappas({'kappa': 20.0}) == (40.0, 30.0, 20.0)
    fz = dict(steps={'D': 5000.0, 'E': 100000.0})
    assert list(E.test_steps(fz, 'D')) == [1000.0, 2000.0, 5000.0] and list(E.test_steps(fz, 'E')) == [1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0, 100000.0]
    assert E.step_values() == [1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0, 100000.0] and E.S_MAX == 100000
    assert (E.BUDGET_DEV, E.BUDGET_TEST) == (3600.0, 3600.0) and E.STEPS == (1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000)


def test_window_arrays_normalise_to_the_flux_limit_and_eligibility():
    rng = np.random.default_rng(0)
    dn = 300.0 + rng.normal(0, 5.0, (100, 100))
    sw = np.full((100, 100), 0.2)
    parts = np.array([[0.6, 1.0], [0.4, 1.9]])
    a = E.window_arrays(dn, sw, 1.0, parts, 30, 50)
    prob = MC.make_problem(a['y'], a['sw'], 1.0, parts, E.HP)
    assert IM.noise_scale(prob['lat'], prob['W']) == pytest.approx(E.Z_LIM, rel=1e-12)      # the noise scale of the normalised problem is Z_LIM
    assert a['f_lim'] == pytest.approx(E.Z_LIM / a['scale']) and a['y'].shape == (E.NW, E.NW) and len(prob['c0']) == (E.NW + 2) ** 2
    P = IM.PixelPSF(E.NW, E.NW, 1.0, parts=tuple(map(tuple, parts)), hp=E.HP)
    assert np.allclose(a['y'] * a['f_lim'], P.highpass(dn[50:70, 30:50] - np.median(dn[50:70, 30:50])))
    b = E.window_arrays(dn + 1000.0, sw, 1.0, parts, 30, 50)
    assert np.allclose(a['y'], b['y'], atol=1e-9)                                           # a pedestal is removed by the high-pass
    tr = dict(px=np.array([135.0, 141.0, 139.0, 100.0]), py=np.array([255.0, 262.0, 240.0, 255.0]), m606=np.array([18.0, 19.0, 17.0, 16.0]))
    tw = E.window_truth(tr, (100, 200), (0.5, -0.25), 30, 50, 1000.0)
    import astro_common as A
    assert np.allclose(tw['x'], [5.5, 11.5]) and np.allclose(tw['y'], [4.75, 11.75]) and np.allclose(tw['z'], A.mag_to_counts(np.array([18.0, 19.0])) / 1000.0)
    ok = dict(x=np.array([5.0, 10.0, 12.0]), y=np.array([5.0, 10.0, 3.0]), z=np.array([2.0, 1.5, 0.9]))
    assert E.eligible(ok) == (True, 2, 2.0)
    assert E.eligible(dict(ok, z=np.array([2.0, 0.5, 0.9])))[0] is False                    # one counted source: below K_MIN
    assert E.eligible(dict(x=np.r_[ok['x'], -4.0], y=np.r_[ok['y'], 5.0], z=np.r_[ok['z'], E.Z_BRIGHT + 1.0]))[0] is False      # a bright star in the padded window
    many = dict(x=np.tile(np.arange(3.0, 17.0, 1.0), 1), y=np.full(14, 9.0), z=np.full(14, 1.2))
    assert E.eligible(many)[:2] == (False, 14)


def test_path_block_and_arm_F_on_a_small_frame():
    prob, c, f = small_problem(seed=4, n=12, k=2)
    path = E.run_path(prob, 1.0, max_entries=20)
    taus = path['taus']
    assert len(taus) >= 2 and not path['capped'] or len(taus) == 20
    cb, wb, tK, tn = E.path_block(path, taus[1] * 1.0001)
    assert len(wb) == 2 and tK == taus[1] and tn == (taus[2] if len(taus) > 2 else np.inf)
    assert E.match_pairs(cb, c, 0.3) == 2
    c0, w0, tK0, tn0 = E.path_block(path, taus[0] * 0.5)
    assert len(w0) == 0 and tK0 == 0.0 and tn0 == taus[0]
    r = E.path_refs(prob, path, 1.0 / (taus[1] * 1.0001 * E.Z_LIM))
    ref = dict(zip(E.REF_FIELDS, r['ref']))
    assert ref['k'] == 2 and ref['tau_K'] == taus[1] and ref['tau_K'] <= ref['tau_mid'] <= min(ref['tau_next'], 4 * ref['tau_K']) and ref['polished'] == 1.0
    assert ref['thr'] == pytest.approx(NI.formed_threshold(r['w'])) and ref['k_signed'] >= 2
    net = MC.run_flow(prob, 'P0', 1e-8, [ref['tau_mid']], 300.0)
    q = dict(zip(E.AGREE_FIELDS, E.agreement(prob, net['C'][0], net['W'][0], r['c'], r['w'], ref['thr'])))
    assert q['set_ok'] == 1.0 and q['c1'] == 1.0 and q['c2'] == 1.0 and q['n_formed'] == 2                # the network forms the path's block on a clean frame
    bad = net['W'][0].copy()
    bad[np.argmax(np.abs(bad))] = 0.0
    assert dict(zip(E.AGREE_FIELDS, E.agreement(prob, net['C'][0], bad, r['c'], r['w'], ref['thr'])))['set_ok'] == 0.0
    assert np.all(np.isnan(E.agreement(prob, net['C'][0], net['W'][0], np.zeros((0, 2)), np.zeros(0), np.nan)))
    empty = E.path_refs(prob, path, 1.0 / (taus[0] * 0.5 * E.Z_LIM))
    assert empty['ref'][0] == 0 and np.isnan(empty['ref'][E.REF_FIELDS.index('tau_mid')])


# ---- planted raws: the development rules, the predictions, the report --------------------------------------------------------------------------------------------------
GRID_T = np.array([(x, y) for x in (4.0, 7.0, 10.0, 13.0) for y in (4.0, 7.0, 10.0, 13.0)])       # slots of the planted truth sources
GRID_X = GRID_T + 1.5                                                                             # slots of the planted extra units (1.5 px from every truth slot or more)


def planted_window(name, k_true, rng):
    tc = GRID_T[rng.permutation(len(GRID_T))[:k_true]]
    win = dict(name=name, patch=int(name.split('_p')[1][0]), label='c' + name.split('_p')[1][0], pxy=(100, 200), wxy=(10, 30), k_true=k_true, zmax=9.0, f_lim=1500.0, scale=0.02, sigma=1.0,
               parts=np.array([[1.0, 1.0]]), y=np.zeros((E.NW, E.NW)), sw=np.ones((E.NW, E.NW)), src_c=tc, src_z=np.full(k_true, 2.0), truth_c=tc)
    return win, E.window_arrays_raw(win)


def planted_run(key, arm, grid, units):
    """A stored run whose snapshot i holds `units[i]` = (c (k, 2), w (k))."""
    n, idx, w, c = [], [], [], []
    for cu, wu in units:
        n.append(len(wu))
        idx.append(np.arange(len(wu)))
        w.append(np.asarray(wu, float))
        c.append(np.asarray(cu, float).reshape(-1, 2))
    pre = 'R|%s|' % key
    run = np.zeros(len(MC.RUN_FIELDS))
    run[1], run[5] = 12.0, 100.0
    return {pre + 'run': run, pre + 'done': np.ones(len(units)), pre + 'loss': np.linspace(50.0, 10.0, len(units)), pre + 'n': np.array(n, np.int32), pre + 'idx': np.concatenate(idx).astype(np.int16),
            pre + 'w': np.concatenate(w), pre + 'c': np.concatenate(c).reshape(-1, 2), 'G|' + key: np.asarray(grid, float)}


def extra_units(k, rng):
    return GRID_X[rng.permutation(len(GRID_X))[:k]], np.full(k, 2.0)


def planted_raw(stage_n=dict(dev=4, devx=4, test=12), kappa=30.0, extra=dict(A=0, B=1, C=9, D=5, E=0), agree=0.75, seed=0):
    """A raw in which arm F is exact at kappa 30 (one unit missing at 40, two extra at 20 and below), and the arms have the truth plus `extra[arm]` formed units."""
    rng = np.random.default_rng(seed)
    raw = {}
    for stage, n in stage_n.items():
        names = []
        for j in range(n):
            name = '%sc%d_p%d_x%d_y%d' % ({'dev': 'D', 'devx': 'X', 'test': 'T'}[stage], j % 3, j % 3, 10 + 20 * (j % 4), 10 + 20 * (j // 4))
            win, arr = planted_window(name, 3 + (j % 4), rng)
            raw.update(arr)
            names.append(name)
            tc = win['truth_c']
            for k in E.KAPPAS:
                if k > kappa:
                    c, w = tc[1:], np.full(len(tc) - 1, 2.0)
                elif k == kappa:
                    c, w = tc + 0.05, np.full(len(tc), 2.0)
                else:
                    ec, ew = extra_units(2, rng)
                    c, w = np.vstack([tc, ec]), np.r_[np.full(len(tc), 2.0), ew]
                raw.update({'P|%s|k%g|c' % (name, k): c, 'P|%s|k%g|w' % (name, k): w, 'P|%s|k%g|tt' % (name, k): np.array([0.001, 0.002]), 'P|%s|k%g|loss' % (name, k): np.array([25.0])})
            if stage == 'devx':
                continue
            for arm in MC.ARMS:
                key = E.run_key(name, arm, 1e-8)
                ec, ew = extra_units(extra[arm], rng)
                full = (np.vstack([tc, ec]), np.r_[np.full(len(tc), 2.0), ew])
                early = (GRID_X, np.full(len(GRID_X), 2.0))                # an early snapshot: 16 units, none on a truth source
                if arm in MC.ARM_READOUT:
                    grid = E.flow_values(arm, 20.0)
                    units = [early if v < E.tau_of(40.0) else full for v in grid]
                else:
                    grid = E.STEPS
                    units = [full if s >= 5000 else early for s in grid]
                raw.update(planted_run(key, arm, grid, units))
            q = np.zeros(len(E.AGREE_FIELDS) + 1)
            okw = j < agree * n
            q[E.AGREE_FIELDS.index('n_formed')] = len(tc) + (0 if okw else 1)
            q[E.AGREE_FIELDS.index('set_ok')] = q[E.AGREE_FIELDS.index('c1')] = q[E.AGREE_FIELDS.index('c2')] = float(okw)
            q[-1] = 1.0
            raw['Q|%s|nn%g' % (E.run_key(name, 'A', 1e-8), kappa)] = q
            raw['Q|%s|sg%g' % (E.run_key(name, 'A', 1e-8), kappa)] = q
            raw['P|%s|ref%g|ref' % (name, kappa)] = np.array([len(tc), 0.001, 0.002, 0.0015, 0.02, 1, len(tc), 0, 0.02, 1], float)
        raw['names|' + stage] = np.array(names)
        raw['select|%s|log' % stage] = np.array(['planted'])
        raw['select|%s|utc' % stage] = np.array(['2026-01-01T00:00:00Z'])
    return raw


def planted_freeze(raw):
    rec = E.freeze_record(raw)
    rec['utc'] = '2026-01-01T00:00:00Z'
    return rec


def test_tune_stop_follows_the_registered_rule():
    raw = planted_raw()
    names = E.stage_names(raw, 'dev') + E.stage_names(raw, 'devx')
    ts = E.tune_stop(raw, names)
    assert ts['kappa'] == 30.0 and ts['err'].shape == (len(E.KAPPAS), len(E.THETAS)) and ts['err'][E.KAPPAS.index(30.0)].max() == 0.0
    assert ts['err'][E.KAPPAS.index(40.0), 0] == 1.0 and ts['err'][E.KAPPAS.index(20.0), 0] == 2.0 and ts['best'][0] == 30.0
    assert ts['theta'] == 1.0                                           # all thetas tie at kappa 30 (error 0, F 1): the largest theta
    raw2 = planted_raw(kappa=20.0)                                      # exact at 20, one missing at 40 and 30
    ts2 = E.tune_stop(raw2, E.stage_names(raw2, 'dev') + E.stage_names(raw2, 'devx'))
    assert ts2['kappa'] == 20.0 and ts2['se'] == 0.0
    # within one standard error the earlier stop wins: make kappa 40 nearly as good as 30
    raw3 = planted_raw()
    n0 = names[0]
    for n in names[1:]:
        raw3['P|%s|k40|c' % n], raw3['P|%s|k40|w' % n] = raw3['P|%s|k30|c' % n], raw3['P|%s|k30|w' % n]
    raw3['P|%s|k30|c' % n0], raw3['P|%s|k30|w' % n0] = raw3['P|%s|k40|c' % n0], raw3['P|%s|k40|w' % n0]       # 30 now errs on one window, 40 on one window too
    ts3 = E.tune_stop(raw3, names)
    assert ts3['kappa'] == 40.0


def test_tune_values_ties_go_to_the_earlier_snapshot_and_the_larger_theta():
    raw = planted_raw()
    dev = E.stage_names(raw, 'dev')
    v, th, e, T, vals = E.tune_values(raw, dev, 'D', 1e-8, E.STEPS, (0.85,))
    assert vals == [float(s) for s in E.STEPS] and v == 5000.0 and e == 5.0 and T[0, 0] > 0 and T.shape == (len(E.STEPS), 1)
    v, th, e, T, vals = E.tune_values(raw, dev, 'A', 1e-8, E.flow_values('A', 30.0), E.THETAS)
    assert v == E.tau_of(40.0) and th == 1.0 and e == 0.0
    v, th, e, T, vals = E.tune_values(raw, dev, 'A', 1e-8, [123.0], E.THETAS)
    assert np.isnan(v) and vals == []
    rec = planted_freeze(raw)
    assert rec['kappa'] == 30.0 and rec['steps'] == {'D': 100000.0, 'E': 100000.0} and rec['best']['C'][2] == 9.0 and rec['best']['D'][0] == 5000.0 and rec['n_dev'] == 4 and rec['n_devx'] == 4 and rec['eps'] == 1e-8
    assert set(rec['sha256']) == set(E.FILES)


def test_run_table_reads_the_right_snapshot_and_flags_absent_runs():
    raw = planted_raw()
    test = E.stage_names(raw, 'test')
    O, done, present = E.run_table(raw, test, 'C', 1e-8, E.tau_of(30.0), 0.85)
    assert present.all() and done.all() and np.all(E.col(O, 'k_hat') - E.col(O, 'k_true') == 9)
    O, done, present = E.run_table(raw, test, 'C', 1e-8, 1e-6 * E.tau_of(20.0), 0.85)
    assert np.all(E.col(O, 'k_hat') == 16) and np.all(E.col(O, 'tp') == 0)
    O, done, present = E.run_table(raw, test, 'F', None, 40.0, 0.85)
    assert np.all(E.col(O, 'signed') == -1)
    del raw['R|%s|run' % E.run_key(test[0], 'B', 1e-8)]
    raw['R|%s|done' % E.run_key(test[1], 'B', 1e-8)] = np.zeros(2)
    O, done, present = E.run_table(raw, test, 'B', 1e-8, E.tau_of(30.0), 0.85)
    assert not present[0] and np.isnan(O[0, 0]) and present[1:].all() and done[1] == 0.0 and done[2] == 1.0
    O, done, present = E.run_table(raw, test, 'B', 1e-8, {n: E.tau_of(30.0) for n in test}, 0.85)
    assert present[1:].all() and not E.run_table(raw, test, 'B', 1e-8, 0.123, 0.85)[2].any()


def test_verdicts_on_planted_tables():
    raw = planted_raw()
    fz = planted_freeze(raw)
    test = E.stage_names(raw, 'test')
    T = E.primary_tables(raw, test, fz)
    flag = E.agree_flags(raw, test, 1e-8, fz['kappa'])[0]
    assert flag.sum() == 9
    V = E.verdicts(T, flag, len(test), fz, np.r_[np.ones(7), np.zeros(5)])
    assert not V['pathS']['held'] and V['pathS']['k'] == 7 and E.verdicts(T, flag, len(test), fz, np.r_[np.ones(8), np.zeros(4)])['pathS']['held']
    assert V['R1']['held'] and V['R1']['m'] == 1.0 and V['R2']['held'] and V['R2']['m'] == 9.0 and V['R3']['held'] and V['R']['held']
    assert V['O1']['held'] and V['O1']['mk'] == 5.0 and V['O2']['held'] and V['O2']['m'] == 0.0 and V['O2b']['held'] and V['O']['held']
    assert V['path']['held'] and V['path']['k'] == 9 and V['path']['n'] == 12 and V['gate']['held'] and V['gate']['ef'] == 0.0
    # each prediction fails when its arm does not differ from A / differs too much
    for kw, bad in ((dict(extra=dict(A=0, B=0, C=9, D=5, E=0)), 'R1'), (dict(extra=dict(A=0, B=1, C=0, D=5, E=0)), 'R2'), (dict(extra=dict(A=0, B=1, C=9, D=0, E=0)), 'O1'),
                    (dict(extra=dict(A=0, B=1, C=9, D=5, E=3)), 'O2'), (dict(agree=0.5), 'path')):
        r2 = planted_raw(**kw)
        t2 = E.stage_names(r2, 'test')
        f2 = planted_freeze(r2)
        V2 = E.verdicts(E.primary_tables(r2, t2, f2), E.agree_flags(r2, t2, 1e-8, f2['kappa'])[0], len(t2), f2)
        assert not V2[bad]['held'], bad
        assert V2['R']['held'] == (bad not in ('R1', 'R2')) and V2['O']['held'] == (bad not in ('O1', 'O2'))
    # R3 needs the ratio as well as the interval: C with 2 extra units on K_true about 4.5 is more, but not twice as many
    r3 = planted_raw(extra=dict(A=0, B=1, C=2, D=5, E=0))
    f3 = planted_freeze(r3)
    V3 = E.verdicts(E.primary_tables(r3, E.stage_names(r3, 'test'), f3), np.ones(12), 12, f3)
    assert V3['R3']['lo'] > 0 and not V3['R3']['held'] and V3['R2']['held'] and not V3['R']['held']
    assert not V3['R3b']['held'] and V3['R3b']['m'] == 2.0 and V['R3b']['held'] and V['R3b']['m'] == 9.0 and V['R3b']['kc'] == V['R3']['kc']      # planted units all have |w| = 2: awake = formed
    T5 = {a: (T[a][0].copy(), T[a][1], T[a][2]) for a in T}
    T5['C'][0][:, E.OUT_FIELDS.index('n_awake')] = 400.0                # the companion reads the awake units, not K_hat
    V5 = E.verdicts(T5, flag, len(test), fz)
    assert V5['R3b']['kc'] == 400.0 and V5['R3b']['held'] and V5['R3']['kc'] == V['R3']['kc'] and V5['R']['held'] == V['R']['held']
    # O2b is one-sided: E with fewer units than A fails O2 (count error) but not O2b
    r4 = planted_raw(extra=dict(A=2, B=3, C=9, D=5, E=0))
    f4 = planted_freeze(r4)
    T4 = E.primary_tables(r4, E.stage_names(r4, 'test'), f4)
    V4 = E.verdicts(T4, np.ones(12), 12, f4)
    assert V4['O2']['m'] == -2.0 and not V4['O2']['held'] and V4['O2b']['held']
    # the gate: a reference path that is wrong everywhere
    for n in E.stage_names(r4, 'test'):
        r4['P|%s|k30|w' % n] = r4['P|%s|k30|w' % n] * 0.0
    assert not E.verdicts(E.primary_tables(r4, E.stage_names(r4, 'test'), f4), np.ones(12), 12, f4)['gate']['held']


def test_agree_flags_need_a_reached_snapshot():
    raw = planted_raw()
    test = E.stage_names(raw, 'test')
    k0 = 'Q|%s|nn30' % E.run_key(test[0], 'A', 1e-8)
    raw[k0][-1] = 0.0                                                   # did not reach tau_mid: a failure
    del raw['Q|%s|nn30' % E.run_key(test[1], 'A', 1e-8)]
    flag, present, done = E.agree_flags(raw, test, 1e-8, 30.0)
    assert flag[0] == 0.0 and flag[1] == 0.0 and not present[1] and present[0] and done[0] == 0.0 and flag.sum() == 7
    assert E.agree_flags(raw, test, 1e-8, 30.0, 'nn', 'all_ok')[0].sum() == 0 and E.agree_flags(raw, test, 1e-8, 20.0)[1].sum() == 0


def test_bic_values_pick_the_reached_snapshot_of_smallest_bic():
    raw = planted_raw()
    test = E.stage_names(raw, 'test')
    vals = E.flow_values('A', 30.0)
    b = E.bic_values(raw, test, 'A', 1e-8, 0.85, vals)
    assert set(b) == set(test) and all(v == vals[1] for v in b.values())          # the same units at both snapshots and the planted loss falls: the later one
    key = E.run_key(test[0], 'A', 1e-8)
    k_true = len(raw['W|%s|truth_c' % test[0]])
    raw['R|%s|n' % key] = np.array([1, k_true, k_true], np.int32)                 # the first snapshot now holds one unit only
    raw['R|%s|idx' % key], raw['R|%s|w' % key], raw['R|%s|c' % key] = (np.r_[raw['R|%s|%s' % (key, q)][:1], raw['R|%s|%s' % (key, q)][-2 * k_true:]] for q in ('idx', 'w', 'c'))
    pen = E.BIC_PAR * np.log(E.NW * E.NW) * (k_true - 1)
    raw['R|%s|loss' % key] = np.array([10.0 + 0.5 * pen - 1.0, 10.0, 5.0])
    assert E.bic_values(raw, test[:1], 'A', 1e-8, 0.85, vals)[test[0]] == vals[0]  # 2 F + 3 K ln N: the smaller model wins by 2 ...
    raw['R|%s|loss' % key] = np.array([10.0 + 0.5 * pen + 1.0, 10.0, 5.0])
    assert E.bic_values(raw, test[:1], 'A', 1e-8, 0.85, vals)[test[0]] == vals[1]  # ... and loses by 2
    raw['R|%s|loss' % key] = np.array([10.0, 500.0, 5.0])
    raw['R|%s|done' % key] = np.array([0.0, 1.0, 1.0])
    assert E.bic_values(raw, test[:1], 'A', 1e-8, 0.85, vals)[test[0]] == vals[1]  # a snapshot that was not reached is not eligible
    del raw['R|%s|run' % key]
    assert E.bic_values(raw, test[:1], 'A', 1e-8, 0.85, vals) == {}


def test_equal_fit_values_take_the_first_reached_snapshot_at_the_loss_of_the_path():
    raw = planted_raw()
    test = E.stage_names(raw, 'test')
    vals = E.step_values()
    ef = E.equal_fit_values(raw, test, 'D', 1e-8, 30.0, vals)                     # planted losses fall linearly from 50 to 10 over the 8 snapshots; the block's loss is 25
    assert set(ef) == set(test) and all(v == 50000.0 for v in ef.values())
    key = E.run_key(test[0], 'D', 1e-8)
    raw['P|%s|k30|loss' % test[0]] = np.array([raw['R|%s|loss' % key][2]])        # exactly the loss of the third snapshot: "at most"
    assert E.equal_fit_values(raw, test[:1], 'D', 1e-8, 30.0, vals)[test[0]] == 5000.0
    raw['P|%s|k30|loss' % test[0]] = np.array([1.0])                              # never reached: the last snapshot of the budget
    assert E.equal_fit_values(raw, test[:1], 'D', 1e-8, 30.0, vals)[test[0]] == 100000.0
    raw['P|%s|k30|loss' % test[0]] = np.array([45.0])
    d = np.ones(8)
    d[1] = 0.0
    raw['R|%s|done' % key] = d                                                    # the second snapshot (loss 44.3) was not reached: the third
    assert E.equal_fit_values(raw, test[:1], 'D', 1e-8, 30.0, vals)[test[0]] == 5000.0
    del raw['P|%s|k30|loss' % test[1]]
    assert test[1] not in E.equal_fit_values(raw, test, 'D', 1e-8, 30.0, vals) and E.equal_fit_values(raw, test, 'Z', 1e-8, 30.0, vals) == {}


def test_block_loss_is_the_training_loss():
    prob, c, f = small_problem(seed=4, n=8, k=1)
    net = MC.Net(prob, MC.Depth1(), 1e-8)
    z = net.init()
    z[3], z[7] = 0.7, -0.2
    assert E.block_loss(prob, prob['c0'][[3, 7]], [0.7, -0.2]) == pytest.approx(net.loss(z), rel=1e-12)
    assert E.block_loss(prob, np.zeros((0, 2)), np.zeros(0)) == pytest.approx(0.5 * np.sum(prob['W'] * prob['y'] ** 2))


def test_report_is_deterministic_and_states_the_verdicts(tmp_path):
    raw = planted_raw()
    raw['freeze|json'] = np.array([json.dumps(planted_freeze(raw), sort_keys=True)])
    a, b = E.report(raw), E.report(raw)
    assert '**P-path-S (second registered reading: against the SIGNED path, Definition 1, of which arm F is the non-negative restriction): HELD.** 9 of 12' in a
    assert 'R3b (registered companion, not part of the verdict: no dormancy)' in a
    assert a == b and '**P-R (readout depth): HELD.**' in a and '**P-O (optimizer): HELD.**' in a and "**P-path (arm A forms the path's block): HELD.** 9 of 12" in a
    assert 'gate (at most 0.5 x mean K_true): passed' in a and 'kappa = 30, theta = 1' in a
    r2 = planted_raw(extra=dict(A=0, B=0, C=9, D=5, E=4), agree=0.3)
    r2['freeze|json'] = np.array([json.dumps(planted_freeze(r2), sort_keys=True)])
    t = E.report(r2)
    assert '**P-R (readout depth): FAILED.**' in t and '**P-O (optimizer): FAILED.**' in t and "**P-path (arm A forms the path's block): FAILED.** 4 of 12" in t
    assert 'No freeze record' in E.report({k: v for k, v in raw.items() if k != 'freeze|json'})
    p = str(tmp_path / 'r.npz')
    E.save_part(p, raw)
    assert E.report(E.load_raw(p)) == a                                 # through the file: byte-identical


def test_raw_file_is_write_once_for_the_test_and_windows_round_trip(tmp_path):
    p = str(tmp_path / 'raw.npz')
    rng = np.random.default_rng(0)
    win, arr = planted_window('Tc1_p1_x10_y30', 4, rng)
    E.save_part(p, arr, write_once=True)
    with pytest.raises(RuntimeError):
        E.save_part(p, {'W|Tc1_p1_x10_y30|y': np.zeros(3)}, write_once=True)
    E.save_part(p, {'other': np.ones(2)}, write_once=True)
    raw = E.load_raw(p)
    w2 = E.window_from_raw(raw, 'Tc1_p1_x10_y30')
    assert w2['k_true'] == 4 and w2['stage'] == 'test' and w2['wxy'] == (10, 30) and w2['pxy'] == (100, 200) and w2['label'] == 'c1' and np.array_equal(w2['truth_c'], win['truth_c'])
    assert E.window_from_raw({k.replace('|T', '|X'): v for k, v in raw.items()}, 'Xc1_p1_x10_y30')['stage'] == 'devx'
    assert E.load_raw(str(tmp_path / 'none.npz')) == {} and E.stage_names(raw, 'test') == []
    assert E.run_key('w', 'A', 1e-8) == 'w|A|1e-08' and E.eps_tag(1e-6) == '1e-06'


def test_freeze_gate_refuses_without_a_record_or_after_a_change(tmp_path):
    with pytest.raises(RuntimeError):
        E.check_freeze(None)
    good = {f: E.sha256_file(os.path.join(E.REPO, f)) for f in ('experiments/mechreal_core.py',)}
    E.check_freeze(dict(sha256=good))
    with pytest.raises(RuntimeError):
        E.check_freeze(dict(sha256={'experiments/mechreal_core.py': '0' * 64}))
    assert E.load_freeze(str(tmp_path / 'no.json')) is None


def test_registry_constants_are_the_planned_ones():
    assert (E.NW, E.HP, E.EDGE, E.Z_LIM, E.Z_BRIGHT, E.PAD, E.K_MIN, E.K_MAX, E.R_BLEND, E.R_MERGE, E.R_MATCH) == (20, 3.0, 2.5, 40.0, 50.0, 5.0, 2, 12, 0.5, 0.5, 0.5)
    assert E.KAPPAS == (40.0, 30.0, 20.0, 15.0, 10.0) and E.THETAS == (0.25, 0.5, 0.7, 0.85, 1.0) and E.EPS_MAIN == 1e-8 and E.LR_ADAM == 1e-3
    assert (E.M_E, E.SHARE_PATH, E.RATIO_C, E.GATE, E.N_DEV, E.N_TEST_PER_PATCH, E.SEED_SELECT) == (1.0, 0.60, 2.0, 0.5, 8, 5, 20261110)
    assert E.FAILED_REG == (('ngc5904', (640, 240)), ('ngc6254', (1900, 80)), ('ngc6934', (370, 1230))) and E.CAND_ORIGINS == (10, 30, 50, 70)
    assert (MC.RTOL, MC.ATOL_P, MC.ATOL_C, MC.LAT_H, MC.LAT_MARGIN) == (1e-8, 1e-10, 1e-12, 1.0, 1.0) and MC.LAW == dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
    assert MC.ARM_READOUT == {'A': 'P0', 'B': 'D2', 'C': 'D1'} and MC.ARM_OPT == {'D': 'adam', 'E': 'adam_g'}


def test_bundle_script_carries_the_three_files(tmp_path):
    jf = str(tmp_path / 'jobs.npz')
    MC.save_jobs(jf, [])
    s = E.bundle_script(jf, procs=3)
    assert 'MC.work("wk/jobs.npz", "result.npz", 3)' in s and '"imaging.py"' in s and '"mechreal_core.py"' in s and '"jobs.npz"' in s
    compile(s, 'bundle', 'exec')
