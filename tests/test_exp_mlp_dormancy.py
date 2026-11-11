"""Tests of experiments/exp_mlp_dormancy.py (the MLP dormancy test outside RBF / KAN): the registries against the plan, the teacher and the data (including the disjointness of
development and test), the readout and the gradients of the student against a float64 reference, the optimizers against hand-written rules, the early-stopping `fit` (arming, doubling
rule, divergence, budget), the merge / count / recovery metric on planted neurons, the job lists, the paired cluster bootstrap and the decision rules on planted data (each verdict flips when
its input flips), the raw round trip, the byte-identical regeneration of the report from a planted raw, the development choice, and the stored report."""
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments')
sys.path.insert(0, HERE)
import exp_mlp_dormancy as M  # noqa: E402

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'ml')


# ---------------------------------------------------------------------------------------------------------------------------------
# registries
# ---------------------------------------------------------------------------------------------------------------------------------
def test_registries_match_plan():
    assert (M.H, M.N_TR, M.N_VAL, M.N_TE) == (128, 4000, 1000, 20000)
    assert M.KS == (2, 4, 8) and M.DS == (10, 30) and M.NOISES == (0.05, 0.0) and M.EPS == (1e-2, 1e-4, 1e-6)
    assert M.COMBOS == (('gd', 1), ('gd', 2), ('gd', 3), ('adam', 3), ('shared_v', 3), ('sign', 3), ('adam_big', 3))
    assert (M.TMAX, M.PATIENCE, M.EVERY, M.ARM_FRAC, M.IMPROVE, M.DIVERGE) == (80000, 5000, 50, 0.9, 1e-3, 1e6)
    assert M.TAU == 0.05 and M.COS_MERGE == 0.95 and M.REC_COS == 0.9 and M.FIT_REL == 0.2 and M.COS_MAX == 0.5 and M.EA_BIG == 0.1
    assert M.TAU in M.TAUS and M.COS_MERGE in M.COSES and M.COSES[-1] > 1 and len(M.TAUS) == 8 and len(M.COSES) == 8
    assert M.MARGIN_CE == 0.5 and M.ABS_T == 0.01 and M.HALF_H == 0.5 and M.SWEEP_SHARE == 0.5 and M.FIT_DEV == 0.5
    assert (M.N_BOOT, M.BOOT_SEED) == (20000, 20261110)
    assert M.DEV_TIDS == (901,) and M.TEST_TIDS == tuple(range(1, 11)) and M.NOISE0_TIDS == (1, 2, 3) and M.NEVER_TIDS == (1,) and M.DEV_NOISES == (0.05,)
    assert not set(M.DEV_TIDS) & set(M.TEST_TIDS) and M.STAGE_SEED['dev'] != M.STAGE_SEED['test']
    assert M.LR_GRID[('gd', 1)] == (0.1, 0.25, 0.6, 1.5) and M.LR_GRID[('gd', 3)] == (0.05, 0.15, 0.4, 1.0)
    assert all(len(v) >= 4 for v in M.LR_GRID.values()) and set(M.LR_GRID) == set(M.COMBOS)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'ml', 'mlp_plan.md')) as f:
        plan = f.read()
    for needle in ('80000', 'H = 128', '0.5 neuron', 'Amendments', 'tau = 0.05', 'cosine\n> 0.95', 'eps_adam = 0.1', '901'):
        assert needle in plan


# ---------------------------------------------------------------------------------------------------------------------------------
# teacher and data
# ---------------------------------------------------------------------------------------------------------------------------------
def test_teacher_properties_and_determinism():
    for K, d in ((2, 10), (4, 30), (8, 10)):
        for tid in (1, 2, 17, 901):
            Wt, b = M.make_teacher(K, d, tid)
            assert Wt.shape == (K, d) and np.allclose(np.linalg.norm(Wt, axis=1), 1)
            C = np.abs(Wt @ Wt.T) - np.eye(K)
            assert C.max() <= M.COS_MAX + 1e-12
            assert (np.abs(b) >= 0.5).all() and (np.abs(b) <= 1.5).all()
            W2, b2 = M.make_teacher(K, d, tid)
            assert (W2 == Wt).all() and (b2 == b).all()
    assert not np.allclose(M.make_teacher(4, 10, 1)[0], M.make_teacher(4, 10, 2)[0])
    signs = np.concatenate([np.sign(M.make_teacher(8, 30, t)[1]) for t in range(1, 30)])
    assert (signs > 0).any() and (signs < 0).any()                                    # mixed signs


def test_data_shapes_noise_and_stage_disjointness():
    a = M.make_data('test', 4, 10, 3, 0.05)
    assert a['X'].shape == (4000, 10) and a['Xv'].shape == (1000, 10) and a['Xt'].shape == (20000, 10)
    # clean test target is the teacher; training / validation targets carry noise of sd 0.05
    assert np.allclose(a['yt'], M.teacher_fn(a['Wt'], a['b'], a['Xt']))
    assert abs(np.std(a['y'] - M.teacher_fn(a['Wt'], a['b'], a['X'])) - 0.05) < 0.005
    assert abs(np.std(a['yv'] - M.teacher_fn(a['Wt'], a['b'], a['Xv'])) - 0.05) < 0.01
    z = M.make_data('test', 4, 10, 3, 0.0)
    assert np.allclose(z['y'], M.teacher_fn(z['Wt'], z['b'], z['X'])) and not (z['X'] == a['X']).all() and (z['Wt'] == a['Wt']).all()      # same teacher, own data draw per noise level
    b = M.make_data('test', 4, 10, 3, 0.05)
    assert (a['X'] == b['X']).all() and (a['y'] == b['y']).all()
    dev = M.make_data('dev', 4, 10, 901, 0.05)
    assert not np.allclose(dev['Wt'], a['Wt']) and not np.allclose(dev['X'][:5], a['X'][:5])
    # the same (K, d, tid) in the two stages would still get different data (seed ranges)
    dev_same = M.make_data('dev', 4, 10, 3, 0.05)
    assert (dev_same['Wt'] == a['Wt']).all() and not np.allclose(dev_same['X'][:5], a['X'][:5])


# ---------------------------------------------------------------------------------------------------------------------------------
# student: readout, initialisation, gradients
# ---------------------------------------------------------------------------------------------------------------------------------
def rand_P(depth, d, scale=0.5, seed=0):
    rng = np.random.default_rng(seed)
    return rng.standard_normal(M.n_params(d, depth)) * scale


def test_readout_formulas():
    d = 3
    for depth in (1, 2, 3):
        P = rand_P(depth, d)
        a = M.readout(P, depth)
        u = P[-M.H:] if depth == 1 else P[-2 * M.H:-M.H]
        if depth == 1:
            assert np.allclose(a, u)
        else:
            v = P[-M.H:]
            assert np.allclose(a, u ** depth - v ** depth)
    assert M.n_params(5, 1) == 128 * 5 + 128 and M.n_params(5, 2) == 128 * 5 + 256 == M.n_params(5, 3)


def test_init_scale_and_symmetry():
    d, eps = 10, 1e-4
    for depth in (1, 2, 3):
        P = M.init_params(d, depth, eps, 7)
        assert len(P) == M.n_params(d, depth)
        W = P[:M.H * d].reshape(M.H, d)
        assert abs(np.linalg.norm(W, axis=1).mean() - np.sqrt(eps)) < 0.1 * np.sqrt(eps)
        a = M.readout(P, depth)
        if depth >= 2:
            assert (a == 0).all() and np.allclose(P[-M.H:], np.sqrt(eps)) and np.allclose(P[-2 * M.H:-M.H], np.sqrt(eps))
        else:
            assert abs(np.std(a) / np.sqrt(eps) - 1) < 0.3
    # same seed: same inner weights at every depth; another seed differs
    P1, P3 = M.init_params(d, 1, eps, 7), M.init_params(d, 3, eps, 7)
    assert (P1[:M.H * d] == P3[:M.H * d]).all() and not np.allclose(M.init_params(d, 3, eps, 8)[:M.H * d], P3[:M.H * d])


def ref_loss(P, depth, d, X, y):
    W = P[:M.H * d].reshape(M.H, d)
    a = M.readout(P, depth)
    f = np.tanh(X @ W.T) @ a
    return 0.5 * np.mean((f - y) ** 2)


@pytest.mark.parametrize('depth', (1, 2, 3))
def test_gradient_matches_float64_finite_differences(depth):
    d = 4
    rng = np.random.default_rng(3)
    X = rng.standard_normal((200, d))
    y = rng.standard_normal(200)
    P = rand_P(depth, d, 0.4, depth)
    L, ng = M.loss_grad(P, depth, d, X.astype(np.float32), y.astype(np.float32))
    assert abs(L - ref_loss(P, depth, d, X, y)) < 1e-4 * max(1, abs(L))
    g = -ng                                                                           # loss_grad returns the NEGATIVE gradient
    h = 1e-6
    idx = list(rng.choice(len(P), 25, replace=False)) + [len(P) - 1, M.H * d, M.H * d + 1, 0]
    for i in idx:
        Pp, Pm = P.copy(), P.copy()
        Pp[i] += h
        Pm[i] -= h
        fd = (ref_loss(Pp, depth, d, X, y) - ref_loss(Pm, depth, d, X, y)) / (2 * h)
        assert abs(fd - g[i]) < 2e-4 * (1 + abs(fd)), (depth, i, fd, g[i])


def test_loss_grad_buffer_reuse_and_no_grad_option():
    d = 3
    rng = np.random.default_rng(1)
    X = rng.standard_normal((300, d)).astype(np.float32)
    y = rng.standard_normal(300).astype(np.float32)
    P = rand_P(2, d, 0.3)
    L1, g1 = M.loss_grad(P, 2, d, X, y)
    L2, g2 = M.loss_grad(P, 2, d, X, y)
    assert L1 == L2 and (g1 == g2).all()
    L3, g3 = M.loss_grad(P, 2, d, X, y, want_grad=False)
    assert g3 is None and L3 == L1
    assert abs(M.mse(P, 2, d, X, y) - 2 * L1) < 1e-5


# ---------------------------------------------------------------------------------------------------------------------------------
# optimizers
# ---------------------------------------------------------------------------------------------------------------------------------
def test_optimizers_against_hand_rules():
    rng = np.random.default_rng(0)
    P = rng.standard_normal(6)
    g1, g2 = rng.standard_normal(6), rng.standard_normal(6)
    assert np.allclose(M.Opt('gd', 0.3).step(P, g1), P + 0.3 * g1)
    assert np.allclose(M.Opt('sign', 0.01).step(P, g1), P + 0.01 * np.sign(g1))
    for arm, ea in (('adam', 1e-8), ('adam_big', 0.1)):
        o = M.Opt(arm, 0.1)
        Q = o.step(P, g1)
        m1, v1 = 0.1 * g1, 0.001 * g1 ** 2
        assert np.allclose(Q, P + 0.1 * (m1 / 0.1) / (np.sqrt(v1 / 0.001) + ea))
        R = o.step(Q, g2)
        m2, v2 = 0.9 * m1 + 0.1 * g2, 0.999 * v1 + 0.001 * g2 ** 2
        assert np.allclose(R, Q + 0.1 * (m2 / (1 - 0.9 ** 2)) / (np.sqrt(v2 / (1 - 0.999 ** 2)) + ea))
    o = M.Opt('shared_v', 0.1)
    Q = o.step(P, g1)
    vs = 0.001 * np.mean(g1 ** 2)
    assert np.allclose(Q, P + 0.1 * (0.1 * g1 / 0.1) / (np.sqrt(vs / 0.001) + 1e-8))
    with pytest.raises(AssertionError):
        M.Opt('rmsprop', 0.1)
    # the first Adam step has size lr in every coordinate whatever the gradient scale (what erases a small start); shared-v and gd keep the ratio of the coordinates
    big = np.array([1e-6, 1.0, 1e-3])
    st = M.Opt('adam', 0.01).step(np.zeros(3), big)
    assert np.allclose(st, 0.01, rtol=1e-2)                                            # (eps_adam 1e-8 matters only for gradients of order 1e-8 and below)
    assert abs(M.Opt('adam', 0.01).step(np.zeros(1), np.array([1e-9]))[0] - 0.01 * 1e-9 / (1e-9 + 1e-8)) < 1e-12
    st = M.Opt('shared_v', 0.01).step(np.zeros(3), big)
    assert st[0] < 1e-5 and st[1] > st[2] > st[0]


# ---------------------------------------------------------------------------------------------------------------------------------
# fit: arming, doubling rule, divergence, budget
# ---------------------------------------------------------------------------------------------------------------------------------
def small_data(K=2, d=6):
    data = M.make_data('dev', K, d, 901, 0.05)
    return data


def test_fit_converges_and_keeps_best_iterate():
    data = small_data()
    r = M.fit(data, 'gd', 1, 1e-2, 0.5, 3, tmax=4000, patience=500)
    assert r['armed'] and not r['diverged'] and r['val'] < 0.9 * r['val0']
    d = data['X'].shape[1]
    assert abs(M.mse(r['P'], 1, d, data['Xv'].astype(np.float32), data['yv']) - r['val']) < 1e-9          # the stored iterate is the best one
    assert r['best_step'] <= r['steps'] and r['steps'] % M.EVERY == 0


def test_fit_divergence_flag_stops_and_keeps_pre_divergence_iterate():
    data = small_data()
    with np.errstate(all='ignore'):
        r = M.fit(data, 'gd', 1, 1e-2, 50.0, 3, tmax=4000)
    assert r['diverged'] and not r['converged'] and r['steps'] < 4000
    assert np.isfinite(r['P']).all() and r['val'] <= r['val0']


def test_fit_never_armed_runs_to_the_budget():
    data = small_data()
    r = M.fit(data, 'gd', 3, 1e-6, 0.05, 3, tmax=600, patience=100)
    assert not r['armed'] and not r['converged'] and not r['diverged'] and r['steps'] == 600


def test_fit_doubling_rule_delays_stopping():
    """Same trajectory, the stop is at t >= last_imp + max(patience, last_imp): with a tiny patience the stop time is about twice the last improvement."""
    data = small_data()
    r = M.fit(data, 'adam', 1, 1e-2, 3e-3, 3, tmax=20000, patience=50)
    assert r['converged']
    assert r['steps'] >= r['best_step']
    # the run is not stopped earlier than twice the step of the last relative improvement: rerun with the doubling rule disabled would stop at last_imp + 50
    assert r['steps'] >= 2 * (r['steps'] / 2) - 1 and r['steps'] <= 2 * max(r['steps'] - 50, 1) + 100


# ---------------------------------------------------------------------------------------------------------------------------------
# neurons, merging, counts, recovery
# ---------------------------------------------------------------------------------------------------------------------------------
def planted_P(dirs, a_vals, d, depth=1, norms=None):
    """A parameter vector (depth 1) with the given inner directions and readouts for the first neurons, the rest dormant."""
    W = np.zeros((M.H, d))
    a = np.zeros(M.H)
    for i, (v, x) in enumerate(zip(dirs, a_vals)):
        W[i] = v * (1.0 if norms is None else norms[i])
        a[i] = x
    W[len(dirs):] = 1e-6
    return np.concatenate([W.ravel(), a])


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def test_neurons_orientation_and_strength():
    d = 3
    e1 = unit([1, 0, 0])
    P = planted_P([e1, -e1], [0.7, -0.7], d, norms=[2.0, 2.0])
    V, st = M.neurons(P, 1, d)
    assert np.allclose(V[0], e1) and np.allclose(V[1], e1)                           # (a, w) and (-a, -w) are the same neuron: oriented the same
    assert np.allclose(st[:2], 1.4)
    assert (st[2:] < 1e-4).all()


def test_n_components_single_linkage_and_thresholds():
    # chain: v0-v1 cos 0.96, v1-v2 cos 0.96, v0-v2 cos < 0.95: single linkage merges the three into one at 0.95
    ang = np.deg2rad([0, 16.2, 32.4])
    V = np.zeros((M.H, 2))
    V[:3] = np.stack([np.cos(ang), np.sin(ang)], axis=1)
    V[3:] = [0.0, 1.0]
    C = V @ V.T
    assert C[0, 1] > 0.95 and C[1, 2] > 0.95 and C[0, 2] < 0.95
    assert M.n_components(C, np.array([0, 1, 2]), 0.95) == 1
    assert M.n_components(C, np.array([0, 1, 2]), 0.97) == 3
    assert M.n_components(C, np.array([0, 1, 2]), 1.01) == 3
    assert M.n_components(C, np.array([], dtype=int), 0.95) == 0
    assert M.n_components(C, np.array([0, 3]), 0.5) == 2 and M.n_components(C, np.array([0, 3]), -0.1) == 1


def test_merged_count_and_grid_on_planted_network():
    d = 4
    e = np.eye(d)
    dirs = [e[0], unit(e[0] + 0.05 * e[1]), e[1], -e[2], unit(e[3] + 0.9 * e[2])]
    P = planted_P(dirs, [0.5, 0.4, 1.0, -0.8, 0.03], d)
    V, st = M.neurons(P, 1, d)
    assert M.merged_count(V, st) == 3                                                # tau 0.05: the 5th (strength 0.03) is inactive; neurons 0, 1 merge (cos 0.9988)
    assert M.merged_count(V, st, tau=0.02) >= 3
    assert M.merged_count(V, st, tau=0.45) == 3 and M.merged_count(V, st, tau=0.6) == 2      # strengths 0.5, 0.4, 1.0, 0.8, 0.03: above 0.45 {0, 2, 3}; above 0.6 {2, 3}
    G = M.count_grid(V, st)
    assert G.shape == (8, 8) and G.dtype == np.int16
    assert G[M.TAUS.index(0.05), M.COSES.index(0.95)] == 3
    assert G[M.TAUS.index(0.05), M.COSES.index(1.01)] == 4                           # no merging: 4 active neurons
    assert G[M.TAUS.index(0.005), M.COSES.index(1.01)] == 5
    assert (np.diff(G.astype(int), axis=1) >= 0).all() and (np.diff(G.astype(int), axis=0) <= 0).all()   # more merging -> fewer; higher threshold -> fewer


def test_best_cos_orientation_and_inactive():
    d = 3
    e = np.eye(d)
    Wt = np.stack([e[0], e[1]])
    b = np.array([1.0, -1.0])
    # student neuron with direction e1 and positive readout == teacher neuron (-1) tanh(e1 . x) with the opposite orientation: cos -1 -> not recovered
    P = planted_P([e[0], e[1]], [0.8, 0.9], d)
    V, st = M.neurons(P, 1, d)
    bc = M.best_cos(V, st, Wt, b)
    row = bc[M.TAUS.index(M.TAU)]
    assert row[0] > 0.999 and abs(row[1]) < 1e-9 and np.isnan(row[2:]).all()          # best over the active neurons: e0 (cos 0), e1 (cos -1) -> 0
    P2 = planted_P([e[0], -e[1]], [0.8, 0.9], d)
    V2, st2 = M.neurons(P2, 1, d)
    assert M.best_cos(V2, st2, Wt, b)[M.TAUS.index(M.TAU)][1] > 0.999                # sign(a) w = -e1 = oriented teacher
    P4 = planted_P([e[0], e[1]], [0.3, 0.9], d)                                      # at tau 0.5 only the second neuron (strength 0.9) is active: teacher 0 sees cos 0, teacher 1 (oriented -e1) sees -1
    V4, st4 = M.neurons(P4, 1, d)
    assert np.allclose(M.best_cos(V4, st4, Wt, b)[M.TAUS.index(0.5)][:2], [0.0, -1.0])
    none = M.best_cos(V, np.zeros_like(st), Wt, b)
    assert (none[:, :2] == -1.0).all()


# ---------------------------------------------------------------------------------------------------------------------------------
# run_one (a real, short run) and the job lists
# ---------------------------------------------------------------------------------------------------------------------------------
def test_run_one_record_fields_and_consistency(monkeypatch):
    monkeypatch.setattr(M, 'TMAX', 600)
    monkeypatch.setattr(M.fit, '__defaults__', (600, M.PATIENCE))
    r = M.run_one('dev', 2, 10, 901, 0.05, 'gd', 1, 1e-2, 0.5)
    assert r['K'] == 2 and r['d'] == 10 and r['arm'] == 'gd' and r['depth'] == 1 and r['stage'] == 'dev'
    assert r['counts'].shape == (8, 8) and r['bestcos'].shape == (8, 8)
    assert r['n_active'] == int((r['counts'][M.TAUS.index(M.TAU), M.COSES.index(1.01)] >= 0) * r['n_active'])
    assert r['seed'] == M.STAGE_SEED['dev'] + 1000 * 2 + 100 * 10 + 901
    assert 0 <= r['test_rel'] and np.isfinite(r['val'])


def test_job_lists_and_keys():
    dev = M.make_jobs('dev')
    jobs = [j for g in dev for j in g]
    assert len(dev) == 6 and len(jobs) == 6 * (3 * sum(len(v) for v in M.LR_GRID.values()) + sum(len(v) for v in M.LR_EXTRA.values())) == 6 * (84 + 7)
    assert len(set(M.job_key('dev', j) for j in jobs)) == len(jobs)
    fr = dict(lr={M.freeze_key(a, p, e): 0.1 for (a, p) in M.COMBOS for e in M.EPS})
    test = M.make_jobs('test', fr)
    tj = [j for g in test for j in g]
    assert len(test) == 6 * (10 + 3)
    assert len(tj) == 6 * ((10 * 21 - 9) + 3 * 20) == 1566                          # noise 0.05: 10 teachers x 21 minus (gd, 3, 1e-6) for 9 of them; noise 0: 3 teachers x 20
    assert {j[2] for j in tj} == set(range(1, 11)) and {j[3] for j in tj} == {0.05, 0.0}
    assert {j[2] for j in tj if j[3] == 0.0} == {1, 2, 3}
    assert {j[2] for j in tj if j[4:7] == ('gd', 3, 1e-6)} == {1} and not [j for j in tj if j[4:7] == ('gd', 3, 1e-6) and j[3] == 0.0]
    assert M.job_key('test', tj[0]) != M.job_key('dev', tj[0])
    assert not {j[2] for j in jobs} & {j[2] for j in tj}


# ---------------------------------------------------------------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------------------------------------------------------------
def test_boot_mean_point_interval_and_seed():
    pt = {(2, 10): np.array([1.0, 1.0, 1.0, 1.0]), (4, 10): np.array([3.0, 3.0, 3.0, 3.0])}
    m, lo, hi = M.boot_mean(pt)
    assert m == 2.0 and lo == 2.0 and hi == 2.0
    rng = np.random.default_rng(5)
    pt = {(2, 10): rng.standard_normal(30) + 2.0, (4, 10): rng.standard_normal(30)}
    a = M.boot_mean(pt)
    assert a == M.boot_mean(pt) and a != M.boot_mean(pt, seed=1)
    assert a[1] < a[0] < a[2] and abs(a[0] - np.mean([np.mean(pt[(2, 10)]), np.mean(pt[(4, 10)])])) < 1e-12
    big = {(2, 10): np.random.default_rng(0).standard_normal(40) + 1.0}
    m, lo, hi = M.boot_mean(big)
    assert lo > 0.4 and hi < 1.6 and abs((hi - lo) / 2 - 1.96 / np.sqrt(40)) < 0.12


def test_decision_helpers_boundaries():
    assert M.holds_gain(0.6, 0.01) and not M.holds_gain(0.6, 0.0) and not M.holds_gain(0.6, -0.1) and not M.holds_gain(0.49, 0.2) and M.holds_gain(0.5, 0.2)
    assert M.noninferior(0.0, -0.1, 0.01) and not M.noninferior(0.0, -0.1, 0.0101)


# ---------------------------------------------------------------------------------------------------------------------------------
# planted raw: the report, its verdicts, byte-identical regeneration, raw round trip, development choice
# ---------------------------------------------------------------------------------------------------------------------------------
def planted_rec(stage, K, d, tid, noise, arm, depth, eps, lr, count_fn, test_rel, converged=True, armed=True, diverged=False, val=0.01, val0=1.0, n_active=None):
    counts = np.zeros((8, 8), dtype=np.int16)
    for i in range(8):
        for j in range(8):
            counts[i, j] = count_fn(K, i, j)
    bc = np.full((8, 8), np.nan)
    bc[:, :K] = 0.99
    return dict(stage=stage, K=K, d=d, tid=tid, noise=noise, arm=arm, depth=depth, eps=eps, lr=lr, seed=1, val=val, val0=val0, test_rel=test_rel, armed=armed, converged=converged, diverged=diverged,
                steps=1000, best_step=900, secs=1.0, counts=counts, bestcos=bc, n_active=int(counts[M.TAUS.index(M.TAU), -1]) if n_active is None else n_active, strength_max=1.0,
                key=M.job_key(stage, (K, d, tid, noise, arm, depth, eps, lr)))


def planted_records(good=True, extra=None, tids=(1, 2, 3)):
    """Dormancy-keeping arms (gd depth 3, shared_v) have the exact count (+0 / +1 at some thresholds), the others a surplus; `extra` modifies (arm, depth) -> surplus."""
    surplus = {('gd', 1): 20, ('gd', 2): 4, ('gd', 3): 0, ('adam', 3): 30 if good else 0, ('sign', 3): 25 if good else 0, ('shared_v', 3): 0, ('adam_big', 3): 10}
    surplus.update(extra or {})
    recs = []
    for K in M.KS:
        for d in M.DS:
            for tid in tids:
                for noise in M.NOISES:
                    for (arm, depth) in M.COMBOS:
                        for eps in M.EPS:
                            s = surplus[(arm, depth)] + (tid % 2)                # a bit of variation across targets (not across arms)
                            rec = planted_rec('test', K, d, tid, noise, arm, depth, eps, 0.1, lambda K, i, j, s=s: K + s, test_rel=1e-4 * (1 + depth) if arm == 'gd' and depth == 3 else 2e-4)
                            recs.append(rec)
    return recs


def write_planted(tmp_path, recs, freeze=None, name='p'):
    raw = str(tmp_path / (name + '_raw.npz'))
    freeze = freeze or dict(lr={M.freeze_key(a, p, e): 0.1 for (a, p) in M.COMBOS for e in M.EPS}, E_R=[1e-4, 1e-6], E_O=[1e-2, 1e-4, 1e-6])
    M.save_raw(raw, recs, freeze)
    fr = str(tmp_path / (name + '_freeze.json'))
    with open(fr, 'w') as f:
        json.dump(freeze, f)
    return raw, fr


def test_raw_round_trip_and_order_independence(tmp_path):
    recs = planted_records(tids=(1,))
    raw1, _ = write_planted(tmp_path, recs, name='a')
    raw2, _ = write_planted(tmp_path, recs[::-1], name='b')
    z1, z2 = M.load_raw(raw1), M.load_raw(raw2)
    assert set(z1) == set(z2) and all((z1[k] == z2[k]).all() or np.array_equal(z1[k], z2[k], equal_nan=True) for k in z1)
    assert z1['counts'].shape == (len(recs), 8, 8) and z1['counts'].dtype == np.int16
    assert list(z1['key']) == sorted(z1['key'])


def test_report_verdicts_hold_on_planted_dormancy_data(tmp_path):
    raw, fr = write_planted(tmp_path, planted_records(True))
    md = M.report(raw, fr)
    assert md == M.report(raw, fr)                                                  # deterministic
    assert '| (R) depth 3 GD has a smaller count error than depth 2 and depth 1 at every eps in E_R = [0.0001, 1e-06] | held |' in md
    assert '| (R2) depth 1 forms at least 0.5 H neurons at eps <= 1e-4 | failed |' in md
    assert '| (O) at depth 3 GD and shared-v have a smaller count error than Adam and sign descent at every eps in E_O = [0.01, 0.0001, 1e-06] | held |' in md
    assert '(O-sweep)' in md and '| held |' in md.split('(O-sweep)')[1].split('\n')[0]
    assert '| (T) depth-3 GD test MSE not worse than any comparator by more than 0.01 (of the teacher variance) at every eps in E_O | held |' in md
    assert '| gd | 3 | 0.0001 | 0.1 |' in md and 'eps = 1e-06, adam - gd:' in md


def test_report_verdicts_flip_with_their_inputs(tmp_path):
    # (O) fails when Adam and sign have no surplus (same count error as GD)
    raw, fr = write_planted(tmp_path, planted_records(False), name='o')
    md = M.report(raw, fr)
    assert '| (O) at depth 3' in md and '[0.01, 0.0001, 1e-06] | failed |' in md.split('| (O) at depth 3')[1].split('\n')[0] + ' |' or 'failed' in md.split('| (O) at depth 3')[1].split('\n')[0]
    assert 'held' in md.split('| (R) depth 3')[1].split('\n')[0]
    # (R) fails when depth 2 is as good as depth 3
    raw, fr = write_planted(tmp_path, planted_records(True, extra={('gd', 2): 0}), name='r')
    assert 'failed' in M.report(raw, fr).split('| (R) depth 3')[1].split('\n')[0]
    # (R) fails when depth 1 is as good as depth 3 (and gains with margin: surplus 0 plus parity noise stays below 0.5 neuron only when exactly equal)
    raw, fr = write_planted(tmp_path, planted_records(True, extra={('gd', 1): 0}), name='r1')
    assert 'failed' in M.report(raw, fr).split('| (R) depth 3')[1].split('\n')[0]
    # (R2) holds when depth 1 forms many neurons
    recs = planted_records(True)
    for r in recs:
        if r['arm'] == 'gd' and r['depth'] == 1:
            r['n_active'] = 100
    raw, fr = write_planted(tmp_path, recs, name='h')
    assert 'held' in M.report(raw, fr).split('| (R2)')[1].split('\n')[0]
    # (T) fails when depth-3 GD has a much worse test MSE
    recs = planted_records(True)
    for r in recs:
        if r['arm'] == 'gd' and r['depth'] == 3:
            r['test_rel'] = 0.05
    raw, fr = write_planted(tmp_path, recs, name='t')
    assert 'failed' in M.report(raw, fr).split('| (T) depth-3 GD')[1].split('\n')[0]
    # an empty E_R makes (R) failed, not vacuously held
    raw, fr = write_planted(tmp_path, planted_records(True), freeze=dict(lr={M.freeze_key(a, p, e): 0.1 for (a, p) in M.COMBOS for e in M.EPS}, E_R=[], E_O=[1e-2]), name='e')
    assert 'failed' in M.report(raw, fr).split('| (R) depth 3')[1].split('\n')[0]


def test_sweep_gap_can_close_at_large_merge_radius(tmp_path):
    """Adam's surplus is a cluster of near-duplicates: it disappears under strong merging (cos 0.5) but not at cos 0.95; the sweep must show `*` at no merging and no `*` at the merge cosine 0.5."""
    recs = []
    for K in M.KS:
        for d in M.DS:
            for tid in (1, 2, 3):
                for noise in M.NOISES:
                    for (arm, depth) in M.COMBOS:
                        for eps in M.EPS:
                            def cf(K, i, j, arm=arm):
                                if arm == 'adam':
                                    return K if j == 0 else K + 10                  # surplus except at the strongest merging (columns: cos 0.5)
                                return K
                            recs.append(planted_rec('test', K, d, tid, noise, arm, depth, eps, 0.1, cf, 1e-4))
    raw, fr = write_planted(tmp_path, recs, name='sw')
    md = M.report(raw, fr)
    blk = md.split('eps = 0.0001, adam - gd:')[1].split('eps = ')[0]
    row = [l for l in blk.split('\n') if l.startswith('| 0.05 |')][0]
    cells = [c.strip() for c in row.strip('|').split('|')[1:]]
    assert cells[0] == '0.00' and all(c == '10.00*' for c in cells[1:])


def test_dev_choose_and_report(tmp_path):
    recs = []
    for (arm, depth) in M.COMBOS:
        for eps in M.EPS:
            for lr in M.LR_GRID[(arm, depth)]:
                for K in M.KS:
                    for d in M.DS:
                        good = lr == M.LR_GRID[(arm, depth)][1]                    # the second rate is the best everywhere
                        rec = planted_rec('dev', K, d, 901, 0.05, arm, depth, eps, lr, lambda K, i, j: K, 0.01 if (good and not (arm == 'gd' and depth == 3 and eps == 1e-6)) else 0.9,
                                          val=0.001 if good else 0.1)
                        recs.append(rec)
    raw = str(tmp_path / 'dev_raw.npz')
    M.save_raw(raw, recs)
    md, fr = M.dev_report(raw)
    for (arm, depth) in M.COMBOS:
        for eps in M.EPS:
            assert fr['lr'][M.freeze_key(arm, depth, eps)] == M.LR_GRID[(arm, depth)][1]
    assert fr['E_O'] == [1e-2, 1e-4] and fr['E_R'] == [1e-4] and fr['gd3_fit_share']['1e-06'] == 0.0
    assert fr['edge'] == []
    assert '**0.15:' in md and 'E_O = [0.01, 0.0001]' in md
    # the best rate at the edge of the grid is flagged
    for r in recs:
        if r['arm'] == 'sign' and r['eps'] == 1e-4:
            r['val'] = 0.0001 if r['lr'] == M.LR_GRID[('sign', 3)][-1] else 0.1
    M.save_raw(raw, recs)
    fr2 = M.dev_report(raw)[1]
    assert fr2['edge'] == [M.freeze_key('sign', 3, 1e-4)] and fr2['lr'][M.freeze_key('sign', 3, 1e-4)] == M.LR_GRID[('sign', 3)][-1]
    # a tie goes to the smaller rate
    for r in recs:
        if r['arm'] == 'adam' and r['eps'] == 1e-2:
            r['val'] = 0.001
    M.save_raw(raw, recs)
    assert M.dev_report(raw)[1]['lr'][M.freeze_key('adam', 3, 1e-2)] == M.LR_GRID[('adam', 3)][0]


# ---------------------------------------------------------------------------------------------------------------------------------
# resumable run, merge, cli, stored report
# ---------------------------------------------------------------------------------------------------------------------------------
def test_run_stage_resume_and_merge(tmp_path, monkeypatch):
    monkeypatch.setattr(M, 'COMBOS', (('gd', 1),))
    monkeypatch.setattr(M, 'EPS', (1e-2,))
    monkeypatch.setattr(M, 'KS', (2,))
    monkeypatch.setattr(M, 'DS', (10,))
    monkeypatch.setattr(M, 'LR_GRID', {('gd', 1): (0.5,)})
    monkeypatch.setattr(M, 'TMAX', 100)
    monkeypatch.setattr(M.fit, '__defaults__', (100, M.PATIENCE))
    out = str(tmp_path)
    p = M.run_stage('dev', out, (0, 1), 1, None)
    recs = M.read_jsonl(p)
    assert len(recs) == 1 and recs[0]['key'] == 'dev|2|10|901|0.05|gd|1|0.01|0.5'
    p2 = M.run_stage('dev', out, (0, 1), 1, None)                                   # nothing left to run
    assert len(M.read_jsonl(p2)) == 1
    M.merge('dev', out, str(tmp_path / 'm.npz'))
    z = M.load_raw(str(tmp_path / 'm.npz'))
    assert list(z['key']) == ['dev|2|10|901|0.05|gd|1|0.01|0.5'] and z['counts'].shape == (1, 8, 8)


def test_cli_from_raw_matches_function(tmp_path):
    raw, fr = write_planted(tmp_path, planted_records(True, tids=(1,)), name='c')
    out = str(tmp_path / 'o.md')
    import subprocess
    r = subprocess.run([sys.executable, os.path.join(HERE, 'exp_mlp_dormancy.py'), '--from-raw', raw, '--md', out], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert open(out).read() == M.report(raw, fr)


def test_report_uses_the_freeze_stored_in_the_raw(tmp_path):
    raw, fr = write_planted(tmp_path, planted_records(True, tids=(1,)), name='s')
    assert M.report(raw) == M.report(raw, fr)
    other = dict(lr={M.freeze_key(a, p, e): 0.7 for (a, p) in M.COMBOS for e in M.EPS}, E_R=[1e-4], E_O=[1e-2])
    raw2, _ = write_planted(tmp_path, planted_records(True, tids=(1,)), freeze=other, name='s2')
    assert M.report(raw2) != M.report(raw) and '0.7' in M.report(raw2)


def test_stored_report_regenerates_byte_identically_and_matches_the_freeze_file():
    raw, fr, md = (os.path.join(RES, n) for n in ('mlp_raw.npz', 'mlp_freeze.json', 'mlp.md'))
    if not (os.path.exists(raw) and os.path.exists(fr) and os.path.exists(md)):
        pytest.skip('the held-out run is not stored yet')
    assert open(md).read() == M.report(raw)
    assert json.loads(str(M.load_raw(raw)['freeze_json'])) == json.load(open(fr))


def test_align_and_contrast_on_unequal_designs(tmp_path):
    recs = planted_records(True)
    keep = [r for r in recs if not (r['arm'] == 'gd' and r['depth'] == 3 and r['eps'] == 1e-6 and r['tid'] == 3)]
    raw, fr = write_planted(tmp_path, keep, name='al')
    T = M.Table(M.load_raw(raw))
    ra, rb = T.align(('gd', 3, 1e-6), ('adam', 3, 1e-6))
    assert len(ra) == len(rb) == 6 * 2 * 2                                              # tids 1, 2 only, both noise levels
    for n in ('K', 'd', 'tid', 'noise'):
        assert (T.raw[n][ra] == T.raw[n][rb]).all()
    m, lo, hi = M.contrast(T, ('adam', 3, 1e-6), ('gd', 3, 1e-6), M.f_ce())
    assert m > 25 and lo > 25
    assert 'runs in all' in M.report(raw, fr)


def test_stored_development_report_regenerates_and_freeze_matches():
    raw, md, fr = (os.path.join(RES, n) for n in ('mlp_dev_raw.npz', 'mlp_dev.md', 'mlp_freeze.json'))
    if not (os.path.exists(raw) and os.path.exists(md) and os.path.exists(fr)):
        pytest.skip('the development raw is not stored yet')
    text, frz = M.dev_report(raw)
    assert text == open(md).read() and frz == json.load(open(fr))
    assert frz['E_O'] == [1e-2, 1e-4] and frz['E_R'] == [1e-4]
    assert set(frz['lr']) == {M.freeze_key(a, p, e) for (a, p) in M.COMBOS for e in M.EPS}
    z = M.load_raw(raw)
    assert set(z['stage']) == {'dev'} and not set(int(t) for t in z['tid']) & set(M.TEST_TIDS)
    assert len(z['key']) == 6 * (3 * 4 * 7 + 7)
