"""Tests of experiments/exp_optdepth.py (the optimizer test P1 and the depth test P2): the registries against the plan and the frozen files, the freshness check (with a planted
overlap), the numpy optimizers against theory3f_escape's, the gradients of the model by finite differences, the early-stopping `fit`, the merge/exact-share metric on planted records,
the control's selection rules, the paired cluster bootstrap and the decision rules on planted data (each flips when its input flips), the closed-form unit (scalar integration, inversion,
exponent fit with exclusions), the formed-set agreement cell against `exp_parametrization.cell`, the byte-identical regeneration of the report from a raw, and the stored report."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments')
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import exp_optdepth as O  # noqa: E402
import exp_kan_spike as KS  # noqa: E402

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'ml')
SCRIPT = O.__file__


def tmp_out(tmp_path):
    return str(tmp_path)


# ---------------------------------------------------------------------------------------------------------------------------------
# registries
# ---------------------------------------------------------------------------------------------------------------------------------
def test_registries_match_plan_and_frozen_files():
    assert O.DEV_TIDS == (401, 402, 403) and O.DEV_SEEDS == (40, 41)
    assert O.TEST_TIDS == tuple(range(411, 431)) and O.TEST_SEEDS == (42, 43, 44) and O.NOISES == (0.1, 0.4)
    assert not (set(O.DEV_TIDS) | set(O.TEST_TIDS)) & set(O.USED_TIDS)
    assert not (set(O.DEV_SEEDS) | set(O.TEST_SEEDS)) & set(O.USED_SEEDS)
    assert O.MERGE_R == 0.3 and 0.3 in O.RADII and O.RADII[0] == 0.0
    assert (O.GATE_GAP, O.GATE_NEAR) == (0.20, 0.10)
    assert O.ARMS == ('gd', 'adam', 'adam_big', 'shared_v', 'sign')
    assert (O.M_UNITS, O.EPS0, O.STEPS, O.PATIENCE, O.EVERY, O.FT_STEPS) == (10, 1e-3, 12000, 2000, 25, 600)
    assert O.TOL_EXP == 0.05 and O.TOL_SCALAR == 0.01 and O.EPS_I == (1e-4, 1e-6, 1e-8, 1e-10) and O.DEPTHS_I == (3, 4, 5)
    assert O.EPS_II == (1e-6, 1e-8) and O.DEPTHS_II == (2, 3, 4, 5)
    assert O.RAND_NET == 0.5 and O.M_NET == 8 and O.MARGIN_PICK == 0.15
    with open(os.path.join(RES, 'kan_freeze.json')) as f:
        fr = json.load(f)
    o = fr['ours_synth']
    assert o['kind'] == 'gd' and o['lr'] == 1.0 and o['eps'] == O.EPS0 and o['nu_ratio'] == 0.05
    assert O.GD_CFG == dict(lr=o['lr']) and O.GD_RATES == dict(a=1.0, n=0.05, c=1.0, ls=1.0, b=0.1)
    assert O.CTRL_LR == fr['lr_synth']['ctrl'] and O.LAMS == KS.LAMS


def test_check_unused_passes_and_detects_a_planted_overlap(monkeypatch):
    lines, ok = O.check_unused()
    assert ok and len(lines) >= 5
    monkeypatch.setattr(O, 'USED_TIDS', O.USED_TIDS + (411,))
    assert O.check_unused()[1] is False
    monkeypatch.setattr(O, 'USED_TIDS', O.USED_TIDS[:-1])
    monkeypatch.setattr(O, 'USED_SEEDS', (0, 1, 2, 3, 42))
    assert O.check_unused()[1] is False


def test_earlier_ids_and_seeds_are_those_of_the_stored_raws():
    for raw, tids, seeds in (('kan_raw.npz', set(range(1, 11)), {0, 1, 2}), ('kan_dev_raw.npz', {101, 102, 103}, {0, 1})):
        _, recs = KS.load_raw(os.path.join(RES, raw))
        syn = [r for r in recs if r['kind'] == 'synth']
        assert {r['tid'] for r in syn} == tids and {r['seed'] for r in syn} == seeds
    assert set(O.USED_TIDS) >= set(range(1, 11)) | {101, 102, 103}


# ---------------------------------------------------------------------------------------------------------------------------------
# optimizers and model
# ---------------------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize('kind,kw', [('gd', {}), ('sign', {}), ('adam', {}), ('adam', dict(eps_adam=0.5)), ('adam_g', {}), ('adam_g', dict(pre=True))])
def test_opt_equals_theory3f(kind, kw):
    import theory3f_escape as T
    rng = np.random.default_rng(3)
    n = 12
    pre = rng.uniform(0.1, 2.0, n) if kw.get('pre') else None
    kw2 = {k: v for k, v in kw.items() if k != 'pre'}
    lr = rng.uniform(0.01, 0.1, n)
    a = O.Opt(kind, lr, pre=pre, **kw2)
    b = T.Opt(kind, lr, pre=pre, **kw2)
    pa = pb = rng.standard_normal(n)
    for _ in range(40):
        g = rng.standard_normal(n) * 10.0 ** rng.uniform(-3, 1, n)
        pa, pb = a.step(pa, g), b.step(pb, g)
        assert np.allclose(pa, pb, rtol=1e-12, atol=1e-14)


def test_opt_sign_adam_shared_properties():
    g = np.array([1e-6, -1.0, 3.0])
    o = O.Opt('sign', 0.1)
    assert np.allclose(o.step(np.zeros(3), g), [0.1, -0.1, 0.1])
    o = O.Opt('adam', 0.1, eps_adam=0.0)
    assert np.allclose(o.step(np.zeros(3), g), [0.1, -0.1, 0.1])          # Adam's first step is the sign step: the scale of the gradient is erased
    o = O.Opt('adam_g', 0.1, eps_adam=0.0)
    st = o.step(np.zeros(3), g)
    assert abs(st[0]) < 1e-6 and abs(st[2]) > abs(st[1])                  # one shared scale: the relative sizes of the gradients are kept
    o = O.Opt('adam', 0.1 * 1.0, eps_adam=1e6)
    assert np.allclose(o.step(np.zeros(3), g), 0.1 * g / 1e6, rtol=1e-6)  # a large eps_adam: a gradient-descent step
    with pytest.raises(AssertionError):
        O.Opt('lbfgs', 1.0)


def _fd_check(model, X, y, l1):
    def loss(z):
        old = {k: v.copy() for k, v in model.p.items()}
        model.setvec(z)
        r = y - model.fwd(X)[0]
        v = 0.5 * np.mean(r ** 2) + l1 * np.abs(model.weights()).sum()
        model.p = old
        return v
    z0 = model.vec(model.p).copy()
    g = model.vec(model.grad(X, y, l1))
    rng = np.random.default_rng(0)
    for k in rng.choice(z0.size, size=25, replace=False):
        h = 1e-6
        e = np.zeros_like(z0)
        e[k] = h
        fd = -(loss(z0 + e) - loss(z0 - e)) / (2 * h)
        assert abs(fd - g[k]) <= 1e-5 * max(1.0, abs(fd)), (k, fd, g[k])


def test_gradients_by_finite_differences_fact_and_free():
    rng = np.random.default_rng(1)
    X = rng.uniform(-2.5, 2.5, (60, 3))
    y = rng.standard_normal(60)
    m = O.AddNet(3, 5, 'fact', 0.05)
    for k in ('a', 'c', 'ls'):
        m.p[k] = m.p[k] + 0.3 * rng.standard_normal(m.p[k].shape)
    m.p['b'] = np.array([0.2])
    _fd_check(m, X, y, 0.0)
    f = O.AddNet(3, 5, 'free', seed=2)
    f.p['b'] = np.array([-0.1])
    f.mask[0, 1] = 0.0
    f.mask[2, 3] = 0.0
    _fd_check(f, X, y, 0.0)
    f.p['w'] = f.p['w'] + 0.2                       # all weights away from 0 (the sign is defined) for the L1 term
    _fd_check(f, X, y, 0.05)
    g = f.grad(X, y, 0.05)
    assert g['w'][0, 1] == 0.0 and g['c'][0, 1] == 0.0 and g['ls'][0, 1] == 0.0      # a pruned unit does not move


def test_model_start_is_the_dormant_small_start_and_matches_theory3f():
    import theory3f_escape as T
    m = O.AddNet(5, 10, 'fact', 1e-3)
    assert np.all(m.weights() == 0.0) and np.allclose(m.p['n'] ** 2, 1e-3) and np.all(m.p['a'] == 0)
    assert np.allclose(m.p['c'][0], np.linspace(-2.5, 2.5, 10)) and np.allclose(np.exp(m.p['ls']), 5.0 / 9.0)
    t = T.AddKan(5, 10, 1e-3)
    rng = np.random.default_rng(5)
    X = rng.uniform(-2.5, 2.5, (80, 5))
    y = rng.standard_normal(80)
    for k in ('a', 'c', 'ls'):
        d = 0.2 * rng.standard_normal(m.p[k].shape)
        m.p[k] = m.p[k] + d
        t.p[k] = t.p[k] + d
    ga, gb = m.grad(X, y), t.grad(X, y)
    assert all(np.allclose(ga[k], gb[k], rtol=1e-12, atol=1e-15) for k in gb)
    assert np.allclose(m.fwd(X)[0], t.fwd(X)[0]) and np.array_equal(m.counts(X), t.counts(X))


def small_data(seed=0, n=300):
    rng = np.random.default_rng(seed)
    X = rng.uniform(-2.5, 2.5, (n + 100 + 50, 2))
    f = 1.2 * np.exp(-0.5 * ((X[:, 0] - 0.8) / 0.4) ** 2)
    y = f + 0.05 * rng.standard_normal(f.shape)
    return dict(Xtr=X[:n], ytr=y[:n], Xva=X[n:n + 100], yva=y[n:n + 100], Xte=X[n + 100:], yte=f[n + 100:])


def test_fit_matches_theory3f_fit_add_and_keeps_best_iterate():
    import theory3f_escape as T
    d = small_data()
    ma = O.AddNet(2, 6, 'fact', 1e-2)
    mb = T.AddKan(2, 6, 1e-2)
    rv = ma.rates_vec(O.GD_RATES)
    sa = O.fit(ma, d, O.Opt('gd', rv), 400)
    sb = T.fit_add(mb, d, T.Opt('gd', rv), 400)
    assert sa[0] == sb and all(np.array_equal(ma.p[k], mb.p[k]) for k in mb.p)
    # the returned model is the best-validation iterate: its validation MSE equals the returned best
    assert abs(np.mean((ma.fwd(d['Xva'])[0] - d['yva']) ** 2) - sa[1]) < 1e-15


def test_fit_stops_on_patience_and_on_nonfinite():
    d = small_data()
    m = O.AddNet(2, 6, 'free', seed=0)
    st = O.fit(m, d, O.Opt('adam', 0.03), 12000, patience=100)
    assert st[0] < 12000 and st[0] % O.EVERY == 0
    m2 = O.AddNet(2, 6, 'free', seed=0)
    st2 = O.fit(m2, d, O.Opt('gd', 1e6), 500)                       # diverges: stops at the first non-finite validation MSE, returns the best finite iterate
    assert np.isfinite(st2[1]) and st2[0] < 500


def test_gd_arm_leaves_dormant_start_and_units_are_few(monkeypatch=None):
    d = small_data()
    m = O.AddNet(2, 8, 'fact', 1e-3)
    O.fit(m, d, O.make_opt('gd', O.GD_CFG, m), 3000)
    c = m.counts(d['Xtr'])
    assert c[0] >= 1 and c[1] == 0                                  # the relevant edge forms, the irrelevant edge stays dormant


def test_make_opt_arms():
    m = O.AddNet(2, 4, 'fact', 1e-3)
    rv = m.rates_vec(O.GD_RATES)
    assert rv.size == m.vec(m.p).size
    assert O.make_opt('gd', dict(lr=1.0), m).kind == 'gd' and np.array_equal(O.make_opt('gd', dict(lr=2.0), m).lr, 2.0 * rv)
    a = O.make_opt('adam', dict(lr=0.01), m)
    assert a.kind == 'adam' and a.lr == 0.01 and a.ea == 1e-8
    b = O.make_opt('adam_big', dict(ea=0.1), m)
    assert b.kind == 'adam' and b.ea == 0.1 and np.allclose(b.lr, 0.1 * rv)
    s = O.make_opt('shared_v', dict(lr=0.003, pre=True), m)
    assert s.kind == 'adam_g' and np.array_equal(s.pre, rv)
    assert O.make_opt('shared_v', dict(lr=0.003, pre=False), m).pre is None
    assert O.make_opt('sign', dict(lr=0.001), m).kind == 'sign'
    with pytest.raises(ValueError):
        O.make_opt('nope', {}, m)


# ---------------------------------------------------------------------------------------------------------------------------------
# metrics on planted records
# ---------------------------------------------------------------------------------------------------------------------------------
def rec(arm, tid, noise, seed, centres, truth, conf=None, val=1.0, P=None, steps=100, rmse=0.1, n_val=500):
    return dict(kind='p1', tid=tid, noise=noise, seed=seed, arm=arm, conf=conf or json.dumps(dict(lr=1.0), sort_keys=True), steps=steps, val_mse=val, test_rmse=rmse, n_val=n_val,
                counts=[len(c) for c in centres], centres=centres, P=P if P is not None else 3 * sum(len(c) for c in centres) + 1, truth=truth, truth_c=[[0.0] * k for k in truth], secs=1.0)


def test_merge_exact_share_planted():
    truth = [0, 1, 2, 3]
    r = rec('gd', 1, 0.1, 0, [[], [0.0, 0.2], [-1.0, 1.0], [-2.0, 0.0, 0.1, 2.0]], truth)
    assert list(O.merged_counts(r, 0.0)) == [0, 2, 2, 4]
    assert list(O.merged_counts(r, 0.3)) == [0, 1, 2, 3]
    assert O.share_pos(r, 0.0) == pytest.approx(1 / 3) and O.share_pos(r, 0.3) == 1.0
    assert list(O.exact_vec(r, 0.3)) == [True] * 4
    r0 = rec('gd', 1, 0.1, 0, [[], []], [0, 0])
    assert np.isnan(O.share_pos(r0, 0.3))
    # a pair exactly one radius apart is not merged (single linkage, strictly closer)
    assert O.merged_counts(rec('gd', 1, 0.1, 0, [[0.0, 0.3]], [1]), 0.3)[0] == 2


def test_pick_ctrl_one_se_and_oracle():
    truth = [2, 1]
    lam = lambda x: json.dumps(dict(lam=x))
    recs = [rec('ctrl', 1, 0.1, 0, [[0.0, 0.1, 1.0], [0.0, 0.1]], truth, lam(0.0), val=1.00, P=19),
            rec('ctrl', 1, 0.1, 0, [[0.0, 1.0], [0.5]], truth, lam(0.01), val=1.01, P=10),      # exact, fewer parameters, within the one-SE band
            rec('ctrl', 1, 0.1, 0, [[0.0, 1.0], [0.5]], truth, lam(0.03), val=1.60, P=10)]
    assert O.pick_ctrl(recs)['conf'] == lam(0.01)                                                # the fewest parameters inside the band
    assert O.pick_ctrl(recs, 0.3)['conf'] == lam(0.01)                                           # oracle: highest merged exact share, then fewer parameters, then the smaller penalty
    tie = [rec('ctrl', 1, 0.1, 0, [[0.0, 1.0], [0.5]], truth, lam(0.003), val=1.0, P=10), rec('ctrl', 1, 0.1, 0, [[0.0, 1.0], [0.5]], truth, lam(0.01), val=1.0, P=10)]
    assert O.pick_ctrl(tie, 0.3)['conf'] == lam(0.003)


def fake_runs(share_by_arm, tids=range(1, 6), noises=(0.1, 0.4), seeds=(0, 1, 2), delta=None):
    """Records for every arm; the arm's K>0 exact share on target t is share_by_arm[arm][t - 1] (realised by making that fraction of the 6 runs exact on the single K>0 edge)."""
    truth = [0, 2]
    recs = []
    for t in tids:
        runs = [(n, s) for n in noises for s in seeds]
        for arm in O.ARMS:
            sh = share_by_arm[arm][t - 1]
            nex = int(round(sh * len(runs)))
            for i, (n, s) in enumerate(runs):
                good = i < nex
                recs.append(rec(arm, t, n, s, [[], [0.0, 1.0] if good else [0.0, 1.0, 2.0]], truth))
        for i, (n, s) in enumerate(runs):
            for lam in (0.0, 0.01):
                recs.append(rec('ctrl', t, n, s, [[], [0.0, 1.0] if lam else [0.0, 0.1, 1.0]], truth, json.dumps(dict(lam=lam)), val=1.0 + lam, P=7 if lam else 10))
    return recs


def test_arm_table_and_contrast_and_verdict_on_planted_data():
    hi, lo = [1.0] * 5, [0.0] * 5
    base = dict(gd=hi, shared_v=hi, adam_big=hi, adam=lo, sign=lo)
    runs = O.complete_runs(O.group_runs(fake_runs(base)), list(O.ARMS) + ['ctrl'])
    assert len(runs) == 5 * 6
    freeze = {a: dict(lr=1.0) for a in O.ARMS}
    t = O.arm_table(runs, 'gd', 0.3, freeze)
    assert t['share_pos'] == 1.0 and t['share_all'] == 1.0 and t['share_empty'] == 1.0 and t['surplus_pos'] == 0.0 and t['all5'] == 1.0
    t = O.arm_table(runs, 'adam', 0.3, freeze)
    assert t['share_pos'] == 0.0 and t['surplus_pos'] == 1.0 and t['share_empty'] == 1.0 and t['missing_pos'] == 0.0
    m, l, h, d = O.contrast(runs, 'gd', 'adam', 0.3, freeze)
    assert m == 1.0 and l == 1.0 and h == 1.0 and d == [1.0] * 5
    V = O.p1_verdict(runs, 0.3, freeze)
    assert V['a_ok'] and V['b_ok'] and V['d']['ok'] and len(V['a']) == 4
    # merging at 0.3 makes no difference here (the surplus unit is at distance 1.0); a merged duplicate would
    # flip (a): lower the gap between GD and Adam
    mid = dict(base, adam=[5 / 6] * 5, sign=[5 / 6] * 5)                    # gap 1/6 < 0.20
    V = O.p1_verdict(O.complete_runs(O.group_runs(fake_runs(mid)), list(O.ARMS) + ['ctrl']), 0.3, freeze)
    assert V['a_ok'] is False and all(abs(x['mean'] - 1 / 6) < 1e-12 for x in V['a']) and V['a'][0]['lo_ci'] > 0
    edge = dict(base, adam=[4 / 6] * 5, sign=[4 / 6] * 5)                   # gap 1/3 >= 0.20, interval above 0
    assert O.p1_verdict(O.complete_runs(O.group_runs(fake_runs(edge)), list(O.ARMS) + ['ctrl']), 0.3, freeze)['a_ok'] is True
    # (b) flips with the GD - shared-v gap
    far = dict(base, shared_v=[0.5] * 5)
    V = O.p1_verdict(O.complete_runs(O.group_runs(fake_runs(far)), list(O.ARMS) + ['ctrl']), 0.3, freeze)
    assert V['b_ok'] is False and V['b']['mean'] == pytest.approx(0.5)
    # (a) needs the interval above 0 too: a gap of 0.5 on a single noisy target out of five is not enough when the others are 0
    one = dict(base, gd=[1.0, 0.0, 0.0, 0.0, 0.0], shared_v=[1.0, 0.0, 0.0, 0.0, 0.0])
    V = O.p1_verdict(O.complete_runs(O.group_runs(fake_runs(one)), list(O.ARMS) + ['ctrl']), 0.3, freeze)
    assert V['a'][0]['mean'] == pytest.approx(0.2) and V['a'][0]['lo_ci'] <= 0.0 and V['a_ok'] is False


def test_merge_changes_arm_table_in_planted_duplicates():
    truth = [0, 1]
    recs = []
    for arm in O.ARMS:
        recs.append(rec(arm, 1, 0.1, 0, [[], [0.0, 0.2]] if arm == 'adam' else [[], [0.0]], truth))
    recs += [rec('ctrl', 1, 0.1, 0, [[], [0.0, 0.2]], truth, json.dumps(dict(lam=0.0)))]
    runs = O.complete_runs(O.group_runs(recs), list(O.ARMS) + ['ctrl'])
    freeze = {a: dict(lr=1.0) for a in O.ARMS}
    assert O.arm_table(runs, 'adam', 0.0, freeze)['share_pos'] == 0.0 and O.arm_table(runs, 'adam', 0.3, freeze)['share_pos'] == 1.0
    assert O.arm_table(runs, 'ctrl', 0.3, freeze)['share_pos'] == 1.0
    assert O.arm_table(runs, 'ctrl_or', 0.3, freeze)['share_pos'] == 1.0


def test_freeze_from_dev_picks_best_merged_share_with_unmerged_tiebreak():
    truth = [0, 1]
    recs = []
    cfgs = O.DEV_GRID
    for arm in O.ARMS[1:]:
        for i, cfg in enumerate(cfgs[arm]):
            for t, seed in ((401, 40), (402, 41)):
                if i == 1:
                    cs = [[], [0.0]]                       # the second configuration is exact
                elif i == 0:
                    cs = [[], [0.0, 0.1]]                  # duplicate: exact only after merging
                else:
                    cs = [[], [0.0, 1.0]]                  # over-count
                recs.append(rec(arm, t, 0.1, seed, cs, truth, json.dumps(cfg, sort_keys=True)))
    frz, table = O.freeze_from_dev(recs)
    assert frz['gd'] == O.GD_CFG
    for arm in O.ARMS[1:]:
        assert frz[arm] == cfgs[arm][1]
    # tie in the merged share: the unmerged share decides
    recs2 = [rec('adam', 401, 0.1, 40, [[], [0.0]] if i == 0 else [[], [0.0, 0.1]], truth, json.dumps(c, sort_keys=True)) for i, c in enumerate(cfgs['adam'])]
    recs2[1] = rec('adam', 401, 0.1, 40, [[], [0.0, 0.1]], truth, json.dumps(cfgs['adam'][1], sort_keys=True))       # (the same merged share as the first, a smaller unmerged one)
    for arm in ('adam_big', 'shared_v', 'sign'):
        for cfg in cfgs[arm]:
            recs2.append(rec(arm, 401, 0.1, 40, [[], [0.0]], truth, json.dumps(cfg, sort_keys=True)))
    frz2, _ = O.freeze_from_dev(recs2)
    assert frz2['adam'] == cfgs['adam'][0]


# ---------------------------------------------------------------------------------------------------------------------------------
# P2
# ---------------------------------------------------------------------------------------------------------------------------------
def test_closed_form_wmax_and_margin_factor():
    assert O.closed_form_wmax(3, 1e-6, 0.5) == pytest.approx((1e-6) ** 1.5 * (0.5 ** -3 - 1.5 ** -3))
    assert O.closed_form_wmax(4, 1e-6, 0.25) == pytest.approx((1e-6) ** 2 * (0.25 ** -2 - 1.75 ** -2))
    assert O.closed_form_wmax(2, 1e-6, 0.5) == pytest.approx(1e-6 ** 0.5 - 1e-6 ** 1.5, rel=1e-9)
    assert O.closed_form_wmax(1, 1e-6, 0.5) == O.PZ_LRC
    assert O.margin_factor(3, 1e-6, 0.5) == pytest.approx(0.5 ** -3 - 1.5 ** -3)
    assert O.closed_form_wmax(5, 1e-8, 0.3) == pytest.approx((1e-8) ** 2.5 * O.margin_factor(5, 1e-8, 0.3))
    # the upper bound of the plan: (eps b)^(p/2) mu^(-p/(p-2)) is not attained
    for p in (3, 4, 5):
        assert O.closed_form_wmax(p, 1e-6, 0.4) < (1e-6) ** (p / 2) * 0.4 ** (-p / (p - 2.0))


@pytest.mark.parametrize('p', [2, 3, 4, 5])
def test_x_of_w_inverts_the_closed_form(p):
    for s2 in (1e-4, 1e-10):
        for mu in (0.05, 0.3, 0.9):
            w = O.closed_form_wmax(p, s2, mu)
            assert O.x_of_w(p, w, s2) == pytest.approx(1 - mu, abs=1e-6)
    assert O.x_of_w(3, 0.0, 1e-6) == 0.0


def test_slope_helper():
    e = np.array([1e-4, 1e-6, 1e-8, 1e-10])
    assert O.slope(e, e ** 1.5 * 7) == pytest.approx(1.5)
    assert O.slope(e, 3.0 * np.ones(4)) == pytest.approx(0.0, abs=1e-12)


def planted_cells(p, margin_by_eps, b=1.0, n=6, form_unit_at=None, status=None, wrong_thr=0.009):
    """Cells (kind 'p2n') of ONE target whose largest dormant coefficient is exactly the closed-form sup at the given margins; unit 2 carries it."""
    cells = []
    for k, e in enumerate(O.EPS_I):
        w = np.zeros(n)
        w[2] = O.closed_form_wmax(p, e * b, margin_by_eps[k])
        w[0] = 0.5                                      # a formed unit
        wend = w.copy()
        if form_unit_at is not None and form_unit_at == k:
            wend[3] = 0.02
            w[3] = 0.02
        cells.append(dict(kind='p2n', eps=e, status=(status or {}).get(k, 'ok'), wmax=w.tolist(), wend=wend.tolist(), b=[b] * n, thr=wrong_thr))
    return cells


def test_fit_exponents_exact_on_planted_cells_and_exclusions():
    for p in (3, 4, 5):
        out = O.fit_exponents(planted_cells(p, [0.3] * 4), p, 'un')
        assert out['excluded'] is None and out['slope'] == pytest.approx(p / 2.0, abs=1e-9) and out['slope_corr'] == pytest.approx(p / 2.0, abs=1e-6)
        assert all(m == pytest.approx(0.3, abs=1e-6) for m in out['margin'])
    # a drifting margin shifts the raw exponent but not the exponent at the realised margin
    out = O.fit_exponents(planted_cells(3, [0.3, 0.3, 0.3, 0.1]), 3, 'un')
    assert abs(out['slope'] - 1.5) > 0.05 and out['slope_corr'] == pytest.approx(1.5, abs=1e-6)
    # margin below the floor excludes the target (but the raw exponent is still returned for the sensitivity table)
    out = O.fit_exponents(planted_cells(3, [0.3, 0.3, 0.3, 0.04]), 3, 'un')
    assert out['excluded'] is not None and np.isfinite(out['slope'])
    assert O.fit_exponents(planted_cells(3, [0.3, 0.3, 0.3, 0.04]), 3, 'un', margin_min=0.0)['excluded'] is None
    # a failed cell is dropped; two failed cells leave too few
    out = O.fit_exponents(planted_cells(5, [0.3] * 4, status={3: 'failed: Requ'}), 5, 'un')
    assert out['n_valid'] == 3 and out['slope'] == pytest.approx(2.5, abs=1e-9)
    out = O.fit_exponents(planted_cells(5, [0.3] * 4, status={2: 'failed', 3: 'failed'}), 5, 'un')
    assert out['excluded'] and np.isnan(out['slope'])
    # a unit formed at some eps is not in the dormant population N*: here unit 3 formed at the largest eps only and is excluded at all eps
    cells = planted_cells(3, [0.3] * 4, form_unit_at=0)
    out = O.fit_exponents(cells, 3, 'un')
    assert out['slope'] == pytest.approx(1.5, abs=1e-9) and set(out['idx']) == {2}
    # no dormant coefficient at all: excluded, not a crash
    cells = planted_cells(3, [0.3] * 4)
    for c in cells:
        c['wmax'] = [0.5, 0, 0, 0, 0, 0]
    assert O.fit_exponents(cells, 3, 'un')['excluded'] == 'zero coefficient'
    # everything formed: no dormant unit
    cells = planted_cells(3, [0.3] * 4)
    for c in cells:
        c['wend'] = [1.0] * 6
    assert O.fit_exponents(cells, 3, 'un')['excluded'] == 'no dormant unit'


def test_scalar_unit_exponents_are_p_over_2():
    for p, mu in ((3, 0.5), (4, 0.25), (5, 0.5)):
        rr = [O.scalar_unit(p, e, mu) for e in O.EPS_I]
        assert all(r['status'] == 'ok' for r in rr)
        for r, e in zip(rr, O.EPS_I):
            assert r['wmax'] == pytest.approx(O.closed_form_wmax(p, e, mu), rel=1e-7)
        assert O.slope(O.EPS_I, [r['wmax'] for r in rr]) == pytest.approx(p / 2.0, abs=1e-6)
    rr = [O.scalar_unit(2, e, 0.5) for e in O.EPS_I]
    assert O.slope(O.EPS_I, [r['wmax'] for r in rr]) == pytest.approx(0.5, abs=1e-3)
    rr = [O.scalar_unit(1, e, 0.0) for e in O.EPS_I]
    assert O.slope(O.EPS_I, [r['wmax'] for r in rr]) == pytest.approx(0.0, abs=1e-9) and rr[0]['wmax'] == pytest.approx(O.PZ_LRC)


def test_scalar_unit_reaches_the_margin_three_times_with_both_signs():
    # the drive h(t) = sin(3 pi t) reaches |h| = 1 at 1/6, 1/2, 5/6
    assert O.scalar_h(1 / 6) == pytest.approx(1.0) and O.scalar_h(0.5) == pytest.approx(-1.0) and O.scalar_h(5 / 6) == pytest.approx(1.0)


def test_mu_theta_values_and_ordering():
    assert O.mu_theta(3, 1e-6, 0.01) == pytest.approx(0.01 ** (-1 / 3) * 1e-3)
    assert O.mu_theta(4, 1e-6, 0.01) == pytest.approx(0.01 ** (-0.5) * 1e-6)
    assert O.mu_theta(5, 1e-8, 0.01) == pytest.approx(0.01 ** (-0.6) * 1e-12)
    assert O.mu_theta(2, 1e-6, 0.01) == pytest.approx(np.log(100) / np.log(1e6))
    for e in (1e-6, 1e-8):
        assert O.mu_theta(2, e, 0.01) > O.mu_theta(3, e, 0.01) > O.mu_theta(4, e, 0.01) > O.mu_theta(5, e, 0.01)
    assert O.mu_theta(3, 1e-8, 0.01) < O.mu_theta(3, 1e-6, 0.01)


def test_pick_i_targets_rule():
    refs = [dict(margin_never=m) for m in (0.5, 0.12, float('nan'), 0.2, 0.9, 0.15, 0.3, 0.16, 0.4, 0.6, 0.25, 0.33, 0.7)]
    pick = O.pick_i_targets([None] * len(refs), refs)
    cand = sorted((r['margin_never'], i) for i, r in enumerate(refs) if np.isfinite(r['margin_never']) and r['margin_never'] >= 0.15)
    assert len(pick) == O.M_NET and pick[0] == cand[0][1] and pick[-1] == cand[-1][1]
    ms = [refs[i]['margin_never'] for i in pick]
    assert ms == sorted(ms) and min(ms) >= 0.15
    few = O.pick_i_targets([None] * 3, [dict(margin_never=0.5), dict(margin_never=0.2), dict(margin_never=0.01)])
    assert few == [1, 0]


def test_dev_target_and_dev_pool_are_not_in_the_sample():
    cs, As = O.dev_target(7001)
    assert cs.shape == (3,) and np.min(np.diff(cs)) >= 0.35 and np.all((0.5 <= np.abs(As)) & (np.abs(As) <= 1.0))
    assert np.array_equal(O.dev_target(7001)[0], cs) and not np.array_equal(O.dev_target(7002)[0], cs)
    smp = O.p2_sample()
    pool = O.p2_dev_pool(5)
    keys = {(tuple(np.round(t['cs'], 9)), tuple(np.round(t['As'], 9))) for t in smp}
    assert len(smp) == 200 and len(pool) == 5
    assert all((tuple(np.round(t['cs'], 9)), tuple(np.round(t['As'], 9))) not in keys for t in pool)


def test_net_cell_design_target_closed_form_and_masses():
    import exp_parametrization as PZ
    ref, tau_end = O.net_target_ref(list(PZ.CSTAR), list(PZ.ASTAR))
    assert ref['set'] and tau_end > 1.0
    c = O.net_cell(list(PZ.CSTAR), list(PZ.ASTAR), 3, 1e-4, 0.0, tau_end)
    assert c['status'] == 'ok' and len(c['wmax']) == 21 and np.allclose(c['b'], 1.0)
    wend = np.array(c['wend'])
    assert (wend > c['thr']).sum() == len(ref['set'])               # the formed set is the path's
    cu = O.net_cell(list(PZ.CSTAR), list(PZ.ASTAR), 3, 1e-4, O.RAND_NET, tau_end)
    b = np.array(cu['b'])
    assert cu['status'] == 'ok' and not np.allclose(b, 1.0) and np.array_equal(b, O.net_cell(list(PZ.CSTAR), list(PZ.ASTAR), 3, 1e-6, O.RAND_NET, tau_end)['b'])
    # the largest dormant coefficient respects the closed-form sup at its own realised margin (identity of the unit)
    out = O.fit_exponents([dict(c, eps=1e-4)], 3, 'eq')
    assert out['n_valid'] == 1 and out['excluded'] == 'fewer than 3 valid eps'


def test_p2ii_cell_equals_exp_parametrization_cell():
    import exp_parametrization as PZ
    L, K, c0 = PZ.setup()
    t = O.p2_sample()[0]
    f = K.synth(t['cs'][:, None], t['As'])
    ref = PZ.reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0
    for p in (3, 4):
        a = O.p2ii_cell(K, f, c0, t['cs'], t['As'], ref, tau_end, p, 1e-6)
        b = PZ.cell(PZ.make_readout('D%d' % p), K, f, c0, t['cs'], t['As'], 1e-6, ref, tau_end)
        assert a['set_equal'] == bool(b['set_equal']) and a['order_equal'] == bool(b['order_equal']) and a['n_formed'] == b['n_formed'] and a['p'] == p


# ---------------------------------------------------------------------------------------------------------------------------------
# report, drivers
# ---------------------------------------------------------------------------------------------------------------------------------
def fake_p2_records():
    recs = []
    for p in (1, 2) + O.DEPTHS_I:
        for mu in ((0.25, 0.5) if p > 1 else (0.0,)):
            for e in O.EPS_I:
                recs.append(dict(kind='p2s', arm=repr(mu), p=p, eps=e, mu=mu, wmax=O.closed_form_wmax(p, e, mu), status='ok', exact=O.closed_form_wmax(p, e, mu)))
    for ti in range(5):
        recs.append(dict(kind='p2t', tgt=ti, name='t%d' % ti, margin_never=0.2 + 0.05 * ti, margin=0.2, order=[1, 2], tau_end=3.0, cs=[0.0], As=[1.0], sample_index=None))
        for arm in ('un', 'eq'):
            for p in (1, 2) + O.DEPTHS_I:
                cells = planted_cells(p if p > 1 else 3, [0.3] * 4, wrong_thr=0.5 if p <= 2 else 0.009)
                for c in cells:
                    c.update(tgt=ti, arm=arm, p=p)
                    if p == 5 and c['eps'] == 1e-10:
                        c.update(status='failed: Requ')
                    recs.append(c)
    n = 10
    for i in range(n):
        recs.append(dict(kind='p2r', tgt=i, fam='A', cs=[0.0], As=[1.0], margin=0.3, margin_never=0.1 + 0.05 * i, thr=0.01, order=[0], set=[0], relaxed=True))
        for p in O.DEPTHS_II:
            for e in O.EPS_II_ALL:
                ok = (p >= 3 and i < 5 + p) or (p == 2 and i < 1)
                recs.append(dict(kind='p2c', tgt=i, p=p, eps=e, set_equal=bool(ok), order_equal=bool(ok), n_formed=1, exact=bool(ok), status='ok', secs=0.1, tau_rel=0.0))
    return recs


def fake_all_records():
    base = dict(gd=[0.9] * 4, shared_v=[0.9] * 4, adam_big=[0.9] * 4, adam=[0.1] * 4, sign=[0.1] * 4)
    return fake_runs(base, tids=range(1, 5)) + fake_p2_records()


def test_report_is_deterministic_and_contains_the_verdicts(tmp_path):
    meta = dict(stage='test', n_records=0)
    recs = fake_all_records()
    freeze = {a: dict(lr=1.0) for a in O.ARMS}
    t1, p1, p2 = O.report(meta, recs, freeze)
    t2, _, _ = O.report(meta, recs, freeze)
    assert t1 == t2
    assert p1['V'][0.3]['a_ok'] and p1['V'][0.3]['b_ok']
    assert p2['i_scalar'] is True and p2['i_net']['un'] == 'HELD' and p2['ii']['ok'] is True
    assert p2['ii']['sec_order'] is True and p2['ii']['sec_mono'] is True
    for needle in ('## Verdicts', 'P1 (a)', 'P2.1 scalar unit', 'P2.2 network, unequal masses', 'P2.3 formed-set agreement', 'POST HOC', 'saturated' if False else 'P2.3'):
        assert needle in t1


def test_report_flips_with_the_predictions():
    meta = dict(stage='test')
    freeze = {a: dict(lr=1.0) for a in O.ARMS}
    recs = fake_runs(dict(gd=[0.5] * 4, shared_v=[0.9] * 4, adam_big=[0.9] * 4, adam=[0.5] * 4, sign=[0.5] * 4), tids=range(1, 5)) + fake_p2_records()
    _, p1, p2 = O.report(meta, recs, freeze)
    assert p1['V'][0.3]['a_ok'] is False and p1['V'][0.3]['b_ok'] is False
    # P2.3 ordering broken: deeper readout worse
    r2 = fake_p2_records()
    for r in r2:
        if r['kind'] == 'p2c' and r['p'] == 5:
            r['set_equal'] = False
    _, _, p2 = O.report(meta, fake_runs(dict(gd=[0.9] * 4, shared_v=[0.9] * 4, adam_big=[0.9] * 4, adam=[0.1] * 4, sign=[0.1] * 4), tids=range(1, 5)) + r2, freeze)
    assert p2['ii']['ok'] is False
    # secondary: break the ordering at eps = 1e-3 only (D4 beyond D5) and the monotonicity in eps for D3
    r4 = fake_p2_records()
    for r in r4:
        if r['kind'] == 'p2c' and r['p'] == 4 and r['eps'] == 1e-3:
            r['set_equal'] = True
        if r['kind'] == 'p2c' and r['p'] == 5 and r['eps'] == 1e-3:
            r['set_equal'] = False
    _, _, p2b = O.report(meta, fake_runs(dict(gd=[0.9] * 4, shared_v=[0.9] * 4, adam_big=[0.9] * 4, adam=[0.1] * 4, sign=[0.1] * 4), tids=range(1, 5)) + r4, freeze)
    assert p2b['ii']['ok'] is True and p2b['ii']['sec_order'] is False
    r5 = fake_p2_records()
    for r in r5:
        if r['kind'] == 'p2c' and r['p'] == 3 and r['eps'] == 1e-8:
            r['set_equal'] = False
    _, _, p2c = O.report(meta, fake_runs(dict(gd=[0.9] * 4, shared_v=[0.9] * 4, adam_big=[0.9] * 4, adam=[0.1] * 4, sign=[0.1] * 4), tids=range(1, 5)) + r5, freeze)
    assert p2c['ii']['ok'] is False and p2c['ii']['sec_mono'] is False
    # P2.2: wrong exponent at p = 4
    r3 = fake_p2_records()
    for r in r3:
        if r['kind'] == 'p2n' and r['p'] == 4:
            r['wmax'] = [w * (r['eps'] ** 0.3) if w < 0.4 else w for w in r['wmax']]
    _, _, p2 = O.report(meta, fake_runs(dict(gd=[0.9] * 4, shared_v=[0.9] * 4, adam_big=[0.9] * 4, adam=[0.1] * 4, sign=[0.1] * 4), tids=range(1, 5)) + r3, freeze)
    assert p2['i_net']['un'] == 'FAILED'


def test_inconclusive_when_too_few_targets_pass_the_margin():
    recs = fake_p2_records()
    for r in recs:
        if r['kind'] == 'p2n' and r['p'] == 3 and r['tgt'] >= 1:
            r['status'] = 'failed: Requ'
    _, p2 = O.p2_report(recs, dict(stage='test'))
    assert p2['i_net']['un'] == 'INCONCLUSIVE'


def test_raw_roundtrip_and_cli_regeneration_is_byte_identical(tmp_path):
    meta = dict(stage='test', n_records=0)
    recs = fake_all_records()
    raw = str(tmp_path / 'optdepth_raw.npz')
    O.save_raw(raw, meta, recs)
    m2, r2 = O.load_raw(raw)
    assert m2 == meta and r2 == json.loads(json.dumps(recs, sort_keys=True))
    freeze = {a: dict(lr=1.0) for a in O.ARMS}
    fz = str(tmp_path / 'freeze.json')
    with open(fz, 'w') as f:
        json.dump(dict(p1=freeze), f)
    outs = []
    for k in range(2):
        md = str(tmp_path / ('o%d.md' % k))
        subprocess.run([sys.executable, SCRIPT, '--from-raw', raw, '--md', md, '--freeze', fz], check=True, env=dict(os.environ, OD_OUT=str(tmp_path)))
        outs.append(open(md).read())
    assert outs[0] == outs[1] and outs[0] == O.report(m2, r2, freeze)[0]


def test_dev_report_runs_on_planted_dev_records():
    recs = []
    truth = [0, 1]
    for arm in O.ARMS[1:]:
        for cfg in O.DEV_GRID[arm]:
            for t, seed in ((401, 40), (402, 41)):
                for noise in O.NOISES:
                    recs.append(rec(arm, t, noise, seed, [[], [0.0]], truth, json.dumps(cfg, sort_keys=True)))
    for t, seed in ((401, 40), (402, 41)):
        for noise in O.NOISES:
            recs.append(rec('gd', t, noise, seed, [[], [0.0]], truth, json.dumps(O.GD_CFG, sort_keys=True)))
            recs.append(rec('ctrl', t, noise, seed, [[], [0.0]], truth, json.dumps(dict(lam=0.0))))
    txt = O.dev_report(dict(stage='dev'), recs)
    assert '# Optimizer and depth: development report' in txt and 'Chosen:' in txt


def test_jsonl_helpers_and_jsonable(tmp_path):
    p = str(tmp_path / 'x.jsonl')
    O.append_jsonl(p, [dict(a=np.float64(1.5), b=np.int64(2), c=np.array([1, 2]), d=np.bool_(True), e=(1, 2))])
    O.append_jsonl(p, [dict(a=2)])
    r = O.read_jsonl(p)
    assert r[0] == dict(a=1.5, b=2, c=[1, 2], d=True, e=[1, 2]) and len(r) == 2 and O.read_jsonl(str(tmp_path / 'none.jsonl')) == []


def test_p1_job_smoke_records_are_complete_and_deterministic():
    arms = {a: [O.DEV_GRID[a][0]] for a in O.DEV_GRID}
    arms['gd'] = [O.GD_CFG]
    old = O.STEPS
    try:
        O.STEPS = 150
        recs, _ = O.p1_job((401, 0.1, 40, {'gd': [O.GD_CFG], 'sign': [O.DEV_GRID['sign'][0]]}, False))
        recs2, _ = O.p1_job((401, 0.1, 40, {'gd': [O.GD_CFG], 'sign': [O.DEV_GRID['sign'][0]]}, False))
    finally:
        O.STEPS = old
    assert [r['arm'] for r in recs] == ['gd', 'sign']
    for r, r2 in zip(recs, recs2):
        assert r['counts'] == r2['counts'] and r['centres'] == r2['centres'] and r['val_mse'] == r2['val_mse']
        assert [len(c) for c in r['centres']] == r['counts'] and len(r['truth']) == 5 and r['P'] == 3 * sum(r['counts']) + 1
    assert recs[0]['truth'] == KS.make_target(401)['K'].tolist()


def test_control_runs_prune_and_fine_tune(monkeypatch):
    d = small_data()
    monkeypatch.setattr(O, 'STEPS', 200)
    monkeypatch.setattr(O, 'FT_STEPS', 50)
    monkeypatch.setattr(O, 'LAMS', (0.0, 0.01))
    monkeypatch.setattr(O.KS, 'D_SYN', 2)
    tg = dict(K=np.array([1, 0]), bumps=[np.array([[1.2, 0.8, 0.4]]), np.zeros((0, 3))])
    recs = O.control_runs(d, tg, (999, 0.1, 0), 0)
    assert len(recs) == 2 and [json.loads(r['conf'])['lam'] for r in recs] == [0.0, 0.01] and all(r['arm'] == 'ctrl' for r in recs)
    assert all(len(r['centres']) == 2 for r in recs)


def test_stored_report_regenerates(tmp_path):
    raw, md = os.path.join(RES, 'optdepth_raw.npz'), os.path.join(RES, 'optdepth.md')
    if not (os.path.exists(raw) and os.path.exists(md)):
        pytest.skip('the held-out report is not stored yet')
    out = str(tmp_path / 'o.md')
    subprocess.run([sys.executable, SCRIPT, '--from-raw', raw, '--md', out], check=True)
    assert open(out).read() == open(md).read()


def test_stored_dev_report_regenerates(tmp_path):
    raw, md = os.path.join(RES, 'optdepth_dev_raw.npz'), os.path.join(RES, 'optdepth_dev.md')
    if not (os.path.exists(raw) and os.path.exists(md)):
        pytest.skip('the development report is not stored yet')
    out = str(tmp_path / 'o.md')
    subprocess.run([sys.executable, SCRIPT, '--from-raw', raw, '--md', out], check=True)
    assert open(out).read() == open(md).read()


# ---------------------------------------------------------------------------------------------------------------------------------
# tests added after the first mutation run (survivors of the first pass)
# ---------------------------------------------------------------------------------------------------------------------------------
def test_free_init_scale_and_start():
    m = O.AddNet(5, 10, 'free', seed=3)
    w = m.p['w']
    assert 0.07 < w.std() < 0.13 and abs(w.mean()) < 0.05
    assert np.array_equal(w, O.AddNet(5, 10, 'free', seed=3).p['w']) and not np.array_equal(w, O.AddNet(5, 10, 'free', seed=4).p['w'])


def test_fit_matches_theory3f_through_early_stopping():
    import theory3f_escape as T
    d = small_data(1)
    ma, mb = O.AddNet(2, 6, 'fact', 1e-2), T.AddKan(2, 6, 1e-2)
    sa = O.fit(ma, d, O.Opt('adam', 0.03), 12000)
    sb = T.fit_add(mb, d, T.Opt('adam', 0.03), 12000)
    assert 0 < sa[0] < 12000 and sa[0] == sb and all(np.array_equal(ma.p[k], mb.p[k]) for k in mb.p)      # same stopping step and best iterate


def test_activity_default_threshold_is_the_spike_rule():
    m = O.AddNet(3, 6, 'fact', 1e-2)
    rng = np.random.default_rng(2)
    m.p['a'] = rng.uniform(-1, 1, m.p['a'].shape) * 10 ** rng.uniform(-4, 0, m.p['a'].shape)
    X = rng.uniform(-2.5, 2.5, (200, 3))
    ref = KS.layer_activity(m.weights()[None], m.p['c'][None], m.widths()[None], X, KS.TAU)
    assert np.array_equal(m.activity(X)[2], ref[2]) and np.array_equal(m.counts(X), KS.edge_counts(ref[2]))
    assert not np.array_equal(KS.layer_activity(m.weights()[None], m.p['c'][None], m.widths()[None], X, 0.2)[2], ref[2])     # the threshold matters on this model
    d = small_data()
    tg = dict(K=np.array([1, 0]), bumps=[np.array([[1.2, 0.8, 0.4]]), np.zeros((0, 3))])
    r = O.summarise(m.__class__(2, 6, 'fact', 1e-2), d, tg, (1, 0.5), 'gd', '{}', (1, 0.1, 0))
    assert r['counts'] == [0, 0] and r['P'] == 1


def test_summarise_records_only_active_centres_and_count_matches():
    d = small_data()
    m = O.AddNet(2, 6, 'fact', 1e-2)
    m.p['a'][0, 2] = 2.0
    m.p['a'][1, 4] = -2.0
    m.p['n'] = m.p['n'] * 0 + 0.5
    tg = dict(K=np.array([1, 1]), bumps=[np.array([[1.2, 0.8, 0.4]]), np.array([[1.0, -1.0, 0.4]])])
    r = O.summarise(m, d, tg, (7, 0.25), 'gd', '{"lr": 1.0}', (5, 0.1, 3))
    assert r['counts'] == [len(c) for c in r['centres']] and r['P'] == 3 * sum(r['counts']) + 1
    assert all(len(c) <= 6 for c in r['centres']) and r['steps'] == 7 and r['val_mse'] == 0.25 and r['tid'] == 5 and r['seed'] == 3
    assert r['truth'] == [1, 1] and r['truth_c'] == [[0.8], [-1.0]]
    assert sum(r['counts']) < 12                                     # inactive units are not recorded


def test_control_runs_call_sequence(monkeypatch):
    calls = []
    real_fit = O.fit

    def spy(model, data, opt, steps, l1=0.0, patience=O.PATIENCE):
        calls.append(dict(kind=opt.kind, lr=opt.lr, l1=l1, steps=steps, mask=model.mask.copy()))
        return real_fit(model, data, opt, 60, l1=l1, patience=patience)
    monkeypatch.setattr(O, 'fit', spy)
    monkeypatch.setattr(O, 'LAMS', (0.0, 0.01))
    monkeypatch.setattr(O.KS, 'D_SYN', 2)
    d = small_data()
    tg = dict(K=np.array([1, 0]), bumps=[np.array([[1.2, 0.8, 0.4]]), np.zeros((0, 3))])
    O.control_runs(d, tg, (999, 0.1, 0), 0)
    assert len(calls) == 4
    c0, c1, c2, c3 = calls
    assert (c0['kind'], c0['lr'], c0['l1'], c0['steps']) == ('adam', O.CTRL_LR, 0.0, O.STEPS) and (c2['l1'], c2['lr']) == (0.01, O.CTRL_LR)
    assert (c1['kind'], c1['l1'], c1['steps']) == ('adam', 0.0, O.FT_STEPS) and c1['lr'] == pytest.approx(O.CTRL_LR / 3.0)
    assert np.all(c0['mask'] == 1.0) and c1['mask'].sum() < c1['mask'].size                  # pruned before the fine-tune


def test_oracle_control_uses_the_radius_and_differs_from_validation_pick():
    truth = [1]
    lam = lambda x: json.dumps(dict(lam=x))
    A = rec('ctrl', 1, 0.1, 0, [[0.0, 0.1]], truth, lam(0.0), val=0.9, P=7)           # exact after merging at 0.3 only
    B = rec('ctrl', 1, 0.1, 0, [[0.0, 1.0]], truth, lam(0.01), val=1.0, P=5)          # never exact
    assert O.pick_ctrl([A, B], 0.3)['conf'] == lam(0.0)
    assert O.pick_ctrl([A, B], 0.0)['conf'] == lam(0.01)                               # unmerged: a tie at 0, the fewer parameters win
    runs = O.complete_runs(O.group_runs([rec(a, 1, 0.1, 0, [[0.0]], truth) for a in O.ARMS] + [A, B]), list(O.ARMS) + ['ctrl'])
    fr = {a: dict(lr=1.0) for a in O.ARMS}
    assert O.arm_table(runs, 'ctrl_or', 0.3, fr)['share_pos'] == 1.0
    assert O.arm_table(runs, 'ctrl', 0.3, fr)['share_pos'] == 1.0                      # one-SE: A is the best validation
    A2 = dict(A, val_mse=5.0)
    runs2 = O.complete_runs(O.group_runs([rec(a, 1, 0.1, 0, [[0.0]], truth) for a in O.ARMS] + [A2, B]), list(O.ARMS) + ['ctrl'])
    assert O.arm_table(runs2, 'ctrl', 0.3, fr)['share_pos'] == 0.0 and O.arm_table(runs2, 'ctrl_or', 0.3, fr)['share_pos'] == 1.0


def test_complete_runs_drops_runs_missing_an_arm():
    truth = [1]
    recs = [rec(a, 1, 0.1, 0, [[0.0]], truth) for a in O.ARMS] + [rec(a, 2, 0.1, 0, [[0.0]], truth) for a in O.ARMS[:-1]]
    recs += [rec('ctrl', 1, 0.1, 0, [[0.0]], truth), rec('ctrl', 2, 0.1, 0, [[0.0]], truth)]
    runs = O.complete_runs(O.group_runs(recs), list(O.ARMS) + ['ctrl'])
    assert list(runs) == [(1, 0.1, 0)]


def test_arm_table_surplus_and_missing_per_edge():
    truth = [0, 2, 1]
    recs = [rec(a, 1, 0.1, 0, [[], [0.0, 1.0, 2.0], [0.0, 1.0]], truth) for a in O.ARMS] + [rec('ctrl', 1, 0.1, 0, [[], [0.0], [0.0]], truth)]
    recs += [rec(a, 1, 0.1, 1, [[5.0], [0.0], [0.0]], truth) for a in O.ARMS] + [rec('ctrl', 1, 0.1, 1, [[5.0], [0.0], [0.0]], truth)]
    runs = O.complete_runs(O.group_runs(recs), list(O.ARMS) + ['ctrl'])
    t = O.arm_table(runs, 'gd', 0.0, {a: dict(lr=1.0) for a in O.ARMS})
    # run 0: K>0 edges: counts 3 vs 2 (surplus 1), 2 vs 1 (surplus 1); run 1: 1 vs 2 (missing 1), 1 vs 1: per K>0 edge = 4 edges
    assert t['surplus_pos'] == pytest.approx(2 / 4) and t['missing_pos'] == pytest.approx(1 / 4) and t['share_empty'] == pytest.approx(0.5)
    assert t['share_all'] == pytest.approx(2 / 6) and t['share_pos'] == pytest.approx(1 / 4)


def test_contrast_uses_the_declared_bootstrap_and_pairs_by_target():
    base = dict(gd=[1.0, 0.5, 0.0, 1.0], shared_v=[1.0] * 4, adam_big=[1.0] * 4, adam=[0.0, 0.5, 0.0, 0.5], sign=[0.0] * 4)
    runs = O.complete_runs(O.group_runs(fake_runs(base, tids=range(1, 5))), list(O.ARMS) + ['ctrl'])
    fr = {a: dict(lr=1.0) for a in O.ARMS}
    m, lo, hi, d = O.contrast(runs, 'gd', 'adam', 0.3, fr)
    assert d == pytest.approx([1.0, 0.0, 0.0, 0.5]) and (m, lo, hi) == KS.boot_ci(d, n_boot=O.N_BOOT, seed=O.BOOT_SEED)
    assert O.N_BOOT == 20000 and O.BOOT_SEED == 20261110


def test_decision_rule_boundaries_are_inclusive_and_arms_are_the_registered_ones(monkeypatch):
    hi, lo = [1.0] * 5, [0.0] * 5
    base = dict(gd=hi, shared_v=[4 / 6] * 5, adam_big=[5 / 6] * 5, adam=[2 / 6] * 5, sign=[0.0] * 5)
    runs = O.complete_runs(O.group_runs(fake_runs(base)), list(O.ARMS) + ['ctrl'])
    fr = {a: dict(lr=1.0) for a in O.ARMS}
    V = O.p1_verdict(runs, 0.3, fr)
    assert [(x['hi'], x['lo']) for x in V['a']] == [('gd', 'adam'), ('gd', 'sign'), ('shared_v', 'adam'), ('shared_v', 'sign')]
    assert V['a'][0]['mean'] == pytest.approx(4 / 6) and V['a'][2]['mean'] == pytest.approx(2 / 6) and V['a'][3]['mean'] == pytest.approx(4 / 6)
    assert V['b']['mean'] == pytest.approx(2 / 6)                                       # GD - shared_v
    assert V['d']['mean'] == pytest.approx(-1 / 6)                                      # adam_big - GD
    m_a = V['a'][2]['mean']
    monkeypatch.setattr(O, 'GATE_GAP', m_a)                                             # a gap exactly at the gate passes (>=)
    assert O.p1_verdict(runs, 0.3, fr)['a'][2]['ok'] is True
    monkeypatch.setattr(O, 'GATE_GAP', m_a + 1e-9)
    assert O.p1_verdict(runs, 0.3, fr)['a'][2]['ok'] is False
    monkeypatch.setattr(O, 'GATE_NEAR', V['b']['mean'])                                 # |GD - shared_v| exactly at the bound passes (<=)
    assert O.p1_verdict(runs, 0.3, fr)['b_ok'] is True
    monkeypatch.setattr(O, 'GATE_NEAR', V['b']['mean'] - 1e-9)
    assert O.p1_verdict(runs, 0.3, fr)['b_ok'] is False
    # (b) is the absolute difference: shared_v ahead of GD by the same amount passes the same way
    swap = dict(base, gd=[4 / 6] * 5, shared_v=hi)
    monkeypatch.setattr(O, 'GATE_NEAR', 2 / 6 + 1e-9)
    assert O.p1_verdict(O.complete_runs(O.group_runs(fake_runs(swap)), list(O.ARMS) + ['ctrl']), 0.3, fr)['b_ok'] is True


def test_verdict_table_rows_follow_their_flags():
    V = dict(a_ok=False, b_ok=True, d=dict(ok=True))
    p1 = dict(V={O.MERGE_R: V})
    p2 = dict(i_scalar=True, i_net=dict(un='INCONCLUSIVE', eq='HELD'), ii=dict(ok=True, sec_order=True, sec_mono=True))
    rows = O.verdict_table(p1, p2)
    row = lambda key: [r for r in rows if r.startswith('| ' + key)][0]
    assert row('P1 (a)').endswith('| FAILED |') and row('P1 (b)').endswith('| HELD |')
    assert row('P2.2 network, unequal masses').endswith('| INCONCLUSIVE |') and row('P2.2 network, equal masses').endswith('| HELD |')
    assert row('P2.3 formed-set').endswith('| HELD |') and row('P2.1').endswith('| HELD |')


def test_net_cell_wmax_is_over_time_and_masses_are_the_registered_draw():
    import exp_parametrization as PZ
    ref, tau_end = O.net_target_ref(list(PZ.CSTAR), list(PZ.ASTAR))
    c = O.net_cell(list(PZ.CSTAR), list(PZ.ASTAR), 3, 1e-4, O.RAND_NET, tau_end)
    wm, we = np.array(c['wmax']), np.array(c['wend'])
    assert np.all(wm >= we - 1e-300) and np.any(wm > 1.001 * we)                      # a dormant unit's dual peaks and falls back: the maximum over time exceeds the end value
    assert np.allclose(c['b'], np.exp(O.RAND_NET * np.random.default_rng(O.MASS_SEED).standard_normal(21)), rtol=1e-12)


def planted_p2c(set_equal_fn, status_fn=lambda p, e, i: 'ok', n=10):
    recs = []
    for i in range(n):
        recs.append(dict(kind='p2r', tgt=i, fam='A', cs=[0.0], As=[1.0], margin=0.3, margin_never=0.3, thr=0.01, order=[0], set=[0], relaxed=True))
        for p in O.DEPTHS_II:
            for e in O.EPS_II_ALL:
                ok = set_equal_fn(p, e, i)
                recs.append(dict(kind='p2c', tgt=i, p=p, eps=e, set_equal=bool(ok), order_equal=bool(ok), n_formed=1, exact=bool(ok), status=status_fn(p, e, i), secs=0.1, tau_rel=0.0))
    return recs


def test_p2_ordering_rule_strictness_and_monotonicity():
    meta = dict(stage='test')
    ii = lambda recs: O.p2_report(recs, meta)[1]['ii']
    base = lambda p, e, i: {2: i < 1, 3: i < 6, 4: i < 8, 5: i < 9}[p]
    assert ii(planted_p2c(base))['ok'] is True
    # f2 == f3 at a registered eps: strict inequality fails
    assert ii(planted_p2c(lambda p, e, i: (i < 6) if p in (2, 3) else base(p, e, i)))['ok'] is False
    # f3 == f4 is allowed (non-strict), f4 > f5 is not
    assert ii(planted_p2c(lambda p, e, i: (i < 6) if p in (3, 4) else base(p, e, i)))['ok'] is True
    assert ii(planted_p2c(lambda p, e, i: (i < 9) if p == 4 else ((i < 8) if p == 5 else base(p, e, i))))['ok'] is False
    # only the monotonicity in eps breaks (D5 slightly lower at 1e-8): the ordering is intact
    assert ii(planted_p2c(lambda p, e, i: (i < 8 if (p == 5 and e == 1e-8) else base(p, e, i))))['ok'] is False
    # a failed integration is a failure even if the stored set equals the path's
    r = ii(planted_p2c(base, status_fn=lambda p, e, i: 'failed: Requ' if (p == 5 and e == 1e-8 and i == 0) else 'ok'))
    assert r['ok'] is False and r['F']['5_1e-08'] == pytest.approx(0.8)


def test_p2_exponent_rule_width_and_scalar_rule_tolerance():
    meta = dict(stage='test')
    r = fake_p2_records()
    for c in r:
        if c['kind'] == 'p2n' and c['p'] == 4 and c['arm'] == 'un':
            c['wmax'] = [w * (c['eps'] ** 0.07) if w < 0.4 else w for w in c['wmax']]       # exponent 2.07: outside 0.05, inside twice that
    assert O.p2_report(r, meta)[1]['i_net']['un'] == 'FAILED'
    r = fake_p2_records()
    for c in r:
        if c['kind'] == 'p2s' and c['p'] == 3:
            c['wmax'] = c['wmax'] * c['eps'] ** 0.05
    assert O.p2_report(r, meta)[1]['i_scalar'] is False
    r = fake_p2_records()
    for c in r:
        if c['kind'] == 'p2s' and c['p'] == 4:
            c['exact'] = c['exact'] * c['eps'] ** 0.5                                       # only the 'predicted' column (closed form) moves; the rule uses p/2
            c['wmax'] = c['wmax'] * c['eps'] ** 0.5
    assert O.p2_report(r, meta)[1]['i_scalar'] is False                                     # exponent 2.5 for p = 4: not p/2 = 2


def test_stopping_is_armed_only_after_a_ten_percent_drop():
    rng = np.random.default_rng(0)
    n = 400
    X = rng.uniform(-2.5, 2.5, (n + 200, 2))
    y = 0.23 + rng.standard_normal(n + 200)                                  # a constant offset: only the bias can reduce the validation MSE, by about 5 percent
    d = dict(Xtr=X[:n], ytr=y[:n], Xva=X[n:n + 100], yva=y[n:n + 100], Xte=X[n + 100:], yte=y[n + 100:])
    m = O.AddNet(2, 4, 'fact', 1e-10)
    st = O.fit(m, d, O.make_opt('gd', O.GD_CFG, m), 800, patience=100)
    v0 = np.mean(y[n:n + 100] ** 2)
    assert st[1] > 0.9 * v0 and st[1] < v0                                   # the validation MSE improved, but never below 0.9 of its start ...
    assert st[0] == 800                                                      # ... so the run is never armed, and is not stopped by the patience
    # (an offset of 0.8 with a small noise drops the validation MSE below 0.9 of its start: armed, then stopped after `patience` steps without improvement)
    d2 = dict(d, yva=0.8 + rng.standard_normal(100) * 0.3, ytr=0.8 + rng.standard_normal(n) * 0.3)
    m2 = O.AddNet(2, 4, 'fact', 1e-10)
    st2 = O.fit(m2, d2, O.make_opt('gd', O.GD_CFG, m2), 800, patience=100)
    assert st2[0] < 800 and st2[0] % O.EVERY == 0


def test_select_arm_rec_oracle_default_radius_is_the_declared_one():
    truth = [1]
    lam = lambda x: json.dumps(dict(lam=x))
    A = rec('ctrl', 1, 0.1, 0, [[0.0, 0.1]], truth, lam(0.0), val=5.0, P=7)
    B = rec('ctrl', 1, 0.1, 0, [[0.0, 1.0]], truth, lam(0.01), val=1.0, P=5)
    assert O.select_arm_rec([A, B], 'ctrl')['conf'] == lam(0.01)
    assert O.select_arm_rec([A, B], 'ctrl_or')['conf'] == lam(0.0) and O.select_arm_rec([A, B], 'ctrl_or', None, 0.3)['conf'] == lam(0.0)
    assert O.select_arm_rec([A, B], 'ctrl_or', None, 0.0)['conf'] == lam(0.01)
