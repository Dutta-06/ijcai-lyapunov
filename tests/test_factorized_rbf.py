"""Tests of src/core/factorized_rbf.py: the free law as the Euclidean gradient flow of a factorised RBF network, the other readouts, the
trainers (flow, fixed-step descent, Adam) and the comparison with the sliding path."""
import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.core import factorized_rbf as FR
from src.core import rich_transit as RT
from src.core import rich_transit_nd as RN
from src.core.mirror_free import FreeLaw

warnings.filterwarnings('ignore')
LAW = FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
LR_C = 12.0
SIG = 0.12


@pytest.fixture(scope='module')
def design():
    xs = np.linspace(-1, 1, 61)
    K = RN.BumpsND(xs[:, None], SIG)
    c0 = np.linspace(-1, 1, 21)
    cs, As = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
    f = K.synth(cs[:, None], As)
    R = RN.sliding_iss_nd(LAW, K, f, c0[:, None], np.ones(21), polish=True)
    ent = [(t, j) for t, j, k, c in R['events'] if k == 'enter']
    ref = dict(order=[j for t, j in ent], tau=[t for t, j in ent], set=sorted(int(j) for j in np.flatnonzero(R['act'])))
    return K, c0, cs, As, f, ref


READOUTS = [FR.Linear(2.0), FR.Product2(0.7), FR.Cubic3(1.3), FR.PowerN(4, 0.9), FR.PowerN(5, 1.1), FR.PaperLaw(LAW), FR.Unsaturated(LAW)]


def F_of(readout, K, f, p, c):
    r = f - K.synth(c, readout.w(p))
    return 0.5 * np.mean(r ** 2)


@pytest.mark.parametrize('ro', READOUTS, ids=lambda r: r.name)
def test_dw_is_the_derivative_of_w(ro):
    rng = np.random.default_rng(1)
    p = rng.uniform(0.2, 1.1, (6, ro.q)) * rng.choice([-1, 1], (6, ro.q))
    h = 1e-6
    for l in range(ro.q):
        e = np.zeros_like(p)
        e[:, l] = h
        fd = (ro.w(p + e) - ro.w(p - e)) / (2 * h)
        assert np.allclose(ro.dw(p)[:, l], fd, rtol=1e-6, atol=1e-9)


@pytest.mark.parametrize('ro', READOUTS, ids=lambda r: r.name)
def test_the_flow_is_the_negative_gradient_of_the_fit_times_the_learning_rates(ro, design):
    """ṗ = −lr·∂F/∂p and ċ = −lr_c·∂F/∂c, checked against central differences of F itself."""
    K, c0, cs, As, f, ref = design
    rng = np.random.default_rng(2)
    d = 5
    p = rng.uniform(0.3, 0.9, (d, ro.q))
    c = np.sort(rng.uniform(-0.8, 0.8, d))[:, None]
    z = np.r_[p.ravel(), c.ravel()]
    got = FR.make_rhs(ro, K, f, LR_C, d, 1)(0.0, z)
    h = 1e-6
    fd = np.zeros(z.size)
    for k in range(z.size):
        e = np.zeros(z.size)
        e[k] = h
        pp, cc = (z + e)[:d * ro.q].reshape(d, ro.q), (z + e)[d * ro.q:].reshape(d, 1)
        pm, cm = (z - e)[:d * ro.q].reshape(d, ro.q), (z - e)[d * ro.q:].reshape(d, 1)
        fd[k] = (F_of(ro, K, f, pp, cc) - F_of(ro, K, f, pm, cm)) / (2 * h)
    rates = np.r_[np.tile(ro.lr, d), np.full(d, LR_C)]
    assert np.allclose(got, -rates * fd, rtol=2e-6, atol=1e-10)


def test_the_paper_law_is_the_euclidean_flow_in_alpha_nu_c(design):
    """The brief's reformulation: `moving_rhs` in (α, log β, c) equals the engine's flow in (α, ν = √β, c) with the learning rates
    (Gm, γm/4, Gm), exactly (an identity, not a numerical agreement)."""
    K, c0, cs, As, f, ref = design
    Kd = RT.Bumps(K.xs[:, 0], K.sig)
    repo = RT.moving_rhs(LAW, Kd, f, 'paradigm')
    ro = FR.PaperLaw(LAW)
    assert np.allclose(ro.lr, [12.0, 0.3], rtol=1e-14, atol=0)  # (Gm, γm/4)
    assert np.isclose(ro.lr[1] / ro.lr[0] * LAW.A ** 2, LAW.c, rtol=1e-14, atol=0)  # v_s = a_max² lr_ν/lr_α
    rng = np.random.default_rng(3)
    d = 8
    for _ in range(4):
        al, beta, c = rng.normal(0, 0.5, d), np.exp(rng.normal(-1, 1, d)), np.sort(rng.uniform(-0.9, 0.9, d))
        mine = FR.make_rhs(ro, K, f, LR_C, d, 1)(0.0, np.r_[np.stack([al, np.sqrt(beta)], 1).ravel(), c])
        theirs = repo(0.0, np.r_[al, np.log(beta), c])
        assert np.allclose(mine[:2 * d:2], theirs[:d], rtol=1e-12, atol=0)
        assert np.allclose(2 * mine[1:2 * d:2] / np.sqrt(beta), theirs[d:2 * d], rtol=1e-12, atol=0)
        assert np.allclose(mine[2 * d:], theirs[2 * d:], rtol=1e-12, atol=0)


def test_ratio_of_learning_rates_fixes_vs():
    for g in (0.1, 1.0, 4.0):
        L = FreeLaw(m=12.0, G=1.0, gamma=g, a_max=1.0)
        ro = FR.PaperLaw(L)
        assert np.isclose(ro.lr[1] / ro.lr[0], g / 4.0, rtol=1e-14, atol=0)
        assert np.isclose(L.c, L.A ** 2 * ro.lr[1] / ro.lr[0], rtol=1e-14, atol=0)
    assert np.allclose(FR.PaperLaw(FreeLaw(m=12.0, G=1.0, gamma=4.0, a_max=1.0)).lr, [12.0, 12.0], rtol=0, atol=0)   # equal rates


@pytest.mark.parametrize('name', ['P0', 'P0u', 'D2', 'D3', 'D4'])
def test_conserved_quantities_along_the_flow(name, design):
    K, c0, cs, As, f, ref = design
    ro = {'P0': FR.PaperLaw(LAW), 'P0u': FR.Unsaturated(LAW), 'D2': FR.Product2(1.0), 'D3': FR.Cubic3(1.0), 'D4': FR.PowerN(4, 1.0)}[name]
    eps, d = 1e-3, 21
    sc = ro.scale(eps)
    from scipy.integrate import solve_ivp
    z0 = np.r_[ro.init(d, eps).ravel(), c0]
    sol = solve_ivp(FR.make_rhs(ro, K, f, LR_C, d, 1), (0, 6 * sc), z0, method='LSODA', rtol=1e-11, atol=1e-16, t_eval=np.linspace(0, 6 * sc, 40))
    P = sol.y[:d * ro.q].reshape(d, ro.q, -1)
    if name == 'P0':
        q = P[:, 1] ** 2 - LAW.c * np.cosh(2 * P[:, 0] / LAW.A)                       # Q = β − v_s cosh(2α/a_max)
        scale = LAW.c
    elif name == 'P0u':
        q = P[:, 1] ** 2 - LAW.gamma / (2 * LAW.G) * P[:, 0] ** 2                      # ν² − (γ/2G)α²
        scale = eps
    elif name == 'D2':
        q, scale = P[:, 0] * P[:, 1], eps                                              # uv = s²
    elif name == 'D3':
        q, scale = 1 / P[:, 0] + 1 / P[:, 1], 2 / np.sqrt(eps)                         # 1/u + 1/v = 2/s
    else:
        q, scale = P[:, 0] ** -2 + P[:, 1] ** -2, 2 / eps                              # u^{−2} + v^{−2} = 2/s²
    assert np.abs(q - q[:, :1]).max() <= 1e-7 * scale


def test_walls():
    """`scale`: D2's wall is where w = 1 (to O(ε²)), D3's coefficient at the margin μ is s³/μ³ − s³/8 (u, v from the conserved 1/u + 1/v)."""
    eps, lr = 1e-6, 0.7
    r2 = FR.Product2(lr)
    s = np.sqrt(eps)
    eta = r2.scale(eps)
    assert np.isclose(2 * s ** 2 * np.sinh(4 * lr * eta), 1.0, rtol=1e-11, atol=0)       # w = 2s² sinh(4 lr η)
    assert np.isclose(r2.scale(eps ** 2), 2 * r2.scale(eps), rtol=1e-14, atol=0)          # logarithmic in ε
    r3 = FR.Cubic3(lr)
    assert np.isclose(r3.scale(eps), 1 / (3 * lr * s), rtol=1e-14, atol=0)
    assert np.isclose(r3.scale(eps / 100), 10 * r3.scale(eps), rtol=1e-14, atol=0)        # a power law ε^{−1/2}
    mu = 0.3
    e = (1 - mu) * r3.scale(eps)
    u, v = 1 / (1 / s - 3 * lr * e), 1 / (1 / s + 3 * lr * e)            # = s/μ and s/(2 − μ)
    assert np.isclose(u, s / mu, rtol=1e-12, atol=0) and np.isclose(v, s / (2 - mu), rtol=1e-12, atol=0)
    assert np.isclose(r3.w(np.array([[u, v]]))[0], s ** 3 * (1 / mu ** 3 - 1 / (2 - mu) ** 3), rtol=1e-9, atol=0)


def test_formation_times_interpolate_log_linearly():
    taus = np.linspace(0, 10, 101)
    W = np.stack([np.exp(taus - 5.0) * 1e-3, np.zeros_like(taus), 0.5 * np.ones_like(taus)], 1)      # exponential, never, always
    ft = FR.formation_times(taus, W, thr=1e-2)
    assert np.isclose(ft[0], 5.0 + np.log(10.0), rtol=1e-12, atol=0)       # exp(τ − 5)/1000 = 1/100 at τ = 5 + ln 10 (log-linear: exact)
    assert np.isnan(ft[1])
    assert ft[2] == 0.0


def test_match_features_and_analyse_on_a_synthetic_record(design):
    K, c0, cs, As, f, ref = design
    T, d = 5, 21
    taus = np.linspace(0, 1, T)
    W = np.zeros((T, d))
    C = np.tile(c0[None, :, None], (T, 1, 1)).astype(float)
    W[2:, 5], W[3:, 16], W[4:, 11] = As[0], As[2], As[1]          # units 5, 16, 11 form at τ ≈ 0.5, 0.75, 1 and sit on their features at the end
    for j, c in zip([5, 11, 16], [cs[0], cs[1], cs[2]]):
        C[-1, j, 0] = c
    W[:, 3] = 2e-4                                    # a dormant unit: above the 1e-4 sweep level (9e-5), below the 1e-3 one, displaced by 0.3σ
    C[-1, 3, 0] += 0.3 * SIG
    res = dict(taus=taus, W=W, C=C)
    a = FR.analyse(res, K, f, c0, cs, As, ref=dict(order=[5, 16, 11], tau=[0.5, 0.75, 1.0], set=[5, 11, 16]))
    assert a['n_formed'] == 3 and a['exact'] and a['set_equal']
    assert np.isclose(a['dormant_disp'], 0.3, rtol=1e-9, atol=0)
    assert np.isclose(a['dormant_wmax'], 2e-4, rtol=1e-12, atol=0)
    assert np.isclose(a['gap'], np.log10(0.6 / 2e-4), rtol=1e-9, atol=0)
    assert a['form_order'][:3] == [5, 16, 11] and a['order_equal']
    assert a['n_formed_sweep'][1e-4] == 4 and a['n_formed_sweep'][1e-3] == 3 and a['n_formed_sweep'][1e-1] == 3
    # one formed unit displaced by more than 0.05σ makes the fit inexact
    C2 = C.copy()
    C2[-1, 11, 0] += 0.2 * SIG
    assert not FR.analyse(dict(taus=taus, W=W, C=C2), K, f, c0, cs, As)['exact']
    # a fourth formed unit: no longer exact, the set differs from the path's
    W2 = W.copy()
    W2[-1, 3] = 0.2
    b = FR.analyse(dict(taus=taus, W=W2, C=C), K, f, c0, cs, As, ref=dict(order=[5, 16, 11], tau=[0.5, 0.75, 1.0], set=[5, 11, 16]))
    assert b['n_formed'] == 4 and not b['exact'] and not b['set_equal']


def test_match_features_requires_distinct_units():
    ok, dist, aerr = FR.match_features(np.array([[0.0], [1.0]]), np.array([1.0, 2.0]), [0.1, 0.9], np.array([1.0, 2.0]))
    assert ok and np.allclose(dist, [0.1, 0.1], rtol=1e-12, atol=0) and np.allclose(aerr, 0.0, rtol=0, atol=0)
    ok, _, _ = FR.match_features(np.array([[0.0]]), np.array([1.0]), [0.1, -0.1], np.array([1.0, 1.0]))
    assert not ok                                   # both features nearest to the same unit


# -- the claims: which readouts reproduce the path (pinned at the design, eps = 1e-6) --------------------------------------------
def _run(ro, design, eps, tau_end=20.0, **kw):
    K, c0, cs, As, f, ref = design
    res = FR.train_flow(ro, K, f, c0, eps, tau_end, LR_C, n_t=801, **kw)
    return FR.analyse(res, K, f, c0, cs, As, ref=ref), res


@pytest.mark.parametrize('name', ['P0', 'P0u', 'D3'])
def test_depth_three_readouts_reproduce_the_sliding_path(name, design):
    ro = {'P0': FR.PaperLaw(LAW), 'P0u': FR.Unsaturated(LAW), 'D3': FR.Cubic3(1.0)}[name]
    a, res = _run(ro, design, 1e-6)
    assert res['status'] == 'ok'
    assert a['n_formed'] == 3 and a['set_equal'] and a['order_equal'] and a['exact']
    assert a['gap'] > 4.5                       # five decades between the formed and the largest dormant coefficient
    assert a['dormant_disp'] < 0.05             # the dormant units stay put (kernel widths)
    assert a['fit_rms'] < 1e-6
    assert a['tau_rel'] < 0.02                  # the entry times are the path's


def test_depth_two_does_not_select_at_a_practical_scale_and_converges_logarithmically(design):
    a6, _ = _run(FR.Product2(1.0), design, 1e-6)
    a12, _ = _run(FR.Product2(1.0), design, 1e-12)
    assert a6['n_formed'] >= 6 and not a6['set_equal'] and not a6['exact']
    assert a12['n_formed'] >= 4 and not a12['exact']
    assert a6['gap'] < 2 and a12['gap'] < 2     # no spectral gap between the formed and the dormant coefficients
    a96, _ = _run(FR.Product2(1.0), design, 1e-96)
    assert a96['n_formed'] == 3 and a96['set_equal'] and a96['exact']
    # τ error ∝ 1/ln(1/ε): doubling ln(1/ε) halves it (6 → 12 → 96 digits: ln ratio 1 : 2 : 16)
    assert 1.6 < a6['tau_rel'] / a12['tau_rel'] < 2.4
    assert 6.0 < a12['tau_rel'] / a96['tau_rel'] < 10.0


def test_depth_one_control_forms_many_units(design):
    K, c0, cs, As, f, ref = design
    sc = RT.transit_scale(LAW, 1e-6)
    a, res = _run(FR.Linear(LR_C), design, 1e-6, t_scale=sc)
    assert a['n_formed'] >= 8 and not a['set_equal']


def test_fixed_step_descent_follows_the_flow(design):
    K, c0, cs, As, f, ref = design
    ro = FR.Cubic3(1.0)
    eps, tau_end = 1e-4, 18.0                          # all three entries (the path's τ = 10.7, 12.8, 16.4) happen before τ_end
    flow = FR.train_flow(ro, K, f, c0, eps, tau_end, LR_C, n_t=451, rtol=1e-10)
    gd = FR.train_gd(ro, K, f, c0, eps, tau_end, LR_C, h=0.02, n_t=451)
    assert gd['status'] == 'ok'
    assert np.abs(gd['W'][-1] - flow['W'][-1]).max() < 1e-2 * np.abs(flow['W'][-1]).max()      # the final state (not inside the entry layers)
    assert np.abs(gd['C'][-1] - flow['C'][-1]).max() < 1e-2 * SIG
    thr = FR.THETA * np.abs(As).max()
    ft_gd, ft_fl = FR.formation_times(gd['taus'], gd['W'], thr), FR.formation_times(flow['taus'], flow['W'], thr)
    assert not np.isnan(ft_fl[[5, 11, 16]]).any() and not np.isnan(ft_gd[[5, 11, 16]]).any()
    assert np.allclose(ft_gd[[5, 11, 16]], ft_fl[[5, 11, 16]], rtol=2e-2, atol=0)


def test_gd_flags_divergence(design):
    """A depth-1 readout with the step h·lr·⟨φ²⟩ ≫ 2 (⟨φ²⟩ ≈ 0.3) overshoots; the trainer stops at the threshold and says so."""
    K, c0, cs, As, f, ref = design
    # (the overshooting kernels run off the samples and stop contributing, so the iteration saturates near 3e4 instead of overflowing:
    # the flag is exercised through its threshold `blow`)
    gd = FR.train_gd(FR.Linear(1e3), K, f, c0, 1e-4, 1.0, LR_C, h=1.0, n_t=101, t_scale=50.0, blow=1e3)
    assert gd['status'] == 'diverged'
    ok = FR.train_gd(FR.Linear(1.0), K, f, c0, 1e-4, 1.0, LR_C, h=0.05, n_t=101, t_scale=50.0)
    assert ok['status'] == 'ok' and np.isfinite(ok['W']).all()


def test_adam_first_step_is_the_sign_of_the_gradient_times_the_step(design):
    """At step 1, m̂/(√v̂ + ε_adam) = g/(|g| + ε_adam), so every parameter moves by lr_adam·sign(−∂F/∂θ) (to ε_adam/|g|). For D3 from (s, s),
    s = 1e-3, lr_adam ≪ s: Δu = +lr·sign(g_i), Δv = −lr·sign(g_i), hence w_i = u³ − v³ moves by 6s²·lr·sign(g_i) (to O(lr/s)). (With
    Adam's default ε_adam = 1e-8 the step is lr·|g′|/(|g′| + ε_adam), g′ = 3s²g ≈ 1e-8: not a sign step — the dormant parameters' gradients are the
    size of Adam's ε.)"""
    K, c0, cs, As, f, ref = design
    ro = FR.Cubic3(1.0)
    lr_adam, eps = 1e-9, 1e-6
    res = FR.train_adam(ro, K, f, c0, eps, n_steps=1, lr_adam=lr_adam, n_t=2, eps_adam=1e-40)   # Adam's default ε = 1e-8 is the size of 3s²g here
    g = K.corr(c0[:, None], f)[0]                      # at the start w = 0, r = f*
    s = np.sqrt(eps)
    big = np.abs(g) > 1e-3 * np.abs(g).max()
    assert big.sum() >= 10
    assert np.allclose(res['W'][-1][big], 6 * s ** 2 * lr_adam * np.sign(g[big]), rtol=1e-3, atol=0)


def test_timeout_is_reported_and_the_record_padded(design):
    K, c0, cs, As, f, ref = design
    res = FR.train_flow(FR.PaperLaw(LAW), K, f, c0, 1e-6, 20.0, LR_C, n_t=101, max_secs=0.0)
    assert res['status'] == 'timeout'
    assert res['W'].shape == (101, 21) and np.all(res['W'] == res['W'][0])


def test_no_wall_requires_a_time_scale(design):
    K, c0, cs, As, f, ref = design
    with pytest.raises(ValueError):
        FR.train_flow(FR.Linear(1.0), K, f, c0, 1e-6, 1.0, LR_C)


def test_random_masses_change_the_start_only_through_the_masses():
    rng = np.random.default_rng(5)
    a = FR.Cubic3(1.0).init(6, 1e-6, np.random.default_rng(5), 0.5)
    b = 1e-6 * np.exp(0.5 * rng.standard_normal(6))
    assert np.allclose(a[:, 0], np.sqrt(b), rtol=1e-14, atol=0) and np.array_equal(a[:, 0], a[:, 1])
    z = FR.PaperLaw(LAW).init(6, 1e-6, np.random.default_rng(5), 0.5)
    assert np.all(z[:, 0] == 0.0) and np.allclose(z[:, 1], np.sqrt(b), rtol=1e-14, atol=0)
    za = FR.PaperLaw(LAW, amp0=1e-2).init(6, 1e-6, np.random.default_rng(5))
    assert np.all(za[:, 0] != 0.0) and np.allclose(za[:, 1], 1e-3, rtol=1e-14, atol=0)

def test_formed_threshold_follows_the_smallest_feature_on_wide_ranges(design):
    """thr = min(0.01·max|A*|, 0.1·min|A*|): Paper 17's convention where the amplitudes are within a factor 10, lower where a feature is
    10³ times smaller than the largest (otherwise the path's own unit for it would not count as formed)."""
    K, c0, cs, As, f, ref = design
    taus, d = np.linspace(0, 1, 3), 21
    C = np.tile(c0[None, :, None], (3, 1, 1)).astype(float)
    for j, c in zip([5, 11, 16], cs):
        C[:, j, 0] = c

    def record(w5, w11, w16):
        W = np.zeros((3, d))
        W[-1, 5], W[-1, 11], W[-1, 16] = w5, w11, w16
        return dict(taus=taus, W=W, C=C)

    # amplitudes 1, 10⁻³, −0.5 (a range of 10³): thr = min(0.01, 10⁻⁴) = 10⁻⁴, so the small feature's unit (10⁻³) has formed
    assert FR.analyse(record(1.0, 1e-3, -0.5), K, f, c0, cs, np.array([1.0, 1e-3, -0.5]))['n_formed'] == 3
    # amplitudes 1, 0.5, −0.5: thr = min(0.01, 0.05) = 0.01 (the paper's), so a unit at 5·10⁻³ has not formed
    assert FR.analyse(record(1.0, 5e-3, -0.5), K, f, c0, cs, np.array([1.0, 0.5, -0.5]))['n_formed'] == 2
    assert FR.analyse(record(1.0, 2e-2, -0.5), K, f, c0, cs, np.array([1.0, 0.5, -0.5]))['n_formed'] == 3


def test_a_piece_over_the_lsoda_budget_is_redone_by_radau(design, monkeypatch):
    """`STIFF_BUDGET`: LSODA's stiffness detection can fail on long horizons (`rich_transit.relax_chunk`); the piece is then redone by Radau and
    the record agrees with the unconstrained one."""
    K, c0, cs, As, f, ref = design
    ro = FR.Cubic3(1.0)
    base = FR.train_flow(ro, K, f, c0, 1e-4, 12.0, LR_C, n_t=121, n_seg=4)
    monkeypatch.setattr(FR, 'STIFF_BUDGET', 30)                      # every LSODA piece exceeds it
    redo = FR.train_flow(ro, K, f, c0, 1e-4, 12.0, LR_C, n_t=121, n_seg=4)
    assert base['status'] == 'ok' and redo['status'] == 'ok'
    assert base['radau_pieces'] == 0 and base['aborted_nfev'] == 0 and redo['radau_pieces'] >= 1 and redo['aborted_nfev'] >= 30
    assert np.abs(redo['W'] - base['W']).max() < 1e-4 * np.abs(base['W']).max()
    assert np.abs(redo['C'] - base['C']).max() < 1e-4 * SIG


@pytest.mark.parametrize('consts', [(12.0, 1.0, 0.1, 1.0), (5.0, 2.3, 0.37, 0.7), (0.8, 0.45, 3.1, 1.9)])
def test_the_paper_law_is_the_euclidean_flow_for_any_constants(consts, design):
    """The identity for constants with G ≠ 1 and a_max ≠ 1 (the constants of the experiments are G = a_max = 1, which hides a dropped G or a_max), and a
    control: a wrong centre rate or a wrong ν rate breaks it, so the test is not vacuous."""
    m, G, gamma, a_max = consts
    L = FreeLaw(m=m, G=G, gamma=gamma, a_max=a_max)
    K, c0, cs, As, f, ref = design
    repo = RT.moving_rhs(L, RT.Bumps(K.xs[:, 0], K.sig), f, 'paradigm')
    ro = FR.PaperLaw(L)
    assert np.allclose(ro.lr, [G * m, gamma * m / 4], rtol=1e-14, atol=0)
    rng = np.random.default_rng(11)
    d = 7
    al, beta, c = rng.normal(0, 0.4, d), np.exp(rng.normal(-1, 0.8, d)), np.sort(rng.uniform(-0.9, 0.9, d))
    z = np.r_[np.stack([al, np.sqrt(beta)], 1).ravel(), c]
    theirs = repo(0.0, np.r_[al, np.log(beta), c])

    def mine(lr_c, lr_nu_scale=1.0):
        r2 = FR.PaperLaw(L)
        r2.lr = r2.lr * np.array([1.0, lr_nu_scale])
        return FR.make_rhs(r2, K, f, lr_c, d, 1)(0.0, z)
    ok = mine(G * m)
    assert np.allclose(ok[:2 * d:2], theirs[:d], rtol=1e-12, atol=0)
    assert np.allclose(2 * ok[1:2 * d:2] / np.sqrt(beta), theirs[d:2 * d], rtol=1e-12, atol=0)
    assert np.allclose(ok[2 * d:], theirs[2 * d:], rtol=1e-12, atol=0)
    assert not np.allclose(mine(G * m / 2)[2 * d:], theirs[2 * d:], rtol=1e-3)                     # a wrong centre rate
    assert not np.allclose(2 * mine(G * m, 4.0)[1:2 * d:2] / np.sqrt(beta), theirs[d:2 * d], rtol=1e-3)   # a wrong ν rate


def test_power_readouts_walls_and_depth_three_is_cubic3():
    eps, lr = 1e-6, 0.8
    assert np.isclose(FR.PowerN(3, lr).scale(eps), FR.Cubic3(lr).scale(eps), rtol=1e-14, atol=0)
    assert np.isclose(FR.PowerN(4, lr).scale(eps), 1 / (8 * lr * eps), rtol=1e-14, atol=0)              # η = s^{−(N−2)}/(N(N−2) lr)
    assert np.isclose(FR.PowerN(5, lr).scale(eps), eps ** -1.5 / (15 * lr), rtol=1e-14, atol=0)
    p = np.array([[0.3, 0.2], [0.5, 0.1]])
    assert np.allclose(FR.PowerN(3, lr).w(p), FR.Cubic3(lr).w(p), rtol=1e-14, atol=0)
    assert FR.PowerN(4, 1.0).name == 'D4' and FR.PowerN(5, 1.0).name == 'D5'


def _synthetic(design, W_rows):
    K, c0, cs, As, f, ref = design
    taus = np.linspace(0, 1, len(W_rows))
    C = np.tile(c0[None, :, None], (len(W_rows), 1, 1)).astype(float)
    for j, c in zip([5, 11, 16], cs):
        C[-1, j, 0] = c
    return dict(taus=taus, W=np.array(W_rows), C=C)


def test_analyse_semantics_tau_rel_is_a_maximum_order_ignores_spurious_units_dormant_peak_is_over_time(design):
    K, c0, cs, As, f, ref = design
    d = 21
    rows = [np.zeros(d) for _ in range(5)]
    for r, (j, a) in zip((2, 3, 4), ((5, As[0]), (16, As[2]), (11, As[1]))):          # entries at τ = 0.5, 0.75, 1 (unit 5, 16, 11)
        for rr in range(r, 5):
            rows[rr][j] = a
    rows[1][3], rows[2][3], rows[3][3], rows[4][3] = 0.05, 0.05, 1e-5, 1e-5            # a dormant unit that peaks mid-record, tiny at the end
    rows[1][9] = 0.5                                                                   # a spurious unit that forms FIRST (τ = 0.25) and decays
    for rr in range(2, 5):
        rows[rr][9] = 0.0
    res = _synthetic(design, rows)
    thr = FR.THETA * np.abs(As).max()
    ft = FR.formation_times(res['taus'], res['W'], thr)
    rt = np.array([0.5 * 1.0, 0.75 * 1.2, 1.0 * 0.9])                                   # reference times 0%, 20% and 10% off for units 5, 16, 11
    a = FR.analyse(res, K, f, c0, cs, As, ref=dict(order=[5, 16, 11], tau=list(rt), set=[5, 11, 16]))
    rel = np.abs(ft[[5, 16, 11]] - rt) / rt
    assert np.isclose(a['tau_rel'], rel.max(), rtol=1e-12, atol=0) and not np.isclose(a['tau_rel'], rel.mean(), rtol=1e-3)   # max, not mean
    assert a['order_equal']                                                            # the spurious unit 9 (formed first) is ignored
    assert not FR.analyse(res, K, f, c0, cs, As, ref=dict(order=[16, 5, 11], tau=list(rt), set=[5, 11, 16]))['order_equal']
    assert np.isclose(a['dormant_wmax'], 0.5, rtol=1e-12, atol=0)                      # the peak over time (unit 9's 0.5), not the end value
    assert a['n_formed'] == 3 and a['set_equal']


def test_threshold_free_gap_and_top_k_exactness(design):
    K, c0, cs, As, f, ref = design
    d = 21
    rows = [np.zeros(d) for _ in range(3)]
    rows[-1][5], rows[-1][11], rows[-1][16], rows[-1][3] = As[0], As[1], As[2], 2e-4
    res = _synthetic(design, rows)
    a = FR.analyse(res, K, f, c0, cs, As)
    assert np.isclose(a['gap_k'], np.log10(0.6 / 2e-4), rtol=1e-12, atol=0) and a['exact_topk']
    res['C'][-1, 11, 0] += 0.5 * SIG                                                    # a top unit 0.5σ off its feature
    assert not FR.analyse(res, K, f, c0, cs, As)['exact_topk']
    # D2-like: the K largest units are not on the features' amplitudes although the count at a high threshold could be right
    rows[-1][5] = 0.6 * As[0]
    assert not FR.analyse(_synthetic(design, rows), K, f, c0, cs, As)['exact_topk']


def test_adam_first_step_with_the_default_epsilon_uses_the_gradient_of_the_fit(design):
    """P0 from α = 0: only α has a gradient at step 1; with Adam's default ε_adam = 1e-8 (the size of ν²g/m here) the step is lr·d/(|d| + 1e-8), d = ν²a′g/m
    the *unit-rate* direction (the gradient of F), so ΔW = (ε/m)·Δα. A direction scaled by the group rates (Gm = 12) changes it."""
    K, c0, cs, As, f, ref = design
    ro = FR.PaperLaw(LAW)
    eps, lr_adam = 1e-6, 1e-9
    res = FR.train_adam(ro, K, f, c0, eps, n_steps=1, lr_adam=lr_adam, n_t=2)
    g = K.corr(c0[:, None], f)[0]
    dirn = g * eps / LAW.m
    expect = (eps / LAW.m) * lr_adam * dirn / (np.abs(dirn) + 1e-8)
    assert np.allclose(res['W'][-1], expect, rtol=1e-6, atol=0)


def test_a_failed_solve_with_an_empty_result_is_a_recorded_failure(design, monkeypatch):
    from types import SimpleNamespace
    K, c0, cs, As, f, ref = design
    monkeypatch.setattr(FR, 'solve_ivp', lambda *a, **k: SimpleNamespace(status=-1, y=[], message='Unexpected istate', nfev=7))
    res = FR.train_flow(FR.Cubic3(1.0), K, f, c0, 1e-4, 5.0, LR_C, n_t=41, n_seg=4)
    assert res['status'].startswith('failed') and res['W'].shape == (41, 21) and np.all(res['W'] == 0.0)
    assert res['radau_pieces'] == 1 and res['aborted_nfev'] == 7                      # LSODA abandoned, Radau failed too


def test_the_default_budget_and_the_counting_wrapper():
    assert FR.STIFF_BUDGET == 150000
    g = FR._counted(lambda t, y: y, 3)
    g(0, 1), g(0, 1), g(0, 1)
    with pytest.raises(FR._OverBudget):
        g(0, 1)



# -- branch pairs (Paper 17, Proposition 'starts on branches'): two units per atom, α(0) = ±a₀ --------------------------------------
def test_branch_pairs_start_and_logarithmic_time_scale():
    ro = FR.PaperPairs(LAW, 0.3)
    p = ro.init(42, 1e-6)
    assert p.shape == (42, 2) and np.all(p[:21, 0] == 0.3) and np.all(p[21:, 0] == -0.3) and np.allclose(p[:, 1], 1e-3, rtol=1e-14, atol=0)
    assert np.array_equal(ro.atoms(42), np.r_[np.arange(21), np.arange(21)])             # unit j and unit j + 21 are the pair of atom j
    assert np.isclose(ro.scale(1e-6), np.log(1e6) / (0.1 * np.tanh(0.3)), rtol=1e-14, atol=0)   # L_ε ρ, ρ = 1/(γ|a(0)|), a(0) = a_max tanh(a₀/a_max)
    assert np.isclose(ro.scale(1e-12), 2 * ro.scale(1e-6), rtol=1e-14, atol=0)           # logarithmic in ε
    assert np.isclose(FR.PaperPairs(LAW, 0.1).scale(1e-6) / ro.scale(1e-6), np.tanh(0.3) / np.tanh(0.1), rtol=1e-14, atol=0)
    p0 = FR.PaperLaw(LAW)
    assert np.array_equal(ro.lr, p0.lr) and np.array_equal(ro.w(p), p0.w(p)) and np.array_equal(ro.dw(p), p0.dw(p))     # P0's unit and rates
    assert ro.a0 == 0.3 and ro.amp0 == 0.0 and ro.name == 'P0p'                           # (no random amplitudes: the start is the pair's, not `amp0`'s)
    assert np.allclose(ro.w(p[:21]), -ro.w(p[21:]), rtol=1e-14, atol=0)                  # the pair's coefficients start opposite (and tiny)
    with pytest.raises(AssertionError):
        FR.PaperPairs(LAW, 0.0)
    with pytest.raises(AssertionError):
        ro.init(21, 1e-6)
    with pytest.raises(AssertionError):
        FR.PaperPairs.atoms(21)
    pr = ro.init(42, 1e-6, np.random.default_rng(5), 0.5)                                # random masses: the paper's pair has equal ones (rand = 0)
    assert np.allclose(pr[:, 1], np.sqrt(1e-6 * np.exp(0.5 * np.random.default_rng(5).standard_normal(42))), rtol=1e-14, atol=0)


@pytest.mark.parametrize('a0,eps', [(0.3, 1e-6), (0.3, 1e-12), (0.1, 1e-6)])
def test_form_ratio_is_the_formation_time_of_one_kernel_in_units_of_the_wall(a0, eps):
    """`form_ratio`: for a kernel on a branch under a constant drive g = 1 (η = t) the time at which |w| = θ over the wall's time L_ε/(γ|a(0)|),
    from dη = dβ/(γaβ), compared with the law integrated directly in (α, ln β)."""
    ro = FR.PaperPairs(LAW, a0)
    thr = 0.009

    def rhs(t, y):
        al, lb = y
        a = LAW.A * np.tanh(al / LAW.A)
        return [LAW.G * np.exp(lb) * (1 - (a / LAW.A) ** 2), LAW.gamma * a]

    def hit(t, y):
        return np.exp(y[1]) * LAW.A * np.tanh(y[0] / LAW.A) / LAW.m - thr
    hit.terminal = True
    from scipy.integrate import solve_ivp
    sol = solve_ivp(rhs, (0, 3 * ro.scale(eps)), [a0, np.log(eps)], method='LSODA', rtol=1e-11, atol=1e-13, events=hit)
    assert sol.status == 1
    assert np.isclose(sol.t_events[0][0] / ro.scale(eps), ro.form_ratio(eps, thr), rtol=1e-6, atol=0)
    assert 0 < ro.form_ratio(eps, thr) < 1


def test_form_ratio_tends_to_one_like_one_over_the_log_scale():
    ro = FR.PaperPairs(LAW, 0.3)
    c = [(1 - ro.form_ratio(e, 0.009)) * np.log(1 / e) for e in (1e-6, 1e-12, 1e-24, 1e-48)]
    assert max(c) / min(c) < 1.01 and 3.9 < c[0] < 4.2                                    # c(a₀) = 4.05 for a₀ = 0.3, θ = 0.009
    assert FR.PaperPairs(LAW, 0.1).form_ratio(1e-6, 0.009) < ro.form_ratio(1e-6, 0.009)  # a smaller a₀ is earlier by more
    assert ro.form_ratio(1e-6, 0.001) < ro.form_ratio(1e-6, 0.009) < ro.form_ratio(1e-6, 0.09)   # a lower threshold is crossed earlier


def test_formed_threshold_is_the_one_of_analyse():
    assert np.isclose(FR.formed_threshold([0.9, -0.6, 0.75]), 0.009, rtol=1e-14, atol=0)                    # 0.01 max|A*|
    assert np.isclose(FR.formed_threshold([1.0, 1e-3, -0.5]), 1e-4, rtol=1e-14, atol=0)                     # 0.1 min|A*| binds
    assert np.isclose(FR.formed_threshold([1.0, 0.5], thr_rel=0.3), 0.05, rtol=1e-14, atol=0)               # thr_rel 0.3: 0.3 but 0.1·0.5 binds
    assert np.isclose(FR.formed_threshold([1.0, -0.5], thr_rel=0.04), 0.04, rtol=1e-14, atol=0)


def _pair_record(design, entries, extra=()):
    """A synthetic record of 42 units (the 21 atoms twice, T = 5 points of τ): entries = [(unit, index of τ at which it reaches its coefficient, coefficient)]."""
    K, c0, cs, As, f, ref = design
    T = 5
    W = np.zeros((T, 42))
    for u, k, w in list(entries) + list(extra):
        W[k:, u] = w
    C = np.tile(np.r_[c0, c0][None, :, None], (T, 1, 1)).astype(float)
    return dict(taus=np.linspace(0, 1, T), W=W, C=C)


def _atom_times(res, As):
    """The formation times of the three atoms 5, 16, 11 of the design in a record of `_pair_record` (an atom's is its earlier unit's)."""
    ft = FR.formation_times(res['taus'], res['W'], FR.formed_threshold(As))
    return {a: np.nanmin([ft[a], ft[a + 21]]) for a in (5, 16, 11)}


def test_atom_summary_forms_an_atom_with_either_unit_of_its_pair(design):
    K, c0, cs, As, f, ref = design
    atom_of = FR.PaperPairs.atoms(42)
    # atom 5 by its + unit (5), atom 16 by its - unit (16 + 21 = 37), atom 11 by its - unit (32): at the 3rd, 4th and 5th point of the record
    res = _pair_record(design, [(5, 2, As[0]), (37, 3, -As[2]), (32, 4, As[1])])
    t = _atom_times(res, As)
    path = dict(order=[5, 16, 11], tau=[t[5], t[16], t[11]], set=[5, 11, 16])
    b = FR.atom_summary(res, As, atom_of, ref=path)
    assert b['n_atoms'] == 3 and b['atoms'] == [5, 11, 16] and b['one_each'] and b['set_equal'] and b['order_equal']
    assert b['form_order'] == [5, 16, 11] and t[5] < t[16] < t[11]
    assert b['tau_rel'] == 0.0 and b['tau_ratio'] == 0.0
    assert set(FR.atom_summary(res, As, atom_of)) == {'n_atoms', 'atoms', 'one_each', 'form_order'}            # without a reference: the counts only
    # both units of atom 5 above the threshold: still one atom, no longer one unit each
    two = FR.atom_summary(_pair_record(design, [(5, 2, As[0]), (26, 3, -0.2), (37, 3, -As[2]), (32, 4, As[1])]), As, atom_of, ref=path)
    assert two['n_atoms'] == 3 and not two['one_each'] and two['set_equal']
    # a fourth atom (a spurious one: unit 30 = atom 9) forming FIRST: the set differs, the path's atoms keep their order (the spurious one is ignored)
    sp = FR.atom_summary(_pair_record(design, [(5, 2, As[0]), (37, 3, -As[2]), (32, 4, As[1])], extra=[(30, 1, 0.5)]), As, atom_of, ref=path)
    assert sp['n_atoms'] == 4 and not sp['set_equal'] and sp['order_equal'] and sp['form_order'][0] == 9
    # an atom of the path's set that never forms: no times, the set differs
    miss = FR.atom_summary(_pair_record(design, [(5, 2, As[0]), (37, 3, -As[2])]), As, atom_of, ref=path)
    assert miss['n_atoms'] == 2 and not miss['set_equal'] and np.isnan(miss['tau_rel']) and np.isnan(miss['tau_ratio'])
    # the wrong order: atom 16 before atom 5
    wr = _pair_record(design, [(5, 3, As[0]), (37, 2, -As[2]), (32, 4, As[1])])
    w = FR.atom_summary(wr, As, atom_of, ref=path)
    assert w['set_equal'] and not w['order_equal'] and w['form_order'] == [16, 5, 11]
    # |w| must EXCEED the formed threshold: a unit exactly at it has not formed
    thr = FR.formed_threshold(As)
    at = FR.atom_summary(_pair_record(design, [(5, 2, As[0]), (37, 3, -As[2]), (32, 4, thr)]), As, atom_of, ref=path)
    assert at['n_atoms'] == 2 and at['atoms'] == [5, 16] and not at['set_equal']
    assert at['form_order'] == [5, 16, 11]                                  # (it did reach it: the formation time is the first ≥, as in `formation_times`)


def test_atom_summary_times_an_atom_by_its_earlier_unit_and_tau_ratio_ignores_a_common_factor(design):
    K, c0, cs, As, f, ref = design
    atom_of = FR.PaperPairs.atoms(42)
    base = _pair_record(design, [(5, 2, As[0]), (37, 3, -As[2]), (32, 4, As[1])])
    t = _atom_times(base, As)
    # the pair of atom 5 forms at the 2nd point (unit 26) and again at the 3rd (unit 5): the earlier counts, so atom 5 is early by the gap between them
    both = _pair_record(design, [(26, 1, -As[0]), (5, 2, As[0]), (37, 3, -As[2]), (32, 4, As[1])])
    tb = _atom_times(both, As)
    assert tb[5] < t[5] and tb[16] == t[16] and tb[11] == t[11]
    b = FR.atom_summary(both, As, atom_of, ref=dict(order=[5, 16, 11], tau=[t[5], t[16], t[11]], set=[5, 11, 16]))
    assert np.isclose(b['tau_rel'], (t[5] - tb[5]) / t[5], rtol=1e-12, atol=0) and b['tau_rel'] > 0.1
    assert b['form_order'] == [5, 16, 11]
    # every entry early by the factor 1/2 (a common factor): the maximum of the relative errors is 1/2, the ratios are exact
    c = FR.atom_summary(base, As, atom_of, ref=dict(order=[5, 16, 11], tau=[2 * t[5], 2 * t[16], 2 * t[11]], set=[5, 11, 16]))
    assert np.isclose(c['tau_rel'], 0.5, rtol=1e-12, atol=0) and np.isclose(c['tau_ratio'], 0.0, atol=1e-12, rtol=0)
    # atom 16 off by 20% (the path's time 1.25 times ours): a maximum, not a mean; the ratio error (to the first entrant 5) is that atom's
    d = FR.atom_summary(base, As, atom_of, ref=dict(order=[5, 16, 11], tau=[t[5], t[16] / 0.8, t[11]], set=[5, 11, 16]))
    assert np.isclose(d['tau_rel'], 0.2, rtol=1e-12, atol=0)
    assert np.isclose(d['tau_ratio'], max(abs((t[16] / t[5]) / ((t[16] / 0.8) / t[5]) - 1), abs((t[11] / t[5]) / (t[11] / t[5]) - 1)), rtol=1e-12, atol=0)
    # the first entrant of the path (the smallest τ) is the reference of the ratios: atom 5 late by 25% makes the others' ratios off by 25%
    e = FR.atom_summary(base, As, atom_of, ref=dict(order=[5, 16, 11], tau=[t[5] / 0.8, t[16], t[11]], set=[5, 11, 16]))
    assert np.isclose(e['tau_rel'], 0.2, rtol=1e-12, atol=0)
    assert np.isclose(e['tau_ratio'], 0.25, rtol=1e-12, atol=0)
    # the formed threshold of `analyse` is used: with amplitudes 1, 0.5, -0.5 a unit at 5e-3 has not formed (theta = 0.01), a unit at 2e-2 has
    small = FR.atom_summary(_pair_record(design, [(5, 2, 1.0), (37, 3, -0.5), (32, 4, 5e-3)]), np.array([1.0, 0.5, -0.5]), atom_of)
    assert small['n_atoms'] == 2 and small['atoms'] == [5, 16]
    big = FR.atom_summary(_pair_record(design, [(5, 2, 1.0), (37, 3, -0.5), (32, 4, 2e-2)]), np.array([1.0, 0.5, -0.5]), atom_of)
    assert big['n_atoms'] == 3
    # the identity grouping (one unit per atom) agrees with `analyse` on the formed set and the order
    res = _pair_record(design, [(5, 2, As[0]), (16, 3, As[2]), (11, 4, As[1])])
    res = dict(taus=res['taus'], W=res['W'][:, :21], C=res['C'][:, :21])
    a = FR.analyse(res, K, f, c0, cs, As, ref=dict(order=[5, 16, 11], tau=[0.5, 0.75, 1.0], set=[5, 11, 16]))
    g = FR.atom_summary(res, As, np.arange(21), ref=dict(order=[5, 16, 11], tau=[0.5, 0.75, 1.0], set=[5, 11, 16]))
    assert g['atoms'] == [5, 11, 16] and g['n_atoms'] == a['n_formed'] and g['form_order'] == a['form_order'] and g['set_equal'] == a['set_equal']
    assert g['order_equal'] == a['order_equal'] and g['tau_rel'] == a['tau_rel']


def _run_pairs(design, a0, eps, tau_end=20.0):
    K, c0, cs, As, f, ref = design
    ro = FR.PaperPairs(LAW, a0)
    c0u = np.r_[c0, c0]
    res = FR.train_flow(ro, K, f, c0u, eps, tau_end, LR_C, n_t=801)
    return FR.analyse(res, K, f, c0u, cs, As), FR.atom_summary(res, As, ro.atoms(42), ref=ref), res


def test_branch_pairs_on_the_design_a0_decides_like_the_depth(design):
    """The referee's validation cells: a₀ = 0.3 at ε = 10⁻⁶ forms 6 units (depth-2-like, no gap), at 10⁻¹² three; a₀ = 0.1 forms three at 10⁻⁶ (gap
    about 1.9 decades, against 5.3 for D3). The sign of the amplitude picks the unit of the pair; on this separated target the entries are the path's up to a common factor (not so on coupled
    targets: part K, `test_the_first_entrant_is_early_by_the_derived_factor...`)."""
    a, b, res = _run_pairs(design, 0.3, 1e-6)
    assert res['status'] == 'ok' and a['n_formed'] == 6 and b['n_atoms'] == 6 and not b['set_equal'] and not a['exact']
    assert a['gap_k'] < 1.0 and b['order_equal'] and b['atoms'] == [4, 5, 10, 11, 16, 17]
    a, b, res = _run_pairs(design, 0.3, 1e-12)
    assert a['n_formed'] == 3 and b['n_atoms'] == 3 and b['set_equal'] and b['order_equal'] and b['one_each'] and a['exact']
    assert 2.0 < a['gap_k'] < 2.4
    assert np.isclose(b['tau_rel'], 1 - FR.PaperPairs(LAW, 0.3).form_ratio(1e-12, FR.formed_threshold(design[3])), atol=5e-3, rtol=0)   # the derived common factor
    assert b['tau_ratio'] < 5e-3
    a, b, res = _run_pairs(design, 0.1, 1e-6)
    assert a['n_formed'] == 3 and b['set_equal'] and a['exact'] and 1.7 < a['gap_k'] < 2.1
    thr = FR.formed_threshold(design[3])
    assert sorted(np.flatnonzero(np.abs(res['W'][-1]) > thr).tolist()) == [5, 16, 32]       # the +0.9 and +0.75 features by + units, the −0.6 one by a − unit
    assert b['tau_ratio'] < 5e-3 and abs(b['tau_rel'] - (1 - FR.PaperPairs(LAW, 0.1).form_ratio(1e-6, thr))) < 5e-3


# -- the survivors of round 156's mutation run: solver bookkeeping, descent steps, threshold edges, tolerances, the rms -------------
def _fake_solver(script):
    """A stand-in for `solve_ivp`: script(call index, method, y0, t_eval) -> (status, list of state columns, nfev)."""
    from types import SimpleNamespace
    calls = []

    def fake(fun, span, y0, method=None, t_eval=None, **kw):
        status, cols, nfev = script(len(calls), method, np.asarray(y0, float), t_eval)
        calls.append(method)
        return SimpleNamespace(status=status, y=np.array(cols).T if len(cols) else [], message='scripted', nfev=nfev)
    fake.calls = calls
    return fake


def test_the_evaluations_of_accepted_and_abandoned_solves_are_counted_separately(design, monkeypatch):
    """LSODA fails on every piece (7 evaluations, no result), Radau accepts (5): per piece one Radau piece, 7 abandoned and 5 accepted evaluations."""
    K, c0, cs, As, f, ref = design
    fake = _fake_solver(lambda i, m, y0, te: (-1, [], 7) if m == 'LSODA' else (0, [y0] * len(te), 5))
    monkeypatch.setattr(FR, 'solve_ivp', fake)
    res = FR.train_flow(FR.Cubic3(1.0), K, f, c0, 1e-4, 5.0, LR_C, n_t=41, n_seg=4)
    assert res['status'] == 'ok' and fake.calls == ['LSODA', 'Radau'] * 4
    assert res['radau_pieces'] == 4 and res['aborted_nfev'] == 28 and res['nfev'] == 20


def test_a_failure_after_a_good_piece_pads_the_record_with_the_last_state_and_counts_the_failed_solve(design, monkeypatch):
    K, c0, cs, As, f, ref = design
    d1, d2 = np.zeros(21 * 2 + 21), np.zeros(21 * 2 + 21)
    d1[0], d2[0] = 0.5, 0.7                                                    # unit 0's u moves, so that its coefficient u³ − v³ is not zero

    def script(i, m, y0, te):
        if i == 0:
            return 0, [y0 + d1] * len(te), 3                                   # the first piece is accepted
        return -1, [y0 + d2] * 3, (7 if m == 'LSODA' else 4)                   # the second: LSODA and then Radau return three columns and fail
    fake = _fake_solver(script)
    monkeypatch.setattr(FR, 'solve_ivp', fake)
    res = FR.train_flow(FR.Cubic3(1.0), K, f, c0, 1e-4, 5.0, LR_C, n_t=41, n_seg=4)
    assert res['status'].startswith('failed') and fake.calls == ['LSODA', 'LSODA', 'Radau']
    s = np.sqrt(1e-4)
    w1, w2 = (s + 0.5) ** 3 - s ** 3, (s + 0.5 + 0.7) ** 3 - s ** 3                  # (the second piece starts from the first one's last state)
    W0 = res['W'][:, 0]
    assert W0[0] == 0.0 and np.allclose(W0[1:11], w1, rtol=1e-12, atol=0)       # the accepted piece
    assert np.allclose(W0[11:], w2, rtol=1e-12, atol=0)                         # the failed piece's three columns, then padded with the last of them
    assert res['radau_pieces'] == 1 and res['aborted_nfev'] == 7 and res['nfev'] == 3 + 4     # accepted piece + the failed Radau solve


def test_a_timeout_pads_the_record_with_the_last_state_not_with_zeros(design, monkeypatch):
    K, c0, cs, As, f, ref = design
    clock = iter([0.0, 0.0, 100.0, 100.0, 100.0])                              # t0, the check before piece 1 (passes), before piece 2 (times out), secs
    monkeypatch.setattr(FR, 'time', type('T', (), {'time': staticmethod(lambda: next(clock))}))
    res = FR.train_flow(FR.Cubic3(1.0), K, f, c0, 1e-4, 5.0, LR_C, n_t=41, n_seg=4, max_secs=1.0)
    assert res['status'] == 'timeout'
    assert np.abs(res['W'][10]).max() > 0 and np.all(res['W'][11:] == res['W'][10]) and np.all(res['C'][11:] == res['C'][10])
    assert np.abs(res['C'][10] - c0[:, None]).max() < 1e-3 and np.all(res['C'][11:] != 0.0)       # the centres are the last state's, not zero


def test_gd_takes_the_ceiling_of_the_steps(design):
    """⌈τ_end·scale/h⌉ steps: D3 at ε = 10⁻⁴ has scale 100/3, so τ_end = 0.1 and h = 0.02 give 166.67 → 167."""
    K, c0, cs, As, f, ref = design
    gd = FR.train_gd(FR.Cubic3(1.0), K, f, c0, 1e-4, 0.1, LR_C, h=0.02, n_t=11)
    assert gd['steps'] == 167 and gd['nfev'] == 167 and gd['status'] == 'ok'
    assert np.isclose(gd['taus'][-1], 0.1, rtol=1e-14, atol=0)


def test_formation_time_of_a_unit_that_only_touches_the_threshold():
    taus = np.linspace(0, 1, 5)
    assert FR.formation_times(taus, np.full((5, 1), 1e-2), thr=1e-2)[0] == 0.0         # touching it counts as reaching it (≥)
    W = np.zeros((5, 1))
    W[2:, 0] = 1e-2
    assert np.isclose(FR.formation_times(taus, W, thr=1e-2)[0], taus[2], rtol=1e-14, atol=0)


def test_match_features_amplitude_error_is_relative():
    ok, dist, aerr = FR.match_features(np.array([[0.0]]), np.array([22.0]), [0.0], np.array([20.0]))
    assert np.isclose(aerr[0], 0.1, rtol=1e-14, atol=0)                                # |22 − 20|/20, not 2
    ok, dist, aerr = FR.match_features(np.array([[0.0]]), np.array([0.012]), [0.0], np.array([0.01]))
    assert np.isclose(aerr[0], 0.2, rtol=1e-12, atol=0)                                # a small feature 20% off is 20% off


def test_exactness_needs_the_amplitude_within_five_percent_and_the_unit_strictly_above_the_threshold(design):
    K, c0, cs, As, f, ref = design
    taus, d = np.linspace(0, 1, 3), 21
    C = np.tile(c0[None, :, None], (3, 1, 1)).astype(float)
    for j, c in zip([5, 11, 16], cs):
        C[:, j, 0] = c

    def record(w5, w11, w16):
        W = np.zeros((3, d))
        W[-1, 5], W[-1, 11], W[-1, 16] = w5, w11, w16
        return dict(taus=taus, W=W, C=C)
    good = FR.analyse(record(As[0], As[1], As[2]), K, f, c0, cs, As)
    assert good['exact'] and good['pos_err'] == 0.0 and good['amp_err'] == 0.0
    ten = FR.analyse(record(As[0], 1.10 * As[1], As[2]), K, f, c0, cs, As)             # position exact, one amplitude 10% off: not exact
    assert not ten['exact'] and ten['n_formed'] == 3 and np.isclose(ten['amp_err'], 0.10, rtol=1e-9, atol=0) and ten['pos_err'] == 0.0
    four = FR.analyse(record(As[0], 1.04 * As[1], As[2]), K, f, c0, cs, As)             # 4% off: within the 5%
    assert four['exact'] and np.isclose(four['amp_err'], 0.04, rtol=1e-9, atol=0)
    # amplitudes 1, 0.5, −0.5: the threshold is 0.01, and a unit exactly at it has not formed (|w| > thr)
    A2 = np.array([1.0, 0.5, -0.5])
    assert FR.analyse(record(1.0, 1e-2, -0.5), K, f, c0, cs, A2)['n_formed'] == 2
    assert FR.analyse(record(1.0, np.nextafter(1e-2, 1.0), -0.5), K, f, c0, cs, A2)['n_formed'] == 3


def test_the_sweep_is_relative_to_the_largest_amplitude(design):
    K, c0, cs, As, f, ref = design
    taus, d = np.linspace(0, 1, 3), 21
    C = np.tile(c0[None, :, None], (3, 1, 1)).astype(float)
    for j, c in zip([5, 11, 16], cs):
        C[:, j, 0] = c
    W = np.zeros((3, d))
    W[-1, 5], W[-1, 11], W[-1, 16] = As[0], As[1], As[2]
    W[-1, 3], W[-1, 7] = 2e-4, 7.5e-4              # between 1e-4·0.9 and 1e-3·0.6 (the levels with the smallest amplitude would count 7.5e-4 at 1e-3)
    sw = FR.analyse(dict(taus=taus, W=W, C=C), K, f, c0, cs, As)['n_formed_sweep']
    assert sw == {1e-4: 5, 1e-3: 3, 1e-2: 3, 1e-1: 3}


def test_fit_value_is_half_the_mean_square_and_the_root_mean_square(design):
    K, c0, cs, As, f, ref = design
    w = np.zeros(21)
    w[5], w[11] = 0.5, -0.2
    F, rms = FR.fit_value(K, f, w, c0)
    r = f - K.synth(c0[:, None], w)
    assert np.isclose(F, 0.5 * np.mean(r ** 2), rtol=1e-14, atol=0) and np.isclose(rms, np.sqrt(np.mean(r ** 2)), rtol=1e-14, atol=0)
    assert not np.isclose(rms, np.mean(np.abs(r)), rtol=1e-3)                          # (the residual is not flat, so rms ≠ mean |r|)
    assert np.isclose(rms ** 2, 2 * F, rtol=1e-14, atol=0)


def test_atom_summary_covers_the_last_atom(design):
    K, c0, cs, As, f, ref = design
    atom_of = FR.PaperPairs.atoms(42)
    res = _pair_record(design, [(5, 2, As[0]), (41, 3, -As[2])])                    # unit 41 is the − unit of atom 20, the last atom
    b = FR.atom_summary(res, As, atom_of, ref=dict(order=[5, 20], tau=[0.5, 0.75], set=[5, 20]))
    assert b['atoms'] == [5, 20] and b['n_atoms'] == 2 and b['set_equal'] and b['order_equal'] and b['form_order'] == [5, 20]
    assert np.isfinite(b['tau_rel'])


# -- the PZ2 round: PowerN's start, the first entrant's earliness, the atoms' formation times ---------------------------------------------
def test_power_readouts_start_at_s_equal_sqrt_eps_like_cubic3():
    """`PowerN.init`: (u, v)(0) = (s, s) with s² = ε b_i, the start of `Cubic3` (D3's), for every N: D4 and D5 of part B start here, and no engine
    test had pinned it (a start at s = ε, a mutant, passed every test: the conserved-quantity test is blind to the start)."""
    eps = 1e-6
    for N in (3, 4, 5):
        p = FR.PowerN(N, 0.9).init(6, eps)
        assert p.shape == (6, 2) and np.allclose(p, 1e-3, rtol=1e-14, atol=0)                  # s = ε^(1/2), not ε
        assert np.array_equal(p, FR.Cubic3(1.1).init(6, eps))
        pr = FR.PowerN(N, 0.9).init(6, eps, np.random.default_rng(5), 0.5)                      # random masses: s² = ε b_i
        assert np.array_equal(pr, FR.Cubic3(1.1).init(6, eps, np.random.default_rng(5), 0.5))
        assert np.allclose(pr[:, 0], np.sqrt(eps * np.exp(0.5 * np.random.default_rng(5).standard_normal(6))), rtol=1e-14, atol=0)
        assert np.array_equal(pr[:, 0], pr[:, 1])
    assert FR.PowerN(4, 1.0).w(FR.PowerN(4, 1.0).init(3, 1e-6)).tolist() == [0.0, 0.0, 0.0]        # w(0) = 0: u^N − v^N at u = v


def test_atom_summary_earliness_is_one_minus_formation_over_entry_in_the_paths_entry_order(design):
    """`form_times`, `earliness`: the atoms of the path's set in the order of the path's entry times (not by index, not in the order of `order`'s
    list), 1 − t_j/τ_j for each, nan only for the atom that never formed."""
    K, c0, cs, As, f, ref = design
    atom_of = FR.PaperPairs.atoms(42)
    res = _pair_record(design, [(5, 2, As[0]), (37, 3, -As[2]), (32, 4, As[1])])
    t = _atom_times(res, As)                                       # the formation times of atoms 5, 16, 11 (t5 < t16 < t11)
    # the path enters 16, then 5, then 11, and our atoms are early by 20%, 50% and 10%
    path = dict(order=[16, 5, 11], tau=[t[16] / 0.8, t[5] / 0.5, t[11] / 0.9], set=[5, 11, 16])
    assert path['tau'][0] < path['tau'][1] < path['tau'][2]
    b = FR.atom_summary(res, As, atom_of, ref=path)
    assert np.allclose(b['form_times'], [t[16], t[5], t[11]], rtol=1e-14, atol=0)
    assert np.allclose(b['earliness'], [0.2, 0.5, 0.1], rtol=1e-12, atol=0)
    assert np.isclose(b['tau_rel'], 0.5, rtol=1e-12, atol=0)       # the largest earliness
    # a path listing its entries out of the order of its times: the atoms are still sorted by the path's times
    path2 = dict(order=[5, 11, 16], tau=[t[5] / 0.5, t[11] / 0.9, t[16] / 0.8], set=[5, 11, 16])
    b2 = FR.atom_summary(res, As, atom_of, ref=path2)
    assert path2['tau'][2] < path2['tau'][0] < path2['tau'][1]
    assert np.allclose(b2['form_times'], [t[16], t[5], t[11]], rtol=1e-14, atol=0) and np.allclose(b2['earliness'], [0.2, 0.5, 0.1], rtol=1e-12, atol=0)
    # an atom that never formed has nan; the others keep theirs
    miss = FR.atom_summary(_pair_record(design, [(5, 2, As[0]), (37, 3, -As[2])]), As, atom_of, ref=path)
    assert np.isnan(miss['tau_rel']) and np.isnan(miss['tau_ratio'])
    assert len(miss['earliness']) == 3 and np.isnan(miss['earliness'][2]) and np.isfinite(miss['earliness'][:2]).all()
    assert np.isnan(miss['form_times'][2]) and np.isfinite(miss['form_times'][:2]).all()
    # a late atom (formed after the path's entry) is negatively early
    late = dict(order=[5, 16, 11], tau=[t[5] / 0.5, t[16] / 0.8, t[11] * 0.5], set=[5, 11, 16])
    assert np.isclose(FR.atom_summary(res, As, atom_of, ref=late)['earliness'][0], -1.0, rtol=1e-12, atol=0)
    assert 'earliness' not in FR.atom_summary(res, As, atom_of)                                       # without a reference: the counts only


def test_the_first_entrant_is_early_by_the_derived_factor_a_later_atom_only_if_its_dual_is_straight(design):
    """`form_ratio` is the first entrant's: on the design (a₀ = 0.1, ε = 10⁻⁶) the first atom of the path's order is early by 1 − form_ratio, to 5·10⁻³.
    The later atoms of this separated target are early by nearly the same (their duals are nearly straight), which is what the part K sample shows;
    the coupled targets of part K (tests/test_parametrization.py) are where it fails."""
    a, b, res = _run_pairs(design, 0.1, 1e-6)
    thr = FR.formed_threshold(design[3])
    pred = 1 - FR.PaperPairs(LAW, 0.1).form_ratio(1e-6, thr)
    assert len(b['earliness']) == 3 and abs(b['earliness'][0] - pred) < 5e-3
    assert all(abs(x - pred) < 0.02 for x in b['earliness'])
