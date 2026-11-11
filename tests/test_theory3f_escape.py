"""Normalised optimizers and the two-layer small start (experiments/theory3f_escape.py): the optimizers, the theorems' closed forms, the checks against the simulators, the report."""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'experiments'))
sys.path.insert(0, os.path.join(HERE, '..'))
import theory3f_escape as T      # noqa: E402


# ---- the optimizers -------------------------------------------------------------------------------------------------------------
def test_adam_first_step_is_sign_and_matches_the_formula():
    g = np.array([3.0, -2.0, 0.0])
    p = T.Opt('adam', 0.1, eps_adam=0.0).step(np.zeros(3), g)
    assert np.allclose(p, [0.1, -0.1, 0.0])                                   # mhat/sqrt(vhat) = sign(g) at t = 1; the zero gradient does not move (0/0 := 0)
    p = T.Opt('sign', 0.1).step(np.zeros(3), g)
    assert np.allclose(p, [0.1, -0.1, 0.0])
    o = T.Opt('adam', 1.0, eps_adam=0.5)
    assert np.allclose(o.step(np.zeros(1), np.array([2.0])), [2.0 / 2.5])    # g/(|g| + eps_adam)


def test_adam_second_step_by_hand():
    o = T.Opt('adam', 1.0, beta1=0.9, beta2=0.999, eps_adam=0.0)
    o.step(np.zeros(1), np.array([1.0]))
    p = o.step(np.zeros(1), np.array([3.0]))
    m, v = 0.9 * 0.1 + 0.1 * 3.0, 0.999 * 0.001 + 0.001 * 9.0
    assert np.allclose(p, (m / (1 - 0.9 ** 2)) / np.sqrt(v / (1 - 0.999 ** 2)))


def test_gd_pre_and_ngd_and_shared_second_moment():
    g = np.array([3.0, 4.0])
    assert np.allclose(T.Opt('gd', 0.5, pre=np.array([2.0, 1.0])).step(np.zeros(2), g), [3.0, 2.0])
    assert np.allclose(T.Opt('ngd', 1.0).step(np.zeros(2), g), [0.6, 0.8])
    assert np.allclose(T.Opt('ngd', 1.0).step(np.zeros(2), np.zeros(2)), 0.0)
    p = T.Opt('adam_g', 1.0, eps_adam=0.0).step(np.zeros(2), g)                # one scalar v: the ratio of the steps is the ratio of the gradients
    assert np.allclose(p[0] / p[1], 3.0 / 4.0) and np.allclose(p, g / np.sqrt(np.mean(g ** 2)))


def test_opt_adam_reproduces_the_repo_adam():
    P, FR, L, K, c0, cs, As, f = T.design()
    ro = P.make_readout('D3')
    ref = FR.train_adam(ro, K, f, c0, 1e-6, 60, 1e-3, n_t=7)
    mine = T.train_opt(FR, ro, K, f, c0, 1e-6, 60, T.Opt('adam', 1e-3, eps_adam=1e-8), n_t=7)
    assert np.allclose(ref['W'], mine['W'], rtol=1e-10, atol=1e-18) and np.allclose(ref['C'], mine['C'], rtol=1e-10, atol=1e-18)
    sg = FR.train_adam(ro, K, f, c0, 1e-6, 60, 1e-3, n_t=7, betas=(0.0, 0.0), eps_adam=1e-40)
    mine = T.train_opt(FR, ro, K, f, c0, 1e-6, 60, T.Opt('sign', 1e-3), n_t=7)
    assert np.allclose(sg['W'], mine['W'], rtol=1e-10, atol=1e-18)


def test_ratio_bounds():
    c, C, tstar = T.adam_ratio_bounds(0.9, 0.999)
    assert abs(c - 0.1) < 1e-15 and abs(C - 7.27) < 0.01 and tstar > 1000
    assert abs(T.adam_ratio_bounds(0.0, 0.999)[1] - 1 / np.sqrt(1 - 0.999)) < 1e-9          # beta1 = 0: R_t = g/sqrt(vhat) <= 1/sqrt(1 - beta2)
    # the upper bound holds for random real sequences, the lower bound for positive nondecreasing ones
    rng = np.random.default_rng(1)
    for _ in range(30):
        G = rng.standard_normal(500) * rng.choice([1.0, 100.0], 500)
        m = v = 0.0
        for t in range(1, 501):
            m = 0.9 * m + 0.1 * G[t - 1]
            v = 0.999 * v + 0.001 * G[t - 1] ** 2
            R = (m / (1 - 0.9 ** t)) / np.sqrt(v / (1 - 0.999 ** t))
            assert abs(R) <= C + 1e-12
        G = np.cumsum(rng.exponential(1.0, 500)) ** 3
        m = v = 0.0
        for t in range(1, 501):
            m = 0.9 * m + 0.1 * G[t - 1]
            v = 0.999 * v + 0.001 * G[t - 1] ** 2
            assert (m / (1 - 0.9 ** t)) / np.sqrt(v / (1 - 0.999 ** t)) >= c - 1e-12


# ---- F1: sign descent, Adam ------------------------------------------------------------------------------------------------------
def test_gd_escape_closed_form_and_ordering():
    assert abs(T.gd_escape_time(1.0, 0.05, 0.5) - 6.0) < 1e-12                          # (20 - 2)/3
    assert abs(T.gd_escape_time(0.1, 0.05, 0.5) / T.gd_escape_time(1.0, 0.05, 0.5) - 10.0) < 1e-12
    u, n = 0.05, 0
    while u < 0.5:
        u += 0.001 * 3 * 0.1 * u * u
        n += 1
    assert abs(n * 0.001 - T.gd_escape_time(0.1, 0.05, 0.5)) < 0.3                      # explicit Euler of the flow


def test_sign_descent_theorem_equal_steps_for_every_drive_and_start():
    G = np.eye(5)
    g = np.array([1.0, 1e-1, 1e-2, 1e-3, 1e-4])
    thw, lr = 1e-6, 1e-4
    for s in (1.3e-3, 1e-5, 1e-9):
        W, Ts = T.run_units(T.Opt('sign', lr), G, g, s, 300, thw)
        assert np.all(Ts == Ts[0]) and Ts[0] == T.sign_escape_step(s, thw, lr)           # Theorem 1: all units at T_env, exactly
    # the envelope: no unit can be earlier; a unit whose hypothesis fails (|g| < theta) is not required to cross
    G2, c2 = np.eye(2), np.array([1.0, 1e-8])
    W, Ts = T.run_units(T.Opt('sign', lr), G2, c2, 1.3e-3, 300, thw)
    assert Ts[0] == T.sign_escape_step(1.3e-3, thw, lr) and (np.isnan(Ts[1]) or Ts[1] >= Ts[0])
    # independence of the start: the escape step changes by at most s^2/(a0 lr) + 1 between s and 0
    a0 = (thw / 2) ** (1 / 3)
    assert 0 <= T.sign_escape_step(0.0, thw, lr) - T.sign_escape_step(1.3e-3, thw, lr) <= 1.3e-3 ** 2 / (a0 * lr) + 1


def test_sign_descent_coupled_gram_still_simultaneous():
    rng = np.random.default_rng(3)
    A = rng.standard_normal((6, 6)) * 0.1
    G = np.eye(6) + A @ A.T * 0.2
    g = np.array([1.0, 0.5, 0.2, 0.1, -0.3, -0.05])
    thw, lr = 1e-6, 1e-4
    ok = np.abs(g) > thw * np.abs(G).sum(1)
    assert ok.all()
    W, Ts = T.run_units(T.Opt('sign', lr), G, G @ g, 1.3e-3, 300, thw)
    assert np.all(Ts == Ts[0])                                                         # both signs of the drive: the D3 unit is two-sided


def test_adam_scale_invariance_and_window():
    steps = [T.frozen_adam_run(1e-3, 0.1, 1e-4, 0.9, 0.999, 0.0, sc)[0] for sc in (1e-8, 1e-3, 1.0, 1e5)]
    assert len(set(steps)) == 1 and steps[0] > 0                                       # exact invariance for eps_adam = 0
    t_big = T.frozen_adam_run(1e-3, 0.1, 1e-4, 0.9, 0.999, 1e-8, 1e-8, max_steps=3000)[0]
    assert t_big == -1                                                                 # eps_adam = 1e-8 against gradients ~ 1e-14: the SGD regime, not reached
    c, C, _ = T.adam_ratio_bounds()
    t, rmin, rmax = T.frozen_adam_run(1e-3, 0.1, 1e-4, 0.9, 0.999, 0.0, 1.0)
    assert c <= rmin and rmax <= C and 0.99 <= rmin <= 1.01 and 1.5 < rmax < 3.0
    assert (0.1 - 1e-3) / (1e-4 * C) <= t <= (0.1 - 1e-3) / (1e-4 * c)               # Prop 2b(iii)
    # an eps_adam far above the gradients: momentum SGD with rate lr/eps_adam, ordered by the drive
    t1 = T.frozen_adam_run(0.05, 0.5, 0.5, 0.9, 0.999, 10.0, 1.0)[0]
    t2 = T.frozen_adam_run(0.05, 0.5, 0.5, 0.9, 0.999, 10.0, 0.1)[0]
    assert 8.0 < t2 / t1 < 12.0                                                        # time ~ 1/g


def test_the_pair_peak_and_the_modifications():
    kinds = dict(gd=dict(lr=0.05), sign=dict(lr=1e-3), adam=dict(lr=1e-3, eps_adam=T.EA_TINY), adam_g=dict(lr=1e-3), ngd=dict(lr=1e-3))
    out = {}
    for k, kw in kinds.items():
        n = 6000 if k == 'gd' else 3000
        out[k] = T.run_pair(k, kw, 0.5, 1.31e-2, n, 1e-3)
    assert abs(out['sign'][2] - 0.5 / 1.5) < 0.01 and out['sign'][0] == out['sign'][1]
    assert out['adam'][2] > 0.2 and out['adam'][0] == out['adam'][1]
    for k in ('gd', 'adam_g', 'ngd'):                                                  # the mechanism is kept: atom 2 stays dormant
        assert out[k][2] < 1e-3 and np.isnan(out[k][1]) and out[k][0] > 0
    assert abs(out['sign'][3][0] - 1.0) < 0.01 and out['sign'][3][1] < 0.01            # the end state is the sparse fit up to the chatter


def test_eps_adam_restores_the_selection():
    lo = T.run_pair('adam', dict(lr=5e-4, eps_adam=1e-8), 0.5, 1.31e-2, 4000, 1e-3)
    hi = T.run_pair('adam', dict(lr=0.05, eps_adam=1.0), 0.5, 1.31e-2, 27000, 1e-3)
    assert lo[2] > 0.2 and lo[0] == lo[1]
    assert hi[2] < 1e-3 and np.isnan(hi[1]) and hi[0] > 400                            # GD-like: late, ordered, atom 2 never forms


def test_duplicates_share_under_sign_descent_and_gd_picks_the_larger_start():
    gd = T.run_duplicates('gd', dict(lr=0.05), 1.31e-2, 0.7, 2700)
    sg = T.run_duplicates('sign', dict(lr=1e-3), 1.31e-2, 0.7, 4000)
    ad = T.run_duplicates('adam', dict(lr=1e-3, eps_adam=T.EA_TINY), 1.31e-2, 0.7, 4000)
    assert gd[1] / gd[0] < 1e-3 and gd[0] > 0.99
    assert abs(sg[0] + sg[1] - 1.0) < 0.01                                                # the fit is shared
    assert 0.98 < sg[1] / sg[0] < 1.0 and 0.98 < ad[1] / ad[0] < 1.02                  # w_2/w_1 = (3 s_2^2 + a^2)/(3 s_1^2 + a^2) < 1
    a = (sg[0] / 2.0) ** (1 / 3)
    assert abs(sg[1] / sg[0] - (3 * (0.7 * 1.31e-2) ** 2 + a * a) / (3 * 1.31e-2 ** 2 + a * a)) < 5e-3     # the closed form, with a from w_1 = 6 s^2 a + 2 a^3 ~ 2a^3


def test_ngd_has_the_orbits_of_gradient_descent():
    G, c = T.pair_gram(0.5)

    def orbit(opt, n):
        W, _ = T.run_units(opt, G, c, 0.3, n, 1e-3)
        return W[int(np.argmax(W[:, 0] >= 0.3))]
    a, b = orbit(T.Opt('gd', 0.002), 20000), orbit(T.Opt('ngd', 2e-4), 20000)
    assert abs(a[0] - 0.3) < 0.005 and abs(b[0] - 0.3) < 0.005 and abs(a[1] / b[1] - 1) < 0.02


# ---- F2 --------------------------------------------------------------------------------------------------------------------------
def test_chain_quadrature_matches_the_ode_and_the_leading_term():
    for D, s in ((2, 0.2), (3, 0.2), (4, 0.1), (5, 0.1)):
        z_t = 0.5 ** (1.0 / D)
        assert abs(T.chain_time_quad(D, 1.0, s, z_t) / T.chain_time_ode(D, 1.0, s, z_t) - 1) < 1e-7
    assert abs(T.chain_time_quad(4, 1.0, 0.02, 0.5 ** 0.25) / T.chain_time_leading(4, 1.0, 0.02) - 1) < 0.01
    assert abs(T.chain_time_quad(2, 1.0, 1e-6, 0.5 ** 0.5) / T.chain_time_leading(2, 1.0, 1e-6) - 1) < 0.1
    e = [np.log(T.chain_time_quad(5, 1.0, s, 0.5 ** 0.2)) for s in (0.02, 0.01)]
    assert abs((e[1] - e[0]) / np.log(0.5) + 3.0) < 0.02                               # T ~ s^-(D-2)


def test_diagonal_two_layer_closed_form():
    t = diag_t = T.diag_two_layer_time(1.0, 1e-4)
    assert abs(T.diag_two_layer_w(t, 1.0, 1e-4) - 0.5) < 1e-12 and abs(diag_t - np.log(1e4 - 1) / 2) < 1e-12
    assert abs(T.diag_gd_check(0.5, 0.1, 1e-3, 4.0) - float(T.diag_two_layer_w(4.0, 0.5, 0.01))) < 1e-4
    g1, g2 = T.diag_two_layer_time(1.0, 1e-16), T.diag_two_layer_time(0.5, 1e-16)
    assert g2 / g1 < 2.0 and (g2 - g1) < 20                                            # the gap is of order ln(1/eps) transits, not 1/s


def test_the_paper_start_is_stationary_in_two_layers():
    assert np.all(np.array(T.pair_rhs(0.0, [0.0, 0.3, 0.0, 0.3], 1.0)) == 0.0)
    kt = T.KanToy('fact')
    P, bias = kt.init(1e-3)
    G, gb, rms, W1, W2 = kt.grads(P, bias)
    assert max(float(np.abs(g).max()) for g in G.values()) < 1e-15 and abs(gb) < 1e-15
    r = kt.train(1e-3, 200)
    assert r['first_grad'] < 1e-15 and r['amax'] < 1e-12 and r['n_units1'] == 0 and abs(r['rms'] - 1.0) < 1e-9
    r = kt.train(1e-3, 200, kind='adam', lr=0.01)
    assert r['amax'] < 1e-8 and r['n_units1'] == 0
    h = T.KanToy('free').train(1e-3, 3000)
    assert h['n_units1'] > 0 and h['rms'] < 0.2                                          # the hybrid leaves the saddle


def test_pair_time_law_and_symmetric_quadrature():
    t1, t2 = T.pair_time(0.4, 1e-3), T.pair_time(0.4, 1e-6)
    assert abs((t2 - t1) * 0.4 ** 4 - np.log(1e3)) < 0.1                                # T s^4 = ln(s/a_0) + const
    assert abs(T.pair_time(0.4, 1e-6) / T.pair_time_sym(0.4, 0.5e-6) - 1) < 1e-3       # alpha = (0, a_0) behaves as the symmetric (a_0/2, a_0/2)
    assert T.pair_time_sym(0.4, 1e-6) < T.pair_time(0.4, 1e-6)


def test_hybrid_conservation_and_bounds():
    for s, b0 in ((0.1, 1.0), (0.1, 0.3), (0.1, np.sqrt(1 / 3) * 0.1)):
        Tm, cons, delta, b, g0, npos = T.hybrid_run(s, b0, 0.2)
        lo, up = T.hybrid_bounds(s, b0, 0.2, g0, delta, npos)
        assert cons < 1e-8 and lo <= Tm <= up and delta < 0.05
        assert np.all(b >= b0 - 1e-12)
    Tm1 = T.hybrid_run(0.1, 1.0, 0.2)[0]
    assert abs(Tm1 / ((1 / 0.1 - 1 / 0.2) / (3 * T.hybrid_run(0.1, 1.0, 0.2)[4] * 1.0)) - 1) < 0.05      # the one-layer law with the drive b_0 g
    pure = T.hybrid_run(0.05, np.sqrt(1 / 3) * 0.05, 0.2)[0]
    assert pure > 10 * T.hybrid_run(0.05, 1.0, 0.2)[0]                                   # the pure small start is much slower


# ---- the additive KAN (numpy) ------------------------------------------------------------------------------------------------------
def _loss(m, X, y):
    return 0.5 * float(np.mean((y - m.fwd(X)[0]) ** 2))


def test_additive_kan_gradient_is_the_negative_gradient_of_the_loss():
    rng = np.random.default_rng(5)
    X, y = rng.uniform(-2.5, 2.5, (40, 3)), rng.standard_normal(40)
    m = T.AddKan(3, 4, 1e-3)
    for k in ('a', 'c', 'ls'):
        m.p[k] = m.p[k] + 0.05 * rng.standard_normal(m.p[k].shape)
    m.p['n'] = m.p['n'] + 0.3 * rng.random(m.p['n'].shape)
    m.p['b'] = np.array([0.1])
    G = m.grad(X, y)
    for k in m.keys:
        idx = np.ndindex(*m.p[k].shape)
        for ij in list(idx)[:5]:
            h = 1e-6
            m.p[k][ij] += h
            lp = _loss(m, X, y)
            m.p[k][ij] -= 2 * h
            lm = _loss(m, X, y)
            m.p[k][ij] += h
            assert abs(G[k][ij] - (-(lp - lm) / (2 * h))) < 1e-6 * max(1.0, abs(G[k][ij])), (k, ij)


def test_additive_kan_vec_roundtrip_counts_and_fit_keeps_the_best_iterate():
    m = T.AddKan(2, 5, 1e-3)
    z = m.vec(m.p)
    m.setvec(z + 1.0)
    assert np.allclose(m.vec(m.p), z + 1.0) and m.p['a'].shape == (2, 5) and m.p['b'].shape == (1,)
    rng = np.random.default_rng(2)
    X = rng.uniform(-2.5, 2.5, (200, 2))
    y = 1.5 * np.exp(-0.5 * ((X[:, 0] - 0.5) / 0.4) ** 2)
    data = dict(Xtr=X[:120], ytr=y[:120], Xva=X[120:], yva=y[120:])
    m = T.AddKan(2, 6, 1e-2)
    rates = m.vec({k: np.full_like(m.p[k], 1.0 if k not in ('n', 'b') else 0.05) for k in m.keys})
    v0 = float(np.mean((m.fwd(data['Xva'])[0] - data['yva']) ** 2))
    t = T.fit_add(m, data, T.Opt('gd', rates), 1500)
    v1 = float(np.mean((m.fwd(data['Xva'])[0] - data['yva']) ** 2))
    assert t > 0 and v1 < 0.5 * v0
    cnt = m.counts(data['Xtr'])
    assert cnt.shape == (2,) and cnt[0] >= 1 and cnt[1] == 0                         # the relevant input has active units, the irrelevant none


# ---- regression values and the stubs that pin the run-time rules ------------------------------------------------------------------
def test_default_eps_adam_and_pre_rules():
    assert T.Opt('adam', 1.0).ea == 1e-8
    assert abs(float(T.Opt('adam', 1.0).step(np.zeros(1), np.array([1e-8]))[0]) - 0.5) < 1e-12        # g/(|g| + 1e-8)
    g, pre = np.array([3.0, 4.0]), np.array([2.0, 1.0])
    assert np.allclose(T.Opt('ngd', 1.0, pre=pre).step(np.zeros(2), g), np.array([6.0, 4.0]) / np.hypot(6.0, 4.0))
    assert np.allclose(T.Opt('adam', 1.0, eps_adam=0.5, pre=pre).step(np.zeros(2), g), T.Opt('adam', 1.0, eps_adam=0.5).step(np.zeros(2), g))   # pre is not used by adam
    a = T.Opt('adam_g', 1.0, eps_adam=0.0, pre=pre).step(np.zeros(2), g)
    assert np.allclose(a, (pre * g) / np.sqrt(np.mean((pre * g) ** 2)))


def test_d3_gradient_step_and_regression_values():
    P = np.array([[0.3, 0.2]])
    assert np.allclose(T.d3_dir(P, np.array([2.0])), [[3 * 0.09 * 2.0, -3 * 0.04 * 2.0]])
    assert np.allclose(T.Opt('gd', 0.1).step(P.ravel(), T.d3_dir(P, np.array([2.0])).ravel()), [0.3 + 0.1 * 0.54, 0.2 - 0.1 * 0.24])
    assert T.frozen_adam_run(1e-3, 0.1, 1e-4, 0.9, 0.999, 0.0, 1.0)[0] == 511
    assert abs(T.pair_time(0.4, 1e-6) / 497.90210320570174 - 1) < 1e-6
    assert abs(T.pair_time_sym(0.4, 0.5e-6) / 497.9021032289769 - 1) < 1e-6 and abs(T.pair_time_sym(0.4, 1e-6) / 470.8260414887198 - 1) < 1e-6
    Tm, cons, delta, b, g0, npos = T.hybrid_run(0.1, 1.0, 0.2)
    assert abs(Tm / 17.90099279008669 - 1) < 1e-6 and abs(g0 / 0.09453087204829419 - 1) < 1e-9 and abs(delta - 0.035207654255412835) < 1e-6
    r = T.KanToy('free').train(1e-3, 1500)
    assert r['n_units1'] == 6 and abs(r['rms'] - 0.05523924073076357) < 1e-9 and abs(r['w2max'] - 0.9103010468721006) < 1e-9


def test_pair_equals_the_symmetric_start_up_to_a0_squared():
    # alpha_2^2 - alpha_1^2 = a_0^2 is conserved by the two-layer pair, so (0, a_0) is the symmetric start a_0/2 up to O(a_0^2)
    for a0 in (1e-3, 1e-5):
        assert abs(T.pair_time(0.4, a0) / T.pair_time_sym(0.4, a0 / 2) - 1) < 50 * a0 ** 2 + 1e-9
    from scipy.integrate import solve_ivp
    sol = solve_ivp(T.pair_rhs, [0, 300], [0.0, 0.4, 1e-3, 0.4], args=(1.0,), rtol=1e-11, atol=1e-22, dense_output=True)
    z = sol.sol(np.linspace(0, 300, 7))
    assert np.allclose(z[2] ** 2 - z[0] ** 2, 1e-6, rtol=1e-6)


def test_kan_toy_gradient_matches_finite_differences_and_the_unit_forms():
    rng = np.random.default_rng(8)
    for second in ('fact', 'free'):
        kt = T.KanToy(second, K=2, n=60)
        P, bias = kt.init(1e-2)
        for k in P:
            P[k] = P[k] + 0.3 * rng.standard_normal(P[k].shape) * (1.0 if k in ('a1', 'a2', 'w2') else 0.2)
        st = kt.stats(P)                                                                  # the gradient treats the standardisation as constants (as the spike's)
        G = kt.grads(P, 0.1, stats=st)[0]

        def loss(PP):
            return 0.5 * kt.grads(PP, 0.1, stats=st)[2] ** 2
        for k in P:
            for ij in list(np.ndindex(*P[k].shape))[:4]:
                h = 1e-6
                Pp = {kk: v.copy() for kk, v in P.items()}
                Pm = {kk: v.copy() for kk, v in P.items()}
                Pp[k][ij] += h
                Pm[k][ij] -= h
                fd = -(loss(Pp) - loss(Pm)) / (2 * h)
                assert abs(G[k][ij] - fd) < 2e-5 * max(1.0, abs(fd)), (second, k, ij)


def test_additive_kan_clamp_amplitude_and_masks():
    m = T.AddKan(1, 3, 1e-2)
    m.p['n'][:] = 1.0
    m.p['a'][:] = np.array([[3.0, 0.0, 0.0]])
    X = np.array([[m.p['c'][0, 0]]])
    f = m.fwd(X)[0][0]
    sp = m.p['c'][0, 1] - m.p['c'][0, 0]
    assert abs(f - 3.0 * np.tanh(1.0) * 1.0) < 1e-12 and abs(sp - 2.5) < 1e-12           # w = nu^2 A tanh(alpha/A), A = 3
    m.p['ls'][:] = np.log(2.5)                                                           # between exp(.) = 2 and 3: not clamped
    assert abs(m.fwd(np.array([[m.p['c'][0, 0] + 2.5]]))[1][0, 0, 0] - np.exp(-0.5)) < 1e-12
    m.p['ls'][:] = np.log(10.0)
    assert abs(m.fwd(np.array([[m.p['c'][0, 0] + 3.0]]))[1][0, 0, 0] - np.exp(-0.5)) < 1e-12      # clamped at 3
    g = m.grad(np.array([[0.0], [1.0]]), np.array([1.0, 2.0]))
    assert np.all(g['ls'] == 0.0)                                                        # outside the clamp: no gradient to the width
    m.p['ls'][:] = np.log(0.01)
    assert np.all(m.grad(np.array([[0.0], [1.0]]), np.array([1.0, 2.0]))['ls'] == 0.0)


class _Stub:
    """A model whose validation MSE after k steps is vs[k // 25] (the parameter counts the steps)."""
    keys = ('x',)

    def __init__(self, vs):
        self.vs, self.p = vs, dict(x=np.array([0.0]))

    def fwd(self, X):
        return (np.full(3, np.sqrt(self.vs[min(int(self.p['x'][0]) // 25, len(self.vs) - 1)])),)

    def grad(self, X, y):
        return dict(x=np.array([1.0]))

    def vec(self, d):
        return d['x'].copy()

    def setvec(self, z):
        self.p['x'] = z.copy()


def test_fit_add_stopping_rule_and_best_iterate():
    data = dict(Xtr=None, ytr=None, Xva=None, yva=np.zeros(3))
    # best at step 50, then worse: armed (0.85 < 0.9), stops 2000 steps after the best, returns the iterate of step 50
    m = _Stub([1.0, 0.95, 0.85] + [0.9] * 200)
    t = T.fit_add(m, data, T.Opt('gd', 1.0), 10000)
    assert m.p['x'][0] == 50.0 and 2050 <= t <= 2075
    # never below 0.9 of the start: never armed, runs to the end; best (step 25) kept
    m = _Stub([1.0, 0.95] + [0.96] * 500)
    t = T.fit_add(m, data, T.Opt('gd', 1.0), 3000)
    assert t == 3000 and m.p['x'][0] == 25.0
    # validation is looked at every 25 steps: a minimum at step 25 is seen
    m = _Stub([1.0, 0.5] + [0.9] * 500)
    T.fit_add(m, data, T.Opt('gd', 1.0), 3000)
    assert m.p['x'][0] == 25.0


def test_train_opt_matches_the_repo_for_the_paper_unit_and_the_clock():
    P, FR, L, K, c0, cs, As, f = T.design()
    ro = P.make_readout('P0')
    ref = FR.train_adam(ro, K, f, c0, 1e-6, 40, 1e-3, n_t=5)
    mine = T.train_opt(FR, ro, K, f, c0, 1e-6, 40, T.Opt('adam', 1e-3, eps_adam=1e-8), n_t=5)
    assert np.allclose(ref['W'], mine['W'], rtol=1e-10, atol=1e-18) and np.allclose(ref['taus'], mine['taus'])
    assert np.allclose(ref['C'], mine['C'], rtol=1e-10, atol=1e-18) and mine['status'] == 'ok'


def test_registries_agree_with_the_stored_raw():
    if not os.path.exists(T.RAW):
        pytest.skip('no stored raw')
    raw = T.load_raw()
    assert np.allclose(raw['A_g'], T.G_A) and np.allclose(raw['A_s'], T.S_A) and np.allclose(raw['B_scales'], T.SCALES_B) and np.allclose(raw['B_ea'], T.EA_B)
    assert np.allclose(raw['D2_grid'], T.GRID_D2) and np.allclose(raw['H_ea'], T.EA_H) and raw['F_labels'].tolist().count('sign') == 2
    assert raw['K_counts'].shape == (len(T.TARGETS_K) * len(T.NOISES_K), len(T.OPTS_K), 5) and np.allclose(raw['K_cells'][::2, 0], T.TARGETS_K)
    assert [round(float(x), 6) for x in raw['E_names'][:2, 1]] == [T.LRS_E[0]] * 2


def test_parts_regenerate_the_stored_arrays():
    if not os.path.exists(T.RAW):
        pytest.skip('no stored raw')
    raw = T.load_raw()
    for part in (T.part_A, T.part_G):
        out = part()
        for k, v in out.items():
            assert np.allclose(np.asarray(v, float), np.asarray(raw[k], float), rtol=1e-9, atol=1e-300, equal_nan=True), k


# ---- the design ------------------------------------------------------------------------------------------------------------------
def test_design_sign_descent_forms_all_units_together():
    P, FR, L, K, c0, cs, As, f = T.design()
    thr = FR.formed_threshold(As)
    res = FR.train_adam(P.make_readout('D3'), K, f, c0, 1e-6, 400, 1e-3, n_t=401, betas=(0.0, 0.0), eps_adam=1e-40)
    ft = T.formation_steps(FR, res, thr, 400)
    assert np.sum(~np.isnan(ft)) == 21 and np.nanmax(ft) - np.nanmin(ft) < 1.5
    g0 = K.corr(FR.pts(c0), f)[0]
    assert np.abs(g0).max() / np.abs(g0).min() > 30                                    # drives differ by > 30, the formation steps by < 1.5


def test_formation_steps_scales_the_clock():
    P, FR, L, K, c0, cs, As, f = T.design()
    res = dict(taus=np.linspace(0, 1, 5), W=np.array([[0.0], [0.001], [0.01], [0.1], [1.0]]))
    assert abs(T.formation_steps(FR, res, 0.01, 40)[0] - 20.0) < 1e-9


# ---- the report ------------------------------------------------------------------------------------------------------------------
def test_report_regenerates_the_stored_md():
    if not os.path.exists(T.RAW):
        pytest.skip('no stored raw')
    raw = T.load_raw()
    md = T.report(raw)
    with open(T.MD) as f:
        assert f.read() == md


def test_report_content_follows_the_raw():
    if not os.path.exists(T.RAW):
        pytest.skip('no stored raw')
    raw = T.load_raw()
    md = T.report(raw)
    A = raw['A_T_sign']
    assert np.all(A == A[:, :1]) and np.all(raw['A_T_gd'][:, 0] * 1e5 / raw['A_T_gd'][:, -1] > 0)
    assert 'sign descent 1, 1, 1' in md and 'gradient descent 1e+05' in md
    assert np.all(raw['B_steps'][0] == raw['B_steps'][0][0])
    assert np.all((raw['G_hybrid'][:, 3] <= raw['G_hybrid'][:, 2]) & (raw['G_hybrid'][:, 2] <= raw['G_hybrid'][:, 4]))
    F = raw['F_rows']
    labs = [str(x) for x in raw['F_labels']]
    assert F.shape == (len(labs), 5) and len(labs) == 24
    n = len(labs) // 2
    for off in (0, n):                                                                  # sign descent and default-eps Adam form all 21 units, with either readout
        assert F[labs.index('sign', off)][0] == 21 and F[labs.index('adam 1e-08', off)][0] == 21 and F[labs.index('adam 1e-08', off)][2] == 0.0
        assert F[labs.index('adam shared v', off)][0] == 3 and F[labs.index('adam shared v', off)][2] == 1.0
    d3 = [int(F[labs.index('adam %g' % e)][0]) for e in T.EA_F]
    assert d3 == [21, 10, 8, 5, 3, 3, 3, 3, 3] and F[labs.index('adam 1')][2] == 1.0    # D3: graded restoration of the sliding path's set
    Kc, truth = raw['K_counts'], raw['K_truth']
    over = [np.maximum(Kc[:, i, :] - truth, 0).sum(1).mean() for i in range(Kc.shape[1])]
    assert over[1] > 3 * over[0] and over[4] < over[1] and over[2] < over[1] and over[3] < over[1]   # GD against Adam on the additive KAN
    ex = [(Kc[:, i, :] == truth).mean() for i in range(Kc.shape[1])]
    assert ex[0] > 0.6 and ex[1] < 0.3 and ex[4] > 0.6
    assert np.all(raw['H_rows'][:, 0] == 3) and np.all(raw['H_rows'][:, 1] == 1) and np.all(raw['H_rows'][:, 2] == 1)     # P0, long horizon
    E, nm = raw['E_T'], raw['E_names']
    for T_, n_ in zip(E, nm):
        assert np.sum(~np.isnan(T_)) == 21
        if n_[2] == 1:
            assert np.nanmax(T_) / np.nanmin(T_) < 1.01                            # sign descent: one step for all 21 units (times interpolated between records)
        else:
            assert np.nanmax(T_) / np.nanmin(T_) < 1.3
    assert float(np.abs(raw['E_g0']).max() / np.abs(raw['E_g0']).min()) > 50
    peaks = raw['D2_peaks']
    cross = []
    for pk in peaks:
        bad = np.flatnonzero(pk >= T.THW_C)
        cross.append(float(raw['D2_grid'][0]) if bad.size == 0 else (float(raw['D2_grid'][bad.max() + 1]) if bad.max() + 1 < len(pk) else np.inf))
    assert cross[0] < cross[2] < cross[4] and cross[1] < cross[3] < cross[5]            # the crossover eps*/s grows with rho (0.2, 0.5, 0.8)
    Cr = raw['C_res']
    sg = Cr[(Cr[:, 2] == 1)]
    assert np.all(np.abs(sg[:, 5] - sg[:, 0] / (1 + sg[:, 0])) < 0.02 * (1 + sg[:, 0])) and np.all(sg[:, 3] == sg[:, 4])
    gd = Cr[(Cr[:, 2] == 0)]
    assert np.all(gd[:, 5] < T.THW_C) and np.all(np.isnan(gd[:, 4]))
