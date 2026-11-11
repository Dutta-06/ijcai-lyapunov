"""
Proposition D‴.40 (PARADIGM.md §5d): the gauge under general noise. Checked here: the noise mask of the dithered integrator dithers only
the chosen coordinates (the amplitudes then move by their drift alone); the averaged fitness of the trace formula is, for isotropic noise
2τ diag(G), minus Proposition D‴.38's correction at κ_F = 1 (−(τ/2)∂_βℓ along M), the same for noise on the amplitudes only, and zero to
leading order for noise on the centres and widths only.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
from src.core import gauge_noise as GN                                   # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def test_the_noise_mask_dithers_only_the_chosen_coordinates():
    # Proposition D‴.40's mask: with (0, 1, 1) the amplitudes move by their drift alone, the centres and widths as without the mask
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    T0 = E.exact_rep(net, np.ones(3) / 3) + 0.01
    out = {}
    for mask in (None, (0.0, 1.0, 1.0)):
        rec = []
        loop.run_dithered(T0, np.ones(3), 0.25, 0.25, 1e-3, np.random.default_rng(3), rec_every=1, callback=lambda s_, T, b: rec.append(T.copy()),
                          mask=mask)
        out[mask] = rec[-1]
    drift = loop.drift(T0, np.ones(3))[0]
    assert np.allclose(out[(0.0, 1.0, 1.0)][:, 0], T0[:, 0] + 0.25 * drift[:, 0], rtol=0, atol=1e-15)
    assert np.array_equal(out[(0.0, 1.0, 1.0)][:, 1:], out[None][:, 1:]) and not np.allclose(out[None][:, 0], out[(0.0, 1.0, 1.0)][:, 0])


def _instance(W):
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    beta = 3 * np.asarray(W, float)
    G9 = np.tile(E.G, 3)
    gradJ = lambda t: -loop.drift(t.reshape(3, 3), beta)[0].ravel() / G9
    u = lambda t: loop.fitness(loop.residual(t.reshape(3, 3), beta), t.reshape(3, 3))
    return net, base, loop, beta, G9, gradJ, u


def test_the_trace_formula_reduces_to_the_log_determinant_for_isotropic_noise_and_vanishes_without_amplitude_noise():
    tau = 1e-4
    W = np.array([0.30, 0.35, 0.35])
    net, base, loop, beta, G9, gradJ, u = _instance(W)
    t = E.exact_rep(net, W).ravel()
    iso = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(G9), G9)[0]
    fix = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=1.0, fixman_tau=tau).fixman_correction(t.reshape(3, 3), beta)
    flow = lambda ub: GN.replicator_flow(W, ub)
    assert np.abs(flow(iso) + flow(fix)).max() < 1e-5 * np.abs(flow(fix)).max()            # ū = −(τ/2)∂_βℓ up to a constant (2·10⁻⁷)
    amp = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(np.tile([1.0, 0.0, 0.0], 3)), G9)[0]
    assert np.abs(flow(amp) - flow(iso)).max() < 1e-5 * np.abs(flow(iso)).max()
    cs = GN.averaged_fitness(gradJ, u, t, 2 * tau * np.diag(np.tile([0.0, 0.25, 0.25], 3)), G9)[0]
    assert np.abs(flow(cs)).max() < 1e-5 * np.abs(flow(iso)).max()


def test_only_each_amplitudes_own_noise_variance_drives_its_class():
    # the generator identity (Proposition D‴.40(i)) at leading order: ū_k = −τ_k g′_k at the fit, τ_k = Σ_{α_kα_k}/(2G_α), g_k = a_k/(β_k a′_k),
    # i.e. ū_k = −(τ_k/m)(1 + (a_k/a_max)²)/(W_k a′_k) — for a correlated noise on all nine coordinates, the correlations and the centres'
    # and widths' noise dropping out at this order (beyond it E[g′_k] depends on all of Σ)
    W = np.array([0.30, 0.35, 0.35])
    net, base, loop, beta, G9, gradJ, u = _instance(W)
    t = E.exact_rep(net, W).ravel()
    a = base.features(E.exact_rep(net, W))[0]
    q = (a / net.a_max) ** 2
    rng = np.random.default_rng(11)
    for _ in range(2):
        L = rng.standard_normal((9, 9))
        S = 2e-4 * 0.3 * (L @ L.T)
        tau_k = np.diag(S)[0::3] / (2 * G9[0])
        ub = GN.averaged_fitness(gradJ, u, t, S, G9)[0]
        assert np.abs(ub / (-(tau_k / 3) * (1 + q) / (W * (1 - q))) - 1).max() < 1e-5


def test_the_weighted_analytic_centre():
    import exp_gauge_fixman as Q
    import exp_gauge_noise as N
    assert np.abs(N.w_sharp_weighted([1.0, 1.0, 1.0]) - Q.w_sharp()).max() < 1e-12
    tk = N.levels(N.noises()['correlated 1'])
    w = N.w_sharp_weighted(tk)
    Wm = np.abs(E.A_BUMP) / 4.0
    grad = tk * (1 + Wm ** 2 / w ** 2) / (w - Wm ** 2 / w)                   # Lagrange: the weighted partials equal
    assert abs(w.sum() - 1) < 1e-12 and np.ptp(grad) < 1e-9 * grad.mean()


def test_the_correlated_noise_path_of_the_dithered_integrator():
    """`run_dithered` with `noise_L`: one step is Θ + Δt·drift + √Δt·L ξ with ξ the generator's next standard normal draws (referee round 31:
    the path was untested)."""
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    T0 = E.exact_rep(net, np.ones(3) / 3) + 0.01
    L = np.linalg.cholesky(np.eye(9) + 0.3 * np.ones((9, 9)))
    rec = []
    loop.run_dithered(T0, np.ones(3), 0.25, 0.25, 1e-3, np.random.default_rng(5), rec_every=1, callback=lambda s_, T, b: rec.append(T.copy()),
                      noise_L=L)
    xi = np.random.default_rng(5).standard_normal(9)
    expect = T0 + 0.25 * loop.drift(T0, np.ones(3))[0] + (np.sqrt(0.25) * L @ xi).reshape(3, 3)
    assert np.allclose(rec[-1], expect, rtol=0, atol=1e-14)


def test_the_levels_and_the_per_class_correction():
    import exp_gauge_noise as N
    tau = 1e-4
    assert np.allclose(N.levels(2 * tau * np.diag(np.tile(E.G, 3))), tau)                     # τ_k = Σ_{α_kα_k}/(2G_α)
    assert np.allclose(N.levels(N.noises(tau)['correlated E']), np.array(N.LEVELS_E) * tau)
    cE = N.corr_E()
    assert np.allclose(np.diag(cE)[1::3], E.G[1]) and np.allclose(np.diag(cE)[2::3], E.G[2]) and np.linalg.eigvalsh(cE).min() > 0
    net, base = E.setup()
    W = np.array([0.30, 0.35, 0.35]); Th = E.exact_rep(net, W); beta = 3 * W
    A0 = base.features(Th)[0] * W
    one = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=tau, fixman_A=A0).fixman_correction(Th, beta)
    per = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, fixman=2.0, fixman_tau=np.array(N.LEVELS_E) * tau,
                                    fixman_A=A0).fixman_correction(Th, beta)
    assert np.allclose(per, np.array(N.LEVELS_E) * one, rtol=1e-14)
    a = base.features(Th)[0]; q = (a / net.a_max) ** 2
    assert np.allclose(N.closed(W, [1.0, 2.0, 3.0]), -(np.array([1.0, 2.0, 3.0]) / 3) * (1 + q) / (W * (1 - q)))
    ub = np.array([0.3, -0.1, 0.2])
    assert np.allclose(GN.replicator_flow(W, ub), W * (ub - W @ ub)) and abs(GN.replicator_flow(W, ub).sum()) < 1e-15


def test_the_gauge_identity_pointwise_and_the_generator_applied_to_psi():
    """u_k = −g_k∂_{α_k}J_fit at every state (g_k = a_k/(β_ka′_k)); Part D's Lψ_k = u_k − g_k∂_{α_k}(J_coh + V) + τ_kg′_k equals
    ψ′b + (Σ_αα/2)ψ″ for ψ_k = (a_max²/4β_kG_α)cosh(2α_k/a_max); the chain's conditional expectation of ψ(α + bΔt + ν) in closed form."""
    import exp_gauge_plateau as P
    net, base = E.setup()
    beta = np.array([0.9, 1.05, 1.05])
    fit = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    rng = np.random.default_rng(2)
    Th = E.exact_rep(net, beta / 3) + 0.05 * rng.standard_normal((3, 3)); Th[0, 0] = 5.6        # α_1 outside the box
    a = base.features(Th)[0]; am = net.a_max; ap = 1 - (a / am) ** 2
    g = a / (beta * ap); gp = (1 + (a / am) ** 2) / (beta * ap)
    r = fit.residual(Th, beta); u = fit.fitness(r, Th)
    Jfit = lambda T: 0.5 * 3 * np.mean(fit.residual(T, beta) ** 2)                                # J = m·½⟨r²⟩, m = 3
    h = 1e-6
    dJ = np.array([(Jfit(Th + h * np.eye(3)[k][:, None] * np.eye(3)[0][None, :]) - Jfit(Th - h * np.eye(3)[k][:, None] * np.eye(3)[0][None, :]))
                   / (2 * h) for k in range(3)])
    assert np.allclose(u, -g * dJ, rtol=1e-6, atol=1e-12)
    Sig = 2e-3 * np.diag(np.tile(E.G, 3)); tk = np.diag(Sig)[0::3] / (2 * E.G[0])
    dT = loop.drift(Th, beta)[0]
    nonfit = -(dT[:, 0] - beta * E.G[0] * base.grad_weighted(r, Th)[:, 0]) / E.G[0]
    lpsi = u - g * nonfit + tk * gp
    c = am ** 2 / (4 * beta * E.G[0]); kk = 2 / am
    psi1 = c * kk * np.sinh(kk * Th[:, 0]); psi2 = c * kk * kk * np.cosh(kk * Th[:, 0])
    assert np.allclose(lpsi, psi1 * dT[:, 0] + 0.5 * np.diag(Sig)[0::3] * psi2, rtol=1e-10)
    assert np.allclose(psi1, g / E.G[0], rtol=1e-12) and np.allclose(psi2, gp / E.G[0], rtol=1e-12)
    x, wq = np.polynomial.hermite_e.hermegauss(60)                                                # E over ν ~ N(0, s²)
    s2, dt = np.diag(Sig)[0::3] * 0.25, 0.25
    num = np.array([np.sum(wq * c[k] * np.cosh(kk * (Th[k, 0] + dt * dT[k, 0] + np.sqrt(s2[k]) * x))) / np.sqrt(2 * np.pi) for k in range(3)])
    assert np.allclose(num, c * np.cosh(kk * (Th[:, 0] + dt * dT[:, 0])) * np.exp(0.5 * kk * kk * s2), rtol=1e-12)


# ---- referee round 147 (D3ep): the bracket of w_sharp_weighted, Part E's predictions (F), the early phase (G), coloured noise (H) ----

def _N():
    import exp_gauge_noise as N
    return N


def _raw():
    N = _N()
    return N.E.load(os.path.join(N.OUT, 'gauge_noise_raw.npz'))


def _direct_centre(tk):
    from scipy.optimize import minimize
    Wm = np.abs(E.A_BUMP) / 4.0; tk = np.asarray(tk, float)

    def F(x):
        W = np.r_[x, 1 - x.sum()]
        return np.inf if np.any(W <= Wm) else -np.sum(tk * np.log(W - Wm ** 2 / W))
    r = minimize(F, [0.4, 0.3], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-15, maxiter=20000))
    return np.r_[r.x, 1 - r.x.sum()]


def test_the_weighted_centre_at_small_levels_and_outside_the_box():
    # D3ep-5: the levels (2, 1, 0.2) raised a sign error in brentq; their centre lies below the box's wall W_B,3 (outside P_B)
    import exp_gauge_plateau as P
    N = _N()
    for lv in ((2.0, 1.0, 0.2), (1.0, 1.0, 0.1), (2.0, 1.0, 0.5), (1.0, 3.0, 0.05)):
        w = N.w_sharp_weighted(lv)
        assert np.abs(w - _direct_centre(lv)).max() < 1e-7 and abs(w.sum() - 1) < 1e-12
    assert N.w_sharp_weighted((2.0, 1.0, 0.2))[2] < P.W_BOX[2] and np.all(N.w_sharp_weighted((2.0, 1.0, 0.5)) > P.W_BOX)
    for bad in ((1.0, 1.0, 0.0), (1.0, -1.0, 1.0)):
        try:
            N.w_sharp_weighted(bad)
        except ValueError as e:
            assert 'positive' in str(e)
        else:
            raise AssertionError('a zero or negative level must be refused')


def test_the_walls_of_P_B():
    import exp_gauge_plateau as P
    N = _N()
    WB = P.W_BOX
    W = np.array([WB[0], 0.40, 1 - WB[0] - 0.40])
    b = np.array([-2.0, 0.5, 1.5])
    out = N.wall_project(W, b)                                     # the wall's component cancelled along e_1 − W, λ ≥ 0
    lam = (out - b)[1] / -W[1]
    assert abs(out[0]) < 1e-15 and abs(out.sum()) < 1e-15 and lam > 0 and np.allclose(out - b, lam * (np.eye(3)[0] - W), atol=1e-15)
    assert np.array_equal(N.wall_project(W, -b), -b) and np.array_equal(N.wall_project(np.ones(3) / 3, b), b)
    V = np.array([WB[0], WB[1], 1 - WB[0] - WB[1]])                  # a vertex: both walls hold
    out = N.wall_project(V, np.array([-1.0, -1.0, 2.0]))
    assert np.abs(out).max() < 1e-12
    X_ = np.array([[WB[0] - 0.01, 0.40, 0.0]]); X_[0, 2] = 1 - X_[0, :2].sum()
    c, r = N.wall_push(X_, 1.0)[0], N.wall_push(X_, 2.0)[0]
    assert abs(c[0] - WB[0]) < 1e-15 and abs(r[0] - (WB[0] + 0.01)) < 1e-15 and abs(c.sum() - 1) < 1e-15 and abs(r.sum() - 1) < 1e-15
    assert np.allclose((r - X_[0]) / 0.02, (np.eye(3)[0] - X_[0]) / (1 - X_[0, 0]))
    g = N.grid_PB(4)
    assert len(g) == 15 and np.allclose(g.sum(axis=1), 1) and np.all(g >= WB - 1e-15) and np.isclose(g[:, 0].max(), 1 - WB[1] - WB[2])


def test_the_slow_limits_leading_order_fitness():
    N = _N()
    W = np.array([0.31, 0.36, 0.33]); tk = N.levels(N.noises()['correlated E'])
    assert np.allclose(N.lead_fitness('correlated E', W), N.closed(W, tk))
    assert np.allclose(N.lead_fitness('correlated E + class levels', W), -N.closed(W, tk))           # κ_F = 2: the noise reversed
    assert np.allclose(N.lead_fitness('nonlinear only + barrier', W), -N.closed(W, np.full(3, N.TAU)))
    assert np.abs(N.lead_fitness('nonlinear only', W)).max() == 0
    assert np.allclose(N.e_target('correlated E + class levels'), N.w_sharp_weighted(tk)) and N.e_target('correlated E') is None


def test_part_F_the_slow_limit_recomputed():
    N = _N(); F = _raw()['F']
    for v in ('nonlinear only + barrier', 'correlated E', 'correlated E + class levels'):
        o, r = N.f_ode(v), F['var'][v]['ode']
        assert o['t_wall'] == r['t_wall'] and np.array_equal(o['half'], r['half'])
        assert all(np.array_equal(o['rec'][t], r['rec'][t]) for t in N.T_F)
    assert F['var']['correlated E']['ode']['t_wall'] == 22180.0 and F['var']['correlated E + class levels']['ode']['t_wall'] is None


def test_part_F_katzenbergers_field_and_paths_recomputed():
    N = _N(); F = _raw()['F']
    g = N.grid_PB()
    for v in ('correlated E', 'correlated E + class levels', 'nonlinear only + barrier'):
        fld = F['var'][v]['field']
        assert np.array_equal(fld[:, :2], g[:, :2])
        Kv = N.katz(v, N.GAMMA_E)
        for i in (0, 417, len(g) - 1):
            b, D = N.k_drift(Kv, g[i])
            assert np.allclose(fld[i, 2:], np.r_[b, D[0, 0], D[0, 1], D[1, 1]], rtol=1e-9, atol=1e-25)
        b13 = N.k_drift(Kv, np.ones(3) / 3)[0]
        assert np.allclose(F['var'][v]['b13'], np.r_[b13, -b13.sum()] / (N.GAMMA_E * N.TAU), rtol=1e-9, atol=1e-12)
        at = N.at_third(v)
        assert np.allclose(at['b13'], F['var'][v]['b13'], rtol=1e-12, atol=1e-15) and np.array_equal(at['lead13'], F['var'][v]['lead13'])
    lead = N.GN.replicator_flow(np.ones(3) / 3, N.closed(np.ones(3) / 3, N.levels(N.noises()['correlated E']))) / N.TAU
    assert np.allclose(F['var']['correlated E']['lead13'], lead)
    r = F['var']['correlated E']
    p = N.f_paths(r['field'], S_end=2e4)
    assert all(np.array_equal(p['rec'][t], r['sde']['rec'][t]) for t in (5e3, 1e4, 2e4))
    q = N.f_paths(r['field'], diffuse=False)                                                       # the drift's path, whole
    assert all(np.array_equal(q['rec'][t], r['kpath']['rec'][t]) for t in N.T_F) and np.array_equal(q['half'], r['kpath']['half'])


def test_part_F_the_corrected_drifts_zero_and_the_centres_in_P_B():
    N = _N(); F = _raw()['F']
    z = [r for r in F['zeros'] if r['variant'] == 'correlated E + class levels' and r['gamma'] == 1e-2][0]
    Kv = N.katz('correlated E + class levels', 1e-2)
    assert np.abs(N.k_drift(Kv, z['zero'])[0]).max() < 1e-9 * 1e-2 * N.TAU
    tg = N.e_target('correlated E + class levels')
    assert np.allclose(N.k_drift(Kv, tg)[0] / (1e-2 * N.TAU), z['at_target'], rtol=1e-9)
    assert np.abs(z['zero'] - tg).max() > 5e-3                                                       # the finite-γ shift
    z2 = N.f_zero.__wrapped__(('correlated E + class levels', 1e-2))
    assert np.abs(z2['zero'] - z['zero']).max() < 1e-9 and z2['resid'] < 1e-8
    for p in F['pb']:
        assert np.allclose(p['w'], N.w_sharp_weighted(p['levels']), rtol=0, atol=1e-14)
    assert [bool(np.all(p['margin'] > 0)) for p in F['pb']] == [True, False, False]


def test_part_G_the_early_phase_reproduced():
    N = _N(); raw = _raw(); G = raw['G']
    assert G['replay_diff'] == 0.0
    Er = np.asarray(raw['E_rec']['correlated E'])
    for r, s in zip(G['replay'], N.REPLAY_G):
        assert np.array_equal(r, Er[N.SEEDS_E.index(s), :len(r)]) and r[-1, 0] == N.S_G
    for v in N.VARIANTS_G:
        R = np.asarray(G['rec'][v])
        assert R.shape[0] == len(N.SEEDS_G) and R[0, -1, 0] == N.S_G
    g = N.g_series.__wrapped__((N.VARIANTS_G[1], N.SEEDS_G[0], 2e3))
    assert np.array_equal(g['rec'], np.asarray(G['rec'][N.VARIANTS_G[1]])[0, :len(g['rec'])])


def test_the_coloured_noise_cross_covariance():
    rng = np.random.default_rng(4)
    A = rng.standard_normal((5, 5)); GH = A @ A.T + 0.1 * np.eye(5)
    L = rng.standard_normal((5, 5)); S = L @ L.T
    for tc in (0.3, 10.0, 1e3):
        assert np.allclose(GN.coloured_cross(GH, S, tc), GN.coloured_cross_augmented(GH, S, tc), rtol=1e-9, atol=1e-12)
    assert np.allclose(GN.coloured_cross(GH, S, 0.0), S / 2)


def test_part_H_the_linear_response_recomputed():
    N = _N(); H = _raw()['H']
    lin = N.part_H_linear()
    assert np.array_equal(lin['rates'], H['linear']['rates']) and lin['check'] < 1e-10
    for a, b in zip(lin['rows'], H['linear']['rows']):
        assert a['name'] == b['name'] and a['tc'] == b['tc'] and np.array_equal(a['ratio'], b['ratio']) and np.array_equal(a['corr'], b['corr'])
    iso = [r for r in lin['rows'] if r['name'] == 'amplitudes only']
    assert all(np.all(np.diff([r['ratio'][k] for r in iso]) < 0) for k in range(3))                 # the level falls with τ_c
    assert np.abs(lin['rows'][0]['corr']).max() == 0                                                 # isotropic: no correlations
    W = np.ones(3) / 3; tk = N.levels(N.noises(N.TAU_H)['amplitudes only'])
    assert np.allclose(H['lead_col'], N.closed(W, H['ratio_col'] * tk)) and np.allclose(H['lead_white'], N.closed(W, tk))
    ph = N.part_H(H['mc'])
    assert np.array_equal(ph['ratio_col'], H['ratio_col']) and np.array_equal(ph['lead_col'], H['lead_col'])
    assert np.array_equal(H['ratio_col'], [r for r in lin['rows'] if r['name'] == 'amplitudes only' and r['tc'] == 100.0][0]['ratio'])


def test_part_H_the_coloured_simulation_step_by_step():
    """h_mc's coloured run against an independent integration over a few steps (one step per batch)."""
    import exp_gauge_plateau as P
    N = _N()
    tc, seed, n = 7.0, 3, 6
    h = N.h_mc.__wrapped__((tc, seed, n * N.DT, n))
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    beta = np.ones(3); sd = np.sqrt(np.diag(N.noises(N.TAU_H)['amplitudes only']))
    rng = np.random.default_rng(seed); Th = E.exact_rep(net, np.ones(3) / 3).copy()
    eta = sd / np.sqrt(2 * tc) * rng.standard_normal(9)
    for k in range(n):
        a = base.features(Th)[0]; g = a / (beta * (1 - (a / 4.0) ** 2))
        assert np.allclose(h['geta'][k], -g * eta.reshape(3, 3)[:, 0] / E.G[0], rtol=1e-12, atol=0)
        assert np.allclose(h['u'][k], loop.fitness(loop.residual(Th, beta), Th), rtol=1e-12)
        Th = Th + N.DT * (loop.drift(Th, beta)[0] + eta.reshape(3, 3))
        eta = np.exp(-N.DT / tc) * eta + np.sqrt(1 - np.exp(-2 * N.DT / tc)) * sd / np.sqrt(2 * tc) * rng.standard_normal(9)
    w = N.h_mc.__wrapped__((0.0, seed, n * N.DT, n))                                                # white: Part D's step
    assert np.array_equal(w['geta'], np.zeros((n, 3))) and np.allclose(w['gp'][0], N.closed(np.ones(3) / 3, N.levels(N.noises(N.TAU_H)['amplitudes only'])), rtol=1e-9)


def test_parts_F_G_H_of_the_md_regenerate_from_the_raw():
    N = _N(); raw = _raw()
    md = N.report(*[raw.get(k) for k in 'ABCDE'], **{k: raw.get(k) for k in ('E_rec', 'F', 'G', 'H', 'I', 'J', 'K', 'M', 'N')})
    with open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8') as f:
        assert f.read() == md
    for h in ('## F. Part E against its predictions', '## G. The early phase', '## H. Coloured noise'):
        assert h in md


def test_part_I_the_batch_is_run_dithered_and_the_decomposition_is_exact():
    """i_runs (one seed, 2·10³ slow units) against run_dithered's records (g_series) and against the raw's batched run; the gauge identity's
    decomposition sums to the log-masses' ratios; the no-box and equilibrated-start runs against the raw."""
    N = _N(); I = _raw()['I']
    sd = N.SEEDS_I[0]
    r = N.i_runs.__wrapped__(('box', (sd,), 2e3))
    g = N.g_series.__wrapped__(('correlated E', sd, 2e3))['rec']
    assert np.allclose(r['al'][0], g[:, 1:4], rtol=0, atol=1e-12) and np.allclose(r['W'][0], g[:, 4:7], rtol=0, atol=1e-13)
    n = r['W'].shape[1]
    for kind in N.KINDS_I:
        q = r if kind == 'box' else N.i_runs.__wrapped__((kind, (sd,), 2e3))
        R = I['runs'][kind]
        assert np.allclose(q['W'][0], R['W'][0, :n], rtol=0, atol=1e-13) and np.allclose(q['acc'][0], R['acc'][0, :n], rtol=1e-9, atol=1e-12)
        assert np.abs(N.i_identity(q)).max() < 1e-12
    assert not np.allclose(I['runs']['burn-in']['al'][:, 0], I['runs']['box']['al'][:, 0])                # the equilibrated start
    assert np.array_equal(I['runs']['burn-in']['W'][:, 0], np.full((len(N.SEEDS_I), 3), 1 / 3))
    for R in I['runs'].values():
        assert np.abs(N.i_identity(R)).max() < 1e-10 and R['W'].shape[1] == int(N.S_I / 1e3) + 1
    raw = _raw(); old = I['runs']['old']; nG = len(N.SEEDS_G)
    assert np.array_equal(old['W'][:nG], np.asarray(raw['G']['rec']['correlated E'])[:, :, 4:7]) or \
        np.abs(old['W'][:nG] - np.asarray(raw['G']['rec']['correlated E'])[:, :, 4:7]).max() < 1e-12


def test_part_I_the_numbers_recomputed_from_the_raw():
    N = _N(); raw = _raw(); I = raw['I']
    smp = N.i_samples(I)
    lab = f"the {len(N.SEEDS_I)} new seeds ({N.SEEDS_I[0]}–{N.SEEDS_I[-1]})"
    Rr = smp[lab][0]; i = int(np.argmin(np.abs(Rr[0, :, 0] - 2e4)))
    m = Rr[:, i, 4].mean(); e = Rr[:, i, 4].std(ddof=1) / np.sqrt(len(Rr))
    sde = I['sde']['reflected at W_B, step 50']['rec'][2e4][:, 0]
    assert np.array_equal(sde, raw['F']['var']['correlated E']['sde']['rec'][2e4][:, 0])                # Part F's ensemble
    z = (m - sde.mean()) / max(e, sde.std(ddof=1) / np.sqrt(len(Rr)))
    md = open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8').read()
    assert f"| {lab} | 20000 | {m:.4f} ± {e:.4f}," in md and f"| {z:+.1f} |" in md
    c = N.i_contrib(smp[lab][1][:, i])                                                                   # the terms add up to W_1 − 1/3
    assert abs(c.sum(axis=1).mean() - (Rr[:, i, 4] - 1 / 3).mean()) < 0.01
    assert np.allclose(N.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])
    assert abs(N._p2(1.959963985) - 0.05) < 1e-9
    fx = I['fixed']; cb = N.i_cbox(fx)
    mb = fx['box'][:, 0].reshape(len(fx['w1']), fx['ns']).mean(axis=1)
    ok = mb > 0
    assert np.allclose(cb(fx['w1'][ok]), mb[ok], rtol=1e-12) and np.all(np.diff(np.log(cb(np.linspace(0.29, 0.33, 50)))) < 0)
    assert np.allclose(I['transient'], N.i_transient(), rtol=1e-9, atol=1e-15)
    m, e, c, b, sm, ss = N.i_regressed(Rr[:, i, 4], N.i_contrib(smp[lab][1][:, i])[:, 1], I['sde']['reflected at W_B, step 2']['rec'][2e4][:, 0])
    w = Rr[:, i, 4]; mm = N.i_contrib(smp[lab][1][:, i])[:, 1]
    assert abs(b - np.polyfit(mm, w, 1)[0]) < 1e-10 and e < 0.75 * Rr[:, i, 4].std(ddof=1) / np.sqrt(len(Rr))   # the variance reduced
    assert f"| {lab} | 20000 | {c:+.3f} | {b:.3f} | {m:.4f} ± {e:.4f} | {sm:.4f} ± {ss:.4f} | {(m - sm) / np.hypot(e, ss):+.1f} |" in md


def test_part_I_the_sde_variants_and_the_field_recomputed():
    N = _N(); raw = _raw(); I = raw['I']
    fe = I['field']; g = N.grid_ext()
    assert np.array_equal(fe[:, :2], g[:, :2]) and np.all(g >= N.L_EXT_I - 1e-15)
    for i in (0, 431, len(g) - 1):
        assert np.allclose(N.i_field_row(N.katz('correlated E', N.GAMMA_E), g[i]), fe[i], rtol=1e-9, atol=1e-25)
    W = np.array([1 / 3, 1 / 3, 1 / 3]); Pi = N.i_field_row(N.katz('correlated E', 1e-6), W)[7:].reshape(2, 3)
    assert np.allclose(Pi, (1e-6 * W[:, None] * (np.eye(3) - W[None, :]))[:2], rtol=2e-3)                 # the slow limit's projection
    p = N.i_sde(fe, N.L_EXT_I, 10.0)
    assert all(np.array_equal(p['rec'][t], I['sde']['no walls, step 10']['rec'][t]) for t in N.T_I)
    p = N.i_sde(raw['F']['var']['correlated E']['field'], N.P.W_BOX, 10.0)
    assert all(np.array_equal(p['rec'][t], I['sde']['reflected at W_B, step 10']['rec'][t]) for t in N.T_I)
    p = N.i_sde(fe, N.L_EXT_I, 2.0, cbox=N.i_cbox(I['fixed']), n_paths=200)
    q = N.i_sde(fe, N.L_EXT_I, 2.0, n_paths=200)                                                          # the box raises W_1
    assert p['rec'][2e4][:, 0].mean() > q['rec'][2e4][:, 0].mean()


def test_part_I_the_fixed_masses_box_term_step_by_step():
    import exp_gauge_plateau as P
    N = _N()
    fx = N.i_fixed.__wrapped__(((0.2947,), 1, 2 * N.DT, 3 * N.DT))
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    W = np.array([0.2947, (1 - 0.2947) / 2, (1 - 0.2947) / 2]); beta = 3 * W; Th = E.exact_rep(net, W).copy()
    Ld = np.sqrt(N.DT) * np.linalg.cholesky(N.noises()['correlated E']); rng = np.random.default_rng(7000)
    Z = rng.standard_normal((N.REC, 9)); box = np.zeros(3)
    for k in range(5):
        if k >= 2:
            a = base.features(Th)[0]; box += a / (beta * (1 - (a / 4.0) ** 2)) * P.SIGMA_B * loop.box_excess(Th)[:, 0] / E.G[0]
        Th = Th + N.DT * loop.drift(Th, beta)[0] + (Ld @ Z[k]).reshape(3, 3)
    assert np.allclose(fx['box'][0], box / 3, rtol=1e-9, atol=1e-15) and box[0] > 0


def test_part_I_of_the_md_regenerates_from_the_raw():
    N = _N(); raw = _raw()
    md = N.report(*[raw.get(k) for k in 'ABCDE'], **{k: raw.get(k) for k in ('E_rec', 'F', 'G', 'H', 'I', 'J', 'K', 'M', 'N')})
    with open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8') as f:
        assert f.read() == md
    assert "## I. The early phase's gap" in md and md.index('## H. Coloured noise') < md.index("## I. The early phase's gap")


def test_part_I_the_loop_without_the_box_leaves_it():
    """Seed 507 without the box: its first amplitude passes the box's edge by 5·10³ slow units (the raw), and the recomputed run (no box term)
    reproduces it; with the box the same draws stay nearer the edge."""
    N = _N(); I = _raw()['I']
    j = N.SEEDS_I.index(507); R = I['runs']['no box']
    assert np.abs(R['al'][j, 5, 0]) > N.P.BOX[0] and np.all(R['acc'][:, :, 3] == 0)
    q = N.i_runs.__wrapped__(('no box', (507,), 5e3))
    assert np.allclose(q['al'][0], R['al'][j, :6], rtol=0, atol=1e-11) and np.allclose(q['W'][0], R['W'][j, :6], rtol=0, atol=1e-13)
    b = N.i_runs.__wrapped__(('box', (507,), 5e3))
    assert np.abs(b['al'][0, 5, 0]) < np.abs(q['al'][0, 5, 0]) and np.any(b['acc'][0, 5, 3] != 0)


# ---- referee round 153 (D3ev): Parts J, K and L ----
def _single_terms(N, loop, base, Th, beta, tk, var, h):
    """d_job's formulas at one state (the reference for j_runs): the identity's terms and the chain's bias at the step h."""
    import exp_gauge_plateau as P
    a = base.features(Th)[0]; ap = 1 - (a / N.A_MAX) ** 2
    g = a / (beta * ap); gp = (1 + (a / N.A_MAX) ** 2) / (beta * ap)
    r = loop.residual(Th, beta); u = loop.fitness(r, Th); dT = loop.drift(Th, beta)[0]
    nonfit = -(dT[:, 0] - beta * E.G[0] * base.grad_weighted(r, Th)[:, 0]) / E.G[0]
    boxt = P.SIGMA_B * loop.box_excess(Th)[:, 0] / E.G[0]
    lpsi = u - g * nonfit + tk * gp
    kk, cpsi = 2.0 / N.A_MAX, N.A_MAX ** 2 / (4 * beta * E.G[0])
    Ep = lambda hh: cpsi * np.cosh(kk * (Th[:, 0] + hh * dT[:, 0])) * np.exp(0.5 * kk * kk * var * hh)
    psi = cpsi * np.cosh(kk * Th[:, 0])
    return dT, dict(u=u, gp=-tk * gp, coh=g * (nonfit - boxt), box=g * boxt, lpsi=lpsi, b=(Ep(h) - psi) / h - lpsi), Ep


def test_part_J_the_coupled_chains_step_by_step():
    """j_runs over ten full steps (one per batch) against d_job's single-state formulas: both chains driven by one Brownian path (the
    half steps' increments √(Δt/2)Lz₁, √(Δt/2)Lz₂, the full step's √(Δt/2)L(z₁ + z₂)), the bias at the own and the other step."""
    import exp_gauge_plateau as P
    N = _N(); S = N.N_BATCH_J * N.DT
    J = N.j_runs.__wrapped__((S, 5))
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    B = len(N.RUNS_J); rng = np.random.default_rng(5); Z = rng.standard_normal((N.REC, 2, B, 9)); h = N.DT / 2; beta = np.ones(3)
    for j in (0, 3, 6):
        nm, tau = N.RUNS_J[j]; Sig = N.noises(tau)[nm]; tk = N.levels(Sig); var = np.diag(Sig)[0::3]
        lam, V = np.linalg.eigh(Sig); Lc = V * np.sqrt(np.clip(lam, 0, None))
        Tf = E.exact_rep(net, np.ones(3) / 3).copy(); Th = Tf.copy()
        for k in range(N.N_BATCH_J):
            z1, z2 = Z[k, 0, j], Z[k, 1, j]
            dT, t, Ep = _single_terms(N, loop, base, Tf, beta, tk, var, N.DT)
            for kn in ('u', 'gp', 'coh', 'box', 'lpsi', 'b'):
                assert np.allclose(J['full'][kn][k, j], t[kn], rtol=1e-9, atol=1e-16), (kn, k, j)
            Tf = Tf + N.DT * dT + (np.sqrt(h) * (Lc @ (z1 + z2))).reshape(3, 3)
            acc, alt = {kn: 0.0 for kn in t}, 0.0
            for z in (z1, z2):
                dT, t, Ep = _single_terms(N, loop, base, Th, beta, tk, var, h)
                for kn in t:
                    acc[kn] = acc[kn] + t[kn] / 2
                alt = alt + ((Ep(N.DT) - N.A_MAX ** 2 / (4 * E.G[0]) * np.cosh(2 * Th[:, 0] / N.A_MAX)) / N.DT - t['lpsi']) / 2
                Th = Th + h * dT + (np.sqrt(h) * (Lc @ z)).reshape(3, 3)
            for kn in acc:
                assert np.allclose(J['half'][kn][k, j], acc[kn], rtol=1e-9, atol=1e-16), (kn, k, j)
            assert np.allclose(J['half']['b_alt'][k, j], alt, rtol=1e-9, atol=1e-16)
    assert J['dt'] == (N.DT, h) and J['full']['out_box'].shape == (B, 3)


def test_part_J_the_bias_is_linear_in_the_step_and_the_md():
    """From the raw: the chains' bias halves with the step (the same states: 2 to 10⁻²; between the chains to 0.1), the right side without
    its Itô term leaves a remainder that does not (ratio 1 to 0.02), the Δt chain's bias agrees with Part D's in order of magnitude; the md
    lines recomputed."""
    N = _N(); raw = _raw(); J = raw['J']; D = raw['D']
    assert J['runs'] == N.RUNS_D and J['S'] == N.S_J and J['seed'] == N.SEED_J
    md = open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8').read()
    for j, (nm, tau) in enumerate(J['runs']):
        bf, bh, ba = J['full']['b'][:, j], J['half']['b'][:, j], J['half']['b_alt'][:, j]
        r = bf.mean(axis=0) / bh.mean(axis=0); rs = ba.mean(axis=0) / bh.mean(axis=0)
        e = (ba - rs * bh).std(axis=0, ddof=1) / np.sqrt(len(bh)) / np.abs(bh.mean(axis=0))
        assert np.all(np.abs(rs - 2) < 1e-2) and np.all(np.abs(r - 2) < 0.15), (nm, tau, r, rs)
        assert np.allclose(N.ratio_se(ba, bh)[1], e, rtol=1e-12)
        assert f"| {nm} | {tau:g} | {', '.join(f'{x:+.2e}' for x in bf.mean(axis=0))} | {', '.join(f'{x:+.2e}' for x in bh.mean(axis=0))} |" in md
        if np.all(J['tk'][j] > 0):
            wf, wh = bf - J['full']['gp'][:, j], bh - J['half']['gp'][:, j]
            assert np.all(np.abs(wf.mean(axis=0) / wh.mean(axis=0) - 1) < 0.02)
            d = [x for x in D if x['name'] == nm and x['tau'] == tau][0]
            assert np.all(np.abs(np.log(bf.mean(axis=0) / d['bias'].mean(axis=0))) < np.log(3))     # Part D's run, another seed
    assert "## J. The identity's own check" in md and md.index("## I. The early phase's gap") < md.index("## J. The identity's own check")


def test_part_K_the_replay_is_part_E_and_the_decomposition_is_exact():
    """k_runs (one seed each, 2·10³ slow units) reproduces Part E's records and the raw's replay; the decomposition sums to the
    log-masses' ratios; the chain's term at the first step is −(φ(α + Δt b_α) − φ(α) − Δtφ′(α)b_α)/β (the fit's fitness is Δtφ′b_α/β plus
    the box's and the cohesion's terms, by the gauge identity)."""
    N = _N(); raw = _raw(); K = raw['K']
    q = N.k_runs.__wrapped__((N.VARIANTS_K, N.SEEDS_E[:1], 2e3))
    for v in N.VARIANTS_K:
        R = np.asarray(raw['E_rec'][v])[0]; r = q['var'][v]; n = r['W'].shape[1]
        assert np.allclose(r['W'][0], R[:n, 4:7], rtol=0, atol=1e-14) and np.allclose(r['al'][0], R[:n, 1:4], rtol=0, atol=1e-13)
        assert np.allclose(r['acc'][0], K['var'][v]['acc'][0, :n], rtol=1e-9, atol=1e-15)
        assert np.allclose(r['Th'][0], K['var'][v]['Th'][0, :n], rtol=0, atol=1e-13)
        assert np.abs(N.k_identity(r)).max() < 1e-13
        assert (v == 'nonlinear only') == np.all(r['acc'][:, :, 4] == 0)
    from src.core import gauge_batch as GB
    loop, kw, Th = N.e_loop('nonlinear only'); bd = GB.BatchDither(loop)
    T3 = (Th + 0.05 * np.array([[1.0, -0.5, 0.3], [-0.4, 0.2, 0.6], [0.7, 0.1, -0.2]]))[None]; b3 = np.array([[0.9, 1.05, 1.05]])
    d3 = bd.drift(T3, b3)[0]; a, u, g, gp, nonfit, boxt = N.gauge_terms(bd, T3, b3, d3)
    kk, c = 2.0 / N.A_MAX, N.A_MAX ** 2 / (4 * E.G[0]); al = T3[0, :, 0]; bA = d3[0, :, 0]
    phi = lambda x: c * np.cosh(kk * x)
    ch = -(phi(al + N.DT * bA) - phi(al) - N.DT * c * kk * np.sinh(kk * al) * bA) / b3[0]
    assert np.allclose(u[0] * N.DT - (phi(al + N.DT * bA) - phi(al)) / b3[0] - g[0] * nonfit[0] * N.DT, ch, rtol=1e-7, atol=1e-17)
    assert np.allclose(u[0], loop.fitness(loop.residual(T3[0], b3[0]), T3[0]), rtol=1e-12, atol=1e-18)
    assert np.all(ch < 0) and np.abs(ch).max() > 1e-9


def test_part_K_the_offset_and_the_response():
    """k_offset vanishes at an exact fit and is Proposition D‴.40's bounded term about it; k_response accumulates a constant fitness where
    the slow limit has no drift (γW∘(du − W·du)·t at the masses 1/3) and settles at −J⁻¹γW∘(du − W·du) where it has one."""
    N = _N(); net, base = E.setup()
    for W in (np.ones(3) / 3, np.array([0.40, 0.32, 0.28])):
        al = E.exact_rep(net, W)[:, 0]
        assert np.abs(N.k_offset(al, W)).max() < 1e-14
        al2 = al + np.array([0.1, -0.05, 0.02]); c = N.A_MAX ** 2 / (4 * E.G[0])
        d = c * (np.cosh(2 * al2 / N.A_MAX) - np.cosh(2 * al / N.A_MAX))
        assert np.allclose(N.k_offset(al2, W), N.GAMMA_E / 3 * (d - W * d.sum()), rtol=1e-12)
    s = np.arange(0, 5e4 + 1, 1e3); du = np.tile(np.array([1e-6, -2e-6, 0.5e-6]), (2, len(s) - 1, 1))
    x = N.k_response('nonlinear only', s, du)
    W = np.ones(3) / 3; pred = N.GAMMA_E * W * (du[0, 0] - du[0, 0] @ W) * s[-1]
    assert np.allclose(x[0, -1], pred, rtol=1e-9) and np.allclose(x[1], x[0])
    s = np.arange(0, 2e6 + 1, 1e4); du = np.tile(np.array([1e-6, -2e-6, 0.5e-6]), (1, len(s) - 1, 1))
    x = N.k_response('nonlinear only + barrier', s, du)[0, -1]
    co = N.lead_coeffs('nonlinear only + barrier'); f = lambda W: N.GAMMA_E * N.GN.replicator_flow(W, N.lead_fitness('nonlinear only + barrier', W, co))
    from scipy.optimize import fsolve
    Wz = fsolve(lambda w: f(np.r_[w, 1 - w.sum()])[:2], [0.37, 0.30], xtol=1e-13); Wz = np.r_[Wz, 1 - Wz.sum()]
    e = 1e-6; Jc = np.array([(f(Wz + e * v) - f(Wz - e * v))[:2] / (2 * e) for v in (np.array([1.0, 0, -1]), np.array([0, 1.0, -1]))]).T
    st = -np.linalg.solve(Jc, (N.GAMMA_E * Wz * (du[0, 0] - du[0, 0] @ Wz))[:2])
    assert np.allclose(x[:2], st, rtol=0.02) and abs(x.sum()) < 1e-15


def test_part_K_numbers_recomputed_from_the_raw():
    N = _N(); raw = _raw(); K = raw['K']; F = raw['F']
    assert K['replay'] < 1e-12 and K['S'] == N.S_E and K['seeds'] == N.SEEDS_E
    md = open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8').read()
    pe = lambda M: ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(M.mean(axis=0), M.std(axis=0, ddof=1) / np.sqrt(len(M))))
    for v in N.VARIANTS_K:
        R = np.asarray(raw['E_rec'][v]); hm = R[:, 500:, 4:7].mean(axis=1); sh = F['var'][v]['sde']['half']
        x = N.k_effects(K, raw['E_rec'], F, v)
        assert np.allclose(x['obs'], hm - sh.mean(axis=0), rtol=0, atol=1e-15)
        z0 = np.abs(x['obs'].mean(axis=0)) / N.zden(x['obs'].std(axis=0, ddof=1) / np.sqrt(12), sh.std(axis=0, ddof=1), 12)
        assert f"| {v} | {pe(x['obs'])} | {z0.max():.1f} | {pe(x['off'])} |" in md
        assert abs(z0.max() - {'nonlinear only': 21.8, 'nonlinear only + barrier': 8.3}[v]) < 0.05           # Part F's |z|
    assert "## K. The Euler chain's bias along Part E" in md


def test_part_L_errors_and_extremes_recomputed():
    N = _N(); raw = _raw(); H, F, I, C = raw['H'], raw['F'], raw['I'], raw['C']
    md = open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8').read()
    mc = {r['tc']: r for r in H['mc']}; c, w = mc[N.TC_H]['u'] - mc[N.TC_H]['box'], mc[0.0]['u'] - mc[0.0]['box']
    d = c.mean(axis=0) / w.mean(axis=0)
    var = d ** 2 * (np.var(c, axis=0, ddof=1) / len(c) / c.mean(axis=0) ** 2 + np.var(w, axis=0, ddof=1) / len(w) / w.mean(axis=0) ** 2)
    assert ', '.join(f'{a:.3f} ± {b:.3f}' for a, b in zip(d, np.sqrt(var))) in md
    zr = [z for z in F['zeros'] if z['variant'] == 'correlated E + class levels' and z['gamma'] == 1e-2][0]['zero']
    assert ', '.join(f'{x:+.4f}' for x in zr - N.w_sharp_weighted(N.levels(N.noises()['correlated E']))) in md
    assert f"at most {max(np.abs(r['centres and widths only']['D']).max() for r in C):.1e} over all γ and W" in md
    smp = N.i_samples(I); lab = f"the {len(N.SEEDS_I)} new seeds ({N.SEEDS_I[0]}–{N.SEEDS_I[-1]})"
    Rr, acc = smp[lab]; s = I['runs']['old']['s']; i = int(np.argmin(np.abs(s - 2e4)))
    y = Rr[:, i, 4]; m = N.i_contrib(acc[:, i])[:, 1]; X = np.c_[np.ones_like(m), m]
    coef, rss = np.linalg.lstsq(X, y, rcond=None)[:2]
    se_int = np.sqrt(rss[0] / (len(y) - 2) * np.linalg.inv(X.T @ X)[0, 0])                    # OLS: the intercept's s.e.
    assert abs(coef[0] - N.i_regressed(y, m, np.zeros(2))[0]) < 1e-14
    assert f"| {lab} | 20000 | {N.i_regressed(y, m, np.zeros(2))[1]:.5f} | {se_int:.5f} |" in md
    assert "## L. Errors and extremes" in md


def test_part_K_the_halved_step_follows_part_Es_brownian_path():
    """k_half against an independent integration with ReplicatorCohesionLoop.drift: two half steps per full step, the increments
    sq·(z ± w)/2 (z run_dithered's draw for Part E's seed, w from the seed + W_OFF_KH generator), so that their sum is Part E's full
    increment; the decomposition exact; the raw's half run reproduced at its first records, its chain term about half the Δt run's."""
    N = _N(); rec = N.REC
    try:
        N.REC = 1
        q = N.k_half.__wrapped__((N.VARIANTS_K, N.SEEDS_E[:1], 3 * N.DT))
    finally:
        N.REC = rec
    for v in N.VARIANTS_K:
        loop, kw, Th = N.e_loop(v); b = np.ones(3); lb = np.zeros(3)
        sq = np.sqrt(2 * N.TAU * N.DT * np.asarray(loop.G)) * np.asarray(kw['mask'])
        g, gw = np.random.default_rng(N.SEEDS_E[0]), np.random.default_rng(N.SEEDS_E[0] + N.W_OFF_KH)
        for k in range(3):
            z, w = g.standard_normal((3, 3)), gw.standard_normal((3, 3))
            assert np.allclose(0.5 * sq * (z + w) + 0.5 * sq * (z - w), sq * z, rtol=0, atol=1e-15)
            for dz in (z + w, z - w):
                b = np.exp(lb); dT, db = loop.drift(Th, b)
                Th = Th + N.DT / 2 * dT + 0.5 * sq * dz
                lb = lb + N.DT / 2 * db / b; lb = lb - np.log(np.exp(lb).sum() / 3)
            r = q['var'][v]
            assert np.allclose(r['al'][0, k + 1], Th[:, 0], rtol=0, atol=1e-13) and np.allclose(r['W'][0, k + 1], np.exp(lb) / 3, rtol=0, atol=1e-15)
        assert np.abs(N.k_identity(r)).max() < 1e-13
    raw = _raw(); H = raw['K']['half']
    assert H['S'] == N.S_KH and H['seeds'] == N.SEEDS_E
    q = N.k_half.__wrapped__((N.VARIANTS_K, N.SEEDS_E[:1], 2e3))
    for v in N.VARIANTS_K:
        assert np.allclose(q['var'][v]['W'][0], H['var'][v]['W'][0, :3], rtol=0, atol=1e-14)
        assert np.allclose(q['var'][v]['acc'][0], H['var'][v]['acc'][0, :3], rtol=1e-10, atol=1e-16)        # every term, not only their sum
        c_full, c_half = raw['K']['var'][v]['acc'][:, 200, 3], H['var'][v]['acc'][:, 200, 3]
        assert np.all(np.abs(c_half / c_full - 0.5) < 0.05)
        assert np.abs(N.k_identity(H['var'][v])).max() < 1e-12


def test_part_K_the_step_effect_recomputed():
    """k_forcing is 2(acc(Δt) − acc(Δt/2)) summed over the identity's terms but the masses form's, fitted through the origin; its response
    reproduces the paired runs' 2(W(Δt) − W(Δt/2)) at 2·10⁵ to 15%; the md's rows recomputed."""
    N = _N(); raw = _raw(); K = raw['K']; F = raw['F']
    md = open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8').read()
    pe = lambda M: ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(M.mean(axis=0), M.std(axis=0, ddof=1) / np.sqrt(len(M))))
    for v in N.VARIANTS_K:
        s = K['half']['s']; n = len(s)
        f = 2 * (K['var'][v]['acc'][:, :n, :4] - K['half']['var'][v]['acc'][:, :, :4]).sum(axis=2)
        sl = np.array([[np.polyfit(np.r_[s, -s], np.r_[f[r, :, k], -f[r, :, k]], 1)[0] for k in range(3)] for r in range(len(f))])
        assert np.allclose(N.k_forcing(K, v), sl, rtol=1e-6, atol=1e-16)
        x = N.k_effects(K, raw['E_rec'], F, v)
        m, r = x['at_half'].mean(axis=0), x['resp_half'].mean(axis=0)
        assert np.all(np.abs(r / m - 1) < 0.15), (v, m, r)
        a, b = x['obs'] - x['eff'], x['obs'] - x['off'] - x['eff']
        assert f"| {pe(x['eff'])} | {pe(a)} |" in md and f"| {pe(b)} |" in md
    x = N.k_effects(K, raw['E_rec'], F, 'nonlinear only')
    own = x['obs'] + F['var']['nonlinear only']['sde']['half'].mean(axis=0) - 1 / 3 - x['off'] - x['eff']
    assert np.allclose(own.mean(axis=0), (np.asarray(raw['E_rec']['nonlinear only'])[:, 500:, 4:7].mean(axis=1) - 1 / 3 - x['off'] - x['eff']).mean(axis=0), rtol=1e-9)
    assert f"= {pe(own)}, is the loop's own secular displacement" in md


# ---- Part M: the missing O(γ²τ) drift of the deterministic-amplitude loop ----

def test_part_M_the_closed_form_of_the_antisymmetric_correlation():
    """m_ou against an independent computation: the Lyapunov equation solved as a linear system (Kronecker), the Euler chain's law as
    the doubled sum Σ_j(I + ΔtB)ʲΔtΣ(I + ΔtB)ʲᵀ; α carries no noise, so (CBᵀ)_αα is antisymmetric; at τ and 3τ the drift per γ²τ is the same
    (the linearised law scales with τ); it is Part C's `limit_drift` at γ = 10⁻⁵ to 3%."""
    N = _N(); o = N.m_ou(N.TAU); o3 = N.m_ou(3 * N.TAU)
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    th0 = E.exact_rep(net, np.ones(3) / 3).ravel(); f = lambda t: loop.drift(t.reshape(3, 3), np.ones(3))[0].ravel()
    Bm = np.array([(f(th0 + 1e-5 * e) - f(th0 - 1e-5 * e)) / 2e-5 for e in np.eye(9)]).T
    Sig = N.noises(N.TAU)['centres and widths only']; I9 = np.eye(9)
    C = np.linalg.solve(np.kron(I9, Bm) + np.kron(Bm, I9), -Sig.ravel()).reshape(9, 9)
    X = C @ Bm.T
    assert np.allclose(X[0::3, 0::3], -X[0::3, 0::3].T, rtol=0, atol=1e-12 * np.abs(X).max())
    dphi = np.sinh(th0[0::3] / 2) * 2.0                                                       # (a_max/2G_α)sinh(2α/a_max), a_max = 4, G_α = 1
    d = (dphi[:, None] * dphi[None, :] * X[0::3, 0::3]).sum(axis=0) / (9 * N.TAU)
    assert np.allclose(o['d'], d, rtol=1e-6, atol=0) and np.allclose(o3['d'], o['d'], rtol=1e-6, atol=0)
    Cd = N.DT * Sig; A = I9 + N.DT * Bm
    for _ in range(40):                                                                     # doubling: Σ_{j<2ⁿ} AʲΔtΣAʲᵀ
        Cd = Cd + A @ Cd @ A.T; A = A @ A
    Xd = dphi[:, None] * dphi[None, :] * (Cd @ Bm.T)[0::3, 0::3]
    assert np.allclose(o['chain'][0], 0.5 * (Xd - Xd.T).sum(axis=0) / (9 * N.TAU), rtol=1e-4, atol=0)
    assert np.allclose(o['chain_sym'][0], 0.5 * (Xd + Xd.T).sum(axis=0) / (9 * N.TAU), rtol=1e-4, atol=0)
    c5 = [r for r in _raw()['C'] if r['gamma'] == 1e-5 and np.allclose(r['w'], 1 / 3)][0]['centres and widths only']['b'] / 1e-5
    assert np.all(np.abs(c5 / o['d'] - 1) < 0.05) and np.abs(c5 - o['d']).max() < 4e-4


def test_part_M_the_limit_map_and_the_offset_with_the_fits_motion():
    """m_limit leaves an exact fit where it is, agrees with gauge_drift.limit_map (LSODA on Katzenberger's flow, no box) from a displaced
    state, and W − Y is (I − M)⁻¹ times Part K's first-order offset up to O(γ²δ²) — the first-order offset itself misses it by far more."""
    from src.core import gauge_drift as GDR
    import exp_gauge_drift as X
    N = _N(); net, _ = E.setup()
    W0 = np.array([0.40, 0.32, 0.28]); Th0 = E.exact_rep(net, W0)
    d = np.array([[0.04, -0.03, 0.10], [-0.03, 0.02, -0.15], [0.05, 0.01, 0.08]])
    Y, chk = N.m_limit(np.array([Th0, Th0 + d]), np.log(3 * np.array([W0, W0])), N.GAMMA_E)
    assert np.abs(Y[0] - W0).max() < 1e-13 and np.abs(Y - chk).max() < 1e-8
    F, chart, coords, _ = X.instance(N.GAMMA_E)
    y, _ = GDR.limit_map(F, np.r_[(Th0 + d).ravel(), np.log(W0[:2] / W0[2])], 4e4)
    e = np.exp(np.r_[y[9:], 0.0])
    assert np.abs(Y[1] - e / e.sum()).max() < 1e-9
    lin, Mm = N.m_offset_lin((Th0 + d)[:, 0], W0)
    off = N.k_offset((Th0 + d)[:, 0], W0)
    a = 1 / W0 * np.abs(E.A_BUMP); h = a / (1 - (a / 4.0) ** 2)
    assert np.allclose(np.diag(Mm), N.GAMMA_E / 3 * (1 - W0) * (-h ** 2 / W0), rtol=1e-12)
    assert np.abs(W0 - Y[1] - lin).max() < 0.1 * np.abs(W0 - Y[1] - off).max()
    assert np.abs(lin - np.linalg.solve(np.eye(3) - Mm, off)).max() < 1e-18


def test_part_M_the_coupled_loop_is_part_Ks_two_chains():
    """m_loop's Δt chain is Part E's run (Part K's replay) and its Δt/2 chain is k_half's, for the same seed; the states' limit maps are
    m_limit's; the raw's runs carry the stated arguments."""
    N = _N(); raw = _raw(); K = raw['K']; ty = N.T_YM
    try:
        N.T_YM = (1000.0, 2000.0)
        q = N.m_loop.__wrapped__((N.GAMMA_E, (N.SEEDS_E[0],), 2e3))
    finally:
        N.T_YM = ty
    r, hv = K['var']['nonlinear only'], K['half']['var']['nonlinear only']
    assert np.array_equal(q['W'][0], r['W'][0, :3]) and np.array_equal(q['al'][0], r['al'][0, :3])
    assert np.array_equal(q['W'][1], hv['W'][0, :3]) and np.array_equal(q['al'][1], hv['al'][0, :3])
    Th = r['Th'][0, 2]; Y = N.m_limit(Th[None], np.log(3 * r['W'][0, 2])[None], N.GAMMA_E)[0]
    assert np.allclose(q['Y'][1, 0], Y[0], rtol=0, atol=1e-15)
    M = raw['M']
    for g in N.GAMMAS_M:
        ch = [M['loop'][g][c] for c in sorted(M['loop'][g])]
        assert sum((c['seeds'] for c in ch), ()) == N.SEEDS_M and all(c['S'] == N.S_M and c['gamma'] == g for c in ch)
        assert all(tuple(c['t']) == N.T_YM for c in ch)


def test_part_M_the_fixed_masses_correlation_step_by_step():
    """m_fixed over a few full steps against an independent integration with ReplicatorCohesionLoop.drift (γ = 0, the box): the two coupled
    chains and the sums of φ_l, φ̇_k = φ′_kb_α_k and φ_lφ̇_k; m_anti's split."""
    import exp_gauge_plateau as P
    N = _N(); rec = N.REC
    try:
        N.REC = 2
        q = N.m_fixed.__wrapped__((3e-4, (900,), 4 * N.DT, 2 * N.DT))
    finally:
        N.REC = rec
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B); b = np.ones(3)
    sq = np.sqrt(2 * 3e-4 * N.DT * E.G) * np.array([0.0, 1.0, 1.0])
    g, gw = np.random.default_rng(900), np.random.default_rng(900 + N.W_OFF_KH)
    Tf = E.exact_rep(net, np.ones(3) / 3).copy(); Th = Tf.copy()
    sums = {c: dict(n=0, p=0, pd=0, ppd=0) for c in ('full', 'half')}
    ph = lambda T, dT: (4.0 * np.cosh(T[:, 0] / 2), 2.0 * np.sinh(T[:, 0] / 2) * dT[:, 0])
    Z, Zw = g.standard_normal((2, 3, 3)), gw.standard_normal((2, 3, 3))
    Z, Zw = np.r_[Z, g.standard_normal((2, 3, 3))], np.r_[Zw, gw.standard_normal((2, 3, 3))]
    Z, Zw = np.r_[Z, g.standard_normal((2, 3, 3))], np.r_[Zw, gw.standard_normal((2, 3, 3))]
    for k in range(6):
        z, w = Z[k], Zw[k]
        dT = loop.drift(Tf, b)[0]
        if k >= 2:
            p, pd = ph(Tf, dT); s = sums['full']; s['n'] += 1; s['p'] += p; s['pd'] += pd; s['ppd'] += p[:, None] * pd[None, :]
        Tf = Tf + N.DT * dT + sq * z
        for dz in (z + w, z - w):
            dT = loop.drift(Th, b)[0]
            if k >= 2:
                p, pd = ph(Th, dT); s = sums['half']; s['n'] += 1; s['p'] += p; s['pd'] += pd; s['ppd'] += p[:, None] * pd[None, :]
            Th = Th + N.DT / 2 * dT + 0.5 * sq * dz
    for c in ('full', 'half'):
        for k_ in ('n', 'p', 'pd', 'ppd'):
            assert np.allclose(q[c][k_][0], sums[c][k_], rtol=1e-10, atol=1e-15), (c, k_)
    a, s = N.m_anti(q['full'], 3e-4)
    n = q['full']['n'][0]; cc = q['full']['ppd'][0] / n - np.outer(q['full']['p'][0] / n, q['full']['pd'][0] / n)
    assert np.allclose(a[0], (cc - cc.T).sum(axis=0) / 2 / (9 * 3e-4), rtol=1e-12) and np.allclose(a[0] + s[0], cc.sum(axis=0) / (9 * 3e-4), rtol=1e-12)


def test_part_M_numbers_recomputed_from_the_raw():
    """From the raw: Part K's records' limit map reproduced at two records; the rates, Richardson's extrapolation, the z against
    Katzenberger's drift and the fixed-mass correlation recomputed and found in the md; Φ's second differences agree with `limit_drift`."""
    N = _N(); raw = _raw(); M = raw['M']; K = raw['K']
    md = open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8').read()
    assert "## M. The missing O(γ²τ) drift" in md and md.index("## L. Errors and extremes") < md.index("## M. The missing O(γ²τ) drift")
    W, al, Y, chk, idx = N.m_records(M, K)
    assert W.shape == (len(N.SEEDS_E), 101, 3) and idx[0] == 500 and idx[-1] == 1000 and np.all(np.diff(idx) == N.EVERY_M)
    ev, ck = N.EVERY_M, N.CHUNK_MK
    try:
        N.EVERY_M, N.CHUNK_MK = 500, 1
        q = N.m_lim_K(K, 5)
    finally:
        N.EVERY_M, N.CHUNK_MK = ev, ck
    assert list(q['idx']) == [500, 1000] and np.allclose(q['Y'][0], Y[5, [0, 100]], rtol=0, atol=1e-14)
    assert np.allclose(q['chk'][0], chk[5, [0, 100]], rtol=0, atol=1e-14)
    se = lambda X: X.std(axis=0, ddof=1) / np.sqrt(len(X))
    pe = lambda X: ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(X.mean(axis=0), se(X)))
    eff = N.k_effects(K, raw['E_rec'], raw['F'], 'nonlinear only')['eff']
    yd = (Y - 1 / 3).mean(axis=1) - eff
    assert f"| Y − 1/3 − (2) | {pe(yd)} |" in md and f"| W − Y, Y the exact limit map | {pe((W - Y).mean(axis=1))} |" in md
    for g in N.GAMMAS_M:
        ch = [M['loop'][g][c] for c in sorted(M['loop'][g])]
        Yl = np.concatenate([c['Y'].reshape(3, 2, -1, 3) for c in ch], axis=2)
        rt = (Yl[2] - Yl[0]) / (N.T_YM[2] - N.T_YM[0]) / (g * g * N.TAU); ri = 2 * rt[1] - rt[0]
        kz = M['katz'][g]; z = (ri.mean(axis=0) - kz) / se(ri)
        p5 = lambda X: ', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(X.mean(axis=0), se(X)))
        assert f"| {g:g} | {p5(rt[0])} | {p5(rt[1])} | {p5(ri)} | {', '.join(f'{x:+.5f}' for x in kz)} | {', '.join(f'{v:+.1f}' for v in z)} |" in md
        b = N.k_drift(N.katz('nonlinear only', g), np.ones(3) / 3)[0] if g == N.GAMMA_E else None
        if b is not None:
            assert np.allclose(np.r_[b, -b.sum()] / (g * g * N.TAU), kz, rtol=1e-9)
    for t in N.TAUS_MF:
        ch = [M['fixed'][t][c] for c in sorted(M['fixed'][t])]
        S, ns = N.S_MF[t]
        assert sum((c['seeds'] for c in ch), ()) == N.SEEDS_MF[:ns] and all(c['S'] == S and c['burn'] == N.BURN_MF for c in ch)
        cat = lambda c_: {k: np.concatenate([x[c_][k] for x in ch]) for k in ('n', 'p', 'pd', 'ppd')}
        out = []
        for c_ in ('full', 'half'):
            a = cat(c_); n = a['n'][:, None]
            cc = a['ppd'] / n[:, :, None] - np.einsum('rl,rk->rlk', a['p'] / n, a['pd'] / n)
            out.append(np.einsum('rlk->rk', cc - cc.transpose(0, 2, 1)) / 2 / (9 * t))
        ri = 2 * out[1] - out[0]
        assert f"| {t:g} | {ns} × {S:g} | " in md and ', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(ri.mean(axis=0), se(ri))) in md
    for r_ in M['fd']:
        assert np.allclose(r_['rich'], (4 * r_['fd'][1] - r_['fd'][0]) / 3, rtol=1e-12)
        # at γ = 10⁻³ the masses' second-order response to h = 0.025 is ~10⁻¹² and the limit map's own error shows (class 3: 2.4%)
        assert np.abs(r_['rich'] / r_['kz'] - 1).max() < (0.002 if r_['gamma'] == 1e-2 else 0.03) and np.abs(r_['at_fit']).max() < 1e-13


def test_part_M_of_the_md_regenerates_from_the_raw():
    N = _N(); raw = _raw()
    md = N.report(*[raw.get(k) for k in 'ABCDE'], **{k: raw.get(k) for k in ('E_rec', 'F', 'G', 'H', 'I', 'J', 'K', 'M', 'N')})
    with open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8') as f:
        assert f.read() == md


def test_part_M_the_shortfall_split_and_the_comparisons_recomputed():
    """Part M(e) and the lines after (c) and (d), recomputed from the raw: Katzenberger's path split into Part K's secular displacement,
    the first-order offset's error, the step correction's excess over (c)'s step effect and the rest (summing to the path); the first
    transient minus 1/3; (d)'s z against the closed form; (c) against (d); the rates' γ dependence."""
    N = _N(); raw = _raw(); M = raw['M']; K = raw['K']
    md = open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8').read()
    se = lambda X: X.std(axis=0, ddof=1) / np.sqrt(len(X))
    W, al, Y, chk, idx = N.m_records(M, K)
    x = N.k_effects(K, raw['E_rec'], raw['F'], 'nonlinear only'); kp = raw['F']['var']['nonlinear only']['kpath']['half'][0] - 1 / 3
    off = N.k_offset(al, W); hs = K['s'][len(K['s']) // 2:].mean()
    assert hs == 750000.0 and len(W) == 12
    rates = {}
    for g in N.GAMMAS_M:
        ch = [M['loop'][g][c] for c in sorted(M['loop'][g])]
        Yl = np.concatenate([c['Y'].reshape(3, 2, -1, 3) for c in ch], axis=2)
        rt = (Yl[2] - Yl[0]) / (N.T_YM[2] - N.T_YM[0]) / (g * g * N.TAU)
        rates[g] = (rt[0], 2 * rt[1] - rt[0])
        y0 = 2 * Yl[0, 1] - Yl[0, 0] - 1 / 3
        assert f"γ = {g:g} " + ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(y0.mean(axis=0), se(y0))) in md
        assert np.all(np.abs(y0.mean(axis=0)) < 1e-5)
    full, rich = rates[1e-2]; ec = (full - rich) * hs * 1e-4 * N.TAU; assert len(ec) == 48
    kres = (W - 1 / 3 - off).mean(axis=1) - x['eff']; err1 = (off - (W - Y)).mean(axis=1)
    y2 = (Y - 1 / 3).mean(axis=1) - ec.mean(axis=0)
    terms = [(kres.mean(axis=0), se(kres)), (err1.mean(axis=0), se(err1)),
             (x['eff'].mean(axis=0) - ec.mean(axis=0), np.hypot(se(x['eff']), se(ec))), (kp - y2.mean(axis=0), np.hypot(se(y2), se(ec)))]
    assert np.allclose(sum(t[0] for t in terms), kp, rtol=1e-12, atol=0)
    for m_, e_ in terms:
        row = f"| {', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(m_, e_))} | {', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(m_ / kp, e_ / np.abs(kp)))} |"
        assert row in md, row
    assert md.index("### (e) Part K's shortfall") > md.index("### (d) The positions' antisymmetric")
    r2, e2 = y2.mean(axis=0) / kp, np.hypot(se(y2), se(ec)) / np.abs(kp)
    kz = M['katz'][1e-2]; rc, ecr = rich.mean(axis=0) / kz, se(rich) / np.abs(kz)
    assert f"over Katzenberger's path {', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(r2, e2))}; " in md
    assert f"at γ = 0.01, z: {', '.join(f'{v:+.1f}' for v in (r2 - rc) / np.hypot(e2, ecr))}." in md
    fx = {}
    for t in N.TAUS_MF:
        ch = [M['fixed'][t][c] for c in sorted(M['fixed'][t])]
        out = []
        for c_ in ('full', 'half'):
            a = {k: np.concatenate([x_[c_][k] for x_ in ch]) for k in ('n', 'p', 'pd', 'ppd')}; n = a['n'][:, None]
            cc = a['ppd'] / n[:, :, None] - np.einsum('rl,rk->rlk', a['p'] / n, a['pd'] / n)
            out.append(np.einsum('rlk->rk', cc - cc.transpose(0, 2, 1)) / 2 / (9 * t))
        fx[t] = 2 * out[1] - out[0]
    d = {o['tau']: o['d'] for o in M['ou']}
    zs = '; '.join(f"τ = {t:g} {', '.join(f'{v:+.1f}' for v in (fx[t].mean(axis=0) - d[t]) / se(fx[t]))}" for t in N.TAUS_MF)
    assert f"Richardson against the closed form, z: {zs}. " in md
    rd, erd = fx[N.TAU].mean(axis=0) / d[N.TAU], se(fx[N.TAU]) / np.abs(d[N.TAU])
    for g in N.GAMMAS_M:
        kz = M['katz'][g]; ri = rates[g][1]
        assert f"γ = {g:g} {', '.join(f'{v:+.1f}' for v in (ri.mean(axis=0) / kz - rd) / np.hypot(se(ri) / np.abs(kz), erd))}" in md
    a, b = rates[3e-3][1], rates[1e-2][1]; r = a.mean(axis=0) / b.mean(axis=0)
    er = (a - r * b).std(axis=0, ddof=1) / np.sqrt(len(a)) / np.abs(b.mean(axis=0))
    assert f"Richardson {', '.join(f'{u:.2f} ± {v:.2f}' for u, v in zip(r, er))}, Katzenberger's " in md


def test_part_M_the_second_differences_and_richardson_on_a_known_map():
    """m_fd with the limit map replaced by a known quartic of the state, Y_k = 1/3 + Σ_i(c_ik δ_i² + e_ik δ_i⁴): the second differences are
    2c_ik + 2e_ik h², Richardson's (4d(h/2) − d(h))/3 removes the h² term exactly, and the drift is ½Σ_iΣ_ii·2c_ik/(γ²τ)."""
    N = _N(); net, _ = E.setup()
    x0 = E.exact_rep(net, np.ones(3) / 3).ravel(); Sig = N.noises()['centres and widths only']
    rng = np.random.default_rng(3); c = rng.normal(size=(9, 3)) * 1e-6; e = rng.normal(size=(9, 3)) * 1e-4
    def fake(Th, lb, gamma):
        d = np.asarray(Th).reshape(len(Th), 9) - x0
        Y = 1 / 3 + (d ** 2) @ c + (d ** 4) @ e
        return Y, Y
    ml = N.m_limit
    try:
        N.m_limit = fake
        r = N.m_fd.__wrapped__(1e-2)
    finally:
        N.m_limit = ml
    want = 0.5 * np.einsum('i,ik->k', np.diag(Sig), 2 * c) / (1e-4 * N.TAU)
    assert np.allclose(r['rich'], want, rtol=1e-6, atol=0)
    assert not np.allclose(r['fd'][0], want, rtol=1e-3, atol=0)


# ---- the finite-τ question (D‴.40 after Part M): Part N ----
def test_part_N_the_design_the_seeds_new_and_the_raw_as_designed():
    """Part N's seeds and run lengths as stated: fresh seeds (none used in Parts A–M, nor their bridges' seed + W_OFF_KH), Part M(d)'s run
    lengths at 10⁻⁴ and 3·10⁻⁴ (so that the seeds pool), and the raw's chunks as designed (seeds in order, τ, S, the burn-in)."""
    N = _N(); raw = _raw(); Nn = raw['N']
    assert N.SEEDS_N == {1e-4: tuple(range(1100, 1132)), 1e-3: tuple(range(1300, 1316)), 3e-4: tuple(range(1400, 1544))}
    assert list(N.SEEDS_N) == [1e-4, 1e-3, 3e-4] and N.CHUNK_N == 16 and N.TAUS_FIT_N == (3e-5, 1e-4, 3e-4)
    assert N.S_N == {1e-4: 2e5, 1e-3: 1e5, 3e-4: 1e5} and all(N.S_N[t] == N.S_MF[t][0] for t in N.S_N if t in N.S_MF)
    old = set(N.SEEDS_MF) | set(N.SEEDS_M) | set(N.SEEDS_E) | set(N.SEEDS_G) | set(N.REPLAY_G) | set(N.SEEDS_I) | set(N.NOISE_SEEDS)
    old |= {N.SEED_D, N.SEED_J, N.SEED_H, N.SEED_F, N.SEED_SDE_I, 21} | set(range(7000, 7000 + len(N.W1_FIX_I) * N.SEEDS_FIX_I))
    old |= {s + N.BURN_OFF_I for s in N.SEEDS_I}
    new = [s for sd in N.SEEDS_N.values() for s in sd]
    assert len(set(new)) == len(new) and not (set(new) & old) and not ({s + N.W_OFF_KH for s in new} & (old | {s + N.W_OFF_KH for s in old}))
    assert set(Nn) == {'fixed'} and set(Nn['fixed']) == set(N.SEEDS_N)
    for t, sd in N.SEEDS_N.items():
        ch = [Nn['fixed'][t][c] for c in sorted(Nn['fixed'][t])]
        assert sorted(Nn['fixed'][t]) == list(range(0, len(sd), N.CHUNK_N)) and sum((c['seeds'] for c in ch), ()) == sd
        assert all(c['tau'] == t and c['S'] == N.S_N[t] and c['burn'] == N.BURN_MF for c in ch)
        n = round((N.S_N[t]) / N.DT)
        assert all(np.all(c['full']['n'] == n) and np.all(c['half']['n'] == 2 * n) for c in ch)
    assert raw['M']['fixed'].keys() == set(N.TAUS_MF)                 # Part M's key untouched


def test_part_N_resumes_and_calls_m_fixed_as_designed():
    """part_N with m_fixed replaced by a recorder: from a raw holding some chunks it computes only the missing ones, with m_fixed's
    arguments (τ, the chunk's seeds, S_N[τ], BURN_MF), in SEEDS_N's order, saving after each; the stored chunks are kept."""
    N = _N(); calls, saves = [], []
    def fake(args):
        calls.append(args); return dict(tau=args[0], seeds=args[1])
    mf = N.m_fixed
    try:
        N.m_fixed = fake
        res = {'N': {'fixed': {3e-4: {16: 'kept'}}}, 'M': 'm'}
        out = N.part_N(res, lambda: saves.append(len(calls)))
    finally:
        N.m_fixed = mf
    want = [(t, sd[c:c + 16], N.S_N[t], N.BURN_MF) for t, sd in N.SEEDS_N.items() for c in range(0, len(sd), 16) if not (t == 3e-4 and c == 16)]
    assert calls == want and saves == list(range(1, len(want) + 1)) and out['fixed'][3e-4][16] == 'kept' and res['N'] is out
    assert all(len(a[1]) == 16 for a in calls) and res['M'] == 'm'


def test_part_N_the_weighted_line_against_numpy():
    """n_fit against numpy's weighted polyfit (its unscaled covariance) and an exact line (χ² zero), and the constant's weighted mean."""
    N = _N(); rng = np.random.default_rng(5)
    x = np.array([3e-5, 1e-4, 3e-4, 1e-3]); e = rng.uniform(0.05, 0.5, size=(4, 3)); y = rng.normal(1, 0.3, size=(4, 3))
    f = N.n_fit(x, y, e)
    for k in range(3):
        p, cov = np.polyfit(x, y[:, k], 1, w=1 / e[:, k], cov='unscaled')
        assert np.allclose([f['b'][k], f['c'][k]], p, rtol=1e-9, atol=0)
        assert np.allclose([f['eb'][k] ** 2, f['ec'][k] ** 2, f['cov'][k]], [cov[0, 0], cov[1, 1], cov[0, 1]], rtol=1e-8, atol=0)
        w = 1 / e[:, k] ** 2; m = (w * y[:, k]).sum() / w.sum()
        assert np.isclose(f['m'][k], m, rtol=1e-12, atol=0) and np.isclose(f['chim'][k], (w * (y[:, k] - m) ** 2).sum(), rtol=1e-12, atol=0)
        assert np.isclose(f['chi'][k], (w * (y[:, k] - np.polyval(p, x)) ** 2).sum(), rtol=1e-9, atol=0)
    g = N.n_fit(x, 0.9 - 200.0 * x[:, None] + 0 * e, e)
    assert np.allclose(g['c'], 0.9, rtol=1e-12) and np.allclose(g['b'], -200.0, rtol=1e-9) and np.all(g['chi'] < 1e-20) and f['dof'] == 2


def test_part_N_numbers_recomputed_from_the_raw():
    """From the raw, independently of m_fixed_all and n_fit: the per-seed Richardson ratio over the closed form, Part M(d)'s seeds and
    Part N's, pooled; the line through the pooled points (numpy's polyfit) and the constant; each row of the md."""
    from scipy.stats import chi2
    N = _N(); raw = _raw(); M = raw['M']; Nn = raw['N']
    md = open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8').read()
    ou = {o['tau']: o for o in M['ou']}; d = ou[N.TAU]['d']
    assert all(np.allclose(o['d'], d, rtol=1e-11, atol=0) for o in M['ou'])            # the closed form per τ is τ-independent
    se = lambda X: X.std(axis=0, ddof=1) / np.sqrt(len(X))
    r2 = lambda m, e: ', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(m, e))
    z1 = lambda v: ', '.join(f'{x:+.1f}' for x in v)

    def ratio(chunks, t):
        out = []
        for c_ in ('full', 'half'):
            a = {k: np.concatenate([x[c_][k] for x in chunks]) for k in ('n', 'p', 'pd', 'ppd')}; n = a['n'][:, None]
            cc = a['ppd'] / n[:, :, None] - np.einsum('rl,rk->rlk', a['p'] / n, a['pd'] / n)
            out.append(np.einsum('rlk->rk', cc - cc.transpose(0, 2, 1)) / 2 / (9 * t))
        return (2 * out[1] - out[0]) / d
    pts = {}
    for t in sorted(set(N.TAUS_FIT_N) | set(N.SEEDS_N)):
        rm = ratio([M['fixed'][t][c] for c in sorted(M['fixed'][t])], t) if t in M['fixed'] else None
        rn = ratio([Nn['fixed'][t][c] for c in sorted(Nn['fixed'][t])], t) if t in Nn['fixed'] else None
        rp = np.concatenate([r for r in (rm, rn) if r is not None]); pts[t] = (rp.mean(axis=0), se(rp), len(rp))
        assert len(rp) == (len(rm) if rm is not None else 0) + (len(rn) if rn is not None else 0)
        row = f"| {t:g} | "
        assert row in md
        line = [l_ for l_ in md.split('\n') if l_.startswith(row) and l_.count('|') == 10][-1]
        assert f" {r2(rp.mean(axis=0), se(rp))} | {z1((rp.mean(axis=0) - 1) / se(rp))} | {', '.join(f'{v:.2f}' for v in rp.std(axis=0, ddof=1))} |" in line
        if rm is not None and rn is not None:
            dz = (rn.mean(axis=0) - rm.mean(axis=0)) / np.sqrt(se(rm) ** 2 + se(rn) ** 2)
            assert f"| {r2(rm.mean(axis=0), se(rm))} | {r2(rn.mean(axis=0), se(rn))} | {z1(dz)} |" in line
    x = np.array(N.TAUS_FIT_N); Y = np.array([pts[t][0] for t in x]); Ee = np.array([pts[t][1] for t in x])
    c, ec, b, eb, cv, chi = [], [], [], [], [], []
    for k in range(3):
        p, cov = np.polyfit(x, Y[:, k], 1, w=1 / Ee[:, k], cov='unscaled')
        b.append(p[0]); c.append(p[1]); eb.append(np.sqrt(cov[0, 0])); ec.append(np.sqrt(cov[1, 1])); cv.append(cov[0, 1])
        chi.append(((Y[:, k] - np.polyval(p, x)) / Ee[:, k]) ** 2 @ np.ones(len(x)))
    c, ec, b, eb, cv, chi = map(np.array, (c, ec, b, eb, cv, chi))
    assert (f"| the line c + bτ | {r2(c, ec)} | {z1((c - 1) / ec)} | {r2(b * 1e-4, eb * 1e-4)} | {z1(b / eb)} | "
            f"{', '.join(f'{v:.2f}' for v in chi)} (1; p = {', '.join(f'{v:.2f}' for v in chi2.sf(chi, 1))}) |") in md
    w = 1 / Ee ** 2; m = (w * Y).sum(axis=0) / w.sum(axis=0); em = 1 / np.sqrt(w.sum(axis=0)); cm = (w * (Y - m) ** 2).sum(axis=0)
    assert f"| a constant | {r2(m, em)} | {z1((m - 1) / em)} | — | — | {', '.join(f'{v:.2f}' for v in cm)} (2; p = " in md
    tx = 1e-3; pr = c + b * tx; pe = np.sqrt(ec ** 2 + tx ** 2 * eb ** 2 + 2 * tx * cv)
    assert (f"At τ = {tx:g} ({pts[tx][2]} seeds, outside the fit): the ratio {r2(pts[tx][0], pts[tx][1])}, the line's prediction {r2(pr, pe)}, "
            f"z = {z1((pts[tx][0] - pr) / np.sqrt(pts[tx][1] ** 2 + pe ** 2))}.") in md


def test_part_N_of_the_md_regenerates_from_the_raw_and_leaves_A_to_M_unchanged():
    N = _N(); raw = _raw(); ks = ('E_rec', 'F', 'G', 'H', 'I', 'J', 'K', 'M')
    md = N.report(*[raw.get(k) for k in 'ABCDE'], **{k: raw.get(k) for k in ks + ('N',)})
    with open(os.path.join(N.OUT, 'gauge_noise.md'), encoding='utf-8') as f:
        assert f.read() == md
    old = N.report(*[raw.get(k) for k in 'ABCDE'], **{k: raw.get(k) for k in ks})
    i = md.index("## N. The finite-τ question"); tail = "Numbers are produced by `experiments/exp_gauge_noise.py`; nothing is edited by hand.\n"
    assert md.index("## M. The missing O(γ²τ) drift") < i and md[:i] + tail == old and md.endswith(tail)
