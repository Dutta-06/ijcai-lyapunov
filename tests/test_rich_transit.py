"""Proposition D‴.48: the rich transit with moving centres — the sliding inverse scale space, the frozen inactive kernels, the
Wasserstein–Fisher–Rao metric's condensation."""
import numpy as np
from scipy.integrate import solve_ivp

from src.core import mirror_free as MF
from src.core import rich_transit as RT

L = MF.FreeLaw(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
K = RT.Bumps(np.linspace(-1, 1, 61), 0.12)
CSTAR, ASTAR = np.array([-0.52, 0.07, 0.61]), np.array([0.9, -0.6, 0.75])
FSTAR = K.phi(CSTAR) @ ASTAR
C0 = np.linspace(-1, 1, 21)
D = C0.size


def flow(eps, tau_end, metric, n_t=221):
    Lam = RT.transit_scale(L, eps)
    taus = np.linspace(0, tau_end, n_t)
    sol = solve_ivp(RT.moving_rhs(L, K, FSTAR, metric), (0, Lam * tau_end), np.r_[np.zeros(D), np.log(np.full(D, eps)), C0],
                    method='LSODA', rtol=1e-9, atol=1e-12, t_eval=Lam * taus)
    al, lb, cc = sol.y[:D].T, sol.y[D:2 * D].T, sol.y[2 * D:].T
    return taus, np.exp(lb) * L.amp(al)[0] / L.m, cc, al, Lam


def test_critical_curve_and_the_depth_three_unit():
    """On Q = −c the coefficient is FreeLaw's, and w ≈ 2cα³/(ma_max²) at small α."""
    al = np.array([-1.3, -0.2, 0.01, 0.4, 2.0])
    assert np.allclose(RT.critical_w(L, al), [L.w(a, -L.c) for a in al], rtol=1e-12, atol=0)
    a = 1e-3
    assert abs(RT.critical_w(L, a) / (2 * L.c * a ** 3 / (L.m * L.A ** 2)) - 1) < 1e-5


def test_sliding_iss_recovers_the_bumps():
    """The sliding inverse scale space on three separated bumps: three entries (the kernels started nearest to the bumps, the largest
    amplitude first), each sliding onto its bump; the fit exact."""
    R = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D))
    assert [e[1] for e in R['events']] == [5, 16, 11] and all(e[2] == 'enter' for e in R['events'])
    order = np.argsort(R['c'][R['act']])
    assert np.abs(R['c'][R['act']][order] - CSTAR).max() < 1e-6
    assert np.abs(R['w'][R['act']][order] - ASTAR).max() < 1e-6
    assert R['hist'][-1][3] < 1e-9


def test_paradigm_flow_follows_the_sliding_path():
    """Proposition D‴.48(iii): at ε = 10⁻⁵ the loop forms the same three kernels, entering at the sliding path's times, at the bumps;
    the others stay where they started (their mass negligible)."""
    R = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D))
    t_ev = {e[1]: e[0] for e in R['events']}
    taus, W, cc, al, Lam = flow(1e-5, 1.3 * R['events'][-1][0], 'paradigm', n_t=441)
    act = np.abs(W[-1]) > 1e-2 * np.abs(ASTAR).max()
    assert sorted(np.flatnonzero(act)) == [5, 11, 16]
    for i in (5, 11, 16):
        k = np.flatnonzero(np.abs(W[:, i]) > 1e-2 * np.abs(ASTAR).max())[0]
        assert abs(taus[k] - t_ev[i]) < 0.15
    order = np.argsort(cc[-1][act])
    assert np.abs(cc[-1][act][order] - CSTAR).max() < 1e-3
    assert np.abs(cc[-1][~act] - C0[~act]).max() < 0.05 and np.abs(W[-1][~act]).max() < 1e-3


def test_inactive_kernels_obey_the_frozen_bound():
    """Proposition D‴.48(iii)(a): before the first entry each kernel's coefficient is the core's profile
    a_max(εb)^{3/2}cos(πμ/2)/(m(2c)^{1/2}sin³(πμ/2)) (μ its dual coordinate's margin; within 3% where μ ≥ 10ε^{1/2} — round 37: the μ^{−3}
    bound is 4 times looser), so its centre hardly moves."""
    eps = 1e-5
    taus, W, cc, al, Lam = flow(eps, 9.0, 'paradigm', n_t=91)
    q = L.Q(np.zeros(D), np.full(D, eps))
    lam1 = np.pi * L.A * np.sqrt((L.c - q) / (q + L.c)) / (4 * L.G * L.c)
    n_checked = 0
    for k in range(1, len(taus)):
        mu = 1 - np.abs(L.eta(al[k], q)) / lam1
        assert mu.min() > 0
        prof = L.core_tail(eps, mu)
        assert np.all(np.abs(W[k]) <= L.A * eps ** 1.5 / (L.m * np.sqrt(2 * L.c)) / mu ** 3)
        ok = mu >= 10 * np.sqrt(eps)
        assert np.abs(np.abs(W[k][ok]) / prof[ok] - 1).max() < 0.03
        n_checked += ok.sum()
    assert n_checked > 1000
    assert np.abs(cc[-1] - C0).max() < 1e-2


def test_wfr_metric_condenses_onto_the_limit():
    """Proposition D‴.48(iv): with the Wasserstein–Fisher–Rao metric the inactive kernels' centres follow the limit
    dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ⟩, du/dτ = ⟨f*φ⟩ before the first entry, gathering at three positions."""
    tau = 10.0
    lim = solve_ivp(RT.wfr_inactive_rhs(L, K, FSTAR, np.ones(D)), (0, tau), np.r_[np.zeros(D), C0], method='LSODA', rtol=1e-10,
                    atol=1e-12)
    c_lim = np.sort(lim.y[D:, -1])
    assert np.sum(np.diff(c_lim) > 1e-3) + 1 == 3
    taus, W, cc, al, Lam = flow(1e-6, tau, 'wfr', n_t=11)
    assert np.abs(np.sort(cc[-1]) - c_lim).max() < 1e-5
    assert np.abs(W[-1]).max() < 1e-3


def test_orbit_filter_is_the_periodic_error():
    """Proposition D‴.49: OrbitBumps.filt is the periodic steady state of ė = −k_ee + r on the reference's orbit — for
    r = sin(ωt), e = (k_e sin ωt − ω cos ωt)/(k_e² + ω²)."""
    ke, om = 10.0, 0.2
    K = RT.OrbitBumps(lambda t: np.sin(om * t), 2 * np.pi / om, 256, 0.12, ke)
    ts = (2 * np.pi / om) * np.arange(256) / 256
    e = K.filt(np.sin(om * ts))
    assert np.abs(e - (ke * np.sin(om * ts) - om * np.cos(om * ts)) / (ke ** 2 + om ** 2)).max() < 1e-12
    assert np.allclose(K.xs, np.sin(om * ts))


def test_moving_widths_identify_wider_bumps():
    """Proposition D‴.48 with moving widths: three bumps of width 0.2 (wider than the kernels' 0.12) are formed at ε = 10⁻⁶ by exactly
    three kernels, at the bumps' centres, widths and amplitudes."""
    fs = np.exp(-(K.xs[:, None] - CSTAR[None, :]) ** 2 / (2 * 0.2 ** 2)) @ ASTAR
    eps, tau = 1e-6, 30.0
    Lam = RT.transit_scale(L, eps)
    sol = solve_ivp(RT.moving_rhs_widths(L, K, fs), (0, Lam * tau), np.r_[np.zeros(D), np.log(np.full(D, eps)), C0, np.full(D, 0.12)],
                    method='LSODA', rtol=1e-9, atol=1e-12)
    z = sol.y[:, -1]
    W = np.exp(z[D:2 * D]) * L.amp(z[:D])[0] / L.m
    act = np.abs(W) > 1e-2 * np.abs(ASTAR).max()
    assert act.sum() == 3
    o = np.argsort(z[2 * D:3 * D][act])
    assert np.abs(z[2 * D:3 * D][act][o] - CSTAR).max() < 1e-3
    assert np.abs(z[3 * D:][act][o] - 0.2).max() < 1e-3
    assert np.abs(W[act][o] - ASTAR).max() < 1e-3


def test_position_laws_and_mobility():
    """The drifts at a random state (round 37: the metric's mass factor and a′ survived): α̇ = Gβa′g, ċ = Gβa⟨r∂_cφ⟩ (the paradigm's metric)
    or Ga⟨r∂_cφ⟩ (WFR); on the curves (critical_rhs with b0) the same with β = b0 + 2c sinh²(α/a_max)."""
    rng = np.random.default_rng(4)
    al, c = rng.uniform(-2, 2, D), C0 + rng.uniform(-0.03, 0.03, D)
    b0 = rng.uniform(0.5, 2.0, D) * 1e-3
    b = b0 + 2 * L.c * np.sinh(al / L.A) ** 2
    a = L.A * np.tanh(al / L.A); da = 1 - np.tanh(al / L.A) ** 2
    r = FSTAR - K.phi(c) @ (b * a / L.m)
    g = K.phi(c).T @ r / K.xs.size
    h = K.dphi(c).T @ r / K.xs.size
    for metric, speed in (('paradigm', L.G * b * a), ('wfr', L.G * a)):
        z = RT.moving_rhs(L, K, FSTAR, metric)(0.0, np.r_[al, np.log(b), c])
        assert np.allclose(z[:D], L.G * b * da * g, rtol=1e-12, atol=1e-15)
        assert np.allclose(z[D:2 * D], L.gamma * a * g, rtol=1e-12, atol=1e-15)
        assert np.allclose(z[2 * D:], speed * h, rtol=1e-12, atol=1e-15)
        zc = RT.critical_rhs(L, K, FSTAR, metric, b0=b0)(0.0, np.r_[al, c])
        assert np.allclose(zc[:D], L.G * b * da * g, rtol=1e-11, atol=1e-15)
        assert np.allclose(zc[D:], speed * h, rtol=1e-12, atol=1e-15)
    assert np.allclose(RT.critical_w(L, al, b0), b * a / L.m, rtol=1e-13)


def test_creep_speed_is_the_core_profile():
    """Proposition D‴.48(iii)(a), round 37: an inactive kernel with the mass εb whose dual coordinate sits at (1 − μ)λ₁ creeps at
    dc/dτ = (πa_max²/(4c))εb^{3/2}cos(πμ/2)/sin³(πμ/2)·s⟨r∂_cφ⟩ — not at the bound's μ^{−3} form (larger by π³/8 as μ → 0)."""
    eps, bj, mu = 1e-8, 2.5, np.array([0.3, 0.1, 0.03])
    Lam = RT.transit_scale(L, eps)
    q = L.Q(0.0, eps * bj)
    lam1 = np.pi * L.A * np.sqrt((L.c - q) / (q + L.c)) / (4 * L.G * L.c)
    al = L.alpha_of_eta((1 - mu) * lam1, np.full(3, q))
    c = np.array([-0.4, 0.05, 0.5])
    f = RT.critical_rhs(L, K, FSTAR, 'paradigm', b0=np.full(3, eps * bj))
    dc = Lam * f(0.0, np.r_[al, c])[3:]
    r = FSTAR - K.phi(c) @ RT.critical_w(L, al, eps * bj)
    h = K.dphi(c).T @ r / K.xs.size
    prof = np.pi * L.A ** 2 / (4 * L.c) * eps * bj ** 1.5 * np.cos(np.pi * mu / 2) / np.sin(np.pi * mu / 2) ** 3 * h
    assert np.abs(dc / prof - 1).max() < 2e-3
    old = np.pi * L.A ** 2 / (4 * L.c) * eps * mu ** -3 * np.abs(h)
    assert np.abs(np.abs(dc[-1]) / old[-1] * np.pi ** 3 / 8 / bj ** 1.5 - 1) < 0.01


def test_wfr_limit_along_the_trajectory_with_unequal_masses():
    """Proposition D‴.48(iv), round 37: with the masses εb_i (ρ_i = b_i^{−1/2} ≠ 1) the loop's centres under the WFR metric follow the
    limit dc/dτ = (πa_max²/(4cρ))tan(πu/(2ρ))⟨f*∂_cφ⟩ along the way (τ = 3, 5, 7), where the limit with the masses ignored (ρ = 1) does not."""
    b = np.exp(np.random.default_rng(31).uniform(np.log(0.25), np.log(4.0), D))
    eps = 1e-6
    Lam = RT.transit_scale(L, eps)
    taus = np.array([3.0, 5.0, 7.0])
    sol = solve_ivp(RT.critical_rhs(L, K, FSTAR, 'wfr', b0=eps * b), (0, Lam * taus[-1]), np.r_[np.zeros(D), C0], method='LSODA',
                    rtol=1e-10, atol=1e-13, t_eval=Lam * taus)
    lim = solve_ivp(RT.wfr_inactive_rhs(L, K, FSTAR, b ** -0.5), (0, taus[-1]), np.r_[np.zeros(D), C0], method='LSODA', rtol=1e-11,
                    atol=1e-13, t_eval=taus)
    lim1 = solve_ivp(RT.wfr_inactive_rhs(L, K, FSTAR, np.ones(D)), (0, taus[-1]), np.r_[np.zeros(D), C0], method='LSODA', rtol=1e-11,
                     atol=1e-13, t_eval=taus)
    assert np.abs(RT.critical_w(L, sol.y[:D], (eps * b)[:, None])).max() < 1e-3          # before any entry
    assert np.abs(sol.y[D:] - lim.y[D:]).max() < 1e-5
    assert np.abs(sol.y[D:] - lim1.y[D:]).max() > 1e-3


def test_sliding_iss_ties_and_seeds():
    """Round 37: a bump midway between two kernels is a tie — both enter together and share it (0.4 each at its centre), as the loop does
    at every ε; the sliding path does not depend on the entrants' seed (10⁻³ … 10⁻¹¹) or on the fast flow's first horizon."""
    R = RT.sliding_iss(L, K, K.phi(np.array([0.05])) @ np.array([0.8]), C0, np.ones(D))
    assert [(e[1], e[2]) for e in R['events']] == [(10, 'enter'), (11, 'enter')] and R['events'][0][0] == R['events'][1][0]
    assert np.allclose(R['c'][R['act']], 0.05, atol=1e-6) and np.allclose(R['w'][R['act']], 0.4, atol=1e-6)
    ref = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D))
    for kw in (dict(seed=1e-3), dict(seed=1e-11), dict(t_fast=1e3)):
        R = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D), **kw)
        assert [(e[1], e[2]) for e in R['events']] == [(e[1], e[2]) for e in ref['events']] and all(R['relaxed'])
        assert max(abs(e1[0] - e2[0]) for e1, e2 in zip(R['events'], ref['events'])) < 1e-6
        assert np.abs(np.sort(R['c'][R['act']]) - CSTAR).max() < 1e-6


def test_polish_reproduces_the_sliding_path_and_resolves_close_bumps():
    """`sliding_iss(polish=True)` (the block's sign-constrained least-squares fit after the first t_fast) gives the same path as the fast
    flow on random three-bump targets, and relaxes two bumps 0.05 apart — where the flow crawls along the block's Δ⁶ mode — to the exact
    two-kernel fit at once."""
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
    import exp_mirror_transit as X
    import exp_transit_mass as T
    L, K, c0 = T.setup()
    for cs, As in T.random_targets(21, 3):
        f = X.bumps_target(K, cs, As)
        R0 = RT.sliding_iss(L, K, f, c0, np.ones(21))
        R1 = RT.sliding_iss(L, K, f, c0, np.ones(21), polish=True)
        assert [(j, k) for t, j, k, c in R0['events']] == [(j, k) for t, j, k, c in R1['events']]
        assert max(abs(a[0] - b[0]) for a, b in zip(R0['events'], R1['events'])) < 1e-9
        assert np.abs(R0['w'] - R1['w']).max() < 1e-9 and np.abs(R0['c'] - R1['c'])[R0['act']].max() < 1e-9
    cs = (0.032 - 0.05 * 0.6 / 1.4, 0.032 + 0.05 * 0.8 / 1.4)
    f = X.bumps_target(K, cs, (0.8, 0.6))
    R = RT.sliding_iss(L, K, f, c0, np.ones(21), polish=True)
    a = np.flatnonzero(R['act'])
    assert a.size == 2 and all(R['relaxed']) and R['hist'][-1][3] < 1e-12
    assert np.abs(np.sort(R['c'][a]) - np.array(cs)).max() < 1e-9 and np.abs(R['w'][a] - np.array([0.8, 0.6])).max() < 1e-9


def test_ascent_target_climbs_to_the_residuals_nearest_peak():
    """`ascent_target`: a near-seed kernel's centre climbs s⟨rφ(c)⟩ — from 0.3 in the residual φ(·; 0) to 0 (within the step), and for the
    sign −1 away from it; with two bumps of r it stops at the nearer peak's basin."""
    K = RT.Bumps(np.linspace(-1, 1, 61), 0.12)
    r = K.phi(np.array([0.0]))[:, 0]
    assert abs(RT.ascent_target(K, r, 0.3, 1.0)[0]) <= 0.12 / 40 + 1e-12
    assert RT.ascent_target(K, r, 0.3, -1.0)[0] > 0.3
    r2 = K.phi(np.array([-0.4]))[:, 0] + K.phi(np.array([0.5]))[:, 0]
    assert abs(RT.ascent_target(K, r2, 0.3, 1.0)[0] - 0.5) < 0.02 and abs(RT.ascent_target(K, r2, -0.2, 1.0)[0] + 0.4) < 0.02


def test_relax_chunk_is_lsoda_below_its_budget_and_radau_above():
    """`relax_chunk`: below the budget the LSODA call itself (bit for bit), above it the chunk redone from its start by Radau."""
    rhs = RT.critical_rhs(L, K, FSTAR)
    z0 = np.r_[np.array([0.3, -0.2, 1e-3]), np.array([-0.5, 0.1, 0.6])]
    ref = solve_ivp(rhs, (0, 50.0), z0, method='LSODA', rtol=1e-10, atol=1e-13).y[:, -1]
    z, redone = RT.relax_chunk(rhs, z0, 50.0)
    assert not redone and np.array_equal(z, ref)
    rad = solve_ivp(rhs, (0, 50.0), z0, method='Radau', rtol=1e-10, atol=1e-13).y[:, -1]
    z, redone = RT.relax_chunk(rhs, z0, 50.0, budget=10)
    assert redone and np.array_equal(z, rad)


def test_sliding_paths_through_a_block_that_stalls_lsoda():
    """Three features at −0.587, 0.013, 0.613 with amplitudes 10⁻³, 10⁻⁶, 1 (referee round 46's wide range): the third block's coefficients
    span six decades and LSODA integrates it non-stiffly (6.7·10⁶ evaluations, 610 s). With the budget at 2·10⁴ the chunk is redone by
    Radau: both paths enter the three features' kernels — the largest first — and relax to the exact fit."""
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
    import exp_mirror_transit as X
    from src.core import rich_transit_nd as RN
    c0 = np.linspace(-1, 1, 21)
    f = X.bumps_target(K, np.array([-0.6, 0.0, 0.6]) + 0.013, np.array([1e-3, 1e-6, 1.0]))
    budget = RT.STIFF_BUDGET
    RT.STIFF_BUDGET = 20000
    try:
        for R in (RT.sliding_iss(L, K, f, c0, np.ones(21), polish=True, exit_tol=1e-12, max_events=12),
                  RN.sliding_iss_nd(L, RN.BumpsND(K.xs[:, None], K.sig), f, c0[:, None], np.ones(21), polish=True, exit_tol=1e-12,
                                    max_events=12)):
            assert [(j, k) for t, j, k, c in R['events']] == [(16, 'enter'), (4, 'enter'), (10, 'enter')]
            assert all(R['relaxed']) and R['stiff_redone'] >= 1 and R['hist'][-1][3] < 1e-14
    finally:
        RT.STIFF_BUDGET = budget


def test_relax_chunk_redoes_a_failed_lsoda_solve():
    """From a relaxed block LSODA fails on a chunk of 3·10¹¹ ("unexpected istate", its last state the start) — `relax_chunk` used to
    return that start as if the chunk had run (referee round 48); it redoes the chunk by Radau and says so (how = 2)."""
    import os
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
    import exp_noisy_certificate as NC
    L_, K_, c0, f, cert, R0, P = NC.noiseless(np.array([-0.52, 0.02, 0.61]), np.array([0.8, -0.35, 0.9]))
    w, c = R0['hist'][-1][2], R0['hist'][-1][1]
    i = np.flatnonzero(w)
    z0 = np.r_[np.sign(w[i]) * RT._alpha_of(L_, np.abs(w[i]), RT.coef_crit), c[i]]
    rhs = RT.critical_rhs(L_, K_, f)
    sol = solve_ivp(rhs, (0, 3e11), z0, method='LSODA', rtol=1e-10, atol=1e-13)
    assert not sol.success and sol.t[-1] == 0.0
    z, how = RT.relax_chunk(rhs, z0, 3e11)
    rad = solve_ivp(rhs, (0, 3e11), z0, method='Radau', rtol=1e-10, atol=1e-13)
    assert how == 2 and rad.success and np.array_equal(z, rad.y[:, -1])
    z, how = RT.relax_chunk(rhs, z0, 1e10)
    assert how == 0


def test_leaves_with_the_kernels_own_contribution_removed():
    """A small coefficient with an inward correlation leaves; under exit_kkt only if the correlation stays inward once the kernel's own
    contribution w_i⟨φ_iφ_i⟩ is removed — else the removal turns it outward and it re-enters at the same τ (referee round 48)."""
    si = np.array([1.0, 1.0, -1.0, 1.0])
    g = np.array([-1e-9, -1e-9, 1e-9, 1e-9])
    wi = np.array([1e-7, 1e-12, -1e-12, 1e-12])
    big = np.array([False, False, False, False])
    own = np.full(4, 0.5)
    assert RT.leaves(si, g, wi, big).tolist() == [True, True, True, False]
    assert RT.leaves(si, g, wi, big, own).tolist() == [False, True, True, False]
    assert not RT.leaves(si, g, wi, np.ones(4, bool), own).any()


def test_the_nd_path_stops_at_its_rhs_budget():
    """`sliding_iss_nd(max_rhs=…)`: a budget of right-hand-side evaluations for the whole path, deterministic; spent, the path stops before
    the block that spent it is relaxed (out_of_budget, that block unrelaxed); a budget the path never reaches leaves it unchanged."""
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
    import exp_mirror_transit as X
    from src.core import rich_transit_nd as RN
    c0 = np.linspace(-1, 1, 21)
    f = X.bumps_target(K, np.array([-0.5, 0.6]) + 0.013, np.array([0.8, -0.5]))
    KN = RN.BumpsND(K.xs[:, None], K.sig)
    R = RN.sliding_iss_nd(L, KN, f, c0[:, None], np.ones(21), polish=True, max_events=6)
    assert not R['out_of_budget'] and R['nfev'] == 0 and all(R['relaxed'])
    Rb = RN.sliding_iss_nd(L, KN, f, c0[:, None], np.ones(21), polish=True, max_events=6, max_rhs=10 ** 7)
    assert not Rb['out_of_budget'] and 0 < Rb['nfev'] <= 10 ** 7
    assert [e[:3] for e in Rb['events']] == [e[:3] for e in R['events']] and np.array_equal(Rb['w'], R['w'])
    Rs = RN.sliding_iss_nd(L, KN, f, c0[:, None], np.ones(21), polish=True, max_events=6, max_rhs=50)
    assert Rs['out_of_budget'] and Rs['relaxed'] == [False] and Rs['nfev'] == 51


def test_a_failed_chunk_leaves_its_block_unrelaxed_without_exits(monkeypatch):
    """Referee round 54: a chunk that fails (`relax_chunk`'s how = 3) at the first event's first chunk or at a later one ends that
    block's relaxation unrelaxed, counts once, and no kernel exits on its correlations — with and without exit_kkt, in one dimension and
    in ℝⁿ (before, the exit rule read undefined or the previous block's arrays); the path goes on from the unrelaxed block."""
    from src.core import rich_transit_nd as RN
    real = RT.relax_chunk
    calls = [0]

    def failing(at):
        def rc(rhs, z, chunk, **kw):
            calls[0] += 1
            return (z, 3) if calls[0] == at else real(rhs, z, chunk, **kw)
        return rc
    KN = RN.BumpsND(K.xs[:, None], K.sig)
    paths = (lambda kkt: RT.sliding_iss(L, K, FSTAR, C0, np.ones(D), exit_kkt=kkt, max_events=8),
             lambda kkt: RN.sliding_iss_nd(L, KN, FSTAR, C0[:, None], np.ones(D), exit_kkt=kkt, max_events=8))
    for path in paths:
        for at in (1, 3):
            for kkt in (False, True):
                calls[0] = 0
                monkeypatch.setattr(RT, 'relax_chunk', failing(at))
                R = path(kkt)
                assert R['solver_failed'] == 1 and False in R['relaxed'] and calls[0] > at
                assert (R['stiff_redone'], R['solver_redone']) == (0, 0)             # a failure is not a redo (round 60)
                tf = R['hist'][R['relaxed'].index(False) + 1][0]                     # the failed block's τ (the path goes on)
                assert not any(e[2] == 'exit' and e[0] == tf for e in R['events'])
    monkeypatch.setattr(RT, 'relax_chunk', real)


def test_relax_chunk_halves_what_radau_cannot_solve(monkeypatch):
    """Referee round 54: `relax_chunk` against a fake solver — LSODA always failing, Radau solving ż = −z exactly on the spans and starts
    it accepts and otherwise failing with its start as its last state. Accepted whole, the chunk is redone (how = 2); accepted only in
    quarters, it is halved twice, each half from the state the first reached, and ends at the exact state; never accepted, it is halved
    along its first half alone to the depth RELAX_SPLITS and fails at its start (how = 3); a first half failing partway — its first quarter
    solved, its second never — fails the chunk at the state the first quarter reached, though the second half would be accepted from
    there (the time between not integrated)."""
    from types import SimpleNamespace
    calls = []

    def fake(ok):
        def solve_ivp(fun, span, y0, method, **kw):
            T = span[1] - span[0]
            good = method == 'Radau' and ok(T, y0[0])
            calls.append((method, T))
            return SimpleNamespace(success=good, y=np.c_[y0, y0 * np.exp(-T) if good else y0])
        return solve_ivp
    rhs = lambda s, z: -z
    z0, C = np.array([1.0]), 8.0
    monkeypatch.setattr(RT, 'solve_ivp', fake(lambda T, y: True))
    assert RT.relax_chunk(rhs, z0, C) == (np.exp(-C), 2)
    monkeypatch.setattr(RT, 'solve_ivp', fake(lambda T, y: T <= C / 4))
    z, how = RT.relax_chunk(rhs, z0, C)
    assert how == 2 and abs(z[0] / np.exp(-C) - 1) < 1e-14
    calls.clear()
    monkeypatch.setattr(RT, 'solve_ivp', fake(lambda T, y: False))
    z, how = RT.relax_chunk(rhs, z0, C)
    assert how == 3 and z[0] == 1.0
    assert [T for m, T in calls if m == 'Radau'] == [C / 2 ** d for d in range(RT.RELAX_SPLITS + 1)] and RT.RELAX_SPLITS == 6
    za = np.exp(-C / 4)
    at_za = lambda y: abs(y - za) < 1e-12
    monkeypatch.setattr(RT, 'solve_ivp', fake(lambda T, y: (T <= C / 4 and not at_za(y)) or (T == C / 2 and at_za(y))))
    z, how = RT.relax_chunk(rhs, z0, C)
    assert how == 3 and z[0] == za


def test_the_path_counts_its_redone_and_failed_chunks(monkeypatch):
    """Referee round 54: the path counts the chunks `relax_chunk` redid after a stall (how = 1) as 'stiff_redone' and after a failed solve
    (how = 2) as 'solver_redone', separately; the flags alone change nothing else on the path."""
    real = RT.relax_chunk
    calls = [0]
    ref = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D))

    def rc(rhs, z, chunk, **kw):
        calls[0] += 1
        z, how = real(rhs, z, chunk, **kw)
        return z, {1: 1, 2: 1, 3: 2}.get(calls[0], how)
    monkeypatch.setattr(RT, 'relax_chunk', rc)
    R = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D))
    assert calls[0] >= 3 and (R['stiff_redone'], R['solver_redone'], R['solver_failed']) == (2, 1, 0)
    assert (ref['stiff_redone'], ref['solver_redone'], ref['solver_failed']) == (0, 0, 0)
    assert R['events'] == ref['events'] and all(R['relaxed']) and np.array_equal(R['w'], ref['w'])


def test_the_path_applies_the_exit_rule_to_its_relaxed_block(monkeypatch):
    """Referee round 54: the exits the path takes from a given relaxed block. On 0.8φ(· + 0.5) − 0.6φ(· − 0.5) the kernel at −0.5 enters
    (sign +); a fake relaxation moves it near the target's correlation's zero with a coefficient 5·10⁻⁷, below the exit tolerance. Where
    the target's correlation there is outward (2·10⁻⁸) and the residual's inward only through the kernel's own contribution w⟨φφ⟩, the
    plain rule lets it leave — and it would re-enter at once — while exit_kkt keeps it; where the target's correlation is inward
    (−10⁻⁵), both let it leave."""
    from scipy.optimize import brentq
    f = K.phi(np.array([-0.5, 0.5])) @ np.array([0.8, -0.6])
    nx = K.xs.size
    g0 = lambda c: float(K.phi(np.array([c]))[:, 0] @ f) / nx
    al = RT._alpha_of(L, np.array([5e-7]), RT.coef_crit)[0]
    for target, kept in ((2e-8, True), (-1e-5, False)):
        cx = brentq(lambda c: g0(c) - target, -0.5, 0.5, xtol=1e-14)
        monkeypatch.setattr(RT, 'relax_chunk', lambda rhs, z, chunk, **kw: (np.array([al, cx]), 0))
        for kkt in (False, True):
            R = RT.sliding_iss(L, K, f, C0, np.ones(D), exit_kkt=kkt, max_events=1, w_scale=1.0)
            assert R['events'][0][1:3] == (5, 'enter') and R['relaxed'] == [True]
            ex = [e for e in R['events'] if e[2] == 'exit']
            assert ex == ([] if kkt and kept else [(R['events'][0][0], 5, 'exit', cx)])


def _wide_path(r_, **kw):
    """D‴.58's wide feature 0.8φ_σ*(· − 0.03), σ* = σ(1 + (τ/σ)²)^{1/2}, on its 321 samples, from D‴.48's grid (polished, KKT exits)."""
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
    import exp_transit_mass as T
    L_, _, c0 = T.setup()
    xs = np.linspace(-1.6, 1.6, 321)
    Kw = RT.Bumps(xs, 0.12)
    f = 0.8 * np.exp(-(xs - 0.03) ** 2 / (2 * 0.0144 * (1 + r_ ** 2)))
    return RT.sliding_iss(L_, Kw, f, c0, np.ones(c0.size), polish=True, max_events=20, tol=0.0, max_doublings=24, exit_kkt=True, **kw)


def test_the_polish_budget_gradient_test_and_repolish(monkeypatch):
    """Referees' rounds 54 and 60, on D‴.58's wide feature (τ/σ = ½), with a loose stationarity test so that blocks pass it unpolished:
    the path hands its polish_nfev and polish_gtol (None: no gradient test) to every polish, the first and the repeated, and no solve
    spends more; a repeated polish is accepted only if it lowers the fit's rms by more than a relative 10⁻⁹, at most `repolish` times a
    block, and after one the flow resumes for as long as the block has run — every chunk the block's elapsed time; repolish = 0 repeats
    none."""
    real_p, real_c = RT.polish_block, RT.relax_chunk
    kws, chunks = [], []

    def pb(*a, **kw):
        kws.append((kw['max_nfev'], kw['gtol']))
        return real_p(*a, **kw)

    def rc(rhs, z, chunk, **kw):
        chunks.append((z.size, chunk))
        return real_c(rhs, z, chunk, **kw)
    monkeypatch.setattr(RT, 'polish_block', pb)
    monkeypatch.setattr(RT, 'relax_chunk', rc)
    R = _wide_path(0.5, max_entries=3, polish_step=2.0, polish_nfev=3, repolish=3, gtol=1e-2, polish_gtol=None)
    assert kws and set(kws) == {(3, None)} and all(p['nfev'] <= 3 for pl in R['polish_log'] for p in pl)
    assert [e[2] for e in R['events']] == ['enter'] * 3                              # no exits: the blocks' sizes tell them apart
    rep = [[p for p in pl if 'e0' in p] for pl in R['polish_log']]
    assert all(len(x) <= 3 for x in rep) and sum(p['accepted'] for x in rep for p in x) >= 2
    assert all(p['e1'] < p['e0'] * (1 - 1e-9) for x in rep for p in x if p['accepted'])
    blocks = [[c for n, c in chunks if n == m] for m in sorted({n for n, c in chunks})]
    assert max(len(b) for b in blocks) >= 4
    for b in blocks:                                                                  # t_fast, then the elapsed time each chunk
        assert b[0] == 1e5 and all(b[k] == sum(b[:k]) for k in range(1, len(b)))
    kws.clear()
    R0 = _wide_path(0.5, max_entries=3, polish_step=2.0, polish_nfev=3, repolish=0, gtol=1e-2)
    assert set(kws) == {(3, 1e-15)} and not any('e0' in p for pl in R0['polish_log'] for p in pl)


def test_a_non_improving_repolish_changes_nothing(monkeypatch):
    """Referee round 60: with a polish that returns its start, every block that passes the stationarity test is polished once more,
    the repeat is rejected (its rms exactly its start's) and the block ends there — the path equals the one without repeats."""
    monkeypatch.setattr(RT, 'polish_block', lambda K_, f_, w, c, s, **kw: (w.copy(), c.copy(), dict(status=1, nfev=1)))
    R3 = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D), polish=True, repolish=3)
    R0 = RT.sliding_iss(L, K, FSTAR, C0, np.ones(D), polish=True, repolish=0)
    for pl, ok in zip(R3['polish_log'], R3['relaxed']):
        rep = [p for p in pl if 'e0' in p]
        assert len(rep) == ok and all(not p['accepted'] and p['e1'] == p['e0'] for p in rep)
    assert R3['events'] == R0['events'] and np.array_equal(R3['w'], R0['w']) and np.array_equal(R3['c'], R0['c'])
    assert not any('e0' in p for pl in R0['polish_log'] for p in pl)


def test_the_polish_log_and_the_rejection_count():
    """Referee round 60: the log's first polishes (status, evaluations, accepted) — their rejections are the path's 'polish_rejected',
    each followed by another try, the block's last first polish accepted where the block was polished."""
    R = _wide_path(1.0, max_entries=4, polish_step=0.1)
    first = [[p for p in pl if 'e0' not in p] for pl in R['polish_log']]
    assert R['polish_rejected'] == sum(not p['accepted'] for x in first for p in x) >= 3
    assert all(x[-1]['accepted'] and all(not p['accepted'] for p in x[:-1]) for x in first if x)
    assert all(p['status'] in (0, 1, 2, 3, 4) and 1 <= p['nfev'] <= 5000 for x in first for p in x)


def test_polish_block_gradient_test():
    """Referee round 60: `polish_block`'s gtol is TRF's absolute gradient test — a loose one (10³) stops the solve at its start (status 1,
    one evaluation); None disables it, and the solve ends on its step or decrease tests (status 2–4) at the exact fit."""
    w0, c0_, s = np.array([0.5, -0.4, 0.6]), np.array([-0.5, 0.1, 0.6]), np.array([1, -1, 1])
    w, c, info = RT.polish_block(K, FSTAR, w0, c0_, s, info=True, gtol=1e3)
    assert info == dict(status=1, nfev=1) and np.array_equal(c, c0_)
    w, c, info = RT.polish_block(K, FSTAR, w0, c0_, s, info=True, gtol=None)
    assert info['status'] in (2, 3, 4) and np.abs(c - CSTAR).max() < 1e-12 and np.abs(w - ASTAR).max() < 1e-12


def test_entry_info_is_the_duals_rebuilt_from_the_path():
    """Referee round 60: `sliding_iss`'s entry_info against the duals rebuilt from the path's states — u_j the sum over the earlier
    iterations i of (τ_{i+1} − τ_i)⟨r_iφ_j⟩ (r_i the residual of hist[i]; every entrant inactive until its entry), its rate ⟨rφ_j⟩, the
    largest rate over the inactive kernels alone, its time to the wall (1 − sign(v_j)u_j)/|v_j|, the least — on D‴.58's wide feature
    (polished), on the three bumps (a negative dual) and on them with a loose stationarity test (blocks whose own correlations exceed
    every inactive rate)."""
    xs = np.linspace(-1.6, 1.6, 321)
    Kw = RT.Bumps(xs, 0.12)
    fw = 0.8 * np.exp(-(xs - 0.03) ** 2 / (2 * 0.0144 * 1.25))
    paths = ((_wide_path(0.5, max_entries=5, polish_step=2.0), Kw, fw, 5),
             (RT.sliding_iss(L, K, FSTAR, C0, np.ones(D)), K, FSTAR, 3),
             (RT.sliding_iss(L, K, FSTAR, C0, np.ones(D), gtol=0.3, t_fast=100.0, max_events=6), K, FSTAR, 6))
    over = 0
    for R, Kp, f, n in paths:
        assert [e[2] for e in R['events']] == ['enter'] * n and len(R['entry_info']) == n
        vs = []
        for tau, c, w, rms in R['hist']:
            v = Kp.phi(c).T @ (f - Kp.phi(c) @ w) / Kp.xs.size
            vs.append((v, w == 0))
            over += bool((w != 0).any()) and np.abs(v[w != 0]).max() > np.abs(v[w == 0]).max()
        taus = [h[0] for h in R['hist']]
        for e in R['entry_info']:
            i, j = e['it'], e['j']
            u = sum((taus[k + 1] - taus[k]) * vs[k][0][j] for k in range(i))
            v, ina = vs[i]
            sc = np.abs(v).max()
            assert ina[j] and abs(e['u0'] - u) <= 1e-9 * max(1.0, abs(u)) and abs(e['rate'] - v[j]) <= 1e-12 * sc
            assert abs(e['rate_max'] - np.abs(v[ina]).max()) <= 1e-12 * sc
            assert abs(e['t_wall'] - (1 - np.sign(v[j]) * e['u0']) / abs(v[j])) <= 1e-9 * e['t_wall'] <= 1e-9 * e['t_min'] * (1 + 1e-9)
        assert sum(abs(e['u0']) > 0.1 for e in R['entry_info']) >= 2                  # later entrants' duals charged before they enter
    assert min(e['u0'] for R, Kp, f, n in paths[1:2] for e in R['entry_info']) < -0.5 and over >= 2
