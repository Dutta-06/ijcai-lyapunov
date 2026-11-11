"""
Checks for PARADIGM.md Theorem C′ (the OT leak toward a prior with a different number of atoms: barycentric projection
of the optimal plan, dissipation up to the within-cell variance, centroidal fixed points, the induced coupling bound and
the Lyapunov inequality), using src/core/transport.py.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core import transport as T

BOX = dict(a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)


def test_plan_reduces_to_assignment_and_pythagorean_split():
    rng = np.random.default_rng(0)
    Ginv = 1.0 / np.array([2.0, 1.0, 0.5])
    for _ in range(10):
        P = rng.normal(size=(5, 3)); Q = rng.normal(size=(5, 3))
        g, W2 = T.optimal_plan(P, Q, Ginv); pi, W2b = T.optimal_assignment(P, Q, Ginv)
        assert abs(W2 - W2b) < 1e-10 * (1 + W2b)
        assert np.allclose(np.sort(g.ravel()), np.sort(np.eye(5)[pi].ravel() / 5), atol=1e-9)   # a permutation / m
        Q2 = rng.normal(size=(17, 3))
        vb, Tbar, gamma, W2, Var = T.centroid_defect(P, Q2, Ginv)
        assert abs(vb + Var - W2) < 1e-10 * (1 + W2)                             # W² = V_bar² + Var_γ
        assert np.allclose(gamma.sum(1), 1 / 5) and np.allclose(gamma.sum(0), 1 / 17) and gamma.min() > -1e-12


def test_pure_barycentric_leak_dissipates_to_a_centroidal_configuration():
    """Along ϑ̇_i = −σ_m(ϑ_i − T̄_i): D⁺(½W²) ≤ −σ_m V_bar² (so W² is non-increasing, strictly while not centroidal), the
    flow converges to an equal-mass centroidal configuration of the prior (V_bar → 0), and that configuration is a
    fixed point of the equal-mass Lloyd iteration with the same W²."""
    rng = np.random.default_rng(1)
    Ginv = 1.0 / np.array([2.0, 1.0, 0.5])
    m, m0, sigma = 4, 24, 1.0
    Q = rng.normal(size=(m0, 3)) * np.array([1.0, 1.5, 0.4]) + np.array([0, 0, 1.5])
    P = rng.normal(size=(m, 3)) + np.array([0, 0, 1.5])
    dt = 1e-2
    W2_prev = None
    for k in range(4000):
        vb, Tbar, gamma, W2, Var = T.centroid_defect(P, Q, Ginv)
        if W2_prev is not None:
            assert W2 <= W2_prev + 1e-9, (k, W2, W2_prev)                             # non-increasing
        # the Dini bound at this state, checked over one explicit Euler step (first order in dt)
        P_next = P - dt * sigma * (P - Tbar)
        W2_next = T.optimal_plan(P_next, Q, Ginv)[1]
        assert (W2_next - W2) / dt <= -2 * sigma * vb + 5e-2 * sigma * vb + 1e-9, (k, (W2_next - W2) / dt, -2 * sigma * vb)
        P, W2_prev = P_next, W2
    vb_end = T.centroid_defect(P, Q, Ginv)[0]
    assert vb_end < 1e-8, vb_end
    P_lloyd, W2_lloyd = T.lloyd_equal_mass(P, Q, Ginv)
    assert np.max(np.abs(P_lloyd - P)) < 1e-6 and abs(W2_lloyd - W2_prev) < 1e-6


def test_induced_coupling_bound_and_plan_stability():
    """√C_bar ≤ b + 2η with C_bar := (1/m)Σ‖ϑ*_{π*(i)} − T̄_i‖²_{G⁻¹} and b := W(μ*, ν); near the ideal, the optimal plan
    from μ to ν is the one from μ* to ν (LP stability), so C_bar equals the centroid defect b_bar² of the ideal."""
    rng = np.random.default_rng(2)
    Ginv = 1.0 / np.array([2.0, 1.0, 0.5])
    m, m0 = 4, 12
    for trial in range(30):
        P_star = rng.normal(size=(m, 3)); Q = rng.normal(size=(m0, 3)) * 1.5
        b2 = T.optimal_plan(P_star, Q, Ginv)[1]; b = np.sqrt(b2)
        bbar2 = T.centroid_defect(P_star, Q, Ginv)[0]
        for scale in (1.0, 0.3, 1e-3):
            P = P_star[rng.permutation(m)] + scale * rng.normal(size=(m, 3))
            pi, eta2 = T.optimal_assignment(P, P_star, Ginv); eta = np.sqrt(eta2)
            Tbar = T.barycentric_projection(P, Q, Ginv)[0]
            d = P_star[pi] - Tbar
            C_bar = float(np.einsum('ik,k,ik->', d, Ginv, d) / m)
            assert np.sqrt(C_bar) <= b + 2 * eta + 1e-9, (trial, scale, np.sqrt(C_bar), b + 2 * eta)
            if scale == 1e-3:
                assert abs(C_bar - bbar2) < 1e-6 * (1 + bbar2), (C_bar, bbar2)
                assert bbar2 <= b2 + 1e-12


def test_lyapunov_inequality_along_the_barycentric_leak_closed_loop():
    """Theorem C′(iii) along the closed loop with the barycentric leak toward a 24-atom prior:
    D⁺V_OT ≤ −k_e‖e‖² + ‖e‖(ε_max + Kλ_Gη²) + V_*η − σ_mη² + σ_mη√C_bar at every recorded state (ε_max = 0, V_* = 0),
    with C_bar from the plan the law uses."""
    rng = np.random.default_rng(3)
    m, n, m0 = 3, 1, 24
    ctrl0 = T.TransportController(m, n, ke=3.0, g_a=2.0, g_c=1.0, g_s=0.5, **BOX)
    P_star = T.sample_lifted_ideal(m, n, ctrl0, rng, min_sep=0.8)
    prior = np.vstack([P_star[i] + rng.normal(0, 0.15, size=(m0 // m, 3)) for i in range(m)])
    prior = ctrl0.hard_clip(prior)
    assert ctrl0.admissible_ideal(prior) and ctrl0.admissible_ideal(P_star)
    ctrl = T.TransportController(m, n, ke=3.0, g_a=2.0, g_c=1.0, g_s=0.5, leak='ot_bary', sigma_m=0.5, prior=prior, **BOX)
    plant = T.LiftedRepresentablePlant(P_star, amplitude=1.5)
    P0 = ctrl.hard_clip(P_star[[1, 2, 0]] + rng.normal(0, 0.3, size=P_star.shape))
    res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0), P0, T=12.0, dt=2e-3, record_every=25)
    worst = -np.inf
    for t, x, P in zip(res['t'], res['x'], res['P']):
        rep = ctrl.lyapunov_report(x, t, plant, P, P_star)
        pi = rep['pi']; eta = rep['W']; e = rep['e_norm']
        Tbar = T.barycentric_projection(P, prior, ctrl.Ginv_diag)[0]
        d = P_star[pi] - Tbar
        C_bar = float(np.einsum('ik,k,ik->', d, ctrl.Ginv_diag, d) / m)
        bound = -ctrl.ke * e ** 2 + e * ctrl.K * ctrl.lam_G * eta ** 2 - ctrl.sigma_m * eta ** 2 + ctrl.sigma_m * eta * np.sqrt(C_bar)
        worst = max(worst, rep['Vdot'] - bound)
        assert rep['Vdot'] <= bound + 1e-9 * (1 + abs(bound)), (t, rep['Vdot'], bound)
    assert worst <= 1e-9
