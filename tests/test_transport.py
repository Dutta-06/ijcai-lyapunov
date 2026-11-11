"""Machine checks of PARADIGM.md Theorems A–E (see the verification map there)."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core import transport as T
from src.core.controller import AdaptiveController

BOX = dict(a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)


def _ctrl(m, n, **kw):
    base = dict(ke=3.0, g_a=2.0, g_c=1.0, g_s=0.5, **BOX)
    base.update(kw)
    return T.TransportController(m, n, **base)


def _oracle_setup(seed, m=6, n=2, noise=0.3, permute=True, **kw):
    rng = np.random.default_rng(seed)
    ctrl = _ctrl(m, n, **kw)
    P_star = T.sample_lifted_ideal(m, n, ctrl, rng)
    assert ctrl.admissible_ideal(P_star)
    plant = T.LiftedRepresentablePlant(P_star)
    P0, pi0 = T.perturb_and_permute(P_star, rng, noise, permute)
    P0 = ctrl.hard_clip(P0)
    x0 = rng.uniform(-1, 1, size=n)
    return ctrl, plant, P_star, P0, x0


# ---------------------------------------------------------------- Theorem A
def test_gaussian_constants_are_suprema():
    s_min, a_max = 0.4, 5.0
    cst = T.gaussian_constants(s_min, a_max)
    rng = np.random.default_rng(0)
    n = 2
    worst_H = 0.0
    for _ in range(3000):
        x = rng.normal(size=n) * 3; c = rng.normal(size=n) * 3; s = rng.uniform(s_min, 3.0)
        # numerical Hessian of φ in (c, s)
        def phi(cs):
            return np.exp(-np.sum((x - cs[:n]) ** 2) / (2 * cs[n] ** 2))
        z = np.concatenate([c, [s]]); h = 1e-4; H = np.zeros((n + 1, n + 1))
        for i in range(n + 1):
            for j in range(n + 1):
                ei = np.zeros(n + 1); ej = np.zeros(n + 1); ei[i] = h; ej[j] = h
                H[i, j] = (phi(z + ei + ej) - phi(z + ei - ej) - phi(z - ei + ej) + phi(z - ei - ej)) / (4 * h * h)
        worst_H = max(worst_H, np.linalg.norm(H, 2))
    assert worst_H <= cst['H_phi'] * (1 + 1e-3)
    assert abs(cst['h_cc'] - 1.0) < 1e-6 and 0.8 < cst['h_cs'] < 0.9 and 1.2 < cst['h_ss'] < 1.25


def test_identification_with_classical_law():
    rng = np.random.default_rng(1)
    m, n = 5, 2
    gW, gc, gs = 3.0, 0.7, 0.4
    cl = AdaptiveController(m=m, n=n, ke=2.0, gamma_W=gW, gamma_c=gc, gamma_sigma=gs,
                            W_max=1e6, c_max=1e6, sigma_min=0.3, sigma_max=1e3, proj_epsilon=0.1)
    tr = T.TransportController(m, n, ke=2.0, g_a=m * gW, g_c=gc / m, g_s=gs / m,
                               a_max=1e6, c_max=1e6, s_min=0.3, s_max=1e3, eps=0.1)
    for _ in range(100):
        x = rng.normal(size=n); e = rng.normal(size=n)
        W = rng.normal(size=(m, n)); c = rng.normal(size=(m, n)); s = rng.uniform(0.5, 2.0, size=m)
        Wd, cd, sd = cl.adaptation_laws(x, e, W, c, s)
        P = np.hstack([m * W, c, s[:, None]])
        v, _ = tr.raw_velocity(x, e, P)
        assert np.allclose(v[:, :n], m * Wd, atol=1e-10)
        assert np.allclose(v[:, n:2 * n], cd.reshape(m, n), atol=1e-10)
        assert np.allclose(v[:, 2 * n], sd, atol=1e-10)
        assert np.allclose(T.f_hat(x, P), cl.kernel.f_hat(x, W, c, s))


def test_identification_projected_rates_per_kernel():
    """
    Theorem A(i), corrected: the PROJECTED rates of (L) coincide with the classical rates only if the
    classical law projects each kernel separately, W_i onto ‖W_i‖ ≤ a_max/m with layer ε/m, c_i onto
    ‖c_i‖ ≤ c_max and σ_i onto [s_min, s_max] with layer ε.  Checked with active projections.
    """
    from src.core.projection import ProjectionOperator, BandwidthProjection
    rng = np.random.default_rng(21)
    m, n = 4, 2
    a_max, c_max, s_min, s_max, eps = 2.0, 1.5, 0.4, 2.0, 0.3
    gW, gc, gs = 1.5, 0.8, 0.5
    tr = T.TransportController(m, n, ke=2.0, g_a=m * gW, g_c=gc / m, g_s=gs / m,
                               a_max=a_max, c_max=c_max, s_min=s_min, s_max=s_max, eps=eps)
    cl = AdaptiveController(m=m, n=n, ke=2.0, gamma_W=gW, gamma_c=gc, gamma_sigma=gs,       # all classical projections
                            W_max=1e6, c_max=1e6, sigma_min=1e-3, sigma_max=1e3, proj_epsilon=0.1)  # inactive: raw rates
    pW = ProjectionOperator(a_max / m, eps / m); pc = ProjectionOperator(c_max, eps); ps = BandwidthProjection(s_min, s_max, eps)
    pW_bad = ProjectionOperator(a_max / m, eps)                     # layer not rescaled
    worst, worst_bad, active = 0.0, 0.0, 0
    for _ in range(400):
        x = rng.normal(size=n) * 1.5; e = rng.normal(size=n)
        a = rng.normal(size=(m, n)); a *= (a_max - rng.uniform(0, eps, size=m))[:, None] / np.linalg.norm(a, axis=1)[:, None]
        c = rng.normal(size=(m, n)); c *= (c_max - rng.uniform(0, eps, size=m))[:, None] / np.linalg.norm(c, axis=1)[:, None]
        s = rng.uniform(s_min, s_max, size=m)
        P = np.hstack([a, c, s[:, None]]); W = a / m
        v, _ = tr.velocity(x, e, P)
        Wd_raw, cd_raw, sd_raw = cl.adaptation_laws(x, e, W, c, s)       # huge radii: raw classical rates
        Wd = np.array([pW.project_vector(W[i], Wd_raw[i]) for i in range(m)])
        Wd_bad = np.array([pW_bad.project_vector(W[i], Wd_raw[i]) for i in range(m)])
        cd = np.array([pc.project_vector(c[i], cd_raw.reshape(m, n)[i]) for i in range(m)])
        sd = ps.project_vector(s, sd_raw)
        worst = max(worst, np.abs(v[:, :n] - m * Wd).max(), np.abs(v[:, n:2 * n] - cd).max(), np.abs(v[:, 2 * n] - sd).max())
        worst_bad = max(worst_bad, np.abs(v[:, :n] - m * Wd_bad).max())
        active += int(np.any(np.abs(v[:, :n] - m * Wd_raw) > 1e-9))
    assert active > 100, "projection must actually be active in this test"
    assert worst < 1e-10, worst
    assert worst_bad > 1e-3, "without rescaling the layer the rates differ (the corrected statement is sharp)"
    # 'only if' (m ≥ 2): for EVERY joint Frobenius ball (W_max, ε_W) some configuration distinguishes it from (L)
    # (two-case construction, review round 2): r_in := (a_max − ε_K)/m, ρ := W_max − ε_W.
    #   ρ < √2 r_in : two equal rows of norm r with √2 r ∈ (ρ, min(√2 r_in, W_max)] — per-kernel inactive, joint active;
    #   ρ ≥ √2 r_in : one row with norm in (r_in, min(a_max/m, ρ))            — per-kernel active, joint inactive.
    m2, amax2, epsK = 2, 2.0, 0.3
    tr2 = T.TransportController(m2, 1, ke=1.0, g_a=2.0, g_c=1.0, g_s=1.0, a_max=amax2, c_max=3.0, s_min=0.5, s_max=2.0, eps=epsK)
    r_in = (amax2 - epsK) / m2
    for Wmax, epsW in ((0.707, 0.05), (0.707, 0.3), (1.0, 0.6), (0.9, 0.8), (3.0, 0.5)):
        rho = Wmax - epsW
        if rho < np.sqrt(2) * r_in:
            hi = min(np.sqrt(2) * r_in, Wmax); r = 0.5 * (rho + hi) / np.sqrt(2)
            Wrows = np.array([[r], [r]])
        else:
            Wrows = np.array([[0.5 * (r_in + min(amax2 / m2, rho))], [0.0]])
        P2 = np.hstack([m2 * Wrows, np.zeros((m2, 1)), np.ones((m2, 1))])     # c = 0, s = 1
        v2, _ = tr2.velocity(np.array([0.0]), np.array([1.0]), P2)            # x = c: φ = 1, outward update
        raw_W = (2.0 / m2) * np.ones((m2, 1))                                  # γ_W φ e, γ_W = g_a/m
        joint = ProjectionOperator(Wmax, epsW).project_matrix(Wrows, raw_W)
        assert np.abs(m2 * joint[:, 0] - v2[:, 0]).max() > 1e-6, (Wmax, epsW, m2 * joint[:, 0], v2[:, 0])


# ---------------------------------------------------------------- Theorem B
def test_theorem_B_inequality_along_trajectories():
    for seed in range(3):
        ctrl, plant, P_star, P0, x0 = _oracle_setup(seed)
        reports = []
        T.rk4_closed_loop(plant, ctrl, x0, P0, T=6.0, dt=1e-3, record_every=5,
                          callback=lambda x, P, t: reports.append(ctrl.lyapunov_report(x, t, plant, P, P_star)))
        for r in reports:
            assert abs(r['Vdot_identity_check']) < 1e-9 * (1 + abs(r['Vdot']))
            assert r['cross'] <= 1e-10
            assert np.all(np.abs(r['rem']) <= r['rem_bound'] * (1 + 1e-6) + 1e-12)
            assert r['Vdot'] <= r['bound_B'] + 1e-9
            assert r['inside']


def test_dini_derivative_of_V_OT_bounded_by_active_branch():
    ctrl, plant, P_star, P0, x0 = _oracle_setup(4, noise=0.25)
    res = T.rk4_closed_loop(plant, ctrl, x0, P0, T=4.0, dt=5e-4, record_every=1)
    Vs = np.array([0.5 * res['e_norm'][k] ** 2 + 0.5 * ctrl.ot_distance(res['P'][k], P_star)[1] for k in range(len(res['t']))])
    dt = res['t'][1] - res['t'][0]
    fd = (Vs[1:] - Vs[:-1]) / dt                                   # forward difference ≈ D⁺V_OT
    bounds = np.array([ctrl.lyapunov_report(res['x'][k], res['t'][k], plant, res['P'][k], P_star)['bound_B'] for k in range(len(res['t']) - 1)])
    viol = fd - bounds
    assert np.max(viol) < 2e-3 * (1 + np.max(np.abs(bounds))), f"D+V_OT exceeded bound by {np.max(viol)}"


class _MovingIdealPlant:
    """
    f(x,t) = (1/m)Σ a*_j φ(x; c*_j(t), s*_j) + ε(t), with the ideal centres rotating rigidly
    (c*_j(t) = R(ωt) c*_j(0), so ‖c*_j(t)‖ is constant and the ideal stays in the inner set) and
    ε(t) = ε₀ sin(3t)·1 (non-representable part, ‖ε‖ ≤ ε₀√n).
    """

    def __init__(self, P0_star, omega=0.6, eps0=0.05):
        self.P0 = np.array(P0_star, float); self.n = (self.P0.shape[1] - 1) // 2
        self.omega, self.eps0 = omega, eps0

    def _R(self, t):
        cth, sth = np.cos(self.omega * t), np.sin(self.omega * t)
        return np.array([[cth, -sth], [sth, cth]])

    def P_star(self, t):
        P = self.P0.copy(); n = self.n
        P[:, n:2 * n] = self.P0[:, n:2 * n] @ self._R(t).T
        return P

    def P_star_dot(self, t):
        D = np.zeros_like(self.P0); n = self.n; w = self.omega
        Rd = w * np.array([[-np.sin(w * t), -np.cos(w * t)], [np.cos(w * t), -np.sin(w * t)]])
        D[:, n:2 * n] = self.P0[:, n:2 * n] @ Rd.T
        return D

    def f(self, x, t):
        return T.f_hat(x, self.P_star(t)) + self.eps0 * np.sin(3 * t) * np.ones(self.n)

    def x_ref(self, t):
        return np.array([np.sin(t), np.cos(0.7 * t)])

    def x_ref_dot(self, t):
        return np.array([np.cos(t), -0.7 * np.sin(0.7 * t)])


def test_theorem_B_moving_ideal_nonrepresentable_switches_and_active_projection():
    """
    Theorem B in the cases the basic tests do not reach: time-varying ideal (drift term V_*·W),
    approximation error ε ≠ 0, optimal-assignment switches, and an active projection (network weights
    started in the boundary layer).  At every sample: chain-rule identity, cross ≤ 0, per-particle
    remainder bound, V̇_π ≤ bound (B) for the active π; forward-difference Dini quotient of V_OT ≤ (B).
    """
    m, n = 4, 2
    box = dict(a_max=2.0, c_max=2.5, s_min=0.4, s_max=2.0, eps=0.2)
    ctrl = T.TransportController(m, n, ke=3.0, g_a=2.0, g_c=0.3, g_s=0.3, **box)
    P0_star = np.array([[1.6, 0.0, 1.2, 0.0, 0.9], [-1.5, 0.3, -1.2, 0.0, 1.1],
                        [0.0, 1.6, 0.0, 1.0, 0.8], [0.4, -1.5, 0.0, -1.0, 1.0]])
    plant = _MovingIdealPlant(P0_star, omega=0.8, eps0=0.05)
    assert ctrl.admissible_ideal(plant.P_star(0.0))
    rng = np.random.default_rng(31)
    P0 = P0_star.copy(); P0[:, :n] *= (box['a_max'] - 0.05) / np.linalg.norm(P0[:, :n], axis=1)[:, None]   # weights in the layer
    P0[:, n:2 * n] += rng.normal(0, 0.2, size=(m, n))
    P0 = ctrl.hard_clip(P0)
    x0 = plant.x_ref(0.0) + np.array([0.3, -0.2])
    reps, pis, active = [], [], []

    def cb(x, P, t):
        r = ctrl.lyapunov_report(x, t, plant, P, plant.P_star(t), P_star_dot=plant.P_star_dot(t))
        reps.append(r); pis.append(tuple(r['pi']))
        v_raw, _ = ctrl.raw_velocity(x, x - plant.x_ref(t), P)
        active.append(bool(np.any(np.abs(ctrl.project(P, v_raw) - v_raw) > 1e-12)))
    res = T.rk4_closed_loop(plant, ctrl, x0, P0, T=16.0, dt=5e-4, record_every=4, callback=cb)
    switches = sum(pis[k] != pis[k - 1] for k in range(1, len(pis)))
    assert switches >= 1, "the optimal assignment must switch in this test"
    assert sum(active) >= 10, "the projection must be active in this test"
    for r in reps:
        assert abs(r['Vdot_identity_check']) < 1e-9 * (1 + abs(r['Vdot']))
        assert r['cross'] <= 1e-10
        assert np.all(np.abs(r['rem']) <= r['rem_bound'] * (1 + 1e-6) + 1e-12)
        assert r['Vdot'] <= r['bound_B'] + 1e-9
        assert r['V_star_t'] > 0
    # forward-difference Dini quotient of V_OT along the run (integrator error budget 2%)
    V = np.array([r['V'] for r in reps]); tt = np.array([r['t'] for r in reps]); B = np.array([r['bound_B'] for r in reps])
    fd = np.diff(V) / np.diff(tt)
    assert np.max(fd - B[:-1]) < 2e-2 * (1 + np.max(np.abs(B))), np.max(fd - B[:-1])


def test_permuted_ideal_labelled_vs_OT():
    """Network initialised as a permutation of the ideal (+ tiny noise): labelled error large, OT error tiny,
    disturbance tiny; the OT bound is predictive, the labelled one is not."""
    ctrl, plant, P_star, P0, x0 = _oracle_setup(5, noise=0.01, permute=True)
    r = ctrl.lyapunov_report(x0, 0.0, plant, P0, P_star)
    assert r['identity_err'] > 10 * r['W'], (r['identity_err'], r['W'])
    actual_excess = abs(r['rem'].sum()) / ctrl.m
    ot_excess_bound = ctrl.excess_bound(r['e_norm'], r['W'])
    assert actual_excess <= ot_excess_bound * (1 + 1e-6) + 1e-12
    labelled_excess_bound = ctrl.excess_bound(r['e_norm'], r['identity_err'])
    assert labelled_excess_bound > 50 * ot_excess_bound


# ---------------------------------------------------------------- Theorem C
def test_theorem_C_i_exact_dissipation_pure_leak():
    rng = np.random.default_rng(6)
    m, n = 6, 2
    ctrl = _ctrl(m, n, g_a=1.0, g_c=1.0, g_s=1.0, leak='ot', sigma_m=0.8)
    prior = T.sample_lifted_ideal(m, n, ctrl, rng); ctrl.prior = prior
    P = ctrl.hard_clip(T.perturb_and_permute(prior, rng, 0.8)[0])
    W0 = T.w2(P, prior, ctrl.Ginv_diag)
    dt, steps = 1e-3, 3000
    switches = 0; prev = None
    for k in range(steps):
        v, tau = ctrl.raw_velocity(np.zeros(n), np.zeros(n), P)       # e = 0 → gradient term vanishes
        if prev is not None and np.any(tau != prev):
            switches += 1
        prev = tau
        P = P + dt * v                                                  # Euler on a linear ODE between switches
    Wt = T.w2(P, prior, ctrl.Ginv_diag)
    assert abs(Wt / W0 - np.exp(-0.8 * dt * steps)) < 2e-3, (Wt / W0, np.exp(-0.8 * dt * steps), switches)
    assert switches == 0, "no-switching lemma: the active assignment stays optimal along the pure leak"


def test_lemma_C_ii_induced_coupling_bound():
    rng = np.random.default_rng(8)
    m, n = 5, 2
    ctrl = _ctrl(m, n)
    R = ctrl.diameter()
    exact_hits = sharp_hits = 0
    for trial in range(300):
        P_star = T.sample_lifted_ideal(m, n, ctrl, rng)
        if trial % 2 == 0:
            P0 = T.sample_lifted_ideal(m, n, ctrl, rng)                       # generic prior
        else:
            P0 = ctrl.hard_clip(P_star[rng.permutation(m)] + rng.normal(0, 0.05, size=P_star.shape))  # prior ≈ permuted ideal (large gap)
        gap, b2 = T.assignment_gap(P_star, P0, ctrl.Ginv_diag)
        noise = rng.uniform(0.001, 0.5) if trial % 2 == 0 else rng.uniform(1e-4, 2e-3)
        P = ctrl.hard_clip(P_star[rng.permutation(m)] + rng.normal(0, noise, size=P_star.shape))
        pi_s, eta2 = T.optimal_assignment(P, P_star, ctrl.Ginv_diag); eta = np.sqrt(eta2)
        tau, _ = T.optimal_assignment(P, P0, ctrl.Ginv_diag)
        C = T.assignment_cost(P_star[pi_s], P0, tau, ctrl.Ginv_diag)
        assert C <= b2 + 4 * eta * R + 1e-9
        b = np.sqrt(b2)
        assert np.sqrt(C) <= b + 2 * eta + 1e-9                              # R-free form: √C ≤ b + 2η
        if 4 * eta * R < gap:
            assert abs(C - b2) < 1e-9; exact_hits += 1
        if 4 * eta * (b + eta) < gap:                                        # sharper gap condition of C(ii)
            assert abs(C - b2) < 1e-9; sharp_hits += 1
    assert exact_hits > 20
    assert sharp_hits >= exact_hits


def test_theorem_C_iv_exponential_estimate():
    """
    An instance satisfying EVERY hypothesis of Theorem C(iv), asserted explicitly:
    small constants (s_min = 1, a_max = 1 ⇒ K ≈ 1.5), m = 2 well-separated ideal particles, prior = ideal
    + small noise (b small, large assignment gap), start = RELABELLED ideal + noise (W(0) ≈ 0.15).
    Checks V_OT(t) ≤ e^{−αt}V_OT(0) + β/α with the theorem's α, β at every sample.
    """
    m, n = 2, 1
    box = dict(a_max=1.0, c_max=2.0, s_min=1.0, s_max=2.0, eps=0.15)
    ctrl = T.TransportController(m, n, ke=4.0, g_a=1.0, g_c=1.0, g_s=1.0, leak='ot', sigma_m=1.0, **box)
    P_star = np.array([[0.5, -1.0, 1.3], [-0.6, 1.0, 1.6]])        # (a, c, s) per particle
    assert ctrl.admissible_ideal(P_star)
    rng = np.random.default_rng(3)
    prior = ctrl.hard_clip(P_star + rng.normal(0, 0.02, size=P_star.shape)); ctrl.prior = prior
    assert ctrl.admissible_ideal(prior)
    plant = T.LiftedRepresentablePlant(P_star, amplitude=1.0)
    gap, b2 = T.assignment_gap(P_star, prior, ctrl.Ginv_diag)
    R = ctrl.diameter(); eta_bar = 0.9 * gap / (4 * R)
    alpha = min(ctrl.ke, ctrl.sigma_m / 2); beta = (ctrl.sigma_m / 2) * b2           # ε = 0, V_* = 0
    assert 4 * eta_bar * (np.sqrt(b2) + eta_bar) < gap, "gap condition of version 3: 4η̄(b + η̄) < Δ_gap"
    P0 = ctrl.hard_clip(P_star[[1, 0]] + rng.normal(0, 0.08, size=P_star.shape))      # relabelled start
    x0 = plant.x_ref(0.0) + 0.05
    V0 = 0.5 * float(np.sum((x0 - plant.x_ref(0.0)) ** 2)) + 0.5 * ctrl.ot_distance(P0, P_star)[1]
    # --- hypotheses of Theorem C(iv), asserted (no escape hatch) ---
    assert 8 * ctrl.K ** 2 * ctrl.lam_G ** 2 * eta_bar ** 2 <= ctrl.ke * ctrl.sigma_m, "gain condition"
    assert max(V0, beta / alpha) < eta_bar ** 2 / 2, "initial level inside the gap-protected set"
    res = T.rk4_closed_loop(plant, ctrl, x0, P0, T=12.0, dt=1e-3, record_every=10)
    Vs = np.array([0.5 * res['e_norm'][k] ** 2 + 0.5 * ctrl.ot_distance(res['P'][k], P_star)[1] for k in range(len(res['t']))])
    env = np.exp(-alpha * res['t']) * Vs[0] + beta / alpha
    assert np.all(Vs <= env * (1 + 1e-6) + 1e-12), np.max(Vs - env)
    assert Vs[-1] < 0.2 * Vs[0], "the estimate is not vacuous: V_OT actually contracts"
    assert np.all(Vs < eta_bar ** 2 / 2), "trajectory stays in the gap-protected set"


# ---------------------------------------------------------------- Theorem F (transport cost of adaptation)
def test_lemma_P1_projection_is_nonexpansive_and_tight_velocity_bound():
    """
    Lemma P.1: ‖Proj_K(ϑ, Gy)‖_{G⁻¹} ≤ ‖y‖_G (block-diagonal G, ramps in [0,1]), with equality when no factor
    is active; and the velocity bound ‖G∇g‖_{G⁻¹} ≤ λ_G^{1/2} L₁ ‖e‖ of Theorem F(i), which is nearly attained.
    """
    rng = np.random.default_rng(41)
    m, n = 5, 2
    ctrl = _ctrl(m, n, a_max=2.0, c_max=2.0, s_min=0.3, s_max=2.0, eps=0.3)
    N = lambda v: float(np.sqrt(np.einsum('ik,k,ik->', v, ctrl.Ginv_diag, v) / m))
    for _ in range(2000):
        P = np.hstack([rng.normal(size=(m, n)), rng.normal(size=(m, n)), rng.uniform(0.3, 2.0, size=(m, 1))])
        P[:, :n] *= (2.0 - rng.uniform(0, 0.4, size=m))[:, None] / np.linalg.norm(P[:, :n], axis=1)[:, None]
        P[:, n:2 * n] *= (2.0 - rng.uniform(0, 0.4, size=m))[:, None] / np.linalg.norm(P[:, n:2 * n], axis=1)[:, None]
        x = rng.normal(size=n); e = rng.normal(size=n)
        raw, _ = ctrl.raw_velocity(x, e, P)
        proj = ctrl.project(P, raw)
        assert N(proj) <= N(raw) * (1 + 1e-12)
        if np.allclose(proj, raw):
            assert abs(N(proj) - N(raw)) < 1e-12
        assert N(raw) <= np.sqrt(ctrl.lam_G) * ctrl.L1 * np.linalg.norm(e) * (1 + 1e-12)
    # tight instance: m = n = 1, small a_max, x at the centre, g_a = λ_G
    t = T.TransportController(1, 1, ke=1.0, g_a=3.0, g_c=1.0, g_s=1.0, a_max=0.05, c_max=2.0, s_min=1.0, s_max=2.0, eps=0.01)
    Pt = np.array([[0.0, 0.0, 1.0]]); e = np.array([1.0])
    v, _ = t.raw_velocity(np.array([0.0]), e, Pt)
    ratio = np.sqrt(np.sum(v ** 2 * t.Ginv_diag)) / (np.sqrt(t.lam_G) * t.L1 * 1.0)
    assert ratio >= 0.99, ratio


def test_theorem_F_leak_clause_along_ot_sigma_mod():
    """Along (L_σ) with the OT leak: W(μ_t1, μ_t0) ≤ ∫ (λ_G^{1/2}L₁‖e‖ + σ_m W(μ_r, μ⁰)) dr on a grid of pairs."""
    rng = np.random.default_rng(42)
    m, n = 4, 2
    ctrl = _ctrl(m, n, leak='ot', sigma_m=0.6)
    P_star = T.sample_lifted_ideal(m, n, ctrl, rng)
    ctrl.prior = ctrl.hard_clip(P_star[rng.permutation(m)] + rng.normal(0, 0.3, size=P_star.shape))
    plant = T.LiftedRepresentablePlant(P_star)
    P0 = ctrl.hard_clip(P_star + rng.normal(0, 0.4, size=P_star.shape))
    res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0) + 0.3, P0, T=6.0, dt=1e-3, record_every=5)
    t, e, P = res['t'], res['e_norm'], res['P']
    Wp = np.array([ctrl.ot_distance(Pk, ctrl.prior)[1] ** 0.5 for Pk in P])
    integrand = np.sqrt(ctrl.lam_G) * ctrl.L1 * e + ctrl.sigma_m * Wp
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(t))])
    for i in range(0, len(t), 60):
        for j in range(i + 30, len(t), 90):
            assert ctrl.ot_distance(P[j], P[i])[1] ** 0.5 <= (cum[j] - cum[i]) * (1 + 1e-3) + 1e-9


def test_theorem_F_wasserstein_length_bounded_by_integrated_error():
    """Particle path length (≥ Wasserstein length) ≤ λ_G^{1/2} L₁ ∫₀ᵗ ‖e‖ ds along the projected law (L)."""
    for seed in range(2):
        ctrl, plant, P_star, P0, x0 = _oracle_setup(seed, noise=0.4)
        res = T.rk4_closed_loop(plant, ctrl, x0, P0, T=6.0, dt=1e-3, record_every=1)
        t, e, P = res['t'], res['e_norm'], res['P']
        steps = np.sqrt(np.einsum('kij,j,kij->k', np.diff(P, axis=0), ctrl.Ginv_diag, np.diff(P, axis=0)) / ctrl.m)
        length = np.concatenate([[0.0], np.cumsum(steps)])
        int_e = np.concatenate([[0.0], np.cumsum(0.5 * (e[1:] + e[:-1]) * np.diff(t))])
        bound = np.sqrt(ctrl.lam_G) * ctrl.L1 * int_e
        assert np.all(length <= bound * (1 + 1e-3) + 1e-9), np.max(length - bound)
        for k in range(0, len(t), len(t) // 6):
            assert ctrl.ot_distance(P[k], P0)[1] ** 0.5 <= length[k] + 1e-9


# ---------------------------------------------------------------- Theorem E
def test_theorem_E_linearisation_and_lyapunov_identity():
    """
    The proof of Theorem E linearises the IMPLEMENTED closed loop at (e, ϑ) = (0, ϑ*∘π).  Check, by
    finite differences of the implemented vector field at several times and for a relabelled
    equilibrium, that the Jacobian is A(t) = [[−k_e I, −(1/m)Ψ_π], [(I_m⊗G)Ψ_πᵀ, 0]], and that
    Q = diag(I/2, (I_m⊗G⁻¹)/(2m)) gives QA + AᵀQ = diag(−k_e I, 0).
    """
    rng = np.random.default_rng(12)
    m, n = 3, 2
    ctrl = _ctrl(m, n, ke=6.0, g_a=30.0, g_c=30.0, g_s=15.0)
    P_star = T.sample_lifted_ideal(m, n, ctrl, rng, min_sep=0.9)
    assert ctrl.admissible_ideal(P_star, strict=True)
    plant = T.RichReferencePlant(P_star, amplitude=2.5, n_tones=6)
    perm = np.array([1, 2, 0]); Peq = P_star[perm]                  # a relabelled equilibrium
    p = 2 * n + 1

    def field(z, t):
        e = z[:n]; P = Peq + z[n:].reshape(m, p)
        x = plant.x_ref(t) + e
        edot = plant.f(x, t) - ctrl.ke * e - T.f_hat(x, P)
        v, _ = ctrl.velocity(x, e, P)
        return np.concatenate([edot, v.reshape(-1)])

    for t in (0.0, 3.7, 11.2):
        z0 = np.zeros(n + m * p)
        assert np.abs(field(z0, t)).max() < 1e-12, "relabelled ideal is an equilibrium"
        h = 1e-6
        Jfd = np.column_stack([(field(z0 + h * np.eye(len(z0))[k], t) - field(z0 - h * np.eye(len(z0))[k], t)) / (2 * h)
                               for k in range(len(z0))])
        Psi = np.concatenate([g for g in T.grad_Phi(plant.x_ref(t), Peq)], axis=1)          # n × mp
        IG = np.kron(np.eye(m), np.diag(ctrl.G_diag))
        A = np.block([[-ctrl.ke * np.eye(n), -Psi / m], [IG @ Psi.T, np.zeros((m * p, m * p))]])
        assert np.abs(Jfd - A).max() < 1e-5 * (1 + np.abs(A).max()), np.abs(Jfd - A).max()
        Q = np.block([[np.eye(n) / 2, np.zeros((n, m * p))], [np.zeros((m * p, n)), np.kron(np.eye(m), np.diag(ctrl.Ginv_diag)) / (2 * m)]])
        S = Q @ A + A.T @ Q
        target = np.zeros_like(S); target[:n, :n] = -ctrl.ke * np.eye(n)
        assert np.abs(S - target).max() < 1e-12 * (1 + np.abs(A).max())


def test_theorem_E_matched_error_converges_labelled_does_not():
    # Qualitative check only. The rate depends on Ψ(·) itself, not only on its PE level (it is not monotone in G or
    # k_e, and it falls like ω⁻² under time-rescaling of the reference); here the smallest Gramian eigenvalue is
    # ≈ 2e-3, so convergence is slow.
    # Criterion: decay of the MATCHED error to a fraction of its initial value while the LABELLED error stays O(1).
    rng = np.random.default_rng(12)
    m, n = 3, 2
    ctrl = _ctrl(m, n, ke=6.0, g_a=30.0, g_c=30.0, g_s=15.0)
    P_star = T.sample_lifted_ideal(m, n, ctrl, rng, min_sep=0.9)
    plant = T.RichReferencePlant(P_star, amplitude=2.5, n_tones=6)
    P0 = ctrl.hard_clip(P_star[[1, 2, 0]] + rng.normal(0, 0.05, size=P_star.shape))   # cyclic permutation + noise
    x0 = plant.x_ref(0.0)
    res = T.rk4_closed_loop(plant, ctrl, x0, P0, T=150.0, dt=1e-3, record_every=1000)
    W = np.array([ctrl.ot_distance(P, P_star)[1] ** 0.5 for P in res['P']])
    lab = np.array([np.sqrt(np.einsum('ik,k,ik->', P - P_star, ctrl.Ginv_diag, P - P_star) / m) for P in res['P']])
    assert W[-1] < 0.25 * W[0], (W[0], W[-1])
    k = len(W) // 3
    slope = np.polyfit(res['t'][k:], np.log(W[k:] + 1e-14), 1)[0]
    assert slope < 0
    assert lab[-1] > 0.5 * lab[0]




# ---------------------------------------------------------------- Corollary B.5 along (L) vs (L_σ)
def test_corollary_B5_fails_along_lsigma_and_its_leak_variant_holds():
    """
    Corollary B.5 integrates (B), which holds along (L) only.  Along (L_σ) with every hypothesis of Theorem C(iv)
    asserted and sup W ≤ η̄ measured, the leak keeps e = O(b) and the plain B.5 bound fails; the (L_σ) variant
    with the extra (σ_m b²/4)/k_e term holds (review round 2, AB.1).
    """
    box = dict(a_max=1.0, c_max=2.0, s_min=1.0, s_max=2.0, eps=0.15)
    P_star = np.array([[0.5, -0.3, 1.4]])
    b = 0.01
    prior = P_star.copy(); prior[:, 0] -= b                                   # W(μ*, μ⁰) = b  (G = I)
    ctrl = T.TransportController(1, 1, ke=4.0, g_a=1.0, g_c=1.0, g_s=1.0, leak='ot', sigma_m=1.0, prior=prior, **box)
    assert ctrl.admissible_ideal(P_star) and ctrl.admissible_ideal(prior)
    gap, b2 = T.assignment_gap(P_star, prior, ctrl.Ginv_diag)                 # single atom: gap = ∞
    plant = T.LiftedRepresentablePlant(P_star, amplitude=1.0)                 # ε_max = 0, V_* = 0
    alpha = min(ctrl.ke, ctrl.sigma_m / 2); beta = (ctrl.sigma_m / 2) * b2
    eta_bar = 1.02 * np.sqrt(2 * beta / alpha)
    assert (not np.isfinite(gap)) or 4 * eta_bar * (np.sqrt(b2) + eta_bar) < gap
    assert 8 * ctrl.K ** 2 * ctrl.lam_G ** 2 * eta_bar ** 2 <= ctrl.ke * ctrl.sigma_m
    assert beta / alpha < eta_bar ** 2 / 2
    res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0), P_star.copy(), T=60.0, dt=2e-3, record_every=5)
    t, e = res['t'], res['e_norm']
    W = np.array([np.sqrt(ctrl.ot_distance(P, P_star)[1]) for P in res['P']])
    assert W.max() <= eta_bar
    k1 = np.searchsorted(t, 20.0)
    Tw = t[-1] - t[k1]
    rms = np.sqrt(np.trapezoid(e[k1:] ** 2, t[k1:]) / Tw)
    V1 = 0.5 * e[k1] ** 2 + 0.5 * W[k1] ** 2
    ke, lam = ctrl.ke, ctrl.lam_G
    Dbar = min(ctrl.K * lam * eta_bar ** 2, ctrl.L0 + ctrl.L1 * np.sqrt(lam) * eta_bar)
    plain = Dbar / (2 * ke) + np.sqrt(Dbar ** 2 / (4 * ke ** 2) + V1 / (ke * Tw))
    variant = Dbar / (2 * ke) + np.sqrt(Dbar ** 2 / (4 * ke ** 2) + (ctrl.sigma_m * b2 / 4) / ke + V1 / (ke * Tw))
    assert rms > 3 * plain, (rms, plain)                                      # B.5 is false along (L_σ)
    assert rms <= variant, (rms, variant)                                     # the leak variant holds


def test_theorem_F_leak_clause_is_tight_when_the_leak_pays_for_all_motion():
    """
    Theorem F (ii′) with e ≡ 0 (review rounds 2 and 3): two identical ideal kernels, prior weights (1.5, 0.5), start
    a = (3, −1).  The OT leak conserves a₁ + a₂, hence f̂ ≡ f and e ≡ 0, while the measure moves by W ≈ 1.05.  The
    leak clause must hold and be nearly tight; the (L) reading (motion paid by tracking error) must fail.
    """
    ctrl = _ctrl(2, 1, leak='ot', sigma_m=0.5)
    P_star = np.array([[1.0, 0.0, 1.0], [1.0, 0.0, 1.0]])
    ctrl.prior = np.array([[1.5, 0.0, 1.0], [0.5, 0.0, 1.0]])
    plant = T.LiftedRepresentablePlant(P_star)
    P0 = np.array([[3.0, 0.0, 1.0], [-1.0, 0.0, 1.0]])
    res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0), P0, T=10.0, dt=1e-3, record_every=10)
    t, e, P = res['t'], res['e_norm'], res['P']
    assert e.max() < 1e-10, e.max()
    W_move = ctrl.ot_distance(P[-1], P[0])[1] ** 0.5
    assert W_move > 1.05, W_move
    Wp = np.array([ctrl.ot_distance(Pk, ctrl.prior)[1] ** 0.5 for Pk in P])
    integrand = np.sqrt(ctrl.lam_G) * ctrl.L1 * e + ctrl.sigma_m * Wp
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(t))])
    ratios = [ctrl.ot_distance(P[j], P[0])[1] ** 0.5 / cum[j] for j in range(10, len(t), 10)]
    assert 0.99 <= max(ratios) <= 1 + 1e-3, max(ratios)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "passed")
