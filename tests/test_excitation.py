"""
Checks for PARADIGM.md Lemma E.0 (explicit strict Lyapunov function for the skew-symmetric mixed error model) and
for the explicit constants of Theorem E (revised), using src/core/excitation.py.

Lemma E.0 is a pointwise statement: at every (t, e, ζ), p_lo‖w‖² ≤ V ≤ p_hi‖w‖² and V̇ ≤ −γ‖w‖² along
ė = −k e − Ψζ, ζ̇ = Ψᵀe.  The tests evaluate V̇ by the chain rule with the analytic Q̇, on random states and times,
for several regressors; then they check the envelope along trajectories, the ordering explicit rate ≤ exact Floquet
rate, the ω⁻² scaling of Remark E.1, the slow-adaptation limit, and the nonlinear envelope of Theorem E on the
implemented closed loop.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core import excitation as X
from src.core import transport as T

PER = 2 * np.pi


def _sinusoid_regressor(rng, n, q, n_tones=3, period=PER):
    freqs = (2 * np.pi / period) * np.arange(1, n_tones + 1)
    A = rng.normal(size=(n, q, n_tones)); B = rng.normal(size=(n, q, n_tones))

    def Psi(t):
        return (A * np.sin(freqs * t)[None, None, :] + B * np.cos(freqs * t)[None, None, :]).sum(-1)

    def Psid(t):
        return (A * freqs * np.cos(freqs * t)[None, None, :] - B * freqs * np.sin(freqs * t)[None, None, :]).sum(-1)
    return Psi, Psid


def _remark_E1_regressor(omega):
    """m = n = 1, ϑ* = (1.5, 0.3, 0.8), G = I, x_m = 1.5 sin ωt (PARADIGM.md Remark E.1)."""
    P_star = np.array([[1.5, 0.3, 0.8]])
    reg = X.TheoremERegressor(P_star, np.ones(3), lambda t: np.array([1.5 * np.sin(omega * t)]),
                              lambda t: np.array([1.5 * omega * np.cos(omega * t)]))
    return reg


def test_lemma_E0_pointwise_inequalities_hold_on_random_states():
    """p_lo‖w‖² ≤ V ≤ p_hi‖w‖² and V̇ ≤ −γ‖w‖² at random (t, w), for the closed-form ε and for the optimised ε,
    on four random PE regressors of different shapes and gains."""
    rng = np.random.default_rng(3)
    for trial, (n, q, scale, k) in enumerate([(2, 3, 1.0, 3.0), (1, 3, 0.1, 2.0), (3, 2, 5.0, 10.0), (2, 5, 0.3, 0.5)]):
        Psi0, Psid0 = _sinusoid_regressor(rng, n, q)
        Psi = lambda t, s=scale, f=Psi0: s * f(t)
        Psid = lambda t, s=scale, f=Psid0: s * f(t)
        tg = np.linspace(0, PER, 120)
        Psi_M, Psi_D = X.sup_norms(Psi, Psid, tg)
        Tw = PER
        alpha = X.pe_level(Psi, Tw, tg, n_nodes=120)
        assert alpha > 0
        for consts in (X.strict_lyapunov_constants(k, alpha, Tw, Psi_M, Psi_D),
                       X.optimise_eps(k, alpha, Tw, Psi_M, Psi_D, n_grid=60)):
            SL = X.StrictLyapunov(Psi, Psid, k, Tw, consts, n_nodes=120)
            assert consts['gamma'] > 0 and consts['p_lo'] > 0
            if not np.isnan(consts['gamma_cf']):
                assert consts['gamma'] >= consts['gamma_cf'] * (1 - 1e-12)            # the closed-form lower bound
                assert consts['p_lo'] >= 0.25 - 1e-12 and consts['p_hi'] <= 0.75 + 1e-12
            for _ in range(60):
                t = rng.uniform(0, 3 * PER); e = rng.normal(size=n); z = rng.normal(size=q)
                V, Vd, g, lo, hi = SL.check_point(t, e, z)
                assert lo - 1e-12 <= V <= hi + 1e-12, (trial, V, lo, hi)
                assert Vd <= -g + 1e-10 * (1 + abs(Vd)), (trial, Vd, -g)


def test_lemma_E0_envelope_and_floquet_ordering():
    """Along (LTV): ‖w(t)‖ ≤ overshoot·e^{−rate t}‖w(0)‖ at every sample, and the explicit rate never exceeds the
    exact Floquet rate (it is a valid lower bound); the ratio is reported in results/explicit_rate/."""
    rng = np.random.default_rng(5)
    Psi, Psid = _sinusoid_regressor(rng, 2, 3)
    k = 3.0
    tg = np.linspace(0, PER, 120)
    rate, best, _ = X.explicit_rate(Psi, Psid, k, [PER / 2, PER, 2 * PER], tg, n_nodes=120)
    c = X.strict_lyapunov_constants(k, best['alpha'], best['T'], best['Psi_M'], best['Psi_D'])
    fl = X.floquet_rate(Psi, k, PER, n_steps=2000)
    assert 0 < c['rate'] <= fl, (c['rate'], fl)
    for w0 in (np.array([0, 0, 1.0, 0, 0]), np.array([1.0, -1.0, 0.3, 0.2, 0.1])):
        ts, nrm = X.simulate_ltv(Psi, k, w0, 60.0, 1e-2)
        env = c['overshoot'] * np.exp(-c['rate'] * ts) * np.linalg.norm(w0)
        assert np.all(nrm <= env * (1 + 1e-9))


def test_remark_E1_omega_scaling_of_the_explicit_rate():
    """Remark E.1's instance: under x_m = 1.5 sin ωt the exact Floquet rate falls like ω⁻² (review round 3).  The
    explicit bound stays below the exact rate at every ω, and for large ω it falls like ω⁻² as well, because Ψ_D ∝ ω
    dominates c₁ there: the bound explains the ω⁻² law.  The ratio exact/explicit decreases with ω and stays bounded."""
    k = 3.0
    oms = (16.0, 32.0, 64.0, 128.0)
    fl, ex = [], []
    for om in oms:
        reg = _remark_E1_regressor(om)
        per = PER / om
        tg = np.linspace(0, per, 80)
        _, best, _ = X.explicit_rate(reg.Psi, reg.Psidot, k, [per, 2 * per], tg, n_nodes=80)
        c = X.optimise_eps(k, best['alpha'], best['T'], best['Psi_M'], best['Psi_D'], n_grid=80)
        f = X.floquet_rate(reg.Psi, k, per, n_steps=1000)
        assert c['rate'] <= f, (om, c['rate'], f)
        fl.append(f); ex.append(c['rate'])
    fl, ex = np.array(fl), np.array(ex)
    slope_fl = np.polyfit(np.log(oms), np.log(fl), 1)[0]
    slope_ex = np.polyfit(np.log(oms), np.log(ex), 1)[0]
    assert -2.3 < slope_fl < -1.7, slope_fl
    assert -2.3 < slope_ex < -1.7, slope_ex
    ratio = fl / ex
    assert np.all(np.diff(ratio) <= 0) and ratio[-1] < 30, ratio


def test_slow_adaptation_limit_ratio_is_bounded():
    """In the mean-field normalisation Ψ_G = m^{−1/2}Ψ(I⊗G^{1/2}) is small for large m.  Scaling G → 0 at fixed
    ideal, the ratio (Floquet rate)/(explicit rate) decreases towards a constant of order 3κ_Ψ(1 + Ψ_D/(kΨ_M))² (Remark
    E.0a), while in the fast regime it is large.  Checked on the m = 3, n = 2 ideal of the Theorem E tests, driven by a
    two-tone reference of period 20π, with G = g·diag(2,2,1); this is a no-growth check of the ratio."""
    rng = np.random.default_rng(12)
    m, n = 3, 2
    ctrl = T.TransportController(m, n, ke=6.0, g_a=30.0, g_c=30.0, g_s=15.0, a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)
    P_star = T.sample_lifted_ideal(m, n, ctrl, rng, min_sep=0.9)
    plant = T.LiftedRepresentablePlant(P_star, amplitude=2.5)                    # period 20π (tones 1 and 1.3)
    per = 20 * np.pi
    ratios = []
    for g in (10.0, 1.0, 0.1, 0.01):
        G = g * np.array([2.0, 2.0, 1.0])
        G_diag = np.concatenate([np.full(n, G[0]), np.full(n, G[1]), [G[2]]])
        reg = X.TheoremERegressor(P_star, G_diag, plant.x_ref, plant.x_ref_dot)
        tg = np.linspace(0, per, 60)
        _, best, _ = X.explicit_rate(reg.Psi, reg.Psidot, ctrl.ke, [per], tg, n_nodes=800)     # converged PE level over 20π
        c = X.optimise_eps(ctrl.ke, best['alpha'], best['T'], best['Psi_M'], best['Psi_D'], n_grid=80)
        f = X.floquet_rate(reg.Psi, ctrl.ke, per, n_steps=3000)
        assert c['rate'] <= f
        ratios.append(f / c['rate'])
    ratios = np.array(ratios)
    assert ratios[0] > ratios[-1]                           # fast regime is more conservative than slow
    assert ratios[-1] < 1.5 * ratios[-2]                    # no growth of the ratio as the gains shrink (a no-growth check, not a saturation check)


def test_theorem_E_explicit_nonlinear_envelope_on_the_implemented_loop():
    """
    Theorem E (revised) with explicit constants on the implemented closed loop (Remark E.1's single-kernel instance
    with G = 0.3·I, where the Gramian is well conditioned and the constants are representable): for a start with
    ‖w(t₀)‖ ≤ r₀ in the w-coordinates, with L_J estimated numerically on D, the trajectory satisfies
    ‖w(t)‖ ≤ overshoot·e^{−rate_nl t}‖w(0)‖ at every sample over a horizon of min(1/rate_nl, 8000) s, which is about
    0.3 e-folds of the explicit envelope here (the exact rate is far larger, so the state itself decays by more).  The explicit constants are
    conservative (r₀ is of order 1e-5 here), so the check is on the inequality, not on sharpness.
    """
    rng = np.random.default_rng(7)
    m, n = 1, 1
    g = 0.3
    ctrl = T.TransportController(m, n, ke=3.0, g_a=g, g_c=g, g_s=g, a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)
    P_star = np.array([[1.5, 0.3, 0.8]])
    assert ctrl.admissible_ideal(P_star, strict=True)
    plant = T.LiftedRepresentablePlant(P_star, amplitude=1.5, freqs=[1.0], phases=[0.0])        # x_m = 1.5 sin t
    per = PER
    p = 2 * n + 1
    reg = X.TheoremERegressor(P_star, ctrl.G_diag, plant.x_ref, plant.x_ref_dot)
    tg = np.linspace(0, per, 80)
    _, best, _ = X.explicit_rate(reg.Psi, reg.Psidot, ctrl.ke, [per], tg, n_nodes=160)
    delta = min(ctrl.a_max - ctrl.eps - abs(P_star[0, 0]), ctrl.c_max - ctrl.eps - abs(P_star[0, 1]),
                P_star[0, 2] - ctrl.s_min - ctrl.eps, ctrl.s_max - ctrl.eps - P_star[0, 2])
    r = float(delta)                                                                         # ρ = ∞ (ε ≡ 0)

    def field(z, t):
        e = z[:n]; P = P_star + z[n:].reshape(m, p)
        x = plant.x_ref(t) + e
        edot = plant.f(x, t) - ctrl.ke * e - T.f_hat(x, P)
        v, _ = ctrl.velocity(x, e, P)
        return np.concatenate([edot, v.reshape(-1)])

    L_J = X.jacobian_lipschitz_estimate(field, n + m * p, r, tg, rng, n_samples=150)
    eps = X.optimise_eps(ctrl.ke, best['alpha'], best['T'], best['Psi_M'], best['Psi_D'], n_grid=120)['eps']
    cst = X.theorem_E_constants(ctrl.ke, best['alpha'], best['T'], best['Psi_M'], best['Psi_D'], m, ctrl.G_diag, r, L_J, eps=eps)
    assert cst['r0_z'] > 1e-8 and cst['rate_nonlinear'] > 1e-6, cst
    S = np.concatenate([np.ones(n), np.tile(1 / np.sqrt(m * ctrl.G_diag), m)])                # w = S z (diagonal)
    z0 = rng.normal(size=n + m * p); z0 *= 0.9 * cst['r0_z'] / np.linalg.norm(z0)
    w0 = S * z0
    assert np.linalg.norm(w0) <= cst['r0_w'] * (1 + 1e-9)
    P0 = P_star + z0[n:].reshape(m, p); x0 = plant.x_ref(0.0) + z0[:n]
    T_end = min(1.0 / cst['rate_nonlinear'], 8.0e3)          # the explicit rate is slow (≈ 3.7e-5): 8000 s is ≈ 0.3 e-folds of the envelope
    res = T.rk4_closed_loop(plant, ctrl, x0, P0, T=T_end, dt=0.02, record_every=50)
    norms = []
    for t, x, P in zip(res['t'], res['x'], res['P']):
        z = np.concatenate([x - plant.x_ref(t), (P - P_star).reshape(-1)])
        w = S * z
        env = cst['overshoot'] * np.exp(-cst['rate_nonlinear'] * t) * np.linalg.norm(w0)
        assert np.linalg.norm(w) <= env * (1 + 1e-6), (t, np.linalg.norm(w), env)
        norms.append(np.linalg.norm(w))
    assert norms[-1] < 0.5 * norms[0]                       # it does converge (the exact rate is far larger)


# ---------------------------------------------------------------- Lemma E.0′ (the slow-adaptation regime)
def test_lemma_E0prime_sinusoidal_regressor_certificate_and_sharpness():
    """
    Ψ(t) = (sin t, cos t), n = 1, q = 2, T = 2π: M̄ = I/2, so λ_min = 1/2 and the peak factors are ρ = ρ_D = ρ_DD = 2.
    The second-order averaged matrix N₂ vanishes (identically, Lemma E.0′); the conditions (A)–(D) hold from k ≈ 10 on;
    the certificate sup_{t,x} V̇/V ≤ −γ holds at k = 20, 50 and 200; the guaranteed rate γ/2 is below the exact Floquet
    rate and within a factor 2 of it, approaching 1/(1 − δ) as k grows.
    """
    Psi = lambda t: np.array([[np.sin(t), np.cos(t)]])
    Psid = lambda t: np.array([[np.cos(t), -np.sin(t)]])
    S = X.SlowRegimeLyapunov(Psi, Psid, 2 * np.pi, N=512)
    assert abs(S.lam_min - 0.5) < 1e-12 and abs(S.rho - 2) < 1e-9 and abs(S.rho_D - 2) < 1e-9 and abs(S.rho_DD - 2) < 1e-6
    assert np.linalg.norm(S.N2) < 1e-12 and S.p1_meas <= 0.5 * 2 * np.pi * S.Psi_M ** 2 + 1e-12
    assert S.threshold(0.2) < 12
    ratios = []
    for k in (20.0, 50.0, 200.0):
        c = S.constants(k, 0.2)
        assert c['ok'], c['conds']
        sup, minus_gamma = S.certificate(k, 0.2, n_t=128)
        assert sup <= minus_gamma * (1 - 1e-9)
        exact = X.floquet_rate(Psi, k, 2 * np.pi)
        assert c['rate'] <= exact
        ratios.append(exact / c['rate'])
    assert ratios[0] < 2.0 and ratios[-1] < 1.3 and ratios[0] > ratios[1] > ratios[2]        # → 1/(1 − δ) = 1.25


def test_lemma_E0prime_pays_for_ill_conditioning_once_on_the_theorem_E_ideal():
    """
    The slow-adaptation ideal of the Theorem E tests (m = 3, n = 2, two-tone reference of period 20π, G = 0.01·diag(2,2,1);
    κ(M̄) ≈ 8·10³, peak factors ρ ≈ 70, ρ_D ≈ 4·10³, ρ_DD ≈ 4·10⁵).  At k = 400 and δ = 0.2 the conditions of Lemma E.0′
    hold, the certificate holds on a 128-point grid, the guaranteed rate is within a factor 1.5 of the exact Floquet
    rate, and it exceeds the optimised rate of Lemma E.0 by more than four orders of magnitude.
    """
    rng = np.random.default_rng(12)
    m, n = 3, 2
    ctrl = T.TransportController(m, n, ke=6.0, g_a=30.0, g_c=30.0, g_s=15.0, a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)
    P_star = T.sample_lifted_ideal(m, n, ctrl, rng, min_sep=0.9)
    plant = T.LiftedRepresentablePlant(P_star, amplitude=2.5)
    per = 20 * np.pi
    g = 0.01
    G_diag = np.concatenate([np.full(n, 2 * g), np.full(n, 2 * g), [g]])
    reg = X.TheoremERegressor(P_star, G_diag, plant.x_ref, plant.x_ref_dot)
    S = X.SlowRegimeLyapunov(reg.Psi, reg.Psidot, per, N=2048)
    assert S.lam_max / S.lam_min > 1e3 and S.rho > 10 and np.linalg.norm(S.N2) < 1e-15 * S.lam_max
    k, delta = 400.0, 0.2
    c = S.constants(k, delta)
    assert c['ok'], c['conds']
    sup, minus_gamma = S.certificate(k, delta, n_t=128)
    assert sup <= minus_gamma * (1 - 1e-9)
    exact = X.floquet_rate(reg.Psi, k, per)
    assert c['rate'] <= exact <= 1.5 * c['rate']
    tg = np.linspace(0, per, 60)
    _, best, _ = X.explicit_rate(reg.Psi, reg.Psidot, k, [per], tg, n_nodes=800)
    e0 = X.optimise_eps(k, best['alpha'], best['T'], best['Psi_M'], best['Psi_D'], n_grid=80)
    assert c['rate'] > 1e4 * e0['rate']


def test_corollary_E0prime_2_envelope_on_the_implemented_loop():
    """
    Corollary E.0′.2 (Theorem E in the slow regime) on Remark E.1's single-kernel instance with G = 0.3·I at k_e = 40,
    δ = 0.2 (above Lemma E.0′'s threshold k₀ ≈ 25 for this regressor): from a start with ‖w(0)‖ ≤ r′₀ the implemented
    closed loop satisfies ‖w(t)‖ ≤ C′e^{−λ′t}‖w(0)‖ at every sample over 1.5 e-folds of the explicit envelope, and the
    explicit nonlinear rate λ′ exceeds Theorem E's rate with Lemma E.0 by more than a factor 10 (it is within a factor 4
    of the exact linear rate).  The step is dt = 0.005: the fixed-step RK4 loop has a discretisation floor of
    ‖w‖ ≈ 6·10⁻⁶ at dt = 0.02 (∝ dt⁴), above the explicit ball r′₀ ≈ 1.4·10⁻⁶, and ≈ 2·10⁻⁸ at dt = 0.005.
    """
    rng = np.random.default_rng(7)
    m, n = 1, 1
    g, ke, delta = 0.3, 40.0, 0.2
    ctrl = T.TransportController(m, n, ke=ke, g_a=g, g_c=g, g_s=g, a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)
    P_star = np.array([[1.5, 0.3, 0.8]])
    plant = T.LiftedRepresentablePlant(P_star, amplitude=1.5, freqs=[1.0], phases=[0.0])
    per = PER
    p = 2 * n + 1
    reg = X.TheoremERegressor(P_star, ctrl.G_diag, plant.x_ref, plant.x_ref_dot)
    S = X.SlowRegimeLyapunov(reg.Psi, reg.Psidot, per, N=1024)
    assert S.constants(ke, delta)['ok']
    r = float(min(ctrl.a_max - ctrl.eps - abs(P_star[0, 0]), ctrl.c_max - ctrl.eps - abs(P_star[0, 1]),
                  P_star[0, 2] - ctrl.s_min - ctrl.eps, ctrl.s_max - ctrl.eps - P_star[0, 2]))
    tg = np.linspace(0, per, 80)

    def field(z, t):
        e = z[:n]; P = P_star + z[n:].reshape(m, p); x = plant.x_ref(t) + e
        edot = plant.f(x, t) - ctrl.ke * e - T.f_hat(x, P); v, _ = ctrl.velocity(x, e, P)
        return np.concatenate([edot, v.reshape(-1)])
    L_J = X.jacobian_lipschitz_estimate(field, n + m * p, r, tg, rng, n_samples=150)
    cst = X.slow_regime_theorem_E_constants(S, ke, delta, m, ctrl.G_diag, r, L_J)
    _, best, _ = X.explicit_rate(reg.Psi, reg.Psidot, ke, [per], tg, n_nodes=160)
    e0 = X.optimise_eps(ke, best['alpha'], best['T'], best['Psi_M'], best['Psi_D'], n_grid=120)
    old = X.theorem_E_constants(ke, best['alpha'], best['T'], best['Psi_M'], best['Psi_D'], m, ctrl.G_diag, r, L_J, eps=e0['eps'])
    assert cst['rate_nonlinear'] > 10 * old['rate_nonlinear']
    assert X.floquet_rate(reg.Psi, ke, per) < 4 * cst['rate_nonlinear']
    Sdiag = np.concatenate([np.ones(n), np.tile(1 / np.sqrt(m * ctrl.G_diag), m)])             # w = S z
    z0 = rng.normal(size=n + m * p); z0 *= 0.9 * cst['r0_z'] / np.linalg.norm(z0)
    w0 = Sdiag * z0
    assert np.linalg.norm(w0) <= cst['r0_w'] * (1 + 1e-9)
    P0 = P_star + z0[n:].reshape(m, p); x0 = plant.x_ref(0.0) + z0[:n]
    T_end = 1.5 / cst['rate_nonlinear']                                                        # 1.5 e-folds of the envelope
    res = T.rk4_closed_loop(plant, ctrl, x0, P0, T=T_end, dt=0.005, record_every=200)
    worst = 0.0
    for t, x, P in zip(res['t'], res['x'], res['P']):
        z = np.concatenate([x - plant.x_ref(t), (P - P_star).reshape(-1)])
        w = Sdiag * z
        worst = max(worst, np.linalg.norm(w) / (cst['C'] * np.exp(-cst['rate_nonlinear'] * t) * np.linalg.norm(w0)))
    assert worst <= 1.0, worst
    assert np.linalg.norm(Sdiag * np.concatenate([res['x'][-1] - plant.x_ref(res['t'][-1]), (res['P'][-1] - P_star).reshape(-1)])) < 0.5 * np.linalg.norm(w0)



def test_corollary_E0prime_3_dither_covariance_bound_against_the_periodic_lyapunov_solution():
    """Corollary E.0′.3 on Ψ = (sin t, cos t), Σ = √(2·10⁻³) I: the explicit bound on limsup E‖ζ‖² dominates the exact
    periodic stationary covariance of (LTV_τ) at k = 20, 50, 200 for δ = 0.1 and 0.2, the ratio decreases toward
    1/(1 − δ) (within 6% at k = 200), and the exact value agrees with the averaged Langevin covariance k·tr X (M̄X + XM̄ = ΣΣᵀ)
    within 2/k² (the deviation is ≈ 1.25/k²)."""
    Psi = lambda t: np.array([[np.sin(t), np.cos(t)]])
    Psid = lambda t: np.array([[np.cos(t), -np.sin(t)]])
    S = X.SlowRegimeLyapunov(Psi, Psid, 2 * np.pi, N=1024)
    Sig = np.eye(2) * np.sqrt(2e-3)
    for delta in (0.2, 0.1):
        ratios = []
        for k in (20.0, 50.0, 200.0):
            b = X.dither_covariance_bound(S, k, delta, Sig)
            eps = 1.0 / k
            L_fun = lambda t, eps=eps: -eps * Psi(t) + eps ** 2 * Psid(t) - eps ** 3 * Psi(t) @ (Psi(t).T @ Psi(t))   # the slaving map of E.0′.2
            ex = X.exact_stationary_covariance(Psi, k, Sig, 2 * np.pi, L_fun=L_fun)
            assert b['bound_zeta'] >= ex['zeta_mean'] and b['bound_zeta'] >= ex['zeta_max']
            assert b['bound_w'] >= ex['w_max'] > 0                                   # the exact stationary E‖w‖² = E‖e − Lζ‖² is below its bound
            assert b['bound_w'] > 1e3 * ex['w_max']                                  # ... by more than three orders of magnitude (the bound is coarse, factor 1/β)
            assert abs(ex['e_mean'] / (b['Sigma_F2'] / (2 * k)) - 1) < 1e-4            # exact E‖e‖² = ‖Σ‖_F²/(2k) = tr(M̄X)/k (round-7 review)
            assert abs(ex['zeta_mean'] / b['averaged_zeta'] - 1) < 2.0 / k ** 2                 # the deviation is ≈ 1.25/k²
            assert b['nu_grid'] <= b['nu_explicit'] * (1 + 1e-9)
            ratios.append(b['bound_zeta'] / ex['zeta_mean'])
        assert ratios[0] > ratios[1] > ratios[2] > 1 / (1 - delta)
        assert ratios[2] < 1.06 / (1 - delta)


def test_lemma_E0dprime_quasi_periodic_regressor_certificate_and_rate():
    """Lemma E.0″ on the quasi-periodic regressor Ψ = (sin t, cos √2 t) (not periodic; PE), window T = 4π: the conditions hold
    at k = 20, 50, 200 for δ = 0.1 and 0.2; the certificate sup V̇/V ≤ −γ holds on a grid of 400 times in [0, 300]; the
    guaranteed rate γ/2 is below the measured decay rate of (LTV) (the top Lyapunov exponent over 200 time units) and within
    a factor 2.6 of it, the ratio improving with k; the exponential averages obey TṀ̄ = M̄ − M and Ṗ₁ = 2(M − M̄) on the grid;
    and the residual second-order term N̄₂ costs δ₂ = ε‖N̄₂‖/λ ≤ 5% at k ≥ 20."""
    r2 = np.sqrt(2.0)
    Psi = lambda t: np.array([[np.sin(t), np.cos(r2 * t)]])
    Psid = lambda t: np.array([[np.cos(t), -r2 * np.sin(r2 * t)]])
    E = X.ExpAverageLyapunov(Psi, Psid, 4 * np.pi, T_end=300.0, dt=0.005)
    assert 0.3 < E.lam_min < 0.5 and E.lam_max < 1.0                           # the long-time mean of M is ½I; the window average is below it
    # the differential identities of the exponential averages on the grid (central differences)
    i = slice(2000, 20000)
    dMbar = (E.Mbar[2001:20001] - E.Mbar[1999:19999]) / (2 * E.dt)
    assert np.max(np.abs(E.T * dMbar - (E.Mbar[i] - E.Ms[i]))) < 1e-4
    dP1 = (E.P1[2001:20001] - E.P1[1999:19999]) / (2 * E.dt)
    assert np.max(np.abs(dP1 - 2 * (E.Ms[i] - E.Mbar[i]))) < 2e-4
    prev = None
    for k in (20.0, 50.0, 200.0):
        exact = X.lyapunov_exponent_rate(Psi, k, 200.0)
        for delta in (0.2, 0.1):
            c = E.constants(k, delta)
            assert c['ok'], (k, delta, c['conds'])
            assert c['delta2'] <= 0.05
            sup, mg = E.certificate(k, delta, n_t=400)
            assert sup <= mg, (k, delta, sup, mg)
            assert c['rate'] < exact < 2.6 * c['rate'], (k, delta, c['rate'], exact)
        ratio = exact / E.constants(k, 0.1)['rate']
        if prev is not None:
            assert ratio < prev
        prev = ratio
