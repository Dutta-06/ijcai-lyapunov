"""
Checks for PARADIGM.md Proposition H.7 (frozen-measure structure of the exactly averaged closed loop) and
Remark H.0 (the reduction to (𝔄) has an O(1/k_e) relative error), using src/core/averaging.py.

All quantities are computed at a FROZEN measure μ: the exact periodic error e_μ of ė = r_μ(x_m + e) − k_e e is
integrated (RK4), and the averaged potential P_μ(ϑ) = ⟨e_μᵀΦ̃(x_m + e_μ; ϑ)⟩ is compared with the expansion
(1/k_e)⟨r_μᵀΦ̃_ϑ⟩_ρ + (1/k_e²)(U_even + U_odd).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core import langevin as L
from src.core import averaging as A

PER = 2 * np.pi


def _plant_1d():
    f = lambda X: np.sin(2 * X) + 0.3 * X
    J = lambda X: (2 * np.cos(2 * X) + 0.3)[:, :, None]
    xr = lambda t: np.array([1.5 * np.sin(t) + 0.5 * np.sin(2 * t)])
    xd = lambda t: np.array([1.5 * np.cos(t) + 1.0 * np.cos(2 * t)])
    return f, J, xr, xd


def _plant_2d(orient=+1):
    """f(x) = (sin x₂ + 0.3x₁, cos x₁ − 0.2x₂); elliptic reference (1.2 cos t, ±0.8 sin t), area 0.96π."""
    def f(X):
        return np.stack([np.sin(X[:, 1]) + 0.3 * X[:, 0], np.cos(X[:, 0]) - 0.2 * X[:, 1]], 1)

    def J(X):
        out = np.zeros((len(X), 2, 2))
        out[:, 0, 0] = 0.3; out[:, 0, 1] = np.cos(X[:, 1]); out[:, 1, 0] = -np.sin(X[:, 0]); out[:, 1, 1] = -0.2
        return out
    xr = lambda t: np.array([1.2 * np.cos(t), orient * 0.8 * np.sin(t)])
    xd = lambda t: np.array([-1.2 * np.sin(t), orient * 0.8 * np.cos(t)])
    return f, J, xr, xd


def _setup(n, seed=0):
    rng = np.random.default_rng(seed)
    net = L.SaturatedLiftedNetwork(n, a_max=3.0, c_max=2.0, s_min=0.3, s_max=1.5)
    p = 2 * n + 1
    return net, rng.normal(0, 1.0, (5, p)), rng.normal(0, 1.0, (4, p))


def test_frozen_measure_expansion_converges_at_rate_one_over_ke():
    """k_e²(P_μ − P₀) → U_even + U_odd with error O(1/k_e), for n = 1 and for both orientations of an n = 2 orbit."""
    cases = [(1, +1)] + [(2, o) for o in (+1, -1)]
    for n, orient in cases:
        net, Th, TE = _setup(n)
        f, J, xr, xd = _plant_1d() if n == 1 else _plant_2d(orient)
        orb = A.Orbit(xr, xd, PER, N=2000, sub=2)
        fm = A.FrozenMeasure(net, Th)
        Ue, Uo = A.expansion_terms(f, J, fm, orb, TE)
        errs = []
        for ke in (40.0, 80.0):
            corr = ke ** 2 * (A.averaged_potential(f, fm, orb, ke, TE) - A.leading_potential(f, fm, orb, ke, TE))
            errs.append(np.max(np.abs(corr - (Ue + Uo))) / np.max(np.abs(Ue + Uo)))
        assert errs[1] < 0.6 * errs[0], (n, orient, errs)                    # halves when k_e doubles
        assert errs[1] < (0.3 if n == 1 else 0.03), (n, orient, errs)
        if n == 1:
            # n = 1 converges slowly (small U_even); a Richardson step at k_e = 80, 160 makes the check sharp
            orb4 = A.Orbit(xr, xd, PER, N=2000, sub=4)
            C = {ke: ke ** 2 * (A.averaged_potential(f, fm, orb4, ke, TE) - A.leading_potential(f, fm, orb4, ke, TE))
                 for ke in (80.0, 160.0)}
            rich = 2 * C[160.0] - C[80.0]
            assert np.max(np.abs(rich - (Ue + Uo))) < 0.03 * np.max(np.abs(Ue + Uo)), np.max(np.abs(rich - (Ue + Uo)))


def test_circulation_kernel_antisymmetric_vanishes_in_1d_and_equals_the_enclosed_area_integral():
    """K is antisymmetric; K ≡ 0 for n = 1; for n = 2, K = (1/T₀) Σ_k ∫_D det[∇Φ̃_k(ϑ), ∇Φ̃_k(ϑ')] dA (Green),
    and reversing the orbit reverses K."""
    net1, Th1, _ = _setup(1)
    f, J, xr, xd = _plant_1d()
    K1 = A.circulation_kernel(net1, A.Orbit(xr, xd, PER, N=4000, sub=1), Th1, Th1)
    assert np.max(np.abs(K1)) < 1e-12
    net, Th, _ = _setup(2)
    _, _, xr, xd = _plant_2d(+1)
    K = A.circulation_kernel(net, A.Orbit(xr, xd, PER, N=4000, sub=1), Th, Th)
    assert np.max(np.abs(K + K.T)) < 1e-12 and np.max(np.abs(K)) > 1e-3
    _, _, xr_, xd_ = _plant_2d(-1)
    K_rev = A.circulation_kernel(net, A.Orbit(xr_, xd_, PER, N=4000, sub=1), Th, Th)
    assert np.max(np.abs(K + K_rev)) < 1e-12
    # area integral over the ellipse {(x/1.2)² + (y/0.8)² ≤ 1}: x = 1.2 ρ cos θ, y = 0.8 ρ sin θ, dA = 0.96 ρ dρ dθ
    g, wg = np.polynomial.legendre.leggauss(120)
    rho, wr = 0.5 * (g + 1), 0.5 * wg
    th = np.arange(400) * 2 * np.pi / 400
    R_, T_ = np.meshgrid(rho, th, indexing='ij')
    X = np.stack([1.2 * R_.ravel() * np.cos(T_.ravel()), 0.8 * R_.ravel() * np.sin(T_.ravel())], 1)
    w = (wr[:, None] * np.full(len(th), 2 * np.pi / 400)[None, :] * 0.96 * R_).ravel()
    _, dPhi = A.features(net, X, Th)                                          # (N, k, n, n): ∂_l Φ̃_k
    det = np.einsum('tikx,tjky->tijxy', dPhi, dPhi)
    area = np.einsum('t,tij->ij', w, det[..., 0, 1] - det[..., 1, 0]) / PER
    assert np.max(np.abs(area - K)) < 1e-6 * np.max(np.abs(K)), np.max(np.abs(area - K))


def test_shift_correction_is_even_and_lag_correction_is_odd_under_orbit_reversal():
    """Same μ and the same occupation measure ρ, orbit traversed both ways: U_even is unchanged, U_odd flips
    (U_odd ≡ 0 for n = 1)."""
    net, Th, TE = _setup(2)
    fm = A.FrozenMeasure(net, Th)
    out = {}
    for o in (+1, -1):
        f, J, xr, xd = _plant_2d(o)
        out[o] = A.expansion_terms(f, J, fm, A.Orbit(xr, xd, PER, N=2000, sub=1), TE)
    (Ue_p, Uo_p), (Ue_m, Uo_m) = out[+1], out[-1]
    scale = np.max(np.abs(np.concatenate([Ue_p, Uo_p])))
    assert np.max(np.abs(Ue_p - Ue_m)) < 1e-10 * scale
    assert np.max(np.abs(Uo_p + Uo_m)) < 1e-10 * scale and np.max(np.abs(Uo_p)) > 1e-3
    net1, Th1, TE1 = _setup(1)
    f, J, xr, xd = _plant_1d()
    _, Uo1 = A.expansion_terms(f, J, A.FrozenMeasure(net1, Th1), A.Orbit(xr, xd, PER, N=2000, sub=1), TE1)
    assert np.max(np.abs(Uo1)) < 1e-12


def test_exactly_averaged_dynamics_is_non_reciprocal_for_n1_and_n2():
    """
    The first-order correction is not a first variation: the cyclic sum of the antisymmetric part S − Sᵀ of the
    second variation is non-zero.  With S the exact second variation (finite differences on the exact averaged
    potential, k_e = 80) and S₀ = ⟨Φ̃_ϑᵀΦ̃_ϑ'⟩_ρ, the full matrix k_e²S + k_eS₀ is within 5% of the predicted coefficient
    S_even + K, and the cyclic sums agree to 3%, for n = 1 (where K ≡ 0, so the non-reciprocity comes from the state
    shift alone) and for n = 2.
    """
    for n in (1, 2):
        net, Th, TE = _setup(n)
        f, J, xr, xd = _plant_1d() if n == 1 else _plant_2d(+1)
        orb = A.Orbit(xr, xd, PER, N=2000, sub=2)
        fm = A.FrozenMeasure(net, Th)
        S_e, K = A.second_variation_formula(f, J, fm, orb, TE[:3])
        pred = A.cyclic_defect(S_e + K)
        assert abs(pred) > 0.1 * np.max(np.abs(S_e)), (n, pred)              # clearly non-reciprocal
        ke = 80.0
        Sx = A.exact_second_variation(f, fm, orb, ke, TE[:3])
        Phi, _ = A.features(net, orb.xm, TE[:3])
        S0 = np.einsum('tik,tjk->ij', Phi, Phi) / orb.N                       # leading order: −S₀/k_e (symmetric)
        assert np.max(np.abs(ke ** 2 * Sx + ke * S0 - (S_e + K))) < 0.05 * np.max(np.abs(S_e + K)), n
        exact = ke ** 2 * A.cyclic_defect(Sx)
        assert abs(exact - pred) < 0.03 * abs(pred), (n, exact, pred)


def test_averaged_drift_error_decays_like_one_over_ke_for_an_orbit_enclosing_area():
    """
    Remark H.0 for n = 2: at frozen μ, the exactly averaged drift G∇_ϑP_μ at the particles differs from the
    (𝔄) drift (1/k_e)G∇_ϑ⟨r_μᵀΦ̃_ϑ⟩ by a relative error that halves when k_e doubles; it does not vanish at finite k_e.
    """
    net, Th, _ = _setup(2, seed=3)
    f, J, xr, xd = _plant_2d(+1)
    orb = A.Orbit(xr, xd, PER, N=2000, sub=2)
    fm = A.FrozenMeasure(net, Th)
    m, p = Th.shape
    h = 1e-5
    E = np.vstack([Th + h * np.eye(p)[k][None, :] for k in range(p)] + [Th - h * np.eye(p)[k][None, :] for k in range(p)])
    errs = []
    for ke in (10.0, 20.0, 40.0, 80.0):
        P = A.averaged_potential(f, fm, orb, ke, E).reshape(2, p, m)
        P0 = A.leading_potential(f, fm, orb, ke, E).reshape(2, p, m)
        g_exact = (P[0] - P[1]).T / (2 * h)                                   # (m, p)
        g_lead = (P0[0] - P0[1]).T / (2 * h)
        errs.append(np.max(np.abs(g_exact - g_lead)) / np.max(np.abs(g_lead)))
    errs = np.array(errs)
    assert np.all(np.diff(errs) < 0), errs
    assert np.all(errs[1:] / errs[:-1] < 0.62) and np.all(errs[1:] / errs[:-1] > 0.4), errs
    assert errs[0] > 1e-3, errs                                              # non-zero at finite k_e


def test_oscillation_correction_is_odd_under_orbit_reversal():
    """The second-order averaging term ⟨D_ΘH·u⟩ of the intra-period parameter oscillation (PARADIGM.md Remark H.7a)
    changes sign when the orbit is reversed (same occupation measure), and is not small for the ellipse."""
    net, Th, _ = _setup(2)
    G = np.array([4.0, 4.0, 1.0, 1.0, 1.0])
    out = {}
    for o in (+1, -1):
        f, J, xr, xd = _plant_2d(o)
        out[o] = A.oscillation_correction(net, f, A.Orbit(xr, xd, PER, N=1000, sub=1), Th, G)
    scale = np.max(np.abs(out[+1]))
    assert scale > 1e-3
    assert np.max(np.abs(out[+1] + out[-1])) < 1e-6 * scale, np.max(np.abs(out[+1] + out[-1]))


def test_oscillation_correction_matches_a_one_period_flow():
    """
    Value check of the oscillation term (Remark H.7a), independent of the formula: integrate Θ̇ = H(Θ,t)/k_e over one
    period (RK4) from the near-identity preimage of Y, map back, and compare the measured mean drift beyond the
    leading term, k_e²[(P₁ − P₀)/T₀ − ⟨H⟩/k_e], with ⟨D_ΘH·u⟩ (k_e = 80; review round 3).
    """
    net, Th, _ = _setup(2, seed=3)
    f, J, xr, xd = _plant_2d(+1)
    G = np.ones(5)
    ke = 80.0

    def H(T, t):
        x = xr(t)
        r = f(x[None, :])[0] - net.f_hat(x, T)
        return G[None, :] * net.grad_g(x, r, T)

    def u0_and_mean(T, Nu=2000):
        Hs = np.array([H(T, t) for t in np.arange(Nu) * PER / Nu])
        hs = Hs - Hs.mean(0)
        I = np.concatenate([np.zeros((1,) + hs.shape[1:]), np.cumsum(0.5 * (hs[1:] + hs[:-1]), 0) * PER / Nu])
        return -I.mean(0), Hs.mean(0)                                         # u(T, 0) and ⟨H⟩(T)

    Th0 = Th + u0_and_mean(Th)[0] / ke
    X = Th0.copy()
    dt = PER / 2000
    for k in range(2000):
        t = k * dt
        k1 = H(X, t) / ke; k2 = H(X + dt / 2 * k1, t + dt / 2) / ke
        k3 = H(X + dt / 2 * k2, t + dt / 2) / ke; k4 = H(X + dt * k3, t + dt) / ke
        X = X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    P0 = Th0 - u0_and_mean(Th0)[0] / ke
    P1 = X - u0_and_mean(X)[0] / ke
    Ym = 0.5 * (P0 + P1)
    D = ke ** 2 * ((P1 - P0) / PER - u0_and_mean(Ym)[1] / ke)
    osc = A.oscillation_correction(net, f, A.Orbit(xr, xd, PER, N=2000, sub=1), Ym, G)
    assert np.max(np.abs(D - osc)) < 0.03 * np.max(np.abs(osc)), (np.max(np.abs(D - osc)), np.max(np.abs(osc)))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "passed")
