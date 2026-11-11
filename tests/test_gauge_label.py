"""
Proposition D‴.41 (PARADIGM.md §5d): measurement noise on the gauge. Checked here: the measured error's noise enters the loop along the
fields U_x of `label_noise.fields` (the drift with r + ξ minus the drift with r is (1/N)Σ_xξ(x)U_x, on the positions and the log-masses);
on M the fields are orthogonal to the gauge's tangent in the flow's metric (no diffusion along M); the Itô drift and the Wong–Zakai term
cancel at a point of M (the Stratonovich drift vanishes); the closed form's stationary point W_ℓ minimises Ψ_ℓ and is interior; the AR(1)
noise at ρ = 0 is the Euler scheme bit for bit, and Heun's scheme without noise is the explicit trapezoidal rule; the draws of the AR(1)
noise (per-step variance vN(1 − ρ)/((1 + ρ)dt), stationary innovations) and one Heun step with correlated noise off M against hand-coded
steps (one draw per step, the corrector reusing it, the predictor renormalised); the weighted and the covariance Wong–Zakai terms against
the field derivatives; the linear-response factor against the exact stationary moments of the residual–noise system, its continuum limits,
and the filtered covariance against the pairwise factors in an eigenbasis.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
from src.core import label_noise as LN                                   # noqa: E402
from src.core import mass_transfer as MT                                 # noqa: E402


def _loop(gamma=1e-2):
    net, base = E.setup()
    return net, base, MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)


def test_the_noise_enters_along_the_fields():
    net, base, loop = _loop()
    W = np.array([0.40, 0.25, 0.35]); b = 3 * W
    Th = E.exact_rep(net, W) + 0.03 * np.random.default_rng(0).standard_normal((3, 3))
    xi = np.random.default_rng(1).standard_normal(base.N)
    dT0, db0 = loop.drift(Th, b)
    dT1, db1 = loop.drift(Th, b, xi)
    UT, Ul = LN.fields(loop, Th, b)
    assert np.allclose(dT1 - dT0, np.einsum('x,xij->ij', xi, UT) / base.N, rtol=1e-12, atol=1e-15)
    assert np.allclose((db1 - db0) / b, xi @ Ul / base.N, rtol=1e-12, atol=1e-15)


def test_the_fields_are_orthogonal_to_the_gauge_in_the_flow_metric():
    # t ∈ T: δα_k = −(a_k/a′_k)δℓ_k, Σβ_kδℓ_k = 0; the metric Σ‖δϑ‖²_{G⁻¹} + Σ(β_k/γ)δℓ_k²
    net, base, loop = _loop()
    W = np.array([0.30, 0.35, 0.35]); b = 3 * W
    Th = E.exact_rep(net, W)
    a = base.features(Th)[0]; ap = 1 - (a / net.a_max) ** 2
    UT, Ul = LN.fields(loop, Th, b)
    for dl in (np.array([1.0, -1.0, 0.0]), np.array([0.5, 0.5, -1.0])):
        dl = dl - (b @ dl) / b.sum()
        t_alpha = -(a / ap) * dl
        inner = UT[:, :, 0] @ (t_alpha / E.G[0]) + Ul @ (b / loop.gamma * dl)
        scale = np.abs(UT[:, :, 0]).max() * np.abs(t_alpha).max() / E.G[0]
        assert np.abs(inner).max() < 1e-12 * scale


def test_the_wong_zakai_term_cancels_the_ito_drift():
    import exp_gauge_label as GL
    r = GL.a_job.__wrapped__((1e-2, (1 / 3, 1 / 3)))                     # uncached
    assert np.abs(r['ito']).max() > 1e-6
    assert np.abs(r['ito'] + r['wz']).max() < 1e-6 * np.abs(r['ito']).max()
    assert np.abs(r['D']).max() < 1e-12 * np.abs(r['ito']).max()


def test_the_closed_forms_target_is_interior_and_minimises_psi():
    net, base, loop = _loop()
    Th = E.exact_rep(net, np.ones(3) / 3)
    P = E.G[0] * np.mean(base.features(Th)[4] ** 2, axis=0)
    Wm = np.abs(E.A_BUMP) / net.a_max
    Wl = LN.w_label(P, Wm)
    assert abs(Wl.sum() - 1) < 1e-10 and np.all(Wl > Wm + 0.05)
    u = LN.selection(loop, E.exact_rep(net, Wl), 3 * Wl)
    assert np.ptp(u) < 1e-9 * np.abs(u).max()
    rng = np.random.default_rng(2)
    for _ in range(200):
        W = Wl + 0.05 * rng.standard_normal(3); W = W - (W.sum() - 1) / 3
        if np.all(W > Wm):
            assert LN.psi(W, P, Wm) >= LN.psi(Wl, P, Wm)
    # the closed form is −∂_βΨ_ℓ along M (Ψ_ℓ in β through W = β/3)
    W = np.array([0.40, 0.25, 0.35]); h = 1e-6
    u = LN.selection(loop, E.exact_rep(net, W), 3 * W)
    grad = np.array([(LN.psi(W + h * e / 3, P, Wm) - LN.psi(W - h * e / 3, P, Wm)) / (2 * h) for e in np.eye(3)])
    assert np.allclose(u, -grad, rtol=1e-6)


def test_the_schemes_reduce_as_they_should():
    net, base, loop = _loop()
    W = np.array([0.40, 0.25, 0.35]); T0 = E.exact_rep(net, W)
    run = lambda **kw: loop.run_dithered(T0, 3 * W, 2.5, 0.25, 0.0, np.random.default_rng(5), **kw)
    T_ito, b_ito = run(label_v=1e-4, label_scheme='ito')
    T_ar, b_ar = run(label_v=1e-4, label_scheme='ar1', label_rho=0.0)
    assert np.array_equal(T_ito, T_ar) and np.array_equal(b_ito, b_ar)
    T_h, b_h = run(label_v=1e-30, label_scheme='heun')                     # no noise to speak of: the explicit trapezoidal rule
    Th, lb = T0.copy(), np.log(3 * W)
    for _ in range(10):
        d1, g1 = loop.drift(Th, np.exp(lb))
        Tp, lp = Th + 0.25 * d1, lb + 0.25 * g1 / np.exp(lb)
        d2, g2 = loop.drift(Tp, np.exp(lp))
        Th, lb = Th + 0.125 * (d1 + d2), lb + 0.125 * (g1 / np.exp(lb) + g2 / np.exp(lp))
        lb = lb - np.log(np.exp(lb).sum() / 3)
    assert np.allclose(T_h, Th, rtol=0, atol=1e-12) and np.allclose(b_h, np.exp(lb), rtol=0, atol=1e-12)
    assert LN.ar1_factor(0.0) == 1.0 and abs(LN.ar1_factor(1 / 3) - 0.5) < 1e-15 and abs(LN.ar1_factor(-1 / 3) - 2.0) < 1e-15


def test_the_scalar_factorisation_example():
    # f = uv, label noise: an Euler step multiplies Q = u² − v² by exactly 1 − (ξ − r)²h² for any draw; Heun's step keeps Q to O(h³)
    rng = np.random.default_rng(6)
    u, v = 1.0 + rng.random(50), 0.3 + rng.random(50)
    xi, h = 3.0 * rng.standard_normal(50), 1e-2
    c = xi - (u * v - 1)
    un, vn = u + h * c * v, v + h * c * u
    assert np.allclose(un ** 2 - vn ** 2, (u ** 2 - v ** 2) * (1 - (c * h) ** 2), rtol=1e-12, atol=1e-14)
    up, vp = u + h * c * v, v + h * c * u
    cp = xi - (up * vp - 1)
    uh, vh = u + 0.5 * h * (c * v + cp * vp), v + 0.5 * h * (c * u + cp * up)
    assert np.abs((uh ** 2 - vh ** 2) - (u ** 2 - v ** 2)).max() < 50 * h ** 3 * np.abs(c).max() ** 3


def test_the_diagonal_network_example():
    # w = u⊙u − v⊙v with label noise: an Euler step multiplies each P_i = u_iv_i by exactly 1 − 4h²(g_i − ξ̃_i)², g = ∇_wL with the noisy labels
    import exp_gauge_label as GL
    X, y, ws = GL.diag_problem()
    rng = np.random.default_rng(8)
    u, v = 0.3 + rng.random(GL.D_E), 0.3 + rng.random(GL.D_E)
    xi, h, n = rng.standard_normal(GL.N_E), 1e-2, GL.N_E
    g = ((u * u - v * v) @ X.T - y - xi) @ X / n
    un, vn = u - h * 2 * g * u, v + h * 2 * g * v
    assert np.allclose(un * vn, u * v * (1 - 4 * h * h * g * g), rtol=1e-12)


class _Draws:
    """An rng that records the label noise's draws (shape (N,)) and returns zeros for any other shape (no dither)."""
    def __init__(self, N, seed):
        self.N, self.rng, self.out = N, np.random.default_rng(seed), []

    def standard_normal(self, shape):
        if tuple(np.atleast_1d(shape)) == (self.N,):
            self.out.append(self.rng.standard_normal(self.N)); return self.out[-1]
        return np.zeros(shape)


def _recording(loop):
    seen, drift = [], loop.drift

    def rec(Th, b, xi=None):
        seen.append(None if xi is None else np.array(xi)); return drift(Th, b, xi)
    loop.drift = rec
    return seen


def test_the_ar1_draws():
    net, base, loop = _loop()
    W = np.array([0.40, 0.25, 0.35])
    v, dt, rho = 3e-4, 0.25, 0.6
    seen, rng = _recording(loop), _Draws(base.N, 3)
    loop.run_dithered(E.exact_rep(net, W), 3 * W, 5 * dt, dt, 0.0, rng, label_v=v, label_scheme='ar1', label_rho=rho)
    s = np.sqrt(v * base.N / dt * (1 - rho) / (1 + rho))                 # the per-step std at the intensity v, started stationary
    assert len(seen) == 5 and len(rng.out) == 5
    xi = s * rng.out[0]
    for k in range(5):
        if k:
            xi = rho * xi + np.sqrt(1 - rho ** 2) * s * rng.out[k]
        assert np.allclose(seen[k], xi, rtol=1e-13, atol=0)


def test_heun_with_correlated_noise_off_the_manifold():
    net, base, loop = _loop(gamma=0.3)
    W = np.array([0.40, 0.25, 0.35]); b0 = 3 * W
    T0 = E.exact_rep(net, W) + 0.05 * np.random.default_rng(4).standard_normal((3, 3))
    v, dt, rho = 1e-3, 0.25, 0.5
    rng = _Draws(base.N, 9)
    T1, b1 = loop.run_dithered(T0, b0, 2 * dt, dt, 0.0, rng, label_v=v, label_scheme='heun', label_rho=rho)
    assert len(rng.out) == 2                                              # one draw per step: the corrector reuses the predictor's
    s = np.sqrt(v * base.N / dt * (1 - rho) / (1 + rho))
    xi = s * rng.out[0]
    Th, lb = T0.copy(), np.log(b0)
    for k in range(2):
        if k:
            xi = rho * xi + np.sqrt(1 - rho ** 2) * s * rng.out[1]
        d1, g1 = loop.drift(Th, np.exp(lb), xi)
        Tp, lp = Th + dt * d1, lb + dt * g1 / np.exp(lb)
        lp = lp - np.log(np.exp(lp).sum() / 3)                            # the predictor renormalised to Σβ = m
        d2, g2 = loop.drift(Tp, np.exp(lp), xi)
        Th, lb = Th + 0.5 * dt * (d1 + d2), lb + 0.5 * dt * (g1 / np.exp(lb) + g2 / np.exp(lp))
        lb = lb - np.log(np.exp(lb).sum() / 3)
    assert np.allclose(T1, Th, rtol=0, atol=1e-13) and np.allclose(b1, np.exp(lb), rtol=0, atol=1e-13)
    assert np.abs(T1 - T0).max() > 1e-3                                   # and it moved


def test_the_weighted_and_the_covariance_wong_zakai_terms():
    import exp_gauge_drift as X
    import exp_gauge_label as GL
    F, chart, coords, S0 = X.instance(1e-2)
    net, base, loop = _loop()
    x0 = chart(np.array([0.40, 0.25]))
    Dd = LN.field_derivatives(loop, x0, GL.beta_of)                        # D[j, i] = ∂_{U_j}U_i
    N = base.N
    wt = 1 + 0.8 * np.sin(2 * np.pi * np.arange(N) / N)
    direct = np.einsum('x,xxk->k', wt, Dd) / N
    assert np.allclose(LN.wong_zakai(loop, x0, GL.beta_of, weights=wt), direct, rtol=1e-10, atol=0)
    assert np.allclose(LN.wong_zakai(loop, x0, GL.beta_of, cov=np.diag(wt)), direct, rtol=1e-12, atol=0)
    assert np.allclose(LN.wong_zakai(loop, x0, GL.beta_of), np.einsum('xxk->k', Dd) / N, rtol=1e-10, atol=0)
    C = np.exp(-np.abs(np.subtract.outer(np.arange(N), np.arange(N))) / 5.0)
    assert np.allclose(LN.wong_zakai(loop, x0, GL.beta_of, cov=C), np.einsum('ij,jik->k', C, Dd) / N, rtol=1e-12, atol=0)


def test_the_linear_response_factor():
    from scipy.linalg import solve_discrete_lyapunov
    for lam, h, rho in ((3.8, 0.01, 0.9), (0.5, 0.25, 0.6), (2.0, 0.1, -0.4), (1.0, 0.05, 0.0), (0.2, 0.25, 0.999)):
        a = 1 - h * lam                                                   # (r, ξ): ξ' = ρξ + √(1 − ρ²)z, r' = ar + hλξ; e = r − ξ
        S = solve_discrete_lyapunov(np.array([[a, h * lam], [0.0, rho]]), np.diag([0.0, 1 - rho ** 2]))
        e2 = S[0, 0] - 2 * S[0, 1] + S[1, 1]
        assert abs(LN.lr_factor(lam, h, ('ar1', rho) if rho else ('fresh',)) - e2) < 1e-10 * e2
    for lam, h, kb in ((1.3, 0.1, 4), (6.0, 0.05, 3), (0.4, 0.25, 7)):      # a draw kept kb steps, averaged over the phase: the
        a, S = 1 - h * lam, np.zeros((2, 2))                              # periodic covariance of (r, ξ), E[e²] over a period
        for it in range(400 * kb):
            if it % kb == 0:
                S[0, 1] = S[1, 0] = 0.0; S[1, 1] = 1.0
            if it >= 399 * kb:
                e2 = (e2 if it > 399 * kb else 0.0) + (S[0, 0] - 2 * S[0, 1] + S[1, 1]) / kb
            S = np.array([[a * a * S[0, 0] + 2 * a * h * lam * S[0, 1] + (h * lam) ** 2 * S[1, 1], a * S[0, 1] + h * lam * S[1, 1]],
                          [a * S[0, 1] + h * lam * S[1, 1], S[1, 1]]])
        assert abs(LN.lr_factor(lam, h, ('block', kb)) - e2) < 1e-10 * e2
    assert abs(LN.lr_factor(1.3, 0.1, ('fresh',)) - (1 + 0.13 / (2 - 0.13))) < 1e-14
    assert abs(LN.lr_factor(1.3, 0.1, ('block', 1)) - LN.lr_factor(1.3, 0.1, ('fresh',))) < 1e-14
    for corr in (('fresh',), ('ar1', 0.9), ('block', 4), ('gauss', 3.0)):  # a mode that is not restored passes the noise unchanged,
        assert np.allclose(LN.lr_factor(np.array([0.0, 1e-18, 1e-12]), 0.5, corr), 1.0, atol=1e-10)   # also when 1 − hλ rounds to 1
    tc, h = 2.0, 1e-5                                                     # finely sampled AR(1): 1/(1 + λτ_c), Ψ = log(1 + λτ_c)/τ_c,
    for lam in (0.1, 1.5, 20.0):                                          # up to O(hλ)
        assert abs(LN.lr_factor(lam, h, ('ar1', np.exp(-h / tc))) * (1 + lam * tc) - 1) < 2 * h * lam + 1e-6
        assert abs(LN.lr_potential(lam, h, ('ar1', np.exp(-h / tc)))[0] * tc / np.log(1 + lam * tc) - 1) < 2 * h * lam + 1e-6


def test_the_filtered_covariance_weights_the_pairs_of_modes():
    rng = np.random.default_rng(1)
    n, k, h = 7, 2, 0.05
    V = np.linalg.qr(rng.standard_normal((n, n)))[0] + 0.3 * rng.standard_normal((n, n))   # a non-orthogonal eigenbasis
    lam = np.array([0.0, 0.0, 0.3, 1.0, 2.5, 2.5, 10.0])                  # a kernel of dimension k, a double rate, hλ_max = 0.5
    Vi = np.linalg.inv(V); A = V @ np.diag(-lam) @ Vi; P = V[:, :k] @ Vi[:k, :]
    Ub = (np.eye(n) - P) @ rng.standard_normal((n, 5)); Sig = Ub @ Ub.T
    for corr in (('fresh',), ('ar1', 0.6), ('ar1', -0.4), ('ar1', 0.95), ('block', 4)):
        Cf = Vi @ LN.filtered_covariance(A, P, Sig, h, corr) @ Vi.T
        C = Vi @ Sig @ Vi.T
        a = 1 - h * lam[k:]; G, R = LN.corr_sums(a, corr); L = lam[k:]
        f = (1 - h * (L[:, None] * R[:, None] + L[None, :] * R[None, :])
             + h * h * L[:, None] * L[None, :] * (1 + G[:, None] + G[None, :]) / (1 - a[:, None] * a[None, :]))
        assert np.abs(Cf[k:, k:] - f * C[k:, k:]).max() < 1e-12 * np.abs(C).max()
        assert np.abs(Cf[:k]).max() < 1e-12 * np.abs(C).max()


def test_smooth_correlations_and_the_slow_limit_at_a_finite_step():
    from scipy.linalg import solve_discrete_lyapunov
    k = 30.0                                                              # Gaussian step correlations exp(−½(d/k)²): the sums,
    d = np.arange(1, 3000); rd = np.exp(-0.5 * (d / k) ** 2)              # scalar and matrix (through an eigendecomposition)
    for a in (0.95, -0.4, 0.999):
        G, R = LN.corr_sums(a, ('gauss', k))
        assert abs(G - (rd * a ** d).sum()) < 1e-12 and abs(R - (rd * a ** (d - 1)).sum()) < 1e-12
    Q = np.linalg.qr(np.random.default_rng(2).standard_normal((4, 4)))[0]
    B = Q @ np.diag([0.9, 0.5, -0.3, 0.99]) @ Q.T
    Gm = LN.corr_sums(B, ('gauss', k))[0]
    assert np.abs(Gm - Q @ np.diag([LN.corr_sums(x, ('gauss', k))[0] for x in (0.9, 0.5, -0.3, 0.99)]) @ Q.T).max() < 1e-12
    h = 0.25                                                              # a smooth slow noise: f ≈ (λτ_c)⁻², an AR(1) one 1/(1 + λτ_c)
    for lam, tc in ((1.0, 400.0), (0.1, 4000.0)):
        assert abs(LN.lr_factor(lam, h, ('gauss', tc / h)) * (lam * tc) ** 2 - 1) < 2e-2
    rng = np.random.default_rng(1)                                        # AR(1) at τ_c → ∞ and a finite step: τ_cΣ_f → X_h,
    n, kk, tc = 6, 2, 1e8                                                 # X_h − BX_hBᵀ = 2hS, the pair (i, l) weighted by
    V = np.linalg.qr(rng.standard_normal((n, n)))[0] + 0.3 * rng.standard_normal((n, n))   # 2/(λ_i + λ_l − hλ_iλ_l)
    lam = np.array([0.0, 0.0, 0.3, 1.0, 2.5, 6.0]); Vi = np.linalg.inv(V)
    A = V @ np.diag(-lam) @ Vi; P = V[:, :kk] @ Vi[:kk, :]
    Ub = (np.eye(n) - P) @ rng.standard_normal((n, 4)); Sig = Ub @ Ub.T
    Cf = Vi @ (tc * LN.filtered_covariance(A, P, Sig, h, ('ar1', float(np.exp(-h / tc))))) @ Vi.T
    C = Vi @ Sig @ Vi.T; L = lam[kk:]
    pair = 2 / (L[:, None] + L[None, :] - h * L[:, None] * L[None, :])
    assert np.abs(Cf[kk:, kk:] - pair * C[kk:, kk:]).max() < 1e-5 * np.abs(pair * C[kk:, kk:]).max()
    Bn = np.diag(1 - h * L)
    assert np.abs(solve_discrete_lyapunov(Bn, 2 * h * C[kk:, kk:]) - pair * C[kk:, kk:]).max() < 1e-12 * np.abs(C).max()


def test_the_amplitude_coordinates_target():
    from scipy.optimize import minimize
    net, base, loop = _loop()
    Wm = np.abs(E.A_BUMP) / net.a_max
    P = E.G[0] * np.mean(base.features(E.exact_rep(net, np.ones(3) / 3))[4] ** 2, axis=0)
    Wa = LN.w_label(P, Wm, coord='a')
    c = (3 * np.abs(E.A_BUMP) / net.a_max) ** 2                            # the strictly convex potential, minimised directly
    pot = lambda W: float(np.sum(P * ((3 * W) ** 2 / 2 + 2 * c * np.log(3 * W) + 3 * c ** 2 / (2 * (3 * W) ** 2))))
    res = minimize(lambda w: pot(np.r_[w, 1 - w.sum()]), np.array([0.3, 0.45]), method='Nelder-Mead',
                   options=dict(xatol=1e-12, fatol=1e-16, maxiter=20000))
    assert np.abs(np.r_[res.x, 1 - res.x.sum()] - Wa).max() < 1e-6
    assert np.abs(Wa - LN.w_label(P, Wm)).max() > 1e-2                    # not the α-coordinates' target
    u = LN.selection(loop, E.exact_rep(net, Wa), 3 * Wa, coord='a')       # the closed form's replicator flow is stationary there
    assert np.abs(u - u.mean()).max() < 1e-9 * np.abs(u).max()
    u0 = LN.selection(loop, E.exact_rep(net, Wa), 3 * Wa)
    assert np.abs(u0 - u0.mean()).max() > 1e-3 * np.abs(u0).max()


def test_heun_with_dither_and_label_noise():
    net, base, loop = _loop(gamma=0.3)
    W = np.array([0.40, 0.25, 0.35]); b0 = 3 * W
    T0 = E.exact_rep(net, W) + 0.02 * np.random.default_rng(5).standard_normal((3, 3))
    v, dt, tau = 1e-3, 0.25, 1e-3

    class Both(_Draws):
        def standard_normal(self, shape):
            z = self.rng.standard_normal(shape); self.out.append(z); return z
    rng = Both(base.N, 11)
    T1, b1 = loop.run_dithered(T0, b0, dt, dt, tau, rng, label_v=v, label_scheme='heun')
    xi = np.sqrt(v * base.N / dt) * rng.out[0]; dW = np.sqrt(2 * tau * dt * loop.G)[None, :] * rng.out[1]
    d1, g1 = loop.drift(T0, b0, xi)
    Tp, lp = T0 + dt * d1 + dW, np.log(b0) + dt * g1 / b0                 # the dither in the predictor and in the step, once each
    lp = lp - np.log(np.exp(lp).sum() / 3)
    d2, g2 = loop.drift(Tp, np.exp(lp), xi)
    Th, lb = T0 + 0.5 * dt * (d1 + d2) + dW, np.log(b0) + 0.5 * dt * (g1 / b0 + g2 / np.exp(lp))
    lb = lb - np.log(np.exp(lb).sum() / 3)
    assert np.allclose(T1, Th, rtol=0, atol=1e-13) and np.allclose(b1, np.exp(lb), rtol=0, atol=1e-13)


def test_the_zero_search_resolves_zeros_near_a_saturation():
    """Part G's flow and zero search on a field with a stable zero 6·10⁻⁵ above a saturation (W₃ = 0.20006, W_min,3 = 0.2): the flow's
    steps are capped by the distance to the guard, the zero accepted relative to the field's scale (referee round 31)."""
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
    import exp_gauge_label as GL
    Wm = np.array([0.25, 0.175, 0.2])
    target = np.array([0.6249, 0.17504])
    orig = GL.g_dr
    GL.g_dr = lambda gamma, w, corr, Wm_: (np.full(2, np.nan) if np.any(np.r_[w, 1 - np.sum(w)] <= GL.GUARD_G * Wm_)
                                         else -(np.asarray(w, float) - target))
    try:
        w, why, lam = GL.g_flow(0.0, np.array([0.40, 0.25]), None, Wm, ref=1.0)
        assert why == 'zero' and np.abs(w - target).max() < 1e-9 and lam.max() < 0
        z = GL.g_zero(0.0, target + 1e-7, None, Wm, ref=1.0)
        assert z is not None and np.abs(z[0] - target).max() < 1e-10 and np.allclose(z[1], -1.0)
    finally:
        GL.g_dr = orig
