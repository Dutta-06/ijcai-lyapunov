"""
Projection operator tests.

The only property of the projection that the Lyapunov proof uses is

    (P)   θ̃ᵀ Γ⁻¹ Proj(θ, Γy) ≤ θ̃ᵀ y     for all θ ∈ Ω, all y,  θ̃ = θ − θ*,

valid whenever θ* lies in the inner set.  These tests check (P) directly
(randomised, isotropic and non-isotropic Γ, vector / matrix / box), check that
the boundary is never crossed in continuous time, and document that (P) can
fail when θ* is outside the inner set (which is why Assumption A3 is needed).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.projection import ProjectionOperator, BandwidthProjection, CenterProjection


def _spd(rng, d, cond=20.0):
    Q, _ = np.linalg.qr(rng.normal(size=(d, d)))
    lam = np.exp(rng.uniform(0, np.log(cond), size=d))
    return Q @ np.diag(lam) @ Q.T


def test_norm_ball_projection_basic():
    proj = ProjectionOperator(theta_max=10.0, epsilon=1.0)
    theta = np.array([5.0, 0.0]); y = np.array([1.0, 1.0])
    assert np.allclose(proj.project_vector(theta, y), y)
    theta = np.array([10.0, 0.0]); y = np.array([2.0, 1.0])
    out = proj.project_vector(theta, y)
    assert np.isclose(out[0], 0.0) and np.isclose(out[1], 1.0)
    y_in = np.array([-2.0, 1.0])
    assert np.allclose(proj.project_vector(theta, y_in), y_in)


def test_property_P_vector_isotropic_and_spd():
    rng = np.random.default_rng(0)
    d, theta_max, eps = 6, 3.0, 0.5
    proj = ProjectionOperator(theta_max, eps)
    worst = -np.inf
    for trial in range(3000):
        # θ* in the inner set, θ anywhere in the ball (biased towards the layer)
        theta_star = rng.normal(size=d); theta_star *= rng.uniform(0, theta_max - eps) / np.linalg.norm(theta_star)
        theta = rng.normal(size=d); theta *= rng.uniform(theta_max - 1.5 * eps, theta_max) / np.linalg.norm(theta)
        grad = rng.normal(size=d) * rng.uniform(0.1, 10)
        if trial % 2 == 0:
            gamma_val = rng.uniform(0.1, 10)
            G = gamma_val * np.eye(d)
            out = proj.project_vector(theta, G @ grad)            # isotropic fast path
        else:
            G = _spd(rng, d)
            out = proj.project_vector(theta, G @ grad, G)         # Γ-weighted path
        lhs = (theta - theta_star) @ np.linalg.solve(G, out)
        rhs = (theta - theta_star) @ grad
        worst = max(worst, lhs - rhs)
    assert worst <= 1e-10, f"property (P) violated by {worst}"


def test_property_P_matrix():
    rng = np.random.default_rng(1)
    m, n, theta_max, eps = 4, 3, 5.0, 0.7
    proj = ProjectionOperator(theta_max, eps)
    worst = -np.inf
    for trial in range(1500):
        Ws = rng.normal(size=(m, n)); Ws *= rng.uniform(0, theta_max - eps) / np.linalg.norm(Ws, 'fro')
        W = rng.normal(size=(m, n)); W *= rng.uniform(theta_max - 1.5 * eps, theta_max) / np.linalg.norm(W, 'fro')
        grad = rng.normal(size=(m, n))
        if trial % 2 == 0:
            g = rng.uniform(0.1, 10); G = g * np.eye(m * n)
            out = proj.project_matrix(W, g * grad)
        else:
            G = _spd(rng, m * n)
            out = proj.project_matrix(W, (G @ grad.reshape(-1)).reshape(m, n), G)
        lhs = (W - Ws).reshape(-1) @ np.linalg.solve(G, out.reshape(-1))
        rhs = np.sum((W - Ws) * grad)
        worst = max(worst, lhs - rhs)
    assert worst <= 1e-10, f"matrix property (P) violated by {worst}"


def test_property_P_box_bandwidths():
    rng = np.random.default_rng(2)
    s_min, s_max, eps = 0.3, 3.0, 0.2
    proj = BandwidthProjection(s_min, s_max, eps)
    worst = -np.inf
    for _ in range(3000):
        m = 5
        s_star = rng.uniform(s_min + eps, s_max - eps, size=m)
        s = rng.uniform(s_min, s_max, size=m)
        grad = rng.normal(size=m) * rng.uniform(0.1, 10)
        gdiag = rng.uniform(0.1, 10, size=m)
        out = proj.project_vector(s, gdiag * grad)
        worst = max(worst, (s - s_star) @ (out / gdiag) - (s - s_star) @ grad)
    assert worst <= 1e-10


def test_boundary_never_crossed_in_continuous_time():
    """On the boundary ‖θ‖ = θ_max the projected update has no outward component."""
    rng = np.random.default_rng(3)
    proj = ProjectionOperator(2.0, 0.3)
    for _ in range(500):
        theta = rng.normal(size=4); theta *= 2.0 / np.linalg.norm(theta)
        y = rng.normal(size=4) * 5
        out = proj.project_vector(theta, y)
        assert theta @ out <= 1e-12
    bp = BandwidthProjection(0.5, 4.0, 0.1)
    s = np.array([0.5, 4.0]); y = np.array([-3.0, 3.0])
    assert np.all(bp.project_vector(s, y) == 0.0)


def test_property_P_can_fail_when_ideal_outside_inner_set():
    """Documents Assumption A3: with θ* outside the inner set, (P) is not guaranteed."""
    proj = ProjectionOperator(theta_max=1.0, epsilon=0.5)
    theta = np.array([1.0, 0.0])          # on the boundary
    theta_star = np.array([2.0, 0.0])     # outside the set
    y = np.array([1.0, 0.0])              # outward
    out = proj.project_vector(theta, y)   # = 0
    assert (theta - theta_star) @ out > (theta - theta_star) @ y  # (P) violated


def test_center_projection_default_is_pure_norm_ball():
    rng = np.random.default_rng(4)
    cp = CenterProjection(c_max=3.0, delta_min=0.5, epsilon=0.2)
    norm = ProjectionOperator(3.0, 0.2)
    for _ in range(200):
        c = rng.normal(size=6); y = rng.normal(size=6)
        assert np.allclose(cp.project_vector(c, y, 2), norm.project_vector(c, y))


def test_legacy_separation_breaks_property_P():
    """The legacy repulsion is outside the proof: exhibit a violation of (P)."""
    # The repulsion changes θ̃ᵀ(update) by  corr · (dist − (cᵢ* − cⱼ*)·dir),
    # which is positive whenever the ideal centres are less separated ALONG THE
    # CURRENT direction than the current ones — e.g. ideal separation orthogonal
    # to it.  Such θ* is admissible (well separated, inside the ball).
    cp = CenterProjection(c_max=10.0, delta_min=0.5, epsilon=0.2, enforce_separation=True)
    c = np.array([[0.0, 0.0], [0.55, 0.0]]).reshape(-1)          # within δ+ε of each other
    c_star = np.array([[0.0, 0.0], [0.0, 0.8]]).reshape(-1)      # ideal: separated by 0.8, orthogonally
    y = np.array([[0.0, 0.0], [-1.0, 0.0]]).reshape(-1)          # compressive raw update
    out = cp.project_vector(c, y, 2)
    assert (c - c_star) @ out > (c - c_star) @ y + 1e-9
    # and without the legacy option the same inputs satisfy (P)
    cp0 = CenterProjection(c_max=10.0, delta_min=0.5, epsilon=0.2)
    assert (c - c_star) @ cp0.project_vector(c, y, 2) <= (c - c_star) @ y + 1e-12


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "passed")
