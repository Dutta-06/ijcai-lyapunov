"""The fast (scalar Γ) and dense (matrix Γ) paths must compute the same law."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.controller import AdaptiveController


def _rand_state(rng, m, n, ctrl):
    x = rng.normal(size=n)
    e = rng.normal(size=n)
    W = rng.normal(size=(m, n)); W *= rng.uniform(0.2, 1.0) * ctrl.W_max / np.linalg.norm(W, 'fro')
    c = rng.normal(size=(m, n)); c *= rng.uniform(0.2, 1.0) * ctrl.c_max / np.linalg.norm(c)
    s = rng.uniform(ctrl.sigma_min, ctrl.sigma_max, size=m)
    return x, e, W, c, s


def test_scalar_and_matrix_gamma_agree():
    rng = np.random.default_rng(0)
    m, n = 6, 3
    kw = dict(m=m, n=n, ke=4.0, W_max=5.0, c_max=4.0, sigma_min=0.2, sigma_max=3.0, proj_epsilon=0.3)
    fast = AdaptiveController(gamma_W=2.0, gamma_c=0.7, gamma_sigma=0.4, **kw)
    dense = AdaptiveController(gamma_W=2.0 * np.eye(m * n), gamma_c=0.7 * np.eye(m * n),
                               gamma_sigma=0.4 * np.eye(m), **kw)
    for _ in range(300):
        x, e, W, c, s = _rand_state(rng, m, n, fast)
        a = fast.adaptation_laws(x, e, W, c, s)
        b = dense.adaptation_laws(x, e, W, c, s)
        for u, v in zip(a, b):
            assert np.allclose(u, v, atol=1e-12, rtol=1e-12)


def test_dense_path_uses_gamma_weighted_projection():
    """With non-isotropic Γ the projected direction is Γθ (checked via property P)."""
    rng = np.random.default_rng(1)
    m, n = 4, 2
    Q, _ = np.linalg.qr(rng.normal(size=(m * n, m * n)))
    GW = Q @ np.diag(rng.uniform(0.5, 8.0, size=m * n)) @ Q.T
    ctrl = AdaptiveController(m=m, n=n, ke=3.0, gamma_W=GW, gamma_c=1.0, gamma_sigma=0.5,
                              W_max=3.0, c_max=4.0, sigma_min=0.2, sigma_max=3.0, proj_epsilon=0.4)
    worst = -np.inf
    for _ in range(500):
        x, e, W, c, s = _rand_state(rng, m, n, ctrl)
        W *= (ctrl.W_max - rng.uniform(0, 0.3)) / np.linalg.norm(W, 'fro')   # in the layer
        Ws = rng.normal(size=(m, n)); Ws *= rng.uniform(0, ctrl.W_max - ctrl.proj_epsilon) / np.linalg.norm(Ws, 'fro')
        W_dot, _, _ = ctrl.adaptation_laws(x, e, W, c, s)
        phi = ctrl.kernel.phi(x, c, s)
        lhs = (W - Ws).reshape(-1) @ ctrl.gamma_W_inv @ W_dot.reshape(-1)
        rhs = np.sum((W - Ws) * np.outer(phi, e))
        worst = max(worst, lhs - rhs)
    assert worst <= 1e-10


def test_zero_learning_rate_freezes_group():
    rng = np.random.default_rng(2)
    m, n = 5, 2
    ctrl = AdaptiveController(m=m, n=n, gamma_W=3.0, gamma_c=0.0, gamma_sigma=0.0,
                              W_max=5.0, c_max=5.0, sigma_min=0.2, sigma_max=3.0)
    x, e, W, c, s = _rand_state(rng, m, n, ctrl)
    W_dot, c_dot, s_dot = ctrl.adaptation_laws(x, e, W, c, s)
    assert np.all(c_dot == 0.0) and np.all(s_dot == 0.0) and np.any(W_dot != 0.0)


def test_rejects_non_diagonal_gamma_sigma_and_indefinite_gamma():
    m, n = 3, 2
    bad = np.ones((m, m))
    try:
        AdaptiveController(m=m, n=n, gamma_sigma=bad)
        assert False, "should reject non-diagonal gamma_sigma"
    except ValueError:
        pass
    try:
        AdaptiveController(m=m, n=n, gamma_W=-np.eye(m * n))
        assert False, "should reject indefinite gamma_W"
    except ValueError:
        pass


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "passed")
