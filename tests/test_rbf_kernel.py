import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.rbf_kernel import RBFKernel


def test_rbf_kernel_values():
    kernel = RBFKernel(m=2, n=2)
    phi = kernel.phi(np.array([0.0, 0.0]), np.array([[0.0, 0.0], [1.0, 0.0]]), np.array([1.0, 1.0]))
    assert np.isclose(phi[0], 1.0) and np.isclose(phi[1], np.exp(-0.5))


def test_jacobians_against_finite_differences():
    rng = np.random.default_rng(0)
    m, n = 5, 3
    kernel = RBFKernel(m, n)
    x = rng.normal(size=n); c = rng.normal(size=(m, n)); s = rng.uniform(0.5, 2.0, size=m)
    jc = kernel.phi_jacobian_c(x, c, s)
    js = kernel.phi_jacobian_sigma(x, c, s)
    h = 1e-6
    for i in range(m):
        for k in range(n):
            cp, cm = c.copy(), c.copy(); cp[i, k] += h; cm[i, k] -= h
            num = (kernel.phi(x, cp, s) - kernel.phi(x, cm, s)) / (2 * h)
            assert np.allclose(jc[:, i * n + k], num, atol=1e-6)
        sp, sm = s.copy(), s.copy(); sp[i] += h; sm[i] -= h
        num = (kernel.phi(x, c, sp) - kernel.phi(x, c, sm)) / (2 * h)
        assert np.allclose(js[:, i], num, atol=1e-6)


def test_block_jacobian_matches_dense():
    rng = np.random.default_rng(1)
    m, n = 6, 2
    kernel = RBFKernel(m, n)
    x = rng.normal(size=n); c = rng.normal(size=(m, n)); s = rng.uniform(0.5, 2.0, size=m)
    blocks = kernel.phi_jacobian_c_blocks(x, c, s)
    dense = kernel.phi_jacobian_c(x, c, s)
    for i in range(m):
        assert np.allclose(dense[i, i * n:(i + 1) * n], blocks[i])
        assert np.allclose(np.delete(dense[i], np.arange(i * n, (i + 1) * n)), 0.0)


def test_jacobian_sup_norms_are_true_suprema():
    """Analytic sup of ‖φ'_c‖₂, ‖φ'_σ‖₂ over σ ≥ σ_min vs brute-force maximisation."""
    s_min = 0.4
    kernel = RBFKernel(1, 1, sigma_min=s_min)
    J_c, J_s = kernel.jacobian_sup_norms()
    r = np.linspace(0, 10, 200001)
    best_c = best_s = 0.0
    for s in [s_min, 0.5, 0.8, 1.3, 2.0]:
        phi = np.exp(-r ** 2 / (2 * s ** 2))
        best_c = max(best_c, np.max(phi * r / s ** 2))
        best_s = max(best_s, np.max(phi * r ** 2 / s ** 3))
    assert best_c <= J_c + 1e-9 and np.isclose(best_c, J_c, rtol=1e-6)
    assert best_s <= J_s + 1e-9 and np.isclose(best_s, J_s, rtol=1e-6)
    # spectral norm of the dense Jacobian is bounded by the row-block sup
    rng = np.random.default_rng(2)
    k2 = RBFKernel(7, 3, sigma_min=s_min)
    for _ in range(200):
        x = rng.normal(size=3) * 3; c = rng.normal(size=(7, 3)) * 3; s = rng.uniform(s_min, 3.0, size=7)
        assert np.linalg.norm(k2.phi_jacobian_c(x, c, s), 2) <= J_c + 1e-12
        assert np.linalg.norm(k2.phi_jacobian_sigma(x, c, s), 2) <= J_s + 1e-12


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "passed")
