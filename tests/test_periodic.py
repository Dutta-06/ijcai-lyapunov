"""
Checks for PARADIGM.md Theorem J (periodic orbits of the finite-m closed loop near the critical points of J_m), using
src/core/periodic.py on the n = 1, m = 2 two-bump instance of experiments/exp_periodic_orbit.py.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core import langevin as L
from src.core import periodic as PO

PER = 2 * np.pi


def _instance():
    net = L.SaturatedLiftedNetwork(1, a_max=3.0, c_max=2.5, s_min=0.2, s_max=1.5)
    f = lambda X: 1.6 * np.exp(-(X - 1.0) ** 2 / 0.4) - 1.2 * np.exp(-(X + 1.0) ** 2 / 0.4) + 0.1 * X

    def J(X):
        d = (1.6 * np.exp(-(X - 1.0) ** 2 / 0.4) * (-2 * (X - 1.0) / 0.4) - 1.2 * np.exp(-(X + 1.0) ** 2 / 0.4) * (-2 * (X + 1.0) / 0.4) + 0.1)
        return d[:, :, None]
    xr = lambda t: np.array([2.0 * np.sin(t)]); xd = lambda t: np.array([2.0 * np.cos(t)])
    Th0 = net.unconstrained(np.array([[1.0], [-1.0]]), np.array([[0.8], [-0.8]]), np.array([0.6, 0.6]))
    return net, f, J, xr, xd, np.array([2.0, 1.0, 1.0]), 0.05, Th0


def test_gradient_of_J_m_matches_finite_differences_and_critical_point_is_nondegenerate():
    net, f, J, xr, xd, G, s0, Th0 = _instance()
    cl = PO.ClosedLoop(net, f, xr, xd, PER, 20.0, G, s0)
    orb = cl.orbit(N=1000)
    rng = np.random.default_rng(0)
    Th = Th0 + 0.1 * rng.normal(size=Th0.shape)
    g = cl.grad_J_m(Th, orb)
    h = 1e-6
    for i in range(Th.shape[0]):
        for k in range(Th.shape[1]):
            E = np.zeros_like(Th); E[i, k] = h
            fd = (cl.J_m(Th + E, orb) - cl.J_m(Th - E, orb)) / (2 * h)
            assert abs(fd - g[i, k]) < 1e-7 * (1 + abs(fd)), (i, k, fd, g[i, k])
    Th_star, H = cl.find_critical_point(Th0, orb)
    assert np.linalg.norm(cl.grad_J_m(Th_star, orb)) < 1e-10
    assert np.linalg.eigvalsh(H)[0] > 1e-3                                   # non-degenerate minimum
    a, c, s = net.physical(Th_star)
    assert abs(c[0, 0] - c[1, 0]) > 1.0 and a[0, 0] * a[1, 0] < 0            # two distinct kernels of opposite sign


def test_theorem_J_periodic_orbit_offset_and_floquet_exponents():
    """
    For k_e = 20 and 40: Newton shooting finds the periodic orbit of the closed loop (residual < 1e-10) near Θ*; the
    period-mean offset Θ̄_p − Θ* has size O(1/k_e) and is predicted by Θ* − J⁻¹(b_frozen + b_osc) with a vector error that
    halves when k_e doubles (the remainder is O(k_e⁻²) against an O(1/k_e) offset); the parameter Floquet exponents match
    the averaged prediction (1/k_e)λ(G_m^{1/2}∇²J_mG_m^{1/2}) with a maximal relative error that halves as well.
    """
    net, f, J, xr, xd, G, s0, Th0 = _instance()
    m = Th0.shape[0]
    w = 1 / np.tile(G, m)
    nrm = lambda v: float(np.sqrt(np.sum(w * np.ravel(v) ** 2)))
    out = {}
    for ke in (20.0, 40.0):
        cl = PO.ClosedLoop(net, f, xr, xd, PER, ke, G, s0)
        orb = cl.orbit(N=1500)
        Th_star, H = cl.find_critical_point(Th0, orb)
        z0 = np.concatenate([xr(0.0), Th_star.ravel()])
        zp, DP, res = cl.shoot(z0, m)
        assert res < 1e-10, res
        Thbar, Th_t, E_t, _ = cl.period_mean(zp, m)
        full, frozen = cl.offset_prediction(Th_star, orb, J)
        off = Thbar - Th_star
        fast, slow = PO.floquet_exponents(DP, PER, 1)
        pred = PO.predicted_parameter_exponents(H, G, m, ke)
        out[ke] = dict(off=nrm(off), verr=nrm(off - full) / nrm(full), rel=np.max(np.abs(slow - pred) / pred),
                       emax=float(np.abs(E_t).max()))
    assert out[40.0]['off'] < 0.65 * out[20.0]['off']                        # ‖Θ̄_p − Θ*‖ = O(1/k_e)
    # this instance has ω_eff ≈ 10 (fast sweep, narrow bumps), so the O(k_e⁻²) remainder of the first-order location is
    # large at these gains: the vector error halves with k_e (0.55 at k_e = 40 in the table) but is not yet small
    assert out[40.0]['verr'] < 0.6 * out[20.0]['verr'] and out[40.0]['verr'] < 0.6
    assert out[40.0]['rel'] < 0.6 * out[20.0]['rel'] and out[40.0]['rel'] < 0.02
    assert out[40.0]['emax'] < 0.65 * out[20.0]['emax']                      # the periodic error is O(1/k_e)


# ---------------------------------------------------------------- Proposition K (collapsed configurations)
def _instance_A():
    """Instance A of experiments/exp_orientation_stability.py (n = 2, m = 3, elliptic reference, collapsed minimiser)."""
    N = 2
    net = L.SaturatedLiftedNetwork(N, a_max=3.0, c_max=2.0, s_min=0.3, s_max=1.5)
    f = lambda X: np.stack([np.sin(X[:, 1]) + 0.3 * X[:, 0], np.cos(X[:, 0]) - 0.2 * X[:, 1]], 1)
    xr = lambda t: np.array([1.2 * np.cos(t), 0.8 * np.sin(t)])
    xd = lambda t: np.array([-1.2 * np.sin(t), 0.8 * np.cos(t)])
    rng = np.random.default_rng(11)
    guess = net.unconstrained(rng.normal(0, 0.4, (1, N)), rng.uniform(-1, 1, (1, N)), rng.uniform(0.5, 1.0, 1))[0]
    return net, f, xr, xd, np.array([4.0, 4.0, 1.0, 1.0, 1.0]), guess


def test_proposition_K_hessian_splitting_and_critical_leak():
    """
    Proposition K(i): at a collapsed critical point Θ* = (ϑ̄, …, ϑ̄) the Hessian of J_m equals
    I_m ⊗ (σ₀G⁻¹ − R₂) + (1/m)(11ᵀ) ⊗ M₂ with R₂ the Hessian of the first-variation potential at the atom (frozen residual)
    and M₂ = ∫∇_ϑΦ̃ᵀ∇_ϑΦ̃ dρ; the separation block is σ₀G⁻¹ − R₂, the collective block ∇²J₁(ϑ̄).
    Proposition K(ii): σ_crit solves σ = λ_max(G^{1/2}R₂(ϑ̄(σ))G^{1/2}); at σ_crit the separation block of ∇²J_m is singular
    with the (m − 1)-fold zero of the Kronecker structure, while the collective block stays positive definite.
    """
    net, f, xr, xd, G, guess = _instance_A()
    m, p = 3, 5
    cl = PO.ClosedLoop(net, f, xr, xd, PER, 20.0, G, 0.5)
    orb = cl.orbit(N=1000)
    tb, _ = PO.collapsed_critical_point(cl, orb, guess)
    Th_star = np.tile(tb, (m, 1))
    assert np.linalg.norm(cl.grad_J_m(Th_star, orb)) < 1e-9                      # a collapsed critical point of J_m
    H = cl.hess_J_m(Th_star, orb)
    R2 = PO.first_variation_hessian(cl, orb, tb)
    M2 = np.zeros((p, p))
    for x in orb.xm:
        Jx = np.stack([net.grad_g(x, np.eye(net.n)[k], tb[None, :])[0] for k in range(net.n)])     # ∇_ϑΦ̃(x;ϑ̄), (n, p)
        M2 += Jx.T @ Jx
    M2 /= orb.N
    sep_block = 0.5 * np.diag(1.0 / G) - R2
    pred = np.kron(np.eye(m), sep_block) + np.kron(np.ones((m, m)) / m, M2)
    assert np.max(np.abs(H - pred)) < 1e-5 * (1 + np.max(np.abs(H)))
    # the two isotypic blocks
    Psep = PO.separation_projector(m, p)
    U, _, _ = np.linalg.svd(Psep)
    B_sep, B_col = U[:, :(m - 1) * p], U[:, (m - 1) * p:]
    e_sep = np.sort(np.linalg.eigvalsh(B_sep.T @ H @ B_sep))
    e_col = np.sort(np.linalg.eigvalsh(B_col.T @ H @ B_col))
    assert np.allclose(e_sep, np.sort(np.tile(np.linalg.eigvalsh(sep_block), m - 1)), atol=1e-5)   # (m − 1)-fold multiplets
    assert np.allclose(e_col, np.sort(np.linalg.eigvalsh(sep_block + M2)), atol=1e-5)
    assert e_sep[0] > 0 and e_col[0] > 0                                           # σ₀ = 0.5 > σ_crit: a strict minimum
    # (ii) the critical leak
    s_crit, tb_c, R2c = PO.sigma_crit(cl, orb, guess)
    sq = np.sqrt(G)
    assert abs(s_crit - np.linalg.eigvalsh(sq[:, None] * R2c * sq[None, :])[-1]) < 1e-8
    assert 0 < s_crit < 0.5
    Hc = cl.hess_J_m(np.tile(tb_c, (m, 1)), orb)                                    # cl.sigma0 is σ_crit after sigma_crit()
    ec_sep = np.sort(np.abs(np.linalg.eigvalsh(B_sep.T @ Hc @ B_sep)))
    ec_col = np.sort(np.linalg.eigvalsh(B_col.T @ Hc @ B_col))
    scale = np.max(np.abs(Hc))
    assert ec_sep[0] < 1e-5 * scale and ec_sep[1] < 1e-5 * scale and ec_sep[2] > 1e-3 * scale    # two-fold zero (m − 1 = 2)
    assert ec_col[0] > 1e-3 * scale                                                 # the collective block stays positive definite
