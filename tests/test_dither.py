"""
Theorem J′ (PARADIGM.md §6f): the stationary law of the dithered finite-m closed loop.

1. The single-kernel instance's averaged data (∇J, b_fr, b_osc) agrees with the general implementation of Proposition H.7
   and Remark H.7a (averaging.expansion_terms, averaging.oscillation_correction, periodic.ClosedLoop.grad_J_m) on the
   saturated network with the width frozen.
2. The weighted Poisson solver reproduces the closed form of the Ornstein–Uhlenbeck case (quadratic J, constant b):
   g = G⁻¹b·Θ/τ₀, first-order mean shift (Gλ)⁻¹b, and the Green identity ∫φ g dπ = −π(b·∇ψ_φ).
3. The real loop at k_e = 40: the stationary means of the real loop are closer to the first-order prediction than to the
   Gibbs value, for the observables whose first-order shift is resolved.
"""
import os
import sys
import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import dither as Dt
from src.core import averaging as A
from src.core import langevin as L
from src.core import periodic as PO


def _instance(orientation=+1, N_theta=2000):
    """The instance of exp_dither_stationary.py: f(x) = e^{−(x−0.5)²/0.5} + 0.4x, x_m = 1.5[sin t′ + (1/3) sin 2t′],
    a_max = 3, c_max = 2.5, s = 0.6 frozen, G = diag(2, 2), σ₀ = 0.4, ϑ⁰ = 0, τ₀ = 0.1."""
    f = lambda X: np.exp(-(X - 0.5) ** 2 / 0.5) + 0.4 * X
    df = lambda X: np.exp(-(X - 0.5) ** 2 / 0.5) * (-2 * (X - 0.5) / 0.5) + 0.4
    return Dt.SingleKernelInstance(f, df, a_max=3.0, c_max=2.5, s=0.6, G=(2.0, 2.0), sigma0=0.4, theta0=(0.0, 0.0), tau0=0.1,
                                   X=1.5, omega=1.0, r2=1.0 / 3.0, orientation=orientation, N_theta=N_theta)


def test_single_kernel_averaged_data_matches_the_general_H7_implementation():
    inst = _instance(N_theta=2000)
    net = L.SaturatedLiftedNetwork(1, a_max=3.0, c_max=2.5, s_min=0.2, s_max=1.5)
    kappa = float(np.log((0.6 - 0.2) / (1.5 - 0.6)))                    # s = 0.6
    f_vec = lambda X: inst.f(X)
    jac_vec = lambda X: inst.df(X)[:, :, None]
    G3 = np.array([2.0, 2.0, 1e-12])                                      # width frozen (gain 0; 1e-12 avoids 1/G in the general code)
    cl = PO.ClosedLoop(net, f_vec, lambda t: np.array([inst.x_ref(t)]), lambda t: np.array([inst.x_ref_dot(t)]), inst.T0, 20.0, G3, inst.sigma0,
                       theta0=np.array([0.0, 0.0, kappa]))                 # leak target at the frozen width
    orb = cl.orbit(N=2000, sub=1)
    for (al, ga) in [(0.7, 0.3), (-0.4, 1.1), (1.5, -0.8)]:
        Th = np.array([[al, ga, kappa]])
        d = inst.averaged(np.array(al), np.array(ga))
        # ∇J (α, γ components) against the general gradient
        gJ = cl.grad_J_m(Th, orb)[0, :2]
        assert np.allclose(gJ, d['gradJ'], rtol=1e-8, atol=1e-10)
        # b_fr: G ∂_{ϑ'}(U_even + U_odd) by central differences in the general code
        fm = A.FrozenMeasure(net, Th)
        h = 1e-6
        E = np.vstack([Th + h * np.eye(3)[k][None, :] for k in range(3)] + [Th - h * np.eye(3)[k][None, :] for k in range(3)])
        v = sum(A.expansion_terms(f_vec, jac_vec, fm, orb, E)).reshape(2, 3)
        b_fr_general = G3 * (v[0] - v[1]) / (2 * h)
        assert np.allclose(b_fr_general[:2], d['b_fr'], rtol=2e-5, atol=1e-8)
        b_osc_general = A.oscillation_correction(net, f_vec, orb, Th, G3)[0]
        assert np.allclose(b_osc_general[:2], d['b_osc'], rtol=2e-4, atol=1e-8)
        # J itself against the general J_m (m = 1)
        assert abs(cl.J_m(Th, orb) - d['J']) < 1e-10


class _OU:
    """A mock instance with quadratic J = ½(λ_a α² + λ_g γ²) and constant b."""
    def __init__(self, lam, b, G, tau0):
        self.lam, self.bconst, self.G, self.tau0 = np.asarray(lam, float), np.asarray(b, float), np.asarray(G, float), float(tau0)

    def averaged(self, al, ga, with_hessian=False):
        J = 0.5 * (self.lam[0] * al ** 2 + self.lam[1] * ga ** 2)
        z = np.zeros_like(J)
        return dict(J=J, Lbar=J, gradJ=np.stack([self.lam[0] * al, self.lam[1] * ga]),
                    b_fr=np.stack([z + self.bconst[0], z + self.bconst[1]]), b_osc=np.stack([z, z]), uu=z)


def test_poisson_solver_reproduces_the_ornstein_uhlenbeck_closed_forms():
    lam, b, G, tau0 = (2.0, 0.5), (0.3, -0.2), (1.0, 2.0), 0.1
    inst = _OU(lam, b, G, tau0)
    sd = np.sqrt(tau0 / np.array(lam))
    grid = Dt.GibbsGrid(inst, (-7 * sd[0], 7 * sd[0]), (-7 * sd[1], 7 * sd[1]), n_al=201, n_ga=201)
    # g = G⁻¹b·Θ/τ₀ (zero π-mean)
    g_exact = (b[0] / G[0] * grid.AL + b[1] / G[1] * grid.GA) / tau0
    mask = grid.pi > 1e-3 * grid.pi.max()
    assert np.max(np.abs(grid.g - g_exact)[mask]) < 5e-3 * np.max(np.abs(g_exact)[mask])
    # first-order shift of the mean: (Gλ)⁻¹b
    for k, phi in enumerate((grid.AL, grid.GA)):
        assert abs(grid.first_order(phi) - b[k] / (G[k] * lam[k])) < 2e-3 * abs(b[k] / (G[k] * lam[k]))
        # Green identity: ∫φ g dπ = −π(b·∇ψ_φ)
        psi = grid.poisson(phi)
        gpsi = grid.grad(psi)
        rhs = -grid.expect(b[0] * gpsi[0] + b[1] * gpsi[1])
        assert abs(grid.first_order(phi) - rhs) < 2e-3 * abs(rhs)
    # constant b is a G-gradient: no current, zero entropy production, zero circulation coefficient
    assert grid.epr_coefficient() < 1e-4 * np.sum(np.array(b) ** 2 / np.array(G)) / tau0
    j1, _ = grid.current()
    assert np.max(np.abs(j1)) < 1e-3 * grid.pi.max() * np.linalg.norm(b)


def test_first_order_stationary_density_and_current_of_the_averaged_sde():
    """Proposition J′.2 on the grid: the stationary density of dΘ = (−G∇J + εb)ds + √(2τ₀)G^{1/2}dW is π(1 + εg) + O(ε²)
    and its current is εj₁ + O(ε²) (the discrete ε = 0 residual subtracted); OU case exact."""
    lam, b, G, tau0 = (2.0, 0.5), (0.3, -0.2), (1.0, 2.0), 0.1
    ou = Dt.GibbsGrid(_OU(lam, b, G, tau0), (-2.5, 2.5), (-4.5, 4.5), n_al=161, n_ga=161)
    m, j = ou.stationary_density(0.2)
    assert abs(np.sum(ou.AL * m) * ou.cell - 0.2 * b[0] / (G[0] * lam[0])) < 1e-4
    assert abs(np.sum(ou.GA * m) * ou.cell - 0.2 * b[1] / (G[1] * lam[1])) < 1e-4
    assert np.abs(j).max() < 1e-9 * m.max()                                   # a gradient drift carries no current
    inst = _instance(N_theta=512)
    grid = Dt.GibbsGrid(inst, (-3.5, 3.5), (-3.5, 3.5), n_al=161, n_ga=161)
    _, wedge1, _ = grid.circulation_prediction()
    assert abs(wedge1) > 1e-3                                                  # the instance carries a current
    m0, j0 = grid.stationary_density(0.0)
    w0 = float(np.sum(grid.AL * j0[1] - grid.GA * j0[0]) * grid.cell)
    ratios, l1s = [], []
    for eps in (0.1, 0.05, 0.025):
        m, j = grid.stationary_density(eps)
        w = float(np.sum(grid.AL * j[1] - grid.GA * j[0]) * grid.cell)
        ratios.append((w - w0) / (eps * wedge1))
        l1s.append(float(np.sum(np.abs(m - m0 * (1 + eps * grid.g))) * grid.cell) / eps ** 2)
    assert abs(ratios[-1] - 1) < 0.1 and abs(ratios[-1] - 1) < abs(ratios[0] - 1)     # → 1 linearly in ε
    assert max(l1s) / min(l1s) < 1.5                                                     # ‖m − π(1 + εg)‖₁ ∝ ε²


def test_real_loop_means_are_closer_to_the_first_order_prediction_than_to_gibbs():
    """Theorem J′(b)–(d) on the real loop at k_e = 20: the stationary means of α and γ are closer to π(φ) + ε∫φg dπ than
    to π(φ) (the shift of γ, 0.021, is resolved at ≥ 3 standard errors from the spread of the per-replica time averages),
    and the stationary tracking error k_e²E e² is within 10% of 2π(L̄)."""
    inst = _instance(N_theta=512)
    grid = Dt.GibbsGrid(inst, (-3.5, 3.5), (-3.5, 3.5), n_al=141, n_ga=141)
    ke = 20.0
    rng = np.random.default_rng(3)
    K = 8000
    al0, ga0 = grid.sample(K, rng, eps=1.0 / ke)
    mean, se, _ = Dt.simulate_stationary(inst, ke, K, T_slow=20.0, burn_slow=10.0, seed=7, al0=al0, ga0=ga0)
    obs = dict(alpha=grid.AL, gamma=grid.GA)
    for k, name in enumerate(('alpha', 'gamma')):
        p0 = grid.expect(obs[name]); p1 = grid.first_order(obs[name])
        assert abs(mean[k] - (p0 + p1 / ke)) < abs(mean[k] - p0), (name, mean[k], se[k], p0, p1 / ke)
    assert abs(grid.first_order(grid.GA) / ke) > 3 * se[1], ("the γ shift is not resolved", se[1])
    assert abs(mean[5] - 2 * grid.expect(grid.data['Lbar'])) < 0.10 * 2 * grid.expect(grid.data['Lbar'])


def test_corollary_J3_mean_field_minimiser_and_gibbs_sampler_consistency():
    """Corollary J′.3's numerics: the damped fixed-point iteration converges to the first-variation condition of F_τ
    (μ* ∝ exp(−V[μ*]/τ) on the grid, residual < 1e-9), and the MALA sampler of π_m reproduces the quadrature of π₁
    (m = 1) within three batch standard errors for E[α], E[γ]."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("exp_gibbs_m", os.path.join(os.path.dirname(__file__), '..', 'experiments', 'exp_gibbs_m.py'))
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
    E = importlib.util.module_from_spec(spec); spec.loader.exec_module(E)
    inst = _instance(N_theta=128)
    grid = Dt.GibbsGrid(inst, (-3.5, 3.5), (-3.5, 3.5), n_al=121, n_ga=121)
    mu, fh, hist = E.mean_field_minimiser(inst, grid)
    assert hist[-1] < 1e-9 and abs(np.sum(mu) * grid.cell - 1) < 1e-10
    # the fixed point: f̂ = ∫Φ̃ dμ* and μ* ∝ exp(−[V_L[μ*] + leak]/τ₀)
    AL, GA = grid.AL[..., None], grid.GA[..., None]
    a, c, _, _, _, _ = inst.phys(AL, GA)
    phi = a * np.exp(-(inst.xm - c) ** 2 / (2 * inst.s ** 2))
    assert np.max(np.abs(np.sum(mu[..., None] * phi, axis=(0, 1)) * grid.cell - fh)) < 1e-8
    # m = 1: MALA against quadrature
    rng = np.random.default_rng(11)
    mean, se, acc = E.mala(inst, 1, 2000, 4000, 0.02, rng)
    assert 0.5 < acc <= 1.0
    for k, phi_g in enumerate((grid.AL, grid.GA)):
        assert abs(mean[k] - grid.expect(phi_g)) < 3 * max(se[k], 2e-3), (k, mean[k], se[k], grid.expect(phi_g))


def test_time_symmetric_reference_has_no_oscillation_correction():
    """Remark J′.1: a reference with x_m(θ₀ + θ) = x_m(θ₀ − θ) (here a cosine series) has b_osc ≡ 0 and ⟨u∧u̇⟩ ≡ 0, and in
    one state dimension U_odd ≡ 0 for every reference, so its first-order drift is the performative term alone and the two
    orientations give identical averaged data; the sine series of the experiment's instance has b_osc ≠ 0 and odd."""
    f = lambda X: np.exp(-(X - 0.5) ** 2 / 0.5) + 0.4 * X
    df = lambda X: np.exp(-(X - 0.5) ** 2 / 0.5) * (-2 * (X - 0.5) / 0.5) + 0.4

    class Cosine(Dt.SingleKernelInstance):
        def x_ref(self, t):
            w = self.omega * self.orientation * np.asarray(t, float)
            return self.X * (np.cos(w) + self.r2 * np.cos(2 * w))

        def x_ref_dot(self, t):
            w = self.omega * self.orientation * np.asarray(t, float)
            return -self.X * self.omega * self.orientation * (np.sin(w) + 2 * self.r2 * np.sin(2 * w))

    al, ga = np.array([0.7, -0.4, 1.2]), np.array([0.3, 1.1, -0.5])
    data = {}
    for o in (+1, -1):
        inst = Cosine(f, df, a_max=3.0, c_max=2.5, s=0.6, G=(2.0, 2.0), sigma0=0.4, tau0=0.1, X=1.5, omega=1.0, r2=1.0 / 3.0, orientation=o, N_theta=512)
        inst.xm = inst.x_ref(inst.ths); inst.xmd = inst.x_ref_dot(inst.ths)
        data[o] = inst.averaged(al, ga)
        assert np.abs(data[o]['b_osc']).max() < 1e-12 and np.abs(data[o]['uu']).max() < 1e-12
    assert np.allclose(data[+1]['b_fr'], data[-1]['b_fr'], atol=1e-12)
    d_plus = _instance(orientation=+1, N_theta=512).averaged(al, ga)
    d_minus = _instance(orientation=-1, N_theta=512).averaged(al, ga)
    assert np.abs(d_plus['b_osc']).max() > 1e-2
    assert np.allclose(d_plus['b_osc'], -d_minus['b_osc'], atol=1e-12) and np.allclose(d_plus['b_fr'], d_minus['b_fr'], atol=1e-12)
