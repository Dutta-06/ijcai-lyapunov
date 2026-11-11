"""Checks for PARADIGM.md Theorem H (Langevin adaptive control)."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core import langevin as L
from src.core.plant_models import PlantModel


class _Scalar(PlantModel):
    def __init__(self):
        super().__init__(n=1, name='scalar')

    def f(self, x, t):
        return np.array([np.sin(2 * x[0]) + 0.3 * x[0]])

    def x_ref(self, t):
        return np.array([1.5 * np.sin(t) + 0.5 * np.sin(2 * t)])

    def x_ref_dot(self, t):
        return np.array([1.5 * np.cos(t) + np.cos(2 * t)])


def test_saturated_gradient_matches_finite_differences():
    rng = np.random.default_rng(0)
    for n in (1, 2):
        net = L.SaturatedLiftedNetwork(n, a_max=3.0, c_max=2.0, s_min=0.2, s_max=1.8)
        Th = rng.normal(size=(5, 2 * n + 1)) * 1.5
        x = rng.normal(size=n); e = rng.normal(size=n)
        G = net.grad_g(x, e, Th)
        h = 1e-6
        for i in range(5):
            for k in range(2 * n + 1):
                Tp, Tm = Th.copy(), Th.copy(); Tp[i, k] += h; Tm[i, k] -= h

                def feat(T):                                   # feature of particle i alone
                    a, c, s = net.physical(T)
                    return a[i] * np.exp(-np.sum((x - c[i]) ** 2) / (2 * s[i] ** 2))
                fd = e @ (feat(Tp) - feat(Tm)) / (2 * h)
                assert abs(G[i, k] - fd) < 1e-6 * (1 + abs(fd)), (n, i, k, G[i, k], fd)


def test_averaged_drift_is_the_mean_field_langevin_drift():
    """
    Theorem H, identification step: with the error channel at its quasi-steady value e = r(x)/k_e
    (Tikhonov) and averaging over one period of the reference, the drift of the adaptive law equals
    −(1/k_e) G ∇_ϑ V_L[μ](ϑ_i) − σ(ϑ_i − ϑ⁰), where V_L[μ](ϑ) = −⟨r_μᵀΦ̃(·;ϑ)⟩_ρ is the first variation of L̄.
    """
    rng = np.random.default_rng(1)
    plant = _Scalar(); net = L.SaturatedLiftedNetwork(1)
    ctrl = L.LangevinAdaptiveController(net, ke=10.0, G_diag=[2.0, 0.5, 0.5], sigma=0.05, theta0=[0.0, 0.0, 0.0])
    Th = rng.normal(size=(8, 3))
    period, N = 2 * np.pi, 400
    ts = np.arange(N) * period / N
    avg = np.zeros_like(Th)
    for t in ts:
        x = plant.x_ref(t)
        r = plant.f(x, t) - net.f_hat(x, Th)
        avg += ctrl.drift(x, r / ctrl.ke, Th) / N
    h = 1e-6
    gradV = np.zeros_like(Th)
    for k in range(3):
        Tp, Tm = Th.copy(), Th.copy(); Tp[:, k] += h; Tm[:, k] -= h
        gradV[:, k] = (L.first_variation(plant, net, Th, Tp, period, N) - L.first_variation(plant, net, Th, Tm, period, N)) / (2 * h)
    predicted = -(1 / ctrl.ke) * ctrl.G[None, :] * gradV - ctrl.sigma * (Th - ctrl.theta0[None, :])
    assert np.abs(avg - predicted).max() < 1e-7 * (1 + np.abs(predicted).max()), np.abs(avg - predicted).max()


def test_noiseless_symmetric_initialisation_never_separates_but_dither_does():
    """
    Exact permutation symmetry on the finite-m REAL closed loop: if all particles start at the same ϑ and
    τ = 0, the configuration stays a (moving) Dirac mass forever — every row is bit-identical at every
    snapshot.  With τ > 0 the rows separate within the first snapshot interval.
    """
    plant = _Scalar(); net = L.SaturatedLiftedNetwork(1)
    Th0 = np.tile(np.array([[0.2, 2.0, -1.0]]), (6, 1))
    ctrl0 = L.LangevinAdaptiveController(net, ke=10.0, G_diag=[2.0, 0.5, 0.5], sigma=0.05, tau=0.0)
    res0 = L.simulate(plant, ctrl0, plant.x_ref(0.0), Th0, T=5.0, dt=2e-3, rng=np.random.default_rng(2), record_every=50)
    assert all(np.all(S == S[0]) for S in res0['Th']), "rows must stay exactly equal without dither"
    assert np.abs(res0['Th'][-1][0] - Th0[0]).max() > 1e-3, "the common particle moves (invariant, not stationary)"
    ctrl1 = L.LangevinAdaptiveController(net, ke=10.0, G_diag=[2.0, 0.5, 0.5], sigma=0.05, tau=0.01)
    res1 = L.simulate(plant, ctrl1, plant.x_ref(0.0), Th0, T=0.2, dt=2e-3, rng=np.random.default_rng(2), record_every=50)
    assert np.max(np.ptp(res1['Th'][1], axis=0)) > 1e-3, "dither separates the rows within 0.1 s"


def test_first_variation_bound_used_in_holley_stroock():
    """|V_L[μ](ϑ)| ≤ K_Φ (F_ρ + K_Φ) for random and adversarial μ, ϑ (the oscillation bound of Theorem H(ii))."""
    rng = np.random.default_rng(5)
    plant = _Scalar(); net = L.SaturatedLiftedNetwork(1, a_max=4.0)
    period = 2 * np.pi
    F_rho = max(abs(plant.f(plant.x_ref(t), 0.0)[0]) for t in np.linspace(0, period, 400))
    bound = net.K_Phi * (F_rho + net.K_Phi)
    worst = 0.0
    for trial in range(60):
        scale = 1.0 if trial % 2 else 20.0                                   # large α ⇒ saturated weights (adversarial)
        Th = rng.normal(size=(16, 3)) * scale
        Te = rng.normal(size=(50, 3)) * scale
        worst = max(worst, float(np.max(np.abs(L.first_variation(plant, net, Th, Te, period, 200)))))
    assert worst <= bound, (worst, bound)


def test_averaging_error_of_the_formal_reduction_decreases_like_one_over_ke():
    """
    The reduction from the closed loop to (𝔄) is FORMAL (Remark H.0).  At frozen μ, integrate the true periodic
    error ė = r_μ(x_m + e) − k_e e (with x = x_m + e) and compare the period-averaged adaptive drift with the
    (𝔄) drift: the relative error must decay roughly like 1/k_e (it does not vanish at finite k_e).
    """
    rng = np.random.default_rng(6)
    plant = _Scalar(); net = L.SaturatedLiftedNetwork(1)
    Th = rng.normal(size=(8, 3))
    period, N = 2 * np.pi, 4000
    dt = period / N
    errs = []
    for ke in (10.0, 20.0, 40.0, 80.0):
        ctrl = L.LangevinAdaptiveController(net, ke=ke, G_diag=[4.0, 1.0, 1.0], sigma=0.0)
        e = np.zeros(1)
        for _ in range(3):                                                   # converge to the periodic solution
            avg = np.zeros_like(Th)
            for k in range(N):
                t = k * dt
                x = plant.x_ref(t) + e
                r = plant.f(x, t) - net.f_hat(x, Th)
                avg += ctrl.drift(x, e, Th) / N
                e = e + dt * (r - ke * e)                                    # ė = r(x_m + e) − k_e e
        ts = np.arange(N) * dt
        ref = np.zeros_like(Th)
        for t in ts:
            x = plant.x_ref(t); r = plant.f(x, t) - net.f_hat(x, Th)
            ref += ctrl.drift(x, r / ke, Th) / N                             # (𝔄) drift (algebraic identity)
        errs.append(np.abs(avg - ref).max() / np.abs(ref).max())
    errs = np.array(errs)
    assert np.all(np.diff(errs) < 0), errs
    r = errs[1:] / errs[:-1]
    assert np.all(r > 0.4) and np.all(r < 0.62), errs                        # halves when k_e doubles (≈ 1/k_e)
    assert errs[0] > 1e-3, errs                                              # non-zero at finite k_e


def test_lsi_lower_bound_formula():
    net = L.SaturatedLiftedNetwork(1, a_max=4.0)
    ctrl = L.LangevinAdaptiveController(net, ke=10.0, sigma=0.05, tau=0.01)
    rho = ctrl.lsi_constant_lower_bound(F_max=2.0)
    K = 4.0
    assert np.isclose(rho, (0.05 / 0.01) * np.exp(-2 * K * (2.0 + K) / (10.0 * 0.01)))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "passed")
