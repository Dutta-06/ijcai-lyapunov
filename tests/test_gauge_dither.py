"""
Proposition D‴.36 (PARADIGM.md §5d): the gauge under dither. Checked here: the entropic leak on the masses is the Fisher–Rao gradient of
(ε/γ)Σβ log β (the weighted loop with it is still a gradient flow); the predicted diffusion of the masses for two classes equals the closed
form 2τT_α/(T_α + T_W)²; at the exact representation with W = w the leak keeps the configuration stationary and turns the gauge from a
zero mode into a stable one, its stationary spread the Ornstein–Uhlenbeck closed form √(D/(2κ)) for two classes; with γ = 0 the dithered
run leaves the weights at 1.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'experiments')))
import exp_cohesion_leak as E                                            # noqa: E402
import exp_gauge_dither as GD                                             # noqa: E402
from src.core import mass_transfer as MT                                  # noqa: E402


def test_the_entropic_leak_is_a_fisher_rao_gradient():
    net, base = E.setup()
    rng = np.random.default_rng(4)
    m = 5
    Th = net.unconstrained(rng.uniform(-2, 2, (m, 1)), rng.uniform(-2, 2, (m, 1)), rng.uniform(0.4, 0.9, m))
    beta = rng.uniform(0.5, 1.5, m); beta *= m / beta.sum()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.2, eps=0.03)
    dT, db = loop.drift(Th, beta)
    h = 1e-6
    gb = np.array([(loop.J(Th, beta + h * np.eye(m)[i]) - loop.J(Th, beta - h * np.eye(m)[i])) / (2 * h) for i in range(m)])
    assert np.abs(db + 0.2 * beta * (gb - beta @ gb / beta.sum())).max() < 1e-7 * np.abs(db).max()
    assert abs(db.sum()) < 1e-14 * m


def test_the_two_class_diffusion_has_the_closed_form():
    A = np.array([1.0, -0.7]); W = np.array([0.45, 0.55]); gamma, tau, amax = 0.05, 1e-4, 4.0
    D = GD.gauge_diffusion(W, gamma, tau, A=A, a_max=amax)
    a = A / W; ap = 1 - (a / amax) ** 2
    Ta = np.sum((A / (W ** 2 * ap)) ** 2) / E.G[0]
    Tw = (2 / gamma) * np.sum(1 / W)
    assert abs(D[0, 0] / (2 * tau * Ta / (Ta + Tw) ** 2) - 1) < 1e-12 and abs(D[0, 1] + D[0, 0]) < 1e-20
    assert GD.gauge_diffusion(W, 1e-6, tau, A=A)[0, 0] < 1e-6 * D[0, 0]  # no replicator, no gauge motion


def test_the_leak_pins_the_gauge_without_bias():
    net, base = E.setup()
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    A = E.A_BUMP; a = A * 3; ap = 1 - (a / 4.0) ** 2
    v = np.zeros(12); v[9] = 1.0; v[11] = -1.0                              # δlog β: class 1 up, class 3 down (Σδβ = 0)
    v[0] = -3 * A[0] / ap[0]; v[6] = +3 * A[2] / ap[2]                       # δα = −A δW/(W²a′), δW = δβ/3 = ±1/3 at β = 1
    out = {}
    for eps in (0.0, 1e-3):
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, eps=eps)

        def F(x):
            T = x[:9].reshape(3, 3); b = np.exp(x[9:])
            dT, db = loop.drift(T, b)
            return np.concatenate([dT.ravel(), db / b])
        x0 = np.concatenate([Th0.ravel(), np.zeros(3)])
        assert np.abs(F(x0)).max() < 1e-15                                  # stationary with and without the leak (W = w representable)
        out[eps] = (F(x0 + 1e-6 * v) - F(x0 - 1e-6 * v)) / 2e-6
    assert np.abs(out[0.0]).max() < 1e-9                                    # the gauge is a zero mode without the leak
    assert out[1e-3][9:] @ v[9:] < -1e-4                                    # and restored by it


def test_without_the_replicator_the_dithered_weights_stay_at_one():
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    Th, beta = loop.run_dithered(E.exact_rep(net, np.ones(3) / 3), np.ones(3), 500.0, 0.25, 1e-4, np.random.default_rng(0))
    assert np.abs(beta - 1).max() < 1e-14


def test_the_leaks_spread_has_the_two_class_closed_form():
    A = np.array([1.0, -0.7]); W = np.array([0.45, 0.55]); gamma, tau, eps, amax = 0.05, 1e-4, 2e-3, 4.0
    a = A / W; ap = 1 - (a / amax) ** 2
    Ta = np.sum((A / (W ** 2 * ap)) ** 2) / E.G[0]
    Tw = (2 / gamma) * np.sum(1 / W)
    D = 2 * tau * Ta / (Ta + Tw) ** 2
    kappa = eps * Tw / (Ta + Tw)                                            # the leak's stiffness along M, slowed by the positions' share
    sd = GD.leak_spread(W, gamma, tau, eps, A=A, a_max=amax)
    assert np.abs(sd / np.sqrt(D / (2 * kappa)) - 1).max() < 1e-10


def test_the_diffusion_is_the_projection_of_the_codes_own_drift():
    # round 26 (D3y-9): D against the spectral projector of the drift's Jacobian (src/core/gauge_drift.py), which ties the metric g to the code
    import exp_gauge_drift as X
    W = np.array([0.30, 0.25, 0.45])
    d, D, ev = X.drift_at(W[:2])
    assert np.abs(D - GD.gauge_diffusion(W, X.GAMMA, X.TAU)[:2, :2]).max() < 1e-4 * np.abs(D).max()


def test_run_dithered_has_the_noise_two_tau_g():
    # round 26 (D3y-9): one Euler–Maruyama step from the exact representation (where the drift vanishes): the positions' increments have
    # the covariance 2τG dt per particle, the weights move only at second order
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    tau, dt, n = 1e-4, 0.25, 4000
    rng = np.random.default_rng(1)
    inc = []
    for _ in range(n):
        Th, b = loop.run_dithered(Th0, np.ones(3), dt, dt, tau, rng)
        inc.append((Th - Th0).ravel())
        assert abs(b.sum() - 3) < 1e-12
    v = np.var(np.array(inc), axis=0) / (2 * tau * dt * np.tile(E.G, 3))
    assert np.abs(v - 1).max() < 0.12
