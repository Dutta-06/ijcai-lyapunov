"""Second-order plants (src/core/sliding.py): the sliding-variable reduction. The averaged phase-plane loop is the gradient flow of L̄ in
the metric (1/m)G⁻¹ (finite differences); isotropic Gaussian kernels on a circular phase-plane orbit are over-parametrised — the
Gauss–Newton block has one zero eigenvalue per kernel (the restriction to the circle is a von Mises function of three effective
parameters); and the real servo under the sliding controller follows the averaged loop along the position orbit in the slow time t/k_s,
the error decreasing with k_s, with the position error at the prediction (d/dt + λ)⁻¹(d/dt + k_s)⁻¹r(x_m) and at its slaved form with
r/k_s (which drops the error filter's attenuation of the residual's harmonics, small here; `exp_servo_cogging.py` Part C shows where it
matters)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import sliding as SL
from src.core import langevin as L
from src.core import ot_leak as OT

PER = 2 * np.pi


def test_the_phase_plane_loop_is_the_gradient_flow_of_the_loss():
    net = SL.PhasePlaneNetwork(2, a_max=4.0, c_max=(3.0, 3.0), s_min=0.15, s_max=2.0)
    G = np.array([1.0, 0.25, 0.25, 0.25])
    loop = SL.AveragedSlidingLoop(net, lambda Z: -np.sin(Z[:, 0]) - 0.2 * Z[:, 1], lambda t: 1.2 * np.sin(t), lambda t: 1.2 * np.cos(t), PER, G, 0.01, n_quad=200)
    rng = np.random.default_rng(0)
    Th = net.unconstrained(rng.normal(0, 1, 5), rng.uniform(-1, 1, (5, 2)), rng.uniform(0.4, 1.0, 5))
    m = Th.shape[0]
    gL = np.zeros_like(Th); h = 1e-6
    for i in range(m):
        for c in range(4):
            E = np.zeros_like(Th); E[i, c] = h
            gL[i, c] = (loop.L_bar(Th + E) - loop.L_bar(Th - E)) / (2 * h)
    T = Th + 0.1
    assert np.abs(loop.drift(Th, T) - (-m * G[None, :] * gL - 0.01 * (Th - T))).max() < 1e-8


def test_phase_plane_kernels_on_a_circular_orbit_are_over_parametrised():
    net = SL.PhasePlaneNetwork(2, a_max=4.0, c_max=(3.0, 3.0), s_min=0.15, s_max=2.0)
    G = np.array([1.0, 0.25, 0.25, 0.25])
    loop = SL.AveragedSlidingLoop(net, lambda Z: -np.sin(Z[:, 0]), lambda t: 1.2 * np.sin(t), lambda t: 1.2 * np.cos(t), PER, G, 0.01, n_quad=400)
    rng = np.random.default_rng(1)
    K = 5
    atoms = net.unconstrained(rng.normal(0, 1, K), rng.uniform(-1, 1, (K, 2)), rng.uniform(0.4, 1.0, K))
    H, ev, V = loop.gauss_newton(np.ones(K) / K, atoms)
    assert np.sum(ev < 1e-15 * ev[-1]) == K and ev[K] > 1e-11 * ev[-1]         # one unexcited direction per kernel (the next 10⁴ times above)


def test_the_servo_follows_the_averaged_loop_along_the_position_orbit():
    tau = lambda X: 0.4 * np.sin(3 * np.asarray(X, float) + 0.3) + 0.1 * np.sin(7 * np.asarray(X, float))
    f = lambda X: -tau(X)
    G = np.array([1.0, 0.25, 0.25]); sig0 = 0.01
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=1.6, s_min=0.05, s_max=1.0)
    Th0 = net.unconstrained(np.array([[1.0], [1.4], [-1.9], [-0.2]]), np.array([[-0.9], [-0.33], [0.3], [0.78]]), np.array([0.25, 0.25, 0.27, 0.1]))
    T = Th0.copy()
    loop = OT.OTPriorLoop(net, f, lambda t: np.array([np.sin(t)]), PER, G, sig0, Th0, n_quad=800)
    errs = {}
    for ks in (25.0, 50.0):
        servo = SL.SlidingServo(net, f, np.sin, np.cos, lambda t: -np.sin(t), ks, 5.0, G, sig0)
        Tr, ts, e2 = servo.run(Th0, lambda Th: T, 10 * PER, dt=2e-3)
        Ta = Th0.copy(); dt = 0.02
        for _ in range(int(round(10 * PER / ks / dt))):
            k1 = loop.drift(Ta, T); k2 = loop.drift(Ta + 0.5 * dt * k1, T); Ta = Ta + dt * k2
        errs[ks] = np.linalg.norm((Tr - Th0) - (Ta - Th0)) / np.linalg.norm(Ta - Th0)
        if ks == 50.0:
            s, e = SL.filtered_error(loop.F - loop.f_hat(Tr), PER, ks, 5.0, exact=True)   # the prediction at the parameters reached
            s1, e1 = SL.filtered_error(loop.F - loop.f_hat(Tr), PER, ks, 5.0)                # with the slaved s = r/k_s
            late = ts > 9 * PER
            ratio = np.sqrt(np.mean(e2[late])) / np.sqrt(np.mean(e ** 2))
            ratio1 = np.sqrt(np.mean(e2[late])) / np.sqrt(np.mean(e1 ** 2))
    assert errs[25.0] < 0.03 and errs[50.0] < 0.6 * errs[25.0]                  # the averaged loop, the error O(1/k_s)
    assert abs(ratio - 1) < 0.01 and abs(ratio1 - 1) < 0.01  # the position error at both predictions (they part where r is high-frequency)
