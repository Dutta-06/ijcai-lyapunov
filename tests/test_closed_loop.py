"""Closed-loop invariants on a real (non-representable) plant."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import VanDerPolPlant
from src.core.controller import AdaptiveController
from src.core.integrator import AdaptiveRK4Integrator, SimulationState


def test_closed_loop_boundedness_and_hard_bounds():
    plant = VanDerPolPlant(mu=1.5, ref_amplitude=1.0, ref_omega=1.2)
    n, m, dt = 2, 9, 2e-3
    rng = np.random.default_rng(0)
    g = np.linspace(-1.5, 1.5, 3)
    c0 = np.array(np.meshgrid(g, g)).T.reshape(-1, 2)
    init = SimulationState(np.array([2.0, -1.0]), rng.normal(0, 3.0, size=(m, n)), c0, np.ones(m))
    ctrl = AdaptiveController(m=m, n=n, ke=5.0, gamma_W=20.0, gamma_c=5.0, gamma_sigma=1.0,
                              W_max=6.0, c_max=4.0, sigma_min=0.2, sigma_max=3.0)
    res = AdaptiveRK4Integrator(plant, ctrl, dt=dt).simulate(init, 0.0, 10.0, record_interval=5)
    d = res.to_arrays()
    assert np.all(np.isfinite(d['e_norm']))
    assert d['e_norm'][-len(d['e_norm']) // 4:].max() < 1.0, "no convergence to a small ball"
    # hard-bound invariants at every recorded step (initial W is outside the ball on purpose)
    for W, c, s in zip(d['W'][1:], d['c'][1:], d['sigma'][1:]):
        assert ctrl.params_inside_sets(W, c, s)


if __name__ == "__main__":
    test_closed_loop_boundedness_and_hard_bounds(); print("ok")
