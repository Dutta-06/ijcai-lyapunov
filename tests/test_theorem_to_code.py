"""
Theorem-to-code gate.

On an exactly representable plant with known constant θ* (oracle setting) the
two facts the proof rests on are evaluated along simulated trajectories:

    identity:   V̇ = −kₑ‖e‖² + eᵀd + cross         (must hold to fp precision)
    property P: cross ≤ 0                            (must hold at every sample)

and their consequence, the RMS bound with the measured disturbance:

    RMS_{[T₁,T]}(‖e‖) ≤ D/(2kₑ) + √(D²/(4kₑ²) + V(T₁)/(kₑ (T−T₁))),   D = sup ‖d‖.

The legacy centre-repulsion option is also exercised to document that it is
outside the proof (it violates property P).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.controller import AdaptiveController
from src.core.integrator import SimulationState
from src.core.oracle import (RepresentablePlant, sample_admissible_ideal, run_oracle,
                             rms_bound, rms_over_window)

BOUNDS = dict(W_max=10.0, c_max=6.0, sigma_min=0.3, sigma_max=3.0, proj_epsilon=0.2)


def _setup(seed, m=8, n=2, **ctrl_extra):
    rng = np.random.default_rng(seed)
    Ws, cs, ss = sample_admissible_ideal(m, n, BOUNDS['W_max'], BOUNDS['c_max'], BOUNDS['sigma_min'],
                                         BOUNDS['sigma_max'], BOUNDS['proj_epsilon'], rng)
    plant = RepresentablePlant(Ws, cs, ss)
    ctrl = AdaptiveController(m=m, n=n, ke=3.0, gamma_W=5.0, gamma_c=2.0, gamma_sigma=1.0,
                              delta_min=0.5, **BOUNDS, **ctrl_extra)
    assert ctrl.ideal_params_admissible(Ws, cs, ss)
    x0 = rng.uniform(-1, 1, size=n)
    c0 = cs + rng.normal(0, 0.4, size=(m, n))
    s0 = np.clip(ss + rng.normal(0, 0.3, size=m), BOUNDS['sigma_min'] + BOUNDS['proj_epsilon'],
                 BOUNDS['sigma_max'] - BOUNDS['proj_epsilon'])
    return ctrl, plant, SimulationState(x0, np.zeros((m, n)), c0, s0)


def test_identity_and_property_P_hold_for_proved_algorithm():
    for seed in range(3):
        ctrl, plant, st0 = _setup(seed)
        res = run_oracle(ctrl, plant, st0, T=8.0, dt=1e-3, every=5)
        ident = np.max(np.abs(res['Vdot'] - res['rhs']) / (1.0 + np.abs(res['Vdot'])))
        assert ident < 1e-12, f"seed {seed}: V̇ decomposition identity broken ({ident})"
        assert np.max(res['cross']) <= 1e-10, f"seed {seed}: property P violated ({np.max(res['cross'])})"
        assert np.all(res['inside']), "hard parameter bounds violated"


def test_rms_bound_with_measured_disturbance():
    for seed in range(3):
        ctrl, plant, st0 = _setup(seed)
        T, T1 = 12.0, 4.0
        res = run_oracle(ctrl, plant, st0, T=T, dt=1e-3, every=5)
        rms, mask = rms_over_window(res, T1)
        D = float(np.max(res['d_norm'][mask]))
        V_T1 = float(res['V'][mask][0])
        bound = rms_bound(D, ctrl.ke, V_T1, T - T1)
        assert rms <= bound * (1 + 1e-6), f"seed {seed}: RMS {rms} > bound {bound}"


def test_legacy_separation_violates_property_P():
    ctrl, plant, st0 = _setup(0, enforce_separation=True)
    res = run_oracle(ctrl, plant, st0, T=4.0, dt=1e-3, every=5)
    assert np.max(res['cross']) > 1e-8, "expected the legacy repulsion to break property P"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "passed")
