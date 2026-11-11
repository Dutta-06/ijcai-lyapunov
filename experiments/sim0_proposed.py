"""
Simulation 0: Proposed Self-Adaptive RBF Controller
=====================================================

Runs the proposed projection-based adaptive RBF controller on the EXACT
same 1D benchmark used by all four kernel baselines for a fair comparison.

Same plant dynamics, same Lorenz reference, same ke=5.0, same m=27,
same dt=0.001, same T=30s.
"""

import os
import sys
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
from tqdm import tqdm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.controller import AdaptiveController
from src.core.projection import ProjectionOperator, BandwidthProjection, CenterProjection


def lorenz(t, state, sigma=10, rho=28, beta=8/3):
    x, y, z = state
    return [sigma*(y-x), x*(rho-z)-y, x*y-beta*z]

def f_true(x, t):
    """Same 1D plant dynamics as all baselines."""
    return -0.5*x[0] + 2*np.sin(x[0]) + 0.3*np.cos(2*t)*x[0]

def generate_lorenz_reference(T=30, dt=0.001):
    t_span = [0, T]
    t_eval = np.arange(0, T, dt)
    sol = solve_ivp(lorenz, t_span, [1.0, 1.0, 1.0], t_eval=t_eval, method='RK45')
    return t_eval, sol.y[0]


def run_proposed(T=30, dt=0.001, m=27):
    t_eval, xm = generate_lorenz_reference(T, dt)
    N = len(t_eval)

    n = 1  # 1D state for this benchmark

    # Controller with same ke=5.0 and m=27 as all baselines
    ctrl = AdaptiveController(
        m=m, n=n, ke=5.0,
        gamma_W=50.0, gamma_c=20.0, gamma_sigma=5.0,
        W_max=50.0, c_max=20.0,
        sigma_min=0.2, sigma_max=10.0, delta_min=0.05
    )

    # Initial conditions
    x = np.array([1.0])
    W = np.zeros((m, n))
    c = np.linspace(-3.0, 3.0, m).reshape(m, n)  # Same uniform grid start
    sigma = np.ones(m) * 1.5

    e_history = []
    t_history = []
    comp_time = []

    for i in tqdm(range(N-1), desc="Proposed (Proj)"):
        t = t_eval[i]

        # Reference
        x_ref = np.array([xm[i]])
        xm_dot = np.array([(xm[i+1] - xm[i]) / dt]) if i < N-2 else np.array([0.0])
        e = x - x_ref

        # Control law: u = x_ref_dot - ke*e - f_hat
        u_vec = ctrl.control_law(x, x_ref, xm_dot, W, c, sigma)
        u_scalar = u_vec[0]

        # Plant step (RK4) -- identical to all baselines
        k1 = f_true(x, t) + u_scalar
        k2 = f_true(x + dt/2*k1, t + dt/2) + u_scalar
        k3 = f_true(x + dt/2*k2, t + dt/2) + u_scalar
        k4 = f_true(x + dt*k3, t + dt) + u_scalar
        x = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        # Adaptation laws (continuous-time, Euler integration for params)
        t0 = time.time()
        W_dot, c_dot, sigma_dot = ctrl.adaptation_laws(x, e, W, c, sigma)
        W = W + dt * W_dot
        c = c + dt * c_dot.reshape(m, n)
        sigma = sigma + dt * sigma_dot

        # Apply projection clips (safety net)
        W = ctrl.proj_W.hard_clip_matrix(W)
        c_flat = ctrl.proj_c.hard_clip(c.flatten(), n)
        c = c_flat.reshape(m, n)
        sigma = ctrl.proj_sigma.hard_clip(sigma)
        comp_time.append(time.time() - t0)

        e_val = x[0] - xm[i+1]
        e_history.append(abs(e_val))
        t_history.append(t)

    return np.array(t_history), np.array(e_history), np.array(comp_time)


if __name__ == "__main__":
    print("Running Proposed (Projection-based) simulation on 1D benchmark...")
    os.makedirs('results', exist_ok=True)
    t, e, ct = run_proposed(T=30, dt=0.001, m=27)
    np.save("results/proposed_error.npy", e)
    np.save("results/proposed_time.npy", t)
    np.save("results/proposed_comptime.npy", ct)
    print(f"Done. Mean |e|: {e[5000:].mean():.6f}")
    plt.semilogy(t, e, label='Proposed (Proj)', color='green')
    plt.xlabel('Time (s)'); plt.ylabel('||e(t)||')
    plt.title('Proposed Architecture -- Tracking Error')
    plt.legend(); plt.savefig('results/sim0_proposed.png', dpi=150)
    print("Plot saved to results/sim0_proposed.png")
