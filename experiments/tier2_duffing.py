"""
Tier 2: Duffing Oscillator Tracking
====================================

Uses the Duffing oscillator with time-varying stiffness to stress the
non-autonomous drift bounds of the Lyapunov convergence proof.

This demonstrates that the controller can handle plant dynamics f(x,t)
that change over time, forcing the ideal parameters W*, c*, σ* to migrate.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import DuffingPlant
from src.core.controller import AdaptiveController
from src.core.integrator import AdaptiveRK4Integrator, SimulationState

def run_experiment():
    print("Setting up Duffing oscillator tracking experiment...")
    
    # 1. Setup the Plant (Time-varying stiffness)
    plant = DuffingPlant(
        delta=0.3, alpha0=1.0, beta_coeff=1.0, gamma=1.0, omega=1.2,
        delta_alpha=0.5, omega_drift=0.5, # High drift to stress the controller
        ref_amplitude=2.0, ref_omega=1.0
    )
    
    n = 2
    m = 25 # 5x5 grid
    dt = 2e-3
    t_start = 0.0
    t_end = 40.0
    
    # Initial conditions
    x0 = np.array([0.5, -0.5])
    
    grid_1d = np.linspace(-3.0, 3.0, 5)
    c0 = np.array(np.meshgrid(grid_1d, grid_1d)).T.reshape(-1, 2)
    W0 = np.zeros((m, n))
    sigma0 = np.ones(m) * 1.5
    
    init_state = SimulationState(x0, W0, c0, sigma0)
    
    # 2. Controller
    ke = 15.0
    ctrl = AdaptiveController(
        m=m, n=n, ke=ke, 
        gamma_W=100.0, gamma_c=20.0, gamma_sigma=5.0,
        W_max=50.0, c_max=20.0, 
        sigma_min=0.2, sigma_max=10.0, delta_min=0.1
    )
    
    # 3. Simulate
    print("Running simulation (this may take a moment)...")
    integrator = AdaptiveRK4Integrator(plant, ctrl, dt=dt)
    res = integrator.simulate(
        init_state, t_start, t_end, record_interval=10,
        use_tqdm=True
    )
    print()
    
    # 4. Extract data
    data = res.to_arrays()
    t = data['t']
    x = data['x']
    x_ref = np.array([plant.x_ref(ti) for ti in t])
    e_norm = data['e_norm']
    
    # Calculate true time-varying stiffness for plotting
    alpha_t = np.array([plant._alpha(ti) for ti in t])
    
    # 5. Plot Results
    print("Generating plots...")
    os.makedirs('results', exist_ok=True)
    
    fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    
    # Plot 1: Tracking Performance
    axs[0].plot(t, x_ref[:, 0], 'k--', label='Reference x₁', alpha=0.7)
    axs[0].plot(t, x[:, 0], 'b-', label='Actual x₁', alpha=0.9)
    axs[0].set_ylabel('Position x₁')
    axs[0].legend(loc='upper right')
    axs[0].grid(True, alpha=0.3)
    axs[0].set_title('Duffing Oscillator Tracking under Time-Varying Dynamics')
    
    # Plot 2: Time-varying Stiffness (The Disturbance)
    axs[1].plot(t, alpha_t, 'r-', linewidth=2, label='Stiffness α(t)')
    axs[1].set_ylabel('Stiffness α')
    axs[1].legend(loc='upper right')
    axs[1].grid(True, alpha=0.3)
    
    # Plot 3: Error Norm
    axs[2].plot(t, e_norm, 'g-', linewidth=2, label='||e(t)||')
    axs[2].set_ylabel('Tracking Error Norm')
    axs[2].set_xlabel('Time (s)')
    axs[2].set_yscale('log')
    axs[2].legend(loc='upper right')
    axs[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('results/tier2_duffing.png', dpi=300)
    print("Plot saved to results/tier2_duffing.png")

if __name__ == "__main__":
    run_experiment()
