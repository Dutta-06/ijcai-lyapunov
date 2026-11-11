"""
Tier 3: Robustness Stress Test
===============================

Tests the graceful degradation of the projection-based architecture when
the fundamental assumptions (e.g., D_max bound) are deliberately violated.

This is done by injecting an unmodeled, high-frequency, high-amplitude
disturbance into the plant dynamics. The goal is to show that while tracking
error increases, the system remains bounded and does not catastrophically
diverge.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import PlantModel
from src.core.controller import AdaptiveController
from src.core.integrator import AdaptiveRK4Integrator, SimulationState

class DisturbedVanDerPol(PlantModel):
    """Van Der Pol plant with an injected unmodeled disturbance."""
    def __init__(self, mu=1.0, dist_amp=0.0, dist_freq=10.0, dist_start=0.0):
        super().__init__(n=2, name="Disturbed VDP")
        self.mu = mu
        self.dist_amp = dist_amp
        self.dist_freq = dist_freq
        self.dist_start = dist_start
        
    def f(self, x, t):
        # Base VDP dynamics
        dx = np.array([
            x[1],
            self.mu * (1 - x[0]**2) * x[1] - x[0]
        ])
        
        # Inject disturbance after dist_start
        if t >= self.dist_start:
            # High frequency disturbance that RBFs can't easily model
            dist = np.array([0, self.dist_amp * np.sin(self.dist_freq * t)])
            dx += dist
            
        return dx
        
    def x_ref(self, t):
        w = 1.0
        return np.array([np.sin(w*t), w*np.cos(w*t)])
        
    def x_ref_dot(self, t):
        w = 1.0
        return np.array([w*np.cos(w*t), -w**2*np.sin(w*t)])


def run_experiment():
    print("Setting up Robustness Stress Test...")
    
    # 1. Setup Disturbed Plant
    dist_start = 10.0 # Inject disturbance halfway through
    plant_normal = DisturbedVanDerPol(mu=1.0, dist_amp=0.0, dist_start=dist_start)
    plant_disturbed = DisturbedVanDerPol(mu=1.0, dist_amp=20.0, dist_freq=25.0, dist_start=dist_start)
    
    n = 2
    m = 16
    dt = 1e-3
    t_start = 0.0
    t_end = 20.0
    
    x0 = np.array([0.5, 0.5])
    grid_1d = np.linspace(-2.0, 2.0, 4)
    c0 = np.array(np.meshgrid(grid_1d, grid_1d)).T.reshape(-1, 2)
    W0 = np.zeros((m, n))
    sigma0 = np.ones(m) * 1.0
    
    init_state = SimulationState(x0, W0, c0, sigma0)
    
    ke = 10.0
    # Controller uses the SAME parameters for both cases
    def create_ctrl():
        return AdaptiveController(
            m=m, n=n, ke=ke, 
            gamma_W=50.0, gamma_c=10.0, gamma_sigma=2.0,
            W_max=30.0, c_max=15.0, sigma_min=0.2, sigma_max=5.0, delta_min=0.1
        )
    
    print("Running Baseline (No Disturbance)...")
    integrator_normal = AdaptiveRK4Integrator(plant_normal, create_ctrl(), dt=dt)
    res_normal = integrator_normal.simulate(init_state.copy(), t_start, t_end, record_interval=10, use_tqdm=True)
    
    print("Running Disturbed (Massive Disturbance after t=10)...")
    integrator_dist = AdaptiveRK4Integrator(plant_disturbed, create_ctrl(), dt=dt)
    res_dist = integrator_dist.simulate(init_state.copy(), t_start, t_end, record_interval=10, use_tqdm=True)
    
    # 5. Plot Results
    print("Generating plots...")
    os.makedirs('results', exist_ok=True)
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    data_normal = res_normal.to_arrays()
    data_dist = res_dist.to_arrays()
    t = data_normal['t']
    
    # Plot 1: Tracking Error
    axs[0].plot(t, data_normal['e_norm'], 'b-', label='Normal Operation')
    axs[0].plot(t, data_dist['e_norm'], 'r-', alpha=0.8, label='With Massive Unmodeled Disturbance')
    axs[0].axvline(dist_start, color='k', linestyle='--', label='Disturbance Injected')
    axs[0].set_yscale('log')
    axs[0].set_ylabel('Tracking Error Norm ||e(t)||')
    axs[0].set_title('Robustness to Violated Disturbance Bounds')
    axs[0].legend()
    axs[0].grid(True, alpha=0.3)
    
    # Plot 2: Parameter Bounds (Weights)
    W_norm_dist = np.linalg.norm(data_dist['W'], axis=(1,2), ord='fro')
    axs[1].plot(t, W_norm_dist, 'g-', label='||W(t)||_F under Disturbance')
    axs[1].axhline(30.0, color='r', linestyle='--', label='Projection Bound Ω_W (Hard Cap)')
    axs[1].axvline(dist_start, color='k', linestyle='--')
    axs[1].set_ylabel('Frobenius Norm of Weights')
    axs[1].set_xlabel('Time (s)')
    axs[1].set_title('Demonstration of Hard Parameter Boundedness Despite Disturbance')
    axs[1].legend()
    axs[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('results/tier3_robustness.png', dpi=300)
    print("Plot saved to results/tier3_robustness.png")

if __name__ == "__main__":
    run_experiment()
