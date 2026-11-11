"""
Tier 2: Parameter Migration Visualization
=======================================

Visualizes the evolution of adaptive parameters (weights, centers, bandwidths)
over time, showing how they migrate and reshape to learn the plant dynamics,
bounded by the projection sets.

This is the "wow" figure showing the kernels physically moving.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import VanDerPolPlant
from src.core.controller import AdaptiveController
from src.core.integrator import AdaptiveRK4Integrator, SimulationState

def run_experiment():
    print("Setting up Parameter Evolution experiment...")
    
    # 1. Setup Plant (Van der Pol - 2D state space is easy to visualize)
    plant = VanDerPolPlant(mu=2.0, ref_amplitude=1.5, ref_omega=1.0)
    
    n = 2
    m = 9 # 3x3 grid for clear visualization
    dt = 2e-3
    t_start = 0.0
    t_end = 15.0
    
    # Initial conditions
    x0 = np.array([0.0, 0.0])
    
    # Initial centers clustered near origin
    grid_1d = np.linspace(-0.5, 0.5, 3)
    c0 = np.array(np.meshgrid(grid_1d, grid_1d)).T.reshape(-1, 2)
    W0 = np.random.normal(0, 0.1, size=(m, n))
    sigma0 = np.ones(m) * 0.5 # Start with small bandwidths
    
    init_state = SimulationState(x0, W0, c0, sigma0)
    
    # 2. Controller
    ke = 5.0
    c_max = 3.0 # Set a tight bound to see projection activate
    ctrl = AdaptiveController(
        m=m, n=n, ke=ke, 
        gamma_W=10.0, gamma_c=10.0, gamma_sigma=2.0,
        W_max=20.0, c_max=c_max, 
        sigma_min=0.2, sigma_max=2.0, delta_min=0.1
    )
    
    # 3. Simulate
    print("Running simulation...")
    integrator = AdaptiveRK4Integrator(plant, ctrl, dt=dt)
    res = integrator.simulate(
        init_state, t_start, t_end, record_interval=50, # Record frequently for smooth animation/plotting
        use_tqdm=True
    )
    print()
    
    # 4. Extract data
    data = res.to_arrays()
    t = data['t']
    x = data['x']
    c_history = data['c'] # Shape (time_steps, m, n)
    sigma_history = data['sigma'] # Shape (time_steps, m)
    
    # 5. Plot Results
    print("Generating parameter evolution plots...")
    os.makedirs('results', exist_ok=True)
    
    # We'll plot snapshots at a few key times
    times_to_plot = [0.0, 3.0, 7.0, 15.0]
    indices = [np.searchsorted(t, tp) for tp in times_to_plot]
    
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))
    axs = axs.flatten()
    
    for i, idx in enumerate(indices):
        idx = min(idx, len(t)-1) # Safety bound
        current_time = t[idx]
        current_x = x[:idx+1] # Trajectory up to this point
        current_c = c_history[idx]
        current_sigma = sigma_history[idx]
        
        ax = axs[i]
        
        # Plot state trajectory
        if len(current_x) > 0:
            ax.plot(current_x[:, 0], current_x[:, 1], 'k-', alpha=0.3, label='State Trajectory')
            ax.plot(current_x[-1, 0], current_x[-1, 1], 'ko', markersize=8) # Current position
        
        # Plot projection boundary
        boundary = Circle((0, 0), c_max, color='r', fill=False, linestyle='--', label='Projection Boundary Ω_c')
        ax.add_patch(boundary)
        
        # Plot kernels (center as dot, bandwidth as circle)
        for j in range(m):
            ax.plot(current_c[j, 0], current_c[j, 1], 'bo', markersize=6)
            kernel_circle = Circle(current_c[j], current_sigma[j], color='blue', fill=True, alpha=0.1)
            ax.add_patch(kernel_circle)
            
        ax.set_xlim(-c_max-0.5, c_max+0.5)
        ax.set_ylim(-c_max-0.5, c_max+0.5)
        ax.set_title(f'Time t = {current_time:.1f}s')
        ax.set_xlabel('State x₁')
        ax.set_ylabel('State x₂')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        
        if i == 0:
            ax.legend(loc='upper right')
            
    plt.suptitle('RBF Kernel Migration in State Space over Time', fontsize=16)
    plt.tight_layout()
    plt.savefig('results/tier2_parameter_evolution.png', dpi=300)
    print("Plot saved to results/tier2_parameter_evolution.png")

if __name__ == "__main__":
    run_experiment()
