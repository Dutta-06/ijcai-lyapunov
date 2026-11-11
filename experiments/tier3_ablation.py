"""
Tier 3: Ablation Study
=======================

Evaluates the benefit of fully adaptive kernels (W, c, σ) versus partially
adaptive architectures (only W, or W and c).
Specifically, it seeks to answer: how many fewer kernels are needed to achieve
a target R_UUB when full adaptation is enabled?
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import LorenzPlant
from src.core.controller import AdaptiveController
from src.core.integrator import AdaptiveRK4Integrator, SimulationState

def run_experiment():
    print("Setting up Ablation Study experiment...")
    
    plant = LorenzPlant(ref_duration=20.0, amplitude_scale=0.1)
    
    n = 3
    dt = 1e-3
    t_start = 0.0
    t_end = 15.0
    
    # We will test a small grid that will perform poorly without adaptation,
    # but well with it.
    m = 27 # 3x3x3 grid
    
    x0 = np.array([0.1, 0.1, 0.1])
    grid_1d = np.linspace(-1.0, 1.0, 3)
    c0 = np.array(np.meshgrid(grid_1d, grid_1d, grid_1d)).T.reshape(-1, 3)
    W0 = np.zeros((m, n))
    sigma0 = np.ones(m) * 1.5
    
    init_state = SimulationState(x0, W0, c0, sigma0)
    
    ke = 10.0
    
    # Configuration 1: Fixed Centers & Bandwidths (gamma_c=0, gamma_sigma=0)
    ctrl_fixed = AdaptiveController(
        m=m, n=n, ke=ke, 
        gamma_W=50.0, gamma_c=0.0, gamma_sigma=0.0,
        W_max=50.0, c_max=20.0, sigma_min=0.5, sigma_max=5.0, delta_min=0.05
    )
    
    # Configuration 2: Adaptive Centers only (gamma_sigma=0)
    ctrl_adapt_c = AdaptiveController(
        m=m, n=n, ke=ke, 
        gamma_W=50.0, gamma_c=20.0, gamma_sigma=0.0,
        W_max=50.0, c_max=20.0, sigma_min=0.5, sigma_max=5.0, delta_min=0.05
    )
    
    # Configuration 3: Full Adaptation
    ctrl_full = AdaptiveController(
        m=m, n=n, ke=ke, 
        gamma_W=50.0, gamma_c=20.0, gamma_sigma=5.0,
        W_max=50.0, c_max=20.0, sigma_min=0.5, sigma_max=5.0, delta_min=0.05
    )
    
    controllers = {
        'Fixed (W only)': ctrl_fixed,
        'Adaptive Centers (W, c)': ctrl_adapt_c,
        'Full Adaptation (W, c, σ)': ctrl_full
    }
    
    results = {}
    
    for name, ctrl in controllers.items():
        print(f"Running simulation for {name}...")
        integrator = AdaptiveRK4Integrator(plant, ctrl, dt=dt)
        # Using a fresh copy of init_state each time is crucial
        res = integrator.simulate(
            initial_state=init_state.copy(),
            t_start=t_start, t_end=t_end, record_interval=20,
            use_tqdm=True
        )
        print()
        results[name] = res.to_arrays()
        
    print("Generating plots...")
    os.makedirs('results', exist_ok=True)
    
    plt.figure(figsize=(10, 6))
    
    colors = ['gray', 'orange', 'blue']
    
    for (name, data), color in zip(results.items(), colors):
        plt.plot(data['t'], data['e_norm'], label=name, color=color, linewidth=2, alpha=0.8)
        
    plt.yscale('log')
    plt.xlabel('Time (s)')
    plt.ylabel('Tracking Error Norm ||e(t)||')
    plt.title(f'Ablation Study: Architecture Capability with m={m} kernels')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('results/tier3_ablation.png', dpi=300)
    print("Plot saved to results/tier3_ablation.png")

if __name__ == "__main__":
    run_experiment()
