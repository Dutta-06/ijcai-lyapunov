"""
Tier 1: Chaotic Reference Tracking
===================================

Compares the proposed projection-based fully adaptive architecture
against classical σ-modification and e-modification baselines.
Uses the Lorenz system to generate a chaotic reference trajectory.

The primary goal is to demonstrate that the proposed method achieves
Uniform Ultimate Boundedness (UUB) without the persistent steady-state
bias characteristic of the leaky modification baselines.
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
from src.core.lyapunov_monitor import LyapunovMonitor
from src.baselines.sigma_mod import SigmaModController
from src.baselines.e_mod import EModController

def run_experiment():
    print("Setting up Lorenz chaotic tracking experiment...")
    
    # 1. Setup the Plant
    # The plant has slightly different parameters from the reference generator
    # to create an unknown dynamics tracking problem.
    plant = LorenzPlant(
        sigma=10.0, rho=28.0, beta=8.0/3.0,
        ref_sigma=10.5, ref_rho=28.5, ref_beta=3.0,
        ref_duration=50.0, ref_dt=1e-3, amplitude_scale=0.1
    )
    
    # 2. Setup common simulation parameters
    n = 3 # State dimension
    m = 27 # Number of RBF kernels (3x3x3 grid)
    dt = 1e-3
    t_start = 0.0
    t_end = 30.0
    
    # Initial conditions
    x0 = np.array([0.1, 0.1, 0.1])
    
    # Create an initial grid of centers
    grid_1d = np.linspace(-2.0, 2.0, 3)
    c0 = np.array(np.meshgrid(grid_1d, grid_1d, grid_1d)).T.reshape(-1, 3)
    W0 = np.zeros((m, n))
    sigma0 = np.ones(m) * 1.5
    
    init_state = SimulationState(x0, W0, c0, sigma0)
    
    # 3. Setup Controllers
    print("Initializing controllers...")
    
    ke = 10.0
    gamma_W = 50.0
    
    # Proposed Projection Controller
    ctrl_proj = AdaptiveController(
        m=m, n=n, ke=ke, 
        gamma_W=gamma_W, gamma_c=10.0, gamma_sigma=2.0,
        W_max=50.0, c_max=20.0, 
        sigma_min=0.5, sigma_max=5.0, delta_min=0.1
    )
    
    # Sigma-mod Controller
    ctrl_sigma = SigmaModController(
        m=m, n=n, ke=ke, gamma_W=gamma_W, sigma_mod=0.05
    )
    
    # E-mod Controller
    ctrl_emod = EModController(
        m=m, n=n, ke=ke, gamma_W=gamma_W, sigma_mod=0.05
    )
    
    # 4. Run Simulations
    results = {}
    controllers = {
        'Proposed (Proj)': ctrl_proj,
        'Baseline (sigma-mod)': ctrl_sigma,
        'Baseline (e-mod)': ctrl_emod
    }
    
    for name, ctrl in controllers.items():
        print(f"Running simulation for {name}...")
        integrator = AdaptiveRK4Integrator(plant, ctrl, dt=dt)
        
        # We don't use LyapunovMonitor here for the baselines since their V_dot is different,
        # but we can log e_norm naturally.
        res = integrator.simulate(
            initial_state=init_state,
            t_start=t_start, t_end=t_end,
            record_interval=10,
            use_tqdm=True
        )
        print()
        results[name] = res.to_arrays()
        
    # 5. Plot Results
    print("Generating plots...")
    plt.figure(figsize=(10, 6))
    
    colors = {'Proposed (Proj)': 'blue', 'Baseline (sigma-mod)': 'red', 'Baseline (e-mod)': 'orange'}
    
    for name, data in results.items():
        plt.plot(data['t'], data['e_norm'], label=name, color=colors[name], alpha=0.8)
        
    plt.yscale('log')
    plt.xlabel('Time (s)')
    plt.ylabel('Tracking Error Norm ||e(t)|| (log scale)')
    plt.title('Chaotic Reference Tracking (Lorenz System)')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    os.makedirs('results', exist_ok=True)
    plt.savefig('results/tier1_lorenz_tracking.png', dpi=300, bbox_inches='tight')
    print("Plot saved to results/tier1_lorenz_tracking.png")
    
if __name__ == "__main__":
    run_experiment()
