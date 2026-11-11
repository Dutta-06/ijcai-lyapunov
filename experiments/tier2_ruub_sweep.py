"""
Tier 2: Empirical vs. Theoretical R_UUB Comparison
===================================================

Sweeps over different values of feedback gain (k_e) and learning rates (Γ)
to compare the theoretical UUB radius to the empirically observed maximum
tracking error.

This confirms the scaling properties derived in Section 7 of the manuscript:
R_UUB = D_max / (2kₑ) + √(D²_max / (4k²ₑ) + C_Ω / kₑ)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
from time import time
from tqdm import tqdm

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import VanDerPolPlant
from src.core.controller import AdaptiveController
from src.core.integrator import AdaptiveRK4Integrator, SimulationState

def run_single_config(args):
    ke, gamma_scale = args
    
    # 1. Setup Plant
    plant = VanDerPolPlant(mu=1.0, ref_amplitude=1.0, ref_omega=1.0)
    
    n = 2
    m = 16
    dt = 2e-3
    t_start = 0.0
    t_end = 25.0
    
    # Fixed Initial Conditions for comparability
    x0 = np.array([1.0, -1.0])
    grid_1d = np.linspace(-2.0, 2.0, 4)
    c0 = np.array(np.meshgrid(grid_1d, grid_1d)).T.reshape(-1, 2)
    W0 = np.zeros((m, n))
    sigma0 = np.ones(m) * 1.0
    
    init_state = SimulationState(x0, W0, c0, sigma0)
    
    # Base projection bounds
    W_max = 20.0
    c_max = 10.0
    sigma_max = 5.0
    
    # 2. Controller
    ctrl = AdaptiveController(
        m=m, n=n, ke=ke,
        gamma_W=20.0 * gamma_scale, 
        gamma_c=5.0 * gamma_scale, 
        gamma_sigma=1.0 * gamma_scale,
        W_max=W_max, c_max=c_max,
        sigma_min=0.2, sigma_max=sigma_max, delta_min=0.1
    )
    
    # Estimate theoretical bound
    # These are very rough conservative estimates for D_max and C_omega
    # In practice, D_max should be derived from the specific function space
    epsilon_max_est = 0.5 
    D_max_est = ctrl.estimate_D_max(epsilon_max_est)
    
    # V_W, V_c, V_sigma are IFT velocity bounds. Assume some constants.
    V_W_est, V_c_est, V_sigma_est = 5.0, 2.0, 0.5
    C_omega_est = ctrl.estimate_C_omega(V_W_est, V_c_est, V_sigma_est)
    
    theoretical_ruub = ctrl.compute_theoretical_ruub(D_max_est, C_omega_est)
    
    # 3. Simulate
    integrator = AdaptiveRK4Integrator(plant, ctrl, dt=dt)
    res = integrator.simulate(
        init_state, t_start, t_end, record_interval=25
    )
    
    # 4. Extract metric
    data = res.to_arrays()
    ss_start = int(0.7 * len(data['t'])) # Last 30%
    max_ss_e = np.max(data['e_norm'][ss_start:])
    
    return ke, gamma_scale, theoretical_ruub, max_ss_e

def run_experiment():
    print("Starting R_UUB Sweep Experiment...")
    
    ke_values = [2.0, 5.0, 10.0, 20.0]
    gamma_scales = [0.1, 0.5, 1.0, 2.0, 5.0]
    
    configs = [(k, g) for k in ke_values for g in gamma_scales]
    
    start_time = time()
    pool = mp.Pool(processes=max(1, mp.cpu_count() - 1))
    results = list(tqdm(pool.imap_unordered(run_single_config, configs), total=len(configs), desc="R_UUB Sweep"))
    pool.close()
    pool.join()
    print(f"Completed in {time() - start_time:.2f} seconds.")
    
    # Process results into matrices for plotting
    empirical_ruub = np.zeros((len(ke_values), len(gamma_scales)))
    theoretical_ruub_mat = np.zeros((len(ke_values), len(gamma_scales)))
    
    for k_val, g_val, t_ruub, e_ruub in results:
        i = ke_values.index(k_val)
        j = gamma_scales.index(g_val)
        empirical_ruub[i, j] = e_ruub
        theoretical_ruub_mat[i, j] = t_ruub
        
    ratio_mat = empirical_ruub / theoretical_ruub_mat
    
    # Plotting
    os.makedirs('results', exist_ok=True)
    
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Empirical vs Theoretical scaling with ke (for fixed gamma_scale)
    mid_g_idx = len(gamma_scales) // 2
    g_val = gamma_scales[mid_g_idx]
    
    axs[0].plot(ke_values, empirical_ruub[:, mid_g_idx], 'o-', linewidth=2, label='Empirical max ||e||_ss')
    axs[0].plot(ke_values, theoretical_ruub_mat[:, mid_g_idx], 's--', linewidth=2, label='Theoretical R_UUB Estimate')
    axs[0].set_xlabel('Feedback Gain k_e')
    axs[0].set_ylabel('Tracking Error')
    axs[0].set_title(f'Error Scaling with k_e (Gamma Scale = {g_val})')
    axs[0].grid(True, alpha=0.3)
    axs[0].legend()
    
    # Plot 2: Heatmap of Empirical Error across Gamma and ke
    im = axs[1].imshow(empirical_ruub, origin='lower', aspect='auto', cmap='viridis')
    axs[1].set_xticks(range(len(gamma_scales)))
    axs[1].set_xticklabels([str(g) for g in gamma_scales])
    axs[1].set_yticks(range(len(ke_values)))
    axs[1].set_yticklabels([str(k) for k in ke_values])
    axs[1].set_xlabel('Learning Rate Scale factor (Γ)')
    axs[1].set_ylabel('Feedback Gain (k_e)')
    axs[1].set_title('Empirical Steady-State Error ||e||_ss')
    plt.colorbar(im, ax=axs[1])
    
    plt.tight_layout()
    plt.savefig('results/tier2_ruub_sweep.png', dpi=300)
    print("Plot saved to results/tier2_ruub_sweep.png")

if __name__ == "__main__":
    run_experiment()
