"""
Tier 3: Computational Overhead Benchmarking
===========================================

Benchmarks the wall-clock time per integration step of the proposed
fully adaptive architecture versus the classical fixed-grid architecture.

The goal is to demonstrate that while analytical Jacobians and projection
operators add overhead, the cost is marginal and suitable for real-time
control applications.
"""

import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import LorenzPlant
from src.core.controller import AdaptiveController
from src.baselines.sigma_mod import SigmaModController

def benchmark_controller(ctrl, plant, x_state, n_steps=1000):
    """Benchmark the raw control and adaptation law evaluation."""
    x = x_state
    W = np.zeros((ctrl.m, ctrl.n))
    c = np.zeros((ctrl.m, ctrl.n))
    sigma = np.ones(ctrl.m)
    
    t_val = 0.0
    x_ref = plant.x_ref(t_val)
    x_ref_dot = plant.x_ref_dot(t_val)
    e = x - x_ref
    
    start_time = time.time()
    for _ in tqdm(range(n_steps), desc="Benchmarking", leave=False):
        # In a real step, e, x_ref, etc. change, but for profiling the overhead
        # of the math operations, static values are fine.
        u = ctrl.control_law(x, x_ref, x_ref_dot, W, c, sigma)
        W_dot, c_dot, sigma_dot = ctrl.adaptation_laws(x, e, W, c, sigma)
        
        # Add a tiny mutation to avoid aggressive compiler optimizations
        W[0, 0] += 1e-10
        
    end_time = time.time()
    return (end_time - start_time) / n_steps * 1e6 # microseconds per step


def run_experiment():
    print("Setting up Computational Benchmarking...")
    
    plant = LorenzPlant()
    n = 3
    x_state = np.array([1.0, 1.0, 1.0])
    
    # We will test scaling with the number of kernels m
    m_values = [8, 27, 64, 125, 216] # 2^3, 3^3, 4^3, 5^3, 6^3 grids
    
    times_proposed = []
    times_baseline = []
    
    for m in m_values:
        print(f"Benchmarking m = {m} kernels...")
        
        # Proposed Controller
        ctrl_proj = AdaptiveController(
            m=m, n=n, ke=5.0, 
            gamma_W=10.0, gamma_c=10.0, gamma_sigma=1.0,
            W_max=50.0, c_max=20.0, sigma_min=0.5, sigma_max=5.0, delta_min=0.1
        )
        
        # Baseline Controller (σ-mod, fixed grid)
        ctrl_sigma = SigmaModController(
            m=m, n=n, ke=5.0, gamma_W=10.0, sigma_mod=0.1
        )
        
        # Warmup
        benchmark_controller(ctrl_proj, plant, x_state, n_steps=100)
        benchmark_controller(ctrl_sigma, plant, x_state, n_steps=100)
        
        # Benchmark
        t_proj = benchmark_controller(ctrl_proj, plant, x_state, n_steps=2000)
        t_sig = benchmark_controller(ctrl_sigma, plant, x_state, n_steps=2000)
        
        times_proposed.append(t_proj)
        times_baseline.append(t_sig)
        
    # Plot Results
    print("Generating plots...")
    os.makedirs('results', exist_ok=True)
    
    plt.figure(figsize=(8, 6))
    
    plt.plot(m_values, times_proposed, 'bo-', linewidth=2, label='Proposed (Full Adaptation + Projection)')
    plt.plot(m_values, times_baseline, 'rs--', linewidth=2, label='Baseline (Weights Only, No Projection)')
    
    plt.xlabel('Number of Kernels (m)')
    plt.ylabel('Computation Time per Step (μs)')
    plt.title('Computational Overhead Scaling')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Add an annotation about the relative cost
    avg_ratio = np.mean(np.array(times_proposed) / np.array(times_baseline))
    plt.annotate(f'Average overhead factor: {avg_ratio:.1f}x', 
                 xy=(0.05, 0.8), xycoords='axes fraction',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
                 
    plt.tight_layout()
    plt.savefig('results/tier3_computational.png', dpi=300)
    print("Plot saved to results/tier3_computational.png")

if __name__ == "__main__":
    run_experiment()
