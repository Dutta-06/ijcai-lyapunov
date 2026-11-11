"""
Tier 3: Optimized Computational Overhead Benchmarking
=====================================================

Compares wall-clock time per integration step across three configurations:

  1. Baseline   — weights-only adaptation (γ_c=0, γ_σ=0)
  2. Dense      — full Γ matrices + dense Jacobian construction (original path)
  3. Optimized  — diagonal Γ exploit + inline Jacobian-vector products

The Dense and Optimized paths produce *mathematically identical* results;
the difference is purely computational.

Expected outcome: the Optimized path should collapse the ~286× overhead
factor (measured in tier3_computational.py) down to single-digit overhead
relative to the baseline.
"""

import os
import sys
import time
import numpy as np
from tqdm import tqdm

# Add the project root to the path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.controller import AdaptiveController
from src.core.plant_models import LorenzPlant

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def benchmark_controller(controller, plant, n_warmup=500, n_steps_per_trial=2000, n_trials=10, dt=0.01):
    """
    Benchmark the wall-clock time per integration step.

    Runs n_warmup steps to warm caches. Then runs n_trials trials, where
    each trial times a block of n_steps_per_trial steps. Returns the median
    time per step in microseconds across trials.
    """
    n = controller.n
    m = controller.m

    # Deterministic initial conditions
    rng = np.random.RandomState(42)
    x = rng.randn(n) * 0.5
    x_ref = np.zeros(n)
    x_ref_dot = np.zeros(n)
    W = rng.randn(m, n) * 0.1
    c = rng.randn(m, n) * 0.5
    sigma = np.ones(m) * 1.0

    # Warm-up (cache, JIT compilation, branch prediction)
    for _ in range(n_warmup):
        e = x - x_ref
        u = controller.control_law(x, x_ref, x_ref_dot, W, c, sigma)
        W_dot, c_dot, sigma_dot = controller.adaptation_laws(
            x, e, W, c, sigma)
        W = W + dt * W_dot
        c = (c.flatten() + dt * c_dot).reshape(m, n)
        sigma = np.clip(sigma + dt * sigma_dot, 0.1, 5.0)
        x = x + dt * (plant.f(x, 0) + u)

    trial_times = []
    
    # Run multiple trials to get stable timing
    for _ in range(n_trials):
        t0 = time.perf_counter()
        
        # Block of steps to amortize timing overhead
        for _ in range(n_steps_per_trial):
            e = x - x_ref
            
            # The core code we are profiling
            u = controller.control_law(x, x_ref, x_ref_dot, W, c, sigma)
            W_dot, c_dot, sigma_dot = controller.adaptation_laws(x, e, W, c, sigma)
            
            # Euler step to keep variables moving
            W = W + dt * W_dot
            c = (c.flatten() + dt * c_dot).reshape(m, n)
            sigma = np.clip(sigma + dt * sigma_dot, 0.1, 5.0)
            x = x + dt * (plant.f(x, 0) + u)
            
        t1 = time.perf_counter()
        
        # Microseconds per step for this trial
        trial_times.append((t1 - t0) * 1e6 / n_steps_per_trial)

    # Return the median time across trials
    return np.median(trial_times)


def run_experiment():
    print("Setting up Optimized Computational Overhead benchmark...")
    print("=" * 68)

    n = 3  # Lorenz system dimension
    plant = LorenzPlant()

    kernel_counts = [5, 10, 25, 50, 75, 125, 216]

    times_baseline = []
    times_dense = []
    times_optimized = []

    for m in tqdm(kernel_counts, desc="Benchmarking kernel counts"):
        # ── Baseline: weights only (scalar gammas, c & σ frozen) ──
        ctrl_base = AdaptiveController(
            m=m, n=n, ke=5.0,
            gamma_W=1.0, gamma_c=0.0, gamma_sigma=0.0,
            W_max=10.0, c_max=10.0,
            sigma_min=0.1, sigma_max=5.0,
            delta_min=0.01, proj_epsilon=0.05
        )
        t_base = benchmark_controller(ctrl_base, plant)
        times_baseline.append(t_base)

        # ── Dense: pass full Γ MATRICES → forces dense code path ──
        ctrl_dense = AdaptiveController(
            m=m, n=n, ke=5.0,
            gamma_W=1.0 * np.eye(m * n),
            gamma_c=0.5 * np.eye(m * n),
            gamma_sigma=0.3 * np.eye(m),
            W_max=10.0, c_max=10.0,
            sigma_min=0.1, sigma_max=5.0,
            delta_min=0.01, proj_epsilon=0.05
        )
        t_dense = benchmark_controller(ctrl_dense, plant)
        times_dense.append(t_dense)

        # ── Optimized: pass SCALAR gammas → triggers fast path ──
        ctrl_opt = AdaptiveController(
            m=m, n=n, ke=5.0,
            gamma_W=1.0, gamma_c=0.5, gamma_sigma=0.3,
            W_max=10.0, c_max=10.0,
            sigma_min=0.1, sigma_max=5.0,
            delta_min=0.01, proj_epsilon=0.05
        )
        t_opt = benchmark_controller(ctrl_opt, plant)
        times_optimized.append(t_opt)

        speedup = t_dense / t_opt if t_opt > 0 else float('inf')
        overhead = t_opt / t_base if t_base > 0 else float('inf')
        print(f"  m={m:>4d}:  baseline={t_base:>8.1f}us  "
              f"dense={t_dense:>8.1f}us  optimized={t_opt:>8.1f}us  "
              f"speedup={speedup:>6.1f}×  overhead={overhead:>5.1f}×")

    # ── Plotting ──────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Left panel: absolute wall-clock times
    ax = axes[0]
    ax.plot(kernel_counts, times_baseline, 'r--s',
            label='Baseline (W only)', linewidth=2, markersize=8)
    ax.plot(kernel_counts, times_dense, 'b--o',
            label='Dense (Original)', linewidth=2, markersize=8)
    ax.plot(kernel_counts, times_optimized, 'g-^',
            label='Optimized (Diagonal Γ + JVP)', linewidth=2, markersize=8)
    ax.set_xlabel('Number of Kernels (m)', fontsize=13)
    ax.set_ylabel('Computation Time per Step (μs)', fontsize=13)
    ax.set_title('Computational Overhead Scaling', fontsize=15)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Right panel: speedup and remaining overhead factors
    ax2 = axes[1]
    speedups = [d / o if o > 0 else 1.0
                for d, o in zip(times_dense, times_optimized)]
    overheads = [o / b if b > 0 else 1.0
                 for o, b in zip(times_optimized, times_baseline)]

    ax2.plot(kernel_counts, speedups, 'g-^',
             label='Speedup vs Dense', linewidth=2, markersize=8)
    ax2.plot(kernel_counts, overheads, 'm-d',
             label='Remaining Overhead vs Baseline',
             linewidth=2, markersize=8)
    ax2.axhline(y=1.0, color='k', linestyle='--', alpha=0.4, label='Parity')

    avg_speedup = np.mean(speedups)
    ax2.axhline(y=avg_speedup, color='g', linestyle=':',
                alpha=0.5, label=f'Avg speedup: {avg_speedup:.1f}×')

    ax2.set_xlabel('Number of Kernels (m)', fontsize=13)
    ax2.set_ylabel('Factor (×)', fontsize=13)
    ax2.set_title('Speedup & Remaining Overhead', fontsize=15)
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    results_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
    os.makedirs(results_dir, exist_ok=True)
    save_path = os.path.join(results_dir, 'tier3_computational_optimized.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nPlot saved → {save_path}")

    # ── Summary ───────────────────────────────────────────────────
    avg_overhead = np.mean(overheads)
    max_speedup = np.max(speedups)
    best_m = kernel_counts[np.argmax(speedups)]
    print(f"\n{'=' * 68}")
    print(f"  Average speedup over dense path  : {avg_speedup:>7.1f}×")
    print(f"  Peak speedup (m={best_m:>3d})             : {max_speedup:>7.1f}×")
    print(f"  Average remaining overhead       : {avg_overhead:>7.1f}×")
    print(f"  Original overhead (tier3)        : ~286.7×")
    print(f"{'=' * 68}")


if __name__ == '__main__':
    run_experiment()
