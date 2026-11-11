"""
Master Runner: All Kernel Baseline Comparisons
================================================

Runs all 5 simulations sequentially and generates the final comparison plot.
Recommended order: Proposed -> NTK -> Exact GP -> Sparse GP -> Deep GP

Expected total runtime (CPU): 90-180 minutes
Expected total runtime (GPU): 15-30 minutes
"""

import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
os.makedirs('results', exist_ok=True)

print("=" * 70)
print("  KERNEL BASELINE COMPARISON SUITE")
print("=" * 70)

overall_start = time.time()

# --- Sim 0: Proposed Architecture ---
print("\n[1/5] Running Proposed (Projection-based) controller...")
t0 = time.time()
from sim0_proposed import run_proposed
import numpy as np
t, e, ct = run_proposed(T=30, dt=0.001, m=27)
np.save("results/proposed_error.npy", e)
np.save("results/proposed_time.npy", t)
np.save("results/proposed_comptime.npy", ct)
print(f"  Done in {time.time()-t0:.1f}s. Mean SS |e|: {e[5000:].mean():.6f}")

# --- Sim 4: NTK (fastest baseline) ---
print("\n[2/5] Running NTK-based control...")
t0 = time.time()
from sim4_ntk_control import run_ntk_control
t, e, ct = run_ntk_control(T=30, dt=0.001, width=256)
np.save("results/ntk_error.npy", e)
np.save("results/ntk_time.npy", t)
np.save("results/ntk_comptime.npy", ct)
print(f"  Done in {time.time()-t0:.1f}s. Mean SS |e|: {e[5000:].mean():.6f}")

# --- Sim 1: Exact GP ---
print("\n[3/5] Running Exact GP control...")
t0 = time.time()
from sim1_exact_gp import run_exact_gp
t, e, ct = run_exact_gp(T=30, dt=0.001, max_gp_points=500)
np.save("results/exact_gp_error.npy", e)
np.save("results/exact_gp_time.npy", t)
np.save("results/exact_gp_comptime.npy", ct)
print(f"  Done in {time.time()-t0:.1f}s. Mean SS |e|: {e[5000:].mean():.6f}")

# --- Sim 2: Sparse GP (Fixed) ---
print("\n[4/5] Running Sparse GP (Fixed Inducing Points)...")
t0 = time.time()
from sim2_sparse_gp import run_sparse_gp_fixed
t, e, ct = run_sparse_gp_fixed(T=30, dt=0.001, m=27)
np.save("results/sparse_gp_fixed_error.npy", e)
np.save("results/sparse_gp_fixed_time.npy", t)
np.save("results/sparse_gp_fixed_comptime.npy", ct)
print(f"  Done in {time.time()-t0:.1f}s. Mean SS |e|: {e[5000:].mean():.6f}")

# --- Sim 3: Deep GP (slowest) ---
print("\n[5/5] Running Deep GP (2-layer)...")
t0 = time.time()
from sim3_deep_gp import run_deep_gp
t, e, ct = run_deep_gp(T=30, dt=0.001)
np.save("results/deep_gp_error.npy", e)
np.save("results/deep_gp_time.npy", t)
np.save("results/deep_gp_comptime.npy", ct)
print(f"  Done in {time.time()-t0:.1f}s. Mean SS |e|: {e[5000:].mean():.6f}")

# --- Final Comparison Plot ---
print("\n" + "=" * 70)
print("  Generating comparison plot...")
print("=" * 70)
plot_script_path = os.path.join(os.path.dirname(__file__), "plot_comparison.py")
exec(open(plot_script_path).read())

total_time = time.time() - overall_start
print(f"\nTotal suite runtime: {total_time/60:.1f} minutes")
