"""
Master Runner: Higher-Dimensional Benchmark Suite
=================================================

Runs all 15 simulations (n=2, n=3, n=6) sequentially and 
generates the final higher-dimensional comparison plot.

Note: Deep GP for n=6 may take several hours.
"""

import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

print("=" * 70)
print("  HIGHER-DIMENSIONAL KERNEL BASELINE SUITE (HD-1, HD-2, HD-3)")
print("=" * 70)

overall_start = time.time()

# ── HD-1: n=2 (Van der Pol) ──────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  [HD-1] n=2: Van der Pol Oscillator (m=27)")
print("=" * 70)
os.makedirs('results/n2', exist_ok=True)

print("\n[1/5] Proposed (Proj) n=2...")
from sim0_proposed_n2 import run as run_p2; t, e, ct = run_p2(); np.save('results/n2/proposed_error.npy', e); np.save('results/n2/proposed_time.npy', t); np.save('results/n2/proposed_comptime.npy', ct)

print("\n[2/5] NTK Control n=2...")
from sim4_ntk_n2 import run as run_ntk2; t, e, ct = run_ntk2(); np.save('results/n2/ntk_error.npy', e); np.save('results/n2/ntk_time.npy', t); np.save('results/n2/ntk_comptime.npy', ct)

print("\n[3/5] Exact GP n=2...")
from sim1_exact_gp_n2 import run as run_ex2; t, e, ct = run_ex2(); np.save('results/n2/exact_gp_error.npy', e); np.save('results/n2/exact_gp_time.npy', t); np.save('results/n2/exact_gp_comptime.npy', ct)

print("\n[4/5] Sparse GP n=2...")
from sim2_sparse_gp_n2 import run as run_sp2; t, e, ct = run_sp2(); np.save('results/n2/sparse_gp_error.npy', e); np.save('results/n2/sparse_gp_time.npy', t); np.save('results/n2/sparse_gp_comptime.npy', ct)

print("\n[5/5] Deep GP n=2...")
from sim3_deep_gp_n2 import run as run_dp2; t, e, ct = run_dp2(); np.save('results/n2/deep_gp_error.npy', e); np.save('results/n2/deep_gp_time.npy', t); np.save('results/n2/deep_gp_comptime.npy', ct)


# ── HD-2: n=3 (Lorenz 3D) ────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  [HD-2] n=3: 3D Lorenz Plant (m=64)")
print("=" * 70)
os.makedirs('results/n3', exist_ok=True)

print("\n[1/5] Proposed (Proj) n=3...")
from sim0_proposed_n3 import run as run_p3; t, e, ct = run_p3(); np.save('results/n3/proposed_error.npy', e); np.save('results/n3/proposed_time.npy', t); np.save('results/n3/proposed_comptime.npy', ct)

print("\n[2/5] NTK Control n=3...")
from sim4_ntk_n3 import run as run_ntk3; t, e, ct = run_ntk3(); np.save('results/n3/ntk_error.npy', e); np.save('results/n3/ntk_time.npy', t); np.save('results/n3/ntk_comptime.npy', ct)

print("\n[3/5] Exact GP n=3...")
from sim1_exact_gp_n3 import run as run_ex3; t, e, ct = run_ex3(); np.save('results/n3/exact_gp_error.npy', e); np.save('results/n3/exact_gp_time.npy', t); np.save('results/n3/exact_gp_comptime.npy', ct)

print("\n[4/5] Sparse GP n=3...")
from sim2_sparse_gp_n3 import run as run_sp3; t, e, ct = run_sp3(); np.save('results/n3/sparse_gp_error.npy', e); np.save('results/n3/sparse_gp_time.npy', t); np.save('results/n3/sparse_gp_comptime.npy', ct)

print("\n[5/5] Deep GP n=3...")
from sim3_deep_gp_n3 import run as run_dp3; t, e, ct = run_dp3(); np.save('results/n3/deep_gp_error.npy', e); np.save('results/n3/deep_gp_time.npy', t); np.save('results/n3/deep_gp_comptime.npy', ct)


# ── HD-3: n=6 (Coupled Duffing) ──────────────────────────────────────────────
print("\n" + "=" * 70)
print("  [HD-3] n=6: Coupled Duffing Oscillators (m=125)")
print("=" * 70)
os.makedirs('results/n6', exist_ok=True)

print("\n[1/5] Proposed (Proj) n=6...")
from sim0_proposed_n6 import run as run_p6; t, e, ct = run_p6(); np.save('results/n6/proposed_error.npy', e); np.save('results/n6/proposed_time.npy', t); np.save('results/n6/proposed_comptime.npy', ct)

print("\n[2/5] NTK Control n=6...")
from sim4_ntk_n6 import run as run_ntk6; t, e, ct = run_ntk6(); np.save('results/n6/ntk_error.npy', e); np.save('results/n6/ntk_time.npy', t); np.save('results/n6/ntk_comptime.npy', ct)

print("\n[3/5] Exact GP n=6...")
from sim1_exact_gp_n6 import run as run_ex6; t, e, ct = run_ex6(); np.save('results/n6/exact_gp_error.npy', e); np.save('results/n6/exact_gp_time.npy', t); np.save('results/n6/exact_gp_comptime.npy', ct)

print("\n[4/5] Sparse GP n=6...")
from sim2_sparse_gp_n6 import run as run_sp6; t, e, ct = run_sp6(); np.save('results/n6/sparse_gp_error.npy', e); np.save('results/n6/sparse_gp_time.npy', t); np.save('results/n6/sparse_gp_comptime.npy', ct)

print("\n[5/5] Deep GP n=6 (WARNING: Overnight expected)...")
from sim3_deep_gp_n6 import run as run_dp6; t, e, ct = run_dp6(); np.save('results/n6/deep_gp_error.npy', e); np.save('results/n6/deep_gp_time.npy', t); np.save('results/n6/deep_gp_comptime.npy', ct)


# ── Final Plotting ───────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("  Generating Higher-Dimensional Comparison Plot...")
print("=" * 70)
plot_script_path = os.path.join(os.path.dirname(__file__), "plot_comparison_hd.py")
exec(open(plot_script_path).read())

total_time = time.time() - overall_start
print(f"\nTotal suite runtime: {total_time/3600:.2f} hours")
