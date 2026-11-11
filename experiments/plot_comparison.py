"""
Composite Comparison Plot: All Kernel Baselines vs Proposed
=============================================================

Generates the single comparison figure combining all five methods
on the same 1D Lorenz tracking benchmark, including computation time.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats

results = {
    'Proposed (Proj)':       ('results/proposed_error.npy',         'results/proposed_time.npy',         'results/proposed_comptime.npy',         '#4CAF50'),
    'Exact GP':              ('results/exact_gp_error.npy',         'results/exact_gp_time.npy',         'results/exact_gp_comptime.npy',         '#2196F3'),
    'Sparse GP (Fixed)':     ('results/sparse_gp_fixed_error.npy',  'results/sparse_gp_fixed_time.npy',  'results/sparse_gp_fixed_comptime.npy',  '#FF9800'),
    'Deep GP (2-layer)':     ('results/deep_gp_error.npy',          'results/deep_gp_time.npy',          'results/deep_gp_comptime.npy',          '#9C27B0'),
}

fig, axes = plt.subplots(1, 3, figsize=(22, 6))

# Panel 1: Tracking error over time
ax1 = axes[0]
for label, (e_path, t_path, ct_path, color) in results.items():
    if os.path.exists(e_path) and os.path.exists(t_path):
        e = np.load(e_path)
        t = np.load(t_path)
        # Smooth with moving average for clarity
        window = 500
        if len(e) > window:
            e_smooth = np.convolve(e, np.ones(window)/window, mode='valid')
            t_smooth = t[:len(e_smooth)]
        else:
            e_smooth = e
            t_smooth = t
        lw = 2.5 if 'Proposed' in label else 1.5
        ax1.semilogy(t_smooth, e_smooth, label=label, color=color, linewidth=lw)

ax1.set_xlabel('Time (s)', fontsize=12)
ax1.set_ylabel('Tracking Error ||e(t)|| (log scale)', fontsize=12)
ax1.set_title('Chaotic Reference Tracking -- All Kernel Methods', fontsize=13)
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Panel 2: Steady-state error distribution (box plot)
ax2 = axes[1]
labels_list = []
data_list = []
colors_list = []
for label, (e_path, t_path, ct_path, color) in results.items():
    if os.path.exists(e_path):
        e = np.load(e_path)
        ss_start = int(0.6 * len(e))
        data_list.append(e[ss_start:])
        labels_list.append(label)
        colors_list.append(color)

if data_list:
    bp = ax2.boxplot(data_list, tick_labels=labels_list, patch_artist=True,
                     showfliers=False, medianprops=dict(color='black', linewidth=2))
    for patch, color in zip(bp['boxes'], colors_list):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

ax2.set_ylabel('Steady-State ||e||', fontsize=12)
ax2.set_title('Steady-State Error Distribution', fontsize=13)
ax2.set_yscale('log')
ax2.grid(True, alpha=0.3, axis='y')
ax2.tick_params(axis='x', rotation=15)

# Panel 3: Computation time per update step over trajectory length
ax3 = axes[2]
for label, (e_path, t_path, ct_path, color) in results.items():
    if ct_path and os.path.exists(ct_path) and os.path.exists(t_path):
        ct = np.load(ct_path)
        t = np.load(t_path)
        # Note: some scripts record ct every N steps. We'll map ct length across t.
        # Create an approximate x-axis for ct array length
        t_ct = np.linspace(0, t[-1], len(ct))
        lw = 2.5 if 'Proposed' in label else 1.5
        
        # Plot raw computation time or smoothed
        window = min(50, len(ct)//10)
        if window > 0:
            ct_smooth = np.convolve(ct, np.ones(window)/window, mode='valid')
            t_ct_smooth = t_ct[:len(ct_smooth)]
            ax3.semilogy(t_ct_smooth, ct_smooth, label=label, color=color, linewidth=lw)
        else:
            ax3.semilogy(t_ct, ct, label=label, color=color, linewidth=lw)
    else:
        print(f"WARNING: Missing comptime data for {label}: {ct_path}")

ax3.set_xlabel('Time (s)', fontsize=12)
ax3.set_ylabel('Computation Time per Step (s) (log scale)', fontsize=12)
ax3.set_title('Computational Cost Over Time', fontsize=13)
ax3.grid(True, alpha=0.3)
ax3.legend(fontsize=10)

plt.tight_layout()
os.makedirs('results', exist_ok=True)
plt.savefig('results/comparison_all_baselines.png', dpi=200, bbox_inches='tight')
print("Comparison plot saved to results/comparison_all_baselines.png")

# Print summary table with Wilcoxon p-values
print("\n--- Summary ---")
print(f"{'Method':<20} {'Mean SS Error':<18} {'Wilcoxon vs Proposed':<25} {'Mean Step Time'}")
print("-" * 80)

# Load proposed data for statistical test
prop_e_path = results['Proposed (Proj)'][0]
if os.path.exists(prop_e_path):
    prop_e = np.load(prop_e_path)
    prop_ss = prop_e[int(0.6*len(prop_e)):]
else:
    prop_ss = None

for label, (e_path, t_path, ct_path, color) in results.items():
    if os.path.exists(e_path):
        e = np.load(e_path)
        ss = e[int(0.6*len(e)):]
        mean_ss = ss.mean()
        
        # Calculate p-value
        p_val_str = "---"
        if prop_ss is not None and label != 'Proposed (Proj)':
            # Ensure lengths match for paired test
            min_len = min(len(prop_ss), len(ss))
            stat, p = stats.wilcoxon(prop_ss[:min_len], ss[:min_len])
            if p < 0.001:
                p_val_str = "p < 0.001"
            elif p > 0.05:
                p_val_str = "p > 0.05 (n.s.)"
            else:
                p_val_str = f"p = {p:.4f}"
                
        # Calculate mean step time in ms
        ct_str = "N/A"
        if ct_path and os.path.exists(ct_path):
            ct = np.load(ct_path)
            mean_ct_ms = ct.mean() * 1000
            ct_str = f"{mean_ct_ms:.2f} ms"
            
        print(f"{label:<20} {mean_ss:<18.4f} {p_val_str:<25} {ct_str}")
