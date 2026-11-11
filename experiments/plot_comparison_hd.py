# plot_comparison_hd.py
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import wilcoxon

DIMS = {
    'n=1 (Lorenz ref)': 'results',
    'n=2 (Van der Pol)': 'results/n2',
    'n=6 (Duffing)': 'results/n6',
}

METHODS = {
    'Proposed (Proj)': ('proposed_error.npy', 'proposed_comptime.npy', '#4CAF50'),
    'Exact GP':        ('exact_gp_error.npy', 'exact_gp_comptime.npy', '#2196F3'),
    'Sparse GP':       ('sparse_gp_fixed_error.npy' if 'results' else 'sparse_gp_error.npy', 'sparse_gp_fixed_comptime.npy' if 'results' else 'sparse_gp_comptime.npy', '#FF9800'),
    'Deep GP':         ('deep_gp_error.npy',  'deep_gp_comptime.npy',  '#9C27B0'),
    'NTK Control':     ('ntk_error.npy',      'ntk_comptime.npy',      '#F44336'),
}
# wait, the file names for n=1 Sparse GP is sparse_gp_fixed_error.npy 
# and for others it is sparse_gp_error.npy. I need to dynamically select the name.

def get_filenames(rdir, mname):
    if mname == 'Proposed (Proj)': return 'proposed_error.npy', 'proposed_comptime.npy', '#4CAF50'
    if mname == 'Exact GP':        return 'exact_gp_error.npy', 'exact_gp_comptime.npy', '#2196F3'
    if mname == 'Deep GP':         return 'deep_gp_error.npy', 'deep_gp_comptime.npy', '#9C27B0'
    if mname == 'NTK Control':     return 'ntk_error.npy', 'ntk_comptime.npy', '#F44336'
    if mname == 'Sparse GP':
        if rdir == 'results':
            return 'sparse_gp_fixed_error.npy', 'sparse_gp_fixed_comptime.npy', '#FF9800'
        else:
            return 'sparse_gp_error.npy', 'sparse_gp_comptime.npy', '#FF9800'

def load(path):
    return np.load(path) if os.path.exists(path) else None

def ss_mean(e, frac=0.6):
    return e[int(frac*len(e)):].mean() if e is not None else None

def wilcoxon_str(a, b, frac=0.6):
    if a is None or b is None: return 'N/A'
    n_len = min(len(a), len(b))
    ss = int(frac * n_len)
    try:
        _, p = wilcoxon(a[ss:n_len], b[ss:n_len], alternative='two-sided')
        return 'p > 0.05 (n.s.)' if p > 0.05 else ('p < 0.05' if p > 0.001 else 'p < 0.001')
    except:
        return 'N/A'

# ── Figure: 2 rows x 4 cols (row=metric, col=dimension) ──────────────────────
fig, axes = plt.subplots(2, len(DIMS), figsize=(20, 10))
fig.suptitle('Higher-Dimensional Validation: All Kernel Methods', fontsize=14)

print(f"\n{'Dim':<25} {'Method':<20} {'Mean SS |e|':<15} {'Wilcoxon vs Prop':<22} {'Mean Step (ms)':<15}")
print("-" * 97)

for col, (dim_label, rdir) in enumerate(DIMS.items()):
    e_prop = load(f'{rdir}/proposed_error.npy')
    ct_prop = load(f'{rdir}/proposed_comptime.npy')

    ax_err = axes[0, col]
    ax_cost = axes[1, col]
    ax_err.set_title(dim_label, fontsize=11)

    for mname in ['Proposed (Proj)', 'Exact GP', 'Sparse GP', 'Deep GP']:
        efile, ctfile, color = get_filenames(rdir, mname)
        e = load(f'{rdir}/{efile}')
        ct = load(f'{rdir}/{ctfile}') if ctfile else None
        if e is None:
            continue

        # Panel 1: steady-state error bar
        ss = ss_mean(e)
        ax_err.bar(mname, ss, color=color, alpha=0.8)
        ax_err.set_yscale('log')
        ax_err.set_ylabel('Mean SS ||e||' if col == 0 else '')
        ax_err.tick_params(axis='x', rotation=30)

        # Panel 2: computation time
        if ct is not None:
            ax_cost.plot(np.arange(len(ct)) * 50 * 0.001,
                        ct * 1000, color=color, label=mname, linewidth=1.2)
        ax_cost.set_xlabel('Time (s)')
        ax_cost.set_ylabel('Step time (ms)' if col == 0 else '')
        ax_cost.set_yscale('log')

        # Summary table
        w_str = wilcoxon_str(e_prop, e) if mname != 'Proposed (Proj)' else '---'
        ct_str = f"{ct.mean()*1000:.2f}" if ct is not None else 'N/A'
        print(f"{dim_label:<25} {mname:<20} {ss:<15.4f} {w_str:<22} {ct_str:<15}")

    if col == len(DIMS) - 1:
        axes[1, col].legend(fontsize=8, loc='upper right')

plt.tight_layout()
plt.savefig('results/comparison_higher_dims.png', dpi=200, bbox_inches='tight')
print("\nSaved: results/comparison_higher_dims.png")
