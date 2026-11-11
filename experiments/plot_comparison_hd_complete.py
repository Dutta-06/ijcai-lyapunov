"""
Complete higher-dimensional comparison from the SAVED result arrays
===================================================================

The published figure (results/comparison_higher_dims.png) omitted the n=3
benchmark (where the proposed method has the LARGEST steady-state error of the
five methods except Exact/Deep GP) and the NTK baseline (which beats the proposed
method at n=2 and n=6).  Rule: report all preregistered configurations and all
baselines.  This script regenerates the figure from the raw .npy files with
nothing excluded, and prints the full table.

Caveats carried from the original protocol (see ASSESSMENT.md):
  * one trajectory per cell — no seed-level uncertainty is available;
  * the GP baselines only refresh their prediction every 50–100 steps (a handicap not shared by the
    proposed method); NTK refreshes every 50 steps at n=1 and every step at n=2, 3, 6;
  * every baseline regresses on noise-free samples of the true f(x, t) (equivalent to measuring the state
    derivative exactly), an advantage not shared by the proposed method, which uses only the tracking error;
  * the Wilcoxon p-values of the original scripts were computed on ~12 000
    autocorrelated time samples (lag-1 autocorrelation 0.9999) and are not
    valid; they are not reproduced here.
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), '..', 'results')
DIMS = {'n=1 (1-D plant, Lorenz ref)': ROOT, 'n=2 (Van der Pol)': os.path.join(ROOT, 'n2'),
        'n=3 (Lorenz plant)': os.path.join(ROOT, 'n3'), 'n=6 (coupled Duffing)': os.path.join(ROOT, 'n6')}
METHODS = [('Proposed (Proj)', 'proposed', '#4CAF50'), ('Exact GP', 'exact_gp', '#2196F3'),
           ('Sparse GP', 'sparse_gp', '#FF9800'), ('Deep GP', 'deep_gp', '#9C27B0'), ('NTK', 'ntk', '#F44336')]


def load(rdir, key, what):
    cands = [f"{key}_{what}.npy", f"{key}_fixed_{what}.npy"]
    for c in cands:
        p = os.path.join(rdir, c)
        if os.path.exists(p):
            return np.load(p)
    return None


def main():
    fig, axes = plt.subplots(2, len(DIMS), figsize=(22, 9))
    lines = ["# Higher-dimensional comparison — complete (nothing excluded)", "",
             "Single trajectory per cell; steady state = last 40% of 30 s. "
             "Step time = mean wall-clock of the method's update call (GP: per update, every 50–100 steps; "
             "NTK: per update, every 50 steps at n=1 and every step at n=2, 3, 6; proposed: per step). "
             "Every baseline regresses on noise-free samples of the true f(x, t); the proposed method uses only the "
             "tracking error. No valid significance test is available (one run per cell).", "",
             "| benchmark | method | mean SS ‖e‖ | median SS | max SS | rank (mean SS) | mean update time (ms) |",
             "|---|---|---|---|---|---|---|"]
    for col, (label, rdir) in enumerate(DIMS.items()):
        rows = []
        for name, key, color in METHODS:
            e = load(rdir, key, 'error'); ct = load(rdir, key, 'comptime')
            if e is None:
                continue
            k = int(0.6 * len(e)); ss = e[k:]
            rows.append((name, key, color, ss.mean(), np.median(ss), ss.max(), ct.mean() * 1e3 if ct is not None else np.nan, e))
        order = sorted(range(len(rows)), key=lambda i: rows[i][3])
        rank = {rows[i][0]: r + 1 for r, i in enumerate(order)}
        for name, key, color, mean_ss, med, mx, ct_ms, e in rows:
            lines.append(f"| {label} | {name} | {mean_ss:.4f} | {med:.4f} | {mx:.3f} | {rank[name]}/{len(rows)} | {ct_ms:.2f} |")
        ax = axes[0, col]
        ax.bar([r[0] for r in rows], [r[3] for r in rows], color=[r[2] for r in rows], alpha=0.85)
        ax.set_yscale('log'); ax.set_title(label, fontsize=11); ax.tick_params(axis='x', rotation=30)
        if col == 0:
            ax.set_ylabel('mean steady-state ‖e‖ (single run)')
        ax2 = axes[1, col]
        for name, key, color, mean_ss, med, mx, ct_ms, e in rows:
            w = 300
            sm = np.convolve(e, np.ones(w) / w, mode='valid')
            ax2.semilogy(np.arange(len(sm)) * 1e-3, sm, color=color, label=name, lw=1.2)
        ax2.set_xlabel('Time (s)'); ax2.grid(True, alpha=0.3)
        if col == 0:
            ax2.set_ylabel('‖e(t)‖ (300-step moving avg)')
        if col == len(DIMS) - 1:
            ax2.legend(fontsize=8)
    fig.suptitle('All benchmarks, all baselines (regenerated from the saved raw arrays)')
    plt.tight_layout()
    out = os.path.join(ROOT, 'comparison_higher_dims_complete.png')
    plt.savefig(out, dpi=160)
    with open(os.path.join(ROOT, 'comparison_complete_table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines)); print("saved", out)


if __name__ == "__main__":
    main()
