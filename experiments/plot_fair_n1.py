"""Aggregate results/fair_n1/<method>_seed<k>.npz into a seed-level table and figure."""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'fair_n1')
METHODS = [('proposed', 'Proposed (Proj)', '#4CAF50'), ('exact_gp', 'Exact GP', '#2196F3'),
           ('sparse_gp', 'Sparse GP (9×3 grid)', '#FF9800'), ('deep_gp', 'Deep GP', '#9C27B0'), ('ntk', 'NTK', '#F44336')]


def main():
    data = {}
    for key, label, color in METHODS:
        runs = []
        for f in sorted(glob.glob(os.path.join(OUT, f"{key}_seed*.npz"))):
            r = np.load(f); e = r['e']; k = int(0.6 * len(e))
            runs.append(dict(seed=int(f.split('seed')[-1].split('.')[0]), mean_ss=e[k:].mean(), max_ss=e[k:].max(),
                             step_ms=np.median(r['step_time']) * 1e3, step_mean_ms=r['step_time'].mean() * 1e3, e=e, t=r['t']))
        data[key] = runs
    seeds = sorted(set(r['seed'] for r in data['proposed']))
    lines = ["# Fair 1-D kernel comparison (per-step prediction for every method, covering inducing grid, seed-level statistics)", "",
             f"Seeds: {len(seeds)} (different reference / plant initial conditions). Steady state = last 40% of 30 s. "
             "Step time = median full online step (prediction + amortised training) after warm-up, same process type, 2 BLAS threads.", "",
             "| method | mean SS ‖e‖ (mean ± 95% CI over seeds) | max SS ‖e‖ (mean) | median step (ms) | mean step (ms) | seeds where proposed has lower mean SS |",
             "|---|---|---|---|---|---|"]
    prop = {r['seed']: r['mean_ss'] for r in data['proposed']}
    for key, label, color in METHODS:
        rr = data[key]
        if not rr:
            continue
        ms = np.array([r['mean_ss'] for r in rr])
        ci = 1.96 * ms.std(ddof=1) / np.sqrt(len(ms)) if len(ms) > 1 else 0.0
        wins = sum(prop[r['seed']] < r['mean_ss'] for r in rr if r['seed'] in prop)
        lines.append(f"| {label} | {ms.mean():.4f} ± {ci:.4f} | {np.mean([r['max_ss'] for r in rr]):.3f} | "
                     f"{np.mean([r['step_ms'] for r in rr]):.3f} | {np.mean([r['step_mean_ms'] for r in rr]):.3f} | "
                     f"{'—' if key == 'proposed' else f'{wins}/{len(rr)}'} |")
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    for key, label, color in METHODS:
        rr = data[key]
        if not rr:
            continue
        E = np.array([r['e'] for r in rr]); w = 300
        med = np.median(E, axis=0); sm = np.convolve(med, np.ones(w) / w, mode='valid')
        axs[0].semilogy(rr[0]['t'][:len(sm)], sm, color=color, label=label, lw=1.5)
        axs[1].bar(label, np.mean([r['mean_ss'] for r in rr]), color=color, alpha=0.85,
                   yerr=1.96 * np.std([r['mean_ss'] for r in rr], ddof=1) / np.sqrt(len(rr)) if len(rr) > 1 else 0)
    axs[0].set_xlabel('Time (s)'); axs[0].set_ylabel('‖e(t)‖ (median over seeds, 300-step avg)'); axs[0].legend(fontsize=8); axs[0].grid(True, alpha=0.3)
    axs[1].set_ylabel('mean steady-state ‖e‖ (± 95% CI over seeds)'); axs[1].tick_params(axis='x', rotation=20)
    fig.suptitle('1-D benchmark, fair protocol'); plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'figure.png'), dpi=160); print("saved figure")


if __name__ == "__main__":
    main()
