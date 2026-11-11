"""
Factorial comparison: robustness mechanism × network geometry
=============================================================

The manuscript's central empirical claim is that projection removes the
steady-state bias of leaky modifications.  The original Tier-1 figure
confounded that claim with a second, different effect (adaptive centres vs a
fixed grid): the baselines there were fixed-grid, weights-only networks whose
grid did not cover the operating region, so they failed representationally
(error ≈ 15, larger than open loop), which says nothing about bias.

This experiment separates the two factors on the identical Lorenz tracking
task used by the original figure:

    mechanism  ∈ {projection (proposed), σ-modification, e-modification}
    geometry   ∈ {fixed grid (W only), adaptive W,c, adaptive W,c,σ}

All nine cells share the same plant, reference, initial geometry (a 3×3×3
grid that covers the bounding box of the reference), learning rates, and the
same K random seeds (random x₀ and centre jitter).  Metrics are per-run scalars;
statistics are computed over seeds (not over autocorrelated time samples).

Outputs (results/factorial/): raw.npz, table.md, figure.png.
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import json
import time
import numpy as np
import multiprocessing as mp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import LorenzPlant
from src.core.controller import AdaptiveController
from src.core.integrator import AdaptiveRK4Integrator, SimulationState
from src.baselines.leaky_geometry import LeakyGeometryController
from src.core.projection import min_center_separation

T_END, DT, M, N, KE = 30.0, 1e-3, 27, 3, 10.0
GW, GC, GS = 50.0, 10.0, 2.0
SIGMA_MOD = 0.05
MECHS = ['projection', 'sigma-mod', 'e-mod']
GEOMS = ['fixed (W)', 'adaptive (W,c)', 'adaptive (W,c,σ)']
OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'factorial')


def make_plant():
    return LorenzPlant(sigma=10.0, rho=28.0, beta=8.0 / 3.0, ref_sigma=10.5, ref_rho=28.5, ref_beta=3.0,
                       ref_duration=50.0, ref_dt=1e-3, amplitude_scale=0.1)


def covering_grid(plant):
    ts = np.arange(0.0, T_END, 0.01)
    xr = np.array([plant.x_ref(t) for t in ts])
    lo, hi = xr.min(axis=0), xr.max(axis=0)
    axes = [np.linspace(lo[k], hi[k], 3) for k in range(N)]
    c0 = np.array(np.meshgrid(*axes, indexing='ij')).reshape(N, -1).T
    spacing = float(np.mean((hi - lo) / 2.0))
    return c0, spacing, lo, hi


def run_cell(args):
    mech, geom, seed = args
    plant = make_plant()
    c0, spacing, lo, hi = covering_grid(plant)
    rng = np.random.default_rng(1000 + seed)
    c0 = c0 + rng.normal(0, 0.1 * spacing, size=c0.shape)
    sigma0 = np.ones(M) * spacing
    x0 = plant.x_ref(0.0) + rng.uniform(-0.5, 0.5, size=N)
    W0 = np.zeros((M, N))
    gc = GC if geom != 'fixed (W)' else 0.0
    gs = GS if geom == 'adaptive (W,c,σ)' else 0.0
    sigma_min = 0.25 * spacing
    if mech == 'projection':
        ctrl = AdaptiveController(m=M, n=N, ke=KE, gamma_W=GW, gamma_c=gc, gamma_sigma=gs,
                                  W_max=50.0, c_max=20.0, sigma_min=sigma_min, sigma_max=5.0 * spacing,
                                  delta_min=0.1)
    else:
        ctrl = LeakyGeometryController(m=M, n=N, ke=KE, gamma_W=GW, gamma_c=gc, gamma_sigma=gs,
                                       sigma_mod=SIGMA_MOD, mode='sigma' if mech == 'sigma-mod' else 'e',
                                       c0=c0, sigma0=sigma0, sigma_min=sigma_min)
    t0 = time.time()
    res = AdaptiveRK4Integrator(plant, ctrl, dt=DT).simulate(SimulationState(x0, W0, c0, sigma0), 0.0, T_END,
                                                             record_interval=10)
    d = res.to_arrays()
    e = d['e_norm']
    k = int(0.6 * len(e))
    diverged = bool(not np.all(np.isfinite(e)) or np.nanmax(e) > 100.0)
    return dict(mech=mech, geom=geom, seed=seed, diverged=diverged,
                mean_ss=float(np.nanmean(e[k:])), max_ss=float(np.nanmax(e[k:])),
                median_ss=float(np.nanmedian(e[k:])), max_all=float(np.nanmax(e)),
                max_W=float(np.nanmax(np.linalg.norm(d['W'], axis=(1, 2)))),
                min_sep=float(min(min_center_separation(c) for c in d['c'][::20])),
                wall=time.time() - t0, t=d['t'], e=e)


def ci95(v):
    v = np.asarray(v, float)
    if len(v) < 2:
        return 0.0
    return 1.96 * v.std(ddof=1) / np.sqrt(len(v))


def main(n_seeds=5, procs=None):
    os.makedirs(OUT, exist_ok=True)
    jobs = [(mech, geom, s) for geom in GEOMS for mech in MECHS for s in range(n_seeds)]
    procs = procs or max(1, mp.cpu_count() - 1)
    t0 = time.time()
    with mp.Pool(procs) as pool:
        rows = pool.map(run_cell, jobs)
    print(f"{len(rows)} runs in {time.time() - t0:.0f}s")
    # raw arrays
    np.savez(os.path.join(OUT, 'raw.npz'),
             meta=json.dumps(dict(T_END=T_END, DT=DT, M=M, N=N, KE=KE, GW=GW, GC=GC, GS=GS,
                                  SIGMA_MOD=SIGMA_MOD, n_seeds=n_seeds, mechs=MECHS, geoms=GEOMS)),
             **{f"e__{r['mech']}__{r['geom']}__{r['seed']}": r['e'] for r in rows},
             t=rows[0]['t'])
    scal = [{k: v for k, v in r.items() if k not in ('t', 'e')} for r in rows]
    with open(os.path.join(OUT, 'per_run.json'), 'w') as f:
        json.dump(scal, f, indent=1)
    # table
    lines = ["# Factorial: robustness mechanism × geometry (Lorenz tracking, m=27, k_e=10)", "",
             f"Seeds: {n_seeds}; metric = steady-state ‖e‖ over the last 40% of a {T_END:.0f}s run; "
             "mean ± 95% CI over seeds (max over window in brackets); 'div' = diverged runs.", "",
             "| geometry | mechanism | mean SS ‖e‖ | median SS | max SS | max ‖W‖_F | min sep | div |",
             "|---|---|---|---|---|---|---|---|"]
    summary = {}
    for geom in GEOMS:
        for mech in MECHS:
            rr = [r for r in rows if r['mech'] == mech and r['geom'] == geom]
            ms = np.array([r['mean_ss'] for r in rr]); md = np.array([r['median_ss'] for r in rr])
            mx = np.array([r['max_ss'] for r in rr])
            summary[(geom, mech)] = ms
            lines.append(f"| {geom} | {mech} | {ms.mean():.4f} ± {ci95(ms):.4f} | {md.mean():.4f} | "
                         f"{mx.mean():.3f} | {np.mean([r['max_W'] for r in rr]):.1f} | "
                         f"{np.mean([r['min_sep'] for r in rr]):.3f} | {sum(r['diverged'] for r in rr)} |")
    lines += ["", "## Paired (same seed) differences, projection minus baseline, mean SS ‖e‖", "",
              "| geometry | baseline | mean diff | 95% CI | projection better in |", "|---|---|---|---|---|"]
    for geom in GEOMS:
        for mech in MECHS[1:]:
            diff = summary[(geom, 'projection')] - summary[(geom, mech)]
            lines.append(f"| {geom} | {mech} | {diff.mean():+.4f} | ±{ci95(diff):.4f} | "
                         f"{int(np.sum(diff < 0))}/{len(diff)} seeds |")
    lines += ["", "## Geometry effect within the projection mechanism (paired, mean SS ‖e‖)", "",
              "| comparison | mean diff | 95% CI | adaptive better in |", "|---|---|---|---|"]
    for a, b in [(GEOMS[1], GEOMS[0]), (GEOMS[2], GEOMS[1])]:
        diff = summary[(a, 'projection')] - summary[(b, 'projection')]
        lines.append(f"| {a} − {b} | {diff.mean():+.4f} | ±{ci95(diff):.4f} | {int(np.sum(diff < 0))}/{len(diff)} |")
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    # figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    colors = {'projection': '#1f77b4', 'sigma-mod': '#d62728', 'e-mod': '#ff7f0e'}
    t = rows[0]['t']
    for ax, geom in zip(axs, GEOMS):
        for mech in MECHS:
            E = np.array([r['e'] for r in rows if r['mech'] == mech and r['geom'] == geom])
            med = np.nanmedian(E, axis=0); lo = np.nanmin(E, axis=0); hi = np.nanmax(E, axis=0)
            ax.plot(t, med, color=colors[mech], label=f"{mech} (median of {n_seeds})")
            ax.fill_between(t, lo, hi, color=colors[mech], alpha=0.15)
        ax.set_yscale('log'); ax.set_title(geom); ax.set_xlabel('Time (s)'); ax.grid(True, alpha=0.3)
    axs[0].set_ylabel('‖e(t)‖'); axs[0].legend(fontsize=8)
    fig.suptitle('Mechanism × geometry on Lorenz tracking (same seeds, same covering grid, m=27)')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'figure.png'), dpi=160)
    print("saved", os.path.join(OUT, 'figure.png'))


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    main(n_seeds=n)
