"""
Theorem H, quantitatively: the speed–bias trade-off of parameter dither
========================================================================

Same plant, reference and network as exp_langevin_global.py (m = 32, σ = 0.002 fixed).  For each
temperature τ (constant), and for an annealed schedule that starts HOT and cools like 1/log t, run the
real closed loop from the 'collapsed' (symmetric, spurious) and the 'spread' start, 2 seeds, T = 1500 s.

Theory (Theorem H, Remark H.2): the averaged mean-field dynamics mixes globally at a rate ≥ 2τρ_τ with
ρ_τ ≥ (σ/τ)exp(−α*/τ) — exponentially small in 1/τ — while the bias of the limit grows with τ.  So:
low τ ⇒ small bias but slow escape from the symmetric configuration (metastability); high τ ⇒ fast
escape but a biased limit.  Reported per τ: final occupation loss for both starts, the escape time
(first t at which the collapsed run's loss falls below 2× the spread run's final loss; ∞ if never),
and W₂ between the final measures of the two starts.

Outputs: results/langevin_tradeoff/{table.md, figure.png}
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from exp_langevin_global import ScalarPlant, init, M, PERIOD, DT, GAINS, KE, w2_phys

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'langevin_tradeoff')
SIGMA = 0.002
T_END = 1500.0
TAUS = [0.0, 0.0005, 0.001, 0.002, 0.004, 0.008]
TAU_HOT = 0.008


def schedule(kind):
    if kind == 'annealed':
        return TAU_HOT, (lambda t: TAU_HOT * np.log(np.e + 10.0) / np.log(np.e + 10.0 + t))
    return float(kind), None


def run(args):
    kind, start, seed = args
    rng = np.random.default_rng(300 + seed)
    plant = ScalarPlant(); net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.15, s_max=2.0)
    tau, sched = schedule(kind)
    ctrl = L.LangevinAdaptiveController(net, ke=KE, G_diag=GAINS, sigma=SIGMA, theta0=np.zeros(3), tau=tau, tau_schedule=sched)
    Th0 = init(start, net, np.random.default_rng(7 + seed))
    res = L.simulate(plant, ctrl, plant.x_ref(0.0), Th0, T=T_END, dt=DT, rng=rng, record_every=2500)
    Ls = np.array([L.occupation_loss(plant, net, Th, PERIOD, 200) for Th in res['Th']])
    phys = np.hstack([*net.physical(res['Th'][-1])[:2], net.physical(res['Th'][-1])[2][:, None]])
    t, e = res['t'], res['e_norm']
    return dict(kind=kind, start=start, seed=seed, t=t, L=Ls, phys=phys,
                rms=float(np.sqrt(np.mean(e[t >= T_END - 10 * PERIOD] ** 2))))



def escape_fit_lines(esc_by_label):
    """Fit t_esc = A + B·log(1/τ) and t_esc = a·exp(b/τ) to the seed-averaged escape times and to the per-seed values
    (τ values at which every seed escapes).  The reading is computed, not asserted."""
    pts = [(float(lab), np.asarray(es, float)) for lab, es in esc_by_label.items()
           if lab not in ('0.0', 'annealed (hot start)') and np.all(np.isfinite(es))]
    if len(pts) < 3:
        return []

    def fits(tau_v, esc_v):
        X1 = np.log(1 / tau_v); c1 = np.polyfit(X1, esc_v, 1); r1 = esc_v - np.polyval(c1, X1)
        X2 = 1 / tau_v; c2 = np.polyfit(X2, np.log(esc_v), 1); r2 = esc_v - np.exp(np.polyval(c2, X2))
        return c1, c2, float(np.sqrt(np.mean((r1 / esc_v) ** 2))), float(np.sqrt(np.mean((r2 / esc_v) ** 2)))
    tau_m = np.array([p[0] for p in pts]); esc_m = np.array([p[1].mean() for p in pts])
    tau_s = np.concatenate([np.full(len(p[1]), p[0]) for p in pts]); esc_s = np.concatenate([p[1] for p in pts])
    c1, c2, rm1, rm2 = fits(tau_m, esc_m)
    _, _, rs1, rs2 = fits(tau_s, esc_s)
    out = ["", "## Escape-time scaling (τ values at which every seed escapes)", "",
           "| model | fitted form (seed-averaged) | RMS relative residual, seed-averaged | RMS relative residual, per seed |",
           "|---|---|---|---|",
           f"| t_esc = A + B·log(1/τ) | {c1[1]:.0f} + {c1[0]:.0f}·log(1/τ) | {rm1:.3f} | {rs1:.3f} |",
           f"| t_esc = a·exp(b/τ) (Arrhenius) | {np.exp(c2[1]):.0f}·exp({c2[0]:.2e}/τ) | {rm2:.3f} | {rs2:.3f} |", ""]
    if max(rs1, rs2) > 2 * min(rm1, rm2):
        out.append("Reading: the seed-averaged fits differ, but on the per-seed values both laws fit poorly, so these data do not "
                   "identify the escape mechanism (PARADIGM.md Remark H.2).")
    else:
        out.append("Reading: the comparison is consistent across seed-averaged and per-seed values.")
    return out


def refit_from_raw():
    """Regenerate only the escape-time section of table.md from raw.npz (the loss paths); the rest of the table is kept."""
    z = np.load(os.path.join(OUT, 'raw.npz'))
    t = z['t']
    kinds = sorted({k.split('__')[1] for k in z.files if k.startswith('L__')}, key=lambda k: (k == 'annealed', k))
    seeds = sorted({int(k.split('__')[3]) for k in z.files if k.startswith('L__')})
    esc_by = {}
    for k in kinds:
        es = []
        for sd in seeds:
            Lc, Ls = z[f'L__{k}__collapsed__{sd}'], z[f'L__{k}__spread__{sd}']
            hit = np.where(Lc <= 2 * Ls[-1])[0]
            es.append(t[hit[0]] if len(hit) else np.inf)
        esc_by['annealed (hot start)' if k == 'annealed' else k] = es
    path = os.path.join(OUT, 'table.md')
    head = open(path, encoding='utf-8').read().split("\n## Escape-time scaling")[0].rstrip("\n")
    lines = head.split("\n") + escape_fit_lines(esc_by)
    with open(path, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

def main(n_seeds=2):
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    kinds = [str(t) for t in TAUS] + ['annealed']
    jobs = [(k, s0, sd) for k in kinds for s0 in ('collapsed', 'spread') for sd in range(n_seeds)]
    with mp.Pool(min(len(jobs), 12)) as pool:
        rows = pool.map(run, jobs)
    lines = [f"# Speed–bias trade-off of parameter dither (m = 32, σ = {SIGMA}, k_e = {KE:g}, T = {T_END:.0f} s, {n_seeds} seeds)", "",
             "| τ | L̄ final, collapsed | L̄ final, spread | escape time of the collapsed start (s) | W₂(final collapsed, final spread) | tracking RMS, spread |",
             "|---|---|---|---|---|---|"]
    summary = []
    for k in kinds:
        Lc, Ls, esc, w2s, rms = [], [], [], [], []
        for sd in range(n_seeds):
            rc = [r for r in rows if r['kind'] == k and r['start'] == 'collapsed' and r['seed'] == sd][0]
            rs = [r for r in rows if r['kind'] == k and r['start'] == 'spread' and r['seed'] == sd][0]
            Lc.append(rc['L'][-1]); Ls.append(rs['L'][-1]); rms.append(rs['rms'])
            thr = 2 * rs['L'][-1]
            hit = np.where(rc['L'] <= thr)[0]
            esc.append(rc['t'][hit[0]] if len(hit) else np.inf)
            w2s.append(w2_phys(rc['phys'], rs['phys']))
        label = 'annealed (hot start)' if k == 'annealed' else k
        esc_str = ', '.join('∞' if np.isinf(x) else f"{x:.0f}" for x in esc)
        lines.append(f"| {label} | {np.mean(Lc):.2e} | {np.mean(Ls):.2e} | {esc_str} | {np.mean(w2s):.2f} | {np.mean(rms):.4f} |")
        summary.append((label, np.mean(Lc), np.mean(Ls), esc))
    lines += ["", "Escape time = first sample at which the collapsed run's loss ≤ 2 × the spread run's final loss (same seed); "
                  "∞ = not within the horizon.  Per-seed values listed."]
    lines += escape_fit_lines({lab: es for lab, _, _, es in summary})
    np.savez(os.path.join(OUT, 'raw.npz'), **{f"L__{r['kind']}__{r['start']}__{r['seed']}": r['L'] for r in rows}, t=rows[0]['t'])
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(kinds), figsize=(4 * len(kinds), 4), sharey=True)
    for ax, k in zip(axs, kinds):
        for r in [r for r in rows if r['kind'] == k]:
            ax.semilogy(r['t'], r['L'], 'r' if r['start'] == 'collapsed' else 'b', alpha=.7)
        ax.set_title('annealed' if k == 'annealed' else f'τ = {k}'); ax.set_xlabel('t'); ax.grid(alpha=.3)
    axs[0].set_ylabel('L̄(μ_t)  (red: collapsed, blue: spread)')
    plt.tight_layout(); plt.savefig(os.path.join(OUT, 'figure.png'), dpi=140); print('saved figure')


if __name__ == "__main__":
    refit_from_raw() if os.environ.get('TR_REFIT') == '1' else main()
