"""
Theorem C in numbers: OT σ-modification vs labelled σ-modification vs projection
=================================================================================

Oracle plant (exactly representable, lifted form).  Base configuration Q; prior μ⁰ = Q + noise
(the a-priori guess); ideal μ* = Q.  Two starting configurations with IDENTICAL geometry and
IDENTICAL represented function, differing only in particle labels:

  (a) 'aligned start'    P0 = μ⁰                 (labels agree with the prior's labels)
  (b) 'relabelled start' P0 = μ⁰ with its particles cyclically relabelled

Case (b) is what happens after kernel migration, warm starts or any relabelling of a network:
the a-priori grid's labels no longer correspond to the network's.  Label-free methods
(projection; OT σ-modification) cannot see the difference.  The labelled σ-modification pulls
particle i toward prior particle i, which in (b) is a different kernel: its bias becomes the
labelled distance b_lab(b) ≫ b.

Mechanisms: projection only; labelled σ-mod (+projection); OT σ-mod (+projection), σ_m = 0.3.
Metrics: steady-state RMS ‖e‖ (second half), final OT error W(μ_T, μ*), final labelled error,
bias constants b = W(μ*, μ⁰) (label-free) and b_lab = labelled distance between the ideal that each
particle converges to and the prior particle its leak pulls it to.  5 seeds, same x₀ per seed.

Outputs: results/ot_sigma_mod/{table.md, figure.png}
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import transport as T

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'ot_sigma_mod')
BOX = dict(a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)
GAINS = dict(ke=3.0, g_a=2.0, g_c=1.0, g_s=0.5)
MECHS = ['projection', 'labelled sigma-mod', 'OT sigma-mod']
SIGMA_M = 0.3


def run(mech, start, seed, m=6, n=2, T_end=40.0):
    rng = np.random.default_rng(seed)
    base = T.TransportController(m, n, **GAINS, **BOX)
    Q = T.sample_lifted_ideal(m, n, base, rng, c_spread=1.5, min_sep=0.6)
    prior = base.hard_clip(Q + rng.normal(0, 0.15, size=Q.shape))
    x0 = Q[0, n:2 * n] * 0 + rng.uniform(-0.5, 0.5, size=n)
    perm = np.roll(np.arange(m), 1) if start == 'relabelled' else np.arange(m)
    P0 = prior[perm]                                            # same geometry & function, different labels
    plant = T.LiftedRepresentablePlant(Q)
    leak = None if mech == 'projection' else ('ot' if mech.startswith('OT') else 'labelled')
    ctrl = T.TransportController(m, n, leak=leak, sigma_m=0.0 if leak is None else SIGMA_M, prior=prior, **GAINS, **BOX)
    res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0) + x0, P0, T=T_end, dt=1e-3, record_every=20)
    k = len(res['t']) // 2
    b = T.w2(Q, prior, ctrl.Ginv_diag)
    # labelled bias: particle i starts at prior[perm[i]], its natural ideal is Q[perm[i]], its labelled leak target is prior[i]
    b_lab = float(np.sqrt(np.einsum('ik,k,ik->', Q[perm] - prior, ctrl.Ginv_diag, Q[perm] - prior) / m))
    lab_end = float(np.sqrt(np.einsum('ik,k,ik->', res['P'][-1] - Q[perm], ctrl.Ginv_diag, res['P'][-1] - Q[perm]) / m))
    return dict(mech=mech, start=start, seed=seed, rms=float(np.sqrt(np.mean(res['e_norm'][k:] ** 2))),
                W_end=float(ctrl.ot_distance(res['P'][-1], Q)[1] ** 0.5), lab_end=lab_end, b=b, b_lab=b_lab,
                t=res['t'], e=res['e_norm'])


def _job(args):
    return run(*args)


def main(n_seeds=5):
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    jobs = [(mech, st, s) for st in ('aligned', 'relabelled') for mech in MECHS for s in range(n_seeds)]
    with mp.Pool(min(len(jobs), 10)) as pool:
        rows = pool.map(_job, jobs)
    lines = ["# OT σ-modification vs labelled σ-modification (oracle plant, m=6, 5 seeds, σ_m = 0.3)", "",
             "Same plant, same prior, same geometry and same x₀ in both starts; only the particle labels differ. "
             "RMS ‖e‖ over the second half of 40 s; W_end = OT error at T; lab_end = labelled error at T; "
             "b = W(μ*, μ⁰); b_lab = labelled bias seen by the labelled leak.", "",
             "| start | mechanism | RMS ‖e‖ (mean ± 95% CI) | W_end | lab_end | b | b_lab |", "|---|---|---|---|---|---|---|"]
    for st in ('aligned', 'relabelled'):
        for mech in MECHS:
            rr = [r for r in rows if r['start'] == st and r['mech'] == mech]
            v = np.array([r['rms'] for r in rr])
            lines.append(f"| {st} | {mech} | {v.mean():.4f} ± {1.96 * v.std(ddof=1) / np.sqrt(len(v)):.4f} | "
                         f"{np.mean([r['W_end'] for r in rr]):.4f} | {np.mean([r['lab_end'] for r in rr]):.4f} | "
                         f"{np.mean([r['b'] for r in rr]):.3f} | {np.mean([r['b_lab'] for r in rr]):.3f} |")
    lines += ["", "## Label invariance (same seed, aligned vs relabelled start): |RMS_aligned − RMS_relabelled|", "",
              "| mechanism | mean | max | relabelled worse in |", "|---|---|---|---|"]
    for mech in MECHS:
        pairs = [(a['rms'], b['rms']) for a in rows for b in rows
                 if a['mech'] == mech == b['mech'] and a['seed'] == b['seed'] and a['start'] == 'aligned' and b['start'] == 'relabelled']
        d = [abs(a - b) for a, b in pairs]
        lines.append(f"| {mech} | {np.mean(d):.2e} | {np.max(d):.2e} | {sum(b > a for a, b in pairs)}/{len(pairs)} |")
    lines += ["", "Reading: projection and OT σ-modification are exactly label-invariant; the labelled σ-modification is not — "
                  "with a relabelled start its leak drags every kernel toward a different kernel's prior position, "
                  "the bias it pays is b_lab instead of b, and tracking degrades accordingly."]
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    col = {'projection': 'k', 'labelled sigma-mod': 'r', 'OT sigma-mod': 'b'}
    for ax, st in zip(axs, ('aligned', 'relabelled')):
        for mech in MECHS:
            E = np.array([r['e'] for r in rows if r['start'] == st and r['mech'] == mech])
            ax.semilogy(rows[0]['t'], np.median(E, axis=0), color=col[mech], label=mech)
        ax.set_title(f'{st} start'); ax.set_xlabel('t'); ax.grid(alpha=.3)
    axs[0].set_ylabel('‖e(t)‖ (median of seeds)'); axs[0].legend()
    plt.tight_layout(); plt.savefig(os.path.join(OUT, 'figure.png'), dpi=160); print('saved figure')


if __name__ == "__main__":
    main()
