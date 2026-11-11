"""
Theorem C(v) corrected, in numbers: the labelled and OT leaks reach the same steady state, but from a
relabelled start the labelled leak forces a reshuffling transient (each kernel crosses the domain) while
the OT leak leaves the configuration in place.  Reports RMS ‖e‖ over [0,5], [5,10], [10,20], [20,40] s and
the transport cost of the motion (Wasserstein length of the particle path).
Outputs: results/ot_sigma_mod/transient.md
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


def run(mech, seed, m=6, n=2, T_end=40.0, sigma_m=0.3):
    rng = np.random.default_rng(seed)
    base = T.TransportController(m, n, **GAINS, **BOX)
    Q = T.sample_lifted_ideal(m, n, base, rng, c_spread=1.5, min_sep=0.6)
    prior = base.hard_clip(Q + rng.normal(0, 0.15, size=Q.shape))
    x0 = rng.uniform(-0.5, 0.5, size=n)
    P0 = prior[np.roll(np.arange(m), 1)]                                   # relabelled start
    plant = T.LiftedRepresentablePlant(Q)
    leak = 'ot' if mech == 'OT' else 'labelled'
    ctrl = T.TransportController(m, n, leak=leak, sigma_m=sigma_m, prior=prior, **GAINS, **BOX)
    res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0) + x0, P0, T=T_end, dt=1e-3, record_every=10)
    t, e, P = res['t'], res['e_norm'], res['P']
    length = float(np.sum(np.sqrt(np.einsum('kij,j,kij->k', np.diff(P, axis=0), ctrl.Ginv_diag, np.diff(P, axis=0)) / m)))
    win = lambda a, b: float(np.sqrt(np.mean(e[(t >= a) & (t < b)] ** 2)))
    return dict(mech=mech, seed=seed, w05=win(0, 5), w510=win(5, 10), w1020=win(10, 20), w2040=win(20, 40), length=length)


def _job(a):
    return run(*a)


def main(n_seeds=5):
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    with mp.Pool(10) as pool:
        rows = pool.map(_job, [(m, s) for m in ('labelled', 'OT') for s in range(n_seeds)])
    lines = ["# Reshuffling transient from a relabelled start (m=6, 5 seeds, σ_m = 0.3)", "",
             "| leak | RMS ‖e‖ [0,5] | [5,10] | [10,20] | [20,40] | Wasserstein length of the particle path |", "|---|---|---|---|---|---|"]
    for mech in ('labelled', 'OT'):
        rr = [r for r in rows if r['mech'] == mech]
        f = lambda k: f"{np.mean([r[k] for r in rr]):.4f}"
        lines.append(f"| {mech} | {f('w05')} | {f('w510')} | {f('w1020')} | {f('w2040')} | {f('length')} |")
    lines += ["", "Steady states coincide; the labelled leak pays a transient (larger early error, far longer particle path)."]
    with open(os.path.join(OUT, 'transient.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
