"""
Theorem E in numbers: matched (OT) parameter error converges under tangential excitation; the
identity-labelled error does not (particles converge to a relabelling of the ideal).
Also reports the tangent-regressor Gramian ∫_t^{t+60}Ψ(s)ᵀΨ(s)ds along the reference: its smallest
eigenvalue on the first window and the MINIMUM over all 60 s windows inside the horizon (the uniform
level that the (PE) hypothesis of Theorem E needs), its largest eigenvalue, and the fitted exponential
decay rate of the matched error.  The PE level enters the guaranteed rate only through constants that
also depend on the gains and on the tube (Remark E.1); the table shows that it does not by itself
predict the observed rate.

Outputs: results/pe_convergence/{table.md, figure.png}
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

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'pe_convergence')
BOX = dict(a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)


def gramian(plant, P_star, T_win=60.0, dt=0.02):
    G = 0
    for t in np.arange(0, T_win, dt):
        Psi = np.concatenate([g for g in T.grad_Phi(plant.x_ref(t), P_star)], axis=1)
        G = G + Psi.T @ Psi * dt
    return np.linalg.eigvalsh(G)


def sliding_gramian_min(plant, P_star, T_end, T_win=60.0, stride=2.0, dt=0.02):
    """min over t ∈ [0, T_end − T_win] (step `stride`) of λ_min ∫_t^{t+T_win} ΨᵀΨ ds, via cumulative sums."""
    ts = np.arange(0, T_end, dt)
    cum = [0.0]
    for t in ts:
        Psi = np.concatenate([g for g in T.grad_Phi(plant.x_ref(t), P_star)], axis=1)
        cum.append(cum[-1] + Psi.T @ Psi * dt)
    w, s_ = int(round(T_win / dt)), int(round(stride / dt))
    lam = [float(np.linalg.eigvalsh(cum[k + w] - cum[k]).min()) for k in range(0, len(ts) - w + 1, s_)]
    return float(min(lam)), float(max(lam))


def fitted_rate(t, W, t0):
    """least-squares slope of −log W on t ≥ t0."""
    mk = t >= t0
    return float(-np.polyfit(t[mk], np.log(W[mk]), 1)[0])


def run(seed, rich=True, m=3, n=2, T_end=200.0):
    rng = np.random.default_rng(seed)
    ctrl = T.TransportController(m, n, ke=6.0, g_a=30.0, g_c=30.0, g_s=15.0, **BOX)
    P_star = T.sample_lifted_ideal(m, n, ctrl, rng, min_sep=0.9)
    plant = T.RichReferencePlant(P_star, amplitude=2.5, n_tones=6) if rich else T.LiftedRepresentablePlant(P_star, amplitude=2.5)
    perm = np.roll(np.arange(m), 1)
    P0 = ctrl.hard_clip(P_star[perm] + rng.normal(0, 0.05, size=P_star.shape))
    res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0), P0, T=T_end, dt=1e-3, record_every=500)
    W = np.array([ctrl.ot_distance(P, P_star)[1] ** 0.5 for P in res['P']])
    lab = np.array([np.sqrt(np.einsum('ik,k,ik->', P - P_star, ctrl.Ginv_diag, P - P_star) / m) for P in res['P']])
    ev = gramian(plant, P_star)
    sw_min, sw_max = sliding_gramian_min(plant, P_star, T_end)
    return dict(seed=seed, rich=rich, t=res['t'], W=W, lab=lab, e=res['e_norm'], ev_min=float(ev.min()), ev_max=float(ev.max()),
                sw_min=sw_min, sw_max=sw_max, rate=fitted_rate(res['t'], W, T_end / 2))


def _job(args):
    return run(*args)


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    jobs = [(s, rich) for rich in (True, False) for s in range(3)]
    with mp.Pool(len(jobs)) as pool:
        rows = pool.map(_job, jobs)
    lines = ["# Tangential excitation and label-free parameter convergence (m=3, n=2, permuted start, T = 200 s)", "",
             "| reference | seed | λ_min(Gramian, first 60 s) | min over 60 s windows | λ_max (first window) | W(0) | W(T)/W(0) | "
             "fitted rate of W, t ≥ T/2 (1/s) | labelled(T)/labelled(0) | max ‖e‖ |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {'rich (6 tones)' if r['rich'] else 'single sine'} | {r['seed']} | {r['ev_min']:.2e} | {r['sw_min']:.2e} | "
                     f"{r['ev_max']:.1f} | {r['W'][0]:.4f} | {r['W'][-1] / r['W'][0]:.3f} | {r['rate']:.2e} | "
                     f"{r['lab'][-1] / r['lab'][0]:.3f} | {r['e'].max():.4f} |")
    rates, levels = np.array([r['rate'] for r in rows]), np.array([r['sw_min'] for r in rows])
    rho_s = float(np.corrcoef(np.argsort(np.argsort(levels)), np.argsort(np.argsort(rates)))[0, 1])
    lines += ["", "The matched (OT) error decays in every run; the identity-labelled error stays O(1) because the kernels converge "
                  "to a relabelling of the ideal configuration.",
              f"Rank correlation between the uniform excitation level (min over windows) and the fitted rate over the 6 runs: "
              f"{rho_s:.2f}.  The excitation level is a hypothesis of Theorem E, not a rate predictor (Remark E.1)."]
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
    for ax, rich in zip(axs, (True, False)):
        for r in [x for x in rows if x['rich'] == rich]:
            ax.semilogy(r['t'], r['W'], 'b', alpha=.7, label='OT error W(μ_t, μ*)' if r['seed'] == 0 else None)
            ax.semilogy(r['t'], r['lab'], 'r', alpha=.7, label='identity-labelled error' if r['seed'] == 0 else None)
        ax.set_title('rich reference (PE)' if rich else 'single-sine reference'); ax.set_xlabel('t'); ax.grid(alpha=.3); ax.legend()
    plt.tight_layout(); plt.savefig(os.path.join(OUT, 'figure.png'), dpi=160); print('saved figure')


if __name__ == "__main__":
    main()
