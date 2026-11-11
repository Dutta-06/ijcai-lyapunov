"""
Monte Carlo over random initial conditions — honest version
===========================================================

Replaces tier1_monte_carlo.py, whose "R_UUB threshold" was the literal
constant 0.5 chosen by hand ("A chosen arbitrary bound for the test") and then
reported as "100% of trials converged strictly below the calculated R_UUB".

What CAN be verified on a real (non-representable) plant:
  1. the hard parameter bounds (W, c, σ) ∈ Ω_W × Ω_c × Ω_σ at every step — exact;
  2. boundedness and the distribution of the steady-state tracking error
     across random initial states, weights, centres and bandwidths;
  3. an a-posteriori ISS bound that follows from the error dynamics
     ė = −kₑe + (f − f̂) alone:  for t ≥ T₁,
        ‖e(t)‖ ≤ e^{−kₑ(t−T₁)}‖e(T₁)‖ + sup_{[T₁,t]}‖f − f̂‖ / kₑ,
     evaluated with the measured sup of ‖f − f̂‖ (this is the quantity the
     a-priori D_max over-estimates by orders of magnitude).

What CANNOT be verified here: the a-priori R_UUB, because ε_max, V_W, V_c, V_σ
(and hence D_max, C_Ω) are not computable for a plant outside the network's
span.  No threshold is therefore drawn; see exp_oracle_ruub.py for the setting
in which every constant is computable.

Outputs: results/monte_carlo_honest/{summary.md, raw.npz, figure.png}
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import json
import numpy as np
import multiprocessing as mp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from src.core.plant_models import VanDerPolPlant
from src.core.controller import AdaptiveController
from src.core.integrator import AdaptiveRK4Integrator, SimulationState

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'monte_carlo_honest')
KE, DT, T_END, M, N = 5.0, 2e-3, 20.0, 16, 2


def run_single_trial(seed):
    rng = np.random.default_rng(seed)
    plant = VanDerPolPlant(mu=1.5, ref_amplitude=1.0, ref_omega=1.2)
    x0 = rng.uniform(-3.0, 3.0, size=N)
    g = np.linspace(-2.0, 2.0, 4)
    c0 = np.array(np.meshgrid(g, g)).T.reshape(-1, 2) + rng.normal(0, 0.2, size=(M, N))
    W0 = rng.uniform(-2.0, 2.0, size=(M, N))
    s0 = rng.uniform(0.8, 2.5, size=M)
    ctrl = AdaptiveController(m=M, n=N, ke=KE, gamma_W=20.0, gamma_c=5.0, gamma_sigma=1.0,
                              W_max=20.0, c_max=10.0, sigma_min=0.2, sigma_max=5.0, delta_min=0.1)
    res = AdaptiveRK4Integrator(plant, ctrl, dt=DT).simulate(SimulationState(x0, W0, c0, s0), 0.0, T_END,
                                                             record_interval=5)
    d = res.to_arrays()
    t, e, x, W, c, s = d['t'], d['e_norm'], d['x'], d['W'], d['c'], d['sigma']
    inside = all(ctrl.params_inside_sets(W[i], c[i], s[i]) for i in range(len(t)))
    # measured model error along the trajectory
    ftil = np.array([np.linalg.norm(plant.f(x[i], t[i]) - ctrl.kernel.f_hat(x[i], W[i], c[i], s[i]))
                     for i in range(len(t))])
    k1 = int(0.5 * len(t)); k2 = int(0.75 * len(t))
    T1 = t[k1]
    iss_bound = np.exp(-KE * (t[k1:] - T1)) * e[k1] + np.maximum.accumulate(ftil[k1:]) / KE
    return dict(seed=seed, finite=bool(np.all(np.isfinite(e))), inside=bool(inside),
                max_ss=float(e[k2:].max()), mean_ss=float(e[k2:].mean()), max_all=float(e.max()),
                sup_ftil_ss=float(ftil[k2:].max()), iss_ok=bool(np.all(e[k1:] <= iss_bound * (1 + 1e-9))),
                iss_ratio=float(np.max(e[k1:] / iss_bound)), init_dist=float(np.linalg.norm(x0)),
                min_sigma=float(s.min()), max_W=float(np.max(np.linalg.norm(W, axis=(1, 2)))))


def main(n_trials=200):
    os.makedirs(OUT, exist_ok=True)
    with mp.Pool(max(1, mp.cpu_count() - 1)) as pool:
        rows = pool.map(run_single_trial, range(n_trials))
    np.savez(os.path.join(OUT, 'raw.npz'), **{k: np.array([r[k] for r in rows]) for k in rows[0]})
    max_ss = np.array([r['max_ss'] for r in rows]); mean_ss = np.array([r['mean_ss'] for r in rows])
    lines = [f"# Monte Carlo over random initial conditions (Van der Pol, m={M}, k_e={KE}, {n_trials} trials)", "",
             "| check | result |", "|---|---|",
             f"| all trajectories finite | {sum(r['finite'] for r in rows)}/{n_trials} |",
             f"| hard parameter bounds held at every recorded step | {sum(r['inside'] for r in rows)}/{n_trials} |",
             f"| a-posteriori ISS bound ‖e(t)‖ ≤ e^{{-kₑ(t-T₁)}}‖e(T₁)‖ + sup‖f−f̂‖/kₑ held | {sum(r['iss_ok'] for r in rows)}/{n_trials} (max ratio {max(r['iss_ratio'] for r in rows):.3f}) |",
             f"| steady-state max ‖e‖ (last 25%): median / 95th pct / max | {np.median(max_ss):.4f} / {np.percentile(max_ss, 95):.4f} / {max_ss.max():.4f} |",
             f"| steady-state mean ‖e‖: median / 95th pct / max | {np.median(mean_ss):.4f} / {np.percentile(mean_ss, 95):.4f} / {mean_ss.max():.4f} |",
             f"| sup‖f − f̂‖ in steady state / kₑ: median / max | {np.median([r['sup_ftil_ss'] for r in rows]) / KE:.4f} / {max(r['sup_ftil_ss'] for r in rows) / KE:.4f} |",
             f"| transient max ‖e‖: max over trials | {max(r['max_all'] for r in rows):.3f} |",
             "", "No a-priori R_UUB threshold is drawn: D_max and C_Ω are not computable for this plant "
             "(the 0.5 threshold of the original figure was an arbitrary constant, not a calculated bound)."]
    with open(os.path.join(OUT, 'summary.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    axs[0].hist(max_ss, bins=30, color='skyblue', edgecolor='k')
    axs[0].set_xlabel('max steady-state ‖e‖ (last 25%)'); axs[0].set_ylabel('trials'); axs[0].set_title('Steady-state error across random initial conditions')
    axs[1].scatter([r['init_dist'] for r in rows], max_ss, alpha=0.6, edgecolors='none')
    axs[1].set_xlabel('‖x(0)‖'); axs[1].set_ylabel('max steady-state ‖e‖'); axs[1].set_title('No a-priori threshold is available for this plant')
    plt.tight_layout(); plt.savefig(os.path.join(OUT, 'figure.png'), dpi=160)
    print("saved", os.path.join(OUT, 'figure.png'))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200)
