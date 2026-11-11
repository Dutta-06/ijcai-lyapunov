"""
Theorem H in numbers: parameter dither buys global kernel placement
====================================================================

Scalar plant ẋ = sin(2x) + 0.3x + u, periodic reference x_m = 1.5 sin t + 0.5 sin 2t (occupation measure ρ),
saturated lifted network with m = 32 particles, σ-modification toward ϑ⁰ = 0 in the G⁻¹ metric,
dither temperature τ ∈ {0, τ₁, annealed}.  Two initialisations of the kernel measure:

  'collapsed' — all 32 particles at one point, a narrow kernel at the edge of the orbit (a symmetric,
                spurious configuration: without noise the particles can never separate);
  'spread'    — centres spread over the orbit, random small weights.

Prediction (Theorem H + Chizat 2022 / Nitanda–Wu–Suzuki 2022): with dither, both initialisations reach the
same occupation loss L̄ and the same tracking RMS (the averaged mean-field dynamics has a unique global
minimiser of F_τ and converges to it exponentially); without dither the collapsed start is stuck.

Metrics (3 seeds): occupation loss L̄(μ_T), tracking RMS over the last 10 periods, W₂ between the final
measures of the two starts (same seed), all in the physical coordinates (a, c, s).
Outputs: results/langevin_global/{table.md, figure.png}
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
from src.core.plant_models import PlantModel

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'langevin_global')
TAG = os.environ.get('LG_TAG', '')                     # output suffix, e.g. 'sigma0.02_tau0.004'
M, PERIOD, DT = 32, 2 * np.pi, 2e-3
T_END = float(os.environ.get('LG_T', 900.0))
GAINS = [4.0, 1.0, 1.0]
SIGMA = float(os.environ.get('LG_SIGMA', 0.002))
TAU1 = float(os.environ.get('LG_TAU', 0.0005))
KE = 10.0


class ScalarPlant(PlantModel):
    def __init__(self):
        super().__init__(n=1, name='scalar')

    def f(self, x, t):
        return np.array([np.sin(2 * x[0]) + 0.3 * x[0]])

    def x_ref(self, t):
        return np.array([1.5 * np.sin(t) + 0.5 * np.sin(2 * t)])

    def x_ref_dot(self, t):
        return np.array([1.5 * np.cos(t) + np.cos(2 * t)])


def init(kind, net, rng):
    if kind == 'collapsed':
        a = np.full((M, 1), 0.2); c = np.full((M, 1), 2.6); s = np.full(M, 0.25)
    else:
        a = rng.normal(0, 0.2, size=(M, 1)); c = rng.uniform(-2.0, 2.0, size=(M, 1)); s = rng.uniform(0.4, 1.2, size=M)
    return net.unconstrained(a, c, s)


def schedule(kind):
    if kind == 'none':
        return 0.0, None
    if kind == 'constant':
        return TAU1, None
    return TAU1, (lambda t: TAU1 * np.log(np.e + 20.0) / np.log(np.e + 20.0 + t))     # annealed, ∝ 1/log t


def run(args):
    dither, start, seed = args
    rng = np.random.default_rng(100 + seed)
    plant = ScalarPlant(); net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.15, s_max=2.0)
    tau, sched = schedule(dither)
    ctrl = L.LangevinAdaptiveController(net, ke=KE, G_diag=GAINS, sigma=SIGMA, theta0=np.zeros(3), tau=tau, tau_schedule=sched)
    Th0 = init(start, net, np.random.default_rng(7 + seed))
    res = L.simulate(plant, ctrl, plant.x_ref(0.0), Th0, T=T_END, dt=DT, rng=rng, record_every=500)
    t, e = res['t'], res['e_norm']
    mask = t >= T_END - 10 * PERIOD
    Ls = np.array([L.occupation_loss(plant, net, Th, PERIOD, 200) for Th in res['Th'][::10]])
    return dict(dither=dither, start=start, seed=seed, rms=float(np.sqrt(np.mean(e[mask] ** 2))),
                Lbar=float(L.occupation_loss(plant, net, res['Th'][-1], PERIOD)), Th_end=res['Th'][-1],
                L_path=Ls, t_L=t[::10], phys_end=np.hstack([*net.physical(res['Th'][-1])[:2], net.physical(res['Th'][-1])[2][:, None]]))


def w2_phys(P, Q):
    C = np.sum((P[:, None, :] - Q[None, :, :]) ** 2, axis=2)
    r, c = linear_sum_assignment(C)
    return float(np.sqrt(C[r, c].mean()))


def offline_best_loss(n_restarts=5):
    """Reference: unregularised min of L̄ over m = 32 kernels in the same box (L-BFGS from spread starts)."""
    from scipy.optimize import minimize
    plant = ScalarPlant(); net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.15, s_max=2.0)
    best = np.inf
    for k in range(n_restarts):
        Th0 = init('spread', net, np.random.default_rng(500 + k))
        res = minimize(lambda v: L.occupation_loss(plant, net, v.reshape(M, 3), PERIOD, 200), Th0.ravel(),
                       method='L-BFGS-B', options=dict(maxiter=3000))
        best = min(best, res.fun)
    return float(best)


def main(n_seeds=3):
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    jobs = [(d, s0, k) for d in ('none', 'constant', 'annealed') for s0 in ('collapsed', 'spread') for k in range(n_seeds)]
    with mp.Pool(min(len(jobs), 12)) as pool:
        rows = pool.map(run, jobs)
    best = offline_best_loss()
    plant = ScalarPlant()
    F_rho = max(abs(plant.f(plant.x_ref(t), 0.0)[0]) for t in np.linspace(0, PERIOD, 400))
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.15, s_max=2.0)
    ctrl = L.LangevinAdaptiveController(net, ke=KE, G_diag=GAINS, sigma=SIGMA, tau=TAU1)
    log_rate = ctrl.log_rate_lower_bound(F_rho)                         # log(2τρ_τ); the rate itself underflows
    lines = ["# Langevin adaptive control: dither and global kernel placement (scalar plant, m = 32, 3 seeds)", "",
             f"σ = {SIGMA} (leak toward ϑ⁰ = 0 in the G⁻¹ metric), k_e = {KE:g}, G = diag{tuple(GAINS)}, τ₁ = {TAU1}, T = {T_END:.0f} s. "
             "L̄ = occupation loss ½⟨(f − f̂)²⟩_ρ at T; RMS over the last 10 periods; W₂(collapsed, spread) between the final measures "
             "in (a, c, s) for the same seed.", "",
             f"Reference: unregularised offline optimum of L̄ over 32 kernels in the same box (L-BFGS, 5 restarts): **{best:.2e}**. "
             f"Guaranteed rate of Theorem H(ii) at τ₁ (Holley–Stroock lower bound on the averaged dynamics): "
             f"2τρ_τ ≥ e^({log_rate:.0f}) s⁻¹, i.e. vacuous at this temperature; the table does not test that rate.", "",
             "| dither | start | L̄(μ_T) (mean over seeds) | tracking RMS | W₂(final collapsed, final spread) |", "|---|---|---|---|---|"]
    for d in ('none', 'constant', 'annealed'):
        for s0 in ('collapsed', 'spread'):
            rr = [r for r in rows if r['dither'] == d and r['start'] == s0]
            w = [w2_phys(a['phys_end'], b['phys_end']) for a in rows for b in rows
                 if a['dither'] == d == b['dither'] and a['seed'] == b['seed'] and a['start'] == 'collapsed' and b['start'] == 'spread']
            lines.append(f"| {d} | {s0} | {np.mean([r['Lbar'] for r in rr]):.2e} | {np.mean([r['rms'] for r in rr]):.4f} | "
                         f"{np.mean(w):.3f} |")
    lines += ["", "Reading (computed from the table, not asserted):"]
    for d in ('none', 'constant', 'annealed'):
        Lc = np.mean([r['Lbar'] for r in rows if r['dither'] == d and r['start'] == 'collapsed'])
        Ls = np.mean([r['Lbar'] for r in rows if r['dither'] == d and r['start'] == 'spread'])
        same = abs(Lc - Ls) <= 0.2 * max(Lc, Ls)
        lines.append(f"- dither = {d}: collapsed/spread loss ratio {Lc / Ls:.2f} — "
                     + ("both starts reach the same loss (within 20%)." if same else
                        "the starts have NOT equalised within the horizon."))
    suffix = f"_{TAG}" if TAG else ""
    with open(os.path.join(OUT, f'table{suffix}.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.5), sharey=True)
    for ax, d in zip(axs, ('none', 'constant', 'annealed')):
        for s0, col in (('collapsed', 'r'), ('spread', 'b')):
            for r in [r for r in rows if r['dither'] == d and r['start'] == s0]:
                ax.semilogy(r['t_L'], r['L_path'], col, alpha=.6, label=s0 if r['seed'] == 0 else None)
        ax.set_title(f'dither: {d}'); ax.set_xlabel('t'); ax.grid(alpha=.3); ax.legend()
    axs[0].set_ylabel('occupation loss L̄(μ_t)')
    plt.tight_layout(); plt.savefig(os.path.join(OUT, f'figure{suffix}.png'), dpi=160); print('saved figure')


if __name__ == "__main__":
    main()
