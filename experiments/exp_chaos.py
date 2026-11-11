"""
Theorem D′: propagation of chaos for the dithered closed loop, checked by synchronous coupling
=============================================================================================

An M-particle system (M = 4096) stands in for the McKean–Vlasov closed loop. For m ∈ {16, 32, …, 1024} the m-particle
system is driven by the first m initial particles and the first m Brownian motions of the M-particle system (the
synchronous coupling of Theorem D′(b)), and the coupled distance

    D_m(t) := ‖x^m_t − x^M_t‖ + (1/m)Σ_{i≤m}‖ϑ^m_i(t) − ϑ^M_i(t)‖

is recorded. Theorem D′(b) and the per-particle bound of D′(c) give E sup_t D_m ≤ C(T)(m^{−1/2} + M^{−1/2}) (the second term
from replacing the limit by the M-system). Because the m-system shares its particles with the reference, the fluctuation of the mean field
that drives D_m has variance proportional to 1/m − 1/M, so the expected sup is compared with (1/m − 1/M)^{1/2} as well as
with m^{−1/2}. sup_t D_m varies by up to a factor 10 between seeds at fixed m (the initial fluctuation of the empirical
measure is amplified anisotropically by the slow dynamics; the diagnostic of the round-5 review found no particle deviating
by more than 0.17, so no basin switching), which is why the expectation is estimated from N_SEEDS seeds with its standard
error, the slope carries a bootstrap standard error over seeds, and the median is reported alongside. Two temperatures (τ = 0 and τ = 2·10⁻³), n = 1, the plant and reference of
the Langevin experiments.

Output: results/chaos/table.md
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
PER = 2 * np.pi
M_REF = 4096
MS = (16, 32, 64, 128, 256, 512, 1024)
T_END = float(os.environ.get('CH_T', 20.0))
DT = 1e-3
KE, SIGMA, GAINS = 10.0, 0.02, np.array([4.0, 1.0, 1.0])
N_SEEDS = int(os.environ.get('CH_SEEDS', 24))
N_WORKERS = int(os.environ.get('CH_WORKERS', 2))


class Plant:
    def f(self, x, t):
        return np.sin(2 * x) + 0.3 * x

    def x_ref(self, t):
        return np.array([1.5 * np.sin(t) + 0.5 * np.sin(2 * t)])

    def x_ref_dot(self, t):
        return np.array([1.5 * np.cos(t) + 1.0 * np.cos(2 * t)])


def coupled_run(seed, tau):
    """Simulate the M-system and all m-systems with shared noise and initial particles; return sup_t D_m for each m."""
    rng = np.random.default_rng(seed)
    net = L.SaturatedLiftedNetwork(1, a_max=4.0, c_max=3.0, s_min=0.15, s_max=2.0)
    ctrl = L.LangevinAdaptiveController(net, ke=KE, G_diag=GAINS, sigma=SIGMA, theta0=np.zeros(3), tau=tau)
    plant = Plant()
    Th_all = net.unconstrained(rng.normal(0, 1.0, (M_REF, 1)), rng.uniform(-2, 2, (M_REF, 1)), rng.uniform(0.3, 1.0, M_REF))
    systems = {m: dict(x=plant.x_ref(0.0).copy(), Th=Th_all[:m].copy()) for m in MS + (M_REF,)}
    sup_D = {m: 0.0 for m in MS}
    sqG = np.sqrt(ctrl.G)
    n_steps = int(round(T_END / DT))
    t = 0.0
    for k in range(n_steps):
        xr, xrd = plant.x_ref(t), plant.x_ref_dot(t)
        noise = rng.standard_normal((M_REF, 3)) if tau > 0 else None
        new = {}
        for m, S in systems.items():
            x, Th = S['x'], S['Th']
            e = x - xr
            u = ctrl.control(x, xr, xrd, Th)
            k1 = plant.f(x, t) + u
            k2 = plant.f(x + 0.5 * DT * k1, t + 0.5 * DT) + u
            k3 = plant.f(x + 0.5 * DT * k2, t + 0.5 * DT) + u
            k4 = plant.f(x + DT * k3, t + DT) + u
            dTh = ctrl.drift(x, e, Th) * DT
            if tau > 0:
                dTh += np.sqrt(2 * tau * DT) * sqG[None, :] * noise[:m]
            new[m] = dict(x=x + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4), Th=Th + dTh)
        systems = new
        t = (k + 1) * DT
        if k % 20 == 0:
            ref = systems[M_REF]
            for m in MS:
                D = abs(float(systems[m]['x'][0] - ref['x'][0])) + float(np.mean(np.linalg.norm(systems[m]['Th'] - ref['Th'][:m], axis=1)))
                sup_D[m] = max(sup_D[m], D)
    return sup_D


def _job(args):
    seed, tau = args
    r = coupled_run(seed, tau)
    return seed, tau, [r[m] for m in MS]


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    seeds = tuple(range(N_SEEDS))
    taus = (0.0, 2e-3)
    with mp.Pool(N_WORKERS) as pool:
        out = pool.map(_job, [(s, tau) for tau in taus for s in seeds])
    res = {tau: np.array([[v for (s, t, v) in out if t == tau and s == seed][0] for seed in seeds]) for tau in taus}   # (seeds, MS)
    lines = ["# Theorem D′: synchronous coupling of the m-particle dithered loop to an M = 4096 reference (`exp_chaos.py`)", "",
             f"n = 1, k_e = {KE:g}, σ = {SIGMA:g}, G = diag{tuple(float(g) for g in GAINS)}, T = {T_END:g} s, dt = {DT:g}, Euler–Maruyama for the particles and RK4 "
             "for the plant (as `langevin.simulate`). The m-system shares its initial particles and Brownian motions with the first m particles of "
             f"the reference (M = {M_REF}). sup_t D_m is averaged over {N_SEEDS} seeds (mean ± standard error of the mean); the median and the spread over "
             "seeds are shown because single seeds differ by up to an order of magnitude at fixed m. The shared-particle coupling makes the "
             "driving fluctuation proportional to (1/m − 1/M)^{1/2}, so the mean is also divided by that factor (a constant row means the "
             "predicted law).", ""]
    ms = np.array(MS, float)
    xi = np.sqrt(1.0 / ms - 1.0 / M_REF)
    for tau in taus:
        R = res[tau]
        mean, sem, med = R.mean(0), R.std(0, ddof=1) / np.sqrt(len(seeds)), np.median(R, 0)
        sel = ms <= 256
        slope_m = np.polyfit(np.log(ms[sel]), np.log(mean[sel]), 1)[0]
        slope_xi = np.polyfit(np.log(xi), np.log(mean), 1)[0]
        slope_med = np.polyfit(np.log(ms[sel]), np.log(med[sel]), 1)[0]
        rng = np.random.default_rng(0)                                   # bootstrap over seeds (1000 resamples)
        bs_m, bs_xi = [], []
        for _ in range(1000):
            idx = rng.integers(0, len(seeds), len(seeds))
            mb = R[idx].mean(0)
            bs_m.append(np.polyfit(np.log(ms[sel]), np.log(mb[sel]), 1)[0])
            bs_xi.append(np.polyfit(np.log(xi), np.log(mb), 1)[0])
        se_m, se_xi = float(np.std(bs_m)), float(np.std(bs_xi))
        lines += [f"## τ = {tau:g}", "", "| m | mean sup_t D_m ± s.e. | median | min, max over seeds | mean / (1/m − 1/M)^{1/2} | mean · m^{1/2} |", "|---|---|---|---|---|---|"]
        for j, m in enumerate(MS):
            lines.append(f"| {m} | {mean[j]:.3e} ± {sem[j]:.1e} | {med[j]:.3e} | {R[:, j].min():.2e}, {R[:, j].max():.2e} | {mean[j] / xi[j]:.3f} | {mean[j] * np.sqrt(m):.3f} |")
        lines += ["", f"Log-log slope of the mean against m over m ≤ 256: {slope_m:.2f} ± {se_m:.2f} (bootstrap s.e. over seeds); of the median: {slope_med:.2f} (Theorem D′: −0.5). "
                  f"Slope of the mean against (1/m − 1/M)^(1/2) over all m: {slope_xi:.2f} ± {se_xi:.2f} (predicted 1).", ""]
    np.savez(os.path.join(OUT, 'raw_sup_D.npz'), ms=np.array(MS), seeds=np.array(seeds), **{f"tau_{tau:g}": res[tau] for tau in taus})
    lines.append("Numbers are produced by `experiments/exp_chaos.py` (raw per-seed values in `raw_sup_D.npz`); nothing is edited by hand.")
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == '__main__':
    main()
