"""
Theorem J′ in the real closed loop: the stationary law of the dithered single-kernel loop against the Gibbs measure
===================================================================================================================

Instance (n = 1, m = 1): plant f(x) = e^{−(x−0.5)²/0.5} + 0.4x, reference x_m(t) = 1.5[sin t′ + (1/3) sin 2t′] with
t′ = ±t (the two orientations have the same occupation measure), saturated amplitude and centre ϑ = (α, γ) with
a_max = 3, c_max = 2.5, frozen width s = 0.6, gains G = diag(2, 2), leak σ₀ = 0.4 toward ϑ⁰ = 0, temperature τ₀ = 0.1;
σ = σ₀/k_e, τ = τ₀/k_e, k_e ∈ {10, 20, 40, 80}.

Grid (201 × 201 on [−3.5, 3.5]²): the Gibbs measure π ∝ e^{−J/τ₀}, the first-order drift b = b_fr + b_osc, the correction
g (τ₀∇·(πG∇g) = ∇·(πb)), the current j₁ = π(b − τ₀G∇g), the entropy-production coefficient τ₀⁻¹π(|b_⊥|²_{G⁻¹}) and the
circulation coefficient ∫Θ∧j₁ + π⟨u∧u̇⟩.
Monte Carlo of the real loop: 20 000 independent replicas started from the first-order law π(1 + εg), Strang splitting with
k_e h = 1/2, 10 slow-time units of burn-in and 40 of averaging (whole numbers of reference periods), both orientations;
time averages of α, γ, α², γ², αγ, (k_e e)² and the circulation α V_γ − γ V_α (V = k_e × drift), with standard errors from
the spread of the per-replica time averages (the replicas are independent; batch means in time are correlated through the
slowest mode of L̄₀, whose relaxation time is 3.5 slow-time units on this instance — referee round 6).

Tables: (1) m_ε(φ) against π(φ) and π(φ) + ε∫φg dπ, with the decay of both errors in k_e; (2) orientation-odd part of
the stationary means against its first-order prediction; (3) the stationary tracking error k_e²E‖e‖² against 2π(L̄);
(4) the circulation observable against ε[∫Θ∧j₁ + π⟨u∧u̇⟩], and the current machinery checked deterministically on the
first-order averaged SDE (stationary Fokker–Planck solve on the grid: density π(1 + εg) + O(ε²), current εj₁ + O(ε²)).

Output: results/dither/table.md (and raw.npz).  Set DS_WORKERS, DS_K, DS_T to change the pool, the replicas and the horizon.
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import dither as Dt

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'dither')
KES = (10.0, 20.0, 40.0, 80.0)
K_REP = int(os.environ.get('DS_K', 20000))
T_SLOW = float(os.environ.get('DS_T', 40.0))
BURN = 10.0
N_WORKERS = int(os.environ.get('DS_WORKERS', 3))
LIM = (-3.5, 3.5)
NG = 201


def make_instance(orientation):
    f = lambda X: np.exp(-(X - 0.5) ** 2 / 0.5) + 0.4 * X
    df = lambda X: np.exp(-(X - 0.5) ** 2 / 0.5) * (-2 * (X - 0.5) / 0.5) + 0.4
    return Dt.SingleKernelInstance(f, df, a_max=3.0, c_max=2.5, s=0.6, G=(2.0, 2.0), sigma0=0.4, theta0=(0.0, 0.0), tau0=0.1,
                                   X=1.5, omega=1.0, r2=1.0 / 3.0, orientation=orientation, N_theta=512)


def mc_run(args):
    ke, orientation, seed = args
    inst = make_instance(orientation)
    grid = Dt.GibbsGrid(inst, LIM, LIM, n_al=121, n_ga=121)                 # grid only to sample the initial law π(1 + εg)
    rng = np.random.default_rng(seed)
    al0, ga0 = grid.sample(K_REP, rng, eps=1.0 / ke)                         # start from the first-order law: the transient is O(ε²)
    mean, se, n_steps = Dt.simulate_stationary(inst, ke, K_REP, T_SLOW, BURN, seed + 1, al0, ga0)
    return dict(ke=ke, orientation=orientation, mean=mean, se=se, n_steps=n_steps)


def first_order_sde_check(grid, epss=(0.2, 0.1, 0.05, 0.025)):
    """
    The current machinery checked on the first-order averaged SDE dΘ = (−G∇J + εb)ds + √(2τ₀)G^{1/2}dW, whose stationary
    density m_ε and current j_ε are computed on the grid (GibbsGrid.stationary_density): Proposition J′.2 predicts
    m_ε = π(1 + εg) + O(ε²) and ∫Θ∧j_ε = ε∫Θ∧j₁ + O(ε²).  The discrete current at ε = 0 (a discretisation residual of the
    central scheme, −6·10⁻⁵ on the 201² grid) is subtracted.  Returns rows (ε, ∫Θ∧(j_ε − j_0)/ε, ratio to ∫Θ∧j₁, ‖m_ε − m_0(1 + εg)‖₁/ε²).
    """
    m0, j0 = grid.stationary_density(0.0)
    w0 = float(np.sum(grid.AL * j0[1] - grid.GA * j0[0]) * grid.cell)
    _, wedge1, _ = grid.circulation_prediction()
    rows = []
    for eps in epss:
        m, j = grid.stationary_density(eps)
        w = float(np.sum(grid.AL * j[1] - grid.GA * j[0]) * grid.cell)
        l1 = float(np.sum(np.abs(m - m0 * (1 + eps * grid.g))) * grid.cell)
        rows.append((eps, (w - w0) / eps, (w - w0) / (eps * wedge1), l1 / eps ** 2))
    return rows, w0


def main():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    grids = {o: Dt.GibbsGrid(make_instance(o), LIM, LIM, n_al=NG, n_ga=NG) for o in (+1, -1)}
    g = grids[+1]
    obs = dict(alpha=g.AL, gamma=g.GA, alpha2=g.AL ** 2, gamma2=g.GA ** 2, alpha_gamma=g.AL * g.GA)
    names = list(obs)
    pi0 = {k: g.expect(v) for k, v in obs.items()}
    p1 = {o: {k: grids[o].first_order(v) for k, v in obs.items()} for o in (+1, -1)}
    Lbar2 = 2 * g.expect(g.data['Lbar'])
    circ = {o: grids[o].circulation_prediction() for o in (+1, -1)}
    epr = {o: grids[o].epr_coefficient() for o in (+1, -1)}
    b2 = {o: grids[o].expect(np.sum(grids[o].b ** 2 / grids[o].inst.G[:, None, None], axis=0)) for o in (+1, -1)}
    # the antisymmetric part of D(G⁻¹b): curl of G⁻¹b (zero iff b is a G-gradient); compare with the symmetric part's scale
    def curl_stats(gr):
        bg = gr.b / gr.inst.G[:, None, None]
        dA = gr.grad(bg[0]); dG = gr.grad(bg[1])
        curl = dG[0] - dA[1]
        sym = np.sqrt(dA[0] ** 2 + dG[1] ** 2 + 0.5 * (dA[1] + dG[0]) ** 2)
        return np.sqrt(gr.expect(curl ** 2)), np.sqrt(gr.expect(sym ** 2))
    curls = {o: curl_stats(grids[o]) for o in (+1, -1)}
    bnd = max(g.pi[0].max(), g.pi[-1].max(), g.pi[:, 0].max(), g.pi[:, -1].max()) / g.pi.max()
    jobs = [(ke, o, 100 + 7 * i) for i, (ke, o) in enumerate([(ke, o) for ke in KES for o in (+1, -1)])]
    jobs.sort(key=lambda j: -j[0])
    with mp.Pool(N_WORKERS) as pool:
        res = pool.map(mc_run, jobs)
    R = {(r['ke'], r['orientation']): r for r in res}
    # the deterministic check of the current machinery on the first-order averaged SDE
    sde_rows, w0 = first_order_sde_check(g)

    L = ["# Theorem J′: the stationary law of the dithered single-kernel loop (`exp_dither_stationary.py`)", "",
         "Instance: n = 1, m = 1, f(x) = e^{−(x−0.5)²/0.5} + 0.4x, x_m(t) = 1.5[sin t′ + (1/3) sin 2t′] (t′ = ±t), a_max = 3, c_max = 2.5, s = 0.6 frozen, "
         f"G = diag(2, 2), σ₀ = 0.4, ϑ⁰ = 0, τ₀ = 0.1; σ = σ₀/k_e, τ = τ₀/k_e. Grid {NG}×{NG} on [{LIM[0]}, {LIM[1]}]² (boundary density / peak density = {bnd:.1e}); "
         f"Monte Carlo: {K_REP} replicas started from the first-order law π(1 + εg), Strang splitting with k_e·h = 1/2, burn-in ≥ {BURN:g} and averaging ≥ {T_SLOW:g} slow-time units "
         "(whole numbers of reference periods), standard errors from the spread of the per-replica time averages (one standard error quoted).", "",
         "## Averaged data", "",
         f"Gibbs expectations: " + ", ".join(f"π({k}) = {pi0[k]:+.5f}" for k in names) + f"; 2π(L̄) = {Lbar2:.5f}.", "",
         "| orientation | rms b_fr | rms b_osc | π(‖b‖²_{G⁻¹}) | non-gradient fraction π(‖b_⊥‖²)/π(‖b‖²) | EPR coefficient τ₀⁻¹π(‖b_⊥‖²_{G⁻¹}) | rms curl(G⁻¹b) / rms sym D(G⁻¹b) | circulation coefficient: ∫Θ∧j₁ + π⟨u∧u̇⟩ |",
         "|---|---|---|---|---|---|---|---|"]
    for o in (+1, -1):
        gr = grids[o]
        L.append(f"| {o:+d} | {np.sqrt(gr.expect(np.sum(gr.data['b_fr'] ** 2, 0))):.4f} | {np.sqrt(gr.expect(np.sum(gr.data['b_osc'] ** 2, 0))):.4f} | {b2[o]:.4f} | "
                 f"{epr[o] * gr.inst.tau0 / b2[o]:.3f} | {epr[o]:.4f} | {curls[o][0] / curls[o][1]:.3f} | {circ[o][1]:+.5f} + {circ[o][2]:+.5f} = {circ[o][0]:+.5f} |")
    L += ["", "The first-order coefficients ∫φ g dπ of the stationary means (m_ε(φ) = π(φ) + ε∫φg dπ + O(ε²)):", "",
          "| φ | orientation +1 | orientation −1 | even part | odd part |", "|---|---|---|---|---|"]
    for k in names:
        L.append(f"| {k} | {p1[+1][k]:+.5f} | {p1[-1][k]:+.5f} | {(p1[+1][k] + p1[-1][k]) / 2:+.5f} | {(p1[+1][k] - p1[-1][k]) / 2:+.5f} |")
    # ---- table 1
    L += ["", "## Table 1 — stationary means of the real loop (orientation +1) against the Gibbs value and the first-order prediction", ""]
    for k, name in enumerate(names):
        L += [f"### φ = {name}", "", "| k_e | m_ε(φ) ± s.e. | π(φ) | π(φ) + ε∫φg dπ | m_ε − π | m_ε − first order | (m_ε − π)·k_e | (m_ε − first order)·k_e² |", "|---|---|---|---|---|---|---|---|"]
        e0, e1 = [], []
        for ke in KES:
            r = R[(ke, +1)]
            m, s = r['mean'][k], r['se'][k]
            d0 = m - pi0[name]; d1 = m - (pi0[name] + p1[+1][name] / ke)
            e0.append(d0); e1.append(d1)
            L.append(f"| {ke:g} | {m:+.5f} ± {s:.5f} | {pi0[name]:+.5f} | {pi0[name] + p1[+1][name] / ke:+.5f} | {d0:+.5f} | {d1:+.5f} | {d0 * ke:+.4f} | {d1 * ke ** 2:+.4f} |")
        sl0 = np.polyfit(np.log(KES), np.log(np.abs(e0)), 1)[0]
        ok1 = [abs(e1[i]) > 2 * R[(KES[i], +1)]['se'][k] for i in range(len(KES))]
        txt1 = f"{np.polyfit(np.log([KES[i] for i in range(len(KES)) if ok1[i]]), np.log([abs(e1[i]) for i in range(len(KES)) if ok1[i]]), 1)[0]:.2f} over the k_e where it is resolved ({sum(ok1)} of {len(KES)})" if sum(ok1) >= 2 else "not resolved beyond the Monte Carlo error"
        L += ["", f"Slope of log|m_ε − π| against log k_e: {sl0:.2f} (first order: −1); slope of log|m_ε − first order|: {txt1} (second order: −2).", ""]
    # ---- table 2: orientation-odd parts
    L += ["## Table 2 — orientation-odd part of the stationary means, (m₊ − m₋)/2, against its first-order prediction (p₁₊ − p₁₋)/(2k_e)", "",
          "| k_e | φ | measured odd part ± s.e. | predicted | measured/predicted | ratio·(1 ± s.e./measured) |", "|---|---|---|---|---|---|"]
    for ke in KES:
        for k, name in enumerate(names):
            rp, rm = R[(ke, +1)], R[(ke, -1)]
            odd = 0.5 * (rp['mean'][k] - rm['mean'][k]); se = 0.5 * np.hypot(rp['se'][k], rm['se'][k])
            pred = 0.5 * (p1[+1][name] - p1[-1][name]) / ke
            L.append(f"| {ke:g} | {name} | {odd:+.5f} ± {se:.5f} | {pred:+.5f} | {odd / pred:.3f} | ± {se / abs(pred):.3f} |")
    # ---- table 3: tracking error
    L += ["", "## Table 3 — stationary tracking error: k_e²·E‖e‖² against 2π(L̄) (Theorem J′(d): equal to O(1/k_e))", "",
          "| k_e | k_e²E e² (orientation +1) ± s.e. | (−1) | 2π(L̄) | ratio (+1) | ratio (−1) |", "|---|---|---|---|---|---|"]
    for ke in KES:
        rp, rm = R[(ke, +1)], R[(ke, -1)]
        L.append(f"| {ke:g} | {rp['mean'][5]:.5f} ± {rp['se'][5]:.5f} | {rm['mean'][5]:.5f} ± {rm['se'][5]:.5f} | {Lbar2:.5f} | {rp['mean'][5] / Lbar2:.4f} | {rm['mean'][5] / Lbar2:.4f} |")
    # ---- table 4: circulation
    L += ["", "## Table 4 — the circulation observable Λ_ε = E[Θ ∧ k_e·drift] against ε[∫Θ∧j₁ + π⟨u∧u̇⟩]", "",
          "| k_e | orientation | measured Λ_ε ± s.e. | predicted ε·coefficient | (measured − predicted)/s.e. |", "|---|---|---|---|---|"]
    for ke in KES:
        for o in (+1, -1):
            r = R[(ke, o)]
            pred = circ[o][0] / ke
            L.append(f"| {ke:g} | {o:+d} | {r['mean'][6]:+.5f} ± {r['se'][6]:.5f} | {pred:+.5f} | {(r['mean'][6] - pred) / r['se'][6]:+.2f} |")
    L += ["", "The predicted circulation is below the Monte Carlo resolution of these runs wherever |predicted| < 2 s.e.; the table is a consistency check, "
          "not a measurement of the current. The current machinery (Hodge decomposition of b, j₁ = π(b − τ₀G∇g), Proposition J′.2) is checked "
          "deterministically on the first-order averaged SDE dΘ = (−G∇J + εb)ds + √(2τ₀)G^{1/2}dW, whose stationary density m_ε and current j_ε are "
          f"computed on the grid (orientation +1; the discrete current at ε = 0, ∫Θ∧j₀ = {w0:+.2e}, a residual of the central scheme, is subtracted): "
          f"Proposition J′.2 predicts ∫Θ∧j_ε = ε∫Θ∧j₁ + O(ε²) with ∫Θ∧j₁ = {circ[+1][1]:+.5f}, and m_ε = π(1 + εg) + O(ε²).", "",
          "| ε | ∫Θ∧(j_ε − j₀)/ε | ratio to ∫Θ∧j₁ | ‖m_ε − m₀(1 + εg)‖₁/ε² |", "|---|---|---|---|"]
    L += [f"| {eps:g} | {w:+.5f} | {r:.4f} | {l1:.3f} |" for (eps, w, r, l1) in sde_rows] + [""]
    L += ["Numbers are produced by `experiments/exp_dither_stationary.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as fh:
        fh.write("\n".join(L) + "\n")
    np.savez(os.path.join(OUT, 'raw.npz'), **{f"mean_{int(r['ke'])}_{'p' if r['orientation'] > 0 else 'm'}": r['mean'] for r in res},
             **{f"se_{int(r['ke'])}_{'p' if r['orientation'] > 0 else 'm'}": r['se'] for r in res},
             pi0=np.array([pi0[k] for k in names]), p1_plus=np.array([p1[+1][k] for k in names]), p1_minus=np.array([p1[-1][k] for k in names]),
             Lbar2=Lbar2, circ_plus=np.array(circ[+1]), circ_minus=np.array(circ[-1]), sde=np.array(sde_rows))
    print("\n".join(L))


if __name__ == '__main__':
    main()
