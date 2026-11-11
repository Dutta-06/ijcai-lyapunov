"""
Theorem C′: the OT leak toward a prior with a different number of atoms (continuous-time equal-mass Lloyd)
=========================================================================================================

Two tables.

1. Pure leak (no plant): m = 6 particles leak toward an m₀ = 48-atom prior in the lifted space (n = 1, p = 3).  W²(μ_t, ν)
   decreases monotonically to the value of an equal-mass centroidal configuration (a Lloyd fixed point), the centroid
   defect V_bar² → 0, and the final W² is compared with the best of 20 Lloyd runs from random starts.  With m₀ = m the
   same code reproduces the exact decay of Theorem C(i).

2. Closed loop: a representable plant with an m* = 4 ideal (ε ≡ 0) is controlled by an m = 4 network under three leaks of
   the same gain toward three priors: (A) the ideal itself (Theorem C, b = 0); (B) a 32-atom cloud around the ideal
   (8 jittered copies per ideal atom, barycentric leak, Theorem C′: b = W(μ*, ν) is the jitter scale but the centroid
   defect b_bar of the ideal is small); (C) the labelled leak toward a random 4-subset of the cloud (b_lab is the
   sub-sampling error).  Steady-state RMS‖e‖ and W(μ_T, μ*) over 3 seeds.

Output: results/barycentric_leak/table.md
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

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'barycentric_leak')
BOX = dict(a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)


def pure_leak_table():
    rng = np.random.default_rng(0)
    Ginv = 1.0 / np.array([2.0, 1.0, 0.5])
    m, m0, sigma, dt = 6, 48, 1.0, 1e-2
    Q = np.column_stack([rng.normal(0, 1.5, m0), rng.uniform(-3, 3, m0), rng.uniform(0.5, 2.5, m0)])
    rows = []
    for seed in range(3):
        r2 = np.random.default_rng(10 + seed)
        P = np.column_stack([r2.normal(0, 1.5, m), r2.uniform(-3, 3, m), r2.uniform(0.5, 2.5, m)])
        W2_0 = T.optimal_plan(P, Q, Ginv)[1]
        mono = True; W2_prev = W2_0; nsw = 0; gamma_prev = None
        for k in range(3000):
            vb, Tbar, gamma, W2, Var = T.centroid_defect(P, Q, Ginv)
            if W2 > W2_prev + 1e-9:
                mono = False
            if gamma_prev is not None and np.max(np.abs(gamma - gamma_prev)) > 1e-6:
                nsw += 1
            gamma_prev, W2_prev = gamma, W2
            P = P - dt * sigma * (P - Tbar)
        vb_end, _, _, W2_end, Var_end = T.centroid_defect(P, Q, Ginv)
        best = min(T.lloyd_equal_mass(np.column_stack([r3.normal(0, 1.5, m), r3.uniform(-3, 3, m), r3.uniform(0.5, 2.5, m)]), Q, Ginv)[1]
                   for r3 in (np.random.default_rng(100 + j) for j in range(20)))
        rows.append(dict(seed=seed, W2_0=W2_0, W2_end=W2_end, vb_end=vb_end, Var_end=Var_end, mono=mono, nsw=nsw, best=best))
    # m₀ = m: exact decay of Theorem C(i)
    Qm = Q[:m]; P = np.column_stack([rng.normal(0, 1.5, m), rng.uniform(-3, 3, m), rng.uniform(0.5, 2.5, m)])
    W0 = np.sqrt(T.optimal_plan(P, Qm, Ginv)[1]); worst = 0.0
    for k in range(500):
        Tbar = T.barycentric_projection(P, Qm, Ginv)[0]
        P = P + dt * (-sigma * (P - Tbar))                      # explicit Euler
        W = np.sqrt(T.optimal_plan(P, Qm, Ginv)[1])
        worst = max(worst, abs(W - (1 - sigma * dt) ** (k + 1) * W0) / W0)
    return rows, worst


def exchange_gap_bound(P, Q, Ginv):
    """Upper bound on Δ_gap at P: the smallest cost increase when two full atoms are swapped between two cells of the vertex
    optimal plan (the swapped plan is another vertex, so its cost is at least the next vertex cost)."""
    gamma, _ = T.optimal_plan(P, Q, Ginv)
    C = T.assignment_cost_matrix(P, Q, Ginv)
    m, m0 = gamma.shape
    full = gamma > 1.0 / m0 - 1e-9
    best = np.inf
    for i in range(m):
        for ip in range(i + 1, m):
            for j in np.flatnonzero(full[i]):
                for jp in np.flatnonzero(full[ip]):
                    inc = (C[i, jp] + C[ip, j] - C[i, j] - C[ip, jp]) / m0
                    if inc > 1e-12:
                        best = min(best, inc)
    return float(best)


def closed_loop_table():
    rows = []
    for seed in range(3):
        rng = np.random.default_rng(20 + seed)
        m, n = 4, 1
        base = T.TransportController(m, n, ke=4.0, g_a=2.0, g_c=1.0, g_s=0.5, **BOX)
        P_star = T.sample_lifted_ideal(m, n, base, rng, min_sep=0.8)
        plant = T.LiftedRepresentablePlant(P_star, amplitude=1.5)
        cloud = base.hard_clip(np.vstack([P_star[i] + rng.normal(0, 0.2, size=(8, 3)) for i in range(m)]))
        sub = cloud[rng.choice(len(cloud), size=m, replace=False)]
        P0 = base.hard_clip(P_star[rng.permutation(m)] + rng.normal(0, 0.4, size=P_star.shape))
        Ginv = base.Ginv_diag
        b_cloud = np.sqrt(T.optimal_plan(P_star, cloud, Ginv)[1]); bbar_cloud = np.sqrt(T.centroid_defect(P_star, cloud, Ginv)[0])
        b_sub = np.sqrt(T.optimal_assignment(P_star, sub, Ginv)[1])
        d = P_star - sub; b_lab = np.sqrt(np.einsum('ik,k,ik->', d, Ginv, d) / m)
        gap = exchange_gap_bound(P_star, cloud, Ginv)                                    # ≥ Δ_gap at μ*
        eta0 = np.sqrt(T.optimal_assignment(P0, P_star, Ginv)[1])                        # η(0) = W(μ₀, μ*)
        out = {}
        for name, leak, prior in (('A: OT leak, prior = ideal (m₀ = 4)', 'ot', P_star), ('B: barycentric leak, 32-atom cloud', 'ot_bary', cloud),
                                  ('C: labelled leak, random 4-subset of the cloud', 'labelled', sub)):
            ctrl = T.TransportController(m, n, ke=4.0, g_a=2.0, g_c=1.0, g_s=0.5, leak=leak, sigma_m=0.5, prior=prior, **BOX)
            res = T.rk4_closed_loop(plant, ctrl, plant.x_ref(0.0), P0, T=40.0, dt=2e-3, record_every=25)
            k1 = np.searchsorted(res['t'], 20.0)
            rms = float(np.sqrt(np.mean(res['e_norm'][k1:] ** 2)))
            W_end = float(np.sqrt(ctrl.ot_distance(res['P'][-1], P_star)[1]))
            out[name] = (rms, W_end)
        rows.append(dict(seed=seed, b_cloud=b_cloud, bbar_cloud=bbar_cloud, b_sub=b_sub, b_lab=b_lab, out=out, gap=gap, eta0=eta0))
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    rows, worst_exact = pure_leak_table()
    L = ["# Theorem C′: the barycentric OT leak (`exp_barycentric_leak.py`)", "",
         "## 1. Pure leak toward a 48-atom prior (m = 6, n = 1, lifted dimension 3, G = diag(2, 1, 0.5), σ_m = 1, dt = 0.01, 3000 steps)", "",
         "| seed | W²(μ₀, ν) | W²(μ_T, ν) | centroid defect V_bar² at T | within-cell variance Var_γ at T | W² non-increasing | plan switches | best of 20 equal-mass Lloyd runs |",
         "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['seed']} | {r['W2_0']:.4f} | {r['W2_end']:.4f} | {r['vb_end']:.1e} | {r['Var_end']:.4f} | {r['mono']} | {r['nsw']} | {r['best']:.4f} |")
    L += ["", "The flow stops exactly where W² = Var_γ (every particle at the barycentre of its cell): an equal-mass centroidal "
          "configuration of the prior. Like Lloyd's algorithm it finds a local minimiser of the quantisation energy, not always the "
          "global one: the best of 20 Lloyd runs lies between the end points of the three starts (last two columns).",
          f"With m₀ = m the same leak gives the exact decay of Theorem C(i): max relative deviation of W(μ_t, μ⁰) from (1 − σ_m dt)^k W₀ over 500 Euler steps: {worst_exact:.1e}.", "",
          "## 2. Closed loop: three leaks of equal gain (σ_m = 0.5) toward three priors; representable plant, m = 4, k_e = 4", "",
          "| seed | b = W(μ*, cloud) | b_bar (centroid defect of the ideal w.r.t. the cloud) | W(μ*, subset) | labelled bias to the subset | A: RMS‖e‖ [20, 40], W(μ_T, μ*) | B | C | exchange bound on Δ_gap at μ* (cloud) | η(0) = W(μ₀, μ*) | 4η(0)(b + η(0)) |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    ratios = []
    for r in closed_loop_table():
        o = r['out']; A, B, C = list(o.values())
        ratios.append(4 * r['eta0'] * (r['b_cloud'] + r['eta0']) / r['gap'])
        L.append(f"| {r['seed']} | {r['b_cloud']:.3f} | {r['bbar_cloud']:.3f} | {r['b_sub']:.3f} | {r['b_lab']:.3f} | {A[0]:.4f}, {A[1]:.3f} | {B[0]:.4f}, {B[1]:.3f} | {C[0]:.4f}, {C[1]:.3f} | "
                 f"{r['gap']:.4f} | {r['eta0']:.3f} | {4 * r['eta0'] * (r['b_cloud'] + r['eta0']):.3f} |")
    L += ["", "Reading: each loop settles at the configuration its leak anchors; for (B) the distance W(μ_T, μ*) equals b_bar, the "
          "distance of the ideal from the centroids of the prior's cells, not b = W(μ*, ν) (Theorem C′(i), (v)). The gap condition "
          "4η(b + η) < Δ_gap of Theorem C′(ii)–(iv) fails at t = 0 in every run: 4η(0)(b + η(0)) exceeds the exchange bound on Δ_gap "
          "(itself an upper bound on Δ_gap) by the factors " + ", ".join(f"{x:.0f}" for x in ratios) + ", so the local estimate (iv) is not "
          "what these runs test. The comparison with (C) shows the cost of a prior of the right size at the wrong positions; a labelled "
          "leak toward the cloud's own equal-mass 4-quantiser would carry the same bias b_bar (Theorem C′(v)).", "",
          "Numbers are produced by `experiments/exp_barycentric_leak.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == '__main__':
    main()
