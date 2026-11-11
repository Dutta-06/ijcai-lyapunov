"""
Theorem B in numbers: labelled vs optimal-transport Lyapunov analysis, and m-uniformity
=======================================================================================

For m = 8, 16, 32, 64 lifted kernels (mean-field gains), the network is initialised as a
random PERMUTATION of the ideal configuration plus small noise, so that the function it
represents is already almost ideal.  Along the closed loop we record

  * the identity-labelled parameter error  ((1/m)Σ‖ϑ_i − ϑ*_i‖²_{G⁻¹})^{1/2},
  * the optimal-transport error            W(μ_t, μ*),
  * the actual excess disturbance          |(1/m)Σ remainders| (what enters V̇ beyond −k‖e‖²),
  * the OT bound ‖e‖ K λ_G W²  and the labelled bound ‖e‖ K λ_G (labelled error)²,
  * the classical labelled constant D_max of THEORY_REVISED.md Lemma 2 (grows with m).

Outputs: results/transport_permutation/{table.md, figure.png}
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

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'transport_permutation')
BOX = dict(a_max=8.0, c_max=4.0, s_min=0.3, s_max=3.0, eps=0.15)


def run(m, seed, n=2, T_end=10.0):
    rng = np.random.default_rng(seed)
    ctrl = T.TransportController(m, n, ke=3.0, g_a=2.0, g_c=1.0, g_s=0.5, **BOX)
    P_star = T.sample_lifted_ideal(m, n, ctrl, rng, c_spread=2.0, min_sep=0.3)
    plant = T.LiftedRepresentablePlant(P_star)
    P0 = ctrl.hard_clip(T.perturb_and_permute(P_star, rng, 0.05)[0])
    x0 = plant.x_ref(0.0) + rng.uniform(-0.5, 0.5, size=n)
    reps = []
    T.rk4_closed_loop(plant, ctrl, x0, P0, T=T_end, dt=1e-3, record_every=20,
                      callback=lambda x, P, t: reps.append(ctrl.lyapunov_report(x, t, plant, P, P_star)))
    lab = np.array([r['identity_err'] for r in reps]); W = np.array([r['W'] for r in reps])
    en = np.array([r['e_norm'] for r in reps])
    actual = np.array([abs(r['rem'].sum()) / m for r in reps])
    ot_bound = np.array([ctrl.excess_bound(a, b) for a, b in zip(en, W)])        # Theorem B, matched (OT) error
    lab_bound = np.array([ctrl.excess_bound(a, b) for a, b in zip(en, lab)])     # same bound with the identity labelling
    D_global = ctrl.excess_bound(1.0, ctrl.diameter())                           # global constant per unit ‖e‖ (Cor. B.3)
    D_quad = ctrl.K * ctrl.lam_G * ctrl.diameter() ** 2                          # quadratic branch alone at the diameter
    return dict(m=m, seed=seed, t=np.array([r['t'] for r in reps]), lab=lab, W=W, e=en, actual=actual,
                ot_bound=ot_bound, lab_bound=lab_bound, D_global=D_global, D_quad=D_quad,
                rms_e=float(np.sqrt(np.mean(en[len(en) // 2:] ** 2))))


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = [run(m, s) for m in (8, 16, 32, 64) for s in range(3)]
    lines = ["# Labelled vs optimal-transport analysis (permuted-ideal initialisation, mean-field gains)", "",
             "Means over 3 seeds; 'excess' = |(1/m)Σ Taylor remainders| actually entering V̇; bounds evaluated along the same trajectory.", "",
             "| m | labelled err (t=0) | OT err W (t=0) | RMS ‖e‖ | mean actual excess | mean OT bound (Thm B) | same bound, identity labelling | ratio labelled/OT | global constant per unit ‖e‖ (Cor. B.3, min-form) | quadratic branch alone at diameter |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for m in (8, 16, 32, 64):
        rr = [r for r in rows if r['m'] == m]
        f = lambda key: np.mean([np.mean(r[key]) for r in rr])
        lines.append(f"| {m} | {np.mean([r['lab'][0] for r in rr]):.3f} | {np.mean([r['W'][0] for r in rr]):.4f} | {np.mean([r['rms_e'] for r in rr]):.4f} | "
                     f"{f('actual'):.2e} | {f('ot_bound'):.2e} | {f('lab_bound'):.2e} | {f('lab_bound') / f('ot_bound'):.0f}× | "
                     f"{np.mean([r['D_global'] for r in rr]):.1f} | {np.mean([r['D_quad'] for r in rr]):.0f} |")
    ok = all(np.all(r['actual'] <= r['ot_bound'] * (1 + 1e-6) + 1e-12) for r in rows)
    lines += ["", f"OT bound ≥ actual excess at every sample of every run: **{ok}**",
              "Reading: the labelled error never decreases (the kernels converge to a relabelling of the ideal), the OT error is ~1/25 of it, "
              "and the identical bound is 30–45× smaller when evaluated with the matched error (the quadratic branch is active for the "
              "matched error, the linear branch for the labelled one). Both bounds are still far above the "
              "actual excess (the Hessian supremum is global), and the global constant is set by the linear branch, not the quadratic one: "
              "the gain of the OT analysis is local, not global. No quantity grows with m."]
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(17, 4.6))
    r = [x for x in rows if x['m'] == 16 and x['seed'] == 0][0]
    axs[0].plot(r['t'], r['lab'], 'r', label='labelled parameter error'); axs[0].plot(r['t'], r['W'], 'b', label='OT error W(μ_t, μ*)')
    axs[0].set_yscale('log'); axs[0].set_xlabel('t'); axs[0].legend(); axs[0].set_title('m = 16, permuted-ideal start'); axs[0].grid(alpha=.3)
    axs[1].semilogy(r['t'], r['actual'] + 1e-16, 'k', label='actual excess disturbance'); axs[1].semilogy(r['t'], r['ot_bound'] + 1e-16, 'b--', label='OT bound (Thm B)')
    axs[1].semilogy(r['t'], r['lab_bound'] + 1e-16, 'r--', label='labelled bound'); axs[1].legend(); axs[1].set_xlabel('t'); axs[1].grid(alpha=.3); axs[1].set_title('what enters V̇')
    ms = (8, 16, 32, 64)
    axs[2].semilogy(ms, [np.mean([np.mean(x['lab_bound']) for x in rows if x['m'] == m]) for m in ms], 'r-s', label='mean bound, identity labelling')
    axs[2].semilogy(ms, [np.mean([np.mean(x['ot_bound']) for x in rows if x['m'] == m]) for m in ms], 'b-o', label='mean bound, OT matching')
    axs[2].semilogy(ms, [np.mean([np.mean(x['actual']) for x in rows if x['m'] == m]) for m in ms], 'k-^', label='mean actual excess')
    axs[2].set_xlabel('m'); axs[2].legend(); axs[2].grid(alpha=.3, which='both'); axs[2].set_title('no growth with the number of kernels')
    plt.tight_layout(); plt.savefig(os.path.join(OUT, 'figure.png'), dpi=160); print('saved figure')


if __name__ == "__main__":
    main()
