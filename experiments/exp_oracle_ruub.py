"""
Oracle R_UUB study: what the theoretical bound actually says
=============================================================

Replaces tier2_ruub_sweep.py, whose "theoretical R_UUB" was built from invented
constants (ε_max = 0.5, V_W = 5, …) and a made-up formula
D_max = ε_max(1 + 2 W_max c_max).

Here the plant is EXACTLY representable with a known constant θ* (ε = 0,
C_Ω = 0), so every constant in the manuscript's bound is computable:

  D_theory  = F_max + 3W_max√m + W_max(J_c·2c_max + J_σ√m(σ_max−σ_min))   (Lemma 2)
  R_theory  = D_theory / kₑ                                                (C_Ω = 0)
  R_full    = √(R_theory² + 2 V_max,Θ)                                     (pointwise, Sec. 7.2)
  D_meas    = sup_t ‖d(t)‖ along the trajectory (oracle measurement)
  R_meas    = D_meas / kₑ        and the finite-window RMS bound (oracle.rms_bound)

and compared with the observed RMS and max steady-state error for a sweep over
kₑ and the learning-rate scale.  Expected (and honest) outcome: the a-priori
bounds are orders of magnitude loose; the RMS bound with the measured
disturbance holds and is within a small factor; the bound scales as 1/kₑ.

Outputs: results/oracle_ruub/{table.md, raw.json, figure.png}
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

from src.core.controller import AdaptiveController
from src.core.integrator import SimulationState
from src.core.oracle import RepresentablePlant, sample_admissible_ideal, run_oracle, rms_bound, rms_over_window

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'oracle_ruub')
BOUNDS = dict(W_max=10.0, c_max=6.0, sigma_min=0.3, sigma_max=3.0, proj_epsilon=0.2)
M, N, T, T1, DT = 8, 2, 30.0, 10.0, 1e-3


def run_one(args):
    ke, gscale, seed = args
    rng = np.random.default_rng(seed)
    Ws, cs, ss = sample_admissible_ideal(M, N, BOUNDS['W_max'], BOUNDS['c_max'], BOUNDS['sigma_min'],
                                         BOUNDS['sigma_max'], BOUNDS['proj_epsilon'], rng)
    plant = RepresentablePlant(Ws, cs, ss)
    ctrl = AdaptiveController(m=M, n=N, ke=ke, gamma_W=5.0 * gscale, gamma_c=2.0 * gscale,
                              gamma_sigma=1.0 * gscale, delta_min=0.5, **BOUNDS)
    x0 = rng.uniform(-1, 1, size=N)
    c0 = cs + rng.normal(0, 0.5, size=(M, N))
    s0 = np.clip(ss + rng.normal(0, 0.3, size=M), BOUNDS['sigma_min'] + BOUNDS['proj_epsilon'],
                 BOUNDS['sigma_max'] - BOUNDS['proj_epsilon'])
    res = run_oracle(ctrl, plant, SimulationState(x0, np.zeros((M, N)), c0, s0), T=T, dt=DT, every=5)
    rms, mask = rms_over_window(res, T1)
    D_meas = float(np.max(res['d_norm'][mask]))
    D_theory = ctrl.theoretical_D_max(plant.F_max())
    V_T1 = float(res['V'][mask][0])
    return dict(ke=ke, gscale=gscale, seed=seed,
                rms=rms, max_ss=float(np.max(res['e_norm'][mask])),
                D_meas=D_meas, D_theory=float(D_theory),
                R_meas_rms_bound=float(rms_bound(D_meas, ke, V_T1, T - T1)),
                R_meas=D_meas / ke, R_theory=float(D_theory / ke),
                R_full=float(np.sqrt((D_theory / ke) ** 2 + 2 * ctrl.parameter_energy_cap())),
                P_violations=int(np.sum(res['cross'] > 1e-10)), identity_err=float(np.max(np.abs(res['Vdot'] - res['rhs']))),
                inside=bool(np.all(res['inside'])))


def main():
    os.makedirs(OUT, exist_ok=True)
    kes = [2.0, 5.0, 10.0, 20.0]; gsc = [0.5, 1.0, 2.0, 5.0]; seeds = [0, 1, 2]
    jobs = [(k, g, s) for k in kes for g in gsc for s in seeds]
    with mp.Pool(max(1, mp.cpu_count() - 1)) as pool:
        rows = pool.map(run_one, jobs)
    with open(os.path.join(OUT, 'raw.json'), 'w') as f:
        json.dump(rows, f, indent=1)
    lines = ["# Oracle R_UUB study (exactly representable plant, constant θ*, m=8, n=2)", "",
             f"Window [{T1:.0f}, {T:.0f}] s; mean over {len(seeds)} seeds. All bounds in the same units as ‖e‖.",
             "`P viol` = samples with cross > 0 (must be 0); `inside` = hard bounds held.", "",
             "| k_e | Γ scale | RMS ‖e‖ | max SS ‖e‖ | RMS bound (D_meas, finite window) | D_meas/k_e | D_theory/k_e | full-state bound | ratio theory/RMS | P viol | inside |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for k in kes:
        for g in gsc:
            rr = [r for r in rows if r['ke'] == k and r['gscale'] == g]
            mean = lambda key: float(np.mean([r[key] for r in rr]))
            lines.append(f"| {k:g} | {g:g} | {mean('rms'):.4f} | {mean('max_ss'):.4f} | {mean('R_meas_rms_bound'):.4f} | "
                         f"{mean('R_meas'):.4f} | {mean('R_theory'):.1f} | {mean('R_full'):.1f} | "
                         f"{mean('R_theory') / mean('rms'):.0f}× | {sum(r['P_violations'] for r in rr)} | "
                         f"{all(r['inside'] for r in rr)} |")
    ok = all(r['rms'] <= r['R_meas_rms_bound'] * (1 + 1e-6) for r in rows)
    lines += ["", f"RMS bound (measured disturbance) satisfied in all {len(rows)} runs: **{ok}**",
              f"Max identity error |V̇ − (−kₑ‖e‖² + eᵀd + cross)|: {max(r['identity_err'] for r in rows):.2e}"]
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    g_mid = 1.0
    for key, lab, st in [('rms', 'observed RMS ‖e‖', 'o-'), ('max_ss', 'observed max SS ‖e‖', 's-'),
                         ('R_meas_rms_bound', 'RMS bound, measured D', '^--'), ('R_theory', 'a-priori R_UUB (Lemma 2)', 'v--'),
                         ('R_full', 'a-priori full-state bound', 'x:')]:
        ys = [np.mean([r[key] for r in rows if r['ke'] == k and r['gscale'] == g_mid]) for k in kes]
        axs[0].loglog(kes, ys, st, label=lab)
    axs[0].set_xlabel('k_e'); axs[0].set_ylabel('error / bound'); axs[0].set_title(f'Scaling with k_e (Γ scale = {g_mid})')
    axs[0].grid(True, which='both', alpha=0.3); axs[0].legend(fontsize=8)
    k_mid = 5.0
    for key, lab, st in [('rms', 'observed RMS ‖e‖', 'o-'), ('R_meas_rms_bound', 'RMS bound, measured D', '^--'),
                         ('R_theory', 'a-priori R_UUB', 'v--')]:
        ys = [np.mean([r[key] for r in rows if r['ke'] == k_mid and r['gscale'] == g]) for g in gsc]
        axs[1].loglog(gsc, ys, st, label=lab)
    axs[1].set_xlabel('learning-rate scale'); axs[1].set_title(f'Scaling with Γ (k_e = {k_mid})')
    axs[1].grid(True, which='both', alpha=0.3); axs[1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(os.path.join(OUT, 'figure.png'), dpi=160)
    print("saved", os.path.join(OUT, 'figure.png'))


if __name__ == "__main__":
    main()
