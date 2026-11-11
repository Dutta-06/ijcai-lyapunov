"""
Theorem D in numbers: convergence of the m-particle closed loop as m grows (mean-field gains)
==============================================================================================

A fixed smooth initial measure μ_0 on the lifted box is quantised deterministically with
m = 16, 32, 64, 128, 256 particles (stratified quantile sampling, so W(μ^m_0, μ_0) decreases
with m).  The plant is the same Van der Pol oscillator for every m (a genuine, non-representable
plant).  We record the tracking-error trajectory e^m(t) and the function f̂^m on a fixed grid,
and report sup_t ‖e^m − e^{2m}‖ and the W-distance between successive particle measures at
t = T, which Theorem D predicts to be controlled by their distance at t = 0.

Outputs: results/mean_field_m/{table.md, figure.png}
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
from src.core.plant_models import VanDerPolPlant

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'mean_field_m')
BOX = dict(a_max=6.0, c_max=3.0, s_min=0.3, s_max=2.5, eps=0.1)
GAINS = dict(ke=4.0, g_a=3.0, g_c=1.5, g_s=0.5)


def quantise_initial_measure(m, n=2):
    """
    Deterministic quantisation of a smooth one-parameter initial measure μ_0 = (ϑ(u))#Uniform[0,1],
    ϑ(u) = (a(u), c(u), s(u)) with a(u) = (sin 2πu, cos 4πu), c(u) = 2(cos 2πu, sin 2πu), s(u) = 1 + 0.4 sin 6πu,
    at the midpoints u_i = (i+½)/m.  Then W(μ^m_0, μ_0) = O(1/m), so the m-dependence is visible.
    """
    u = (np.arange(m) + 0.5) / m
    P = np.empty((m, 2 * n + 1))
    P[:, 0] = np.sin(2 * np.pi * u); P[:, 1] = np.cos(4 * np.pi * u)
    P[:, 2] = 2 * np.cos(2 * np.pi * u); P[:, 3] = 2 * np.sin(2 * np.pi * u)
    P[:, 4] = 1.0 + 0.4 * np.sin(6 * np.pi * u)
    return P


def run(m, T_end=15.0, dt=1e-3):
    plant = VanDerPolPlant(mu=1.5, ref_amplitude=1.2, ref_omega=1.0)
    ctrl = T.TransportController(m, 2, **GAINS, **BOX)
    P0 = quantise_initial_measure(m)
    res = T.rk4_closed_loop(plant, ctrl, np.array([0.5, -0.5]), P0, T=T_end, dt=dt, record_every=10)
    return ctrl, res


def main():
    os.makedirs(OUT, exist_ok=True)
    ms = [16, 32, 64, 128, 256]
    runs = {m: run(m) for m in ms}
    lines = ["# Mean-field consistency in m (Van der Pol, mean-field gains, smooth initial measure quantised at m points)", "",
             "Theorem D predicts sup_t ‖e^m − e^{2m}‖ + W(μ^m_t, μ^{2m}_t) ≤ e^{Λt} W(μ^m_0, μ^{2m}_0); the initial distance decays like 1/m.", "",
             "| m | sup_t ‖e^m − e^{2m}‖ | W(μ^m_0, μ^{2m}_0) | W(μ^m_T, μ^{2m}_T) | ratio W_T / W_0 | RMS ‖e‖ (last half) |", "|---|---|---|---|---|---|"]
    for m in ms[:-1]:
        c1, r1 = runs[m]; c2, r2 = runs[2 * m]
        sup_e = float(np.max(np.linalg.norm(r1['x'] - r2['x'], axis=1)))
        # W between measures of different size: use the 2m-particle configuration duplicated vs doubled m-configuration
        P1 = np.repeat(r1['P'][0], 2, axis=0); P1T = np.repeat(r1['P'][-1], 2, axis=0)
        W0 = T.w2(P1, r2['P'][0], c2.Ginv_diag); WT = T.w2(P1T, r2['P'][-1], c2.Ginv_diag)
        k = len(r1['e_norm']) // 2
        lines.append(f"| {m} vs {2 * m} | {sup_e:.4f} | {W0:.4f} | {WT:.4f} | {WT / W0:.2f} | {np.sqrt(np.mean(r1['e_norm'][k:] ** 2)):.4f} |")
    c, r = runs[ms[-1]]; k = len(r['e_norm']) // 2
    lines.append(f"| {ms[-1]} | — | — | — | — | {np.sqrt(np.mean(r['e_norm'][k:] ** 2)):.4f} |")
    with open(os.path.join(OUT, 'table.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for m in ms:
        ax.semilogy(runs[m][1]['t'], runs[m][1]['e_norm'], label=f'm = {m}', lw=1)
    ax.set_xlabel('t'); ax.set_ylabel('‖e(t)‖'); ax.legend(); ax.grid(alpha=.3); ax.set_title('closed loop converges as m grows (mean-field gains)')
    plt.tight_layout(); plt.savefig(os.path.join(OUT, 'figure.png'), dpi=160); print('saved figure')


if __name__ == "__main__":
    main()
