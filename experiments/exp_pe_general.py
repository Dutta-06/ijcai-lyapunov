"""
Lemma E.0″: the slow-regime strict Lyapunov function for non-periodic persistently exciting regressors
======================================================================================================

Lemma E.0′'s construction with the period mean replaced by the forward exponential average over a window T,
M̄ = ⟨M⟩_T (TṀ̄ = M̄ − M), P₁ = −2T(M̄ − M̂), S = −T(F̄ − F̂), P₂ = −(M − M̄) + S, which leaves a residual second-order term
N̄₂ = ½(F̄ − Ṃ̄) and costs δ₂ = ε sup‖N̄₂‖/λ in the rate: γ/2 = (1 − δ − δ₂)λ/(k(1 + εp₁ + ε²p₂)), λ = inf λ_min(M̄_T).

Two regressors that are not periodic:
  1. quasi-periodic: Ψ(t) = (sin t, cos √2 t)  (n = 1, q = 2; the long-time mean of M is ½I);
  2. chaotic: Ψ(t) = (x₁(t/5)/20, (x₃(t/5) − 25)/20) along a Lorenz orbit (σ = 10, r = 28, b = 8/3) traversed five times slower
     than standard Lorenz time, with Ψ̇ by the chain rule through the Lorenz vector field (n = 1, q = 2).
For each: λ_T and the other constants at windows T, the conditions (A)–(D), δ₂, the guaranteed rate γ/2, the certificate
sup V̇/V on the time grid, and the exact decay rate of (LTV) (top Lyapunov exponent over the horizon) at k = 20, 50, 100, 200, 400.

Output: results/explicit_rate/general_pe.md
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import excitation as X

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'explicit_rate')
KES = (20.0, 50.0, 100.0, 200.0, 400.0)
T_END = 300.0


def quasi_periodic():
    r2 = np.sqrt(2.0)
    Psi = lambda t: np.array([[np.sin(t), np.cos(r2 * t)]])
    Psid = lambda t: np.array([[np.cos(t), -r2 * np.sin(r2 * t)]])
    return dict(name='quasi-periodic Ψ = (sin t, cos √2 t)', Psi=Psi, Psid=Psid, windows=(2 * np.pi, 4 * np.pi, 8 * np.pi), lam_inf=0.5)


def lorenz_regressor(slow=5.0):
    """Ψ(t) = (x₁(t/τ_s)/20, (x₃(t/τ_s) − 25)/20) along a Lorenz orbit (σ = 10, r = 28, b = 8/3) traversed τ_s = 5 times slower than
    the standard Lorenz time (the reference's time scale is a design choice: the derivative peak factors scale like τ_s⁻², τ_s⁻⁴);
    Ψ̇ by the chain rule through the Lorenz vector field.  n = 1, q = 2."""
    sig, r, b = 10.0, 28.0, 8.0 / 3.0
    f = lambda t, x: np.array([sig * (x[1] - x[0]), x[0] * (r - x[2]) - x[1], x[0] * x[1] - b * x[2]])
    T_lor = (T_END + 12.0 * 20.0 + 20.0) / slow + 1.0            # covers the padded grid of the longest window (T_pad = 12T, T ≤ 20)
    sol = solve_ivp(f, (0.0, T_lor), np.array([1.0, 1.0, 20.0]), method='DOP853', rtol=1e-11, atol=1e-12, t_eval=np.arange(0.0, T_lor, 5e-4))
    sp = CubicSpline(sol.t, sol.y.T, extrapolate=False)
    t_max = float(sol.t[-1])

    def _state(t):
        assert 0.0 <= t / slow <= t_max, f'the Lorenz orbit was integrated to Lorenz time {t_max}; t/τ_s = {t / slow} is outside it'
        return sp(t / slow)

    def Psi(t):
        x = _state(t)
        return np.array([[x[0] / 20.0, (x[2] - 25.0) / 20.0]])

    def Psid(t):
        x = _state(t)
        v = f(0.0, x) / slow
        return np.array([[v[0] / 20.0, v[2] / 20.0]])
    return dict(name=f'chaotic Ψ = (x₁/20, (x₃ − 25)/20) along a Lorenz orbit traversed {slow:g} times slower than standard time', Psi=Psi, Psid=Psid,
                windows=(5.0, 10.0, 20.0), lam_inf=None)


def analyse(inst):
    rows = []
    exact = {k: X.lyapunov_exponent_rate(inst['Psi'], k, 200.0) for k in KES}
    for T in inst['windows']:
        E = X.ExpAverageLyapunov(inst['Psi'], inst['Psid'], T, T_end=T_END, dt=0.005)
        head = dict(T=T, lam=E.lam_min, Lam=E.lam_max, rho=E.rho, rho_D=E.rho_D, rho_DD=E.rho_DD, Psi_M=E.Psi_M, Psi_D=E.Psi_D, Psi_DD=E.Psi_DD,
                    p1=E.p1_meas, p2=E.p2_meas, pi1=E.pi1_meas, pi2=E.pi2_meas, n2=E.N2_norm)
        for k in KES:
            for delta in (0.2, 0.1):
                c = E.constants(k, delta)
                cx = E.constants(k, delta, explicit=True)
                sup, mg = E.certificate(k, delta, n_t=400) if c['ok'] else (float('nan'), float('nan'))
                rows.append(dict(T=T, k=k, delta=delta, c=c, cx=cx, sup=sup, mg=mg, exact=exact[k], head=head))
    return rows, exact


def main():
    os.makedirs(OUT, exist_ok=True)
    L = ["# Lemma E.0″: the slow-regime strict Lyapunov function on non-periodic PE regressors (`exp_pe_general.py`)", "",
         f"Exponential averages over a window T on [0, {T_END:g}] (grid step 0.005, the averages formed from the extended grid); the certificate "
         "is the supremum of V̇/V over 400 times and all states; the exact decay rate of (LTV) is −(1/200)log‖Φ(200, 0)‖₂ (top Lyapunov exponent "
         "over 200 time units). Constants: 'measured' uses the grid suprema of ‖P₁‖, ‖P₂‖, ‖P₁v‖/‖v‖_{M̄}, ‖P₂v‖/‖v‖_{M̄}, ‖N̄₂‖ (computable constants "
         "of the regressor, like the peak factors); 'explicit' uses the closed-form bounds π₁ = 2T(√Λ + ‖M̂‖/√λ) etc.", ""]
    for inst in (quasi_periodic(), lorenz_regressor()):
        rows, exact = analyse(inst)
        L += [f"## {inst['name']}", "",
              "Exact decay rates of (LTV): " + ", ".join(f"k = {k:g}: {exact[k]:.5f}" for k in KES)
              + (f" (the long-time mean of M has λ_min = {inst['lam_inf']}, so λ_min/k = " + ", ".join(f"{inst['lam_inf'] / k:.5f}" for k in KES) + ")" if inst['lam_inf'] else "") + ".", ""]
        for T in inst['windows']:
            h = [r for r in rows if r['T'] == T][0]['head']
            L += [f"### window T = {T:.4g}: λ = inf λ_min(M̄_T) = {h['lam']:.4f}, Λ = {h['Lam']:.4f}, ρ = {h['rho']:.2f}, ρ_D = {h['rho_D']:.2f}, ρ_DD = {h['rho_DD']:.2f}, "
                  f"Ψ_M = {h['Psi_M']:.3f}, Ψ_D = {h['Psi_D']:.3f}, Ψ_DD = {h['Psi_DD']:.3f}, p₁ = {h['p1']:.2f}, p₂ = {h['p2']:.2f}, π₁ = {h['pi1']:.2f}, π₂ = {h['pi2']:.2f}, ‖N̄₂‖ = {h['n2']:.3f}", "",
                  "| k | δ | conditions (A)–(D), measured | δ₂ | guaranteed rate γ/2 | exact rate | exact / guaranteed | λ_T/k | sup V̇/V on the grid, −γ | explicit constants: (A)–(D), rate |",
                  "|---|---|---|---|---|---|---|---|---|---|"]
            for r in [r for r in rows if r['T'] == T]:
                c, cx = r['c'], r['cx']
                conds = "".join(k_ for k_, v in c['conds'].items() if v) or "none"
                condx = "".join(k_ for k_, v in cx['conds'].items() if v) or "none"
                L.append(f"| {r['k']:g} | {r['delta']} | {conds} | {c['delta2']:.3f} | {c['rate']:.5f} | {r['exact']:.5f} | "
                         f"{r['exact'] / c['rate'] if c['ok'] else float('nan'):.2f} | {c['averaged_rate']:.5f} | {r['sup']:.5f}, {r['mg']:.5f} | {condx}, {cx['rate']:.5f} |")
            L.append("")
    L += ["Reading: wherever the lemma's conditions (A)–(D) hold, the certificate sup V̇/V ≤ −γ holds on the grid, and the guaranteed rate is "
          "(1 − δ − δ₂)λ_T/k up to the factor 1/(1 + εp₁ + ε²p₂), with λ_T the exponentially weighted PE level of the window and δ₂ the price of the "
          "residual second-order term, which decreases like 1/k. The window trades λ_T (closer to the long-time excitation for longer windows) against "
          "p₁ and p₂ (which grow with T) and the residual ‖N̄₂‖ (no monotone rule: it decreases with T on the quasi-periodic regressor and is "
          "smallest at T = 10 on the chaotic one); the ratio exact/guaranteed is the product of the window loss λ_∞/λ_T, the factor 1/(1 − δ − δ₂) and "
          "the O(1/k) factors. Where a condition fails the row says which: (A) fails when the derivative peak factors ρ_D, ρ_DD are large relative "
          "to the excitation (a chaotic regressor whose components are nearly collinear has a tiny λ_T and astronomical peak factors, and the "
          "threshold in k moves beyond the range), (C) fails when δ₂ = εn₂/λ_T is not small. The explicit closed-form constants are looser than the "
          "measured suprema and move the threshold to larger k.", "",
          "Numbers are produced by `experiments/exp_pe_general.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'general_pe.md'), 'w', encoding='utf-8') as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == '__main__':
    main()
