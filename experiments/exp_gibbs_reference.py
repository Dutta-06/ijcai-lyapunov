"""
Reference values for Theorem H and Conjecture H.6 (n = 1 Langevin experiment)
==============================================================================

Adapted from the round-3 reviewer's scripts (review/REVIEW_LOG.md, round 3), so that the numbers quoted in
PARADIGM.md §6c are produced by a project script.

1. μ*_τ, the minimiser of F_τ = (1/k_e)L̄ + τH(·|γ_τ) + const, γ_τ = N(ϑ⁰ = 0, (τ/σ)G), computed on a grid through the
   concave dual
        D(λ) = ⟨λ, f⟩_ρ − (k_e/2)‖λ‖²_ρ − τ log ∫ exp(⟨λ, Φ̃_ϑ⟩_ρ/τ) dγ_τ(ϑ),
   whose maximiser λ* = r_{μ*}/k_e gives μ*_τ ∝ γ_τ exp(⟨λ*, Φ̃_ϑ⟩_ρ/τ).  The product structure
   Φ̃(x;ϑ) = a(α)φ(x; c(γ), s(κ)) lets the α-integral be done exactly on a 1-D grid.
2. The stationary state of the EXACTLY averaged loop at fixed k_e (the small-gain limit in Conjecture H.6's
   discussion): μ ∝ γ_τ exp(P_μ/τ), P_μ(ϑ) = a·⟨e_μ(t) φ(x_m + e_μ; c, s)⟩_t, e_μ the periodic solution of
   ė = r_μ(x_m + e) − k_e e; solved as a fixed point in f̂ with Anderson mixing.

Setting: f(x) = sin 2x + 0.3x, x_m = 1.5 sin t + 0.5 sin 2t, saturated network a_max = 4, c_max = 3, s_min = 0.15,
s_max = 2, G = diag(4, 1, 1), ϑ⁰ = 0.

Outputs: results/langevin_global/gibbs_reference.md
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, expit

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'langevin_global')


def setup(sigma=0.02, tau=0.004, G=(4.0, 1.0, 1.0), ke=10.0, Nt=256, Ng=241, Nk=241, Na=1201,
          ga_rng=6.0, ka_rng=6.0, al_rng=14.0):
    P = dict(sigma=sigma, tau=tau, ke=ke)
    var = np.array(G) * tau / sigma
    t = np.arange(Nt) * 2 * np.pi / Nt
    xm = 1.5 * np.sin(t) + 0.5 * np.sin(2 * t)
    f = lambda x: np.sin(2 * x) + 0.3 * x
    al = np.linspace(-al_rng, al_rng, Na); ga = np.linspace(-ga_rng, ga_rng, Ng); ka = np.linspace(-ka_rng, ka_rng, Nk)
    a = 4.0 * np.tanh(al / 4.0)
    c = 3.0 * np.tanh(ga / 3.0)
    s = 0.15 + 1.85 * expit(ka)
    logNa = -al ** 2 / (2 * var[0]) + np.log(al[1] - al[0]) - 0.5 * np.log(2 * np.pi * var[0])
    logNg = -ga ** 2 / (2 * var[1]) + np.log(ga[1] - ga[0]) - 0.5 * np.log(2 * np.pi * var[1])
    logNk = -ka ** 2 / (2 * var[2]) + np.log(ka[1] - ka[0]) - 0.5 * np.log(2 * np.pi * var[2])
    C, S = np.meshgrid(c, s, indexing='ij')
    P.update(t=t, xm=xm, f=f, fx=f(xm), a=a, logNa=logNa, Cf=C.ravel(), Sf=S.ravel(),
             logw_gk=(logNg[:, None] + logNk[None, :]).ravel())
    return P


def phi_mat(P, X):
    return np.exp(-(X[:, None] - P['Cf'][None, :]) ** 2 / (2 * P['Sf'][None, :] ** 2))


def gibbs_from_h(P, h):
    """μ ∝ γ_τ exp(a h(c,s)/τ): marginal weights on the (γ, κ) grid, E[a | γ, κ], log-normaliser."""
    E = P['a'][None, :] * h[:, None] / P['tau'] + P['logNa'][None, :]
    logZ = logsumexp(E, axis=1)
    Ea = np.exp(E - logZ[:, None]) @ P['a']
    lw = logZ + P['logw_gk']
    LZ = logsumexp(lw)
    return np.exp(lw - LZ), Ea, LZ


def solve_mu_star(P):
    ke, tau = P['ke'], P['tau']
    Phi_m = phi_mat(P, P['xm'])
    Nt = len(P['xm'])
    fx = P['fx']

    def negD(lam):
        h = Phi_m.T @ lam / Nt
        w, Ea, LZ = gibbs_from_h(P, h)
        D = lam @ fx / Nt - 0.5 * ke * lam @ lam / Nt - tau * LZ
        g = (fx - ke * lam - Phi_m @ (w * Ea)) / Nt
        return -D, -g
    res = minimize(negD, fx / ke, jac=True, method='L-BFGS-B', options=dict(maxiter=5000, gtol=1e-14, ftol=1e-16))
    lam = res.x
    w, Ea, _ = gibbs_from_h(P, Phi_m.T @ lam / Nt)
    r = fx - Phi_m @ (w * Ea)
    return 0.5 * np.mean(r ** 2), float(np.max(np.abs(r / ke - lam))), w, Ea, lam


def exact_avg_stationary(P, w, Ea, Nsub=8, eta=0.3, depth=6, tol=1e-10, maxit=150):
    """Stationary state of the exactly averaged loop, starting at μ*_τ (Anderson-mixed fixed point in f̂)."""
    ke = P['ke']
    xg = np.linspace(-2.8, 2.8, 1401)
    PhiG = phi_mat(P, xg)
    Nt = len(P['t']); M = Nt * Nsub; dt = 2 * np.pi / M
    tt = np.arange(M + 1) * dt; xr = 1.5 * np.sin(tt) + 0.5 * np.sin(2 * tt)
    th = (np.arange(M) + 0.5) * dt; xh = 1.5 * np.sin(th) + 0.5 * np.sin(2 * th)
    f = P['f']

    def per_err(fg):
        r = lambda x: f(x) - np.interp(x, xg, fg)
        e = 0.0; out = np.zeros(Nt)
        for rep in range(3):
            for k in range(M):
                if rep == 2 and k % Nsub == 0:
                    out[k // Nsub] = e
                k1 = r(xr[k] + e) - ke * e
                k2 = r(xh[k] + e + dt / 2 * k1) - ke * (e + dt / 2 * k1)
                k3 = r(xh[k] + e + dt / 2 * k2) - ke * (e + dt / 2 * k2)
                k4 = r(xr[k + 1] + e + dt * k3) - ke * (e + dt * k3)
                e = e + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return out

    def Tmap(fg):
        e = per_err(fg)
        hn = phi_mat(P, P['xm'] + e).T @ e / Nt
        w2, Ea2, _ = gibbs_from_h(P, hn)
        return PhiG @ (w2 * Ea2)

    fg = PhiG @ (w * Ea)
    X, F = [], []
    nr = np.inf
    for it in range(maxit):
        Tf = Tmap(fg)
        res = Tf - fg
        nr = float(np.max(np.abs(res)))
        X.append(fg.copy()); F.append(res.copy())
        X, F = X[-depth:], F[-depth:]
        if nr < tol:
            break
        if len(F) > 1:
            dF = np.array([F[i + 1] - F[i] for i in range(len(F) - 1)]).T
            dX = np.array([X[i + 1] - X[i] for i in range(len(X) - 1)]).T
            gam, *_ = np.linalg.lstsq(dF, res, rcond=None)
            fg = fg + eta * res - (dX + eta * dF) @ gam
        else:
            fg = fg + eta * res
    fxm = np.interp(P['xm'], xg, Tf)
    return 0.5 * np.mean((P['fx'] - fxm) ** 2), nr, it, fxm


def main():
    os.makedirs(OUT, exist_ok=True)
    lines = ["# Reference values for Theorem H and Conjecture H.6 (n = 1 experiment)", "",
             "Computed by `experiments/exp_gibbs_reference.py` (concave dual on a grid; Anderson-mixed fixed point for the "
             "exactly averaged loop).  G = diag(4, 1, 1), ϑ⁰ = 0, f = sin 2x + 0.3x, x_m = 1.5 sin t + 0.5 sin 2t.", "",
             "| quantity | (σ, τ) | k_e | grid (Nγ, Nκ, Nα, N_t) | value | residual |", "|---|---|---|---|---|---|"]
    t0 = time.time()
    for (sig, tau) in ((0.02, 0.004), (0.002, 0.0005)):
        for grid in (dict(), dict(Ng=301, Nk=301, Na=1601, ga_rng=7.0, ka_rng=7.0, al_rng=16.0, Nt=320)):
            P = setup(sigma=sig, tau=tau, **grid)
            L, resid, w, Ea, lam = solve_mu_star(P)
            g = (grid.get('Ng', 241), grid.get('Nk', 241), grid.get('Na', 1201), grid.get('Nt', 256))
            lines.append(f"| L̄(μ*_τ) | ({sig:g}, {tau:g}) | 10 | {g} | {L:.5f} | {resid:.1e} |")
            print(lines[-1], f"({time.time() - t0:.0f}s)", flush=True)
    for ke in (10.0, 20.0, 40.0):
        P = setup(sigma=0.02 * 10 / ke, tau=0.004 * 10 / ke, ke=ke, Ng=161, Nk=161, Na=801, Nt=256)
        Lstar, _, w, Ea, _ = solve_mu_star(P)
        Lavg, nr, it, fx_avg = exact_avg_stationary(P, w, Ea)
        fx_star = phi_mat(P, P['xm']) @ (w * Ea)
        d = float(np.sqrt(np.mean((fx_avg - fx_star) ** 2)))
        lines.append(f"| L̄(μ*_τ), scaling σ = 0.2/k_e, τ = 0.04/k_e | ({0.2 / ke:g}, {0.04 / ke:g}) | {ke:g} | (161, 161, 801, 256) | {Lstar:.5f} | — |")
        lines.append(f"| L̄ of the exactly averaged stationary state, same scaling | ({0.2 / ke:g}, {0.04 / ke:g}) | {ke:g} | (161, 161, 801, 256) | "
                     f"{Lavg:.5f} | fixed point {nr:.1e} after {it} iterations |")
        lines.append(f"| ‖f̂(exactly averaged) − f̂(μ*_τ)‖ in L²(ρ), and k_e × it | ({0.2 / ke:g}, {0.04 / ke:g}) | {ke:g} | (161, 161, 801, 256) | "
                     f"{d:.5f} | k_e·d = {ke * d:.4f} |")
        print(lines[-3]); print(lines[-2]); print(lines[-1], f"({time.time() - t0:.0f}s)", flush=True)
    lines += ["", "In the scaling of Conjecture H.6 (σ = σ₀/k_e, τ = τ₀/k_e) μ*_τ does not depend on k_e, so the first column of "
                  "each pair should be constant.  The offset of the exactly averaged stationary state is measured by the L²(ρ) "
                  "distance between the two learned functions (the scalar loss difference can change sign); it is the O(1/k_e) "
                  "offset that makes the fixed-k_e statement false."]
    with open(os.path.join(OUT, 'gibbs_reference.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
