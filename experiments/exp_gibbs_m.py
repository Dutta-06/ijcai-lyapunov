"""
Corollary J′.3: the Gibbs measures π_m ∝ exp(−J_m/τ₀) of Theorem J′ are chaotic, with limit the minimiser μ*_τ of Theorem H's
free energy (Sanov's theorem and Varadhan's lemma): the one-particle marginal of π_m converges to μ*_τ as m → ∞.

Instance of `exp_dither_stationary.py` (n = 1, ϑ = (α, γ), frozen width, G = diag(2, 2), σ₀ = 0.4, τ₀ = 0.1).
1. μ*_τ on the 201² grid by the damped fixed-point iteration f̂ ← (1 − η)f̂ + η∫Φ̃ dμ[f̂], μ[f̂] ∝ exp(−[−⟨(f − f̂)Φ̃⟩_θ +
   (σ₀/2)‖ϑ‖²_{G⁻¹}]/τ₀) (the first variation of F_τ vanishes at μ*_τ: μ*_τ ∝ exp(−V[μ*_τ]/τ), Theorem H(ii)).
2. π_m for m = 1, 2, 4, 8, 16, 32 by MALA (exact stationary law) with K replicas, in the 2m parameter coordinates; the
   one-particle marginal means E[α₁], E[γ₁], E[α₁²], E[γ₁²] (averaged over the m exchangeable particles), the Gibbs average of
   the occupation loss E[L̄(μ_Θ)], with standard errors from the independent replicas; m = 1 is cross-checked against the quadrature of `dither.GibbsGrid`.
3. The deviations from the μ*_τ values against m (expected O(1/m)).

Output: results/dither/gibbs_m.md
"""
import os
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(__file__))
from src.core import dither as Dt
from exp_dither_stationary import make_instance, LIM

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'dither')
MS = (1, 2, 4, 8, 16, 32)
K_REP = int(os.environ.get('GM_K', 4000))
N_STEPS = int(os.environ.get('GM_STEPS', 1500))
DT = 0.05
N_THETA = 96


def mean_field_minimiser(inst, grid, eta=0.3, iters=400, tol=1e-11):
    """μ*_τ on the grid: fixed point of f̂ ↦ ∫Φ̃ dμ[f̂].  Returns (density on the grid, f̂ on the θ-grid, residual history)."""
    AL, GA = grid.AL[..., None], grid.GA[..., None]
    a, c, _, _, _, _ = inst.phys(AL, GA)
    phi = a * np.exp(-(inst.xm - c) ** 2 / (2 * inst.s ** 2))                   # Φ̃(x_m(θ); ϑ) on grid × θ, (nα, nγ, N)
    f_orb = inst.f(inst.xm)
    leak = 0.5 * (inst.sigma0 / inst.G[0] * grid.AL ** 2 + inst.sigma0 / inst.G[1] * grid.GA ** 2)
    fh = np.zeros_like(f_orb)
    hist = []
    for it in range(iters):
        VL = -np.mean((f_orb - fh)[None, None, :] * phi, axis=-1)              # V_L[μ](ϑ) = −⟨r_μ Φ̃_ϑ⟩
        logw = -(VL + leak) / inst.tau0
        w = np.exp(logw - logw.max()); mu = w / (w.sum() * grid.cell)
        fh_new = np.sum(mu[..., None] * phi, axis=(0, 1)) * grid.cell          # ∫Φ̃ dμ on the θ-grid
        res = float(np.max(np.abs(fh_new - fh)))
        hist.append(res)
        fh = (1 - eta) * fh + eta * fh_new
        if res < tol:
            break
    return mu, fh, hist


def J_m_and_grad(inst, A, Gm, m):
    """J_m(Θ) and ∇J_m for a batch: A, Gm of shape (K, m) (α and γ of each particle).  Returns J (K,), gA, gG (K, m)."""
    a, c, a1, _, c1, _ = inst.phys(A[..., None], Gm[..., None])                  # (K, m, 1)
    d = inst.xm[None, None, :] - c
    phi = np.exp(-d ** 2 / (2 * inst.s ** 2))
    fh = np.mean(a * phi, axis=1)                                                  # (K, N)
    r = inst.f(inst.xm)[None, :] - fh                                              # residual (K, N)
    Lbar = 0.5 * np.mean(r ** 2, axis=1)
    lk = inst.sigma0 / inst.G
    J = m * Lbar + 0.5 * (lk[0] * np.sum((A - inst.theta0[0]) ** 2, axis=1) + lk[1] * np.sum((Gm - inst.theta0[1]) ** 2, axis=1))
    # ∇_{ϑ_i}(m L̄_m) = −⟨r ∇_ϑΦ̃(x_m; ϑ_i)⟩ (the 1/m of f̂ cancels the m)
    gA = -np.mean(r[:, None, :] * a1 * phi, axis=2) + lk[0] * (A - inst.theta0[0])
    gG = -np.mean(r[:, None, :] * a * phi * d / inst.s ** 2 * c1, axis=2) + lk[1] * (Gm - inst.theta0[1])
    return J, gA, gG, Lbar


def mala(inst, m, K, n_steps, dt, rng, burn_frac=0.25):
    """MALA for π_m ∝ exp(−J_m/τ₀) in the G⁻¹ geometry: proposal Θ′ = Θ − dt·G∇J + √(2τ₀dt)G^{1/2}ξ.  Returns the
    one-particle observables (averaged over the m exchangeable particles of each replica and over the steps after the
    burn-in; the mean over the K independent replicas, with the standard error from their spread) and the acceptance rate."""
    tau0, G = inst.tau0, inst.G
    A = rng.normal(0.5, 0.5, (K, m)); Gm = rng.normal(0.3, 0.5, (K, m))
    J, gA, gG, Lb = J_m_and_grad(inst, A, Gm, m)
    acc = []
    n_burn = int(burn_frac * n_steps)
    per_rep = np.zeros((5, K)); cnt = 0
    for k in range(n_steps):
        xiA, xiG = rng.standard_normal((K, m)), rng.standard_normal((K, m))
        A2 = A - dt * G[0] * gA + np.sqrt(2 * tau0 * dt * G[0]) * xiA
        G2 = Gm - dt * G[1] * gG + np.sqrt(2 * tau0 * dt * G[1]) * xiG
        J2, gA2, gG2, Lb2 = J_m_and_grad(inst, A2, G2, m)
        # log q(Θ|Θ′) − log q(Θ′|Θ) with q(·|Θ) = N(Θ − dt G∇J(Θ), 2τ₀dt G)
        fwd = (np.sum((A2 - A + dt * G[0] * gA) ** 2, axis=1) / G[0] + np.sum((G2 - Gm + dt * G[1] * gG) ** 2, axis=1) / G[1]) / (4 * tau0 * dt)
        bwd = (np.sum((A - A2 + dt * G[0] * gA2) ** 2, axis=1) / G[0] + np.sum((Gm - G2 + dt * G[1] * gG2) ** 2, axis=1) / G[1]) / (4 * tau0 * dt)
        log_alpha = -(J2 - J) / tau0 - bwd + fwd
        accept = np.log(rng.random(K)) < log_alpha
        A[accept] = A2[accept]; Gm[accept] = G2[accept]; J[accept] = J2[accept]; gA[accept] = gA2[accept]; gG[accept] = gG2[accept]; Lb[accept] = Lb2[accept]
        acc.append(float(accept.mean()))
        if k >= n_burn:
            per_rep += np.array([A.mean(1), Gm.mean(1), (A ** 2).mean(1), (Gm ** 2).mean(1), Lb]); cnt += 1
    per_rep /= cnt
    return per_rep.mean(1), per_rep.std(1, ddof=1) / np.sqrt(K), float(np.mean(acc[n_burn:]))


def main():
    os.makedirs(OUT, exist_ok=True)
    inst = make_instance(+1)
    inst.N_theta = N_THETA; inst.ths = np.arange(N_THETA) * inst.T0 / N_THETA; inst.xm = inst.x_ref(inst.ths); inst.xmd = inst.x_ref_dot(inst.ths)
    grid = Dt.GibbsGrid(inst, LIM, LIM, n_al=201, n_ga=201)
    mu, fh, hist = mean_field_minimiser(inst, grid)
    ex = lambda v: float(np.sum(v * mu) * grid.cell)
    star = dict(alpha=ex(grid.AL), gamma=ex(grid.GA), alpha2=ex(grid.AL ** 2), gamma2=ex(grid.GA ** 2))
    Lbar_star = 0.5 * float(np.mean((inst.f(inst.xm) - fh) ** 2))
    pi1 = dict(alpha=grid.expect(grid.AL), gamma=grid.expect(grid.GA), alpha2=grid.expect(grid.AL ** 2), gamma2=grid.expect(grid.GA ** 2), Lbar=grid.expect(grid.data['Lbar']))
    rng = np.random.default_rng(5)
    rows = []
    for m in MS:
        mean, se, acc = mala(inst, m, K_REP, N_STEPS, DT, rng)
        rows.append((m, mean, se, acc))
    names = ['alpha', 'gamma', 'alpha2', 'gamma2']
    L = ["# Corollary J′.3: the Gibbs measures π_m are chaotic with limit μ*_τ (`exp_gibbs_m.py`)", "",
         f"Instance of `exp_dither_stationary.py`; θ-grid of {N_THETA} points; μ*_τ on the 201² grid by {len(hist)} damped fixed-point iterations "
         f"(final residual {hist[-1]:.1e}); π_m by MALA with {K_REP} replicas, {N_STEPS} steps of dt = {DT:g} (first {int(0.25 * N_STEPS)} discarded), "
         "standard errors from the spread of the per-replica averages (independent replicas). The m = 1 row is cross-checked against the quadrature of `dither.GibbsGrid`.", "",
         f"μ*_τ: ∫α dμ* = {star['alpha']:+.5f}, ∫γ dμ* = {star['gamma']:+.5f}, ∫α² dμ* = {star['alpha2']:+.5f}, ∫γ² dμ* = {star['gamma2']:+.5f}, L̄(μ*_τ) = {Lbar_star:.5f}. "
         f"Quadrature of π₁: E[α] = {pi1['alpha']:+.5f}, E[γ] = {pi1['gamma']:+.5f}, E[α²] = {pi1['alpha2']:+.5f}, E[γ²] = {pi1['gamma2']:+.5f}, E[L̄] = {pi1['Lbar']:.5f}.", "",
         "| m | E[α₁] ± s.e. | E[γ₁] | E[α₁²] | E[γ₁²] | E[L̄(μ_Θ)] | acceptance | (E[α₁] − ∫α dμ*)·m | (E[γ₁] − ∫γ dμ*)·m | (E[L̄(μ_Θ)] − L̄(μ*))·m |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for m, mean, se, acc in rows:
        L.append(f"| {m} | {mean[0]:+.5f} ± {se[0]:.5f} | {mean[1]:+.5f} ± {se[1]:.5f} | {mean[2]:+.5f} ± {se[2]:.5f} | {mean[3]:+.5f} ± {se[3]:.5f} | {mean[4]:.5f} ± {se[4]:.5f} | {acc:.3f} | "
                 f"{(mean[0] - star['alpha']) * m:+.4f} | {(mean[1] - star['gamma']) * m:+.4f} | {(mean[4] - Lbar_star) * m:+.4f} |")
    dev = np.array([[abs(r[1][0] - star['alpha']), abs(r[1][1] - star['gamma']), abs(r[1][4] - Lbar_star)] for r in rows])
    ms = np.array([r[0] for r in rows], float)
    sl = [np.polyfit(np.log(ms[1:]), np.log(dev[1:, k] + 1e-16), 1)[0] for k in range(3)]
    L += ["", f"Slopes of log|E[·] − (μ*_τ value)| against log m over m = 2–32: {sl[0]:.2f} (α), {sl[1]:.2f} (γ), {sl[2]:.2f} (L̄) (Corollary J′.3: the deviations vanish; "
          "the standard rate for mean-field Gibbs measures is 1/m, slope −1). The last three columns are the rescaled deviations; a constant means a clean 1/m law. "
          "The Gibbs average of the occupation loss is the cleanest observable: it is a smooth functional of the empirical measure and its 1/m correction is the "
          "finite-size (fluctuation) energy of the Gibbs measure.", "",
          "Numbers are produced by `experiments/exp_gibbs_m.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'gibbs_m.md'), 'w', encoding='utf-8') as fh_:
        fh_.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == '__main__':
    main()
