"""
Proposition D‴.46: a leak in the mirror coordinates — the free mass–amplitude law made robust, its bias a Bregman-regularised fit (a ridge
from large starting masses, a LASSO from small ones, which selects the target's atoms under the irrepresentable condition)
====================================================================================================================================

The free law of Proposition D‴.45 with the leak −σ(η_i − η_p,i) added to each particle's scalar drive (Q_i kept, η explicit):
dα_i = Gβ_ia′_i[g_i − σ(η_i − η_p,i)], dβ_i = γβ_ia_i[g_i − σ(η_i − η_p,i)] (`mirror_free.FreeLaw.real_rhs_leak`, `leak_equilibrium`,
`filtered_gram`).

A. Identities: V̇ = −k_ee² + ed − n(f − f̂) − σΣ[B(w*, w) + B(w, w_p) − B(w*, w_p)] in the real loop with a disturbance d and a
   measurement noise n in the adaptation, at random states; the averaged equilibrium against the ridge (an O(1) ridge effect, large starting
   masses, and at the saturation's wall) and the LASSO (small masses) on random positive definite H, k_e = 2.
B. (1) The real loop on a persistently exciting orbit under an approximation error and a biased disturbance (six Gaussian atoms,
   x_m = sin t + 0.4 sin 2.7t): without the leak the coefficients converge slowly to the biased fit w* + H_L⁻¹d̄_L of the averaged loop
   with the error filter's Gram (its slowest time constant given), with σ = 0.01, 0.1 they settle at once; V against the comparison bound
   V(0)e^{−λt} + (C/λ)(1 − e^{−λt}), λ = min(k_e, σ), C = d̄²/(2k_e) + σB(w*, w_p), asserted along the runs. (2) A drift without
   persistent excitation: the error model ė = −k_ee + (w* − w)φ(t) + d(t), φ(t) = (1 + t)^{−1/2}, one free particle, a bounded d built so
   that without the leak η = η₀ + κ log(1 + t) and w grows like (1 + t)^{1/4} without bound; the leak bounds it.
C. Slow adaptation: the real loop started at the averaged equilibrium w_σ^L (the filtered Gram), its period means (over 20π) against w_σ^L
   and against the idealised w_σ (H/k_e), the last periods' drift, the within-period oscillation; a lazy and an intermediate start.
D. On an orbit's Gram (twelve atoms of width 0.25): (1) the regularisation paths for one truth (the ridge from a lazy start, the LASSO from
   a rich one); (2) every 3-sparse support and sign pattern: the irrepresentable quantity ‖H_{S^cS}H_SS⁻¹s‖∞, the rich equilibrium's support
   at μ = 10⁻³ against the LASSO's prediction (the truth's support exactly when the quantity is < 1 and no sign flips); (3) a biased
   disturbance correlated with the regressors: w_σ → w* + H⁻¹d̄ as σ → 0; (4) the averaged loop's time from a rich start to its
   equilibrium.
E. The 1760 patterns of D(2) against the weighted LASSO itself (round 141): per pattern the LASSO on the same Gram, solved exactly
   (coordinate descent, then the active-set equations on its support, KKT checked), the equilibrium–LASSO distance, the LASSO's exact
   support and its support under D's threshold (|w| > 10⁻³ max|w|) against the truth's, the prediction and the equilibrium's; at the
   LASSO's inactive coordinates the equilibrium against the core's tangent profile at the LASSO's relative slack.

Output: results/chaos/mirror_leak.md and mirror_leak_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D,E`).
"""
import itertools
import os
import sys
import time
import warnings
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mirror_free as MF
import exp_cohesion_leak as E

OUT = os.environ.get('ML_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
KE, M6, G6, GAM6, AMAX = 5.0, 6.0, 1.0, 0.5, 1.0
XR = (lambda t: np.sin(t) + 0.4 * np.sin(2.7 * t))
XRD = (lambda t: np.cos(t) + 1.08 * np.cos(2.7 * t))
PERIOD = 20 * np.pi                                                             # x_m's period (2.7 = 27/10)
CEN6, WID6 = np.linspace(-1.2, 1.2, 6), np.full(6, 0.4)
SIGMAS_B, T_B, T_B0 = (0.0, 0.01, 0.1), 3000.0, 2e4                            # T_B0: the run without the leak
DIST_B = (lambda t: 0.2 * np.sin(3.1 * t) + 0.1)                             # a disturbance with a bias
APPROX_B = 0.1                                                                  # f = Σw*φ + 0.1 sin(4x): not representable
LAW_B2, KE_B2, ETA0_B2, T_B2, SIGMAS_B2 = dict(m=1.0, G=1.0, gamma=0.5, a_max=1.0), 1.0, 0.0, 1e6, (0.0, 0.01, 0.1)
KE_C, G_C, GAM_C = 20.0, 0.5, 0.05
CASES_C = (('lazy', 3.0, 5e-4, 2e4), ('intermediate', 0.3, 1e-4, 6e4))         # (start, β⁰, σ, T): T ≥ 6/σ
D12, WID12, SIGMAS_D, MU_IC = 12, 0.25, (1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2), 1e-3


def phi_of(x, cen, wid):
    return np.exp(-(np.asarray(x)[..., None] - cen) ** 2 / (2 * wid ** 2))


def orbit_gram(cen, wid, T=PERIOD, n=200000):
    """H = ⟨φφᵀ⟩ over one period of x_m (a uniform grid, endpoint excluded), and the grid's φ."""
    ts = np.arange(n) * T / n
    P = phi_of(XR(ts), cen, wid)
    return P.T @ P / n, ts, P


def comparison(V, ts, ke, sig, dbar, Bp):
    """V(t) against V(0)e^{−λt} + (C/λ)(1 − e^{−λt}), λ = min(k_e, σ), C = d̄²/(2k_e) + σB(w*, w_p): the largest ratio."""
    lam = min(ke, sig)
    C = dbar ** 2 / (2 * ke) + sig * Bp
    bound = V[0] * np.exp(-lam * ts) + (C / lam) * (1 - np.exp(-lam * ts))
    return float((V / bound).max()), float(C / lam)


def part_A():
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(M6, G6, GAM6, AMAX)
    rng = np.random.default_rng(4)
    al0, b0 = 0.1 * rng.standard_normal(6), rng.uniform(0.5, 1.5, 6); q = L.Q(al0, b0)
    wstar = rng.uniform(-0.15, 0.15, 6); al_star = L.alpha_of_w(wstar, q)
    eta_p = L.eta(al0, q)
    ftrue = lambda x: phi_of(x, CEN6, WID6) @ wstar
    noise = lambda t: 0.05 * np.cos(7.3 * t + 0.4)
    errs = []
    for sig in (0.0, 0.3):
        F = L.real_rhs_leak(CEN6, WID6, ftrue, XR, XRD, KE, q, eta_p, sig, dist=DIST_B, noise=noise)
        for _ in range(6):
            al = al0 + 0.3 * rng.standard_normal(6); z = np.r_[rng.uniform(-1, 1), al, np.log(L.beta(al, q))]; t = rng.uniform(0, 6)
            dz = F(t, z); hh = 1e-6
            V = lambda s: L.real_lyapunov(z + s * dz, t + s, XR, al_star, q)
            dV = (V(hh) - V(-hh)) / (2 * hh)
            e = z[0] - XR(t)
            B = lambda a1, a2: L.bregman_div(a1, a2, q)
            ff = phi_of(z[0], CEN6, WID6) @ (wstar - L.w(al, q))
            pred = -KE * e ** 2 + e * DIST_B(t) - noise(t) * ff - sig * (B(al_star, al) + B(al, al0) - B(al_star, al0))
            errs.append(abs(dV - pred) / max(1.0, abs(pred)))
    Lz = MF.FreeLaw(12.0, 1.0, 0.1, 1.0)
    ke, rho = 2.0, 0.5                                                          # the ridge parameter ρ = k_eσm/(Gβ⁰²): an O(1) effect
    lim = []
    for rep in range(4):
        A = rng.standard_normal((30, 10)); H = A.T @ A / 30
        ws = np.zeros(10); ws[rng.choice(10, 2, replace=False)] = rng.choice([-1, 1], 2) * rng.uniform(0.15, 0.3, 2)
        for b0 in (1e4, 1e2, 3.0):
            q = Lz.Q(np.zeros(10), np.full(10, b0))
            sig = rho * Lz.G * b0 ** 2 / (ke * Lz.m)
            w, _, it, res = Lz.leak_equilibrium(H, ws, q, np.zeros(10), sig, ke)
            ridge = np.linalg.solve(H + rho * np.eye(10), H @ ws)
            lim.append(dict(kind='ridge', b0=b0, d=float(np.abs(w - ridge).max() / np.abs(ridge).max()), res=res,
                            effect=float(np.abs(ridge - ws).max() / np.abs(ws).max()), wall=float(b0 * Lz.A / Lz.m),
                            wmax=float(np.abs(ws).max())))
        for b0 in (1e-4, 1e-6, 1e-8):
            q = Lz.Q(np.zeros(10), np.full(10, b0))
            kap = np.sqrt((2 * Lz.c - b0) / b0); lam1 = np.pi * Lz.A * kap / (4 * Lz.G * Lz.c)
            sig = 0.02 / lam1                                                     # the LASSO penalty σλ₁ = 0.02 at every β⁰
            w, _, it, res = Lz.leak_equilibrium(H, ws, q, np.zeros(10), sig, ke)
            ls = MF.lasso(H, ws, np.full(10, sig * lam1), ke)
            lim.append(dict(kind='lasso', b0=b0, d=float(np.abs(w - ls).max() / np.abs(ls).max()), res=res,
                            same_support=bool(np.array_equal(np.abs(w) > 1e-3 * np.abs(w).max(), np.abs(ls) > 1e-12)),
                            lam1_lead=float(np.pi * Lz.A / (4 * Lz.G) * np.sqrt(2 / (Lz.c * b0)) / lam1 - 1)))
    return dict(dV=float(max(errs)), lim=lim, ke=ke, rho=rho)


def part_B():
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(M6, G6, GAM6, AMAX)
    rng = np.random.default_rng(11)
    al0, b0 = 0.1 * rng.standard_normal(6), rng.uniform(0.5, 1.5, 6); q = L.Q(al0, b0)
    wstar = rng.uniform(-0.15, 0.15, 6)
    ftrue = lambda x: phi_of(x, CEN6, WID6) @ wstar + APPROX_B * np.sin(4 * x)
    eta_p = L.eta(al0, q)
    al_star = L.alpha_of_w(wstar, q)
    dbar = 0.2 + 0.1 + APPROX_B
    Bp = float(L.bregman_div(al_star, al0, q))
    # the averaged loop with the error filter: H_L, d̄_L over one period of x_m (the approximation error and the disturbance as d)
    _, ts_o, P = orbit_gram(CEN6, WID6, n=20000)
    HL, dL = MF.filtered_gram(P, PERIOD, KE, APPROX_B * np.sin(4 * XR(ts_o)) + DIST_B(ts_o))
    w_lim = wstar + np.linalg.solve(HL, dL)                                       # σ = 0: the biased fit
    al_lim = L.alpha_of_w(w_lim, q)
    rates = np.abs(np.linalg.eigvals(HL * (L.m * L.relevance(al_lim, q))[None, :] / KE))
    rows = []
    for sig in SIGMAS_B:
        F = L.real_rhs_leak(CEN6, WID6, ftrue, XR, XRD, KE, q, eta_p, sig, dist=DIST_B)
        T = T_B0 if sig == 0 else T_B
        z0 = np.r_[XR(0.0) + 0.3, al0, np.log(b0)]
        ts = np.linspace(0, T, int(T) + 1)
        sol = solve_ivp(F, (0, T), z0, method='DOP853', rtol=1e-10, atol=1e-12, t_eval=ts)
        Z = sol.y
        W = np.exp(Z[7:]) * L.amp(Z[1:7])[0] / L.m
        V = np.array([L.real_lyapunov(Z[:, j], t, XR, al_star, q) for j, t in enumerate(ts)])
        e = Z[0] - XR(ts)
        qs = np.array([L.Q(Z[1:7, j], np.exp(Z[7:, j])) for j in range(len(ts))])
        n4 = len(ts) // 4
        quart = [slice(k * n4, (k + 1) * n4 + 1) for k in range(4)]
        row = dict(sig=sig, T=T, ok=bool(sol.success), wmax=[float(np.abs(W[:, s]).max()) for s in quart],
                   Vmax=[float(V[s].max()) for s in quart], erms=[float(np.sqrt(np.mean(e[s] ** 2))) for s in quart],
                   dQ=float(np.abs(qs - q).max()), V0=float(V[0]), Bp=Bp)
        if sig > 0:
            w_eq, _, it, res = L.leak_equilibrium(HL, wstar, q, eta_p, sig, KE, dphi=dL)
            ratio, ult = comparison(V, ts, KE, sig, dbar, Bp)
            assert ratio <= 1 + 1e-6, ('the comparison bound is violated', sig, ratio)
            late = W[:, ts > 0.75 * T].mean(axis=1)
            row.update(ratio=ratio, ultimate=ult, w_eq_max=float(np.abs(w_eq).max()), d_eq=float(np.abs(late - w_eq).max()), eq_res=res)
        rows.append(row)
    b1 = dict(rows=rows, dbar=dbar, ke=KE, w_lim_max=float(np.abs(w_lim).max()), tau_slow=float(1 / rates.min()),
              eigH=np.linalg.eigvalsh(0.5 * (HL + HL.T)).tolist(),
              V_lim=float(0.5 * 0 + L.bregman_div(al_star, al_lim, q)))
    # (2) a drift without persistent excitation (the error model with one free particle)
    L2 = MF.FreeLaw(**LAW_B2)
    kap = 1.0 / (4 * L2.gamma * L2.A)                                             # w ~ (1 + t)^{γa_max κ} = (1 + t)^{1/4}
    alp0, bp0 = 0.5, 1.0
    q2 = float(L2.Q(alp0, bp0))
    e0 = float(L2.eta(alp0, q2)[0])
    phi = lambda t: (1 + t) ** -0.5
    w_of_t = lambda t: float(L2.w(L2.alpha_of_eta(e0 + kap * np.log1p(t), q2)[0], q2))
    dist2 = lambda t: -0.5 * kap * (1 + t) ** -1.5 + KE_B2 * kap * (1 + t) ** -0.5 + w_of_t(t) * phi(t)      # w* = 0
    tg = np.r_[0.0, np.geomspace(1e-2, T_B2, 400)]
    dbar2 = float(max(abs(dist2(t)) for t in tg))
    rows2 = []
    for sig in SIGMAS_B2:
        def F2(t, z, sig=sig):
            e, al, lb = z
            b = np.exp(lb); a, da = L2.amp(al)
            w = b * a / L2.m
            u = e * phi(t) - (sig * (float(L2.eta(al, L2.Q(al, b))[0]) - e0) if sig else 0.0)
            return [-KE_B2 * e - w * phi(t) + dist2(t), L2.G * b * da * u, L2.gamma * a * u]
        ts2 = np.r_[0.0, np.geomspace(1.0, T_B2, 25)]
        sol = solve_ivp(F2, (0, T_B2), [kap, alp0, np.log(bp0)], method='LSODA', rtol=1e-10, atol=1e-12, t_eval=ts2)
        W2 = np.exp(sol.y[2]) * L2.amp(sol.y[1])[0] / L2.m
        V2 = 0.5 * sol.y[0] ** 2 + np.array([L2.bregman_div(np.zeros(1), np.array([a_]), q2) for a_ in sol.y[1]])
        row = dict(sig=sig, t=sol.t.tolist(), w=W2.tolist(), ok=bool(sol.success), V=V2.tolist())
        if sig > 0:
            Bp2 = float(L2.bregman_div(np.zeros(1), np.array([alp0]), q2))
            ratio, ult = comparison(V2, sol.t, KE_B2, sig, dbar2, Bp2)
            assert ratio <= 1 + 1e-6, ('the comparison bound is violated (B2)', sig, ratio)
            row.update(ratio=ratio, ultimate=ult)
        rows2.append(row)
    b2 = dict(rows=rows2, dbar=dbar2, kap=kap, w_pred=[w_of_t(t) for t in rows2[0]['t']])
    return dict(b1=b1, b2=b2)


def part_C():
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(M6, G_C, GAM_C, AMAX)
    rng = np.random.default_rng(13)
    wstar = rng.uniform(-0.15, 0.15, 6)
    ftrue = lambda x: phi_of(x, CEN6, WID6) @ wstar
    H, ts_o, P = orbit_gram(CEN6, WID6, n=20000)
    HL = MF.filtered_gram(P, PERIOD, KE_C)
    rows = []
    for lab, b0v, sig, T in CASES_C:
        al0 = np.zeros(6); b0 = np.full(6, b0v); q = L.Q(al0, b0); eta_p = L.eta(al0, q)
        w_id, _, _, res_id = L.leak_equilibrium(H, wstar, q, eta_p, sig, KE_C)
        w_L, al_L, _, res_L = L.leak_equilibrium(HL, wstar, q, eta_p, sig, KE_C)
        F = L.real_rhs_leak(CEN6, WID6, ftrue, XR, XRD, KE_C, q, eta_p, sig)
        n_per = int(T // PERIOD)
        T = n_per * PERIOD
        ts = np.linspace(0, T, n_per * 400 + 1)
        z0 = np.r_[XR(0.0), al_L, np.log(L.beta(al_L, q))]
        sol = solve_ivp(F, (0, T), z0, method='LSODA', rtol=1e-9, atol=1e-12, t_eval=ts)
        W = np.exp(sol.y[7:]) * L.amp(sol.y[1:7])[0] / L.m
        means = np.array([W[:, k * 400:(k + 1) * 400].mean(axis=1) for k in range(n_per)])
        last = means[-1]
        rows.append(dict(lab=lab, b0=b0v, sig=sig, T=float(T), ok=bool(sol.success), n_per=n_per,
                         d_L=float(np.abs(last - w_L).max()), d_id=float(np.abs(last - w_id).max()),
                         lag=float(np.abs(w_L - w_id).max()), bias=float(np.abs(w_L - wstar).max()),
                         drift=float(np.abs(means[-1] - means[-2]).max()), osc=float((W[:, -400:].max(axis=1) - W[:, -400:].min(axis=1)).max()),
                         d_L_first=float(np.abs(means[0] - w_L).max()), res=max(res_id, res_L)))
    return dict(rows=rows, HL_rel=float(np.abs(HL - H).max() / np.abs(H).max()), HL_asym=float(np.abs(HL - HL.T).max()))


def part_D():
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(12.0, 1.0, 0.1, AMAX)
    cen = np.linspace(-1.4, 1.4, D12); wid = np.full(D12, WID12)
    H, _, _ = orbit_gram(cen, wid)
    ws = np.zeros(D12); idx = np.array([2, 6, 9]); ws[idx] = [0.25, -0.2, 0.3]
    out = dict(eig=np.linalg.eigvalsh(H).tolist(), truth=idx.tolist(), rows=[],
               ic_truth=MF.irrepresentable(H, idx, np.sign(ws[idx])))
    lam1_of = lambda b0: np.pi * L.A * np.sqrt(abs(2 * L.c - b0) / b0) / (4 * L.G * L.c)
    for lab, b0 in (('lazy', 30.0), ('rich', 1e-6)):
        q = L.Q(np.zeros(D12), np.full(D12, b0))
        for mu in SIGMAS_D:
            sig = mu / lam1_of(b0) if lab == 'rich' else mu
            w, _, it, res = L.leak_equilibrium(H, ws, q, np.zeros(D12), sig, 1.0)
            if lab == 'rich':
                ref = MF.lasso(H, ws, np.full(D12, mu), 1.0)
            else:
                ref = np.linalg.solve(H + sig * L.m / (L.G * b0 ** 2) * np.eye(D12), H @ ws)
            sup = np.flatnonzero(np.abs(w) > 1e-3 * max(np.abs(w).max(), 1e-12))
            out['rows'].append(dict(lab=lab, mu=mu, d=float(np.abs(w - ref).max()), sup=sup.tolist(), res=res,
                                    truth_kept=bool(set(idx) <= set(sup)), extra=int(len(set(sup) - set(idx))),
                                    err=float(np.abs(w - ws).max())))
    # (2) every 3-sparse support and sign pattern: the irrepresentable quantity and the rich equilibrium's support at μ = MU_IC
    b0 = 1e-6
    q = L.Q(np.zeros(D12), np.full(D12, b0)); sig = MU_IC / lam1_of(b0)
    rng = np.random.default_rng(23)
    pats = []
    for S in itertools.combinations(range(D12), 3):
        for sg in itertools.product((1.0, -1.0), repeat=3):
            S_ = np.array(S); sg_ = np.array(sg)
            mag = rng.uniform(0.15, 0.3, 3)
            wt = np.zeros(D12); wt[S_] = sg_ * mag
            ic = MF.irrepresentable(H, S_, sg_)
            wS = wt[S_] - MU_IC * np.linalg.solve(H[np.ix_(S_, S_)], sg_)            # the LASSO on S (k_e = 1)
            pred = bool(ic < 1 and np.all(np.sign(wS) == sg_))
            w, _, it, res = L.leak_equilibrium(H, wt, q, np.zeros(D12), sig, 1.0)
            sup = np.flatnonzero(np.abs(w) > 1e-3 * np.abs(w).max())
            pats.append(dict(ic=float(ic), pred=pred, exact=bool(set(sup) == set(S)), extra=int(len(set(sup) - set(S))),
                             missing=int(len(set(S) - set(sup))), res=res))
    out['ic'] = dict(n=len(pats), ic_lt1=int(sum(p['ic'] < 1 for p in pats)), pred=int(sum(p['pred'] for p in pats)),
                     exact=int(sum(p['exact'] for p in pats)), agree=int(sum(p['pred'] == p['exact'] for p in pats)),
                     extra_when_fail=float(np.mean([p['extra'] for p in pats if not p['exact']])) if any(not p['exact'] for p in pats) else 0.0,
                     max_res=float(max(p['res'] for p in pats)))
    # (3) a biased disturbance: w_σ → w* + H⁻¹d̄ as σ → 0 (a lazy start)
    v = 0.05 * np.random.default_rng(29).standard_normal(D12)
    dbar = H @ v                                                                  # a disturbance correlated with the regressors: H⁻¹d̄ = v
    qz = L.Q(np.zeros(D12), np.full(D12, 30.0))
    shift = ws + v
    out['bias'] = []
    for sig in (1e-2, 1e-4, 1e-6, 1e-8):
        w, _, it, res = L.leak_equilibrium(H, ws, qz, np.zeros(D12), sig, 1.0, dphi=dbar)
        out['bias'].append(dict(sig=sig, d_star=float(np.abs(w - ws).max()), d_shift=float(np.abs(w - shift).max()), res=res))
    # (4) the averaged loop from a rich start (β⁰ = 10⁻⁶, μ = 10⁻³): the time to its equilibrium, for the truth above and for a pattern
    # failing the irrepresentable condition
    out['dyn'] = []
    for lab, S_, sg_ in (('the truth above', idx, np.sign(ws[idx])), ('(0, 5, 11; +, −, +)', np.array([0, 5, 11]), np.array([1.0, -1.0, 1.0]))):
        wt = np.zeros(D12); wt[S_] = sg_ * np.array([0.25, 0.2, 0.3])
        w_eq, _, _, res = L.leak_equilibrium(H, wt, q, np.zeros(D12), sig, 1.0)

        def Fa(t, al):
            b = L.beta(al, q); w = b * L.amp(al)[0] / L.m
            return L.G * b * L.amp(al)[1] * (H @ (wt - w) - sig * L.eta(al, q))
        ts = np.geomspace(1.0, 1e8, 81)
        sol = solve_ivp(Fa, (0, 1e8), np.zeros(D12), method='LSODA', rtol=1e-9, atol=1e-13, t_eval=ts)
        dist = [float(np.abs(L.w(sol.y[:, j], q) - w_eq).max() / np.abs(w_eq).max()) for j in range(sol.y.shape[1])]
        reach = next((t for t, d_ in zip(sol.t, dist) if d_ < 1e-3), None)
        out['dyn'].append(dict(lab=lab, ic=MF.irrepresentable(H, S_, sg_), reach=reach, ok=bool(sol.success),
                               sup=[int(j) for j in np.flatnonzero(np.abs(w_eq) > 1e-3 * np.abs(w_eq).max())]))
    return out


def lasso_exact(H, ws, mu, ke, it_max=60):
    """The LASSO min (1/(2k_e))(w − w*)ᵀH(w − w*) + μ|w|₁ (H positive definite, so the minimiser is unique): coordinate descent
    (`mirror_free.lasso`), then the active-set equations H_SS w_S = (Hw*)_S − k_eμs_S on the current support and signs, the support
    corrected (a sign flip leaves it, the largest KKT violator enters) until the KKT conditions hold. Returns (w, the KKT residual
    relative to k_eμ, the LASSO's relative slacks 1 − |(H(w* − w))_j|/(k_eμ) on its inactive coordinates, NaN on the support)."""
    H, ws = np.asarray(H, float), np.asarray(ws, float)
    n = len(ws)
    w = MF.lasso(H, ws, np.full(n, mu), ke)
    S = set(np.flatnonzero(w).tolist())
    sg = np.sign(w)
    for _ in range(it_max):
        Sl = np.array(sorted(S), int)
        w = np.zeros(n)
        if Sl.size:
            w[Sl] = np.linalg.solve(H[np.ix_(Sl, Sl)], (H @ ws)[Sl] - ke * mu * sg[Sl])
        flip = [j for j in Sl if np.sign(w[j]) != sg[j]]
        if flip:
            S -= set(flip)
            continue
        c = H @ (ws - w) / (ke * mu)                                              # the KKT multipliers: s on S, |c| ≤ 1 off S
        off = np.array([j for j in range(n) if j not in S], int)
        if off.size and np.abs(c[off]).max() > 1 + 1e-12:
            j = int(off[np.argmax(np.abs(c[off]))])
            S.add(j); sg[j] = np.sign(c[j])
            continue
        break
    c = H @ (ws - w) / (ke * mu)
    on = np.flatnonzero(w)
    kkt = max(float(np.abs(c[on] - np.sign(w[on])).max()) if on.size else 0.0,
              float(max(np.abs(c[w == 0]).max() - 1, 0.0)) if (w == 0).any() else 0.0)
    slack = np.where(w == 0, 1 - np.abs(c), np.nan)
    return w, kkt, slack


def tangent_profile(L, b0, slack, sgn):
    """The core's tangent profile (Derived statement 4.4 of the paper; the rich limit of the leak's equilibrium): an inactive coordinate
    whose mirror coordinate sits at (1 − μ)λ₁ has w = sign·a_max β⁰^{3/2}/(m(2c)^{1/2})·cos(πμ/2)/sin³(πμ/2)."""
    x = np.pi * np.asarray(slack, float) / 2
    return sgn * L.A * b0 ** 1.5 / (L.m * np.sqrt(2 * L.c)) * np.cos(x) / np.sin(x) ** 3


def pattern_record(L, H, wt, S_, sg_, q, sig, mu, b0, thr=1e-3):
    """One pattern of Part E: the leak's equilibrium (rich start, penalty μ = σλ₁) against the exact LASSO of penalty μ (k_e = 1)."""
    n = len(wt)
    ic = MF.irrepresentable(H, S_, sg_)
    wS = wt[S_] - mu * np.linalg.solve(H[np.ix_(S_, S_)], sg_)
    pred = bool(ic < 1 and np.all(np.sign(wS) == sg_))
    w, _, it, res = L.leak_equilibrium(H, wt, q, np.zeros(n), sig, 1.0)
    wl, kkt, slack = lasso_exact(H, wt, mu, 1.0)
    sup_eq = set(np.flatnonzero(np.abs(w) > thr * np.abs(w).max()).tolist())
    sup_l0 = set(np.flatnonzero(wl).tolist())
    sup_lt = set(np.flatnonzero(np.abs(wl) > thr * np.abs(wl).max()).tolist())
    S = set(int(i) for i in S_)
    inact = np.flatnonzero(wl == 0)
    c = H @ (wt - wl) / mu
    tp = tangent_profile(L, b0, slack[inact], np.sign(c[inact]))
    tp_err = np.abs(w[inact] - tp) / np.abs(tp)
    return dict(ic=float(ic), pred=pred, res=res, kkt=kkt,
                eq_exact=sup_eq == S, l0_exact=sup_l0 == S, lt_exact=sup_lt == S,
                eq_l0=sup_eq == sup_l0, eq_lt=sup_eq == sup_lt, lt_l0=sup_lt == sup_l0,
                eq_extra_vs_l0=int(len(sup_eq - sup_l0)), eq_missing_vs_l0=int(len(sup_l0 - sup_eq)),
                dist=float(np.abs(w - wl).max()), rel=float(np.abs(w - wl).max() / np.abs(wl).max()),
                n_l0=len(sup_l0), slack_min=float(np.nanmin(slack)) if inact.size else float('nan'),
                l_small=float(np.abs(wl[wl != 0]).min() / np.abs(wl).max()),
                eq_off_max=float(np.abs(w[inact]).max() / np.abs(w).max()) if inact.size else 0.0,
                tp_off_max=float(np.abs(tp).max() / np.abs(w).max()) if inact.size else 0.0,
                tp_err_med=float(np.median(tp_err)) if inact.size else float('nan'),
                tp_err_max=float(tp_err.max()) if inact.size else float('nan'))


def part_E():
    """D(2)'s 1760 patterns (the same Gram, magnitudes and rich start; rng seed 23) against the weighted LASSO itself."""
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(12.0, 1.0, 0.1, AMAX)
    cen = np.linspace(-1.4, 1.4, D12); wid = np.full(D12, WID12)
    H, _, _ = orbit_gram(cen, wid)
    b0 = 1e-6
    lam1 = np.pi * L.A * np.sqrt(abs(2 * L.c - b0) / b0) / (4 * L.G * L.c)
    q = L.Q(np.zeros(D12), np.full(D12, b0)); sig = MU_IC / lam1
    rng = np.random.default_rng(23)
    pats = []
    for S in itertools.combinations(range(D12), 3):
        for sg in itertools.product((1.0, -1.0), repeat=3):
            S_ = np.array(S); sg_ = np.array(sg)
            mag = rng.uniform(0.15, 0.3, 3)
            wt = np.zeros(D12); wt[S_] = sg_ * mag
            r = pattern_record(L, H, wt, S_, sg_, q, sig, MU_IC, b0)
            r.update(S=list(S), sg=list(sg))
            pats.append(r)
    return dict(pats=pats, b0=b0, mu=MU_IC, core=float(L.A * b0 ** 1.5 / (L.m * np.sqrt(2 * L.c))))


def report_E(Es):
    fmt = E.fmt
    P = Es['pats']
    n = len(P)
    cnt = lambda f: int(sum(1 for p in P if f(p)))
    rng_ = lambda xs: f"{fmt(min(xs))}–{fmt(max(xs))}" if xs else "—"
    med = lambda xs: fmt(float(np.median(xs))) if xs else "—"
    mm = [p for p in P if p['pred'] != p['eq_exact']]
    mt, mf = [p for p in mm if p['pred']], [p for p in mm if not p['pred']]
    eqd = [p for p in P if not p['eq_l0']]
    ltd = [p for p in P if not p['lt_l0']]
    L_ = [f"## E. The patterns of D(2) against the weighted LASSO itself ({n} patterns; the same Gram, magnitudes and rich start, "
          f"β⁰ = {Es['b0']:g}, μ = σλ₁ = {Es['mu']:g}, k_e = 1)", "",
          "The LASSO min ½(w − w*)ᵀH(w − w*) + μ|w|₁ solved exactly (coordinate descent, then the active-set equations on its support, "
          "the KKT conditions checked); its exact support (w_j ≠ 0) and its support under D's threshold (|w| > 10⁻³ max|w|) against the "
          "truth's, the prediction of D(2) and the equilibrium's (threshold support, as in D):", "",
          "| quantity | count |", "|---|---|",
          f"| the LASSO's exact support is the truth's | {cnt(lambda p: p['l0_exact'])} |",
          f"| — and agrees with the prediction | {cnt(lambda p: p['l0_exact'] == p['pred'])} of {n} |",
          f"| the LASSO's threshold support is the truth's | {cnt(lambda p: p['lt_exact'])} |",
          f"| — and agrees with the prediction | {cnt(lambda p: p['lt_exact'] == p['pred'])} of {n} |",
          f"| the LASSO's threshold support equals its exact support | {cnt(lambda p: p['lt_l0'])} of {n} |",
          f"| the equilibrium's support is the truth's (D(2)) | {cnt(lambda p: p['eq_exact'])} |",
          f"| — and agrees with the prediction (D(2)) | {cnt(lambda p: p['eq_exact'] == p['pred'])} of {n} |",
          f"| the equilibrium's support equals the LASSO's exact support | {cnt(lambda p: p['eq_l0'])} of {n} |",
          f"| the equilibrium's support equals the LASSO's threshold support | {cnt(lambda p: p['eq_lt'])} of {n} |", "",
          f"The equilibrium against the LASSO, max|w_eq − w_LASSO| over the patterns: median {med([p['dist'] for p in P])}, largest "
          f"{fmt(max(p['dist'] for p in P))} (relative to max|w_LASSO|: median {med([p['rel'] for p in P])}, largest "
          f"{fmt(max(p['rel'] for p in P))}); the LASSO's largest KKT residual {fmt(max(p['kkt'] for p in P))} (relative to μ), the "
          f"equilibria's largest solver residual {fmt(max(p['res'] for p in P))}.", "",
          f"Where the LASSO's threshold support differs from its exact support ({len(ltd)} patterns): the LASSO's smallest nonzero coefficient "
          f"relative to its largest {rng_([p['l_small'] for p in ltd])}; irrepresentable quantity {min(p['ic'] for p in ltd):.4f}–{max(p['ic'] for p in ltd):.4f}.", "",
          f"Where the equilibrium's support differs from the LASSO's exact support ({len(eqd)} patterns): atoms the equilibrium keeps beyond the "
          f"LASSO's {sum(p['eq_extra_vs_l0'] for p in eqd)}, LASSO atoms it drops {sum(p['eq_missing_vs_l0'] for p in eqd)}; the LASSO's "
          f"smallest relative slack {rng_([p['slack_min'] for p in eqd])} (all patterns: median {med([p['slack_min'] for p in P if np.isfinite(p['slack_min'])])}); "
          f"the equilibrium's largest coefficient on the LASSO's inactive atoms, relative to its largest, {rng_([p['eq_off_max'] for p in eqd])}, "
          f"and the core's tangent profile there {rng_([p['tp_off_max'] for p in eqd])}; the LASSO's smallest nonzero coefficient, relative, "
          f"{rng_([p['l_small'] for p in eqd])}.", "",
          f"The equilibrium on the LASSO's inactive atoms against the core's tangent profile a_max β⁰^{{3/2}}/(m(2c)^{{1/2}})·cos(πμ_j/2)/sin³(πμ_j/2) "
          f"(a_max β⁰^{{3/2}}/(m(2c)^{{1/2}}) = {fmt(Es['core'])}) at the LASSO's relative slacks μ_j, the relative error per pattern: median of the "
          f"medians {med([p['tp_err_med'] for p in P if np.isfinite(p['tp_err_med'])])}, largest {fmt(max(p['tp_err_max'] for p in P if np.isfinite(p['tp_err_max'])))}.", "",
          f"Where the prediction and the equilibrium disagree ({len(mm)} patterns), irrepresentable quantity "
          f"{min(p['ic'] for p in mm):.4f}–{max(p['ic'] for p in mm):.4f}:", "",
          f"- the prediction holds ({len(mt)}): the LASSO's exact support is the truth's for {sum(p['l0_exact'] for p in mt)}; the equilibrium "
          f"keeps {sum(p['eq_extra_vs_l0'] for p in mt)} atoms beyond it, its largest coefficient on the LASSO's inactive atoms relative to its "
          f"largest {rng_([p['eq_off_max'] for p in mt])} (above the threshold 10⁻³), at the LASSO's smallest relative slacks "
          f"{rng_([p['slack_min'] for p in mt])};",
          f"- the prediction fails ({len(mf)}): the LASSO's exact support has {', '.join(str(n_) for n_ in sorted({p['n_l0'] for p in mf}))} atoms, "
          f"its smallest coefficient relative to its largest {rng_([p['l_small'] for p in mf])} (below the threshold), its threshold support is "
          f"the truth's for {sum(p['lt_exact'] for p in mf)}, and the equilibrium's support equals the LASSO's threshold support for "
          f"{sum(p['eq_lt'] for p in mf)}.", ""]
    return L_


def report(A=None, B=None, C=None, D=None, Es=None):
    fmt = E.fmt
    L_ = ["# Proposition D‴.46: a leak in the mirror coordinates (`exp_mirror_leak.py`)", "",
          "The free law with the leak −σ(η_i − η_p,i) on each particle's scalar drive (Q kept); η, B the mirror coordinate and Bregman "
          "divergence of Proposition D‴.45 (`mirror_free`).", ""]
    if A is not None:
        L_ += ["## A. Identities", "", f"V̇ against −k_ee² + ed − n(f − f̂) − σΣ[B(w*, w) + B(w, w_p) − B(w*, w_p)] in the real loop (k_e = {KE:g}, "
               f"a disturbance and a measurement noise in the adaptation, σ = 0 and 0.3, 12 random states): largest relative difference {fmt(A['dV'])}.", "",
               f"The averaged equilibrium (k_e = {A['ke']:g}) against its limits on random positive definite H (10 × 10, 4 draws), relative "
               "distances, the largest over the draws: the ridge with ρ = k_eσm/(Gβ⁰²) = "
               f"{A['rho']:g} (its effect max|w_ridge − w*|/max|w*| given; at β⁰ = 3 the saturation's wall w_⊓ = β⁰a_max/m lies within the "
               "coefficients' range), the LASSO with penalty σλ₁(β⁰) = 0.02 (λ₁ exact; its leading term's relative error given):", "",
               "| the limit | β⁰ | distance | the solver's residual | the ridge's effect, or the supports equal | w_⊓ against max|w*|, or λ₁'s leading term |",
               "|---|---|---|---|---|---|"]
        for kind in ('ridge', 'lasso'):
            for b0 in sorted({r['b0'] for r in A['lim'] if r['kind'] == kind}, reverse=True):
                rs = [r for r in A['lim'] if r['kind'] == kind and r['b0'] == b0]
                if kind == 'ridge':
                    L_.append(f"| ridge | {b0:g} | {fmt(max(r['d'] for r in rs))} | {fmt(max(r['res'] for r in rs))} | "
                              f"{fmt(min(r['effect'] for r in rs))}–{fmt(max(r['effect'] for r in rs))} | {fmt(rs[0]['wall'])} against "
                              f"{fmt(max(r['wmax'] for r in rs))} |")
                else:
                    L_.append(f"| LASSO | {b0:g} | {fmt(max(r['d'] for r in rs))} | {fmt(max(r['res'] for r in rs))} | "
                              f"{sum(r['same_support'] for r in rs)} of {len(rs)} | {fmt(rs[0]['lam1_lead'])} |")
        L_.append("")
    if B is not None:
        b1 = B['b1']
        L_ += [f"## B. (1) The real loop under an approximation error ({APPROX_B:g} sin 4x) and a disturbance (0.2 sin 3.1t + 0.1) (k_e = {KE:g}, "
               f"m = {M6:g}, G = {G6:g}, γ = {GAM6:g}; six atoms of width 0.4; x_m's period 20π)", "",
               "The orbit is persistently exciting: the filtered Gram H_L = k_e⟨φL[φ]ᵀ⟩ (L = 1/(s + k_e)) has eigenvalues "
               f"{fmt(min(b1['eigH']))} to {fmt(max(b1['eigH']))}. Without the leak the averaged loop converges to the biased fit "
               f"w* + H_L⁻¹d̄_L (d̄_L = k_e⟨φL[ε + d]⟩): max|w| = {fmt(b1['w_lim_max'])}, ΣB(w*, w) = {fmt(b1['V_lim'])}, its slowest time "
               f"constant {fmt(b1['tau_slow'])}. By quarters of each run: the largest |w|, the largest V (V = ½e² + ΣB(w*, w), w* the "
               "representable part), the rms error; the largest |Q − Q(0)|; for σ > 0 the late mean's distance to the averaged equilibrium "
               "w_σ^L, the comparison bound's largest ratio V/(V(0)e^{−λt} + (C/λ)(1 − e^{−λt})) (asserted ≤ 1) and the ultimate bound C/λ "
               f"(d̄ = {b1['dbar']:g}):", "",
               "| σ | T | max |w| by quarters | max V by quarters | rms e by quarters | max |Q − Q(0)| | late mean − w_σ^L | V/bound (largest) | C/λ |",
               "|---|---|---|---|---|---|---|---|---|"]
        for r in b1['rows']:
            L_.append(f"| {r['sig']:g} | {r['T']:g} | " + ', '.join(f"{x:.3f}" for x in r['wmax']) + " | " + ', '.join(fmt(x) for x in r['Vmax'])
                      + " | " + ', '.join(f"{x:.4f}" for x in r['erms']) + f" | {fmt(r['dQ'])} | "
                      + (f"{fmt(r['d_eq'])} | {r['ratio']:.4f} | {fmt(r['ultimate'])} |" if r['sig'] > 0 else "— | — | — |"))
        b2 = B['b2']
        L_ += ["", "(2) A drift without persistent excitation: ė = −k_ee + (w* − w)φ(t) + d(t), φ(t) = (1 + t)^{−1/2}, w* = 0, one free particle "
               f"(m = {LAW_B2['m']:g}, G = {LAW_B2['G']:g}, γ = {LAW_B2['gamma']:g}, α(0) = 0.5, β(0) = 1), k_e = {KE_B2:g}; d built so that "
               f"without the leak η = η(0) + κ log(1 + t), κ = {b2['kap']:g}, and w(t) grows like (1 + t)^{{1/4}}; |d| ≤ {fmt(b2['dbar'])}. "
               "w at t = 1, 10², 10⁴, 10⁶ (and the construction's w(t) for σ = 0); for σ > 0 the comparison bound's largest ratio and C/λ:", "",
               "| σ | w(1), w(10²), w(10⁴), w(10⁶) | V/bound (largest) | C/λ |", "|---|---|---|---|"]
        pick = lambda r: [int(np.argmin(np.abs(np.array(r['t']) - tt))) for tt in (1.0, 1e2, 1e4, 1e6)]
        for r in b2['rows']:
            ks = pick(r)
            L_.append(f"| {r['sig']:g} | " + ', '.join(fmt(r['w'][k]) for k in ks)
                      + (f" | {r['ratio']:.4f} | {fmt(r['ultimate'])} |" if r['sig'] > 0 else " | — | — |"))
        ks = pick(b2['rows'][0])
        L_ += [f"| the construction (σ = 0) | " + ', '.join(fmt(b2['w_pred'][k]) for k in ks) + " | | |", ""]
    if C is not None:
        L_ += [f"## C. Slow adaptation (k_e = {KE_C:g}, G = {G_C:g}, γ = {GAM_C:g}, from α = 0 and equal masses β⁰, the leak toward the start)", "",
               f"The filtered Gram against the idealised one: max|H_L − H|/max|H| = {fmt(C['HL_rel'])}, asymmetry {fmt(C['HL_asym'])}. The real loop "
               "started at the averaged equilibrium w_σ^L (with H_L), its means over x_m's periods (20π): the last period's mean against "
               "w_σ^L and against the idealised w_σ (H/k_e), the first period's, the change between the last two periods' means, the last "
               "period's peak-to-peak oscillation; the lag |w_σ^L − w_σ| and the leak's bias |w_σ^L − w*|:", "",
               "| start | β⁰ | σ | T (periods) | last mean − w_σ^L | last mean − w_σ | first mean − w_σ^L | the last periods' drift | oscillation | lag | bias |",
               "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in C['rows']:
            L_.append(f"| {r['lab']} | {r['b0']:g} | {r['sig']:g} | {r['T']:.0f} ({r['n_per']}) | {fmt(r['d_L'])} | {fmt(r['d_id'])} | "
                      f"{fmt(r['d_L_first'])} | {fmt(r['drift'])} | {fmt(r['osc'])} | {fmt(r['lag'])} | {fmt(r['bias'])} |")
        L_.append("")
    if D is not None:
        L_ += [f"## D. On an orbit's Gram ({D12} atoms of width {WID12:g} on [−1.4, 1.4]; eigenvalues {fmt(min(D['eig']))} to {fmt(max(D['eig']))})", "",
               f"(1) A 3-sparse truth on atoms {D['truth']} (signs +, −, +; its irrepresentable quantity {D['ic_truth']:.3f}): the equilibrium "
               "against the ridge (lazy start, β⁰ = 30) and the LASSO of penalty μ = σλ₁(β⁰) (rich start, β⁰ = 10⁻⁶); its support "
               "(|w| > 10⁻³ max|w|), whether it keeps the truth's atoms, the atoms beyond them, max |w − w*| (absolute distances):", "",
               "| start | μ (rich) or σ (lazy) | distance to the reference | support | the truth kept | extra atoms | max |w − w*| |",
               "|---|---|---|---|---|---|---|"]
        for r in D['rows']:
            L_.append(f"| {r['lab']} | {r['mu']:g} | {fmt(r['d'])} | {r['sup']} | {'yes' if r['truth_kept'] else 'no'} | {r['extra']} | {fmt(r['err'])} |")
        ic = D['ic']
        L_ += ["", f"(2) Every 3-sparse support and sign pattern ({ic['n']}, magnitudes 0.15–0.3), the rich start at μ = {MU_IC:g}: the "
               f"irrepresentable quantity ‖H_{{S^cS}}H_SS⁻¹s‖∞ is below 1 for {ic['ic_lt1']}; the LASSO's prediction (below 1 and no sign "
               f"flip on S) holds for {ic['pred']}; the equilibrium's support is exactly the truth's for {ic['exact']} — the prediction and "
               f"the equilibrium agree on {ic['agree']} of {ic['n']}; where the support is not the truth's it carries {ic['extra_when_fail']:.2f} "
               f"extra atoms on average (largest solver residual {fmt(ic['max_res'])}).", "",
               "(3) A biased disturbance correlated with the regressors (d̄ = Hv, v of size 0.05, a lazy start): the equilibrium's distance "
               "to w* and to w* + H⁻¹d̄ as σ → 0 — "
               + '; '.join(f"σ = {r['sig']:g}: {fmt(r['d_star'])}, {fmt(r['d_shift'])}" for r in D['bias']) + ".", "",
               "(4) The averaged loop from the rich start (β⁰ = 10⁻⁶, μ = 10⁻³): the time at which max|w − w_σ| < 10⁻³ max|w_σ| — "
               + '; '.join(f"{r['lab']} (irrepresentable quantity {r['ic']:.3f}, the equilibrium's support {r['sup']}): "
                           + (f"t ≈ {r['reach']:.2g}" if r['reach'] else "not by t = 10⁸") for r in D['dyn']) + ".", ""]
    if Es is not None:
        L_ += report_E(Es)
    L_.append("Numbers are produced by `experiments/exp_mirror_leak.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mirror_leak_raw.npz')
    md = os.path.join(OUT, 'mirror_leak.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDE']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDE')
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('D', part_D), ('B', part_B), ('C', part_C), ('E', part_E)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDE']))


if __name__ == '__main__':
    main()
