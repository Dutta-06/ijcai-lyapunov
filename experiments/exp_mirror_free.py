"""
Proposition D‴.45: the free mass–amplitude law is a mirror flow — a conservation law per particle, the implicit bias as a Bregman
projection, the noise's invariance in continuous time, the sampled-data drift, the real loop's Lyapunov function
====================================================================================================================================

The replicator–amplitude law of Proposition D‴.34 with free masses (λ ≡ 0), each particle driven by one scalar signal g_i:
dα_i = Gβ_ia′_ig_i dt, dβ_i = γβ_ia_ig_i dt; Q_i = β_i − c cosh(2α_i/a_max) (c = γa_max²/(4G)) is conserved for every drive, and on a
dictionary w = β⊙a/m is the mirror flow dη = Φᵀ(r + ξ)dt/N of the potential Ψ_q, ψ′ = η = ∫dα/(Gβa′), ψ* = (β − β⁰)/(γm)
(`mirror_free`).

A. Identities: η's closed form against quadrature (the three cases and the branch); the conjugate and the metric; Q along the loop,
   noiseless and driven by a smooth node-correlated coloured noise; the conserved law's shift dQ_i = −γβ_iλ dt, integrated (ΔQ_i along
   the conserved law's trajectory against the quadrature of −γβ_iλ, λ = Σ_jw_jg_j from the state); Q's drift along Proposition D‴.44's
   reduced flow (free masses, white noise, the paradigm's update) against the Euler step's closed form.
B. The implicit bias: from α = 0 and equal masses s·m/d, the Bregman projection (its residual recorded; a projection whose residual
   exceeds 10⁻⁹(1 + |y|) counts as a failure and is left out of every aggregate) against the minimum-norm fit (the lazy regime) and basis
   pursuit (the rich regime), for s = 10 … 10⁻⁴ (N × d = 5 × 12, 2- and 3-sparse truths), medians and largest values; the noiseless
   loop's end for s ≥ 10⁻², and once at s = 10⁻³ (the time to reach the projection); unequal masses (weighted basis pursuit, weighted
   minimum norm); and from small random amplitudes at tiny masses (the branches Q < −c, Proposition D‴.47): whether the frozen signs
   still admit an exact fit, the coordinates vanishing on the sign-constrained fits, the truth's atoms on a wrong-signed branch, the
   distance to basis pursuit — as distributions.
C. Noise: (1) white label noise at a fixed intensity as the step h → 0: Euler–Maruyama on (α, log β) (Itô), Heun (Stratonovich) and the
   mirror update — Q's change against the Itô correction's time integral (Itô keeps a drift as h → 0, Stratonovich's vanishes like h, the
   mirror update keeps Q exactly); (2) the Euler step's exact mean change of Q (Gauss–Hermite over each particle's Gaussian drive) against
   the closed form at h = 0.5, 0.05, 0.005; (3) along the paradigm's reduced flow (D‴.44's selection, free masses, white noise) from the
   lazy start, on instances where basis pursuit is the truth and where it is not, from s = 10² to 10⁸: ‖w‖₁, the distance to basis
   pursuit, the end support against basis pursuit's, Q; (4) the sampled-data loop with label noise: the paradigm's Euler update moves Q and
   the projection, the mirror update on the curves moves neither.
D. The real loop (a scalar plant, Gaussian atoms at fixed centres, two particles on one atom without cohesion between them,
   x_m = sin t + 0.4 sin 2.7t): V = ½e² + Σ Bregman divergences monotone with V̇ = −k_ee²; the pair's mirror coordinates move together
   (η_i − η_i⁰ = η_j − η_j⁰) under any common drive — noiseless, with a smooth noise in the adaptation only, and with a sensor noise in
   the feedback, the regressor and the adaptation.
E. The conserved law (Σβ = m) from the same start: Q moves, and its end differs from the free law's projection.
F. The rich limit's leading order (unequal masses εr_i, the instances of B): where weighted basis pursuit is a vertex with N atoms, its
   certificate's slacks μ_j, the projection against the leading order computed from the certificate (`rich_leading_order`) and against
   the closed form a_max(εr_j)^{3/2}cos(πμ_j/2)/(m(2c)^{1/2}sin³(πμ_j/2)) off the support (`core_tail`); the distance's ratio between
   ε = 10⁻⁷ and 10⁻⁹ (ε^{3/2}: 1000); the instances whose distance at ε = 10⁻⁷ exceeds 10⁻³ followed to 10⁻¹¹ — small slacks
   (near-optimal vertices) delay the asymptotic regime until ε^{1/2} ≪ μ_min.

Output: results/chaos/mirror_free.md and mirror_free_raw.npz (MF_WORKERS; `--from-raw` regenerates the md; `--only A,B,C,D,E,F`).
"""
import os
import sys
import time
import warnings
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import mirror_free as MF
from src.core import mass_dictionary as MD
from src.core import colour_selection as CS
from src.core.runcache import cached
import exp_cohesion_leak as E

OUT = os.environ.get('MF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('MF_WORKERS', 2))
N_I, D_I, M_I, G_I, GAM_I, AMAX = 5, 12, 12.0, 1.0, 0.1, 1.0
SCALES_B, KS_B, N_INST_B, T_B = (10.0, 1.0, 0.1, 1e-2, 1e-3, 1e-4), (2, 3), 20, 1e7
EPS_RICH, SCALES_LAZY = (1e-3, 1e-5, 1e-7), (1e2, 1e4)                     # unequal masses: ε·r and scale·r
EPS_F, EPS_SLOW, MU_SLOW = (1e-5, 1e-7, 1e-9), (1e-11,), 0.01                # Part F: the leading order; the slow instances
FLOW_MIN_B = 1e-2                                                           # the flow is integrated for s ≥ this (slower below)
FLOW_S3, T_S3 = 1e-3, 3e7                                                   # one instance's flow at s = 10⁻³ (k = 2, instance 0)
OK_TOL = 1e-9                                                               # a projection fails when |Φw − y| > OK_TOL(1 + |y|)
S_BR = (1e-2, 1e-4, 1e-6)                                                   # the branches' start scales
S_C3 = (0.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
I_C3 = (0, 1, 2, 3, 4, 5, 8, 9, 13, 17)                                     # k = 2; 8, 9, 13, 17: basis pursuit ≠ the truth
H_C1, SIG_C1, T_C1, P_C1 = (0.1, 0.03, 0.01), 0.3, 200.0, 64
KE_D, GAM_D, M_D, T_D = 5.0, 0.5, 6.0, 2000.0
H_C4, SIG_C4, N_C4, P_C4 = 0.5, 0.3, 20000, 16


def law():
    return MF.FreeLaw(M_I, G_I, GAM_I, AMAX)


def instance(k, i):
    """A 5 × 12 Gaussian design with unit columns and a k-sparse truth (signs ±1, magnitudes in [0.15, 0.3])."""
    rng = np.random.default_rng(7000 + 100 * k + i)
    Phi = rng.standard_normal((N_I, D_I)); Phi /= np.linalg.norm(Phi, axis=0)
    ws = np.zeros(D_I); idx = rng.choice(D_I, k, replace=False); ws[idx] = rng.choice([-1, 1], k) * rng.uniform(0.15, 0.3, k)
    return Phi, ws


def smooth_noise(rng, N, K=30, amp=0.3, lo=1e-3, hi=3.0):
    """A smooth noise on the nodes: K sinusoids of log-uniform frequencies in [lo, hi], mixed over the nodes by a random matrix."""
    om = np.exp(rng.uniform(np.log(lo), np.log(hi), K)); ph = rng.uniform(0, 2 * np.pi, K); Mx = amp * rng.standard_normal((N, K))
    return lambda t: Mx @ np.sin(om * t + ph)


def flow_end(L, Phi, y, z0, T, conserved=False, xi=None):
    sol = solve_ivp(L.rhs(Phi, y, xi, conserved), (0.0, T), z0, method='LSODA', rtol=1e-10, atol=1e-14)
    return sol.y[:, -1]


def part_A():
    warnings.filterwarnings('ignore')
    L = law()
    f = lambda x, q: 1.0 / (L.G * L.beta(x, q) * L.amp(x)[1])
    qs = (100 * L.c, 3.0, 0.5, L.c, 0.5 * L.c, 0.0, -0.5 * L.c)
    eta_err = 0.0
    for q in qs:
        for al in np.linspace(-4, 4, 17):
            ref = quad(lambda x: f(x, q), 0.0, al, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
            eta_err = max(eta_err, abs(L.eta(al, q)[0] - ref) / max(1.0, abs(ref)))
    qb = -2 * L.c; a0 = float(L.branch_point(qb))                              # the branch, from a reference point
    for al in a0 + np.array([0.2, 0.5, 1.0, 2.0, 3.5]):
        ref = quad(lambda x: f(x, qb), a0 + 0.1, al, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
        eta_err = max(eta_err, abs((L.eta(al, qb) - L.eta(a0 + 0.1, qb))[0] - ref) / max(1.0, abs(ref)))
    conj_err = metric_err = 0.0
    for q in (0.97, 0.01, -0.02):
        for al in (-1.3, 0.2, 2.2):
            hh = 1e-5
            dw = (L.w(al + hh, q) - L.w(al - hh, q)) / (2 * hh)
            de = (L.eta(al + hh, q) - L.eta(al - hh, q))[0] / (2 * hh)
            dps = (L.psi_star(al + hh, q) - L.psi_star(al - hh, q)) / (2 * hh)
            conj_err = max(conj_err, abs(dps / de - L.w(al, q)) / max(abs(L.w(al, q)), 1e-3))
            metric_err = max(metric_err, abs(dw / de - L.m * L.relevance(al, q)) / (L.m * L.relevance(al, q)))
    cons = []
    for i in range(3):
        Phi, ws = instance(2, i); y = Phi @ ws; d = D_I
        rng = np.random.default_rng(300 + i)
        z0 = np.r_[0.05 * rng.standard_normal(d), np.log(np.full(d, 0.3))]
        q0 = L.Q(z0[:d], np.exp(z0[d:]))
        out = {}
        for lab, xi, T in (('noiseless', None, 1e4), ('smooth noise', smooth_noise(rng, N_I), 3e4)):
            sol = solve_ivp(L.rhs(Phi, y, xi), (0.0, T), z0, method='DOP853', rtol=1e-11, atol=1e-13, dense_output=True)
            ts = np.linspace(0, T, 60)
            out[lab] = float(max(np.abs(L.Q(sol.sol(t)[:d], np.exp(sol.sol(t)[d:])) - q0).max() for t in ts))
        # the conserved law, integrated: ΔQ_i along its trajectory against the quadrature of −γβ_iλ, λ = Σ_jβ_ja_jg_j/Σ_jβ_j from the state
        xi9 = smooth_noise(np.random.default_rng(9), N_I)
        Tc = 200.0
        solc = solve_ivp(L.rhs(Phi, y, xi9, conserved=True), (0.0, Tc), z0, method='DOP853', rtol=1e-12, atol=1e-14, dense_output=True)
        tg = np.linspace(0.0, Tc, 20001)
        Zc = solc.sol(tg)
        al_c, b_c = Zc[:d], np.exp(Zc[d:])
        a_c = L.amp(al_c)[0]
        r_c = y[:, None] - Phi @ (b_c * a_c / L.m) + np.array([xi9(t) for t in tg]).T
        g_c = Phi.T @ r_c / N_I
        lam_c = (b_c * a_c * g_c).sum(axis=0) / b_c.sum(axis=0)
        integrand = -L.gamma * b_c * lam_c[None, :]
        quadQ = np.sum(0.5 * (integrand[:, 1:] + integrand[:, :-1]) * np.diff(tg)[None, :], axis=1)
        dQc = L.Q(al_c[:, -1], b_c[:, -1]) - L.Q(al_c[:, 0], b_c[:, 0])
        out['shift'] = float(np.abs(dQc - quadQ).max() / np.abs(dQc).max())
        cons.append(out)
    # Q's drift along D‴.44's reduced flow (free masses, white noise, the paradigm's update) against the Euler step's closed form
    red = []
    for i in range(3):
        Phi, ws = instance(2, i); y = Phi @ ws
        X = MD.MassDictionary(Phi, y, M_I, a_max=AMAX, G=G_I, gamma=GAM_I, balanced=False)
        al0 = np.zeros(D_I); q = L.Q(al0, np.full(D_I, M_I / D_I))
        _, al, _, _ = L.bregman_projection(Phi, y, q, L.eta(al0, q))
        th = X.theta(al, L.beta(al, q))
        dr = X.drift(th, lambda k: np.ones_like(k), ito=True)
        b = L.beta(al, q); a, da = L.amp(al)
        dQ = -L.c * (2 / L.A) * np.sinh(2 * al / L.A) * dr[:D_I] + 2 * np.sqrt(b) * dr[D_I:]
        pred = (2 * (Phi ** 2).sum(axis=0) / (L.m * N_I)) * L.gamma * b * (L.gamma * a ** 2 - L.G * b * da * (1 + (a / L.A) ** 2))
        red.append(float(np.abs(dQ - pred).max() / np.abs(pred).max()))
    return dict(eta_err=eta_err, conj_err=conj_err, metric_err=metric_err, cons=cons, red=red)


@cached
def b_job(args):
    """Part B, one (k, instance): the projection (its residual recorded) for every start scale, the flow's end for s ≥ FLOW_MIN_B (and
    the time to reach the projection at s = FLOW_S3 for k = 2, instance 0); unequal masses; the branches (Proposition D‴.47)."""
    warnings.filterwarnings('ignore')
    k, i = args
    t0 = time.time()
    L = law()
    Phi, ws = instance(k, i); y = Phi @ ws
    tol_y = OK_TOL * (1 + np.linalg.norm(y))
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    wl2 = np.linalg.pinv(Phi) @ y
    out = dict(k=k, i=i, bp_ok=bool(np.abs(wbp - ws).max() < 1e-4 * np.abs(ws).max()), l1_bp=float(np.abs(wbp).sum()))
    for s in SCALES_B:
        al0 = np.zeros(D_I); b0 = np.full(D_I, s * M_I / D_I)
        q = L.Q(al0, b0)
        w, al, it, res = L.bregman_projection(Phi, y, q, L.eta(al0, q))
        r = dict(conv=bool(res <= tol_y), res=res, it=it, l1=float(np.abs(w).sum()), used=int((np.abs(w) > 1e-3).sum()),
                 d_bp=float(np.abs(w - wbp).max()), d_l2=float(np.abs(w - wl2).max() / np.abs(wl2).max()),
                 ok=bool(np.abs(w - ws).max() < 1e-3 * np.abs(ws).max()), amax=float(np.abs(L.amp(al)[0]).max() / AMAX),
                 bmax=float(L.beta(al, q).max()))
        if s >= FLOW_MIN_B:
            z = flow_end(L, Phi, y, np.r_[al0, np.log(b0)], T_B)
            wT = np.exp(z[D_I:]) * L.amp(z[:D_I])[0] / L.m
            r.update(flow_res=float(np.linalg.norm(Phi @ wT - y)), flow_d=float(np.abs(wT - w).max()))
        out[s] = r
    if (k, i) == (2, 0):                                                         # the loop at s = 10⁻³: the time to reach the projection
        al0 = np.zeros(D_I); b0 = np.full(D_I, FLOW_S3 * M_I / D_I); q = L.Q(al0, b0)
        w, al, it, res = L.bregman_projection(Phi, y, q, L.eta(al0, q))
        ts = np.geomspace(1e2, T_S3, 41)
        sol = solve_ivp(L.rhs(Phi, y), (0.0, T_S3), np.r_[al0, np.log(b0)], method='LSODA', rtol=1e-10, atol=1e-14, t_eval=ts)
        dist = [float(np.abs(np.exp(sol.y[D_I:, j]) * L.amp(sol.y[:D_I, j])[0] / L.m - w).max()) for j in range(sol.y.shape[1])]
        out['flow3'] = dict(t=sol.t.tolist(), d=dist, res=res, status=int(sol.status))
    # the branches (Proposition D‴.47): small random amplitudes at tiny masses — the frozen signs, the sign-constrained fits S, the
    # projection onto S
    for s in S_BR:
        al0 = 0.05 * np.random.default_rng(900 + 100 * k + i).standard_normal(D_I); b0 = np.full(D_I, s * M_I / D_I)
        P = MF.ParticleSet(L, al0, b0)
        yh, Z, _ = MF.sign_ls_set(Phi, y, P.s)
        keep = np.setdiff1d(np.arange(D_I), Z)
        w, it, res = MF.projection(P, Phi, yh, P.eta(P.x0), keep)
        truth = np.flatnonzero(ws)
        out[('branch', s)] = dict(n_branch=int(P.br.sum()), yres=float(np.linalg.norm(y - yh)), exact=bool(np.linalg.norm(y - yh) <= tol_y),
                                  nZ=int(Z.size), wrong=int(sum(bool(P.br[j]) and P.s[j] != np.sign(ws[j]) for j in truth)),
                                  conv=bool(res <= OK_TOL * (1 + np.linalg.norm(yh))), res=res, it=it,
                                  used=int((np.abs(w) > 1e-3).sum()), d_bp=float(np.abs(w - wbp).max()))
    # unequal masses (ratios r in [1/4, 4]): the rich limit against weighted basis pursuit (weights 1/√β⁰), the lazy one against the
    # weighted minimum norm (weights 1/β⁰²)
    r = np.exp(np.random.default_rng(800 + 100 * k + i).uniform(np.log(0.25), np.log(4.0), D_I))
    wl1 = MF.weighted_bp(Phi, y, 1.0 / np.sqrt(r))
    for eps in EPS_RICH:
        q = L.Q(np.zeros(D_I), eps * r)
        w, _, it, res = L.bregman_projection(Phi, y, q, L.eta(np.zeros(D_I), q))
        out[('rich', eps)] = dict(d=float(np.abs(w - wl1).max()), d_unw=float(np.abs(w - wbp).max()), res=res, conv=bool(res <= tol_y))
    for big in SCALES_LAZY:
        b0 = big * r; q = L.Q(np.zeros(D_I), b0)
        w, _, it, res = L.bregman_projection(Phi, y, q, L.eta(np.zeros(D_I), q))
        wl2 = b0 ** 2 * (Phi.T @ np.linalg.solve(Phi @ ((b0 ** 2)[:, None] * Phi.T), y))
        out[('lazy', big)] = dict(d=float(np.abs(w - wl2).max() / np.abs(wl2).max()), res=res, conv=bool(res <= tol_y))
    print(f"[{time.time() - t0:.0f}s] B k {k} instance {i}", flush=True)
    return out


def f_job(args):
    """Part F, one (k, instance): unequal masses εr (Part B's r), the projection from α = 0 (exact fallback) against weighted basis pursuit,
    the leading order from its certificate and the closed form off the support, where the vertex has N atoms."""
    warnings.filterwarnings('ignore')
    k, i = args
    L = law()
    Phi, ws = instance(k, i); y = Phi @ ws
    r = np.exp(np.random.default_rng(800 + 100 * k + i).uniform(np.log(0.25), np.log(4.0), D_I))
    rho = 1.0 / np.sqrt(r)
    wl1 = MF.weighted_bp(Phi, y, rho)
    S = np.flatnonzero(np.abs(wl1) > 1e-10)
    out = dict(k=k, i=i, nS=int(S.size), rows=[])
    if S.size != N_I:
        return out
    mu, xi = MF.certificate_slacks(Phi, wl1, S, rho)
    off = np.setdiff1d(np.arange(D_I), S)
    out['mu_min'] = float(mu[off].min())
    for eps in EPS_F + EPS_SLOW:
        if eps in EPS_SLOW and out['rows'][1]['d'] <= 1e-3:
            break
        P = MF.ParticleSet(L, np.zeros(D_I), eps * r)
        w, it, res = MF.projection(P, Phi, y, P.eta(P.x0), exact=True)
        w1 = MF.rich_leading_order(P, Phi, y, wl1, S)
        cf = np.sign(Phi.T @ xi) * L.core_tail(eps * r, mu)
        out['rows'].append(dict(eps=eps, it=it, res=res, d=float(np.abs(w - wl1).max()), d1=float(np.abs(w1 - wl1).max()),
                                agree=float(np.abs(w - w1).max()), cf=float(np.max(np.abs(w[off] - cf[off]) / np.abs(cf[off]))),
                                supp=[int(j) for j in np.flatnonzero(np.abs(w) > 1e-6)]))
    return out


def part_C12():
    """C1: white label noise at a fixed intensity as h → 0 — Euler–Maruyama on (α, log β) (Itô), Heun (Stratonovich) and the mirror
    update: Q's change against the time integral of the Itô correction (γβ_i/2)(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²))σ²|φ_i|²/N² along each
    path. C2: the Euler step's exact mean change of Q (Gauss–Hermite over each particle's Gaussian drive) against the closed form."""
    warnings.filterwarnings('ignore')
    L = law()
    c1 = []
    for i in range(3):
        Phi, ws = instance(2, i); y = Phi @ ws; d = D_I
        rng0 = np.random.default_rng(400 + i)
        al_s = 0.3 * rng0.standard_normal(d); b_s = rng0.uniform(0.2, 2.0, d)
        q = L.Q(al_s, b_s)
        phi2 = (Phi ** 2).sum(axis=0)
        F = lambda al, b: (L.G * b * L.amp(al)[1], L.gamma * L.amp(al)[0])
        gdet = lambda al, b: (y[None, :] - (b * L.amp(al)[0] / L.m) @ Phi.T) @ Phi / N_I
        rate = lambda al, b: 0.5 * L.gamma * b * (L.gamma * L.amp(al)[0] ** 2 - L.G * b * L.amp(al)[1] * (1 + (L.amp(al)[0] / L.A) ** 2)) \
            * SIG_C1 ** 2 * phi2 / N_I ** 2
        for h in H_C1:
            n = int(round(T_C1 / h))
            for scheme in ('Euler (Itô)', 'Heun (Stratonovich)', 'mirror'):
                rng = np.random.default_rng(4100 + i)                           # the same noise paths for every scheme
                al = np.tile(al_s, (P_C1, 1)); b = np.tile(b_s, (P_C1, 1)); pred = np.zeros((P_C1, d))
                for _ in range(n):
                    pred += rate(al, b) * h
                    dxi = (SIG_C1 * np.sqrt(h) / N_I) * (rng.standard_normal((P_C1, N_I)) @ Phi)
                    de1 = gdet(al, b) * h + dxi
                    if scheme == 'Euler (Itô)':
                        fa, fb = F(al, b)
                        al, b = al + fa * de1, b * np.exp(fb * de1)
                    elif scheme == 'Heun (Stratonovich)':
                        fa, fb = F(al, b)
                        al1, b1 = al + fa * de1, b * np.exp(fb * de1)
                        de2 = gdet(al1, b1) * h + dxi
                        fa1, fb1 = F(al1, b1)
                        al, b = al + 0.5 * (fa * de1 + fa1 * de2), b * np.exp(0.5 * (fb * de1 + fb1 * de2))
                    else:
                        al = L.mirror_step(al, q, de1 / h, h); b = L.beta(al, q)
                dQ = L.Q(al, b) - q[None, :]
                c1.append(dict(i=i, h=h, scheme=scheme, dQ=float(dQ.mean()), pred=float(pred.mean()),
                               ratio=float(dQ.mean() / pred.mean()), dQmax=float(np.abs(dQ).max())))
    c2 = []
    xg, wg = np.polynomial.hermite.hermgauss(80)
    for h in (0.5, 0.05, 0.005):
        for i in range(3):
            Phi, ws = instance(2, i); y = Phi @ ws; d = D_I
            rng = np.random.default_rng(500 + i)
            al = 0.4 * rng.standard_normal(d); b = rng.uniform(0.2, 2.0, d)
            sig = 0.3
            a, da = L.amp(al)
            r0 = y - Phi @ (b * a / L.m)
            mu = Phi.T @ r0 / N_I; sd = sig * np.sqrt((Phi ** 2).sum(axis=0)) / N_I
            g = mu[:, None] + np.sqrt(2) * sd[:, None] * xg[None, :]
            EdQ = (L.dQ_exact(al[:, None], b[:, None], g, h) * wg[None, :]).sum(axis=1) / np.sqrt(np.pi)
            pred = L.dQ_euler(al, b, mu ** 2 + sd ** 2, h)
            c2.append(dict(h=h, rel=float(np.abs(EdQ - pred).max() / np.abs(pred).max())))
    return dict(c1=c1, c2=c2)


def part_C4():
    """C4: the sampled-data loop at the step H_C4 with label noise drawn afresh (σ = SIG_C4, P_C4 paths, N_C4 steps) from the lazy
    projection (α = 0, β = m/d): the paradigm's Euler update on (α, log β) against the mirror update η⁺ = η + hg on each particle's curve."""
    warnings.filterwarnings('ignore')
    L = law()
    rows = []
    for i in range(2):
        Phi, ws = instance(2, i); y = Phi @ ws
        wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
        al00 = np.zeros(D_I); q0 = L.Q(al00, np.full(D_I, M_I / D_I))
        w0, al0, _, _ = L.bregman_projection(Phi, y, q0, L.eta(al00, q0))
        rng = np.random.default_rng(600 + i)
        out = dict(i=i, d_bp0=float(np.abs(w0 - wbp).max()), l1_0=float(np.abs(w0).sum()))
        for upd in ('Euler', 'mirror'):
            al = np.tile(al0, (P_C4, 1)); b = np.tile(L.beta(al0, q0), (P_C4, 1))
            for n in range(N_C4):
                a, da = L.amp(al)
                r = y[None, :] - (b * a / L.m) @ Phi.T + SIG_C4 * rng.standard_normal((P_C4, N_I))
                g = r @ Phi / N_I
                if upd == 'Euler':
                    al, b = al + H_C4 * L.G * b * da * g, b * np.exp(H_C4 * L.gamma * a * g)
                else:
                    al = L.mirror_step(al, q0, g, H_C4); b = L.beta(al, q0)
            Qp = L.Q(al, b)
            proj = []
            for p_ in range(P_C4):
                wp = L.bregman_projection(Phi, y, Qp[p_], L.eta(al[p_], Qp[p_]), sgn=np.sign(al[p_]))[0]
                proj.append(wp)
            proj = np.array(proj)
            out[upd] = dict(dQ_mean=float((Qp - q0).mean()), dQ_max=float(np.abs(Qp - q0).max()),
                            d_proj=float(np.abs(proj - w0).max(axis=1).mean()), d_bp=float(np.abs(proj - wbp).max(axis=1).mean()),
                            l1=float(np.abs(proj).sum(axis=1).mean()), res=float(np.linalg.norm((b * L.amp(al)[0] / L.m) @ Phi.T - y, axis=1).mean()),
                            proj_mean=proj.mean(axis=0).tolist())
        # Proposition D‴.44's reduced flow (the step's filter) to the same s = (m v/4) t, v = hσ²/N
        X = MD.MassDictionary(Phi, y, M_I, a_max=AMAX, G=G_I, gamma=GAM_I, balanced=False)
        s_end = M_I * (H_C4 * SIG_C4 ** 2 / N_I) / 4 * N_C4 * H_C4
        filt = lambda k: MD.LN.lr_factor(k, H_C4, ('fresh',))
        th = X.flow(X.theta(al0, L.beta(al0, q0)), s_end, filt, ito=True, chunks=20)[1][-1]
        alr, _, br, _, _, _ = X.parts(th); wr = X.w(th)
        out['reduced'] = dict(s=float(s_end), l1=float(np.abs(wr).sum()), dQ_mean=float((L.Q(alr, br) - q0).mean()),
                              d_euler=float(np.abs(np.array(out['Euler']['proj_mean']) - wr).max()))
        rows.append(out)
    return rows


@cached
def c3_job(i):
    """Part C3, one instance: D‴.44's reduced flow (free masses, white noise, the paradigm's update) from the lazy projection (s = 1)."""
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L = law()
    Phi, ws = instance(2, i); y = Phi @ ws
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    X = MD.MassDictionary(Phi, y, M_I, a_max=AMAX, G=G_I, gamma=GAM_I, balanced=False)
    al0 = np.zeros(D_I); q0 = L.Q(al0, np.full(D_I, M_I / D_I))
    _, al, _, _ = L.bregman_projection(Phi, y, q0, L.eta(al0, q0))
    th0 = X.theta(al, L.beta(al, q0))
    _, Th = X.flow(th0, S_C3[-1], lambda k: np.ones_like(k), ito=True, t_eval=np.array(S_C3), chunks=60)
    rows = []
    for th in Th:
        alx, _, b, a, _, _ = X.parts(th); w = X.w(th); Qx = L.Q(alx, b)
        rows.append(dict(l1=float(np.abs(w).sum()), d_bp=float(np.abs(w - wbp).max()), used=int((np.abs(w) > 1e-3).sum()),
                         qmean=float(Qx.mean()), qmin=float(Qx.min()), qmax=float(Qx.max()), dead=int((Qx <= -L.c).sum()),
                         supp=[int(j) for j in np.flatnonzero(np.abs(w) > 1e-3)]))
    supp_bp = [int(j) for j in np.flatnonzero(np.abs(wbp) > 1e-3)]
    print(f"[{time.time() - t0:.0f}s] C3 instance {i}", flush=True)
    return dict(i=i, rows=rows, bp_ok=bool(np.abs(wbp - ws).max() < 1e-4 * np.abs(ws).max()), supp_bp=supp_bp,
                l1_bp=float(np.abs(wbp).sum()), end_is_bp=bool(rows[-1]['supp'] == supp_bp))


def part_D():
    """The real loop: Gaussian atoms at fixed centres (six, width 0.4) plus a second particle on the third atom (no cohesion between
    them, σ₀ = 0); x_m = sin t + 0.4 sin 2.7t."""
    warnings.filterwarnings('ignore')
    L = MF.FreeLaw(M_D, G_I, GAM_D, AMAX)
    rng = np.random.default_rng(11)
    cen6 = np.linspace(-1.2, 1.2, 6)
    cen = np.r_[cen6, cen6[2]]; wid = np.full(7, 0.4)
    al0 = 0.1 * rng.standard_normal(7); b0 = rng.uniform(0.5, 1.5, 7)
    q = L.Q(al0, b0)
    wtrue6 = rng.uniform(-0.15, 0.15, 6)
    # the pair's split: the Bregman split of wtrue6[2] given the start (η_i − η_i⁰ = η_j − η_j⁰ and w_i + w_j = wtrue6[2])
    from scipy.optimize import brentq
    e0 = L.eta(al0, q)
    pair_w = lambda dl: L.w(L.alpha_of_eta(np.array([e0[2] + dl]), q[2:3])[0], q[2]) + \
        L.w(L.alpha_of_eta(np.array([e0[6] + dl]), q[6:7])[0], q[6])
    dl = brentq(lambda x: pair_w(x) - wtrue6[2], -200, 200, xtol=1e-14)
    wi = L.w(L.alpha_of_eta(np.array([e0[2] + dl]), q[2:3])[0], q[2]); wj = wtrue6[2] - wi
    wstar = np.r_[wtrue6[:2], wi, wtrue6[3:], wj]                     # the target's coefficients, the pair at its Bregman split
    al_star = L.alpha_of_w(wstar, q)
    xr, xrd = (lambda t: np.sin(t) + 0.4 * np.sin(2.7 * t)), (lambda t: np.cos(t) + 1.08 * np.cos(2.7 * t))
    out = {}
    nz = smooth_noise(np.random.default_rng(12), 1, K=20, amp=0.1, lo=0.05, hi=20.0)
    for lab, noise in (('noiseless', None), ('a smooth noise in the adaptation only', nz), ('a sensor noise (feedback, regressor, adaptation)', nz)):
        F0 = L.real_rhs(cen, wid, wstar, xr, xrd, KE_D)
        if noise is None:
            F = F0
        elif lab.startswith('a smooth'):
            def F(t, z, F0=F0, noise=noise):                                 # the adaptation sees e + n(t); the plant does not
                x, al, lb = z[0], z[1:8], z[8:]
                b = np.exp(lb); a, da = L.amp(al)
                phi = np.exp(-(x - cen) ** 2 / (2 * wid ** 2))
                em = x - xr(t) + noise(t)[0]
                return np.r_[F0(t, z)[0], L.G * b * da * em * phi, L.gamma * a * em * phi]
        else:
            def F(t, z, noise=noise):                                         # the measured state x + n in the feedback, f̂ and the adaptation
                x, al, lb = z[0], z[1:8], z[8:]
                b = np.exp(lb); a, da = L.amp(al)
                xm_ = x + noise(t)[0]
                phi = np.exp(-(xm_ - cen) ** 2 / (2 * wid ** 2))
                em = xm_ - xr(t)
                f_true = np.exp(-(x - cen) ** 2 / (2 * wid ** 2)) @ wstar
                return np.r_[f_true + xrd(t) - KE_D * em - phi @ (b * a / L.m), L.G * b * da * em * phi, L.gamma * a * em * phi]
        z0 = np.r_[xr(0.0) + 0.3, al0, np.log(b0)]
        ts = np.linspace(0, T_D, 8001)
        sol = solve_ivp(F, (0, T_D), z0, method='DOP853', rtol=1e-12, atol=1e-14, t_eval=ts)
        Z = sol.y
        V = np.array([L.real_lyapunov(Z[:, j], t, xr, al_star, q) for j, t in enumerate(ts)])
        e = Z[0] - xr(ts)
        dVerr = np.nan
        if noise is None:
            dV = np.empty_like(ts)
            for j, t in enumerate(ts):
                zd = F0(t, Z[:, j]); hh = 1e-6
                dV[j] = (L.real_lyapunov(Z[:, j] + hh * zd, t + hh, xr, al_star, q) - L.real_lyapunov(Z[:, j] - hh * zd, t - hh, xr, al_star, q)) / (2 * hh)
            dVerr = float(np.abs(dV + KE_D * e ** 2).max() / np.abs(KE_D * e ** 2).max())
        etas = np.array([L.eta(Z[1:8, j], q) for j in range(len(ts))])
        pair = np.abs((etas[:, 2] - e0[2]) - (etas[:, 6] - e0[6])).max() / max(1.0, np.abs(etas[:, 2] - e0[2]).max())
        Qs = np.array([L.Q(Z[1:8, j], np.exp(Z[8:, j])) for j in range(len(ts))])
        wend = np.exp(Z[8:, -1]) * L.amp(Z[1:8, -1])[0] / L.m
        out[lab] = dict(V0=float(V[0]), Vend=float(V[-1]), monotone=bool(np.all(np.diff(V) <= 1e-12 * max(1.0, V[0]))), dVerr=dVerr,
                        pair=float(pair), dQ=float(np.abs(Qs - q).max()), e_end=float(np.abs(e[-400:]).max()),
                        d6=float(np.abs(np.r_[wend[:2], wend[3:6]] - np.r_[wtrue6[:2], wtrue6[3:]]).max()),
                        dsum=float(abs(wend[2] + wend[6] - wtrue6[2])))
    return dict(out=out, wi=float(wi), wj=float(wj), wsum=float(wtrue6[2]))


def part_E():
    """The conserved law (Σβ = m) against the free law from α = 0, β = m/d."""
    warnings.filterwarnings('ignore')
    L = law()
    rows = []
    for i in range(6):
        Phi, ws = instance(2, i); y = Phi @ ws
        al0 = np.zeros(D_I); b0 = np.full(D_I, M_I / D_I); z0 = np.r_[al0, np.log(b0)]
        q0 = L.Q(al0, b0)
        wf = L.bregman_projection(Phi, y, q0, L.eta(al0, q0))[0]
        z = flow_end(L, Phi, y, z0, 1e6, conserved=True)
        wc = np.exp(z[D_I:]) * L.amp(z[:D_I])[0] / L.m
        rows.append(dict(res=float(np.linalg.norm(Phi @ wc - y)), d=float(np.abs(wc - wf).max()), dQ=float(np.abs(L.Q(z[:D_I], np.exp(z[D_I:])) - q0).max()),
                         l1c=float(np.abs(wc).sum()), l1f=float(np.abs(wf).sum()), mass=float(np.exp(z[D_I:]).sum())))
    return rows


def report(A=None, B=None, C=None, D=None, Ee=None, F=None):
    fmt = E.fmt
    L = law()
    L_ = ["# Proposition D‴.45: the free mass–amplitude law is a mirror flow (`exp_mirror_free.py`)", "",
          f"Free masses (λ ≡ 0): dα_i = Gβ_ia′_ig_i, dβ_i = γβ_ia_ig_i, a = a_max tanh(α/a_max), w = β⊙a/m; m = {M_I:g}, G = {G_I:g}, "
          f"γ = {GAM_I:g}, a_max = {AMAX:g} (c = γa_max²/(4G) = {L.c:g}) unless said otherwise; Q_i = β_i − c cosh(2α_i/a_max); on a "
          "dictionary the mirror flow dη = Φᵀ(r + ξ)dt/N of Ψ_q, ψ′ = η = ∫dα/(Gβa′) (`mirror_free`).", ""]
    if A is not None:
        L_ += ["## A. Identities", "",
               f"η's closed form against quadrature (q = 100c, 3, 0.5, c, c/2, 0, −c/2 and the branch q = −2c, 17 and 5 points): largest "
               f"relative difference {fmt(A['eta_err'])}; (ψ*)′ = w against finite differences {fmt(A['conj_err'])}; dw/dη = mD "
               f"{fmt(A['metric_err'])}.", "",
               "Q along the loop (5 × 12, 2-sparse truths, β(0) = 0.3, α(0) ~ N(0, 0.05²)), largest |Q(t) − Q(0)| over 60 times; the "
               "conserved law's shift: ΔQ_i along its integrated trajectory (t ≤ 200, the smooth noise) against the quadrature of −γβ_iλ, "
               "λ from the state (relative):", "",
               "| instance | noiseless, t ≤ 10⁴ | a smooth noise (30 sinusoids, frequencies 10⁻³–3, mixed over the nodes, amplitude 0.3), t ≤ 3·10⁴ "
               "| the conserved law's shift |", "|---|---|---|---|"]
        for j, r in enumerate(A['cons']):
            L_.append(f"| {j} | {fmt(r['noiseless'])} | {fmt(r['smooth noise'])} | {fmt(r['shift'])} |")
        L_ += ["", "Q's drift along Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's update; at the lazy projection) "
               "against the Euler step's closed form (2|φ_i|²/(mN))γβ_i(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²)) per unit s: largest relative "
               "difference " + ', '.join(fmt(x) for x in A['red']) + ".", ""]
    if B is not None:
        nfail = lambda xs: sum(1 for x in xs if not x['conv'])
        med = lambda v: fmt(float(np.median(v))) if len(v) else '—'
        mx = lambda v: fmt(float(np.max(v))) if len(v) else '—'
        L_ += [f"## B. The implicit bias: the start's masses set the regime ({N_I} × {D_I}, {N_INST_B} instances per k)", "",
               "From α = 0 and equal masses β = s·m/d: the Bregman projection (Newton on the dual; a projection whose residual exceeds "
               f"{OK_TOL:g}(1 + |y|) is a failure, counted and left out of every aggregate) — ‖w‖₁ over basis pursuit's, the atoms used "
               "(|w| > 10⁻³), the fraction recovering the truth (10⁻³·max|w*|), the median and largest distances to basis pursuit, the "
               f"median relative distance to the minimum-norm fit Φ⁺y, the largest |a|/a_max and β; for s ≥ {FLOW_MIN_B:g} the noiseless "
               f"loop's end at t = {T_B:.0e} against the projection (its largest residual and distance).", "",
               "| k | s | failures | ‖w‖₁/‖w_bp‖₁ | used | the truth | ‖w − w_bp‖∞: median, largest | median ‖w − Φ⁺y‖∞/‖Φ⁺y‖∞ | max |a|/a_max, max β | "
               "the loop's end: residual, distance to the projection |", "|---|---|---|---|---|---|---|---|---|---|"]
        for k in KS_B:
            rs = [r for r in B if r['k'] == k]
            L_.append(f"| {k} | basis pursuit | | 1 | | {np.mean([r['bp_ok'] for r in rs]):.2f} | | | | |")
            for s_ in SCALES_B:
                allq = [r[s_] for r in rs]; q = [x for x in allq if x['conv']]; rq = [r for r in rs if r[s_]['conv']]
                fl = (f"{fmt(max(x['flow_res'] for x in allq))}, {fmt(max(x['flow_d'] for x in allq))}" if 'flow_d' in allq[0] else "—")
                L_.append(f"| {k} | {s_:g} | {nfail(allq)} | {np.mean([x['l1'] / r['l1_bp'] for x, r in zip(q, rq)]):.4f} | "
                          f"{np.mean([x['used'] for x in q]):.1f} | {np.mean([x['ok'] for x in q]):.2f} | {med([x['d_bp'] for x in q])}, "
                          f"{mx([x['d_bp'] for x in q])} | {med([x['d_l2'] for x in q])} | {np.max([x['amax'] for x in q]):.3f}, "
                          f"{np.max([x['bmax'] for x in q]):.3g} | {fl} |")
        f3 = [r['flow3'] for r in B if 'flow3' in r]
        if f3:
            f3 = f3[0]
            reach = next((t for t, d_ in zip(f3['t'], f3['d']) if d_ < 1e-8), None)
            L_ += ["", f"The loop at s = {FLOW_S3:g} (k = 2, instance 0): distance to its projection "
                   + ", ".join(f"{fmt(d_)} at t = {t:.0e}" for t, d_ in list(zip(f3['t'], f3['d']))[::8])
                   + (f"; below 10⁻⁸ from t ≈ {reach:.1e}." if reach else "; not below 10⁻⁸ by the end.")]
        L_ += ["", "Small random amplitudes (α(0) ~ N(0, 0.05²)) at tiny masses: the particles with Q < −c sit on branches (Proposition "
               "D‴.47) — the sign-constrained fits S, the coordinates vanishing on S (|Z|), the instances with a truth atom on a "
               "wrong-signed branch, whether the frozen signs still fit y exactly, and the distribution of the projection's distance to "
               "basis pursuit:", "",
               "| k | s | branches (mean) | failures | exact fit kept | a truth atom on a wrong-signed branch | |Z| (mean) | ‖w − w_bp‖∞: "
               "≤ 10⁻³, median, largest |", "|---|---|---|---|---|---|---|---|"]
        for k in KS_B:
            rs = [r for r in B if r['k'] == k]
            for s_ in S_BR:
                allq = [r[('branch', s_)] for r in rs]; q = [x for x in allq if x['conv']]
                db = [x['d_bp'] for x in q]
                L_.append(f"| {k} | {s_:g} | {np.mean([x['n_branch'] for x in allq]):.1f} | {nfail(allq)} | "
                          f"{sum(x['exact'] for x in allq)}/{len(allq)} | {sum(x['wrong'] > 0 for x in allq)}/{len(allq)} | "
                          f"{np.mean([x['nZ'] for x in allq]):.1f} | {sum(d_ <= 1e-3 for d_ in db)}/{len(db)}, {med(db)}, {mx(db)} |")
        L_ += ["", "Unequal masses β⁰_i = s·r_i (r_i log-uniform in [1/4, 4]), α = 0: the projection against weighted basis pursuit "
               "(min Σ|w_i|/√r_i) as s → 0 — and against the unweighted one — and against the weighted minimum norm (min Σw_i²/r_i²) as "
               "s → ∞ (relative); median and largest over the converged projections (failures in brackets):", "",
               "| k | " + " | ".join(f"s = {e:g}: weighted (median, largest), unweighted (median)" for e in EPS_RICH)
               + " | " + " | ".join(f"s = {b_:g}: weighted ℓ² (median, largest)" for b_ in SCALES_LAZY) + " |",
               "|---|" + "---|" * (len(EPS_RICH) + len(SCALES_LAZY))]
        for k in KS_B:
            rs = [r for r in B if r['k'] == k]
            cells = []
            for e in EPS_RICH:
                q = [r[('rich', e)] for r in rs if r[('rich', e)]['conv']]
                cells.append(f"{med([x['d'] for x in q])}, {mx([x['d'] for x in q])}, {med([x['d_unw'] for x in q])} [{len(rs) - len(q)}]")
            for b_ in SCALES_LAZY:
                q = [r[('lazy', b_)] for r in rs if r[('lazy', b_)]['conv']]
                cells.append(f"{med([x['d'] for x in q])}, {mx([x['d'] for x in q])} [{len(rs) - len(q)}]")
            L_.append(f"| {k} | " + " | ".join(cells) + " |")
        L_.append("")
    if C is not None:
        L_ += ["## C. Noise", "", f"(1) White label noise at the intensity σ = {SIG_C1:g} (dη_i = (Φᵀr)_idt/N + σ(ΦᵀdW)_i/N), {P_C1} paths "
               f"to t = {T_C1:g} from random states (α ~ N(0, 0.3²), β ~ U(0.2, 2)): Q's mean change against the mean time integral of the "
               "Itô correction (γβ_i/2)(γa_i² − Gβ_ia′_i(1 + a_i²/a_max²))σ²|φ_i|²/N² along the same paths, for three schemes as the step h → 0:", "",
               "| instance | h | scheme | mean ΔQ | the Itô integral | ratio | max |ΔQ| |", "|---|---|---|---|---|---|---|"]
        for r in C['c1']:
            L_.append(f"| {r['i']} | {r['h']:g} | {r['scheme']} | {r['dQ']:+.4e} | {r['pred']:+.4e} | {r['ratio']:.4f} | {fmt(r['dQmax'])} |")
        L_ += ["", "(2) The Euler step on (α, log β) with label noise (σ = 0.3, at random states): its exact mean change of Q (Gauss–Hermite, "
               "80 nodes, over each particle's Gaussian drive) against the closed form (h²E[g²]γβ/2)(γa² − Gβa′(1 + a²/a_max²)) — largest "
               "relative difference over the particles, three instances per h: "
               + '; '.join(f"h = {h_:g}: " + ', '.join(fmt(r['rel']) for r in C['c2'] if r['h'] == h_) for h_ in sorted({r['h'] for r in C['c2']}, reverse=True))
               + ".", ""]
        if 'c3' in C:
            L_ += ["(3) The selection in the mirror coordinates: Proposition D‴.44's reduced flow (free masses, white noise, the paradigm's "
                   "update) from the lazy projection (s = 1) — ‖w‖₁, the distance to basis pursuit, the atoms used, the dead atoms (Q ≤ −c), "
                   "Q's mean and range; at the end the support against basis pursuit's:", "",
                   "| instance (basis pursuit = truth) | s | ‖w‖₁ | ‖w − w_bp‖∞ | used | dead | Q: mean, min, max |", "|---|---|---|---|---|---|---|"]
            for r in C['c3']:
                for s_, x in zip(S_C3, r['rows']):
                    L_.append(f"| {r['i']} ({'yes' if r['bp_ok'] else 'no'}) | {s_:.0e} | {x['l1']:.4f} | {fmt(x['d_bp'])} | {x['used']} | "
                              f"{x.get('dead', '—')} | {x['qmean']:.4g}, {x['qmin']:.4g}, {x['qmax']:.4g} |")
            L_ += ["", "The end supports: " + "; ".join(
                f"instance {r['i']}: {r['rows'][-1].get('supp')} against basis pursuit's {r.get('supp_bp')} (‖w‖₁ {r['rows'][-1]['l1']:.4f} "
                f"against {r.get('l1_bp', float('nan')):.4f}) — {'equal' if r.get('end_is_bp') else 'different'}" for r in C['c3']) + "."]
            L_.append("")
        if 'c4' in C:
            L_ += [f"(4) The sampled-data loop at the step h = {H_C4:g} with label noise drawn afresh (σ = {SIG_C4:g}, {P_C4} paths, "
                   f"{N_C4} steps) from the lazy projection (α = 0, β = m/d): the paradigm's Euler update on (α, log β) against the mirror "
                   "update η⁺ = η + hg on each particle's curve — Q's change (mean over the particles and paths, largest), the projection of "
                   "the end state: its distance to the start's (mean over the paths), to basis pursuit's, its ‖w‖₁; the end's residual:", "",
                   "| instance | the start: ‖w‖₁, ‖w − w_bp‖∞ | the update | Q − Q(0): mean, max | ‖P − P(0)‖∞ | ‖P − w_bp‖∞ | ‖P‖₁ | |Φw − y| |",
                   "|---|---|---|---|---|---|---|---|"]
            for r in C['c4']:
                for upd in ('Euler', 'mirror'):
                    x = r[upd]
                    L_.append(f"| {r['i']} | {r['l1_0']:.4f}, {fmt(r['d_bp0'])} | {upd} | {x['dQ_mean']:+.3e}, {fmt(x['dQ_max'])} | "
                              f"{fmt(x['d_proj'])} | {fmt(x['d_bp'])} | {x['l1']:.4f} | {x['res']:.3f} |")
            L_ += ["", "Proposition D‴.44's reduced flow (the step's filter) to the same s = (m v/4)t: its ‖w‖₁ and mean change of Q, and the "
                   "Euler paths' mean projection against its state: " + '; '.join(
                       f"instance {r['i']} (s = {r['reduced']['s']:g}): {r['reduced']['l1']:.4f}, {r['reduced']['dQ_mean']:+.3e}, {fmt(r['reduced']['d_euler'])}"
                       for r in C['c4'] if 'reduced' in r) + "."]
            L_.append("")
    if D is not None:
        L_ += [f"## D. The real loop (k_e = {KE_D:g}, γ = {GAM_D:g}, m = {M_D:g}, six Gaussian atoms of width 0.4 on [−1.2, 1.2] and a second "
               f"particle on the third, without cohesion between them; x_m = sin t + 0.4 sin 2.7t; {T_D:g} time units)", "",
               f"The target's coefficient on the shared atom w*₃ = {D['wsum']:.5f}, at the pair's Bregman split {D['wi']:.5f} + {D['wj']:.5f} "
               "(V's reference).", "",
               "| the noise | V(0) | V(T) | V monotone | max |V̇ + k_ee²|/max k_ee² | the pair's mirror coordinates "
               "(max |Δη_i − Δη_j|, relative) | max |Q − Q(0)| | max |e| (last 10%) | the distinct atoms' |w − w*| | the pair's sum's error |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for lab, r in D['out'].items():
            L_.append(f"| {lab} | {fmt(r['V0'])} | {fmt(r['Vend'])} | {'yes' if r['monotone'] else 'no'} | {fmt(r['dVerr'])} | {fmt(r['pair'])} | "
                      f"{fmt(r['dQ'])} | {fmt(r['e_end'])} | {fmt(r['d6'])} | {fmt(r['dsum'])} |")
        L_.append("")
    if Ee is not None:
        L_ += ["## E. The conserved law (Σβ = m) from the same start (α = 0, β = m/d; t = 10⁶)", "",
               "| instance | its residual | max |Q − Q(0)| | ‖w‖₁: conserved, free (the projection) | the ends' distance | Σβ |", "|---|---|---|---|---|---|"]
        for j, r in enumerate(Ee):
            L_.append(f"| {j} | {fmt(r['res'])} | {fmt(r['dQ'])} | {r['l1c']:.4f}, {r['l1f']:.4f} | {fmt(r['d'])} | {r['mass']:.6f} |")
        L_.append("")
    if F is not None:
        L_ += ["## F. The rich limit's leading order (unequal masses εr_i, the instances of B; Proposition D‴.45(iii))", "",
               "Where weighted basis pursuit (weights r^{−1/2}) is a vertex with N atoms: its certificate's smallest relative slack μ_min off "
               "the support; the projection from α = 0 (exact fallback) — its distance |w − w_bp|∞ — against the leading order computed from "
               "the certificate (|w − w₁|∞, absolute: the rounding floor is ≈ 10⁻¹²) and against the closed form a_max(εr_j)^{3/2}cos(πμ_j/2)/"
               f"(m(2c)^{{1/2}}sin³(πμ_j/2)) off the support (its largest relative error); medians and largest values over the instances with "
               f"μ_min ≥ {MU_SLOW:g}:", ""]
        for k in KS_B:
            q = [f for f in F if f['k'] == k]
            v = [f for f in q if f['nS'] == N_I]
            g = [f for f in v if f['mu_min'] >= MU_SLOW]
            L_ += [f"k = {k}: {len(v)} of {len(q)} instances with an N-atom vertex, {len(g)} with μ_min ≥ {MU_SLOW:g}; μ_min median "
                   f"{fmt(np.median([f['mu_min'] for f in v]))}, smallest {fmt(min(f['mu_min'] for f in v))}.", "",
                   "| ε | distance: median, largest | |w − w₁|∞: median, largest | the closed form: median, largest | projection: largest residual |",
                   "|---|---|---|---|---|"]
            for j, eps in enumerate(EPS_F):
                rows = [f['rows'][j] for f in g]
                L_.append(f"| {eps:g} | {fmt(np.median([rr['d'] for rr in rows]))}, {fmt(max(rr['d'] for rr in rows))} | "
                          f"{fmt(np.median([rr['agree'] for rr in rows]))}, {fmt(max(rr['agree'] for rr in rows))} | "
                          f"{fmt(np.median([rr['cf'] for rr in rows]))}, {fmt(max(rr['cf'] for rr in rows))} | {fmt(max(rr['res'] for rr in rows))} |")
            ratio = [f['rows'][1]['d'] / f['rows'][2]['d'] for f in g]
            L_ += ["", f"The distance's ratio between ε = 10⁻⁷ and 10⁻⁹ (ε^{{3/2}}: 1000): median {np.median(ratio):.0f}, range "
                   f"{min(ratio):.0f}–{max(ratio):.0f}.", ""]
        slow = [f for f in F if f['nS'] == N_I and f['mu_min'] < MU_SLOW]
        if slow:
            L_ += [f"The instances with μ_min < {MU_SLOW:g} (near-optimal second vertices; those whose distance at ε = 10⁻⁷ exceeds 10⁻³ "
                   "followed to 10⁻¹¹):", "",
                   "| instance | μ_min | ε | distance | the leading order's distance | |w − w₁|∞ | support |", "|---|---|---|---|---|---|---|"]
            for f in slow:
                for rr in f['rows']:
                    L_.append(f"| k = {f['k']}, {f['i']} | {fmt(f['mu_min'])} | {rr['eps']:g} | {fmt(rr['d'])} | {fmt(rr['d1'])} | "
                              f"{fmt(rr['agree'])} | {rr['supp']} |")
            L_.append("")
    L_.append("Numbers are produced by `experiments/exp_mirror_free.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'mirror_free_raw.npz')
    md = os.path.join(OUT, 'mirror_free.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEF']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEF')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(k, i) for k in KS_B for i in range(N_INST_B)], chunksize=1)
        if 'C' in parts:
            h['C3'] = pool.map_async(c3_job, list(I_C3), chunksize=1)
        if 'F' in parts:
            h['F'] = pool.map_async(f_job, [(k, i) for k in KS_B for i in range(N_INST_B)], chunksize=1)
        if 'A' in parts:
            res['A'] = part_A(); E.save(raw, **res); log('A')
        if 'D' in parts:
            res['D'] = part_D(); E.save(raw, **res); log('D')
        if 'E' in parts:
            res['E'] = part_E(); E.save(raw, **res); log('E')
        if 'C' in parts:
            C = part_C12(); C['c4'] = part_C4(); C['c3'] = h['C3'].get(); res['C'] = C; E.save(raw, **res); log('C')
        if 'B' in parts:
            res['B'] = h['B'].get(); E.save(raw, **res); log('B')
        if 'F' in parts:
            res['F'] = h['F'].get(); E.save(raw, **res); log('F')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEF']))


if __name__ == '__main__':
    main()
