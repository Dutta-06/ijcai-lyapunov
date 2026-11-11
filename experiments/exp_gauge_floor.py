"""
Proposition D‴.28: the adaptive prior's floor removal is a gauge property — it survives an approximation error
==========================================================================================================

The three-bump instance with a fourth bump that three kernels cannot represent, f_ε = f + ε·exp(−(x − 0.8)²/(2·0.3²)), ε ∈ EPS; K₀ = 3, the
masses q = 1/3, σ₀ = 0.15σ_crit. ϑ^b(q): the best weighted K₀-atom approximant (a nondegenerate critical point of the unregularised weighted
loss L_q, continued in ε from the exact representation of f_0). The amplitude–mass gauge G_w (α_k ↦ a⁻¹(q_ka(α_k)/w_k)) leaves f̂ unchanged,
so L_w∘G_w = L_q.

A. No dynamics, per ε: the approximant's fit error; the whitened unregularised Hessian Ã = H_L − R̃ (its smallest eigenvalue, ‖R̃‖/‖H_L‖);
   the gauge set W_g = {w_k > q_k|a^b_k|/a_max} and its probability P(w^m ∈ W_g) at m = 15, 60, 240; the gauge check at sampled weights of
   W_g — the unregularised weighted gradient at G_wϑ^b(q), and the critical point reached by Newton continued in w from ϑ^b(q) (not started
   at the gauge transform): its distance to G_wϑ^b(q) and its f̂ against f̂_b(q); plain Newton from ϑ^b(q) at w, which may land elsewhere;
   the function-level sensitivity of the adaptive prior's stationary point (critical points reached from ϑ^b(q), central differences in w)
   and of the fixed matching prior's (ϑ⁰ = ϑ^b(q)) against σ₀Ĵ(Ã + σ₀I)⁻¹g̃ and against D‴.24's form with H_L in place of Ã; the fixed
   prior's floor law tr(FC̃Fᵀ)/N; the residual curvature at ϑ^b(q) — its amplitude row (zero at stationarity) and λ⁺_k := λ_max(R^G_k) —
   and the adaptive prior's splitting margins σ₀ − (q_k/w_k)λ⁺_k against the margins computed at G_wϑ^b(q), over all count vectors of
   m = 15, 60 with the probability of a negative margin, and against the separation eigenvalues of the dynamics (the Jacobian of the
   averaged adaptive drift at the replicated configuration, counts (6, 4, 5)); the branch ϑ^b(q) continued to ε = 0.68, near the
   amplitude's saturation; at the Part B samples outside W_g, the minimum of L_w over the closed amplitude box from the clipped gauge
   transform (the systems' limit) and the existence of a gauge transform with the classes permuted.
B. Dynamics at ε = 0.1 and 0.3: the nearest-atom leak toward the fixed matching prior ϑ^b(q) (η = 0) and the adaptive prior (η = σ₀, ∞) from
   warm clouds around ϑ^b(q) (width 0.02), m = 15, 60, 240 against M = 510, four seeds, S_B slow units: ‖f̂_m − f̂_M‖² and the fit error at
   s = 600, 3000, S_B, the counts, the class means' distance from G_wϑ^b(q), and the end diameters of the classes.

Output: results/chaos/gauge_floor.md and gauge_floor_raw.npz (GF_WORKERS, GF_SB, GF_OUT; `--only-B` computes Part B, saving after
every job and skipping the jobs already saved; `--from-raw` regenerates the md).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import ot_leak as OT
from src.core import prior_design as PD
from src.core import adaptive_prior as AP
from src.core import floor_bias as FB
from src.core.ot_leak import _cost_matrix

OUT = os.environ.get('GF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
EXTRA = (0.8, 0.3)                                                         # the fourth bump's centre and width
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q = np.ones(3) / 3
EPS = (0.0, 0.05, 0.1, 0.2, 0.3, 0.5)
EPS_B = (0.1, 0.3)
EPS_END = (0.6, 0.65, 0.68)
MS_WG = (15, 60, 240)
MS, M_REF = (15, 60, 240), 510
JITTER = 0.02
DT, PLAN_EVERY, REC = 0.25, 4, 200
S_B = float(os.environ.get('GF_SB', 6000.0))
N_SEEDS = 4
N_WORKERS = int(os.environ.get('GF_WORKERS', 2))
H_FD = 1e-5


def target(eps):
    def f(X):
        X = np.asarray(X, float)
        return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS) + eps * np.exp(-(X - EXTRA[0]) ** 2 / (2 * EXTRA[1] ** 2))
    return f


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def net_():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact0(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def best_approximant(eps, n_steps=None):
    """ϑ^b(q) at ε: the critical point of the unregularised weighted loss at q, continued in ε from f_0's exact representation."""
    net = net_()
    th = exact0(net)
    n_steps = n_steps or max(1, int(np.ceil(eps / 0.025)))
    for k in range(1, n_steps + 1):
        e = eps * k / n_steps
        loop = OT.OTPriorLoop(net, target(e), x_ref, PER, G, SIGMA0, th, n_quad=N_Q)
        th, res, fit = FB.exact_weighted_atoms(loop, Q, th, n_cont=1, w0=Q, fit_tol=1.0)
    return net, th, res, fit


def setup(eps):
    """The loop with the fixed matching prior ϑ⁰ = ϑ^b(q) and the approximant."""
    net, tb, res, fit = best_approximant(eps)
    loop = OT.OTPriorLoop(net, target(eps), x_ref, PER, G, SIGMA0, tb, n_quad=N_Q)
    return net, tb, loop, res, fit


def gauge(net, tb, w):
    """G_wϑ^b(q): the amplitudes q_ka^b_k/w_k (None outside W_g)."""
    a, c, s = net.physical(tb)
    an = Q[:, None] * a / np.asarray(w, float)[:, None]
    if np.any(np.abs(an) >= net.a_max):
        return None
    return net.unconstrained(an, c, s)


def fhat(loop, T, w):
    a, c, s, d, phi = loop.features(T)
    return phi @ (np.asarray(w, float) * a)


def hessians(loop, tb, h=1e-6):
    """At ϑ^b(q): the whitened weighted unregularised Hessian Ã, the Gauss–Newton block H_L, Ĵ and the whitened gauge displacement g̃ along
    e₁ − e₃, e₂ − e₃ (g̃ := −D^{1/2}∂_w(G_wϑ^b), so that f̂'s direct weight effect is Ĵg̃)."""
    K, p = tb.shape; d = K * p
    gradL = lambda v: (Q[:, None] * FB.unregularised_grad(loop, Q, v.reshape(K, p))).ravel()          # ∇_{ϑ_k}L_q = q_k∇V_L(ϑ_k)
    A = np.zeros((d, d))
    for j in range(d):
        E = np.zeros(d); E[j] = h
        A[:, j] = (gradL(tb.ravel() + E) - gradL(tb.ravel() - E)) / (2 * h)
    A = 0.5 * (A + A.T)
    Dh = np.kron(np.sqrt(Q), np.sqrt(1.0 / G))                              # D^{1/2}, D = diag(q_kG⁻¹)
    At = A / Dh[:, None] / Dh[None, :]                                      # whitened: D^{−1/2}AD^{−1/2}
    Jh = FB.whitened_jacobian(loop, Q, tb)
    HL = Jh.T @ Jh / loop.N
    net = loop.net
    Vd = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
    g = np.column_stack([((gauge(net, tb, Q + h * Vd[:, c]) - gauge(net, tb, Q - h * Vd[:, c])) / (2 * h)).ravel() for c in range(2)])
    gt = -Dh[:, None] * g                                                   # g̃ := −D^{1/2}∂_w(G_wϑ^b): f̂'s direct weight effect is Ĵg̃
    return dict(A=A, At=At, HL=HL, Jh=Jh, gt=gt, Vd=Vd)


def jac(fun, v, h=1e-7):
    """The central-difference Jacobian of fun at v."""
    d = v.size; J = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h * max(1.0, abs(v[j]))
        J[:, j] = (fun(v + e) - fun(v - e)) / (2 * e[j])
    return J


def sep_basis(lab, K=3, p=3, extra=0):
    """An orthonormal basis of the separations (class means fixed, the prior's offsets zero) in ℝ^{mp + extra}: Helmert contrasts per class."""
    m = len(lab); cols = []
    for k in range(K):
        idx = np.where(lab == k)[0]
        for a in range(1, len(idx)):
            c = np.zeros(m); c[idx[:a]] = 1.0; c[idx[a]] = -a; c /= np.linalg.norm(c)
            for j in range(p):
                v = np.zeros(m * p + extra); v[np.arange(m) * p + j] = c; cols.append(v)
    return np.array(cols).T


def prob_Wg(lo, m):
    """P(w^m ∈ W_g) for the multinomial sample w^m = n/m of the masses q: the count vectors with n_k/m > lo_k for every k."""
    from scipy.stats import multinomial
    tot = 0.0
    for n1 in range(m + 1):
        for n2 in range(m + 1 - n1):
            n = np.array([n1, n2, m - n1 - n2])
            if np.all(n / m > lo):
                tot += multinomial.pmf(n, m, Q)
    return float(tot)


def part_A_eps(eps):
    t0 = time.time()
    net, tb, loop, res, fit = setup(eps)
    Hs = hessians(loop, tb)
    At, HL, Jh, gt, Vd = Hs['At'], Hs['HL'], Hs['Jh'], Hs['gt'], Hs['Vd']
    Rt = HL - At
    ev = np.linalg.eigvalsh(At)
    a_b = net.physical(tb)[0][:, 0]
    lo = Q * np.abs(a_b) / net.a_max                                        # W_g: w_k > lo_k
    # the gauge check at sampled weights
    rng = np.random.default_rng(1)
    gchk = []
    for _ in range(6):
        while True:
            w = rng.dirichlet(30 * np.ones(3))
            if np.all(w > lo + 0.01):
                break
        tg = gauge(net, tb, w)
        gnorm = float(np.abs(FB.unregularised_grad(loop, w, tg)).max())
        tc = None
        for nc in (20, 80, 320):                                                              # continued in w from ϑ^b(q), refined on a fold
            try:
                tc, rc, fc = FB.exact_weighted_atoms(loop, w, tb, n_cont=nc, w0=Q, fit_tol=1.0)
                break
            except np.linalg.LinAlgError:
                continue
        if tc is None:
            gchk.append(dict(w=w, grad=gnorm, dist=float('inf'), dfh=float('inf'), res=float('nan'), dist_plain=float('nan'), res_plain=float('nan')))
            continue
        try:
            tp, rp, fp = FB.exact_weighted_atoms(loop, w, tb, n_cont=1, w0=w, fit_tol=1.0)   # plain Newton from ϑ^b(q) at w
            dpl = float(np.abs(tp - tg).max())
        except np.linalg.LinAlgError:                                                         # a singular Jacobian on the way: no limit
            dpl, rp = float('inf'), float('nan')
        gchk.append(dict(w=w, grad=gnorm, dist=float(np.abs(tc - tg).max()), dfh=float(np.abs(fhat(loop, tc, w) - fhat(loop, tb, Q)).max()), res=rc,
                         dist_plain=dpl, res_plain=rp))
    # sensitivities: adaptive (critical points of L_w reached by Newton from ϑ^b(q)) and fixed (Θ*(w) toward ϑ⁰ = ϑ^b(q)), central differences
    Fa, Ff = [], []
    for c in range(2):
        o = {}
        for sg in (+1, -1):
            w = Q + sg * H_FD * Vd[:, c]
            ta = FB.exact_weighted_atoms(loop, w, tb, n_cont=1, w0=w, fit_tol=1.0)[0]
            tf, rf = PD.weighted_stationary_atoms_prior(loop.nearest, w, tb)
            o[sg] = (fhat(loop, ta, w), fhat(loop, tf, w))
        Fa.append((o[1][0] - o[-1][0]) / (2 * H_FD)); Ff.append((o[1][1] - o[-1][1]) / (2 * H_FD))
    Fa, Ff = np.array(Fa).T, np.array(Ff).T
    F_gauge = SIGMA0 * Jh @ np.linalg.solve(At + SIGMA0 * np.eye(9), gt)    # σ₀Ĵ(Ã + σ₀I)⁻¹g̃
    F_gn = SIGMA0 * Jh @ np.linalg.solve(HL + SIGMA0 * np.eye(9), gt)       # D‴.24's form with H_L
    Ct = (np.diag(Q) - np.outer(Q, Q))[:2, :2]
    law = lambda F: float(np.trace(F @ Ct @ F.T)) / loop.N
    # the residual curvature at ϑ^b(q): its amplitude row and λ⁺_k; the adaptive prior's margins
    a, c_, s_, d_, phi = loop.features(tb)
    r = loop.F - phi @ (Q * a)
    Rk = loop.residual_curvature(r, tb)                                     # (K, 3, 3), per particle
    sqG = np.sqrt(G)
    RG = sqG[None, :, None] * Rk * sqG[None, None, :]
    amp_row = float(np.abs(Rk[:, 0, :]).max())
    lam = np.array([np.linalg.eigvalsh(0.5 * (RG[k] + RG[k].T))[-1] for k in range(3)])
    # (iii) at the level of the dynamics: the separation eigenvalues of the averaged nearest-atom adaptive drift at the replicated
    # configuration G_wϑ^b(q) (counts (6, 4, 5), m = 15) against −σ₀ + (q_k/w_k)eig(G^{1/2}R_k(q)G^{1/2}), each n_k − 1 times
    sepchk = {}
    counts = np.array([6, 4, 5]); m = 15; w = counts / m; lab = np.repeat(np.arange(3), counts)
    tg = gauge(net, tb, w)
    if tg is not None:
        ev_q = [np.linalg.eigvalsh(0.5 * (RG[k] + RG[k].T)) for k in range(3)]
        pred = np.sort(np.concatenate([np.repeat(-SIGMA0 + (Q[k] / w[k]) * ev_q[k], counts[k] - 1) for k in range(3)]))
        for eta in (SIGMA0, np.inf):
            apl = AP.AdaptivePriorLoop(loop, eta, plan='nearest')
            Th = tg[lab].copy(); Z = tg.copy()
            gamma, _ = apl.plan(Th, Z)
            if np.isinf(eta):
                fun = lambda v: apl.drift(v.reshape(m, 3), Z, gamma)[0].ravel()
                v0 = Th.ravel(); extra = 0
            else:
                def fun(v):
                    d1, d2 = apl.drift(v[:3 * m].reshape(m, 3), v[3 * m:].reshape(3, 3), gamma)
                    return np.concatenate([d1.ravel(), d2.ravel()])
                v0 = np.concatenate([Th.ravel(), Z.ravel()]); extra = 9
            J = jac(fun, v0)
            S = sep_basis(lab, extra=extra)
            JS = J @ S
            es = np.sort(np.linalg.eigvals(S.T @ JS).real)
            sepchk[eta] = dict(dmax=float(np.abs(es - pred).max()), largest=float(es[-1]), pred_largest=float(pred[-1]),
                               inv=float(np.abs(JS - S @ (S.T @ JS)).max() / np.abs(JS).max()), drift0=float(np.abs(fun(v0)).max()))
    pW = {m_: prob_Wg(lo, m_) for m_ in MS_WG}
    print(f"[{time.time() - t0:.0f}s] A eps {eps:g}", flush=True)
    return dict(eps=eps, fit=fit, res=res, ev=ev, R_rel=float(np.linalg.norm(Rt) / np.linalg.norm(HL)), lo=lo, a_b=a_b, gauge=gchk,
                Fa=float(np.linalg.norm(Fa) / np.linalg.norm(Ff)), law_fixed=law(Ff), law_gauge=law(F_gauge), law_gn=law(F_gn),
                rel_gauge=float(np.linalg.norm(F_gauge - Ff) / np.linalg.norm(Ff)), rel_gn=float(np.linalg.norm(F_gn - Ff) / np.linalg.norm(Ff)),
                amp_row=amp_row, lam=lam, R_norm=float(np.abs(Rk).max()), sepchk=sepchk, pW=pW)


def branch_end():
    """ϑ^b(q) continued in ε beyond Part A's range, toward the amplitude's saturation: the fit error, max|a^b_k|/a_max, the smallest
    eigenvalue of Ã, max_kλ⁺_k/σ₀ and W_g's thresholds q_k|a^b_k|/a_max."""
    rows = []
    for eps in EPS_END:
        net, tb, res, fit = best_approximant(eps)
        loop = OT.OTPriorLoop(net, target(eps), x_ref, PER, G, SIGMA0, tb, n_quad=N_Q)
        Hs = hessians(loop, tb)
        a_b = net.physical(tb)[0][:, 0]
        a, c_, s_, d_, phi = loop.features(tb)
        Rk = loop.residual_curvature(loop.F - phi @ (Q * a), tb)
        sqG = np.sqrt(G)
        lam = max(np.linalg.eigvalsh(0.5 * (sqG[:, None] * Rk[k] * sqG[None, :] + (sqG[:, None] * Rk[k] * sqG[None, :]).T))[-1] for k in range(3))
        rows.append(dict(eps=eps, fit=fit, res=res, amax=float(np.abs(a_b).max() / net.a_max), ev=float(np.linalg.eigvalsh(Hs['At'])[0]),
                         lam=float(lam), lo=Q * np.abs(a_b) / net.a_max))
    print("[A end]", flush=True)
    return rows


def outside_job(args):
    """At a weight vector w outside W_g: the minimum of L_w over the closed box (|a| ≤ a_max, |c| ≤ c_max, s ∈ [s_min, s_max]; L-BFGS-B in the
    physical coordinates) from the clipped gauge transform — the limit of systems that keep ϑ^b(q)'s branch — its relative fit error and
    the amplitudes at the bound; and the permutations of the classes under which a gauge transform exists (then the fit error is the
    approximant's)."""
    from scipy.optimize import minimize
    from itertools import permutations
    eps, w = args
    w = np.asarray(w, float)
    net, tb, loop, res, fit = setup(eps)
    a_b, c_b, s_b = (x.copy() for x in net.physical(tb))
    X, F = loop.X, loop.F

    def fun(v):
        a, c, s = v[:3], v[3:6], v[6:9]
        dx = X[:, None] - c[None, :]
        E = np.exp(-dx ** 2 / (2 * s[None, :] ** 2))
        r = F - E @ (w * a)
        g = np.concatenate([-(r @ (E * w[None, :])), -(r @ (E * (w * a)[None, :] * dx / s[None, :] ** 2)),
                            -(r @ (E * (w * a)[None, :] * dx ** 2 / s[None, :] ** 3))]) / len(X)
        return 0.5 * float(np.mean(r ** 2)), g
    a0 = np.clip(Q * a_b[:, 0] / w, -net.a_max * (1 - 1e-9), net.a_max * (1 - 1e-9))
    v0 = np.concatenate([a0, c_b[:, 0], s_b])
    bnds = [(-net.a_max, net.a_max)] * 3 + [(-net.c_max, net.c_max)] * 3 + [(net.s_min, net.s_max)] * 3
    r_ = minimize(fun, v0, jac=True, method='L-BFGS-B', bounds=bnds, options=dict(maxiter=50000, ftol=1e-16, gtol=1e-13))
    fit_box = float(np.sqrt(2 * r_.fun / np.mean(F ** 2)))
    lo = Q * np.abs(a_b[:, 0]) / net.a_max
    perms = [pp for pp in permutations(range(3)) if np.all(w[list(pp)] > lo)]
    return dict(eps=eps, w=w, fit_box=fit_box, fit_b=fit, at_bound=int(np.sum(np.abs(r_.x[:3]) > net.a_max * (1 - 1e-6))), perms=perms)


def margins_job(args):
    """The adaptive prior's splitting margins over all count vectors of m: predicted σ₀ − (q_k/w_k)λ⁺_k against computed at G_wϑ^b(q)."""
    eps, m, lam = args
    from scipy.stats import multinomial
    net, tb, loop, res, fit = setup(eps)
    rows = []
    for n1 in range(1, m - 1):
        for n2 in range(1, m - n1):
            n = np.array([n1, n2, m - n1 - n2]); w = n / m
            tg = gauge(net, tb, w)
            if tg is None:
                continue
            a, c, s, d, phi = loop.features(tg)
            sep = SIGMA0 - loop.kappa_R(loop.F - phi @ (w * a), tg)
            pred = SIGMA0 - (Q / w) * lam
            rows.append(dict(n=n, p=float(multinomial.pmf(n, m, Q)), sep=sep, pred=pred))
    return (eps, m), rows


# ---------------------------------------------------------------------------- B. dynamics
def class_diameter(Th, lab, K=3):
    out = 0.0
    for k in range(K):
        idx = lab == k
        if idx.sum() >= 2:
            out = max(out, float(np.sqrt(_cost_matrix(Th[idx], Th[idx], G).max())))
    return out


def run_job(args):
    """All sizes in lockstep (midpoint rule, the cell plan every PLAN_EVERY steps): the nearest-atom leak toward ϑ^b(q), fixed or adaptive."""
    eps, seed, eta = args
    t0 = time.time()
    net, tb, loop, res, fit = setup(eps)
    clouds = CL.mixture_start(tb, M_REF, np.random.default_rng(seed), JITTER, G)
    lab0 = np.argmin(_cost_matrix(clouds, tb, G), axis=1)
    sizes = MS + (M_REF,)
    Th = {m: clouds[:m].copy() for m in sizes}
    Z = {m: tb.copy() for m in sizes}
    apl = {m: AP.AdaptivePriorLoop(loop, eta, plan='nearest') for m in sizes}
    w0 = {m: np.bincount(lab0[:m], minlength=3) / m for m in sizes}
    rec = dict(s=[], h2={m: [] for m in MS}, fit={m: [] for m in sizes}, counts={m: [] for m in sizes})
    gamma = {}
    n = int(round(S_B / DT))
    for k in range(n + 1):
        if k % PLAN_EVERY == 0:
            for m in sizes:
                gamma[m], _ = apl[m].plan(Th[m], Z[m])
                if np.isinf(eta):
                    Z[m] = apl[m].bary(Th[m], gamma[m], Z[m]); gamma[m], _ = apl[m].plan(Th[m], Z[m])
        if k % REC == 0:
            rec['s'].append(k * DT)
            fM = loop.f_hat(Th[M_REF])
            for m in sizes:
                rec['fit'][m].append(loop.fit_error(Th[m]))
                rec['counts'][m].append(np.bincount(np.argmin(_cost_matrix(Th[m], Z[m], G), axis=1), minlength=3))
            for m in MS:
                rec['h2'][m].append(float(np.mean((loop.f_hat(Th[m]) - fM) ** 2)))
        if k == n:
            break
        for m in sizes:
            a = apl[m]
            k1T, k1Z = a.drift(Th[m], Z[m], gamma[m])
            k2T, k2Z = a.drift(Th[m] + 0.5 * DT * k1T, Z[m] + 0.5 * DT * k1Z, gamma[m])
            Th[m] = Th[m] + DT * k2T; Z[m] = Z[m] + DT * k2Z
            if np.isinf(eta):
                Z[m] = a.bary(Th[m], gamma[m], Z[m])
    out = dict(eps=eps, seed=seed, eta=eta, s=np.array(rec['s']), h2={m: np.array(v) for m, v in rec['h2'].items()},
               fit={m: np.array(v) for m, v in rec['fit'].items()}, counts={m: np.array(v) for m, v in rec['counts'].items()}, w0=w0,
               final={m: Th[m].astype(np.float32) for m in sizes}, approx_fit=fit)
    out['dist'], out['diam'] = {}, {}
    for m in sizes:
        tg = gauge(net, tb, w0[m])
        lab = lab0[:m]
        means = np.array([Th[m][lab == j].mean(axis=0) for j in range(3)])
        out['dist'][m] = float(np.linalg.norm(AP.whitened(means - tg, w0[m], G))) if tg is not None else float('nan')
        out['diam'][m] = class_diameter(Th[m], lab)
    print(f"[{time.time() - t0:.0f}s] B eps {eps:g} seed {seed} eta {eta:g}", flush=True)
    return out


# ---------------------------------------------------------------------------- report
def report(A, M, B, t0=None, End=None, Out=None):
    es = lambda e: '∞' if np.isinf(e) else f'{e:.5f}'
    lines = ["# Proposition D‴.28: the adaptive prior's floor removal is a gauge property (`exp_gauge_floor.py`)", "",
             f"The three-bump instance with a fourth bump ε·exp(−(x − {EXTRA[0]:g})²/(2·{EXTRA[1]:g}²)) that K₀ = 3 kernels cannot represent; masses q = 1/3, "
             f"σ₀ = 0.15σ_crit = {SIGMA0:.5f}; ϑ^b(q) the best weighted approximant (a critical point of the unregularised loss, continued in ε).", "",
             "## A. The approximant, the gauge and the sensitivities (no dynamics)", "",
             "Per ε: the approximant's relative fit error and gradient residual; the smallest eigenvalue of the whitened unregularised Hessian Ã = H_L − R̃ and "
             "‖R̃‖/‖H_L‖; the gauge check at six weight vectors 0.01 inside W_g — the largest unregularised weighted gradient at G_wϑ^b(q), the critical point "
             "reached by Newton continued in w from ϑ^b(q) (20 steps; its largest distance from G_wϑ^b(q), and of its f̂ from f̂_b(q)), and in how many of the "
             "six cases plain Newton from ϑ^b(q) at w lands at G_wϑ^b(q) (else at another critical point of L_w, or at none); the adaptive prior's "
             "function-level sensitivity (critical points reached by Newton from ϑ^b(q)) relative to the fixed prior's (central differences); the fixed "
             "matching prior's sensitivity against σ₀Ĵ(Ã + σ₀I)⁻¹g̃ and against D‴.24's form with H_L; the fixed prior's floor law per unit of 1/m − 1/M.", "",
             "| ε | fit error | residual | min eig Ã | ‖R̃‖/‖H_L‖ | ∇L_w at G_wϑ^b | continued from ϑ^b(q): distance, Δf̂ | plain Newton at G_wϑ^b | ‖F_ad‖/‖F_fixed‖ | "
             "F_fixed vs σ₀Ĵ(Ã + σ₀I)⁻¹g̃ | vs H_L form | law (fixed) |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in A:
        gd = max(x['dist'] for x in r['gauge']); gf = max(x['dfh'] for x in r['gauge'])
        gg = max(x.get('grad', float('nan')) for x in r['gauge'])
        npl = sum(1 for x in r['gauge'] if x.get('dist_plain', np.inf) < 1e-8)
        lines.append(f"| {r['eps']:g} | {r['fit']:.3e} | {r['res']:.1e} | {r['ev'][0]:.3e} | {r['R_rel']:.2e} | {gg:.1e} | {gd:.1e}, {gf:.1e} | {npl} of 6 | "
                     f"{r['Fa']:.1e} | {r['rel_gauge']:.1e} | {r['rel_gn']:.1e} | {r['law_fixed']:.4e} |")
    lines += ["", "The residual curvature at ϑ^b(q): the largest entry of its amplitude row and column (zero at a critical point), and λ⁺_k = λ_max(G^{1/2}R_kG^{1/2}) per "
              "class; W_g's thresholds q_k|a^b_k|/a_max.", "", "| ε | amplitude row | largest entry of R | λ⁺_k | W_g: w > |", "|---|---|---|---|---|"]
    for r in A:
        lines.append(f"| {r['eps']:g} | {r['amp_row']:.1e} | {r['R_norm']:.2e} | " + ", ".join(f"{x:.2e}" for x in r['lam']) + " | (" + ", ".join(f"{x:.4f}" for x in r['lo']) + ") |")
    if 'pW' in A[0]:
        lines += ["", "The probability that the sampled weights fall in W_g, P(w^m ∈ W_g) for w^m multinomial(m, q)/m:", "",
                  "| ε | m = " + " | m = ".join(str(m_) for m_ in MS_WG) + " |", "|---|" + "---|" * len(MS_WG)]
        for r in A:
            lines.append(f"| {r['eps']:g} | " + " | ".join(f"{r['pW'][m_]:.4f}" for m_ in MS_WG) + " |")
    if any(r.get('sepchk') for r in A):
        lines += ["", "(iii) at the level of the dynamics: the separation eigenvalues of the averaged nearest-atom adaptive drift (its central-difference Jacobian) at the "
                  "replicated configuration G_wϑ^b(q), counts (6, 4, 5) (w = (0.4, 0.267, 0.333)), against −σ₀ + (q_k/w_k)eig(G^{1/2}R_k(q)G^{1/2}), each n_k − 1 times:", "",
                  "| ε | η | drift there | separations invariant to | largest |difference| | largest eigenvalue (−margin) | predicted |", "|---|---|---|---|---|---|---|"]
        for r in A:
            for eta, c_ in sorted(r.get('sepchk', {}).items(), key=lambda kv: kv[0]):
                lines.append(f"| {r['eps']:g} | {es(eta)} | {c_['drift0']:.1e} | {c_['inv']:.1e} | {c_['dmax']:.1e} | {c_['largest']:+.6e} | {c_['pred_largest']:+.6e} |")
    if End:
        lines += ["", "The branch ϑ^b(q) continued toward the amplitude's saturation (no dynamics): beyond ε ≈ 0.68 the largest amplitude reaches a_max, the smallest "
                  "eigenvalue of Ã vanishes and W_g's threshold for that class reaches q_k = 1/3, so w = q leaves W_g.", "",
                  "| ε | fit error | max_k|a^b_k|/a_max | min eig Ã | max_kλ⁺_k/σ₀ | W_g: w > |", "|---|---|---|---|---|---|"]
        for r in End:
            lines.append(f"| {r['eps']:g} | {r['fit']:.4e} | {r['amax']:.4f} | {r['ev']:.3e} | {r['lam'] / SIGMA0:.4f} | (" + ", ".join(f"{x:.4f}" for x in r['lo']) + ") |")
    if M:
        lines += ["", "The adaptive prior's splitting margins at G_wϑ^b(q) over every count vector of m in W_g (with no empty class): the largest difference from the "
                  "prediction σ₀ − (q_k/w_k)λ⁺_k, the probability of a negative margin (some class can split), and the smallest margin.", "",
                  "| ε | m | count vectors | largest |computed − predicted| | P(negative margin) | smallest margin |", "|---|---|---|---|---|---|"]
        for (eps, m) in sorted(M):
            rows = M[(eps, m)]
            dmax = max(float(np.abs(x['sep'] - x['pred']).max()) for x in rows)
            pneg = sum(x['p'] for x in rows if np.min(x['sep']) <= 0)
            smin = min(float(np.min(x['sep'])) for x in rows)
            lines.append(f"| {eps:g} | {m} | {len(rows)} | {dmax:.1e} | {pneg:.4f} | {smin:+.2e} |")
    if B:
        lines += ["", f"## B. The nearest-atom leak from warm clouds around ϑ^b(q): fixed (η = 0) against adaptive priors ({S_B:g} slow units)", "",
                  "Per ε, η and m: ‖f̂_m − f̂_M‖² at s = 600, 3000, " + f"{S_B:g} (min–max over the seeds), the fit error at the end (min–max) against the approximant's, "
                  "the class means' whitened weighted distance from G_wϑ^b(q) at the end (max), the largest class diameter at the end (max), and the counts conserved.", "",
                  "| ε | η | m | ‖f̂_m − f̂_M‖² at those s | fit error at the end (approximant) | distance from G_wϑ^b | class diameter | counts conserved |",
                  "|---|---|---|---|---|---|---|---|"]
        for eps in sorted(set(r['eps'] for r in B)):
            for eta in (0.0, SIGMA0, np.inf):
                rs = [r for r in B if r['eps'] == eps and r['eta'] == eta]
                if not rs:
                    continue
                for m in MS:
                    h2 = [", ".join(f"{r['h2'][m][int(np.argmin(np.abs(r['s'] - t)))]:.1e}" for t in (600.0, 3000.0, S_B)) for r in rs]
                    h2v = np.array([[r['h2'][m][int(np.argmin(np.abs(r['s'] - t)))] for t in (600.0, 3000.0, S_B)] for r in rs])
                    fits = np.array([r['fit'][m][-1] for r in rs])
                    dist = np.array([r['dist'][m] for r in rs]); diam = np.array([r['diam'][m] for r in rs])
                    cons = all(np.all(r['counts'][m] == r['counts'][m][0]) for r in rs)
                    lines.append(f"| {eps:g} | {es(eta)} | {m} | " + ", ".join(f"{h2v[:, j].min():.1e}–{h2v[:, j].max():.1e}" for j in range(3)) + " | "
                                 f"{fits.min():.4e}–{fits.max():.4e} ({rs[0]['approx_fit']:.4e}) | {np.nanmax(dist):.1e} | {diam.max():.1e} | {'yes' if cons else 'no'} |")
        # the samples by their margin to the boundary of W_g
        Aeps = {r['eps']: r for r in A}
        near, excess = [], {m: [] for m in MS}
        for eps in sorted(set(r['eps'] for r in B)):
            lo = np.asarray(Aeps[eps]['lo'], float)
            for seed in sorted(set(r['seed'] for r in B)):
                byeta = {r['eta']: r for r in B if r['eps'] == eps and r['seed'] == seed}
                for m in MS:
                    w = np.asarray(byeta[0.0]['w0'][m], float)
                    marg = float(np.min(w - lo))
                    ad = [byeta[e] for e in (SIGMA0, np.inf) if e in byeta]
                    if marg < 0.02:
                        near.append((eps, m, seed, w, marg, [byeta[e]['fit'][m][-1] for e in (0.0, SIGMA0, np.inf) if e in byeta],
                                     [r['dist'][m] for r in ad], byeta[0.0]['approx_fit']))
                    else:
                        excess[m] += [r['fit'][m][-1] / r['approx_fit'] - 1 for r in ad]
        lines += ["", "The samples within 0.02 of the boundary of W_g or outside it (the margin min_k(w_k − q_k|a^b_k|/a_max); negative: outside, no gauge transform): "
                  "the fit error at the end under the fixed prior and the adaptive priors (η = σ₀, ∞) against the approximant's, and the adaptive class means' distance "
                  "from G_wϑ^b(q).", "",
                  "| ε | m | seed | w | margin | fit error: fixed, η = σ₀, η = ∞ (approximant) | distance from G_wϑ^b: η = σ₀, ∞ |", "|---|---|---|---|---|---|---|"]
        for eps, m, seed, w, marg, fits, dists, af in near:
            lines.append(f"| {eps:g} | {m} | {seed} | ({', '.join(f'{x:.3f}' for x in w)}) | {marg:+.4f} | " + ", ".join(f"{x:.4e}" for x in fits)
                         + f" ({af:.4e}) | " + ", ".join('—' if not np.isfinite(x) else f"{x:.2e}" for x in dists) + " |")
        if Out:
            lines += ["", "At the samples outside W_g: the minimum of L_w over the closed amplitude box from the clipped gauge transform — the limit of the systems "
                      "that keep ϑ^b(q)'s branch, approached as the saturated amplitude's coordinate runs off — and the gauge transforms with the classes permuted.", "",
                      "| ε | w | the box minimum's fit error | relative to the approximant's | amplitudes at the bound | permuted gauge transforms |",
                      "|---|---|---|---|---|---|"]
            for o in Out:
                lines.append(f"| {o['eps']:g} | ({', '.join(f'{x:.3f}' for x in o['w'])}) | {o['fit_box']:.4e} | {o['fit_box'] / o['fit_b']:.3f} | {o['at_bound']} | "
                             + (", ".join('(' + ', '.join(str(i + 1) for i in pp) + ')' for pp in o['perms']) if o['perms'] else 'none') + " |")
        ratios = {}
        for r in B:
            if r['eta'] != 0.0:
                continue
            i600 = int(np.argmin(np.abs(r['s'] - 600.0)))
            for m in MS:
                ratios.setdefault((r['eps'], m), []).append(r['h2'][m][-1] / r['h2'][m][i600])
        lines += ["", "The fixed prior's floor after s = 600: ‖f̂_m − f̂_M‖²(" + f"{S_B:g}" + ")/‖f̂_m − f̂_M‖²(600) over the seeds, "
                  + "; ".join(f"ε = {e:g}, m = {m}: {min(v):.3f}–{max(v):.3f}" for (e, m), v in sorted(ratios.items())) + "."]
        if any(excess.values()):
            lines += ["", "At the other samples (margin ≥ 0.02) the adaptive priors' fit error at the end exceeds the approximant's by at most "
                      + ", ".join(f"{max(v):.1e} (m = {m}, median {np.median(v):.1e})" for m, v in excess.items() if v) + " (relative): still converging toward the gauge "
                      "transform, slowest at small m, where the samples lie nearest the boundary."]
    lines += ["", "Numbers are produced by `experiments/exp_gauge_floor.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'gauge_floor.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def save_raw(raw, **arrays):
    tmp = raw[:-4] + '.tmp.npz'
    np.savez(tmp, **arrays)
    os.replace(tmp, raw)


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'gauge_floor_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            report(list(z['A']), z['M'].item(), list(z['B']) if 'B' in z.files else None, t0,
                   End=list(z['E']) if 'E' in z.files else None, Out=list(z['O']) if 'O' in z.files else None)
        return
    if '--only-A' in sys.argv:                                              # Part A, the branch's end and the outside samples; M and B kept
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        B = list(save['B'])
        with mp.Pool(N_WORKERS) as pool:
            A = pool.map(part_A_eps, list(EPS), chunksize=1)
            Aeps = {r['eps']: r for r in A}
            samples = []
            for r in B:
                if r['eta'] != 0.0:
                    continue
                for m in MS:
                    w = np.asarray(r['w0'][m], float)
                    if np.min(w - np.asarray(Aeps[r['eps']]['lo'], float)) < 0 and not any(np.allclose(w, w2) and r['eps'] == e2 for e2, w2 in samples):
                        samples.append((r['eps'], w))
            Out = pool.map(outside_job, samples, chunksize=1)
        save.update(A=np.array(A, dtype=object), E=np.array(branch_end(), dtype=object), O=np.array(Out, dtype=object))
        save_raw(raw, **save)
        report(A, save['M'].item(), B, t0, End=list(save['E']), Out=Out)
        return
    jobsB = [(eps, seed, eta) for eps in EPS_B for eta in (0.0, SIGMA0, np.inf) for seed in range(N_SEEDS)]
    with mp.Pool(N_WORKERS) as pool:
        if '--only-B' in sys.argv:                                         # Parts A and the margins from the raw file
            with np.load(raw, allow_pickle=True) as z:
                A, M = list(z['A']), z['M'].item()
                done = list(z['B']) if 'B' in z.files else []
        else:
            A = pool.map(part_A_eps, list(EPS), chunksize=1)
            lam = {r['eps']: r['lam'] for r in A}
            M = dict(pool.map(margins_job, [(eps, m, lam[eps]) for eps in EPS if eps > 0 for m in (15, 60)], chunksize=1))
            save_raw(raw, A=np.array(A, dtype=object), M=np.array(M, dtype=object))
            done = []
        B = None
        if '--no-B' not in sys.argv:
            have = {(r['eps'], r['seed'], r['eta']) for r in done}
            B = list(done)
            for r in pool.imap_unordered(run_job, [j for j in jobsB if j not in have], chunksize=1):
                B.append(r)                                                # saved after every job: an interrupted run keeps its results
                save_raw(raw, A=np.array(A, dtype=object), M=np.array(M, dtype=object), B=np.array(B, dtype=object))
            B.sort(key=lambda r: (r['eps'], r['eta'], r['seed']))
            save_raw(raw, A=np.array(A, dtype=object), M=np.array(M, dtype=object), B=np.array(B, dtype=object))
    report(A, M, B, t0)


if __name__ == '__main__':
    main()
