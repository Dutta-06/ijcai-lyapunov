"""
Proposition D‴.24: the finite-size floor is the leak's bias at the sampled weights; an adaptive prior removes it without a mass constraint
====================================================================================================================================

The averaged slow loop on the three-bump instance (σ₀ = 0.15σ_crit unless stated), the warm clouds of Proposition D‴.12 (m = 15, 60, 240
against M = 510, the m-system the first m rows, four seeds), midpoint rule dt = 0.25, the plan (cell or OT) every 4 steps.

A. The floor as a bias (no dynamics), matching prior, σ₀ = 0.05, 0.15, 0.5 σ_crit: Proposition D‴.9's function-level sensitivity F (central
   differences of the fixed-prior nearest-atom stationary point) against the bias formula Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰ (src/core/floor_bias.py), the
   floor law per unit of 1/m − 1/M, its terms per mode h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j and the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ); the predicted floors at
   m = 15, 60, 240. The admissible weights of the instance (an exact representation needs |A_k|/w_k < a_max) and, per seed and size, the
   sampled weights' margin min_k(w_k − |A_k|/a_max) and the exact representation's residual.
B. Dynamics at the matching prior, the nearest-atom leak toward the fixed prior (η = 0) and toward the adaptive prior (η = σ₀, ∞; the cell plan,
   the cell masses), 6000 slow units: ‖f̂_m − f̂_M‖² and the fit error every 50 slow units, the counts, the class atoms' whitened weighted
   distance from the exact representation of the sampled weights, and at the end W²(μ_m, μ_M) against the W² between the exact
   representations of the two systems' weights (the measure-level distance, dominated by the cross-class transport ½‖w^m − w^M‖₁Δ²).
C. The design table at the misplaced prior (centres +0.15), 3000 slow units: the nearest-atom and the OT leak, each toward the fixed and the
   adaptive prior (η = σ₀; the OT one with the settled gate of Theorem D‴.23): the fit error (the bias) and ‖f̂_m − f̂_M‖² (the floor) at the end.
D. The real closed loop (k_e = 25, σ = σ₀/k_e, the prior's rate η/k_e; RK4 plant with the parameters frozen over the step, Euler parameters,
   dt = 10⁻³, the cells every 100 steps), matching prior, nearest-atom leak, m = 15 and 60 at the sampled weights of seed 0: from the lift of
   the fixed prior's stationary point Θ*(w), toward the fixed prior and toward the adaptive prior (η = ∞), 3000 time units; and from the exact
   representation of w with the adaptive prior at its atoms. The rms tracking error over windows of five periods against (1/k_e)‖f − f̂‖_ρ of
   the averaged loop at s = t/k_e from the same start (`--only D` computes it and merges it).
Checks at the sampled weights (`--only P` computes them from the stored finals and merges them; referee round 18): per seed and size the
   fixed prior's Θ*(w) with Proposition D‴.5's splitting margins (matching and misplaced prior) and ϑ^e(w); per system the end diameter
   of the classes, the coupled distance D of the shared particles and W² split into cross- and within-class transport; the exact bias at the
   sampled weights and the delta method's F(w^m − w^M) against the observed floor; Part C's gate and fit errors at s = 1250, 2000.
E. The instance's probabilities (exact multinomial sums): representability, the fixed prior's saddles, Hoeffding's bound, the adaptive prior's
   expected limiting floor (the constrained minima on the branch of ϑ^e(q) outside W).
F. Monte Carlo (AF_NMC = 600 nested samples): the replicated stationary points' floors and coupled distances against the laws.

Output: results/chaos/adaptive_floor.md and adaptive_floor_raw.npz (AF_SEEDS, AF_WORKERS, AF_SB, AF_SC, AF_TREAL, AF_OUT; `--from-raw` regenerates
the md).
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
from src.core.ot_leak import _cost_matrix, transport_plan

OUT = os.environ.get('AF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
A_BUMP = np.array([b[0] for b in BUMPS])
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q = np.ones(3) / 3
MS, M_REF = (15, 60, 240), 510
JITTER = 0.02
DT, PLAN_EVERY, REC = 0.25, 4, 200                                        # records every 50 slow units
S_B, S_C = float(os.environ.get('AF_SB', 6000.0)), float(os.environ.get('AF_SC', 3000.0))
SHIFT_C = 0.15
N_SEEDS = int(os.environ.get('AF_SEEDS', 4))
N_WORKERS = int(os.environ.get('AF_WORKERS', 2))


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def representation(net, shift=0.0):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def setup(shift=0.0, sigma0=SIGMA0):
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    exact = representation(net, 0.0)
    prior = representation(net, shift)
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, sigma0, prior, n_quad=N_Q)
    return net, exact, prior, loop


# ---------------------------------------------------------------------------- A. the floor as a bias
def part_A():
    out = dict(sig=[])
    for frac in (0.05, 0.15, 0.5):
        net, exact, prior, loop = setup(0.0, frac * SIGMA_CRIT)
        law = FB.bias_floor_law(loop, Q, exact, frac * SIGMA_CRIT)
        S = PD.weight_sensitivity_prior(loop.nearest, Q, exact)
        F2 = PD.function_sensitivity(loop.nearest, Q, exact, S) @ np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
        numeric = float(np.trace(F2 @ law['C2'] @ F2.T)) / F2.shape[0]
        # the parameter level: the fixed prior's sensitivity S (whitened weighted) against H_L(H_L + σ₀I)⁻¹S^e, S^e = ∂_wϑ^e (the gauge's);
        # Proposition D‴.5's coupled-distance law per unit of 1/m − 1/M for both
        V2 = np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T
        Sy = np.column_stack([(np.sqrt(Q)[:, None] * np.einsum('kcj,j->kc', S, V2[:, c]) / np.sqrt(G)[None, :]).ravel() for c in range(2)])
        Sey = -law['dx0']
        split = law['HL'] @ np.linalg.solve(law['HL'] + frac * SIGMA_CRIT * np.eye(9), Sey)
        out['sig'].append(dict(frac=frac, sigma0=frac * SIGMA_CRIT, rel=float(np.linalg.norm(law['F_bias'] - F2) / np.linalg.norm(F2)),
                               law=law['law'], numeric=numeric, bound=law['bound'], h=law['h'], terms=law['terms'],
                               split_rel=float(np.linalg.norm(Sy - split) / np.linalg.norm(Sy)),
                               D_fixed=float(np.trace(Sy @ law['C2'] @ Sy.T)), D_gauge=float(np.trace(Sey @ law['C2'] @ Sey.T))))
    # admissibility of the sampled weights and the exact representations
    net, exact, prior, loop = setup(0.0)
    adm = {}
    for seed in range(N_SEEDS):
        clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
        for m in MS + (M_REF,):
            counts = np.bincount(np.argmin(_cost_matrix(clouds[:m], exact, G), axis=1), minlength=3)
            w = counts / m
            margin = float(np.min(w - np.abs(A_BUMP) / net.a_max))
            if margin > 0:
                te, res, fit = FB.exact_weighted_atoms(loop, w, exact)
            else:
                res, fit = float('nan'), float('nan')
            adm[(seed, m)] = dict(counts=counts, margin=margin, res=res, fit=fit)
    out['adm'] = adm
    return out


# ---------------------------------------------------------------------------- B, C. dynamics
def class_deviation(Th, Z, w_exact_atoms, w):
    """The whitened weighted norm of the cell means' deviation from the exact representation of the weights w (cells of Z)."""
    lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
    means = np.array([Th[lab == k].mean(axis=0) if np.any(lab == k) else w_exact_atoms[k] for k in range(3)])
    return float(np.linalg.norm(AP.whitened(means - w_exact_atoms, w, G)))


def run_job(args):
    part, seed, plan, eta, gate, shift, S_end = args
    t0 = time.time()
    net, exact, prior, loop = setup(shift)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
    sysm = {m: clouds[:m].copy() for m in MS + (M_REF,)}
    out = dict(part=part, seed=seed, plan=plan, eta=eta, gate=gate, shift=shift, rec={}, final={}, t_open={})
    for m, Th0 in sysm.items():
        counts0 = np.bincount(np.argmin(_cost_matrix(Th0, exact, G), axis=1), minlength=3)
        w0 = counts0 / m
        te = FB.exact_weighted_atoms(loop, w0, exact)[0] if np.min(w0 - np.abs(A_BUMP) / net.a_max) > 0 else exact
        apl = AP.AdaptivePriorLoop(loop, eta, plan=plan)
        rec = dict(s=[], fit=[], dev=[], counts=[])

        def cb(s, Th, Z, adapting, rec=rec, te=te, w0=w0):
            rec['s'].append(s); rec['fit'].append(loop.fit_error(Th)); rec['dev'].append(class_deviation(Th, Z, te, w0))
            rec['counts'].append(np.bincount(np.argmin(_cost_matrix(Th, Z, G), axis=1), minlength=3))
        Th, Z, t_open = apl.run(Th0, prior.copy(), S_end, DT, plan_every=PLAN_EVERY, rec_every=REC, callback=cb, gate=gate)
        out['rec'][m] = {k: np.array(v) for k, v in rec.items()}
        out['final'][m] = Th.astype(np.float32)
        out['t_open'][m] = t_open
    # the function-level distance to the reference, along the records: rerun is not needed — the f̂'s at the records are not stored, so the
    # distance is computed from the stored finals at the end and, along the way, from the fits (both systems converge to f when exact)
    fM = loop.f_hat(np.asarray(out['final'][M_REF], float))
    out['h2_end'] = {m: float(np.mean((loop.f_hat(np.asarray(out['final'][m], float)) - fM) ** 2)) for m in MS}
    print(f"[{time.time() - t0:.0f}s] {part} seed {seed} {plan} eta {eta:g} gate {gate}", flush=True)
    return out


def run_job_h2(args):
    """As run_job, but recording ‖f̂_m − f̂_M‖² along the way (all sizes integrated in lockstep)."""
    part, seed, plan, eta, gate, shift, S_end = args
    t0 = time.time()
    net, exact, prior, loop = setup(shift)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
    sizes = MS + (M_REF,)
    Th = {m: clouds[:m].copy() for m in sizes}
    Z = {m: prior.copy() for m in sizes}
    apl = {m: AP.AdaptivePriorLoop(loop, eta, plan=plan) for m in sizes}
    te, w0 = {}, {}
    for m in sizes:
        c0 = np.bincount(np.argmin(_cost_matrix(Th[m], exact, G), axis=1), minlength=3); w0[m] = c0 / m
        te[m] = FB.exact_weighted_atoms(loop, w0[m], exact)[0] if np.min(w0[m] - np.abs(A_BUMP) / net.a_max) > 0 else exact
    rec = dict(s=[], h2={m: [] for m in MS}, fit={m: [] for m in sizes}, dev={m: [] for m in sizes}, counts={m: [] for m in sizes})
    gamma = {}
    n = int(round(S_end / DT))
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
                rec['fit'][m].append(loop.fit_error(Th[m])); rec['dev'][m].append(class_deviation(Th[m], Z[m], te[m], w0[m]))
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
    out = dict(part=part, seed=seed, plan=plan, eta=eta, gate=gate, shift=shift, s=np.array(rec['s']),
               h2={m: np.array(v) for m, v in rec['h2'].items()}, fit={m: np.array(v) for m, v in rec['fit'].items()},
               dev={m: np.array(v) for m, v in rec['dev'].items()}, counts={m: np.array(v) for m, v in rec['counts'].items()},
               final={m: Th[m].astype(np.float32) for m in sizes}, w0=w0)
    # the parameter level at the end: W²(μ_m, μ_M) against W² between the exact representations of the two weight vectors
    out['W2_end'] = {m: OT.W2_empirical(Th[m], Th[M_REF], G) for m in MS}
    out['W2_exact'] = {m: float(np.sum(transport_plan(_cost_matrix(te[m], te[M_REF], G), w0[m], w0[M_REF]) * _cost_matrix(te[m], te[M_REF], G))) for m in MS}
    print(f"[{time.time() - t0:.0f}s] {part} seed {seed} {plan} eta {eta:g}", flush=True)
    return out


# ---------------------------------------------------------------------------- D. the real closed loop
KE, DT_REAL, T_REAL = 25.0, 1e-3, float(os.environ.get('AF_TREAL', 3000.0))


def x_ref_dot(t):
    return np.array([2.0 * np.cos(t)])


def real_loop(net, Th0, Z0, eta, T, adaptive):
    """The real closed loop with the nearest-atom leak toward Z (fixed, or adaptive at the rate η/k_e; η = ∞: the cell means): e² every 0.1 time
    units, the final state."""
    x = x_ref(0.0).copy(); Th = np.array(Th0, float); Z = np.array(Z0, float)
    sig = SIGMA0 / KE
    n = int(round(T / DT_REAL)); t = 0.0
    ts, e2 = [], []
    for i in range(n):
        if i % 100 == 0:
            lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
            if adaptive and np.isinf(eta):
                for k in range(Z.shape[0]):
                    if np.any(lab == k):
                        Z[k] = Th[lab == k].mean(axis=0)
            tgt = Z[lab]
            ts.append(t); e2.append(float((x - x_ref(t))[0] ** 2))
        e = x - x_ref(t)
        xdot = lambda y, tau: f_vec(y) + x_ref_dot(tau) - KE * (y - x_ref(tau)) - net.f_hat(y, Th)
        k1 = xdot(x, t); k2 = xdot(x + 0.5 * DT_REAL * k1, t + 0.5 * DT_REAL); k3 = xdot(x + 0.5 * DT_REAL * k2, t + 0.5 * DT_REAL)
        k4 = xdot(x + DT_REAL * k3, t + DT_REAL)
        Thn = Th + DT_REAL * (G[None, :] * net.grad_g(x, e, Th) - sig * (Th - tgt))
        if adaptive and np.isfinite(eta) and eta > 0:
            means = np.array([Th[lab == k].mean(axis=0) if np.any(lab == k) else Z[k] for k in range(Z.shape[0])])
            Z = Z - DT_REAL * (eta / KE) * (Z - means)
        x = x + DT_REAL / 6 * (k1 + 2 * k2 + 2 * k3 + k4); Th = Thn
        t = (i + 1) * DT_REAL
    return np.array(ts), np.array(e2), Th, Z


def part_D_job(args):
    m, case = args
    t0 = time.time()
    net, exact, prior, loop = setup(0.0)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(0), JITTER, G)
    counts = np.bincount(np.argmin(_cost_matrix(clouds[:m], exact, G), axis=1), minlength=3); w = counts / m
    Tw, res = PD.weighted_stationary_atoms_prior(loop.nearest, w, exact, w0=Q, n_cont=4)    # the fixed prior's stationary point at w
    te, rese, fite = FB.exact_weighted_atoms(loop, w, exact)                                  # the exact representation at w
    if case == 'exact':
        Th0, Z0, adaptive, eta = np.repeat(te, counts, axis=0), te.copy(), True, np.inf
    else:
        Th0, Z0, adaptive, eta = np.repeat(Tw, counts, axis=0), exact.copy(), case == 'adaptive', (np.inf if case == 'adaptive' else 0.0)
    ts, e2, Th, Z = real_loop(net, Th0, Z0, eta, T_REAL, adaptive)
    # the averaged loop from the same start, at s = t/k_e
    apl = AP.AdaptivePriorLoop(loop, eta if adaptive else 0.0, plan='nearest')
    rec = dict(s=[], fit=[])
    apl.run(Th0, Z0, T_REAL / KE, 0.05, plan_every=1, rec_every=20, callback=lambda s_, T_, Z_, ad: (rec['s'].append(s_), rec['fit'].append(loop.fit_error(T_))))
    normF = float(np.sqrt(np.mean(loop.F ** 2)))
    print(f"[{time.time() - t0:.0f}s] D m {m} {case}", flush=True)
    return dict(m=m, case=case, counts=counts, ts=ts, e2=e2, s_avg=np.array(rec['s']), fit_avg=np.array(rec['fit']), normF=normF,
                res_fixed=res, res_exact=rese, fit_exact=fite, final=Th.astype(np.float32))


def post_B(B):
    """Part B's parameter-level columns from the stored finals: the exact representation of each system's sampled weights (none on the
    boundary of W), W² between the exact representations of the m-system's and the reference's weights, and the end
    distance of the class means (classes by the exact atoms' cells) from the exact representation."""
    net, exact, prior, loop = setup(0.0)
    cache = {}
    for r in B:
        te = {}
        for m in MS + (M_REF,):
            w = np.asarray(r['w0'][m], float)
            key = tuple(np.round(w * m).astype(int))
            if key not in cache:
                cache[key] = FB.exact_weighted_atoms(loop, w, exact) if np.min(w - np.abs(A_BUMP) / net.a_max) > 0 else None
            te[m] = cache[key]
        r['W2_exact'], r['dev_end'] = {}, {}
        for m in MS:
            ok = te[m] is not None and te[M_REF] is not None and te[m][2] < 1e-10 and te[M_REF][2] < 1e-10
            if ok:
                C = _cost_matrix(te[m][0], te[M_REF][0], G)
                r['W2_exact'][m] = float(np.sum(transport_plan(C, r['w0'][m], r['w0'][M_REF]) * C))
                r['dev_end'][m] = class_deviation(np.asarray(r['final'][m], float), exact, te[m][0], np.asarray(r['w0'][m], float))
            else:
                r['W2_exact'][m] = float('nan'); r['dev_end'][m] = float('nan')
    return B


# ---------------------------------------------------------------------------- the checks at the sampled weights (referee round 18)
P_MIN = {15: 0.0, 60: 1e-12, 240: 1e-9}                                    # the count vectors enumerated: multinomial probability above P_MIN
N_MC = int(os.environ.get('AF_NMC', 600))
THR = np.array([10, 7, 8])                                                 # |A_k|/a_max = THR_k/40 exactly


def representability(counts):
    """+1 if w = counts/m lies in W (w_k > |A_k|/a_max for every k: an exact representation on the branch of ϑ^e(q)), 0 on its boundary,
    −1 outside; exact integer arithmetic (40n_k against THR_k·m)."""
    n = np.asarray(counts, int); g = 40 * n - THR * int(n.sum())
    return 1 if np.all(g > 0) else (0 if np.all(g >= 0) else -1)


def k_atom_hessian_min(nloop, w, Tw, h=1e-6):
    """The smallest eigenvalue of the whitened weighted Hessian of the K-atom energy 𝒢_w at Tw (Proposition D‴.5's collective block)."""
    K = Tw.shape[0]; d = 3 * K
    g = lambda v: PD.weighted_energy_grad_prior(nloop, w, v.reshape(K, 3)).ravel()
    J = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        J[:, j] = (g(Tw.ravel() + e) - g(Tw.ravel() - e)) / (2 * h)
    sw = np.repeat(np.sqrt(w), 3) * np.tile(np.sqrt(G), K); sv = np.tile(np.sqrt(G), K) / np.repeat(np.sqrt(w), 3)
    H = sw[:, None] * J * sv[None, :]
    return float(np.linalg.eigvalsh(0.5 * (H + H.T))[0])


def n_continuation(w):
    return max(4, int(np.ceil(np.abs(np.asarray(w, float) - Q).sum() / 0.02)))


def fixed_point(loop, counts, start):
    """The fixed prior's replicated stationary point at the counts (nearest-atom leak toward loop's prior): Θ*(w) by Newton continued from q
    at `start`, its residual, Proposition D‴.5's splitting margins σ₀ − λ_max R^G_k (+∞ for a class with fewer than two particles: it has no
    separation mode), the smallest eigenvalue of the K-atom Hessian, whether every atom lies in its own cell, f̂ on the grid, and stability
    (D‴.5-admissibility: all of these positive). Every class must be non-empty."""
    counts = np.asarray(counts, int); w = counts / counts.sum()
    Tw, res = PD.weighted_stationary_atoms_prior(loop.nearest, w, start, w0=Q, n_cont=n_continuation(w))
    a, c, s_, d_, phi = loop.features(Tw)
    sep = np.where(counts >= 2, loop.sigma0 - loop.kappa_R(loop.F - phi @ (w * a), Tw), np.inf)
    hmin = k_atom_hessian_min(loop.nearest, w, Tw)
    own = bool(np.all(loop.nearest.assignment(Tw) == np.arange(Tw.shape[0])))
    return dict(Tw=Tw, res=float(res), sep=sep, hess_min=hmin, own=own, fhat=phi @ (w * a),
                stable=bool(np.all(sep > 0) and hmin > 0 and own and res < 1e-9))


def exact_point(loop, counts, start):
    """ϑ^e(w) on the branch of ϑ^e(q), with f̂ on the grid; None outside W or if the continuation does not reach an exact representation."""
    counts = np.asarray(counts, int); w = counts / counts.sum()
    if representability(counts) <= 0:
        return None
    te, res, fit = FB.exact_weighted_atoms(loop, w, start, n_cont=n_continuation(w))
    if fit > 1e-10:
        return None
    a, c, s_, d_, phi = loop.features(te)
    return dict(te=te, res=res, fit=fit, fhat=phi @ (w * a))


def cloud_labels(seed):
    """The classes of the seed's warm cloud: the cells of the exact atoms at s = 0 (the m-system is the first m rows)."""
    net, exact, prior, loop = setup(0.0)
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
    return np.argmin(_cost_matrix(clouds, exact, G), axis=1)


def class_spread(Th, lab, K=3):
    """The largest G⁻¹-diameter of a class, and the sizes of that class's sub-clusters (single linkage at 0.05 in the G⁻¹ metric)."""
    from scipy.cluster.hierarchy import linkage, fcluster
    best, sizes = -1.0, ()
    for k in range(K):
        idx = lab == k
        if idx.sum() < 2:
            continue
        dmax = float(np.sqrt(_cost_matrix(Th[idx], Th[idx], G).max()))
        if dmax > best:
            cl = fcluster(linkage(Th[idx] / np.sqrt(G)[None, :], 'single'), t=0.05, criterion='distance')
            best, sizes = dmax, tuple(int(x) for x in sorted(np.bincount(cl)[1:], reverse=True))
    return max(best, 0.0), sizes


def coupled_D(Th_m, Th_M):
    """Proposition D‴.5's coupled distance (1/m)Σ_{i≤m}‖ϑ_i^m − ϑ_i^M‖²_{G⁻¹} of the shared particles."""
    m = Th_m.shape[0]
    return float(np.mean(np.sum((Th_m - Th_M[:m]) ** 2 / G[None, :], axis=1)))


def W2_parts(Th_m, lab_m, Th_M, lab_M):
    """W²(μ_m, μ_M) by an optimal plan, split into the transport between different classes and within classes; the cross-class mass
    against its minimum ½‖w^m − w^M‖₁."""
    m, M = len(Th_m), len(Th_M)
    C = _cost_matrix(Th_m, Th_M, G)
    gam = transport_plan(C, np.full(m, 1.0 / m), np.full(M, 1.0 / M))
    cross = lab_m[:, None] != lab_M[None, :]
    wm, wM = np.bincount(lab_m, minlength=3) / m, np.bincount(lab_M, minlength=3) / M
    return dict(W2=float(np.sum(gam * C)), cross=float(np.sum(gam * C * cross)), within=float(np.sum(gam * C * ~cross)),
                cross_mass=float(np.sum(gam * cross)), half_l1=0.5 * float(np.abs(wm - wM).sum()))


def replicated_D(P_m, P_M, w_m, key):
    """The coupled distance of two replicated configurations with the same classes: Σ_kw^m_k‖ϑ_k(w^m) − ϑ_k(w^M)‖²_{G⁻¹}."""
    return float(np.sum(w_m * np.sum((P_m[key] - P_M[key]) ** 2 / G[None, :], axis=1)))


def function_sensitivity_q():
    """Proposition D‴.9's F at q (σ₀ = 0.15σ_crit) along e₁ − e₃, e₂ − e₃ (N, 2)."""
    net, exact, prior, loop = setup(0.0)
    S = PD.weight_sensitivity_prior(loop.nearest, Q, exact)
    return PD.function_sensitivity(loop.nearest, Q, exact, S) @ np.array([[1.0, 0.0, -1.0], [0.0, 1.0, -1.0]]).T


def post_checks(B, C, D):
    """From the stored finals: per seed and size the fixed prior's Θ*(w) and its D‴.5 margins at the matching prior (Parts B, D) and at the
    misplaced one (Part C), and ϑ^e(w); per system of Parts B and C the end spread, the relabelled particles (end cells of the exact atoms
    against the initial classes), the coupled distance D and W² split by classes; in Part B the replicated stationary points' D and
    ‖f̂_m − f̂_M‖² (the exact bias at the sampled weights) and the delta method's F(w^m − w^M); in Part C the gate's opening time and the fit
    errors at s = 1250 and 2000."""
    net, exact, prior, loop = setup(0.0)
    loopC = setup(SHIFT_C)[3]
    F2 = function_sensitivity_q()
    labels = {seed: cloud_labels(seed) for seed in range(N_SEEDS)}
    fp, fpC, ex = {}, {}, {}
    for seed in range(N_SEEDS):
        for m in MS + (M_REF,):
            counts = np.bincount(labels[seed][:m], minlength=3)
            fp[(seed, m)] = fixed_point(loop, counts, exact)
            fpC[(seed, m)] = fixed_point(loopC, counts, exact)
            ex[(seed, m)] = exact_point(loop, counts, exact)
    outB, outC = {}, {}
    for r in B:
        seed, lab = r['seed'], labels[r['seed']]
        ThM = np.asarray(r['final'][M_REF], float)
        cellsM = np.argmin(_cost_matrix(ThM, exact, G), axis=1)
        for m in MS:
            Thm = np.asarray(r['final'][m], float)
            cells = np.argmin(_cost_matrix(Thm, exact, G), axis=1)
            wm, wM = np.bincount(lab[:m], minlength=3) / m, np.bincount(lab, minlength=3) / M_REF
            spread, sizes = class_spread(Thm, lab[:m])
            o = dict(spread=spread, sizes=sizes, relabelled=int(np.sum(cells != lab[:m])), D=coupled_D(Thm, ThM), W2=W2_parts(Thm, cells, ThM, cellsM),
                     h2_obs=float(r['h2'][m][-1]), D_fixed_star=replicated_D(fp[(seed, m)], fp[(seed, M_REF)], wm, 'Tw'),
                     h2_exact=float(np.mean((fp[(seed, m)]['fhat'] - fp[(seed, M_REF)]['fhat']) ** 2)), h2_delta=float(np.mean((F2 @ (wm - wM)[:2]) ** 2)))
            o['D_exact_star'] = replicated_D(ex[(seed, m)], ex[(seed, M_REF)], wm, 'te') if ex[(seed, m)] is not None and ex[(seed, M_REF)] is not None else float('nan')
            outB[(r['eta'], seed, m)] = o
        print(f"checks B eta {r['eta']:g} seed {seed}", flush=True)
    for r in C:
        seed, lab = r['seed'], labels[r['seed']]
        ThM = np.asarray(r['final'][M_REF], float)
        cellsM = np.argmin(_cost_matrix(ThM, exact, G), axis=1)
        for m in MS:
            Thm = np.asarray(r['final'][m], float)
            cells = np.argmin(_cost_matrix(Thm, exact, G), axis=1)
            spread, sizes = class_spread(Thm, cells)
            rc = r['rec'][m]
            outC[(r['plan'], r['eta'], r['gate'], seed, m)] = dict(
                spread=spread, sizes=sizes, relabelled=int(np.sum(cells != lab[:m])), D=coupled_D(Thm, ThM), W2=W2_parts(Thm, cells, ThM, cellsM),
                t_open=r['t_open'][m], fit_at={t: float(rc['fit'][int(np.argmin(np.abs(rc['s'] - t)))]) for t in (1250.0, 2000.0)})
        print(f"checks C {r['plan']} eta {r['eta']:g} seed {seed}", flush=True)
    strip = lambda d: {k: v for k, v in d.items() if k != 'fhat'} if d is not None else None
    return dict(fp={k: strip(v) for k, v in fp.items()}, fpC={k: strip(v) for k, v in fpC.items()}, ex={k: strip(v) for k, v in ex.items()},
                B=outB, C=outC)


def constrained_floor(loop, w):
    """min ⟨(f − Σ_kw_ka_kφ(·;c_k, s_k))²⟩_ρ over |a_k| ≤ a_max, |c_k| ≤ c_max, s_k ∈ [s_min, s_max] (L-BFGS-B from the bumps, the amplitudes
    A_k/w_k clipped): at weights outside W the constrained minimum of the unregularised weighted loss on the branch of ϑ^e(q), taken as the
    adaptive prior's limiting error (the dynamics from the warm clouds follows the branch). It is a local minimum, not the infimum: where a
    permutation of the classes over the bumps is representable the global infimum is 0 (referee round 19, D3n-10)."""
    from scipy.optimize import minimize
    X, F = loop.X, loop.F
    kw = NET_KW
    w = np.asarray(w, float)

    def obj(v):
        a, c, s = v[:3], v[3:6], v[6:]
        dx = X[:, None] - c[None, :]
        phi = np.exp(-dx ** 2 / (2 * s[None, :] ** 2))
        r = F - phi @ (w * a)
        wa = (w * a)[None, :]
        g = -2 * np.concatenate([np.mean(r[:, None] * phi * w[None, :], axis=0), np.mean(r[:, None] * wa * phi * dx / s[None, :] ** 2, axis=0),
                                 np.mean(r[:, None] * wa * phi * dx ** 2 / s[None, :] ** 3, axis=0)])
        return float(np.mean(r ** 2)), g
    a0 = np.clip(A_BUMP / np.maximum(w, 1e-12), -kw['a_max'], kw['a_max'])
    v0 = np.concatenate([a0, [b[1] for b in BUMPS], [b[2] for b in BUMPS]])
    bounds = [(-kw['a_max'], kw['a_max'])] * 3 + [(-kw['c_max'], kw['c_max'])] * 3 + [(kw['s_min'], kw['s_max'])] * 3
    res = minimize(obj, v0, jac=True, method='L-BFGS-B', bounds=bounds, options=dict(ftol=1e-15, gtol=1e-12, maxiter=5000))
    return float(res.fun)


def enumerate_job(m):
    """Every count vector of size m with multinomial probability above P_MIN[m]: its probability, its representability (also with the classes
    permuted over the bumps), the fixed matching prior's D‴.5 stability of Θ*(w), and outside W the constrained minimum of the loss on the
    branch."""
    from scipy.stats import multinomial
    from itertools import permutations
    t0 = time.time()
    net, exact, prior, loop = setup(0.0)
    rows = []
    for n1 in range(m + 1):
        for n2 in range(m + 1 - n1):
            n = np.array([n1, n2, m - n1 - n2]); p = float(multinomial.pmf(n, m, Q))
            if p <= P_MIN[m]:
                continue
            rep = representability(n)
            perm = any(representability(n[list(pi)]) > 0 for pi in permutations(range(3)))
            row = dict(n=n, p=p, rep=rep, perm=perm, empty=bool(np.any(n == 0)), stable=None, res=float('nan'), sep_min=float('nan'), floor=0.0)
            if not row['empty']:
                f = fixed_point(loop, n, exact)
                row.update(stable=f['stable'], res=f['res'], sep_min=float(np.min(f['sep'])), hess_min=f['hess_min'])
            if rep < 0:
                row['floor'] = constrained_floor(loop, n / m)
            rows.append(row)
    print(f"[{time.time() - t0:.0f}s] enumerated m = {m}: {len(rows)} count vectors", flush=True)
    return m, rows


def mc_job(args):
    """Nested samples i0, …, i0 + n − 1 (M = 510 labels i.i.d. uniform from the seed 10⁴ + i, the m-system the first m): the replicated
    stationary points' ‖f̂_m − f̂_M‖² and coupled distance under the fixed prior (Θ*(w)) and the adaptive prior (ϑ^e(w)), with the flags."""
    i0, n = args
    t0 = time.time()
    net, exact, prior, loop = setup(0.0)
    cfp, cex, out = {}, {}, []
    for i in range(i0, i0 + n):
        lab = np.random.default_rng(10_000 + i).integers(0, 3, M_REF)
        cnt = {m: tuple(int(x) for x in np.bincount(lab[:m], minlength=3)) for m in MS + (M_REF,)}
        for c in cnt.values():
            if c not in cfp:
                cfp[c] = fixed_point(loop, c, exact) if min(c) > 0 else None
                cex[c] = exact_point(loop, c, exact)
        rec = {}
        cM = cnt[M_REF]
        for m in MS:
            c = cnt[m]; wm = np.array(c) / m
            Pf, PfM, Pe, PeM = cfp[c], cfp[cM], cex[c], cex[cM]
            rec[m] = dict(counts=c, ok_f=Pf is not None and Pf['res'] < 1e-9 and PfM['res'] < 1e-9,
                          stable=bool(Pf is not None and Pf['stable'] and PfM['stable']), rep=representability(c), rep_M=representability(cM))
            if rec[m]['ok_f']:
                rec[m].update(h2=float(np.mean((Pf['fhat'] - PfM['fhat']) ** 2)), D_f=replicated_D(Pf, PfM, wm, 'Tw'))
            if Pe is not None and PeM is not None:
                rec[m].update(D_e=replicated_D(Pe, PeM, wm, 'te'))
        out.append(rec)
    print(f"[{time.time() - t0:.0f}s] Monte Carlo samples {i0}–{i0 + n - 1}", flush=True)
    return out


def checks(B, C, D, pool):
    """The checks above: post_checks, the enumeration for m = 15, 60, 240 and N_MC nested Monte Carlo samples."""
    enum = dict(pool.map(enumerate_job, list(MS)))
    chunks = [(i, min(50, N_MC - i)) for i in range(0, N_MC, 50)]
    mc = [r for part in pool.map(mc_job, chunks, chunksize=1) for r in part]
    return dict(post=post_checks(B, C, D), enum=enum, mc=mc)


def _rng(x, f='{:.2e}'):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    return '—' if x.size == 0 else (f.format(x.min()) + '–' + f.format(x.max()) if x.min() != x.max() else f.format(x.min()))


def report_A_checks(A, P):
    po = P['post']
    lines = ["", "The weight h_jσ₀²/(h_j + σ₀)² is a band-pass centred at h = σ₀: it is at least half its peak σ₀/4 for h/σ₀ ∈ [3 − 2√2, 3 + 2√2] = [0.17, 5.83]. "
             "The shares of the modes with the weight and without it (z_jᵀC̃z_j alone: the bound's distribution, independent of σ₀):", "",
             "| σ₀/σ_crit | modes 2–6: weighted share | unweighted | modes with h_j/σ₀ ∈ [1/5, 5]: weighted share | unweighted | the weight relative to its peak, modes 1–6 |",
             "|---|---|---|---|---|---|"]
    for r in A['sig']:
        h, sg = np.asarray(r['h']), r['sigma0']
        wgt = h * sg ** 2 / (h + sg) ** 2
        un = r['terms'] / wgt; un = un / un.sum(); tw = r['terms'] / r['terms'].sum()
        near = (h / sg >= 0.2) & (h / sg <= 5.0)
        lines.append(f"| {r['frac']:g} | {tw[1:6].sum():.3f} | {un[1:6].sum():.3f} | {tw[near].sum():.3f} | {un[near].sum():.3f} | "
                     + ", ".join(f"{x:.2f}" for x in (4 * wgt / sg)[:6]) + " |")
    lines += ["", "Since S^e = −∂_wx⁰, the bound is σ₀/4 times the adaptive prior's coupled-distance law per unit of 1/m − 1/M, tr(∂_wx⁰C̃∂_wx⁰ᵀ) = "
              f"{A['sig'][0]['D_gauge']:.4f} (independent of σ₀): the fixed prior's function-level law is at most σ₀/4 times the adaptive prior's parameter-level law."]
    lines += ["", "The fixed matching prior's replicated stationary point Θ*(w) at the sampled weights (σ₀ = 0.15σ_crit): Proposition D‴.5's splitting margin "
              "min_k(σ₀ − λ_max R^G_k) over the classes with at least two particles (the class attaining it), and the smallest eigenvalue of the whitened weighted "
              "K-atom Hessian. Θ*(w) is a strict local minimiser (w ∈ W_adm, \"D‴.5-admissible\") iff both are positive; otherwise it is a saddle and a cluster can split. "
              "This is a different set from W, where an exact representation exists (\"representable\").", "",
              "| m | " + " | ".join(f"seed {s}" for s in range(N_SEEDS)) + " |", "|---|" + "---|" * N_SEEDS]
    for m in MS + (M_REF,):
        cells = []
        for s in range(N_SEEDS):
            f = po['fp'][(s, m)]
            k = int(np.argmin(f['sep']))
            cells.append(f"{f['sep'][k]:+.2e} (class {k + 1}), {f['hess_min']:.1e}" + ("" if f['stable'] else " **saddle**"))
        lines.append(f"| {m} | " + " | ".join(cells) + " |")
    return lines


def report_B_checks(B, P):
    po = P['post']
    lines = ["", "### The fixed prior at the sampled weights: stability, splitting, and the floor as the exact bias", "",
             "Per seed and size at s = 6000: the D‴.5 margin of the fixed prior's Θ*(w) (its class), the largest G⁻¹-diameter of a class at the end with the "
             "sub-cluster sizes of that class (single linkage at 0.05) for the fixed prior and the diameters for the adaptive priors; the observed ‖f̂_m − f̂_M‖² of the "
             "fixed prior against the exact bias at the sampled weights, ‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖² (Proposition D‴.24(i) without the linearisation), and against "
             "the delta method's ‖F(w^m − w^M)‖² (Proposition D‴.9's law per sample).", "",
             "| seed | m | counts | D‴.5 margin of Θ*(w) | end diameter, fixed (sub-clusters) | adaptive η = σ₀, ∞ | ‖f̂_m − f̂_M‖², fixed | exact bias (observed/exact) | delta method (observed/delta) |",
             "|---|---|---|---|---|---|---|---|---|"]
    rs = {(r['eta'], r['seed']): r for r in B}
    ratios = dict(stable=[], saddle=[], d_stable=[], d_saddle=[])
    for seed in range(N_SEEDS):
        for m in MS:
            f = po['fp'][(seed, m)]; k = int(np.argmin(f['sep']))
            o0, o1, o2 = po['B'][(0.0, seed, m)], po['B'][(SIGMA0, seed, m)], po['B'][(np.inf, seed, m)]
            c = tuple(int(x) for x in np.round(rs[(0.0, seed)]['w0'][m] * m))
            ex_r = o0['h2_obs'] / o0['h2_exact']; de_r = o0['h2_obs'] / o0['h2_delta']
            stable = f['stable'] and po['fp'][(seed, M_REF)]['stable']
            ratios['stable' if stable else 'saddle'].append(ex_r); ratios['d_stable' if stable else 'd_saddle'].append(de_r)
            lines.append(f"| {seed} | {m} | {c} | {f['sep'][k]:+.2e} (class {k + 1}){'' if f['stable'] else ' saddle'} | {o0['spread']:.2f} {o0['sizes']} | "
                         f"{o1['spread']:.1e}, {o2['spread']:.1e} | {o0['h2_obs']:.3e} | {o0['h2_exact']:.3e} ({ex_r:.4f}) | {o0['h2_delta']:.3e} ({de_r:.2f}) |")
    n_sad = len(ratios['saddle'])
    lines += ["", f"{n_sad} of {n_sad + len(ratios['stable'])} fixed-prior m-systems sit at weights where Θ*(w) is a saddle; their observed floors are "
              f"{_rng(ratios['saddle'], '{:.3f}')} of the replicated point's (they have split), against {_rng(ratios['stable'], '{:.4f}')} at the stable ones. "
              f"The delta method's prediction misses the observed floor by the factors {_rng(ratios['d_stable'], '{:.2f}')} at the stable samples "
              f"({_rng(ratios['d_saddle'], '{:.2f}')} at the saddles).",
              "", "### The parameter level: the coupled distance D and W²", "",
              "Per seed and size at s = 6000: Proposition D‴.5's coupled distance D of the shared particles, observed and at the replicated stationary points with the "
              "same classes (the fixed prior's Θ*(w); the adaptive prior's exact representations ϑ^e(w), —: one of the two weight vectors outside W); W²(μ_m, μ_M) by an "
              "optimal plan split into the transport between different classes (the cells of the exact atoms) and within classes, and ½‖w^m − w^M‖₁ (the cross-class mass "
              "of every optimal plan here).", "",
              "| seed | m | ½‖w^m − w^M‖₁ | D fixed: observed / Θ*(w) | D adaptive η = σ₀, ∞: observed / ϑ^e(w) | D adaptive (∞) / fixed | W² fixed: cross + within | W² adaptive η = σ₀: cross + within | W² adaptive η = ∞: cross + within |",
              "|---|---|---|---|---|---|---|---|---|"]
    stats = dict(rf=[], re=[], ratio=[], cross_f=[], cross_a=[], cross_a_out=[], w2_ratio=[], within_ratio=[])
    for seed in range(N_SEEDS):
        for m in MS:
            o0, o1, o2 = po['B'][(0.0, seed, m)], po['B'][(SIGMA0, seed, m)], po['B'][(np.inf, seed, m)]
            ex_ok = np.isfinite(o0['D_exact_star'])
            stats['rf'].append(o0['D'] / o0['D_fixed_star'])
            if ex_ok:
                stats['re'] += [o1['D'] / o0['D_exact_star'], o2['D'] / o0['D_exact_star']]
                stats['w2_ratio'].append(o2['W2']['W2'] / o0['W2']['W2']); stats['within_ratio'].append(o2['W2']['within'] / max(o0['W2']['within'], 1e-300))
            stats['ratio'].append(o2['D'] / o0['D'])
            stats['cross_f'].append(o0['W2']['cross'] / o0['W2']['W2'])
            stats['cross_a' if ex_ok else 'cross_a_out'] += [o1['W2']['cross'] / o1['W2']['W2'], o2['W2']['cross'] / o2['W2']['W2']]
            dex = f"{o0['D_exact_star']:.3f}" if ex_ok else "—"
            w2 = lambda o: f"{o['W2']['cross']:.3f} + {o['W2']['within']:.3f}"
            lines.append(f"| {seed} | {m} | {o0['W2']['half_l1']:.4f} | {o0['D']:.3f} / {o0['D_fixed_star']:.3f} | {o1['D']:.3f}, {o2['D']:.3f} / {dex} | {o2['D'] / o0['D']:.2f} | "
                         f"{w2(o0)} | {w2(o1)} | {w2(o2)} |")
    lines += ["", f"The fixed prior's D is {_rng(stats['rf'], '{:.3f}')} of the replicated Θ*(w)'s; the adaptive prior's is {_rng(stats['re'], '{:.3f}')} of the exact representations' "
              f"(both weight vectors in W); the adaptive/fixed ratio of D is {_rng(stats['ratio'], '{:.2f}')}. The cross-class transport is {_rng(stats['cross_f'], '{:.2f}')} of W² "
              f"under the fixed prior and {_rng(stats['cross_a'], '{:.2f}')} under the adaptive ones where both weight vectors are in W ({_rng(stats['cross_a_out'], '{:.2f}')} otherwise, the parameters "
              f"running off); there the adaptive (η = ∞)/fixed ratio of W² is "
              f"{_rng(stats['w2_ratio'], '{:.2f}')}, of its within-class part {_rng(stats['within_ratio'], '{:.2f}')}."]
    return lines


def report_C_checks(C, P):
    po = P['post']
    designs = (('nearest', 0.0, False), ('nearest', SIGMA0, False), ('ot', 0.0, False), ('ot', SIGMA0, 'settled'))
    name = lambda plan, eta, gate: f"{'OT' if plan == 'ot' else 'nearest atom'}, {'adaptive' if eta > 0 else 'fixed'}"
    lines = ["", "The coupling level and the measure level at the end (means over the seeds, min–max): the coupled distance D of the shared particles and W²(μ_m, μ_M), "
             "with the cross-class part's share (classes: the cells of the exact atoms); the particles whose end cell differs from their initial class.", "",
             "| leak, prior | m | D | W² | cross-class share of W² (—: W² below 10⁻¹²) | relabelled particles (all seeds) |", "|---|---|---|---|---|---|"]
    for d in designs:
        for m in MS:
            o = [po['C'][d + (s, m)] for s in range(N_SEEDS)]
            D_ = np.array([x['D'] for x in o]); W_ = np.array([x['W2']['W2'] for x in o])
            cs = np.array([x['W2']['cross'] / x['W2']['W2'] if x['W2']['W2'] > 1e-12 else np.nan for x in o])
            lines.append(f"| {name(*d)} | {m} | {D_.mean():.3f} ({_rng(D_, '{:.3f}')}) | {W_.mean():.2e} ({_rng(W_, '{:.1e}')}) | {_rng(cs, '{:.2f}')} | {sum(x['relabelled'] for x in o)} |")
    lines += ["", "The fixed misplaced prior's Θ*(w) at the sampled weights (Part C's nearest-atom fixed design): the D‴.5 margin (its class) and the end diameter "
              "(sub-clusters).", "", "| m | " + " | ".join(f"seed {s}" for s in range(N_SEEDS)) + " |", "|---|" + "---|" * N_SEEDS]
    for m in MS:
        cells = []
        for s in range(N_SEEDS):
            f = po['fpC'][(s, m)]; k = int(np.argmin(f['sep'])); o = po['C'][('nearest', 0.0, False, s, m)]
            cells.append(f"{f['sep'][k]:+.2e} (class {k + 1}){'' if f['stable'] else ' saddle'}; {o['spread']:.2f} {o['sizes']}")
        lines.append(f"| {m} | " + " | ".join(cells) + " |")
    to = [po['C'][('ot', SIGMA0, 'settled', s, m)]['t_open'] for s in range(N_SEEDS) for m in MS]
    to = [float(x) for x in to if x is not None and np.isfinite(float(x))]
    lines += ["", f"The ranking depends on the horizon: the OT leak toward the adaptive prior equals the fixed one until its settled gate opens (s = {_rng(to, '{:.0f}')}). "
              "The fit errors at s = 1250 and 2000 (means over the seeds, min–max):", "", "| leak, prior | m | fit error at s = 1250 | at s = 2000 |", "|---|---|---|---|"]
    for d in designs:
        for m in MS:
            o = [po['C'][d + (s, m)]['fit_at'] for s in range(N_SEEDS)]
            a, b = np.array([x[1250.0] for x in o]), np.array([x[2000.0] for x in o])
            lines.append(f"| {name(*d)} | {m} | {a.mean():.2e} ({_rng(a, '{:.1e}')}) | {b.mean():.2e} ({_rng(b, '{:.1e}')}) |")
    return lines


def report_D_checks(P):
    po = P['post']
    out = []
    for m in (15, 60):
        f = po['fp'][(0, m)]; k = int(np.argmin(f['sep']))
        out.append(f"m = {m}: D‴.5 margin {f['sep'][k]:+.2e} (class {k + 1}){', a saddle: the fixed-prior run sits on a saddle orbit (Proposition D‴.15), held there by the exact replication of its start' if not f['stable'] else ', stable'}")
    return ["", "The fixed prior's Θ*(w) at seed 0's sampled weights: " + "; ".join(out) + "."]


def report_E(A, P):
    en, mc = P['enum'], P['mc']
    law = {r['frac']: r['numeric'] for r in A['sig']}[0.15]
    lines = ["", "## E. The instance's probabilities: representability, the fixed prior's saddles, and the adaptive prior's expected limiting floor", "",
             "Exact sums of the multinomial probabilities (q = 1/3) over the count vectors of probability above " + ", ".join(f"{P_MIN[m]:g} (m = {m})" for m in MS)
             + "; the neglected mass is given. Representable: w ∈ W (an exact representation on the branch of ϑ^e(q)); boundary: some w_k = |A_k|/a_max; with a permutation: "
             "outside W, but a permutation of the classes over the bumps is representable (an exact representation exists off the branch). Hoeffding: "
             "Σ_kexp(−2m(q_k − |A_k|/a_max)²) bounds P(outside W ∪ boundary). The fixed matching prior's saddle: Θ*(w) not D‴.5-admissible (classes with at least "
             "two particles; a vector with an empty class has no K₀-cluster configuration). The adaptive prior's expected limiting floor: Σ p(w)·min_{|a| ≤ a_max} "
             "⟨(f − f̂_w)²⟩_ρ over the vectors outside W (the constrained minimum on the branch of ϑ^e(q), by L-BFGS-B from the bumps, taken as the dynamics' limit "
             "— the dynamics from the warm clouds follows the branch; a local minimum, not the infimum, which is 0 where a permutation is representable; 0 on W and on "
             "its boundary, the reference taken exact), a lower bound on that limit since the sum is truncated; against the fixed prior's law (1/m − 1/M)·tr(FC̃Fᵀ)/N and its Monte Carlo "
             "expectation (Part F).", "",
             "| m | vectors | neglected mass | P(outside W) | P(boundary) | Hoeffding | share of P(outside W) with a permutation | P(fixed-prior saddle) | P(empty class) | Newton failed | expected limiting floor, adaptive | fixed: law | fixed: Monte Carlo | adaptive/law | adaptive/Monte Carlo |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for m in MS:
        rows = en[m]
        S = lambda f: float(sum(r['p'] for r in rows if f(r)))
        tot = S(lambda r: True)
        pin = S(lambda r: r['rep'] < 0); pb = S(lambda r: r['rep'] == 0)
        hoef = float(np.sum(np.exp(-2 * m * (Q - THR / 40.0) ** 2)))
        perm = S(lambda r: r['rep'] < 0 and r['perm']) / pin if pin > 0 else float('nan')
        sad = S(lambda r: r['stable'] is False); emp = S(lambda r: r['empty']); fail = S(lambda r: (not r['empty']) and not (r['res'] < 1e-9))
        Ef = float(sum(r['p'] * r['floor'] for r in rows))
        lawm = (1.0 / m - 1.0 / M_REF) * law
        h2 = np.array([s[m]['h2'] for s in mc if 'h2' in s[m]])
        lines.append(f"| {m} | {len(rows)} | {'none (all vectors)' if P_MIN[m] == 0 else f'{max(1 - tot, 0.0):.1e}'} | {pin:.4f} | {pb:.4f} | {hoef:.3f} | {perm:.2f} | {sad:.4f} | {emp:.1e} | {fail:.1e} | {Ef:.2e} | {lawm:.2e} | "
                     f"{h2.mean():.2e} | {Ef / lawm:.3f} | {Ef / h2.mean():.3f} |")
    return lines


def report_F(A, P):
    mc = P['mc']
    sig = {r['frac']: r for r in A['sig']}[0.15]
    lines = ["", f"## F. Monte Carlo: the exact floors and coupled distances against the laws ({len(mc)} nested samples)", "",
             "M = 510 labels i.i.d. uniform (seeds 10⁴ + i), the m-system the first m. At each sample the replicated stationary points: the fixed prior's Θ*(w) and "
             "the adaptive prior's ϑ^e(w). The ratios to the laws per unit of 1/m − 1/M (Part A, σ₀ = 0.15σ_crit), mean ± standard error: ‖f̂_m − f̂_M‖² of the fixed prior "
             "(all samples, and conditionally on both samples D‴.5-admissible), and the coupled distances (fixed: all samples; adaptive: both samples in W).", "",
             "| m | E‖f̂_m − f̂_M‖²/law, fixed | conditional on D‴.5-admissible | samples | E D/law, fixed | E D/law, adaptive (in W) | samples | E D adaptive / E D fixed (same samples) |",
             "|---|---|---|---|---|---|---|---|"]
    msr = lambda x: f"{np.mean(x):.3f} ± {np.std(x, ddof=1) / np.sqrt(len(x)):.3f}"
    for m in MS:
        c = 1.0 / m - 1.0 / M_REF
        h2 = np.array([s[m]['h2'] for s in mc if 'h2' in s[m]]) / (c * sig['numeric'])
        h2c = np.array([s[m]['h2'] for s in mc if 'h2' in s[m] and s[m]['stable']]) / (c * sig['numeric'])
        Df = np.array([s[m]['D_f'] for s in mc if 'D_f' in s[m]]) / (c * sig['D_fixed'])
        both = [s[m] for s in mc if 'D_e' in s[m] and 'D_f' in s[m]]
        De = np.array([s['D_e'] for s in both]) / (c * sig['D_gauge'])
        rat = np.mean([s['D_e'] for s in both]) / np.mean([s['D_f'] for s in both])
        lines.append(f"| {m} | {msr(h2)} | {msr(h2c)} | {len(h2)}, {len(h2c)} | {msr(Df)} | {msr(De)} | {len(De)} | {rat:.2f} |")
    tail = []
    for m in MS:
        c = 1.0 / m - 1.0 / M_REF
        Df = np.array([s[m]['D_f'] for s in mc if 'D_f' in s[m]]) / (c * sig['D_fixed'])
        De = np.sort(np.array([s[m]['D_e'] for s in mc if 'D_e' in s[m] and 'D_f' in s[m]]) / (c * sig['D_gauge']))
        tail.append((m, Df.max(), De[-1], np.mean(De[:-1])))
    lines += ["", "The tails (the largest sample's D/law; in parentheses the adaptive prior's mean without that sample): "
              + "; ".join(f"m = {m}: fixed {a:.2f}, adaptive {b:.1f} ({c:.3f})" for m, a, b, c in tail)
              + ". The adaptive prior's gauge diverges at the boundary of W, so its coupled distance has a far heavier tail, and its mean over these samples is "
              "not a reliable estimate of its law."]
    return lines


# ---------------------------------------------------------------------------- report
def report(A, B, C, D=None, t0=None, P=None):
    lines = ["# Proposition D‴.24: the finite-size floor is the leak's bias; an adaptive prior removes it without a mass constraint (`exp_adaptive_floor.py`)", "",
             f"The three-bump instance (σ_crit = {SIGMA_CRIT:.5f}); the warm clouds of Proposition D‴.12 (width {JITTER:g}), m = {', '.join(map(str, MS))} against M = {M_REF}, {N_SEEDS} seeds; "
             f"midpoint rule dt = {DT:g}, the plan every {PLAN_EVERY} steps.", "",
             "## A. The floor as the weight-derivative of the leak's bias (matching prior, no dynamics)", "",
             "F: Proposition D‴.9's function-level sensitivity along e₁ − e₃, e₂ − e₃ (central differences of the fixed-prior nearest-atom stationary point); F_bias = Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰. "
             "The law per unit of 1/m − 1/M, tr(FC̃Fᵀ)/N, from F and as the sum of the mode terms h_jσ₀²/(h_j + σ₀)²z_jᵀC̃z_j; the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ).", "",
             "| σ₀/σ_crit | ‖F − F_bias‖/‖F‖ | law from F | sum of the mode terms | bound | law/bound | floor (1/m − 1/M)·law at m = 15, 60, 240 |",
             "|---|---|---|---|---|---|---|"]
    for r in A['sig']:
        fl = ", ".join(f"{(1.0 / m - 1.0 / M_REF) * r['numeric']:.2e}" for m in MS)
        lines.append(f"| {r['frac']:g} | {r['rel']:.1e} | {r['numeric']:.4e} | {r['law']:.4e} | {r['bound']:.4e} | {r['law'] / r['bound']:.3f} | {fl} |")
    lines += ["", "The parameter level: the fixed prior's weight sensitivity S (whitened weighted) against H_L(H_L + σ₀I)⁻¹S^e with S^e = ∂_wϑ^e (per mode the fraction "
              "h_j/(h_j + σ₀) of the gauge's displacement, the rest σ₀/(h_j + σ₀) being the bias), and Proposition D‴.5's coupled-distance law per unit of 1/m − 1/M, "
              "Σ tr(S̃C̃S̃ᵀ), for the fixed prior and for the adaptive one (S^e: the whole displacement in the parameters).", "",
              "| σ₀/σ_crit | ‖S − H_L(H_L + σ₀I)⁻¹S^e‖/‖S‖ | coupled-distance law, fixed prior | adaptive prior | ratio |", "|---|---|---|---|---|"]
    for r in A['sig']:
        lines.append(f"| {r['frac']:g} | {r['split_rel']:.1e} | {r['D_fixed']:.4e} | {r['D_gauge']:.4e} | {r['D_gauge'] / r['D_fixed']:.3f} |")
    lines += ["", "Mode by mode (whitened weighted Gauss–Newton eigenvalues h_j; the share of each term in the law):", "",
              "| σ₀/σ_crit | " + " | ".join(f"mode {j + 1}" for j in range(9)) + " |", "|---|" + "---|" * 9]
    for r in A['sig']:
        lines.append(f"| {r['frac']:g}: h_j/σ₀ | " + " | ".join(f"{x:.2f}" for x in r['h'] / r['sigma0']) + " |")
        lines.append(f"| {r['frac']:g}: share | " + " | ".join(f"{x:.3f}" for x in r['terms'] / r['terms'].sum()) + " |")
    if P:
        lines += report_A_checks(A, P)
    lines += ["", f"Representable weights (the set W): an exact representation at the weights w needs |A_k|/w_k < a_max = 4, i.e. w > ({', '.join(f'{x:.3f}' for x in np.abs(A_BUMP) / 4.0)}). "
              "Per seed and size: the sampled counts, the margin min_k(w_k − |A_k|/a_max), and the exact representation's gradient residual and fit error.", "",
              "| m | " + " | ".join(f"seed {s}" for s in range(N_SEEDS)) + " |", "|---|" + "---|" * N_SEEDS]
    for m in MS + (M_REF,):
        cells = []
        for s in range(N_SEEDS):
            a = A['adm'][(s, m)]
            cells.append(f"{tuple(int(x) for x in a['counts'])}: {a['margin']:+.4f}" + (f", {a['res']:.0e}, {a['fit']:.0e}" if np.isfinite(a['res']) else " (none)"))
        lines.append(f"| {m} | " + " | ".join(cells) + " |")
    # B
    lines += ["", f"## B. The nearest-atom leak at the matching prior: fixed (η = 0) against adaptive priors ({S_B:g} slow units)", "",
              "‖f̂_m − f̂_M‖² (the function-level floor) and the fit error at s = 600, 3000, 6000, per seed; W²(μ_m, μ_M) at the end against the W² between the exact "
              "representations of the two systems' sampled weights (the measure-level distance; its split, and the coupled distance D, below); the class means' whitened "
              "weighted distance from the exact representation of the sampled weights at the end; the counts (all conserved). —: no exact representation at those weights "
              "(outside W or on its boundary).", "",
              "| η | m | ‖f̂_m − f̂_M‖² at s = 600, 3000, 6000 (per seed) | fit error at 6000 (per seed) | W²(μ_m, μ_M) at the end / exact-representation W² (per seed) | distance from ϑ^e(w) at the end | counts conserved |",
              "|---|---|---|---|---|---|---|"]
    for eta in (0.0, SIGMA0, np.inf):
        rs = sorted([r for r in B if r['eta'] == eta], key=lambda r: r['seed'])
        for m in MS:
            h2 = "; ".join(", ".join(f"{r['h2'][m][int(np.argmin(np.abs(r['s'] - t)))]:.1e}" for t in (600.0, 3000.0, S_B)) for r in rs)
            fit = "; ".join(f"{r['fit'][m][-1]:.1e}" for r in rs)
            w2 = "; ".join(f"{r['W2_end'][m]:.2e} / " + (f"{r['W2_exact'][m]:.2e}" if np.isfinite(r['W2_exact'][m]) else "—") for r in rs)
            dev = "; ".join((f"{r['dev_end'][m]:.1e}" if np.isfinite(r['dev_end'][m]) else "—") if 'dev_end' in r else f"{r['dev'][m][-1]:.1e}" for r in rs)
            cons = all(np.all(r['counts'][m] == r['counts'][m][0]) for r in rs)
            lines.append(f"| {'∞' if np.isinf(eta) else f'{eta:.5f}'} | {m} | {h2} | {fit} | {w2} | {dev} | {'yes' if cons else 'no'} |")
    if P:
        lines += report_B_checks(B, P)
    # C
    lines += ["", f"## C. The design table at the misplaced prior (centres +{SHIFT_C:g}), {S_C:g} slow units", "",
              "The fit error (the bias; 0 at an exact representation) and ‖f̂_m − f̂_M‖² (the floor) at the end, means over the seeds (min–max); the OT leak toward the adaptive "
              "prior with the settled gate of Theorem D‴.23.", "",
              "| leak | prior | m | fit error at the end | ‖f̂_m − f̂_M‖² at the end |", "|---|---|---|---|---|"]
    for plan, eta, gate in (('nearest', 0.0, False), ('nearest', SIGMA0, False), ('ot', 0.0, False), ('ot', SIGMA0, 'settled')):
        rs = [r for r in C if r['plan'] == plan and r['eta'] == eta and r['gate'] == gate]
        for m in MS:
            fit = np.array([r['rec'][m]['fit'][-1] for r in rs]); h2 = np.array([r['h2_end'][m] for r in rs])
            lines.append(f"| {'OT' if plan == 'ot' else 'nearest atom'} | {'adaptive (η = σ₀' + (', settled gate)' if gate else ')') if eta > 0 else 'fixed'} | {m} | "
                         f"{fit.mean():.2e} ({fit.min():.1e}–{fit.max():.1e}) | {h2.mean():.2e} ({h2.min():.1e}–{h2.max():.1e}) |")
    if P:
        lines += report_C_checks(C, P)
    if D:
        win = lambda ts, e2, t: float(np.sqrt(np.mean(e2[(ts > t - 5 * PER) & (ts <= t)])))
        times = (5 * PER, 1000.0, 2000.0, T_REAL)
        lines += ["", f"## D. The real closed loop (k_e = {KE:g}, dt = {DT_REAL:g}, {T_REAL:g} time units), matching prior, nearest-atom leak, seed 0's sampled weights", "",
                  "The rms tracking error over the five periods ending at t = " + ", ".join(f"{x:.0f}" for x in times) + ", against (1/k_e)‖f − f̂‖_ρ of the averaged loop "
                  "from the same start at s = t/k_e (the slaved error of Theorem J at leading order); from the exact representation of w with the adaptive prior at its atoms "
                  "the tracking error stays at the integrator's level.", "",
                  "| m | counts | start, prior | e_rms at those times | (1/k_e)‖f − f̂‖_ρ averaged at s = t/k_e | ratio at the end |", "|---|---|---|---|---|---|"]
        name = {'fixed': 'Θ*(w), fixed prior', 'adaptive': 'Θ*(w), adaptive prior (η = ∞)', 'exact': 'ϑ^e(w), adaptive prior at ϑ^e(w)'}
        for d in sorted(D, key=lambda d: (d['m'], ['fixed', 'adaptive', 'exact'].index(d['case']))):
            er = [win(d['ts'], d['e2'], t) for t in times]
            pr = [float(np.interp(t / KE, d['s_avg'], d['fit_avg'])) * d['normF'] / KE for t in times]
            lines.append(f"| {d['m']} | {tuple(int(x) for x in d['counts'])} | {name[d['case']]} | " + ", ".join(f"{x:.3e}" for x in er) + " | "
                         + ", ".join(f"{x:.3e}" for x in pr) + f" | {er[-1] / pr[-1]:.3f} |" if pr[-1] > 1e-12 else
                         f"| {d['m']} | {tuple(int(x) for x in d['counts'])} | {name[d['case']]} | " + ", ".join(f"{x:.1e}" for x in er) + " | " + ", ".join(f"{x:.1e}" for x in pr) + " | — |")
    if P and D:
        lines += report_D_checks(P)
    if P:
        lines += report_E(A, P) + report_F(A, P)
    lines += ["", "Numbers are produced by `experiments/exp_adaptive_floor.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'adaptive_floor.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def save_raw(raw, **arrays):
    """np.savez to a temporary file, then an atomic rename: a crash while pickling cannot truncate the stored results."""
    tmp = raw[:-4] + '.tmp.npz'
    np.savez(tmp, **arrays)
    os.replace(tmp, raw)


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'adaptive_floor_raw.npz')
    if '--from-raw' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        report(z['A'].item(), list(z['B']), list(z['C']), list(z['D']) if 'D' in z.files else None, t0, z['P'].item() if 'P' in z.files else None)
        return
    if '--only' in sys.argv and sys.argv[sys.argv.index('--only') + 1] == 'P':          # the checks at the sampled weights, from the stored finals
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        with mp.Pool(N_WORKERS) as pool:
            P = checks(list(save['B']), list(save['C']), list(save['D']) if 'D' in save else None, pool)
        save['P'] = np.array(P, dtype=object)
        save_raw(raw, **save)
        report(save['A'].item(), list(save['B']), list(save['C']), list(save['D']) if 'D' in save else None, t0, P)
        return
    if '--only' in sys.argv and sys.argv[sys.argv.index('--only') + 1] == 'A':          # Part A and Part B's post-processing
        with np.load(raw, allow_pickle=True) as z:                                   # read everything before the file is rewritten
            save = {k: z[k] for k in z.files}
        A = part_A(); B = post_B(list(save['B']))
        save['A'] = np.array(A, dtype=object); save['B'] = np.array(B, dtype=object)
        save_raw(raw, **save)
        report(A, B, list(save['C']), list(save['D']) if 'D' in save else None, t0, save['P'].item() if 'P' in save else None)
        return
    if '--only' in sys.argv and sys.argv[sys.argv.index('--only') + 1] == 'D':
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        with mp.Pool(N_WORKERS) as pool:
            D = pool.map(part_D_job, [(m, case) for m in (15, 60) for case in ('fixed', 'adaptive', 'exact')], chunksize=1)
        save['D'] = np.array(D, dtype=object)
        save_raw(raw, **save)
        report(save['A'].item(), list(save['B']), list(save['C']), D, t0, save['P'].item() if 'P' in save else None)
        return
    A = part_A()
    print("[A done]", flush=True)
    jobsB = [('B', seed, 'nearest', eta, False, 0.0, S_B) for eta in (0.0, SIGMA0, np.inf) for seed in range(N_SEEDS)]
    jobsC = [('C', seed, plan, eta, gate, SHIFT_C, S_C) for plan, eta, gate in (('nearest', 0.0, False), ('nearest', SIGMA0, False), ('ot', 0.0, False), ('ot', SIGMA0, 'settled'))
             for seed in range(N_SEEDS)]
    with mp.Pool(N_WORKERS) as pool:
        B = post_B(pool.map(run_job_h2, jobsB, chunksize=1))
        C = pool.map(run_job, jobsC, chunksize=1)
        P = checks(B, C, None, pool)
    save_raw(raw, A=np.array(A, dtype=object), B=np.array(B, dtype=object), C=np.array(C, dtype=object), P=np.array(P, dtype=object))
    report(A, B, C, None, t0, P)


if __name__ == '__main__':
    main()
