"""
Theorem D‴.23: adaptive priors — a learned prior removes the leak's bias, at the rate of the measure-space excitation
==================================================================================================================

The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit with the OT leak toward a prior whose centres are shifted by +0.15
(a misplaced prior: the fixed-prior stationary point Θ*(q) is biased), masses q = 1/3, and the adaptive prior of src/core/adaptive_prior.py.

A. The linear theory at the exact representation: the whitened weighted Gauss–Newton spectrum h_j of the loss, the fraction σ₀/(h_j + σ₀)
   of the prior's error that a fixed prior keeps in mode j, and the rates λ_j^− of the collective pairs for η = σ₀/10, σ₀, 10σ₀, ∞; the
   radius bookkeeping of (v): a sampled lower bound on the Lipschitz constant of the collective Jacobian, the resulting upper bound on the
   proof's radius for each η, and the norms of the starts of Part C.
B. The fixed prior's bias: Θ*(q) by Newton against the first-order prediction y* = σ₀(H_L + σ₀I)⁻¹x⁰, mode by mode, and against the shift.
C. The K₀-atom adaptive loop (one particle per class; every replicated lift carries it) from Θ*(q) with Z = the shifted prior, for the four
   η, 6000 slow units: the mode amplitudes of the atoms' deviation from the exact representation against the exact solutions of the pairs,
   the fit error, the energy.
D. m-particle systems: m = 15, 60, 240 and M = 510 from the warm clouds of Proposition D‴.12 (four seeds), η = 0 (the fixed prior), σ₀ and ∞,
   3000 slow units, ungated and with two gates — 'transit' (the prior frozen while some particle is assigned to an atom that is not its
   nearest) and 'settled' (in addition, the first time, until the fixed-prior drift is below 10⁻⁵): the fit error, W²(μ_Θ, μ_exact), the
   class means' whitened weighted deviation from the exact representation against the K₀-atom run of Part C (per seed at s − t_open for
   the settled gate), the prior's distance, the final counts. `--only-gated transit|settled` recomputes one gate and merges it.
E. The ungated η = ∞ system with the largest ‖ȳ‖ at s = 3000, re-run from its seed and continued to 23000 slow units: the fit error, the
   energy, the class-centre spread and the drift; at the end the unregularised loss's Hessian at the collapsed atoms, a BFGS descent of
   that loss from them, and the collapsed η = ∞ dynamics from them — the unregularised K₀-atom gradient flow — with its crossing times
   (`--only-E` recomputes Part E and merges it).
F. The test of (v): from the lift of Θ*(q) plus zero-mean per-particle separations (rms 0.02, 0.1; m = 15, 510; η = σ₀, ∞), ungated, 1750
   slow units, ‖ȳ‖ against the K₀-atom run of Part C (`--only-F`; `--only-lin` recomputes Part A–B's linear parts with the shift table).
G. The limits of (i) (`--only-G`): a non-replicated stationary point held by the leak (a four-bump target, m = 6, one class split), and the
   joint flow from the exact split fit; unbounded solutions of the adaptive prior (a target bump beyond the saturation, K₀-atom system,
   LSODA to s = 10⁶) against the fixed prior's bounded one; the exact representation with the prior's atoms permuted.
H. Proposition D‴.27 (`--only-H`): the projected adaptive prior (the atoms confined to a box in the unconstrained coordinates, amplitudes and
   centres to 90% of their saturations) on Part G(b)'s target (bounded, convergent, the a-priori constant B/σ₀; the projected midpoint rule
   at two steps), on Part C's K₀-atom loop at the misplaced prior (the box never approached: the run identical to Part C's by
   construction), with the box's amplitude bound varied on Part G(b)'s target, and on Part G(a)'s four-bump target (the box stops the
   width's run-off).

Output: results/chaos/adaptive_prior.md and adaptive_prior_raw.npz (AP_SEEDS, AP_WORKERS; `--from-raw` regenerates the md).
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
from src.core import adaptive_prior as AP

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
SHIFT = 0.15
Q = np.ones(3) / 3
ETAS = (SIGMA0 / 10, SIGMA0, 10 * SIGMA0, np.inf)
ETAS_D = (0.0, SIGMA0, np.inf)
M_REF, MS = 510, (15, 60, 240)
JITTER = 0.02
S_C, S_D, DT, PLAN_EVERY = 6000.0, 3000.0, 0.25, 4
REC_C, REC_D = 100, 100                                                   # records every 25 slow units
FLOOR = 1e-6                                                              # mode amplitudes below it are at the second-order remainder's level
N_SEEDS = int(os.environ.get('AP_SEEDS', 4))
N_WORKERS = int(os.environ.get('AP_WORKERS', 3))


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def representation(net, shift=0.0):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def setup():
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    exact = representation(net, 0.0)
    prior = representation(net, SHIFT)
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, prior, n_quad=N_Q)
    return net, exact, prior, loop


def class_means(Th, gamma, q):
    return AP.barycentres(Th, gamma, q)


def linear_parts():
    net, exact, prior, loop = setup()
    HL, h, V = AP.gauss_newton_block(loop, exact, Q)
    star, res, marg = loop.stationary_atoms(prior)                          # the fixed prior's stationary point
    y0 = AP.whitened(prior - exact, Q, G)
    ya = AP.whitened(star - exact, Q, G)
    yp = SIGMA0 * np.linalg.solve(HL + SIGMA0 * np.eye(HL.shape[0]), y0)
    fit_star = loop.fit_error(np.repeat(star, 1, axis=0))
    rates = {eta: AP.pair_rates(h, SIGMA0, eta)[0] for eta in ETAS}
    conv = []                                                               # the first-order formula against the shift
    for sh in (0.15, 0.075, 0.0375, 0.01875):
        pr = representation(net, sh)
        lp = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, pr, n_quad=N_Q)
        st, rs, _ = lp.stationary_atoms(pr)
        x0s = AP.whitened(pr - exact, Q, G); ys = AP.whitened(st - exact, Q, G)
        yps = SIGMA0 * np.linalg.solve(HL + SIGMA0 * np.eye(HL.shape[0]), x0s)
        conv.append(dict(shift=sh, res=rs, rel=float(np.linalg.norm(ys - yps) / np.linalg.norm(ys)), mode_ratio=(V.T @ ys) / (V.T @ yps)))
    return dict(h=h, V=V, HL=HL, star=star, newton_res=res, diag=marg, y0=y0, ya=ya, yp=yp, fit_star=loop.fit_error(star), rates=rates, exact=exact, prior=prior,
                conv=conv, radius=radius_bookkeeping(loop, exact, star, h, HL, y0, ya))


def collective_drift(loop, atoms):
    """The unregularised K₀-atom drift dϑ_k/ds = −G∇V_L[Σ_lq_lδ_{ϑ_l}](ϑ_k) (the collective part of the joint drift without the leak, whose
    part is linear)."""
    a_, c_, s_, d_, phi = loop.features(atoms)
    return -G[None, :] * loop.grad_V_L(loop.F - phi @ (Q * a_), atoms)


def collective_jacobian(loop, atoms, h=1e-6):
    """Its Jacobian in the whitened weighted coordinates y_k = √q_kG^{−1/2}(ϑ_k − ϑ*_k), by central differences (= −H_L at an exact
    representation)."""
    S = np.repeat(np.sqrt(Q), 3) / np.tile(np.sqrt(G), 3)
    d = atoms.size; J = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        J[:, j] = (collective_drift(loop, (atoms.ravel() + e).reshape(atoms.shape)) - collective_drift(loop, (atoms.ravel() - e).reshape(atoms.shape))).ravel() / (2 * h)
    return S[:, None] * J / S[None, :]


def radius_bookkeeping(loop, exact, star, h, HL, y0, ya):
    """(v): a sampled lower bound C on the Lipschitz constant C_N of the collective Jacobian in the whitened weighted norm (at Θ*(q) and at 160
    random points at whitened distances 0.005–0.05 from the exact representation); for each η the upper bound min{σ₀, λ_1^−(η)}/(4C) on the
    proof's radius (an upper bound: C ≤ C_N, and the plan-margin term only lowers it); the norms |||·|||_η of Part C's start (Θ*(q), the shifted
    prior: (‖y*‖² + (σ₀/η)‖x⁰‖²)^{1/2}; at η = ∞ the prior is the class barycentres and the norm has no ζ part)."""
    J0 = collective_jacobian(loop, exact)
    rng = np.random.default_rng(3)
    pts = [star]
    for r in (0.005, 0.01, 0.02, 0.05):
        for _ in range(40):
            y = rng.standard_normal(9); y *= r / np.linalg.norm(y)
            pts.append(exact + y.reshape(3, 3) * np.sqrt(G)[None, :] / np.sqrt(Q)[:, None])
    C = max(np.linalg.norm(collective_jacobian(loop, P) - J0, 2) / np.linalg.norm(AP.whitened(P - exact, Q, G)) for P in pts)
    rows = []
    prior = representation(loop.net, SHIFT)
    clouds = [CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G) for seed in range(N_SEEDS)]
    for eta in ETAS:
        lam1 = AP.pair_rates(h, SIGMA0, eta)[0][0]; lam0 = min(SIGMA0, lam1)
        start = float(np.linalg.norm(ya)) if np.isinf(eta) else float(np.sqrt(np.sum(ya ** 2) + (SIGMA0 / eta) * np.sum(y0 ** 2)))
        d_starts = [start_norm(cl[:m], prior, exact, eta) for cl in clouds for m in list(MS) + [M_REF]]   # Part D's warm clouds
        rows.append(dict(eta=eta, lam1=lam1, lam0=lam0, rho=lam0 / (4 * C), start=start, d_min=min(d_starts), d_max=max(d_starts)))
    return dict(C=C, n_pts=len(pts), J0_err=float(np.linalg.norm(J0 + HL, 2) / np.linalg.norm(HL, 2)), rows=rows)


def k0_job(eta):
    """Part C: the K₀-atom adaptive loop from Θ*(q) with Z = the shifted prior."""
    net, exact, prior, loop = setup()
    lin = linear_parts()
    apl = AP.AdaptivePriorLoop(loop, eta)
    rec = dict(s=[], y=[], z=[], fit=[], J=[])

    def cb(s, Th, Z, adapting):
        rec['s'].append(s); rec['y'].append(AP.whitened(Th - exact, Q, G)); rec['z'].append(AP.whitened(Z - exact, Q, G))
        rec['fit'].append(loop.fit_error(Th)); rec['J'].append(apl.energy(Th, Z))
    apl.run(lin['star'].copy(), prior.copy(), S_C, DT, plan_every=1, rec_every=REC_C, callback=cb)
    return dict(eta=eta, **{k: np.array(v) for k, v in rec.items()})


def m_job(args):
    """Part D: the m-particle systems from the warm clouds."""
    seed, eta, gate = args
    t0 = time.time()
    net, exact, prior, loop = setup()
    rng = np.random.default_rng(seed)
    clouds = CL.mixture_start(exact, M_REF, rng, JITTER, G)
    out = dict(seed=seed, eta=eta, gate=gate, rec={}, t_open={})
    for m in list(MS) + [M_REF]:
        apl = AP.AdaptivePriorLoop(loop, eta if eta > 0 else 0.0)
        rec = dict(s=[], fit=[], W2ex=[], ybar=[], zdist=[], J=[], counts=[], adapting=[])

        def cb(s, Th, Z, adapting, rec=rec, apl=apl):
            gamma, W2 = AP.plan_to(Th, Z, Q, G)
            rec['s'].append(s); rec['fit'].append(loop.fit_error(Th)); rec['W2ex'].append(OT.W2_empirical(Th, exact, G))
            rec['ybar'].append(AP.whitened(class_means(Th, gamma, Q) - exact, Q, G)); rec['zdist'].append(float(np.linalg.norm(AP.whitened(Z - exact, Q, G))))
            rec['J'].append(loop.L_bar(Th) + 0.5 * SIGMA0 * W2); rec['counts'].append(np.bincount(np.argmin(loop.cost_matrix(Th), axis=1), minlength=3))
            rec['adapting'].append(adapting)
        _, _, t_open = apl.run(clouds[:m].copy(), prior.copy(), S_D, DT, plan_every=PLAN_EVERY, rec_every=REC_D, callback=cb, gate=gate)
        out['rec'][m] = {k: np.array(v) for k, v in rec.items()}
        out['t_open'][m] = t_open
    print(f"[{time.time() - t0:.0f}s] seed {seed} eta {eta:g} gate {gate}", flush=True)
    return out


S_E, REC_E = 23000.0, 4000                                                 # Part E: continued to 23000 slow units, records every 1000


def trapped_job(D):
    """Part E: the ungated η = ∞ system with the largest ‖ȳ‖ at s = 3000 in Part D, re-run from its seed (the run is deterministic: the fit
    error at s = 3000 is compared with the stored one) and continued to S_E; every 2000 slow units the fit error, the energy, the largest
    spread of the class centres and the largest drift; at the end the collapsed class atoms, the eigenvalues of the unregularised loss's
    Hessian there (whitened, weighted) and a BFGS minimisation of the unregularised weighted loss from them."""
    from scipy.optimize import minimize
    cand = [(np.linalg.norm(d['rec'][m]['ybar'][int(np.argmin(np.abs(d['rec'][m]['s'] - 3000.0)))]), d['seed'], m, d)
            for d in D if np.isinf(d['eta']) and d.get('gate', False) is False for m in d['rec']]
    yb, seed, m, d = max(cand, key=lambda t: t[0])
    net, exact, prior, loop = setup()
    clouds = CL.mixture_start(exact, M_REF, np.random.default_rng(seed), JITTER, G)
    apl = AP.AdaptivePriorLoop(loop, np.inf)
    rec = dict(s=[], fit=[], J=[], spread=[], drift=[])

    def cb(s, Th, Z, adapting):
        gamma, _ = AP.plan_to(Th, Z, Q, G)
        lab = np.argmax(gamma, axis=1); c = loop.features(Th)[1]
        rec['s'].append(s); rec['fit'].append(loop.fit_error(Th)); rec['J'].append(apl.energy(Th, Z))
        rec['spread'].append(max(float(np.ptp(c[lab == j])) for j in range(3))); rec['drift'].append(float(np.abs(apl.drift(Th, Z, gamma)[0]).max()))
    Th, Z, _ = apl.run(clouds[:m].copy(), prior.copy(), S_E, DT, plan_every=PLAN_EVERY, rec_every=REC_E, callback=cb, gate=False)
    gamma, _ = AP.plan_to(Th, Z, Q, G); lab = np.argmax(gamma, axis=1)
    atoms = np.array([Th[lab == j].mean(axis=0) for j in range(3)])

    def L(v):
        at = v.reshape(3, 3); a_, c_, s_, d_, phi = loop.features(at)
        return 0.5 * float(np.mean((loop.F - phi @ (Q * a_)) ** 2))
    res = minimize(L, atoms.ravel(), method='BFGS', options=dict(gtol=1e-13, maxiter=20000))
    ev = AP.gauss_newton_block(loop, atoms, Q)[1]
    # the collapsed η = ∞ dynamics from there: the unregularised K₀-atom gradient flow dϑ_k/ds = −G∇V_L[Σ_lq_lδ_{ϑ_l}](ϑ_k) (LSODA)
    from scipy.integrate import solve_ivp

    def rhs(s_, v):
        at = v.reshape(3, 3); a_, c_, s__, d_, phi = loop.features(at)
        return (-G[None, :] * loop.grad_V_L(loop.F - phi @ (Q * a_), at)).ravel()
    t_eval = np.concatenate([[0.0], np.minimum(np.geomspace(1e2, 2e5, 200), 2e5)])
    L0 = L(atoms.ravel())
    levels = (0.9 * L0, 0.5 * L0, 1e-12)
    events = []
    for lev in levels:                                                     # the crossings located by events, not read off the grid
        evf = (lambda t, v, lev=lev: L(v) - lev); evf.terminal = False; evf.direction = -1
        events.append(evf)
    sol = solve_ivp(rhs, (0.0, 2e5), atoms.ravel(), method='LSODA', rtol=1e-10, atol=1e-12, t_eval=t_eval, events=events)
    Ls = np.array([L(sol.y[:, k]) for k in range(sol.y.shape[1])])
    first = lambda te: float(te[0]) if len(te) else float('nan')
    flow = dict(t=sol.t, L=Ls, plateau_until=first(sol.t_events[0]), t_half=first(sol.t_events[1]), t_exact=first(sol.t_events[2]), L_end=float(Ls[-1]),
                atoms_end=tuple(np.asarray(x) for x in loop.features(sol.y[:, -1].reshape(3, 3))[:3]))
    stored_fit = d['rec'][m]['fit'][int(np.argmin(np.abs(d['rec'][m]['s'] - 3000.0)))]
    feat = lambda at: tuple(np.asarray(x) for x in loop.features(at)[:3])
    return dict(seed=seed, m=m, ybar3000=yb, stored_fit3000=stored_fit, rec={k: np.array(v) for k, v in rec.items()}, atoms=feat(atoms),
                L_atoms=L(atoms.ravel()), hess_ev=ev, bfgs_L=float(res.fun), bfgs_atoms=feat(res.x.reshape(3, 3)), exact=feat(exact), flow=flow)


S_F, F_RMS, F_MS, F_DRAWS = 1750.0, (0.02, 0.1), (15, 510), 5
F_JOBS = [(m, eta, r, dr) for eta in (SIGMA0, np.inf) for r in F_RMS for m in F_MS for dr in range(F_DRAWS)]


def triple_norm(u_bar, u_perp, zeta, q, eta):
    """|||(u, ζ)|||_η = max{(Σ_kq_k(‖ū_k‖² + (σ₀/η)‖ζ_k‖²))^{1/2}, max_i‖u_i^⊥‖} (G⁻¹ norms; η = ∞: no ζ part)."""
    n2 = lambda X: np.sum(np.asarray(X) ** 2 / G[None, :], axis=1)
    coll = np.sum(q * n2(u_bar)) + (0.0 if np.isinf(eta) else (SIGMA0 / eta) * np.sum(q * n2(zeta)))
    return float(max(np.sqrt(coll), np.sqrt(n2(u_perp).max()) if len(u_perp) else 0.0))


def start_norm(Th, Z, exact, eta):
    """The norm |||·|||_η of a start (Θ, Z) around the exact representation, classes by the plan to ν_Z."""
    gamma, _ = AP.plan_to(Th, Z, Q, G); lab = np.argmax(gamma, axis=1)
    means = np.array([Th[lab == k].mean(axis=0) for k in range(3)])
    return triple_norm(means - exact, Th - means[lab], Z - exact, Q, eta)


def sep_job(args):
    """Part F, the structure behind (v) beyond its ball: from the lift of Θ*(q) plus zero-mean per-particle separations of the given rms (G⁻¹
    norm; the draw `dr` fixes their direction, the same for both rms), the prior at the shifted centres, ungated, S_F slow units: ‖ȳ‖ and the
    fit error at the end, to be compared with the K₀-atom run of Part C at S_F; the start's norm |||·|||_η."""
    m, eta, rms, dr = args
    t0 = time.time()
    net, exact, prior, loop = setup()
    star = loop.stationary_atoms(prior)[0]
    lab = np.repeat(np.arange(3), m // 3)
    rng = np.random.default_rng(7 + dr)
    u = rng.standard_normal((m, 3)) * np.sqrt(G)[None, :]
    for k in range(3):
        u[lab == k] -= u[lab == k].mean(axis=0)
    u *= rms / np.sqrt(np.mean(np.sum(u ** 2 / G[None, :], axis=1)))
    Th0 = star[lab] + u
    apl = AP.AdaptivePriorLoop(loop, eta)
    Th, Z, _ = apl.run(Th0, prior.copy(), S_F, DT, plan_every=PLAN_EVERY, rec_every=10 ** 9, gate=False)
    gamma, _ = AP.plan_to(Th, Z, Q, G)
    yb = float(np.linalg.norm(AP.whitened(class_means(Th, gamma, Q) - exact, Q, G)))
    print(f"[{time.time() - t0:.0f}s] F m {m} eta {eta:g} rms {rms:g} draw {dr}", flush=True)
    return dict(m=m, eta=eta, rms=rms, draw=dr, ybar=yb, fit=loop.fit_error(Th), start=start_norm(Th0, prior, exact, eta))


# ---------------------------------------------------------------------------- G: the limits of (i)
B4 = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.4, 1.3, 0.4), (0.4, 2.3, 0.4)]    # four bumps: six kernels fit them exactly with one class split
B_SAT = [(1.6, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]                   # the first bump's weighted amplitude 3·1.6 = 4.8 > a_max = 4
S_G = 40000.0
T_CONT = (1e6, 1e7, 1e8, 1e9, 1e10)                                        # Part G(a): the collapsed flow continued
T_G = (0.0, 1e2, 1e3, 1e4, 1e5, 1e6)


def _bumps(B):
    return lambda X: sum(A * np.exp(-(X - C_) ** 2 / (2 * S_ ** 2)) for A, C_, S_ in B)


def _fd_jacobian(F, v0, h=1e-6):
    d = v0.size; J = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        J[:, j] = (F(v0 + e) - F(v0 - e)) / (2 * h)
    return J


def boundary_job():
    """Part G. (a) A non-replicated stationary point: the four-bump target B4, K₀ = 3, q = 1/3, m = 6, the third class split over the last two
    bumps. The stationary point of the η = ∞ form (the leak toward the class barycentres, the plan the class partition) by least squares from
    the exact six-kernel fit; checked for η = σ₀, 10σ₀ as a stationary point of the joint flow with the prior at the barycentres (the drift, the
    plan's optimality), with its leak force, the unregularised gradient there and the largest real part of the joint Jacobian (plan frozen,
    central differences). Then the joint flow from the exact split fit (η = ∞, σ₀; S_G slow units): its state at S_G, the collapsed η = ∞
    dynamics (the unregularised K₀-atom gradient flow, LSODA) continued from it to the times T_CONT, and the minimiser of the loss on the face
    where the merged atom's width is s_max, with the loss's outward slope there. (b) Unbounded solutions: the target B_SAT, the K₀-atom system (m = 3) from the instance's exact atoms with the prior there,
    η = 0, σ₀, ∞, LSODA at the times T_G. (c) The exact representation, m = 6, with the prior's atoms permuted cyclically."""
    from scipy.optimize import least_squares
    from scipy.integrate import solve_ivp
    t0 = time.time()
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    out = {}
    # (a)
    m = 6
    Th0 = net.unconstrained(np.array([[3.0], [3.0], [-2.1], [-2.1], [2.4], [2.4]]), np.array([[-1.5], [-1.5], [0.0], [0.0], [1.3], [2.3]]),
                            np.array([0.5, 0.5, 0.4, 0.4, 0.4, 0.4]))
    loop4 = OT.OTPriorLoop(net, _bumps(B4), x_ref, PER, G, SIGMA0, Th0[[0, 2, 4]], n_quad=N_Q)
    lab = np.array([0, 0, 1, 1, 2, 2]); gam = np.zeros((m, 3)); gam[np.arange(m), lab] = 1.0 / m
    split = lambda T_: float(np.sqrt(np.sum((T_[4] - T_[5]) ** 2 / G)))
    unreg = lambda T_: float(np.abs(loop4.drift(T_, T_) / np.sqrt(G)[None, :]).max())

    def F_inf(v, g_=gam):
        T_ = v.reshape(m, 3)
        return (loop4.drift(T_, m * (g_ @ AP.barycentres(T_, g_, Q))) / np.sqrt(G)[None, :]).ravel()

    def examine(T_, g_):
        Z_ = AP.barycentres(T_, g_, Q)
        leak = SIGMA0 * np.sqrt(np.sum((T_ - m * (g_ @ Z_)) ** 2 / G[None, :], axis=1))
        per_eta = []
        for eta in (SIGMA0, 10 * SIGMA0):
            apl = AP.AdaptivePriorLoop(loop4, eta)
            gp, W2 = AP.plan_to(T_, Z_, Q, G)
            dT, dZ = apl.drift(T_, Z_, g_)

            def Fj(v):
                x, y = apl.drift(v[:3 * m].reshape(m, 3), v[3 * m:].reshape(3, 3), g_)
                return np.concatenate([x.ravel(), y.ravel()])
            ev = np.linalg.eigvals(_fd_jacobian(Fj, np.concatenate([T_.ravel(), Z_.ravel()])))
            per_eta.append(dict(eta=eta, plan_is_classes=bool(np.array_equal(np.argmax(gp, axis=1), np.argmax(g_, axis=1))), W2=W2,
                                dT=float(np.abs(dT).max()), dZ=float(np.abs(dZ).max()), max_re=float(ev.real.max())))
        return dict(split=split(T_), leak=leak, unreg=unreg(T_), fit=loop4.fit_error(T_), L=loop4.L_bar(T_), per_eta=per_eta)
    sol = least_squares(F_inf, Th0.ravel(), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
    out['a'] = dict(res=float(np.abs(F_inf(sol.x)).max()), point=examine(sol.x.reshape(m, 3), gam), exact_split=split(Th0), exact_fit=loop4.fit_error(Th0),
                    exact_unreg=unreg(Th0))
    flows = []
    for eta in (np.inf, SIGMA0):
        apl = AP.AdaptivePriorLoop(loop4, eta)
        rec = dict(s=[], fit=[], J=[], split=[])

        def cb(s_, T_, Z_, ad, rec=rec, apl=apl):
            rec['s'].append(s_); rec['fit'].append(loop4.fit_error(T_)); rec['J'].append(apl.energy(T_, Z_)); rec['split'].append(split(T_))
        Th, Z, _ = apl.run(Th0.copy(), Th0[[0, 2, 4]].copy(), S_G, DT, plan_every=PLAN_EVERY, rec_every=16000, callback=cb)
        g_, _ = AP.plan_to(Th, Z, Q, G)
        atoms = AP.barycentres(Th, g_, Q)
        spread = float(np.sqrt(np.sum((Th - m * (g_ @ atoms)) ** 2 / G[None, :], axis=1)).max())
        kk = int(np.argmax(net.physical(atoms)[2]))                          # the merged atom: the widest
        # the collapsed η = ∞ dynamics continued: the unregularised K₀-atom gradient flow from the class atoms (three particles of weight 1/3)
        rhs = lambda t, v: loop4.drift(v.reshape(3, 3), v.reshape(3, 3)).ravel()
        sol = solve_ivp(rhs, (0.0, T_CONT[-1]), atoms.ravel(), method='LSODA', t_eval=np.array(T_CONT), rtol=1e-10, atol=1e-12)
        cont = []
        for k in range(len(sol.t)):
            at = sol.y[:, k].reshape(3, 3)
            cont.append(dict(s=float(sol.t[k]), kappa=float(at[kk, 2]), gap=float(net.s_max - net.physical(at)[2][kk]), fit=loop4.fit_error(at)))
        flows.append(dict(eta=eta, rec={k: np.array(v) for k, v in rec.items()}, classes=tuple(int(x) for x in np.argmax(g_, axis=1)), spread=spread, kk=kk,
                          atoms=tuple(np.asarray(x) for x in net.physical(atoms)), kappa_end=float(atoms[kk, 2]), fit_end=loop4.fit_error(Th),
                          drift_end=float(np.abs(rhs(0.0, atoms.ravel()) / np.tile(np.sqrt(G), 3)).max()), cont=cont, status=int(sol.status)))
    out['a_flow'] = flows
    # the minimiser of the K₀-atom loss on the face where the merged atom's width is s_max (physical coordinates, BFGS), and the outward slope there
    from scipy.optimize import minimize
    kk = flows[0]['kk']
    a0, c0, s0 = (np.asarray(x, float).ravel() for x in flows[0]['atoms'])
    X_ = loop4.X; F_ = loop4.F

    def L_phys(a, c, sw):
        fh = (np.exp(-(X_[:, None] - c[None, :]) ** 2 / (2 * sw[None, :] ** 2)) @ a) / 3.0
        return 0.5 * float(np.mean((F_ - fh) ** 2))
    free = [k for k in range(3) if k != kk]

    def L_face(v):
        a, c = v[:3], v[3:6]; sw = np.empty(3); sw[kk] = net.s_max; sw[free] = v[6:]
        return L_phys(a, c, sw)
    res = minimize(L_face, np.concatenate([a0, c0, s0[free]]), method='BFGS', options=dict(gtol=1e-12, maxiter=20000))
    a, c = res.x[:3], res.x[3:6]; sw = np.empty(3); sw[kk] = net.s_max; sw[free] = res.x[6:]
    hh = 1e-6; sp, sm_ = sw.copy(), sw.copy(); sp[kk] += hh; sm_[kk] -= hh
    slope = (L_phys(a, c, sp) - L_phys(a, c, sm_)) / (2 * hh)
    fit_face = float(np.sqrt(2 * L_phys(a, c, sw) / np.mean(F_ ** 2)))
    out['face'] = dict(fit=fit_face, slope=float(slope), grad=float(np.abs(res.jac).max()), a=a, c=c, s=sw,
                       inside=bool(np.all(np.abs(a) < net.a_max) and np.all(np.abs(c) < net.c_max) and np.all(sw[free] > net.s_min) and np.all(sw[free] < net.s_max)))
    # (b)
    loop2 = OT.OTPriorLoop(net, _bumps(B_SAT), x_ref, PER, G, SIGMA0, representation(net, 0.0), n_quad=N_Q)
    start = representation(net, 0.0); gamma = np.eye(3) / 3.0
    runs = []
    for eta in (0.0, SIGMA0, np.inf):
        apl = AP.AdaptivePriorLoop(loop2, eta)

        def rhs(t, v, apl=apl, eta=eta):
            Th = v[:9].reshape(3, 3); Z = v[9:].reshape(3, 3)
            dT, dZ = apl.drift(Th, Th if np.isinf(eta) else Z, gamma)
            return np.concatenate([dT.ravel(), (np.zeros(9) if np.isinf(eta) else dZ.ravel())])
        sol = solve_ivp(rhs, (0.0, T_G[-1]), np.concatenate([start.ravel(), start.ravel()]), method='LSODA', t_eval=np.array(T_G), rtol=1e-10, atol=1e-12)
        rows = []
        for k in range(len(sol.t)):
            Th = sol.y[:9, k].reshape(3, 3); a, c, s_ = net.physical(Th)
            rows.append(dict(s=float(sol.t[k]), th_a0=float(Th[0, 0]), z_a0=(float(sol.y[9, k]) if not np.isinf(eta) else float('nan')), a0=float(a[0, 0]),
                             fit=loop2.fit_error(Th), drift=float(np.abs(rhs(sol.t[k], sol.y[:, k])).max())))
        runs.append(dict(eta=eta, status=int(sol.status), rows=rows))
    out['b'] = dict(runs=runs, a_max=net.a_max)
    # (c)
    net_, exact, prior, loop = setup()
    apl = AP.AdaptivePriorLoop(loop, SIGMA0)
    Th = np.repeat(exact, 2, axis=0); Zp = exact[[1, 2, 0]]
    g_, W2 = AP.plan_to(Th, Zp, Q, G)
    dT, dZ = apl.drift(Th, Zp, g_)
    out['c'] = dict(classes=tuple(int(x) for x in np.argmax(g_, axis=1)), W2=W2, dT=float(np.abs(dT).max()), dZ=float(np.abs(dZ).max()),
                    offset=float(np.abs(Zp - exact).max()))
    print(f"[{time.time() - t0:.0f}s] G", flush=True)
    return out


A_BOUND, C_BOUND, K_BOUND = 0.9, 0.9, 6.0                                  # Part H: the box, amplitudes and centres to 90% of a_max, c_max
S_H = 1e5                                                                  # Part H: the projected midpoint run at η = σ₀
DT_H = (0.25, 0.125)                                                       # Part H(a): the projected midpoint rule at two steps
A_SCAN = (0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99)                # Part H(c): the box's amplitude bound, a fraction of a_max
DT_TRAJ, DT_REF, S_TRAJ = (0.25, 0.125, 0.0625), 0.015625, 300.0            # Part H(a): the trajectory error against a fine reference


def projected_box(net, a_bound=A_BOUND):
    lo = np.array([-net.a_max * np.arctanh(a_bound), -net.c_max * np.arctanh(C_BOUND), -K_BOUND])
    return lo, -lo


def projected_job():
    """Part H. (a) Part G(b)'s target B_SAT, the K₀-atom system from the instance's exact atoms with the prior there, the prior projected onto
    the box: η = ∞ by LSODA (the prior is the projected class mean: a Lipschitz vector field, with a kink where an atom meets a face) at the
    times T_G; η = σ₀ by the projected midpoint rule of AdaptivePriorLoop.run (positions projected inside the step; the projected flow's
    velocity jumps at the faces) to S_H at the two steps DT_H, the projected flow's drift at the end checked to be ≤ 10⁻¹²; with
    B = (max|f| + a_max)·L_ϑ^G, L_ϑ^G ≤ the closed-form bound of |G^{1/2}∇_ϑΦ̃| (so B is an upper bound), and the largest distance of the
    particles from the box. (b) Part C's K₀-atom loop at the misplaced prior (η = σ₀, S_C slow units) with and without the box: the largest
    difference of the atoms and of the prior, and the smallest margin of the prior and of the particles to a face along the run, relative to
    the box's width in that coordinate. (c) The box's size: (a)'s η = ∞ system (LSODA to T_G[-1]) for the amplitude bounds A_SCAN. (d) Part
    G(a)'s four-bump target with the box at η = ∞: the joint flow (m = 6) from the exact split fit for S_G slow units, then the collapsed
    dynamics with the box (the K₀-atom flow with the leak toward the projection of the atoms, LSODA) to the times T_CONT counted from its end:
    the merged atom's width coordinate κ and its prior atom's, the fit error, the drift."""
    from scipy.integrate import solve_ivp
    t0 = time.time()
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    box = projected_box(net)
    start = representation(net, 0.0)
    loop2 = OT.OTPriorLoop(net, _bumps(B_SAT), x_ref, PER, G, SIGMA0, start, n_quad=N_Q)
    gamma = np.eye(3) / 3.0
    LG = float(np.sqrt(G[0] + G[1] * (net.a_max * np.exp(-0.5) / net.s_min) ** 2 + G[2] * (net.a_max * (net.s_max - net.s_min) / (2 * np.e * net.s_min)) ** 2))
    B = (float(np.max(np.abs(loop2.F))) + net.a_max) * LG
    dist_box = lambda Th, bx=box: float(np.sqrt(np.max(np.sum((Th - np.clip(Th, bx[0][None, :], bx[1][None, :])) ** 2 / G[None, :], axis=1))))

    def lsoda_inf(bx, T_eval):
        apl = AP.AdaptivePriorLoop(loop2, np.inf, box=bx)

        def rhs(t, v):
            Th = v.reshape(3, 3)
            return apl.drift(Th, Th, gamma)[0].ravel()
        t_dense = np.concatenate([[0.0], np.geomspace(1.0, T_eval[-1], 120)])
        sol = solve_ivp(rhs, (0.0, T_eval[-1]), start.ravel(), method='LSODA', t_eval=np.unique(np.concatenate([t_dense, T_eval])), rtol=1e-10, atol=1e-12)
        rows = []
        for t in T_eval:
            k = int(np.argmin(np.abs(sol.t - t)))
            Th = sol.y[:, k].reshape(3, 3); a, c, s_ = net.physical(Th)
            rows.append(dict(s=float(sol.t[k]), th_a0=float(Th[0, 0]), z_a0=float(apl.project(Th)[0, 0]), a0=float(a[0, 0]), fit=loop2.fit_error(Th),
                             drift=float(np.abs(rhs(sol.t[k], sol.y[:, k])).max())))
        return dict(eta=np.inf, status=int(sol.status), rows=rows, max_dist=max(dist_box(sol.y[:, k].reshape(3, 3), bx) for k in range(sol.y.shape[1])))
    runs = [lsoda_inf(box, T_G)]
    # η = σ₀: the projected midpoint rule at two steps
    ends = []
    for dt in DT_H:
        apl = AP.AdaptivePriorLoop(loop2, SIGMA0, box=box)
        rec = dict(s=[], th=[], z=[], fit=[], drift=[], dist=[])

        def cb(s_, T_, Z_, ad, rec=rec, apl=apl):
            g_, _ = AP.plan_to(T_, Z_, Q, G)
            rec['dist'].append(dist_box(T_))
            if any(abs(s_ - t) < 1e-9 for t in T_G):
                dT, dZ = apl.drift(T_, Z_, g_)
                rec['s'].append(s_); rec['th'].append(T_.copy()); rec['z'].append(Z_.copy()); rec['fit'].append(loop2.fit_error(T_))
                rec['drift'].append(float(max(np.abs(dT).max(), np.abs(dZ).max())))
        Th_end, Z_end, _ = apl.run(start.copy(), start.copy(), S_H, dt, plan_every=PLAN_EVERY, rec_every=int(round(100.0 / dt)), callback=cb)
        g_, _ = AP.plan_to(Th_end, Z_end, Q, G)
        dT, dZ = apl.drift(Th_end, Z_end, g_)
        end_drift = float(max(np.abs(dT).max(), np.abs(dZ).max()))
        assert end_drift <= 1e-12, (dt, end_drift)                         # a stall of the prior inside 𝒵 would leave a drift (D3n-2)
        rows = []
        for k in range(len(rec['s'])):
            Th = rec['th'][k]; a, c, s_ = net.physical(Th)
            rows.append(dict(s=rec['s'][k], th_a0=float(Th[0, 0]), z_a0=float(rec['z'][k][0, 0]), a0=float(a[0, 0]), fit=rec['fit'][k], drift=rec['drift'][k]))
        ends.append(dict(dt=dt, Th=Th_end, Z=Z_end, drift=end_drift, gap=float(box[1][0] - Z_end[0, 0]), rows=rows, max_dist=max(rec['dist'])))
    runs.insert(0, dict(eta=SIGMA0, status=0, rows=ends[0]['rows'], max_dist=ends[0]['max_dist']))
    dt_check = dict(dts=DT_H, diff=float(max(np.abs(ends[0]['Th'] - ends[1]['Th']).max(), np.abs(ends[0]['Z'] - ends[1]['Z']).max())),
                    drifts=[e['drift'] for e in ends], gaps=[e['gap'] for e in ends], max_dist2=ends[1]['max_dist'])
    # the scheme's consistency: the trajectory at s = S_TRAJ against a reference at DT_REF (the end states above share the fixed point)
    traj = {}
    for dt in DT_TRAJ + (DT_REF,):
        apl = AP.AdaptivePriorLoop(loop2, SIGMA0, box=box)
        traj[dt] = apl.run(start.copy(), start.copy(), S_TRAJ, dt, plan_every=max(1, int(round(PLAN_EVERY * DT / dt))), rec_every=10 ** 9)[:2]
    dt_check['traj'] = {dt: float(max(np.abs(traj[dt][0] - traj[DT_REF][0]).max(), np.abs(traj[dt][1] - traj[DT_REF][1]).max())) for dt in DT_TRAJ}
    # (b) the misplaced prior: the box never reached
    net_, exact, prior, loop = setup()
    star = loop.stationary_atoms(prior)[0]
    out_b = {}
    width = box[1] - box[0]
    marg = dict(Z=np.inf, Th=np.inf)

    def cb_m(s_, T_, Z_, ad):
        marg['Z'] = min(marg['Z'], float(np.min(np.minimum(Z_ - box[0][None, :], box[1][None, :] - Z_) / width[None, :])))
        marg['Th'] = min(marg['Th'], float(np.min(np.minimum(T_ - box[0][None, :], box[1][None, :] - T_) / width[None, :])))
    for name, bx in (('free', None), ('box', box)):
        apl = AP.AdaptivePriorLoop(loop, SIGMA0, box=bx)
        Th, Z, _ = apl.run(star.copy(), prior.copy(), S_C, DT, plan_every=PLAN_EVERY, rec_every=1 if bx is not None else 10 ** 9,
                           callback=cb_m if bx is not None else None)
        out_b[name] = (Th, Z)
    diff = float(max(np.abs(out_b['free'][0] - out_b['box'][0]).max(), np.abs(out_b['free'][1] - out_b['box'][1]).max()))
    inside = bool(np.all(np.abs(out_b['free'][1]) < box[1][None, :]) and np.all(np.abs(out_b['free'][0]) < box[1][None, :]))
    # (c) the box's size
    scan = []
    for ab in A_SCAN:
        bx = projected_box(net, ab)
        r = lsoda_inf(bx, (T_G[-1],))
        scan.append(dict(a_bound=ab, face=float(bx[1][0]), **{k: r['rows'][-1][k] for k in ('th_a0', 'z_a0', 'fit', 'drift')}, max_dist=r['max_dist'], status=r['status'],
                         start_inside=bool(np.all(start >= bx[0][None, :]) and np.all(start <= bx[1][None, :]))))
    # (d) Part G(a)'s four-bump target with the box (η = ∞)
    m = 6
    Th0 = net.unconstrained(np.array([[3.0], [3.0], [-2.1], [-2.1], [2.4], [2.4]]), np.array([[-1.5], [-1.5], [0.0], [0.0], [1.3], [2.3]]),
                            np.array([0.5, 0.5, 0.4, 0.4, 0.4, 0.4]))
    loop4 = OT.OTPriorLoop(net, _bumps(B4), x_ref, PER, G, SIGMA0, Th0[[0, 2, 4]], n_quad=N_Q)
    apl4 = AP.AdaptivePriorLoop(loop4, np.inf, box=box)
    Th, Z, _ = apl4.run(Th0.copy(), Th0[[0, 2, 4]].copy(), S_G, DT, plan_every=PLAN_EVERY, rec_every=10 ** 9)
    g_, _ = AP.plan_to(Th, Z, Q, G)
    atoms = AP.barycentres(Th, g_, Q)
    spread = float(np.sqrt(np.sum((Th - m * (g_ @ atoms)) ** 2 / G[None, :], axis=1)).max())
    kk = int(np.argmax(net.physical(atoms)[2]))
    rhs4 = lambda t, v: loop4.drift(v.reshape(3, 3), apl4.project(v.reshape(3, 3))).ravel()
    sol = solve_ivp(rhs4, (0.0, T_CONT[-1]), atoms.ravel(), method='LSODA', t_eval=np.array(T_CONT), rtol=1e-10, atol=1e-12)
    cont = []
    for k in range(len(sol.t)):
        at = sol.y[:, k].reshape(3, 3)
        cont.append(dict(s=float(sol.t[k]), kappa=float(at[kk, 2]), kappa_prior=float(apl4.project(at)[kk, 2]), fit=loop4.fit_error(at),
                         drift=float(np.abs(rhs4(0.0, sol.y[:, k]) / np.tile(np.sqrt(G), 3)).max())))
    G4 = dict(classes=tuple(int(x) for x in np.argmax(g_, axis=1)), spread=spread, kk=kk, kappa_start=float(atoms[kk, 2]), fit_start=loop4.fit_error(Th), cont=cont,
              status=int(sol.status), face=float(box[1][2]))
    print(f"[{time.time() - t0:.0f}s] H", flush=True)
    return dict(box=box, LG=LG, B=B, runs=runs, dt_check=dt_check, diff_b=diff, inside_b=inside, margin_b=dict(marg), scan=scan, G4=G4, sel=selection_job())


def selection_job():
    """Part H(e): at η = ∞ the plan is Lloyd-consistent (optimal for its own class means), not jointly optimal. A configuration of m = 6
    particles (K₀ = 3 atoms of mass 1/3; zero amplitudes, the target 0.02 × the three bumps; Part H's box), found by a seeded search over
    planar configurations in the whitened (γ, κ) coordinates for a Lloyd-consistent, non-optimal pairing along whose pure-leak flow the
    jointly optimal pairing's energy rises; the prior started at the class means of that pairing. Along AdaptivePriorLoop's run (plan 'ot',
    η = ∞): the selected energy 𝒥(Θ, Z) (Z the class means of the current plan, optimal for ν_Z) and the reduced energy min_Z𝒥 (the minimum
    over the 15 balanced pairings, with the box's q_j dist²(T̂_j, 𝒵) term)."""
    from scipy.optimize import linear_sum_assignment

    def pairings(rem):
        if not rem:
            yield []
            return
        a = rem[0]
        for b in rem[1:]:
            for pp in pairings([x for x in rem if x not in (a, b)]):
                yield [(a, b)] + pp
    PAIRS = list(pairings(list(range(6))))
    energy = lambda Y, pp: sum(np.sum((Y[a] - Y[b]) ** 2) / 2 for a, b in pp)
    rng = np.random.default_rng(0)
    Y = pair = None
    for trial in range(20000):
        Yt = rng.uniform(-2, 2, (6, 2))
        Es = np.array([energy(Yt, pp) for pp in PAIRS])
        order = np.argsort(Es)
        if Es[order[1]] - Es[order[0]] < 1e-3:
            continue
        act = PAIRS[order[0]]
        for k in order[1:]:
            pp = PAIRS[k]
            Zt = np.array([(Yt[a] + Yt[b]) / 2 for a, b in pp])
            Cm = np.sum((Yt[:, None, :] - Zt[None, :, :]) ** 2, axis=2)
            rep_ = np.repeat(np.arange(3), 2)
            r_, c_ = linear_sum_assignment(Cm[:, rep_])
            if sum(Cm[a, j] + Cm[b, j] for j, (a, b) in enumerate(pp)) > Cm[r_, rep_[c_]].sum() + 1e-12:
                continue                                                   # not Lloyd-consistent
            V = np.zeros_like(Yt)
            for a, b in pp:
                mu = (Yt[a] + Yt[b]) / 2; V[a] = -(Yt[a] - mu); V[b] = -(Yt[b] - mu)
            if sum(np.dot(Yt[a] - Yt[b], V[a] - V[b]) for a, b in act) > 1e-3:
                Y, pair = Yt, pp
            break
        if Y is not None:
            break
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    lo, hi = projected_box(net)
    Th0 = np.column_stack([np.zeros(6), np.sqrt(G[1]) * Y[:, 0], np.sqrt(G[2]) * Y[:, 1]])
    Z0 = np.array([Th0[list(pp)].mean(axis=0) for pp in pair])
    loop = OT.OTPriorLoop(net, _bumps([(0.02 * A, C, S) for A, C, S in BUMPS]), x_ref, PER, G, SIGMA0, Z0, n_quad=N_Q)
    apl = AP.AdaptivePriorLoop(loop, np.inf, plan='ot', box=(lo, hi))

    def reduced(Th):
        best = np.inf
        for pp in PAIRS:
            e = 0.0
            for a, b in pp:
                T = (Th[a] + Th[b]) / 2
                e += (np.sum((Th[a] - T) ** 2 / G) + np.sum((Th[b] - T) ** 2 / G)) / 6 + (1 / 3) * np.sum((T - np.clip(T, lo, hi)) ** 2 / G)
            best = min(best, loop.L_bar(Th) + 0.5 * SIGMA0 * e)
        return best
    rec = dict(s=[], sel=[], red=[], optimal=[])

    def cb(s_, Th, Z, ad):
        g, _ = apl.plan(Th, Z)
        rec['s'].append(s_); rec['sel'].append(apl.energy(Th, Z)); rec['red'].append(reduced(Th))
        rec['optimal'].append(bool(abs(rec['sel'][-1] - rec['red'][-1]) < 1e-12))
    apl.run(Th0.copy(), Z0.copy(), 400.0, 0.25, plan_every=4, rec_every=40, callback=cb)
    out = {k: np.array(v) for k, v in rec.items()}
    out.update(trial=trial, pairing=pair)
    print("[H selection]", flush=True)
    return out


def report(lin, C, D, E=None, t0=None, F=None, Gb=None, H=None):
    h, V = lin['h'], lin['V']
    lines = ["# Theorem D‴.23: adaptive priors — a learned prior removes the leak's bias (`exp_adaptive_prior.py`)", "",
             f"The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit = {SIGMA0:.5f}; the OT leak toward the prior with centres shifted by +{SHIFT:g} (masses 1/3); "
             f"the adaptive prior moves its atoms toward the class barycentres at the rate η. Midpoint rule dt = {DT:g} slow units.", "",
             "## A. The linear theory at the exact representation", "",
             "The whitened weighted Gauss–Newton eigenvalues h_j of the loss (the collective block without the leak), the fraction σ₀/(h_j + σ₀) of the prior's error that a fixed prior "
             "keeps in mode j, and the slower rate λ_j^− of the collective pair for each η (η = ∞: λ_j^− = h_j). The separations decay at σ₀ in every case.", "",
             "| j | h_j | σ₀/(h_j + σ₀) | " + " | ".join(f"λ_j^− at η = {('∞' if np.isinf(e) else f'{e:.5f}')}" for e in ETAS) + " |",
             "|---|---|---|" + "---|" * len(ETAS)]
    for j in range(len(h)):
        lines.append(f"| {j + 1} | {h[j]:.4e} | {SIGMA0 / (h[j] + SIGMA0):.3f} | " + " | ".join(f"{lin['rates'][e][j]:.3e}" for e in ETAS) + " |")
    lines += ["", f"A fixed prior's collective rates are h_j + σ₀ = {h[0] + SIGMA0:.4e}–{h[-1] + SIGMA0:.4e}."]
    if 'radius' in lin:
        rb = lin['radius']
        lines += ["", f"(v)'s radius. A sampled lower bound on the Lipschitz constant C_N of the collective Jacobian in the whitened weighted norm (the unregularised K₀-atom drift by "
                  f"central differences at Θ*(q) and at {rb['n_pts'] - 1} random points at whitened distances 0.005–0.05 from the exact representation, where the Jacobian is −H_L to "
                  f"{rb['J0_err']:.1e}): C_N ≥ {rb['C']:.4f}. The proof's radius is then at most min{{σ₀, λ_1^−}}/(4C_N) (the plan-margin term can only lower it), against the "
                  "norm |||·|||_η of Part C's start (Θ*(q) and the shifted prior; at η = ∞ the prior is the class barycentres and the norm has no ζ part):", "",
                  "| η | λ_1^− | min{σ₀, λ_1^−}/σ₀ | upper bound on ρ | |||start|||_η | start / bound |", "|---|---|---|---|---|---|"]
        for r in rb['rows']:
            en = '∞' if np.isinf(r['eta']) else f"{r['eta']:.5f}"
            lines.append(f"| {en} | {r['lam1']:.3e} | {r['lam0'] / SIGMA0:.4f} | {r['rho']:.2e} | {r['start']:.3e} | {r['start'] / r['rho']:.0f} |")
        en_ = lambda e: '∞' if np.isinf(e) else f"{e:.5f}"
        dd = [f"{r['d_min'] / r['rho']:.1e}–{r['d_max'] / r['rho']:.1e} at η = {en_(r['eta'])}" for r in rb['rows'] if 'd_min' in r]
        ff = []
        if F:
            for r in rb['rows']:
                st = [x['start'] for x in F if x['eta'] == r['eta'] and 'start' in x]
                if st:
                    lo, hi = f"{min(st) / r['rho']:.0f}", f"{max(st) / r['rho']:.0f}"
                    ff.append((lo if lo == hi else f"{lo}–{hi}") + f" at η = {en_(r['eta'])}")
        lines += ["", "The runs start outside the ball that the proof of (v) covers: Part C's start by the factors of the table; Part D's warm clouds (largest separations "
                  "of order 5 at s = 0) by " + "; ".join(dd) + (("; Part F's separated starts by " + "; ".join(ff)) if ff else "") + ". They probe the dynamics beyond it."]
    lines += ["",
              "## B. The fixed prior's bias (shifted prior)", "",
              f"Θ*(q) by Newton (residual {lin['newton_res']:.1e}; diagonal plan ε = {lin['diag'][0]:.3f} < Δ_ν/2 = {lin['diag'][1]:.3f}); fit error at Θ*(q) {lin['fit_star']:.4e} (0 at the exact representation). "
              f"Mode amplitudes (components on the eigenvectors of H_L, whitened weighted coordinates) of the prior's error x⁰, of the actual bias y = Θ*(q) − exact, and of the first-order prediction "
              f"σ₀(H_L + σ₀I)⁻¹x⁰; ‖y − prediction‖/‖y‖ = {np.linalg.norm(lin['ya'] - lin['yp']) / np.linalg.norm(lin['ya']):.3e}.", "",
              "| j | x⁰_j | y_j (actual) | predicted | actual/x⁰ | actual/predicted |", "|---|---|---|---|---|---|"]
    a0, aa, ap = V.T @ lin['y0'], V.T @ lin['ya'], V.T @ lin['yp']
    for j in range(len(h)):
        lines.append(f"| {j + 1} | {a0[j]:+.4e} | {aa[j]:+.4e} | {ap[j]:+.4e} | {aa[j] / a0[j]:.3f} | {aa[j] / ap[j]:.3f} |")
    if 'conv' in lin:
        lines += ["", "The first-order formula against the shift of the prior's centres (Newton residual of Θ*(q); the relative error in norm; actual/predicted in the modes 1–4):", "",
                  "| shift | Newton residual | ‖y − prediction‖/‖y‖ | actual/predicted, modes 1, 2, 3, 4 |", "|---|---|---|---|"]
        for c in lin['conv']:
            lines.append(f"| {c['shift']:g} | {c['res']:.1e} | {c['rel']:.3e} | " + ", ".join(f"{x:.3f}" for x in c['mode_ratio'][:4]) + " |")
    # C
    lines += ["", "## C. The K₀-atom adaptive loop from Θ*(q) with the prior at the shifted centres", "",
              "Mode amplitudes of the atoms' deviation from the exact representation at s = 0, 1000, 3000, 6000 against the exact solution of the linear pairs (started from the actual Θ*(q) "
              "and the actual prior); the fit error and the energy.", "",
              "| η | mode | amplitude at s = 0, 1000, 3000, 6000 | linear prediction at the same times | ratio at s = 1000, 3000, 6000 | log-slope on [3000, 6000] / λ_j^− | fit error at s = 0, 1000, 3000, 6000 | 𝒥 at s = 6000 |",
              "|---|---|---|---|---|---|---|---|"]
    low = []                                                               # amplitudes below FLOOR at s ≥ 3000 (the remainder's level)
    slopes, ratios, drift = [], [], {}                                     # every printed log-slope and ratio; the ratios' change from s = 1000 on, per η
    for c in C:
        eta = c['eta']
        A = c['y'] @ V; Zm = c['z'] @ V
        late = np.abs(A[c['s'] >= 3000.0]); low += list(late[late < FLOOR])
        idx = [int(np.argmin(np.abs(c['s'] - t))) for t in (0.0, 1000.0, 3000.0, 6000.0)]
        lam = AP.pair_rates(h, SIGMA0, eta)[0]
        for j in range(len(h)):
            pred_y, _ = AP.pair_solution(h[j], SIGMA0, eta, A[0, j], Zm[0, j], c['s'][idx])
            ratio = ", ".join(f"{A[i, j] / p_:.3f}" if min(abs(A[i, j]), abs(p_)) >= FLOOR else "—" for i, p_ in zip(idx[1:], pred_y[1:]))
            slope = (-(np.log(abs(A[idx[3], j])) - np.log(abs(A[idx[2], j]))) / (c['s'][idx[3]] - c['s'][idx[2]]) / lam[j]
                     if min(abs(A[idx[2], j]), abs(A[idx[3], j])) >= FLOOR else float('nan'))
            rv = [A[i, j] / p_ if min(abs(A[i, j]), abs(p_)) >= FLOOR else np.nan for i, p_ in zip(idx[1:], pred_y[1:])]
            ratios += [x for x in rv if np.isfinite(x)]
            if np.isfinite(slope):
                slopes.append(slope)
            later = [abs(x / rv[0] - 1) for x in rv[1:] if np.isfinite(x)]
            if np.isfinite(rv[0]) and later:
                drift[eta] = max(drift.get(eta, 0.0), max(later))
            lines.append(f"| {'∞' if np.isinf(eta) else f'{eta:.5f}'} | {j + 1} | " + ", ".join(f"{A[i, j]:+.3e}" for i in idx) + " | " + ", ".join(f"{x:+.3e}" for x in pred_y)
                         + f" | {ratio} | " + (f"{slope:.3f}" if np.isfinite(slope) else "—")
                         + " | " + (", ".join(f"{c['fit'][i]:.3e}" for i in idx) if j == 0 else "") + " | " + (f"{c['J'][-1]:.3e}" if j == 0 else "") + " |")
    lines += ["", "Ratio: the amplitude over the linear prediction (constant in time when the rate is the predicted one; its offset from 1 is set in the nonlinear first phase, "
              f"the shift being finite). Log-slope: the decay rate of the amplitude on [3000, 6000] over λ_j^−. Both only where the amplitudes are above {FLOOR:.0e}: below it the "
              f"faster modes stop at the level of the second-order remainder driven by the slow modes (the amplitudes below {FLOOR:.0e} at s ≥ 3000, every mode and η: "
              f"{min(low):.1e}–{max(low):.1e}).", "",
              f"Over every mode and η above {FLOOR:.0e}: log-slopes {min(slopes):.3f}–{max(slopes):.3f} of λ_j^− ({len(slopes)} modes), ratios {min(ratios):.3f}–{max(ratios):.3f}; "
              "the largest change of a ratio from s = 1000 to s = 3000 or 6000 (modes with a ratio at s = 1000 and a later one): "
              + ", ".join(f"{100 * drift[e]:.2f}% at η = {('∞' if np.isinf(e) else f'{e:.5f}')}" for e in sorted(drift)) + "."]
    # D
    lines += ["", "## D. m-particle systems from the warm clouds (shifted prior)", "",
              f"Sizes m = {', '.join(map(str, MS))} and M = {M_REF}, {N_SEEDS} seeds, {S_D:g} slow units, plan every {PLAN_EVERY} steps. Means over the seeds (min–max). ‖ȳ‖: the whitened weighted norm of "
              "the class means' deviation from the exact representation (classes by the plan); W²: to the exact representation's measure. Gated: the prior is frozen while mass is in transit "
              "(some particle assigned to an atom that is not its nearest), so that the class barycentres it moves to contain no particle of another cluster ('transit'), and in addition, the first time, "
              "until the fixed-prior drift is below 10⁻⁵ ('settled': the adaptation starts from the fixed-prior stationary point); t_open is the first time it adapts. "
              "In parentheses: the K₀-atom run of Part C at s (ungated) or at s − t_open (gated, mean t_open).", "",
              "| η | gated | m | t_open (seeds) | fit error at s = 0, 500, 3000 | W²(μ_Θ, μ_exact) at s = 0, 500, 3000 | ‖ȳ‖ at s = 500, 3000 | prior's distance at s = 3000 | final counts |",
              "|---|---|---|---|---|---|---|---|---|"]
    Cmap = {c['eta']: c for c in C}
    combos = [(0.0, False), (SIGMA0, False), (np.inf, False), (SIGMA0, True), (np.inf, True), (SIGMA0, 'settled'), (np.inf, 'settled')]
    for eta, gate in combos:
        for m in list(MS) + [M_REF]:
            rs = [d for d in D if d['eta'] == eta and d.get('gate', False) == gate]
            if not rs:
                continue
            recs = [d['rec'][m] for d in rs]
            at = lambda key, t: np.array([r[key][int(np.argmin(np.abs(r['s'] - t)))] for r in recs])
            fit = [at('fit', t) for t in (0.0, 500.0, 3000.0)]
            w2 = [at('W2ex', t) for t in (0.0, 500.0, 3000.0)]
            yb = [np.array([np.linalg.norm(v) for v in at('ybar', t)]) for t in (500.0, 3000.0)]
            tops = [d.get('t_open', {}).get(m, 0.0) for d in rs] if gate else [0.0 for d in rs]
            ref = ""
            if eta in Cmap:
                c = Cmap[eta]; sh = float(np.mean(tops))
                ref = " (" + ", ".join(f"{np.linalg.norm(c['y'][int(np.argmin(np.abs(c['s'] - max(t - sh, 0.0))))]):.3e}" for t in (500.0, 3000.0)) + ")"
            elif eta == 0.0:
                ref = f" (Θ*(q): {np.linalg.norm(lin['ya']):.3e})"
            zd = at('zdist', 3000.0)
            cnt = "; ".join(str(tuple(int(x) for x in r['counts'][-1])) for r in recs)
            to = ", ".join(f"{x:g}" for x in tops) if gate else "—"
            gname = {False: 'no', True: 'transit', 'settled': 'settled'}[gate]
            lines.append(f"| {('∞' if np.isinf(eta) else f'{eta:.5f}')} | {gname} | {m} | {to} | " + ", ".join(f"{x.mean():.3e}" for x in fit) + " | " + ", ".join(f"{x.mean():.3e}" for x in w2)
                         + " | " + ", ".join(f"{x.mean():.3e} ({x.min():.2e}–{x.max():.2e})" for x in yb) + ref + f" | {zd.mean():.3e} | {cnt} |")
    # the settled gate against the K₀-atom run of Part C, per seed (the K₀-atom records interpolated log-linearly at 3000 − t_open)
    sl = []
    for eta in (SIGMA0, np.inf):
        rs = [d for d in D if d['eta'] == eta and d.get('gate', False) == 'settled']
        if not rs or eta not in Cmap:
            continue
        c = Cmap[eta]
        ly = np.log([np.linalg.norm(v) for v in c['y']]); lf = np.log(c['fit'])
        for m in list(MS) + [M_REF]:
            ry, rf = [], []
            for d in rs:
                r = d['rec'][m]; t_o = d['t_open'][m]
                i3 = int(np.argmin(np.abs(r['s'] - 3000.0)))
                ry.append(np.linalg.norm(r['ybar'][i3]) / np.exp(np.interp(r['s'][i3] - t_o, c['s'], ly)))
                rf.append(r['fit'][i3] / np.exp(np.interp(r['s'][i3] - t_o, c['s'], lf)))
            sl.append(f"η = {'∞' if np.isinf(eta) else f'{eta:.5f}'}, m = {m}: ‖ȳ‖ {min(ry):.3f}–{max(ry):.3f}, fit error {min(rf):.3f}–{max(rf):.3f}")
    if sl:
        lines += ["", "The settled gate against the K₀-atom run of Part C: per seed, ‖ȳ‖ and the fit error at s = 3000 over the K₀-atom run's at s = 3000 − t_open (min–max over the seeds): "
                  + "; ".join(sl) + "."]
        # the state in which the settled systems open: separations (W² to the exact atoms minus the class-mean part) and the class means' offset from Θ*(q)
        sep, off = [], []
        for d in [d for d in D if d.get('gate', False) == 'settled']:
            for m in list(MS) + [M_REF]:
                r = d['rec'][m]; kb = int(np.flatnonzero(r['s'] < d['t_open'][m])[-1])
                sep.append(np.sqrt(max(float(r['W2ex'][kb] - np.sum(r['ybar'][kb] ** 2)), 0.0))); off.append(float(np.linalg.norm(r['ybar'][kb] - lin['ya'])))
        lines += ["", f"The settled systems open replicated: at the last record before the opening ({len(sep)} runs) the rms separation is {min(sep):.1e}–{max(sep):.1e} and the class "
                  f"means lie {min(off):.1e}–{max(off):.1e} from Θ*(q) (‖y*(q)‖ = {np.linalg.norm(lin['ya']):.3e}). Each system is then the K₀-atom system started from its own class "
                  "means, so the ratios above measure the invariance of the replicated configurations and the opening offset in the weakest modes, not the uniformity of (v); "
                  "Part F puts the separations in at the start."]
        # the two gates system by system at s = 3000
        cnt = dict(all3=0, total=0)
        for eta in (SIGMA0, np.inf):
            for d_s in [d for d in D if d['eta'] == eta and d.get('gate', False) == 'settled']:
                d_t = next((d for d in D if d['eta'] == eta and d.get('gate', False) is True and d['seed'] == d_s['seed']), None)
                if d_t is None:
                    continue
                for m in list(MS) + [M_REF]:
                    rs, rt = d_s['rec'][m], d_t['rec'][m]
                    i3s, i3t = int(np.argmin(np.abs(rs['s'] - 3000.0))), int(np.argmin(np.abs(rt['s'] - 3000.0)))
                    lower = [np.linalg.norm(rt['ybar'][i3t]) < np.linalg.norm(rs['ybar'][i3s]), rt['W2ex'][i3t] < rs['W2ex'][i3s], rt['fit'][i3t] < rs['fit'][i3s]]
                    cnt['total'] += 1; cnt['all3'] += int(all(lower))
        lines += ["", f"System by system at s = 3000 the transit gate ends lower than the settled one in all three of ‖ȳ‖, W² and the fit error in {cnt['all3']} of the "
                  f"{cnt['total']} systems: the comparison of the table's seed means does not hold for every system."]
    if E:
        r = E['rec']
        i3 = int(np.flatnonzero(r['s'] == 3000.0)[0])
        fa = lambda t: "a = (" + ", ".join(f"{x:+.3f}" for x in t[0]) + "), c = (" + ", ".join(f"{x:+.3f}" for x in t[1]) + "), s = (" + ", ".join(f"{x:.3f}" for x in t[2]) + ")"
        lines += ["", "## E. The ungated η = ∞ system with the largest ‖ȳ‖ at s = 3000, continued", "",
                  f"Seed {E['seed']}, m = {E['m']} (‖ȳ‖ = {E['ybar3000']:.3f} at s = 3000), re-run from its seed (fit error at s = 3000: {r['fit'][i3]:.6e}, stored {E['stored_fit3000']:.6e}) "
                  f"and continued to s = {S_E:g}.", "",
                  "| s | fit error | 𝒥 | largest spread of the class centres | largest drift |", "|---|---|---|---|---|"]
        for k in range(len(r['s'])):
            if r['s'][k] <= 3000.0 or (r['s'][k] - 3000.0) % 4000.0 == 0.0 or k == len(r['s']) - 1:
                lines.append(f"| {r['s'][k]:.0f} | {r['fit'][k]:.4e} | {r['J'][k]:.6e} | {r['spread'][k]:.4f} | {r['drift'][k]:.2e} |")
        lines += ["", f"At s = {S_E:g} the classes are collapsed onto the atoms {fa(E['atoms'])} (exact representation: {fa(E['exact'])}); the unregularised weighted loss "
                  f"there is {E['L_atoms']:.4e}, the eigenvalues of its Hessian (whitened, weighted) are " + ", ".join(f"{x:.2e}" for x in E['hess_ev'][:3]) + " (smallest three) — a plateau, not "
                  f"a minimum — and a BFGS minimisation of the unregularised loss from those atoms ends at {E['bfgs_L']:.1e} with {fa(E['bfgs_atoms'])}."]
        if 'flow' in E:
            fl = E['flow']
            lines += ["", f"The collapsed η = ∞ dynamics from the same atoms, the unregularised K₀-atom gradient flow (LSODA, rtol 10⁻¹⁰, flow time counted from the "
                      f"s = {S_E:g} state, the crossings located as events): the loss stays above 0.9 of its "
                      f"{fl['L'][0]:.4e} until {fl['plateau_until']:.4e}, halves at {fl['t_half']:.4e}, falls below 10⁻¹² at {fl['t_exact']:.4e} and "
                      f"ends at {fl['L_end']:.1e} with {fa(fl['atoms_end'])}: the same basin, crossed on the time scale of the negative curvature "
                      f"(1/{abs(E['hess_ev'][0]):.2e} = {1.0 / abs(E['hess_ev'][0]):.2e})."]
    if F:
        Cmap2 = {c['eta']: c for c in C}
        draws = sorted(set(x.get('draw', 0) for x in F))
        lines += ["", f"## F. Separated starts: the structure behind (v) beyond its ball", "",
                  f"From the lift of Θ*(q) plus zero-mean per-particle separations (rms in the G⁻¹ norm; {len(draws)} draws of their direction, the same for both rms), the "
                  f"prior at the shifted centres, ungated, {S_F:g} slow units: ‖ȳ‖ at the end over the K₀-atom run's at the same time (Part C) − 1, the mean over the draws "
                  "and their range. A zero-separation system is the K₀-atom system by the invariance of the replicated configurations, so the departure is the separations' "
                  "second-order effect, which Theorem D‴.18's facts (1)–(2) bound uniformly in m; the starts lie outside (v)'s ball (Part A), so (v) itself asserts nothing here.", "",
                  "| η | " + " | ".join(f"rms {r:g}, m = {m}" for r in F_RMS for m in F_MS) + " |", "|---|" + "---|" * (len(F_RMS) * len(F_MS))]
        growth, zs, sds = [], [], []
        for eta in (SIGMA0, np.inf):
            c = Cmap2[eta]; ref = float(np.linalg.norm(c['y'][int(np.argmin(np.abs(c['s'] - S_F)))]))
            cells, dep = [], {}
            for r in F_RMS:
                for m in F_MS:
                    v = np.array([x['ybar'] / ref - 1 for x in sorted([x for x in F if x['eta'] == eta and x['rms'] == r and x['m'] == m], key=lambda x: x.get('draw', 0))])
                    dep[(r, m)] = v
                    cells.append(f"{v.mean():+.4f} ({v.min():+.4f} to {v.max():+.4f})")
            lines.append(f"| {'∞' if np.isinf(eta) else f'{eta:.5f}'} | " + " | ".join(cells) + " |")
            growth += [float(a / b) for m in F_MS for a, b in zip(dep[(F_RMS[1], m)], dep[(F_RMS[0], m)])]
            for r in F_RMS:
                a, b = dep[(r, F_MS[0])], dep[(r, F_MS[1])]
                se = lambda v: float(np.std(v, ddof=1) / np.sqrt(len(v)))
                zs.append(float((a.mean() - b.mean()) / np.hypot(se(a), se(b)))); sds.append(float(np.std(a, ddof=1) / np.std(b, ddof=1)))
        st = [x['start'] for x in F if 'start' in x]
        lines += ["", f"From rms {F_RMS[0]:g} to {F_RMS[1]:g} the departure grows, draw by draw, by the factors {min(growth):.2f}–{max(growth):.2f} (the square of the ratio "
                  f"of the rms, {(F_RMS[1] / F_RMS[0]) ** 2:.0f}, for a second-order effect). At second order it is linear in the separations' class covariances, whose "
                  "expectation after the normalisation is the same at every m (classes of equal size, isotropic whitened directions): the prediction is a departure of the "
                  f"same expected size at both m, with a spread that shrinks as the covariances concentrate. Across the draws the standard deviation at m = {F_MS[0]} is "
                  f"{min(sds):.1f}–{max(sds):.1f} times that at m = {F_MS[1]}, and the means at the two sizes differ by {min(zs):+.1f} to {max(zs):+.1f} combined standard errors "
                  f"({len(draws)} draws; the four columns share their draws, so these are not independent tests)." + (f" The starts' norms |||·|||_η: {min(st):.3f}–{max(st):.3f}." if st else "")]
    if Gb:
        A_, P_ = Gb['a'], Gb['a']['point']
        fb = lambda B: ", ".join(f"({A:g}, {C_:g}, {S_:g})" for A, C_, S_ in B)
        eta_s = lambda e: '∞' if np.isinf(e) else f'{e:.5f}'
        pe = lambda P: "; ".join(f"η = {eta_s(d['eta'])}: the plan the classes {d['plan_is_classes']}, drift {d['dT']:.1e} (Θ) and {d['dZ']:.1e} (Z), largest real part of the joint "
                                 f"Jacobian {d['max_re']:+.2e}" for d in P['per_eta'])
        lines += ["", "## G. The limits of (i)", "",
                  f"(a) A non-replicated stationary point. The four-bump target (A, C, S) = {fb(B4)}, which six kernels of weight 1/6 fit exactly with the third class split over the "
                  f"last two bumps (fit error {A_['exact_fit']:.1e}, G⁻¹-distance of the split pair {A_['exact_split']:.4f}; K₀ = 3, q = 1/3, m = 6). The stationary point of the η = ∞ "
                  f"form (the leak toward the class barycentres) by least squares from that fit (residual {A_['res']:.1e}): the split pair at distance {P_['split']:.4f}, the leak force "
                  f"σ₀‖ϑ_i − T̂_π(i)‖ on its particles {P_['leak'][4]:.4f}, {P_['leak'][5]:.4f} (largest on the others {max(P_['leak'][:4]):.1e}), the unregularised gradient there "
                  f"{P_['unreg']:.2e} (at the exact fit {A_['exact_unreg']:.1e}), fit error {P_['fit']:.3e}. With the prior at the class barycentres it is stationary for the joint "
                  f"flow — {pe(P_)} — a saddle: the leak holds the split class against the loss.", "",
                  f"The joint flow from the exact split fit ({S_G:g} slow units): the split distance and the fit error at the records; at the end the classes' largest "
                  "spread, the merged atom (the widest) — its width coordinate κ and width s — and the largest whitened unregularised drift at the class atoms:", "",
                  "| η | split distance at s = " + ", ".join(f"{x:g}" for x in Gb['a_flow'][0]['rec']['s']) + " | fit error at the same s | at the end: spread, κ, s, drift |",
                  "|---|---|---|---|"]
        for fl in Gb['a_flow']:
            r_ = fl['rec']; sw_end = float(np.ravel(fl['atoms'][2])[fl['kk']])
            lines.append(f"| {eta_s(fl['eta'])} | " + ", ".join(f"{x:.4f}" for x in r_['split']) + " | " + ", ".join(f"{x:.3e}" for x in r_['fit'])
                         + f" | {fl['spread']:.1e}, {fl['kappa_end']:+.3f}, {sw_end:.3f}, {fl['drift_end']:.1e} |")
        lines += ["", f"The collapsed η = ∞ dynamics continued from the end state (the unregularised K₀-atom gradient flow, LSODA, rtol 10⁻¹⁰; s counted from the end of the "
                  "run): the merged atom's width coordinate κ, the gap s_max − s of its width, the fit error.", "",
                  "| η | s | κ | s_max − s | fit error |", "|---|---|---|---|---|"]
        for fl in Gb['a_flow']:
            for c_ in fl['cont']:
                lines.append(f"| {eta_s(fl['eta'])} | {c_['s']:.0e} | {c_['kappa']:+.3f} | {c_['gap']:.2e} | {c_['fit']:.6e} |")
        incs = []
        for fl in Gb['a_flow']:
            kap = [c_['kappa'] for c_ in fl['cont']]
            incs.append(f"η = {eta_s(fl['eta'])}: " + ", ".join(f"{(kap[k + 1] - kap[k]) / np.log(10):.2f}" for k in range(len(kap) - 1)))
        fc = Gb['face']
        lines += ["", "The increments of κ per decade over ln 10: " + "; ".join(incs) + " — κ grows like ln s (the logistic saturation of the width) and s tends to "
                  "s_max like 1/s: the solution leaves every compact set, as in (b). Its physical image approaches the minimiser of the loss on the face where the merged "
                  f"atom's width is s_max (BFGS over the other parameters: fit error {fc['fit']:.6e}, gradient {fc['grad']:.1e}, inside the other bounds: {fc['inside']}), "
                  f"where the loss still decreases outward (∂L/∂s = {fc['slope']:.3e}): a Karush–Kuhn–Tucker point of the box-constrained loss, not a critical point, "
                  "reached only as s → ∞.", "",
                  "The exact six-kernel fit is not stationary: the joint flow collapses the split class onto a replicated configuration with a larger fit error and then "
                  "runs off — the K₀-cluster bias of (iv)'s regulariser, and a second instance of the lost a-priori bound, at a target the network fits exactly.", "",
                  f"(b) Unbounded solutions. The target (A, C, S) = {fb(B_SAT)}: the first bump's weighted amplitude 3·{B_SAT[0][0]:g} = {3 * B_SAT[0][0]:g} exceeds a_max = "
                  f"{Gb['b']['a_max']:g}, so no exact representation with q = 1/3 exists. The K₀-atom system from the instance's exact atoms with the prior there (LSODA, rtol 10⁻¹⁰): the "
                  "first atom's amplitude coordinate ϑ_a and its prior atom's z_a (unconstrained), its amplitude a, the fit error and the largest drift.", "",
                  "| η | s | ϑ_a | z_a | a | fit error | largest drift |", "|---|---|---|---|---|---|---|"]
        for run in Gb['b']['runs']:
            for r_ in run['rows']:
                lines.append(f"| {eta_s(run['eta'])} | {r_['s']:.0e} | {r_['th_a0']:.4f} | " + ("—" if not np.isfinite(r_['z_a0']) else f"{r_['z_a0']:.4f}")
                             + f" | {r_['a0']:.6f} | {r_['fit']:.4e} | {r_['drift']:.2e} |")
        incs = []
        for run in Gb['b']['runs']:
            if run['eta'] > 0:
                th = [r_['th_a0'] for r_ in run['rows']]
                incs.append(f"η = {eta_s(run['eta'])}: " + ", ".join(f"{(th[k + 1] - th[k]) / np.log(10):.3f}" for k in range(2, len(th) - 1)))
        lines += ["", "The increments of ϑ_a per decade of s over ln 10, from s = 10³ on: " + "; ".join(incs) + " — ϑ_a grows like 2 ln s, and a approaches a_max like 1/s: under an "
                  "adaptive prior the solution leaves every compact set, while the fixed prior's converges.", "",
                  f"(c) The relabelling. The exact representation (m = 6) with the prior's atoms permuted cyclically (largest offset {Gb['c']['offset']:.3f}, unconstrained "
                  f"coordinates): the plan maps the classes to the atoms {Gb['c']['classes']}, W² = {Gb['c']['W2']:.1e}, drift {Gb['c']['dT']:.1e} (Θ) and {Gb['c']['dZ']:.1e} (Z) — "
                  "stationary, with z_σ(k) = ϑ_k for the plan's bijection σ."]
    if H:
        eta_s = lambda e: '∞' if np.isinf(e) else f'{e:.5f}'
        lo, hi = H['box']
        dc = H['dt_check']
        lines += ["", "## H. Proposition D‴.27: the projected adaptive prior", "",
                  f"The box 𝒵 for the prior's atoms (unconstrained coordinates): α ∈ [{lo[0]:.4f}, {hi[0]:.4f}], γ ∈ [{lo[1]:.4f}, {hi[1]:.4f}], κ ∈ [{lo[2]:g}, {hi[2]:g}] "
                  f"(amplitudes and centres to {A_BOUND:g} of a_max and c_max). (a) Part G(b)'s target (the first bump's weighted amplitude 4.8 > a_max), the K₀-atom system "
                  f"from the instance's exact atoms with the prior there, the prior projected onto 𝒵 (η = σ₀: the projected midpoint rule — positions projected inside the "
                  f"step — at dt = {dc['dts'][0]:g} (the rows) and {dc['dts'][1]:g}, to s = {S_H:.0e}; η = ∞: LSODA, rtol 10⁻¹⁰, the prior the projected class mean): the first "
                  "atom's amplitude coordinate ϑ_a, its prior atom's z_a, the amplitude a, the fit error and the largest drift (of the projected flow).", "",
                  "| η | s | ϑ_a | z_a | a | fit error | largest drift |", "|---|---|---|---|---|---|---|"]
        for run in H['runs']:
            for r_ in run['rows']:
                lines.append(f"| {eta_s(run['eta'])} | {r_['s']:.0e} | {r_['th_a0']:.4f} | {r_['z_a0']:.4f} | {r_['a0']:.6f} | {r_['fit']:.4e} | {r_['drift']:.2e} |")
        lines += ["", f"At dt = {dc['dts'][1]:g} the η = σ₀ run ends within {dc['diff']:.1e} of the dt = {dc['dts'][0]:g} run (atoms and prior), with the projected flow's "
                  f"drift {dc['drifts'][1]:.1e} at s = {S_H:.0e} ({dc['drifts'][0]:.1e} at dt = {dc['dts'][0]:g}) and the prior's amplitude atom on the face (gaps "
                  + ", ".join(f"{g:.1e}" for g in dc['gaps']) + ")" + (": the shared fixed point. The scheme's consistency: at s = " + f"{S_TRAJ:g} the "
                  "trajectory (atoms and prior) differs from a reference at dt = " + f"{DT_REF:g} by " + ", ".join(f"{v:.1e} (dt = {k:g})" for k, v in dc['traj'].items())
                  + "." if 'traj' in dc else "."),
                  "", "The largest G⁻¹-distance of a particle from 𝒵 along the run: " + "; ".join(f"η = {eta_s(r['eta'])}: {r['max_dist']:.4f}" for r in H['runs'])
                  + f" (η = σ₀ at dt = {dc['dts'][1]:g}: {dc['max_dist2']:.4f}), against the a-priori bound B/σ₀ ≤ {H['B'] / SIGMA0:.1f} (B ≤ {H['B']:.2f} = (max|f| + a_max)·L, "
                  f"L = {H['LG']:.3f} the closed-form bound of L_ϑ^G; 4.044 by search, `learning_table.md`). Both priors converge to the same bounded stationary point, with "
                  "the prior's amplitude atom on the boundary of 𝒵; in this run the fit error lies between the fixed prior's and the unprojected adaptive prior's limit "
                  "(Part G(b)) — (c) shows its dependence on the box.", "",
                  f"(b) Part C's K₀-atom loop at the misplaced prior (η = σ₀, {S_C:g} slow units) with and without 𝒵: the largest difference of the atoms and of the prior "
                  f"at the end is {H['diff_b']:.1e} (the run stays inside 𝒵: {H['inside_b']}). The box is never approached — along the run the prior keeps a margin of at least "
                  f"{H['margin_b']['Z']:.1%} of the box's width to every face, the particles {H['margin_b']['Th']:.1%} — so the two runs coincide by construction (the "
                  "projection is the identity inside 𝒵).", "",
                  "(c) The box's size: (a)'s η = ∞ system (LSODA to s = " + f"{T_G[-1]:.0e}" + ") with the amplitude bound of 𝒵 varied (the centre and width bounds as above).", "",
                  "| amplitude bound (× a_max) | the α-face | the start (the exact atoms) inside 𝒵 | ϑ_a | z_a | fit error | largest drift |", "|---|---|---|---|---|---|---|"]
        for r_ in H['scan']:
            lines.append(f"| {r_['a_bound']:g} | {r_['face']:.4f} | {'yes' if r_.get('start_inside', True) else 'no'} | {r_['th_a0']:.4f} | {r_['z_a0']:.4f} | {r_['fit']:.4e} | {r_['drift']:.1e} |")
        if Gb:
            fx = [r for r in Gb['b']['runs'] if r['eta'] == 0.0][0]['rows'][-1]['fit']
            un = [r for r in Gb['b']['runs'] if np.isinf(r['eta'])][0]['rows'][-1]['fit']
            btw = [r_['a_bound'] for r_ in H['scan'] if un <= r_['fit'] <= fx]
            worse = [(r_['a_bound'], r_['fit']) for r_ in H['scan'] if r_['fit'] > fx]
            outs = [r_['a_bound'] for r_ in H['scan'] if not r_.get('start_inside', True)]
            lines += ["", f"The fixed prior's limit is {fx:.4e} and the unprojected adaptive prior's fit error at s = {T_G[-1]:.0e} {un:.4e} (Part G(b)): for the amplitude "
                      f"bounds {min(btw):g}–{max(btw):g} the box's fit error lies between them, decreasing as the box grows"
                      + (f"; the tighter boxes are worse than the fixed prior (" + ", ".join(f"{f_:.4f} at {a_:g}" for a_, f_ in worse) + ")" if worse else "")
                      + (f", and the boxes with the bound ≤ {max(outs):g} exclude the instance's exact atoms (the start is projected at once): a box should contain "
                         "the representation." if outs else ".")]
        if 'sel' in H:
            sl = H['sel']
            i_ = int(np.argmax(np.diff(sl['red'])))
            lines += ["", f"(e) At η = ∞ the plan is Lloyd-consistent — optimal for its own class means — not jointly optimal: a configuration of six particles "
                      f"(zero amplitudes, the target 0.02 × the three bumps, the box; found by a seeded search, trial {int(sl['trial'])}) with the prior at the class means "
                      f"of the consistent pairing {[tuple(int(x) for x in pp) for pp in sl['pairing']]}. Along AdaptivePriorLoop's run (plan 'ot', η = ∞, 400 slow units) "
                      f"the selected energy 𝒥(Θ, Z) falls ({sl['sel'][0]:.4e} → {sl['sel'][-1]:.4e}, its largest increment {np.diff(sl['sel']).max():+.1e}), while the "
                      f"reduced energy min_Z𝒥 (the minimum over the 15 balanced pairings) rises from {sl['red'][0]:.4e} to {sl['red'].max():.4e} (its largest increment "
                      f"{np.diff(sl['red']).max():+.2e}, at s = {sl['s'][i_]:g}–{sl['s'][i_ + 1]:g}) before the plan becomes optimal (at s = "
                      f"{sl['s'][int(np.argmax(sl['optimal']))]:g}): the selected energy is the Lyapunov function, the reduced energy is not."]
        g4 = H['G4']
        lines += ["", f"(d) Part G(a)'s four-bump target with 𝒵 at η = ∞: the joint flow from the exact split fit ({S_G:g} slow units; plan classes {g4['classes']}, largest "
                  f"spread {g4['spread']:.1e}, the merged atom's κ = {g4['kappa_start']:.3f}, fit error {g4['fit_start']:.6e}), then the collapsed dynamics with 𝒵 (the "
                  "leak toward the projection of the atoms; LSODA, rtol 10⁻¹⁰; s counted from the end of the run): the merged atom's width coordinate κ, its prior atom's, "
                  "the fit error and the largest whitened drift.", "",
                  "| s | κ | the prior's κ | fit error | largest drift |", "|---|---|---|---|---|"]
        for r_ in g4['cont']:
            lines.append(f"| {r_['s']:.0e} | {r_['kappa']:+.6f} | {r_['kappa_prior']:+.4f} | {r_['fit']:.6e} | {r_['drift']:.1e} |")
        lines += ["", f"Without 𝒵 (Part G(a)) κ grows like ln s and the solution leaves every compact set; with 𝒵 the merged atom's prior stops on the face κ = {g4['face']:g} and "
                  f"the atom stays within {abs(g4['cont'][-1]['kappa'] - g4['face']):.1e} of it: a bounded stationary point, its fit error within "
                  + (f"{abs(g4['cont'][-1]['fit'] - Gb['face']['fit']):.1e} of the face minimiser's of Part G(a)." if Gb else "that of Part G(a)'s face minimiser.")]
    lines += ["", "Numbers are produced by `experiments/exp_adaptive_prior.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'adaptive_prior.md'), 'w', encoding='utf-8') as f:
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
    raw = os.path.join(OUT, 'adaptive_prior_raw.npz')
    if '--from-raw' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        report(z['lin'].item(), list(z['C']), list(z['D']), z['E'].item() if 'E' in z.files else None, t0, list(z['F']) if 'F' in z.files else None,
               z['G'].item() if 'G' in z.files else None, z['H'].item() if 'H' in z.files else None)
        return
    if '--only-H' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        save['H'] = np.array(projected_job(), dtype=object)
        save_raw(raw, **save)
        report(save['lin'].item(), list(save['C']), list(save['D']), save['E'].item() if 'E' in save else None, t0, list(save['F']) if 'F' in save else None,
               save['G'].item() if 'G' in save else None, save['H'].item())
        return
    if '--only-G' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        save['G'] = np.array(boundary_job(), dtype=object)
        save_raw(raw, **save)
        report(save['lin'].item(), list(save['C']), list(save['D']), save['E'].item() if 'E' in save else None, t0, list(save['F']) if 'F' in save else None, save['G'].item(),
               save['H'].item() if 'H' in save else None)
        return
    if '--only-F' in sys.argv or '--only-lin' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        if '--only-lin' in sys.argv:
            save['lin'] = np.array(linear_parts(), dtype=object)
        if '--only-F' in sys.argv:
            with mp.Pool(N_WORKERS) as pool:
                Fr = pool.map(sep_job, F_JOBS, chunksize=1)
            save['F'] = np.array(Fr, dtype=object)
        save_raw(raw, **save)
        report(save['lin'].item(), list(save['C']), list(save['D']), save['E'].item() if 'E' in save else None, t0, list(save['F']) if 'F' in save else None,
               save['G'].item() if 'G' in save else None, save['H'].item() if 'H' in save else None)
        return
    if '--only-E' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        E = trapped_job(list(save['D']))
        save['E'] = np.array(E, dtype=object)
        save_raw(raw, **save)
        report(save['lin'].item(), list(save['C']), list(save['D']), E, t0, list(save['F']) if 'F' in save else None, save['G'].item() if 'G' in save else None, save['H'].item() if 'H' in save else None)
        return
    if '--only-gated' in sys.argv:
        which = sys.argv[sys.argv.index('--only-gated') + 1] if len(sys.argv) > sys.argv.index('--only-gated') + 1 else 'transit'
        gate = True if which == 'transit' else which
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        D = [d for d in save['D'] if d.get('gate', False) != gate]
        with mp.Pool(N_WORKERS) as pool:
            D += pool.map(m_job, [(seed, eta, gate) for seed in range(N_SEEDS) for eta in ETAS_D if eta > 0])
        save['D'] = np.array(D, dtype=object)
        save_raw(raw, **save)
        report(save['lin'].item(), list(save['C']), D, save['E'].item() if 'E' in save else None, t0, list(save['F']) if 'F' in save else None,
               save['G'].item() if 'G' in save else None, save['H'].item() if 'H' in save else None)
        return
    lin = linear_parts()
    with mp.Pool(N_WORKERS) as pool:
        C = pool.map(k0_job, list(ETAS))
        D = pool.map(m_job, [(seed, eta, False) for seed in range(N_SEEDS) for eta in ETAS_D] + [(seed, eta, g) for g in (True, 'settled') for seed in range(N_SEEDS) for eta in ETAS_D if eta > 0])
        F = pool.map(sep_job, F_JOBS, chunksize=1)
    E = trapped_job(D)
    Gb = boundary_job()
    H = projected_job()
    save_raw(raw, lin=np.array(lin, dtype=object), C=np.array(C, dtype=object), D=np.array(D, dtype=object), E=np.array(E, dtype=object), F=np.array(F, dtype=object),
             G=np.array(Gb, dtype=object), H=np.array(H, dtype=object))
    report(lin, C, D, E, t0, F, Gb, H)


if __name__ == '__main__':
    main()
