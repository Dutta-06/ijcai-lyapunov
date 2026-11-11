"""
An application: cogging-torque compensation in a position servo — kernel adaptive control with transport leaks
=================================================================================================================

The servo ẍ = −τ_cog(x) + u (unit inertia, viscous friction compensated), τ_cog(x) = 0.4 sin(3x + 0.3) + 0.1 sin 7x unknown, tracks
x_m = sin t under the sliding controller u = ẍ_m − λė − k_s s − f̂_Θ(x), s := ė + λe (src/core/sliding.py: the error channel is
ṡ = −k_s s + r_Θ(x), so the averaged slow loop is Theorem J's along the position orbit x_m(t) ∈ [−1, 1]). The network: K₀ = 4 saturated
Gaussian kernels on the position (a_max = 4, c_max = 1.6, widths in [0.05, 1]); τ_cog has no exact four-kernel representation. ϑ^b(q)
(q = 1/4) is the nondegenerate local minimum of the unregularised weighted loss that BFGS reaches from START (a 1.5% fit error), where
Proposition D‴.28's analysis applies; the loss's infimum lies at the saturation (Part A's global picture), so globally the servo is
Proposition D‴.27's case and Parts B–C are statements about the basin of ϑ^b(q) and about an engineer's model.

A. No dynamics: ϑ^b(q), its fit error, the whitened unregularised Hessian Ã (its spectrum: the weakest modes), the residual curvature's
   λ⁺_k and the gauge set W_g; the fixed prior's floor law (D‴.28(ii)); the collective rates — the fixed prior's ã_j + σ₀ against the
   adaptive prior's λ_j^− (η = σ₀, ∞) — and the misplaced prior's bias (an engineer's cogging model with the phase 0.5 in place of 0.3):
   the fixed prior's stationary point and its fit error. The global picture: N_GLOBAL BFGS runs from random starts (their end fit
   errors, the saturated amplitudes) and, from the best end point, the saturated amplitudes held at |α| = ALPHA_HOLD with the rest
   re-optimised — the infimum approached at the saturation.
B. The averaged loop, m = 16, 64, 128 against M = 256, three seeds, S_B slow units: the nearest-atom leak toward the fixed and the adaptive
   prior (η = ∞ and σ₀), the OT leak toward the fixed prior, the classical σ-modification toward ϑ⁰ = 0 at the same σ₀, and no leak
   (σ₀ = 0); from warm clouds (width 0.02) around ϑ^b(q) — the oracle start, computed from the true torque — with the matching prior
   ϑ^b(q) or the misplaced prior (at η = ∞ the prior's initial position enters only through the first cells, so those runs are the
   matching ones), and from warm clouds around the engineer's model with the prior at it (the start an engineer has); the fit error (the
   bias) and ‖f̂_m − f̂_M‖² (the finite-size floor). And the collapsed K₀-atom flows (one particle per class) toward the misplaced prior to
   s = 10⁵, from both starts: the time scale of the adaptive prior's bias removal.
C. The real servo (k_s = 25, λ = 5, RK4 dt = 2·10⁻³, T_C time units, m = 16 and 64, seed 0): from the end states of the oracle-start
   designs, with the design's leak (the finite-η prior integrated in real time, dz/dt = −(η/k_s)(z − T̂)), the position-tracking error's
   rms over five periods (records every 0.02 time units) against the averaged predictions — the residual r(x_m) filtered by
   (d/dt + λ)⁻¹(d/dt + k_s)⁻¹ (exact) and by (d/dt + λ)⁻¹ after the slaved s = r/k_s.
D. Phase-plane kernels (Proposition D‴.29(iv)): the whitened Gauss–Newton spectrum of isotropic kernels along the circle, along ellipses
   and along a two-harmonic orbit, and with the kernels centred at the origin — the number of modes below 10⁻⁶σ₀.

Output: results/servo/servo_cogging.md and servo_cogging_raw.npz (SV_WORKERS, SV_SB, SV_TC; `--from-raw` regenerates the md; `--only-A`
recomputes Parts A and D; `--add-B` runs the Part B jobs and long flows not yet stored; `--only-C` reruns Part C from the stored Part B,
`--only-C <design> …` only the designs named).
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

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import chaos_learning as CL
from src.core import ot_leak as OT
from src.core import prior_design as PD
from src.core import adaptive_prior as AP
from src.core import floor_bias as FB
from src.core.ot_leak import _cost_matrix

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'servo')
NET_KW = dict(a_max=4.0, c_max=1.6, s_min=0.05, s_max=1.0)
G = np.array([1.0, 0.25, 0.25])
PER = 2 * np.pi
N_Q = 400
K0 = 4
Q = np.ones(K0) / K0
SIGMA0 = 0.01
PHASE, PHASE_PRIOR = 0.3, 0.5
MS, M_REF = (16, 64, 128), 256
JITTER = 0.02
DT, PLAN_EVERY, REC = 0.25, 4, 200
N_GLOBAL = 120
ALPHA_HOLD = (6.0, 8.0, 10.0, 14.0)
S_LONG = (600.0, 2000.0, 1e4, 1e5)
REC_C = 10                                                                 # Part C's records: every REC_C real-time steps
S_B = float(os.environ.get('SV_SB', 2000.0))
N_SEEDS = 3
N_WORKERS = int(os.environ.get('SV_WORKERS', 2))
KS, LAM, DT_REAL = 25.0, 5.0, 2e-3
T_C = float(os.environ.get('SV_TC', 1000.0))


def tau_cog(X, phase=PHASE):
    X = np.asarray(X, float)
    return 0.4 * np.sin(3 * X + phase) + 0.1 * np.sin(7 * X)


def f_vec(X, phase=PHASE):
    return -tau_cog(X, phase)


def x_ref(t):
    return np.array([np.sin(t)])


def net_():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


START = (np.array([1.0, 1.4, -1.9, -0.2]), np.array([-0.9, -0.33, 0.3, 0.78]), np.array([0.25, 0.25, 0.27, 0.1]))


def approximant(phase):
    """The best weighted K₀-kernel approximant at q of −τ_cog with the given phase: BFGS on the unregularised weighted loss from START,
    polished by Newton; returns (Θ, gradient residual, relative fit error)."""
    net = net_()
    th0 = net.unconstrained(START[0][:, None], START[1][:, None], START[2])
    loop = OT.OTPriorLoop(net, lambda X: f_vec(X, phase), x_ref, PER, G, SIGMA0, th0, n_quad=N_Q)

    def fun(v):
        T = v.reshape(K0, 3); a, c, s, d, ph = loop.features(T); r = loop.F - ph @ (Q * a)
        return 0.5 * float(np.mean(r ** 2)), (Q[:, None] * loop.grad_V_L(r, T)).ravel()
    res = minimize(fun, th0.ravel(), jac=True, method='BFGS', options=dict(gtol=1e-14, maxiter=40000))
    th, rr, fit = FB.exact_weighted_atoms(loop, Q, res.x.reshape(K0, 3), n_cont=1, w0=Q, fit_tol=1.0)
    return th, rr, fit


def setup(prior_phase=PHASE, sigma0=SIGMA0, zero_prior=False):
    """The loop for the true target with the prior at the approximant of the cogging model with prior_phase (the matching prior when it
    equals PHASE), or at ϑ⁰ = 0 (σ-modification)."""
    net = net_()
    tb = approximant(PHASE)[0]
    prior = np.zeros((1, 3)) if zero_prior else approximant(prior_phase)[0]
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, sigma0, prior, n_quad=N_Q)
    return net, tb, prior, loop


# ---------------------------------------------------------------------------- A
def part_A():
    net, tb, prior, loop = setup(PHASE)
    th, rr, fit = approximant(PHASE)
    gradL = lambda v: (Q[:, None] * FB.unregularised_grad(loop, Q, v.reshape(K0, 3))).ravel()
    A = np.zeros((3 * K0, 3 * K0)); e = 1e-6
    for j in range(3 * K0):
        E = np.zeros(3 * K0); E[j] = e
        A[:, j] = (gradL(tb.ravel() + E) - gradL(tb.ravel() - E)) / (2 * e)
    Dh = np.kron(np.sqrt(Q), np.sqrt(1.0 / G))
    At = 0.5 * (A + A.T) / Dh[:, None] / Dh[None, :]
    ev = np.linalg.eigvalsh(At)
    a, c, s, d, ph = loop.features(tb)
    r = loop.F - ph @ (Q * a)
    lam = loop.kappa_R(r, tb)
    lo = Q * np.abs(a) / net.a_max
    rates = {eta: np.sort(AP.pair_rates(ev, SIGMA0, eta)[0]) for eta in (SIGMA0, np.inf)}
    # the fixed matching prior's floor law: σ₀Ĵ(Ã + σ₀I)⁻¹g̃, the gauge displacement along e_j − e_K
    Jh = FB.whitened_jacobian(loop, Q, tb)
    h = 1e-6
    Vd = np.array([np.eye(K0)[j] - np.eye(K0)[K0 - 1] for j in range(K0 - 1)]).T
    def gauge(w):
        a_, c_, s_ = net.physical(tb)
        return net.unconstrained(Q[:, None] * a_ / np.asarray(w, float)[:, None], c_, s_)
    gt = np.column_stack([-Dh * ((gauge(Q + h * Vd[:, k]) - gauge(Q - h * Vd[:, k])) / (2 * h)).ravel() for k in range(K0 - 1)])
    F = SIGMA0 * Jh @ np.linalg.solve(At + SIGMA0 * np.eye(3 * K0), gt)
    Ct = (np.diag(Q) - np.outer(Q, Q))[:K0 - 1, :K0 - 1]
    law = float(np.trace(F @ Ct @ F.T)) / loop.N
    # the misplaced prior: the fixed prior's stationary point and its fit error
    netp, tbp, priorp, loopp = setup(PHASE_PRIOR)
    Tw, resw = PD.weighted_stationary_atoms_prior(loopp.nearest, Q, priorp)
    fitw = float(np.sqrt(np.mean((loopp.F - loopp.features(Tw)[4] @ (Q * loopp.features(Tw)[0])) ** 2) / np.mean(loopp.F ** 2)))
    fit_prior = float(np.sqrt(np.mean((loopp.F - loopp.features(priorp)[4] @ (Q * loopp.features(priorp)[0])) ** 2) / np.mean(loopp.F ** 2)))
    # the same replicated point from the approximant, and its D‴.5 splitting margin (the fixed prior's replicated point a saddle if negative)
    Tw2, resw2 = PD.weighted_stationary_atoms_prior(loopp.nearest, Q, tbp)
    marg_mis, _, kap_mis = PD.splitting_margin_prior(loopp.nearest, Q, priorp)
    marg_match = float(SIGMA0 - np.max(lam))
    return dict(fit=fit, res=rr, ev=ev, lam=lam, lo=lo, rates=rates, law=law, phys=tuple(x for x in net.physical(tb)), fit_misplaced=fitw,
                res_misplaced=resw, fit_prior=fit_prior, same_point=float(np.abs(Tw2 - Tw).max()), margin_misplaced=marg_mis,
                kappa_misplaced=np.asarray(kap_mis), margin_matching=marg_match)


def perm_distance(T, T0):
    """The largest coordinate difference between two K₀-atom configurations, minimised over the relabellings of the kernels."""
    from itertools import permutations
    return min(float(np.abs(np.asarray(T)[list(pp)] - T0).max()) for pp in permutations(range(T0.shape[0])))


def part_A_global():
    """The global picture of L_q: BFGS from N_GLOBAL random starts (physical amplitudes, centres and widths uniform inside their boxes),
    the end points' fit errors, their largest |α| and their distance from ϑ^b(q) up to relabelling; from the best end point, the saturated
    amplitudes (|a| > 0.999a_max) held at ±ALPHA_HOLD and the rest re-optimised."""
    net = net_()
    tb = approximant(PHASE)[0]
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, tb, n_quad=N_Q)

    def fun(v, fixed=None, vals=None):
        T = v.reshape(K0, 3).copy()
        if fixed is not None:
            T[fixed, 0] = vals
        a, c, s, d, ph = loop.features(T); r = loop.F - ph @ (Q * a)
        g = (Q[:, None] * loop.grad_V_L(r, T))
        if fixed is not None:
            g[fixed, 0] = 0.0
        return 0.5 * float(np.mean(r ** 2)), g.ravel()
    rng = np.random.default_rng(11)
    ends = []
    for k in range(N_GLOBAL):
        a0 = rng.uniform(-0.9, 0.9, K0) * net.a_max; c0 = rng.uniform(-0.9, 0.9, K0) * net.c_max
        s0 = rng.uniform(net.s_min + 0.02, net.s_max - 0.02, K0)
        th0 = net.unconstrained(a0[:, None], c0[:, None], s0)
        res = minimize(fun, th0.ravel(), jac=True, method='BFGS', options=dict(gtol=1e-12, maxiter=20000))
        T = res.x.reshape(K0, 3)
        ends.append(dict(fit=float(loop.fit_error(T)), amax=float(np.abs(T[:, 0]).max()), dist=perm_distance(T, tb), T=T))
    best = min(ends, key=lambda e: e['fit'])
    sat = np.where(np.abs(np.tanh(best['T'][:, 0] / net.a_max)) > 0.999)[0]
    hold = []
    for al in ALPHA_HOLD:
        vals = np.sign(best['T'][sat, 0]) * al
        res = minimize(lambda v: fun(v, sat, vals), best['T'].ravel(), jac=True, method='BFGS', options=dict(gtol=1e-12, maxiter=20000))
        T = res.x.reshape(K0, 3).copy(); T[sat, 0] = vals
        hold.append(dict(alpha=al, a=float(net.a_max * np.tanh(al / net.a_max)), fit=float(loop.fit_error(T))))
    a, c, s_ = net.physical(best['T'])
    print("[A global]", flush=True)
    return dict(fits=np.array([e['fit'] for e in ends]), amax=np.array([e['amax'] for e in ends]), dist=np.array([e['dist'] for e in ends]),
                best_phys=(a[:, 0], c[:, 0], s_), best_fit=best['fit'], n_sat=len(sat), hold=hold, fit_b=float(loop.fit_error(tb)))


def part_D():
    """Proposition D‴.29(iv): the whitened Gauss–Newton spectrum of K = 5 isotropic phase-plane kernels (tests/test_sliding.py's atoms) along
    the circle (1.2 sin t, 1.2 cos t), the ellipses with the velocity axis scaled by 1.1 and 2, and the orbit of x_m = 1.2 sin t + 0.3 sin 2t,
    and with the five centres at the origin on the circle; the counts below 10⁻⁶σ₀ (σ₀ = 0.01) and 10⁻¹⁵ of the largest."""
    from src.core import sliding as SL
    netp = SL.PhasePlaneNetwork(2, a_max=4.0, c_max=(3.0, 3.0), s_min=0.15, s_max=2.0)
    Gp = np.array([1.0, 0.25, 0.25, 0.25]); sig = 0.01
    rng = np.random.default_rng(1)
    K = 5
    atoms = netp.unconstrained(rng.normal(0, 1, K), rng.uniform(-1, 1, (K, 2)), rng.uniform(0.4, 1.0, K))
    centred = netp.unconstrained(rng.normal(0, 1, K), np.zeros((K, 2)), rng.uniform(0.4, 1.0, K))
    orbits = [('circle', lambda t: 1.2 * np.sin(t), lambda t: 1.2 * np.cos(t), atoms),
              ('ellipse, velocity axis ×1.1', lambda t: 1.2 * np.sin(t), lambda t: 1.32 * np.cos(t), atoms),
              ('ellipse, velocity axis ×2', lambda t: 1.2 * np.sin(t), lambda t: 2.4 * np.cos(t), atoms),
              ('x_m = 1.2 sin t + 0.3 sin 2t', lambda t: 1.2 * np.sin(t) + 0.3 * np.sin(2 * t), lambda t: 1.2 * np.cos(t) + 0.6 * np.cos(2 * t), atoms),
              ('circle, the centres at the origin', lambda t: 1.2 * np.sin(t), lambda t: 1.2 * np.cos(t), centred)]
    rows = []
    for name, xr, xrd, at in orbits:
        lp = SL.AveragedSlidingLoop(netp, lambda Z: -np.sin(Z[:, 0]), xr, xrd, PER, Gp, sig, n_quad=400)
        H, ev, V = lp.gauss_newton(np.ones(K) / K, at)
        rows.append(dict(name=name, n_small=int(np.sum(ev < 1e-6 * sig)), n_zero=int(np.sum(ev < 1e-15 * ev[-1])), smallest=float(ev[0] / sig),
                         next_rel=float(ev[K] / ev[-1]), dim=len(ev)))
    print("[D]", flush=True)
    return rows


# ---------------------------------------------------------------------------- B
DESIGNS = ('nearest-fixed', 'nearest-adaptive', 'nearest-adaptive-slow', 'ot-fixed', 'sigma-mod', 'no-leak')


def design_loop(design, prior_phase):
    zero = design == 'sigma-mod'
    net, tb, prior, loop = setup(prior_phase, zero_prior=zero, sigma0=0.0 if design == 'no-leak' else SIGMA0)
    eta = np.inf if design == 'nearest-adaptive' else (SIGMA0 if design == 'nearest-adaptive-slow' else 0.0)
    plan = 'ot' if design == 'ot-fixed' else 'nearest'
    return net, tb, prior, loop, eta, plan


def start_atoms(start):
    """The centre of the warm clouds: ϑ^b(q) (the oracle start, from the true torque) or the engineer's model (the approximant of the
    cogging model with the phase PHASE_PRIOR)."""
    return approximant(PHASE)[0] if start == 'oracle' else approximant(PHASE_PRIOR)[0]


def run_job(args):
    design, prior_phase, seed = args[:3]
    start = args[3] if len(args) > 3 else 'oracle'
    t0 = time.time()
    net, tb, prior, loop, eta, plan = design_loop(design, prior_phase)
    clouds = CL.mixture_start(start_atoms(start), M_REF, np.random.default_rng(seed), JITTER, G)
    out = dict(design=design, prior_phase=prior_phase, seed=seed, start=start, rec={}, final={})
    for m in MS + (M_REF,):
        apl = AP.AdaptivePriorLoop(loop, eta, plan=plan)
        rec = dict(s=[], fit=[], fh=[])
        Z0 = prior.copy()

        def cb(s, Th, Z, adapting, rec=rec):
            rec['s'].append(s); rec['fit'].append(loop.fit_error(Th)); rec['fh'].append(loop.f_hat(Th).astype(np.float32))
        Th, Z, _ = apl.run(clouds[:m].copy(), Z0, S_B, DT, plan_every=PLAN_EVERY, rec_every=REC, callback=cb, gate=False)
        out['rec'][m] = {k: np.array(v) for k, v in rec.items()}
        out['final'][m] = (Th.astype(np.float32), np.asarray(Z, np.float32))
    print(f"[{time.time() - t0:.0f}s] B {design} prior phase {prior_phase} seed {seed} start {start}", flush=True)
    return out


def long_job(args):
    """A collapsed K₀-atom flow (one particle per class, the weights 1/K₀) toward the misplaced prior, from ϑ^b(q) (oracle) or from the
    engineer's model, to s = S_LONG[-1]: the fit error and the distance from ϑ^b(q) (up to relabelling) at S_LONG."""
    design, start = args
    t0 = time.time()
    net, tb, prior, loop, eta, plan = design_loop(design, PHASE_PRIOR)
    apl = AP.AdaptivePriorLoop(loop, eta, plan=plan)
    Th = start_atoms(start).copy(); Z = prior.copy()
    out = dict(design=design, start=start, s=[], fit=[], dist=[], amax=[])
    s_done = 0.0
    for S in S_LONG:
        Th, Z, _ = apl.run(Th, Z, S - s_done, DT, plan_every=PLAN_EVERY, rec_every=10 ** 9, gate=False)
        s_done = S
        out['s'].append(S); out['fit'].append(float(loop.fit_error(Th))); out['dist'].append(perm_distance(Th, tb))
        out['amax'].append(float(np.abs(Th[:, 0]).max()))
    print(f"[{time.time() - t0:.0f}s] long {design} {start}", flush=True)
    return out


# ---------------------------------------------------------------------------- C
def real_job(args):
    """The real servo from the end state of a Part B design (seed 0), the design's leak (a finite-η prior relaxed toward the class means
    at the rate η/k_s in real time), T_C time units: e² every REC_C steps (0.02 time units)."""
    design, prior_phase, m, Th0, Z0 = args
    t0 = time.time()
    net, tb, prior, loop, eta, plan = design_loop(design, prior_phase)
    Th = np.asarray(Th0, float).copy(); Z = np.asarray(Z0, float).copy()
    x, v, t = 0.0, 1.0, 0.0                                                # on the reference: x_m(0) = 0, ẋ_m(0) = 1
    sig = loop.sigma0 / KS                                                 # the design's σ₀ (0 without a leak)

    def targets(Th, Z):
        if design == 'ot-fixed':
            return loop.leak_target(Th)
        lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
        if np.isinf(eta):
            Z = np.array([Th[lab == k].mean(axis=0) if np.any(lab == k) else Z[k] for k in range(Z.shape[0])])
        return Z[lab]

    def relax(Th, Z, dtau):
        """The finite-η prior over dtau time units toward the class means of its nearest-atom cells (exact relaxation, cells frozen)."""
        if eta == 0 or np.isinf(eta) or design == 'ot-fixed':
            return Z
        lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
        cm = np.array([Th[lab == k].mean(axis=0) if np.any(lab == k) else Z[k] for k in range(Z.shape[0])])
        return cm + (Z - cm) * np.exp(-(eta / KS) * dtau)

    def F(t, x, v, Th, T):
        e = x - np.sin(t); ed = v - np.cos(t); s = ed + LAM * e
        fh = float(net.f_hat(np.array([x]), Th)[0])
        u = -np.sin(t) - LAM * ed - KS * s - fh
        vd = float(f_vec(np.array([x]))[0]) + u
        return v, vd, G[None, :] * net.grad_g(np.array([x]), np.array([s]), Th) - sig * (Th - T)
    n = int(round(T_C / DT_REAL)); ts, e2 = [], []
    marks = {int(round(tt / DT_REAL)): tt for tt in (5 * PER, 500.0, T_C)}
    snaps = {}
    T = targets(Th, Z)
    for i in range(n + 1):
        if i in marks:
            snaps[marks[i]] = Th.copy()
        if i == n:
            break
        if i % REC_C == 0:
            if i > 0:
                Z = relax(Th, Z, REC_C * DT_REAL)
            T = targets(Th, Z)
            ts.append(t); e2.append((x - np.sin(t)) ** 2)
        k1 = F(t, x, v, Th, T); k2 = F(t + DT_REAL / 2, x + DT_REAL / 2 * k1[0], v + DT_REAL / 2 * k1[1], Th + DT_REAL / 2 * k1[2], T)
        k3 = F(t + DT_REAL / 2, x + DT_REAL / 2 * k2[0], v + DT_REAL / 2 * k2[1], Th + DT_REAL / 2 * k2[2], T)
        k4 = F(t + DT_REAL, x + DT_REAL * k3[0], v + DT_REAL * k3[1], Th + DT_REAL * k3[2], T)
        x += DT_REAL / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); v += DT_REAL / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        Th = Th + DT_REAL / 6 * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]); t += DT_REAL
    # the averaged prediction at the parameters reached at each window's end: s = r(x_m)/k_s, e = (d/dt + λ)⁻¹s along the orbit
    from src.core.sliding import filtered_error
    pred, pred_ex = {}, {}
    for tt, Ts in snaps.items():
        s_pred, e_pred = filtered_error(loop.F - loop.f_hat(Ts), PER, KS, LAM)
        pred[tt] = float(np.sqrt(np.mean(e_pred ** 2)))
        pred_ex[tt] = float(np.sqrt(np.mean(filtered_error(loop.F - loop.f_hat(Ts), PER, KS, LAM, exact=True)[1] ** 2)))
    print(f"[{time.time() - t0:.0f}s] C {design} prior phase {prior_phase} m {m}", flush=True)
    return dict(design=design, prior_phase=prior_phase, m=m, ts=np.array(ts), e2=np.array(e2), e_pred=pred, e_pred_exact=pred_ex, final=Th.astype(np.float32),
                Zfinal=np.asarray(Z, np.float32), Z0=np.asarray(Z0, np.float32))


# ---------------------------------------------------------------------------- report
NAMES = {'nearest-fixed': 'nearest atom, fixed prior', 'nearest-adaptive': 'nearest atom, adaptive prior (η = ∞)',
         'nearest-adaptive-slow': 'nearest atom, adaptive prior (η = σ₀)', 'ot-fixed': 'OT, fixed prior',
         'sigma-mod': 'σ-modification (toward ϑ⁰ = 0, the same σ₀)', 'no-leak': 'no leak (σ₀ = 0)'}


def class_diameter(Th, Z):
    """The largest whitened diameter of the classes of Θ by the nearest atom of Z."""
    Th = np.asarray(Th, float); Z = np.asarray(Z, float)
    lab = np.argmin(_cost_matrix(Th, Z, G), axis=1)
    out = 0.0
    for k in range(Z.shape[0]):
        X = Th[lab == k]
        if len(X) > 1:
            out = max(out, float(np.sqrt(np.max(np.sum((X[:, None, :] - X[None, :, :]) ** 2 / G[None, None, :], axis=2)))))
    return out


def report(A, B, C=None, t0=None, Lg=None, D=None):
    from src.core.sliding import filtered_error
    a, c, s = A['phys']
    win = 5 * PER
    lines = ["# An application: cogging-torque compensation in a position servo (`exp_servo_cogging.py`)", "",
             "The servo ẍ = −τ_cog(x) + u, τ_cog(x) = 0.4 sin(3x + 0.3) + 0.1 sin 7x, tracking x_m = sin t under the sliding controller (s = ė + λe); "
             f"K₀ = {K0} saturated Gaussian kernels on the position, masses 1/{K0}, σ₀ = {SIGMA0:g}; the averaged slow loop along the position orbit.", "",
             "## A. The local minimum ϑ^b(q), the rates, and the global picture (no dynamics)", "",
             f"ϑ^b(q), the nondegenerate local minimum of the unregularised weighted loss that BFGS reaches from START: a = ({', '.join(f'{x:.3f}' for x in a[:, 0])}), "
             f"c = ({', '.join(f'{x:.3f}' for x in c[:, 0])}), s = ({', '.join(f'{x:.3f}' for x in s)}); gradient residual {A['res']:.1e}, relative fit error "
             f"{A['fit']:.4e}. The whitened unregularised Hessian Ã's spectrum (units of σ₀): " + ", ".join(f"{x / SIGMA0:.2e}" for x in A['ev'])
             + f" — {int(np.sum(A['ev'] < 0.1 * SIGMA0))} modes below 0.1σ₀.", "",
             "The collective rates (units of σ₀; the slowest three): the fixed prior's ã_j + σ₀: " + ", ".join(f"{(x + SIGMA0) / SIGMA0:.4f}" for x in A['ev'][:3])
             + "; the adaptive prior's λ_j^− at η = σ₀: " + ", ".join(f"{x / SIGMA0:.2e}" for x in A['rates'][SIGMA0][:3]) + ", at η = ∞: "
             + ", ".join(f"{x / SIGMA0:.2e}" for x in A['rates'][np.inf][:3]) + ".", "",
             "The residual curvature: λ⁺_k/σ₀ = " + ", ".join(f"{x / SIGMA0:.4f}" for x in A['lam']) + " (Proposition D‴.28(iii): margins σ₀ − (q_k/w_k)λ⁺_k); the gauge set "
             "W_g: w > (" + ", ".join(f"{x:.4f}" for x in A['lo']) + f"). The fixed matching prior's floor law (D‴.28(ii)) per unit of 1/m − 1/M: {A['law']:.4e}.", "",
             f"The misplaced prior (the approximant of the cogging model with the phase {PHASE_PRIOR:g} in place of {PHASE:g}, the engineer's model): its own fit "
             f"error on the true target {A['fit_prior']:.4e}; the fixed prior's replicated stationary point there (residual {A['res_misplaced']:.1e}; the same from "
             f"the prior's atoms and from ϑ^b(q), to {A.get('same_point', float('nan')):.1e}) has the fit error {A['fit_misplaced']:.4e} — the leak's bias — "
             + (f"and Proposition D‴.5's splitting margin {A['margin_misplaced']:+.2e} (κ_R/σ₀ = " + ", ".join(f"{x / SIGMA0:.3f}" for x in A['kappa_misplaced'])
                + "): a saddle, which the m-systems leave by splitting their classes (Part B). At the matching prior the replicated point ϑ^b(q) has the "
                f"margin {A['margin_matching']:+.2e} at w = q." if 'margin_misplaced' in A else "")]
    if 'global' in A:
        gl = A['global']
        near = gl['dist'] < 1e-3
        below = gl['fits'] < gl['fit_b'] - 1e-6
        sat_txt = (f", and those all with a saturated amplitude (the largest |α| of each at least {gl['amax'][below].min():.1f}; a = a_max·tanh(α/a_max))"
                   if below.any() and np.all(gl['amax'][below] > 4.0 * np.arctanh(0.999)) else (f" (the largest |α| of each from {gl['amax'][below].min():.1f})" if below.any() else ""))
        lines += ["", f"The global picture: of {len(gl['fits'])} BFGS runs from random starts (physical parameters uniform in their boxes), {int(np.sum(near))} end at "
                  f"ϑ^b(q) (up to relabelling), {int(np.sum(below))} below its fit error {gl['fit_b']:.4e}" + sat_txt + "; the end fit errors' "
                  f"quartiles " + ", ".join(f"{x:.4f}" for x in np.percentile(gl['fits'], [0, 25, 50, 75])) + f". The best end point (fit error {gl['best_fit']:.4e}): "
                  f"a = (" + ", ".join(f"{x:.3f}" for x in gl['best_phys'][0]) + "), c = (" + ", ".join(f"{x:.3f}" for x in gl['best_phys'][1]) + "), s = ("
                  + ", ".join(f"{x:.3f}" for x in gl['best_phys'][2]) + f"). With its {gl['n_sat']} saturated amplitudes held at |α| = "
                  + ", ".join(f"{h['alpha']:g}" for h in gl['hold']) + " (|a| = " + ", ".join(f"{h['a']:.3f}" for h in gl['hold']) + ") and the rest re-optimised, the "
                  "fit error is " + ", ".join(f"{h['fit']:.4e}" for h in gl['hold']) + ": the infimum of the loss is approached at the saturation, not attained — "
                  "globally the servo is Proposition D‴.27's case (critical points at infinity), and ϑ^b(q) is a local minimum, where Proposition D‴.28 applies."]
    lines += ["", f"## B. The averaged loop ({S_B:g} slow units, {N_SEEDS} seeds)", "",
              "The fit error at s = 600, 1500 and the end (mean over the seeds), ‖f̂_m − f̂_M‖² at the end (mean, min–max), and the largest whitened class "
              "diameter at the end over the seeds (the classes by the nearest atom of the design's prior at the end; —: σ-modification and no leak). Start: "
              "warm clouds around ϑ^b(q) (oracle — computed from the true torque) or around the engineer's model.", "",
              "| start | prior | design | m | fit error at 600, 1500, end | ‖f̂_m − f̂_M‖² at the end | class diameter |", "|---|---|---|---|---|---|---|"]
    for st in ('oracle', 'engineer'):
        for pp in (PHASE, PHASE_PRIOR):
            for dsg in DESIGNS:
                rs = [r for r in B if r['design'] == dsg and r['prior_phase'] == pp and r.get('start', 'oracle') == st]
                if not rs:
                    continue
                for m in MS:
                    fits = np.array([[r['rec'][m]['fit'][int(np.argmin(np.abs(r['rec'][m]['s'] - t)))] for t in (600.0, 1500.0, S_B)] for r in rs])
                    h2 = np.array([float(np.mean((r['rec'][m]['fh'][-1].astype(float) - r['rec'][M_REF]['fh'][-1].astype(float)) ** 2)) for r in rs])
                    pname = '—' if dsg in ('sigma-mod', 'no-leak') else ('matching' if pp == PHASE else f'misplaced (phase {pp:g})')
                    diam = '—' if dsg in ('sigma-mod', 'no-leak') else f"{max(class_diameter(*r['final'][m]) for r in rs):.2f}"
                    lines.append(f"| {st} | {pname} | {NAMES[dsg]} | {m} | " + ", ".join(f"{x:.3e}" for x in fits.mean(axis=0))
                                 + f" | {h2.mean():.2e} ({h2.min():.1e}–{h2.max():.1e}) | {diam} |")
    if Lg:
        lines += ["", "The collapsed K₀-atom flows toward the misplaced prior (one particle per class), from ϑ^b(q) and from the engineer's model: the fit error and "
                  "the parameter distance from ϑ^b(q) (largest coordinate difference, up to relabelling) at s = " + ", ".join(f"{x:g}" for x in S_LONG) + ".", "",
                  "| design | start | fit error | distance from ϑ^b(q) | largest |α| at the end |", "|---|---|---|---|---|"]
        for r in Lg:
            lines.append(f"| {NAMES[r['design']]} | {r['start']} | " + ", ".join(f"{x:.4e}" for x in r['fit']) + " | " + ", ".join(f"{x:.3f}" for x in r['dist'])
                         + f" | {r['amax'][-1]:.2f} |")
    if C:
        lines += ["", f"## C. The real servo (k_s = {KS:g}, λ = {LAM:g}, dt = {DT_REAL:g}, {T_C:g} time units), seed 0, from the oracle-start end states", "",
                  "The position-tracking error's rms over the five periods ending at t = 31, 500 and the end (records every "
                  f"{REC_C * DT_REAL:g} time units), against the averaged predictions at the parameters reached then: the residual r(x_m) filtered by "
                  "(d/dt + λ)⁻¹(d/dt + k_s)⁻¹ (exact) and by (d/dt + λ)⁻¹ after the slaved s = r/k_s; their ratios.", "",
                  "| prior | design | m | e_rms at t = 31, 500, end | exact-filter prediction | ratios | slaved prediction | ratios |", "|---|---|---|---|---|---|---|---|"]
        for cc in sorted(C, key=lambda c: (c['prior_phase'], DESIGNS.index(c['design']), c['m'])):
            er = []
            for tt in (win, 500.0, T_C):
                msk = (cc['ts'] > tt - win) & (cc['ts'] <= tt + 1e-9)
                er.append(float(np.sqrt(np.mean(cc['e2'][msk]))))
            pname = '—' if cc['design'] in ('sigma-mod', 'no-leak') else ('matching' if cc['prior_phase'] == PHASE else 'misplaced')
            pr = [cc['e_pred'][tt] for tt in (win, 500.0, T_C)]
            if 'e_pred_exact' in cc:
                px = [cc['e_pred_exact'][tt] for tt in (win, 500.0, T_C)]
            else:
                net_, tb_, prior_, loop_, eta_, plan_ = design_loop(cc['design'], cc['prior_phase'])
                px = [float('nan'), float('nan'), float(np.sqrt(np.mean(filtered_error(loop_.F - loop_.f_hat(np.asarray(cc['final'], float)), PER, KS, LAM, exact=True)[1] ** 2)))]
            lines.append(f"| {pname} | {NAMES[cc['design']]} | {cc['m']} | " + ", ".join(f"{x:.4e}" for x in er) + " | " + ", ".join(f"{x:.4e}" for x in px) + " | "
                         + ", ".join(f"{a_ / b_:.4f}" for a_, b_ in zip(er, px)) + " | " + ", ".join(f"{x:.4e}" for x in pr) + " | "
                         + ", ".join(f"{a_ / b_:.4f}" for a_, b_ in zip(er, pr)) + " |")
    if D:
        lines += ["", "## D. Phase-plane kernels: the whitened Gauss–Newton spectrum along the orbit (Proposition D‴.29(iv))", "",
                  "K = 5 isotropic kernels (tests/test_sliding.py's atoms), weights 1/K, σ₀ = 0.01: the number of eigenvalues below 10⁻⁶σ₀ and below 10⁻¹⁵ of the "
                  "largest, the smallest (units of σ₀), and the (K + 1)-th relative to the largest.", "",
                  "| orbit | dimension | below 10⁻⁶σ₀ | below 10⁻¹⁵ of the largest | smallest/σ₀ | (K + 1)-th/largest |", "|---|---|---|---|---|---|"]
        for r in D:
            lines.append(f"| {r['name']} | {r['dim']} | {r['n_small']} | {r['n_zero']} | {r['smallest']:.2e} | {r['next_rel']:.2e} |")
    lines += ["", "Numbers are produced by `experiments/exp_servo_cogging.py`; nothing is edited by hand."]
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'servo_cogging.md'), 'w', encoding='utf-8') as f:
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
    os.makedirs(OUT, exist_ok=True)
    raw = os.path.join(OUT, 'servo_cogging_raw.npz')
    load = lambda z, k: (z[k].item() if z[k].shape == () else list(z[k])) if k in z.files else None
    if any(f in sys.argv for f in ('--from-raw', '--only-A', '--add-B', '--only-C')):
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            report(z['A'].item(), list(z['B']), load(z, 'C'), t0, Lg=load(z, 'L'), D=load(z, 'D'))
        return
    jobs_B = ([(d, pp, seed) for pp in (PHASE, PHASE_PRIOR) for d in DESIGNS if not (d in ('sigma-mod', 'no-leak') and pp != PHASE) for seed in range(N_SEEDS)]
              + [(d, PHASE_PRIOR if d not in ('sigma-mod', 'no-leak') else PHASE, seed, 'engineer') for d in DESIGNS for seed in range(N_SEEDS)])
    jobs_L = [(d, st) for d in ('nearest-fixed', 'nearest-adaptive', 'nearest-adaptive-slow', 'no-leak') for st in ('oracle', 'engineer')]
    key = lambda r: (r['design'], r['prior_phase'], r['seed'], r.get('start', 'oracle'))
    jkey = lambda j: (j[0], j[1], j[2], j[3] if len(j) > 3 else 'oracle')
    with mp.Pool(N_WORKERS) as pool:
        if '--only-A' in sys.argv:                                          # Parts A and D recomputed, merged with the stored B, L and C
            A = part_A(); A['global'] = part_A_global()
            save['A'] = np.array(A, dtype=object); save['D'] = np.array(part_D(), dtype=object)
            save_raw(raw, **save)
        elif '--add-B' in sys.argv:                                         # the Part B jobs and long flows not yet stored
            B = list(save['B']); have = {key(r) for r in B}
            for r in pool.imap_unordered(run_job, [j for j in jobs_B if jkey(j) not in have], chunksize=1):
                B.append(r); save['B'] = np.array(B, dtype=object); save_raw(raw, **save)
            Lg = list(save['L']) if 'L' in save else []
            haveL = {(r['design'], r['start']) for r in Lg}
            for r in pool.imap_unordered(long_job, [j for j in jobs_L if j not in haveL], chunksize=1):
                Lg.append(r); save['L'] = np.array(Lg, dtype=object); save_raw(raw, **save)
        elif '--only-C' in sys.argv:                                        # Part C rerun from the stored oracle-start Part B (all, or the designs named)
            B = list(save['B'])
            only = [a for a in sys.argv[sys.argv.index('--only-C') + 1:] if a in DESIGNS]
            cj = [(r['design'], r['prior_phase'], m, r['final'][m][0], r['final'][m][1]) for r in B if r['seed'] == 0 and r.get('start', 'oracle') == 'oracle'
                  and (not only or r['design'] in only) for m in (16, 64)]
            keep = [c for c in save['C'] if c['design'] not in only] if only and 'C' in save else []
            save['C'] = np.array(keep + pool.map(real_job, cj, chunksize=1), dtype=object)
            save_raw(raw, **save)
        else:
            A = part_A(); A['global'] = part_A_global()
            Dd = part_D()
            B = pool.map(run_job, jobs_B, chunksize=1)
            save = dict(A=np.array(A, dtype=object), B=np.array(B, dtype=object), D=np.array(Dd, dtype=object))
            save_raw(raw, **save)
            save['L'] = np.array(pool.map(long_job, jobs_L, chunksize=1), dtype=object)
            save_raw(raw, **save)
            cj = [(r['design'], r['prior_phase'], m, r['final'][m][0], r['final'][m][1]) for r in B if r['seed'] == 0 and r.get('start', 'oracle') == 'oracle'
                  for m in (16, 64)]
            save['C'] = np.array(pool.map(real_job, cj, chunksize=1), dtype=object)
            save_raw(raw, **save)
    with np.load(raw, allow_pickle=True) as z:
        report(z['A'].item(), list(z['B']), load(z, 'C'), t0, Lg=load(z, 'L'), D=load(z, 'D'))


if __name__ == '__main__':
    main()
