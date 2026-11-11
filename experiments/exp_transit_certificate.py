"""
Theorem D‴.55: the transit's certificate — exact recovery of separated features by their starting kernels of largest correlation
(the nearest in the continuum), the path's skeleton and the fitted skeleton as predictors of the outcome at a finite separation
====================================================================================================================================

The sliding inverse scale space (`rich_transit_nd.sliding_iss_nd`, in one dimension on a 1-D `BumpsND` — the same equations as
`rich_transit.sliding_iss` with the exact Jacobian, which keeps a block of a large kernel and a seed from stalling the integrator
(referee round 46): the limit of the free law with moving centres from small masses) on random targets of K separated Gaussian
features, against the certificate on the starting set (the samples' inner product; e_k the starting point of largest correlation with
the k-th feature, the nearest in the continuum) and the path's skeletons; the unpolished flow (polish off) as the audit of the
least-squares polish on A and D(b):

A. One dimension (D‴.48's dictionary: 21 kernels of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]): three features at
   least DMIN_A apart, amplitudes ±U(0.4, 1): the certificate's margin m = 1 − max_{j∉E}|η(c_j)| against the path's outcome — exact with
   the nearest kernels (the active set E, the residual zero, no exits), exact with others (three kernels, the residual zero), or not.
B. The amplitudes' range: the same positions and signs with |A_k| = 10^{−r·u_k} (u_k uniform in [0, 1]) for r = R_B — the certificate
   does not change; does the outcome?
C. The plane (D‴.54's lattice: 17 × 17 kernels of width 0.15, spacing 0.125, 41² samples of [−1, 1]²): three features at least
   DMIN_C apart, amplitudes ±U(0.4, 1).
D. Designed families (referee round 43): pairs on starting points with amplitude ratios, pairs at s·σ over a grid cell, every order of
   the amplitudes.
E. The two-feature map: the first feature on a starting point, the second at the distance D with the amplitude ±10^ℓ, over a grid of
   (D, ℓ): the fitted skeleton's outcome everywhere (no flow), the threshold separation D*(ℓ), the other margins against it, the linear
   skeleton's error against the fitted one's in the ratio (its order), and the path at sampled points against it.
F. A wide range: three features at 5σ and 3σ, the amplitudes 1, 10⁻³, 10⁻⁶ in all six orders, same and alternating signs (exits only
   below 10⁻¹² of the largest coefficient).
G. The polish's audit on the rest (referee round 140): the unpolished flow (polish off) on B, C, D(a) and D(c), so that with A and D(b)
   every one of the 1154 targets of A–D has it.
H. Baselines that do not follow the path (referee round 140), on every target of A–D: Duval–Peyré's vanishing-derivative pre-certificate
   η_V (`eta_v`: the least-norm p with ⟨pφ(c*_k)⟩ = sign A_k, ⟨p∇φ(c*_k)⟩ = 0 on the samples, η_V = ⟨φ(·)p⟩) — its margin on the starting
   points off E (an amplitude-free predictor of exactness with E, like m) and its nondegeneracy on a fine grid (|η_V| < 1 beyond EX_R·σ
   of the features, s_k∇²η_V(c*_k) negative definite), which certifies that basis pursuit over measures (BLASSO as its regularisation
   vanishes) recovers the features —; basis pursuit itself on a fine grid with the features' positions added (`bp_certificate`: the
   least sup-norm certificate on the grid, by a linear program; < 1 iff the features are its unique solution there); and continuous
   OMP (Elvira–Gribonval–Soussen–Herzet: K steps, each the atom of largest normalised correlation on a fine grid refined locally, the
   amplitudes by least squares at fixed positions), without and with a joint least-squares slide of all selected atoms after each step.

K. Part I's failures by type (referee round 149), from the stored solutions without a new solve (`k_classify`, `part_K`): a feature
   lost, a spurious nonzero (beyond K_RD·σ of every feature), a displaced one (beyond the tolerance, within K_RD·σ), a cluster of the
   wrong sign; and Part G's audit as entry sequences (the unpolished flow's entries against the polished path's).
L. Basis pursuit's least certificate in the plane exactly (after round 149; Part J bounded it by the jets): Part H's linear program in
   full, in the basis of the tensor grid's singular vectors (`_kron_basis`), by a primal–dual interior-point method whose normal matrix is
   formed through the tensor factors, with certified upper and lower bounds (`bp_certificate_kron`), at Part H's truncation of the span
   and at two coarser ones (L_CUTS).
M. The BLASSO in the plane (after round 149): Part I's problem, solver, λ's and recovery rule on the plane targets, the plane's tolerance
   I_TOL·H_BP[1], the Gram matrix through the samples' tensor factors (`KronGram`, `blasso_kron`); the failures by Part K's types.

Usage: python experiments/exp_transit_certificate.py [--only A,B,C,D,E,F,G,H,I,J,K,L,M] [--from-raw] → results/chaos/transit_certificate.md (+ _raw.npz).
"""
import copy
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
from src.core import rich_transit_nd as RN
import exp_cohesion_leak as E
import exp_mirror_transit as X
import exp_transit_mass as TM
import exp_transit_nd as TN

OUT = os.environ.get('TC_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TC_WORKERS', '2'))
K_F = 3
SEED_A, N_A, DMIN_A = 551, 80, (0.25, 0.3, 0.35, 0.45, 0.6)
SEED_B, N_B, DMIN_B, R_B = 552, 40, 0.35, (0, 1, 2, 3)
SEED_C, N_C, DMIN_C = 553, 40, (0.3, 0.45, 0.6)
RMS_TOL = 1e-8
BANDS = (0.5, 0.8, 1.0)                                                # the band |m| ≤ Cq around the certificate's threshold
T_MAX = 600                                                            # D: a path's wall-clock limit (s)
LAT_OFF, LAT_D, RATIOS_D = 1e-6, (0.3, 0.4, 0.5), (1e-3, 0.02, 0.1, 1 / 3, 1.0, 3.0, 10.0, 50.0, 1e3)
SEED_D, SEPS_D, N_OFF_D = 554, (1.5, 2.0, 2.5, 3.0, 3.5, 4.0), 20       # D(b): pairs at s·σ, the midpoint uniform over a grid cell
N_ORD_D, RANGE_D = 30, 10.0                                            # D(c): A's targets outside 0.8q, the amplitudes 1, 10^{−1/2}, 10^{−1} in every order
D_E, L_E = np.round(np.arange(0.12, 0.6001, 0.01), 10), np.round(np.arange(-3.0, 3.0001, 0.125), 10)  # E: the map (D, log₁₀ ratio)
SEED_E, N_EV, LMAX_EV, T_EV = 555, 120, 2.5, 300                      # E: paths per sign checking the map (|log₁₀ ratio| ≤ LMAX_EV)
M0_E = 0.05                                                            # E: the second feature's own near ties (m₀ ≤ M0_E) set apart
F_SEPS, F_MAGS, F_OFF = (0.6, 0.36), (1.0, 1e-3, 1e-6), 0.013          # F: 5σ and 3σ, the range 10⁶, the features off the grid
EX_R = 0.1                                                             # H: η_V's nondegeneracy checked beyond EX_R·σ of the features
H_ETA, H_BP, H_OMP = (0.002, 0.01), (0.005, 0.025), (0.001, 0.01)      # H: the grids' spacings (one dimension, the plane)
H_BOX = ((-1.2, 1.2), (-1.0, 1.0))                                     # H: the grids' extent per coordinate (one dimension, the plane)
I_LAMS, I_TOL, I_KKT = (1e-1, 1e-2, 1e-3, 1e-4), 2, 1e-12              # I: λ/λ_max; recovery within I_TOL grid spacings; the KKT tolerance
I_CHUNK = 100                                                          # I, J: targets per resumable chunk
K_RD = 0.5                                                             # K: beyond the tolerance, within K_RD·σ of its feature a nonzero is displaced, beyond it spurious
L_CUTS = (1e-13, 1e-10, 1e-7)                                          # L: the span's truncation s_is_j > c·s_0² (10⁻¹³ is Part H's)
L_TOL, L_IT, L_STALL, L_REG, L_CHUNK = 1e-9, 80, 8, 1e-15, 20     # L: the bracket's target width, the iteration cap, the stall, the shift, the chunk


def entrants(K, c0, cs):
    """E: for each feature the starting point of largest correlation with it on the samples (the nearest in the continuum; on samples
    they differ near the window's edge — referee round 46), and the correlations ⟨φ(c_j)φ(c*_k)⟩."""
    c0 = np.asarray(c0, float).reshape(len(c0), -1)
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    P0 = K.phi(c0 if c0.shape[1] > 1 else c0[:, 0])
    Pf = K.phi(cs if cs.shape[1] > 1 else cs[:, 0])
    gram = P0.T @ Pf / P0.shape[0]
    return [int(np.argmax(gram[:, k])) for k in range(len(cs))], gram


def certificate(K, c0, cs, As):
    """η at every starting point (the samples' inner product), E (`entrants`), the margin m off E, and the cross-talk q (continuum,
    q = max_k Σ_{l≠k} G((Δ_kl − 2R)₊)/G(R), G(d) = e^{−d²/(4σ²)})."""
    c0 = np.asarray(c0, float).reshape(len(c0), -1)
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    d0 = np.linalg.norm(c0[:, None, :] - cs[None, :, :], axis=2)
    e, gram = entrants(K, c0, cs)
    own = np.array([gram[e[k], k] for k in range(len(cs))])
    eta = (np.sign(As)[None, :] * gram / own[None, :]).sum(1)
    off = np.array([j for j in range(len(c0)) if j not in e])
    sig = K.sig
    G = lambda d: np.exp(-np.asarray(d, float) ** 2 / (4 * sig ** 2))
    R = float(max(d0[e[k], k] for k in range(len(cs))))
    D = np.linalg.norm(cs[:, None, :] - cs[None, :, :], axis=2)
    q = float(max(sum(G(max(D[k, l] - 2 * R, 0.0)) for l in range(len(cs)) if l != k) for k in range(len(cs)))) / float(G(R))
    # the relative margin (referee round 43): the entrants' own values s_kη(c_{e_k}) — 1 up to the cross-talk — against the largest off E;
    # the single-feature (Voronoi) margin: each feature alone, no signs or sums
    m_rel = float(min(np.sign(As[k]) * eta[e[k]] for k in range(len(cs))) - np.abs(eta[off]).max())
    m_vor = float(min(1 - max(gram[j, k] for j in range(len(c0)) if j != e[k]) / gram[e[k], k] for k in range(len(cs))))
    return dict(eta=eta, e=e, m=float(1 - np.abs(eta[off]).max()), argm=int(off[np.argmax(np.abs(eta[off]))]), q=q, R=R,
                unique=len(set(e)) == len(e), m_rel=m_rel, m_vor=m_vor)


def skeleton(K, c0, cs, As, fit=False, atoms=False):
    """The sliding path's skeleton (referee round 43's screening): each block replaced by its exact fit — the fitted features' kernels at
    their centres, the residual between entries the target's projection off their tangent space span{φ(c*_m), ∇_cφ(c*_m)} — the duals
    piecewise linear in τ at ⟨rφ(c_j)⟩ from the starting points. It runs the entries in order (each the first dual to reach ±1) until all K
    features are fitted or a kernel off E enters. Returns the entrants, whether they are E (ok), and the margin: 1 − max |u_j| over the
    entries and the kernels off E, or, if a kernel off E enters, minus the rightful entrants' least distance from their walls then.
    fit: the fitted skeleton — each block at its own sign-constrained least-squares fit started at its features (`polish_block`: the
    block's stationary point near them, the projection's nonlinear completion) — the skeleton without its linearisation. atoms: the
    baseline without the tangent directions — the fitted features' atoms A_mφ(c*_m) subtracted (referee round 46). 'phase': the number
    of features fitted when the outcome is decided (the intruder's entry, or the phase of the least margin)."""
    c0 = np.asarray(c0, float).reshape(len(c0), -1)
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    sq = lambda c: c if c.shape[1] > 1 else c[:, 0]
    nx = K.xs.size if hasattr(K, 'xs') and np.ndim(K.xs) == 1 else K.N
    Pf = K.phi(sq(cs))
    f = Pf @ np.asarray(As, float)
    P0 = K.phi(sq(c0))
    e = entrants(K, c0, cs)[0]
    off = [j for j in range(len(c0)) if j not in e]
    u, act, fitted, ent, worst, wphase = np.zeros(len(c0)), np.zeros(len(c0), bool), [], [], 0.0, 0
    while len(fitted) < len(cs):
        if fitted and atoms:
            r = f - Pf[:, fitted] @ np.asarray(As, float)[fitted]
        elif fitted and fit:
            a_, s_ = np.asarray(As, float)[fitted], np.sign(np.asarray(As, float)[fitted])
            if cs.shape[1] > 1:
                w_, c_ = RN.polish_block_nd(K, f, a_, cs[fitted], s_)
                r = f - K.synth(c_, w_)
            else:
                w_, c_ = RT.polish_block(K, f, a_, cs[fitted, 0], s_)
                r = f - K.phi(c_) @ w_
        elif fitted:
            Dp = K.dphi(sq(cs[fitted])).reshape(nx, -1)
            B = np.c_[Pf[:, fitted], Dp]
            r = f - B @ np.linalg.lstsq(B, f, rcond=None)[0]
        else:
            r = f
        v = P0.T @ r / nx
        ina = ~act & (v != 0)
        if not ina.any():
            break
        dt = np.full(len(c0), np.inf)
        dt[ina] = (1 - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])
        j = int(np.argmin(dt))
        u[ina] += dt[j] * v[ina]
        u[j] = np.sign(u[j])
        act[j] = True
        ent.append(j)
        if j not in e or e.index(j) in fitted:
            # an intruder: the margin is minus the rightful entrants' distance from their walls at its entry (a near tie when small)
            gap = min(1 - abs(u[e[m]]) for m in range(len(cs)) if m not in fitted)
            return dict(entrants=ent, ok=False, margin=float(-gap), phase=len(fitted))
        cur = float(np.abs(u[[i for i in off if not act[i]]]).max()) if any(not act[i] for i in off) else 0.0
        if cur > worst:
            worst, wphase = cur, len(fitted)
        fitted.append(e.index(j))
    return dict(entrants=ent, ok=True, margin=float(1 - worst), phase=wphase)


def outcome(R, e, K):
    act = sorted(int(j) for j in np.flatnonzero(R['act']))
    clean = R['hist'][-1][3] < RMS_TOL and all(ev[2] == 'enter' for ev in R['events'])
    if clean and act == sorted(e) and len(R['events']) == K:
        return 'nearest'
    if R['hist'][-1][3] < RMS_TOL and len(act) == K:
        return 'other'
    return 'not'


def targets(rng, n, dmin, lo=-0.7, hi=0.7, dim=1):
    out = []
    while len(out) < n:
        cs = rng.uniform(lo, hi, (K_F, dim))
        if min(np.linalg.norm(cs[i] - cs[j]) for i in range(K_F) for j in range(i)) < dmin:
            continue
        out.append((cs[:, 0] if dim == 1 else cs, rng.choice([-1.0, 1.0], K_F), rng.uniform(0.4, 1.0, K_F), rng.uniform(0, 1, K_F)))
    return out


def path_1d(L, K, f, c0, polish=True, exit_tol=1e-6):
    """The sliding path in one dimension by `sliding_iss_nd` on the 1-D samples (the exact Jacobian)."""
    return RN.sliding_iss_nd(L, RN.BumpsND(K.xs[:, None], K.sig), f, np.asarray(c0, float)[:, None], np.ones(len(c0)), polish=polish,
                             exit_tol=exit_tol)


def job_1d(args):
    """One target: the certificate and the path (polished); with flow=True also the unpolished flow's outcome (the polish's audit)."""
    cs, As = args[:2]
    flow = len(args) > 2 and args[2]
    exit_tol = args[3] if len(args) > 3 else 1e-6
    warnings.filterwarnings('ignore')
    L, K, c0 = TM.setup()
    cert = certificate(K, c0, cs, As)
    if not cert['unique']:
        return dict(cs=cs, As=As, skip=True)
    f = X.bumps_target(K, cs, As)
    R = path_1d(L, K, f, c0, exit_tol=exit_tol)
    out = dict(cs=cs, As=As, skip=False, m=cert['m'], q=cert['q'], argm=cert['argm'], e=cert['e'], out=outcome(R, cert['e'], len(cs)),
               order=[int(j) for t, j, k, c in R['events']], rms=R['hist'][-1][3], m_rel=cert['m_rel'], m_vor=cert['m_vor'],
               rej=int(R['polish_rejected']))
    if flow:
        t0 = time.time()
        Rf = path_1d(L, K, f, c0, polish=False, exit_tol=exit_tol)
        out.update(out_flow=outcome(Rf, cert['e'], len(cs)), order_flow=[int(j) for t, j, k, c in Rf['events']],
                   flow_relaxed=bool(all(Rf['relaxed'])), flow_secs=time.time() - t0)
    return out


def job_d(args):
    """D: `job_1d` under a wall-clock limit T_MAX (out = 'timeout', the certificate still recorded); the family's tag kept."""
    import signal
    tag, cs, As = args[:3]
    tmax = int(args[3]) if len(args) > 3 else T_MAX
    flow = len(args) > 4 and args[4]
    exit_tol = args[5] if len(args) > 5 else 1e-6

    def alarm(signum, frame):
        raise TimeoutError
    signal.signal(signal.SIGALRM, alarm)
    signal.alarm(tmax)
    t0 = time.time()
    try:
        r = job_1d((cs, As, flow, exit_tol))
    except TimeoutError:
        L, K, c0 = TM.setup()
        cert = certificate(K, c0, cs, As)
        r = dict(cs=cs, As=As, skip=not cert['unique'], m=cert['m'], q=cert['q'], argm=cert['argm'], e=cert['e'], out='timeout', order=[],
                 rms=np.nan, m_rel=cert['m_rel'], m_vor=cert['m_vor'], rej=-1)
    finally:
        signal.alarm(0)
    r.update(tag=tag, secs=time.time() - t0)
    return r


def part_D(A):
    """Designed families (referee round 43): (a) pairs on starting points (LAT_OFF off) at the distances LAT_D, same and opposite signs,
    the second amplitude the first's times RATIOS_D; (b) pairs at s·σ (SEPS_D), equal amplitudes, same and opposite signs, the midpoint
    uniform over a grid cell; (c) N_ORD_D of A's targets with |m| > 0.8q, the amplitudes 1, RANGE_D^{−1/2}, RANGE_D^{−1} in every order."""
    from itertools import permutations
    from multiprocessing import Pool
    jobs = []
    for D in LAT_D:
        for sg in (1.0, -1.0):
            for rt in RATIOS_D:
                jobs.append((('a', D, sg, rt), np.array([-0.1, -0.1 + D]) + LAT_OFF, np.array([1.0, sg * rt])))
    rng = np.random.default_rng(SEED_D)
    for s_ in SEPS_D:
        for sg in (1.0, -1.0):
            for i in range(N_OFF_D):
                mid = rng.uniform(-0.05, 0.05)
                jobs.append((('b', s_, sg, i), np.array([mid - s_ * 0.06, mid + s_ * 0.06]), np.array([1.0, sg]), T_MAX, True))
    out_ = [r for r in A if not r['skip'] and abs(r['m']) > 0.8 * r['q']]
    pick = sorted(rng.choice(len(out_), size=min(N_ORD_D, len(out_)), replace=False))
    mags = RANGE_D ** -np.array([0.0, 0.5, 1.0])
    for t in pick:
        r = out_[t]
        for perm in permutations(range(3)):
            jobs.append((('c', int(t), perm), np.array(r['cs']), np.sign(r['As']) * mags[list(perm)]))
    with Pool(N_WORKERS) as pool:
        return pool.map(job_d, jobs, chunksize=1)


def e_point(args):
    """E: one point of the two-feature map — the margins (the certificate's, the relative, the single-feature) and both skeletons'."""
    D, l, sg = args
    warnings.filterwarnings('ignore')
    L, K, c0 = TM.setup()
    cs = np.array([-0.1, -0.1 + D]) + LAT_OFF
    As = np.array([1.0, sg * 10.0 ** l])
    cert = certificate(K, c0, cs, As)
    if not cert['unique']:
        return dict(D=D, l=l, sg=sg, skip=True)
    sk, fs = skeleton(K, c0, cs, As), skeleton(K, c0, cs, As, fit=True)
    return dict(D=D, l=l, sg=sg, skip=False, m=cert['m'], q=cert['q'], m_rel=cert['m_rel'], m_vor=cert['m_vor'], m_skel=sk['margin'],
                m_fit=fs['margin'], fit_ok=bool(fs['ok']), m_atoms=skeleton(K, c0, cs, As, atoms=True)['margin'], phase=fs['phase'])


def part_E():
    """The two-feature map on the grid D_E × L_E for both signs; the path (under T_EV s) at N_EV sampled points per sign, |ℓ| ≤ LMAX_EV."""
    from multiprocessing import Pool
    pts = [(float(D), float(l), sg) for sg in (1.0, -1.0) for D in D_E for l in L_E]
    rng = np.random.default_rng(SEED_E)
    pick = []
    for sg in (1.0, -1.0):
        ch = [p for p in pts if p[2] == sg and abs(p[1]) <= LMAX_EV]
        pick += [ch[i] for i in sorted(rng.choice(len(ch), size=N_EV, replace=False))]
    with Pool(N_WORKERS) as pool:
        grid = pool.map(e_point, pts, chunksize=40)
        paths = pool.map(job_d, [(('e', D, l, sg), np.array([-0.1, -0.1 + D]) + LAT_OFF, np.array([1.0, sg * 10.0 ** l]), T_EV)
                                 for D, l, sg in pick], chunksize=1)
    return dict(grid=grid, paths=paths)


def job_2d(args):
    cs, As = args
    warnings.filterwarnings('ignore')
    L = X.law()
    K, c0 = TN.full()
    cert = certificate(K, c0, cs, As)
    if not cert['unique']:
        return dict(cs=cs, As=As, skip=True)
    R = RN.sliding_iss_nd(L, K, K.synth(cs, As), c0, np.ones(len(c0)))
    return dict(cs=cs, As=As, skip=False, m=cert['m'], q=cert['q'], argm=cert['argm'], e=cert['e'], out=outcome(R, cert['e'], len(cs)),
                order=[int(j) for t, j, k, c in R['events']], rms=R['hist'][-1][3], rej=int(R['polish_rejected']))


def part_A():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_A)
    jobs, tags = [], []
    for dmin in DMIN_A:
        for cs, s, a, u in targets(rng, N_A, dmin):
            jobs.append((cs, s * a, True)); tags.append(dmin)
    with Pool(N_WORKERS) as pool:
        rs = pool.map(job_1d, jobs, chunksize=4)
    return [dict(dmin=t, **r) for t, r in zip(tags, rs)]


def part_F():
    """Three features at 5σ and 3σ (off the grid by F_OFF), the amplitudes F_MAGS in all six orders, signs + + + and + − +; exits only
    below 10⁻¹² of the largest coefficient."""
    from itertools import permutations
    from multiprocessing import Pool
    jobs = []
    for s_ in F_SEPS:
        for sg in ((1.0, 1.0, 1.0), (1.0, -1.0, 1.0)):
            for perm in permutations(range(3)):
                jobs.append((('f', s_, sg, perm), np.array([-s_, 0.0, s_]) + F_OFF, np.array(sg) * np.array(F_MAGS)[list(perm)], T_MAX,
                             False, 1e-12))
    with Pool(N_WORKERS) as pool:
        return pool.map(job_d, jobs, chunksize=1)


def part_B():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_B)
    jobs, tags = [], []
    for i, (cs, s, a, u) in enumerate(targets(rng, N_B, DMIN_B)):
        for r in R_B:
            jobs.append((cs, s * 10.0 ** (-r * u))); tags.append((i, r))
    with Pool(N_WORKERS) as pool:
        rs = pool.map(job_1d, jobs, chunksize=4)
    return [dict(cfg=t[0], r=t[1], **x) for t, x in zip(tags, rs)]


def part_C():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    jobs, tags = [], []
    for dmin in DMIN_C:
        for cs, s, a, u in targets(rng, N_C, dmin, dim=2):
            jobs.append((cs, s * a)); tags.append(dmin)
    with Pool(N_WORKERS) as pool:
        rs = pool.map(job_2d, jobs, chunksize=2)
    return [dict(dmin=t, **r) for t, r in zip(tags, rs)]


def confusion(rows):
    rows = [r for r in rows if not r['skip']]
    pos = [r for r in rows if r['m'] > 0]
    neg = [r for r in rows if r['m'] <= 0]
    cnt = lambda rs, o: sum(r['out'] == o for r in rs)
    return rows, pos, neg, cnt


def fill_margins(rows, dim=1):
    """The skeletons' margins recomputed from the targets at report time (a stored margin would hide a change of `skeleton`: referee
    round 46) — the linear skeleton, the fitted one (and the phase at which it decides) and the atoms-only baseline —, and the relative and
    single-feature margins where a row lacks them."""
    if dim == 1:
        L, K, c0 = TM.setup()
    else:
        K, c0 = TN.full()
    for r in rows:
        if r['skip']:
            continue
        cs, As = np.asarray(r['cs']), np.asarray(r['As'])
        if 'm_rel' not in r:
            c = certificate(K, c0, cs, As)
            r['m_rel'], r['m_vor'] = c['m_rel'], c['m_vor']
        r['m_skel'] = skeleton(K, c0, cs, As)['margin']
        fs = skeleton(K, c0, cs, As, fit=True)
        r['m_fit'], r['phase'] = fs['margin'], fs['phase']
        r['m_atoms'] = skeleton(K, c0, cs, As, atoms=True)['margin']
    return rows


PREDICTORS = (('m', 'the certificate, 1 − max off E'), ('m_rel', 'the relative margin'), ('m_vor', 'the single feature (no signs or sums)'),
              ('m_atoms', 'the atoms subtracted (no tangent directions)'), ('m_skel', "the path's skeleton"), ('m_fit', 'the fitted skeleton'))
SHORT = dict(m='the certificate', m_rel='relative', m_vor='single-feature', m_atoms='atoms subtracted', m_skel='the skeleton',
             m_fit='the fitted skeleton')                              # D(b)'s column labels


def predictors(rows):
    """How often each predictor's sign (and the baselines 'always E' and 'E iff the path's first entrant is in E') matches exactness with
    E; where computed, how often the unpolished flow's outcome is the polished path's."""
    rs = [r for r in rows if not r['skip'] and r['out'] != 'timeout']
    ok = [r['out'] == 'nearest' for r in rs]
    parts = [f"{lab} {sum((r[k] > 0) == o for r, o in zip(rs, ok))}" for k, lab in PREDICTORS if all(k in r for r in rs)]
    first = [bool(r['order']) and r['order'][0] in r['e'] for r in rs]
    txt = (f"of {len(rs)}: " + ', '.join(parts) + f", always E {sum(ok)}, E iff the first entrant is in E {sum(f == o for f, o in zip(first, ok))} "
           f"(the first entrant off E — the outcome decided at the first entry — on {sum(not f for f in first)}")
    if all('phase' in r for r in rs):
        txt += f"; the fitted skeleton decided at its first phase on {sum(r['phase'] == 0 for r in rs)}"
    txt += ")"
    fl = [r for r in rs if 'out_flow' in r]
    if fl:
        txt += (f"; the unpolished flow's outcome the polished path's on {sum(r['out_flow'] == r['out'] for r in fl)} of {len(fl)} (its "
                f"blocks relaxed on {sum(r['flow_relaxed'] for r in fl)})")
    return txt


def band(rs):
    """Outside the band |m| ≤ Cq: how many targets, and on how many the certificate's sign decides; the paths' ends."""
    parts = []
    for C in BANDS:
        ob = [r for r in rs if abs(r['m']) > C * r['q']]
        parts.append(f"C = {C:g}: {sum((r['m'] > 0) == (r['out'] == 'nearest') for r in ob)} of {len(ob)}")
    return ("Outside the band |m| ≤ Cq the sign decides on " + '; '.join(parts) + f". Every path ends at an exact fit with "
            f"{K_F} kernels on {sum(r['out'] != 'not' for r in rs)} of {len(rs)} targets.")


def job_2d_flow(args):
    """G: a target of C, the polished path and the unpolished flow (polish off)."""
    cs, As = args
    warnings.filterwarnings('ignore')
    L = X.law()
    K, c0 = TN.full()
    cert = certificate(K, c0, cs, As)
    if not cert['unique']:
        return dict(cs=cs, As=As, skip=True)
    f = K.synth(cs, As)
    R = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)))
    t0 = time.time()
    Rf = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)), polish=False)
    return dict(cs=cs, As=As, skip=False, e=cert['e'], out=outcome(R, cert['e'], len(cs)), out_flow=outcome(Rf, cert['e'], len(cs)),
                flow_relaxed=bool(all(Rf['relaxed'])), flow_secs=time.time() - t0)


def part_G(B, C, D):
    """The unpolished flow on B, C, D(a) and D(c) (A and D(b) have it already)."""
    from multiprocessing import Pool
    j1 = [(('B', i), np.asarray(r['cs']), np.asarray(r['As']), T_MAX, True) for i, r in enumerate(B) if not r['skip']]
    j1 += [(('D', i), np.asarray(r['cs']), np.asarray(r['As']), T_MAX, True) for i, r in enumerate(D) if r['tag'][0] in 'ac' and not r['skip']]
    with Pool(N_WORKERS) as pool:
        r1 = pool.map(job_d, j1, chunksize=2)
        r2 = pool.map(job_2d_flow, [(np.asarray(r['cs']), np.asarray(r['As'])) for r in C if not r['skip']], chunksize=2)
    ic = [i for i, r in enumerate(C) if not r['skip']]
    return dict(rows=r1, plane=[dict(src=('C', i), **x) for i, x in zip(ic, r2)])


def _k1(K):
    """The 1-D samples as a `BumpsND` (one code path for both dimensions)."""
    return RN.BumpsND(np.asarray(K.xs, float)[:, None], K.sig)


def _grid(dim, h):
    g = np.arange(H_BOX[dim - 1][0], H_BOX[dim - 1][1] + 1e-9, h)
    return g[:, None] if dim == 1 else np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)


def _corr(K, c, q, chunk=4000):
    """Σ_s φ_s(c)q_s at the points c (chunked)."""
    return np.concatenate([K.phi(c[i:i + chunk]).T @ q for i in range(0, len(c), chunk)])


def eta_v(K, cs, As, c0, e):
    """Duval–Peyré's vanishing-derivative pre-certificate on the samples: q = Γ(ΓᵀΓ)⁻¹(sign A, 0), Γ = [φ(c*_k), ∇φ(c*_k)], η_V(c) =
    Σ_sφ_s(c)q_s (η_V(c*_k) = sign A_k, ∇η_V(c*_k) = 0). Returns its margin on the starting points off E, its largest |η_V| on the fine grid
    beyond EX_R·σ of every feature, and the largest eigenvalue of s_k∇²η_V(c*_k) over k (< 0: nondegenerate at the features)."""
    cs = np.asarray(cs, float).reshape(len(As), -1)
    n = cs.shape[1]
    Gm = np.c_[K.phi(cs), K.dphi(cs).reshape(K.N, -1)]
    q = Gm @ np.linalg.solve(Gm.T @ Gm, np.r_[np.sign(As), np.zeros(len(As) * n)])
    off = [j for j in range(len(c0)) if j not in e]
    m_v = float(1 - np.abs(_corr(K, np.asarray(c0, float).reshape(len(c0), -1)[off], q)).max())
    g = _grid(n, H_ETA[n - 1])
    far = np.min(np.linalg.norm(g[:, None, :] - cs[None, :, :], axis=2), axis=1) > EX_R * K.sig
    sup = float(np.abs(_corr(K, g[far], q)).max())
    far1 = np.min(np.linalg.norm(g[:, None, :] - cs[None, :, :], axis=2), axis=1) > K.sig
    sup1 = float(np.abs(_corr(K, g[far1], q)).max())
    Q = K.corr(cs, q, second=True)[2] * K.N
    hmax = float(max(np.linalg.eigvalsh(np.sign(As[k]) * Q[k]).max() for k in range(len(As))))
    return dict(m_v=m_v, sup=sup, sup1=sup1, hess=hmax, nondeg=bool(sup < 1 and hmax < 0))


def bp_certificate(K, cs, As):
    """Basis pursuit on the fine grid, its points within half a spacing of a feature replaced by the feature: the least t with a p such that ⟨pφ(c*_k)⟩ = sign A_k and
    |⟨pφ(c)⟩| ≤ t at every other grid point (a linear program, HiGHS, with p in an orthonormal basis of the span of the columns,
    singular values below 10⁻¹³ of the largest dropped). t < 1 iff the features (with their amplitudes) are the unique
    minimum-ℓ¹ exact fit on that grid (the columns at the features being independent)."""
    from scipy.optimize import linprog
    cs = np.asarray(cs, float).reshape(len(As), -1)
    n = cs.shape[1]
    g = _grid(n, H_BP[n - 1])
    g = g[np.min(np.linalg.norm(g[:, None, :] - cs[None, :, :], axis=2), axis=1) > 0.5 * H_BP[n - 1]]   # the features replace their grid neighbours
    Pg = np.concatenate([K.phi(g[i:i + 4000]) for i in range(0, len(g), 4000)], axis=1)       # (N, grid)
    Ps = K.phi(cs)
    U, sv, _ = np.linalg.svd(np.c_[Pg, Ps], full_matrices=False)       # p in an orthonormal basis of the columns' span: the
    Ub = U[:, sv > sv[0] * 1e-13]                                       # Gaussian columns are nearly dependent, and in the samples'
    Qg, Qs = Pg.T @ Ub, Ps.T @ Ub                                       # own basis HiGHS reports numerical difficulties
    M = Ub.shape[1]
    A_ub = np.r_[np.c_[Qg, -np.ones(len(g))], np.c_[-Qg, -np.ones(len(g))]]
    sol = linprog(np.r_[np.zeros(M), 1.0], A_ub=A_ub, b_ub=np.zeros(2 * len(g)), A_eq=np.c_[Qs, np.zeros(len(As))],
                  b_eq=np.sign(As), bounds=[(None, None)] * M + [(0, None)], method='highs')
    return float(sol.x[-1]) if sol.status == 0 else np.nan


def omp(K, f, n_steps, dim, slide=False):
    """Continuous OMP: n_steps selections, each the grid point of largest |⟨rφ(c)⟩|/⟨φ(c)²⟩^{1/2} refined by Nelder–Mead, the amplitudes
    by least squares at the selected positions; with slide, all selected atoms then jointly refitted (`polish_block_nd`, the signs of
    their amplitudes kept). Returns (rms residual, centres, coefficients)."""
    from scipy.optimize import minimize
    g = _grid(dim, H_OMP[dim - 1])
    nrm = np.sqrt(np.concatenate([(K.phi(g[i:i + 4000]) ** 2).sum(0) for i in range(0, len(g), 4000)]))
    sel = np.zeros((0, dim))
    w = np.zeros(0)
    r = f.copy()
    for _ in range(n_steps):
        cc = np.abs(_corr(K, g, r)) / nrm
        c1 = g[int(np.argmax(cc))]
        obj = lambda c: -abs(float(K.phi(c[None, :])[:, 0] @ r)) / float(np.sqrt((K.phi(c[None, :])[:, 0] ** 2).sum()))
        c1 = minimize(obj, c1, method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-16, maxiter=4000)).x
        sel = np.r_[sel, c1[None, :]]
        Ph = K.phi(sel)
        w = np.linalg.lstsq(Ph, f, rcond=None)[0]
        if slide:
            w, cflat = RN.polish_block_nd(K, f, w, sel, np.sign(w))
            sel = np.asarray(cflat, float).reshape(len(w), dim)
        r = f - K.phi(sel) @ w
    return float(np.sqrt(np.mean(r ** 2))), sel, w


def h_job(args):
    """H: one target's baselines (η_V, basis pursuit on a grid, continuous OMP without and with the slide)."""
    src, cs, As, dim = args
    warnings.filterwarnings('ignore')
    if dim == 1:
        L, K0, c0 = TM.setup()
        K = _k1(K0)
    else:
        K, c0 = TN.full()
    cs = np.asarray(cs, float).reshape(len(As), -1)
    e = entrants(K0 if dim == 1 else K, c0, cs if dim > 1 else cs[:, 0])[0]
    f = K.phi(cs) @ np.asarray(As, float)
    t0 = time.time()
    ev = eta_v(K, cs, np.asarray(As, float), c0, e)
    t_bp = bp_certificate(K, cs, np.asarray(As, float)) if dim == 1 else np.nan      # the plane's linear program is too large here
    o_rms, o_c, o_w = omp(K, f, len(As), dim)
    s_rms, s_c, s_w = omp(K, f, len(As), dim, slide=True)
    dmax = float(np.max(np.min(np.linalg.norm(o_c[:, None, :] - cs[None, :, :], axis=2), axis=1))) / K.sig
    return dict(src=src, **ev, t_bp=t_bp, omp_rms=o_rms, omp_dmax=dmax, somp_rms=s_rms, secs=time.time() - t0)


def part_H(A, B, C, D):
    from multiprocessing import Pool
    jobs = [(('A', i), r['cs'], r['As'], 1) for i, r in enumerate(A) if not r['skip']]
    jobs += [(('B', i), r['cs'], r['As'], 1) for i, r in enumerate(B) if not r['skip']]
    jobs += [(('D', i), r['cs'], r['As'], 1) for i, r in enumerate(D) if not r['skip']]
    jobs += [(('C', i), r['cs'], r['As'], 2) for i, r in enumerate(C) if not r['skip']]
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, jobs, chunksize=4)


# ---- I and J (after round 140): BLASSO at a finite regularisation; basis pursuit in the plane -------------------------------------------

def bp_grid(K, cs, h):
    """Part H's basis-pursuit grid of spacing h with the features added: its points within half a spacing of a feature dropped, the
    features appended last (positions (M, n))."""
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    g = _grid(cs.shape[1], h)
    g = g[np.min(np.linalg.norm(g[:, None, :] - cs[None, :, :], axis=2), axis=1) > 0.5 * h]
    return np.r_[g, cs]


def lasso_fs(G, b, lam, x0=None, tol=I_KKT, it_max=10000):
    """The LASSO min_x ½xᵀGx − bᵀx + λ‖x‖₁ (G = ⟨φ_iφ_j⟩, b = ⟨φ_if*⟩) by feature-sign search (Lee, Battle, Raina, Ng 2007): on the current
    support and signs the exact Newton step, an exact line search over the step's zero crossings, then the inactive coordinate of largest
    |∂| added while it exceeds λ; warm-started from x0. Returns x, the largest KKT violation relative to λ (off the support |c_j|/λ − 1,
    on it |c_j − λ sign x_j|/λ; c = b − Gx) and the iterations."""
    n = len(b)
    x = np.zeros(n) if x0 is None else np.array(x0, float)
    obj = lambda z: 0.5 * z @ G @ z - b @ z + lam * np.abs(z).sum()
    it, full = 0, False
    while it < it_max:
        it += 1
        gr = G @ x - b
        S = np.flatnonzero(x)
        th = np.sign(x)
        if full or not (len(S) and np.abs(gr[S] + lam * th[S]).max() > tol * lam):   # a full step solved the support's equations
            z = np.flatnonzero(x == 0)
            if len(z) == 0:
                break
            i = z[np.argmax(np.abs(gr[z]))]
            if abs(gr[i]) <= lam * (1 + tol):
                break
            th[i] = -np.sign(gr[i])
            S = np.r_[S, i]
        xn = np.zeros(n)
        xn[S] = np.linalg.solve(G[np.ix_(S, S)], b[S] - lam * th[S])
        dl = xn - x
        with np.errstate(divide='ignore', invalid='ignore'):
            tc = -x[S] / dl[S]
        ts = [1.0] + [t for t in tc if 0 < t < 1]
        bt = min(ts, key=lambda t: obj(x + t * dl))
        x = x + bt * dl
        full = bt == 1 and bool(np.all(np.sign(x[S]) == th[S]))
        if bt < 1:                                                     # the coordinate whose zero crossing the search stopped at
            x[S[np.argmin(np.where(np.isfinite(tc), np.abs(tc - bt), np.inf))]] = 0.0
    c = b - G @ x
    S = np.flatnonzero(x)
    off = np.ones(n, bool)
    off[S] = False
    kkt = max(float(np.abs(c[off]).max() / lam - 1) if off.any() else -1.0,
              float(np.abs(c[S] - lam * np.sign(x[S])).max() / lam) if len(S) else 0.0)
    return x, kkt, it


def blasso(K, cs, As, h=None):
    """The BLASSO on Part H's basis-pursuit grid (spacing h, the features added) at λ = ℓλ_max for ℓ in I_LAMS (largest first, each warm
    started from the previous), λ_max = max_j|⟨φ(c_j)f*⟩|: per ℓ the support's positions and coefficients, the KKT violation, the
    iterations."""
    cs = np.asarray(cs, float).reshape(len(As), -1)
    h = H_BP[cs.shape[1] - 1] if h is None else h
    pos = bp_grid(K, cs, h)
    P = K.phi(pos)
    f = K.phi(cs) @ np.asarray(As, float)
    G, b = P.T @ P / K.N, P.T @ f / K.N
    lmax = float(np.abs(b).max())
    x, sols = None, []
    for l in I_LAMS:
        x, kkt, it = lasso_fs(G, b, l * lmax, x)
        S = np.flatnonzero(x)
        sols.append(dict(l=l, pos=pos[S].copy(), x=x[S].copy(), kkt=kkt, it=it))
    return dict(lmax=lmax, sols=sols)


def blasso_verdict(sol, cs, As, tol):
    """A BLASSO solution against the features: every nonzero within tol of a feature, each feature's cluster (its nonzeros within tol)
    nonempty with the sum of A_k's sign — recovered —; the support exactly the features; the features lost (empty cluster); the largest
    relative amplitude error max_k|Σ_cluster x − A_k|/|A_k| over the clusters."""
    cs = np.asarray(cs, float).reshape(len(As), -1)
    pos = np.asarray(sol['pos'], float).reshape(len(sol['x']), cs.shape[1])
    x = np.asarray(sol['x'], float)
    As = np.asarray(As, float)
    d = np.linalg.norm(pos[:, None, :] - cs[None, :, :], axis=2) if len(x) else np.zeros((0, len(As)))
    near = d.argmin(1) if len(x) else np.zeros(0, int)
    far = int(np.sum(d.min(1) > tol)) if len(x) else 0
    amp = [float(x[(near == k) & (d.min(1) <= tol)].sum()) if len(x) else 0.0 for k in range(len(As))]
    cnt = [int(np.sum((near == k) & (d.min(1) <= tol))) if len(x) else 0 for k in range(len(As))]
    lost = sum(c == 0 for c in cnt)
    signs = all(c > 0 and np.sign(a) == np.sign(A) for c, a, A in zip(cnt, amp, As))
    exact = len(x) == len(As) and bool(np.all(d.min(1) < 1e-12))
    err = float(max(abs(a - A) / abs(A) for a, A in zip(amp, As)))
    return dict(rec=bool(far == 0 and signs), exact=exact, lost=int(lost), far=far, err=err, n=len(x))


def i_job(args):
    """I: one one-dimensional target's BLASSO."""
    src, cs, As = args
    warnings.filterwarnings('ignore')
    L, K0, c0 = TM.setup()
    t0 = time.time()
    r = blasso(_k1(K0), np.asarray(cs, float).reshape(len(As), 1), np.asarray(As, float))
    return dict(src=src, **r, secs=time.time() - t0)


def i_jobs(A, B, D):
    jobs = [(('A', i), r['cs'], r['As']) for i, r in enumerate(A) if not r['skip']]
    jobs += [(('B', i), r['cs'], r['As']) for i, r in enumerate(B) if not r['skip']]
    jobs += [(('D', i), r['cs'], r['As']) for i, r in enumerate(D) if not r['skip']]
    return jobs


def jets(K, cs):
    """φ, ∇_cφ and the second derivatives ∂²_{ab}φ (a ≤ b) at the features, on the samples: (N, K(1 + n + n(n + 1)/2))."""
    cs = np.asarray(cs, float).reshape(len(cs), -1)
    n = cs.shape[1]
    P = K.phi(cs)
    dx = (K.xs[:, None, :] - cs[None, :, :]) / K.sig ** 2
    cols = [P] + [P * dx[:, :, a] for a in range(n)]
    cols += [P * (dx[:, :, a] * dx[:, :, b] - (a == b) / K.sig ** 2) for a in range(n) for b in range(a, n)]
    return np.concatenate(cols, 1)


def bp_jet_certificate(K, cs, As, h=None):
    """Basis pursuit on the grid of Part H (spacing h, the features replacing their grid neighbours) certified in a subspace: the least
    t such that a p in the span of the features' second-order jets (`jets`, orthonormalised, singular values below 10⁻¹³ of the largest
    dropped) has ⟨pφ(c*_k)⟩ = sign A_k and |⟨pφ(c)⟩| ≤ t at every other grid point (a linear program, HiGHS). Any such p is a certificate,
    so t_J ≥ the least certificate's t of Part H and t_J < 1 proves that the features are the unique minimum-ℓ¹ fit on the grid. Also
    η_V's largest |η_V| on the same grid points (η_V is a p of that span: t_J ≤ it)."""
    from scipy.optimize import linprog
    cs = np.asarray(cs, float).reshape(len(As), -1)
    n = cs.shape[1]
    h = H_BP[n - 1] if h is None else h
    g = bp_grid(K, cs, h)[:-len(As)]
    Pg = np.concatenate([K.phi(g[i:i + 4000]) for i in range(0, len(g), 4000)], axis=1)
    U, sv, _ = np.linalg.svd(jets(K, cs), full_matrices=False)
    Ub = U[:, sv > sv[0] * 1e-13]
    Qg, Qs = Pg.T @ Ub, K.phi(cs).T @ Ub
    M = Ub.shape[1]
    sol = linprog(np.r_[np.zeros(M), 1.0], A_ub=np.r_[np.c_[Qg, -np.ones(len(g))], np.c_[-Qg, -np.ones(len(g))]],
                  b_ub=np.zeros(2 * len(g)), A_eq=np.c_[Qs, np.zeros(len(As))], b_eq=np.sign(As), bounds=[(None, None)] * M + [(0, None)],
                  method='highs')
    Gm = np.c_[K.phi(cs), K.dphi(cs).reshape(K.N, -1)]
    q = Gm @ np.linalg.solve(Gm.T @ Gm, np.r_[np.sign(As), np.zeros(len(As) * n)])
    return dict(t_jet=float(sol.x[-1]) if sol.status == 0 else np.nan, dim=int(M), t_etav=float(np.abs(Pg.T @ q).max()))


def j_job(args):
    """J: one target's jet certificate (one dimension or the plane)."""
    src, cs, As, dim = args
    warnings.filterwarnings('ignore')
    K = _k1(TM.setup()[1]) if dim == 1 else TN.full()[0]
    t0 = time.time()
    r = bp_jet_certificate(K, np.asarray(cs, float).reshape(len(As), dim), np.asarray(As, float))
    return dict(src=src, **r, secs=time.time() - t0)


def j_jobs(A, B, C, D):
    jobs = [(('A', i), r['cs'], r['As'], 1) for i, r in enumerate(A) if not r['skip']]
    jobs += [(('B', i), r['cs'], r['As'], 1) for i, r in enumerate(B) if not r['skip']]
    jobs += [(('D', i), r['cs'], r['As'], 1) for i, r in enumerate(D) if not r['skip']]
    jobs += [(('C', i), r['cs'], r['As'], 2) for i, r in enumerate(C) if not r['skip']]
    return jobs


def run_chunks(fn, jobs, done, save, chunk=None):
    """Resumable: the jobs whose src is not in `done` (a list of rows, extended in place), in chunks of `chunk`, `save()` after each."""
    from multiprocessing import Pool
    have = {tuple(r['src']) for r in done}
    todo = [j for j in jobs if tuple(j[0]) not in have]
    chunk = chunk or I_CHUNK
    for i in range(0, len(todo), chunk):
        with Pool(N_WORKERS) as pool:
            done.extend(pool.map(fn, todo[i:i + chunk], chunksize=4))
        save()
    order = {tuple(j[0]): n for n, j in enumerate(jobs)}
    done.sort(key=lambda r: order.get(tuple(r['src']), len(order)))
    return done


def families(A, B, C, D):
    """The rows of A–D by family, as in the paper's table: A, B, C, D(a), D(b), D(c) (the skipped rows dropped); keys (part, index)."""
    fam = dict(A=[('A', i) for i, r in enumerate(A) if not r['skip']], B=[('B', i) for i, r in enumerate(B) if not r['skip']],
               C=[('C', i) for i, r in enumerate(C) if not r['skip']])
    for t in 'abc':
        fam[f'D({t})'] = [('D', i) for i, r in enumerate(D) if r['tag'][0] == t and not r['skip']]
    return fam


# ---- K (after round 149): Part I's failures by type, from the stored solutions (no re-solve); Part G's entry sequences -------------------

def k_classify(sol, cs, As, tol, rd):
    """One stored BLASSO solution against the features, by type. Each nonzero belongs to its nearest feature; within `tol` it is in that
    feature's cluster (as in `blasso_verdict`), beyond `tol` but within `rd` it is *displaced*, beyond `rd` *spurious*. A feature is *lost*
    when no nonzero within `rd` has it as its nearest feature; a cluster has the *wrong sign* when its sum's sign is not A_k's. The class of
    a target that `blasso_verdict` does not recover is the first of lost, spurious, displaced, wrong sign that applies (the four cover every
    failure); 'rec' otherwise. Also: the empty clusters within `tol` (`blasso_verdict`'s 'lost'), the largest displaced and the smallest
    spurious distance, and the lost features' |A_k|/max|A|."""
    cs = np.asarray(cs, float).reshape(len(As), -1)
    As = np.asarray(As, float)
    x = np.asarray(sol['x'], float)
    pos = np.asarray(sol['pos'], float).reshape(len(x), cs.shape[1])
    d = np.linalg.norm(pos[:, None, :] - cs[None, :, :], axis=2) if len(x) else np.zeros((0, len(As)))
    near, dm = (d.argmin(1), d.min(1)) if len(x) else (np.zeros(0, int), np.zeros(0))
    lost = [k for k in range(len(As)) if not np.any((near == k) & (dm <= rd))]
    disp, spur = dm[(dm > tol) & (dm <= rd)], dm[dm > rd]
    cnt = [int(np.sum((near == k) & (dm <= tol))) for k in range(len(As))]
    sign = any(c > 0 and np.sign(x[(near == k) & (dm <= tol)].sum()) != np.sign(A) for k, (c, A) in enumerate(zip(cnt, As)))
    rec = blasso_verdict(sol, cs, As, tol)['rec']
    cls = 'rec' if rec else 'lost' if lost else 'spurious' if len(spur) else 'displaced' if len(disp) else 'sign' if sign else '?'
    return dict(cls=cls, lost_tol=sum(c == 0 for c in cnt), lost=len(lost), disp=int(len(disp)), spur=int(len(spur)), sign=bool(sign),
                ddisp=float(disp.max()) if len(disp) else np.nan, dspur=float(spur.min()) if len(spur) else np.nan,
                lost_rel=[float(abs(As[k]) / np.abs(As).max()) for k in lost])


def part_K(A, B, D, G, I):
    """K: every stored Part I solution classified by `k_classify` (tolerance I_TOL spacings as in Part I, the displacement radius K_RD·σ);
    and Part G's audit as entry sequences: the unpolished flow's entries against the polished path's (A and D(b) from their own rows, B,
    D(a), D(c) from Part G's), and in the plane, where no sequence of the unpolished flow is stored, the outcome classes."""
    P = dict(A=A, B=B, D=D)
    tol, rd = I_TOL * H_BP[0], K_RD * TM.setup()[1].sig
    rows = []
    for r in I:
        t = P[r['src'][0]][r['src'][1]]
        rows.append(dict(src=tuple(r['src']), types=[k_classify(s, t['cs'], t['As'], tol, rd) for s in r['sols']]))
    seq = [(('A', i), list(r['order']), list(r['order_flow'])) for i, r in enumerate(A) if not r['skip'] and 'order_flow' in r]
    seq += [(('D', i), list(r['order']), list(r['order_flow'])) for i, r in enumerate(D) if not r['skip'] and 'order_flow' in r]
    seq += [(tuple(g['tag']), list(P[g['tag'][0]][g['tag'][1]]['order']), list(g['order_flow'])) for g in G['rows'] if not g.get('skip')]
    return dict(tol=tol, rd=rd, rows=rows, seq=[dict(src=s, same=o == f) for s, o, f in seq],
                plane=[dict(src=tuple(g['src']), same=g['out'] == g['out_flow']) for g in G['plane'] if not g.get('skip')])


# ---- L and M (after round 149): basis pursuit's least certificate in the plane exactly; the BLASSO in the plane -----------------------

def _kron_basis(K, h, cut=1e-13):
    """L: the plane's grid columns on the tensor samples are G1 ⊗ G1 (G1: the one-dimensional Gaussians at the grid's abscissae on the
    samples' abscissae); with G1 = U1 S1 V1ᵀ their span's orthonormal basis is U1 ⊗ U1 with the singular values s_is_j. Returns
    (U1, s1, V1, ti, tj), the pairs (ti, tj) those with s_is_j > 10⁻¹³ s_0² (Part H's truncation)."""
    assert np.array_equal(K.x1, K.x2)
    g1 = np.arange(H_BOX[1][0], H_BOX[1][1] + 1e-9, h)
    U1, s1, V1t = np.linalg.svd(np.exp(-(K.x1[:, None] - g1[None, :]) ** 2 / (2 * K.sig ** 2)), full_matrices=False)
    ti, tj = np.nonzero(np.outer(s1, s1) > s1[0] ** 2 * cut)
    return U1, s1, V1t.T, ti, tj


def _dual_lower(BT, Bw, solve, keep, l1, l2, nu, Cf, s):
    """L: a dual iterate made exactly feasible and its bound. For min t subject to Cw = s and −t ≤ Bw ≤ t on the kept grid, (l1, l2 ≥ 0,
    ν) is dual feasible when B_keepᵀ(l1 − l2) + Cᵀν = 0 and Σ(l1 + l2) = 1, and then −s·ν ≤ the least t. The residual r of the first
    condition is removed by δ = B_keep(B_keepᵀB_keep)⁻¹(−r) (`solve` applies (B_keepᵀB_keep)⁻¹; BT, Bw apply Bᵀ, B), δ's positive part
    added to l1 and its negative part to l2 (both stay ≥ 0), and the mass is renormalised: returns −s·ν/mass and the corrected
    (unnormalised) l1, l2."""
    r = BT(np.where(keep, l1 - l2, 0.0)) + Cf.T @ nu
    dlt = Bw(solve(-r)) * keep
    a, b = l1 + np.maximum(dlt, 0), l2 + np.maximum(-dlt, 0)
    return float(-(s @ nu) / (a.sum() + b.sum())), a, b


def bp_certificate_kron(K, cs, As, h=None, tol=L_TOL, it_max=L_IT, cut=1e-13):
    """L: Part H's least certificate in the plane, the full linear program: p = Σ_{(i,j)} w_ij u_i ⊗ u_j/(s_is_j) over the pairs kept at
    the truncation `cut` (`_kron_basis`; 10⁻¹³ is Part H's), so that the values at the grid points are η = (V1 ⊗ V1)w (orthonormal
    columns); minimise t subject to η(c*_k) = sign A_k and |η| ≤ t at every grid point farther than half a spacing from the features. A
    primal–dual interior-point method (Mehrotra's predictor–corrector; the normal matrix (V1 ⊗ V1)ᵀD(V1 ⊗ V1) formed through the tensor
    factors, D the diagonal weight on the grid). Each iterate gives two certified bounds: max|η(w)| over the grid (w satisfies the
    equalities to `re`: an upper bound), and the dual iterate made exactly feasible (the stationarity residual removed by a grid vector
    δ = B_keep(B_keepᵀB_keep)⁻¹(−r), its positive and negative parts added to the two multipliers, rescaled to unit mass: a lower bound by
    weak duality). Stops when the best bracket is narrower than `tol`, or has not narrowed for L_STALL iterations, or the normal matrix
    is not numerically definite. Returns t (the best upper bound), t_low, the iterations, the dimension of p, the constraint pairs, the
    grid points within 10⁻⁷ of t, and the certificate p evaluated directly on the samples: its sup on the grid t_dir and its largest
    |⟨pφ(c*_k)⟩ − sign A_k| feq (the rounding of p's large components)."""
    from scipy.linalg import cho_factor, cho_solve
    cs = np.asarray(cs, float).reshape(len(As), 2)
    s = np.sign(np.asarray(As, float))
    h = H_BP[1] if h is None else h
    U1, s1, V1, ti, tj = _kron_basis(K, h, cut)
    n1, r1 = V1.shape
    M = len(ti)
    g = _grid(2, h)
    keep = np.min(np.linalg.norm(g[:, None, :] - cs[None, :, :], axis=2), axis=1) > 0.5 * h
    km = keep.astype(float)
    n = int(keep.sum())
    W2 = (V1[:, :, None] * V1[:, None, :]).reshape(n1, r1 * r1)
    idx = ti * r1 + tj

    def Bw(w):
        Wm = np.zeros((r1, r1))
        Wm[ti, tj] = w
        return (V1 @ Wm @ V1.T).ravel()

    def BT(q):
        return (V1.T @ q.reshape(n1, n1) @ V1)[ti, tj]

    def BDB(d):
        F = (W2.T @ (d.reshape(n1, n1) @ W2)).reshape(r1, r1, r1, r1).transpose(0, 2, 1, 3).reshape(r1 * r1, r1 * r1)
        return F[np.ix_(idx, idx)]

    fac = lambda x: np.exp(-(K.x1 - x) ** 2 / (2 * K.sig ** 2)) @ U1 / s1
    Cf = np.array([fac(c[0])[ti] * fac(c[1])[tj] for c in cs])
    Ef = np.c_[Cf, np.zeros(len(s))]
    w = np.linalg.lstsq(Cf, s, rcond=None)[0]
    eta = Bw(w)
    t = float(np.abs(eta[keep]).max()) * 1.1 + 0.1
    u1, u2 = (t - eta) * km, (t + eta) * km
    l1, l2 = km / (2 * n), km / (2 * n)
    nu = np.zeros(len(s))
    safe = lambda a: np.where(keep, a, 1.0)
    amax = lambda v, dv: min(1.0, float((-v[dv < 0] / dv[dv < 0]).min())) if np.any(dv < 0) else 1.0
    cfK = cho_factor(BDB(km))                                          # B_keepᵀB_keep, for the dual's feasibility correction
    best = dict(up=np.inf, low=-np.inf, w=w, it=0)
    t0 = time.time()
    for it in range(1, it_max + 1):
        eta = Bw(w)
        re = Cf @ w - s
        up = float(np.abs(eta[keep]).max())                            # w is a certificate (re at rounding): an upper bound
        rd = np.r_[BT(np.where(keep, l1 - l2, 0.0)) + Cf.T @ nu, 1 - l1.sum() - l2.sum()]
        low = _dual_lower(BT, Bw, lambda v: cho_solve(cfK, v), keep, l1, l2, nu, Cf, s)[0]
        if up < best['up']:
            best.update(up=up, w=w.copy(), re=float(np.abs(re).max()))
        if low > best['low']:
            best.update(low=low)
        if best['up'] - best['low'] < best.get('width', np.inf):
            best.update(width=best['up'] - best['low'], it=it)
        if best['width'] <= tol or it - best['it'] >= L_STALL:
            break
        rp1, rp2 = (eta - t + u1) * km, (-eta - t + u2) * km
        mu = (l1 @ u1 + l2 @ u2) / (2 * n)
        D1, D2 = np.where(keep, l1 / safe(u1), 0.0), np.where(keep, l2 / safe(u2), 0.0)
        H = np.zeros((M + 1, M + 1))
        H[:M, :M] = BDB(D1 + D2)
        H[:M, M] = H[M, :M] = BT(D2 - D1)
        H[M, M] = (D1 + D2).sum()
        try:
            cf = cho_factor(H + np.eye(M + 1) * (L_REG * np.abs(np.diag(H)).max()))
        except np.linalg.LinAlgError:                                  # the normal matrix lost definiteness in rounding: stop here
            break
        HiE = cho_solve(cf, Ef.T)
        Sch = Ef @ HiE

        def step(rc1, rc2):
            q1 = np.where(keep, (l1 * rp1 - rc1) / safe(u1), 0.0)
            q2 = np.where(keep, (l2 * rp2 - rc2) / safe(u2), 0.0)
            Hr = cho_solve(cf, -rd - np.r_[BT(q1 - q2), -(q1.sum() + q2.sum())])
            dnu = np.linalg.solve(Sch, Ef @ Hr + re)
            dx = Hr - HiE @ dnu
            de = Bw(dx[:M])
            du1, du2 = (-rp1 - de + dx[M]) * km, (-rp2 + de + dx[M]) * km
            return (dx, dnu, du1, du2, np.where(keep, (-rc1 - l1 * du1) / safe(u1), 0.0),
                    np.where(keep, (-rc2 - l2 * du2) / safe(u2), 0.0))

        dx, dnu, du1, du2, dl1, dl2 = step(l1 * u1, l2 * u2)          # the affine (predictor) step
        ap = min(amax(u1[keep], du1[keep]), amax(u2[keep], du2[keep]))
        ad = min(amax(l1[keep], dl1[keep]), amax(l2[keep], dl2[keep]))
        sig = (((l1 + ad * dl1) @ (u1 + ap * du1) + (l2 + ad * dl2) @ (u2 + ap * du2)) / (2 * n) / mu) ** 3
        dx, dnu, du1, du2, dl1, dl2 = step(l1 * u1 + du1 * dl1 - sig * mu * km, l2 * u2 + du2 * dl2 - sig * mu * km)
        ap = 0.995 * min(amax(u1[keep], du1[keep]), amax(u2[keep], du2[keep]))
        ad = 0.995 * min(amax(l1[keep], dl1[keep]), amax(l2[keep], dl2[keep]))
        w, t, u1, u2 = w + ap * dx[:M], t + ap * dx[M], u1 + ap * du1, u2 + ap * du2
        l1, l2, nu = l1 + ad * dl1, l2 + ad * dl2, nu + ad * dnu
    secs = time.time() - t0
    w = best['w']
    a = np.abs(Bw(w)[keep])
    Wm = np.zeros((r1, r1))
    Wm[ti, tj] = w / (s1[ti] * s1[tj])
    p = (U1 @ Wm @ U1.T).ravel()                                        # the certificate on the samples (index a·n + b)
    return dict(t=best['up'], t_low=best['low'], it=it, dim=M, n_con=n, nact=int(np.sum(a >= best['up'] - 1e-7)), re=best['re'],
                t_dir=float(np.abs(_corr(K, g[keep], p)).max()), feq=float(np.abs(K.phi(cs).T @ p - s).max()), secs=secs)


def l_job(args):
    """L: one plane target's least certificate at each truncation of L_CUTS."""
    src, cs, As = args
    warnings.filterwarnings('ignore')
    K = TN.full()[0]
    return dict(src=src, cuts=[dict(cut=c, **bp_certificate_kron(K, cs, np.asarray(As, float), cut=c)) for c in L_CUTS])


def l_jobs(C):
    return [(('C', i), r['cs'], r['As']) for i, r in enumerate(C) if not r['skip']]


class KronGram:
    """M: the Gram matrix ⟨φ_cφ_d⟩ (the samples' mean) of Gaussian columns on a tensor grid of samples, from the columns' one-dimensional
    factors F1, F2 (n1 × m, n2 × m): G = (F1ᵀF1) ∘ (F2ᵀF2)/N, never formed. `G @ x` uses only x's nonzeros, `G[np.ix_(I, J)]` gives a
    block; enough for `lasso_fs`."""
    __array_ufunc__ = None                                              # z @ G with an ndarray z goes to __rmatmul__

    def __init__(self, F1, F2, N):
        self.F1, self.F2, self.N = F1, F2, N
        self.shape = (F1.shape[1], F1.shape[1])

    def __matmul__(self, x):
        S = np.flatnonzero(x)
        R = self.F1[:, S] @ (np.asarray(x)[S, None] * self.F2[:, S].T)
        return np.einsum('am,ab,bm->m', self.F1, R, self.F2, optimize=True) / self.N

    def __rmatmul__(self, z):                                          # G is symmetric
        return self @ z

    def __getitem__(self, ij):
        I, J = (np.asarray(a).ravel() for a in ij)
        return (self.F1[:, I].T @ self.F1[:, J]) * (self.F2[:, I].T @ self.F2[:, J]) / self.N


def blasso_kron(K, cs, As, h=None):
    """M: `blasso` in the plane (K a `BumpsGrid2`): the same problem, grid (Part H's basis-pursuit grid of spacing h with the features
    added), λ = ℓλ_max for ℓ in I_LAMS with warm starts and feature-sign search, the Gram matrix as a `KronGram`."""
    cs = np.asarray(cs, float).reshape(len(As), 2)
    h = H_BP[1] if h is None else h
    pos = bp_grid(K, cs, h)
    G = KronGram(np.exp(-(K.x1[:, None] - pos[None, :, 0]) ** 2 / (2 * K.sig ** 2)),
                 np.exp(-(K.x2[:, None] - pos[None, :, 1]) ** 2 / (2 * K.sig ** 2)), K.N)
    xf = np.zeros(len(pos))
    xf[-len(As):] = As                                                  # f* = Σ_k A_kφ(c*_k): b = ⟨φ_jf*⟩ = (G x_f)_j
    b = G @ xf
    lmax = float(np.abs(b).max())
    x, sols = None, []
    for l in I_LAMS:
        x, kkt, it = lasso_fs(G, b, l * lmax, x)
        S = np.flatnonzero(x)
        sols.append(dict(l=l, pos=pos[S].copy(), x=x[S].copy(), kkt=kkt, it=it))
    return dict(lmax=lmax, sols=sols)


def m_job(args):
    """M: one plane target's BLASSO."""
    src, cs, As = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    r = blasso_kron(TN.full()[0], cs, np.asarray(As, float))
    return dict(src=src, **r, secs=time.time() - t0)


def report(A=None, B=None, C=None, D=None, E_=None, F=None, G=None, H=None, I=None, J=None, K=None, L=None, M=None):
    PA = dict(A=A, B=B, C=C, D=D)                                      # G, H: the Parts' rows (the names A–D are reused below)
    L_ = ["# Theorem D‴.55: the transit's certificate (`exp_transit_certificate.py`)", "",
          "The sliding inverse scale space (the free law's limit from small masses, moving centres) on random targets of three Gaussian "
          "features of the kernels' width; the certificate η(c_j) = Σ_k sign(A_k)⟨φ(c*_k)φ(c_j)⟩/⟨φ(c*_k)φ(c_{e_k})⟩ on the starting points "
          "(the samples' inner product, e_k the starting point of largest correlation with feature k — the nearest in the continuum —, E "
          "their set), its margin m = 1 − max_{j∉E}|η(c_j)| and "
          "the cross-talk q = max_k Σ_{l≠k}G((Δ_kl − 2R)₊)/G(R). The path's outcome: exact with E (its active set E, the residual below "
          f"{RMS_TOL:g}, no exits), exact with other kernels (three kernels, the residual below {RMS_TOL:g}), or not exact. The paths by "
          "`sliding_iss_nd` (in one dimension on the 1-D samples, with the exact Jacobian), the block polished by its least-squares fit "
          "unless a kernel below 10⁻² of the block's largest coefficient would land beyond one kernel width of both its start and its "
          "correlation's ascent's end (a rejected polish; the bound does not check single-kernel blocks or the block's largest kernel); the "
          "skeletons recomputed from the targets.", ""]

    def table(rows, key, label):
        out = [f"| {label} | targets | m > 0: exact with E, with others, not | m ≤ 0: exact with E, with others, not | the median q | "
               "the disagreements' m (q) |", "|---|---|---|---|---|---|"]
        for v in sorted(set(r[key] for r in rows)):
            rs, pos, neg, cnt = confusion([r for r in rows if r[key] == v])
            bad = [r for r in rs if (r['m'] > 0) != (r['out'] == 'nearest')]
            out.append(f"| {v:g} | {len(rs)} | {cnt(pos, 'nearest')}, {cnt(pos, 'other')}, {cnt(pos, 'not')} | {cnt(neg, 'nearest')}, "
                       f"{cnt(neg, 'other')}, {cnt(neg, 'not')} | {np.median([r['q'] for r in rs]):.3f} | "
                       + (', '.join(f"{r['m']:+.4f} ({r['q']:.2f})" for r in bad) or '—') + " |")
        return out
    if A is not None:
        L_ += ["## A. One dimension", "", f"{N_A} targets per minimal separation (seed {SEED_A}; σ = 0.12, the grid's spacing 0.1); "
               "amplitudes ±U(0.4, 1):", ""] + table(A, 'dmin', 'Δ ≥')
        rs, pos, neg, cnt = confusion(A)
        bad = [r for r in rs if (r['m'] > 0) != (r['out'] == 'nearest')]
        L_ += ["", f"Over all {len(rs)} targets the certificate's sign predicts exactness with E on {len(rs) - len(bad)}; the "
               f"disagreements' margins |m| ≤ {max([abs(r['m']) for r in bad], default=0):.4f}, at most "
               f"{max([abs(r['m']) / r['q'] for r in bad], default=0):.4f} times their cross-talk q. " + band(rs),
               "", "Exactness with E predicted by the sign of each margin (the relative margin min_k s_kη(c_{e_k}) − max_{j∉E}|η(c_j)|; the "
               "single-feature margin min_k(1 − max_{j≠e_k}⟨φ_kφ_j⟩/⟨φ_kφ_{e_k}⟩)), " + predictors(fill_margins(A)) + ".", ""]
    if B is not None:
        cfgs = sorted(set(r['cfg'] for r in B))
        same, n_ok, flip = 0, 0, []
        for c in cfgs:
            rs = [r for r in B if r['cfg'] == c and not r['skip']]
            outs = [r['out'] for r in sorted(rs, key=lambda r: r['r'])]
            if len(set(outs)) == 1:
                same += 1
            else:
                flip.append((c, outs, rs[0]['m']))
            n_ok += int(rs[0]['m'] > 0)
        orders = sum(len(set(tuple(r['order']) for r in B if r['cfg'] == c and not r['skip'])) > 1 for c in cfgs)
        ratio = {c: min(abs(r['m']) / r['q'] for r in B if r['cfg'] == c and not r['skip']) for c in cfgs}
        rng_r = {rr: [float(np.abs(r['As']).max() / np.abs(r['As']).min()) for r in B if r['r'] == rr and not r['skip']] for rr in R_B}
        keep = [c for c in cfgs if ratio[c] > 0.8 and c not in [f[0] for f in flip]]
        L_ += ["## B. The amplitudes' range", "",
               f"{N_B} configurations (seed {SEED_B}; positions at least {DMIN_B:g} apart, signs) with |A_k| = 10^{{−r·u_k}}, r = "
               f"{', '.join(str(r) for r in R_B)} (the amplitudes' range up to 10^{max(R_B)}); the certificate depends on the positions and "
               f"signs only ({n_ok} of {len(cfgs)} with m > 0). The outcome is the same at every range on {same} of {len(cfgs)} "
               f"configurations, while the order of the entries changes with the range on {orders}"
               + (f"; it changes on {len(flip)}: " + '; '.join(f"configuration {c} (m = {m:+.4f}): {', '.join(o)}" for c, o, m in flip)
                  if flip else '') + (f". The configurations whose outcome changes have |m|/q ≤ {max(ratio[c] for c, o, m in flip):.2f}; of "
                                      f"those with |m| > 0.8q, {len(keep)} of {sum(ratio[c] > 0.8 for c in cfgs)} keep their outcome at every "
                                      "range" if flip else '') + ". The ranges realised, max|A|/min|A|, median [max]: " +
               ', '.join(f"r = {rr}: {np.median(v):.3g} [{max(v):.3g}]" for rr, v in rng_r.items()) + "; the amplitudes' order is the same "
               "for every r > 0 (u_k fixed), so this Part varies the range, not the order (D(c) varies the order).",
               "", "Exactness with E predicted by the margins' signs, " + predictors(fill_margins(B)) + ".", ""]
    if C is not None:
        L_ += ["## C. The plane", "", f"{N_C} targets per minimal separation (seed {SEED_C}; σ = 0.15, the lattice's spacing 0.125; "
               "centres in [−0.7, 0.7]²), amplitudes ±U(0.4, 1):", ""] + table(C, 'dmin', 'Δ ≥')
        rs, pos, neg, cnt = confusion(C)
        bad = [r for r in rs if (r['m'] > 0) != (r['out'] == 'nearest')]
        L_ += ["", f"Over all {len(rs)} targets the certificate's sign predicts exactness with E on {len(rs) - len(bad)}. " + band(rs),
               "", "Exactness with E predicted by the margins' signs, " + predictors(fill_margins(C, dim=2)) + ".", ""]
    if D is not None:
        code = {'nearest': 'E', 'other': 'O', 'not': '×', 'timeout': 't'}
        D = fill_margins(D)                                            # the skeleton's margin (and older rows' relative margins)
        da = [r for r in D if r['tag'][0] == 'a']
        db = [r for r in D if r['tag'][0] == 'b']
        dc = [r for r in D if r['tag'][0] == 'c']
        L_ += ["## D. Designed families (referee round 43)", "",
               f"(a) Two features on starting points ({LAT_OFF:g} off) at −0.1 and −0.1 + D, the first's amplitude 1, the second's ±ρ_A: the "
               "margins (amplitude-free) and the outcome at each ρ_A (E exact with E, O exact with others, × not exact, t over "
               f"{T_MAX} s):", "",
               "| D | signs | m | m/q | relative | single-feature | " + ' | '.join(f"{x:.3g}" for x in RATIOS_D) + " |",
               '|---|---|---|---|---|---|' + '---|' * len(RATIOS_D)]
        for Dd in LAT_D:
            for sg in (1.0, -1.0):
                rs = sorted([r for r in da if r['tag'][1] == Dd and r['tag'][2] == sg], key=lambda r: r['tag'][3])
                r0 = rs[0]
                L_.append(f"| {Dd:g} | {'+ +' if sg > 0 else '+ −'} | {r0['m']:+.4f} | {r0['m'] / r0['q']:+.3f} | {r0['m_rel']:+.4f} | "
                          f"{r0['m_vor']:+.4f} | " + ' | '.join(code[r['out']] for r in rs) + " |")
        L_ += ["", f"(b) Two features at s·σ, equal amplitudes, the midpoint uniform over a grid cell ({N_OFF_D} per row): the outcomes and "
               "how often each margin's sign is right:", "",
               "| s | signs | E, O, ×, t | m > 0 | " + ' | '.join(f"{SHORT[k]} right" for k, _ in PREDICTORS) + " | always E right |",
               '|---' * (5 + len(PREDICTORS)) + '|']
        for s_ in SEPS_D:
            for sg in (1.0, -1.0):
                rs = [r for r in db if r['tag'][1] == s_ and r['tag'][2] == sg and not r['skip']]
                rr = [r for r in rs if r['out'] != 'timeout']
                ok = [r['out'] == 'nearest' for r in rr]
                L_.append(f"| {s_:g} | {'+ +' if sg > 0 else '+ −'} | " + ', '.join(str(sum(r['out'] == o for r in rs)) for o in code) +
                          f" | {sum(r['m'] > 0 for r in rs)} | " + ' | '.join(str(sum((r[k] > 0) == o for r, o in zip(rr, ok)))
                                                                         for k, _ in PREDICTORS) + f" | {sum(ok)} |")
        tg = sorted(set(r['tag'][1] for r in dc))
        vary = [t for t in tg if len(set(r['out'] for r in dc if r['tag'][1] == t)) > 1]
        L_ += ["", f"(c) {len(tg)} of A's targets with |m| > 0.8q, the amplitude magnitudes 1, {RANGE_D:g}^{{−1/2}}, {RANGE_D:g}^{{−1}} in all six "
               f"orders (the signs kept): the outcome changes with the order on {len(vary)} of them" +
               (': ' + '; '.join(f"m/q = {[r for r in dc if r['tag'][1] == t][0]['m'] / [r for r in dc if r['tag'][1] == t][0]['q']:+.2f}: "
                                 + ''.join(code[r['out']] for r in sorted([r for r in dc if r['tag'][1] == t], key=lambda r: r['tag'][2]))
                                 for t in vary) if vary else '') + ".", "",
               "Over D, exactness with E predicted by the margins' signs: (a) " + predictors(da) + "; (b) " + predictors(db) + "; (c) " +
               predictors(dc) + f". Paths over {T_MAX} s: {sum(r['out'] == 'timeout' for r in D)}; polishes rejected: "
               f"{sum(max(r.get('rej', 0), 0) for r in D)} on {sum(r.get('rej', 0) > 0 for r in D)} paths.", ""]
        bad_b = [r for r in db if not r['skip'] and r['out'] != 'nearest' and r['m_vor'] > M0_E]
        by_s = [(s_, sg, [r for r in bad_b if r['tag'][1] == s_ and r['tag'][2] == sg]) for s_ in SEPS_D for sg in (1.0, -1.0)]
        L_ += [f"D(b)'s failures that are not the second feature's own near ties (m₀ > {M0_E:g}), per row: " + ', '.join(
            f"{s_:g}σ {'+ +' if sg > 0 else '+ −'} {len(rs)}" + (" (m₀ = " + ', '.join(f"{r['m_vor']:.3f}" for r in rs) + ")"
                                                              if 0 < len(rs) <= 4 else '') for s_, sg, rs in by_s) + ".", ""]
        dd = [r for r in D if not r['skip'] and r['out'] != 'timeout']
        sent = []
        for C_ in BANDS:                                               # the band of A and C, tested on the designed families
            ob = [r for r in dd if abs(r['m']) > C_ * r['q']]
            wr = [r for r in ob if (r['m'] > 0) != (r['out'] == 'nearest')]
            sent.append(f"|m| > {C_:g}q: wrong on {len(wr)} of {len(ob)}" + (f" (|m|/q up to {max(abs(r['m']) / r['q'] for r in wr):.3f}, "
                                                                          f"{sum(r['m'] > 0 for r in wr)} of them with m > 0)" if wr else ""))
        L_ += ["Outside the band, the certificate's sign on the designed targets: " + '; '.join(sent) + ".", ""]
    if E_ is not None:
        g = [r for r in E_['grid'] if not r['skip']]
        sig = 0.12
        L_ += ["## E. The two-feature map", "",
               f"The first feature (amplitude 1) on the starting point −0.1 ({LAT_OFF:g} off), the second at the distance D = {D_E[0]:g} … "
               f"{D_E[-1]:g} (step {D_E[1] - D_E[0]:g}; σ = {sig:g}) with the amplitude ±10^ℓ, ℓ = {L_E[0]:g} … {L_E[-1]:g} (step "
               f"{L_E[1] - L_E[0]:g}): the fitted skeleton's outcome at every point (E: exact with the nearest kernels, ·: not), shown at "
               "every 0.06 in D and every ½ in ℓ:", ""]
        lv = [l for l in L_E if abs(l * 2 - round(l * 2)) < 1e-9]
        dv = [D for D in D_E if abs(round((D - D_E[0]) / 0.06) * 0.06 - (D - D_E[0])) < 1e-9]
        at = {(round(r['D'], 9), round(r['l'], 9), r['sg']): r for r in g}
        for sg in (1.0, -1.0):
            L_ += [f"{'Same' if sg > 0 else 'Opposite'} signs:", "", "| D/σ \\ ℓ | " + ' | '.join(f"{l:g}" for l in lv) + " |",
                   '|---|' + '---|' * len(lv)]
            for D in dv:
                L_.append(f"| {D / sig:.2f} | " + ' | '.join(('E' if at[(round(D, 9), round(l, 9), sg)]['fit_ok'] else '·')
                                                        if (round(D, 9), round(l, 9), sg) in at else ' ' for l in lv) + " |")
            L_.append("")
        thr, thr0 = {}, {}
        for sg in (1.0, -1.0):
            for l in L_E:
                col = sorted([r for r in g if r['sg'] == sg and abs(r['l'] - l) < 1e-9], key=lambda r: r['D'])
                bad = [r['D'] for r in col if not r['fit_ok']]
                thr[(sg, round(l, 9))] = (max(bad) + (D_E[1] - D_E[0])) if bad else col[0]['D']
                # without the second feature's own near ties (its single-feature margin m₀ ≤ M0_E: near a midpoint of the grid)
                bad0 = [r['D'] for r in col if not r['fit_ok'] and r['m_vor'] > M0_E]
                thr0[(sg, round(l, 9))] = (max(bad0) + (D_E[1] - D_E[0])) if bad0 else col[0]['D']
        show = [-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0]
        far = max(r['D'] for r in g if not r['fit_ok'] and r['m_vor'] > M0_E)  # the farthest failure that is not a near tie of its own
        tie_D = sorted(set(round(r['D'], 3) for r in g if not r['fit_ok'] and r['D'] > far))
        L_ += ["The fitted skeleton's threshold D*(ℓ) on this grid (the first feature on a starting point, |ℓ| ≤ 3: not (iii)'s Δ₀) — every "
               "D of the grid at or beyond it exact with E —, in units of σ: " + '; '.join(
            f"{'same' if sg > 0 else 'opposite'} signs " + ', '.join(f"ℓ = {l:g}: {thr[(sg, l)] / sig:.2f}" for l in show) for sg in (1.0, -1.0))
               + f"; over all ℓ at most {max(thr.values()) / sig:.2f}σ. Beyond {far / sig:.2f}σ every failure is the second feature's own near "
               f"tie (its single-feature margin m₀ ≤ {M0_E:g}, near a midpoint of the grid: D = {', '.join(f'{d:g}' for d in tie_D)}); without "
               "them D*(ℓ) is " + '; '.join(f"{'same' if sg > 0 else 'opposite'} signs " + ', '.join(
                   f"ℓ = {l:g}: {thr0[(sg, l)] / sig:.2f}" for l in show) for sg in (1.0, -1.0)) + f", at most {max(thr0.values()) / sig:.2f}σ."]
        ends = []                                                      # (v): the outcome at both ends of the ratio
        for sg in (1.0, -1.0):
            same_lo = same_hi = differ = 0
            for D in D_E:
                col = {round(r['l'], 9): r['fit_ok'] for r in g if r['sg'] == sg and abs(r['D'] - D) < 1e-9}
                lo = {col[k] for k in col if k <= -2 + 1e-9}
                hi = {col[k] for k in col if k >= 2 - 1e-9}
                same_lo += len(lo) == 1
                same_hi += len(hi) == 1
                differ += len(lo) == 1 and len(hi) == 1 and lo != hi
            ends.append(f"{'same' if sg > 0 else 'opposite'} signs: constant for ℓ ≤ −2 at {same_lo} and for ℓ ≥ 2 at {same_hi} of {len(D_E)} "
                        f"separations, the two ends' outcomes different at {differ}")
        L_ += ["", "The two ends of the ratio ((v)): " + '; '.join(ends) + "."]
        agree = {k: sum((r[k] > 0) == r['fit_ok'] for r in g) for k, _ in PREDICTORS if k != 'm_fit' and all(k in r for r in g)}
        L_ += ["", f"The margins' signs against the fitted skeleton on the map's {len(g)} points: " + ', '.join(
            f"{lab} {agree[k]}" for k, lab in PREDICTORS if k != 'm_fit' and all(k in r for r in g)) +
               f"; exact with E on {sum(r['fit_ok'] for r in g)}; the linear skeleton's |m_sk| > 0.01 where its sign differs from the fitted "
               f"skeleton's at {sum((r['m_skel'] > 0) != r['fit_ok'] and abs(r['m_skel']) > 0.01 for r in g)} points (|m_sk| up to "
               f"{max([abs(r['m_skel']) for r in g if (r['m_skel'] > 0) != r['fit_ok']], default=0):.3f})."]
        # (ii)'s order: the linear skeleton's margin against the fitted one's along the ratio, the tail being the smaller feature
        slopes = []
        for sg in (1.0, -1.0):
            for D in D_E:
                if not any(abs(D / sig - x) < 1e-6 for x in (2.5, 3.0, 3.5, 4.0)):
                    continue
                for side in (-1, 1):
                    col = [r for r in g if r['sg'] == sg and abs(r['D'] - D) < 1e-9 and 1 - 1e-9 <= side * r['l'] <= 3 + 1e-9
                           and r['fit_ok'] and r['m_skel'] > 0 and abs(r['m_skel'] - r['m_fit']) > 1e-13]
                    if len(col) >= 6:
                        x = np.array([abs(r['l']) for r in col]) * np.log(10)
                        y = np.log([abs(r['m_skel'] - r['m_fit']) for r in col])
                        slopes.append((D / sig, sg, side, float(-np.polyfit(x, y, 1)[0])))
        if slopes:
            L_ += ["", "(ii)'s order: |m_sk − m_fit| against the ratio of the smaller feature's amplitude to the larger's, its log–log "
                   "slope over 10⁻³ … 10⁻¹ (1: first order in the tail, 2: second): " + ', '.join(
                       f"{d:.2f}σ {'+ +' if sg > 0 else '+ −'} {'second small' if side < 0 else 'first small'} {sl:.2f}"
                       for d, sg, side, sl in slopes) + "."]
        P = [r for r in E_['paths'] if not r['skip']]
        done = [r for r in P if r['out'] != 'timeout']
        key = lambda r: (round(r['tag'][1], 9), round(r['tag'][2], 9), r['tag'][3])
        ok = sum((r['out'] == 'nearest') == at[key(r)]['fit_ok'] for r in done if key(r) in at)
        L_ += ["", f"The path at {len(P)} sampled points (|ℓ| ≤ {LMAX_EV:g}; under {T_EV} s each): it agrees with the fitted skeleton at {ok} of "
               f"{len(done)} ({len(P) - len(done)} over the limit); exact with E at {sum(r['out'] == 'nearest' for r in done)}; polishes "
               f"rejected {sum(max(r.get('rej', 0), 0) for r in done)} on {sum(r.get('rej', 0) > 0 for r in done)} paths.", ""]
    if F is not None:
        code = {'nearest': 'E', 'other': 'O', 'not': '×', 'timeout': 't'}
        L_ += ["## F. A wide range", "", f"Three features at −s, 0, s (+{F_OFF:g}; s = {', '.join(f'{x:g}' for x in F_SEPS)}: "
               f"{', '.join(f'{x / 0.12:g}σ' for x in F_SEPS)}), the amplitudes {', '.join(f'{x:g}' for x in F_MAGS)} in all six orders, "
               "the signs + + + and + − +; exits only below 10⁻¹² of the largest coefficient: the outcome per order (E exact with E, O with "
               "others, × not, t over the limit), the polishes rejected, the skeletons' signs:", "",
               "| s | signs | outcomes (orders 012, 021, 102, 120, 201, 210) | polishes rejected | the skeleton right | the fitted skeleton right |",
               "|---|---|---|---|---|---|"]
        F = fill_margins(list(F))
        for s_ in F_SEPS:
            for sg in ((1.0, 1.0, 1.0), (1.0, -1.0, 1.0)):
                rs = sorted([r for r in F if r['tag'][1] == s_ and tuple(r['tag'][2]) == sg], key=lambda r: tuple(r['tag'][3]))
                rr = [r for r in rs if not r['skip'] and r['out'] != 'timeout']
                L_.append(f"| {s_ / 0.12:g}σ | {' '.join('+' if v > 0 else '−' for v in sg)} | " + ''.join(code[r['out']] for r in rs) +
                          f" | {sum(max(r.get('rej', 0), 0) for r in rs)} | {sum((r['m_skel'] > 0) == (r['out'] == 'nearest') for r in rr)} of "
                          f"{len(rr)} | {sum((r['m_fit'] > 0) == (r['out'] == 'nearest') for r in rr)} of {len(rr)} |")
        L_.append("")
    if G is not None and all(x is not None for x in PA.values()):
        src = PA
        flow = {}
        for r in G['rows']:
            flow[tuple(r['tag'])] = r
        for r in G['plane']:
            flow[tuple(r['src'])] = r
        for i, r in enumerate(PA['A']):
            if not r['skip'] and 'out_flow' in r:
                flow[('A', i)] = r
        for i, r in enumerate(PA['D']):
            if not r['skip'] and 'out_flow' in r:
                flow[('D', i)] = r
        L_ += ["## G. The polish's audit on every target of A–D (referee round 140)", "",
               "The unpolished flow (polish off; the blocks relaxed by the fast flow alone) against the polished path, per family (A and D(b) from "
               "their Parts, the rest computed here): its outcome the polished path's, its blocks relaxed, and the polished path's outcome as "
               "recomputed here the stored one:", "",
               "| family | targets | the flow's outcome the path's | blocks relaxed | the path recomputed as stored | polishes rejected (paths) |",
               "|---|---|---|---|---|---|"]
        tot = [0, 0, 0]
        for fam, keys in families(PA['A'], PA['B'], PA['C'], PA['D']).items():
            rs = [(k, flow[k]) for k in keys if k in flow]
            same = sum(x['out_flow'] == src[k[0]][k[1]]['out'] for k, x in rs)
            rel = sum(bool(x['flow_relaxed']) for k, x in rs)
            rec = sum(x['out'] == src[k[0]][k[1]]['out'] for k, x in rs)
            rej = [max(src[k[0]][k[1]].get('rej', 0), 0) for k in keys]
            tot = [tot[0] + len(rs), tot[1] + same, tot[2] + rel]
            L_.append(f"| {fam} | {len(rs)} of {len(keys)} | {same} | {rel} | {rec} | {sum(rej)} ({sum(x > 0 for x in rej)}) |")
        L_ += ["", f"Over A–D: the unpolished flow's outcome the polished path's on {tot[1]} of {tot[0]} targets, its blocks relaxed on {tot[2]}.", ""]
    if H is not None and all(x is not None for x in PA.values()):
        src = PA
        hb = {tuple(h['src']): h for h in H}
        L_ += ["## H. Baselines off the path: η_V, basis pursuit, continuous OMP (referee round 140)", "",
               f"Per family: the targets; the path's outcomes (exact with E, exact with other kernels — both fit the features exactly with {K_F} "
               f"kernels, at their centres —, not exact); η_V nondegenerate (|η_V| < 1 on the grid of spacing {H_ETA[0]:g} (one dimension) or "
               f"{H_ETA[1]:g} (the plane) over [{H_BOX[0][0]:g}, {H_BOX[0][1]:g}] or [{H_BOX[1][0]:g}, {H_BOX[1][1]:g}]² beyond {EX_R:g}σ of the "
               "features, and s_k∇²η_V(c*_k) negative definite); basis pursuit on the grid of spacing "
               f"{H_BP[0]:g} (one dimension only; not computed in the plane), its points within half a spacing of a feature replaced by the "
               "feature, recovering the features (its least certificate's sup-norm t < 1); continuous OMP exact "
               f"after {K_F} steps (the rms residual below {RMS_TOL:g}; its selection grid {H_OMP[0]:g} or {H_OMP[1]:g}, refined), and with the "
               "joint slide after each step; and the targets on which the sign of η_V's margin on the starting points off E, ν_V = 1 − "
               "max_{j∉E}|η_V(c_j)|, and of the certificate's m predict exactness with E:", "",
               "| family | targets | path: E, other, not | η_V nondegenerate | basis pursuit exact | OMP exact | OMP with slide exact | "
               "ν_V's sign right | m's sign right |", "|---|---|---|---|---|---|---|---|---|"]
        tot = np.zeros(9, int)
        n1 = 0
        for fam, keys in families(PA['A'], PA['B'], PA['C'], PA['D']).items():
            rs = [(src[k[0]][k[1]], hb[k]) for k in keys if k in hb and src[k[0]][k[1]]['out'] != 'timeout']
            ok = [r['out'] == 'nearest' for r, h in rs]
            row = [len(rs), sum(r['out'] == 'nearest' for r, h in rs), sum(r['out'] == 'other' for r, h in rs),
                   sum(r['out'] == 'not' for r, h in rs), sum(h['nondeg'] for r, h in rs), sum(h['t_bp'] < 1 for r, h in rs),
                   sum(h['omp_rms'] < RMS_TOL for r, h in rs), sum(h['somp_rms'] < RMS_TOL for r, h in rs),
                   sum((h['m_v'] > 0) == o for (r, h), o in zip(rs, ok))]
            mright = sum((r['m'] > 0) == o for (r, h), o in zip(rs, ok))
            tot += np.array(row)
            n1 += row[0] if fam != 'C' else 0
            L_.append(f"| {fam} | {row[0]} | {row[1]}, {row[2]}, {row[3]} | {row[4]} | {row[5] if fam != 'C' else '—'} | {row[6]} | {row[7]} | "
                      f"{row[8]} | {mright} |")
        allm = sum((src[k[0]][k[1]]['m'] > 0) == (src[k[0]][k[1]]['out'] == 'nearest') for keys in families(PA['A'], PA['B'], PA['C'], PA['D']).values() for k in keys
                   if k in hb and src[k[0]][k[1]]['out'] != 'timeout')
        L_.append(f"| all | {tot[0]} | {tot[1]}, {tot[2]}, {tot[3]} | {tot[4]} | {tot[5]} of {n1} (one dimension) | {tot[6]} | {tot[7]} | {tot[8]} | {allm} |")
        allh = [h for h in H]
        dm = [h['omp_dmax'] for h in allh]
        L_ += ["", f"Continuous OMP without the slide: the largest distance of a selected atom from its nearest feature, median "
               f"{np.median(dm):.3g}σ, range {min(dm):.3g}σ–{max(dm):.3g}σ; its rms residual after {K_F} steps median "
               f"{E.fmt(float(np.median([h['omp_rms'] for h in allh])))}; the largest |η_V| beyond 1σ of the features: median {np.median([h['sup1'] for h in allh]):.4f}, range "
               f"{min(h['sup1'] for h in allh):.4f}–{max(h['sup1'] for h in allh):.4f}; "
               f"basis pursuit's t (one dimension): median {np.nanmedian([h['t_bp'] for h in allh]):.4f}, the linear program not solved on "
               f"{sum(not np.isfinite(h['t_bp']) for h in allh if h['src'][0] != 'C')}. Where the path is not exact with E but fits the features with other kernels, η_V is "
               "nondegenerate on " + str(sum(hb[k]['nondeg'] for keys in families(PA['A'], PA['B'], PA['C'], PA['D']).values() for k in keys
                                              if k in hb and src[k[0]][k[1]]['out'] == 'other'))
               + " of " + str(sum(src[k[0]][k[1]]['out'] == 'other' for keys in families(PA['A'], PA['B'], PA['C'], PA['D']).values() for k in keys if k in hb))
               + ", and basis pursuit exact on "
               + str(sum(hb[k]['t_bp'] < 1 for keys in families(PA['A'], PA['B'], PA['C'], PA['D']).values() for k in keys if k in hb and src[k[0]][k[1]]['out'] == 'other'))
               + " of the " + str(sum(src[k[0]][k[1]]['out'] == 'other' for keys in families(PA['A'], PA['B'], PA['C'], PA['D']).values() for k in keys
                                      if k in hb and k[0] != 'C')) + " in one dimension.", ""]
    if I is not None and all(x is not None for x in PA.values()):
        src = PA
        ib = {tuple(r['src']): r for r in I}
        tol = I_TOL * H_BP[0]
        fam1 = {k: v for k, v in families(PA['A'], PA['B'], PA['C'], PA['D']).items() if k != 'C'}
        rows1 = [k for v in fam1.values() for k in v if k in ib]
        for k in rows1:                                                # the fitted skeleton's margin where the Parts' sections did not fill it
            if 'm_fit' not in src[k[0]][k[1]]:
                fill_margins([src[k[0]][k[1]]])
        V = {k: [blasso_verdict(sol, src[k[0]][k[1]]['cs'], src[k[0]][k[1]]['As'], tol) for sol in ib[k]['sols']] for k in rows1}
        kk = [s_['kkt'] for k in rows1 for s_ in ib[k]['sols']]
        L_ += ["## I. BLASSO at a finite regularisation (referee round 140)", "",
               "The BLASSO on Part H's basis-pursuit grid (spacing " f"{H_BP[0]:g}" " over " f"[{H_BOX[0][0]:g}, {H_BOX[0][1]:g}]" ", its points "
               "within half a spacing of a feature replaced by the feature), min_x ½⟨(Σ_jx_jφ(c_j) − f*)²⟩ + λ‖x‖₁ (the samples' mean, no noise), "
               "at λ = ℓλ_max, λ_max = max_j|⟨φ(c_j)f*⟩| the least λ at which x = 0 is optimal, ℓ = " + ', '.join(f"{l:g}" for l in I_LAMS) +
               f", on every one-dimensional target of A, B and D ({len(rows1)}); solved by feature-sign search (exact Newton steps on the "
               "support and signs with an exact line search over their zero crossings, warm started from the larger λ) to the KKT conditions: "
               f"the largest violation over all {len(kk)} solves {E.fmt(max(kk))} relative to λ. Recovered: every nonzero coefficient within "
               f"{I_TOL} spacings ({tol:g}) of a feature, and each feature's cluster nonempty with its sum of A_k's sign; the support exactly "
               "the features; features lost (an empty cluster); the amplitude error max_k|Σ_cluster x − A_k|/|A_k| on the recovered targets, "
               "and that error over ℓ; and how recovery splits by the path's outcome (exact with E, exact with other kernels):", "",
               "| family | ℓ | targets | recovered | support exactly the features | targets with a feature lost | amplitude error median [max] "
               "| error/ℓ median | recovered: path E, path other |", "|---|---|---|---|---|---|---|---|---|"]
        for fam, keys in list(fam1.items()) + [('all', rows1)]:
            ks = [k for k in keys if k in ib]
            for n_l, l in enumerate(I_LAMS):
                vs = [(k, V[k][n_l]) for k in ks]
                er = [v['err'] for k, v in vs if v['rec']]
                L_.append(f"| {fam} | {l:g} | {len(vs)} | {sum(v['rec'] for k, v in vs)} | {sum(v['exact'] for k, v in vs)} | "
                          f"{sum(v['lost'] > 0 for k, v in vs)} | " + (f"{E.fmt(float(np.median(er)))} [{E.fmt(max(er))}]" if er else '—') +
                          " | " + (f"{np.median(er) / l:.3g}" if er else '—') + f" | {sum(v['rec'] for k, v in vs if src[k[0]][k[1]]['out'] == 'nearest')} of "
                          f"{sum(src[k[0]][k[1]]['out'] == 'nearest' for k, v in vs)}, {sum(v['rec'] for k, v in vs if src[k[0]][k[1]]['out'] == 'other')} "
                          f"of {sum(src[k[0]][k[1]]['out'] == 'other' for k, v in vs)} |")
        hb = {tuple(h['src']): h for h in H} if H is not None else {}
        lo = len(I_LAMS) - 1
        nd = [k for k in rows1 if k in hb and hb[k]['nondeg']]
        dg = [k for k in rows1 if k in hb and not hb[k]['nondeg']]
        fit_ok = [k for k in rows1 if src[k[0]][k[1]]['m_fit'] > 0]
        L_ += ["", f"At ℓ = {I_LAMS[lo]:g}: recovered on {sum(V[k][lo]['rec'] for k in nd)} of the {len(nd)} targets with η_V nondegenerate "
               f"and on {sum(V[k][lo]['rec'] for k in dg)} of the {len(dg)} without; on {sum(V[k][lo]['rec'] for k in fit_ok)} of the "
               f"{len(fit_ok)} the fitted skeleton predicts exact with E and on {sum(V[k][lo]['rec'] for k in rows1 if k not in fit_ok)} of the "
               f"{len(rows1) - len(fit_ok)} it does not; the support's size beyond the number of features (the grid neighbours of Duval–Peyré's "
               "extended support), median [max] per ℓ: " + ', '.join(
                   f"{l:g}: {np.median([V[k][n_l]['n'] - len(src[k[0]][k[1]]['As']) for k in rows1]):g} "
                   f"[{max(V[k][n_l]['n'] - len(src[k[0]][k[1]]['As']) for k in rows1)}]"
                   for n_l, l in enumerate(I_LAMS)) +
               f"; feature-sign iterations per solve, median [max]: {np.median([s_['it'] for k in rows1 for s_ in ib[k]['sols']]):g} "
               f"[{max(s_['it'] for k in rows1 for s_ in ib[k]['sols'])}].", ""]
    if J is not None and all(x is not None for x in PA.values()):
        src = PA
        hb = {tuple(h['src']): h for h in H} if H is not None else {}
        jb = {tuple(r['src']): r for r in J}
        fam = families(PA['A'], PA['B'], PA['C'], PA['D'])
        K2 = TN.full()[0]
        g2 = _grid(2, H_BP[1])
        sv2 = np.linalg.svd(np.concatenate([K2.phi(g2[i:i + 4000]) for i in range(0, len(g2), 4000)], axis=1), compute_uv=False)
        L_ += ["## J. Basis pursuit in the plane by a certificate in the features' jets (referee round 140)", "",
               "Part H's linear program for basis pursuit's least certificate has, in the plane, a p of dimension "
               f"{int(np.sum(sv2 > sv2[0] * 1e-13))} (the span of the {len(g2)} grid columns on the {K2.N} samples, singular values below 10⁻¹³ "
               f"of the largest dropped) and {len(g2) - K_F} pairs of constraints; here p is restricted to the span of "
               "the features' second-order jets (φ, ∇_cφ, ∇²_cφ at the features: "
               "three dimensions per feature in one dimension, six in the plane), and the same linear program (HiGHS) is solved over every grid "
               f"point (spacing {H_BP[0]:g} or {H_BP[1]:g}, the features replacing their neighbours). Every such p is a certificate, so its t_J "
               "bounds Part H's t from above, and t_J < 1 proves that the features are the unique minimum-ℓ¹ fit on the grid; η_V lies in "
               "that span, and its largest |η_V| on the same points bounds t_J from above:", "",
               "| family | targets | t_J < 1 (basis pursuit exact, certified) | t_J median [max] | max\\|η_V\\| on the grid median [max] | "
               "Part H's t median | t_J − t median [range] | path: E, other |", "|---|---|---|---|---|---|---|---|"]
        allr = []
        for f_, keys in fam.items():
            rs = [(k, jb[k]) for k in keys if k in jb]
            allr += rs
            tj = [x['t_jet'] for k, x in rs]
            tv = [x['t_etav'] for k, x in rs]
            th = [hb[k]['t_bp'] for k, x in rs if k in hb and np.isfinite(hb[k]['t_bp'])]
            gap = [x['t_jet'] - hb[k]['t_bp'] for k, x in rs if k in hb and np.isfinite(hb[k]['t_bp'])]
            L_.append(f"| {f_} | {len(rs)} | {sum(t < 1 for t in tj)} | {np.median(tj):.5f} [{max(tj):.5f}] | {np.median(tv):.5f} [{max(tv):.5f}] | "
                      + (f"{np.median(th):.4f} | {np.median(gap):.4f} [{min(gap):.4f}, {max(gap):.4f}]" if gap else "— | —") +
                      f" | {sum(src[k[0]][k[1]]['out'] == 'nearest' for k, x in rs)}, {sum(src[k[0]][k[1]]['out'] == 'other' for k, x in rs)} |")
        one = [(k, x) for k, x in allr if k[0] != 'C' and k in hb and np.isfinite(hb[k]['t_bp'])]
        pl = [(k, x) for k, x in allr if k[0] == 'C']
        L_ += ["", f"In one dimension ({len(one)} targets) t_J < 1 exactly where Part H's t < 1 on "
               f"{sum((x['t_jet'] < 1) == (hb[k]['t_bp'] < 1) for k, x in one)}, and t_J ≥ t − 10⁻⁹ on {sum(x['t_jet'] >= hb[k]['t_bp'] - 1e-9 for k, x in one)}. "
               f"In the plane basis pursuit on the grid is certified exact on {sum(x['t_jet'] < 1 for k, x in pl)} of {len(pl)} targets, on "
               f"{sum(x['t_jet'] < 1 for k, x in pl if src[k[0]][k[1]]['out'] == 'other')} of the "
               f"{sum(src[k[0]][k[1]]['out'] == 'other' for k, x in pl)} where the path fits the features with other kernels than E; the full "
               f"linear program in the plane (Part H's t) was not solved here (Part L solves it at three truncations of the span). Wall time per target (one process), median: one dimension {np.median([x['secs'] for k, x in one]):.2g} s, "
               f"the plane {np.median([x['secs'] for k, x in pl]):.2g} s.", ""]
    if K is not None and I is not None and all(x is not None for x in PA.values()):
        src = PA
        fam1 = {k: v for k, v in families(PA['A'], PA['B'], PA['C'], PA['D']).items() if k != 'C'}
        kb = {tuple(r['src']): r['types'] for r in K['rows']}
        hb = {tuple(h['src']): h for h in H} if H is not None else {}
        allk = [k for v in fam1.values() for k in v if k in kb]
        L_ += ["## K. Part I's failures by type (referee round 149)", "",
               "Every stored solution of Part I classified, without a new solve: each nonzero belongs to its nearest feature, within the "
               f"tolerance {K['tol']:g} to its cluster (Part I's rule), beyond it but within {K_RD:g}σ = {K['rd']:g} *displaced*, beyond "
               f"{K['rd']:g} *spurious*. A feature is *lost* when no nonzero within {K['rd']:g} has it as its nearest feature. A target that Part I "
               "does not recover is counted in the first of lost, spurious, displaced, wrong sign (a cluster's sum not of A_k's sign) that "
               "applies; the four cover every failure. The last columns: Part I's column \"a feature lost\" (a cluster empty within the "
               "tolerance) restricted to the failures, the displaced nonzeros' largest distance, and the lost features' |A_k|/max|A|:", "",
               "| family | ℓ | targets | not recovered | a feature lost | else spurious | else displaced | else wrong sign | a cluster empty "
               "within the tolerance | displaced: largest distance | lost: \\|A_k\\|/max\\|A\\| median [max] |", "|---|---|---|---|---|---|---|---|---|---|---|"]
        for fam, keys in list(fam1.items()) + [('all', allk)]:
            ks = [k for k in keys if k in kb]
            for n_l, l in enumerate(I_LAMS):
                ts = [kb[k][n_l] for k in ks]
                bad = [t for t in ts if t['cls'] != 'rec']
                dd = [t['ddisp'] for t in bad if np.isfinite(t['ddisp'])]
                lr = [a for t in bad for a in t['lost_rel']]
                L_.append(f"| {fam} | {l:g} | {len(ts)} | {len(bad)} | " + ' | '.join(str(sum(t['cls'] == c for t in bad))
                                                                            for c in ('lost', 'spurious', 'displaced', 'sign')) +
                          f" | {sum(t['lost_tol'] > 0 for t in bad)} | " + (f"{max(dd):.4f}" if dd else '—') + " | " +
                          (f"{np.median(lr):.3g} [{max(lr):.3g}]" if lr else '—') + " |")
        ts = [t for k in allk for t in kb[k] if t['cls'] != 'rec']
        dd = [t['ddisp'] for t in ts if np.isfinite(t['ddisp'])]
        ds = [t['dspur'] for t in ts if np.isfinite(t['dspur'])]
        sp = [k for k in allk if any(t['cls'] == 'spurious' for t in kb[k])]
        unk = sum(t['cls'] == '?' for t in ts)
        sq, pl = K['seq'], K['plane']
        L_ += ["", f"Over all failures the displaced nonzeros lie between {K['tol']:g} and {max(dd):.4f} of their features and the spurious ones "
               f"at least {min(ds):.4f} from every feature, so the split does not depend on the radius between those two distances; "
               f"unclassified failures: {unk}. The targets with a spurious nonzero at some ℓ: {len(sp)}, of which η_V is degenerate on "
               f"{sum(k in hb and not hb[k]['nondeg'] for k in sp)}.", "",
               f"Part G's audit as entry sequences: the unpolished flow enters the same kernels in the same order as the polished path on "
               f"{sum(x['same'] for x in sq)} of the {len(sq)} one-dimensional targets ({sum(x['src'][0] == 'A' or (x['src'][0] == 'D' and src['D'][x['src'][1]]['tag'][0] == 'b') for x in sq)} "
               f"of A and D(b) from their own rows, the others from Part G's); in the plane, where the unpolished flow's sequence is not stored, "
               f"its outcome class is the polished path's on {sum(x['same'] for x in pl)} of {len(pl)}.", ""]
    if L is not None and all(x is not None for x in PA.values()):
        src = PA
        jb = {tuple(r['src']): r for r in J} if J is not None else {}
        lb = [(tuple(r['src']), r) for r in L]
        n_c = len(L_CUTS)
        L_ += ["## L. Basis pursuit's least certificate in the plane, at three truncations of the span (after round 149)", "",
               "Part H's linear program in the plane, in full: p in the span of the grid's columns (spacing " f"{H_BP[1]:g}" ", "
               f"{len(_grid(2, H_BP[1]))} points over the samples' square), every grid point farther than half a spacing from the features "
               "constrained. On the tensor samples these columns are G1 ⊗ G1 (G1 the one-dimensional Gaussians), so with G1 = U1S1V1ᵀ the span's "
               "orthonormal basis is U1 ⊗ U1 with the singular values s_is_j, truncated at s_is_j > c·s_0² (c = 10⁻¹³ is Part H's truncation; "
               "the dimension below); in the coordinates w = (s_is_j)(u_i ⊗ u_j)ᵀp the grid values are η = (V1 ⊗ V1)w. The program, min t "
               "subject to η(c*_k) = sign A_k and |η| ≤ t on the grid, is solved by a primal–dual interior-point method (Mehrotra's "
               "predictor–corrector) whose normal matrix (V1 ⊗ V1)ᵀD(V1 ⊗ V1) is formed through the tensor factors. Each iterate gives a "
               "certified bracket: t, the sup of |η| over the grid of a w that satisfies the equalities, from above, and t_low, the dual "
               "iterate made exactly feasible (its stationarity residual removed on the grid and its mass renormalised), from below by "
               f"weak duality; the iteration stops when the bracket is below {L_TOL:g} or has not narrowed for {L_STALL} iterations. t_J is "
               "Part J's jet certificate (an upper bound on t, up to the truncation: its p need not lie in the truncated span). |t_dir − t| "
               "is the certificate's p, evaluated directly on the samples, against t, and feq its largest |⟨pφ(c*_k)⟩ − sign A_k|: the "
               "rounding of p's large components in the directions of small s_is_j.", "",
               "| truncation c | dimension of p | targets | t < 1 | t median [min, max] | t ≤ t_J | t_J − t median [min, max] | bracket t − t_low "
               "max | iterations median [max] | \\|t_dir − t\\| max | feq max | seconds median |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for n_l, c in enumerate(L_CUTS):
            xs = [(k, r['cuts'][n_l]) for k, r in lb]
            t = [x['t'] for k, x in xs]
            gap = [jb[k]['t_jet'] - x['t'] for k, x in xs if k in jb]
            L_.append(f"| {c:g} | {xs[0][1]['dim']} | {len(xs)} | {sum(v < 1 for v in t)} | {np.median(t):.5f} [{min(t):.5f}, {max(t):.5f}] | "
                      f"{sum(x['t'] <= jb[k]['t_jet'] for k, x in xs if k in jb)} | {np.median(gap):.4f} [{min(gap):.4f}, {max(gap):.4f}] | "
                      f"{E.fmt(max(x['t'] - x['t_low'] for k, x in xs))} | {np.median([x['it'] for k, x in xs]):g} [{max(x['it'] for k, x in xs)}] | "
                      f"{E.fmt(max(abs(x['t_dir'] - x['t']) for k, x in xs))} | {E.fmt(max(x['feq'] for k, x in xs))} | "
                      f"{np.median([x['secs'] for k, x in xs]):.2g} |")
        mono = sum(all(r['cuts'][i]['t'] <= r['cuts'][i + 1]['t'] + 1e-9 for i in range(n_c - 1)) for k, r in lb)
        oth = [k for k, r in lb if src[k[0]][k[1]]['out'] == 'other']
        g2 = _grid(2, H_BP[1])
        rep = {k: all(np.min(np.linalg.norm(g2 - c[None, :], axis=1)) <= 0.5 * H_BP[1] for c in np.asarray(src[k[0]][k[1]]['cs'], float).reshape(-1, 2))
               for k, r in lb}
        t_r, t_n = [r['cuts'][0]['t'] for k, r in lb if rep[k]], [r['cuts'][0]['t'] for k, r in lb if not rep[k]]
        L_ += ["", f"At Part H's truncation ({L_CUTS[0]:g}) basis pursuit on the grid is exact (t < 1) on {sum(r['cuts'][0]['t'] < 1 for k, r in lb)} "
               f"of {len(lb)} plane targets, on {sum(r['cuts'][0]['t'] < 1 for k, r in lb if k in oth)} of the {len(oth)} where the path fits the "
               "features with other kernels than E; t median by the path's outcome (E, other): "
               f"{np.median([r['cuts'][0]['t'] for k, r in lb if k not in oth]):.5f}, {np.median([r['cuts'][0]['t'] for k, r in lb if k in oth]):.5f}. "
               f"t grows as the span is truncated more (a smaller span, a larger least t) on {mono} of {len(lb)}. t follows the grid's geometry "
               "at the features: where every feature replaced a grid point (one within half a spacing) t has the median [range] "
               f"{np.median(t_r):.5f} [{min(t_r):.5f}, {max(t_r):.5f}] ({len(t_r)} targets), where some feature replaced none "
               f"{np.median(t_n):.5f} [{min(t_n):.5f}, {max(t_n):.5f}] ({len(t_n)}). The points within 10⁻⁷ of t "
               f"(the program's active set), median per truncation: " + ', '.join(
                   f"{c:g}: {np.median([r['cuts'][n_l]['nact'] for k, r in lb]):g}" for n_l, c in enumerate(L_CUTS)) +
               f"; the equalities' largest residual in w: {E.fmt(max(x['re'] for k, r in lb for x in r['cuts']))}.", ""]
        t0 = [r['cuts'][0]['t'] for k, r in lb]
        g0 = [jb[k]['t_jet'] - r['cuts'][0]['t'] for k, r in lb if k in jb]
        L_ += ["At three digits (t at Part H's truncation is meaningful to about 10⁻³, see \\|t_dir − t\\|), the quantities the paper quotes: "
               f"t at {L_CUTS[0]:g} min, median, max {min(t0):.3f}, {np.median(t0):.3f}, {max(t0):.3f}; t median per truncation "
               + ', '.join(f"{c:g}: {np.median([r['cuts'][n_l]['t'] for k, r in lb]):.3f}" for n_l, c in enumerate(L_CUTS)) +
               f"; t_J − t at {L_CUTS[0]:g} min, max {min(g0):.3f}, {max(g0):.3f}; t median where every feature replaced a grid point, where "
               f"some feature replaced none: {np.median(t_r):.3f}, {np.median(t_n):.3f}; t median by the path's outcome (E, other): "
               f"{np.median([r['cuts'][0]['t'] for k, r in lb if k not in oth]):.3f}, {np.median([r['cuts'][0]['t'] for k, r in lb if k in oth]):.3f}.", ""]
    if M is not None and all(x is not None for x in PA.values()):
        src = PA
        mb = [(tuple(r['src']), r) for r in M]
        hb = {tuple(h['src']): h for h in H} if H is not None else {}
        K2 = TN.full()[0]
        tol, rd = I_TOL * H_BP[1], K_RD * K2.sig
        V = {k: [blasso_verdict(s_, src[k[0]][k[1]]['cs'], src[k[0]][k[1]]['As'], tol) for s_ in r['sols']] for k, r in mb}
        Ty = {k: [k_classify(s_, src[k[0]][k[1]]['cs'], src[k[0]][k[1]]['As'], tol, rd) for s_ in r['sols']] for k, r in mb}
        kk = [s_['kkt'] for k, r in mb for s_ in r['sols']]
        L_ += ["## M. The BLASSO in the plane (after round 149)", "",
               "Part I's BLASSO on the plane targets: the same problem on Part H's basis-pursuit grid in the plane (spacing " f"{H_BP[1]:g}" ", "
               "its points within half a spacing of a feature replaced by the feature), λ = ℓλ_max for ℓ = " + ', '.join(f"{l:g}" for l in I_LAMS) +
               ", feature-sign search warm started from the larger λ, the Gram matrix applied through the samples' tensor factors (never "
               f"formed); the largest KKT violation over all {len(kk)} solves {E.fmt(max(kk))} relative to λ. Part I's recovery rule with the "
               f"plane's tolerance: every nonzero within {I_TOL} spacings ({tol:g}) of a feature, each feature's cluster nonempty with its "
               f"sum of A_k's sign. Failures by Part K's types with the radius {K_RD:g}σ = {rd:g} (σ = {K2.sig:g}): a feature lost (no "
               f"nonzero within {rd:g} has it as its nearest feature), else a spurious nonzero (beyond {rd:g} of every feature), else a "
               "displaced one (beyond the tolerance), else a cluster of the wrong sign:", "",
               "| ℓ | targets | recovered | support exactly the features | a cluster empty within the tolerance | failures: lost, spurious, "
               "displaced, wrong sign | amplitude error median [max] | error/ℓ median | recovered: path E, path other | support beyond the "
               "features median [max] |", "|---|---|---|---|---|---|---|---|---|---|"]
        for n_l, l in enumerate(I_LAMS):
            vs = [(k, V[k][n_l]) for k, r in mb]
            er = [v['err'] for k, v in vs if v['rec']]
            bad = [Ty[k][n_l]['cls'] for k, r in mb if Ty[k][n_l]['cls'] != 'rec']
            ex = [v['n'] - len(src[k[0]][k[1]]['As']) for k, v in vs]
            L_.append(f"| {l:g} | {len(vs)} | {sum(v['rec'] for k, v in vs)} | {sum(v['exact'] for k, v in vs)} | {sum(v['lost'] > 0 for k, v in vs)} | "
                      + ', '.join(str(bad.count(c)) for c in ('lost', 'spurious', 'displaced', 'sign')) + " | "
                      + (f"{E.fmt(float(np.median(er)))} [{E.fmt(max(er))}]" if er else '—') + " | " + (f"{np.median(er) / l:.3g}" if er else '—') +
                      f" | {sum(v['rec'] for k, v in vs if src[k[0]][k[1]]['out'] == 'nearest')} of {sum(src[k[0]][k[1]]['out'] == 'nearest' for k, v in vs)}, "
                      f"{sum(v['rec'] for k, v in vs if src[k[0]][k[1]]['out'] == 'other')} of {sum(src[k[0]][k[1]]['out'] == 'other' for k, v in vs)} | "
                      f"{np.median(ex):g} [{max(ex)}] |")
        lo = len(I_LAMS) - 1
        nd = [k for k, r in mb if k in hb and hb[k]['nondeg']]
        dd = [t['ddisp'] for k, r in mb for t in Ty[k] if t['cls'] != 'rec' and np.isfinite(t['ddisp'])]
        ds = [t['dspur'] for k, r in mb for t in Ty[k] if t['cls'] != 'rec' and np.isfinite(t['dspur'])]
        L_ += ["", f"At ℓ = {I_LAMS[lo]:g}: recovered on {sum(V[k][lo]['rec'] for k in nd)} of the {len(nd)} targets with η_V nondegenerate "
               f"(Part H) and on {sum(V[k][lo]['rec'] for k, r in mb if k not in nd)} of the {len(mb) - len(nd)} without. Over the failures the "
               "displaced nonzeros' largest distance from their feature: " + (f"{max(dd):.4f}" if dd else '—') + "; the spurious ones' "
               "smallest distance from every feature: " + (f"{min(ds):.4f}" if ds else '—') + "; unclassified failures: "
               f"{sum(t['cls'] == '?' for k, r in mb for t in Ty[k])}. Feature-sign iterations per solve, median [max]: "
               f"{np.median([s_['it'] for k, r in mb for s_ in r['sols']]):g} [{max(s_['it'] for k, r in mb for s_ in r['sols'])}]; wall time "
               f"per target (four λ), median {np.median([r['secs'] for k, r in mb]):.2g} s.", ""]
    L_.append("Numbers are produced by `experiments/exp_transit_certificate.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'transit_certificate_raw.npz')
    md = os.path.join(OUT, 'transit_certificate.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in 'ABCDEFGHIJKLM']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFGHIJKLM')
    res = dict(old)
    t0 = time.time()
    save = lambda: E.save(raw, **res)
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', lambda: part_D(res['A'])), ('E', part_E), ('F', part_F),
                  ('G', lambda: part_G(res['B'], res['C'], res['D'])), ('H', lambda: part_H(res['A'], res['B'], res['C'], res['D'])),
                  ('I', lambda: run_chunks(i_job, i_jobs(res['A'], res['B'], res['D']), res.setdefault('I', []), save)),
                  ('J', lambda: run_chunks(j_job, j_jobs(res['A'], res['B'], res['C'], res['D']), res.setdefault('J', []), save)),
                  ('K', lambda: part_K(res['A'], res['B'], res['D'], res['G'], res['I'])),
                  ('L', lambda: run_chunks(l_job, l_jobs(res['C']), res.setdefault('L', []), save, L_CHUNK)),
                  ('M', lambda: run_chunks(m_job, l_jobs(res['C']), res.setdefault('M', []), save))):
        if k in parts:                                                 # I, J, L and M resume from the rows already in the raw
            res[k] = fn(); save(); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[copy.deepcopy(res.get(k)) for k in 'ABCDEFGHIJKLM']))   # the report fills margins in place


if __name__ == '__main__':
    main()
