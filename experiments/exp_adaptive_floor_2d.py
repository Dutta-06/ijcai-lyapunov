"""
Proposition D‴.24 in two dimensions: the floor is the leak's bias, and an adaptive prior removes it without a mass constraint
=============================================================================================================================

The planar instance of `exp_ot_2d.py` — ẋ = f(x) + u on ℝ², f four Gaussian bumps with vector amplitudes, the Lissajous reference
(2 sin t, 1.5 sin 2t), the saturated network with p = 5 parameters per kernel, G = diag(1, 1, 0.25, 0.25, 0.25), σ₀ = 0.05, the matching
prior = the exact representation with masses 1/4 — and the averaged slow loop of `src/core/planar_loop.py` (ρ: 400 times of one period).

A. The law, no dynamics: the excitation spectrum h_j at the exact representation; the function-level sensitivity F of the fixed prior's
   weighted stationary point (central differences in w along e_j − e_4) against the bias form Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰; the law per unit of
   1/m − 1/M from F and as the sum of the mode terms, against the bound (σ₀/4)tr(∂_wx⁰C̃∂_wx⁰ᵀ); the representable weights W (an exact
   representation at w needs |A_{k,c}|/w_k < a_max for both components).
B. The nearest-atom leak from the warm clouds of `exp_ot_2d.py` (width 0.02), m = 16, 32, 64 against M = 128 (the first m rows), six seeds,
   fixed (η = 0) against adaptive priors (η = σ₀ and ∞; η = ∞: the atoms reset to the cell means at every refresh of the cells and frozen in
   between, a lagged prior), S_B slow units, midpoint rule dt = 0.25, the cells every 4 steps (one slow unit): the sampled
   weights and their representability margin, ‖f̂_m − f̂_M‖² and the fit error at three times, the class means' distance from the exact
   representation of the sampled weights at the end, W²(μ_m, μ_M) at the end against the W² between the two systems' exact
   representations, the counts.
C. The real closed loop (k_e = 25, σ = σ₀/k_e, the planar RK4 plant of `exp_ot_2d.py`, dt = 2·10⁻³, the cells every 100 steps; η = ∞: the
   atoms are the cell means at every refresh) for one seed and size whose sample is representable: fixed and adaptive priors from the warm cloud,
   and the adaptive prior from the exact representation of the sampled weights with the atoms there; the rms tracking error over windows of
   five periods against (1/k_e)‖f − f̂‖_ρ of the averaged loop from the same start at s = t/k_e (`--only C` recomputes it and merges it).
The fixed prior at the sampled weights (`--only P` computes it from the stored records and merges it): per seed and size its replicated
   stationary point Θ*(w), Proposition D‴.5's splitting margins and the collective block's largest eigenvalue, and the exact bias at the
   sampled weights against the observed floor.

Output: results/chaos/adaptive_floor_2d.md and adaptive_floor_2d_raw.npz (AF2_WORKERS, AF2_SB, AF2_TC, AF2_OUT; `--from-raw` regenerates the md).
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
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import planar_loop as PL
from src.core.ot_leak import transport_plan
import exp_ot_2d as P2

OUT = os.environ.get('AF2_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
G = P2.G
SIGMA0 = P2.SIGMA0
Q = P2.Q
K0 = 4
A_B = np.array([b[0] for b in P2.BUMPS]); C_B = np.array([b[1] for b in P2.BUMPS]); S_B_ = np.array([b[2] for b in P2.BUMPS])
MS, M_REF = (16, 32, 64), 128
SEEDS = (0, 1, 2, 3, 4, 5)
ETAS = (0.0, SIGMA0, np.inf)
S_B = float(os.environ.get('AF2_SB', 6000.0))
DT, CELL_EVERY = 0.25, 4
REC_TIMES = tuple(t for t in (300.0, 1000.0, 3000.0) if t < S_B) + (S_B,)
N_WORKERS = int(os.environ.get('AF2_WORKERS', 3))
DELTA_W = 1e-4
T_C = float(os.environ.get('AF2_TC', 2500.0))
KE, DT_REAL, CELL_REAL = P2.KE, P2.DT, 100


def setup():
    net = P2.make_net()
    exact = P2.exact_representation(net)
    loop = PL.AveragedLoopN(net, P2.f_vec, P2.x_ref, P2.PER, G, SIGMA0, n_quad=400)
    return net, exact, loop


def exact_at(net, w):
    """The exact representation at the weights w (the amplitudes A_k/w_k), or None if some |A_{k,c}|/w_k ≥ a_max."""
    w = np.asarray(w, float)
    if margin(net, w) <= 0:
        return None
    return net.unconstrained(A_B / w[:, None], C_B, S_B_)


def margin(net, w):
    """min_k(w_k − max_c|A_{k,c}|/a_max): positive iff the weights are representable (in W)."""
    return float(np.min(np.asarray(w, float) - np.abs(A_B).max(axis=1) / net.a_max))


def directions():
    """e_j − e_4, j = 1, 2, 3."""
    D = np.zeros((K0 - 1, K0))
    for j in range(K0 - 1):
        D[j, j], D[j, -1] = 1.0, -1.0
    return D


# ---------------------------------------------------------------------------- A
def part_A():
    net, exact, loop = setup()
    H, h, V = PL.gauss_newton(loop, Q, exact)
    J = PL.whitened_jacobian(loop, Q, exact)
    D = directions()
    F = np.zeros((loop.N * loop.n, K0 - 1)); dx0 = np.zeros((K0 * loop.p, K0 - 1)); res = []
    for j in range(K0 - 1):
        fh, x0 = [], []
        for sgn in (+1, -1):
            w = Q + sgn * DELTA_W * D[j]
            at, r, own = PL.weighted_stationary(loop, w, exact)
            res.append((r, own))
            fh.append(loop.f_hat(at, w).ravel())
            x0.append(PL.whitened(exact - exact_at(net, w), w, G))
        F[:, j] = (fh[0] - fh[1]) / (2 * DELTA_W)
        dx0[:, j] = (x0[0] - x0[1]) / (2 * DELTA_W)
    Fb = J @ (SIGMA0 * np.linalg.solve(H + SIGMA0 * np.eye(H.shape[0]), dx0))
    Ct = np.diag(Q[:-1]) - np.outer(Q[:-1], Q[:-1])
    law = float(np.trace(F @ Ct @ F.T)) / loop.N
    z = V.T @ dx0
    terms = h * SIGMA0 ** 2 / (h + SIGMA0) ** 2 * np.einsum('jk,kl,jl->j', z, Ct, z)
    bound = 0.25 * SIGMA0 * float(np.trace(dx0 @ Ct @ dx0.T))
    at_q, r_q, own_q = PL.weighted_stationary(loop, Q, exact)
    return dict(h=h, rel=float(np.linalg.norm(F - Fb) / np.linalg.norm(F)), law=law, terms=terms, bound=bound,
                newton=max(r for r, _ in res), own=all(o for _, o in res), fit_q=loop.fit_error(at_q, Q), res_q=r_q,
                lo=np.abs(A_B).max(axis=1) / net.a_max)


# ---------------------------------------------------------------------------- B
def run_job(args):
    seed, eta = args
    t0 = time.time()
    net, exact, loop = setup()
    clouds = P2.mixture_start(exact, M_REF, np.random.default_rng(seed), P2.JITTER)
    out = dict(seed=seed, eta=eta, rec={})
    for m in MS + (M_REF,):
        Th0 = clouds[:m].copy()
        lab0 = PL.cells(Th0, exact, G)
        rec = dict(s=[], fh=[], fit=[], counts=[])

        def cb(s, Th, Z, rec=rec):
            if any(abs(s - t) < 1e-9 for t in (0.0,) + REC_TIMES):
                rec['s'].append(s); rec['fh'].append(loop.f_hat(Th).astype(np.float32)); rec['fit'].append(loop.fit_error(Th))
                rec['counts'].append(np.bincount(PL.cells(Th, Z, G), minlength=K0))
        Th, Z = PL.run(loop, Th0, exact.copy(), S_B, DT, eta=eta, cell_every=CELL_EVERY, rec_every=int(round(100.0 / DT)), callback=cb)
        out['rec'][m] = dict(s=np.array(rec['s']), fh=np.array(rec['fh']), fit=np.array(rec['fit']), counts=np.array(rec['counts']),
                             counts0=np.bincount(lab0, minlength=K0), final=Th, Z=Z)
    print(f"[{time.time() - t0:.0f}s] B seed {seed} eta {eta:g}", flush=True)
    return out


def class_means(Th, lab):
    return np.array([Th[lab == k].mean(axis=0) for k in range(K0)])


def w2_atoms(A1, w1, A2, w2):
    C = PL.cost_matrix(A1, A2, G)
    return float(np.sum(transport_plan(C, np.asarray(w1, float), np.asarray(w2, float)) * C))


def w2_empirical(T1, T2):
    m1, m2 = T1.shape[0], T2.shape[0]
    C = PL.cost_matrix(T1, T2, G)
    return float(np.sum(transport_plan(C, np.full(m1, 1.0 / m1), np.full(m2, 1.0 / m2)) * C))


def post_B(B):
    """Per system: the sampled weights, the margin, the distances at the record times, the end W² against the exact representations'."""
    net, exact, loop = setup()
    rows = []
    for r in B:
        R = r['rec']
        wM = R[M_REF]['counts0'] / M_REF; eM = exact_at(net, wM)
        for m in MS + (M_REF,):
            o = R[m]
            w = o['counts0'] / m; e = exact_at(net, w)
            d = [float(np.mean(np.sum((o['fh'][k].astype(float) - R[M_REF]['fh'][k].astype(float)) ** 2, axis=1))) for k in range(len(o['s']))]
            lab = PL.cells(o['final'], o['Z'], G)
            dev = float(np.linalg.norm(PL.whitened(class_means(o['final'], lab) - e, w, G))) if e is not None and np.all(np.bincount(lab, minlength=K0) > 0) else float('nan')
            w2 = w2_empirical(o['final'], R[M_REF]['final']) if m != M_REF else 0.0
            w2e = w2_atoms(e, w, eM, wM) if (e is not None and eM is not None and m != M_REF) else float('nan')
            rows.append(dict(seed=r['seed'], eta=r['eta'], m=m, counts=tuple(int(x) for x in o['counts0']), margin=margin(net, w), s=o['s'], dist=d,
                             fit=o['fit'], conserved=bool(np.all(o['counts'] == o['counts0'][None, :])), dev=dev, w2=w2, w2e=w2e))
    return rows


# ---------------------------------------------------------------------------- the fixed prior at the sampled weights (referee round 18's lesson)
def fixed_point_2d(loop, counts, prior, h=1e-6):
    """The fixed prior's replicated stationary point at the counts: Θ*(w) (hybrid Powell from the prior), its residual, whether every atom is in
    its own cell, Proposition D‴.5's splitting margins σ₀ − λ_max(G^{1/2}R_kG^{1/2}) (R_k = ⟨rᵀ∇²_ϑΦ̃(ϑ*_k)⟩_ρ, the learning force's Jacobian at
    the frozen residual; +∞ for a class with fewer than two particles), the largest eigenvalue of the weighted K-atom drift's Jacobian (the
    collective block; negative iff it is stable), f̂ on the grid, and stability (D‴.5-admissibility)."""
    counts = np.asarray(counts, int); w = counts / counts.sum()
    K, p = prior.shape
    atoms, res, own = PL.weighted_stationary(loop, w, prior)
    r = loop.residual(atoms, w)
    sqG = np.sqrt(G)
    sep = np.full(K, np.inf)
    for k in range(K):
        if counts[k] < 2:
            continue
        Rk = np.zeros((p, p))
        for c in range(p):
            e = np.zeros(p); e[c] = h
            Rk[:, c] = (loop.learning((atoms[k] + e)[None, :], r=r)[0] - loop.learning((atoms[k] - e)[None, :], r=r)[0]) / (2 * h)
        Rk = 0.5 * (Rk + Rk.T)
        sep[k] = SIGMA0 - float(np.linalg.eigvalsh(sqG[:, None] * Rk * sqG[None, :])[-1])
    J = np.zeros((K * p, K * p))
    for j in range(K * p):
        e = np.zeros(K * p); e[j] = h
        J[:, j] = ((loop.drift((atoms.ravel() + e).reshape(K, p), prior, w) - loop.drift((atoms.ravel() - e).reshape(K, p), prior, w)) / (2 * h)).ravel()
    coll = float(np.max(np.linalg.eigvals(J).real))
    return dict(atoms=atoms, res=res, own=own, sep=sep, coll=coll, fhat=loop.f_hat(atoms, w),
                stable=bool(np.all(sep > 0) and coll < 0 and own and res < 1e-9))


def checks_2d(rows):
    """Per seed and size (references included) the fixed prior's Θ*(w) at the sampled counts; per m-system the exact bias at the sampled
    weights, ‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖², against the observed floor of the fixed prior at the end."""
    net, exact, loop = setup()
    fp = {}
    for r in rows:
        if r['eta'] == 0.0:
            fp[(r['seed'], r['m'])] = fixed_point_2d(loop, r['counts'], exact)
    out = {}
    for r in rows:
        if r['eta'] == 0.0 and r['m'] != M_REF:
            a, b = fp[(r['seed'], r['m'])], fp[(r['seed'], M_REF)]
            ex = float(np.mean(np.sum((a['fhat'] - b['fhat']) ** 2, axis=1)))
            out[(r['seed'], r['m'])] = dict(h2_exact=ex, h2_obs=float(r['dist'][-1]))
    strip = lambda d: {k: v for k, v in d.items() if k != 'fhat'}
    return dict(fp={k: strip(v) for k, v in fp.items()}, floor=out)


def report_checks_2d(rows, P):
    fp, fl = P['fp'], P['floor']
    lines = ["", "### The fixed prior at the sampled weights: stability and the floor as the exact bias", "",
             "Per seed and size: the fixed prior's replicated stationary point Θ*(w) at the sampled counts, Proposition D‴.5's splitting margin "
             "min_k(σ₀ − λ_max(G^{1/2}R_kG^{1/2})) over the classes with at least two particles (its class), the largest eigenvalue of the weighted K-atom "
             "drift's Jacobian (negative: the collective block is stable); Θ*(w) is a strict local minimiser (D‴.5-admissible) iff the margin is positive and "
             "the eigenvalue negative. For the m-systems: the fixed prior's observed ‖f̂_m − f̂_M‖² at the end against the exact bias at the sampled weights, "
             "‖f̂_{Θ*(w^m)} − f̂_{Θ*(w^M)}‖² (ratio 1 for a system that stays at its replicated point).", "",
             "| seed | m | counts | splitting margin (class) | collective: largest eigenvalue | observed / exact bias |", "|---|---|---|---|---|---|"]
    stab, sad = [], []
    for (seed, m) in sorted(fp):
        f = fp[(seed, m)]; k = int(np.argmin(f['sep']))
        cnt = next(r['counts'] for r in rows if r['seed'] == seed and r['m'] == m and r['eta'] == 0.0)
        cell = "—"
        if (seed, m) in fl:
            q_ = fl[(seed, m)]['h2_obs'] / fl[(seed, m)]['h2_exact']
            cell = f"{fl[(seed, m)]['h2_obs']:.3e} / {fl[(seed, m)]['h2_exact']:.3e} ({q_:.4f})"
            (stab if (f['stable'] and fp[(seed, M_REF)]['stable']) else sad).append(q_)
        lines.append(f"| {seed} | {m} | {cnt} | {f['sep'][k]:+.2e} (class {k + 1}){'' if f['stable'] else ' saddle'} | {f['coll']:+.2e} | {cell} |")
    rng_ = lambda x: f"{min(x):.4f}–{max(x):.4f}" if x else "—"
    n_sad = sum(1 for (s_, m_), f in fp.items() if m_ != M_REF and not f['stable'])
    lines += ["", f"{n_sad} of {sum(1 for (s_, m_) in fp if m_ != M_REF)} m-systems (and {sum(1 for (s_, m_), f in fp.items() if m_ == M_REF and not f['stable'])} "
              f"of the {sum(1 for (s_, m_) in fp if m_ == M_REF)} references) sit at weights where the fixed prior's replicated point is a saddle. At the stable pairs the "
              f"observed floor is {rng_(stab)} of the exact bias at the sampled weights; at the others {rng_(sad)}."]
    return lines


# ---------------------------------------------------------------------------- C
def real_loop(net, Th0, Z0, eta_inf, T):
    """The planar real loop with the nearest-atom leak toward Z (fixed, or the cell means at every refresh): e² every 0.1 time units."""
    x = P2.x_ref(0.0).copy(); Th = np.array(Th0, float); Z = np.array(Z0, float)
    sig = SIGMA0 / KE
    n = int(round(T / DT_REAL)); t = 0.0
    ts, e2 = [], []
    tgt = None
    for i in range(n):
        if i % CELL_REAL == 0:
            lab = PL.cells(Th, Z, G)
            if eta_inf:
                Z = PL.cell_means(Th, lab, Z)
            tgt = Z[lab]
        e = x - P2.x_ref(t)
        if i % 50 == 0:
            ts.append(t); e2.append(float(e @ e))
        xdot = lambda y, tau: P2.f_vec(y)[0] + P2.x_ref_dot(tau) - KE * (y - P2.x_ref(tau)) - net.f_hat(y, Th)
        k1 = xdot(x, t); k2 = xdot(x + 0.5 * DT_REAL * k1, t + 0.5 * DT_REAL); k3 = xdot(x + 0.5 * DT_REAL * k2, t + 0.5 * DT_REAL)
        k4 = xdot(x + DT_REAL * k3, t + DT_REAL)
        Thn = Th + DT_REAL * (G[None, :] * net.grad_g(x, e, Th) - sig * (Th - tgt))
        x = x + DT_REAL / 6 * (k1 + 2 * k2 + 2 * k3 + k4); Th = Thn
        t = (i + 1) * DT_REAL
    return np.array(ts), np.array(e2), Th, Z


def part_C_job(args):
    seed, m, case = args
    t0 = time.time()
    net, exact, loop = setup()
    clouds = P2.mixture_start(exact, M_REF, np.random.default_rng(seed), P2.JITTER)
    Th0 = clouds[:m].copy()
    lab = PL.cells(Th0, exact, G); counts = np.bincount(lab, minlength=K0); w = counts / m
    e = exact_at(net, w)
    if case == 'exact':
        Th0 = e[np.repeat(np.arange(K0), counts)]; Z0 = e.copy(); eta_inf = True
    else:
        Z0 = exact.copy(); eta_inf = case == 'adaptive'
    ts, e2, Th, Z = real_loop(net, Th0, Z0, eta_inf, T_C)
    rec = dict(s=[], fit=[])
    PL.run(loop, Th0, Z0, T_C / KE, 0.05, eta=(np.inf if eta_inf else 0.0), cell_every=4, rec_every=20,
           callback=lambda s_, T_, Z_: (rec['s'].append(s_), rec['fit'].append(loop.fit_error(T_))))
    print(f"[{time.time() - t0:.0f}s] C seed {seed} m {m} {case}", flush=True)
    return dict(seed=seed, m=m, case=case, counts=counts, margin=margin(net, w), ts=ts, e2=e2, s_avg=np.array(rec['s']), fit_avg=np.array(rec['fit']),
                normF=float(np.sqrt(loop.F2)), final=Th.astype(np.float32))


def pick_C(rows):
    """The first (seed, m) with m ≥ 32 whose sample is representable."""
    for r in rows:
        if r['eta'] == 0.0 and r['m'] in (32, 64) and r['margin'] > 0:
            return r['seed'], r['m']
    return None


# ---------------------------------------------------------------------------- report
def report(A, rows, C=None, t0=None, P=None):
    eta_s = lambda e: '∞' if np.isinf(e) else f'{e:g}'
    lines = ["# Proposition D‴.24 in two dimensions: the floor is the leak's bias, and an adaptive prior removes it (`exp_adaptive_floor_2d.py`)", "",
             f"The planar instance of `exp_ot_2d.py` (four bumps with vector amplitudes, the Lissajous reference, p = 5, G = diag{tuple(float(g) for g in G)}, "
             f"σ₀ = {SIGMA0:g}, the matching prior = the exact representation with masses 1/4) and the averaged slow loop of `src/core/planar_loop.py` "
             f"(ρ: 400 times of one period).", "",
             "## A. The floor as the weight-derivative of the leak's bias (matching prior, no dynamics)", "",
             "The whitened weighted Gauss–Newton spectrum at the exact representation (h_j/σ₀): " + ", ".join(f"{x / SIGMA0:.3f}" for x in A['h']) + ".", "",
             f"The fixed prior's weighted stationary points at q ± {DELTA_W:g}(e_j − e_4) by Newton (largest residual {A['newton']:.1e}; every atom in its own cell: "
             f"{A['own']}); at q itself the fit error is {A['fit_q']:.1e} (residual {A['res_q']:.1e}). F against the bias form Ĵσ₀(H_L + σ₀I)⁻¹∂_wx⁰: "
             f"‖F − F_bias‖/‖F‖ = {A['rel']:.1e}. The law per unit of 1/m − 1/M: {A['law']:.4e} from F, {np.sum(A['terms']):.4e} as the sum of the mode terms, "
             f"against the bound {A['bound']:.4e} (law/bound {A['law'] / A['bound']:.3f}).", "",
             "| mode j | h_j/σ₀ | share of the law |", "|---|---|---|"]
    tot = float(np.sum(A['terms']))
    for j in range(len(A['h'])):
        lines.append(f"| {j + 1} | {A['h'][j] / SIGMA0:.3f} | {A['terms'][j] / tot:.3f} |")
    lines += ["", f"Representable weights (the set W): an exact representation at w needs |A_k,c|/w_k < a_max, i.e. w > ({', '.join(f'{x:.3f}' for x in A['lo'])}).", "",
              f"## B. The nearest-atom leak from the warm clouds: fixed (η = 0) against adaptive priors ({S_B:g} slow units)", "",
              f"Midpoint rule dt = {DT:g}, the cells every {CELL_EVERY} steps; η = ∞: the atoms reset to the cell means at every refresh of the cells (every "
              f"{CELL_EVERY * DT:g} slow unit{'' if CELL_EVERY * DT == 1 else 's'}) and frozen in between — a lagged form of the η = ∞ prior.", "",
              "‖f̂_m − f̂_M‖² (the function-level floor) and the fit error at s = " + ", ".join(f"{x:g}" for x in REC_TIMES) + "; the sampled counts and their margin "
              "min_k(w_k − max_c|A_k,c|/a_max) (representable iff positive); at the end the class means' whitened weighted distance from the exact representation of the "
              "sampled weights, and W²(μ_m, μ_M) against the W² between the two systems' exact representations (—: no exact representation); the counts conserved.", "",
              "| seed | m | counts (margin) | η | ‖f̂_m − f̂_M‖² at those s | fit error at those s | distance from ϑ^e(w) | W² / exact-representation W² | conserved |",
              "|---|---|---|---|---|---|---|---|---|"]
    idx = lambda r: [int(np.argmin(np.abs(r['s'] - t))) for t in REC_TIMES]
    for r in sorted(rows, key=lambda r: (r['seed'], r['m'], r['eta'])):
        if r['m'] == M_REF:
            continue
        ii = idx(r)
        lines.append(f"| {r['seed']} | {r['m']} | {r['counts']} ({r['margin']:+.3f}) | {eta_s(r['eta'])} | " + ", ".join(f"{r['dist'][k]:.1e}" for k in ii) + " | "
                     + ", ".join(f"{r['fit'][k]:.1e}" for k in ii) + " | " + ("—" if not np.isfinite(r['dev']) else f"{r['dev']:.1e}") + " | "
                     + f"{r['w2']:.3f} / " + ("—" if not np.isfinite(r['w2e']) else f"{r['w2e']:.3f}") + f" | {'yes' if r['conserved'] else 'no'} |")
    refs = [r for r in rows if r['m'] == M_REF and r['eta'] == 0.0]
    lines += ["", "The reference systems (M = 128): " + "; ".join(f"seed {r['seed']}: {r['counts']} ({r['margin']:+.3f})" for r in refs) + "."]
    # the reading, computed
    def last(r):
        return r['dist'][idx(r)[-1]]
    adm = lambda r: r['margin'] > 0 and next(x for x in rows if x['seed'] == r['seed'] and x['m'] == M_REF and x['eta'] == r['eta'])['margin'] > 0
    fixed = {(r['seed'], r['m']): last(r) for r in rows if r['eta'] == 0.0}
    sub = [r for r in rows if r['m'] != M_REF and r['eta'] > 0]
    fx_a = [fixed[(r['seed'], r['m'])] for r in rows if r['m'] != M_REF and r['eta'] == 0.0 and adm(r)]
    fx_n = [fixed[(r['seed'], r['m'])] for r in rows if r['m'] != M_REF and r['eta'] == 0.0 and not adm(r)]
    ratio = lambda r: last(r) / fixed[(r['seed'], r['m'])]
    falling = lambda r: r['dist'][idx(r)[-1]] < r['dist'][idx(r)[-2]]
    parts = []
    for e in ETAS:
        if e == 0:
            continue
        A_ = [r for r in sub if r['eta'] == e and adm(r)]; N_ = [r for r in sub if r['eta'] == e and not adm(r)]
        if A_:
            parts.append(f"η = {eta_s(e)}, representable ({len(A_)} systems): {min(last(r) for r in A_):.1e}–{max(last(r) for r in A_):.1e}, "
                         f"{min(ratio(r) for r in A_):.3f}–{max(ratio(r) for r in A_):.3f} of the fixed prior's, falling between the last two records in "
                         f"{sum(falling(r) for r in A_)} of {len(A_)}")
        if N_:
            parts.append(f"η = {eta_s(e)}, not representable ({len(N_)} systems): {min(last(r) for r in N_):.1e}–{max(last(r) for r in N_):.1e}, "
                         f"{min(ratio(r) for r in N_):.3f}–{max(ratio(r) for r in N_):.3f} of the fixed prior's, falling between the last two records in "
                         f"{sum(falling(r) for r in N_)} of {len(N_)}")
    lines += ["", f"Reading (at s = {S_B:g}; a system counts as representable when its sample and its reference's are). The fixed prior's ‖f̂_m − f̂_M‖²: "
              f"{min(fx_a):.1e}–{max(fx_a):.1e} at the representable systems, {min(fx_n):.1e}–{max(fx_n):.1e} at the others. The adaptive prior's: " + "; ".join(parts) + "."]
    if P:
        lines += report_checks_2d(rows, P)
    if C:
        lines += ["", f"## C. The real closed loop (k_e = {KE:g}, dt = {DT_REAL:g}, {T_C:g} time units), nearest-atom leak, seed {C[0]['seed']}, m = {C[0]['m']}, "
                  f"counts {tuple(int(x) for x in C[0]['counts'])} (representability margin {C[0]['margin']:+.3f})", "",
                  "The rms tracking error over the five periods ending at the times below, against (1/k_e)‖f − f̂‖_ρ of the averaged loop from the same start at s = t/k_e.", "",
                  "| start, prior | t | e_rms | (1/k_e)‖f − f̂‖_ρ averaged | ratio |", "|---|---|---|---|---|"]
        win = 5 * P2.PER
        for c in C:
            name = {'fixed': 'warm cloud, fixed prior', 'adaptive': 'warm cloud, adaptive prior (η = ∞)', 'exact': 'ϑ^e(w), adaptive prior at ϑ^e(w)'}[c['case']]
            for tt in (win, 500.0, 1000.0, T_C):
                msk = (c['ts'] > tt - win) & (c['ts'] <= tt + 1e-9)
                er = float(np.sqrt(np.mean(c['e2'][msk])))
                pr = float(np.interp(tt / KE, c['s_avg'], c['fit_avg'])) * c['normF'] / KE
                lines.append(f"| {name} | {tt:.0f} | {er:.3e} | {pr:.3e} | " + (f"{er / pr:.3f}" if pr > 1e-9 else "—") + " |")
    lines += ["", "Numbers are produced by `experiments/exp_adaptive_floor_2d.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'adaptive_floor_2d.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def save_raw(raw, **arrays):
    """np.savez to a temporary file, then an atomic rename."""
    tmp = raw[:-4] + '.tmp.npz'
    np.savez(tmp, **arrays)
    os.replace(tmp, raw)


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'adaptive_floor_2d_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        report(save['A'].item(), list(save['rows']), list(save['C']) if 'C' in save else None, t0, save['P'].item() if 'P' in save else None)
        return
    if '--only' in sys.argv and sys.argv[sys.argv.index('--only') + 1] == 'P':          # the fixed prior at the sampled weights
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        P = checks_2d(list(save['rows']))
        save['P'] = np.array(P, dtype=object)
        save_raw(raw, **save)
        report(save['A'].item(), list(save['rows']), list(save['C']) if 'C' in save else None, t0, P)
        return
    if '--only' in sys.argv and sys.argv[sys.argv.index('--only') + 1] == 'C':
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        rows = list(save['rows'])
        sm = pick_C(rows)
        with mp.Pool(N_WORKERS) as pool:
            C = pool.map(part_C_job, [(sm[0], sm[1], case) for case in ('fixed', 'adaptive', 'exact')], chunksize=1)
        save['C'] = np.array(C, dtype=object)
        save_raw(raw, **save)
        report(save['A'].item(), rows, C, t0, save['P'].item() if 'P' in save else None)
        return
    A = part_A()
    with mp.Pool(N_WORKERS) as pool:
        B = pool.map(run_job, [(seed, eta) for seed in SEEDS for eta in ETAS], chunksize=1)
    rows = post_B(B)
    P = checks_2d(rows)
    save_raw(raw, A=np.array(A, dtype=object), rows=np.array(rows, dtype=object), P=np.array(P, dtype=object))
    sm = pick_C(rows)
    C = None
    if sm is not None and '--no-C' not in sys.argv:
        with mp.Pool(N_WORKERS) as pool:
            C = pool.map(part_C_job, [(sm[0], sm[1], case) for case in ('fixed', 'adaptive', 'exact')], chunksize=1)
        save_raw(raw, A=np.array(A, dtype=object), rows=np.array(rows, dtype=object), P=np.array(P, dtype=object), C=np.array(C, dtype=object))
    report(A, rows, C, t0, P)


if __name__ == '__main__':
    main()
