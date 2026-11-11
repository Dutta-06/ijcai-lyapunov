"""
Proposition D‴.54: the transit in ℝⁿ — one feature on a lattice of kernels (the k-fold near tie and its shares, the formation map at a
vertex) and several features in the plane
====================================================================================================================================

The free law with moving centres (Proposition D‴.48(iii); m = 12, G = 1, γ = 0.1, a_max = 1, c = 0.025), Gaussian kernels of width SIG in
the plane (`rich_transit_nd`, the samples on a tensor grid), masses ε, from α = 0.

A. The k-fold near tie. One feature of amplitude A near a tie of the kernels' square lattice of spacing H: an edge (the kernels at (0, 0)
   and (H, 0), the feature near (H/2, 0): a 2-fold tie) and a vertex ((0, 0), (H, 0), (0, H), (H, H), the feature near (H/2, H/2): a
   4-fold tie), placed so that the runner-ups' margins when the entrant enters are μ_j = ε^{1/2}x_j for the given gaps (the diagonal's
   gap then follows: x₁₁ ≈ x₁₀ + x₀₁). The loop to TAU_A τ_entry at ε = EPS_A; the shares |w_j|/Σ|w| against the ε-free inner problem
   with the k kernels at the same gaps (the feature at the tie; `inner_k`), the softmax e^{−κx_j}/Σe^{−κx_i} (κ = πρ(2c)^{1/2}) and, at
   the vertex, the product of the two edges' inner problems (the factorisation of the 4-fold split).
B. The formation map at a vertex: the number of kernels formed (|w| > TH_REL·|A|) over gaps (x₁₀, x₀₁) on a grid, the loop at EPS_B
   against the 4-kernel inner problem and the independent rule #{j: x_j < x_θ} (x_θ: the 2-fold inner problem's crossing of TH_REL);
   the probabilities for a feature placed uniformly in a Voronoi cell, P(N ≥ 2) ∝ ε^{1/2}, P(N ≥ 3) and P(N = 4) ∝ ε, from the map's areas.
C. Several features in the plane: random three-feature targets on the lattice over [−1, 1]², the sliding path against the loop at ε =
   EPS_C: the kernels formed, the never-entering kernels' largest coefficient (∝ ε^{3/2}) and displacement (∝ ε), the block against the
   path's.
T. The 2-kernel (edge) inner problem's runner-up share p(x; t) on a grid of gaps x and of the feature's offsets t along the face
   (referee round 42: the split depends on where the feature sits on the face, not only on the gap) — x_θ(t) and its average over the
   face; the product rule: at a vertex, where each edge's offset is t = H/2, the kernel e + h1_S has the share
   Π_{i∈S}p(x_i; H/2)Π_{i∉S}(1 − p(x_i; H/2)), formed when it exceeds TH_REL —, used by A, B's map and D.
B4. The 4-kernel inner problem on a finer grid of gaps (GRID_B4): the areas of the regions N ≥ 3 and N = 4 at a vertex, against the
   product rule's.
D. The probabilities: features placed uniformly in a Voronoi cell (N_D positions) at ε = EPS_D, the loop's count against the product
   rule at each position's gaps and offsets and against the asymptotic P(N ≥ 2) = (4/h)·∫w(t)dt/h, w(t) = 2σ²x_θ(t)ε^{1/2}/h.
H. A hexagonal lattice (spacing H; a Voronoi vertex meets three cells, so k ≤ n + 1 = 3): the 3-kernel inner problem at a vertex's exact
   tie, its formation region N ≥ 3 in the two runner-ups' gaps and P(N ≥ 3) ∝ ε from it, the fourth kernel's gap there, and the loop at
   the exact tie against the inner problem.
T2. x_θ(t) at the offsets T_FAR beyond H/2 (by bisection): the faces of hℤⁿ for n ≥ 3 reach |t| = (n − 1)^{1/2}H/2, and in the
   continuum the 2-kernel problem with the feature at the transverse offset t reduces to the plane of the axis and t — x_θ depends on
   |t| only —, so the face averages x̄_θ⁽ⁿ⁾ over the (n − 1)-cube of offsets give E[N − 1] on hℤⁿ.
W. Where the selection is the Voronoi cell (referee round 42): W1 the samples' window (Part C's grid: the faces between neighbours on
   the axis x₂ = 0 and the margins' error against the continuum's toward the window's edge); W2 unequal masses (the kernel at (H, 0) of
   mass βε: the entrant the kernel of the least entry time ρ_j/v_j — a power diagram —, the loop at EPS_W); W3 coarse samples along x₂
   (spacing DS_W: the face's shift against the Poisson summation's harmonics); W4 two coinciding kernels with the feature off their
   centre (the inner problem against the profile split); W5 the profile split against the softmax at the formation threshold.
C2. Designed targets in the plane: a near tie (one feature near a face of the lattice, its runner-up's margin MU_C2 at the first entry:
   the gap MU_C2/ε^{1/2} crosses x_θ between the masses EPS_C) and a feature at the window's edge (the samples' entrant is not the
   nearest kernel) — the loop against the path and the near-tie rule.

Usage: python experiments/exp_transit_nd.py [--only A,B,B4,C,T,D,H,T2,W,C2] [--from-raw] → results/chaos/transit_nd.md (+ _raw.npz).
"""
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

OUT = os.environ.get('TN_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TND_WORKERS', '2'))
SIG, H, DS = 0.15, 0.125, 0.05                                        # the kernels' width, the lattice's spacing, the samples' spacing
TH_REL = 1e-2                                                          # formed: |w| > TH_REL·|A|
TAU_A = 1.32                                                           # the loop to 1.32 τ_entry (the relaxation after the entry)
AMPS_A = (0.3, 1.0)
EPS_A = (1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
EDGE_X = (0.0, 0.5, 1.0, 2.0, 4.0, 8.0)
VERTEX_X = ((0.0, 0.0), (0.5, 0.5), (1.0, 1.0), (1.0, 3.0), (2.0, 2.0), (2.0, 5.0), (4.0, 4.0))
A_B, EPS_B, GRID_B = 0.8, (1e-6, 1e-8, 1e-10), np.arange(0.0, 10.5, 1.0)
SEED_C, N_C, EPS_C = 54, 6, (1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
X_T = np.round(np.arange(0.0, 14.05, 0.1), 10)                         # the 2-kernel (edge) inner problem's share, tabulated
SEED_D, N_D, EPS_D = 541, 600, (1e-5, 1e-6)
T_OFFS = (0.0, H / 8, H / 4, 3 * H / 8, H / 2)                          # T: the feature's offsets along the edge's face
XI0_A = (800.0, 1600.0)                                                # A: the inner problems' starts, extrapolated (the tail ∝ 1/xi0)
GRID_B4 = np.arange(0.0, 10.01, 0.25)                                  # B4: the 4-kernel inner problem's finer grid
GRID_H, EPS_H = np.arange(0.0, 12.01, 0.5), (1e-8, 1e-10)              # H: the hexagonal vertex's gaps; the loop's masses
SHARE_TOL = 1e-7                                                       # the inner problems' shares settled over a doubling
INNER_SECS = 300.0                                                     # an inner problem's wall-clock budget (then 'not stationary')
T_FAR = (0.625 * H, 0.75 * H, 0.875 * H, H)                            # T2: the offsets beyond H/2 (n ≥ 3)
EPS_W = 1e-8                                                           # W1/W2: the margins' scale ε^{1/2}x_θ and the loop's mass
BETAS_W = (1.0, 1.001, 1.01, 1.5)                                      # W2: the mass ratio of the kernel at (H, 0)
DS_W = (0.05, 0.1, 0.15, 0.2, 0.3)                                     # W3: the samples' spacing along x₂ (offset by 0.3 of a spacing)
OFF_W, XR_W = (0.0, 0.07, 0.2, 0.4), (2.0, 4.0)                        # W4: the feature's offset (units of σ), the runner-up's gap
AMPS_W = (0.3, 0.8, 1.0)                                               # W5
MU_C2 = 3e-4                                                           # C2: the near tie's margin (the gap 3 at ε = 10⁻⁸)
D_V, AP_V, D_V2 = (5.0, 5.5, 6.0, 6.5, 7.0), 0.7, 6.0                   # V: the unfitted feature's distances (σ), amplitude; the planar check's
W_V = np.geomspace(0.3, 1e-4, 49)                                      # V: the tail's coefficients, descending (continuation)
EPS_C3 = tuple(float(f'{e:.3g}') for e in np.geomspace(1.5e-10, 1e-8, 13)) + (3e-9, 3.13e-10)  # C3: the capture scan's masses
RTOL_C3 = 1e-12                                                        # C3: the loop's tolerance (captures sit on a saddle: round 47)
T3_X = ((1.0, 1.0, 1.0), (2.0, 2.0, 2.0), (1.0, 2.0, 3.0), (0.5, 3.0, 5.0), (3.0, 3.0, 3.0), (2.0, 4.0, 6.0), (4.0, 4.0, 4.0),
        (5.9, 9.0, 9.0), (5.95, 9.0, 9.0), (6.0, 9.0, 9.0))                # the last three: a corner between x_θ(H/√2) and x_θ(H/2)
X_C4 = (3.0, 4.9, 7.0, 10.0, 17.3, 24.5)                                 # C4: the near tie's gap held fixed while ε varies (round 49)
EPS_C4 = (1e-8, 1e-9, 1e-10, 1e-11, 1e-12)
RTOL_C4B = 1e-13                                                        # C4b: the race again at a tighter tolerance (a saddle's passage)
N_OUT_C4 = 160                                                           # C4: outputs from feature 0's entry to feature 2's, log-spaced
MU_R, XI0_R, RTOL_R, S_END_R, SIDE_R = 1e-9, 400.0, 1e-13, 3e8, 0.05                  # R: the ε-free race — the geometry's margin (feature 0 on the face
                                                                       # as ε → 0 at a fixed gap), the start, the tolerance, the horizon in s
X_R_EXTRA = (2.5, 4.0, 4.5, 30.0)                                       # R: gaps beyond C3's and C4's (30: C2's far tail)
S_MARKS_R = (1e4, 1e5, 1e6, 1e7, 1e8, 3e8)                              # R: the runner-up's coefficient and distance at these s
X_RM = (3.0, 4.9, 7.0, 10.0, 17.3)                                      # R: the mirrored target's gaps (the far side)
X_RS, RTOLS_RS, S_END_RS = (2.5, 3.0, 3.5, 4.0, 4.5, 4.9, 5.5, 6.0), (1e-9, 1e-10, 1e-11, 1e-12, 1e-13), 3e6   # R: the side's tolerance scan
C4M = ((7.0, 1e-13), (7.0, 1e-10), (7.0, 1e-17), (7.0, 2e-12), (7.0, 4.5e-12), (10.0, 2.5e-13), (10.0, 6e-13), (7.0, 6.5e-17),
       (7.0, 3.7e-17))  # C4M: the loop on the
# mirrored target — the far rule's boundary s₂ ≈ s_θ and beyond s_{d,e} (the swap, referee round 58)
R_LONG = ((7.0, 3e9), (10.0, 3e9))                                    # R: the mirror past its horizon — the swap and the entrant's race
X_RP, TOPS_RP = (3.0, 4.5, 4.9, 5.07, 5.48, 6.04, 7.0), (0.95, 1.0, 1.5)  # RP: the gaps (held and formed runner-ups, a tail); the top edge
S_END_RP, S0_RP, S_MARKS_RP = 5e7, 400.0, (500.0, 700.0, 1000.0)        # RP: the run's end; the law from s₀; D reported at these s
WIN_FAR, WIN_RIGHT = (-2.0, 2.0, -1.5, 2.0), (-1.0, 0.8, -1.0, 1.5)      # RP (referee round 62): (x₋, x₊, y₋, y₊) — every edge ≥ 11.67σ
                                                                        # from feature 0 (the continuum); the right edge 4.92σ from it
D_CLIMB, S_JAC_RP = 1e-2, (450.0, 5e3)                                  # RP: the climb's end |D| = 1e-2σ; the Jacobian's checks at s
X_RC, S1_RC, SUB_RC, INJ_RC = X_RP, 400.0, 4, 1e-9                         # RC (referee round 66): the gaps; the merged block from s₁; the matrix
S_CONT_RC, N_CONT_RC, S_END_RC = 2e6, 800, 2e8                           # exponentials per output; the injection's |D| (σ); the continuation's
                                                                        # end and outputs, the injected race's end (WIN_FAR's time)
RC_CHECKS = (('near tie', 4.9, 16, 'y', 1e-9), ('near tie', 4.9, 4, 'zero', 1e-9), ('near tie', 4.9, 4, 'y', 1e-7),
             ('near tie', 4.9, 4, 'y', 1e-11))                         # RC: 16 substeps, y = 0, the injection at 10⁻⁷ and 10⁻¹¹σ
THETAS_CLIMB = (1e-3, 1e-2, 1e-1)                                       # RC, RP (referee round 72): the climb's end at these |D| (σ)
RD_JOBS, S_RD, RD_PREC = (('near tie', 3.0), ('near tie mirrored', 3.0)), (400.0, 4e3, 4e4), 80   # RD (referee round 72): RC's merged
DS_DYADIC = 3 / 64                                                      # block at these s; exact sums to RD_PREC digits; a dyadic grid
W_FC = np.geomspace(0.06, 0.003, 120)                                   # FC: the configuration's own fold (continuation in w)
D_V3, AP_V3 = (5.0, 6.0, 7.0), (0.35, 0.7, 1.4, 2.8)                     # V3: the fold against A′/A (A = A_B)
MAX_DOUBLINGS_B4X = 24                                                   # B4x: B4's unsettled problems continued
EPS_D2, T_D2 = (1e-5, 1e-6, 1e-8), (0.0, 3 * H / 8)                       # D2: the loop's threshold at the offsets 0 and 3H/8


def law():
    warnings.filterwarnings('ignore')
    return X.law()


def local():
    """One feature: the samples symmetric about the vertex (H/2, H/2) (so that the ties are exact), the kernels on the lattice H·(i, j),
    i, j = −3 … 4, around it."""
    xs = H / 2 + DS * np.arange(-18, 19)
    K = RN.BumpsGrid2(xs, xs, SIG)
    g = H * np.arange(-3, 5)
    return K, np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)


def kidx(c0, p):
    return int(np.argmin(((c0 - np.asarray(p, float)) ** 2).sum(1)))


TIES = {'edge': ([(0, 0), (H, 0)], np.array([H / 2, 0.0])),
        'vertex': ([(0, 0), (H, 0), (0, H), (H, H)], np.array([H / 2, H / 2]))}


def margins(K, c0, A, cstar, idx):
    """The margins μ_j = 1 − v_j/v_e (v = |⟨fφ_j⟩| at the start, e the largest) of the kernels idx (the first the entrant)."""
    v = np.abs(K.corr(c0[idx], A * K.phi(np.atleast_2d(cstar))[:, 0])[0])
    return 1 - v / v[0], v[0]


def place(K, c0, A, kind, x, eps, t=0.0):
    """The feature's centre near the tie with the runner-ups' margins ε^{1/2}x (edge: x₁₀, the feature at the offset t along the face;
    vertex: x₁₀, x₀₁ — the diagonal's follows), by Newton's method on the offset from the tie toward the entrant."""
    pts, T = TIES[kind]
    T = T + (np.array([0.0, t]) if kind == 'edge' else 0.0)
    idx = [kidx(c0, p) for p in pts]
    m = 1 if kind == 'edge' else 2
    tgt = np.sqrt(eps) * np.asarray(x, float)[:m]
    d = np.full(m, 1e-12) + tgt * 2 * SIG ** 2 / H                    # the continuum's linear map as the start
    for _ in range(30):
        def mu(dd):
            off = np.zeros(2); off[:m] = dd
            return margins(K, c0, A, T - off, idx)[0][1:1 + m]
        f0 = mu(d) - tgt
        J = np.empty((m, m))
        for i in range(m):
            e = np.zeros(m); e[i] = 1e-8
            J[:, i] = (mu(d + e) - mu(d - e)) / (2 * e[i])
        step = np.linalg.solve(J, f0)
        d = d - step
        if np.all(np.abs(step) <= 1e-15 + 1e-12 * np.abs(d)):
            break
    off = np.zeros(2); off[:m] = d
    return T - off, idx


def softmax(x, L):
    kap = np.pi * (2 * L.c) ** 0.5
    s = np.exp(-kap * (np.asarray(x, float) - np.min(x)))
    return s / s.sum()


def a_job(args):
    """One loop: (kind, A, x, ε) → the shares of the tie's kernels, their margins in units of ε^{1/2}, the others' largest |w|."""
    kind, A, x, eps = args
    t0 = time.time()
    L = law(); K, c0 = local()
    cstar, idx = place(K, c0, A, kind, x, eps)
    mu, v_e = margins(K, c0, A, cstar, idx)
    f = A * K.phi(np.atleast_2d(cstar))[:, 0]
    taus = np.linspace(0, TAU_A / v_e, 23)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus)
    w = np.abs(W[-1][idx])
    rest = np.delete(np.abs(W[-1]), idx)
    return dict(kind=kind, A=A, x=list(x), eps=eps, cstar=cstar, x_act=mu / np.sqrt(eps), shares=w / w.sum(), w=W[-1][idx],
                c=C[-1][idx], rest=float(rest.max()), n_formed=int((np.abs(W[-1]) > TH_REL * abs(A)).sum()), secs=time.time() - t0)


def inner_job(args):
    """The ε-free inner problem at the tie: (kind, A, gaps of the tie's kernels[, the feature's offset t along the edge's face[, xi0]])
    → shares (the edge's feature at (H/2, t): still on the face, so still a tie)."""
    kind, A, xs = args[:3]
    t = float(args[3]) if len(args) > 3 else 0.0
    xi0 = float(args[4]) if len(args) > 4 else 400.0
    L = law(); K, c0 = local()
    pts, T = TIES[kind]
    T = T + (np.array([0.0, t]) if kind == 'edge' else 0.0)
    f = A * K.phi(np.atleast_2d(T))[:, 0]
    out = RN.inner_k(L, K, f, np.array(pts, float), np.asarray(xs, float), xi0=xi0, share_tol=SHARE_TOL, max_secs=INNER_SECS)
    return dict(kind=kind, A=A, xs=list(xs), t=t, xi0=xi0, shares=out['shares'], stationary=out['stationary'], s=out['s'], c=out['c'],
                drift=out['drift'], stop=out['stop'], remaining=out['remaining'])


def inner_ext(args):
    """The inner problem's shares extrapolated in its start xi0 (Richardson from XI0_A: the far field's tail ∝ 1/xi0; referee round 42)."""
    kind, A, xs = args[:3]
    t = args[3] if len(args) > 3 else 0.0
    r = [inner_job((kind, A, xs, t, x0)) for x0 in XI0_A]
    out = dict(r[-1])
    out['shares'] = 2 * r[-1]['shares'] - r[-2]['shares']
    out['spread'] = float(np.abs(r[-1]['shares'] - r[-2]['shares']).max())
    out['stationary'] = all(x['stationary'] for x in r)
    return out


def part_A():
    from multiprocessing import Pool
    loops = [('edge', A, (x,), e) for A in AMPS_A for x in EDGE_X for e in EPS_A]
    loops += [('vertex', A, x, e) for A in AMPS_A for x in VERTEX_X for e in EPS_A]
    with Pool(N_WORKERS) as pool:
        rl = pool.map(a_job, loops, chunksize=1)
        # the inner problem at the loop's gaps at the smallest ε (the limit's), extrapolated in xi0, and for the factorisation the edges'
        # problems in the vertex's geometry: the feature at the vertex sits H/2 off each edge's axis (referee round 42)
        jobs = []
        for r in rl:
            if r['eps'] == min(EPS_A):
                jobs.append((r['kind'], r['A'], tuple(np.r_[0.0, r['x_act'][1:]]), 0.0))
                if r['kind'] == 'vertex':
                    jobs += [('edge', r['A'], (0.0, r['x_act'][1]), H / 2), ('edge', r['A'], (0.0, r['x_act'][2]), H / 2)]
        ri = pool.map(inner_ext, jobs, chunksize=1)
    return dict(loops=rl, inner=ri)


def crossing(share, lo, hi, it):
    """The gap in [lo, hi] at which a share falling with the gap crosses TH_REL, by bisection (it halvings)."""
    for _ in range(it):
        m_ = 0.5 * (lo + hi)
        if share(m_) > TH_REL:
            lo = m_
        else:
            hi = m_
    return 0.5 * (lo + hi)


def threshold_x(L, K, A, lo=0.0, hi=14.0, it=30):
    """x_θ: the gap at which the 2-fold (edge) inner problem's runner-up share crosses TH_REL."""
    return crossing(lambda x: inner_job(('edge', A, (0.0, x)))['shares'][1], lo, hi, it)


def b_job(args):
    x1, x2, eps = args
    r = a_job(('vertex', A_B, (x1, x2), eps))
    return dict(x=(x1, x2), eps=eps, n=r['n_formed'], x_act=r['x_act'], shares=r['shares'])


def b_inner(args):
    x1, x2 = args
    L = law(); K, c0 = local()
    xd = margins_diag(K, c0, x1, x2)
    r = inner_job(('vertex', A_B, (0.0, x1, x2, xd)))
    return dict(x=(x1, x2), xd=xd, n=int((r['shares'] > TH_REL).sum()), shares=r['shares'], stationary=r['stationary'])


def margins_diag(K, c0, x1, x2, eps=1e-10):
    """The diagonal's gap at the vertex for the edge gaps (x₁₀, x₀₁) (at a small ε, where the map is linear)."""
    cstar, idx = place(K, c0, A_B, 'vertex', (x1, x2), eps)
    return float(margins(K, c0, A_B, cstar, idx)[0][3] / np.sqrt(eps))


def part_B():
    from multiprocessing import Pool
    pts = [(a, b) for a in GRID_B for b in GRID_B if a <= b]
    with Pool(N_WORKERS) as pool:
        rl = pool.map(b_job, [(a, b, e) for e in EPS_B for a, b in pts], chunksize=1)
        ri = pool.map(b_inner, pts, chunksize=1)
    L = law(); K, c0 = local()
    return dict(loops=rl, inner=ri, x_theta=threshold_x(L, K, A_B))


def b4_job(args):
    """B4: the 4-kernel inner problem at a vertex at the gaps (x₁₀, x₀₁) (the diagonal's from the samples' margins)."""
    x1, x2 = args
    L = law(); K, c0 = local()
    xd = margins_diag(K, c0, x1, x2)
    t0 = time.time()
    r = inner_job(('vertex', A_B, (0.0, x1, x2, xd)))
    return dict(x=(x1, x2), xd=xd, n=int((r['shares'] > TH_REL).sum()), shares=r['shares'], stationary=r['stationary'],
                stop=r.get('stop'), drift=r.get('drift'), remaining=r.get('remaining'), s=r.get('s'), secs=time.time() - t0)


def part_B4u(B4):
    """B4u: B4's problems that did not settle, rerun with their stop's reason recorded (referee round 47: the doubling cap, not the
    clock)."""
    from multiprocessing import Pool
    pts = [tuple(r['x']) for r in B4['rows'] if not r['stationary']]
    with Pool(N_WORKERS) as pool:
        return dict(rows=pool.map(b4_job, pts, chunksize=4))


def part_B4():
    from multiprocessing import Pool
    pts = [(a, b) for a in GRID_B4 for b in GRID_B4 if a <= b]
    with Pool(N_WORKERS) as pool:
        return dict(grid=GRID_B4, rows=pool.map(b4_job, pts, chunksize=4))


def local_hex():
    """H: the kernels on a hexagonal lattice of spacing H around the vertex V of the cells of (0, 0), (H, 0), (H/2, √3H/2) — the triangle's
    circumcentre, H/√3 from each —, the samples' tensor grid of spacing DS centred at V (no threefold symmetry: the tie is exact in the
    continuum, approximately on the samples)."""
    V = np.array([H / 2, H / (2 * np.sqrt(3))])
    xs1, xs2 = V[0] + DS * np.arange(-18, 19), V[1] + DS * np.arange(-18, 19)
    K = RN.BumpsGrid2(xs1, xs2, SIG)
    pts = np.array([[H * (i + j / 2), H * j * np.sqrt(3) / 2] for i in range(-5, 6) for j in range(-5, 6)])
    c0 = pts[np.linalg.norm(pts - V, axis=1) <= 4 * H]
    tri = [kidx(c0, p_) for p_ in ((0, 0), (H, 0), (H / 2, H * np.sqrt(3) / 2))]
    return K, c0, V, tri


def hex_tie(K, c0, V, tri):
    """The point near V where the three kernels' correlations with the feature are equal on the samples (Newton on the two margins)."""
    c = V.copy()
    for _ in range(30):
        def mu(p_):
            v = np.abs(K.corr(c0[tri], K.phi(np.atleast_2d(p_))[:, 0])[0])
            return np.array([v[1] / v[0] - 1, v[2] / v[0] - 1])
        f0 = mu(c)
        J = np.empty((2, 2))
        for i in range(2):
            e = np.zeros(2); e[i] = 1e-7
            J[:, i] = (mu(c + e) - mu(c - e)) / 2e-7
        step = np.linalg.solve(J, f0)
        c = c - step
        if np.abs(step).max() < 1e-14:
            break
    return c


def h_inner(args):
    """H: the 3-kernel inner problem at the hexagonal vertex (the samples' exact tie): the runner-ups x₁, x₂ behind the entrant."""
    x1, x2 = args
    L = law(); K, c0, V, tri = local_hex()
    f = A_B * K.phi(np.atleast_2d(hex_tie(K, c0, V, tri)))[:, 0]
    out = RN.inner_k(L, K, f, c0[tri], np.array([0.0, x1, x2]), share_tol=SHARE_TOL, max_secs=INNER_SECS)
    return dict(x=(x1, x2), shares=out['shares'], n=int((out['shares'] > TH_REL).sum()), stationary=out['stationary'])


def h_loop(eps):
    """H: the loop with the feature at the hexagonal vertex (the exact tie in the continuum): the three kernels' shares, the others' largest
    |w| and the margins of the tie's runner-ups and of the fourth-nearest kernel, in units of ε^{1/2}."""
    L = law(); K, c0, V, tri = local_hex()
    T0 = hex_tie(K, c0, V, tri)
    f = A_B * K.phi(np.atleast_2d(T0))[:, 0]
    v = np.abs(K.corr(c0, f)[0])
    order = np.argsort(-v)
    taus = np.linspace(0, TAU_A / v[order[0]], 23)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus)
    w = np.abs(W[-1][tri])
    return dict(eps=eps, tie=T0.tolist(), V=V.tolist(), shares=w / w.sum(), rest=float(np.delete(np.abs(W[-1]), tri).max()),
                n=int((np.abs(W[-1]) > TH_REL * A_B).sum()), gaps=((1 - v[order[1:4]] / v[order[0]]) / np.sqrt(eps)).tolist(),
                top3=sorted(int(j) for j in order[:3]) == sorted(tri))


def part_H():
    from multiprocessing import Pool
    pts = [(a, b) for a in GRID_H for b in GRID_H if a <= b]
    with Pool(N_WORKERS) as pool:
        ri = pool.map(h_inner, pts, chunksize=4)
        rl = pool.map(h_loop, list(EPS_H), chunksize=1)
    return dict(inner=ri, loops=rl, grid=GRID_H)


def targets_C(seed, n):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < n:
        cs = rng.uniform(-0.7, 0.7, (3, 2))
        if min(np.linalg.norm(cs[i] - cs[j]) for i in range(3) for j in range(i)) < 3 * SIG:
            continue
        out.append((cs, rng.choice([-1.0, 1.0], 3) * rng.uniform(0.4, 1.0, 3)))
    return out


def full():
    xs = np.arange(-1, 1 + 1e-9, DS)
    K = RN.BumpsGrid2(xs, xs, SIG)
    g = np.arange(-1, 1 + 1e-9, H)
    return K, np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)


def c_job(args):
    q, eps = args
    t0 = time.time()
    L = law(); K, c0 = full()
    cs, As = targets_C(SEED_C, N_C)[q]
    f = K.synth(cs, As)
    R = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)))
    path = sorted(int(j) for j in np.flatnonzero(R['act']))
    t_last = R['events'][-1][0]
    taus = np.linspace(0, 1.3 * t_last, 27)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus)
    near = np.argmin(np.linalg.norm(c0[:, None, :] - cs[None, :, :], axis=2), axis=1)
    big = np.abs(W[-1]) > TH_REL * np.abs(As)[near]                    # formed: against its nearest feature's amplitude (round 42)
    never = np.array([j for j in range(len(c0)) if j not in path])
    rms = float(np.sqrt(np.mean((f - K.synth(C[-1], W[-1])) ** 2)))
    # the creep as a vector law (D‴.48(iii)(a) in ℝⁿ): the never-entering kernels' displacements against ε·C_j(τ_end), C_j integrated
    # along the path (`creep_constants_nd`); kernels whose predicted displacement is above 10⁻³ of the largest
    P = RN.path_duals_nd(K, f, R)
    Cj = np.array([RN.creep_constants_nd(L, K, P, j, c0[j], [taus[-1]])[0][0] for j in never])
    disp = C[-1][never] - c0[never]
    nrm = np.linalg.norm(Cj, axis=1)
    sel = nrm > 1e-3 * nrm.max()
    rel = float((np.linalg.norm(disp[sel] - eps * Cj[sel], axis=1) / (eps * nrm[sel])).max())
    cos = float(((disp[sel] * Cj[sel]).sum(1) / (np.linalg.norm(disp[sel], axis=1) * nrm[sel])).min())
    return dict(q=q, eps=eps, path=path, events=[(float(t), int(j), k) for t, j, k, c in R['events']], path_rms=R['hist'][-1][3],
                formed=sorted(int(j) for j in np.flatnonzero(big)), w_never=float(np.abs(W[-1][never]).max()),
                disp_never=float(np.linalg.norm(C[-1][never] - c0[never], axis=1).max()),
                dw=float(np.abs(W[-1][path] - R['w'][path]).max()), dc=float(np.abs(C[-1][path] - R['c'][path]).max()), rms=rms,
                creep_rel=rel, creep_cos=cos, secs=time.time() - t0)


def part_C():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, [(q, e) for q in range(N_C) for e in EPS_C], chunksize=1)


def t_job(args):
    x, t = args
    return float(inner_job(('edge', A_B, (0.0, float(x)), float(t)))['shares'][1])


def part_T():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        p = np.array(pool.map(t_job, [(x, t) for t in T_OFFS for x in X_T], chunksize=4)).reshape(len(T_OFFS), len(X_T))
    return dict(x=X_T, t=np.array(T_OFFS), pt=p, p=p[0])


def p_edge(T, x, t=H / 2):
    """The edge's runner-up share p(x; t) from the table T (linear in x and in |t| ∈ [0, H/2]; the old table: t = 0 only)."""
    if 'pt' not in T:
        return np.interp(x, T['x'], T['p'], right=0.0)
    ta = np.clip(np.abs(t), 0.0, H / 2)
    rows = np.array([np.interp(x, T['x'], T['pt'][i], right=0.0) for i in range(len(T['t']))])
    return np.interp(ta, T['t'], rows) if np.ndim(x) == 0 else np.array([np.interp(ta, T['t'], rows[:, i]) for i in range(len(x))])


def x_theta_t(T):
    """x_θ(t): the gap at which p(x; t) crosses TH_REL, per offset of the table."""
    return np.array([float(np.interp(-TH_REL, -T['pt'][i], T['x'])) for i in range(len(T['t']))])


def product_count(T, x1, x2, t1=H / 2, t2=H / 2):
    """The product rule's count with the two faces' gaps x1, x2 and the feature's offsets t1, t2 along them (at a vertex both H/2)."""
    p1, p2 = float(p_edge(T, x1, t1)), float(p_edge(T, x2, t2))
    shares = [(1 - p1) * (1 - p2), p1 * (1 - p2), (1 - p1) * p2, p1 * p2]
    return int(sum(s > TH_REL for s in shares))


def d_job(args):
    """A feature at a uniform position in the Voronoi cell of (0, 0): the loop's count and the gaps to the two nearest faces."""
    q, eps = args
    cstar = np.random.default_rng(SEED_D + q).uniform(-H / 2, H / 2, 2)
    L = law(); K, c0 = local()
    f = A_B * K.phi(np.atleast_2d(cstar))[:, 0]
    v = np.abs(K.corr(c0, f)[0])
    e = int(np.argmax(v))
    taus = np.linspace(0, TAU_A / v[e], 23)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus)
    gaps = []
    for ax in range(2):                                                # the nearer face along each axis
        nb = [kidx(c0, c0[e] + s * H * np.eye(2)[ax]) for s in (-1, 1)]
        gaps.append(min((1 - v[j] / v[e]) / np.sqrt(eps) for j in nb))
    return dict(q=q, eps=eps, cstar=cstar, e=e, n=int((np.abs(W[-1]) > TH_REL * A_B).sum()), gaps=gaps)


def part_D():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(d_job, [(q, e) for e in EPS_D for q in range(N_D)], chunksize=4)


def t2_job(t):
    """T2: x_θ(t), the gap at which the edge's runner-up share (A_B) crosses TH_REL, by bisection, the feature at the offset t."""
    return crossing(lambda x: inner_job(('edge', A_B, (0.0, x), t))['shares'][1], 3.0, 9.0, 16)


def part_T2():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return dict(t=np.array(T_FAR), x=np.array(pool.map(t2_job, list(T_FAR), chunksize=1)))


def face_average(tt, xt, n, m=24):
    """x̄_θ⁽ⁿ⁾: x_θ(|t|) (interpolated in |t| on the offsets tt) averaged over a face of hℤⁿ, the (n − 1)-cube [−H/2, H/2]^{n−1} of
    transverse offsets (the midpoint rule, m points per side)."""
    if n == 1:
        return float(xt[0])
    g = (np.arange(m) + 0.5) / m * H - H / 2
    pts = np.stack(np.meshgrid(*([g] * (n - 1)), indexing='ij'), -1).reshape(-1, n - 1)
    return float(np.interp(np.linalg.norm(pts, axis=1), tt, xt).mean())


def w_edge():
    """W1: the selection near the samples' window (Part C's grid, the samples on [−1, 1]²). On the axis x₂ = 0: the faces between lattice
    neighbours (where their correlations with a feature are equal on the samples; in the continuum the midpoints), searched up to 2H
    beyond the outer kernel; and for features on the axis the largest error of the margins of the nearest kernel's neighbours against the
    continuum's, and whether the samples' entrant is the nearest kernel."""
    K, c0 = full()
    row = [kidx(c0, (x, 0.0)) for x in np.arange(0.0, 1.0 + 1e-9, H)]

    def gap(x, a, b):
        v = np.abs(K.corr(c0[[a, b]], K.phi(np.array([[x, 0.0]]))[:, 0])[0])
        return v[0] - v[1]
    faces = []
    for a, b in zip(row[:-1], row[1:]):
        g = np.linspace(c0[a][0], c0[b][0] + 2 * H, 241)
        s = np.array([gap(x, a, b) for x in g])
        k = np.flatnonzero((s[:-1] > 0) & (s[1:] <= 0))
        face = np.nan
        if k.size:
            lo, hi = g[k[0]], g[k[0] + 1]
            for _ in range(50):
                m_ = 0.5 * (lo + hi)
                lo, hi = (m_, hi) if gap(m_, a, b) > 0 else (lo, m_)
            face = 0.5 * (lo + hi)
        faces.append(dict(a=float(c0[a][0]), b=float(c0[b][0]), face=float(face)))
    rows = []
    for x in np.round(np.arange(0.0, 0.985, 0.01), 10):
        cs = np.array([x, 0.0])
        v = np.abs(K.corr(c0, K.phi(np.atleast_2d(cs))[:, 0])[0])
        en = int(np.argmin(((c0 - cs) ** 2).sum(1)))
        nb = [j for j in range(len(c0)) if j != en and np.linalg.norm(c0[j] - c0[en]) < 1.5 * H]
        mu_s = 1 - v[nb] / v[en]
        mu_c = 1 - np.exp(-(((cs - c0[nb]) ** 2).sum(1) - ((cs - c0[en]) ** 2).sum()) / (4 * SIG ** 2))
        rows.append(dict(x=float(x), nearest=bool(int(np.argmax(v)) == en), err=float(np.abs(mu_s - mu_c).max())))
    return dict(faces=faces, rows=rows)


def power_face(beta):
    """W2: the point on the axis x₂ = 0 where the kernels at (0, 0) and (H, 0) (the latter of mass βε) enter together (equal entry times
    ρ_j/v_j on the local samples, ρ = b^{−1/2}); in the continuum H/2 − σ²ln β/H."""
    K, c0 = local()
    e, j = kidx(c0, (0, 0)), kidx(c0, (H, 0))

    def d(x):
        v = np.abs(K.corr(c0[[e, j]], K.phi(np.array([[x, 0.0]]))[:, 0])[0])
        return v[1] / beta ** -0.5 - v[0]                              # > 0: the heavier kernel enters first
    lo, hi = -2 * H, H
    for _ in range(60):
        m_ = 0.5 * (lo + hi)
        lo, hi = (lo, m_) if d(m_) > 0 else (m_, hi)
    return 0.5 * (lo + hi)


def w_mass(args):
    """W2: one loop (ε = EPS_W) on the local lattice, the kernel at (H, 0) of mass βε and the rest ε, the feature (A_B) on the axis at x₁:
    the kernels formed, the shares of (0, 0) and (H, 0), the entrant by the least entry time ρ_j/v_j (the power diagram) and the nearest
    kernel (Voronoi), and the runner-up's gap (1 − τ_e/τ_r)/ε^{1/2}."""
    beta, x1 = args
    L = law(); K, c0 = local()
    e, j = kidx(c0, (0, 0)), kidx(c0, (H, 0))
    b = np.ones(len(c0)); b[j] = beta
    cs = np.array([x1, 0.0])
    f = A_B * K.phi(np.atleast_2d(cs))[:, 0]
    tau = b ** -0.5 / np.abs(K.corr(c0, f)[0])
    o = np.argsort(tau)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, EPS_W, np.linspace(0, TAU_A * tau[o[0]], 23), b=b)
    w = np.abs(W[-1])
    return dict(beta=beta, x1=x1, e=e, j=j, formed=[int(i) for i in np.flatnonzero(w > TH_REL * A_B)], shares=(w[[e, j]] / w.sum()).tolist(),
                power=int(o[0]), runner=int(o[1]), voronoi=int(np.argmin(((c0 - cs) ** 2).sum(1))),
                gap=float((1 - tau[o[0]] / tau[o[1]]) / np.sqrt(EPS_W)))


def w_alias(ds):
    """W3: the samples' spacing ds along x₂ (the grid offset by 0.3ds; along x₁ DS): the face between the kernels at (0, 0) and (0, H) for
    features at (0, x₂) — its shift from H/2 against the Poisson summation (the x₂-sum of a product of Gaussians is (√πσ/ds)(1 +
    2Σ_k q^{k²}cos(2πk(m − o)/ds)), q = e^{−π²σ²/ds²}, m the product's centre, o the grid's offset) — and the margin's slope there against
    the continuum's H/(2σ²)."""
    n2 = int(np.ceil(1.0 / ds))
    o = 0.3 * ds
    K = RN.BumpsGrid2(DS * np.arange(-18, 19), o + ds * np.arange(-n2, n2 + 1), SIG)
    cc = np.array([[0.0, 0.0], [0.0, H]])
    q = float(np.exp(-np.pi ** 2 * SIG ** 2 / ds ** 2))

    def lr(x2):                                                        # log(v₀/v₁) on the samples
        v = np.abs(K.corr(cc, K.phi(np.array([[0.0, x2]]))[:, 0])[0])
        return float(np.log(v[0] / v[1]))

    def lr_p(x2):                                                      # the same from the Poisson sum
        S = lambda m: 1 + 2 * sum(q ** (k * k) * np.cos(2 * np.pi * k * (m - o) / ds) for k in (1, 2, 3))
        return (H * H - 2 * H * x2) / (4 * SIG ** 2) + np.log(S(x2 / 2)) - np.log(S((x2 + H) / 2))

    def root(fn):
        lo, hi = -H, 2 * H
        for _ in range(60):
            m_ = 0.5 * (lo + hi)
            lo, hi = (m_, hi) if fn(m_) > 0 else (lo, m_)
        return 0.5 * (lo + hi)
    face, face_p = root(lr), root(lr_p)
    h_ = 1e-6
    mu = lambda x2: 1 - np.exp(-lr(x2))                                # the runner-up (0, H)'s margin below the face
    slope = abs(mu(face + h_) - mu(face - h_)) / (2 * h_) / (H / (2 * SIG ** 2))
    return dict(ds=ds, q=q, shift=(face - H / 2) / H, shift_p=(face_p - H / 2) / H, slope=float(slope))


def w_coinc(args):
    """W4: two kernels coinciding at the local samples' centre of symmetry V, the runner-up x_r behind, the feature (A_B) at V + (δσ, 0):
    the inner problem's runner-up share (extrapolated in xi0) against the profile split's (exact for the feature at their centre; off it
    the kernels slide onto the feature at speeds ∝ their coefficients, apart on the way and together again at its centre)."""
    dlt, xr = args
    L = law(); K, c0 = local()
    V = TIES['vertex'][1]
    f = A_B * K.phi(np.atleast_2d(V + np.array([dlt * SIG, 0.0])))[:, 0]
    r = [RN.inner_k(L, K, f, np.array([V, V]), np.array([0.0, xr]), xi0=x0, share_tol=SHARE_TOL) for x0 in XI0_A]
    sh = 2 * r[-1]['shares'] - r[-2]['shares']
    prof = RN.profile_split(L, A_B, np.array([0.0, xr]))
    return dict(delta=dlt, xr=xr, share=float(sh[1]), profile=float(prof[1]), stationary=all(x['stationary'] for x in r),
                end=float(np.linalg.norm(r[-1]['c'] - (V + np.array([dlt * SIG, 0.0])), axis=1).max()) / SIG)


def w_profile():
    """W5: the 2-fold profile split's runner-up share p(x) (coinciding kernels at the feature's centre) against the softmax's
    1/(1 + e^{κx}): p's crossing of TH_REL against the softmax's ln(1/TH_REL − 1)/κ, the softmax's relative error at the crossing and at
    x = 8, and e^{−2α/a} at the runner-up's coefficient TH_REL·|A| (the size of the exponential branch's correction, D‴.51(i))."""
    L = law()
    kap = np.pi * (2 * L.c) ** 0.5
    soft = lambda x: 1 / (1 + np.exp(kap * x))
    out = []
    for A in AMPS_W:
        p = lambda x: float(RN.profile_split(L, A, np.array([0.0, x]))[1])
        lo, hi = 0.0, 14.0
        for _ in range(50):
            m_ = 0.5 * (lo + hi)
            lo, hi = (m_, hi) if p(m_) > TH_REL else (lo, m_)
        xc = 0.5 * (lo + hi)
        al = float(RT._alpha_of(L, np.array([TH_REL * A]), RT.coef_crit)[0])
        out.append(dict(A=A, x_cross=xc, x_soft=float(np.log(1 / TH_REL - 1) / kap), err_cross=soft(xc) / p(xc) - 1,
                        err8=soft(8.0) / p(8.0) - 1, alpha_a=al / L.A, corr=float(np.exp(-2 * al / L.A))))
    return out


def part_W(T=None):
    from multiprocessing import Pool
    L = law()
    xbar = x_theta_avg(T) if T is not None and 'pt' in T else threshold_x(L, None, A_B)
    w_ = 2 * SIG ** 2 * xbar * np.sqrt(EPS_W) / H                      # the near-tie layer's width at EPS_W (positions only)
    jobs = []
    for beta in BETAS_W:
        p = power_face(beta)
        xs = [H / 2, p - 2 * w_, p + 2 * w_] + ([p] if beta != 1.0 else [])
        jobs += [(beta, float(x)) for x in xs]
    with Pool(N_WORKERS) as pool:
        mass = pool.map(w_mass, jobs, chunksize=1)
        coinc = pool.map(w_coinc, [(d, x) for x in XR_W for d in OFF_W], chunksize=1)
        alias = pool.map(w_alias, list(DS_W), chunksize=1)
    return dict(edge=w_edge(), mass=mass, faces={b: power_face(b) for b in BETAS_W}, coinc=coinc, alias=alias, profile=w_profile(),
                layer=w_)


def x_theta_avg(T):
    """x̄_θ over a face of hℤ² from the T table (x_θ(t) at its offsets 0 … H/2, the trapezoid rule)."""
    tt = np.asarray(T['t'])
    return float(np.trapezoid(x_theta_t(T), tt) / tt[-1])


def targets_C2(mu_t=MU_C2):
    """C2's designed targets on Part C's lattice: 'near tie' — feature 0 (A_B) near the face between (0, 2H) and (H, 2H), on their axis,
    placed so that its runner-up's margin at the first entry (the target's correlations at the start: feature 0 enters first) is mu_t;
    'window edge' — feature 0 at (0.96, 0), nearest (1, 0), the samples' entrant (0.875, 0)."""
    K, c0 = full()
    oth = np.array([[-0.5, -0.5], [0.5, -0.55]]), np.array([-0.6, 0.7])
    e, r = kidx(c0, (0, 2 * H)), kidx(c0, (H, 2 * H))

    def mu(d):
        cs = np.r_[[[H / 2 + d, 2 * H]], oth[0]]
        v = np.abs(K.corr(c0[[e, r]], K.synth(cs, np.r_[A_B, oth[1]]))[0])
        return 1 - v[1] / v[0]
    d = mu_t * 2 * SIG ** 2 / H
    for _ in range(30):
        g = (mu(d + 1e-8) - mu(d - 1e-8)) / 2e-8
        step = (mu(d) - mu_t) / g
        d -= step
        if abs(step) < 1e-15:
            break
    c_nt = np.array([H / 2 + d, 2 * H])
    # the variant: the same-sign feature (0.7) moved to 5σ from feature 0 at −60°, feature 0 re-placed for the same margin
    oth5 = np.array([[-0.5, -0.5], c_nt + 5 * SIG * np.array([0.5, -np.sqrt(3) / 2])]), oth[1]

    def mu5(d5):
        cs = np.r_[[[H / 2 + d5, 2 * H]], oth5[0]]
        v = np.abs(K.corr(c0[[e, r]], K.synth(cs, np.r_[A_B, oth5[1]]))[0])
        return 1 - v[1] / v[0]
    d5 = d
    for _ in range(30):
        g = (mu5(d5 + 1e-8) - mu5(d5 - 1e-8)) / 2e-8
        step = (mu5(d5) - mu_t) / g
        d5 -= step
        if abs(step) < 1e-15:
            break
    # the mirror image of features 1 and 2 across the face x = H/2 — feature 2 on the entrant's side, so a runner-up displaced toward
    # the entrant's side of feature 0 sits on the far side of feature 2 (referee round 53) —, feature 0 re-placed for the same margin
    othm = np.array([[H - q[0], q[1]] for q in oth[0]]), oth[1]

    def mum(dm):
        cs = np.r_[[[H / 2 + dm, 2 * H]], othm[0]]
        v = np.abs(K.corr(c0[[e, r]], K.synth(cs, np.r_[A_B, othm[1]]))[0])
        return 1 - v[1] / v[0]
    dm = d
    for _ in range(30):
        g = (mum(dm + 1e-8) - mum(dm - 1e-8)) / 2e-8
        step = (mum(dm) - mu_t) / g
        dm -= step
        if abs(step) < 1e-15:
            break
    return {'near tie': (np.r_[[c_nt], oth[0]], np.r_[A_B, oth[1]]),
            'near tie 5σ': (np.r_[[[H / 2 + d5, 2 * H]], oth5[0]], np.r_[A_B, oth5[1]]),
            'near tie mirrored': (np.r_[[[H / 2 + dm, 2 * H]], othm[0]], np.r_[A_B, othm[1]]),
            'window edge': (np.array([[0.96, 0.0], [-0.4, 0.5], [-0.3, -0.5]]), np.array([0.8, -0.7, 0.9]))}


def formed_by_feature(c0, cs, As, w, c):
    """The bookkeeping of a loop's end: each kernel's feature at the start (the nearest to its starting point) and at the end (the
    nearest to its final centre), formed (|w| above TH_REL of the start's feature's amplitude), the counts formed per feature by
    the start and by the end, and the formed kernels that changed feature (captured: (kernel, start's feature, end's))."""
    near = np.argmin(np.linalg.norm(c0[:, None, :] - cs[None, :, :], axis=2), axis=1)
    fin = np.argmin(np.linalg.norm(np.asarray(c)[:, None, :] - cs[None, :, :], axis=2), axis=1)
    big = np.abs(w) > TH_REL * np.abs(As)[near]
    return dict(near=near, fin=fin, big=big, n_start=[int(big[near == k].sum()) for k in range(len(cs))],
                n_fin=[int(big[fin == k].sum()) for k in range(len(cs))],
                captured=[(int(j), int(near[j]), int(fin[j])) for j in np.flatnonzero(big) if near[j] != fin[j]])


def c2_job(args):
    """C2: one designed target at one ε: the path, the loop's kernels formed (against their nearest feature's amplitude), feature 0's
    entrant and runner-up at the start (and its nearest kernel), the runner-up's margin there, its share in the loop and the number of
    kernels formed for feature 0, the block against the path's."""
    name, eps = args[:2]
    rtol = args[2] if len(args) > 2 else 1e-9
    t0 = time.time()
    L = law(); K, c0 = full()
    cs, As = targets_C2()[name]
    f = K.synth(cs, As)
    R = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)))
    path = sorted(int(j) for j in np.flatnonzero(R['act']))
    t_e = R['events'][0][0]
    early = t_e * (1 + np.array([0.005, 0.02, 0.1]))                     # the runner-up just after feature 0's entry (output points only)
    taus = np.unique(np.r_[np.linspace(0, 1.3 * R['events'][-1][0], 27), early])
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus, rtol=rtol, atol=1e-12 if rtol >= 1e-10 else 1e-14)
    bk = formed_by_feature(c0, cs, As, W[-1], C[-1])
    near, big, fin = bk['near'], bk['big'], bk['fin']
    v = np.abs(K.corr(c0, f)[0])
    o = np.argsort(-np.where(near == 0, v, -1.0))
    e0, r0 = int(o[0]), int(o[1])
    w0 = np.abs(W[-1])[near == 0]
    # the runner-up's coefficient at the first output after feature 0's entry
    i1 = int(np.argmax(taus > R['events'][0][0]))
    dist = np.linalg.norm(C[:, r0, :] - cs[0], axis=1) / SIG                  # the runner-up's distance from feature 0 (referee round 53)
    return dict(name=name, eps=eps, path=path, d_max=float(dist[taus > t_e].max()), d_end=float(dist[-1]), events=[(float(t), int(j), k) for t, j, k, c in R['events']],
                formed=sorted(int(j) for j in np.flatnonzero(big)), e0=e0, r0=r0, nearest0=int(np.argmin(np.linalg.norm(c0 - cs[0], axis=1))),
                c_e0=c0[e0].tolist(), c_r0=c0[r0].tolist(), mu=float(1 - v[r0] / v[e0]), share_r=float(abs(W[-1][r0]) / w0.sum()),
                n0=bk['n_start'][0], n0_fin=bk['n_fin'][0], t_off=float(cs[0][1] - c0[e0][1]), captured=bk['captured'],
                fin=[(int(j), float(C[-1][j][0]), float(C[-1][j][1]), float(W[-1][j])) for j in np.flatnonzero(big)],
                w_tail=float(abs(W[i1][r0])), tau_tail=float(taus[i1]),
                dw=float(np.abs(W[-1][path] - R['w'][path]).max()), rms=float(np.sqrt(np.mean((f - K.synth(C[-1], W[-1])) ** 2))),
                n_fin=bk['n_fin'], rtol=rtol,
                w_early=[float(abs(W[int(np.argmin(np.abs(taus - te)))][r0])) for te in early], secs=time.time() - t0)


def part_C2():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c2_job, [(n, e) for n in ('near tie', 'window edge') for e in EPS_C], chunksize=1)


def v_residual(n, dd, A=A_B, Ap=AP_V):
    """V: the fit's residual on a fine trapezoid grid in ℝⁿ (n = 1, 2: the continuum to many digits), the feature A at the origin, an
    unfitted one A′ at dd·σ along e (the axis; in the plane the diagonal); p = (w_e, c_e, c_t) with the tail's coefficient w fixed:
    the residual vector (scaled by h^{n/2}) and e."""
    h = SIG / (20 if n == 1 else 8)
    e = np.ones(n) / np.sqrt(n)
    lo, hi = -6 * SIG, dd * SIG / np.sqrt(n) + 6 * SIG
    g = np.arange(lo, hi, h)
    x = g[:, None] if n == 1 else np.stack(np.meshgrid(g, g, indexing='ij'), -1).reshape(-1, 2)
    phi = lambda c: np.exp(-((x - c) ** 2).sum(1) / (2 * SIG ** 2))
    f = A * phi(np.zeros(n)) + Ap * phi(dd * SIG * e)

    def res(p, w):
        return (f - p[0] * phi(p[1:1 + n]) - w * phi(p[1 + n:])) * h ** (n / 2)
    return res, e


def v_job(args):
    """V: the tail's local minimum near its feature — the entrant (w_e, c_e) and the tail's centre free, its coefficient w fixed — for the
    coefficients W_V by continuation (the tail started at the feature), its displacement along e and across it (σ units); where it runs
    off (beyond 2σ) the fold, located by bisection in log w between the last bound and the first unbound coefficient (from the last
    bound state). If the tail is unbound at the grid's top (V3's large A′), the continuation starts nearer A (12 coefficients from 0.99A
    down to the grid); `top` the first coefficient, and no rows if the tail is unbound even there (the fold above it). The continuation
    stops if the entrant leaves its feature (|c_e| > 2σ: off V's branch, for A′ > A the entrant goes to A′). Each rest records its energy
    E*(w) (the minimum over the other parameters, insensitive to their error at second order): E* increasing in w along the branch means
    that with w free the coefficient decays (the hold of (v), referee round 49). A fourth argument, side = −1, starts the tail on the far
    side of its feature (away from A′: the mirror branch, referee round 53). Where the cold start at the grid's top is not a rest of the
    branch (at large A′ its small entrant is drawn to A′), fresh starts move the entrant along ∓e (referee round 53); `start` says which."""
    n, dd = args[:2]
    Ap = args[2] if len(args) > 2 else AP_V
    side = args[3] if len(args) > 3 else 1
    from scipy.optimize import least_squares
    res, e = v_residual(n, dd, Ap=Ap)

    def settle(p0, w):
        s_ = least_squares(res, p0, args=(w,), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
        return s_.x, float(np.linalg.norm(s_.x[1 + n:])) <= 2 * SIG
    ws = W_V if n == 1 else W_V[::4]
    t0_ = side * 1e-3 * SIG * e
    if not settle(np.r_[A_B - ws[0], np.zeros(n), t0_], ws[0])[1]:
        ws = np.r_[np.geomspace(0.99 * A_B, ws[0], 13)[:-1], ws]
    p, start = np.r_[A_B - ws[0], np.zeros(n), t0_], 'cold'
    q0, b0 = settle(p, ws[0])
    if not b0 or np.linalg.norm(q0[1:1 + n]) > 2 * SIG:
        for te in (-0.5, -1.0, -1.5, -2.0, 0.5, 1.0):
            q2, b2 = settle(np.r_[A_B - ws[0], te * SIG * e, t0_], ws[0])
            if b2 and np.linalg.norm(q2[1:1 + n]) <= 2 * SIG:
                p, start = q2, f'fresh, the entrant at {te:+g}σ along e'
                break
    rows, fold, along_fold, left = [], None, None, False
    for w in ws:
        q, bound = settle(p, w)
        if np.linalg.norm(q[1:1 + n]) > 2 * SIG:
            left = True
            break
        if not bound:
            if not rows:
                break
            lo, hi, pl = float(w), float(rows[-1]['w']), p.copy()
            for _ in range(14):
                m_ = float(np.sqrt(lo * hi))
                qm, bm = settle(pl, m_)
                if bm:
                    hi, pl = m_, qm
                else:
                    lo = m_
            fold = (lo, hi)
            along_fold = float(pl[1 + n:] @ e) / SIG
            break
        ct = q[1 + n:]
        rows.append(dict(w=float(w), we=float(q[0]), along=float(ct @ e) / SIG, across=float(np.linalg.norm(ct - (ct @ e) * e)) / SIG,
                         E=float(0.5 * np.sum(res(q, w) ** 2))))
        p = q
    return dict(n=n, d=dd, Ap=Ap, rows=rows, fold=fold, along_fold=along_fold, top=float(ws[0]), entrant_left=left, side=side, start=start)


def part_V():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(v_job, [(1, d) for d in D_V] + [(2, D_V2)], chunksize=1)


def part_VF():
    """VF: V's continuation on the far side (the tail started away from A′; referee round 53)."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(v_job, [(1, d, AP_V, -1) for d in D_V] + [(2, D_V2, AP_V, -1)], chunksize=1)


def v_law(dd, w, we, Ap=AP_V):
    """(v)'s displacement to leading order: Δ²/σ² = ℓ(d)·w_e/(w(w_e + w)), ℓ(d) = A′(d/σ)²e^{−d²/(4σ²)} — the binding (∝ M₂²) and the
    curvature pull (∝ M₂, M₂ = w(1 + w/w_e)|δ|² the block's second moment), the pull projected off the entrant's φ as the binding is (the
    entrant's coefficient slaved: ⟨r′, (∂²φ)_⊥⟩ = A′π^{n/2}σ^{n−2}(d/σ)²e^{−d²/(4σ²)}/4; referee round 49 — the unprojected pull gave
    (d/σ)² − 2)."""
    return float(np.sqrt(max(Ap * dd ** 2 * np.exp(-dd ** 2 / 4) * we / (w * (we + w)), 0.0)))


def toward_branches(V, V3=None):
    """V's and V3's one-dimensional branches toward A′, each (d, A′) once (V3's A′/A = AP_V/A_B column repeats V's: referee round 53)."""
    seen, out = set(), []
    for r in list(V) + (list(V3) if V3 is not None else []):
        k = (float(r['d']), float(r.get('Ap', AP_V)))
        if r['n'] == 1 and r.get('side', 1) == 1 and k not in seen:
            seen.add(k)
            out.append(r)
    return out


def e3_law(w, dd, side, A=A_B, Ap=AP_V):
    """(v)'s hold at the third moment (referee round 53), one dimension, the continuum: at the second moment the rests' energy
    E₀ − P²/(2B) does not depend on w; the third, M₃ = w(1 − w²/w_e²)δ³ (w_e = A − w), adds −(1/6)M₃P₃ with
    P₃ = ⟨r′, (∂³_cφ)_⊥⟩ = A′√π(d/σ)³e^{−d²/(4σ²)}/(8σ²) ((∂³_cφ) projected off φ and ∂_cφ), so at the rest δ = ±(M₂*/(w(1 + w/w_e)))^{1/2},
    M₂* = σ²ℓ(d): E*₃ = ∓(1/6)P₃(M₂*)^{3/2}(1 − w/w_e)(w(1 + w/w_e))^{−1/2} — increasing in w toward A′ (side +1), decreasing away (−1)."""
    l = Ap * dd ** 2 * np.exp(-dd ** 2 / 4)
    P3 = Ap * np.sqrt(np.pi) * dd ** 3 * np.exp(-dd ** 2 / 4) / (8 * SIG ** 2)
    we = A - w
    return -side * P3 * (SIG ** 2 * l) ** 1.5 * (1 - w / we) / np.sqrt(w * (1 + w / we)) / 6


def e3_ratios(r, ws=(0.28, 0.1, 0.037)):
    """A branch's secant dE*/dw over the third-moment law's (`e3_law`) at the steps nearest ws, and the displacement there (σ)."""
    out = []
    pairs = list(zip(r['rows'][:-1], r['rows'][1:]))
    for wq in ws:
        if not pairs:
            break
        a, b = min(pairs, key=lambda q: abs(np.log(0.5 * (q[0]['w'] + q[1]['w']) / wq)))
        if abs(np.log(0.5 * (a['w'] + b['w']) / wq)) > 0.2:
            continue
        side = r.get('side', 1)
        pred = (e3_law(b['w'], r['d'], side, Ap=r.get('Ap', AP_V)) - e3_law(a['w'], r['d'], side, Ap=r.get('Ap', AP_V))) / (b['w'] - a['w'])
        out.append((0.5 * (a['w'] + b['w']), ((b['E'] - a['E']) / (b['w'] - a['w'])) / pred, b['along']))
    return out


def c3_job(args):
    name, eps = args
    return c2_job((name, eps, RTOL_C3))


def c3_outcome(r):
    """C3's outcome from a C2 or C4 row: a capture tie (feature 2 two kernels by their final positions), captured, bound (feature 0 two
    kernels), one kernel (the fit exact) — or, with one kernel per feature and the residual above 10⁻⁸ at the end (the tail not settled),
    by the runner-up's recorded distance from feature 0 (referee round 53): not departed (it never left 2σ: still held, it can be
    recalled), stranded (it left and is at neither feature at the end), unsettled (no distance recorded). The bound loops' residual at
    the end (4–7·10⁻⁸) is the finite mass's relaxation, not a tail."""
    nf = r['n_fin']
    if nf[2] >= 2:
        return 'capture tie'
    if r['captured']:
        return 'captured'
    if nf[0] >= 2:
        return 'bound'
    if r['rms'] <= 1e-8:
        return 'one kernel'
    d = r.get('d_held', r.get('d_max'))
    return 'unsettled' if d is None else ('not departed' if d < 2.0 else 'stranded')


def part_C3():
    """C3: the capture scan — the near tie and its 5σ variant over EPS_C3 at the loop's tolerance RTOL_C3."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c3_job, [(n, e) for n in ('near tie', 'near tie 5σ') for e in EPS_C3], chunksize=1)


def fold_level(d2, A2, V):
    """The fold c(d)ℓ(d) of a same-sign feature A2 at the distance d2 (σ) from V's one-dimensional folds at A′ = AP_V, c(d) interpolated in
    d (None without V)."""
    vd = {r['d']: r for r in (V if V is not None else []) if r['n'] == 1 and r['fold'] is not None and r.get('Ap', AP_V) == AP_V}
    if not vd:
        return None
    ds_ = sorted(vd)
    c_ = np.interp(d2, ds_, [np.sqrt(vd[q]['fold'][0] * vd[q]['fold'][1]) / (AP_V * q ** 2 * np.exp(-q ** 2 / 4)) for q in ds_])
    return float(c_ * A2 * d2 ** 2 * np.exp(-d2 ** 2 / 4))


def c4_atol(eps):
    """C4's absolute tolerance: 1e-14, or 10⁻⁵ε below ε = 10⁻⁹ (the loop's state before the entry is at the masses' scale; referee round
    66: the absolute tolerance sets the step on every component, the dormant kernels' α included, and at ε = 10⁻¹¹ the step's
    error exceeded a formed runner-up's pin, which decides its side at 10⁻¹⁴σ; referee round 72: at ε = 10⁻¹¹)."""
    return min(1e-14, 1e-5 * eps)


def c4_job(args):
    """C4: C2's near tie with the runner-up's gap x held fixed (its margin μ = xε^{1/2}) at the mass ε: the loop (rtol RTOL_C3) with
    N_OUT_C4 outputs log-spaced from feature 0's entry τ_e to feature 2's own entry τ₂ on the path; the runner-up's distance from feature 0
    and its coefficient; its departure (the first output at which it is 2σ from feature 0) in s = (τ − τ_e)/ε^{1/2} against τ₂'s, its
    coefficient there and its largest after the entry; the kernels formed per feature at the end and the outcome (`c3_outcome`); a third
    argument sets the tolerance (C4b)."""
    x, eps = args[:2]
    rtol = args[2] if len(args) > 2 else RTOL_C3
    name = args[3] if len(args) > 3 else 'near tie'
    win = tuple(args[4]) if len(args) > 4 and args[4] is not None else None
    early = bool(args[5]) if len(args) > 5 else False
    t0 = time.time()
    L = law(); K, c0 = full()
    if win is not None:                                                # referee round 62: the samples' window (x₋, x₊, y₋, y₊)
        K = window(win[3], win[0], win[1], win[2])
    cs, As = targets_C2(x * np.sqrt(eps))[name]
    f = K.synth(cs, As)
    R = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)))
    near = np.argmin(np.linalg.norm(c0[:, None, :] - cs[None, :, :], axis=2), axis=1)
    ent = [e for e in R['events'] if e[2] == 'enter']
    t_e = ent[0][0]
    t2 = next(e[0] for e in ent if near[e[1]] == 2)
    v = np.abs(K.corr(c0, f)[0])
    o = np.argsort(-np.where(near == 0, v, -1.0))
    r0 = int(o[1])
    dense = t_e + (t2 - t_e) * np.geomspace(1e-7, 1.0, N_OUT_C4)
    if early:                       # referee round 62: from s = −100 about the path's entry, then from s = 0.1 (the loop's own layer)
        se = np.sqrt(eps)
        dense = t_e + se * np.r_[np.linspace(-100.0, 0.0, 51)[:-1], np.geomspace(0.1, (t2 - t_e) / se, N_OUT_C4)]
    taus = np.unique(np.r_[np.linspace(0, 1.3 * R['events'][-1][0], 27), dense])
    # the absolute tolerance below the masses (referee round 62's C4M at 10⁻¹⁷: at atol 1e-14 > ε the pre-entry state is unresolved
    # and the loop enters before s = −2000; referee round 66: 10⁻⁵ε wherever that is below 1e-14, the pin's side resolved)
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, taus, rtol=rtol, atol=c4_atol(eps))
    dist = np.linalg.norm(C[:, r0, :] - cs[0], axis=1) / SIG
    after = np.flatnonzero(taus > t_e)
    dep = next((int(i) for i in after if dist[i] > 2.0), None)
    held = after if dep is None else after[after < dep]
    bk = formed_by_feature(c0, cs, As, W[-1], C[-1])
    # the side (referee round 53): the runner-up's displacement along e (feature 0 → feature 2) where it first passes SIDE_R·σ, its
    # displacement along e at feature 2's entry, and the pair's closest approach after the entry
    ev = (cs[2] - cs[0]) / np.linalg.norm(cs[2] - cs[0])
    e0 = int(o[0])
    far = [int(i) for i in after if dist[i] > SIDE_R]
    side = None if not far else ('toward' if (C[far[0], r0] - cs[0]) @ ev > 0 else 'far')
    gap = np.linalg.norm(C[after, r0] - C[after, e0], axis=1) / SIG if after.size else np.array([np.nan])
    r = dict(name=name, side=side, closest=float(np.nanmin(gap)), x=x, eps=eps, mu=float(x * np.sqrt(eps)), tau_e=float(t_e), tau2=float(t2), s2=float((t2 - t_e) / np.sqrt(eps)),
             s_dep=float((taus[dep] - t_e) / np.sqrt(eps)) if dep is not None else None,
             w_dep=float(abs(W[dep - 1][r0])) if dep is not None else None,
             w_max=float(np.abs(W[held, r0]).max()) if held.size else 0.0, d_held=float(dist[held].max()) if held.size else 0.0,
             d2=float(np.linalg.norm(cs[2] - cs[0]) / SIG), A2=float(As[2]), n_fin=bk['n_fin'], captured=bk['captured'],
             rms=float(np.sqrt(np.mean((f - K.synth(C[-1], W[-1])) ** 2))), d_end=float(dist[-1]), w_end=float(abs(W[-1][r0])),
             rtol=rtol, atol=c4_atol(eps), win=None if win is None else list(win), early=early, secs=time.time() - t0)
    print(f"  C4 {name} x {x} eps {eps:g} rtol {rtol:g}: {c3_outcome(r)}, side {side}, s_dep {r['s_dep']}, s2 {r['s2']:.3g} "
          f"[{r['secs']:.0f}s]", flush=True)
    return r


def part_C4M():
    """C4M: the loop on the mirrored target (feature 2 across the face: the runner-up on its far side; referee round 53)."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c4_job, [(x, e, RTOL_C3, 'near tie mirrored', None, True) for x, e in C4M], chunksize=1)


def part_C4W():
    """C4W (referee round 62): C4's loop at (4.9, 10⁻¹⁰) on the window with the right edge 4.92σ from the block, on both targets."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c4_job, [(4.9, 1e-10, RTOL_C3, nm, WIN_RIGHT) for nm in ('near tie', 'near tie mirrored')], chunksize=1)


def x_R():
    """R's gaps: C3's (μ = MU_C2 at its masses), C4's and X_R_EXTRA."""
    return sorted({round(float(MU_C2 / np.sqrt(e)), 4) for e in EPS_C3} | set(X_C4) | set(X_R_EXTRA))


def race_job(args):
    """R: the ε-free race (referee round 53) — C2's near tie at μ = MU_R → 0 (or its mirror image, 'near tie mirrored') with only feature
    0's entrant (0, 2H) and runner-up (H, 2H) active and features 1 and 2 left in the residual: `inner_k`'s inner problem (the entrant from
    XI0_R behind its wall, the runner-up x behind it) on the full target, integrated in s = (τ − τ_e)/ε^{1/2} to s_end beyond the
    entrant's wall — C4 with s₂ = ∞. Records the pair's closest approach |c_r − c_e| after the wall (the inner layer's collapse), the side
    the runner-up leaves feature 0 to (its displacement along e, feature 0 → feature 2, where it first passes SIDE_R·σ after the wall — a
    formed runner-up's rest can be 0.15σ away; None if it never does), its displacement along e at the end,
    its departure s_d (2σ from feature 0), its coefficient after the layer (s = 10³), at the departure (the last output before) and at
    S_MARKS_R, its distance there, and the first s after the layer at which a coefficient below θ|A| has grown past it."""
    name, x = args[:2]
    rtol = args[2] if len(args) > 2 else RTOL_R
    s_end = args[3] if len(args) > 3 else S_END_R
    from scipy.integrate import solve_ivp
    t0 = time.time()
    L = law(); K, c0 = full()
    cs, As = targets_C2(MU_R)[name]
    f = K.synth(cs, As)
    e, r = kidx(c0, (0, 2 * H)), kidx(c0, (H, 2 * H))
    lam0 = RT.transit_scale(L, 1.0)
    z0 = np.r_[RT._alpha_of(L, -lam0 * (XI0_R + np.array([0.0, x])), RT.eta_crit), c0[[e, r]].ravel()]
    f0, j0 = RN.critical_rhs_nd(L, K, f), RN.critical_jac_nd(L, K, f)
    sw0 = XI0_R / abs(float(K.corr(c0[[e]], f)[0][0]))
    ss = np.unique(np.r_[np.linspace(0, 2 * sw0, 200), sw0 + np.geomspace(1e2, s_end, 400), sw0 + np.array(S_MARKS_R)])
    ss = ss[ss <= sw0 + s_end]
    sol = solve_ivp(lambda s_, z: lam0 * f0(s_, z), (0.0, ss[-1]), z0, method='LSODA', jac=lambda s_, z: lam0 * j0(s_, z), rtol=rtol,
                    atol=min(1e-13, rtol * 1e-3), t_eval=ss, dense_output=True)
    al, C = sol.y[:2].T, sol.y[2:].T.reshape(-1, 2, 2)
    w = RT.critical_w(L, al)
    xi_e = -RT.eta_crit(L, al[:, 0]) / lam0
    k = int(np.argmax(xi_e < 0))
    t_wall = float(np.interp(0.0, -xi_e[k - 1:k + 1], sol.t[k - 1:k + 1]))
    s = sol.t - t_wall

    def crossing(i, g):
        """The s at which g(z) changes sign between the outputs i − 1 and i, by bisection on the dense output (referee round 58: not
        the output grid's point)."""
        a_, b_ = sol.t[i - 1], sol.t[i]
        for _ in range(80):
            m_ = 0.5 * (a_ + b_)
            if g(sol.sol(m_)) > 0:
                b_ = m_
            else:
                a_ = m_
        return 0.5 * (a_ + b_) - t_wall
    ev = (cs[2] - cs[0]) / np.linalg.norm(cs[2] - cs[0])
    dist = np.linalg.norm(C[:, 1] - cs[0], axis=1) / SIG
    post = np.flatnonzero(s > 0)
    gap = np.linalg.norm(C[post, 1] - C[post, 0], axis=1) / SIG
    far = [int(i) for i in post if dist[i] > SIDE_R]
    side = None if not far else ('toward' if (C[far[0], 1] - cs[0]) @ ev > 0 else 'far')
    dep = next((int(i) for i in post if dist[i] > 2.0), None)
    s_dep = (crossing(dep, lambda z: np.linalg.norm(z[4:6] - cs[0]) / SIG - 2.0) if dep is not None and dep - 1 in post else
             (float(s[dep]) if dep is not None else None))
    at = lambda q: int(np.argmin(np.abs(s - q)))
    i3 = at(1e3)
    thr = TH_REL * abs(A_B)
    grow = (next((float(s[i]) for i in post if s[i] > 1e3 and (dep is None or i < dep) and abs(w[i, 1]) > thr), None)
            if abs(w[i3, 1]) < thr else None)                          # before the departure only (referee round 58)
    i_cl = int(post[np.argmin(gap)])
    dcl = C[i_cl, 1] - C[i_cl, 0]
    e_cl = float(dcl @ ev / np.linalg.norm(dcl)) if np.linalg.norm(dcl) > 0 else None
    marks = [q for q in S_MARKS_R if q <= s[-1] * (1 + 1e-9)]
    sw = next((int(i) for i in post if abs(w[i, 1]) > abs(w[i, 0])), None)            # the swap: the runner-up outgrows the entrant
    dist_e = np.linalg.norm(C[:, 0] - cs[0], axis=1) / SIG
    dep_e = next((int(i) for i in post if dist_e[i] > 2.0), None)                     # the entrant leaves feature 0
    wz = lambda z: np.abs(RT.critical_w(L, z[:2]))
    s_sw = crossing(sw, lambda z: wz(z)[1] - wz(z)[0]) if sw is not None and sw - 1 in post else None
    s_de = crossing(dep_e, lambda z: np.linalg.norm(z[2:4] - cs[0]) / SIG - 2.0) if dep_e is not None and dep_e - 1 in post else None
    near_e = [float(np.linalg.norm(C[-1, 0] - c_) / SIG) for c_ in cs]
    out = dict(name=name, x=float(x), rtol=rtol, s_end=float(s[-1]), status=int(sol.status), side=side, closest=float(gap.min()),
               e_closest=e_cl,
               s_swap=s_sw, s_de=s_de,
               e_dest=int(np.argmin(near_e)), e_dest_d=float(min(near_e)), w_e_end=float(abs(w[-1, 0])),
               s_d=s_dep, s_d_grid=float(s[dep]) if dep is not None else None, w_dep=float(abs(w[dep - 1, 1])) if dep is not None else None,
               d_dep=float(dist[dep]) if dep is not None else None, d_pre=float(dist[dep - 1]) if dep is not None else None,
               s_layer=float(s[i3]), w_layer=float(abs(w[i3, 1])), marks=marks, w_marks=[float(abs(w[at(q), 1])) for q in marks],
               d_marks=[float(dist[at(q)]) for q in marks], s_theta=grow, d_end=float(dist[-1]), w_end=float(abs(w[-1, 1])),
               e_end=float((C[-1, 1] - cs[0]) @ ev / SIG),
               secs=time.time() - t0)
    print(f"  R {name} x {x:g} rtol {rtol:g}: side {side}, closest {out['closest']:.1e}σ, s_d {out['s_d']}, w after the layer "
          f"{out['w_layer']:.4g}, at the end {out['w_end']:.4g} ({out['d_end']:.2f}σ) [{out['secs']:.0f}s]", flush=True)
    return out


def part_R():
    """R: the race without ε — the near tie at R's gaps, the mirrored target at X_RM (the far side), and the side's dependence on the
    tolerance at X_RS × RTOLS_RS (to S_END_RS)."""
    from multiprocessing import Pool
    jobs = ([('near tie', x) for x in x_R()] + [('near tie mirrored', x) for x in X_RM]
            + [('near tie', x, rt, S_END_RS) for x in X_RS for rt in RTOLS_RS])
    with Pool(N_WORKERS) as pool:
        out = pool.map(race_job, jobs, chunksize=1)
    n1, n2 = len(x_R()), len(X_RM)
    with Pool(N_WORKERS) as pool:
        long_ = pool.map(race_job, [('near tie mirrored', x, RTOL_R, se) for x, se in R_LONG], chunksize=1)
    return dict(main=out[:n1], mirror=out[n1:n1 + n2], scan=out[n1 + n2:], long=long_)


def window(top, xl=-1.0, xr=1.0, yb=-1.0):
    """Part C's samples with the window's top edge at `top` (`full`'s is 1) and its other edges at xl, xr, yb (referee round 62)."""
    return RN.BumpsGrid2(np.arange(xl, xr + 1e-9, DS), np.arange(yb, top + 1e-9, DS), SIG)


def edges(win):
    """The distances (σ) from the near tie's feature 0 to a window's edges (x₋, x₊, y₋, y₊) = win, in that order."""
    p = targets_C2(MU_R)['near tie'][0][0]
    return [(p[0] - win[0]) / SIG, (win[1] - p[0]) / SIG, (p[1] - win[2]) / SIG, (win[3] - p[1]) / SIG]


class PairSums:
    """R's two-kernel right-hand side (`RN.critical_rhs_nd` for the entrant and the runner-up, features 1 and 2 in the residual) with every
    sum over the samples — the target, the residual and its correlations — in the dtype dt: float64, or numpy's long double (the 80-bit
    format on x86-64, unit roundoff 5.4·10⁻²⁰); the state stays float64."""

    def __init__(self, L, K, cs, As, dt):
        self.L, self.dt = L, dt
        self.x1, self.x2 = K.x1.astype(dt), K.x2.astype(dt)
        self.s2, self.N = dt(SIG) ** 2, dt(K.N)
        self.F = sum(dt(a) * np.outer(self.g(self.x1, c[0]), self.g(self.x2, c[1])) for c, a in zip(cs, As))

    def g(self, x, c):
        return np.exp(-(x - self.dt(c)) ** 2 / (2 * self.s2))

    def __call__(self, s_, z):
        dt, L = self.dt, self.L
        w, M = RT.critical_w(L, z[:2]), RT.curve_mobility(L, z[:2], 0.0)
        C = z[2:].reshape(2, 2)
        G, Q = [self.g(self.x1, c[0]) for c in C], [self.g(self.x2, c[1]) for c in C]
        R = self.F - dt(w[0]) * np.outer(G[0], Q[0]) - dt(w[1]) * np.outer(G[1], Q[1])
        out = np.empty(6)
        for i in range(2):
            gp, qp = G[i] * (self.x1 - dt(C[i, 0])) / self.s2, Q[i] * (self.x2 - dt(C[i, 1])) / self.s2
            out[i] = M[i] * float(G[i] @ R @ Q[i] / self.N)
            out[2 + 2 * i:4 + 2 * i] = (dt(L.G * L.m * w[i]) * np.array([gp @ R @ Q[i], G[i] @ R @ qp]) / self.N).astype(float)
        return out

    def moments(self, cb, W):
        """At one centre c̄ with the coefficient W (r = f* − Wφ(c̄)): F = ⟨r, ∇_cφ⟩, ⟨r, ∇²_cφ⟩, ⟨∇_cφ∇_cφᵀ⟩ and ⟨φ, ∇_cφ⟩, in dt."""
        dt, s2, N = self.dt, self.s2, self.N
        u, v = self.x1 - dt(cb[0]), self.x2 - dt(cb[1])
        g, q = np.exp(-u ** 2 / (2 * s2)), np.exp(-v ** 2 / (2 * s2))
        gp, qp, gpp, qpp = g * u / s2, q * v / s2, g * (u ** 2 / s2 - 1) / s2, q * (v ** 2 / s2 - 1) / s2
        R = self.F - dt(W) * np.outer(g, q)
        F = np.array([gp @ R @ q, g @ R @ qp]) / N
        H = np.array([[gpp @ R @ q, gp @ R @ qp], [gp @ R @ qp, g @ R @ qpp]]) / N
        Q2 = np.array([[(gp @ gp) * (q @ q), (gp @ g) * (q @ qp)], [(gp @ g) * (q @ qp), (g @ g) * (qp @ qp)]]) / N
        pd = np.array([(g @ gp) * (q @ q), (g @ g) * (q @ qp)]) / N
        return F, H, Q2, pd

    def block(self, cb, W):
        """At the block's centre c̄ with the coefficient W: ⟨r, ∇²_cφ(c̄)⟩ (2 × 2, r = f* − Wφ(c̄)) and the rest's shift with W,
        ∂c̄*/∂W = (∂F/∂c̄)⁻¹⟨φ, ∇_cφ⟩ (F = ⟨f* − Wφ(c̄), ∇_cφ(c̄)⟩; ∂F/∂c̄ = ⟨r, ∇²_cφ⟩ − W⟨∇_cφ∇_cφᵀ⟩)."""
        _, H, Q2, pd = self.moments(cb, W)
        return H.astype(float), np.linalg.solve((H - self.dt(W) * Q2).astype(float), pd.astype(float))

    def force(self, cb, W):
        """RC (referee round 66): the merged block's force F = ⟨f* − Wφ(c̄), ∇_cφ(c̄)⟩ and its Jacobian ∂F/∂c̄."""
        F, H, Q2, _ = self.moments(cb, W)
        return F.astype(float), (H - self.dt(W) * Q2).astype(float)


def sample_coords(lo, hi, ds, ideal):
    """RD (referee round 72): one axis of a window's samples as Decimals — the committed grid's float64 coordinates
    (`np.arange(lo, hi + 1e-9, ds)`, each converted exactly; np.arange rounds lo + k·fl(ds) at every point) or the ideal grid's exact
    lo + k·ds, as many points."""
    from decimal import Decimal as Dm
    xs = np.arange(lo, hi + 1e-9, ds)
    return [Dm(float(v)) for v in xs] if not ideal else [Dm(repr(float(lo))) + k * Dm(repr(float(ds))) for k in range(xs.size)]


def exact_drift(cb, W, cs, As, x1, x2, prec=RD_PREC):
    """RD (referee round 72): the merged block's rest's shift with W, ∂c̄*/∂W = (∂F/∂c̄)⁻¹⟨φ, ∇_cφ⟩ (`PairSums.block`), every sum over the
    samples x1 × x2 (Decimals) exact to prec digits (Python's decimal; the two-dimensional sums factor into one-dimensional ones); the
    state (cb, W), the target and σ the float64 values, converted exactly. Returns (∂c̄*/∂W, ∂F/∂c̄ as (11, 12, 22), ⟨φ²⟩), floats."""
    from decimal import Decimal as Dm, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        s2 = Dm(float(SIG)) ** 2

        def fac(xs, c):
            u = [x_ - Dm(float(c)) for x_ in xs]
            g = [(-(d * d) / (2 * s2)).exp() for d in u]
            return g, [g_ * d / s2 for g_, d in zip(g, u)], [g_ * (d * d / s2 - 1) / s2 for g_, d in zip(g, u)]
        dot = lambda a, b: sum((p_ * q_ for p_, q_ in zip(a, b)), Dm(0))
        g, gp, gpp = fac(x1, cb[0])
        q, qp, qpp = fac(x2, cb[1])
        fe = [(Dm(float(a)), fac(x1, c[0])[0], fac(x2, c[1])[0]) for c, a in zip(cs, As)]
        N, Wd = Dm(len(x1) * len(x2)), Dm(float(W))
        H11 = sum(a * dot(gf, gpp) * dot(qf, q) for a, gf, qf in fe) - Wd * dot(g, gpp) * dot(q, q)
        H12 = sum(a * dot(gf, gp) * dot(qf, qp) for a, gf, qf in fe) - Wd * dot(g, gp) * dot(q, qp)
        H22 = sum(a * dot(gf, g) * dot(qf, qpp) for a, gf, qf in fe) - Wd * dot(g, g) * dot(q, qpp)
        J11 = (H11 - Wd * dot(gp, gp) * dot(q, q)) / N
        J12 = (H12 - Wd * dot(gp, g) * dot(q, qp)) / N
        J22 = (H22 - Wd * dot(g, g) * dot(qp, qp)) / N
        p1, p2 = dot(g, gp) * dot(q, q) / N, dot(g, g) * dot(q, qp) / N
        det = J11 * J22 - J12 * J12
        return ((float((J22 * p1 - J12 * p2) / det), float((J11 * p2 - J12 * p1) / det)), (float(J11), float(J12), float(J22)),
                float(dot(g, g) * dot(q, q) / N))


def alias_drift(cb, J, phi2, lo, ds):
    """RD (referee round 72): the aliasing's shift on an ideal grid of spacing ds from lo = (x₋, y₋) — Poisson summation's first harmonic,
    ⟨φ, ∂_iφ⟩ = −⟨φ²⟩(2πq/ds)sin(2π(c_i − lo_i)/ds) up to O(q²), q = e^{−π²σ²/ds²} ((ii)(c)'s), through ∂F/∂c̄ = J (11, 12, 22)."""
    q = np.exp(-np.pi ** 2 * SIG ** 2 / ds ** 2)
    pd = [-phi2 * 2 * np.pi * q / ds * np.sin(2 * np.pi * (cb[i] - lo[i]) / ds) for i in range(2)]
    return np.linalg.solve(np.array([[J[0], J[1]], [J[1], J[2]]]), pd)


def relative_law(P, s, C, w, W, i0, iend, D0, lam0, Gm, slaved=True):
    """RP: the relative offset D = c_r − c_e (e the entrant, r the runner-up) integrated from the output i0 with D0 over the run's outputs
    to iend, the coefficients frozen per output interval (exponential steps; stopped where D stops being finite). With the block's
    centroid c_b = (w_e c_e + w c_r)/W slaved — the midpoint c̄ = c_b + ρD, ρ = (w_e − w)/(2W) — the exact first order
        dD/ds = (2(w − w_e)/W)dc̄/ds + Λ₀(2Gm ww_e/W)⟨r, ∇²_cφ(c̄)⟩D
    becomes dD/ds = κ[Λ₀(2Gm ww_e/W)⟨r, ∇²_cφ(c̄)⟩D + (2(w − w_e)/W)(ρ′D + dc_b/ds)], κ = W²/(2(w² + w_e²)) (referee round 62: the rate is
    κ times the coefficient; slaved=False: κ = 1, ρ′ = 0, round 58's law), dc_b/ds the block's rest moving with W, (∂c̄*/∂W)dW/ds
    (`PairSums.block`): the pin, about the run's own pair (referee round 66: below the floor the relative mode is continued about a
    merged block instead, RC). Returns (the outputs' indices, D at each (σ), ∂c̄*/∂W at each (σ per unit W))."""
    dWds = np.gradient(W, s)
    rho = (w[:, 0] - w[:, 1]) / (2 * W)
    drho = np.gradient(rho, s)
    Dl, idx, Ds, dc = np.asarray(D0, float).copy(), [i0], [np.asarray(D0, float) / SIG], [None]
    for i in range(i0, iend):
        Hb, dcdW = P.block(0.5 * (C[i, 0] + C[i, 1]), W[i])
        A = lam0 * 2 * Gm * w[i, 0] * w[i, 1] / W[i] * Hb
        k_ = 2 * (w[i, 1] - w[i, 0]) / W[i]
        g = k_ * dcdW * dWds[i]
        if slaved:
            kap = W[i] ** 2 / (2 * (w[i, 0] ** 2 + w[i, 1] ** 2))
            A, g = kap * (A + k_ * drho[i] * np.eye(2)), kap * (g + 0.0)
        lv, V = np.linalg.eig(A)
        Vi = np.linalg.inv(V)
        h_ = s[i + 1] - s[i]
        E1 = (V @ np.diag(np.exp(lv * h_)) @ Vi).real
        P1 = (V @ np.diag(np.where(np.abs(lv * h_) > 1e-12, np.expm1(lv * h_) / np.where(lv == 0, 1, lv), h_)) @ Vi).real
        Dn = E1 @ Dl + P1 @ g
        if not np.all(np.isfinite(Dn)):                                 # referee round 66: no NaN carried on
            break
        Dl = Dn
        idx.append(i + 1); Ds.append(Dl / SIG); dc.append(dcdW / SIG)
    return np.array(idx), np.array(Ds), dc


def rp_job(args):
    """RP (referee rounds 58 and 62): R's race (C2's near tie at μ = MU_R, the runner-up x behind the entrant) on Part C's window with
    the top edge at `top`, or on a window (x₋, x₊, y₋, y₊) (WIN_FAR: every edge ≥ 11.7σ from the block, the continuum; WIN_RIGHT: the
    right edge 4.92σ from it), the sums in float64 (prec 64) or the long double (80) (`PairSums`): the relative offset D = c_r − c_e at
    S_MARKS_RP against `relative_law` from D(S0_RP) — the slaved-centroid law, and round 58's for comparison —; the Jacobian's
    eigenvalues nearest the relative mode's at S_JAC_RP against the law's rates; the least |D| after the wall before the departure and its
    share along e; the side (SIDE_R), the departure s_d (2σ, bisected on the dense output) and the climb's end (|D| first back at
    D_CLIMB·σ after its least; referee round 66: s_d is the climb's end plus the hold after it); the growth rate of |D| on
    s ∈ [3·10³, 10⁴] (`growth_rate`); the samples' count N. Below the floating-point floor the relative mode is RC's."""
    name, x, top, prec = args[:4]
    win = tuple(args[4]) if len(args) > 4 and args[4] is not None else None
    from scipy.integrate import solve_ivp
    t0 = time.time()
    L = law(); c0 = full()[1]
    K = window(top) if win is None else window(win[3], win[0], win[1], win[2])
    cs, As = targets_C2(MU_R)[name]
    P = PairSums(L, K, cs, As, np.longdouble if prec == 80 else np.float64)
    e, r = kidx(c0, (0, 2 * H)), kidx(c0, (H, 2 * H))
    lam0 = RT.transit_scale(L, 1.0)
    f = K.synth(cs, As)
    j0 = RN.critical_jac_nd(L, K, f)
    z0 = np.r_[RT._alpha_of(L, -lam0 * (XI0_R + np.array([0.0, x])), RT.eta_crit), c0[[e, r]].ravel()]
    sw0 = XI0_R / abs(float(K.corr(c0[[e]], f)[0][0]))
    ss = np.unique(np.r_[np.linspace(0, 2 * sw0, 200), sw0 + np.geomspace(1e2, S_END_RP, 1200)])
    sol = solve_ivp(lambda s_, z: lam0 * P(s_, z), (0.0, ss[-1]), z0, method='LSODA', jac=lambda s_, z: lam0 * j0(s_, z), rtol=RTOL_R,
                    atol=min(1e-13, RTOL_R * 1e-3), t_eval=ss, dense_output=True)
    al, C = sol.y[:2].T, sol.y[2:].T.reshape(-1, 2, 2)
    w = RT.critical_w(L, al)
    xi_e = -RT.eta_crit(L, al[:, 0]) / lam0
    k = int(np.argmax(xi_e < 0))
    t_wall = float(np.interp(0.0, -xi_e[k - 1:k + 1], sol.t[k - 1:k + 1]))
    s = sol.t - t_wall
    D, W = C[:, 1] - C[:, 0], w.sum(1)
    ev = (cs[2] - cs[0]) / np.linalg.norm(cs[2] - cs[0])
    dist = np.linalg.norm(C[:, 1] - cs[0], axis=1) / SIG
    post = np.flatnonzero(s > 0)
    dep = next((int(i) for i in post if dist[i] > 2.0), None)
    far = [int(i) for i in post if dist[i] > SIDE_R]
    side = None if not far else ('toward' if (C[far[0], 1] - cs[0]) @ ev > 0 else 'far')
    s_d = None
    if dep is not None and dep - 1 in post:
        a_, b_ = sol.t[dep - 1], sol.t[dep]
        for _ in range(80):
            m_ = 0.5 * (a_ + b_)
            zm = sol.sol(m_)
            if np.linalg.norm(zm[4:6] - cs[0]) / SIG > 2.0:
                b_ = m_
            else:
                a_ = m_
        s_d = 0.5 * (a_ + b_) - t_wall
    held = post if dep is None else post[post < dep]
    i_min = int(held[np.argmin(np.linalg.norm(D[held], axis=1))])
    Gm = L.G * L.m
    i0 = int(np.argmin(np.abs(s - S0_RP)))
    i1 = int(np.argmin(np.abs(s - max(S_MARKS_RP)))) + 1
    at = lambda q: int(np.argmin(np.abs(s - q)))

    def marks(slaved):
        idx, Ds, dc = relative_law(P, s, C, w, W, i0, i1, D[i0], lam0, Gm, slaved)
        return [Ds[int(np.argmin(np.abs(s[idx] - q)))].tolist() for q in S_MARKS_RP], dc[int(np.argmin(np.abs(s[idx] - S_MARKS_RP[1])))]
    law_, dcdW_at = marks(True)
    law0 = marks(False)[0]
    def climb_end(th):                                                 # the climb's end, bisected on the dense output (referee round 66)
        up = [int(i) for i in held if i > i_min and np.linalg.norm(D[i]) >= th * SIG]
        if not up:
            return None
        a_, b_ = sol.t[up[0] - 1], sol.t[up[0]]
        for _ in range(80):
            m_ = 0.5 * (a_ + b_)
            zm = sol.sol(m_)
            if np.linalg.norm(zm[4:6] - zm[2:4]) >= th * SIG:
                b_ = m_
            else:
                a_ = m_
        return 0.5 * (a_ + b_) - t_wall
    s_climb = climb_end(D_CLIMB)
    s_climbs = [climb_end(th) for th in THETAS_CLIMB]                   # referee round 72: s_d − s_climb(θ) against θ
    jac = []
    for q in S_JAC_RP:
        i = at(q)
        if not (0 < s[i] and (dep is None or i < dep)):
            continue
        Hb, _ = P.block(0.5 * (C[i, 0] + C[i, 1]), W[i])
        lamD = lam0 * 2 * Gm * w[i, 0] * w[i, 1] / W[i]
        kap = W[i] ** 2 / (2 * (w[i, 0] ** 2 + w[i, 1] ** 2))
        pr = np.sort(np.linalg.eigvals(lamD * Hb).real)
        ej = np.linalg.eigvals(lam0 * j0(0.0, sol.y[:, i]))
        near = sorted(float(ej[q].real) for q in np.argsort(np.min(np.abs(ej[:, None] - kap * pr[None, :]), axis=1))[:2])
        # referee round 72: the centroid's offset y unslaved — the (D, y) system on the merged block's moments, Q = ⟨r, ∇²_cφ⟩,
        # R = Q − W⟨∇_cφ∇_cφᵀ⟩: dD/ds = Λ₀Gm[(2ww_e/W)QD + (w − w_e)Ry], dy/ds = Λ₀Gm[((w² + w_e²)/W)Ry + (ww_e(w − w_e)/W²)QD];
        # slaving y gives the κ-law, and its first correction is λ ≈ κλ_D(1 − ε_b), ε_b = (w − w_e)²ww_e/(w² + w_e²)²·R⁻¹Q
        _, Hm, Q2m, _ = P.moments(0.5 * (C[i, 0] + C[i, 1]), W[i])
        we_, wr_, Rm = w[i, 0], w[i, 1], (Hm - P.dt(W[i]) * Q2m).astype(float)
        Hm = Hm.astype(float)
        cD, cy = we_ * wr_ * (wr_ - we_) / W[i] ** 2, (we_ ** 2 + wr_ ** 2) / W[i]
        M4 = lam0 * Gm * np.block([[2 * we_ * wr_ / W[i] * Hm, (wr_ - we_) * Rm], [cD * Hm, cy * Rm]])
        e4 = np.linalg.eigvals(M4)
        slow = sorted(float(e4[q].real) for q in np.argsort(np.min(np.abs(e4[:, None] - kap * pr[None, :]), axis=1))[:2])
        eb = np.sort(np.linalg.eigvals((wr_ - we_) * cD / cy ** 2 * np.linalg.solve(Rm, Hm)).real)
        jac.append(dict(s=float(s[i]), jac=near, law=(kap * pr).tolist(), law0=pr.tolist(), kappa=float(kap), law_y=slow,
                        eps_b=eb.tolist()))
    grow = growth_rate(s, D, post)
    out = dict(name=name, x=float(x), top=float(top), win=None if win is None else list(win), N=int(K.N), prec=int(prec),
               status=int(sol.status), nfev=int(sol.nfev), D=[(D[at(q)] / SIG).tolist() for q in S_MARKS_RP], law=law_, law0=law0,
               dcdW=None if dcdW_at is None else np.asarray(dcdW_at).tolist(), s_climb=s_climb, s_climbs=s_climbs, jac=jac,
               least=float(np.linalg.norm(D[i_min]) / SIG), s_least=float(s[i_min]),
               e_least=float(D[i_min] @ ev / np.linalg.norm(D[i_min])) if np.linalg.norm(D[i_min]) > 0 else None, side=side, s_d=s_d,
               grow=grow, w_layer=float(abs(w[at(1e3), 1])), ulp=float(np.spacing(abs(C[at(S_MARKS_RP[1]), 1, 1])) / SIG),
               secs=time.time() - t0)
    print(f"  RP {name} x {x:g} top {top:g} win {win} prec {prec}: D(700) {out['D'][1]}, law {out['law'][1]}, side {side}, s_d {s_d} "
          f"[{out['secs']:.0f}s]", flush=True)
    return out


def growth_rate(s, D, idx, lo=3e3, hi=1e4):
    """RP: |D|'s growth rate, the slope of log|D| against s fitted over the outputs idx with s in [lo, hi] (None below three outputs)."""
    gi = [i for i in idx if lo <= s[i] <= hi and np.linalg.norm(D[i]) > 0]
    return float(np.polyfit(s[gi], np.log(np.linalg.norm(D[gi], axis=1)), 1)[0]) if len(gi) > 2 else None


def merged_rhs(P, lam0):
    """RC (referee round 66): the merged block's flow (α_e, α_r, c_m) — both kernels at one centre c_m: the coefficients' laws as in the
    race, the centre at the pair's centroid velocity Λ₀Gm((w_e² + w²)/W)⟨r, ∇_cφ(c_m)⟩ (each kernel's velocity ∝ its coefficient)."""
    def rhs(s_, y):
        out = P(s_, np.r_[y[:4], y[2:4]])
        w = RT.critical_w(P.L, y[:2])
        return lam0 * np.r_[out[:2], (w[0] * out[2:4] + w[1] * out[4:6]) / w.sum()]
    return rhs


def merged_block(name, x):
    """RC, RD (referee rounds 66, 72): R's race on WIN_FAR (`rp_job`'s run) to s₁ = S1_RC, and there the merged block — both kernels at
    the centroid c_b = (w_e c_e + w c_r)/W less its offset y = (∂F/∂c_b)⁻¹F(c_b) off the merged rest (`PairSums.force`) — run on to
    S_CONT_RC (`merged_rhs`): the setting, the race's state at s₁ (w, W, D = c_r − c_e, y) and the merged run; s in WIN_FAR's time from
    the entrant's wall at t_wall."""
    from scipy.integrate import solve_ivp
    L = law(); c0 = full()[1]
    K = window(WIN_FAR[3], WIN_FAR[0], WIN_FAR[1], WIN_FAR[2])
    cs, As = targets_C2(MU_R)[name]
    P = PairSums(L, K, cs, As, np.float64)
    e, r = kidx(c0, (0, 2 * H)), kidx(c0, (H, 2 * H))
    lam0 = RT.transit_scale(L, 1.0)
    f = K.synth(cs, As)
    j0 = RN.critical_jac_nd(L, K, f)
    z0 = np.r_[RT._alpha_of(L, -lam0 * (XI0_R + np.array([0.0, x])), RT.eta_crit), c0[[e, r]].ravel()]
    sw0 = XI0_R / abs(float(K.corr(c0[[e]], f)[0][0]))
    tol = dict(method='LSODA', rtol=RTOL_R, atol=min(1e-13, RTOL_R * 1e-3), dense_output=True)
    race = lambda s_, z: lam0 * P(s_, z)
    jac = lambda s_, z: lam0 * j0(s_, z)
    ss = np.linspace(0, 2 * sw0, 200)                                   # rp_job's outputs to the wall: the same t_wall
    sol = solve_ivp(race, (0.0, 2 * sw0 + 2 * S1_RC), z0, jac=jac, t_eval=ss, **tol)
    xi_e = -RT.eta_crit(L, sol.y[0]) / lam0
    k = int(np.argmax(xi_e < 0))
    t_wall = float(np.interp(0.0, -xi_e[k - 1:k + 1], sol.t[k - 1:k + 1]))
    z1 = sol.sol(t_wall + S1_RC)
    w1 = RT.critical_w(L, z1[:2])
    C1 = z1[2:].reshape(2, 2)
    W1 = w1.sum()
    cb, D1 = (w1[0] * C1[0] + w1[1] * C1[1]) / W1, C1[1] - C1[0]
    F1, dF1 = P.force(cb, W1)
    y = np.linalg.solve(dF1, F1)
    cm = cb - y
    msol = solve_ivp(merged_rhs(P, lam0), (t_wall + S1_RC, t_wall + S_CONT_RC), np.r_[z1[:2], cm], **tol)
    assert msol.status == 0, msol.message
    return dict(L=L, K=K, cs=cs, As=As, P=P, lam0=lam0, j0=j0, tol=tol, race=race, jac=jac, t_wall=t_wall, w1=w1, W1=W1, D1=D1, y=y,
                msol=msol)


def rc_job(args):
    """RC (referee round 66): the relative mode in the continuum, continued about a merged block (`merged_block`): the centres'
    perturbation δc — δc_e = y − (w/W)D, δc_r = y + (w_e/W)D, D the race's c_r − c_e at s₁ (ymode 'zero': y = 0) — by the race's exact
    Jacobian's centre block about the merged block, δc′ = J_cc δc (`critical_jac_nd`; the coefficients' coupling vanishes on a merged
    block at rest), in midpoint matrix exponentials, `sub` per output interval (N_CONT_RC outputs, geometric to S_CONT_RC), the norm
    carried in logarithms: no floor. Past its least, once the relative offset D = δc_r − δc_e has regrown to inj·σ and to twice its least,
    δc is injected into the float64 race — the merged block's coefficients, its centre plus δc — run on (to S_END_RC) to the departure
    (2σ). Returns the least |D| (log₁₀ σ), its share along e and its s, the contraction's e-folds from s₁ to it, the injection's s, the
    side (SIDE_R), the climb's end (|D| = D_CLIMB·σ) and s_d — s in WIN_FAR's time. Referee round 72: (1) the rest's drift as a forcing —
    the merged base is no solution of the race: its kernels' speeds go as their coefficients, so a block moving at dc_b/ds forces D by
    κk·dc_b/ds (κ = W²/(2(w² + w_e²)), k = 2(w − w_e)/W; δc_e by −(w/W) of it, δc_r by +(w_e/W)), here dc_b/ds = (∂c̄*/∂W)dW/ds for
    the drift of the float64 sums (`PairSums.block`), of exact sums on the committed samples (their float64 coordinates) and on the ideal
    grid (the aliasing; `exact_drift` at s₁, checked at S_CONT_RC): the forced D at the least and at the injection, |D| (log₁₀ σ) and its
    share along e; (2) the climb's end at every THETAS_CLIMB; (3) the runner-up's coefficient after the injection — its change at each
    climb's end, its peak (with |D| and s there) and where it has fallen 1% below its value at the injection."""
    name, x = args[:2]
    sub = int(args[2]) if len(args) > 2 else SUB_RC
    ymode = args[3] if len(args) > 3 else 'y'
    inj = float(args[4]) if len(args) > 4 else INJ_RC
    from scipy.integrate import solve_ivp
    from scipy.linalg import expm
    t0 = time.time()
    m = merged_block(name, x)
    L, K, cs, As, P, lam0, j0, tol, race, jac, t_wall, msol = (m[k_] for k_ in ('L', 'K', 'cs', 'As', 'P', 'lam0', 'j0', 'tol', 'race',
                                                                                 'jac', 't_wall', 'msol'))
    w1, W1, D1, y = m['w1'], m['W1'], m['D1'], m['y']
    yv = y if ymode == 'y' else np.zeros(2)
    dc = np.r_[yv - w1[1] / W1 * D1, yv + w1[0] / W1 * D1]
    ev = (cs[2] - cs[0]) / np.linalg.norm(cs[2] - cs[0])
    lnn, u = np.log(np.linalg.norm(dc)), dc / np.linalg.norm(dc)
    lg0 = np.log10(np.linalg.norm(D1) / SIG)
    least, s_inj = (np.inf, None, None), None
    xc = [sample_coords(WIN_FAR[0], WIN_FAR[1], DS, False), sample_coords(WIN_FAR[2], WIN_FAR[3], DS, False)]
    xi = [sample_coords(WIN_FAR[0], WIN_FAR[1], DS, True), sample_coords(WIN_FAR[2], WIN_FAR[3], DS, True)]

    def exact_pair(s_):
        z_ = msol.sol(t_wall + s_)
        c_, W_ = z_[2:4], float(RT.critical_w(L, z_[:2]).sum())
        return np.array(exact_drift(c_, W_, cs, As, *xc)[0]), np.array(exact_drift(c_, W_, cs, As, *xi)[0])
    d_coord, d_alias = exact_pair(S1_RC)
    d_coord1, d_alias1 = exact_pair(S_CONT_RC)
    Wat = lambda s_: float(RT.critical_w(L, msol.sol(t_wall + s_)[:2]).sum())
    Fc, forced = np.zeros((4, 3)), {}

    def snap():
        Du = Fc[2:] - Fc[:2]
        out_ = {}
        for k_, v_ in (('coord', Du[:, :2] @ d_coord), ('alias', Du[:, :2] @ d_alias), ('f64', Du[:, 2])):
            nv = float(np.linalg.norm(v_))
            out_[k_] = (float(np.log10(nv / SIG)), float(v_ @ ev / nv)) if nv > 0 else (None, None)
        return out_
    grid = np.geomspace(S1_RC, S_CONT_RC, N_CONT_RC)
    for a_, b_ in zip(grid[:-1], grid[1:]):
        hs = (b_ - a_) / sub
        for j in range(sub):
            zm = msol.sol(t_wall + a_ + (j + 0.5) * hs)
            Jcc = lam0 * j0(0.0, np.r_[zm[:4], zm[2:4]])[2:, 2:]
            u = expm(Jcc * hs) @ u
            nu = np.linalg.norm(u)
            lnn, u = lnn + np.log(nu), u / nu
            wm = RT.critical_w(L, zm[:2])                               # the forcing, constant over the substep (referee round 72)
            Wm = wm.sum()
            fk = Wm ** 2 / (2 * (wm[0] ** 2 + wm[1] ** 2)) * 2 * (wm[1] - wm[0]) / Wm * (Wat(a_ + (j + 1) * hs) - Wat(a_ + j * hs)) / hs
            Bm = np.column_stack([np.r_[-wm[1] / Wm * fk * d_, wm[0] / Wm * fk * d_]
                                  for d_ in (np.array([1.0, 0.0]), np.array([0.0, 1.0]), P.block(zm[2:4], Wm)[1])])
            Aug = np.zeros((7, 7))
            Aug[:4, :4], Aug[:4, 4:] = Jcc * hs, Bm * hs
            Ea = expm(Aug)
            Fc = Ea[:4, :4] @ Fc + Ea[:4, 4:]
        Dv = u[2:] - u[:2]
        lg = (lnn + np.log(np.linalg.norm(Dv)) - np.log(SIG)) / np.log(10)
        if lg < least[0]:
            least = (float(lg), float(Dv @ ev / np.linalg.norm(Dv)), float(b_))
            forced['least'] = snap()
        if lg >= np.log10(inj) and lg > least[0] + np.log10(2.0):
            s_inj = float(b_)
            forced['inj'] = snap()
            break
    out = dict(name=name, x=float(x), sub=sub, ymode=ymode, inj=inj, N=int(K.N), lg_D1=float(lg0), least=least[0], e_least=least[1],
               s_least=least[2], efolds=float((lg0 - least[0]) * np.log(10)), y=float(np.linalg.norm(y) / SIG), s_inj=s_inj,
               side=None, s_climb=None, s_d=None, forced=forced, d_coord=(d_coord / SIG).tolist(), d_alias=(d_alias / SIG).tolist(),
               d_change=float(max(np.linalg.norm(d_coord1 - d_coord) / np.linalg.norm(d_coord),
                                  np.linalg.norm(d_alias1 - d_alias) / np.linalg.norm(d_alias))),
               s_climbs=[None] * len(THETAS_CLIMB), w_climbs=[None] * len(THETAS_CLIMB), w_peak=None, s_wdecay=None)
    if s_inj is not None:
        zm = msol.sol(t_wall + s_inj)
        dcv = np.exp(lnn) * u
        zi = np.r_[zm[:2], zm[2:4] + dcv[:2], zm[2:4] + dcv[2:]]

        def ev_dep(t_, z):
            return np.linalg.norm(z[4:6] - cs[0]) / SIG - 2.0
        ev_dep.terminal, ev_dep.direction = True, 1

        def ev_side(t_, z):
            return np.linalg.norm(z[4:6] - cs[0]) / SIG - SIDE_R
        ev_side.direction = 1

        def ev_climb(t_, z):
            return np.linalg.norm(z[4:6] - z[2:4]) / SIG - D_CLIMB
        ev_climb.direction = 1

        def ev_far(t_, z):                                              # on the far side the race ends at its rests: stop there
            return (z[4:6] - cs[0]) @ ev / SIG + SIDE_R
        ev_far.terminal, ev_far.direction = True, -1

        def climb_at(th):                                               # referee round 72: the climb's end at each THETAS_CLIMB
            def ev_(t_, z):
                return np.linalg.norm(z[4:6] - z[2:4]) / SIG - th
            ev_.direction = 1
            return ev_
        sol2 = solve_ivp(race, (t_wall + s_inj, t_wall + S_END_RC), zi, jac=jac,
                         events=(ev_dep, ev_side, ev_climb, ev_far) + tuple(climb_at(th) for th in THETAS_CLIMB), **tol)
        if sol2.t_events[1].size:
            out['side'] = 'toward' if (sol2.y_events[1][0][4:6] - cs[0]) @ ev > 0 else 'far'
        out['s_climb'] = float(sol2.t_events[2][0] - t_wall) if sol2.t_events[2].size else None
        out['s_d'] = float(sol2.t_events[0][0] - t_wall) if sol2.t_events[0].size else None
        out['status'] = int(sol2.status)
        out['s_climbs'] = [float(sol2.t_events[4 + i][0] - t_wall) if sol2.t_events[4 + i].size else None
                           for i in range(len(THETAS_CLIMB))]
        # the runner-up's coefficient along the injected race (referee round 72): its change from the injection, at each climb's end,
        # at its peak (with |D| there) and where it has fallen 1% below its value at the injection
        ts = t_wall + s_inj + np.geomspace(1.0, max(sol2.t[-1] - t_wall - s_inj, 2.0), 4000)
        Z = sol2.sol(ts)
        rel = RT.critical_w(L, Z[:2].T)[:, 1] / RT.critical_w(L, zi[:2])[1] - 1
        Dn = np.linalg.norm(Z[4:6] - Z[2:4], axis=0) / SIG
        out['w_climbs'] = [None if sc is None else float(np.interp(sc + t_wall, ts, rel)) for sc in out['s_climbs']]
        ip = int(np.argmax(rel))
        out['w_peak'] = (float(rel[ip]), float(ts[ip] - t_wall), float(Dn[ip]))
        dn = np.flatnonzero((np.arange(rel.size) > ip) & (rel < -0.01))
        out['s_wdecay'] = float(ts[dn[0]] - t_wall) if dn.size else None
    out['secs'] = time.time() - t0
    print(f"  RC {name} x {x:g} sub {sub} {ymode} inj {inj:g}: least 10^{least[0]:.2f} ({least[1]:+.2f}) at {least[2]:.4g}, "
          f"side {out['side']}, climb {out['s_climb']}, s_d {out['s_d']} [{out['secs']:.0f}s]", flush=True)
    return out


def rd_job(args):
    """RD (referee round 72): what moves the block's rest with W on samples without an edge within reach — RC's merged block
    (`merged_block`: the target `name`, the gap x, on WIN_FAR) at S_RD: ∂c̄*/∂W with the sums in float64 and in the long double
    (`PairSums.block`) and exact (`exact_drift`, RD_PREC digits) on the committed samples (np.arange's float64 coordinates), exact on
    the ideal grid (lo + k·ds exactly) against the aliasing's first harmonic (`alias_drift`), and all four on a dyadic grid (DS_DYADIC,
    its float64 coordinates exact) over the same window; the committed and the dyadic coordinates' largest errors."""
    name, x = args
    t0 = time.time()
    m = merged_block(name, x)
    L, cs, As = m['L'], m['cs'], m['As']
    lo = (WIN_FAR[0], WIN_FAR[2])
    P80 = PairSums(L, m['K'], cs, As, np.longdouble)
    Kd = RN.BumpsGrid2(np.arange(WIN_FAR[0], WIN_FAR[1] + 1e-9, DS_DYADIC), np.arange(WIN_FAR[2], WIN_FAR[3] + 1e-9, DS_DYADIC), SIG)
    Pd64, Pd80 = PairSums(L, Kd, cs, As, np.float64), PairSums(L, Kd, cs, As, np.longdouble)
    ax = lambda ds, ideal: [sample_coords(WIN_FAR[0], WIN_FAR[1], ds, ideal), sample_coords(WIN_FAR[2], WIN_FAR[3], ds, ideal)]
    xc, xi, xdc, xd = ax(DS, False), ax(DS, True), ax(DS_DYADIC, False), ax(DS_DYADIC, True)
    v = lambda d: (np.asarray(d, float) / SIG).tolist()
    rows = []
    for s_ in S_RD:
        z = m['msol'].sol(m['t_wall'] + s_)
        W, c = float(RT.critical_w(L, z[:2]).sum()), z[2:4]
        ex = exact_drift(c, W, cs, As, *xc)[0]
        ei, Ji, p2i = exact_drift(c, W, cs, As, *xi)
        ed, Jd, p2d = exact_drift(c, W, cs, As, *xd)
        rows.append(dict(s=float(s_), W=W, c=c.tolist(), f64=v(m['P'].block(c, W)[1]), ld=v(P80.block(c, W)[1]), exact=v(ex),
                         ideal=v(ei), alias=v(alias_drift(c, Ji, p2i, lo, DS)), dy_f64=v(Pd64.block(c, W)[1]),
                         dy_ld=v(Pd80.block(c, W)[1]), dy_exact=v(ed), dy_alias=v(alias_drift(c, Jd, p2d, lo, DS_DYADIC))))
    err = lambda a, b: [max(abs(float(p_ - q_)) for p_, q_ in zip(a[i], b[i])) for i in range(2)]
    out = dict(name=name, x=float(x), rows=rows, coord_err=err(xc, xi), dy_coord_err=err(xdc, xd), N=int(m['K'].N), N_dy=int(Kd.N),
               secs=time.time() - t0)
    print(f"  RD {name} x {x:g}: exact {rows[1]['exact']}, ideal {rows[1]['ideal']} [{out['secs']:.0f}s]", flush=True)
    return out


def part_RD():
    """RD: RD_JOBS."""
    return [rd_job(a) for a in RD_JOBS]


def part_RC():
    """RC: every gap on both targets, and RC_CHECKS (the substeps, y, the injection's level) at x = 4.9."""
    from multiprocessing import Pool
    jobs = [(nm, x) for nm in ('near tie', 'near tie mirrored') for x in X_RC] + list(RC_CHECKS)
    with Pool(N_WORKERS) as pool:
        return pool.map(rc_job, jobs, chunksize=1)


def part_RP():
    """RP: the relative mode's pin — the formed gaps and a tail on Part C's window in float64, x = 4.9 in the long double and on the
    mirrored target, the formed gaps 3, 4.9, 6.04 with the window's top edge nearer (0.95) and farther (1.5, also in the long double);
    referee round 62: every gap on both targets on WIN_FAR (the continuum), 3 and 4.9 on WIN_RIGHT."""
    from multiprocessing import Pool
    jobs = ([('near tie', x, 1.0, 64) for x in X_RP] + [('near tie', 4.9, 1.0, 80), ('near tie mirrored', 4.9, 1.0, 64)]
            + [('near tie', x, top, 64) for top in TOPS_RP if top != 1.0 for x in (3.0, 4.9, 6.04)] + [('near tie', 4.9, 1.5, 80)]
            + [(nm, x, 1.0, 64, WIN_FAR) for nm in ('near tie', 'near tie mirrored') for x in X_RP]
            + [(nm, x, 1.0, 64, WIN_RIGHT) for nm in ('near tie', 'near tie mirrored') for x in (3.0, 4.9)])
    with Pool(N_WORKERS) as pool:
        return pool.map(rp_job, jobs, chunksize=1)


def fold_job(args):
    """FC: the static fold in C4's own geometry (referee round 53) — Part C's samples, C2's near tie at μ = MU_R with features 1 (−0.6,
    opposite sign) and 2 (0.7) unfitted, or feature 1 removed ('alone') —: the entrant (w_e, c_e) and the runner-up's centre free, its
    coefficient w held, the rest continued in w along W_FC from the side toward feature 2 (side +1) or away (−1); the fold — the rest
    beyond 2σ from feature 0 — by bisection. Returns the fold (None: none down to W_FC's end), the rest's distance and its displacement
    along e at the fold (or at the end)."""
    name, side = args
    from scipy.optimize import least_squares
    K, c0 = full()
    cs, As = targets_C2(MU_R)['near tie']
    if name == 'alone':
        cs, As = cs[[0, 2]], As[[0, 2]]
    f = K.synth(cs, As)
    ev = (cs[-1] - cs[0]) / np.linalg.norm(cs[-1] - cs[0])

    def res(p, w):
        return (f - p[0] * K.phi(p[1:3][None, :])[:, 0] - w * K.phi(p[3:5][None, :])[:, 0]) / np.sqrt(K.N)

    settle = lambda p0, w: least_squares(res, p0, args=(w,), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000).x
    inside = lambda q: np.linalg.norm(q[3:5] - cs[0]) / SIG <= 2.0 and np.linalg.norm(q[1:3] - cs[0]) / SIG <= 2.0
    where = lambda q: (float(np.linalg.norm(q[3:5] - cs[0]) / SIG), float((q[3:5] - cs[0]) @ ev / SIG))
    p, last, fold, at_ = np.r_[A_B - W_FC[0], cs[0], cs[0] + side * 0.05 * SIG * ev], None, None, None
    for w in W_FC:
        q = settle(p, w)
        if not inside(q):
            if last is not None:
                lo, hi, pl = float(w), float(last[0]), last[1]
                for _ in range(30):
                    m_ = float(np.sqrt(lo * hi))
                    qm = settle(pl, m_)
                    if inside(qm):
                        hi, pl = m_, qm
                    else:
                        lo = m_
                fold, at_ = (lo, hi), where(pl)
            break
        last, p = (float(w), q), q
    return dict(name=name, side=side, fold=fold, at=at_ if fold is not None else (where(last[1]) if last is not None else None),
                w_end=float(W_FC[-1]))


def part_FC():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(fold_job, [('near tie', 1), ('near tie', -1), ('alone', 1)], chunksize=1)


def part_C4():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c4_job, [(x, e) for x in X_C4 for e in EPS_C4], chunksize=1)


def part_C4b():
    """C4b: C4's points at the tolerance RTOL_C4B (an outcome decided near a saddle can depend on the integration's error)."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c4_job, [(x, e, RTOL_C4B) for x in X_C4 for e in EPS_C4], chunksize=1)


def b4x_job(x):
    """B4x: one of B4's unsettled problems continued to MAX_DOUBLINGS_B4X doublings (the drifts per doubling, the stop, the shares)."""
    x1, x2 = x
    L = law(); K, c0 = local()
    xd = margins_diag(K, c0, x1, x2)
    pts, T = TIES['vertex']
    f = A_B * K.phi(np.atleast_2d(T))[:, 0]
    t0 = time.time()
    out = RN.inner_k(L, K, f, np.array(pts, float), np.array([0.0, x1, x2, xd]), share_tol=SHARE_TOL, max_doublings=MAX_DOUBLINGS_B4X,
                     max_secs=8 * INNER_SECS)
    return dict(x=(x1, x2), shares=out['shares'], n=int((out['shares'] > TH_REL).sum()), stop=out['stop'], drifts=list(out['drifts']),
                secs=time.time() - t0)


def part_B4x(B4):
    from multiprocessing import Pool
    pts = [tuple(r['x']) for r in B4['rows'] if not r['stationary']]
    with Pool(N_WORKERS) as pool:
        return dict(rows=pool.map(b4x_job, pts, chunksize=4))


def part_V3():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(v_job, [(1, d, a) for d in D_V3 for a in AP_V3], chunksize=1)


def d2_share(x, eps, t):
    """D2: the edge loop's runner-up share at the gap x, the mass ε, the feature at the offset t along the face (Part A's loop)."""
    L = law(); K, c0 = local()
    cstar, idx = place(K, c0, A_B, 'edge', (x,), eps, t)
    mu, v_e = margins(K, c0, A_B, cstar, idx)
    f = A_B * K.phi(np.atleast_2d(cstar))[:, 0]
    T_, W, C, Lam, AL = RN.loop_nd(L, K, f, c0, eps, np.linspace(0, TAU_A / v_e, 23))
    w = np.abs(W[-1][idx])
    return float(w[1] / w.sum())


def d2_job(args):
    """D2: the loop's formation threshold in the gap (the runner-up share crossing TH_REL, by bisection in [4.5, 8]) at (ε, t)."""
    eps, t = args
    t0 = time.time()
    return dict(eps=eps, t=t, x=crossing(lambda x: d2_share(x, eps, t), 4.5, 8.0, 12), secs=time.time() - t0)


def part_D2():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(d2_job, [(e, t) for e in EPS_D2 for t in T_D2], chunksize=1)


def local_3d():
    """T3: the 8 kernels at the corners of the cube {0, H}³ and a tensor grid of samples (spacing DS) about it in ℝ³."""
    g = DS * np.arange(-9, 12)
    xs = np.stack(np.meshgrid(g, g, g, indexing='ij'), -1).reshape(-1, 3)
    corners = np.array([[a, b, c] for a in (0, H) for b in (0, H) for c in (0, H)], float)
    return RN.BumpsND(xs, SIG), corners


def t3_job(x):
    """T3: the 8-kernel inner problem at a vertex of hℤ³ (the feature A_B at the cube's centre, the entrant the corner 0, the corner
    1_S behind it by Σ_{i∈S}x_i — the margins' product rule to leading order) against the product of the three edges' 2-kernel shares at
    the vertex's transverse offset (n − 1)^{1/2}H/2 = H/√2."""
    L = law(); K, corners = local_3d()
    f = A_B * K.phi(np.full((1, 3), H / 2))[:, 0]
    S = [tuple(int(v) for v in np.round(c / H)) for c in corners]
    gaps = [sum(xi for xi, si in zip(x, s_) if si) for s_ in S]
    o = np.argsort(gaps, kind='stable')
    t0 = time.time()
    r = RN.inner_k(L, K, f, corners[o], np.array(gaps)[o], share_tol=SHARE_TOL, max_secs=4 * INNER_SECS)
    sh = np.empty(8); sh[o] = r['shares']
    return dict(x=x, S=S, shares=sh.tolist(), stop=r['stop'], drift=r['drift'], secs=time.time() - t0)


def t3_edge(x):
    """T3: the edge's 2-kernel inner problem's runner-up share at the gap x with the feature at the transverse offset H/√2."""
    return float(inner_job(('edge', A_B, (0.0, x), H / np.sqrt(2)))['shares'][1])


def part_T3():
    from multiprocessing import Pool
    xs = sorted({v for x in T3_X for v in x})
    with Pool(N_WORKERS) as pool:
        edge = dict(zip(xs, pool.map(t3_edge, xs, chunksize=1)))
        cube = pool.map(t3_job, list(T3_X), chunksize=1)
    return dict(edge=edge, cube=cube)


def slope(x, y):
    return float(np.polyfit(np.log10(x), np.log10(y), 1)[0])


def areas(N, grid):
    """The areas of {N ≥ 3} and {N = 4} in the gaps' quarter-plane from counts on a symmetric grid (trapezoid weights: half at 0)."""
    g = np.asarray(grid, float)
    d = g[1] - g[0]
    w = np.where(g == 0, d / 2, d)
    W2 = np.outer(w, w)
    return float((W2 * (N >= 3)).sum()), float((W2 * (N == 4)).sum())


def vertex_counts(B4):
    """B4's counts on its grid of gaps, symmetric (the problems are posed at x₁₀ ≤ x₀₁; the vertex's count is symmetric in its gaps)."""
    g = np.asarray(B4['grid'], float)
    idx = {v: i for i, v in enumerate(np.round(g, 9))}
    N = np.zeros((g.size, g.size), int)
    for r in B4['rows']:
        i, j = idx[round(r['x'][0], 9)], idx[round(r['x'][1], 9)]
        N[i, j] = N[j, i] = r['n']
    return g, N


def p_layers(xbar, eps):
    """P(N ≥ 2) for a feature uniform in a cell of hℤ² at leading order: the four faces' layers of width w = 2σ²x̄_θε^{1/2}/H less
    their overlaps at the corners, 4w/H − 4w²/H²."""
    w = 2 * SIG ** 2 / H ** 2 * H * xbar * np.sqrt(eps) / H
    return 4 * w - 4 * w * w


def e_extra(xbar_n, n):
    """The coefficient of ε^{1/2} in E[N − 1] on hℤⁿ: 2n faces, each a layer of width (2σ²/H)x̄_θ⁽ⁿ⁾ε^{1/2} over the (n − 1)-cube."""
    return 2 * n * (2 * SIG ** 2 / H ** 2) * xbar_n


def p_vertices(a):
    """The coefficient of ε in P(N ≥ 3) from B4's area a of {N ≥ 3} (or in P(N = 4) from the area of {N = 4}): four vertices per
    cell, each vertex's sector mapped to the gaps by (2σ²ε^{1/2}/H)·(x₁₀, x₀₁)."""
    k2 = 2 * SIG ** 2 / H ** 2
    return 4 * k2 * k2 * a


def b4u_note(B4, B4u, B4x=None):
    """B4's caption: why the unsettled problems stopped (rerun, B4u), whether their counts changed, and their continuation to
    MAX_DOUBLINGS_B4X doublings (B4x: the counts and the shares' change after the cap — a margin over the last drift bounds nothing while
    the drifts still grow: referee round 49)."""
    old = {tuple(r['x']): r for r in B4['rows']}
    by = {}
    for r in B4u['rows']:
        by.setdefault(r['stop'], []).append(r)
    same = sum(int(r['n'] == old[tuple(r['x'])]['n']) for r in B4u['rows'])
    parts = []
    for k, rs in sorted(by.items()):
        lab = dict(cap='at the doubling cap', clock='at the clock', converged='settled on rerun').get(k, 'with no reason recorded')
        if k == 'converged':
            parts.append(f"{len(rs)} {lab}")
            continue
        d = [r['drift'] for r in rs]
        parts.append(f"{len(rs)} {lab} (the last drift {min(d):.1e}" + (f" … {max(d):.1e}" if max(d) > min(d) else '') + ")")
    chg = max(float(np.abs(np.asarray(r['shares']) - np.asarray(old[tuple(r['x'])]['shares'])).max()) for r in B4u['rows'])
    note = ''
    if B4x is not None:
        u = {tuple(r['x']): r for r in B4u['rows']}
        xs_ = B4x['rows']
        same_x = sum(int(r['n'] == old[tuple(r['x'])]['n']) for r in xs_)
        mv = max(float(np.abs(np.asarray(r['shares']) - np.asarray(u[tuple(r['x'])]['shares'])).max()) for r in xs_ if tuple(r['x']) in u)
        pk = [int(np.argmax(r['drifts'])) + 1 for r in xs_ if r['drifts']]
        stops = {k: sum(r['stop'] == k for r in xs_) for k in ('converged', 'cap', 'clock')}
        note = (f"; continued to {MAX_DOUBLINGS_B4X} doublings (B4x): the counts unchanged at {same_x} of {len(xs_)}, the shares moved "
                f"after the first cap by at most {mv:.1e}, the drifts peaking at doublings {min(pk)}–{max(pk)}, "
                f"{stops['converged']} settled, {stops['cap']} at the new cap, {stops['clock']} at the clock")
    return (f"; the {len(B4u['rows'])} others rerun with the stop's reason recorded: " + ', '.join(parts) +
            f", their counts unchanged at {same} of {len(B4u['rows'])}, their shares moved by at most {chg:.1e}" + note)


def t3_product(pv, S):
    """T3's product rule: per corner 1_S the product of the edges' runner-up shares p_i (i ∈ S) and 1 − p_i (i ∉ S)."""
    return np.array([np.prod([pv[i] if s_[i] else 1 - pv[i] for i in range(len(pv))]) for s_ in S])


def fc_fold(FC):
    """The configuration's own fold toward feature 2 (FC's 'near tie', side +1), its geometric mean; None without it."""
    r = next((q for q in (FC if FC is not None else []) if q['name'] == 'near tie' and q['side'] == 1 and q['fold'] is not None), None)
    return None if r is None else float(np.sqrt(r['fold'][0] * r['fold'][1]))


def vf_lines(VF):
    """VF's report: the far branches' rests and the secant of their energy (referee round 53)."""
    if VF is None:
        return []
    out = ["", "VF. The far side (the tail started away from A′): the rest's displacement along e (σ; negative: away) at w = 0.3, 0.03, "
           "0.003, 0.001, the fold, and the secant dE*/dw along the branch:", "",
           "| n | d/σ | displacement at w = 0.3, 0.03, 0.003, 0.001 | fold | dE*/dw: steps < 0 of all, range |", "|---|---|---|---|---|"]
    for r in VF:
        cells = []
        for wq in (0.3, 0.03, 0.003, 0.001):
            rr = min(r['rows'], key=lambda q: abs(np.log(q['w'] / wq))) if r['rows'] else None
            cells.append(f"{rr['along']:+.4f}" if rr is not None and abs(np.log(rr['w'] / wq)) < 0.1 else '—')
        gr = [(b['E'] - a['E']) / (b['w'] - a['w']) for a, b in zip(r['rows'][:-1], r['rows'][1:])]
        fold = (f"{np.sqrt(r['fold'][0] * r['fold'][1]):.4g}" if r['fold'] is not None else
                f"none down to {r['rows'][-1]['w']:g}" if r['rows'] else 'no rest')
        out.append(f"| {r['n']} | {r['d']:g} | {', '.join(cells)} | {fold} | "
                   + (f"{sum(g < 0 for g in gr)} of {len(gr)}, {min(gr):.2e} … {max(gr):.2e}" if gr else '—') + " |")
    return out


def side_of(r, sep=''):
    """A C4 row's recorded side (referee round 53; rows from before it carry none)."""
    return f"{r['side']}{sep}" if 'side' in r else ''


def c4m_lines(C4M):
    """C4M's report: the loop on the mirrored target (referee round 53); referee round 62: the outputs from s = −100 about the path's
    entry (s = (τ − τ_e)/ε^{1/2}) — before, the grid began at s = 10⁻⁷(τ₂ − τ_e)/ε^{1/2}, past the inner layer at the smallest ε —
    and the far rule's upper boundary (ε = 6.5·10⁻¹⁷, 3.7·10⁻¹⁷)."""
    if C4M is None:
        return []
    out = ["", "C4M. The loop on the mirrored target (features 1 and 2 reflected across the face: feature 2 on the entrant's side, the "
           "runner-up displaced away from it), rtol as C4's, the absolute tolerance min(10⁻¹⁴, 10⁻⁵ε) (`c4_atol`: "
           "at 10⁻¹⁴ the state before the entry, at the masses' scale, is unresolved at ε = 10⁻¹⁷ and the loop enters thousands of s "
           "early), the "
           "outputs from s = −100 about the path's entry (s = (τ − τ_e)/ε^{1/2}; a row marked † has the grid of C4, from s = 10⁻⁷s₂, "
           "which can begin past the inner layer):", "",
           "| x | ε | side | the pair's closest approach (σ) | max w held | departure s | s₂ | kernels per feature (final) | the runner-up "
           "at the end: w, distance from feature 0 (σ) | outcome |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(C4M, key=lambda r: (r['x'], -r['eps'])):
        out.append(f"| {r['x']:g} | {r['eps']:g}{'' if r.get('early') else '†'} | {r['side']} | {r['closest']:.1e} | {r['w_max']:.4g} | "
                   + (f"{r['s_dep']:.3g}" if r['s_dep'] is not None else 'none') + f" | {r['s2']:.3g} | {', '.join(map(str, r['n_fin']))} | "
                   f"{r['w_end']:.4g}, {r['d_end']:.3f} | {c4m_outcome(r)} |")
    return out


def c4w_lines(C4W):
    """C4W's report (referee round 62): C4's loop on the window with the right edge 4.92σ from the block."""
    if C4W is None or not len(C4W):
        return []
    out = ["", f"C4W. C4's loop at (x, ε) = (4.9, 10⁻¹⁰) on the window x ∈ [{WIN_RIGHT[0]:g}, {WIN_RIGHT[1]:g}], y ∈ "
           f"[{WIN_RIGHT[2]:g}, {WIN_RIGHT[3]:g}] (the right edge {edges(WIN_RIGHT)[1]:.2f}σ from feature 0, the top {edges(WIN_RIGHT)[3]:.2f}σ): "
           "the side the formed runner-up "
           "takes, its closest approach, largest held w, departure, the kernels per feature and the outcome, on both targets:", "",
           "| target | x | ε | side | the pair's closest approach (σ) | max w held | departure s | s₂ | kernels per feature (final) | "
           "outcome |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in C4W:
        out.append(f"| {r['name']} | {r['x']:g} | {r['eps']:g} | {r['side']} | {r['closest']:.2e} | {r['w_max']:.4g} | "
                   + (f"{r['s_dep']:.3g}" if r['s_dep'] is not None else 'none') + f" | {r['s2']:.3g} | {', '.join(map(str, r['n_fin']))} | "
                   f"{c4m_outcome(r) if 'mirrored' in r['name'] else c3_outcome(r)} |")
    return out


def c4m_outcome(r):
    """C4M's outcome: `c3_outcome`, an unformed runner-up that never left 2σ labelled by its side whatever the residual (referee round 58:
    'not departed' and 'one kernel' told the far side's held tails apart only by the residual's cut)."""
    o = c3_outcome(r)
    if o in ('not departed', 'one kernel') and r.get('d_held', 2.0) < 2.0 and r.get('w_end', 1.0) < TH_REL * A_B:
        return f"held, unformed ({r['side']})"
    return o


def rp_lines(RP):
    """RP's report: the relative offset against the slaved-centroid law (and round 58's), per window and arithmetic; the Jacobian's rates
    against the law's at every row (referee round 66); s_d and the climb's end in the samples' time (referee round 62)."""
    if RP is None or not len(RP):
        return []
    N0 = window(1.0).N
    edge = lambda top: (top - 2 * H) / SIG
    wname = lambda r: (f"top {r['top']:g} ({edge(r['top']):.2f}σ)" if r.get('win') is None else
                       ('far edges' if list(r['win']) == list(WIN_FAR) else 'right edge' if list(r['win']) == list(WIN_RIGHT) else
                        str(r['win'])))
    out = ["## RP. The relative mode: the pin (referee rounds 58 and 62)", "",
           "R's race (the near tie at μ → 0, the runner-up x behind the entrant) with the sums over the samples in float64 or in the 80-bit "
           "long double (unit roundoff 5.4·10⁻²⁰; the state float64), on Part C's window with its top edge at 1 (the block "
           f"{edge(1.0):.2f}σ below it), nearer or farther, on the window with every edge at least {min(edges(WIN_FAR)):.2f}σ from "
           f"feature 0 (x ∈ "
           f"[{WIN_FAR[0]:g}, {WIN_FAR[1]:g}], y ∈ [{WIN_FAR[2]:g}, {WIN_FAR[3]:g}]: no edge within reach, the block's rest moving with W "
           "only by the samples' own float64 coordinates — the same in the long double and in exact sums — and by the aliasing, RD) and "
           f"with the right edge {edges(WIN_RIGHT)[1]:.2f}σ from it (x ∈ [{WIN_RIGHT[0]:g}, {WIN_RIGHT[1]:g}], y ∈ [{WIN_RIGHT[2]:g}, "
           f"{WIN_RIGHT[3]:g}]); N the samples (⟨·⟩ averages over them, so a window of N samples runs at N₀/N the speed of Part C's, "
           f"N₀ = {N0}). The relative offset D = c_r − c_e (σ) at s = " + ', '.join(f"{q:g}" for q in S_MARKS_RP)
           + f" against the law integrated from s = {S0_RP:g} with the block's centroid c_b = (w_e c_e + w c_r)/W slaved to its rest — "
           "dD/ds = κ[Λ₀(2Gm ww_e/W)⟨r, ∇²_cφ(c̄)⟩D + (2(w − w_e)/W)(ρ′D + (∂c̄*/∂W)dW/ds)], κ = W²/(2(w² + w_e²)), ρ = (w_e − w)/(2W) "
           "(the midpoint c̄ = c_b + ρD), ∂c̄*/∂W = (∂F/∂c̄)⁻¹⟨φ, ∇_cφ⟩, F = ⟨f* − Wφ(c̄), ∇_cφ(c̄)⟩ —; ∂c̄*/∂W (σ per unit W) at s = 700; "
           "the least |D| after the wall before the departure and its share along e; the side, the climb's end (|D| back at "
           f"{D_CLIMB:g}σ) and the departure s_d (in brackets times N₀/N: in Part C's time) and |D|'s growth rate on s ∈ [3·10³, 10⁴]:", "",
           "| target | x | window | N | sums | D at s = " + ', '.join(f"{q:g}" for q in S_MARKS_RP) + " (x, y) | the law | ∂c̄*/∂W (x, y) | "
           "least \\|D\\| (its share along e) | side | climb's end (×N₀/N) | s_d (×N₀/N) | growth rate |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    vec = lambda v: '—' if v is None else f"({v[0]:+.2e}, {v[1]:+.2e})"
    key = lambda r: (r['name'], r.get('win') is not None, str(r.get('win')), r['top'], r['prec'], r['x'])
    for r in sorted(RP, key=key):
        sd = 'none' if r['s_d'] is None else f"{r['s_d']:.4g} ({r['s_d'] * N0 / r.get('N', N0):.4g})"
        sc = '—' if r.get('s_climb') is None else f"{r['s_climb']:.4g} ({r['s_climb'] * N0 / r.get('N', N0):.4g})"
        out.append(f"| {r['name']} | {r['x']:g} | {wname(r)} | {r.get('N', N0)} | {'long double' if r['prec'] == 80 else 'float64'} | "
                   + '; '.join(vec(v) for v in r['D']) + " | " + '; '.join(vec(v) for v in r['law']) + f" | {vec(r['dcdW'])} | "
                   + f"{r['least']:.2e} (" + ('—' if r['e_least'] is None else f"{r['e_least']:+.2f}") + f") | {r['side']} | {sc} | {sd} | "
                   + ('—' if r['grow'] is None else f"{r['grow']:.3g}") + " |")
    std = sorted([r for r in RP if r['top'] == 1.0 and r['prec'] == 64 and r['name'] == 'near tie' and r.get('win') is None],
                 key=lambda r: r['x'])
    rat = lambda r, k, i: (f"{r[k][i][1] / r['D'][i][1]:.3f}" if r.get(k) and r[k][i] is not None and r['D'][i][1] != 0 else '—')
    if std:
        out += ["", "On Part C's window, the law's D_y over the run's at s = " + ', '.join(f"{q:g}" for q in S_MARKS_RP) + " — the slaved "
                "centroid's law (round 58's, κ = 1 and ρ′ = 0, in brackets): " + '; '.join(
                    f"x = {r['x']:g}: " + ', '.join(f"{rat(r, 'law', i)} ({rat(r, 'law0', i)})" for i in range(len(S_MARKS_RP)))
                    for r in std) + "."]
        if all('ulp' in r for r in std):
            out += ["The run's D is the difference of float64 centres, so D_y is quantised in steps of the runner-up's c_y's last place, "
                    + ', '.join(f"{r['ulp']:.3g}σ" for r in std) + f" (x = {', '.join(f'{r[chr(120)]:g}' for r in std)}): the pin "
                    "at x = 3 is " + ', '.join(f"{abs(r['D'][1][1]) / r['ulp']:.0f}" for r in std if r['x'] == 3.0) + " such steps, the "
                    "law there within " + ', '.join(f"{abs(r['law'][1][1] - r['D'][1][1]) / r['ulp']:.1f}" for r in std if r['x'] == 3.0)
                    + " of them at s = 700."]
    jr = [r for r in RP if r.get('jac') and r['prec'] == 64]
    if jr:
        js = [(r, j) for r in sorted(jr, key=key) for j in r['jac']]
        err = [max(abs(a - b) for a, b in zip(sorted(j['law']), sorted(j['jac']))) / max(abs(b) for b in j['jac']) for r, j in js]
        erry = [max(abs(a - b) for a, b in zip(sorted(j['law_y']), sorted(j['jac']))) / max(abs(b) for b in j['jac']) if 'law_y' in j
                else None for r, j in js]
        out += ["", "The rate (referee rounds 62 and 66), at every float64 row: the Jacobian's eigenvalues nearest the relative mode's (the "
                "full two-kernel system, Λ₀ times `critical_jac_nd`) against the slaved law's κλ_D·eig⟨r, ∇²_cφ⟩ — its error over the largest rate — "
                "and round 58's λ_D·eig⟨r, ∇²_cφ⟩ (λ_D = Λ₀·2Gm ww_e/W): " + '; '.join(
                    f"{r['name']}, x = {r['x']:g}, {wname(r)}, s = {j['s']:.4g}: " + ', '.join(f"{a:+.4g}" for a in j['jac'])
                    + " against " + ', '.join(f"{a:+.4g}" for a in j['law']) + f" ({e_:.1%}; κ = {j['kappa']:.3f}"
                    + (f"; ε_b {', '.join(f'{b_:.3f}' for b_ in j['eps_b'])}; the (D, y) system " + ', '.join(f"{a:+.4g}" for a in j['law_y'])
                       + f" ({ey:.1e})" if 'law_y' in j else '') + "; round 58: "
                    + ', '.join(f"{a:+.4g}" for a in j['law0']) + ")" for (r, j), e_, ey in zip(js, err, erry)) + ".",
                f"Over these rows κ = {min(j['kappa'] for r, j in js):.3f}–{max(j['kappa'] for r, j in js):.3f} and the slaved law's error "
                f"{min(err):.1%}–{max(err):.1%} of the largest rate: the slaving holds while the centroid's modes are faster than the "
                "relative mode, its error growing as they near it"
                + ("" if not all('law_y' in j for r, j in js) else
                   f". Referee round 72: the slaving drops the centroid's offset y, slaved at first order in D (y = −(ww_e(w − w_e)/"
                   "(W(w² + w_e²)))R⁻¹QD, Q = ⟨r, ∇²_cφ⟩, R = Q − W⟨∇_cφ∇_cφᵀ⟩), and its first correction is λ ≈ κλ_D(1 − ε_b), ε_b = "
                   "(w − w_e)²ww_e/(w² + w_e²)²·R⁻¹Q = (1 − κ)λ_D/λ_b (λ_b the centroid's rate): ε_b (the larger) "
                   f"{min(max(j['eps_b']) for r, j in js):.1%}–{max(max(j['eps_b']) for r, j in js):.1%} over these rows against the "
                   f"slaved law's error {min(err):.1%}–{max(err):.1%}; the (D, y) system on the block's moments, y unslaved, within "
                   f"{max(erry):.1e} of the Jacobian's rates (relative to the largest)") + "."]
    for top in [t for t in TOPS_RP if t != 1.0]:
        pairs = [(a, b) for a in std for b in RP if b['top'] == top and b['prec'] == 64 and b['name'] == 'near tie' and b['x'] == a['x']
                 and b.get('win') is None]
        if pairs:
            out.append(f"The top edge at {top:g} ({edge(top):.2f}σ; N = {pairs[0][1].get('N', N0)}): ∂c̄*/∂W along y {pairs[0][1]['dcdW'][1]:.3g}σ against "
                       f"{pairs[0][0]['dcdW'][1]:.3g}σ; " + '; '.join(
                           f"x = {a['x']:g}: D_y at s = 700 {b['D'][1][1]:.3g} (the law {b['law'][1][1]:.3g}) against {a['D'][1][1]:.3g}, "
                           f"the side {b['side']}, s_d×N₀/N " + ('none' if b['s_d'] is None else f"{b['s_d'] * N0 / b.get('N', N0):.4g}")
                           + " against " + ('none' if a['s_d'] is None else f"{a['s_d']:.4g}") for a, b in pairs) + ".")
    p80 = [(a, b) for a in RP for b in RP if a['top'] == b['top'] and a['x'] == b['x'] and a['name'] == b['name'] and a['prec'] == 64
           and b['prec'] == 80 and a.get('win') is None and b.get('win') is None]
    out += [f"Long double against float64 (x = {a['x']:g}, top {a['top']:g}): D at s = 700 {vec(b['D'][1])} against {vec(a['D'][1])}, "
            f"the side {b['side']} against {a['side']}, s_d " + ('none' if b['s_d'] is None else f"{b['s_d']:.4g}") + " against "
            + ('none' if a['s_d'] is None else f"{a['s_d']:.4g}") + "." for a, b in p80]
    out.append("")
    return out


def rc_lines(RC, RP=None):
    """RC's report: the relative mode in the continuum (continued about a merged block, injected into the race) per gap and target, its
    checks, and s_d against Part C's window's (RP) — the climb's end and the hold (referee round 66)."""
    if RC is None or not len(RC):
        return []
    N0 = window(1.0).N
    main_ = [r for r in RC if (r['sub'], r['ymode'], r['inj']) == (SUB_RC, 'y', INJ_RC)]
    nf = main_[0]['N']
    pc = lambda v: '—' if v is None else f"{v:.4g} ({v * N0 / nf:.4g})"
    out = ["## RC. The relative mode in the continuum: continued about a merged block, injected into the race (referee rounds 66, 72)",
           "",
           "The continuum here is the far window's Jacobian with the rest's drift dropped (referee round 72: the drift as a forcing, "
           "below and RD).", "",
           f"On the window with every edge at least {min(edges(WIN_FAR)):.2f}σ from the block (N = {nf}; s in its time, in brackets "
           f"times N₀/N: Part C's, N₀ = {N0}) R's race (`rp_job`'s) to s₁ = {S1_RC:g}; there the merged block — both kernels at the "
           "centroid c_b less its offset y = (∂F/∂c_b)⁻¹F(c_b) off the merged rest — runs on (`merged_rhs`), and the centres' "
           "perturbation (δc_e = y − (w/W)D, δc_r = y + (w_e/W)D) by the race's exact Jacobian's centre block about it, in midpoint "
           f"matrix exponentials ({SUB_RC} per output, the norm in logarithms: no floor), until D = δc_r − δc_e, past its least, has "
           f"regrown to {INJ_RC:g}σ (and to twice its least); there it is injected into the float64 race, run to the departure (2σ) or, "
           f"on the far side, to its rests ({SIDE_R:g}σ beyond feature 0 away from feature 2). The least |D| (log₁₀ σ), its share along e "
           "and its s; the contraction's e-folds from s₁; the injection's s; the side; the climb's end (|D| = "
           f"{D_CLIMB:g}σ) and s_d:", "",
           "| target | x | least \\|D\\|: log₁₀ σ (share along e), at s | e-folds | injected at s | side | climb's end | s_d |",
           "|---|---|---|---|---|---|---|---|"]
    for r in sorted(main_, key=lambda r: (r['name'], r['x'])):
        out.append(f"| {r['name']} | {r['x']:g} | {r['least']:.2f} ({r['e_least']:+.2f}), {r['s_least']:.4g} | {r['efolds']:.0f} | "
                   + ('—' if r['s_inj'] is None else f"{r['s_inj']:.4g}") + f" | {r['side']} | {pc(r['s_climb'])} | {pc(r['s_d'])} |")
    fz = [r for r in main_ if r.get('forced', {}).get('least')]
    if fz:
        lab = lambda t: '—' if t is None or t[0] is None else f"{t[0]:.2f} ({t[1]:+.2f})"
        out += ["", "The rest's drift as a forcing (referee round 72): the merged base is no solution of the race — its kernels' speeds go "
                "as their coefficients —, so a block whose rest moves with W forces the relative mode by κk(∂c̄*/∂W)dW/ds, k = 2(w − w_e)/W "
                "(δc_e by −w/W of it, δc_r by +w_e/W). The forced D, log₁₀|D| σ (its share along e), at the continuum's least and at its "
                "injection (the forced part's size there over the homogeneous one's, inj·σ: the linear response, its sign deciding the "
                "side wherever it is the larger), for three drifts: the float64 sums' (`PairSums.block` along the merged run), exact sums on "
                "these samples (their "
                f"float64 coordinates: ∂c̄*/∂W = ({fz[0]['d_coord'][0]:+.3e}, {fz[0]['d_coord'][1]:+.3e})σ per unit W) and on the ideal "
                f"grid (the aliasing: ({fz[0]['d_alias'][0]:+.3e}, {fz[0]['d_alias'][1]:+.3e})σ) — the exact drifts at s₁ and "
                f"{S_CONT_RC:g} within {max(r['d_change'] for r in fz):.1e} relative —, against the continuum's least (the first column):",
                "", "| target | x | least \\|D\\| (share along e) | forced at the least: float64 sums | exact (these samples) | ideal grid "
                "(the aliasing) | at the injection, the forced over the homogeneous (dex; share): float64 sums | exact | ideal grid |",
                "|---|---|---|---|---|---|---|---|---|"]
        for r in sorted(fz, key=lambda r: (r['name'], r['x'])):
            fl, fi = r['forced']['least'], r['forced'].get('inj', {})
            over = lambda t: '—' if t is None or t[0] is None else f"{t[0] - np.log10(r['inj']):+.1f} ({t[1]:+.2f})"
            out.append(f"| {r['name']} | {r['x']:g} | {r['least']:.2f} ({r['e_least']:+.2f}) | {lab(fl['f64'])} | {lab(fl['coord'])} | "
                       f"{lab(fl['alias'])} | {over(fi.get('f64'))} | {over(fi.get('coord'))} | {over(fi.get('alias'))} |")
        dom = lambda k: sorted(r['x'] for r in fz if r['name'] == 'near tie' and r['forced']['least'][k][0] is not None
                               and r['forced']['least'][k][0] > r['least'])
        out += ["", "On the near tie the forced D exceeds the continuum's least at x = " + ', '.join(f"{x:g}" for x in dom('coord'))
                + " for the samples' coordinates and at x = " + (', '.join(f"{x:g}" for x in dom('alias')) or 'none') + " for the aliasing."]
    ck = [r for r in RC if r not in main_]
    ref = {(r['name'], r['x']): r for r in main_}
    if ck:
        out += ["", "Checks (the same gap; 16 substeps; y = 0 in the perturbation; the injection at another level): " + '; '.join(
            f"{r['name']}, x = {r['x']:g}, {r['sub']} substeps, " + ('y' if r['ymode'] == 'y' else 'y = 0') + f", {r['inj']:g}σ: least "
            f"{r['least']:.2f} against {ref[(r['name'], r['x'])]['least']:.2f}, s_d " + ('—' if r['s_d'] is None else f"{r['s_d']:.6g}")
            + " against " + ('—' if ref[(r['name'], r['x'])]['s_d'] is None else f"{ref[(r['name'], r['x'])]['s_d']:.6g}") + f", the side "
            f"{r['side']}" for r in ck if (r['name'], r['x']) in ref) + "."]
    if RP is not None and len(RP):
        far = {(q['name'], q['x']): q for q in RP if q.get('win') is not None and list(q['win']) == list(WIN_FAR)}
        std = {q['x']: q for q in RP if q.get('win') is None and q['top'] == 1.0 and q['prec'] == 64 and q['name'] == 'near tie'}
        val = [(r, far[(r['name'], r['x'])]) for r in main_ if (r['name'], r['x']) in far and far[(r['name'], r['x'])]['least'] > INJ_RC
               and r['s_d'] is not None and far[(r['name'], r['x'])]['s_d'] is not None]
        if val:
            out += ["", "Where float64 resolves the offset itself (RP's least |D| on this window above the injection's level) the injected "
                    "race's s_d against the plain run's: " + '; '.join(f"{r['name']}, x = {r['x']:g}: {r['s_d']:.5g} against {q['s_d']:.5g}"
                                                                       for r, q in val) + "."]
        cmp_ = [(r, std[r['x']]) for r in main_ if r['name'] == 'near tie' and r['x'] in std and r['s_d'] is not None
                and std[r['x']]['s_d'] is not None and std[r['x']].get('s_climb') is not None and r['s_climb'] is not None]
        if cmp_:
            out += ["", "s_d is the climb's end plus the hold after it: on Part C's window (RP, top edge at 1) the climb starts from the pin, "
                    "in the continuum from the layer's offset. In Part C's time — the climb's end, s_d and their shifts, the continuum's "
                    "less the pin's: " + '; '.join(
                        f"x = {r['x']:g}: climb {q['s_climb']:.4g} → {r['s_climb'] * N0 / nf:.4g}, s_d {q['s_d']:.4g} → "
                        f"{r['s_d'] * N0 / nf:.4g} ({(r['s_d'] * N0 / nf) / q['s_d'] - 1:+.1%}; shifts {r['s_d'] * N0 / nf - q['s_d']:.3g} "
                        f"and {r['s_climb'] * N0 / nf - q['s_climb']:.3g})" for r, q in sorted(cmp_, key=lambda t: t[0]['x'])) + "."]
        th = [(r, q) for r, q in cmp_ if r.get('s_climbs') and q.get('s_climbs') and None not in r['s_climbs'] + q['s_climbs']]
        if th:
            out += ["", "The split is a convention, the shift an identity (referee round 72): s_d less the climb's end at |D| = θσ, θ = "
                    + ', '.join(f"{t:g}" for t in THETAS_CLIMB) + ", on Part C's window and in the continuum (Part C's time), and their "
                    "relative difference: " + '; '.join(
                        f"x = {r['x']:g}: " + ', '.join(
                            f"{q['s_d'] - q['s_climbs'][i]:.4g} and {(r['s_d'] - r['s_climbs'][i]) * N0 / nf:.4g} "
                            f"({((r['s_d'] - r['s_climbs'][i]) * N0 / nf) / (q['s_d'] - q['s_climbs'][i]) - 1:+.1e})"
                            for i in range(len(THETAS_CLIMB))) for r, q in sorted(th, key=lambda t: t[0]['x'])) + "."]
    wr = [r for r in main_ if r.get('w_peak') is not None and r['side'] == 'toward']
    if wr:
        out += ["", "The runner-up's coefficient along the injected race toward feature 2 (referee round 72), relative to its value at the "
                "injection: at the climb's ends θ = " + ', '.join(f"{t:g}" for t in THETAS_CLIMB) + ", its peak (|D| there; s in Part C's "
                "time) and where it is 1% below: " + '; '.join(
                    f"{r['name']}, x = {r['x']:g}: " + ', '.join('—' if v is None else f"{v:+.2e}" for v in r['w_climbs'])
                    + f"; {r['w_peak'][0]:+.2e} at |D| = {r['w_peak'][2]:.3f}σ (s = {r['w_peak'][1] * N0 / nf:.4g}); "
                    + ('none' if r['s_wdecay'] is None else f"{r['s_wdecay'] * N0 / nf:.4g}")
                    for r in sorted(wr, key=lambda r: (r['name'], r['x']))) + "."]
    out.append("")
    return out


def rd_lines(RD):
    """RD's report (referee round 72): ∂c̄*/∂W at RC's merged block in each arithmetic and on each grid."""
    if RD is None or not len(RD):
        return []
    q = lambda ds: np.exp(-np.pi ** 2 * SIG ** 2 / ds ** 2)
    vec = lambda v: f"({v[0]:+.3e}, {v[1]:+.3e})"
    out = ["## RD. What moves the block's rest with W without an edge (referee round 72)", "",
           f"RC's merged block (the gap x on the window x ∈ [{WIN_FAR[0]:g}, {WIN_FAR[1]:g}], y ∈ [{WIN_FAR[2]:g}, {WIN_FAR[3]:g}]) at "
           f"s = {', '.join(f'{v:g}' for v in S_RD)}: the rest's shift ∂c̄*/∂W = (∂F/∂c̄)⁻¹⟨φ, ∇_cφ⟩ (σ per unit W; in the continuum 0, "
           "⟨φ, ∇_cφ⟩ = ½∇_c‖φ‖²) with the sums in float64, in the long double and exact "
           f"({RD_PREC} digits) on the committed samples (spacing {DS:g}, `np.arange`'s float64 coordinates, off the ideal grid by up to "
           f"{max(RD[0]['coord_err']):.2e}), exact on the ideal grid against the aliasing's first harmonic (Poisson summation: "
           f"⟨φ, ∂_iφ⟩ = −⟨φ²⟩(2πq/ds)sin(2π(c_i − x₀,ᵢ)/ds), q = e^{{−π²σ²/ds²}} = {q(DS):.3e}), and the same on a dyadic grid "
           f"(spacing {DS_DYADIC:g}, its float64 coordinates exact: off by {max(RD[0]['dy_coord_err']):g}; q = {q(DS_DYADIC):.3e}):", "",
           "| target | x | s | W | float64 | long double | exact | ideal grid, exact | the aliasing | dyadic: float64 | long double | exact | "
           "the aliasing |", "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in RD:
        for w in r['rows']:
            out.append(f"| {r['name']} | {r['x']:g} | {w['s']:g} | {w['W']:.4f} | " + ' | '.join(
                vec(w[k]) for k in ('f64', 'ld', 'exact', 'ideal', 'alias', 'dy_f64', 'dy_ld', 'dy_exact', 'dy_alias')) + " |")
    ex = [w for r in RD for w in r['rows']]
    rel = lambda a, b: max(abs(x_ - y_) for x_, y_ in zip(a, b)) / max(abs(y_) for y_ in b)
    out += ["", f"On the committed samples the long double's drift is the exact sums' to {max(rel(w['ld'], w['exact']) for w in ex):.1e} "
            "relative — the drift is the samples' coordinates', in any arithmetic — and float64's differs from it by up to "
            f"{max(max(abs(a - b) for a, b in zip(w['f64'], w['exact'])) for w in ex):.1e}σ, its sums' rounding; on the ideal grid the exact drift is the aliasing's to {max(rel(w['alias'], w['ideal']) for w in ex):.1e}, "
            f"on the dyadic grid to {max(rel(w['dy_alias'], w['dy_exact']) for w in ex):.1e}, where the long double's "
            f"(|·| ≤ {max(max(abs(v) for v in w['dy_ld']) for w in ex):.1e}σ) and float64's (≤ {max(max(abs(v) for v in w['dy_f64']) for w in ex):.1e}σ) "
            "are their sums' rounding.", ""]
    return out


def race_rule(sd, s2):
    """The race's prediction (v): a runner-up displaced toward feature 2 is captured iff it departs before feature 2's own entry."""
    return 'captured' if sd is not None and sd < s2 else 'not departed'


def r_lines(R, C3=None, C4=None, C4b=None, FC=None, C4M=None):
    """R's report: the ε-free race on both sides, the rule against C3's and C4's outcomes, the side's tolerance scan, FC's folds."""
    if R is None:
        return []
    out = ["## R. The race without ε (referee round 53)", "",
           f"C2's near tie at μ → 0 ({MU_R:g}) with only feature 0's entrant and runner-up active and features 1 and 2 in the residual: the "
           f"inner problem (the entrant {XI0_R:g} behind its wall, the runner-up x behind it) in s = (τ − τ_e)/ε^{{1/2}} to {S_END_R:g} — C4 "
           f"without its horizon s₂ —, rtol {RTOL_R:g}: the side the runner-up leaves feature 0 to (its displacement along e where it first "
           f"passes {SIDE_R:g}σ from feature 0: toward feature 2, or away), "
           "the pair's closest approach in the inner layer, its coefficient after the layer (s = 10³), its departure s_d (2σ from feature 0) "
           "and its coefficient there:", "",
           "| x | side | closest approach (σ; its share along e) | w after the layer | s_d | w at s_d | at s = 10⁵: w, distance (σ) |",
           "|---|---|---|---|---|---|---|"]
    for r in sorted(R['main'], key=lambda r: r['x']):
        i5 = r['marks'].index(1e5) if 1e5 in r['marks'] else None
        m5 = '—' if i5 is None else f"{r['w_marks'][i5]:.4g}, {r['d_marks'][i5]:.3f}"
        out.append(f"| {r['x']:g} | {r['side']} | {r['closest']:.1e}" + (f" ({r['e_closest']:+.2f})" if r.get('e_closest') is not None else '')
                   + f" | {r['w_layer']:.4g} | "
                   + (f"{r['s_d']:.4g} | {r['w_dep']:.4g}" if r['s_d'] is not None else f"none by {r['s_end']:.2g} | —") + f" | {m5} |")
    sd = {round(r['x'], 3): r['s_d'] for r in R['main']}
    rows, agree, ties = [], 0, 0
    tau21 = (np.median([r['s2'] * np.sqrt(r['eps']) for r in C4]) if C4 is not None else None)
    for lab, rs in (('C4', C4), ('C4b', C4b), ('C3', [r for r in (C3 if C3 is not None else []) if r['name'] == 'near tie'])):
        for r in (rs if rs is not None else []):
            x = float(r['x']) if 'x' in r else float(r['mu'] / np.sqrt(r['eps']))
            s2 = r['s2'] if 's2' in r else (tau21 / np.sqrt(r['eps']) if tau21 is not None else None)
            k = round(x, 3)
            if s2 is None or k not in sd:
                continue
            pred, act = race_rule(sd[k], s2), c3_outcome(r)
            if act == 'capture tie':
                ties += 1
            elif (pred == 'captured') == (act == 'captured'):
                agree += 1
            rows.append((lab, x, r['eps'], s2, sd[k], pred, act))
    if rows:
        n = len(rows) - ties
        out += ["", f"The rule — captured iff s_d(x) < s₂ = (τ₂ − τ_e)ε^{{−1/2}} (C3's s₂ from C4's τ₂ − τ_e = "
                    f"{tau21:.4g}) — against the loops' outcomes: {agree} of {n} agree, {ties} capture ties:", "",
                "| loops | x | ε | s₂ | s_d | predicted | the loop's outcome |", "|---|---|---|---|---|---|---|"]
        out += [f"| {lab} | {x:.3g} | {e:g} | {s2:.3g} | " + (f"{d:.3g}" if d is not None else 'none') + f" | {p} | {a} |"
                for lab, x, e, s2, d, p, a in sorted(rows, key=lambda q: (q[0], q[1], -q[2]))]
    if R.get('mirror') is not None and len(R['mirror']):
        out += ["", "The mirrored target (feature 2 across the face: the runner-up held on the far side): the side, the closest approach, "
                    "the coefficient after the layer and at s = " + ', '.join(f"{q:g}" for q in S_MARKS_R) + ", the distance from feature 0 "
                    "at the end, the first s at which a coefficient below θ|A| has grown past it, the departure:", "",
                "| x | side | closest (σ) | w after the layer | w at the marks | distance at the end (σ) | s at θ\\|A\\| | departure |",
                "|---|---|---|---|---|---|---|---|"]
        for r in sorted(R['mirror'], key=lambda r: r['x']):
            out.append(f"| {r['x']:g} | {r['side']} | {r['closest']:.1e} | {r['w_layer']:.4g} | {', '.join(f'{w:.3g}' for w in r['w_marks'])} | "
                       f"{r['d_end']:.3f} | " + (f"{r['s_theta']:.3g}" if r['s_theta'] is not None else '—') + " | "
                       + (f"{r['s_d']:.3g}" if r['s_d'] is not None else 'none') + " |")
    if R.get('long') is not None and len(R['long']):
        out += ["", "The mirror past R's horizon (referee round 58): the coefficients at the end, the s at which the runner-up outgrows "
                    "the entrant (the swap), the entrant's departure (2σ from feature 0) and the feature it ends nearest (its distance), "
                    "the runner-up's distance from feature 0 at the end:", "",
                "| x | s_end | w: runner-up, entrant | the swap s | the entrant's departure s | the entrant ends at feature (σ) | the "
                "runner-up's distance (σ) |", "|---|---|---|---|---|---|---|"]
        for r in sorted(R['long'], key=lambda r: r['x']):
            out.append(f"| {r['x']:g} | {r['s_end']:.3g} | {r['w_end']:.4g}, {r['w_e_end']:.4g} | "
                       + (f"{r['s_swap']:.3g}" if r.get('s_swap') is not None else 'none') + ' | '
                       + (f"{r['s_de']:.3g}" if r.get('s_de') is not None else 'none') + f" | {r['e_dest']} ({r['e_dest_d']:.3f}) | "
                       f"{r['d_end']:.3f} |")
        th = {r['x']: r.get('s_theta') for r in (R.get('mirror') if R.get('mirror') is not None else [])}
        dts = [r['s2'] * np.sqrt(r['eps']) for r in (C4M if C4M is not None else [])]
        win = [(r['x'], th.get(r['x']), r.get('s_de')) for r in sorted(R['long'], key=lambda r: r['x'])]
        if dts and any(a is not None and b is not None for x, a, b in win):
            dt = float(np.median(dts))
            out += ["", f"The far rule's window in ε (referee round 62): s₂ = (τ₂ − τ_e)ε^{{−1/2}} with τ₂ − τ_e = {dt:.4g} (C4M's path, "
                    f"{min(dts):.4g}–{max(dts):.4g}), so s_θ < s₂ < s_{{d,e}} holds for ((τ₂ − τ_e)/s_{{d,e}})² < ε < ((τ₂ − τ_e)/s_θ)²: "
                    + '; '.join(f"x = {x:g}: {(dt / b) ** 2:.2g} < ε < {(dt / a) ** 2:.2g}" for x, a, b in win if a is not None and b is not None)
                    + "."]
    if R.get('scan') is not None and len(R['scan']):
        rts = sorted({r['rtol'] for r in R['scan']}, reverse=True)
        out += ["", f"The side against the tolerance (to s = {S_END_RS:g}): per gap the closest approach at the tightest tolerance and the "
                    "side at rtol " + ', '.join(f"{q:g}" for q in rts) + f" (None: the pair merged — the runner-up never {SIDE_R:g}σ "
                    "from feature 0):", "",
                "| x | closest approach (σ) | sides |", "|---|---|---|"]
        for x in sorted({r['x'] for r in R['scan']}):
            rs = {r['rtol']: r for r in R['scan'] if r['x'] == x}
            out.append(f"| {x:g} | {rs[min(rts)]['closest']:.1e} | " + ', '.join(str(rs[q]['side']) for q in rts) + " |")
        out.append("")
        out.append("Runs not toward, per tolerance: " + '; '.join(f"{q:g}: {sum(r['side'] != 'toward' for r in R['scan'] if r['rtol'] == q)} of "
                                                                f"{sum(1 for r in R['scan'] if r['rtol'] == q)}" for q in rts) + ".")
    if FC is not None:
        def fc_line(r):
            lab = {('near tie', 1): 'toward feature 2', ('near tie', -1): 'away from it', ('alone', 1): 'toward, feature 1 removed'}
            if r['fold'] is not None:
                return (f"{lab.get((r['name'], r['side']))}: the fold at w = {np.sqrt(r['fold'][0] * r['fold'][1]):.5g}, the rest there "
                        f"{r['at'][0]:.3f}σ from feature 0 ({r['at'][1]:+.3f}σ along e)")
            return (f"{lab.get((r['name'], r['side']))}: no fold down to w = {r['w_end']:g}"
                    + (f" (the rest {r['at'][0]:.3f}σ, {r['at'][1]:+.3f}σ along e)" if r['at'] is not None else ''))
        out += ["", "FC. The static fold in C4's own geometry (Part C's samples; the entrant and the runner-up's centre free, its coefficient "
                    "held): " + '; '.join(fc_line(r) for r in FC) + "."]
    out.append("")
    return out


def report(A=None, B=None, C=None, T=None, D=None, B4=None, H_=None, T2=None, W=None, C2=None, V=None, C3=None, T3=None, B4u=None,
           C4=None, B4x=None, V3=None, D2=None, C4b=None, R=None, VF=None, FC=None, C4M=None, RP=None, C4W=None, RC=None, RD=None):
    L = law()
    kap = np.pi * (2 * L.c) ** 0.5
    L_ = ["# Proposition D‴.54: the transit in ℝⁿ (`exp_transit_nd.py`)", "",
          f"The free law with moving centres (m = {L.m:g}, G = {L.G:g}, γ = {L.gamma:g}, a_max = {L.A:g}, c = {L.c:g}), Gaussian kernels of width "
          f"{SIG:g} in the plane on a square lattice of spacing {H:g} (H: a hexagonal one), the samples on a grid of spacing {DS:g}; masses ε; "
          f"formed: |w| > {TH_REL:g}|A|; κ = πρ(2c)^{{1/2}} = {kap:.4f}.", ""]
    if A is not None:
        L_ += ["## A. The k-fold near tie", "",
               f"One feature near an edge (the kernels at (0, 0), ({H:g}, 0)) or a vertex (and ({H:g}, 0), (0, {H:g}), ({H:g}, {H:g})) of the "
               f"lattice — 8 × 8 kernels around it, the samples symmetric about the vertex —, its margins ε^{{1/2}}x at the entry; the loop to "
               f"{TAU_A:g} τ_entry. The gaps x (the runner-ups', as placed; at the vertex the diagonal's as it follows), the shares of the "
               "tie's kernels in the loop against the inner problem's at the gaps of the smallest ε (the feature at the tie; extrapolated in "
               f"its start from xi0 = {', '.join(f'{v:g}' for v in XI0_A)}, the change between those in the last column), the softmax at those "
               "gaps and, at the vertex, the product of the two edges' inner problems posed in the vertex's geometry (the feature H/2 off each "
               "edge's axis):", "",
               "| tie | A | x | shares: loop (ε = " + ', '.join(f'{e:g}' for e in EPS_A) + ") | inner problem | softmax | the edges' product | "
               "the others' largest \\|w\\| (ε = " + ', '.join(f'{e:g}' for e in EPS_A) + ") | the loop's successive changes | xi0's change |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        inner = {(r['kind'], r['A'], tuple(np.round(r['xs'], 9)), round(float(r.get('t', 0.0)), 9)): r for r in A['inner']}
        groups = {}
        for r in A['loops']:
            groups.setdefault((r['kind'], r['A'], tuple(r['x'])), []).append(r)
        worst, succ, prod_a = {}, [], {}
        for key, rs in groups.items():
            rs = sorted(rs, key=lambda r: -r['eps'])
            kind, Aa, x = key
            xa = np.r_[0.0, rs[-1]['x_act'][1:]]
            ri = inner.get((kind, Aa, tuple(np.round(xa, 9)), 0.0))
            sm = softmax(xa, L)
            prod = '—'
            if kind == 'vertex':
                e1 = inner.get(('edge', Aa, tuple(np.round([0.0, xa[1]], 9)), round(H / 2, 9)))
                e2 = inner.get(('edge', Aa, tuple(np.round([0.0, xa[2]], 9)), round(H / 2, 9)))
                if e1 is not None and e2 is not None:
                    p1, p2 = e1['shares'][1], e2['shares'][1]
                    pv = np.array([(1 - p1) * (1 - p2), p1 * (1 - p2), (1 - p1) * p2, p1 * p2])
                    prod = ', '.join(f'{v:.4g}' for v in pv)
                    if ri is not None:
                        worst.setdefault(('prod', kind), []).append(float(np.abs(pv - ri['shares']).max()))
                        prod_a.setdefault(Aa, []).append((float(np.abs(pv - ri['shares']).max()), float(np.abs(pv - ri['shares'])[1:3].max())))
            if ri is not None:
                worst.setdefault(('loop', kind), []).append(float(np.abs(rs[-1]['shares'] - ri['shares']).max()))
                worst.setdefault(('soft', kind), []).append(float(np.abs(sm - ri['shares']).max()))
            dif = [float(np.abs(a_['shares'] - b_['shares']).max()) for a_, b_ in zip(rs[:-1], rs[1:])]
            succ.append(dif)
            spread = f"{ri['spread']:.1e}" if ri is not None and 'spread' in ri else '—'
            L_.append(f"| {kind} | {Aa:g} | {', '.join(f'{round(float(v), 3) + 0.0:.3f}' for v in xa[1:])} | "
                      + '; '.join(', '.join(f'{v:.4g}' for v in r['shares']) for r in rs) + " | "
                      + (', '.join(f'{v:.4g}' for v in ri['shares']) + ('' if ri['stationary'] else ' (not stationary)') if ri is not None else '—')
                      + f" | {', '.join(f'{v:.4g}' for v in sm)} | {prod} | " + ', '.join(f"{r['rest']:.2e}" for r in rs) + " | "
                      + ', '.join(f'{v:.1e}' for v in dif) + f" | {spread} |")
        rat = [a_ / b_ for d in succ for a_, b_ in zip(d[:-1], d[1:]) if b_ > 0]
        L_ += ["", "The largest difference of the shares from the inner problem's (at the smallest ε): " + '; '.join(
            f"{k[1]} — {dict(loop='the loop', soft='the softmax', prod='the edges’ product')[k[0]]} {max(v):.2e}" for k, v in sorted(worst.items()))
               + f". The loop's successive changes fall by {np.median(rat):.2f} per decade of ε (median; √10 = 3.16), their last {max(d[-1] for d in succ):.1e} "
                 f"at most" + (f"; the inner problem's change between its starts at most {max(r['spread'] for r in A['inner']):.1e}"
                               if all('spread' in r for r in A['inner']) else '') + "." +
               ("" if not prod_a else " The edges' product against the 4-kernel problem, by amplitude (all shares; the edges' runner-ups only — the rest is the diagonal's "
                "frozen tail, which does not factorise, and the entrant's balance): " +
                '; '.join(f"A = {a:g}: {max(v[0] for v in prod_a[a]):.1e} ({max(v[1] for v in prod_a[a]):.1e})" for a in sorted(prod_a)) + "."), ""]
    if T is not None and 'pt' in T:
        xt = x_theta_t(T)
        tt = np.asarray(T['t'])
        xbar = x_theta_avg(T)
        L_ += ["## T. The edge's split along the face", "",
               f"The 2-kernel inner problem (A = {A_B:g}) with the feature on the face at the offset t from the edge's axis: the runner-up's share "
               f"p(x; t) at x = 2, 4, 6, 8 and the crossing x_θ(t) of {TH_REL:g}:", "",
               "| t/H | p(2) | p(4) | p(6) | p(8) | x_θ(t) |", "|---|---|---|---|---|---|"]
        for i, t_ in enumerate(tt):
            L_.append(f"| {t_ / H:g} | " + ' | '.join(f"{np.interp(x, T['x'], T['pt'][i]):.5f}" for x in (2, 4, 6, 8)) + f" | {xt[i]:.3f} |")
        L_ += ["", f"x_θ averaged over the face (t uniform in [−H/2, H/2]): {xbar:.3f} (at t = 0: {xt[0]:.3f}; at the vertex, t = H/2: {xt[-1]:.3f}). "
                   f"For a feature uniform in a Voronoi cell the offsets from a face are (2σ²/H)ε^{{1/2}}x at leading order, so P(N ≥ 2) = "
                   f"4·(2σ²/H²)·x̄_θ·ε^{{1/2}} = {e_extra(xbar, 2):.1f}·ε^{{1/2}} (the four faces' layers, less O(ε) where they overlap).", ""]
        if T2 is not None:
            ta = np.r_[tt, T2['t']]
            xa = np.r_[xt, T2['x']]
            fa = {n: face_average(ta, xa, n) for n in (2, 3, 4, 5)}
            L_ += [f"Beyond H/2 (T2, by bisection): x_θ = " + ', '.join(f"{x:.3f} (t = {t / H:g}H)" for t, x in zip(T2['t'], T2['x'])) +
                   ". In the continuum x_θ depends on |t| only, so on hℤⁿ (σ/H as here) a face is the (n − 1)-cube of offsets and E[N − 1] = "
                   "2n(2σ²/H²)x̄_θ⁽ⁿ⁾ε^{1/2}: " + '; '.join(f"n = {n}: x̄_θ = {fa[n]:.3f}, {e_extra(fa[n], n):.1f}·ε^{{1/2}}" for n in fa) +
                   " (the n = 2 average by the midpoint rule; the trapezoid rule above).", ""]
    if B is not None:
        xt, where = B['x_theta'], "the edge's inner problem on its axis"
        if T is not None and 'pt' in T:                                # at a vertex the feature sits H/2 off each edge's axis (round 42)
            xt, where = float(x_theta_t(T)[-1]), "the edge's inner problem at the vertex's offset t = H/2 (T)"
        L_ += ["## B. The formation map at a vertex", "",
               f"A = {A_B:g}; the gaps (x₁₀, x₀₁) of the edge runner-ups on a grid, x₁₀ ≤ x₀₁, the diagonal's x₁₁ = x₁₀ + x₀₁ − ε^{{1/2}}x₁₀x₀₁ "
               "exactly (the tensor grid's correlations multiply); the number of kernels formed by the loop at ε = "
               f"{', '.join(f'{e:g}' for e in EPS_B)}, by the 4-kernel inner problem, by the independent rule #{{j: x_j < x_θ}} + 1 with "
               f"x_θ = {xt:.3f} ({where}), and by the product rule with the edges' factors at the vertex's offset "
               "t = H/2 (T):", "",
               "| x₁₀ | x₀₁ | x₁₁ | loop | inner problem | the independent rule | the product rule (T) |", "|---|---|---|---|---|---|---|"]
        byx = {}
        for r in B['loops']:
            byx.setdefault(tuple(r['x']), {})[r['eps']] = r['n']
        inn = {tuple(r['x']): r for r in B['inner']}
        if B4 is not None:                                             # the 4-kernel problem settled to SHARE_TOL (B4) at B's points
            b4 = {(round(r['x'][0], 9), round(r['x'][1], 9)): r for r in B4['rows']}
            if all((round(x[0], 9), round(x[1], 9)) in b4 for x in inn):
                inn = {x: b4[(round(x[0], 9), round(x[1], 9))] for x in inn}
                L_[-4] = L_[-4].replace("by the 4-kernel inner problem,", "by the 4-kernel inner problem (B4's),")
        agree = {e: 0 for e in EPS_B}; agree_rule = 0; agree_prod = 0; tot = 0
        for x in sorted(byx):
            ri = inn[x]
            rule = 1 + int(x[0] < xt) + int(x[1] < xt) + int(ri['xd'] < xt)
            tot += 1; agree_rule += int(rule == ri['n'])
            for e in EPS_B:
                agree[e] += int(byx[x].get(e) == ri['n'])
            pr = product_count(T, *x) if T is not None else None
            agree_prod += int(pr == ri['n'])
            L_.append(f"| {x[0]:g} | {x[1]:g} | {ri['xd']:.2f} | {', '.join(str(byx[x].get(e, '—')) for e in EPS_B)} | {ri['n']} | {rule} | "
                      f"{pr if pr is not None else '—'} |")
        L_ += ["", f"The loop's count equals the inner problem's at {', '.join(f'{agree[e]} of {tot} (ε = {e:g})' for e in EPS_B)}; the "
               f"independent rule at {agree_rule} of {tot}" + (f", the product rule at {agree_prod} of {tot}" if T is not None else '') + ".", ""]
    if B4 is not None:
        g, N = vertex_counts(B4)
        a3, a4 = areas(N, g)
        line = (f"## B4. The areas at a vertex", "",
                f"The 4-kernel inner problem on the grid of gaps 0 … {g[-1]:g} in steps of {g[1] - g[0]:g} ({len(B4['rows'])} problems, x₁₀ ≤ x₀₁): the "
                f"region N ≥ 3 has the area {a3:.1f} and N = 4 {a4:.1f} in the gaps' quarter-plane (trapezoid weights; "
                f"{sum(bool(r['stationary']) for r in B4['rows'])} of {len(B4['rows'])} problems settled to the shares' tolerance "
                f"{SHARE_TOL:g}" + (b4u_note(B4, B4u, B4x) if B4u is not None else '') + ")")
        L_ += list(line[:2])
        if T is not None and 'pt' in T:
            Np = np.array([[product_count(T, a, b_) for b_ in g] for a in g])
            p3, p4 = areas(Np, g)
            agree = int((Np == N).sum())
            gf = np.arange(0.0, float(g[-1]) + 1e-9, 0.05)
            f3, f4 = areas(np.array([[product_count(T, a, b_) for b_ in gf] for a in gf]), gf)
            L_.append(line[2] + f"; the product rule with the factors at t = H/2 on the same grid: {p3:.1f} and {p4:.1f} (on a grid of step 0.05: "
                      f"{f3:.1f} and {f4:.1f}, the coarse grid's error), its count equal to the "
                      f"inner problem's at {agree} of {N.size} grid points. Each vertex's sector in a cell maps to the gaps by (2σ²ε^{{1/2}}/H)·(x₁₀, "
                      f"x₀₁) and a cell has four: P(N ≥ 3) = 4(2σ²/H²)²·{a3:.1f}·ε = {p_vertices(a3):.0f}·ε and P(N = 4) = {p_vertices(a4):.0f}·ε.")
        else:
            L_.append(line[2] + f": P(N ≥ 3) = {p_vertices(a3):.0f}·ε and P(N = 4) = {p_vertices(a4):.0f}·ε.")
        L_.append("")
    if C is not None:
        L_ += ["## C. Several features in the plane", "",
               f"Random three-feature targets (seed {SEED_C}; centres in [−0.7, 0.7]² at least 3σ apart, amplitudes ±U(0.4, 1)) on the "
               f"lattice over [−1, 1]² ({int(round(2 / H)) + 1}² kernels; the samples' grid {int(round(2 / DS)) + 1}²); the sliding path against "
               "the loop to 1.3 times the path's last entry; the never-entering kernels' displacement vectors against the creep's εC_j "
               "integrated along the path (`creep_constants_nd`; where |εC_j| is above 10⁻³ of the largest: the largest relative error and the "
               "smallest cosine):", "",
               "| target | the path (entries) | ε | formed = the path's | the never-entering kernels' largest \\|w\\|·ε^{−3/2} | their largest "
               "displacement/ε | creep: relative error, smallest cosine | the block: max \\|w − w_path\\|, max \\|c − c_path\\| | rms residual |",
               "|---|---|---|---|---|---|---|---|---|"]
        for q in sorted(set(r['q'] for r in C)):
            rs = sorted([r for r in C if r['q'] == q], key=lambda r: -r['eps'])
            for r in rs:
                cr = f"{r['creep_rel']:.2e}, {r['creep_cos']:.6f}" if 'creep_rel' in r else '—'
                L_.append(f"| {q} | {r['path']} ({', '.join(f'{t_:.2f}' for t_, j, k in r['events'])}) | {r['eps']:g} | "
                          f"{'yes' if r['formed'] == r['path'] else 'no: ' + str(r['formed'])} | {r['w_never'] / r['eps'] ** 1.5:.4g} | "
                          f"{r['disp_never'] / r['eps']:.4g} | {cr} | {r['dw']:.2e}, {r['dc']:.2e} | {r['rms']:.2e} |")
            es = [r['eps'] for r in rs]
            L_.append(f"| {q} | slopes in log ε (10⁻⁸ … 10⁻¹⁰) | | | {slope(es[-3:], [r['w_never'] for r in rs][-3:]):.3f} | "
                      f"{slope(es[-3:], [r['disp_never'] for r in rs][-3:]):.3f} | | {slope(es[-3:], [r['dw'] for r in rs][-3:]):.3f}, "
                      f"{slope(es[-3:], [r['dc'] for r in rs][-3:]):.3f} | |")
        L_.append("")
    if D is not None and T is not None:
        K, c0 = local()
        tt = np.asarray(T.get('t', [0.0]))
        xbar = x_theta_avg(T) if 'pt' in T else float(np.interp(-TH_REL, -T['p'], T['x']))
        L_ += ["## D. The probabilities over uniformly placed features", "",
               f"A = {A_B:g} at {N_D} uniform positions in the Voronoi cell of the kernel at (0, 0) (seed {SEED_D}; the same positions at both "
               "masses, so the two rows are not independent); the loop's count N against the product rule at each position's gaps to its two "
               "nearest faces and its offsets along them (T), and the fractions against the rule's at the sampled positions and the asymptotics "
               f"(P(N ≥ 2) = 4w/H − 4w²/H², w = 2σ²x̄_θε^{{1/2}}/H with the face average x̄_θ = {xbar:.3f}; P(N ≥ 3) from B4's area):", "",
               "| ε | N = 1, 2, 3, 4 (loop) | the product rule agrees | P(N ≥ 2): loop (± s.e.) | the rule at the positions | asymptotic | "
               "P(N ≥ 3): loop | the rule at the positions | asymptotic |", "|---|---|---|---|---|---|---|---|---|"]
        a3 = None
        if B4 is not None:
            g, N4 = vertex_counts(B4)
            a3 = areas(N4, g)[0]
        blind = []
        for e in EPS_D:
            rs = [r for r in D if r['eps'] == e]
            n = np.array([r['n'] for r in rs])
            pr = np.array([product_count(T, r['gaps'][0], r['gaps'][1], r['cstar'][1] - c0[r['e']][1], r['cstar'][0] - c0[r['e']][0])
                           for r in rs])
            p2, p3 = float(np.mean(n >= 2)), float(np.mean(n >= 3))
            asy3 = f"{p_vertices(a3) * e:.4f}" if a3 is not None else '—'
            pr0 = np.array([product_count(T, r['gaps'][0], r['gaps'][1], 0.0, 0.0) for r in rs])
            blind.append(f"{int((n == pr0).sum())} of {len(rs)} (ε = {e:g})")
            L_.append(f"| {e:g} | {', '.join(str(int((n == j).sum())) for j in (1, 2, 3, 4))} | {int((n == pr).sum())} of {len(rs)} | "
                      f"{p2:.4f} ± {np.sqrt(p2 * (1 - p2) / len(rs)):.4f} | {np.mean(pr >= 2):.4f} | {p_layers(xbar, e):.4f} | {p3:.4f} | "
                      f"{np.mean(pr >= 3):.4f} | {asy3} |")
        if 'pt' in T:                                                  # the misses: the loop's count, the rule's, the gaps, x_θ at the offsets
            xt_ = x_theta_t(T)
            miss = []
            for e in EPS_D:
                for r in [r for r in D if r['eps'] == e]:
                    t1, t2 = r['cstar'][1] - c0[r['e']][1], r['cstar'][0] - c0[r['e']][0]
                    pr_ = product_count(T, r['gaps'][0], r['gaps'][1], t1, t2)
                    if pr_ != r['n']:
                        th = [float(np.interp(abs(t_), T['t'], xt_)) for t_ in (t1, t2)]
                        miss.append(f"{e:g}: {r['n']}/{pr_} at ({r['gaps'][0]:.2f}, {r['gaps'][1]:.2f}) [{th[0]:.2f}, {th[1]:.2f}]")
            L_ += ["", "The misses — ε: the loop's count/the rule's at the gaps (x₁₀, x₀₁) [x_θ at their offsets]: " + '; '.join(miss) + "."
                   if miss else "", "",
                   "The offset-blind rule (every factor at t = 0, x_θ(0)) agrees with the loop at " + ', '.join(blind) + ": at these "
                   "masses the loop's own threshold is shifted by more than x_θ varies along a face (D2), so D measures the shift, not "
                   "x_θ(t) (referee round 49).", ""]
        L_.append("")
    if H_ is not None:
        g = np.asarray(H_['grid'], float)
        idx = {v: i for i, v in enumerate(np.round(g, 9))}
        N = np.zeros((g.size, g.size), int)
        ex = None
        for r in H_['inner']:
            i, j = idx[round(r['x'][0], 9)], idx[round(r['x'][1], 9)]
            N[i, j] = N[j, i] = r['n']
            if r['x'] == (0.0, 0.0) or tuple(r['x']) == (0.0, 0.0):
                ex = r
        a3 = areas(N, g)[0]
        xth = [float(g[np.argmax(N[i] < 3)]) if (N[i] < 3).any() else np.nan for i in range(g.size)]
        L_ += ["## H. A hexagonal lattice", "",
               f"Kernels on a hexagonal lattice of spacing {H:g}, the feature (A = {A_B:g}) at a Voronoi vertex — three cells meet there, the "
               "kernels (0, 0), (H, 0), (H/2, √3H/2) at H/√3 —, placed where the three correlate equally with it on the samples. The 3-kernel "
               "inner problem at the exact tie: shares " + (', '.join(f'{v:.5f}' for v in ex['shares']) if ex is not None else '—') +
               f"; on the gaps' grid 0 … {g[-1]:g} (step {g[1] - g[0]:g}) the region where all three form, N = 3, has the area {a3:.1f} (the "
               f"square vertex's N ≥ 3: B4). A sector of a hexagonal vertex maps to the gaps with the Jacobian (2σ²ε^{{1/2}})²/((√3/2)H²) and a "
               f"cell has six: P(N = 3) = 32σ⁴/H⁴·{a3:.1f}·ε = {32 * SIG ** 4 / H ** 4 * a3:.0f}·ε. The loop at the exact tie:", "",
               "| ε | shares | formed | the others' largest \\|w\\| | the gaps of the runner-ups and of the fourth kernel (units of ε^{1/2}) |",
               "|---|---|---|---|---|"]
        for r in H_['loops']:
            L_.append(f"| {r['eps']:g} | {', '.join(f'{v:.5f}' for v in r['shares'])} | {r['n']} | {r['rest']:.2e} | "
                      f"{', '.join(f'{v:.4g}' for v in r['gaps'])} |")
        L_ += ["", "The fourth kernel's margin at the vertex is ε-free, 1 − e^{−H²/(4σ²)} in the continuum, so its gap grows like ε^{−1/2}: "
                   "a fourth kernel never forms at a simple vertex (k ≤ n + 1).", ""]
    if C2 is not None:
        Kf, cf = full()
        pt = lambda i: f"({cf[i][0]:.3f}, {cf[i][1]:.3f})"
        xth = (lambda t: float(np.interp(abs(t), T['t'], x_theta_t(T)))) if T is not None and 'pt' in T else None
        L_ += ["## C2. Designed targets in the plane", "",
               f"Part C's lattice; 'near tie': feature 0 (A = {A_B:g}) near the face between (0, {2 * H:g}) and ({H:g}, {2 * H:g}), on their "
               f"axis, its runner-up's margin {MU_C2:g} when it enters (first), two features far from it (−0.6 at (−0.5, −0.5), 0.7 at (0.5, "
               "−0.55)); 'window edge': feature 0 (0.8) at (0.96, 0) with −0.7 at (−0.4, 0.5) and 0.9 at (−0.3, −0.5). The loop to 1.3 "
               "times the path's last entry; feature 0's entrant and runner-up (the target's correlations at the start) and its nearest "
               "kernel; the runner-up's share in the loop (of the kernels started nearest feature 0) against the edge's 2-kernel inner "
               "problem p(x; t) at its gap (T), and its coefficient at the first output after the entry against the fold c(d)ℓ(d) of the "
               "same-sign features after it ((v); V's c(d) at A′/A = 0.875, interpolated in d and held at its last value beyond 7σ — there a "
               "lower estimate: c(d) grows with d), the largest; the kernels formed for "
               "feature 0 (|w| > 1% of its amplitude; by their final positions) against the rule 1 + [x < x_θ(t)], and the formed kernels "
               "that ended at another feature than the one they started nearest (captured):", "",
               "| target | ε | entrant, runner-up, nearest | margin μ, gap x = μ/ε^{1/2} | runner-up's share: loop, p(x; t) | the runner-up "
               "after the entry: \\|w\\| (τ), the fold (d) | formed for feature 0: loop, rule | captured (kernel: start's feature → "
               "end's) | formed = the path's | max\\|w − w_path\\| | rms residual |", "|---|---|---|---|---|---|---|---|---|---|---|"]
        cs2 = targets_C2()

        def c2_fold(name):
            cs_, As_ = cs2[name]
            lv = [(fold_level(float(np.linalg.norm(cs_[j] - cs_[0]) / SIG), abs(float(As_[j])), V), float(np.linalg.norm(cs_[j] - cs_[0]) / SIG))
                  for j in range(1, len(cs_)) if As_[j] * As_[0] > 0]
            lv = [q for q in lv if q[0] is not None]
            return max(lv) if lv else None
        for r in sorted(C2, key=lambda r: (r['name'], -r['eps'])):
            x = r['mu'] / np.sqrt(r['eps'])
            pr = f"{float(p_edge(T, x, r['t_off'])):.4g}" if T is not None and x <= T['x'][-1] else '—'
            rule = str(1 + int(x < xth(r['t_off']))) if xth is not None else '—'
            fo = c2_fold(r['name'])
            tail = (f"{r['w_tail']:.3g} ({r['tau_tail']:.1f}), " + (f"{fo[0]:.2g} ({fo[1]:.2f}σ)" if fo else '—') if 'w_tail' in r else '—')
            cap = (', '.join(f"{pt(j)}: {a} → {b}" for j, a, b in r['captured']) or 'none') if 'captured' in r else '—'
            gone = 'captured' in r and any(j == r['r0'] for j, a, b in r['captured'])
            L_.append(f"| {r['name']} | {r['eps']:g} | {pt(r['e0'])}, {pt(r['r0'])}, {pt(r['nearest0'])} | {r['mu']:.3e}, {x:.3g} | "
                      + ("— (captured)" if gone else f"{r['share_r']:.4g}") + f", {pr} | {tail} | {r.get('n0_fin', r['n0'])}, {rule} | {cap} | "
                      f"{'yes' if r['formed'] == r['path'] else 'no: ' + str(r['formed'])} | {r['dw']:.2e} | {r['rms']:.2e} |")
        L_.append("")
    if V is not None:
        L_ += ["## V. A runner-up held at its feature: the static law of (v)", "",
               f"The feature (A = {A_B:g}) at the origin and an unfitted one (A′ = {AP_V:g}) at the distance d (the continuum: a trapezoid of "
               "step σ/20 in one dimension, σ/8 in the plane along the diagonal); the entrant (w_e, c_e) and the runner-up's centre free, its "
               "coefficient w fixed: the local minimum's displacement Δ of the runner-up toward the unfitted feature (by continuation in w) "
               "against (v)'s law Δ² = σ²ℓ(d)·w_e/(w(w_e + w)), ℓ(d) = A′(d/σ)²e^{−d²/(4σ²)} (the pull projected off the entrant's φ, as "
               "the binding is: referee round 49); the fold — below it no minimum within 2σ — by bisection, the fold over ℓ(d), c(d), the "
               "minimum's displacement at the fold, and the law's there:", "",
               "| d/σ | Δ/σ (the law) at w = 0.3, 0.1, 0.03, 0.01, 0.003, 0.001 | the fold w_f | ℓ(d) | c(d) = w_f/ℓ | Δ/σ at w_f: the "
               "minimum's (the law's) |", "|---|---|---|---|---|---|"]
        for r in [r for r in V if r['n'] == 1]:
            cells = []
            for wq in (0.3, 0.1, 0.03, 0.01, 0.003, 0.001):
                rr = min(r['rows'], key=lambda q: abs(np.log(q['w'] / wq))) if r['rows'] else None
                cells.append(f"{rr['along']:.4f} ({v_law(r['d'], rr['w'], rr['we']):.4f})" if rr is not None and abs(np.log(rr['w'] / wq)) < 0.1
                             else '—')
            lev = AP_V * r['d'] ** 2 * np.exp(-r['d'] ** 2 / 4)
            if r['fold'] is not None:
                wf = float(np.sqrt(r['fold'][0] * r['fold'][1]))
                wef = r['rows'][-1]['we']
                fold = (f"{wf:.4g} | {lev:.3g} | {wf / lev:.2f} | " + (f"{r['along_fold']:.3f}" if r.get('along_fold') is not None else '—')
                        + f" ({v_law(r['d'], wf, wef):.3f})")
            else:
                fold = f"none above {W_V[-1]:g} | {lev:.3g} | — | —"
            L_.append(f"| {r['d']:g} | {', '.join(cells)} | {fold} |")
        gr = [(b['E'] - a['E']) / (b['w'] - a['w']) for r in toward_branches(V, V3) for a, b in zip(r['rows'][:-1], r['rows'][1:])
              if 'E' in a]
        if gr:
            L_ += ["", f"Along V's and V3's distinct one-dimensional branches (the runner-up toward A′; V3's A′/A = {AP_V / A_B:g} column "
                       f"repeats V's) the rests' energy E*(w) (the minimum over the other parameters) "
                       f"{'increases with w at every step' if min(gr) > 0 else 'does NOT increase with w at every step'} "
                       f"({len(gr)} steps; the secant dE*/dw from {min(gr):.2e} to {max(gr):.2e}): with w free the runner-up's coefficient "
                       "decays toward the fold on this side — the hold of (v) (referee round 49; the far side: VF)."]
        e3 = [(r, e3_ratios(r)) for r in list(V) + (list(VF) if VF is not None else []) if r['n'] == 1 and r.get('Ap', AP_V) == AP_V]
        if any(q for _, q in e3):
            L_ += ["", "The hold's third-moment law (`e3_law`: E*₃ = ∓(1/6)P₃(M₂*)^{3/2}(1 − w/w_e)(w(1 + w/w_e))^{−1/2}, P₃ = A′√π(d/σ)³"
                       "e^{−d²/(4σ²)}/(8σ²)): the measured secant over the law's at w ≈ 0.28, 0.1, 0.037 (the displacement there, σ) — toward "
                       "A′ (V) and away (VF):", "", "| d/σ | side | secant / the law's (displacement) |", "|---|---|---|"]
            L_ += [f"| {r['d']:g} | {'toward' if r.get('side', 1) == 1 else 'away'} | "
                   + ', '.join(f"{q:.3f} ({al:+.3f})" for _, q, al in qs) + " |" for r, qs in e3 if qs]
            two = {}
            for r, qs in e3:
                if qs and abs(qs[0][0] - 0.28) < 0.05:
                    two.setdefault(r['d'], {})[r.get('side', 1)] = qs[0][1]
            mean2 = [(d_, (v[1] + v[-1]) / 2) for d_, v in sorted(two.items()) if 1 in v and -1 in v]
            if mean2:
                L_ += ["", "The two sides' mean at w ≈ 0.28 (an even next order enters the two ratios with opposite signs): "
                       + ', '.join(f"d = {d_:g}σ: {m_:.5f}" for d_, m_ in mean2) + "."]
        L_ += vf_lines(VF)
        if V3 is not None:
            L_ += ["", "V3. The fold's constant against the unfitted feature's amplitude: c = w_f/ℓ(d) (one dimension, A = "
                   f"{A_B:g}) at A′/A = " + ', '.join(f"{a / A_B:.3g}" for a in AP_V3) + " (where the tail is unbound at w = "
                   f"{W_V[0]:g} the continuation starts at 0.99A; '>' the tail unbound from there on, the fold above it; '<' bound down to "
                   f"w = {W_V[-1]:g}; 'entrant to A′' where the entrant leaves its feature for A′ — off the branch):", "",
                   "| d/σ | " + ' | '.join(f"A′/A = {a / A_B:.3g}" for a in AP_V3) + " |", "|---|" + "---|" * len(AP_V3)]
            for d_ in D_V3:
                cells = []
                for a in AP_V3:
                    r3 = next((q for q in V3 if q['d'] == d_ and q['Ap'] == a), None)
                    lev = a * d_ ** 2 * np.exp(-d_ ** 2 / 4)
                    if r3 is None:
                        cells.append('—')
                    elif r3.get('entrant_left') and r3['fold'] is None:
                        cells.append('entrant to A′')
                    elif r3['fold'] is None:
                        cells.append(f"> {r3.get('top', W_V[0]) / lev:.3f}" if not r3['rows'] else f"< {W_V[-1] / lev:.3g}")
                    else:
                        cells.append(f"{np.sqrt(r3['fold'][0] * r3['fold'][1]) / lev:.3f}"
                                     + (' (fresh start)' if str(r3.get('start', 'cold')).startswith('fresh') else ''))
                L_.append(f"| {d_:g} | " + ' | '.join(cells) + " |")
        for r in [r for r in V if r['n'] == 2]:
            L_ += ["", f"The plane, d = {r['d']:g}σ along the diagonal: " + '; '.join(
                f"w = {q['w']:.3g}: along {q['along']:.4f}, across {q['across']:.1e} (the law {v_law(r['d'], q['w'], q['we']):.4f})"
                for q in r['rows'][::2]) + (f"; the fold between w = {r['fold'][0]:.4g} and {r['fold'][1]:.4g}" if r['fold'] else '') + "."]
        L_.append("")
    if C3 is not None:
        Kf, cf = full()
        cs3 = targets_C2()
        d_c3 = {k: float(np.linalg.norm(cs3[k][0][2] - cs3[k][0][0]) / SIG) for k in ('near tie', 'near tie 5σ')}
        L_ += ["## C3. The capture window", "",
               f"C2's near tie (feature 2, 0.7, at {d_c3['near tie']:.2f}σ from feature 0) and its variant with feature 2 at "
               f"{d_c3['near tie 5σ']:.2f}σ (−60°), at {len(EPS_C3)} masses, "
               f"the loop at the tolerance {RTOL_C3:g} (captures sit on a saddle): the runner-up's gap x, its coefficient after the entry "
               "against the fold at feature 2's distance (V: c(d)ℓ(d), interpolated in d), the kernels formed per feature by their final "
               "positions, the captured kernels, and the outcome — bound (feature 0 two kernels, the runner-up formed there), captured "
               "(feature 0 one kernel, feature 2 fitted by the migrant), a capture tie (feature 2 two kernels), stranded at the horizon (one "
               "kernel per feature, the fit not exact at the end: the tail away from both features then — it can still be recalled to its "
               "own), one kernel (the runner-up not formed, the fit exact). The runner-up's coefficient is read at the first output after "
               "the entry (on the 5σ variant the runner-up has left by then: the coefficient is the captured kernel's). At a fixed μ the "
               "gap x and the mass vary together; C4 separates them:", "",
               "| target | ε | gap x | the runner-up after the entry | the fold c(d)ℓ(d) | kernels per feature (final) | captured | rms | outcome |",
               "|---|---|---|---|---|---|---|---|---|"]
        vd = {r['d']: r for r in (V if V is not None else []) if r['n'] == 1 and r['fold'] is not None}
        fc_own = fc_fold(FC)
        for r in sorted(C3, key=lambda r: (r['name'], -r['eps'])):
            cs_, As_ = cs3[r['name']]
            d2 = float(np.linalg.norm(cs_[2] - cs_[0]) / SIG)
            if r['name'] == 'near tie' and fc_own is not None:
                foldw = f"{fc_own:.4g} (FC)"
            elif vd:
                ds_ = sorted(vd)
                cf_ = np.interp(d2, ds_, [np.sqrt(vd[q]['fold'][0] * vd[q]['fold'][1]) / (AP_V * q ** 2 * np.exp(-q ** 2 / 4)) for q in ds_])
                foldw = f"{cf_ * As_[2] * d2 ** 2 * np.exp(-d2 ** 2 / 4):.3g} (d = {d2:.2f}σ)"
            else:
                foldw = '—'
            x = r['mu'] / np.sqrt(r['eps'])
            nf = r['n_fin']
            out = c3_outcome(r)
            cap = ', '.join(f"({cf[j][0]:.3f}, {cf[j][1]:.3f}): {a} → {b}" for j, a, b in r['captured']) or 'none'
            L_.append(f"| {r['name']} | {r['eps']:.3g} | {x:.3g} | {r['w_tail']:.3g} | {foldw} | {', '.join(map(str, nf))} | {cap} | "
                      f"{r['rms']:.1e} | {out} |")
        L_.append("")
    if C4 is not None:
        L_ += ["## C4. The race: the gap held, the mass varied", "",
               f"C2's near tie with the runner-up's gap x held fixed (its margin μ = xε^{{1/2}}) at the masses {', '.join(f'{e:g}' for e in EPS_C4)} "
               f"(rtol {RTOL_C3:g}, the absolute tolerance min(10⁻¹⁴, 10⁻⁵ε): `c4_atol`): the runner-up's largest coefficient while within 2σ of feature 0 (after its entry) against the fold "
               "c(d)ℓ(d) at feature 2's distance; its departure (2σ from feature 0) in s = (τ − τ_e)/ε^{1/2} against feature 2's own entry "
               "s₂ = (τ₂ − τ_e)/ε^{1/2}, its coefficient at the departure; the final kernels per feature and the outcome; at the tolerance "
               f"{RTOL_C4B:g} (C4b) the outcome, the largest coefficient held and the departure:", "",
               "| x | ε | side (closest approach, σ) | max w held (/ fold) | departure s | s₂ | w at the departure (/ fold) | kernels per "
               f"feature (final) | outcome | rtol {RTOL_C4B:g}: side, outcome, max w held, departure s |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        b_ = {(r["x"], r["eps"]): r for r in (C4b if C4b is not None else [])}
        agree, differ = 0, []
        fc_own = fc_fold(FC)
        for r in sorted(C4, key=lambda r: (r['x'], -r['eps'])):
            fl = fc_own if fc_own is not None else fold_level(r['d2'], r['A2'], V)
            rel = lambda w: f"{w:.3g}" + (f" ({w / fl:#.3g})" if fl else '')
            rb = b_.get((r['x'], r['eps']))
            if rb is not None:
                same = c3_outcome(rb) == c3_outcome(r)
                agree += same
                if not same:
                    differ.append(f"x = {r['x']:g}, ε = {r['eps']:g}: {c3_outcome(r)} at {RTOL_C3:g} (" + side_of(r, ', ') + f"max w held "
                                  f"{r['w_max']:.4g}, the runner-up {r['d_end']:.2g}σ from feature 0 at the end), {c3_outcome(rb)} at "
                                  f"{RTOL_C4B:g} (" + side_of(rb, ', ') + f"{rb['w_max']:.4g})")
            tb = ('—' if rb is None else side_of(rb, ', ') + f"{c3_outcome(rb)}, {rb['w_max']:.3g}, "
                  + (f"{rb['s_dep']:.3g} (w {rb['w_dep']:.3g})" if rb['s_dep'] is not None else 'none'))
            sd_ = f"{r['side']} ({r['closest']:.1e})" if 'side' in r else '—'
            L_.append(f"| {r['x']:g} | {r['eps']:g} | {sd_} | {rel(r['w_max'])} | " + (f"{r['s_dep']:.3g}" if r['s_dep'] is not None else 'none')
                      + f" | {r['s2']:.3g} | " + (rel(r['w_dep']) if r['w_dep'] is not None else '—') + f" | "
                      f"{', '.join(map(str, r['n_fin']))} | {c3_outcome(r)} | {tb} |")
        if b_:
            L_ += ["", f"The two tolerances give the same outcome at {agree} of {len(C4)} points" + ("; they differ at " + '; '.join(differ)
                   if differ else '') + "."]
        L_ += ["", "The fold divided by: " + ("the configuration's own (FC), " + f"{fc_own:.4g}" if fc_own is not None else
                                                "V's one-dimensional c(d)ℓ(d) at feature 2's distance") + "."]
        L_ += c4m_lines(C4M)
        L_ += c4w_lines(C4W)
        L_.append("")
    L_ += r_lines(R, C3, C4, C4b, FC, C4M)
    L_ += rp_lines(RP)
    L_ += rc_lines(RC, RP)
    L_ += rd_lines(RD)
    if T3 is not None:
        L_ += ["## T3. The product rule at a vertex of hℤ³", "",
               f"The 8 kernels at the corners of the cube {{0, H}}³ (samples on a tensor grid of spacing {DS:g}), the feature (A = {A_B:g}) at "
               "the cube's centre — each edge's axis at the transverse offset H/√2 ((n − 1)^{1/2}H/2; H/2 only in the plane) —; the corner "
               "1_S behind the entrant by Σ_{i∈S}x_i: the 8-kernel inner problem against the product of the edges' 2-kernel shares at the "
               "offset H/√2 (the planar edge problem: in the continuum it depends on |t| only) and at H/2 (the planar offset); the last "
               "triples put a corner between x_θ(H/√2) and x_θ(H/2), where only the true offset's product gets the count:", "",
               "| gaps (x₁, x₂, x₃) | 8-kernel: formed (stop) | the product at H/√2: formed, max \\|Δ share\\|, the formed corners' relative "
               "errors | the unformed corners with shares ≥ 5·10⁻⁴: the product at H/√2's relative errors (corners) | at H/2: formed, "
               "max \\|Δ share\\|, the formed corners' relative errors |",
               "|---|---|---|---|---|"]
        for r in T3['cube']:
            sh = np.array(r['shares'])
            row, tails = [], '—'
            for pf in ((lambda v: T3['edge'][v]), ((lambda v: float(p_edge(T, v, H / 2))) if T is not None and 'pt' in T else None)):
                if pf is None:
                    row.append('—')
                    continue
                pv = [pf(v) for v in r['x']]
                prod = t3_product(pv, r['S'])
                big = sh > TH_REL
                rel = (prod[big] - sh[big]) / sh[big]
                row.append(f"{int((prod > TH_REL).sum())}, {np.abs(prod - sh).max():.2e}, {', '.join(f'{q:+.2%}' for q in rel)}")
                small = (~big) & (sh >= 5e-4)
                if not row[1:] and small.any():                      # the frozen tails, at the true offset (referee round 53)
                    rs = (prod[small] - sh[small]) / sh[small]
                    tails = f"{rs.min():+.1%} … {rs.max():+.1%} ({int(small.sum())})" if small.sum() > 1 else f"{rs[0]:+.1%} (1)"
            L_.append(f"| ({', '.join(f'{v:g}' for v in r['x'])}) | {int((sh > TH_REL).sum())} ({r['stop']}) | {row[0]} | {tails} | {row[1]} |")
        L_.append("")
    if D2 is not None and T is not None and 'pt' in T:
        xt_ = x_theta_t(T)
        L_ += ["## D2. The loop's threshold at a finite mass", "",
               f"The edge loop (Part A's, A = {A_B:g}) with the feature at the offset t along the face: the gap at which its runner-up's share "
               f"crosses {TH_REL:g} (bisection in [4.5, 8]) against the inner problem's x_θ(t) (T); the shift is the finite mass's:", "",
               "| ε | t/H | the loop's crossing | x_θ(t) | shift |", "|---|---|---|---|---|"]
        for r in sorted(D2, key=lambda r: (-r['eps'], r['t'])):
            th = float(np.interp(abs(r['t']), T['t'], xt_))
            L_.append(f"| {r['eps']:g} | {r['t'] / H:.3f} | {r['x']:.3f} | {th:.3f} | {r['x'] - th:+.3f} |")
        L_.append("")
    if W is not None:
        Kl, cl = local()
        lp = lambda i: f"({cl[i][0] / H:+.0f}, {cl[i][1] / H:+.0f})H"
        xb = W['layer'] * H / (2 * SIG ** 2 * np.sqrt(EPS_W))
        L_ += ["## W. Where the selection is the Voronoi cell", "",
               "W1. The samples' window (Part C's grid: the samples on [−1, 1]², the lattice up to ±1). The faces between neighbours on the "
               "axis x₂ = 0 (where their correlations with a feature are equal on the samples; the continuum's: the midpoints):", "",
               "| pair | the midpoint | the samples' face |", "|---|---|---|"]
        for fc in W['edge']['faces']:
            L_.append(f"| {fc['a']:g} \\| {fc['b']:g} | {0.5 * (fc['a'] + fc['b']):.4f} | "
                      + (f"{fc['face']:.4f}" if np.isfinite(fc['face']) else f"none up to {fc['b'] + 2 * H:g}") + " |")
        rows = W['edge']['rows']
        thr = {e: np.sqrt(e) * xb for e in (EPS_W, 1e-10)}
        first = {e: next((r['x'] for r in rows if r['err'] > t), np.nan) for e, t in thr.items()}
        wrong = [r['x'] for r in rows if not r['nearest']]
        L_ += ["", "The largest error of the nearest kernel's neighbours' margins against the continuum's, features on the axis at x₁: "
               + ', '.join(f"{r['x']:g}: {r['err']:.1e}" for r in rows if round(r['x'] * 100) in (0, 30, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95))
               + f"; it first exceeds ε^{{1/2}}x̄_θ ({xb:.3f}) at x₁ = " + ', '.join(f"{first[e]:g} for ε = {e:g} ({(1 - first[e]) / SIG:.2f}σ "
                                                                              "from the window's edge)" for e in thr)
               + "; the samples' entrant is not the nearest kernel at " + (f"{len(wrong)} of {len(rows)} positions, from x₁ = {min(wrong):g}"
                                                                          if wrong else "no position") + ".", "",
               f"W2. Unequal masses: the kernel at (H, 0) of mass βε, the rest ε (the local lattice, ε = {EPS_W:g}); the entrant is the kernel "
               "of the least entry time ρ_j/v_j (ρ = b^{−1/2}): a power diagram, the face between (0, 0) and (H, 0) at H/2 − σ²ln β/H in the "
               "continuum. Features (A = 0.8) on the axis at the Voronoi face H/2, at the power face p (the samples') and 2 layer widths "
               f"(w = 2σ²x̄_θε^{{1/2}}/H = {W['layer']:.2e}) to either side:", "",
               "| β | p: samples, continuum | x₁ | entrant: power, Voronoi | runner-up, its gap | formed | shares of (0, 0), (H, 0) |",
               "|---|---|---|---|---|---|---|"]
        for r in sorted(W['mass'], key=lambda r: (r['beta'], r['x1'])):
            p = W['faces'][r['beta']]
            L_.append(f"| {r['beta']:g} | {p:.5f}, {H / 2 - SIG ** 2 * np.log(r['beta']) / H:.5f} | {r['x1']:.5f} | {lp(r['power'])}, "
                      f"{lp(r['voronoi'])} | {lp(r['runner'])}, {r['gap']:.3g} | {', '.join(lp(i) for i in r['formed'])} | "
                      f"{r['shares'][0]:.4f}, {r['shares'][1]:.4f} |")
        L_ += ["", "W3. Coarse samples along x₂ (spacing ds, the grid offset by 0.3ds; along x₁ the spacing 0.05): the face between (0, 0) "
               "and (0, H) for features at (0, x₂), its shift from H/2 on the samples and from the Poisson sum (harmonics k ≤ 3 of "
               "q = e^{−π²σ²/ds²}), the margin's slope there against the continuum's H/(2σ²), and the shift in layer widths at ε = "
               f"{EPS_W:g}:", "", "| ds | ds/σ | q | shift/H: samples, Poisson | slope ratio | shift/w |", "|---|---|---|---|---|---|"]
        for r in W['alias']:
            L_.append(f"| {r['ds']:g} | {r['ds'] / SIG:.3g} | {r['q']:.3e} | {r['shift']:+.3e}, {r['shift_p']:+.3e} | {r['slope']:.4f} | "
                      f"{abs(r['shift']) * H / W['layer']:.3g} |")
        L_ += ["", "W4. Two kernels coinciding at the samples' centre of symmetry, the runner-up x_r behind, the feature (A = 0.8) δ off "
               "their centre: the inner problem's runner-up share (extrapolated in xi0) against the profile split's (iii)(b) — off the "
               "centre the kernels slide onto the feature at speeds ∝ their coefficients, apart on the way —, and the kernels' largest "
               "distance from the feature at the end:", "",
               "| x_r | δ/σ | inner problem | profile split | difference | at the end (σ) |", "|---|---|---|---|---|---|"]
        for r in sorted(W['coinc'], key=lambda r: (r['xr'], r['delta'])):
            L_.append(f"| {r['xr']:g} | {r['delta']:g} | {r['share']:.6f}{'' if r['stationary'] else ' (not stationary)'} | "
                      f"{r['profile']:.6f} | {r['share'] - r['profile']:+.2e} | {r['end']:.1e} |")
        L_ += ["", "W5. The 2-fold profile split's runner-up share p(x) against the softmax 1/(1 + e^{κx}) at the formation threshold: "
               f"p's crossing of {TH_REL:g} (x_θ of coinciding kernels at the feature's centre) and the softmax's, the softmax's relative "
               "error at p's crossing and at x = 8, and e^{−2α/a} at the runner-up's coefficient (its amplitude coordinate α; the exponential "
               "branch's correction, D‴.51(i)):", "",
               "| A | x_θ: profile split, softmax | the softmax's error: at the crossing, at x = 8 | α/a | e^{−2α/a} |", "|---|---|---|---|---|"]
        for r in W['profile']:
            L_.append(f"| {r['A']:g} | {r['x_cross']:.3f}, {r['x_soft']:.3f} | {r['err_cross']:+.3f}, {r['err8']:+.3f} | {r['alpha_a']:.3f} | "
                      f"{r['corr']:.3f} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_transit_nd.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


KEYS = ('A', 'B', 'C', 'T', 'D', 'B4', 'H', 'T2', 'W', 'C2', 'V', 'C3', 'T3', 'B4u', 'C4', 'B4x', 'V3', 'D2', 'C4b', 'R', 'VF', 'FC', 'C4M',
        'RP', 'C4W', 'RC', 'RD')                                                     # the raw's parts, in `report`'s order


def main():
    raw = os.path.join(OUT, 'transit_nd_raw.npz')
    md = os.path.join(OUT, 'transit_nd.md')
    old = E.load(raw)
    keys = KEYS
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    if res.get('C3') is not None:                                      # round 47's w_unbind, never reported (referee round 58)
        res['C3'] = [{k: v for k, v in r.items() if k != 'w_unbind'} for r in res['C3']]
    for k_ in ('C4', 'C4b'):                                           # e_t2 and w_t2, dropped from the code in round 58 (round 62)
        if res.get(k_) is not None:
            res[k_] = [{k: v for k, v in r.items() if k not in ('e_t2', 'w_t2')} for r in res[k_]]
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('T', part_T), ('D', part_D), ('B4', part_B4), ('H', part_H),
                  ('T2', part_T2), ('W', lambda: part_W(res.get('T'))), ('C2', part_C2), ('V', part_V), ('C3', part_C3),
                  ('T3', part_T3), ('B4u', lambda: part_B4u(res['B4'])), ('C4', part_C4), ('B4x', lambda: part_B4x(res['B4'])),
                  ('V3', part_V3), ('D2', part_D2), ('C4b', part_C4b), ('R', part_R), ('VF', part_VF), ('FC', part_FC),
                  ('C4M', part_C4M), ('RP', part_RP), ('C4W', part_C4W), ('RC', part_RC), ('RD', part_RD)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
