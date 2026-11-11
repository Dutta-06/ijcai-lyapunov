"""
Proposition D‴.52: close features — the transit counts them at the number-detection limit and its block places them at the support limit
====================================================================================================================================

The sliding inverse scale space (Proposition D‴.48(iii); `rich_transit.sliding_iss` with the block's least-squares polish) on two Gaussian
bumps of the kernels' width σ = 0.12 a distance Δ apart (61 samples of [−1, 1], 21 kernels started on a grid of spacing 0.1, equal masses).
With A = A₁ + A₂ ≠ 0 the pair is a monopole at the centroid c̄ = (A₁c₁ + A₂c₂)/A (kept at 0.032) with the moments V = A₁A₂Δ²/A,
T₃ = A₁A₂(A₁ − A₂)Δ³/A²: f* = Aφ(c̄) + ½V∂²_cφ(c̄) + ⅙T₃∂³_cφ(c̄) + O(Δ⁴); the first entrant fits Aφ(c̄) and leaves ½Vψ,
ψ = ∂²_cφ − (⟨∂²_cφ φ⟩/⟨φφ⟩)φ, whose correlations ½V⟨ψφ_j⟩ drive the second entry: Δτ₂ = min_j (s − u_j)/v_j over the walls s = ±1
the duals u_j(τ₁) = ⟨f*φ_j⟩/|⟨f*φ_{j₁}⟩| reach at the rates v_j = ½V⟨ψφ_j⟩ — ∝ Δ⁻²; a dipole (A = 0) has its first entry and the second
∝ Δ⁻¹.

A. The second entry against Δ (no noise): A = (0.8, 0.6), (0.8, 0.8), (0.8, −0.4), (0.8, −0.8): τ₁, the entrants, Δτ₂ against the
   leading-order prediction (for the dipole the transit of its derivative shape −A₁∂_cφ, slowed by 1/Δ), the scaled delays Δτ₂Δ² (Δτ₂Δ
   for the dipole), and the path's end (two kernels at the features?). The expansion parameter is |A₁ − A₂|Δ/(|A|σ).
B. The two-kernel block at its exact fit: the Gauss–Newton Hessian's eigenvalues and the fast flow's relaxation rates (the critical
   curve's flow linearised at the fit) against Δ — ∝ 1, 1, Δ², Δ⁶.
C. Number detection under label noise: y = f* + ξ, ξ i.i.d. N(0, s²) on the samples, A = (0.8, 0.6); under H1 the pair, under H0 the
   one-bump fit closest to it (the single-kernel least-squares fit of f*); R_C draws each; the second entry's delay Δτ₂ as the statistic:
   the power at the false-alarm rate 5% — the threshold the H0 draws' quantile at the true offset: an oracle-calibrated simple null
   (referee round 55) — against Δ for each s, beside the Neyman–Pearson bound Φ(‖f₁ − f₀‖/s − z₀.₉₅) of the same simple hypotheses; the
   separation Δ₅₀ where the power crosses ½ (with bootstrap intervals), its slope in log s (½ predicted: Δ² ∝ s).
D. Placement under label noise: the two-kernel block after the second entry on the H1 draws — the rms errors of the positions and of the
   amplitude ratio against s/Δ² and s/Δ³ (the moments' singular values Δ, Δ³) and against their Cramér–Rao bounds at the truth; how many
   blocks hold both kernels within σ of the centroid, and how many second entrants are the first's neighbours.
E. The null law and its calibration across the grid offset (referee round 55): at the centroids OFFSETS_E within one grid cell and the
   cells CELLS_E near Δ₅₀, R0_E draws under H0 and R1_E under H1; after the first fit the duals move at ⟨r₁φ_j⟩, r₁ = Rξ to first order
   (R the projection off the fit's tangent span(φ, w∂_cφ)), so the null delay is min_j (1 − sign(n_j)u_j)/|n_j|, n ~ N(0, s²ΦᵀRΦ/N²),
   over the whole dictionary (`knot_law`). Its quantile at the truth against the H0 draws'; the test calibrated per draw by a parametric
   bootstrap of that law at the draw's own first fit (ĉ, ŵ): its size and power at every offset, against the oracle-calibrated test, the
   fixed threshold of the offset CBAR, and the Neyman–Pearson bound.

Output: results/chaos/close_features.md and close_features_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D,E,S,F,U,G`).
"""
import math
import os
import sys
import time
import warnings
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('CF_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('CF_WORKERS', '3'))
CBAR = 0.032
PAIRS_A = ((0.8, 0.6), (0.8, 0.8), (0.8, -0.4), (0.8, -0.8))
DELTAS_A = (0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625)
DELTAS_B = (0.2, 0.1, 0.05, 0.025, 0.0125)
PAIR_C, SIGMAS_C, R_C, ALPHA_C = (0.8, 0.6), (0.001, 0.003, 0.01, 0.03), 200, 0.05
R_D, RHS_D = 40, 400_000          # D: the placement draws per (Δ, s) under H1 and their budget of right-hand-side evaluations each
DELTAS_C = (0.3, 0.2, 0.14, 0.1, 0.07, 0.05, 0.035, 0.025, 0.0175, 0.0125)
OFFSETS_E, CELLS_E = (0.0, 0.02, 0.032, 0.045, 0.049), ((0.0175, 0.003), (0.025, 0.01), (0.05, 0.03))
R0_E, R1_E, M_E, M_LAW = 1000, 400, 2000, 20000   # E: draws under H0 and H1, the bootstrap per draw, the law's samples at the truth
N_BOOT = 2000                                    # C: bootstrap resamples for Δ₅₀'s intervals
PAIRS_F = ((0.8, 0.6), (0.8, -0.4))               # F (referee round 59): both signs, the centroid across the whole grid cell
OFFSETS_F = (0.0, 0.02, 0.032, 0.045, 0.049, 0.055, 0.068, 0.08)
SIGMAS_F = (0.003, 0.01, 0.03)
DELTAS_F = {0.003: (0.007, 0.01, 0.014, 0.02, 0.028, 0.04), 0.01: (0.0125, 0.0175, 0.025, 0.035, 0.05, 0.07),
            0.03: (0.025, 0.035, 0.05, 0.07, 0.1, 0.14)}
R1_F, R0_F, H0_F = 200, 400, 2                   # F: H1 draws per cell; H0 draws at the separation DELTAS_F[s][H0_F] per (pair, offset, s)
G_CELLS = (((0.8, 0.6), 0.032, 0.1, 0.1), ((0.8, 0.6), 0.04, 0.06, 0.07), ((0.8, -0.4), 0.0, 0.05, 0.07),
           ((0.8, 0.6), 0.0, 0.1, 0.1), ((0.8, 0.6), 0.032, 0.06, 0.07), ((0.8, -0.4), 0.02, 0.02, 0.035))
R0_G, BAND_G, N_PERM = 3000, (1.2, 5.0), 2000    # G (referee round 82): H0 nulls (pair, c̄, s, Δ) chosen by μ₀/sd(μ̂) before they were run —
                                                  # the first three in the band BAND_G (the first two the round-82 referee's), the others
                                                  # outside it —, R0_G draws each; the permutations of the null-as-unit tests


def setup():
    warnings.filterwarnings('ignore')
    L = X.law()
    K = RT.Bumps(np.linspace(-1, 1, X.NS_C), X.SIG_C)
    return L, K, np.linspace(-1, 1, X.D_C)


def pair_centres(A1, A2, D, cbar=None):
    """The pair's centres: the centroid at cbar (default CBAR; A ≠ 0), the midpoint at cbar (a dipole)."""
    cb = CBAR if cbar is None else cbar
    A = A1 + A2
    if abs(A) < 1e-12:
        return np.array([cb - D / 2, cb + D / 2])
    return np.array([cb - A2 * D / A, cb + A1 * D / A])


def d2phi(K, c):
    c = np.atleast_1d(c)
    dx = K.xs[:, None] - c[None, :]
    return K.phi(c) * (dx ** 2 / K.sig ** 4 - 1 / K.sig ** 2)


def taylor_basis(K, c):
    """e_n = ∂ⁿ_cφ(c)/n!, n = 0 … 3, on the samples."""
    dx = K.xs - c
    ph, s2 = np.exp(-dx ** 2 / (2 * K.sig ** 2)), K.sig ** 2
    return [ph, ph * dx / s2, ph * (dx ** 2 / s2 ** 2 - 1 / s2) / 2, ph * (dx ** 3 / s2 ** 3 - 3 * dx / s2 ** 2) / 6]


def perp_norm2(vs, k):
    """⟨(P⊥e_k)²⟩, P⊥ the projection orthogonal to e_0 … e_{k−1} (sample means)."""
    B = np.array(vs[:k]).T
    e = vs[k]
    r = e - B @ np.linalg.lstsq(B, e, rcond=None)[0]
    return float(np.mean(r ** 2))


def hessian_prediction(K, A1, A2, D):
    """The block's two small Hessian eigenvalues to leading order: the split v = (1, −1, Δ/(2w₁), Δ/(2w₂)) leaves Jv = (Δ³/2)e₃,
    λ₁ = Δ⁶⟨(P⊥e₃)²⟩/8, amplitude-free; the separation (0, 0, 1/w₁, −1/w₂) leaves Jv = −2Δe₂, λ₂ = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²)."""
    cs = pair_centres(A1, A2, D)
    vs = taylor_basis(K, 0.5 * (cs[0] + cs[1]))
    return D ** 6 * perp_norm2(vs, 3) / 8, 4 * D ** 2 * perp_norm2(vs, 2) * A1 ** 2 * A2 ** 2 / (A1 ** 2 + A2 ** 2)


def predicted_delay(K, f, c0, A1, A2, D, order=2, cbar=None):
    """The leading-order second entry: the duals at the first entry u_j = ⟨fφ_j⟩/|⟨fφ_{j₁}⟩| moving at v_j = ½V⟨ψφ_j⟩ (ψ at the centroid)
    to the wall ±1 they meet first; returns (Δτ₂, j₂, the wall's sign). order = 3 (referee round 55): the residual
    P⊥(½V∂²_cφ + ⅙T₃∂³_cφ), P⊥ off the first fit's tangent span(φ, ∂_cφ) at the centroid, T₃ = A₁A₂(A₁ − A₂)Δ³/A² — the O(ϑ) term."""
    nx = K.xs.size
    A = A1 + A2
    cb = (A1 * pair_centres(A1, A2, D, cbar)[0] + A2 * pair_centres(A1, A2, D, cbar)[1]) / A
    V = A1 * A2 * D ** 2 / A
    ph, p2 = K.phi(np.array([cb]))[:, 0], d2phi(K, cb)[:, 0]
    if order == 3:
        e = taylor_basis(K, cb)                                       # e_n = ∂ⁿ_cφ/n!
        T = np.c_[e[0], e[1]]
        r = V * e[2] + A1 * A2 * (A1 - A2) * D ** 3 / A ** 2 * e[3]   # ½V∂²φ + ⅙T₃∂³φ
        psi = (r - T @ np.linalg.lstsq(T, r, rcond=None)[0]) / (0.5 * V)
    else:
        psi = p2 - (p2 @ ph) / (ph @ ph) * ph
    g = K.phi(c0).T @ f / nx
    j1 = int(np.argmax(np.abs(g)))
    u = g / abs(g[j1])
    v = 0.5 * V * (K.phi(c0).T @ psi) / nx                            # the duals' rates after the first fit (ρ = 1, per unit τ)
    best = (np.inf, None, 0)
    for j in range(c0.size):
        if j == j1 or v[j] == 0:
            continue
        wall = 1.0 if v[j] > 0 else -1.0
        t = (wall - u[j]) / v[j]
        if 0 < t < best[0]:
            best = (float(t), j, int(wall))
    return best


def a_job(args):
    A1, A2, D = args
    L, K, c0 = setup()
    cs = pair_centres(A1, A2, D)
    f = X.bumps_target(K, cs, (A1, A2))
    t0 = time.time()
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True)
    ent = [(float(t), int(j)) for t, j, k, c in R['events'] if k == 'enter']
    act = np.flatnonzero(R['act'])
    out = dict(A=(A1, A2), D=D, cs=cs.tolist(), entries=ent, n_act=int(act.size), c_end=np.sort(R['c'][act]).tolist(),
               w_end=R['w'][act][np.argsort(R['c'][act])].tolist(), rms=float(R['hist'][-1][3]), relaxed=all(R['relaxed']),
               secs=time.time() - t0)
    if abs(A1 + A2) > 1e-12:
        out['pred'] = predicted_delay(K, f, c0, A1, A2, D)
    else:                                                              # a dipole: the transit of its derivative shape, slowed by 1/Δ
        shape = -A1 * K.dphi(np.array([CBAR]))[:, 0]
        Rs = RT.sliding_iss(L, K, shape, c0, np.ones(c0.size), max_events=2, polish=True)
        es = [(float(t), int(j)) for t, j, k, c in Rs['events'] if k == 'enter']
        out['shape'] = dict(tau1=es[0][0] / D, dtau2=(es[1][0] - es[0][0]) / D, kernels=[j for t, j in es])
    return out


def part_A():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, [(a1, a2, D) for a1, a2 in PAIRS_A for D in DELTAS_A])


def b_job(args):
    A1, A2, D = args
    L, K, c0 = setup()
    nx = K.xs.size
    cs = pair_centres(A1, A2, D)
    w = np.array([A1, A2])
    J = np.c_[K.phi(cs), K.dphi(cs) * w[None, :]] / np.sqrt(nx)
    hess = np.sort(np.linalg.eigvalsh(J.T @ J))
    f = X.bumps_target(K, cs, (A1, A2))
    rhs = RT.critical_rhs(L, K, f)
    al = np.sign(w) * RT._alpha_of(L, np.abs(w), RT.coef_crit)
    z = np.r_[al, cs]
    Jf = np.zeros((4, 4))
    for k in range(4):
        h = 1e-6 * max(1.0, abs(z[k]))
        e = np.zeros(4); e[k] = h
        Jf[:, k] = (rhs(0, z + e) - rhs(0, z - e)) / (2 * h)
    ev, V = np.linalg.eig(Jf)
    rates = np.sort(-np.real(ev))
    l1, l2 = hessian_prediction(K, A1, A2, D)
    # the flow's split and separation modes, identified by their eigenvectors in z = (α₁, α₂, c₁, c₂), against their leading order:
    # ẇ = −μ∂E/∂w (μ = w′(α)·curve_mobility, the amplitudes' mobility), ċ = −Gm∂E/∂c; the split's rate λ₁·2/(1/μ₁ + 1/μ₂) (its
    # Rayleigh quotient in the mobility's metric, the centres' share O(Δ²)); the separation's at frozen amplitudes (Gm Δ² ≫ μ)
    # Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) — the projection orthogonal to the translation only.
    h_ = 1e-7
    wp = np.array([(RT.coef_crit(L, abs(a) + h_) - RT.coef_crit(L, abs(a) - h_)) / (2 * h_) for a in al])
    mu = wp * RT.curve_mobility(L, al)
    vs = taylor_basis(K, 0.5 * (cs[0] + cs[1]))
    B1 = vs[1][:, None]
    e2p = vs[2] - B1 @ np.linalg.lstsq(B1, vs[2], rcond=None)[0]
    pred_split = l1 * 2 / (1 / mu[0] + 1 / mu[1])
    pred_sep = L.G * L.m * 4 * D ** 2 * float(np.mean(e2p ** 2)) * A1 ** 2 * A2 ** 2 / (A1 ** 2 + A2 ** 2)

    def mode(s_):
        Vr = np.real(V)
        cos = np.abs(Vr.T @ s_) / (np.linalg.norm(Vr, axis=0) * np.linalg.norm(s_))
        k = int(np.argmax(cos))
        return float(-np.real(ev[k])), float(cos[k]), Vr[:, k]
    split = mode(np.r_[1 / wp[0], -1 / wp[1], D / (2 * A1), D / (2 * A2)])
    sep = mode(np.r_[0.0, 0.0, 1 / A1, -1 / A2])
    v = split[2]
    a = (v[0] * wp[0] - v[1] * wp[1]) / 2                             # the weight moved, w = ±a; the centres against Δ/(2w_k)·a
    comp = [float(v[2] / (D / (2 * A1) * a)), float(v[3] / (D / (2 * A2) * a))]
    return dict(A=(A1, A2), D=D, hess=hess.tolist(), rates=rates.tolist(), pred=(l1, l2), mu=mu.tolist(),
                split=split[:2], sep=sep[:2], pred_flow=(pred_split, pred_sep), split_c=comp)


DELTAS_SEP = (0.025, 0.0125, 0.00625, 0.003, 0.0015, 0.0005)


def sep_regimes(A1, A2, D):
    """(ii)'s separation rate in its two regimes (referee round 55): the exact linearisation −D·JᵀJ at the fit (D = diag(μ₁, μ₂, Gm, Gm)
    in (w, c), μ_k the amplitudes' mobility), its rate on the mode nearest the separation direction (0, 0, 1/w₁, −1/w₂), over the
    frozen-amplitude rate Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) and over the relaxed-amplitude one Gm·λ_sep (P⊥ off e₀ and e₁)."""
    L, K, c0 = setup()
    cs = pair_centres(A1, A2, D)
    w = np.array([A1, A2])
    J = np.c_[K.phi(cs), K.dphi(cs) * w[None, :]] / np.sqrt(K.xs.size)
    al = np.sign(w) * RT._alpha_of(L, np.abs(w), RT.coef_crit)
    h_ = 1e-7
    wp = np.array([(RT.coef_crit(L, abs(a) + h_) - RT.coef_crit(L, abs(a) - h_)) / (2 * h_) for a in al])
    d = np.r_[wp * RT.curve_mobility(L, al), L.G * L.m, L.G * L.m]
    lam, U = np.linalg.eigh(np.sqrt(d)[:, None] * (J.T @ J) * np.sqrt(d)[None, :])
    V = np.sqrt(d)[:, None] * U                                       # the eigenvectors of D·JᵀJ
    sv = np.r_[0.0, 0.0, 1 / A1, -1 / A2]
    k = int(np.argmax(np.abs(V.T @ sv) / np.linalg.norm(V, axis=0)))
    vs = taylor_basis(K, 0.5 * (cs[0] + cs[1]))
    B1 = vs[1][:, None]
    e2p = vs[2] - B1 @ np.linalg.lstsq(B1, vs[2], rcond=None)[0]
    wf = 4 * D ** 2 * A1 ** 2 * A2 ** 2 / (A1 ** 2 + A2 ** 2) * L.G * L.m
    return float(lam[k] / (wf * np.mean(e2p ** 2))), float(lam[k] / (wf * perp_norm2(vs, 2))), float(np.min(d[:2] * np.diag(J.T @ J)[:2]))


def part_B():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(a1, a2, D) for a1, a2 in PAIRS_A for D in DELTAS_B])


def one_fit(K, f):
    """The single-kernel least-squares fit of f (the closest one-bump target): returns its centre and amplitude."""
    nx = K.xs.size
    g = K.phi(np.linspace(-1, 1, 201)).T @ f / nx
    c = np.linspace(-1, 1, 201)[int(np.argmax(np.abs(g)))]
    w = 0.0
    for _ in range(100):
        ph = K.phi(np.array([c]))[:, 0]
        w = (ph @ f) / (ph @ ph)
        dp = K.dphi(np.array([c]))[:, 0]
        r = f - w * ph
        gr = w * (dp @ r)
        hs = w * w * (dp @ dp)
        c = c + gr / hs
    return float(c), float(w)


def tangent_projector(K, c, w):
    """R = I − T(TᵀT)⁻¹Tᵀ, T = [φ(c), w∂_cφ(c)]: the projection off the one-kernel fit's tangent span."""
    T = np.c_[K.phi(np.array([c]))[:, 0], w * K.dphi(np.array([c]))[:, 0]]
    return np.eye(K.xs.size) - T @ np.linalg.solve(T.T @ T, T.T)


def knot_of(K, c0, c, w, Xi):
    """The second knot of the one-bump target wφ(c) under the noise rows Xi, first order in the noise, in the sliding limit (referee round
    55): the first entrant j₁ the largest |⟨yφ_j⟩| (equal masses: the duals at its entry u_j = ⟨yφ_j⟩/|⟨yφ_{j₁}⟩|), then the duals' rates
    n_j = ⟨Rξ, φ_j⟩ (r₁ = Rξ) to their walls: Δτ₂ = min_{j ≠ j₁} (1 − sign(n_j)u_j)/|n_j|."""
    N, M = K.xs.size, Xi.shape[0]
    Phi = K.phi(c0)
    v0 = (w * K.phi(np.array([c]))[:, 0][None, :] + Xi) @ Phi / N
    j1 = np.argmax(np.abs(v0), axis=1)
    u = v0 / np.abs(v0[np.arange(M), j1])[:, None]
    n = Xi @ (tangent_projector(K, c, w) @ Phi) / N
    with np.errstate(divide='ignore'):
        t = (1 - np.sign(n) * u) / np.abs(n)
    t[np.arange(M), j1] = np.inf
    return t.min(axis=1)


def knot_law(K, c0, c, w, s, M, rng):
    """The second knot's null law at the one-bump fit (c, w): `knot_of` on M draws of i.i.d. N(0, s²) noise."""
    return knot_of(K, c0, c, w, s * rng.standard_normal((M, K.xs.size)))


def np_power(K, f1, f0, s, alpha=0.05):
    """The Neyman–Pearson power of the simple hypotheses f₀ against f₁ under i.i.d. N(0, s²) labels: Φ(‖f₁ − f₀‖/s − z_{1−α})."""
    from scipy.stats import norm
    return float(norm.cdf(np.linalg.norm(f1 - f0) / s - norm.ppf(1 - alpha)))


def e_job(args):
    """E: at the centroid cbar, R0_E draws under H0 (the one-bump fit of the pair) and R1_E under H1: the second knot Δτ₂, the entrants,
    the first fit (ĉ, ŵ) and the knot's p-value under the law bootstrapped at that fit (M_E draws); the law at the truth (M_LAW draws),
    the Neyman–Pearson power and the projected noise's std over the unprojected one per kernel, (⟨(Rφ_j)²⟩/⟨φ_j²⟩)^{1/2}."""
    D, s, cbar, seed = args
    t0 = time.time()
    L, K, c0 = setup()
    N = K.xs.size
    f1 = X.bumps_target(K, pair_centres(*PAIR_C, D, cbar), PAIR_C)
    c1, w1 = one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    rows = {}
    for hyp, f, R_ in ((0, f0, R0_E), (1, f1, R1_E)):
        out = []
        for r in range(R_):
            y = f + s * rng.standard_normal(N)
            R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
            ent = [(t, j) for t, j, k, c in R['events'] if k == 'enter']
            dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.inf
            h = R['hist'][1]
            a = int(np.flatnonzero(h[2])[0])
            law = knot_law(K, c0, float(h[1][a]), float(h[2][a]), s, M_E, np.random.default_rng([seed, hyp, r]))
            out.append(dict(dt=float(dt), j=[j for t, j in ent], c1=float(h[1][a]), w1=float(h[2][a]), p=float(np.mean(law <= dt))))
        rows[hyp] = out
    law0 = knot_law(K, c0, c1, w1, s, M_LAW, np.random.default_rng([seed, 9]))
    Rm = tangent_projector(K, c1, w1)
    Phi = K.phi(c0)
    ratio = np.sqrt(np.mean((Rm @ Phi) ** 2, axis=0) / np.mean(Phi ** 2, axis=0))
    print(f"  E Δ {D:g} s {s:g} c̄ {cbar:g}: H0 q·s {np.quantile([q['dt'] for q in rows[0]], ALPHA_C) * s:.3f} (law "
          f"{np.quantile(law0, ALPHA_C) * s:.3f}) [{time.time() - t0:.0f}s]", flush=True)
    return dict(D=D, s=s, cbar=cbar, rows0=rows[0], rows1=rows[1], q_law=float(np.quantile(law0, ALPHA_C)),
                np=np_power(K, f1, f0, s), ratio=ratio.tolist(), c1=c1, w1=w1)


def e_jobs():
    """E's jobs (Δ, s, c̄, seed), one seed per cell."""
    return [(D, s, cb, 55000 + 100 * i + k) for i, (D, s) in enumerate(CELLS_E) for k, cb in enumerate(OFFSETS_E)]


def part_E():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(e_job, e_jobs(), chunksize=1)


def c_job(args):
    """C and D: R_C draws of y = f + ξ under H0 (the one-bump fit of the pair) or H1 (the pair): the path to its second entry, not
    relaxed after it (`max_entries` = 2: the delay Δτ₂ only); under H1 the first R_D draws once more to the two-kernel block's relaxation
    within RHS_D right-hand-side evaluations — its centres and coefficients if it relaxed, 'oob' if the budget ran out first (a block
    near coalescence relaxes at a rate ∝ Δ⁶, Part B; unbounded, such relaxations exhausted the memory)."""
    D, s, hyp, seed = args
    t0 = time.time()
    L, K, c0 = setup()
    cs = pair_centres(*PAIR_C, D)
    f1 = X.bumps_target(K, cs, PAIR_C)
    if hyp == 1:
        f = f1
    else:
        c1, w1 = one_fit(K, f1)
        f = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    out = []
    for r in range(R_C):
        y = f + s * rng.standard_normal(K.xs.size)
        R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
        ent = [(t, j) for t, j, k, c in R['events'] if k == 'enter']
        dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.inf
        row = dict(dt=float(dt), j=[j for t, j in ent])
        if hyp == 1 and r < R_D and len(ent) > 1:
            R2 = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_rhs=RHS_D)
            act = np.flatnonzero(R2['act'])                          # the block after the second entry's relaxation
            row['oob'] = bool(R2['out_of_budget'])
            if not R2['out_of_budget'] and len(R2['relaxed']) > 1 and R2['relaxed'][1] and act.size == 2:
                o = np.argsort(R2['c'][act])
                row.update(c=R2['c'][act][o].tolist(), w=R2['w'][act][o].tolist())
        out.append(row)
    print(f"  C Δ {D:g} s {s:g} H{hyp}: second entries {sum(np.isfinite(q['dt']) for q in out)}/{len(out)}, placed "
          f"{sum('c' in q for q in out)}, out of budget {sum(q.get('oob', False) for q in out)} [{time.time() - t0:.0f}s]", flush=True)
    return dict(D=D, s=s, hyp=hyp, rows=out)


def c_jobs():
    """C's jobs (Δ, s, hypothesis, seed), one seed per cell and hypothesis."""
    return [(D, s, h, 1000 * k + 10 * h + i) for i, (D, s) in enumerate([(D, s) for D in DELTAS_C for s in SIGMAS_C])
            for h, k in ((0, 1), (1, 2))]


def part_C():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, c_jobs(), chunksize=1)


def knot_walls(K, c0, c, w, Xi):
    """`knot_of` split by wall (referee round 59): the first knot at the wall of the fit's sign (sign w: where a same-sign pair's
    runner-up goes) and at the other wall (an opposite-sign pair's flank kernel), first order in the noise; their minimum is `knot_of`."""
    N, M = K.xs.size, Xi.shape[0]
    Phi = K.phi(c0)
    v0 = (w * K.phi(np.array([c]))[:, 0][None, :] + Xi) @ Phi / N
    j1 = np.argmax(np.abs(v0), axis=1)
    u = v0 / np.abs(v0[np.arange(M), j1])[:, None]
    n = Xi @ (tangent_projector(K, c, w) @ Phi) / N
    with np.errstate(divide='ignore'):
        t = (1 - np.sign(n) * u) / np.abs(n)
    t[np.arange(M), j1] = np.inf
    own = np.sign(n) == np.sign(w)
    return np.where(own, t, np.inf).min(axis=1), np.where(own, np.inf, t).min(axis=1)


def walls_at_fit(K, c0, y, c, w):
    """The second knot of the labels y per wall from the path's own first fit (c, w), exact in the sliding limit: the duals at the first
    entry u_j = ⟨yφ_j⟩/|⟨yφ_{j₁}⟩| move at ⟨r₁φ_j⟩, r₁ = y − wφ(c), to the wall sign⟨r₁φ_j⟩ — the minimum at the wall of the fit's sign and
    at the other (their minimum the path's Δτ₂)."""
    N = K.xs.size
    Phi = K.phi(c0)
    g = Phi.T @ y / N
    j1 = int(np.argmax(np.abs(g)))
    u = g / abs(g[j1])
    v = Phi.T @ (y - w * K.phi(np.array([c]))[:, 0]) / N
    with np.errstate(divide='ignore'):
        t = (1 - np.sign(v) * u) / np.abs(v)
    t[j1] = np.inf
    own = np.sign(v) == np.sign(w)
    return float(np.where(own, t, np.inf).min()), float(np.where(own, np.inf, t).min())


def score_stats(K, y, c, w, s):
    """The score statistic of one bump against two at the one-bump fit (c, w) (referee round 59): z = ⟨r, ψ̂⟩/(s‖ψ̂‖), r = y − wφ(c),
    ψ̂ = R∂²_cφ(c) (R off the fit's tangent span(φ, w∂_cφ); at a fit ⟨r, ψ̂⟩ = ⟨r, ∂²_cφ⟩) — N(0, 1) under H0 to first order, mean
    ½V‖ψ̂‖/s under H1 —, and its t-version with ŝ² = ‖r‖²/(N − 2)."""
    N = K.xs.size
    r = y - w * K.phi(np.array([c]))[:, 0]
    psi = tangent_projector(K, c, w) @ d2phi(K, c)[:, 0]
    q = float(psi @ r / np.linalg.norm(psi))
    return q / s, q / (np.linalg.norm(r) / np.sqrt(N - 2))


def width_fit(K, y, c, w):
    """The one-bump least-squares fit with a free width (a moving-width kernel, D‴.59) from (w, c, σ) (Levenberg–Marquardt): (σ̂, the
    standard deviation of σ̂ per unit noise, [(JᵀJ)⁻¹]_σσ^{1/2} at the fit)."""
    from scipy.optimize import least_squares
    x, s0 = K.xs, K.sig

    def jac(th):
        dx = x - th[1]
        g = np.exp(-dx ** 2 / (2 * th[2] ** 2))
        return np.c_[g, th[0] * g * dx / th[2] ** 2, th[0] * g * dx ** 2 / th[2] ** 3]

    res = least_squares(lambda th: th[0] * np.exp(-(x - th[1]) ** 2 / (2 * th[2] ** 2)) - y, np.array([w, c, s0]), jac=jac,
                        method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=2000)
    J = jac(res.x)
    try:
        return float(res.x[2]), float(np.sqrt(np.linalg.inv(J.T @ J)[2, 2]))
    except np.linalg.LinAlgError:                                  # a width run off to 0: no information
        return float(res.x[2]), np.inf


def width_z(K, y, c, w, s):
    """z_σ = (σ̂ − σ)/sd(σ̂) from `width_fit`, sd(σ̂) at the noise s — to first order sign(w)·z (∂_σφ = σ∂²_cφ + φ/σ); ±inf if the width
    leaves (σ/5, 5σ)."""
    s0 = K.sig
    sg, sd = width_fit(K, y, c, w)
    if not s0 / 5 < sg < 5 * s0:
        return float(np.sign(sg - s0) * np.inf)
    return float((sg - s0) / (s * sd))


def draw_stats(K, c0, y, c, w, s):
    """Per draw at the one-bump fit (c, w): the score z, its t-version, the free-width z_σ and the knot per wall (`walls_at_fit`)."""
    z, t = score_stats(K, y, c, w, s)
    tA, tB = walls_at_fit(K, c0, y, c, w)
    return dict(z=z, t=t, zs=width_z(K, y, c, w, s), tA=tA, tB=tB)


def s_job(args):
    """S (referee round 59): C's and E's draws regenerated from their seeds (the jobs' rng draws only the labels, y = f + sξ in order; E's
    H0 draws first, then H1's): per draw `draw_stats` at the one-bump fit — E's path's own first fit (stored; `one_fit` of the labels is
    checked against it: 'dfit'), C's `one_fit` —; for E also the runner-up's filter z₂ = wall·⟨r, Rφ_{j₂}⟩/(s‖Rφ_{j₂}‖) (j₂ the noiseless
    second entrant, its wall), whether the first entrant differs from the noiseless one, the runner-up's margin 1 − wall·u_{j₂} (noiseless
    'mu_true'), the cosine of Rφ_{j₂} with f₁ − f₀; for a C cell under H0 the knot's law at the truth (`knot_walls`, M_LAW samples): the α
    quantile of the minimum and the α/2 quantile per wall."""
    kind, D, s, cbar, hyp, seed, ns, fits = args
    L, K, c0 = setup()
    N = K.xs.size
    f1 = X.bumps_target(K, pair_centres(*PAIR_C, D, cbar), PAIR_C)
    c1, w1 = one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    xi = [rng.standard_normal(N) for _ in range(sum(ns))]
    hyps = (hyp,) if kind == 'C' else (0, 1)
    pd = predicted_delay(K, f1, c0, *PAIR_C, D, cbar=cbar)
    j1n = int(np.argmax(np.abs(K.phi(c0).T @ f1)))
    out = dict(kind=kind, D=D, s=s, cbar=cbar, hyp=hyp)
    k = 0
    for h, n in zip(hyps, ns):
        rows = []
        for r in range(n):
            y = (f1 if h == 1 else f0) + s * xi[k]
            cf, wf = one_fit(K, y)
            row = dict(dfit=0.0)
            if fits is not None:
                row['dfit'] = abs(cf - fits[k][0])
                cf, wf = fits[k]
            row.update(draw_stats(K, c0, y, cf, wf, s))
            if kind == 'E':
                Rp = tangent_projector(K, cf, wf) @ K.phi(c0[[pd[1]]])[:, 0]
                g = K.phi(c0).T @ y
                row.update(z2=float(pd[2] * (y - wf * K.phi(np.array([cf]))[:, 0]) @ Rp / (s * np.linalg.norm(Rp))),
                           swap=bool(int(np.argmax(np.abs(g))) != j1n), mu=float(1 - pd[2] * g[pd[1]] / np.abs(g).max()),
                           mu_u=float(1 - pd[2] * g[pd[1]] / abs(g[j1n])))
            rows.append(row)
            k += 1
        out[f'rows{h}'] = rows
    if kind == 'C' and hyp == 0:
        la, lb = knot_walls(K, c0, c1, w1, s * np.random.default_rng([seed, 9]).standard_normal((M_LAW, N)))
        out.update(q_law=float(np.quantile(np.minimum(la, lb), ALPHA_C)), qA=float(np.quantile(la, ALPHA_C / 2)),
                   qB=float(np.quantile(lb, ALPHA_C / 2)))
    if kind == 'E':
        g = K.phi(c0).T @ f1
        Rp = tangent_projector(K, c1, w1) @ K.phi(c0[[pd[1]]])[:, 0]
        out.update(j2=pd[1], wall2=pd[2], mu_true=float(1 - pd[2] * g[pd[1]] / np.abs(g).max()),
                   cos2=float(Rp @ (f1 - f0) / (np.linalg.norm(Rp) * np.linalg.norm(f1 - f0))))
    zz = [q['z'] for h in hyps for q in out[f'rows{h}']]
    print(f"  S {kind} Δ {D:g} s {s:g} c̄ {cbar:g} H{hyp}: z mean {np.mean(zz):+.3f}, sd {np.std(zz):.3f}", flush=True)
    return out


def part_S(C, E_):
    """S's jobs: C's (seeds as `part_C`'s) and E's (seeds as `part_E`'s; the path's first fits from E's rows)."""
    from multiprocessing import Pool
    cseed = {(D, s, h): sd for D, s, h, sd in c_jobs()}
    eseed = {(D, s, cb): sd for D, s, cb, sd in e_jobs()}
    jobs = [('C', r['D'], r['s'], CBAR, r['hyp'], cseed[(r['D'], r['s'], r['hyp'])], (len(r['rows']),), None) for r in C]
    for e in E_:
        fits = [(q['c1'], q['w1']) for q in e['rows0']] + [(q['c1'], q['w1']) for q in e['rows1']]
        jobs.append(('E', e['D'], e['s'], e['cbar'], 2, eseed[(e['D'], e['s'], e['cbar'])], (len(e['rows0']), len(e['rows1'])), fits))
    with Pool(N_WORKERS) as pool:
        return pool.map(s_job, jobs, chunksize=1)


def f_seed(A, cbar, s, D, hyp):
    return 700000 + 10000 * PAIRS_F.index(A) + 1000 * OFFSETS_F.index(cbar) + 100 * SIGMAS_F.index(s) + 10 * DELTAS_F[s].index(D) + hyp


def f_job(args):
    """F (referee round 59): the pair A at the centroid cbar, noise s, separation D; R_ draws under hyp (H1 the pair, H0 its one-bump fit):
    per draw the path's second knot Δτ₂ (to the second entry) and entrants, its first fit (ĉ, ŵ), `draw_stats` there, and the knot's
    bootstrap p-values at that fit (M_E samples of `knot_walls`): of Δτ₂ against the minimum, of each wall's knot against its own law;
    per cell the law at the truth (M_LAW samples: the α quantile of the minimum, the α/2 quantile per wall) and the Neyman–Pearson power."""
    A, cbar, s, D, hyp, R_ = args[:6]
    seed = args[6] if len(args) > 6 else f_seed(A, cbar, s, D, hyp)            # G's cells carry their own seed (referee round 82)
    t0 = time.time()
    L, K, c0 = setup()
    N = K.xs.size
    f1 = X.bumps_target(K, pair_centres(*A, D, cbar), A)
    c1, w1 = one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    rows = []
    for r in range(R_):
        y = (f1 if hyp == 1 else f0) + s * rng.standard_normal(N)
        R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), max_events=2, polish=True, max_entries=2)
        ent = [(t, j) for t, j, k, c in R['events'] if k == 'enter']
        dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.inf
        h = R['hist'][1]
        a = int(np.flatnonzero(h[2])[0])
        ch, wh = float(h[1][a]), float(h[2][a])
        row = dict(dt=float(dt), j=[j for t, j in ent], c1=ch, w1=wh)
        row.update(draw_stats(K, c0, y, ch, wh, s))
        la, lb = knot_walls(K, c0, ch, wh, s * np.random.default_rng([seed, 1, r]).standard_normal((M_E, N)))
        row.update(p=float(np.mean(np.minimum(la, lb) <= dt)), pA=float(np.mean(la <= row['tA'])), pB=float(np.mean(lb <= row['tB'])))
        rows.append(row)
    la, lb = knot_walls(K, c0, c1, w1, s * np.random.default_rng([seed, 0]).standard_normal((M_LAW, N)))
    out = dict(A=A, cbar=cbar, s=s, D=D, hyp=hyp, rows=rows, c1=c1, w1=w1, np=np_power(K, f1, f0, s, ALPHA_C),
               q_law=float(np.quantile(np.minimum(la, lb), ALPHA_C)), qA=float(np.quantile(la, ALPHA_C / 2)),
               qB=float(np.quantile(lb, ALPHA_C / 2)))
    print(f"  F {A} c̄ {cbar:g} s {s:g} Δ {D:g} H{hyp}: knot < law {np.mean([q['dt'] < out['q_law'] for q in rows]):.3f}, |z| "
          f"{np.mean([abs(q['z']) > 1.959964 for q in rows]):.3f}, NP {out['np']:.3f} [{time.time() - t0:.0f}s]", flush=True)
    return out


def f_jobs():
    """F's jobs: the H0 draws at the third separation of each (pair, offset, s), then the H1 cells."""
    return ([(A, cb, s, DELTAS_F[s][H0_F], 0, R0_F) for A in PAIRS_F for cb in OFFSETS_F for s in SIGMAS_F]
            + [(A, cb, s, D, 1, R1_F) for A in PAIRS_F for cb in OFFSETS_F for s in SIGMAS_F for D in DELTAS_F[s]])


def part_F():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(f_job, f_jobs(), chunksize=1)


def g_jobs():
    """G's jobs (referee round 82): the H0 null of each G_CELLS pair at its offset, noise and separation, R0_G draws (`f_job`), seeds
    820000 + 10i."""
    return [(A, cb, s, D, 0, R0_G, 820000 + 10 * i) for i, (A, cb, s, D) in enumerate(G_CELLS)]


def part_G():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(f_job, g_jobs(), chunksize=1)


def crossing(Ds, pw, level=0.5):
    """The separation where the power, increasing in Δ, first crosses `level` (log-linear interpolation), NaN if it does not."""
    for k in range(len(Ds) - 1):
        if pw[k] < level <= pw[k + 1]:
            t = (level - pw[k]) / (pw[k + 1] - pw[k])
            return float(np.exp(np.log(Ds[k]) + t * (np.log(Ds[k + 1]) - np.log(Ds[k]))))
    return np.nan


def detection(C):
    """Per (Δ, s): the H0 quantile of Δτ₂ at ALPHA_C and the H1 power; Δ₅₀(s) by log-linear interpolation of the power."""
    tab = {}
    for r in C:
        tab.setdefault((r['D'], r['s']), {})[r['hyp']] = np.array([q['dt'] for q in r['rows']])
    out = {}
    for (D, s), d in tab.items():
        q = np.quantile(d[0], ALPHA_C)
        out[(D, s)] = dict(q0=float(q), power=float(np.mean(d[1] < q)), med1=float(np.median(d[1])), med0=float(np.median(d[0])))
    Ds = sorted(DELTAS_C)
    d50 = {s: crossing(Ds, [out[(D, s)]['power'] for D in Ds]) for s in SIGMAS_C}
    return out, d50


def d50_bootstrap(C, n=N_BOOT, seed=5500):
    """Δ₅₀'s and its slope's 95% intervals: both hypotheses' draws resampled per cell (the H0 quantile recomputed each time)."""
    tab = {}
    for r in C:
        tab.setdefault((r['D'], r['s']), {})[r['hyp']] = np.array([q['dt'] for q in r['rows']])
    rng = np.random.default_rng(seed)
    Ds = sorted(DELTAS_C)
    pw = {}
    for (D, s), d in sorted(tab.items()):
        q = np.quantile(d[0][rng.integers(0, d[0].size, (n, d[0].size))], ALPHA_C, axis=1)
        pw[(D, s)] = np.mean(d[1][rng.integers(0, d[1].size, (n, d[1].size))] < q[:, None], axis=1)
    d50 = {s: np.array([crossing(Ds, [pw[(D, s)][b] for D in Ds]) for b in range(n)]) for s in SIGMAS_C}
    sl = []
    for b in range(n):
        x = [(np.log(s), np.log(d50[s][b])) for s in SIGMAS_C if np.isfinite(d50[s][b])]
        sl.append(np.polyfit(*np.array(x).T, 1)[0] if len(x) >= 2 else np.nan)
    ci = lambda v: ((float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))) if np.isfinite(v).any() else (np.nan, np.nan)) + \
        (int(np.sum(~np.isfinite(v))),)
    return {s: ci(d50[s]) for s in SIGMAS_C}, ci(np.array(sl))


def crb(K, cs, w, s):
    """The Cramér–Rao bounds of the two-kernel fit at the truth, s²(JᵀJ)⁻¹ (J the samples' Jacobian in (w₁, w₂, c₁, c₂)): the positions'
    rms ((Var c₁ + Var c₂)/2)^{1/2} and the split w₁/(w₁ + w₂)'s std."""
    J = np.c_[K.phi(cs), K.dphi(cs) * w[None, :]]
    Cov = s ** 2 * np.linalg.inv(J.T @ J)
    W = w.sum()
    g = np.array([w[1] / W ** 2, -w[0] / W ** 2, 0.0, 0.0])
    return float(np.sqrt((Cov[2, 2] + Cov[3, 3]) / 2)), float(np.sqrt(g @ Cov @ g))


def crb_sep(K, cs, w, s, known=False):
    """The Cramér–Rao bound of the separation ĉ₂ − ĉ₁ (referee round 64): its standard deviation from s²(JᵀJ)⁻¹, the amplitudes
    unknown (J in (w₁, w₂, c₁, c₂)) or known (J in (c₁, c₂))."""
    J = K.dphi(cs) * w[None, :] if known else np.c_[K.phi(cs), K.dphi(cs) * w[None, :]]
    Cov = s ** 2 * np.linalg.inv(J.T @ J)
    g = np.array([-1.0, 1.0]) if known else np.array([0.0, 0.0, -1.0, 1.0])
    return float(np.sqrt(g @ Cov @ g))


SMITH_PAIRS, SMITH_D = ((0.8, 0.6), (0.8, -0.4), (0.8, 0.8)), (0.0125, 0.00625)


def smith_line(K):
    """Smith's criterion (referee round 64): the separation's bound sd(ĉ₂ − ĉ₁) = c·s/Δ^p (p the local slope between SMITH_D), the
    amplitudes unknown or known, and the noise s = Δ^{p+1}/c at which sd = Δ."""
    cells = []
    for A in SMITH_PAIRS:
        for known in (False, True):
            sd = [crb_sep(K, pair_centres(*A, D), np.array(A), 1.0, known) for D in SMITH_D]
            p = float(np.log(sd[1] / sd[0]) / np.log(SMITH_D[0] / SMITH_D[1]))
            c = sd[0] * SMITH_D[0] ** p
            pi = round(p)                                               # referee round 67: the constant at the integer power too
            ci = [q * D ** pi for q, D in zip(sd, SMITH_D)]
            cells.append(f"({A[0]:g}, {A[1]:g}), {'known' if known else 'unknown'}: {c:.3g}·s/Δ^{p:.2f}, sd = Δ at s ≈ {1 / c:.3g}Δ^{p + 1:.2f} "
                         f"(with sd ∝ s/Δ^{pi} exactly: s ≈ {1 / ci[0]:.3g}Δ^{pi + 1}, {1 / ci[1]:.3g}Δ^{pi + 1} from Δ = {SMITH_D[0]:g}, "
                         f"{SMITH_D[1]:g})")
    return ("Smith's criterion (the separation's Cramér–Rao standard deviation equal to Δ; referee round 64), the amplitudes unknown or known: "
            + '; '.join(cells) + '.')


Z95, Z975 = 1.6448536269514722, 1.959963984540054   # the standard normal's 95% and 97.5% quantiles
Z99, Z995 = 2.3263478740408408, 2.5758293035489004  # and 99%, 99.5% (referee round 64: the sizes at α = 0.01)


def cp_int(k, n, a=0.05):
    """The Clopper–Pearson interval of the proportion k/n."""
    from scipy.stats import beta
    lo = float(beta.ppf(a / 2, k, n - k + 1)) if k > 0 else 0.0
    hi = float(beta.ppf(1 - a / 2, k + 1, n - k)) if k < n else 1.0
    return lo, hi


def d50_text(Ds, pw):
    """Δ₅₀ from the powers pw at the separations Ds (increasing): '< Δ_min' if the power is ½ already there, '> Δ_max' if it never is."""
    if pw[0] >= 0.5:
        return f"< {Ds[0]:g}"
    x = crossing(list(Ds), list(pw))
    return f"{x:.4g}" if np.isfinite(x) else f"> {Ds[-1]:g}"


def d50_boot(Ds, rej, n=N_BOOT, seed=5600):
    """Δ₅₀ by `crossing` and its bootstrap 95% interval, the draws (rejections per separation) resampled per separation."""
    rng = np.random.default_rng(seed)
    pt = crossing(list(Ds), [float(np.mean(rej[D])) for D in Ds])
    bs = np.array([crossing(list(Ds), [float(np.mean(rej[D][rng.integers(0, rej[D].size, rej[D].size)])) for D in Ds])
                   for _ in range(n)])
    if not np.isfinite(bs).any():
        return pt, (np.nan, np.nan)
    return pt, (float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5)))


def rejections(rows, sgn, q_law=None, qA=None, qB=None, dts=None, alloc=None):
    """Per draw (referee round 59), at the level ALPHA_C: the knot against the law at the truth (q_law), the calibrated knot (its bootstrap
    p, if recorded), both walls against their laws (qA, qB at α/2 each), the score z one-sided toward the pair's sign sgn = sign(A₁A₂)
    (z times the fit's sign), |z| and the width's |z_σ| two-sided; referee round 64: each wall's knot one-sided, calibrated at the fit
    (p_A < α: the wall of the fit's sign, a same-sign pair's; p_B < α: the other, an opposite-sign pair's). dts: the knots, if not in the
    rows (C's draws, scored by S). alloc (referee round 67): the score spending α as the calibrated knot does — the share alloc of α
    toward the fit's sign (ŵz > z_{1−α·alloc}), the rest away from it."""
    z = np.array([q['z'] for q in rows])
    zw = np.array([q['zs'] for q in rows])
    w = np.array([np.sign(q['w1']) if 'w1' in q else 1.0 for q in rows])
    out = dict(z=sgn * w * z > Z95, z2=np.abs(z) > Z975, zs=np.abs(zw) > Z975)
    if q_law is not None:
        out['knot'] = (np.array([q['dt'] for q in rows]) if dts is None else np.asarray(dts, float)) < q_law
    if qA is not None:
        out['walls'] = (np.array([q['tA'] for q in rows]) < qA) | (np.array([q['tB'] for q in rows]) < qB)
    if rows and 'p' in rows[0]:
        out['cal'] = np.array([q['p'] for q in rows]) < ALPHA_C
    if rows and 'pA' in rows[0]:                                       # referee round 64: each wall one-sided, calibrated at the fit
        out['own'] = np.array([q['pA'] for q in rows]) < ALPHA_C
        out['other'] = np.array([q['pB'] for q in rows]) < ALPHA_C
    if alloc is not None:
        from scipy.stats import norm
        out['alloc'] = (w * z > norm.ppf(1 - ALPHA_C * alloc)) | (-w * z > norm.ppf(1 - ALPHA_C * (1 - alloc)))
    return out


def c_prime_lines(C, S, npw):
    """C′ (referee round 59): C's draws with the knot calibrated by its law at the truth (M_LAW samples, `knot_walls`) instead of the 200
    H0 draws' quantile, and the score z at the one-bump fit (S), against the Neyman–Pearson bound."""
    Sc = {(r['D'], r['s'], r['hyp']): r for r in S if r['kind'] == 'C'}
    dts = {(r['D'], r['s'], r['hyp']): [q['dt'] for q in r['rows']] for r in C}
    Ds = sorted(DELTAS_C)
    out = ['### C′. The same draws: the knot calibrated by its law at the truth, and the score', '',
           f"Per (Δ, s): the knot's power at its law's {ALPHA_C:g} quantile at the truth (`knot_walls`, {M_LAW} samples) with that "
           "threshold's false-alarm rate on C's H0 draws in brackets; the score z = ⟨r, R∂²_cφ⟩/(s‖R∂²_cφ‖) at the draw's one-bump fit "
           f"(one-sided at z₀.₉₅; its false-alarm rate in brackets); the Neyman–Pearson power:", '',
           '| Δ | ' + ' | '.join(f's = {s:g}: knot (law), z, NP' for s in SIGMAS_C) + ' |', '|---|' + '---|' * len(SIGMAS_C)]
    rk, rz, r2 = {}, {}, {}
    for D in Ds:
        cells = []
        for s_ in SIGMAS_C:
            h0, h1 = Sc[(D, s_, 0)], Sc[(D, s_, 1)]
            a0 = rejections(h0['rows0'], 1.0, h0['q_law'], dts=dts[(D, s_, 0)])
            a1 = rejections(h1['rows1'], 1.0, h0['q_law'], dts=dts[(D, s_, 1)])
            rk[(D, s_)], rz[(D, s_)], r2[(D, s_)] = a1['knot'], a1['z'], a1['z2']
            cells.append(f"{np.mean(a1['knot']):.3f} ({np.mean(a0['knot']):.3f}), {np.mean(a1['z']):.3f} ({np.mean(a0['z']):.3f}), "
                         f"{npw[(D, s_)]:.3f}")
        out.append(f"| {D:g} | " + ' | '.join(cells) + ' |')
    sizes_k = np.concatenate([rejections(Sc[(D, s_, 0)]['rows0'], 1.0, Sc[(D, s_, 0)]['q_law'], dts=dts[(D, s_, 0)])['knot']
                              for D in Ds for s_ in SIGMAS_C])
    sizes_z = np.concatenate([rejections(Sc[(D, s_, 0)]['rows0'], 1.0)['z'] for D in Ds for s_ in SIGMAS_C])
    d50k = {s_: d50_boot(Ds, {D: rk[(D, s_)] for D in Ds}) for s_ in SIGMAS_C}
    d50z = {s_: d50_boot(Ds, {D: rz[(D, s_)] for D in Ds}) for s_ in SIGMAS_C}
    d502 = {s_: crossing(Ds, [float(np.mean(r2[(D, s_)])) for D in Ds]) for s_ in SIGMAS_C}
    npd = {s_: crossing(Ds, [npw[(D, s_)] for D in Ds]) for s_ in SIGMAS_C}
    ck, cz = cp_int(int(sizes_k.sum()), sizes_k.size), cp_int(int(sizes_z.sum()), sizes_z.size)
    out += ['', "Δ₅₀ (bootstrap 95% intervals, the H1 draws resampled): the knot calibrated by its law " +
            ', '.join(f"s = {s:g}: {fmt_(d50k[s][0])} [{fmt_(d50k[s][1][0])}, {fmt_(d50k[s][1][1])}]" for s in SIGMAS_C) +
            "; the score z " + ', '.join(f"s = {s:g}: {fmt_(d50z[s][0])} [{fmt_(d50z[s][1][0])}, {fmt_(d50z[s][1][1])}]" for s in SIGMAS_C) +
            "; |z| (two-sided) " + ', '.join(f"{fmt_(d502[s])}" for s in SIGMAS_C) + "; the Neyman–Pearson bound " +
            ', '.join(f"{fmt_(npd[s])}" for s in SIGMAS_C) + f". Pooled over C's {sizes_k.size} H0 draws, the law's threshold's false-alarm "
            f"rate {np.mean(sizes_k):.4f} [{ck[0]:.4f}, {ck[1]:.4f}] and the score's {np.mean(sizes_z):.4f} [{cz[0]:.4f}, {cz[1]:.4f}] "
            f"(nominal {ALPHA_C:g}; Clopper–Pearson). The knots at the one-bump fit (`walls_at_fit`) reproduce C's Δτ₂ to "
            f"{max(np.nanmax(np.abs(np.array([min(q['tA'], q['tB']) for q in Sc[k_]['rows' + str(k_[2])]]) / np.array(dts[k_]) - 1)[np.isfinite(dts[k_])]) for k_ in Sc):.1e} "
            "(relative).", '']
    return out


def fmt_(x):
    return '—' if not np.isfinite(x) else f"{x:.4g}"


def e_prime_lines(E_, S, U=None):
    """E′ (referee round 59): E's draws scored at the path's own first fit (S): the score z (one-sided), |z|, its t-version with
    ŝ² = ‖r‖²/(N − 2), the free-width z_σ, the runner-up's filter z₂, against the knot (oracle, calibrated) and the Neyman–Pearson bound;
    the knot's loss mechanisms: the first entrant swapped under H0, the runner-up's margin and its noise."""
    Se = {(r['D'], r['s'], r['cbar']): r for r in S if r['kind'] == 'E'}
    out = ['### E′. The same draws: the score at the first fit, and what the knot adds to the runner-up\'s filter', '',
           "Per cell: the score z = ⟨r̂₁, R̂∂²_cφ(ĉ)⟩/(s‖R̂∂²_cφ(ĉ)‖) at the path's own first fit (ĉ, ŵ) — N(0, 1) under H0 to first order —: "
           "its false-alarm rate (Clopper–Pearson 95%) and power one-sided at z₀.₉₅, the power of |z| at z₀.₉₇₅, the false-alarm rate of its "
           "t-version (ŝ² = ‖r̂₁‖²/(N − 2)) and of the free-width fit's z_σ = (σ̂ − σ)/sd(σ̂) with its power; the runner-up's filter "
           "z₂ = ⟨r̂₁, R̂φ_{j₂}⟩/(s‖R̂φ_{j₂}‖) (j₂ the noiseless second entrant, toward its wall; one-sided) with the cosine of R̂φ_{j₂} with "
           "f₁ − f₀; the knot's powers (oracle, calibrated) and, against its law at the truth (as the runner-up's own knot), its size and "
           "power; the Neyman–Pearson bound; under H0 the share of draws whose first entrant is not the noiseless one, and the runner-up's "
           "margin μ = 1 − wall·u_{j₂} noiseless at the null f₀ (referee round 71: the null's, as its standard deviation is the H0 "
           "draws'; H1's before; over f₀'s first entrant's |⟨f₀φ_{j}⟩|), its standard deviation over s — uncensored (referee round 67): "
           "u_{j₂} = ⟨yφ_{j₂}⟩ over |⟨yφ_{j₁}⟩|, j₁ f₁'s first entrant, not over the largest, which censors μ̂ at 0 where j₂ leads — and "
           "the margin over that standard deviation (the two first entrants coincide where U records them: the summary):", '',
           '| Δ, s | c̄ | z: size [CP], power | \\|z\\| power | t size | z_σ: size, power | z₂: size, power (cos) | knot: oracle, calibrated; '
           "against its law at the truth: size, power | the runner-up's own knot: size, power | Neyman–Pearson | first entrant swapped "
           "(H0) | margin μ, sd(μ̂)/s, μ/sd(μ̂) |",
           '|---|---|---|---|---|---|---|---|---|---|---|---|']
    Ue = {(u['D'], u['s'], u['cbar']): u for u in (U if U is not None else [])}
    summ = []
    for D, s_ in CELLS_E:
        for e in sorted([e for e in E_ if e['D'] == D and e['s'] == s_], key=lambda e: e['cbar']):
            r = Se[(D, s_, e['cbar'])]
            a0, a1 = rejections(r['rows0'], 1.0), rejections(r['rows1'], 1.0)
            t0 = np.array([q['t'] for q in r['rows0']]) > Z95
            z20 = np.array([q['z2'] for q in r['rows1']]) > Z95
            z2s = float(np.mean(np.array([q['z2'] for q in r['rows0']]) > Z95))
            d0 = np.array([q['dt'] for q in e['rows0']]); d1 = np.array([q['dt'] for q in e['rows1']])
            p1 = np.array([q['p'] for q in e['rows1']])
            po, pc = float(np.mean(d1 < np.quantile(d0, ALPHA_C))), float(np.mean(p1 < ALPHA_C))
            k0 = int(a0['z'].sum()); ci = cp_int(k0, a0['z'].size)
            mu = np.array([q['mu_u'] for q in r['rows0']])
            uq = Ue.get((D, s_, e['cbar']))
            mu_n = uq['mu0'] if uq is not None and 'mu0' in uq else None           # the null's margin from U (— without it)
            kl0, kl1 = float(np.mean(d0 < e['q_law'])), float(np.mean(d1 < e['q_law']))
            summ.append(dict(np=e['np'], z=float(np.mean(a1['z'])), z2=float(np.mean(z20)), knot=max(po, pc), cal=pc, k0=k0, n0=a0['z'].size,
                             t0=int(t0.sum()), zs0=int(a0['zs'].sum()), zs1=float(np.mean(np.array([q['zs'] for q in r['rows1']]) > Z95))))
            out.append(f"| {D:g}, {s_:g} | {e['cbar']:g} | {np.mean(a0['z']):.3f} [{ci[0]:.3f}, {ci[1]:.3f}], {np.mean(a1['z']):.3f} | "
                       f"{np.mean(a1['z2']):.3f} | {np.mean(t0):.3f} | "
                       f"{np.mean(np.array([q['zs'] for q in r['rows0']]) > Z95):.3f}, {summ[-1]['zs1']:.3f} | {z2s:.3f}, {np.mean(z20):.3f} "
                       f"({r['cos2']:.3f}) | {po:.3f}, {pc:.3f}; {kl0:.3f}, {kl1:.3f} | {own_knot(Ue.get((D, s_, e['cbar'])))} | {e['np']:.3f} | "
                       f"{np.mean([q['swap'] for q in r['rows0']]):.3f} | "
                       + (f"—, {np.std(mu) / s_:.3f}, — |" if mu_n is None else f"{mu_n:.4f}, {np.std(mu) / s_:.3f}, {mu_n / np.std(mu):.2f} |"))
            u = Ue.get((D, s_, e['cbar']))
            summ[-1].update(law1=kl1, cb=e['cbar'], D=D, s=s_, swap=float(np.mean([q['swap'] for q in r['rows0']])))
            if u is not None:
                own1 = np.array(u['t2_1']) < u['q2']
                dd = (d1 < e['q_law']).astype(float) - own1.astype(float)      # referee round 73: the part's paired standard error
                summ[-1].update(u0=float(np.mean(np.array(u['t2_0']) < u['q2'])), u1=float(np.mean(own1)),
                                first=float(np.mean(u['first_0'])), se=float(np.std(dd, ddof=1) / np.sqrt(dd.size)),
                                nd=int(np.sum(dd != 0)), nn=int(dd.size))
                if 'q_full' in u:                                       # referee round 71: the two laws from the same samples
                    t2a, dta = np.r_[u['t2_0'], u['t2_1']].astype(float), np.r_[d0, d1].astype(float)
                    ok = np.isfinite(t2a) & (t2a > 0)                    # referee round 75: an infinite path knot beside a finite own
                    exc = np.where(np.isfinite(dta[ok]), dta[ok] / t2a[ok] - 1, np.inf)   # one is the worst exceedance, counted
                    summ[-1].update(dq=abs(u['q_full'] - e['q_law']), eq=bool(u['q2'] == u['q_full']),
                                    betw=int(np.sum((d1 < u['q2']) & ~own1)), n1=int(own1.size),
                                    exc=(int(np.sum(exc > 0)), float(exc.max()) if exc.size else 0.0, int(t2a.size), int(ok.sum())))
    dz = np.array([q['z'] - q['np'] for q in summ]); dk = np.array([q['z'] - q['cal'] for q in summ])
    dn = np.array([q['np'] - q['knot'] for q in summ])
    k0, n0 = sum(q['k0'] for q in summ), sum(q['n0'] for q in summ)
    ci = cp_int(k0, n0)
    t0 = sum(q['t0'] for q in summ); zs0 = sum(q['zs0'] for q in summ)
    fit_dev = max(max(q['dfit'] for q in r['rows0'] + r['rows1']) for r in Se.values())
    out += ['', f"Over the {len(summ)} cells: the score's power minus the Neyman–Pearson bound {dz.min():+.3f} … {dz.max():+.3f}, minus the "
                f"calibrated knot's {dk.min():+.3f} … {dk.max():+.3f} (positive in {int(np.sum(dk > 0))} of {len(summ)}); the bound minus "
                f"the larger of the knot's two powers {dn.min():.3f} … {dn.max():.3f} (within 0.05 in {int(np.sum(dn <= 0.05 + 1e-9))}); the "
                f"runner-up's filter's power {min(q['z2'] for q in summ):.3f} … {max(q['z2'] for q in summ):.3f}. The score's pooled "
                f"false-alarm rate {k0 / n0:.4f} [{ci[0]:.4f}, {ci[1]:.4f}] over {n0} H0 draws, its t-version's {t0 / n0:.4f}, the width's "
                f"{zs0 / n0:.4f}; the first fit by `one_fit` equals the path's to {fit_dev:.1e}.", '']
    uu = [q for q in summ if 'u1' in q]
    if uu:                                                              # referee round 64: the knot's losses decomposed
        far = [q for q in uu if q['cb'] <= 0.032]
        tie = [q for q in uu if q['cb'] == max(q2['cb'] for q2 in uu)]
        near = [q for q in uu if q['cb'] >= 0.045]
        out += [f"The runner-up's own knot (the tie pair's second knot against its own law at the truth): its size {min(q['u0'] for q in uu):.3f}"
                f"–{max(q['u0'] for q in uu):.3f}; at c̄ ≤ 0.032 its power minus the runner-up's filter z₂'s "
                f"{min(q['u1'] - q['z2'] for q in far):+.3f} … {max(q['u1'] - q['z2'] for q in far):+.3f}, and the larger of the knot's two "
                f"powers below it by {min(q['u1'] - q['knot'] for q in far):.3f} … {max(q['u1'] - q['knot'] for q in far):.3f} (the "
                f"dictionary's minimum); at c̄ = {tie[0]['cb']:g} its power {', '.join(f'{q['u1']:.3f}' for q in tie)} against z₂'s "
                f"{', '.join(f'{q['z2']:.3f}' for q in tie)} (the margin's noise).",
                "The minimum's part (referee round 67) — the full knot against its own law at the truth less the runner-up's own knot, per cell "
                "(c̄, Δ, s), with its standard error over the paired H1 draws (referee round 73; the laws' own Monte Carlo error not included — "
                "0.003–0.006 per cell by referee round 78's estimate, not reproduced here): "
                + '; '.join(f"({q['cb']:g}, {q['D']:g}, {q['s']:g}) {q['law1'] - q['u1']:+.4f}" + part_se(q) for q in uu)
                + f": at c̄ ≤ 0.032 {min(q['law1'] - q['u1'] for q in far):+.4f} … {max(q['law1'] - q['u1'] for q in far):+.4f} (pooled per "
                "offset over its cells, with its standard error: " + '; '.join(
                    f"c̄ = {cb:g} {np.mean([q['law1'] - q['u1'] for q in far if q['cb'] == cb]):+.4f} ("
                    f"{np.sqrt(sum(q['se'] ** 2 for q in far if q['cb'] == cb)) / sum(q['cb'] == cb for q in far):.4f})"
                    for cb in sorted({q['cb'] for q in far}) if all('se' in q for q in far if q['cb'] == cb))
                + (lambda pp: f"; the two lowest offsets' pooled parts differ by z = {(pp[0][0] - pp[1][0]) / np.hypot(pp[0][1], pp[1][1]):+.2f}"
                   if len(pp) >= 2 else '')([(np.mean([q['law1'] - q['u1'] for q in far if q['cb'] == cb]),
                                              np.sqrt(sum(q['se'] ** 2 for q in far if q['cb'] == cb)) / sum(q['cb'] == cb for q in far))
                                             for cb in sorted({q['cb'] for q in far}) if all('se' in q for q in far if q['cb'] == cb)]) + ")"
                + (f", at c̄ ≥ 0.045 {min(q['law1'] - q['u1'] for q in near):+.4f} … {max(q['law1'] - q['u1'] for q in near):+.4f}"
                   if near else '') + f". U's share of H0 draws whose first entrant is j₂ against S's swaps: "
                f"{max(abs(q['first'] - q['swap']) for q in uu):.4f} at most apart."
                + (f" f₁'s and f₀'s first entrants (S's normalisation of μ̂, U's of μ₀) coincide in "
                   f"{sum(u_['j1'] == u_['j10'] for u_ in U)} of {len(U)} cells (referee round 75; for E's pair at offsets below half the "
                   "spacing the two centres lie in the same kernel's half-cell, so they cannot differ: referee rounds 78, 82)." if U is not None and len(U)
                   and all('j10' in u_ for u_ in U) else '')
                + (f" Both laws from the same samples (referee round 71): the full knot's α quantile from U's samples equals E's at the truth "
                   f"to {max(q['dq'] for q in uu if 'dq' in q):.1e}, and the own knot's equals it in "
                   f"{sum(q['eq'] for q in uu if 'eq' in q)} of {sum('eq' in q for q in uu)} cells (c̄ = "
                   f"{', '.join(sorted({f'{q['cb']:g}' for q in uu if q.get('eq')}))}: the pair sets the whole tail). The full knot is "
                   "never later than the own knot — the minimum over the kernels' knots, the own knot among them (draw by draw, over all "
                   f"cells, the path's exceeds it in {sum(q['exc'][0] for q in uu if 'exc' in q)} of {sum(q['exc'][2] for q in uu if 'exc' in q)} "
                   f"draws ({sum(q['exc'][3] for q in uu if 'exc' in q)} with a finite own knot), by at most "
                   f"{max(q['exc'][1] for q in uu if 'exc' in q):.1e} relative: round-off) —, so at a common threshold the minimum's part is "
                   "the share of H1 draws whose two knots straddle it, ≥ 0: where the quantiles coincide the minimum cannot lose, and it gains "
                   f"nothing there because the threshold falls between the two knots of no H1 draw ({sum(q['betw'] for q in uu if q.get('eq'))} "
                   f"of {sum(q['n1'] for q in uu if q.get('eq'))} in those cells"
                   + (f", no discordant draw: the share at most {1 - 0.05 ** (1 / sum(q['n1'] for q in uu if q.get('eq'))):.4f} (one-sided 95%; "
                      f"{1 - 0.025 ** (1 / sum(q['n1'] for q in uu if q.get('eq'))):.4f} two-sided) for the cells' mean share"
                      if all(q.get('nd') == 0 for q in uu if q.get('eq')) else '') + "; referee rounds 73, 75)."
                   if any('dq' in q for q in uu) else ''), '']
    return out


def part_se(q):
    """The minimum's part's paired standard error — or, with no discordant draw, the count and the one-sided 95% bound 1 − 0.05^{1/n} on
    the share between the knots (referee round 75: a zero count's standard error is degenerate)."""
    if 'se' not in q:
        return ''
    if q.get('nd') == 0:
        return f" (0 discordant of {q['nn']}: ≤ {1 - 0.05 ** (1 / q['nn']):.4f} one-sided, {1 - 0.025 ** (1 / q['nn']):.4f} two-sided at 95%)"
    return f" ({q['se']:.3f})"


def own_knot(u):
    """U's cell: the runner-up's own knot against its own law at the truth, its size and power ('—' without U)."""
    if u is None:
        return '—'
    return f"{np.mean(np.array(u['t2_0']) < u['q2']):.3f}, {np.mean(np.array(u['t2_1']) < u['q2']):.3f}"


def null_distance(c, h):
    """The distance of a centre to the nearest kernel of the grid of spacing h (mirror-symmetric about the tie at h/2)."""
    return min(c % h, h - c % h)


def trend_z(groups):
    """The Cochran–Armitage trend statistic of rejections over groups (x, k, n) — score x, k rejections of n — and its two-sided p-value
    (referee round 78)."""
    x, k, n = (np.array([g[i] for g in groups], float) for i in range(3))
    N, K_ = n.sum(), k.sum()
    if N == 0 or K_ in (0, N) or len(groups) < 2:
        return None, None
    pb, xb = K_ / N, float(np.sum(n * x) / N)
    v = pb * (1 - pb) * float(np.sum(n * (x - xb) ** 2))
    if v <= 0:
        return None, None
    z = float(np.sum(k * (x - xb)) / np.sqrt(v))
    return z, float(math.erfc(abs(z) / np.sqrt(2)))


def logit_fit(X, k, n, it=100):
    """The binomial logistic regression of k of n on the design X by iteratively reweighted least squares: the coefficients, their
    covariance and the deviance (referee round 78); np.linalg.LinAlgError if the iteration does not converge (referee round 82)."""
    X, k, n = np.asarray(X, float), np.asarray(k, float), np.asarray(n, float)
    b = np.zeros(X.shape[1])
    conv = False
    for _ in range(it):
        with np.errstate(over='ignore'):
            p = 1 / (1 + np.exp(-(X @ b)))
        W = n * p * (1 - p)
        H = X.T @ (W[:, None] * X)
        step = np.linalg.solve(H, X.T @ (k - n * p))
        b = b + step
        if np.max(np.abs(step)) < 1e-12:
            conv = True
            break
    if not conv:                                                        # referee round 82: a separated design diverges without a
        raise np.linalg.LinAlgError('logit_fit: no convergence')       # singular matrix — reported as one
    p = 1 / (1 + np.exp(-(X @ b)))
    W = n * p * (1 - p)
    cov = np.linalg.inv(X.T @ (W[:, None] * X))
    with np.errstate(divide='ignore', invalid='ignore'):
        dv = 2 * np.sum(np.where(k > 0, k * np.log(k / (n * p)), 0.0) + np.where(n - k > 0, (n - k) * np.log((n - k) / (n * (1 - p))), 0.0))
    return b, cov, float(dv)


def quasi_logit(X, k, n):
    """`logit_fit` with the quasi-binomial dispersion (referee round 82: nulls at one distance disagree beyond binomial noise): φ = Pearson
    χ²/(m − p) over the m groups and the p coefficients; the standard errors scaled by √φ. Returns (b, se, φ, deviance)."""
    X, k, n = np.asarray(X, float), np.asarray(k, float), np.asarray(n, float)
    b, cov, dv = logit_fit(X, k, n)
    p = 1 / (1 + np.exp(-(X @ b)))
    phi = float(np.sum((k - n * p) ** 2 / (n * p * (1 - p))) / max(len(k) - X.shape[1], 1))
    return b, np.sqrt(np.diag(cov) * phi), phi, dv


def perm_trend(groups, n_perm=None, seed=8200):
    """The Cochran–Armitage z over groups (x, k, n) and its permutation p-value with the groups as the units (referee round 82): the
    scores x permuted over the groups n_perm times (default N_PERM), p = (1 + #{|z_π| ≥ |z|})/(n_perm + 1). (None, None) without a z."""
    z, _ = trend_z(groups)
    if z is None:
        return None, None
    n_perm = N_PERM if n_perm is None else n_perm
    rng = np.random.default_rng(seed)
    xs = np.array([g[0] for g in groups], float)
    hits = 0
    for _ in range(n_perm):
        zp, _ = trend_z([(x, g[1], g[2]) for x, g in zip(rng.permutation(xs), groups)])
        hits += zp is not None and abs(zp) >= abs(z) - 1e-12
    return z, (1 + hits) / (n_perm + 1)


def null_margin(K, c0, c, w, Phi=None):
    """The runner-up's margin at the one-bump null wφ(c) and its noise (referee round 82): the duals at the first entry u_j = ⟨fφ_j⟩/
    |⟨fφ_{j₁}⟩| (j₁ = argmax|⟨fφ_j⟩|), the runner-up j₂ the other kernel nearest its wall, its margin μ₀ = 1 − |u_{j₂}|, and the
    standard deviation per unit noise of μ̂ = 1 − sign(u_{j₂})⟨yφ_{j₂}⟩/|⟨yφ_{j₁}⟩| by the delta method, ‖φ_{j₂} − (⟨fφ_{j₂}⟩/⟨fφ_{j₁}⟩)
    φ_{j₁}‖/|⟨fφ_{j₁}⟩| (sums over the samples). For the pair (0.8, 0.6) j₂ is (i)'s noiseless second entrant and μ₀ U's mu0.
    Returns (j₂, μ₀, sd(μ̂)/s)."""
    Phi = K.phi(c0) if Phi is None else Phi
    g = Phi.T @ (w * K.phi(np.array([c]))[:, 0])
    j1 = int(np.argmax(np.abs(g)))
    m = 1 - np.abs(g / abs(g[j1]))
    m[j1] = np.inf
    j2 = int(np.argmin(m))
    return j2, float(m[j2]), float(np.linalg.norm(Phi[:, j2] - (g[j2] / g[j1]) * Phi[:, j1]) / abs(g[j1]))


def dist_bins():
    """The bins of the nulls' own distances to the nearest kernel (referee rounds 75, 78): the distances of OFFSETS_F's centroids
    (mirrors about the tie pooled) and the edges halfway between them, from 0 to the half spacing."""
    h = 2.0 / (X.D_C - 1)
    grid = sorted({round(null_distance(cb, h), 6) for cb in OFFSETS_F})
    return grid, [0.0] + [0.5 * (a + b) for a, b in zip(grid[:-1], grid[1:])] + [h / 2]


def bin_index(d, edges):
    """The bin of a distance d: lo ≤ d < hi, the last bin closed."""
    for i in range(len(edges) - 1):
        if edges[i] <= d and (d < edges[i + 1] or (i == len(edges) - 2 and d <= edges[i + 1])):
            return i
    return None


def size_nulls(F, E_=None):
    """Every H0 null of F (per row) and of E (per cell, the strong pair PAIR_C) with its own distance to the nearest kernel and its
    calibrated p-values: (source, pair, distance, the p-values, the null's |ŵ|) (referee rounds 73, 75, 78, 82)."""
    h = 2.0 / (X.D_C - 1)
    out = [('F', tuple(r['A']), null_distance(r.get('c1', r['cbar']), h), np.array([q['p'] for q in r['rows']]), abs(r.get('w1', np.nan)))
           for r in F if r['hyp'] == 0]
    out += [('E', tuple(PAIR_C), null_distance(e['c1'], h), np.array([q['p'] for q in e['rows0']]), abs(e.get('w1', np.nan)))
            for e in (E_ if E_ is not None else [])]
    return out


def size_by_distance(F, E_=None):
    """The calibrated knot's size by the H0 null's own distance to the nearest kernel (referee rounds 73, 75, 78), for the weak null (F's
    (0.8, −0.4), |ŵ| ≈ 0.4) and the strong one (0.8, 0.6) — F's draws, E's (the same null and calibration) and both pooled —: per bin of
    distances (the interval of distances nearest each offset's {d(c̄)}, mirrors about the tie pooled; its nulls' distances) the size at α
    (Clopper–Pearson 95%) and at 0.01; per null the trend over the nulls' own distances (Cochran–Armitage), and the two nulls' slopes in
    one logistic model (distance, pair and their interaction; the interaction's likelihood-ratio test, and the common slope's z)."""
    grid, edges = dist_bins()
    nulls = size_nulls(F, E_)
    out, data = [], []
    for A in (tuple(A) for A in PAIRS_F):
        for src in (('F',), ('E',), ('F', 'E')):
            ns = [x for x in nulls if x[1] == A and x[0] in src]
            if not ns or (src != ('F',) and not any(x[0] == 'E' for x in ns)):
                continue
            cells = []
            for i, g in enumerate(grid):
                lo, hi = edges[i], edges[i + 1]
                last = i == len(grid) - 1
                inb = [x for x in ns if lo <= x[2] and (x[2] < hi or (last and x[2] <= hi))]
                if not inb:
                    continue
                pc = np.concatenate([x[3] for x in inb])
                k5, k1 = int(np.sum(pc < ALPHA_C)), int(np.sum(pc < 0.01))
                ci = cp_int(k5, pc.size)
                ds = [x[2] for x in inb]
                cells.append(f"[{lo:.4g}, {hi:.4g}{']' if last else ')'} (nulls at {min(ds):.4f}–{max(ds):.4f}): "
                             f"{k5 / pc.size:.4f} [{ci[0]:.4f}, {ci[1]:.4f}], {k1 / pc.size:.4f}")
            z, p = trend_z([(x[2], int(np.sum(x[3] < ALPHA_C)), x[3].size) for x in ns])
            nd = sum(x[3].size for x in ns)
            ws = [x[4] for x in ns if np.isfinite(x[4])]
            out.append(f"Pair {A}, {'+'.join(src)} ({nd} H0 draws of {len(ns)} nulls"
                       + (f", the nulls' |ŵ| {min(ws):.4f}–{max(ws):.4f}" if ws else '') + "), the size by the null's own distance to the "
                       "nearest kernel (referee rounds 73, 75, 78): " + '; '.join(cells)
                       + ('' if z is None else f"; the trend over the nulls' own distances (Cochran–Armitage, binomial variance: the draws as "
                                               f"the units) z = {z:+.2f} (p = {p:.3f})") + '.')
            if src == ('F', 'E') or (src == ('F',) and not any(x[0] == 'E' for x in nulls if x[1] == A)):
                data += [(x[2], int(np.sum(x[3] < ALPHA_C)), x[3].size, A) for x in ns]
    prs = sorted({d[3] for d in data})
    if len(prs) == 2:
        g = np.array([1.0 if d[3] == prs[1] else 0.0 for d in data])
        dd = np.array([d[0] for d in data]); kk = np.array([d[1] for d in data]); nn = np.array([d[2] for d in data])
        try:
            bf, cf, dvf = logit_fit(np.c_[np.ones_like(dd), dd, g, g * dd], kk, nn)
            br, cr, dvr = logit_fit(np.c_[np.ones_like(dd), dd, g], kk, nn)
        except np.linalg.LinAlgError:                                   # a separated design: no finite fit
            return out
        lr = max(dvr - dvf, 0.0)
        out.append(f"The two nulls in one logistic model of rejection at α over the nulls' own distances (the strong null's draws E and F "
                   f"pooled): the slopes' difference (pair {prs[1]} against {prs[0]}) z = {bf[3] / np.sqrt(cf[3, 3]):+.2f}, its likelihood-"
                   f"ratio p = {math.erfc(np.sqrt(lr / 2)):.3f}; the common slope (without the interaction) z = {br[1] / np.sqrt(cr[1, 1]):+.2f}.")
    return out


def margin_nulls(F, E_=None, G=None):
    """Every H0 null of F, E (the pair PAIR_C) and G with its runner-up's margin (referee round 82): source, pair, s, c̄, Δ, the null's own
    distance to the nearest kernel and |ŵ|, `null_margin` at its one-bump fit (j₂, μ₀, sd(μ̂)/s) and μ₀/sd(μ̂); its draws' calibrated
    p-values, knots and the law at the truth's α quantile, and per draw the margin at the draw's own first fit (ĉ, ŵ) — the null the
    bootstrap draws from."""
    L_, K, c0 = setup()
    Phi = K.phi(c0)
    h = 2.0 / (X.D_C - 1)
    src = ([('F', r, tuple(r['A']), r['rows']) for r in F if r['hyp'] == 0]
           + [('E', e, tuple(PAIR_C), e['rows0']) for e in (E_ if E_ is not None else [])]
           + [('G', r, tuple(r['A']), r['rows']) for r in (G if G is not None else []) if r['hyp'] == 0])
    out = []
    for nm, r, A, rows in src:
        j2, mu0, sd = null_margin(K, c0, r['c1'], r['w1'], Phi)
        out.append(dict(src=nm, A=A, s=r['s'], cbar=r['cbar'], D=r['D'], dist=null_distance(r['c1'], h), w1=abs(r['w1']), j2=j2, mu0=mu0,
                        sd=sd, ratio=mu0 / (sd * r['s']), p=np.array([q['p'] for q in rows]), dt=np.array([q['dt'] for q in rows]),
                        q_law=r['q_law'], fit=np.array([null_margin(K, c0, q['c1'], q['w1'], Phi)[1] for q in rows])))
    return out


def in_band(x):
    """A null in the band BAND_G of μ₀/sd(μ̂) (referee round 82)."""
    return BAND_G[0] < x['ratio'] < BAND_G[1]


def pooled(ns, key='p'):
    """Over the nulls ns: (rejections, draws) of the calibrated test (key 'p': p < α) or of the knot against the law at the truth ('t')."""
    k = sum(int(np.sum(x['p'] < ALPHA_C)) if key == 'p' else int(np.sum(x['dt'] < x['q_law'])) for x in ns)
    return k, sum(x['p'].size for x in ns)


def fmt_rate(k, n):
    lo, hi = cp_int(k, n)
    return f"{k / n:.4f} [{lo:.4f}, {hi:.4f}]"


def margin_lines(F, E_=None, G=None):
    """F″ (referee round 82): the calibrated knot's size against the runner-up's margin at the null over its noise — per null (`margin_nulls`)
    its μ₀/sd(μ̂), its size and the same draws against the law at the truth, and the draws the law at the truth rejects (how many the
    calibrated test rejects too, and the margin at their own first fits); pooled in and outside the band BAND_G; with the nulls as the
    units, the trend in the distance (permutation, quasi-binomial) and the band against it; G's nulls, chosen by μ₀/sd(μ̂) before they were
    run, against the band's and the distance's predictions."""
    if F is None or len(F) == 0:
        return []
    ns = margin_nulls(F, E_, G)
    base = [x for x in ns if x['src'] != 'G']
    gn = [x for x in ns if x['src'] == 'G']
    out = ['', "### F″. The calibrated knot's size against the runner-up's margin at the null (referee round 82)", '',
           "Every H0 null — F's, E's (the pair (0.8, 0.6)) and G's (below) — at its one-bump fit: the runner-up j₂ (the kernel, other than the "
           "first entrant, whose dual at the first entry is nearest its wall), its margin μ₀ = 1 − |u_{j₂}|, the standard deviation of its "
           "estimate μ̂ per unit noise by the delta method (`null_margin`; for the pair (0.8, 0.6) j₂ is (i)'s noiseless second entrant and μ₀ "
           "E′'s), and μ₀/sd(μ̂); the calibrated knot's size at α (Clopper–Pearson 95%) and the same draws' knots against the law at the "
           "truth; the draws the law at the truth rejects, how many of them the calibrated test rejects, and the median of the margin at their "
           "own first fits (the null the bootstrap draws from) against all draws':", '',
           "| pair | null | s | c̄ | distance | j₂ | μ₀ | sd(μ̂)/s | μ₀/sd(μ̂) | draws | calibrated size [CP] | against the law at the truth | "
           "rejected by the law at the truth: calibrated too; their fitted margin, all draws' |",
           '|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for x in sorted(ns, key=lambda x: (x['A'], x['src'] == 'G', x['ratio'])):
        k, n = pooled([x])
        tr = x['dt'] < x['q_law']
        out.append(f"| {x['A']} | {x['src']} | {x['s']:g} | {x['cbar']:g} | {x['dist']:.4f} | {x['j2']} | {x['mu0']:.4f} | {x['sd']:.4f} | "
                   f"{x['ratio']:.2f} | {n} | {fmt_rate(k, n)} | {np.mean(tr):.4f} | {int(np.sum(tr & (x['p'] < ALPHA_C)))} of {int(tr.sum())}; "
                   + (f"{np.median(x['fit'][tr]):.4f}" if tr.any() else '—') + f", {np.median(x['fit']):.4f} |")
    parts = []
    for A in sorted({x['A'] for x in base}) + ['both']:
        sel = [x for x in base if A == 'both' or x['A'] == A]
        for lab, g in (('in the band', [x for x in sel if in_band(x)]), ('outside it', [x for x in sel if not in_band(x)])):
            if g:
                parts.append(f"{'both pairs' if A == 'both' else f'pair {A}'} {lab}: {len(g)} nulls, {fmt_rate(*pooled(g))}, against the "
                             f"law at the truth {fmt_rate(*pooled(g, 't'))}")
    out += ['', f"Pooled by the band {BAND_G[0]:g} < μ₀/sd(μ̂) < {BAND_G[1]:g} (the round-82 referee's, read off E's two lowest cells) over "
                "F's and E's nulls — the calibrated size [CP], and the same draws against the law at the truth: " + '; '.join(parts) + '.']
    mech = []
    for lab, g in (('in the band', [x for x in base if in_band(x)]), ('outside it', [x for x in base if not in_band(x)])):
        tr = [(x['dt'] < x['q_law'], x) for x in g]
        kt = sum(int(t.sum()) for t, _ in tr)
        kc = sum(int(np.sum(t & (x['p'] < ALPHA_C))) for t, x in tr)
        fr = [np.median(x['fit'][t]) / np.median(x['fit']) for t, x in tr if t.any() and np.median(x['fit']) > 0]
        if g:
            mech.append(f"{lab} the calibrated test makes {kc} of the law at the truth's {kt} rejections, and their fitted margin is "
                        f"{np.median(fr):.2f} of all draws' (the median over the nulls of the ratio of the medians; {min(fr):.2f}–{max(fr):.2f})")
    if mech:
        out += ['', "Where the calibration loses its rejections (referee round 82): " + '; '.join(mech) + '.']
    tests = []
    for A in sorted({x['A'] for x in base}):
        sel = [x for x in base if x['A'] == A]
        grp = lambda g: [(x['dist'], int(np.sum(x['p'] < ALPHA_C)), x['p'].size) for x in g]
        z, pp = perm_trend(grp(sel))
        dd = np.array([x['dist'] for x in sel]); bb = np.array([1.0 if in_band(x) else 0.0 for x in sel])
        kk = np.array([pooled([x])[0] for x in sel]); nn = np.array([x['p'].size for x in sel])
        try:
            bd, sd_d, phd, dvd = quasi_logit(np.c_[np.ones_like(dd), dd], kk, nn)
            bf, sd_f, phf, dvf = quasi_logit(np.c_[np.ones_like(dd), bb, dd], kk, nn)
            bo, sd_o, pho, dvo = quasi_logit(np.c_[np.ones_like(dd), bb], kk, nn)
        except np.linalg.LinAlgError:
            continue
        out_ = [x for x in sel if not in_band(x)]
        zo, po = perm_trend(grp(out_))
        try:
            bq, sd_q, phq, _ = quasi_logit(np.c_[np.ones(len(out_)), [x['dist'] for x in out_]], [pooled([x])[0] for x in out_],
                                           [x['p'].size for x in out_])
            outside = (f"; outside the band ({len(out_)} nulls) the distance's z = {bq[1] / sd_q[1]:+.2f} (φ = {phq:.2f}), the permutation "
                       f"p = {po:.3f}")
        except np.linalg.LinAlgError:
            outside = ''
        m = len(sel)
        tests.append(f"pair {A} ({m} nulls): the trend in the distance — Cochran–Armitage z = {z:+.2f} with the draws as the units, its "
                     f"permutation p = {pp:.3f} with the nulls as the units; the quasi-binomial slope z = {bd[1] / sd_d[1]:+.2f} (dispersion "
                     f"φ = {phd:.2f}, deviance {dvd:.1f} on {m - 2} df) —; the band alone: deviance {dvo:.1f} on {m - 2} df (φ = {pho:.2f}), its "
                     f"z = {bo[1] / sd_o[1]:+.2f}; both: the band's z = {bf[1] / sd_f[1]:+.2f}, the distance's {bf[2] / sd_f[2]:+.2f} (φ = "
                     f"{phf:.2f})" + outside)
    if tests:
        out += ['', "With the nulls as the units (referee round 82): " + '; '.join(tests) + '.']
    if gn:
        grid, edges = dist_bins()
        rows = []
        for x in gn:
            A = x['A']
            same = [y for y in base if y['A'] == A]
            band = [y for y in same if in_band(y) == in_band(x)]
            i = bin_index(x['dist'], edges)
            dist_ = [y for y in same if bin_index(y['dist'], edges) == i]
            kb, nb = pooled(band); kd, nd = pooled(dist_)
            rows.append(f"| {A} | {x['s']:g} | {x['cbar']:g} | {x['D']:g} | {x['dist']:.4f} | {x['ratio']:.2f} | "
                        f"{'yes' if in_band(x) else 'no'} | " + (f"{kb / nb:.4f}" if nb else '—') + ' | ' + (f"{kd / nd:.4f}" if nd else '—')
                        + f" | {fmt_rate(*pooled([x]))} | {fmt_rate(*pooled([x], 't'))} |")
        gi, go = [x for x in gn if in_band(x)], [x for x in gn if not in_band(x)]
        out += ['', "G — H0 nulls chosen by μ₀/sd(μ̂) before they were run (the first two the round-82 referee's), R0_G draws each, the same "
                    "test (`f_job`): the band's prediction (F's and E's nulls of the pair on the same side of the band, pooled) and the "
                    "distance's (those in the same distance bin as `size_by_distance`'s), against the calibrated size and the law at the "
                    "truth:", '', '| pair | s | c̄ | Δ | distance | μ₀/sd(μ̂) | in the band | the band predicts | the distance predicts | '
                                  'calibrated size [CP] | against the law at the truth [CP] |', '|---|---|---|---|---|---|---|---|---|---|---|'] + rows
        out += ['', "G pooled: " + '; '.join(f"{lab} {len(g)} nulls, {fmt_rate(*pooled(g))}, against the law at the truth "
                                               f"{fmt_rate(*pooled(g, 't'))}" for lab, g in (('in the band', gi), ('outside it', go)) if g) + '.']
    return out + ['']


def f_alloc_lines(F, least, small, names, keys, E_=None):
    """F's allocation (referee round 64): under H0 the share of the calibrated knot's rejections at the wall of the fit's sign (tA < tB)
    per pair, pooled and per offset; each test's least power over Δ per cell (the cells below α) and its power pooled over the cells at
    the two smallest separations (Clopper–Pearson); the sizes at α = 0.01 pooled over the H0 draws."""
    out = []
    for A in PAIRS_F:
        h0 = [r for r in F if tuple(r['A']) == tuple(A) and r['hyp'] == 0]
        rej = [(q['tA'] < q['tB'], r['cbar']) for r in h0 for q in r['rows'] if q['p'] < ALPHA_C]
        per = {cb: [a for a, c in rej if c == cb] for cb in OFFSETS_F}
        out.append(f"Pair {A} under H0: of the calibrated knot's {len(rej)} rejections {sum(a for a, _ in rej)} at the wall of the fit's "
                   f"sign ({sum(a for a, _ in rej) / max(len(rej), 1):.3f}; per offset "
                   + ', '.join(f"{cb:g}: {sum(v)}/{len(v)}" for cb, v in per.items() if v) + ").")
    out.append('')
    for A in PAIRS_F:
        cells = []
        for k in keys:
            lw = least[(A, k)]
            if not lw:
                continue
            k0, n0 = small[(A, k, 0)]; k1, n1 = small[(A, k, 1)]
            c0, c1 = cp_int(k0, n0), cp_int(k1, n1)
            cells.append(f"{names[k]}: least {min(lw):.3f}–{max(lw):.3f}, below α in {sum(x < ALPHA_C for x in lw)} of {len(lw)}; "
                         f"pooled {k0 / max(n0, 1):.3f} [{c0[0]:.3f}, {c0[1]:.3f}], {k1 / max(n1, 1):.3f} [{c1[0]:.3f}, {c1[1]:.3f}]")
        out.append(f"Pair {A}, each test's least power over Δ per cell and its power pooled over the cells at the two smallest separations: "
                   + '; '.join(cells) + '.')
    out.append('')
    h0 = [q for r in F if r['hyp'] == 0 for q in r['rows']]
    if h0:
        z = np.array([q['z'] for q in h0]); w = np.sign([q['w1'] for q in h0]); p = np.array([q['p'] for q in h0])
        n = len(h0)
        sz = [('knot (calibrated)', int(np.sum(p < 0.01))), ('z, one-sided toward the same sign', int(np.sum(w * z > Z99))),
              ('z, one-sided toward the opposite sign', int(np.sum(-w * z > Z99))), ('\\|z\\|', int(np.sum(np.abs(z) > Z995)))]
        out.append(f"The sizes at α = 0.01 over the {n} H0 draws (Clopper–Pearson 95%): "
                   + '; '.join(f"{nm} {k / n:.4f} [{cp_int(k, n)[0]:.4f}, {cp_int(k, n)[1]:.4f}]" for nm, k in sz) + '.')
        per = []                                                       # referee round 67: the calibrated knot's size per offset
        for cb in OFFSETS_F:
            pc = np.array([q['p'] for r in F if r['hyp'] == 0 and r['cbar'] == cb for q in r['rows']])
            if pc.size:
                k5, k1 = int(np.sum(pc < ALPHA_C)), int(np.sum(pc < 0.01))
                per.append(f"c̄ = {cb:g}: {k5 / pc.size:.4f} [{cp_int(k5, pc.size)[0]:.4f}, {cp_int(k5, pc.size)[1]:.4f}], "
                           f"{k1 / pc.size:.4f}")
        out.append(f"The calibrated knot's size per offset c̄ — by the label, the two pairs' H0 draws pooled (referee round 78: by the null's "
                   f"own distance and per null below) — (α = {ALPHA_C:g} with Clopper–Pearson 95%, then α = 0.01): "
                   + '; '.join(per) + '.')
        out += size_by_distance(F, E_)
        out.append('')
    comp = []                                                           # referee round 67: the walls test's own components at α/2
    for A in PAIRS_F:
        for h in (0, 1):
            rs = [r for r in F if tuple(r['A']) == tuple(A) and r['hyp'] == h and (h == 0 or r['D'] in DELTAS_F[r['s']][:2])]
            a = [q['tA'] < r['qA'] for r in rs for q in r['rows']]
            b = [q['tB'] < r['qB'] for r in rs for q in r['rows']]
            if a:
                comp.append((A, h, float(np.mean(a)), float(np.mean(b))))
    if comp:
        out.append(f"The walls test's components, each at α/2 against its law at the truth — the fit's-sign wall and the other — under H0 and "
                   "under H1 at the two smallest separations: " + '; '.join(
                       f"pair {A}, H{h}: {a:.4f}, {b:.4f}" for A, h, a, b in comp) + '.')
        out.append('')
    return out


def equal_size_line(F, A, d50z):
    """Referee round 71: the one-sided half at equal size — the wall the pair's runner-up goes to (the fit's sign for a same-sign pair,
    p_A; the other for an opposite-sign pair, p_B), calibrated at the fit (p < α′) at the level α′ (a grid of step 0.0005) whose H0 size,
    pooled over the pair's H0 draws, is nearest the one-sided z's; its Δ₅₀ over the Neyman–Pearson bound's per cell and paired against
    z's (d50z)."""
    sgn = float(np.sign(A[0] * A[1]))
    h0 = [q for r in F if tuple(r['A']) == tuple(A) and r['hyp'] == 0 for q in r['rows']]
    if not h0:
        return f"Pair {A} at equal size: —"
    key, wall = ('pA', "the fit's-sign wall") if sgn > 0 else ('pB', 'the other wall')
    pA = np.array([q[key] for q in h0])
    size_z = float(np.mean(sgn * np.sign([q['w1'] for q in h0]) * np.array([q['z'] for q in h0]) > Z95))
    grid = np.round(np.arange(0.030, 0.0801, 0.0005), 4)
    a_eq = float(grid[int(np.argmin([abs(float(np.mean(pA < a)) - size_z) for a in grid]))])
    d50 = {}
    for s_ in SIGMAS_F:
        Ds = DELTAS_F[s_]
        for cb in OFFSETS_F:
            h1 = {r['D']: r for r in F if tuple(r['A']) == tuple(A) and r['cbar'] == cb and r['s'] == s_ and r['hyp'] == 1}
            if len(h1) < len(Ds):
                continue
            x = crossing(list(Ds), [float(np.mean([q[key] < a_eq for q in h1[D]['rows']])) for D in Ds])
            xn = crossing(list(Ds), [h1[D]['np'] for D in Ds])
            if np.isfinite(x) and np.isfinite(xn):
                d50[(s_, cb)] = x / xn
    q = np.array(list(d50.values()))
    return (f"Pair {A} at equal size (referee round 71): {wall} at p < α′ = {a_eq:.4f}, its H0 size {np.mean(pA < a_eq):.4f} "
            f"against the one-sided z's {size_z:.4f} (at α = {ALPHA_C:g}: {np.mean(pA < ALPHA_C):.4f}) — Δ₅₀ over the Neyman–Pearson bound's "
            + (f"{np.median(q):.3f} [{q.min():.3f}, {q.max():.3f}]; " if q.size else '—; ')
            + paired(d50, d50z, "paired against z, one-sided") + '.')


def knot_allocation(F, A):
    """Referee round 67: the share of the calibrated knot's rejections under H0 at the wall of the fit's sign (tA < tB), pooled over the
    pair A's H0 draws — how the knot spends α."""
    rej = [q['tA'] < q['tB'] for r in F if tuple(r['A']) == tuple(A) and r['hyp'] == 0 for q in r['rows'] if q['p'] < ALPHA_C]
    return sum(rej) / len(rej) if rej else 1.0


def paired(a, b, name):
    """Two tests' Δ₅₀ (dicts by cell) paired over the cells where both cross: the median [range] of a/b and the cells where a < b."""
    ks = [k for k in a if k in b]
    if not ks:
        return f"{name} —"
    q = np.array([a[k] / b[k] for k in ks])
    return f"{name} {np.median(q):.3f} [{q.min():.3f}, {q.max():.3f}], smaller in {int(np.sum(q < 1))} of {len(ks)}"


def knot_of_kernel(K, c0, c, w, Xi, j, alt):
    """The second knot of the tie pair {alt, j} (referee round 64: the runner-up's own knot), first order in the noise at the one-bump fit
    (c, w): of j — the noiseless runner-up —, or of alt — the noiseless first entrant — where j enters first; its dual at the first entry
    u_k and its rate n_k = ⟨Rξ, φ_k⟩, the time (1 − sign(n_k)u_k)/|n_k| to the wall it moves toward."""
    N, M = K.xs.size, Xi.shape[0]
    Phi = K.phi(c0)
    v0 = (w * K.phi(np.array([c]))[:, 0][None, :] + Xi) @ Phi / N
    k = np.where(np.argmax(np.abs(v0), axis=1) == j, alt, j)
    u = v0[np.arange(M), k] / np.abs(v0).max(axis=1)
    n = (Xi @ (tangent_projector(K, c, w) @ Phi))[np.arange(M), k] / N
    with np.errstate(divide='ignore'):
        return (1 - np.sign(n) * u) / np.abs(n)


def u_job(args):
    """U (referee round 64): E's draws regenerated from their seeds (as `s_job`), the path's own first fits from E's rows: per draw the
    noiseless second entrant j₂'s own knot from the fit (exact in the sliding limit: its dual at the first entry, moving at ⟨r₁φ_{j₂}⟩; the
    noiseless first entrant's where j₂ enters first), against its own law at the truth (`knot_of_kernel`, M_LAW samples: the α quantile),
    and the share of the draws whose first entrant is j₂. Referee round 71: the law from E's own law samples (rng [seed, 9], `e_job`'s
    law0), so that the full knot's law and the own knot's share their samples — at the tie the pair sets the whole tail and the two α
    quantiles coincide (with independent samples their Monte Carlo errors made a spurious difference) —; the full knot's quantile from
    the same samples (q_full: E's q_law) and the null's margin mu0 = 1 − wall·u_{j₂} at f₀ (u over the null's own first entrant's dual,
    argmax|Φᵀf₀|: referee round 73)."""
    D, s, cbar, seed, ns, fits = args
    L, K, c0 = setup()
    N = K.xs.size
    f1 = X.bumps_target(K, pair_centres(*PAIR_C, D, cbar), PAIR_C)
    c1, w1 = one_fit(K, f1)
    f0 = w1 * K.phi(np.array([c1]))[:, 0]
    rng = np.random.default_rng(seed)
    xi = [rng.standard_normal(N) for _ in range(sum(ns))]
    pd = predicted_delay(K, f1, c0, *PAIR_C, D, cbar=cbar)
    j2 = pd[1]
    Phi = K.phi(c0)
    j1 = int(np.argmax(np.abs(Phi.T @ f1)))
    Xi = s * np.random.default_rng([seed, 9]).standard_normal((M_LAW, N))     # E's law samples (referee round 71)
    q2 = float(np.quantile(knot_of_kernel(K, c0, c1, w1, Xi, j2, j1), ALPHA_C))
    g0 = Phi.T @ f0
    j10 = int(np.argmax(np.abs(g0)))                                   # the null's own first entrant (referee round 73; j₁ in every cell)
    out = dict(D=D, s=s, cbar=cbar, j2=int(j2), q2=q2, q_full=float(np.quantile(knot_of(K, c0, c1, w1, Xi), ALPHA_C)),
               mu0=float(1 - pd[2] * g0[j2] / abs(g0[j10])), j1=j1, j10=j10)
    k = 0
    for h, n in zip((0, 1), ns):
        t2, first = [], []
        for r in range(n):
            y = (f1 if h == 1 else f0) + s * xi[k]
            cf, wf = fits[k]
            g = Phi.T @ y / N
            jj = j1 if int(np.argmax(np.abs(g))) == j2 else j2
            v = Phi[:, jj] @ (y - wf * K.phi(np.array([cf]))[:, 0]) / N
            u = g[jj] / np.abs(g).max()
            t2.append(float((1 - np.sign(v) * u) / abs(v)) if v != 0 else np.inf)
            first.append(jj != j2)
            k += 1
        out[f't2_{h}'], out[f'first_{h}'] = t2, first
    print(f"  U Δ {D:g} s {s:g} c̄ {cbar:g}: size {np.mean(np.array(out['t2_0']) < q2):.3f}, power {np.mean(np.array(out['t2_1']) < q2):.3f}",
          flush=True)
    return out


def part_U(E_):
    """U's jobs: E's cells (seeds as `part_E`'s; the path's first fits from E's rows)."""
    from multiprocessing import Pool
    eseed = {(D, s, cb): sd for D, s, cb, sd in e_jobs()}
    jobs = [(e['D'], e['s'], e['cbar'], eseed[(e['D'], e['s'], e['cbar'])], (len(e['rows0']), len(e['rows1'])),
             [(q['c1'], q['w1']) for q in e['rows0']] + [(q['c1'], q['w1']) for q in e['rows1']]) for e in E_]
    with Pool(N_WORKERS) as pool:
        return pool.map(u_job, jobs, chunksize=1)


def f_lines(F, E_=None):
    """F (referee round 59): per pair, noise and centroid offset the separation of power ½ of each test (law-calibrated knot, the knot
    calibrated per draw, both walls, the score z toward the pair's sign, |z|, the width's |z_σ|) against the Neyman–Pearson bound's, and
    the tests' false-alarm rates on the H0 draws; the opposite-sign pair's powers against Δ at c̄ = CBAR."""
    out = ['## F. Resolution across the grid cell, for both signs', '',
           f"The pairs {', '.join(f'({a:g}, {b:g})' for a, b in PAIRS_F)} at the centroids c̄ = {', '.join(f'{c:g}' for c in OFFSETS_F)} (both "
           f"halves of a grid cell: c̄ = 0.05 is a tie), noise s = {', '.join(f'{s:g}' for s in SIGMAS_F)}, the separations "
           + '; '.join(f"s = {s:g}: {', '.join(f'{D:g}' for D in DELTAS_F[s])}" for s in SIGMAS_F)
           + f"; {R1_F} H1 draws per cell, {R0_F} H0 draws at the third separation. The tests at the level {ALPHA_C:g}: the knot against its "
           f"law at the truth ({M_LAW} samples), the knot calibrated per draw by its law at the first fit ({M_E} samples), both walls' knots "
           "against their laws at α/2 each, the score z one-sided toward the pair's sign sign(A₁A₂), |z| and the width's |z_σ| two-sided; Δ₅₀ "
           "where the power crosses ½ (log-linear interpolation; '<' or '>' the grid's ends), and the false-alarm rates at the third "
           "separation:", '']
    keys = ('knot', 'cal', 'own', 'other', 'walls', 'z', 'z2', 'zs', 'alloc')
    names = dict(knot='knot (law)', cal='knot (calibrated)', own="the fit's-sign wall (one-sided)", other='the other wall (one-sided)',
                 walls='both walls', z='z, one-sided', z2='\\|z\\|', zs='\\|z_σ\\|', alloc="z with the knot's allocation")
    alloc = {A: knot_allocation(F, A) for A in PAIRS_F}
    d50s = {(A, k): {} for A in PAIRS_F for k in keys}
    pool = {k: [0, 0] for k in keys}
    ratios = {(A, k): [] for A in PAIRS_F for k in keys}
    least = {(A, k): [] for A in PAIRS_F for k in keys}
    small = {(A, k, i): [0, 0] for A in PAIRS_F for k in keys for i in (0, 1)}
    for A in PAIRS_F:
        sgn = float(np.sign(A[0] * A[1]))
        out += [f"Pair {A}:", '', '| s | c̄ | Δ₅₀: ' + ', '.join(names[k] for k in keys) + ', Neyman–Pearson | false alarms: '
                + ', '.join(names[k] for k in keys) + " | the knot's least power (Δ) |", '|---|---|---|---|---|']
        for s_ in SIGMAS_F:
            Ds = DELTAS_F[s_]
            for cb in OFFSETS_F:
                h1 = {r['D']: r for r in F if tuple(r['A']) == tuple(A) and r['cbar'] == cb and r['s'] == s_ and r['hyp'] == 1}
                h0 = [r for r in F if tuple(r['A']) == tuple(A) and r['cbar'] == cb and r['s'] == s_ and r['hyp'] == 0]
                if len(h1) < len(Ds):
                    continue
                rej = {D: rejections(h1[D]['rows'], sgn, h1[D]['q_law'], h1[D]['qA'], h1[D]['qB'], alloc=alloc[A]) for D in Ds}
                for k in keys:
                    least[(A, k)].append(min(float(np.mean(rej[D][k])) for D in Ds))
                    for i, D in enumerate(Ds[:2]):
                        small[(A, k, i)][0] += int(rej[D][k].sum()); small[(A, k, i)][1] += rej[D][k].size
                pw = {k: [float(np.mean(rej[D][k])) for D in Ds] for k in keys}
                npw = [h1[D]['np'] for D in Ds]
                d50np = crossing(list(Ds), npw)
                for k in keys:
                    x = crossing(list(Ds), pw[k])
                    if np.isfinite(x) and np.isfinite(d50np):
                        ratios[(A, k)].append(x / d50np)
                        d50s[(A, k)][(s_, cb)] = x / d50np
                fa = ''
                if h0:
                    r0 = rejections(h0[0]['rows'], sgn, h0[0]['q_law'], h0[0]['qA'], h0[0]['qB'], alloc=alloc[A])
                    fa = ', '.join(f"{np.mean(r0[k]):.3f}" for k in keys)
                    for k in keys:
                        pool[k][0] += int(r0[k].sum()); pool[k][1] += r0[k].size
                j = int(np.argmin(pw['knot']))
                out.append(f"| {s_:g} | {cb:g} | " + ', '.join(d50_text(Ds, pw[k]) for k in keys) + f", {d50_text(Ds, npw)} | {fa} | "
                           f"{pw['knot'][j]:.3f} ({Ds[j]:g}) |")
        out.append('')
        out.append(f"Pair {A}, Δ₅₀ over the Neyman–Pearson bound's (median [range] over the cells where both cross): "
                   + '; '.join(f"{names[k]} {np.median(ratios[(A, k)]):.3f} [{min(ratios[(A, k)]):.3f}, {max(ratios[(A, k)]):.3f}]"
                               if ratios[(A, k)] else f"{names[k]} —" for k in keys) + '.')
        from scipy.stats import norm
        out.append(f"Pair {A}, like for like (referee round 67; the score with the knot's allocation spends α·{alloc[A]:.3f} toward the fit's "
                   f"sign, at z > {norm.ppf(1 - ALPHA_C * alloc[A]):.3f}, and α·{1 - alloc[A]:.3f} away, at −z > "
                   f"{norm.ppf(1 - ALPHA_C * (1 - alloc[A])):.3f}): paired over the cells where both cross, the first's Δ₅₀ over the "
                   "second's (median [range]) and the cells where the first is the smaller — " + '; '.join(
                       paired(d50s[(A, a)], d50s[(A, b)], f"{names[a]} against {names[b]}")
                       for a, b in (('cal', 'alloc'), ('walls', 'z2'), ('own', 'z'), ('cal', 'z2'))) + '.')
        out.append(equal_size_line(F, A, d50s[(A, 'z')]))
        out.append('')
    rng_ = [r['s'] / abs(q['w1']) for r in F for q in r['rows']]
    out.append(f"The noise over the one-bump fit's coefficient, s/|ŵ|, over F's draws (each at its own fit; referee rounds 67 and 71): "
               f"{min(rng_):.4f}–{max(rng_):.4f} (the range in which the score's first order is tested).")
    out.append("False-alarm rates pooled over the H0 draws of both pairs (Clopper–Pearson 95%): "
               + '; '.join(f"{names[k]} {pool[k][0] / max(pool[k][1], 1):.4f} [{cp_int(*pool[k])[0]:.4f}, {cp_int(*pool[k])[1]:.4f}]"
                           for k in keys if pool[k][1]) + '.')
    out += [''] + f_alloc_lines(F, least, small, names, keys, E_)
    A = PAIRS_F[-1]
    sgn = float(np.sign(A[0] * A[1]))
    out += ['', f"The pair {A} at c̄ = {CBAR:g}: the powers against Δ (knot (law), both walls, z one-sided, |z|; Neyman–Pearson) — with "
                "the null's centre (the pair's one-bump fit) against c̄, the law's α quantile and the knot's median under H1, each times s "
                "(referee round 64: the null moves across the tie as Δ grows):", '']
    for s_ in SIGMAS_F:
        cells = []
        for D in DELTAS_F[s_]:
            h = [r for r in F if tuple(r['A']) == tuple(A) and r['cbar'] == CBAR and r['s'] == s_ and r['D'] == D and r['hyp'] == 1]
            if h:
                rj = rejections(h[0]['rows'], sgn, h[0]['q_law'], h[0]['qA'], h[0]['qB'])
                cells.append(f"Δ = {D:g}: {np.mean(rj['knot']):.3f}, {np.mean(rj['walls']):.3f}, {np.mean(rj['z']):.3f}, "
                             f"{np.mean(rj['z2']):.3f}; {h[0]['np']:.3f} (the null's ĉ − c̄ {h[0]['c1'] - CBAR:+.4f}, q_law·s "
                             f"{h[0]['q_law'] * s_:.2f}, the knot's median·s {np.median([q['dt'] for q in h[0]['rows']]) * s_:.2f})")
        out.append(f"- s = {s_:g}: " + '; '.join(cells))
    out.append('')
    return out


def report(A=None, B=None, C=None, E_=None, S=None, F=None, U=None, G=None):
    fmt = E.fmt
    L_ = ['# Proposition D‴.52: close features (`exp_close_features.py`)', '',
          f"The sliding inverse scale space (polished) on two bumps of the kernels' width {X.SIG_C} a distance Δ apart, the centroid (the "
          f"midpoint for a dipole) at {CBAR}; 61 samples of [−1, 1], 21 kernels on a grid of spacing 0.1, equal masses.", '']
    if A is not None:
        L_ += ['## A. The second entry', '',
               '| A₁, A₂ | Δ | entries (τ, kernel) | Δτ₂ | predicted (kernel, wall) | ratio | with T₃ | Δτ₂Δ² | Δτ₂Δ | the end: kernels, centres '
               '(truth), rms |', '|---|---|---|---|---|---|---|---|---|---|']
        La, Ka, c0a = setup()
        for r in A:
            ent = r['entries']
            dt = ent[1][0] - ent[0][0] if len(ent) > 1 else np.nan
            p, sh = r.get('pred'), r.get('shape')
            r3 = '—'
            if p:
                pr, ratio = f"{p[0]:.4g} ({p[1]}, {p[2]:+d})", f"{dt / p[0]:.4f}"
                p3 = predicted_delay(Ka, X.bumps_target(Ka, pair_centres(*r['A'], r['D']), r['A']), c0a, *r['A'], r['D'], order=3)
                r3 = f"{dt / p3[0]:.4f}" + (f" ({p3[1]}, {p3[2]:+d})" if p3[1] != p[1] else '')
            elif sh:
                pr = f"shape: τ₁ {sh['tau1']:.4g}, Δτ₂ {sh['dtau2']:.4g} ({sh['kernels']})"
                ratio = f"{ent[0][0] / sh['tau1']:.4f}, {dt / sh['dtau2']:.4f}"
            else:
                pr = ratio = '—'
            L_.append(f"| {r['A'][0]:g}, {r['A'][1]:g} | {r['D']:g} | {', '.join(f'({t:.4g}, {j})' for t, j in ent[:3])} | {dt:.4g} | {pr} | "
                      f"{ratio} | {r3} | {dt * r['D'] ** 2:.4g} | {dt * r['D']:.4g} | {r['n_act']}, {[round(v, 5) for v in r['c_end']]} "
                      f"({[round(v, 5) for v in r['cs']]}), {fmt(r['rms'])} |")
        L_.append('')
    if B is not None:
        L_ += ['## B. The two-kernel block at its exact fit', '',
               'The Gauss–Newton Hessian JᵀJ (J = ∂(Φ(c)w)/∂(w, c)/N^{1/2}) and the fast flow on the critical curve linearised at the fit '
               '(its relaxation rates): eigenvalues, the two smallest scaled by Δ² and Δ⁶ against their leading order (the split mode: '
               'λ₁ = Δ⁶⟨(P⊥e₃)²⟩/8, amplitude-free; the separation: λ₂ = 4Δ²⟨(P⊥e₂)²⟩w₁²w₂²/(w₁² + w₂²); e_n = ∂ⁿ_cφ/n! at the pair, P⊥ '
               'orthogonal to the lower e_m), and the local slopes of the smallest two in log Δ. The flow\'s split and separation modes are '
               'identified by their eigenvectors\' overlap with (1/w′₁, −1/w′₂, Δ/(2w₁), Δ/(2w₂)) and (0, 0, 1/w₁, −1/w₂) in (α, c), and '
               'compared with λ₁·2/(1/μ₁ + 1/μ₂) (μ_k the amplitudes\' mobility, ẇ = −μ∂E/∂w) and, at frozen amplitudes, '
               'Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²).', '',
               '| A₁, A₂ | Δ | Hessian eigenvalues | λ₂/Δ², λ₁/Δ⁶ | predicted | slopes | flow rates | split: r/Δ⁶ (predicted; '
               'the mode\'s overlap; its centre components over Δ/(2w_k) per unit weight moved) | separation: r/Δ² (predicted; overlap) '
               '| slopes |',
               '|---|---|---|---|---|---|---|---|---|---|']
        prev = {}
        for r in B:
            key = tuple(r['A'])
            h, q = r['hess'], r['rates']
            sl = sq = '—'
            sp, se = r.get('split'), r.get('sep')
            fq = (sp[0], se[0]) if sp else (q[0], q[1])
            if key in prev:
                pD, ph, pq = prev[key]
                sl = ', '.join(f"{np.log(h[k] / ph[k]) / np.log(r['D'] / pD):.3f}" for k in (0, 1))
                sq = ', '.join(f"{np.log(fq[k] / pq[k]) / np.log(r['D'] / pD):.3f}" for k in (0, 1))
            prev[key] = (r['D'], h, fq)
            pl, pf = r.get('pred'), r.get('pred_flow')
            pt = f"{pl[1] / r['D'] ** 2:.4g}, {pl[0] / r['D'] ** 6:.4g}" if pl else '—'
            D6, D2 = r['D'] ** 6, r['D'] ** 2
            if sp:
                sc = r.get('split_c')
                fs = (f"{sp[0] / D6:.4g} ({pf[0] / D6:.4g}; {sp[1]:.3f}" + (f"; centres {sc[0]:.4f}, {sc[1]:.4f}" if sc is not None else '')
                      + f") | {se[0] / D2:.4g} ({pf[1] / D2:.4g}; {se[1]:.3f})")
            else:
                fs = f"{q[0] / D6:.4g} | {q[1] / D2:.4g}"
            L_.append(f"| {key[0]:g}, {key[1]:g} | {r['D']:g} | {', '.join(f'{v:.3e}' for v in h)} | {h[1] / D2:.4g}, "
                      f"{h[0] / D6:.4g} | {pt} | {sl} | {', '.join(f'{v:.3e}' for v in q)} | {fs} | {sq} |")
        L_ += ['', "The separation's rate in its two regimes (the exact linearisation −D·JᵀJ at the fit, the mode nearest the separation): "
               "over the frozen-amplitude rate Gm·4Δ²⟨(P_{⊥e₁}e₂)²⟩w₁²w₂²/(w₁² + w₂²) and over the relaxed-amplitude rate Gm·λ_sep "
               "(P⊥ off e₀ and e₁) — the first while the separation outpaces the amplitudes' mode μ_k⟨φ_k²⟩, the second once it is slower:", '',
               '| A₁, A₂ | ' + ' | '.join(f'Δ = {D:g}' for D in DELTAS_SEP) + ' | the amplitudes\' mode min μ_k⟨φ_k²⟩ |',
               '|---|' + '---|' * (len(DELTAS_SEP) + 1)]
        for A1, A2 in [p_ for p_ in PAIRS_A if abs(p_[0] + p_[1]) > 1e-12]:
            sr = [sep_regimes(A1, A2, D) for D in DELTAS_SEP]
            L_.append(f"| {A1:g}, {A2:g} | " + ' | '.join(f"{a:.3f}, {b:.3f}" for a, b, m in sr) + f" | {sr[0][2]:.3g} |")
        rr = [(r['A'], r['sep'][0] / sorted(r['rates'])[1]) for r in B if r['D'] == 0.0125 and r.get('sep') is not None
              and abs(r['A'][0] + r['A'][1]) > 1e-12]
        if rr:
            L_ += ['', "At Δ = 0.0125 the flow's separation rate over its second-smallest rate (the amplitudes' mode): "
                   + ', '.join(f"({a[0]:g}, {a[1]:g}) {x:.1f}" for a, x in rr) + '.']
        L_.append('')
    if C is not None:
        det, d50 = detection(C)
        L_ += ['## C. Number detection under label noise', '',
               f"A = {PAIR_C}; {R_C} draws per hypothesis; H0 the closest one-bump target; the power of Δτ₂ < its H0 quantile "
               f"{ALPHA_C:g} at the true offset (an oracle-calibrated simple null; E calibrates it), in brackets the Neyman–Pearson power "
               "Φ(‖f₁ − f₀‖/s − z₀.₉₅) of the same simple hypotheses (the medians of Δτ₂ under H1 and H0):", '',
               '| Δ | ' + ' | '.join(f's = {s:g}' for s in SIGMAS_C) + ' |', '|---|' + '---|' * len(SIGMAS_C)]
        Lc, Kc, c0c = setup()
        npw = {}
        for D in DELTAS_C:
            f1 = X.bumps_target(Kc, pair_centres(*PAIR_C, D), PAIR_C)
            c1_, w1_ = one_fit(Kc, f1)
            for s_ in SIGMAS_C:
                npw[(D, s_)] = np_power(Kc, f1, w1_ * Kc.phi(np.array([c1_]))[:, 0], s_, ALPHA_C)
        for D in DELTAS_C:
            L_.append(f"| {D:g} | " + ' | '.join(f"{det[(D, s)]['power']:.2f} [{npw[(D, s)]:.2f}] ({det[(D, s)]['med1']:.3g}, "
                                                f"{det[(D, s)]['med0']:.3g})" for s in SIGMAS_C) + ' |')
        ss = np.array([s for s in SIGMAS_C if np.isfinite(d50[s])])
        dd = np.array([d50[s] for s in ss])
        sl = np.polyfit(np.log(ss), np.log(dd), 1)[0] if ss.size >= 2 else np.nan
        # the prediction: under H1 the second entry comes after ≈ C₁/Δ² (C₁ the leading order's Δτ₂Δ² at the smallest Δ, (i)), so the
        # power is ½ where C₁/Δ² meets the H0 quantile q_α (∝ √N/s): Δ₅₀ = (C₁/q_α)^{1/2}, q_α interpolated in log Δ at the crossing
        Dm = min(DELTAS_C)
        C1 = predicted_delay(Kc, X.bumps_target(Kc, pair_centres(*PAIR_C, Dm), PAIR_C), c0c, *PAIR_C, Dm)[0] * Dm ** 2
        pred50 = {}
        for s_ in SIGMAS_C:
            Ds = np.array(sorted(DELTAS_C))
            gap = np.log(C1 / Ds ** 2) - np.log([det[(D, s_)]['q0'] for D in Ds])
            pred50[s_] = np.nan
            for k in range(len(Ds) - 1):
                if gap[k] > 0 >= gap[k + 1] or gap[k] <= 0 < gap[k + 1]:
                    t = gap[k] / (gap[k] - gap[k + 1])
                    pred50[s_] = float(np.exp(np.log(Ds[k]) + t * (np.log(Ds[k + 1]) - np.log(Ds[k]))))
                    break
        qs = {s_: float(np.median([det[(D, s_)]['q0'] for D in DELTAS_C]) * s_) for s_ in SIGMAS_C}
        ci50, cisl = d50_bootstrap(C)
        h0 = [(q['dt'], abs(q['j'][1] - q['j'][0]) <= 1) for r in C if r['hyp'] == 0 for q in r['rows'] if len(q['j']) > 1]
        early = [a for r in C if r['hyp'] == 0 for q, a in [(q, abs(q['j'][1] - q['j'][0]) <= 1) for q in r['rows'] if len(q['j']) > 1]
                 if q['dt'] <= np.quantile([q_['dt'] for q_ in r['rows']], ALPHA_C)]
        npd = {s_: crossing(sorted(DELTAS_C), [npw[(D, s_)] for D in sorted(DELTAS_C)]) for s_ in SIGMAS_C}
        L_ += ['', 'Δ₅₀ (the power ½): ' + ', '.join(f's = {s:g}: {fmt(d50[s])}' for s in SIGMAS_C) +
               f"; the slope of log Δ₅₀ in log s {sl:.3f} (½ predicted). Predicted from (i) and the H0 quantile, Δ₅₀ = (C₁/q_α)^{{1/2}} "
               f"(C₁ = {C1:.4f}, the leading order's Δτ₂Δ² at Δ = {Dm:g}): " + ', '.join(f's = {s:g}: {fmt(pred50[s])}' for s in SIGMAS_C) +
               "; the H0 quantile times s (its median over Δ; constant if q_α ∝ 1/s): " + ', '.join(f'{qs[s]:.4g}' for s in SIGMAS_C) + '. '
               f"The bootstrap's 95% intervals ({N_BOOT} resamples of both hypotheses' draws per cell): Δ₅₀ "
               + ', '.join(f"s = {s:g}: " + ('— (no crossing)' if ci50[s][2] == N_BOOT else f"[{ci50[s][0]:.4g}, {ci50[s][1]:.4g}]"
                                                 + (f" ({ci50[s][2]} resamples without a crossing)" if ci50[s][2] else ''))
                           for s in SIGMAS_C)
               + f"; the slope [{cisl[0]:.3f}, {cisl[1]:.3f}]. The Neyman–Pearson bound's Δ₅₀: "
               + ', '.join(f's = {s:g}: {fmt(npd[s])}' for s in SIGMAS_C) + '. '
               f"Under H0 the second entrant is the first's neighbour in {np.mean([a for t, a in h0]):.1%} of the {len(h0)} draws, in "
               f"{np.mean(early):.1%} of the earliest {ALPHA_C:.0%} of each cell (which set the quantile).", '']
        if S is not None:
            L_ += c_prime_lines(C, S, npw)
        L_ += ['## D. Placement under label noise (the two-kernel block after the second entry, H1)', '',
               f"The first {R_D} H1 draws of each (Δ, s) relaxed after the second entry within {RHS_D} right-hand-side evaluations (a "
               "block near coalescence relaxes at a rate ∝ Δ⁶, B): the placed blocks (two kernels, relaxed) and those out of budget; "
               "the rms errors over the placed ones, against the Cramér–Rao bounds at the truth (s²(JᵀJ)⁻¹, `crb`); the placed blocks "
               "with both kernels within σ of their features, and the draws whose second entrant is the first's neighbour:", '',
               '| Δ | s | placed, out of budget (of the draws) | rms position error | × Δ²/s | / CRB | rms error of w₁/(w₁ + w₂) | × Δ³/s | '
               '/ CRB | both within σ | second entrant adjacent |',
               '|---|---|---|---|---|---|---|---|---|---|---|']
        cs_truth = {D: pair_centres(*PAIR_C, D) for D in DELTAS_C}
        for r in C:
            if r['hyp'] != 1:
                continue
            ok = [q for q in r['rows'] if 'c' in q]
            tried = [q for q in r['rows'] if 'oob' in q]
            adj = sum(len(q['j']) > 1 and abs(q['j'][1] - q['j'][0]) <= 1 for q in r['rows'][:R_D])
            if not ok:
                if tried:
                    L_.append(f"| {r['D']:g} | {r['s']:g} | 0, {sum(q['oob'] for q in tried)} ({len(tried)}) | — | — | — | — | — | — | — | "
                              f"{adj} |")
                continue
            pe = np.sqrt(np.mean([np.mean((np.array(q['c']) - cs_truth[r['D']]) ** 2) for q in ok]))
            fr = np.sqrt(np.mean([(q['w'][0] / (q['w'][0] + q['w'][1]) - PAIR_C[0] / sum(PAIR_C)) ** 2 for q in ok]))
            cp, cf_ = crb(Kc, cs_truth[r['D']], np.array(PAIR_C), r['s'])
            both = sum(bool(np.all(np.abs(np.array(q['c']) - cs_truth[r['D']]) < X.SIG_C)) for q in ok)
            L_.append(f"| {r['D']:g} | {r['s']:g} | {len(ok)}, {sum(q['oob'] for q in tried)} ({len(tried)}) | {pe:.3e} | "
                      f"{pe * r['D'] ** 2 / r['s']:.3g} | {pe / cp:.3g} | {fr:.3e} | {fr * r['D'] ** 3 / r['s']:.3g} | {fr / cf_:.3g} | "
                      f"{both} | {adj} |")
        L_.append('')
        D0 = min(DELTAS_C)
        cp, cf_ = crb(Kc, pair_centres(*PAIR_C, D0), np.array(PAIR_C), 1.0)
        good, big = [], None
        for r in C:
            ok = [q for q in r['rows'] if 'c' in q] if r['hyp'] == 1 else []
            if not ok:
                continue
            cs_ = cs_truth[r['D']]
            cp_, cf2 = crb(Kc, cs_, np.array(PAIR_C), r['s'])
            if cf2 <= 0.2:
                pe = np.sqrt(np.mean([np.mean((np.array(q['c']) - cs_) ** 2) for q in ok]))
                fr = np.sqrt(np.mean([(q['w'][0] / (q['w'][0] + q['w'][1]) - PAIR_C[0] / sum(PAIR_C)) ** 2 for q in ok]))
                good.append((pe / cp_, fr / cf2))
            for q in ok:
                if big is None or max(np.abs(q['w'])) > max(np.abs(big[2]['w'])):
                    big = (r['D'], r['s'], q)
        cells = []
        for r in C:
            ok = [q for q in r['rows'] if 'c' in q] if r['hyp'] == 1 else []
            if ok:
                cs_ = cs_truth[r['D']]
                pe = np.sqrt(np.mean([np.mean((np.array(q['c']) - cs_) ** 2) for q in ok]))
                fr = np.sqrt(np.mean([(q['w'][0] / (q['w'][0] + q['w'][1]) - PAIR_C[0] / sum(PAIR_C)) ** 2 for q in ok]))
                cells.append((r['s'] / r['D'] ** 2, r['s'] / r['D'] ** 3, pe, fr / crb(Kc, cs_, np.array(PAIR_C), r['s'])[1]))
        jumped = [c for c in cells if c[2] >= 0.2]
        kept = [c for c in cells if c[2] < 0.2]
        spl = sorted((c[1], c[3]) for c in cells)
        L_ += [f"The positions' rms error against s/Δ²: at least 0.2 in {len(jumped)} cells, all with s/Δ² ≥ "
               f"{min(c[0] for c in jumped) if jumped else np.nan:.3g} ({min(c[2] for c in jumped) if jumped else np.nan:.3f}–"
               f"{max(c[2] for c in jumped) if jumped else np.nan:.3f}), below it in {len(kept)} cells, all with s/Δ² ≤ "
               f"{max(c[0] for c in kept) if kept else np.nan:.3g} (at most {max(c[2] for c in kept) if kept else np.nan:.3f}). The split's "
               "rms error over its bound against s/Δ³: " + ', '.join(f"{a:.3g}: {b:.2f}" for a, b in spl) + '.']
        L_ += [f"The Cramér–Rao constants at Δ = {D0:g}: the positions' rms × Δ²/s {cp * D0 ** 2:.4f}, the split's × Δ³/s "
               f"{cf_ * D0 ** 3:.4f} — its std reaches 0.1 at s ≈ {0.1 / (cf_ * D0 ** 3):.1f}Δ³. In the {len(good)} cells where the split's "
               f"bound is at most 0.2 the measured errors are {min(g[0] for g in good):.2f}–{max(g[0] for g in good):.2f} (positions) and "
               f"{min(g[1] for g in good):.2f}–{max(g[1] for g in good):.2f} (split) of the bounds. The largest coefficient of a placed "
               f"block: Δ = {big[0]:g}, s = {big[1]:g}, w = ({big[2]['w'][0]:.4g}, {big[2]['w'][1]:.4g}) (their sum "
               f"{sum(big[2]['w']):.4g}), the centres {abs(big[2]['c'][1] - big[2]['c'][0]):.2g} apart, the entrants {big[2]['j']}, its own "
               f"split error {abs(big[2]['w'][0] / sum(big[2]['w']) - PAIR_C[0] / sum(PAIR_C)):.4g}.", '', smith_line(Kc), '']
    if E_ is not None:
        L_ += ['## E. The null law and its calibration across the grid offset', '',
               f"At the centroids c̄ = {', '.join(f'{c:g}' for c in OFFSETS_E)} (the kernels on a grid of spacing 0.1: c̄ = 0.05 is a tie) and "
               f"the cells (Δ, s) = {', '.join(f'({D:g}, {s:g})' for D, s in CELLS_E)}, {R0_E} draws under H0 (the one-bump fit) and "
               f"{R1_E} under H1. The H0 quantile of Δτ₂ at {ALPHA_C:g} times s from the draws and from the null law at the truth "
               f"(`knot_law`, {M_LAW} samples: min_j (1 − sign(n_j)u_j)/|n_j|, n ~ N(0, s²ΦᵀRΦ/N²), R the projection off the fit's "
               "tangent span) and that quantile's false-alarm rate on the draws (Clopper–Pearson 95%); the projected noise's std over the "
               "unprojected one, (⟨(Rφ_j)²⟩/⟨φ_j²⟩)^{1/2}, for the kernel j₁ nearest the fit and its neighbours j₁ ∓ 1; the share of the H0 second entrants next to the first, over all draws and over "
               f"the earliest {ALPHA_C:g} (the draws below the H0 quantile); the false-alarm rate of the threshold set at "
               f"c̄ = {CBAR:g} and of the test calibrated per draw by the law bootstrapped at its own first fit ({M_E} samples); the power of "
               "the oracle-calibrated test (the offset's own H0 quantile), of the calibrated test, and the Neyman–Pearson bound:", '',
               '| Δ, s | c̄ | noiseless Δτ₂Δ² (leading order; kernel) | H0 quantile × s: draws, law | the law\'s false alarms [CP] | '
               f'std ratio: j₁; j₁ ∓ 1 | H0 second entrants adjacent: all, earliest | false alarms: c̄ = {CBAR:g} threshold, calibrated [CP] | power: '
               'oracle, calibrated, Neyman–Pearson |',
               '|---|---|---|---|---|---|---|---|---|']
        Le, Ke, c0e = setup()
        summ = []                                                      # per cell: the law over the draws, the false alarms, the powers
        for D, s_ in CELLS_E:
            es = sorted([e for e in E_ if e['D'] == D and e['s'] == s_], key=lambda e: e['cbar'])
            ref = next((e for e in es if abs(e['cbar'] - CBAR) < 1e-12), None)
            qref = np.quantile([q['dt'] for q in ref['rows0']], ALPHA_C) if ref is not None else np.nan
            for e in es:
                d0 = np.array([q['dt'] for q in e['rows0']]); d1 = np.array([q['dt'] for q in e['rows1']])
                p0 = np.array([q['p'] for q in e['rows0']]); p1 = np.array([q['p'] for q in e['rows1']])
                q0 = np.quantile(d0, ALPHA_C)
                summ.append(dict(cbar=e['cbar'], law=e['q_law'] / q0, fa_fix=np.mean(d0 < qref), fa_cal=np.mean(p0 < ALPHA_C),
                                 pw_or=np.mean(d1 < q0), pw_cal=np.mean(p1 < ALPHA_C), np=e['np'], k_law=int(np.sum(d0 < e['q_law'])),
                                 k_cal=int(np.sum(p0 < ALPHA_C)), n0=d0.size))
                cg = np.linspace(-1, 1, len(e['ratio']))
                j1 = int(np.argmin(np.abs(cg - e['c1'])))
                cl, cc = cp_int(summ[-1]['k_law'], d0.size), cp_int(summ[-1]['k_cal'], d0.size)
                adj = np.mean([abs(q['j'][1] - q['j'][0]) <= 1 for q in e['rows0'] if len(q['j']) > 1])
                adj5 = np.mean([abs(q['j'][1] - q['j'][0]) <= 1 for q in e['rows0'] if len(q['j']) > 1 and q['dt'] < q0])
                pd = predicted_delay(Ke, X.bumps_target(Ke, pair_centres(*PAIR_C, D, e['cbar']), PAIR_C), c0e, *PAIR_C, D,
                                     cbar=e['cbar'])
                L_.append(f"| {D:g}, {s_:g} | {e['cbar']:g} | {pd[0] * D ** 2:.3f} ({pd[1]}) | {q0 * s_:.3f}, {e['q_law'] * s_:.3f} | "
                          f"{np.mean(d0 < e['q_law']):.3f} [{cl[0]:.3f}, {cl[1]:.3f}] | {e['ratio'][j1]:.3f}; "
                          + ', '.join(f"{e['ratio'][k]:.3f}" for k in (j1 - 1, j1 + 1) if 0 <= k < len(e['ratio']))
                          + f" | {adj:.2f}, {adj5:.2f} | {np.mean(d0 < qref):.3f}, {np.mean(p0 < ALPHA_C):.3f} [{cc[0]:.3f}, {cc[1]:.3f}] | "
                          f"{np.mean(d1 < q0):.3f}, {np.mean(p1 < ALPHA_C):.3f}, {e['np']:.3f} |")
        if summ:
            law = np.array([q['law'] for q in summ])
            near = [q for q in summ if q['cbar'] > CBAR + 1e-12]
            far_ = [q for q in summ if q['cbar'] < CBAR - 1e-12]                # out of sample: the threshold is the CBAR cells' own
            kl, kc, n_ = sum(q['k_law'] for q in summ), sum(q['k_cal'] for q in summ), sum(q['n0'] for q in summ)
            dp = np.array([q['pw_cal'] - q['pw_or'] for q in summ])
            L_ += ['', f"Over the {len(summ)} cells: the law's H0 quantile over the draws' {law.min():.3f}–{law.max():.3f} (median "
                       f"{np.median(law):.3f}); the false alarms of the threshold fixed at c̄ = {CBAR:g} "
                       + (f"{min(q['fa_fix'] for q in far_):.3f}–{max(q['fa_fix'] for q in far_):.3f}" if far_ else '—')
                       + f" at c̄ < {CBAR:g} (out of sample) and "
                       f"{min(q['fa_fix'] for q in near):.3f}–{max(q['fa_fix'] for q in near):.3f} beyond, of the calibrated test "
                       f"{min(q['fa_cal'] for q in summ):.3f}–{max(q['fa_cal'] for q in summ):.3f} (nominal {ALPHA_C:g}); the calibrated "
                       f"test's power minus the oracle's {dp.min():+.3f} … {dp.max():+.3f}, within 0.05 in {int(np.sum(np.abs(dp) <= 0.05 + 1e-9))} "
                       f"of {len(summ)}; the Neyman–Pearson bound minus the larger of the two powers "
                       f"{min(q['np'] - max(q['pw_or'], q['pw_cal']) for q in summ):.3f}–{max(q['np'] - max(q['pw_or'], q['pw_cal']) for q in summ):.3f}. "
                       f"Pooled over the {n_} H0 draws: the law at the truth's false-alarm rate {kl / n_:.4f} [{cp_int(kl, n_)[0]:.4f}, "
                       f"{cp_int(kl, n_)[1]:.4f}], the calibrated test's {kc / n_:.4f} [{cp_int(kc, n_)[0]:.4f}, {cp_int(kc, n_)[1]:.4f}]."]
        L_.append('')
        if S is not None:
            L_ += e_prime_lines(E_, S, U)
    if F is not None:
        L_ += f_lines(F, E_)
        L_ += margin_lines(F, E_, G)
    L_.append('Numbers are produced by `experiments/exp_close_features.py`; nothing is edited by hand.')
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'close_features_raw.npz')
    md = os.path.join(OUT, 'close_features.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'E', 'S', 'F', 'U', 'G')
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else {'A', 'B', 'C', 'D', 'E', 'S', 'F', 'U', 'G'}
    res = dict(old)
    t0 = time.time()
    for k, fn in (('B', part_B), ('A', part_A), ('C', part_C), ('E', part_E), ('S', lambda: part_S(res['C'], res['E'])), ('F', part_F),
                  ('U', lambda: part_U(res['E'])), ('G', part_G)):
        if k in parts or (k == 'C' and 'D' in parts):
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
