"""
Proposition D‴.36: the gauge under dither — the replicator's masses diffuse along the manifold of exact fits; a leak on the masses pins them
=========================================================================================================================================

Proposition D‴.32's three-bump instance with one particle per class at the exact representation ϑ^e(1/3) (m = 3, the compact cohesion
inactive between the far classes), the replicator weight law of Proposition D‴.34 at rate γ, the positions dithered √(2τG) per particle
(Euler–Maruyama, dt = 0.25), and optionally the entropic leak on the masses, energy (ε/γ)Σβ log β.

A. The prediction: the manifold M of exact fits {W_ka(α_k) = A_k, the centres and widths the bumps'} through the configuration, its tangent T
   (the K − 1 gauges, δα_k = −A_k δW_k/(W_k²a′_k)), the flow's metric g (G⁻¹ on the positions, (m/γ)/W_k on the masses) and the position
   noise Σ = 2τG: the masses' diffusion D = (TᵀgT)⁻¹Tᵀ g Σ g T(TᵀgT)⁻¹ per unit time (the g-orthogonal projection of the noise onto M: the
   linearised gradient flow's limit map), at W = 1/3 for several γ.
B. The diffusion measured: 64 seeds, τ = 10⁻⁴, γ = 10⁻², 5·10⁴ slow units; the masses' displacement variance and its slope over
   s ∈ [10⁴, 5·10⁴] (a bootstrap error), against A; the masses' mean at the end (the limit's drift term, which A does not include).
C. No leak over long times (10 seeds, 2·10⁵ slow units): the masses' spread, the amplitude coordinates' excursions toward the saturation,
   the fit error; γ = 0 (equal weights, no gauge) as the control (4 seeds).
D. The leak on the masses (ε = 10⁻⁴, 10⁻³; 4 seeds each, 2·10⁵ slow units): the masses' spread and mean over the second half (the
   stationary point at W = w, no bias where w is representable), against the leak's Ornstein–Uhlenbeck limit on M (`leak_spread`), the
   fit error.

Output: results/chaos/gauge_dither.md and gauge_dither_raw.npz (GD_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`).
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
from src.core import mass_transfer as MT
import exp_cohesion_leak as E

OUT = os.environ.get('GD_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GD_WORKERS', 3))
TAU, GAMMA, DT = 1e-4, 1e-2, 0.25
SEEDS_B, S_B = tuple(range(100, 164)), 50000.0
SEEDS_C, S_C = tuple(range(10)), 200000.0
EPS_D, SEEDS_D = (1e-4, 1e-3), (0, 1, 2, 3)
REC = 4000                                                                    # steps between records (10³ slow units)


def gauge_diffusion(W, gamma, tau, A=None, a_max=4.0):
    """D (K × K, rows summing to zero): the per-unit-time covariance of the masses' displacement along M at the masses W."""
    A = E.A_BUMP if A is None else np.asarray(A, float)
    K = len(A); m = K; Ga = E.G[0]
    W = np.asarray(W, float)
    a = A / W; ap = 1 - (a / a_max) ** 2
    T = np.zeros((2 * K, K - 1))
    for j in range(K - 1):
        dW = np.zeros(K); dW[j] = 1.0; dW[K - 1] = -1.0
        T[:K, j] = -(A / (W ** 2 * ap)) * dW
        T[K:, j] = dW
    g = np.diag(np.r_[np.full(K, 1 / Ga), (m / gamma) / W])
    Sig = np.diag(np.r_[np.full(K, 2 * tau * Ga), np.zeros(K)])
    M = np.linalg.solve(T.T @ g @ T, T.T @ g)
    C = M @ Sig @ M.T
    D = np.zeros((K, K)); D[:K - 1, :K - 1] = C
    D[K - 1, :K - 1] = D[:K - 1, K - 1] = -C.sum(axis=1); D[K - 1, K - 1] = C.sum()
    return D


def leak_spread(W, gamma, tau, eps, A=None, a_max=4.0):
    """The stationary standard deviations of the masses under the entropic leak, linearised at W (the Ornstein–Uhlenbeck limit on M): in the
    coordinates δW = Σ_j c_j(e_j − e_K) the leak's stiffness κ = (TᵀgT)⁻¹ H, H = Tᵀ diag(0, εm/(γW)) T (the Hessian of (ε/γ)Σβ log β), the
    covariance C from κC + Cκᵀ = (TᵀgT)⁻¹TᵀgΣgT(TᵀgT)⁻¹."""
    from scipy.linalg import solve_continuous_lyapunov
    A = E.A_BUMP if A is None else np.asarray(A, float)
    K = len(A); m = K; Ga = E.G[0]
    W = np.asarray(W, float)
    a = A / W; ap = 1 - (a / a_max) ** 2
    T = np.zeros((2 * K, K - 1))
    for j in range(K - 1):
        dW = np.zeros(K); dW[j] = 1.0; dW[K - 1] = -1.0
        T[:K, j] = -(A / (W ** 2 * ap)) * dW
        T[K:, j] = dW
    g = np.diag(np.r_[np.full(K, 1 / Ga), (m / gamma) / W])
    Sig = np.diag(np.r_[np.full(K, 2 * tau * Ga), np.zeros(K)])
    TgT = T.T @ g @ T
    M = np.linalg.solve(TgT, T.T @ g)
    kappa = np.linalg.solve(TgT, T.T @ np.diag(np.r_[np.zeros(K), eps * m / (gamma * W)]) @ T)
    C = solve_continuous_lyapunov(kappa, M @ Sig @ M.T)
    Tw = T[K:]
    return np.sqrt(np.diag(Tw @ C @ Tw.T))


def part_A():
    return dict(rows=[dict(gamma=g, D=gauge_diffusion(np.ones(3) / 3, g, TAU)) for g in (1e-3, 1e-2, 1e-1, 1.0, 10.0)])


def series(args):
    gamma, eps, seed, S = args
    t0 = time.time()
    net, base = E.setup()
    Th = E.exact_rep(net, np.ones(3) / 3).copy()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, eps=eps)
    rec = []
    loop.run_dithered(Th, np.ones(3), S, DT, TAU, np.random.default_rng(seed), rec_every=REC,
                      callback=lambda s_, T, b: rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b)]))
    print(f"[{time.time() - t0:.0f}s] gamma {gamma:g} eps {eps:g} seed {seed}", flush=True)
    return dict(gamma=gamma, eps=eps, seed=seed, rec=np.array(rec))


def msd(runs):
    R = np.array([r['rec'] for r in runs])
    W = R[:, :, 4:7]
    return R[0, :, 0], np.mean((W - W[:, :1, :]) ** 2, axis=0), W


def report(A=None, B=None, C=None, D=None):
    fmt = E.fmt
    L_ = ["# Proposition D‴.36: the gauge under dither (`exp_gauge_dither.py`)", "",
          f"Proposition D‴.32's instance, one particle per class at ϑ^e(1/3) (m = 3), the replicator at rate γ, the positions dithered √(2τG), τ = {TAU:g}, "
          f"Euler–Maruyama dt = {DT:g}; the masses W_k = β_k/3.", ""]
    if A is not None:
        L_ += ["## A. The predicted diffusion of the masses along the manifold of exact fits (W = 1/3, per unit slow time)", "",
               "| γ | D_11, D_22, D_33 | D_12, D_13, D_23 |", "|---|---|---|"]
        for r in A['rows']:
            D_ = r['D']
            L_.append(f"| {r['gamma']:g} | {', '.join(fmt(D_[k, k]) for k in range(3))} | {fmt(D_[0, 1])}, {fmt(D_[0, 2])}, {fmt(D_[1, 2])} |")
        L_.append("")
    if B is not None:
        s, ms, W = msd(B)
        dvar = lambda W_: np.var(W_ - W_[:, :1, :], axis=0, ddof=1)      # round 26 (D3y-5): the variance, not the mean square (the drift removed)
        vs = dvar(W)
        i0 = int(np.argmin(np.abs(s - 1e4))); i1 = len(s) - 1
        slope = (vs[i1] - vs[i0]) / (s[i1] - s[i0])
        rng = np.random.default_rng(0); bs = []
        for _ in range(1000):
            idx = rng.integers(0, len(W), len(W)); v_ = dvar(W[idx])
            bs.append((v_[i1] - v_[i0]) / (s[i1] - s[i0]))
        err = np.std(bs, axis=0)
        pred = np.diag(A['rows'][[r['gamma'] for r in A['rows']].index(GAMMA)]['D']) if A is not None else np.full(3, np.nan)
        L_ += [f"## B. The diffusion measured ({len(W)} seeds, γ = {GAMMA:g}, {S_B:g} slow units)", "",
               "| s | mean-square displacement of W_1, W_2, W_3 | the displacement's variance |", "|---|---|---|"]
        for t in (5e3, 1e4, 2e4, 3e4, 5e4):
            i = int(np.argmin(np.abs(s - t)))
            L_.append(f"| {s[i]:g} | {', '.join(fmt(v) for v in ms[i])} | {', '.join(fmt(v) for v in vs[i])} |")
        L_ += ["", f"The variance's slope over s ∈ [{s[i0]:g}, {s[i1]:g}]: " + ', '.join(f"{fmt(v)} ± {fmt(e)}" for v, e in zip(slope, err))
               + f" (bootstrap standard errors); predicted (A): {', '.join(fmt(v) for v in pred)}; ratios " + ', '.join(f"{v / p:.2f}" for v, p in zip(slope, pred)) + ".", ""]
        Wend = W[:, -1, :]
        L_ += [f"The masses' mean at s = {s[-1]:g} (from 1/3): " + ', '.join(f"{v:.4f} ± {e:.4f}" for v, e in zip(Wend.mean(axis=0), Wend.std(axis=0, ddof=1) / np.sqrt(len(Wend))))
               + " (standard errors): the drift term of the limit, which (A) does not include.", ""]
    if C is not None or D is not None:
        L_ += ["## C–D. Long runs: no leak, the control γ = 0, and the leak on the masses (2·10⁵ slow units)", "",
               "| γ | ε | seeds | the masses over the second half: mean | std | predicted std (the leak's Ornstein–Uhlenbeck limit at W = 1/3) | the amplitude coordinates' range (α_1, α_2, α_3) | the largest |a|/a_max | fit error (second-half mean) |",
               "|---|---|---|---|---|---|---|---|---|"]
        groups = {}
        for r in list(C or []) + list(D or []):
            groups.setdefault((r['gamma'], r['eps']), []).append(r)
        for (g, e), rs in sorted(groups.items(), key=lambda kv: (kv[0][0] == 0, kv[0][1])):
            R = np.array([r['rec'] for r in rs]); h = R.shape[1] // 2
            W = R[:, h:, 4:7]; al = R[:, :, 1:4]
            amax = np.abs(np.tanh(al / 4.0)).max()
            L_.append(f"| {g:g} | {e:g} | {len(rs)} | {', '.join(f'{v:.4f}' for v in W.mean(axis=(0, 1)))} | {', '.join(f'{v:.4f}' for v in W.std(axis=(0, 1)))} | "
                      + (', '.join(f'{v:.4f}' for v in leak_spread(np.ones(3) / 3, g, TAU, e)) if e > 0 else '—') + " | "
                      + ', '.join(f"[{al[:, :, k].min():+.2f}, {al[:, :, k].max():+.2f}]" for k in range(3)) + f" | {amax:.4f} | {fmt(R[:, h:, 7].mean())} |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_gauge_dither.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_dither_raw.npz')
    md = os.path.join(OUT, 'gauge_dither.md')
    old = E.load(raw)
    get = lambda d, k: (list(d[k]) if k in ('B', 'C', 'D') else d[k]) if k in d else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(get(old, 'A'), get(old, 'B'), get(old, 'C'), get(old, 'D')))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCD')
    from multiprocessing import Pool
    res = dict(old)
    if 'A' in parts:
        res['A'] = part_A(); E.save(raw, **res)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'B' in parts:
            h['B'] = pool.map_async(series, [(GAMMA, 0.0, s, S_B) for s in SEEDS_B], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(series, [(GAMMA, 0.0, s, S_C) for s in SEEDS_C] + [(0.0, 0.0, s, S_C) for s in SEEDS_D], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(series, [(GAMMA, e, s, S_C) for e in EPS_D for s in SEEDS_D], chunksize=1)
        for k in 'BCD':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res)
    open(md, 'w').write(report(get(res, 'A'), get(res, 'B'), get(res, 'C'), get(res, 'D')))


if __name__ == '__main__':
    main()
