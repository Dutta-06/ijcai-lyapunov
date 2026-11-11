"""
Proposition D‴.39: dissipativity of the dithered gauge — the escape to saturated approximate fits, and a dead-zone confinement
==============================================================================================================================

Proposition D‴.36's instance (Proposition D‴.32's three bumps, one particle per class, m = 3, the cohesion leak), the positions dithered at
τ. Under the cohesion leak J is bounded on the unconstrained coordinates, so the positions' Gibbs law e^{−J/τ} is not normalisable; the
dither leaves the exact fit's well for approximate fits with saturated amplitudes (the plateau), on which the saturated coordinates diffuse
freely. The dead zone V = (σ_B/2)Σ_i‖ϑ_i − Π_B(ϑ_i)‖²_{G⁻¹} (`ReplicatorCohesionLoop(box=L, box_sigma=σ_B)`, B = Π_j[−L_j, L_j] in the
unconstrained coordinates (α, γ, κ)) makes it normalisable and leaves everything inside B unchanged.

A. The plateau: at the masses 1/3 and W♯ (Proposition D‴.38), for each class k the lowest configuration with α_k fixed at ±20 (the sign of
   A_k; L-BFGS over the other eight coordinates from N_A perturbed exact fits), its level ΔJ/τ above the exact fit (τ = 10⁻⁴) and its fit
   error; with α_k at ±ESC (the escape set's boundary, Part C's ΔJ); and with the box, α_k at its edge ±L_α.
B. The minimum-energy path from the exact fit to the plateau (the string method in the G⁻¹ metric, N_B images, the endpoint the lowest
   plateau of A with its saturated amplitude set to 10 and the rest relaxed): the profile ΔJ/τ along it.
C. The escapes at fixed masses (γ = 0, W = 1/3, τ = 1.5, 2, 3·10⁻⁴, 24 seeds, at most S_C slow units): the first time max_k|α_k| > ESC; the
   escape rate (escapes over the total time of the runs, censored runs included) against e^{−ΔJ/τ}, ΔJ the lowest level of J on the escape
   set's boundary {max_k|α_k| = ESC} (A; Freidlin–Wentzell: τ log of the mean exit time tends to it).
D. The dead zone at fixed masses (γ = 0, W = 1/3, τ = 2·10⁻⁴, 16 seeds, S_D slow units): without and with the box L = BOX (σ_B = 1): the
   fraction of the records with max|α| > ESC and the fit error by quarters, the largest |α|.
F. The two channels to the saturation in the stored runs (no new runs): in exp_gauge_drift.py's Part H (no leak) and in exp_gauge_fixman.py's
   Part F without the correction (the masses held at 1/3 by the leak), the records with some |α_k| > 10 — their mass's distance from the
   saturation |A_k|/a_max and their fit error: at the saturation with the exact fit's error, the gauge's channel (Proposition D‴.37); away from
   it with the plateau's error, the positions' own (ii).
E. The full loop with the box (γ = 10⁻², τ = 10⁻⁴, the 12 seeds of exp_gauge_drift.py's Part H, 10⁶ slow units): κ_F = 0, and the masses
   form at κ_F = 2: the masses at 10⁵, 3·10⁵, 10⁶ against the vertices of the box's representable masses {W_k ≥ |A_k|/(a_max tanh(L_α/a_max))}
   (κ_F = 0) and W♯ (κ_F = 2); κ_F = 0 without the box is exp_gauge_drift.py's Part H.

G. (Referee round 143, D3el-6.) The escapes at fixed masses 1/3 (γ = 0) from a bounded domain with a smooth boundary, bounded in all
   nine coordinates: 𝒪 = {Σ_{k,j}((ϑ_kj − ϑ*_kj)/L_j)^p < 1}, ϑ* the exact fit, L = DOM_L in (α, γ, κ), p = DOM_P (a rounded box). Its
   barrier ΔJ_𝒪 = min of J − J(exact) on ∂𝒪 (L-BFGS over the boundary from N_DOM random directions); the deterministic flow from N_FLOW
   points of ∂𝒪 (the minimiser and random directions) over S_FLOW slow units (it should end at the exact fit). Then the dithered runs at
   the τ of TAUS_G (SEEDS_G, each up to S_G slow units; an exit: the first step after which Σ((ϑ − ϑ*)/L)^p ≥ 1), integrated together by
   src/core/gauge_batch.py (the noise of each seed is run_dithered's), checkpointed per seed in GP_CKPT and stopped at a block boundary
   after GP_WALL seconds (rerun to resume). The escape rate per τ (exits over the total time, censored runs included), and the weighted
   least-squares fits of log(rate) + c·log τ against 1/τ for the prefactors τ^0, τ^{−1/2}, τ^{−1} and with c free: the fitted barrier over
   ΔJ_𝒪 with its standard error, χ², and the mean exit time extrapolated to τ = 10⁻⁴ against 1/(γτ) = 10⁶ (γ = 10⁻²).

H. (Referee round 150, D3es-1, -4, -9, -10; no new dithered runs.) From Part G's stored exit times: the prefactor with the barrier fixed
   at ΔJ_𝒪 (the weighted fit of log(rate) + ΔJ_𝒪/τ = a − c·log τ: ĉ, its standard error, χ² and the fixed prefactors c = 0, ½, 1 with
   their χ², p-values and z); the lag (the first LAGS_H time units of every run discarded: B/ΔJ_𝒪 for c = 0, ½, 1 and ĉ); the
   extrapolated mean exit times at τ = 10⁻⁴ with the standard error of the fitted log rate there (Part G's fits, the fit with the barrier
   fixed, and Part C's fit with τ^{−1/2}). The local minima of J on ∂𝒪: domain_job's N_DOM L-BFGS runs (the same starts), all kept,
   clustered (CLUSTER_H), and the deterministic flow from each over S_FLOW (the largest Σ((ϑ − ϑ*)/L)^p after the first step, the final
   distance and energy). The exit points: Part G's runs keep their state at the exit step in the checkpoint (GP_CKPT/plateau_G.pkl;
   read with `--only H`, checked against the stored exit times and stored as G_exits); by τ, the coordinate whose term dominates
   Σ((ϑ − ϑ*)/L)^p at the exit, and the second width's sign there (broadened: the side of the lowest minimum; narrowed: the other).

Output: results/chaos/gauge_plateau.md and gauge_plateau_raw.npz (GP_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`;
`--only G` with GP_CKPT (a directory for the checkpoints) and GP_WALL (seconds) runs Part G resumably; `--only H` computes H from the
raw, and its exit points from GP_CKPT's checkpoint when it is there).
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
from src.core.runcache import cached
import exp_cohesion_leak as E
import exp_gauge_dither as GD
import exp_gauge_drift as X
import exp_gauge_fixman as Q

OUT = os.environ.get('GP_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GP_WORKERS', 4))
TAU, A_MAX, DT, REC = GD.TAU, 4.0, GD.DT, GD.REC
BOX, SIGMA_B = (5.0, 3.0, 4.0), 1.0                                         # L_α, L_γ, L_κ (|a| ≤ 3.39, |c| ≤ 2.28, s ∈ [0.33, 1.97])
AL_SAT, N_A = 20.0, 200
N_B, IT_B, DT_B, AL_B = 40, 30000, 0.5, 10.0
ESC = 6.0                                                                    # |a|/a_max > 0.905; in the well |α_k| < 4.5
S2_PLATEAU = 0.6                                                             # the second width on the plateaus' family (exact fit 0.4)
TAUS_C, SEEDS_C, S_C, REC_C = (1.5e-4, 2e-4, 3e-4), tuple(range(500, 524)), 4e5, 40
TAU_D, SEEDS_D, S_D = 2e-4, tuple(range(600, 616)), 4e5
GAMMA_E, SEEDS_E, S_E, T_E = 1e-2, X.SEEDS_H, 1e6, X.T_H
W_MIN = np.abs(E.A_BUMP) / A_MAX
W_BOX = np.abs(E.A_BUMP) / (A_MAX * np.tanh(BOX[0] / A_MAX))
DOM_L, DOM_P, N_DOM, N_FLOW, S_FLOW = (2.5, 0.8, 1.2), 4, 200, 64, 4e4      # Part G (referee round 143)
TAUS_G = (2.5e-4, 2e-4, 1.6e-4, 1.4e-4, 1.25e-4, 1.1e-4)                      # ΔJ_𝒪/τ = 2.5 … 5.7 (to 8 is beyond the budget)
SEEDS_G = tuple(range(700, 748))
S_G = (6e4, 1e5, 1.5e5, 2e5, 2.5e5, 3e5)                                       # the caps, per τ
BLOCK_G, GAMMA_REF = 4000, 1e-2
LAGS_H = (1650.0, 3300.0, 6600.0)                                             # ½, 1, 2 times 1/(3.06·10⁻⁴), the slowest relaxation time
C_GRID_H = (0.0, 0.5, 1.0)
CLUSTER_H = (1e-3, 0.05)                    # one minimum: (J − J(exact))/τ within 10⁻³, (ϑ − ϑ*)/L within 0.05 in every coordinate
CKPT = os.environ.get('GP_CKPT')
WALL = float(os.environ.get('GP_WALL', 'inf'))


def masses(which):
    return np.ones(3) / 3 if which == 'uniform' else Q.w_sharp()


@cached
def plateau_job(args):
    """The lowest J with α_k fixed at al (with the box: J + V), over N_A L-BFGS runs from the exact fit perturbed (seeded)."""
    from scipy.optimize import minimize
    which, k, at = args
    boxed = at == 'box'
    net, base = E.setup()
    W = masses(which); beta = 3 * W
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2, box=(BOX if boxed else None), box_sigma=SIGMA_B)
    Th0 = E.exact_rep(net, W); J0 = loop.J(Th0, beta)
    al = np.sign(E.A_BUMP[k]) * dict(sat=AL_SAT, esc=ESC, box=BOX[0])[at]
    f = lambda x: loop.J(np.insert(x, 3 * k, al).reshape(3, 3), beta)
    rng = np.random.default_rng(1000 + 10 * k + ('sat', 'esc', 'box').index(at))
    best = None
    for _ in range(N_A):
        x0 = np.delete(Th0.ravel(), 3 * k) + rng.normal(0, 1.5, 8) * rng.uniform(0.2, 1.5)
        r = minimize(f, x0, method='L-BFGS-B', options=dict(maxiter=5000, ftol=1e-15, gtol=1e-12))
        if best is None or r.fun < best.fun:
            best = r
    T = np.insert(best.x, 3 * k, al).reshape(3, 3)
    a, c, s = net.physical(T)
    return dict(which=which, k=k, at=at, dJ=(best.fun - J0) / TAU, fit=loop.fit_error(T, beta), a=a[:, 0], c=c[:, 0], s=s, Th=T)


@cached
def string_job(args):
    """The string method (the simplified string method of E–Ren–Vanden-Eijnden) in the G⁻¹ metric between the exact fit and the plateau's
    point `end`, the end relaxed with its saturated amplitude fixed at AL_B: the profile ΔJ/τ and the images."""
    from scipy.optimize import minimize
    from scipy.interpolate import interp1d
    which, k, end = args
    net, base = E.setup()
    W = masses(which); beta = 3 * W; G = E.G
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    A = E.exact_rep(net, W); J0 = loop.J(A, beta)
    f = lambda x: loop.J(np.insert(x, 3 * k, np.sign(E.A_BUMP[k]) * AL_B).reshape(3, 3), beta)
    r = minimize(f, np.delete(np.asarray(end, float).ravel(), 3 * k), method='L-BFGS-B', options=dict(maxiter=20000, ftol=1e-15, gtol=1e-12))
    B = np.insert(r.x, 3 * k, np.sign(E.A_BUMP[k]) * AL_B).reshape(3, 3)
    S = np.array([A + (B - A) * t for t in np.linspace(0, 1, N_B)])
    metric = np.sqrt(1 / G)[None, :]

    def repar(S):
        d = np.r_[0, np.cumsum(np.sqrt((((S[1:] - S[:-1]) * metric) ** 2).sum(axis=(1, 2))))]
        return interp1d(d / d[-1], S.reshape(N_B, -1), axis=0)(np.linspace(0, 1, N_B)).reshape(S.shape)
    for _ in range(IT_B):
        for i in range(1, N_B - 1):
            S[i] = S[i] + DT_B * loop.drift(S[i], beta)[0]
        S = repar(S)
    prof = np.array([(loop.J(x, beta) - J0) / TAU for x in S])
    fits = np.array([loop.fit_error(x, beta) for x in S])
    return dict(which=which, k=k, prof=prof, fits=fits, images=S)


class Escaped(Exception):
    pass


@cached
def escape_job(args):
    """The first record (every REC_C steps) with max|α_k| > ESC at fixed masses 1/3 (γ = 0), or None by S_C."""
    tau, seed = args
    t0 = time.time()
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    out = dict(tau=tau, seed=seed, t=None)

    def cb(s_, T, b):
        if np.abs(T[:, 0]).max() > ESC:
            out.update(t=s_, alpha=T[:, 0].copy(), fit=loop.fit_error(T, b)); raise Escaped
    try:
        loop.run_dithered(E.exact_rep(net, np.ones(3) / 3), np.ones(3), S_C, DT, tau, np.random.default_rng(seed), rec_every=REC_C, callback=cb)
    except Escaped:
        pass
    print(f"[{time.time() - t0:.0f}s] C tau {tau:g} seed {seed}: {out['t']}", flush=True)
    return out


@cached
def series(args):
    """Records (s, α_1..3, W_1..3, the fit error, c_1..3, s_1..3) every REC steps: (γ, τ, κ_F, form, boxed, seed, S)."""
    gamma, tau, kF, form, boxed, seed, S = args
    t0 = time.time()
    net, base = E.setup()
    Th = E.exact_rep(net, np.ones(3) / 3)
    A0 = base.features(Th)[0] / 3 if form == 'masses' else None
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, fixman=kF, fixman_tau=tau, fixman_A=A0,
                                     box=(BOX if boxed else None), box_sigma=SIGMA_B)
    rec = []

    def cb(s_, T, b):
        a, c, s = net.physical(T)
        rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b), c[:, 0], s])
    loop.run_dithered(Th, np.ones(3), S, DT, tau, np.random.default_rng(seed), rec_every=REC, callback=cb)
    print(f"[{time.time() - t0:.0f}s] gamma {gamma:g} tau {tau:g} {form} kF {kF:g} box {boxed} seed {seed}", flush=True)
    return dict(gamma=gamma, tau=tau, kF=kF, form=form, boxed=boxed, seed=seed, rec=np.array(rec))


def part_C(runs, dJ):
    out = []
    for tau in TAUS_C:
        rr = [r for r in runs if r['tau'] == tau]
        t = np.array([r['t'] if r['t'] is not None else S_C for r in rr]); n = sum(r['t'] is not None for r in rr)
        rate = n / t.sum()
        fits = [r['fit'] for r in rr if r['t'] is not None]
        out.append(dict(tau=tau, n=n, N=len(rr), rate=rate, rate_se=rate / np.sqrt(max(n, 1)), med=np.median(t),
                        scaled=rate * np.exp(dJ * TAU / tau), fit=(np.mean(fits) if fits else np.nan)))
    ok = [r for r in out if r['n'] > 0]
    slope = np.polyfit([1 / r['tau'] for r in ok], np.log([r['rate'] for r in ok]), 1)[0] if len(ok) > 1 else np.nan
    return dict(rows=out, slope=slope, dJ=dJ * TAU)


def part_D(runs):
    out = []
    for boxed in (False, True):
        R = np.array([r['rec'] for r in runs if r['boxed'] == boxed]); q = np.array_split(np.arange(R.shape[1]), 4)
        al = np.abs(R[:, :, 1:4]).max(axis=2); pl = R[:, :, 12] > S2_PLATEAU
        out.append(dict(boxed=boxed, frac=[float(np.mean(al[:, ix] > ESC)) for ix in q], plat=[float(np.mean(pl[:, ix])) for ix in q],
                        fit=[float(R[:, ix, 7].mean()) for ix in q], almax=[float(al[:, ix].max()) for ix in q],
                        seeds=int(np.sum(pl.any(axis=1))), N=len(R)))
    return out


def summarise_E(R, ref=None):
    """The masses' statistics; with `ref` (the same seeds without the box: records (s, α, W, the fit error, …)) the first record at which the
    two differ — inside B the dynamics are identical, so this is the box's first action, to within a record."""
    s = R[0, :, 0]; W = R[:, :, 4:7]; al = R[:, :, 1:4]
    out = dict(times=[], half=W[:, R.shape[1] // 2:].mean(axis=1), fit=R[:, R.shape[1] // 2:, 7].mean(), almax=float(np.abs(al).max()),
               end=W[:, -1], split=None)
    if ref is not None:
        ref = np.asarray(ref, float)[:, :R.shape[1], :8]
        diff = np.abs(R[:, :ref.shape[1], :8] - ref).max(axis=2) > 0
        out['split'] = [float(s[np.argmax(d)]) if d.any() else None for d in diff]
    for t in T_E:
        i = int(np.argmin(np.abs(s - t)))
        out['times'].append(dict(s=s[i], mean=W[:, i].mean(axis=0), se=W[:, i].std(axis=0, ddof=1) / np.sqrt(len(W)), std=W[:, i].std(axis=0, ddof=1)))
    return out


def part_F():
    """The records with some |α_k| > 10 in the stored no-leak runs (exp_gauge_drift.py's H) and leak runs without the correction
    (exp_gauge_fixman.py's F, 'none'): per class, their number, how many have W_k − W_min,k < 0.01 and > 0.05, their mean fit error."""
    out = []
    zh = E.load(os.path.join(OUT, 'gauge_drift_raw.npz')).get('H_rec')
    zf = E.load(os.path.join(OUT, 'gauge_fixman_raw.npz')).get('F_rec')
    srcs = [('no leak (exp_gauge_drift.py H)', zh), ('the leak, no correction (exp_gauge_fixman.py F)', None if zf is None else zf.get(('none', 0.0)))]
    for name, R in srcs:
        if R is None:
            continue
        R = np.asarray(R, float)
        al, W, fit = np.abs(R[:, :, 1:4]), R[:, :, 4:7], R[:, :, 7]
        rows = []
        for k in range(3):
            big = al[:, :, k] > 10; d = W[:, :, k] - W_MIN[k]
            rows.append(dict(n=int(big.sum()), near=int((big & (d < 0.01)).sum()), far=int((big & (d > 0.05)).sum()),
                             fit=float(fit[big].mean()) if big.any() else np.nan, seeds=int(big.any(axis=1).sum())))
        out.append(dict(src=name, rows=rows, N=R.shape[0], fit_all=float(fit[:, R.shape[1] // 2:].mean())))
    return out


def vertices():
    return np.array([[1 - W_BOX[1] - W_BOX[2], W_BOX[1], W_BOX[2]], [W_BOX[0], 1 - W_BOX[0] - W_BOX[2], W_BOX[2]],
                     [W_BOX[0], W_BOX[1], 1 - W_BOX[0] - W_BOX[1]]])


def dom_phi(Th, Th0):
    """Σ_{k,j}((ϑ_kj − ϑ*_kj)/L_j)^p for Θ (B, 3, 3): 𝒪 = {dom_phi < 1}."""
    return np.sum(((np.asarray(Th, float) - Th0[None]) / np.asarray(DOM_L)[None, None, :]) ** DOM_P, axis=(1, 2))


def domain_job():
    """Part G's domain: ΔJ_𝒪 = min of J − J(exact) on ∂𝒪 (L-BFGS from N_DOM random directions, the boundary point Θ* + L∘u/‖u‖_p), and the
    deterministic flow (Euler, dt = DT) from N_FLOW boundary points over S_FLOW: the largest final G⁻¹-distance from the exact fit and the
    largest final (J − J(exact))/τ."""
    from scipy.optimize import minimize
    from src.core.gauge_batch import BatchDither
    net, base = E.setup()
    beta = np.ones(3)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    Th0 = E.exact_rep(net, np.ones(3) / 3); J0 = loop.J(Th0, beta)
    Lv = np.tile(np.asarray(DOM_L, float), 3)
    pt = lambda u: (Th0.ravel() + Lv * u / np.sum(np.abs(u) ** DOM_P) ** (1 / DOM_P)).reshape(3, 3)
    f = lambda u: (loop.J(pt(u), beta) - J0) / TAU
    rng = np.random.default_rng(1700)
    best, vals = None, []
    for _ in range(N_DOM):
        r = minimize(f, rng.standard_normal(9), method='L-BFGS-B', options=dict(maxiter=5000, ftol=1e-15, gtol=1e-12))
        vals.append(r.fun)
        if best is None or r.fun < best.fun:
            best = r
    Tb = pt(best.x)
    a, c, sw = net.physical(Tb)
    starts = [Tb] + [pt(rng.standard_normal(9)) for _ in range(N_FLOW - 1)]
    bd = BatchDither(loop)
    T = np.array(starts); lb = np.zeros((len(starts), 3))
    for _ in range(int(round(S_FLOW / DT))):
        T, lb = bd.step(T, lb, DT, np.zeros(3), 0.0)
    dist = np.sqrt(np.sum((T - Th0[None]) ** 2 / E.G[None, None, :], axis=(1, 2)))
    dJ_end = np.array([(loop.J(x, beta) - J0) / TAU for x in T])
    dJ_start = np.array([(loop.J(x, beta) - J0) / TAU for x in starts])
    return dict(dJ=best.fun * TAU, at=Tb - Th0, a=a[:, 0], c=c[:, 0], s=sw, fit=loop.fit_error(Tb, beta),
                n_best=int(np.sum(np.array(vals) < best.fun + 1e-3)), dist_end=dist, dJ_end=dJ_end, dJ_start=dJ_start)


def escape_runs(dom):
    """Part G's dithered runs at fixed masses 1/3 (resumable; see the docstring). Returns the per-run results, or None while some run is
    open (the checkpoint keeps the state)."""
    from src.core.gauge_batch import BatchDither, run_batch
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    runs = [dict(key=(tau, seed), seed=seed, n=int(round(cap / DT)), sq=np.sqrt(2 * tau * DT * E.G), Th0=Th0, lb0=np.zeros(3))
            for tau, cap in zip(TAUS_G, S_G) for seed in SEEDS_G]
    ck = None if CKPT is None else os.path.join(CKPT, 'plateau_G.pkl')
    st = run_batch(BatchDither(loop), runs, DT, BLOCK_G, ckpt=ck, wall=WALL, inside=lambda T: dom_phi(T, Th0) < 1.0,
                   log=lambda x: print(x, flush=True))
    if not all(v['done'] for v in st.values()):
        return None
    return [dict(tau=r['key'][0], seed=r['seed'], cap=r['n'] * DT, t=(None if st[r['key']]['exit'] is None else st[r['key']]['exit'] * DT))
            for r in runs]


def arrhenius(rows, dJ, c=None):
    """The weighted (escape counts) least-squares fit of log(rate) + c·log τ against 1/τ (c None: log(rate) = a − c·log τ − B/τ with c
    free): the barrier B, its standard error, c (and its standard error), χ² = Σn(residual)², the fitted log rate at τ = 10⁻⁴."""
    ok = [r for r in rows if r['n'] > 0]
    x = np.array([1 / r['tau'] for r in ok]); lt = np.log([r['tau'] for r in ok]); y = np.log([r['rate'] for r in ok])
    w = np.array([r['n'] for r in ok], float)
    X = np.c_[np.ones_like(x), -x] if c is not None else np.c_[np.ones_like(x), -x, -lt]
    yy = y + c * lt if c is not None else y
    cov = np.linalg.inv(X.T @ (w[:, None] * X)); beta = cov @ X.T @ (w * yy)
    res = yy - X @ beta
    x4 = np.r_[1.0, -1e4] if c is not None else np.r_[1.0, -1e4, -np.log(1e-4)]
    lr4 = float(x4 @ beta) - (c * np.log(1e-4) if c is not None else 0.0)
    return dict(c=(c if c is not None else float(beta[2])), c_se=(None if c is not None else float(np.sqrt(cov[2, 2]))), B=float(beta[1]),
                B_se=float(np.sqrt(cov[1, 1])), ratio=float(beta[1] / dJ), ratio_se=float(np.sqrt(cov[1, 1]) / dJ), chi2=float(w @ res ** 2),
                dof=len(ok) - X.shape[1], t4=float(np.exp(-lr4)))


def part_G(dom, runs):
    rows = []
    for tau in TAUS_G:
        rr = [r for r in runs if r['tau'] == tau]
        t = np.array([r['t'] if r['t'] is not None else r['cap'] for r in rr]); n = sum(r['t'] is not None for r in rr)
        rate = n / t.sum()
        rows.append(dict(tau=tau, n=n, N=len(rr), cap=rr[0]['cap'], T=float(t.sum()), rate=rate, rate_se=rate / np.sqrt(max(n, 1)),
                         x=dom['dJ'] / tau, med=float(np.median(t)), scaled=rate * np.exp(dom['dJ'] / tau),
                         mean_exit=(float(np.mean([r['t'] for r in rr if r['t'] is not None])) if n else np.nan)))
    fits = [arrhenius(rows, dom['dJ'], c) for c in (0.0, 0.5, 1.0)] + [arrhenius(rows, dom['dJ'])]
    return dict(dom=dom, rows=rows, fits=fits)


def report_G(Gg):
    fmt = E.fmt
    dom = Gg['dom']
    v3 = lambda v, p=2: ', '.join(f'{x:.{p}f}' for x in v)
    L_ = [f"## G. The escapes from a bounded domain with a smooth boundary (referee round 143; fixed masses 1/3, γ = 0, {len(SEEDS_G)} seeds per τ)", "",
          f"The domain 𝒪 = {{Σ_k,j ((ϑ_kj − ϑ*_kj)/L_j)^{DOM_P} < 1}}, ϑ* the exact fit, L = ({v3(DOM_L, 1)}) in (α, γ, κ): bounded in the amplitudes', "
          f"centres' and widths' coordinates, its boundary smooth. The lowest J on ∂𝒪 (L-BFGS from {N_DOM} random directions; {dom['n_best']} reach it "
          f"to 10⁻³τ): ΔJ_𝒪 = {fmt(dom['dJ'])}, ΔJ_𝒪/τ = {dom['dJ'] / TAU:.3f} at τ = {TAU:g}, at the deviation from ϑ* "
          f"{'; '.join(v3(r) for r in dom['at'])} (a = {v3(dom['a'])}, c = {v3(dom['c'])}, s = {v3(dom['s'])}; the fit error {dom['fit']:.4f}).", "",
          f"The deterministic flow from {len(dom['dist_end'])} points of ∂𝒪 (the minimiser and random directions; (J − J(exact))/τ between "
          f"{dom['dJ_start'].min():.2f} and {dom['dJ_start'].max():.1f} there, τ = {TAU:g}) over {S_FLOW:g} slow units: the largest final G⁻¹-distance from "
          f"the exact fit {fmt(float(dom['dist_end'].max()))}, the largest final (J − J(exact))/τ {fmt(float(dom['dJ_end'].max()))}.", "",
          "An exit: the first step (dt = 0.25) after which Σ((ϑ − ϑ*)/L)^p ≥ 1; the runs are integrated together by src/core/gauge_batch.py with each "
          "seed's noise as in run_dithered.", "",
          "| τ | ΔJ_𝒪/τ | exits | the cap (slow units) | the total time | the rate (per slow unit, ± s.e.) | the mean exit time of the exits | "
          "the median time | the rate × e^{ΔJ_𝒪/τ} |", "|---|---|---|---|---|---|---|---|---|"]
    for r in Gg['rows']:
        L_.append(f"| {r['tau']:g} | {r['x']:.2f} | {r['n']}/{r['N']} | {r['cap']:g} | {r['T']:.4g} | {fmt(r['rate'])} ± {fmt(r['rate_se'])} | "
                  f"{fmt(r['mean_exit'])} | {r['med']:.3g} | {fmt(r['scaled'])} |")
    L_ += ["", f"Weighted (by the exit counts; var log(rate) ≈ 1/n) least-squares fits of log(rate) = a − c·log τ − B/τ; the mean exit time "
           f"extrapolated to τ = 10⁻⁴ against 1/(γτ) = {1 / (GAMMA_REF * 1e-4):g} at γ = {GAMMA_REF:g}:", "",
           "| prefactor τ^(−c) | B (± s.e.) | B/ΔJ_𝒪 (± s.e.) | c | χ² (dof) | the extrapolated mean exit time at τ = 10⁻⁴ | over 1/(γτ) |",
           "|---|---|---|---|---|---|---|"]
    for f in Gg['fits']:
        cc = f"{f['c']:.2f} ± {f['c_se']:.2f} (free)" if f['c_se'] is not None else f"{f['c']:g}"
        lab = 'free' if f['c_se'] is not None else f"{f['c']:g}"
        L_.append(f"| {lab} | {fmt(f['B'])} ± {fmt(f['B_se'])} | {f['ratio']:.2f} ± {f['ratio_se']:.2f} | "
                  f"{cc} | {f['chi2']:.2f} ({f['dof']}) | {f['t4']:.3g} | {f['t4'] * GAMMA_REF * 1e-4:.2f} |")
    L_.append("")
    return L_


def rows_G(runs, t0=0.0):
    """Part G's exits per τ with the first t0 time units of every run discarded: an exit counts when it comes after t0, and each run
    contributes max(t − t0, 0) to the total time (t its exit time or cap). t0 = 0 gives Part G's n, T and rate."""
    out = []
    for tau in TAUS_G:
        rr = [r for r in runs if r['tau'] == tau]
        t = np.array([r['t'] if r['t'] is not None else r['cap'] for r in rr]); ex = np.array([r['t'] is not None for r in rr])
        n = int(np.sum(ex & (t > t0))); T = float(np.sum(np.maximum(t - t0, 0.0)))
        out.append(dict(tau=tau, n=n, T=T, rate=n / T, first=(float(t[ex].min()) if ex.any() else np.nan)))
    return out


def fixed_barrier(rows, dJ):
    """The barrier fixed at ΔJ_𝒪: the weighted (exit counts) least-squares fit of y = log(rate) + ΔJ_𝒪/τ = a − c·log τ. Returns ĉ, its
    standard error (χ² is quadratic in c, so the profile's Δχ² = 1 and 4 intervals are ĉ ± s.e. and ĉ ± 2 s.e.), χ²(ĉ) on n − 2 degrees of
    freedom with its p-value, the covariance of (a, c) and the fitted log rate at τ = 10⁻⁴ with its standard error; and for the prefactors
    c = 0, ½, 1 (only a fitted) χ² on n − 1 degrees of freedom, the p-value and z = (c − ĉ)/s.e."""
    from scipy.stats import chi2 as C2
    ok = [r for r in rows if r['n'] > 0]
    tau = np.array([r['tau'] for r in ok]); lt = np.log(tau); w = np.array([r['n'] for r in ok], float)
    y = np.log([r['rate'] for r in ok]) + dJ / tau
    X = np.c_[np.ones_like(lt), -lt]
    cov = np.linalg.inv(X.T @ (w[:, None] * X)); beta = cov @ X.T @ (w * y)
    c, se = float(beta[1]), float(np.sqrt(cov[1, 1]))

    def chi2(cc):
        yy = y + cc * lt
        return float(w @ (yy - w @ yy / w.sum()) ** 2)
    x4 = np.r_[1.0, -np.log(1e-4)]
    lr4 = float(x4 @ beta) - dJ / 1e-4
    fixed = [dict(c=cc, chi2=chi2(cc), dof=len(ok) - 1, p=float(C2.sf(chi2(cc), len(ok) - 1)), z=(cc - c) / se) for cc in C_GRID_H]
    return dict(c=c, se=se, chi2=chi2(c), dof=len(ok) - 2, p=float(C2.sf(chi2(c), len(ok) - 2)), int1=(c - se, c + se),
                int2=(c - 2 * se, c + 2 * se), t4=float(np.exp(-lr4)), lse4=float(np.sqrt(x4 @ cov @ x4)), fixed=fixed)


def extrap(rows, dJ, c):
    """arrhenius with the prefactor τ^{−c} fixed, and the standard error of its fitted log rate at τ = 10⁻⁴ (the fit's covariance of the
    intercept and the barrier, var log(rate) ≈ 1/n) and of the barrier's part alone (10⁴ times its standard error)."""
    f = arrhenius(rows, dJ, c)
    ok = [r for r in rows if r['n'] > 0]
    x = np.array([1 / r['tau'] for r in ok]); w = np.array([r['n'] for r in ok], float)
    X = np.c_[np.ones_like(x), -x]
    cov = np.linalg.inv(X.T @ (w[:, None] * X)); x4 = np.r_[1.0, -1e4]
    return dict(c=c, ratio=f['ratio'], ratio_se=f['ratio_se'], t4=f['t4'], lse4=float(np.sqrt(x4 @ cov @ x4)), lseB=float(1e4 * f['B_se']))


def cluster_minima(vals, U):
    """The distinct minima among L-BFGS results with the values vals ((J − J(exact))/τ) and the normalised deviations U ((ϑ − ϑ*)/L), in
    increasing value: one minimum when the values agree within CLUSTER_H[0] and the deviations within CLUSTER_H[1] in every coordinate.
    Returns [dict(dJ, u, n, i)] (i: the first result of the cluster, n: its size)."""
    cl = []
    for i in np.argsort(vals, kind='stable'):
        for c in cl:
            if abs(vals[i] - c['dJ']) < CLUSTER_H[0] and np.abs(U[i] - c['u']).max() < CLUSTER_H[1]:
                c['n'] += 1
                break
        else:
            cl.append(dict(dJ=float(vals[i]), u=np.asarray(U[i]), n=1, i=int(i)))
    return cl


def boundary_minima():
    """The local minima of J on ∂𝒪: domain_job's N_DOM L-BFGS runs from the same starts (the generator 1700), all kept and clustered
    (CLUSTER_H), by (J − J(exact))/τ; and the deterministic flow (Euler, dt = DT) from each over S_FLOW: the largest Σ((ϑ − ϑ*)/L)^p after
    the first step (below 1: the flow enters 𝒪 at once and stays in it), the final G⁻¹-distance from the exact fit and (J − J(exact))/τ."""
    from scipy.optimize import minimize
    from src.core.gauge_batch import BatchDither
    net, base = E.setup()
    beta = np.ones(3)
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    Th0 = E.exact_rep(net, np.ones(3) / 3); J0 = loop.J(Th0, beta)
    Lv = np.tile(np.asarray(DOM_L, float), 3)
    pt = lambda u: (Th0.ravel() + Lv * u / np.sum(np.abs(u) ** DOM_P) ** (1 / DOM_P)).reshape(3, 3)
    f = lambda u: (loop.J(pt(u), beta) - J0) / TAU
    rng = np.random.default_rng(1700)
    res = []
    for _ in range(N_DOM):
        r = minimize(f, rng.standard_normal(9), method='L-BFGS-B', options=dict(maxiter=5000, ftol=1e-15, gtol=1e-12))
        res.append((float(r.fun), pt(r.x)))
    cl = cluster_minima([v for v, _ in res], [(T - Th0) / np.asarray(DOM_L)[None] for _, T in res])
    for c in cl:
        c['Th'] = res[c.pop('i')][1]
    bd = BatchDither(loop)
    T = np.array([c['Th'] for c in cl]); lb = np.zeros((len(cl), 3)); phi = np.zeros(len(cl))
    for _ in range(int(round(S_FLOW / DT))):
        T, lb = bd.step(T, lb, DT, np.zeros(3), 0.0)
        phi = np.maximum(phi, dom_phi(T, Th0))
    for c, x, ph in zip(cl, T, phi):
        c.update(phi_max=float(ph), dist_end=float(np.sqrt(np.sum((x - Th0) ** 2 / E.G[None, :]))), dJ_end=float((loop.J(x, beta) - J0) / TAU))
        del c['Th']
    return dict(minima=cl, n=N_DOM)


def exit_points(runs, ckpt):
    """Part G's exit points from the runs' checkpoint (run_batch keeps each run's state at its exit step: the first state outside 𝒪), each
    exit step checked against the stored exit time. Returns [dict(tau, seed, dev = Θ − Θ*)] for the runs that exit, or None without the
    checkpoint."""
    from src.core.gauge_batch import load_states
    path = None if ckpt is None else os.path.join(ckpt, 'plateau_G.pkl')
    if path is None or not os.path.exists(path):
        return None
    st = load_states(path)
    net, _ = E.setup()
    Th0 = E.exact_rep(net, np.ones(3) / 3)
    out = []
    for r in runs:
        s = st[(r['tau'], r['seed'])]
        assert (None if s['exit'] is None else s['exit'] * DT) == r['t'], (r['tau'], r['seed'])
        if r['t'] is not None:
            out.append(dict(tau=r['tau'], seed=r['seed'], dev=np.array(s['Th'], float) - Th0))
    return out


def exits_by_tau(ex):
    """By τ: the exits, the coordinate (class, j) whose term dominates Σ((ϑ − ϑ*)/L)^p at the exit point (counts), and the exits through
    the second width broadened and narrowed (its term dominant, its deviation positive or negative)."""
    out = []
    for tau in TAUS_G:
        U = np.array([e['dev'] / np.asarray(DOM_L)[None] for e in ex if e['tau'] == tau])
        top = [np.unravel_index(int(np.argmax(u ** DOM_P)), (3, 3)) for u in U]
        cnt = {}
        for k, j in top:
            cnt[(int(k), int(j))] = cnt.get((int(k), int(j)), 0) + 1
        w2 = np.array([(k, j) == (1, 2) for k, j in top])
        out.append(dict(tau=tau, n=len(U), top=sorted(cnt.items()), broad=int(np.sum(w2 & (U[:, 1, 2] > 0))),
                        narrow=int(np.sum(w2 & (U[:, 1, 2] < 0)))))
    return out


def part_H(dom, runs, C, mins, ex):
    dJ = dom['dJ']
    lags = []
    for t0 in (0.0,) + LAGS_H:
        rr = rows_G(runs, t0)
        lags.append(dict(t0=t0, n=[r['n'] for r in rr], fits=[arrhenius(rr, dJ, c) for c in C_GRID_H], fb=fixed_barrier(rr, dJ)))
    rr = rows_G(runs)
    first = [dict(tau=r['tau'], first=r['first'], before=int(sum(1 for x in runs if x['tau'] == r['tau'] and x['t'] is not None
                                                                  and x['t'] < LAGS_H[1]))) for r in rr]
    ext = [extrap(rr, dJ, c) for c in C_GRID_H]
    extC = extrap(C['rows'], C['dJ'], 0.5) if C is not None else None
    return dict(fb=fixed_barrier(rr, dJ), lags=lags, first=first, ext=ext, extC=extC, mins=mins,
                exits=(None if ex is None else exits_by_tau(ex)))


def report_H(Hh):
    fmt = E.fmt
    fb = Hh['fb']
    L_ = ["## H. The prefactor with the barrier fixed, the lag, the extrapolations' uncertainty, the minima on ∂𝒪 and the exit points "
          "(referee round 150; from Part G's stored runs)", "",
          "With the barrier fixed at ΔJ_𝒪: the weighted (exit counts) least-squares fit of log(rate) + ΔJ_𝒪/τ = a − c·log τ (χ² is quadratic "
          "in c, so the profile's Δχ² = 1 and 4 intervals are ĉ ± s.e. and ĉ ± 2 s.e.):", "",
          "| prefactor τ^(−c) | c (± s.e.) | Δχ² ≤ 1 | Δχ² ≤ 4 | χ² (dof) | p | z = (c − ĉ)/s.e. |", "|---|---|---|---|---|---|---|",
          f"| fitted | {fb['c']:.2f} ± {fb['se']:.2f} | {fb['int1'][0]:.2f} to {fb['int1'][1]:.2f} | {fb['int2'][0]:.2f} to {fb['int2'][1]:.2f} | "
          f"{fb['chi2']:.2f} ({fb['dof']}) | {fb['p']:.2f} | — |"]
    for f in fb['fixed']:
        L_.append(f"| {f['c']:g} | {f['c']:g} | — | — | {f['chi2']:.2f} ({f['dof']}) | {f['p']:.3f} | {f['z']:+.2f} |")
    L_ += ["", f"The lag: the first t₀ time units of every run discarded (an exit counts after t₀; each run contributes max(t − t₀, 0)); "
           f"t₀ = {', '.join(f'{v:g}' for v in LAGS_H)} is ½, 1, 2 times the slowest relaxation time 1/(3.06·10⁻⁴):", "",
           "| t₀ | exits by τ | B/ΔJ_𝒪 (± s.e.), c = 0 | c = ½ | c = 1 | ĉ with the barrier ΔJ_𝒪 (± s.e.) |", "|---|---|---|---|---|---|"]
    for g in Hh['lags']:
        L_.append(f"| {g['t0']:g} | {', '.join(str(v) for v in g['n'])} | "
                  + ' | '.join(f"{f['ratio']:.3f} ± {f['ratio_se']:.3f}" for f in g['fits']) + f" | {g['fb']['c']:.2f} ± {g['fb']['se']:.2f} |")
    L_ += ["", f"The first exit by τ (and the number of exits before {LAGS_H[1]:g}): "
           + '; '.join(f"{r['tau']:g}: {r['first']:g} ({r['before']})" for r in Hh['first']) + ".", "",
           f"The mean exit time extrapolated to τ = 10⁻⁴ (1/(γτ) = {1 / (GAMMA_REF * 1e-4):g} at γ = {GAMMA_REF:g}), with the standard error of the "
           "fitted log rate there (the fit's covariance, var log(rate) ≈ 1/n), the range e^{∓s.e.} and the barrier's part alone (10⁴ × its s.e.):", "",
           "| fit | the mean exit time at τ = 10⁻⁴ | s.e. of its log | the range (one s.e.) | over 1/(γτ) | the barrier's part of the s.e. |",
           "|---|---|---|---|---|---|"]
    t = 1 / (GAMMA_REF * 1e-4)
    rws = [(f"G, c = {e['c']:g}", e['t4'], e['lse4'], e['lseB']) for e in Hh['ext']]
    rws.append((f"G, the barrier ΔJ_𝒪, ĉ = {fb['c']:.2f}", fb['t4'], fb['lse4'], None))
    if Hh['extC'] is not None:
        rws.append(("C (the escape set max|α_k| < 6), c = ½", Hh['extC']['t4'], Hh['extC']['lse4'], Hh['extC']['lseB']))
    for lab, t4, s, sB in rws:
        L_.append(f"| {lab} | {t4:.3g} | {s:.2f} | {t4 * np.exp(-s):.3g}–{t4 * np.exp(s):.3g} | {t4 / t:.2f} ({t4 * np.exp(-s) / t:.2f}–{t4 * np.exp(s) / t:.2f}) | "
                  f"{'—' if sB is None else f'{sB:.2f}'} |")
    m = Hh['mins']
    L_ += ["", f"The local minima of J on ∂𝒪 (Part G's {m['n']} L-BFGS runs from the same starts, clustered: (J − J(exact))/τ within "
           f"{CLUSTER_H[0]:g} and (ϑ − ϑ*)/L within {CLUSTER_H[1]:g}), and the deterministic flow from each over {S_FLOW:g} slow units "
           "(1 minus the largest Σ((ϑ − ϑ*)/L)^p after the first step: positive, the flow enters 𝒪 at once and stays in it):", "",
           "| (J − J(exact))/τ | starts | (ϑ − ϑ*)/L by class (α, γ, κ) | 1 − the largest Σ((ϑ − ϑ*)/L)^p on the flow | the final G⁻¹-distance | "
           "the final (J − J(exact))/τ |", "|---|---|---|---|---|---|"]
    for c in m['minima']:
        L_.append(f"| {c['dJ']:.3f} | {c['n']} | {'; '.join(', '.join(f'{v:+.2f}' for v in r) for r in c['u'])} | {fmt(1 - c['phi_max'])} | "
                  f"{fmt(c['dist_end'])} | {fmt(c['dJ_end'])} |")
    L_.append("")
    if Hh['exits'] is None:
        L_ += ["The exit points were not recorded (no checkpoint of Part G's runs).", ""]
    else:
        L_ += ["The exit points (the first state outside 𝒪, read from Part G's checkpoint and checked against the stored exit times): the term "
               "of Σ((ϑ − ϑ*)/L)^p that dominates there, by (class, coordinate), and the exits through the second width (its term dominant) "
               "broadened, on the side of the lowest minimum, and narrowed:", "",
               "| τ | exits | the dominant term: (class, coordinate) counts | the second width broadened | narrowed | broadened, fraction |",
               "|---|---|---|---|---|---|"]
        for r in Hh['exits']:
            L_.append(f"| {r['tau']:g} | {r['n']} | {'; '.join(f'({k + 1}, {chr(945) if j == 0 else (chr(947) if j == 1 else chr(954))}): {v}' for (k, j), v in r['top'])} | "
                      f"{r['broad']} | {r['narrow']} | {r['broad'] / r['n']:.2f} |")
        L_.append("")
    return L_


def report(A=None, B=None, C=None, D=None, Ee=None, Ff=None, Gg=None, Hh=None):
    fmt = E.fmt
    v3 = lambda v, p=3: ', '.join(f'{x:.{p}f}' for x in v)
    L_ = ["# Proposition D‴.39: dissipativity of the dithered gauge (`exp_gauge_plateau.py`)", "",
          f"Proposition D‴.36's instance (three bumps, one particle per class, m = 3, the cohesion leak); energies over τ = {TAU:g}; the box "
          f"L = ({v3(BOX, 1)}) in (α, γ, κ), σ_B = {SIGMA_B:g}: |a| ≤ {A_MAX * np.tanh(BOX[0] / A_MAX):.3f}, the representable masses "
          f"W_k ≥ {v3(W_BOX, 4)} (the saturation's {v3(W_MIN, 4)}).", ""]
    if A is not None:
        L_ += ["## A. The plateau: the lowest configuration with one amplitude fixed far out", "",
               "| masses | class | α_k | ΔJ/τ above the exact fit | fit error | a | c | s |", "|---|---|---|---|---|---|---|---|"]
        for r in A:
            L_.append(f"| {r['which']} | {r['k'] + 1} | {dict(sat=f'±{AL_SAT:g}', esc=f'±{ESC:g}', box='±L_α (box)')[r['at']]} | {r['dJ']:.3f} | {r['fit']:.4f} | {v3(r['a'], 2)} | "
                      f"{v3(r['c'], 2)} | {v3(r['s'], 2)} |")
        L_.append("")
    if B is not None:
        L_ += ["## B. The minimum-energy path from the exact fit to the plateau (string method)", ""]
        for r in B:
            p = r['prof']
            s2 = [E.setup()[0].physical(T)[2][1] for T in r['images'][::4]]
            L_ += [f"Masses {r['which']}, toward class {r['k'] + 1}'s saturation: the largest ΔJ/τ along the path {p.max():.3f} (image {int(p.argmax())} of "
                   f"{len(p) - 1}); at every fourth image the profile {', '.join(f'{v:.2f}' for v in p[::4])}, the fit error "
                   f"{', '.join(f'{v:.4f}' for v in r['fits'][::4])}, the second width {', '.join(f'{v:.2f}' for v in s2)}.", ""]
    if C is not None:
        L_ += [f"## C. The escapes at fixed masses 1/3 (γ = 0; {len(SEEDS_C)} seeds, at most {S_C:g} slow units; an escape: max|α_k| > {ESC:g})", "",
               "| τ | escapes | the rate (per slow unit, ± s.e.) | the median time | the rate × e^{ΔJ/τ} | the fit error at the escape |", "|---|---|---|---|---|---|"]
        for r in C['rows']:
            L_.append(f"| {r['tau']:g} | {r['n']}/{r['N']} | {fmt(r['rate'])} ± {fmt(r['rate_se'])} | {r['med']:.3g} | {fmt(r['scaled'])} | {fmt(r['fit'])} |")
        ok = [r for r in C['rows'] if r['n'] > 0]
        x = np.array([1 / r['tau'] for r in ok]); yv = np.log([r['rate'] for r in ok]); w = np.array([r['n'] for r in ok], float)
        xb, yb = w @ x / w.sum(), w @ yv / w.sum()
        sl, se_sl = (w @ ((x - xb) * (yv - yb))) / (w @ (x - xb) ** 2), 1 / np.sqrt(w @ (x - xb) ** 2)
        L_ += ["", f"The slope of log(rate) against 1/τ: {fmt(-C['slope'])} (the lowest J on the escape set's boundary, ΔJ = {fmt(C['dJ'])}); "
               f"weighted by the escape counts (log(rate)'s variance ≈ 1/n): {fmt(-sl)} ± {fmt(se_sl)}, {-sl / C['dJ']:.2f} ± {se_sl / C['dJ']:.2f} of ΔJ "
               "(referee round 31).", ""]
    if D is not None:
        L_ += [f"## D. The dead zone at fixed masses 1/3 (γ = 0, τ = {TAU_D:g}, {len(SEEDS_D)} seeds, {S_D:g} slow units)", "",
               f"A record is on the plateau when its second width exceeds {S2_PLATEAU:g} (A and B: the plateaus' configurations and the path's).", "",
               f"| box | seeds that reach the plateau | fraction of the records on it (s_2 > {S2_PLATEAU:g}), by quarters | fraction with max|α| > {ESC:g} | "
               "the largest |α|, by quarters | the fit error, by quarters |", "|---|---|---|---|---|---|"]
        for r in D:
            L_.append(f"| {'L = (' + v3(BOX, 1) + ')' if r['boxed'] else 'none'} | {r['seeds']}/{r['N']} | {', '.join(f'{v:.3f}' for v in r['plat'])} | "
                      f"{', '.join(f'{v:.3f}' for v in r['frac'])} | "
                      f"{', '.join(f'{v:.1f}' for v in r['almax'])} | {', '.join(f'{v:.4f}' for v in r['fit'])} |")
        L_.append("")
    if Ee is not None:
        ws = Q.w_sharp(); V = vertices()
        L_ += [f"## E. The full loop with the box (γ = {GAMMA_E:g}, τ = {TAU:g}, {len(SEEDS_E)} seeds, {S_E:g} slow units)", "",
               f"The box's vertices: {'; '.join(v3(v, 4) for v in V)}; W♯ = {v3(ws, 4)}.", "",
               "| κ_F | s | mean of W_1, W_2, W_3 (± s.e.) | std |", "|---|---|---|---|"]
        for r in Ee:
            for t in r['sum']['times']:
                L_.append(f"| {r['kF']:g} | {t['s']:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(t['mean'], t['se']))} | {v3(t['std'], 4)} |")
        L_ += ["", "| κ_F | the masses at the end, by seed (the nearest vertex for κ_F = 0, its distance) | the largest |α| | the fit error (second half) |",
               "|---|---|---|---|"]
        for r in Ee:
            sm = r['sum']
            if r['kF'] == 0:
                d = np.linalg.norm(sm['end'][:, None, :] - V[None, :, :], axis=2)
                cell = '; '.join(f"{v3(w, 3)} (v{int(np.argmin(dd)) + 1}, {dd.min():.3f})" for w, dd in zip(sm['end'], d))
            else:
                cell = '; '.join(f"{v3(w, 3)} ({np.abs(w - ws).max():.3f})" for w in sm['end'])
            L_.append(f"| {r['kF']:g} | {cell} | {sm['almax']:.2f} | {fmt(sm['fit'])} |")
        for r in Ee:
            sp = r['sum']['split']
            if sp is not None:
                src = "exp_gauge_drift.py's Part H" if r['kF'] == 0 else "exp_gauge_fixman.py's Part E"
                L_ += ["", f"κ_F = {r['kF']:g}: the records coincide with {src} (the same seeds without the box) up to the first record at which a "
                       f"coordinate has left B: {', '.join('never' if v is None else f'{v:g}' for v in sp)}."]
        L_.append("")
    if Ff is not None:
        L_ += ["## F. The two channels to the saturation in the stored runs", "",
               "| runs | class | seeds with records at |α_k| > 10 | such records | of them with W_k − |A_k|/a_max < 0.01 | > 0.05 | their mean fit error |",
               "|---|---|---|---|---|---|---|"]
        for r in Ff:
            for k, q in enumerate(r['rows']):
                L_.append(f"| {r['src']} | {k + 1} | {q['seeds']}/{r['N']} | {q['n']} | {q['near']} | {q['far']} | {fmt(q['fit'])} |")
        L_.append("")
    if Gg is not None:
        L_ += report_G(Gg)
    if Hh is not None:
        L_ += report_H(Hh)
    L_.append("Numbers are produced by `experiments/exp_gauge_plateau.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_plateau_raw.npz')
    md = os.path.join(OUT, 'gauge_plateau.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        txt = report(*[get(k) for k in 'ABCDEF'], Gg=get('G'), Hh=get('H'))       # before the md is opened
        open(md, 'w').write(txt)
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEF')
    if 'F' in parts:
        old['F'] = part_F(); E.save(raw, **old)
    if 'G' in parts:                                                         # referee round 143: resumable, no worker pool
        if 'G_dom' not in old:
            old['G_dom'] = domain_job(); E.save(raw, **old)
            print(f"G: ΔJ_O/τ = {old['G_dom']['dJ'] / TAU:.3f}", flush=True)
        runs = escape_runs(old['G_dom'])
        if runs is None:
            print("G: open runs remain (checkpointed); rerun to resume", flush=True)
            return
        old['G_runs'] = runs; old['G'] = part_G(old['G_dom'], runs); E.save(raw, **old)
        parts = parts - {'G'}
        if not parts:
            txt = report(*[old.get(k) for k in 'ABCDEF'], Gg=old.get('G'), Hh=old.get('H'))
            open(md, 'w').write(txt)
            return
    if 'H' in parts:                                                         # referee round 150: from the stored runs, no worker pool
        if 'H_min' not in old:
            old['H_min'] = boundary_minima(); E.save(raw, **old)
        if old.get('G_exits') is None:
            ex = exit_points(old['G_runs'], CKPT)
            if ex is not None:
                old['G_exits'] = ex
        old['H'] = part_H(old['G_dom'], old['G_runs'], old.get('C'), old['H_min'], old.get('G_exits')); E.save(raw, **old)
        parts = parts - {'H'}
        if not parts:
            txt = report(*[old.get(k) for k in 'ABCDEF'], Gg=old.get('G'), Hh=old.get('H'))
            open(md, 'w').write(txt)
            return
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'A' in parts:
            res['A'] = pool.map(plateau_job, [(w, k, at) for at in ('sat', 'esc', 'box') for w in ('uniform', 'sharp') for k in range(3)], chunksize=1)
            E.save(raw, **res); log('A')
        if 'E' in parts:
            h['E'] = pool.map_async(series, [(GAMMA_E, TAU, k, f, True, s, S_E) for k, f in ((0.0, 'none'), (2.0, 'masses')) for s in SEEDS_E], chunksize=1)
        if 'B' in parts:
            ends = []
            for w in ('uniform', 'sharp'):
                rr = [r for r in res['A'] if r['which'] == w and r['at'] == 'sat']
                best = min(rr, key=lambda r: r['dJ']); ends.append((w, best['k'], best['Th']))
            h['B'] = pool.map_async(string_job, ends, chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(series, [(0.0, TAU_D, 0.0, 'none', b, s, S_D) for b in (False, True) for s in SEEDS_D], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(escape_job, [(t, s) for t in TAUS_C for s in SEEDS_C], chunksize=1)
        if 'B' in h:
            res['B'] = h['B'].get(); E.save(raw, **res); log('B')
        if 'C' in h:
            dJ = min(r['dJ'] for r in res['A'] if r['which'] == 'uniform' and r['at'] == 'esc')
            res['C_runs'] = h['C'].get(); res['C'] = part_C(res['C_runs'], dJ); E.save(raw, **res); log('C')
        if 'D' in h:
            res['D_rec'] = h['D'].get(); res['D'] = part_D(res['D_rec']); E.save(raw, **res); log('D')
        if 'E' in h:
            runs = h['E'].get()
            res['E_rec'] = {k: np.array([r['rec'] for r in runs if r['kF'] == k]) for k in (0.0, 2.0)}
            refs = {0.0: E.load(os.path.join(OUT, 'gauge_drift_raw.npz')).get('H_rec'),
                    2.0: (E.load(os.path.join(OUT, 'gauge_fixman_raw.npz')).get('E_rec') or {}).get(2.0)}
            res['E'] = [dict(kF=k, sum=summarise_E(res['E_rec'][k], refs[k])) for k in (0.0, 2.0)]
            E.save(raw, **res); log('E')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEF'], Gg=res.get('G'), Hh=res.get('H')))


if __name__ == '__main__':
    main()
