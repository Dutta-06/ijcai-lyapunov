"""
Proposition D‴.40: the gauge under general noise — the replicator's averaged fitness is the generator's image of one coordinate
==========================================================================================================================

Proposition D‴.36's instance (Proposition D‴.32's three bumps, one particle per class, m = 3). With white noise of constant covariance Σ (per
unit slow time) on the positions' unconstrained coordinates (referee round 147: the noise class and the chart), the gauge identity u_k = −g_k∂_{α_k}J_fit (g_k = a_k/(β_k a′_k)) makes the masses' fitness the generator's image of
ψ_k = ∫g_k dα_k/G_α: in a stationary law of the positions at fixed masses (with Proposition D‴.39's box) E[u_k] = −τ_kE[g′_k(α_k)] +
E[g_k∂_{α_k}(J_coh + V)], τ_k := Σ_{α_kα_k}/(2G_α); at leading order ū_k = −(τ_k/m)(1 + (a_k/a_max)²)/(W_k a′_k), only each amplitude's own
level entering — the replicator flow of Σ_kτ_k log(W_k a′_k).

A. The trace formula (`gauge_noise.averaged_fitness`: ū = ∇u·μ + ½tr(∇²u C), C the Lyapunov solution, μ the mean shift) against the
   closed form, for isotropic noise, noise on the amplitudes only, on the centres and widths only, per-class amplitude levels and three
   random correlated covariances, at three points of M.
B. Katzenberger's limit drift at γ = 10⁻⁵ (`gauge_drift.limit_drift` with the noise Σ) against the replicator flow of the closed form,
   for the same noises at two points; and its diffusion.
C. Noise on the centres and widths only at γ = 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵: the masses' diffusion (zero: the gauge's tangent lies in the
   amplitude–mass coordinates and the metric has no block coupling them to the centres and widths) and drift (O(γ)), against the isotropic noise's.
D. The identity in simulation at fixed masses 1/3 (γ = 0), inside the box of Proposition D‴.39 (so the stationary law exists): the time
   averages of u_k and of each term of the generator identity E[u_k] = −τ_kE[g′_k] + E[g_k∂_{α_k}J_coh] + E[g_k∂_{α_k}V] + E[Lψ_k] (ψ_k =
   ∫g_k dα_k/G_α; the last term the Euler chain's bias, measured exactly from the Gaussian increment's closed-form expectation), with
   (i)'s box-free residual E[u_k + τ_kg′_k] and its z, the share of steps outside the box, and the whole identity's residual against its
   martingale error; isotropic noise, the amplitudes only (the same τ_k, no centre or width noise), 'correlated 1' at τ = 10⁻³ and 10⁻²,
   and the centres and widths only at 10⁻²; the amplitudes' noise intensity from their quadratic variation, with and without the drift's
   part (referee round 31).
E. The full loop inside the box at γ = 10⁻², τ = 10⁻⁴ (the 12 seeds of exp_gauge_drift.py's Part H, 10⁶ slow units): the centres and
   widths dithered only (the masses frozen); the same with the barrier of Proposition D‴.38's masses form at κ_F = 1 (no amplitude noise
   to cancel: the masses ascend ℓ to W♯); a correlated noise on all coordinates whose amplitudes have the levels (2, 1, 0.5)τ and whose
   centres and widths stay at the dither's level ('correlated E'), uncorrected, and with the masses form at the classes' levels τ_k and
   κ_F = 2 (the masses ascend Σ_kτ_k log(W_k a′_k) to its maximiser W♯_τ, 0.084 from W♯).
F. (referee round 147) Part E against its predictions: the slow limit's leading-order ODE (`f_ode`: the replicator flow of the closed form
   with the masses form, P_B's walls enforced as the box does, `wall_project`, `wall_push`) and Katzenberger's limit at γ = 10⁻² (`katz`,
   `f_field` on `grid_PB`; `f_paths`: its drift's path and its SDE's ensemble, reflected at the walls); the zeros of Katzenberger's corrected
   drift at γ = 10⁻², 10⁻³, 10⁻⁵ (`f_zero`); W♯_τ against P_B for the levels (2, 1, 0.5), (2, 1, 0.2), (1, 1, 0.1).
G. (referee round 147) The early phase: 'correlated E' with and without the class levels, 36 new seeds over 2·10⁴ slow units (`g_series`);
   two of Part E's seeds replayed against its records; the paired difference against Katzenberger's.
H. (referee round 147) Coloured noise: the linear response's levels τ_k(τ_c) = [(I + τ_cGH)⁻¹Σ]_{α_kα_k}/(2G_α) at the masses 1/3
   (`part_H_linear`; `gauge_noise.coloured_cross` against the augmented system's Lyapunov equation), and the identity in simulation inside the
   box (`h_mc`): the amplitudes only at τ = 10⁻³, an Ornstein–Uhlenbeck noise at τ_c = 100 against white, one seed.
I. (round 147's open item) The early phase's gap: 'correlated E' (uncorrected) over 2·10⁴ slow units for 192 new seeds, the same without
   the box and from the positions' equilibrated state (`i_runs`, batched; each run with the gauge identity's decomposition: the boundary
   term, the martingale (mean zero), the noise's drift, the box's, the cohesion's, the chain's bias); Part G's and Part E's seeds replayed and
   decomposed; the box's term at fixed masses against W_1 (`i_fixed`); Katzenberger's SDE at smaller steps, without walls (`i_field` on
   `grid_ext`) and with the measured box term through Katzenberger's projection (`i_sde`); the multiple-time tests (Holm).
J. (referee round 153) The identity's own check: Part D's runs as two Euler–Maruyama chains at the steps Δt and Δt/2 driven by one Brownian
   path (`j_runs`), the chain's exact per-state bias at both steps (O(Δt) if and only if the identity's right side is the generator's
   image), against the same with the right side lacking its Itô term. (Part D's χ² holds by construction and checks the averaging.)
K. (referee round 153) The Euler chain's bias along Part E's two deterministic-amplitude variants: the runs replayed with the gauge
   identity's decomposition (`k_runs`), the positions' offset from the fit (`k_offset`) and the chain's term through the masses' linear
   response (`k_response`) against the second half's offsets from Katzenberger's SDE.
L. (referee round 153) From the raw: Part H's coloured-over-white ratios with delta-method errors, Part F's zero shift per class, Part C's
   and H's extremes, Part I's regressed s.e. with the regression coefficient's estimation.
M. (D‴.40 after round 153: the missing O(γ²τ) drift of the deterministic-amplitude loop) Katzenberger's drift at the masses 1/3 from the
   linearised positions' antisymmetric correlation (`m_ou`) and from the limit map's second differences (`m_limit`, `m_fd`) against
   `limit_drift`; Part K's records with the exact limit map (`m_lim_K`) against its first-order offset and with the fit's motion
   (`m_offset_lin`); new coupled runs at Δt and Δt/2 at γ = 10⁻², 3·10⁻³ (`m_loop`), the limit map's rate with Richardson's extrapolation;
   the antisymmetric correlation at fixed masses at three levels τ (`m_fixed`, `m_anti`). Resumable: each piece is saved to the raw as it
   completes (`part_M`).
N. (D‴.40 after Part M: the finite-τ question) More seeds of M(d)'s fixed-mass runs (`m_fixed`, unchanged) at τ = 10⁻⁴ and 3·10⁻⁴ (Part
   M(d)'s run lengths, so that the seeds pool) and 16 at 10⁻³ (`part_N`, resumable as M); the pooled ratio to the closed form per τ
   (`n_pool`, `n_ratios`) and its weighted line against τ (`n_fit`): the slope and the limit τ → 0.

Output: results/chaos/gauge_noise.md and gauge_noise_raw.npz (GN_OUT, GN_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`; F needs E's records,
I needs F's and G's, K needs E's and F's, M needs C's, E's, F's and K's, N's report needs M's).
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
from src.core import gauge_drift as GDR
from src.core import gauge_noise as GN
from src.core.runcache import cached
import exp_cohesion_leak as E
import exp_gauge_dither as GD
import exp_gauge_drift as X
import exp_gauge_fixman as Q
import exp_gauge_plateau as P

OUT = os.environ.get('GN_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GN_WORKERS', 4))
TAU, A_MAX, DT, REC = GD.TAU, 4.0, GD.DT, GD.REC
G9 = np.tile(E.G, 3)
W_A = ((1 / 3, 1 / 3), (0.30, 0.35), (0.40, 0.25))
W_B = ((1 / 3, 1 / 3), (0.30, 0.35))
GAMMAS_C = (1e-2, 1e-3, 1e-4, 1e-5)
TAUS_D, S_D, SEED_D, N_BATCH, N_BATCH_LONG = (1e-3, 1e-2), 2e5, 700, 40, 10
RUNS_D = tuple((n, t) for n in ('isotropic', 'amplitudes only', 'correlated 1') for t in TAUS_D) + (('centres and widths only', 1e-2),)
GAMMA_E, SEEDS_E, S_E, T_E = 1e-2, X.SEEDS_H, 1e6, X.T_H
NOISE_SEEDS = (11, 12, 13)


def noises(tau=TAU):
    """name → Σ (9 × 9, per unit slow time) on the flattened positions (α, γ, κ per particle)."""
    out = {'isotropic': 2 * tau * np.diag(G9),
           'amplitudes only': 2 * tau * np.diag(np.tile([1.0, 0.0, 0.0], 3)),
           'centres and widths only': 2 * tau * np.diag(np.tile([0.0, E.G[1], E.G[2]], 3)),
           'amplitude levels (2, 1, 0.5)τ': 2 * tau * np.diag(np.repeat([2.0, 1.0, 0.5], 3) * np.tile([1.0, 0.0, 0.0], 3))}
    for i, sd in enumerate(NOISE_SEEDS):
        Lr = np.random.default_rng(sd).standard_normal((9, 9))
        out[f'correlated {i + 1}'] = 2 * tau * 0.3 * (Lr @ Lr.T)
    out['correlated E'] = 2 * tau * corr_E()
    return out


LEVELS_E = (2.0, 1.0, 0.5)


def corr_E():
    """Part E's correlated noise over 2τ: the amplitudes at the levels LEVELS_E·G_α, the centres and widths at the dither's level G_γ, G_κ,
    correlated through the mixture ½I + ½R of a random correlation matrix R (seed 21) — so the amplitudes' levels differ and every
    coordinate's noise is correlated, while the centres and widths are not heated beyond the dither (referee round 31)."""
    Lr = np.random.default_rng(21).standard_normal((9, 9))
    C = Lr @ Lr.T; d = np.sqrt(np.diag(C))
    R = 0.5 * np.eye(9) + 0.5 * C / np.outer(d, d)
    sd = np.sqrt(np.tile(E.G, 3) * np.ravel([[lv, 1.0, 1.0] for lv in LEVELS_E]))
    return R * np.outer(sd, sd)


def levels(Sig):
    """τ_k = Σ_{α_kα_k}/(2G_α)."""
    return np.diag(Sig)[0::3] / (2 * E.G[0])


def closed(W, tau_k):
    """−(τ_k/m)(1 + (a_k/a_max)²)/(W_k a′_k) at the exact fit with the masses W."""
    W = np.asarray(W, float); q = (np.abs(E.A_BUMP) / W / A_MAX) ** 2
    return -(np.asarray(tau_k, float) / 3) * (1 + q) / (W * (1 - q))


def w_sharp_weighted(tk):
    """argmax of Σ_kτ_k log(W_k − W_min,k²/W_k) on the simplex (Lagrange: τ_k(1 + W_min,k²/W_k²)/(W_k − W_min,k²/W_k) equal) — the
    maximiser of the barrier weighted by the classes' noise levels (not an analytic centre: the barrier's terms are not logarithms of affine
    functions; referee round 153); Q.w_sharp() for equal levels. Every level must be positive (with a
    zero level F_τ is flat in that class's mass and its maximiser lies on the saturation: referee round 147). The multiplier's bracket ends at
    the smallest of the g_k at their lower ends (it was 10¹², above that for a small level: a sign error at the levels (2, 1, 0.2))."""
    from scipy.optimize import brentq
    tk = np.asarray(tk, float)
    if np.any(tk <= 0):
        raise ValueError('every level τ_k must be positive')
    Wm = np.abs(E.A_BUMP) / A_MAX; tk = tk / np.mean(tk)
    g = lambda W, k: tk[k] * (1 + Wm[k] ** 2 / W ** 2) / (W - Wm[k] ** 2 / W)        # decreasing in W
    lo = Wm * (1 + 1e-12)
    W_of = lambda mu, k: brentq(lambda W: g(W, k) - mu, lo[k], 1.0)
    mu = brentq(lambda mu: sum(W_of(mu, k) for k in range(3)) - 1, max(g(1.0, k) for k in range(3)), min(g(lo[k], k) for k in range(3)))
    return np.array([W_of(mu, k) for k in range(3)])


@cached
def a_job(args):
    name, w = args
    W = np.r_[w, 1 - sum(w)]
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    beta = 3 * W
    gradJ = lambda t: -loop.drift(t.reshape(3, 3), beta)[0].ravel() / G9
    u = lambda t: loop.fitness(loop.residual(t.reshape(3, 3), beta), t.reshape(3, 3))
    Sig = noises()[name]
    ub = GN.averaged_fitness(gradJ, u, E.exact_rep(net, W).ravel(), Sig, G9)[0]
    cf = closed(W, levels(Sig))
    flow = lambda x: GN.replicator_flow(W, x)
    return dict(name=name, w=w, ubar=ub, closed=cf, rel=np.abs(flow(ub) - flow(cf)).max() / max(np.abs(flow(cf)).max(), 1e-30),
                flow=flow(ub) / TAU)


@cached
def b_job(args):
    name, w = args
    Sig = noises()[name]
    F, chart, coords, S = X.instance(1e-5)
    S11 = np.zeros((11, 11)); S11[:9, :9] = Sig
    b, D, _ = GDR.limit_drift(F, chart, coords, np.asarray(w, float), S11, h=X.H_ID, h2=X.H2_ID)
    W = np.r_[w, 1 - sum(w)]
    pred = GN.replicator_flow(W, closed(W, levels(Sig))) * 1e-5
    return dict(name=name, w=w, b=np.r_[b, -b.sum()], pred=pred, D=D)


@cached
def c_job(args):
    gamma, w = args
    F, chart, coords, S = X.instance(gamma)
    out = {}
    for name in ('isotropic', 'centres and widths only'):
        S11 = np.zeros((11, 11)); S11[:9, :9] = noises()[name]
        b, D, _ = GDR.limit_drift(F, chart, coords, np.asarray(w, float), S11, h=X.H_ID, h2=X.H2_ID)
        out[name] = dict(b=np.r_[b, -b.sum()] / (gamma * TAU), D=D / (gamma ** 2 * TAU))
    return dict(gamma=gamma, w=w, **out)


@cached
def d_job(args):
    """Fixed masses 1/3, the box, the noise Σ(τ): per batch the time averages of u_k and of the identity's terms −τ_kg′_k(α_k),
    g_k∂_{α_k}J_coh, g_k∂_{α_k}V (the box), Lψ_k (ψ_k = ∫g_k dα_k/G_α = (a_max²/4β_kG_α)cosh(2α_k/a_max), the generator applied to it: its
    stationary mean is the Euler chain's bias) and the chain's exact per-step bias E_n[ψ_k(α_{n+1}) − ψ_k(α_n)]/Δt − Lψ_k (the Gaussian
    increment's expectation in closed form); the martingale increments' sum of squares (the exact error of E[Lψ_k]); the share of steps with
    α_k outside the box; the amplitudes' quadratic variation with and without the drift's part."""
    name, tau = args
    t0 = time.time()
    net, base = E.setup()
    W = np.ones(3) / 3; beta = 3 * W
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    Sig = noises(tau)[name]; tk = levels(Sig)
    lam, V = np.linalg.eigh(Sig); Lc = V * np.sqrt(np.clip(lam, 0, None))            # Lc Lcᵀ = Σ
    n = int(round(S_D / DT)); nb = n // N_BATCH
    keys = ('u', 'gp', 'coh', 'box', 'lpsi', 'bias')
    acc = {k: np.zeros((N_BATCH, 3)) for k in keys}
    m2, out_box, qv, qv_noise = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)
    Th = E.exact_rep(net, W).copy(); rng = np.random.default_rng(SEED_D)
    sq = np.sqrt(DT) * Lc
    kk, cpsi, s2 = 2.0 / A_MAX, A_MAX ** 2 / (4 * beta * E.G[0]), np.diag(Sig)[0::3] * DT
    psi0 = cpsi * np.cosh(kk * Th[:, 0])
    psi_start = psi0.copy()
    for k in range(nb * N_BATCH):
        a = base.features(Th)[0]
        ap = 1 - (a / A_MAX) ** 2
        g = a / (beta * ap); gp = (1 + (a / A_MAX) ** 2) / (beta * ap)
        r = loop.residual(Th, beta)
        u = loop.fitness(r, Th)
        dT = loop.drift(Th, beta)[0]
        nonfit = -(dT[:, 0] - beta * E.G[0] * base.grad_weighted(r, Th)[:, 0]) / E.G[0]   # ∂_{α_k}(J_coh + V)
        boxt = P.SIGMA_B * loop.box_excess(Th)[:, 0] / E.G[0]                            # ∂_{α_k}V
        lpsi = u - g * nonfit + tk * gp                                                   # Lψ_k = (g_k/G_α)b_k + τ_kg′_k
        Epsi = cpsi * np.cosh(kk * (Th[:, 0] + DT * dT[:, 0])) * np.exp(0.5 * kk * kk * s2)
        i = k // nb
        acc['u'][i] += u; acc['gp'][i] += -tk * gp; acc['coh'][i] += g * (nonfit - boxt); acc['box'][i] += g * boxt
        acc['lpsi'][i] += lpsi; acc['bias'][i] += (Epsi - psi0) / DT - lpsi
        out_box += np.abs(Th[:, 0]) > P.BOX[0]
        dW = (sq @ rng.standard_normal(9)).reshape(3, 3)
        qv += (DT * dT[:, 0] + dW[:, 0]) ** 2; qv_noise += dW[:, 0] ** 2
        Th = Th + DT * dT + dW
        psi1 = cpsi * np.cosh(kk * Th[:, 0])
        m2 += (psi1 - Epsi) ** 2
        psi0 = psi1
    for k_ in keys:
        acc[k_] /= nb
    T = nb * N_BATCH * DT
    print(f"[{time.time() - t0:.0f}s] D {name} tau {tau:g}", flush=True)
    # Σ_n(E_n[ψ_k(α_(n+1))] − ψ_k(α_n)) = ψ_k(end) − ψ_k(start) − (the martingale's sum): the run's E[Lψ_k] + E[b_k] is its boundary term
    # (ψ_k(end) − ψ_k(start))/T, O(1/T), less the martingale's mean, whose error se_mart is
    return dict(name=name, tau=tau, tk=tk, **acc, se_mart=np.sqrt(m2) / T, boundary=(psi0 - psi_start) / T, out_box=out_box / (nb * N_BATCH),
                tk_hat=qv / (2 * E.G[0] * T), tk_hat_noise=qv_noise / (2 * E.G[0] * T), closed=closed(W, tk))


@cached
def e_series(args):
    """The full loop inside the box of Proposition D‴.39 (so that a stationary law exists): (variant, seed) → records (s, α_1..3, W_1..3,
    the fit error)."""
    variant, seed = args
    t0 = time.time()
    loop, kw, Th = e_loop(variant)
    rec = []
    loop.run_dithered(Th, np.ones(3), S_E, DT, TAU, np.random.default_rng(seed), rec_every=REC,
                      callback=lambda s_, T, b: rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b)]), **kw)
    print(f"[{time.time() - t0:.0f}s] E {variant} seed {seed}", flush=True)
    return dict(variant=variant, seed=seed, rec=np.array(rec))


def e_loop(variant):
    """Part E's boxed loop for a variant, run_dithered's noise keywords and the start (the exact fit at the masses 1/3) — Parts E and G."""
    net, base = E.setup()
    Th = E.exact_rep(net, np.ones(3) / 3).copy()
    A0 = base.features(Th)[0] / 3
    Sig = noises()['correlated E']
    kw = dict(mask=None, noise_L=None)
    bx = dict(box=P.BOX, box_sigma=P.SIGMA_B)
    if variant == 'nonlinear only':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, **bx); kw['mask'] = (0.0, 1.0, 1.0)
    elif variant == 'nonlinear only + barrier':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, fixman=1.0, fixman_tau=TAU, fixman_A=A0, **bx)
        kw['mask'] = (0.0, 1.0, 1.0)
    elif variant == 'correlated E':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, **bx); kw['noise_L'] = np.linalg.cholesky(Sig)
    elif variant == 'correlated E + class levels':
        loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E, fixman=2.0, fixman_tau=levels(Sig), fixman_A=A0, **bx)
        kw['noise_L'] = np.linalg.cholesky(Sig)
    else:
        raise ValueError(variant)
    return loop, kw, Th


VARIANTS_E = ('nonlinear only', 'nonlinear only + barrier', 'correlated E', 'correlated E + class levels')


def summarise_E(R):
    s = R[0, :, 0]; W = R[:, :, 4:7]
    out = dict(times=[], half=W[:, R.shape[1] // 2:].mean(axis=1), fit=R[:, R.shape[1] // 2:, 7].mean(), almax=float(np.abs(R[:, :, 1:4]).max()))
    for t in T_E:
        i = int(np.argmin(np.abs(s - t)))
        out['times'].append(dict(s=s[i], mean=W[:, i].mean(axis=0), se=W[:, i].std(axis=0, ddof=1) / np.sqrt(len(W)), std=W[:, i].std(axis=0, ddof=1)))
    return out


# ---- referee round 147 (D3ep): Part E against its predictions (F), the early phase (G), coloured noise (H) ----
T_F = (5e3, 1e4, 2e4) + tuple(T_E)
GAMMAS_F, GRID_F, DS_ODE_F, SEED_F = (1e-2, 1e-3, 1e-5), 40, 10.0, 0
LEVELS_PB = ((2.0, 1.0, 0.5), (2.0, 1.0, 0.2), (1.0, 1.0, 0.1))
SEEDS_G, S_G, REPLAY_G = tuple(range(400, 436)), 2e4, (300, 301)
VARIANTS_G = ('correlated E', 'correlated E + class levels')
TCS_H, TC_H, TAU_H, S_H, SEED_H, N_BATCH_H = (1.0, 10.0, 100.0, 1e3, 1e4), 100.0, 1e-3, 6e5, 900, 20
NAMES_H = ('isotropic', 'amplitudes only', 'correlated E')


def e_noise(variant):
    """The positions' noise Σ of a Part E variant."""
    return noises()['centres and widths only' if variant.startswith('nonlinear only') else 'correlated E']


def e_corr(variant):
    """(κ_F, the correction's levels) of a Part E variant's masses form, or None."""
    return {'nonlinear only + barrier': (1.0, np.full(3, TAU)),
            'correlated E + class levels': (2.0, levels(noises()['correlated E']))}.get(variant)


def e_target(variant):
    return {'nonlinear only': np.ones(3) / 3, 'nonlinear only + barrier': Q.w_sharp(), 'correlated E': None,
            'correlated E + class levels': w_sharp_weighted(levels(noises()['correlated E']))}[variant]


def lead_coeffs(variant):
    """(the noise's levels τ_k, (κ_F, the correction's levels τ_F,k)) of a Part E variant; κ_F = 0 without a correction."""
    c = e_corr(variant)
    return levels(e_noise(variant)), ((0.0, np.zeros(3)) if c is None else c)


def lead_fitness(variant, W, co=None):
    """The slow limit's leading-order fitness of a Part E variant at the exact fit with the masses W: the noise's −τ_kg′_k plus the masses
    form's κ_Fτ_F,kg′_k (`co`: lead_coeffs(variant), precomputed)."""
    tk, (kF, tF) = lead_coeffs(variant) if co is None else co
    return closed(W, tk) - kF * closed(W, tF)


def wall_project(W, b, tol=1e-12):
    """The mass drift b with P_B's walls (Proposition D‴.39(iii)(c)) enforced as the box does it: at a wall W_k = W_B,k a drift out of P_B
    is cancelled by raising class k's fitness, which adds λ_k(e_k − W), λ_k ≥ 0, to dW/ds."""
    W = np.asarray(W, float); b0 = np.asarray(b, float); act = []; out = b0
    while True:
        new = [k for k in range(3) if k not in act and W[k] <= P.W_BOX[k] + tol and out[k] < 0]
        if not new:
            return out
        act.append(new[0])
        lam = np.linalg.solve(np.eye(3)[np.ix_(act, act)] - W[act][:, None], -b0[act])
        out = b0 + (np.eye(3)[:, act] - W[:, None]) @ lam


def wall_push(W, f=1.0):
    """Points (n, 3) below a wall moved back along e_k − W (the box's push on class k, the others in proportion to their masses): onto the
    wall (f = 1) or reflected (f = 2)."""
    W = np.array(W, float, ndmin=2)
    for k in range(3):
        e = np.maximum(P.W_BOX[k] - W[:, k], 0.0)[:, None]
        W = W + f * e * (np.eye(3)[k][None, :] - W) / (1 - W[:, k:k + 1])
    return W


def f_ode(variant):
    """The slow limit at leading order: the replicator flow of lead_fitness from the masses 1/3, P_B's walls enforced (Euler, step
    DS_ODE_F): W at T_F, the second half's mean over the records' times (every 10³ slow units), the first time at a wall."""
    W = np.ones(3) / 3; n = int(round(S_E / DS_ODE_F)); every = int(round(1e3 / DS_ODE_F)); co = lead_coeffs(variant)
    rec, half, t_wall = {}, [], None
    for i in range(n + 1):
        s = i * DS_ODE_F
        for t in T_F:
            if abs(s - t) < DS_ODE_F / 2:
                rec[t] = W.copy()
        if s >= S_E / 2 and i % every == 0:
            half.append(W.copy())
        if t_wall is None and np.any(W <= P.W_BOX + 1e-12):
            t_wall = s
        if i == n:
            break
        b = GAMMA_E * GN.replicator_flow(W, lead_fitness(variant, W, co))
        W = wall_push(W + DS_ODE_F * wall_project(W, b), 1.0)[0]
    return dict(rec=rec, half=np.mean(half, axis=0), t_wall=t_wall)


def katz(variant, gamma):
    """Katzenberger's limit on M at the rate γ for a Part E variant (no box: inside B it changes nothing, Proposition D‴.39(iii)(b)):
    (F, chart, coords, the noise (11 × 11), F1 = the masses form's drift or None)."""
    F, chart, coords, _ = X.instance(gamma)
    S11 = np.zeros((11, 11)); S11[:9, :9] = e_noise(variant)
    c = e_corr(variant)
    if c is None:
        return F, chart, coords, S11, None
    net, base = E.setup()
    A0 = base.features(E.exact_rep(net, np.ones(3) / 3))[0] / 3
    loopc = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, fixman=c[0], fixman_tau=c[1], fixman_A=A0)

    def Fc(x):
        e = np.exp(np.r_[x[9:], 0.0]); b = 3 * e / e.sum()
        dT, db = loopc.drift(x[:9].reshape(3, 3), b); dl = db / b
        return np.concatenate([dT.ravel(), dl[:2] - dl[2]])
    return F, chart, coords, S11, (lambda x: Fc(x) - F(x))


def k_drift(Kv, W):
    F, chart, coords, S11, F1 = Kv
    return GDR.limit_drift(F, chart, coords, np.asarray(W, float)[:2], S11, h=X.H_ID, h2=X.H2_ID, F1=F1)[:2]


def grid_PB(n=GRID_F):
    """A barycentric grid on P_B (its vertices at the walls)."""
    span = 1 - P.W_BOX.sum()
    return np.array([P.W_BOX + span * np.array([i, j, n - i - j]) / n for i in range(n + 1) for j in range(n + 1 - i)])


@cached
def f_field(variant):
    """Katzenberger's drift and diffusion of (W_1, W_2) at γ = GAMMA_E on grid_PB: rows (W_1, W_2, b_1, b_2, D_11, D_12, D_22)."""
    t0 = time.time()
    Kv = katz(variant, GAMMA_E)
    out = []
    for W in grid_PB():
        b, D = k_drift(Kv, W)
        out.append(np.r_[W[:2], b, D[0, 0], D[0, 1], D[1, 1]])
    print(f"[{time.time() - t0:.0f}s] F field {variant}", flush=True)
    return np.array(out)


def f_paths(field, n_paths=None, S_end=S_E, seed=SEED_F, t_rec=T_F, diffuse=True):
    """Katzenberger's limit SDE on M from the masses 1/3 (Euler–Maruyama, step X.DS_SDE, the field interpolated linearly on grid_PB),
    reflected at P_B's walls along e_k − W; with diffuse=False the drift's path (one path, the walls enforced). Returns the records at t_rec
    (paths × 3) and each path's mean over the second half at the records' times."""
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    ip, ipn = LinearNDInterpolator(field[:, :2], field[:, 2:7]), NearestNDInterpolator(field[:, :2], field[:, 2:7])
    n_paths = (X.N_PATHS if diffuse else 1) if n_paths is None else n_paths
    rng = np.random.default_rng(seed); ds = X.DS_SDE
    W = np.tile(np.ones(3) / 3, (n_paths, 1)); rec = {}; half = np.zeros((n_paths, 3)); nh = 0
    n = int(round(S_end / ds)); every = int(round(1e3 / ds))
    for i in range(n + 1):
        s = i * ds
        for t in t_rec:
            if abs(s - t) < ds / 2:
                rec[t] = W.copy()
        if s >= S_E / 2 and i % every == 0:
            half += W; nh += 1
        if i == n:
            break
        v = ip(W[:, :2]); bad = np.isnan(v[:, 0])
        if bad.any():
            v[bad] = ipn(W[bad, :2])
        b = np.c_[v[:, :2], -v[:, :2].sum(axis=1)]
        if not diffuse:
            W = wall_push(W + ds * wall_project(W[0], b[0])[None, :], 1.0)
            continue
        L11 = np.sqrt(np.maximum(v[:, 2], 0)); L21 = v[:, 3] / np.where(L11 > 0, L11, 1); L22 = np.sqrt(np.maximum(v[:, 4] - L21 ** 2, 0))
        z = rng.standard_normal((n_paths, 2)) * np.sqrt(ds)
        dw = v[:, :2] * ds + np.c_[L11 * z[:, 0], L21 * z[:, 0] + L22 * z[:, 1]]
        W = wall_push(wall_push(np.c_[W[:, :2] + dw, 1 - W[:, :2].sum(axis=1) - dw.sum(axis=1)], 2.0), 1.0)
    return dict(rec=rec, half=half / max(nh, 1))


@cached
def f_zero(args):
    """The zero of Katzenberger's corrected drift (the masses form included) near the variant's target, at the rate γ."""
    variant, gamma = args
    from scipy.optimize import fsolve
    Kv = katz(variant, gamma); tg = e_target(variant)
    f = lambda w: k_drift(Kv, np.r_[w, 1 - sum(w)])[0] / (gamma * TAU)
    z = fsolve(f, tg[:2], xtol=1e-11)
    return dict(variant=variant, gamma=gamma, zero=np.r_[z, 1 - z.sum()], resid=float(np.abs(f(z)).max()), at_target=f(tg[:2]))


def at_third(variant):
    """The slow limit's and Katzenberger's (γ = GAMMA_E) mass drifts at W = 1/3, per γτ."""
    b13 = k_drift(katz(variant, GAMMA_E), np.ones(3) / 3)[0]
    return dict(lead13=GN.replicator_flow(np.ones(3) / 3, lead_fitness(variant, np.ones(3) / 3)) / TAU,
                b13=np.r_[b13, -b13.sum()] / (GAMMA_E * TAU))


def part_F(Erec):
    """Part E's variants against the slow limit's leading-order ODE and Katzenberger's limit at γ = GAMMA_E (the drift's path and the SDE's
    ensemble); the corrected drift's zeros; W♯_τ against P_B."""
    out = dict(var={}, zeros=[], pb=[])
    for v in VARIANTS_E:
        fld = f_field(v)
        out['var'][v] = dict(ode=f_ode(v), kpath=f_paths(fld, diffuse=False), sde=f_paths(fld), field=fld, **at_third(v))
        print(f"F {v}", flush=True)
    out['zeros'] = [f_zero((v, g)) for v in ('nonlinear only + barrier', 'correlated E + class levels') for g in GAMMAS_F]
    for lv in LEVELS_PB:
        w = w_sharp_weighted(lv)
        out['pb'].append(dict(levels=lv, w=w, margin=w - P.W_BOX))
    return out


@cached
def g_series(args):
    """A Part E variant's full loop from the masses 1/3 over S slow units with the generator's seed (Part G: the early phase)."""
    variant, seed, S = args
    loop, kw, Th = e_loop(variant)
    rec = []
    loop.run_dithered(Th, np.ones(3), S, DT, TAU, np.random.default_rng(seed), rec_every=REC,
                      callback=lambda s_, T, b: rec.append(np.r_[s_, T[:, 0], b / 3, loop.fit_error(T, b)]), **kw)
    return dict(variant=variant, seed=seed, rec=np.array(rec))


def lin_GH(W=np.ones(3) / 3, hH=1e-5):
    """G∇²_ΘJ at the exact fit with the masses W (the positions' linearisation at fixed masses)."""
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=1e-2)
    beta = 3 * np.asarray(W, float)
    gradJ = lambda t: -loop.drift(t.reshape(3, 3), beta)[0].ravel() / G9
    th = E.exact_rep(net, np.asarray(W, float)).ravel()
    H = np.array([(gradJ(th + hH * e) - gradJ(th - hH * e)) / (2 * hH) for e in np.eye(9)]).T
    return G9[:, None] * (0.5 * (H + H.T))


def part_H_linear():
    """The coloured noise's levels at linear response (fixed masses 1/3): τ_k(τ_c) = [(I + τ_cGH)⁻¹Σ]_{α_kα_k}/(2G_α) over the white τ_k,
    with Σ and with Σ's diagonal only; the closed form against the augmented system's Lyapunov solution."""
    GH = lin_GH()
    out = dict(rates=np.sort(np.linalg.eigvals(GH).real), rows=[], check=0.0)
    for name in NAMES_H:
        S = noises()[name]
        for tc in TCS_H:
            full = np.diag(GN.coloured_cross(GH, S, tc))[0::3] / E.G[0]
            diag = np.diag(GN.coloured_cross(GH, np.diag(np.diag(S)), tc))[0::3] / E.G[0]
            aug = np.diag(GN.coloured_cross_augmented(GH, S, tc))[0::3] / E.G[0]
            out['check'] = max(out['check'], float(np.abs(aug / full - 1).max()))
            out['rows'].append(dict(name=name, tc=tc, ratio=full / levels(S), corr=np.abs(full - diag) / full))
    return out


def part_H(mc):
    """The linear response's levels and the simulations' leading-order comparison (amplitudes only at TAU_H, masses 1/3)."""
    lin = part_H_linear()
    W = np.ones(3) / 3; tk = levels(noises(TAU_H)['amplitudes only'])
    ratio = [r['ratio'] for r in lin['rows'] if r['name'] == 'amplitudes only' and r['tc'] == TC_H][0]
    return dict(linear=lin, mc=list(mc), lead_white=closed(W, tk), lead_col=closed(W, ratio * tk), ratio_col=ratio)


@cached
def h_mc(args):
    """Fixed masses 1/3 inside the box, the amplitudes-only noise at the level TAU_H: white (tc = 0, Euler–Maruyama as in Part D) or an
    Ornstein–Uhlenbeck noise η of correlation time tc with the same integrated covariance Σ (dϑ = −G∇(J + V)ds + η ds, η's covariance
    Σe^{−|s−s′|/tc}/(2tc); the exact AR(1) step, started stationary). Per batch the time averages of u_k, g_k∂_{α_k}J_coh, g_k∂_{α_k}V,
    −τ_kg′_k and −g_kη_{α_k}/G_α (dψ_k/ds = u_k − g_k∂_{α_k}(J_coh + V) + g_kη_{α_k}/G_α for the coloured noise: E[u_k] = −E[g_kη_{α_k}]/G_α +
    E[g_k∂_{α_k}(J_coh + V)] in its stationary law); the amplitudes' quadratic variation over the run."""
    tc, seed, S, nbatch = args
    t0 = time.time()
    net, base = E.setup()
    W = np.ones(3) / 3; beta = 3 * W
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B)
    Sig = noises(TAU_H)['amplitudes only']; tk = levels(Sig); sd = np.sqrt(np.diag(Sig))
    n = int(round(S / DT)); nb = n // nbatch
    keys = ('u', 'coh', 'box', 'gp', 'geta')
    acc = {k: np.zeros((nbatch, 3)) for k in keys}
    qv = np.zeros(3)
    Th = E.exact_rep(net, W).copy(); rng = np.random.default_rng(seed)
    if tc > 0:
        eta = sd / np.sqrt(2 * tc) * rng.standard_normal(9); ar = np.exp(-DT / tc); s_ar = np.sqrt(1 - ar * ar) * sd / np.sqrt(2 * tc)
    else:
        eta = np.zeros(9)
    for k in range(nb * nbatch):
        a = base.features(Th)[0]
        ap = 1 - (a / A_MAX) ** 2
        g = a / (beta * ap); gp = (1 + (a / A_MAX) ** 2) / (beta * ap)
        r = loop.residual(Th, beta)
        u = loop.fitness(r, Th)
        dT = loop.drift(Th, beta)[0]
        nonfit = -(dT[:, 0] - beta * E.G[0] * base.grad_weighted(r, Th)[:, 0]) / E.G[0]   # ∂_{α_k}(J_coh + V)
        boxt = P.SIGMA_B * loop.box_excess(Th)[:, 0] / E.G[0]                            # ∂_{α_k}V
        i = k // nb
        acc['u'][i] += u; acc['coh'][i] += g * (nonfit - boxt); acc['box'][i] += g * boxt; acc['gp'][i] += -tk * gp
        acc['geta'][i] += -g * eta.reshape(3, 3)[:, 0] / E.G[0]
        xi = rng.standard_normal(9)
        if tc > 0:
            dTh = DT * (dT + eta.reshape(3, 3))
            eta = ar * eta + s_ar * xi
        else:
            dTh = DT * dT + (np.sqrt(DT) * sd * xi).reshape(3, 3)
        qv += dTh[:, 0] ** 2
        Th = Th + dTh
    for k_ in keys:
        acc[k_] /= nb
    print(f"[{time.time() - t0:.0f}s] H tc {tc:g}", flush=True)
    return dict(tc=tc, seed=seed, S=S, tk=tk, **acc, tk_hat=qv / (2 * E.G[0] * nb * nbatch * DT))


# ---- the early phase's gap (round 147's open item): Part I ----
SEEDS_I, B_I, S_I, BURN_I, BURN_OFF_I = tuple(range(500, 692)), 48, 2e4, 1e4, 10 ** 5
KINDS_I = ('box', 'no box', 'burn-in')
T_I = (5e3, 1e4, 2e4)
ACC_I = ('bdry', 'mart', 'noise', 'box', 'coh', 'bias')
W1_FIX_I, SEEDS_FIX_I, BURN_FIX_I, S_FIX_I = (0.2947, 0.2975, 0.3, 0.3025, 0.305, 0.3075, 0.31, 0.3125, 0.315, 0.32), 12, 5e3, 2e4
L_EXT_I = np.abs(E.A_BUMP) / A_MAX + 0.012
DS_I, SEED_SDE_I = (50.0, 10.0, 2.0), 0
DS_SOFT_I = (2.0, 0.5)
JAC13 = (np.eye(3) - 1 / 3) * (1 / 3)                         # ∂W/∂log β at the masses 1/3: W_k(δ_kl − W_l)


def i_bd(box, gamma=GAMMA_E):
    """Part E's 'correlated E' loop (uncorrected) as a gauge_batch.BatchDither, with or without the box."""
    from src.core import gauge_batch as GB
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma, **(dict(box=P.BOX, box_sigma=P.SIGMA_B) if box else {}))
    return net, GB.BatchDither(loop)


@cached
def i_runs(args):
    """Part E's 'correlated E' loop (uncorrected, γ = GAMMA_E, τ = TAU) for a batch of seeds over S slow units, integrated together
    (`gauge_batch.BatchDither`'s drift; each run draws run_dithered's noise √dt·Lξ from its own generator, so a run is run_dithered's run
    with that seed) from the exact fit at the masses 1/3. kind 'box': Part E's loop; 'no box': without the box; 'burn-in': the positions
    first equilibrated over BURN_I slow units at the fixed masses 1/3 (γ = 0, the box, a generator seeded seed + BURN_OFF_I), then the loop
    with the seed's generator — the same draws as 'box' from the start of the masses' motion. Every REC steps: s, α_1..3, W_1..3, and the
    gauge identity's decomposition accumulated per class: with φ_k = (a_max²/4G_α)cosh(2α_k/a_max) and the Euler chain's closed-form
    E_n[φ_k(α_(n+1))], u_kΔt = (φ_k(α_(n+1)) − φ_k(α_n))/β_k − (φ_k(α_(n+1)) − E_nφ_k(α_(n+1)))/β_k + [−τ_kg′_k + g_k∂_{α_k}V + g_k∂_{α_k}J_coh]Δt
    + the chain's bias Δt: ACC_I = the boundary term, the martingale (mean zero exactly), the noise's drift, the box's, the cohesion's, the
    bias; Δlog(β_k/β_3) = γΣ(u_k − u_3)Δt exactly."""
    kind, seeds, S = args
    t0 = time.time()
    net, bd = i_bd(kind != 'no box')
    Sig = noises()['correlated E']; Ld = np.sqrt(DT) * np.linalg.cholesky(Sig); tk = levels(Sig)
    B = len(seeds); kk = 2.0 / A_MAX; cphi = A_MAX ** 2 / (4 * E.G[0]); s2 = np.diag(Sig)[0::3] * DT
    Th = np.tile(E.exact_rep(net, np.ones(3) / 3), (B, 1, 1)); lb = np.zeros((B, 3))
    if kind == 'burn-in':
        _, b0 = i_bd(True, 0.0); g0 = [np.random.default_rng(s + BURN_OFF_I) for s in seeds]
        for k in range(int(round(BURN_I / DT))):
            if k % REC == 0:
                Z = np.array([g.standard_normal((REC, 9)) for g in g0])
            Th = Th + DT * b0.drift(Th, np.ones((B, 3)))[0] + (Z[:, k % REC] @ Ld.T).reshape(B, 3, 3)
    gens = [np.random.default_rng(s) for s in seeds]
    n = int(round(S / DT)); acc = np.zeros((B, len(ACC_I), 3)); s_, al, Wr, ac = [], [], [], []
    for k in range(n + 1):
        b = np.exp(lb)
        if k % REC == 0:
            s_.append(k * DT); al.append(Th[:, :, 0].copy()); Wr.append(b / 3); ac.append(acc.copy())
        if k == n:
            break
        if k % REC == 0:
            Z = np.array([g.standard_normal((REC, 9)) for g in gens])
        dT, db = bd.drift(Th, b)
        a, c, sw, diff, phi = bd.features(Th)
        r = bd.F[None, :] - np.matmul((b * a)[:, None, :], phi)[:, 0, :] / 3
        pr = np.matmul(phi, r[:, :, None])[:, :, 0]
        q = (a / A_MAX) ** 2
        u = a * pr / bd.N; g = a / (b * (1 - q)); gp = (1 + q) / (b * (1 - q))
        nonfit = -(dT[:, :, 0] - b * E.G[0] * pr * (1 - q) / bd.N) / E.G[0]                     # ∂_{α_k}(J_coh + V)
        boxt = (P.SIGMA_B * (Th[:, :, 0] - np.clip(Th[:, :, 0], -P.BOX[0], P.BOX[0])) / E.G[0]) if bd.box is not None else 0.0 * nonfit
        phi0 = cphi * np.cosh(kk * Th[:, :, 0]); Ephi = cphi * np.cosh(kk * (Th[:, :, 0] + DT * dT[:, :, 0])) * np.exp(0.5 * kk * kk * s2)
        Th = Th + DT * dT + (Z[:, k % REC] @ Ld.T).reshape(B, 3, 3)
        lb = lb + DT * db / b
        lb = lb - np.log(np.exp(lb).sum(axis=1, keepdims=True) / 3)
        phi1 = cphi * np.cosh(kk * Th[:, :, 0])
        drift = (-tk * gp, g * boxt, g * (nonfit - boxt))
        acc[:, 0] += (phi1 - phi0) / b; acc[:, 1] += (phi1 - Ephi) / b
        for j, d in enumerate(drift):
            acc[:, 2 + j] += d * DT
        acc[:, 5] += u * DT - (Ephi - phi0) / b - sum(drift) * DT
    print(f"[{time.time() - t0:.0f}s] I {kind} seeds {seeds[0]}–{seeds[-1]}", flush=True)
    return dict(kind=kind, seeds=tuple(seeds), s=np.array(s_), al=np.array(al).transpose(1, 0, 2), W=np.array(Wr).transpose(1, 0, 2),
                acc=np.array(ac).transpose(1, 0, 2, 3))


@cached
def i_fixed(args):
    """The box's term at fixed masses (γ = 0, the box, 'correlated E'): for each W_1 in w1s (W_2 = W_3 = (1 − W_1)/2) `ns` runs (seeds
    7000 + i) from the exact fit, burn-in, then the time averages of g_k∂_{α_k}V (= g_kσ_B(α_k − Π_Bα_k)/G_α), of −τ_kg′_k, and the share
    of steps with α_k outside the box."""
    w1s, ns, burn, S = args
    t0 = time.time()
    net, bd = i_bd(True, 0.0)
    Sig = noises()['correlated E']; Ld = np.sqrt(DT) * np.linalg.cholesky(Sig); tk = levels(Sig)
    Ws = np.array([[w, (1 - w) / 2, (1 - w) / 2] for w in w1s for _ in range(ns)]); B = len(Ws); beta = 3 * Ws
    Th = np.array([E.exact_rep(net, W) for W in Ws]); gens = [np.random.default_rng(7000 + i) for i in range(B)]
    nb, nm = int(round(burn / DT)), int(round(S / DT)); box, gpa, out = np.zeros((B, 3)), np.zeros((B, 3)), np.zeros((B, 3))
    for k in range(nb + nm):
        if k % REC == 0:
            Z = np.array([g.standard_normal((REC, 9)) for g in gens])
        dT = bd.drift(Th, beta)[0]
        if k >= nb:
            al = Th[:, :, 0]; a = A_MAX * np.tanh(al / A_MAX); q = (a / A_MAX) ** 2; ex = al - np.clip(al, -P.BOX[0], P.BOX[0])
            box += a / (beta * (1 - q)) * P.SIGMA_B * ex / E.G[0]; gpa += -tk * (1 + q) / (beta * (1 - q)); out += ex != 0
        Th = Th + DT * dT + (Z[:, k % REC] @ Ld.T).reshape(B, 3, 3)
    print(f"[{time.time() - t0:.0f}s] I fixed masses", flush=True)
    return dict(w1=np.array(w1s), ns=ns, W=Ws, box=box / nm, gp=gpa / nm, out=out / nm)


def grid_ext(n=GRID_F):
    """A barycentric grid on {W_k ≥ L_EXT_I,k} ⊃ P_B (0.012 above the saturations W_min)."""
    span = 1 - L_EXT_I.sum()
    return np.array([L_EXT_I + span * np.array([i, j, n - i - j]) / n for i in range(n + 1) for j in range(n + 1 - i)])


@cached
def i_field(variant):
    """Katzenberger's drift and diffusion at γ = GAMMA_E on grid_ext (rows W_1, W_2, b_1, b_2, D_11, D_12, D_22) and the projection Π (2 × 3,
    row-major) of the box's term: Π_{·j} the drift of (W_1, W_2) when class j's fitness is raised by 1 through a force on α_j alone
    (−G_α/g_j at the fit: limit_drift's F1 with Σ = 0); in the slow limit Π_{kj} = γW_k(δ_kj − W_j)."""
    t0 = time.time()
    Kv = katz(variant, GAMMA_E)
    out = np.array([i_field_row(Kv, W) for W in grid_ext()])
    print(f"[{time.time() - t0:.0f}s] I field", flush=True)
    return out


def i_field_row(Kv, W):
    """One row of i_field at the masses W (Kv = katz(variant, γ))."""
    F, chart, coords = Kv[:3]
    net, base = E.setup()
    b, D = k_drift(Kv, W)
    a = base.features(E.exact_rep(net, W))[0]; gk = a / (3 * W * (1 - (a / A_MAX) ** 2)); Pi = []
    for j in range(3):
        f = np.zeros(11); f[3 * j] = -E.G[0] / gk[j]
        Pi.append(GDR.limit_drift(F, chart, coords, W[:2], np.zeros((11, 11)), h=X.H_ID, h2=X.H2_ID, F1=lambda x, f=f: f)[0])
    return np.r_[W[:2], b, D[0, 0], D[0, 1], D[1, 1], np.array(Pi).T.ravel()]


def i_cbox(fx):
    """W_1 ↦ class 1's box term from Part I's fixed-mass runs: log-linear between the measured W_1 (the runs' means), extrapolated with
    the end slopes; classes 2 and 3 zero (their amplitudes stay far inside the box)."""
    m = fx['box'][:, 0].reshape(len(fx['w1']), fx['ns']).mean(axis=1); ok = m > 0
    x, y = fx['w1'][ok], np.log(m[ok])

    def c(W1):
        W1 = np.asarray(W1, float)
        lo, hi = (y[1] - y[0]) / (x[1] - x[0]), (y[-1] - y[-2]) / (x[-1] - x[-2])
        v = np.interp(W1, x, y) + np.where(W1 < x[0], (W1 - x[0]) * lo, 0.0) + np.where(W1 > x[-1], (W1 - x[-1]) * hi, 0.0)
        return np.exp(v)
    return c


def i_sde(field, walls, ds, cbox=None, n_paths=None, seed=SEED_SDE_I, S_end=S_I, t_rec=T_I):
    """Katzenberger's SDE (as f_paths: Euler–Maruyama, the field interpolated linearly, from the masses 1/3) reflected at `walls` (P_B's
    W_B, or grid_ext's lower ends: no walls inside it); with `cbox` (W_1 ↦ class 1's box term) the drift adds Π(W)·(cbox(W_1), 0, 0), Π
    interpolated from field's columns 7–12. Returns the records at t_rec (paths × 3) and the share of paths that reached a wall."""
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    cols = field[:, 2:] if cbox is not None else field[:, 2:7]
    ip, ipn = LinearNDInterpolator(field[:, :2], cols), NearestNDInterpolator(field[:, :2], cols)
    n_paths = X.N_PATHS if n_paths is None else n_paths
    rng = np.random.default_rng(seed); W = np.tile(np.ones(3) / 3, (n_paths, 1)); rec = {}; hit = np.zeros(n_paths, bool)
    walls = np.asarray(walls, float); n = int(round(S_end / ds))
    for i in range(n + 1):
        for t in t_rec:
            if abs(i * ds - t) < ds / 2:
                rec[t] = W.copy()
        if i == n:
            break
        v = ip(W[:, :2]); bad = np.isnan(v[:, 0])
        if bad.any():
            v[bad] = ipn(W[bad, :2])
        b = v[:, :2]
        if cbox is not None:
            b = b + v[:, 5:11].reshape(-1, 2, 3)[:, :, 0] * cbox(W[:, 0])[:, None]
        L11 = np.sqrt(np.maximum(v[:, 2], 0)); L21 = v[:, 3] / np.where(L11 > 0, L11, 1); L22 = np.sqrt(np.maximum(v[:, 4] - L21 ** 2, 0))
        z = rng.standard_normal((n_paths, 2)) * np.sqrt(ds)
        dw = b * ds + np.c_[L11 * z[:, 0], L21 * z[:, 0] + L22 * z[:, 1]]
        W = np.c_[W[:, :2] + dw, 1 - W[:, :2].sum(axis=1) - dw.sum(axis=1)]
        hit |= np.any(W < walls, axis=1)
        for f in (2.0, 1.0):                                   # wall_push at these walls: reflected (f = 2), then onto the wall (f = 1)
            for k in range(3):
                e = np.maximum(walls[k] - W[:, k], 0.0)[:, None]
                W = W + f * e * (np.eye(3)[k][None, :] - W) / (1 - W[:, k:k + 1])
    return dict(rec=rec, hit=float(hit.mean()))


def i_transient():
    """The positions' spreading at fixed masses 1/3 (Part E's start, the exact fit) as a mass offset: γ(E[φ_k] − φ_k(fit))/β_k at the
    Gaussian law's mean shift μ and covariance C (referee round 147's estimate), mapped to W by JAC13."""
    net, base = E.setup(); W = np.ones(3) / 3; beta = 3 * W
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA_E)
    gradJ = lambda t: -loop.drift(t.reshape(3, 3), beta)[0].ravel() / G9
    u = lambda t: loop.fitness(loop.residual(t.reshape(3, 3), beta), t.reshape(3, 3))
    th = E.exact_rep(net, W).ravel()
    _, C, mu = GN.averaged_fitness(gradJ, u, th, noises()['correlated E'], G9)
    al = th[0::3]; cphi = A_MAX ** 2 / (4 * E.G[0])
    dphi = cphi * ((2 / A_MAX) * np.sinh(2 * al / A_MAX) * mu[0::3] + 0.5 * (2 / A_MAX) ** 2 * np.cosh(2 * al / A_MAX) * np.diag(C)[0::3])
    return JAC13 @ (GAMMA_E * dphi / beta)


def part_I(F):
    """The early phase's gap: the loop's runs (I_runs), the fixed-mass box term, the extended field, the SDE variants."""
    out = dict(runs={}, replay={})
    fr = [SEEDS_I[i:i + B_I] for i in range(0, len(SEEDS_I), B_I)]
    for kind in KINDS_I:
        rs = [i_runs((kind, sd, S_I)) for sd in fr]
        out['runs'][kind] = {k: np.concatenate([r[k] for r in rs]) for k in ('al', 'W', 'acc')}
        out['runs'][kind]['s'] = rs[0]['s']
    old = i_runs(('box', SEEDS_G + SEEDS_E, S_I))
    out['runs']['old'] = {k: old[k] for k in ('s', 'al', 'W', 'acc')}
    out['fixed'] = i_fixed((W1_FIX_I, SEEDS_FIX_I, BURN_FIX_I, S_FIX_I))
    out['field'] = i_field('correlated E')
    out['transient'] = i_transient()
    fp = F['var']['correlated E']['field']; cb = i_cbox(out['fixed'])
    out['sde'] = {f'reflected at W_B, step {ds:g}': i_sde(fp, P.W_BOX, ds) for ds in DS_I}
    out['sde']['no walls, step 10'] = i_sde(out['field'], L_EXT_I, 10.0)
    for ds in DS_SOFT_I:
        out['sde'][f'the box measured, step {ds:g}'] = i_sde(out['field'], L_EXT_I, ds, cbox=cb)
    return out


# ---- referee round 153 (D3ev): the identity's own check (J), the Euler chain's bias along Part E (K), errors and extremes (L) ----
S_J, SEED_J, N_BATCH_J = 2e5, 710, 10
RUNS_J = RUNS_D
TERMS_J = ('u', 'gp', 'coh', 'box', 'lpsi', 'b', 'b_alt')
VARIANTS_K = ('nonlinear only', 'nonlinear only + barrier')
ACC_K = ('bdry', 'box', 'coh', 'chain', 'corr')


def gauge_terms(bd, Th, b, dT):
    """The gauge identity's terms at a batch (Θ (B, 3, 3), β (B, 3), its drift dΘ), as in i_runs: a_k, u_k, g_k, g′_k, ∂_{α_k}(J_coh + V),
    ∂_{α_k}V."""
    a, c, sw, diff, phi = bd.features(Th)
    r = bd.F[None, :] - np.matmul((b * a)[:, None, :], phi)[:, 0, :] / 3
    pr = np.matmul(phi, r[:, :, None])[:, :, 0]
    q = (a / A_MAX) ** 2
    u = a * pr / bd.N; g = a / (b * (1 - q)); gp = (1 + q) / (b * (1 - q))
    nonfit = -(dT[:, :, 0] - b * E.G[0] * pr * (1 - q) / bd.N) / E.G[0]
    boxt = (P.SIGMA_B * (Th[:, :, 0] - np.clip(Th[:, :, 0], -P.BOX[0], P.BOX[0])) / E.G[0]) if bd.box is not None else 0.0 * nonfit
    return a, u, g, gp, nonfit, boxt


@cached
def j_runs(args):
    """Part D's per-state check at two steps (referee round 153): the runs RUNS_J at the fixed masses 1/3 with the box, each as two
    Euler–Maruyama chains driven by one Brownian path — the step Δt = DT ('full') and Δt/2 ('half': the half steps' increments
    √(Δt/2)Lz₁, √(Δt/2)Lz₂; the full step's √(Δt/2)L(z₁ + z₂)). At every state of each chain the identity's terms as in d_job and the
    chain's exact per-step bias b_k(h) = (E_n[ψ_k(α_k + hb_α + ν)] − ψ_k(α_k))/h − Lψ_k (ν the amplitude's Gaussian increment, of variance
    Σ_{α_kα_k}h; the expectation in closed form) at the chain's own step ('b') and at the other chain's ('b_alt': Δt at the half chain's
    states, Δt/2 at the full chain's); per chain the batch means over N_BATCH_J batches, the martingale increments' sum of squares, the
    boundary term and the share of states with α_k outside the box."""
    S, seed = args
    t0 = time.time()
    from src.core import gauge_batch as GB
    net, base = E.setup()
    bd = GB.BatchDither(MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0, box=P.BOX, box_sigma=P.SIGMA_B))
    B = len(RUNS_J); Sig = [noises(t)[nm] for nm, t in RUNS_J]
    tk = np.array([levels(s) for s in Sig]); var = np.array([np.diag(s)[0::3] for s in Sig])
    Lc = np.array([V * np.sqrt(np.clip(lam, 0, None)) for lam, V in (np.linalg.eigh(s) for s in Sig)])     # Lc Lcᵀ = Σ, per run
    h = DT / 2; kk, cpsi = 2.0 / A_MAX, A_MAX ** 2 / (4 * E.G[0])                       # ψ_k = cpsi·cosh(kk α_k)/β_k, β_k = 1
    tk2, var2, beta2 = np.r_[tk, tk], np.r_[var, var], np.ones((2 * B, 3))               # rows: the full chain, then the half chain
    own2 = np.r_[np.full(B, DT), np.full(B, h)][:, None]; alt2 = np.r_[np.full(B, h), np.full(B, DT)][:, None]
    F_, H_ = slice(0, B), slice(B, 2 * B)

    def terms(T, sl):
        dT = bd.drift(T, beta2[sl])[0]
        _, u, g, gp, nonfit, boxt = gauge_terms(bd, T, beta2[sl], dT)
        lpsi = u - g * nonfit + tk2[sl] * gp
        psi = cpsi * np.cosh(kk * T[:, :, 0])
        Ep = lambda hh: cpsi * np.cosh(kk * (T[:, :, 0] + hh * dT[:, :, 0])) * np.exp(0.5 * kk * kk * var2[sl] * hh)
        Eo = Ep(own2[sl])
        return dT, dict(u=u, gp=-tk2[sl] * gp, coh=g * (nonfit - boxt), box=g * boxt, lpsi=lpsi, b=(Eo - psi) / own2[sl] - lpsi,
                        b_alt=(Ep(alt2[sl]) - psi) / alt2[sl] - lpsi), Eo

    n = int(round(S / DT)); nb = n // N_BATCH_J
    Th = np.tile(E.exact_rep(net, np.ones(3) / 3), (2 * B, 1, 1))
    acc = {k: np.zeros((N_BATCH_J, 2 * B, 3)) for k in TERMS_J}
    m2, out = np.zeros((2 * B, 3)), np.zeros((2 * B, 3))
    psi_start = cpsi * np.cosh(kk * Th[:, :, 0])
    rng = np.random.default_rng(seed); sq = np.sqrt(h)
    noise = lambda z: sq * np.einsum('bij,bj->bi', Lc, z).reshape(B, 3, 3)
    for k in range(nb * N_BATCH_J):
        if k % REC == 0:
            Z = rng.standard_normal((REC, 2, B, 9))
        i = k // nb; z1, z2 = Z[k % REC]
        dT, t, Eo = terms(Th, slice(None))                                   # both chains (the half chain at its first half step)
        for kn in TERMS_J:
            acc[kn][i] += t[kn]
        out += np.abs(Th[:, :, 0]) > P.BOX[0]
        Th = Th.copy()
        Th[F_] = Th[F_] + DT * dT[F_] + noise(z1 + z2)
        Th[H_] = Th[H_] + h * dT[H_] + noise(z1)
        m2 += (cpsi * np.cosh(kk * Th[:, :, 0]) - Eo) ** 2
        dT, t, Eo = terms(Th[H_], H_)                                       # the half chain's second half step
        for kn in TERMS_J:
            acc[kn][i, H_] += t[kn]
        out[H_] += np.abs(Th[H_, :, 0]) > P.BOX[0]
        Th[H_] = Th[H_] + h * dT + noise(z2)
        m2[H_] += (cpsi * np.cosh(kk * Th[H_, :, 0]) - Eo) ** 2
    T = nb * N_BATCH_J * DT; cnt = np.r_[np.full(B, nb), np.full(B, 2 * nb)][:, None]
    for kn in TERMS_J:
        acc[kn] /= cnt[None]
    psi_end = cpsi * np.cosh(kk * Th[:, :, 0])
    print(f"[{time.time() - t0:.0f}s] J", flush=True)
    ch = lambda x, sl: x[:, sl] if x.ndim == 3 else x[sl]
    return dict(runs=RUNS_J, seed=seed, S=S, tk=tk, dt=(DT, h),
                **{c: dict(**{kn: ch(acc[kn], sl) for kn in TERMS_J}, se_mart=np.sqrt(m2[sl]) / T, boundary=(psi_end - psi_start)[sl] / T,
                           out_box=out[sl] / (cnt[sl] * N_BATCH_J)) for c, sl in (('full', F_), ('half', H_))})


def ratio_se(x, y):
    """The ratio of the means of the paired batch means x, y (batches × …) and its delta-method s.e. from the batches' x − r·y."""
    r = x.mean(axis=0) / y.mean(axis=0)
    return r, (x - r * y).std(axis=0, ddof=1) / np.sqrt(len(x)) / np.abs(y.mean(axis=0))


@cached
def k_runs(args):
    """Part E's deterministic-amplitude variants replayed over S slow units with the gauge identity's decomposition (referee round 153;
    `gauge_batch.BatchDither`'s drift, the variants in one batch through a per-run gain κ_F, each seed run_dithered's draws, so that the
    records are Part E's): with φ_k = (a_max²/4G_α)cosh(2α_k/a_max) and the amplitudes noise-free (α_(n+1) = α_n + Δt b_α exactly),
    u_kΔt = (φ_k(α_(n+1)) − φ_k(α_n))/β_k + (g_k∂_{α_k}V + g_k∂_{α_k}J_coh)Δt + the chain's term, which is −(φ_k(α_n + Δt b_α) − φ_k(α_n) −
    Δtφ′_k(α_n)b_α)/β_k, minus Δt times the chain's per-step bias of Part D; the masses form's fitness u^F_kΔt apart. ACC_K per class every
    REC steps, with the positions; Δlog(β_k/β_3) = γΣ_n((u + u^F)_k − (u + u^F)_3)Δt exactly."""
    variants, seeds, S = args
    t0 = time.time()
    from src.core import gauge_batch as GB
    loop, kw, Th0 = e_loop('nonlinear only + barrier')
    for v in variants:
        kv = e_loop(v)[1]
        assert kv['noise_L'] is None and kv['mask'] == kw['mask'] and kv['mask'][0] == 0.0
    runs = [(v, s) for v in variants for s in seeds]; B = len(runs)
    bd = GB.BatchDither(loop, fixman=np.array([{'nonlinear only': 0.0, 'nonlinear only + barrier': 1.0}[v] for v, _ in runs]))
    sq = np.sqrt(2.0 * TAU * DT * np.asarray(loop.G, float)) * np.asarray(kw['mask'], float)
    kk, cphi = 2.0 / A_MAX, A_MAX ** 2 / (4 * E.G[0])
    Th = np.tile(Th0, (B, 1, 1)); lb = np.zeros((B, 3)); gens = [np.random.default_rng(s) for _, s in runs]
    n = int(round(S / DT)); acc = np.zeros((B, len(ACC_K), 3)); s_, al, Wr, ac, Tr = [], [], [], [], []
    for k in range(n + 1):
        b = np.exp(lb)
        if k % REC == 0:
            s_.append(k * DT); al.append(Th[:, :, 0].copy()); Wr.append(b / 3); ac.append(acc.copy()); Tr.append(Th.copy())
        if k == n:
            break
        if k % REC == 0:
            Z = np.array([g_.standard_normal((REC, 3, 3)) for g_ in gens])
        dT, db = bd.drift(Th, b)
        a, u, g, gp, nonfit, boxt = gauge_terms(bd, Th, b, dT)
        uF = bd.correction(a, b)
        phi0 = cphi * np.cosh(kk * Th[:, :, 0])
        Th = Th + DT * dT + sq[None, None, :] * Z[:, k % REC]
        lb = lb + DT * db / b
        lb = lb - np.log(np.exp(lb).sum(axis=1, keepdims=True) / 3)
        dphi = (cphi * np.cosh(kk * Th[:, :, 0]) - phi0) / b; bx, co = g * boxt * DT, g * (nonfit - boxt) * DT
        acc[:, 0] += dphi; acc[:, 1] += bx; acc[:, 2] += co; acc[:, 3] += u * DT - dphi - bx - co; acc[:, 4] += uF * DT
    print(f"[{time.time() - t0:.0f}s] K", flush=True)
    tr = lambda x: np.array(x).swapaxes(0, 1)
    out = dict(seeds=tuple(seeds), S=S, s=np.array(s_), var={})
    for v in variants:
        ix = [i for i, (w, _) in enumerate(runs) if w == v]
        out['var'][v] = dict(al=tr(al)[ix], W=tr(Wr)[ix], acc=tr(ac)[ix], Th=tr(Tr)[ix])
    return out


S_KH, W_OFF_KH = 2e5, 10 ** 6


@cached
def k_half(args):
    """Part K's runs at the step Δt/2 over S slow units, coupled to Part E's (referee round 153): each full step's increment of Part E's
    run, sq·z (z run_dithered's draw), is split by the Brownian bridge into the half steps' increments sq·(z + w)/2 and sq·(z − w)/2, w
    standard normal from a second generator seeded seed + W_OFF_KH, so that the two chains share the Brownian path at the full steps. The
    decomposition of k_runs, with the step Δt/2, every REC full steps."""
    variants, seeds, S = args
    t0 = time.time()
    from src.core import gauge_batch as GB
    loop, kw, Th0 = e_loop('nonlinear only + barrier')
    runs = [(v, s) for v in variants for s in seeds]; B = len(runs)
    bd = GB.BatchDither(loop, fixman=np.array([{'nonlinear only': 0.0, 'nonlinear only + barrier': 1.0}[v] for v, _ in runs]))
    sq = np.sqrt(2.0 * TAU * DT * np.asarray(loop.G, float)) * np.asarray(kw['mask'], float)
    kk, cphi, h = 2.0 / A_MAX, A_MAX ** 2 / (4 * E.G[0]), DT / 2
    Th = np.tile(Th0, (B, 1, 1)); lb = np.zeros((B, 3))
    gens = [np.random.default_rng(s) for _, s in runs]; gw = [np.random.default_rng(s + W_OFF_KH) for _, s in runs]
    n = int(round(S / DT)); acc = np.zeros((B, len(ACC_K), 3)); s_, al, Wr, ac = [], [], [], []
    for k in range(n + 1):
        if k % REC == 0:
            s_.append(k * DT); al.append(Th[:, :, 0].copy()); Wr.append(np.exp(lb) / 3); ac.append(acc.copy())
        if k == n:
            break
        if k % REC == 0:
            Z = np.array([g_.standard_normal((REC, 3, 3)) for g_ in gens]); Zw = np.array([g_.standard_normal((REC, 3, 3)) for g_ in gw])
        z, w = Z[:, k % REC], Zw[:, k % REC]
        for dz in (z + w, z - w):
            b = np.exp(lb)
            dT, db = bd.drift(Th, b)
            a, u, g, gp, nonfit, boxt = gauge_terms(bd, Th, b, dT)
            uF = bd.correction(a, b)
            phi0 = cphi * np.cosh(kk * Th[:, :, 0])
            Th = Th + h * dT + 0.5 * sq[None, None, :] * dz
            lb = lb + h * db / b
            lb = lb - np.log(np.exp(lb).sum(axis=1, keepdims=True) / 3)
            dphi = (cphi * np.cosh(kk * Th[:, :, 0]) - phi0) / b; bx, co = g * boxt * h, g * (nonfit - boxt) * h
            acc[:, 0] += dphi; acc[:, 1] += bx; acc[:, 2] += co; acc[:, 3] += u * h - dphi - bx - co; acc[:, 4] += uF * h
    print(f"[{time.time() - t0:.0f}s] K half", flush=True)
    tr = lambda x: np.array(x).swapaxes(0, 1)
    out = dict(seeds=tuple(seeds), S=S, s=np.array(s_), var={})
    for v in variants:
        ix = [i for i, (x, _) in enumerate(runs) if x == v]
        out['var'][v] = dict(al=tr(al)[ix], W=tr(Wr)[ix], acc=tr(ac)[ix])
    return out


def k_identity(r):
    """The decomposition's residual: Δlog(β_k/β_3) − γΣ(the terms), all runs and records."""
    x = np.log(r['W']); x = x - x[:, :, 2:3]
    t = GAMMA_E * r['acc'].sum(axis=2); t = t - t[:, :, 2:3]
    return x - t


def k_offset(al, W):
    """The positions' offset from the fit as a mass offset, (γ/m)[(Υ_k − Υ_k^fit) − W_kΣ_l(Υ_l − Υ_l^fit)] (Proposition 7.3(a)'s bounded
    term about the exact fit at the masses W, Υ^fit_k = Υ(α_k) at |a_k| = |A_k|/W_k): W less the masses of the fit that the
    noise-free flow reaches, at first order."""
    c = A_MAX ** 2 / (4 * E.G[0]); W = np.asarray(W, float)
    afit = A_MAX * np.arctanh(np.abs(E.A_BUMP) / (W * A_MAX))
    d = c * (np.cosh(2 * np.asarray(al, float) / A_MAX) - np.cosh(2 * afit / A_MAX))
    return GAMMA_E / 3 * (d - W * d.sum(axis=-1, keepdims=True))


def k_response(variant, s, du, ds=100.0):
    """The masses' linear response to an added fitness du (seeds × (records − 1) × 3, constant between the records s) along the slow
    limit's leading-order path from the masses 1/3 (lead_fitness's replicator flow; no wall is reached): δW′ = J(W)δW + γW∘(du − W·du),
    J the flow's Jacobian in (W_1, W_2) by central differences, once per record interval; Euler at the step ds. Returns δW at the records (seeds × records × 3)."""
    co = lead_coeffs(variant)
    f = lambda W: GAMMA_E * GN.replicator_flow(W, lead_fitness(variant, W, co))
    W = np.ones(3) / 3; dW = np.zeros((du.shape[0], 3)); out = [dW.copy()]; eps = 1e-6
    for j in range(len(s) - 1):
        Jc = np.array([(f(W + eps * e) - f(W - eps * e))[:2] / (2 * eps) for e in (np.array([1.0, 0, -1]), np.array([0, 1.0, -1]))]).T
        d = du[:, j]
        for _ in range(int(round((s[j + 1] - s[j]) / ds))):
            dW2 = dW[:, :2] + ds * (dW[:, :2] @ Jc.T + GAMMA_E * (W * (d - (d @ W)[:, None]))[:, :2])
            dW = np.c_[dW2, -dW2.sum(axis=1)]
            W = W + ds * f(W)
        out.append(dW.copy())
    return np.array(out).swapaxes(0, 1)


def k_forcing(K, v):
    """The O(Δt) part of the identity's terms (boundary, box, cohesion, chain: all but the masses form's, whose change is the masses'
    response), from the paired runs at Δt and Δt/2: 2(acc(Δt) − acc(Δt/2)) per class, fitted through the origin in s, per seed — the
    fitness the step adds, per unit time."""
    hv, fv = K['half']['var'][v], K['var'][v]; s = K['half']['s']; n = len(s)
    f = 2 * (fv['acc'][:, :n, :4] - hv['acc'][:, :, :4]).sum(axis=2)
    return (f * s[None, :, None]).sum(axis=1) / (s ** 2).sum()


def k_effects(K, E_rec, F, v):
    """Per seed, the second half's means (Part E's half): W − Katzenberger's SDE mean (Part F); the positions' offset (k_offset at the
    records); the O(Δt) effect, the masses' linear response (k_response) to k_forcing held constant over the run; the SDE's std of a path's
    second-half mean; and the paired runs' 2(W(Δt) − W(Δt/2)) at the half run's end with the response's value there."""
    r = K['var'][v]; s = K['s']; R = np.asarray(E_rec[v]); hi = R.shape[1] // 2
    du = np.repeat(k_forcing(K, v)[:, None, :], len(s) - 1, axis=1)
    eff = k_response(v, s, du)
    sh = F['var'][v]['sde']['half']; n = len(K['half']['s'])
    return dict(obs=R[:, hi:, 4:7].mean(axis=1) - sh.mean(axis=0), off=k_offset(r['al'], r['W'])[:, hi:].mean(axis=1),
                eff=eff[:, hi:].mean(axis=1), sde_std=sh.std(axis=0, ddof=1), at_half=2 * (r['W'][:, n - 1] - K['half']['var'][v]['W'][:, -1]),
                resp_half=eff[:, n - 1], kpath=F['var'][v]['kpath']['half'][0])


def part_K(E_rec, old=None):
    """k_runs (reused from a stored K of the same arguments: it is deterministic in them), k_half, and the replay's agreement with
    Part E's records."""
    same = old is not None and old.get('S') == S_E and tuple(old.get('seeds', ())) == tuple(SEEDS_E) and set(old['var']) == set(VARIANTS_K)
    K = {k: v for k, v in old.items() if k in ('seeds', 'S', 's', 'var')} if same else k_runs((VARIANTS_K, SEEDS_E, S_E))
    K['half'] = k_half((VARIANTS_K, SEEDS_E, S_KH))
    K['replay'] = float(max(max(np.abs(K['var'][v]['W'] - np.asarray(E_rec[v])[:, :, 4:7]).max(),
                                np.abs(K['var'][v]['al'] - np.asarray(E_rec[v])[:, :, 1:4]).max()) for v in VARIANTS_K))
    return K


# ---- the missing O(γ²τ) drift (D‴.40 after round 153): Part M ----
VARIANT_M = 'nonlinear only'
LIM_STEPS = ((1.0, 4000), (8.0, 4500))       # the limit map: RK4 at the step 1 to s = 4000 (the fast modes), then 8 to s = 4·10⁴
S_CHK_M = 3e4                                # the limit map's convergence: its masses at 3·10⁴ against 4·10⁴
EVERY_M, CHUNK_MK = 5, 4                     # Part K's second-half records: every fifth; seeds per chunk
GAMMAS_M, SEEDS_M, S_M, T_YM, CHUNK_M = (1e-2, 3e-3), tuple(range(800, 848)), 2e5, (2e4, 1e5, 2e5), 16
TAUS_MF, SEEDS_MF, BURN_MF, CHUNK_MF = (3e-5, 1e-4, 3e-4), tuple(range(900, 932)), 2e4, 16
S_MF = {3e-5: (1e5, 16), 1e-4: (2e5, 32), 3e-4: (1e5, 16)}      # slow units after the burn-in and seeds, per level τ
GAMMAS_FD_M, H_FD_M = (1e-2, 1e-3), (0.05, 0.025)


def m_limit(Th, lb, gamma, steps=LIM_STEPS, s_check=S_CHK_M):
    """The limit map's masses: the exact fit that the noise-free 'nonlinear only' loop (the box, the rate γ) reaches from the states
    (Θ (B, 3, 3), log β (B, 3)), by RK4 on (Θ, log(β_1/β_3), log(β_2/β_3)) at the steps `steps` ((step, count) pairs: the step 1 while the
    fast modes decay — RK4 at the step 4 misses the masses' quadratic response by 10% —, then 8); returns W at the end and at s_check."""
    from src.core import gauge_batch as GB
    Th = np.array(Th, float); B = len(Th)
    bd = GB.BatchDither(e_loop(VARIANT_M)[0], gamma=float(gamma))
    lr = np.asarray(lb, float)[:, :2] - np.asarray(lb, float)[:, 2:3]

    def f(T, l):
        e = np.exp(np.c_[l, np.zeros(B)]); b = 3 * e / e.sum(axis=1, keepdims=True)
        dT, db = bd.drift(T, b); dl = db / b
        return dT, dl[:, :2] - dl[:, 2:3]
    w_of = lambda l: np.exp(np.c_[l, np.zeros(B)]) / np.exp(np.c_[l, np.zeros(B)]).sum(axis=1, keepdims=True)
    s, chk = 0.0, None
    for ds, n in steps:
        for _ in range(n):
            k1 = f(Th, lr); k2 = f(Th + ds / 2 * k1[0], lr + ds / 2 * k1[1])
            k3 = f(Th + ds / 2 * k2[0], lr + ds / 2 * k2[1]); k4 = f(Th + ds * k3[0], lr + ds * k3[1])
            Th = Th + ds / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); lr = lr + ds / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            s += ds
            if chk is None and s >= s_check - 1e-9:
                chk = w_of(lr)
    return w_of(lr), chk


def m_offset_lin(al, W, gamma=GAMMA_E):
    """Part K's offset (1) with the fit's motion: W − Y = (I − M)⁻¹(1) at first order, M = (γ/m)(I − W1ᵀ)diag(dΥ_k^fit/dW_k), dΥ_k^fit/dW_k
    = −h_k²/(G_αW_k), h_k = a_k/a′_k at the fit (|a_k| = |A_k|/W_k): the masses' change moves the fit, whose amplitudes' change moves the
    masses again (the gauge identity along the noise-free flow)."""
    W = np.asarray(W, float); al = np.asarray(al, float)
    off = k_offset(al, W) * gamma / GAMMA_E
    a = np.abs(E.A_BUMP) / W; h = a / (1 - (a / A_MAX) ** 2); D = -h ** 2 / (E.G[0] * W)
    Mm = (gamma / 3) * (np.eye(3) - W[..., :, None]) * D[..., None, :]
    return np.linalg.solve(np.eye(3) - Mm, off[..., None])[..., 0], Mm


def m_lim_K(K, c0):
    """The exact limit map at Part K's 'nonlinear only' records (the second half, every EVERY_M-th record) for the seeds c0 … c0 + CHUNK_MK − 1."""
    r = K['var'][VARIANT_M]; n = r['W'].shape[1]; idx = np.arange(n // 2, n, EVERY_M); sl = slice(c0, c0 + CHUNK_MK)
    Th = r['Th'][sl][:, idx]; W = r['W'][sl][:, idx]; sh = W.shape
    Y, chk = m_limit(Th.reshape(-1, 3, 3), np.log(3 * W.reshape(-1, 3)), GAMMA_E)
    return dict(c0=c0, idx=idx, Y=Y.reshape(sh), chk=chk.reshape(sh))


@cached
def m_loop(args):
    """Part E's 'nonlinear only' loop (the box, τ = TAU) at the rate γ for the seeds, from the exact fit at the masses 1/3, as two Euler
    chains on one Brownian path (k_half's coupling: the step Δt with run_dithered's draws sq·z, the step Δt/2 with sq·(z + w)/2 and
    sq·(z − w)/2, w from the generator seed + W_OFF_KH); rows: the Δt chain's, then the Δt/2 chain's. W and α every REC full steps; the
    states at T_YM and their limit maps (m_limit)."""
    gamma, seeds, S = args
    t0 = time.time()
    from src.core import gauge_batch as GB
    loop, kw, Th0 = e_loop(VARIANT_M)
    B = len(seeds); h = DT / 2; F_, H_ = slice(0, B), slice(B, 2 * B)
    bd = GB.BatchDither(loop, gamma=float(gamma))
    sq = np.sqrt(2.0 * TAU * DT * np.asarray(loop.G, float)) * np.asarray(kw['mask'], float)
    gens = [np.random.default_rng(s) for s in seeds]; gw = [np.random.default_rng(s + W_OFF_KH) for s in seeds]
    st = np.r_[np.full(B, DT), np.full(B, h)][:, None]
    Th = np.tile(Th0, (2 * B, 1, 1)); lb = np.zeros((2 * B, 3))
    n = int(round(S / DT)); s_, Wr, ar, states = [], [], [], {}
    renorm = lambda x: x - np.log(np.exp(x).sum(axis=1, keepdims=True) / 3)
    for k in range(n + 1):
        if k % REC == 0:
            s_.append(k * DT); Wr.append(np.exp(lb) / 3); ar.append(Th[:, :, 0].copy())
        for t in T_YM:
            if abs(k * DT - t) < DT / 2:
                states[t] = (Th.copy(), lb.copy())
        if k == n:
            break
        if k % REC == 0:
            Z = np.array([g_.standard_normal((REC, 3, 3)) for g_ in gens]); Zw = np.array([g_.standard_normal((REC, 3, 3)) for g_ in gw])
        z, w = Z[:, k % REC], Zw[:, k % REC]
        b = np.exp(lb); dT, db = bd.drift(Th, b)
        Th = Th + st[:, :, None] * dT + np.concatenate([sq * z, 0.5 * sq * (z + w)])
        lb = renorm(lb + st * db / b)
        b = np.exp(lb[H_]); dT, db = bd.drift(Th[H_], b)
        Th[H_] = Th[H_] + h * dT + 0.5 * sq * (z - w)
        lb[H_] = renorm(lb[H_] + h * db / b)
    ts = [t for t in T_YM if t in states]
    Y, chk = m_limit(np.concatenate([states[t][0] for t in ts]), np.concatenate([states[t][1] for t in ts]), gamma)
    print(f"[{time.time() - t0:.0f}s] M loop gamma {gamma:g} seeds {seeds[0]}-{seeds[-1]}", flush=True)
    tr = lambda x: np.array(x).swapaxes(0, 1)
    return dict(gamma=gamma, seeds=tuple(seeds), S=S, s=np.array(s_), W=tr(Wr), al=tr(ar), t=np.array(ts),
                Y=Y.reshape(len(ts), 2 * B, 3), chk=chk.reshape(len(ts), 2 * B, 3))


@cached
def m_fixed(args):
    """The positions at the fixed masses 1/3 (γ = 0, the box), the centres and widths dithered at the level τ, from the exact fit, as the two
    coupled Euler chains of m_loop; after BURN_MF slow units, per seed and chain the sums over the states (each half step's for the Δt/2
    chain) of φ_l, φ̇_k = φ′_k(α_k)b_α_k (the amplitudes move by their drift alone) and φ_lφ̇_k, and their count."""
    tau, seeds, S, burn = args
    t0 = time.time()
    from src.core import gauge_batch as GB
    loop, kw, Th0 = e_loop(VARIANT_M)
    B = len(seeds); h = DT / 2; H_ = slice(B, 2 * B)
    bd = GB.BatchDither(loop, gamma=0.0)
    sq = np.sqrt(2.0 * tau * DT * np.asarray(loop.G, float)) * np.asarray(kw['mask'], float)
    gens = [np.random.default_rng(s) for s in seeds]; gw = [np.random.default_rng(s + W_OFF_KH) for s in seeds]
    st = np.r_[np.full(B, DT), np.full(B, h)][:, None, None]
    kk, cphi = 2.0 / A_MAX, A_MAX ** 2 / (4 * E.G[0]); beta = np.ones((2 * B, 3))
    Th = np.tile(Th0, (2 * B, 1, 1))
    acc = dict(n=np.zeros(2 * B), p=np.zeros((2 * B, 3)), pd=np.zeros((2 * B, 3)), ppd=np.zeros((2 * B, 3, 3)))

    def add(T, dT, rows):
        p = cphi * np.cosh(kk * T[:, :, 0]); pd = cphi * kk * np.sinh(kk * T[:, :, 0]) * dT[:, :, 0]
        acc['n'][rows] += 1; acc['p'][rows] += p; acc['pd'][rows] += pd; acc['ppd'][rows] += p[:, :, None] * pd[:, None, :]
    nb, n = int(round(burn / DT)), int(round(S / DT))
    for k in range(nb + n):
        if k % REC == 0:
            Z = np.array([g_.standard_normal((REC, 3, 3)) for g_ in gens]); Zw = np.array([g_.standard_normal((REC, 3, 3)) for g_ in gw])
        z, w = Z[:, k % REC], Zw[:, k % REC]
        dT = bd.drift(Th, beta)[0]
        if k >= nb:
            add(Th, dT, slice(None))
        Th = Th + st * dT + np.concatenate([sq * z, 0.5 * sq * (z + w)])
        dT = bd.drift(Th[H_], beta[H_])[0]
        if k >= nb:
            add(Th[H_], dT, H_)
        Th[H_] = Th[H_] + h * dT + 0.5 * sq * (z - w)
    print(f"[{time.time() - t0:.0f}s] M fixed tau {tau:g} seeds {seeds[0]}-{seeds[-1]}", flush=True)
    return dict(tau=tau, seeds=tuple(seeds), S=S, burn=burn, **{c: {k_: v[sl] for k_, v in acc.items()} for c, sl in (('full', slice(0, B)), ('half', H_))})


def m_anti(a, tau):
    """Per seed, the positions' correlation E[φ_lφ̇_k] − E[φ_l]E[φ̇_k] (m_fixed's sums) → the drift it gives the masses per γ²τ,
    (1/m²)Σ_l of its antisymmetric part (Proposition D‴.40(ii): (γ/m)²Σ_lE[φ_l dφ_k/ds]), and the same of its symmetric part (zero for the
    diffusion, whose E[φ_lφ_k] is stationary: the chain's bias)."""
    n = a['n'][:, None]; c = a['ppd'] / n[:, :, None] - (a['p'] / n)[:, :, None] * (a['pd'] / n)[:, None, :]
    return 0.5 * (c - c.swapaxes(1, 2)).sum(axis=1) / (9 * tau), 0.5 * (c + c.swapaxes(1, 2)).sum(axis=1) / (9 * tau)


def m_ou(tau=TAU, dts=(DT, DT / 2), hJ=1e-5):
    """The linearised positions at the fixed masses 1/3 (Ornstein–Uhlenbeck, drift matrix B = the loop's Jacobian at the fit, the centres'
    and widths' noise Σ): C from BC + CBᵀ + Σ = 0, E[δα_l dδα_k/ds] = (CBᵀ)_{α_lα_k} (antisymmetric: α carries no noise), the drift per γ²τ
    (1/m²)Σ_lφ′_lφ′_k(CBᵀ)_{α_lα_k}/τ (Katzenberger's at γ → 0); the Euler chains' at the steps dts (C_Δt = (I + ΔtB)C_Δt(I + ΔtB)ᵀ +
    ΔtΣ, the states' law), the antisymmetric and the symmetric part of φ′_lφ′_k(C_ΔtBᵀ)_{α_lα_k}; the antisymmetric correlation's share, E[δα_lδα̇_k]/√(E[δα_l²]E[δα̇_k²]), largest over l ≠ k; B's rates."""
    from scipy.linalg import solve_continuous_lyapunov, solve_discrete_lyapunov
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=0.0)
    th0 = E.exact_rep(net, np.ones(3) / 3).ravel(); f = lambda t: loop.drift(t.reshape(3, 3), np.ones(3))[0].ravel()
    Bm = np.array([(f(th0 + hJ * e) - f(th0 - hJ * e)) / (2 * hJ) for e in np.eye(9)]).T
    Sig = noises(tau)['centres and widths only']
    dphi = (A_MAX / (2 * E.G[0])) * np.sinh(2 * th0[0::3] / A_MAX)
    def drift(C, sym=False):
        X = (dphi[:, None] * dphi[None, :]) * (C @ Bm.T)[0::3, 0::3]
        return 0.5 * (X + (1 if sym else -1) * X.T).sum(axis=0) / (9 * tau)
    C = solve_continuous_lyapunov(Bm, -Sig)
    CB = (C @ Bm.T)[0::3, 0::3]; V = (Bm @ C @ Bm.T)[0::3, 0::3]
    rho = CB / np.sqrt(np.diag(C)[0::3][:, None] * np.diag(V)[None, :])
    Cd = [solve_discrete_lyapunov(np.eye(9) + dt * Bm, dt * Sig) for dt in dts]
    return dict(tau=tau, d=drift(C), chain=np.array([drift(c) for c in Cd]), chain_sym=np.array([drift(c, True) for c in Cd]), dts=dts,
                rho=np.abs(rho - np.diag(np.diag(rho))).max(), sd_alpha=np.sqrt(np.diag(C)[0::3]),
                rates=np.sort(-np.linalg.eigvals(Bm).real), dphi=dphi)


@cached
def m_fd(gamma):
    """Katzenberger's drift ½∂²Φ[Σ] of the masses, independently: the limit map Φ (m_limit) at the exact fit at the masses 1/3 displaced by
    ±h along each dithered coordinate (the centres' and widths'), its second differences at h ∈ H_FD_M, Richardson's extrapolation; with
    `limit_drift`'s (Part C's) at the same γ. Per γ²τ."""
    net, _ = E.setup()
    Th0 = E.exact_rep(net, np.ones(3) / 3).ravel(); Sig = noises()['centres and widths only']; idx = [i for i in range(9) if Sig[i, i] > 0]
    xs = [Th0]
    for hh in H_FD_M:
        for i in idx:
            for sg in (1.0, -1.0):
                x = Th0.copy(); x[i] += sg * hh; xs.append(x)
    Y, chk = m_limit(np.array(xs).reshape(-1, 3, 3), np.zeros((len(xs), 3)), gamma)
    Yp = Y[1:].reshape(len(H_FD_M), len(idx), 2, 3)
    d2 = (Yp[:, :, 0] + Yp[:, :, 1] - 2 * Y[0]) / np.asarray(H_FD_M)[:, None, None] ** 2
    dr = 0.5 * np.einsum('i,hik->hk', np.diag(Sig)[idx], d2) / (gamma ** 2 * TAU)
    b = k_drift(katz(VARIANT_M, gamma), np.ones(3) / 3)[0]
    return dict(gamma=gamma, h=H_FD_M, at_fit=Y[0] - 1 / 3, d2=d2, fd=dr, rich=(4 * dr[1] - dr[0]) / 3,
                kz=np.r_[b, -b.sum()] / (gamma ** 2 * TAU), chk=float(np.abs(Y - chk).max()))


def part_M(res, save):
    """Part M, resumable: each piece is put into res['M'] and saved (`save`) as it completes, and is not recomputed when present."""
    M = dict(res.get('M') or {})

    def put(k, v):
        M[k] = v; res['M'] = M; save(); print(f"M {k}", flush=True)
    if 'ou' not in M:
        put('ou', [m_ou(t) for t in TAUS_MF])
    if 'fd' not in M:
        put('fd', [m_fd(g) for g in GAMMAS_FD_M])
    if 'katz' not in M:
        put('katz', {g: np.r_[b, -b.sum()] / (g ** 2 * TAU) for g in GAMMAS_M for b in [k_drift(katz(VARIANT_M, g), np.ones(3) / 3)[0]]})
    lk = dict(M.get('lim_K') or {})
    for c0 in range(0, len(SEEDS_E), CHUNK_MK):
        if c0 not in lk:
            lk[c0] = m_lim_K(res['K'], c0); put('lim_K', dict(lk))
    lp = {g: dict(v) for g, v in (M.get('loop') or {}).items()}
    for g in GAMMAS_M:
        for c0 in range(0, len(SEEDS_M), CHUNK_M):
            if c0 not in lp.setdefault(g, {}):
                lp[g][c0] = m_loop((g, SEEDS_M[c0:c0 + CHUNK_M], S_M)); put('loop', {k: dict(v) for k, v in lp.items()})
    fx = {t: dict(v) for t, v in (M.get('fixed') or {}).items()}
    for t in TAUS_MF:
        S, ns = S_MF[t]
        for c0 in range(0, ns, CHUNK_MF):
            if c0 not in fx.setdefault(t, {}):
                fx[t][c0] = m_fixed((t, SEEDS_MF[c0:min(c0 + CHUNK_MF, ns)], S, BURN_MF)); put('fixed', {k: dict(v) for k, v in fx.items()})
    return M


# ---- the finite-τ question (D‴.40 after Part M): Part N ----
SEEDS_N = {1e-4: tuple(range(1100, 1132)), 1e-3: tuple(range(1300, 1316)), 3e-4: tuple(range(1400, 1544))}   # in the order computed
S_N = {1e-4: 2e5, 1e-3: 1e5, 3e-4: 1e5}           # slow units after BURN_MF: Part M(d)'s at 10⁻⁴ and 3·10⁻⁴, so that the seeds pool
CHUNK_N = 16
TAUS_FIT_N = (3e-5, 1e-4, 3e-4)                    # the line's points (10⁻³: 16 seeds, against the line's prediction only)


def part_N(res, save):
    """Part N, resumable: more seeds of m_fixed (Part M(d)'s runs, the same burn-in, step and coupling; new seeds) at the levels SEEDS_N,
    each chunk of CHUNK_N seeds put into res['N']['fixed'][τ][c0] and saved (`save`) as it completes, and not recomputed when present."""
    Nn = dict(res.get('N') or {}); fx = {t: dict(v) for t, v in (Nn.get('fixed') or {}).items()}
    for t, sd in SEEDS_N.items():
        for c0 in range(0, len(sd), CHUNK_N):
            if c0 not in fx.setdefault(t, {}):
                fx[t][c0] = m_fixed((t, sd[c0:c0 + CHUNK_N], S_N[t], BURN_MF))
                Nn['fixed'] = {k: dict(v) for k, v in fx.items()}; res['N'] = Nn; save(); print(f"N fixed {t:g} {c0}", flush=True)
    return Nn


def n_pool(M, Nn, t):
    """m_fixed_all over Part M(d)'s seeds and Part N's at the level τ (either may be absent), the chunks in that order; None if neither."""
    ch = {(0, c): v for c, v in ((M.get('fixed') or {}).get(t) or {}).items()}
    ch.update({(1, c): v for c, v in ((Nn.get('fixed') or {}).get(t) or {}).items()})
    return m_fixed_all({'fixed': {t: ch}}, t) if ch else None


def n_fit(tau, y, e):
    """Per class, the weighted least-squares line y = c + bτ through the points (τ_i, y_i ± e_i): the intercept c (the limit τ → 0) and the
    slope b with their standard errors and χ² (len(τ) − 2 degrees of freedom); and the constant (the weighted mean, its s.e., χ² with
    len(τ) − 1)."""
    x = np.asarray(tau, float)[:, None]; y = np.asarray(y, float); w = 1 / np.asarray(e, float) ** 2
    S, Sx, Sxx, Sy, Sxy = w.sum(axis=0), (w * x).sum(axis=0), (w * x * x).sum(axis=0), (w * y).sum(axis=0), (w * x * y).sum(axis=0)
    D = S * Sxx - Sx ** 2; c = (Sxx * Sy - Sx * Sxy) / D; b = (S * Sxy - Sx * Sy) / D; m = Sy / S
    return dict(c=c, ec=np.sqrt(Sxx / D), b=b, eb=np.sqrt(S / D), cov=-Sx / D, chi=(w * (y - c - b * x) ** 2).sum(axis=0),
                m=m, em=1 / np.sqrt(S), chim=(w * (y - m) ** 2).sum(axis=0), dof=len(x) - 2)


def report_J(J):
    se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M))
    f3 = lambda v, p='+.2e': ', '.join(f'{x:{p}}' for x in v)
    fu, ha = J['full'], J['half']; dt, h = J['dt']
    L_ = [f"## J. The identity's own check: the chain's per-state bias at the steps Δt and Δt/2 (referee round 153)", "",
          f"Part D's runs at the fixed masses 1/3 with the box ({J['S']:g} slow units, seed {J['seed']}), each as two Euler–Maruyama chains "
          f"driven by one Brownian path, at the step Δt = {dt:g} and Δt/2 = {h:g}. At every state the chain's exact per-step bias "
          "b_k(h) = (E_n[ψ_k(α_(n+1))] − ψ_k(α_n))/h − Lψ_k (closed form) is computed at both steps. If Lψ_k = u_k − g_k∂_{α_k}(J_coh + V) + "
          "τ_kg′_k (the identity) holds, b_k(h) = O(h) at every state and halves with the step; were a term of the identity's right side wrong, "
          "b_k would keep that term's O(1) error at every step — as it does below for the right side without its Itô term τ_kg′_k. "
          f"(Part D's χ² test of E[Lψ_k] + E[b_k] holds by construction: it checks the averaging, not the identity.) Errors from "
          f"{N_BATCH_J} batch means; the ratios' errors by the delta method on the paired batches.", "",
          "| noise | τ | E[b_k] at Δt (the Δt chain) | E[b_k] at Δt/2 (the Δt/2 chain) | their ratio (± s.e.) | at the Δt/2 chain's states: "
          "E[b_k(Δt)]/E[b_k(Δt/2)] (± s.e.) | without the Itô term: E[b_k] at Δt, Δt/2 and their ratio | E[b_k]/τ_kE[g′_k] at Δt, Δt/2 | "
          "E[u_k] at Δt − at Δt/2 (± s.e.) |", "|---|---|---|---|---|---|---|---|---|"]
    rat_c, rat_s, rat_w, rel_f, rel_h = [], [], [], [], []
    for j, (nm, tau) in enumerate(J['runs']):
        bf, bh, ba = fu['b'][:, j], ha['b'][:, j], ha['b_alt'][:, j]
        r1, e1 = ratio_se(bf, bh); r2, e2 = ratio_se(ba, bh); rat_c += list(r1); rat_s += list(r2)
        if np.all(J['tk'][j] > 0):
            wf, wh = bf - fu['gp'][:, j], bh - ha['gp'][:, j]; rw = wf.mean(axis=0) / wh.mean(axis=0)
            wrong = '; '.join(f'{a:+.2e}, {b:+.2e} ({c:.3f})' for a, b, c in zip(wf.mean(axis=0), wh.mean(axis=0), rw))
            rf, rh = bf.mean(axis=0) / -fu['gp'][:, j].mean(axis=0), bh.mean(axis=0) / -ha['gp'][:, j].mean(axis=0)
            rel = f3(rf, '+.4f') + '; ' + f3(rh, '+.4f'); rat_w += list(rw); rel_f += list(rf); rel_h += list(rh)
        else:
            wrong, rel = '— (τ_k = 0)', '— (τ_k = 0)'
        du = fu['u'][:, j] - ha['u'][:, j]
        L_.append(f"| {nm} | {tau:g} | {f3(bf.mean(axis=0))} | {f3(bh.mean(axis=0))} | {', '.join(f'{a:.3f} ± {b:.3f}' for a, b in zip(r1, e1))} | "
                  f"{', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(r2, e2))} | {wrong} | {rel} | "
                  f"{', '.join(f'{a:+.1e} ± {b:.1e}' for a, b in zip(du.mean(axis=0), se(du)))} |")
    rat_c, rat_s = np.array(rat_c), np.array(rat_s)
    zf = lambda c: ((c['lpsi'].mean(axis=0) + c['b'].mean(axis=0)) / np.where(c['se_mart'] > 0, c['se_mart'], np.nan)).ravel()
    L_ += ["", f"The two chains' bias ratios: {rat_c.min():.3f} to {rat_c.max():.3f} between the chains, {rat_s.min():.4f} to "
           f"{rat_s.max():.4f} at the same states (2 for a bias exactly linear in the step); without the Itô term {min(rat_w):.3f} to "
           f"{max(rat_w):.3f}. The bias over τ_kE[g′_k], runs with amplitude noise: {min(rel_f):.4f} to {max(rel_f):.4f} at Δt, "
           f"{min(rel_h):.4f} to {max(rel_h):.4f} at Δt/2. The share of states with α_k outside the box, "
           f"largest over the runs: {fu['out_box'].max():.3f} (Δt), {ha['out_box'].max():.3f} (Δt/2). The bookkeeping check of Part D in each "
           f"chain, (E[Lψ_k] + E[b_k]) over its martingale error, runs with amplitude noise: largest |z| {np.nanmax(np.abs(zf(fu))):.2f} (Δt), "
           f"{np.nanmax(np.abs(zf(ha))):.2f} (Δt/2).", ""]
    return L_


def report_K(K, E_rec, F):
    pe = lambda M: ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(*_ms(M)))
    zz = lambda M, sd: (np.abs(M.mean(axis=0)) / zden(_ms(M)[1], sd, len(M))).max()
    H_ = K['half']
    L_ = [f"## K. The Euler chain's bias along Part E's deterministic-amplitude runs (referee round 153)", "",
          f"Part E's variants without amplitude noise replayed over {K['S']:g} slow units with the gauge identity's decomposition (`k_runs`, "
          f"gauge_batch's drift, each seed run_dithered's draws: Part E's records reproduced to {K['replay']:.1e}; the decomposition's largest "
          f"residual {max(np.abs(k_identity(r)).max() for r in K['var'].values()):.1e}), and over {H_['S']:g} slow units at the step Δt/2 = "
          f"{DT / 2:g} (`k_half`: each full step's increment of Part E's run split by the Brownian bridge, so that both chains follow one "
          "Brownian path). With the amplitudes noise-free the gauge identity reads u_kΔt = Δφ_k/β_k + (g_k∂_{α_k}V + g_k∂_{α_k}J_coh)Δt + the "
          "chain's term, −(φ_k(α_n + Δt b_α) − φ_k(α_n) − Δtφ′_kb_α)/β_k, O(Δt²) per step. Against Katzenberger's SDE (Part F), which describes "
          "the masses of the fit that the noise-free flow reaches, the measured masses carry (1) the positions' offset from that fit, "
          "(γ/m)[(Υ_k − Υ_k^fit) − W_kΣ_l(Υ_l − Υ_l^fit)] at first order (Proposition D‴.40's bounded term about the fit at the current "
          "masses), and (2) the step's O(Δt) effect: the masses' linear response (along the slow limit's path) to the fitness the step adds, "
          "2(the identity's terms at Δt − at Δt/2) per unit time, fitted through the origin over the paired runs and held constant. Where "
          "the masses form is the only restoring drift ('nonlinear only + barrier') it is evaluated at the measured masses, so that its zero "
          "holds them, not the fit's masses: (1) shifts the fit's masses by its opposite and cancels from the measured ones. Second-half "
          "means (Part E's half), per seed; z as in Part F.", "",
          "| variant | measured − the SDE's mean (± s.e.) | largest \\|z\\| | (1) the positions' offset (± s.e.) | (2) the O(Δt) effect (± s.e.) | "
          "measured − SDE − (2) (± s.e.) | largest \\|z\\| | measured − SDE − (1) − (2) (± s.e.) | largest \\|z\\| |",
          "|---|---|---|---|---|---|---|---|---|"]
    X = {v: k_effects(K, E_rec, F, v) for v in VARIANTS_K}
    for v, x in X.items():
        a, b = x['obs'] - x['eff'], x['obs'] - x['off'] - x['eff']
        L_.append(f"| {v} | {pe(x['obs'])} | {zz(x['obs'], x['sde_std']):.1f} | {pe(x['off'])} | {pe(x['eff'])} | {pe(a)} | "
                  f"{zz(a, x['sde_std']):.1f} | {pe(b)} | {zz(b, x['sde_std']):.1f} |")
    x = X[VARIANTS_K[0]]
    L_ += ["", f"'{VARIANTS_K[0]}' has no restoring drift, so (1) and (2) add to its masses; what is left, measured − 1/3 − (1) − (2) = "
           f"{pe(x['obs'] + F['var'][VARIANTS_K[0]]['sde']['half'].mean(axis=0) - 1 / 3 - x['off'] - x['eff'])}, is the loop's own secular "
           f"displacement, against Katzenberger's drift path's {', '.join(f'{c:+.2e}' for c in x['kpath'] - 1 / 3)} (its O(γ²τ) drift).", "",
           f"The paired runs at Δt and Δt/2: W(Δt) − W(Δt/2), half the step's O(Δt) effect, with its parts from the chain's term and from the "
           "boundary term (their log-mass differences mapped by ∂W/∂log β at the Δt run's masses), and at the end the response model's value "
           "of 2(W(Δt) − W(Δt/2)) against the runs'.", "",
           "| variant | s | W(Δt) − W(Δt/2) (± s.e.) | the chain's term's part | the boundary term's part |", "|---|---|---|---|---|"]
    jac = lambda W: W[..., :, None] * (np.eye(3) - W[..., None, :])
    for v in VARIANTS_K:
        hv, fv = H_['var'][v], K['var'][v]; s = H_['s']; n = len(s); Wf = fv['W'][:, :n]
        lin = lambda A: np.einsum('rtkl,rtl->rtk', jac(Wf), GAMMA_E * (A - A[..., 2:3]))
        dacc = fv['acc'][:, :n] - hv['acc']
        for t in (5e4, 1e5, 2e5):
            i = int(np.argmin(np.abs(s - t)))
            L_.append(f"| {v} | {s[i]:g} | {pe(Wf[:, i] - hv['W'][:, i])} | {pe(lin(dacc[:, :, 3])[:, i])} | {pe(lin(dacc[:, :, 0])[:, i])} |")
    L_ += ["", "At s = " + f"{H_['s'][-1]:g}, 2(W(Δt) − W(Δt/2)) measured and from the response model: "
           + '; '.join(f"'{v}' {', '.join(f'{c:+.2e}' for c in X[v]['at_half'].mean(axis=0))} and "
                       f"{', '.join(f'{c:+.2e}' for c in X[v]['resp_half'].mean(axis=0))}" for v in VARIANTS_K) + ".", ""]
    v = VARIANTS_K[0]; r = K['var'][v]; hi = r['W'].shape[1] // 2
    c = GAMMA_E * (r['acc'] - r['acc'][:, :, :, 2:3]); lin = np.einsum('kl,rtjl->rtjk', JAC13, c)       # each term's W, linearised at 1/3
    off = k_offset(r['al'], r['W']); off0 = off - off[:, :1]
    half = lambda M: M[:, hi:].mean(axis=1)
    rows = [("W − 1/3 (measured)", half(r['W'] - 1 / 3)),
            ("the boundary term ΣΔφ_k/β_k", half(lin[:, :, 0])),
            ("of which the positions' offset (1), from its start", half(off0)),
            ("of which the rest (secular)", half(lin[:, :, 0] - off0)),
            ("the box's term", half(lin[:, :, 1])), ("the cohesion's", half(lin[:, :, 2])), ("the chain's term", half(lin[:, :, 3])),
            ("the terms' sum − (W − 1/3): the linearisation", half(lin[:, :, :4].sum(axis=2) - (r['W'] - 1 / 3)))]
    L_ += [f"'{v}' at the step Δt: W − 1/3 split by the identity (each term's contribution linearised at the masses 1/3, ∂W/∂log β = "
           "W_k(δ_kl − W_l); second-half means, per seed); the step's effect (2) is the chain's term with the boundary term's own O(Δt) part.",
           "", "| term | W_1, W_2, W_3 (± s.e.) |", "|---|---|"]
    L_ += [f"| {lab} | {pe(M)} |" for lab, M in rows]
    L_.append("")
    return L_


def report_L(F, H, I, C):
    """Errors and extremes that the paper quotes (referee round 153): Part H's ratios with their delta-method errors, Part F's zero shift per
    class, Part C's and H's largest values, Part I's regressed s.e. with b's uncertainty."""
    L_ = ["## L. Errors and extremes of Parts C, F, H and I (referee round 153)", ""]
    mc = {r['tc']: r for r in H['mc']}; w, c = mc[0.0], mc[TC_H]
    xs, xw = (c['u'] - c['box']), (w['u'] - w['box']); n = len(xs)
    d = xs.mean(axis=0) / xw.mean(axis=0)
    ed = np.abs(d) * np.sqrt((xs.std(axis=0, ddof=1) / np.sqrt(n) / xs.mean(axis=0)) ** 2 + (xw.std(axis=0, ddof=1) / np.sqrt(n) / xw.mean(axis=0)) ** 2)
    L_ += [f"Part H: coloured (τ_c = {TC_H:g}) over white, E[u_k] − E[g_k∂_{{α_k}}V], with the delta method's errors (the runs as independent): "
           + ', '.join(f'{a:.3f} ± {b:.3f}' for a, b in zip(d, ed)) + "; from the linear response's τ_k(τ_c)/τ_k "
           + ', '.join(f'{x:.3f}' for x in H['ratio_col']) + " by " + ', '.join(f'{x:+.1f}' for x in (d - H['ratio_col']) / ed)
           + " errors. The coloured run's τ̂_k/τ_k: largest " + f"{(c['tk_hat'] / c['tk']).max():.4f}.", ""]
    zr = [z for z in F['zeros'] if z['variant'] == 'correlated E + class levels' and z['gamma'] == GAMMA_E][0]
    tg = e_target('correlated E + class levels')
    L_ += [f"Part F, 'correlated E + class levels': the zero of Katzenberger's corrected drift at γ = {GAMMA_E:g} minus W♯_τ, per class: "
           + ', '.join(f'{x:+.4f}' for x in zr['zero'] - tg) + " (the measured second half minus W♯_τ is in Part E's table).", ""]
    cw = [(r['gamma'], np.abs(r['centres and widths only']['b']).max(), np.abs(r['centres and widths only']['D']).max()) for r in C]
    gs = sorted(set(g for g, _, _ in cw), reverse=True)
    L_ += ["Part C, the centres and widths dithered only, largest over the points W: |drift|/γτ "
           + ', '.join(f"{max(b for g_, b, _ in cw if g_ == g):.2e} (γ = {g:g})" for g in gs)
           + f"; |D|/γ²τ at most {max(D for _, _, D in cw):.1e} over all γ and W.", ""]
    smp = i_samples(I); s = I['runs']['old']['s']
    L_ += ["Part I, the martingale-regressed first mass: its s.e. from the residuals alone and with b's estimation, √(s²(1/n + m̄²/S_mm)) "
           "(m̄ the sample's mean of m, S_mm its sum of squares about it):", "",
           "| sample | s | s.e. (residuals) | s.e. (with b's estimation) | their ratio |", "|---|---|---|---|---|"]
    rr = {}
    for lab, (Rr, acc) in smp.items():
        if len(Rr) < len(SEEDS_G):
            continue
        for t in T_I:
            i = int(np.argmin(np.abs(s - t))); wv = Rr[:, i, 4]; m = i_contrib(acc[:, i])[:, 1]; nn = len(wv)
            e = i_regressed(wv, m, np.zeros(2))[1]; e2 = e * np.sqrt(1 + nn * m.mean() ** 2 / ((m - m.mean()) ** 2).sum())
            L_.append(f"| {lab} | {s[i]:g} | {e:.5f} | {e2:.5f} | {e2 / e:.3f} |"); rr.setdefault(len(Rr), []).append(e2 / e)
    L_ += ["", "Their ratio over the samples of " + '; '.join(f"{n} seeds: {min(v):.3f} to {max(v):.3f}" for n, v in sorted(rr.items())) + ".", ""]
    return L_


def m_records(M, K):
    """Part K's 'nonlinear only' records at the limit map's indices (seeds × records): W, α, the exact Y, the convergence check."""
    r = K['var'][VARIANT_M]; ch = [M['lim_K'][c] for c in sorted(M['lim_K'])]; idx = ch[0]['idx']
    return r['W'][:, idx], r['al'][:, idx], np.concatenate([c['Y'] for c in ch]), np.concatenate([c['chk'] for c in ch]), idx


def m_loop_rates(M, g):
    """The coupled loops at the rate γ: per seed the limit map's rate over [T_YM[0], T_YM[-1]] per γ²τ, for the Δt chain, the Δt/2 chain and
    Richardson's 2r(Δt/2) − r(Δt); the limit map at T_YM[0] − 1/3; the convergence check."""
    ch = [M['loop'][g][c] for c in sorted(M['loop'][g])]
    Y = np.concatenate([c['Y'].reshape(len(c['t']), 2, -1, 3) for c in ch], axis=2)        # times × chains × seeds × 3
    t = ch[0]['t']; rt = (Y[-1] - Y[0]) / (t[-1] - t[0]) / (g ** 2 * TAU)
    return dict(full=rt[0], half=rt[1], rich=2 * rt[1] - rt[0], y0=Y[0], chk=max(float(np.abs(c['Y'] - c['chk']).max()) for c in ch),
                n=Y.shape[2], t=t, Y=Y)


def m_fixed_all(M, t):
    """The fixed-mass runs at the level τ, their seeds concatenated: m_anti for each chain and Richardson's."""
    ch = [M['fixed'][t][c] for c in sorted(M['fixed'][t])]
    cat = lambda c_: {k: np.concatenate([x[c_][k] for x in ch]) for k in ch[0][c_]}
    (af, sf), (ah, sh) = m_anti(cat('full'), t), m_anti(cat('half'), t)
    return dict(full=af, half=ah, rich=2 * ah - af, sym_full=sf, sym_half=sh, n=len(af), S=ch[0]['S'])


def report_M(M, K, E_rec, F, C):
    pe = lambda X: ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(*_ms(X)))
    p5 = lambda X: ', '.join(f'{a:+.5f} ± {b:.5f}' for a, b in zip(*_ms(X)))
    v5 = lambda v: ', '.join(f'{x:+.5f}' for x in v)
    L_ = ["## M. The missing O(γ²τ) drift of the deterministic-amplitude loop (D‴.40 after round 153)", "",
          f"'{VARIANT_M}' (the centres and widths dithered, the amplitudes noise-free; γ = {GAMMA_E:g}, τ = {TAU:g}): Katzenberger's limit "
          "gives the masses an O(γ²τ) drift, (γ/m)²Σ_lE[φ_l dφ_k/ds] at leading order in γ (Proposition D‴.40(ii), φ_k = "
          "(a_max²/4G_α)cosh(2α_k/a_max)), and Part K found the loop's secular displacement far below its drift path. The tests: (a) the "
          "prediction, from the linearised positions' antisymmetric correlation in closed form and from the limit map Φ's second differences "
          "(`m_limit`: the noise-free loop integrated by RK4 at the step 1 to s = 4000, then 8 to 4·10⁴; `m_fd`), against `limit_drift`; (b) "
          f"Part K's measurement with the exact limit map Y = Φ(state) at its records (every {EVERY_M}th of the second half) in place of its "
          "first-order offset (1); (c) the loop's own rate: new runs from the exact fit at two γ, each as two Euler chains at Δt and Δt/2 on "
          "one Brownian path (`m_loop`), Y's rate between "
          f"{T_YM[0]:g} and {T_YM[-1]:g} slow units, Richardson's 2r(Δt/2) − r(Δt); (d) the antisymmetric correlation at fixed masses at "
          "three levels τ (`m_fixed`, the same coupling). Drifts per γ²τ unless stated.", "",
          "### (a) The prediction at the masses 1/3", "",
          "| source | γ | the masses' drift per γ²τ |", "|---|---|---|"]
    ou = {o['tau']: o for o in M['ou']}; o = ou[TAU]
    L_.append(f"| the linearised positions' antisymmetric correlation (closed form; the limit γ → 0) | 0 | {v5(o['d'])} |")
    c5 = [r for r in C if r['gamma'] == min(GAMMAS_C) and np.allclose(r['w'], 1 / 3)][0]
    L_.append(f"| `limit_drift` (Part C) | {c5['gamma']:g} | {v5(c5['centres and widths only']['b'] / c5['gamma'])} |")
    for r in sorted(M['fd'], key=lambda r: r['gamma']):
        L_.append(f"| `limit_drift` | {r['gamma']:g} | {v5(r['kz'])} |")
        L_.append(f"| Φ's second differences, h = {', '.join(f'{x:g}' for x in r['h'])}; Richardson | {r['gamma']:g} | "
                  f"{'; '.join(v5(x) for x in r['fd'])}; {v5(r['rich'])} |")
    for g in sorted(M['katz']):
        if g not in [r['gamma'] for r in M['fd']]:
            L_.append(f"| `limit_drift` | {g:g} | {v5(M['katz'][g])} |")
    L_ += ["", "Φ's second differences (Richardson) against `limit_drift`, the largest relative difference: "
           + '; '.join(f"γ = {r['gamma']:g} {np.abs(r['rich'] / r['kz'] - 1).max():.4f}" for r in sorted(M['fd'], key=lambda r: r['gamma']))
           + "; the limit map at the exact fit moves the masses by at most "
           + f"{max(np.abs(r['at_fit']).max() for r in M['fd']):.1e}, and between s = {S_CHK_M:g} and the end by at most "
           + f"{max(r['chk'] for r in M['fd']):.1e}. The linearised positions (masses 1/3, τ = {TAU:g}): their rates "
           + f"{o['rates'][0]:.3e} to {o['rates'][-1]:.4f}, the amplitudes' standard deviations {', '.join(f'{x:.4f}' for x in o['sd_alpha'])}; "
           + f"the antisymmetric correlation's share E[δα_lδα̇_k]/√(E[δα_l²]E[δα̇_k²]) at most {o['rho']:.4f} (zero for a reversible law).", ""]
    # (b) Part K's records with the exact limit map
    W, al, Y, chk, idx = m_records(M, K)
    x = k_effects(K, E_rec, F, VARIANT_M); s = K['s']; kp = F['var'][VARIANT_M]['kpath']['half'][0] - 1 / 3
    off = k_offset(al, W); lin, Mm = m_offset_lin(al, W)
    hm = lambda A: A.mean(axis=1)
    yd = hm(Y - 1 / 3) - x['eff']; se = _ms(yd)[1]
    rows = [("W − 1/3 (measured)", hm(W - 1 / 3)), ("(1) Part K's positions' offset (first order, about the fit at W)", hm(off)),
            ("W − Y, Y the exact limit map", hm(W - Y)), ("(I − M)⁻¹(1): the fit's motion with the masses included", hm(lin)),
            ("Y − 1/3", hm(Y - 1 / 3)), ("(2) Part K's O(Δt) effect", x['eff']), ("Y − 1/3 − (2)", yd),
            ("for comparison, Part K's measured − 1/3 − (1) − (2)", hm(W - 1 / 3 - off) - x['eff'])]
    L_ += ["### (b) Part K's measurement with the exact limit map", "",
           f"Part K's {len(W)} seeds, the records {s[idx[0]]:g}, {s[idx[1]]:g}, …, {s[idx[-1]]:g} ({len(idx)} per seed); means over them, "
           "per seed. Katzenberger's SDE describes Y; Part K subtracted from W the first-order offset (1) about the fit at the current masses, "
           "but the masses' own change moves that fit, and its amplitudes' change moves the masses again: W − Y = (I − M)⁻¹(1), M = "
           "(γ/m)(I − W1ᵀ)diag(dΥ_k^fit/dW_k), dΥ_k^fit/dW_k = −h_k²/(G_αW_k), h_k = a_k/a′_k. At the masses 1/3 −diag M = "
           + ', '.join(f'{-v:.4f}' for v in np.diag(m_offset_lin(E.exact_rep(E.setup()[0], np.ones(3) / 3)[:, 0], np.ones(3) / 3)[1]))
           + " (the first class's amplitude near its saturation).", "",
           "| quantity | W_1, W_2, W_3 (± s.e.) |", "|---|---|"]
    L_ += [f"| {lab} | {pe(X)} |" for lab, X in rows]
    L_ += ["", f"Katzenberger's drift path (Part F) over the second half: {', '.join(f'{c:+.2e}' for c in kp)}; Y − 1/3 − (2) against it: "
           f"z = {', '.join(f'{v:+.1f}' for v in (yd.mean(axis=0) - kp) / se)}, their ratio {', '.join(f'{v:.2f}' for v in yd.mean(axis=0) / kp)} "
           f"(± {', '.join(f'{v:.2f}' for v in se / np.abs(kp))}). Per record, |W − Y − (I − M)⁻¹(1)| at most "
           f"{np.abs(W - Y - lin).max():.1e} (the O(γ²) remainder), |W − Y − (1)| up to {np.abs(W - Y - off).max():.1e}; the limit map between "
           f"s = {S_CHK_M:g} and the end at most {np.abs(Y - chk).max():.1e}.", ""]
    # (c) the coupled loops
    L_ += ["### (c) The loop's rate at two steps and two γ", "",
           f"New seeds {SEEDS_M[0]}–{SEEDS_M[-1]}, {S_M:g} slow units from the exact fit at the masses 1/3 (the box, τ = {TAU:g}); Y at s = "
           f"{', '.join(f'{t:g}' for t in T_YM)}; the rate (Y({T_YM[-1]:g}) − Y({T_YM[0]:g}))/{T_YM[-1] - T_YM[0]:g} per γ²τ, per seed.", "",
           "| γ | rate per γ²τ: the Δt chain | the Δt/2 chain | Richardson | Katzenberger's drift at 1/3 | z | Richardson over Katzenberger |",
           "|---|---|---|---|---|---|---|"]
    rr = {}
    for g in GAMMAS_M:
        q = m_loop_rates(M, g); rr[g] = q; kz = M['katz'][g]; m_, e_ = _ms(q['rich'])
        L_.append(f"| {g:g} | {p5(q['full'])} | {p5(q['half'])} | {p5(q['rich'])} | {v5(kz)} | "
                  f"{', '.join(f'{v:+.1f}' for v in (m_ - kz) / e_)} | {', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(m_ / kz, e_ / np.abs(kz)))} |")
    q = rr[GAMMA_E]; hs = s[len(s) // 2:].mean()
    L_ += ["", f"{', '.join(str(rr[g]['n']) for g in GAMMAS_M)} seeds. At γ = {GAMMA_E:g} the step's effect on the rate, r(Δt) − Richardson: "
           f"{p5(q['full'] - q['rich'])}; Part K's (2) as a rate (its second-half mean over the half's mean time {hs:g}): "
           f"{v5(x['eff'].mean(axis=0) / hs / (GAMMA_E ** 2 * TAU))}. Y({T_YM[0]:g}) − 1/3 (the first transient), Richardson: "
           + '; '.join(f"γ = {g:g} {pe(2 * rr[g]['y0'][1] - rr[g]['y0'][0] - 1 / 3)}" for g in GAMMAS_M)
           + f". The limit map between s = {S_CHK_M:g} and the end: at most {max(rr[g]['chk'] for g in GAMMAS_M):.1e}.", ""]
    # (d) fixed masses
    L_ += ["### (d) The positions' antisymmetric correlation at the fixed masses 1/3", "",
           f"γ = 0, the box; after {BURN_MF:g} slow units of burn-in; (1/m²)Σ_l of the antisymmetric part of E[φ_lφ̇_k] − E[φ_l]E[φ̇_k] over the "
           "states, per τ (the drift it gives the masses per γ²τ at γ → 0); the linearised positions (closed form) for the diffusion and for "
           "the Euler chains' laws.", "",
           "| τ | seeds × slow units | the Δt chain | the Δt/2 chain | Richardson | closed form: the diffusion | the Δt chain | the Δt/2 chain | "
           "Richardson over the closed form |", "|---|---|---|---|---|---|---|---|---|"]
    for t in TAUS_MF:
        q = m_fixed_all(M, t); o = ou[t]; m_, e_ = _ms(q['rich'])
        L_.append(f"| {t:g} | {q['n']} × {q['S']:g} | {p5(q['full'])} | {p5(q['half'])} | {p5(q['rich'])} | {v5(o['d'])} | {v5(o['chain'][0])} | "
                  f"{v5(o['chain'][1])} | {', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(m_ / o['d'], e_ / np.abs(o['d'])))} |")
    q = m_fixed_all(M, TAU); o = ou[TAU]
    L_ += ["", f"The symmetric part (the chain's bias; zero for the diffusion), τ = {TAU:g}: the Δt chain {p5(q['sym_full'])}, the Δt/2 chain "
           f"{p5(q['sym_half'])}; closed form {v5(o['chain_sym'][0])} and {v5(o['chain_sym'][1])}.", ""]
    zf = lambda m1, e1, m2, e2: ', '.join(f'{v:+.1f}' for v in (m1 - m2) / np.sqrt(e1 ** 2 + e2 ** 2))
    fr = {t: _ms(m_fixed_all(M, t)['rich']) for t in TAUS_MF}
    cr = {g: (_ms(rr[g]['rich'])[0] / M['katz'][g], _ms(rr[g]['rich'])[1] / np.abs(M['katz'][g])) for g in GAMMAS_M}
    L_ += ["Richardson against the closed form, z: " + '; '.join(f"τ = {t:g} {zf(fr[t][0], fr[t][1], ou[t]['d'], 0.0)}" for t in TAUS_MF)
           + f". (c)'s Richardson over Katzenberger against (d)'s Richardson over the closed form at τ = {TAU:g} "
           f"({', '.join(f'{v:.2f}' for v in fr[TAU][0] / ou[TAU]['d'])}), z: "
           + '; '.join(f"γ = {g:g} {zf(cr[g][0], cr[g][1], fr[TAU][0] / ou[TAU]['d'], fr[TAU][1] / np.abs(ou[TAU]['d']))}" for g in GAMMAS_M)
           + f". The γ dependence, the rate per γ²τ at γ = {GAMMAS_M[1]:g} over that at γ = {GAMMAS_M[0]:g}: Richardson "
           + ', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(*ratio_se(rr[GAMMAS_M[1]]['rich'], rr[GAMMAS_M[0]]['rich'])))
           + f", Katzenberger's {', '.join(f'{v:.2f}' for v in M['katz'][GAMMAS_M[1]] / M['katz'][GAMMAS_M[0]])}.", ""]
    # (e) Part K's shortfall, term by term
    q = rr[GAMMA_E]; ec = (q['full'] - q['rich']) * hs * GAMMA_E ** 2 * TAU
    ecm, ece = _ms(ec); kres = hm(W - 1 / 3 - off) - x['eff']; y2 = hm(Y - 1 / 3) - ecm; y2e = _ms(y2)[1]
    parts = [("Part K's secular displacement, measured − 1/3 − (1) − (2) (at (b)'s records)", *_ms(kres)),
             ("the first-order offset's error, (1) − (W − Y)", *_ms(off.mean(axis=1) - hm(W - Y))),
             ("the step correction's excess, (2) − (c)'s step effect", _ms(x['eff'])[0] - ecm, np.sqrt(_ms(x['eff'])[1] ** 2 + ece ** 2)),
             ("the rest, Katzenberger's path − (Y − 1/3 − (c)'s step effect)", kp - y2.mean(axis=0), np.sqrt(y2e ** 2 + ece ** 2))]
    tot = sum(p_[1] for p_ in parts)
    L_ += ["### (e) Part K's shortfall, term by term", "",
           f"Katzenberger's drift path over the second half (Part F) is the sum of four terms: Part K's secular displacement; the error of "
           "its first-order offset (1) against the exact W − Y ((b)); the excess of its step correction (2) over the loop's own step effect "
           f"(c) at γ = {GAMMA_E:g}, r(Δt) − Richardson per seed times γ²τ and the half's mean time {hs:g}; and the rest, which "
           f"(c) and (d) measure on their own. Each as a share of the path (± s.e.; the {len(W)} seeds of (b) and the {q['n']} of (c) independent).", "",
           "| term | W_1, W_2, W_3 (± s.e.) | share of Katzenberger's path (± s.e.) |", "|---|---|---|"]
    L_ += [f"| {lab} | {', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(m_, e_))} | "
           f"{', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(m_ / kp, e_ / np.abs(kp)))} |" for lab, m_, e_ in parts]
    L_ += ["", f"The four terms' sum − Katzenberger's path: at most {np.abs(tot - kp).max():.1e}. Y − 1/3 − (c)'s step effect: "
           f"{', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(y2.mean(axis=0), np.sqrt(y2e ** 2 + ece ** 2)))}, over Katzenberger's path "
           f"{', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(y2.mean(axis=0) / kp, np.sqrt(y2e ** 2 + ece ** 2) / np.abs(kp)))}; "
           f"against (c)'s Richardson over Katzenberger at γ = {GAMMA_E:g}, z: "
           f"{zf(y2.mean(axis=0) / kp, np.sqrt(y2e ** 2 + ece ** 2) / np.abs(kp), cr[GAMMA_E][0], cr[GAMMA_E][1])}.", ""]
    return L_


def n_ratios(M, Nn):
    """Per level τ: the per-seed ratio of Richardson's antisymmetric correlation to the closed form (Part M(a)'s at τ = TAU; it is
    τ-independent per τ) for Part M(d)'s seeds ('M', None if none), Part N's ('N') and both ('pool'); with the seeds × slow units."""
    d = {o['tau']: o for o in M['ou']}[TAU]['d']; out = {}
    for t in sorted(set(TAUS_FIT_N) | set(SEEDS_N)):
        q = {k: n_pool(*a, t) for k, a in (('M', (M, {})), ('N', ({}, Nn)), ('pool', (M, Nn)))}
        out[t] = {k: (None if v is None else dict(r=v['rich'] / d, n=v['n'], S=v['S'])) for k, v in q.items()}
    return out


def report_N(M, Nn):
    from scipy.stats import chi2
    R = n_ratios(M, Nn)
    r2 = lambda m, e: ', '.join(f'{a:.2f} ± {b:.2f}' for a, b in zip(m, e))
    z1 = lambda v: ', '.join(f'{x:+.1f}' for x in v)
    ns = lambda q: '—' if q is None else f"{q['n']} × {q['S']:g}"
    rs = lambda q: '—' if q is None else r2(*_ms(q['r']))
    L_ = ["## N. The finite-τ question: more seeds at fixed masses (D‴.40 after Part M)", "",
          "Part M(d)'s fixed-mass antisymmetric correlation (Richardson over the closed form) fell short of 1 at τ = 0.0001 and was "
          "undetermined at 0.0003. Is the loop's shortfall (Part M(c)) a finite-τ correction of the positions' law, a ratio falling with τ, "
          "or noise? `m_fixed` (unchanged: the masses 1/3, γ = 0, the box, two Euler chains at Δt and Δt/2 on one Brownian path, "
          f"{BURN_MF:g} slow units of burn-in) at new seeds: "
          + '; '.join(f"τ = {t:g} seeds {sd[0]}–{sd[-1]} × {S_N[t]:g} slow units" for t, sd in sorted(SEEDS_N.items()))
          + ". At 0.0001 and 0.0003 the run lengths are Part M(d)'s, so its seeds and the new ones pool with equal weight. Per seed the "
          "ratio of Richardson's 2r(Δt/2) − r(Δt) to the closed form (Part M(a), the limit τ → 0 of the linearised positions); means ± s.e. "
          "over the seeds.", "",
          "| τ | seeds × slow units: Part M(d) | Part N | the ratio: Part M(d)'s seeds | Part N's | their difference, z | pooled | "
          "pooled against 1, z | the seeds' spread of the ratio (standard deviation) |", "|---|---|---|---|---|---|---|---|---|"]
    R = {t: q for t, q in R.items() if q['pool'] is not None}
    for t in sorted(R):
        q = R[t]; mp, ep = _ms(q['pool']['r'])
        dz = '—'
        if q['M'] is not None and q['N'] is not None:
            (m1, e1), (m2, e2) = _ms(q['M']['r']), _ms(q['N']['r']); dz = z1((m2 - m1) / np.sqrt(e1 ** 2 + e2 ** 2))
        L_.append(f"| {t:g} | {ns(q['M'])} | {ns(q['N'])} | {rs(q['M'])} | {rs(q['N'])} | {dz} | {r2(mp, ep)} | {z1((mp - 1) / ep)} | "
                  f"{', '.join(f'{v:.2f}' for v in q['pool']['r'].std(axis=0, ddof=1))} |")
    pts = [_ms(R[t]['pool']['r']) for t in TAUS_FIT_N]
    f = n_fit(TAUS_FIT_N, [p[0] for p in pts], [p[1] for p in pts]); u = 1e-4
    L_ += ["", f"The ratio against τ, weighted least squares through the pooled points at τ = {', '.join(f'{t:g}' for t in TAUS_FIT_N)} (per class "
           "W_1, W_2, W_3; the third class's closed form is small, its ratio undetermined):", "",
           "| fit | the limit τ → 0 | its z against 1 | the slope, per 10⁻⁴ of τ | its z against 0 | χ² (degrees of freedom; p) |",
           "|---|---|---|---|---|---|",
           f"| the line c + bτ | {r2(f['c'], f['ec'])} | {z1((f['c'] - 1) / f['ec'])} | {r2(f['b'] * u, f['eb'] * u)} | {z1(f['b'] / f['eb'])} | "
           f"{', '.join(f'{x:.2f}' for x in f['chi'])} ({f['dof']}; p = {', '.join(f'{x:.2f}' for x in chi2.sf(f['chi'], f['dof']))}) |",
           f"| a constant | {r2(f['m'], f['em'])} | {z1((f['m'] - 1) / f['em'])} | — | — | "
           f"{', '.join(f'{x:.2f}' for x in f['chim'])} ({f['dof'] + 1}; p = {', '.join(f'{x:.2f}' for x in chi2.sf(f['chim'], f['dof'] + 1))}) |", ""]
    tx = max(t for t in R if t not in TAUS_FIT_N) if any(t not in TAUS_FIT_N for t in R) else None
    if tx is not None and R[tx]['pool'] is not None:
        pr = f['c'] + f['b'] * tx; pe_ = np.sqrt(f['ec'] ** 2 + (tx * f['eb']) ** 2 + 2 * tx * f['cov']); mx, ex = _ms(R[tx]['pool']['r'])
        L_ += [f"At τ = {tx:g} ({R[tx]['pool']['n']} seeds, outside the fit): the ratio {r2(mx, ex)}, the line's prediction {r2(pr, pe_)}, "
               f"z = {z1((mx - pr) / np.sqrt(ex ** 2 + pe_ ** 2))}.", ""]
    return L_


def report(A=None, B=None, C=None, D=None, Ee=None, E_rec=None, F=None, G=None, H=None, I=None, J=None, K=None, M=None, N=None):
    fmt = E.fmt
    v3 = lambda v, p=4: ', '.join(f'{x:.{p}f}' for x in v)
    L_ = ["# Proposition D‴.40: the gauge under general noise (`exp_gauge_noise.py`)", "",
          f"Proposition D‴.36's instance (three bumps, one particle per class, m = 3); the noises at τ = {TAU:g}: isotropic 2τ diag(G); the "
          "amplitudes only; the centres and widths only; amplitude levels (2, 1, 0.5)τ; three correlated covariances 2τ·0.3·LLᵀ (L standard "
          f"normal, seeds {', '.join(str(s) for s in NOISE_SEEDS)}). τ_k = Σ_{{α_kα_k}}/(2G_α); the closed form −(τ_k/m)(1 + (a_k/a_max)²)/(W_k a′_k).", ""]
    if A is not None:
        L_ += ["## A. The trace formula against the closed form", "",
               "| noise | W | τ_k / τ | ū (trace formula) | the closed form | their replicator flows' largest relative difference |", "|---|---|---|---|---|---|"]
        for r in A:
            zero = not np.any(levels(noises()[r['name']]))
            L_.append(f"| {r['name']} | {v3(np.r_[r['w'], 1 - sum(r['w'])], 3)} | {', '.join(f'{v:.3f}' for v in levels(noises()[r['name']]) / TAU)} | "
                      f"{', '.join(fmt(v) for v in r['ubar'])} | {', '.join(fmt(v) for v in r['closed'])} | "
                      f"{(fmt(np.abs(r['flow']).max() * TAU) + ' (the flow itself; the closed form is 0)') if zero else fmt(r['rel'])} |")
        L_.append("")
    if B is not None:
        L_ += ["## B. Katzenberger's limit drift at γ = 10⁻⁵ against the replicator flow of the closed form", "",
               "| noise | W | the drift per unit slow time | the closed form's flow (γ = 10⁻⁵) | largest relative difference | the diffusion's diagonal |",
               "|---|---|---|---|---|---|"]
        for r in B:
            den = max(np.abs(r['pred']).max(), 1e-30)
            L_.append(f"| {r['name']} | {v3(np.r_[r['w'], 1 - sum(r['w'])], 3)} | {', '.join(fmt(v) for v in r['b'])} | {', '.join(fmt(v) for v in r['pred'])} | "
                      f"{fmt(np.abs(r['b'] - r['pred']).max() / den) if np.abs(r['pred']).max() > 0 else fmt(np.abs(r['b']).max())} | "
                      f"{', '.join(fmt(v) for v in np.diag(r['D']))} |")
        L_.append("")
    if C is not None:
        L_ += ["## C. The centres and widths dithered only, at finite γ", "",
               "| γ | W | drift/γτ, isotropic | drift/γτ, centres and widths only | D/γ²τ, isotropic (diagonal) | D/γ²τ, centres and widths only |",
               "|---|---|---|---|---|---|"]
        for r in C:
            iso, cs = r['isotropic'], r['centres and widths only']
            L_.append(f"| {r['gamma']:g} | {v3(np.r_[r['w'], 1 - sum(r['w'])], 3)} | {', '.join(f'{v:+.5f}' for v in iso['b'])} | {', '.join(f'{v:+.2e}' for v in cs['b'])} | "
                      f"{', '.join(f'{v:.4f}' for v in np.diag(iso['D']))} | {np.abs(cs['D']).max():.1e} |")
        L_.append("")
    if D is not None:
        se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M))
        lng = lambda M: M.reshape(N_BATCH_LONG, -1, 3).mean(axis=1)
        pm = lambda M: ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(lng(M).mean(axis=0), se(lng(M))))
        L_ += [f"## D. The identity in simulation (fixed masses 1/3, the box of Proposition D‴.39, {S_D:g} slow units)", "",
               "The generator applied to ψ_k = ∫g_k dα_k/G_α = (a_max²/4β_kG_α)cosh(2α_k/a_max) gives, in the stationary law, "
               "E[u_k] = −τ_kE[g′_k] + E[g_k∂_{α_k}J_coh] + E[g_k∂_{α_k}V] + E[Lψ_k], the last the Euler chain's bias (zero for the diffusion); "
               "(i)'s box-free form leaves out the cohesion's and the box's terms. Errors: batch means over "
               f"{N_BATCH_LONG} batches of {S_D / N_BATCH_LONG:g} slow units (the plateau gives a slow tail: referee round 31).", "",
               "| noise | τ | τ_k/τ | E[u_k] | E[−τ_kg′_k] | E[g_k∂_{α_k}J_coh] | E[g_k∂_{α_k}V] (the box) | (i)'s residual E[u_k + τ_kg′_k] (z) | "
               "α_k outside the box (share of steps) |", "|---|---|---|---|---|---|---|---|---|"]
        zi = []
        for r in D:
            res = lng(r['u'] - r['gp']); z = res.mean(axis=0) / se(res); zi += list(z)
            L_.append(f"| {r['name']} | {r['tau']:g} | {', '.join(f'{v:.3f}' for v in r['tk'] / r['tau'])} | {pm(r['u'])} | {pm(r['gp'])} | "
                      f"{pm(r['coh'])} | {pm(r['box'])} | {', '.join(f'{a:+.2e} ({b:+.1f})' for a, b in zip(res.mean(axis=0), z))} | "
                      f"{', '.join(f'{v:.3f}' for v in r['out_box'])} |")
        L_ += ["", "The whole identity and the Euler chain's bias (the chain's exact per-step bias b_k = E_n[ψ_k(α_(n+1)) − ψ_k(α_n)]/Δt − Lψ_k, "
               "the Gaussian increment's expectation in closed form, so that E[Lψ_k] = −E[b_k] exactly; the error of E[Lψ_k] from the "
               "martingale increments ψ_k(α_(n+1)) − E_n[ψ_k(α_(n+1))]):", "",
               "| noise | τ | E[Lψ_k] (± its martingale error) | −E[b_k] | E[Lψ_k] + E[b_k] over that error | of which the boundary term "
               "(ψ_k(end) − ψ_k(start))/T | the bias E[b_k] over τ_kE[g′_k] | "
               "τ_k/τ̂_k: the amplitudes' quadratic variation, with the drift's part, without |", "|---|---|---|---|---|---|---|---|"]
        zf, zb, det = [], [], []
        for r in D:
            lp, bi = r['lpsi'].mean(axis=0), r['bias'].mean(axis=0)
            bd = r.get('boundary', np.full(3, np.nan))
            if np.all(r['se_mart'] > 0):
                z = (lp + bi) / r['se_mart']; zf += list(z); zb += list(bd / r['se_mart'])
                zs, bs = ', '.join(f'{v:+.2f}' for v in z), ', '.join(f'{v:+.2f}' for v in bd / r['se_mart'])
            else:                                                     # no noise on the amplitudes: the martingale is zero
                det.append((r['name'], np.abs(lp + bi - bd).max()))
                zs = '— (no amplitude noise: ' + ', '.join(f'{v:+.2e}' for v in lp + bi) + ')'
                bs = ', '.join(f'{v:+.2e}' for v in bd)
            L_.append(f"| {r['name']} | {r['tau']:g} | {', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(lp, r['se_mart']))} | "
                      f"{', '.join(f'{-v:+.2e}' for v in bi)} | {zs} | {bs} | "
                      f"{', '.join(f'{v:+.4f}' for v in bi / -r['gp'].mean(axis=0)) if np.all(r['tk'] > 0) else '— (τ_k = 0)'} | "
                      f"{', '.join(f'{a / b:.4f}, {a / c:.4f}' for a, b, c in zip(r['tk'], r['tk_hat'], r['tk_hat_noise'])) if np.all(r['tk'] > 0) else '—'} |")
        zi, zf, zb = np.array(zi), np.array(zf), np.array(zb)
        L_ += ["", f"The chain's bookkeeping, a check of the averaging and not of the identity (the runs with noise on the amplitudes; Part J "
               f"checks the identity): χ² = {float(zf @ zf):.1f} for {len(zf)} values (largest |z| "
               f"{np.abs(zf).max():.2f}; {int((zf > 0).sum())} of them positive, their mean {zf.mean():+.2f}); E[Lψ_k] + E[b_k] is exactly the boundary "
               "term (ψ_k(end) − ψ_k(start))/T less the martingale's mean, so that without the boundary term χ² = "
               f"{float((zf - zb) @ (zf - zb)):.1f} (largest |z| {np.abs(zf - zb).max():.2f}, mean {(zf - zb).mean():+.2f})"
               + (''.join(f"; '{n}': E[Lψ_k] + E[b_k] is the boundary term to {e:.1e}" for n, e in det)) + ". (i)'s box-free form: "
               f"the largest |z| {np.abs(zi).max():.1f}, its residual the box's term.", ""]
    if Ee is not None:
        ws = Q.w_sharp(); wt = w_sharp_weighted(levels(noises()['correlated E']))
        target = {'nonlinear only': np.ones(3) / 3, 'nonlinear only + barrier': ws, 'correlated E': None, 'correlated E + class levels': wt}
        L_ += [f"## E. The full loop inside the box (γ = {GAMMA_E:g}, τ = {TAU:g}, {len(SEEDS_E)} seeds, {S_E:g} slow units)", "",
               f"Targets: the initial masses 1/3 (frozen), W♯ = {v3(ws)}, W♯_τ = {v3(wt)} (the maximiser of Σ_kτ_k log(W_k a′_k) for "
               f"'correlated E', τ_k/τ = {', '.join(f'{v:.3f}' for v in levels(noises()['correlated E']) / TAU)}, its centres and widths at the "
               f"dither's level); the box's representable masses W_B = {v3(P.W_BOX)}.", "",
               "| variant | s | mean of W_1, W_2, W_3 (± s.e.) | std |", "|---|---|---|---|"]
        for r in Ee:
            for t in r['sum']['times']:
                L_.append(f"| {r['variant']} | {t['s']:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(t['mean'], t['se']))} | {v3(t['std'])} |")
        L_ += ["", "| variant | the second half's mean (± s.e.) | its target | the distance from it (largest component) | the distance from W♯ | "
               "the largest |α| | the fit error (second half) |", "|---|---|---|---|---|---|---|"]
        for r in Ee:
            hm = r['sum']['half']; se = hm.std(axis=0, ddof=1) / np.sqrt(len(hm)); tg = target[r['variant']]
            dist = '—' if tg is None else f"{np.abs(hm.mean(axis=0) - tg).max():.4f}"
            L_.append(f"| {r['variant']} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(hm.mean(axis=0), se))} | "
                      f"{'none' if tg is None else v3(tg)} | {dist} | {np.abs(hm.mean(axis=0) - ws).max():.4f} | {r['sum']['almax']:.2f} | "
                      f"{fmt(r['sum']['fit'])} |")
        L_.append("")
    if F is not None and E_rec is not None:
        L_ += report_F(F, E_rec)
    if G is not None and E_rec is not None and F is not None:
        L_ += report_G(G, E_rec, F)
    if H is not None:
        L_ += report_H(H)
    if I is not None and E_rec is not None and F is not None and G is not None:
        L_ += report_I(I, E_rec, F, G)
    if J is not None:
        L_ += report_J(J)
    if K is not None and E_rec is not None and F is not None:
        L_ += report_K(K, E_rec, F)
    if F is not None and H is not None and I is not None and C is not None:
        L_ += report_L(F, H, I, C)
    if M is not None and K is not None and E_rec is not None and F is not None and C is not None:
        L_ += report_M(M, K, E_rec, F, C)
    if N is not None and M is not None:
        L_ += report_N(M, N)
    L_.append("Numbers are produced by `experiments/exp_gauge_noise.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def _ms(W):
    """seeds × 3 → (mean, s.e.)."""
    return W.mean(axis=0), W.std(axis=0, ddof=1) / np.sqrt(len(W))


def zden(se, sde_std, n):
    """The error of an n-seed mean against Katzenberger's SDE: the larger of its own s.e. and the SDE's std/√n."""
    d = np.maximum(se, sde_std / np.sqrt(n))
    return np.where(d > 0, d, np.inf)


def report_F(F, E_rec):
    v3 = lambda v, p=4: ', '.join(f'{x:.{p}f}' for x in v)
    pm = lambda m, e: ', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(m, e))
    L_ = ["## F. Part E against its predictions (referee round 147)", "",
          "Each variant of Part E beside (a) the slow limit at leading order — the replicator flow of the closed form −τ_kg′_k plus the masses "
          "form's κ_Fτ_F,kg′_k, from the masses 1/3 —, (b) Katzenberger's limit at γ = "
          f"{GAMMA_E:g}: its drift's path, and the ensemble of its SDE ({X.N_PATHS} paths from 1/3, Euler–Maruyama at the step {X.DS_SDE:g}; "
          f"the drift and diffusion computed at the {len(grid_PB())} points of a barycentric grid of P_B and interpolated linearly). P_B's walls "
          "W_B are enforced as the box enforces them: at a wall a drift out of P_B is cancelled by raising that class's fitness (the paths "
          "reflected there). The measured columns are Part E's 12 seeds; z = (measured − the SDE's mean) over the larger of the measured s.e. "
          "and the SDE's std/√12 (a 12-seed mean's s.e. if the seeds were the SDE's paths), the largest |z| of the three classes.", "",
          "| variant | the slow limit's drift at W = 1/3 per γτ | Katzenberger's drift at W = 1/3 per γτ (γ = 0.01) | their ratio (class 1) | "
          "the slow limit first at a wall |", "|---|---|---|---|---|"]
    for v in VARIANTS_E:
        r = F['var'][v]
        rat = f"{r['b13'][0] / r['lead13'][0]:.3f}" if abs(r['lead13'][0]) > 0 else '—'
        tw = '—' if r['ode']['t_wall'] is None else f"s = {r['ode']['t_wall']:g}"
        L_.append(f"| {v} | {', '.join(f'{x:+.4f}' for x in r['lead13'])} | {', '.join(f'{x:+.4f}' for x in r['b13'])} | {rat} | {tw} |")
    L_ += ["", "| variant | s | measured mean (± s.e.) | the slow limit (ODE) | Katzenberger's drift path | Katzenberger's SDE: mean (std) | "
           "largest \\|z\\| |", "|---|---|---|---|---|---|---|"]
    for v in VARIANTS_E:
        r = F['var'][v]; R = np.asarray(E_rec[v]); s = R[0, :, 0]
        for t in T_F:
            i = int(np.argmin(np.abs(s - t))); m, e = _ms(R[:, i, 4:7])
            sm, sd = r['sde']['rec'][t].mean(axis=0), r['sde']['rec'][t].std(axis=0, ddof=1)
            z = np.abs(m - sm) / zden(e, sd, len(R))
            L_.append(f"| {v} | {s[i]:g} | {pm(m, e)} | {v3(r['ode']['rec'][t])} | {v3(r['kpath']['rec'][t][0])} | "
                      f"{', '.join(f'{a:.4f} ({b:.4f})' for a, b in zip(sm, sd))} | {z.max():.1f} |")
        hm = R[:, R.shape[1] // 2:, 4:7].mean(axis=1); m, e = _ms(hm); sh = r['sde']['half']
        z = np.abs(m - sh.mean(axis=0)) / zden(e, sh.std(axis=0, ddof=1), len(R))
        L_.append(f"| {v} | the second half | {pm(m, e)} | {v3(r['ode']['half'])} | {v3(r['kpath']['half'][0])} | "
                  f"{', '.join(f'{a:.4f} ({b:.4f})' for a, b in zip(sh.mean(axis=0), sh.std(axis=0, ddof=1)))} | {z.max():.1f} |")
    L_ += ["", "The zeros of Katzenberger's corrected drift (the masses form included; Newton from the target):", "",
           "| variant | γ | the zero | its target | the distance (largest component) | the drift at the target per γτ | the residual per γτ |",
           "|---|---|---|---|---|---|---|"]
    for z in F['zeros']:
        tg = e_target(z['variant'])
        L_.append(f"| {z['variant']} | {z['gamma']:g} | {v3(z['zero'])} | {v3(tg)} | {np.abs(z['zero'] - tg).max():.4f} | "
                  f"{', '.join(f'{x:+.2e}' for x in z['at_target'])} | {z['resid']:.1e} |")
    L_ += ["", f"The weighted-barrier maximiser W♯_τ against P_B (W_B = {v3(P.W_BOX)}): (iii)'s target under the box needs W♯_τ ∈ P_B.", "",
           "| τ_k/τ | W♯_τ | W♯_τ − W_B | in P_B |", "|---|---|---|---|"]
    for p in F['pb']:
        L_.append(f"| {', '.join(f'{x:g}' for x in p['levels'])} | {v3(p['w'])} | {', '.join(f'{x:+.4f}' for x in p['margin'])} | "
                  f"{'yes' if np.all(p['margin'] > 0) else 'no'} |")
    L_ += ["", "Part D's last column, τ_k/τ̂_k, lists per class the pair (with the drift's part, without it).", ""]
    return L_


def report_G(G, E_rec, F):
    pm = lambda m, e: ', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(m, e))
    L_ = [f"## G. The early phase: {len(SEEDS_G)} new seeds over {S_G:g} slow units (referee round 147)", "",
          f"The full loop of Part E (the box, γ = {GAMMA_E:g}, τ = {TAU:g}) from the masses 1/3 with the seeds {SEEDS_G[0]}–{SEEDS_G[-1]}; "
          f"Part E's seeds {', '.join(str(s) for s in REPLAY_G)} replayed over the same span reproduce its records to "
          f"{G['replay_diff']:.1e}. z = (the new seeds' mean − Katzenberger's SDE mean (Part F)) over the larger of their s.e. and the SDE's "
          f"std/√{len(SEEDS_G)}.", "",
          "| variant | s | the new seeds: mean (± s.e.) | std | Part E's 12 seeds | all | Katzenberger's SDE mean | the slow limit (ODE) | z |",
          "|---|---|---|---|---|---|---|---|---|"]
    v3 = lambda v, p=4: ', '.join(f'{x:.4f}' for x in v)
    times = [t for t in T_F if t <= S_G]
    for v in VARIANTS_G:
        R = np.asarray(G['rec'][v]); s = R[0, :, 0]; Ro = np.asarray(E_rec[v]); so = Ro[0, :, 0]
        for t in times:
            i = int(np.argmin(np.abs(s - t))); io = int(np.argmin(np.abs(so - t)))
            m, e = _ms(R[:, i, 4:7]); mo, eo = _ms(Ro[:, io, 4:7]); ma, ea = _ms(np.r_[R[:, i, 4:7], Ro[:, io, 4:7]])
            sm, sd = F['var'][v]['sde']['rec'][t].mean(axis=0), F['var'][v]['sde']['rec'][t].std(axis=0, ddof=1)
            L_.append(f"| {v} | {s[i]:g} | {pm(m, e)} | {v3(R[:, i, 4:7].std(axis=0, ddof=1))} | {pm(mo, eo)} | {pm(ma, ea)} | {v3(sm)} | "
                      f"{v3(F['var'][v]['ode']['rec'][t])} | {', '.join(f'{x:+.1f}' for x in (m - sm) / zden(e, sd, len(R)))} |")
    a, b = (np.asarray(G['rec'][v]) for v in VARIANTS_G)
    L_ += ["", f"The paired difference ('{VARIANTS_G[1]}' − '{VARIANTS_G[0]}', the same seeds and noise):", "",
           "| s | measured (± s.e.) | Katzenberger's SDE means' difference | the slow limit's difference |", "|---|---|---|---|"]
    s = a[0, :, 0]
    for t in times:
        i = int(np.argmin(np.abs(s - t))); m, e = _ms(b[:, i, 4:7] - a[:, i, 4:7])
        ds = F['var'][VARIANTS_G[1]]['sde']['rec'][t].mean(axis=0) - F['var'][VARIANTS_G[0]]['sde']['rec'][t].mean(axis=0)
        do = F['var'][VARIANTS_G[1]]['ode']['rec'][t] - F['var'][VARIANTS_G[0]]['ode']['rec'][t]
        L_.append(f"| {s[i]:g} | {pm(m, e)} | {', '.join(f'{x:+.4f}' for x in ds)} | {', '.join(f'{x:+.4f}' for x in do)} |")
    L_ += ["", "The share of the new seeds' records (every 10³ slow units) with an amplitude outside the box, per class: "
           + '; '.join(f"'{v}' " + ', '.join(f'{x:.4f}' for x in (np.abs(np.asarray(G['rec'][v])[:, :, 1:4]) > P.BOX[0]).mean(axis=(0, 1)))
                       for v in VARIANTS_G) + ".", ""]
    return L_


def report_H(H):
    lin = H['linear']
    L_ = ["## H. Coloured noise (referee round 147)", "",
          "An Ornstein–Uhlenbeck noise η of correlation time τ_c with the white noise's integrated covariance Σ (η's covariance "
          "Σe^{−|s−s′|/τ_c}/(2τ_c)). The positions then move by dϑ = −G∇(J + V)ds + η ds and the gauge identity gives, in a stationary law, "
          "E[u_k] = −E[g_kη_{α_k}]/G_α + E[g_k∂_{α_k}(J_coh + V)] — no Itô term —; at linear response E[δϑ ηᵀ] = (I + τ_cGH)⁻¹Σ/2, so the level "
          "τ_k is replaced by τ_k(τ_c) = [(I + τ_cGH)⁻¹Σ]_{α_kα_k}/(2G_α). At the masses 1/3 the rates of GH run from "
          f"{lin['rates'][0]:.3e} to {lin['rates'][-1]:.3f}; the closed form agrees with the Lyapunov solution of the system augmented with η "
          f"to {lin['check']:.1e}.", "",
          "| noise | τ_c | τ_k(τ_c)/τ_k | the correlations' share \\|τ_k(Σ) − τ_k(diag Σ)\\|/τ_k(Σ) |", "|---|---|---|---|"]
    for r in lin['rows']:
        L_.append(f"| {r['name']} | {r['tc']:g} | {', '.join(f'{x:.4f}' for x in r['ratio'])} | {', '.join(f'{x:.4f}' for x in r['corr'])} |")
    mc = {r['tc']: r for r in H['mc']}
    w, c = mc[0.0], mc[TC_H]
    se = lambda M: M.std(axis=0, ddof=1) / np.sqrt(len(M))
    pm = lambda M: ', '.join(f'{a:+.2e} ± {b:.1e}' for a, b in zip(M.mean(axis=0), se(M)))
    L_ += ["", f"In simulation (fixed masses 1/3, the box, the amplitudes only at τ = {TAU_H:g}, {w['S']:g} slow units, errors from "
           f"{len(w['u'])} batch means; the white run and the coloured run at τ_c = {TC_H:g} from the same seed {w['seed']}). The leading "
           "order: −τ_kg′_k at the fit (white), −τ_k(τ_c)g′_k (coloured).", "",
           "| noise | E[u_k] − E[g_k∂_{α_k}V] | the leading order | the identity's residual (z) | τ̂_k/τ_k (the quadratic variation at Δt = "
           f"{DT:g}) |", "|---|---|---|---|---|"]
    for lab, r, lead in (('white', w, H['lead_white']), (f'coloured, τ_c = {TC_H:g}', c, H['lead_col'])):
        res = r['u'] - r['box'] - r['coh'] - (r['gp'] if r['tc'] == 0 else r['geta'])
        z = res.mean(axis=0) / se(res)
        L_.append(f"| {lab} | {pm(r['u'] - r['box'])} | {', '.join(f'{x:+.2e}' for x in lead)} | "
                  f"{', '.join(f'{a:+.1e} ({b:+.1f})' for a, b in zip(res.mean(axis=0), z))} | {', '.join(f'{x:.4f}' for x in r['tk_hat'] / r['tk'])} |")
    d = (c['u'] - c['box']).mean(axis=0) / (w['u'] - w['box']).mean(axis=0)
    zd = ((c['u'] - c['box']).mean(axis=0) - (w['u'] - w['box']).mean(axis=0)) / np.sqrt(se(c['u'] - c['box']) ** 2 + se(w['u'] - w['box']) ** 2)
    L_ += ["", f"The white run's residual is the Euler chain's O(Δt) bias with −τ_kE[g′_k] (Part D); the coloured run's identity is "
           "−E[g_kη_{α_k}]/G_α. Coloured over white, E[u_k] − E[g_k∂_{α_k}V]: " + ', '.join(f'{x:.3f}' for x in d)
           + " (the linear response's τ_k(τ_c)/τ_k: " + ', '.join(f'{x:.3f}' for x in H['ratio_col']) + "); their difference over the errors "
           "combined as if independent: " + ', '.join(f'{x:+.1f}' for x in zd) + ".", ""]
    return L_


def _p2(z):
    """The two-sided normal p-value of z."""
    import math
    return math.erfc(abs(z) / math.sqrt(2))


def holm(p):
    """Holm's step-down adjustment of the p-values p."""
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); out = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, min(1.0, (m - r) * p[i])); out[i] = run
    return out


def i_samples(I):
    """Part I's samples of the loop at T_I: label → (records (seeds, times, ≥ 7: s, α, W), accumulators or None)."""
    R = I['runs']; mk = lambda r: np.concatenate([np.broadcast_to(r['s'][None, :, None], r['W'].shape[:2] + (1,)), r['al'], r['W']], axis=2)
    old = mk(R['old']); nG = len(SEEDS_G)
    return {f"Part G's {nG} seeds ({SEEDS_G[0]}–{SEEDS_G[-1]})": (old[:nG], R['old']['acc'][:nG]),
            f"Part E's {len(SEEDS_E)} seeds ({SEEDS_E[0]}–{SEEDS_E[-1]})": (old[nG:], R['old']['acc'][nG:]),
            f"the {len(SEEDS_I)} new seeds ({SEEDS_I[0]}–{SEEDS_I[-1]})": (mk(R['box']), R['box']['acc']),
            f"the {len(SEEDS_I)} new seeds, no box": (mk(R['no box']), R['no box']['acc']),
            f"the {len(SEEDS_I)} new seeds, equilibrated start": (mk(R['burn-in']), R['burn-in']['acc'])}


def report_I(I, E_rec, F, G):
    pm = lambda m, e: ', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(m, e))
    v3 = lambda v: ', '.join(f'{x:.4f}' for x in v)
    R = I['runs']; old = R['old']; nG = len(SEEDS_G)
    Gr = np.asarray(G['rec']['correlated E']); Er = np.asarray(E_rec['correlated E'])[:, :old['W'].shape[1]]
    rep = max(np.abs(old['W'][:nG] - Gr[:, :, 4:7]).max(), np.abs(old['al'][:nG] - Gr[:, :, 1:4]).max(),
              np.abs(old['W'][nG:] - Er[:, :, 4:7]).max(), np.abs(old['al'][nG:] - Er[:, :, 1:4]).max())
    S = I['sde']; ref = 'reflected at W_B, step 50'; soft = f'the box measured, step {DS_SOFT_I[0]:g}'
    smp = i_samples(I); s = old['s']
    L_ = [f"## I. The early phase's gap (round 147's open item)", "",
          f"Part E's 'correlated E' loop (uncorrected, the box, γ = {GAMMA_E:g}, τ = {TAU:g}) over {S_I:g} slow units, integrated in batches "
          f"of {B_I} seeds (`i_runs`, gauge_batch's drift, each seed run_dithered's draws: Part G's {nG} and Part E's {len(SEEDS_E)} seeds "
          f"reproduce their records to {rep:.1e}); {len(SEEDS_I)} new seeds, the same seeds without the box, and the same seeds and draws "
          f"from the positions' equilibrated state (first {BURN_I:g} slow units at the fixed masses 1/3, the box, γ = 0). Katzenberger's SDE "
          f"at γ = {GAMMA_E:g} ({X.N_PATHS} paths from 1/3, seed {SEED_SDE_I}): reflected at P_B's walls W_B at the steps "
          f"{', '.join(f'{d:g}' for d in DS_I)} (step 50 is Part F's ensemble); without walls (the field on a grid of {len(grid_ext())} points "
          f"on W_k ≥ W_min,k + 0.012, reflected only at its edge); and with the box's term measured at fixed masses (the table at the end) "
          "added to the drift through Katzenberger's projection at γ (a force on α_1 alone), without walls. z = (the sample's mean − the "
          "model's mean) over the larger of the sample's s.e. and the model's std/√n.", "",
          "| sample or model | s | W_1, W_2, W_3: mean (± s.e.) | std of W_1 | W_1 < 0.30 (share) | z of W_1 against Part F's SDE | "
          "against the SDE with the box measured |", "|---|---|---|---|---|---|---|"]
    zf = {}
    for lab, (Rr, _) in smp.items():
        n = len(Rr)
        for t in T_I:
            i = int(np.argmin(np.abs(s - t))); m, e = _ms(Rr[:, i, 4:7])
            z = [(m[0] - S[k]['rec'][t][:, 0].mean()) / zden(e[0], S[k]['rec'][t][:, 0].std(ddof=1), n) for k in (ref, soft)]
            zf[(lab, t)] = z
            L_.append(f"| {lab} | {s[i]:g} | {pm(m, e)} | {Rr[:, i, 4].std(ddof=1):.4f} | {(Rr[:, i, 4] < 0.30).mean():.4f} | "
                      f"{z[0]:+.1f} | {z[1]:+.1f} |")
    for lab, r in S.items():
        for t in T_I:
            w = r['rec'][t]
            L_.append(f"| Katzenberger's SDE, {lab} | {t:g} | {v3(w.mean(axis=0))} | {w[:, 0].std(ddof=1):.4f} | {(w[:, 0] < 0.30).mean():.4f} | — | — |")
    L_ += ["", "The share of the SDE's paths that reached a wall by s = 2·10⁴: " + '; '.join(f"{k} {r['hit']:.4f}" for k, r in S.items()) + ".", ""]
    # multiple testing
    labs = list(smp)[:3]
    L_ += ["The first mass against Part F's SDE, three recorded times per sample (their z are positively correlated, the same paths; "
           "W_3's mirror W_1's, the masses summing to 1): two-sided p and Holm's adjustment over the three times.", "",
           "| sample | z at 5000, 10000, 20000 | p | Holm-adjusted p | the smallest adjusted p |", "|---|---|---|---|---|"]
    for lab in labs:
        z = [zf[(lab, t)][0] for t in T_I]; p = [_p2(x) for x in z]; h = holm(p)
        L_.append(f"| {lab} | {', '.join(f'{x:+.2f}' for x in z)} | {', '.join(f'{x:.4f}' for x in p)} | {', '.join(f'{x:.4f}' for x in h)} | {h.min():.4f} |")
    a, b = smp[labs[0]][0], smp[labs[2]][0]; i = int(np.argmin(np.abs(s - S_I)))
    (ma, ea), (mb, eb) = _ms(a[:, i, 4:7]), _ms(b[:, i, 4:7])
    L_ += ["", f"At s = {S_I:g} Part G's {len(a)} seeds and the {len(b)} new seeds differ in W_1 by {ma[0] - mb[0]:+.4f} ± "
           f"{np.hypot(ea[0], eb[0]):.4f} (z = {(ma[0] - mb[0]) / np.hypot(ea[0], eb[0]):+.2f}).", ""]
    # the decomposition
    L_ += ["The gauge identity's decomposition of each run's log-masses, Δlog(β_k/β_3) = γΣ_n(u_k − u_3)Δt with u_kΔt = the boundary term "
           "Δφ_k/β_k − the martingale term (φ_k(α_(n+1)) − E_nφ_k(α_(n+1)))/β_k + (−τ_kg′_k + g_k∂_{α_k}V + g_k∂_{α_k}J_coh + the chain's bias)Δt "
           "(φ_k = (a_max²/4G_α)cosh(2α_k/a_max); the martingale's mean is zero exactly), each term's contribution to W_1 linearised at the "
           "masses 1/3 (∂W/∂log β = W_k(δ_kl − W_l)); the identity's largest residual over all runs and records "
           f"{max(np.abs(i_identity(r)).max() for r in R.values()):.1e}.", "",
           "| sample | s | W_1 − 1/3 | the martingale term | the boundary term | the noise's drift −τ_kg′_k | the box's term | the cohesion's | "
           "the chain's bias |", "|---|---|---|---|---|---|---|---|---|"]
    for lab, (Rr, acc) in smp.items():
        for t in T_I:
            i = int(np.argmin(np.abs(s - t))); c = i_contrib(acc[:, i])
            m0, e0 = _ms(Rr[:, i, 4:5] - 1 / 3)
            L_.append(f"| {lab} | {s[i]:g} | {m0[0]:+.4f} ± {e0[0]:.4f} | " + ' | '.join(
                f"{c[:, j].mean():+.4f} ± {c[:, j].std(ddof=1) / np.sqrt(len(c)):.4f}" if j < 4 else f"{c[:, j].mean():+.1e}" for j in (1, 0, 2, 3, 4, 5)) + " |")
    # the martingale term regressed out
    L_ += ["", f"The first mass with the martingale term m (its contribution to W_1 above) regressed out: W_1 − b·m, b the sample's regression "
           "coefficient of W_1 on m — the same expectation as W_1, m's being zero, up to b's O(1/n) estimation; its s.e. from the residuals "
           f"(n − 2); samples of at least {len(SEEDS_G)} seeds — against Katzenberger's SDE (reflected at W_B, the step 2, for the loops with "
           "the box; without walls, the step 10, for the loop without it); z over that s.e. and the SDE's Monte Carlo error (std/√paths) combined; the "
           "displacements' ratio (W_1 − b·m − 1/3)/(the SDE's W_1 − 1/3).", "",
           "| sample | s | corr(W_1, m) | b | W_1 − b·m (± s.e.) | the SDE (± its error) | z | the displacements' ratio (± s.e.) |",
           "|---|---|---|---|---|---|---|---|"]
    for lab, (Rr, acc) in smp.items():
        if len(Rr) < len(SEEDS_G):
            continue
        ref = S['no walls, step 10' if 'no box' in lab else 'reflected at W_B, step 2']
        for t in T_I:
            i = int(np.argmin(np.abs(s - t))); m, e, c, b_, sm, ss = i_regressed(Rr[:, i, 4], i_contrib(acc[:, i])[:, 1], ref['rec'][t][:, 0])
            L_.append(f"| {lab} | {s[i]:g} | {c:+.3f} | {b_:.3f} | {m:.4f} ± {e:.4f} | {sm:.4f} ± {ss:.4f} | {(m - sm) / np.hypot(e, ss):+.1f} | "
                      f"{(m - 1 / 3) / (sm - 1 / 3):.2f} ± {e / abs(sm - 1 / 3):.2f} |")
    # paired differences
    tr = I['transient']
    L_ += ["", "Paired differences (the same seeds and draws) against the models' differences:", "",
           "| s | the box − no box: W_1 (± s.e.) | the SDE: reflected at W_B (step 2) − no walls | with the box measured (step 2) − no walls | "
           "the equilibrated start − the exact fit: W_1 (± s.e.) |", "|---|---|---|---|---|"]
    lab3, labnb, labeq = list(smp)[2:5]
    for t in T_I:
        i = int(np.argmin(np.abs(s - t)))
        d1 = _ms(smp[lab3][0][:, i, 4:5] - smp[labnb][0][:, i, 4:5]); d2 = _ms(smp[labeq][0][:, i, 4:5] - smp[lab3][0][:, i, 4:5])
        nw = S['no walls, step 10']['rec'][t][:, 0].mean()
        L_.append(f"| {s[i]:g} | {d1[0][0]:+.4f} ± {d1[1][0]:.4f} | {S['reflected at W_B, step 2']['rec'][t][:, 0].mean() - nw:+.4f} | "
                  f"{S[soft]['rec'][t][:, 0].mean() - nw:+.4f} | {d2[0][0]:+.4f} ± {d2[1][0]:.4f} |")
    L_ += ["", f"The positions' spreading from the exact fit as a mass offset (referee round 147's estimate: γ(E[φ_k] − φ_k(fit))/β_k at the "
           f"masses 1/3, the Gaussian law's mean shift and covariance): W = 1/3 + ({', '.join(f'{x:+.1e}' for x in tr)}).", ""]
    fx = I['fixed']; nW = len(fx['w1']); bx = fx['box'].reshape(nW, fx['ns'], 3); gp = fx['gp'].reshape(nW, fx['ns'], 3)
    ot = fx['out'].reshape(nW, fx['ns'], 3)
    L_ += [f"The box's term at fixed masses (W_2 = W_3 = (1 − W_1)/2; γ = 0, the box, 'correlated E'; {fx['ns']} runs each, {BURN_FIX_I:g} "
           f"slow units of burn-in, then {S_FIX_I:g}): class 1's time averages.", "",
           "| W_1 | E[g_1∂_{α_1}V] (± s.e.) | E[−τ_1g′_1] | their ratio | α_1 outside the box (share of steps) |", "|---|---|---|---|---|"]
    for j, w in enumerate(fx['w1']):
        mb, eb = bx[j, :, 0].mean(), bx[j, :, 0].std(ddof=1) / np.sqrt(fx['ns'])
        L_.append(f"| {w:.4f} | {mb:.3e} ± {eb:.1e} | {gp[j, :, 0].mean():.3e} | {mb / -gp[j, :, 0].mean():.3f} | {ot[j, :, 0].mean():.4f} |")
    L_.append("")
    return L_


def i_regressed(w, m, sde):
    """w (runs) regressed on the zero-mean m: (the mean of w − b·m, its s.e. (residuals, n − 2), corr(w, m), b, the SDE's mean and Monte
    Carlo error)."""
    b = np.cov(w, m)[0, 1] / m.var(ddof=1); r = w - b * m; n = len(w)
    return (r.mean(), r.std(ddof=2) / np.sqrt(n), np.corrcoef(w, m)[0, 1], b, sde.mean(), sde.std(ddof=1) / np.sqrt(len(sde)))


def i_contrib(acc):
    """Runs × (ACC_I) × 3 accumulators → each term's contribution to W_1 (runs × len(ACC_I)), linearised at 1/3; the martingale's sign
    as it enters u (−)."""
    d = GAMMA_E * (acc - acc[:, :, 2:3]); sg = np.array([1, -1, 1, 1, 1, 1])[None, :]
    return sg * np.einsum('l,rjl->rj', JAC13[0], d)


def i_identity(r):
    """The decomposition's residual: Δlog(β_k/β_3) − γΣ(terms), all runs and records."""
    x = np.log(r['W']); x = x - x[:, :, 2:3]
    sg = np.array([1, -1, 1, 1, 1, 1])[None, None, :, None]
    t = GAMMA_E * (sg * r['acc']).sum(axis=2); t = t - t[:, :, 2:3]
    return x - t


def main():
    raw = os.path.join(OUT, 'gauge_noise_raw.npz')
    md = os.path.join(OUT, 'gauge_noise.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    extra = lambda src: {k: src.get(k) for k in ('E_rec', 'F', 'G', 'H', 'I', 'J', 'K', 'M', 'N')}
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDE'], **extra(old)))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFGHIJKMN')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    names = list(noises())
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'A' in parts:
            h['A'] = pool.map_async(a_job, [(n, w) for n in names for w in W_A], chunksize=1)
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(n, w) for n in names for w in W_B], chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(g, w) for g in GAMMAS_C for w in W_B], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(d_job, list(RUNS_D), chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(e_series, [(v, s) for v in VARIANTS_E for s in SEEDS_E], chunksize=1)
        for k in 'ABCD':
            if k in h:
                res[k] = h[k].get(); E.save(raw, **res); log(k)
        if 'E' in h:
            runs = h['E'].get()
            res['E_rec'] = {v: np.array([r['rec'] for r in runs if r['variant'] == v]) for v in VARIANTS_E}
            res['E'] = [dict(variant=v, sum=summarise_E(res['E_rec'][v])) for v in VARIANTS_E]
            E.save(raw, **res); log('E')
        if 'G' in parts:                                                     # referee round 147: the early phase
            h['G'] = pool.map_async(g_series, [(v, s, S_G) for v in VARIANTS_G for s in SEEDS_G]
                                    + [(VARIANTS_G[0], s, S_G) for s in REPLAY_G], chunksize=1)
        if 'H' in parts:                                                     # coloured noise in simulation
            h['H'] = pool.map_async(h_mc, [(tc, SEED_H, S_H, N_BATCH_H) for tc in (0.0, TC_H)], chunksize=1)
        if 'F' in parts:                                                     # Part E's predictions (in this process)
            res['F'] = part_F(res['E_rec']); E.save(raw, **res); log('F')
        if 'G' in h:
            runs = h['G'].get(); new = runs[:len(VARIANTS_G) * len(SEEDS_G)]
            rp = [r['rec'] for r in runs[len(new):]]
            Er = np.asarray(res['E_rec'][VARIANTS_G[0]])
            res['G'] = dict(rec={v: np.array([r['rec'] for r in new if r['variant'] == v]) for v in VARIANTS_G}, replay=rp,
                            replay_diff=float(max(np.abs(r - Er[SEEDS_E.index(s), :len(r)]).max() for r, s in zip(rp, REPLAY_G))))
            E.save(raw, **res); log('G')
        if 'H' in h:
            res['H'] = part_H(h['H'].get()); E.save(raw, **res); log('H')
    if 'I' in parts:                                                         # the early phase's gap: one process (batched runs)
        res['I'] = part_I(res['F']); E.save(raw, **res); log('I')
    if 'J' in parts:                                                         # referee round 153: one process each
        res['J'] = j_runs((S_J, SEED_J)); E.save(raw, **res); log('J')
    if 'K' in parts:
        res['K'] = part_K(res['E_rec'], res.get('K')); E.save(raw, **res); log('K')
    if 'M' in parts:                                                         # the missing O(γ²τ) drift: one process, resumable
        res['M'] = part_M(res, lambda: E.save(raw, **res)); E.save(raw, **res); log('M')
    if 'N' in parts:                                                         # the finite-τ question: one process, resumable
        res['N'] = part_N(res, lambda: E.save(raw, **res)); E.save(raw, **res); log('N')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDE'], **extra(res)))


if __name__ == '__main__':
    main()
