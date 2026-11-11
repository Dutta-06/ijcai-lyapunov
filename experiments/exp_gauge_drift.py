"""
Proposition D‴.37: the dither's drift along the gauge — in the slow-replicator limit the masses descend the log-determinant of the
positions' Hessian, toward the saturation, and stop in a boundary layer of width O(√γ)
=================================================================================================================================

Proposition D‴.36's instance (Proposition D‴.32's three bumps, one particle per class at an exact representation, the replicator at γ, the
positions dithered at τ = 10⁻⁴): the drift of Katzenberger's limit on the manifold M of exact fits, from the second-order identities of the
limit map (`src/core/gauge_drift.py`), in the masses' coordinates (W_1, W_2), W_3 = 1 − W_1 − W_2.

A. The identities against the limit map computed by integration (LSODA to 1.2·10⁵ slow units) and its second differences along the nine
   noise coordinates, at W = 1/3 and at (0.30, 0.35, 0.35), γ = 10⁻².
B. The drift split by the noise coordinates (each class's amplitude, centre, width) at the same points.
C. The slow-replicator limit: the drift over γτ for γ = 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵ against the replicator flow of the entropic free energy,
   −½W_k(∂_kℓ − Σ_lW_l∂_lℓ), ℓ(β) = log det ∇²_ϑJ at the exact representation with the masses β/3 (finite differences).
D. The saturation layer: the zeros of the first class's drift on the line W_2 = W_3 for γ = 10⁻²–10⁻⁴, their distance from the first
   class's saturation |A_1|/a_max = 0.25, against √γ and the rates' crossing δ*; deep inside, the drift against the amplitude's free
   diffusion carried along M, τG_α∂²_αW_1.
E. The drift field at γ = 10⁻² on a grid of the representable masses (W_k ≥ |A_k|/a_max + 0.004); its zeros (Newton from the grid's
   twelve smallest drifts and from the three points 0.035 inside two saturations) and the eigenvalues of its Jacobian there.
F. The ensemble mean of exp_gauge_dither.py's Part B (64 seeds, 5·10⁴ slow units, from its raw file) against the mean path dW/ds = b(W).
G. The entropic leak's pinned state (ε = 10⁻⁴, 3·10⁻⁴, 10⁻³; 16 seeds each, 4·10⁵ slow units, the means over the second half) against
   the shift κ⁻¹b(w) (κ the leak's stiffness along M, Proposition D‴.36(iii)) and the nonlinear fixed point b(W) = κ(W)(W − w).
H. Without a leak over 10⁶ slow units (12 seeds): the masses' means and spreads at 10⁵, 3·10⁵, 10⁶ against the limit SDE on M
   (Euler–Maruyama, 4000 paths, ds = 50, the drift and the diffusion interpolated linearly on E's grid, steps leaving the grid refused).
I. The time scales (referee round 143, D3el-1, -4, -5; no new runs of the loop): the slowest normal rates of the linearisation at W = 1/3
   (γ = 10⁻³, 10⁻²) and the slowest mode's share of the masses in its g-norm; the leak's stiffness along M at G's three ε against them, and
   the leak's linearised spread over the distance 1/3 − |A_k|/a_max to the faces; the normal rates along W_2 = W_3 down to 5·10⁻⁴ from the
   first class's face (γ = 10⁻², 10⁻⁴); the Péclet number |b_1|(1/3 − |A_1|/a_max)/D_11 at W = 1/3 and the free spread √(D_11/(γτ)) at the
   drift's time 1/(γτ) (γ = 10⁻², 10⁻³, 10⁻⁴); and F's 64-seed ensemble against the limit SDE's ensemble (as in H, 4000 paths, its own seed)
   rather than against the mean path, with the ensembles' mean drift over [0, s].

Output: results/chaos/gauge_drift.md and gauge_drift_raw.npz (GDR_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`; `--only E,S`
recomputes H's limit SDE on a new field from the stored full-loop records; `--only I` computes I from the stored E and F).
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
import exp_cohesion_leak as E
import exp_gauge_dither as GD

OUT = os.environ.get('GDR_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('GDR_WORKERS', 3))
GAMMA, TAU, A_MAX = GD.GAMMA, GD.TAU, 4.0
W_A = ((1 / 3, 1 / 3), (0.30, 0.35))
H_BF, S_BF = 0.02, 1.2e5
GAMMAS_C, W_C = (1e-2, 1e-3, 1e-4, 1e-5), ((1 / 3, 1 / 3), (0.30, 0.35), (0.40, 0.25), (0.28, 0.36))
GAMMAS_D = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4)
GRID_N, GRID_PAD = 24, 0.004
EPS_G, SEEDS_G, S_G = (1e-4, 3e-4, 1e-3), tuple(range(200, 216)), 4e5
SEEDS_H, S_H = tuple(range(300, 312)), 1e6
T_H = (1e5, 3e5, 1e6)
N_PATHS, DS_SDE = 4000, 50.0
GAMMAS_I, DELTAS_I, GAMMAS_IF, GAMMAS_PE = (1e-3, 1e-2), (0.08, 0.02, 0.005, 0.002, 0.0005), (1e-2, 1e-4), (1e-2, 1e-3, 1e-4)
T_I, SEED_I = (1e4, 2e4, 3e4, 5e4), 1
W_MIN = np.abs(E.A_BUMP) / A_MAX
H_ID, H2_ID = 3e-6, 1e-4                                                       # the identities' steps (near the saturation the chart curves on 10⁻³)
NAMES = [f'{q}_{k + 1}' for k in range(3) for q in ('α', 'c', 's')]


def instance(gamma=GAMMA):
    """The flow in (Θ, log(β_1/β_3), log(β_2/β_3)), the chart of M by (W_1, W_2), the coordinates, the noise covariance."""
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=gamma)

    def beta_of(z):
        e = np.exp(np.r_[z, 0.0]); return 3 * e / e.sum()

    def F(x):
        T = x[:9].reshape(3, 3); b = beta_of(x[9:]); dT, db = loop.drift(T, b); dl = db / b
        return np.concatenate([dT.ravel(), dl[:2] - dl[2]])

    def chart(w):
        W = np.r_[w, 1 - w.sum()]; return np.concatenate([E.exact_rep(net, W).ravel(), np.log(W[:2] / W[2])])

    coords = lambda x: beta_of(x[9:])[:2] / 3
    S = np.zeros((11, 11)); S[np.arange(9), np.arange(9)] = 2 * TAU * np.tile(E.G, 3)
    return F, chart, coords, S


def drift_at(w, gamma=GAMMA, by_noise=False):
    F, chart, coords, S = instance(gamma)
    return GDR.limit_drift(F, chart, coords, np.asarray(w, float), S, h=H_ID, h2=H2_ID, by_noise=by_noise)


def log_det(beta, h=1e-4):
    """ℓ(β) = log det ∇²_ϑJ at the exact representation with the masses β/3 (the nine position coordinates; central differences)."""
    net, base = E.setup()
    loop = MT.ReplicatorCohesionLoop(base, 1.0, 'compact', gamma=GAMMA)
    Th = E.exact_rep(net, np.asarray(beta) / 3).ravel()
    J = lambda t: loop.J(t.reshape(3, 3), np.asarray(beta, float))
    H = np.zeros((9, 9)); I = np.eye(9) * h
    for i in range(9):
        for j in range(i, 9):
            H[i, j] = H[j, i] = (J(Th + I[i] + I[j]) - J(Th + I[i] - I[j]) - J(Th - I[i] + I[j]) + J(Th - I[i] - I[j])) / (4 * h * h)
    sgn, ld = np.linalg.slogdet(H)
    assert sgn > 0
    return ld


def log_det_law(W, hb=1e-3):
    """The slow-replicator limit of the drift over γτ: −½W_k(∂_kℓ − Σ_lW_l∂_lℓ), ∂_k = ∂/∂β_k (the replicator flow of (τ/2)ℓ)."""
    W = np.asarray(W, float); beta = 3 * W
    g = np.array([(log_det(beta + hb * e) - log_det(beta - hb * e)) / (2 * hb) for e in np.eye(3)])
    return -0.5 * W * (g - W @ g)


def bf_job(args):
    """½∂²(W∘Φ)[Σ] by second differences of the integrated limit map along one noise coordinate (± h)."""
    w, i = args
    F, chart, coords, S = instance()
    x0 = chart(np.asarray(w, float)); e = np.zeros(11); e[i] = H_BF
    yp, rp = GDR.limit_map(F, x0 + e, S_BF); ym, rm = GDR.limit_map(F, x0 - e, S_BF)
    return dict(w=w, i=i, val=0.5 * S[i, i] * (coords(yp) - 2 * np.asarray(w) + coords(ym)) / H_BF ** 2, res=max(rp, rm))


def c_job(args):
    w, gamma = args
    d = drift_at(w, gamma)[0]
    return dict(w=w, gamma=gamma, d=np.r_[d, -d.sum()] / (gamma * TAU))


def d_job(gamma):
    from scipy.optimize import brentq
    b1 = lambda W1: drift_at(np.array([W1, (1 - W1) / 2]), gamma)[0][0]
    xs = W_MIN[0] + np.geomspace(2e-3, 0.08, 22)
    v = [b1(x) for x in xs]
    roots = [brentq(b1, xs[i], xs[i + 1], xtol=1e-7) for i in range(len(xs) - 1) if np.sign(v[i]) != np.sign(v[i + 1])]
    return dict(gamma=gamma, roots=np.array(roots), signs=np.sign(v), xs=xs, b=np.array(v), ito=np.array([ito_inside(x) for x in xs]))


def ito_inside(W1):
    """Inside the layer: τG_α d²W_1/dα_1² for W_1 = A_1/a(α_1), a = a_max tanh(α/a_max) — the amplitude's free diffusion carried along M."""
    A = abs(E.A_BUMP[0]); a = A / W1; ap = 1 - (a / A_MAX) ** 2
    return TAU * E.G[0] * A * (2 * ap ** 2 / a ** 3 + 2 * ap / (a * A_MAX ** 2))


def grid_points():
    pts = []
    for W1 in np.linspace(W_MIN[0] + GRID_PAD, 1 - W_MIN[1] - W_MIN[2] - 2 * GRID_PAD, GRID_N):
        for W2 in np.linspace(W_MIN[1] + GRID_PAD, 1 - W_MIN[0] - W_MIN[2] - 2 * GRID_PAD, GRID_N):
            if 1 - W1 - W2 >= W_MIN[2] + GRID_PAD - 1e-12:
                pts.append((W1, W2))
    return pts


def grid_job(w):
    d, D, ev = drift_at(w)
    return np.r_[w, d, D[0, 0], D[0, 1], D[1, 1], ev.max()]


def zeros(field):
    from scipy.optimize import fsolve

    def safe(w):                                                             # Newton's trial points outside the representable masses
        if not np.all(np.r_[w, 1 - np.sum(w)] > W_MIN + 1e-3):
            return np.full(2, 1e3)
        try:
            return drift_at(w)[0] * 1e7
        except ValueError:
            return np.full(2, 1e3)
    out = []
    W3 = np.c_[field[:, :2], 1 - field[:, 0] - field[:, 1]] - W_MIN               # distances from the saturations
    corners = []                                                             # two classes 0.035 inside their saturation, the third the rest
    for i, j in ((0, 1), (0, 2), (1, 2)):
        c = np.zeros(3); c[i] = W_MIN[i] + 0.035; c[j] = W_MIN[j] + 0.035; c[3 - i - j] = 1 - c[i] - c[j]
        corners.append(c[:2])
    for w0 in list(field[np.argsort(np.hypot(field[:, 2], field[:, 3]))[:12], :2]) + corners:
        z, info, ier, msg = fsolve(safe, w0, full_output=True)
        if ier == 1 and np.all(np.r_[z, 1 - z.sum()] > W_MIN + 1e-3) and not any(np.hypot(*(z - o['w'])) < 1e-4 for o in out):
            h = 1e-4
            Jb = np.array([(drift_at(z + h * e)[0] - drift_at(z - h * e)[0]) / (2 * h) for e in np.eye(2)]).T
            out.append(dict(w=z, eig=np.linalg.eigvals(Jb), D=drift_at(z)[1]))
    return out


def e_job(args):
    return GD.series(args)


def kappa(W, eps):
    """The leak's stiffness along M in the coordinates (W_1, W_2) (Proposition D‴.36(iii))."""
    A = E.A_BUMP; K = 3; m = 3; Ga = E.G[0]; W = np.asarray(W, float)
    a = A / W; ap = 1 - (a / A_MAX) ** 2
    T = np.zeros((2 * K, K - 1))
    for j in range(K - 1):
        dW = np.zeros(K); dW[j] = 1.0; dW[-1] = -1.0
        T[:K, j] = -(A / (W ** 2 * ap)) * dW; T[K:, j] = dW
    g = np.diag(np.r_[np.full(K, 1 / Ga), (m / GAMMA) / W])
    return np.linalg.solve(T.T @ g @ T, T.T @ np.diag(np.r_[np.zeros(K), eps * m / (GAMMA * W)]) @ T)


def part_F(field):
    z = E.load(os.path.join(OUT, 'gauge_dither_raw.npz'))
    R = np.array([r['rec'] for r in z['B']]); s = R[0, :, 0]; W = R[:, :, 4:7]
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda t, w: drift_at(w)[0], (0, s[-1]), [1 / 3, 1 / 3], t_eval=[1e4, 2e4, 3e4, 5e4], rtol=1e-6, atol=1e-9)
    out = dict(n=len(W), s=[], mean=[], se=[], path=[], spread=[])
    Dk = np.diag(GD.gauge_diffusion(np.ones(3) / 3, GAMMA, TAU))
    for j, t in enumerate(sol.t):
        i = int(np.argmin(np.abs(s - t)))
        out['s'].append(s[i]); out['mean'].append(W[:, i].mean(axis=0)); out['se'].append(W[:, i].std(axis=0, ddof=1) / np.sqrt(len(W)))
        p = sol.y[:, j]; out['path'].append(np.r_[p, 1 - p.sum()]); out['spread'].append(np.sqrt(Dk * s[i]))
    return out


def part_G(runs):
    from scipy.optimize import fsolve
    w = np.ones(3) / 3
    b = drift_at(w[:2])[0]
    out = []
    for eps in EPS_G:
        rs = [r for r in runs if r['eps'] == eps]
        Rr = np.array([r['rec'] for r in rs]); h = Rr.shape[1] // 2
        means = Rr[:, h:, 4:7].mean(axis=1)
        c = np.linalg.solve(kappa(w, eps), b)
        f = lambda c_: (drift_at(w[:2] + c_)[0] - kappa(np.r_[w[:2] + c_, 1 - (w[:2] + c_).sum()], eps) @ c_) * 1e7
        cn = fsolve(f, c)
        out.append(dict(eps=eps, shift=means.mean(axis=0) - 1 / 3, se=means.std(axis=0, ddof=1) / np.sqrt(len(means)), lin=np.r_[c, -c.sum()],
                        nl=np.r_[cn, -cn.sum()], std=Rr[:, h:, 4:7].std(axis=(0, 1)), std_pred=GD.leak_spread(w, GAMMA, TAU, eps)))
    return out


def sde_paths(field, S_end, t_rec, seed=0):
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import Delaunay
    P = field[:, :2]
    ib = LinearNDInterpolator(P, field[:, 2:4]); iD = LinearNDInterpolator(P, field[:, 4:7]); tri = Delaunay(P)
    rng = np.random.default_rng(seed)
    W = np.tile([1 / 3, 1 / 3], (N_PATHS, 1))
    rec = {}
    n = int(round(S_end / DS_SDE))
    for step in range(n + 1):
        s = step * DS_SDE
        for t in t_rec:
            if abs(s - t) < DS_SDE / 2:
                rec[t] = W.copy()
        if step == n:
            break
        b = ib(W); D = iD(W)
        L11 = np.sqrt(np.maximum(D[:, 0], 0)); L21 = D[:, 1] / np.where(L11 > 0, L11, 1); L22 = np.sqrt(np.maximum(D[:, 2] - L21 ** 2, 0))
        z = rng.standard_normal((N_PATHS, 2)) * np.sqrt(DS_SDE)
        Wn = W + b * DS_SDE + np.c_[L11 * z[:, 0], L21 * z[:, 0] + L22 * z[:, 1]]
        out = tri.find_simplex(Wn) < 0
        Wn[out] = W[out]                                                     # a step that leaves the grid is refused
        W = Wn
    return {t: np.c_[v, 1 - v.sum(axis=1)] for t, v in rec.items()}


def part_H(R, field):
    """R: the full loop's records (seeds × records × (s, α_1..3, W_1..3, fit))."""
    sde = sde_paths(field, S_H, T_H)
    s = R[0, :, 0]; W = R[:, :, 4:7]; al = R[:, :, 1:4]
    out = []
    for t in T_H:
        i = int(np.argmin(np.abs(s - t)))
        out.append(dict(s=s[i], mean=W[:, i].mean(axis=0), se=W[:, i].std(axis=0, ddof=1) / np.sqrt(len(W)), std=W[:, i].std(axis=0, ddof=1),
                        sde_mean=sde[t].mean(axis=0), sde_std=sde[t].std(axis=0, ddof=1),
                        sat=np.abs(np.tanh(al[:, i] / A_MAX)).max(axis=0), near=np.mean(W[:, i] - W_MIN < 0.02, axis=0),
                        sde_near=np.mean(sde[t] - W_MIN < 0.02, axis=0)))
    q = np.array_split(np.arange(R.shape[1]), 4)                             # round 27 (D3z-1): the excursions into the saturation by quarters
    out.append(dict(quarters=np.array([np.abs(al[:, ix]).max(axis=1).mean(axis=0) for ix in q]),
                    big=np.array([np.mean(np.abs(al[:, ix]) > 10, axis=(0, 1)) for ix in q]),
                    below=np.mean(W[:, q[-1]] < W_MIN, axis=(0, 1))))
    return out


def normal_rates(w, gamma, n=4):
    """The n slowest normal rates of the linearisation at the exact fit with the masses (w, 1 − Σw) (the eigenvalues of −∂F off the tangent),
    and the slowest normal mode's shares of its g-norm: the masses' (Σδβ_i²/(γβ_i)) and the positions' by coordinate (α, c, s)."""
    from src.core.gauge_drift import _jac
    F, chart, coords, S = instance(gamma)
    w = np.asarray(w, float)
    A = _jac(F, chart(w), H_ID)
    lam, V = np.linalg.eig(A)
    o = np.argsort(np.abs(lam)); lam, V = lam[o], V[:, o]
    v = V[:, 2].real
    W = np.r_[w, 1 - w.sum()]; beta = 3 * W
    dz = np.r_[v[9:], 0.0]; db = beta * (dz - W @ dz)                         # δβ from δ log(β_k/β_3)
    pos = (v[:9].reshape(3, 3) ** 2 / E.G[None, :]).sum(axis=0); mas = np.sum(db ** 2 / (gamma * beta))
    return dict(w=W, gamma=gamma, rates=-lam[2:2 + n].real, mass=mas / (pos.sum() + mas), pos=pos / (pos.sum() + mas))


def part_I(field, Ff):
    w = np.ones(3) / 3; dist = w - W_MIN
    out = dict(rates=[normal_rates(w[:2], g) for g in GAMMAS_I], stiff=[], face=[], pe=[], sde=[])
    for eps in EPS_G:
        sp = GD.leak_spread(w, GAMMA, TAU, eps)
        out['stiff'].append(dict(eps=eps, eig=np.sort(np.linalg.eigvals(kappa(w, eps)).real), spread=sp, rel=sp / dist))
    for g in GAMMAS_IF:
        for d in DELTAS_I:
            W1 = W_MIN[0] + d; wf = np.array([W1, (1 - W1) / 2])
            r = normal_rates(wf, g); r['delta'] = d; r['alpha1'] = float(instance(g)[1](wf)[0])
            out['face'].append(r)
    for g in GAMMAS_PE:
        b, D, _ = drift_at(w[:2], g)
        out['pe'].append(dict(gamma=g, b1=b[0], D11=D[0, 0], pe=abs(b[0]) * dist[0] / D[0, 0], spread=np.sqrt(D[0, 0] / (g * TAU))))
    if Ff is not None and field is not None:
        sde = sde_paths(field, max(T_I), T_I, seed=SEED_I)
        for s_, m_, se, p_ in zip(Ff['s'], Ff['mean'], Ff['se'], Ff['path']):
            t = min(T_I, key=lambda u: abs(u - s_)); Y = sde[t]
            ys = Y.std(axis=0, ddof=1) / np.sqrt(len(Y))
            out['sde'].append(dict(s=s_, loop=m_, loop_se=se, sde=Y.mean(axis=0), sde_se=ys, sde_std=Y.std(axis=0, ddof=1), path=p_,
                                   z=(m_ - Y.mean(axis=0)) / np.hypot(se, ys), z_path=(m_ - p_) / se,
                                   rate_loop=(m_ - 1 / 3) / (s_ * GAMMA * TAU), rate_loop_se=se / (s_ * GAMMA * TAU),
                                   rate_sde=(Y.mean(axis=0) - 1 / 3) / (s_ * GAMMA * TAU), rate_sde_se=ys / (s_ * GAMMA * TAU),
                                   rate_path=(p_ - 1 / 3) / (s_ * GAMMA * TAU)))
        b0 = drift_at(w[:2])[0]; out['b0'] = np.r_[b0, -b0.sum()] / (GAMMA * TAU)
    return out


def report_I(Ii):
    fmt = E.fmt
    w3 = lambda v, p=4: ', '.join(f'{x:.{p}f}' for x in v)
    L_ = ["## I. The time scales: the normal rates, the leak's window, the Péclet number, and the limit SDE's ensemble", "",
          "The normal rates: the eigenvalues of −∂F off the tangent of M (per unit slow time); the slowest normal mode's shares of its g-norm "
          "(g = G⁻¹ on the positions, Σδβ²/(γβ) on the masses).", "",
          "| W | γ | the four slowest normal rates | the slowest mode's share: masses | positions α, c, s |", "|---|---|---|---|---|"]
    for r in Ii['rates']:
        L_.append(f"| {w3(r['w'], 3)} | {r['gamma']:g} | {', '.join(fmt(v) for v in r['rates'])} | {fmt(r['mass'])} | {', '.join(fmt(v) for v in r['pos'])} |")
    L_ += ["", f"The leak's stiffness along M at W = 1/3 (γ = {GAMMA:g}; the eigenvalues of κ = ε(TᵀgT)⁻¹H_W, per unit slow time) against those rates, "
           "and the leak's linearised spread (Proposition D‴.36(iii)) over the distance 1/3 − |A_k|/a_max to the faces:", "",
           "| ε | the stiffness's eigenvalues | the spread (std of W_1, W_2, W_3) | over the distance to the faces |", "|---|---|---|---|"]
    for r in Ii['stiff']:
        L_.append(f"| {r['eps']:g} | {', '.join(fmt(v) for v in r['eig'])} | {w3(r['spread'])} | {w3(r['rel'], 3)} |")
    L_ += ["", "The normal rates along W_2 = W_3 toward the first class's face (α_1 the exact fit's amplitude coordinate):", "",
           "| γ | W_1 − 0.25 | α_1 | the three slowest normal rates |", "|---|---|---|---|"]
    for r in Ii['face']:
        L_.append(f"| {r['gamma']:g} | {r['delta']:g} | {r['alpha1']:.2f} | {', '.join(fmt(v) for v in r['rates'][:3])} |")
    L_ += ["", "The Péclet number at W = 1/3 over the distance to the first face, Pe = |b_1|(1/3 − 0.25)/D_11, and the free spread at the drift's time "
           "1/(γτ), √(D_11/(γτ)) (the distance to the face 0.0833):", "",
           "| γ | b_1 | D_11 | Pe | √(D_11/(γτ)) |", "|---|---|---|---|---|"]
    for r in Ii['pe']:
        L_.append(f"| {r['gamma']:g} | {fmt(r['b1'])} | {fmt(r['D11'])} | {r['pe']:.2f} | {r['spread']:.3f} |")
    if Ii['sde']:
        L_ += ["", f"F's ensemble (64 seeds, the full loop) against the limit SDE's ensemble ({N_PATHS} paths from 1/3, seed {SEED_I}, as in H) and "
               "against the mean path; the ensembles' mean drift over [0, s], (mean(s) − 1/3)/(sγτ), against Katzenberger's drift at 1/3 over γτ, "
               + ', '.join(f'{v:+.4f}' for v in Ii['b0']) + ":", "",
               "| s | full loop: mean (± s.e.) | limit SDE: mean (± s.e.) | its std | the distance in s.e. (combined) | the mean path | the distance in s.e. | "
               "mean drift over [0, s] / γτ: full loop (± s.e.) | limit SDE (± s.e.) | the mean path |", "|---|---|---|---|---|---|---|---|---|---|"]
        for r in Ii['sde']:
            L_.append(f"| {r['s']:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(r['loop'], r['loop_se']))} | "
                      f"{', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(r['sde'], r['sde_se']))} | {w3(r['sde_std'])} | "
                      f"{', '.join(f'{v:+.2f}' for v in r['z'])} | {w3(r['path'])} | {', '.join(f'{v:+.2f}' for v in r['z_path'])} | "
                      f"{', '.join(f'{a:+.3f} ± {b:.3f}' for a, b in zip(r['rate_loop'], r['rate_loop_se']))} | "
                      f"{', '.join(f'{a:+.3f} ± {b:.3f}' for a, b in zip(r['rate_sde'], r['rate_sde_se']))} | {', '.join(f'{v:+.3f}' for v in r['rate_path'])} |")
    L_.append("")
    return L_


def report(A=None, B=None, C=None, D=None, Ef=None, Ff=None, Gg=None, Hh=None, Ii=None):
    fmt = E.fmt
    w3 = lambda w: ', '.join(f'{v:.3f}' for v in np.r_[w, 1 - sum(w)])
    L_ = ["# Proposition D‴.37: the dither's drift along the gauge (`exp_gauge_drift.py`)", "",
          f"Proposition D‴.36's instance (three bumps, one particle per class at an exact representation, m = 3, the replicator at γ, "
          f"the positions dithered at τ = {TAU:g}); the masses W_k = β_k/3, the first class saturating at W_1 = {W_MIN[0]:g}; drifts and diffusions "
          "per unit slow time; γ = 10⁻² unless stated.", ""]
    if A is not None:
        L_ += ["## A. The second-order identities against the integrated limit map", "",
               f"The limit map by LSODA to {S_BF:g} slow units, its second differences at ± {H_BF:g} along the nine noise coordinates.", "",
               "| W | drift of W_1, W_2, W_3: identities | integrated limit map | relative difference | the flow's residual at the end |", "|---|---|---|---|---|"]
        for r in A:
            d3 = np.r_[r['id'], -r['id'].sum()]; b3 = np.r_[r['bf'], -r['bf'].sum()]
            L_.append(f"| {w3(r['w'])} | {', '.join(fmt(v) for v in d3)} | {', '.join(fmt(v) for v in b3)} | "
                      f"{np.abs(d3 - b3).max() / np.abs(b3).max():.1e} | {fmt(r['res'])} |")
        L_.append("")
    if B is not None:
        L_ += ["## B. The drift split by the noise coordinates", "", "| W | noise on | drift of W_1, W_2, W_3 |", "|---|---|---|"]
        for r in B:
            for nm, p in zip(NAMES, r['parts']):
                L_.append(f"| {w3(r['w'])} | {nm} | {', '.join(fmt(v) for v in np.r_[p, -p.sum()])} |")
            L_.append(f"| {w3(r['w'])} | all | {', '.join(fmt(v) for v in np.r_[r['total'], -r['total'].sum()])} |")
        L_.append("")
    if C is not None:
        L_ += ["## C. The slow-replicator limit: the drift over γτ against the replicator flow of the entropic free energy (τ/2)ℓ", "",
               "| W | " + ' | '.join(f'γ = {g:g}' for g in GAMMAS_C) + " | −½W_k(∂_kℓ − Σ_lW_l∂_lℓ) | ℓ |", "|---|" + "---|" * (len(GAMMAS_C) + 2)]
        for r in C:
            L_.append(f"| {w3(r['w'])} | " + ' | '.join(', '.join(f'{v:+.3f}' for v in r['d'][g]) for g in GAMMAS_C)
                      + f" | {', '.join(f'{v:+.3f}' for v in r['law'])} | {r['ell']:.3f} |")
        L_.append("")
    if D is not None:
        L_ += [f"## D. The saturation layer: the zeros of the first class's drift on W_2 = W_3 (W_1 − {W_MIN[0]:g} from 0.002 to 0.08)", "",
               "δ* = ½√(γa_max|A_1|/(mG_α)): where the saturating amplitude's relaxation rate mG_αΓ_11W_1²a′_1² ≈ 4mG_αΓ_11(W_1 − 0.25)² meets "
               "the mass channel's γa_1²W_1Γ_11 ≈ γa_max|A_1|Γ_11.", "",
               "| γ | the zeros W_1° − 0.25 | over √γ | over δ* | the sign of b_1 from the saturation outward | b_1 / τG_α∂²_αW_1 at W_1 − 0.25 = 0.002, 0.0024, 0.0028 |", "|---|---|---|---|---|---|"]
        for r in D:
            sg = ''.join('+' if v > 0 else '−' for v in r['signs'])
            ds = 0.5 * np.sqrt(r['gamma'] * A_MAX * abs(E.A_BUMP[0]) / (3 * E.G[0]))
            L_.append(f"| {r['gamma']:g} | {', '.join(f'{v - W_MIN[0]:.5f}' for v in r['roots']) or '—'} | "
                      f"{', '.join(f'{(v - W_MIN[0]) / np.sqrt(r['gamma']):.3f}' for v in r['roots']) or '—'} | "
                      f"{', '.join(f'{(v - W_MIN[0]) / ds:.3f}' for v in r['roots']) or '—'} | {sg} | "
                      f"{', '.join(f'{b_ / i_:.3f}' for b_, i_ in zip(r['b'][:3], r['ito'][:3]))} |")
        L_.append("")
    if Ef is not None:
        f = Ef['field']
        L_ += [f"## E. The drift field on the representable masses ({len(f)} grid points, W_k ≥ |A_k|/a_max + {GRID_PAD:g})", "",
               f"Largest |drift|: {fmt(np.abs(f[:, 2:4]).max())}; D_11 from {fmt(f[:, 4].min())} to {fmt(f[:, 4].max())}; the slowest transverse rate from "
               f"{fmt(-f[:, 7].max())} to {fmt(-f[:, 7].min())}.", "",
               "| the zero W° of the drift | eigenvalues of the drift's Jacobian (in W_1, W_2) | type | D there (D_11, D_12, D_22) |", "|---|---|---|---|"]
        for z in Ef['zeros']:
            ev = np.sort(z['eig'].real); typ = 'attracting' if ev.max() < 0 else ('repelling' if ev.min() > 0 else 'saddle')
            L_.append(f"| {', '.join(f'{v:.4f}' for v in np.r_[z['w'], 1 - z['w'].sum()])} | {', '.join(fmt(v) for v in z['eig'])} | {typ} | "
                      f"{fmt(z['D'][0, 0])}, {fmt(z['D'][0, 1])}, {fmt(z['D'][1, 1])} |")
        L_.append("")
    if Ff is not None:
        L_ += [f"## F. The ensemble mean of exp_gauge_dither.py's Part B ({Ff['n']} seeds) against the mean path dW/ds = b(W)", "",
               "| s | ensemble mean of W_1, W_2, W_3 (± standard error) | the mean path | the free diffusion's spread √(D_kk s) |", "|---|---|---|---|"]
        for s_, m_, se, p_, sp in zip(Ff['s'], Ff['mean'], Ff['se'], Ff['path'], Ff['spread']):
            L_.append(f"| {s_:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(m_, se))} | {', '.join(f'{v:.4f}' for v in p_)} | {', '.join(f'{v:.3f}' for v in sp)} |")
        L_.append("")
    if Gg is not None:
        L_ += [f"## G. The entropic leak's pinned state ({len(SEEDS_G)} seeds each, {S_G:g} slow units, the second half)", "",
               "| ε | the masses' mean − 1/3 (± standard error over the seeds) | predicted κ⁻¹b(1/3) | the nonlinear fixed point − 1/3 | std (measured) | std (Proposition D‴.36(iii)) |",
               "|---|---|---|---|---|---|"]
        for r in Gg:
            L_.append(f"| {r['eps']:g} | {', '.join(f'{a:+.4f} ± {b:.4f}' for a, b in zip(r['shift'], r['se']))} | {', '.join(f'{v:+.4f}' for v in r['lin'])} | "
                      f"{', '.join(f'{v:+.4f}' for v in r['nl'])} | {', '.join(f'{v:.4f}' for v in r['std'])} | {', '.join(f'{v:.4f}' for v in r['std_pred'])} |")
        L_.append("")
    if Hh is not None:
        L_ += [f"## H. Without a leak: the full loop ({len(SEEDS_H)} seeds) against the limit SDE on M ({N_PATHS} paths)", "",
               "| s | full loop: mean of W_1, W_2, W_3 (± s.e.) | std | within 0.02 of the saturation (fraction of seeds) | the largest |a|/a_max | limit SDE: mean | std | within 0.02 |",
               "|---|---|---|---|---|---|---|---|"]
        for r in [r for r in Hh if 's' in r]:
            L_.append(f"| {r['s']:g} | {', '.join(f'{a:.4f} ± {b:.4f}' for a, b in zip(r['mean'], r['se']))} | {', '.join(f'{v:.4f}' for v in r['std'])} | "
                      f"{', '.join(f'{v:.2f}' for v in r['near'])} | {', '.join(f'{v:.4f}' for v in r['sat'])} | {', '.join(f'{v:.4f}' for v in r['sde_mean'])} | "
                      f"{', '.join(f'{v:.4f}' for v in r['sde_std'])} | {', '.join(f'{v:.2f}' for v in r['sde_near'])} |")
        ex = [r for r in Hh if 'quarters' in r]
        if ex:
            ex = ex[0]
            L_ += ["", "The full loop's excursions into the saturation, by quarters of the run (the seeds' mean of the largest |α_k| in the quarter; the "
                   "fraction of records with |α_k| > 10):", "", "| quarter | mean largest |α_1|, |α_2|, |α_3| | fraction |α_k| > 10 |", "|---|---|---|"]
            for i_, (qa, qb) in enumerate(zip(ex['quarters'], ex['big'])):
                L_.append(f"| {i_ + 1} | {', '.join(f'{v:.1f}' for v in qa)} | {', '.join(f'{v:.2f}' for v in qb)} |")
            L_ += ["", "In the last quarter the masses are below their saturation (off the manifold of exact fits, the fit then not exact) in "
                   + ', '.join(f'{v:.2f}' for v in ex['below']) + " of the records."]
        L_.append("")
    if Ii is not None:
        L_ += report_I(Ii)
    L_.append("Numbers are produced by `experiments/exp_gauge_drift.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'gauge_drift_raw.npz')
    md = os.path.join(OUT, 'gauge_drift.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEFGHI']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFGHI')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}                                                               # the cheap parts first: the pool's queue is FIFO
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(w, g) for w in W_C for g in GAMMAS_C], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(d_job, GAMMAS_D, chunksize=1)
        if parts & set('EFHS'):
            h['E'] = pool.map_async(grid_job, grid_points(), chunksize=4)
        if 'A' in parts:
            h['A'] = pool.map_async(bf_job, [(w, i) for w in W_A for i in range(9)], chunksize=1)
        if 'H' in parts:
            h['H'] = pool.map_async(e_job, [(GAMMA, 0.0, s, S_H) for s in SEEDS_H], chunksize=1)
        if 'G' in parts:
            h['G'] = pool.map_async(e_job, [(GAMMA, e, s, S_G) for e in EPS_G for s in SEEDS_G], chunksize=1)
        if 'B' in parts:
            res['B'] = [dict(w=np.array(w), total=d, parts=p) for w in W_A for d, _, _, p in [drift_at(w, by_noise=True)]]
            E.save(raw, **res); log('B')
        if 'C' in parts:
            cr = h['C'].get()
            res['C'] = [dict(w=np.array(w), d={r['gamma']: r['d'] for r in cr if r['w'] == w}, law=log_det_law(np.r_[w, 1 - sum(w)]),
                             ell=log_det(3 * np.r_[w, 1 - sum(w)])) for w in W_C]
            E.save(raw, **res); log('C')
        if 'D' in parts:
            res['D'] = h['D'].get(); E.save(raw, **res); log('D')
        if 'E' in h:
            field = np.array(h['E'].get())
            res['E'] = dict(field=field, zeros=zeros(field)); E.save(raw, **res); log('E')
        else:
            field = res['E']['field'] if 'E' in res else None
        if 'A' in parts:
            bf = h['A'].get()
            res['A'] = [dict(w=np.array(w), id=drift_at(w)[0], bf=sum(r['val'] for r in bf if r['w'] == w), res=max(r['res'] for r in bf if r['w'] == w)) for w in W_A]
            E.save(raw, **res); log('A')
        if 'F' in parts:
            res['F'] = part_F(field); E.save(raw, **res); log('F')
        if 'G' in parts:
            res['G'] = part_G(h['G'].get()); E.save(raw, **res); log('G')
        if 'H' in parts:
            res['H_rec'] = np.array([r['rec'] for r in h['H'].get()])
            res['H'] = part_H(res['H_rec'], field); E.save(raw, **res); log('H')
        elif 'S' in parts:                                                   # the limit SDE again, on the stored full-loop records
            res['H'] = part_H(np.asarray(res['H_rec'], float), field); E.save(raw, **res); log('H (SDE)')
        if 'I' in parts:
            res['I'] = part_I(field, res.get('F')); E.save(raw, **res); log('I')
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEFGHI']))


if __name__ == '__main__':
    main()
