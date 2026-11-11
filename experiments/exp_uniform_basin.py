"""
Theorem D‴.18: Theorem J uniformly in m at a replicated configuration — the per-particle basin and the m-independent rate
=====================================================================================================================

The real closed loop (CL) of Theorem J on the three-bump instance (k_e = 25, σ = σ₀/k_e, σ₀ = 0.15σ_crit; RK4 plant with the control and the
reference evaluated at the stages and the parameters frozen over the step, Euler parameters, dt = 5·10⁻³), integrated for many systems at
once (all particles in one array, the plant states in another; each system's f̂ by a segment sum).

A. The matching prior at the matching masses q = (1/3, 1/3, 1/3): the lift is the exact representation, a constant orbit with zero tracking
   error. Sizes m = 3, 15, 60, 240; every particle displaced from its class atom by ρ in a random direction of the G⁻¹-unit sphere,
   ρ = 0.3, 1.0, 2.5 (inside the Voronoi cells, whose distance from the atoms is ≥ Δ_ν/2 = 3.13) and 4.5 (across a cell boundary for some
   particles); the same displacements under the nearest-atom leak and under the OT leak. Records: d_∞ = max_i min_k‖ϑ_i − ϑ⁰_k‖_{G⁻¹}, the
   class counts by nearest atom, W²_{G⁻¹}(μ_Θ, ν), the rms tracking error.
B. A moving orbit: the OT leak toward the prior with centres shifted by +0.15 (the residual is not zero, so the weighted 3-atom orbit Γ_q moves
   and the tracking error is not zero), masses q; the reference is the 3-atom loop started at Θ*(q), integrated alongside; sizes m = 15, 60,
   240 started at the lift of Θ*(q) displaced per particle by ρ = 0.3, 1.0; d_∞ = max_i‖ϑ_i − ϑ^{ref}_{π₀(i)}‖_{G⁻¹} against the lift of the
   reference trajectory.

C. Remark D‴.19's mechanism: the lift of q at the matching prior with ONE particle of class 0 moved 60% of the way to atom 1 (into cell 1),
   m = 60, 240, both leaks — a vanishing mean-square displacement.
F. The exponents that test the rates: the 3-atom orbit of the real loop by shooting (Parts A and B), its parameter Floquet exponents against
   λ_j(H_q)/k_e, the ghosts' exponents per class against (σ₀ − λ_j(R^G_k))/k_e (Proposition D‴.15); and a sampled lower bound on the
   |||·|||-Lipschitz constant C_N of D⟨H⟩ at the replicated configuration (m = 15), hence an upper bound λ₀/(4C_N) on the proof's radius ρ.

Theorem D‴.18(b) predicts that every d_∞ started in the per-particle ball of radius ρ/2 decays at the rate (1 − κ)min{λ_coll, λ_sep}/k_e — the
same for every m — with an m-independent constant. The proof's ρ on this instance is small (Part F); the runs of Parts A and B use per-particle
radii far beyond it and probe how far the basin extends, which the theorem does not claim.

Output: results/chaos/uniform_basin.md and uniform_basin_raw.npz (UB_T for the horizon). `--from-raw` regenerates the report; `--only C,F`
runs the listed parts and merges them into the existing raw file (the other parts' records are kept).
"""
import os
import sys
import time
try:
    sys.stdout.reconfigure(encoding='utf-8')
except Exception:
    pass
import numpy as np
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.core import langevin as L
from src.core import prior_design as PD
from src.core import ot_leak as OT
from src.core import lifted_orbit as LO

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
KE = 25.0
SIG = SIGMA0 / KE
DT = 5e-3
T_END = float(os.environ.get('UB_T', 8000.0))
PLAN_EVERY, REC_EVERY, E_EVERY = 200, 1000, 20
MS_A, RHOS_A = (3, 15, 60, 240), (0.3, 1.0, 2.5, 4.5)
MS_B, RHOS_B = (15, 60, 240), (0.3, 1.0)
MS_C, T_MOVE = (60, 240), 0.6
Q = np.ones(3) / 3
A_MAX, C_MAX, S_MIN, S_MAX = 4.0, 3.0, 0.3, 2.0


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact_representation(net, c_shift=0.0):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + c_shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def sphere(rng, n):
    v = rng.standard_normal((n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


class Batch:
    """Many real-loop systems integrated together: particles Th (P, 3) with owner (P,), plant states X (S,), per-system leak ('nearest' or
    'ot') and prior atoms (S, K, 3); the step is the scheme of exp_ot_real_loop.py."""

    def __init__(self, systems):
        self.S = len(systems)
        self.Th = np.vstack([s['Th0'] for s in systems])
        self.owner = np.concatenate([np.full(s['Th0'].shape[0], k) for k, s in enumerate(systems)])
        self.m = np.array([s['Th0'].shape[0] for s in systems])
        self.inv_m = 1.0 / self.m[self.owner]
        self.leak = [s['leak'] for s in systems]
        self.prior = np.array([s['prior'] for s in systems])
        self.X = np.full(self.S, 0.0)                               # x_m(0) = 0: e(0) = 0
        self.slices = np.concatenate([[0], np.cumsum(self.m)])
        self.T = self.targets()

    def physical(self, Th):
        a = A_MAX * np.tanh(Th[:, 0] / A_MAX)
        c = C_MAX * np.tanh(Th[:, 1] / C_MAX)
        s = S_MIN + (S_MAX - S_MIN) / (1.0 + np.exp(-np.clip(Th[:, 2], -500, 500)))
        return a, c, s

    def f_hat(self, Y, a, c, s):
        y = Y[self.owner]
        phi = np.exp(-(y - c) ** 2 / (2 * s ** 2))
        return np.bincount(self.owner, weights=self.inv_m * a * phi, minlength=self.S)

    def targets(self):
        T = np.empty_like(self.Th)
        for k in range(self.S):
            sl = slice(self.slices[k], self.slices[k + 1])
            Th = self.Th[sl]; P = self.prior[k]
            C = np.sum((Th[:, None, :] - P[None, :, :]) ** 2 / G[None, None, :], axis=2)
            if self.leak[k] == 'nearest':
                T[sl] = P[np.argmin(C, axis=1)]
            else:
                m = Th.shape[0]; rep = np.repeat(np.arange(P.shape[0]), m // P.shape[0])
                rows, cols = linear_sum_assignment(C[:, rep])
                lab = np.empty(m, int); lab[rows] = rep[cols]
                T[sl] = P[lab]
        return T

    def step(self, t, refresh):
        if refresh:
            self.T = self.targets()
        a, c, s = self.physical(self.Th)
        X = self.X
        xr = lambda tau: 2.0 * np.sin(tau); xrd = lambda tau: 2.0 * np.cos(tau)
        xdot = lambda Y, tau: f_vec(Y) + xrd(tau) - KE * (Y - xr(tau)) - self.f_hat(Y, a, c, s)
        k1 = xdot(X, t); k2 = xdot(X + 0.5 * DT * k1, t + 0.5 * DT); k3 = xdot(X + 0.5 * DT * k2, t + 0.5 * DT); k4 = xdot(X + DT * k3, t + DT)
        e = X - xr(t)
        # grad_g per particle with its system's (x, e)
        x = X[self.owner]; ee = e[self.owner]
        diff = x - c
        phi = np.exp(-diff ** 2 / (2 * s ** 2))
        da = 1.0 - (a / A_MAX) ** 2; dc = 1.0 - (c / C_MAX) ** 2
        sk = 1.0 / (1.0 + np.exp(-np.clip(self.Th[:, 2], -500, 500))); ds = (S_MAX - S_MIN) * sk * (1 - sk)
        g = np.empty_like(self.Th)
        g[:, 0] = phi * ee * da
        g[:, 1] = a * ee * phi / s ** 2 * diff * dc
        g[:, 2] = a * ee * phi * diff ** 2 / s ** 3 * ds
        self.Th = self.Th + DT * (G[None, :] * g - SIG * (self.Th - self.T))
        self.X = X + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return e

    def system(self, k):
        return self.Th[self.slices[k]:self.slices[k + 1]]


def w2_to_prior(Th, P):
    """W²_{G⁻¹}(μ_Θ, ν) for ν uniform on the rows of P (m divisible by K): the assignment problem on the replicated atoms."""
    m = Th.shape[0]; rep = np.repeat(np.arange(P.shape[0]), m // P.shape[0])
    C = np.sum((Th[:, None, :] - P[None, rep, :]) ** 2 / G[None, None, :], axis=2)
    r, cidx = linear_sum_assignment(C)
    return float(C[r, cidx].mean())


def run_batch(part, systems, meta):
    B = Batch(systems)
    n = int(round(T_END / DT))
    rec = {k: dict(t=[], dinf=[], counts=[], w2=[]) for k in range(B.S)}
    e2acc = np.zeros(B.S); e2n = 0
    t0 = time.time()
    for i in range(n):
        t = i * DT
        if i % REC_EVERY == 0:
            for k in range(B.S):
                Th = B.system(k); P = B.prior[k]
                d = np.sqrt(np.sum((Th[:, None, :] - P[None, :, :]) ** 2 / G[None, None, :], axis=2))
                lab = np.argmin(d, axis=1)
                rec[k]['t'].append(t); rec[k]['counts'].append(np.bincount(lab, minlength=P.shape[0]))
                if part in ('A', 'C'):
                    rec[k]['dinf'].append(float(d.min(axis=1).max()))
                    rec[k]['w2'].append(w2_to_prior(Th, P))
                else:
                    ref = B.system(meta['ref_index'])
                    lab0 = systems[k]['labels0']
                    rec[k]['dinf'].append(float(np.sqrt(np.sum((Th - ref[lab0]) ** 2 / G[None, :], axis=1)).max()))
                    rec[k]['w2'].append(float('nan'))
        e = B.step(t, i % PLAN_EVERY == 0)
        if i % E_EVERY == 0 and t > T_END - 10 * PER:
            e2acc += e ** 2; e2n += 1
    out = []
    for k, s in enumerate(systems):
        r = {kk: np.array(v) for kk, v in rec[k].items()}
        r.update(dict(leak=s['leak'], m=int(B.m[k]), rho=s['rho'], seed=s['seed'], e_rms=float(np.sqrt(e2acc[k] / max(e2n, 1))), final=B.system(k).copy()))
        out.append(r)
    print(f"[part {part}] {B.S} systems, {B.Th.shape[0]} particles, {time.time() - t0:.0f}s", flush=True)
    return out


def predictions(net, prior, w=Q, full=False):
    """λ_sep = min_k λ_min(σ₀I − R^G_k) and λ_coll = λ_min(H_w) at the weighted stationary point Θ*(w) of the nearest-atom loop toward `prior`;
    with full=True also the eigenpairs of H_w (whitened, weighted coordinates √w_k G^{−1/2}u_k) and the separation spectra per class."""
    f_x = lambda X: f_vec(X)
    loop = PD.MixturePriorLoop(net, f_x, lambda t: np.array([2.0 * np.sin(t)]), PER, G, SIGMA0, prior, n_quad=256)
    atoms, res = PD.weighted_stationary_atoms_prior(loop, w, prior)
    a, c, s, diff, phi = loop.features(atoms); r = loop.F - phi @ (w * a)
    R = loop.residual_curvature(r, atoms); Rw = SQG[None, :, None] * R * SQG[None, None, :]
    lam_sep = min(float(np.linalg.eigvalsh(SIGMA0 * np.eye(3) - Rw[k]).min()) for k in range(3))
    K = 3; d = 3 * K; v = atoms.ravel(); Jm = np.zeros((d, d)); h = 1e-6
    gf = lambda V: PD.weighted_energy_grad_prior(loop, w, V.reshape(K, 3)).ravel()
    for j in range(d):
        ej = np.zeros(d); ej[j] = h
        Jm[:, j] = (gf(v + ej) - gf(v - ej)) / (2 * h)
    sw = np.repeat(np.sqrt(w), 3); sg = np.tile(SQG, K)
    Hs = (sw * sg)[:, None] * Jm * (sg / sw)[None, :]; Hs = 0.5 * (Hs + Hs.T)
    lam_coll = float(np.linalg.eigvalsh(Hs).min())
    if full:
        ev, V = np.linalg.eigh(Hs)
        sep = [np.sort(np.linalg.eigvalsh(SIGMA0 * np.eye(3) - Rw[k])) for k in range(3)]
        return atoms, res, lam_sep, lam_coll, ev, V, sep
    return atoms, res, lam_sep, lam_coll


def norm3(u, lab, w):
    """|||u||| = max{(Σ_k w_k‖ū_k‖²_{G⁻¹})^{1/2}, max_i‖u_i − ū_{lab(i)}‖_{G⁻¹}}."""
    K = len(w)
    ubar = np.array([u[lab == k].mean(axis=0) for k in range(K)])
    coll = np.sqrt(np.sum(np.asarray(w) * np.sum(ubar ** 2 / G[None, :], axis=1)))
    sep = np.sqrt(np.sum((u - ubar[lab]) ** 2 / G[None, :], axis=1)).max()
    return max(coll, sep)


def jacobian_lipschitz_lower_bound(net, prior, atoms, m=15, delta=0.05, n_samples=40, seed=7, h=1e-5):
    """A sampled lower bound on C_N := the |||·|||-Lipschitz constant of D⟨H⟩ at Θ*_m(q): max over displacements u (every particle moved by
    δ in a random direction, plus one-particle and collective displacements) of |||(D⟨H⟩(Θ* + u) − D⟨H⟩(Θ*))u||| / |||u|||², with ⟨H⟩ the
    averaged slow drift and directional derivatives by central differences."""
    loop = PD.MixturePriorLoop(net, f_vec, lambda t: np.array([2.0 * np.sin(t)]), PER, G, SIGMA0, prior, n_quad=256)
    lab = np.repeat(np.arange(3), m // 3)
    Th0 = atoms[lab]
    dirder = lambda Y, u: (loop.drift(Y + h * u) - loop.drift(Y - h * u)) / (2 * h)
    rng = np.random.default_rng(seed)
    us = []
    for _ in range(n_samples):
        u = rng.standard_normal((m, 3)) * SQG[None, :]
        us.append(delta * u / np.sqrt(np.sum(u ** 2 / G[None, :], axis=1))[:, None])
    u = np.zeros((m, 3)); u[0] = delta * SQG * np.array([1.0, 0.0, 0.0]); us.append(u)
    u = np.zeros((m, 3)); u[0] = delta * SQG * np.array([0.0, 1.0, 0.0]); u[1] = -u[0]; us.append(u)
    v = rng.standard_normal((3, 3)) * SQG[None, :]; us.append(delta * (v / np.sqrt(np.sum(v ** 2 / G[None, :], axis=1))[:, None])[lab])
    best = 0.0
    for u in us:
        d = dirder(Th0 + u, u) - dirder(Th0, u)
        best = max(best, norm3(d, lab, Q) / norm3(u, lab, Q) ** 2)
    return best


def floquet_job(net, prior, atoms):
    """The 3-atom orbit of the real loop with class targets `prior` (the OT leak at the masses q), shot from Θ*(q): its parameter Floquet
    exponents and each class's ghost exponents (Proposition D‴.15), with the residuals."""
    xr = lambda t: np.array([2.0 * np.sin(t)]); xrd = lambda t: np.array([2.0 * np.cos(t)])
    clK = LO.ClassLeakLoop(net, f_vec, xr, xrd, PER, KE, G, SIGMA0, prior)
    zK, DPK, res = clK.shoot(np.concatenate([xr(0.0), atoms.ravel()]), 3, max_iter=12, tol=1e-10)
    ex = np.sort(-np.log(np.abs(np.linalg.eigvals(DPK))) / PER)
    ex = ex[ex < 1.0]                                                      # the error channel's exponent (≈ k_e) is at the floor of the differences
    ghosts = {}
    for k in range(3):
        DPg, gex, gres = LO.ghost_monodromy(clK, zK, 3, k, prior[k])
        ghosts[k] = (np.sort(np.asarray(gex)), float(gres))
    return dict(res=float(res), offset=float(np.abs(zK[1:] - atoms.ravel()).max()), ex=ex, ghosts=ghosts)


def job(part):
    net = make_net()
    if part == 'A':
        prior = exact_representation(net)
        systems = []
        for m in MS_A:
            for rho in RHOS_A:
                seed = 1000 * m + int(10 * rho)
                rng = np.random.default_rng(seed)
                Th0 = np.repeat(prior, m // 3, axis=0) + rho * sphere(rng, m) * SQG[None, :]
                for leak in ('nearest', 'ot'):
                    systems.append(dict(Th0=Th0.copy(), leak=leak, prior=prior, rho=rho, seed=seed))
        return dict(part='A', runs=run_batch('A', systems, {}), prior=prior)
    if part == 'C':
        prior = exact_representation(net)
        systems = []
        for m in MS_C:
            Th0 = np.repeat(prior, m // 3, axis=0)
            Th0[0] = (1.0 - T_MOVE) * prior[0] + T_MOVE * prior[1]                # one particle of class 0 moved into cell 1
            for leak in ('nearest', 'ot'):
                systems.append(dict(Th0=Th0.copy(), leak=leak, prior=prior, rho=float('nan'), seed=0))
        return dict(part='C', runs=run_batch('C', systems, {}), prior=prior)
    if part == 'F':
        priorA = exact_representation(net)
        priorB = exact_representation(net, c_shift=0.15)
        atomsB = predictions(net, priorB)[0]
        out = dict(part='F')
        for key, pr, at in (('A', priorA, priorA), ('B', priorB, atomsB)):
            t1 = time.time()
            out[key] = floquet_job(net, pr, at)
            out[key]['CN_lower'] = jacobian_lipschitz_lower_bound(net, pr, at)
            print(f"[part F/{key}] {time.time() - t1:.0f}s", flush=True)
        return out
    prior = exact_representation(net, c_shift=0.15)
    atoms, res, lam_sep, lam_coll = predictions(net, prior)
    systems = [dict(Th0=atoms.copy(), leak='ot', prior=prior, rho=0.0, seed=0, labels0=np.arange(3))]
    for m in MS_B:
        for rho in RHOS_B:
            seed = 2000 * m + int(10 * rho)
            rng = np.random.default_rng(seed)
            lab0 = np.repeat(np.arange(3), m // 3)
            systems.append(dict(Th0=atoms[lab0] + rho * sphere(rng, m) * SQG[None, :], leak='ot', prior=prior, rho=rho, seed=seed, labels0=lab0))
    return dict(part='B', runs=run_batch('B', systems, dict(ref_index=0)), prior=prior, atoms=atoms, res=res, lam_sep=lam_sep, lam_coll=lam_coll)


def decay_rate(t, d, frac=0.5):
    msk = (t > frac * t[-1]) & (d > 0)
    if msk.sum() < 5:
        return float('nan')
    return float(-np.polyfit(t[msk], np.log(d[msk]), 1)[0])


def C_of_T(t, d, mu, Ts):
    """C(T) = max_{t ≤ T} d(t)e^{μt}/d(0) for the horizons Ts, and whether the maximum over the whole run sits at the last record."""
    g = d * np.exp(mu * t) / d[0]
    return [float(g[t <= T].max()) for T in Ts], bool(np.argmax(g) == len(g) - 1)


def count_history(r):
    """Number of record-to-record changes of the nearest-atom counts and the time of the last one (nan if none)."""
    c = np.asarray(r['counts'])
    ch = np.where(np.any(c[1:] != c[:-1], axis=1))[0]
    return int(ch.size), (float(r['t'][ch[-1] + 1]) if ch.size else float('nan'))


def nearest_labels(Th, P):
    return np.argmin(np.sum((Th[:, None, :] - P[None, :, :]) ** 2 / G[None, None, :], axis=2), axis=1)


def end_weights_analysis(loop, prior, final):
    """The realised weights at the end (by nearest atom), the weighted stationary point of the nearest-atom loop for them (Newton, by
    continuation from the uniform weights on the non-empty classes), its residual and splitting margin σ₀ − max κ_R, the largest distance of
    the final configuration from the lift of that point, the lift's W² to ν, and the within-class spread of the final configuration."""
    Th = np.asarray(final, float); m = Th.shape[0]
    lab = nearest_labels(Th, prior)
    counts = np.bincount(lab, minlength=3); nz = np.where(counts > 0)[0]
    w = counts[nz] / m
    Tw, res = PD.weighted_stationary_atoms_prior(loop, w, prior[nz], w0=np.ones(nz.size) / nz.size, n_cont=4)
    a, c, s, diff, phi = loop.features(Tw)
    kR = loop.kappa_R(loop.F - phi @ (w * a), Tw)
    idx = {int(k): j for j, k in enumerate(nz)}
    dist = max(float(np.sqrt(np.sum((Th[i] - Tw[idx[int(lab[i])]]) ** 2 / G))) for i in range(m))
    lifted = np.vstack([np.repeat(Tw[j:j + 1], counts[k], axis=0) for j, k in enumerate(nz)])
    spread = []
    for k in nz:
        X = Th[lab == k]
        dd = np.sqrt(np.sum((X - X.mean(axis=0)) ** 2 / G[None, :], axis=1))
        spread.append((float(np.sqrt(np.mean(dd ** 2))), float(dd.max())))
    return dict(counts=counts, res=float(res), margin=float(SIGMA0 - kR.max()), dist=dist, w2_lift=w2_to_prior(lifted, prior), spread=spread)


def initial_A(prior, m, rho, seed):
    rng = np.random.default_rng(seed)
    return np.repeat(prior, m // 3, axis=0) + rho * sphere(rng, m) * SQG[None, :]


def report(A, B, C=None, F=None, t0=None):
    net = make_net()
    priorA = exact_representation(net)
    _, _, lsepA, lcollA, evA, VA, sepA = predictions(net, priorA, full=True)
    sepB = predictions(net, B['prior'], full=True)[6]
    muA = min(lsepA, lcollA) / KE; muB = min(B['lam_sep'], B['lam_coll']) / KE
    mu3 = lcollA / KE                                                        # m = 3: one particle per class, no separation modes
    Ts = (2000.0, 4000.0, 6000.0, 8000.0)
    loopA = PD.MixturePriorLoop(net, f_vec, lambda t: np.array([2.0 * np.sin(t)]), PER, G, SIGMA0, priorA, n_quad=256)
    lines = ["# Theorem D‴.18: Theorem J uniformly in m at a replicated configuration — the per-particle basin and the m-independent rate (`exp_uniform_basin.py`)", "",
             f"The real closed loop of Theorem J on the three-bump instance: k_e = {KE:g}, σ = σ₀/k_e with σ₀ = 0.15σ_crit = {SIGMA0:.5f}, RK4 plant with the control and the reference at the stages and the parameters "
             f"frozen over the step, Euler parameters, dt = {DT:g}, horizon {T_END:g} time units ({T_END / PER:.0f} periods, {T_END / KE:.0f} slow units); the leak target refreshed every {PLAN_EVERY} steps; records every "
             f"{REC_EVERY * DT:g} time units. Fitted rates: least-squares slope of log d_∞ over the second half of the run; they are compared with the exponents of Part F, which are the test of the rate.", "",
             f"## A. The matching prior at the matching masses (the lift is the exact representation, a constant orbit)", "",
             f"Prediction: λ_sep = {lsepA:.5f} (= σ₀ at zero residual), λ_coll = {lcollA:.5f}, rate min{{λ_coll, λ_sep}}/k_e = {muA:.3e} per time unit for m ≥ 15; for m = 3 (one particle per class, no separation modes) "
             f"the slowest mode is the collective one, λ_coll/k_e = {mu3:.3e}, followed by λ_j(H_q)/k_e = " + ", ".join(f"{x / KE:.3e}" for x in evA[1:3]) + ". Every particle is displaced from its class atom by ρ "
             "in a random direction of the G⁻¹-unit sphere (the same displacements under both leaks); the Voronoi cells extend at least Δ_ν/2 = 3.13 from the atoms. "
             "The rate and C are printed only for runs that end at the counts q with constant counts over the second half.", "",
             "| leak | ρ | m | final counts by nearest atom | d_∞(0) → d_∞(end) | fitted rate of d_∞ (second half) | C(T) = max_{t≤T} d_∞(t)e^{μt}/d_∞(0) at T = 2000, 4000, 6000, 8000 (μ the prediction for that m) | W²(μ_Θ, ν): start → end | rms tracking error (last 10 periods) |",
             "|---|---|---|---|---|---|---|---|---|"]
    growing = []
    for leak in ('nearest', 'ot'):
        for rho in RHOS_A:
            for m in MS_A:
                r = next(x for x in A['runs'] if x['leak'] == leak and x['rho'] == rho and x['m'] == m)
                c = np.asarray(r['counts']); half = r['t'] > 0.5 * r['t'][-1]
                ok = bool(np.all(c[-1] == m // 3) and np.all(c[half] == c[-1]))
                mu = mu3 if m == 3 else muA
                if ok:
                    rate = f"{decay_rate(r['t'], r['dinf']):.3e}"
                    cT, last = C_of_T(r['t'], r['dinf'], mu, Ts)
                    cs = ", ".join(f"{x:.3f}" for x in cT) + (" (maximum at the last record)" if last else "")
                    if last and rho < 3.13:
                        growing.append(f"{leak}, ρ = {rho:g}, m = {m}")
                else:
                    rate, cs = "— (does not end at q)", "—"
                lines.append(f"| {leak} | {rho:g} | {m} | {tuple(int(x) for x in r['counts'][-1])} | {r['dinf'][0]:.3f} → {r['dinf'][-1]:.2e} | {rate} | {cs} | {r['w2'][0]:.3e} → {r['w2'][-1]:.2e} | {r['e_rms']:.1e} |")
    inside = [x for x in A['runs'] if x['rho'] < 3.13]
    big = [x for x in inside if x['m'] >= 15]; one = [x for x in inside if x['m'] == 3]
    rb = [decay_rate(x['t'], x['dinf']) for x in big]; r3 = [decay_rate(x['t'], x['dinf']) for x in one]
    Cb = [C_of_T(x['t'], x['dinf'], muA, Ts)[0][-1] for x in big]; C3 = [C_of_T(x['t'], x['dinf'], mu3, Ts)[0][-1] for x in one]
    ident = max(float(np.abs(np.asarray(x['dinf']) - np.asarray(next(y for y in inside if y['leak'] == 'ot' and y['rho'] == x['rho'] and y['m'] == x['m'])['dinf'])).max())
                for x in inside if x['leak'] == 'nearest')
    lines += ["", f"Inside the cells (ρ = 0.3, 1.0, 2.5; both leaks): the two leaks give identical runs (largest difference of d_∞ {ident:.1e}); for m ≥ 15 fitted rates {min(rb):.3e}–{max(rb):.3e} against σ₀/k_e = {muA:.3e} "
              f"(ratios {min(rb) / muA:.3f}–{max(rb) / muA:.3f}), C(8000) = {min(Cb):.3f}–{max(Cb):.3f}; for m = 3 fitted rates {min(r3):.3e}–{max(r3):.3e} against λ_coll/k_e = {mu3:.3e} (ratios {min(r3) / mu3:.3f}–{max(r3) / mu3:.3f}), "
              f"C(8000) = {min(C3):.3f}–{max(C3):.3f}. C is still growing at the horizon (maximum of d_∞e^{{μt}} at the last record) for: " + ("; ".join(growing) if growing else "none") + "."]
    # m = 3: the final deviation on the collective eigenmodes (whitened, weighted coordinates)
    amp = []
    for rho in RHOS_A[:3]:
        r = next(x for x in A['runs'] if x['leak'] == 'ot' and x['rho'] == rho and x['m'] == 3)
        y = (np.sqrt(1.0 / 3.0) * (np.asarray(r['final'], float) - priorA) / SQG[None, :]).ravel()
        a3 = np.abs(VA[:, :3].T @ y)
        amp.append(f"ρ = {rho:g}: " + ", ".join(f"{x:.4f}" for x in a3) + f" (fitted rate {decay_rate(r['t'], r['dinf']):.3e})")
    lines += ["", "m = 3: the final deviation from the exact representation on the three slowest collective modes (amplitudes in the whitened, weighted coordinates "
              "√w_k G^{−1/2}u_k; modes at λ_j(H_q)/k_e = " + ", ".join(f"{x / KE:.3e}" for x in evA[:3]) + "): " + "; ".join(amp) + ". A fitted rate is that of the mode that carries the end of the run: "
              "the runs are mixtures of modes, and only a run whose end is carried by the slowest mode shows λ_coll/k_e."]
    # A′: the cell bookkeeping at ρ = 4.5
    lines += ["", "### A′. Across a cell boundary (ρ = 4.5): the cell bookkeeping", "",
              "Classes by nearest prior atom. For the nearest-atom leak: the weighted stationary point for the realised end weights (Newton by continuation), its residual and splitting margin σ₀ − max κ_R "
              "(negative: the lift is a saddle), the largest distance of the final configuration from the lift of that point and the lift's W² to ν. For the OT leak: the particles that end outside the cell of their "
              "initial class (the plan relabels) and the log-slope of W² over t > 6000 per slow unit, against the leak's 2σ₀ = " + f"{2 * SIGMA0:.5f}.", "",
              "| leak | m | particles outside their class cell at t = 0 | counts at t = 0 → end | count changes (records) / last change at t | within-class spread at the end, per class (rms / max) | end weights: Newton residual, splitting margin | max distance from the lift of the end weights | W²(μ_Θ, ν) end / W² of that lift | particles ending outside their initial class's cell | log-slope of W² over t > 6000 (per slow unit) |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    for leak in ('nearest', 'ot'):
        for m in MS_A:
            r = next(x for x in A['runs'] if x['leak'] == leak and x['rho'] == 4.5 and x['m'] == m)
            Th0 = initial_A(priorA, m, 4.5, r['seed']); lab0 = np.repeat(np.arange(3), m // 3)
            outside0 = int(np.sum(nearest_labels(Th0, priorA) != lab0))
            nch, tlast = count_history(r)
            fin = np.asarray(r['final'], float)
            moved = int(np.sum(nearest_labels(fin, priorA) != lab0))
            tl = f"{tlast:.0f}" if np.isfinite(tlast) else "—"
            if leak == 'nearest':
                ea = end_weights_analysis(loopA, priorA, fin)
                sp = "; ".join(f"{a:.2f} / {b:.2f}" for a, b in ea['spread'])
                conv = ea['res'] < 1e-8
                ew = f"{ea['res']:.1e}, {ea['margin']:+.5f}" if conv else f"no stationary point found (residual {ea['res']:.1e})"
                dl = f"{ea['dist']:.3f}" if conv else "—"
                wl = f"{ea['w2_lift']:.3f}" if conv else "—"
                lines.append(f"| {leak} | {m} | {outside0} | {tuple(int(x) for x in r['counts'][0])} → {tuple(int(x) for x in r['counts'][-1])} | {nch} / {tl} | {sp} | {ew} | {dl} | {r['w2'][-1]:.3f} / {wl} | {moved} | — |")
            else:
                msk = r['t'] > 6000.0
                sl = float(np.polyfit(r['t'][msk], np.log(np.asarray(r['w2'])[msk]), 1)[0]) * KE
                lab = nearest_labels(fin, priorA)
                sp = []
                for k in range(3):
                    X = fin[lab == k]; dd = np.sqrt(np.sum((X - X.mean(axis=0)) ** 2 / G[None, :], axis=1))
                    sp.append(f"{np.sqrt(np.mean(dd ** 2)):.3f} / {dd.max():.3f}")
                lines.append(f"| {leak} | {m} | {outside0} | {tuple(int(x) for x in r['counts'][0])} → {tuple(int(x) for x in r['counts'][-1])} | {nch} / {tl} | {'; '.join(sp)} | — | — | {r['w2'][-1]:.3f} / — | {moved} | {sl:.5f} |")
    # B
    lines += ["", "## B. A moving orbit: the OT leak toward the prior with centres shifted by +0.15 (non-zero residual), masses q", "",
              f"Θ*(q) by Newton (residual {B['res']:.1e}); prediction λ_sep = {B['lam_sep']:.5f}, λ_coll = {B['lam_coll']:.5f}, rate {muB:.3e} per time unit (the separation exponents (σ₀ − λ_j(R^G_k))/k_e: "
              + "; ".join(f"class {k}: " + ", ".join(f"{x / KE:.3e}" for x in sepB[k]) for k in range(3)) + "). The reference is the 3-atom loop from Θ*(q); d_∞ is measured against the lift of the reference trajectory.", "",
              "| ρ | m | d_∞(0) → d_∞(end) | fitted rate | C(T) at T = 2000, 4000, 6000, 8000 | rms tracking error (m-system / reference) |", "|---|---|---|---|---|---|"]
    ref = B['runs'][0]
    for r in B['runs'][1:]:
        rate = decay_rate(r['t'], r['dinf'])
        cT, last = C_of_T(r['t'], r['dinf'], muB, Ts)
        lines.append(f"| {r['rho']:g} | {r['m']} | {r['dinf'][0]:.3f} → {r['dinf'][-1]:.2e} | {rate:.3e} | " + ", ".join(f"{x:.3f}" for x in cT) + (" (maximum at the last record)" if last else "") + f" | {r['e_rms']:.2e} / {ref['e_rms']:.2e} |")
    ratesB = [decay_rate(r['t'], r['dinf']) for r in B['runs'][1:]]
    lines += ["", f"Fitted rates {min(ratesB):.3e}–{max(ratesB):.3e} against the prediction {muB:.3e} (ratios {min(ratesB) / muB:.3f}–{max(ratesB) / muB:.3f}): the second-half fits are mixtures of the separation "
              "modes, whose gaps (inverse gaps above 10⁴ time units) exceed the fit window; they are consistent with the bound (1 − κ)min{λ_coll, λ_sep}/k_e, and the exponents of Part F are the test of the rate."]
    # C
    if C is not None:
        loopC = loopA
        lines += ["", f"## C. One particle across a boundary (Remark D‴.19): the lift of q with one particle of class 0 moved {T_MOVE:g} of the way to atom 1", "",
                  "| leak | m | mean-square displacement (1/m)Σ‖u_i‖² at t = 0 | counts at the end | W²(μ_Θ, ν): start → end | d_∞ at the end | end weights: Newton residual, splitting margin | max distance from the lift of the end weights | W² of that lift |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for r in C['runs']:
            m = r['m']
            Th0 = np.repeat(C['prior'], m // 3, axis=0); Th0[0] = (1.0 - T_MOVE) * C['prior'][0] + T_MOVE * C['prior'][1]
            msd = float(np.mean(np.sum((Th0 - np.repeat(C['prior'], m // 3, axis=0)) ** 2 / G[None, :], axis=1)))
            fin = np.asarray(r['final'], float)
            if r['leak'] == 'nearest':
                ea = end_weights_analysis(loopC, C['prior'], fin)
                conv = ea['res'] < 1e-8
                ew = f"{ea['res']:.1e}, {ea['margin']:+.5f}" if conv else f"no stationary point found (residual {ea['res']:.1e})"
                tail = f"{ew} | {ea['dist']:.3f} | {ea['w2_lift']:.4f} |" if conv else f"{ew} | — | — |"
            else:
                tail = "— | — | — |"
            lines.append(f"| {r['leak']} | {m} | {msd:.4f} | {tuple(int(x) for x in r['counts'][-1])} | {r['w2'][0]:.4f} → {r['w2'][-1]:.2e} | {r['dinf'][-1]:.2e} | " + tail)
    # F
    if F is not None:
        lines += ["", "## F. The exponents that test the rates, and the proof's radius on this instance", ""]
        for key, ev, sep, lam0, name in (('A', evA, sepA, min(lsepA, lcollA), 'A (matching prior)'), ('B', predictions(net, B['prior'], full=True)[4], sepB, min(B['lam_sep'], B['lam_coll']), 'B (shifted prior)')):
            fr = F[key]
            ex = np.asarray(fr['ex'])
            pred = np.sort(ev) / KE
            nn = min(len(ex), len(pred))
            rel = np.abs(ex[:nn] - pred[:nn]) / pred[:nn]
            lines.append(f"**{name}.** The 3-atom orbit of the real loop by shooting: residual {fr['res']:.1e}, period-start offset from Θ*(q) {fr['offset']:.2e}; its parameter Floquet exponents "
                         + ", ".join(f"{x:.4e}" for x in ex[:nn]) + " against λ_j(H_q)/k_e = " + ", ".join(f"{x:.4e}" for x in pred[:nn]) + f" (relative errors {rel.min() * 100:.2f}–{rel.max() * 100:.2f}%). "
                         "The ghosts' exponents (separation exponents of every lift, Proposition D‴.15): " + "; ".join(f"class {k}: " + ", ".join(f"{x:.4e}" for x in fr['ghosts'][k][0]) + f" (prediction " + ", ".join(f"{x / KE:.4e}" for x in sep[k]) + f"; return residual {fr['ghosts'][k][1]:.1e})" for k in range(3))
                         + f". The slowest exponent overall: {min(ex.min(), min(fr['ghosts'][k][0].min() for k in range(3))):.4e}.")
            rho_p = lam0 / (4.0 * fr['CN_lower'])
            radii = RHOS_A[:3] if key == 'A' else RHOS_B
            lines.append("")
            lines.append(f"The proof's radius: a sampled lower bound C_N ≥ {fr['CN_lower']:.4f} on the |||·|||-Lipschitz constant of D⟨H⟩ at Θ*_15(q) (displacements of per-particle size 0.05) gives "
                         f"ρ ≤ λ₀/(4C_N) = {rho_p:.4f} with λ₀ = {lam0:.5f}; the runs start at per-particle radii {', '.join(f'{x:g}' for x in radii)}, i.e. {min(radii) / (rho_p / 2):.0f}–{max(radii) / (rho_p / 2):.0f} times the proved start radius ρ/2. "
                         "The contraction observed there is evidence that the basin extends to the cell margins, which the theorem does not claim.")
            lines.append("")
    lines += ["Numbers are produced by `experiments/exp_uniform_basin.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'uniform_basin.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'uniform_basin_raw.npz')
    if '--from-raw' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        report(z['A'].item(), z['B'].item(), z['C'].item() if 'C' in z.files else None, z['F'].item() if 'F' in z.files else None, t0=t0)
        return
    parts = ['A', 'B', 'C', 'F']
    old = {}
    if '--only' in sys.argv:
        parts = sys.argv[sys.argv.index('--only') + 1].split(',')
        z = np.load(raw, allow_pickle=True)
        old = {k: z[k].item() for k in z.files if k in ('A', 'B', 'C', 'F')}
    with mp.Pool(min(len(parts), 3)) as pool:
        res = pool.map(job, parts)
    new = dict(old)
    for p, part in zip(parts, res):
        for r in part.get('runs', []):
            r['final'] = r['final'].astype(np.float32)
        new[p] = part
    np.savez(raw, **{k: np.array(v, dtype=object) for k, v in new.items()})
    report(new['A'], new['B'], new.get('C'), new.get('F'), t0=t0)


if __name__ == '__main__':
    main()
