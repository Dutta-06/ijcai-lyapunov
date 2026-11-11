"""
Proposition D‴.50: the rich transit under label noise — the dual is a Langevin diffusion
====================================================================================================================================

The free mass–amplitude law on the 6 × 14 Gaussian designs of Proposition D‴.48 (`exp_mirror_transit`: 3-sparse truths, masses εb_i, b_i
log-uniform in [1/4, 4]) with an isotropic white label noise of intensity σ: dζ = (y − Φw)dt/N + (σ/N)dW (`dual_langevin`).

A. The Gibbs law on a small well-conditioned dictionary (4 × 8, orthonormal rows scaled by √8, masses 2b_i with b_i ∈ U(1, 2) — the
   lazy side, where every particle moves and the dual mixes in a few time units — and σ = 1, far from Gaussian): the mirror update (the
   dual's Euler–Maruyama, at two steps) against an independent random-walk Metropolis sampler of π ∝ exp(−F/T) — the coefficients' means
   and standard deviations, E[Φw] − y; the paradigm's law on its curves read in Stratonovich's sense by the stochastic Heun scheme in α
   follows the mirror update pathwise (the same noise; the strong error ∝ h) but drifts along the dual's conserved directions
   P_{ker Φ}η over long times (a discretisation that does not keep them samples another law); read in Itô's sense (Euler on
   (α, log β), the paradigm's update) the conserved quantities drift (the Heun scheme's drift along them ∝ h√t: a discretisation
   artefact, at fixed h no stationary law). The exact integration-by-parts identity Cov_π(Φw) = T·ΦE_π[D]Φᵀ on the mirror update's run
   (against the Laplace value at the mode); a small-T run (σ = 0.1) against (i)(b)'s T·D^{1/2}PD^{1/2} (P of rank N < d).
B. The fluctuation–dissipation relation on the support (ε = 10⁻⁶, three designs, σ = 0.03, 0.1, 0.3, four independent noise paths each):
   from the noiseless transit's end, over at least NRELAX_B relaxation times of the slowest mode, Var(w_i) against T·E[dw_i/dη] (E over
   the run) and T·dw_i/dη at w_BP, the hypothesis' size (T·D_i)^{1/2}/|w_i|, the correlations against their sampling error (0 predicted),
   the mean against weighted basis pursuit, the inactive coefficients.
C. The noisy transit (the first design, ε = 10⁻⁴, 10⁻⁶, 10⁻⁸, σ = 0.1, 0.3, 1, four noise paths, each one Brownian path for every ε and
   σ): the noise's share sup_τ|u − u_det| (u_det the noiseless run on the same grid) against κ = σ/(NΛ^{1/2}) and against its first-order
   response along the inverse scale space path (δν ← P_T(δν + κdW): for entries only κP_{T(τ)}W(τ)) — its error over all τ, outside fixed
   windows around the path's and the run's events and outside windows sized by the two mechanisms (the noise's entry shift 3κ/|u̇_j| and
   the entrant's finite-ε absorption 5N/(Λ·D_j|P_{T⁻}φ_j|²)), inside them (σ-independent: finite-ε, linear in the noise), and in L²; the
   residual outside the windows against the active block's Ornstein–Uhlenbeck component κ²Σ1/(2Λλ_i); the entries' order and times.
D. The local quasi-equilibrium in the rich regime (the first design, ε = 10⁻⁴, 10⁻⁶, 10⁻⁸, σ = 0.3, 1, four independent paths): over a
   window that is a small fraction of the face's crossing time Λ²N²/σ², the support's mean against weighted basis pursuit with its
   standard error, the inactive coefficients — not a sample of π, which (iii) concerns.

Output: results/chaos/transit_noise.md and transit_noise_raw.npz (`--from-raw` regenerates the md; `--only A,B,C,D`).
"""
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
from src.core import mirror_free as MF
from src.core import rich_transit as RT
from src.core import dual_langevin as DL
from src.core.runcache import cached
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('TN_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('TN_WORKERS', '2'))
NA_A, DA_A, EPS_A, SIG_A = 4, 8, 2.0, 1.0                                   # A: the small dictionary, masses 2b_i, the noise
H_A, T_A, N_RWM, BURN_A = (0.01, 0.005), 2.0e4, 1_000_000, 0.1               # A: mirror-update steps, horizon; Metropolis steps; burn-in
H_PATH_A, T_PATH_A, T_DRIFT_A = (0.01, 0.0025), 50.0, 2000.0                  # A: Heun against the mirror update, pathwise and long
SIG_A_SMALL = 0.1                                                        # A: (i)(b)'s Laplace expansion at small T
SEEDS_B, EPS_B, SIGS_B, H_B, TW_B = (3, 4, 5), 1e-6, (0.03, 0.1, 0.3), 0.1, 4.0e4
PATHS_B, NRELAX_B, H_MAX_B = 4, 30.0, 0.5                                # B: independent paths; the window ≥ NRELAX_B slowest relaxations
SEED_CD = 3                                                              # C, D: the first design of exp_mirror_transit
EPS_C, SIGS_C, PATHS_C, H_C, REC_C = (1e-4, 1e-6, 1e-8), (0.1, 0.3, 1.0), 4, 0.25, 20
WIN_C = (0.05, 0.1)                                                       # C: fixed windows around the path's and the run's events (in τ)
WIN_FAC_C = (1.0, 2.0, 4.0)                                               # C: the mechanisms' windows max(3κ/|u̇_j|, 5N/(Λ·D_j|P_{T⁻}φ_j|²)) ×
EPS_D, SIGS_D, H_D, TW_D, PATHS_D = (1e-4, 1e-6, 1e-8), (0.3, 1.0), 0.1, 4.0e4, 4
NB = 20                                                                  # batch means


def setup(seed, eps):
    L = X.law()
    Phi, ws, b = X.design(seed)
    y = Phi @ ws
    rho = b ** -0.5
    path = MF.iss_path(Phi, y, rho)
    q = L.Q(np.zeros(b.size), eps * b)
    return L, Phi, y, b, rho, path, q


def batch_se(x, nb=NB):
    """Standard error of the mean of a time series by batch means."""
    x = np.asarray(x)
    n = (x.shape[0] // nb) * nb
    bm = x[:n].reshape(nb, n // nb, *x.shape[1:]).mean(1)
    return bm.std(0, ddof=1) / np.sqrt(nb)


def plateau(L, Phi, y, q, path, eps):
    """The noiseless transit's end (1.5τ*), as the start of a noisy window."""
    F = DL.FastFree(L, q)
    Lam = RT.transit_scale(L, eps)
    ts, Z, W, al = DL.mirror_update(F, Phi, y, 0.0, 0.25, 1.5 * path['taus'][-1] * Lam, rec_every=10 ** 9)
    return F, Z[-1], al, W[-1]


# --- A ---------------------------------------------------------------------------------------------------------------------------------
def small_dictionary():
    """4 × 8, orthonormal rows scaled by √8 (ΦΦᵀ = 8I), y = Φ(0.6, −0.4, 0, …), masses 2b_i, b_i ∈ U(1, 2)."""
    L = X.law()
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.standard_normal((DA_A, DA_A)))
    Phi = np.sqrt(DA_A) * Q[:NA_A]
    y = Phi @ np.r_[0.6, -0.4, np.zeros(DA_A - 2)]
    b = rng.uniform(1.0, 2.0, DA_A)
    return L, Phi, y, b, L.Q(np.zeros(DA_A), EPS_A * b)


@cached
def a_job(args):
    kind, h, seed = args
    warnings.filterwarnings('ignore')
    L, Phi, y, b, q = small_dictionary()
    N, d = Phi.shape
    T = SIG_A ** 2 / (2 * N)
    F = DL.FastFree(L, q)
    t0 = time.time()
    if kind == 'rwm':
        ts, Z, W, al = DL.mirror_update(F, Phi, y, 0.0, 0.05, 2000.0, rec_every=10 ** 9)   # the mode: the Bregman projection's dual
        D = F.dw(al)
        cov = (2.38 ** 2 / N) * T * np.linalg.inv(Phi @ (D[:, None] * Phi.T))
        Zs, Ws, acc = DL.gibbs_rwm(F, Phi, y, T, Z[-1], cov, N_RWM, seed=seed, thin=10)
        Ws = Ws[int(BURN_A * len(Ws)):]
        return dict(kind=kind, h=h, mean=Ws.mean(0), sd=Ws.std(0), se=batch_se(Ws), fit=(Ws @ Phi.T).mean(0) - y, acc=acc,
                    secs=time.time() - t0, w_proj=W[-1], rates=np.linalg.eigvalsh(Phi @ (D[:, None] * Phi.T) / N))
    if kind == 'mirror':
        ts, Z, Ws, AL = DL.mirror_update(F, Phi, y, SIG_A, h, T_A, seed=seed, rec_every=10, rec_al=True)
        k0 = int(BURN_A * len(Ws))
        Ws, AL = Ws[k0:], AL[k0:]
        # the exact identity E_π[∇F∇Fᵀ] = T·E_π[∇²F], i.e. Cov_π(Φw) = T·ΦE_π[D]Φᵀ, against the Laplace value T·ΦD(w*)Φᵀ; the identity's
        # sampling error from its spread over the batches
        _, _, _, al_star = DL.mirror_update(F, Phi, y, 0.0, 0.05, 2000.0, rec_every=10 ** 9)
        Dt = F.dw(AL)
        rel = lambda Ws_, Dt_: float(np.linalg.norm(np.cov((Ws_ @ Phi.T).T) - T * Phi @ (Dt_.mean(0)[:, None] * Phi.T))
                                     / np.linalg.norm(T * Phi @ (Dt_.mean(0)[:, None] * Phi.T)))
        n = (len(Ws) // NB) * NB
        rb = [rel(a, b) for a, b in zip(Ws[:n].reshape(NB, -1, d), Dt[:n].reshape(NB, -1, d))]
        Pm = T * Phi @ (F.dw(al_star)[:, None] * Phi.T)
        return dict(kind=kind, h=h, mean=Ws.mean(0), sd=Ws.std(0), se=batch_se(Ws), fit=(Ws @ Phi.T).mean(0) - y, secs=time.time() - t0,
                    ident=rel(Ws, Dt), ident_batches=float(np.mean(rb)), ident_se=float(np.std(rb, ddof=1) / np.sqrt(NB)),
                    ident_mode=float(np.linalg.norm(np.cov((Ws @ Phi.T).T) - Pm) / np.linalg.norm(Pm)))
    if kind == 'laplace':                                                  # (i)(b) at small T: Laplace's expansion at the mode
        Ts = SIG_A_SMALL ** 2 / (2 * N)
        _, Zm, Wm, al_star = DL.mirror_update(F, Phi, y, 0.0, 0.05, 2000.0, rec_every=10 ** 9)
        ts, Z, Ws, al = DL.mirror_update(F, Phi, y, SIG_A_SMALL, h, T_A, seed=seed, rec_every=10, zeta0=Zm[-1], al0=al_star)
        Ws = Ws[int(BURN_A * len(Ws)):]
        Cp = DL.fd_covariance(Phi, F.dw(al_star), Ts)
        Cw = np.cov(Ws.T)
        return dict(kind=kind, h=h, T=Ts, cov_rel=float(np.linalg.norm(Cw - Cp) / np.linalg.norm(Cp)),
                    sd_ratio=np.sqrt(np.diag(Cw) / np.diag(Cp)), shift=(Ws.mean(0) - Wm[-1]) / Ws.std(0), se_rel=batch_se(Ws) / Ws.std(0),
                    rank=int(np.linalg.matrix_rank(Cp, tol=1e-10 * np.abs(Cp).max())), secs=time.time() - t0)
    if kind == 'path':                                                     # Heun on α against the mirror update, the same noise
        Pk = np.eye(d) - np.linalg.pinv(Phi) @ Phi                          # the projector onto ker Φ: P_{ker Φ}η is conserved
        out = dict(kind=kind, h=h)
        for name, T_end in (('short', T_PATH_A), ('long', T_DRIFT_A)):
            rec = max(1, int(round(1.0 / h)))
            ts, Z, Wm, al = DL.mirror_update(F, Phi, y, SIG_A, h, T_end, seed=seed, rec_every=rec)
            ts2, A, B = DL.primal_sde(L, Phi, y, np.zeros(d), np.log(EPS_A * b), SIG_A, h, T_end, seed=seed, rec_every=rec,
                                      scheme='stratonovich')
            eta_h = F.eta(np.abs(A)) * np.sign(A)
            out[name] = dict(dw=float(np.abs(np.exp(B) * L.amp(A)[0] / L.m - Wm).max()), deta=float(np.abs(eta_h - Z @ Phi).max()),
                             inv=np.abs(eta_h @ Pk).max(1), ts=ts)
        out['secs'] = time.time() - t0
        return out
    # 'ito': the paradigm's Euler update on (α, log β)
    ts, A, B = DL.primal_sde(L, Phi, y, np.zeros(d), np.log(EPS_A * b), SIG_A, h, T_A / 10, seed=seed, rec_every=10, scheme='ito')
    Ws = np.exp(B) * L.amp(A)[0] / L.m
    Qs = L.Q(A, np.exp(B))
    k1 = max(1, int(round(1.0 / (10 * h))))                                # t = 1
    return dict(kind=kind, h=h, mean=Ws[len(Ws) // 2:].mean(0), sd=Ws[len(Ws) // 2:].std(0), dq_end=Qs[-1] - q,
                rate1=(Qs[k1] - q) / ts[k1], drift0=DL.ito_q_drift(L, np.zeros(d), EPS_A * b, Phi, SIG_A), t_end=float(ts[-1]),
                secs=time.time() - t0)


def part_A(pool):
    jobs = [('rwm', 0.0, 0)] + [('mirror', h, 2) for h in H_A] + [('path', h, 5) for h in H_PATH_A] + [('ito', H_A[0], 3),
                                                                                                         ('laplace', H_A[0], 4)]
    out = {}
    for r in pool.map(a_job, jobs):
        out[f"{r['kind']} {r['h']}"] = r
        print(f"  A {r['kind']} h {r['h']}: {r['secs']:.0f}s", flush=True)
    return out


# --- B ---------------------------------------------------------------------------------------------------------------------------------
@cached
def b_job(args):
    seed, sig, k = args
    warnings.filterwarnings('ignore')
    L, Phi, y, b, rho, path, q = setup(seed, EPS_B)
    N, d = Phi.shape
    S = np.flatnonzero(path['w_end'])
    F, z0, al0, w0 = plateau(L, Phi, y, q, path, EPS_B)
    T = sig ** 2 / (2 * N)
    D0 = F.dw(al0)                                                         # at the noiseless end, w_BP
    G = Phi[:, S].T @ Phi[:, S]
    rates = np.sort(np.linalg.eigvals((G / N) @ np.diag(D0[S])).real)
    tw = max(TW_B, NRELAX_B / rates[0])                                    # at least NRELAX_B relaxation times of the slowest mode
    h = min(H_B * tw / TW_B, H_MAX_B)
    rec = max(1, int(round(tw / h / 4e5)))
    t0 = time.time()
    ts, Z, Ws, AL = DL.mirror_update(F, Phi, y, sig, h, tw, seed=100000 * seed + 100 * int(round(1000 * sig)) + k, rec_every=rec,
                                     zeta0=z0, al0=al0, rec_al=True)
    k0 = len(Ws) // 20
    Ws, AL = Ws[k0:], AL[k0:]
    ED = F.dw(AL).mean(0)                                                  # E_π[D] over the run
    C = np.cov(Ws[:, S].T)
    sd = np.sqrt(np.diag(C))
    iu = np.triu_indices(S.size, 1)
    n = (len(Ws) // NB) * NB
    cb = []
    for Wb in Ws[:n].reshape(NB, -1, d):
        Cb = np.cov(Wb[:, S].T)
        sb = np.sqrt(np.diag(Cb))
        cb.append((Cb / np.outer(sb, sb))[iu])
    off = np.setdiff1d(np.arange(d), S)
    return dict(seed=seed, sig=sig, k=k, h=h, tw=tw, S=S, ratio=np.diag(C) / (T * ED[S]), ratio_bp=np.diag(C) / (T * D0[S]),
                ED_rel=ED[S] / D0[S], sd_rel=np.sqrt(T * D0[S]) / np.abs(path['w_end'][S]), corr=(C / np.outer(sd, sd))[iu],
                corr_se=np.std(cb, axis=0, ddof=1) / np.sqrt(NB), bias=Ws[:, S].mean(0) - path['w_end'][S], se=batch_se(Ws[:, S]),
                w_bp=path['w_end'][S], inactive=float(np.abs(Ws[:, off]).max()) if off.size else 0.0,
                start=float(np.abs(w0 - path['w_end']).max()), rates=rates, n_relax=float(tw * rates[0]), secs=time.time() - t0)


def part_B(pool):
    jobs = [(s, g, k) for s in SEEDS_B for g in SIGS_B for k in range(PATHS_B)]
    rows = pool.map(b_job, jobs)
    for r in rows:
        print(f"  B design {r['seed']} σ {r['sig']} path {r['k']} (h {r['h']}, t {r['tw']:.3g}): Var/(T·E[dw/dη]) {np.round(r['ratio'], 3)} "
              f"({r['secs']:.0f}s)", flush=True)
    return rows


# --- C ---------------------------------------------------------------------------------------------------------------------------------
def entries_interp(W, taus, th):
    """The first τ at which |w_i| exceeds th, linearly interpolated between the records (inf if never)."""
    out = np.full(W.shape[1], np.inf)
    for i in range(W.shape[1]):
        k = np.flatnonzero(np.abs(W[:, i]) > th)
        if k.size:
            j = k[0]
            if j == 0:
                out[i] = taus[0]
            else:
                a, b = abs(W[j - 1, i]), abs(W[j, i])
                out[i] = taus[j - 1] + (th - a) / (b - a) * (taus[j] - taus[j - 1])
    return out


@cached
def c_job(args):
    eps, sig, k = args                                                     # sig = 0: the noiseless run on the same grid
    warnings.filterwarnings('ignore')
    L, Phi, y, b, rho, path, q = setup(SEED_CD, eps)
    N, d = Phi.shape
    F = DL.FastFree(L, q)
    Lam = RT.transit_scale(L, eps)
    r = int(round(RT.transit_scale(L, min(EPS_C)) / Lam))                  # the finest ε's steps per step: one Brownian path for every ε
    t_end = 1.5 * path['taus'][-1] * Lam
    nst = int(round(t_end / H_C))
    rng = np.random.default_rng(1000 + k)
    xi = np.empty((nst + 1, N))
    for j in range(nst + 1):
        xi[j] = rng.standard_normal((r, N)).sum(0) / np.sqrt(r)
    t0 = time.time()
    ts, Z, W, AL = DL.mirror_update(F, Phi, y, sig, H_C, t_end, rec_every=REC_C, noise=lambda j: xi[j], rec_al=True)
    taus = ts / Lam
    U = Z @ Phi / Lam
    th = 1e-3 * np.abs(path['w_end']).max()
    out = dict(eps=eps, sig=sig, k=k, taus=taus, U=U, nu=Z / Lam, entries=entries_interp(W, taus, th), secs=time.time() - t0, Lam=Lam,
               w_end=W[-1], iss=float(np.abs(U - np.array([MF.iss_eval(path, t, Phi, y)[0] for t in taus])).max()))
    if sig == 0:
        out['D'], out['W'] = F.dw(AL), W
    else:
        # the first-order response on the record grid: increments of W in τ, summed over each record block
        blocks = np.add.reduceat(xi[:nst], np.arange(0, nst, REC_C), axis=0)[:taus.size - 1] * np.sqrt(H_C / Lam)
        out['nu_lin'] = DL.linear_response(path, Phi, taus, blocks, sig / (N * np.sqrt(Lam)), return_nu=True)
    return out


def windows_C(path, Phi, y, det, kap, Lam):
    """Per event of the path (all entries): the noise's entry-shift width 3κ/|u̇_j| (u̇_j the entrant's rate before it) and the entrant's
    finite-ε absorption width 5N/(Λ·D_j|P_{T⁻}φ_j|²) (D_j = dw/dη of the noiseless run where its coefficient first reaches half its
    final value; P_{T⁻} the projector onto the face before the entry)."""
    N = Phi.shape[0]
    P = DL.tangent_projectors(path, Phi)
    out = []
    for q, (te, j, kind) in enumerate(path['events']):
        seg = int(np.searchsorted(path['taus'], te, side='left'))          # the segment before the event
        w_prev = path['W'][max(seg - 1, 0)]
        udot = abs(Phi[:, j] @ (y - Phi @ w_prev)) / N
        wj = np.abs(det['W_end'][j])
        i_half = np.flatnonzero(np.abs(det['Wd'][:, j]) >= 0.5 * wj)
        Dj = det['D'][i_half[0], j] if i_half.size else det['D'][-1, j]
        pphi = P[max(seg - 1, 0)] @ Phi[:, j]
        out.append((te, 3 * kap / udot, 5 * N / (Lam * Dj * (pphi @ pphi))))
    return out


def part_C(pool):
    jobs = [(e, 0.0, 0) for e in EPS_C] + [(e, s, k) for e in EPS_C for s in SIGS_C for k in range(PATHS_C)]
    res = pool.map(c_job, jobs)
    det = {r['eps']: r for r in res if r['sig'] == 0}
    rows = []
    for e, dd in det.items():                                              # the noiseless run's coefficients, for the windows
        dd['Wd'], dd['W_end'] = dd['W'], dd['w_end']
    for r in res:
        if r['sig'] == 0:
            continue
        dd = det[r['eps']]
        L, Phi, y, b, rho, path, q = setup(SEED_CD, r['eps'])
        N = Phi.shape[0]
        kap = r['sig'] / (N * np.sqrt(r['Lam']))
        dU = r['U'] - dd['U']
        nu_lin = r['nu_lin']
        dev = dU - nu_lin @ Phi
        err = np.abs(dev).max(1)
        scale = np.abs(dU).max()
        taus = r['taus']
        ev_path = [t for t, i, kind in path['events']]
        ev_run = [t for t in dd['entries'] if np.isfinite(t)]
        dist = lambda evs: np.min(np.abs(taus[:, None] - np.array(evs)[None, :]), axis=1)
        out = {w: (float(err[dist(ev_path + ev_run) >= w].max() / scale) if np.any(dist(ev_path + ev_run) >= w) else np.nan) for w in WIN_C}
        wins = windows_C(path, Phi, y, dd, kap, r['Lam'])
        inside = np.zeros(taus.size, bool)
        mech = {}
        for f in WIN_FAC_C:
            keep = np.ones(taus.size, bool)
            for te, wn, wa in wins:
                keep &= np.abs(taus - te) >= f * max(wn, wa)
            mech[f] = float(err[keep].max() / scale) if keep.any() else np.nan
            if f == WIN_FAC_C[0]:
                inside = ~keep
        keep2 = np.ones(taus.size, bool)
        for te, wn, wa in wins:
            keep2 &= np.abs(taus - te) >= WIN_FAC_C[1] * max(wn, wa)
        # the residual outside the windows against the active block's Ornstein–Uhlenbeck component
        R = (r['nu'] - dd['nu']) - nu_lin
        seg = np.clip(np.searchsorted(path['taus'], taus, side='right') - 1, 0, len(path['W']) - 1)
        pred = np.array([DL.normal_ou_var(Phi, dd['D'][i], path['W'][seg[i]] != 0, kap, r['Lam']) for i in range(taus.size)])
        ou = float(np.sqrt(np.mean((R[keep2] ** 2).sum(1)) / np.mean(pred[keep2]))) if keep2.any() else np.nan
        en_r, en_d = r['entries'], dd['entries']
        row = dict(eps=r['eps'], sig=r['sig'], k=r['k'], share=float(scale), kappa=float(kap), lin_err=float(err.max() / scale), lin_out=out,
                   lin_mech=mech, lin_in=float(err[inside].max() / scale) if inside.any() else 0.0,
                   widths=[(float(wn), float(wa)) for te, wn, wa in wins], ou_ratio=ou,
                   ou_rms=float(np.sqrt(np.mean((R[keep2] ** 2).sum(1)))) if keep2.any() else np.nan, lin_l2=float(np.sqrt(np.mean(dev ** 2) / np.mean(dU ** 2))),
                   iss=r['iss'], iss_det=dd['iss'],
                   order=[i for i in np.argsort(en_r) if np.isfinite(en_r[i])] == [i for i in np.argsort(en_d) if np.isfinite(en_d[i])],
                   dt_entry=float(np.nanmax(np.abs(np.where(np.isfinite(en_r), en_r, np.nan) - np.where(np.isfinite(en_d), en_d, np.nan)))),
                   run_vs_path=float(np.nanmax(np.abs(np.sort(np.array(ev_run)) - np.sort(np.array(ev_path))))) if len(ev_run) == len(ev_path)
                   else np.nan, end=float(np.abs(r['w_end'] - dd['w_end']).max()))
        rows.append(row)
        print(f"  C ε {r['eps']:g} σ {r['sig']} path {r['k']}: share {row['share']:.2e} = {row['share'] / kap:.2f}κ, linear response "
              f"error {row['lin_err']:.3f} (outside the fixed windows {', '.join(f'{v:.3f}' for v in out.values())}; the mechanisms' "
              f"{', '.join(f'{v:.3f}' for v in mech.values())}; inside {row['lin_in']:.3f}), OU {ou:.3f} ({r['secs']:.0f}s)", flush=True)
    return dict(rows=rows, det={e: dict(iss=det[e]['iss'], entries=det[e]['entries']) for e in det})


# --- D ---------------------------------------------------------------------------------------------------------------------------------
@cached
def d_job(args):
    eps, sig, k = args
    warnings.filterwarnings('ignore')
    L, Phi, y, b, rho, path, q = setup(SEED_CD, eps)
    N, d = Phi.shape
    S = np.flatnonzero(path['w_end'])
    off = np.setdiff1d(np.arange(d), S)
    F, z0, al0, w0 = plateau(L, Phi, y, q, path, eps)
    Lam = RT.transit_scale(L, eps)
    t0 = time.time()
    ts, Z, Ws, al = DL.mirror_update(F, Phi, y, sig, H_D, TW_D, seed=7 + 1000 * k + int(round(100 * sig)), rec_every=1, zeta0=z0, al0=al0)
    Ws = Ws[len(Ws) // 20:]
    th = 1e-3 * np.abs(path['w_end']).max()
    return dict(eps=eps, sig=sig, k=k, bias=Ws[:, S].mean(0) - path['w_end'][S], se=batch_se(Ws[:, S]), sd=Ws[:, S].std(0),
                inactive=float(np.abs(Ws[:, off]).max()), frac=float(np.mean(np.abs(Ws[:, off]).max(1) > th)),
                fit=float(np.abs(Ws.mean(0) @ Phi.T - y).max()), start=float(np.abs(w0 - path['w_end']).max()),
                crossing=float(Lam ** 2 * N ** 2 / sig ** 2), secs=time.time() - t0)


def part_D(pool):
    rows = pool.map(d_job, [(e, s, k) for e in EPS_D for s in SIGS_D for k in range(PATHS_D)])
    for r in rows:
        print(f"  D ε {r['eps']:g} σ {r['sig']} path {r['k']}: bias/SE {np.round(r['bias'] / r['se'], 2)}, inactive max {r['inactive']:.1e} "
              f"({r['secs']:.0f}s)", flush=True)
    return rows


# --- report ----------------------------------------------------------------------------------------------------------------------------
def fmt(v, p=3):
    return '(' + ', '.join(f'{x:.{p}g}' for x in np.atleast_1d(v)) + ')'


def report(A=None, B=None, C=None, D=None):
    s = ['# Proposition D‴.50: the rich transit under label noise — the dual is a Langevin diffusion', '',
         'Generated by `experiments/exp_transit_noise.py` (the designs of `exp_mirror_transit`; m = 12, G = 1, γ = 0.1, a_max = 1, '
         'c = 0.025; the label noise dζ = (y − Φw)dt/N + (σ/N)dW, T = σ²/(2N)).', '']
    if A is not None:
        L, Phi, y, b, q = small_dictionary()
        r = A['rwm 0.0']
        s += ['## A. The Gibbs law (a 4 × 8 dictionary with orthonormal rows, masses 2b_i, σ = 1, T = 0.125)', '',
              f"Random-walk Metropolis on π ∝ exp(−F/T) ({N_RWM} steps, acceptance {r['acc']:.2f}); the mirror update over t = {T_A:g} after "
              f"a burn-in of {BURN_A:.0%}; the linearised dual's rates at the projection {fmt(r['rates'])}. The coefficients' means and "
              'standard deviations, the largest deviation from the sampler\'s in units of the combined standard error (batch means) and '
              'relatively, the mean fit.', '',
              '| dynamics | E[w] | sd(w) | max_i |ΔE[w_i]|/se | max_i |Δsd_i|/sd_i | |E[Φw] − y|_∞ |', '|---|---|---|---|---|---|',
              f"| Metropolis (π) | {fmt(r['mean'])} | {fmt(r['sd'])} | — | — | {np.abs(r['fit']).max():.1e} |"]
        for h in H_A:
            x = A[f'mirror {h}']
            z = np.abs(x['mean'] - r['mean']) / np.sqrt(x['se'] ** 2 + r['se'] ** 2)
            s.append(f"| mirror update, h = {h} | {fmt(x['mean'])} | {fmt(x['sd'])} | {z.max():.1f} | {np.abs(x['sd'] / r['sd'] - 1).max():.3f} | "
                     f"{np.abs(x['fit']).max():.1e} |")
        s += ['', f"The noiseless projection: {fmt(r['w_proj'])}.", '',
              'The exact identity E_π[∇F∇Fᵀ] = T·E_π[∇²F], i.e. Cov_π(Φw) = T·ΦE_π[D]Φᵀ (D = diag dw/dη, E over the run): its relative '
              'Frobenius error on the mirror update\'s run (± its batch standard error), against the Laplace value T·ΦD(w*)Φᵀ at the mode: ' +
              '; '.join(f"h = {h}: {A[f'mirror {h}']['ident']:.3f} (± {A[f'mirror {h}']['ident_se']:.3f}) against "
                        f"{A[f'mirror {h}']['ident_mode']:.3f}" for h in H_A) + '.', '']
        if any(k.startswith('laplace') for k in A):
            x = [v for k, v in A.items() if k.startswith('laplace')][0]
            s += [f"(i)(b) at small T (σ = {SIG_A_SMALL}, T = {x['T']:.3g}; from the mode): Cov(w) against T·D^{{1/2}}PD^{{1/2}} (P of rank "
                  f"{x['rank']} in ℝ^{DA_A}) — relative Frobenius error {x['cov_rel']:.3f}, the standard deviations' ratios "
                  f"{fmt(x['sd_ratio'])}; the mean's shift from the projection in standard deviations {fmt(x['shift'], 2)} "
                  f"(the standard errors {fmt(x['se_rel'], 2)} of a standard deviation).", '']
        s += ['The paradigm\'s law read in Stratonovich\'s sense (stochastic Heun on α, β on its curve) against the mirror update with the '
              'same noise: the largest differences in w and η, and the conserved directions\' drift |P_{ker Φ}η|_∞ at the end (∝ h√t: the '
              'scheme\'s artefact; at a fixed step it has no stationary law) —', '',
              f'| h | t ≤ {T_PATH_A:g}: max|Δw|, max|Δη| | t ≤ {T_DRIFT_A:g}: max|Δw|, max|Δη| | |P_{{ker Φ}}η| at t = {T_PATH_A:g}, '
              f'{T_DRIFT_A:g} |', '|---|---|---|---|']
        for h in H_PATH_A:
            x = A[f'path {h}']
            s.append(f"| {h} | {x['short']['dw']:.1e}, {x['short']['deta']:.1e} | {x['long']['dw']:.1e}, {x['long']['deta']:.1e} | "
                     f"{x['short']['inv'][-1]:.1e}, {x['long']['inv'][-1]:.1e} |")
        x = A[f'ito {H_A[0]}']
        s += ['', f"Read in Itô's sense (Euler on (α, log β), h = {x['h']}): Q(t) − Q(0) at t = {x['t_end']:g} from {x['dq_end'].min():.3g} to "
              f"{x['dq_end'].max():.3g}; the rate over the first time unit against (γβ/2)(γa² − Gβa′(1 + a²/a_max²))σ²|φ|²/N² at the start: "
              f"ratios {fmt(x['rate1'] / x['drift0'], 3)} (one path over one time unit while the state moves: an illustration, the formula "
              f"itself is the unit test's); the coefficients' means {fmt(x['mean'])}, standard deviations {fmt(x['sd'])} (a non-stationary "
              "window).", '']
    if B is not None:
        s += ['## B. Fluctuation–dissipation on the support (ε = 10⁻⁶)', '',
              f'From the noiseless transit\'s end, over max({TW_B:g}, {NRELAX_B:g} relaxation times of the slowest mode) (the first 5 % dropped), '
              f'{PATHS_B} independent noise paths per design and σ: Var(w_i)/(T·E[dw_i/dη]) (E over the run; the exact identity\'s '
              'prediction) and Var(w_i)/(T·dw_i/dη at w_BP), as mean ± sd over the paths; E[dw/dη]/dw/dη(w_BP); the hypothesis\' size '
              '(T·D_i)^{1/2}/|w_i|; the correlations (0 predicted) against their sampling error (batch means), the largest |corr|/se; the '
              'mean\'s offset from weighted basis pursuit in standard errors (the largest over the paths); the inactive coefficients\' '
              'largest |w|; the window in relaxation times.', '',
              '| design | σ | h, t | support | Var/(T·E[D]) | Var/(T·D(w_BP)) | E[D]/D(w_BP) | (T·D)^{1/2}/\\|w\\| | max \\|corr\\| (max \\|corr\\|/se) | '
              'max \\|E[w] − w_BP\\|/se | max inactive \\|w\\| | relaxation times |',
              '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for sd_ in SEEDS_B:
            for g in SIGS_B:
                rr = [r for r in B if r['seed'] == sd_ and r['sig'] == g]
                if not rr:
                    continue
                ra = np.array([r['ratio'] for r in rr]); rb = np.array([r['ratio_bp'] for r in rr])
                pm = lambda a: '(' + ', '.join(f'{m:.2f} ± {v:.2f}' for m, v in zip(a.mean(0), a.std(0, ddof=1))) + ')'
                cz = max(float(np.max(np.abs(r['corr']) / r['corr_se'])) for r in rr) if rr[0]['corr'].size else 0.0
                cm = max(float(np.max(np.abs(r['corr']))) for r in rr) if rr[0]['corr'].size else 0.0
                bz = max(float(np.max(np.abs(r['bias']) / r['se'])) for r in rr)
                s.append(f"| {sd_} | {g} | {rr[0]['h']:g}, {rr[0]['tw']:.3g} | {rr[0]['S'].tolist()} | {pm(ra)} | {pm(rb)} | "
                         f"{fmt(np.mean([r['ED_rel'] for r in rr], 0), 3)} | {fmt(rr[0]['sd_rel'], 2)} | {cm:.3f} ({cz:.1f}) | {bz:.1f} | "
                         f"{max(r['inactive'] for r in rr):.1e} | {rr[0]['n_relax']:.0f} |")
        s.append('')
    if C is not None:
        rows = C['rows']
        s += ['## C. The noisy transit (the first design)', '',
              f'Per ε and σ over {PATHS_C} noise paths (each one Brownian path for every ε and σ — the finest ε\'s increments summed for '
              'the coarser): the noise\'s share sup_τ|u − u_det| (mean, and divided by κ = σ/(NΛ^{1/2})); the first-order response\'s '
              '(δν ← P_T(δν + κdW) along the path) relative error |δu − δu_lin|/sup|δu| — its sup over all τ; outside the fixed windows '
              f'±{WIN_C[0]:g}, ±{WIN_C[1]:g} around the path\'s and the noiseless run\'s events; outside the mechanisms\' windows '
              '(× ' + ', '.join(f'{f:g}' for f in WIN_FAC_C) + ' the larger of the noise\'s entry shift 3κ/|u̇_j| and the entrant\'s '
              'finite-ε absorption 5N/(Λ·D_j|P_{T⁻}φ_j|²)); inside the ×1 windows (σ-independent: finite-ε, linear in the noise); its L² '
              'over τ (the largest over the paths) —; the residual outside the ×2 windows over the active block\'s Ornstein–Uhlenbeck '
              'component, (E|δν − δν_lin|²/(κ²Σ1/(2Λλ_i)))^{1/2} (the mean over the paths); the entries\' order kept and their largest '
              'shift (interpolated between records); the noiseless run\'s sup|u_det − u_ISS| and its entries\' largest distance from the '
              'path\'s.', '',
              f'| ε | σ | κ | share (mean) | share/κ | sup | outside ±{WIN_C[0]:g}, ±{WIN_C[1]:g} | outside the mechanisms\' windows | inside | '
              'L² | OU ratio | widths: shift, absorption (max) | order kept | max entry shift | sup\\|u_det − u_ISS\\| | run − path |',
              '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
        for e in EPS_C:
            for g in SIGS_C:
                rr = [r for r in rows if r['eps'] == e and r['sig'] == g]
                if not rr:
                    continue
                sh = np.array([r['share'] for r in rr])
                wn = max(w[0] for r in rr for w in r['widths']); wa = max(w[1] for r in rr for w in r['widths'])
                s.append(f"| {e:g} | {g} | {rr[0]['kappa']:.2e} | {sh.mean():.2e} | {np.mean(sh / rr[0]['kappa']):.2f} | "
                         f"{max(r['lin_err'] for r in rr):.3f} | {max(r['lin_out'][WIN_C[0]] for r in rr):.3f}, "
                         f"{max(r['lin_out'][WIN_C[1]] for r in rr):.3f} | "
                         + ', '.join(f"{max(r['lin_mech'][f] for r in rr):.3f}" for f in WIN_FAC_C) + " | "
                         f"{max(r['lin_in'] for r in rr):.3f} | {max(r['lin_l2'] for r in rr):.3f} | {np.mean([r['ou_ratio'] for r in rr]):.3f} | "
                         f"{wn:.2e}, {wa:.2e} | {sum(r['order'] for r in rr)}/{len(rr)} | {max(r['dt_entry'] for r in rr):.2e} | "
                         f"{rr[0]['iss_det']:.2e} | {rr[0]['run_vs_path']:.3f} |")
        # the share's slope in ε: per path (one Brownian path for every ε), mean ± sd over the paths
        sl = {}
        for g in SIGS_C:
            per = []
            for k in range(PATHS_C):
                v = [next((r['share'] for r in rows if r['eps'] == e and r['sig'] == g and r['k'] == k), np.nan) for e in EPS_C]
                if np.all(np.isfinite(v)):
                    per.append(np.polyfit(np.log(EPS_C), np.log(v), 1)[0])
            sl[g] = (np.mean(per), np.std(per, ddof=1) if len(per) > 1 else np.nan)
        s += ['', 'Log-log slope of the share in ε, per path (mean ± sd over the paths): ' +
              ', '.join(f'σ = {g}: {sl[g][0]:.3f} ± {sl[g][1]:.3f}' for g in SIGS_C) + ' (κ ∝ ε^{1/4}).', '']
    if D is not None:
        s += ['## D. The local quasi-equilibrium in the rich regime (the first design)', '',
              f'From the noiseless transit\'s end, a window t = {TW_D:g} on {PATHS_D} independent noise paths — a fraction of the face\'s '
              'crossing time Λ²N²/σ² (the column "window/crossing"), so not a sample of π: the support\'s mean against weighted basis pursuit '
              'in standard errors (mean ± sd over the paths of the largest |z|, and its largest), the standard deviations, the inactive '
              'coefficients (largest |w| and the fraction of the time any exceeds 10⁻³max|w_BP|), the time-averaged fit.', '',
              '| ε | σ | window/crossing | max\\|z\\| per path: mean ± sd (max) | sd(w_S) | max inactive \\|w\\| | fraction | \\|Φw̄ − y\\|_∞ |',
              '|---|---|---|---|---|---|---|---|']
        for e in EPS_D:
            for g in SIGS_D:
                rr = [r for r in D if r['eps'] == e and r['sig'] == g]
                if not rr:
                    continue
                zz = np.array([np.max(np.abs(r['bias']) / r['se']) for r in rr])
                s.append(f"| {e:g} | {g} | {TW_D / rr[0]['crossing']:.1e} | {zz.mean():.2f} ± {zz.std(ddof=1) if zz.size > 1 else 0:.2f} "
                         f"({zz.max():.2f}) | {fmt(np.mean([r['sd'] for r in rr], 0))} | {max(r['inactive'] for r in rr):.1e} | "
                         f"{max(r['frac'] for r in rr):.3f} | {max(r['fit'] for r in rr):.1e} |")
        s.append('')
    return '\n'.join(s)


def main():
    raw = os.path.join(OUT, 'transit_noise_raw.npz')
    md = os.path.join(OUT, 'transit_noise.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCD']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCD')
    res = dict(old)
    t0 = time.time()
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D)):
            if k in parts:
                res[k] = fn(pool); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCD']))


if __name__ == '__main__':
    main()
