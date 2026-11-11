"""
Proposition D‴.43: the colour of the measurement noise selects the kernels — from the trace (a weighted ℓ¹ norm) to the log-determinant
(a support count)
=====================================================================================================================================

A kernel dictionary sampled at N < d points (the design Φ), the coefficients w = u⊙u − v⊙v (the unbalanced Fisher–Rao law of a signed
kernel measure's two parts on a fixed dictionary), fitted by gradient descent with label noise; along the manifold M of exact fits the
parameters follow the reduced flow of `colour_selection` (R = tr Ψ(K): the trace for white draws, (1/τ_c)log det(I + τ_cK) for a slow
AR(1) noise).

A. The potential's identities at random points: det(I + τK) against its subset (Cauchy–Binet) expansion; ∂R/∂D_i = (4/N)ω_i; R's second
   differences along random segments in D (zero for white draws, negative for coloured ones: concave); on a k-sparse exact fit
   τ_cR − k log(4τ_c/N) → log(det(Φ_SᵀΦ_S)Π_{i∈S}|w_i|) as τ_c grows (the support's size times log τ_c).
B. Label-noise gradient descent (P_B paths, step H_B, per-step variance σ_B²) against the reduced flow with the scheme's filter at the step,
   on N_B × D_B Gaussian designs with a K_B-sparse truth, the first eight seeds: fresh draws and AR(1) draws of correlation time τ_c = 100
   (in time units); the paths' mean against the flow at nine times of the rescaled time s = 4v_ht/N, and the interpolant the noiseless
   gradient flow reaches from each path's end state — the fraction of paths recovering the truth, against the flow's own success and
   basis pursuit's; and, per seed, that fraction at σ = 0.3 and 0.5 (the noise's size: the reduced flow is its small-σ limit). The
   start: the noiseless flow's interpolant from u = v = α (exactly: the hyperbolic entropy's minimiser on M, P = α²).
C. Recovery against sparsity (the reduced flows, continuum filters): the fraction of instances whose flow ends within 10⁻⁴ (relative)
   of the truth, for weighted basis pursuit, the white flow (run to convergence), slow flows (τ_c = 10, 100, 1000), and the white flow
   followed by the slow one (τ_c = 100) — N × d = 10 × 40 and 20 × 80; where a slow flow fails, whether its end has R below the truth's.
D. Gaussian kernel dictionaries on [0, 1] (D_D centres, N_D equispaced samples, two widths — neighbouring atoms nearly collinear),
   random k-atom truths (signs ±1, magnitudes in [0.5, 1.5]): the recovery rates of basis pursuit and of the white and slow flows, and,
   on the first instances where basis pursuit fails, label-noise descent (fresh and AR(1) draws, τ_c = 100) from the same start, with
   the reduced flow at the same s and run to convergence.
E. Other colours at the step H_B: AR(1), a smooth noise (Gaussian step correlations) and a draw kept k steps, at comparable
   correlation times, against fresh draws (the reduced flow with the scheme's filter), on part C's 10 × 40 instances at k = 3, 4; a smooth
   noise at τ_c = 100 over the slow horizon and ten times it (its selection runs at (λτ_c)⁻²).
F. The start: the gradient flow's interpolant from u = v = α for α = 0.1, 0.5, 2, the flows run to convergence (s ≤ S_MAX), on part C's
   10 × 40 instances at k = 3, 4.
G. References on part C's 10 × 40 instances: reweighted ℓ¹ (Candès–Wakin–Boyd: weights 1/(|w_i| + ε), ε = 0.01, 0.1, 1 times max|w*|, from
   basis pursuit's solution and from the flows' start), the majorisation–minimisation of the slow noise's own potential (weighted ℓ¹ with
   the weights ω(|w|), τ_c = 100 — Fazel–Hindi–Boyd's log-det heuristic) and sparse Bayesian learning (Tipping's EM), against the slow
   flow: the same kind of penalty minimised by other algorithms instead of descended by the noise.
H. Which exact fits the selection holds: at a k-sparse exact fit (P = 0, support S) the inactive atoms decay when the weighted dual
   certificate |c_i| < ω_i holds, c = Φ_Iᵀ(Φ_S⁺)ᵀ(ω_S⊙sign w_S) (sufficient, not necessary; for white draws Fuchs's least-squares
   certificate); the largest |c_i|/ω_i at the truth and at basis pursuit's solution, for white draws and slow ones (τ_c = 10, 100, 1000),
   on part C's 10 × 40 instances, beside basis pursuit's success; and on kernel dictionaries (widths 0.03, 0.05, 0.08) for τ_c up to 10⁵,
   beside whether the slow flow returns to the truth from P = 10⁻³.
I. A fit floor: y = Φw* + ε·e — whether the slow flow's and basis pursuit's k largest atoms are the truth's support.
J. The selection is the explicit (Euler) scheme's: label-noise descent at the same physical noise and time with Euler at the step H_B and
   H_B/4 and with Heun's scheme (the same draw in both stages) — the mean change of log P on the atoms off the truth's support.
K. Statistics and the competing fits (referee round 142). (K1) From Parts C and G's stored outcomes, no new run: exact (Clopper–Pearson)
   95% intervals of the recovery rates, and the paired comparisons on the same instances — the discordant counts and the exact
   (binomial) McNemar p-value — of the slow flow (τ_c = 100) against basis pursuit, sparse Bayesian learning and the slow flow at
   τ_c = 1000. (K2) Where the slow flow (τ_c = 100) fails on Part C's 10 × 40 instances (k = 3 … 6), the flow from the same start run to
   convergence (s ≤ S_MAX): whether it agrees with Part C's fixed-horizon outcome, its end's support and whether it is a basic exact fit
   (at most N atoms, Φ_S of full column rank), the certificate max_{i∉S}|c_i|/ω_i at the end for τ_c = 100 … 10⁵ and its τ_c → ∞ limit at
   an N-atom end (|Σ_jα_ij/w_j| over Σ_jα_ij²/|w_j|, α_i = Φ_S⁻¹φ_i: at an N-atom fit both sides of the certificate are O(1/τ_c) and the
   ratio tends to this τ_c-free number), and R at the end against R at the truth as τ_c varies (τ_c = 1 … 10⁵).

Output: results/chaos/colour_selection.md and colour_selection_raw.npz (CS_WORKERS; `--from-raw` regenerates the md; `--only A,B,…`;
`--only K` needs Parts C and G in the raw file).
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
from src.core import colour_selection as CS
from src.core.runcache import cached
import exp_cohesion_leak as E

OUT = os.environ.get('CS_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('CS_WORKERS', 4))
ALPHA, TOL = 0.5, 1e-4                                                        # the start u = v = α; success within TOL·max|w*|
N_B, D_B, K_B, H_B, SIG_B, P_B, SEEDS_B = 10, 40, 4, 0.5, 0.1, 16, tuple(range(8))
TAUS_B, S_B = (0.0, 100.0), {0.0: 300.0, 10.0: 150.0, 100.0: 400.0}            # the rescaled time each colour runs
SIGMAS_B = (0.3, 0.5)                                                         # B's σ sweep (AR(1), τ_c = 100)
SIZES_C, KS_C, N_INST_C = ((10, 40), (20, 80)), {(10, 40): (1, 2, 3, 4, 5, 6), (20, 80): (2, 4, 6, 8, 10)}, {(10, 40): 40, (20, 80): 20}
TAUS_C, S_WHITE, S_SLOW = (10.0, 100.0, 1000.0), 2000.0, lambda tc: 400.0 + 60.0 * np.sqrt(1 + tc)
D_D, N_D, WIDTHS_D, KS_D, N_INST_D, N_NOISY_D = 40, 20, (0.03, 0.05), (3, 4, 5, 6), 40, 3
KS_E, N_INST_E, TC_E = (3, 4), 20, 10.0
ALPHAS_F = (0.1, 0.5, 2.0)
S_MAX = 1e5                                                                   # the convergence criterion's cap (`flow_to_convergence`)
KS_G, N_INST_G, EPS_G, IT_G, TC_G = (2, 3, 4, 5, 6), 40, 0.1, 12, 100.0
EPS_GS = (0.01, 0.1, 1.0)                                                     # reweighted ℓ¹'s ε, relative to max|w*|
TAUS_H = (0.0, 10.0, 100.0, 1000.0)
WIDTHS_H, KS_H, N_INST_H, TAUS_HK = (0.03, 0.05, 0.08), (4, 6), 20, (1e2, 1e3, 1e4, 1e5)
EPS_I, KS_I, N_INST_I = (1e-3, 1e-2, 5e-2), (3, 4), 20                 # the fit floor y = Φw* + ε·e
T_J, TC_J, P_J = 2e5, 50.0, 16                                          # Euler against Heun: physical time, OU time, paths
KS_K, TAUS_KC, TAUS_KR = (3, 4, 5, 6), (100.0, 1e3, 1e4, 1e5), (1.0, 10.0, 100.0, 1e3, 1e4, 1e5)   # K2: sparsities, certificate and R's τ_c


def instance(N, d, k, seed):
    """A Gaussian design (columns of norm ≈ 1) and a k-sparse truth with signs ±1 and magnitudes in [0.5, 1.5]."""
    rng = np.random.default_rng(seed)
    Phi = rng.standard_normal((N, d)) / np.sqrt(N)
    ws = np.zeros(d); idx = rng.choice(d, k, replace=False); ws[idx] = rng.choice([-1, 1], k) * rng.uniform(0.5, 1.5, k)
    return Phi, ws


def start(Phi, y, alpha):
    """The noiseless gradient flow's interpolant from u = v = α (exactly, `colour_selection.flow_interpolant`): P = α²."""
    d = Phi.shape[1]
    return CS.flow_interpolant(Phi, y, np.full(d, alpha), np.full(d, alpha))


def seed_of(N, k, i):
    return 100000 * N + 1000 * k + i


def ok(w, ws):
    return bool(np.abs(w - ws).max() < TOL * np.abs(ws).max())


def corr_ct(tc):
    return ('fresh',) if tc == 0 else ('ou', tc)


def part_A():
    rng = np.random.default_rng(1)
    rows = []
    for N, d in ((4, 7), (5, 9)):
        Phi = rng.standard_normal((N, d)) / np.sqrt(N); D = rng.uniform(0.05, 1.5, d)
        for tau in (0.3, 10.0, 1e3):
            K = CS.gram(Phi, D)
            exact = np.linalg.det(np.eye(N) + tau * K)
            e = 1e-6
            g = np.array([(CS.potential(Phi, D + e * np.eye(d)[i], 0, ('ou', tau)) - CS.potential(Phi, D - e * np.eye(d)[i], 0, ('ou', tau)))
                          / (2 * e) for i in range(d)])
            sec = []
            for _ in range(20):                                               # second differences along random segments
                D1, D2 = rng.uniform(0.05, 1.5, d), rng.uniform(0.05, 1.5, d)
                f = [CS.potential(Phi, D1 + t * (D2 - D1), 0, ('ou', tau)) for t in (0.0, 0.5, 1.0)]
                sec.append(f[0] + f[2] - 2 * f[1])
            fw = [CS.potential(Phi, D1 + t * (D2 - D1), 0, ('fresh',)) for t in (0.0, 0.5, 1.0)]
            rows.append(dict(N=N, d=d, tau=tau, cb=abs(CS.det_subsets(Phi, D, tau) / exact - 1),
                             grad=float(np.abs(g - (4 / N) * CS.weights(Phi, D, 0, ('ou', tau))).max() / np.abs(g).max()),
                             sec_max=float(max(sec)), sec_min=float(min(sec)), sec_white=float(abs(fw[0] + fw[2] - 2 * fw[1]))))
    lim = []
    for N, d, k in ((10, 40, 3), (10, 40, 6), (20, 80, 5)):
        Phi, ws = instance(N, d, k, 7 + k)
        S = np.flatnonzero(ws)
        ref = np.log(np.linalg.det(Phi[:, S].T @ Phi[:, S]) * np.prod(np.abs(ws[S])))
        lim.append(dict(N=N, d=d, k=k, ref=ref, vals=[(tc, tc * CS.potential(Phi, np.abs(ws), 0, ('ou', tc)) - k * np.log(4 * tc / N))
                                                      for tc in (1e2, 1e4, 1e6)]))
    return dict(rows=rows, lim=lim)


@cached
def b_job(args):
    """Part B, one (seed, τ_c, σ): the paths' mean against the reduced flow at nine times, the end interpolants, basis pursuit."""
    seed, tc, sig = args
    t0 = time.time()
    Phi, ws = instance(N_B, D_B, K_B, seed)
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    u0, v0 = CS.uv_of(w0, P0)
    corr = ('fresh',) if tc == 0 else ('ar1', float(np.exp(-H_B / tc)))
    rho = 0.0 if tc == 0 else float(np.exp(-H_B / tc))
    vh = H_B * sig ** 2 / N_B
    S = S_B[tc]
    n = int(round(S * N_B / (4 * vh) / H_B)); every = n // 8
    s_rec = np.arange(9) * every * H_B * 4 * vh / N_B
    _, Wf, _ = CS.reduced_flow(Phi, w0, P0, s_rec[-1], H_B, corr, t_eval=s_rec)
    u, v, rec = CS.label_noise_gd(Phi, y, u0, v0, H_B, sig, rho, n, np.random.default_rng(seed + 77), P_B, rec_every=every)
    ends, conv = CS.end_interpolants(Phi, y, u, v)
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    print(f"[{time.time() - t0:.0f}s] B seed {seed} tau_c {tc:g} sigma {sig:g}", flush=True)
    return dict(seed=seed, tc=tc, sig=sig, s=s_rec, dev=np.abs(rec - Wf.T).max(axis=1), flow_err=np.abs(Wf.T - ws[None, :]).max(axis=1),
                path_err=np.abs(rec - ws[None, :]).max(axis=1), end_err=np.abs(ends - ws[None, :]).max(axis=1),
                end_nnz=np.where(conv, (np.abs(ends) > 1e-3).sum(axis=1), -1), bp_err=float(np.abs(wbp - ws).max()), n=n, vh=vh,
                wmax=float(np.abs(ws).max()), raw_err=np.abs(u * u - v * v - ws[None, :]).max(axis=1))


@cached
def c_job(args):
    """Part C, one instance: success of basis pursuit, the white flow, the slow flows and the annealed one; the supports' sizes."""
    N, d, k, i = args
    Phi, ws = instance(N, d, k, seed_of(N, k, i))
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    nnz = lambda w: int((np.abs(w) > 1e-3 * np.abs(ws).max()).sum())
    out = {}
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0)); out['bp'] = (ok(wbp, ws), nnz(wbp))
    Ww, Pw, s_w, conv_w = CS.flow_to_convergence(Phi, w0, P0, S_WHITE, 0, ('fresh',), S_max=S_MAX)
    out['white'] = (ok(Ww, ws), nnz(Ww))
    R = {}
    for tc in TAUS_C:
        _, Wt, Pt = CS.reduced_flow(Phi, w0, P0, S_SLOW(tc), 0, ('ou', tc)); out[f'slow {tc:g}'] = (ok(Wt[:, -1], ws), nnz(Wt[:, -1]))
        D_end = np.sqrt(Wt[:, -1] ** 2 + 4 * Pt[:, -1] ** 2)
        R[tc] = (CS.potential(Phi, D_end, 0, ('ou', tc)), CS.potential(Phi, np.abs(ws), 0, ('ou', tc)))
    _, Wa, _ = CS.reduced_flow(Phi, Ww, np.maximum(Pw, 1e-300), S_SLOW(100.0), 0, ('ou', 100.0))
    out['white, then slow 100'] = (ok(Wa[:, -1], ws), nnz(Wa[:, -1]))
    return dict(N=N, d=d, k=k, i=i, out=out, white_conv=conv_w, white_s=s_w, white_bp=float(np.abs(Ww - wbp).max()), R=R)


def dict_instance(width, k, i):
    """A Gaussian kernel dictionary (D_D centres on [0, 1], N_D equispaced samples, unit columns) and a random k-atom truth."""
    Phi = CS.rbf_dictionary(np.linspace(0.0, 1.0, N_D), np.linspace(0.0, 1.0, D_D), width)
    rng = np.random.default_rng(9000 + 100 * k + i)
    ws = np.zeros(D_D); idx = rng.choice(D_D, k, replace=False); ws[idx] = rng.choice([-1, 1], k) * rng.uniform(0.5, 1.5, k)
    return Phi, ws


@cached
def d_job(args):
    """Part D, one (width, k, instance): success of basis pursuit and of the white and slow flows; the supports' sizes."""
    width, k, i = args
    Phi, ws = dict_instance(width, k, i)
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    nnz = lambda w: int((np.abs(w) > 1e-3 * np.abs(ws).max()).sum())
    out = {}
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0)); out['bp'] = (ok(wbp, ws), nnz(wbp))
    _, Ww, _ = CS.reduced_flow(Phi, w0, P0, S_WHITE, 0, ('fresh',)); out['white'] = (ok(Ww[:, -1], ws), nnz(Ww[:, -1]))
    for tc in (10.0, 100.0):
        _, Wt, _ = CS.reduced_flow(Phi, w0, P0, S_SLOW(tc), 0, ('ou', tc)); out[f'slow {tc:g}'] = (ok(Wt[:, -1], ws), nnz(Wt[:, -1]))
    return dict(width=width, k=k, i=i, out=out)


@cached
def d_noisy_job(args):
    """Part D, label-noise descent on one instance (fresh draws and AR(1) τ_c = 100, P_B paths, step H_B, σ_B): the end interpolants'
    errors and supports, against the reduced flows'."""
    width, k, i = args
    t0 = time.time()
    Phi, ws = dict_instance(width, k, i)
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    u0, v0 = CS.uv_of(w0, P0)
    lam = np.linalg.eigvalsh(CS.gram(Phi, np.sqrt(w0 ** 2 + 4 * P0 ** 2))).max()
    h = min(H_B, 0.5 / lam)                                                   # hλ_max ≤ 1/2 at the start
    vh = h * SIG_B ** 2 / N_D
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    out = dict(width=width, k=k, i=i, h=h)
    for tc in (0.0, 100.0):
        rho = 0.0 if tc == 0 else float(np.exp(-h / tc))
        S = S_B[tc]
        n = int(round(S * N_D / (4 * vh) / h))
        u, v, _ = CS.label_noise_gd(Phi, y, u0, v0, h, SIG_B, rho, n, np.random.default_rng(1000 * k + i), P_B)
        ends, conv = CS.end_interpolants(Phi, y, u, v)
        raw = u * u - v * v
        corr = ('fresh',) if tc == 0 else ('ar1', rho)
        flow = CS.reduced_flow(Phi, w0, P0, S, h, corr)[1][:, -1]
        wc, _, s_c, conv_c = CS.flow_to_convergence(Phi, w0, P0, S, h, corr, S_max=S_MAX)
        out[tc] = dict(n=n, err=np.abs(ends - ws[None, :]).max(axis=1), conv_err=float(np.abs(wc - ws).max()), conv_s=s_c, conv=conv_c,
                       conv_nnz=int((np.abs(wc) > 1e-3 * np.abs(ws).max()).sum()),
                       nnz=np.where(conv, (np.abs(ends) > 1e-3 * np.abs(ws).max()).sum(axis=1), -1),
                       raw_err=np.abs(raw - ws[None, :]).max(axis=1), raw_bp=np.abs(raw - wbp[None, :]).max(axis=1),
                       alive=((u * u + v * v) > 1e-8 * (u * u + v * v).max(axis=1, keepdims=True)).sum(axis=1),
                       flow_err=float(np.abs(flow - ws).max()), flow_nnz=int((np.abs(flow) > 1e-3 * np.abs(ws).max()).sum()),
                       flow_bp=float(np.abs(flow - wbp).max()))
    print(f"[{time.time() - t0:.0f}s] D noisy width {width:g} k {k} instance {i}", flush=True)
    return out


def d_noisy_cases():
    """The first N_NOISY_D instances of width WIDTHS_D[1], k = 4 where basis pursuit fails."""
    cases, i = [], 0
    while len(cases) < N_NOISY_D and i < N_INST_D:
        Phi, ws = dict_instance(WIDTHS_D[1], 4, i)
        if not ok(CS.basis_pursuit(Phi, Phi @ ws, (Phi ** 2).sum(axis=0)), ws):
            cases.append((WIDTHS_D[1], 4, i))
        i += 1
    return cases


def e_corrs():
    k = int(round(TC_E / H_B))
    return {'fresh': ('fresh',), f'AR(1), τ_c = {TC_E:g}': ('ar1', float(np.exp(-H_B / TC_E))),
            f'smooth (Gaussian), τ_c = {TC_E:g}': ('gauss', TC_E / H_B), f'kept {2 * k} steps': ('block', 2 * k)}


@cached
def e_job(args):
    """Part E, one instance: success of the reduced flow with each colour's filter at the step H_B."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    out = {}
    for name, corr in e_corrs().items():
        _, W, _ = CS.reduced_flow(Phi, w0, P0, S_WHITE if corr[0] == 'fresh' else S_SLOW(TC_E), H_B, corr)
        out[name] = ok(W[:, -1], ws)
    for mult in (1, 10):                                                      # a smooth noise selects at (λτ_c)⁻²: its horizon
        _, W, _ = CS.reduced_flow(Phi, w0, P0, mult * S_SLOW(100.0), H_B, ('gauss', 100.0 / H_B))
        out[f'smooth, τ_c = 100, {mult}× the slow horizon'] = ok(W[:, -1], ws)
    out['basis pursuit'] = ok(CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0)), ws)
    return dict(k=k, i=i, out=out)


@cached
def f_job(args):
    """Part F, one instance: success of the white and slow (τ_c = 100) flows from the interpolants of u = v = α."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    out, conv = {}, {}
    for a in ALPHAS_F:
        w0, P0 = start(Phi, y, a)
        Ww, _, sw, cw = CS.flow_to_convergence(Phi, w0, P0, S_WHITE, 0, ('fresh',), S_max=S_MAX)
        Wt, _, st, ct = CS.flow_to_convergence(Phi, w0, P0, S_SLOW(100.0), 0, ('ou', 100.0), S_max=S_MAX)
        out[a] = (ok(w0, ws), ok(Ww, ws), ok(Wt, ws))
        conv[a] = (cw, sw, ct, st)
    out['bp'] = ok(CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0)), ws)
    return dict(k=k, i=i, out=out, conv=conv)


def mm_logdet(Phi, y, w, tc, iters=IT_G):
    """The majorisation–minimisation of R = (1/τ_c)log det(I + τ_cK) on M at P = 0: w ← argmin Σ_iω_i(|w|)|w_i| subject to Φw = y."""
    for _ in range(iters):
        w = CS.basis_pursuit(Phi, y, CS.weights(Phi, np.abs(w) + 1e-12, 0, ('ou', tc)))
    return w


def reweighted_l1(Phi, y, w, eps=EPS_G, iters=IT_G):
    """Candès–Wakin–Boyd: w ← argmin Σ_i|w_i|/(|w_i^{old}| + ε) subject to Φw = y."""
    for _ in range(iters):
        w = CS.basis_pursuit(Phi, y, 1.0 / (np.abs(w) + eps))
    return w


@cached
def g_job(args):
    """Part G, one instance: success of reweighted ℓ¹ and of the log-det's majorisation–minimisation (from basis pursuit and from the
    flows' start), with their support sizes."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    w0, _ = start(Phi, y, ALPHA)
    nnz = lambda w: int((np.abs(w) > 1e-3 * np.abs(ws).max()).sum())
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    out = {'bp': (ok(wbp, ws), nnz(wbp))}
    for e_ in EPS_GS:
        eps = e_ * np.abs(ws).max()
        w = reweighted_l1(Phi, y, wbp, eps); out[f'reweighted ℓ¹ (ε = {e_:g}), from basis pursuit'] = (ok(w, ws), nnz(w))
        w = reweighted_l1(Phi, y, w0, eps); out[f'reweighted ℓ¹ (ε = {e_:g}), from the start'] = (ok(w, ws), nnz(w))
    w = mm_logdet(Phi, y, wbp, TC_G); out['log-det MM from basis pursuit'] = (ok(w, ws), nnz(w))
    w = mm_logdet(Phi, y, w0, TC_G); out['log-det MM from the start'] = (ok(w, ws), nnz(w))
    w = CS.sbl_em(Phi, y); out['sparse Bayesian learning (EM)'] = (ok(w, ws), nnz(w))
    return dict(k=k, i=i, out=out)


def certificate(Phi, w, tc):
    """The largest |c_i|/ω_i over the atoms off w's support, c = Φ_Iᵀ(Φ_S⁺)ᵀ(ω_S⊙sign w_S), ω at D = |w| (white: tc = 0)."""
    S = np.flatnonzero(np.abs(w) > 1e-9 * np.abs(w).max()); I = np.setdiff1d(np.arange(Phi.shape[1]), S)
    om = CS.weights(Phi, np.abs(w) + 1e-300, 0, corr_ct(tc))
    c = Phi[:, I].T @ np.linalg.pinv(Phi[:, S]).T @ (om[S] * np.sign(w[S]))
    return float(np.max(np.abs(c) / om[I])), len(S)


@cached
def h_job(args):
    """Part H, one instance: the certificate ratios at the truth and at basis pursuit's solution, for each colour; basis pursuit's success."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    wbp = np.where(np.abs(wbp) > 1e-9 * np.abs(wbp).max(), wbp, 0.0)
    return dict(k=k, i=i, bp=ok(wbp, ws), truth={tc: certificate(Phi, ws, tc)[0] for tc in TAUS_H},
                at_bp={tc: certificate(Phi, wbp, tc)[0] for tc in TAUS_H}, bp_nnz=certificate(Phi, wbp, 0.0)[1])


@cached
def hk_job(args):
    """Part H on kernel dictionaries, one (width, k, instance): the certificate at the truth for τ_c up to 10⁵, basis pursuit's success,
    and whether the slow flow (τ_c = 100, 1000) started at the truth with every P = 10⁻³ returns to it (s ≤ 2·10⁴)."""
    width, k, i = args
    Phi, ws = dict_instance(width, k, i)
    y = Phi @ ws
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    hold = {}
    for tc in (100.0, 1000.0):
        w_, _, _, _ = CS.flow_to_convergence(Phi, ws, np.full(D_D, 1e-3), 1e3, 0, ('ou', tc), S_max=2e4)
        hold[tc] = ok(w_, ws)
    return dict(width=width, k=k, i=i, bp=ok(wbp, ws), truth={tc: certificate(Phi, ws, tc)[0] for tc in (0.0,) + TAUS_HK}, hold=hold)


@cached
def i_job(args):
    """Part I, one (ε, k, instance): the fit floor y = Φw* + ε·e (|e| = 1): the slow flow's (τ_c = 100) and basis pursuit's k largest atoms
    against the truth's support."""
    eps, k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    e = np.random.default_rng(31 + 100 * k + i).standard_normal(10); e /= np.linalg.norm(e)
    y = Phi @ ws + eps * e
    w0, P0 = start(Phi, y, ALPHA)
    _, Wt, _ = CS.reduced_flow(Phi, w0, P0, S_SLOW(100.0), 0, ('ou', 100.0))
    wbp = CS.basis_pursuit(Phi, y, (Phi ** 2).sum(axis=0))
    S = set(np.flatnonzero(ws))
    top = lambda w: set(np.argsort(-np.abs(w))[:k])
    nnz = lambda w: int((np.abs(w) > 1e-3 * np.abs(ws).max()).sum())
    return dict(eps=eps, k=k, i=i, slow=top(Wt[:, -1]) == S, bp=top(wbp) == S, slow_nnz=nnz(Wt[:, -1]), bp_nnz=nnz(wbp))


def j_case():
    """Part J's instance: the first 10 × 40, k = 4 instance (Part C's seeds) where basis pursuit fails."""
    for i in range(N_INST_C[(10, 40)]):
        Phi, ws = instance(10, 40, 4, seed_of(10, 4, i))
        if not ok(CS.basis_pursuit(Phi, Phi @ ws, (Phi ** 2).sum(axis=0)), ws):
            return i
    return 0


@cached
def j_job(args):
    """Part J, one scheme: label-noise descent at the same physical OU noise (correlation time TC_J, per-step variance σ_B²) over the same
    physical time T_J — Euler at h, or Heun (the same draw in both stages) — the mean change of log P on the atoms off the truth's support,
    and the end interpolants' success."""
    scheme, h, i = args
    t0 = time.time()
    Phi, ws = instance(10, 40, 4, seed_of(10, 4, i))
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    u0, v0 = CS.uv_of(w0, P0)
    rho, n = float(np.exp(-h / TC_J)), int(round(T_J / h))
    rng = np.random.default_rng(4242)
    if scheme == 'Euler':
        u, v, _ = CS.label_noise_gd(Phi, y, u0, v0, h, SIG_B, rho, n, rng, P_J)
    else:
        u, v = CS.label_noise_heun(Phi, y, u0, v0, h, SIG_B, rho, n, rng, P_J)
    off = ws == 0
    dlogP = np.log(np.maximum(u * v, 1e-300))[:, off] - np.log(P0[off])[None, :]
    ends, conv = CS.end_interpolants(Phi, y, u, v)
    print(f"[{time.time() - t0:.0f}s] J {scheme} h {h:g}", flush=True)
    return dict(scheme=scheme, h=h, n=n, dlogP=float(dlogP.mean()), dlogP_se=float(dlogP.mean(axis=1).std(ddof=1) / np.sqrt(P_J)),
                success=float(np.mean([ok(x, ws) for x in ends])))


def clopper_pearson(x, n, level=0.95):
    """The exact (Clopper–Pearson) interval of a binomial proportion x/n."""
    from scipy.stats import beta
    a = 1 - level
    lo = 0.0 if x == 0 else float(beta.ppf(a / 2, x, n - x + 1))
    hi = 1.0 if x == n else float(beta.ppf(1 - a / 2, x + 1, n - x))
    return lo, hi


def mcnemar_exact(b, c):
    """The exact two-sided McNemar test on the discordant counts b, c: the binomial test of min(b, c) in b + c trials at 1/2."""
    from scipy.stats import binomtest
    return 1.0 if b + c == 0 else float(binomtest(min(b, c), b + c, 0.5).pvalue)


def vertex_limit(Phi, w):
    """At an exact fit on N atoms (Φ_S square, invertible) the slow certificate's τ_c → ∞ limit: max over the atoms i off S of
    |Σ_jα_ij/w_j| / Σ_jα_ij²/|w_j|, α_i = Φ_S⁻¹φ_i (ω_i and c_i both ≈ (N/(4τ_c))× these sums); nan unless |S| = N."""
    S = np.flatnonzero(np.abs(w) > 1e-9 * np.abs(w).max()); I = np.setdiff1d(np.arange(Phi.shape[1]), S)
    if len(S) != Phi.shape[0]:
        return float('nan')
    al = np.linalg.solve(Phi[:, S], Phi[:, I])
    return float(np.max(np.abs((al / w[S][:, None]).sum(axis=0)) / (al ** 2 / np.abs(w[S])[:, None]).sum(axis=0)))


def k_stats(C, G):
    """Part K1 from Parts C and G's stored outcomes: the rates with exact intervals, and paired discordant counts with McNemar's p."""
    cs = {(r['N'], r['d'], r['k'], r['i']): r for r in C}
    gs = {(r['k'], r['i']): r for r in G}
    rates, pairs = [], []
    for N, d in SIZES_C:
        for k in KS_C[(N, d)]:
            ii = [i for i in range(N_INST_C[(N, d)]) if (N, d, k, i) in cs]
            out = {m: [cs[(N, d, k, i)]['out'][m][0] for i in ii] for m in ('bp', 'slow 10', 'slow 100', 'slow 1000')}
            if (N, d) == (10, 40) and all((k, i) in gs for i in ii):
                out['sparse Bayesian learning (EM)'] = [gs[(k, i)]['out']['sparse Bayesian learning (EM)'][0] for i in ii]
            for m, v in out.items():
                x = int(np.sum(v)); rates.append(dict(N=N, d=d, k=k, method=m, x=x, n=len(v), ci=clopper_pearson(x, len(v))))
            for m in ('bp', 'sparse Bayesian learning (EM)', 'slow 1000'):
                if m in out:
                    a, b_ = np.asarray(out['slow 100']), np.asarray(out[m])
                    nb, nc = int(np.sum(a & ~b_)), int(np.sum(~a & b_))
                    pairs.append(dict(N=N, d=d, k=k, other=m, n=len(a), slow_only=nb, other_only=nc, p=mcnemar_exact(nb, nc)))
    return dict(rates=rates, pairs=pairs)


@cached
def k_job(args):
    """Part K2, one 10 × 40 instance where Part C's slow flow (τ_c = 100) fails: the flow run to convergence from the same start, its end
    (support, basic or not), the certificate there for each τ_c and its τ_c → ∞ limit, and R(end) − R(truth) for each τ_c."""
    k, i = args
    Phi, ws = instance(10, 40, k, seed_of(10, k, i))
    y = Phi @ ws
    w0, P0 = start(Phi, y, ALPHA)
    w, P, s, conv = CS.flow_to_convergence(Phi, w0, P0, S_SLOW(100.0), 0, ('ou', 100.0), S_max=S_MAX)
    we = np.where(np.abs(w) > 1e-9 * np.abs(w).max(), w, 0.0)
    S = np.flatnonzero(we)
    basic = bool(len(S) <= Phi.shape[0] and np.linalg.matrix_rank(Phi[:, S]) == len(S))
    De = np.sqrt(w ** 2 + 4 * P ** 2)
    dR = {tc: CS.potential(Phi, De, 0, ('ou', tc)) - CS.potential(Phi, np.abs(ws), 0, ('ou', tc)) for tc in TAUS_KR}
    return dict(k=k, i=i, w=w, P_max=float(P.max()), s=float(s), conv=bool(conv), recovers=ok(w, ws), nnz=int(len(S)), basic=basic,
                cert={tc: certificate(Phi, we, tc)[0] for tc in TAUS_KC} if basic else {}, limit=vertex_limit(Phi, we) if basic else float('nan'),
                dR=dR)


def fmt_end(x):
    return 'n.c.' if not np.isfinite(x) else f'{x:.0e}'


def report(A=None, B=None, C=None, D=None, Ee=None, Ff=None, G=None, H=None, I=None, J=None, K=None):
    L_ = ["# Proposition D‴.43: the colour of the measurement noise selects the kernels (`exp_colour_selection.py`)", "",
          "w = u⊙u − v⊙v on a fixed dictionary Φ (N × d), fitted by gradient descent with label noise; along the exact fits the reduced flow "
          "dP/ds = −ω⊙P, dw/ds = −Π_D(ω⊙w), ω_i = φ_iᵀf(K)φ_i, the gradient flow of R = tr Ψ(K) (`colour_selection`). Success: within "
          f"{TOL:g}·max|w*| of the truth. The start: the noiseless gradient flow's interpolant from u = v = {ALPHA:g} (the hyperbolic "
          "entropy's minimiser on M, P = α²); the end interpolants likewise from each path's end state.", ""]
    if A is not None:
        L_ += ["## A. The potential's identities", "",
               "| N × d | τ | det(I + τK) against its subset expansion (relative) | ∂R/∂D against (4/N)ω (relative) | R's second differences "
               "along 20 random segments (max, min) | white: |second difference| |", "|---|---|---|---|---|---|"]
        for r in A['rows']:
            L_.append(f"| {r['N']} × {r['d']} | {r['tau']:g} | {r['cb']:.1e} | {r['grad']:.1e} | {r['sec_max']:+.2e}, {r['sec_min']:+.2e} | "
                      f"{r['sec_white']:.1e} |")
        L_ += ["", "On a k-sparse exact fit (P = 0, D = |w|): τ_cR − k log(4τ_c/N) against log(det(Φ_SᵀΦ_S)Π_{i∈S}|w_i|) —", "",
               "| N × d | k | the limit | τ_c = 10², 10⁴, 10⁶ |", "|---|---|---|---|"]
        for r in A['lim']:
            L_.append(f"| {r['N']} × {r['d']} | {r['k']} | {r['ref']:+.6f} | {', '.join(f'{v:+.6f}' for _, v in r['vals'])} |")
        L_.append("")
    if B is not None:
        ok_ = lambda e, r: e < TOL * r['wmax']
        B1 = [r for r in B if r['sig'] == SIG_B]
        L_ += [f"## B. Label-noise descent against the reduced flow ({N_B} × {D_B}, a {K_B}-sparse truth, step {H_B:g}, {P_B} paths, the first "
               f"{len(SEEDS_B)} seeds)", "", f"σ = {SIG_B:g}: the paths' mean against the flow (largest deviation over w and the nine times "
               "s = 0, S/8, …, S), the flow's error at S, the fraction of paths whose end interpolant recovers the truth (within "
               f"{TOL:g}·max|w*|), the end interpolants' errors (median over the paths) and supports' sizes (median; n.c. not converged), "
               "the raw ends' error (median), basis pursuit's error:", "",
               "| seed | the noise | steps | s at the end | the paths' mean against the flow | the flow's error at S | paths recovering the "
               "truth | the end interpolants' error | their support | the raw ends' error | basis pursuit's error |",
               "|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in B1:
            noise = 'fresh' if r['tc'] == 0 else f"AR(1), τ_c = {r['tc']:g}"
            nnz = [x for x in r['end_nnz'] if x >= 0]
            L_.append(f"| {r['seed']} | {noise} | {r['n']} | {r['s'][-1]:.0f} | {r['dev'].max():.3f} | {r['flow_err'][-1]:.1e} | "
                      f"{np.mean([ok_(e, r) for e in r['end_err']]):.2f} | {fmt_end(np.median(r['end_err']))} | "
                      f"{('%g' % np.median(nnz)) if nnz else 'n.c.'} | {np.median(r['raw_err']):.2f} | {r['bp_err']:.2f} |")
        L_ += ["", f"Per seed (AR(1), τ_c = 100): the fraction of paths recovering the truth at σ = {SIG_B:g}, " + ', '.join(f"{x:g}" for x in SIGMAS_B)
               + " (the same rescaled time), against the reduced flow's own success and basis pursuit's:", "",
               "| seed | the reduced flow recovers | basis pursuit recovers | " + ' | '.join(f"σ = {x:g}" for x in (SIG_B,) + SIGMAS_B) + " |",
               "|---|---|---|" + "---|" * (1 + len(SIGMAS_B))]
        for seed in SEEDS_B:
            rs = {r['sig']: r for r in B if r['seed'] == seed and r['tc'] == 100.0}
            r0 = rs[SIG_B]
            L_.append(f"| {seed} | {'yes' if ok_(r0['flow_err'][-1], r0) else 'no'} | {'yes' if ok_(r0['bp_err'], r0) else 'no'} | "
                      + ' | '.join(f"{np.mean([ok_(e, rs[x]) for e in rs[x]['end_err']]):.2f}" for x in (SIG_B,) + SIGMAS_B) + " |")
        L_.append("")
    if C is not None:
        meth = list(C[0]['out'])
        L_ += ["## C. Recovery against sparsity (the reduced flows, continuum filters)", "",
               "| N × d | k | instances | " + ' | '.join(meth) + " | the slow flow's (τ_c = 100) mean support size | basis pursuit's |",
               "|---|---|---|" + "---|" * len(meth) + "---|---|"]
        keys = sorted({(r['N'], r['d'], r['k']) for r in C})
        for N, d, k in keys:
            rs = [r for r in C if (r['N'], r['d'], r['k']) == (N, d, k)]
            L_.append(f"| {N} × {d} | {k} | {len(rs)} | " + ' | '.join(f"{np.mean([r['out'][m][0] for r in rs]):.2f}" for m in meth) +
                      f" | {np.mean([r['out']['slow 100'][1] for r in rs]):.1f} | {np.mean([r['out']['bp'][1] for r in rs]):.1f} |")
        L_ += ["", f"The white flow run to convergence (max P < 10⁻⁸·max|w|, w still to 10⁻¹⁰, s ≤ {S_MAX:.0e}): converged in "
               f"{sum(r['white_conv'] for r in C)} of {len(C)} instances, its largest distance to basis pursuit's solution "
               f"{max(r['white_bp'] for r in C):.1e}. Where a slow flow fails, the fraction of its ends with R below the truth's (the truth not "
               "the global minimiser) — τ_c = 10, 100:", "", "| N × d | k | failures (10, 100) | R(end) < R(truth) |", "|---|---|---|---|"]
        for N, d, k in keys:
            rs = [r for r in C if (r['N'], r['d'], r['k']) == (N, d, k)]
            cells = []
            for tc in (10.0, 100.0):
                f_ = [r for r in rs if not r['out'][f'slow {tc:g}'][0]]
                cells.append((len(f_), (np.mean([r['R'][tc][0] < r['R'][tc][1] for r in f_]) if f_ else float('nan'))))
            L_.append(f"| {N} × {d} | {k} | {cells[0][0]}, {cells[1][0]} | " + ', '.join('—' if not np.isfinite(x) else f"{x:.2f}" for _, x in cells) + " |")
        L_.append("")
    if D is not None:
        Dr, Dn = D['rates'], D['noisy']
        meth = list(Dr[0]['out'])
        L_ += [f"## D. Gaussian kernel dictionaries on [0, 1] ({D_D} centres, {N_D} equispaced samples), random k-atom truths ({N_INST_D} per row)", "",
               "| width | the largest off-diagonal |φ_iᵀφ_j| | k | " + ' | '.join(meth) + " | mean support size: basis pursuit, slow 100 |",
               "|---|---|---|" + "---|" * len(meth) + "---|"]
        for width in WIDTHS_D:
            Phi = dict_instance(width, 3, 0)[0]
            coh = float(np.abs(Phi.T @ Phi - np.eye(D_D)).max())
            for k in KS_D:
                rs = [r for r in Dr if r['width'] == width and r['k'] == k]
                L_.append(f"| {width:g} | {coh:.2f} | {k} | " + ' | '.join(f"{np.mean([r['out'][m][0] for r in rs]):.2f}" for m in meth) +
                          f" | {np.mean([r['out']['bp'][1] for r in rs]):.1f}, {np.mean([r['out']['slow 100'][1] for r in rs]):.1f} |")
        L_ += ["", f"Label-noise descent on the first instances where basis pursuit fails (width {WIDTHS_D[1]:g}, k = 4; {P_B} paths, σ = {SIG_B:g}, "
               "the step min(0.5, 1/(2λ_max)) at the start):", "",
               "| instance | step | the noise | steps | the end interpolants' errors (per path) | their supports' sizes | the ends' atoms "
               "alive (D > 10⁻⁸ max D) | the raw ends' distance to the truth, to basis pursuit's solution (max over paths) | the reduced "
               "flow at the same s: its error, support, distance to basis pursuit's solution | the reduced flow run to convergence: its error, "
               "support, the s reached |", "|---|---|---|---|---|---|---|---|---|---|"]
        for r in Dn:
            for tc in (0.0, 100.0):
                q = r[tc]
                L_.append(f"| {r['i']} | {r['h']:.3g} | {'fresh' if tc == 0 else 'AR(1), τ_c = 100'} | {q['n']} | {', '.join(fmt_end(x) for x in q['err'])} | "
                          f"{', '.join(str(int(x)) if x >= 0 else 'n.c.' for x in q['nnz'])} | {', '.join(str(int(x)) for x in q['alive'])} | "
                          f"{q['raw_err'].max():.2f}, {q['raw_bp'].max():.2f} | {q['flow_err']:.1e}, {q['flow_nnz']}, {q['flow_bp']:.1e} | "
                          f"{q['conv_err']:.1e}, {q['conv_nnz']}, {q['conv_s']:.0f}{'' if q['conv'] else ' (not converged)'} |")
        L_ += ["", "n.c.: the end's flow interpolant does not converge — the exact fit needs atoms whose masses have underflowed (fewer than N "
               "alive, a residual outside their span), its maximiser is at infinity; the raw end (residual O(σ)) is reported instead.", ""]
    if Ee is not None:
        names = list(Ee[0]['out'])
        L_ += [f"## E. Other colours at the step {H_B:g} (the reduced flow with each scheme's filter; 10 × 40, {N_INST_E} instances)", "",
               "| k | " + ' | '.join(names) + " |", "|---|" + "---|" * len(names)]
        for k in KS_E:
            rs = [r for r in Ee if r['k'] == k]
            L_.append(f"| {k} | " + ' | '.join(f"{np.mean([r['out'][n_] for r in rs]):.2f}" for n_ in names) + " |")
        L_.append("")
    if Ff is not None:
        L_ += [f"## F. The start: the interpolant of u = v = α (10 × 40), the flows run to convergence (s ≤ {S_MAX:.0e})", "",
               "| k | α | the interpolant itself | the white flow | the slow flow (τ_c = 100) | converged (white, slow) | median s reached "
               "(white, slow) | basis pursuit |", "|---|---|---|---|---|---|---|---|"]
        for k in KS_E:
            rs = [r for r in Ff if r['k'] == k]
            for a in ALPHAS_F:
                cv = [r['conv'][a] for r in rs]
                L_.append(f"| {k} | {a:g} | " + ' | '.join(f"{np.mean([r['out'][a][j] for r in rs]):.2f}" for j in range(3)) +
                          f" | {np.mean([c[0] for c in cv]):.2f}, {np.mean([c[2] for c in cv]):.2f} | {np.median([c[1] for c in cv]):.3g}, "
                          f"{np.median([c[3] for c in cv]):.3g} | {np.mean([r['out']['bp'] for r in rs]):.2f} |")
        L_.append("")
    if G is not None:
        meth = list(G[0]['out'])
        L_ += [f"## G. References on the same instances (10 × 40, {N_INST_G} per k): reweighted ℓ¹ (Candès–Wakin–Boyd; {IT_G} iterations, ε "
               "relative to max|w*|, from basis pursuit's solution or from the flows' start), the log-det's majorisation–minimisation "
               f"(τ_c = {TC_G:g}), sparse Bayesian learning (Tipping's EM, λ = 10⁻¹⁰) — the fraction recovering the truth, and the slow "
               "flow's (τ_c = 100, Part C) beside them", "", "| method | " + ' | '.join(f"k = {k}" for k in KS_G) + " | mean support (k = 4) |",
               "|---|" + "---|" * (len(KS_G) + 1)]
        for m in meth:
            L_.append(f"| {m} | " + ' | '.join(f"{np.mean([r['out'][m][0] for r in G if r['k'] == k]):.2f}" for k in KS_G) +
                      f" | {np.mean([r['out'][m][1] for r in G if r['k'] == 4]):.1f} |")
        if C is not None:
            cs = {(r['k'], r['i']): r for r in C if (r['N'], r['d']) == (10, 40)}
            L_.append("| the slow flow (τ_c = 100) | " + ' | '.join(
                f"{np.mean([cs[(k, r['i'])]['out']['slow 100'][0] for r in G if r['k'] == k and (k, r['i']) in cs]):.2f}" for k in KS_G)
                + f" | {np.mean([cs[(4, r['i'])]['out']['slow 100'][1] for r in G if r['k'] == 4 and (4, r['i']) in cs]):.1f} |")
        L_.append("")
    if H is not None:
        L_ += ["## H. Which exact fits the selection holds: the dual certificate max_{i∉S}|c_i|/ω_i (below 1: the fit is stable)", "",
               "| k | instances | at the truth: the fraction certified, white; τ_c = 10, 100, 1000 | the median ratio at the truth, white; τ_c = 10, 100, 1000 | "
               "where basis pursuit fails: its solution's support size, and the fraction of its solutions certified for τ_c = 10, 100, 1000 | "
               "basis pursuit's success |", "|---|---|---|---|---|---|"]
        Hk = H['kernel'] if isinstance(H, dict) else None
        H = H['gauss'] if isinstance(H, dict) else H
        for k in KS_G:
            rs = [r for r in H if r['k'] == k]
            fails = [r for r in rs if not r['bp']]
            frac_truth = '; '.join(f"{np.mean([r['truth'][tc] < 1 for r in rs]):.2f}" for tc in TAUS_H)
            med_truth = '; '.join(f"{np.median([r['truth'][tc] for r in rs]):.3g}" for tc in TAUS_H)
            at_bp = (f"{np.mean([r['bp_nnz'] for r in fails]):.1f}; " + '; '.join(f"{np.mean([r['at_bp'][tc] < 1 for r in fails]):.2f}" for tc in TAUS_H[1:])
                     if fails else '—')
            L_.append(f"| {k} | {len(rs)} | {frac_truth} | {med_truth} | {at_bp} | {np.mean([r['bp'] for r in rs]):.2f} |")
        L_ += ["", "The certificate is sufficient, not necessary (for white draws it is Fuchs's least-squares certificate): basis pursuit's own "
               "success is in the last column.", ""]
        if Hk is not None:
            L_ += [f"On the kernel dictionaries ({D_D} centres, {N_D} samples; {N_INST_H} instances per row): the fraction of truths certified, "
                   "white and τ_c = " + ', '.join(f"{t:.0e}" for t in TAUS_HK) + "; the fraction to which the slow flow returns from P = 10⁻³ "
                   "(τ_c = 100, 1000); basis pursuit's success:", "",
                   "| width | k | certified: white; " + '; '.join(f"{t:.0e}" for t in TAUS_HK) + " | returns (100, 1000) | basis pursuit |",
                   "|---|---|---|---|---|"]
            for width in WIDTHS_H:
                for k in KS_H:
                    rs = [r for r in Hk if r['width'] == width and r['k'] == k]
                    L_.append(f"| {width:g} | {k} | " + '; '.join(f"{np.mean([r['truth'][t] < 1 for r in rs]):.2f}" for t in (0.0,) + TAUS_HK)
                              + f" | {np.mean([r['hold'][100.0] for r in rs]):.2f}, {np.mean([r['hold'][1000.0] for r in rs]):.2f} | "
                              f"{np.mean([r['bp'] for r in rs]):.2f} |")
            L_.append("")
    if I is not None:
        L_ += ["## I. A fit floor: y = Φw* + ε·e (|e| = 1; 10 × 40, the instances of Part C) — the fraction whose k largest atoms are the truth's "
               "support, and the support sizes (|w| > 10⁻³·max|w*|)", "",
               "| k | ε | the slow flow (τ_c = 100) | basis pursuit | support sizes (slow, basis pursuit) |", "|---|---|---|---|---|"]
        for k in KS_I:
            for eps in EPS_I:
                rs = [r for r in I if r['k'] == k and r['eps'] == eps]
                L_.append(f"| {k} | {eps:g} | {np.mean([r['slow'] for r in rs]):.2f} | {np.mean([r['bp'] for r in rs]):.2f} | "
                          f"{np.mean([r['slow_nnz'] for r in rs]):.1f}, {np.mean([r['bp_nnz'] for r in rs]):.1f} |")
        L_.append("")
    if J is not None:
        L_ += [f"## J. The selection is the explicit scheme's: label-noise descent at the same physical noise (AR(1) of correlation time "
               f"{TC_J:g} time units, per-step σ = {SIG_B:g}) over the same physical time {T_J:g} ({P_J} paths; Part C's first 10 × 40, k = 4 "
               "instance where basis pursuit fails)", "", "| scheme | step | steps | mean Δlog P off the truth's support (± s.e.) | "
               "the end interpolants recovering the truth |", "|---|---|---|---|---|"]
        for r in J:
            L_.append(f"| {r['scheme']} | {r['h']:g} | {r['n']} | {r['dlogP']:+.4f} ± {r['dlogP_se']:.1e} | {r['success']:.2f} |")
        L_.append("")
    if K is not None:
        st, cp = K['stats'], K['competitors']
        L_ += ["## K. Statistics and the competing fits (referee round 142)", "",
               "K1. The recovery rates of Parts C and G with exact (Clopper–Pearson) 95% intervals — successes/instances, rate [interval]:", "",
               "| N × d | k | " + ' | '.join(m for m in ('bp', 'slow 10', 'slow 100', 'slow 1000', 'sparse Bayesian learning (EM)')) + " |",
               "|---|---|---|---|---|---|---|"]
        keys = sorted({(r['N'], r['d'], r['k']) for r in st['rates']})
        for N, d, k in keys:
            cell = {r['method']: r for r in st['rates'] if (r['N'], r['d'], r['k']) == (N, d, k)}
            L_.append(f"| {N} × {d} | {k} | " + ' | '.join(
                (f"{cell[m]['x']}/{cell[m]['n']}, {cell[m]['x'] / cell[m]['n']:.2f} [{cell[m]['ci'][0]:.2f}, {cell[m]['ci'][1]:.2f}]" if m in cell else '—')
                for m in ('bp', 'slow 10', 'slow 100', 'slow 1000', 'sparse Bayesian learning (EM)')) + " |")
        L_ += ["", "Paired on the same instances: the slow flow (τ_c = 100) against another method — the instances only the slow flow recovers, "
               "those only the other recovers, and the exact two-sided McNemar p-value (the binomial test on the discordant pairs):", "",
               "| N × d | k | instances | against | slow only | the other only | p |", "|---|---|---|---|---|---|---|"]
        for r in st['pairs']:
            L_.append(f"| {r['N']} × {r['d']} | {r['k']} | {r['n']} | {r['other']} | {r['slow_only']} | {r['other_only']} | {r['p']:.2g} |")
        L_ += ["", f"K2. Where the slow flow (τ_c = 100) fails at its fixed horizon (Part C, 10 × 40), the same flow run to convergence (s ≤ {S_MAX:.0e}): "
               "the failures, how many converged and still fail (agreeing with Part C), the end's support (median, range), the fraction of ends "
               "that are basic exact fits on N = 10 atoms, the fraction certified at the end (max_{i∉S}|c_i|/ω_i < 1) for τ_c = 100, 10³, 10⁴, 10⁵ "
               "and in the τ_c → ∞ limit, that limit's median:", "",
               "| k | failures | converged, failing | the end's support | basic on N atoms | certified at the end: τ_c = 100; 10³; 10⁴; 10⁵; ∞ | "
               "the limit's median ratio |", "|---|---|---|---|---|---|---|"]
        for k in KS_K:
            rs = [r for r in cp if r['k'] == k]
            if not rs:
                continue
            nb = [r for r in rs if r['basic'] and r['nnz'] == 10]
            fr = lambda tc: f"{np.mean([r['cert'][tc] < 1 for r in nb]):.2f}" if nb else '—'
            lim = [r['limit'] for r in nb]
            L_.append(f"| {k} | {len(rs)} | {sum(r['conv'] and not r['recovers'] for r in rs)} | {int(np.median([r['nnz'] for r in rs]))} "
                      f"({min(r['nnz'] for r in rs)}–{max(r['nnz'] for r in rs)}) | {len(nb) / len(rs):.2f} | "
                      + '; '.join(fr(tc) for tc in TAUS_KC) + f"; {np.mean([x < 1 for x in lim]):.2f} | {np.median(lim):.3g} |" if nb else
                      f"| {k} | {len(rs)} | {sum(r['conv'] and not r['recovers'] for r in rs)} | — | 0.00 | — | — |")
        L_ += ["", "R at the converged end against R at the truth, the same two points evaluated with the potential of each τ_c (the slow flow "
               "ran at τ_c = 100): the fraction of the failures whose end has the lower R, and the median of τ_c(R(end) − R(truth)):", "",
               "| k | " + ' | '.join(f"τ_c = {tc:g}" for tc in TAUS_KR) + " |", "|---|" + "---|" * len(TAUS_KR)]
        for k in KS_K:
            rs = [r for r in cp if r['k'] == k]
            if rs:
                L_.append(f"| {k} | " + ' | '.join(f"{np.mean([r['dR'][tc] < 0 for r in rs]):.2f}; {np.median([tc * r['dR'][tc] for r in rs]):+.3g}"
                                                 for tc in TAUS_KR) + " |")
        L_.append("")
    L_.append("Numbers are produced by `experiments/exp_colour_selection.py`; nothing is edited by hand.")
    return '\n'.join(L_) + '\n'


def main():
    raw = os.path.join(OUT, 'colour_selection_raw.npz')
    md = os.path.join(OUT, 'colour_selection.md')
    old = E.load(raw)
    get = lambda k: old[k] if k in old else None
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[get(k) for k in 'ABCDEFGHIJK']))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set('ABCDEFGHIJK')
    from multiprocessing import Pool
    res = dict(old)
    t0 = time.time()
    log = lambda k: print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
    with Pool(N_WORKERS) as pool:
        h = {}
        if 'B' in parts:
            h['B'] = pool.map_async(b_job, [(s_, tc, SIG_B) for s_ in SEEDS_B for tc in TAUS_B]
                                    + [(s_, 100.0, sg) for s_ in SEEDS_B for sg in SIGMAS_B], chunksize=1)
        if 'D' in parts:
            h['D'] = pool.map_async(d_job, [(w_, k, i) for w_ in WIDTHS_D for k in KS_D for i in range(N_INST_D)], chunksize=1)
            h['Dn'] = pool.map_async(d_noisy_job, d_noisy_cases(), chunksize=1)
        if 'C' in parts:
            h['C'] = pool.map_async(c_job, [(N, d, k, i) for N, d in SIZES_C for k in KS_C[(N, d)] for i in range(N_INST_C[(N, d)])], chunksize=1)
        if 'E' in parts:
            h['E'] = pool.map_async(e_job, [(k, i) for k in KS_E for i in range(N_INST_E)], chunksize=1)
        if 'F' in parts:
            h['F'] = pool.map_async(f_job, [(k, i) for k in KS_E for i in range(N_INST_E)], chunksize=1)
        if 'G' in parts:
            h['G'] = pool.map_async(g_job, [(k, i) for k in KS_G for i in range(N_INST_G)], chunksize=1)
        if 'H' in parts:
            h['H'] = pool.map_async(h_job, [(k, i) for k in KS_G for i in range(N_INST_G)], chunksize=1)
            h['Hk'] = pool.map_async(hk_job, [(w_, k, i) for w_ in WIDTHS_H for k in KS_H for i in range(N_INST_H)], chunksize=1)
        if 'I' in parts:
            h['I'] = pool.map_async(i_job, [(e_, k, i) for e_ in EPS_I for k in KS_I for i in range(N_INST_I)], chunksize=1)
        if 'K' in parts:
            Cc = h['C'].get() if 'C' in h else res['C']
            fails = [(r['k'], r['i']) for r in Cc if (r['N'], r['d']) == (10, 40) and r['k'] in KS_K and not r['out']['slow 100'][0]]
            h['K'] = pool.map_async(k_job, sorted(fails), chunksize=1)
        if 'J' in parts:
            ij = j_case()
            h['J'] = pool.map_async(j_job, [('Euler', H_B, ij), ('Heun', H_B, ij), ('Euler', H_B / 4, ij)], chunksize=1)
        if 'A' in parts:
            res['A'] = part_A(); E.save(raw, **res); log('A')
        for k in 'BDCEFGHIJK':
            if k in h:
                res[k] = (dict(rates=h['D'].get(), noisy=h['Dn'].get()) if k == 'D' else
                          dict(gauss=h['H'].get(), kernel=h['Hk'].get()) if k == 'H' else
                          dict(stats=k_stats(res['C'], res['G']), competitors=h['K'].get()) if k == 'K' else h[k].get())
                E.save(raw, **res); log(k)
    open(md, 'w').write(report(*[res.get(k) for k in 'ABCDEFGHIJK']))


if __name__ == '__main__':
    main()
