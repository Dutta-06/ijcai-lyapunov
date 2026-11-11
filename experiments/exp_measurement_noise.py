"""
Proposition D‴.30: measurement noise at perfect tracking — a coherent, high-passed noise; the floors of the fixed and the adaptive prior
=======================================================================================================================================

The three-bump instance (σ₀ = 0.15σ_crit, x_m = 2 sin t), the real closed loop with the nearest-atom leak toward a fixed (η = 0) or an
adaptive prior (η > 0), at the lift of the exact representation of the matching weights q = 1/3 with the prior at its atoms (perfect
tracking), and the measurement y = x + ν used by the controller and the adaptive law (src/core/measurement_noise.py).

A. The linear theory (no dynamics): H = ⟨ĴᵀĴ⟩ and Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩ in the eigenbasis of H (a = ∂_xf̂ on the orbit); per mode the
   function-level floor's share (N₀/2k_e)Q_jj(h_j + η)/(h_j + η + σ₀) and the parameter variance at η = 0, σ₀, 10σ₀, ∞, with the
   bias-removal rate λ_j^−(η); the second-order (errors-in-variables) force per unit of the noise's variance and the stationary mean
   deviations it produces, (H + σ₀I)⁻¹ against H⁻¹ times k_e times the force. The full second-order mean (`second_order_mean`): the force's
   noise part against the quasi-static one, the parameters' fluctuation's part, the biases of the mean parameters and of the mean function,
   and the exact periodic mean of the linearised loop; the bias's growth in k_e at B = 100, 10³, 10⁴ (local exponents); the separations'
   second-order drift S_k per class and the noise level at which a class's separations stop decaying.
B. The exact linearised loop driven by white measurement noise: its periodic stationary covariance against the averaged prediction, at
   k_e = 25, …, 800 and η = 0, σ₀, ∞ (the function-level floor ⟨Ĵ P_yy Ĵᵀ⟩ over a period, the weakest mode's variance against its slow
   and fast parts); a control with the slope a removed from the noise input (−k_eν), where the deviation is even in 1/k_e; and the exact
   stationary second-order mean of the linearised loop with the OU noise (MR) against the averaged bias, k_e = 25–200.
C. The real nonlinear loop (the weighted K₀-atom loop; RK4, dt = 10⁻³, k_e = 25) with an Ornstein–Uhlenbeck measurement noise (rate
   B = 100, exact at the steps), N_MC independent runs per case, T_C time units, at two noise levels: E‖f̂ − f‖²_ρ(t) against the linear
   prediction — the transient covariance of the exactly linearised loop with the OU noise, plus the second-order mean deviation — for the
   fixed and the adaptive prior (η = ∞). And one run of the lift m = 15 driven by the same noise path as a K₀-atom run: the same f̂ to
   rounding (the noise is coherent: a replicated start stays replicated; f̂ is stored in float32).
D. Invariance is not attraction: the float64 difference between the lift and the K₀-atom run over 25 time units; and the lift m = 15 from
   a separation of class 3 of size 10⁻⁶ along S₃'s top eigenvector, N_D antithetic pairs of noise paths (ν and −ν) at the stds STDS_D: the
   separation's growth rate over t ∈ [50, T_D] (the pair means cancel its first-order part) and the noise-induced rate
   λ_eff := (rate + σ₀/k_e)/Eν², against the lift's drift λ_max(S₃) and against the drift along the second-order transient from the lift,
   whose mean error carries the mean tracking error of the developing errors-in-variables deviation (round 22, ID-1): S₃(t) by
   `second_order_transient` and `separation_drift_series`, its window averages, and the predicted slope of log‖u‖ for u̇ = (Eν²S₃(t) − σ₀/k_e)u.

Output: results/chaos/measurement_noise.md and measurement_noise_raw.npz (MN_WORKERS, MN_RUNS, MN_TC; `--from-raw` regenerates the md,
`--only-A`, `--only-B`, `--only-D` recompute that part and merge it).
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
from src.core import langevin as L
from src.core import adaptive_prior as AP
from src.core import adaptive_floquet as AF
from src.core import measurement_noise as MN

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q3 = np.ones(3) / 3
KE = 25.0
ETAS_A = (0.0, SIGMA0, 10 * SIGMA0, np.inf)
KES_B = (25.0, 50.0, 100.0, 200.0, 400.0, 800.0)
KES_B2 = (25.0, 50.0, 100.0, 200.0)
KES_GROW = (12.5, 25.0, 50.0, 100.0, 200.0, 400.0)
BS_GROW = (100.0, 1000.0, 10000.0)
STDS_D = (0.01, 0.03, 0.05)
N_D = 8                                                                    # antithetic pairs per std
T_D = 300.0
ETAS_B = (0.0, SIGMA0, np.inf)
B_OU, DT_C = 100.0, 1e-3
STDS_C = (0.003, 0.01)
N_MC = int(os.environ.get('MN_RUNS', 12))
T_C = float(os.environ.get('MN_TC', 300.0))
REC_C = 25.0
N_WORKERS = int(os.environ.get('MN_WORKERS', 3))
N_GRID = 256


def f_vec(X):
    X = np.asarray(X, float)
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def fprime(x):
    return sum(-A * (x - C) / S ** 2 * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def fsecond(x):
    return sum(A * ((x - C) ** 2 / S ** 4 - 1 / S ** 2) * np.exp(-(x - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


xr = lambda t: np.array([2.0 * np.sin(t)])
xrd = lambda t: np.array([2.0 * np.cos(t)])


def setup():
    net = L.SaturatedLiftedNetwork(1, **NET_KW)
    A = np.array([b[0] for b in BUMPS]); C = np.array([b[1] for b in BUMPS]); S = np.array([b[2] for b in BUMPS])
    atoms = net.unconstrained((A / Q3)[:, None], C[:, None], S)
    return net, atoms


def save_raw(path, **arrays):
    tmp = path + '.tmp.npz'
    np.savez(tmp, **arrays)
    os.replace(tmp, path)


# ---------------------------------------------------------------------------- A
def part_A():
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, np.inf, atoms, Q3)
    H, Q = MN.noise_matrices(cv, fprime)
    h, V = np.linalg.eigh(H)
    qd = np.einsum('ij,ik,kj->j', V, Q, V)
    t = np.arange(1024) * PER / 1024
    J = np.array([cv.J(s)[0] for s in t]); k = np.fft.fftfreq(1024, d=1.0 / 1024)
    Jd = np.real(np.fft.ifft((1j * k)[:, None] * np.fft.fft(J, axis=0), axis=0))
    a = np.array([fprime(2 * np.sin(s)) for s in t])
    parts = dict(deriv=float(np.trace(Jd.T @ Jd) / 1024), cross=float(2 * np.trace(Jd.T @ (a[:, None] * J)) / 1024),
                 slope=float(np.trace((a[:, None] * J).T @ (a[:, None] * J)) / 1024))
    rows = {}
    for eta in ETAS_A:
        fl, sh, pv = MN.floor_closed_form(H, Q, SIGMA0, eta, 1.0, KE)
        lm = (h + SIGMA0) * 0 if eta == 0 else AP.pair_rates(h, SIGMA0, eta)[0]
        rows[eta] = dict(floor=fl, shares=sh, pvar=pv, rate=lm)
    F1 = MN.eiv_force(cv, fprime, fsecond, KE, B_OU, 1.0)          # per unit of the noise's variance (real time, whitened)
    Fs = KE * F1                                                    # slow time
    fv = V.T @ Fs
    bias = {0.0: float(np.sum(fv ** 2 * h / (h + SIGMA0) ** 2)), np.inf: float(np.sum(fv ** 2 / h))}
    trH2, J4 = MN.fast_parts(cv)
    # the full second-order mean: the noise part against the quasi-static force, the parameters' fluctuation, the exact periodic mean
    so = {}
    for eta in (0.0, np.inf):
        cl = MN.NoisyRealLoop(net, f_vec, xr, xrd, PER, KE, G, SIGMA0, np.arange(3), 3, eta, weights=Q3)
        out = MN.second_order_mean(cl, cv, fprime, B_OU, 1.0, periodic=True)
        bs = MN.mean_biases(H, out, SIGMA0 if eta == 0 else None)
        my = out['m'][:, 1:1 + cv.d]
        so[eta] = dict(bs, rel=float(np.linalg.norm(out['F_noise'] - Fs) / np.linalg.norm(Fs)), exact=float(np.mean(np.einsum('ti,ij,tj->t', my, H, my))),
                       per=max(out['per_err'], out['m_per_err']))
    grow = {B: [MN.equivalent_perturbation(H, ke * MN.eiv_force(cv, fprime, fsecond, ke, B, 1.0)) for ke in KES_GROW] for B in BS_GROW}
    S = MN.separation_drift(cv, fprime, fsecond, KE, B_OU, 1.0)
    sep = [float(np.linalg.eigvalsh(Sk).max()) for Sk in S]
    print("[A]", flush=True)
    return dict(H=H, Q=Q, h=h, V=V, qd=qd, parts=parts, rows=rows, F1=F1, fv=fv, bias=bias, trH2=trH2, J4=J4, so=so, grow=grow, sep=sep)


# ---------------------------------------------------------------------------- B
def job_B(args):
    ke, eta = args[:2]
    slope = len(args) < 3 or args[2] == 'slope'                         # 'noslope': the noise input −k_eν (the control)
    t0 = time.time()
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, eta if eta > 0 else np.inf, atoms, Q3)
    if eta == 0:
        gen = MN.fixed_generator(cv, ke); D = 1 + cv.d
    else:
        gen = lambda t: cv.generator(t, ke); D = cv.dim
    b0 = MN.collapsed_noise_input(cv, fprime, ke, compensation='measured' if slope else 'reference')
    bv = (lambda t: b0(t)[:D])
    P0 = MN.periodic_covariance(gen, bv, PER, D, 1.0)
    # the function-level floor over a period: ⟨Ĵ(t)P_yy(t)Ĵ(t)ᵀ⟩ along the periodic solution
    ts = np.linspace(0, PER, 65)
    tt, Ps = MN.transient_covariance(gen, bv, D, 1.0, PER, ts, P0=P0, rtol=1e-10)
    H, Q = MN.noise_matrices(cv, fprime if slope else (lambda x: 0.0))
    trH2, J4 = MN.fast_parts(cv)
    point = [cv.J(t)[0] @ P[1:1 + cv.d, 1:1 + cv.d] @ cv.J(t)[0] for t, P in zip(tt, Ps)]
    orbit = [np.trace(H @ P[1:1 + cv.d, 1:1 + cv.d]) for P in Ps]
    exact_point = float(np.trapezoid(point, tt) / PER); exact_orbit = float(np.trapezoid(orbit, tt) / PER)
    slow = MN.floor_closed_form(H, Q, SIGMA0, eta, 1.0, ke)[0]
    pred_orbit, pred_point = MN.floors(H, Q, trH2, J4, SIGMA0, eta, 1.0, ke)
    # per mode: the period-averaged variance against the slow part (N₀/2k_e)Q_jj(h_j + η)/(h_j(h_j + η + σ₀)) plus the fast part (N₀/2k_e)h_j
    h, V = np.linalg.eigh(H)
    Pm = np.trapezoid(np.array([P[1:1 + cv.d, 1:1 + cv.d] for P in Ps]), tt, axis=0) / PER
    mode_exact = np.einsum('ij,ik,kj->j', V, Pm, V)
    mode_slow = MN.floor_closed_form(H, Q, SIGMA0, eta, 1.0, ke)[2]
    print(f"[{time.time() - t0:.0f}s] B ke {ke:g} eta {eta:g} {'slope' if slope else 'noslope'}", flush=True)
    return dict(ke=ke, eta=eta, slope=slope, exact_orbit=exact_orbit, exact_point=exact_point, slow=slow, pred_orbit=pred_orbit,
                pred_point=pred_point, mode_exact=mode_exact, mode_slow=mode_slow, mode_fast=h / (2 * ke))


def job_B2(args):
    """The exact stationary second-order mean of the linearised loop with the OU noise (MR, B = B_OU) — its mean parameters' function-level
    bias over a period — against the averaged bias of the quasi-static force."""
    ke, eta = args
    t0 = time.time()
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, np.inf, atoms, Q3)
    H, Q = MN.noise_matrices(cv, fprime)
    cl = MN.NoisyRealLoop(net, f_vec, xr, xrd, PER, ke, G, SIGMA0, np.arange(3), 3, eta, weights=Q3)
    out = MN.second_order_mean(cl, cv, fprime, B_OU, 1.0, periodic=True)
    my = out['m'][:, 1:1 + cv.d]
    Fq = ke * MN.eiv_force(cv, fprime, fsecond, ke, B_OU, 1.0)
    avg = MN.response_size(H, Fq, SIGMA0) if eta == 0 else MN.equivalent_perturbation(H, Fq)
    print(f"[{time.time() - t0:.0f}s] B2 ke {ke:g} eta {eta:g}", flush=True)
    return dict(ke=ke, eta=eta, exact=float(np.mean(np.einsum('ti,ij,tj->t', my, H, my))), averaged=avg,
                per=max(out['per_err'], out['m_per_err']))


# ---------------------------------------------------------------------------- C
def job_C(args):
    kind, eta, std, seed = args
    t0 = time.time()
    net, atoms = setup()
    grid = 2 * np.sin(np.arange(N_GRID) * PER / N_GRID); Fg = f_vec(grid)
    if kind == 'lift':
        labels, weights, m = np.repeat(np.arange(3), 5), None, 15
    else:
        labels, weights, m = np.arange(3), Q3, 3
    cl = MN.NoisyRealLoop(net, f_vec, xr, xrd, PER, KE, G, SIGMA0, labels, 3, eta, weights=weights)
    z0 = cl.lift(atoms)

    def rec(t, z):
        Th = z[1:1 + m * 3].reshape(m, 3)
        fh = cl.f_hat(grid[:, None], Th) if False else np.array([cl.f_hat(np.array([x]), Th)[0] for x in grid])
        return np.concatenate([[np.mean((fh - Fg) ** 2)], fh])
    out = MN.simulate_noisy(cl, z0, m, T_C, DT_C, B_OU, std ** 2, np.random.default_rng(1000 + seed), int(round(REC_C / DT_C)), rec)
    out = np.array(out)
    print(f"[{time.time() - t0:.0f}s] C {kind} eta {eta:g} std {std:g} seed {seed}", flush=True)
    return dict(kind=kind, eta=eta, std=std, seed=seed, err=out[:, 0], fh=out[:, 1:].astype(np.float32) if kind == 'lift' or seed == 0 else None)


def predictions_C(std, eta):
    """The linear prediction of E‖f̂ − f‖²_ρ at the records: the transient covariance of the exactly linearised loop with the OU noise (from
    the stationary noise and P = 0) plus the second-order mean deviation."""
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, np.inf, atoms, Q3)
    var = std ** 2; N0 = 2 * var / B_OU
    if eta == 0:
        gen = MN.fixed_generator(cv, KE); D = 1 + cv.d
    else:
        gen = lambda t: cv.generator(t, KE); D = cv.dim
    b0 = MN.collapsed_noise_input(cv, fprime, KE)
    Ma, cfun = MN.ou_augmented(gen, lambda t: b0(t)[:D], D, B_OU, N0)
    P0 = np.zeros((D + 1, D + 1)); P0[-1, -1] = var
    tev = np.arange(0.0, T_C + 1e-9, REC_C)
    tt, Ps = MN.transient_covariance(Ma, cfun, D + 1, N0, T_C, tev, P0=P0)
    H, Q = MN.noise_matrices(cv, fprime)
    cov = np.array([np.trace(H @ P[1:1 + cv.d, 1:1 + cv.d]) for P in Ps])
    F = KE * MN.eiv_force(cv, fprime, fsecond, KE, B_OU, var)
    mean = np.array([(lambda y: y @ H @ y)(MN.mean_response(H, F, SIGMA0, eta, t / KE)) for t in tev])
    return dict(t=tev, cov=cov, mean=mean)


# ---------------------------------------------------------------------------- D
def job_D0(args):
    """The lift m = 15 and the K₀-atom loop driven by one noise path (std 0.01, η = ∞) over 25 time units, in float64: the largest difference
    of f̂ on the grid and the lift's largest within-class spread."""
    t0 = time.time()
    net, atoms = setup()
    grid = 2 * np.sin(np.arange(N_GRID) * PER / N_GRID)
    res = {}
    for kind in ('lift', 'K0'):
        labels, weights, m = (np.repeat(np.arange(3), 5), None, 15) if kind == 'lift' else (np.arange(3), Q3, 3)
        cl = MN.NoisyRealLoop(net, f_vec, xr, xrd, PER, KE, G, SIGMA0, labels, 3, np.inf, weights=weights)
        out = MN.simulate_noisy(cl, cl.lift(atoms), m, 25.0, DT_C, B_OU, 0.01 ** 2, np.random.default_rng(1000), int(round(25.0 / DT_C)),
                                lambda t, z: z.copy())
        Th = out[-1][1:1 + m * 3].reshape(m, 3)
        res[kind] = (np.array([cl.f_hat(np.array([x]), Th)[0] for x in grid]), Th)
    spread = max(float(np.abs(res['lift'][1][np.repeat(np.arange(3), 5) == k] - res['lift'][1][np.repeat(np.arange(3), 5) == k].mean(axis=0)).max())
                 for k in range(3))
    print(f"[{time.time() - t0:.0f}s] D0", flush=True)
    return dict(df=float(np.abs(res['lift'][0] - res['K0'][0]).max()), spread=spread)


def job_D(args):
    """The lift m = 15 with class 3's particles (index 2) separated by c_i·10⁻⁶ along G^{1/2}u (u the top eigenvector of S₃, Σc_i = 0), one noise path
    (std, seed, sign: ±1 the antithetic pair): the separation's norm every 5 time units and its growth rate (least squares over t ≥ 50)."""
    std, seed, sign = args
    t0 = time.time()
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, np.inf, atoms, Q3)
    S2 = MN.separation_drift(cv, fprime, fsecond, KE, B_OU, 1.0)[2]
    u = np.linalg.eigh(S2)[1][:, -1]
    labels = np.repeat(np.arange(3), 5); m = 15
    cl = MN.NoisyRealLoop(net, f_vec, xr, xrd, PER, KE, G, SIGMA0, labels, 3, np.inf)
    z0 = cl.lift(atoms)
    for i, c in zip(range(10, 15), (2.0, 1.0, 0.0, -1.0, -2.0)):
        z0[1 + 3 * i:4 + 3 * i] += 1e-6 * c * np.sqrt(G) * u

    def rec(t, z):
        Th = z[1:1 + m * 3].reshape(m, 3)[10:15]
        return np.array([t, np.sqrt(np.sum(((Th - Th.mean(axis=0)) / np.sqrt(G)) ** 2))])
    out = np.array(MN.simulate_noisy(cl, z0, m, T_D, DT_C, B_OU, std ** 2, np.random.default_rng(3000 + seed), int(round(5.0 / DT_C)), rec,
                                     sign=sign))
    msk = out[:, 0] >= 50.0
    rate = float(np.polyfit(out[msk, 0], np.log(out[msk, 1]), 1)[0])
    print(f"[{time.time() - t0:.0f}s] D std {std:g} seed {seed} sign {sign:+g}", flush=True)
    return dict(std=std, seed=seed, sign=sign, t=out[:, 0], sep=out[:, 1], rate=rate)


def job_D_pred(args):
    """The drift along the second-order transient from the lift (Eν² = 1; the drift is linear in Eν² at this order): S_k(t) on the grid, its
    50-unit window averages' largest eigenvalues, the mean error's rms per window, and per std the predicted rate — the least-squares slope of
    log‖u‖ over t ≥ 50 for u̇ = (Eν²S₃(t) − σ₀/k_e)u from S₃'s top eigenvector at the lift — and the std at which it vanishes."""
    t0 = time.time()
    net, atoms = setup()
    cv = AF.CollapsedVariational(net, xr, PER, G, SIGMA0, np.inf, atoms, Q3)
    cl = MN.NoisyRealLoop(net, f_vec, xr, xrd, PER, KE, G, SIGMA0, np.arange(3), 3, np.inf, weights=Q3)
    o = MN.second_order_transient(cl, cv, fprime, B_OU, 1.0, T_D, dt=0.05)
    S = MN.separation_drift_series(cv, o['t'], o['Ee'], o['Eem2'])
    Slift = MN.separation_drift(cv, fprime, fsecond, KE, B_OU, 1.0)
    u0 = np.linalg.eigh(Slift[2])[1][:, -1]
    edges = np.arange(0.0, T_D + 1e-9, 50.0)
    win = [dict(a=a, b=b, lam=[float(np.linalg.eigvalsh(S[(o['t'] >= a) & (o['t'] < b)].mean(axis=0)[k])[-1]) for k in range(3)],
                rmsE=float(np.sqrt(np.mean(o['Ee'][(o['t'] >= a) & (o['t'] < b)] ** 2)))) for a, b in zip(edges[:-1], edges[1:])]

    def pred(std):
        u = u0.copy(); ln = [0.0]; dt = o['t'][1] - o['t'][0]
        for k in range(len(o['t']) - 1):
            A1 = std ** 2 * S[k, 2] - SIGMA0 / KE * np.eye(3); A2 = std ** 2 * S[k + 1, 2] - SIGMA0 / KE * np.eye(3); Am = 0.5 * (A1 + A2)
            k1 = A1 @ u; k2 = Am @ (u + dt / 2 * k1); k3 = Am @ (u + dt / 2 * k2); k4 = A2 @ (u + dt * k3)
            u = u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            nrm = np.linalg.norm(u); ln.append(ln[-1] + np.log(nrm)); u = u / nrm
        ln = np.array(ln); msk = o['t'] >= 50.0
        return float(np.polyfit(o['t'][msk], ln[msk], 1)[0])
    rates = {std: pred(std) for std in STDS_D}
    lo, hi = 0.01, 0.06
    for _ in range(40):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if pred(mid) < 0 else (lo, mid)
    print(f"[{time.time() - t0:.0f}s] D pred", flush=True)
    return dict(win=win, rates=rates, lift=[float(np.linalg.eigvalsh(Slift[k])[-1]) for k in range(3)], std_star=0.5 * (lo + hi))


# ---------------------------------------------------------------------------- report
def report(A, Bres, Cres, Cpred, t0=None, B2=None, D=None):
    es = lambda e: '∞' if np.isinf(e) else ('0 (fixed)' if e == 0 else f'{e:.5f}')
    h, qd = A['h'], A['qd']
    lines = ["# Proposition D‴.30: measurement noise at perfect tracking (`exp_measurement_noise.py`)", "",
             f"The three-bump instance at σ₀ = 0.15σ_crit = {SIGMA0:.5f}, x_m = 2 sin t, the nearest-atom leak toward a fixed (η = 0) or an adaptive prior, at the "
             "lift of the exact representation of q = 1/3 with the prior at its atoms; the measurement y = x + ν in the controller and the adaptive law.", "",
             "## A. The linear theory and the second-order force (no dynamics)", "",
             f"k_e = {KE:g}, N₀ = 1. Q = ⟨(Ĵ̇ + aĴ)ᵀ(Ĵ̇ + aĴ)⟩, a = ∂_xf̂ on the orbit: tr Q = {np.trace(A['Q']):.4f} = {A['parts']['deriv']:.4f} (⟨Ĵ̇ᵀĴ̇⟩) "
             f"{A['parts']['cross']:+.4f} (cross) {A['parts']['slope']:+.4f} (⟨a²ĴᵀĴ⟩). Per mode (the eigenbasis of H): h_j, Q_jj, Q_jj/h_j, and the "
             "function-level shares of the floor, (N₀/2k_e)Q_jj(h_j + η)/(h_j + η + σ₀), for each η.", "",
             "| j | h_j | Q_jj | Q_jj/h_j | " + " | ".join(f"share at η = {es(e)}" for e in ETAS_A) + " |", "|---|---|---|---|" + "---|" * len(ETAS_A)]
    for j in range(len(h)):
        lines.append(f"| {j + 1} | {h[j]:.4e} | {qd[j]:.4e} | {qd[j] / h[j]:.2f} | " + " | ".join(f"{A['rows'][e]['shares'][j]:.4e}" for e in ETAS_A) + " |")
    lines += ["", "| η | floor (N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) | relative to η = ∞ | the weakest mode's parameter variance | its bias-removal rate λ_1^− |",
              "|---|---|---|---|---|"]
    finf = A['rows'][np.inf]['floor']
    for e in ETAS_A:
        r = A['rows'][e]
        lines.append(f"| {es(e)} | {r['floor']:.5e} | {r['floor'] / finf:.4f} | {r['pvar'][0]:.4e} | {r['rate'][0]:.3e} |")
    lines += ["", f"At η = ∞ the floor is (N₀/2k_e)tr Q = {np.trace(A['Q']) / (2 * KE):.5e} whatever H. The second-order force (errors in variables; OU noise of rate "
              f"B = {B_OU:g}), slow time, per unit of the noise's variance, in the eigenbasis of H: " + ", ".join(f"{x:+.3e}" for x in A['fv']) + ". The stationary "
              f"mean deviation's function-level size per unit of the variance squared: fixed prior Σ_jF_j²h_j/(h_j + σ₀)² = {A['bias'][0.0]:.4e}, adaptive prior "
              f"(every η > 0) Σ_jF_j²/h_j = {A['bias'][np.inf]:.4e} (ratio {A['bias'][np.inf] / A['bias'][0.0]:.2f})."]
    if 'trH2' in A:
        xo = {e: A['rows'][e]['floor'] * 2.0 / B_OU / A['bias'][e] for e in (0.0, np.inf)}
        lines += ["", f"The fast parts: tr H² = {A['trH2']:.4f} (the function-level floor) and ⟨|Ĵ|⁴⟩ = {A['J4']:.4f} (the compensation error at the current point), "
                  f"each times N₀/(2k_e). With the OU noise (N₀ = 2Eν²/B) the errors-in-variables bias equals the floor of the parameters at Eν² = "
                  f"{xo[0.0]:.2e} (std {np.sqrt(xo[0.0]):.4f}) under the fixed prior and {xo[np.inf]:.2e} (std {np.sqrt(xo[np.inf]):.4f}) under the adaptive one, "
                  "and dominates above."]
    if 'so' in A:
        so = A['so']
        lines += ["", "The full second-order mean (`second_order_mean`: the periodic covariance of the linearised loop with the OU noise as a state, contracted "
                  "with the second derivatives of the true vector field). Its noise part — the noise's own second moments at fixed parameters — equals the "
                  f"quasi-static force above to {so[np.inf]['rel']:.1e} (relative); the parameters' first-order fluctuation adds a force of the same order in "
                  "Eν² (the compensation's curvature and slope–parameter terms in the error equation, ∇²_ϑΦ̃E[ε_mδϑ] in the law). Function-level biases per unit "
                  "(Eν²)²: of the mean parameters (‖Ĵȳ‖²_ρ, the averaged response) from the noise part and from the full force, of the mean function "
                  "(‖Ĵȳ + ½D²_Θf̂[Σ_ΘΘ]‖²_ρ), and of the exact periodic second-order mean of the linearised loop (mean parameters).", "",
                  "| prior | noise part | full force | mean function | exact periodic mean | periodicity |", "|---|---|---|---|---|---|"]
        for e in (0.0, np.inf):
            r = so[e]
            lines.append(f"| {'fixed' if e == 0 else 'adaptive (η = ∞; the noise part for every η > 0)'} | {r['noise']:.6g} | {r['full']:.6g} | {r['function']:.6g} | {r['exact']:.6g} | {r['per']:.1e} |")
        kg = KES_GROW
        lines += ["", "The adaptive prior's averaged bias per unit (Eν²)² against k_e at the OU rates B, with its local exponents d log/d log k_e between "
                  "successive k_e (→ 2 while k_e ≪ B, saturating for k_e ≫ B).", "",
                  "| B | " + " | ".join(f"k_e = {k:g}" for k in kg) + " | local exponents |", "|---|" + "---|" * (len(kg) + 1)]
        for B, vals in A['grow'].items():
            ex = [np.log(vals[i + 1] / vals[i]) / np.log(kg[i + 1] / kg[i]) for i in range(len(kg) - 1)]
            lines.append(f"| {B:g} | " + " | ".join(f"{v_:.4g}" for v_ in vals) + " | " + ", ".join(f"{x:.2f}" for x in ex) + " |")
        marg = SIGMA0 / KE
        lines += ["", f"Invariance is not attraction: a separation of class k decays at σ₀/k_e = {marg:.3e} without noise and at σ₀/k_e − Eν²λ_max(S_k) with it "
                  "(S_k the averaged second-order drift, real time, `separation_drift`). λ_max(S_k) per unit Eν²: " + ", ".join(f"{x:.4f}" for x in A['sep'])
                  + "; the separations stop decaying at the std " + ", ".join(f"{np.sqrt(marg / x):.3f}" for x in A['sep']) + " (classes 1, 2, 3); at Part C's "
                  "stds the margin is used to " + ", ".join(f"{100 * std ** 2 * max(A['sep']) / marg:.1f}%" for std in STDS_C) + " (the weakest class)."]
    if Bres:
        lines += ["", "## B. The exactly linearised loop with white measurement noise against the averaged prediction", "",
                  "The periodic stationary covariance P(t) (the monodromy and a discrete Lyapunov equation; N₀ = 1), averaged over a period: the function-level floor "
                  "E‖f̂ − f‖²_ρ = ⟨tr(HP_yy)⟩ and the compensation error at the plant's current point ⟨ĴP_yyĴᵀ⟩, against the prediction — the slow (parameter) part "
                  "(N₀/2k_e)Σ_jQ_jj(h_j + η)/(h_j + η + σ₀) plus the fast part (N₀/2k_e)tr H² and (N₀/2k_e)⟨|Ĵ|⁴⟩ respectively.", "",
                  "| k_e | η | E‖f̂ − f‖²_ρ: exact | predicted | ratio | at the current point: exact | predicted | ratio | the slow part alone |",
                  "|---|---|---|---|---|---|---|---|---|"]
        main_ = [r for r in Bres if r.get('slope', True)]
        for r in sorted(main_, key=lambda r: (r['eta'] if np.isfinite(r['eta']) else 1e9, r['ke'])):
            lines.append(f"| {r['ke']:g} | {es(r['eta'])} | {r['exact_orbit']:.5e} | {r['pred_orbit']:.5e} | {r['exact_orbit'] / r['pred_orbit']:.4f} | "
                         f"{r['exact_point']:.5e} | {r['pred_point']:.5e} | {r['exact_point'] / r['pred_point']:.4f} | {r['slow']:.5e} |")
        if 'mode_exact' in main_[0]:
            lines += ["", "The weakest mode's variance (period average, N₀ = 1) against its slow part (N₀/2k_e)Q₁₁(h₁ + η)/(h₁(h₁ + η + σ₀)) plus its fast part "
                      "(N₀/2k_e)h₁, and the ratio of the adaptive (η = ∞) to the fixed prior's variance, exact and averaged.", "",
                      "| k_e | η | exact | slow + fast | ratio |", "|---|---|---|---|---|"]
            for r in sorted(main_, key=lambda r: (r['eta'] if np.isfinite(r['eta']) else 1e9, r['ke'])):
                pr = r['mode_slow'][0] + r['mode_fast'][0]
                lines.append(f"| {r['ke']:g} | {es(r['eta'])} | {r['mode_exact'][0]:.5e} | {pr:.5e} | {r['mode_exact'][0] / pr:.4f} |")
            byk = {(r['ke'], r['eta']): r for r in main_}
            lines += ["", "Adaptive over fixed, the weakest mode: exact " + ", ".join(
                f"{byk[(k, np.inf)]['mode_exact'][0] / byk[(k, 0.0)]['mode_exact'][0]:.2f}" for k in KES_B if (k, 0.0) in byk and (k, np.inf) in byk)
                + f" (k_e = {', '.join(f'{k:g}' for k in KES_B)}); averaged (slow parts) {byk[(KES_B[0], np.inf)]['mode_slow'][0] / byk[(KES_B[0], 0.0)]['mode_slow'][0]:.2f}."]
        ctrl = [r for r in Bres if not r.get('slope', True)]
        if ctrl:
            lines += ["", "The deviation of the floor ratio from 1 (η = 0 and ∞) with the noise input −(k_e + a)ν and in the control with −k_eν (the slope a removed, "
                      "the generator unchanged; its prediction with Q = ⟨Ĵ̇ᵀĴ̇⟩), and the deviation's successive ratios per doubling of k_e (→ 4 for an even "
                      "expansion in 1/k_e, → 2 for a first-order term).", "",
                      "| η | noise input | ratio − 1 at k_e = " + ", ".join(f"{k:g}" for k in KES_B) + " | successive ratios |", "|---|---|---|---|"]
            for e in (0.0, np.inf):
                for lab, rs in (('−(k_e + a)ν', main_), ('−k_eν (control)', ctrl)):
                    dv = [next(r for r in rs if r['ke'] == k and r['eta'] == e) for k in KES_B]
                    dv = [r['exact_orbit'] / r['pred_orbit'] - 1 for r in dv]
                    lines.append(f"| {es(e)} | {lab} | " + ", ".join(f"{x:+.3e}" for x in dv) + " | " + ", ".join(f"{dv[i] / dv[i + 1]:.2f}" for i in range(len(dv) - 1)) + " |")
    if B2:
        lines += ["", "The exact stationary second-order mean of the linearised loop with the OU noise (B = " f"{B_OU:g}" ", the measured regressor): the function-level "
                  "bias of its mean parameters over a period against the averaged bias of the quasi-static force, per unit (Eν²)².", "",
                  "| k_e | prior | averaged | exact | exact/averaged | periodicity |", "|---|---|---|---|---|---|"]
        for r in sorted(B2, key=lambda r: (r['eta'] if np.isfinite(r['eta']) else 1e9, r['ke'])):
            lines.append(f"| {r['ke']:g} | {es(r['eta'])} | {r['averaged']:.6g} | {r['exact']:.6g} | {r['exact'] / r['averaged']:.4f} | {r['per']:.1e} |")
    if Cres:
        lines += ["", f"## C. The real loop with an Ornstein–Uhlenbeck measurement noise (k_e = {KE:g}, B = {B_OU:g}, RK4 dt = {DT_C:g}, {T_C:g} time units)", "",
                  "E‖f̂ − f‖²_ρ at the records (mean ± standard error over the runs) against the linear prediction: the transient covariance of the exactly linearised "
                  "loop with the OU noise plus the second-order mean deviation (their sum, and the two parts).", ""]
        for std in STDS_C:
            for eta in (0.0, np.inf):
                rs = [r for r in Cres if r['kind'] == 'K0' and r['std'] == std and r['eta'] == eta]
                if not rs:
                    continue
                E = np.array([r['err'] for r in rs])
                pr = Cpred[(std, eta)]
                lines += [f"Noise std {std:g} (variance {std ** 2:.1e}), η = {es(eta)}, {len(rs)} runs:", "",
                          "| t | Monte Carlo | prediction | covariance part | mean part |", "|---|---|---|---|---|"]
                for k_, t in enumerate(pr['t']):
                    if k_ == 0:
                        continue
                    lines.append(f"| {t:g} | {E[:, k_].mean():.4e} ± {E[:, k_].std(ddof=1) / np.sqrt(len(rs)):.1e} | {pr['cov'][k_] + pr['mean'][k_]:.4e} | "
                                 f"{pr['cov'][k_]:.4e} | {pr['mean'][k_]:.4e} |")
                lines += [""]
        dd = []
        for std in STDS_C:
            E0 = np.array([r['err'] for r in Cres if r['kind'] == 'K0' and r['std'] == std and r['eta'] == 0.0])
            Ei = np.array([r['err'] for r in Cres if r['kind'] == 'K0' and r['std'] == std and np.isinf(r['eta'])])
            if len(E0) and len(Ei):
                p0, pi = Cpred[(std, 0.0)], Cpred[(std, np.inf)]
                mc = Ei.mean(axis=0) / E0.mean(axis=0) - 1
                pr = (pi['cov'] + pi['mean']) / (p0['cov'] + p0['mean']) - 1
                dd.append(f"std {std:g}: " + ", ".join(f"{100 * mc[k]:+.1f}% ({100 * pr[k]:+.1f}%)" for k in (1, len(pr) // 2, len(pr) - 1))
                          + f" at t = {p0['t'][1]:g}, {p0['t'][len(pr) // 2]:g}, {p0['t'][-1]:g}")
        if dd:
            lines += ["The adaptive prior's E‖f̂ − f‖²_ρ over the fixed prior's, Monte Carlo (predicted): " + "; ".join(dd) + ". All four cases run on the same "
                      f"{N_MC} noise paths (the seeds 1000 + run, scaled by the std), so their agreements are not independent.", ""]
        lift = [r for r in Cres if r['kind'] == 'lift']
        if lift:
            r = lift[0]
            twin = [x for x in Cres if x['kind'] == 'K0' and x['seed'] == r['seed'] and x['std'] == r['std'] and x['eta'] == r['eta']][0]
            lines.append(f"The lift m = 15 (counts 5, 5, 5) driven by the same noise path as the K₀-atom run (seed {r['seed']}, std {r['std']:g}, η = {es(r['eta'])}): "
                         f"the largest difference of the stored f̂ (float32) on the grid over the records {np.abs(r['fh'].astype(float) - twin['fh'].astype(float)).max():.1e}, "
                         f"of E‖f̂ − f‖² {np.abs(r['err'] - twin['err']).max():.1e} — the noise is coherent: a replicated start stays replicated (in float64: Part D).")
    if D:
        d0 = D['D0']
        lines += ["", "## D. Invariance is not attraction", "",
                  f"The lift and the K₀-atom run on one noise path (std 0.01, η = ∞, 25 time units) in float64: the largest difference of f̂ on the grid "
                  f"{d0['df']:.1e} (a mean of five against a weighted atom: rounding), the lift's largest within-class spread {d0['spread']:.1e} (identical "
                  "arithmetic per particle, whether or not the replicated state attracts).", "",
                  f"The lift m = 15 with class 3's particles separated by 10⁻⁶·(2, 1, 0, −1, −2) along S₃'s top eigenvector, {N_D} antithetic pairs of noise "
                  f"paths per std (ν and −ν), {T_D:g} time units: the separation's growth rate (least squares on log‖separation‖ over t ≥ 50) — its pair "
                  "means cancel the first-order part — and λ_eff = (rate + σ₀/k_e)/Eν², against the lift's drift (−σ₀/k_e + λ_max(S₃)Eν²) and against the "
                  "drift along the second-order transient from the lift (the slope predicted for u̇ = (Eν²S₃(t) − σ₀/k_e)u).", ""]
        pr = D.get('pred')
        if pr is not None:
            lines += ["The drift along the transient (Eν² = 1): per 50-unit window, the mean error's rms and the largest eigenvalues of the window-averaged "
                      "S_k (classes 1, 2, 3); at the lift: " + ", ".join(f"{v:.4f}" for v in pr['lift']) + ".", "",
                      "| window | rms E[e] | λ_max(S₁), λ_max(S₂), λ_max(S₃) |", "|---|---|---|"]
            for w_ in pr['win']:
                lines.append(f"| [{w_['a']:g}, {w_['b']:g}) | {w_['rmsE']:.3f} | " + ", ".join(f"{v:.4f}" for v in w_['lam']) + " |")
            lines += ["", f"The separation's rate vanishes at the std {np.sqrt(SIGMA0 / KE / max(pr['lift'])):.4f} by the lift's drift and at {pr['std_star']:.4f} "
                      f"by the transient's over t ∈ [50, {T_D:g}].", ""]
        lines += ["| std | lift's prediction | transient's prediction | measured (pair means ± s.e.) | λ_eff measured | λ_eff: lift; transient |", "|---|---|---|---|---|---|"]
        lm = max(A['sep']) if 'sep' in A else np.nan
        fit = []
        for std in STDS_D:
            rs = [r for r in D['runs'] if r['std'] == std]
            if not rs:
                continue
            seeds = sorted(set(r['seed'] for r in rs))
            pm = np.array([np.mean([r['rate'] for r in rs if r['seed'] == sd]) for sd in seeds])
            se = pm.std(ddof=1) / np.sqrt(len(pm))
            le = (pm + SIGMA0 / KE) / std ** 2
            tp = pr['rates'][std] if pr is not None else np.nan
            fit.append((std ** 2, le.mean(), le.std(ddof=1) / np.sqrt(len(le))))
            lines.append(f"| {std:g} | {-SIGMA0 / KE + lm * std ** 2:+.3e} | {tp:+.3e} | {pm.mean():+.4e} ± {se:.1e} | {le.mean():.4f} ± {fit[-1][2]:.4f} | "
                         f"{lm:.4f}; {(tp + SIGMA0 / KE) / std ** 2:.4f} |")
        if len(fit) >= 3:                                                          # λ_eff = c₀ + c₁Eν², weighted by the standard errors
            v, le, se = (np.array(c) for c in zip(*fit))
            (c1, c0), cov = np.polyfit(v, le, 1, w=1.0 / se, cov='unscaled')
            root = (-c0 + np.sqrt(c0 ** 2 + 4 * c1 * SIGMA0 / KE)) / (2 * c1) if c1 != 0 else SIGMA0 / KE / c0
            lines += ["", f"The weighted fit λ_eff = c₀ + c₁Eν²: c₀ = {c0:.4f} ± {np.sqrt(cov[1, 1]):.4f} (the small-noise limit; the lift's drift "
                      f"{lm:.4f}, the transient's {(pr['rates'][STDS_D[0]] + SIGMA0 / KE) / STDS_D[0] ** 2 if pr is not None else np.nan:.4f} at std "
                      f"{STDS_D[0]:g}), c₁ = {c1:.2f} ± {np.sqrt(cov[0, 0]):.2f} (the next order); the measured rate vanishes at the std {np.sqrt(root):.4f}."]
    lines += ["", "Numbers are produced by `experiments/exp_measurement_noise.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'measurement_noise.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'measurement_noise_raw.npz')
    jobs_B = [(ke, eta) for ke in KES_B for eta in ETAS_B] + [(ke, eta, 'noslope') for ke in KES_B for eta in (0.0, np.inf)]
    jobs_B2 = [(ke, eta) for ke in KES_B2 for eta in (0.0, np.inf)]
    jobs_D = [(std, seed, sign) for std in STDS_D for seed in range(N_D) for sign in (1.0, -1.0)]

    def rep_(save):
        report(save['A'].item(), list(save['B']), list(save['C']), save['P'].item(), t0,
               B2=list(save['B2']) if 'B2' in save else None, D=save['D'].item() if 'D' in save else None)
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            rep_({k: z[k] for k in z.files})
        return
    if any(a in sys.argv for a in ('--only-A', '--only-B', '--only-D')):    # the part recomputed, merged with the stored ones
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        if '--only-A' in sys.argv:
            save['A'] = np.array(part_A(), dtype=object)
        with mp.Pool(N_WORKERS) as pool:
            if '--only-B' in sys.argv:
                save['B'] = np.array(pool.map(job_B, jobs_B, chunksize=1), dtype=object)
                save['B2'] = np.array(pool.map(job_B2, jobs_B2, chunksize=1), dtype=object)
            if '--only-D' in sys.argv:
                pr = pool.apply_async(job_D_pred, (None,))
                runs = pool.map(job_D, jobs_D, chunksize=1)
                save['D'] = np.array(dict(D0=job_D0(None), runs=runs, pred=pr.get()), dtype=object)
        save_raw(raw, **save)
        rep_(save)
        return
    A = part_A()
    with mp.Pool(N_WORKERS) as pool:
        Bres = pool.map(job_B, jobs_B, chunksize=1)
        B2 = pool.map(job_B2, jobs_B2, chunksize=1)
        save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bres, dtype=object), B2=np.array(B2, dtype=object))
        jobs = [('K0', eta, std, seed) for std in STDS_C for eta in (0.0, np.inf) for seed in range(N_MC)] + [('lift', np.inf, STDS_C[1], 0)]
        Cres = pool.map(job_C, jobs, chunksize=1)
        Dd = dict(D0=job_D0(None), runs=pool.map(job_D, jobs_D, chunksize=1), pred=job_D_pred(None))
    Cpred = {(std, eta): predictions_C(std, eta) for std in STDS_C for eta in (0.0, np.inf)}
    save = dict(A=np.array(A, dtype=object), B=np.array(Bres, dtype=object), B2=np.array(B2, dtype=object), C=np.array(Cres, dtype=object),
                P=np.array(Cpred, dtype=object), D=np.array(Dd, dtype=object))
    save_raw(raw, **save)
    rep_(save)


if __name__ == '__main__':
    main()
