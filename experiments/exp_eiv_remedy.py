"""
Proposition D‴.31: the adaptive prior is a washout leak, and the desired regressor turns the noise's persistent force into a least-squares one
==============================================================================================================================================

The instance of Proposition D‴.30 (`exp_measurement_noise.py`: the three bumps, σ₀ = 0.15σ_crit, x_m = 2 sin t, the lift of the exact
representation of q = 1/3 with the prior at its atoms) and three versions of the controller and the adaptive law under the measurement
y = x + ν: the measured regressor (MR, Proposition D‴.30: ϑ̇ ∋ G∇_ϑ((y − x_m)Φ̃(y; ϑ)), the compensation f̂(y)), the desired regressor (DR: the
regressor at the reference, ∇_ϑΦ̃(x_m; ϑ), the compensation f̂(y)) and the desired compensation adaptation law (DR*: both at the reference).

A. The second-order force and its stationary responses (no dynamics), at k_e = 25, 50, 100, 200 and the Ornstein–Uhlenbeck rates B = 100,
   1000: the function-level stationary bias per unit (Eν²)² under the fixed prior, Fᵀ(H + σ₀I)⁻¹H(H + σ₀I)⁻¹F, and under the adaptive prior,
   FᵀH⁻¹F — the squared norm of the smallest target perturbation whose learning force is F — for the force of the noise's own second
   moments (the quasi-static `eiv_force`) and for the full second-order force (`second_order_mean`: the parameters' fluctuation
   included; the bias of the mean parameters and of the mean function); for DR and DR* the bound ‖g‖²_ρ, g = k_eE[e]; the rescaled
   projected mean tracking error k_e²‖P E[e]‖²_ρ = F₁ᵀH⁻¹F₁ under the adaptive prior and the mean tracking error itself; the noise level at
   which the bias equals the floor of the parameters.
B. The stationary identity RH + B = I along the prior family dz/ds = −η(z − T̂) − κ(z − z_a) (η = σ₀; κ/η = 0, 0.1, 1, 10, and the fixed
   prior): the leak's zero-frequency gain S₀, the weakest mode's response and bias factor, the errors-in-variables bias of MR and DR, the
   bias of an anchor misplaced by a unit along the weakest and along the strongest mode, and the identity's residual (from the block); and
   priors whose zero-frequency leak S₀ does not commute with H, R and B again from their dynamic blocks (`anchored_block`: a prior fixed on
   the amplitudes or on the centres and adaptive at η = σ₀ elsewhere, one active face of the projected adaptive prior at η = ∞, a matrix
   anchor gain K = ηP_amp, a random anchor gain with K + Kᵀ ⪰ 0), against (H + S₀)⁻¹ and (H + S₀)⁻¹S₀ with S₀ by elimination: RH + B = I
   holds and HR + B = I does not (round 22 of the review, ID-5: the identity derived from the block, not from the definitions).
C. The real loop (RK4, dt = 10⁻³, k_e = 25, B = 100, the noise's std 0.01, η = ∞, N_RUNS runs of T_C time units): MR, DR and DR* driven by
   the same noise paths (common random numbers). E‖f̂ − f‖²_ρ against the prediction — the transient covariance of the linearised loop
   (MR and DR share it: their linearisations coincide) plus the second-order mean deviation — and the run mean of f̂_MR − f̂_DR on the
   orbit's grid against its prediction Ĵ(ȳ_MR − ȳ_DR) (the first-order responses cancel path by path).
D. The separations under the desired regressor (Corollary D‴.33): the separations' second-order drift S_k at the lift for MR, DR and DR*
   (with the regressor at x_m only the mean error's term remains), the noise thresholds under the adaptive prior and the cohesion leak, the
   drift at DR's and DR*'s stationary points (k_eE[e] = (I − P)g) and along DR's second-order transient; the real loop — the lift m = 15
   with class 3 separated along MR's top eigenvector, N_D antithetic pairs per law and std (0.05, 0.1) on common paths — the separation's
   rate under MR and DR.

Output: results/chaos/eiv_remedy.md and eiv_remedy_raw.npz (ER_WORKERS, ER_RUNS, ER_TC; `--from-raw` regenerates the md; `--only-AB`
recomputes Parts A and B and merges them, `--only-B` Part B, `--only-D` Part D, `--only-Dpred` Part D's predictions; ER_ND the pairs of
Part D).
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
from src.core import adaptive_floquet as AF
from src.core import measurement_noise as MN
import exp_measurement_noise as X                                   # the instance of Proposition D‴.30

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
VARIANTS = ('measured', 'reference', 'reference_both')
NAMES = {'measured': 'MR', 'reference': 'DR', 'reference_both': 'DR*'}
LAW = {'measured': ('measured', 'measured'), 'reference': ('reference', 'measured'), 'reference_both': ('reference', 'reference')}
KES_A = (25.0, 50.0, 100.0, 200.0)
BS_A = (100.0, 1000.0)
KAPPAS_B = (0.0, 0.1, 1.0, 10.0)                                   # κ/η at η = σ₀; then the fixed prior
KE, B_OU, DT_C, STD_C = X.KE, X.B_OU, X.DT_C, 0.01
N_RUNS = int(os.environ.get('ER_RUNS', 12))
T_C = float(os.environ.get('ER_TC', 300.0))
REC_C = 25.0
N_WORKERS = int(os.environ.get('ER_WORKERS', 3))
STDS_D = (0.05, 0.1)
N_D = int(os.environ.get('ER_ND', 4))                                 # antithetic pairs per (law, std)


def collapsed():
    net, atoms = X.setup()
    return net, atoms, AF.CollapsedVariational(net, X.xr, X.PER, X.G, X.SIGMA0, np.inf, atoms, X.Q3)


# ---------------------------------------------------------------------------- A
def job_A(args):
    """The full second-order mean (`second_order_mean`) of variant v at (k_e, B) under the adaptive (η = ∞) or the fixed prior (η = 0): its
    biases (`mean_biases`) per unit (Eν²)²."""
    ke, B, v, eta = args
    t0 = time.time()
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, ke, X.G, X.SIGMA0, np.arange(3), 3, eta, weights=X.Q3, regressor=reg, compensation=comp)
    out = MN.second_order_mean(cl, cv, X.fprime, B, 1.0)
    res = MN.mean_biases(H, out, X.SIGMA0 if eta == 0 else None)
    print(f"[{time.time() - t0:.0f}s] A {NAMES[v]} ke {ke:g} B {B:g} eta {eta:g}", flush=True)
    return dict(ke=ke, B=B, v=v, eta=eta, per=out['per_err'], **res)


def part_A(pool=None):
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    jobs = [(ke, B, v, eta) for ke in KES_A for B in BS_A for v in VARIANTS for eta in (np.inf, 0.0)]
    full = {j[:4]: r for j, r in zip(jobs, (pool.map(job_A, jobs, chunksize=1) if pool is not None else map(job_A, jobs)))}
    rows = []
    for ke in KES_A:
        for B in BS_A:
            for v in VARIANTS:
                F, g = MN.eiv_force(cv, X.fprime, X.fsecond, ke, B, 1.0, variant=v, return_g=True)
                Fs = ke * F                                           # slow time, per unit of the noise's variance
                Fl = ke * MN.eiv_force(cv, X.fprime, X.fsecond, ke, B, 1.0, variant='reference' if v == 'measured' else v)   # ⟨Ĵᵀg⟩
                F1 = Fs - Fl                                          # the regressor's errors in variables proper (MR only)
                adapt = MN.equivalent_perturbation(H, Fs)
                fixed = MN.response_size(H, Fs, X.SIGMA0)
                merr = MN.equivalent_perturbation(H, F1)              # ‖P(g − Ĵȳ_∞)‖²_ρ under the adaptive prior: F₁ᵀH⁻¹F₁
                floor = MN.floor_closed_form(H, Q, X.SIGMA0, np.inf, 2.0 / B, ke)[0]   # the parameters' floor per unit Eν² (N₀ = 2Eν²/B)
                fa, ff = full[(ke, B, v, np.inf)], full[(ke, B, v, 0.0)]
                rows.append(dict(ke=ke, B=B, v=v, adapt=adapt, fixed=fixed, g2=float(np.mean(g ** 2)), merr=merr, floor=floor,
                                 cross=floor / fa['full'], full_adapt=fa, full_fixed=ff))
    print("[A]", flush=True)
    return dict(rows=rows, h=np.linalg.eigvalsh(H))


# ---------------------------------------------------------------------------- B
def part_B():
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    h, V = np.linalg.eigh(H)
    F = {v: KE * MN.eiv_force(cv, X.fprime, X.fsecond, KE, B_OU, 1.0, variant=v) for v in ('measured', 'reference')}
    eta = X.SIGMA0
    cases = [(eta, k * eta) for k in KAPPAS_B] + [(0.0, 0.0)]
    rows = []
    for e_, kap in cases:
        R, Bm = MN.stationary_operators(H, X.SIGMA0, e_, kap)
        S0 = MN.zero_frequency_leak(X.SIGMA0, e_, kap)
        resid = float(np.abs(H @ R + Bm - np.eye(H.shape[0])).max())
        r1 = float(V[:, 0] @ R @ V[:, 0]); b1 = float(V[:, 0] @ Bm @ V[:, 0])
        eiv = {v: float((R @ F[v]) @ H @ (R @ F[v])) for v in F}
        anc = [float((Bm @ V[:, j]) @ H @ (Bm @ V[:, j])) for j in (0, len(h) - 1)]
        rows.append(dict(eta=e_, kappa=kap, S0=S0, resid=resid, resid_RH=float(np.abs(R @ H + Bm - np.eye(H.shape[0])).max()), r1=r1, b1=b1,
                         eiv=eiv, anchor=anc))
    # priors whose zero-frequency leak does not commute with H, R and B from their dynamic blocks (whitened coordinates (k, l) ↦ 3k + l:
    # l = 0 the amplitude, 1 the centre, 2 the width), against the definitions with S₀ by elimination (round 22, ID-5)
    d = H.shape[0]
    mask = lambda idx: np.isin(np.arange(d), idx)
    Mr = np.random.default_rng(7).standard_normal((d, d))
    Kr = eta * (Mr @ Mr.T / np.linalg.norm(Mr @ Mr.T, 2) + 0.5 * (Mr - Mr.T) / np.linalg.norm(Mr - Mr.T, 2))   # K + Kᵀ ⪰ 0, not symmetric
    cases = [('fixed on the amplitudes, η = σ₀ elsewhere', eta, None, mask([0, 3, 6])),
             ('fixed on the centres, η = σ₀ elsewhere', eta, None, mask([1, 4, 7])),
             ("one active face (atom 1's amplitude), η = ∞", np.inf, None, mask([0])),
             ('anchor gain K = ηP_amp', eta, eta * np.diag(mask([0, 3, 6]).astype(float)), None),
             ('a random anchor gain, K + Kᵀ ⪰ 0, ‖K‖ ≈ η', eta, Kr, None)]
    gen_rows = []
    I = np.eye(d)
    for name, e_, K, fixed in cases:
        R, Bm = MN.anchored_operators(H, X.SIGMA0, e_, K, fixed)
        S0m = MN.anchored_leak(H, X.SIGMA0, e_, K, fixed)
        R0, B0 = MN.zero_frequency_operators(H, S0m)
        per_mode = np.abs(h * np.einsum('ij,ik,kj->j', V, R, V) + np.einsum('ij,ik,kj->j', V, Bm, V) - 1).max()
        eiv = {v: float((R @ F[v]) @ H @ (R @ F[v])) for v in F}
        gen_rows.append(dict(name=name, RH=float(np.abs(R @ H + Bm - I).max()), HR=float(np.abs(H @ R + Bm - I).max()), per_mode=float(per_mode),
                             comm=float(np.abs(H @ S0m - S0m @ H).max()), eiv=eiv,
                             vs_def=float(max(np.abs(R - R0).max() / np.abs(R0).max(), np.abs(Bm - B0).max())),
                             dim=int(MN.anchored_block(H, X.SIGMA0, e_, K, fixed)[0].shape[0])))
    print("[B]", flush=True)
    return dict(rows=rows, h=h, general=gen_rows, n_weak=int(np.sum(h < X.SIGMA0)))


# ---------------------------------------------------------------------------- C
def job_C(args):
    v, seed = args
    t0 = time.time()
    net, atoms = X.setup()
    grid = 2 * np.sin(np.arange(X.N_GRID) * X.PER / X.N_GRID); Fg = X.f_vec(grid)
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, KE, X.G, X.SIGMA0, np.arange(3), 3, np.inf, weights=X.Q3,
                          regressor=reg, compensation=comp)
    z0 = cl.lift(atoms)

    def rec(t, z):
        Th = z[1:1 + 9].reshape(3, 3)
        fh = np.array([cl.f_hat(np.array([x]), Th)[0] for x in grid])
        return np.concatenate([[np.mean((fh - Fg) ** 2)], fh - Fg])
    out = np.array(MN.simulate_noisy(cl, z0, 3, T_C, DT_C, B_OU, STD_C ** 2, np.random.default_rng(5000 + seed), int(round(REC_C / DT_C)), rec))
    print(f"[{time.time() - t0:.0f}s] C {NAMES[v]} seed {seed}", flush=True)
    return dict(v=v, seed=seed, err=out[:, 0], dev=out[:, 1:])


def predictions_C():
    """The linear prediction at the records for each variant: the transient covariance part tr(HP_yy) (the linearised loop with the OU noise,
    MR and DR identical; DR* with the error channel −(k_e − a) and the noise input −k_e), the mean part ȳᵀHȳ, and the mean deviation Ĵȳ on
    the grid."""
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    var = STD_C ** 2; N0 = 2 * var / B_OU
    tev = np.arange(0.0, T_C + 1e-9, REC_C)
    tg = np.arange(X.N_GRID) * X.PER / X.N_GRID
    Jg = np.array([cv.J(s)[0] for s in tg])                              # the whitened Jacobian on the grid (x = 2 sin t_k)
    out = {}
    for v in VARIANTS:
        comp = LAW[v][1]
        gen = lambda t: cv.generator(t, KE)
        if comp == 'reference':
            gen = MN.reference_compensation(gen, cv, X.fprime, KE)
        D = cv.dim
        b0 = MN.collapsed_noise_input(cv, X.fprime, KE, compensation=comp)
        Ma, cfun = MN.ou_augmented(gen, b0, D, B_OU, N0)
        P0 = np.zeros((D + 1, D + 1)); P0[-1, -1] = var
        tt, Ps = MN.transient_covariance(Ma, cfun, D + 1, N0, T_C, tev, P0=P0)
        cov = np.array([np.trace(H @ P[1:1 + cv.d, 1:1 + cv.d]) for P in Ps])
        F = KE * MN.eiv_force(cv, X.fprime, X.fsecond, KE, B_OU, var, variant=v)
        ys = np.array([MN.mean_response(H, F, X.SIGMA0, np.inf, t / KE) for t in tev])
        out[v] = dict(t=tev, cov=cov, mean=np.einsum('ti,ij,tj->t', ys, H, ys), dev=ys @ Jg.T)
    return out


# ---------------------------------------------------------------------------- report
# ---------------------------------------------------------------------------- D (Corollary D‴.33)
KES_D = (25.0, 50.0, 100.0, 200.0, 400.0, 800.0)
N_LONG = int(os.environ.get('ER_NLONG', 82000))                      # periods of the long transient (≈ 5.2·10⁵ time units)
N_GRID = 504                                                         # grid points per period of the long transient (0.3% at the stationary point)
DT_TR = 0.0125                                                       # the transient over the real loop's 300 time units
REC_D = 5.0                                                          # the real loop's record interval (the rate's least squares over t ≥ 50)
TARGETS_D = {'instance': X.BUMPS,
             'two opposite bumps': [(1.0, -1.0, 0.5), (-0.8, 1.0, 0.45)],
             'three same-sign bumps': [(0.8, -1.5, 0.5), (0.8, 0.0, 0.5), (0.8, 1.5, 0.5)],
             'wide bumps (S = 1)': [(1.0, -1.5, 1.0), (-0.7, 0.0, 1.0), (0.8, 1.5, 1.0)],
             'bumps at the orbit\'s edge': [(1.0, -2.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 2.5, 0.5)],
             'two close same-sign bumps': [(0.8, -0.4, 0.5), (0.8, 0.4, 0.5)]}
EXTRA_D = [('reference', std, k, 'top') for std in STDS_D for k in range(3)] + [('reference_both', std, 1, 'top') for std in STDS_D]
FAMILIES = {'measured': [(2, 'u')], 'reference': [(2, 'u'), (0, 'top'), (1, 'top'), (2, 'top')], 'reference_both': [(1, 'top')]}


def target(bumps):
    """The exact representation of q = 1/K₀ for a sum of Gaussian bumps (A, C, S), with f′ and f″ in closed form."""
    net = X.L.SaturatedLiftedNetwork(1, **X.NET_KW)
    K0 = len(bumps)
    q = np.ones(K0) / K0
    A = np.array([b_[0] for b_ in bumps]); C = np.array([b_[1] for b_ in bumps]); S = np.array([b_[2] for b_ in bumps])
    atoms = net.unconstrained((A / q)[:, None], C[:, None], S)
    cv = AF.CollapsedVariational(net, X.xr, X.PER, X.G, X.SIGMA0, np.inf, atoms, q)
    fp = lambda x: sum(-A_ * (x - C_) / S_ ** 2 * np.exp(-(x - C_) ** 2 / (2 * S_ ** 2)) for A_, C_, S_ in bumps)
    fpp = lambda x: sum(A_ * ((x - C_) ** 2 / S_ ** 4 - 1 / S_ ** 2) * np.exp(-(x - C_) ** 2 / (2 * S_ ** 2)) for A_, C_, S_ in bumps)
    return net, atoms, cv, fp, fpp


def lift_drift(v, ke=KE, B=B_OU, cv=None, fp=X.fprime, fpp=X.fsecond, N=512):
    cv = collapsed()[2] if cv is None else cv
    return MN.separation_drift(cv, fp, fpp, ke, B, 1.0, N=N, variant=v)


def directions():
    """The real loop's start directions (whitened unit vectors): MR's top eigenvector of S₃ at the lift ('u', Part D of
    `exp_measurement_noise.py`) and each law's own top eigenvector of S_k ('top')."""
    S = {v: lift_drift(v) for v in VARIANTS}
    out = {('measured', 2, 'u'): np.linalg.eigh(S['measured'][2])[1][:, -1]}
    out[('reference', 2, 'u')] = out[('measured', 2, 'u')]
    for v in ('reference', 'reference_both'):
        for k in range(3):
            out[(v, k, 'top')] = np.linalg.eigh(S[v][k])[1][:, -1]
    return out, S


def pred_lift():
    """The drift at the lift for MR, DR, DR* (extreme eigenvalues per class, the Rayleigh quotients along the start directions); the
    scans in k_e at B = 4k_e and at B = 100 (the largest eigenvalue per law, ‖S_k‖ against (ii)'s bound (Eν²/2k_e)‖∂²_xf‖_∞⟨‖M_k‖⟩, and the
    k_e-uniform floor of the threshold std that the bound gives); and DR's (with DR*'s and MR's largest) drift at the lift for other
    exactly represented targets, with the centre entries."""
    t0 = time.time()
    net, atoms, cv = collapsed()
    dirs, S = directions()
    lift = {v: dict(lo=[float(np.linalg.eigvalsh(sk)[0]) for sk in S[v]], hi=[float(np.linalg.eigvalsh(sk)[-1]) for sk in S[v]]) for v in VARIANTS}
    ts = np.arange(4096) * X.PER / 4096
    xs = 2 * np.sin(ts)
    fpp_inf = float(np.abs(X.fsecond(xs)).max())
    g = np.sqrt(X.G)
    tm = np.arange(512) * X.PER / 512
    Mn = np.array([np.mean([np.linalg.norm(g[:, None] * MN._feature_hessian(net, 2 * np.sin(s_), th) * g[None, :], 2) for s_ in tm])
                   for th in cv.atoms])
    scan = {}
    for ke in KES_D:
        for Bk in ('4ke', 100.0):
            Bv = 4 * ke if Bk == '4ke' else Bk
            for v in VARIANTS:
                Sv = lift_drift(v, ke, Bv, cv)
                scan[(ke, Bk, v)] = dict(lam=float(max(np.linalg.eigvalsh(sk)[-1] for sk in Sv)),
                                         lam_lo=float(min(np.linalg.eigvalsh(sk)[0] for sk in Sv)),
                                         ratio=float(max(np.linalg.norm(sk, 2) / (fpp_inf * Mn[k] / (2 * ke)) for k, sk in enumerate(Sv))))
    floor = float(np.sqrt(2 * X.SIGMA0 / (fpp_inf * Mn.max())))
    targets = {}
    for name, bumps in TARGETS_D.items():
        net_t, atoms_t, cv_t, fp_t, fpp_t = target(bumps)
        r = {}
        for v in VARIANTS:
            Sv = MN.separation_drift(cv_t, fp_t, fpp_t, KE, B_OU, 1.0, N=512, variant=v)
            r[v] = dict(lo=[float(np.linalg.eigvalsh(sk)[0]) for sk in Sv], hi=[float(np.linalg.eigvalsh(sk)[-1]) for sk in Sv],
                        centre=[float(sk[1, 1]) for sk in Sv])
        targets[name] = dict(bumps=[list(b_) for b_ in bumps], **r)
    print(f"[{time.time() - t0:.0f}s] D lift", flush=True)
    return dict(lift=lift, scan=scan, floor=floor, fpp_inf=fpp_inf, Mnorm=Mn.tolist(), targets=targets, w=X.Q3.tolist(),
                dirs={f'{k[0]}:{k[1]}{k[2]}': v.tolist() for k, v in dirs.items()})


def pred_stationary(args):
    """The drift at the stationary point (the exact periodic second-order mean, `second_order_mean(periodic=True)`, referee round 23,
    D3v-1) under the adaptive (η = ∞) or the fixed prior (η = 0), for DR or DR*: the mean error's rms (total, noise part, parameters' part;
    its tangent and orthogonal parts) and the drift's extreme eigenvalues per class from each part — the noise part of E[e], its parameters'
    part, the cross term — and in total; with the averaged noise-only balance k_eE[e] = (I − P)g (the round-22 statement) for comparison."""
    v, eta = args
    t0 = time.time()
    net, atoms, cv = collapsed()
    H, Q = MN.noise_matrices(cv, X.fprime)
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, KE, X.G, X.SIGMA0, np.arange(3), 3, eta, weights=X.Q3, regressor=reg, compensation=comp)
    o = MN.second_order_mean(cl, cv, X.fprime, B_OU, 1.0, N=504, periodic=True)
    t = o['t']; J = o['J']
    parts = dict(noise=o['m_noise'][:, 0], param=o['m_param'][:, 0])
    Eall = parts['noise'] + parts['param']
    z = np.zeros_like(Eall)
    eig = lambda Sm: [[float(e[0]), float(e[-1])] for e in (np.linalg.eigvalsh(sk) for sk in Sm.mean(axis=0))]
    Sn = MN.separation_drift_series(cv, t, parts['noise'], z); Sp = MN.separation_drift_series(cv, t, parts['param'], z)
    C = MN.separation_cross_term(cv, t, o['Eeth'])
    tang = lambda e: float(np.sqrt((J.T @ e / len(e)) @ np.linalg.solve(H, J.T @ e / len(e))))
    rms = lambda e: float(np.sqrt(np.mean(e ** 2)))
    F, g = MN.eiv_force(cv, X.fprime, X.fsecond, KE, B_OU, 1.0, N=504, variant=v, return_g=True)
    Ea = (g - J @ np.linalg.solve(H, KE * F)) / KE
    out = dict(v=v, eta=eta, rms=rms(Eall), rms_noise=rms(parts['noise']), rms_param=rms(parts['param']), tang=tang(Eall),
               orth=float(np.sqrt(max(rms(Eall) ** 2 - tang(Eall) ** 2, 0.0))), tang_param=tang(parts['param']),
               eig_noise=eig(Sn), eig_param=eig(Sp), eig_cross=eig(C), eig_tot=eig(Sn + Sp + C), eig_nocross=eig(Sn + Sp),
               rms_avg=rms(Ea), eig_avg=eig(MN.separation_drift_series(cv, t, Ea, z)), rms_lift=rms(g / KE),
               cross_norm=float(max(np.linalg.norm(c_, 2) for c_ in C.mean(axis=0))))
    print(f"[{time.time() - t0:.0f}s] D stationary {NAMES[v]} eta {eta:g}", flush=True)
    return out


def pred_long(v):
    """The second-order transient from the lift over N_LONG periods (`long_transient`, referee round 23, D3v-2) for DR or DR*."""
    t0 = time.time()
    net, atoms, cv = collapsed()
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, KE, X.G, X.SIGMA0, np.arange(3), 3, np.inf, weights=X.Q3, regressor=reg, compensation=comp)
    o = MN.long_transient(cl, cv, X.fprime, B_OU, 1.0, N_LONG, N=N_GRID,
                          keep=np.unique(np.concatenate([np.arange(48), np.rint(np.logspace(np.log10(48), np.log10(N_LONG - 1), 240)).astype(int)])))
    lt = o['lam_tot'].astype(float); ln_ = o['lam_noise'].astype(float)
    first = lambda lam: [float(o['t'][np.argmax(lam[:, k] > 0)]) if np.any(lam[:, k] > 0) else None for k in range(lam.shape[1])]
    last_neg = lambda lam: [float(o['t'][np.where(lam[:, k] < 0)[0][-1]]) if np.any(lam[:, k] < 0) else None for k in range(lam.shape[1])]
    imax = lt.argmax(axis=0)
    rec = [dict(n=int(n), t=float(n * X.PER), eig_tot=[[float(e[0]), float(e[-1])] for e in np.linalg.eigvalsh(r_['S_noise'] + r_['S_param'] + r_['C'])],
                eig_noise=[[float(e[0]), float(e[-1])] for e in np.linalg.eigvalsh(r_['S_noise'])],
                eig_param=[[float(e[0]), float(e[-1])] for e in np.linalg.eigvalsh(r_['S_param'])],
                eig_cross=[[float(e[0]), float(e[-1])] for e in np.linalg.eigvalsh(r_['C'])],
                rms=r_['rms'], rms_noise=r_['rms_noise'], rms_param=r_['rms_param'], S=(r_['S_noise'] + r_['S_param'] + r_['C']).tolist())
           for n, r_ in zip(o['keep'], o['rec'])]
    Sp = o['S_per']
    out = dict(v=v, n_periods=N_LONG, N=N_GRID, first_pos=first(lt), first_pos_noise=first(ln_), last_neg=last_neg(lt),
               max_tot=[float(lt[imax[k], k]) for k in range(3)], t_max=[float(o['t'][imax[k]]) for k in range(3)],
               max_noise=[float(ln_.max(axis=0)[k]) for k in range(3)], t_max_noise=[float(o['t'][ln_.argmax(axis=0)[k]]) for k in range(3)],
               min_tot=[float(lt[:, k].min()) for k in range(3)], rec=rec, lam_tot=o['lam_tot'], lam_noise=o['lam_noise'],
               stat=dict(eig_tot=[[float(e[0]), float(e[-1])] for e in np.linalg.eigvalsh(Sp['S_noise'] + Sp['S_param'] + Sp['C'])],
                         rms=o['rms_per'], rms_noise=o['rms_per_noise'], rms_param=o['rms_per_param']))
    print(f"[{time.time() - t0:.0f}s] D long transient {NAMES[v]}", flush=True)
    return out


def pred_runs(v):
    """The second-order prediction of the real loop's rate estimator for law v along each start direction used in Part D (referee round 23,
    D3v-6): along the transient from the lift (`second_order_transient`, dt = DT_TR, over T_D) E log(‖u(t)‖/‖u(0)‖) = −σ₀t/k_e +
    Eν²[∫₀ᵗθᵀ(S_k(τ) + C_k(τ))θdτ + V(t)] — S_k the mean error's term (and for MR the regressor's E[ε_m²] term), C_k the cross term, V the
    variance term (`separation_variance`) — and its least-squares slope over the record times t ≥ 50, per unit Eν²: the noise-induced part of
    the predicted rate, split into the three terms, against the lift's θᵀS_kθ."""
    t0 = time.time()
    net, atoms, cv = collapsed()
    dirs, S = directions()
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, KE, X.G, X.SIGMA0, np.arange(3), 3, np.inf, weights=X.Q3, regressor=reg, compensation=comp)
    o = MN.second_order_transient(cl, cv, X.fprime, B_OU, 1.0, int(round(X.T_D / DT_TR)) * DT_TR, dt=DT_TR)
    Sd = MN.separation_drift_series(cv, o['t'], o['Ee'], o['Eem2'] if v == 'measured' else np.zeros_like(o['Ee']))
    Cd = MN.separation_cross_term(cv, o['t'], o['Eeth'])
    trec = np.arange(0.0, X.T_D + 1e-9, REC_D)
    idx = np.rint(trec / DT_TR).astype(int)
    fit = trec >= 50.0
    slope = lambda y: float(np.polyfit(trec[fit], y[fit], 1)[0])
    cum = lambda f: np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * DT_TR)])
    out = {}
    for k, d in FAMILIES[v]:
        th = dirs[(v, k, d)]
        fS = np.einsum('i,tij,j->t', th, Sd[:, k], th); fC = np.einsum('i,tij,j->t', th, Cd[:, k], th)
        V = MN.separation_variance(cl, cv, X.fprime, B_OU, 1.0, k, th, trec)
        out[f'{k}{d}'] = dict(k=k, d=d, mean=slope(cum(fS)[idx]), cross=slope(cum(fC)[idx]), var=slope(V), lift=float(th @ S[v][k] @ th),
                              window=[float(np.mean(fS[(o['t'] >= a_) & (o['t'] < a_ + 50.0)])) for a_ in np.arange(0.0, X.T_D, 50.0)])
    win = [dict(a=float(a_), b=float(a_ + 50.0), lam=[float(np.linalg.eigvalsh((Sd + Cd)[(o['t'] >= a_) & (o['t'] < a_ + 50.0)].mean(axis=0)[k])[-1]) for k in range(3)],
                lam_nocross=[float(np.linalg.eigvalsh(Sd[(o['t'] >= a_) & (o['t'] < a_ + 50.0)].mean(axis=0)[k])[-1]) for k in range(3)],
                rmsE=float(np.sqrt(np.mean(o['Ee'][(o['t'] >= a_) & (o['t'] < a_ + 50.0)] ** 2)))) for a_ in np.arange(0.0, X.T_D, 50.0)]
    print(f"[{time.time() - t0:.0f}s] D run predictions {NAMES[v]}", flush=True)
    return dict(v=v, dirs=out, win=win)


def part_D_pred(pool):
    """Part D's predictions: the lift and the scans, the stationary points (DR, DR*; adaptive and fixed prior), the long transients (DR, DR*)
    and the run predictions (MR, DR, DR*), as pool jobs."""
    a_lift = pool.apply_async(pred_lift)
    a_stat = pool.map_async(pred_stationary, [(v, eta) for v in ('reference', 'reference_both') for eta in (np.inf, 0.0)], chunksize=1)
    a_long = pool.map_async(pred_long, ['reference', 'reference_both'], chunksize=1)
    a_runs = pool.map_async(pred_runs, list(VARIANTS), chunksize=1)
    return a_lift, a_stat, a_long, a_runs


def collect_D_pred(handles):
    a_lift, a_stat, a_long, a_runs = handles
    out = a_lift.get()
    out['stat'] = {(r['v'], 'adaptive' if np.isinf(r['eta']) else 'fixed'): r for r in a_stat.get()}
    out['long'] = {r['v']: r for r in a_long.get()}
    out['runs'] = {r['v']: r for r in a_runs.get()}
    return out


def job_D(args):
    """The lift m = 15 with class k's particles separated by c_i·10⁻⁶ along G^{1/2}θ (θ: 'u' MR's top eigenvector of S₃ at the lift, Part D
    of `exp_measurement_noise.py`; 'top' the law's own top eigenvector of S_k), one noise path under the law v (std, seed, sign: ±1 the
    antithetic pair; the paths common to the laws and directions): the separation's growth rate, least squares on log‖separation‖ over t ≥ 50."""
    v, std, seed, sign = args[:4]
    k, d = (args[4], args[5]) if len(args) > 4 else (2, 'u')
    t0 = time.time()
    net, atoms, cv = collapsed()
    th = directions()[0][(v, k, d)]
    labels = np.repeat(np.arange(3), 5); m = 15
    reg, comp = LAW[v]
    cl = MN.NoisyRealLoop(net, X.f_vec, X.xr, X.xrd, X.PER, KE, X.G, X.SIGMA0, labels, 3, np.inf, regressor=reg, compensation=comp)
    z0 = cl.lift(atoms)
    for i, c in zip(range(5 * k, 5 * k + 5), (2.0, 1.0, 0.0, -1.0, -2.0)):
        z0[1 + 3 * i:4 + 3 * i] += 1e-6 * c * np.sqrt(X.G) * th

    def rec(t, z):
        Th = z[1:1 + m * 3].reshape(m, 3)[5 * k:5 * k + 5]
        return np.array([t, np.sqrt(np.sum(((Th - Th.mean(axis=0)) / np.sqrt(X.G)) ** 2))])
    out = np.array(MN.simulate_noisy(cl, z0, m, X.T_D, DT_C, B_OU, std ** 2, np.random.default_rng(5000 + seed), int(round(REC_D / DT_C)), rec,
                                     sign=sign))
    msk = out[:, 0] >= 50.0
    rate = float(np.polyfit(out[msk, 0], np.log(out[msk, 1]), 1)[0])
    print(f"[{time.time() - t0:.0f}s] D {NAMES[v]} std {std:g} seed {seed} sign {sign:+g} class {k + 1} {d}", flush=True)
    return dict(v=v, std=std, seed=seed, sign=sign, k=k, d=d, rate=rate, t=out[:, 0], sep=out[:, 1])


def report_D(D):
    """Part D's section of the report (Corollary D‴.33, as revised in round 23 of the review)."""
    pr = D['pred']
    s0 = X.SIGMA0 / KE
    thr = lambda lam, stiff=X.SIGMA0: np.sqrt(stiff / KE / lam) if lam > 0 else np.inf
    fs = lambda x: '—' if (x is None or np.isinf(x)) else f'{x:.3f}'
    mm = lambda e: f"[{e[0]:+.2e}, {e[1]:+.2e}]"
    L = ["", "## D. The separations under the desired regressor (Corollary D‴.33)", "",
         "### D1. At the lift", "",
         f"The separations' second-order drift S_k at the lift (real time, per unit Eν²; k_e = {KE:g}, B = {B_OU:g}, the noise's own quasi-static moments): its "
         "extreme eigenvalues per class and the noise std at which its largest eigenvalue reaches the separations' stiffness — σ₀/k_e under the adaptive prior, "
         "σ₀w_k/k_e under the cohesion leak (Proposition D‴.32, w_k = 1/3); '—' where every eigenvalue is negative.", "",
         "| law | class 1: [min, max] | class 2 | class 3 | threshold std, adaptive prior | threshold std, cohesion leak |", "|---|---|---|---|---|---|"]
    for v in VARIANTS:
        r = pr['lift'][v]; lam = max(r['hi'])
        L.append(f"| {NAMES[v]} | " + " | ".join(f"[{lo:+.5f}, {hi:+.5f}]" for lo, hi in zip(r['lo'], r['hi'])) + f" | {fs(thr(lam))} | "
                 f"{fs(thr(lam, X.SIGMA0 * pr['w'][0]))} |")
    L += ["", "### D2. Against k_e", "",
          "The drift's extreme eigenvalues over the classes (λ_min the smallest, λ_max the largest) against k_e, with the noise's bandwidth scaled with the gain "
          "(B = 4k_e) and fixed (B = 100); (ii)'s bound: the largest ratio ‖S_k‖/((1/2k_e)‖∂²_xf‖_∞⟨‖M_k‖⟩) over the classes (M_k = G^{1/2}∇²_ϑΦ̃G^{1/2}; "
          f"‖∂²_xf‖_∞ = {pr['fpp_inf']:.4f} on the orbit, ⟨‖M_k‖⟩ = " + ", ".join(f'{x:.4f}' for x in pr['Mnorm']) + "), and the threshold std per law (adaptive "
          f"prior). The bound gives the threshold a floor uniform in k_e, std ≥ (2σ₀/(‖∂²_xf‖_∞max_k⟨‖M_k‖⟩))^{{1/2}} = {pr['floor']:.4f}, for DR and DR*.", "",
          "| B | k_e | MR: λ_max | DR: λ_min, λ_max; k_e·λ_max | DR*: λ_min, λ_max; k_e·λ_max | bound ratio: DR; DR* | threshold std: MR; DR; DR* |",
          "|---|---|---|---|---|---|---|"]
    for Bk in ('4ke', 100.0):
        for ke in KES_D:
            sc = {v: pr['scan'][(ke, Bk, v)] for v in VARIANTS}
            L.append(f"| {'4k_e' if Bk == '4ke' else f'{Bk:g}'} | {ke:g} | {sc['measured']['lam']:+.4e} | {sc['reference']['lam_lo']:+.2e}, {sc['reference']['lam']:+.3e}; "
                     f"{ke * sc['reference']['lam']:+.4f} | {sc['reference_both']['lam_lo']:+.2e}, {sc['reference_both']['lam']:+.3e}; {ke * sc['reference_both']['lam']:+.4f} | "
                     f"{sc['reference']['ratio']:.3f}; {sc['reference_both']['ratio']:.3f} | "
                     + "; ".join(fs(np.sqrt(X.SIGMA0 / ke / sc[v]['lam']) if sc[v]['lam'] > 0 else np.inf) for v in VARIANTS) + " |")
    L += ["", "### D3. Other targets (DR at the lift)", "",
          f"Exactly represented targets with q = 1/K₀ (k_e = {KE:g}, B = {B_OU:g}): DR's drift per class (the centre entry, the (γ, γ) entry of S_k, and the extreme "
          "eigenvalues), DR*'s and MR's largest eigenvalues.", "",
          "| target | bumps (A, C, S) | DR: centre entries | DR: [min, max] per class | DR*: λ_max | MR: λ_max |", "|---|---|---|---|---|---|"]
    for name, t_ in pr['targets'].items():
        r = t_['reference']
        L.append(f"| {name} | " + "; ".join(f"({A_:g}, {C_:g}, {S_:g})" for A_, C_, S_ in t_['bumps']) + " | " + ", ".join(f"{c_:+.2e}" for c_ in r['centre'])
                 + " | " + "; ".join(mm((lo, hi)) for lo, hi in zip(r['lo'], r['hi'])) + f" | {max(t_['reference_both']['hi']):+.2e} | {max(t_['measured']['hi']):+.3e} |")
    L += ["", "### D4. At the stationary point", "",
          "The exact periodic second-order mean of the K₀-atom loop (`second_order_mean(periodic=True)`, the forcing split into the noise part — the noise's "
          "own second moments — and the parameters' part), per unit Eν²: the mean error's rms (the lift's; the round-22 averaged noise-only balance (I − P)g/k_e; "
          "exact: total, noise part, parameters' part; the tangent part ‖PE[e]‖ and the orthogonal part), and the drift's [min, max] per class from the noise "
          "part of E[e], from its parameters' part, from the cross term, and in total.", "",
          "| law | prior | rms E[e]: lift; averaged (I − P)g/k_e | exact: total; noise; parameters' | tangent; orthogonal | drift, noise part of E[e] | "
          "parameters' part of E[e] | cross term | total | threshold std (total) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for (v, prior), r in sorted(pr['stat'].items(), key=lambda kv: (kv[0][0] != 'reference', kv[0][1])):
        lam = max(e[1] for e in r['eig_tot'])
        L.append(f"| {NAMES[v]} | {prior} | {r['rms_lift']:.3e}; {r['rms_avg']:.3e} | {r['rms']:.3e}; {r['rms_noise']:.3e}; {r['rms_param']:.3e} | "
                 f"{r['tang']:.3e}; {r['orth']:.3e} | " + " | ".join("; ".join(mm(e) for e in r[key]) for key in ('eig_noise', 'eig_param', 'eig_cross', 'eig_tot'))
                 + f" | {fs(thr(lam))} |")
    L += ["", "### D5. Along the second-order transient from the lift", ""]
    for v in ('reference', 'reference_both'):
        if v not in pr['long']:
            continue
        g_ = pr['long'][v]
        tpos = lambda xs: ", ".join('—' if x is None else f'{x:.0f}' for x in xs)
        L += [f"{NAMES[v]}: the exact second-order transient (`long_transient`, {g_['n_periods']} periods ≈ {g_['n_periods'] * X.PER:.3g} time units, {g_['N']} grid "
              "points per period), per unit Eν²: the period-averaged drift's extreme eigenvalues per class (all second-order terms: the mean error's noise and "
              "parameters' parts, the cross term) and the mean error's rms (total; noise; parameters'). The largest eigenvalue first positive at t = "
              f"{tpos(g_['first_pos'])} (classes 1, 2, 3; with the noise part of E[e] alone at t = {tpos(g_['first_pos_noise'])}), its maximum "
              + ", ".join(f"{x:+.3e} (t ≈ {t:.3g})" for x, t in zip(g_['max_tot'], g_['t_max'])) + " (the noise part alone "
              + ", ".join(f"{x:+.3e} (t ≈ {t:.3g})" for x, t in zip(g_['max_noise'], g_['t_max_noise'])) + "); at the stationary point (the period maps' limit) "
              + "; ".join(mm(e) for e in g_['stat']['eig_tot']) + f", rms E[e] {g_['stat']['rms']:.3e}.", "",
              "| t | class 1: [min, max] | class 2 | class 3 | rms E[e]: total; noise; parameters' |", "|---|---|---|---|---|"]
        tsel = [0.0, 50.0, 300.0, 1257.0, 2700.0, 4100.0, 6283.0, 1.26e4, 3.2e4, 1e5, 2.5e5, 5.1e5]
        recs = g_['rec']
        for ts_ in tsel:
            r = min(recs, key=lambda r_: abs(r_['t'] - ts_))
            L.append(f"| {r['t']:.4g} | " + " | ".join(mm(e) for e in r['eig_tot']) + f" | {r['rms']:.3e}; {r['rms_noise']:.3e}; {r['rms_param']:.3e} |")
        L.append("")
    if 'runs' in pr and 'reference' in pr['runs']:
        w = pr['runs']['reference']['win']
        L += ["DR over the real loop's 300 time units (`second_order_transient`, dt = " f"{DT_TR:g}): per 50-unit window the mean error's rms and the largest eigenvalue "
              "of the window-averaged drift per class, with the cross term and without it.", "",
              "| window | rms E[e] | λ_max(S₁), λ_max(S₂), λ_max(S₃) | without the cross term |", "|---|---|---|---|"]
        for w_ in w:
            L.append(f"| [{w_['a']:g}, {w_['b']:g}) | {w_['rmsE']:.4f} | " + ", ".join(f"{x:+.3e}" for x in w_['lam']) + " | "
                     + ", ".join(f"{x:+.3e}" for x in w_['lam_nocross']) + " |")
        if 'reference_both' in pr['runs']:
            w0 = pr['runs']['reference_both']['win'][0]
            l1, l0 = max(w0['lam']), max(w0['lam_nocross'])
            L += ["", f"DR* over [0, 50): the largest eigenvalue {l1:+.4e} with the cross term ({l0:+.4e} without it), the threshold std {fs(thr(l1))} "
                  f"({fs(thr(l0))}; at the lift's quasi-static moments {fs(thr(max(pr['lift']['reference_both']['hi'])))})."]
    L += ["", "### D6. The real loop", "",
          f"The lift m = 15 with class k separated by 10⁻⁶·(2, 1, 0, −1, −2) along a start direction θ — u (MR's top eigenvector of S₃ at the lift) or the law's own "
          f"top eigenvector of S_k ('top') —, {N_D} antithetic pairs per family on paths common to the laws, {X.T_D:g} time units, the separation's rate (least squares "
          f"over t ≥ 50) against −σ₀/k_e = {-s0:.4e} plus Eν² times: the lift's θᵀS_kθ, and the transient's prediction — the least-squares slope of "
          "∫₀ᵗθᵀ(S_k + C_k)θ + E‖P_⊥v₁(t)‖² (the mean error's term along the transient, the cross term, the variance term). The noise-induced part is the rate "
          "plus σ₀/k_e; its ratio to the prediction is given with its sampling error.", "",
          "| law | std | class, θ | measured rate (pair means ± s.e.) | lift prediction | transient prediction (mean-error, cross, variance terms per unit Eν²) | "
          "noise part: measured / transient (± s.e.) | runs growing | runs slower than σ₀/k_e |", "|---|---|---|---|---|---|---|---|---|"]
    fam = sorted({(r['v'], r['std'], r.get('k', 2), r.get('d', 'u')) for r in D['runs']}, key=lambda f_: (VARIANTS.index(f_[0]), f_[3] != 'u', f_[2], f_[1]))
    for v, std, k, d in fam:
        rs = [r for r in D['runs'] if (r['v'], r['std'], r.get('k', 2), r.get('d', 'u')) == (v, std, k, d)]
        seeds = sorted(set(r['seed'] for r in rs))
        pm = np.array([np.mean([r['rate'] for r in rs if r['seed'] == sd]) for sd in seeds])
        se = pm.std(ddof=1) / np.sqrt(len(pm)) if len(pm) > 1 else np.nan
        q = pr['runs'].get(v, {}).get('dirs', {}).get(f'{k}{d}') if 'runs' in pr else None
        if q is None:
            continue
        tr = q['mean'] + q['cross'] + q['var']
        meas = pm.mean() + s0
        L.append(f"| {NAMES[v]} | {std:g} | {k + 1}, {d} | {pm.mean():+.5e} ± {se:.1e} | {-s0 + std ** 2 * q['lift']:+.5e} | {-s0 + std ** 2 * tr:+.5e} "
                 f"({q['mean']:+.3e}, {q['cross']:+.2e}, {q['var']:+.2e}) | {meas / (std ** 2 * tr):.3f} ± {se / abs(std ** 2 * tr):.3f} | "
                 f"{sum(r['rate'] > 0 for r in rs)} of {len(rs)} | {sum(r['rate'] > -s0 for r in rs)} of {len(rs)} |")
    return L


def report(A, Bp, Cres, Cpred, t0=None, D=None):
    es = lambda e: '∞' if np.isinf(e) else ('0 (fixed)' if e == 0 else f'{e:.5f}')
    L = ["# Proposition D‴.31: the adaptive prior is a washout leak; the desired regressor (`exp_eiv_remedy.py`)", "",
         f"The instance of Proposition D‴.30 (σ₀ = {X.SIGMA0:.5f}, x_m = 2 sin t, the lift of the exact representation of q = 1/3 with the prior at its atoms) "
         "under the measurement y = x + ν, with an Ornstein–Uhlenbeck noise of rate B. MR: the measured regressor and compensation (Proposition D‴.30); DR: the "
         "regressor at the reference, the compensation f̂(y); DR*: both at the reference (the desired compensation adaptation law).", "",
         "## A. The second-order force and its stationary responses (no dynamics)", "",
         "Per unit (Eν²)², function level, the averaged stationary response: the bias of the mean parameters under the fixed prior and under the adaptive prior "
         "(η = ∞, the noise part's for every η > 0; = FᵀH⁻¹F, the squared norm of the smallest target perturbation exerting F) for the force of the noise's own second moments (quasi-static, "
         "`eiv_force`) and for the full second-order force (the parameters' fluctuation included, `second_order_mean`), the bias of the mean function under the "
         "adaptive prior with the full force (‖Ĵȳ + ½D²_Θf̂[Σ_ΘΘ]‖²_ρ), the least-squares bound ‖g‖²_ρ (g = k_eE[e] of the noise part; DR, DR*), the rescaled "
         "projected mean tracking error at the adaptive prior's stationary point, k_e²‖P E[e]‖²_ρ = F₁ᵀH⁻¹F₁ (F₁ the regressor's errors in variables; zero by "
         "construction for DR and DR*), and the mean tracking error itself ‖P E[e]‖²_ρ = F₁ᵀH⁻¹F₁/k_e². The last column is the noise variance at which the adaptive "
         "prior's full bias equals its floor of the parameters (N₀/2k_e)tr Q, N₀ = 2Eν²/B (std above 0.1 is outside the small-noise regime).", "",
         "| k_e | B | law | fixed: noise part | fixed: full | adaptive: noise part | adaptive: full | adaptive: mean function | ‖g‖²_ρ | k_e²‖P E[e]‖²_ρ | ‖P E[e]‖²_ρ | Eν² at bias = floor (std) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in A['rows']:
        g2 = f"{r['g2']:.4e}" if r['v'] != 'measured' else '—'
        fa, ff = r['full_adapt'], r['full_fixed']
        L.append(f"| {r['ke']:g} | {r['B']:g} | {NAMES[r['v']]} | {r['fixed']:.4e} | {ff['full']:.4e} | {r['adapt']:.4e} | {fa['full']:.4e} | "
                 f"{fa['function']:.4e} | {g2} | {r['merr']:.3e} | {r['merr'] / r['ke'] ** 2:.3e} | {r['cross']:.2e} ({np.sqrt(r['cross']):.4f}) |")
    rel = lambda r: r['full_adapt']['full'] / r['adapt'] - 1
    L += ["", "The parameters' fluctuation's share of the adaptive prior's bias (full over noise part − 1): " + "; ".join(
        f"{NAMES[v]} " + ", ".join(f"{rel(r):+.2%}" for r in A['rows'] if r['v'] == v) for v in VARIANTS) + " (the rows' order); of the mean function's: "
        + "; ".join(f"{NAMES[v]} " + ", ".join(f"{r['full_adapt']['function'] / r['adapt'] - 1:+.2%}" for r in A['rows'] if r['v'] == v) for v in VARIANTS) + "."]
    if 'param' in A['rows'][0]['full_adapt']:
        L += ["", "The parameters' part's own function-level response under the adaptive prior (and its learning-law part's), per unit (Eν²)²: " + "; ".join(
            f"{NAMES[v]} " + ", ".join(f"{r['full_adapt']['param']:.2e} ({r['full_adapt']['law']:.2e})" for r in A['rows'] if r['v'] == v) for v in VARIANTS)
            + " (the rows' order) — set by the parameters' covariance, falling with k_e."]
    mr = {(r['ke'], r['B']): r['full_adapt']['full'] for r in A['rows'] if r['v'] == 'measured'}
    dr = {(r['ke'], r['B']): r['full_adapt']['full'] for r in A['rows'] if r['v'] == 'reference'}
    mr25, mr200 = mr[(25.0, 1000.0)], mr[(200.0, 1000.0)]
    law = lambda ke, B: (ke * B / (ke + B)) ** 2
    L += ["", "The adaptive prior's full bias, MR over DR: " + ", ".join(f"{mr[k] / dr[k]:.3g} (k_e = {k[0]:g}, B = {k[1]:g})" for k in sorted(mr)) + ". "
          f"MR's grows from k_e = 25 to 200 at B = 1000 by {mr200 / mr25:.1f} — the law (k_eB/(k_e + B))² gives {law(200.0, 1000.0) / law(25.0, 1000.0):.1f}, k_e² 64; "
          "the noise parts of DR's and DR*'s stay below ¼‖∂²_xf‖²_ρ = "
          f"{0.25 * np.mean(X.fsecond(2 * np.sin(np.arange(4096) * X.PER / 4096)) ** 2):.4e} for every k_e and B, and their full biases at most "
          f"{max(r['full_adapt']['full'] for r in A['rows'] if r['v'] != 'measured'):.4f} here."]
    L += ["", "## B. The stationary identity RH + B = I along the prior family dz/ds = −η(z − T̂) − κ(z − z_a)", "",
          f"η = σ₀; R and B from the averaged collective block at zero frequency (k_e = {KE:g}, B = {B_OU:g}). h₁ = {Bp['h'][0]:.4e}, h_max = {Bp['h'][-1]:.4e}. "
          "The errors-in-variables biases are function-level per unit (Eν²)²; the anchor's biases are function-level for a unit misplacement along the weakest and "
          "the strongest mode.", "",
          "| prior | S₀ (zero-frequency leak) | weakest mode: response r₁ | bias factor b₁ | h₁r₁ + b₁ | EIV bias, MR | EIV bias, DR | anchor bias, weakest | anchor bias, strongest | max|RH + B − I| |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in Bp['rows']:
        name = 'fixed (η = 0)' if r['eta'] == 0 else (f"adaptive, κ = 0" if r['kappa'] == 0 else f"κ = {r['kappa'] / r['eta']:g}η")
        L.append(f"| {name} | {r['S0']:.4e} | {r['r1']:.4e} | {r['b1']:.4f} | {Bp['h'][0] * r['r1'] + r['b1']:.6f} | {r['eiv']['measured']:.4e} | "
                 f"{r['eiv']['reference']:.4e} | {r['anchor'][0]:.3e} | {r['anchor'][1]:.3e} | {r.get('resid_RH', r['resid']):.1e} |")
    ratios = [r['eiv']['measured'] / r['eiv']['reference'] for r in Bp['rows']]
    L += ["", "Along the family the zero-frequency leak S₀ = σ₀κ/(η + κ) trades the anchor's bias for the response to persistent forces, mode by mode, at "
          f"h_jr_j + b_j = 1; the trade binds on the modes with h_j < σ₀ ({Bp.get('n_weak', '—')} of {len(Bp['h'])} here), where 1/h_j exceeds the fixed "
          "prior's 1/σ₀. The desired regressor changes the force instead: MR's bias over DR's (noise parts) along the family "
          + ", ".join(f"{x:.3g}" for x in ratios) + "."]
    if 'general' in Bp:
        if 'vs_def' in Bp['general'][0]:
            L += ["", "Priors whose zero-frequency leak S₀ does not commute with H: R and B from the dynamic block of (y, z_f) at zero frequency "
                  "(`anchored_block`; the block's dimension), against (H + S₀)⁻¹ and (H + S₀)⁻¹S₀ with S₀ = σ₀[E_xE_xᵀ + E_f(ηI + K)⁻¹KE_fᵀ] by "
                  "elimination of z_f (the largest relative difference of R, the largest difference of B): RH + B = I holds, HR + B = I does not; "
                  "per mode h_jr_jj + b_jj = 1 holds whatever S₀.", "",
                  "| prior | block dimension | block against the definitions | max|HS₀ − S₀H| | max|RH + B − I| | max|HR + B − I| | max_j|h_jr_jj + b_jj − 1| | EIV bias, MR | EIV bias, DR |",
                  "|---|---|---|---|---|---|---|---|---|"]
            for r in Bp['general']:
                L.append(f"| {r['name']} | {r['dim']} | {r['vs_def']:.1e} | {r['comm']:.2e} | {r['RH']:.1e} | {r['HR']:.3f} | {r['per_mode']:.1e} | "
                         f"{r['eiv']['measured']:.4e} | {r['eiv']['reference']:.4e} |")
        else:
            L += ["", "Zero-frequency leaks S₀ that do not commute with H (R = (H + S₀)⁻¹, B = (H + S₀)⁻¹S₀): RH + B = I holds, HR + B = I does not; per mode "
                  "h_jr_jj + b_jj = 1 holds whatever S₀.", "",
                  "| S₀ | max|HS₀ − S₀H| | max|RH + B − I| | max|HR + B − I| | max_j|h_jr_jj + b_jj − 1| | EIV bias, MR | EIV bias, DR |", "|---|---|---|---|---|---|---|"]
            for r in Bp['general']:
                L.append(f"| {r['name']} | {r['comm']:.2e} | {r['RH']:.1e} | {r['HR']:.3f} | {r['per_mode']:.1e} | {r['eiv']['measured']:.4e} | {r['eiv']['reference']:.4e} |")
    if Cres:
        L += ["", f"## C. The real loop: MR, DR and DR* driven by the same noise paths (k_e = {KE:g}, B = {B_OU:g}, std {STD_C:g}, η = ∞, RK4 dt = {DT_C:g})", "",
              f"E‖f̂ − f‖²_ρ at the records (mean ± standard error over {N_RUNS} runs) against the prediction (covariance part + mean part); MR and DR share the "
              "covariance part (their linearisations coincide).", "",
              "| t | MR: Monte Carlo | MR: predicted | DR: Monte Carlo | DR: predicted | DR*: Monte Carlo | DR*: predicted | covariance part (MR, DR) | mean part: MR | DR | DR* |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
        E = {v: np.array([r['err'] for r in sorted([r for r in Cres if r['v'] == v], key=lambda r: r['seed'])]) for v in VARIANTS}
        Dv = {v: np.array([r['dev'] for r in sorted([r for r in Cres if r['v'] == v], key=lambda r: r['seed'])]) for v in VARIANTS}
        for k_, t in enumerate(Cpred['measured']['t']):
            if k_ == 0:
                continue
            cells = []
            for v in VARIANTS:
                p = Cpred[v]
                cells += [f"{E[v][:, k_].mean():.4e} ± {E[v][:, k_].std(ddof=1) / np.sqrt(len(E[v])):.1e}", f"{p['cov'][k_] + p['mean'][k_]:.4e}"]
            L.append(f"| {t:g} | " + " | ".join(cells) + f" | {Cpred['measured']['cov'][k_]:.4e} | {Cpred['measured']['mean'][k_]:.4e} | "
                     f"{Cpred['reference']['mean'][k_]:.4e} | {Cpred['reference_both']['mean'][k_]:.4e} |")
        diff = Dv['measured'] - Dv['reference']                           # path by path: the first-order responses cancel
        pd = Cpred['measured']['dev'] - Cpred['reference']['dev']
        L += ["", "The run mean of f̂_MR − f̂_DR on the orbit's grid (common random numbers) against Ĵ(ȳ_MR − ȳ_DR): its mean square (with the run-to-run "
              "standard error of the mean function, in the same units), the prediction's, and the relative L²_ρ error of the prediction.", "",
              "| t | ⟨(run mean)²⟩ | predicted | relative error ‖run mean − Ĵ(ȳ_MR − ȳ_DR)‖/‖Ĵ(ȳ_MR − ȳ_DR)‖ | sampling error of the run mean (rms) | relative error over the relative sampling error |",
              "|---|---|---|---|---|---|"]
        for k_, t in enumerate(Cpred['measured']['t']):
            if k_ == 0:
                continue
            m_ = diff[:, k_].mean(axis=0); se = np.sqrt(np.mean(diff[:, k_].var(axis=0, ddof=1) / len(diff)))
            re_ = np.sqrt(np.mean((m_ - pd[k_]) ** 2) / np.mean(pd[k_] ** 2))
            L.append(f"| {t:g} | {np.mean(m_ ** 2):.4e} | {np.mean(pd[k_] ** 2):.4e} | {re_:.4f} | {se:.2e} | {re_ / (se / np.sqrt(np.mean(pd[k_] ** 2))):.2f} |")
    if D:
        L += report_D(D)
    L += ["", "Numbers are produced by `experiments/exp_eiv_remedy.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'eiv_remedy.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'eiv_remedy_raw.npz')
    if '--from-raw' in sys.argv:
        with np.load(raw, allow_pickle=True) as z:
            report(z['A'].item(), z['B'].item(), list(z['C']), z['P'].item(), t0, D=z['D'].item() if 'D' in z.files else None)
        return
    if '--only-Dpred' in sys.argv or '--only-D2' in sys.argv:                # Part D's predictions recomputed (and with --only-D2 the runs
        with np.load(raw, allow_pickle=True) as z:                            # along the laws' own top eigenvectors added), its runs kept
            save = {k: z[k] for k in z.files}
        Dd = save['D'].item()
        with mp.Pool(N_WORKERS) as pool:
            h = part_D_pred(pool)
            if '--only-D2' in sys.argv:
                extra = pool.map(job_D, [(v, std, seed, sign, k, d) for v, std, k, d in EXTRA_D for seed in range(N_D) for sign in (1.0, -1.0)],
                                 chunksize=1)
                keyf = lambda r: (r['v'], r['std'], r['seed'], r['sign'], r.get('k', 2), r.get('d', 'u'))
                new_keys = {keyf(r) for r in extra}
                Dd['runs'] = [r for r in Dd['runs'] if keyf(r) not in new_keys] + extra
            Dd['pred'] = collect_D_pred(h)
        save['D'] = np.array(Dd, dtype=object)
        X.save_raw(raw, **save)
        report(save['A'].item(), save['B'].item(), list(save['C']), save['P'].item(), t0, D=Dd)
        return
    if '--only-D' in sys.argv:                                               # Part D (Corollary D‴.33) computed, merged with the stored parts
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        with mp.Pool(N_WORKERS) as pool:
            h = part_D_pred(pool)
            runs = pool.map(job_D, [(v, std, seed, sign) for v in ('measured', 'reference') for std in STDS_D for seed in range(N_D)
                                    for sign in (1.0, -1.0)]
                            + [(v, std, seed, sign, k, d) for v, std, k, d in EXTRA_D for seed in range(N_D) for sign in (1.0, -1.0)], chunksize=1)
            save['D'] = np.array(dict(pred=collect_D_pred(h), runs=runs), dtype=object)
        X.save_raw(raw, **save)
        report(save['A'].item(), save['B'].item(), list(save['C']), save['P'].item(), t0, D=save['D'].item())
        return
    if '--only-B' in sys.argv:                                               # Part B recomputed, merged with the stored A, C and P
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        save['B'] = np.array(part_B(), dtype=object)
        X.save_raw(raw, **save)
        report(save['A'].item(), save['B'].item(), list(save['C']), save['P'].item(), t0, D=save['D'].item() if 'D' in save else None)
        return
    if '--only-AB' in sys.argv:                                              # Parts A and B recomputed, merged with the stored C and P
        with np.load(raw, allow_pickle=True) as z:
            save = {k: z[k] for k in z.files}
        with mp.Pool(N_WORKERS) as pool:
            save['A'] = np.array(part_A(pool), dtype=object)
        save['B'] = np.array(part_B(), dtype=object)
        X.save_raw(raw, **save)
        report(save['A'].item(), save['B'].item(), list(save['C']), save['P'].item(), t0, D=save['D'].item() if 'D' in save else None)
        return
    with mp.Pool(N_WORKERS) as pool:
        A, Bp = part_A(pool), part_B()
        X.save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bp, dtype=object))
        Cres = pool.map(job_C, [(v, seed) for seed in range(N_RUNS) for v in VARIANTS], chunksize=1)
        h = part_D_pred(pool)
        Druns = pool.map(job_D, [(v, std, seed, sign) for v in ('measured', 'reference') for std in STDS_D for seed in range(N_D)
                                 for sign in (1.0, -1.0)]
                         + [(v, std, seed, sign, k, d) for v, std, k, d in EXTRA_D for seed in range(N_D) for sign in (1.0, -1.0)], chunksize=1)
        Dd = dict(pred=collect_D_pred(h), runs=Druns)
    Cpred = predictions_C()
    X.save_raw(raw, A=np.array(A, dtype=object), B=np.array(Bp, dtype=object), C=np.array(Cres, dtype=object), P=np.array(Cpred, dtype=object),
               D=np.array(Dd, dtype=object))
    report(A, Bp, Cres, Cpred, t0, D=Dd)


if __name__ == '__main__':
    main()
