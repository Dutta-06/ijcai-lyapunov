"""
Theorem D‴.21: a measure-level basin for the barycentric OT leak — the Polyak–Łojasiewicz inequality along the runs
================================================================================================================

The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit, as in `exp_ot_leak.py`: commensurate sizes m = 15, 60, 240 and
M = 510 from the warm clouds of Proposition D‴.12 (the jittered mixture of width 0.02 around the exact representation, one cloud per seed,
the m-system its first m rows), four seeds, 600 slow units, midpoint rule with dt = 0.02 and the leak target recomputed every 10 steps, under
the OT leak and — for the contrast — the nearest-atom leak. Two priors: the matching prior (the exact representation; r* = 0, μ* = ν) and the
prior with centres shifted by +0.15 (μ* = Θ*(q) by Newton, r* ≠ 0).

Records every 5 slow units, for every system, against μ* (src/core/ot_basin.py): the excess energy E = 𝒥(μ_Θ) − 𝒥* (𝒥 := L̄ + (σ₀/2)W²(·, ν),
the OT functional, evaluated along both leaks' trajectories), the class-respecting distance D_a, W²(μ_Θ, μ*), ‖Δf̂‖², the exact effective
margins κ_E and κ_Q of the identities (E) and (P) and their residuals, the observed PL ratio |∇J_a|²_G/(mE) (= −d log E/ds along the OT flow)
and its lower bound Q²/(DE), and the plan changes; every 25 slow units the plan-segment curvatures κ*, κ*_−, κ^Θ (9 points per segment),
the margin α = σ₀ − (κ* + κ^Θ)/2 and the theorem's rate λ = α·min{2, 2α/(σ₀ + κ*_−)}.

The report adds: the a-posteriori bound E(0)exp(−∫λ₊) against E(s) (λ interpolated linearly between the segment records, trapezoid rule);
the local PL constant 2min{λ_sep, λ_coll} at μ*; the a-priori basin E₀ = κ²/(2(L₂^G)²) with κ := σ₀ − sup λ_max R^G[μ*] over the whole
parameter space (a grid over the saturation ranges with Nelder–Mead polish) and over the box of the end configurations, the priors and the
atoms (random points with L-BFGS-B polish; the trajectories happen to lie inside it, but it is not forward-invariant), with the largest
fraction of particles outside their cells that {E ≤ E₀} can hold; for the nearest-atom runs, the end energy against the energy of the
replicated configuration of the sampled weights, the admissibility of those weights (the separation margins and the whitened weighted
K-atom Hessian), the end state's class spreads and E at s = 300, 450, 600; and the measure-level basin (particles moved into another cell,
the segment suprema on 129 points per segment), with one particle moved: the OT run with ∫λ₊ evaluated every step, and the nearest-atom
run continued to 3000 slow units against the spurious minimum of its counts.

Output: results/chaos/ot_basin.md and ot_basin_raw.npz (OTB_SEEDS, OTB_S_END, OTB_WORKERS to change the pool; `--from-raw` regenerates the md;
`--only consts,adm,moved` recomputes those parts from the stored runs and merges them into the raw file).
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
from src.core import chaos_learning as CL
from src.core import prior_design as PD
from src.core import ot_leak as OT
from src.core import ot_basin as OB

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
M_REF, MS = 510, (15, 60, 240)
S_END = float(os.environ.get('OTB_S_END', 600.0))
DT, PLAN_EVERY, REC, SEG = 0.02, 10, 250, 1250                          # records every 5 slow units, segment curvatures every 25
JITTER = 0.02
N_SEEDS = int(os.environ.get('OTB_SEEDS', 4))
N_WORKERS = int(os.environ.get('OTB_WORKERS', 2))
Q = np.ones(3) / 3
SHIFTS = {'matching': 0.0, 'shifted': 0.15}
TIMES = (0.0, 25.0, 100.0, 300.0, 600.0)


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def representation(net, shift=0.0):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C + shift] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def setup(prior_name):
    net = make_net()
    prior = representation(net, SHIFTS[prior_name])
    loop = OT.OTPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, prior, n_quad=N_Q)
    atoms, res, marg = loop.stationary_atoms(prior)
    return net, loop, atoms, res, marg


def record(B, Th, prev_lab):
    t = B.terms(Th)
    W2s = OT.W2_empirical(Th, B.atoms, G)
    sw = int(np.sum(t['lab'] != prev_lab)) if prev_lab is not None else 0
    rel = lambda a, b: abs(a - b) / max(abs(a), 1e-300)
    return dict(E=t['E'], D=t['D'], W2s=W2s, ndf2=t['ndf2'], kE=t['kappa_E'], kQ=t['kappa_Q'], ratio=t['ratio'], lam_eff=t['lam_eff'],
                idE=rel(t['E'], t['E_ident']), idQ=rel(t['Q'], t['Q_ident']), switches=sw), t['lab'], t['r']


def job(args):
    seed, prior_name, leak = args
    t0 = time.time()
    net, loop, atoms, res, marg = setup(prior_name)
    B = OB.BasinDiagnostics(loop, atoms)
    exact = representation(net, 0.0)
    rng = np.random.default_rng(seed)
    clouds = CL.mixture_start(exact, M_REF, rng, JITTER, G)
    sys_ = {m: clouds[:m].copy() for m in list(MS) + [M_REF]}
    n_steps = int(round(S_END / DT))
    keys = ('E', 'D', 'W2s', 'ndf2', 'kE', 'kQ', 'ratio', 'lam_eff', 'idE', 'idQ', 'switches')
    rec = {m: dict(s=[], **{k: [] for k in keys}) for m in sys_}
    seg = {m: dict(s=[], kstar=[], kstar_minus=[], ktheta=[], alpha=[], lam=[]) for m in sys_}
    prev = {m: None for m in sys_}
    T = {}
    for k in range(n_steps + 1):
        if k % REC == 0:
            for m, Th in sys_.items():
                d, lab, r = record(B, Th, prev[m]); prev[m] = lab
                rec[m]['s'].append(k * DT)
                for kk in keys:
                    rec[m][kk].append(d[kk])
                if k % SEG == 0:
                    sc = B.segment_curvatures(Th, lab, n_t=9, r=r)
                    al, lam = B.theorem_rate(**sc)
                    seg[m]['s'].append(k * DT); seg[m]['kstar'].append(sc['kappa_star']); seg[m]['kstar_minus'].append(sc['kappa_star_minus'])
                    seg[m]['ktheta'].append(sc['kappa_theta']); seg[m]['alpha'].append(al); seg[m]['lam'].append(lam)
        if k == n_steps:
            break
        if k % PLAN_EVERY == 0:
            for m, Th in sys_.items():
                T[m] = loop.leak_target(Th) if leak == 'ot' else loop.nearest.leak_target(Th)
        for m, Th in list(sys_.items()):
            k1 = loop.drift(Th, T[m])
            sys_[m] = Th + DT * loop.drift(Th + 0.5 * DT * k1, T[m])
    out = dict(seed=seed, prior=prior_name, leak=leak, atoms=atoms, newton_res=res, diag_margin=marg, Jstar=B.Jstar, stationarity=B.stationarity,
               rec={m: {k: np.array(v) for k, v in r.items()} for m, r in rec.items()},
               seg={m: {k: np.array(v) for k, v in r.items()} for m, r in seg.items()},
               final={m: Th.astype(np.float32) for m, Th in sys_.items()})
    # the nearest-atom contrast: the replicated configuration of the sampled weights (a critical point of the nearest-atom energy)
    if leak == 'nearest':
        spur = {}
        for m, Th in sys_.items():
            lab = np.argmin(loop.cost_matrix(Th), axis=1)
            counts = np.bincount(lab, minlength=3)
            if np.all(counts > 0):
                w = counts / m
                Tw, resw = PD.weighted_stationary_atoms_prior(loop.nearest, w, atoms, w0=Q, n_cont=4)
                lift = np.repeat(Tw, counts, axis=0)
                d0 = float(np.abs(loop.nearest.drift(lift)).max())
                dOT = float(np.abs(loop.drift(lift)).max())
                spur[m] = dict(counts=counts, newton_res=resw, E=B.terms(lift)['E'], drift_nearest=d0, drift_ot=dOT)
        out['spurious'] = spur
    print(f"[{time.time() - t0:.0f}s] seed {seed} {prior_name} {leak}", flush=True)
    return out


def k_atom_hessian(nloop, w, Tw, h=1e-6):
    """The whitened weighted Hessian of the K-atom energy 𝒢_w (nearest-atom leak, each atom in its own cell) at Tw: eigenvalues ascending."""
    K = Tw.shape[0]; d = 3 * K
    g = lambda v: PD.weighted_energy_grad_prior(nloop, w, v.reshape(K, 3)).ravel()
    J = np.zeros((d, d))
    for j in range(d):
        e = np.zeros(d); e[j] = h
        J[:, j] = (g(Tw.ravel() + e) - g(Tw.ravel() - e)) / (2 * h)
    sw = np.repeat(np.sqrt(w), 3) * np.tile(SQG, K); sv = np.tile(SQG, K) / np.repeat(np.sqrt(w), 3)
    H = sw[:, None] * J * sv[None, :]
    return np.linalg.eigvalsh(0.5 * (H + H.T))


def weight_admissibility(loop, atoms, counts):
    """At the replicated configuration of the counts N (weights w = N/m) under the nearest-atom leak: Θ*(w) by continuation from q, its
    Newton residual, each class's separation margin σ₀ − λ_max R^G_k(w), the smallest eigenvalue of the K-atom Hessian, and whether every atom
    lies in its own cell. Admissible iff all margins and the Hessian are positive."""
    w = np.asarray(counts, float) / np.sum(counts)
    Tw, res = PD.weighted_stationary_atoms_prior(loop.nearest, w, atoms, w0=Q, n_cont=4)
    a, c, s_, diff, phi = loop.features(Tw)
    sep = SIGMA0 - loop.kappa_R(loop.F - phi @ (w * a), Tw)
    hmin = float(k_atom_hessian(loop.nearest, w, Tw)[0])
    own = bool(np.all(loop.nearest.assignment(Tw) == np.arange(Tw.shape[0])))
    return dict(Tw=Tw, res=res, sep=sep, hess_min=hmin, own_cells=own, admissible=bool(np.all(sep > 0) and hmin > 0 and own))


def admissibility(runs):
    """For every nearest-atom run and size: the admissibility of the sampled counts, the end state's class spreads (rms G⁻¹-distance to the class
    mean, classes by nearest prior atom) and its largest G⁻¹-distance from the class atom of Θ*(w), and E at s = 300, 450, 600."""
    out = {}
    cache = {}
    for r in runs:
        if r['leak'] != 'nearest':
            continue
        net, loop, atoms, res, marg = setup(r['prior'])
        for m, sp in r['spurious'].items():
            key = (r['prior'], tuple(int(x) for x in sp['counts']))
            if key not in cache:
                cache[key] = weight_admissibility(loop, atoms, sp['counts'])
            A = cache[key]
            Th = np.asarray(r['final'][m], float)
            lab = loop.nearest.assignment(Th)
            spread = [float(np.sqrt(np.mean(np.sum((Th[lab == k] - Th[lab == k].mean(axis=0)) ** 2 / G[None, :], axis=1)))) for k in range(3)]
            far = float(np.sqrt(np.sum((Th - A['Tw'][lab]) ** 2 / G[None, :], axis=1)).max())
            rc = r['rec'][m]
            Es = [at(rc['s'], rc['E'], t) for t in (300.0, 450.0, 600.0)]
            out[(r['seed'], r['prior'], m)] = dict(counts=sp['counts'], sep=A['sep'], hess_min=A['hess_min'], own_cells=A['own_cells'],
                                                  admissible=A['admissible'], res=A['res'], spread=spread, far=far, E_t=Es, E_spur=sp['E'])
    return out


def local_pl_constants():
    """At μ* for each prior: λ_sep = min_k λ_min(σ₀I − R^G_k[μ*]) and λ_coll = λ_min of the whitened weighted K-atom Hessian of the OT energy
    (the diagonal plan), so that the local PL constant near μ* is 2min{λ_sep, λ_coll}."""
    out = {}
    for name in SHIFTS:
        net, loop, atoms, res, marg = setup(name)
        B = OB.BasinDiagnostics(loop, atoms)
        lam_sep = float(np.min(SIGMA0 - loop.kappa_R(B.r_s, atoms)))
        lam_coll = float(k_atom_hessian(loop.nearest, Q, atoms)[0])
        out[name] = dict(lam_sep=lam_sep, lam_coll=lam_coll, local_pl=2 * min(lam_sep, lam_coll))
    return out


def cell_distance(loop, atoms):
    """The smallest G⁻¹-distance from a stationary atom ϑ*_j to the boundary of the prior's Voronoi cell of ϑ⁰_j (its bisectors with the other
    prior atoms): a particle of class j outside that cell is at least this far from ϑ*_j."""
    P, d = loop.prior, np.inf
    for j in range(P.shape[0]):
        for k in range(P.shape[0]):
            if k != j:
                n2 = float(np.sum((P[k] - P[j]) ** 2 / G))
                d = min(d, (float(np.sum((atoms[j] - P[k]) ** 2 / G)) - float(np.sum((atoms[j] - P[j]) ** 2 / G))) / (2 * np.sqrt(n2)))
    return d


def sup_R_global(loop, r, n_a=61, n_c=241, n_s=49, t_max=6.0, ts_max=12.0, n_polish=16):
    """sup over the whole parameter space of λ_max R^G[μ](ϑ): a grid over the saturation ranges (a = a_max tanh(t_a), |t_a| ≤ t_max; c on n_c
    equispaced values of (−c_max, c_max); s = s_min + (s_max − s_min)σ(t_s), |t_s| ≤ ts_max), a scan along c at the four saturation corners
    (t_a = ±3t_max, t_s = ±3ts_max), then Nelder–Mead (unconstrained) from the best n_polish grid points at least 0.5 apart in (t_a, c/c_max, t_s).
    The supremum may be approached only at infinity (a saturation corner); the value returned is the largest found."""
    from scipy.optimize import minimize
    net = loop.net
    ta = np.linspace(-t_max, t_max, n_a); cs = np.linspace(-1.0, 1.0, n_c + 2)[1:-1] * net.c_max; ts = np.linspace(-ts_max, ts_max, n_s)
    A, C, S = np.meshgrid(ta, cs, ts, indexing='ij')
    X = np.column_stack([net.a_max * A.ravel(), net.c_max * np.arctanh(C.ravel() / net.c_max), S.ravel()])   # unconstrained coordinates
    corners = [np.column_stack([np.full(n_c, net.a_max * sa * 3 * t_max), net.c_max * np.arctanh(cs / net.c_max), np.full(n_c, ss * 3 * ts_max)])
               for sa in (-1, 1) for ss in (-1, 1)]
    X = np.vstack([X] + corners)
    f = lambda Y: np.linalg.eigvalsh(SQG[None, :, None] * loop.residual_curvature(r, Y) * SQG[None, None, :])[:, -1]
    vals = np.concatenate([f(X[i:i + 4000]) for i in range(0, X.shape[0], 4000)])
    best = float(vals.max()); x_best = X[int(np.argmax(vals))]
    key = np.column_stack([X[:, 0] / net.a_max, np.tanh(X[:, 1] / net.c_max), X[:, 2]])
    starts = []
    for i in np.argsort(vals)[::-1]:
        if all(np.abs(key[i] - key[j]).max() > 0.5 for j in starts):
            starts.append(i)
        if len(starts) == n_polish:
            break
    for i in starts:
        res = minimize(lambda v: -float(f(v[None, :])[0]), X[i], method='Nelder-Mead', options=dict(maxiter=800, xatol=1e-7, fatol=1e-13))
        if -res.fun > best:
            best, x_best = float(-res.fun), res.x
    a, c, s_, _, _ = loop.features(x_best[None, :])
    return best, (float(a[0]), float(c[0]), float(s_[0]))


def a_priori_constants(runs):
    """L₂^G (search, as in exp_chaos_learning.py) and, for each prior, κ := σ₀ − sup λ_max R^G[μ*] (a) over the whole parameter space and (b)
    over the bounding box of the end configurations, the prior and the stationary atoms, enlarged by 10% (20000 random points, L-BFGS-B polish
    from the best 12); E₀ = κ²/(2(L₂^G)²), and the largest fraction p of the particles outside their cells that {E ≤ E₀} can hold,
    p ≤ 2E₀/(κd²) with d the cell distance (E ≥ (κ/2)D_π under (H*))."""
    from scipy.optimize import minimize
    net = make_net()
    Lc = CL.whitened_lipschitz_constants(net, G, 2.5, n_samples=60000)
    out = dict(L2G=Lc['L_2_G'])
    for name in SHIFTS:
        net, loop, atoms, res, marg = setup(name)
        B = OB.BasinDiagnostics(loop, atoms)
        pts = [atoms, loop.prior] + [np.asarray(r['final'][m], float) for r in runs if r['prior'] == name for m in r['final']]
        P = np.vstack(pts)
        lo, hi = P.min(axis=0), P.max(axis=0); pad = 0.1 * (hi - lo)
        rng = np.random.default_rng(11)
        X = rng.uniform(lo - pad, hi + pad, (20000, 3))
        f = lambda Y: np.linalg.eigvalsh(SQG[None, :, None] * loop.residual_curvature(B.r_s, Y) * SQG[None, None, :])[:, -1]
        vals = f(X)
        sup_box = float(vals.max())
        for i in np.argsort(vals)[::-1][:12]:
            rr = minimize(lambda v: -float(f(v[None, :])[0]), X[i], method='L-BFGS-B', bounds=list(zip(lo - pad, hi + pad)))
            sup_box = max(sup_box, float(-rr.fun))
        sup_glob, where = sup_R_global(loop, B.r_s)
        d_cell = cell_distance(loop, atoms)
        entry = dict(box=(lo - pad, hi + pad), d_cell=d_cell, sup_glob_at=where)
        for tag, sup in (('box', sup_box), ('glob', sup_glob)):
            kappa = SIGMA0 - sup
            E0 = kappa ** 2 / (2 * Lc['L_2_G'] ** 2) if kappa > 0 else float('nan')
            pmax = 2 * E0 / (kappa * d_cell ** 2) if kappa > 0 else float('nan')
            entry[tag] = dict(sup=sup, kappa=kappa, E0=E0, p_max=pmax, m_min=(1.0 / pmax if kappa > 0 else float('nan')))
        out[name] = entry
    out['local'] = local_pl_constants()
    return out


def moved_fraction(ms=(30, 300), ks=(1, 2, 3, 5), t_move=0.6, s_run=300.0, s_long=3000.0, ds=0.5, n_t_table=129, n_t_run=17):
    """The measure-level basin at the matching prior: the lift of q with k particles of class 0 (k = 1, 2, 3, 5 at m = 30, ten times as many at
    m = 300: the same fractions, the same measures) moved t_move of the way to atom 1, into cell 1. At s = 0: the plan's labels of the moved
    particles, κ^Θ, the margin α and λ, the segment suprema on n_t_table points per segment. With one particle moved: the OT run over s_run slow
    units (midpoint rule, plan every step) with λ evaluated at every step (n_t_run points per segment) and ∫λ₊ by the smaller endpoint value
    on each step — a numerical evaluation, not a bound; the nearest-atom run over s_long, E at s = 300, 1000, s_long, against the replicated
    configuration of its counts (the spurious minimum: its E and admissibility)."""
    net, loop, atoms, res, marg = setup('matching')
    B = OB.BasinDiagnostics(loop, atoms)
    rows, runs = [], {}
    for m in ms:
        lab0 = np.repeat(np.arange(3), m // 3)
        for k0 in ks:
            k = k0 * m // ms[0]
            Th = atoms[lab0].copy()
            Th[:k] = (1.0 - t_move) * atoms[0] + t_move * atoms[1]
            t = B.terms(Th)
            sc = B.segment_curvatures(Th, t['lab'], n_t=n_t_table)
            al, lam = B.theorem_rate(**sc)
            rows.append(dict(m=m, k=k, p=k / m, kept=bool(np.all(t['lab'][:k] == 0)), E=t['E'], D=t['D'], ktheta=sc['kappa_theta'], alpha=al, lam=lam, ratio=t['ratio']))
        for leak in ('ot', 'nearest'):
            Th = atoms[lab0].copy(); Th[0] = (1.0 - t_move) * atoms[0] + t_move * atoms[1]
            t = B.terms(Th); al, lam = B.theorem_rate(**B.segment_curvatures(Th, t['lab'], n_t=n_t_run))
            E0, integral, lam_prev, lam_min = t['E'], 0.0, lam, lam
            n = int(round((s_run if leak == 'ot' else s_long) / ds))
            E_at = {}
            for j in range(n):
                T = loop.leak_target(Th) if leak == 'ot' else loop.nearest.leak_target(Th)
                k1 = loop.drift(Th, T); Th = Th + ds * loop.drift(Th + 0.5 * ds * k1, T)
                if leak == 'ot':
                    t = B.terms(Th); al, lam = B.theorem_rate(**B.segment_curvatures(Th, t['lab'], n_t=n_t_run))
                    integral += ds * min(lam_prev, lam); lam_prev = lam; lam_min = min(lam_min, lam)
                s_now = (j + 1) * ds
                if leak == 'nearest' and s_now in (300.0, 1000.0, s_long):
                    E_at[s_now] = B.terms(Th)['E']
            t = B.terms(Th)
            out = dict(E0=E0, E=t['E'], ratio=t['ratio'], lab0=int(t['lab'][0]), cell=int(np.argmin(loop.cost_matrix(Th[:1])[0])))
            if leak == 'ot':
                out.update(bound=E0 * np.exp(-integral), lam_min=lam_min)
            else:
                counts = np.bincount(loop.nearest.assignment(Th), minlength=3)
                A = weight_admissibility(loop, atoms, counts)
                lift = np.repeat(A['Tw'], counts, axis=0)
                out.update(E_at=E_at, counts=counts, E_spur=B.terms(lift)['E'], admissible=A['admissible'], sep_min=float(np.min(A['sep'])), hess_min=A['hess_min'])
            runs[(m, leak)] = out
    return rows, runs


# ---------------------------------------------------------------------------- report
def at(s, v, t):
    return float(v[int(np.argmin(np.abs(s - t)))])


def fmt_ms(vals, p=2):
    vals = np.asarray(vals, float)
    return f"{vals.mean():.{p}e}" + (f" ({vals.min():.1e}–{vals.max():.1e})" if len(vals) > 1 else "")


def report(runs, consts, adm=None, moved=None, t0=None):
    lines = ["# Theorem D‴.21: a measure-level basin for the OT leak — the Polyak–Łojasiewicz inequality along the runs (`exp_ot_basin.py`)", "",
             f"The averaged slow loop on the three-bump instance at σ₀ = 0.15σ_crit = {SIGMA0:.5f} (2σ₀ = {2 * SIGMA0:.5f}); sizes m = {', '.join(map(str, MS))} and M = {M_REF} from the warm clouds "
             f"of Proposition D‴.12 (width {JITTER:g} around the exact representation; the m-system the first m rows); {N_SEEDS} seeds; {S_END:g} slow units, midpoint rule dt = {DT:g}, leak target every "
             f"{PLAN_EVERY} steps. E = 𝒥(μ_Θ) − 𝒥* with 𝒥 = L̄ + (σ₀/2)W²(·, ν) the OT functional, evaluated along both leaks' trajectories; D_a the class-respecting distance to μ*; the observed PL ratio "
             "|∇J_a|²_G/(mE), which along the OT flow is −d log E/ds; λ the theorem's rate from the plan-segment curvatures (9 points per segment, every 25 slow units). "
             "Entries: mean over the seeds (min–max).", ""]
    for name in SHIFTS:
        rr = [r for r in runs if r['prior'] == name]
        r0 = rr[0]
        lin = None
        lines += [f"## Prior: {name} (centres shifted by {SHIFTS[name]:g})", "",
                  f"μ*: Newton residual {r0['newton_res']:.1e}, stationarity of the atoms {r0['stationarity']:.1e}, 𝒥* = {r0['Jstar']:.4e}; diagonal plan: ε = {r0['diag_margin'][0]:.3f} < Δ_ν/2 = {r0['diag_margin'][1]:.3f}. "
                  f"Identities (E), (P) along every run: largest relative residual {max(float(np.max(r['rec'][m]['idE'])) for r in rr for m in r['rec']):.1e} (E) and {max(float(np.max(r['rec'][m]['idQ'])) for r in rr for m in r['rec']):.1e} (P).", "",
                  "| leak | m | E at s = " + ", ".join(f"{t:g}" for t in TIMES) + " | D_a at s = 0, 100, 600 | W²(μ_Θ, μ*) at s = 0, 100, 600 | plan changes (particle relabellings) |",
                  "|---|---|---|---|---|---|"]
        for leak in ('ot', 'nearest'):
            for m in list(MS) + [M_REF]:
                rs = [r for r in rr if r['leak'] == leak]
                Es = [" / ".join(f"{at(r['rec'][m]['s'], r['rec'][m]['E'], t):.2e}" for r in rs) for t in TIMES]
                Em = [np.mean([at(r['rec'][m]['s'], r['rec'][m]['E'], t) for r in rs]) for t in TIMES]
                Dm = [np.mean([at(r['rec'][m]['s'], r['rec'][m]['D'], t) for r in rs]) for t in (0.0, 100.0, 600.0)]
                Wm = [np.mean([at(r['rec'][m]['s'], r['rec'][m]['W2s'], t) for r in rs]) for t in (0.0, 100.0, 600.0)]
                sw = [int(np.sum(r['rec'][m]['switches'])) for r in rs]
                lines.append(f"| {leak} | {m} | " + ", ".join(f"{x:.2e}" for x in Em) + " | " + ", ".join(f"{x:.2e}" for x in Dm) + " | " + ", ".join(f"{x:.2e}" for x in Wm)
                             + f" | {min(sw)}–{max(sw)} |")
        # the PL inequality along the OT runs
        lines += ["", "The PL inequality along the OT runs. s₊: the first segment record after which the margin α stays positive; the observed ratio and the bounds over the records "
                  "with α > 0 (the chain ratio ≥ Q²/(DE) ≥ λ must hold); the late values: means over s ∈ [500, 600] (records with E > 10⁻²⁶).", "",
                  "| m | s₊ (seeds) | E(s₊)/E(0) | min over α > 0 of ratio/λ and Q²/(DE)/λ | late observed ratio | late λ | late α | late κ*, κ^Θ | a-posteriori bound E(0)e^{−∫λ₊} at s = 600 against E(600) |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for m in list(MS) + [M_REF]:
            rs = [r for r in rr if r['leak'] == 'ot']
            splus, Erat, mins, late_ratio, late_lam, late_al, late_k, bounds = [], [], [], [], [], [], [], []
            for r in rs:
                sg, rc = r['seg'][m], r['rec'][m]
                pos = sg['alpha'] > 0
                idx = None
                for j in range(len(pos)):
                    if np.all(pos[j:]):
                        idx = j; break
                sp = float(sg['s'][idx]) if idx is not None else float('nan')
                splus.append(sp)
                Erat.append(at(rc['s'], rc['E'], sp) / rc['E'][0] if idx is not None else float('nan'))
                # ratio and lam_eff at the segment records
                ok = []
                for j, s_ in enumerate(sg['s']):
                    if sg['alpha'][j] > 0 and sg['lam'][j] > 0:
                        i = int(np.argmin(np.abs(rc['s'] - s_)))
                        if rc['E'][i] > 1e-26:
                            ok.append((rc['ratio'][i] / sg['lam'][j], rc['lam_eff'][i] / sg['lam'][j]))
                mins.append((min(x[0] for x in ok), min(x[1] for x in ok)) if ok else (float('nan'), float('nan')))
                late = (rc['s'] >= 500.0) & (rc['E'] > 1e-26)
                late_ratio.append(float(np.mean(rc['ratio'][late])) if late.any() else float('nan'))
                lateg = sg['s'] >= 500.0
                late_lam.append(float(np.mean(sg['lam'][lateg]))); late_al.append(float(np.mean(sg['alpha'][lateg])))
                late_k.append((float(np.mean(sg['kstar'][lateg])), float(np.mean(sg['ktheta'][lateg]))))
                # a-posteriori bound: λ interpolated linearly between the segment records, ∫ by the trapezoid rule
                lam_t = np.interp(rc['s'], sg['s'], np.maximum(sg['lam'], 0.0))
                integral = float(np.trapezoid(lam_t, rc['s'])) if hasattr(np, 'trapezoid') else float(np.trapz(lam_t, rc['s']))
                bounds.append((rc['E'][0] * np.exp(-integral), rc['E'][-1]))
            lines.append(f"| {m} | " + ", ".join(f"{x:g}" for x in splus) + " | " + ", ".join(f"{x:.1e}" for x in Erat) + " | "
                         + f"{min(x[0] for x in mins):.2f}, {min(x[1] for x in mins):.2f} | {np.nanmean(late_ratio):.5f} | {np.mean(late_lam):.5f} | {np.mean(late_al):.5f} | "
                         + f"{np.mean([k[0] for k in late_k]):.1e}, {np.mean([k[1] for k in late_k]):.1e} | " + "; ".join(f"{b:.1e} ≥ {e:.1e}" for b, e in bounds) + " |")
        # nearest-atom spurious minima, with the admissibility of the sampled counts
        lines += ["", "The nearest-atom contrast: E at the end of the nearest-atom runs against the E of the replicated configuration of their sampled weights (a critical point of the "
                  "nearest-atom energy; its drift under the nearest-atom leak and under the OT leak), and the admissibility of those weights: the smallest separation margin σ₀ − λ_max R^G_k(w) "
                  "(the class), the smallest eigenvalue of the whitened weighted K-atom Hessian; inadmissible weights (a margin ≤ 0) have a saddle as their replicated configuration. "
                  "E at s = 300, 450, 600; the end state's largest class spread (rms G⁻¹-distance to the class mean) and its largest distance from the atom of Θ*(w).", "",
                  "| m | seed | sampled counts | admissible | smallest separation margin (class) | K-atom Hessian min | E at s = 300, 450, 600 | E of the replicated configuration | E(600)/E_spur − 1 | class spread, largest distance | its drift: nearest / OT |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]
        for m in list(MS) + [M_REF]:
            for r in sorted([r for r in rr if r['leak'] == 'nearest'], key=lambda r: r['seed']):
                if m not in r['spurious']:
                    continue
                sp = r['spurious'][m]; ad = adm.get((r['seed'], name, m)) if adm else None
                if ad is None:
                    lines.append(f"| {m} | {r['seed']} | {tuple(int(x) for x in sp['counts'])} | — | — | — | — | {sp['E']:.4e} | — | — | {sp['drift_nearest']:.0e} / {sp['drift_ot']:.1e} |")
                    continue
                kmin = int(np.argmin(ad['sep']))
                lines.append(f"| {m} | {r['seed']} | {tuple(int(x) for x in sp['counts'])} | {'yes' if ad['admissible'] else 'no (saddle)'} | {ad['sep'][kmin]:+.5f} ({kmin}) | {ad['hess_min']:.5f} | "
                             + ", ".join(f"{x:.4e}" for x in ad['E_t']) + f" | {sp['E']:.4e} | {ad['E_t'][-1] / sp['E'] - 1:+.4f} | {max(ad['spread']):.3f}, {ad['far']:.3f} | "
                             + f"{sp['drift_nearest']:.0e} / {sp['drift_ot']:.1e} |")
        c = consts[name]
        lam_late = float(np.mean([np.mean(r['seg'][m]['lam'][r['seg'][m]['s'] >= 500.0]) for r in rr if r['leak'] == 'ot' for m in r['seg']]))
        if 'glob' in c:
            g_, b_ = c['glob'], c['box']
            lines += ["", f"The a-priori basin of Theorem D‴.21(iv). Over the whole parameter space sup λ_max R^G[μ*] = {g_['sup']:.4e} (approached at a = {c['sup_glob_at'][0]:.3f}, "
                      f"c = {c['sup_glob_at'][1]:.3f}, s = {c['sup_glob_at'][2]:.3f})" + (f", so (H*) holds on the whole space with κ = {g_['kappa']:.3e} and E₀ = κ²/(2(L₂^G)²) = {g_['E0']:.2e} "
                      f"(L₂^G = {consts['L2G']:.3f})" if np.isfinite(g_['E0']) else ": (H*) fails on the whole space") + f". Over the box of the end configurations, the priors and the atoms "
                      f"(not forward-invariant) sup = {b_['sup']:.4e}, κ = {b_['kappa']:.5f}, E₀ = {b_['E0']:.2e}. A particle outside its cell is at G⁻¹-distance ≥ {c['d_cell']:.4f} from its class "
                      f"atom, so {{E ≤ E₀}} holds a fraction p ≤ {g_['p_max']:.1e} (whole space) or {b_['p_max']:.1e} (box) of the particles in other cells: one such particle needs m ≥ "
                      f"{g_['m_min']:.1e} or {b_['m_min']:.1e}. Local PL constant at μ*: 2min{{λ_sep, λ_coll}} = {2 * min(consts['local'][name]['lam_sep'], consts['local'][name]['lam_coll']):.5f} "
                      f"(λ_sep = {consts['local'][name]['lam_sep']:.5f}, λ_coll = {consts['local'][name]['lam_coll']:.5f}); the theorem's late λ (mean over the OT runs, s ≥ 500) "
                      f"is {lam_late / (2 * SIGMA0):.3f} of 2σ₀ = {lam_late / (2 * SIGMA0) / (consts['local'][name]['local_pl'] / (2 * SIGMA0)):.3f} (the slack of (ii)) × "
                      f"{consts['local'][name]['local_pl'] / (2 * SIGMA0):.3f} (the curvature's cost, the local constant over 2σ₀).", ""]
        else:
            lines += ["", f"The a-priori basin of Theorem D‴.21(iv): sup λ_max R^G[μ*] over the box of the runs = {c['sup_kappa_R']:.4e}, so κ = σ₀ − sup = {c['kappa']:.5f}"
                      + (f" and E₀ = κ²/(2(L₂^G)²) = {c['E0']:.2e} (L₂^G = {consts['L2G']:.3f})." if np.isfinite(c['E0']) else " ≤ 0: hypothesis (H*) fails on that box, and only the a-posteriori form applies."), ""]
    rows, mruns = moved if moved is not None else moved_fraction()
    lines += ["## The measure-level basin: particles moved into another cell (matching prior)", "",
              "The lift of q with k particles of class 0 moved 0.6 of the way to atom 1 (into cell 1); the same fractions p = k/m at m = 30 and 300 (the same measures, so the same values). At s = 0: "
              "whether the plan keeps the moved particles in class 0, E, D_π, the current-residual curvature κ^Θ along the plan segments (129 points per segment), the margin α and the rate λ "
              "(κ* = 0 at the matching prior).", "",
              "| m | k | p | kept in class 0 | E | D_π | κ^Θ | α | λ | observed ratio |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['m']} | {r['k']} | {r['p']:.4f} | {r['kept']} | {r['E']:.3e} | {r['D']:.3e} | {r['ktheta']:.3e} | {r['alpha']:+.5f} | {r['lam']:.5f} | {r['ratio']:.5f} |")
    lines += ["", "One particle moved (midpoint rule, ds = 0.5, plan every step): ", ""]
    for (m, leak), r in mruns.items():
        if leak == 'ot':
            lines.append(f"- m = {m}, OT leak, 300 slow units: E {r['E0']:.3e} → {r['E']:.3e} (E(0)exp(−∫λ₊) = {r['bound']:.3e}, λ evaluated at every step, smallest λ {r['lam_min']:.5f}); "
                         f"observed ratio at the end {r['ratio']:.5f}; the moved particle labelled {r['lab0']} by the OT plan, in cell {r['cell']} at the end.")
        else:
            lines.append(f"- m = {m}, nearest-atom leak, {max(r['E_at']):g} slow units (E is the OT functional's excess, not this leak's Lyapunov function): E {r['E0']:.3e} → "
                         + ", ".join(f"{v:.4e} (s = {k:g})" for k, v in sorted(r['E_at'].items())) + f"; the particle in cell {r['cell']}; the counts {tuple(int(x) for x in r['counts'])}, whose "
                         f"replicated configuration has E = {r['E_spur']:.4e} ({'admissible' if r['admissible'] else 'inadmissible'}: smallest separation margin {r['sep_min']:+.5f}, K-atom Hessian min {r['hess_min']:.5f}).")
    lines += ["", "Numbers are produced by `experiments/exp_ot_basin.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'ot_basin.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'ot_basin_raw.npz')
    if '--from-raw' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        report(list(z['runs']), z['consts'].item(), z['adm'].item() if 'adm' in z.files else None,
               tuple(z['moved']) if 'moved' in z.files else None, t0)
        return
    if '--only' in sys.argv:
        parts = sys.argv[sys.argv.index('--only') + 1].split(',')
        z = np.load(raw, allow_pickle=True)
        data = {k: z[k] for k in z.files}
        runs = list(z['runs'])
        if 'consts' in parts:
            data['consts'] = np.array(a_priori_constants(runs), dtype=object)
        if 'adm' in parts:
            data['adm'] = np.array(admissibility(runs), dtype=object)
        if 'moved' in parts:
            mv = moved_fraction(); arr = np.empty(2, dtype=object); arr[0], arr[1] = mv
            data['moved'] = arr
        np.savez(raw, **data)
        report(runs, data['consts'].item(), data['adm'].item() if 'adm' in data else None, tuple(data['moved']) if 'moved' in data else None, t0)
        return
    jobs = [(seed, name, leak) for seed in range(N_SEEDS) for name in SHIFTS for leak in ('ot', 'nearest')]
    with mp.Pool(N_WORKERS) as pool:
        runs = pool.map(job, jobs)
    consts = a_priori_constants(runs)
    adm = admissibility(runs)
    mv = moved_fraction(); moved = np.empty(2, dtype=object); moved[0], moved[1] = mv
    np.savez(raw, runs=np.array(runs, dtype=object), consts=np.array(consts, dtype=object), adm=np.array(adm, dtype=object), moved=moved)
    report(runs, consts, adm, mv, t0)


if __name__ == '__main__':
    main()
