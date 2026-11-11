"""
Which factorised readouts of a Gaussian-kernel network, trained by plain gradient descent from a small start, reproduce the sliding
transit? (the generality of Paper 17's mechanism)
====================================================================================================================================

The free law is the Euclidean gradient flow in (α, ν, c) with the learning rates (Gm, γm/4, Gm) of a standard RBF network with the
readout w = ν²a_max tanh(α/a_max)/m (`factorized_rbf.PaperLaw`). Here the same network is trained with other readouts, by the same
gradient flow, fixed-step gradient descent and Adam, and compared with the sliding inverse scale space (`rich_transit_nd.sliding_iss_nd`,
the ε → 0 limit): the formed set (|w| > 0.01 max|A*|), the order and the times of the entries, the gap between formed and dormant
coefficients, the dormant units' displacement, exactness (one unit within 0.05σ and 5% of each feature). D‴.48's dictionary (21 kernels
of width 0.12 on a grid of spacing 0.1, 61 samples of [−1, 1]), centre rate Gm = 12 in every readout.

A. The reformulation: the law's right-hand side against the Euclidean gradient flow's (an identity), the raw flow against the reduced
   one (`rich_transit.critical_rhs` on the curves), the conserved quantities along the raw flows.
B. The design of D‴.48 (features −0.52, 0.07, 0.61 with amplitudes 0.9, −0.6, 0.75): every readout over ε = 10⁻² … 10⁻¹⁰ — the paper
   law P0, without saturation P0u, with equal learning rates P4, depth 3 D3, depth 2 D2, depth 1 D1 (the control).
C. 200 of the 1034 one-dimensional targets of `exp_transit_certificate` (a seeded stratified sample of its families A, B, D): P0, D3
   and D2 against the sliding path of each target.
D. Depth 2's convergence in ε: the logarithmic scale — ε = 10⁻⁶ … 10⁻¹⁴⁴ on the design and the first targets of C.
E. The learning rates: P0 at lr_ν/lr_α = γ/(4G) = 0.025, 0.25, 1, and D3 at lr = 0.1, 1, 10.
F. Fixed-step gradient descent (explicit Euler on the same flow) on the design and a few targets.
G. Random starts: random masses ε b_i (b_i log-normal; the path with ρ_i = b_i^{−1/2}) and random initial amplitudes α_i(0).
H. Adam (reported, not claimed).
I. The centres' rate: 0 (fixed centres), 0.1, 1 and 10 times Gm — is sliding what makes K units enough?
J. The plane: 12 of the 120 planar targets of `exp_transit_certificate` (D‴.54's lattice: 17 × 17 kernels of width 0.15, 41² samples).
K. Branch pairs: the paper's logarithmic regime — two kernels per atom with α(0) = +a₀ and −a₀ (`factorized_rbf.PaperPairs`), the design
   and the first targets of C's sample, a₀ = 0.1, 0.3 and ε = 10⁻⁶ … 10⁻⁹⁶ (P0, D3 only at 10⁻⁶, 10⁻¹²), against the sliding path of the 21 atoms,
   with P0, D3, D2 as controls; the scan of a₀ on the design; the coupled targets (the 16 targets of C's sample whose two features are at most
   2σ apart: where the pairs' common early factor holds for the first entrant only).
L. The dictionary's spacing: 21, 16, 13, 11, 9, 7 kernels (spacing 0.1 … 0.33) on the design's off-grid features and, as a control that removes
   the sliding, with the features moved onto kernel centres; P0, D3, D2 at ε = 10⁻⁶, 10⁻¹², each against the path of its own dictionary; and
   the scan of ε = 10⁻⁶ … 10⁻⁹⁰ of D2 on six dictionaries (11 … 31 kernels): where does depth 2 start to pass, against the margin rule?

Usage: python experiments/exp_parametrization.py [--only A,B,C,D,E,F,G,H,I,J,K,L] [--from-raw] → results/chaos/parametrization.md (+ _raw.npz).
Environment: PZ_OUT (output directory), PZ_WORKERS (processes, default 3).
"""
import copy
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import factorized_rbf as FR
from src.core import rich_transit as RT
from src.core import rich_transit_nd as RN
from src.core.mirror_free import FreeLaw
import exp_cohesion_leak as E

OUT = os.environ.get('PZ_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
CERT_RAW = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos', 'transit_certificate_raw.npz')
N_WORKERS = int(os.environ.get('PZ_WORKERS', '3'))

LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
NS, SIG, ND = 61, 0.12, 21
CSTAR, ASTAR = (-0.52, 0.07, 0.61), (0.9, -0.6, 0.75)
LR_C = LAW['G'] * LAW['m']                                              # the centres' rate Gm, in every readout
EPS_B = (1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
TAU_B = 20.0
SEED_C, N_C = 901, dict(A_nearest=35, A_other=35, B=30, Da=20, Db=50, Dc=30)
CELLS_C = (('P0', 1e-6), ('P0', 1e-8), ('D3', 1e-6), ('D3', 1e-8), ('D2', 1e-6), ('D2', 1e-12))
EPS_D = (1e-6, 1e-12, 1e-24, 1e-48, 1e-96, 1e-144)
N_D = 20                                                                # D: the first N_D targets of C's sample
GAMMAS_E, LRS_E, N_E, EPS_E = (0.1, 1.0, 4.0), (0.1, 1.0, 10.0), 10, 1e-6
CELLS_F = (('P0', 1e-4, 0.02), ('P0', 1e-4, 0.05), ('P0', 1e-6, 0.05), ('D3', 1e-4, 0.02), ('D3', 1e-6, 0.02), ('D3', 1e-6, 0.05))
N_F = 4                                                                 # F: also this many targets of C's sample (D3 at 1e-6 and P0 at 1e-4, h = 0.02)
RAND_G, AMP0_G, N_G, EPS_G, SEED_G = 0.5, (1e-3, 1e-2), 10, 1e-6, 903
LRS_H, STEPS_H, EPS_H = (1e-3, 1e-4, 1e-5, 1e-6), 200000, 1e-6
LRC_I, EPS_I, N_I = (0.0, 0.1, 1.0, 10.0), 1e-6, 10                    # I: the centres' rate in units of Gm; ε; extra targets of C's sample
N_J, SEED_J, CELLS_J = 12, 905, (('P0', 1e-6), ('D3', 1e-6), ('D2', 1e-6))        # J: planar targets sampled from the raw's family C
A0_K, EPS_K, N_K, NAMES_K = (0.1, 0.3), (1e-6, 1e-12), 10, ('P0', 'D3', 'D2')       # K: the pairs' a₀ and ε; extra targets of C's sample; the controls
EPS_K_FAR = (1e-24, 1e-48, 1e-96)                                       # K: further down, for the pairs and D2 (P0, D3 on the power scale ε^(−1/2): not integrated there)
SEP_KC, A0_KC, EPS_KC = 2.0 * SIG, 0.1, (1e-6, 1e-12, 1e-24, 1e-48)    # K: the coupled targets (features at most 2σ apart); their pairs' a₀ and ε
A0_SCAN_K, EPS_SCAN_K = (0.003, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0, 3.0), (1e-6, 1e-12)    # K: the scan of a₀ on the design
ND_L, ND_ON_L, EPS_L, NAMES_L = (21, 16, 13, 11, 9, 7), (21, 16, 13, 11, 9, 7), (1e-6, 1e-12), ('P0', 'D3', 'D2')   # L: kernels on [−1, 1]; on-grid arm at all
ND_SCAN_L, K_SCAN_L = (11, 13, 16, 21, 25, 31), tuple(range(6, 91, 2))  # L: the scan of D2 over ε = 10^(−k): dictionaries (kernels) and k
MAX_SECS = 900                                                          # a safety net (s), checked between the 8 pieces of a run; above every cell's time
TOL_TIE = 1e-9


def law(**kw):
    return FreeLaw(**dict(LAW, **kw))


def setup():
    warnings.filterwarnings('ignore')
    xs = np.linspace(-1, 1, NS)
    return law(), RN.BumpsND(xs[:, None], SIG), np.linspace(-1, 1, ND)


def make_readout(name):
    """P0 the paper's unit, P0u without saturation, P4 with equal learning rates (γ = 4G, v_s = a_max²), D3 depth 3, D2 depth 2, D1 depth 1."""
    if name == 'P0':
        return FR.PaperLaw(law())
    if name == 'P0u':
        return FR.Unsaturated(law())
    if name == 'P4':
        return FR.PaperLaw(law(gamma=4 * LAW['G']))
    if name == 'D3':
        return FR.Cubic3(1.0)
    if name in ('D4', 'D5'):
        return FR.PowerN(int(name[1:]), 1.0)
    if name == 'D2':
        return FR.Product2(1.0)
    if name == 'D1':
        return FR.Linear(LR_C)
    raise ValueError(name)


def margins(K, f, R, d):
    """The runner-up margins of the sliding path R (equal masses): the smallest 1 − |u_j| over the kernels not yet active at an entry (the later
    entrants included: `margin`) and over the kernels that never enter (`margin_never`), from the duals at the entries (`path_duals_nd`).
    A limit statement about the loop at a mass ε has to be read against it (Paper 17: the limit is not uniform in the margin)."""
    P = RN.path_duals_nd(K, f, R)
    taus, U, ev = P['taus'], np.abs(P['U']), R['events']
    ever = np.zeros(d, bool)
    ever[P['ever']] = True
    act, e, m_all, m_nev = np.zeros(d, bool), 0, np.inf, np.inf
    for r in range(1, len(taus)):
        while e < len(ev) and ev[e][0] <= taus[r] + 1e-12:
            act[ev[e][1]] = ev[e][2] == 'enter'
            e += 1
        inact = ~act
        if inact.any():
            m_all = min(m_all, float((1 - U[r][inact]).min()))
            if (inact & ~ever).any():
                m_nev = min(m_nev, float((1 - U[r][inact & ~ever]).min()))
    nan = float('nan')
    return (m_all if np.isfinite(m_all) else nan), (m_nev if np.isfinite(m_nev) else nan)


def rule_margin(m):
    """The margin the depth-2 rule reads: μ never, the smallest 1 − |u_j| over the kernels that never enter. NaN (every kernel enters, so no kernel
    is left to violate the rule; `margins`) is read as +∞: the rule then predicts that depth 2 forms the path's set."""
    return float('inf') if m is None or np.isnan(m) else float(m)


def rule_threshold(eps):
    """ln(1/θ)/ln(1/ε), θ = `FR.THETA`: the margin that the rule asks a kernel that never enters to exceed (a coefficient at margin μ is ε^μ)."""
    return float(np.log(1 / FR.THETA) / np.log(1 / eps))


def rule_passes(m_never, eps):
    """The depth-2 rule: depth 2 forms the path's set iff μ never > ln(1/θ)/ln(1/ε) (NaN as +∞)."""
    return rule_margin(m_never) > rule_threshold(eps)


def reference(L, K, c0, f, rho=None):
    """The sliding path of the target f: the entrants in order with their times, the final active set, exits, the rms residual, its solver flags
    and (equal masses) the runner-up margins."""
    R = RN.sliding_iss_nd(L, K, f, np.asarray(c0, float)[:, None], np.ones(len(c0)) if rho is None else rho, polish=True)
    ent = [(float(t), int(j)) for t, j, k, c in R['events'] if k == 'enter']
    m_all, m_nev = margins(K, f, R, len(c0)) if rho is None else (float('nan'), float('nan'))
    return dict(order=[j for t, j in ent], tau=[t for t, j in ent], set=sorted(int(j) for j in np.flatnonzero(R['act'])),
                exits=int(sum(k == 'exit' for t, j, k, c in R['events'])), rms=float(R['hist'][-1][3]), margin=m_all, margin_never=m_nev,
                relaxed=bool(all(R['relaxed'])), solver_failed=int(R.get('solver_failed', 0)), out_of_budget=bool(R.get('out_of_budget', False)))


def cell(ro, K, f, c0, cs, As, eps, ref, tau_end, kind='flow', lr_c=None, **kw):
    """One training run and its analysis against the path (a flat dict of numbers and short lists). lr_c: the centres' rate (default Gm)."""
    lr_c = LR_C if lr_c is None else lr_c
    try:
        if kind == 'flow':
            t_scale = kw.pop('t_scale', None)
            res = FR.train_flow(ro, K, f, c0, eps, tau_end, lr_c, t_scale=t_scale, max_secs=MAX_SECS, **kw)
        elif kind == 'gd':
            res = FR.train_gd(ro, K, f, c0, eps, tau_end, lr_c, max_secs=MAX_SECS, **kw)
        else:
            res = FR.train_adam(ro, K, f, c0, eps, max_secs=MAX_SECS, **kw)
        a = FR.analyse(res, K, f, c0, cs, As, ref=ref)
    except Exception as e:                                   # a failed run is a recorded result (never counted as a success), not a crash of the part
        nan = float('nan')
        return dict(n_formed=0, n_formed_sweep=[0, 0, 0, 0], gap=nan, dormant_disp=nan, dormant_wmax=nan, fit_rms=nan, exact=False, set_equal=False,
                    order_equal=False, tau_rel=nan, pos_err=nan, amp_err=nan, form_order=[], secs=0.0, status=('error: ' + type(e).__name__)[:12],
                    nfev=0, eps=float(eps), readout=ro.name, gap_k=nan, exact_topk=False, radau=0, aborted=0)
    a['n_formed_sweep'] = [a['n_formed_sweep'][s] for s in (1e-4, 1e-3, 1e-2, 1e-1)]
    a.update(secs=float(res['secs']), status=str(res['status'])[:12], nfev=int(res['nfev']), eps=float(eps), readout=ro.name,
             radau=int(res.get('radau_pieces', 0)), aborted=int(res.get('aborted_nfev', 0)))
    return a


def setup_nd(nd):
    """The design's dictionary with nd kernels of width SIG at spacing 2/(nd − 1) on [−1, 1] (`setup` has nd = ND = 21)."""
    L, K, _ = setup()
    return L, K, np.linspace(-1, 1, nd)


def cell_atoms(ro, K, f, c0, cs, As, eps, ref, tau_end, lr_c=None):
    """One run with the units grouped by atom and its analysis against the sliding path `ref` of the atoms' dictionary (centres c0). For the
    branch pairs (`FR.PaperPairs`) the atoms' centres are used twice (2·len(c0) units, unit j and j + len(c0) on atom j); any other readout has one
    unit per atom. *n_formed*: the units above the formed threshold; *n_atoms*: the atoms with at least one such unit; *set_equal*, *order_equal*,
    *tau_rel* (and *tau_ratio*, the entry times up to a common factor) are the atoms'; *exact*: K atoms formed, one unit each, every feature
    within the study's tolerances of a unit (K units formed, in K atoms); *gap_k*, *exact_topk*: threshold-free, over the units; *form_times*,
    *earliness*: the atoms' formation times and 1 − t_j/τ_j against the path's entry times, in the order of the path's entries (`FR.atom_summary`)."""
    lr_c = LR_C if lr_c is None else lr_c
    c0u = np.concatenate([FR.pts(c0)] * 2) if isinstance(ro, FR.PaperPairs) else FR.pts(c0)
    atom_of = ro.atoms(len(c0u)) if isinstance(ro, FR.PaperPairs) else np.arange(len(c0u))
    try:
        res = FR.train_flow(ro, K, f, c0u, eps, tau_end, lr_c, max_secs=MAX_SECS)
        a = FR.analyse(res, K, f, c0u, cs, As)                     # the units (no reference: the path's indices are the atoms')
        b = FR.atom_summary(res, As, atom_of, ref=ref)
    except Exception as e:                                       # a failed run is a recorded result (never counted as a success), not a crash
        nan = float('nan')
        return dict(n_formed=0, n_formed_sweep=[0, 0, 0, 0], n_atoms=0, atoms=[], one_each=False, gap=nan, dormant_disp=nan, dormant_wmax=nan,
                    fit_rms=nan, exact=False, set_equal=False, order_equal=False, tau_rel=nan, tau_ratio=nan, form_times=[], earliness=[],
                    pos_err=nan, amp_err=nan, form_order=[], secs=0.0, status=('error: ' + type(e).__name__)[:12], nfev=0, eps=float(eps),
                    readout=ro.name, gap_k=nan, exact_topk=False, radau=0, aborted=0)
    a['n_formed_sweep'] = [a['n_formed_sweep'][s] for s in (1e-4, 1e-3, 1e-2, 1e-1)]
    a.update(n_atoms=b['n_atoms'], atoms=b['atoms'], one_each=b['one_each'], form_order=b['form_order'], set_equal=b['set_equal'],
             order_equal=b['order_equal'], tau_rel=b['tau_rel'], tau_ratio=b['tau_ratio'], form_times=b['form_times'], earliness=b['earliness'],
             exact=bool(a['exact'] and b['n_atoms'] == len(As)),
             secs=float(res['secs']), status=str(res['status'])[:12], nfev=int(res['nfev']), eps=float(eps), readout=ro.name,
             radau=int(res.get('radau_pieces', 0)), aborted=int(res.get('aborted_nfev', 0)))
    return a


# --------------------------------------------------------------------------------------------------------------------------------
# A. the reformulation
# --------------------------------------------------------------------------------------------------------------------------------
def part_A():
    L, K, c0 = setup()
    cs, As = np.array(CSTAR), np.array(ASTAR)
    f = K.synth(cs[:, None], As)
    rng = np.random.default_rng(7)
    d = 9
    out = dict(identity=[], flows=[], conserved=[])
    ro = FR.PaperLaw(L)
    rhs = FR.make_rhs(ro, K, f, LR_C, d, 1)
    Kd = RT.Bumps(K.xs[:, 0], K.sig)
    repo = RT.moving_rhs(L, Kd, f, 'paradigm')
    for _ in range(5):
        al = rng.normal(0, 0.5, d); beta = np.exp(rng.normal(-1.0, 1.0, d)); c = np.sort(rng.uniform(-0.9, 0.9, d))
        z = np.r_[np.stack([al, np.sqrt(beta)], 1).ravel(), c]                    # the engine's state is unit-major: (α_i, ν_i)_i, c
        mine = rhs(0.0, z)
        theirs = repo(0.0, np.r_[al, np.log(beta), c])
        dal, dnu, dc = mine[:2 * d:2], mine[1:2 * d:2], mine[2 * d:]
        out['identity'].append(float(max(np.abs(dal - theirs[:d]).max() / np.abs(theirs[:d]).max(),
                                         np.abs(2 * dnu / np.sqrt(beta) - theirs[d:2 * d]).max() / np.abs(theirs[d:2 * d]).max(),
                                         np.abs(dc - theirs[2 * d:]).max() / np.abs(theirs[2 * d:]).max())))
    for eps in (1e-3, 1e-4, 1e-6):
        tau_end = 20.0
        raw = FR.train_flow(ro, K, f, c0, eps, tau_end, LR_C, n_t=401, rtol=1e-10)
        # the reduced flow: each kernel on its curve (round 37), `critical_rhs` with b0 = ε
        from scipy.integrate import solve_ivp
        Lam = RT.transit_scale(L, eps)
        taus = np.linspace(0, tau_end, 401)
        fr = RT.critical_rhs(L, Kd, f, 'paradigm', b0=eps * np.ones(ND))
        al_core = L.A * np.sqrt(eps / (2 * L.c))
        tol = np.r_[1e-12 * al_core * np.ones(ND), np.full(ND, 1e-14)]
        sol = solve_ivp(fr, (0, Lam * tau_end), np.r_[np.zeros(ND), c0], method='LSODA', rtol=1e-10, atol=tol, t_eval=Lam * taus)
        Wr = RT.critical_w(L, sol.y[:ND].T, eps)
        Cr = sol.y[ND:].T
        out['flows'].append(dict(eps=eps, dw=float(np.abs(raw['W'] - Wr).max() / np.abs(np.array(ASTAR)).max()),
                                 dc=float(np.abs(raw['C'][:, :, 0] - Cr).max() / SIG)))
    eps = 1e-4
    for name in ('P0', 'D2', 'D3'):
        ro_ = make_readout(name)
        d0 = ND
        z0 = np.r_[ro_.init(d0, eps).ravel(), c0]
        sc = ro_.scale(eps)
        from scipy.integrate import solve_ivp
        f_ = FR.make_rhs(ro_, K, f, LR_C, d0, 1)
        sol = solve_ivp(f_, (0, sc * 12.0), z0, method='LSODA', rtol=1e-11, atol=1e-16, t_eval=np.linspace(0, sc * 12.0, 61))
        P = sol.y[:d0 * ro_.q].reshape(d0, ro_.q, -1)
        if name == 'P0':
            Lq = law()
            q = P[:, 1] ** 2 - Lq.c * np.cosh(2 * P[:, 0] / Lq.A)
            drift = float(np.abs(q - q[:, :1]).max() / Lq.c)
            what = 'Q = β − v_s cosh(2α/a_max), relative to v_s'
        elif name == 'D2':
            q = P[:, 0] * P[:, 1]
            drift = float(np.abs(q - q[:, :1]).max() / eps)
            what = 'uv, relative to s²'
        else:
            q = 1 / P[:, 0] + 1 / P[:, 1]
            drift = float(np.abs(q - q[:, :1]).max() / np.abs(q[:, :1]).max())
            what = '1/u + 1/v, relative'
        out['conserved'].append(dict(readout=name, what=what, drift=drift))
    return out


# --------------------------------------------------------------------------------------------------------------------------------
# B. the design
# --------------------------------------------------------------------------------------------------------------------------------
def job_b(args):
    name, eps = args
    L, K, c0 = setup()
    cs, As = np.array(CSTAR), np.array(ASTAR)
    f = K.synth(cs[:, None], As)
    ref = reference(L, K, c0, f)
    ro = make_readout(name)
    kw = dict(t_scale=RT.transit_scale(L, eps)) if name == 'D1' else {}
    a = cell(ro, K, f, c0, cs, As, eps, ref, TAU_B, **kw)
    if name == 'D1':
        a['tau_rel'] = float('nan')                                       # no wall, no time scale to compare with the path's
    return dict(name=name, **a)


def part_B():
    from multiprocessing import Pool
    jobs = [(n, e) for n in ('P0', 'P0u', 'P4', 'D3', 'D4', 'D5', 'D2', 'D1') for e in EPS_B]
    with Pool(N_WORKERS) as pool:
        rows = pool.map(job_b, jobs, chunksize=1)
    L, K, c0 = setup()
    ref = reference(L, K, c0, K.synth(np.array(CSTAR)[:, None], np.array(ASTAR)))
    return dict(rows=rows, ref=ref)


# --------------------------------------------------------------------------------------------------------------------------------
# C. the benchmark targets
# --------------------------------------------------------------------------------------------------------------------------------
def sample_targets():
    """A seeded stratified sample of D‴.55's one-dimensional targets (families A, B, D of `exp_transit_certificate`'s raw)."""
    raw = E.load(CERT_RAW)
    rng = np.random.default_rng(SEED_C)
    out = []

    def take(rows, n, fam):
        rows = [r for r in rows if not r.get('skip')]
        for i in sorted(rng.choice(len(rows), size=min(n, len(rows)), replace=False)):
            r = rows[int(i)]
            out.append(dict(fam=fam, cs=np.asarray(r['cs'], float), As=np.asarray(r['As'], float), out=str(r['out']),
                            order=[int(j) for j in r['order']], e=[int(j) for j in r['e']]))
    A = list(raw['A'])
    take([r for r in A if r['out'] == 'nearest'], N_C['A_nearest'], 'A')
    take([r for r in A if r['out'] != 'nearest'], N_C['A_other'], 'A')
    take(list(raw['B']), N_C['B'], 'B')
    D = list(raw['D'])
    for tag, key in (('a', 'Da'), ('b', 'Db'), ('c', 'Dc')):
        take([r for r in D if r['tag'][0] == tag], N_C[key], 'D' + tag)
    return out


def job_c(args):
    i, tgt, cells = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    t0 = time.time()
    ref = reference(L, K, c0, f)
    out = dict(i=i, fam=tgt['fam'], path_out=tgt['out'], cs=tgt['cs'], As=tgt['As'], ref=ref, ref_secs=time.time() - t0,
               order_stored=tgt['order'], cells=[])
    tau_end = 1.05 * max(ref['tau']) + 1.0 if ref['tau'] else TAU_B
    for name, eps in cells:
        out['cells'].append(dict(name=name, **cell(make_readout(name), K, f, c0, tgt['cs'], tgt['As'], eps, ref, tau_end)))
    return out


def part_C():
    from multiprocessing import Pool
    tg = sample_targets()
    jobs = [(i, t, CELLS_C) for i, t in enumerate(tg)]
    with Pool(N_WORKERS) as pool:
        return pool.map(job_c, jobs, chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# D. depth 2's logarithmic scale
# --------------------------------------------------------------------------------------------------------------------------------
def job_d(args):
    i, tgt = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    ref = reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0
    ro = make_readout('D2')
    # rows_long: the first EPS_LONG_D at three times the horizon (is depth 2's failure a horizon too short for its later entries?)
    return dict(i=i, margin=ref['margin'], margin_never=ref['margin_never'],
                rows=[dict(**cell(ro, K, f, c0, tgt['cs'], tgt['As'], e, ref, tau_end)) for e in EPS_D],
                rows_long=[dict(**cell(ro, K, f, c0, tgt['cs'], tgt['As'], e, ref, 3 * tau_end)) for e in EPS_D[:2]])


def part_D():
    from multiprocessing import Pool
    tg = [dict(fam='design', cs=np.array(CSTAR), As=np.array(ASTAR), out='nearest', order=[], e=[])] + sample_targets()[:N_D]
    with Pool(N_WORKERS) as pool:
        return pool.map(job_d, list(enumerate(tg)), chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# E. the learning rates
# --------------------------------------------------------------------------------------------------------------------------------
def job_e(args):
    i, tgt = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    ref = reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0
    rows = []
    for g in GAMMAS_E:
        ro = FR.PaperLaw(law(gamma=g))
        rows.append(dict(kind='P0', param=g, ratio=g / (4 * LAW['G']), **cell(ro, K, f, c0, tgt['cs'], tgt['As'], EPS_E, ref, tau_end)))
    for lr in LRS_E:
        rows.append(dict(kind='D3', param=lr, ratio=lr, **cell(FR.Cubic3(lr), K, f, c0, tgt['cs'], tgt['As'], EPS_E, ref, tau_end)))
    return dict(i=i, rows=rows)


def part_E():
    from multiprocessing import Pool
    tg = [dict(cs=np.array(CSTAR), As=np.array(ASTAR))] + sample_targets()[:N_E]
    with Pool(N_WORKERS) as pool:
        return pool.map(job_e, list(enumerate(tg)), chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# F. fixed-step gradient descent
# --------------------------------------------------------------------------------------------------------------------------------
def job_f(args):
    i, tgt, name, eps, h = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    ref = reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0 if i else TAU_B
    # explicit Euler is stable only for h < 2/λ_max; at the exact fit the centres' curvature of a unit of amplitude A is lr_c A² ⟨(∂_cφ)²⟩
    dph = K.dphi(np.asarray(tgt['cs'], float).reshape(len(tgt['cs']), -1))
    lam = float(max(LR_C * A ** 2 * np.mean(dph[:, k, :] ** 2, axis=0).max() for k, A in enumerate(np.asarray(tgt['As'], float))))
    return dict(i=i, name=name, h=h, lam_max=lam, **cell(make_readout(name), K, f, c0, tgt['cs'], tgt['As'], eps, ref, tau_end, kind='gd', h=h))


def part_F():
    from multiprocessing import Pool
    design = dict(cs=np.array(CSTAR), As=np.array(ASTAR))
    tg = sample_targets()
    jobs = [(0, design, n, e, h) for n, e, h in CELLS_F]
    for i, t in enumerate(tg[:N_F]):
        jobs += [(i + 1, t, 'D3', 1e-6, 0.02), (i + 1, t, 'P0', 1e-4, 0.02)]       # h = 0.02 < 2/λ_max = 0.045 for |A| = 1
    with Pool(N_WORKERS) as pool:
        return pool.map(job_f, jobs, chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# G. random starts
# --------------------------------------------------------------------------------------------------------------------------------
def job_g(args):
    i, tgt = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    rows = []
    b = np.exp(RAND_G * np.random.default_rng(SEED_G + i).standard_normal(ND))        # the masses the readouts draw below
    ref_eq = reference(L, K, c0, f)
    ref_b = reference(L, K, c0, f, rho=b ** -0.5)
    tau_end = 1.05 * max(max(ref_eq['tau']), max(ref_b['tau'])) + 1.0
    for name in ('P0', 'D3'):
        rows.append(dict(kind='masses', name=name, ref_same=ref_b['set'] == ref_eq['set'] and ref_b['order'] == ref_eq['order'],
                         **cell(make_readout(name), K, f, c0, tgt['cs'], tgt['As'], EPS_G, ref_b, tau_end,
                                rng=np.random.default_rng(SEED_G + i), rand=RAND_G)))
    for a0 in AMP0_G:
        ro = FR.PaperLaw(law(), amp0=a0)
        rows.append(dict(kind='amp', name='P0', amp0=a0,
                         **cell(ro, K, f, c0, tgt['cs'], tgt['As'], EPS_G, ref_eq, tau_end, rng=np.random.default_rng(SEED_G + i))))
    return dict(i=i, rows=rows)


def part_G():
    from multiprocessing import Pool
    tg = [dict(cs=np.array(CSTAR), As=np.array(ASTAR))] + sample_targets()[:N_G]
    with Pool(N_WORKERS) as pool:
        return pool.map(job_g, list(enumerate(tg)), chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# H. Adam
# --------------------------------------------------------------------------------------------------------------------------------
def job_h(args):
    name, lr = args
    L, K, c0 = setup()
    cs, As = np.array(CSTAR), np.array(ASTAR)
    f = K.synth(cs[:, None], As)
    ref = reference(L, K, c0, f)
    return dict(name=name, lr=lr, **cell(make_readout(name), K, f, c0, cs, As, EPS_H, ref, 1.0, kind='adam', n_steps=STEPS_H, lr_adam=lr))


def part_H():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(job_h, [(n, lr) for n in ('P0', 'D3', 'D2') for lr in LRS_H], chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# I. the centres' rate
# --------------------------------------------------------------------------------------------------------------------------------
def job_i(args):
    i, tgt = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    ref = reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0
    rows = []
    for name in ('P0', 'D3', 'D2'):
        for m in LRC_I:
            rows.append(dict(name=name, mult=m, **cell(make_readout(name), K, f, c0, tgt['cs'], tgt['As'], EPS_I, ref, tau_end, lr_c=m * LR_C)))
    return dict(i=i, rows=rows)


def part_I():
    from multiprocessing import Pool
    tg = [dict(cs=np.array(CSTAR), As=np.array(ASTAR))] + sample_targets()[:N_I]
    with Pool(N_WORKERS) as pool:
        return pool.map(job_i, list(enumerate(tg)), chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# J. the plane
# --------------------------------------------------------------------------------------------------------------------------------
def job_j(args):
    i, cs, As = args
    import exp_transit_nd as TN
    warnings.filterwarnings('ignore')
    L = law()
    K, c0 = TN.full()
    f = K.synth(cs, As)
    R = RN.sliding_iss_nd(L, K, f, c0, np.ones(len(c0)))
    ent = [(float(t), int(j)) for t, j, k, c in R['events'] if k == 'enter']
    ref = dict(order=[j for t, j in ent], tau=[t for t, j in ent], set=sorted(int(j) for j in np.flatnonzero(R['act'])),
               exits=int(sum(k == 'exit' for t, j, k, c in R['events'])), rms=float(R['hist'][-1][3]))
    ref['margin'], ref['margin_never'] = margins(K, f, R, len(c0))
    tau_end = 1.05 * max(ref['tau']) + 1.0
    cells = [dict(name=nm, **cell(make_readout(nm), K, f, c0, cs, As, eps, ref, tau_end)) for nm, eps in CELLS_J]
    return dict(i=i, cs=cs, As=As, ref=ref, cells=cells)


def part_J():
    from multiprocessing import Pool
    raw = E.load(CERT_RAW)
    rows = [r for r in raw['C'] if not r.get('skip')]
    rng = np.random.default_rng(SEED_J)
    pick = sorted(int(i) for i in rng.choice(len(rows), size=N_J, replace=False))
    jobs = [(k, np.asarray(rows[i]['cs'], float), np.asarray(rows[i]['As'], float)) for k, i in enumerate(pick)]
    with Pool(N_WORKERS) as pool:
        return pool.map(job_j, jobs, chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# K. branch pairs
# --------------------------------------------------------------------------------------------------------------------------------
def job_k(args):
    """The pairs (a₀ ∈ A0_K, ε ∈ EPS_K + EPS_K_FAR) and the controls of one target of the main sample: P0, D3 at EPS_K, D2 also at EPS_K_FAR."""
    i, tgt = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    ref = reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0
    rows = []
    thr = FR.formed_threshold(tgt['As'])
    for a0 in A0_K:
        for eps in EPS_K + EPS_K_FAR:
            ro = FR.PaperPairs(law(), a0)
            rows.append(dict(start=f"pairs a₀ = {a0:g}", a0=a0, tau_pred=1.0 - ro.form_ratio(eps, thr),
                             **cell_atoms(ro, K, f, c0, tgt['cs'], tgt['As'], eps, ref, tau_end)))
    for nm in NAMES_K:
        for eps in EPS_K + (EPS_K_FAR if nm == 'D2' else ()):
            rows.append(dict(start=nm + (', α(0) = 0' if nm == 'P0' else ''), a0=0.0, tau_pred=float('nan'),
                             **cell_atoms(make_readout(nm), K, f, c0, tgt['cs'], tgt['As'], eps, ref, tau_end)))
    cs = np.sort(np.asarray(tgt['cs'], float))
    return dict(i=i, kind='std', margin=ref['margin'], sep=float(np.min(np.diff(cs)) / SIG) if cs.size > 1 else float('nan'), rows=rows)


def coupled_targets():
    """(index, target) of C's sample for the coupled test: the K = 2 targets whose features are at most 2σ apart (family D(b): 1.5σ and 2σ),
    eight with amplitudes of the same sign and eight of opposite signs."""
    return [(i, t) for i, t in enumerate(sample_targets()) if t['cs'].size == 2 and abs(t['cs'][1] - t['cs'][0]) <= SEP_KC + 1e-9]


def job_kc(args):
    """One coupled target: the pairs a₀ = A0_KC at EPS_KC, with the path's entries, the features' separation and sign, and the initial slopes of
    the entrants' duals (u_j = τ⟨f*, φ_j⟩ until the first entry: the straight line to the wall)."""
    i, tgt = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    ref = reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0
    thr = FR.formed_threshold(tgt['As'])
    g0 = np.abs(K.corr(c0[:, None], f)[0])
    ro = FR.PaperPairs(law(), A0_KC)
    rows = [dict(start=f"pairs a₀ = {A0_KC:g}", a0=A0_KC, tau_pred=1.0 - ro.form_ratio(eps, thr),
                 **cell_atoms(ro, K, f, c0, tgt['cs'], tgt['As'], eps, ref, tau_end)) for eps in EPS_KC]
    return dict(i=i, kind='coupled', margin=ref['margin'], same_sign=bool(tgt['As'][0] * tgt['As'][1] > 0),
                sep=float(abs(tgt['cs'][1] - tgt['cs'][0]) / SIG), path_order=[int(j) for j in ref['order']], path_tau=[float(t) for t in ref['tau']],
                slope=[float(g0[j]) for j in ref['order']], rows=rows)


def job_ka(args):
    """The scan of a₀ on the design: the pairs at A0_SCAN_K × EPS_SCAN_K (the controls P0, D3, D2 of the design are in the main sample's first job)."""
    i, tgt = args
    L, K, c0 = setup()
    f = K.synth(tgt['cs'][:, None], tgt['As'])
    ref = reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0
    thr = FR.formed_threshold(tgt['As'])
    rows = []
    for a0 in A0_SCAN_K:
        for eps in EPS_SCAN_K:
            ro = FR.PaperPairs(law(), a0)
            rows.append(dict(start=f"pairs a₀ = {a0:g}", a0=a0, tau_pred=1.0 - ro.form_ratio(eps, thr),
                             **cell_atoms(ro, K, f, c0, tgt['cs'], tgt['As'], eps, ref, tau_end)))
    return dict(i=i, kind='a0scan', margin=ref['margin'], rows=rows)


def run_k(job):
    kind, args = job
    return dict(std=job_k, coupled=job_kc, a0scan=job_ka)[kind](args)


def part_K():
    from multiprocessing import Pool
    design = dict(cs=np.array(CSTAR), As=np.array(ASTAR))
    tg = [design] + sample_targets()[:N_K]
    jobs = [('std', (i, t)) for i, t in enumerate(tg)] + [('coupled', (i, t)) for i, t in coupled_targets()] + [('a0scan', (0, design))]
    with Pool(N_WORKERS) as pool:
        return pool.map(run_k, jobs, chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# L. the dictionary's spacing
# --------------------------------------------------------------------------------------------------------------------------------
def on_grid(c0):
    """The design's features moved onto the nearest kernel centres of the dictionary c0 (the amplitudes unchanged): a control that removes the
    sliding, since the dictionary then contains the features."""
    c0 = np.asarray(c0, float)
    return np.array([c0[int(np.argmin(np.abs(c0 - c)))] for c in CSTAR])


def job_l(args):
    """One dictionary: variant 'off' (the design's features), 'on' (moved onto kernel centres: a control that removes the sliding) or 'scan' (D2 over
    ε = 10^(−k), k ∈ K_SCAN_L, on the design's features: where does depth 2 start to pass?)."""
    nd, variant = args
    L, K, c0 = setup_nd(nd)
    cs, As = (on_grid(c0) if variant == 'on' else np.array(CSTAR)), np.array(ASTAR)
    f = K.synth(cs[:, None], As)
    ref = reference(L, K, c0, f)
    tau_end = 1.05 * max(ref['tau']) + 1.0
    if variant == 'scan':
        rows = [dict(name='D2', **cell(make_readout('D2'), K, f, c0, cs, As, float(f'1e-{k}'), ref, tau_end)) for k in K_SCAN_L]
        return dict(nd=nd, variant=variant, spacing=float(2.0 / (nd - 1)), cs=cs, ref=ref, rows=rows)
    rows = [dict(name=nm, **cell(make_readout(nm), K, f, c0, cs, As, eps, ref, tau_end)) for nm in NAMES_L for eps in EPS_L]
    # rows_long: D2 at three times the horizon (is a depth-2 amplitude that has not converged at the horizon what makes a cell inexact?)
    rows_long = [dict(name='D2', **cell(make_readout('D2'), K, f, c0, cs, As, eps, ref, 3 * tau_end)) for eps in EPS_L]
    return dict(nd=nd, variant=variant, spacing=float(2.0 / (nd - 1)), cs=cs, ref=ref, rows=rows, rows_long=rows_long)


def part_L():
    from multiprocessing import Pool
    jobs = [(nd, 'off') for nd in ND_L] + [(nd, 'on') for nd in ND_ON_L] + [(nd, 'scan') for nd in ND_SCAN_L]
    with Pool(N_WORKERS) as pool:
        return pool.map(job_l, jobs, chunksize=1)


# --------------------------------------------------------------------------------------------------------------------------------
# the report
# --------------------------------------------------------------------------------------------------------------------------------
def g(x, fmt='.3g'):
    return 'nan' if x is None or (isinstance(x, float) and np.isnan(x)) else format(x, fmt)


def yn(b):
    return 'yes' if b else 'no'


def span(vals, fmt='.3g'):
    """'min – max' of the finite values (nan if there is none)."""
    v = [float(x) for x in vals if not np.isnan(x)]
    return f"{format(min(v), fmt)} – {format(max(v), fmt)}" if v else 'nan'


def row_cells(a):
    """The columns shared by the cell tables: formed (the sweep 1e-4…1e-1 in brackets), set, order, exact, gap, displacement, rms."""
    sw = '/'.join(str(v) for v in a['n_formed_sweep'])
    return [f"{a['n_formed']} ({sw})", yn(a['set_equal']), yn(a['order_equal']), yn(a['exact']), g(a['tau_rel']), g(a['gap'], '.1f'),
            g(a['dormant_disp']), g(a['dormant_wmax']), g(a['fit_rms']), g(a['secs'], '.1f'), a['status']]


HEAD = ['formed (1e-4/1e-3/1e-2/1e-1)', 'set = path', 'order = path', 'exact', 'τ rel. err.', 'gap (dec.)', 'dormant disp. (σ)',
        'dormant max\\|w\\|', 'rms', 's', 'status']


def table(head, rows):
    return ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)] + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows] + ['']


def pct(n, m):
    return f"{n} of {m}"


def report(A=None, B=None, C=None, D=None, Ee=None, F=None, G=None, H=None, I=None, J=None, Kb=None, Ls=None):
    o = ["# Which factorised readouts reproduce the sliding transit? (WP1)", "",
         "Generated by `experiments/exp_parametrization.py` from its raw; every number below is printed by the script. D‴.48's dictionary "
         f"(21 kernels, width {SIG}, {NS} samples of [−1, 1]), the centres' rate Gm = {LR_C:g} in every readout. A unit has *formed* when "
         "|w| > min(0.01 max\\|A*\\|, 0.1 min\\|A*\\|) (Paper 17's 0.01 max\\|A*\\|, lowered wherever a feature is more than 10 times smaller than the largest); *set*, *order* and *τ rel. err.* compare with the sliding path (`sliding_iss_nd`, the ε → 0 limit; τ = t/scale, "
         "scale the readout's wall); *gap* = log₁₀(smallest formed \\|w\\| / largest unformed \\|w\\|) at the end; *dormant disp.* the largest "
         "centre displacement of an unformed unit in kernel widths; *exact*: K units, one within 0.05σ and 5% of each feature.", ""]
    # the engine and the status of every cell: a cell that is not `ok` is excluded from the tables' denominators and never counted as a success
    tallies, radau = {}, {}

    def tally(part, cells):
        cs_ = [c for c in cells if isinstance(c, dict)]
        n0, b0 = tallies.get(part, (0, 0))
        tallies[part] = (n0 + len(cs_), b0 + sum(c.get('status', 'ok') != 'ok' for c in cs_))
        r0, m0, rec0 = radau.get(part, (0, 0, True))                       # cells with a piece redone by Radau; whether every cell records the field
        radau[part] = (r0 + sum(c.get('radau', 0) > 0 for c in cs_), m0 + len(cs_), rec0 and all('radau' in c for c in cs_))
    if B is not None:
        tally('B', B['rows'])
    for part, val in (('C', C), ('D', D), ('E', Ee), ('G', G), ('I', I), ('J', J), ('K', Kb), ('L', Ls)):
        for t in ([] if val is None else val):
            tally(part, t.get('cells', t.get('rows', [])) + t.get('rows_long', []))
    for part, val in (('F', F), ('H', H)):
        tally(part, [] if val is None else val)
    o += [f"Engine: gradient flows are integrated by LSODA (rtol 1e-9); a piece of a run is redone by Radau after {FR.STIFF_BUDGET} right-hand-side "
          f"evaluations or on a failed solve (a cell records the pieces redone and the abandoned evaluations); a wall-clock safety limit of {MAX_SECS} s "
          "per run is checked between the 8 pieces of a run. *order = path* compares the entry order of the units of the path's set (other units are "
          "ignored). Cells run and cells whose status is not `ok` (excluded from the denominators, never counted as a success), by part: "
          + '; '.join(f"{k}: {tallies[k][0]} / {tallies[k][1]}" for k in sorted(tallies)) + ". Cells with a piece redone by Radau (the raw's `radau` field), "
          "by part (the parts whose cells record it): " + '; '.join(f"{k}: {radau[k][0]} of {radau[k][1]}" for k in sorted(radau) if radau[k][2]) + ".", ""]
    if A is not None:
        o += ["## A. The reformulation", "",
              f"Largest relative difference between the law's right-hand side (`moving_rhs`) and the Euclidean gradient flow in (α, ν, c) with "
              f"learning rates (Gm, γm/4, Gm), over 5 random states: {g(max(A['identity']), '.2e')}.", ""]
        o += table(['ε', 'max \\|ΔW\\| / max\\|A*\\| (raw vs reduced)', 'max \\|ΔC\\| / σ'],
                   [[f"{r['eps']:.0e}", g(r['dw'], '.2e'), g(r['dc'], '.2e')] for r in A['flows']])
        o += table(['readout', 'conserved quantity', 'largest drift along the raw flow (ε = 1e-4, τ ≤ 12)'],
                   [[r['readout'], r['what'], g(r['drift'], '.2e')] for r in A['conserved']])
    if B is not None:
        ref = B['ref']
        o += ["## B. The design of D‴.48", "",
              f"Features at {CSTAR} with amplitudes {ASTAR}; the sliding path enters the kernels {ref['order']} at τ = "
              f"{', '.join(g(t, '.3f') for t in ref['tau'])} and ends with the set {ref['set']}.", ""]
        names = {'P0': 'P0, the paper law (lr_ν/lr_α = 0.025)', 'P0u': 'P0u, without saturation', 'P4': 'P4, equal learning rates (lr_ν = lr_α)',
                 'D3': 'D3, depth 3: w = u³ − v³', 'D4': 'D4, depth 4: w = u⁴ − v⁴', 'D5': 'D5, depth 5: w = u⁵ − v⁵', 'D2': 'D2, depth 2: w = u² − v²', 'D1': 'D1, depth 1 (the control; same time as P0 at each ε)'}
        for n in ('P0', 'P0u', 'P4', 'D3', 'D4', 'D5', 'D2', 'D1'):
            o += [f"**{names[n]}**", ""]
            o += table(['ε'] + HEAD, [[f"{r['eps']:.0e}"] + row_cells(r) for r in B['rows'] if r['name'] == n])
    if C is not None:
        o += ["## C. 200 targets of D‴.55's one-dimensional benchmark", "",
              "A seeded stratified sample (seed %d) of the families A (%d + %d), B (%d), D(a) (%d), D(b) (%d), D(c) (%d); the path of each "
              "target recomputed (`sliding_iss_nd`, polished)." % (SEED_C, N_C['A_nearest'], N_C['A_other'], N_C['B'], N_C['Da'], N_C['Db'],
                                                                   N_C['Dc']), ""]
        n = len(C)
        ref_flags = (f"reference paths: blocks relaxed on {sum(r['ref'].get('relaxed', False) for r in C)}, solver failures on "
                     f"{sum(r['ref'].get('solver_failed', 0) > 0 for r in C)}, out of budget on {sum(r['ref'].get('out_of_budget', False) for r in C)}"
                     if all('relaxed' in r['ref'] for r in C) else 'reference flags not stored')
        outs = {k: sum(r['path_out'] == k for r in C) for k in ('nearest', 'other', 'not')}
        o += [f"Targets: {n}; the path is exact with the nearest kernels on {outs['nearest']}, with other kernels on {outs['other']}, not exact on "
              f"{outs['not']}; its entry order equals the stored one on {sum(r['ref']['order'] == r['order_stored'] for r in C)} of {n}; "
              f"{sum(r['ref']['exits'] > 0 for r in C)} paths have an exit; {ref_flags}.", ""]
        rows = []
        for k, (name, eps) in enumerate(CELLS_C):
            cs_ = [r['cells'][k] for r in C]
            done = [c for c in cs_ if c['status'] == 'ok']
            tr = np.array([c['tau_rel'] for c in done], float)
            rows.append([name, f"{eps:.0e}", pct(len(done), n), pct(sum(c['set_equal'] for c in done), len(done)),
                         pct(sum(c['exact'] for c in done), len(done)), pct(sum(c['order_equal'] for c in done), len(done)),
                         g(np.nanmedian(tr)), g(np.nanpercentile(tr, 95) if len(tr) else float('nan')), g(np.median([c['gap'] for c in done if not np.isnan(c['gap'])]), '.1f'),
                         g(np.median([c['dormant_disp'] for c in done])), g(np.max([c['dormant_disp'] for c in done])),
                         g(np.median([c['secs'] for c in cs_]), '.1f')])
        o += table(['readout', 'ε', 'runs finished', 'set = path', 'exact', 'order = path', 'median τ err.', '95th pct τ err.', 'median gap',
                    'median dormant disp. (σ)', 'max dormant disp. (σ)', 'median s'], rows)
        o += ["τ err. is measured on a record of 1201 points of τ: values below about 10⁻³ are limited by that grid, and a target whose first entry is at "
              "τ < 0.1 has a grid floor above 1% (the maximum is therefore not reported).", ""]
        fams = sorted(set(r['fam'] for r in C))
        o += ["By family, the number of targets on which the formed set is the path's:", ""]
        rows = []
        for k, (name, eps) in enumerate(CELLS_C):
            rows.append([name, f"{eps:.0e}"] + [pct(sum(r['cells'][k]['set_equal'] for r in C if r['fam'] == f_ and r['cells'][k]['status'] == 'ok'),
                                                     sum(r['fam'] == f_ for r in C)) for f_ in fams])
        o += table(['readout', 'ε'] + fams, rows)
        bad_t = {}
        for r in C:
            fl = [c for c in r['cells'] if c['name'] in ('P0', 'D3') and not c['set_equal']]
            if fl:
                bad_t[r['i']] = fl
        o += [f"Targets on which P0 or D3 does not form the path's set in at least one cell: {len(bad_t)} of {n} "
              f"({sum(len(v) for v in bad_t.values())} cells).", ""]
        if all('margin' in r['ref'] for r in C):
            edges = (0.0, 0.002, 0.005, 0.01, 0.02, 0.05, float('inf'))
            c4 = [(k, nm, ep) for k, (nm, ep) in enumerate(CELLS_C) if nm in ('P0', 'D3')]
            rows = []
            for lo, hi in zip(edges[:-1], edges[1:]):
                sel = [r for r in C if lo <= r['ref']['margin'] < hi]
                rows.append([f"{lo:g} – {hi:g}" if hi < float('inf') else f"≥ {lo:g}", len(sel)]
                            + [f"{sum(not r['cells'][k]['set_equal'] for r in sel)} / {sum(not r['cells'][k]['exact'] for r in sel)}" for k, nm, ep in c4])
            rows.append(['all', n] + [f"{sum(not r['cells'][k]['set_equal'] for r in C)} / {sum(not r['cells'][k]['exact'] for r in C)}" for k, nm, ep in c4])
            o += ["The same against the path's **runner-up margin** μ — the smallest 1 − |u_j| over the kernels not yet active at an entry, later entrants "
                  "included (`path_duals_nd`); ε^½ is 10⁻³ at ε = 10⁻⁶ and 10⁻⁴ at 10⁻⁸. Each cell: the number of targets whose formed set is not the path's / "
                  "that are not exact (K units on the features).", ""]
            o += table(['μ', 'targets'] + [f"{nm} {ep:.0e}: set ≠ path / not exact" for k, nm, ep in c4], rows)
        rows = []
        for i, fl in sorted(bad_t.items()):
            r = C[i]
            rows.append([i, r['fam'], g(r['ref'].get('margin', float('nan'))),
                         ', '.join(f"{c['name']} {c['eps']:.0e}: {c['n_formed']} formed, {'exact' if c['exact'] else 'not exact'}" for c in fl)])
        o += table(['target', 'family', 'margin μ', "cells whose formed set differs from the path's"], rows)
    if D is not None:
        thr_d = FR.formed_threshold(ASTAR)
        c_d2 = ''
        if C is not None and ('D2', 1e-6) in CELLS_C:
            tr_c = np.array([c['tau_rel'] for c in (r['cells'][CELLS_C.index(('D2', 1e-6))] for r in C) if c['status'] == 'ok'], float)
            c_d2 = (f" (Part C, D2 at ε = 1e-06 on {len(C)} targets: median τ err. {g(np.nanmedian(tr_c))}, 95th percentile {g(np.nanpercentile(tr_c, 95))})")
        o += ["## D. Depth 2's convergence in ε", "",
              "D2 (w = u² − v², s² = ε) on the design (target 0) and the first targets of C's sample (all of family A, the easy family). Its wall is at "
              "η = ln(1/ε)/(4 lr). The formed threshold θ is crossed before the wall (at |w| ≈ 1), so the **first entrant** is early by the factor "
              "ln(1/θ)/ln(1/ε) to leading order: its *τ err.* ≈ ln(1/θ)/ln(1/ε), and *τ err. × ln(1/ε)* is a constant set by θ "
              f"(ln(1/{thr_d:g}) = {np.log(1 / thr_d):.2f} for the design's θ), not a fitted law. A later entrant is early by the same factor only while its "
              "dual is straight until its wall (separated features, as on these targets, whose median follows the law); where a dual bends before its wall "
              f"the earliness is larger and the entry times are not the path's up to one common factor{c_d2}; Part K's coupled targets show the mechanism "
              "for the pairs. What decides the outcome is the margin μ (below).", ""]
        rows = []
        for j, e in enumerate(EPS_D):
            rs = [t['rows'][j] for t in D]
            ok = [r for r in rs if r['status'] == 'ok']
            tr = np.array([r['tau_rel'] for r in ok], float)
            rows.append([f"{e:.0e}", pct(sum(r['set_equal'] for r in ok), len(ok)), pct(sum(r['exact'] for r in ok), len(ok)),
                         g(np.median([r['n_formed'] for r in ok])), g(np.median([r['gap'] for r in ok if not np.isnan(r['gap'])]), '.1f'),
                         g(np.nanmedian(tr)), g(np.nanmedian(tr) * np.log(1 / e)), g(np.median([r['dormant_disp'] for r in ok])),
                         g(rs[0]['n_formed']), g(rs[0]['tau_rel']), g(rs[0]['gap'], '.1f')])
        o += table(['ε', 'set = path', 'exact', 'median formed', 'median gap', 'median τ err.', 'median τ err. × ln(1/ε)',
                    'median dormant disp. (σ)', 'design: formed', 'design: τ err.', 'design: gap'], rows)
        if all('margin' in t for t in D):
            has_nev = all('margin_never' in t for t in D)
            rows = []
            for j, e in enumerate(EPS_D):
                rs = [t['rows'][j] for t in D]
                oki = [i for i, r in enumerate(rs) if r['status'] == 'ok']
                thr_mu = rule_threshold(e)
                agree = sum((D[i]['margin'] > thr_mu) == rs[i]['set_equal'] for i in oki)
                agree_nev = sum(rule_passes(D[i]['margin_never'], e) == rs[i]['set_equal'] for i in oki) if has_nev else None
                rows.append([f"{e:.0e}", pct(len(oki), len(rs)), g(thr_mu)] + ([pct(agree_nev, len(oki))] if has_nev else []) + [pct(agree, len(oki)),
                             pct(sum(rs[i]['exact_topk'] for i in oki), len(oki)), g(np.nanmedian([rs[i]['gap_k'] for i in oki]), '.2f')])
            mus = np.array([t['margin'] for t in D], float)
            mus_nev = np.array([rule_margin(t['margin_never']) for t in D], float) if has_nev else None
            o += [f"The rule behind these counts, with no fitted parameter: depth 2 forms the path's set when every kernel that never enters stays at a margin "
                  f"μ never > ln(1/θ)/ln(1/ε) from its wall (θ = {FR.THETA:g} of the largest amplitude; a coefficient at margin μ is ε^μ; a path on which "
                  "every kernel enters has no such kernel, and its μ never, NaN, is read as +∞). The same margin as in Part L: the kernels that never enter "
                  "are the ones whose walls the dynamics must not cross, a later entrant's wall is crossed on purpose. *rule agrees*: the targets on which "
                  "that prediction equals the outcome" + (" (*μ never*; *μ*: with the smallest margin of all kernels not yet active at an entry, the later "
                                                           "entrants included, as in the first version of this part)" if has_nev else "") + "; *top-K exact*: those on which "
                  "the K largest units are exact whatever the threshold; the last "
                  f"column is the threshold-free gap (K-th over (K+1)-th largest |w|) in decades. The margins of these {len(D)} targets: μ median "
                  f"{g(np.nanmedian(mus))}, maximum {g(np.nanmax(mus))}" + (f"; μ never median {g(np.median(mus_nev))}, maximum {g(np.max(mus_nev))}" if has_nev else "") + ".", ""]
            o += table(['ε', 'runs ok', 'rule: margin must exceed'] + (['rule agrees (μ never)'] if has_nev else []) + ['rule agrees' + (' (μ)' if has_nev else ''),
                        'top-K exact', 'median K-th/(K+1)-th gap'], rows)
        if all('rows_long' in t for t in D):
            rows = []
            for j, e in enumerate(EPS_D[:2]):
                rs = [t['rows_long'][j] for t in D]
                ok = [r for r in rs if r['status'] == 'ok']
                rows.append([f"{e:.0e}", pct(len(ok), len(rs)), pct(sum(r['set_equal'] for r in ok), len(ok)), pct(sum(r['exact'] for r in ok), len(ok)),
                             g(np.median([r['n_formed'] for r in ok])), g(np.median([r['fit_rms'] for r in ok])),
                             g(np.median([r['dormant_disp'] for r in ok]))])
            o += ['At three times the horizon (the entries of D2 come later or earlier than the path\'s by the τ error above):', '']
            o += table(['ε', 'runs finished', 'set = path', 'exact', 'median formed', 'median rms', 'median dormant disp. (σ)'], rows)
    if Ee is not None:
        o += ["## E. The learning rates", "",
              f"ε = {EPS_E:.0e}; P0 at lr_ν/lr_α = γ/(4G) (0.025 is the paper's), D3 at the learning rate lr of u and v (the centres' rate stays "
              f"{LR_C:g}); the design and the first {N_E} targets of C's sample.", ""]
        rows = []
        kinds = [(r['kind'], r['param']) for r in Ee[0]['rows']]
        for j, (kd, pa) in enumerate(kinds):
            rs = [t['rows'][j] for t in Ee]
            ok = [r for r in rs if r['status'] == 'ok']
            rows.append([kd, g(rs[0]['ratio']) if kd == 'P0' else g(pa), pct(sum(r['set_equal'] for r in ok), len(ok)),
                         pct(sum(r['exact'] for r in ok), len(ok)), g(np.nanmedian([r['tau_rel'] for r in ok])),
                         g(np.median([r['gap'] for r in ok if not np.isnan(r['gap'])]), '.1f'), g(np.median([r['dormant_disp'] for r in ok]))])
        o += table(['readout', 'lr_ν/lr_α (P0) or lr (D3)', 'set = path', 'exact', 'median τ err.', 'median gap', 'median dormant disp. (σ)'], rows)
    if F is not None:
        o += ["## F. Fixed-step gradient descent", "",
              "Explicit Euler on the same flow (every parameter by h·lr·(−∂F/∂p) per step); target 0 is the design.", ""]
        o += ["h·λ_max is the explicit Euler iteration's stability margin at the exact fit (λ_max = Gm·A²·⟨(∂_cφ)²⟩ for the largest amplitude A; stable "
              "only below 2): a run above 2 is an unstable iteration, not evidence about the flow.", ""]
        o += table(['target', 'readout', 'ε', 'h', 'h·λ_max', 'steps'] + HEAD,
                   [[r['i'], r['name'], f"{r['eps']:.0e}", g(r['h']), g(r['h'] * r['lam_max'], '.2f'), f"{r['nfev']:.2e}"] + row_cells(r) for r in F])
    if G is not None:
        o += ["## G. Random starts", "",
              f"ε = {EPS_G:.0e}. *masses*: ε b_i with b_i = exp({RAND_G} N(0, 1)); the path is computed with ρ_i = b_i^(−1/2). *amp*: the paper "
              "law with α_i(0) = amp0·N(0, 1) (the equal-mass path as reference; the limit is a sign-constrained path on a logarithmic scale, "
              "Paper 17 Proposition 3.2, so only the count and exactness are meaningful).", ""]
        rows = []
        keys = [(r['kind'], r['name'], r.get('amp0')) for r in G[0]['rows']]
        for j, (kd, nm, a0) in enumerate(keys):
            rs = [t['rows'][j] for t in G]
            ok = [r for r in rs if r['status'] == 'ok']
            rows.append([kd, nm, g(a0) if a0 else '—', pct(sum(r['set_equal'] for r in ok), len(ok)), pct(sum(r['exact'] for r in ok), len(ok)),
                         g(np.median([r['n_formed'] for r in ok])), g(np.median([r['gap'] for r in ok if not np.isnan(r['gap'])]), '.1f'),
                         g(np.median([r['dormant_disp'] for r in ok])), g(np.median([r['fit_rms'] for r in ok]))])
        o += table(['start', 'readout', 'amp0', 'set = path', 'exact', 'median formed', 'median gap', 'median dormant disp. (σ)', 'median rms'], rows)
    if H is not None:
        o += ["## H. Adam", "",
              f"One step size for every parameter and centre, {STEPS_H} steps, ε = {EPS_H:.0e}, the design. Reported, not claimed: the time scale "
              "is the step count, so *τ rel. err.* is not comparable and omitted.", ""]
        o += table(['readout', 'step size', 'formed (1e-4/1e-3/1e-2/1e-1)', 'exact', 'gap (dec.)', 'dormant disp. (σ)', 'rms', 's', 'status'],
                   [[r['name'], g(r['lr']), f"{r['n_formed']} ({'/'.join(str(v) for v in r['n_formed_sweep'])})", yn(r['exact']),
                     g(r['gap'], '.1f'), g(r['dormant_disp']), g(r['fit_rms']), g(r['secs'], '.1f'), r['status']] for r in H])
    if I is not None:
        o += ["## I. The centres' rate", "",
              f"ε = {EPS_I:.0e}; the centres' rate is the multiple shown of Gm = {LR_C:g} (0: fixed centres, the network is then a dictionary and "
              f"the limit the fixed-centre inverse scale space, which cannot fit an off-grid feature with one unit); the design and {N_I} targets "
              "of C's sample. *set = path* compares with the *sliding* path.", ""]
        rows = []
        keys = [(r['name'], r['mult']) for r in I[0]['rows']]
        for j, (nm, mu) in enumerate(keys):
            rs = [t['rows'][j] for t in I]
            ok = [r for r in rs if r['status'] == 'ok']
            rows.append([nm, g(mu), pct(sum(r['set_equal'] for r in ok), len(ok)), pct(sum(r['exact'] for r in ok), len(ok)),
                         g(np.median([r['n_formed'] for r in ok])), g(np.median([r['gap'] for r in ok if not np.isnan(r['gap'])]), '.1f'),
                         g(np.median([r['dormant_disp'] for r in ok])), g(np.median([r['fit_rms'] for r in ok])),
                         g(rs[0]['n_formed']), g(rs[0]['fit_rms'])])
        o += table(['readout', 'centres\' rate / Gm', 'set = path', 'exact', 'median formed', 'median gap', 'median dormant disp. (σ)', 'median rms',
                    'design: formed', 'design: rms'], rows)
    if J is not None:
        o += ["## J. The plane", "",
              f"{len(J)} of the 120 planar targets (seed {SEED_J}) on D‴.54's lattice (17 × 17 kernels of width 0.15, 41² samples of [−1, 1]²); "
              f"ε = {CELLS_J[0][1]:.0e}; the path of each target recomputed.", ""]
        rows = []
        for k, (nm, eps) in enumerate(CELLS_J):
            cs_ = [t['cells'][k] for t in J]
            ok = [c for c in cs_ if c['status'] == 'ok']
            rows.append([nm, f"{eps:.0e}", pct(len(ok), len(J)), pct(sum(c['set_equal'] for c in ok), len(ok)), pct(sum(c['exact'] for c in ok), len(ok)),
                         g(np.median([c['n_formed'] for c in ok])), g(np.median([c['gap'] for c in ok if not np.isnan(c['gap'])]), '.1f'),
                         g(np.median([c['dormant_disp'] for c in ok])), g(np.median([c['fit_rms'] for c in ok])),
                         g(np.median([c['secs'] for c in cs_]), '.1f')])
        o += table(['readout', 'ε', 'runs finished', 'set = path', 'exact', 'median formed', 'median gap', 'median dormant disp. (σ)', 'median rms',
                    'median s'], rows)
    if Kb is not None:
        Ks = [t for t in Kb if t.get('kind', 'std') == 'std']
        Kc = [t for t in Kb if t.get('kind') == 'coupled']
        Ka = [t for t in Kb if t.get('kind') == 'a0scan']
        pred = {(r['a0'], r['eps']): r['tau_pred'] for r in Ks[0]['rows'] if r['start'].startswith('pairs')}
        e_lo, e_hi = EPS_K[-1], EPS_K_FAR[-1]
        sep_min = min(t['sep'] for t in Ks) if all('sep' in t for t in Ks) else float('nan')
        o += ["## K. Branch pairs: the paper's logarithmic regime", "",
              "Paper 17's Proposition (starts on branches): a kernel started at α(0) ≠ 0 sits on a branch of its curve, its mass growing as "
              "ε·exp(γ a(0)(η − η(0))), with the weight ρ = 1/(γ\\|a(0)\\|) and the logarithmic scale L_ε = ln(1/ε), τ = t/L_ε; two kernels of opposite signs on the "
              "same atom give the two-sided box, the limit being the path of equal weights. **That proposition is proved for fixed centres. Its sliding "
              "extension (the centres move, the limit is the sliding path of the atoms), which is what is run here, is tested, not proved:** Paper 17's "
              "theorem on the sliding limit is proved for starts with α(0) = 0. Part G's random α_i(0) do not test it (their weights differ "
              "by decades); here every one of the 21 centres carries a *pair*, α(0) = +a₀ and −a₀, masses ε (42 units, `PaperPairs`), "
              "τ = t/scale with scale = ln(1/ε)/(γ a_max tanh(a₀/a_max)) = L_ε ρ, the time at which the path of unit weights has τ = 1, "
              "against the sliding path of the 21-atom dictionary. *formed*: units above the formed threshold (in brackets, the atoms that have at least "
              "one such unit); *set*, *order*, *τ rel. err.* are the atoms' (an atom forms when its first unit does); *exact*: K atoms formed, one unit "
              "each, every feature within 0.05σ and 5% of a unit; *τ ratio err.*: the entry times up to a common factor (t_j/t_first against the "
              "path's τ_j/τ_first); *gap K*: log₁₀ of the K-th over the (K+1)-th largest \\|w\\| over all units (threshold-free). *predicted*: "
              "the factor 1 − c/L_ε by which the **first entrant** of the path is early, derived from the branch law with no fitted parameter "
              "(`PaperPairs.form_ratio`: along the curve dη = dβ/(γaβ), so a kernel has \\|w\\| = θ after the mirror coordinate has grown by "
              "(1/γ)∫dβ/(aβ) from β = ε, against the wall's L_ε/(γ\\|a(0)\\|); a dual that is a straight line to its wall reaches the formed level at the "
              "fraction 1 − c/L_ε of its entry time; c is the constant of the column τ rel. err. × ln(1/ε)). A later atom is early by the same factor only while its dual is straight until its wall (separated "
              "features); an earlier entry changes its dual's slope and the formed level 1 − c/L_ε is then crossed on another part of the dual, so its "
              "earliness differs at the same order 1/L_ε. The tables on the design and on the "
              f"{len(Ks) - 1} targets of family A below (features at least {g(sep_min, '.2g')}σ apart) are where it looks exact; the coupled "
              "targets at the end show where it fails. The common factor is therefore claimed for the first entrant, and for a later atom only while its "
              "dual is straight. "
              f"ε = {e_lo:.0e} is pre-asymptotic for the pairs: L_ε = {np.log(1 / e_lo):.1f} and the early factor is {g(pred.get((0.1, e_lo)))} (a₀ = 0.1), "
              f"{g(pred.get((0.3, e_lo)))} (a₀ = 0.3) on the design, against {g(pred.get((0.1, e_hi)))} and {g(pred.get((0.3, e_hi)))} at ε = {e_hi:.0e} "
              f"(L_ε = {np.log(1 / e_hi):.0f}); the pairs are therefore run down to ε = {e_hi:.0e}, with D2 as the control there (P0 and D3 run on the power "
              f"scale ε^(−1/2), which this engine does not integrate in a practical time at those ε; they are run at ε = {', '.join(f'{e:.0e}' for e in EPS_K)} only). Controls: the same "
              "targets with P0 started at α(0) = 0, D3 and D2.", ""]
        names_k = [(r['start'], r['eps']) for r in Ks[0]['rows']]
        o += ["**The design** (target 0: part B's target).", ""]
        rows = []
        for j, (st, e) in enumerate(names_k):
            r = Ks[0]['rows'][j]
            rows.append([st, f"{e:.0e}", f"{r['n_formed']} ({r['n_atoms']})", yn(r['set_equal']), yn(r['order_equal']), yn(r['exact']), g(r['tau_rel']),
                         g(r['tau_pred']), g(r['tau_rel'] * np.log(1 / e)), g(r['tau_ratio']), g(r['gap'], '.2f'), g(r['gap_k'], '.2f'),
                         g(r['dormant_disp']), g(r['fit_rms']), r['status']])
        o += table(['start', 'ε', 'formed units (atoms)', 'set = path', 'order = path', 'exact', 'τ rel. err.', 'predicted τ rel. err.',
                    'τ rel. err. × ln(1/ε)', 'τ ratio err.', 'gap (dec.)', 'gap K (dec.)', 'dormant disp. (σ)', 'rms', 'status'], rows)
        o += [f"**The design and the first {len(Ks) - 1} targets of C's sample** (all of family A, the easy family; the features of every target are at least "
              f"{g(sep_min, '.2g')}σ apart): the number of targets on which each case holds, and medians; the last two columns are the median and the "
              "maximum over the targets of the largest displacement of an unformed unit (a statement about dormant units needs the maximum, not only the median).", ""]
        rows = []
        for j, (st, e) in enumerate(names_k):
            rs = [t['rows'][j] for t in Ks]
            ok = [r for r in rs if r['status'] == 'ok']
            rows.append([st, f"{e:.0e}", pct(len(ok), len(rs)), pct(sum(r['set_equal'] for r in ok), len(ok)), pct(sum(r['order_equal'] for r in ok), len(ok)),
                         pct(sum(r['exact'] for r in ok), len(ok)), g(np.median([r['n_formed'] for r in ok])), g(np.median([r['n_atoms'] for r in ok])),
                         g(np.median([r['gap_k'] for r in ok if not np.isnan(r['gap_k'])]), '.2f'), g(np.nanmedian([r['tau_rel'] for r in ok])),
                         g(np.nanmedian([r['tau_pred'] for r in ok])) if st.startswith('pairs') else '—',
                         g(np.nanmedian([r['tau_rel'] for r in ok]) * np.log(1 / e)), g(np.nanmedian([r['tau_ratio'] for r in ok])),
                         g(np.median([r['dormant_disp'] for r in ok])), g(np.max([r['dormant_disp'] for r in ok]))])
        o += table(['start', 'ε', 'runs finished', 'set = path', 'order = path', 'exact', 'median formed units', 'median atoms', 'median gap K (dec.)',
                    'median τ rel. err.', 'median predicted τ rel. err.', 'median τ rel. err. × ln(1/ε)', 'median τ ratio err.',
                    'median dormant disp. (σ)', 'max dormant disp. (σ)'], rows)
        # the ordering of the gap K by start
        chain = ('D2', 'pairs a₀ = 0.3', 'pairs a₀ = 0.1', 'D3')
        by_target = [{(r['start'], r['eps']): r['gap_k'] for r in t['rows']} for t in Ks]
        if all((c, e) in by_target[0] for c in chain for e in EPS_K):
            med = {k: np.median([bt[k] for bt in by_target]) for k in by_target[0]}
            links = list(zip(chain[:-1], chain[1:]))
            per_link = {e: [sum(bt[(a, e)] < bt[(b, e)] for bt in by_target) for a, b in links] for e in EPS_K}
            per_chain = {e: sum(all(bt[(a, e)] < bt[(b, e)] for a, b in links) for bt in by_target) for e in EPS_K}
            o += ["The ordering of the gap K by start, D2 < pairs a₀ = 0.3 < pairs a₀ = 0.1 < D3, holds for the medians at "
                  + ' and '.join(f"ε = {e:.0e}: {yn(all(med[(a, e)] < med[(b, e)] for a, b in links))}" for e in EPS_K)
                  + ". It is an ordering of medians: target by target, " + '; '.join(
                      f"at ε = {e:.0e} the three links hold on {', '.join(str(n) for n in per_link[e])} of {len(Ks)} targets and the whole chain on {per_chain[e]}"
                      for e in EPS_K) + ".", ""]
        rows = []
        for a0 in A0_K:
            for e in EPS_K + EPS_K_FAR:
                rs = [r for t in Ks for r in t['rows'] if r['start'] == f"pairs a₀ = {a0:g}" and r['eps'] == e]
                fin = [r for r in rs if r['status'] == 'ok' and not np.isnan(r['tau_rel'])]
                rows.append([f"pairs a₀ = {a0:g}", f"{e:.0e}", g(np.log(1 / e), '.1f'), g(np.median([r['tau_pred'] for r in rs])), pct(len(fin), len(rs)),
                             g(max(abs(r['tau_rel'] - r['tau_pred']) for r in fin)) if fin else 'nan', g(max(r['tau_ratio'] for r in fin)) if fin else 'nan'])
        o += ["The pairs' entry times against the derived common factor, by a₀ and ε (the runs in which every atom of the path's set formed; the same "
              f"{len(Ks)} targets, features at least {g(sep_min, '.2g')}σ apart):", ""]
        o += table(['start', 'ε', 'L_ε', 'predicted τ rel. err. (median)', 'every atom formed', 'largest \\|measured − predicted τ rel. err.\\|',
                    'largest τ ratio err.'], rows)
        orig = [r for t in Ks for r in t['rows'] if r['start'].startswith('pairs') and r['eps'] in EPS_K]
        dev = [r for r in orig if r['status'] == 'ok' and not np.isnan(r['tau_rel'])]
        o += [f"At the first two ε alone (ε = {', '.join(f'{e:.0e}' for e in EPS_K)}), the pair runs in which every atom of the path's set formed are "
              f"{len(dev)} of {len(orig)}; the largest difference between the measured and the predicted τ rel. err. among them is "
              f"{g(max(abs(r['tau_rel'] - r['tau_pred']) for r in dev) if dev else float('nan'))} and the largest τ ratio err. "
              f"{g(max(r['tau_ratio'] for r in dev) if dev else float('nan'))}: both figures belong to this sample (the design and family-A targets whose features are "
              f"at least {g(sep_min, '.2g')}σ apart), not to every target (see the coupled targets below).", ""]
        if Ka:
            sc = Ka[0]['rows']
            ctl = {(r['start'], r['eps']): r for r in Ks[0]['rows'] if not r['start'].startswith('pairs')}
            rows = []
            for e in EPS_SCAN_K:
                for r in [r for r in sc if r['eps'] == e] + [ctl[k] for k in ctl if k[1] == e]:
                    rows.append([r['start'], f"{e:.0e}", f"{r['n_formed']} ({r['n_atoms']})", yn(r['set_equal']), yn(r['order_equal']), yn(r['exact']),
                                 g(r['gap_k'], '.2f'), g(r['tau_rel']), g(r['tau_pred']), g(r['dormant_disp']), g(r['fit_rms']), r['status']])
            series = {e: [(r['a0'], r['gap_k']) for r in sc if r['eps'] == e] for e in EPS_SCAN_K}
            rises = {e: [(b[0], a[1], b[1]) for a, b in zip(series[e][:-1], series[e][1:]) if b[1] >= a[1]] for e in EPS_SCAN_K}
            steps = ' and '.join(f"at ε = {e:.0e} at {len(series[e]) - 1 - len(rises[e])} of {len(series[e]) - 1} steps"
                                 + (" (the exception: " if len(rises[e]) == 1 else " (the exceptions: " if rises[e] else '')
                                 + '; '.join(f"a₀ → {a0b:g}: {ga:.3f} → {gb:.3f}" for a0b, ga, gb in rises[e]) + (")" if rises[e] else '') for e in EPS_SCAN_K)
            o += ["**The pairs' a₀ on the design** (the scan; at each ε the controls of the design table follow the pairs): the threshold-free gap K falls with a₀ "
                  + steps + ", from the value of P0 (α(0) = 0) as a₀ → 0 to that of D2 for a₀ ≳ 1. a₀ is the interpolation parameter between P0 (an effective "
                  "depth-3 unit at small mass) and D2, and the pairs of a₀ = 0.1 and 0.3 lie between depth 2 and depth 3 because of their a₀, not by coincidence "
                  "of two values.", ""]
            o += table(['start', 'ε', 'formed units (atoms)', 'set = path', 'order = path', 'exact', 'gap K (dec.)', 'τ rel. err.', 'predicted τ rel. err.',
                        'dormant disp. (σ)', 'rms', 'status'], rows)
        if Kc:
            o += [f"**Coupled targets: where the common factor holds and where it fails.** The {len(Kc)} targets of C's sample whose two features are at "
                  f"most 2σ apart (family D(b), K = 2): {sum(t['same_sign'] for t in Kc)} with amplitudes of the same sign (A = (1, 1)) and "
                  f"{sum(not t['same_sign'] for t in Kc)} of opposite signs, the pairs a₀ = {A0_KC:g}. *Earliness* of an atom: 1 − t_j/τ_j, its formation time "
                  "over the path's entry time (both in units of the pair scale); the *first atom* is the path's first entrant, the *second atom* the other. "
                  "The prediction 1 − form_ratio is the first entrant's earliness; the second atom is early by it only if its dual is straight until its "
                  "wall.", ""]
            rows = []
            for name, same in (('opposite signs', False), ('the same sign', True)):
                ts = [t for t in Kc if t['same_sign'] == same]
                for e in EPS_KC:
                    rs = [next(r for r in t['rows'] if r['eps'] == e) for t in ts]
                    both = [r for r in rs if r['status'] == 'ok' and len(r['earliness']) == 2 and not np.any(np.isnan(r['earliness']))]
                    pr = np.array([r['tau_pred'] for r in rs])
                    f1 = np.array([r['earliness'][0] for r in both])
                    f2 = np.array([r['earliness'][1] for r in both])
                    pb = np.array([r['tau_pred'] for r in both])
                    rows.append([name, f"{e:.0e}", len(ts), pct(len(both), len(rs)), g(np.median(pr)), span(f1), span(f2),
                                 g(np.median(f2 / pb)) if len(both) else 'nan', span([r['tau_ratio'] for r in both]), pct(sum(r['set_equal'] for r in rs), len(rs)),
                                 pct(sum(r['exact'] for r in rs), len(rs))])
            o += table(['features', 'ε', 'targets', 'both atoms formed', 'predicted earliness', 'first atom: earliness (min – max)',
                        'second atom: earliness (min – max)', 'second atom: earliness / predicted (median)', 'τ ratio err. (min – max)', 'set = path', 'exact'], rows)
            e_lo1 = EPS_KC[0]
            rows = []
            for t in Kc:
                r = next(x for x in t['rows'] if x['eps'] == e_lo1)
                ft = r['form_times'] if len(r['form_times']) == 2 else [float('nan')] * 2
                eo = r['earliness'] if len(r['earliness']) == 2 else [float('nan')] * 2
                rows.append([t['i'], 'same' if t['same_sign'] else 'opposite', g(t['sep'], '.2g'), ' / '.join(g(x, '.2f') for x in t['path_tau']),
                             ' / '.join(g(x, '.2f') for x in ft), g(t['slope'][1] * t['path_tau'][0], '.2f'), g(t['slope'][1] * t['path_tau'][1], '.2f'),
                             g(r['tau_pred']), ' / '.join(g(x) for x in eo), g(r['tau_ratio']), yn(r['set_equal'])])
            o += [f"The same targets one by one at ε = {e_lo1:.0e}: the path's entry times and the atoms' formation times (τ); the second atom's dual: its value "
                  "at the first entry on the straight line u = τ⟨f*, φ⟩ that it follows until then, and the ratio of the path's entry time τ₂ to the time at which that "
                  "straight line would reach the wall (1: a dual that is straight to its wall, more: it has slowed down after the first entry); and the earliness "
                  "of the two atoms:", ""]
            o += table(['target', 'amplitudes', 'separation (σ)', 'path τ₁ / τ₂', 'formed at τ₁ / τ₂', 'second dual at the first entry',
                        'τ₂ / straight-line time to the wall', 'predicted earliness', 'earliness 1 / 2', 'τ ratio err.', 'set = path'], rows)
            allr = [(t, next(x for x in t['rows'] if x['eps'] == e)) for t in Kc for e in EPS_KC]
            ok_ = [(t, r) for t, r in allr if r['status'] == 'ok' and len(r['earliness']) == 2 and not np.any(np.isnan(r['earliness']))]
            d_first = max(abs(r['earliness'][0] - r['tau_pred']) for t, r in ok_)
            d_opp = max(abs(r['earliness'][1] - r['tau_pred']) for t, r in ok_ if not t['same_sign'])
            ratio_same = [r['earliness'][1] / r['tau_pred'] for t, r in ok_ if t['same_sign']]
            o += [f"On these {len(Kc)} targets and {len(EPS_KC)} values of ε the first entrant's earliness differs from the prediction by at most "
                  f"{g(d_first)} ({len(ok_)} of {len(allr)} runs with both atoms formed); the second atom of an opposite-sign pair by at most {g(d_opp)}; the "
                  f"second atom of a same-sign pair is early by {g(min(ratio_same))} to {g(max(ratio_same))} times the prediction. "
                  "So the factor 1 − c/L_ε is exact for the first entrant and for later atoms whose duals are straight, and not for the second atom of coupled "
                  "same-sign features; there \"formed\" is not \"entered\".", ""]
    if Ls is not None:
        Lm = [t for t in Ls if t['variant'] in ('off', 'on')]
        Lsc = [t for t in Ls if t['variant'] == 'scan']
        o += ["## L. The dictionary's spacing", "",
              f"Dictionaries of n kernels of width {SIG} at spacing 2/(n − 1) on [−1, 1] ({NS} samples), each readout against the sliding path of its own "
              "dictionary. *off-grid*: the design's features (the part B target); *on-grid*: the features moved onto the nearest kernel centres, at every "
              "spacing (neighbouring kernels overlap by exp(−spacing²/(4σ²)), which is above 0.5 below spacing 0.2). The on-grid arm is a control "
              "that removes the sliding — the dictionary contains the features, so no unit has to move onto one — and not another spacing sweep of the "
              "same problem; D2 passes there iff the same margin rule holds. The path's runner-up margins (`margins`): μ, the smallest "
              "1 − \\|u_j\\| of a kernel not yet active at an entry (the later entrants included), and μ never, the same over the kernels that never enter "
              "(NaN when every kernel enters: read as +∞, nothing is left to violate the rule). "
              "*rule* (D2 rows): depth 2 forms the path's set iff the margin exceeds ln(1/θ)/ln(1/ε) (Part D; the threshold is printed, then the rule's "
              "prediction of *set = path* from μ never); *gap K*: log₁₀ of the K-th over the (K+1)-th largest \\|w\\|; *pos err*, *amp err*: the largest "
              "distance of a feature from its unit (σ) and the largest relative amplitude error, when K units have formed.", ""]
        for variant, title in (('off', 'off-grid features (the design)'), ('on', 'on-grid features (the nearest kernel centres)')):
            blocks = [t for t in Lm if t['variant'] == variant]
            if not blocks:
                continue
            o += [f"**{title}**", ""]
            rows = []
            for t in blocks:
                for r in t['rows']:
                    rule = f"{g(rule_threshold(r['eps']))}: {yn(rule_passes(t['ref']['margin_never'], r['eps']))}" if r['readout'] == 'D2' else '—'
                    rows.append([t['nd'], g(t['spacing']), g(t['ref']['margin']), g(t['ref']['margin_never']), r['readout'], f"{r['eps']:.0e}",
                                 f"{r['n_formed']} ({'/'.join(str(v) for v in r['n_formed_sweep'])})", yn(r['set_equal']), yn(r['order_equal']), yn(r['exact']),
                                 g(r['pos_err']), g(r['amp_err']), g(r['gap'], '.2f'), g(r['gap_k'], '.2f'), g(r['dormant_disp']), g(r['fit_rms']), rule, r['status']])
            o += table(['kernels', 'spacing', 'μ', 'μ never', 'readout', 'ε', 'formed (1e-4/1e-3/1e-2/1e-1)', 'set = path', 'order = path', 'exact',
                        'pos err (σ)', 'amp err', 'gap (dec.)', 'gap K (dec.)', 'dormant disp. (σ)', 'rms', 'rule: threshold, prediction', 'status'], rows)
        if Lm and all('rows_long' in t for t in Lm):
            o += ["**D2 at three times the horizon** (the same cells with τ_end tripled):", ""]
            rows = []
            for t in Lm:
                for r in t['rows_long']:
                    rows.append([t['nd'], t['variant'] + '-grid', f"{r['eps']:.0e}", f"{r['n_formed']} ({'/'.join(str(v) for v in r['n_formed_sweep'])})", yn(r['set_equal']),
                                 yn(r['exact']), g(r['pos_err']), g(r['amp_err']), g(r['gap'], '.2f'), g(r['gap_k'], '.2f'), g(r['fit_rms']), r['status']])
            o += table(['kernels', 'features', 'ε', 'formed (1e-4/1e-3/1e-2/1e-1)', 'set = path', 'exact', 'pos err (σ)', 'amp err', 'gap (dec.)', 'gap K (dec.)',
                        'rms', 'status'], rows)
        d2 = [(t, r) for t in Lm for r in t['rows'] if r['readout'] == 'D2' and r['status'] == 'ok']
        o += [f"The depth-2 rule against *set = path* on the {len(d2)} D2 cells of the two tables ({len(Lm)} dictionaries × {len(EPS_L)} ε): with μ never "
              f"(NaN as +∞) it agrees on {sum(rule_passes(t['ref']['margin_never'], r['eps']) == r['set_equal'] for t, r in d2)}; with μ (the later entrants "
              f"included, as in Part D's first version) on {sum((t['ref']['margin'] > rule_threshold(r['eps'])) == r['set_equal'] for t, r in d2)}.", ""]
        fin = []
        if Lsc:
            o += ["**The rule out of sample: where does depth 2 start to pass?** D2 on the design's off-grid features over ε = 10^(−k), "
                  f"k = {K_SCAN_L[0]}, {K_SCAN_L[1]} … {K_SCAN_L[-1]}, on dictionaries of {', '.join(str(t['nd']) for t in Lsc)} kernels (those of {' and '.join(str(t['nd']) for t in Lsc if t['nd'] not in ND_L)} kernels are "
                  "not in the tables above). *ε from which D2 passes*: the largest ε such that D2 forms the path's set at it and at every smaller ε of the scan; "
                  "L_obs = ln(1/ε) there. The rule's value is L_rule = ln(1/θ)/μ never (the ε at which μ never = ln(1/θ)/ln(1/ε)); the rule has no fitted "
                  "parameter. *pattern*: one character per k, Y where D2 forms the path's set.", ""]
            rows, ratios = [], []
            for t in Lsc:
                ks = [int(round(-np.log10(r['eps']))) for r in t['rows']]
                okk = [bool(r['set_equal']) and r['status'] == 'ok' for r in t['rows']]
                kstar = next((k for j, k in enumerate(ks) if all(okk[j:])), None)
                m_nev = rule_margin(t['ref']['margin_never'])
                l_rule = float(np.log(1 / FR.THETA) / m_nev) if np.isfinite(m_nev) else 0.0
                l_obs = kstar * np.log(10) if kstar else float('nan')
                ratios.append(l_obs / l_rule if kstar and l_rule > 0 else float('nan'))
                rows.append([t['nd'], g(t['spacing']), g(t['ref']['margin_never']), g(l_rule, '.1f'), f"1e-{kstar}" if kstar else f"none to 1e-{ks[-1]}", g(l_obs, '.1f'),
                             g(ratios[-1], '.2f'), ''.join('Y' if x else '.' for x in okk)])
            o += table(['kernels', 'spacing', 'μ never', 'L_rule', 'ε from which D2 passes', 'L_obs', 'L_obs / L_rule', 'pattern'], rows)
            fin = [x for x in ratios if not np.isnan(x)]
            o += [f"L_obs/L_rule is between {g(min(fin), '.2f')} and {g(max(fin), '.2f')} on the {len(fin)} dictionaries: the rule's \"iff\" is right to about "
                  "that factor in L, at no fitted parameter." if fin else "", ""]
        # the rule's validation, all together (D's 21 targets × 6 ε, L's dictionaries, the scan)
        in_sample = []
        if D is not None and all('margin_never' in t for t in D):
            n_d = sum(1 for t in D for r in t['rows'] if r['status'] == 'ok')
            a_d = sum(rule_passes(t['margin_never'], r['eps']) == r['set_equal'] for t in D for r in t['rows'] if r['status'] == 'ok')
            in_sample.append(f"Part D's {len(D)} targets × {len(EPS_D)} ε ({a_d} of {n_d} D2 cells)")
        in_sample.append(f"the {len(Lm)} dictionaries of this part ({sum(rule_passes(t['ref']['margin_never'], r['eps']) == r['set_equal'] for t, r in d2)} of {len(d2)} D2 cells)")
        out_sample = ''
        if Lsc and fin:
            out_sample = (f", and out of sample on the scan above ({len(Lsc)} dictionaries × {len(K_SCAN_L)} values of ε, {sum(t['nd'] not in ND_L for t in Lsc)} of the dictionaries not in the tables: the transition at "
                          f"L_obs/L_rule between {g(min(fin), '.2f')} and {g(max(fin), '.2f')})")
        o += ["The margin rule (μ never > ln(1/θ)/ln(1/ε), NaN as +∞; no fitted parameter) is checked, in this report, on " + '; on '.join(in_sample) + out_sample + ".", ""]
    return '\n'.join(o) + '\n'


PARTS = 'ABCDEFGHIJKL'


def main():
    raw = os.path.join(OUT, 'parametrization_raw.npz')
    md = os.path.join(OUT, 'parametrization.md')
    old = E.load(raw)
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in PARTS]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(PARTS)
    res = dict(old)
    t0 = time.time()
    for k, fn in zip(PARTS, (part_A, part_B, part_C, part_D, part_E, part_F, part_G, part_H, part_I, part_J, part_K, part_L)):
        if k in parts:
            res[k] = fn()
            E.save(raw, **res)
            print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[copy.deepcopy(res.get(x)) for x in PARTS]))


if __name__ == '__main__':
    main()
