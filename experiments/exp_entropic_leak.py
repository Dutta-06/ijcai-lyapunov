"""
Proposition D‴.20: the entropic OT leak — a soft mass constraint with a critical excess
======================================================================================

The averaged slow loop of Theorem D‴ on the three-bump instance (σ₀ = 0.15σ_crit), with the entropic OT leak toward the matching prior
(three atoms at the exact representation, masses q = 1/3) at regularisations ε, against the hard OT leak (ε → 0) of Proposition D‴.12. The
entropic plan by Newton's method on the semi-dual (three columns), with ε-scaling and a minimum-norm fallback.

1. The critical excess (no dynamics). For the weights w_δ = q + δ(1.5e₁ − 0.5) (the excess of class 1 taken from both other atoms) and
   w′_δ = q + δ(e₁ − e₃) (taken from atom 3 alone): the weighted stationary point Θ*_ε(w) continued in ε from 256 down to 1 (four steps per
   octave; every Newton solve checked to 10⁻⁹), the separation block of the over-weighted class σ₀(I − (2/ε)Cov^G₁) − R^G₁ along the branch,
   ε_c by bisection inside the bracket of the sign change (each solve started from the nearer converged end), and at ε_c the conditional plan,
   the rank of Cov^G₁, λ_min(σ₀I − R^G_k), the collective block, the exact threshold 2σ₀λ_max(A^{−1/2}Cov^G₁A^{−1/2}) and the rank-one values
   2σ₀p(1 − p)Δ̂_jᵀA⁻¹Δ̂_j with all of the outflow p sent to atom j. Also δ_c(ε), the critical excess at fixed ε, by continuation in δ.
2. Dynamics. m = 15, 60, 240 against M = 510 from the warm clouds of `exp_ot_leak.py` (the same seeds), ε = 1, 4, 16, 64, 600 slow units
   (dt = 0.02): the cluster weights by nearest atom, the fit error, W²_{G⁻¹}(μ_m, ν) and the function-level distance to the reference, per seed.
3. Continuation of the ε = 4, 16, 64 runs from their stored s = 600 states (float32) to 6000 slow units with dt = 0.1: the count changes,
   the end counts, the final drift, and the margins of Θ*_ε(w) at the realised weights of s = 600 and of s = 6000 (stable iff every
   separation block and the collective block are positive).
4. The floors: Proposition D‴.5's function-level law at the entropic stationary point Θ*_ε(q) (ε = 16, 64) and at the nearest-atom one, against
   the observed floors per seed and the nearest-atom leak's floors on the same clouds (`ot_leak_raw.npz`).

Output: results/chaos/entropic_leak.md and entropic_leak_raw.npz (EL_S_END, EL_SEEDS, CLW_WORKERS; `--from-raw` regenerates the md; `--only
crit,cont,margins,sens` recomputes the listed parts and merges each into the raw file as soon as it is finished, keeping the 600-unit runs;
EL_CACHE=<dir> caches every finished continuation job there, so that an interrupted run resumes without repeating them).
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
from src.core import entropic_leak as EL
from src.core import ot_leak as OT

OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos')
NET_KW = dict(a_max=4.0, c_max=3.0, s_min=0.3, s_max=2.0)
G = np.array([1.0, 0.25, 0.25])
SQG = np.sqrt(G)
BUMPS = [(1.0, -1.5, 0.5), (-0.7, 0.0, 0.4), (0.8, 1.5, 0.5)]
PER = 2 * np.pi
N_Q = 256
SIGMA_CRIT = 0.047162567584665846
SIGMA0 = 0.15 * SIGMA_CRIT
Q = np.ones(3) / 3
EPSS = (1.0, 4.0, 16.0, 64.0)
MS, M_REF = (15, 60, 240), 510
S_END, DT_S, REC, PLAN_EVERY = float(os.environ.get('EL_S_END', 600.0)), 0.02, 50, 10
N_SEEDS = int(os.environ.get('EL_SEEDS', 4))
N_WORKERS = int(os.environ.get('CLW_WORKERS', 3))
JITTER = 0.02


def f_vec(X):
    return sum(A * np.exp(-(X - C) ** 2 / (2 * S ** 2)) for A, C, S in BUMPS)


def x_ref(t):
    return np.array([2.0 * np.sin(t)])


def make_net():
    return L.SaturatedLiftedNetwork(1, **NET_KW)


def exact_representation(net):
    return net.unconstrained(np.array([[3 * A] for A, C, S in BUMPS]), np.array([[C] for A, C, S in BUMPS]), np.array([S for A, C, S in BUMPS]))


def make_loop(net, eps):
    return EL.EntropicPriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact_representation(net), eps, n_quad=N_Q)


# ---------------------------------------------------------------------------- 1. the critical excess (continuation in ε, checked convergence)
FAMILIES = {'w_δ': lambda d: Q + d * (1.5 * np.eye(3)[0] - 0.5),             # the excess taken from both other atoms
            'w′_δ': lambda d: Q + d * (np.eye(3)[0] - np.eye(3)[2])}          # the excess taken from atom 3 alone
DELTAS = {'w_δ': (0.02, 0.05, 0.1, 0.2), 'w′_δ': (0.02, 0.05, 0.1)}
EPS_BRANCH = tuple(float(2.0 ** (8 - k / 4)) for k in range(33))             # 256 … 1, four steps per octave
RES_TOL = 1e-9


def solve(net, w, eps, start):
    """Θ*_ε(w) by Newton from `start`; raises if it does not converge (the D‴.20 table of round 14 was built on unconverged solves)."""
    lp = make_loop(net, eps)
    atoms, res = EL.weighted_stationary_atoms_entropic(lp, w, start, eps, Q)
    if res > RES_TOL:
        raise RuntimeError(f"Newton did not converge: w = {w}, ε = {eps}, residual {res:.1e}")
    return atoms, res, lp


def branch(net, w, eps_list=EPS_BRANCH, start=None):
    """Continuation of Θ*_ε(w) along eps_list (decreasing), from the exact representation at the first ε."""
    atoms = exact_representation(net) if start is None else start
    out = []
    for eps in eps_list:
        atoms, res, lp = solve(net, w, eps, atoms)
        out.append((eps, atoms.copy(), res))
    return out


def collective_min(lp, w, atoms, eps, h=1e-6):
    """λ_min of the whitened weighted Hessian of the K-atom entropic energy at Θ*_ε(w) (Proposition D‴.5(i)'s collective block)."""
    K = atoms.shape[0]; d = 3 * K
    v = atoms.ravel(); Jm = np.zeros((d, d))
    gf = lambda V: EL.weighted_energy_grad_entropic(lp, w, V.reshape(K, 3), eps, Q).ravel()
    for j in range(d):
        e = np.zeros(d); e[j] = h
        Jm[:, j] = (gf(v + e) - gf(v - e)) / (2 * h)
    sw = np.repeat(np.sqrt(w), 3); sg = np.tile(SQG, K)
    Hs = (sw * sg)[:, None] * Jm * (sg / sw)[None, :]
    return float(np.linalg.eigvalsh(0.5 * (Hs + Hs.T)).min())


def analyse(lp, w, atoms, eps, k=0):
    """At Θ*_ε(w): every class's separation minimum, λ_min of A_k = σ₀I − R^G_k, the collective minimum, the conditional plan of class k,
    the exact threshold 2σ₀λ_max(A^{−1/2}Cov^G_kA^{−1/2}) (when A ≻ 0) and the rank-one values with all of p sent to each other atom."""
    blocks, mins, cov = EL.separation_blocks_entropic(lp, w, atoms, eps, Q)
    a_, c_, s_, diff_, phi_ = lp.features(atoms); r = lp.F - phi_ @ (w * a_)
    Rw = SQG[None, :, None] * lp.residual_curvature(r, atoms) * SQG[None, None, :]
    A = SIGMA0 * np.eye(3)[None] - Rw
    Amin = np.array([np.linalg.eigvalsh(Ak).min() for Ak in A])
    pi, T, _ = EL.weighted_targets(lp.prior, w, atoms, Q, eps, G)
    cond = pi[k] / pi[k].sum()
    p = float(1.0 - cond[k])
    exact_thr, rank1 = float('nan'), {}
    if Amin[k] > 0:
        ev, U = np.linalg.eigh(A[k]); Ah = U @ np.diag(ev ** -0.5) @ U.T
        exact_thr = float(2 * SIGMA0 * np.linalg.eigvalsh(Ah @ cov[k] @ Ah).max())
        for j in range(3):
            if j != k:
                dh = (lp.prior[j] - lp.prior[k]) / SQG
                rank1[j] = float(2 * SIGMA0 * p * (1 - p) * dh @ np.linalg.solve(A[k], dh))
    return dict(sep=mins, Amin=Amin, coll=collective_min(lp, w, atoms, eps), cond=cond, p=p, exact_thr=exact_thr, rank1=rank1,
                cov_ev=np.linalg.eigvalsh(cov[k]))


def critical_row(fam, delta):
    net = make_net()
    w = FAMILIES[fam](delta)
    br = branch(net, w)
    margins = [EL.separation_blocks_entropic(make_loop(net, e), w, at, e, Q)[1][0] for e, at, r in br]
    grid = [(e, mg, r) for (e, at, r), mg in zip(br, margins) if e in (256.0, 64.0, 16.0, 4.0, 1.0)]
    eps_c, info = float('nan'), None
    signs = [i for i in range(len(br) - 1) if margins[i] > 0 >= margins[i + 1]]
    if signs:
        i = signs[0]
        a, b = np.log(br[i][0]), np.log(br[i + 1][0])                         # margin > 0 at a (larger ε), ≤ 0 at b
        at_a, at_b = br[i][1], br[i + 1][1]
        for _ in range(40):
            c = 0.5 * (a + b)
            start = at_a if abs(c - a) <= abs(c - b) else at_b
            atoms, res, lp = solve(net, w, float(np.exp(c)), start)
            mg = EL.separation_blocks_entropic(lp, w, atoms, float(np.exp(c)), Q)[1][0]
            if mg > 0:
                a, at_a = c, atoms
            else:
                b, at_b = c, atoms
        eps_c = float(np.exp(0.5 * (a + b)))
        atoms, res, lp = solve(net, w, eps_c, at_a)
        info = analyse(lp, w, atoms, eps_c)
        info['res'] = res
    return dict(family=fam, delta=delta, w=w, grid=grid, n_sign_changes=len(signs), eps_c=eps_c, info=info,
                branch_res=max(r for e, at, r in br))


def critical_excess_at(fam, eps, d_step=0.0025, d_max=0.3):
    """δ_c(ε): continuation in δ at fixed ε from δ = 0 (Θ*_ε(q) from the ε-branch at w = q) to the first sign change of the over-weighted
    class's separation margin, then bisection; returns δ_c, the conditional fraction p_c there and ε/(2‖Δ̂‖²) with Δ̂ to the atom of the
    larger share."""
    net = make_net()
    br_q = branch(net, Q, eps_list=[e for e in EPS_BRANCH if e >= eps] + ([eps] if eps not in EPS_BRANCH else []))
    atoms = br_q[-1][1]
    prev_d, prev_at, prev_mg = 0.0, atoms, None
    d = 0.0
    while d < d_max:
        d = round(d + d_step, 10)
        w = FAMILIES[fam](d)
        atoms, res, lp = solve(net, w, eps, prev_at)
        mg = EL.separation_blocks_entropic(lp, w, atoms, eps, Q)[1][0]
        if mg <= 0:
            a, b, at_a = prev_d, d, prev_at
            for _ in range(30):
                c = 0.5 * (a + b)
                at, rr, lp = solve(net, FAMILIES[fam](c), eps, at_a)
                if EL.separation_blocks_entropic(lp, FAMILIES[fam](c), at, eps, Q)[1][0] > 0:
                    a, at_a = c, at
                else:
                    b = c
            dc = 0.5 * (a + b)
            at, rr, lp = solve(net, FAMILIES[fam](dc), eps, at_a)
            info = analyse(lp, FAMILIES[fam](dc), at, eps)
            j = 1 + int(np.argmax(info['cond'][1:]))
            dh = (lp.prior[j] - lp.prior[0]) / SQG
            return dict(family=fam, eps=eps, delta_c=dc, p_c=info['p'], scale=eps / (2 * float(dh @ dh)), res=rr)
        prev_d, prev_at, prev_mg = d, atoms, mg
    return dict(family=fam, eps=eps, delta_c=float('nan'), p_c=float('nan'), scale=float('nan'), res=float('nan'))


# ---------------------------------------------------------------------------- 2. dynamics
def job(args):
    seed, eps = args
    net = make_net(); exact = exact_representation(net)
    rng = np.random.default_rng(seed)
    Th_all = CL.mixture_start(exact, M_REF, rng, JITTER, G)                 # the clouds of exp_ot_leak.py (same seeds, same generator calls)
    sysm = {m: Th_all[:m].copy() for m in MS + (M_REF,)}
    loops = {m: make_loop(net, eps) for m in sysm}
    T = {}
    n = int(round(S_END / DT_S))
    rec = dict(s=[], w={m: [] for m in sysm}, fit={m: [] for m in sysm}, W2nu={m: [] for m in sysm}, h2={m: [] for m in MS})
    t0 = time.time()
    for k in range(n):
        if k % PLAN_EVERY == 0:
            for m, Th in sysm.items():
                T[m] = loops[m].leak_target(Th)
        new = {}
        for m, Th in sysm.items():
            lp = loops[m]
            k1 = lp.drift(Th, T[m]); k2 = lp.drift(Th + 0.5 * DT_S * k1, T[m]); new[m] = Th + DT_S * k2
        sysm = new
        if k % REC == 0 or k == n - 1:
            rec['s'].append((k + 1) * DT_S)
            fM = loops[M_REF].f_hat(sysm[M_REF])
            for m, Th in sysm.items():
                rec['w'][m].append(CL.cluster_weights(Th, exact, G)[1])
                rec['fit'][m].append(loops[m].fit_error(Th))
                rec['W2nu'][m].append(_w2nu(Th, exact))
            for m in MS:
                rec['h2'][m].append(float(np.mean((loops[m].f_hat(sysm[m]) - fM) ** 2)))
    out = dict(seed=seed, eps=eps, s=np.array(rec['s']))
    for key in ('w', 'fit', 'W2nu', 'h2'):
        out[key] = {m: np.array(v) for m, v in rec[key].items()}
    out['final'] = {m: Th.astype(np.float32) for m, Th in sysm.items()}
    print(f"seed {seed}, ε = {eps}: {time.time() - t0:.0f}s", flush=True)
    return out


def _w2nu(Th, prior):
    """W²_{G⁻¹}(μ_Θ, ν) for ν uniform on the three prior atoms (the hard transport cost, a diagnostic for every leak)."""
    m = Th.shape[0]
    C = OT._cost_matrix(Th, prior, G)
    if m % 3 == 0:
        rep = np.repeat(np.arange(3), m // 3)
        from scipy.optimize import linear_sum_assignment
        r, c = linear_sum_assignment(C[:, rep]); return float(C[r, rep[c]].mean())
    gamma = OT.transport_plan(C, np.full(m, 1.0 / m), Q); return float(np.sum(gamma * C))


# ---------------------------------------------------------------------------- 3. continuation of the runs at ε = 4, 16, 64
EPS_CONT, S_CONT, DT_CONT, PLAN_CONT, REC_CONT = (4.0, 16.0, 64.0), 6000.0, 0.1, 2, 500     # records every 50 slow units


def drift_norm(lp, Th, T):
    d = lp.drift(Th, T)
    return float(np.sqrt(np.sum(d ** 2 / G[None, :], axis=1)).max())


def job_continue(args):
    """Continue the stored s = 600 states (float32) of one (seed, ε) job to S_CONT with dt = 0.1 (the midpoint rule, the target every
    2 steps): records of the counts by nearest atom, the fit error, W²(μ, ν) and the function-level distance to the continued reference; at the
    end the states, the largest whitened drift and the within-class spreads."""
    seed, eps, finals = args
    t0 = time.time()
    net = make_net(); exact = exact_representation(net)
    sysm = {m: np.asarray(finals[m], float) for m in MS + (M_REF,)}
    loops = {m: make_loop(net, eps) for m in sysm}
    T = {}
    n = int(round((S_CONT - S_END) / DT_CONT))
    rec = dict(s=[], counts={m: [] for m in sysm}, fit={m: [] for m in sysm}, W2nu={m: [] for m in sysm}, h2={m: [] for m in MS})
    for k in range(n + 1):
        if k % PLAN_CONT == 0:
            for m, Th in sysm.items():
                T[m] = loops[m].leak_target(Th)
        if k % REC_CONT == 0:
            rec['s'].append(S_END + k * DT_CONT)
            fM = loops[M_REF].f_hat(sysm[M_REF])
            for m, Th in sysm.items():
                rec['counts'][m].append(np.bincount(np.argmin(OT._cost_matrix(Th, exact, G), axis=1), minlength=3))
                rec['fit'][m].append(loops[m].fit_error(Th)); rec['W2nu'][m].append(_w2nu(Th, exact))
            for m in MS:
                rec['h2'][m].append(float(np.mean((loops[m].f_hat(sysm[m]) - fM) ** 2)))
        if k == n:
            break
        for m, Th in list(sysm.items()):
            lp = loops[m]
            k1 = lp.drift(Th, T[m]); sysm[m] = Th + DT_CONT * lp.drift(Th + 0.5 * DT_CONT * k1, T[m])
    out = dict(seed=seed, eps=eps, s=np.array(rec['s']))
    for key in ('counts', 'fit', 'W2nu', 'h2'):
        out[key] = {m: np.array(v) for m, v in rec[key].items()}
    out['final'] = {m: Th.astype(np.float32) for m, Th in sysm.items()}
    out['drift'] = {m: drift_norm(loops[m], Th, loops[m].leak_target(Th)) for m, Th in sysm.items()}
    print(f"[continue] seed {seed}, ε = {eps}: {time.time() - t0:.0f}s", flush=True)
    return out


# ---------------------------------------------------------------------------- 4. the margins at the realised weights
def end_margins(eps, counts_list):
    """For each count vector: Θ*_ε(w) by continuation in the weights from Θ*_ε(q) (itself from the ε-branch at w = q), its residual, every
    class's separation minimum and the collective minimum (stable iff all positive)."""
    net = make_net()
    br_q = branch(net, Q, eps_list=[e for e in EPS_BRANCH if e >= eps] + ([eps] if eps not in EPS_BRANCH else []))
    atoms_q = br_q[-1][1]
    out = {}
    for counts in counts_list:
        key = tuple(int(x) for x in counts)
        if key in out:
            continue
        w = np.asarray(counts, float) / np.sum(counts)
        at = atoms_q
        try:
            for j in range(1, 9):
                at, res, lp = solve(net, Q + (w - Q) * j / 8, eps, at)
            blocks, mins, cov = EL.separation_blocks_entropic(lp, w, at, eps, Q)
            out[key] = dict(res=res, sep=mins, coll=collective_min(lp, w, at, eps), excess=float(np.abs(w - Q).max()))
        except RuntimeError:
            out[key] = dict(res=float('nan'), sep=np.full(3, np.nan), coll=float('nan'), excess=float(np.abs(w - Q).max()))
    return out


# ---------------------------------------------------------------------------- 5. the weight-sensitivity law at the entropic stationary point
def sensitivity_law(eps, h=1e-4):
    """Proposition D‴.5's function-level law at Θ*_ε(q): the derivatives of f̂_w = Σ_kw_kΦ̃(·;ϑ*_k(w)) along the zero-sum directions
    e₁ − e₃, e₂ − e₃ (central differences of the entropic weighted stationary point; eps = 0: the nearest-atom leak toward the matching prior),
    and the predicted E‖f̂_m − f̂_M‖² = (1/m − 1/M)tr(F₂C₂F₂ᵀ)/N for m = 15, 60, 240 against M = 510, C₂ the multinomial covariance of (w₁, w₂)."""
    net = make_net(); exact = exact_representation(net)
    if eps == 0:
        from src.core import prior_design as PD
        lp = PD.MixturePriorLoop(net, f_vec, x_ref, PER, G, SIGMA0, exact, n_quad=N_Q)
        stat = lambda w, start: PD.weighted_stationary_atoms_prior(lp, w, start)
        atoms_q, _ = stat(Q, exact)
    else:
        lp = make_loop(net, eps)
        br_q = branch(net, Q, eps_list=[e for e in EPS_BRANCH if e >= eps])
        atoms_q = br_q[-1][1]
        stat = lambda w, start: EL.weighted_stationary_atoms_entropic(lp, w, start, eps, Q)
    fh = lambda w, at: lp.features(at)[4] @ (w * lp.features(at)[0])
    F2 = []
    for j in (0, 1):
        v = np.zeros(3); v[j] = 1.0; v[2] = -1.0
        ap, rp = stat(Q + h * v, atoms_q); am, rm = stat(Q - h * v, atoms_q)
        assert rp < RES_TOL and rm < RES_TOL
        F2.append((fh(Q + h * v, ap) - fh(Q - h * v, am)) / (2 * h))
    F2 = np.array(F2).T                                                        # (N, 2)
    C2 = np.array([[Q[0] * (1 - Q[0]), -Q[0] * Q[1]], [-Q[0] * Q[1], Q[1] * (1 - Q[1])]])
    tr = float(np.trace(F2 @ C2 @ F2.T)) / F2.shape[0]
    return {m: (1.0 / m - 1.0 / M_REF) * tr for m in MS}


def nearest_floors():
    """The nearest-atom leak's function-level floors per seed on the same clouds: mean of ‖f̂_m − f̂_M‖² over (450, 600] in exp_ot_leak.py's runs."""
    path = os.path.join(OUT, 'ot_leak_raw.npz')
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    out = {}
    for r in z['runs']:
        if r['set'] == 'commensurate' and r['leak'] == 'nearest':
            for m in MS:
                out[(int(r['seed']), m)] = float(np.mean(r[f'h2{m}'][r['s'] > S_END - 150]))
    return out


def fmt_counts(c):
    return '(' + ', '.join(str(int(x)) for x in c) + ')'


def report(crit, crit_eps, runs, cont, margins, sens, t0=None):
    lines = ["# Proposition D‴.20: the entropic OT leak — a soft mass constraint with a critical excess (`exp_entropic_leak.py`)", "",
             f"The three-bump instance with the matching prior (masses q = 1/3) at σ₀ = 0.15σ_crit = {SIGMA0:.5f}; the entropic plan by Newton's method on the semi-dual "
             "(three columns), every weighted stationary point checked to a residual of 10⁻⁹.", "",
             "## 1. The critical excess: the separation margin of the over-weighted class of Θ*_ε(w) along the ε-branch", "",
             "w_δ = q + δ(1.5e₁ − 0.5) takes the excess of class 1 from both other atoms, w′_δ = q + δ(e₁ − e₃) from atom 3 alone. The branch Θ*_ε(w) is continued from ε = 256 "
             "down to 1 (four steps per octave); ε_c by bisection inside the bracket of the sign change. At ε_c: the conditional plan of class 1, the eigenvalues of Cov^G₁, "
             "λ_min(A) with A = σ₀I − R^G₁ (the hypothesis of (iv)), the collective block, the exact threshold 2σ₀λ_max(A^{−1/2}Cov^G₁A^{−1/2}) and the rank-one values "
             "2σ₀p(1 − p)Δ̂_jᵀA⁻¹Δ̂_j with all of the outflow p sent to atom j = 2 or 3.", "",
             "| family | δ | margin at ε = 256, 64, 16, 4, 1 | max Newton residual on the branch | sign changes | ε_c | conditional plan of class 1 at ε_c | eigenvalues of Cov^G₁ | λ_min(A) | collective min | other classes' margins | exact threshold | rank-one (atom 2 / atom 3) and ratios to ε_c |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in crit:
        i = c['info']
        if i is None:
            lines.append(f"| {c['family']} | {c['delta']:g} | " + ", ".join(f"{mg:+.5f}" for e, mg, r in c['grid']) + f" | {c['branch_res']:.1e} | {c['n_sign_changes']} | — | — | — | — | — | — | — | — |")
            continue
        r1 = i['rank1']
        lines.append(f"| {c['family']} | {c['delta']:g} | " + ", ".join(f"{mg:+.5f}" for e, mg, r in c['grid']) + f" | {c['branch_res']:.1e} | {c['n_sign_changes']} | {c['eps_c']:.3f} | "
                     + "(" + ", ".join(f"{x:.4f}" for x in i['cond']) + f") | " + ", ".join(f"{x:.3f}" for x in i['cov_ev']) + f" | {i['Amin'][0]:+.5f} | {i['coll']:+.5f} | "
                     + ", ".join(f"{x:+.5f}" for x in i['sep'][1:]) + f" | {i['exact_thr']:.3f} | {r1[1]:.2f} / {r1[2]:.2f} ({r1[1] / c['eps_c']:.2f} / {r1[2] / c['eps_c']:.2f}) |")
    lines += ["", "The critical excess at fixed ε: δ_c(ε), the first δ along each family at which the over-weighted class's margin vanishes (continuation in δ from the w = q point of the "
              "ε-branch), the conditional fraction p_c sent away there and the small-ε scale ε/(2‖Δ̂‖²) (Δ̂ to the atom of the larger share).", "",
              "| family | ε | δ_c | p_c | ε/(2‖Δ̂‖²) | p_c/(ε/(2‖Δ̂‖²)) | w₁p_c |", "|---|---|---|---|---|---|---|"]
    for c in crit_eps:
        w1 = FAMILIES[c['family']](c['delta_c'])[0] if np.isfinite(c['delta_c']) else float('nan')
        lines.append(f"| {c['family']} | {c['eps']:g} | {c['delta_c']:.4f} | {c['p_c']:.4f} | {c['scale']:.4f} | {c['p_c'] / c['scale']:.3f} | {w1 * c['p_c']:.4f} |")
    # 2. dynamics per seed
    lines += ["", f"## 2. Dynamics over {S_END:g} slow units: m = {', '.join(str(m) for m in MS)} against M = {M_REF} from the warm clouds of `exp_ot_leak.py` (dt = {DT_S:g}), per seed", "",
              "| ε | m | counts by nearest atom: start → s = 600 (per seed) | fit error at s = 600 (mean) | W²(μ_m, ν) at s = 600 (mean) | ‖f̂_m − f̂_M‖² on (450, 600] (per seed) |",
              "|---|---|---|---|---|---|"]
    for eps in EPSS:
        R = sorted([r for r in runs if r['eps'] == eps], key=lambda r: r['seed'])
        for m in MS + (M_REF,):
            cs = "; ".join(fmt_counts(np.round(r['w'][m][0] * m)) + " → " + fmt_counts(np.round(r['w'][m][-1] * m)) for r in R)
            fit = np.mean([r['fit'][m][-1] for r in R]); W2 = np.mean([r['W2nu'][m][-1] for r in R])
            h2 = "; ".join(f"{float(np.mean(r['h2'][m][r['s'] > S_END - 150])):.2e}" for r in R) if m in MS else "—"
            lines.append(f"| {eps:g} | {m}{' (reference)' if m == M_REF else ''} | {cs} | {fit:.3e} | {W2:.3e} | {h2} |")
    # 3. continuation
    if cont:
        lines += ["", f"## 3. The ε = 4, 16, 64 runs continued to {S_CONT:g} slow units (from the stored s = 600 states, dt = {DT_CONT:g})", "",
                  "Per system: the counts at s = 600 → 6000, the number of 50-unit records at which the counts changed and the last one, the largest whitened drift at the end, "
                  "and the stability of Θ*_ε(w) at the realised weights of s = 600 and of s = 6000 (the smallest of the separation and collective minima; negative: a saddle).", "",
                  "| ε | m | seed | counts s = 600 → 6000 | count changes / last at s | drift at the end | smallest margin at w(600) / w(6000) | ‖f̂_m − f̂_M‖² on (5850, 6000] |",
                  "|---|---|---|---|---|---|---|---|"]
        for eps in EPS_CONT:
            Cs = sorted([c for c in cont if c['eps'] == eps], key=lambda c: c['seed'])
            mg = margins.get(eps, {})
            for m in MS + (M_REF,):
                for c in Cs:
                    cnt = c['counts'][m]
                    ch = np.where(np.any(cnt[1:] != cnt[:-1], axis=1))[0]
                    last = f"{c['s'][ch[-1] + 1]:.0f}" if ch.size else "—"
                    m0 = mg.get(tuple(int(x) for x in cnt[0])); m1 = mg.get(tuple(int(x) for x in cnt[-1]))
                    f = lambda d: "—" if d is None or not np.isfinite(d['coll']) else f"{min(float(np.min(d['sep'])), d['coll']):+.5f}"
                    h2 = f"{float(np.mean(c['h2'][m][c['s'] > S_CONT - 150])):.2e}" if m in MS else "—"
                    lines.append(f"| {eps:g} | {m} | {c['seed']} | {fmt_counts(cnt[0])} → {fmt_counts(cnt[-1])} | {ch.size} / {last} | {c['drift'][m]:.1e} | {f(m0)} / {f(m1)} | {h2} |")
    # 4. floors
    if sens:
        nf = nearest_floors()
        lines += ["", "## 4. The floors: Proposition D‴.5's function-level law at the stationary point of each leak, and the observed floors on the same clouds", "",
                  "Predicted E‖f̂_m − f̂_M‖² = (1/m − 1/M)tr(F₂C₂F₂ᵀ)/N with F₂ the derivatives of f̂_w along e₁ − e₃, e₂ − e₃ at Θ*(q) of the leak (ε = 0: the nearest-atom leak); observed: the "
                  "per-seed means on (450, 600] (and on (5850, 6000] after the continuation), and their ratios to the nearest-atom leak's floors of `exp_ot_leak.py` on the same clouds.", "",
                  "| leak | m | predicted | observed per seed on (450, 600] | ratio to the nearest-atom floor per seed | observed per seed on (5850, 6000] |", "|---|---|---|---|---|---|"]
        for eps in (0.0, 16.0, 64.0):
            for m in MS:
                pred = sens[eps][m]
                if eps == 0:
                    obs = [nf.get((s_, m), float('nan')) for s_ in range(N_SEEDS)]
                    lines.append(f"| nearest atom | {m} | {pred:.2e} | " + "; ".join(f"{x:.2e}" for x in obs) + " | 1 | — |")
                    continue
                R = sorted([r for r in runs if r['eps'] == eps], key=lambda r: r['seed'])
                obs = [float(np.mean(r['h2'][m][r['s'] > S_END - 150])) for r in R]
                rat = [o / nf[(r['seed'], m)] if (r['seed'], m) in nf else float('nan') for o, r in zip(obs, R)]
                Cs = sorted([c for c in cont if c['eps'] == eps], key=lambda c: c['seed']) if cont else []
                late = "; ".join(f"{float(np.mean(c['h2'][m][c['s'] > S_CONT - 150])):.2e}" for c in Cs) if Cs else "—"
                lines.append(f"| ε = {eps:g} | {m} | {pred:.2e} | " + "; ".join(f"{x:.2e}" for x in obs) + " | " + "; ".join(f"{x:.2f}" for x in rat) + f" | {late} |")
        lines += ["", "Ratios of the laws, entropic/nearest: " + "; ".join(f"ε = {eps:g}: " + ", ".join(f"{sens[eps][m] / sens[0.0][m]:.2f}" for m in MS) for eps in (16.0, 64.0)) + "."]
    lines += ["", "Numbers are produced by `experiments/exp_entropic_leak.py`; nothing is edited by hand."]
    with open(os.path.join(OUT, 'entropic_leak.md'), 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if t0 is not None:
        print(f"[{time.time() - t0:.0f}s]")


def _cont_indexed(ia):
    i, a = ia
    return i, job_continue(a)


def _continue_all(args, n_workers):
    """The continuation jobs, each result cached as soon as it is finished when EL_CACHE names a directory (a restart then skips the
    finished jobs; the cache is keyed by (seed, ε) and is valid only for the same stored states and code)."""
    import multiprocessing as mp
    import pickle
    cache = os.environ.get('EL_CACHE')
    path = lambda a: os.path.join(cache, f"cont_seed{a[0]}_eps{a[1]:g}.pkl")
    out = [None] * len(args)
    if cache:
        os.makedirs(cache, exist_ok=True)
        for i, a in enumerate(args):
            if os.path.exists(path(a)):
                with open(path(a), 'rb') as f:
                    out[i] = pickle.load(f)
                print(f"[continue] seed {a[0]}, ε = {a[1]}: cached", flush=True)
    todo = [(i, a) for i, a in enumerate(args) if out[i] is None]
    if todo:
        with mp.Pool(n_workers) as pool:
            for i, res in pool.imap_unordered(_cont_indexed, todo, chunksize=1):
                out[i] = res
                if cache:
                    with open(path(args[i]), 'wb') as f:
                        pickle.dump(res, f)
    return out


def compute_parts(parts, runs, n_workers, save=None):
    """The listed parts; `save(out)` is called after each one, so that an interrupted run keeps the finished parts."""
    out = {}
    if 'crit' in parts:
        out['crit'] = [critical_row(fam, d) for fam in FAMILIES for d in DELTAS[fam]]
        out['crit_eps'] = [critical_excess_at(fam, eps) for fam in FAMILIES for eps in (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)]
        print("[critical excess done]", flush=True)
        if save:
            save(out)
    if 'cont' in parts:
        args = [(r['seed'], r['eps'], r['final']) for r in runs if r['eps'] in EPS_CONT]
        out['cont'] = _continue_all(args, n_workers)
        if save:
            save(out)
    if 'margins' in parts:
        cont = out.get('cont')
        mg = {}
        for eps in EPS_CONT:
            counts = []
            for c in (cont or []):
                if c['eps'] == eps:
                    for m in MS + (M_REF,):
                        counts += [c['counts'][m][0], c['counts'][m][-1]]
            mg[eps] = end_margins(eps, counts)
        out['margins'] = mg
        print("[margins done]", flush=True)
        if save:
            save(out)
    if 'sens' in parts:
        out['sens'] = {eps: sensitivity_law(eps) for eps in (0.0, 16.0, 64.0)}
        print("[sensitivity laws done]", flush=True)
        if save:
            save(out)
    return out


def _save_raw(raw, data):
    np.savez(raw, **{k: np.array(v, dtype=object) for k, v in data.items()})


def main():
    import multiprocessing as mp
    t0 = time.time()
    raw = os.path.join(OUT, 'entropic_leak_raw.npz')
    keys = ('crit', 'crit_eps', 'runs', 'cont', 'margins', 'sens')
    if '--from-raw' in sys.argv:
        z = np.load(raw, allow_pickle=True)
        get = lambda k, d: (z[k].item() if z[k].dtype == object and z[k].shape == () else list(z[k])) if k in z.files else d
        report(get('crit', []), get('crit_eps', []), get('runs', []), get('cont', []), get('margins', {}), get('sens', {}), t0)
        return
    if '--only' in sys.argv:
        parts = sys.argv[sys.argv.index('--only') + 1].split(',')
        z = np.load(raw, allow_pickle=True)
        data = {k: (z[k].item() if z[k].shape == () else list(z[k])) for k in z.files if k in keys}
        if 'margins' in parts and 'cont' not in parts:
            parts.append('cont')
        merge = lambda out: (data.update(out), _save_raw(raw, data))          # each finished part is merged into the raw file at once
        compute_parts(parts, data['runs'], N_WORKERS, save=merge)
    else:
        crit_parts = compute_parts(['crit'], [], N_WORKERS)
        jobs = [(seed, eps) for eps in EPSS for seed in range(N_SEEDS)]
        with mp.Pool(N_WORKERS) as pool:
            runs = pool.map(job, jobs, chunksize=1)
        data = dict(runs=runs, **crit_parts)
        _save_raw(raw, data)
        merge = lambda out: (data.update(out), _save_raw(raw, data))
        compute_parts(['cont', 'margins', 'sens'], runs, N_WORKERS, save=merge)
    report(data.get('crit', []), data.get('crit_eps', []), data['runs'], data.get('cont', []), data.get('margins', {}), data.get('sens', {}), t0)


if __name__ == '__main__':
    main()
