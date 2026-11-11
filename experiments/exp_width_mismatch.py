"""
Proposition D‴.58: the transit under width mismatch — a feature wider than the kernels is formed as a sequence of Gaussian quadratures,
the error falling by t = (σ* − σ)/(σ* + σ) per kernel
====================================================================================================================================

Kernels of width σ = 0.12 (the program's), a feature f* = Aφ_{σ*}(· − c*) with σ* = σ(1 + (τ/σ)²)^{1/2} (`src/core/wide_feature.py`):

A. The rules on ℝ: for τ/σ in TAUS_A and K ≤ K_A, the half-squeeze rule's L² error (trapezoid on a fine grid) against the exact Mehler
   series and its first term; the Gauss–Hermite rule of the mixing law (moment matching); the half-squeeze nodes with least-squares
   amplitudes; the optimal K-kernel fit (variable projection from the half-squeeze rule) — the ratios per kernel against 1/t and 1/λ;
   the optimum's nodes against the half-squeeze nodes. A2: the series' tail — the ratios D_{2K+2j,K}/D_{2K,K} for K ≤ K_A2 (integers;
   their j-th differences in K against j!·binom(1/2, j)(−2)^j) and the tail's share against ₂F₁(−½, −½; 1; t²) − 1. A3: the scaled
   Gauss–Hermite rules (`wide_feature.scaled_rule`) at the natural Fock and the Mercer-zero scalings against the half-squeeze one (their
   ratios per kernel against (2 − λ)/λ and 1/q), and the variance of least error (`best_scale`): one parameter's share of the optimum's
   gain over the rule (referee round 48). A4: (v), the frame continued to features narrower than the kernels (σ*/σ in NARROW_A): the
   Gauss–Hermite cosine sum's error (`narrow_rule`) against the series at |t|. A5: the Mehler frames and their rules' exact series (the
   Mercer frame's rule further out, at τ/σ = 1/4 too). A6: the scaled rules' exact L² errors against one conjectured rate; A6E: the
   defect generating function's asymptotics that rate rests on (referee round 84); A6X: A6 beyond K = 160 at seven fast-branch rules
   (referee rounds 101, 123); A6C: five of them at every K = 159 … 400, the exact errors in Arb (referee round 129).
B. The transit: the sliding path (polished; the KKT exit rule) on f* (A_B, c* = C_B off the kernels' grid) sampled on X_B, in six
   numerical variants (VARIANTS_B: the polish step, the starting grid's offset, the polish's budget of evaluations, its repetition
   after the stationarity test and its gradient test; referee rounds 48, 54 and 60): after each entry the relaxed block's relative error against the optimal fit
   with as many kernels on the same samples (variable projection from the block and from the half-squeeze rule, the better), whether
   the relaxation passed the stationarity test, the first non-optimal block's polishes (status, evaluations, accepted), the block's
   nodes against the optimum's, the continuous half-squeeze error; the entries (the entrant's rate against the largest, its dual), the
   chunks redone after a stall or a failed solve, the polish's rejections. B2: variable projection from every non-optimal block. B3: each
   path's first non-optimal block, its first polish restarted from its recorded start one setting at a time (the budget, the gradient
   test). B4: the floor — the same Levenberg–Marquardt solve with the residual in double and in long double (referee round 70), about an
   anchor, every end evaluated exactly. B5: the floors' order — the double's rounding scaled by κ (referee round 84).
C. Under noise: y = f* + ξ (i.i.d. N(0, s²) on X_B), the path to K_C + 1 entries (its centres in the window ± BOX_C·σ, the thresholds
   on max|y|, the KKT exit rule): the discrepancy rule's stop (the first block with rms ≤ (1 + DELTA)s) against its prediction from (i)
   — the first K with e_K²‖f*‖²_rms ≤ s²((1 + δ)² − 1 + 2K/N), e_K the half-squeeze rule's L² error — and the probability of each stop
   from the rms's expectation and scatter (`stop_law`), the stopped fit's error against (e_K²‖f*‖²_rms + 2Ks²/N)^{1/2} (e_K B's optima
   on the samples); the kernels at the feature and elsewhere; the largest relative gap's K.

Usage: python experiments/exp_width_mismatch.py [--only A,B,C,D,B2,A5,B3,B4,B5,A6,A6E,A6X,A6C] [--variants i,j] [--rules i,j] [--from-raw]
→ results/chaos/width_mismatch.md (+ _raw.npz). --variants: B for those variants only (indices into VARIANTS_B), merged into the raw's B;
--rules: A6X for those rules only (indices into A6X_RULES), merged into the raw's A6X. WM_WORKERS (default 2). --from-raw computes
nothing (A6C, like every part, is computed only by its --only run).
"""
import os
import sys
import time
import warnings
from fractions import Fraction
from math import comb, factorial

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
from src.core import rich_transit as RT
from src.core import wide_feature as W
import exp_cohesion_leak as E
import exp_transit_mass as TM

OUT = os.environ.get('WM_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('WM_WORKERS', '2'))
SIG = 0.12
TAUS_A, K_A, K_A2 = (0.25, 0.5, 1.0, 2.0), 10, 30
NARROW_A = (0.5, 0.8)                                                 # A4: σ*/σ < 1
E_FLOOR, G_FLOOR = 1e-10, 1e-13                                      # A: the optimum not computed below (its solve stops short), the grid's floor
TAUS_B, A_B, C_B, K_B = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0), 0.8, 0.03, 12
X_B = np.linspace(-1.6, 1.6, 321)
LEVELS_C, N_DRAW_C, TAU_C, K_C, DELTA = (1e-4, 1e-3, 1e-2), 10, 1.0, 10, 0.1
KS_D, J_MAX_D = (10, 20, 40, 80, 160, 320), 7                        # D: the tail ratio's K and the ratios' degree checked
MAX_DOUBLINGS = 24                                                    # chunks to 10⁵·2²⁴ ≈ 1.7·10¹² (beyond ~10¹¹ LSODA fails from relaxed
                                                                      # states and `relax_chunk` redoes the chunk by Radau: round 48)
VARIANTS_B = (('base', 1.0, 0.0, 5000, 0, 1e-15), ('polish 2σ', 2.0, 0.0, 5000, 0, 1e-15), ('grid +0.05', 1.0, 0.05, 5000, 0, 1e-15),
              ('polish 2σ, 5·10⁴ evaluations', 2.0, 0.0, 50_000, 0, 1e-15), ('polish 2σ, 5·10⁴, repeated', 2.0, 0.0, 50_000, 3, 1e-15),
              ('polish 2σ, 5·10⁴, repeated, no gradient test', 2.0, 0.0, 50_000, 3, None))
# B: (name, the polish step, the starting grid's shift, the polish's budget of evaluations, its repetitions after the stationarity test,
# its absolute gradient test — None: none; the fourth and fifth variants: referee round 54, the sixth: round 60)
OPT_TOL = 1e-3                                                        # B: a block is the optimal fit if its error is within 0.1% of it
UNCONV = 1e-9                                                         # A3: an optimum the one-parameter rule beats by more is unconverged
BOX_C = 2.0                                                           # C: the noisy path's centres in the samples' window ± BOX_C·σ


def sigma_star(r_):
    return SIG * np.sqrt(1 + r_ ** 2)


def part_A():
    rows = []
    for r_ in TAUS_A:
        ss = sigma_star(r_)
        P = W.params(SIG, ss)
        for K in range(1, K_A + 1):
            y, a = W.half_squeeze(K, SIG, ss)
            e = W.l2_error(y, a, SIG, ss)
            s2, f2 = W.mehler_series(K, P['t'])
            if np.sqrt(s2) < G_FLOOR:                                    # the trapezoid rule's double-precision floor
                break
            ym, am = W.moment_matching(K, SIG, ss)
            V = W.scalings(SIG, ss)
            vb, eb = W.best_scale(K, SIG, ss)
            row = dict(r=r_, K=K, t=P['t'], lam=P['lam'], e_hs=e, e_series=float(np.sqrt(s2)), e_first=float(np.sqrt(f2)),
                       e_mm=W.l2_error(ym, am, SIG, ss), e_ls=W.l2_error(y, W.ls_amplitudes(y, SIG, ss), SIG, ss), y_hs=y,
                       e_fock=W.l2_error(*W.scaled_rule(K, SIG, ss, V['fock']), SIG, ss),
                       e_mercer=W.l2_error(*W.scaled_rule(K, SIG, ss, V['mercer']), SIG, ss), e_best=eb, v_best=vb / V['hs'])
            if e > E_FLOOR:
                yo, ao, eo = W.optimal_fit(y, SIG, ss)
                row.update(e_opt=eo, y_opt=yo, a_opt=ao)
            rows.append(row)
    # A2: the tail's integer ratios and the tail's share
    ratios = {K: [W.gh_defect(2 * K + 2 * j, K) / W.gh_defect(2 * K, K) for j in range(1, 5)] for K in range(1, K_A2 + 1)}
    share = []
    for r_ in TAUS_A:
        t = W.params(SIG, sigma_star(r_))['t']
        for K in (5, 10, 20, 30):
            s2, f2 = W.mehler_series(K, t)
            share.append(dict(r=r_, K=K, share=s2 / f2 - 1))
    narrow = []
    for q in NARROW_A:
        ss, t = q * SIG, (1 - q) / (1 + q)
        for K in range(1, K_A + 1):
            s2, f2 = W.mehler_series(K, t)
            if np.sqrt(s2) < G_FLOOR:
                break
            narrow.append(dict(q=q, K=K, t=t, e=W.narrow_error(*W.narrow_rule(K, SIG, ss), SIG, ss), e_series=float(np.sqrt(s2))))
    return dict(rows=rows, ratios=ratios, share=share, narrow=narrow)


def blocks(R):
    """The relaxed blocks along a path: (k kernels, τ, centres, coefficients, rms, whether the relaxation became stationary) after each
    event (hist[i + 1] is iteration i's block, relaxed[i] its flag)."""
    out = []
    for i, h in enumerate(R['hist'][1:]):
        act = np.flatnonzero(h[2] != 0)
        o = np.argsort(h[1][act])
        out.append(dict(k=int(act.size), tau=float(h[0]), c=h[1][act][o], w=h[2][act][o], rms=float(h[3]),
                        relaxed=bool(R['relaxed'][i]) if i < len(R['relaxed']) else False))
    return out


def sample_fit(Phi_of, f, y0):
    """The least-squares fit with free centres on the samples by variable projection from y0: (centres, coefficients, rms)."""
    from scipy.optimize import least_squares

    def res(c):
        P = Phi_of(c)
        return f - P @ np.linalg.lstsq(P, f, rcond=None)[0]
    sol = least_squares(res, np.asarray(y0, float), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    c = np.sort(sol.x)
    P = Phi_of(c)
    a = np.linalg.lstsq(P, f, rcond=None)[0]
    return c, a, float(np.sqrt(np.mean((f - P @ a) ** 2)))


def vp_fit(Phi_of, f, y0, max_nfev=20000):
    """`sample_fit` with its solve's stop (referee round 65): (centres, coefficients, rms, status, evaluations) — MINPACK's status 1 the
    gradient test, 2 the decrease test, 3 the step test, 4 both, 0 the budget."""
    from scipy.optimize import least_squares

    def res(c):
        P = Phi_of(c)
        return f - P @ np.linalg.lstsq(P, f, rcond=None)[0]
    sol = least_squares(res, np.asarray(y0, float), method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=max_nfev)
    c = np.sort(sol.x)
    P = Phi_of(c)
    a = np.linalg.lstsq(P, f, rcond=None)[0]
    return c, a, float(np.sqrt(np.mean((f - P @ a) ** 2))), int(sol.status), int(sol.nfev)


def vp_jacobian(K, f, c):
    """Variable projection's Jacobian in the centres, analytic (Golub–Pereyra): with P = Φ(c), the least-squares amplitudes a = P⁺f, the
    residual r = f − Pa and R = I − PP⁺, ∂r/∂c_i = −a_iR∂_cφ_i − (∂_cφ_i·r)(P⁺)ᵀe_i (referee round 65: central differences lose the
    small singular values)."""
    c = np.asarray(c, float)
    P, D = K.phi(c), K.dphi(c)
    Pp = np.linalg.pinv(P)
    a = Pp @ f
    r = f - P @ a
    return -(D - P @ (Pp @ D)) * a[None, :] - Pp.T * (D.T @ r)[None, :]


def b_job(args):
    r_, var = args
    name, step, shift, nfev, rep_, gt = next(v for v in VARIANTS_B if v[0] == var)
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L, _, c0 = TM.setup()
    c0 = c0 + shift
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(r_)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    R = RT.sliding_iss(L, K, f, c0, np.ones(c0.size), polish=True, max_events=3 * K_B, tol=0.0, max_entries=K_B + 1,
                       max_doublings=MAX_DOUBLINGS, polish_step=step, exit_kkt=True, polish_nfev=nfev, repolish=rep_,
                       polish_gtol=gt)
    rf = float(np.sqrt(np.mean(f ** 2)))
    out = []
    for i, b in enumerate(blocks(R)):
        if b['k'] == 0:
            continue
        yh, _ = W.half_squeeze(b['k'], SIG, ss)
        cands = [sample_fit(K.phi, f, b['c']), sample_fit(K.phi, f, C_B + yh)]
        c_o, a_o, rms_o = min(cands, key=lambda z: z[2])
        out.append(dict(k=b['k'], tau=b['tau'], c=b['c'], w=b['w'], e=b['rms'] / rf, e_opt=rms_o / rf, c_opt=c_o, relaxed=b['relaxed'],
                        dc=float(np.abs(b['c'] - c_o).max()) if c_o.size == b['c'].size else np.nan,
                        dc_hs=float(np.abs(b['c'] - (C_B + np.sort(yh))).max()),
                        e_hs=W.l2_error(*W.half_squeeze(b['k'], SIG, ss), SIG, ss),
                        polish=R['polish_log'][i] if i < len(R['polish_log']) else []))
    print(f"  B tau/sigma {r_} ({var}): {len(out)} blocks, errors " + ', '.join(f"{o['e']:.2e}" for o in out) +
          f" [{time.time() - t0:.0f}s]", flush=True)
    return dict(r=r_, var=var, blocks=out, events=[(float(e[0]), int(e[1]), e[2], float(e[3])) for e in R['events']],
                relaxed=R['relaxed'], stiff=R['stiff_redone'], redone=R['solver_redone'], failed=R['solver_failed'],
                rejected=R['polish_rejected'], entries=R['entry_info'], secs=time.time() - t0)


def part_B(idx=None):
    """B's paths for the variants of VARIANTS_B with indices idx (all if None), in the order of VARIANTS_B, then TAUS_B."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b_job, [(r_, v[0]) for i, v in enumerate(VARIANTS_B) if idx is None or i in idx for r_ in TAUS_B], chunksize=1)


def run_B(old, idx):
    """B's part: all paths (idx None), or the variants with indices idx rerun and merged into the stored ones (`merge_B`; referee round
    70: a partial rerun must not drop the others)."""
    return part_B() if idx is None else merge_B(old, part_B(idx))


def merge_B(old, new):
    """B's paths: the new ones replace the old of the same width and variant; in the order of VARIANTS_B, then TAUS_B."""
    names = [v[0] for v in VARIANTS_B]
    key = lambda b: (b['r'], b['var'])
    got = {key(b): b for b in (old if old is not None else [])}
    got.update({key(b): b for b in new})
    return sorted(got.values(), key=lambda b: (names.index(b['var']), TAUS_B.index(b['r'])))


def b2_job(args):
    """B2 (referee round 54): from one of B's non-optimal blocks, variable projection (Levenberg–Marquardt on the centres) started at the
    block's own centres and at the half-squeeze nodes — whether the block lies in the basin of the best optimum found (its relaxation
    unfinished) or at another stationary point; dc, the block's largest distance from the better optimum's nodes; and at the block's
    centres (least-squares amplitudes) the condition numbers of the two problems' Jacobians — the joint fit in (w, c), as the polish
    solves it, and variable projection's in c alone (`vp_jacobian`, analytic), both with the centres in units of σ, and the ratio of their
    smallest singular values (referee rounds 60, 65) — and the same in the polish's own units (the centres in x) and with every column
    scaled to unit norm (referee round 70: the comparison depends on the units); where variable projection from the block stops: its
    status, its closest nodes (σ) and the largest cosine between its residual and a column of its Jacobian (zero at a stationary point)."""
    r_, var, k, c = args
    warnings.filterwarnings('ignore')
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(r_)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    rf = float(np.sqrt(np.mean(f ** 2)))
    yh, _ = W.half_squeeze(k, SIG, ss)
    cb, ab, eb, st, nf = vp_fit(K.phi, f, c)
    cr, _, er, _, _ = vp_fit(K.phi, f, C_B + yh)
    c = np.sort(np.asarray(c, float))
    P = K.phi(c)
    wl = np.linalg.lstsq(P, f, rcond=None)[0]
    sv = lambda J: np.linalg.svd(J, compute_uv=False)
    eq = lambda J: J / np.where(np.linalg.norm(J, axis=0) > 0, np.linalg.norm(J, axis=0), 1.0)[None, :]
    Jj, Jv = np.c_[P, K.dphi(c) * wl[None, :]], vp_jacobian(K, f, c)  # the polish's own units (the centres in x)
    sj, sp = sv(np.c_[P, K.dphi(c) * wl[None, :] * SIG]), sv(Jv * SIG)
    sjx, spx, sje, spe = sv(Jj), sv(Jv), sv(eq(Jj)), sv(eq(Jv))
    Je, re = vp_jacobian(K, f, cb), f - K.phi(cb) @ ab
    cos = float(np.max(np.abs(Je.T @ re) / (np.linalg.norm(Je, axis=0) * np.linalg.norm(re)))) if np.linalg.norm(re) > 0 else 0.0
    return dict(r=r_, var=var, k=k, e_block=eb / rf, e_rule=er / rf, dc=float(np.abs(c - (cb if eb <= er else cr)).max()),
                cond_joint=float(sj[0] / sj[-1]), cond_vp=float(sp[0] / sp[-1]), smin_ratio=float(sp[-1] / sj[-1]),
                cond_joint_x=float(sjx[0] / sjx[-1]), cond_vp_x=float(spx[0] / spx[-1]), smin_ratio_x=float(spx[-1] / sjx[-1]),
                smin_ratio_eq=float(spe[-1] / sje[-1]), vp_status=st, vp_nfev=nf,
                vp_gap=float(np.diff(cb).min() / SIG) if cb.size > 1 else np.nan, vp_cos=cos)


def seed_rule(K, f, w, c, s, cp, step, wfrac=1e-2):
    """`sliding_iss`'s acceptance of a polish from (w, c) to the centres cp: every kernel whose |w| is at most wfrac·max|w| (near its seed)
    lands within step kernel widths of its start or of its correlation's ascent target (`RT.ascent_target`). Returns whether it is
    accepted and the largest such kernel's move (σ; the nearer of the two)."""
    seedy = np.flatnonzero(np.abs(w) <= wfrac * np.abs(w).max())
    r0 = f - K.phi(c) @ w
    mv = [min(abs(cp[i] - c[i]), abs(cp[i] - RT.ascent_target(K, r0, c[i], s[i])[0])) / K.sig for i in seedy]
    return dict(accepted=bool(all(m <= step for m in mv)), move=float(max(mv)) if mv else 0.0)


B3_BIG = 500_000                                                      # B3 (referee round 70): the larger budget
B3_KEYS = ('path', 'budget', 'no test', 'budget, no test')            # B3's restarts: the path's settings, then one or both changed


class _StopPath(Exception):
    """`path_polishes`' recorder: the call asked for is recorded, the rest of the path is not needed."""


def path_polishes(r_, var, stop_after=None):
    """B's path (r_, var) rerun as `b_job`'s with every `polish_block` call recorded (referee round 70): its start (w, c, s) and settings,
    in the path's order (R['polish_log'] flattened: each call appends one record to its block's log). With stop_after = n the path stops
    once call n (from 0) is recorded — before that call's solve — and R is None. Returns (R, K, f, calls)."""
    name, step, shift, nfev, rep_, gt = next(v for v in VARIANTS_B if v[0] == var)
    L, _, c0 = TM.setup()
    K = RT.Bumps(X_B, SIG)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))
    calls, real = [], RT.polish_block

    def rec(K_, fstar, w, c, s, cbox=None, max_nfev=5000, info=False, gtol=1e-15):
        calls.append(dict(w=np.array(w, float), c=np.array(c, float), s=np.array(s, float), max_nfev=max_nfev, gtol=gtol))
        if stop_after is not None and len(calls) > stop_after:
            raise _StopPath
        return real(K_, fstar, w, c, s, cbox=cbox, max_nfev=max_nfev, info=info, gtol=gtol)
    RT.polish_block = rec
    try:
        R = RT.sliding_iss(L, K, f, c0 + shift, np.ones(c0.size), polish=True, max_events=3 * K_B, tol=0.0, max_entries=K_B + 1,
                           max_doublings=MAX_DOUBLINGS, polish_step=step, exit_kkt=True, polish_nfev=nfev, repolish=rep_,
                           polish_gtol=gt)
    except _StopPath:
        R = None
    finally:
        RT.polish_block = real
    return R, K, f, calls


def first_nonoptimal(b):
    """B's path b (its stored blocks): the first block with polishes whose error is at least 1 + OPT_TOL times the optimum's, as
    dict(n0 — the number of polish calls before its first —, k, e, e_opt, n_polish, first — its first polish's stored stop —), or None
    (the stored blocks are those with kernels: a block without kernels has no polish, so n0 counts every call)."""
    n0 = 0
    for o in b['blocks']:
        n_i = len(o['polish']) if 'polish' in o else 0
        if o['k'] > 0 and n_i > 0 and abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL:
            return dict(n0=n0, k=int(o['k']), e=float(o['e']), e_opt=float(o['e_opt']), n_polish=n_i, first=dict(o['polish'][0]))
        n0 += n_i
    return None


def b3_job(args):
    """B3 (referee rounds 65, 70): the path's first non-optimal block (its error over the optimum's at least 1 + OPT_TOL), its first
    polish restarted from that polish's recorded start (`path_polishes`) by `polish_block`, one setting at a time — the path's budget and
    gradient test; the budget B3_BIG with the test; the path's budget without it; both changed —: each the solve's stop (TRF's status: 0
    the budget, 1 the gradient test, 2 the decrease test, 3 the step test), its evaluations, its error over the optimum's and the path's
    acceptance of its end (`seed_rule` with the path's polish step: the first polish has no improvement test)."""
    r_, var = args[:2]
    _, step, _, nfev, _, gt = next(v for v in VARIANTS_B if v[0] == var)
    warnings.filterwarnings('ignore')
    t0 = time.time()
    if len(args) > 2:                                                    # the block known from B's raw (`first_nonoptimal`): the path
        blk = args[2]                                                    # rerun only to its first polish
        _, K, f, calls = path_polishes(r_, var, stop_after=blk['n0'])
        if len(calls) <= blk['n0'] or calls[blk['n0']]['c'].size != blk['k']:
            raise RuntimeError(f"B3 {r_} {var}: the rerun's polish {blk['n0']} is not the stored block's ({len(calls)} calls)")
        return b3_restart(r_, var, K, f, calls[blk['n0']], blk['k'], blk['e_opt'], blk['e'] / blk['e_opt'], blk['n_polish'], step,
                          nfev, gt, t0, first=blk.get('first'))
    R, K, f, calls = path_polishes(r_, var)
    rf, ss, n0 = float(np.sqrt(np.mean(f ** 2))), sigma_star(r_), 0
    for i, b in enumerate(blocks(R)):
        n_i = len(R['polish_log'][i]) if i < len(R['polish_log']) else 0
        if b['k'] > 0 and n_i > 0:
            yh, _ = W.half_squeeze(b['k'], SIG, ss)
            e_opt = min(sample_fit(K.phi, f, b['c'])[2], sample_fit(K.phi, f, C_B + yh)[2]) / rf
            if abs(b['rms'] / rf / e_opt - 1) >= OPT_TOL:
                return b3_restart(r_, var, K, f, calls[n0], b['k'], e_opt, b['rms'] / rf / e_opt, n_i, step, nfev, gt, t0)
        n0 += n_i
    return dict(r=r_, var=var, k=None, secs=time.time() - t0)


def b3_restart(r_, var, K, f, st, k, e_opt, e_block, n_polish, step, nfev, gt, t0, first=None):
    """B3's restarts of one recorded polish start st, one setting at a time (B3_KEYS), each judged by `seed_rule` with the path's
    polish step. first: B's stored first polish of the block — the 'path' restart must reproduce its stop, evaluations and acceptance
    (referee rounds 77, 81: checked before the larger restarts run; a RuntimeError otherwise)."""
    rf = float(np.sqrt(np.mean(f ** 2)))
    out = dict(r=r_, var=var, k=k, e_opt=e_opt, e_block=e_block, n_polish=n_polish,
               e_start=float(np.sqrt(np.mean((f - K.phi(st['c']) @ st['w']) ** 2))) / rf / e_opt,
               settings=dict(nfev=st['max_nfev'], gtol=st['gtol']))
    for key, nf, g in zip(B3_KEYS, (nfev, B3_BIG, nfev, B3_BIG), (gt, gt, None, None)):
        wp, cp, info = RT.polish_block(K, f, st['w'], st['c'], st['s'], max_nfev=nf, info=True, gtol=g)
        out[key] = dict(info, ratio=float(np.sqrt(np.mean((f - K.phi(cp) @ wp) ** 2))) / rf / e_opt,
                        **seed_rule(K, f, st['w'], st['c'], st['s'], cp, step))
        fs = first or {}
        if key == 'path' and any(fs.get(x) is not None and out['path'][x] != fs[x] for x in ('status', 'nfev', 'accepted')):
            raise RuntimeError(f"B3 {r_} {var}: the 'path' restart {out['path']['status']}/{out['path']['nfev']}/"
                               f"{out['path']['accepted']} is not the stored first polish {fs.get('status')}/{fs.get('nfev')}/"
                               f"{fs.get('accepted')}")
    out['secs'] = time.time() - t0
    print(f"  B3 tau/sigma {r_} ({var}) K {k}: start {out['e_start']:.4g}, " + ', '.join(
        f"{key} {out[key]['status']}/{out[key]['nfev']} {out[key]['ratio']:.5f}" for key in B3_KEYS)
          + f" [{out['secs']:.0f}s]", flush=True)
    return out


def b3_jobs(B):
    """B3's jobs: each of B's paths with a non-optimal block (its error at least 1 + OPT_TOL times the optimum's) — with the first such
    block that was polished (`first_nonoptimal`) where the stored blocks record their polishes, so that the rerun stops at its first
    polish."""
    out = []
    for b in B:
        if any(abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL for o in b['blocks']):
            blk = first_nonoptimal(b) if all('polish' in o and 'k' in o for o in b['blocks']) else None
            out.append((b['r'], b['var']) if blk is None else (b['r'], b['var'], blk))
    return out


def part_B3(B):
    """B3: `b3_job` for `b3_jobs(B)`."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b3_job, b3_jobs(B), chunksize=1)


def b3_binds(q, ref=1.0):
    """What bounds a B3 block's range end (referee rounds 70, 77), from the restarts' ends that reach the optimum — judged against ref, the
    optimum's error over the stored optimum's (below B4_MAX B4's long-double end) — and that the path's seed rule accepts: the fewest
    settings whose change gives such an end (the budget, the gradient test, either alone, or both); where every end that reaches it is
    rejected, the fewest settings that reach it "(and the radius)", and the radius alone where the path's own settings reach it; "neither"
    with the both-changed restart's stop (its budget spent, or the solve's own stop); an inconsistency where the path's own polish reaches
    it, accepted."""
    reach = lambda k: q[k]['ratio'] / ref < 1 + OPT_TOL
    ok = lambda k: reach(k) and q[k]['accepted']
    name = lambda ks: {('budget',): 'the budget', ('no test',): 'the gradient test'}.get(tuple(ks), 'the budget or the gradient test')
    if reach('path'):
        return 'the radius' if not q['path']['accepted'] else "nothing (inconsistent: the path's own polish reaches it, accepted)"
    for test in (ok, reach):
        one = [k for k in ('budget', 'no test') if test(k)]
        tail = '' if test is ok else ' (and the radius)'
        if one:
            return name(one) + tail
        if test('budget, no test'):
            return 'the budget and the gradient test' + tail
    st = q['budget, no test']['status']
    return 'neither (' + ('the budget spent' if st == 0 else f"the solve's stop, status {st}") + ')'


def b3_lines(B3, B4=None, B=None):
    """B3's table and summary (referee rounds 70, 77, 81): per block the restarts' stops, errors and the seed rule, and what bounds its
    end — judged, below B4_MAX, against B4's best end (`b4_best`: the lowest exact error any B4 solve reached there, over the row's own
    optimum), the better optimum there; with B, the 'path' restarts against the first polishes B stored (their stops and evaluations)."""
    if B3 is None or len(B3) == 0:                                    # raw-loaded: an object array
        return []
    rows = [q for q in B3 if q.get('k') is not None]
    refs = b4_best(B4)                                                 # referee round 81: B4's best end, absolute, over each row's
    ref = lambda q: refs[(q['r'], q['k'])]['err'] / q['e_opt'] if (q['r'], q['k']) in refs else 1.0   # own optimum
    cell = lambda x: f"{x['status']}/{x['nfev']}, {x['ratio']:.5f}, {x['move']:.2f}σ {'a' if x['accepted'] else 'r'}"
    out = ['', "### B3. The range ends: the polish restarted from its recorded start, one setting at a time (referee rounds 65, 70, 77)", '',
           f"Each path's first non-optimal block (its error at least {1 + OPT_TOL:g} times the optimum's): the path rerun with every polish's "
           "start recorded, and from the block's first polish's start — its error over the optimum's — the polish's own solve "
           f"(`polish_block`, TRF) with the path's budget and gradient test, with {B3_BIG:g} evaluations and the test, with the path's budget "
           "without the test, and with both changed: each its stop (status 0: the budget; 1: the gradient test; 2: the decrease test; 3: the "
           "step test)/evaluations, its error over the optimum's and the path's seed rule on its end (the largest move of a kernel near its "
           f"seed, σ; a: accepted within the path's polish step, r: rejected); the block the path ended with; † the optimum below {E_FLOOR:g}:",
           '', "| τ/σ | variant | K | optimum | start | the path's block | " + ' | '.join(B3_KEYS) + ' | bounded by |',
           '|---|---|---|---|---|---|' + '---|' * len(B3_KEYS) + '---|']
    for q in rows:
        out.append(f"| {q['r']:g} | {q['var']} | {q['k']} | {q['e_opt']:.2e}{'†' if q['e_opt'] < E_FLOOR else ''} | {q['e_start']:.4g} | "
                   f"{q['e_block']:.4g} | " + ' | '.join(cell(q[k]) for k in B3_KEYS) + f" | {b3_binds(q, ref(q))} |")
    up = [q for q in rows if q['e_opt'] >= E_FLOOR]
    cnt = {}
    for q in up:
        cnt.setdefault(b3_binds(q, ref(q)), []).append(f"{q['r']:g} ({q['var']})")
    lo = [q for q in rows if q['e_opt'] < E_FLOOR]
    both = lambda q: q['budget, no test']
    if B is not None and len(B):
        first = {(b['r'], b['var']): (first_nonoptimal(b) or {}).get('first') for b in B}
        cmp = [(q, first.get((q['r'], q['var']))) for q in rows]
        cmp = [(q, fs) for q, fs in cmp if fs is not None]
        same = [q for q, fs in cmp if (q['path']['status'], q['path']['nfev']) == (fs.get('status'), fs.get('nfev'))]
        out += ['', f"The 'path' restart reproduces the first polish B stored — its stop and evaluations — at {len(same)} of {len(cmp)} blocks"
                + (" (not at " + '; '.join(f"{q['r']:g} ({q['var']})" for q, _ in cmp if q not in same) + ")" if len(same) < len(cmp) else '')
                + "."]
    out += ['', f"Above {E_FLOOR:g} ({len(up)} blocks: τ/σ (variant)) the end is bounded by", '']
    out += [f"- {k}: {len(v)} — {'; '.join(v)}" + (';' if i < len(cnt) - 1 else '.') for i, (k, v) in enumerate(sorted(cnt.items()))]
    if lo:
        out += ['', f"Below it ({len(lo)} blocks) the restart with both changed stops at " + ', '.join(
            f"{both(q)['ratio']:.4f} ({both(q)['status']}/{both(q)['nfev']}: {b3_binds(q, ref(q))})" for q in lo) + " times the stored "
            "optimum" + (", " + ', '.join(f"{both(q)['ratio'] / ref(q):.4f}" for q in lo) + " times B4's best end" if all(
                (q['r'], q['k']) in refs for q in lo) else '') + "."]
    by = {}
    for q in rows:
        by.setdefault(float(f"{q['e_opt']:.3g}"), []).append(both(q)['ratio'] / ref(q) - 1)
    out += ['', "The restart with both changed (the polish's TRF): its end over the optimum, minus 1 — below " + f"{B4_MAX:g} over B4's "
            "best end —, by the optimum (one solver from one kind of start; above B4_MAX a negative value is the end below the stored "
            "optimum, converged only to that): " + '; '.join(
                f"{e:.3g}: {min(v):+.1e}" + (f" … {max(v):+.1e}" if max(v) > min(v) else '') for e, v in sorted(by.items(), reverse=True))
            + '.', '']
    return out



B4_CAP, B4_MAX = 400_000, 2e-9                                        # B4 (referee round 70): the solve's iterations; the blocks below
B4_PREC = (('double', np.float64), ('long double', np.longdouble))    # B4: the parameters' and the residual's arithmetic
B4_ANCHOR, B4_FLAT = 2000, 1e-10                                       # B4 (referee rounds 77, 81): the anchored solve's period (accepted
                                                                       # steps, counted across re-anchors) and its stop: the exact rms over a
                                                                       # period falls by less than B4_FLAT of itself
B4_MARKS = (100_000, 200_000, 300_000)                                 # B4 (round 77): the anchored double's end after these iterations
EXACT_PREC = 60                                                        # B4 (referee round 81): the exact residual's digits


def _dec(v):
    """A double, or a long double as its double-double hi + lo (the 64-bit significand splits exactly into two doubles), as an exact
    Decimal (the sum formed at 200 digits, whatever the caller's context)."""
    from decimal import Decimal, localcontext
    if isinstance(v, np.longdouble):
        hi = float(v)
        with localcontext() as ctx:
            ctx.prec = 200
            return +(Decimal(hi) + Decimal(float(v - np.longdouble(hi))))
    return Decimal(float(v))


def exact_resid(f, p, n, prec=EXACT_PREC, dec=False):
    """The residual f − Φ(c)w on X_B, evaluated exactly to its double rounding (referee round 81): the samples, f, σ and p = (w, c) (n
    kernels, double or long double) taken exactly, the arithmetic in prec-digit decimal (exp correctly rounded). Returns (the residual,
    each sample's rounded to double; its rms, rounded to double) and, if dec, the rms itself, a Decimal (referee round 84: the ends'
    differences formed exactly)."""
    from decimal import Decimal, localcontext
    p = np.asarray(p)
    with localcontext() as ctx:
        ctx.prec = prec
        ws, cs = [_dec(v) for v in p[:n]], [_dec(v) for v in p[n:]]
        s2 = 2 * Decimal(SIG) ** 2
        r, tot = [], Decimal(0)
        for xi, fi in zip(X_B, np.asarray(f, float)):
            x_ = Decimal(float(xi))
            ri = Decimal(float(fi)) - sum((w_ * (-(x_ - c_) ** 2 / s2).exp() for w_, c_ in zip(ws, cs)), Decimal(0))
            r.append(float(ri))
            tot += ri * ri
        rms = (tot / len(X_B)).sqrt()
        return (np.array(r), float(rms)) + ((rms,) if dec else ())


def exact_rms(f, p, n):
    """`exact_resid`'s rms (referee round 81: B4's instrument, finer than the floors it measures)."""
    return exact_resid(f, p, n)[1]


GN_PREC = 80                                                           # B4 (referee round 84): the Gauss–Newton projections' digits


def exact_parts(f, p, n, prec=GN_PREC):
    """The residual f − Φ(c)w on X_B and its Jacobian in (w, c) — the rows (−φ_k(x), −w_kφ_k(x)(x − c_k)/σ²) —, exactly to prec digits
    (referee round 84): the samples, f, σ and p (n kernels, double or long double) taken exactly (`_dec`). Returns (r, J), lists of
    Decimals (J by samples)."""
    from decimal import Decimal, localcontext
    p = np.asarray(p)
    with localcontext() as ctx:
        ctx.prec = prec
        ws, cs = [_dec(v) for v in p[:n]], [_dec(v) for v in p[n:]]
        sg2 = Decimal(SIG) ** 2
        s2 = 2 * sg2
        r, J = [], []
        for xi, fi in zip(X_B, np.asarray(f, float)):
            x_ = Decimal(float(xi))
            ph = [(-(x_ - c_) ** 2 / s2).exp() for c_ in cs]
            r.append(Decimal(float(fi)) - sum((w_ * e_ for w_, e_ in zip(ws, ph)), Decimal(0)))
            J.append([-e_ for e_ in ph] + [-w_ * e_ * (x_ - c_) / sg2 for w_, e_, c_ in zip(ws, ph, cs)])
        return r, J


def gn_share(J, r, vs, prec=GN_PREC):
    """vᵀP_Jv/(2‖r‖²) for each v in vs, P_J the orthogonal projection on J's range — at v = r the exact rms's relative decrease that a
    Gauss–Newton step predicts —, from the normal equations' Cholesky factor in prec-digit decimal (referee round 84: J's condition
    number, ≲ 10⁹, costs its square ≲ 18 of the digits): vᵀP_Jv = ‖L⁻¹Jᵀv‖², LLᵀ = JᵀJ."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        cols = list(zip(*J))
        m = len(cols)
        L = []
        for i in range(m):
            row = []
            for j in range(i + 1):
                Lj = row if j == i else L[j]
                s = sum((a * b for a, b in zip(cols[i], cols[j])), Decimal(0)) - sum((row[k] * Lj[k] for k in range(j)), Decimal(0))
                row.append(s.sqrt() if j == i else s / L[j][j])
            L.append(row)
        rr = sum((x * x for x in r), Decimal(0))
        out = []
        for v in vs:
            z = []
            for i in range(m):
                b = sum((a * x for a, x in zip(cols[i], v)), Decimal(0))
                z.append((b - sum((L[i][k] * z[k] for k in range(i)), Decimal(0))) / L[i][i])
            out.append(float(sum((y * y for y in z), Decimal(0)) / (2 * rr)))
        return out


def b4_gap(f, p, k):
    """The exact rms's relative decrease that a Gauss–Newton step predicts at p (referee round 84): ‖P_Jr‖²/(2‖r‖²), r and J exact
    (`exact_parts`, `gn_share`) — a first-order distance to the local minimum."""
    r, J = exact_parts(f, p, k)
    return gn_share(J, r, [r])[0]


def direct_double(f, p, k):
    """f − Φ(c)w by direct subtraction in double, as `lm_fit` evaluates it at an anchor (the same operations on the same values)."""
    x, fd, s2 = X_B.astype(np.float64), np.asarray(f, float).astype(np.float64), 2 * np.float64(SIG) ** 2
    q = np.asarray(p, float)
    w, c = q[:k].copy(), q[k:].copy()
    return fd - np.exp(-(x[:, None] - c[None, :]) ** 2 / s2) @ w


def b4_floors(f, p, k, pa, r0):
    """The anchored double's floor, predicted (referee rounds 81, 84): a rounding e₀ frozen in the residual moves the minimiser by −J⁺e₀
    and raises the exact rms by ‖P_Je₀‖²/(2‖r‖²) — P_J the projection on the Jacobian's range, r the exact residual, both at the end p,
    exactly (`exact_parts`, `gn_share`) —, e₀ the direct subtraction in double less the exact residual: at the end ('end', round 81's
    `anchor_floor`, an anchor the solve never used), at the last anchor pa, whose r₀ (as the solve froze it) the solve's floor is
    ('last'), and its part from 2σ² rounded to double alone at pa ('s2': −Σ_kw_k(e^{−u_k²/fl(2σ²)} − e^{−u_k²/(2σ²)}), u_k = x − c_k,
    the same at every anchor); and the gap ‖P_Jr‖²/(2‖r‖²) at p ('gap', `b4_gap`). Returns a dict of the four."""
    from decimal import Decimal, localcontext
    r, J = exact_parts(f, p, k)
    ra = exact_parts(f, pa, k)[0]
    pa = np.asarray(pa, float)
    with localcontext() as ctx:
        ctx.prec = GN_PREC
        e_end = [Decimal(float(a)) - b for a, b in zip(direct_double(f, p, k), r)]
        e_last = [Decimal(float(a)) - b for a, b in zip(np.asarray(r0, float), ra)]
        ws, cs = [Decimal(float(v)) for v in pa[:k]], [Decimal(float(v)) for v in pa[k:]]
        s2x, s2d = 2 * Decimal(SIG) ** 2, Decimal(float(2 * np.float64(SIG) ** 2))
        e_s2 = []
        for xi in X_B:
            x_ = Decimal(float(xi))
            e_s2.append(-sum((w_ * ((-(x_ - c_) ** 2 / s2d).exp() - (-(x_ - c_) ** 2 / s2x).exp()) for w_, c_ in zip(ws, cs)),
                             Decimal(0)))
    gap, end, last, s2 = gn_share(J, r, [r, e_end, e_last, e_s2])
    return dict(gap=gap, end=end, last=last, s2=s2)


def kappa_blend(f, p, n, rd, kappa):
    """r_κ = r_ld + κ(r_d − r_ld), rounded to double (referee round 84, D3ce-1): rd = r_d, the residual f − Φ(c)w by direct subtraction in
    double at p, r_ld the same subtraction in long double at the same parameters, κ the scale of r_d's rounding — κ = 1 gives rd exactly
    (r_d − r_ld is exact by Sterbenz's lemma, and so is the sum), κ = 0 r_ld rounded to double."""
    LD = np.longdouble
    q = np.asarray(p).astype(LD)
    rl = np.asarray(f, float).astype(LD) - np.exp(-(X_B.astype(LD)[:, None] - q[None, n:]) ** 2 / (2 * LD(SIG) ** 2)) @ q[:n]
    return (rl + LD(kappa) * (np.asarray(rd).astype(LD) - rl)).astype(np.float64)


def lm_fit(f, w, c, dtype, max_iter=B4_CAP, mu0=1e-20, anchored=False, anchor_every=B4_ANCHOR, marks=None, stall_reanchor=True,
           kappa=None, anchor_log=None):
    """B4 (referee round 70): Levenberg–Marquardt on the joint fit in (w, c) on X_B — the parameters and the residual f − Φ(c)w in dtype
    (np.float64 or np.longdouble), the Jacobian and the damped step in double: the least-squares problem [J; μ^{1/2}D]d ≈ [−r; 0], D the
    columns' largest norms so far (its condition is the Jacobian's, not its square's). A step is accepted when its predicted decrease is
    positive and the objective ½‖r‖² computed in dtype decreases; Nielsen's damping from μ = mu0, near Gauss–Newton — the slow modes' curvature is ~10⁻¹⁷ of the largest,
    and a larger damping shrinks their steps' predicted decrease below the objective's rounding. Stops after max_iter iterations ('cap')
    or once μ exceeds 10¹⁰⁰ without a decrease ('stall': no step the solve computes decreases the objective as its residual evaluates
    it). Referee round 77, anchored: the residual about an anchor (w₀, c₀), r = r₀ − Σ_k φ₀_k[(w_k − w₀_k) + w_k·expm1(Δ_k)], Δ_k =
    (c_k − c₀_k)(2x − c_k − c₀_k)/(2σ²), r₀ = f − Φ(c₀)w₀ by direct subtraction — its rounding frozen until the next anchor, the change
    computed with an error proportional to the parameters' displacement from the anchor —, all in dtype (referee round 81: double or long
    double). At a stall the solve re-anchors at its point and resets μ to mu0 (stall_reanchor; without it, it stops there), unless no
    step was accepted since the last anchor; every anchor_every accepted steps, counted across those re-anchors, it re-anchors too, and
    stops ('flat', referee round 81) if the exact rms (`exact_rms`) fell by less than B4_FLAT of itself over that period. Referee round 84
    (D3ce-1, D3ce-4): kappa — the direct subtraction (the residual, or an anchor's r₀) with its rounding scaled by κ (`kappa_blend`;
    dtype np.float64 only); anchor_log — a dict given the anchors' counts by kind ('start', 'stall', 'period') and the last anchor's
    parameters ('p'), frozen r₀ ('r0'), iteration ('it') and kind ('kind'). marks: a dict whose keys are iteration counts, given the
    parameters after that many iterations. Returns (the parameters (w, c) in dtype, iterations, accepted steps, stop)."""
    if kappa is not None and dtype is not np.float64:
        raise ValueError("lm_fit: kappa scales the double's rounding (dtype np.float64)")
    mix = (lambda rd, p_: rd) if kappa is None else (lambda rd, p_: kappa_blend(f, p_, len(w), rd, kappa))
    x, fd = X_B.astype(dtype), np.asarray(f, float).astype(dtype)
    n = len(w)
    p = np.r_[np.asarray(w, float), np.asarray(c, float)].astype(dtype)
    s2 = 2 * dtype(SIG) ** 2
    anc, prog = {}, {}
    if anchored:
        def set_anchor(p_, kind, it_):
            anc['w'], anc['c'] = p_[:n].copy(), p_[n:].copy()
            anc['P'] = np.exp(-(x[:, None] - anc['c'][None, :]) ** 2 / s2)
            anc['r'] = mix(fd - anc['P'] @ anc['w'], p_)
            anc['acc'] = 0
            if anchor_log is not None:
                anchor_log[kind] = anchor_log.get(kind, 0) + 1
                anchor_log.update(p=p_.copy(), r0=anc['r'].copy(), it=it_, kind=kind)

        def resid(p_):
            dl = (p_[None, n:] - anc['c'][None, :]) * (2 * x[:, None] - p_[None, n:] - anc['c'][None, :]) / s2
            return anc['r'] - (anc['P'] * ((p_[:n] - anc['w'])[None, :] + p_[None, :n] * np.expm1(dl))).sum(axis=1)
        set_anchor(p, 'start', 0)
        prog.update(acc=0, rms=exact_rms(f, p, n))
    else:
        resid = lambda p_: mix(fd - np.exp(-(x[:, None] - p_[None, n:]) ** 2 / s2) @ p_[:n], p_)

    def jac(p_):
        q = p_.astype(float)
        u = X_B[:, None] - q[None, n:]
        P = np.exp(-u ** 2 / (2 * SIG ** 2))
        return -np.c_[P, P * u / SIG ** 2 * q[None, :n]]
    r = resid(p)
    F, J = r @ r / 2, jac(p)
    D = np.linalg.norm(J, axis=0)
    mu, nu, it, acc = mu0, 2.0, 0, 0
    while it < max_iter:
        if marks is not None and it in marks:
            marks[it] = p.copy()
        it += 1
        r64 = r.astype(float)
        d = np.linalg.lstsq(np.r_[J, np.sqrt(mu) * np.diag(D)], np.r_[-r64, np.zeros(2 * n)], rcond=None)[0]
        Jd = J @ d
        pred = -(Jd @ r64) - 0.5 * (Jd @ Jd)
        pn = p + d.astype(dtype)
        rn = resid(pn)
        Fn = rn @ rn / 2
        if pred > 0 and Fn < F:
            rho = float(F - Fn) / pred
            p, r, F, acc = pn, rn, Fn, acc + 1
            mu, nu = max(mu * max(1 / 3, 1 - (2 * rho - 1) ** 3), 1e-40), 2.0
            J = jac(p)
            D = np.maximum(D, np.linalg.norm(J, axis=0))
            if anchored:
                anc['acc'] += 1
                prog['acc'] += 1
                if prog['acc'] >= anchor_every:                         # referee round 81: the period's exact decrease; re-anchor
                    e_ = exact_rms(f, p, n)
                    if prog['rms'] - e_ < B4_FLAT * e_:
                        return p, it, acc, 'flat'
                    prog.update(acc=0, rms=e_)
                    set_anchor(p, 'period', it)
                    r = resid(p)
                    F = r @ r / 2
        else:
            mu, nu = mu * nu, 2 * nu
            if mu > 1e100:
                if anchored and stall_reanchor and anc['acc'] > 0:      # a stall since the last anchor: re-anchor, μ reset
                    set_anchor(p, 'stall', it)
                    r = resid(p)
                    F, mu, nu = r @ r / 2, mu0, 2.0
                    continue
                return p, it, acc, 'stall'
    return p, it, acc, 'cap'


def lm_rms(f, p, n):
    """The rms of f − Φ(c)w on X_B in long double, p = (w, c) with n kernels."""
    x, fd, p = X_B.astype(np.longdouble), np.asarray(f, float).astype(np.longdouble), np.asarray(p).astype(np.longdouble)
    r = fd - np.exp(-(x[:, None] - p[None, n:]) ** 2 / (2 * np.longdouble(SIG) ** 2)) @ p[:n]
    return float(np.sqrt(np.mean(r * r)))


def b4_end(f, p, k, rf, e_opt, gap=None, **kw):
    """One B4 end (referee rounds 81, 84): its error over the optimum's evaluated exactly (`exact_resid`; the rms also as a decimal string,
    for exact differences) and as `lm_rms` evaluates it, the exact rms's decrease that a Gauss–Newton step predicts there (`b4_gap`,
    unless given), its parameters (a long double as the doubles hi and lo, hi + lo exact) and kw."""
    p = np.asarray(p)
    hi = p.astype(float)
    lo = (p - hi.astype(p.dtype)).astype(float) if p.dtype == np.longdouble else np.zeros_like(hi)
    _, e, d = exact_resid(f, p, k, dec=True)
    return dict(ratio=e / rf / e_opt, rms=str(d), gap=b4_gap(f, p, k) if gap is None else gap, ratio_ld=lm_rms(f, p, k) / rf / e_opt,
                p=hi.tolist(), p_lo=lo.tolist(), **kw)


def anchor_info(log, p, it):
    """An anchored solve's anchors (referee round 84, `lm_fit`'s anchor_log): its re-anchors at stalls and by period, and its last
    anchor's kind, the iterations from it to the end and its largest parameter displacement from the end p."""
    return dict(stall=log.get('stall', 0), period=log.get('period', 0), kind=log['kind'], before=int(it - log['it']),
                dist=float(np.abs(np.asarray(p, float) - np.asarray(log['p'], float)).max()))


def b4_job(args):
    """B4 (referee round 70): from one start — the path's block ('block') or the stored optimum's centres with least-squares amplitudes
    ('optimum') — `lm_fit` with the residual in double and in long double (B4_PREC): the start's error and each end's over the stored
    optimum's, its iterations, accepted steps and stop, and long double's ε (an error if not extended); referee round 77: the same solve
    in double with the residual about an anchor ('anchored'), and the polish's own TRF with B3_BIG evaluations and no gradient test
    ('trf': its status and evaluations); the anchored end's error also after B4_MARKS iterations ('marks'). Referee round 81: every error
    exact (`b4_end`: also as `lm_rms` evaluates it, and the end's parameters); the anchored solve without re-anchoring at stalls
    ('nostall'), the long double's anchored solve from the anchored double's end ('anchored ld'). Referee round 84: every end's exact rms
    as a decimal string and its Gauss–Newton gap (`b4_end`); the anchored double's floor predicted at its end ('floor'), at its last anchor
    ('floor_last') and from 2σ²'s rounding alone ('floor_s2') (`b4_floors`), its anchors and the 'nostall' solve's (`anchor_info`); the
    long double's anchored end with its parameters rounded to double ('rms_d': a scale, not the double parameters' floor — referee
    rounds 88, 92)."""
    r_, var, k, start, w0, c0, e_opt = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))
    rf = float(np.sqrt(np.mean(f ** 2)))
    out = dict(r=r_, var=var, k=k, start=start, e_opt=e_opt, e_start=exact_rms(f, np.r_[w0, c0], k) / rf / e_opt,
               eps=float(np.finfo(np.longdouble).eps))
    if out['eps'] >= np.finfo(float).eps:
        raise RuntimeError("B4: long double is not extended on this platform")
    for key, dt in B4_PREC:
        p, it, acc, stop = lm_fit(f, w0, c0, dt)
        out[key] = b4_end(f, p, k, rf, e_opt, it=it, acc=acc, stop=stop)
    marks, alog = dict.fromkeys(B4_MARKS), {}
    p, it, acc, stop = lm_fit(f, w0, c0, np.float64, anchored=True, marks=marks, anchor_log=alog)  # referee round 77: about an anchor
    fl = b4_floors(f, p, k, alog['p'], alog['r0'])
    out['anchored'] = b4_end(f, p, k, rf, e_opt, gap=fl['gap'], it=it, acc=acc, stop=stop, floor=fl['end'], floor_last=fl['last'],
                             floor_s2=fl['s2'], anchors=anchor_info(alog, p, it),
                             marks=[[m, exact_rms(f, q, k) / rf / e_opt] for m, q in marks.items() if q is not None])
    nlog = {}
    pn, it, acc, stop = lm_fit(f, w0, c0, np.float64, anchored=True, stall_reanchor=False, anchor_log=nlog)  # round 81: no stall re-anchoring
    out['nostall'] = b4_end(f, pn, k, rf, e_opt, it=it, acc=acc, stop=stop, anchors=anchor_info(nlog, pn, it))
    pl, it, acc, stop = lm_fit(f, p[:k], p[k:], np.longdouble, anchored=True)                # and the long double's, from that end
    out['anchored ld'] = b4_end(f, pl, k, rf, e_opt, it=it, acc=acc, stop=stop, rms_d=str(exact_resid(f, pl.astype(float), k, dec=True)[2]))
    K = RT.Bumps(X_B, SIG)                                              # the polish's own solve from the same start: TRF, B3's budget,
    wp, cp, info = RT.polish_block(K, f, np.asarray(w0, float), np.asarray(c0, float), np.ones(k), max_nfev=B3_BIG, info=True,
                                   gtol=None)                          # no gradient test
    out['trf'] = b4_end(f, np.r_[wp, cp], k, rf, e_opt, status=int(info['status']), nfev=int(info['nfev']))
    out['secs'] = time.time() - t0
    print(f"  B4 tau/sigma {r_} ({var}) K {k} from the {start}: {out['e_start']:.6g} → " + ', '.join(
        f"{key} {out[key]['ratio']:.12f} ({out[key]['it']}, {out[key]['stop']})" for key in [k for k, _ in B4_PREC] + ['anchored', 'nostall', 'anchored ld'])
          + f", TRF {out['trf']['ratio']:.12f} ({out['trf']['status']}/{out['trf']['nfev']}) [{out['secs']:.0f}s]",
          flush=True)
    return out


def b4_jobs(B, B3):
    """B4's jobs: each distinct block (width, K) of B3 whose optimum is below B4_MAX — B3's first row for it, its block in B's path (the
    block with that K whose optimum is the row's) —, from the path's block and from the stored optimum."""
    out, seen = [], set()
    K = RT.Bumps(X_B, SIG)
    for q in (B3 if B3 is not None else []):
        if q.get('k') is None or q['e_opt'] >= B4_MAX or (q['r'], q['k']) in seen:
            continue
        seen.add((q['r'], q['k']))
        b = next(b for b in B if b['r'] == q['r'] and b['var'] == q['var'])
        o = min((o for o in b['blocks'] if o['k'] == q['k']), key=lambda o: abs(o['e_opt'] / q['e_opt'] - 1))
        f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(q['r']) ** 2))
        wo = np.linalg.lstsq(K.phi(o['c_opt']), f, rcond=None)[0]
        out += [(q['r'], q['var'], q['k'], 'block', np.asarray(o['w'], float), np.asarray(o['c'], float), q['e_opt']),
                (q['r'], q['var'], q['k'], 'optimum', wo, np.asarray(o['c_opt'], float), q['e_opt'])]
    return out


def part_B4(B, B3):
    """B4: `b4_job` for `b4_jobs(B, B3)`."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b4_job, b4_jobs(B, B3), chunksize=1)


B4_ENDS = (('double', 'double'), ('long double', 'long double'), ('anchored', 'double, anchored'),
           ('nostall', 'double, anchored, no stall re-anchoring'), ('anchored ld', "long double, anchored, from the anchored double's end"),
           ('trf', 'TRF'))                                             # B4's solves (referee rounds 70, 77, 81)


def b4_best(B4):
    """B4's best end per block (τ/σ, K) (referee rounds 81, 84): the lowest exact rms any of its solves reached from either start —
    {(τ/σ, K): dict(rms, a Decimal; err, its error relative to the target's rms (the ratio times its optimum's); start; end, the solve)}."""
    from decimal import Decimal
    best = {}
    for q in (B4 if B4 is not None else []):
        for k_, _ in B4_ENDS:
            if k_ in q:
                d = Decimal(q[k_]['rms'])
                b = best.get((q['r'], q['k']))
                if b is None or d < b['rms']:
                    best[(q['r'], q['k'])] = dict(rms=d, err=q[k_]['ratio'] * q['e_opt'], start=q['start'], end=k_)
    return best


def b4_excess(d, ref):
    """d/ref − 1 for two exact rms (decimal strings or Decimals), formed in 60-digit decimal (referee round 84, D3ce-5: B4's excesses,
    not differences of double ratios)."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = 60
        return float(Decimal(d) / Decimal(ref) - 1)


def b4_lines(B4):
    """B4's table and summaries (referee rounds 70, 77, 81, 84): from the blocks near the floor and from the stored optima,
    Levenberg–Marquardt with the residual in double and in long double by direct subtraction, in double about an anchor (with and without
    re-anchoring at stalls), the long double's anchored solve from the anchored double's end, and the polish's own TRF — every end's
    error exact (`exact_resid`) and against its block's best end (`b4_best`), the excesses formed exactly (`b4_excess`); `lm_rms`'s error
    at every end; the best end's solve, its two starts' agreement and its Gauss–Newton gap; the direct subtraction's floor in double over
    its floor in long double from the same start, against ε_d/ε_ld; the anchored double's floor predicted at its end, at its last anchor
    and from 2σ²'s rounding alone, against its gap, and its anchors; the best end with its parameters rounded to double."""
    if B4 is None or len(B4) == 0:                                    # raw-loaded: an object array
        return []
    from decimal import Decimal, localcontext
    key = lambda q: (q['r'], q['k'])
    best = b4_best(B4)
    ex = lambda q, k_: b4_excess(q[k_]['rms'], best[key(q)]['rms'])
    cell = lambda q, k_: f"{ex(q, k_):.1e}, " + (f"{q[k_]['status']}/{q[k_]['nfev']}" if k_ == 'trf'
                                                  else f"{q[k_]['it']}/{q[k_]['acc']}, {q[k_]['stop']}")
    out = ['', "### B4. The floor: the same solve with the residual in double, in long double and about an anchor, evaluated exactly "
               "(referee rounds 70, 77, 81, 84)", '',
           f"B3's blocks whose optimum is below {B4_MAX:g}, from the path's block and from the stored optimum (its centres with least-squares "
           "amplitudes): Levenberg–Marquardt on the joint fit in (w, c) (`lm_fit`: the Jacobian and the step in double; a step accepted when "
           "its predicted decrease is positive and the objective, as the solve computes its residual, decreases) with the parameters and the "
           f"residual f − Φ(c)w by direct subtraction in double or in long double (its ε {B4[0]['eps']:.3g}), or about an anchor (w₀, c₀) — "
           "r = r₀ − Σ_kφ₀_k[(w_k − w₀_k) + w_k·expm1(Δ_k)], Δ_k = (c_k − c₀_k)(2x − c_k − c₀_k)/(2σ²), r₀ by direct subtraction: in double "
           f"re-anchored at each stall (μ reset) and every {B4_ANCHOR:g} accepted steps, or not at stalls, and in long double from the "
           f"anchored double's end —, and the polish's own TRF (`polish_block`, {B3_BIG:g} evaluations, no gradient test). Every error is "
           f"exact (`exact_resid`: the residual in {EXACT_PREC}-digit decimal arithmetic, every double and long double taken exactly); the "
           "start's error over the stored optimum's; the block's best end — the lowest error any solve reached from either start — over it; "
           "and each end's excess over that best (formed exactly, `b4_excess`), its iterations/accepted steps and stop (stall: no step the "
           f"solve computes decreases its objective; flat: the exact rms fell by less than {B4_FLAT:g} of itself over {B4_ANCHOR:g} accepted "
           f"steps; cap: {B4_CAP:g} iterations) or TRF's status/evaluations; † the optimum below {E_FLOOR:g}:", '',
           "| τ/σ | variant | K | optimum | start | its error | the best end | " + ' | '.join(lab for _, lab in B4_ENDS) + ' |',
           '|---|---|---|---|---|---|---|' + '---|' * len(B4_ENDS)]
    for q in B4:
        out.append(f"| {q['r']:g} | {q['var']} | {q['k']} | {q['e_opt']:.2e}{'†' if q['e_opt'] < E_FLOOR else ''} | {q['start']} | "
                   f"{q['e_start']:.6g} | {best[key(q)]['err'] / q['e_opt']:.12f} | " + ' | '.join(cell(q, k_) for k_, _ in B4_ENDS) + ' |')
    blk = sorted((q for q in B4 if q['start'] == 'block'), key=lambda q: -q['e_opt'])
    opt = {key(q): q for q in B4 if q['start'] == 'optimum'}
    both = [q for q in blk if key(q) in opt]
    lab = lambda q: f"{q['e_opt']:.3g} (τ/σ = {q['r']:g}, K = {q['k']})"
    rng = lambda v: f"{min(v):+.1e} … {max(v):+.1e}"
    pos = lambda v: f"{min(v):.1e} … {max(v):.1e}"
    ev = [q[k_]['ratio_ld'] / q[k_]['ratio'] - 1 for q in B4 for k_, _ in B4_ENDS]
    evl = [q['long double']['ratio_ld'] / q['long double']['ratio'] - 1 for q in B4]
    out += ['', f"The long-double rms by direct subtraction (`lm_rms`, B4's instrument until round 81) against the exact one, minus 1: "
                f"{rng(ev)} over every end, {rng(evl)} at the long-double solve's ends ({sum(v < 0 for v in evl)} of {len(evl)} below it: that "
                "solve minimised the same evaluation)."]
    per = lambda k_: [ex(q, k_) for q in B4]
    out += ['', "Each solve's excess over the best end, minus 1, from both starts: " + '; '.join(
        f"{lab_} {rng(per(k_))}" for k_, lab_ in B4_ENDS) + ". By the optimum, from the path's blocks — " + '; '.join(
        f"{lab(q)}: double {ex(q, 'double'):.1e}, long double {ex(q, 'long double'):.1e}" for q in blk) + '.']
    names = dict(B4_ENDS)
    cells = [key(q) for q in blk]
    who = [best[c_]['end'] for c_ in cells]
    spread = {c_: abs(b4_excess(q[best[c_]['end']]['rms'], opt[c_][best[c_]['end']]['rms'])) for q in blk for c_ in [key(q)] if c_ in opt}
    gaps = [next(q for q in B4 if key(q) == c_ and q['start'] == best[c_]['start'])[best[c_]['end']]['gap'] for c_ in cells]
    out += ['', "The best end is reached by " + ' and by '.join(f"{names[e_]} at {who.count(e_)}" for e_ in dict.fromkeys(who))
            + f" of the {len(cells)} blocks" + (f"; that solve's ends from the two starts agree within {max(spread.values()):.1e} (exactly)"
                                                 if spread else '')
            + f", and a Gauss–Newton step from the best end's exact residual predicts at most {max(gaps):.1e} of further decrease "
            "(`b4_gap`: ‖P_Jr‖²/(2‖r‖²), J and r exact)."]
    eg = [ex(q, k_) / q[k_]['gap'] for q in B4 for k_, _ in B4_ENDS if ex(q, k_) >= 1e3 * max(gaps)]
    if eg:
        out[-1] += (f" The {len(eg)} of the {len(B4) * len(B4_ENDS)} ends whose excess is at least 10³ times that lie at their own Gauss–Newton "
                    f"gap from it: their excess over their gap is {min(eg):.3f}–{max(eg):.3f} (median {np.median(eg):.3f}) — the best end the "
                    "minimum the solves approach, each excess a first-order distance to it.")
    eps_d, eps_l = float(np.finfo(float).eps), float(B4[0]['eps'])
    dl = [(q, ex(q, 'double'), ex(q, 'long double')) for q in B4]
    use = [(q, a / b) for q, a, b in dl if b > spread.get(key(q), 0.0)]
    oth = [(q, b) for q, a, b in dl if b <= spread.get(key(q), 0.0)]
    if use:
        rs = sorted(x for _, x in use)
        ex_ = [np.log(x) / np.log(eps_d / eps_l) for x in rs]
        out += ['', "The direct subtraction's floor in double over its floor in long double, from the same start (the two ends' excesses over "
                    f"the best end), where the long double's exceeds that spread ({len(use)} of {len(B4)} rows): " + ', '.join(f"{x:.2g}" for x in rs)
                + f" — median {np.median(rs):.3g}, as (ε_d/ε_ld)^a with a = {min(ex_):.2f}–{max(ex_):.2f} (median {np.median(ex_):.2f}); a floor "
                  f"linear in the rounding (a = 1) gives ε_d/ε_ld = {eps_d / eps_l:.4g}"
                + ('; at the other ' + ', '.join(f"{lab(q)} from the {q['start']}: {b:.1e}" for q, b in oth) if oth else '') + '.']
    A_ = [q['anchored'] for q in B4]
    if all('floor_last' in a for a in A_):
        with localcontext() as ctx:
            ctx.prec = 60
            d2 = float((Decimal(float(2 * np.float64(SIG) ** 2)) - 2 * Decimal(SIG) ** 2) / (2 * Decimal(SIG) ** 2))
        rl = [a['floor_last'] / a['gap'] for a in A_]
        re_ = [a['floor'] / a['gap'] for a in A_]
        rs2 = [a['floor_s2'] / a['floor_last'] for a in A_]
        an = [a['anchors'] for a in A_]
        fl_ = [x for a, x in zip(A_, an) if a['stop'] == 'flat']
        ns = [q['nostall']['anchors']['period'] for q in blk]
        out += ['', "The anchored double's floor, the rise of the exact rms that a rounding e₀ frozen in its residual predicts, ‖P_Je₀‖²/(2‖r‖²) "
                    "at its end (`b4_floors`, exactly): with e₀ at the end " + pos([a['floor'] for a in A_]) + "; with e₀ at its last anchor, the "
                    "rounding it froze, " + pos([a['floor_last'] for a in A_]) + "; its gap — the decrease a Gauss–Newton step predicts there — "
                + pos([a['gap'] for a in A_]) + f": the last anchor's prediction {min(rl):.2f}–{max(rl):.2f} times the gap, the end's "
                  f"{min(re_):.2f}–{max(re_):.2f}; 2σ² rounded to double (relative {d2:.3g}), the same at every anchor, alone predicts "
                  f"{min(rs2):.2f}–{max(rs2):.2f} of the last anchor's. Its last anchor: " + ', '.join(
                      f"{sum(x['kind'] == k_ for x in an)} {k_}" for k_ in ('stall', 'period', 'start') if any(x['kind'] == k_ for x in an))
                + f", {min(x['before'] for x in an)}–{max(x['before'] for x in an)} iterations before the end"
                + (f" ({min(x['dist'] for x in fl_):.1e} … {max(x['dist'] for x in fl_):.1e} from it at the flat stops)" if fl_ else '')
                + f"; per solve {min(x['stall'] for x in an)}–{max(x['stall'] for x in an)} re-anchors at stalls and "
                  f"{min(x['period'] for x in an)}–{max(x['period'] for x in an)} by period; without re-anchoring at stalls, from the blocks, "
                  f"{min(ns)}–{max(ns)} by period."]
    if all('rms_d' in q['anchored ld'] for q in B4):
        rd_ = {}
        for c_ in cells:
            qb = next(q for q in B4 if key(q) == c_ and q['start'] == best[c_]['start'])
            rd_[c_] = (b4_excess(qb['anchored ld']['rms_d'], best[c_]['rms']), [ex(q, 'anchored') for q in B4 if key(q) == c_])
        out += ['', "The long double's anchored end of the best end's start with its parameters rounded to double, over the best end, minus 1 "
                    "(a scale, not the double parameters' floor: B5; the anchored double's excesses in brackets): " + '; '.join(
                        f"{lab(q)}: {rd_[key(q)][0]:.1e} (" + ', '.join(f"{x:.1e}" for x in rd_[key(q)][1]) + ')' for q in blk) + '.']
    if both:
        agree = [abs(b4_excess(q['anchored']['rms'], opt[key(q)]['anchored']['rms'])) for q in both]
        out += ['', "The anchored double's ends from the block and from the stored optimum agree within " f"{max(agree):.1e}; its stops: "
                + ', '.join(f"{k} {sum(q['anchored']['stop'] == k for q in B4)}" for k in ('flat', 'stall', 'cap')
                            if any(q['anchored']['stop'] == k for q in B4))
                + f", after {min(q['anchored']['it'] for q in B4)} … {max(q['anchored']['it'] for q in B4)} iterations; without re-anchoring "
                f"at stalls it stops after {min(q['nostall']['it'] for q in B4)} … {max(q['nostall']['it'] for q in B4)}. The long double's "
                "anchored solve from the anchored double's end gains on it, minus 1, "
                + rng([b4_excess(q['anchored ld']['rms'], q['anchored']['rms']) for q in B4]) + ". The stored optimum over the best end, minus 1: "
                + '; '.join(f"{lab(q)}: {opt[key(q)]['e_start'] * q['e_opt'] / best[key(q)]['err'] - 1:+.1e}" for q in both) + '.']
    if any(q['anchored'].get('marks') for q in B4):
        out += ['', "The anchored double's error after " + ', '.join(f"{m:g}" for m in B4_MARKS) + " iterations, over the best end, minus 1: "
                + '; '.join(f"{lab(q)} from the {q['start']}: " + ', '.join(f"{x * q['e_opt'] / best[key(q)]['err'] - 1:+.1e}" for _, x in
                                                                           q['anchored']['marks']) for q in B4 if q['anchored'].get('marks')) + '.']
    out.append('')
    return out


B5_CELLS = ((1.0, 12), (0.25, 6), (0.75, 11))                          # B5 (referee round 84, D3ce-1): the κ-scan's blocks (τ/σ, K)
B5_KAPPA = (0.0, 1 / 16, 1 / 8, 1 / 4, 1 / 2, 1.0, 2.0, 4.0)           # B5: the direct subtraction's rounding scaled by κ
B5_KAPPA_ANC = (0.0, 1 / 4, 1 / 2, 1.0, 2.0, 4.0)                      # B5: the anchor's, from the stored optimum
B5_THETA = (0.0, 0.25, 0.5, 0.75, 1.0)                                 # B5's starts: the stored optimum (0), the path's block (1), between


def b5_jobs(B, B3):
    """B5's jobs (referee round 84): for each block of B5_CELLS among B4's (`b4_jobs`), the direct solve from B4's two starts and from
    the starts between them — the centres c_opt + θ(c_block − c_opt) with least-squares amplitudes, θ in B5_THETA (θ = 0 the stored
    optimum's start, θ = 1 the path's block, as in B4) —, and the anchored solve from the stored optimum's start."""
    K = RT.Bumps(X_B, SIG)
    jobs = {(j[0], j[2], j[3]): j for j in b4_jobs(B, B3)}
    out = []
    for r_, k in B5_CELLS:
        if (r_, k, 'block') not in jobs or (r_, k, 'optimum') not in jobs:
            continue
        jb, jo = jobs[(r_, k, 'block')], jobs[(r_, k, 'optimum')]
        f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))
        for th in B5_THETA:
            if th == 0:
                w0, c0 = jo[4], jo[5]
            elif th == 1:
                w0, c0 = jb[4], jb[5]
            else:
                c0 = jo[5] + th * (jb[5] - jo[5])
                w0 = np.linalg.lstsq(K.phi(c0), f, rcond=None)[0]
            out.append((r_, jo[1], k, th, 'direct', w0, c0, jo[6]))
        out.append((r_, jo[1], k, 0.0, 'anchored', jo[4], jo[5], jo[6]))
    return out


def b5_job(args):
    """B5 (referee round 84, D3ce-1): B4's solve in double from one start with the direct subtraction's rounding scaled by κ (`lm_fit`'s
    kappa) — in the residual ('direct', κ in B5_KAPPA) or in the anchor's r₀ ('anchored', κ in B5_KAPPA_ANC, B4's re-anchoring and flat
    stop) —: each end's exact rms (a decimal string), its ratio over the optimum's, iterations, accepted steps, stop and parameters."""
    r_, var, k, th, mode, w0, c0, e_opt = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))
    rf = float(np.sqrt(np.mean(f ** 2)))
    rows = []
    for kap in (B5_KAPPA if mode == 'direct' else B5_KAPPA_ANC):
        p, it, acc, stop = lm_fit(f, w0, c0, np.float64, anchored=mode == 'anchored', kappa=kap)
        _, e, d = exact_resid(f, p, k, dec=True)
        rows.append(dict(kappa=kap, rms=str(d), ratio=e / rf / e_opt, it=it, acc=acc, stop=stop, p=p.tolist()))
    out = dict(r=r_, var=var, k=k, theta=th, mode=mode, e_opt=e_opt, e_start=exact_rms(f, np.r_[w0, c0], k) / rf / e_opt, rows=rows,
               secs=time.time() - t0)
    print(f"  B5 tau/sigma {r_} K {k} theta {th:g} {mode}: " + ', '.join(f"{x['kappa']:g}: {x['ratio']:.12f} ({x['it']}, {x['stop']})"
                                                                       for x in rows) + f" [{out['secs']:.0f}s]", flush=True)
    return out


def part_B5(B, B3):
    """B5: `b5_job` for `b5_jobs(B, B3)`."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(b5_job, b5_jobs(B, B3), chunksize=1)


def log_slope(xs, ys):
    """The least-squares slope of log y against log x over the pairs with x, y > 0 (None with fewer than two)."""
    pts = [(np.log(x), np.log(y)) for x, y in zip(xs, ys) if x > 0 and y > 0]
    if len(pts) < 2:
        return None
    a, b = np.array(pts).T
    return float(np.polyfit(a, b, 1)[0])


def b5_lines(B5, B4):
    """B5's tables (referee round 84, D3ce-1): each end's excess over B4's best end of its block (exact, `b4_excess`) by start and κ; the
    slope of log(excess) against log κ (`log_slope`) over κ in [1/16, 1] and in [1/4, 4], per start and for the geometric mean over the
    starts; the anchored solve's excess by κ, its slope over [1/2, 4], and at κ = 0 against B4's best end with its parameters rounded to
    double; whether the ends at κ = 1 are B4's."""
    if B5 is None or len(B5) == 0 or B4 is None or len(B4) == 0:
        return []
    best = b4_best(B4)
    fr = lambda k_: '0' if k_ == 0 else (f"1/{round(1 / k_)}" if k_ < 1 else f"{k_:g}")
    ex = lambda x, y: b4_excess(y['rms'], best[(x['r'], x['k'])]['rms'])
    sl = lambda ks, es, lo, hi: log_slope([k_ for k_ in ks if lo <= k_ <= hi], [e for k_, e in zip(ks, es) if lo <= k_ <= hi])
    fs = lambda s: '—' if s is None else f"{s:.2f}"
    out = ['', "### B5. The floors' order: the double's rounding scaled (referee round 84)", '',
           "B4's solve in double at three blocks with the residual's rounding scaled by κ — r_κ = r_ld + κ(r_d − r_ld), formed in long double "
           "from the direct subtractions in double (r_d) and in long double (r_ld) at the same double parameters and rounded to double "
           "(`kappa_blend`; κ = 1 the double's residual exactly, κ = 0 the long double's) —, by direct subtraction from five starts (the stored "
           "optimum's, θ = 0; the path's block, θ = 1; the centres c_opt + θ(c_block − c_opt) with least-squares amplitudes between) and, "
           "from the stored optimum's, about an anchor whose r₀ is r_κ (B4's re-anchoring and flat stop): each end's excess over B4's best "
           "end of its block, minus 1 (exact), and the slope of log(excess) against log κ (least squares) over κ in [1/16, 1] and [1/4, 4]:", '',
           '| τ/σ | K | start θ (its error) | ' + ' | '.join(f"κ = {fr(k_)}" for k_ in B5_KAPPA) + ' | slope on [1/16, 1] | on [1/4, 4] |',
           '|---|---|---|' + '---|' * len(B5_KAPPA) + '---|---|']
    dirs = [x for x in B5 if x['mode'] == 'direct']
    for x in dirs:
        ks, es = [y['kappa'] for y in x['rows']], [ex(x, y) for y in x['rows']]
        out.append(f"| {x['r']:g} | {x['k']} | {x['theta']:g} ({x['e_start']:.4g}) | "
                   + ' | '.join(f"{e:.1e}" + ('*' if y['acc'] <= 2 else '') for e, y in zip(es, x['rows']))
                   + f" | {fs(sl(ks, es, 1 / 16, 1))} | {fs(sl(ks, es, 1 / 4, 4))} |")
    if any(y['acc'] <= 2 for x in dirs for y in x['rows']):              # referee round 88 (D3ci-5c)
        out += ['', "* at most two accepted steps: the end is the start's, censored (the slopes include it)."]
    gm = []
    for c_ in dict.fromkeys((x['r'], x['k']) for x in dirs):
        xs = [x for x in dirs if (x['r'], x['k']) == c_]
        ks = [y['kappa'] for y in xs[0]['rows']]
        es = [[ex(x, y) for y in x['rows']] for x in xs]
        g = [float(np.exp(np.mean(np.log(col)))) if all(v > 0 for v in col) else 0.0 for col in zip(*es)]
        gm.append((c_, len(xs), g, sl(ks, g, 1 / 16, 1), sl(ks, g, 1 / 4, 4), [max(col) / min(col) if min(col) > 0 else np.inf
                                                                                 for col in zip(*es)]))
    out += ['', "The geometric mean over the starts: " + '; '.join(
        f"τ/σ = {c_[0]:g}, K = {c_[1]} ({n} starts): " + ', '.join(f"{fr(k_)}: {v:.1e}" for k_, v in zip(B5_KAPPA, g))
        + f" — slopes {fs(a)} on [1/16, 1], {fs(b)} on [1/4, 4]; the starts' largest over least at one κ {min(sp):.3g}–{max(sp):.3g}"
        for c_, n, g, a, b, sp in gm) + '.']
    per = [[sl([y['kappa'] for y in x['rows']], [ex(x, y) for y in x['rows']], lo, hi) for x in dirs] for lo, hi in ((1 / 16, 1), (1 / 4, 4))]
    per = [[s for s in p_ if s is not None] for p_ in per]
    sp_all = [(v, c_, k_) for c_, n, g, a, b, sp in gm for v, k_ in zip(sp, B5_KAPPA)]
    if per[0] and per[1] and sp_all:
        v, c_, k_ = max(sp_all, key=lambda z: z[0])
        out += ['', f"Per start ({len(dirs)}): the slopes " + ', '.join(f"{min(p_):.2f}–{max(p_):.2f} on {iv} (median {float(np.median(p_)):.2f})"
                                                                    for p_, iv in zip(per, ('[1/16, 1]', '[1/4, 4]')))
                + f"; the starts' largest over least at one κ at most {v:.3g} (τ/σ = {c_[0]:g}, K = {c_[1]}, κ = {fr(k_)})."]
    anc = [x for x in B5 if x['mode'] == 'anchored']
    if anc:
        out += ['', '| τ/σ | K | about an anchor, from the stored optimum: ' + ' | '.join(f"κ = {fr(k_)}" for k_ in B5_KAPPA_ANC)
                + ' | slope on [1/2, 4] | over the best end rounded to double, κ = ' + ', '.join(fr(k_) for k_ in B5_KAPPA_ANC) + ' |',
                '|---|---|---|' + '---|' * (len(B5_KAPPA_ANC) - 1) + '---|---|']
        sa = []
        for x in anc:
            ks, es = [y['kappa'] for y in x['rows']], [ex(x, y) for y in x['rows']]
            c_ = (x['r'], x['k'])
            qb = next((q for q in B4 if (q['r'], q['k']) == c_ and q['start'] == best[c_]['start']), None)
            r0 = (', '.join(f"{e / b4_excess(qb['anchored ld']['rms_d'], best[c_]['rms']):.3g}" for e in es)
                  if qb is not None and 'rms_d' in qb['anchored ld'] else '—')
            sa.append(sl(ks, es, 1 / 2, 4))
            out.append(f"| {x['r']:g} | {x['k']} | " + ' | '.join(f"{e:.1e} ({y['it']}, {y['stop']})" for e, y in zip(es, x['rows']))
                       + f" | {fs(sa[-1])} | {r0} |")
        sa = [s for s in sa if s is not None]
        out += ['', "The anchored solves' slopes on [1/2, 4] " + (f"{min(sa):.2f}–{max(sa):.2f}" if sa else '—') + "; their flat stops by κ: "
                + ', '.join(f"{fr(k_)}: {sum(1 for x in anc for y in x['rows'] if y['kappa'] == k_ and y['stop'] == 'flat')} of "
                            f"{sum(1 for x in anc for y in x['rows'] if y['kappa'] == k_)}" for k_ in B5_KAPPA_ANC) + '.']
    same, n_ = 0, 0
    for x in B5:
        y = next((y for y in x['rows'] if y['kappa'] == 1.0), None)
        st = {0.0: 'optimum', 1.0: 'block'}.get(x['theta'])
        q = next((q for q in B4 if (q['r'], q['k'], q['start']) == (x['r'], x['k'], st)), None)
        if y is None or q is None or (x['mode'] == 'anchored' and st != 'optimum'):
            continue
        z = q['double' if x['mode'] == 'direct' else 'anchored']
        n_ += 1
        same += y['rms'] == z['rms'] and y['it'] == z['it'] and y['acc'] == z['acc'] and y['stop'] == z['stop'] and list(y['p']) == list(z['p'])
    out += ['', f"At κ = 1 the ends are B4's (the double's from the stored optimum and the block, the anchored double's from the stored optimum) "
                f"in {same} of {n_}: rms, iterations, accepted steps, stop and parameters.", '']
    return out


def gh_defect_v(n, K, v):
    """D_{n,K}(v) = E He_n(√v Z) − Σ_k w_k He_n(√v z_k), exactly (v rational, n even; referee round 65): He_n(√v z) = Σ_m c_m v^{m/2}z^m
    reduced modulo He_K(z) — the rule integrates the remainder, of degree < K, exactly, and He_K vanishes at its nodes —, and
    E He_n(√v Z) = (v − 1)^{n/2}(n − 1)!!. At v = 2 it is `wide_feature.gh_defect`."""
    v = Fraction(v)
    c = W._he(n)
    p = [Fraction(c[m]) * v ** (m // 2) if m % 2 == 0 else Fraction(0) for m in range(len(c))]
    q = W._he(K)
    for d in range(len(p) - 1, K - 1, -1):
        f_ = p[d]
        if f_:
            for i in range(K + 1):
                p[d - K + i] -= f_ * q[i]
    QK = sum(p[m] * W._gauss_moment(m) for m in range(min(K, len(p))))
    return (v - 1) ** (n // 2) * W._gauss_moment(n) - QK


def frame(rho, sig, tau2):
    """The Mehler frame of ratio ρ of the kernel e^{−(y−y′)²/(4σ²)} (referee round 65): k(y, y′) = g(y)g(y′)Σ_n (ρⁿ/n!)He_n(y/ℓ)He_n(y′/ℓ)
    with ℓ² = 2σ²ρ/(1 − ρ²); the Gauss–Hermite rule of the frame's own Gaussian N(0, v_ρ) — g times N(0, τ²) —, v_ρ = 2σ²τ²/(2σ² +
    (1 − ρ)τ²), has no defect below 2K in it; the measure's variance in the frame's units v′ = v_ρ/ℓ². Returns (v_ρ, v′)."""
    v = 2 * sig ** 2 * tau2 / (2 * sig ** 2 + (1 - rho) * tau2)
    return v, v / (2 * sig ** 2 * rho / (1 - rho ** 2))


A5_KMAX, A5_PREC = 40, 110                                            # A5 (referee round 81): the frame series exactly to K = A5_KMAX
A5_FAR, A5_FAR_PREC = ((60, 80, 100, 120), (0.5, 1.0, 2.0)), 140       # A5: the Mercer frame's rule further, at these K and τ/σ
A5_FAR_Q, A5_FAR_CHECK = (40, 60, 80, 100, 120, 160), 40              # A5 (referee round 84, D3ce-3): at τ/σ = 1/4 at these K in 4K + 60
                                                                       # digits; each τ/σ's last K again with A5_FAR_CHECK more digits


def he_all(x, nmax):
    """He_0 … He_nmax at x (He_{n+1} = xHe_n − nHe_{n−1})."""
    h = [x * 0 + 1, x]
    for n in range(1, nmax):
        h.append(x * h[n] - n * h[n - 1])
    return h


def gh_nodes(K, prec=A5_PREC):
    """The K-point Gauss–Hermite rule of N(0, 1) in prec-digit Decimal (referee round 81): numpy's nodes refined by Newton on He_K
    (He_K′ = KHe_{K−1}), the weights K!/(K²He_{K−1}(z_k)²)."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        zs, ws = [], []
        for z in np.polynomial.hermite_e.hermegauss(K)[0]:
            x_ = Decimal(float(z))
            for _ in range(8):
                h = he_all(x_, K)
                x_ -= h[K] / (K * h[K - 1])
            h = he_all(x_, K)
            zs.append(x_)
            ws.append(Decimal(factorial(K)) / (K * K * h[K - 1] ** 2))
        return zs, ws


def frame_series(K, rho, v, prec=A5_PREC, tol_exp=40):
    """A K-point rule's error series in a Mehler frame, exactly (referee round 81): S(K) = Σ_{n ≥ 2K even} (ρⁿ/n!)D_n(v)², D_n(v) =
    E He_n(√v Z) − Σ_k w_k He_n(√v z_k) the defect on He_n of the Gauss–Hermite rule of N(0, v) in the frame's units (E He_n(√v Z) =
    (v − 1)^{n/2}(n − 1)!!), summed until the last two terms fall below 10^{−tol_exp} of the sum, one of them nonzero (referee round 84,
    D3ce-6: at K = 1 in the half-squeeze frame D_n = 0 for n ≡ 0 mod 4); ρ and v Decimals. Returns (S, the terms, the largest defect below
    2K relative to its parts — zero but for the rounding where the rule is exact; it bounds that rounding, not the tail's)."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        zs, ws = gh_nodes(K, prec)
        xs = [v.sqrt() * z for z in zs]
        h_prev, h_cur = [Decimal(0)] * K, [Decimal(1)] * K
        S, terms, low = Decimal(0), [], Decimal(0)
        dfact, nfact, rhon, tol, n = Decimal(1), Decimal(1), Decimal(1), Decimal(10) ** -tol_exp, 0
        while True:
            if n % 2 == 0:
                E_ = (v - 1) ** (n // 2) * dfact
                Q = sum(w * h for w, h in zip(ws, h_cur))
                D = E_ - Q
                if n < 2 * K:
                    low = max(low, abs(D) / (abs(E_) + abs(Q) + 1))
                else:
                    terms.append(rhon / nfact * D * D)
                    S += terms[-1]
                    if len(terms) > 10 and 0 < max(terms[-2:]) < tol * S:
                        return S, terms, low
            h_prev, h_cur = h_cur, [x_ * hc - n * hp for x_, hc, hp in zip(xs, h_cur, h_prev)]
            n += 1
            nfact *= n
            rhon *= rho
            if n % 2 == 0:
                dfact *= n - 1


def a5_series(r_):
    """A5's exact series at τ/σ = r_ (referee round 81): the Mercer frame's rule (ρ = q, v′ = 1 + q) and the half-squeeze rule (ρ = t,
    v′ = 2), K = 1 … A5_KMAX: S(K), the first term's share, the tail over the first term, and the largest low defect."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = A5_PREC
        b = (1 + 2 * Decimal(r_) ** 2).sqrt()
        q = (b - 1) / (b + 1)
        st = (1 + Decimal(r_) ** 2).sqrt()                              # σ*/σ
        t = (st - 1) / (st + 1)
        out = {}
        for nm, rho, v in (('Mercer', q, 1 + q), ('half-squeeze', t, Decimal(2))):
            rows = []
            for K in range(1, A5_KMAX + 1):
                S_, terms, low = frame_series(K, rho, v)
                rows.append(dict(K=K, S=float(S_), first=float(terms[0] / S_), tail=float((S_ - terms[0]) / terms[0]), low=float(low)))
            out[nm] = rows
        return dict(r=r_, q=float(q), t=float(t), frames=out)


def a5_far_row(r_, K, prec):
    """The Mercer frame's rule's ratio per kernel (S(K − 1)/S(K))^{1/2} at τ/σ = r_ in prec-digit decimal arithmetic, with the largest low
    defect (referee rounds 81, 84)."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        b = (1 + 2 * Decimal(r_) ** 2).sqrt()
        q = (b - 1) / (b + 1)
        S0 = frame_series(K - 1, q, 1 + q, prec=prec)
        S1 = frame_series(K, q, 1 + q, prec=prec)
        return dict(r=r_, K=K, prec=prec, ratio=float((S0[0] / S1[0]).sqrt()), low=float(max(S0[2], S1[2])))


def a5_far():
    """The Mercer frame's rule's ratio per kernel further out (`a5_far_row`; referee rounds 81, 84): at τ/σ in A5_FAR[1] and K in A5_FAR[0]
    in A5_FAR_PREC digits, at τ/σ = 1/4 and K in A5_FAR_Q in 4K + 60 (its cancellation grows like K digits); each τ/σ's last K again
    with A5_FAR_CHECK more digits ('check')."""
    out = [a5_far_row(r_, K, A5_FAR_PREC) for r_ in A5_FAR[1] for K in A5_FAR[0]] + [a5_far_row(0.25, K, 4 * K + 60) for K in A5_FAR_Q]
    last = {y['r']: y for y in out}
    for y in last.values():
        y['check'] = a5_far_row(y['r'], y['K'], y['prec'] + A5_FAR_CHECK)['ratio']
    return out


def part_A5(A, series=True, far=True):
    """A5 (referee round 65): per τ/σ the frames at ρ = t (the half-squeeze rule), q (the Mercer expansion) and λ (v′ = 1): v_ρ over the
    half-squeeze variance and v′; the Mercer frame's rule's (N(0, 2τ²/(β² + 1))) errors at A's K; and the first tail ratio in a frame of
    variance ratio v′, D_{2K+2}/D_{2K} against −v′K(K + 2) + (v′ − 1)(K + 1)(2K + 1), exactly; the exact series (`a5_series`) and the
    Mercer frame's rule further out (`a5_far`) unless series, far are false."""
    rows = []
    for r_ in TAUS_A:
        ss = sigma_star(r_)
        P = W.params(SIG, ss)
        tau2 = P['tau'] ** 2
        b2 = np.sqrt(1 + 2 * tau2 / SIG ** 2)
        q = (b2 - 1) / (b2 + 1)
        vhs = tau2 * SIG / ss
        fr = {nm: (lambda vv: dict(rho=rho, v=vv[0] / vhs, vp=vv[1]))(frame(rho, SIG, tau2)) for nm, rho in (('t', P['t']), ('q', q),
                                                                                                          ('λ', P['lam']))}
        Ks = [x['K'] for x in A['rows'] if x['r'] == r_]
        vM = 2 * tau2 / (b2 + 1)
        rows.append(dict(r=r_, t=P['t'], q=q, frames=fr, vM=vM / vhs, Ks=Ks,
                         e_mframe=[W.l2_error(*W.scaled_rule(K, SIG, ss, vM), SIG, ss) for K in Ks]))
    chk = [dict(v=str(v), K=K, ok=bool(gh_defect_v(2 * K + 2, K, v) / gh_defect_v(2 * K, K, v)
                                       == -v * K * (K + 2) + (v - 1) * (K + 1) * (2 * K + 1)))
           for v in (Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3), Fraction(7, 4)) for K in range(1, 13)]
    return dict(rows=rows, check=chk, series=[a5_series(r_) for r_ in TAUS_A] if series else [], far=a5_far() if far else [])


def a5_lines(A5, A=None):
    """A5's table, its exact frame series (`a5_series_lines`) and the tail-ratio check."""
    if A5 is None:
        return []
    out = ['', '### A5. The frames: which Gauss–Hermite rule a Mehler expansion singles out', '',
           "The Mehler frames of the kernel e^{−(y−y′)²/(4σ²)} (ratio ρ, Hermite scale ℓ² = 2σ²ρ/(1 − ρ²)): the Gauss–Hermite rule of the "
           "frame's own Gaussian N(0, v_ρ), v_ρ = 2σ²τ²/(2σ² + (1 − ρ)τ²), has no defect below 2K in the frame; v′ = v_ρ/ℓ² the measure's "
           "variance in the frame's units. At ρ = t (the half-squeeze rule), q (the Mercer expansion) and λ: v_ρ over the half-squeeze "
           "variance, and v′; the Mercer frame's rule N(0, 2τ²/(β² + 1)) — not the Mercer-zero rule — and its ratio per kernel at the last "
           "three K against 1/t and 1/q:", '',
           '| τ/σ | ρ = t: v/v_hs, v′ | ρ = q: v/v_hs, v′ | ρ = λ: v/v_hs, v′ | the Mercer frame\'s rule: ratio per kernel (1/t, 1/q) |',
           '|---|---|---|---|---|']
    for x in A5['rows']:
        f_ = x['frames']
        e = x['e_mframe']
        out.append(f"| {x['r']:g} | " + ' | '.join(f"{f_[k]['v']:.6f}, {f_[k]['vp']:.6f}" for k in ('t', 'q', 'λ'))
                   + f" | {', '.join(f'{e[i - 1] / e[i]:.4f}' for i in range(len(e) - 3, len(e)))} at K = {x['Ks'][-3]}–{x['Ks'][-1]} "
                   f"({1 / x['t']:.4g}, {1 / x['q']:.4g}) |")
    out += a5_series_lines(A5, A)
    ok = [c_ for c_ in A5['check'] if c_['ok']]
    out += ['', f"The first tail ratio D_{{2K+2}}/D_{{2K}} = −v′K(K + 2) + (v′ − 1)(K + 1)(2K + 1) (`gh_defect_v`, exact) at "
                f"v′ = {', '.join(sorted({c_['v'] for c_ in A5['check']}, key=lambda s: Fraction(s)))}, K = 1 … "
                f"{max(c_['K'] for c_ in A5['check'])}: " + (f"exact in all {len(ok)}" if len(ok) == len(A5['check'])
                                                              else f"FAILS in {len(A5['check']) - len(ok)} of {len(A5['check'])}")
            + " — −(K − 1) at v′ = 2, growing like (v′ − 2)K² at any other v′.", '']
    return out


def a5_series_lines(A5, A=None):
    """A5's exact frame series (referee round 81): per τ/σ and frame the tail over the first term at K = 5, 10, 20, 40 and the first K from
    which it exceeds 1, the ratio per kernel (S(K − 1)/S(K))^{1/2} at K = 10, 20, 30, 40 against 1/t and 1/q, and the series' ratio per
    kernel against the actual errors' (the Mercer frame's rule's of A5's table, the half-squeeze rule's of A) where both exist."""
    if A5 is None or not A5.get('series'):
        return []
    rows = {x['r']: x for x in A5['rows']}
    act = {}
    for x in (A['rows'] if A is not None else []):
        act.setdefault(x['r'], {})[x['K']] = x['e_hs']
    out = ['', f"The rules' error series in their frames, exactly (`frame_series`: {A5_PREC} digits, K = 1 … {A5_KMAX}; referee round 81) — the "
               "Mercer frame's rule (ρ = q, v′ = 1 + q) and the half-squeeze rule (ρ = t, v′ = 2): the tail after the first term over the "
               "first at K = 5, 10, 20, 40 and the first K from which it exceeds it; the ratio per kernel (S(K − 1)/S(K))^{1/2} at K = 10, "
               "20, 30, 40 against 1/t and 1/q; and the series' ratio per kernel against the actual errors' where both are tabulated (the "
               "largest relative difference):", '',
           '| τ/σ | frame | tail/first at K = 5, 10, 20, 40 | the tail above the first from K | ratio per kernel at K = 10, 20, 30, 40 (1/t, 1/q) '
           '| against the actual errors |', '|---|---|---|---|---|---|']
    for x in A5['series']:
        for nm, rs in x['frames'].items():
            by = {y['K']: y for y in rs}
            rat = lambda K: float(np.sqrt(by[K - 1]['S'] / by[K]['S']))
            over = next((y['K'] for y in rs if y['tail'] > 1), None)
            if nm == 'Mercer':
                e = rows.get(x['r'], {}).get('e_mframe') or []
                ks = rows.get(x['r'], {}).get('Ks') or []
                ea = dict(zip(ks, e))
            else:
                ea = act.get(x['r'], {})
            dif = [abs(rat(K) / (ea[K - 1] / ea[K]) - 1) for K in sorted(ea) if K - 1 in ea and K in by and K > 1]
            out.append(f"| {x['r']:g} | {nm} | " + ', '.join(f"{by[K]['tail']:.3g}" for K in (5, 10, 20, 40)) + ' | '
                       + ('—' if over is None else str(over)) + ' | ' + ', '.join(f"{rat(K):.4f}" for K in (10, 20, 30, 40))
                       + f" ({1 / x['t']:.4g}, {1 / x['q']:.4g}) | " + (f"{max(dif):.1e} (K ≤ {max(ea)})" if dif else '—') + ' |')
    low = max(y['low'] for x in A5['series'] for rs in x['frames'].values() for y in rs)
    out += ['', f"The defects below 2K, zero exactly, are at most {low:.1e} of their parts (the rounding of the {A5_PREC}-digit arithmetic where "
                "the rule is exact; it bounds that rounding, not the tail's)."]
    if A5.get('far'):
        rs = {}
        for y in A5['far']:
            rs.setdefault(y['r'], []).append(y)
        tq = {x['r']: x['t'] for x in A5['series']}
        parts, lims = [], []
        for r_, ys in sorted(rs.items(), reverse=True):
            Ks, rr = [y['K'] for y in ys], [y['ratio'] for y in ys]
            L = float(np.linalg.solve(np.array([[1.0, K ** -2.0, K ** -3.0] for K in Ks[-3:]]), rr[-3:])[0])
            d = np.diff(rr)
            dd = [d[i + 1] / d[i] for i in range(len(d) - 1)]
            q2 = [(Ks[i + 1] ** -2.0 - Ks[i + 2] ** -2.0) / (Ks[i] ** -2.0 - Ks[i + 1] ** -2.0) for i in range(len(Ks) - 2)]
            b = np.sqrt(1 + 2 * r_ ** 2)
            conj = (b * b + b + 2) / (b * (b - 1))
            ck = [abs(y['check'] / y['ratio'] - 1) for y in ys if 'check' in y]
            lims.append((r_, L, conj, 1 / tq[r_]))
            parts.append(f"τ/σ = {r_:g} ({', '.join(str(y['prec']) for y in ys) if len({y['prec'] for y in ys}) > 1 else ys[0]['prec']} digits): "
                         + ', '.join(f"{K}: {x:.6f}" for K, x in zip(Ks, rr)) + "; the successive differences' ratios "
                         + ', '.join(f"{x:.3f}" for x in dd) + " (K⁻²'s " + ', '.join(f"{x:.3f}" for x in q2) + f"); the limit by L + c/K² + d/K³ "
                         f"through the last three {L:.6f}, the formula's (b² + b + 2)/(b(b − 1)), b = (1 + 2τ²/σ²)^{{1/2}}, {conj:.6f} (A6), "
                         f"1/t {1 / tq[r_]:.4f}" + ((f"; the last K with {A5_FAR_CHECK} more digits the same in double" if max(ck) == 0 else
                                                     f"; the last K with {A5_FAR_CHECK} more digits {max(ck):.1e} apart") if ck else ''))
        tm = up_f(max(1 / c for _, _, c, _ in lims), 4)                 # the frame's largest t(v), rounded up
        st, rs = d61_step(float(tm))                                    # round 121: the first of D‴.61's tiers covering it
        out += ['', "The Mercer frame's rule further out (`a5_far`): its ratio per kernel by K — " + '; '.join(parts)
                + f" — the low defects at most {max(y['low'] for y in A5['far']):.1e}. Its limit lies "
                + ('below 1/t at every τ/σ' if all(L < it for _, L, _, it in lims) else 'NOT below 1/t at every τ/σ')
                + f", within {up_e(max(abs(L / c - 1) for _, L, c, _ in lims))} of the formula's — D‴.61 (iv)'s limit there (the frame's t(v) at "
                  f"most {tm}, below {st}'s r* = {f'{rs:.4f}' if rs < 1 else '1'}" + (" already, and with VA (iv) holds at every frame; " if st != 'VA' else "; ")
                + f"{D61_STATUS}; referee rounds 96, 101, 107)."]
    return out


A6_TAUS = (0.25, 0.5, 1.0, 2.0)                                         # A6 (referee round 84, D3ce-3): the scaled rules' exact errors
A6_VR = (0.5, 0.8, 1.0, 1.25, 2.0)                                     # A6: the variances v/v_hs beside the named scalings
A6_K = (20, 40, 60, 80, 120, 160)                                      # A6: the ratio per kernel e_{K−1}/e_K at these K
A6_MARGIN = 60                                                         # A6: digits beyond the 2K·log₁₀(1/t) the cancellation costs
A6_NEAR = 0.05                                                         # A6: |ln(v/v_hs)| below which a rule is summarised as near the cusp


def a6_scalings(r_):
    """A6's rules at τ/σ = r_ (σ = 1), their variances as Decimals in the current context: moment matching τ², the Fock scaling
    2τ²/(1 + σ*²), the Mercer-zero scaling τ²/b and the Mercer frame's rule 2τ²/(b + 1), b = (1 + 2τ²)^{1/2}, and v/v_hs in A6_VR,
    v_hs = τ²/σ*."""
    from decimal import Decimal
    tau2 = Decimal(repr(r_)) ** 2
    ss = (1 + tau2).sqrt()
    b = (1 + 2 * tau2).sqrt()
    vhs = tau2 / ss
    return ([('moment matching', tau2), ('Fock', 2 * tau2 / (2 + tau2)), ('Mercer-zero', tau2 / b), ('Mercer frame', 2 * tau2 / (b + 1))]
            + [(f"{x:g}·v_hs", Decimal(repr(x)) * vhs) for x in A6_VR])


D61_RSTAR = 0.8052                                                       # D‴.61 (iv)'s r*, rounded down (its round 93; pinned by a test)
D61_RSTAR_V = 0.9900                                                     # D‴.61's r* with its part V (the valley certificate to |η| ≤ 1.995)
D61_RSTAR_VA = 1                                                         # with its part VA (1.995 ≤ |η| < 2): (iv) at every frame, r < 1
D61_STATUS = ("by hand since referee round 148 (the valley condition on 1 < |η| < 2 along a three-segment path: "
              "review/RESEARCH_note_annulus.md, Theorem A); until then computer-assisted, D‴.61's two certificates — the branch and bound "
              "of its part D to |η| ≤ 1.8927 and its part VA from 1.89 (VA: referee rounds 120, 121; its adoption reviewed in round 122; D "
              "and VA alone since round 131, its part V an independent cross-check) — rigorous in ball arithmetic, the arithmetic of D‴.61's "
              "Appendix A.0 assumed correct, now independent cross-checks")
                                                                         # rounds 107, 121, 129 (D3dx-5, as D58 states it): pinned by a test
                                                                         # to D‴.61's raw
# D‴.61 (iv)'s tiers until referee round 148 (since then (iv) holds at every frame by hand; the tiers are printed as history): frames with r < D61_RSTAR by the branch and bound (|η| ≤ 1.8927), r < D61_RSTAR_V with V (to |η| ≤ 1.995), every
# r < 1 with VA (the open disc |η| < 2: c_need(r)² = 4r ln(1/r)/(1 − r) < 4 as ln x < x − 1) — so a scaled rule's rate is proved iff it
# has a Mehler frame (|ρ| < 1; its t(v) < 1 always); moment matching (ρ = 1) and |ρ| > 1 are D‴.62 (b)'s, by hand (referee round 132).


def d61_step(t):
    """The first of D‴.61's steps whose tier covers a frame of ratio t < 1, with that tier's r*: its branch and bound, V or VA (t < 1
    asserted, round 123: VA's tier is the open disc, and a frame's t(v) is below λ < 1)."""
    assert t < 1, t
    return ('its branch and bound', D61_RSTAR) if t < D61_RSTAR else ('V', D61_RSTAR_V) if t < D61_RSTAR_V else ('VA', D61_RSTAR_VA)


def up_e(v, d=1):
    """|v| in e-notation with d decimals, rounded up at the last shown digit: an upper bound printed outward (referee round 96, D3cq-4)."""
    import math
    v = abs(float(v))
    if v == 0:
        return f"{0:.{d}e}"
    e = math.floor(math.log10(v))
    m = math.ceil(v / 10.0 ** e * 10 ** d * (1 - 1e-12)) / 10 ** d
    if m >= 10:
        m, e = m / 10, e + 1
    return f"{m:.{d}f}e{e:+03d}"


def out_e(v, up, d=1):
    """A signed value in e-notation rounded outward at its last digit — toward +∞ (up) or −∞ — for a range's ends (referee round 101,
    D3cv-3)."""
    import math
    v = float(v)
    if v == 0:
        return f"{0:+.{d}e}"
    a = abs(v)
    e = math.floor(math.log10(a))
    m = a / 10.0 ** e * 10 ** d
    m = (math.ceil(m * (1 - 1e-12)) if (v > 0) == up else math.floor(m * (1 + 1e-12))) / 10 ** d
    if m >= 10:
        m, e = m / 10, e + 1
    return f"{'-' if v < 0 else '+'}{m:.{d}f}e{e:+03d}"


def dn_f(v, d=3):
    """v rounded down at d decimals: a threshold claimed to be stayed below (referee round 96)."""
    import math
    return f"{math.floor(v * 10 ** d) / 10 ** d:.{d}f}"


def up_f(v, d=3):
    """v rounded up at d decimals: a threshold claimed to be exceeded (referee round 96)."""
    import math
    return f"{math.ceil(v * 10 ** d) / 10 ** d:.{d}f}"


def rate_conj(r_, v):
    """The conjectured rate per kernel of `scaled_rule` with variance v (referee round 84, D3ce-3): max{h(τ²/v), h(vσ*²/(σ²τ²))},
    h(y) = (y − 1)/(y + 1), σ = 1, τ = r_ — the arguments are (σ*/σ)(v_hs/v) and (σ*/σ)(v/v_hs), v_hs = τ²σ/σ*, so the rate is
    tanh((ln(σ*/σ) + |ln(v/v_hs)|)/2), at least t with equality only at v_hs. Returns (the rate, the argument attaining it: 1 for τ²/v,
    2 for vσ*²/(σ²τ²))."""
    tau2, ss2 = r_ ** 2, 1 + r_ ** 2
    h = lambda y: (y - 1) / (y + 1)
    a, b = h(tau2 / v), h(v * ss2 / tau2)
    return (a, 1) if a >= b else (b, 2)


def a6_error(r_, v, K, prec):
    """`scaled_rule`'s squared relative L² error at σ = 1, τ = r_ (the rule of N(0, v), v a Decimal), exactly in prec-digit decimal
    arithmetic (referee round 84): E = 1 − 2√2(1 + σ*²)^{−1/2}Σ_ka_ke^{−y_k²/(2(1 + σ*²))} + σ*⁻¹Σ_{jk}a_ja_ke^{−(y_j − y_k)²/4} — the
    Gaussian inner products in closed form —, y_k = √v z_k, a_k = σ*w_k(v/τ²)^{1/2}e^{y_k²(1/v − 1/τ²)/2}, (z_k, w_k) `gh_nodes`."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        v = +v
        tau2 = Decimal(repr(r_)) ** 2
        ss2 = 1 + tau2
        ss = ss2.sqrt()
        zs, ws = gh_nodes(K, prec)
        sv = v.sqrt()
        ys = [sv * z for z in zs]
        g = (1 / v - 1 / tau2) / 2
        a = [ss * w * (v / tau2).sqrt() * (y * y * g).exp() for w, y in zip(ws, ys)]
        cross = sum((a_ * (-(y * y) / (2 * (ss2 + 1))).exp() for a_, y in zip(a, ys)), Decimal(0))
        gram = Decimal(0)
        for j in range(K):
            gram += a[j] * a[j]
            for k in range(j):
                gram += 2 * a[j] * a[k] * (-((ys[j] - ys[k]) ** 2) / 4).exp()
        return 1 - 2 * Decimal(2).sqrt() / (ss2 + 1).sqrt() * cross + gram / ss


def a6_job(args):
    """A6 (referee round 84, D3ce-3): at τ/σ = r_ the rule i of `a6_scalings` — its exact relative L² error (`a6_error`) at K − 1 and K for
    K in A6_K, in 2K·log₁₀(1/t) + A6_MARGIN digits, and the ratio per kernel e_{K−1}/e_K —, its variance and the conjecture's rate
    (`rate_conj`); a third argument replaces A6_K (A6X, referee round 101)."""
    from decimal import Decimal, localcontext
    import math
    r_, i = args[:2]
    Ks_ = args[2] if len(args) > 2 else A6_K
    t0 = time.time()
    ss = math.sqrt(1 + r_ ** 2)
    t = (ss - 1) / (ss + 1)
    rows = []
    for K in Ks_:
        prec = int(math.ceil(2 * K * math.log10(1 / t))) + A6_MARGIN
        with localcontext() as ctx:
            ctx.prec = prec
            name, v = a6_scalings(r_)[i]
            E0, E1 = a6_error(r_, v, K - 1, prec), a6_error(r_, v, K, prec)
            ok = E0 > 0 and E1 > 0
            rows.append(dict(K=K, prec=prec, e0=format(E0.sqrt(), '.25e') if E0 > 0 else None, e1=format(E1.sqrt(), '.25e') if E1 > 0 else None,
                             ratio=float((E0 / E1).sqrt()) if ok else None))
            tau2 = Decimal(repr(r_)) ** 2
            vv, vr = float(v), float(v / (tau2 / (1 + tau2).sqrt()))
    rate, arg = rate_conj(r_, vv)
    out = dict(r=r_, name=name, v=vv, v_hs=vr, rate=rate, arg=arg, t=t, rows=rows, secs=time.time() - t0)
    print(f"  A6 tau/sigma {r_} {name} (v/v_hs {vr:.6f}): " + ', '.join(f"{y['K']}: {y['ratio']}" for y in rows)
          + f" against 1/r {1 / rate:.6f} [{out['secs']:.0f}s]", flush=True)
    return out


def part_A6():
    """A6: `a6_job` for every τ/σ of A6_TAUS and rule of `a6_scalings`."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a6_job, [(r_, i) for r_ in A6_TAUS for i in range(4 + len(A6_VR))], chunksize=1)


A6X_RULES = ((0.25, 8), (0.5, 0), (1.0, 0), (2.0, 0), (2.0, 8), (1.0, 3), (2.0, 3))   # A6X (referee round 101, D3cv-1): (τ/σ,
A6X_K = (200, 240, 320)                                                 # `a6_scalings` index) — the fast branch's outlier (¼, 2·v_hs),
                                                                        # moment matching at ½, 1, 2, and the proved (2, 2·v_hs) and, since
                                                                        # round 123 (D3dr-2), the proved Mercer frame's rule at τ/σ = 1, 2:
                                                                        # the known cases, beyond K = 160


def part_A6X(old=None, idx=None):
    """A6X (referee round 101, D3cv-1): `a6_job` at A6X_K for the rules of A6X_RULES — with idx (indices into A6X_RULES; round 123) for
    those rules only, merged into old's records, every record in A6X_RULES' order."""
    from multiprocessing import Pool
    todo = [k for k in range(len(A6X_RULES)) if idx is None or k in idx]
    with Pool(N_WORKERS) as pool:
        new = pool.map(a6_job, [(*A6X_RULES[k], A6X_K) for k in todo], chunksize=1)
    if idx is None:
        return new
    have = {(x['r'], x['name']): x for x in (list(old) if old is not None else []) + list(new)}
    return [have[k] for k in [(r_, a6_scalings(r_)[i][0]) for r_, i in A6X_RULES] if k in have]


A6C_RULES = ((2.0, 3), (1.0, 3), (2.0, 8), (2.0, 0), (0.25, 8))       # A6C (referee round 129, D3dx-1, -4): (τ/σ, `a6_scalings` index) —
                                                                        # the proved Mercer frame's rule at τ/σ = 2, 1 (A6X's failing known
                                                                        # cases), the controls (2, 2·v_hs), proved, and (2, MM), frameless,
                                                                        # and the frameless outlier (¼, 2·v_hs)
A6C_K = (160, 400)                                                      # A6C: the ratio e_{K−1}/e_K at every K = 160 … 400
A6C_BITS = 450                                                          # A6C: Arb's bits beyond the 2K·log₂(1/t(v)) the cancellation costs
A6C_DIGITS = 70                                                         # A6C: the stored errors' significant digits
A6C_TERMS = tuple(range(3, 16))                                         # A6C: the orders n = 3 … 15 (four … sixteen terms)
A6C_FLOOR = (6, 10, 14, 18, 22, 26)                                     # A6C: the residual floor's numbers of powers N
A6C_COEF = ((160, 20), (160, 24), (160, 28), (200, 20), (200, 24), (200, 28))   # A6C: the coefficient fits (K₀, N powers) on [K₀, 400]
A6C_SPREAD = 1e-2                                                       # A6C: a coefficient printed while its fits agree within this
A6C_KMAX = 15                                                           # A6C: … and at most to c_15
A6C_TR_K0 = 200                                                         # A6C: the transient's rule on [A6C_TR_K0, 400], q^K carried


def a6c_v(r_, i):
    """`a6_scalings`' variance i at τ/σ = r_ (σ = 1) as an Arb ball in flint's current precision (referee round 129): τ², 2τ²/(2 + τ²),
    τ²/b, 2τ²/(b + 1) and x·τ²/σ* (x in A6_VR), b = (1 + 2τ²)^{1/2}, σ* = (1 + τ²)^{1/2}."""
    from flint import arb
    t2 = arb(repr(r_)) ** 2
    b = (1 + 2 * t2).sqrt()
    return [t2, 2 * t2 / (2 + t2), t2 / b, 2 * t2 / (b + 1)][i] if i < 4 else arb(repr(A6_VR[i - 4])) * t2 / (1 + t2).sqrt()


def a6c_nodes(K, bits):
    """The K-point Gauss–Hermite rule of N(0, 1) in Arb at bits (referee round 129, independent of `gh_nodes`): LAPACK's tridiagonal
    eigenvalues refined by Newton on He_K in Arb (the iterate's midpoint kept) to a step below 2^{20 − bits}, the middle node of an odd
    rule exactly 0, the rule symmetric; the weights K!/(K·He_{K−1}(z))²."""
    import math
    from scipy.linalg import eigh_tridiagonal
    from flint import arb, ctx
    ctx.prec = bits
    z0 = eigh_tridiagonal(np.zeros(K), np.sqrt(np.arange(1, K, dtype=float)), eigvals_only=True)
    fK, tol = arb(math.factorial(K)), arb(2) ** (20 - bits)

    def he(x):                                                          # (He_{K−1}(x), He_K(x)) by the three-term recurrence
        hm, h = arb(1), x
        for n in range(1, K):
            hm, h = h, x * h - n * hm
        return hm, h
    zs, ws = [], []
    for j, x0 in enumerate(z0[:(K + 1) // 2]):
        x = arb(float(x0))
        if K % 2 and j == (K - 1) // 2:
            x = arb(0)
        else:
            for it in range(60):
                hm, h = he(x)
                dx = h / (K * hm)
                x = arb((x - dx).mid())
                if it > 2 and abs(dx.mid()) < tol * (1 + abs(x.mid())):
                    break
            else:
                raise RuntimeError(f"a6c_nodes: Newton at K = {K}")
        zs.append(x)
        ws.append(fK / (K * he(x)[0]) ** 2)
    if K % 2:
        return zs + [-z for z in zs[-2::-1]], ws + ws[-2::-1]
    return zs + [-z for z in zs[::-1]], ws + ws[::-1]


def a6c_error(r_, i, K, bits):
    """`scaled_rule`'s relative L² error e_K at τ/σ = r_ (σ = 1), the rule i of `a6_scalings`, as an Arb ball at bits (referee round 129:
    python-flint, independent of `a6_error`'s decimal path and of `gh_nodes`) — e_K² = 1 − 2√2(1 + σ*²)^{−1/2}Σ_k a_ke^{−y_k²/(2(1 + σ*²))}
    + σ*⁻¹Σ_{jk}a_ja_ke^{−(y_j − y_k)²/4} (the Gaussian inner products), y_k = √v z_k, a_k = σ*w_k(v/τ²)^{1/2}e^{z_k²(1 − v/τ²)/2}, the
    nodes `a6c_nodes`. Leaves flint's precision at bits (`a6c_job` restores it)."""
    from flint import arb
    zs, ws = a6c_nodes(K, bits)
    t2 = arb(repr(r_)) ** 2
    v = a6c_v(r_, i)
    ss2 = 1 + t2
    ss, sv, g = ss2.sqrt(), v.sqrt(), (1 - v / t2) / 2
    ys = [sv * z for z in zs]
    a = [ss * w * (v / t2).sqrt() * (z * z * g).exp() for w, z in zip(ws, zs)]
    cross = sum((a_ * (-(y * y) / (2 * (ss2 + 1))).exp() for a_, y in zip(a, ys)), arb(0))
    gram = arb(0)
    for j in range(K):
        s = arb(0)
        for k in range(j):
            d = ys[j] - ys[k]
            s += a[k] * (-(d * d) / 4).exp()
        gram += a[j] * (a[j] + 2 * s)
    return (1 - 2 * arb(2).sqrt() / (ss2 + 1).sqrt() * cross + gram / ss).sqrt()


def a6c_str(x, digits=A6C_DIGITS):
    """An Arb ball's midpoint correctly rounded to digits significant digits (one rounding, the decimal context's), and its radius
    relative to it."""
    from decimal import Decimal, localcontext
    m, e = (int(u) for u in x.mid().man_exp())
    num, den = (m * 2 ** e, 1) if e >= 0 else (m, 2 ** -e)
    with localcontext() as c:
        c.prec = digits
        d = Decimal(num) / Decimal(den)
    rm, re_ = (int(u) for u in x.rad().mid().man_exp())
    return format(d, f'.{digits - 1}e'), float(Fraction(rm * 2 ** re_ if re_ >= 0 else Fraction(rm, 2 ** -re_)) / abs(Fraction(num, den)))


def a6c_job(args):
    """A6C (referee round 129, D3dx-1, -4): the rule (r_, i) of `a6_scalings` — e_K (`a6c_error`) at every K = A6C_K[0] − 1 … A6C_K[1]
    in ⌈2K·log₂(1/t(v))⌉ + A6C_BITS bits, stored to A6C_DIGITS significant digits (the balls' midpoints; their relative radii asserted
    below 10^{−A6C_DIGITS − 5}) —, its variance and the conjecture's rate (`rate_conj`); flint's precision restored."""
    import math
    import flint
    r_, i = args
    t0 = time.time()
    name, vd = a6_scalings(r_)[i]
    rate, arg = rate_conj(r_, float(vd))
    p0 = flint.ctx.prec
    Ks, es, bits, rad = list(range(A6C_K[0] - 1, A6C_K[1] + 1)), [], [], 0.0
    try:
        for K in Ks:
            bits.append(int(math.ceil(2 * K * math.log2(1 / rate))) + A6C_BITS)
            s, rr = a6c_str(a6c_error(r_, i, K, bits[-1]))
            es.append(s)
            rad = max(rad, rr)
    finally:
        flint.ctx.prec = p0
    assert rad < 10.0 ** (-A6C_DIGITS - 5), (r_, i, rad)
    out = dict(r=r_, name=name, i=i, v=float(vd), rate=rate, arg=arg, K0=Ks[0], K1=Ks[-1], e=es, bits=bits, rad=rad,
               flint=flint.__version__, secs=time.time() - t0)
    print(f"  A6C tau/sigma {r_} {name}: K = {Ks[0]} … {Ks[-1]}, the radii at most {rad:.1e} [{out['secs']:.0f}s]", flush=True)
    return out


def part_A6C():
    """A6C: `a6c_job` for the rules of A6C_RULES, in their order."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a6c_job, list(A6C_RULES), chunksize=1)


A6X_DIGITS = 26                                                        # A6X: the stored errors' significant digits ('.25e', `a6_job`)
A6X_CONS = 10                                                           # A6X (round 118): the last order within this many times its change
A6X_FLAT = 1e-3                                                         # A6X (round 118): B_K within this of B below the transient fit
A6X_TERMS = {3: 'four', 4: 'five', 5: 'six', 6: 'seven', 7: 'eight'}   # A6X: n + 1 terms, the order n of the successive fits


def a6x_resolvable(q, K=None, digits=A6X_DIGITS):
    """Referee rounds 118 (D3dm-4), 123 (D3dr-1): a transient Bq^K is resolvable in the ratios at K from 26-digit errors only if q^K ≫
    10^{−digits} — q^K > 10⁻²⁰: at A6's last K (160, the default) q ≳ 0.75, at the smallest K of the orders of `a6x_lines` (160, 120, 80,
    60, 40 for four … eight terms) q ≳ 0.75 … 0.32 —, not merely q > ½."""
    K = A6_K[-1] if K is None else K
    return q > 0 and q ** K > 10.0 ** (6 - digits)


def a6x_decreasing(o, b):
    """Rounds 118 (D3dm-1), 123 (D3dr-3): successive orders o decreasing in magnitude through the last pair, an order below its
    propagated input bound (b, the same index) counting as converged."""
    return all(abs(y) <= e_ or abs(y) < abs(x) for x, y, e_ in zip(o, o[1:], b[1:]))


def a6x_resolution(o, b):
    """Round 123 (D3dr-3): the orders' resolution — A6X_CONS times the last change, or the last order's input bound if that is larger;
    the limit is 1/t(v) to that resolution if the last order lies within it."""
    return max(A6X_CONS * abs(o[-1] - o[-2]), b[-1])


def a6x_transient_ok(tr, po):
    """Rounds 118 (D3dm-1, -2, -3), 123 (D3dr-4): the transient accounts for the residual — B > 0 (the two branches' contributions add),
    B_K within A6X_FLAT of B below the transient fit (no power of K beside q^K), the transient fit's deviation within its next order's
    change (or its input bound, if larger), and the power fits (po) worsening with more terms."""
    return (tr['B'] > 0 and tr['flat'] < A6X_FLAT and abs(tr['dev']) <= max(abs(tr['dev5'] - tr['dev']), tr['bnd'])
            and abs(po[0]) < abs(po[1]) < abs(po[2]))


def a6x_branches(r_, v):
    """The formula's two arguments at τ/σ = r_ (σ = 1) for the variance v (a Fraction or a Decimal, in its own arithmetic): h(τ²/v) and
    h(vσ*²/τ²), h(y) = (y − 1)/(y + 1). Returns (t(v), h_other, fast): t(v) the larger, h_other the other, fast whether t(v) is
    attained at vσ*²/τ² (the fast branch, `rate_conj`'s argument 2; a tie counted as the slow branch, as there)."""
    t2 = type(v)(repr(r_)) ** 2
    h = lambda y: (y - 1) / (y + 1)
    a, b = h(t2 / v), h(v * (1 + t2) / t2)
    return (a, b, False) if a >= b else (b, a, True)


def a6x_c2(r_, v):
    """Referee round 129 (D3dx-3): the K⁻² coefficient of the ratio per kernel off the cusp, e_{K−1}/e_K = (1/t(v))(1 + c₂/K² + O(K⁻³)),
    c₂ = s₁/2, s₁ = 1/16 − ϱ/(4(1 + ϱ)²) − 2ϱ(1 − ϱ)/(1 + ϱ)³ + √q(1 − ϱ)²/(8(1 − √q)(1 + ϱ)²), ϱ = t(v) on the fast branch and −t(v)
    on the slow, √q = h_other/t(v) signed (`a6x_branches`). The referee's derivation from D‴.61 (ii)'s moments — the cumulants' next
    coefficients d_j = −c_j·j(8j² + j + 1)/4, the 1/K term F₁(y) = y(1 + y)/(1 − y)³ + y/(8(1 − y)²) + y/(8(1 − y)) of log ĝ_K —
    by a two-dimensional Laplace step at the mass point: formal at next order (D‴.61 proves only the leading limit; the next-order
    uniformity and the Laplace remainder are not proved); applied formally to the frameless rules too. In v's own arithmetic (exact for
    a Fraction); None at the cusp (√q = 1 to 10⁻²⁰, where the ratio has a K⁻¹ term)."""
    t, ho, fast = a6x_branches(r_, v)
    one = type(v)(1)
    sq = ho / t
    if abs(one - sq) < one / 10 ** 20:
        return None
    g = t if fast else -t
    return (one / 16 - g / (4 * (1 + g) ** 2) - 2 * g * (1 - g) / (1 + g) ** 3 + sq * (1 - g) ** 2 / (8 * (1 - sq) * (1 + g) ** 2)) / 2


def a6x_carry(q, K, rem):
    """Round 129 (D3dx-2): the other branch's transient Bq^K enters an order (in place of its top power, on the order's own K — never
    below its smallest K) only if q^K at the order's smallest K is resolvable in the 26-digit errors (`a6x_resolvable`) and exceeds the
    order's own remainder there, rem: its next change — the power the carried order drops (K^{−n−1} for n + 1 terms) at that K, as the
    next carried order (Bq^K in it too, through one more K) fits it; None: no next order, not carried."""
    return a6x_resolvable(q, K) and rem is not None and q ** K > rem


def a6x_lines(A6X, A6, A6C=None):
    """A6X's line (referee rounds 101, 107, 113, 118, 123, 129): per rule the two fits' deviations from 1/t(v) — L + c/K² + d/K³ and
    L + c/K + d/K² — through A6's last three K and through A6X_K, the ratio's excess over 1/t(v) times K², and, from the stored 26-digit
    errors in exact rationals with t(v) exact (`rate_conj` on `a6_scalings`' Decimal variance), the limits by successive orders —
    L + Σ_{k=2..n+1} c_k/K^k through the last n + 1 K for n = 3 … 7 (four … eight terms, K from 160, 120, 80, 60, 40 to 320). The other
    branch's transient Bq^K, q = (h_other/t(v))² of the formula's other argument, replaces an order's top power on the order's own K
    where `a6x_carry` admits it (round 129, D3dx-2: q^K resolvable and above the order's own remainder at its smallest K, the dropped
    power there as the next carried order fits it); where it is resolvable but not carried, the reason is printed. Where it is resolvable at every order's smallest K (at A6's last K: the outlier) it
    is carried in every order and the rule is judged by the transient's own criterion (D3dx-4; round 118, 123): the five-term carried
    order L + c/K² + d/K³ + e/K⁴ + Bq^K through K = 120 … 320, its next order (with K⁻⁵, from K = 80), B_K = r_K/q^K at the two K of A6
    below that fit (at its own K B_K is B by construction), the residual's decay per kernel against the four-term power fit at A6's K, the
    power fits worsening; every carried order printed. The shrink factors between windows against a pure K⁻⁴ (first fit) and K⁻³ (second
    fit) remainder's — window geometry, the same for any rule. The verdict, rule by rule (round 123, D3dr-2): the orders decreasing in
    magnitude — an order below its propagated input bound (the 26-digit errors' rounding through the fit's weights) counting as converged
    (D3dr-3) — and the last within A6X_CONS times its last change (that, the resolution); else inconclusive — the limit, the premise, or
    the window, the expansion not yet asymptotic at the order's smallest K —, and where A6C's consecutive K show the orders decreasing on
    K ≥ A6C_K[0] (`a6c_analysis`), the window (round 129, D3dx-1). The next term (D3dx-3; `a6x_c2`, formal): c₂ per rule beside a fitted
    one (A6C's where computed, else the stored errors' with the limit exact, the eight-term order's seven powers through the last seven
    K), and the check with c₂ fixed (R_K − c₂/(t(v)K²), the same top powers through one K fewer). The proved rules are known cases of the
    check, and a proved rule's bias calibrates nothing (round 113, D3dh-1: the proved rules' biases differ, by the factor printed)."""
    from fractions import Fraction as Fr
    from decimal import Decimal, localcontext
    import math
    if A6X is None or len(A6X) == 0 or A6 is None:
        return []
    by = {(x['r'], x['name']): x for x in A6}
    byc = {(x['r'], x['name']): x for x in (list(A6C) if A6C is not None else [])}
    fit = lambda Ks, rr, p, q: float(np.linalg.solve(np.array([[1.0, K ** p, K ** q] for K in Ks], float), np.array(rr, float))[0])

    def solve(A, b):                                                    # exact Gauss–Jordan elimination
        A, b, n = [r_[:] for r_ in A], b[:], len(b)
        for c in range(n):
            piv = next(r_ for r_ in range(c, n) if A[r_][c] != 0)
            A[c], A[piv], b[c], b[piv] = A[piv], A[c], b[piv], b[c]
            for r_ in range(n):
                if r_ != c and A[r_][c] != 0:
                    f = A[r_][c] / A[c][c]
                    A[r_] = [a - f * cc for a, cc in zip(A[r_], A[c])]
                    b[r_] -= f * b[c]
        return [b[k] / A[k][k] for k in range(n)]
    pw = lambda e: (lambda K: Fr(1, int(K) ** e))
    one = lambda K: Fr(1)
    efit = lambda R, Ks, basis: solve([[f(K) for f in basis] for K in Ks], [R[K] for K in Ks])
    ratio_hp = lambda y: Fr(y['e0']) / Fr(y['e1']) if y.get('e0') and y.get('e1') else None   # e_{K−1}/e_K from the 26-digit errors

    def rel(s):                                                         # round 123 (D3dr-3): a stored error's rounding, relative
        d = Decimal(s)
        return Fr(1, 2 * 10 ** (A6X_DIGITS - 1)) / abs(Fr(d.scaleb(-d.adjusted())))

    def order(R, dR, Ks, basis):                                        # the fit and its input bound: L = w·R, w = A^{−T}e₀, |δL| ≤ Σ|w_k|δR_k
        A = [[f(K) for f in basis] for K in Ks]
        w = solve([list(c_) for c_ in zip(*A)], [Fr(int(j == 0)) for j in range(len(Ks))])
        return solve(A, [R[K] for K in Ks]), sum(abs(w_) * dR[K] for w_, K in zip(w, Ks))

    def geom(Ks0, Ks1, p, ps):                                          # a pure K^{−p} remainder: the error's shrink between windows
        e_ = lambda Ks: solve([[Fr(1)] + [Fr(1, int(K) ** q) for q in ps] for K in Ks], [Fr(1, int(K) ** p) for K in Ks])[0]
        return float(abs(e_(Ks0) / e_(Ks1)))

    def dec18(d):                                                       # a Decimal to 18 significant digits, trailing zeros dropped
        with localcontext() as c_:
            c_.prec = 18
            return format((+d).normalize(), 'f')
    rows_ = []
    for x in A6X:
        a = by.get((x['r'], x['name']))
        if a is None:
            continue
        Ka = [float(K) for K in A6_K]
        ra = [y['ratio'] for y in a['rows']]
        K1, r1 = [float(y['K']) for y in x['rows']], [y['ratio'] for y in x['rows']]
        if any(v is None for v in ra + r1):
            continue
        Kall, rall = Ka + K1, ra + r1
        win = lambda k0: (Kall[k0:k0 + 3], rall[k0:k0 + 3])
        i0, iL = len(Ka) - 3, len(Kall) - 3                            # 80–160 and 200–320
        d = {k0: (fit(*win(k0), -2.0, -3.0) * x['rate'] - 1, fit(*win(k0), -1.0, -2.0) * x['rate'] - 1) for k0 in (i0, iL - 1, iL)}
        c2d = None
        with localcontext() as ctx:                                     # t(v) exact from the Decimal variance (round 113)
            ctx.prec = 60
            vv = dict(a6_scalings(x['r'])).get(x['name'])
            tx = Fr(rate_conj(Decimal(repr(x['r'])), vv)[0]) if vv is not None else Fr(x['rate'])
            hs = [Fr(v_) for v_ in ((Decimal(repr(x['r'])) ** 2 / vv - 1) / (Decimal(repr(x['r'])) ** 2 / vv + 1),
                                    (vv * (1 + Decimal(repr(x['r'])) ** 2) / Decimal(repr(x['r'])) ** 2 - 1)
                                    / (vv * (1 + Decimal(repr(x['r'])) ** 2) / Decimal(repr(x['r'])) ** 2 + 1))] if vv is not None else []
            if vv is not None:
                c2d = a6x_c2(x['r'], vv)                                 # round 129 (D3dx-3): the next term, formal
        ys = list(a['rows']) + list(x['rows'])
        R = {int(y['K']): ratio_hp(y) for y in ys}
        cx = byc.get((x['r'], x['name']))
        z = dict(x=x, d=d, orders=None, tr=None, q=0.0, c2=None, a6c=a6c_analysis(cx) if cx is not None else None)
        if all(v is not None for v in R.values()):
            Ks = sorted(R)
            dR = {int(y['K']): R[int(y['K'])] * (rel(y['e0']) + rel(y['e1'])) for y in ys}
            ho = (hs[0] if abs(hs[1] - tx) < abs(hs[0] - tx) else hs[1]) if len(hs) == 2 else None   # the other argument's h
            q = float((ho / tx) ** 2) if ho is not None else 0.0
            qf = Fr(q)
            qK = lambda K: qf ** int(K)
            pcb = []                                                    # the power orders, and the next one beyond the last (its top term)
            for n in list(A6X_TERMS) + [max(A6X_TERMS) + 1]:
                if n + 1 <= len(Ks):
                    pcb.append(order(R, dR, Ks[-(n + 1):], [one] + [pw(e) for e in range(2, n + 2)]))
            po = [float(c[0] * tx - 1) for c, _ in pcb[:len(A6X_TERMS)]]
            pb = [float(b_ * tx) for _, b_ in pcb[:len(A6X_TERMS)]]
            trr = a6x_resolvable(q, Ks[-(min(A6X_TERMS) + 1)])          # resolvable at every order's smallest K: the transient's rule
            cc = {}                                                     # the carried orders: L + Σ_{k=2..n} c_k/K^k + Bq^K, own n + 1 K
            for n in list(A6X_TERMS) + [max(A6X_TERMS) + 1]:
                if n + 1 <= len(Ks) and (trr or a6x_resolvable(q, Ks[-n]) or a6x_resolvable(q, Ks[-(n + 1)])):
                    cc[n] = order(R, dR, Ks[-(n + 1):], [one] + [pw(e) for e in range(2, n + 1)] + [qK])
            ow, ob, cw, rz, rem, qk, k0s = [], [], [], [], [], [], []
            for j, n in enumerate(A6X_TERMS):                           # round 129 (D3dx-2): carried in place of the top power, own K
                k0s.append(Ks[-(n + 1)])                                # the remainder: the power it drops, K^{−n−1}, at its smallest K,
                rem.append(float(abs(cc[n + 1][0][-2] * tx) / k0s[-1] ** (n + 1)) if n + 1 in cc else None)   # by the next carried order
                qk.append(q ** k0s[-1])
                rz.append(trr or a6x_carry(q, k0s[-1], rem[-1]))
                c, b_ = cc[n] if rz[-1] else pcb[j]
                cw.append(c)
                ow.append(float(c[0] * tx - 1))
                ob.append(float(b_ * tx))
            z.update(orders=po, pbnd=pb, ow=ow, bnd=ob, rz=rz, rem=rem, qk=qk, k0s=k0s, q=q,
                     Bw=[float(c[-1]) if r_ else None for c, r_ in zip(cw, rz)])
            if trr:                                                     # the transient's own analysis (its criterion; round 118, 123)
                c4 = efit(R, Ks[-4:], [one, pw(2), pw(3), pw(4)])
                ct = cw[1]                                              # the five-term carried order L + c/K² + d/K³ + e/K⁴ + Bq^K
                res = [(K, R[K] - sum(c * f(K) for c, f in zip(c4, [one, pw(2), pw(3), pw(4)]))) for K in A6_K[:-1]]
                dec = [float(r2 / r1_) ** (1 / (k2 - k1)) for (k1, r1_), (k2, r2) in zip(res, res[1:]) if r1_ * r2 > 0]
                pp = lambda K: ct[0] + ct[1] / K ** 2 + ct[2] / K ** 3 + ct[3] / K ** 4
                kb = [K for K in Ks if K < Ks[-5]][-2:]                 # round 123 (D3dr-4): B_K below the fit (at its K it is B)
                Bk = [float((R[K] - pp(Fr(K))) / qf ** int(K)) for K in kb]
                tr = dict(q=q, B=float(ct[4]), dev=ow[1], dev5=ow[2], bnd=ob[1], dec=dec, K0=Ks[-5], K05=Ks[-6], kb=kb,
                          flat=max(abs(b_ / float(ct[4]) - 1) for b_ in Bk), mono=a6x_decreasing(ow, ob))
                tr['ok'] = a6x_transient_ok(tr, po)
                z['tr'] = tr
            else:                                                       # round 123 (D3dr-2, -3): the orders' verdict, rule by rule
                z.update(mono=a6x_decreasing(ow, ob), pmono=a6x_decreasing(po, pb), res=a6x_resolution(ow, ob))
                z['cons'] = abs(ow[-1]) <= z['res']
            if c2d is not None:                                         # round 129 (D3dx-3): c₂ fitted and the check with c₂ fixed
                c2 = Fr(c2d)
                m = max(A6X_TERMS)
                eps = {K: R[K] * tx - 1 for K in Ks}
                cf = efit(eps, Ks[-m:], [pw(e) for e in range(2, m + 2)])[0]
                z.update(c2=c2d, c2fit=float(cf), c2rel=float(abs(cf / c2 - 1)), c2src='the stored errors')
                if z['a6c'] is not None and z['a6c'].get('c2rel') is not None:
                    z.update(c2fit=z['a6c']['c2fit'], c2rel=z['a6c']['c2rel'], c2src='A6C')
                Rm = {K: R[K] - c2 / (tx * K * K) for K in Ks}
                o2, b2 = [], []
                for n in A6X_TERMS:
                    c, b_ = order(Rm, dR, Ks[-n:], [one] + [pw(e) for e in range(3, n + 2)])
                    o2.append(float(c[0] * tx - 1))
                    b2.append(float(b_ * tx))
                z.update(o2=o2, b2=b2, mono2=a6x_decreasing(o2, b2), res2=a6x_resolution(o2, b2))
                z['cons2'] = abs(o2[-1]) <= z['res2']
        frame = abs(1 + 2 / x['r'] ** 2 - 2 / float(x['v'])) < 1 - 1e-12
        fr = frame                                                      # proved iff a Mehler frame: D‴.61 (iv) at every frame with VA
        ek = [(ra[-1] * x['rate'] - 1) * Ka[-1] ** 2] + [(v * x['rate'] - 1) * K ** 2 for K, v in zip(K1, r1)]
        z.update(fr=fr, frame=frame, ek=ek, i0=i0, iL=iL, prec=max(y['prec'] for y in x['rows']))
        rows_.append(z)
    if not rows_:
        return []
    shr = lambda z, a_, b_, j: abs(z['d'][a_][j]) / max(abs(z['d'][b_][j]), 1e-300)
    pf, nf = [z for z in rows_ if z['fr']], [z for z in rows_ if not z['fr']]
    sig3 = lambda v, up: (lambda e: (math.ceil if up else math.floor)(v / 10 ** e) * 10 ** e)(math.floor(math.log10(abs(v))) - 2) if v else 0.0
    rng3 = lambda vs: (f"{sig3(min(vs), False):.3g}" if len(vs) == 1 or min(vs) == max(vs) else
                       f"{sig3(min(vs), False):.3g}–{sig3(max(vs), True):.3g}")
    nrule = lambda n: f"{n} rule" + ('' if n == 1 else 's')
    rl = lambda zs: ', '.join(f"({z['x']['r']:g}, {z['x']['name']})" for z in zs)
    terms = lambda ns: (', '.join(A6X_TERMS[n] for n in ns[:-1]) + ' and ' if len(ns) > 1 else '') + f"{A6X_TERMS[ns[-1]]}-term"
    Kw = lambda k0: [int(K) for K in (list(A6_K) + list(A6X_K))[k0:k0 + 3]]
    conv = lambda z: z.get('a6c') is not None and z['a6c'].get('mono') and z['a6c'].get('cons')   # A6C: the orders converge on K ≥ 160
    jn = lambda z, key, ns, f: ', '.join(f.format(z[key][n - min(A6X_TERMS)]) for n in ns)
    z0 = rows_[0]
    parts = []
    for z in rows_:
        o_, t_ = z['orders'], z['tr']
        hy = [n for n, r_ in zip(A6X_TERMS, z.get('rz') or []) if r_]
        nc = [n for n, r_, k0 in zip(A6X_TERMS, z.get('rz') or [], z.get('k0s') or []) if not r_ and a6x_resolvable(z['q'], k0)]
        parts.append(f"τ/σ = {z['x']['r']:g}, {z['x']['name']} ({'a frame' if z['frame'] else 'no frame'}"
                     f"; {z['prec']} digits): "
                     f"{z['d'][z['i0']][0]:+.2e} → {z['d'][z['iL']][0]:+.2e} and {z['d'][z['i0']][1]:+.2e} → {z['d'][z['iL']][1]:+.2e}, the excess "
                     f"times K² at K = {A6_K[-1]}, " + ', '.join(str(int(y['K'])) for y in z['x']['rows']) + ' ' + ', '.join(f"{v:+.4f}" for v in z['ek'])
                     + ('' if o_ is None else "; by successive orders (four … eight terms) " + ', '.join(f"{v:+.1e}" for v in o_))
                     + ('' if t_ is not None or not hy else
                        f"; with the other branch's transient Bq^K, q = (h_other/t(v))² = {z['q']:.6f} per kernel, carried in the {terms(hy)} orders "
                        f"(at their smallest K, {jn(z, 'k0s', hy, '{}')}, q^K {jn(z, 'qk', hy, '{:.1e}')} above the order's own remainder there, "
                        f"its next change — the power it drops, by the next order with Bq^K — {jn(z, 'rem', hy, '{:.1e}')}; Bq^K in place of the "
                        f"top power, on the order's own K): "
                        + jn(z, 'ow', hy, '{:+.1e}') + f" (B {jn(z, 'Bw', hy, '{:.3g}')})")
                     + ('' if t_ is not None or not nc else
                        f"; the other branch's transient Bq^K, q = (h_other/t(v))² = {z['q']:.6f} per kernel, not carried in the {terms(nc)} orders, "
                        f"where q^K > 10⁻²⁰ at their smallest K ({jn(z, 'k0s', nc, '{}')}): q^K {jn(z, 'qk', nc, '{:.1e}')} there is below the order's "
                        f"own remainder, its next change (the power the order would drop for Bq^K, at that K, by the next order with Bq^K) " + ', '.join(
                            'none' if z['rem'][n - min(A6X_TERMS)] is None else f"{z['rem'][n - min(A6X_TERMS)]:.1e}" for n in nc)
                        + (f"; on K ≥ {A6C_K[0]}, where the orders converge (A6C), q^{A6C_K[0]} = {z['a6c']['qK0']:.1e}" if conv(z) else ''))
                     + ('' if t_ is None else f"; with the other branch's transient Bq^K, q = (h_other/t(v))² = {t_['q']:.6f} per kernel, "
                        f"through K = {t_['K0']} … {A6X_K[-1]}: {t_['dev']:+.1e} (B {t_['B']:.4f}; B_K = r_K/q^K at K = {', '.join(str(K) for K in t_['kb'])}, "
                        f"outside the fit, within {t_['flat']:.1e} of B — at the fit's own K it is B by construction —: no power of K beside q^K to that "
                        f"resolution), its next "
                        f"order (with K⁻⁵, from K = {t_['K05']}) {t_['dev5']:+.1e}, the residual against the four-term fit "
                        f"falling by " + ', '.join(f"{v:.5f}" for v in t_['dec']) + f" per kernel over K = {A6_K[0]} … {A6_K[-2]}; its carried orders "
                        f"(four … eight terms, Bq^K in place of the top power, on each order's K) " + ', '.join(f"{v:+.1e}" for v in z['ow'])))
    out = ['', f"A6 beyond K = {A6_K[-1]} on the fast branch (A6X, referee rounds 101, 107, 113, 118, 123, 129): the two fits' deviations from 1/t(v) "
               f"through K = {', '.join(str(K) for K in A6_K[-3:])} → through {', '.join(str(K) for K in A6X_K)}, and the limits from the "
               f"{A6X_DIGITS}-digit errors (t(v) exact) — " + '; '.join(parts) + " —"]
    g = [geom(Kw(z0['i0']), Kw(z0['iL']), 4, (2, 3)), geom(Kw(z0['i0']), Kw(z0['iL']), 3, (1, 2)),
         geom(Kw(z0['iL'] - 1), Kw(z0['iL']), 4, (2, 3)), geom(Kw(z0['iL'] - 1), Kw(z0['iL']), 3, (1, 2))]
    nlab = "without a frame"                                            # round 121: every unproved rule is one without a frame
    for lab, zs in ((nlab, nf), ("proved", pf)):
        if zs:
            out[-1] += (f" at the {nrule(len(zs))} {lab}" if lab == "without a frame" else f" at the {'' if len(zs) == 1 else f'{len(zs)} '}{lab} "
                        f"{'rule' if len(zs) == 1 else 'rules'}") + (f" the deviations shrink by {rng3([shr(z, z['i0'], z['iL'], 0) for z in zs])} and "
                        f"{rng3([shr(z, z['i0'], z['iL'], 1) for z in zs])} times, in the later window "
                        f"({', '.join(str(K) for K in Kw(z0['iL'] - 1))} → {', '.join(str(K) for K in Kw(z0['iL']))}) by "
                        f"{rng3([shr(z, z['iL'] - 1, z['iL'], 0) for z in zs])} and {rng3([shr(z, z['iL'] - 1, z['iL'], 1) for z in zs])};")
    out[-1] += (f" a pure K⁻⁴ and K⁻³ remainder gives {g[0]:.4g} and {g[1]:.4g}, {g[2]:.4g} and {g[3]:.4g} (window geometry, any rule);"
                + (" every deviation shrinks as the window moves out (a necessary condition)" if all(
                    shr(z, z['i0'], z['iL'], j) > 1 for z in rows_ for j in (0, 1)) else " NOT every deviation shrinks as the window moves out"))
    nofit = [z for z in rows_ if z['tr'] is None and z['orders'] is not None]
    trs = [z for z in rows_ if z['tr'] is not None]
    nfast = sum(1 for y in A6 if y['arg'] == 2 and not abs(1 + 2 / y['r'] ** 2 - 2 / float(y['v'])) < 1 - 1e-12)   # frameless fast-branch rules
    if nofit:                                                           # round 123 (D3dr-2): each rule's own flags
        ok_ = [z for z in nofit if z['mono'] and z['cons']]
        mo_ = [z for z in nofit if z['mono'] and not z['cons']]
        no_ = [z for z in nofit if not z['mono']]
        oth = [z for z in ok_ if z.get('a6c') is not None and z['a6c'].get('g4') is not None]   # A6C's rules where the check passes
        k8 = (list(A6_K) + list(A6X_K))[-(max(A6X_TERMS) + 1)]           # the eight-term order's smallest K

        def why(zs):                                                    # round 129 (D3dx-1): the window, where A6C shows the orders converge
            w_ = [z for z in zs if conv(z) and z['a6c'].get('g4') is not None]
            return (" — inconclusive: the limit, the premise, or the window (the expansion not yet asymptotic at the order's smallest K)"
                    + ('' if not w_ else f"; at {rl(w_)} the window: on K ≥ {A6C_K[0]} the orders decrease (A6C), to "
                       + ', '.join(f"{z['a6c']['orders'][-1]:+.1e}" for z in w_) + f" at {a6c_nw(max(A6C_TERMS) + 1)} terms, the last within ten "
                       "times its last change, the excess's coefficients growing " + ', '.join(
                           f"{z['a6c']['g4'][0]:.1f}–{z['a6c']['g4'][1]:.1f}" for z in w_) + " times per order (the last four ratios c_k/c_{k−1}"
                       + (f"; at most {max(z['a6c']['g4'][1] for z in oth):.1f} at A6C's rules where the check passes" if oth else '')
                       + f"), so that at the eight-term order's smallest K, {k8}, the terms c_k/K^k shrink by only " + ', '.join(
                           f"{dn_f(z['a6c']['g4'][0] / k8, 2)}–{up_f(z['a6c']['g4'][1] / k8, 2)}" for z in w_) + " per order: the expansion is not yet "
                       "asymptotic there"))
        out[-1] += (f"; by the orders (four … eight terms, with Bq^K wherever it is carried; an order below its "
                    f"propagated input bound — the {A6X_DIGITS}-digit errors' rounding through the fit, here at most "
                    f"{up_e(max(max(z['bnd']) for z in nofit))} — counting as converged), rule by rule at the {nrule(len(nofit))} without a transient "
                    f"resolvable at K = {A6_K[-1]}, (τ/σ, rule) = " + '; '.join(
                        ([f"{rl(ok_)}: decreasing in magnitude through eight terms, the last within ten times its last change (that, the resolution: "
                          + ', '.join(up_e(z['res']) for z in ok_) + ") — 1/t(v) to that resolution"] if ok_ else [])
                        + ([f"{rl(mo_)}: decreasing in magnitude through eight terms, NOT within ten times the last change" + why(mo_)] if mo_ else [])
                        + ([f"{rl(no_)}: NOT decreasing in magnitude through eight terms" + why(no_)] if no_ else [])))
        kp = [z for z in nofit if z['fr']]
        if kp:
            np_ = sum(1 for z in kp if z['mono'] and z['cons'])
            out[-1] += (f"; the known cases, the {nrule(len(kp))} proved: the check passes at {np_} of {len(kp)}"
                        + (" — a failure is inconclusive, it does not show the limit off" if np_ < len(kp) else ''))
        neg = [z for z in nofit if any(z['rz']) and not z['mono']]       # rounds 123, 129 (D3dr-1, D3dx-2): the other branch's transient
        pos = [z for z in nofit if any(z['rz']) and z['mono'] and not z['pmono']]
        nrm = [z for z in nofit if not any(z['rz']) and not (z['mono'] and z['cons']) and any(a6x_resolvable(z['q'], k0) for k0 in z['k0s'])]
        if neg:
            out[-1] += (f"; at {rl(neg)}, with the other branch's transient carried where it exceeds the order's own remainder, the orders still "
                        "do NOT decrease: the transient does not account for the stall (a negative result)")
        if pos:
            out[-1] += (f"; at {rl(pos)} the orders decrease only with the other branch's transient carried: the power fits' stall is that "
                        "unmodelled transient (a positive result)")
        if nrm:
            cv = [z for z in nrm if conv(z)]
            out[-1] += (f"; at {rl(nrm)} the other branch's transient, resolvable in the {A6X_DIGITS}-digit errors at some orders' smallest K, "
                        "lies below the order's own remainder at each: not carried, it does not account for the stall"
                        + ('' if not cv else f", and where the orders converge, on K ≥ {A6C_K[0]} (A6C), q^{A6C_K[0]} is at most "
                           f"{up_e(max(z['a6c']['qK0'] for z in cv))}") + " (a negative result)")
    c2z = [z for z in rows_ if z.get('c2') is not None]
    if c2z:                                                             # round 129 (D3dx-3): the next term, formal
        out[-1] += ("; the next term (referee round 129, derived from D‴.61 (ii)'s moments by a two-dimensional Laplace step at the mass point; "
                    "formal at next order — D‴.61 proves only the leading limit, the next-order uniformity and the Laplace remainder are not "
                    "proved): e_{K−1}/e_K = (1/t(v))(1 + c₂/K² + O(K⁻³)), c₂ = s₁/2, s₁ = 1/16 − ϱ/(4(1 + ϱ)²) − 2ϱ(1 − ϱ)/(1 + ϱ)³ + "
                    "√q(1 − ϱ)²/(8(1 − √q)(1 + ϱ)²), ϱ = t(v) on the fast branch (−t(v) on the slow), √q = h_other/t(v) signed; per rule c₂ "
                    "and a fitted one (the limit exact; A6C's where computed, else from the stored errors, the eight-term order's seven powers "
                    "through the last seven K), relative: " + '; '.join(
                        f"({z['x']['r']:g}, {z['x']['name']}) {dec18(z['c2'])} and {z['c2fit']:+.13e} ({z['c2src']}), {z['c2rel']:.1e}"
                        for z in c2z))
        cz = [z for z in c2z if z['tr'] is None]
        if cz:
            kp2 = [z for z in cz if z['fr']]
            out[-1] += ("; with c₂ fixed (R_K − c₂/(t(v)K²) = L + Σ_{k=3..n+1} c_k/K^k through the last n K: the same top powers, one K fewer) "
                        "the orders " + '; '.join(
                            f"({z['x']['r']:g}, {z['x']['name']}) " + ', '.join(f"{v:+.2e}" for v in z['o2'])
                            + (f", decreasing, the last within ten times its last change ({up_e(z['res2'])})" if z['mono2'] and z['cons2'] else
                               ", decreasing, NOT within ten times its last change" if z['mono2'] else ", NOT decreasing") for z in cz)
                        + ('' if not kp2 else f" — the known cases, the {nrule(len(kp2))} proved: with c₂ fixed the check passes at "
                           f"{sum(1 for z in kp2 if z['mono2'] and z['cons2'])} of {len(kp2)}"))
    if trs:
        tok, tno = [z for z in trs if z['tr']['ok']], [z for z in trs if not z['tr']['ok']]
        if tok:
            out[-1] += (f"; at the transient's {nrule(len(tok))} ({rl(tok)}: resolvable at every order's smallest K) the power fits worsen with "
                        f"more terms, B > 0, B_K flat (within {A6X_FLAT:g} of B outside the fit), and the transient fit is within "
                        + ', '.join(f"{up_e(z['tr']['dev'])} (rounded up)" for z in tok) + " of 1/t(v), inside its next order's change (or its "
                        "input bound, if larger) " + ', '.join(up_e(max(abs(z['tr']['dev5'] - z['tr']['dev']), z['tr']['bnd'])) for z in tok)
                        + " (that, the resolution; the transient's own q^{2K} term unmodelled) — the transient's own criterion, not the orders' "
                        "decrease: its carried orders " + ', '.join("decrease in magnitude" if z['tr']['mono'] else "do NOT decrease in magnitude"
                                                                   for z in tok))
        if tno:
            out[-1] += f"; at {rl(tno)} the transient fit does NOT account for the residual"
    pb_ = {}
    for z in pf:
        pb_.setdefault(z['x']['r'], []).append(z['d'][z['iL']][0])
    rat = max(((max(abs(v) for v in vs) / min(abs(v) for v in vs), r_, vs) for r_, vs in pb_.items()
               if len(vs) > 1 and min(abs(v) for v in vs) > 0), default=None, key=lambda w: w[0])
    out[-1] += (" — numerical, consistent with an expansion in integer powers of 1/K from K⁻² (no K⁻¹ term) and the transient (a non-integer "
                "power can pass this check), not a proof; a proved rule's bias calibrates nothing here, since the proved rules' biases differ"
                + ('' if rat is None else f" (by {rat[0]:.1f} times at τ/σ = {rat[1]:g}: " + ' and '.join(f"{v:+.2e}" for v in rat[2])
                   + f", the first fit through K = {', '.join(str(K) for K in A6X_K)})")
                + " and at τ/σ = 0.25 no far "
                f"fast-branch rule has a frame (A6X extends {len(nf)} of the {nfast} fast-branch rules without a frame).")
    return out


def a6c_nw(n):
    """A number of terms in words (A6C's orders: four … sixteen terms)."""
    return {4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 16: 'sixteen'}.get(n, str(n))


def a6c_win(K0, K1, N):
    """N integers spread uniformly in 1/K over [K0, K1] (the referee's windows, round 129), asserted distinct."""
    Ks = sorted({int(round(1 / u)) for u in np.linspace(1 / K0, 1 / K1, N)})
    assert len(Ks) == N, (K0, K1, N)
    return Ks


_A6C_MEMO = {}


def a6c_analysis(x):
    """A6C's analysis of one record (referee round 129, D3dx-1, -3, -4; memoised on its errors), from its A6C_DIGITS-digit errors in
    exact rationals (python-flint's fmpq), t(v), q and c₂ (`a6x_branches`, `a6x_c2`) from `a6_scalings`' variance at 100 digits — the
    ratios R_K = e_{K−1}/e_K at every K of A6C_K. A rule without a transient resolvable at K₀ = A6C_K[0] (`a6x_resolvable`): the orders by
    windows on [K₀, A6C_K[1]] (n + 1 K spread uniformly in 1/K, `a6c_win`; limit free; L + Σ_{k=2..n+1} c_k/K^k for n in A6C_TERMS) with
    their input bounds (the stored digits' half unit and the balls' radius, through the fit's weights), their decrease and resolution
    (`a6x_decreasing`, `a6x_resolution`); the residual floor — N powers K⁻² … K^{−N−1} of the excess ε_K = t(v)R_K − 1, the limit exact,
    through N such K: the largest |ε_K − fit| over every K (N in A6C_FLOOR) —; a K⁻¹ term beside 15 powers (its coefficient times t(v));
    the coefficients c_k of ε_K (the limit exact) — the median of the fits A6C_COEF, kept while their relative spread is at most
    A6C_SPREAD (k ≤ A6C_KMAX) —, the last four ratios' range in magnitude (g4) and c₂ fitted (the median). The transient's rule (the
    outlier): on [A6C_TR_K0, A6C_K[1]] the power orders and those with q^K(B₀ + … + B_J/K^J), J = 1, 2, in the basis (B₀ per order, B₁/B₀
    at the last), q^{2K} at A6C_TR_K0 (the transient's own unmodelled next term) and c₂ fitted (the limit exact, 20 powers and
    q^K(B₀ + B₁/K) on that window)."""
    import statistics
    from decimal import Decimal, localcontext
    from flint import fmpq, fmpq_mat
    key = (x['r'], x['name'], x['K0'], tuple(x['e']), x.get('rad'))
    if key in _A6C_MEMO:
        return _A6C_MEMO[key]
    fq = lambda f: fmpq(f.numerator, f.denominator)
    fl = lambda v: float(Fraction(int(v.p), int(v.q)))
    K0, K1 = A6C_K
    Ke = list(range(x['K0'], x['K1'] + 1))
    with localcontext() as c_:
        c_.prec = 100
        vv = dict(a6_scalings(x['r']))[x['name']]
        t_, ho, fast = a6x_branches(x['r'], vv)
        qd = (ho / t_) ** 2
        c2 = a6x_c2(x['r'], vv)
        qp = {K: fq(Fraction(qd ** K)) for K in Ke}
        t = fq(Fraction(t_))
    E_ = {K: Fraction(s) for K, s in zip(Ke, x['e'])}
    dE = {K: Fraction(1, 2 * 10 ** (A6C_DIGITS - 1)) / abs(Fraction(Decimal(s).scaleb(-Decimal(s).adjusted()))) + Fraction(x.get('rad', 0.0))
          for K, s in zip(Ke, x['e'])}
    R = {K: fq(E_[K - 1] / E_[K]) for K in Ke if K - 1 in E_ and K0 <= K <= K1}
    dR = {K: fl(R[K]) * float(dE[K - 1] + dE[K]) for K in R}
    eps = {K: R[K] * t - 1 for K in R}
    P = lambda K, k: fmpq(1, K ** k)

    def solve(rows, rhs):
        X = fmpq_mat(len(rows), len(rows[0]), [v for r_ in rows for v in r_]).solve(fmpq_mat(len(rhs), 1, rhs))
        return [X[i, 0] for i in range(len(rows[0]))]

    def order(Ks, row, rhs):                                            # the fit and its input bound (the weights w = A^{−T}e₀)
        A = [row(K) for K in Ks]
        w = solve([list(c_) for c_ in zip(*A)], [fmpq(int(j == 0)) for j in range(len(Ks))])
        return solve(A, [rhs[K] for K in Ks]), sum(abs(fl(w_)) * dR[K] for w_, K in zip(w, Ks))
    z = dict(q=float(qd), qK0=float(qd ** K0), qK0s=format(qd ** K0, '.1e') if qd else '0', c2=c2, tr=None, c2fit=None, c2rel=None, g4=None)
    pr = lambda n: (lambda K: [fmpq(1)] + [P(K, k) for k in range(2, n + 2)])
    if not a6x_resolvable(float(qd), K0):
        ob = [order(a6c_win(K0, K1, n + 1), pr(n), R) for n in A6C_TERMS]
        o, b = [fl(c[0] * t - 1) for c, _ in ob], [bb * fl(t) for _, bb in ob]
        z.update(orders=o, bnd=b, mono=a6x_decreasing(o, b), res=a6x_resolution(o, b))
        z['cons'] = abs(o[-1]) <= z['res']
        fls = []
        with localcontext() as c_:
            c_.prec = 120
            dq = lambda v: Decimal(int(v.p)) / Decimal(int(v.q))
            ed = {K: dq(eps[K]) for K in R}
            for N in A6C_FLOOR:
                Ks = a6c_win(K0, K1, N)
                cd = [dq(v) for v in solve([[P(K, k) for k in range(2, N + 2)] for K in Ks], [eps[K] for K in Ks])][::-1]
                mx = Decimal(0)
                for K in R:
                    u, p = Decimal(1) / K, Decimal(0)
                    for cc in cd:
                        p = p * u + cc
                    mx = max(mx, abs(ed[K] - p * u * u))
                fls.append(float(mx))
        Ks = a6c_win(K0, K1, 16)
        z.update(floor=fls, k1=fl(solve([[fmpq(1)] + [P(K, k) for k in range(1, 16)] for K in Ks], [R[K] for K in Ks])[1] * t))
        fq_ = [solve([[P(K, k) for k in range(2, N + 2)] for K in a6c_win(a_, K1, N)], [eps[K] for K in a6c_win(a_, K1, N)])
               for a_, N in A6C_COEF]
        fits = [[fl(v) for v in f_] for f_ in fq_]
        c2s = sorted(Fraction(int(f_[0].p), int(f_[0].q)) for f_ in fq_)
        c2q = (c2s[(len(c2s) - 1) // 2] + c2s[len(c2s) // 2]) / 2       # c₂'s median in exact rationals
        coef = []
        for k in range(2, A6C_KMAX + 1):
            vs = [f_[k - 2] for f_ in fits]
            med = statistics.median(vs)
            sp = (max(vs) - min(vs)) / abs(med)
            if sp > A6C_SPREAD:
                break
            coef.append((k, med, sp))
        rat = [abs(b_[1] / a_[1]) for a_, b_ in zip(coef, coef[1:])]
        z.update(coef=coef, c2fit=float(c2q), c2rel=float(abs(c2q / Fraction(c2) - 1)) if c2 is not None else None,
                 g4=(min(rat[-4:]), max(rat[-4:])) if len(rat) >= 4 else None)
    else:                                                               # the transient's rule: q^K(B₀ + … + B_J/K^J) carried
        k0 = A6C_TR_K0
        po = [fl(solve([pr(n)(K) for K in a6c_win(k0, K1, n + 1)], [R[K] for K in a6c_win(k0, K1, n + 1)])[0] * t - 1) for n in A6C_TERMS]
        tj = {}
        for J in (1, 2):
            oj, bj = [], []
            for n in A6C_TERMS:
                Ks = a6c_win(k0, K1, n + J + 2)
                c = solve([pr(n)(K) + [qp[K] * P(K, j) for j in range(J + 1)] for K in Ks], [R[K] for K in Ks])
                oj.append(fl(c[0] * t - 1))
                bj.append([fl(c[n + 1 + j]) for j in range(J + 1)])
            tj[J] = dict(orders=oj, B0=[b_[0] for b_ in bj], b10=bj[-1][1] / bj[-1][0])
        Ks = a6c_win(k0, K1, 22)
        c = solve([[P(K, k) for k in range(2, 22)] + [qp[K], qp[K] * P(K, 1)] for K in Ks], [eps[K] for K in Ks])
        c2q = Fraction(int(c[0].p), int(c[0].q))
        z.update(tr=dict(K0=k0, orders=po, J=tj, q2K=float(qd ** (2 * k0))), c2fit=float(c2q),
                 c2rel=float(abs(c2q / Fraction(c2) - 1)) if c2 is not None else None)
    _A6C_MEMO[key] = z
    return z


def a6c_lines(A6C, A6=None, A6X=None):
    """A6C's line (referee round 129, D3dx-1, -4): the exact errors' provenance (Arb, the bits, the stored digits, the balls' radii) and
    their agreement with A6's and A6X's 26-digit errors at the K they share; per rule `a6c_analysis` — the orders on K ≥ A6C_K[0] with
    their verdict, the residual floor, the K⁻¹ coefficient, the excess's coefficients and their ratios c_k/c_{k−1}, the other branch's
    q^{A6C_K[0]}; at the transient's rule the power orders and those with q^K(B₀ + B₁/K) and with B₂/K² too, B₀ and B₁/B₀."""
    if A6C is None or len(A6C) == 0:
        return []
    from decimal import Decimal
    stored, rnd = {}, 0.0
    for xs in (A6 if A6 is not None else [], A6X if A6X is not None else []):
        for x in xs:
            for y in x['rows']:
                for K, s in ((int(y['K']) - 1, y.get('e0')), (int(y['K']), y.get('e1'))):
                    if s:
                        stored[(x['r'], x['name'], K)] = s
    xd, kx = 0.0, set()
    for x in A6C:
        for K, s in zip(range(x['K0'], x['K1'] + 1), x['e']):
            s0 = stored.get((x['r'], x['name'], K))
            if s0 is not None:
                xd = max(xd, abs(float(Fraction(s) / Fraction(s0) - 1)))
                rnd = max(rnd, 0.5 * 10.0 ** (1 - A6X_DIGITS) / abs(float(Decimal(s0).scaleb(-Decimal(s0).adjusted()))))
                kx.add(K)
    kr = sorted(K for K in kx if K - 1 in kx)
    bits = [b for x in A6C for b in x['bits']]
    parts = []
    for x in A6C:
        z = a6c_analysis(x)
        fr = abs(1 + 2 / x['r'] ** 2 - 2 / float(x['v'])) < 1 - 1e-12
        head = f"({x['r']:g}, {x['name']}; {'a frame' if fr else 'no frame'}; q = {z['q']:.6g}, q^{A6C_K[0]} = {z['qK0s']})"
        if z['tr'] is None:
            parts.append(head + f": the orders {', '.join(f'{v:+.1e}' for v in z['orders'])} — "
                         + ("decreasing in magnitude" if z['mono'] else "NOT decreasing in magnitude")
                         + (f", the last within ten times its last change ({up_e(z['res'])}; the input bound at most {up_e(max(z['bnd']))})"
                            if z['cons'] else ", the last NOT within ten times its last change")
                         + f"; the residual floor {', '.join(f'{v:.1e}' for v in z['floor'])}; the K⁻¹ coefficient {z['k1']:+.1e}; "
                         + (f"c_2 … c_{z['coef'][-1][0]} " + ', '.join(f"{c:+.6g}" for _, c, _ in z['coef']) + ", the ratios c_k/c_{k−1} "
                            + ', '.join(f"{b_[1] / a_[1]:+.1f}" for a_, b_ in zip(z['coef'], z['coef'][1:])) if z['coef'] else "no coefficient stable"))
        else:
            t_ = z['tr']
            parts.append(head + f": on [{t_['K0']}, {A6C_K[1]}] the power orders {', '.join(f'{v:+.1e}' for v in t_['orders'])}; with q^K(B₀ + B₁/K) "
                         f"in the basis {', '.join(f'{v:+.1e}' for v in t_['J'][1]['orders'])} (B₀ {', '.join(f'{v:.6f}' for v in t_['J'][1]['B0'])}; "
                         f"B₁/B₀ {t_['J'][1]['b10']:+.4f} at {a6c_nw(max(A6C_TERMS) + 1)} terms); with B₂/K² too {', '.join(f'{v:+.1e}' for v in t_['J'][2]['orders'])} "
                         f"(B₀ {', '.join(f'{v:.6f}' for v in t_['J'][2]['B0'])}; B₁/B₀ {t_['J'][2]['b10']:+.4f}: a small 1/K term in B) — the transient's own next term, "
                         f"q^{{2K}} = {t_['q2K']:.1e} at K = {t_['K0']}, unmodelled")
    out = ['', f"A6X's rules at consecutive K (A6C, referee round 129): the exact errors e_K at every K = {A6C[0]['K0']} … {A6C[0]['K1']} in Arb "
               f"(python-flint {A6C[0].get('flint', '?')}: the Gauss–Hermite nodes from LAPACK, refined by Newton in Arb; e_K² from the Gaussian "
               f"inner products, independent of A6's decimal arithmetic; ⌈2K·log₂(1/t(v))⌉ + {A6C_BITS} bits, {min(bits)}–{max(bits)}; stored to "
               f"{A6C_DIGITS} digits, the balls' relative radii at most {up_e(max(x['rad'] for x in A6C))}), within {up_e(xd)} (relative) of A6's and "
               f"A6X's {A6X_DIGITS}-digit errors at K − 1 and K for K = {', '.join(str(K) for K in kr)} (their rounding at most {up_e(rnd)}); from the "
               f"ratios in exact rationals, t(v) exact, the orders by windows on [{A6C_K[0]}, {A6C_K[1]}] (n + 1 K spread uniformly in 1/K, the "
               f"limit free; {a6c_nw(min(A6C_TERMS) + 1)} … {a6c_nw(max(A6C_TERMS) + 1)} terms), the residual floor (N powers K⁻² … K^{{−N−1}} of the "
               f"excess t(v)R_K − 1, the limit exact, through N such K; the largest residual over all {A6C_K[1] - A6C_K[0] + 1} K, N = "
               f"{', '.join(str(N) for N in A6C_FLOOR)}), the K⁻¹ coefficient (times t(v); beside 15 powers) and the excess's coefficients c_k "
               f"(the limit exact; the median of {len(A6C_COEF)} fits, {', '.join(sorted({str(N) for _, N in A6C_COEF}, key=int))} powers on "
               + ' and '.join(f"[{a_}, {A6C_K[1]}]" for a_ in sorted({a_ for a_, _ in A6C_COEF})) + ", printed while they agree within "
               f"{A6C_SPREAD * 100:g}%) — " + '; '.join(parts) + " — numerical, not a proof."]
    return out


def a6_lines(A6, A5=None):
    """A6's table (referee round 84, D3ce-3): per τ/σ and rule its variance over v_hs, the conjecture's 1/r (`rate_conj`) and the argument
    attaining it, the exact ratios per kernel at A6_K, the successive differences' ratios against K⁻²'s and K⁻¹'s, the limit — by
    L + c/K² + d/K³ through the last three K, at v_hs (the cusp) by L + c/K + d/K² — over 1/r, minus 1, and the ratio's excess over 1/r
    at the last K with its decay exponent from the K before; by the argument attaining r, the rules away from the cusp (|ln(v/v_hs)| ≥
    A6_NEAR) summarised and those near it listed; at v_hs the ratio over the sketch's (1 − 1/(4K))/t (e_K ∝ K^{1/4}t^K), minus 1, and
    its successive ratios against K⁻²'s; the Mercer frame's rule's exact squared errors over A5's frame series."""
    if A6 is None or len(A6) == 0:
        return []
    Ks = [float(K) for K in A6_K]
    dq = lambda p: ', '.join(f"{(Ks[i + 1] ** p - Ks[i + 2] ** p) / (Ks[i] ** p - Ks[i + 1] ** p):.3f}" for i in range(len(Ks) - 2))
    out = ['', "### A6. The scaled rules' rate: one formula (referee round 84; proved by D‴.61 (iv) where a rule has a Mehler frame"
           + f" ({D61_STATUS}); by D‴.62 (b) where it has none, without computer assistance (referee round 132))", '',
           "The Gauss–Hermite rule of N(0, v) for the reweighted mixing law (`scaled_rule`), its relative L² error computed exactly "
           f"(`a6_error`: the closed-form Gaussian inner products in decimal arithmetic of 2K·log₁₀(1/t) + {A6_MARGIN} digits, σ = 1), against "
           "the formula t(v) = max{h(τ²/v), h(vσ*²/(σ²τ²))}, h(y) = (y − 1)/(y + 1) — the two arguments are (σ*/σ)(v_hs/v) and "
           "(σ*/σ)(v/v_hs), v_hs = τ²σ/σ*, so t(v) = tanh((r + |ln(v/v_hs)|)/2), r = ln(σ*/σ): a mis-scaling of the variance by e^u adds |u| "
           "to the squeeze, and t(v) ≥ t with equality only at v_hs, a cusp —: the ratio per kernel e_{K−1}/e_K at K = " + ', '.join(f"{K}" for K in A6_K)
           + f"; the successive differences' ratios (as K⁻²: {dq(-2.0)}; as K⁻¹: {dq(-1.0)}); the limit by L + c/K² + d/K³ through the last "
           "three K (at v_hs by L + c/K + d/K²) against 1/t(v):", '',
           '| τ/σ | rule | v/v_hs | 1/t(v) (attained at) | e_{K−1}/e_K at K = ' + ', '.join(f"{K}" for K in A6_K)
           + f" | the differences' ratios | the limit | over 1/t(v), minus 1 | e_{{K−1}}/e_K·t(v) − 1 at K = {A6_K[-1]} (its exponent from "
             f"K = {A6_K[-2]}) |", '|---|---|---|---|---|---|---|---|---|']
    dev, cusp = [], []
    for x in A6:
        rr = [y['ratio'] for y in x['rows']]
        if any(v is None for v in rr):
            out.append(f"| {x['r']:g} | {x['name']} | {x['v_hs']:.6f} | {1 / x['rate']:.6f} | — | — | — | — | — |")
            continue
        at = x['name'] == '1·v_hs'
        p, q = (-1.0, -2.0) if at else (-2.0, -3.0)
        L = float(np.linalg.solve(np.array([[1.0, K ** p, K ** q] for K in Ks[-3:]]), rr[-3:])[0])
        d = np.diff(rr)
        y1, y0 = rr[-1] * x['rate'] - 1, rr[-2] * x['rate'] - 1
        ex = -float(np.log(y1 / y0) / np.log(Ks[-1] / Ks[-2])) if y1 * y0 > 0 else None
        L2 = float(np.linalg.solve(np.array([[1.0, K ** -1.0, K ** -2.0] for K in Ks[-3:]]), rr[-3:])[0])   # referee round 88 (D3ci-1)
        (cusp if at else dev).append((x, L * x['rate'] - 1, y1, ex, L2 * x['rate'] - 1))
        lab = 'both arguments, the cusp' if at else ('τ²/v' if x['arg'] == 1 else 'vσ*²/τ²')         # D3ci-5d: at v_hs the two tie
        out.append(f"| {x['r']:g} | {x['name']} | {x['v_hs']:.6f} | {1 / x['rate']:.6f} ({lab}) | " + ', '.join(f"{v:.6f}" for v in rr) + ' | '
                   + ', '.join(f"{d[i + 1] / d[i]:.3f}" for i in range(len(d) - 1)) + f" | {L:.6f} | {L * x['rate'] - 1:+.1e} | {y1:+.1e} ("
                   + ('—' if ex is None else f"{ex:.2f}") + ') |')
    if dev:
        w = max(dev, key=lambda z: abs(z[1]))
        dist = lambda zs: len({(z[0]['r'], round(float(z[0]['v']), 12)) for z in zs}) if all('v' in z[0] for z in zs) else len(zs)
        dst_ = lambda zs: f" ({dist(zs)} distinct)" if dist(zs) < len(zs) else ''     # round 107 (D3db-4): Fock ≡ Mercer-zero at τ/σ = 2
        out += ['', f"Off the cusp the limit lies within {up_e(w[1])} of the formula's 1/t(v) at all {len(dev)} rules{dst_(dev)} (the farthest: "
                    f"τ/σ = {w[0]['r']:g}, {w[0]['name']})" + (f"; at v_hs within {up_e(max(abs(z[1]) for z in cusp))} of 1/t" if cusp else '')
                + '.' + (lambda w2: f" With a K⁻¹ term instead (L + b/K + c/K² through the same K) off the cusp within {up_e(w2[4])} "
                         f"(the farthest: τ/σ = {w2[0]['r']:g}, {w2[0]['name']}): the second fit's largest deviation; the two fits differ by at most "
                         f"{up_e(max(abs(z[1] - z[4]) for z in dev))} (bounds rounded up: referee round 96).")(max(dev, key=lambda z: abs(z[4])))]
        fr_ = lambda z: abs(1 + 2 / z[0]['r'] ** 2 - 2 / float(z[0]['v'])) < 1 - 1e-12 if 'v' in z[0] else False   # referee round 92 (D3cm-1)
        for lab_, zs in (("with a frame, whose limits D‴.61 (iv) proves (the fits' deviations their bias)",
                          [z for z in dev if fr_(z)]),
                         ("without a frame, whose limits D‴.62 (b) proves by hand (the fits' deviations their bias too; 'the proved rules' "
                          "below are the frame rules, the calibration predating D‴.62)", [z for z in dev if not fr_(z)])):
            if zs:
                out[-1] += (f" Off the cusp, the {len(zs)} rules{dst_(zs)} {lab_}: the first fit within {up_e(max(abs(z[1]) for z in zs))}, the second "
                            f"within {up_e(max(abs(z[4]) for z in zs))}, the two bracketing 1/t(v) at {sum(1 for z in zs if z[1] * z[4] <= 0)} of "
                            f"the {len(zs)}.")
        cells = []                                                      # referee round 96 (D3cq-2): by frame, branch and nearness
        for fl, fz in (('with a frame', True), ('without', False)):
            for a_, bl in ((1, 'τ²/v'), (2, 'vσ*²/τ²')):
                for nl, nf in (('far', True), ('near', False)):
                    zs = [z for z in dev if fr_(z) == fz and z[0]['arg'] == a_ and (abs(np.log(z[0]['v_hs'])) >= A6_NEAR) == nf]
                    if zs:
                        cells.append(f"{fl}, {bl}, {nl} ({len(zs)}" + (f", {dist(zs)} distinct" if dist(zs) < len(zs) else '')
                                     + f"): within {up_e(max(abs(z[1]) for z in zs))} and "
                                     f"{up_e(max(abs(z[4]) for z in zs))}, bracketing at {sum(1 for z in zs if z[1] * z[4] <= 0)}")
        out[-1] += (f" By frame, the branch attaining t(v) and nearness to the cusp (far: |ln(v/v_hs)| ≥ {A6_NEAR:g}; referee round 96): "
                    + '; '.join(cells) + " — the fits' bracket calibrated by the proved limits only where they bracket those.")
        fz = lambda fr: [z for z in dev if fr_(z) == fr and z[0]['arg'] == 2 and abs(np.log(z[0]['v_hs'])) >= A6_NEAR]   # round 101 (D3cv-1)
        if fz(True) and fz(False):
            w1, w2 = (max(fz(False), key=lambda z: abs(z[j])) for j in (1, 4))
            q1, q2 = (abs(w[j]) / max(abs(z[j]) for z in fz(True)) for w, j in ((w1, 1), (w2, 4)))
            nm_ = lambda w: f"τ/σ = {w[0]['r']:g}, {w[0]['name']}"
            out[-1] += (" On the fast branch far from the cusp the largest deviations without a frame, " + (f"both at {nm_(w1)}, are" if w1 is w2 else
                        f"at {nm_(w1)} and {nm_(w2)}, are") + f" {dn_f(q1, 1)} and {dn_f(q2, 1)} times the proved rules' largest (the first fit, "
                        "the second; referee round 101).")
        sf = sorted([z for z in dev if z[0]['arg'] == 1 and abs(np.log(z[0]['v_hs'])) >= A6_NEAR and 't' in z[0]], key=lambda z: 1 / z[0]['rate'])
        if sf:                                                          # round 101 (D3cv-2): the slow branch's rules by t(v), frame or not
            out[-1] += (" The slow branch far from the cusp by t(v) (both fits' deviations): " + '; '.join(
                f"t(v) {z[0]['rate']:.3f} at τ/σ = {z[0]['r']:g}, {z[0]['name']} ({'a frame' if fr_(z) else 'no frame'}): {z[1]:+.2e}, {z[4]:+.2e}" for z in sf) + '.')
        rng = lambda vs, f: f"{f(min(vs))} … {f(max(vs))}" if vs else '—'
        rgo = lambda vs: f"{out_e(min(vs), False)} … {out_e(max(vs), True)}" if vs else '—'   # round 101 (D3cv-3): ranges outward
        far = lambda a_: [z for z in dev if z[0]['arg'] == a_ and abs(np.log(z[0]['v_hs'])) >= A6_NEAR]
        near = [z for z in dev if abs(np.log(z[0]['v_hs'])) < A6_NEAR]
        out += ['', f"By the argument attaining t(v), away from the cusp (|ln(v/v_hs)| ≥ {A6_NEAR:g}): " + '; '.join(
            f"{lab} ({len(far(a_))} rules): the limit over 1/t(v), minus 1, " + rgo([z[1] for z in far(a_)])
            + f"; the ratio's excess at K = {A6_K[-1]} " + rgo([z[2] for z in far(a_)]) + ", its exponent "
            + (lambda vs: f"{np.floor(min(vs) * 100) / 100:.2f} … {np.ceil(max(vs) * 100) / 100:.2f}" if vs else '—')(
                [z[3] for z in far(a_) if z[3] is not None])                # round 107 (D3db-3): outward
            + (f" ({sum(1 for z in far(a_) if z[3] is None)} changing sign)" if any(z[3] is None for z in far(a_)) else '')
            + "; with a K⁻¹ term " + rgo([z[4] for z in far(a_)])
            + f", the two fits bracketing 1/t(v) at {sum(1 for z in far(a_) if z[1] * z[4] <= 0)} of the {len(far(a_))}"
            for a_, lab in ((1, 'τ²/v'), (2, 'vσ*²/τ²'))) + '. Near it: ' + (', '.join(
                f"τ/σ = {z[0]['r']:g}, {z[0]['name']} (v/v_hs {z[0]['v_hs']:.6f}): {z[2]:+.1e}" for z in near) if near else 'none') + '.']
    if cusp:
        cs = []
        for x, *_ in cusp:
            res = [y['ratio'] * x['rate'] / (1 - 1 / (4 * K)) - 1 for y, K in zip(x['rows'], Ks)]
            cs.append(f"τ/σ = {x['r']:g}: " + ', '.join(f"{v:+.1e}" for v in res) + ' (' + ', '.join(f"{res[i + 1] / res[i]:.3f}"
                                                                                                 for i in range(len(res) - 1)) + ')')
        out += ['', "At v_hs the sketch's integral spreads over every angle, so e_K ∝ K^{1/4}t^K and the ratio per kernel is "
                    "(1 − 1/(4K))/t + O(K⁻²): the ratio over (1 − 1/(4K))/t, minus 1, at K = " + ', '.join(f"{K}" for K in A6_K)
                + " (its successive ratios; K⁻²'s " + ', '.join(f"{(Ks[i] / Ks[i + 1]) ** 2:.3f}" for i in range(len(Ks) - 1)) + '): '
                + '; '.join(cs) + '.']
    if all('v' in x for x in A6):                                       # referee round 88 (D3ci-3): the rules D‴.61 (iv) covers
        rho = [(x, 1 + 2 / x['r'] ** 2 - 2 / float(x['v'])) for x in A6]  # σ = 1
        fin, one = [x for x, p in rho if abs(p) < 1 - 1e-12], [x for x, p in rho if abs(p - 1) <= 1e-12]
        dst = lambda xs: len({(x['r'], round(float(x['v']), 12)) for x in xs})   # round 107 (D3db-4): Fock ≡ Mercer-zero at τ/σ = 2
        out += ['', f"A Mehler frame (|ρ| < 1, ρ = 1 + 2σ²/τ² − 2σ²/v: v ∈ (τ²σ²/σ*², τ²); D‴.61) exists for {len(fin)} of the {len(A6)} "
                    f"rules ({dst(fin)} of {dst(A6)} distinct; " + ', '.join(f"{sum(1 for x in fin if x['r'] == r_)} at τ/σ = {r_:g}" for r_ in dict.fromkeys(x['r'] for x in A6))
                + f"); ρ = 1 at {len(one)} (" + ', '.join(dict.fromkeys(x['name'] for x in one)) + f"), |ρ| > 1 at {len(A6) - len(fin) - len(one)}."]
    if all('v' in x for x in A6):                                       # referee round 92 (D3cm-2): the rules with a frame at every τ/σ
        import math
        def cross(f, rs=D61_RSTAR, lo=1e-3, hi=50.0):
            for _ in range(100):
                m = (lo + hi) / 2
                lo, hi = (m, hi) if f(m) < rs else (lo, m)
            return lo
        fock = lambda r_: rate_conj(r_, 2 * r_ * r_ / (2 + r_ * r_))[0]
        mz = lambda r_: rate_conj(r_, r_ * r_ / math.sqrt(1 + 2 * r_ * r_))[0]
        rf, rm = cross(fock), cross(mz)
        rfv, rmv = cross(fock, D61_RSTAR_V, hi=1e3), cross(mz, D61_RSTAR_V, hi=1e3)   # round 101: with D‴.61's part V (an earlier tier)
        out += ['', "The Fock rule (ρ = 0) and the Mercer-zero rule have a frame at every τ/σ, so with D‴.61's part VA (its (iv) at every "
                    f"frame, r* = {D61_RSTAR_VA:g}; {D61_STATUS}; referee rounds 101, 107) their rates — λ/(2 − λ) and q — are proved at "
                    f"every τ/σ, no threshold left. The earlier tiers (thresholds rounded inward: referee round 96): with part V (r* = "
                    f"{D61_RSTAR_V:.4f}) while below that r*, σ*/σ < {dn_f(math.sqrt(1 + rfv * rfv))} (τ/σ < {dn_f(rfv)}) and τ/σ < {dn_f(rmv)} — "
                    f"near the frame interval's ends t(v) → λ, above that r* once σ*/σ > {up_f(math.sqrt((1 + D61_RSTAR_V) / (1 - D61_RSTAR_V)))} —; "
                    f"with the branch and bound alone (r* = {D61_RSTAR:g}): σ*/σ < {dn_f(math.sqrt(1 + rf * rf))} (τ/σ < {dn_f(rf)}), "
                    f"τ/σ < {dn_f(rm)}, the frame ends above it once σ*/σ > {up_f(math.sqrt((1 + D61_RSTAR) / (1 - D61_RSTAR)))}."]
    if A5 is not None and A5.get('series'):
        sr = {y['r']: {z['K']: z['S'] for z in y['frames']['Mercer']} for y in A5['series']}
        cs = []
        for x in A6:
            if x['name'] != 'Mercer frame' or x['r'] not in sr:
                continue
            e2 = {}
            for y in x['rows']:
                for K_, e_ in ((y['K'] - 1, y['e0']), (y['K'], y['e1'])):
                    if e_ is not None and K_ in sr[x['r']]:
                        e2[K_] = float(e_) ** 2 / sr[x['r']][K_]
            if len(e2) > 1:
                k0 = min(e2)
                cs.append((x['r'], max(abs(v / e2[k0] - 1) for v in e2.values()), sorted(e2)))
        if cs:
            out += ['', "The Mercer frame's rule's exact squared error over its frame series (A5) is one constant per τ/σ: within "
                        + ', '.join(f"{up_e(c)} at τ/σ = {r_:g} (K = {', '.join(str(k) for k in ks)})" for r_, c, ks in cs)
                        + ' — the resolution of the stored doubles, not an agreement measured below it.']
    return out


A6E_K, A6E_S = (20, 40, 80), tuple(i / 10 for i in range(1, 21))     # A6E (referee round 84): K, and |ξ|/√K on each axis


def dpi():
    """π in the current decimal context (Machin's 16 arctan(1/5) − 4 arctan(1/239))."""
    from decimal import Decimal, getcontext
    eps = Decimal(10) ** -(getcontext().prec + 2)

    def atan_inv(n):
        x = Decimal(1) / n
        x2, term, tot, k = x * x, x, x, 1
        while abs(term) > eps:
            term *= -x2
            tot += term / (2 * k + 1)
            k += 1
        return tot
    return 16 * atan_inv(5) - 4 * atan_inv(239)


def dcos(t, pi):
    """cos t in the current decimal context (t reduced modulo 2π, then its Taylor series)."""
    from decimal import Decimal, getcontext
    t = t % (2 * pi)
    eps = Decimal(10) ** -(getcontext().prec + 2)
    term, tot, n = Decimal(1), Decimal(1), 0
    while abs(term) > abs(tot) * eps:
        n += 2
        term *= -t * t / (n * (n - 1))
        tot += term
    return tot


def gen_ratio(K, s, axis, prec):
    """The K-point Gauss–Hermite rule's defect generating function ε̃_K(ξ) = 1 − e^{−ξ²/2}Q_K[e^{ξz}] (Σ_n D_n(v′)sⁿ/n! =
    e^{(v′−1)s²/2}ε̃_K(√v′s)) over (K!/(2K)!)ξ^{2K}e^{−ξ²/4}, at ξ = s√K ('real') or ξ = is√K ('imag', where Q_K[e^{ξz}] =
    Q_K[cos(|ξ|z)]), in prec-digit decimal arithmetic (referee round 84, D3ce-3)."""
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = prec
        zs, ws = gh_nodes(K, prec)
        x = Decimal(repr(s)) * Decimal(K).sqrt()
        a = Decimal(factorial(K)) / factorial(2 * K) * x ** (2 * K)
        if axis == 'real':
            e = 1 - (-(x * x) / 2).exp() * sum((w * (x * z).exp() for w, z in zip(ws, zs)), Decimal(0))
            return float(e / (a * (-(x * x) / 4).exp()))
        pi = dpi()
        e = 1 - ((x * x) / 2).exp() * sum((w * dcos(x * z, pi) for w, z in zip(ws, zs)), Decimal(0))
        return float(e / ((-1) ** K * a * ((x * x) / 4).exp()))


def part_A6E():
    """A6E (referee round 84, D3ce-3): `gen_ratio` at K in A6E_K and s in A6E_S on both axes, in 3K + 60 digits."""
    return [dict(K=K, axis=ax, s=list(A6E_S), ratio=[gen_ratio(K, s, ax, 3 * K + 60) for s in A6E_S]) for K in A6E_K for ax in ('real', 'imag')]


def a6e_lines(A6E):
    """A6E's line: per K and axis the generating function's ratio to (K!/(2K)!)ξ^{2K}e^{−ξ²/4} — its range over |ξ| ≤ 2√K and its values
    at |ξ|/√K = 0.5, 1, 1.5, 2."""
    if A6E is None or len(A6E) == 0:
        return []
    at = lambda x, s: x['ratio'][list(x['s']).index(s)]
    return ['', "The step of the rate's Laplace sketch (A6E): the rule's defect generating function ε̃_K(ξ) = 1 − e^{−ξ²/2}Q_K[e^{ξz}] — "
                "Σ_n D_n(v′)sⁿ/n! = e^{(v′−1)s²/2}ε̃_K(√v′s), so that the frame series is its Fock norm, (πρv′)⁻¹∫|ε̃_K(ξ)|²e^{−Ax²−By²}dA, "
                "A = (1 − ρ(v′ − 1))/(ρv′), B = (1 + ρ(v′ − 1))/(ρv′) — over (K!/(2K)!)ξ^{2K}e^{−ξ²/4}, on ξ = s√K and ξ = is√K for "
                f"s = {A6E[0]['s'][0]:g} … {A6E[0]['s'][-1]:g} (`gen_ratio`, 3K + 60 digits): " + '; '.join(
                    f"K = {x['K']}, {'the real axis' if x['axis'] == 'real' else 'the imaginary axis'}: {min(x['ratio']):.3f}–{max(x['ratio']):.3f} "
                    "(at s = 0.5, 1, 1.5, 2: " + ', '.join(f"{at(x, s):.3f}" for s in (0.5, 1.0, 1.5, 2.0)) + ')' for x in A6E)
            + ". Bounded above at these K; on the real axis its minimum falls with K toward D‴.61's limit (1 − s²/4)^{1/2}, which vanishes at "
              "s = 2 (referee round 88), so its K-th root tends to 1 where s < 2; Laplace's method on the Fock integral — |ξ|^{4K}e^{−(A + ½)x² − "
              "(B − ½)y²}, its mass at |ξ| = (2K/min(A + ½, B − ½))^{1/2} — then gives the rate 1/(2 min(A + ½, B − ½)), the formula's "
              f"max{{h(τ²/v), h(vσ*²/(σ²τ²))}} (a sketch here; D‴.61 (iii)–(iv) prove it for every rule with a Mehler frame, {D61_STATUS}, "
              "and D‴.62 (b) for every rule without one, by hand)."]


def part_B2(B):
    """B2: `b2_job` from every non-optimal block of B."""
    from multiprocessing import Pool
    jobs = [(b['r'], b['var'], o['k'], np.asarray(o['c'])) for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL]
    with Pool(N_WORKERS) as pool:
        return pool.map(b2_job, jobs, chunksize=1)


def b2_line(B2):
    """B2's summary: the non-optimal blocks from which variable projection reaches the best optimum found (within OPT_TOL), counted where
    that optimum is at least E_FLOOR — below it the solve stops short of the optimum (A), and a miss says nothing of the block's basin
    (referee round 60); of the misses, those within 1% of the optimum (in its basin, unconverged) and the others with the block's largest
    distance from the optimum's nodes (σ); the same for each path's first non-optimal block."""
    if B2 is None or not len(B2):
        return []
    best = [min(q['e_block'], q['e_rule']) for q in B2]
    hit = [abs(q['e_block'] / b - 1) < OPT_TOL for q, b in zip(B2, best)]
    up = [b >= E_FLOOR for b in best]
    miss = [(q['e_block'] / b, q.get('dc', np.nan) / SIG, q) for q, b, h, u in zip(B2, best, hit, up) if u and not h]
    near = sorted(v for v, d, q in miss if v < 1.01)
    far = [(v, d, q) for v, d, q in miss if v >= 1.01]
    seen, first = set(), []                                            # each path's first non-optimal block (B2 keeps the paths' order)
    for i, q in enumerate(B2):
        if (q['r'], q['var']) not in seen:
            seen.add((q['r'], q['var']))
            first.append(i)
    fu = [i for i in first if up[i]]
    fw = [B2[i]['e_block'] / best[i] for i in fu if not hit[i]]
    ends = lambda n: f"{n} {'ends' if n == 1 else 'end'}"
    span = lambda v, f_: f"{min(v):{f_}}" if min(v) == max(v) else f"{min(v):{f_}}–{max(v):{f_}}"
    nb = len(first) - len(fu)
    sp_ = lambda k, f_: span([B2[i][k] for i in fu], f_)
    md_ = lambda k: f"{np.median([B2[i][k] for i in fu]):.3g}"
    return [f"Variable projection from the non-optimal block's own centres (B2) reaches the best optimum found in "
            f"{sum(h for h, u in zip(hit, up) if u)} of the {sum(up)} blocks whose optimum is at least {E_FLOOR:g}"
            + (f"; {ends(len(near))} within 1% of it ({', '.join(f'{v:.4g}' for v in near)} times: in its basin, unconverged)" if near else '')
            + (f"; {ends(len(far))} {span([v for v, d, q in far], '.4g')} times above it — the blocks' largest node distances from the "
               f"optimum's {span([d for v, d, q in far], '.2f')}σ —" + (
                   f" where variable projection stalls: its two closest nodes {span([q['vp_gap'] for v, d, q in far], '.3f')}σ apart, its "
                   f"end not stationary (the largest cosine of its residual with a Jacobian column "
                   f"{span([q['vp_cos'] for v, d, q in far], '.1e')}), its stop "
                   f"{', '.join(f'{q['vp_status']}/{q['vp_nfev']}' for v, d, q in far)} (status/evaluations: a miss says nothing of the "
                   f"block's basin)" if 'vp_gap' in far[0][2] else '') if far else '')
            + f". The {len(B2) - sum(up)} blocks with an optimum below {E_FLOOR:g}, where the solve stops short of the optimum (A), are not "
            f"counted ({sum(h for h, u in zip(hit, up) if not u)} reached). From the first non-optimal block of each path: "
            f"{sum(hit[i] for i in fu)} of {len(fu)}" + (f" (the others {', '.join(f'{v:.4g}' for v in sorted(fw))} times above it)" if fw else '')
            + (f"; on {nb} {'path' if nb == 1 else 'paths'} the first non-optimal block's optimum is below {E_FLOOR:g}" if nb else '') + "."
            + (f" There the Jacobian's condition number is {sp_('cond_joint_x', '.2g')} for the joint fit in (w, c) with the centres in x, "
               f"the polish's parameters ({sp_('cond_joint', '.2g')} with the centres in units of σ) and {sp_('cond_vp_x', '.2g')} for variable projection's "
               f"in c alone (analytic; the same in either unit of the centres); the ratio of their smallest singular values, variable "
               f"projection's over the joint fit's, is at least 1 in any common scaling of the centres for Kaufman's form of the former, "
               f"the centre columns projected off the amplitudes' (referee round 77: the joint fit's smallest singular value is at most "
               f"theirs — a Schur complement), from which the analytic Golub–Pereyra form used here differs by a term proportional to the "
               f"residual; its size depends on the units: {sp_('smin_ratio', '.3g')} (median {md_('smin_ratio')}) with the "
               f"centres in units of σ, {sp_('smin_ratio_x', '.3g')} (median {md_('smin_ratio_x')}) in x, "
               f"{sp_('smin_ratio_eq', '.3g')} (median {md_('smin_ratio_eq')}) with every column scaled to unit norm, which scales the two "
               f"problems' centre columns differently; TRF itself works in Coleman–Li's variables, the bounded amplitudes scaled by "
               f"v^{{1/2}}." if fu and 'smin_ratio_eq' in B2[fu[0]] else ''),
            '']


def cs_margin(r_, var, o, gtol=1e-9):
    """The stationarity test's threshold over its Cauchy–Schwarz bound — an effective cosine (referee round 65): the test is absolute,
    |⟨rφ_i⟩| and |w_i⟨r∂_cφ_i⟩| at most gtol·max_j|⟨f*φ_j⟩| (the starting grid's correlations, gtol the path's default), and by
    Cauchy–Schwarz |⟨rφ_i⟩| ≤ rms(r)·rms(φ_i), |⟨r∂_cφ_i⟩| ≤ rms(r)·rms(∂_cφ_i): the test passes whatever the state once the ratio is ≥ 1,
    and below 1 it asks the residual's cosine with the block's directions to be below the ratio."""
    shift = next(v for v in VARIANTS_B if v[0] == var)[2]
    K = RT.Bumps(X_B, SIG)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))
    scale = float(np.abs(K.phi(TM.setup()[2] + shift).T @ f).max()) / X_B.size
    c, w = np.asarray(o['c'], float), np.asarray(o['w'], float)
    lim = np.maximum(np.sqrt(np.mean(K.phi(c) ** 2, axis=0)), np.abs(w) * np.sqrt(np.mean(K.dphi(c) ** 2, axis=0)))
    return float(gtol * scale / (o['e'] * float(np.sqrt(np.mean(f ** 2))) * lim.max()))


def cs_vacuous(r_, var, o, gtol=1e-9):
    """Whether a block passes `sliding_iss`'s stationarity test whatever its state (referee round 60): `cs_margin` ≥ 1."""
    return bool(cs_margin(r_, var, o, gtol) >= 1)


def gap_share(x):
    """A3: the share of the half-squeeze rule's gap to the optimal fit that the best one-parameter variance recovers; an optimum the
    one-parameter rule beats is unconverged (referee round 54)."""
    ref = min(x['e_opt'], x['e_best'])
    sh = (x['e_hs'] - x['e_best']) / (x['e_hs'] - ref) if x['e_hs'] > ref else np.nan
    return f"{sh:.1%}" + (', optimum unconverged' if x['e_best'] < x['e_opt'] * (1 - UNCONV) else '')


def gap_range(A):
    """The range of `gap_share` over A's rows with an optimum, and how many optima the variance beats."""
    xs = [x for x in A['rows'] if 'e_opt' in x and x['e_hs'] > min(x['e_opt'], x['e_best'])]
    if not xs:
        return '—'
    sh = [(x['e_hs'] - x['e_best']) / (x['e_hs'] - min(x['e_opt'], x['e_best'])) for x in xs]
    return f"{min(sh):.1%}–{max(sh):.1%} ({sum(x['e_best'] < x['e_opt'] * (1 - UNCONV) for x in xs)} optima unconverged)"


def polish_cell(o, rf=None):
    """A block's polishes: status/evaluations/accepted (status 0: the budget spent first; referee round 54), and for a repeated polish
    (its end's rms recorded) that end's error over the optimum's in brackets (referee round 65: a gradient-test stop at or above it;
    B3 restarts the solve from the block)."""
    return ' '.join(f"{p['status']}/{p['nfev']}/{'a' if p['accepted'] else 'r'}"
                    + (f" ({p['e1'] / (o['e_opt'] * rf):.4g})" if rf and 'e1' in p else '') for p in o.get('polish', [])) or '—'


def entry_lines(B):
    """B's entries (`sliding_iss`'s entry_info) by the block they join (referee round 60: the regimes apart) — none (a path's first entry:
    every dual zero, so the least time to the wall is the largest rate), an optimal block, a non-optimal one (block it − 1, the entries in
    the order of the events) — and over all: whether the entrant had the largest current rate, its rate over the largest, its dual before
    the entry (the entrant is the kernel of the least time to the wall (1 − |u_j|)/|⟨rφ_j⟩|: referee round 54), and where it starts
    against the block it joins: beyond its outermost node (how far) or between its nodes (their duals)."""
    rows = {'none (the first entry)': [], 'an optimal block': [], 'a non-optimal block': []}
    for b in B:
        info = b.get('entries') if b.get('entries') is not None else []
        shift = next(v for v in VARIANTS_B if v[0] == b['var'])[2]
        for e, ev in zip(info, [ev for ev in b['events'] if ev[2] == 'enter']):
            key = (b['r'], shift, e['it'], e['j'])
            if e['it'] == 0:
                rows['none (the first entry)'].append((e, None, key))
            elif e['it'] <= len(b['blocks']):
                o = b['blocks'][e['it'] - 1]
                nodes = np.asarray(o['c'])
                rows['an optimal block' if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL else 'a non-optimal block'].append(
                    (e, max(nodes.min() - ev[3], ev[3] - nodes.max()) / SIG, key))
    rows['all'] = [x for v in list(rows.values()) for x in v]
    if not rows['all']:
        return []
    groups = {}
    for e, d, key in rows['all']:
        groups.setdefault(key, []).append(abs(e['u0']))
    cid, spread = {}, 0.0                                              # an event: one key, its duals within 1e-4 (the variants part)
    for key, us in groups.items():
        us_, c_, lo_ = sorted(set(us)), 0, None
        for i, u in enumerate(us_):
            if i and u - us_[i - 1] > 1e-4:
                spread, c_, lo_ = max(spread, us_[i - 1] - lo_), c_ + 1, None
            lo_ = u if lo_ is None else lo_
            cid[(key, u)] = c_
        spread = max(spread, us_[-1] - lo_)
    out = ['', "The entries over the variants that record them, by the block they join — one event per width, starting grid, iteration and "
           "kernel whose duals agree within 10⁻⁴ — the same entrant at the same iteration, whether or not the variants' paths have parted "
           f"(the entries over the variants in brackets; the duals of one event agree to {spread:.1e}; referee rounds 65, 70) —: the entrant "
           "had the largest current rate |⟨rφ_j⟩|; its rate over "
           "the largest; its dual |u_j| at the block's start; whether it starts beyond the block's outermost node (how far) or between its "
           "nodes (their duals). The entrant is the kernel of the least time to the wall (1 − |u_j|)/|⟨rφ_j⟩| (`sliding_iss`'s rule), its "
           "dual charged during the earlier blocks.", '',
           '| joining | events (entries) | largest rate | rate / largest: median (range) | dual: median (least) | beyond: entrants, median '
           '(range) | between: entrants (duals) |', '|---|---|---|---|---|---|---|']
    counts = {}
    for name, xs_all in rows.items():
        if not xs_all:
            continue
        seen_, xs = set(), []
        for e, d, key in xs_all:
            k_ = (key, cid[(key, abs(e['u0']))])
            if k_ not in seen_:
                seen_.add(k_)
                xs.append((e, d))
        rel = np.array([abs(e['rate']) / e['rate_max'] for e, d in xs])
        u0 = np.array([abs(e['u0']) for e, d in xs])
        top = int(sum(abs(abs(e['rate']) - e['rate_max']) <= 1e-12 * e['rate_max'] for e, d in xs))
        out_ = [d for e, d in xs if d is not None and d > 0]
        ins = [abs(e['u0']) for e, d in xs if d is not None and d <= 0]
        counts[name] = len(xs)
        out.append(f"| {name} | {len(xs)} ({len(xs_all)}) | {top} | {np.median(rel):.2f} ({rel.min():.3f}–{rel.max():.2f}) | {np.median(u0):.3f} "
                   f"({u0.min():.3f}) | " + (f"{len(out_)}, {np.median(out_):.2f}σ ({min(out_):.2f}–{max(out_):.2f}σ)" if out_ else '—')
                   + ' | ' + (f"{len(ins)} ({min(ins):.3f}–{max(ins):.3f})" if ins else '—') + ' |')
    ov = sum(v for k, v in counts.items() if k != 'all') - counts['all']
    if ov:                                                             # one event, an optimal block in some variants, not in others
        out += ['', f"{ov} events join an optimal block in some variants and a non-optimal one in others: counted in both rows, once in "
                    "'all'."]
    return out + ['']


def c_job(args):
    s, q = args
    warnings.filterwarnings('ignore')
    t0 = time.time()
    L, _, c0 = TM.setup()
    K = RT.Bumps(X_B, SIG)
    ss = sigma_star(TAU_C)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    xi = s * np.random.default_rng(5800 + 100 * LEVELS_C.index(s) + q).standard_normal(X_B.size)
    y = f + xi
    R = RT.sliding_iss(L, K, y, c0, np.ones(c0.size), polish=True, max_events=3 * K_C, tol=0.0, max_entries=K_C + 1,
                       max_doublings=MAX_DOUBLINGS, box=(X_B[0] - BOX_C * SIG, X_B[-1] + BOX_C * SIG),
                       w_scale=float(np.abs(y).max()), exit_kkt=True)
    rf = float(np.sqrt(np.mean(f ** 2)))
    bl = [b for b in blocks(R) if b['k'] > 0]
    ent = [e for e in R['events'] if e[2] == 'enter']
    taus = [e[0] for e in ent]
    out = dict(s=s, q=q, k=[b['k'] for b in bl], rms=[b['rms'] for b in bl],
               err=[float(np.sqrt(np.mean((K.phi(b['c']) @ b['w'] - f) ** 2)) / rf) for b in bl],
               at=[int(np.sum(np.abs(b['c'] - C_B) < 4 * ss)) for b in bl],
               gaps=[float((taus[i + 1] - taus[i]) / taus[i]) for i in range(len(taus) - 1)], exits=sum(e[2] == 'exit' for e in R['events']),
               stiff=R['stiff_redone'], runaway=R['runaway'], relaxed=[b['relaxed'] for b in bl], rf=rf, secs=time.time() - t0)
    return out


def part_C():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, [(s, q) for s in LEVELS_C for q in range(N_DRAW_C)], chunksize=1)


def d_job(args):
    """D: the tail ratio at (K, t), with enough terms for t^{2j} to fall below 10⁻¹⁷ (40/ln(1/t) of them, at least 40)."""
    K, t = args
    n = max(40, int(np.ceil(40 / np.log(1 / t))))
    return dict(K=K, t=t, ratio=W.tail_ratio(K, t, n_terms=n), n_terms=n)


def lead_coeffs(j_max):
    """The ratios P_j(K) = D_{2K+2j,K}/D_{2K,K} as polynomials in K: the j-th divided difference over K = 1 … j + 1 (the leading
    coefficient, exact) and the (j + 1)-th over K = 1 … j + 2 (zero for degree j)."""
    from fractions import Fraction

    def dd(ps):
        return ps[0][1] if len(ps) == 1 else (dd(ps[1:]) - dd(ps[:-1])) / (ps[-1][0] - ps[0][0])
    out = []
    for j in range(1, j_max + 1):
        pts = [(K, Fraction(W.gh_defects(K, 2 * K + 2 * j)[2 * K + 2 * j], W.gh_defects(K, 2 * K)[2 * K])) for K in range(1, j + 3)]
        lead = dd(pts[:j + 1])
        out.append(dict(j=j, lead=(lead.numerator, lead.denominator), next_dd=float(dd(pts))))
    return out


def part_D():
    from multiprocessing import Pool
    ts = tuple(float(W.params(SIG, sigma_star(r_))['t']) for r_ in (1.0, 2.0)) + (0.5, 0.8)
    with Pool(N_WORKERS) as pool:
        rows = pool.map(d_job, [(K, t) for t in ts for K in KS_D], chunksize=1)
    return dict(rows=rows, lead=lead_coeffs(J_MAX_D))


def flank(r_, c, w):
    """B: the residual's correlation with a kernel at c, C(c) = ⟨r, φ_σ(· − c)⟩ on the samples, after a block (c, w) on the feature of
    τ/σ = r_, on a fine grid of c: min C/max C (zero or positive where C ≥ 0) and its largest values' distances beyond the outermost
    nodes, left and right (in σ)."""
    ss = sigma_star(r_)
    f = A_B * np.exp(-(X_B - C_B) ** 2 / (2 * ss ** 2))
    r = f - np.exp(-(X_B[:, None] - np.asarray(c)[None, :]) ** 2 / (2 * SIG ** 2)) @ np.asarray(w)
    cg = np.linspace(X_B[0], X_B[-1], 3201)
    C = np.exp(-(X_B[:, None] - cg[None, :]) ** 2 / (2 * SIG ** 2)).T @ r / X_B.size
    lo, hi = float(np.min(c)), float(np.max(c))
    L, R = cg < lo, cg > hi
    return float(C.min() / C.max()), float((lo - cg[L][np.argmax(C[L])]) / SIG), float((cg[R][np.argmax(C[R])] - hi) / SIG)


def predict_stop(s, rf, N):
    """(i)'s prediction of the discrepancy stop: the first K with e_K²rf² ≤ s²((1 + δ)² − 1 + 2K/N), and the error there."""
    ss = sigma_star(TAU_C)
    for K in range(1, 40):
        e = W.l2_error(*W.half_squeeze(K, SIG, ss), SIG, ss)
        if (e * rf) ** 2 <= s ** 2 * ((1 + DELTA) ** 2 - 1 + 2 * K / N):
            return K, float(np.sqrt(e ** 2 + 2 * K * s ** 2 / (N * rf ** 2)))
    return -1, np.nan


def stop_law(s, rf, N, k_max=12):
    """(iv)'s stop as a margin (referee round 54). The rms of the K-kernel fit to y = f* + ξ has the expectation
    m_K = s(b_K² + 1 − 2K/N)^{1/2}, b_K = e_K‖f*‖_N/s (2K parameters), and the scatter sd(rms²)/s² = (2(N − 2K)/N² + 4b_K²/N)^{1/2} (the
    noise's chi-square and its cross term with the bias); the rms falls along the path, so the rule stops at K with probability
    P(rms_K ≤ (1 + δ)s) − P(rms_{K−1} ≤ (1 + δ)s), in the normal approximation; e_K the half-squeeze rule's L² error, as in
    `predict_stop`. Also the first-term formula's K(s): the root of K = log((1 − t²)^{1/4}(πK)^{1/4}‖f*‖_N/(s(2δ + δ²)^{1/2}))/log(1/t)
    → (K(s), {K: P(stop at K)})."""
    from scipy.stats import norm
    ss = sigma_star(TAU_C)
    t = W.params(SIG, ss)['t']
    Kr = 1.0
    for _ in range(100):
        Kr = np.log((1 - t * t) ** 0.25 * (np.pi * Kr) ** 0.25 * rf / (s * np.sqrt(2 * DELTA + DELTA ** 2))) / np.log(1 / t)
    cdf = [0.0]
    for K in range(1, k_max + 1):
        b2 = (W.l2_error(*W.half_squeeze(K, SIG, ss), SIG, ss) * rf / s) ** 2
        m2 = b2 + 1 - 2 * K / N
        sd2 = np.sqrt(2 * (N - 2 * K) / N ** 2 + 4 * b2 / N)
        cdf.append(float(norm.cdf(((1 + DELTA) ** 2 - m2) / sd2)))
    return float(Kr), {K: cdf[K] - cdf[K - 1] for K in range(1, k_max + 1)}


def report(A=None, B=None, C=None, D=None, B2=None, A5=None, B3=None, B4=None, B5=None, A6=None, A6E=None, A6X=None, A6C=None):
    L_ = ['# Proposition D‴.58: the transit under width mismatch (`exp_width_mismatch.py`)', '',
          f"Kernels of width σ = {SIG}; a feature Aφ_σ*(· − c*), σ* = σ(1 + (τ/σ)²)^{{1/2}}; t = (σ* − σ)/(σ* + σ), λ = τ²/(σ² + σ*²) "
          "(`src/core/wide_feature.py`).", '']
    if A is not None:
        L_ += ['## A. The rules on ℝ', '',
               "Relative L² errors: the half-squeeze rule (nodes τ(σ/σ*)^{1/2}z_k, amplitudes A(σ*/σ)^{1/2}w_k exp(y_k²/(2σ(σ + σ*)))) by "
               "the trapezoid rule and by the Mehler series (their relative difference), the series' first term's share of the squared "
               "error; moment matching (nodes τz_k, amplitudes A(σ*/σ)w_k); the half-squeeze nodes with least-squares amplitudes; the "
               f"optimal K-kernel fit (variable projection; not computed below {E_FLOOR:g}, where its solve stops short of the optimum — A3 marks "
               "the optimum unconverged where the one-parameter rule beats it). Rows stop "
               f"where the series falls below {G_FLOOR:g} (the trapezoid rule's floor). Ratios per kernel e_{{K−1}}/e_K; the optimum's "
               "nodes over the half-squeeze nodes (range over the nodes).", '']
        for r_ in TAUS_A:
            rows = [x for x in A['rows'] if x['r'] == r_]
            t, lam = rows[0]['t'], rows[0]['lam']
            L_ += [f"### τ/σ = {r_:g} (1/t = {1 / t:.4f}, 1/λ = {1 / lam:.4f})", '',
                   '| K | half-squeeze | series (rel. diff.) | first term (share) | moment matching | least-squares amplitudes | optimal | '
                   'ratios: half-squeeze, optimal, moment matching | optimal nodes / half-squeeze |', '|---|---|---|---|---|---|---|---|---|']
            for i, x in enumerate(rows):
                rat = (f"{rows[i - 1]['e_hs'] / x['e_hs']:.4f}, " + (f"{rows[i - 1]['e_opt'] / x['e_opt']:.4f}"
                       if 'e_opt' in x and 'e_opt' in rows[i - 1] else '—') + f", {rows[i - 1]['e_mm'] / x['e_mm']:.4f}") if i else '—'
                nod = '—'
                if 'y_opt' in x and x['K'] > 1:
                    m = np.abs(x['y_hs']) > 1e-9
                    q = np.sort(x['y_opt'])[m] / np.sort(x['y_hs'])[m]
                    nod = f"{q.min():.4f} … {q.max():.4f}"
                L_.append(f"| {x['K']} | {x['e_hs']:.4e} | {x['e_series']:.4e} ({abs(x['e_hs'] / x['e_series'] - 1):.0e}) | "
                          f"{(x['e_first'] / x['e_series']) ** 2:.5f} | {x['e_mm']:.3e} | {x['e_ls']:.4e} | "
                          + (f"{x['e_opt']:.4e}" if 'e_opt' in x else '—') + f" | {rat} | {nod} |")
            L_.append('')
        L_ += ['### A2. The series\' tail', '',
               "D_{2K+2j,K}/D_{2K,K} (integers) for j = 1, …, 4 and K = 1, …, " + f"{K_A2}" + "; their j-th differences in K, "
               "constant if they are polynomials of degree j, against j!·binom(1/2, j)(−2)^j = −1, −1, −3, −15.", '',
               '| j | K = 1, 2, 3, 4, 5, …, ' + f"{K_A2}" + ' | j-th differences (all K) | j!·binom(1/2, j)(−2)^j |', '|---|---|---|---|']
        ks = sorted(A['ratios'])
        for j in range(1, 5):
            v = np.array([A['ratios'][K][j - 1] for K in ks], float)
            dif = np.diff(v, n=j)
            lead = factorial(j) * comb(1, 1) * float(np.prod([0.5 - i for i in range(j)]) / factorial(j)) * (-2) ** j
            L_.append(f"| {j} | {', '.join(f'{int(x)}' for x in v[:5])}, …, {int(v[-1])} | "
                      f"{', '.join(sorted({f'{d:g}' for d in dif}))} | {lead:g} |")
        if 'e_fock' in A['rows'][0]:
            L_ += ['', '### A3. Scaled rules', '',
                   "The Gauss–Hermite rule of N(0, v) for the reweighted mixing law (`scaled_rule`) at the half-squeeze variance v_hs = τ²σ/σ*, "
                   "the natural Fock scaling v_hs/cosh r and the Mercer-zero scaling τ²(1 + 2τ²/σ²)^{−1/2}: their ratios per kernel at the last "
                   "K against 1/t, (2 − λ)/λ = (1 − t + t²)/t and 1/q, q = (β² − 1)/(β² + 1), β² = (1 + 2τ²/σ²)^{1/2}; the variance of least "
                   "error (`best_scale`): v/v_hs, its error and the optimal fit's over the rule's, at K = 1 and at the last K with an optimum "
                   "(in brackets the share of the gap to the optimum that the variance recovers, (e_hs − e_best)/(e_hs − min(e_opt, e_best)); "
                   f"'optimum unconverged' where e_best < e_opt by more than a relative {UNCONV:g}); over all K with an optimum the shares "
                   + gap_range(A) + ":", '',
                   '| τ/σ | half-squeeze (1/t) | Fock ((2 − λ)/λ) | Mercer-zero (1/q) | K = 1: v/v_hs, e_best/e_hs, e_opt/e_hs | the last K: K, '
                   'v/v_hs, e_best/e_hs, e_opt/e_hs |', '|---|---|---|---|---|---|']
            for r_ in TAUS_A:
                rows = [x for x in A['rows'] if x['r'] == r_]
                t, lam = rows[0]['t'], rows[0]['lam']
                b2 = np.sqrt(1 + 2 * r_ ** 2)
                q = (b2 - 1) / (b2 + 1)
                a, b = rows[-2], rows[-1]
                opt = [x for x in rows if 'e_opt' in x]
                cell = lambda x: (f"{x['v_best']:.5f}, {x['e_best'] / x['e_hs']:.5f}, " + (f"{x['e_opt'] / x['e_hs']:.5f}" if 'e_opt' in x else '—')
                                  + (f" ({gap_share(x)})" if 'e_opt' in x else ''))
                L_.append(f"| {r_:g} | {a['e_hs'] / b['e_hs']:.3f} ({1 / t:.3f}) | {a['e_fock'] / b['e_fock']:.3f} ({(2 - lam) / lam:.3f}) | "
                          f"{a['e_mercer'] / b['e_mercer']:.3f} ({1 / q:.3f}) | {cell(rows[0])} | "
                          + (f"{opt[-1]['K']}, {cell(opt[-1])}" if opt else '—') + ' |')
        if 'narrow' in A:
            L_ += ['', '### A4. Narrower features: the continued frame', '',
                   "σ* < σ: the nodes imaginary in conjugate pairs, the K-term sum e^{−x²/(2σ²)}Σ_k b_k cos(η_kx/σ²) at the Hermite-zero "
                   "frequencies η_k = (2σ²|sinh r|)^{1/2}z_k (`narrow_rule`); its relative L² error (trapezoid) against the Mehler series at "
                   "|t| = (σ − σ*)/(σ + σ*), and the ratios per kernel against 1/|t|:", '',
                   '| σ*/σ | 1/|t| | K: error (relative difference from the series) | ratio per kernel at the last K |', '|---|---|---|---|']
            for q in NARROW_A:
                rs = [x for x in A['narrow'] if x['q'] == q]
                L_.append(f"| {q:g} | {1 / rs[0]['t']:.4f} | " + ', '.join(f"{x['K']}: {x['e']:.3e} ({abs(x['e'] / x['e_series'] - 1):.0e})"
                                                                         for x in rs) + f" | {rs[-2]['e'] / rs[-1]['e']:.4f} |")
        L_ += a5_lines(A5, A)
        L_ += a6_lines(A6, A5)
        L_ += a6x_lines(A6X, A6, A6C)
        L_ += a6c_lines(A6C, A6, A6X)
        L_ += a6e_lines(A6E)
        from scipy.special import hyp2f1
        L_ += ['', "The tail's share Σ_{n>2K}/(first term) against ₂F₁(−½, −½; 1; t²) − 1 (the share's limit if the ratios' leading "
               "coefficients govern the tail):", '', '| τ/σ | t | K = 5 | 10 | 20 | 30 | ₂F₁(−½, −½; 1; t²) − 1 |', '|---|---|---|---|---|---|---|']
        for r_ in TAUS_A:
            sh = [x['share'] for x in A['share'] if x['r'] == r_]
            t = W.params(SIG, sigma_star(r_))['t']
            L_.append(f"| {r_:g} | {t:.5f} | " + ' | '.join(f"{v:.5f}" for v in sh) + f" | {hyp2f1(-0.5, -0.5, 1, t ** 2) - 1:.5f} |")
    if B is not None:
        L_ += ['', '## B. The transit', '',
               f"The sliding path (polished; a kernel leaves only if its correlation stays inward with its own contribution removed) on "
               f"f* = {A_B}φ_σ*(· − {C_B}) sampled at {X_B.size} points of [{X_B[0]:g}, {X_B[-1]:g}], 21 kernels started at the grid "
               "−1, −0.9, …, 1, to " + f"{K_B + 1} entries, the fast flow's chunks up to 10⁵·2^{MAX_DOUBLINGS}: after each entry the relaxed "
               "block's relative rms error against the optimal fit with as many kernels on the same samples (the better of variable "
               f"projection from the block and from the half-squeeze nodes: the best local optimum found); a block is optimal within "
               f"{OPT_TOL:g}. The numerical variants: " + ', '.join(f"{n} (polish step {st:g}σ, grid shifted by {sh:g}, the polish's budget "
                                                                     f"{nf:g} evaluations" + (f", repeated up to {rp} times" if rp else '')
                                                                     + (f", its gradient test {gt:g}" if gt is not None else ', no gradient test')
                                                                     + ')' for n, st, sh, nf, rp, gt in VARIANTS_B)
               + ". 'Relaxed': the block passed the stationarity test. The first non-optimal block's polishes as status/evaluations/accepted "
               f"(status 0: the budget was spent first; 1: the gradient test), a repeated polish's end over the optimum in brackets. 'none': every block optimal, to the path's cap of {K_B} kernels "
               "(the range censored).", '',
               '| τ/σ | variant | optimal blocks (kernels) | the first non-optimal: kernels, its optimum\'s error, block / optimum, relaxed | its '
               'polishes | '
               'relaxed blocks | chunks redone (stall, failure), failed | polish rejected | exits |', '|---|---|---|---|---|---|---|---|---|']
        for r_ in TAUS_B:
            for v in VARIANTS_B:
                b = next((x for x in B if x['r'] == r_ and x['var'] == v[0]), None)
                if b is None:
                    continue
                opt = [o['k'] for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL]
                bad = next((o for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL), None)
                L_.append(f"| {r_:g} | {v[0]} | {', '.join(map(str, opt)) or '—'} | "
                          + (f"{bad['k']}, {bad['e_opt']:.2e}, {bad['e'] / bad['e_opt']:.4g}, {bad['relaxed']}" if bad else 'none') + f" | "
                          + (polish_cell(bad, float(np.sqrt(np.mean((A_B * np.exp(-(X_B - C_B) ** 2 / (2 * sigma_star(r_) ** 2))) ** 2))))
                             if bad else '—') + f" | "
                          f"{sum(o['relaxed'] for o in b['blocks'])}/{len(b['blocks'])} | {b['stiff']}, {b['redone']}, {b['failed']} | "
                          f"{b['rejected']} | {sum(e[2] == 'exit' for e in b['events'])} |")
        ok = [o['e_opt'] for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL]
        no = [o for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL]
        vac = [o['e'] for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL and o['relaxed']
               and cs_vacuous(b['r'], b['var'], o)]
        mg = [cs_margin(b['r'], b['var'], o) for b in B for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) >= OPT_TOL and o['relaxed']]
        mg = [m for m in mg if m < 1]
        L_ += ['', f"Over every width and variant: the optimal blocks' optimum errors down to {min(ok):.2e}; the non-optimal blocks' "
               + (f"from {max(o['e_opt'] for o in no):.2e} down ({sum(not o['relaxed'] for o in no)} of {len(no)} failed the stationarity "
                  "test" + (f"; {len(vac)} of those that passed it pass it at any state by Cauchy–Schwarz, the test being absolute — their "
                            f"errors up to {max(vac):.1e}" if vac else '')
                  + (f"; for the other {len(mg)} its threshold over the Cauchy–Schwarz bound — the residual's cosine with the block's "
                     f"directions it allows — is {min(mg):.1e}–{max(mg):.2f} (median {np.median(mg):.1e}; above 0.1 in "
                     f"{sum(m > 0.1 for m in mg)}): a weak test" if mg else '') + "); "
                  f"non-optimal blocks with an optimum above 10⁻⁷: {sum(o['e_opt'] > 1e-7 for o in no)}; the non-optimal blocks' errors "
                  f"{min(o['e'] for o in no):.2e} … {max(o['e'] for o in no):.2e}." if no else "none."), '']
        L_ += b2_line(B2)
        L_ += b3_lines(B3, B4, B)
        L_ += b4_lines(B4)
        L_ += b5_lines(B5, B4)
        L_ += entry_lines(B)
        for b in [x for x in B if x['var'] == VARIANTS_B[0][0]]:
            L_ += [f"### τ/σ = {b['r']:g}, {b['var']} (events: " + ', '.join(f"{j}{'+' if k == 'enter' else '−'}" for t_, j, k, c in b['events'])
                   + ")", '',
                   '| kernels | τ | block | optimal fit | block / optimal | relaxed | centres: block − optimal (σ) | block − half-squeeze (σ) | '
                   'half-squeeze on ℝ | ratio per kernel |', '|---|---|---|---|---|---|---|---|---|---|']
            for i, o in enumerate(b['blocks']):
                rat = f"{b['blocks'][i - 1]['e'] / o['e']:.3f}" if i and b['blocks'][i - 1]['k'] == o['k'] - 1 else '—'
                L_.append(f"| {o['k']} | {o['tau']:.4g} | {o['e']:.4e} | {o['e_opt']:.4e} | {o['e'] / o['e_opt']:.6f} | {o['relaxed']} | "
                          f"{o['dc'] / SIG:.1e} | {o['dc_hs'] / SIG:.3f} | {o['e_hs']:.4e} | {rat} |")
            opt = [o for o in b['blocks'] if abs(o['e'] / o['e_opt'] - 1) < OPT_TOL]
            fl = [flank(b['r'], o['c'], o['w']) for o in opt]
            if fl:
                L_ += ['', "The residual's correlation C(c) = ⟨r, φ_σ(· − c)⟩ after the optimal blocks (K = "
                       + ', '.join(str(o['k']) for o in opt) + f"): min C/max C within ±{max(abs(x[0]) for x in fl):.0e}, its peaks beyond "
                       "the outermost nodes at (left, right) " + ', '.join(f"({x[1]:.2f}, {x[2]:.2f})σ" for x in fl) + '.', '']
    if C is not None:
        N = X_B.size
        L_ += ['## C. Under noise', '',
               f"τ/σ = {TAU_C:g}, {N_DRAW_C} draws per level, the path to {K_C + 1} entries: the discrepancy rule's stop (kernels), its "
               "prediction from (i) and the predicted error there, the stopped fit's relative error (median [range]), the kernels within "
               "4σ* of the feature at the stop, the largest relative gap's kernel count (in how many draws it is the stop's).", '',
               '| s | stop: kernels (median [range]) | predicted | error at the stop: median [range] | predicted | kernels at the feature '
               '| largest gap at |', '|---|---|---|---|---|---|---|']
        L2 = ['', "At the stop (K kernels): the error over the oracle's s(2K/N)^{1/2} (the variance alone, the right support), over the "
              "bias–variance prediction (e_K²‖f*‖²_N + 2Ks²/N)^{1/2} — e_K the noiseless optimum's error with K kernels (B's base path "
              f"at τ/σ = {TAU_C:g}) — and over the discrepancy rule's bound s((1 + δ)² − 1 + 4K/N)^{{1/2}}; paths with a runaway kernel "
              "before the stop:", '', '| s | error / oracle: median [range] | / bias–variance | / bound: max | runaway |', '|---|---|---|---|---|']
        L3 = ['', "The stop as a margin (`stop_law`): the first-term formula's root K*(s) of K = log((1 − t²)^{1/4}(πK)^{1/4}‖f*‖_N/"
              "(s(2δ + δ²)^{1/2}))/log(1/t), whose ceiling is the predicted stop, and the probability of the stop at each K — the rms's "
              "expectation s(b_K² + 1 − 2K/N)^{1/2} and its scatter from the noise's chi-square and its cross term with the bias, "
              "b_K = e_K‖f*‖_N/s, e_K the half-squeeze rule's error, in the normal approximation (the rms falls along the path) —; the "
              "expected number of the draws stopping at the predicted K against the observed:", '',
              '| s | K*(s) | P(stop at K) | at the predicted K | expected | observed |', '|---|---|---|---|---|---|']
        for s in LEVELS_C:
            ds = [d for d in C if d['s'] == s]
            if not ds:
                continue
            ks, es, at, kg, same = [], [], [], [], 0
            for d in ds:
                i = next((i for i, r_ in enumerate(d['rms']) if r_ <= (1 + DELTA) * s), None)
                if i is None:
                    continue
                ks.append(d['k'][i]); es.append(d['err'][i]); at.append(d['at'][i])
                if d['gaps']:
                    kg.append(int(np.argmax(d['gaps'])) + 1)
                    same += kg[-1] == d['k'][i]
            Kp, ep = predict_stop(s, ds[0]['rf'], N)
            fmt = lambda v, f_: f"{np.median(v):{f_}} [{np.min(v):{f_}}, {np.max(v):{f_}}]" if len(v) else '—'
            L_.append(f"| {s:g} | {fmt(ks, 'g')} ({len(ks)}/{len(ds)} stopped) | {Kp} | {fmt(es, '.3e')} | {ep:.3e} | {fmt(at, 'g')} | "
                      f"{fmt(kg, 'g')}" + (f" (the stop in {same})" if kg else '') + " |")
            Kr, P = stop_law(s, ds[0]['rf'], N)
            L3.append(f"| {s:g} | {Kr:.2f} | " + ', '.join(f"{K}: {p:.3f}" for K, p in P.items() if p >= 5e-4) + f" | {P.get(Kp, 0.0):.3f} | "
                      f"{len(ds) * P.get(Kp, 0.0):.1f} of {len(ds)} | {sum(k == Kp for k in ks)} |")
            rf = ds[0]['rf']
            eopt = ({o['k']: o['e_opt'] for b in B for o in b['blocks'] if b['r'] == TAU_C and b['var'] == VARIANTS_B[0][0]}
                    if B is not None else {})
            r1 = [e / (s * np.sqrt(2 * k / N) / rf) for k, e in zip(ks, es)]
            r2 = [e / np.sqrt(eopt[k] ** 2 + 2 * k * s ** 2 / (N * rf ** 2)) for k, e in zip(ks, es) if k in eopt]
            r3 = [e / (s * np.sqrt((1 + DELTA) ** 2 - 1 + 4 * k / N) / rf) for k, e in zip(ks, es)]
            L2.append(f"| {s:g} | {fmt(r1, '.2f')} | {fmt(r2, '.3f')} | {max(r3):.3f} | {sum(d.get('runaway', False) for d in ds)} |"
                      if ks else f"| {s:g} | — | — | — | — |")
        L_ += L2 + L3
    if D is not None:
        from scipy.special import hyp2f1
        ts = sorted({r['t'] for r in D['rows']})
        Ks = sorted({r['K'] for r in D['rows']})
        by = {(r['t'], r['K']): r for r in D['rows']}
        L_ += ['', f"## D. The tail of the exact series (the sharp constant: proved by D‴.61 (iv) at every t < 1 — every σ*/σ —, {D61_STATUS}; "
                   f"before VA for t < {D61_RSTAR_V:.4f} with V, for t < {D61_RSTAR:.4f} with the branch and bound alone)", '',
               "The half-squeeze rule's squared error over the series' first term, Σ_j t^{2j}(2K)!/(2K + 2j)!·(D_{2K+2j,K}/D_{2K,K})² (the "
               "defects exact, by the remainders modulo He_K), against ₂F₁(−½, −½; 1; t²), the sum of the terms' limits binom(1/2, j)²t^{2j}; "
               "K times the gap in parentheses:", '',
               '| t | ₂F₁ | ' + ' | '.join(f'K = {K}' for K in Ks) + ' |', '|---|---|' + '---|' * len(Ks)]
        for t in ts:
            F = float(hyp2f1(-0.5, -0.5, 1.0, t * t))
            L_.append(f"| {t:.4g} | {F:.6f} | " + ' | '.join(f"{by[(t, K)]['ratio']:.6f} ({K * (by[(t, K)]['ratio'] - F):+.3f})"
                                                           for K in Ks) + ' |')
        cat = lambda j: comb(2 * (j - 1), j - 1) // j
        L_ += ['', 'The ratios D_{2K+2j,K}/D_{2K,K} are polynomials in K of degree j (the (j + 1)-th divided difference over K = 1 … j + 2: '
               + ', '.join(f"{r['next_dd']:g}" for r in D['lead']) + ') with the leading coefficients '
               + ', '.join(f"{r['lead'][0]}/{r['lead'][1]}" if r['lead'][1] != 1 else f"{r['lead'][0]}" for r in D['lead']) + ' at j = 1 … ' + str(len(D['lead'])) + ' — −C_{j−1}/2^{j−1} '
               + ('at every j' if all(r['lead'][0] * 2 ** (r['j'] - 1) == -cat(r['j']) * r['lead'][1] for r in D['lead']) else 'NOT at every j')
               + ' (C the Catalan numbers), binom(1/2, j)(−2)^j.']
    L_ += ['', 'Numbers are produced by `experiments/exp_width_mismatch.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'width_mismatch_raw.npz')
    md = os.path.join(OUT, 'width_mismatch.md')
    old = E.load(raw)
    keys = ('A', 'B', 'C', 'D', 'B2', 'A5', 'B3', 'B4', 'B5', 'A6', 'A6E', 'A6X', 'A6C')
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    idx = ({int(i) for i in sys.argv[sys.argv.index('--variants') + 1].split(',')} if '--variants' in sys.argv else None)
    rix = ({int(i) for i in sys.argv[sys.argv.index('--rules') + 1].split(',')} if '--rules' in sys.argv else None)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', lambda: run_B(res.get('B'), idx)), ('C', part_C), ('D', part_D), ('B2', lambda: part_B2(res['B'])),
                  ('A5', lambda: part_A5(res['A'])), ('B3', lambda: part_B3(res['B'])), ('B4', lambda: part_B4(res['B'], res['B3'])),
                  ('B5', lambda: part_B5(res['B'], res['B3'])), ('A6', part_A6), ('A6E', part_A6E), ('A6X', lambda: part_A6X(res.get('A6X'), rix)),
                  ('A6C', part_A6C)):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
