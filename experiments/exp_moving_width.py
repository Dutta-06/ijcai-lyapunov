"""
Proposition D‴.59: the transit with moving widths — one kernel per Gaussian feature, whatever its width (`src/core/moving_width.py`).

The paradigm's law moves each kernel's width beside its centre (σ̇ = Gβa⟨r∂_σφ⟩, ċ = Gβa⟨r∂_cφ⟩); on one-dimensional samples a kernel is
θ = (c, σ), and `rich_transit_nd`'s sliding path, loop and polish run on `moving_width.BumpsCW`. The kernels start on the grid −1, −0.9,
…, 1 with the width σ₀ = 0.12; the samples are 321 points of [−1.6, 1.6]. Against the same path with the widths fixed at σ₀
(`rich_transit_nd.BumpsND`), to MAX_ENT_FIXED entries and within MAX_RHS_FIXED right-hand-side evaluations: on a feature narrower than
the kernels the fixed widths' blocks are stiff and their relaxation runs for hours (a path stopped by the budget is reported so).

A. One feature (A = 0.8 at c* = 0.03) at the width ratios σ*/σ₀ = 1/4 … 4 (D‴.58's 1.031, 1.118, 1.414, 2.236 among them): moving widths
   against fixed widths — the kernels formed, the final (c, σ, w) against (c*, σ*, A), the residual.
B. The one-kernel reduced energy in the continuum: its local minima on a fine (c, σ) grid (one, at the exact fit), the closed-form
   gradient against finite differences.
C. Random three-feature targets (widths log-uniform in [0.05, 0.3], amplitudes ±U(0.4, 1), centres in [−0.7, 0.7] at least Δ apart):
   moving widths' exact recovery against fixed widths' kernel count and residual.
D. A narrow feature (σ = 0.06, A = 0.5) on a wide pedestal (σ = 0.3, A = 0.8) at the distances d: the kernels formed with moving widths.
E. Under noise (C's targets at Δ = 0.6, s = 10⁻³, 10⁻²): the discrepancy rule's stop (the residual's rms within (1 + δ)s), moving widths
   against fixed: the kernels and the error against the truth, against the oracle's s(3K/N)^{1/2}.
F. The real loop at small masses (ε = 10⁻⁶, 10⁻⁸, the free law on its curves, `rich_transit_nd.loop_nd` with moving widths) on C's
   first targets: its formed kernels against the path's, at TAUS_F times the path's last entry (to 150: referee round 61).
G. C's paths stopped by their budget: the stalled block's state when the budget ran out (the path as C's), its flow from the entry on
   the critical curve to s = 1e11 with the dipole law's constant (K = v²(E − E_inf), E_inf the optimal dipole fit's), and the real loop
   at the masses EPS_G on the grid TAUS_G (referee round 57): the dipole pair's half-difference v, centre of mass and net coefficient
   against the features, its three phases — the dipole's growth v³ = 3μ₀KΛ(ε)(τ − τ_d) until the next kernel forms, the landing, the
   collapse — and the landed v against (3μ₀KΛ(ε)(τ₃ − τ_d))^{1/3} ∝ ε^{−1/6} (referee round 61).
H. Two features of one sign at ±0.08 (1.33σ₀ apart), symmetric and perturbed: the kernels of the other sign (referee round 57), the
   end's kernels narrower than the samples' spacing or of negative width; H2: the σ* = 0.06 paths to the tolerance TOL_H2; H3: the
   separation scan SEP_H at σ* = 0.06, equal amplitudes; HL: the real loop on H's two equal targets (referee round 61).
"""
import itertools
import math
import os
import sys
import time
import warnings

import numpy as np
from statistics import NormalDist

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.core import rich_transit as RT
from src.core import rich_transit_nd as RN
from src.core import moving_width as MW
import exp_cohesion_leak as E
import exp_mirror_transit as X

OUT = os.environ.get('MW_OUT', os.path.join(os.path.dirname(__file__), '..', 'results', 'chaos'))
N_WORKERS = int(os.environ.get('MW_WORKERS', '2'))
XS = np.linspace(-1.6, 1.6, 321)
S0, C0 = 0.12, np.linspace(-1.0, 1.0, 21)
A_A, C_A = 0.8, 0.03
RATIOS_A = (0.25, 0.5, 0.75, float(np.sqrt(1 + 0.25 ** 2)), float(np.sqrt(1 + 0.5 ** 2)), float(np.sqrt(2.0)), float(np.sqrt(5.0)), 3.0, 4.0)
MAX_ENT_FIXED = 4                                                       # fixed widths: the path to 4 entries (a narrow feature's slow blocks)
MAX_RHS_FIXED, MAX_RHS_MOVING = 300_000, 1_000_000                      # the paths' budgets of right-hand-side evaluations (a moving-width
                                                                        # path to its exact fit spends 10³–10⁴)
SEED_C, N_C, DELTAS_C, K_C, WID_C = 590, 40, (0.4, 0.6), 3, (0.05, 0.3)
EXACT_TOL = 1e-6                                                        # exact: every feature's (c, σ, w) within it, the residual < 1e-8
D_D = (0.0, 0.1, 0.2, 0.3, 0.45, 0.6)
LEVELS_E, DELTA_E, N_E = (1e-3, 1e-2), 0.1, 20
EPS_F, N_F, TAUS_F = (1e-6, 1e-8), 4, (1.3, 6.0, 30.0, 150.0)            # F: the loop to 150 times the last entry (round 61)
EPS_G = (1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10)                         # G: the loop at six masses (the landing law: round 61)
TAUS_G = tuple(float(x) for x in np.unique(np.r_[np.linspace(10.0, 50.0, 21), np.geomspace(50.0, 6400.0, 36)].round(6)))  # 50·2^{k/5}
TAUS_G_MD = (400.0, 1600.0, 6400.0)                                     # G: the loop's states the md tabulates
S_G = (1e3, 1e5, 1e7, 1e9, 3.2e10, 1e11)                                # G: the block's flow to s = 1e11
LAND_TOL = 0.01                                                         # G: landed — the pair's centre of mass within 1% of a feature's
                                                                        # width of it and its net coefficient within 1% of its amplitude
SEP_H, TOL_H2 = (0.05, 0.1, 0.12, 0.15, 0.2), 1e-13                     # H3: σ* = 0.06, equal, at ±SEP_H; H2: the paths to this tolerance
EPS_HL, TAUS_HL = (1e-6, 1e-8), (1600.0, 6400.0, 2.5e4, 1e5)            # HL: the loop on H's targets (referee round 68: to 1e5)
HL_NAMES = ('σ* 0.04, equal', 'σ* 0.06, equal', 'σ* 0.06, A₂ 0.99', 'σ* 0.06, A₂ 0.95', 'σ* 0.06, A₂ 0.9', 'σ* 0.06, c₂ +0.01',
            'σ* 0.06, σ₂ 0.065')                                        # HL: the symmetric targets and the asymmetric ones
HP_NAMES = ('σ* 0.04, equal', 'σ* 0.06, equal', 'σ* 0.06, c₂ +0.01', 'σ* 0.06, A₂ 0.9')   # HP: H's path without the polish
TAUS_FINE = tuple(float(x) for x in np.round(np.arange(10.0, 140.05, 0.1), 6))  # G (referee round 68): the formations on a fine grid
EPS_GX = (1e-11, 1e-12)                                                 # G: two smaller masses (the fine phases only)
P_LAND = 0.28                                                           # G: (ii)(a)'s parameter p = |w_a − w_b||Δc|/(|W|σ) at the landing
SEQ_TOL = 0.5                                                           # G: the loop's intermediate slope over the path's law beyond 1 ± it:
                                                                        # another sequence (referee round 80: 27 at ε ≥ 10⁻⁹, 18 at 10⁻⁷)
CLOSURE_EXACT = 1e-12                                                   # GT: a closure fit whose residual's rms is below it is exact
                                                                        # (referee round 80: target 20's K′ and τ₄ were roundoff)
Q_GL, EPS_GL = (20, 9), 1e-10                                           # GL: two dipole targets' collapse (9: referee round 74), to TAUS_GL
TAUS_GL = tuple(float(x) for x in np.geomspace(400.0, 1e7, 29))
RATIOS_A2, WMUL_A2 = (0.25, 1.0, 4.0), (2.5, 3.0, 5.0, 10.0, 20.0)      # A2: one kernel from (c*, σ*) with w = WMUL_A2·A (rounds 74, 80)
A2_OFF = 1e-9                                                           # A2 (referee round 87, D3ch-1): the start's centre c* + A2_OFF·σ*, off
                                                                        # the sample c* — from c* exactly LSODA holds 3A at 0.25 on the saddle
S_A2 = tuple(float(x) for x in np.geomspace(1.0, 1e13, 27))             # A2: the block's flow's outputs (referee round 80: to 10¹³)
PLATEAU_A2, ESC_A2 = 1e-3, 0.5                                          # A2: on the plateau, E within it of E₂ (`a2_plateau`;
                                                                        # referee round 83); escaped once σ/σ* is back above ESC_A2
S_A2E = tuple(10.0 ** (k / 4) for k in range(-16, 5))                    # A2E (referee round 83): A2's 5A at 0.25 to s = 10
A2E_ROUND = 1e-14                                                       # A2E (referee round 97, D3cr-9): an offset below it is roundoff
A2E_MORE = ((1e13, 0.25, 3.0, 0.0, 10.0),)                                         # A2E (referee round 87, D3ch-1): 3A at 0.25 from c* exactly
SPANS_A2E, S_A2E_STEPS = (10.0, 1e13), (1e-3, 1e-2)                                                # A2E: the spans (A2's is 10¹³)
WMUL_A2P, RTOL_A2P, K_A2P, SEED_A2P = (3.0, 5.0), (1e-10, 1e-12), 100, (131, 132)   # A2P (referee rounds 125, 130: D3dt-1, D3dy-1): the
                                                                        # symmetric samples' exit over the code's order (one draw) and K_A2P
                                                                        # random orders per seed (the distribution; fresh seeds, round 125's
                                                                        # being its referee's own 125)
RTOL_A2H, K_A2H, FRAC_A2H = 1e-10, 100, (0.25, 0.0625)                  # A2H (referee round 130, D3dy-2): A2P's first K_A2H orders of
                                                                        # SEED_A2P[0] at RTOL_A2H with DOP853's step at most f·h₀ (f in
                                                                        # FRAC_A2H, both arithmetics), h₀ the double's mean step over A2P's
                                                                        # draws (`a2p_h0`)
SEED_A2R, ATOL_A2R = (133, 134, 135, 136), 1e-14                         # A2R (referee round 133, D3eb-1, -2): A2P's cell farthest from
CELL_A2R, CAP_A2R = (3.0, 1e-10), (5.0, 0.0625)                         # 11 ln 2 (CELL_A2R) on the further seeds SEED_A2R, and on all its
                                                                        # seeds with the double at the long double's atol ATOL_A2R·σ*; A2H's
                                                                        # cap beyond 3 sem (CAP_A2R: w, f) on SEED_A2P[1]'s first K_A2H orders
EXP_A2P, GM_A2P = (0.75, 1.25), (0.05, 20.0)                            # A2P's rule (referee rounds 133, 137: D3eb-1, D3ef-2): every cell's
                                                                        # exponent Λ/(11 ln 2) with its 99% interval inside EXP_A2P on all its
                                                                        # draws (A2P with A2R: `a2p_gate`), every geometric mean of the seeds
                                                                        # within GM_A2P (units of uσ*). EXP_A2P halfway to the neighbouring
                                                                        # half-integer exponents ½ and 3/2 (at γ = 1 a prefactor ratio within
                                                                        # e^{±0.25·11 ln 2}); round 133's band, set after the data, EXP_A2P_133
EXP_A2P_133 = (0.9, 1.1)
DT_A2T = (('float32', np.float32), ('double', np.float64), ('long double', np.longdouble))   # A2T (referee round 137, D3ef-1, -3): three
                                                                        # arithmetics of the same right-hand side (`_DTBumps`)
S_A2T = ((3.0, (0.02, 0.04, 0.06)), (5.0, (0.001, 0.002, 0.004)))       # A2T's scan: the line's states (w, s on the long-double line flow)
NT_A2T, ULP_A2T, SEED_A2T = 801, 4, 141                                 # the scan of ċ(δ) on NT_A2T points of |δ| ≤ ULP_A2T ulps of the ±1
                                                                        # samples' coordinate, in the code's order and SEED_A2T's first orders
RK_A2T = ((3.0, 1.0, 30, 142), (3.0, 0.5, 30, 142), (3.0, 0.25, 30, 142), (5.0, 0.25, 200, 143), (5.0, 0.0625, 100, 144))
SMAX_A2T = 0.5                                                          # A2T's fixed-step RK4 arms (w, f, orders, seed): step f·h₀, to s = SMAX_A2T
HLL = (('σ* 0.06, equal', 1e-6), ('σ* 0.06, equal', 1e-8), ('σ* 0.06, A₂ 0.99', 1e-8), ('σ* 0.06, A₂ 0.95', 1e-6))
TAUS_HLL = (1e5, 1e6, 1e7, 1e8, 1e9)                                    # HLL: HL's cells to 10⁹ (referee rounds 74, 80)
VIS = 1e-10                                                             # visible on the samples: the largest φ there at least VIS
H_TARGETS = (('σ* 0.04, equal', (-0.08, 0.08), (0.04, 0.04), (1.0, 1.0)), ('σ* 0.06, equal', (-0.08, 0.08), (0.06, 0.06), (1.0, 1.0)),
             ('σ* 0.06, A₂ 0.99', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.99)), ('σ* 0.06, A₂ 0.95', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.95)),
             ('σ* 0.06, A₂ 0.9', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.9)), ('σ* 0.06, c₂ +0.01', (-0.08, 0.09), (0.06, 0.06), (1.0, 1.0)),
             ('σ* 0.06, σ₂ 0.065', (-0.08, 0.08), (0.06, 0.065), (1.0, 1.0)), ('σ* 0.04, A₂ 0.6', (-0.08, 0.08), (0.04, 0.04), (1.0, 0.6)),
             ('σ* 0.06, A₂ 0.6', (-0.08, 0.08), (0.06, 0.06), (1.0, 0.6)))
PAIR_TOL = 0.01                                                         # a coinciding pair: |Δc|/σ + |Δlog σ| below it (a same-sign
                                                                        # pair a little apart fits one Gaussian to O(Δc⁴))


def law():
    return X.law()


def theta0():
    return np.stack([C0, np.full(C0.size, S0)], -1)


def feature_sum(cs, ss, As, x=XS):
    return np.exp(-(x[:, None] - np.asarray(cs)[None, :]) ** 2 / (2 * np.asarray(ss)[None, :] ** 2)) @ np.asarray(As)


def run_paths(f, max_entries=None, max_events=40, fixed=True):
    """The sliding path with moving widths (to max_entries) and with fixed widths (to MAX_ENT_FIXED entries) on the samples' target f."""
    L = law()
    Rm = RN.sliding_iss_nd(L, MW.BumpsCW(XS, S0), f, theta0(), np.ones(C0.size), max_events=max_events, max_entries=max_entries,
                           exit_kkt=True, max_rhs=MAX_RHS_MOVING)
    Rf = (RN.sliding_iss_nd(L, RN.BumpsND(XS[:, None], S0), f, C0[:, None], np.ones(C0.size), max_events=max_events,
                            max_entries=MAX_ENT_FIXED, exit_kkt=True, max_rhs=MAX_RHS_FIXED) if fixed else None)
    return Rm, Rf


def match(R, cs, ss, As):
    """The active kernels against the features: the count, and for each feature the nearest active kernel — in (c, log σ) with moving
    widths, so that features at one centre are told apart by their widths — with its |Δc|, |Δσ| and |Δw|; exact if as many kernels as
    features, each within EXACT_TOL, the residual below 1e-8."""
    act = np.flatnonzero(R['act'])
    th, w = np.atleast_2d(R['c'][act]), R['w'][act]
    errs = []
    for c_, s_, a_ in zip(cs, ss, As):
        if not act.size:
            errs.append((np.inf, np.inf, np.inf))
            continue
        dist = np.abs(th[:, 0] - c_) / S0 + (np.abs(np.log(np.abs(th[:, 1]) / s_)) if th.shape[1] > 1 else 0.0)
        j = int(np.argmin(dist))
        ds = abs(th[j, 1] - s_) if th.shape[1] > 1 else abs(S0 - s_)
        errs.append((abs(th[j, 0] - c_), ds, abs(w[j] - a_)))
    errs = np.array(errs)
    rms = R['hist'][-1][3]
    st = R.get('stalled')
    return dict(n=int(act.size), err=errs.max(0).tolist(), rms=float(rms),
                exact=bool(act.size == len(cs) and errs.max() < EXACT_TOL and rms < 1e-8),
                entries=sum(e[2] == 'enter' for e in R['events']), exits=sum(e[2] == 'exit' for e in R['events']),
                theta=th.tolist(), w=w.tolist(), oob=bool(R.get('out_of_budget', False)), nfev=int(R.get('nfev', 0)),
                stalled=None if st is None else dict(theta=np.atleast_2d(st['c']).tolist(), w=np.asarray(st['w']).tolist()))


def pairs_of(theta, w):
    """The coinciding pairs among active kernels: |Δc|/σ + |Δlog σ| < PAIR_TOL (σ their mean width), with their coefficients."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    out = []
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            sm = 0.5 * (abs(th[i, 1]) + abs(th[j, 1]))
            if abs(th[i, 0] - th[j, 0]) / sm + abs(np.log(abs(th[i, 1]) / abs(th[j, 1]))) < PAIR_TOL:
                out.append((float(w[i]), float(w[j])))
    return out


def closest_pair(theta, w, opposite=False):
    """The closest pair of kernels (|Δc|/σ + |Δlog σ|, σ their mean width; of opposite signs if asked): (distance, i, j) or None."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    best = None
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            if opposite and w[i] * w[j] >= 0:
                continue
            sm = 0.5 * (abs(th[i, 1]) + abs(th[j, 1]))
            d = abs(th[i, 0] - th[j, 0]) / sm + abs(np.log(abs(th[i, 1]) / abs(th[j, 1])))
            if best is None or d < best[0]:
                best = (float(d), i, j)
    return best


def width_split(theta, w, i, j):
    """A pair's width split against the third moment's (referee round 57): with Δc = c_j − c_i, W = w_i + w_j and σ their mean width the
    split σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ) cancels the pair's third moment, so the pair fits one Gaussian to O(Δc⁴). Returns (Δc,
    the observed split, the predicted)."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    if th[i, 0] > th[j, 0]:
        i, j = j, i
    dc = th[j, 0] - th[i, 0]
    sm = 0.5 * (abs(th[i, 1]) + abs(th[j, 1]))
    return float(dc), float(abs(th[j, 1]) - abs(th[i, 1])), float((w[j] - w[i]) * dc ** 2 / (6 * (w[i] + w[j]) * sm))


def moving_state(r):
    """C's moving-width path's formed kernels and residual — at its stalled block (G) where its budget stopped it (referee round 68: the
    table printed the block before it)."""
    m = r['moving']
    st = m.get('stalled')
    if st is None:
        return n_formed(m, r['As']), m['rms']
    w = np.asarray(st['w'], float)
    f = feature_sum(r['cs'], r['ss'], r['As'])
    rms = float(np.sqrt(np.mean((f - MW.BumpsCW(XS, S0).synth(np.atleast_2d(np.asarray(st['theta'], float)), w)) ** 2)))
    return int(np.sum(np.abs(w) > 1e-2 * min(abs(a) for a in r['As']))), rms


def n_formed(m, As):
    """The formed kernels of a path's end: |w| above 1% of the smallest amplitude (as F)."""
    return int((np.abs(np.asarray(m['w'], float)) > 1e-2 * np.abs(np.asarray(As)).min()).sum())


def outcome(m, As, cs=None, ss=None):
    """A moving-width path's outcome: 'exact' (as many formed kernels — |w| above 1% of the smallest amplitude — as features, within
    EXACT_TOL, the residual below 1e-8; idle entrants below 1% allowed), 'exact, null cluster' (each feature's nearest kernel within
    EXACT_TOL and the residual below 1e-8, with more formed kernels: the others sum to nothing on the samples — referee round 61),
    'coinciding pair' (the residual below 1e-8 with two formed kernels at one (c, σ), PAIR_TOL, whose net coefficient and centre of mass
    fit a feature within LAND_TOL — of its amplitude, of its width — when the features (cs, ss) are given: a split or a cancelling pair),
    'coinciding pair off the features' (such a pair fitting no feature; round 103: the null cluster is `null_cluster`'s), 'cancellation' (the residual below 1e-8 with a coefficient beyond 1.5 times the
    largest amplitude or a width below the samples' spacing: several kernels cancelling), 'dipole run-off' (stopped by the budget with an
    opposite-sign pair in the stalled block, |w| beyond 10 times the largest amplitude and 0.05 apart: G), 'stopped by the budget' (another
    stall), or 'not fitted' (referee round 57)."""
    w = np.asarray(m['w'], float)
    As_ = np.asarray(As, float)
    A_ = np.abs(As_)
    formed = np.abs(w) > 1e-2 * A_.min()
    if m.get('oob'):
        st = m.get('stalled')
        if st is not None:
            ws = np.asarray(st['w'], float)
            cp = closest_pair(st['theta'], ws, opposite=True)
            if cp is not None and np.abs(ws).max() > 10 * A_.max() and cp[0] < 0.05:
                return 'dipole run-off'
        return 'stopped by the budget'
    th = np.atleast_2d(m['theta'])
    if m['rms'] < 1e-8 and max(m['err']) < EXACT_TOL:
        if formed.sum() == len(As):
            return 'exact'
        if cs is None or null_cluster(th[formed], w[formed], cs, ss, As, m['rms']):   # referee round 100: one predicate
            return 'exact, null cluster'
    if m['rms'] < 1e-8:
        thf, wf = th[formed], w[formed]
        pr = [p for p in pair_index(thf, wf)]
        if pr:
            fit = cs is None or any(pair_fits(thf, wf, i, j, cs, ss, As_) for i, j in pr)
            return 'coinciding pair' if fit else 'coinciding pair off the features'   # round 103 (D3cx-4): not a null cluster
    if m['rms'] < 1e-8 and formed.any() and (np.abs(w[formed]).max() > 1.5 * A_.max() or np.abs(th[formed, 1]).min() < XS[1] - XS[0]):
        return 'cancellation'
    return 'not fitted'


def pair_index(theta, w):
    """The coinciding pairs among the kernels (`pairs_of`'s rule) by their indices."""
    th = np.atleast_2d(np.asarray(theta, float))
    out = []
    for i in range(len(w)):
        for j in range(i + 1, len(w)):
            sm = 0.5 * (abs(th[i, 1]) + abs(th[j, 1]))
            if abs(th[i, 0] - th[j, 0]) / sm + abs(np.log(abs(th[i, 1]) / abs(th[j, 1]))) < PAIR_TOL:
                out.append((i, j))
    return out


def pair_fits(theta, w, i, j, cs, ss, As):
    """Whether the pair (i, j) fits a feature: its net coefficient W within LAND_TOL of some feature's amplitude and its centre of mass
    within LAND_TOL of that feature's width from its centre (any feature: features may share a centre, as D's at d = 0)."""
    th = np.atleast_2d(np.asarray(theta, float))
    W = w[i] + w[j]
    if W == 0:
        return False
    cm = (w[i] * th[i, 0] + w[j] * th[j, 0]) / W
    return any(abs(cm - c_) <= LAND_TOL * s_ and abs(W - a_) <= LAND_TOL * abs(a_) for c_, s_, a_ in zip(cs, ss, As))


def a_job(ratio):
    t0 = time.time()
    ss = ratio * S0
    f = feature_sum([C_A], [ss], [A_A])
    Rm, Rf = run_paths(f)
    return dict(ratio=ratio, moving=match(Rm, [C_A], [ss], [A_A]), fixed=match(Rf, [C_A], [ss], [A_A]),
                tau=float(Rm['events'][0][0]) if Rm['events'] else np.nan, secs=time.time() - t0)


def part_A():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a_job, RATIOS_A, chunksize=1)


def a2_job(args):
    """A2 (referee round 74): one kernel's block on one feature (A's), started at (c*, σ*) with the coefficient w = wmul·A — E(start) =
    (1 − wmul)²·½⟨f*²⟩, below ½⟨f*²⟩ (i)'s condition) only for 0 < wmul < 2 —: its flow on the critical curve (LSODA, the exact
    Jacobian) to S_A2[-1], at the outputs S_A2 its w/A, σ/σ*, (c − c*)/σ* and E over ½⟨f*²⟩. Referee round 87 (D3ch-1): the centre starts
    at c* + A2_OFF·σ*, off the sample c* — the symmetric state on it is a saddle in c that an implicit solver's damping can hold."""
    ratio, wmul = args
    from scipy.integrate import solve_ivp
    t0 = time.time()
    L, K = law(), MW.BumpsCW(XS, S0)
    ss = ratio * S0
    f = feature_sum([C_A], [ss], [A_A])
    w0 = np.array([wmul * A_A])
    z0 = np.r_[np.sign(w0) * RT._alpha_of(L, np.abs(w0), RT.coef_crit), C_A + A2_OFF * ss, ss]
    sol = solve_ivp(RN.critical_rhs_nd(L, K, f), (0.0, S_A2[-1]), z0, method='LSODA', jac=RN.critical_jac_nd(L, K, f), rtol=1e-10,
                    atol=1e-14, t_eval=S_A2)
    e0 = float(np.mean(f ** 2) / 2)
    rows = []
    for i in range(sol.y.shape[1]):
        w, th = RT.critical_w(L, sol.y[:1, i]), sol.y[1:, i].reshape(1, 2)
        rows.append(dict(s=float(sol.t[i]), w=float(w[0]) / A_A, sig=float(th[0, 1]) / ss, dc=float(th[0, 0] - C_A) / ss,
                         E=float(np.mean((f - K.synth(th, w)) ** 2) / 2) / e0))
    return dict(ratio=ratio, wmul=wmul, E0=(1 - wmul) ** 2, rows=rows, status=int(sol.status), off=A2_OFF, secs=time.time() - t0)


def part_A2():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(a2_job, [(r_, m_) for r_ in RATIOS_A2 for m_ in WMUL_A2], chunksize=1)


def a2e_job(span, ratio=0.25, wmul=5.0, off=0.0, steps_to=None):
    """A2E (referee round 83): A2's start (`a2_job`: LSODA, the exact Jacobian, its tolerances) over the span (0, span), the outputs at
    S_A2E (s ≤ 10; a terminal event past them stops a longer span without changing its steps): w/A, σ/σ*, (c − c*)/σ*, E over
    ½⟨f*²⟩ and the flow's ∂ċ/∂c (its exact Jacobian); and LSODA's own steps on (S_A2E_STEPS[0], S_A2E_STEPS[1]] (the same solver,
    stepped directly): s, the step h, the method (ODEPACK's MUSED: 1 Adams, 2 BDF) and (c − c*)/σ*. c* = C_A is a sample: the collapse
    passes the symmetric state on it (c = c*, J_N's one-sample branch), a saddle in c (∂ċ/∂c > 0), which it leaves once the seed — the
    data's drift on A's samples (referee round 109) — has grown, unless the solver damps the offset. Referee round 87 (D3ch-1, D3ch-2): `off` starts the centre at c* + off·σ*
    (A2's start since round 87); each step also records its order (ODEPACK's NQU), BDF-k damping a positive real eigenvalue only for
    hλ beyond Σ_{j≤k} 2^j/j (2, 4, 6.67, … for k = 1, 2, 3). Referee round 100 (D3cu-1): steps_to — the steps recorded to this s
    instead of S_A2E_STEPS[1] (the exact-start run's BDF steps to s = 10)."""
    from scipy.integrate import solve_ivp
    L, K = law(), MW.BumpsCW(XS, S0)
    ss = ratio * S0
    f = feature_sum([C_A], [ss], [A_A])
    w0 = np.array([wmul * A_A])
    z0 = np.r_[np.sign(w0) * RT._alpha_of(L, np.abs(w0), RT.coef_crit), C_A + off * ss, ss]
    stop = lambda s, z: s - 1.05 * S_A2E[-1]
    stop.terminal = True
    jac = RN.critical_jac_nd(L, K, f)
    sol = solve_ivp(RN.critical_rhs_nd(L, K, f), (0.0, span), z0, method='LSODA', jac=jac, rtol=1e-10, atol=1e-14, t_eval=S_A2E,
                    events=stop)
    e0 = float(np.mean(f ** 2) / 2)
    rows = []
    for i in range(sol.y.shape[1]):
        w, th = RT.critical_w(L, sol.y[:1, i]), sol.y[1:, i].reshape(1, 2)
        rows.append(dict(s=float(sol.t[i]), w=float(w[0]) / A_A, sig=float(th[0, 1]) / ss, dc=float(th[0, 0] - C_A) / ss,
                         E=float(np.mean((f - K.synth(th, w)) ** 2) / 2) / e0, jcc=float(np.asarray(jac(sol.t[i], sol.y[:, i]))[1, 1])))
    from scipy.integrate import LSODA
    so = LSODA(RN.critical_rhs_nd(L, K, f), 0.0, z0, span, rtol=1e-10, atol=1e-14, jac=jac)
    steps, tp = [], 0.0
    while so.status == 'running' and so.t < (S_A2E_STEPS[1] if steps_to is None else steps_to):
        so.step()
        if so.t > S_A2E_STEPS[0]:
            iw = so._lsoda_solver._integrator.iwork
            steps.append((float(so.t), float(so.t - tp), int(iw[18]), float((so.y[1] - C_A) / ss), int(iw[13])))
        tp = so.t
    return dict(span=float(span), ratio=ratio, wmul=wmul, off=off, rows=rows, status=int(sol.status), steps=steps,
                steps_to=S_A2E_STEPS[1] if steps_to is None else float(steps_to))


def part_A2E():
    return [a2e_job(s_) for s_ in SPANS_A2E] + [a2e_job(*x) for x in A2E_MORE]


def bdf_bound(k):
    """The largest hλ > 0 at which BDF-k does not damp a positive real eigenvalue λ: its root −1 at hλ = Σ_{j≤k} 2^j/j (referee round 87)."""
    return sum(2 ** j / j for j in range(1, k + 1))


def a2_invariance(ratio, wmuls=(3.0, 5.0), sigs=(1.0, 0.5, 0.2, 0.1)):
    """Referee round 97 (D3cr-1): is the line c = c* invariant for A2's problem? The samples without a mirror image about c* (their
    number and largest f*/A), and the flow's |ċ| at c = c* (w = wmuls·A, σ/σ* = sigs) on the samples and on samples exactly symmetric
    about the centre (the problem translated to c* = 0 on a symmetric grid of the same spacing and size): roundoff on the latter, and on
    the former unless the unpaired samples carry f*."""
    L, d, ss = law(), XS[1] - XS[0], ratio * S0
    un = [x for x in XS if not np.any(np.abs(XS - (2 * C_A - x)) < d / 4)]
    fu = max([float(np.exp(-(x - C_A) ** 2 / (2 * ss * ss))) for x in un] or [0.0])
    lg = max([-(x - C_A) ** 2 / (2 * ss * ss) / np.log(10) for x in un] or [-np.inf])   # round 100: log₁₀ of that f*/A, exactly
    out = []
    for xs, c0 in ((XS, C_A), (d * np.arange(-(XS.size // 2), XS.size // 2 + 1), 0.0)):
        rhs = RN.critical_rhs_nd(L, MW.BumpsCW(xs, S0), A_A * np.exp(-(xs - c0) ** 2 / (2 * ss * ss)))
        out.append(max(abs(float(rhs(0.0, np.r_[RT._alpha_of(L, np.array([m_ * A_A]), RT.coef_crit), c0, g_ * ss])[1]))
                       for m_ in wmuls for g_ in sigs))
    return len(un), fu, out[0], out[1], float(lg)


def a2_state(ratio, w, sig, dc=0.0):
    """A2E's stored state (w/A, σ/σ*, (c − c*)/σ*) as the block's z = (α, c, σ), with the law and the samples' kernel."""
    L, K = law(), MW.BumpsCW(XS, S0)
    ss = ratio * S0
    al = np.sign(w) * RT._alpha_of(L, np.array([abs(w) * A_A]), RT.coef_crit)
    return L, K, ss, np.r_[al, C_A + dc * ss, sig * ss]


def a2_fixed_point(ratio, w, sig, dc):
    """Referee round 100 (D3cu-1): at a held state, the line's drift d = ċ at c = c* (in double: the samples' roundoff) and the
    unstable rate λ = ∂ċ/∂c there; an implicit step at hλ ≫ 1 holds the offset at its fixed point −d/λ. Returns (d/σ*, λ, −d/(λσ*),
    the held offset (c − c*)/σ*)."""
    L, K, ss, z = a2_state(ratio, w, sig)
    f = feature_sum([C_A], [ss], [A_A])
    dd = float(RN.critical_rhs_nd(L, K, f)(0.0, z)[1])
    lam = float(np.asarray(RN.critical_jac_nd(L, K, f)(0.0, z))[1, 1])
    return dd / ss, lam, -dd / lam / ss, dc


def a2_lemma(ratio, w, sig):
    """Lemma D‴.59a's split (referee round 100, D3cu-2): on samples symmetric about c* the line c = c* is invariant and ∂ċ/∂c =
    Gm·w·⟨r∂²_cφ − w(∂_cφ)²⟩ there, ∂_cφ = φu/σ², ∂²_cφ = φ(u²/σ⁴ − 1/σ²), u = x − c*; per sample offset k from c* (in samples) its
    term, the centre's (Gm w/(Nσ²))(w − f*(c*)). Returns ({'0': centre, '1': ±1, '2': ±2, 'rest'}, their sum, the exact Jacobian's
    ∂ċ/∂c, the sample spacing over σ)."""
    L, K, ss, z = a2_state(ratio, w, sig)
    f = feature_sum([C_A], [ss], [A_A])
    wv, s_ = float(RT.critical_w(L, z[:1])[0]), float(z[2])
    u = XS - C_A
    P = np.exp(-u ** 2 / (2 * s_ * s_))
    r = f - wv * P
    t = L.G * L.m * wv / XS.size * (r * P * (u * u / s_ ** 4 - 1 / s_ ** 2) - wv * (P * u / s_ ** 2) ** 2)
    k = np.rint(u / (XS[1] - XS[0])).astype(int)
    parts = {'0': float(t[k == 0].sum()), '1': float(t[np.abs(k) == 1].sum()), '2': float(t[np.abs(k) == 2].sum()),
             'rest': float(t[np.abs(k) > 2].sum())}
    J = float(np.asarray(RN.critical_jac_nd(L, K, f)(0.0, z))[1, 1])
    return parts, float(t.sum()), J, float((XS[1] - XS[0]) / s_)


def a2_ulps(ratio, w, sig, kmax=40):
    """Referee round 103 (D3cx-1): in double, where the computed problem's own symmetric state lies — the centre sample nearest c* and
    the zero of ċ along c = c* + k·ulp(c*) at the state (w/A, σ/σ*), both in ulps of c* (k interpolated between the sign change).
    Returns (the centre sample's offset in ulps, the zero of ċ in ulps or None, one ulp over σ*)."""
    L, K, ss, z = a2_state(ratio, w, sig)
    f = feature_sum([C_A], [ss], [A_A])
    rhs = RN.critical_rhs_nd(L, K, f)
    u = float(np.spacing(C_A))
    ic = int(np.argmin(np.abs(XS - C_A)))
    v = []
    for k in range(kmax + 1):
        zz = z.copy()
        zz[1] = C_A + k * u
        v.append(float(rhs(0.0, zz)[1]))
    kz = next((k + v[k] / (v[k] - v[k + 1]) for k in range(kmax) if v[k] <= 0 < v[k + 1] or v[k] >= 0 > v[k + 1]), None)
    return float((XS[ic] - C_A) / u), kz, u / ss


def a2_ulps_exact(ratio, w, sig, kmax=20, prec=50):
    """Referee round 106 (D3da-5): `a2_ulps`' zero of ċ along c = c* + k·ulp(c*), in prec-digit decimal arithmetic on the double data
    (XS, the stored f*, the state's w and σ as doubles): ċ ∝ Σ_k (f_k − wφ_k)φ_k(x_k − c) (Gm w/σ² > 0 and the mean's 1/N aside). The
    zero in ulps (interpolated between the sign change), or None."""
    from decimal import Decimal, localcontext
    L, K, ss, z = a2_state(ratio, w, sig)
    f = feature_sum([C_A], [ss], [A_A])
    wv = float(RN.critical_w(L, np.atleast_1d(z[0]), 0.0)[0])
    u = float(np.spacing(C_A))
    with localcontext() as ctx:
        ctx.prec = prec
        xs, fs = [Decimal(float(x)) for x in XS], [Decimal(float(v)) for v in f]
        W, S2 = Decimal(wv), 2 * Decimal(float(z[2])) ** 2
        def S(k):
            c = Decimal(float(C_A)) + k * Decimal(u)
            tot = Decimal(0)
            for x, fk in zip(xs, fs):
                ph = (-(x - c) ** 2 / S2).exp()
                tot += (fk - W * ph) * ph * (x - c)
            return tot
        v = [S(k) for k in range(kmax + 1)]
        k0 = next((k for k in range(kmax) if (v[k] <= 0 < v[k + 1]) or (v[k] >= 0 > v[k + 1])), None)
        return None if k0 is None else float(k0 + v[k0] / (v[k0] - v[k0 + 1]))


class _FsumBumps(MW.BumpsCW):
    """`BumpsCW` with the samples' sums in `corr` exactly rounded (`math.fsum`; referee round 109, D3dd-1): the same data, another
    arithmetic for the sums."""

    def corr(self, th, r, second=False):
        assert not second
        u, s, P = self._parts(th)
        rP = r[:, None] * P
        a, b = u / s ** 2, u ** 2 / s ** 3
        g = np.array([math.fsum(rP[:, j]) for j in range(P.shape[1])]) / self.N
        h = np.array([[math.fsum((rP * a)[:, j]), math.fsum((rP * b)[:, j])] for j in range(P.shape[1])]) / self.N
        return g, h


class _LDBumps:
    """`BumpsCW` for one kernel in long double with the centre's offset δ = c − c_ref as the state (referee round 112, D3dg-1): u = (x −
    c_ref) − δ, the Gaussians, products and sums in np.longdouble (x87's 64-bit mantissa here), the correlations returned in double. The
    data are the double samples; only the arithmetic changes — and the state no longer rounds c near c_ref to c's ulp."""

    def __init__(self, xs, sig, cref):
        LD = np.longdouble
        self.xr = np.asarray(xs, float).astype(LD) - LD(float(cref))
        self.N, self.n, self.sig = self.xr.size, 2, float(sig)

    def _parts(self, th):
        LD = np.longdouble
        th = np.atleast_2d(np.asarray(th, float))
        d, s = th[:, 0].astype(LD), th[:, 1].astype(LD)
        u = self.xr[:, None] - d[None, :]
        return u, s, np.exp(-u ** 2 / (2 * s ** 2))

    def synth(self, th, w):
        return self._parts(th)[2] @ np.asarray(w, float).astype(np.longdouble)

    def corr(self, th, r, second=False):
        assert not second
        u, s, P = self._parts(th)
        rP = np.asarray(r)[:, None] * P
        a, b = u / s ** 2, u ** 2 / s ** 3
        g = rP.sum(0) / self.N
        h = np.stack([(rP * a).sum(0), (rP * b).sum(0)], -1) / self.N
        return g.astype(float), h.astype(float)


def a2_exit_ld(ratio, wmul, sym=False, level=1e-2, smax=1.0, rtol=1e-12, off=0.0, dense=False, perm=None, line=False, max_step=np.inf,
               steps=False):
    """Referee round 112 (D3dg-1): `a2_exit` in long double with the offset δ = c − c* as the state (`_LDBumps`; DOP853, atol 10⁻¹⁴·σ*):
    the signed exit (s, ±1) from c* + off (off in units of σ*; None by smax) on A's samples or (sym) on exactly symmetric ones; with
    dense, also the solution and the right-hand side (for `a2_seed`). Round 125 (D3dt-1): perm — the samples (and f*) in that order, the
    same data summed in another order; line — ċ set to 0, the flow on the line c = c* itself (no exit). Round 130 (D3dy-2): DOP853's
    step at most max_step, and (steps) the exit with DOP853's number of steps to it, ((s, ±1), n)."""
    from scipy.integrate import solve_ivp
    LD = np.longdouble
    L, ss, d = law(), ratio * S0, XS[1] - XS[0]
    xs, c0 = (d * np.arange(-(XS.size // 2), XS.size // 2 + 1), 0.0) if sym else (XS, C_A)
    xs = xs if perm is None else xs[perm]
    K = _LDBumps(xs, S0, c0)
    f = (A_A * np.exp(-(xs - c0) ** 2 / (2 * ss * ss))).astype(LD)   # round 116 (D3dk-1): the stored double data, not re-evaluated
    z0 = np.r_[RT._alpha_of(L, np.array([wmul * A_A]), RT.coef_crit), off * ss, ss]
    ev = lambda t, y: abs(y[1]) - level * ss
    ev.terminal = True
    rhs = RN.critical_rhs_nd(L, K, f)
    fl = (lambda t, z: rhs(t, z) * (1.0, 0.0, 1.0)) if line else rhs
    sol = solve_ivp(fl, (0.0, smax), z0, method='DOP853', rtol=rtol, atol=1e-14 * ss, events=ev, dense_output=dense, max_step=max_step)
    t = (None, 0) if not sol.t_events[0].size else (float(sol.t_events[0][0]), int(np.sign(sol.y_events[0][0][1])))
    return (t, sol, rhs) if dense else (t, sol.t.size - 1) if steps else t


def _route_rates(sol, rhs, ts, ss):
    """Along a long-double run's (α, σ) at the times ts: d(s) = ċ at δ = 0 and λ(s) = ∂ċ/∂c there (a central difference of 10⁻⁹σ* in
    the δ state)."""
    Y, hd = sol.sol(ts), 1e-9 * ss
    dd = np.array([rhs(t, np.r_[y[0], 0.0, y[2]])[1] for t, y in zip(ts, Y.T)])
    lam = np.array([(rhs(t, np.r_[y[0], hd, y[2]])[1] - rhs(t, np.r_[y[0], -hd, y[2]])[1]) / (2 * hd) for t, y in zip(ts, Y.T)])
    return dd, lam


def a2_sym_lag(ratio, wmul, rtol, n=4001):
    """Referee round 116 (D3dk-2): on exactly symmetric samples from their centre, Λ = ∫λ ds along the long-double run's route from the
    double exit (DOP853, the same rtol) to the long-double one — the linear escape's measure of the two seeds' ratio, against
    ln(2⁶³/2⁵²) = 11 ln 2 if each seed scales with its arithmetic's unit roundoff. None if either run does not leave."""
    td = a2_exit(ratio, wmul, sym=True, rtol=rtol, signed=True)[0]
    ex, sol, rhs = a2_exit_ld(ratio, wmul, sym=True, rtol=rtol, dense=True)
    if td is None or ex[0] is None:
        return None
    ts = np.linspace(td, ex[0], n)
    lam = _route_rates(sol, rhs, ts, ratio * S0)[1]
    return float(np.sum((lam[1:] + lam[:-1]) / 2 * np.diff(ts)))


def a2_seed(ratio, wmul, n=8001, rtol=1e-12):
    """Referee round 112 (D3dg-1): the data's seed integrated along the route — on A's samples in long double (`a2_exit_ld`), with
    d(s) = ċ at c = c* and λ(s) = ∂ċ/∂c there (a central difference of 10⁻⁹σ* in the δ state) on the run's (α, σ), Λ = ∫λ, the seed
    S = ∫_0^{s_e} e^{−Λ}d ds (trapezoid on n points): the offset at s = 0 that, carried by the linearisation alone, matches the data's
    drift. Returns (S, the exit (s, ±1), the symmetric problem's long-double exit from δ = S, Λ at the exit (round 116)); if the run does
    not leave, (None, the exit (None, 0), (None, 0), None) (round 125, D3dt-2)."""
    ex, sol, rhs = a2_exit_ld(ratio, wmul, rtol=rtol, dense=True)
    if ex[0] is None:
        return None, ex, (None, 0), None
    ss = ratio * S0
    ts = np.linspace(0.0, ex[0], n)
    dd, lam = _route_rates(sol, rhs, ts, ss)
    Lam = np.r_[0.0, np.cumsum((lam[1:] + lam[:-1]) / 2 * np.diff(ts))]
    g = np.exp(-Lam) * dd
    S = float(np.sum((g[1:] + g[:-1]) / 2 * np.diff(ts)))
    return S, ex, a2_exit_ld(ratio, wmul, sym=True, rtol=rtol, off=S / ss), float(Lam[-1])   # Λ at the exit (round 116, D3dk-2)


def a2_exit(ratio, wmul, sym=False, level=1e-2, smax=1.0, method='DOP853', rtol=1e-10, off=0.0, signed=False, summ=None, perm=None,
            max_step=np.inf, steps=False, atol=1e-14):
    """Referee round 103 (D3cx-2): the flow from c* exactly (w = wmul·A, σ = σ*) — by default with an explicit integrator (DOP853,
    rtol 10⁻¹⁰, atol 10⁻¹⁴) —: the s at which |c − c*| first reaches level·σ* (None by smax) — on A's samples, or (sym) on samples
    exactly symmetric about the centre (the problem translated to c* = 0). Round 106 (D3da-1, D3da-2): another method or rtol, a start
    off·σ* from c*, and (signed) the exit's direction, (s, ±1). Round 109 (D3dd-1): summ 'fsum' sums over the samples exactly rounded,
    'rev' in the reversed order — the data unchanged. Round 125 (D3dt-1): perm — the samples (and f*) in that order. Round 130
    (D3dy-2): the step at most max_step, and (steps) the exit with the solver's number of steps to it, (exit, n). Round 133 (D3eb-1):
    another atol."""
    from scipy.integrate import solve_ivp
    L, ss, d = law(), ratio * S0, XS[1] - XS[0]
    xs, c0 = (d * np.arange(-(XS.size // 2), XS.size // 2 + 1), 0.0) if sym else (XS, C_A)
    xs = xs[::-1] if summ == 'rev' else xs
    xs = xs if perm is None else xs[perm]
    K = _FsumBumps(xs, S0) if summ == 'fsum' else MW.BumpsCW(xs, S0)
    f = A_A * np.exp(-(xs - c0) ** 2 / (2 * ss * ss))
    z0 = np.r_[RT._alpha_of(L, np.array([wmul * A_A]), RT.coef_crit), c0 + off * ss, ss]
    ev = lambda t, y: abs(y[1] - c0) - level * ss
    ev.terminal = True
    sol = solve_ivp(RN.critical_rhs_nd(L, K, f), (0.0, smax), z0, method=method, rtol=rtol, atol=atol, events=ev, max_step=max_step)
    if not sol.t_events[0].size:
        t = (None, 0) if signed else None
    else:
        t = float(sol.t_events[0][0])
        t = (t, int(np.sign(sol.y_events[0][0][1] - c0))) if signed else t
    return (t, sol.t.size - 1) if steps else t


def a2_exits(ratio=0.25, wmuls=(3.0, 5.0), kzs=None):
    """Referee round 106 (D3da-1, D3da-2): per w the signed exit from c* — on A's samples by DOP853 at rtol 10⁻¹⁰ and 10⁻¹², by scipy's
    Radau and BDF at 10⁻¹⁰; on exactly symmetric samples by DOP853 at 10⁻¹⁰ and 10⁻¹², and (kzs[w], the zero of ċ in ulps at the held
    state: c* lies that far below the posed problem's repelling state) the symmetric problem started that far below its own. Returns
    {w: dict}."""
    out = {}
    for wm in wmuls:
        kz = (kzs or {}).get(wm)
        u = float(np.spacing(C_A))
        ss = ratio * S0
        r = dict(a=[a2_exit(ratio, wm, rtol=q, signed=True) for q in (1e-10, 1e-12)],
                 imp=[(m, a2_exit(ratio, wm, method=m, signed=True)) for m in ('Radau', 'BDF')]
                 + [('BDF at rtol 1e-12', a2_exit(ratio, wm, method='BDF', rtol=1e-12, signed=True))],   # round 116 (D3dk-6)
                 sym=[a2_exit(ratio, wm, sym=True, rtol=q, signed=True) for q in (1e-10, 1e-12)],
                 summ={k_: (a2_exit(ratio, wm, summ=k_, signed=True), a2_exit(ratio, wm, sym=True, summ=k_, signed=True))
                       for k_ in ('fsum', 'rev')},                     # round 109 (D3dd-1): the arithmetic changed, the data kept
                 kz=kz, symoff=None if kz is None else a2_exit(ratio, wm, sym=True, off=-kz * u / ss, signed=True))
        S, exl, syS, Le = a2_seed(ratio, wm)                            # round 112 (D3dg-1): long double and the route's seed
        r['ld'] = dict(a=[a2_exit_ld(ratio, wm, rtol=1e-10), exl], sym=[a2_exit_ld(ratio, wm, sym=True, rtol=q) for q in (1e-10, 1e-12)],
                       seed=S, symS=syS, ulp=u, Lam=Le, ss=ss, lag=[a2_sym_lag(ratio, wm, q) for q in (1e-10, 1e-12)])
        out[wm] = r
    return out


def a2_line_lam(ratio, wmul, rtol, smax, n=8001):
    """Referee round 125 (D3dt-1): Λ(0, s) = ∫λ along the line c = c* itself — the long-double flow on exactly symmetric samples with ċ
    set to 0 (`a2_exit_ld`'s line; the runs leave the line by 10⁻²σ*, which moves (α, σ) at second order), λ as `_route_rates`' — on
    n points of [0, smax] (the trapezoid). Returns (s, Λ)."""
    sol, rhs = a2_exit_ld(ratio, wmul, sym=True, smax=smax, rtol=rtol, dense=True, line=True)[1:]
    ts = np.linspace(0.0, smax, n)
    lam = _route_rates(sol, rhs, ts, ratio * S0)[1]
    return ts, np.r_[0.0, np.cumsum((lam[1:] + lam[:-1]) / 2 * np.diff(ts))]


def a2p_job(ratio, wmul, rtol, k=K_A2P, seed=SEED_A2P[0], max_step=None, atol_d=None):
    """Referee round 125 (D3dt-1): the exactly symmetric samples' exit from their centre as a draw of the rounding — over summation orders
    of the same data, in double (`a2_exit`, DOP853) and in long double (`a2_exit_ld`) at rtol: per order the signed exits, the effective
    seeds ln|S| = ln(10⁻²σ*) − Λ(0, s_e) (`a2_line_lam`) and the paired Λ = ln|S_d| − ln|S_ld| = Λ(s_d, s_ld) (`a2_sym_lag`'s in the
    code's order); u, the arithmetics' unit roundoffs. Round 130 (D3dy-1, D3dy-2): the orders are the code's (row 0, one draw) and then
    k permutations of the samples from default_rng(seed) (the distribution; a seed's first orders do not depend on k); per order
    DOP853's steps to each exit (nd, nld); DOP853's step at most max_step in both arithmetics (None: its own steps; the line's Λ at its
    own steps either way). Round 133 (D3eb-1): atol_d, the double run's atol (None: `a2_exit`'s 10⁻¹⁴; the long double's is 10⁻¹⁴σ*)."""
    rng = np.random.default_rng(seed)
    perms = [np.arange(XS.size)] + [rng.permutation(XS.size) for _ in range(k)]
    ms = np.inf if max_step is None else max_step
    d = [a2_exit(ratio, wmul, sym=True, rtol=rtol, signed=True, perm=p, max_step=ms, steps=True,
                 **({} if atol_d is None else dict(atol=atol_d))) for p in perms]
    l = [a2_exit_ld(ratio, wmul, sym=True, rtol=rtol, perm=p, max_step=ms, steps=True) for p in perms]
    ss = ratio * S0
    ts, Lam = a2_line_lam(ratio, wmul, rtol, 1.05 * max((t for (t, _), _ in d + l if t is not None), default=0.1))
    lnS = lambda t: None if t is None else float(math.log(1e-2 * ss) - np.interp(t, ts, Lam))
    rows = [dict(d=a, ld=b, nd=na, nld=nb, lnSd=lnS(a[0]), lnSld=lnS(b[0]),
                 lag=None if a[0] is None or b[0] is None else lnS(a[0]) - lnS(b[0])) for (a, na), (b, nb) in zip(d, l)]
    return dict(ratio=ratio, wmul=wmul, rtol=rtol, k=k, seed=seed, max_step=max_step, atol_d=atol_d, ss=ss,
                u=(float(np.finfo(float).eps) / 2, float(np.finfo(np.longdouble).eps) / 2), rows=rows)


def _a2p_job(args):
    return a2p_job(*args[:3], seed=args[3])


def part_A2P():
    """A2P (referee rounds 125, 130): per w, rtol and seed one record (the code's order and K_A2P random orders)."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(_a2p_job, [(0.25, w_, q_, s_) for w_ in WMUL_A2P for q_ in RTOL_A2P for s_ in SEED_A2P], chunksize=1)


def a2p_h0(A2P, wmul, rtol):
    """Referee round 130 (D3dy-2): h₀, the double runs' mean step s_e/n over A2P's random orders (both seeds) at (wmul, rtol)."""
    h = [q['d'][0] / q['nd'] for r in A2P if (r['wmul'], r['rtol'], r.get('max_step')) == (wmul, rtol, None)
         for q in r['rows'][1:] if q['d'][0] is not None]
    return float(np.mean(h))


def _a2h_job(args):
    return a2p_job(0.25, args[0], RTOL_A2H, k=K_A2H, seed=SEED_A2P[0], max_step=args[1])


def part_A2H(A2P):
    """A2H (referee round 130, D3dy-2): A2P's first K_A2H random orders of SEED_A2P[0] (and the code's) at RTOL_A2H with DOP853's step at
    most f·h₀ (f in FRAC_A2H) in both arithmetics — the seed's prefactor at other steps, the two arithmetics at equal ones."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(_a2h_job, [(w_, f_ * a2p_h0(A2P, w_, RTOL_A2H)) for w_ in WMUL_A2P for f_ in FRAC_A2H], chunksize=1)


def _a2r_job(args):
    return a2p_job(0.25, args[0], args[1], k=K_A2P if args[3] is None else K_A2H, seed=args[2], max_step=args[3], atol_d=args[4])


def a2r_jobs(A2P):
    """A2R's jobs (referee round 133, D3eb-1, -2): (w, rtol, seed, max_step, atol_d) — CELL_A2R on the seeds SEED_A2R at the default
    atols; CELL_A2R on every seed (A2P's and SEED_A2R) with the double run at ATOL_A2R·σ*, the long double's; A2H's cap CAP_A2R on
    SEED_A2P[1]."""
    w, q = CELL_A2R
    return ([(w, q, s_, None, None) for s_ in SEED_A2R] + [(w, q, s_, None, ATOL_A2R * (0.25 * S0)) for s_ in SEED_A2P + SEED_A2R]
            + [(CAP_A2R[0], RTOL_A2H, SEED_A2P[1], CAP_A2R[1] * a2p_h0(A2P, CAP_A2R[0], RTOL_A2H), None)])


def part_A2R(A2P):
    """A2R (referee round 133, D3eb-1, -2): the replications — more seeds at A2P's farthest cell, the atols matched, A2H's cap beyond
    3 sem on A2P's second seed."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(_a2r_job, a2r_jobs(A2P), chunksize=1)


def a2p_stats(rs):
    """A2P's distribution over the random orders (referee rounds 125, 130: D3dt-1, D3dy-1 … D3dy-3), of one record or of several pooled;
    the code's order (row 0) left out — one draw: the paired Λ's mean, sem (sd/√n, n the orders where both leave), t = (mean − ln(u_d/
    u_ld))/sem, median, quartiles and extremes; per arithmetic (double, long double) the seeds in units of uσ* — their geometric mean,
    quartiles and their logarithm's mean and sem —, the exits upward of those leaving, DOP853's mean number of steps to the exit and mean
    step s_e/n; the long double's lateness s_ld/s_d − 1: its median, quartiles and extremes."""
    rs = [rs] if isinstance(rs, dict) else list(rs)
    rows = [(r, q) for r in rs for q in r['rows'][1:]]
    lnr = math.log(rs[0]['u'][0] / rs[0]['u'][1])
    x = [q['lag'] for _, q in rows if q['lag'] is not None]
    lp = [[q[k_] - math.log(r['u'][i] * r['ss']) for r, q in rows if q[k_] is not None] for i, k_ in enumerate(('lnSd', 'lnSld'))]
    late = [q['ld'][0] / q['d'][0] - 1 for _, q in rows if q['d'][0] is not None and q['ld'][0] is not None]
    up = [(sum(1 for _, q in rows if q[k_][0] is not None and q[k_][1] > 0), sum(1 for _, q in rows if q[k_][0] is not None))
          for k_ in ('d', 'ld')]
    nh = [[(q[n_], q[k_][0] / q[n_]) for _, q in rows if q[k_][0] is not None] for k_, n_ in (('d', 'nd'), ('ld', 'nld'))]
    qs = lambda v: None if not v else dict(med=float(np.median(v)), q1=float(np.percentile(v, 25)), q3=float(np.percentile(v, 75)),
                                           lo=float(min(v)), hi=float(max(v)), n=len(v))
    se = lambda v: float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else None
    m, sem = (float(np.mean(x)) if x else None), se(x)
    return dict(n=len(x), mean=m, sem=sem, t=None if sem is None else (m - lnr) / sem, q=qs(x), lo=min(x, default=None),
                hi=max(x, default=None), gm=[math.exp(float(np.mean(v))) if v else None for v in lp],
                pq=[qs([math.exp(y) for y in v]) for v in lp], lp=[(float(np.mean(v)), se(v)) if v else None for v in lp], up=up,
                steps=[(float(np.mean([a for a, _ in v])), float(np.mean([h for _, h in v]))) if v else None for v in nh], late=qs(late))


Z99 = NormalDist().inv_cdf(0.995)                                       # a 99% interval: mean ± Z99·sem (referee round 133, D3eb-1)


def a2p_exp(s, lnr):
    """Referee round 133 (D3eb-1): the exponent Λ/(11 ln 2) of a2p_stats' s and its 99% interval (mean ± Z99·sem, over lnr); None
    without a sem."""
    return None if s['sem'] is None else (s['mean'] / lnr, (s['mean'] - Z99 * s['sem']) / lnr, (s['mean'] + Z99 * s['sem']) / lnr)


def a2p_blocked(A2P):
    """Referee round 133 (D3eb-5): the pooled Λ with its sem over orders — each seed's random order j is the same permutation in every
    cell, so the draws of one order are averaged over the cells where both arithmetics leave and the sem taken over the orders: (mean,
    sem, number of orders)."""
    bl = {}
    for r in A2P:
        for j, q in enumerate(r['rows'][1:]):
            if q['lag'] is not None:
                bl.setdefault((r['seed'], j), []).append(q['lag'])
    v = [float(np.mean(x)) for x in bl.values()]
    return float(np.mean(v)), (float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else None), len(v)


def a2p_gate(A2P, A2R=None):
    """Referee round 137 (D3ef-2): A2P's rule on all the draws of each cell — A2P's records with A2R's further seeds at default steps and
    atols (`a2r_line`'s 'more') at the cells they replicate: per cell (w, rtol) the pooled `a2p_stats`, its exponent with the 99%
    interval (`a2p_exp`) and the records' seeds; and the verdict — every interval inside EXP_A2P and every geometric mean within
    GM_A2P. Returns (verdict, {cell: (stats, exponent, seeds)})."""
    more = [r for r in (A2R if A2R is not None else []) if r['max_step'] is None and r.get('atol_d') is None]
    cells = list(dict.fromkeys((r['wmul'], r['rtol']) for r in A2P))
    lnr = math.log(A2P[0]['u'][0] / A2P[0]['u'][1])
    out = {}
    for c in cells:
        rs = [r for r in A2P if (r['wmul'], r['rtol']) == c] + [r for r in more if (r['wmul'], r['rtol']) == c]
        s = a2p_stats(rs)
        out[c] = (s, a2p_exp(s, lnr), [r['seed'] for r in rs])
    ok = all(e is not None and EXP_A2P[0] <= e[1] and e[2] <= EXP_A2P[1]
             and all(g is not None and GM_A2P[0] < g < GM_A2P[1] for g in s['gm']) for s, e, _ in out.values())
    return ok, out


def a2p_line(A2P, lag=None, A2R=None):
    """A2P's line (referee rounds 125, 130, 133: D3dt-1, D3dy-1 … D3dy-3, D3eb-1, D3eb-5): per w and rtol the paired Λ's mean ± sem per
    seed and pooled over the seeds, against ln(u_d/u_ld) = 11 ln 2 (in sem), the exponent Λ/(11 ln 2) with its 99% interval, the pooled
    draws' median, quartiles and extremes; the seeds' geometric mean and quartiles in units of uσ*, the exits upward, the long double's
    lateness (median, quartiles, extremes), DOP853's mean steps per arithmetic; over all the draws the exponent ± sem (and its sem over
    orders, `a2p_blocked`) and the exits upward per arithmetic; the conclusion only if every cell's exponent's 99% interval on all its
    draws (A2P with A2R's further seeds, `a2p_gate`: round 137, D3ef-2) lies inside EXP_A2P and every geometric mean within GM_A2P
    (round 133's rule: an interval, not a test of Λ = 11 ln 2, which asks the two arithmetics' prefactors to be equal — the cells beyond
    3 sem of it are named), the band's meaning and round 133's band's margin on A2P's draws at CELL_A2R printed. lag: {(w, rtol):
    `a2_sym_lag`'s Λ}, against the code order's."""
    cells = list(dict.fromkeys((r['wmul'], r['rtol']) for r in A2P))
    by = {c: [r for r in A2P if (r['wmul'], r['rtol']) == c] for c in cells}
    st, al = {c: a2p_stats(by[c]) for c in cells}, a2p_stats(A2P)
    dl = [abs(r['rows'][0]['lag'] - lag[(r['wmul'], r['rtol'])]) for r in A2P
          if (lag or {}).get((r['wmul'], r['rtol'])) is not None and r['rows'][0]['lag'] is not None]
    lnr = math.log(A2P[0]['u'][0] / A2P[0]['u'][1])
    ex = {c: a2p_exp(s, lnr) for c, s in st.items()}
    ok, gt = a2p_gate(A2P, A2R)                                         # round 137 (D3ef-2): the rule on all of each cell's draws
    off = [c for c, s in st.items() if s['t'] is not None and abs(s['t']) >= 3]
    tm = max((c for c in cells if st[c]['t'] is not None), key=lambda c: abs(st[c]['t']), default=None)
    bm, bs, bn = a2p_blocked(A2P)
    seeds = list(dict.fromkeys(r['seed'] for r in A2P))
    pw = lambda u: f"2^{{{math.log2(u):.0f}}}"
    gm = lambda g, q: '—' if g is None else f"{g:.2f}" + ('' if q is None else f" (interquartile {q['q1']:.2g}–{q['q3']:.2g})")
    ms = lambda s: '—' if s['sem'] is None else f"{s['mean']:.2f} ± {s['sem']:.2f}"
    ts = lambda s: '—' if s['t'] is None else f"{s['t']:+.1f}"
    xs = lambda e: '' if e is None else f", the exponent {e[0]:.3f} (99% interval {e[1]:.3f} … {e[2]:.3f})"
    late = lambda q: ('' if q is None else f", the long double later by a median {100 * q['med']:.0f}% (interquartile {100 * q['q1']:.0f}–"
                      f"{100 * q['q3']:.0f}%, extremes {100 * q['lo']:.1f}% … {100 * q['hi']:.1f}% of {q['n']})")
    stp = lambda s: ('' if None in s['steps'] else f", DOP853's {s['steps'][0][0]:.0f} and {s['steps'][1][0]:.0f} steps to the exit on "
                     f"average (mean step {s['steps'][0][1]:.2e} and {s['steps'][1][1]:.2e})")
    cell = lambda c: (lambda ps, s: f"{c[0]:g}A at rtol {c[1]:g}: Λ " + ' and '.join(ms(x) for x in ps) + " (seed" + (
        's ' if len(ps) > 1 else ' ') + ' and '.join(str(r['seed']) for r in by[c]) + ", over " + ' and '.join(str(x['n']) for x in ps)
        + ": " + ' and '.join(ts(x) for x in ps) + " sem from 11 ln 2), pooled " + ms(s) + f" (over {s['n']}, {ts(s)} sem" + (
        '' if s['q'] is None else f"; median {s['q']['med']:.2f}, interquartile {s['q']['q1']:.2f} … {s['q']['q3']:.2f}, extremes "
        f"{s['q']['lo']:.2f} … {s['q']['hi']:.2f} of {s['q']['n']}") + ")" + xs(ex[c]) + f", the seeds' geometric mean "
        f"{gm(s['gm'][0], s['pq'][0])} u_dσ* in double and {gm(s['gm'][1], s['pq'][1])} u_ldσ* in long double, upward in {s['up'][0][0]} "
        f"of {s['up'][0][1]} (double) and {s['up'][1][0]} of {s['up'][1][1]} (long double)" + late(s['late']) + stp(s))(
        [a2p_stats(r) for r in by[c]], st[c])
    pc = lambda u: f"{u[0]} of {u[1]}" + ('' if not u[1] else f" ({100 * u[0] / u[1]:.0f}%)")
    pt = max(((a2p_stats(r), r) for r in A2P), key=lambda x: -1 if x[0]['t'] is None else abs(x[0]['t']))   # the per-seed largest offset
    return (f"The symmetric samples' exit as a draw (A2P; referee rounds 125, 130, 133, 137): over summation orders of the same data — the "
            f"exactly symmetric samples in {' and '.join(str(r_) for r_ in sorted({r['k'] for r in A2P}))} random permutations from each "
            f"of the seeds {' and '.join(str(x) for x in seeds)} (drawn in round 130, after round 125's rule was set: round 125's stored "
            "draws, seed 125, were its referee's own) and in the code's order, one draw apart (" + (
                f"its Λ within {max(dl):.1e} of `a2_sym_lag`'s, along the long-double run's own route" if dl else
                "`a2_sym_lag`'s") + ") —, DOP853 from their centre in double and in long double at its own steps, each exit's effective "
            "seed |S| = 10⁻²σ*·e^{−Λ(0, s_e)} (Λ along the line c = c*, `a2_line_lam`) in units of uσ*, u the unit roundoff "
            f"({pw(A2P[0]['u'][0])} in double, {pw(A2P[0]['u'][1])} in long double), and per order the paired Λ = ln(|S_d|/|S_ld|), "
            "from the double exit to the long-double one, and the exponent Λ/(11 ln 2) with its 99% interval (mean ± "
            f"{Z99:.3f} sem): " + '; '.join(cell(c) for c in cells)
            + f" — against ln(u_d/u_ld) = 11 ln 2 = {lnr:.2f}, the unit roundoffs' ratio; over all {al['n']} draws Λ = "
            + ('—' if al['sem'] is None else f"{al['mean']:.3f} ± {al['sem']:.3f}, the exponent Λ/(11 ln 2) = {al['mean'] / lnr:.3f} ± "
               f"{al['sem'] / lnr:.3f}" + ('' if bs is None else f" (± {bs / lnr:.4f} over the {bn} orders, each order's cells "
                                           "averaged: the same orders recur in every cell)"))
            + f", the exits upward in {pc(al['up'][0])} in double and {pc(al['up'][1])} in long double"
            + ('' if pt[0]['t'] is None else f", the per-seed means' largest offset {pt[0]['t']:+.1f} sem ({pt[1]['wmul']:g}A at rtol "
               f"{pt[1]['rtol']:g}, seed {pt[1]['seed']})")
            + "; the rule on all of each cell's draws (" + ("A2P with A2R's further seeds" if any(
                len(sd) > len(by[c]) for c, (_, _, sd) in gt.items()) else "A2P's") + ": " + ', '.join(
                f"{c[0]:g}A at rtol {c[1]:g} over {s['n']}" + ('' if e is None else f", {e[1]:.3f} … {e[2]:.3f}") for c, (s, e, _) in gt.items())
            + f") with the band [{EXP_A2P[0]:g}, {EXP_A2P[1]:g}] — halfway to the neighbouring exponents ½ and 3/2; at γ = 1 the two "
            f"prefactors' ratio within e^{{±{(EXP_A2P[1] - 1) * lnr:.2f}}}, while each geometric mean within {GM_A2P[0]:g}–{GM_A2P[1]:g} "
            "uσ* asks each prefactor O(1); round 133's band [" + f"{EXP_A2P_133[0]:g}, {EXP_A2P_133[1]:g}], set after the data (e^{{±"
            f"{(EXP_A2P_133[1] - 1) * lnr:.2f}}})" + ('' if ex.get(CELL_A2R) is None else ': ' + (lambda m, hw: (
                f"A2P's {st[CELL_A2R]['n']} alone passed at {CELL_A2R[0]:g}A, rtol {CELL_A2R[1]:g} by {m:.4f} at its upper end, against "
                f"the interval's half-width {hw:.4f}" + (", marginal" if m < 0.1 * hw else '')) if m >= 0 else
                f"A2P's {st[CELL_A2R]['n']} alone failed it at {CELL_A2R[0]:g}A, rtol {CELL_A2R[1]:g}")(
                EXP_A2P_133[1] - ex[CELL_A2R][2], ex[CELL_A2R][2] - ex[CELL_A2R][0])) + " —: " + (
                "every cell's interval inside the band and every geometric mean within it, so the seed scales with the arithmetic's "
                "unit roundoff at DOP853's steps (A2H, below) — the two seeds' ratio the unit roundoffs' to within the band, the "
                "exponent and the prefactors' ratio not separable from two arithmetics (A2T, below: three) —, while the direction, Λ "
                "and the lateness of one order are draws" if ok else
                "not every cell inside the band, or a geometric mean outside it")
            + "; Λ = 11 ln 2 itself asks the two arithmetics' prefactors to be equal — " + (
                f"{len(off)} of {len(cells)} cells beyond 3 sem of it (" + ', '.join(f"{c[0]:g}A at rtol {c[1]:g}, {st[c]['t']:+.1f}"
                                                                                      for c in off) + "): there they differ"
                if off else "no cell beyond 3 sem of it" + ('' if tm is None else f" (the largest {st[tm]['t']:+.1f}, {tm[0]:g}A at "
                                                                                   f"rtol {tm[1]:g}: A2R, below)")) + ".")


def a2p_shift(r, base, k_):
    """Referee rounds 130, 133 (D3dy-2, D3eb-3): the shift of ln|S|'s mean (k_ 'lnSd' or 'lnSld') from base's run to r's over the same
    orders (rows 1 … n, the random orders, where both leave): the difference of the means, its sem (the two sems in quadrature: an
    order's draws at different steps are independent, their correlation r returned beside) and that correlation, (shift, sem, r)."""
    v = np.array([(a[k_], b[k_]) for a, b in zip(r['rows'][1:], base['rows'][1:]) if a[k_] is not None and b[k_] is not None])
    se = lambda x: np.var(x, ddof=1) / len(x)
    return float(np.mean(v[:, 0]) - np.mean(v[:, 1])), float(math.sqrt(se(v[:, 0]) + se(v[:, 1]))), float(np.corrcoef(v.T)[0, 1])


def dl_meta(x, e):
    """Referee round 137 (D3ef-5): DerSimonian–Laird random effects of the estimates x with sems e — the fixed-effect mean, Cochran's
    Q on len(x) − 1 dof, the between-estimate sd τ (τ² = max(0, (Q − k + 1)/(Σw − Σw²/Σw)), w = 1/e²), the random-effects mean and
    its sem (weights 1/(e² + τ²)). Returns (fixed mean, Q, dof, τ, mean, sem)."""
    x, e = np.asarray(x, float), np.asarray(e, float)
    w = 1 / e ** 2
    mu = float((w * x).sum() / w.sum())
    Q = float((w * (x - mu) ** 2).sum())
    t2 = max(0.0, (Q - (x.size - 1)) / float(w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (e ** 2 + t2)
    return mu, Q, x.size - 1, math.sqrt(t2), float((wr * x).sum() / wr.sum()), float(wr.sum() ** -0.5)


def a2h_line(A2H, A2P):
    """A2H's line (referee rounds 130, 133: D3dy-2, D3eb-2, D3eb-3): per w DOP853's own steps (A2P's record of that seed and rtol, its
    first orders) and each cap (A2H): the mean step and steps per arithmetic, ln(S/uσ*)'s mean ± sem, at a cap the shift of its mean
    from the same orders at DOP853's own steps (`a2p_shift`) against Brouwer's law S ∝ h^{1/2} (½ ln f at the cap f·h₀: in sem), the
    paired Λ's mean ± sem (in sem from 11 ln 2, and the exponent) and the exits upward; then the shifts' extremes, how many lie beyond 3
    sem of 0 and of Brouwer's law, the correlations of an order's draws, at DOP853's own steps the two arithmetics' mean steps' largest
    relative difference over A2P's cells, the caps beyond 3 sem of 11 ln 2, and the exponent's range over every configuration (A2P's
    cells and the caps) with the estimates' sems and their DerSimonian–Laird heterogeneity (`dl_meta`: Q, τ, the random-effects
    exponent; round 137, D3ef-5) beside the pooled exponent's statistical sem."""
    lnr = math.log(A2H[0]['u'][0] / A2H[0]['u'][1])
    sh, bw, parts, cap = [], [], [], []
    lp = lambda s, i: f"{s['lp'][i][0]:+.2f} ± {s['lp'][i][1]:.2f}"
    lam = lambda s: f"Λ {s['mean']:.2f} ± {s['sem']:.2f} ({s['t']:+.1f} sem from 11 ln 2, the exponent {s['mean'] / lnr:.3f} ± {s['sem'] / lnr:.3f})"
    ar = ('double', 'long double')
    for w in dict.fromkeys(r['wmul'] for r in A2H):
        rs = [r for r in A2H if r['wmul'] == w]
        base = next(x for x in A2P if (x['wmul'], x['rtol'], x['seed'], x.get('max_step')) == (w, rs[0]['rtol'], rs[0]['seed'], None))
        base = dict(base, rows=base['rows'][:len(rs[0]['rows'])])      # the same orders at DOP853's own steps
        b, h0 = a2p_stats(base), a2p_h0(A2P, w, rs[0]['rtol'])
        sub = []
        for r in rs:
            a, d, f = a2p_stats(r), [a2p_shift(r, base, k_) for k_ in ('lnSd', 'lnSld')], r['max_step'] / h0
            sh += d
            bw += [((m - 0.5 * math.log(f)) / e, f"{w:g}A at {f:g}·h₀ in {ar[i]}") for i, (m, e, _) in enumerate(d)]
            cap.append((a, f"{w:g}A at {f:g}·h₀"))
            sub.append(f"at most {f:g}·h₀: mean step {a['steps'][0][1]:.3e} and {a['steps'][1][1]:.3e} ({a['steps'][0][0]:.0f} "
                       f"and {a['steps'][1][0]:.0f} steps), ln(S/uσ*) {lp(a, 0)} and {lp(a, 1)}, its mean moved by {d[0][0]:+.2f} ± "
                       f"{d[0][1]:.2f} and {d[1][0]:+.2f} ± {d[1][1]:.2f} (Brouwer's {0.5 * math.log(f):+.2f}: {bw[-2][0]:+.2f} and "
                       f"{bw[-1][0]:+.2f} sem from it), {lam(a)}, upward in {a['up'][0][0]} and {a['up'][1][0]} of {a['n']}")
        parts.append(f"{w:g}A (h₀ {h0:.2e}): at DOP853's own steps mean step {b['steps'][0][1]:.3e} and {b['steps'][1][1]:.3e} (double, long "
                     f"double; {b['steps'][0][0]:.0f} and {b['steps'][1][0]:.0f} steps), ln(S/uσ*) {lp(b, 0)} and {lp(b, 1)}, {lam(b)}; "
                     + '; '.join(sub))
    cs = [a2p_stats([x for x in A2P if (x['wmul'], x['rtol']) == c]) for c in dict.fromkeys((x['wmul'], x['rtol']) for x in A2P)]
    dh = max(abs(s['steps'][1][1] / s['steps'][0][1] - 1) for s in cs)
    big = sum(1 for m, e, _ in sh if abs(m) > 3 * e)
    nb = [lab for z, lab in bw if abs(z) > 3]
    off = [(a, lab) for a, lab in cap if abs(a['t']) >= 3]
    xs = [s['mean'] / lnr for s in cs] + [a['mean'] / lnr for a, _ in cap]
    es = [s['sem'] / lnr for s in cs] + [a['sem'] / lnr for a, _ in cap]
    dm = dl_meta(xs, es)                                                # round 137 (D3ef-5): the heterogeneity's size
    al = a2p_stats(A2P)
    return (f"The steps (A2H; referee rounds 130, 133, 137): on A2P's first {A2H[0]['k']} random orders of seed {A2H[0]['seed']} at rtol "
            f"{A2H[0]['rtol']:g}, DOP853 with its step at most " + ' and '.join(f"{f_:g}·h₀" for f_ in FRAC_A2H) + " in both "
            "arithmetics, h₀ the double's mean step over A2P's draws at that rtol (the line's Λ at its own steps), ln(S/uσ*) as mean ± sem "
            "and the shift of its mean from DOP853's own steps on the same orders (± the two sems in quadrature: the orders' draws are "
            f"independent at different steps, their correlation {min(c for _, _, c in sh):+.2f} … {max(c for _, _, c in sh):+.2f}), "
            "against Brouwer's law for rounding in the derivative, S ∝ h^{1/2} (a shift ½ ln f at the cap f·h₀): " + '; '.join(parts) + " — "
            + (f"the seed's prefactor in units of uσ* moves with the step sequence (the shifts of the means {min(m for m, _, _ in sh):+.2f} … "
               f"{max(m for m, _, _ in sh):+.2f} in ln, {big} of {len(sh)} beyond 3 sem), the integrator sampling the rounding" if big else
               f"the shifts ({min(m for m, _, _ in sh):+.2f} … {max(m for m, _, _ in sh):+.2f} in ln) within 3 sem: no step dependence resolved")
            + (f", {len(nb)} of {len(bw)} beyond 3 sem of Brouwer's law (" + ', '.join(nb) + ")" if nb else
               f", every shift within 3 sem of Brouwer's law")
            + f"; at DOP853's own steps the two arithmetics' mean steps differ by at most {100 * dh:.0f}% (A2P's cells), under a cap both "
            f"step at it, and there Λ is {min(a['mean'] for a, _ in cap):.2f} … {max(a['mean'] for a, _ in cap):.2f} (the exponent "
            f"{min(a['mean'] for a, _ in cap) / lnr:.3f} … {max(a['mean'] for a, _ in cap) / lnr:.3f}, {len(off)} of {len(cap)} beyond 3 "
            "sem of 11 ln 2" + (": " + ', '.join(f"{lab}, {a['t']:+.1f}" for a, lab in off) + " — the two arithmetics' prefactors "
                                "differ at equal steps too" if off else '') + "), the step sequence moving the two arithmetics' "
            f"prefactors in part alike; over every configuration ({len(cs)} cells at DOP853's own steps, A2P; {len(cap)} caps) the exponent "
            f"Λ/(11 ln 2) spans {min(xs):.3f} … {max(xs):.3f} — a range of {len(xs)} estimates of sems {min(es):.3f}–{max(es):.3f}; their "
            f"heterogeneity Q = {dm[1]:.1f} on {dm[2]} dof, the between-configuration sd τ = {dm[3]:.3f} (DerSimonian–Laird), the "
            f"random-effects exponent {dm[4]:.3f} ± {dm[5]:.3f} over configurations —, beside the pooled exponent's statistical ± "
            f"{al['sem'] / lnr:.3f} (A2P).")


def a2_paired_shift(rs, bases, k_='lag'):
    """Referee round 137 (D3ef-4): the paired shift of k_ from each base record to the record of the same seed in rs over the same
    orders (rows 1 …, where both are set): (mean, sem, n)."""
    bd = {b['seed']: b for b in bases}
    v = np.array([a[k_] - b[k_] for r in rs for a, b in zip(r['rows'][1:], bd[r['seed']]['rows'][1:])
                  if a[k_] is not None and b[k_] is not None])
    return float(np.mean(v)), float(np.std(v, ddof=1) / math.sqrt(v.size)), int(v.size)


def a2r_line(A2R, A2P, A2H):
    """A2R's line (referee rounds 133, 137: D3eb-1, D3eb-2, D3ef-4): at CELL_A2R — selected as A2P's cell farthest from 11 ln 2 — the
    further seeds' Λ (per seed, over them alone, and over every seed of the cell with A2P's selected draws), the exponent with its 99%
    interval, each arithmetic's ln(S/uσ*) and Λ − 11 ln 2, their difference, as mean ± sem (t to two decimals); the same orders with
    the double at the long double's atol (the steps per arithmetic, the further seeds alone, and the paired shift from the double's own
    atol: the long double's runs are the same); A2H's cap CAP_A2R — selected as A2H's farthest — on SEED_A2P[1] alone and pooled with
    A2H's SEED_A2P[0]. The verdicts on the further draws alone: "differ" only beyond 3 sem of 11 ln 2 there; the two replications
    combined by Stouffer's (t₁ + t₂)/√2; "not the atol" if the paired shift lies within 3 sem of 0."""
    lnr = math.log(A2R[0]['u'][0] / A2R[0]['u'][1])
    more = [r for r in A2R if r['max_step'] is None and r.get('atol_d') is None]
    mat = [r for r in A2R if r.get('atol_d') is not None]
    cap = [r for r in A2R if r['max_step'] is not None]
    w, q = more[0]['wmul'], more[0]['rtol']
    base = [r for r in A2P if (r['wmul'], r['rtol']) == (w, q)]
    ms = lambda s: f"{s['mean']:.3f} ± {s['sem']:.3f} ({s['t']:+.2f} sem)"
    xs = lambda s: (lambda e: f"the exponent {e[0]:.3f} (99% interval {e[1]:.3f} … {e[2]:.3f})")(a2p_exp(s, lnr))
    lp = lambda s: f"ln(S/uσ*) {s['lp'][0][0]:+.2f} ± {s['lp'][0][1]:.2f} in double and {s['lp'][1][0]:+.2f} ± {s['lp'][1][1]:.2f} in long double"
    sm, sa, sx = a2p_stats(more), a2p_stats(base + more), a2p_stats(mat)
    sxm = a2p_stats([r for r in mat if r['seed'] not in {b['seed'] for b in base}])
    out = (f"The replications (A2R; referee rounds 133, 137): at {w:g}A, rtol {q:g} — selected as A2P's cell farthest from 11 ln 2 (of "
           f"{len({(r['wmul'], r['rtol']) for r in A2P})}, on its {a2p_stats(base)['n']} draws: the selection's) — {more[0]['k']} random "
           "orders from each of the further seeds " + ', '.join(str(r['seed']) for r in more) + ": Λ " + ', '.join(
               f"{s['mean']:.2f} ± {s['sem']:.2f}" for s in map(a2p_stats, more)) + f", over their {sm['n']} alone {ms(sm)}; with the "
           f"selected {sa['n'] - sm['n']}, over all {sa['n']} Λ = {ms(sa)} from 11 ln 2, {xs(sa)}, {lp(sa)} — their difference Λ − 11 "
           f"ln 2 = {sa['mean'] - lnr:+.2f} ± {sa['sem']:.2f}" + (
               f": on the further draws alone the two arithmetics' runs' prefactors differ at this cell, by about e^{{{sa['mean'] - lnr:+.2f}}}"
               " — at DOP853's step sequences (A2T, below: at 3A the exit is the integrator's escape from ulp-scale traps)"
               if abs(sm['t']) >= 3 else f": on the further draws alone within 3 sem ({sm['t']:+.2f}), not resolved by this replication alone"))
    nd = lambda rs, k_: float(np.mean([x[k_] for r in rs for x in r['rows'][1:]]))
    ps = a2_paired_shift(mat, base + more)
    out += (f"; with the double run at the long double's atol ({mat[0]['atol_d']:.1e} = {ATOL_A2R:g}·σ*, against 1e-14) on the same "
            f"{sx['n']} orders (seeds " + ', '.join(str(r['seed']) for r in mat) + f") Λ = {ms(sx)}, on the further seeds' {sxm['n']} alone "
            f"{ms(sxm)}, DOP853's double taking {nd(mat, 'nd'):.2f} steps against {nd(base + more, 'nd'):.2f} at its own atol (the long "
            f"double {nd(base + more, 'nld'):.2f}), its mean step {sx['steps'][0][1]:.3e} against {sa['steps'][0][1]:.3e} (the long "
            f"double's {sa['steps'][1][1]:.3e}); Λ's paired shift from the double's own atol on the same {ps[2]} orders (the long double's "
            f"runs the same) {ps[0]:+.3f} ± {ps[1]:.3f}" + (": not the atol" if abs(ps[0]) < 3 * ps[1] else ": the atol moves Λ"))
    tr = [sm['t']]
    for r in cap:
        h = next(x for x in A2H if (x['wmul'], x['rtol']) == (r['wmul'], r['rtol']) and abs(x['max_step'] / r['max_step'] - 1) < 1e-12)
        s1, sp = a2p_stats(r), a2p_stats([h, r])
        tr.append(s1['t'])
        f = r['max_step'] / a2p_h0(A2P, r['wmul'], r['rtol'])
        out += (f"; A2H's cap at {f:g}·h₀ at {r['wmul']:g}A — selected as A2H's farthest from 11 ln 2 (of {len(A2H)}, "
                f"{a2p_stats(h)['t']:+.2f} sem on seed {h['seed']}: the selection's) — on seed {r['seed']}'s first {r['k']} orders alone: "
                f"Λ {ms(s1)}, {lp(s1)}; with the selected seed's, over both seeds' {sp['n']} Λ = {ms(sp)}, {xs(sp)}" + (
                    ": replicated, the arithmetics' runs' prefactors differ at equal steps too" if abs(s1['t']) >= 3 else
                    f": on the further draws alone within 3 sem ({s1['t']:+.2f}), not replicated beyond 3 sem by itself"))
    if len(tr) > 1:
        z = sum(tr) / math.sqrt(len(tr))
        out += (f"; the {len(tr)} replications alone combined (Stouffer, Σt/√{len(tr)}) {z:+.2f} sem" + (
            ": the two arithmetics' runs' prefactors differ, on the further draws alone" if abs(z) >= 3 else
            ": not resolved on the further draws alone"))
    return out + "."


class _DTBumps(_LDBumps):
    """`_LDBumps` in any floating dtype (referee round 137, D3ef-3): the samples' offsets, the Gaussians, products and sums in dt, the
    correlations returned in double. At np.longdouble its right-hand side is `_LDBumps`' bitwise, at np.float64 (cref = 0) `BumpsCW`'s
    bitwise; at np.float32 a third arithmetic of the same data (u = 2⁻²⁴)."""

    def __init__(self, xs, sig, cref, dt):
        self.dt = dt
        self.xr = np.asarray(xs, float).astype(dt) - dt(float(cref))
        self.N, self.n, self.sig = self.xr.size, 2, float(sig)

    def _parts(self, th):
        dt = self.dt
        th = np.atleast_2d(np.asarray(th, float))
        d, s = th[:, 0].astype(dt), th[:, 1].astype(dt)
        u = self.xr[:, None] - d[None, :]
        return u, s, np.exp(-u ** 2 / (2 * s ** 2))

    def synth(self, th, w):
        return self._parts(th)[2] @ np.asarray(w, float).astype(self.dt)


def a2t_rhs(wmul, dt, perm=None, ratio=0.25):
    """Referee round 137 (D3ef-1, -3): on the exactly symmetric samples (in perm's order) with the centre's offset δ as the state, the
    block's right-hand side in the arithmetic dt (`_DTBumps`; f* the stored double data cast to dt) and the start (w = wmul·A, δ = 0,
    σ = σ*). Returns (rhs, z0, σ*)."""
    L, ss, d = law(), ratio * S0, XS[1] - XS[0]
    xs = d * np.arange(-(XS.size // 2), XS.size // 2 + 1)
    xs = xs if perm is None else xs[perm]
    f = (A_A * np.exp(-xs ** 2 / (2 * ss * ss))).astype(dt)
    z0 = np.r_[RT._alpha_of(L, np.array([wmul * A_A]), RT.coef_crit), 0.0, ss]
    return RN.critical_rhs_nd(L, _DTBumps(xs, S0, 0.0, dt), f), z0, ss


def a2t_rk4(wmul, dt, perm, h, smax=SMAX_A2T, level=1e-2):
    """Referee round 137 (D3ef-1, -3): classical RK4 at the fixed step h (no error control: the same step sequence in every arithmetic)
    from the symmetric samples' centre in the arithmetic dt: the exit — |δ| first reaching level·σ*, s by linear interpolation of ln|δ|
    over the step — and its direction, (s, ±1), or (None, 0) by smax; and the number of steps taken."""
    F, z, ss = a2t_rhs(wmul, dt, perm)
    t, lev, n = 0.0, level * ss, 0
    while t < smax:
        k1 = F(t, z)
        k2 = F(t + h / 2, z + h / 2 * k1)
        k3 = F(t + h / 2, z + h / 2 * k2)
        k4 = F(t + h, z + h * k3)
        zn = z + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        n += 1
        if abs(zn[1]) >= lev:
            a, b = abs(z[1]), abs(zn[1])
            fr = (math.log(lev) - math.log(a)) / (math.log(b) - math.log(a)) if a > 0 else 1.0
            return (t + fr * h, int(np.sign(zn[1]))), n
        z, t = zn, t + h
    return (None, 0), n


def a2t_states(wmul, ss_, rtol=1e-12):
    """Referee round 137 (D3ef-1): the long-double line flow's (α, σ) at the times ss_ (exactly symmetric samples, ċ set to 0:
    `a2_exit_ld`'s line), the states of the scan."""
    sol = a2_exit_ld(0.25, wmul, sym=True, smax=max(ss_), rtol=rtol, dense=True, line=True)[1]
    return [(float(sol.sol(s)[0]), float(sol.sol(s)[2])) for s in ss_]


def a2t_scan(wmul, al, sg, dt, perm=None, nt=NT_A2T, nu=ULP_A2T):
    """Referee round 137 (D3ef-1): the computed ċ(δ) at the state (α, σ) on nt points of |δ| ≤ nu ulps of the ±1 samples' coordinate
    x₁ (ulp(x₁) in dt) against the linear law λδ (λ by a central difference of ±10⁻⁶σ*, where the rounding is negligible): λ, ulp(x₁)
    over uσ*, the attracting zeros (ċ from + to − as δ grows) and the repelling ones (− to +) in ulps (the midpoint of the grid's
    bracket), the slope of ċ over |δ| < 0.45 ulp in units of λ (a least-squares line), and the rms of ċ − λδ over λuσ*."""
    F, z0, ss = a2t_rhs(wmul, dt, perm)
    u = float(np.finfo(dt).eps) / 2
    cd = lambda x: float(F(0.0, np.array([al, x, sg]))[1])
    big = 1e-6 * ss
    lam = (cd(big) - cd(-big)) / (2 * big)
    ul = float(np.spacing(dt(XS[1] - XS[0])))
    g = np.linspace(-nu * ul, nu * ul, nt)
    v = np.array([cd(x) for x in g])
    at = [float((g[i] + g[i + 1]) / 2 / ul) for i in np.where((v[:-1] > 0) & (v[1:] <= 0))[0]]
    rp = [float((g[i] + g[i + 1]) / 2 / ul) for i in np.where((v[:-1] < 0) & (v[1:] >= 0))[0]]
    m = np.abs(g) < 0.45 * ul
    return dict(lam=float(lam), ulp=ul / (u * ss), att=at, rep=rp, slope0=float(np.polyfit(g[m], v[m], 1)[0] / lam),
                rms=float(np.sqrt(np.mean((v - lam * g) ** 2)) / (lam * u * ss)))


def a2t_job(args):
    """A2T's job (referee round 137, D3ef-1, -3): ('scan', w, the line's times) — the line's states (`a2t_states`) and, at each, the
    scan of ċ(δ) in every arithmetic of DT_A2T on the code's order and SEED_A2T's first two random orders; ('rk', w, f, k, seed, h₀) —
    k random orders from default_rng(seed) at the fixed RK4 step f·h₀ in every arithmetic (`a2t_rk4`), each exit's effective seed
    ln|S| = ln(10⁻²σ*) − Λ(0, s_e) (`a2_line_lam` at rtol 10⁻¹²)."""
    ss = 0.25 * S0
    if args[0] == 'scan':
        w, ts = args[1], args[2]
        rng = np.random.default_rng(SEED_A2T)
        perms = [None] + [rng.permutation(XS.size) for _ in range(2)]
        rows = [dict(s=s_, al=a, sig=sg / ss, scans=[dict(arith=nm, order=j, **a2t_scan(w, a, sg, dt, p)) for nm, dt in DT_A2T
                                                     for j, p in enumerate(perms)]) for s_, (a, sg) in zip(ts, a2t_states(w, ts))]
        return dict(kind='scan', wmul=w, rows=rows)
    w, f, k, seed, h0 = args[1:]
    rng = np.random.default_rng(seed)
    perms = [rng.permutation(XS.size) for _ in range(k)]
    rows = []
    for p in perms:
        r = {}
        for nm, dt in DT_A2T:
            (s_, dr), n = a2t_rk4(w, dt, p, f * h0)
            r[nm] = dict(s=s_, dir=dr, n=n)
        rows.append(r)
    sm = max((x['s'] for r in rows for x in r.values() if x['s'] is not None), default=None)
    if sm is not None:
        ts, Lam = a2_line_lam(0.25, w, 1e-12, 1.05 * sm)
        for r in rows:
            for x in r.values():
                x['lnS'] = None if x['s'] is None else float(math.log(1e-2 * ss) - np.interp(x['s'], ts, Lam))
    return dict(kind='rk', wmul=w, frac=f, h=f * h0, k=k, seed=seed, smax=SMAX_A2T, ss=ss,
                u={nm: float(np.finfo(dt).eps) / 2 for nm, dt in DT_A2T}, rows=rows)


def _a2t_job(args):
    return a2t_job(args)


def a2t_jobs(A2P):
    """A2T's jobs: a scan per w of S_A2T, an RK4 arm per RK_A2T entry (h₀ the double's mean DOP853 step at rtol 10⁻¹⁰, `a2p_h0`)."""
    return ([('scan', w, ts) for w, ts in S_A2T]
            + [('rk', w, f, k, s_, a2p_h0(A2P, w, RTOL_A2H)) for w, f, k, s_ in RK_A2T])


def part_A2T(A2P):
    """A2T (referee round 137, D3ef-1, -3): the computed field at the ulp scale and three arithmetics at identical fixed steps."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(_a2t_job, a2t_jobs(A2P), chunksize=1)


def a2t_rk_stats(r):
    """A2T's RK4 arm r (referee round 137, D3ef-1, -3): per arithmetic the exits by smax of the k orders, the upward ones, the mean s_e,
    ln(S/uσ*)'s mean, sem and median; over the orders where all three leave the three-point slope γ of ln|S| against ln u per order
    (least squares), its mean, sem, 99% interval and median, the pairwise exponents (ln|S_a| − ln|S_b|)/(ln u_a − ln u_b) and the
    prefactors' log-ratios ln(C_a/C_ld) (C = S/(uσ*)), as mean ± sem."""
    nms = [nm for nm, _ in DT_A2T]
    lu = {nm: math.log(r['u'][nm]) for nm in nms}
    se = lambda v: float(np.std(v, ddof=1) / math.sqrt(len(v))) if len(v) > 1 else None
    per = {}
    for nm in nms:
        ex = [x[nm] for x in r['rows'] if x[nm]['s'] is not None]
        v = [x['lnS'] - lu[nm] - math.log(r['ss']) for x in ex]
        per[nm] = dict(n=len(ex), up=sum(1 for x in ex if x['dir'] > 0), s=float(np.mean([x['s'] for x in ex])) if ex else None,
                       mean=float(np.mean(v)) if v else None, sem=se(v), med=float(np.median(v)) if v else None)
    ok = [x for x in r['rows'] if all(x[nm]['s'] is not None for nm in nms)]
    out = dict(per=per, n3=len(ok))
    if len(ok) > 1:
        X = np.array([lu[nm] for nm in nms])
        Xc = X - X.mean()
        g = np.array([float(np.dot(Xc, [x[nm]['lnS'] for nm in nms]) / np.dot(Xc, Xc)) for x in ok])
        m, e = float(g.mean()), se(g)
        pr = lambda a, b: np.array([(x[a]['lnS'] - x[b]['lnS']) / (lu[a] - lu[b]) for x in ok])
        cr = lambda a, b: np.array([(x[a]['lnS'] - lu[a]) - (x[b]['lnS'] - lu[b]) for x in ok])
        out.update(gam=(m, e, m - Z99 * e, m + Z99 * e, float(np.median(g))),
                   pair={(a, b): (float(pr(a, b).mean()), se(pr(a, b))) for a, b in itertools.combinations(nms, 2)},
                   ratio={a: (float(cr(a, nms[-1]).mean()), se(cr(a, nms[-1]))) for a in nms[:-1]})
    return out


def a2t_line(A2T):
    """A2T's line (referee round 137, D3ef-1, -3): the scans of ċ(δ) — per w the states, over every scan (state × arithmetic × order)
    how many have attracting zeros, the nearest ones on either side, every distinct one, ulp(x₁) over uσ*, the slope inside ±0.45 ulp
    over λ, the rms of ċ − λδ —; the RK4 arms — per w and step the exits per arithmetic (`a2t_rk_stats`), at 5A the three-point
    slope γ with its interval, the pairwise exponents and the prefactors' log-ratios —; the verdicts: at a w where some arm holds the
    line in double and long double while another leaves, the exit's existence the integrator's; γ's 99% intervals inside EXP_A2P at every
    arm with three arithmetics at a w whose scans find no attracting zero, the scaling (where there are traps γ is printed, not
    concluded from)."""
    sc = [r for r in A2T if r['kind'] == 'scan']
    rk = [r for r in A2T if r['kind'] == 'rk']
    nms = [nm for nm, _ in DT_A2T]
    parts, traps = [], []
    for r in sc:
        al = [x for row in r['rows'] for x in row['scans']]
        wa = [x for x in al if x['att']]
        neg = [max(a for a in x['att'] if a < 0) for x in wa if any(a < 0 for a in x['att'])]
        pos = [min(a for a in x['att'] if a > 0) for x in wa if any(a > 0 for a in x['att'])]
        dist = sorted({round(a, 3) for x in wa for a in x['att']})
        cnt = {a: sum(1 for x in wa if a in {round(b, 3) for b in x['att']}) for a in dist}
        rg = lambda v, f: f"{min(v):{f}}" if min(v) == max(v) else f"{min(v):{f}} … {max(v):{f}}"
        if wa:
            traps.append(r['wmul'])
        parts.append(f"{r['wmul']:g}A (s = " + ', '.join(f"{row['s']:g}" for row in r['rows']) + ": σ/σ* " + ', '.join(
            f"{row['sig']:.4f}" for row in r['rows']) + f"; λ {rg([x['lam'] for x in al], '.4g')}): attracting zeros in {len(wa)} of "
            f"{len(al)} scans" + ('' if not wa else (
                ", the nearest to 0 at " + ' and '.join(rg(v, '+.3f') for v in (neg, pos) if v) + " ulp, by position (in how many "
                "scans) " + ', '.join(
                    f"{a:+.3f} ({cnt[a]})" for a in dist))) + f", the slope inside ±0.45 ulp {rg([x['slope0'] for x in al], '.2f')}λ"
            f", the rms of ċ − λδ {rg([x['rms'] for x in al], '.2f')} λuσ*")
    ulp = sorted({round(x['ulp'], 3) for r in sc for row in r['rows'] for x in row['scans']})
    arms, held, gam = [], [], []
    for w in dict.fromkeys(r['wmul'] for r in rk):
        rs = [r for r in rk if r['wmul'] == w]
        sub = []
        for r in rs:
            st = a2t_rk_stats(r)
            p = st['per']
            sub.append(f"at {r['frac']:g}·h₀ ({r['h']:.3e}, {r['k']} orders of seed {r['seed']}): exits " + ', '.join(
                str(p[nm]['n']) for nm in nms) + f" of {r['k']} by s = {r['smax']:g}" + ''.join(
                f"; {nm} " + ("—" if p[nm]['n'] == 0 else f"s_e {p[nm]['s']:.4f}, upward {p[nm]['up']}, ln(S/uσ*) {p[nm]['mean']:+.2f}"
                              + ('' if p[nm]['sem'] is None else f" ± {p[nm]['sem']:.2f}") + f" (median {p[nm]['med']:+.2f})")
                for nm in nms) + ('' if 'gam' not in st else (
                    f"; over the {st['n3']} orders where all three leave the three-point slope γ = {st['gam'][0]:.4f} ± {st['gam'][1]:.4f}"
                    f" (99% interval {st['gam'][2]:.4f} … {st['gam'][3]:.4f}, median {st['gam'][4]:.4f}), pairwise " + ', '.join(
                        f"{a}/{b} {m:.4f} ± {e:.4f}" for (a, b), (m, e) in st['pair'].items()) + ", the prefactors' log-ratios " + ', '.join(
                        f"ln(C_{a}/C_{nms[-1]}) {m:+.2f} ± {e:.2f}" for a, (m, e) in st['ratio'].items()))))
            if 'gam' in st and w not in traps:                          # the scaling where the field has no ulp-scale zeros
                gam.append((st['gam'], st['ratio'], r))
        hd = [r for r in rs if all(a2t_rk_stats(r)['per'][nm]['n'] == 0 for nm in nms[1:])]
        lv = [r for r in rs if all(a2t_rk_stats(r)['per'][nm]['n'] > 0 for nm in nms[1:])]
        if hd and lv and w in traps:
            held.append((w, hd, lv))
        arms.append(f"{w:g}A (h₀ {rs[0]['h'] / rs[0]['frac']:.3e}): " + '; '.join(sub))
    okg = bool(gam) and all(EXP_A2P[0] <= g[2] and g[3] <= EXP_A2P[1] for g, _, _ in gam)
    dg = max((max(abs(g[2] - 1), abs(g[3] - 1)) for g, _, _ in gam), default=None)
    cm = max((abs(m) for _, rt, _ in gam for m, _ in rt.values()), default=None)
    return (f"The computed field at the scale of the rounding and three arithmetics at equal steps (A2T; referee round 137): the "
            "exactly symmetric samples with the centre's offset δ as the state and the same right-hand side in " + ', '.join(nms)
            + " (`_DTBumps`: at double and long double bitwise `BumpsCW`'s and `_LDBumps`'"
            ") — (a) ċ(δ) at the long-double line flow's states, on " + f"{NT_A2T} points of |δ| ≤ {ULP_A2T} ulp(x₁) (x₁ = "
            f"{XS[1] - XS[0]:g}, the ±1 samples' coordinate, whose x₁ − δ rounds at half an ulp; ulp(x₁) = " + ', '.join(f"{u_:.3f}" for u_ in ulp)
            + f" uσ* in every arithmetic), in each arithmetic on the code's order and seed {SEED_A2T}'s first two random orders, against "
            "λδ (λ from ±10⁻⁶σ*): " + '; '.join(parts) + " — (b) classical RK4 at the fixed step f·h₀ (h₀ the double's mean DOP853 step "
            f"at rtol {RTOL_A2H:g}, A2P), the same step sequence in every arithmetic and the same orders in each, the exit at |δ| = 10⁻²σ* "
            f"by s = {SMAX_A2T:g} and its effective seed as A2P's (Λ along the line at rtol 1e-12), exits and upward ones per arithmetic ("
            + ', '.join(nms) + "): " + '; '.join(arms) + " — " + '; '.join(
                f"at {w:g}A the computed field's ulp-scale zeros decide: RK4 at " + ' and '.join(f"{r['frac']:g}·h₀" for r in hd)
                + " holds the line in double and long double to s = " + f"{hd[0]['smax']:g}" + ", at " + ' and '.join(
                    f"{r['frac']:g}·h₀" for r in lv) + " it leaves: the exit's existence and its seed are the integrator's, S ∝ u only "
                "through ulp(x₁)" for w, hd, lv in held) + ('; ' if held else '') + (
                "where the scans find no attracting zero (" + ', '.join(f"{w:g}A" for w in dict.fromkeys(r['wmul'] for _, _, r in gam))
                + f") three arithmetics at equal steps put γ's 99% interval inside [{EXP_A2P[0]:g}, {EXP_A2P[1]:g}] at every arm, within "
                f"{dg:.3f} of 1: the seed scales with the unit roundoff, its prefactors a separate O(1) scatter (log-ratios up to "
                f"{cm:.2f} in size)" if okg else "γ not inside the band at every arm without traps" if gam else
                "no arm with three arithmetics where the scans find no attracting zero")
            + ".")


def nh_verdict(lam, tg):
    """Referee round 130 (D3dy-4): normal hyperbolicity at a state on the line from its rates — the normal rate lam = ∂ċ/∂c against the
    tangential ones tg (the (α, σ) block's eigenvalues). Repelling (lam > 0): normally hyperbolic if lam exceeds the largest tangential
    rate, and for every r (lam > k·max tg for k = 0 … r: Hirsch–Pugh–Shub's rate condition) if that rate is negative. Attracting
    (lam < 0): normally hyperbolic if lam lies below the smallest, for every r if that rate is positive. Round 133 (D3eb-4): a zero
    normal rate is not hyperbolic, (False, False). The rates are the frozen Jacobian's eigenvalues: the verdict holds in a metric adapted
    to the frozen linearization (the (α, σ) block is strongly non-normal). Returns (normally hyperbolic, for every r)."""
    if lam == 0:
        return False, False
    if lam > 0:
        return lam > max(tg), lam > max(tg) and max(tg) < 0
    return lam < min(tg), lam < min(tg) and min(tg) > 0


def a2_tangent(ratio, w, sig):
    """Referee round 125 (D3dt-5): at A2E's stored state on the line (c = c*), the exact Jacobian's normal rate ∂ċ/∂c and the
    eigenvalues of its (α, σ) block, the tangential rates (ascending)."""
    L, K, ss, z = a2_state(ratio, w, sig)
    J = np.asarray(RN.critical_jac_nd(L, K, feature_sum([C_A], [ss], [A_A]))(0.0, z))
    return float(J[1, 1]), sorted(float(e.real) for e in np.linalg.eigvals(J[np.ix_([0, 2], [0, 2])]))


def a2_lemma_signs(ratio=0.25, step=0.01, top=1.2, bot=0.15):
    """Referee round 106 (D3da-4): `a2_lemma_scan` on a grid of σ/h (step), and for each w (4A, the least-squares) the largest grid σ/h
    below which (down to bot) the total is positive, the ±1 term positive, and the ±1 term above the centre's. Returns {w: (three)}."""
    rel = tuple(float(x) for x in np.round(np.arange(top, bot - 1e-9, -step), 10))
    sc = a2_lemma_scan(ratio, rel)
    out = {}
    for lab in ('w', 'ls'):
        rows = sorted(sc[lab])                                          # ascending σ/h
        def edge(pred):
            q, nxt = None, None
            for r in rows:
                if not pred(r):
                    nxt = r[0]
                    break
                q = r[0]
            if q is None or nxt is None:                                # round 109 (D3dd-5): the root, by bisection between grid points
                return q
            lo, hi = q, nxt
            for _ in range(40):
                m = (lo + hi) / 2
                lo, hi = (m, hi) if pred(next(x for x in a2_lemma_scan(ratio, (m,))[lab])) else (lo, m)
            return (lo + hi) / 2
        out[lab] = (edge(lambda r: r[3] > 0), edge(lambda r: r[2] > 0), edge(lambda r: r[2] > r[1]))
    return out


def a2_lemma_scan(ratio=0.25, rel=(1.0, 0.8, 0.6, 0.45, 0.3, 0.2), wmul=4.0):
    """Referee round 103 (D3cx-3): Lemma D‴.59a's split on samples exactly symmetric about the centre, at σ = (σ/h)·h for σ/h in rel,
    with w = wmul·A and with the least-squares w on the line (⟨f*φ⟩/⟨φ²⟩): per σ/h (σ/h, the centre term, the ±1 term, the total) for
    each w. The overshoot (the centre) leads only for σ well below h and an O(1) overshoot; at the least-squares w the ±1 term can lead."""
    L, ss, d = law(), ratio * S0, XS[1] - XS[0]
    xs = d * np.arange(-(XS.size // 2), XS.size // 2 + 1)
    f = A_A * np.exp(-xs ** 2 / (2 * ss * ss))
    k = np.rint(xs / d).astype(int)
    out = {}
    for lab in ('w', 'ls'):
        rows = []
        for q in rel:
            sg = q * d
            P = np.exp(-xs ** 2 / (2 * sg * sg))
            wv = wmul * A_A if lab == 'w' else float(np.sum(f * P) / np.sum(P * P))
            r = f - wv * P
            t = L.G * L.m * wv / xs.size * (r * P * (xs * xs / sg ** 4 - 1 / sg ** 2) - wv * (P * xs / sg ** 2) ** 2)
            rows.append((q, float(t[k == 0].sum()), float(t[np.abs(k) == 1].sum()), float(t.sum()), wv / A_A))
        out[lab] = rows
    return out


def mirror_asym(theta, w, c0=0.0):
    """Referee rounds 103, 106 (D3cx-6, D3da-3): the distance of a configuration from its mirror image about c0 (the reflection
    combined with a kernel permutation) — the bottleneck matching: the least, over permutations π, of the largest over
    kernels i of the distance (in c, σ and w, the largest of the three) between kernel i and the mirror image of kernel π(i); 0 exactly
    on the subspace (round 103's nearest image was a lower bound, 0 off it for {a, a, Ra}). It is the distance from the mirror image,
    twice the sup-norm distance to the fixed-point subspace (the midpoint configuration lies on it; referee round 109, D3dd-2)."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import maximum_bipartite_matching
    th, w = np.atleast_2d(np.asarray(theta, float)).reshape(-1, 2), np.asarray(w, float).reshape(-1)
    n = w.size
    if not n:
        return 0.0
    D = np.maximum(np.maximum(np.abs(th[:, None, 0] + th[None, :, 0] - 2 * c0), np.abs(np.abs(th[:, None, 1]) - np.abs(th[None, :, 1]))),
                   np.abs(w[:, None] - w[None, :]))
    vals = np.unique(D)
    lo, hi = 0, vals.size - 1
    while lo < hi:
        m = (lo + hi) // 2
        if (maximum_bipartite_matching(csr_matrix((D <= vals[m]).astype(np.int8)), perm_type='column') >= 0).all():
            hi = m
        else:
            lo = m + 1
    return float(vals[lo])


def mirror_odd(theta, w, c0=0.0, A=1.0, half=1.6, n=3201, xs=None):
    """Referee round 106 (D3da-3): the fitted function's odd part about c0, max_x |f(c0 + x) − f(c0 − x)|/2 over A — at the samples'
    offsets xs ≥ 0 (what the fit sees) or on a fine grid of |x| ≤ half (between the samples too) — gauge-invariant (a coinciding pair's
    parameters may be asymmetric while its function is symmetric)."""
    th, w = np.atleast_2d(np.asarray(theta, float)).reshape(-1, 2), np.asarray(w, float).reshape(-1)
    x = np.linspace(0.0, half, n) if xs is None else np.asarray(xs, float)
    if xs is None:                                                      # round 109 (D3dd-2): resolve every kernel, however narrow
        d_, s_ = np.abs(th[:, 0] - c0), np.abs(th[:, 1])
        x = np.unique(np.concatenate([x] + [np.clip(d_ + k * s_, 0.0, half) for k in np.linspace(-3, 3, 25)]))
    f = lambda y: np.exp(-(y[:, None] - th[None, :, 0]) ** 2 / (2 * th[None, :, 1] ** 2)) @ w
    return float(np.max(np.abs(f(c0 + x) - f(c0 - x))) / 2 / A)


def a2_plateau(ratio):
    """A2's plateau (referee round 83), over ½⟨f*²⟩: lim_{σ→0} min_c J_N — a kernel narrower than the samples' spacing fits at most two
    neighbouring samples, so the limit interpolates the best pair, E₂ = 1 − max_k (f*(x_k)² + f*(x_{k+1})²)/Σf*² — and J_N's σ → 0
    limit at a fixed c on the best sample, E₁ = 1 − max_k f*(x_k)²/Σf*². Returns (E₂, E₁)."""
    q = feature_sum([C_A], [ratio * S0], [A_A]) ** 2
    q = q / q.sum()
    return float(1 - np.max(q[:-1] + q[1:])), float(1 - q.max())


def a2_spike(rows, e2):
    """A2's plateau (referee rounds 80, 83): the outputs on it (E over ½⟨f*²⟩ within PLATEAU_A2 of the two-sample value e2,
    `a2_plateau`), their E, and the first output after them with σ/σ* back above ESC_A2 (the escape; None if none by the last
    output): (first s on it, last s on it, E's range, the escape's s) or None."""
    on = [i for i, x in enumerate(rows) if abs(x['E'] - e2) < PLATEAU_A2]
    if not on:
        return None
    esc = next((rows[i]['s'] for i in range(on[-1] + 1, len(rows)) if rows[i]['sig'] > ESC_A2), None)
    es = [rows[i]['E'] for i in on]
    return rows[on[0]]['s'], rows[on[-1]]['s'], (min(es), max(es)), esc


def a2_lines(A2, A2E=None, A2P=None, A2H=None, A2R=None, A2T=None):
    """A2's table (referee rounds 74, 80, 83): one kernel from (c*, σ*) with a large coefficient — its end, the smallest σ/σ* (at the
    outputs), and the plateau: the outputs with E within PLATEAU_A2 of lim_{σ→0} min_c J_N (`a2_plateau`: the kernel interpolating the
    two samples nearest c*), its E, and the escape —; the plateau's and the one-sample limit's E; A2E's early collapse."""
    if A2 is None or len(A2) == 0:                                      # raw-loaded: an object array
        return []
    out = ['### A2. One kernel from a large coefficient (referee rounds 74, 80, 83)', '',
           f"A's feature at σ*/σ₀ = " + ', '.join(f"{r:g}" for r in RATIOS_A2) + ": one kernel at (c*, σ*) with w = " + ', '.join(
               f"{m:g}" for m in WMUL_A2) + f"·A — E(start) (1 − w/A)²·½⟨f*²⟩, (i)'s condition E(start) < ½‖f*‖² holding for w < 2A only —, "
           f"its block's flow on the critical curve to s = {S_A2[-1]:.0e}: at the end w/A, σ/σ*, (c − c*)/σ* and E over ½⟨f*²⟩; the "
           f"smallest σ/σ* on the way (at the outputs); the plateau — the outputs with E/½⟨f*²⟩ within {PLATEAU_A2:g} of E₂, the kernel "
           "interpolating the two samples nearest c* (below) —, its E over ½⟨f*²⟩, and the first output after it with σ/σ* back above "
           f"{ESC_A2:g} (the escape):", '',
           '| σ*/σ₀ | w(start)/A | E(start)/½⟨f*²⟩ | end: w/A, σ/σ*, (c − c*)/σ*, E/½⟨f*²⟩ | smallest σ/σ* | on the plateau: s, E/½⟨f*²⟩ | '
           'escaped at s | status |', '|---|---|---|---|---|---|---|---|']
    for r in A2:
        e = r['rows'][-1] if r['rows'] else None
        sp = a2_spike(r['rows'], a2_plateau(r['ratio'])[0]) if r['rows'] else None
        out.append(f"| {r['ratio']:g} | {r['wmul']:g} | {r['E0']:g} | " + ('—' if e is None else
                   f"{e['w']:.4f}, {e['sig']:.4g}, {e['dc']:+.1e}, {e['E']:.2e} (s = {e['s']:.0e})") + ' | '
                   + (f"{min(x['sig'] for x in r['rows']):.3g}" if r['rows'] else '—') + ' | '
                   + ('—' if sp is None else f"{sp[0]:.0e} … {sp[1]:.0e}, {sp[2][0]:.3f}–{sp[2][1]:.3f}") + ' | '
                   + ('—' if sp is None else ('not by the end' if sp[3] is None else f"{sp[3]:.0e}")) + f" | {r['status']} |")
    pl = {r_: a2_plateau(r_) for r_ in RATIOS_A2}
    out += ['', "The plateau (referee round 83): a kernel narrower than the samples' spacing fits at most two neighbouring samples, so "
            "lim_{σ→0} min_c J_N interpolates the best pair, E₂ = 1 − max_k (f*(x_k)² + f*(x_{k+1})²)/Σf*² — here the two samples "
            "nearest c* —, and J_N's σ → 0 limit at a fixed c, on the best sample, is E₁ = 1 − max_k f*(x_k)²/Σf*²; over ½⟨f*²⟩ at "
            "σ*/σ₀ = " + ', '.join(f"{r_:g}" for r_ in RATIOS_A2) + ": E₂ " + ', '.join(f"{pl[r_][0]:.6f}" for r_ in RATIOS_A2)
            + ", E₁ " + ', '.join(f"{pl[r_][1]:.4f}" for r_ in RATIOS_A2) + "."]
    inv = {r_: a2_invariance(r_) for r_ in RATIOS_A2}
    out += ['', "The symmetric line (referee rounds 97, 100): on samples exactly symmetric about c* the line c = c* is invariant (Lemma "
            "D‴.59a: the flow commutes with the reflection c ↦ 2c* − c) — the flow's |ċ| on it (w = 3A, 5A; σ/σ* = 1, 0.5, 0.2, 0.1) at most "
            + ', '.join(f"{inv[r_][3]:.1e}" for r_ in RATIOS_A2) + " there, roundoff —; A's samples are symmetric about c* but for "
            f"{inv[RATIOS_A2[0]][0]} at the window's far end, whose f*/A is at most 10^" + ', 10^'.join(f"{{{inv[r_][4]:.1f}}}" for r_ in RATIOS_A2)
            + " (" + ', '.join(f"{inv[r_][1]:.1e}" for r_ in RATIOS_A2) + " in double) at σ*/σ₀ = " + ', '.join(f"{r_:g}" for r_ in RATIOS_A2)
            + ", and in double the samples' and f*'s roundoff (the centre sample lies "
            + f"{a2_ulps(RATIOS_A2[0], 3.0, 1.0)[0]:.0f} ulps of c* above c*) — and where their f* is not negligible the unpaired samples — "
              "drift c on the line: |ċ| there at most " + ', '.join(f"{inv[r_][2]:.1e}" for r_ in RATIOS_A2)
            + ". The problem as posed is a perturbation of the symmetric one, and the symmetric state, for σ well below h, a saddle in c "
            "(Lemma D‴.59a): from c* exactly the flow leaves it too, after a time set by the asymmetry's logarithm over the unstable rate, "
            "and a solver can stay on it by damping the unstable direction (LSODA's steps; below). The samples XS are symmetric about 0 to "
            f"{float(np.max(np.abs(XS + XS[::-1]))):.1e} in double (referee round 106). A2 starts off the line (c* + {A2_OFF:g}σ*): the "
            "generic flow."]
    if A2E is not None and len(A2E):
        more = [r for r in A2E if (r['ratio'], r['wmul']) != (A2E[0]['ratio'], A2E[0]['wmul'])]   # referee round 87: other starts
        A2E = [r for r in A2E if (r['ratio'], r['wmul']) == (A2E[0]['ratio'], A2E[0]['wmul'])]
        e2 = a2_plateau(A2E[0]['ratio'])[0]
        cols = S_A2E[::4]
        out += ['', f"The early collapse (referee round 83): from {A2E[0]['wmul']:g}A at σ*/σ₀ = {A2E[0]['ratio']:g} — c* = {C_A:g} is "
                "a sample — the block's flow (`a2_job`'s, LSODA) over the span (0, S), S = 10¹³ being A2's (a terminal event past s = 10 "
                "leaves its steps alone): (c − c*)/σ*, σ/σ* and E/½⟨f*²⟩ at s ≤ 10, and the first output on the plateau (E within "
                f"{PLATEAU_A2:g} of E₂):", '',
                '| S | ' + ' | '.join(f"s = {s_:.0e}" for s_ in cols) + ' | on the plateau from s |', '|---|' + '---|' * (len(cols) + 1)]
        for r in A2E:
            rs = {x['s']: x for x in r['rows']}
            first = next((x['s'] for x in r['rows'] if abs(x['E'] - e2) < PLATEAU_A2), None)
            out.append(f"| {r['span']:.0e} | " + ' | '.join('—' if s_ not in rs else f"{rs[s_]['dc']:+.1e}, {rs[s_]['sig']:.4f}, "
                                                            f"{rs[s_]['E']:.3f}" for s_ in cols)
                       + ' | ' + ('not by s = 10' if first is None else f"{first:.2g}") + ' |')
        sym = [x for r in A2E for x in r['rows'] if abs(x['dc']) < 1e-12 and 'jcc' in x]
        neg, col = [x for x in sym if x['jcc'] <= 0], [x for x in sym if x['jcc'] > 0]
        if col:
            rg = lambda xs, k: f"{min(x[k] for x in xs):{'.4f' if k == 'sig' else '+.3e'}} … {max(x[k] for x in xs):{'.4f' if k == 'sig' else '+.3e'}}"
            tn = [a2_tangent(A2E[0]['ratio'], x['w'], x['sig'])[1] for x in neg]   # referee round 125 (D3dt-5): not dominated there
            nv = [nh_verdict(x['jcc'], t)[0] for x, t in zip(neg, tn)]   # round 130 (D3dy-4): the verdict from the printed rates
            out += ['', "On the symmetric state (\\|c − c*\\| below 10⁻¹²σ*) the flow's ∂ċ/∂c is " + (
                f"{rg(neg, 'jcc')} while the width falls (σ/σ* {rg(neg, 'sig')}, {len(neg)} outputs: stable in c; the (α, σ) block's "
                f"eigenvalues there {min(t[0] for t in tn):.0f} … {max(t[0] for t in tn):.0f} and {min(t[1] for t in tn):+.2g} … "
                f"{max(t[1] for t in tn):+.2g}, " + ("c contracting faster than both at every output: normally hyperbolic" if all(nv) else
                                                    "the width contracting faster than c (∂ċ/∂c above the smallest at every output): not "
                                                    "normally hyperbolic" if not any(nv) else
                                                    f"normally hyperbolic at {sum(nv)} of {len(nv)} outputs")
                + ", referee rounds 125, 130) and " if neg else '')
                + f"{rg(col, 'jcc')} at σ/σ* {rg(col, 'sig')} ({len(col)} outputs): a saddle in c, an offset growing ×e per "
                f"{1 / max(x['jcc'] for x in col):.1e} … {1 / min(x['jcc'] for x in col):.1e} of s, left once the seed — on A's samples the data's drift, on exactly symmetric ones the sums' rounding (below) — has grown."]
            lam = max(x['jcc'] for x in col)
            for r in A2E:
                st = [x for x in r.get('steps', []) if abs(x[3]) < 1e-6]   # on the symmetric state
                if not st:
                    continue
                ms = [f"{nm} at hλ {min(x[1] * lam for x in st if x[2] == m_):.2g}–{max(x[1] * lam for x in st if x[2] == m_):.2g}"
                      + ('' if m_ == 1 else f" from s = {min(x[0] for x in st if x[2] == m_):.3g}") for m_, nm in ((1, 'Adams'), (2, 'BDF'))
                      if any(x[2] == m_ for x in st)]
                bd = [(x, y) for x, y in zip(st[:-1], st[1:]) if y[2] == 2 and len(y) > 4 and abs(x[3]) > 1e-15
                      and abs(y[3]) > A2E_ROUND]                       # BDF steps by order; round 97: not those into roundoff
                bo = ''.join(f"; its BDF order-{k_} steps at hλ {min(y[1] * lam for _, y in v):.3g}–{max(y[1] * lam for _, y in v):.3g} "
                             f"multiply the offset by {min(y[3] / x[3] for x, y in v):+.3g} … {max(y[3] / x[3] for x, y in v):+.3g}"
                             + (f" (1/(1 − hλ): {min(1 / (1 - y[1] * lam) for _, y in v):+.3g} … {max(1 / (1 - y[1] * lam) for _, y in v):+.3g})"
                                if k_ == 1 else '')
                             for k_, v in sorted({k_: [z for z in bd if z[1][4] == k_] for k_ in {y[4] for _, y in bd}}.items()))
                nr_ = sum(1 for x, y in zip(st[:-1], st[1:]) if y[2] == 2 and len(y) > 4 and abs(x[3]) > 1e-15 and abs(y[3]) <= A2E_ROUND)
                bo += (f"; {nr_} further BDF step{'s' if nr_ > 1 else ''} left out, where the offset falls below {A2E_ROUND:g}σ* — to "
                       "the scale of the drift's fixed point −d/λ (below; referee round 100)" if nr_ else '')
                gr = [x for x in st if x[2] == 1 and 1e-12 < abs(x[3])]
                rate = (np.log(abs(gr[-1][3]) / abs(gr[0][3])) / (gr[-1][0] - gr[0][0])) if len(gr) > 1 and gr[-1][0] > gr[0][0] else None
                out.append(f"LSODA's steps on the symmetric state (s in ({S_A2E_STEPS[0]:g}, {S_A2E_STEPS[1]:g}], the offset below 10⁻⁶σ*; "
                           f"λ = {lam:.4g}), S = {r['span']:.0e}: " + ', '.join(ms) + f"; the offset's largest {max(abs(x[3]) for x in st):.1e}σ*, "
                           f"at the last such step {st[-1][3]:+.1e}σ*" + ('' if rate is None else f"; over the Adams steps above 10⁻¹²σ* it "
                           f"grows at {rate:.4g} per unit s") + bo + '.')
            xl = max(col, key=lambda x: x['jcc'])                       # referee round 100 (D3cu-2): the lemma's split there
            pt, tot, J, hs = a2_lemma(A2E[0]['ratio'], xl['w'], xl['sig'])
            out.append(f"Lemma D‴.59a's split of ∂ċ/∂c = Gm·w·⟨r∂²_cφ − w(∂_cφ)²⟩ at its largest output (s = {xl['s']:.2g}, σ/σ* {xl['sig']:.4f}, "
                       f"the samples' spacing {hs:.2f}σ): the centre sample {pt['0']:+.5g} — (Gm w/(Nσ²))(w − f*(c*)), the kernel overshooting "
                       f"the sample it sits on —, the samples at ±1 {pt['1']:+.5g}, at ±2 {pt['2']:+.4g}, the rest {pt['rest']:+.2g}; their sum "
                       f"{tot:.10g} against the exact Jacobian's {J:.10g}.")
            out.append("(BDF-k damps a positive real eigenvalue only for hλ beyond " + ', '.join(f"{bdf_bound(k_):.4g}" for k_ in range(1, 6))
                       + " (k = 1 … 5; referee round 87): below that an implicit step amplifies the unstable offset; on the Adams steps it "
                       "grows at λ, the flow's own instability.)")
        for r in more:                                                  # referee round 87 (D3ch-1): the saddle holding the exact start
            mx = max(abs(x['dc']) for x in r['rows'])
            e2m = a2_plateau(r['ratio'])[0]
            lam_at = lambda s_: next((x['jcc'] for x in r['rows'] if x['s'] >= s_), r['rows'][-1]['jcc'])
            bq = [(x[0], x[1] * lam_at(x[0]), x[4]) for x in r.get('steps', []) if x[2] == 2 and len(x) > 4]
            first = next((x['s'] for x in r['rows'] if abs(x['E'] - e2m) < PLATEAU_A2), None)
            out.append(f"From {r['wmul']:g}A at σ*/σ₀ = {r['ratio']:g} with the centre at c* exactly (A2's start before referee round 87), "
                       f"over the span {r['span']:.0e}: |c − c*|/σ* at most {mx:.1e} at the outputs to s = {max(x['s'] for x in r['rows']):g}, "
                       f"σ/σ* down to {min(x['sig'] for x in r['rows']):.4f}, " + ('not on the plateau by then' if first is None else
                                                                                  f"on the plateau from s = {first:.2g}")
                       + (" — held on the symmetric state by LSODA's damping" if abs(r['rows'][-1]['dc']) < 1e-12 else " — leaving the symmetric state")
                       + f" (A2 starts off it, at c* + {A2_OFF:g}σ*; referee round 106: the clause from the data)."
                       + ((lambda st: f" Its steps to s = {r['steps_to']:g} (referee round 100): Adams to s = "
                           f"{max([x[0] for x in st if x[2] == 1] or [0]):.3g}, then BDF at hλ "
                           + ', '.join(f"{min(hl for _, hl, o in bq if o == k_):.3g}–{max(hl for _, hl, o in bq if o == k_):.3g} (order {k_})"
                                       for k_ in sorted({o for _, _, o in bq}))
                           + f", λ the outputs' ∂ċ/∂c; the offset at most {max(abs(x[3]) for x in st):.1e}σ*, at the last step {st[-1][3]:+.2e}σ*.")(
                           r['steps']) if r.get('steps') and r.get('steps_to', 0) > S_A2E_STEPS[1] else ''))
        held = [(r, r['rows'][-1]) for r in A2E + more if r['rows'] and abs(r['rows'][-1]['dc']) < 1e-12]   # referee round 100 (D3cu-1)
        if held:
            out.append("The held runs at s = 10 (referee round 100): " + '; '.join(
                (lambda fp, tg, nh: f"{r['wmul']:g}A over the span {r['span']:.0e}: the line's drift ċ at c = c* {fp[0]:+.2e}σ* per unit s, "
                                    f"λ = {fp[1]:.4g} (the (α, σ) block's eigenvalues {tg[1][1]:.2g} and {tg[1][0]:.4g}: " + (
                                        "λ above both, normally hyperbolic" + (", for every r as both are negative" if nh[1] else '')
                                        if nh[0] else "a tangential rate above λ: not normally hyperbolic")
                                    + f"; referee rounds 125, 130), the implicit step's fixed point −d/λ = {fp[2]:+.3e}σ* against the held "
                                    f"offset {fp[3]:+.3e}σ*")(
                    *(lambda fp, tg: (fp, tg, nh_verdict(fp[1], tg[1])))(a2_fixed_point(r['ratio'], x['w'], x['sig'], x['dc']),
                                                                       a2_tangent(r['ratio'], x['w'], x['sig']))) for r, x in held)
                + " — the offset sits within an ulp of the zero of ċ (−d/λ with d summed in double carries the double sum's error in "
                "d), the computed problem's own near-symmetric state (referee rounds 103, 106, 109: "
                + '; '.join((lambda u, ue: f"{r['wmul']:g}A: held {x['dc'] / u[2]:.0f} ulps of c* above c*"
                                           + ('' if not r.get('steps') or r.get('steps_to', 0) < 10 else " (its last 10 steps to s = "
                                              f"{r['steps_to']:g} at " + ', '.join(str(k_) for k_ in sorted({int(round(y[3] / u[2])) for y in r['steps'][-10:]}))
                                              + " ulps)")
                                           + ", the zero of ċ at " + ('—' if ue is None else f"{ue:.2f}") + " ulps in 50-digit arithmetic on the "
                                           "double data (" + ('—' if u[1] is None else f"{u[1]:.1f}") + f" in double), the centre sample {u[0]:.0f} "
                                           "ulps above c*")(
                    a2_ulps(r['ratio'], x['w'], x['sig']), a2_ulps_exact(r['ratio'], x['w'], x['sig'])) for r, x in held)
                + "): a repelling state the flow leaves, which LSODA's steps (BDF beyond its bounds) make attracting — the location the "
                "problem's, the stability LSODA's step sequence's (referee round 106: scipy's Radau and BDF need not hold it, below).")
        kzs = {r['wmul']: a2_ulps_exact(r['ratio'], x['w'], x['sig']) for r, x in held} if held else {}   # round 109: the 50-digit zero
        ex = a2_exits(0.25, (3.0, 5.0), kzs)
        sg = lambda t: '—' if t[0] is None else f"{t[0]:.4f} ({'+' if t[1] > 0 else '−'})"
        out.append("The flow from c* exactly (referee rounds 103, 106, 109): |c − c*| reaches 10⁻²σ* (the exit's direction in brackets) at s = "
                   + '; '.join(f"{wm:g}A: on A's samples {sg(r['a'][0])}, {sg(r['a'][1])} (DOP853, rtol 1e-10, 1e-12), "
                               + ', '.join(f"{sg(t)} ({m})" if t[0] is not None else f"held to s = 1 ({m})" for m, t in r['imp'])
                               + f" (scipy, rtol 1e-10 unless stated), {sg(r['summ']['fsum'][0])} with the samples' sums exactly rounded, "
                                 f"{sg(r['summ']['rev'][0])} summed in reversed order; on exactly symmetric samples from their centre "
                                 f"{sg(r['sym'][0])}, {sg(r['sym'][1])} (DOP853), "
                               + ("no exit by s = 1" if r['summ']['fsum'][1][0] is None else sg(r['summ']['fsum'][1])) + " exactly rounded, "
                               + f"{sg(r['summ']['rev'][1])} reversed"
                               + ('' if r['symoff'] is None else f"; on them started {r['kz']:.2f} ulps of c* below their centre — the "
                                  f"posed problem's c* below its zero of ċ — {sg(r['symoff'])}") for wm, r in sorted(ex.items())))
        ok_ = lambda t: t is not None and t[0] is not None
        rel = lambda a, b: abs(a[0] / b[0] - 1)
        sgl = lambda t, k=5: '—' if t[0] is None else f"{t[0]:#.{k}g} ({'+' if t[1] > 0 else '−'})"   # round 125 (D3dt-3): 5 digits
        lda = lambda a, b: (f"{sgl(((a[0] + b[0]) / 2, a[1]))} at rtol 1e-10 and 1e-12 (apart by {rel(a, b):.1e} relatively)"
                            if ok_(a) and ok_(b) and a[1] == b[1] else f"{sgl(a)}, {sgl(b)} (rtol 1e-10, 1e-12)")
        out[-1] += (". In long double with the offset c − c* as the state (`a2_exit_ld`; referee round 112): " + '; '.join(
            f"{wm:g}A: on A's samples {lda(*r['ld']['a'])}, on the symmetric samples "
            f"{sgl(r['ld']['sym'][0], 3)}, {sgl(r['ld']['sym'][1], 3)}" + ('' if r['ld']['seed'] is None else
            f"; the data's drift carried along the route (`a2_seed`: S = ∫e^{{−Λ}}ċ|_{{c*}}ds to the exit, Λ = ∫∂ċ/∂c|_{{c*}}) S = "
            f"{r['ld']['seed']:.3e} ({r['ld']['seed'] / r['ld'].get('ss', S0 / 4):.2e}σ*)" + ('' if r['kz'] is None else
            f", {r['ld']['seed'] / (-r['kz'] * r['ld']['ulp']):.2f} times c*'s distance below the zero of ċ ({r['kz']:.2f} ulps)")
            + ('' if r['ld'].get('Lam') is None else f", e^{{Λ(s_e)}}|S| {math.exp(r['ld']['Lam']) * abs(r['ld']['seed']) / (0.01 * r['ld'].get('ss', S0 / 4)):.3f} "
               "times the exit's level (the linear response)") + f", and the symmetric problem started at c* + S leaving at {sgl(r['ld']['symS'])}"
            + ('' if not (ok_(r['ld']['symS']) and ok_(r['ld']['a'][1])) else
               (", A's time within 3·10⁻⁵ (a 10⁻³ change of S; the floor of a few 10⁻⁶ is the long double's own rounding of ċ, which "
                "no tolerance moves — referee round 125)" if rel(r['ld']['symS'], r['ld']['a'][1]) < 3e-5 else
                f", A's time to {rel(r['ld']['symS'], r['ld']['a'][1]):.1e}")))
            for wm, r in sorted(ex.items())))
        A_ = lambda r: r['a'] + [m[1] for m in r['imp'] if m[1][0] is not None] + [r['summ']['fsum'][0], r['summ']['rev'][0]]
        data = all(all(t[0] is not None and t[1] < 0 for t in A_(r) + r['ld']['a'])     # round 112: decided in long double
                   and all(abs(t[0] / r['ld']['a'][1][0] - 1) < 0.02 for t in r['a'] + [r['summ']['fsum'][0], r['summ']['rev'][0]])
                   and rel(r['ld']['a'][0], r['ld']['a'][1]) < 1e-4
                   and ok_(r['ld']['symS']) and r['ld']['symS'][1] == r['ld']['a'][1][1] and rel(r['ld']['symS'], r['ld']['a'][1]) < 1e-4
                   and r['symoff'] is not None and r['symoff'][0] is not None and r['symoff'][1] < 0 for r in ex.values())
        flip = lambda r: r['sym'][0][1] != r['sym'][1][1] or r['ld']['sym'][0][1] != r['ld']['sym'][1][1]
        flipa = lambda r: all(ok_(r['sym'][i]) and ok_(r['ld']['sym'][i]) and r['sym'][i][1] != r['ld']['sym'][i][1] for i in (0, 1))
        dr = lambda ts: ', '.join('+' if t[1] > 0 else '−' for t in ts)   # round 125 (D3dt-6): both tolerances' directions
        dev = max(abs(t[0] / r['ld']['a'][1][0] - 1) for r in ex.values() for t in r['a'] + [r['summ']['fsum'][0], r['summ']['rev'][0]]
                  if t[0] is not None) if all(ok_(r['ld']['a'][1]) for r in ex.values()) else None
        lags = [x for r in ex.values() for x in (r['ld'].get('lag') or []) if x is not None]
        late = [r['ld']['sym'][i][0] / r['sym'][i][0] - 1 for r in ex.values() for i in (0, 1)
                if ok_(r['ld']['sym'][i]) and ok_(r['sym'][i])]
        arith = all(r['summ']['fsum'][1][0] is None and r['sym'][0][0] is not None and r['summ']['rev'][1][0] is not None
                    and r['summ']['rev'][1][1] == -r['sym'][0][1] for r in ex.values()) and len(late) == 2 * len(ex) and min(late) > 0.05
        out[-1] += ((" — on A's samples every exit downward, in double its time within "
                     + f"{math.ceil(dev * 1000 - 1e-9) / 10:.1f}% of the long-double one whatever the "
                     "tolerance and the sums' order or rounding (DOP853), in long double the same at both tolerances to 10⁻⁴, and the "
                     "symmetric problem started at the data's seed leaving with A's direction at A's time to 10⁻⁴: the exit is the data's"
                     if data else '')
                    + ("; on exactly symmetric samples from their centre exact sums keep c at the centre (no exit by s = 1) and reversed "
                       "sums mirror the run — identities of the symmetric data, not tests —, "
                       + ("the integrator's tolerance flips the direction (" + ', '.join(f"{wm:g}A" + (" in long double" if r['sym'][0][1] == r['sym'][1][1] else '')
                                                                                for wm, r in sorted(ex.items()) if flip(r)) + "), " if any(flip(r) for r in ex.values()) else '')
                       + ("in the code's order of the sums the arithmetic flips it at both tolerances (" + '; '.join(
                           f"{wm:g}A: {dr(r['sym'])} in double and {dr(r['ld']['sym'])} in long double at rtol 1e-10, 1e-12"
                           for wm, r in sorted(ex.items()) if flipa(r)) + "), " if any(flipa(r) for r in ex.values()) else '')
                       + f"and in long double the exit comes {math.floor(min(late) * 100):.0f}–{math.ceil(max(late) * 100):.0f}% "
                       "later" + ('' if len(lags) != 2 * len(ex) else
                                  f" — by Λ = ∫λ ds {min(lags):.2f} … {max(lags):.2f} along the route (`a2_sym_lag`, from the double exit to the "
                                  "long-double one at each tolerance: the two seeds' ratio e^Λ) —") + ", each in the code's order one draw (A2P, "
                       "below): its seed is the arithmetic's" if arith else '') + ".")
        if A2P is not None and len(A2P):                                # referee round 125 (D3dt-1): the draws' distribution
            out.append(a2p_line(A2P, {(wm, q): x for wm, r in ex.items() for q, x in zip((1e-10, 1e-12), r['ld'].get('lag') or [])}, A2R))
            if A2H is not None and len(A2H):                            # referee round 130 (D3dy-2): the steps
                out.append(a2h_line(A2H, A2P))
                if A2R is not None and len(A2R):                        # referee round 133 (D3eb-1, -2): the replications
                    out.append(a2r_line(A2R, A2P, A2H))
                if A2T is not None and len(A2T):                        # referee round 137 (D3ef-1, -3): the traps, three arithmetics
                    out.append(a2t_line(A2T))
        sc, sg_ = a2_lemma_scan(), a2_lemma_signs()
        out.append("Lemma D‴.59a's terms on exactly symmetric samples at σ*/σ₀ = 0.25 against σ/h (h the spacing; the centre, the ±1 samples, "
                   "the total; referee round 103): " + '; '.join(
                       ("w = 4A" if lab == 'w' else "the least-squares w") + ": " + ', '.join(
                           f"σ/h {q:g}: {a:+.4g}, {b:+.4g}, {t:+.4g}" + (f" (w/A {wa:.4f})" if lab == 'ls' else '') for q, a, b, t, wa in sc[lab])
                       for lab in ('w', 'ls')) + ". On a grid of σ/h step 0.01 from 1.2 down to 0.15 (`a2_lemma_signs`; referee round 106) "
                   "the σ/h below which the total, the ±1 term, and the ±1 term above the centre's hold positive (each a root, by bisection "
                   "between the grid's points; referee round 109): " + '; '.join(
                       ("w = 4A" if lab == 'w' else "the least-squares w") + ": " + ', '.join('—' if v is None else f"{v:.4f}" for v in sg_[lab])
                       for lab in ('w', 'ls')) + " — the scan's states, not A2E's.")
    return out + ['']


def part_B():
    """The reduced energy on a (c, σ) grid around the exact fit for each ratio: its grid-local minima, and the gradient's check."""
    out = []
    cs_, A = 0.0, 1.0
    for ratio in RATIOS_A:
        ss = ratio * S0
        cg = np.linspace(-6 * ss, 6 * ss, 241)
        sg = np.geomspace(ss / 8, 8 * ss, 241)
        J = MW.reduced_energy(cg[:, None], sg[None, :], A, cs_, ss)
        inner = J[1:-1, 1:-1]
        mins = np.argwhere((inner < J[:-2, 1:-1]) & (inner < J[2:, 1:-1]) & (inner < J[1:-1, :-2]) & (inner < J[1:-1, 2:]))
        locs = [(float(cg[i + 1]), float(sg[j + 1])) for i, j in mins]
        rng = np.random.default_rng(SEED_C + int(100 * ratio))
        gerr = 0.0
        q = lambda c, s: MW.overlap(c, s, cs_, ss) ** 2 / (np.sqrt(np.pi) * s)        # J = ½A²(π^{1/2}σ* − q): differenced without
        for _ in range(20):                                                            # the constant, which would swamp a tail's gradient
            c, s = rng.normal(0, 2 * ss), ss * np.exp(rng.normal(0, 0.7))
            g = MW.one_kernel_grad(c, s, A, cs_, ss)
            h = 1e-6 * ss
            fd = -0.5 * A ** 2 * np.array([(q(c + h, s) - q(c - h, s)) / (2 * h), (q(c, s + h) - q(c, s - h)) / (2 * h)])
            gerr = max(gerr, float(np.abs(fd - g).max() / max(np.abs(g).max(), 1e-300)))
        out.append(dict(ratio=ratio, minima=locs, J_min=float(J.min()), grad_err=gerr))
    return out


def targets_C(rng, n, dmin):
    out = []
    while len(out) < n:
        cs = rng.uniform(-0.7, 0.7, K_C)
        if min(abs(cs[i] - cs[j]) for i in range(K_C) for j in range(i)) < dmin:
            continue
        ss = np.exp(rng.uniform(np.log(WID_C[0]), np.log(WID_C[1]), K_C))
        As = rng.choice([-1.0, 1.0], K_C) * rng.uniform(0.4, 1.0, K_C)
        out.append((cs, ss, As))
    return out


def c_job(args):
    q, dmin, (cs, ss, As) = args
    t0 = time.time()
    f = feature_sum(cs, ss, As)
    Rm, Rf = run_paths(f)
    return dict(q=q, dmin=dmin, cs=cs.tolist(), ss=ss.tolist(), As=As.tolist(), moving=match(Rm, cs, ss, As), fixed=match(Rf, cs, ss, As),
                secs=time.time() - t0)


def part_C():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    jobs = [(q, d, t) for d in DELTAS_C for q, t in enumerate(targets_C(rng, N_C, d))]
    with Pool(N_WORKERS) as pool:
        return pool.map(c_job, jobs, chunksize=2)


def d_job(d):
    cs, ss, As = [0.0, d], [0.3, 0.06], [0.8, 0.5]
    f = feature_sum(cs, ss, As)
    Rm, Rf = run_paths(f)
    return dict(d=d, moving=match(Rm, cs, ss, As), fixed=match(Rf, cs, ss, As))


def part_D():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(d_job, D_D, chunksize=1)


def blocks(R):
    """The relaxed blocks along a path: (kernels, the residual's rms) after each event, with the block."""
    out = []
    for tau, c, w, rms in R['hist'][1:]:
        act = np.flatnonzero(w != 0)
        out.append(dict(k=int(act.size), rms=float(rms), c=np.atleast_2d(c)[act], w=w[act]))
    return out


def e_job(args):
    q, s, (cs, ss, As) = args
    f = feature_sum(cs, ss, As)
    xi = s * np.random.default_rng(SEED_C + 1000 * q + int(1e4 * s)).standard_normal(XS.size)
    Rm, _ = run_paths(f + xi, max_entries=K_C + 2, fixed=False)
    out = dict(q=q, s=s, oracle=oracle_error(f, xi, cs, ss, As))
    for name, R, moving in (('moving', Rm, True),):
        bl = [b for b in blocks(R) if b['k'] > 0]
        stop = next((b for b in bl if b['rms'] <= (1 + DELTA_E) * s), None)
        if stop is None:
            out[name] = dict(stopped=False, oob=bool(R.get('out_of_budget', False)), rms=bl[-1]['rms'] if bl else np.nan)
            continue
        th = stop['c']
        P = (np.exp(-(XS[:, None] - th[:, 0][None, :]) ** 2 / (2 * th[:, 1][None, :] ** 2)) if moving else
             np.exp(-(XS[:, None] - th[:, 0][None, :]) ** 2 / (2 * S0 ** 2)))
        out[name] = dict(stopped=True, k=stop['k'], err=float(np.sqrt(np.mean((P @ stop['w'] - f) ** 2))))
    return out


def oracle_error(f, xi, cs, ss, As):
    """The oracle fit: the K-kernel moving-width least-squares fit of f + ξ started at the truth (the right support); its error against f
    (rms on the samples) — per draw, against the stopped fit's (referee round 57)."""
    from scipy.optimize import least_squares
    k = len(As)
    K = MW.BumpsCW(XS, S0)
    y = f + xi
    res = lambda z: y - K.synth(z[k:].reshape(k, 2), z[:k])
    sol = least_squares(res, np.r_[As, np.c_[cs, ss].ravel()], method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    return float(np.sqrt(np.mean((K.synth(sol.x[k:].reshape(k, 2), sol.x[:k]) - f) ** 2)))


def part_E():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    targets_C(rng, N_C, DELTAS_C[0])                                   # C's stream: the Δ = 0.4 targets first, then Δ = 0.6's
    tg = targets_C(rng, N_C, DELTAS_C[1])[:N_E]
    with Pool(N_WORKERS) as pool:
        return pool.map(e_job, [(q, s, t) for s in LEVELS_E for q, t in enumerate(tg)], chunksize=2)


def f_job(args):
    q, (cs, ss, As) = args
    t0 = time.time()
    L = law()
    K = MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    R = RN.sliding_iss_nd(L, K, f, theta0(), np.ones(C0.size), max_rhs=MAX_RHS_MOVING)
    t_last = R['events'][-1][0]
    taus = np.r_[0.0, np.array(TAUS_F) * t_last]
    path = sorted(int(j) for j in np.flatnonzero(R['act']))
    pp = closest_pair(R['c'][path], R['w'][path])
    snaps = []
    for eps in EPS_F:
        _, W, C, Lam, al = RN.loop_nd(L, K, f, theta0(), eps, taus, rtol=1e-9, atol=1e-12)
        for t in range(1, taus.size):
            big = np.abs(W[t]) > 1e-2 * np.abs(As).min()
            fj = np.flatnonzero(big)
            cp = closest_pair(C[t][fj], W[t][fj])
            snaps.append(dict(eps=eps, times=float(TAUS_F[t - 1]), formed=sorted(int(j) for j in fj),
                              dtheta=float(np.abs(C[t][path] - R['c'][path]).max()), dw=float(np.abs(W[t][path] - R['w'][path]).max()),
                              rest=float(np.abs(np.delete(W[t], path)).max()),
                              pair=None if cp is None else dict(d=cp[0], dc=float(abs(C[t][fj[cp[1]], 0] - C[t][fj[cp[2]], 0])),
                                                                ds=float(abs(C[t][fj[cp[1]], 1] - C[t][fj[cp[2]], 1])),
                                                                w=(float(W[t][fj[cp[1]]]), float(W[t][fj[cp[2]]]))),
                              rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2)))))
    return dict(q=q, path=path, path_pair=None if pp is None else pp[0], path_rms=float(R['hist'][-1][3]), snaps=snaps,
                formed=snaps[0]['formed'], dtheta=snaps[0]['dtheta'], dw=snaps[0]['dw'], rest=snaps[0]['rest'], secs=time.time() - t0)


def part_F():
    from multiprocessing import Pool
    rng = np.random.default_rng(SEED_C)
    targets_C(rng, N_C, DELTAS_C[0])
    tg = targets_C(rng, N_C, DELTAS_C[1])[:N_F]
    with Pool(N_WORKERS) as pool:
        return pool.map(f_job, list(enumerate(tg)), chunksize=1)


def dipole_law(L, K, f, th, w):
    """The dipole's law at a block's state (θ, w) (referee round 57): its pair — the largest |w| and the closest kernel of the other sign —
    (a, b), v = |w_a − w_b|/2 (the
    half-difference, moved at the mobility μ₀|w|, μ₀ = 4Gc/A on the free law's exponential branch), the moment p = w_aθ_a + w_bθ_b −
    (w_a + w_b)θ̄ about the midpoint θ̄ (its σ-share |p_σ|/|p|), and the optimal dipole fit — the other kernels, Wφ(θ̄) and p·∇_θφ(θ̄)
    by least squares from the state (`closure_fit`) —: E_inf = ½⟨r²⟩ there. Returns dict(v, E, E_inf, rms_inf, p, mu0) or None."""
    th, w = np.atleast_2d(np.asarray(th, float)), np.asarray(w, float)
    a = int(np.argmax(np.abs(w)))
    opp = [j for j in range(w.size) if w[j] * w[a] < 0]
    if not opp:
        return None
    dist = lambda j: (abs(th[a, 0] - th[j, 0]) / (0.5 * (abs(th[a, 1]) + abs(th[j, 1]))) + abs(np.log(abs(th[a, 1]) / abs(th[j, 1]))))
    b = min(opp, key=dist)
    oth = [i for i in range(w.size) if i not in (a, b)]
    tb = 0.5 * (th[a] + th[b])
    p0 = w[a] * (th[a] - tb) + w[b] * (th[b] - tb)
    cf = closure_fit(K, f, th[oth], w[oth], w[a] + w[b], tb, p0)
    E_inf = float(np.mean(cf['r'] ** 2) / 2)
    return dict(v=float(abs(w[a] - w[b]) / 2), E=float(np.mean((f - K.synth(th, w)) ** 2) / 2), E_inf=E_inf,
                rms_inf=float(np.sqrt(2 * E_inf)), p=cf['p'].tolist(), mu0=float(4 * L.G * L.c / L.A),
                dc=float(abs(th[a, 0] - th[b, 0])), ds=float(abs(th[a, 1] - th[b, 1])))


def closure_fit(K, f, th_o, w_o, W0, tb0, p0, max_nfev=20000):
    """The dipole closure's least-squares fit (referee rounds 57, 74): the other kernels (θ, w), Wφ(θ̄) and p·∇_θφ(θ̄) by
    Levenberg–Marquardt from the given start. Returns dict(r — the fit's residual —, tb, p, th, w — the other kernels, their widths as
    |σ| (referee round 80: φ is even in σ) —, W)."""
    from scipy.optimize import least_squares
    th_o, w_o = np.asarray(th_o, float).reshape(-1, 2), np.asarray(w_o, float).reshape(-1)
    m = w_o.size

    def model(z):
        g = K.synth(z[m:3 * m].reshape(m, 2), z[:m]) if m else 0.0
        t_ = z[3 * m + 1:3 * m + 3]
        return g + z[3 * m] * K.phi(t_[None, :])[:, 0] + K.dphi(t_[None, :])[:, 0, :] @ z[3 * m + 3:]
    z0 = np.r_[w_o, th_o.ravel(), W0, tb0, p0]
    sol = least_squares(lambda z: (f - model(z)) / np.sqrt(f.size), z0, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15,
                        max_nfev=max_nfev)
    z = sol.x
    th = z[m:3 * m].reshape(m, 2).copy()
    th[:, 1] = np.abs(th[:, 1])
    return dict(r=f - model(z), tb=z[3 * m + 1:3 * m + 3], p=z[3 * m + 3:], th=th, w=z[:m], W=float(z[3 * m]))


def g_job(args):
    """G: a path stopped by its budget (C) — C's path again, the stalled block's state when the budget ran out (its largest |w|, the
    opposite-sign pair's distance, its residual) and the polishes of that block; the block from the entry that started it (the path to
    that many entries, the entrant at its seed): its flow (the fast flow on the critical curve, LSODA with the exact Jacobian) to
    s = 1e11, the largest |w|, the closest opposite-sign pair (|Δc|/σ + |Δlog σ|) and the dipole law (`dipole_law`) at S_G; and the
    real loop at the masses EPS_G at the times TAUS_G: its formed kernels, the closest pair's distance, coefficients and midpoint's offset
    from the nearest feature (in its width), the coinciding pairs, the largest |w| over the largest amplitude, the residual."""
    dmin, q, (cs, ss, As), n_ent = args
    from scipy.integrate import solve_ivp
    t0 = time.time()
    L = law()
    K = MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    Rc, _ = run_paths(f, fixed=False)
    stalled = stalled_summary(Rc, f, K)
    R = RN.sliding_iss_nd(L, K, f, theta0(), np.ones(C0.size), max_entries=n_ent, exit_kkt=True)
    idx = np.flatnonzero(R['act'])
    al = np.sign(R['w'][idx]) * RT._alpha_of(L, np.abs(R['w'][idx]), RT.coef_crit)
    z0 = np.r_[al, R['c'][idx].ravel()]
    k = idx.size
    sol = solve_ivp(RN.critical_rhs_nd(L, K, f), (0.0, S_G[-1]), z0, method='LSODA', jac=RN.critical_jac_nd(L, K, f), rtol=1e-10,
                    atol=1e-14, t_eval=S_G)
    flow = []
    for i in range(sol.y.shape[1]):
        w = RT.critical_w(L, sol.y[:k, i])
        th = sol.y[k:, i].reshape(k, 2)
        sep = min((abs(th[a, 0] - th[b, 0]) / (0.5 * (abs(th[a, 1]) + abs(th[b, 1]))) + abs(np.log(abs(th[a, 1]) / abs(th[b, 1])))
                   for a in range(k) for b in range(a + 1, k) if w[a] * w[b] < 0), default=None)
        flow.append(dict(s=float(sol.t[i]), wmax=float(np.abs(w).max()), sep=None if sep is None else float(sep),
                         rms=float(np.sqrt(np.mean((f - K.synth(th, w)) ** 2))), theta=th.tolist(), w=w.tolist()))
    law_ = [dipole_law(L, K, f, x['theta'], x['w']) for x in flow[-2:]] if sol.status == 0 else []
    fine, loops = {}, []
    tg = np.unique(np.r_[TAUS_G, TAUS_FINE])                            # one run per mass: the extra outputs leave the steps alone
    for eps in EPS_G + EPS_GX:
        _, W, C, Lam, _ = RN.loop_nd(L, K, f, theta0(), eps, np.r_[0.0, tg], rtol=1e-9, atol=1e-12)
        sts = []
        for t, tau in enumerate(tg, start=1):
            big = np.flatnonzero(np.abs(W[t]) > 1e-2 * np.abs(As).min())
            sts.append((float(tau), C[t][big], W[t][big]))
            if eps in EPS_G and tau in TAUS_G:
                loops.append(dict(eps=eps, tau=float(tau), Lam=float(Lam), theta=C[t][big].tolist(), w=W[t][big].tolist(),
                                  wmax=float(np.abs(W[t]).max()), rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2)))))
        fine[eps] = fine_phases(sts, n_ent, float(R['events'][-1][0]) if R['events'] else float('nan'), cs, ss, As)
    return dict(dmin=dmin, q=q, cs=list(cs), ss=list(ss), As=list(As), entries=n_ent,
                tau_d=float(R['events'][-1][0]) if R['events'] else float('nan'), stalled=stalled, flow=flow, law=law_,
                status=int(sol.status), loops=loops, fine=fine, secs=time.time() - t0)


def path_duals_at(K, f, R):
    """The path's duals at the dipole block's entry τ_d (referee round 68): `path_duals_nd` gives them at the path's history times, and
    the path stops at τ_d, which has none, so they are advanced from the last history time with that interval's residual. Returns (u,
    act, the dictionary's θ, τ_d)."""
    P = RN.path_duals_nd(K, f, R)
    td = float(R['events'][-1][0])
    act = np.asarray(R['act'], bool)
    cen = np.asarray(P['centres'][-1])
    return P['U'][-1] + np.where(act, 0.0, (td - P['taus'][-1]) * K.corr(cen, P['res'][-1])[0]), act, cen, td


def first_wall(u, act, v):
    """The first inactive kernel whose dual u, moving at the correlations v, reaches its wall ±1: (its index, the time) or (None, inf)."""
    ina = ~np.asarray(act, bool) & (v != 0)
    if not ina.any():
        return None, float('inf')
    dt = np.full(v.size, np.inf)
    dt[ina] = (1 - np.sign(v[ina]) * u[ina]) / np.abs(v[ina])
    j = int(np.argmin(dt))
    return j, float(dt[j])


def t3_from_duals(K, f, R, theta, w):
    """τ₃ from the sliding path alone (referee round 68): the path's duals at the dipole block's entry τ_d (`path_duals_at`) moved at
    the correlations of the block's residual at the end of its flow on the critical curve (θ, w at s = 10¹¹: the dipole's asymptotic
    state, which the loop reaches in Λ(ε)τ); the first inactive kernel to reach its wall ±1 (`first_wall`)."""
    u, act, cen, td = path_duals_at(K, f, R)
    j, dt = first_wall(u, act, K.corr(cen, f - K.synth(np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)))[0])
    return None if j is None else float(td + dt)


def d3phi(tb, p, x=XS):
    """∇³_θφ(x; θ̄):p³, the third derivative of φ(x; θ̄ + tp) in t at 0 (referee round 74; analytic): with u = x − c and
    h = −u²/(2σ²) along the line, h′ = p_cu/σ² + p_σu²/σ³, h″ = −p_c²/σ² − 4p_cp_σu/σ³ − 3p_σ²u²/σ⁴, h‴ = 6p_c²p_σ/σ³ +
    18p_cp_σ²u/σ⁴ + 12p_σ³u²/σ⁵, and (e^h)‴ = (h‴ + 3h′h″ + h′³)e^h."""
    (c, sg), (pc, ps) = np.asarray(tb, float), np.asarray(p, float)
    u = np.asarray(x, float) - c
    h1 = pc * u / sg ** 2 + ps * u ** 2 / sg ** 3
    h2 = -pc ** 2 / sg ** 2 - 4 * pc * ps * u / sg ** 3 - 3 * ps ** 2 * u ** 2 / sg ** 4
    h3 = 6 * pc ** 2 * ps / sg ** 3 + 18 * pc * ps ** 2 * u / sg ** 4 + 12 * ps ** 3 * u ** 2 / sg ** 5
    return (h3 + 3 * h1 * h2 + h1 ** 3) * np.exp(-u ** 2 / (2 * sg ** 2))


def k_third(r, tb, p):
    """The dipole's energy coefficient from the third order (referee round 74): at the optimal dipole fit (residual r, midpoint θ̄,
    moment p) E(v) = E_inf + K/v² + O(v⁻⁴) with K = −⟨r, ∇³φ:p³⟩/24 (⟨·⟩ the samples' mean, as E = ½⟨r²⟩). Referee round 80: O(v⁻⁴),
    not O(v⁻³) — about the midpoint the pair's sum has only even powers of 1/v (W multiplies the even orders of Δθ = p/v, w_a − w_b = 2v
    the odd ones)."""
    return float(-np.mean(np.asarray(r, float) * d3phi(tb, p)) / 24)


def kappa_from_path(K, f, R, theta, w):
    """The dipole's intermediate phase from the sliding path alone (referee round 74). At the dipole block's asymptotic state (θ, w: its
    flow's end at s = 10¹¹) the optimal dipole fit (`closure_fit`) and K (`k_third`); τ₃ and its entrant from the duals at τ_d moved by
    that state's residual (`t3_from_duals`' rule); the entrant added at its seed with the coefficient 10⁻³ of its correlation's sign
    and the closure refitted from the first fit — the three-kernel block's optimal dipole fit and its K′ by the same formula —; the
    duals moved on from τ₃ at that fit's residual: τ₄, the first to reach a wall. While the three-kernel block holds, its pair follows
    the block's law with K′: v³ = 3μ₀Λ(ε)(K(τ₃ − τ_d) + K′(τ − τ₃)), so v_land = κv(τ₃) with κ³ = 1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)).
    Referee round 80: the three-kernel block's closure rms (rms2) — below CLOSURE_EXACT the block fits exactly: K′ = 0, no τ₄, κ = 1
    (exact) —; for K′ < 0 the pair retracts and v³ would reach 0 at t_diss = τ₃ + K(τ₃ − τ_d)/|K′| (the dipole dissolves unless τ₄ comes
    first); κ only where K > 0 and the bracket is positive.
    Referee round 83: the block's kernels at τ_d and the entrants' indices j3, j4 (SQ compares the loop's with them).
    Returns dict(K_third, Kp, t3_path, t4_path, kappa_path, t_diss, rms2, exact, block, j3, j4), None where undefined."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    out = dict(K_third=None, Kp=None, t3_path=None, t4_path=None, kappa_path=None, t_diss=None, rms2=None, exact=False, block=None,
               j3=None, j4=None)
    dp = dipole_pair(th, w)
    if dp is None:
        return out
    a, b = dp[:2]
    oth = [i for i in range(w.size) if i not in (a, b)]
    tb = 0.5 * (th[a] + th[b])
    c1 = closure_fit(K, f, th[oth], w[oth], w[a] + w[b], tb, w[a] * (th[a] - tb) + w[b] * (th[b] - tb))
    out['K_third'] = k1 = k_third(c1['r'], c1['tb'], c1['p'])
    u, act, cen, td = path_duals_at(K, f, R)
    out['block'] = [int(i) for i in np.flatnonzero(act)]
    v1 = K.corr(cen, f - K.synth(th, w))[0]
    j3, d3 = first_wall(u, act, v1)
    if j3 is None:
        return out
    out['t3_path'], out['j3'] = td + d3, int(j3)
    u = u + np.where(act, 0.0, d3 * v1)
    act = act.copy()
    act[j3], u[j3] = True, np.sign(v1[j3])
    c2 = closure_fit(K, f, np.r_[c1['th'], cen[j3][None, :]], np.r_[c1['w'], np.sign(v1[j3]) * 1e-3], c1['W'], c1['tb'], c1['p'],
                     max_nfev=50000)
    out['Kp'] = k2 = k_third(c2['r'], c2['tb'], c2['p'])
    out['rms2'] = float(np.sqrt(np.mean(c2['r'] ** 2)))
    if out['rms2'] < CLOSURE_EXACT:                                    # referee round 80: an exact fit — K′ and τ₄ are roundoff
        out.update(Kp=0.0, exact=True, kappa_path=1.0)
        return out
    if k1 > 0 and k2 < 0:                                              # referee round 80: the pair retracts — v³ = 0 at t_diss
        out['t_diss'] = out['t3_path'] + k1 * d3 / abs(k2)
    j4, d4 = first_wall(u, act, K.corr(cen, c2['r'])[0])
    if j4 is None:
        return out
    out['t4_path'], out['j4'] = out['t3_path'] + d4, int(j4)
    br = 1 + k2 * d4 / (k1 * d3)
    out['kappa_path'] = float(np.cbrt(br)) if k1 > 0 and br > 0 else None   # referee round 80: no κ past the dissolution
    return out


def gt_job(args):
    """GT: `t3_from_duals` and `kappa_from_path` (referee round 74) for one of G's targets, its path rerun as `g_job`'s and its block's
    flow end from G."""
    q, cs, ss, As, n_ent, theta, w = args
    L = law()
    K = MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    R = RN.sliding_iss_nd(L, K, f, theta0(), np.ones(C0.size), max_entries=n_ent, exit_kkt=True)
    return dict(q=q, t3_pred=t3_from_duals(K, f, R, theta, w), **kappa_from_path(K, f, R, theta, w))


def part_GT(G):
    """GT (referee round 68's τ₃ from the path, corrected: G's stored t3_pred used the duals at the last history time, before τ_d)."""
    return [gt_job((r['q'], np.array(r['cs']), np.array(r['ss']), np.array(r['As']), int(r['entries']), r['flow'][-1]['theta'],
                    r['flow'][-1]['w'])) for r in G if r['status'] == 0]


def sq_job(args):
    """SQ (referee round 83): the loop's formation sequence by kernel index — `g_job`'s loop rerun at each of G's masses with G's outputs
    (so G's steps) —: at each output where the formed set (|w| above 1% of the smallest amplitude) changes, the kernels newly formed and
    those no longer formed. Returns dict(q, seq={ε: [(τ, formed, unformed)]})."""
    q, cs, ss, As = args
    L, K = law(), MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    tg = np.unique(np.r_[TAUS_G, TAUS_FINE])
    seq = {}
    for eps in EPS_G + EPS_GX:
        _, W, _, _, _ = RN.loop_nd(L, K, f, theta0(), eps, np.r_[0.0, tg], rtol=1e-9, atol=1e-12)
        prev, ev = set(), []
        for t, tau in enumerate(tg, start=1):
            fm = set(np.flatnonzero(np.abs(W[t]) > 1e-2 * np.abs(As).min()).tolist())
            if fm != prev:
                ev.append((float(tau), sorted(fm - prev), sorted(prev - fm)))
            prev = fm
        seq[eps] = ev
    return dict(q=q, seq=seq)


def sq_args(G):
    """SQ's jobs: G's targets with a block flow (status 0), as GT's."""
    return [(r['q'], np.array(r['cs']), np.array(r['ss']), np.array(r['As'])) for r in G if r['status'] == 0]


def part_SQ(G):
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(sq_job, sq_args(G), chunksize=1)


def sq_check(ev, td, n0, g):
    """SQ against the path (referee round 83): the loop's τ₃ and τ₄ as `fine_phases`' (on G's outputs, the first past τ_d with more than
    n0 formed kernels, then the next with more) and its formed sets there against the path's — the block at τ_d with j3, then with j4
    too (`kappa_from_path` in GT; an exact closure has no τ₄) —. Returns (τ₃, τ₄, the loop's set at τ₃, at τ₄, on the path's sequence)."""
    st, k, s = [], 0, set()
    for tau in np.unique(np.r_[TAUS_G, TAUS_FINE]):
        while k < len(ev) and ev[k][0] <= tau:
            s = (s | set(ev[k][1])) - set(ev[k][2])
            k += 1
        st.append((float(tau), frozenset(s)))
    a3 = next(((t, s_) for t, s_ in st if t > td and len(s_) > n0), None)
    a4 = next(((t, s_) for t, s_ in st if a3 is not None and t > a3[0] and len(s_) > len(a3[1])), None)
    w3 = None if g.get('j3') is None or g.get('block') is None else set(g['block']) | {g['j3']}
    w4 = None if w3 is None or g.get('j4') is None else w3 | {g['j4']}
    ok = (w3 is not None and a3 is not None and set(a3[1]) == w3
          and ((a4 is not None and set(a4[1]) == w4) if w4 is not None else (a4 is None and bool(g.get('exact')))))
    return (None if a3 is None else a3[0], None if a4 is None else a4[0], None if a3 is None else sorted(a3[1]),
            None if a4 is None else sorted(a4[1]), bool(ok))


def pair_p(th, w, a, b):
    """(ii)(a)'s expansion parameter of the pair (a, b): p = |w_a − w_b||Δθ|/(|W|σ̄), |Δθ| = (Δc² + Δσ²)^{1/2} (referee round 74: a
    σ-dipole's Δc alone understates it), σ̄ the pair's mean |σ|; inf at W = 0."""
    th, w = np.atleast_2d(np.asarray(th, float)), np.asarray(w, float)
    W = w[a] + w[b]
    if W == 0:
        return float('inf')
    return float(abs(w[a] - w[b]) * np.hypot(th[a, 0] - th[b, 0], th[a, 1] - th[b, 1]) / (abs(W) * 0.5 * (abs(th[a, 1]) + abs(th[b, 1]))))


def fine_phases(sts, n0, td, cs, ss, As):
    """G on the fine grid (referee rounds 68, 74): from the states (τ, formed θ, formed w), the block's dipole phase and what follows — τ₃
    (the first τ past τ_d with more than n0 formed kernels) and τ₄ (the next formation); the dipole pair's v over the times in
    (τ_d, τ₃) with n0 formed (every tenth, for the slopes of v³) and at the last of them (v3 at t_v3); over the intermediate phase,
    the times in [τ₃, τ₄) with n3 formed, v at every tenth and the last (dip2: the slope of v³ there); the intermediate phase at
    (τ₃ + τ₄)/2 (formed, the pair's offset from the nearest feature (σ*) and (W − A)/|A|); the first landing past τ₃ (`landing`:
    centre of mass and W within LAND_TOL) — its τ, v and (ii)(a)'s p (`pair_p`) —; and over every state with a dipole pair and p < 1
    whether p ≤ P_LAND and the landing agree (agree), and over those past τ₄ (agree4)."""
    ns = [(tau, len(w)) for tau, th, w in sts]
    t3 = next((tau for tau, n in ns if tau > td and n > n0), None)
    n3 = next((n for tau, n in ns if tau == t3), None)
    t4 = next((tau for tau, n in ns if t3 is not None and tau > t3 and n > n3), None)
    dip = [(tau, dipole_pair(th, w)) for tau, th, w in sts if td < tau < (t3 if t3 is not None else np.inf) and len(w) == n0]
    dip = [(tau, d[2]) for tau, d in dip if d is not None]
    mid = [(tau, dipole_pair(th, w)) for tau, th, w in sts if t4 is not None and t3 <= tau < t4 and len(w) == n3]
    mid = [(tau, d[2]) for tau, d in mid if d is not None]
    out = dict(t3=t3, t4=t4, v3=dip[-1][1] if dip else None, t_v3=dip[-1][0] if dip else None, dip=[(a, b) for a, b in dip[::10]],
               dip2=[(a, b) for a, b in mid[::10]] + ([mid[-1]] if mid and (len(mid) - 1) % 10 else []), inter=None, t_land=None,
               v_land=None, p_land=None, agree=None, agree4=None)
    if t3 is not None and t4 is not None:
        mid = next(((th, w) for tau, th, w in sts if tau >= 0.5 * (t3 + t4)), None)
        lp = landing(mid[0], mid[1], cs, ss, As) if mid else None
        out['inter'] = None if lp is None else dict(n=len(mid[1]), off=lp['off'], dW=lp['dW'])
    agree, tot, agree4, tot4 = 0, 0, 0, 0
    for tau, th, w in sts:
        lp = landing(th, w, cs, ss, As)
        if lp is None:
            continue
        a, b = dipole_pair(th, w)[:2]
        p = pair_p(th, w, a, b)
        ok = lp['off'] <= LAND_TOL and abs(lp['dW']) <= LAND_TOL
        if p < 1:
            tot += 1
            agree += (p <= P_LAND) == ok
            if t4 is not None and tau > t4:
                tot4 += 1
                agree4 += (p <= P_LAND) == ok
        if ok and out['t_land'] is None and t3 is not None and tau > t3:
            out.update(t_land=tau, v_land=lp['v'], p_land=float(p))
    out['agree'], out['agree4'] = (agree, tot), (agree4, tot4)
    return out


def gl_job(args):
    """GL (referee rounds 68, 74): the loop on one dipole target at the mass EPS_GL to TAUS_GL — per output the dipole pair among the
    formed kernels: its distance, centre of mass's offset and (W − A)/|A| (`landing`), |Δc|, w_a, w_b, its mean width and the nearest
    feature's amplitude, the residual, the energy ½⟨r²⟩ of all kernels and of the formed ones —, with Λ(ε), Gm and μ₀ for the collapse
    law's constant (Δc⁻⁶ ∝ τ − τ_L, Δc ∝ τ^{−1/6}, the residual ∝ τ^{−2/3}, W − A ∝ τ^{−1/3} as Δc², E ∝ Δc⁸ with the centre of mass
    pinned)."""
    q, cs, ss, As = args
    t0 = time.time()
    L, K = law(), MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    _, W, C, Lam, _ = RN.loop_nd(L, K, f, theta0(), EPS_GL, np.r_[0.0, TAUS_GL], rtol=1e-9, atol=1e-12)
    rows = []
    for t, tau in enumerate(TAUS_GL, start=1):
        big = np.abs(W[t]) > 1e-2 * np.abs(As).min()
        th, w = C[t][big], W[t][big]
        lp = landing(th, w, cs, ss, As)
        r = f - K.synth(C[t], W[t])
        row = dict(tau=float(tau), d=None if lp is None else lp['d'], dW=None if lp is None else lp['dW'], rms=float(np.sqrt(np.mean(r ** 2))),
                   E=float(np.mean(r ** 2) / 2), E_formed=float(np.mean((f - K.synth(th, w)) ** 2) / 2))
        if lp is not None:
            a, b = dipole_pair(th, w)[:2]
            cm = (w[a] * th[a, 0] + w[b] * th[b, 0]) / (w[a] + w[b])
            k = int(np.argmin(np.abs(np.asarray(cs) - cm) / np.asarray(ss)))
            row.update(dc=float(abs(th[a, 0] - th[b, 0])), wa=float(w[a]), wb=float(w[b]), sbar=float(0.5 * (abs(th[a, 1]) + abs(th[b, 1]))),
                       A=float(As[k]))
        rows.append(row)
    return dict(q=q, rows=rows, Lam=float(Lam), Gm=float(L.G * L.m), mu0=float(4 * L.G * L.c / L.A), secs=time.time() - t0)


def part_GL(G):
    """GL: G's targets Q_GL."""
    return [gl_job((q, np.array(r['cs']), np.array(r['ss']), np.array(r['As']))) for q in Q_GL for r in G if r['q'] == q]


def gl_constant(x0, x1, Gm, mu0, Lam):
    """The collapse law's constant over GL's interval (x0, x1) (referee round 74): the measured d(Δc⁻⁶)/dτ against the proof's
    48aGm(W²/(w_a² + w_b²))Λ, a = E/Δc⁸, and against it over 1 + C — C = Gm W²(W − A)²(1/|w_a| + 1/|w_b|)/(μ₀Δc²(w_a² + w_b²)) the
    slaved coefficients' share of the separation's metric (W − A ∝ Δc² moves at the coefficients' mobility μ₀|w|, far below the centres'
    Gm), (Gm/(16μ₀))W²p²/(vσ²) to leading order —, each averaged over the interval's ends. Returns (p at x1, ratio, ratio with C, C at
    x1)."""
    rate = lambda x: 48 * x['E'] / x['dc'] ** 8 * Gm * (x['wa'] + x['wb']) ** 2 / (x['wa'] ** 2 + x['wb'] ** 2) * Lam
    cc = lambda x: (Gm * (x['wa'] + x['wb']) ** 2 * (x['dW'] * abs(x['A'])) ** 2 * (1 / abs(x['wa']) + 1 / abs(x['wb']))
                    / (mu0 * x['dc'] ** 2 * (x['wa'] ** 2 + x['wb'] ** 2)))
    meas = (x1['dc'] ** -6 - x0['dc'] ** -6) / (x1['tau'] - x0['tau'])
    p = abs(x1['wa'] - x1['wb']) * x1['dc'] / (abs(x1['wa'] + x1['wb']) * x1['sbar'])
    return (float(p), float(meas / (0.5 * (rate(x0) + rate(x1)))),
            float(meas / (0.5 * (rate(x0) / (1 + cc(x0)) + rate(x1) / (1 + cc(x1))))), float(cc(x1)))


def stalled_summary(Rc, f, K):
    """G: a budget-stopped path's stalled block — its kernels, largest |w| (the accepted polish's: solver-dependent), the closest
    opposite-sign pair's distance, the residual, the block's polishes — or None."""
    st = Rc.get('stalled')
    if st is None:
        return None
    cp = closest_pair(st['c'], st['w'], opposite=True)
    return dict(wmax=float(np.abs(st['w']).max()), sep=None if cp is None else cp[0], n=int(np.size(st['w'])),
                rms=float(np.sqrt(np.mean((f - K.synth(np.atleast_2d(st['c']), st['w'])) ** 2))),
                polish=Rc['polish_log'][-1] if Rc['polish_log'] else [], nfev=int(Rc['nfev']))


def dipole_exponents(flow):
    """The last two outputs' log-slopes of the largest |w| and of the opposite-sign pair's distance against s (the dipole: +1/3, −1/3 —
    on the free law's exponential branch the coefficient's mobility grows like |w| and the pair's energy excess falls like 1/w² at a
    fixed moment)."""
    a, b = flow[-2], flow[-1]
    if a['sep'] is None or b['sep'] is None:
        return None
    ls = np.log(b['s'] / a['s'])
    return float(np.log(b['wmax'] / a['wmax']) / ls), float(np.log(b['sep'] / a['sep']) / ls)


def dipole_pair(theta, w):
    """The loop's dipole pair among formed kernels (`dipole_law`'s choice): the largest |w| and the closest kernel of the other sign —
    (a, b, v = |w_a − w_b|/2, its distance |Δc|/σ + |Δlog σ|) or None."""
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    if w.size < 2:
        return None
    a = int(np.argmax(np.abs(w)))
    opp = [j for j in range(w.size) if w[j] * w[a] < 0]
    if not opp:
        return None
    dist = lambda j: (abs(th[a, 0] - th[j, 0]) / (0.5 * (abs(th[a, 1]) + abs(th[j, 1]))) + abs(np.log(abs(th[a, 1]) / abs(th[j, 1]))))
    b = min(opp, key=dist)
    return a, b, float(abs(w[a] - w[b]) / 2), float(dist(b))


def landing(theta, w, cs, ss, As):
    """The dipole pair against the features: its centre of mass's offset from the nearest feature (in that feature's width) and its net
    coefficient's relative difference from the feature's amplitude, (v, distance, offset, (W − A)/A) — landed if both within LAND_TOL."""
    dp = dipole_pair(theta, w)
    if dp is None:
        return None
    a, b, v, d = dp
    th, w = np.atleast_2d(np.asarray(theta, float)), np.asarray(w, float)
    W = w[a] + w[b]
    cm = (w[a] * th[a, 0] + w[b] * th[b, 0]) / W if W != 0 else 0.5 * (th[a, 0] + th[b, 0])
    k = int(np.argmin(np.abs(np.asarray(cs) - cm) / np.asarray(ss)))
    return dict(v=v, d=d, off=float(abs(cm - cs[k]) / ss[k]), dW=float((W - As[k]) / abs(As[k])), w=(float(w[a]), float(w[b])))


def phases(r, eps):
    """G: the loop at the mass eps on the dipole's target r, in three phases (referee round 61). The dipole: from the block's entry τ_d
    (the path's) the pair's v grows as v³ = 3μ₀KΛ(ε)(τ − τ_d) (the block's law in its fast time Λ(ε)τ) until the next kernel forms, at
    τ₃ (the first grid time with more formed kernels than the block's); the slopes of v³ over the grid times strictly between τ_d and
    τ₃ against 3μ₀KΛ (median, range). The landing: v at τ = 400 against the law's (3μ₀KΛ(ε)(τ₃ − τ_d))^{1/3}, and whether the pair is
    landed (`landing`) at 400, 1600, 6400. The collapse: over 1600 → 6400 the exponents of the pair's distance and of the residual in τ."""
    xs = sorted([x for x in r['loops'] if x['eps'] == eps], key=lambda x: x['tau'])
    lw = r.get('law') or []
    if not xs or not lw or lw[-1] is None:
        return None
    K_ = lw[-1]['v'] ** 2 * (lw[-1]['E'] - lw[-1]['E_inf'])
    rate = 3 * lw[-1]['mu0'] * K_ * xs[0]['Lam']
    n0, td = r['entries'], r['tau_d']
    t3 = next((x['tau'] for x in xs if x['tau'] > td and len(x['w']) > n0), None)
    at = {x['tau']: x for x in xs}
    win = [x for x in xs if td < x['tau'] < (t3 if t3 is not None else np.inf) and len(x['w']) == n0]
    vs = [(x['tau'], dipole_pair(x['theta'], x['w'])) for x in win]
    vs = [(tau, dp[2]) for tau, dp in vs if dp is not None]
    sl = [(b[1] ** 3 - a[1] ** 3) / (b[0] - a[0]) / rate for a, b in zip(vs[:-1], vs[1:])]
    land = {tau: landing(at[tau]['theta'], at[tau]['w'], r['cs'], r['ss'], r['As']) if tau in at else None for tau in TAUS_G_MD}
    v400 = land[400.0]['v'] if land.get(400.0) else None
    vstar = (rate * (t3 - td)) ** (1 / 3) if t3 is not None else None
    col = None
    a, b = land.get(1600.0), land.get(6400.0)
    if a and b and a['d'] > 0 and b['d'] > 0:
        lt = np.log(6400.0 / 1600.0)
        col = (float(np.log(b['d'] / a['d']) / lt), float(np.log(at[6400.0]['rms'] / at[1600.0]['rms']) / lt))
    return dict(rate=rate, K=K_, t3=t3, slopes=sl, v400=v400, vstar=vstar, land=land, collapse=col)


def part_G(C):
    """G: C's targets whose moving-width path was stopped by its budget."""
    from multiprocessing import Pool
    jobs = [(r['dmin'], r['q'], (np.array(r['cs']), np.array(r['ss']), np.array(r['As'])), int(r['moving']['entries'])) for r in C
            if r['moving'].get('oob')]
    with Pool(N_WORKERS) as pool:
        return pool.map(g_job, jobs, chunksize=1)


def visibility(th):
    """Each kernel's largest φ on the samples (0 at σ = 0): below VIS it is invisible there (referee rounds 68, 74)."""
    th = np.atleast_2d(np.asarray(th, float)).reshape(-1, 2)
    return np.array([float(np.exp(-((XS - t[0]) ** 2).min() / (2 * t[1] ** 2))) if t[1] != 0 else 0.0 for t in th])


def path_stop(Rm, tol, max_events=40):
    """Why a path ended (referee round 74): 'tolerance' (the last block's rms below tol), 'budget' (its evaluations spent), 'dictionary'
    (every kernel active), 'events' (max_events blocks) or 'no entrant'."""
    if Rm['hist'] and Rm['hist'][-1][3] < tol:
        return 'tolerance'
    if Rm.get('out_of_budget'):
        return 'budget'
    if np.all(Rm['act']):
        return 'dictionary'
    return 'events' if len(Rm['hist']) >= max_events else 'no entrant'


def h_stats(blocks, m, As):
    """H's statistics of a path — blocks: each block's active (θ, w) — and of its end m (referee rounds 57, 61, 74), counting only
    kernels visible on the samples (`visibility` ≥ VIS): the most kernels of the other sign (below −1% of the smallest amplitude) in a
    block and the most negative coefficient; the end's formed kernels of the other sign (and all of them, the invisible included) and
    the smallest formed width; over all the end's active kernels those narrower than the samples' spacing and those of negative width
    (the widths are unconstrained), and the formed ones invisible on the samples; the largest |w|."""
    w, th = np.asarray(m['w'], float), np.atleast_2d(np.asarray(m['theta'], float)).reshape(-1, 2)
    amin = min(abs(a) for a in As)
    formed = np.abs(w) > 1e-2 * amin
    vis = visibility(th) >= VIS
    seen = [(np.asarray(wb, float), visibility(tb) >= VIS) for tb, wb in blocks]
    return dict(neg_path=int(max([int(((wb < -1e-2 * amin) & vb).sum()) for wb, vb in seen], default=0)),
                wmin_path=float(min([wb[vb].min() for wb, vb in seen if vb.any()], default=0.0)),
                neg_end=int((formed & vis & (w < 0)).sum()), neg_end_all=int((formed & (w < 0)).sum()),
                smin=float(np.abs(th[formed & vis, 1]).min()) if (formed & vis).any() else float('nan'),
                sub=int((np.abs(th[:, 1]) < XS[1] - XS[0]).sum()) if w.size else 0, negw=int((th[:, 1] < 0).sum()) if w.size else 0,
                sub_formed=int((formed & (np.abs(th[:, 1]) < XS[1] - XS[0])).sum()) if w.size else 0,
                negw_formed=int((formed & (th[:, 1] < 0)).sum()) if w.size else 0,
                inv_formed=int((formed & ~vis).sum()), wmax=float(np.abs(w).max()) if w.size else 0.0)


def clusters(m, cs, ss, As):
    """H2: the end's formed kernels by the feature nearest each (in its width): per feature their count and their mass over the
    feature's, Σw_i|σ_i|/(A_kσ_k) (one kernel fitting it exactly: 1)."""
    w, th = np.asarray(m['w'], float), np.atleast_2d(np.asarray(m['theta'], float))
    formed = np.flatnonzero(np.abs(w) > 1e-2 * min(abs(a) for a in As))
    near = [int(np.argmin(np.abs(np.asarray(cs) - th[i, 0]) / np.asarray(ss))) for i in formed]
    return [(sum(k == j for k in near), float(sum(w[i] * abs(th[i, 1]) for i, k in zip(formed, near) if k == j) / (As[j] * ss[j])))
            for j in range(len(cs))]


def h_job(args):
    """H: two features of one sign — the moving-width path's outcome and `h_stats` (referee round 57); with a tolerance (H2, referee
    round 61) the path run to it instead of 1e-9: the last entry's τ, the last blocks' residuals and the formed kernels per feature."""
    name, cs, ss, As = args[:4]
    tol = args[4] if len(args) > 4 else None
    pol = bool(args[5]) if len(args) > 5 else True
    t0 = time.time()
    f = feature_sum(cs, ss, As)
    if tol is None and pol:
        Rm, _ = run_paths(f, fixed=False)
    else:
        Rm = RN.sliding_iss_nd(law(), MW.BumpsCW(XS, S0), f, theta0(), np.ones(C0.size), max_events=40, exit_kkt=True,
                               max_rhs=MAX_RHS_MOVING, tol=1e-9 if tol is None else tol, polish=pol)
    m = match(Rm, cs, ss, As)
    blocks = [(h[1][h[2] != 0], h[2][h[2] != 0]) for h in Rm['hist'][1:]]
    out = dict(name=name, cs=list(cs), ss=list(ss), As=list(As), moving=m, outcome=outcome(m, As, cs, ss), tol=tol, polish=pol,
               stop=path_stop(Rm, 1e-9 if tol is None else tol), secs=time.time() - t0, **h_stats(blocks, m, As), **stop_anatomy(Rm))
    if tol is not None:
        out.update(last_tau=float(Rm['events'][-1][0]) if Rm['events'] else float('nan'),
                   rms_blocks=[float(h[3]) for h in Rm['hist'][-6:]], clusters=clusters(m, cs, ss, As))
    return out


def stop_anatomy(Rm):
    """Referee round 68: why a path ended and what its polishes did — the active kernels (all of them: the dictionary exhausted), the
    accepted polishes and those of them stopped by their budget (status 0), the active kernels of negative width and those invisible on
    the samples (largest φ below 10⁻¹⁰)."""
    act = np.flatnonzero(Rm['act'])
    th = np.asarray(Rm['c'], float)[act]
    pl = [p_ for b in Rm.get('polish_log', []) for p_ in b if p_.get('accepted')]
    vis = visibility(th) if act.size else np.zeros(0)
    return dict(n_act=int(act.size), exhausted=bool(act.size == C0.size), pol_acc=len(pl),
                pol_capped=sum(p_.get('status') == 0 for p_ in pl), neg_sigma=int(np.sum(th[:, 1] < 0)) if act.size else 0,
                invisible=int(np.sum(vis < VIS)))


def part_HP():
    """HP (referee round 68): H's path without the least-squares polish on the symmetric targets and two asymmetric ones."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, [tuple(x) + (None, False) for x in H_TARGETS if x[0] in HP_NAMES], chunksize=1)


def part_H():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, list(H_TARGETS), chunksize=1)


def part_H2():
    """H2: the σ* = 0.06 targets' paths to TOL_H2 (referee round 61: below 1e-9 the 'exact' ends were the path's stopping tolerance)."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, [tuple(x) + (TOL_H2,) for x in H_TARGETS if x[0].startswith('σ* 0.06')], chunksize=1)


def h3_targets():
    return [(f"σ* 0.06, ±{s:g}", (-s, s), (0.06, 0.06), (1.0, 1.0)) for s in SEP_H]


def part_H3():
    """H3: the separation scan at σ* = 0.06, equal amplitudes (H's ±0.08 besides; referee round 61)."""
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(h_job, h3_targets(), chunksize=1)


def hl_job(args):
    """HL: the real loop on one of H's targets at the masses EPS_HL to the times TAUS_HL: the formed kernels (|w| above 1% of the
    smallest amplitude), their states, the largest |w|, the residual (referee round 61)."""
    name, cs, ss, As = args
    t0 = time.time()
    L, K = law(), MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    out = []
    for eps in EPS_HL:
        _, W, C, Lam, _ = RN.loop_nd(L, K, f, theta0(), eps, np.r_[0.0, TAUS_HL], rtol=1e-9, atol=1e-12)
        for t, tau in enumerate(TAUS_HL, start=1):
            big = np.flatnonzero(np.abs(W[t]) > 1e-2 * min(abs(a) for a in As))
            out.append(dict(eps=eps, tau=float(tau), theta=C[t][big].tolist(), w=W[t][big].tolist(), wmax=float(np.abs(W[t]).max()),
                            rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2))), neg=int(np.sum(W[t][big] < 0)),
                            smin=float(np.abs(C[t][big][:, 1]).min()) if big.size else float('nan')))
    return dict(name=name, cs=list(cs), ss=list(ss), As=list(As), loops=out, secs=time.time() - t0)


def part_HL():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(hl_job, [x for x in H_TARGETS if x[0] in HL_NAMES], chunksize=1)


SPLIT_CS = 1.5                                                          # referee round 97 (D3cr-4): a main kernel's (c, σ) error at most this, w free
NULL_TOL = 1e-8                                                         # referee round 97 (D3cr-3): a null cluster's others' rms over max f* (the number of `outcome`'s residual test)


def fit_split(th, w, cs, ss, As, cap=None):
    """The fit's split (referee rounds 74, 80, 83): the main kernels — distinct formed kernels, one per feature of its sign, the
    assignment minimizing the largest relative error in (c, σ, w) (c over σ*), then the next largest (over all assignments; a feature
    left without one only if its sign's kernels run out) —, each one's error and signed w/A − 1 (inflated, > 0, the others' net weight
    of the other sign: a cancellation; deficient, < 0, the others making up the rest: a split); the other formed kernels: their number,
    Σw and Σ|w|, their sum's largest value on the samples over max f* (a null cluster: small against the residual, beside exact main
    kernels), their widths' range and their centres' largest |c|. th, w: the formed kernels. Referee round 87 (D3ch-3): a kernel is a
    feature's main kernel only within a cap of it (else the feature is left unassigned, its error None: no broad kernel made a main
    one by the tie-break); the others' signed mass √(2π)Σw|σ| (the side of a split or a cancellation, which Σw's sign need not tell) and
    their sum's rms on the samples over max f* (the null cluster's norm). Referee round 97 (D3cr-2, D3cr-4): the cap SPLIT_CS acts on the
    (c, σ) error alone, w free — a cancellation's main kernel is inflated in w, a transient broad kernel is off in (c, σ) —; main_cs the
    main kernels' (c, σ) errors, near_cs an unassigned feature's nearest unused kernel of its sign in (c, σ); rest_win the others' mass on the
    samples' window, Σ_x oth(x)Δx (where the fit lives: the continuum mass also counts tails beyond the window)."""
    th, w = np.atleast_2d(np.asarray(th, float)).reshape(-1, 2), np.asarray(w, float).reshape(-1)
    ecs = lambda i, k: float(max(abs(th[i, 0] - cs[k]) / ss[k], abs(abs(th[i, 1]) / ss[k] - 1)))
    er = lambda i, k: max(ecs(i, k), float(abs(w[i] / As[k] - 1)))
    best = None
    for combo in itertools.product(*[[i for i in range(w.size) if np.sign(w[i]) == np.sign(A_) and ecs(i, k) <= (SPLIT_CS if cap is None else cap)] + [None]
                                     for k, A_ in enumerate(As)]):
        used = [i for i in combo if i is not None]
        if len(set(used)) < len(used):
            continue
        key = (len(combo) - len(used), sorted((er(i, k) for k, i in enumerate(combo) if i is not None), reverse=True))
        if best is None or key < best[0]:
            best = (key, list(combo))
    main = best[1]
    err = [None if i is None else er(i, k) for k, i in enumerate(main)]
    dw = [None if i is None else float(w[i] / As[k] - 1) for k, i in enumerate(main)]
    rest = [i for i in range(w.size) if i not in main]
    oth = MW.BumpsCW(XS, S0).synth(th[rest], w[rest]) if rest else None
    mf = float(np.max(np.abs(feature_sum(cs, ss, As))))
    foot = float(np.max(np.abs(oth))) / mf if rest else 0.0
    rms_ = float(np.sqrt(np.mean(oth ** 2))) / mf if rest else 0.0
    near = [None if i is not None else min([ecs(j, k) for j in rest if np.sign(w[j]) == np.sign(As[k])] or [None],
                                           key=lambda v: np.inf if v is None else v) for k, i in enumerate(main)]
    return dict(main=list(main), main_cs=[None if i is None else ecs(i, k) for k, i in enumerate(main)], near_cs=near,
                rest_win=float(np.sum(oth) * (XS[1] - XS[0])) if rest else 0.0, main_err=err, main_dw=dw, n_rest=len(rest), rest_sw=float(w[rest].sum()), rest_w=float(np.abs(w[rest]).sum()), rest_foot=foot,
                rest_sig=[float(np.abs(th[rest, 1]).min()), float(np.abs(th[rest, 1]).max())] if rest else None,
                rest_c=float(np.abs(th[rest, 0]).max()) if rest else None,
                rest_mass=float(np.sqrt(2 * np.pi) * np.sum(w[rest] * np.abs(th[rest, 1]))) if rest else 0.0, rest_rms=rms_)


def null_cluster(th, w, cs, ss, As, rms):
    """The one null-cluster predicate (referee round 100, D3cu-3), `outcome`'s and `split_summary`'s: the residual's rms below 1e-8, every
    feature's main kernel (`fit_split`) within EXACT_TOL of it in c, σ and w (absolute, as `match`), other formed kernels present, and
    their sum's rms on the samples below NULL_TOL of max f*."""
    th, w = np.atleast_2d(np.asarray(th, float)).reshape(-1, 2), np.asarray(w, float).reshape(-1)
    sp = fit_split(th, w, cs, ss, As)
    if rms >= 1e-8 or not sp['n_rest'] or any(i is None for i in sp['main']):
        return False
    err = max(max(abs(th[i, 0] - cs[k]), abs(abs(th[i, 1]) - ss[k]), abs(w[i] - As[k])) for k, i in enumerate(sp['main']))
    return bool(err < EXACT_TOL and sp['rest_rms'] < NULL_TOL)


def cap_interval(inputs):
    """The caps over which every split's assignment is the one at SPLIT_CS (referee round 100, D3cu-5): an assignment can change only
    where the cap crosses a same-sign kernel's (c, σ) error off a feature; from SPLIT_CS the next such errors up and down are tried in
    turn (`fit_split` at each, just past it) until one changes an assignment. inputs: (th, w, cs, ss, As) per split. Returns (lo, hi)."""
    base = [fit_split(*x)['main'] for x in inputs]
    errs = set()
    for th, w, cs, ss, As in inputs:
        th, w = np.atleast_2d(np.asarray(th, float)).reshape(-1, 2), np.asarray(w, float).reshape(-1)
        for k in range(len(cs)):
            for i in range(w.size):
                if np.sign(w[i]) == np.sign(As[k]):
                    errs.add(float(max(abs(th[i, 0] - cs[k]) / ss[k], abs(abs(th[i, 1]) / ss[k] - 1))))
    same = lambda c: all(fit_split(*x, cap=c)['main'] == b for x, b in zip(inputs, base))
    hi = next((e for e in sorted(v for v in errs if v > SPLIT_CS) if not same(e * (1 + 1e-12))), np.inf)
    lo = next((e for e in sorted((v for v in errs if v <= SPLIT_CS), reverse=True) if not same(e * (1 - 1e-12))), 0.0)
    return lo, hi


def split_summary(H=None, HP=None, H3=None, HL=None, HLL=None):
    """Referee round 97 (D3cr-3, D3cr-4): over every split of H, H3, HP, HL and HLL — the main kernels' largest (c, σ) error and w/A − 1,
    an unassigned feature's nearest unused kernel of its sign in (c, σ) (the cap SPLIT_CS between), and the null clusters (every main kernel
    within EXACT_TOL, others present, their rms below NULL_TOL of max f*)."""
    rows, nl, inp = [], [], []
    for lab, rs in (('H', H), ('H3', H3), ('HP', HP)):
        for r in (rs if rs is not None else []):
            rows.append((f"{lab} {r['name']}", end_split(r)))
            m = r['moving']
            fm = np.abs(np.asarray(m['w'], float)) > 1e-2 * min(abs(a) for a in r['As'])
            inp.append((np.atleast_2d(np.asarray(m['theta'], float))[fm], np.asarray(m['w'], float)[fm], r['cs'], r['ss'], r['As']))
            if null_cluster(np.atleast_2d(np.asarray(m['theta'], float))[fm], np.asarray(m['w'], float)[fm], r['cs'], r['ss'], r['As'], m['rms']):
                nl.append(rows[-1][0])
    for r in (HL if HL is not None else []):
        for x in r['loops']:
            rows.append((f"HL {r['name']} {x['eps']:g} τ {x['tau']:g}", fit_split(x['theta'], x['w'], r['cs'], r['ss'], r['As'])))
            inp.append((x['theta'], x['w'], r['cs'], r['ss'], r['As']))
            if null_cluster(x['theta'], x['w'], r['cs'], r['ss'], r['As'], x['rms']):
                nl.append(rows[-1][0])
    ht = {x[0]: x[1:] for x in H_TARGETS}
    for r in (HLL if HLL is not None else []):
        for x in r['loops']:
            rows.append((f"HLL {r['name']} {r['eps']:g} τ {x['tau']:g}", fit_split(x['theta'], x['w'], *ht[r['name']])))
            inp.append((x['theta'], x['w']) + tuple(ht[r['name']]))
            if null_cluster(x['theta'], x['w'], *ht[r['name']], x['rms']):
                nl.append(rows[-1][0])
    if not rows:
        return []
    mc = [(v, lab) for lab, sp in rows for v in sp['main_cs'] if v is not None]
    md_ = [(d, lab) for lab, sp in rows for d in sp['main_dw'] if d is not None]
    nr = [(v, lab) for lab, sp in rows for v in sp['near_cs'] if v is not None]
    un = sum(v is None for _, sp in rows for v in sp['main_cs'])
    nul = nl
    lo, hi = cap_interval(inp)
    sh = [lab for lab, sp in rows if any(v is None for v in sp['main_cs']) and any(v is None for v, m_ in zip(sp['near_cs'], sp['main_cs'])
                                                                                     if m_ is None)]
    sym = lambda r: abs(r['cs'][0] + r['cs'][1]) < 1e-12 and r['ss'][0] == r['ss'][1] and r['As'][0] == r['As'][1]
    ma = []                                                              # referee rounds 103, 106 (D3cx-6, D3da-3): the ends' mirror symmetry
    for lab, rs in (('H', H), ('H3', H3), ('HP', HP)):
        for r in (rs if rs is not None else []):
            if sym(r):
                m = r['moving']
                th_, w_ = np.atleast_2d(np.asarray(m['theta'], float)).reshape(-1, 2), np.asarray(m['w'], float).reshape(-1)
                fm = np.abs(w_) > 1e-2 * min(abs(a) for a in r['As'])
                ac = w_ != 0
                Am = max(abs(a) for a in r['As'])
                bn = mirror_asym(th_[fm], w_[fm])
                ma.append(f"{lab} {r['name']} {bn:.2g} from the mirror image ({bn / 2:.2g} from the subspace; all {int(ac.sum())} kernels "
                          f"{mirror_asym(th_[ac], w_[ac]):.2g}), the formed kernels' function odd by {mirror_odd(th_[fm], w_[fm], A=Am):.1e} on "
                          f"and between the samples, the whole state's by {mirror_odd(th_, w_, A=Am, xs=XS[XS >= 0]):.1e} on the samples "
                          f"(the residual's: f* is even) and {mirror_odd(th_, w_, A=Am):.1e} between them"
                          + (lambda wd, nw: f"; the {int(wd.sum())} kernels at least as wide as the samples' spacing (σ ≥ h; the "
                             f"narrowest {np.abs(th_[wd, 1]).min():.3g}" + ('' if not nw.any() else f", the widest narrower one "
                             f"{np.abs(th_[nw, 1]).max():.2g}") + f") odd by {mirror_odd(th_[wd], w_[wd], A=Am):.1e} on and between the samples"
                             + f" ({int((wd & fm).sum())} formed" + ('' if not (wd & fm).any() else
                                f", the narrowest of them {np.abs(th_[wd & fm, 1]).min():.3g}") + f"; {int((wd & ~fm).sum())} idle)"
                             + ('' if not (ac & ~fm).any() else f"; the formed threshold's margin: the largest unformed |w| "
                                f"{np.abs(w_[ac & ~fm]).max() / Am:.2g} (over max A" + ('' if min(abs(a) for a in r['As']) == Am else   # round 125
                                f"; over min A {np.abs(w_[ac & ~fm]).max() / min(abs(a) for a in r['As']):.2g}") + ", "
                                f"against the threshold 0.01 of min A), the smallest formed {np.abs(w_[fm]).min() / Am:.2g} (over max A)"))(
                              ac & (np.abs(th_[:, 1]) >= XS[1] - XS[0]), ac & (np.abs(th_[:, 1]) < XS[1] - XS[0])))   # round 112 (D3dg-2)
    v = [mirror_asym(x['theta'], x['w']) for r in (HL if HL is not None else []) if sym(r) for x in r['loops']]
    v += [mirror_asym(x['theta'], x['w']) for r in (HLL if HLL is not None else []) if sym(dict(zip(('cs', 'ss', 'As'), ht[r['name']])))
          for x in r['loops']]
    if v:
        ma.append(f"the loops' formed kernels (HL, HLL) at most {max(v):.1e}")
    return ['', f"The splits' classes (referee rounds 97, 100): over the {len(rows)} splits above, the main kernels' (c, σ) errors at most "
            f"{max(mc)[0]:.3g} ({max(mc)[1]}" + (", where the other feature's only kernel of its sign is this one, left to it by "
                                                 "distinctness" if max(mc)[1] in sh else '') + f") and their w/A − 1 "
            f"{min(md_)[0]:+.3g} … {max(md_)[0]:+.3g}; {un} features unassigned"
            + (f", {len(nr)} of them with an unused kernel of their sign, at least {min(nr)[0]:.3g} off in (c, σ) ({min(nr)[1]})" if nr else '')
            + f" — the cap {SPLIT_CS:g} on the (c, σ) error alone (w free) lies in that gap, and every split's assignment is the same for "
              f"every cap in ({math.ceil(lo * 1e4) / 1e4:.4f}, {math.floor(hi * 1e4) / 1e4 if np.isfinite(hi) else hi:.4f}) (`cap_interval`, "
              "rounded inward; "
              "where a same-sign kernel's (c, σ) error is crossed); null clusters (`null_cluster`, `outcome`'s "
            f"test too: the residual below 1e-08, every main kernel within {EXACT_TOL:g} in c, σ and w, the others' rms below {NULL_TOL:g} of "
            "max f*; referee round 100): " + (', '.join(nul) if nul else 'none') + '.'
            + (" On the targets symmetric about 0 the formed kernels' distance from their mirror image (`mirror_asym`, a bottleneck "
               "matching; half of it from the reflection's fixed-point subspace; referee rounds 103, 106, 112): " + '; '.join(ma) + "." if ma else ''), '']


def split_cells(sp):
    """`fit_split`'s two table cells: the main kernels' errors (their signed w/A − 1); the others' n, Σw (Σ|w|), sum on the samples,
    widths, largest |c|."""
    return (', '.join('—' if v is None else f"{v:.1e} ({d:+.2g})" for v, d in zip(sp['main_err'], sp['main_dw'])) + ' | '
            + f"{sp['n_rest']}, {sp['rest_sw']:+.3f} ({sp['rest_w']:.3f}), mass {sp.get('rest_mass', 0.0):+.3f}"
            + (f" (window {sp['rest_win']:+.1e})" if 'rest_win' in sp else '') + f", {sp['rest_foot']:.1e}"
            + (f" (rms {sp['rest_rms']:.1e})" if 'rest_rms' in sp else '') + ", " + ('—' if sp['rest_sig'] is None else f"{sp['rest_sig'][0]:.3f}–{sp['rest_sig'][1]:.3f}") + ', '
            + ('—' if sp['rest_c'] is None else f"{sp['rest_c']:.2f}"))


def end_split(r):
    """`fit_split` of a path's end (H, HP: r['moving'], its formed kernels)."""
    m = r['moving']
    w, th = np.asarray(m['w'], float), np.atleast_2d(np.asarray(m['theta'], float)).reshape(-1, 2)
    fm = np.abs(w) > 1e-2 * min(abs(a) for a in r['As'])
    return fit_split(th[fm], w[fm], r['cs'], r['ss'], r['As'])


def hll_job(args):
    """HLL (referee rounds 74, 80): one of HL's cells continued to TAUS_HLL[-1]: at TAUS_HLL the formed kernels' states, the largest
    |w|, the residual, and the fit's split (`fit_split`: per feature its main kernel, its largest relative error in (c, σ, w) and its
    signed w/A − 1; the other formed kernels' Σw, Σ|w|, their sum's largest value on the samples over max f*, their widths and their
    centres' largest |c|)."""
    name, eps = args
    _, cs, ss, As = next(x for x in H_TARGETS if x[0] == name)
    t0 = time.time()
    L, K = law(), MW.BumpsCW(XS, S0)
    f = feature_sum(cs, ss, As)
    _, W, C, _, _ = RN.loop_nd(L, K, f, theta0(), eps, np.r_[0.0, TAUS_HLL], rtol=1e-9, atol=1e-12)
    out = []
    for t, tau in enumerate(TAUS_HLL, start=1):
        big = np.flatnonzero(np.abs(W[t]) > 1e-2 * min(abs(a) for a in As))
        th, w = C[t][big], W[t][big]
        out.append(dict(tau=float(tau), theta=th.tolist(), w=w.tolist(), wmax=float(np.abs(W[t]).max()),
                        rms=float(np.sqrt(np.mean((f - K.synth(C[t], W[t])) ** 2))), **fit_split(th, w, cs, ss, As)))
    return dict(name=name, eps=eps, loops=out, secs=time.time() - t0)


def part_HLL():
    from multiprocessing import Pool
    with Pool(N_WORKERS) as pool:
        return pool.map(hll_job, list(HLL), chunksize=1)


def loop_row(name, eps, x):
    """HL's and HLL's row (referee round 74): the formed kernels, the visible ones (`visibility` ≥ VIS) of the other sign (the most
    negative w), the invisible formed ones, the smallest visible width, those centred outside the samples' window, the largest |w|, the
    residual."""
    w, th = np.asarray(x['w'], float), np.atleast_2d(np.asarray(x['theta'], float)).reshape(-1, 2)
    vis = visibility(th) >= VIS if w.size else np.zeros(0, bool)
    neg = (w < 0) & vis
    return (f"| {name} | {eps:g} | {x['tau']:g} | {w.size} | {int(neg.sum())}" + (f" ({w[neg].min():+.3f})" if neg.any() else '')
            + f" | {int((~vis).sum())} | " + (f"{np.abs(th[vis, 1]).min():.3f}" if vis.any() else '—')
            + f" | {int((np.abs(th[:, 0]) > XS[-1]).sum())} | {x['wmax']:.3g} | {x['rms']:.1e} |")


def report(A=None, B=None, C=None, D=None, Ee=None, F=None, G=None, H=None, H2=None, H3=None, HL=None, HP=None, GL=None, GT=None,
           A2=None, HLL=None, SQ=None, A2E=None, A2P=None, A2H=None, A2R=None, A2T=None):
    N = XS.size
    L_ = ['# Proposition D‴.59: the transit with moving widths (`exp_moving_width.py`)', '',
          f"The sliding path with moving widths (θ = (c, σ), `moving_width.BumpsCW`) and with the widths fixed at σ₀ = {S0:g} "
          f"(`rich_transit_nd.BumpsND`); {C0.size} kernels started at −1, −0.9, …, 1; {N} samples of [{XS[0]:g}, {XS[-1]:g}]; exact: as many "
          f"formed kernels (|w| above 1% of the smallest amplitude) as features, each feature's centre, width and amplitude within "
          f"{EXACT_TOL:g}, the residual's rms below 1e-08.", '']
    if A is not None:
        L_ += ['## A. One feature', '',
               f"A = {A_A:g} at c* = {C_A:g}, the width σ* = ratio·σ₀; moving widths: the kernels, the largest error of (c, σ, w), the "
               f"residual; fixed widths (up to {MAX_ENT_FIXED} entries, {MAX_RHS_FIXED:.0e} evaluations — † the budget spent first): the kernels and the residual:", '',
               '| σ*/σ₀ | moving: kernels | exact | max error (c, σ, w) | residual | fixed: kernels | residual |', '|---|---|---|---|---|---|---|']
        for r in A:
            m, f_ = r['moving'], r['fixed']
            L_.append(f"| {r['ratio']:.4g} | {m['n']} | {'yes' if m['exact'] else 'no'} | {', '.join(f'{v:.1e}' for v in m['err'])} | "
                      f"{m['rms']:.1e} | {f_['n']}{'†' if f_.get('oob') else ''} | {f_['rms']:.2e} |")
        L_.append('')
    L_ += a2_lines(A2, A2E, A2P, A2H, A2R, A2T)
    if B is not None:
        L_ += ['## B. The one-kernel reduced energy (the continuum)', '',
               "J(c, σ) = ½‖f*‖² − ½⟨f*, φ_{c,σ}⟩²/‖φ_{c,σ}‖² on a 241 × 241 grid (c within ±6σ*, σ from σ*/8 to 8σ*, A = 1, c* = 0): its "
               "grid-local minima and the closed-form gradient's largest relative error against central differences (20 random points):", '',
               '| σ*/σ₀ | local minima (c, σ/σ*) | min J | gradient error |', '|---|---|---|---|']
        for r in B:
            ss = r['ratio'] * S0
            L_.append(f"| {r['ratio']:.4g} | " + ', '.join(f"({c:+.1e}, {s / ss:.3f})" for c, s in r['minima']) + f" | {r['J_min']:.1e} | "
                      f"{r['grad_err']:.1e} |")
        L_.append('')
    def oc(r):
        return outcome(r['moving'], r['As'], r['cs'], r['ss'])

    def split_text(m, As):
        w = np.asarray(m['w'], float)
        fj = np.flatnonzero(np.abs(w) > 1e-2 * np.abs(np.asarray(As)).min())
        th = np.atleast_2d(m['theta'])[fj]
        out = []
        for i in range(fj.size):
            for j in range(i + 1, fj.size):
                if pairs_of(th[[i, j]], w[fj][[i, j]]):
                    dc, ds, pr = width_split(th, w[fj], i, j)
                    out.append(f"{w[fj][i]:+.3f}/{w[fj][j]:+.3f}, Δc {dc:.2e}, the width split {ds:+.4e} (predicted {pr:+.4e})")
        return '; '.join(out)
    if C is not None:
        L_ += ['## C. Random three-feature targets', '',
               f"{N_C} targets per minimal separation (seed {SEED_C}; widths log-uniform in [{WID_C[0]:g}, {WID_C[1]:g}], amplitudes ±U(0.4, 1); "
               "a moving-width path stopped by its budget counted at its stalled block — G — not at the block before it, referee round 68):", '',
               '| Δ ≥ | moving widths: exact | formed kernels (median [range]) | residual (median [max]) | fixed widths: kernels (median '
               '[range]) | residual (median [max]) |', '|---|---|---|---|---|---|']
        for d in DELTAS_C:
            rs = [r for r in C if r['dmin'] == d]
            km, kf = [moving_state(r)[0] for r in rs], [r['fixed']['n'] for r in rs]
            em, ef = [moving_state(r)[1] for r in rs], [r['fixed']['rms'] for r in rs]
            L_.append(f"| {d:g} | {sum(oc(r) == 'exact' for r in rs)} of {len(rs)} | {np.median(km):g} [{min(km)}, {max(km)}] | "
                      f"{np.median(em):.1e} [{max(em):.1e}] | {np.median(kf):g} [{min(kf)}, {max(kf)}] ({sum(r['fixed'].get('oob', False) for r in rs)}† "
                      f"stopped by the budget) | {np.median(ef):.1e} [{max(ef):.1e}] |")
        cls = {}
        for r in C:
            cls.setdefault((r['dmin'], oc(r)), []).append(r)
        L_ += ['', "The moving widths' outcomes (`outcome`: the formed kernels |w| > 1% of the smallest amplitude; exact with a null "
               "cluster where each feature's nearest kernel is exact and the others sum to nothing; a coinciding pair within "
               f"{PAIR_TOL:g} in |Δc|/σ + |Δlog σ| whose net coefficient and centre of mass fit a feature within {LAND_TOL:g}, a coinciding "
               "pair off the features where it fits none (referee round 106); a budget stop a dipole run-off if the stalled block has an opposite-sign pair under 0.05 "
               "with |w| beyond 10 times the largest amplitude): " + '; '.join(
            f"Δ ≥ {d:g}: " + ', '.join(f"{len(cls.get((d, o), []))} {o}" for o in ('exact', 'exact, null cluster', 'coinciding pair', 'coinciding pair off the features', 'cancellation', 'dipole run-off', 'stopped by the budget', 'not fitted') if cls.get((d, o)) or o == 'exact')
            for d in DELTAS_C) + '.']
        bad = [r for r in C if oc(r) != 'exact']

        def stall(m):
            st = m.get('stalled')
            if not m.get('oob') or st is None:
                return ''
            cp = closest_pair(st['theta'], st['w'], opposite=True)
            return (f"; the stalled block: {len(st['w'])} kernels, max|w| {np.abs(st['w']).max():.4g} (the accepted polish's, "
                    "solver-dependent)" + ('' if cp is None else f", its closest opposite pair {cp[0]:.2e} apart"))
        L_ += ['', "Not exact with moving widths (the end's state; for a budget stop the state before the stalled block, and the "
               "stalled block): " + ('; '.join(
            f"Δ ≥ {r['dmin']:g}, target {r['q']} (widths {', '.join(f'{v:.3f}' for v in r['ss'])}; centres "
            f"{', '.join(f'{v:+.2f}' for v in r['cs'])}; amplitudes {', '.join(f'{v:+.2f}' for v in r['As'])}): "
            f"{oc(r)}, {r['moving']['n']} kernels"
            + (f" (the pair: {split_text(r['moving'], r['As'])})" if split_text(r['moving'], r['As']) else '')
            + f", residual {r['moving']['rms']:.1e}" + stall(r['moving']) for r in bad) or 'none') + '. The width split '
            "(σ_j − σ_i = (w_j − w_i)Δc²/(6Wσ), the pair's third moment cancelled) against the observed.", '']
    if D is not None:
        L_ += ['## D. A narrow feature on a wide pedestal', '',
               "The pedestal 0.8φ(·; 0, 0.3) and the narrow feature 0.5φ(·; d, 0.06):", '',
               '| d | moving: formed kernels | outcome | max error (c, σ, w) | residual | fixed: kernels | residual |', '|---|---|---|---|---|---|---|']
        for r in D:
            m, f_ = r['moving'], r['fixed']
            As_ = [0.8, 0.5]
            sp = split_text(m, As_)
            oc_ = outcome(m, As_, [0.0, r['d']], [0.3, 0.06]) + (f" ({sp})" if sp else '')
            L_.append(f"| {r['d']:g} | {n_formed(m, As_)} | {oc_} | {', '.join(f'{v:.1e}' for v in m['err'])} | {m['rms']:.1e} | "
                      f"{f_['n']}{'†' if f_.get('oob') else ''} | {f_['rms']:.2e} |")
        L_.append('')
    if Ee is not None:
        L_ += ['## E. Under noise: the discrepancy stop', '',
               f"C's first {N_E} targets at Δ ≥ {DELTAS_C[1]:g}, one draw each, the path with moving widths to {K_C + 2} entries (past "
               "the features a noise kernel can shrink onto a sample: the fit has no minimiser); the first block whose residual's rms (against "
               f"the data) is within (1 + {DELTA_E:g})s; its kernels and its error against the truth (rms on the samples), median [range] and "
               "the rms over the draws, against the oracle's s(3K/N)^{1/2} (three parameters per kernel, K = 3, the rms of the errors), and "
               "per draw against the oracle fit (the K-kernel least-squares fit of the data started at the truth):", '',
               '| s | stopped | kernels | error: median [range]; rms | the oracle | the stopped error over the oracle fit\'s: median [range] '
               '(within 10⁻⁴ of it: how many, the largest deviation among them) '
               '| not stopped: the last block\'s rms |', '|---|---|---|---|---|---|---|']
        for s in LEVELS_E:
            rs = [r for r in Ee if r['s'] == s]
            st = [(r['moving'], r.get('oracle')) for r in rs if r['moving']['stopped']]
            ks, es = [x['k'] for x, o in st], [x['err'] for x, o in st]
            qo = [x['err'] / o for x, o in st if o]
            ns = [r['moving'] for r in rs if not r['moving']['stopped']]
            L_.append(f"| {s:g} | {len(st)} of {len(rs)} | " + (f"{np.median(ks):g} [{min(ks)}, {max(ks)}] | {np.median(es):.2e} "
                      f"[{min(es):.2e}, {max(es):.2e}]; {np.sqrt(np.mean(np.square(es))):.3e}" if st else '— | —') + f" | "
                      f"{s * np.sqrt(3 * K_C / N):.3e} | " + (f"{np.median(qo):.4f} [{min(qo):.4f}, {max(qo):.4f}] "
                                                              f"({sum(abs(q - 1) < 1e-4 for q in qo)} of {len(qo)}" + (f", {max(abs(q - 1) for q in qo if abs(q - 1) < 1e-4):.0e})" if any(abs(q - 1) < 1e-4 for q in qo) else ')') if qo else '—') + ' | '
                      + (', '.join(f"{x.get('rms', np.nan):.2e}" + (' (budget)' if x.get('oob') else '') for x in ns) or '—') + ' |')
        L_.append('')
    if F is not None:
        L_ += ['## F. The real loop at a small mass', '',
               f"The free law with moving centres and widths from α = 0, masses ε = {', '.join(f'{e:g}' for e in EPS_F)} (`loop_nd` on `BumpsCW`, rtol 10⁻⁹), to "
               f"{', '.join(f'{t:g}' for t in TAUS_F)} times the path's last entry, on C's first {N_F} targets at Δ ≥ {DELTAS_C[1]:g}: its formed "
               "kernels (|w| above 1% of the smallest amplitude) against the path's, the path's kernels' largest deviation, the others' largest "
               "|w|, the closest pair among the formed kernels (|Δc|/σ + |Δlog σ|; its |Δc|, |Δσ|) and the residual; the path's closest pair "
               "and residual:", '',
               "| target | the path's kernels (closest pair; residual) | ε | times | formed by the loop | max \\|Δθ\\| | max \\|Δw\\| | the others' "
               "largest \\|w\\| | the loop's closest pair: distance (\\|Δc\\|, \\|Δσ\\|; coefficients) | residual |", '|---|---|---|---|---|---|---|---|---|---|']
        for r in F:
            for x in r.get('snaps', [dict(times=None, formed=r['formed'], dtheta=r['dtheta'], dw=r['dw'], rest=r['rest'], pair=None,
                                          rms=np.nan)]):
                pp = r.get('path_pair')
                pa = x.get('pair')
                L_.append(f"| {r['q']} | {r['path']} (" + ('—' if pp is None else f"{pp:.1e}") + f"; {r.get('path_rms', np.nan):.1e}) | "
                          + f"{x.get('eps', EPS_F[0]):g} | " + ('—' if x['times'] is None else f"{x['times']:g}") + f" | {x['formed']} | "
                          f"{x['dtheta']:.1e} | {x['dw']:.1e} | {x['rest']:.1e} | "
                          + ('—' if pa is None else f"{pa['d']:.1e} ({pa['dc']:.1e}, {pa['ds']:.1e}"
                             + (f"; {pa['w'][0]:+.3f}/{pa['w'][1]:+.3f}" if 'w' in pa else '') + ")") + f" | {x['rms']:.1e} |")
        L_.append('')
    if G is not None:
        L_ += ['## G. The blocks that never settle: a dipole', '',
               "C's paths stopped by their budget. The stalled block (C's path again; its state when the budget ran out): its kernels, the "
               "largest |w|, the closest opposite-sign pair's distance (|Δc|/σ + |Δlog σ|), its residual against the optimal dipole fit's "
               "(below), its polishes (status/evaluations/accepted, the largest coefficient over the block's). The block from its entry, "
               f"flowed on the critical curve to s = {S_G[-1]:g}: the largest |w| and the pair's distance at s = {', '.join(f'{v:g}' for v in S_G)}; "
               "over the last factor in s the exponents of the largest |w|, the distance, the pair's |Δc| and |Δσ| and the residual's "
               "excess over the dipole infimum; the dipole law's constant K = v²(E − E_inf) at the last two outputs (v the pair's "
               "half-difference |w_a − w_b|/2, E_inf the optimal dipole fit's — the other kernels, Wφ(θ̄) and p·∇_θφ(θ̄) by least squares), "
               "the growth of v³ over the last interval against 3μ₀K (μ₀ = 4Gc/A: the coefficient's mobility μ₀|w| on the exponential "
               "branch), the moment's σ-share |p_σ|/|p|:", '',
               '| Δ ≥ | target | entries | the stalled block: kernels, max\\|w\\|, distance, residual (infimum), polishes | the flow: largest '
               '\\|w\\| | distance | exponents: \\|w\\|, distance, \\|Δc\\|, \\|Δσ\\|, the residual\'s excess | K | Δv³/Δs, 3μ₀K | \\|p_σ\\|/\\|p\\| |',
               '|---|---|---|---|---|---|---|---|---|---|']
        for r in G:
            fl, lw, st = r['flow'], r.get('law') or [], r.get('stalled')
            ls = np.log(fl[-1]['s'] / fl[-2]['s'])
            ex = dipole_exponents(fl)
            stc = '—'
            if st is not None:
                pol = ' '.join(f"{x['status']}/{x['nfev']}/{'a' if x['accepted'] else 'r'}/{x['wmax_ratio']:.3g}" for x in st['polish'])
                stc = (f"{st['n']}, {st['wmax']:.4g}, " + ('—' if st['sep'] is None else f"{st['sep']:.2e}") + f", {st['rms']:.6e}"
                       + (f" ({lw[-1]['rms_inf']:.6e})" if lw else '') + f", {pol or '—'}")
            exs, kk, vv, ps = '—', '—', '—', '—'
            if ex is not None and len(lw) == 2 and lw[0] is not None and lw[1] is not None:
                a_, b_ = lw
                exs = (f"{ex[0]:+.4f}, {ex[1]:+.4f}, {np.log(b_['dc'] / a_['dc']) / ls:+.4f}, {np.log(b_['ds'] / a_['ds']) / ls:+.4f}, "
                       f"{np.log((fl[-1]['rms'] - b_['rms_inf']) / (fl[-2]['rms'] - b_['rms_inf'])) / ls:+.4f}")
                K2 = [x['v'] ** 2 * (x['E'] - x['E_inf']) for x in lw]
                kk = ', '.join(f"{k_:.5e}" for k_ in K2)
                vv = f"{(b_['v'] ** 3 - a_['v'] ** 3) / (fl[-1]['s'] - fl[-2]['s']):.5e}, {3 * b_['mu0'] * K2[1]:.5e}"
                ps = f"{abs(b_['p'][1]) / np.hypot(*b_['p']):.3f}"
            L_.append(f"| {r['dmin']:g} | {r['q']} | {r['entries']} | {stc} | {', '.join(format(x['wmax'], '.3g') for x in fl)} | "
                      + ', '.join('—' if x['sep'] is None else f"{x['sep']:.2e}" for x in fl) + f" | {exs} | {kk} | {vv} | {ps} |")
        L_ += g_loop_lines(G)
        L_ += g_fine_lines(G, GT, SQ)
    if GL is not None:
        L_ += gl_lines(GL)
    if H is not None:
        L_ += ['## H. Two features of one sign, close', '',
               "Two features φ(·; c, σ*) of the given amplitudes (both positive) at the centres ±0.08 — 0.16 apart, 1.33 times the kernels' "
               f"starting width σ₀ = {S0:g} — unless stated, the moving-width path (as C's: its tolerance 1e-09, its budget "
               f"{MAX_RHS_MOVING:.0e}): the outcome (as C's), the most kernels of the other sign in a block along the path and the most "
               "negative coefficient, the end's formed kernels of the other sign and the smallest formed width — counting only kernels visible on "
               f"the samples (largest φ there at least {VIS:g}; referee round 74), the end's formed ones of the other sign including the "
               "invisible in brackets —, over all the end's active kernels those narrower than the samples' spacing and those of negative "
               "width (the widths are not constrained) and how many of them are formed, the formed kernels invisible on the samples, the "
               "largest |w|, the residual:", '']
        L_ += h_table(H)
    if H2 is not None:
        L_ += ['### H2. The σ* = 0.06 paths to a tolerance of ' + f"{TOL_H2:g}", '',
               "The same paths with the tolerance lowered (the 1e-09 ends were the tolerance's): the outcome, the entries and the last one's "
               "τ, the residuals of the last blocks, the formed kernels by the feature nearest each (their count and mass over the feature's, "
               "Σw_i|σ_i|/(Aσ*)), and "
               "H's columns:", '',
               '| target | outcome | entries (last τ) | the last blocks\' residuals | formed per feature (count: mass) | visible negative along '
               'the path (most negative w) | at the end (all formed) | smallest visible width | sub-sample, negative widths (all active; formed) | '
               'max\\|w\\| |', '|---|---|---|---|---|---|---|---|---|---|']
        for r in H2:
            L_.append(f"| {r['name']} | {r['outcome']} | {r['moving']['entries']} ({r['last_tau']:.2e}) | "
                      + ', '.join(f"{v:.2e}" for v in r['rms_blocks']) + ' | ' + ', '.join(f"{n}: {s_:.4f}" for n, s_ in r['clusters'])
                      + f" | {r['neg_path']} ({r['wmin_path']:+.3f}) | {r['neg_end']} ({r.get('neg_end_all', '—')}) | {r['smin']:.2e} | "
                      f"{r['sub']}, {r['negw']} ({r.get('sub_formed', '—')}, {r.get('negw_formed', '—')}) | {r['wmax']:.3g} |")
        if all('n_act' in r for r in H2):                               # referee round 68: why the paths ended
            L_ += ['', "Why they ended (referee rounds 68, 74): " + '; '.join(
                f"{r['name']}: the stop {r.get('stop', '—')}, {r['n_act']} of {C0.size} kernels active" + (' (the dictionary exhausted)' if r['exhausted'] else '')
                + f", {r['pol_capped']} of the {r['pol_acc']} accepted polishes stopped by their budget (status 0), {r['neg_sigma']} "
                f"active kernels of negative width, {r['invisible']} invisible on the samples (largest φ below 10⁻¹⁰)" for r in H2) + '.']
        L_.append('')
    if HP is not None:
        L_ += ['### HP. H\'s path without the polish (referee round 68)', '',
               "The moving-width path (tolerance 10⁻⁹) with the least-squares polish off, against H's (with it): the outcome, the "
               "negative kernels visible on the samples along the path and at the end (all formed at the end in brackets), the largest |w|, "
               "the residual, the active kernels and why the path stopped (referee round 74: tolerance, budget, dictionary exhausted), and the "
               "fit's split at the end (referee round 80; as HLL's):", '',
               '| target | outcome (H) | visible negative along the path, at the end (H) | max\\|w\\| (H) | residual (H) | active | stop | '
               'the main kernels\' errors (w/A − 1) | the others: n, Σw (Σ\\|w\\|), on the samples, widths, \\|c\\| max |',
               '|---|---|---|---|---|---|---|---|---|']
        Hn = {r['name']: r for r in (H if H is not None else [])}
        for r in HP:
            h = Hn.get(r['name'], {})
            L_.append(f"| {r['name']} | {r['outcome']} ({h.get('outcome', '—')}) | {r['neg_path']}, {r['neg_end']} [{r.get('neg_end_all', '—')}] "
                      f"({h.get('neg_path', '—')}, {h.get('neg_end', '—')} [{h.get('neg_end_all', '—')}]) | {r['wmax']:.3g} ("
                      + (f"{h['wmax']:.3g}" if 'wmax' in h else '—') + f") | {r['moving']['rms']:.1e} ("
                      + (f"{h['moving']['rms']:.1e}" if 'moving' in h else '—') + f") | {r['n_act']} | {r.get('stop', '—')} | "
                      f"{split_cells(end_split(r))} |")
        L_.append('')
    if H3 is not None:
        L_ += ['### H3. The separation scan', '',
               "σ* = 0.06, equal amplitudes 1, at ±c (2c in units of σ₀ in brackets; H's ±0.08 is 1.33σ₀), H's columns:", '']
        L_ += h_table(sorted(list(H3) + [r for r in (H if H is not None else []) if r['name'] == 'σ* 0.06, equal'],
                             key=lambda r: r['cs'][1]), sep=True)
    if HL is not None:
        L_ += ['### HL. The real loop on H\'s targets (referee round 68: the asymmetric ones, to τ = 10⁵)', '',
               "The free law from α = 0 at the masses ε (`loop_nd`, rtol 10⁻⁹) at the times τ: the formed kernels (|w| above 1% of the "
               f"smallest amplitude), those visible on the samples (largest φ there at least {VIS:g}; referee round 74) of the other sign, the "
               "invisible formed ones, the smallest visible width, the formed kernels centred outside the samples' window, the largest |w|, "
               "the residual, and the fit's split (referee round 80; as HLL's):", '',
               '| target | ε | τ | formed | visible of the other sign (most negative w) | invisible | smallest visible width | outside the window | '
               'max\\|w\\| | residual | the main kernels\' errors (w/A − 1) | the others: n, Σw (Σ\\|w\\|), on the samples, widths, \\|c\\| max |',
               '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for r in HL:
            L_ += [loop_row(r['name'], x['eps'], x) + f" {split_cells(fit_split(x['theta'], x['w'], r['cs'], r['ss'], r['As']))} |"
                   for x in r['loops']]
        L_.append('')
    if HLL is not None and len(HLL):
        L_ += ['### HLL. HL\'s cells to τ = 10⁹ (referee rounds 74, 80, 83)', '', "HL's columns, and the fit's split (`fit_split`, referee "
               "round 83: the main kernels the assignment minimizing the largest error): per feature its main kernel's largest relative "
               "error in (c, σ, w) and its signed w/A − 1 (inflated: a cancellation; deficient: a split); the other formed kernels — their "
               "number, Σw and Σ\\|w\\|, their sum's largest value on the samples over max f* (a null cluster where it is small against the "
               "residual), their widths and their centres' largest \\|c\\|:", '',
               '| target | ε | τ | formed | visible of the other sign (most negative w) | invisible | smallest visible width | outside the window | '
               'max\\|w\\| | residual | the main kernels\' errors (w/A − 1) | the others: n, Σw (Σ\\|w\\|), on the samples, widths, \\|c\\| max |',
               '|---|---|---|---|---|---|---|---|---|---|---|---|']
        ht = {x[0]: x[1:] for x in H_TARGETS}
        nc = []
        for r in HLL:
            sps = [fit_split(x['theta'], x['w'], *ht[r['name']]) for x in r['loops']]
            for x, sp in zip(r['loops'], sps):
                L_.append(loop_row(r['name'], r['eps'], x) + f" {split_cells(sp)} |")
            ft, mf = [sp['rest_foot'] for sp in sps], float(np.abs(feature_sum(*ht[r['name']])).max())
            rr = [sp['rest_rms'] / (x['rms'] / mf) for sp, x in zip(sps, r['loops']) if x['rms'] > 0]   # referee round 87: one norm
            fr = [f"{a / b:.3g}" for a, b in zip(ft[:-1], ft[1:]) if b > 0]
            nc.append(f"{r['name']} at {r['eps']:g}: " + (f"×{', '.join(fr)}; " if fr else '')
                      + f"at τ = {r['loops'][-1]['tau']:g} {ft[-1]:.1e} against {r['loops'][-1]['rms'] / mf:.1e}"
                      + (f"; the others' rms over the residual's " + ', '.join(f"{v:.3g}" for v in rr) if rr else '')
                      + (f" (1% of max f* after {np.log(ft[-1] / 0.01) / np.log(ft[-2] / ft[-1]):.1f} more decades at the last decade's "
                         "rate)" if len(ft) > 1 and ft[-1] > 0.01 and ft[-2] > ft[-1] else ''))
        L_ += ['', "The others' sum on the samples (referee round 83): its fall per decade of τ, and at the last τ against the residual's "
               f"rms, both over max f* (a null cluster where, the main kernels exact, the others' rms is below {NULL_TOL:g} of max f*: "
               "referee round 97): " + '; '.join(nc) + '.']
        L_.append('')
    L_ += split_summary(H, HP, H3, HL, HLL)
    L_ += ['Numbers are produced by `experiments/exp_moving_width.py`; nothing is edited by hand.']
    return '\n'.join(L_) + '\n'


def v3_slopes(dp, rate):
    """The slopes of v³ between successive (τ, v) over rate."""
    return [(b[1] ** 3 - a[1] ** 3) / (b[0] - a[0]) / rate for a, b in zip(dp[:-1], dp[1:])]


def kappa_loop(ph, td, t3p=None):
    """The loop's κ (referee round 80): v_land over v(τ₃⁻), the loop's v at its last fine-grid time t_v3 before τ₃ carried by the law
    v³ ∝ τ − τ_d to the path's τ₃ t3p — v(t_v3) itself lags v(τ₃⁻) by the grid's offset ((τ₃ − τ_d)/(t_v3 − τ_d))^{1/3} —; without
    t3p, v_land/v(t_v3). Returns (κ, the offset) or (None, None)."""
    if not ph.get('v_land') or not ph.get('v3'):
        return None, None
    tv = ph.get('t_v3')
    off = ((t3p - td) / (tv - td)) ** (1 / 3) if t3p is not None and tv is not None and tv > td and t3p > td else 1.0
    return ph['v_land'] / (ph['v3'] * off), off


def g_fine_lines(G, GT=None, SQ=None):
    """G's phases on the fine grid (referee rounds 68, 74): per target and mass, τ_d, τ₃, τ₄, the dipole phase's slopes of v³ over the
    law (median), v before τ₃ over the law's at that time, the intermediate phase (formed, offset, (W − A)/|A|; its slopes of v³ over
    3μ₀K′Λ, K′ the path's), the first landing (τ, κ = v_land/v(τ₃⁻) by `kappa_loop`, p); per target the path's predictions
    (`kappa_from_path`, GT) against the loop's; κ against the mass; over every state whether p ≤ P_LAND and the landing agree."""
    if not any(r.get('fine') for r in G):
        return []
    gt = {x['q']: x for x in (GT if GT is not None else [])}
    out = ['', "The phases on a fine grid (τ = 10 … 140 by 0.1, with G's times; referee rounds 68, 74): τ₃ the first formation past the "
           "block's entry τ_d, τ₄ the next; over (τ_d, τ₃) the slopes of v³ over 3μ₀KΛ(ε) (median); v at the last time before τ₃ over the "
           "law's (3μ₀KΛ(ε)(t − τ_d))^{1/3} at that time t; at (τ₃ + τ₄)/2 the formed kernels, the dipole pair's offset from the nearest "
           "feature (σ*) and (W − A)/|A|, and over [τ₃, τ₄) the slopes of v³ over 3μ₀K′Λ(ε), K′ the three-kernel block's from the path "
           f"(below; median); the first landing past τ₃ (centre of mass and W within {LAND_TOL:g}): its τ, κ = v_land/v(τ₃⁻) with "
           "v(τ₃⁻) the loop's v at the last time before τ₃, κ′ = v_land over the law's v at the path's τ₃ (the loop's v before τ₃ "
           "carried by v³ ∝ τ − τ_d: `kappa_loop`, referee rounds 80, 83 — at larger masses not the loop's own κ, whose crossing is "
           "the loop's), and p = |w_a − w_b||Δθ|/(|W|σ̄); † another sequence — "
           + ("the loop's formed kernels at its τ₃ or τ₄ not the path's, by index (SQ, referee round 83)" if any(
               gt.get(x['q'], {}).get('block') is not None for x in (SQ if SQ is not None else []))
              else f"the intermediate slope beyond 1 ± {SEQ_TOL:g} times the path's law") + " —, its κ not the path's:", '',
           '| target | ε | τ_d | τ₃ | τ₄ | dipole: slope | v(τ₃⁻)/law | intermediate: formed, offset, (W − A)/\\|A\\|, slope | '
           "landing: τ, κ, κ′, p |", '|---|---|---|---|---|---|---|---|---|']
    agree, agree4, agree7, kap, offs = [0, 0], [0, 0], [0, 0], {}, {}
    sq = {x['q']: x['seq'] for x in (SQ if SQ is not None else []) if gt.get(x['q'], {}).get('block') is not None}   # SQ needs GT's
    sqr, prox, brk = {}, {}, []
    for r in G:
        lw = r.get('law') or []
        if not r.get('fine') or not lw or lw[-1] is None:
            continue
        Kc, mu0 = lw[-1]['v'] ** 2 * (lw[-1]['E'] - lw[-1]['E_inf']), lw[-1]['mu0']
        Kp = gt.get(r['q'], {}).get('Kp')
        Lam = {x['eps']: x['Lam'] for x in r['loops']}
        for eps, ph in sorted(r['fine'].items(), key=lambda t: -t[0]):
            for acc, a_ in ((agree, ph.get('agree')), (agree4, ph.get('agree4')), (agree7, ph.get('agree') if eps <= 1e-7 else None)):
                if a_:
                    acc[0] += a_[0]
                    acc[1] += a_[1]
            lam = Lam.get(eps, RT.transit_scale(law(), eps))              # `loop_nd`'s Λ(ε) (the masses without G's tables)
            rate = 3 * mu0 * Kc * lam
            sl = v3_slopes(ph['dip'], rate)
            tv = ph.get('t_v3', ph['t3'])
            vlaw = (rate * (tv - r['tau_d'])) ** (1 / 3) if tv is not None and tv > r['tau_d'] else None
            sl2 = v3_slopes(ph.get('dip2') or [], 3 * mu0 * Kp * lam) if Kp else []
            k_, of_ = kappa_loop(ph, r['tau_d'], gt.get(r['q'], {}).get('t3_path'))
            k0 = kappa_loop(ph, r['tau_d'])[0]
            seq = bool(sl2) and abs(np.median(sl2) - 1) > SEQ_TOL         # referee round 80: the loop on another sequence (a proxy)
            prox.setdefault(r['q'], {})[eps] = seq
            if r['q'] in sq and eps in sq[r['q']]:                       # referee round 83: by index
                sqr.setdefault(r['q'], {})[eps] = c_ = sq_check(sq[r['q']][eps], r['tau_d'], r['entries'], gt[r['q']])
                seq = not c_[4]
            kap.setdefault(r['q'], {})[eps] = (k_, seq)
            offs.setdefault(r['q'], {})[eps] = of_
            it = ph['inter']
            out.append(f"| {r['q']} | {eps:g} | {r['tau_d']:.4g} | " + ('—' if ph['t3'] is None else f"{ph['t3']:g}") + ' | '
                       + ('—' if ph['t4'] is None else f"{ph['t4']:g}") + ' | ' + (f"{np.median(sl):.4f}" if sl else '—') + ' | '
                       + (f"{ph['v3'] / vlaw:.4f}" if ph['v3'] and vlaw else '—') + ' | '
                       + ('—' if it is None else f"{it['n']}, {it['off']:.3g}, {it['dW']:+.3f}, " + (f"{np.median(sl2):.4f}" if sl2 else '—'))
                       + ' | ' + ('—' if ph['t_land'] is None else f"{ph['t_land']:g}, " + ('—' if k0 is None else f"{k0:.3f}") + ', '
                                  + ('—' if k_ is None else f"{k_:.3f}") + ('†' if seq else '') + f", {ph['p_land']:.3f}") + ' |')
    if gt:
        out += ['', "The intermediate phase from the path (`kappa_from_path`, referee rounds 74, 80): K from the third order, −⟨r_dip, "
                "∇³φ:p³⟩/24 at the block's optimal dipole fit, beside G's v²(E − E_inf); τ₃ from the duals; K′/K, K′ the three-kernel block's "
                f"(the τ₃ entrant added), with that block's closure rms (below {CLOSURE_EXACT:g} the block fits exactly: K′ = 0, no τ₄, "
                "κ = 1); τ₄ from the duals moved on by that block's residual; for K′ < 0 the dissolution time τ₃ + K(τ₃ − τ_d)/|K′|, at which "
                "v³ would reach 0 (the four stages need τ₄ before it); κ = (1 + K′(τ₄ − τ₃)/(K(τ₃ − τ_d)))^{1/3}; the landing's v by the law "
                "v³ = 3μ₀Λ(ε)(K(τ₃ − τ_d) + K′(τ₄ − τ₃)) with the path's τ₃ and τ₄ — each against the loop's at the smallest mass (its κ as in "
                "the table above) —; and v at τ = 6400 over v_land at ε = 10⁻¹⁰, minus 1 (v after τ₄):", '',
                '| target | K: third order (G\'s) | τ₃: path (loop) | K′/K (the closure\'s rms) | τ₄: path (loop) | dissolution | '
                'κ: path (loop: κ′; its own crossing\'s bracket) | v_land: law (loop) | v(6400)/v_land − 1 |', '|---|---|---|---|---|---|---|---|---|']
        for r in G:
            g_, lw = gt.get(r['q']), r.get('law') or []
            if not g_ or not r.get('fine') or not lw or lw[-1] is None or g_.get('K_third') is None:
                continue
            Kc, mu0, e0 = lw[-1]['v'] ** 2 * (lw[-1]['E'] - lw[-1]['E_inf']), lw[-1]['mu0'], min(r['fine'])
            ph = r['fine'][e0]
            Lam = {x['eps']: x['Lam'] for x in r['loops']}
            lam0 = Lam.get(e0, RT.transit_scale(law(), e0))
            t3p, t4p, K3, Kp, ex = g_['t3_path'], g_.get('t4_path'), g_['K_third'], g_.get('Kp'), bool(g_.get('exact'))
            br = K3 * (t3p - r['tau_d']) + (0.0 if ex else Kp * (t4p - t3p) if t4p is not None and Kp is not None else np.nan)
            vl = (3 * mu0 * lam0 * br) ** (1 / 3) if np.isfinite(br) and br > 0 else None
            p10 = r['fine'].get(1e-10)
            x64 = next((x for x in r['loops'] if x['eps'] == 1e-10 and x['tau'] == 6400.0), None)
            dp64 = dipole_pair(x64['theta'], x64['w']) if x64 is not None else None
            fr = dp64[2] / p10['v_land'] - 1 if dp64 is not None and p10 and p10['v_land'] else None
            k_loop = kappa_loop(ph, r['tau_d'], t3p)[0]
            br_ = (kappa_loop(ph, r['tau_d'], ph['t3'])[0], kappa_loop(ph, r['tau_d'])[0]) if ph['t3'] is not None else (None, None)
            if g_.get('kappa_path') is not None and None not in br_:
                brk.append((r['q'], g_['kappa_path'], br_[0], br_[1]))
            rm = '' if g_.get('rms2') is None else f" ({g_['rms2']:.1e})"
            kpc = f"0, exact{rm}" if ex else ('—' if Kp is None else f"{Kp / K3:+.4f}{rm}")
            t4c = '— (exact)' if ex else ('—' if t4p is None else f"{t4p:.4g}")
            out.append(f"| {r['q']} | {K3:.5e} ({Kc:.5e}) | {t3p:.4g} (" + ('—' if ph['t3'] is None else f"{ph['t3']:g}") + f") | {kpc} | "
                       + t4c + ' (' + ('—' if ph['t4'] is None else f"{ph['t4']:g}") + ') | '
                       + ('—' if g_.get('t_diss') is None else f"{g_['t_diss']:.5g}") + ' | '
                       + ('—' if g_.get('kappa_path') is None else f"{g_['kappa_path']:.4f}")
                       + ' (' + ('—' if k_loop is None else f"{k_loop:.4f}") + ('' if None in br_ else f"; [{br_[0]:.4f}, {br_[1]:.4f}]")
                       + ') | ' + ('—' if vl is None else f"{vl:.4g}") + ' ('
                       + ('—' if not ph['v_land'] else f"{ph['v_land']:.4g}") + ') | ' + ('—' if fr is None else f"{fr:+.1e}") + ' |')
    if brk:                                                             # referee round 83: the path's κ against the loop's own
        ins = [str(q) for q, k, lo, hi in brk if lo <= k <= hi]
        outs = [f"{q} (the bracket {lo / k - 1:+.2%} … {hi / k - 1:+.2%} off it)" for q, k, lo, hi in brk if not lo <= k <= hi]
        out += ['', "The path's κ against the loop's own crossing's bracket at the smallest mass (referee round 83): inside on "
                + (', '.join(ins) or 'none') + ('; outside on ' + ', '.join(outs) if outs else '') + '.']
    out += ['', "κ′ against the mass (ε = " + ', '.join(f"{e:g}" for e in sorted({e for d in kap.values() for e in d}, reverse=True))
            + "; † another sequence): " + '; '.join(f"{q}: " + ', '.join(('—' if d[e][0] is None else f"{d[e][0]:.3f}") + ('†' if d[e][1] else '')
                                                                          for e in sorted(d, reverse=True)) for q, d in kap.items())
            + ". The grid's offset ((τ₃ − τ_d)/(t_v3 − τ_d))^{1/3} − 1 that κ′ removes, at ε = " + f"{min(e for d in offs.values() for e in d):g}: "
            + ', '.join(f"{q}: " + ('—' if d[min(d)] is None else f"{d[min(d)] - 1:+.2%}") for q, d in offs.items())
            + f". Over every fine-grid state with a dipole pair and p < 1, p ≤ {P_LAND:g} and the landing agree in {agree[0]} of {agree[1]} "
            f"({agree7[0]} of {agree7[1]} at ε ≤ 10⁻⁷; past τ₄ {agree4[0]} of {agree4[1]}).", '']
    rt = []
    for q, d in kap.items():                                            # referee round 83: κ′'s convergence rate on the path's sequence
        kp_ = gt.get(q, {}).get('kappa_path')
        es = [e for e in sorted(d, reverse=True) if d[e][0] is not None and not d[e][1]]
        if kp_ is None or len(es) < 2:
            continue
        dd = {e: d[e][0] - kp_ for e in es}
        fr = [abs(dd[a]) / abs(dd[b]) for a, b in zip(es[:-1], es[1:]) if abs(np.log10(a / b) - 1) < 1e-9 and min(abs(dd[a]), abs(dd[b])) >= 1e-3]
        rr = '' if not fr else (f"×{min(fr):.2f}, as ε^{{{np.log10(min(fr)):.2f}}}" if f"{min(fr):.2f}" == f"{max(fr):.2f}" else
                                f"×{min(fr):.2f}–{max(fr):.2f}, as ε^{{{np.log10(min(fr)):.2f}–{np.log10(max(fr)):.2f}}}")
        rt.append(f"{q}: " + ', '.join(f"{dd[e]:+.4f}" for e in es) + (f" ({rr})" if rr else '')
                  + " — over ε^{1/3}: " + ', '.join(f"{dd[e] / e ** (1 / 3):+.1f}" for e in es))   # referee round 87 (D3ch-7)
    if rt:
        out += ["κ′ − κ (the path's) on the path's sequence (referee round 83; ε as above, the masses on another sequence left out), "
                "and |κ′ − κ|'s fall per decade of ε where both are at least 10⁻³ — O(ε^{1/3}) is ×10^{1/3} = ×2.154 —: " + '; '.join(rt)
                + '.', '']
    if sqr:
        tot = sum(len(d) for d in sqr.values())
        on = sum(c_[4] for d in sqr.values() for c_ in d.values())
        same = sum((c_[0], c_[1]) == (r_['fine'][e]['t3'], r_['fine'][e]['t4']) for r_ in G if r_['q'] in sqr for e, c_ in sqr[r_['q']].items())
        fs = lambda t, s_: '—' if t is None else f"{t:g} {{" + ', '.join(map(str, s_)) + '}'
        dis = [f"{q} at {e:g}: τ₃ {fs(c_[0], c_[2])}, τ₄ {fs(c_[1], c_[3])}" for q, d in sqr.items()
               for e, c_ in sorted(d.items(), key=lambda t: -t[0]) if not c_[4]]
        pth = [f"{q}: " + '{' + ', '.join(map(str, gt[q]['block'])) + f"}} + {gt[q]['j3']}" + ('' if gt[q].get('j4') is None else f" + {gt[q]['j4']}")
               for q in sqr]
        mis = [f"{q} at {e:g}" for q, d in sqr.items() for e in d if prox.get(q, {}).get(e) != (not d[e][4])]
        out += ["The loop's sequence by kernel index (SQ, referee round 83: G's loop rerun at G's outputs, the formed sets at the loop's "
                "τ₃ and τ₄ against the path's block at τ_d with its entrants j3 and j4, `kappa_from_path`): the path's " + '; '.join(pth)
                + f". The loop follows it at {on} of {tot} cells (SQ's τ₃ and τ₄ are G's at {same} of {tot}); elsewhere: " + '; '.join(dis)
                + ". The intermediate slope's proxy (beyond 1 ± " + f"{SEQ_TOL:g} times the path's law, round 80) "
                + ("marks the same cells." if not mis else "differs at " + ', '.join(mis) + "."), '']
    return out


def gl_lines(GL):
    """GL's collapse (referee rounds 68, 74): per target the local exponents in τ of the dipole pair's distance, (W − A) and the residual,
    and the collapse law's constant over successive rows (`gl_constant`)."""
    out = []
    for g in GL:
        rs = [x for x in g['rows'] if x['d'] is not None and x['d'] > 0]
        if len(rs) < 3:
            continue
        ex = lambda k, a, b: float(np.log(abs(b[k]) / abs(a[k])) / np.log(b['tau'] / a['tau'])) if a[k] and b[k] else float('nan')
        pts = rs[::4] + ([rs[-1]] if (len(rs) - 1) % 4 else [])
        out += ['', f"### GL. The collapse (referee rounds 68, 74): target {g['q']} at ε = {EPS_GL:g} to τ = {TAUS_GL[-1]:g}", '',
                "The dipole pair's distance, (W − A)/|A| and the residual, and their local exponents in τ against the law's −1/6, −1/3, "
                "−2/3 (E ∝ Δc⁸ with the centre of mass pinned):", '',
                '| τ | distance | (W − A)/\\|A\\| | residual | exponents since the previous row: distance, W − A, residual |',
                '|---|---|---|---|---|']
        prev = None
        for x in pts:
            out.append(f"| {x['tau']:.4g} | {x['d']:.3e} | {x['dW']:+.3e} | {x['rms']:.3e} | "
                       + ('—' if prev is None else f"{ex('d', prev, x):+.4f}, {ex('dW', prev, x):+.4f}, {ex('rms', prev, x):+.4f}") + ' |')
            prev = x
        if 'Gm' in g and all('dc' in x for x in rs):
            cons = [gl_constant(a, b, g['Gm'], g['mu0'], g['Lam']) for a, b in zip(rs[:-1], rs[1:])]
            ef = [x['E_formed'] / x['E'] for x in rs]
            sel = list(range(0, len(cons), 4)) + ([len(cons) - 1] if (len(cons) - 1) % 4 else [])
            out += ['', "The law's constant (referee round 74) over successive rows: the measured d(Δc⁻⁶)/dτ over 48aGm(W²/(w_a² + "
                    "w_b²))Λ(ε), a = E/Δc⁸ (the proof's metric: the centres alone), and over it divided by 1 + C, C the slaved coefficients' "
                    "share of the separation's metric (`gl_constant`); p = |w_a − w_b||Δc|/(|W|σ̄) and C at the interval's end:", '',
                    '| τ | p | ratio | with C | C |', '|---|---|---|---|---|']
            out += [f"| {rs[k + 1]['tau']:.4g} | {cons[k][0]:.4f} | {cons[k][1]:.4f} | {cons[k][2]:.4f} | {cons[k][3]:.4f} |" for k in sel]
            out += ['', f"Over the {len(cons)} intervals the ratio is {min(c[1] for c in cons):.4f}–{max(c[1] for c in cons):.4f}, with C "
                    f"{min(c[2] for c in cons):.4f}–{max(c[2] for c in cons):.4f}; the formed kernels' energy over all kernels' "
                    f"{min(ef):.6f}–{max(ef):.6f}."]
            pr = np.array([abs(x['wa'] - x['wb']) * x['dc'] / (abs(x['wa'] + x['wb']) * x['sbar']) for x in rs])
            av = np.array([x['E'] / x['dc'] ** 8 for x in rs])
            if pr.max() - pr.min() > 0.05:                              # referee round 74: what the constant with C leaves
                sl = [np.log(b['E_formed'] / a['E_formed']) / np.log(b['dc'] / a['dc']) for a, b in zip(rs[:-1], rs[1:])]
                st = [np.log(b['E'] / a['E']) / np.log(b['dc'] / a['dc']) for a, b in zip(rs[:-1], rs[1:])]
                ks = [(x - 8) / (0.5 * (pa + pb)) for x, pa, pb in zip(sl, pr[:-1], pr[1:])]   # referee round 80: the local slope
                kt = [(x - 8) / (0.5 * (pa + pb)) for x, pa, pb in zip(st, pr[:-1], pr[1:])]
                rem = [c[2] / (x / 8) - 1 for c, x in zip(cons, st)]           # referee round 83: a = E/Δc⁸'s E, all kernels
                out += [f"The remainder with C against p: (ratio with C − 1)/p {min((c[2] - 1) / c[0] for c in cons):.3f}–"
                        f"{max((c[2] - 1) / c[0] for c in cons):.3f}; a = E/Δc⁸ itself goes {av[0]:.4g} → {av[-1]:.4g} along the run. "
                        "E's next order in Δc from the local slope (referee round 80; the formed kernels' energy, the idle ones' floor "
                        f"left out): (d ln E/d ln Δc − 8)/p over the intervals (p their mean) {min(ks):.3f}–{max(ks):.3f}, median "
                        f"{np.median(ks):.3f} (all kernels' energy: {min(kt):.3f}–{max(kt):.3f}, the idle kernels' floor in it), which "
                        f"raises the force over 8aΔc⁷ by 1 + {np.median(ks) / 8:.3f}p; with each interval's own force factor (d ln E/d ln "
                        f"Δc)/8 — E all kernels' energy, as in a = E/Δc⁸ (referee round 83) — the ratio with C is 1 {min(rem):+.1e} … "
                        f"{max(rem):+.1e} (the largest |·| at "
                        f"τ = {rs[1 + int(np.argmax(np.abs(rem)))]['tau']:.4g}); the mean width's share of the metric, p²/36, at most "
                        f"{pr.max() ** 2 / 36:.5f}."]
        out.append('')
    return out


def h_table(H, sep=False):
    """H's table (H and H3)."""
    out = [('| separation | ' if sep else '| target | ') + 'outcome | visible negative kernels along the path (most negative w) | at the end '
           '(all formed) | smallest visible width | sub-sample, negative widths (all active; formed) | formed invisible | max\\|w\\| | residual | '
           'the main kernels\' errors (w/A − 1) | the others: n, Σw (Σ\\|w\\|), on the samples, widths, \\|c\\| max |',
           '|---|---|---|---|---|---|---|---|---|---|---|']
    for r in H:
        lab = f"±{r['cs'][1]:g} ({(r['cs'][1] - r['cs'][0]) / S0:.2f}σ₀)" if sep else r['name']
        out.append(f"| {lab} | {r['outcome']} | {r['neg_path']} ({r['wmin_path']:+.3f}) | {r['neg_end']} ({r.get('neg_end_all', '—')}) | "
                   f"{r['smin']:.2e} | {r.get('sub', '—')}, {r.get('negw', '—')} ({r.get('sub_formed', '—')}, {r.get('negw_formed', '—')}) | "
                   f"{r.get('inv_formed', '—')} | {r['wmax']:.3g} | {r['moving']['rms']:.1e} | {split_cells(end_split(r))} |")
    return out + ['']


def g_loop_lines(G):
    """G's loop tables (referee round 61): the states at TAUS_G_MD — the dipole pair (`dipole_pair`) with its distance, coefficients,
    centre of mass's offset from the nearest feature (in its width) and net coefficient against the feature's amplitude, the coinciding
    pairs, the largest |w| over the largest amplitude, the residual —; the three phases (`phases`); the landed v against ε."""
    if any('theta' not in x for r in G for x in r['loops']):
        return ['', "(G's loops of a run before referee round 61, without the loop's states: rerun G.)", '']
    L_ = ['', "The real loop at the masses ε from α = 0 (`loop_nd`, rtol 10⁻⁹) at τ = " + ', '.join(f"{t:g}" for t in TAUS_G_MD) + ": the "
          "formed kernels, the dipole pair (the largest |w| and the closest kernel of the other sign: its distance |Δc|/σ + |Δlog σ|, "
          "coefficients, its centre of mass's offset from the nearest feature in that feature's width and its net coefficient W against the "
          f"feature's amplitude, (W − A)/|A|; landed (L) if both within {LAND_TOL:g}), the coinciding pairs (within {PAIR_TOL:g}), the "
          "largest |w| over the largest amplitude, the residual:", '',
          '| Δ ≥ | target | ε | τ | formed | the dipole pair: distance, coefficients, offset (σ*), (W − A)/\\|A\\| | coinciding pairs | '
          'max\\|w\\|/max\\|A\\| | residual |', '|---|---|---|---|---|---|---|---|---|']
    for r in G:
        for x in sorted(r['loops'], key=lambda x: (x['eps'], x['tau'])):
            if x['tau'] not in TAUS_G_MD:
                continue
            lp = landing(x['theta'], x['w'], r['cs'], r['ss'], r['As'])
            th = np.atleast_2d(np.asarray(x['theta'], float)).reshape(-1, 2)
            L_.append(f"| {r['dmin']:g} | {r['q']} | {x['eps']:g} | {x['tau']:g} | {len(x['w'])} | "
                      + ('—' if lp is None else f"{lp['d']:.2e}, {lp['w'][0]:+.3f}/{lp['w'][1]:+.3f}, {lp['off']:.1e}, {lp['dW']:+.1e}"
                         + (' L' if lp['off'] <= LAND_TOL and abs(lp['dW']) <= LAND_TOL else '')) + ' | '
                      + (', '.join(f"{a:+.3f}/{b:+.3f}" for a, b in pairs_of(th, x['w'])) or 'none')
                      + f" | {x['wmax'] / np.abs(r['As']).max():.2f} | {x['rms']:.1e} |")
    spl = []                                                           # (ii)(a)'s width split on the landed pairs (referee round 61)
    for r in G:
        for x in r['loops']:
            lp = landing(x['theta'], x['w'], r['cs'], r['ss'], r['As']) if x['tau'] in TAUS_G_MD else None
            if lp is None or lp['off'] > LAND_TOL or abs(lp['dW']) > LAND_TOL:
                continue
            th, w = np.atleast_2d(np.asarray(x['theta'], float)).reshape(-1, 2), np.asarray(x['w'], float)
            a, b = dipole_pair(th, w)[:2]
            dc, ds, pr = width_split(th, w, a, b)
            spl.append((ds / pr, abs(w[a] - w[b]) * abs(dc) / (abs(w[a] + w[b]) * 0.5 * (abs(th[a, 1]) + abs(th[b, 1])))))
    if spl:
        rt, pm = np.array([s_[0] for s_ in spl]), np.array([s_[1] for s_ in spl])
        L_ += ['', f"(ii)(a)'s width split on the landed pairs (L, {len(spl)} rows, {len({round(float(x), 6) for x in rt})} distinct to six "
               "decimals: frozen states repeat): the observed σ_b − σ_a over (w_b − w_a)Δc²/(6Wσ) "
               f"{rt.min():.4f}–{rt.max():.4f} (median {np.median(rt):.4f}), at |w_a − w_b||Δc|/(|W|σ) = {pm.min():.2e}–{pm.max():.2e}; "
               f"its largest deviation {np.abs(rt - 1).max():.1e} at the parameter {pm[np.argmax(np.abs(rt - 1))]:.2e}."]
    L_ += ['', "The three phases (`phases`): the dipole block's entry τ_d (the path's) and the next kernel's formation τ₃ in the loop; over "
           "the grid times between them the slopes of v³ (v the dipole pair's half-difference) against 3μ₀KΛ(ε) — the block's law in its "
           "fast time Λ(ε)τ, K and μ₀ from the flow's last output — (median [range]); the landed v (at τ = 400) against the law's "
           "(3μ₀KΛ(ε)(τ₃ − τ_d))^{1/3}; landed at τ = " + ', '.join(f"{t:g}" for t in TAUS_G_MD) + "; over 1600 → 6400 the exponents in τ of "
           "the pair's distance and of the residual:", '',
           '| Δ ≥ | target | ε | τ_d | τ₃ | d(v³)/dτ / 3μ₀KΛ | v at 400 | the law\'s v | landed | exponents: distance, residual |',
           '|---|---|---|---|---|---|---|---|---|---|']
    law_v, phs = {}, {}
    for r in G:
        for eps in EPS_G:
            ph = phases(r, eps)
            if ph is None:
                continue
            sl = ph['slopes']
            law_v[(r['dmin'], r['q'], eps)] = ph['v400']
            phs[(r['dmin'], r['q'], eps)] = ph
            L_.append(f"| {r['dmin']:g} | {r['q']} | {eps:g} | {r['tau_d']:.4g} | " + ('—' if ph['t3'] is None else f"{ph['t3']:g}") + ' | '
                      + (f"{np.median(sl):.3f} [{min(sl):.3f}, {max(sl):.3f}]" if sl else '—') + ' | '
                      + ('—' if ph['v400'] is None else f"{ph['v400']:.4g}") + ' | ' + ('—' if ph['vstar'] is None else f"{ph['vstar']:.4g}")
                      + ' | ' + ', '.join('—' if ph['land'][t] is None else ('yes' if ph['land'][t]['off'] <= LAND_TOL
                                                                                and abs(ph['land'][t]['dW']) <= LAND_TOL else 'no')
                                          for t in TAUS_G_MD) + ' | '
                      + ('—' if ph['collapse'] is None else f"{ph['collapse'][0]:+.3f}, {ph['collapse'][1]:+.3f}") + ' |')
    L_ += ['', "The landed v against the mass: v at τ = 400 and its ratio per decade of ε (10^{1/6} = " + f"{10 ** (1 / 6):.4f} if v ∝ "
           "ε^{−1/6}):", '', '| Δ ≥ | target | v at ε = ' + ', '.join(f"{e:g}" for e in EPS_G) + ' | ratios per decade |', '|---|---|---|---|']
    for r in G:
        vs = [law_v.get((r['dmin'], r['q'], e)) for e in EPS_G]
        rat = [f"{b / a:.4f}" if a and b else '—' for a, b in zip(vs[:-1], vs[1:])]
        L_.append(f"| {r['dmin']:g} | {r['q']} | " + ', '.join('—' if v is None else f"{v:.4g}" for v in vs) + f" | {', '.join(rat)} |")
    L_ += g_phase_summary(phs, law_v)
    return L_ + ['']


def g_phase_summary(phs, law_v):
    """G's three phases summed up (referee round 61): per ε the cells landed at τ = 6400; at the smallest ε the dipole phase's slope over
    the law (the cells' medians) and the landed v over the law's; the last decade's ratio of v; on the cells landed at both 1600 and 6400
    the collapse's exponents (distance, residual) and those not yet collapsing by τ = 6400 (|exponent| below 10⁻³: landed at a small p,
    they follow the same law later — referee round 68, GL)."""
    if not phs:
        return []
    ok = lambda lp: lp is not None and lp['off'] <= LAND_TOL and abs(lp['dW']) <= LAND_TOL
    e0 = min(EPS_G)
    n6 = {e: (sum(ok(ph['land'].get(6400.0)) for k, ph in phs.items() if k[2] == e), sum(k[2] == e for k in phs)) for e in EPS_G}
    sl = [np.median(ph['slopes']) for k, ph in phs.items() if k[2] == e0 and ph['slopes']]
    vr = [ph['v400'] / ph['vstar'] for k, ph in phs.items() if k[2] == e0 and ph['v400'] and ph['vstar']]
    dec = [law_v[k[:2] + (e0,)] / law_v[k[:2] + (10 * e0,)] for k in phs if k[2] == e0 and law_v.get(k[:2] + (10 * e0,)) and law_v.get(k)]
    both = [(k, ph['collapse']) for k, ph in phs.items() if ok(ph['land'].get(1600.0)) and ok(ph['land'].get(6400.0)) and ph['collapse']]
    fr = [k for k, c in both if abs(c[0]) < 1e-3]
    col = np.array([c for k, c in both if abs(c[0]) >= 1e-3])
    span = lambda v, f_: f"{min(v):{f_}}–{max(v):{f_}}" if len(v) else '—'
    return ['', "Over the six targets: landed at τ = 6400 in " + ', '.join(f"{a} of {b} at ε = {e:g}" for e, (a, b) in n6.items())
            + f"; at ε = {e0:g} the dipole phase's slope over the law {span(sl, '.3f')} (the cells' medians, {len(sl)} cells with a "
            f"phase), the landed v over the law's at the grid's τ₃ {span(vr, '.2f')} (on the fine grid κ, below), v's ratio over the last "
            f"decade {span(dec, '.4f')}; on the {len(both)} cells landed at both 1600 and 6400, {len(fr)} not yet collapsing by τ = 6400 "
            "(|exponent| < 10⁻³: landed at a small p, the same law later — GL) (" + ', '.join(f"{k[1]} at {k[2]:g}" for k in fr)
            + ")" + (f", the others collapsing with the exponents {col[:, 0].min():+.3f} … {col[:, 0].max():+.3f} (distance) and "
                     f"{col[:, 1].min():+.3f} … {col[:, 1].max():+.3f} "
                     "(residual)" if len(col) else '') + '.']


KEYS = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'H2', 'H3', 'HL', 'HP', 'GL', 'GT', 'A2', 'HLL', 'SQ', 'A2E', 'A2P', 'A2H', 'A2R', 'A2T')   # report's order


def main():
    warnings.filterwarnings('ignore')
    raw = os.path.join(OUT, 'moving_width_raw.npz')
    md = os.path.join(OUT, 'moving_width.md')
    old = E.load(raw)
    keys = KEYS
    if '--from-raw' in sys.argv:
        open(md, 'w').write(report(*[old.get(k) for k in keys]))
        return
    parts = set(sys.argv[sys.argv.index('--only') + 1].split(',')) if '--only' in sys.argv else set(keys)
    res = dict(old)
    t0 = time.time()
    for k, fn in (('A', part_A), ('B', part_B), ('C', part_C), ('D', part_D), ('E', part_E), ('F', part_F),
                  ('G', lambda: part_G(res['C'])), ('H', part_H), ('H2', part_H2), ('H3', part_H3), ('HL', part_HL), ('HP', part_HP),
                  ('GL', lambda: part_GL(res['G'])), ('GT', lambda: part_GT(res['G'])), ('A2', part_A2), ('HLL', part_HLL),
                  ('SQ', lambda: part_SQ(res['G'])), ('A2E', part_A2E), ('A2P', part_A2P), ('A2H', lambda: part_A2H(res['A2P'])),
                  ('A2R', lambda: part_A2R(res['A2P'])), ('A2T', lambda: part_A2T(res['A2P']))):
        if k in parts:
            res[k] = fn(); E.save(raw, **res); print(f"[{time.time() - t0:.0f}s] {k}", flush=True)
            open(md, 'w').write(report(*[res.get(k) for k in keys]))


if __name__ == '__main__':
    main()
