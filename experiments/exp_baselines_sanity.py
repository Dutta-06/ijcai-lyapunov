"""
Are the order-selection baselines of the 432-trial comparison a strawman?  A post hoc check of forward selection, penalised training and random start + prune
=================================================================================================================================================================

Why.  In `results/chaos/noise_order.md` (the 432-trial comparison of `experiments/exp_noise_order.py`, the source of the supplement's order-selection table) the three
grid-based or network baselines "forward selection on the grid, BIC" (`step`), "penalised network, BIC" (`pen`) and "random start + prune" (`prune`) reach 0-10 % exact order
(`pen` and `prune` 0-1 % at every noise level, even at s = 0.01).  A reader can ask whether (a) an implementation bug, (b) the design (centres fixed on a grid cannot fit an
off-grid feature with one unit, and counting units that coincide or neighbour each other is not what the method estimates), or (c) badly chosen hyper-parameters
(penalty, threshold, pruning rule) is behind it.  Everything here is POST HOC: the frozen `src/core/baselines.py` and `experiments/exp_noise_order.py` are imported and not
edited; the corrected and fair variants live in this file.

The setting is that of `exp_noise_order.py` (x in [-1, 1], 61 samples, a dictionary of 21 Gaussian kernels of width 0.12 at spacing 0.1, K-hat = K scored as "exact").

THE BASELINES AS THE COMPARISON RUNS THEM (frozen code, reproduced here with the same seeds)
  step    `baselines.stepwise_bic`: forward selection by least squares among the 21 grid kernels (centres stay on the grid), at most 10 units, BIC with ONE parameter per unit,
          K-hat = the number of grid units chosen.
  pen     `baselines.penalised_net`: 21 units, centres start on the grid and are trained, coefficients start at 0.1 N(0,1) (seed = target seed), plain gradient descent (step
          1/L, L the largest eigenvalue of the Gram matrix at the start; centres 0.1/L) with a proximal l1 step on the coefficients, 2500 steps, for 7 penalties
          lam = geomspace(3e-4, 3e-1, 7); the penalty by BIC with TWO parameters per nonzero coefficient; K-hat = the NUMBER OF NONZERO COEFFICIENTS (|w| > 1e-6 max|w|), units
          that sit at the same centre are counted one by one.
  prune   `baselines.prune_net`: the same 21-unit network, coefficients start at 0.3 N(0,1), centres trainable, no penalty, 5000 steps, units with |w| < 0.1 max|w| removed;
          K-hat = the number of survivors.
  (The comparison's blasso rows `sfw` and `cpgd` post-process their solution before counting: `threshold_units` drops weak spikes and `cpgd` merges particles closer than
  1.5 grid steps into one spike.  `pen` and `prune` get neither.)

DECLARED ANALYSES (declared before the runs; nothing was added after seeing results except where marked "added after")
  A0  Implementation checks.  (i) Reproduction: the three frozen baselines are re-run on the 432 trials of the main design (same seeds) and their K-hat compared with the stored
      raw `results/chaos/noise_order_raw.npz` (columns step|bic, pen|bic, prune).  (ii) Gradient: the update of `train_net` against finite differences of the loss, in w and in
      c.  (iii) Fixed point: the planted solution is a fixed point of `train_net` without penalty on noiseless data.  (iv) Convergence: `pen` (lam = 3e-3) with 2500 steps
      against 20000 steps (objective and unit count), `prune` with 5000 against 50000 steps (residual and unit count), on three planted cases.  (v) Optimality: on the sanity
      trials, the objective J = 1/2 mean r^2 + lam |w|_1 reached by `pen` at each of its 7 penalties against J of the planted truth (a network that fails to optimise would
      have J above it).
  A1  Sanity cases with a known answer ("sanity"): K = 1, and K = 2 with the centres at least 10 widths (1.2) apart, the largest amplitude 1, the other in [0.5, 1], random signs;
      placement ON the dictionary grid, OFF-grid at the worst case (every centre half a grid step, 0.05, from the grid) and OFF-grid at random (uniform in [-0.85, 0.85]);
      noise s = 0 (noiseless; BIC with the residual mean square floored at (1e-4)^2, declared because log of a vanishing residual is not defined), 0.001, 0.01 of the largest
      amplitude; 20 targets per (K, placement) cell, 1 draw if noiseless and 2 draws otherwise: 600 trials.  For each baseline: exact order (K-hat = K), exact and placed (every
      unit matched within half a width), over- and under-estimation, the failure mode.  Positive controls on the same trials: matching pursuit + sliding with BIC and
      sliding Frank-Wolfe with BIC (the comparison's own methods whose centres are free).
  A2  Three counting conventions for each baseline (the design question (b)): V0 = the comparison's own (units counted one by one); V1 = FAIR COUNT: the units are merged into
      clusters (single linkage on the centres, link 0.12 = one width, which on the 0.1 grid joins neighbouring grid kernels; the coincident units of a trained network have
      distance ~0), a cluster has the sum of the coefficients as its mass and the |w|-weighted mean as its position, clusters with |mass| below 0.1 of the largest are dropped
      (the post-processing of `threshold_units`), K-hat = the number of clusters; V2 = V1 followed by a joint slide of the clusters' centres and masses on the least-squares
      objective (`baselines.refit_block`) and a second merge of units that came to coincide (centre refinement after selection).  The sensitivity of V1 to the link
      (0.06, 0.12, 0.15, 0.2, 0.3) is reported on the sanity trials.
  A3  Hyper-parameter retuning (the design question (c)) on DEVELOPMENT targets disjoint from every target of the main comparison: seeds 1000 K + 100 + t, t = 0..7, K = 1..6,
      the four noise levels, draw 0 (192 trials).  Each baseline's own hyper-parameters are chosen by the largest pooled exact-order rate on the development trials (ties: the
      first of the listed grid, which starts at the comparison's own setting).  Grids -- step: the BIC multiplier c in (1, 2, 3, 4, 6, 8, 12, 16) [BIC = N ln(rss/N) + c k ln N]
      (V1 also: link in (0.12, 0.06, 0.15, 0.2, 0.3), mass threshold in (0.1, 0, 0.05, 0.2)); pen: the rule in (BIC with c in (1, 2, 4, 8) [2 c nnz ln N], or one fixed penalty out
      of the 7), the count threshold theta in (1e-6, 1e-3, 1e-2, 0.05, 0.1, 0.2) [unit if |w| > theta max|w|] (V1: instead link and mass threshold as above, theta = 1e-6);
      prune: the initial scale in (0.3, 0.1), the pruning threshold thr in (0.1, ..., 0.9) (V1: also the link).  The chosen settings are then scored on the MAIN design's 432
      trials (same targets, noise draws and design as the comparison), which no tuning ever saw, and on the sanity trials.  A fresh scale of the initial network (0.03, the
      "small start" of the paper's own method) was NOT tried: it would change the baseline into the method it is compared with.
  A4  The main comparison with the corrected counting: exact order, by noise level and by K, of the three baselines with V0 (the comparison's), V1, V2 and the tuned
      settings, beside the stored rows of the sliding ISS (BIC and discrepancy), matching pursuit + sliding (BIC), sliding Frank-Wolfe (BIC), conic particle gradient descent
      (BIC), forward selection, the merged l1 row, `pen` and `prune`.  Intervals: 95 % Wilson over the trials AND 95 % percentile bootstrap over the 36 targets (2000 resamples,
      seed 0; the trials of a target are resampled together, as in `exp_noise_order.py`); differences from the sliding ISS (BIC) by the original's
      `target_cluster_diffs`.

ADDED AFTER the first look at the sanity trials (not declared in advance, marked as such in the report): the sensitivity of V1 to the link on the main design, and the summary
section (which only re-reads numbers printed below it).  The link 0.12, the mass threshold 0.1 and every grid above were fixed before any run.

Scoring is the comparison's: K-hat = K (exact order), among those the units matched within half a width (`baselines.match_units`, tol 0.06).  "Exact and placed" additionally
needs every one of the K true features matched.

Usage:
    python experiments/exp_baselines_sanity.py --run [--stage checks,sanity,dev,main] [--workers 1]   (writes results/synthetic/baselines_sanity_raw.npz; BS_OUT overrides the directory)
    python experiments/exp_baselines_sanity.py --from-raw                                              (regenerates results/synthetic/baselines_sanity.md, byte-identically)
"""
import argparse
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)
import exp_noise_order as X          # noqa: E402  (frozen: the main comparison)

B = X.B
OUT = os.environ.get('BS_OUT', os.path.join(REPO, 'results', 'synthetic'))
RAW, MD = os.path.join(OUT, 'baselines_sanity_raw.npz'), os.path.join(OUT, 'baselines_sanity.md')
REF_RAW = os.path.join(REPO, 'results', 'chaos', 'noise_order_raw.npz')

SIG, TOL = X.SIG, X.TOL
GRID = np.linspace(-1, 1, X.ND)
ND = X.ND
PEN_LAMS = X.PEN_LAMS
NPEN = len(PEN_LAMS)
STEPS_PEN, STEPS_PRUNE = 2500, 5000
SCALES = (0.3, 0.1)
MAXU = X.MAX_UNITS
FLOOR = 1e-4                                           # residual rms floor of the information criteria (matters only for noiseless data)
SAN_NOISE = (0.0, 0.001, 0.01)
SAN_T = 20
SAN_SEP = 10 * SIG                                     # K = 2: the centres at least ten widths apart
PLC = ('on grid', 'off grid, half a step', 'off grid, random')
DEV_T = 8
LINKS = (0.12, 0.06, 0.15, 0.2, 0.3)                   # the first is the default of V1
THMS = (0.1, 0.0, 0.05, 0.2)
C_STEP = (1, 2, 3, 4, 6, 8, 12, 16)
PEN_RULES = tuple(('bic', c) for c in (1, 2, 4, 8)) + tuple(('lam', j) for j in range(NPEN))
THETAS = (1e-6, 1e-3, 1e-2, 0.05, 0.1, 0.2)
THRS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
METHODS = ('step', 'pen', 'prune')
MLABEL = {'step': 'forward selection on the grid, BIC', 'pen': 'penalised network, BIC', 'prune': 'random start + prune'}
CTRL = ('omp|bic', 'sfw|bic')
CTRL_LABEL = {'omp|bic': 'matching pursuit + sliding, BIC (positive control)', 'sfw|bic': 'sliding Frank-Wolfe, BIC (positive control)'}
REF_ROWS = (('iss|bic', 'sliding ISS, BIC'), ('iss|disc', 'sliding ISS, discrepancy (s)'), ('omp|bic', 'matching pursuit + sliding, BIC'), ('sfw|bic', 'sliding Frank-Wolfe, BIC'),
            ('cpgd|bic', 'conic particle gradient descent, BIC'), ('step|bic', 'forward selection on the grid, BIC (stored)'), ('lasso|merged', 'l1 on the grid, BIC (merged atoms; stored)'),
            ('pen|bic', 'penalised network, BIC (stored)'), ('prune', 'random start + prune (stored)'))


# ---------------------------------------------------------------------------------------------------------------------------
# the fair count: merging, the information criterion, the slide
# ---------------------------------------------------------------------------------------------------------------------------
def merge_units(w, c, link, thm):
    """Clusters of units: single linkage on the centres (consecutive sorted centres farther apart than `link` start a new cluster), mass = the sum of the coefficients, position = the
    |w|-weighted mean of the centres; clusters of |mass| below thm times the largest (or of zero mass) are dropped. Units of zero coefficient are ignored. Returns (masses, positions)."""
    w, c = np.asarray(w, float).ravel(), np.asarray(c, float).ravel()
    m = w != 0
    w, c = w[m], c[m]
    if len(w) == 0:
        return np.zeros(0), np.zeros(0)
    o = np.argsort(c, kind='stable')
    w, c = w[o], c[o]
    lab = np.concatenate([[0], np.cumsum(np.diff(c) > link)]).astype(int)
    mass = np.bincount(lab, weights=w)
    absw = np.bincount(lab, weights=np.abs(w))
    pos = np.bincount(lab, weights=np.abs(w) * c) / absw
    keep = (np.abs(mass) >= thm * np.abs(mass).max()) & (np.abs(mass) > 0)
    return mass[keep], pos[keep]


def ic_value(rss, N, n_params):
    """N ln(max(rss, N FLOOR^2)/N) + n_params ln N: the BIC with the residual mean square floored (noiseless data)."""
    return N * np.log(max(rss, N * FLOOR ** 2) / N) + n_params * np.log(N)


def ic_select(rss, per_unit, N, cmul=1.0):
    """argmin over k of the floored BIC with cmul * per_unit * k parameters; rss[k] the residual sum of squares of the block of k units."""
    vals = [ic_value(r, N, cmul * per_unit * k) for k, r in enumerate(rss)]
    return int(np.argmin(vals))


def units_step(idx, w_pad, rss, cfg, merged, N):
    """Forward selection: the block chosen by the BIC multiplier cfg['c'], then (merged) the clusters. Returns (K-hat, masses, positions)."""
    k = ic_select(rss, 1, N, cfg['c'])
    w = np.asarray(w_pad[k - 1, :k], float) if k > 0 else np.zeros(0)
    c = GRID[np.asarray(idx[:k], int)]
    if merged:
        W, C = merge_units(w, c, cfg['link'], cfg['thm'])
        return len(W), W, C
    return k, w, c


def _nz(w, theta):
    return np.abs(w) > theta * max(np.abs(w).max(), 1e-12)


def units_pen(w, c, rss, cfg, merged, N):
    """Penalised network: the penalty by cfg['rule'] (('bic', c) over the ladder, or ('lam', j) one fixed penalty), the units counted one by one above cfg['theta'] (V0) or merged (V1)."""
    th = cfg['theta'] if not merged else THETAS[0]
    rule = cfg['rule']
    if rule[0] == 'lam':
        j = rule[1]
    else:
        j = int(np.argmin([ic_value(rss[q], N, 2 * rule[1] * int(_nz(w[q], th).sum())) for q in range(len(rss))]))
    if merged:
        W, C = merge_units(w[j][_nz(w[j], th)], c[j][_nz(w[j], th)], cfg['link'], cfg['thm'])
        return len(W), W, C
    nz = _nz(w[j], th)
    return int(nz.sum()), w[j][nz], c[j][nz]


def units_prune(w, c, cfg, merged):
    """Random start + prune: the net of initial scale SCALES[cfg['scale']]; units of |w| >= thr max|w| kept (V0), or the units merged first and the clusters of |mass| >= thr max kept (V1)."""
    ws, cs = w[cfg['scale']], c[cfg['scale']]
    if merged:
        W, C = merge_units(ws, cs, cfg['link'], cfg['thr'])
        return len(W), W, C
    kept = np.abs(ws) >= cfg['thr'] * np.abs(ws).max()
    return int(kept.sum()), ws[kept], cs[kept]


DEFAULT = {'step': ({'c': 1}, {'c': 1, 'link': LINKS[0], 'thm': THMS[0]}),
           'pen': ({'rule': ('bic', 1), 'theta': THETAS[0]}, {'rule': ('bic', 1), 'link': LINKS[0], 'thm': THMS[0]}),
           'prune': ({'scale': 0, 'thr': 0.1}, {'scale': 0, 'link': LINKS[0], 'thr': 0.1})}


def config_grid(method, merged):
    """The development grid of one baseline and counting convention, the comparison's own setting first."""
    if method == 'step':
        L = [{'c': c} for c in C_STEP] if not merged else [{'c': c, 'link': l, 'thm': t} for c in C_STEP for l in LINKS for t in THMS]
    elif method == 'pen':
        L = [{'rule': r, 'theta': t} for r in PEN_RULES for t in THETAS] if not merged else [{'rule': r, 'link': l, 'thm': t} for r in PEN_RULES for l in LINKS for t in THMS]
    else:
        L = [{'scale': s, 'thr': t} for s in range(len(SCALES)) for t in THRS] if not merged else [{'scale': s, 'link': l, 'thr': t} for s in range(len(SCALES)) for l in LINKS for t in THRS]
    d = DEFAULT[method][1 if merged else 0]
    L.sort(key=lambda g: g != d)                       # stable: the default first, the rest in listed order
    return L


def cfg_str(method, cfg):
    if method == 'step':
        s = 'BIC x%g' % cfg['c']
    elif method == 'pen':
        s = ('BIC x%g' % cfg['rule'][1] if cfg['rule'][0] == 'bic' else 'fixed lam = %.3g' % PEN_LAMS[cfg['rule'][1]]) + (', theta = %g' % cfg['theta'] if 'theta' in cfg else '')
    else:
        s = 'scale %g, thr %g' % (SCALES[cfg['scale']], cfg['thr'])
    if 'link' in cfg:
        s += ', link %g' % cfg['link'] + (', mass thr %g' % cfg['thm'] if 'thm' in cfg else '')
    return s


def slide_units(Kd, y, W, C):
    """V2: the clusters' masses and centres refitted jointly by least squares (`baselines.refit_block`), then units that came to coincide merged again. Returns (K-hat, masses, positions)."""
    if len(W) == 0:
        return 0, W, C
    w2, c2 = B.refit_block(Kd, y, np.asarray(W, float), np.asarray(C, float).reshape(-1, 1))
    W2, C2 = merge_units(w2, c2.ravel(), LINKS[0], THMS[0])
    return len(W2), W2, C2


# ---------------------------------------------------------------------------------------------------------------------------
# the designs
# ---------------------------------------------------------------------------------------------------------------------------
def sanity_target(K, plc, seed):
    """A sanity target: K = 1 or 2 (centres at least SAN_SEP = 10 widths apart), placement plc = 0 on the grid, 1 half a grid step off it, 2 uniform at random; the largest amplitude 1."""
    rng = np.random.default_rng(seed)
    step = GRID[1] - GRID[0]
    if K == 1:
        j = int(rng.integers(-8, 9))
        c = np.array([0.1 * j + (0.5 * step if plc == 1 else 0.0)]) if plc < 2 else np.array([rng.uniform(-0.85, 0.85)])
        a = np.array([rng.choice([-1.0, 1.0])])
    else:
        if plc < 2:
            j1 = int(rng.integers(-8, -3))
            j2 = int(rng.integers(j1 + 12, 9))
            c = 0.1 * np.array([j1, j2]) + (0.5 * step if plc == 1 else 0.0)
        else:
            c1 = rng.uniform(-0.85, 0.85 - SAN_SEP)
            c = np.array([c1, rng.uniform(c1 + SAN_SEP, 0.85)])
        a = np.array([1.0, rng.uniform(0.5, 1.0)]) * rng.choice([-1.0, 1.0], 2)
        if rng.random() < 0.5:
            a = a[::-1].copy()
    return np.round(c, 10), a


def sanity_jobs():
    """(stage 0, K, noise index, target seed, draw, placement, noise value) for the 600 sanity trials."""
    jobs = []
    for K in (1, 2):
        for plc in range(3):
            cell = 3 * (K - 1) + plc
            for t in range(SAN_T):
                for ni, s in enumerate(SAN_NOISE):
                    for d in range(1 if s == 0 else 2):
                        jobs.append((0, K, ni, 500000 + 1000 * cell + t, d, plc, s))
    return jobs


def dev_jobs():
    return [(1, K, ni, 1000 * K + 100 + t, 0, -1, X.NOISE[ni]) for K in X.KS for t in range(DEV_T) for ni in range(len(X.NOISE))]


def main_jobs(n_targets=6, n_draws=3):
    return [(2, K, ni, 1000 * K + t, d, -1, X.NOISE[ni]) for K in X.KS for t in range(n_targets) for ni in range(len(X.NOISE)) for d in range(n_draws)]


STAGES = {'sanity': (0, sanity_jobs), 'dev': (1, dev_jobs), 'main': (2, main_jobs)}


def _pad(a, n):
    out = np.full(n, np.nan)
    out[:len(a)] = a
    return out


def _trial(args):
    """One trial (stage code, K, noise index, target seed, draw, placement, noise value) -> dict of arrays: the planted truth, the raw output of the three baselines (so that every
    counting convention and every hyper-parameter setting can be scored afterwards), the V2 scores at the default settings, the controls (sanity only) and the seconds."""
    stg, K_true, ni, seed, draw, plc, nz = args
    Kd, grid = X.setup()
    N = Kd.N
    ct, at = sanity_target(K_true, plc, seed) if stg == 0 else X.make_target(np.random.default_rng(seed), K_true)
    f = Kd.synth(ct[:, None], at)
    s = nz * np.abs(at).max()
    y = f + s * np.random.default_rng(seed * 1000 + draw + 1).standard_normal(N) if s > 0 else f.copy()
    secs = np.zeros(5)
    t0 = time.time()
    S = B.stepwise_bic(Kd, y, grid, max_units=MAXU)
    idx = np.array(S['order'], int)
    srss = np.array([b['rss'] for b in S['blocks']])
    sw = np.full((MAXU, MAXU), np.nan)
    for k in range(1, MAXU + 1):
        sw[k - 1, :k] = S['blocks'][k]['w']
    secs[0] = time.time() - t0
    t0 = time.time()
    net = B.penalised_net(Kd, y, grid, PEN_LAMS, steps=STEPS_PEN, seed=seed)
    pw, pc, prss = np.stack([q['w'] for q in net]), np.stack([q['c'].ravel() for q in net]), np.array([q['rss'] for q in net])
    secs[1] = time.time() - t0
    t0 = time.time()
    qw, qc = np.zeros((len(SCALES), ND)), np.zeros((len(SCALES), ND))
    for i, sc in enumerate(SCALES):
        Z = B.prune_net(Kd, y, grid, steps=STEPS_PRUNE, seed=seed, scale=sc)
        qw[i], qc[i] = Z['w'], Z['c'].ravel()
    secs[2] = time.time() - t0
    ctrl = np.full((2, 4), np.nan)
    if stg == 0:
        t0 = time.time()
        O = B.omp_slide(Kd, y, grid, max_units=MAXU, s=None)
        k = ic_select([b['rss'] for b in O['blocks']], 2, N)
        ctrl[0] = X._score(len(O['blocks'][k]['w']), O['blocks'][k]['w'], O['blocks'][k]['c'], ct, at)
        G = float(np.abs(Kd.corr(grid[:, None], y)[0]).max())
        SF = [B.threshold_units(Kd, y, B.sfw(Kd, y, grid, fr * G, max_units=MAXU), fr * G) for fr in X.LAM_REL]
        k = int(np.argmin([ic_value(q['rss'], N, 2 * len(q['w'])) for q in SF]))
        ctrl[1] = X._score(len(SF[k]['w']), SF[k]['w'], SF[k]['c'], ct, at)
        secs[3] = time.time() - t0
    t0 = time.time()
    v2 = np.full((3, 4), np.nan)
    for i, (m, un) in enumerate((('step', units_step(idx, sw, srss, DEFAULT['step'][1], True, N)), ('pen', units_pen(pw, pc, prss, DEFAULT['pen'][1], True, N)),
                                 ('prune', units_prune(qw, qc, DEFAULT['prune'][1], True)))):
        kh, W, C = slide_units(Kd, y, un[1], un[2])
        v2[i] = X._score(kh, W, C[:, None], ct, at)
    secs[4] = time.time() - t0
    return dict(ct=_pad(ct, 6), at=_pad(at, 6), meta=np.array([K_true, nz, seed, draw, plc, ni, s]), step_idx=idx, step_w=sw, step_rss=srss, pen_w=pw, pen_c=pc, pen_rss=prss,
                pr_w=qw, pr_c=qc, ctrl=ctrl, v2=v2, secs=secs)


def _shapes():
    return dict(ct=(6,), at=(6,), meta=(7,), step_idx=(MAXU,), step_w=(MAXU, MAXU), step_rss=(MAXU + 1,), pen_w=(NPEN, ND), pen_c=(NPEN, ND), pen_rss=(NPEN,),
                pr_w=(len(SCALES), ND), pr_c=(len(SCALES), ND), ctrl=(2, 4), v2=(3, 4), secs=(5,))


def trial(args):
    """`_trial`, an exception recorded as a failed trial (NaN, the message printed) instead of ending the run."""
    try:
        return _trial(args)
    except Exception as e:                               # noqa: BLE001 - recorded, counted in the raw and reported
        print('FAILED trial %r: %s: %s' % (args, type(e).__name__, e), flush=True)
        out = {k: np.full(v, np.nan) for k, v in _shapes().items()}
        out['meta'] = np.array([args[1], args[6], args[3], args[4], args[5], args[2], np.nan])
        return out


def run_stage(name, workers=1, out_path=None, chunk=24, jobs=None):
    """Run the trials of one stage ('sanity', 'dev' or 'main'), saving the raw after every chunk (a stopped run resumes from the partial file)."""
    out_path = out_path or RAW
    code, gen = STAGES[name]
    jobs = jobs if jobs is not None else gen()
    part = out_path + '.%s.part.npz' % name
    done = []
    if os.path.exists(part):
        z = np.load(part, allow_pickle=False)
        n0 = int(z['n'])
        done = [{k: z[k][i] for k in _shapes()} for i in range(n0)]
        print('%s: resuming after %d of %d trials' % (name, n0, len(jobs)), flush=True)
    t0 = time.time()
    pool = None
    if workers > 1:
        from multiprocessing import Pool
        pool = Pool(workers)
    while len(done) < len(jobs):
        batch = jobs[len(done):len(done) + chunk]
        done += (pool.map(trial, batch, chunksize=1) if pool else [trial(a) for a in batch])
        np.savez_compressed(part, n=len(done), **{k: np.stack([d[k] for d in done]) for k in _shapes()})
        print('%s: %d / %d trials, %.0f s' % (name, len(done), len(jobs), time.time() - t0), flush=True)
    if pool:
        pool.close()
    arrs = {'%s_%s' % (name, k): np.stack([d[k] for d in done]) for k in _shapes()}
    arrs['%s_wall' % name] = np.array(time.time() - t0)
    if name == 'main':
        ref = np.load(REF_RAW, allow_pickle=False)
        rm = ref['meta']
        mm = np.stack([d['meta'] for d in done])
        assert np.array_equal(rm[:, :4], mm[:, [0, 5, 2, 3]]), 'the main design is not the stored comparison\'s job list'
        arrs.update(ref_res=ref['res'], ref_cols=ref['cols'], ref_meta=rm)
    old = dict(np.load(out_path, allow_pickle=False)) if os.path.exists(out_path) else {}
    old.update(arrs)
    np.savez_compressed(out_path, **old)
    os.remove(part)
    print('%s: %d trials in %.0f s' % (name, len(jobs), time.time() - t0), flush=True)


# ---------------------------------------------------------------------------------------------------------------------------
# implementation checks (stage 'checks')
# ---------------------------------------------------------------------------------------------------------------------------
def checks():
    """Returns dict of small arrays: the gradient of `train_net` against finite differences, the planted fixed point, the convergence of `pen` and `prune`."""
    Kd, grid = X.setup()
    rng = np.random.default_rng(7)
    out = {}
    # (ii) gradient: one step of train_net (lam = 0) against the finite-difference gradient of L = 1/2 mean r^2 in w and in c
    ct, at = np.array([-0.55, 0.65]), np.array([1.0, -0.7])
    y = Kd.synth(ct[:, None], at) + 0.01 * rng.standard_normal(Kd.N)
    c0 = grid[:, None].copy()
    w0 = 0.3 * rng.standard_normal(ND)

    def loss(w, c):
        r = y - Kd.synth(c, w)
        return 0.5 * float(r @ r) / Kd.N
    h = 1e-6
    gw = np.array([(loss(w0 + h * np.eye(ND)[j], c0) - loss(w0 - h * np.eye(ND)[j], c0)) / (2 * h) for j in range(ND)])
    gc = np.array([(loss(w0, c0 + h * np.eye(ND)[j][:, None]) - loss(w0, c0 - h * np.eye(ND)[j][:, None])) / (2 * h) for j in range(ND)])
    w1, c1 = B.train_net(Kd, y, c0, w0, lam=0.0, steps=1)
    Phi = Kd.phi(c0)
    L = max(float(np.linalg.eigvalsh(Phi.T @ Phi / Kd.N).max()), 1e-12)
    lw, lc = 1.0 / L, 0.1 / L
    out['chk_grad'] = np.array([np.abs((w1 - w0) / (-lw) - gw).max() / np.abs(gw).max(), np.abs((c1 - c0).ravel() / (-lc) - gc).max() / np.abs(gc).max()])
    # (iii) the planted solution is a fixed point (noiseless, no penalty)
    f = Kd.synth(ct[:, None], at)
    wf, cf = B.train_net(Kd, f, ct[:, None], at, lam=0.0, steps=300)
    out['chk_fixed'] = np.array([np.abs(wf - at).max(), np.abs(cf.ravel() - ct).max()])
    # (iv) convergence on three planted noisy cases: pen (lam = PEN_LAMS[2]) 2500 against 20000 steps; prune (scale 0.3) 5000 against 50000 steps
    cv, cp = [], []
    for k, (cc, aa) in enumerate((([0.3], [1.0]), ([-0.55, 0.65], [1.0, 0.8]), ([-0.6, 0.0, 0.6], [1.0, -0.8, 0.9]))):
        cc, aa = np.array(cc), np.array(aa)
        yy = Kd.synth(cc[:, None], aa) + 0.01 * np.random.default_rng(100 + k).standard_normal(Kd.N)
        w0 = 0.1 * np.random.default_rng(k).standard_normal(ND)
        row = []
        for steps in (STEPS_PEN, 8 * STEPS_PEN):
            w, c = B.train_net(Kd, yy, grid[:, None], w0, lam=PEN_LAMS[2], steps=steps)
            r = yy - Kd.synth(c, w)
            row += [0.5 * float(r @ r) / Kd.N + PEN_LAMS[2] * float(np.abs(w).sum()), float(_nz(w, 1e-6).sum())]
        cv.append(row)
        w0 = 0.3 * np.random.default_rng(k).standard_normal(ND)
        row = []
        for steps in (STEPS_PRUNE, 10 * STEPS_PRUNE):
            w, c = B.train_net(Kd, yy, grid[:, None], w0, lam=0.0, steps=steps)
            r = yy - Kd.synth(c, w)
            row += [float(r @ r) / Kd.N, float((np.abs(w) >= 0.1 * np.abs(w).max()).sum())]
        cp.append(row)
    out['chk_conv'], out['chk_prune'] = np.array(cv), np.array(cp)
    return out


def run_checks(out_path=None):
    out_path = out_path or RAW
    t0 = time.time()
    arrs = checks()
    arrs['chk_wall'] = np.array(time.time() - t0)
    old = dict(np.load(out_path, allow_pickle=False)) if os.path.exists(out_path) else {}
    old.update(arrs)
    np.savez_compressed(out_path, **old)
    print('checks in %.0f s' % (time.time() - t0), flush=True)


# ---------------------------------------------------------------------------------------------------------------------------
# scoring a stage from its raw
# ---------------------------------------------------------------------------------------------------------------------------
class Stage:
    """The arrays of one stage ('sanity', 'dev', 'main') of the raw, the scoring of every method / counting convention / setting from them (cached)."""

    def __init__(self, z, name):
        self.name = name
        g = lambda k: z['%s_%s' % (name, k)]                      # noqa: E731
        self.a = {k: g(k) for k in _shapes()}
        self.bad = np.isnan(self.a['meta'][:, 6])
        self.K = self.a['meta'][:, 0].astype(int)
        self.seed = self.a['meta'][:, 2].astype(int)
        self.noise = self.a['meta'][:, 1]
        self.plc = self.a['meta'][:, 4].astype(int)
        self.n = len(self.K)
        self.N = X.NS
        self.cache = {}

    def truth(self, i):
        k = self.K[i]
        return self.a['ct'][i, :k], self.a['at'][i, :k]

    def score(self, method, merged, cfg):
        key = (method, merged, repr(sorted(cfg.items())))
        if key in self.cache:
            return self.cache[key]
        a = self.a
        out = np.full((self.n, 4), np.nan)
        for i in range(self.n):
            if self.bad[i]:
                continue
            if method == 'step':
                kh, W, C = units_step(a['step_idx'][i], a['step_w'][i], a['step_rss'][i], cfg, merged, self.N)
            elif method == 'pen':
                kh, W, C = units_pen(a['pen_w'][i], a['pen_c'][i], a['pen_rss'][i], cfg, merged, self.N)
            else:
                kh, W, C = units_prune(a['pr_w'][i], a['pr_c'][i], cfg, merged)
            ct, at = self.truth(i)
            out[i] = X._score(kh, W, C[:, None], ct, at)
        self.cache[key] = out
        return out

    def exact(self, r):
        """(exact order, exact and placed) as boolean arrays from a score array (n, 4)."""
        e = r[:, 0] == self.K
        return e, e & (r[:, 1] == self.K)


def exact_rate(st, method, merged, cfg):
    r = st.score(method, merged, cfg)
    ok = ~st.bad
    return float((r[ok, 0] == st.K[ok]).mean())


def tune(dev):
    """The development-set choice for every baseline and counting convention: {(method, merged): (cfg, dev rate, default's dev rate)}."""
    out = {}
    for m in METHODS:
        for merged in (False, True):
            best, br = None, -1.0
            for cfg in config_grid(m, merged):
                r = exact_rate(dev, m, merged, cfg)
                if r > br:
                    best, br = cfg, r
            out[(m, merged)] = (best, br, exact_rate(dev, m, merged, DEFAULT[m][1 if merged else 0]))
    return out


# ---------------------------------------------------------------------------------------------------------------------------
# statistics and markdown
# ---------------------------------------------------------------------------------------------------------------------------
def cluster_rate_ci(x, tid, B_=2000, seed=0):
    """Mean over targets of the target's rate and its 95 % percentile bootstrap interval over the targets (the trials of a target resampled together)."""
    ids, inv = np.unique(tid, return_inverse=True)
    means = np.bincount(inv, weights=np.asarray(x, float)) / np.bincount(inv)
    draw = np.random.default_rng(seed).integers(0, len(ids), size=(B_, len(ids)))
    boot = means[draw].mean(axis=1)
    return float(means.mean()), float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def pct(k, n):
    lo, hi = X.wilson(k, n)
    lo, hi = max(lo, 0.0), min(hi, 1.0)                  # rounding must not print "-0"
    return '%.0f %% [%.0f, %.0f]' % (100 * k / n, 100 * lo, 100 * hi) if n else 'n/a'


def table(hdr, rows):
    return '| ' + ' | '.join(hdr) + ' |\n|' + '|'.join(['---'] * len(hdr)) + '|\n' + '\n'.join('| ' + ' | '.join(r) + ' |' for r in rows) + '\n'


def sanity_rows(st, tuned):
    """[(label, (n, 4) score array)] for the sanity tables: the three baselines in V0, V1, V2 at the comparison's settings, the tuned settings, and the controls."""
    rows = []
    for mi, m in enumerate(METHODS):
        rows.append(('%s, V0 (the comparison\'s count)' % MLABEL[m], st.score(m, False, DEFAULT[m][0])))
        rows.append(('%s, V1 (merged count)' % MLABEL[m], st.score(m, True, DEFAULT[m][1])))
        rows.append(('%s, V2 (merged, then slid)' % MLABEL[m], st.a['v2'][:, mi, :]))
    if tuned:
        for m in METHODS:
            for merged in (False, True):
                rows.append(('%s, %s, tuned on the development set' % (MLABEL[m], 'V1' if merged else 'V0'), st.score(m, merged, tuned[(m, merged)][0])))
    for j, c in enumerate(CTRL):
        rows.append((CTRL_LABEL[c], st.a['ctrl'][:, j, :]))
    return rows


def failure_rows(st, sel, rows):
    hdr = ['method and count', 'trials', 'exact order', 'exact and placed', 'over-estimates', 'under-estimates', 'right count, wrong place', 'mean K-hat - K', 'median K-hat']
    out = []
    for lab, r in rows:
        r = r[sel]
        K = st.K[sel]
        n = len(K)
        if n == 0:
            out.append([lab, '0'] + ['n/a'] * 7)
            continue
        e = r[:, 0] == K
        pl = e & (r[:, 1] == K)
        out.append([lab, str(n), pct(int(e.sum()), n), pct(int(pl.sum()), n), pct(int((r[:, 0] > K).sum()), n), pct(int((r[:, 0] < K).sum()), n), pct(int((e & ~pl).sum()), n),
                    '%+.2f' % float(np.mean(r[:, 0] - K)), '%g' % float(np.median(r[:, 0]))])
    return hdr, out


def summary_lines(z, have, SUM):
    """The summary section: the exact-order rate of each baseline under each count on the sanity trials and on the main design, with the readings computed from the same arrays."""
    L = ['## Summary (every number computed from the raw below)\n']
    if 'sanity' in SUM and 'main' in SUM:
        sn, ok_s, tuned, n_over, n_cured, ratio, n_opt, n_J = SUM['sanity']
        mn, ok_m, keyed, cdd = SUM['main']
        on, off = ok_s & (sn.plc == 0), ok_s & (sn.plc > 0)

        def srate(m, merged, msk):
            r = sn.score(m, merged, DEFAULT[m][1 if merged else 0])
            return pct(int((r[msk, 0] == sn.K[msk]).sum()), int(msk.sum()))

        def mrate(key):
            r = keyed[key][0]
            return pct(int(((r[:, 0] == mn.K) & ok_m).sum()), int(ok_m.sum()))
        hdr = ['baseline', 'sanity, on grid: V0', 'sanity, off grid: V0', 'sanity, on grid: V1', 'sanity, off grid: V1', 'main: V0 (stored)', 'main: V1', 'main: V0 tuned', 'main: V1 tuned', 'main: V1 minus sliding ISS (BIC), points']
        body = []
        for m, c in zip(METHODS, ('step|bic', 'pen|bic', 'prune')):
            d, lo, hi = cdd[('V1', m)]
            body.append([MLABEL[m], srate(m, False, on), srate(m, False, off), srate(m, True, on), srate(m, True, off), mrate(('ref', c)), mrate(('V1', m)), mrate(('V0t', m)), mrate(('V1t', m)),
                         '%+.1f [%+.1f, %+.1f]' % (-100 * d, -100 * hi, -100 * lo)])
        body.append(['sliding ISS, BIC (stored)', '', '', '', '', mrate(('ref', 'iss|bic')), '', '', '', ''])
        body.append(['matching pursuit + sliding, BIC (stored; sanity: positive control)', pct(int((sn.a['ctrl'][on, 0, 0] == sn.K[on]).sum()), int(on.sum())),
                     pct(int((sn.a['ctrl'][off, 0, 0] == sn.K[off]).sum()), int(off.sum())), '', '', mrate(('ref', 'omp|bic')), '', '', '', ''])
        L.append(table(hdr, body))
        c0 = z['chk_conv']
        L.append('- **Implementation (a).** The frozen code is reproduced: see the reproduction line below. The training\'s update equals the finite-difference gradient to %.1e; the trained penalised '
                 'network reaches an objective no larger than the planted truth\'s (within 0.1 %%) in %d of %d (trial, penalty) cases; 8 times more steps change its objective by at most '
                 '%.1e (relative) and the unit counts not at all (%d of %d cases).' % (max(float(z['chk_grad'][0]), float(z['chk_grad'][1])), n_opt, n_J,
                                                                                    float(np.max(np.abs(c0[:, 2] - c0[:, 0]) / c0[:, 0])), int((c0[:, 1] == c0[:, 3]).sum()), len(c0)))
        L.append('- **Design (b), penalised network.** On the %d sanity trials it over-estimates in %d (its units sit on top of each other at the right centres) and %d of these are exact once coincident '
                 'units are merged; on the main design the exact-order rate goes from %s (stored) to %s (merged count), against %s for the sliding ISS (BIC).' % (
                     int(ok_s.sum()), n_over, n_cured, mrate(('ref', 'pen|bic')), mrate(('V1', 'pen')), mrate(('ref', 'iss|bic'))))
        L.append('- **Design (b), forward selection.** The right unit count is reached in %s of the off-grid sanity trials (%s on the grid) because it needs several neighbouring kernels per '
                 'feature; merged, the sanity rate is %s and the main-design rate %s against %s for the sliding ISS (BIC): a gap that the merged count does not close.' % (
                     srate('step', False, off), srate('step', False, on), srate('step', True, off), mrate(('V1', 'step')), mrate(('ref', 'iss|bic'))))
        L.append('- **Random start + prune.** The comparison\'s count is right in %s of the sanity trials with the features on the grid, noiseless included, with a median of %g of %d units surviving; '
                 'with the units merged the main-design rate is %s.' % (srate('prune', False, on), float(np.median(sn.score('prune', False, DEFAULT['prune'][0])[ok_s, 0])), ND, mrate(('V1', 'prune'))))
        L.append('- **Hyper-parameters (c).** Retuned on development targets and scored on the main design: V0 count %s / %s / %s (forward selection / penalised / prune), merged count %s / %s / %s.\n' % (
            mrate(('V0t', 'step')), mrate(('V0t', 'pen')), mrate(('V0t', 'prune')), mrate(('V1t', 'step')), mrate(('V1t', 'pen')), mrate(('V1t', 'prune'))))
    return L


def report(path=RAW):
    z = np.load(path, allow_pickle=False)
    have = {k: any(f.startswith(k + '_') for f in z.files) for k in ('sanity', 'dev', 'main', 'chk')}
    L = ['# Are the order-selection baselines a strawman? Forward selection, penalised training and random start + prune, checked on cases with a known answer', '',
         'Generated by `experiments/exp_baselines_sanity.py` from `baselines_sanity_raw.npz`; every number below is printed by the script; the analyses were declared in the script\'s docstring '
         'before the runs. The comparison under test is `results/chaos/noise_order.md` (432 trials, `experiments/exp_noise_order.py`); the frozen code is imported, not edited. '
         'Intervals are 95 %% Wilson intervals over trials unless a bootstrap over targets is named. Counting conventions: V0 = the comparison\'s (units counted one by one); '
         'V1 = units merged into clusters (single linkage, link %g = one width; clusters below %g of the largest mass dropped) and clusters counted; V2 = V1, then the clusters\' '
         'masses and centres slid jointly on least squares and units that coincide merged again. "Exact" is K-hat = K; "exact and placed" also needs all K true features matched '
         'within half a width.\n' % (LINKS[0], THMS[0])]
    dev = Stage(z, 'dev') if have['dev'] else None
    tuned = tune(dev) if dev is not None else None
    SUMMARY = {}
    L.append('## What each baseline does in the comparison\n')
    L.append('- **Forward selection (`step`).** Least squares on the 21 grid kernels (centres stay on the grid), greedy, at most %d units, the block by BIC with one parameter per unit; '
             'K-hat = the number of grid kernels chosen. It cannot place a centre between grid points, so an off-grid feature is fitted by several neighbouring kernels.' % MAXU)
    L.append('- **Penalised network (`pen`).** 21 units, centres trainable (start on the grid), coefficients start at 0.1 N(0,1); proximal gradient descent, %d steps, at %d penalties '
             '%.0e ... %.0e, the penalty by BIC with two parameters per nonzero coefficient; K-hat = the number of nonzero coefficients. The l1 norm of the coefficients does not '
             'change when a unit is split into two at the same centre, so units that sit at one centre are counted one by one.' % (STEPS_PEN, NPEN, PEN_LAMS[0], PEN_LAMS[-1]))
    L.append('- **Random start + prune (`prune`).** The same network without a penalty from a random start of scale %g, %d steps, units with |w| below 0.1 max|w| removed; K-hat = '
             'the survivors. Nothing in the training sparsifies the coefficients.\n' % (SCALES[0], STEPS_PRUNE))
    idx_sum = len(L)
    if have['chk']:
        L.append('## A0. Implementation checks\n')
        g, fx, cv, cp = z['chk_grad'], z['chk_fixed'], z['chk_conv'], z['chk_prune']
        L.append('- Gradient of `train_net` against central finite differences of the loss (relative to the largest component), one step without penalty from a random start: '
                 'coefficients %.1e, centres %.1e.' % (g[0], g[1]))
        L.append('- The planted solution (K = 2, noiseless, no penalty) after 300 steps: largest change of a coefficient %.1e, of a centre %.1e.' % (fx[0], fx[1]))
        rows = [['penalised, lam = %.3g, %d vs %d steps (case %d)' % (PEN_LAMS[2], STEPS_PEN, 8 * STEPS_PEN, i + 1), 'objective %.5g -> %.5g' % (r[0], r[2]), '%d -> %d' % (r[1], r[3])]
                for i, r in enumerate(cv)]
        rows += [['random start + prune, %d vs %d steps (case %d)' % (STEPS_PRUNE, 10 * STEPS_PRUNE, i + 1), 'residual mean square %.3g -> %.3g' % (r[0], r[2]), '%d -> %d' % (r[1], r[3])]
                 for i, r in enumerate(cp)]
        L.append('\nConvergence of the training (three planted cases, K = 1, 2, 3, s = 0.01; unit counts at the comparison\'s rule):\n')
        L.append(table(['training', 'objective or residual', 'units counted'], rows))
    if have['main']:
        mn = Stage(z, 'main')
        cols = [str(c) for c in z['ref_cols']]
        ref = z['ref_res']
        ok = ~mn.bad
        L.append('- Reproduction of the comparison: the frozen baselines were re-run on the %d trials of the main design (%d failed). K-hat identical to the stored raw in: '
                 'forward selection %d / %d, penalised network %d / %d, random start + prune %d / %d trials.\n' % (
                     mn.n, int(mn.bad.sum()), int((mn.score('step', False, DEFAULT['step'][0])[ok, 0] == ref[ok, cols.index('step|bic'), 0]).sum()), int(ok.sum()),
                     int((mn.score('pen', False, DEFAULT['pen'][0])[ok, 0] == ref[ok, cols.index('pen|bic'), 0]).sum()), int(ok.sum()),
                     int((mn.score('prune', False, DEFAULT['prune'][0])[ok, 0] == ref[ok, cols.index('prune'), 0]).sum()), int(ok.sum())))
    if have['sanity']:
        sn = Stage(z, 'sanity')
        ok = ~sn.bad
        L.append('## A1. Sanity cases with a known answer (%d trials, %d failed)\n' % (sn.n, int(sn.bad.sum())))
        L.append('K = 1, and K = 2 with the centres at least ten widths apart (%.2f); placement on the grid, half a grid step off it, or uniform at random; s = 0 (noiseless), 0.001, 0.01 of the '
                 'largest amplitude; %d targets per (K, placement) cell, 1 draw if noiseless, 2 otherwise. Positive controls: the comparison\'s own methods with free centres.\n' % (SAN_SEP, SAN_T))
        rows_all = sanity_rows(sn, tuned)
        cells = [(K, p) for K in (1, 2) for p in range(3)]
        hdr = ['method and count'] + ['K = %d, %s' % (K, PLC[p]) for K, p in cells]
        body = []
        for lab, r in rows_all:
            row = [lab]
            for K, p in cells:
                m = ok & (sn.K == K) & (sn.plc == p)
                row.append(pct(int((r[m, 0] == K).sum()), int(m.sum())))
            body.append(row)
        L.append('### Exact order by case (all noise levels pooled)\n')
        L.append(table(hdr, body))
        nz = sorted(set(sn.noise.tolist()))
        hdr = ['method and count'] + ['s = %g, %s' % (v, g) for v in nz for g in ('on grid', 'off grid')]
        body = []
        for lab, r in rows_all:
            row = [lab]
            for v in nz:
                for g in ('on grid', 'off grid'):
                    m = ok & (sn.noise == v) & ((sn.plc == 0) if g == 'on grid' else (sn.plc > 0))
                    row.append(pct(int((r[m, 0] == sn.K[m]).sum()), int(m.sum())))
            body.append(row)
        L.append('### Exact order by noise level, on-grid against off-grid (K = 1 and 2 pooled; off grid = half step and random)\n')
        L.append(table(hdr, body))
        base_rows = [(lab, r) for lab, r in rows_all if 'tuned' not in lab]
        for title, sel in (('on-grid features', ok & (sn.plc == 0)), ('off-grid features (half step and random)', ok & (sn.plc > 0))):
            hdr, body = failure_rows(sn, sel, base_rows)
            L.append('### Failure modes, %s (all noise levels; K = 1 and 2)\n' % title)
            L.append(table(hdr, body))
        # sensitivity to the link
        hdr = ['method (V1 at the comparison\'s settings)'] + ['link %g' % lk for lk in sorted(LINKS)]
        body = []
        for m in METHODS:
            row = [MLABEL[m]]
            for lk in sorted(LINKS):
                cfg = dict(DEFAULT[m][1], link=lk)
                r = sn.score(m, True, cfg)
                row.append('on %s; off %s' % (pct(int((r[ok & (sn.plc == 0), 0] == sn.K[ok & (sn.plc == 0)]).sum()), int((ok & (sn.plc == 0)).sum())).split(' [')[0],
                                              pct(int((r[ok & (sn.plc > 0), 0] == sn.K[ok & (sn.plc > 0)]).sum()), int((ok & (sn.plc > 0)).sum())).split(' [')[0]))
            body.append(row)
        L.append('### Sensitivity of the merged count to the link (exact order, on-grid / off-grid, all sanity trials)\n')
        L.append(table(hdr, body))
        # optimality of the penalised network: J against the planted truth's
        J, T = [], []
        Kd, _ = X.setup()
        for i in np.flatnonzero(ok):
            ct, at = sn.truth(i)
            f = Kd.synth(ct[:, None], at)
            s = sn.a['meta'][i, 6]
            y = f + s * np.random.default_rng(sn.seed[i] * 1000 + int(sn.a['meta'][i, 3]) + 1).standard_normal(Kd.N) if s > 0 else f
            for j, lam in enumerate(PEN_LAMS):
                w, c = sn.a['pen_w'][i, j], sn.a['pen_c'][i, j]
                r = y - Kd.synth(c[:, None], w)
                rt = y - f
                J.append(0.5 * float(r @ r) / Kd.N + lam * float(np.abs(w).sum()))
                T.append(0.5 * float(rt @ rt) / Kd.N + lam * float(np.abs(at).sum()))
        J, T = np.array(J), np.array(T)
        L.append('### Optimality of the penalised network\n')
        n_opt = int((J <= T * (1 + 1e-3) + 1e-9).sum())
        r0, r1 = sn.score('pen', False, DEFAULT['pen'][0]), sn.score('pen', True, DEFAULT['pen'][1])
        o = ok & (r0[:, 0] > sn.K)                                              # the penalised network over-estimates with the comparison's count
        n_cured = int((o & (r1[:, 0] == sn.K)).sum())
        L.append('For every sanity trial and each of the %d penalties, J = 1/2 mean r^2 + lam |w|_1 of the trained network against J of the planted truth (same data, same lam): '
                 'J(network) <= J(truth) (1 + 1e-3) + 1e-9 in %d of %d cases; the largest ratio J(network)/J(truth) is %.3f.\n' % (
                     NPEN, n_opt, len(J), float((J / np.maximum(T, 1e-12)).max())))
        L.append('Penalised network, over-estimating trials (V0 K-hat > K): %d of %d sanity trials; %d of them are exact once the units are merged (V1), i.e. the excess is units at '
                 'the same centre or next to it.\n' % (int(o.sum()), int(ok.sum()), n_cured))
        SUMMARY['sanity'] = (sn, ok, tuned, int(o.sum()), n_cured, float(np.max((J / np.maximum(T, 1e-12)))), n_opt, len(J))
        L.append('Random start + prune: median number of survivors at the comparison\'s rule (thr 0.1) over the sanity trials: %g of %d units.\n' % (
            float(np.median(sn.score('prune', False, DEFAULT['prune'][0])[ok, 0])), ND))
    if dev is not None:
        L.append('## A3. Hyper-parameters retuned on development targets (%d trials, seeds 1000 K + 100 + t, t < %d)\n' % (dev.n, DEV_T))
        rows = []
        for m in METHODS:
            for merged in (False, True):
                cfg, br, dr = tuned[(m, merged)]
                rows.append([MLABEL[m], 'V1' if merged else 'V0', cfg_str(m, cfg), '%.1f %%' % (100 * dr), '%.1f %%' % (100 * br), str(len(config_grid(m, merged)))])
        L.append(table(['baseline', 'count', 'chosen setting', 'development exact order, comparison\'s setting', 'development exact order, chosen', 'settings tried'], rows))
    if have['main']:
        mn = Stage(z, 'main')
        ok = ~mn.bad
        nl = len(X.NOISE)
        tid = X.target_ids(mn.K, mn.seed)
        ni = mn.a['meta'][:, 5].astype(int)
        rowsM = [(lab, ref[:, cols.index(c), :], ('ref', c)) for c, lab in REF_ROWS]
        for m in METHODS:
            rowsM.append(('%s, V1 (merged count)' % MLABEL[m], mn.score(m, True, DEFAULT[m][1]), ('V1', m)))
        for mi, m in enumerate(METHODS):
            rowsM.append(('%s, V2 (merged, then slid)' % MLABEL[m], mn.a['v2'][:, mi, :], ('V2', m)))
        for m in METHODS:
            for merged in (False, True):
                if tuned:
                    rowsM.append(('%s, %s, tuned on development targets' % (MLABEL[m], 'V1' if merged else 'V0'), mn.score(m, merged, tuned[(m, merged)][0]), ('V1t' if merged else 'V0t', m)))
        L.append('## A4. The main comparison (%d trials, %d failed; the design of `exp_noise_order.py`: K = 1..6, 6 targets per K, 3 draws, four noise levels) with the corrected count\n' % (mn.n, int(mn.bad.sum())))
        L.append('Rows marked (stored) are the comparison\'s own numbers from `noise_order_raw.npz`; the others were computed here on the same trials. The last column has the 95 %% bootstrap '
                 'interval over the %d targets.\n' % len(np.unique(tid)))
        hdr = ['method and count'] + ['s = %g (n = %d)' % (v, int((ni == i).sum())) for i, v in enumerate(X.NOISE)] + ['all, Wilson', 'all, bootstrap over targets']
        body = []
        for lab, r, src in rowsM:
            e = (r[:, 0] == mn.K) & ok
            row = [lab]
            for i in range(nl):
                m = ok & (ni == i)
                row.append(pct(int(e[m].sum()), int(m.sum())))
            row.append(pct(int(e[ok].sum()), int(ok.sum())))
            mu, lo, hi = cluster_rate_ci(e[ok], tid[ok])
            row.append('%.0f %% [%.0f, %.0f]' % (100 * mu, 100 * lo, 100 * hi))
            body.append(row)
        L.append('### Exact order, by noise level\n')
        L.append(table(hdr, body))
        hdr = ['method and count'] + ['K = %d' % k for k in X.KS]
        body = [[lab] + [('%.0f %%' % (100 * ((r[:, 0] == mn.K) & ok)[mn.K == k].mean())) if (mn.K == k).any() else 'n/a' for k in X.KS] for lab, r, src in rowsM]
        L.append('### Exact order, by the true number of kernels (all noise levels)\n')
        L.append(table(hdr, body))
        hdr = ['method and count', 'exact and placed', 'over-estimates', 'under-estimates', 'mean K-hat - K', 'position rms (exact and placed, units of the width)']
        body = []
        for lab, r, src in rowsM:
            e = (r[:, 0] == mn.K) & ok
            pl = e & (r[:, 1] == mn.K)
            d = r[ok, 0] - mn.K[ok]
            pos = float(np.sqrt(np.nanmean(r[pl, 2] ** 2)) / SIG) if pl.any() else float('nan')
            body.append([lab, pct(int(pl.sum()), int(ok.sum())), '%.0f %%' % (100 * (d > 0).mean()), '%.0f %%' % (100 * (d < 0).mean()), '%+.2f' % d.mean(), '%.3f' % pos])
        L.append('### Failure modes and accuracy (all noise levels)\n')
        L.append(table(hdr, body))
        E = np.stack([(r[:, 0] == mn.K) for _, r, _ in rowsM], axis=1)
        jb = [lab for lab, _, _ in rowsM].index('sliding ISS, BIC')
        oth = [j for j in range(len(rowsM)) if j != jb]
        cd = X.target_cluster_diffs(E[ok], tid[ok], jb, oth)
        L.append('### Difference from the sliding ISS (BIC), clustered by target\n')
        L.append('Exact-order rate of the sliding ISS (BIC) minus each row (positive: the ISS is better), points, 95 %% percentile bootstrap over the %d targets (2000 resamples, seed 0).\n' % len(np.unique(tid)))
        L.append(table(['row', 'difference (points)'], [[rowsM[j][0], '%+.1f [%+.1f, %+.1f]' % (100 * d, 100 * lo, 100 * hi)] for j, (d, lo, hi) in zip(oth, cd)]))
        hdr = ['merged count V1 at the comparison\'s other settings'] + ['link %g' % lk for lk in sorted(LINKS)]
        body = []
        for m in METHODS:
            row = [MLABEL[m]]
            for lk in sorted(LINKS):
                r = mn.score(m, True, dict(DEFAULT[m][1], link=lk))
                row.append(pct(int(((r[:, 0] == mn.K) & ok).sum()), int(ok.sum())))
            body.append(row)
        L.append('### Sensitivity of the merged count to the link, main design (exact order, all trials; added after the first look at the sanity trials, the link 0.12 was declared in advance)\n')
        L.append(table(hdr, body))
        keyed = {key: (r, j) for j, (lab, r, key) in enumerate(rowsM)}
        cdd = {rowsM[j][2]: (d, lo, hi) for j, (d, lo, hi) in zip(oth, cd)}
        SUMMARY['main'] = (mn, ok, keyed, cdd)
    sm = summary_lines(z, have, SUMMARY)
    L[idx_sum:idx_sum] = sm
    tm = []
    for nm in ('sanity', 'dev', 'main'):
        if have[nm]:
            s_ = z['%s_secs' % nm]
            tm.append([nm, str(len(s_)), '%.3f / %.2f / %.2f / %.2f / %.2f' % tuple(float(np.nanmedian(s_[:, i])) for i in (0, 1, 2, 3, 4))])
    if tm:
        L.append('## Cost (median seconds per trial: forward selection / penalised network (7 penalties) / random start + prune (2 scales) / controls / V2 slides; one core, a shared machine)\n')
        L.append(table(['stage', 'trials', 'seconds'], tm))
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--stage', default='checks,sanity,dev,main')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('BS_WORKERS', '1')))
    A = ap.parse_args()
    if A.run:
        os.makedirs(OUT, exist_ok=True)
        for st in A.stage.split(','):
            if st == 'checks':
                run_checks()
            else:
                run_stage(st, A.workers)
    if A.from_raw:
        with open(MD, 'w') as f:
            f.write(report())
        print('wrote', MD)


if __name__ == '__main__':
    main()
