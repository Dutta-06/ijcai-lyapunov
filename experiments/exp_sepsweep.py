"""
Where does the TRAINED network deviate from the sliding path, and how does the separation at which it forms exactly K units compare with the computed threshold?
=============================================================================================================================================================

Two reviewers asked for the same experiment: a separation sweep on the trained factorised network itself (not on the path) -- K Gaussian features at separations 1.5 to 7
widths, amplitude ratios 1 to 100, white noise up to 20 % of the smallest amplitude, many targets, the start eps from 1e-4 to 1e-8 -- with the empirical transition set against
the paper's computed sufficient threshold Delta_entry (4.56-6.70 widths for K <= 3 and the ratios 1-10; `theory3e.md`).  The design, the outcomes, the development rule that fixes the
stopping rule / formed-unit threshold / merge radius, and the registered predictions are in `results/synthetic/sepsweep_plan.md` (committed before any test target was run).

Machinery (all existing): the paper's free-law network `factorized_rbf.PaperLaw` with (m, G, gamma, a_max) = (12, 1, 0.1, 1), centre rate Gm, trained by the exact gradient flow
(`factorized_rbf.train_flow`, LSODA with Radau fallback) from the small start nu_i^2(0) = eps, alpha_i(0) = 0 on a lattice dictionary (spacing 0.1 = 0.833 widths, width 0.12, samples
at spacing 1/30), and the sliding inverse scale space `rich_transit_nd.sliding_iss_nd` (the eps -> 0 path) on the SAME data.  Per target: one path run, one flow per eps; every
candidate configuration of the stopping rule, formed threshold and merge radius is scored from the stored snapshots (the raw keeps the outcome bits, not the snapshots).

    python experiments/exp_sepsweep.py --dev  [--resume]      development targets (seeds DEV_SEED0), writes sepsweep_dev_raw.npz
    python experiments/exp_sepsweep.py --test [--resume]      the held-out run (seeds TEST_SEED0), writes sepsweep_raw.npz       (run once, after the freeze)
    python experiments/exp_sepsweep.py --from-raw             regenerates results/synthetic/sepsweep.md from the raws (byte for byte)
    (SS_OUT redirects the raws and the md; SS_WORK is the scratch directory of the resumable checkpoints; SS_CELLS = a comma list of cell indices, for smoke runs)
"""
import argparse
import json
import math
import os
import sys
import time
import warnings

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

OUT = os.environ.get('SS_OUT', os.path.join(REPO, 'results', 'synthetic'))
WORK = os.environ.get('SS_WORK', os.path.join(os.environ.get('TMPDIR', '/tmp'), 'sepsweep_work'))
DEV_RAW, RAW, MD = (os.path.join(OUT, n) for n in ('sepsweep_dev_raw.npz', 'sepsweep_raw.npz', 'sepsweep.md'))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# the declared design (the plan, sections 3-4)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
SIG = 0.12                                    # the kernel width (the repository's lattice: spacing 0.1 = 0.833 widths, 21 kernels on [-1, 1])
H_LAT, HX = 0.1, 1.0 / 30.0                   # dictionary spacing, sample spacing
MARGIN = 0.5                                  # the data window extends this far (4.2 widths) beyond the outer features, rounded out to the dictionary lattice
LAW = dict(m=12.0, G=1.0, gamma=0.1, a_max=1.0)
A_MIN = 0.7                                   # the smallest amplitude (the theory's unit; T1/theory3b scale it to 0.7)
KS = (2, 3, 4)
DELTAS = (1.5, 2.5, 3.5, 4.5, 5.5, 7.0)       # the separation between neighbouring features, in widths (all pairs at least Delta: neighbours exactly Delta)
RS = (1.0, 2.4, 10.0, 100.0)                  # amplitude ratio largest / smallest
NOISES = (0.0, 0.01, 0.05, 0.2)               # white noise sd as a fraction of the SMALLEST amplitude
EPS = (1e-4, 1e-6, 1e-8)
RS_K4 = (1.0, 10.0)                           # the drop rules (fixed beforehand): K = 4 only at these ratios;
RS_MID = (1.0, 10.0)                          # the intermediate noise levels (0.01, 0.05) only at these ratios (R = 2.4 and 100: noise 0 and 0.2 only)
N_TEST, N_DEV = 40, 6                         # targets per cell
DEV_SEED0, TEST_SEED0 = 7_100_000, 8_200_000
N_T = 601                                     # snapshots of the flow
OMEGA_MAX = 0.5                               # the flow runs to (1 + OMEGA_MAX) tau_K + 1
MAX_SECS = 300                                # a flow that takes longer is a recorded failure (never a success)
MAX_RHS = 600_000                             # budget of the path (a stalled block is a recorded failure)
EXTRA_ENTRIES = 2                             # the path is run to K + 2 entries
LOC_TOL, AMP_TOL = 0.5, 0.5                   # a unit is 'located at' its feature: within 0.5 widths and 50 % of the amplitude (success)
EXACT_POS, EXACT_AMP = 0.05, 0.05             # the paper's exactness tolerances (noiseless secondary outcome)
LEVEL = 0.9                                   # the transition: smallest Delta from which the success rate exceeds 0.9
CONF = 0.95                                   # exact (Clopper-Pearson) intervals of the rates
BOOT, BOOT_SEED = 2000, 0                     # bootstrap resamples of the transition's interval
BIC_FLOOR = 1e-10                             # of mean(y^2), in the residual variance of the BIC

# candidate configurations: stop rule x formed threshold x merge radius (widths); fixed order = tie-break order of the development rule
STOPS = ('A05', 'A25', 'A50', 'BIC')          # A: the snapshot at (1 + omega) tau_K, tau_K = the path's K-th entry time (omega = 0.05, 0.25, 0.5);  BIC: over the snapshots
OMEGAS = {'A05': 0.05, 'A25': 0.25, 'A50': 0.5}
THRS = ('P', 'D')                             # P: min(0.01 max|A*|, 0.1 min|A*|) (the paper's);  D: 1e-3 max|w| of the snapshot (no truth)
RHOS = (0.0, 0.25, 0.5, 0.75)                 # merge radius in widths (single linkage on the formed units' centres)
CFGS = [(s, t, r) for s in STOPS for t in THRS for r in RHOS]
N_CFG = len(CFGS)
D_THETA = 1e-3
FROZEN = ('A25', 'P', 0.0)                    # the frozen configuration (the development rule's choice, see the plan's freeze record); (stop, thr, rho)

# computed thresholds Delta_entry(K, R, d0 = 0.4167) for the Gaussian Gram (Proposition G, blockwise spectral norm, grid 0.02): the five of theory3e.md Table A and the rest by the same code
THRESH = {(2, 1.0): 4.56, (2, 2.4): 5.44, (2, 10.0): 6.70, (2, 100.0): 8.22, (3, 1.0): 5.00, (3, 2.4): 5.80, (3, 10.0): 6.96, (3, 100.0): 8.42,
          (4, 1.0): 5.22, (4, 2.4): 5.98, (4, 10.0): 7.10, (4, 100.0): 8.54}
THRESH_TAB = {(2, 1.0), (2, 2.4), (2, 10.0), (3, 1.0), (3, 2.4)}      # those printed in theory3e.md (the others are computed with the same function for this report)


def cells():
    """The cells (K, iDelta, iR, iNoise) in a fixed order (K, R, noise outer, Delta inner), after the drop rule."""
    out = []
    for K in KS:
        for iR, R in enumerate(RS):
            if K == 4 and R not in RS_K4:
                continue
            for iN, nz in enumerate(NOISES):
                if nz in (0.01, 0.05) and R not in RS_MID:
                    continue
                for iD in range(len(DELTAS)):
                    out.append((K, iD, iR, iN))
    return out


CELLS = cells()


def targets(split, n):
    """All targets of a split as rows (K, iDelta, iR, iNoise, idx) in the order of CELLS, idx = 0..n-1."""
    return [(K, iD, iR, iN, i) for (K, iD, iR, iN) in CELLS for i in range(n)]


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# the setup, the targets
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
_ST = {}


def _cached_bumps():
    """`BumpsND` that keeps the last kernel matrix (the flow's right-hand side asks for phi(c) twice at the same c: bit-identical results, about twice as fast)."""
    from src.core import rich_transit_nd as RN

    class CachedBumps(RN.BumpsND):
        _last = None

        def phi(self, c):
            c = np.atleast_2d(np.asarray(c, float))
            if self._last is not None and self._last[0].shape == c.shape and np.array_equal(self._last[0], c):
                return self._last[1]
            P = super().phi(c)
            self._last = (c.copy(), P)
            return P
    return CachedBumps


def setup(win):
    """(law, kernel object on the window [0, win/10] (samples at spacing 1/30), dictionary centres = the lattice 0, 0.1, ..., win/10) -- cached; win in lattice steps."""
    if win not in _ST:
        from src.core.mirror_free import FreeLaw
        xs = np.linspace(0.0, win * H_LAT, 3 * win + 1)
        c0 = np.linspace(0.0, win * H_LAT, win + 1)
        _ST[win] = (FreeLaw(**LAW), _cached_bumps()(xs[:, None], SIG), c0)
    return _ST[win]


def make_target(K, iD, iR, iN, idx, seed0):
    """The target of a cell: K features with neighbours exactly Delta widths apart, the first at 0.5 + (uniform lattice phase) x 0.1 (the window [0, win/10] is the features' span plus
    MARGIN on both sides, rounded out to the lattice: translation invariance makes the window's place irrelevant, so the random position is the lattice phase), the amplitudes
    A_min R^(j/(K-1)) in a random order (all positive), and the noise vector A_min * noise * z, z ~ N(0, I) the same for every noise level of the cell (paired).
    Returns dict(cs sorted, As, y = noisy samples, f = noiseless samples, win = window in lattice steps)."""
    rng = np.random.default_rng(np.random.SeedSequence([seed0, K, iD, iR, idx]))
    c1 = MARGIN + H_LAT * rng.uniform()
    cs = c1 + np.arange(K) * DELTAS[iD] * SIG
    mags = A_MIN * RS[iR] ** (np.arange(K) / (K - 1))
    As = mags[rng.permutation(K)]
    win = int(math.ceil((cs[-1] + MARGIN) / H_LAT - 1e-9))
    _, Kobj, _ = setup(win)
    z = rng.standard_normal(Kobj.N)
    f = Kobj.synth(cs[:, None], As)
    return dict(cs=cs, As=As, y=f + A_MIN * NOISES[iN] * z, f=f, win=win)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# outcomes of one snapshot
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
def formed_threshold(kind, w, astar):
    """P: min(0.01 max|A*|, 0.1 min|A*|); D: D_THETA max|w|."""
    if kind == 'P':
        a = np.abs(np.asarray(astar, float))
        return min(0.01 * float(a.max()), 0.1 * float(a.min()))
    return D_THETA * float(np.abs(w).max())


def merge_units(w, c, thr, rho):
    """Formed units (|w| > thr) merged by single linkage on the line: neighbours within rho (absolute) join; the cluster's centre is the |w|-weighted mean, its coefficient the sum.
    Returns (centres sorted, coefficients, formed index array)."""
    idx = np.flatnonzero(np.abs(w) > thr)
    if idx.size == 0:
        return np.zeros(0), np.zeros(0), idx
    o = idx[np.argsort(c[idx], kind='stable')]
    cs, ws = c[o], w[o]
    brk = np.flatnonzero(np.diff(cs) > rho) + 1
    cen, amp = [], []
    for g in np.split(np.arange(o.size), brk):
        a = np.abs(ws[g])
        cen.append(float(np.sum(a * cs[g]) / np.sum(a)))
        amp.append(float(np.sum(ws[g])))
    return np.array(cen), np.array(amp), idx


def n_merged(w, c, thr, rho):
    return merge_units(w, c, thr, rho)[0].size


def located(cen, amp, cstar, astar, pos_tol, amp_tol):
    """K merged units, one per feature (features and units both sorted by position), each within pos_tol widths and amp_tol (relative) of its feature."""
    K = len(cstar)
    if cen.size != K:
        return False
    o = np.argsort(cstar)
    return bool(np.all(np.abs(cen - cstar[o]) <= pos_tol * SIG) and np.all(np.abs(amp - astar[o]) <= amp_tol * np.abs(astar[o])))


def snapshot_outcome(w, c, cstar, astar, thr_kind, rho):
    """(success, formed index set (tuple), n formed (raw), n merged, exact) of one snapshot; exact = the paper's tolerances (0.05 widths, 5 %) on the merged units."""
    thr = formed_threshold(thr_kind, w, astar)
    cen, amp, idx = merge_units(w, c, thr, rho * SIG)
    ok = located(cen, amp, cstar, astar, LOC_TOL, AMP_TOL)
    ex = located(cen, amp, cstar, astar, EXACT_POS, EXACT_AMP)
    return ok, tuple(int(j) for j in idx), int(idx.size), int(cen.size), ex


def bic_value(rss, N, k, floor):
    return N * math.log(rss / N + floor) + 2.0 * k * math.log(N)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# the path and the network on one target
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
def run_path(win, K, y):
    """The sliding path on y to K + EXTRA_ENTRIES entries: the list of states (tau, c (d,), w (d,), rms) after 0, 1, ... entries (the final state appended), the entry times, flags."""
    from src.core import rich_transit_nd as RN
    L, Kobj, c0 = setup(win)
    R = RN.sliding_iss_nd(L, Kobj, y, c0[:, None], np.ones(len(c0)), polish=True, max_events=K + EXTRA_ENTRIES, max_rhs=MAX_RHS)
    states = [(float(t), np.array(c)[:, 0].copy(), np.array(w).copy(), float(r)) for (t, c, w, r) in R['hist']]
    ent = [float(e[0]) for e in R['events'] if e[2] == 'enter']
    if not states or not (np.array_equal(states[-1][2], R['w']) and np.array_equal(states[-1][1], R['c'][:, 0])):
        r = y - Kobj.synth(R['c'], R['w'])
        states.append((ent[-1] if ent else 0.0, np.array(R['c'][:, 0]).copy(), np.array(R['w']).copy(), float(np.sqrt(np.mean(r ** 2)))))
    return states, ent, dict(oob=bool(R.get('out_of_budget', False)), failed=int(R.get('solver_failed', 0)), nfev=int(R.get('nfev', 0)),
                             exits=int(sum(e[2] == 'exit' for e in R['events'])))


def index_of_kth_entry(states, ent, K):
    """The first state with at least K entries (state i has the entries up to its tau); returns (i, tau_K, found)."""
    for i, (t, c, w, r) in enumerate(states):
        if sum(e <= t * (1 + 1e-12) + 1e-300 for e in ent) >= K:
            return i, t, True
    return len(states) - 1, states[-1][0], False


def run_network(win, y, eps, tau_end):
    from src.core import factorized_rbf as FR
    L, Kobj, c0 = setup(win)
    ro = FR.PaperLaw(L)
    return FR.train_flow(ro, Kobj, y, c0, eps, tau_end, LAW['G'] * LAW['m'], n_t=N_T, max_secs=MAX_SECS, method='Radau')


def net_rss(Kobj, y, C, W):
    """Residual sum of squares of every snapshot: C (T, d, 1), W (T, d)."""
    out = np.empty(W.shape[0])
    for i in range(W.shape[0]):
        r = y - Kobj.synth(C[i], W[i])
        out[i] = float(np.sum(r * r))
    return out


def evaluate_path(win, K, cs, As, y, states, ent):
    """Per configuration (CFGS order): the path's outcome fields (success, n formed, n merged, exact) and the chosen state's formed set, plus tau_K."""
    _, Kobj, _ = setup(win)
    N = Kobj.N
    floor = BIC_FLOOR * float(np.mean(y ** 2))
    iK, tauK, found = index_of_kth_entry(states, ent, K)
    out, sets = np.zeros((N_CFG, 4), np.int16), []
    for ic, (stop, thr, rho) in enumerate(CFGS):
        if stop == 'BIC':
            vals = [bic_value(N * st[3] ** 2, N, n_merged(st[2], st[1], formed_threshold(thr, st[2], As), rho * SIG), floor) for st in states]
            i = int(np.argmin(vals))
        else:
            i = iK
        st = states[i]
        ok, idx, nf, nm, ex = snapshot_outcome(st[2], st[1], cs, As, thr, rho)
        out[ic] = (ok, nf, nm, ex)
        sets.append(idx)
    return out, sets, tauK, found, iK


def evaluate_network(win, cs, As, y, res, path_sets, tauK):
    """Per configuration: (success, agree with the path's set, n formed, n merged, exact) from the flow's snapshots `res`."""
    _, Kobj, _ = setup(win)
    N = Kobj.N
    W, C, taus = res['W'], res['C'], res['taus']
    floor = BIC_FLOOR * float(np.mean(y ** 2))
    rss = None
    out = np.zeros((N_CFG, 5), np.int16)
    for ic, (stop, thr, rho) in enumerate(CFGS):
        if stop == 'BIC':
            if rss is None:
                rss = net_rss(Kobj, y, C, W)
            vals = [bic_value(rss[i], N, n_merged(W[i], C[i, :, 0], formed_threshold(thr, W[i], As), rho * SIG), floor) for i in range(W.shape[0])]
            i = int(np.argmin(vals))
        else:
            i = int(np.argmin(np.abs(taus - (1.0 + OMEGAS[stop]) * tauK)))
        ok, idx, nf, nm, ex = snapshot_outcome(W[i], C[i, :, 0], cs, As, thr, rho)
        out[ic] = (ok, idx == path_sets[ic], nf, nm, ex)
    return out


STATUS = {'ok': 0}


def status_code(s):
    if s == 'ok':
        return 0
    return 1 if s.startswith('timeout') else 2


def run_target(row, seed0):
    """One target: the path and one flow per eps.  Returns dict(path (N_CFG, 4), net (n_eps, N_CFG, 5), diag floats)."""
    K, iD, iR, iN, idx = row
    tg = make_target(K, iD, iR, iN, idx, seed0)
    cs, As, y, win = tg['cs'], tg['As'], tg['y'], tg['win']
    t0 = time.time()
    states, ent, pinfo = run_path(win, K, y)
    pout, psets, tauK, found, iK = evaluate_path(win, K, cs, As, y, states, ent)
    tpath = time.time() - t0
    tau_end = (1.0 + OMEGA_MAX) * tauK + 1.0
    net = np.zeros((len(EPS), N_CFG, 5), np.int16)
    diag = np.zeros(8 + 3 * len(EPS))
    diag[:8] = (tauK, float(found), tpath, pinfo['nfev'], float(pinfo['oob']), pinfo['failed'], pinfo['exits'], len(ent))
    for ie, eps in enumerate(EPS):
        try:
            res = run_network(win, y, eps, tau_end)
            st = status_code(res['status'])
            net[ie] = evaluate_network(win, cs, As, y, res, psets, tauK) if st == 0 else 0
            diag[8 + 3 * ie:11 + 3 * ie] = (st, res['secs'], res['nfev'])
        except Exception:                                       # a failed flow is a recorded result (never a success)
            diag[8 + 3 * ie:11 + 3 * ie] = (3, 0.0, 0)
    return dict(path=pout, net=net, diag=diag)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# the run (resumable: one checkpoint line per target in the scratch directory), the raw files
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
def load_raw(path):
    if not os.path.exists(path):
        return {}
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def save_raw(path, d):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp.npz'
    np.savez_compressed(tmp, **d)
    os.replace(tmp, path)


def _job(args):
    row, seed0 = args
    warnings.filterwarnings('ignore')
    o = run_target(row, seed0)
    return dict(row=list(row), path=o['path'].tolist(), net=o['net'].tolist(), diag=o['diag'].tolist())


def pending_rows(split, only=None, shard=None):
    """The targets of a split (restricted to the cells `only` and to the cells of a shard (i, n): cell index mod n == i)."""
    n = N_DEV if split == 'dev' else N_TEST
    rows = targets(split, n)
    if only is not None:
        keep = set(CELLS[i] for i in only)
        rows = [r for r in rows if r[:4] in keep]
    if shard is not None:
        keep = set(c for i, c in enumerate(CELLS) if i % shard[1] == shard[0])
        rows = [r for r in rows if r[:4] in keep]
    return sorted(rows, key=lambda r: tier(r[:4]))                  # stable: priority tiers, then the cell order


def tier(cell):
    """The run priority of a cell (what is dropped last): 0 noise 0 and 0.2 (the predictions), 1 noise 0.05, 2 noise 0.01."""
    nz = NOISES[cell[3]]
    return 0 if nz in (0.0, 0.2) else (1 if nz == 0.05 else 2)


def run_split(split, resume=False, only=None, shard=None, workers=1, tag=''):
    """Run the targets of a split, one checkpoint line per target in WORK/ckpt_<split><tag>.jsonl (resumable).  Returns nothing; `assemble` builds the raw."""
    seed0 = DEV_SEED0 if split == 'dev' else TEST_SEED0
    os.makedirs(WORK, exist_ok=True)
    ck = os.path.join(WORK, 'ckpt_%s%s.jsonl' % (split, tag))
    done = set()
    if resume and os.path.exists(ck):
        for line in open(ck):
            done.add(tuple(json.loads(line)['row']))
    elif os.path.exists(ck):
        os.remove(ck)
    rows = [r for r in pending_rows(split, only, shard) if tuple(r) not in done]
    warnings.filterwarnings('ignore')
    t0 = time.time()
    with open(ck, 'a') as fh:
        if workers > 1:
            from multiprocessing import Pool
            with Pool(workers) as pool:
                it = pool.imap_unordered(_job, [(r, seed0) for r in rows], chunksize=1)
                for k, rec in enumerate(it):
                    fh.write(json.dumps(rec) + '\n')
                    fh.flush()
                    if k % 50 == 0:
                        print('%s %d/%d  %.0fs' % (split, k, len(rows), time.time() - t0), flush=True)
        else:
            for k, row in enumerate(rows):
                fh.write(json.dumps(_job((row, seed0))) + '\n')
                fh.flush()
                if k % 20 == 0:
                    print('%s %d/%d  %.0fs' % (split, k, len(rows), time.time() - t0), flush=True)


def assemble(split, files):
    """The raw of a split from checkpoint files (duplicates of a target: the first kept); targets in the order of `targets`; cells with no target at all are absent."""
    recs = {}
    for fn in files:
        for line in open(fn):
            r = json.loads(line)
            recs.setdefault(tuple(r['row']), r)
    order = [r for r in targets(split, N_DEV if split == 'dev' else N_TEST) if tuple(r) in recs]
    rr = [recs[tuple(r)] for r in order]
    raw = dict(rows=np.array([r['row'] for r in rr], np.int16), path=np.array([r['path'] for r in rr], np.int16),
               net=np.array([r['net'] for r in rr], np.int16), diag=np.array([r['diag'] for r in rr], float))
    return {split + '_' + k: v for k, v in raw.items()}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
def cp_interval(k, n, conf=CONF):
    """Exact (Clopper-Pearson) interval of a binomial proportion."""
    from scipy.stats import beta
    a = 1.0 - conf
    lo = 0.0 if k == 0 else float(beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def crossing(rates, grid=DELTAS, level=LEVEL):
    """The smallest grid separation from which the rate exceeds `level` at every larger one, linearly interpolated from the previous grid point.
    Returns (value, flag): flag '<' (already above at the first grid point: value = grid[0], a bound), '>' (never: value = inf), '' (interpolated)."""
    r = list(rates)
    j = len(r)
    while j > 0 and r[j - 1] > level:
        j -= 1
    if j == len(r):
        return float('inf'), '>'
    if j == 0:
        return float(grid[0]), '<'
    r0, r1 = r[j - 1], r[j]
    return float(grid[j - 1] + (level - r0) / (r1 - r0) * (grid[j] - grid[j - 1])), ''


def crossing_many(R, grid=DELTAS, level=LEVEL):
    """`crossing` of every row of R (B, len(grid)) at once (the value only: inf = never, grid[0] = already above at the first point)."""
    R = np.asarray(R, float)
    g = np.asarray(grid, float)
    above = R > level
    suff = np.flip(np.cumprod(np.flip(above, 1), 1), 1)            # suff[:, i]: every point from i on is above the level
    j = R.shape[1] - suff.sum(1).astype(int)                       # the first such i (len: none)
    out = np.empty(R.shape[0])
    never, first = j == R.shape[1], j == 0
    mid = ~never & ~first
    jm = np.where(mid, j, 1)
    r0 = R[np.arange(R.shape[0]), jm - 1]
    r1 = R[np.arange(R.shape[0]), np.minimum(jm, R.shape[1] - 1)]
    with np.errstate(divide='ignore', invalid='ignore'):
        interp = g[jm - 1] + (level - r0) / (r1 - r0) * (g[np.minimum(jm, R.shape[1] - 1)] - g[jm - 1])
    out[:] = np.where(mid, interp, 0.0)
    out[never] = float('inf')
    out[first] = g[0]
    return out


def flag_of(v, grid=DELTAS):
    return '>' if v == float('inf') else ('<' if v == grid[0] else '')


def transition(ks, ns, grid=DELTAS):
    """(estimate, flag, lo, lo_flag, hi, hi_flag): the transition of the success curve and its interval, the 2.5th and 97.5th percentile of the transition of BOOT parametric bootstrap
    resamples of the success counts (each grid point binomial(n, k/n), seed BOOT_SEED)."""
    e = crossing([k / n for k, n in zip(ks, ns)], grid)
    rng = np.random.default_rng(BOOT_SEED)
    p = np.array([k / n for k, n in zip(ks, ns)])
    draws = rng.binomial(np.array(ns)[None, :], p[None, :], size=(BOOT, len(ks)))
    vals = np.sort(crossing_many(draws / np.array(ns)[None, :], grid))
    lo, hi = float(vals[int(0.025 * BOOT)]), float(vals[int(math.ceil(0.975 * BOOT)) - 1])
    return e[0], e[1], lo, flag_of(lo, grid), hi, flag_of(hi, grid)


def fmt_t(v, flag):
    if flag == '>':
        return '>%g' % DELTAS[-1]
    if flag == '<':
        return '<=%g' % DELTAS[0]
    return '%.2f' % v


def fmt_tr(t):
    return '%s [%s, %s]' % (fmt_t(t[0], t[1]), fmt_t(t[2], t[3]), fmt_t(t[4], t[5]))


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# aggregation of a raw split
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
def aggregate(raw, split, cfg):
    """Counts per (K, iR, iN, iDelta): n targets; path successes; for each eps the network successes, agreements, exacts; the discordant pairs (network ok / path not, path ok / network not);
    n formed raw summed...  `cfg` = the configuration (stop, thr, rho).  Returns dict keyed (K, iR, iN, iD) -> dict."""
    ic = CFGS.index(tuple(cfg))
    rows, path, net, diag = raw[split + '_rows'], raw[split + '_path'], raw[split + '_net'], raw[split + '_diag']
    agg = {}
    for t in range(rows.shape[0]):
        K, iD, iR, iN, idx = (int(x) for x in rows[t])
        d = agg.setdefault((K, iR, iN, iD), dict(n=0, path=0, path_exact=0, net=[0] * len(EPS), agree=[0] * len(EPS), exact=[0] * len(EPS), only_net=[0] * len(EPS),
                                                 only_path=[0] * len(EPS), bad=[0] * len(EPS), nofound=0, oob=0, exits=0))
        d['n'] += 1
        pk = int(path[t, ic, 0])
        d['path'] += pk
        d['path_exact'] += int(path[t, ic, 3])
        d['nofound'] += int(diag[t, 1] == 0)
        d['oob'] += int(diag[t, 4] > 0 or diag[t, 5] > 0)
        d['exits'] += int(diag[t, 6] > 0)
        for ie in range(len(EPS)):
            nk = int(net[t, ie, ic, 0])
            d['net'][ie] += nk
            d['agree'][ie] += int(net[t, ie, ic, 1])
            d['exact'][ie] += int(net[t, ie, ic, 4])
            d['only_net'][ie] += int(nk and not pk)
            d['only_path'][ie] += int(pk and not nk)
            d['bad'][ie] += int(diag[t, 8 + 3 * ie] > 0)
    return agg


def curve(agg, K, iR, iN, key, ie=None):
    """The counts over DELTAS of a cell group: (ks, ns)."""
    ks, ns = [], []
    for iD in range(len(DELTAS)):
        d = agg[(K, iR, iN, iD)]
        v = d[key] if ie is None else d[key][ie]
        ks.append(int(v))
        ns.append(int(d['n']))
    return ks, ns


def group_keys(agg):
    """The (K, iR, iN) groups that are complete: every separation present, each with as many targets as the fullest cell (an incomplete group is listed in the report, not analysed)."""
    if not agg:
        return []
    nmax = max(d['n'] for d in agg.values())
    g = set((K, iR, iN) for (K, iR, iN, iD) in agg)
    return sorted(x for x in g if all(x + (iD,) in agg and agg[x + (iD,)]['n'] == nmax for iD in range(len(DELTAS))))


def dev_objective(raw):
    """The development objective of every configuration (amendment 2 of the plan): the mean over ALL development targets of [network success at eps = 1e-8 + path success] / 2.
    Returns the list over CFGS of (J, network part, path part)."""
    path, net = raw['dev_path'], raw['dev_net']
    ie = EPS.index(1e-8)
    out = []
    for ic in range(N_CFG):
        a = float(net[:, ie, ic, 0].mean())
        b = float(path[:, ic, 0].mean())
        out.append((0.5 * (a + b), a, b))
    return out


def select_config(raw, tol=0.005):
    """The development rule: the configuration of maximal objective; among those within `tol` of the maximum the first in CFGS order."""
    obj = dev_objective(raw)
    best = max(o[0] for o in obj)
    return next(ic for ic, o in enumerate(obj) if o[0] >= best - tol)


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# predictions
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
def diff(a, b):
    """a - b of two transitions' values, 0 when both are 'never' (>7)."""
    return 0.0 if (a == float('inf') and b == float('inf')) else a - b


def fmt_diff(d):
    return 'censored (>7)' if abs(d) == float('inf') else '%+.2f' % d


def val(t):
    """A transition as a number for the comparisons: a bound counts as its value, 'never' as +inf."""
    return t[0]


def predictions(agg):
    """P1-P4 on the primary configuration's aggregate.  Returns dict with the per-(K, R) rows and the verdicts."""
    i8, i4 = EPS.index(1e-8), EPS.index(1e-4)
    out = dict(rows=[])
    ok1 = ok2 = ok3 = True
    n1 = n3 = 0
    up, strict, n4 = 0, 0, 0
    gk = set(group_keys(agg))
    for (K, iR) in sorted(set((K, iR) for (K, iR, iN) in gk if (K, iR, 0) in gk and (K, iR, 3) in gk)):
        thr = THRESH[(K, RS[iR])]
        r = dict(K=K, R=RS[iR], thr=thr)
        ks, ns = curve(agg, K, iR, 0, 'net', i8)
        r['net8'] = transition(ks, ns)
        ks, ns = curve(agg, K, iR, 0, 'net', i4)
        r['net4'] = transition(ks, ns)
        ks, ns = curve(agg, K, iR, 0, 'path')
        r['path'] = transition(ks, ns)
        ks, ns = curve(agg, K, iR, 0, 'net', EPS.index(1e-6))
        r['net6'] = transition(ks, ns)
        ks, ns = curve(agg, K, iR, 3, 'net', i8)
        r['net8_n02'] = transition(ks, ns)
        r['p1'] = val(r['net8']) <= thr
        r['p1_untestable'] = False
        r['p2_diff'] = diff(val(r['net8']), val(r['path']))
        r['p2'] = r['p2_diff'] <= 1.0
        r['p3_testable'] = thr <= DELTAS[-1]
        r['p3'] = (val(r['net8_n02']) > thr) if r['p3_testable'] else None
        r['p4_diff'] = diff(val(r['net4']), val(r['net8']))
        ok1 &= r['p1']
        ok2 &= r['p2']
        if r['p3_testable']:
            ok3 &= bool(r['p3'])
            n3 += 1
        n4 += 1
        up += r['p4_diff'] >= 0
        strict += r['p4_diff'] > 0.25
        out['rows'].append(r)
    out['P1'], out['P2'], out['P3'] = ok1, ok2, ok3
    out['P4'] = bool(up == n4 and strict >= (2 * n4 + 2) // 3)
    out['P4_counts'] = (up, strict, n4)
    out['P3_n'] = n3
    return out


# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# the report
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
def table(header, rows):
    s = '| ' + ' | '.join(header) + ' |\n|' + '|'.join('---' for _ in header) + '|\n'
    return s + ''.join('| ' + ' | '.join(str(c) for c in r) + ' |\n' for r in rows)


def pct_cell(k, n):
    lo, hi = cp_interval(k, n)
    return '%d/%d [%.2f, %.2f]' % (k, n, lo, hi)


def frac_cell(k, n):
    return '%d/%d' % (k, n)


def cfg_str(c):
    return '%s, threshold %s, merge radius %g' % c


def report(raws):
    lines = []
    w = lines.append
    test, dev = raws.get('test'), raws.get('dev')
    w('# Separation sweep on the trained network: where it deviates from the path, and the empirical transition against the computed threshold\n')
    w('Generated by `experiments/exp_sepsweep.py` (`--from-raw` regenerates this file byte for byte from `sepsweep_raw.npz` and `sepsweep_dev_raw.npz`); plan, '
      'registered predictions and the freeze: `results/synthetic/sepsweep_plan.md`.  Lengths are in kernel widths (sigma = 0.12; the benchmark lattice has spacing 0.833).  '
      'The law is (m, G, gamma, a_max) = (12, 1, 0.1, 1), the centre rate Gm; the start is nu_i^2 = eps, alpha_i = 0 on a lattice dictionary of spacing 0.833 widths '
      '(the data window of a target is its features\' span plus %g on both sides, rounded out to the lattice, samples at spacing 1/30); the smallest amplitude is %g, the noise is white with sd = '
      '(noise) x (smallest amplitude); the neighbours of every target are exactly Delta apart, the lattice phase is random, all amplitudes are positive.\n' % (MARGIN, A_MIN))
    if dev is not None:
        obj = dev_objective(dev)
        ic = select_config(dev)
        w('## 0. Development (declared before the test; used only for the stopping rule, the formed threshold and the merge radius)\n')
        w('%d development targets (%d per cell, seeds from %d).  The objective of a configuration is the mean over ALL development targets of '
          '[network success at eps = 1e-8 + path success]/2 (amendment 2; the objective first registered, over Delta >= 5.5 only, ties every window rule); '
          'the first configuration within 0.005 of the maximum (in the order listed) is selected.  '
          'Stop rules: A05/A25/A50 = the snapshot at (1 + omega) tau_K with omega = 0.05/0.25/0.5 and tau_K the path\'s K-th entry time (the path\'s own state after its K-th entry); '
          'BIC = the snapshot of least BIC = N ln(rss/N) + 2 k ln N over the recorded snapshots (the path: over its states).  Formed threshold: P = min(0.01 max|A*|, 0.1 min|A*|), '
          'D = 1e-3 max|w|.  Merge radius in widths.\n' % (dev['dev_rows'].shape[0], N_DEV, DEV_SEED0))
        order = sorted(range(N_CFG), key=lambda i: (-obj[i][0], i))[:12]
        w(table(['rank', 'stop', 'threshold', 'merge radius', 'objective', 'network part', 'path part'],
                [[r + 1, CFGS[i][0], CFGS[i][1], '%g' % CFGS[i][2], '%.4f' % obj[i][0], '%.4f' % obj[i][1], '%.4f' % obj[i][2]] for r, i in enumerate(order)]))
        w('Selected by the rule: **%s**.  Frozen in the script and the plan: **%s**.\n' % (cfg_str(CFGS[ic]), cfg_str(FROZEN) if FROZEN else 'not yet frozen'))
    if test is None:
        return '\n'.join(lines) + '\n'
    cfg = FROZEN
    agg = aggregate(test, 'test', cfg)
    n_t = int(test['test_rows'].shape[0])
    w('## 1. Design and size\n')
    w('%d cells (K in %s; Delta in %s; R in %s, K = 4 only R in %s, noise 0.01 and 0.05 only R in %s; noise in %s) of %d targets = %d targets, each run through the path and %d flows (eps = %s).  Primary configuration (frozen): '
      '**%s**.  Success = exactly K merged units, one per feature, each within %g widths and %d %% of its feature (the paper\'s exactness, %g widths and %d %%, is reported for the noiseless '
      'cells as `exact`); agreement = the formed unit set equals the path\'s K-th-entry block (unmerged, same threshold).  The success rate is the share of targets; intervals are exact '
      '(Clopper-Pearson, %d %%).  The empirical transition is the smallest Delta from which the rate exceeds %g at every larger grid point, interpolated linearly between grid points; its '
      'interval is the 2.5th to 97.5th percentile of the transition of %d parametric bootstrap resamples of the success counts (seed %d); `<=1.5` and `>7` are bounds.\n'
      % (len(CELLS), ', '.join(map(str, KS)), ', '.join('%g' % d for d in DELTAS), ', '.join('%g' % r for r in RS), ', '.join('%g' % r for r in RS_K4),
         ', '.join('%g' % r for r in RS_MID), ', '.join('%g' % x for x in NOISES), N_TEST, n_t, len(EPS), ', '.join('%g' % e for e in EPS), cfg_str(cfg), LOC_TOL, int(AMP_TOL * 100),
         EXACT_POS, int(EXACT_AMP * 100), int(CONF * 100), LEVEL, BOOT, BOOT_SEED))
    nbad = int(sum(d['bad'][ie] for d in agg.values() for ie in range(len(EPS))))
    w('Flows that timed out or failed (recorded as failures, counted in the rates): %d of %d.  Paths whose budget ran out or whose block solver failed: %d of %d targets; '
      'paths in which the K-th entry was never reached: %d.\n' % (nbad, n_t * len(EPS), sum(d['oob'] for d in agg.values()), n_t, sum(d['nofound'] for d in agg.values())))
    P = predictions(agg)
    w('## 2. Transitions against the computed thresholds (noiseless)\n')
    w('`computed` = Delta_entry(K, R, d0 = 0.4167) for the Gaussian Gram (Proposition G, blockwise spectral norm, grid 0.02), the number of `theory3e.md` Table A where it is printed (marked *) and '
      'the same function\'s value otherwise.  The thresholds of R = 100 lie beyond the largest separation tested (7), so a transition at or below 7 is below them.\n')
    rows = []
    for r in P['rows']:
        rows.append([r['K'], '%g' % r['R'], '%.2f%s' % (r['thr'], '*' if (r['K'], r['R']) in THRESH_TAB else ''), fmt_tr(r['path']), fmt_tr(r['net8']), fmt_tr(r['net4']), fmt_tr(r['net6'])])
    w(table(['K', 'R', 'computed', 'path', 'network eps 1e-8', 'network eps 1e-4', 'network eps 1e-6'], rows))
    w('## 3. Registered predictions\n')
    w('* **P1** (eps = 1e-8, noise 0: the network\'s transition is at or below the computed threshold for every (K, R)): **%s**.  Point estimates; per (K, R) in the table of section 2 '
      '(%d of %d at or below; upper interval end at or below the threshold in %d).\n'
      % ('HELD' if P['P1'] else 'FAILED', sum(r['p1'] for r in P['rows']), len(P['rows']),
         sum((r['net8'][5] == '' and r['net8'][4] <= r['thr']) or r['net8'][5] == '<' for r in P['rows'])))
    w('* **P2** (the network\'s transition is within +1 width of the path\'s, eps = 1e-8, noise 0): **%s**; differences network minus path: %s.\n'
      % ('HELD' if P['P2'] else 'FAILED', ', '.join('(%d, %g): %s' % (r['K'], r['R'], fmt_diff(r['p2_diff'])) for r in P['rows'])))
    w('* **P3** (noise 0.2 moves the transition above the computed threshold, eps = 1e-8; the %d of the %d (K, R) whose threshold is within the tested range): **%s**; '
      'transitions at noise 0.2: %s.\n' % (P['P3_n'], len(P['rows']), 'HELD' if P['P3'] else 'FAILED',
                                           ', '.join('(%d, %g): %s vs %.2f' % (r['K'], r['R'], fmt_tr(r['net8_n02']), r['thr']) for r in P['rows'] if r['p3_testable'])))
    up, strict, n4 = P['P4_counts']
    w('* **P4** (eps = 1e-4 moves the transition up relative to 1e-8, noiseless; rule: not lower in every (K, R) and higher by more than 0.25 width in at least two thirds): **%s**; '
      'not lower in %d of %d, higher by more than 0.25 in %d; differences: %s.\n'
      % ('HELD' if P['P4'] else 'FAILED', up, n4, strict, ', '.join('(%d, %g): %s' % (r['K'], r['R'], fmt_diff(r['p4_diff'])) for r in P['rows'])))
    w('## 4. Transitions at every noise level (network eps 1e-8 / 1e-6 / 1e-4; path)\n')
    rows = []
    for (K, iR, iN) in group_keys(agg):
        row = [K, '%g' % RS[iR], '%g' % NOISES[iN], '%.2f' % THRESH[(K, RS[iR])]]
        for ie in (2, 1, 0):
            ks, ns = curve(agg, K, iR, iN, 'net', ie)
            row.append(fmt_tr(transition(ks, ns)))
        ks, ns = curve(agg, K, iR, iN, 'path')
        row.append(fmt_tr(transition(ks, ns)))
        rows.append(row)
    w(table(['K', 'R', 'noise', 'computed', 'net 1e-8', 'net 1e-6', 'net 1e-4', 'path'], rows))
    w('## 5. Success curves (successes / targets [exact 95 %% interval]) against Delta\n')
    hdr = ['K', 'R', 'noise', 'who'] + ['%g' % d for d in DELTAS]
    rows = []
    for (K, iR, iN) in group_keys(agg):
        ks, ns = curve(agg, K, iR, iN, 'path')
        rows.append([K, '%g' % RS[iR], '%g' % NOISES[iN], 'path'] + [pct_cell(k, n) for k, n in zip(ks, ns)])
        for ie in (2, 1, 0):
            ks, ns = curve(agg, K, iR, iN, 'net', ie)
            rows.append([K, '%g' % RS[iR], '%g' % NOISES[iN], 'net %g' % EPS[ie]] + [pct_cell(k, n) for k, n in zip(ks, ns)])
    w(table(hdr, rows))
    w('## 6. Agreement with the path (formed set = the path\'s block) and the paper\'s exactness (noiseless cells)\n')
    rows = []
    for (K, iR, iN) in group_keys(agg):
        for ie in (2, 1, 0):
            ks, ns = curve(agg, K, iR, iN, 'agree', ie)
            row = [K, '%g' % RS[iR], '%g' % NOISES[iN], 'agree %g' % EPS[ie]] + [frac_cell(k, n) for k, n in zip(ks, ns)]
            rows.append(row)
        if iN == 0:
            ks, ns = curve(agg, K, iR, iN, 'path_exact')
            rows.append([K, '%g' % RS[iR], '0', 'path exact'] + [frac_cell(k, n) for k, n in zip(ks, ns)])
            ks, ns = curve(agg, K, iR, iN, 'exact', 2)
            rows.append([K, '%g' % RS[iR], '0', 'net exact 1e-8'] + [frac_cell(k, n) for k, n in zip(ks, ns)])
    w(table(hdr, rows))
    w('## 7. Where the network and the path part (per target, eps = 1e-8 / 1e-6 / 1e-4)\n')
    w('Targets on which the network succeeds and the path does not (`net only`) and the reverse (`path only`), summed over the noise levels, per (K, R) and separation.\n')
    rows = []
    for (K, iR) in sorted(set((K, iR) for (K, iR, iN) in group_keys(agg))):
        for key, name in (('only_net', 'net only'), ('only_path', 'path only')):
            for ie in (2, 1, 0):
                cellv = []
                for iD in range(len(DELTAS)):
                    cellv.append(str(sum(agg[(K, iR, iN, iD)][key][ie] for iN in range(len(NOISES)) if (K, iR, iN, iD) in agg)))
                rows.append([K, '%g' % RS[iR], '%s %g' % (name, EPS[ie])] + cellv)
    w(table(['K', 'R', 'who'] + ['%g' % d for d in DELTAS], rows))
    w('## 8. Sensitivity (registered): the 32 configurations on the same test targets\n')
    w('The primary configuration is the frozen one; the rest are not verdicts.  `P1` ... `P4` as in section 3, recomputed under each configuration; `P1 count` = the (K, R) whose '
      'transition is at or below the computed threshold; the medians are over the (K, R) of the design (`>7` = never reached).\n')
    w(table(['#', 'stop', 'threshold', 'merge radius', 'P1', 'P2', 'P3', 'P4', 'P1 count', 'median net 1e-8', 'median path'], sensitivity_rows(test)))
    miss = sorted(set((c[0], c[2], c[3]) for c in CELLS) - set(group_keys(agg)))
    if miss:
        w('Cells missing or incomplete (no claim is made for them): %s.\n' % ', '.join('(K=%d, R=%g, noise=%g)' % (K_, RS[iR], NOISES[iN]) for (K_, iR, iN) in miss))
    return '\n'.join(lines) + '\n'


def sensitivity_rows(test):
    """For every configuration: the verdicts of P1-P4 and the median transition (over (K, R)) of the network at eps = 1e-8 and of the path, noiseless."""
    rows = []
    for ic, cfg in enumerate(CFGS):
        agg = aggregate(test, 'test', cfg)
        P = predictions(agg)
        n8 = [val(r['net8']) for r in P['rows']]
        pa = [val(r['path']) for r in P['rows']]
        rows.append([ic, cfg[0], cfg[1], '%g' % cfg[2], 'held' if P['P1'] else 'failed', 'held' if P['P2'] else 'failed', 'held' if P['P3'] else 'failed',
                     'held' if P['P4'] else 'failed', sum(r['p1'] for r in P['rows']), fmt_t(float(np.median(n8)), flag_of(float(np.median(n8)))),
                     fmt_t(float(np.median(pa)), flag_of(float(np.median(pa))))])
    return rows


def load_raws():
    d = {}
    t, dv = load_raw(RAW), load_raw(DEV_RAW)
    if t:
        d['test'] = t
    if dv:
        d['dev'] = dv
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dev', action='store_true')
    ap.add_argument('--test', action='store_true')
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--shard', default=None, help='i/n: only the cells with index mod n == i')
    ap.add_argument('--assemble', default=None, help='a directory of ckpt_<split>*.jsonl files: build the raw of --dev/--test from them')
    A = ap.parse_args()
    only = [int(x) for x in os.environ['SS_CELLS'].split(',')] if os.environ.get('SS_CELLS') else None
    shard = tuple(int(x) for x in A.shard.split('/')) if A.shard else None
    tag = '_s%d' % shard[0] if shard else ''
    for split, flag in (('dev', A.dev), ('test', A.test)):
        if not flag:
            continue
        if A.assemble:
            import glob
            save_raw(DEV_RAW if split == 'dev' else RAW, assemble(split, sorted(glob.glob(os.path.join(A.assemble, 'ckpt_%s*.jsonl' % split)))))
        else:
            assert split == 'dev' or FROZEN is not None, 'the configuration is not frozen'
            run_split(split, A.resume, only, shard, A.workers, tag)
            if not shard:
                save_raw(DEV_RAW if split == 'dev' else RAW, assemble(split, [os.path.join(WORK, 'ckpt_%s.jsonl' % split)]))
    if A.from_raw or A.assemble or (not shard and (A.dev or A.test)):
        os.makedirs(OUT, exist_ok=True)
        with open(MD, 'w') as fh:
            fh.write(report(load_raws()))


if __name__ == '__main__':
    main()
