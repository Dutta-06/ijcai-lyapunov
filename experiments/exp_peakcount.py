"""
Does the small-start sliding Gaussian path count pulses in a real 1-D signal better than tuned classical counters?  (workstream J spike)
=======================================================================================================================================

SPIKE with a kill rule (plan and GO/KILL rule: results/ml/peakcount_plan.md, committed before any test record was scored; dataset scouting:
results/ml/peakcount_scout.md).  Data: the 2000 training records of the CPSC 2019 challenge (single-lead ECG, 10 s at 500 Hz, expert-annotated R
peaks; Zenodo record 23179391).  The target is the NUMBER of R peaks in a record (positions are scored as a secondary quantity).

Front end (the same for every method): band-pass 5-18 Hz, derivative, square (or absolute value), 120 ms moving mean, decimation to 50 Hz, division
by the 95th percentile.  Each beat is then a bump of about 0.1 s; the components of the model are Gaussian bumps of ONE width fitted to that envelope.

Methods (each family gets the same number of development configurations, N_CFG = 42, and its configuration is frozen before the test run):
    fp       scipy find_peaks on the envelope: prominence (fraction of the 95th percentile) x minimum distance
    fp_rel   find_peaks with a prominence relative to the record's own tall peaks (75th percentile of the prominences of a low-threshold pass) x distance
    mf       matched filter (Gaussian template) then find_peaks: template width x prominence x distance
    step     forward selection of grid Gaussians by least squares (centres on the grid, `baselines.stepwise_bic`), stopped by BIC with a penalty multiplier or by a residual-drop rule
    omp      matching pursuit with sliding (`baselines.omp_slide`: the greedy cousin), same stopping rules
    iss      the sliding inverse scale space path (`baselines.iss_path`, the limit form of the factorised small-start network), same stopping rules
The path families (step, omp, iss) run at three kernel widths; their stopping rules are the penalised BIC (7 multipliers of 2 ln N per unit) and the
residual-drop rule (7 thresholds), 3 x 14 = 42 configurations.  Counting is the number of units of the selected block of the path.

    python experiments/exp_peakcount.py --fetch                          (download + checksum of the raw archive into PC_DATA)
    python experiments/exp_peakcount.py --stage fe|dev|freeze|test       (needs PC_DATA; writes peakcount_*.{json,npz,md} into PC_OUT, default results/ml)
    python experiments/exp_peakcount.py --from-raw                       (numpy + scipy only: regenerates peakcount.md byte-identically from peakcount_raw.npz)
"""
import argparse
import hashlib
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, REPO)

OUT = os.environ.get('PC_OUT', os.path.join(REPO, 'results', 'ml'))
DATA = os.environ.get('PC_DATA', '')

ZENODO_URL = 'https://zenodo.org/api/records/23179391/files/CPSC2019-train.zip/content'
ZIP_SHA256 = '1a838bcab24d5b93de5dfebd0c9f5b43b9a621415a6472f2056f6b79c7d5c7fc'

# ------------------------------------------------------------------------------------------------------------------------------------
# registries (everything the plan fixes)
# ------------------------------------------------------------------------------------------------------------------------------------
N_RECORDS, FS, DUR, DEC = 2000, 500, 10.0, 10
EFS = FS // DEC                                   # envelope rate (Hz)
SPLIT_SEED, N_DEV, N_TEST = 20261110, 150, 500    # records are permuted with this seed; the first N_DEV are development, the next N_TEST are test, the rest unused
BAND, WIN_S, P_NORM = (5.0, 18.0), 0.12, 95.0
FRONT_ENDS = ('sq', 'abs')
N_CFG = 42
FP_PROM = (0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35)
FP_DIST = (0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
REL_Q = (0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60)
REL_LOW, REL_PCT = 0.02, 75.0
MF_SIG = (0.03, 0.05, 0.07)
MF_DIST = (0.20, 0.30)
WIDTHS = (0.03, 0.04, 0.05)                       # kernel widths (s) of the path families (wider kernels make the sign-constrained polish of the ISS path 20-100 times slower: computing limit)
BIC_C = (1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)    # multipliers of the BIC penalty 2 ln N per unit
DROP_D = (0.001, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08)   # residual-drop thresholds, fraction of the empty-fit residual
GRID_STEP = 0.1                                   # start grid of candidate centres (s)
MAXM = 44                                         # cap on the number of entries of a path
CAP_PROM, CAP_DIST, CAP_MULT, CAP_ADD = 0.05, 0.14, 1.4, 4   # label-free cap M = min(MAXM, int(CAP_MULT * P) + CAP_ADD), P = peaks of a loose find_peaks pass
TOL_S = 0.10                                      # position tolerance (s) for the secondary F1
PATH_FAMS = ('step', 'omp', 'iss')
CLASSICAL = ('fp', 'fp_rel', 'mf', 'step')        # the "best tuned classical" is the one of these with the lowest development error
FAMS = ('fp', 'fp_rel', 'mf', 'step', 'omp', 'iss')
LABEL = {'fp': 'find_peaks (prominence x distance)', 'fp_rel': 'find_peaks, prominence relative to the record', 'mf': 'matched filter + find_peaks',
         'step': 'forward selection on the grid + stopping rule', 'omp': 'matching pursuit with sliding + stopping rule (greedy control)',
         'iss': 'sliding inverse scale space path + stopping rule (ours)'}
N_BOOT, BOOT_SEED = 10000, 20261111
NI_MARGIN = 0.10                                  # ours may be worse than the greedy sliding control by at most this many beats per record (upper bound of the interval)


# ------------------------------------------------------------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------------------------------------------------------------
def fetch(dest):
    """Download the Zenodo archive into dest/CPSC2019-train.zip, check SHA-256, unpack to dest/data and dest/ref."""
    import urllib.request
    import zipfile
    os.makedirs(dest, exist_ok=True)
    z = os.path.join(dest, 'CPSC2019-train.zip')
    if not os.path.exists(z):
        urllib.request.urlretrieve(ZENODO_URL, z)
    h = hashlib.sha256(open(z, 'rb').read()).hexdigest()
    if h != ZIP_SHA256:
        raise RuntimeError('checksum mismatch: %s' % h)
    with zipfile.ZipFile(z) as f:
        f.extractall(dest)
    return h


def load_record(i, root=None):
    """(ecg (5000,), R-peak sample indices) of record i (1-based)."""
    import scipy.io as sio
    root = root or DATA
    x = sio.loadmat(os.path.join(root, 'data', 'data_%05d.mat' % i))['ecg'].ravel().astype(float)
    r = sio.loadmat(os.path.join(root, 'ref', 'R_%05d.mat' % i))['R_peak'].ravel().astype(int)
    return x, r


def split_ids():
    """(dev ids, test ids): the 1-based record ids permuted by numpy default_rng(SPLIT_SEED); first N_DEV development, next N_TEST test."""
    p = np.random.default_rng(SPLIT_SEED).permutation(np.arange(1, N_RECORDS + 1))
    return np.sort(p[:N_DEV]), np.sort(p[N_DEV:N_DEV + N_TEST])


def envelope(x, fe='sq'):
    """The shared front end: band-pass, derivative, square or absolute value, moving mean, decimation to EFS Hz, division by the 95th percentile."""
    from scipy import signal
    b, a = signal.butter(2, BAND, btype='band', fs=FS)
    d = np.gradient(signal.filtfilt(b, a, x))
    e = d ** 2 if fe == 'sq' else np.abs(d)
    win = int(round(WIN_S * FS))
    e = np.convolve(e, np.ones(win) / win, mode='same')
    e = e[:len(e) // DEC * DEC].reshape(-1, DEC).mean(1)
    s = float(np.percentile(e, P_NORM))
    return e / (s if s > 0 else 1.0)


def path_cap(e):
    """Label-free cap on the number of entries of a path: from a loose find_peaks pass on the envelope."""
    from scipy.signal import find_peaks
    P = len(find_peaks(e, prominence=CAP_PROM, distance=max(1, int(CAP_DIST * EFS)))[0])
    return int(min(MAXM, int(CAP_MULT * P) + CAP_ADD))


# ------------------------------------------------------------------------------------------------------------------------------------
# classical counters (positions in seconds)
# ------------------------------------------------------------------------------------------------------------------------------------
def fp_positions(e, prom, dist):
    from scipy.signal import find_peaks
    pk, _ = find_peaks(e, prominence=prom, distance=max(1, int(round(dist * EFS))))
    return pk / EFS


def fp_rel_positions(e, q, dist):
    """Peaks of a low-threshold pass whose prominence is at least q times the REL_PCT percentile of the pass's prominences."""
    from scipy.signal import find_peaks
    pk, pr = find_peaks(e, prominence=REL_LOW, distance=max(1, int(round(dist * EFS))))
    if len(pk) == 0:
        return pk / EFS
    ref = np.percentile(pr['prominences'], REL_PCT)
    return pk[pr['prominences'] >= q * ref] / EFS


def mf_positions(e, sig, prom, dist):
    """Matched filter: correlate with a Gaussian template of width sig (s), normalise by the 95th percentile of the output, find_peaks."""
    from scipy.signal import find_peaks
    h = int(np.ceil(3 * sig * EFS))
    t = np.arange(-h, h + 1) / EFS
    g = np.exp(-t ** 2 / (2 * sig ** 2))
    y = np.convolve(e, g / g.sum(), mode='same')
    s = float(np.percentile(y, P_NORM))
    y = y / (s if s > 0 else 1.0)
    pk, _ = find_peaks(y, prominence=prom, distance=max(1, int(round(dist * EFS))))
    return pk / EFS


# ------------------------------------------------------------------------------------------------------------------------------------
# paths
# ------------------------------------------------------------------------------------------------------------------------------------
def run_path(e, fam, width, M):
    """One path (fam in PATH_FAMS) of at most M entries on the envelope e with Gaussians of the given width (s) on the grid GRID_STEP.
    Returns rss (MAXM+1,), cnt (MAXM+1,), cen (MAXM+1, MAXM): NaN-padded beyond the blocks produced; block k holds cnt[k] units with sorted centres cen[k, :cnt[k]]."""
    from src.core import baselines as B
    from src.core import rich_transit_nd as RN
    N = len(e)
    t = np.arange(N) / EFS
    K = RN.BumpsND(t[:, None], width)
    grid = np.arange(0.0, N / EFS, GRID_STEP)
    if fam == 'iss':
        R = B.iss_path(K, e, grid, max_entries=M)
    elif fam == 'omp':
        R = B.omp_slide(K, e, grid, max_units=M)
    else:
        R = B.stepwise_bic(K, e, grid, max_units=M)
    rss = np.full(MAXM + 1, np.nan)
    cnt = np.full(MAXM + 1, np.nan)
    cen = np.full((MAXM + 1, MAXM), np.nan)
    for k, b in enumerate(R['blocks'][:MAXM + 1]):
        c = np.sort(np.asarray(b['c'], float).ravel())
        rss[k], cnt[k] = b['rss'], len(c)
        cen[k, :len(c)] = c
    return rss, cnt, cen


def select_block(rss, cnt, N, rule, p):
    """Index of the selected block of a path.  rule 'bic': argmin N ln(rss/N) + p * 2 * cnt * ln N (p the multiplier of the usual penalty);
    rule 'drop': the first block k with rss[k] - rss[k+1] < p * rss[0] (the last block if none)."""
    ok = np.flatnonzero(~np.isnan(rss))
    if rule == 'bic':
        v = N * np.log(np.maximum(rss[ok], 1e-300) / N) + p * 2.0 * cnt[ok] * np.log(N)
        return int(ok[int(np.argmin(v))])
    r0 = rss[ok[0]]
    for a, b in zip(ok[:-1], ok[1:]):
        if rss[a] - rss[b] < p * r0:
            return int(a)
    return int(ok[-1])


def path_positions(P, i, rule, p, N):
    """Positions (s) of the selected block of record i of the stored path arrays P = (rss, cnt, cen)."""
    rss, cnt, cen = P[0][i], P[1][i], P[2][i]
    k = select_block(rss, cnt, N, rule, p)
    return cen[k, :int(cnt[k])]


def path_cfgs():
    """The 42 configurations of a path family: (width index, rule, parameter)."""
    return [(wi, 'bic', c) for wi in range(len(WIDTHS)) for c in BIC_C] + [(wi, 'drop', d) for wi in range(len(WIDTHS)) for d in DROP_D]


def classical_cfgs(fam):
    if fam == 'fp':
        return [(p, d) for p in FP_PROM for d in FP_DIST]
    if fam == 'fp_rel':
        return [(q, d) for q in REL_Q for d in FP_DIST]
    if fam == 'mf':
        return [(s, p, d) for s in MF_SIG for p in FP_PROM for d in MF_DIST]
    raise KeyError(fam)


def classical_positions(fam, e, cfg):
    if fam == 'fp':
        return fp_positions(e, *cfg)
    if fam == 'fp_rel':
        return fp_rel_positions(e, *cfg)
    return mf_positions(e, *cfg)


# ------------------------------------------------------------------------------------------------------------------------------------
# scoring and statistics
# ------------------------------------------------------------------------------------------------------------------------------------
def match_f1(est, ref, tol=TOL_S):
    """One-to-one matching (Hungarian on |difference|, pairs farther than tol are not matches): (tp, fp, fn)."""
    from scipy.optimize import linear_sum_assignment
    est, ref = np.asarray(est, float), np.asarray(ref, float)
    if len(est) == 0 or len(ref) == 0:
        return 0, len(est), len(ref)
    d = np.abs(est[:, None] - ref[None, :])
    r, c = linear_sum_assignment(np.where(d <= tol, d, 1e6))
    tp = int(np.sum(d[r, c] <= tol))
    return tp, len(est) - tp, len(ref) - tp


def micro_f1(tp, fp, fn):
    return 2.0 * tp / max(2.0 * tp + fp + fn, 1)


def boot_ci(d, n_boot=N_BOOT, seed=BOOT_SEED):
    """Percentile 95% interval of the mean of d over records (paired bootstrap): (mean, lo, hi)."""
    d = np.asarray(d, float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=(n_boot, len(d)))
    m = d[idx].mean(1)
    return float(d.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


class Store:
    """The records of one split with everything the methods need: ids, envelopes, reference positions (s), and (after run) the stored paths."""

    def __init__(self, ids, fe):
        self.ids, self.fe = np.asarray(ids), fe
        self.E, self.ref = [], []
        for i in self.ids:
            x, r = load_record(int(i))
            self.E.append(envelope(x, fe))
            self.ref.append(r / FS)
        self.E = np.array(self.E)
        self.K = np.array([len(r) for r in self.ref])
        self.paths = {}                                 # (fam, width index) -> (rss, cnt, cen), arrays over records

    def positions(self, fam, cfg, shift=0.0):
        """Estimated positions (s, shifted by -shift) of every record under the configuration cfg of family fam."""
        out = []
        for j in range(len(self.ids)):
            if fam in PATH_FAMS:
                wi, rule, p = cfg
                pos = path_positions(self.paths[(fam, wi)], j, rule, p, self.E.shape[1])
            else:
                pos = classical_positions(fam, self.E[j], cfg)
            out.append(np.asarray(pos, float) - shift)
        return out

    def counts(self, fam, cfg):
        if fam in PATH_FAMS:
            wi, rule, p = cfg
            rss, cnt, _ = self.paths[(fam, wi)]
            return np.array([cnt[j][select_block(rss[j], cnt[j], self.E.shape[1], rule, p)] for j in range(len(self.ids))])
        return np.array([len(classical_positions(fam, e, cfg)) for e in self.E])

    def f1(self, fam, cfg, shift):
        tp = fp = fn = 0
        for est, ref in zip(self.positions(fam, cfg, shift), self.ref):
            a, b, c = match_f1(est, ref)
            tp, fp, fn = tp + a, fp + b, fn + c
        return micro_f1(tp, fp, fn)


def all_cfgs(fam):
    return path_cfgs() if fam in PATH_FAMS else classical_cfgs(fam)


def cfg_name(fam, cfg):
    if fam in PATH_FAMS:
        wi, rule, p = cfg
        return 'width %.2f s, %s %g' % (WIDTHS[wi], rule, p)
    return ', '.join('%g' % v for v in cfg)


def tune(store, fam):
    """Development tuning of a family: (best cfg, mean abs error per cfg (list), signed mean per cfg).  Lowest mean absolute count error; ties go to the first configuration."""
    cfgs = all_cfgs(fam)
    mae, bias = [], []
    for cfg in cfgs:
        err = store.counts(fam, cfg) - store.K
        mae.append(float(np.mean(np.abs(err))))
        bias.append(float(np.mean(err)))
    k = int(np.argmin(mae))
    return cfgs[k], mae, bias


def position_shift(store, cfg_fp):
    """The shared offset (s) between the envelope peaks and the annotated R peaks: the median of (estimate - reference) over matched pairs of the tuned find_peaks counter."""
    d = []
    for est, ref in zip(store.positions('fp', cfg_fp), store.ref):
        if len(est) and len(ref):
            j = np.abs(est[:, None] - ref[None, :]).argmin(1)
            dd = est - ref[j]
            d += list(dd[np.abs(dd) <= 0.2])
    return float(np.median(d)) if d else 0.0


# ------------------------------------------------------------------------------------------------------------------------------------
# stages
# ------------------------------------------------------------------------------------------------------------------------------------
def stage_fe():
    """Choose the front end on development records by the best find_peaks error (the classical counter picks the front end: favours the baselines)."""
    dev, _ = split_ids()
    res = {}
    for fe in FRONT_ENDS:
        st = Store(dev, fe)
        cfg, mae, _ = tune(st, 'fp')
        res[fe] = dict(cfg=list(cfg), mae=float(min(mae)))
    best = min(FRONT_ENDS, key=lambda f: res[f]['mae'])
    out = dict(front_end=best, dev_fp=res, dev_ids=[int(i) for i in dev])
    json.dump(out, open(os.path.join(OUT, 'peakcount_fe.json'), 'w'), indent=1, sort_keys=True)
    print(json.dumps(out)[:400])
    return out


def compute_paths(store, which, log=None):
    """Run the path families listed in which = [(fam, width index), ...] on every record of the store (resumable by the caller)."""
    caps = np.array([path_cap(e) for e in store.E])
    for fam, wi in which:
        n = len(store.ids)
        rss, cnt, cen = np.full((n, MAXM + 1), np.nan), np.full((n, MAXM + 1), np.nan), np.full((n, MAXM + 1, MAXM), np.nan)
        t0 = time.time()
        for j in range(n):
            rss[j], cnt[j], cen[j] = run_path(store.E[j], fam, WIDTHS[wi], int(caps[j]))
            if log and (j + 1) % 25 == 0:
                print('  %s w%d: %d/%d records, %.0f s' % (fam, wi, j + 1, n, time.time() - t0), file=log, flush=True)
        store.paths[(fam, wi)] = (rss, cnt, cen)
    return caps


def save_raw(store, caps, path, extra=None):
    d = dict(ids=store.ids, fe=np.array(store.fe), E=store.E, K=store.K, caps=caps,
             ref_flat=np.concatenate(store.ref) if len(store.ref) else np.zeros(0), ref_len=np.array([len(r) for r in store.ref]))
    for (fam, wi), (rss, cnt, cen) in store.paths.items():
        d['%s_%d_rss' % (fam, wi)], d['%s_%d_cnt' % (fam, wi)], d['%s_%d_cen' % (fam, wi)] = rss, cnt, cen
    if extra:
        d.update(extra)
    np.savez_compressed(path, **d)


def load_raw(path):
    z = np.load(path, allow_pickle=False)
    st = Store.__new__(Store)
    st.ids, st.fe, st.E, st.K = z['ids'], str(z['fe']), z['E'], z['K']
    flat, lens = z['ref_flat'], z['ref_len']
    off = np.r_[0, np.cumsum(lens)]
    st.ref = [flat[off[j]:off[j + 1]] for j in range(len(lens))]
    st.paths = {}
    for fam in PATH_FAMS:
        for wi in range(len(WIDTHS)):
            if '%s_%d_rss' % (fam, wi) in z.files:
                st.paths[(fam, wi)] = (z['%s_%d_rss' % (fam, wi)], z['%s_%d_cnt' % (fam, wi)], z['%s_%d_cen' % (fam, wi)])
    return st


def stage_dev():
    fe = json.load(open(os.path.join(OUT, 'peakcount_fe.json')))['front_end']
    dev, _ = split_ids()
    st = Store(dev, fe)
    caps = compute_paths(st, [(f, wi) for f in PATH_FAMS for wi in range(len(WIDTHS))], log=sys.stdout)
    save_raw(st, caps, os.path.join(OUT, 'peakcount_dev_raw.npz'))


def stage_freeze():
    """Tune every family on the development raw, write peakcount_freeze.json and the development report."""
    st = load_raw(os.path.join(OUT, 'peakcount_dev_raw.npz'))
    frz = dict(front_end=st.fe, families={})
    rows = []
    for fam in FAMS:
        cfg, mae, bias = tune(st, fam)
        frz['families'][fam] = dict(cfg=list(cfg), dev_mae=float(min(mae)), dev_bias=float(bias[int(np.argmin(mae))]), n_cfg=len(mae))
        rows.append((fam, cfg, min(mae), bias[int(np.argmin(mae))], len(mae)))
    best = min(CLASSICAL, key=lambda f: frz['families'][f]['dev_mae'])
    frz['best_classical'] = best
    frz['shift'] = position_shift(st, tuple(frz['families']['fp']['cfg']))
    json.dump(frz, open(os.path.join(OUT, 'peakcount_freeze.json'), 'w'), indent=1, sort_keys=True)
    lines = ['# Peak counting spike: development summary (150 development records)', '',
             'Front end: `%s`.  Mean absolute count error of the frozen configuration of each family (lowest over its %d configurations):' % (st.fe, N_CFG), '',
             '| family | frozen configuration | dev mean abs error | dev mean signed error |', '|---|---|---|---|']
    for fam, cfg, m, b, n in rows:
        lines.append('| %s | %s | %.3f | %+.3f |' % (LABEL[fam], cfg_name(fam, tuple(cfg)), m, b))
    lines += ['', 'Best tuned classical counter (lowest development error among %s): **%s**.  Shared position offset: %.3f s.' % (', '.join(CLASSICAL), best, frz['shift'])]
    open(os.path.join(OUT, 'peakcount_dev.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def stage_test():
    frz = json.load(open(os.path.join(OUT, 'peakcount_freeze.json')))
    _, test = split_ids()
    st = Store(test, frz['front_end'])
    which = sorted({(f, int(frz['families'][f]['cfg'][0])) for f in PATH_FAMS})
    caps = compute_paths(st, which, log=sys.stdout)
    save_raw(st, caps, os.path.join(OUT, 'peakcount_raw.npz'))


# ------------------------------------------------------------------------------------------------------------------------------------
# verdict and report
# ------------------------------------------------------------------------------------------------------------------------------------
def frozen_cfg(frz, fam):
    c = frz['families'][fam]['cfg']
    if fam in PATH_FAMS:
        return (int(c[0]), str(c[1]), float(c[2]))
    return tuple(float(v) for v in c)


def evaluate(st, frz):
    """Per-record absolute errors, signed errors and F1 counts of every family at its frozen configuration."""
    res = {}
    for fam in FAMS:
        cfg = frozen_cfg(frz, fam)
        if fam in PATH_FAMS and (fam, cfg[0]) not in st.paths:
            continue
        cnt = st.counts(fam, cfg)
        pos = st.positions(fam, cfg, frz['shift'])
        tp = fp = fn = 0
        for est, ref in zip(pos, st.ref):
            a, b, c = match_f1(est, ref)
            tp, fp, fn = tp + a, fp + b, fn + c
        res[fam] = dict(cfg=cfg, err=(cnt - st.K).astype(float), f1=micro_f1(tp, fp, fn), tp=tp, fp=fp, fn=fn)
    return res


def verdict(res, frz):
    """GO if ours has a lower mean absolute count error than the best tuned classical counter (the one frozen on development) with the paired-bootstrap interval of the
    difference excluding 0, AND ours is not worse than the greedy sliding matching pursuit by more than NI_MARGIN (upper end of the interval of ours minus omp < NI_MARGIN).
    Returns dict(go, win, ni, base, win_ci, ni_ci)."""
    base = frz['best_classical']
    a_o, a_b, a_g = np.abs(res['iss']['err']), np.abs(res[base]['err']), np.abs(res['omp']['err'])
    win_ci = boot_ci(a_b - a_o)
    ni_ci = boot_ci(a_o - a_g)
    win = win_ci[1] > 0
    ni = ni_ci[2] < NI_MARGIN
    return dict(go=bool(win and ni), win=bool(win), ni=bool(ni), base=base, win_ci=win_ci, ni_ci=ni_ci)


def render_report(st, frz):
    res = evaluate(st, frz)
    v = verdict(res, frz)
    K = st.K
    n = len(K)
    L = []
    L.append('# Peak counting spike: test report (%d held-out CPSC 2019 records)' % n)
    L.append('')
    L.append('Plan and rule: `results/ml/peakcount_plan.md` (committed before any test record was scored). Data: CPSC 2019 training records, expert R-peak annotations; '
             'target = number of annotated R peaks per 10 s record (mean %.2f, range %d to %d). Front end `%s`. Every configuration below was frozen on the development records '
             '(`results/ml/peakcount_freeze.json`).' % (K.mean(), K.min(), K.max(), st.fe))
    L.append('')
    L.append('## Verdict: **%s**' % ('GO' if v['go'] else 'KILL'))
    L.append('')
    L.append('* Win condition: ours (`iss`) against the best tuned classical counter frozen on development (`%s`): mean absolute count error difference (baseline minus ours) '
             '%+.3f, 95%% bootstrap interval [%+.3f, %+.3f] -> %s.' % (v['base'], v['win_ci'][0], v['win_ci'][1], v['win_ci'][2], 'met' if v['win'] else 'not met'))
    L.append('* Non-inferiority to the greedy sliding matching pursuit (`omp`): ours minus omp %+.3f, 95%% interval [%+.3f, %+.3f], margin %.2f beats per record -> %s.'
             % (v['ni_ci'][0], v['ni_ci'][1], v['ni_ci'][2], NI_MARGIN, 'met' if v['ni'] else 'not met'))
    L.append('')
    L.append('## Count error of every family at its frozen configuration (test records)')
    L.append('')
    L.append('| family | frozen configuration | mean abs error | exact count | mean signed error | secondary: position F1 (tol %.2f s) |' % TOL_S)
    L.append('|---|---|---|---|---|---|')
    for fam in FAMS:
        if fam not in res:
            continue
        r = res[fam]
        L.append('| %s | %s | %.3f | %.3f | %+.3f | %.3f |' % (LABEL[fam], cfg_name(fam, r['cfg']), np.mean(np.abs(r['err'])), np.mean(r['err'] == 0), np.mean(r['err']), r['f1']))
    L.append('')
    L.append('## Paired differences of mean absolute error against ours (baseline minus ours; positive = ours better), 95% bootstrap intervals')
    L.append('')
    L.append('| baseline | difference | 95% interval |')
    L.append('|---|---|---|')
    ao = np.abs(res['iss']['err'])
    for fam in FAMS:
        if fam in res and fam != 'iss':
            m, lo, hi = boot_ci(np.abs(res[fam]['err']) - ao)
            L.append('| %s | %+.3f | [%+.3f, %+.3f] |' % (LABEL[fam], m, lo, hi))
    L.append('')
    L.append('## Where the errors are (mean absolute error by true count)')
    L.append('')
    bins = [(0, 10), (11, 14), (15, 19), (20, 99)]
    L.append('| true beats | records | ' + ' | '.join(f for f in FAMS if f in res) + ' |')
    L.append('|---|---|' + '---|' * len([f for f in FAMS if f in res]))
    for lo, hi in bins:
        m = (K >= lo) & (K <= hi)
        if m.sum():
            L.append('| %d-%d | %d | ' % (lo, hi, m.sum()) + ' | '.join('%.2f' % np.mean(np.abs(res[f]['err'][m])) for f in FAMS if f in res) + ' |')
    L.append('')
    L.append('Notes. Errors mix the model-order error with the front-end error, identically for every method. Positions are scored after removing the shared offset '
             '%.3f s found on development; the F1 is secondary and not part of the verdict. The "ours" family is the limit form of the factorised small-start network '
             '(`baselines.iss_path`), not a trained torch network.' % frz['shift'])
    return '\n'.join(L) + '\n', v


def from_raw(raw=None, md=None):
    raw = raw or os.path.join(OUT, 'peakcount_raw.npz')
    md = md or os.path.join(OUT, 'peakcount.md')
    frz = json.load(open(os.path.join(os.path.dirname(raw), 'peakcount_freeze.json')))
    st = load_raw(raw)
    text, v = render_report(st, frz)
    open(md, 'w').write(text)
    print('verdict', 'GO' if v['go'] else 'KILL')
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fetch', action='store_true')
    ap.add_argument('--stage', choices=['fe', 'dev', 'freeze', 'test'])
    ap.add_argument('--from-raw', nargs='?', const='')
    ap.add_argument('--md')
    a = ap.parse_args()
    if a.fetch:
        print(fetch(DATA or os.path.join(REPO, 'data', 'cpsc2019')))
    elif a.stage:
        {'fe': stage_fe, 'dev': stage_dev, 'freeze': stage_freeze, 'test': stage_test}[a.stage]()
    elif a.from_raw is not None:
        from_raw(a.from_raw or None, a.md)


if __name__ == '__main__':
    main()
