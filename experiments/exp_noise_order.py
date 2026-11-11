"""
Choosing the number of kernels under noise: the sliding inverse scale space against the usual alternatives (WP3), one-dimensional design
=======================================================================================================================================

The setting of Paper 17 / the parametrisation study: 61 samples of [-1, 1], a dictionary of 21 Gaussian kernels of width 0.12 at spacing 0.1,
a target of K = 1 .. 6 Gaussian features (random centres at least 1.5 widths apart, signs random, amplitudes in [0.3, 1] x the largest), white
Gaussian noise of standard deviation s = 0.01, 0.03, 0.1, 0.3 x the largest amplitude (the noise level is KNOWN to the rules that need it). For
each trial every method returns a number of kernels K-hat and the kernels themselves:

  iss    the sliding inverse scale space (`baselines.iss_path`, the limit form reproducing `sliding_iss_nd`), K-hat by the discrepancy rule, by the
         significance of the entrant (4.5), by the largest relative gap of the entry times (no noise level), by BIC (2 parameters per kernel)
  omp    matching pursuit with sliding (same candidates and block fits; the entrant is the instantaneous correlation): discrepancy, significance, BIC
  step   forward selection on the grid by least squares (centres stay on the grid), BIC
  lasso  the l1 path on the grid, penalty by BIC (K-hat = the support size; 'merged': atoms within 1.5 grid steps count as one spike)
  pen    penalised training: the 21-unit network trained from a standard start with an l1 penalty on the coefficients and trainable centres,
         penalty by BIC (K-hat = the number of nonzero coefficients)
  prune  random initialisation: the same network without penalty from a random start of O(1), units below 10 % of the largest |w| removed
  sfw    sliding Frank–Wolfe for the BLASSO (`baselines.sfw`) on a ladder of penalties lam = f x max|<y phi>| (f = 0.5 ... 0.003, geometric, fixed a priori),
         lam by BIC (2 parameters per kernel) or by the discrepancy principle (the largest lam whose residual is within the noise; needs s)
  cpgd   conic particle gradient descent (`baselines.cpgd`, 3000 steps) on the same ladder, the particles clustered into spikes; lam by BIC or discrepancy
         (both BLASSO solvers drop spikes below 10 % of the largest weight, `baselines.threshold_units`: the weak spikes of a small penalty fit noise)

Scored: K-hat = K (exact order), and among those the position and relative amplitude errors of the kernels (matched within half a width), and the
seconds. Nothing is tuned per method beyond the one rule it comes with (the penalty ladders are fixed a priori; BIC counts two parameters for a
kernel whose centre is fitted: pen, sfw, cpgd, omp, iss). Usage:
    python experiments/exp_noise_order.py --run [--workers 4] [--targets 6] [--draws 3]      (writes results/chaos/noise_order_raw.npz)
    python experiments/exp_noise_order.py --from-raw                                          (regenerates results/chaos/noise_order.md)
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
from src.core import baselines as B             # noqa: E402
from src.core import rich_transit_nd as RN      # noqa: E402

OUT = os.environ.get('NO_OUT', os.path.join(REPO, 'results', 'chaos'))
RAW, MD = os.path.join(OUT, 'noise_order_raw.npz'), os.path.join(OUT, 'noise_order.md')

SIG, NS, ND = 0.12, 61, 21
NOISE = (0.01, 0.03, 0.1, 0.3)
KS = (1, 2, 3, 4, 5, 6)
MIN_SEP = 1.5 * SIG
TOL = 0.5 * SIG
KAPPA = 4.5
MAX_UNITS = 10
PEN_LAMS = tuple(np.geomspace(3e-4, 3e-1, 7))
LAM_REL = tuple(np.geomspace(0.5, 0.003, 8))                 # the penalties of sfw and cpgd as fractions of max|<y phi>| over the start grid
COLS = ('iss|disc', 'iss|sig', 'iss|gap', 'iss|bic', 'omp|disc', 'omp|sig', 'omp|bic', 'sfw|disc', 'sfw|bic', 'cpgd|disc', 'cpgd|bic', 'step|bic', 'lasso|bic', 'lasso|merged', 'pen|bic', 'prune')
LABEL = {'iss|disc': 'sliding ISS, discrepancy (s)', 'iss|sig': 'sliding ISS, significance (s)', 'iss|gap': 'sliding ISS, largest gap (no s)',
         'iss|bic': 'sliding ISS, BIC', 'omp|disc': 'matching pursuit + sliding, discrepancy (s)', 'omp|sig': 'matching pursuit + sliding, significance (s)',
         'omp|bic': 'matching pursuit + sliding, BIC', 'sfw|disc': 'sliding Frank-Wolfe, discrepancy (s)', 'sfw|bic': 'sliding Frank-Wolfe, BIC',
         'cpgd|disc': 'conic particle gradient descent, discrepancy (s)', 'cpgd|bic': 'conic particle gradient descent, BIC', 'step|bic': 'forward selection on the grid, BIC', 'lasso|bic': 'l1 on the grid, BIC (support)',
         'lasso|merged': 'l1 on the grid, BIC (merged atoms)', 'pen|bic': 'penalised network, BIC', 'prune': 'random start + prune'}


def make_target(rng, K):
    """K centres in [-0.85, 0.85] at least MIN_SEP apart (rejection), random signs, magnitudes uniform in [0.3, 1], the largest scaled to 1."""
    for _ in range(10000):
        c = np.sort(rng.uniform(-0.85, 0.85, K))
        if K == 1 or np.diff(c).min() >= MIN_SEP:
            break
    a = rng.uniform(0.3, 1.0, K) * rng.choice([-1.0, 1.0], K)
    a[int(np.argmax(np.abs(a)))] = np.sign(a[int(np.argmax(np.abs(a)))])
    return c, a


def setup():
    xs = np.linspace(-1, 1, NS)
    return RN.BumpsND(xs[:, None], SIG), np.linspace(-1, 1, ND)


def _score(Khat, w, c, ct, at):
    if len(w) == 0:                                      # a method that returns no unit (e.g. every coefficient pruned): order 0, nothing matched
        return [Khat, 0, np.nan, np.nan]
    n, pos, amp = B.match_units(np.asarray(c).reshape(len(w), -1), w, ct[:, None], at, TOL)
    return [Khat, n, pos, amp]


def trial(args):
    """One trial: (K, noise level index, target seed, draw), see `_trial`; an exception inside it is recorded as a failed trial (NaN rows, the message printed once) instead of ending the run."""
    try:
        return _trial(args)
    except Exception as e:                               # noqa: BLE001 - recorded, counted in the raw and reported
        print('FAILED trial %r: %s: %s' % (args, type(e).__name__, e), flush=True)
        K_true, ni, seed, draw = args
        return np.full((len(COLS), 4), np.nan), np.full(8, np.nan), np.array([K_true, ni, seed, draw, np.nan, np.nan])


def _trial(args):
    """One trial: (K, noise level index, target seed, draw) -> array (len(COLS), 4) = K-hat, matched, position rms, amplitude rms, and seconds."""
    K_true, ni, seed, draw = args
    K, grid = setup()
    rng = np.random.default_rng(seed)
    ct, at = make_target(rng, K_true)
    f = K.synth(ct[:, None], at)
    s = NOISE[ni] * np.abs(at).max()
    y = f + s * np.random.default_rng(seed * 1000 + draw + 1).standard_normal(K.N)
    N = K.N
    out = np.full((len(COLS), 4), np.nan)
    secs = {}

    def unit(blk):
        return blk['w'], blk['c']

    t0 = time.time()
    R = B.iss_path(K, y, grid, max_entries=MAX_UNITS, s=s)
    secs['iss'] = time.time() - t0
    for col, k in (('iss|disc', B.select_discrepancy(R['blocks'], N, s)), ('iss|sig', B.select_significance(R['sig'], KAPPA)),
                   ('iss|gap', B.select_gap(R['taus'])), ('iss|bic', B.select_bic(R['blocks'], N, 2))):
        w, c = unit(R['blocks'][k])
        out[COLS.index(col)] = _score(len(w), w, c, ct, at)
    t0 = time.time()
    O = B.omp_slide(K, y, grid, max_units=MAX_UNITS, s=s)
    secs['omp'] = time.time() - t0
    for col, k in (('omp|disc', B.select_discrepancy(O['blocks'], N, s)), ('omp|sig', B.select_significance(O['sig'], KAPPA)), ('omp|bic', B.select_bic(O['blocks'], N, 2))):
        w, c = unit(O['blocks'][k])
        out[COLS.index(col)] = _score(len(w), w, c, ct, at)
    G = float(np.abs(K.corr(grid[:, None], y)[0]).max())
    t0 = time.time()
    SF = [B.threshold_units(K, y, B.sfw(K, y, grid, f * G, max_units=MAX_UNITS), f * G) for f in LAM_REL]
    secs['sfw'] = time.time() - t0
    for col, k in (('sfw|disc', B.select_discrepancy(SF, N, s)), ('sfw|bic', int(np.argmin([B.bic_value(q['rss'], N, 2 * len(q['w'])) for q in SF])))):
        out[COLS.index(col)] = _score(len(SF[k]['w']), SF[k]['w'], SF[k]['c'], ct, at)
    t0 = time.time()
    CP = [B.cpgd(K, y, grid, f * G) for f in LAM_REL]
    secs['cpgd'] = time.time() - t0
    for col, k in (('cpgd|disc', B.select_discrepancy(CP, N, s)), ('cpgd|bic', int(np.argmin([B.bic_value(q['rss'], N, 2 * q['n_kept']) for q in CP])))):
        out[COLS.index(col)] = _score(CP[k]['n_kept'], CP[k]['w'], CP[k]['c'], ct, at)
    t0 = time.time()
    S = B.stepwise_bic(K, y, grid, max_units=MAX_UNITS)
    secs['step'] = time.time() - t0
    w, c = unit(S['blocks'][B.select_bic(S['blocks'], N, 1)])
    out[COLS.index('step|bic')] = _score(len(w), w, c, ct, at)
    t0 = time.time()
    P = B.lasso_grid(K, y, grid, n_lam=40)
    secs['lasso'] = time.time() - t0
    bics = [B.bic_value(p['rss'], N, p['nnz']) for p in P]
    p = P[int(np.argmin(bics))]
    nz = p['w'] != 0
    out[COLS.index('lasso|bic')] = _score(int(nz.sum()), p['w'][nz], grid[nz], ct, at)
    out[COLS.index('lasso|merged')] = _score(B.merge_adjacent(p['w'], grid, 1.5 * (grid[1] - grid[0])), p['w'][nz], grid[nz], ct, at)
    t0 = time.time()
    net = B.penalised_net(K, y, grid, PEN_LAMS, steps=2500, seed=seed)
    secs['pen'] = time.time() - t0
    b = int(np.argmin([B.bic_value(q['rss'], N, 2 * q['nnz']) for q in net]))
    q = net[b]
    nz = np.abs(q['w']) > 1e-6 * max(np.abs(q['w']).max(), 1e-12)
    out[COLS.index('pen|bic')] = _score(int(nz.sum()), q['w'][nz], q['c'][nz], ct, at)
    t0 = time.time()
    Z = B.prune_net(K, y, grid, steps=5000, seed=seed)
    secs['prune'] = time.time() - t0
    out[COLS.index('prune')] = _score(Z['n_kept'], Z['w'][Z['kept']], Z['c'][Z['kept']], ct, at)
    return out, np.array([secs[m] for m in ('iss', 'omp', 'sfw', 'cpgd', 'step', 'lasso', 'pen', 'prune')]), np.array([K_true, ni, seed, draw, np.diff(ct).min() if K_true > 1 else np.nan,
                                                                                                 np.abs(at).min() / np.abs(at).max()])


def run(workers, n_targets, n_draws, out_path=RAW):
    from multiprocessing import Pool
    jobs = [(K, ni, 1000 * K + t, d) for K in KS for t in range(n_targets) for ni in range(len(NOISE)) for d in range(n_draws)]
    t0 = time.time()
    with Pool(workers) as pool:
        res = pool.map(trial, jobs, chunksize=2)
    np.savez_compressed(out_path, res=np.stack([r[0] for r in res]), secs=np.stack([r[1] for r in res]), meta=np.stack([r[2] for r in res]),
                        cols=np.array(COLS), noise=np.array(NOISE), n_targets=n_targets, n_draws=n_draws)
    print('%d trials in %.0f s' % (len(jobs), time.time() - t0), flush=True)


def flow_trial(args):
    """The network of Paper 17 itself, trained by its gradient flow from the small start eps = 1e-8 on the noisy data, stopped midway between the
    entry times of the discrepancy block and the next entrant (the plateau of that block): the number of kernels it has formed (|w| above a twentieth of
    the smallest coefficient of the path's block) against the path's block. Returns (K, K-hat of the path by the discrepancy rule, K-hat of the trained
    network, rms difference of the formed coefficients' fit from the path block's, seconds)."""
    from experiments import exp_parametrization as EP
    from src.core import factorized_rbf as FR
    K_true, ni, seed, draw = args
    K, grid = setup()
    rng = np.random.default_rng(seed)
    ct, at = make_target(rng, K_true)
    f = K.synth(ct[:, None], at)
    s = NOISE[ni] * np.abs(at).max()
    y = f + s * np.random.default_rng(seed * 1000 + draw + 1).standard_normal(K.N)
    R = B.iss_path(K, y, grid, max_entries=MAX_UNITS, s=s)
    kd = B.select_discrepancy(R['blocks'], K.N, s)
    taus = R['taus']
    if kd == 0:
        return np.array([K_true, 0, 0, 0.0, 0.0])
    tau_end = 0.5 * (taus[kd - 1] + taus[kd]) if kd < len(taus) else 1.5 * taus[kd - 1]
    t0 = time.time()
    res = FR.train_flow(FR.PaperLaw(EP.law()), K, y, FR.pts(grid), 1e-8, tau_end, EP.LR_C, max_secs=300)
    w, c = res['W'][-1], res['C'][-1]
    wb = R['blocks'][kd]['w']
    formed = np.abs(w) > 0.05 * np.abs(wb).min()
    fit = K.synth(c[formed].reshape(-1, 1), w[formed]) if formed.any() else np.zeros(K.N)
    blk = R['blocks'][kd]
    ref = K.synth(blk['c'], blk['w'])
    return np.array([K_true, kd, int(formed.sum()), float(np.sqrt(np.mean((fit - ref) ** 2)) / max(np.abs(f).max(), 1e-12)), time.time() - t0])


def run_flow(workers, n, out_path=RAW):
    from multiprocessing import Pool
    jobs = [(K, 1 + (i % 3), 1000 * K + i, 0) for i, K in enumerate([2, 3, 4, 3, 5, 2, 3, 4] * ((n + 7) // 8))][:n]
    with Pool(workers) as pool:
        res = np.stack(pool.map(flow_trial, jobs, chunksize=1))
    old = dict(np.load(out_path, allow_pickle=False)) if os.path.exists(out_path) else {}
    old['flow'] = res
    np.savez_compressed(out_path, **old)
    print('flow: %d trials, trained network = path block in %d' % (len(res), int((res[:, 1] == res[:, 2]).sum())), flush=True)


def wilson(k, n, z=1.96):
    if n == 0:
        return float('nan'), float('nan')
    p = k / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return mid - half, mid + half


def target_ids(Kt, seeds):
    """One id per trial for its target (K, seed): trials of the same target share it."""
    return np.asarray(Kt).astype(int) * 10 ** 6 + np.asarray(seeds).astype(int)


def target_cluster_diffs(exact, targets, base, others, B=2000, seed=0):
    """Differences of the exact-order rate, column `base` minus each column of `others`, over the trials (rows of the boolean matrix `exact`), with 95 % percentile intervals of a bootstrap
    over the TARGETS (`targets`: one id per trial; the trials of a target are resampled together, so that the trials of one target are not counted as independent). Returns [(difference, lo, hi)]."""
    ids = np.unique(targets)
    groups = [np.flatnonzero(targets == t) for t in ids]
    draw = np.random.default_rng(seed).integers(0, len(ids), size=(B, len(ids)))
    out = []
    for jo in others:
        d = np.array([exact[g, base].mean() - exact[g, jo].mean() for g in groups])
        boot = d[draw].mean(axis=1)
        out.append((float(d.mean()), float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))))
    return out


def order_strata(pairs):
    """For each target (K, seed): 0 if the noiseless entry sequences (the first K entries) of the sliding inverse scale space (the accumulated dual enters first) and of matching pursuit
    with sliding (the instantaneous correlation) are identical, 1 if they are the same kernels in another order, 2 if the kernels differ. Deterministic: the targets are rebuilt from
    their seeds."""
    K, grid = setup()
    out = {}
    for Kn, seed in sorted(set((int(a), int(b)) for a, b in pairs)):
        ct, at = make_target(np.random.default_rng(seed), Kn)
        f = K.synth(ct[:, None], at)
        a = [int(j) for j in B.iss_path(K, f, grid, max_entries=Kn)['order'][:Kn]]
        b = [int(j) for j in B.omp_slide(K, f, grid, max_units=Kn)['order'][:Kn]]
        out[(Kn, seed)] = 0 if a == b else (1 if sorted(a) == sorted(b) else 2)
    return out


STRATA = {0: 'identical entry sequences', 1: 'same kernels, another order', 2: 'different kernels'}


def report(path=RAW):
    z = np.load(path, allow_pickle=False)
    res, meta, secs = z['res'], z['meta'], z['secs']
    bad = np.isnan(res[:, 0, 0])                                         # a failed trial (an exception, see `trial`) is all NaN: it is left out and counted
    res, meta, secs = res[~bad], meta[~bad], secs[~bad]
    cols = [str(c) for c in z['cols']]
    nl = len(NOISE)
    Kt, ni = meta[:, 0].astype(int), meta[:, 1].astype(int)
    exact = res[:, :, 0] == Kt[:, None]
    L = ['# Choosing the number of kernels under noise: the sliding inverse scale space against the usual alternatives (one-dimensional design)', '',
         'Generated by `experiments/exp_noise_order.py` from `noise_order_raw.npz`; every number below is printed by the script. Setting and methods: the '
         'script\'s docstring. %d trials: K = 1..6 features, %d targets per K, %d noise draws, four noise levels (s as a fraction of the largest amplitude: %s). '
         '"Exact" is K-hat = K; intervals are 95 %% Wilson intervals of the proportion over the trials of that column. Failed trials (an exception): %d.\n' % (
             len(res), int(z['n_targets']), int(z['n_draws']), ', '.join('%g' % v for v in NOISE), int(bad.sum()))]
    L.append('## Exact order, by noise level\n')
    hdr = ['method and rule'] + ['s = %g (n = %d)' % (v, int((ni == i).sum())) for i, v in enumerate(NOISE)] + ['all']
    rows = []
    for j, c in enumerate(cols):
        cells = []
        for i in range(nl):
            m = ni == i
            k = int(exact[m, j].sum())
            lo, hi = wilson(k, int(m.sum()))
            cells.append('%.0f %% [%.0f, %.0f]' % (100 * k / m.sum(), 100 * lo, 100 * hi))
        k = int(exact[:, j].sum())
        lo, hi = wilson(k, len(res))
        cells.append('%.0f %% [%.0f, %.0f]' % (100 * k / len(res), 100 * lo, 100 * hi))
        rows.append([LABEL.get(c, c)] + cells)
    L.append('| ' + ' | '.join(hdr) + ' |\n|' + '|'.join(['---'] * len(hdr)) + '|\n' + '\n'.join('| ' + ' | '.join(r) + ' |' for r in rows))
    tid = target_ids(Kt, meta[:, 2])
    jb = cols.index('iss|bic')
    oth = [j for j in range(len(cols)) if j != jb]
    cd = target_cluster_diffs(exact, tid, jb, oth)
    L.append('\n## Exact order: differences clustered by target\n')
    L.append('The intervals above count the trials of a target as independent. Difference of the exact-order rate of the sliding ISS with the BIC rule minus each other column (positive: the ISS is better), '
             'with 95 %% percentile intervals of a bootstrap over the %d targets (2000 resamples, seed 0; the trials of a target are resampled together).\n' % len(np.unique(tid)))
    L.append('| method and rule | difference (points) |\n|---|---|\n' + '\n'.join('| %s | %+.1f [%+.1f, %+.1f] |' % (LABEL.get(cols[j], cols[j]), 100 * d, 100 * lo, 100 * hi) for j, (d, lo, hi) in zip(oth, cd)))
    L.append('\n## Exact order, by the true number of kernels (all noise levels)\n')
    hdr = ['method and rule'] + ['K = %d' % k for k in KS]
    rows = [[LABEL.get(c, c)] + ['%.0f %%' % (100 * exact[Kt == k, j].mean()) for k in KS] for j, c in enumerate(cols)]
    L.append('| ' + ' | '.join(hdr) + ' |\n|' + '|'.join(['---'] * len(hdr)) + '|\n' + '\n'.join('| ' + ' | '.join(r) + ' |' for r in rows))
    L.append('\n## Bias of the order and accuracy of the kernels\n')
    hdr = ['method and rule', 'mean K-hat - K', 'over-estimates %', 'under-estimates %', 'position rms (exact order, units of the width)', 'amplitude rms (exact order, relative)']
    rows = []
    for j, c in enumerate(cols):
        d = res[:, j, 0] - Kt
        e = exact[:, j] & (res[:, j, 1] == Kt)
        pos = np.sqrt(np.nanmean(res[e, j, 2] ** 2)) / SIG if e.any() else float('nan')
        amp = np.sqrt(np.nanmean(res[e, j, 3] ** 2)) if e.any() else float('nan')
        rows.append([LABEL.get(c, c), '%+.2f' % d.mean(), '%.0f' % (100 * (d > 0).mean()), '%.0f' % (100 * (d < 0).mean()), '%.3f' % pos, '%.3f' % amp])
    L.append('| ' + ' | '.join(hdr) + ' |\n|' + '|'.join(['---'] * len(hdr)) + '|\n' + '\n'.join('| ' + ' | '.join(r) + ' |' for r in rows))
    L.append('\nPosition and amplitude errors are over the trials where K-hat = K and every kernel is matched within half a width.\n')
    strata = order_strata(zip(Kt, meta[:, 2]))
    st = np.array([strata[(int(k), int(sd))] for k, sd in zip(Kt, meta[:, 2])])
    L.append('## Where the two orders differ (accumulated dual against instantaneous correlation)\n')
    L.append('Strata by the noiseless entry sequences of the two methods on the same target (first K entries): %s. Exact-order rates over the trials of a stratum, 95 %% Wilson intervals.\n' % (
        '; '.join('%d = %s' % (k, v) for k, v in STRATA.items())))
    pairs_ = (('discrepancy', 'iss|disc', 'omp|disc'), ('significance', 'iss|sig', 'omp|sig'), ('BIC', 'iss|bic', 'omp|bic'))
    hdr = ['stratum', 'targets', 'trials'] + [h for nm, _, _ in pairs_ for h in ('%s: sliding ISS' % nm, '%s: matching pursuit + sliding' % nm)]
    rows = []
    for g in sorted(STRATA):
        m = st == g
        if not m.any():
            continue
        row = [STRATA[g], str(len({(int(k), int(sd)) for k, sd in zip(Kt[m], meta[m, 2])})), str(int(m.sum()))]
        for _, ci_, co_ in pairs_:
            for c_ in (ci_, co_):
                k = int(exact[m, cols.index(c_)].sum())
                lo, hi = wilson(k, int(m.sum()))
                row.append('%.0f %% [%.0f, %.0f]' % (100 * k / m.sum(), 100 * lo, 100 * hi))
        rows.append(row)
    L.append('| ' + ' | '.join(hdr) + ' |\n|' + '|'.join(['---'] * len(hdr)) + '|\n' + '\n'.join('| ' + ' | '.join(r) + ' |' for r in rows))
    L.append('\n## Cost (seconds per trial, one core, a shared machine)\n')
    names = ('sliding ISS (limit path, 10 entries)', 'matching pursuit + sliding', 'sliding Frank-Wolfe (8 penalties)', 'conic particle gradient descent (8 penalties, 3000 steps)',
             'forward selection', 'l1 path (40 penalties)', 'penalised network (7 penalties, 2500 steps)', 'random start + prune (5000 steps)')
    L.append('| ' + ' | '.join(['method'] + ['median s']) + ' |\n|---|---|\n' + '\n'.join('| %s | %.3f |' % (n, float(np.median(secs[:, i]))) for i, n in enumerate(names)))
    if 'flow' in z.files:
        fl = z['flow']
        L.append('\n## The network itself (gradient flow from eps = 1e-8) against the path, noisy data\n')
        L.append('%d noisy targets (noise levels 0.03, 0.1, 0.3 of the largest amplitude): the network of Paper 17 trained by its gradient flow to the plateau of the '
                 'path\'s discrepancy block has formed the path block\'s kernels in %d of %d trials (the number of formed kernels differs in %d); the relative rms '
                 'difference of the fits is %.2g at the median and %.2g at the maximum (the data fit rms of the same targets: %s of the largest amplitude).\n' % (
                     len(fl), int((fl[:, 1] == fl[:, 2]).sum()), len(fl), int((fl[:, 1] != fl[:, 2]).sum()), float(np.median(fl[:, 3])), float(fl[:, 3].max()),
                     ', '.join('%g' % v for v in NOISE[1:])))
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--flow', type=int, default=0, help='also train the real network on this many noisy targets (slow)')
    ap.add_argument('--workers', type=int, default=int(os.environ.get('NO_WORKERS', '3')))
    ap.add_argument('--targets', type=int, default=6)
    ap.add_argument('--draws', type=int, default=3)
    A = ap.parse_args()
    if A.run:
        run(A.workers, A.targets, A.draws)
    if A.flow:
        run_flow(A.workers, A.flow)
    if A.from_raw:
        with open(MD, 'w') as f:
            f.write(report())
        print('wrote', MD)


if __name__ == '__main__':
    main()
