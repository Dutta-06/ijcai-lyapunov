"""
Picasso baseline LADDER on the simulated 2D microtubule benchmarks
==================================================================

Picasso is the baseline of the localization experiments. A baseline is only as fair as the options it is allowed to use, so it is built as a
ladder, every rung tuned on DEVELOPMENT frames only and scored on held-out TEST frames with the harness of `experiments/smlm_common.py`:

  R0   Picasso's defaults (minimum net gradient 5000, box 7, raw frames, no padding, least-squares fit gausslq)
  R1   the minimum net gradient tuned on the development frames (the old "tuned" rows), by efficiency E (R1E) and by Jaccard J (R1J)
  R2   tuned over ALL the documented options at once (grid below, declared before any run): identification box, fit box, Picasso's Gaussian
       filter, median padding of the frame, fit method (gausslq, gaussmle), then, for every structural setting, a coordinate ascent of the
       post-filters (net gradient, minimum photons, width window); selection by development E (R2E) and by development J (R2J)
  +c   every rung also with a constant position-offset correction, estimated on the development frames of that very configuration (the mean signed
       error of its matched true positives) and subtracted from every detection (development and test) before scoring

Sets (development frames / test frames): MT0.N1.LD (0-999 / 1000-1499), MT0.N1.HD (0-99 / 100-499), MT0.N2.LD and MT0.N2.HD (same protocol; the
N2 sets have about 10 times fewer photons, so the net-gradient and photon grids are scaled). Test frames only SCORE the rungs, nothing is ever chosen
on them. Reported per rung and set, in the primary setting (all ground-truth emitters) and in the challenge setting (emitters within 450 nm of the
border and the dimmest quarter of the whole set's activations dropped after pairing): Jaccard, precision, recall, RMSE, the RMSE about the mean error
(bias removed), the RMSE on the true positives COMMON to every rung (free of the survivor bias of a stricter rung), efficiency, all with 95 %
intervals of a bootstrap in blocks of 10 frames, and paired differences between rungs; the ablation ladder on the development frames; recall by distance
to the border, by photons and by |z|.

The detections of every rung (and of its +c version) on the development and the test frames are saved, packed as offsets + flat arrays, in
`results/smlm/picasso_ladder_raw.npz` so that another method can be scored against exactly these detections (see the README entry in that file).

Picasso runs in an interpreter that can import it: set SMLM_PICASSO_PYTHON or pass --picasso-python.
Output: results/smlm/picasso_baseline.md and results/smlm/picasso_ladder_raw.npz (both regenerate byte for byte; timings go to the console only).
Usage: python experiments/exp_smlm_picasso.py [--sets MT0.N1.LD ...] [--boot 1000] [--jobs 2] [--cache-dir DIR]
"""
import argparse
import ast
import collections
import concurrent.futures
import hashlib
import itertools
import json
import multiprocessing
import os
import platform
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
import smlm_common as C

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SETS = C.SETS
PROTOCOL = {s: dict(dev=(0, 1000), test=(1000, 1500)) if s.endswith('LD') else dict(dev=(0, 100), test=(100, 500)) for s in SETS}
CHALLENGE = dict(border_nm=450.0, keep=0.75)
BOOT = 1000                          # bootstrap resamples (seed 0), in blocks of C.DEFAULT_BLOCK = 10 frames
BLOCK = C.DEFAULT_BLOCK
DEFAULT_NG = 5000.0                  # Picasso's default minimum net gradient (R0), for every set

# ==================================================================================================================================
# THE GRIDS OF THE LADDER. DECLARED HERE BEFORE ANY RUN OF THE SEARCH; NOT CHANGED AFTER SEEING A RESULT.
# ==================================================================================================================================
# R1: the minimum net gradient only, 20 values (the grid of the first tuning; the N2 grid is the one of experiments/exp_smlm_iss.py, PIC_GRID['N2'])
GRID_R1 = {'N1': [1000, 1500, 2000, 2500, 3000, 3250, 3500, 3750, 4000, 4500, 5000, 5500, 6000, 7000, 8000, 9000, 10000, 12000, 15000, 20000],
           'N2': [100, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 8000]}
# R2, structural settings (every combination is run once at the smallest net gradient of its set; 3 x 3 x 3 x 2 x 2 = 108 settings)
STRUCT_GRID = collections.OrderedDict([('box', (5, 7, 9)),                      # identification box (px)
                                       ('fit_box', (5, 7, 9)),                  # fit box (px)
                                       ('gauss_sigma', (None, 1.0, 1.5)),       # Picasso's "Gaussian Filter Sigma" (px); None = off
                                       ('pad', (0, 4)),                         # median padding (px); 0 = off
                                       ('fit', ('gausslq', 'gaussmle'))])       # least squares / Poisson maximum likelihood, both point-sampled Gaussians
DEFAULT_STRUCT = collections.OrderedDict([('box', 7), ('fit_box', 7), ('gauss_sigma', None), ('pad', 0), ('fit', 'gausslq')])
# R2, post-filters (coordinate ascent for every structural setting; the order below is the order of the line searches; at most MAX_ROUNDS rounds)
MAX_ROUNDS = 3
FILTER_GRIDS = {
    'N1': collections.OrderedDict([
        ('min_net_gradient', [300, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750, 4000, 4500, 5000, 5500, 6000, 7000, 8000, 9000,
                              10000, 12000, 15000, 20000]),
        ('min_photons', [None, 100, 200, 300, 400, 600, 800, 1000, 1500]),
        ('min_sigma_nm', [None, 70, 90, 100, 110]),
        ('max_sigma_nm', [None, 200, 250, 300, 350, 400, 500])]),
    'N2': collections.OrderedDict([
        ('min_net_gradient', [30, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000,
                              5000, 6000, 8000]),
        ('min_photons', [None, 10, 20, 30, 40, 60, 80, 100, 150]),
        ('min_sigma_nm', [None, 70, 90, 100, 110]),
        ('max_sigma_nm', [None, 200, 250, 300, 350, 400, 500])])}
RUNGS = ('R0', 'R1E', 'R1J', 'R2E', 'R2J')
VARIANTS = tuple(r + s for r in RUNGS for s in ('', '+c'))          # the ten sets of detections scored: every rung without and with the offset correction
PHOTON_BINS = {'N1': ((0, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 7000), (7000, 1e9)),
               'N2': ((0, 50), (50, 100), (100, 200), (200, 400), (400, 700), (700, 1e9))}
EDGE_BINS = ((0, 300), (300, 450), (450, 700), (700, 1000), (1000, 7000))
Z_BINS = ((0, 150), (150, 300), (300, 450), (450, 1000))


def noise_of(name):
    return 'N2' if '.N2.' in name else 'N1'


def struct_list():
    """Every structural setting of the declared grid, in a fixed order (the order of STRUCT_GRID, last key fastest)."""
    return [collections.OrderedDict(zip(STRUCT_GRID, v)) for v in itertools.product(*STRUCT_GRID.values())]


def skey(struct):
    return '|'.join('%s=%s' % (k, struct[k]) for k in STRUCT_GRID)


def ng_low(noise):
    return float(min(FILTER_GRIDS[noise]['min_net_gradient']))


def iss_pic_grid():
    """PIC_GRID of experiments/exp_smlm_iss.py read from its source text (not imported: that module is heavy and owned by another stage), or None."""
    try:
        tree = ast.parse(open(os.path.join(os.path.dirname(__file__), 'exp_smlm_iss.py')).read())
        for node in tree.body:
            if isinstance(node, ast.Assign) and any(getattr(t, 'id', None) == 'PIC_GRID' for t in node.targets):
                return ast.literal_eval(node.value)
    except Exception:
        pass
    return None


def describe(params):
    """One-line description of a Picasso configuration (a dict of run_picasso options)."""
    p = C._merge_params(params)
    box, fb = int(p['box']), int(p['fit_box'] or p['box'])
    parts = ['ident %d / fit %d px' % (box, fb) if box != fb else 'box %d px' % box,
             'Gaussian filter %g px' % p['gauss_sigma'] if p['gauss_sigma'] else 'raw frames',
             'median pad %d px' % p['pad'] if p['pad'] else 'no padding', p['fit'], 'net gradient > %g' % p['min_net_gradient']]
    if p['min_photons'] is not None:
        parts.append('photons > %g' % p['min_photons'])
    lo, hi = p['min_sigma_nm'], p['max_sigma_nm']
    if lo is not None or hi is not None:
        parts.append('%s < width < %s nm' % ('%g' % lo if lo is not None else '0', '%g' % hi if hi is not None else 'inf'))
    return ', '.join(parts)


def ci_str(v, ci, scale=1.0, nd=1):
    return '%.*f [%.*f, %.*f]' % (nd, scale * v, nd, scale * ci[0], nd, scale * ci[1])


def md_table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    return '\n'.join(out + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows])


# ==================================================================================================================================
# the development search (R1 and R2)
# ==================================================================================================================================
def _compact(res):
    return dict(filters=res['filters'], value=res['value'], metrics=res['metrics'], rounds=res['rounds'], converged=res['converged'], first=res['first'])


def search_structures(session, gt_dev, dev_idx, noise, log=print):
    """R2: every structural setting of STRUCT_GRID, one low-threshold run each, then the coordinate ascent of the post-filters twice (objective
    efficiency and objective Jaccard; one shared cache, so no filter combination is scored twice). Development frames only."""
    grids = FILTER_GRIDS[noise]
    out = collections.OrderedDict()
    structs = struct_list()
    t_all = time.time()
    for i, st in enumerate(structs):
        t0 = time.time()
        cand = session.run(dict(st, min_net_gradient=ng_low(noise)))
        cache = {}
        rE = C.tune_postfilters(cand, gt_dev, dev_idx, grids, objective='efficiency', max_rounds=MAX_ROUNDS, cache=cache)
        rJ = C.tune_postfilters(cand, gt_dev, dev_idx, grids, objective='jaccard', max_rounds=MAX_ROUNDS, cache=cache)
        n_cand = int(sum(len(d) for d in cand))
        n_bad = int(sum(int((~np.isfinite(d[:, :2])).any(axis=1).sum()) for d in cand if len(d)))
        out[skey(st)] = dict(struct=dict(st), n_cand=n_cand, n_nonfinite=n_bad, n_eval=len(cache), E=_compact(rE), J=_compact(rJ))
        log('  [%d/%d] %s: %d candidates, %d configurations scored, dev E %.1f (rounds %d), dev J %.1f (rounds %d)  [%.0f s, total %.0f s]' % (
            i + 1, len(structs), skey(st), n_cand, len(cache), rE['metrics']['efficiency'], rE['rounds'], 100 * rJ['metrics']['jaccard'], rJ['rounds'],
            time.time() - t0, time.time() - t_all))
    return out


def choose(results):
    """The structural settings with the highest development efficiency (E) and Jaccard (J); ties: the earlier one in the grid order."""
    bestE = bestJ = None
    for k, r in results.items():
        if bestE is None or r['E']['metrics']['efficiency'] > results[bestE]['E']['metrics']['efficiency']:
            bestE = k
        if bestJ is None or r['J']['metrics']['jaccard'] > results[bestJ]['J']['metrics']['jaccard']:
            bestJ = k
    return bestE, bestJ


def params_of(struct, filters):
    p = dict(struct)
    p.update({k: v for k, v in filters.items()})
    return p


def ablation_rows(results, key, obj, extra):
    """Marginal effect of every option on the development score at the chosen configuration ``key`` of objective ``obj`` ('E' or 'J'): the best
    development configuration (post-filters re-tuned by the same ascent) with ONE option put back to its default. ``extra``: development scores of
    the chosen structure with the photon / width post-filter switched off, from the candidates of that structure.
    Returns (label, metrics or None, note) rows; None with a note when the option is already at its default (or the structure is not in the grid)."""
    ch = results[key]
    st = dict(ch['struct'])
    rows = [('chosen configuration', ch[obj]['metrics'], '')]
    names = collections.OrderedDict([('box', 'identification box -> 7 px'), ('fit_box', 'fit box -> the identification box'),
                                     ('gauss_sigma', 'Gaussian filter -> off'), ('pad', 'median padding -> 0'), ('fit', 'fit -> gausslq')])
    for k, label in names.items():
        s2 = dict(st)
        s2[k] = s2['box'] if k == 'fit_box' else DEFAULT_STRUCT[k]
        if s2[k] == st[k]:
            rows.append((label, None, 'already the default'))
        elif skey(s2) not in results:
            rows.append((label, None, 'not in the grid'))
        else:
            rows.append((label, results[skey(s2)][obj]['metrics'], ''))
    rows += [('photon filter -> off (the other post-filters kept, not re-tuned)', extra['photons'], ''), ('width window -> off (the other post-filters kept, not re-tuned)', extra['width'], '')]
    rows.append(('every structural option -> default (post-filters re-tuned)', results[skey(DEFAULT_STRUCT)][obj]['metrics'], ''))
    return rows


def forward_rows(results, obj):
    """Each option alone on top of R1 (default structure with the net gradient tuned, no other filter): the structural setting with ONE option
    changed from the default, net gradient tuned by the objective (the first line search of the ascent). Development scores."""
    rows = [('default structure, net gradient tuned (R1 on the finer grid)', results[skey(DEFAULT_STRUCT)][obj]['first']['metrics'])]
    for k, vals in STRUCT_GRID.items():
        for v in vals:
            if v == DEFAULT_STRUCT[k]:
                continue
            s2 = dict(DEFAULT_STRUCT)
            s2[k] = v
            rows.append(('%s = %s' % (k, 'off' if v is None else v), results[skey(s2)][obj]['first']['metrics']))
    return rows


# ==================================================================================================================================
# scoring helpers
# ==================================================================================================================================
def score_both(dets, gt, idx, pmin):
    return C.score(dets, gt, idx), C.score(dets, gt, idx, border_nm=CHALLENGE['border_nm'], photon_min=pmin)


def kept_mask(gt_rows_gt, pmin):
    """Ground-truth rows kept by the challenge setting (record array with x_nm, y_nm, photons)."""
    b = CHALLENGE['border_nm']
    return ((gt_rows_gt['x_nm'] >= b) & (gt_rows_gt['x_nm'] <= C.FIELD_NM - b) & (gt_rows_gt['y_nm'] >= b) & (gt_rows_gt['y_nm'] <= C.FIELD_NM - b)
            & (gt_rows_gt['photons'] >= pmin))


def count_stealing(dets, gt, idx, pmin):
    """Kept ground-truth rows (challenge setting) that stay unmatched although a detection within the tolerance was paired with an EXCLUDED row
    of the frame (the pairing precedes the exclusions): returns (stolen, kept rows)."""
    g_sorted = gt[np.argsort(gt['frame'], kind='stable')]
    lo = np.searchsorted(g_sorted['frame'], idx, side='left')
    hi = np.searchsorted(g_sorted['frame'], idx, side='right')
    stolen = kept_total = 0
    for i in range(len(idx)):
        gf = g_sorted[lo[i]:hi[i]]
        if len(gf) == 0:
            continue
        kp = kept_mask(gf, pmin)
        kept_total += int(kp.sum())
        if kp.all() or len(dets[i]) == 0:
            continue
        gxy = np.stack([gf['x_nm'], gf['y_nm']], axis=1)
        dxy = C._xy(dets[i])
        m = C.match(dxy, gxy)
        paired_gt = np.zeros(len(gf), bool)
        partner = {}
        for di, gj in m.pairs:
            paired_gt[gj] = True
            partner[int(di)] = int(gj)
        with np.errstate(invalid='ignore'):
            d = np.sqrt(((dxy[:, None, :] - gxy[None, :, :]) ** 2).sum(axis=2))
        for j in np.flatnonzero(kp & ~paired_gt):
            near = np.flatnonzero(d[:, j] <= C.TOL_NM)
            if any((int(di) in partner) and not kp[partner[int(di)]] for di in near):
                stolen += 1
    return stolen, kept_total


def common_columns(me, common, n_frames):
    """Per-frame columns (sse, tp) of the squared errors of the rows flagged ``common`` (one rung's ``matched_errors`` record)."""
    sel = common & me['matched']
    sse = np.bincount(me['pos'][sel], weights=me['dx'][sel] ** 2 + me['dy'][sel] ** 2, minlength=n_frames)
    tp = np.bincount(me['pos'][sel], minlength=n_frames).astype(float)
    return dict(sse=sse, tp=tp)


def rmse_stat(s):
    return np.sqrt(s['sse'] / s['tp'])


def argmax_stability(per_frame_by_thr, B=1000, seed=0, block=BLOCK):
    """How well defined is the optimal threshold? Block-bootstrap resamples of the development frames (the same resampling as the intervals); in each, the
    threshold that maximises the efficiency / the Jaccard index. per_frame_by_thr: {threshold: per_frame array of ``score``}. Returns {objective: [(threshold,
    share of the resamples), ...] sorted by decreasing share} (ties in a resample go to the lower threshold)."""
    thr = sorted(per_frame_by_thr)
    arr = {k: np.stack([np.asarray(per_frame_by_thr[t][k], dtype=float) for t in thr]) for k in ('tp', 'fp', 'fn', 'sse')}
    F = arr['tp'].shape[1]
    wins = {'efficiency': np.zeros(len(thr)), 'jaccard': np.zeros(len(thr))}
    for _, idx in C._boot_indices(F, B, seed, block):
        for ix in idx:
            tp, fp, fn, sse = (arr[k][:, ix].sum(axis=1) for k in ('tp', 'fp', 'fn', 'sse'))
            with np.errstate(invalid='ignore', divide='ignore'):
                jac = tp / (tp + fp + fn)
                eff = 100.0 - np.sqrt((100.0 - 100.0 * jac) ** 2 + sse / tp)
            for name, v in (('efficiency', eff), ('jaccard', jac)):
                wins[name][int(np.argmax(np.where(np.isnan(v), -np.inf, v)))] += 1
    return {name: sorted(zip(thr, w / w.sum()), key=lambda tv: (-tv[1], tv[0])) for name, w in wins.items()}


def neighbour_fractions(gt, frame_indices, radii=(250.0, 400.0)):
    """Fraction of the ground-truth rows of the frames that have another ground-truth emitter of the same frame within each radius (nm)."""
    g = gt[np.argsort(gt['frame'], kind='stable')]
    lo = np.searchsorted(g['frame'], frame_indices, side='left')
    hi = np.searchsorted(g['frame'], frame_indices, side='right')
    n = 0
    hit = np.zeros(len(radii))
    for i in range(len(frame_indices)):
        gf = g[lo[i]:hi[i]]
        if len(gf) < 1:
            continue
        n += len(gf)
        if len(gf) > 1:
            d = np.hypot(gf['x_nm'][:, None] - gf['x_nm'][None, :], gf['y_nm'][:, None] - gf['y_nm'][None, :])
            np.fill_diagonal(d, np.inf)
            dmin = d.min(axis=1)
            hit += np.array([(dmin <= r).sum() for r in radii])
    return hit / max(n, 1), n


# ==================================================================================================================================
# one set
# ==================================================================================================================================
def search_signature(name):
    blob = json.dumps(dict(name=name, proto=PROTOCOL[name], struct=[(k, list(v)) for k, v in STRUCT_GRID.items()], filt=[(k, v) for k, v in FILTER_GRIDS[noise_of(name)].items()],
                           rounds=MAX_ROUNDS, defaults=C.PICASSO_DEFAULTS, r1=GRID_R1[noise_of(name)]), sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def run_set(name, args, log=print):
    """Everything for one set. Returns (markdown lines, raw arrays dict, summary dict). Deterministic: no timing enters the lines."""
    noise = noise_of(name)
    (d0, d1), (t0, t1) = PROTOCOL[name]['dev'], PROTOCOL[name]['test']
    py = args.picasso_python
    stack, gt_all = C.data_path(name, 'stack'), C.load_gt(C.data_path(name, 'gt'))
    dev = C.load_frames(stack, d0, d1)
    test = C.load_frames(stack, t0, t1)
    dev_idx, test_idx = np.arange(d0, d1), np.arange(t0, t1)
    gt_dev = gt_all[(gt_all['frame'] >= d0) & (gt_all['frame'] < d1)]
    gt_test = gt_all[(gt_all['frame'] >= t0) & (gt_all['frame'] < t1)]
    pmin = C.photon_threshold(gt_all, CHALLENGE['keep'])
    nT = len(test)
    L = []
    P = L.append
    P('## %s\n' % name)
    rl = C.run_lengths(gt_all)
    P('Development frames %d-%d (%d frames, %d ground-truth activations), test frames %d-%d (%d frames, %d activations). Frame indices are 0-based stack indices. '
      'Ground-truth activations of the whole set: %d; photon threshold of the challenge setting (25 %% quantile of the whole set\'s ground truth): %.0f photons. '
      'Persistence: a fluorophore stays on for %.2f consecutive frames on average (%d runs, %.0f %% of them 3 frames or longer), so neighbouring frames are correlated, which the block bootstrap accounts for.\n'
      % (d0, d1 - 1, len(dev), len(gt_dev), t0, t1 - 1, len(test), len(gt_test), len(gt_all), pmin, rl['mean'], rl['runs'], 100 * rl['ge3']))

    # ---------------------------------------------------------------- R1: the net gradient only (the first tuning, 20 values)
    t_ = time.time()
    best_e, curve = C.tune_picasso(dev, gt_all, GRID_R1[noise], frame_indices=dev_idx, objective='efficiency', python=py, verify_frames=60)
    best_j = max(curve, key=lambda c: (c['jaccard'], -c['min_net_gradient']))['min_net_gradient']
    log('%s: R1 tuning %.0f s, best efficiency at %g, best Jaccard at %g' % (name, time.time() - t_, best_e, best_j))
    P('### R1: tuning curve of the minimum net gradient on the development frames (%d values)\n' % len(GRID_R1[noise]))
    P(md_table(['min net gradient', 'Jaccard %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'detections / frame'],
               [[('**%g**' if c['min_net_gradient'] in (best_e, best_j) else '%g') % c['min_net_gradient'], '%.1f' % (100 * c['jaccard']), '%.1f' % (100 * c['precision']),
                 '%.1f' % (100 * c['recall']), '%.1f' % c['rmse_nm'], '%.1f' % c['efficiency'], '%.2f' % (c['n_det'] / len(dev))] for c in curve]))
    P('\nBold: the value maximising the efficiency (%g; R1E) and the value maximising the Jaccard index (%g; R1J). The curve is computed from one run at the smallest grid '
      'value, filtered at each threshold (identification keeps the maxima with net gradient above the threshold and every spot is fitted independently); this was verified '
      'to reproduce direct runs at %g and %g on the first 60 development frames.\n' % (best_e, best_j, best_e, GRID_R1[noise][-1]))
    curve_at = {c['min_net_gradient']: c for c in curve}
    with C.PicassoSession(dev, py) as s1:
        full = s1.run(dict(min_net_gradient=float(min(GRID_R1[noise]))))
    stab = argmax_stability({t: C.score([d[d[:, C.NG_COL] > t] for d in full], gt_dev, dev_idx)['per_frame'] for t in curve_at}, B=args.boot, seed=0)
    e_best = max(c['efficiency'] for c in curve)
    near = [c['min_net_gradient'] for c in curve if c['efficiency'] >= e_best - 1.0]
    P('How well defined are these two operating points? On the development frames %d of the %d thresholds (from %g to %g) have an efficiency within 1 point of the best (%.1f): the efficiency is a plateau and its '
      'argmax is arbitrary within it. In %d block-bootstrap resamples of the development frames (blocks of %d frames) the efficiency-optimal threshold is %s; the Jaccard-optimal threshold is %s. '
      'A falling RMSE with a rising threshold is largely survivor bias (the hard emitters are lost): compare the column "RMSE common TP" below, where the rows are scored on the same emitters.\n'
      % (len(near), len(curve), min(near), max(near), e_best, args.boot, BLOCK, ', '.join('%g in %.0f %%' % (t, 100 * v) for t, v in stab['efficiency'][:3]),
         ', '.join('%g in %.0f %%' % (t, 100 * v) for t, v in stab['jaccard'][:3])))

    # ---------------------------------------------------------------- R2: the search over all options (development frames only)
    cache_file = os.path.join(args.cache_dir, '%s.search.json' % name) if args.cache_dir else None
    sig = search_signature(name)
    results = None
    if cache_file and os.path.exists(cache_file):
        blob = json.load(open(cache_file))
        if blob.get('signature') == sig:
            results = collections.OrderedDict(blob['results'])
            log('%s: R2 search results read from %s' % (name, cache_file))
    if results is None:
        t_ = time.time()
        with C.PicassoSession(dev, py) as sess:
            results = search_structures(sess, gt_dev, dev_idx, noise, log=lambda s: log('%s %s' % (name, s)))
        log('%s: R2 search %.0f s' % (name, time.time() - t_))
        if cache_file:
            os.makedirs(args.cache_dir, exist_ok=True)
            json.dump(dict(signature=sig, results=results), open(cache_file, 'w'))
    keyE, keyJ = choose(results)
    n_struct = len(results)
    n_eval = sum(r['n_eval'] for r in results.values())
    n_cand = sum(r['n_cand'] for r in results.values())
    n_bad = sum(r['n_nonfinite'] for r in results.values())
    n_conv = sum(int(r[o]['converged']) for r in results.values() for o in ('E', 'J'))
    chosen = collections.OrderedDict([('R2E', (keyE, 'E')), ('R2J', (keyJ, 'J'))])
    base_params = {'R0': {}, 'R1E': dict(min_net_gradient=float(best_e)), 'R1J': dict(min_net_gradient=float(best_j))}
    for r, (k, o) in chosen.items():
        base_params[r] = params_of(results[k]['struct'], results[k][o]['filters'])
        base_params[r]['min_net_gradient'] = float(base_params[r]['min_net_gradient'])

    # direct runs of every rung on the development and the test frames, and the self-check that filtering the low-threshold run equals the direct run
    extra, ceiling = {}, {}
    with C.PicassoSession(dev, py) as sd, C.PicassoSession(test, py) as st_:
        dets_dev = {r: sd.run(base_params[r]) for r in RUNGS}
        dets_test = {r: st_.run(base_params[r]) for r in RUNGS}
        for r, (k, o) in chosen.items():
            cand = sd.run(dict(results[k]['struct'], min_net_gradient=ng_low(noise)))
            for a, b in zip(dets_dev[r], C.filter_detections(cand, base_params[r])):
                if a.shape != b.shape or not np.allclose(a, b, equal_nan=True):
                    raise RuntimeError('%s %s: the filtered low-threshold run differs from the direct run' % (name, r))
            f0 = dict(results[k][o]['filters'])

            def dev_score(f, cand=cand):
                q = C.score(C.filter_detections(cand, f), gt_dev, dev_idx)
                return {kk: q[kk] for kk in C._TUNE_KEYS}
            extra[r] = dict(photons=dev_score(dict(f0, min_photons=None)), width=dev_score(dict(f0, min_sigma_nm=None, max_sigma_nm=None)))
            qc = C.score(cand, gt_dev, dev_idx)
            ceiling[r] = (qc['recall'], qc['n_det_all'] / len(dev))
    log('%s: direct runs and checks done' % name)

    P('### R2: search over all options on the development frames\n')
    P('Declared grid (identical for the sets of one noise class; at the top of `experiments/exp_smlm_picasso.py`): identification box %s px, fit box %s px, Gaussian filter sigma %s px, '
      'median padding %s px, fit %s: **%d structural settings**, each run once at net gradient %g and then post-filtered by a coordinate ascent (at most %d rounds over, in this order: '
      'net gradient (%d values), minimum photons (%d), minimum width (%d), maximum width (%d)) once for the efficiency and once for the Jaccard index. '
      'Development frames only.\n' % ('/'.join(map(str, STRUCT_GRID['box'])), '/'.join(map(str, STRUCT_GRID['fit_box'])), '/'.join('off' if v is None else '%g' % v for v in STRUCT_GRID['gauss_sigma']),
                                    '/'.join(map(str, STRUCT_GRID['pad'])), '/'.join(STRUCT_GRID['fit']), n_struct, ng_low(noise), MAX_ROUNDS,
                                    *[len(v) for v in FILTER_GRIDS[noise].values()]))
    P('**Configurations tried on the development frames: %d structural settings run, %d distinct (structure, post-filter) configurations scored** (R1: %d net-gradient values of the '
      'default structure, all contained in them); %d candidate localizations fitted in total%s. The coordinate ascent converged (a round without change) within the %d rounds in %d of %d runs.\n'
      % (n_struct, n_eval, len(GRID_R1[noise]), n_cand, ', of which %d with a non-finite position' % n_bad if n_bad else '', MAX_ROUNDS, n_conv, 2 * n_struct))
    dev_res = {r: C.score(dets_dev[r], gt_dev, dev_idx) for r in RUNGS}
    off = {r: (dev_res[r]['bias_x_nm'], dev_res[r]['bias_y_nm']) for r in RUNGS}
    P('Chosen configurations and their scores ON THE DEVELOPMENT FRAMES (all ground-truth emitters); the offset is the mean signed error (detection - truth) of the matched pairs, '
      'the correction of the "+c" rows:\n')
    P(md_table(['rung', 'configuration', 'dev J %', 'dev precision %', 'dev recall %', 'dev RMSE nm', 'dev E', 'offset (dx, dy) nm', 'dev TP'],
               [[r, describe(base_params[r]), '%.1f' % (100 * dev_res[r]['jaccard']), '%.1f' % (100 * dev_res[r]['precision']), '%.1f' % (100 * dev_res[r]['recall']),
                 '%.1f' % dev_res[r]['rmse_nm'], '%.1f' % dev_res[r]['efficiency'], '(%.1f, %.1f)' % off[r], dev_res[r]['tp']] for r in RUNGS]))
    P('')

    def abl_cells(m, note=''):
        return ['%.1f' % (100 * m['jaccard']), '%.1f' % m['efficiency'], '%.1f' % m['rmse_nm']] if m is not None else ['(%s)' % note, '', '']

    for r, (k, o) in chosen.items():
        rows = [[label] + abl_cells(m, note) for label, m, note in ablation_rows(results, k, o, extra[r])]
        rows.append(['R1: net gradient only, 20 values'] + abl_cells(curve_at[best_e if o == 'E' else best_j]))
        rows.append(['R0: Picasso defaults'] + abl_cells(curve_at[DEFAULT_NG]))
        P('Ablation of %s on the development frames (one structural option at a time put back to its default, the post-filters re-tuned by the same ascent for the same objective; the photon and width filters switched off one at a time with the rest of the chosen configuration kept; the chosen structure is `%s`):\n'
          % (r, skey(results[k]['struct'])))
        P(md_table(['configuration', 'dev J %', 'dev E', 'dev RMSE nm'], rows))
        P('')
    P('Each option alone on top of R1 (the default structure with only the net gradient tuned, no other filter; the net gradient is tuned for the objective of the column), development frames:\n')
    fE, fJ = forward_rows(results, 'E'), forward_rows(results, 'J')
    P(md_table(['setting', 'dev E (net gradient tuned for E)', 'dev J % (net gradient tuned for J)'],
               [[a[0], '%.1f' % a[1]['efficiency'], '%.1f' % (100 * b[1]['jaccard'])] for a, b in zip(fE, fJ)]))
    P('')

    # ---------------------------------------------------------------- the +c versions
    dd, dt = {}, {}
    for r in RUNGS:
        dd[r], dt[r] = dets_dev[r], dets_test[r]
        dd[r + '+c'], dt[r + '+c'] = C.shift_detections(dets_dev[r], *off[r]), C.shift_detections(dets_test[r], *off[r])

    # ---------------------------------------------------------------- test frames
    res, resc, ci, cic = {}, {}, {}, {}
    for v in VARIANTS:
        res[v], resc[v] = score_both(dt[v], gt_all, test_idx, pmin)
        ci[v] = C.bootstrap_ci(res[v], B=args.boot, seed=0, keys=C.BOOT_KEYS_EXT)
        cic[v] = C.bootstrap_ci(resc[v], B=args.boot, seed=0, keys=C.BOOT_KEYS_EXT)
    me = {v: C.matched_errors(dt[v], gt_all, test_idx) for v in VARIANTS}
    rg = gt_all[me[VARIANTS[0]]['gt_row']]
    kept = kept_mask(rg, pmin)
    common = {'all': np.stack([me[v]['matched'] for v in VARIANTS]).all(axis=0)}
    common['ch'] = common['all'] & kept
    rm, rmci = {'all': {}, 'ch': {}}, {'all': {}, 'ch': {}}
    for s_ in ('all', 'ch'):
        for v in VARIANTS:
            c = common_columns(me[v], common[s_], nT)
            ok = c['tp'].sum() > 0
            rm[s_][v] = float(np.sqrt(c['sse'].sum() / c['tp'].sum())) if ok else float('nan')
            rmci[s_][v] = C.bootstrap_pooled_ci(c, rmse_stat, B=args.boot, seed=0) if ok else (float('nan'), float('nan'))

    def table(resd, cid, rmd, rmcid):
        rows = [[v, ci_str(resd[v]['jaccard'], cid[v]['jaccard'], 100), ci_str(resd[v]['precision'], cid[v]['precision'], 100), ci_str(resd[v]['recall'], cid[v]['recall'], 100),
                 ci_str(resd[v]['rmse_nm'], cid[v]['rmse_nm']), ci_str(resd[v]['rmse_centered_nm'], cid[v]['rmse_centered_nm']), ci_str(rmd[v], rmcid[v]),
                 ci_str(resd[v]['efficiency'], cid[v]['efficiency'])] for v in VARIANTS]
        return md_table(['row', 'J %', 'precision %', 'recall %', 'RMSE nm', 'centred RMSE nm', 'RMSE common TP nm', 'efficiency'], rows)

    def diffs(resd):
        pairs = [('R1E', 'R0'), ('R1J', 'R0'), ('R2E', 'R0'), ('R2J', 'R0'), ('R2E', 'R1E'), ('R2J', 'R1J')]
        pairs = [(a + sfx, b + sfx) for sfx in ('', '+c') for a, b in pairs] + [(r + '+c', r) for r in RUNGS]
        rows = []
        for A, B in pairs:
            d = C.bootstrap_diff_ci(resd[A], resd[B], B=args.boot, seed=0)
            rows.append(['%s minus %s' % (A, B)] + ['%+.*f [%+.*f, %+.*f]' % (1, sc * (resd[A][key] - resd[B][key]), 1, sc * d[key][0], 1, sc * d[key][1])
                                                    for key, sc in (('jaccard', 100.0), ('rmse_nm', 1.0), ('efficiency', 1.0))])
        return md_table(['difference', 'Jaccard (points)', 'RMSE (nm)', 'efficiency (points)'], rows)

    P('### Test frames: all ground-truth emitters (primary setting)\n')
    P('Rungs R0 (Picasso defaults), R1E / R1J (net gradient tuned on the development frames for efficiency / Jaccard), R2E / R2J (all options tuned on the development frames for efficiency / '
      'Jaccard), each without and with the "+c" offset correction (the constant (dx, dy) estimated on the development frames for that configuration, subtracted from every detection). '
      '95 %% intervals of a bootstrap in blocks of %d consecutive frames (%d resamples, seed 0; pooled metrics recomputed in every resample). The centred RMSE is the RMSE about the mean error of '
      'the matched pairs of the test frames. "RMSE common TP": RMSE over the %d ground-truth emitters that EVERY one of the ten rows detects, so that no row is flattered by losing its hard emitters.\n'
      % (BLOCK, args.boot, int(common['all'].sum())))
    P(table(res, ci, rm['all'], rmci['all']))
    P('\nTP / FP / FN: %s. Detected per frame: %s; true per frame: %.3f. Offsets (dx, dy) in nm estimated on the development frames: %s.\n'
      % ('; '.join('%s %d / %d / %d' % (v, res[v]['tp'], res[v]['fp'], res[v]['fn']) for v in RUNGS), '; '.join('%s %.3f' % (v, res[v]['n_det_all'] / nT) for v in RUNGS),
         res['R0']['n_gt_all'] / nT, '; '.join('%s (%.1f, %.1f)' % (r, *off[r]) for r in RUNGS)))
    P('Paired differences between rows on the same test frames (both evaluated on every bootstrap resample, blocks of %d frames; 95 %% interval; an interval that contains 0 is not '
      'distinguishable from frame-sampling noise):\n' % BLOCK)
    P(diffs(res))
    P('')
    P('### Test frames: challenge assessment setting (emitters within 450 nm of the border and the dimmest 25 % of the whole set\'s activations dropped after pairing)\n')
    P('The same ten rows; "RMSE common TP" over the %d kept emitters that every row detects.\n' % int(common['ch'].sum()))
    P(table(resc, cic, rm['ch'], rmci['ch']))
    P('\nGround-truth emitters counted: %d; detections counted per row: %s.\n' % (resc['R0']['n_gt'], '; '.join('%s %d' % (v, resc[v]['n_det']) for v in RUNGS)))
    P('Paired differences (challenge setting):\n')
    P(diffs(resc))
    P('')
    stl = {r: count_stealing(dt[r], gt_all, test_idx, pmin) for r in RUNGS}
    P('Pairing precedes the exclusions, so an excluded (dim or border) emitter can take the detection of a kept one. Kept emitters left unmatched although a detection within 250 nm was '
      'paired with an excluded emitter, of all kept emitters: %s.\n' % '; '.join('%s %d of %d' % (r, *stl[r]) for r in RUNGS))
    rows = []
    for v in VARIANTS:
        cells = [v]
        for kw in ({}, dict(border_nm=CHALLENGE['border_nm']), dict(photon_min=pmin), dict(border_nm=CHALLENGE['border_nm'], photon_min=pmin)):
            q = C.score(dt[v], gt_all, test_idx, **kw)
            cells.append('%.1f / %.1f / %.1f' % (100 * q['jaccard'], q['rmse_nm'], q['efficiency']))
        rows.append(cells)
    P('### Effect of each exclusion on the same detections (Jaccard % / RMSE nm / efficiency, point estimates)\n')
    P(md_table(['row', 'all emitters', 'border 450 nm dropped', 'dimmest 25 % dropped', 'both (challenge setting)'], rows))
    P('\nPicasso never searches the box // 2 pixel rows / columns at the top and left (box // 2 + 1 at the bottom and right: 3 and 4 px for the default 7 px box) when the frame is not padded, '
      'so emitters there are missed by construction; the border setting removes most of that zone from the ground truth (450 nm = 4.5 px). Median padding (an R2 option) removes the blind zone.\n')

    # ---------------------------------------------------------------- recall tables
    edge = np.minimum(np.minimum(rg['x_nm'], C.FIELD_NM - rg['x_nm']), np.minimum(rg['y_nm'], C.FIELD_NM - rg['y_nm']))
    zz = np.abs(rg['z_nm'])
    inner = edge >= CHALLENGE['border_nm']

    def recall_table(header, masks):
        rows = [[lab, int(m.sum())] + ['%.1f' % (100 * me[r]['matched'][m].mean()) if m.sum() else 'nan' for r in RUNGS] for lab, m in masks]
        return md_table([header, 'emitters'] + ['recall %% %s' % r for r in RUNGS], rows)

    P('### Recall by distance to the border, by photons and by |z| (test frames, rungs without the offset correction)\n')
    P('By distance of the emitter to the nearest field edge (all emitters):\n')
    P(recall_table('distance to the edge', [('%d-%d nm' % (lo, hi), (edge >= lo) & (edge < hi)) for lo, hi in EDGE_BINS]))
    P('\nBy emitted photons (emitters at least 450 nm from the edge):\n')
    P(recall_table('photons', [('%d-%s photons' % (lo, '%d' % hi if hi < 1e8 else 'inf'), inner & (rg['photons'] >= lo) & (rg['photons'] < hi)) for lo, hi in PHOTON_BINS[noise]]))
    P('\nBy axial distance |z| from the focal plane (emitters at least 450 nm from the edge and at least %.0f photons; the PSF broadens with |z|, see `calibration.md`):\n' % pmin)
    P(recall_table('|z|', [('%d-%d nm' % (lo, hi), inner & (rg['photons'] >= pmin) & (zz >= lo) & (zz < hi)) for lo, hi in Z_BINS]))
    P('')

    # ---------------------------------------------------------------- photon bias and the structural ceiling
    P('### Picasso photons against the emitted photons, and what a one-peak-per-blob detector cannot do\n')
    rows = []
    for r in ('R0', 'R2E'):
        m = me[r]['matched']
        ratio = np.full(len(rg), np.nan)
        idx_m = np.flatnonzero(m)
        ratio[idx_m] = [dt[r][me[r]['pos'][i]][me[r]['det'][i], C.DET_COLUMNS.index('photons')] for i in idx_m] / rg['photons'][idx_m]
        for lab, mk in (('all matched', m), ('|z| 0-100 nm', m & (zz < 100)), ('|z| 100-200 nm', m & (zz >= 100) & (zz < 200)), ('|z| 200-350 nm', m & (zz >= 200) & (zz < 350)),
                        ('|z| 350-800 nm', m & (zz >= 350))):
            v = ratio[mk]
            rows.append([r, lab, int(mk.sum())] + (['%.2f' % np.median(v), '%.2f .. %.2f' % tuple(np.percentile(v, [25, 75]))] if len(v) else ['nan', '']))
    P('Photons reported by Picasso (QE-corrected) divided by the ground-truth photons of the matched emitter, test frames, median and inter-quartile range (the Gaussian core holds only part of the emitted '
      'photons, `calibration.md` section on the halo, so the photon count of a Picasso fit is biased low and more so out of focus):\n')
    P(md_table(['row', 'emitters', 'n', 'median ratio', 'IQR'], rows))
    nb, nrows = neighbour_fractions(gt_test, test_idx)
    P('\nStructural ceiling of a detector that reports one peak per blob: %.1f %% of the %d ground-truth rows of the test frames have another emitter of the same frame within 250 nm and %.1f %% within 400 nm; '
      'on the development frames, keeping EVERY candidate of the chosen structure (R2E: %.1f per frame, R2J: %.1f per frame) the recall is at most %.1f %% (R2E structure) / %.1f %% (R2J structure) '
      '(pairing with the maximum number of pairs). The development and test frames show the same three microtubules, and `positions.csv` holds the exact emitter sites of the fluorophores '
      '(no frame information), which a method that uses the structure as a prior could exploit; Picasso uses neither.\n'
      % (100 * nb[0], nrows, 100 * nb[1], ceiling['R2E'][1], ceiling['R2J'][1], 100 * ceiling['R2E'][0], 100 * ceiling['R2J'][0]))

    # ---------------------------------------------------------------- raw detections
    raw = {}
    for v in VARIANTS:
        for split, dsrc in (('dev', dd), ('test', dt)):
            raw['%s|%s|%s|off' % (name, v, split)], raw['%s|%s|%s|flat' % (name, v, split)] = C.pack_detections(dsrc[v], columns=(0, 1, 2, 3))
    for r in RUNGS:
        raw['%s|%s|offset_nm' % (name, r)] = np.array(off[r], dtype=np.float64)
        raw['%s|%s|config' % (name, r)] = np.array(json.dumps(C._merge_params(base_params[r]), sort_keys=True, default=C._json_default))
    raw['%s|dev_frames' % name] = np.array([d0, d1], dtype=np.int64)
    raw['%s|test_frames' % name] = np.array([t0, t1], dtype=np.int64)
    raw['%s|photon_threshold' % name] = np.array(pmin)
    summary = dict(res=res, resc=resc, ci=ci, cic=cic, rm=rm, rmci=rmci, off=off, n_struct=n_struct, n_eval=n_eval, keys=(keyE, keyJ), params=base_params)
    return L, raw, summary


# ==================================================================================================================================
# the report
# ==================================================================================================================================
RAW_README = ('Detections of the Picasso ladder (experiments/exp_smlm_picasso.py). Keys "<set>|<row>|<split>|off" (F + 1 int64 offsets) and "<set>|<row>|<split>|flat" '
              '(M x 4 float64: x_nm, y_nm, photons, sx_nm): the detections of frame i (0-based position in the split\'s frame list) are flat[off[i]:off[i + 1]]; '
              'rows R0, R1E, R1J, R2E, R2J as detected, and the same with "+c" (the constant offset of "<set>|<rung>|offset_nm", estimated on the development frames, '
              'subtracted from x_nm, y_nm); split "dev" or "test" (frame ranges in "<set>|dev_frames", "<set>|test_frames", [start, stop) stack indices); '
              '"<set>|<rung>|config" is the JSON of the Picasso options. positions in nm, corner convention; photons QE-corrected. Load with '
              'smlm_common.load_ladder_detections(path, set, row, split) (-> dets, frame_indices) and score with smlm_common.score(dets, gt, frame_indices).')


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--sets', nargs='+', default=list(SETS), choices=sorted(PROTOCOL))
    ap.add_argument('--boot', type=int, default=BOOT)
    ap.add_argument('--jobs', type=int, default=1, help='sets run in parallel (each with one single-threaded Picasso worker); at most 2 heavy processes are wanted')
    ap.add_argument('--picasso-python', default=os.environ.get('SMLM_PICASSO_PYTHON'))
    ap.add_argument('--cache-dir', default=None, help='development only: keep / reuse the result of the R2 search per set here (never used for the committed results)')
    ap.add_argument('--out', default=os.path.join(REPO, 'results', 'smlm', 'picasso_baseline.md'))
    ap.add_argument('--raw-out', default=os.path.join(REPO, 'results', 'smlm', 'picasso_ladder_raw.npz'))
    args = ap.parse_args(argv)
    t_all = time.time()
    info = {}
    C.run_picasso(np.full((2, 64, 64), 700.0), python=args.picasso_python, info=info)
    ck = {}
    try:
        ck = json.load(open(os.path.join(REPO, 'data', 'smlm', 'checksums.json')))
    except Exception:
        pass
    iss = iss_pic_grid()
    if iss is not None and iss.get('N2') != GRID_R1['N2']:
        print('WARNING: GRID_R1[N2] differs from PIC_GRID[N2] in exp_smlm_iss.py')
    out = ['# Picasso baseline ladder on the simulated 2D microtubule sets\n']
    out.append('Produced by `experiments/exp_smlm_picasso.py` with the harness `experiments/smlm_common.py`. Frame indices are 0-based stack indices, positions in nm '
               '(corner convention; Picasso\'s pixel-index coordinates are converted with nm = (p + 0.5) * 100), matching tolerance 250 nm (Hungarian pairing per frame, inclusive), '
               'metrics pooled over frames. Every number below is printed by that script; the file regenerates byte for byte (no timings).\n')
    out.append('**The baseline is a ladder, not one setting.** A baseline is only as fair as the options it may use (referee finding SM1-1, `review/REVIEW_roundSM1.md`): tuning the minimum net gradient alone '
               'leaves Picasso well below what its own documented options reach, and a single operating point is not a baseline (the point of best efficiency and the point of best Jaccard index differ strongly). '
               'Rungs, each tuned on development frames only and scored once on held-out test frames:\n')
    out.append('* **R0**: Picasso\'s defaults (minimum net gradient 5000, box 7 px, raw frames, no padding, least-squares fit `gausslq`).\n'
               '* **R1E / R1J**: the minimum net gradient alone tuned (20 values) for efficiency E / Jaccard J (the earlier "tuned (efficiency)" / "tuned (Jaccard)" rows).\n'
               '* **R2E / R2J**: every documented option tuned together for E / J: identification box, fit box, Picasso\'s Gaussian filter, median padding of each frame (Picasso cannot identify in a border zone '
               'of box // 2 px; the padded frame removes it, positions are shifted back and spots whose centre lies outside the field are dropped), fit method (`gausslq` / `gaussmle`), and post-filters '
               '(net gradient, minimum photons, width window |sx|, |sy|). The grids are declared at the top of the script before any run; the search is described in each set\'s section '
               'and the number of configurations tried is reported there, so that another method can be given the same tuning effort.\n'
               '* **+c** (appended to any row): the same detections with a constant position offset subtracted, estimated on the development frames only as the mean signed error (detection - truth) of the '
               'matched true positives of that configuration (the offset of a symmetric fit against the simulated emitter, `calibration.md`; it depends on the estimator and on |z|).\n')
    out.append('Two scoring settings on the same detections: ALL ground-truth emitters (primary), and the challenge\'s assessment setting (emitters closer than 450 nm to the border and the dimmest quarter of the '
               'whole set\'s activations dropped after pairing; the rules are from memory except the 250 nm tolerance and the metric definitions, see the docstring of `smlm_common.py`). '
               'Intervals are 95 %% percentile intervals of a bootstrap over blocks of %d consecutive frames (an emitter stays on for several consecutive frames, see the run lengths in each set\'s section and `calibration.md`); the earlier i.i.d. frame bootstrap was too narrow. '
               'They cover the sampling of frames only, not other datasets or the simulation, and the development and test frames share the same three microtubules. The tuning '
               'effort on the development frames is part of the result: a method that is tuned less is compared against a better-tuned baseline.\n' % BLOCK)
    out.append('**Exact Picasso setting.** Picasso %s (`picasso.localize.identify` and `fit`, in-memory stack, single thread; `localize_frames` is the same two calls with one box). Identification: local maxima over a box '
               'of 7 px (default) whose net gradient exceeds the minimum net gradient (default 5000, in raw DN, on the frame or, with the Gaussian filter, on the smoothed frame); fit: `gausslq`, Picasso\'s default least-squares fit '
               '(Levenberg-Marquardt, elliptical Gaussian with free sx, sy, amplitude, background, centre; convergence schedule of Picasso: tolerance %s, at most %s iterations), box 7 x 7. '
               'Camera passed to Picasso: baseline 100 DN, EM gain 300, sensitivity 45 e-/ADU, QE 0.9 (Picasso converts to photoelectrons (DN - baseline) x sensitivity / gain and does not divide by QE; the harness '
               'divides the reported photons by QE; the photon count of its Gaussian fit is biased low against the emitted photons because the Gaussian core does not hold the PSF halo, `calibration.md`; table in each set section). '
               '**The Gaussian of Picasso\'s least-squares fit is evaluated at the pixel centres, not integrated over the pixel** (the integrated-Gaussian fit exists only in its deprecated maximum-likelihood module); with the '
               'fitted sigma of about 1.3-1.4 px the pixel-integration correction is a width change of about 0.03 px. R0 and R1 have no post-filtering of the localizations (no sigma, photon or precision filter, no drift '
               'correction, no linking of consecutive frames). Cross-checks done once by hand, not by this script: Picasso\'s own command line (`picasso localize -a lq -b 7 -g 5000 -bl 100 -s 45 -ga 300 -qe 0.9 -px 100`, whose default '
               'fit would otherwise be `mle`) on a TIFF of test frames 1000-1199 gave exactly the localizations of this runner (158 spots, largest position difference 2e-4 nm); and the independent pipeline of the SM1 referee (Picasso identification and fit called separately, as here) gave the same R0 numbers and, for the '
               'configuration "identification box 7, fit box 9, Gaussian filter 1.5, median pad 5, net gradient > 2000" on the MT0.N1.LD test frames, the same J, precision, recall, RMSE and E as this runner (82.5, 99.4, 82.9, 45.7, 51.1).\n' % (info['picasso_version'], info.get('fit_tolerance', '0.01'), info.get('fit_max_iterations', '200')))
    out.append('Data: %s. Python %s, numpy %s.\n' % ('; '.join('%s sha256 %s...' % (k, v['sha256'][:12]) for k, v in sorted(ck.items()) if k.split('/')[0] in args.sets and ('activations' in k or 'stack' in k)) or 'see data/smlm/checksums.json',
                                                      platform.python_version(), np.__version__))
    jobs = max(1, int(args.jobs))
    names = list(args.sets)
    sections = {}
    raws = {}
    summaries = {}
    if jobs == 1 or len(names) == 1:
        for n in names:
            sections[n], raws[n], summaries[n] = run_set(n, args)
    else:
        ctx = multiprocessing.get_context('spawn')
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs, mp_context=ctx) as ex:
            futs = {n: ex.submit(_job, n, vars(args)) for n in names}
            for n in names:
                sections[n], raws[n], summaries[n] = futs[n].result()
    out.append('## Summary: test frames, every row of the ladder\n')
    out.append('Point estimate [95 % block-bootstrap interval] on the held-out test frames; the sections of the sets below give every metric, the paired differences, the development-frame search and the ablations.\n')
    for setting, key, ckey in (('all ground-truth emitters', 'res', 'ci'), ('challenge setting', 'resc', 'cic')):
        for metric, label, scale, nd in (('jaccard', 'Jaccard %', 100.0, 1), ('rmse_nm', 'RMSE nm', 1.0, 1), ('efficiency', 'efficiency E', 1.0, 1)):
            rows = []
            for v in VARIANTS:
                rows.append([v] + [ci_str(summaries[n][key][v][metric], summaries[n][ckey][v][metric], scale, nd) for n in names])
            out.append('%s, %s (test frames):\n' % (label, setting))
            out.append(md_table(['row'] + names, rows))
            out.append('')
    out.append('Number of configurations scored on the development frames (R2): ' + '; '.join('%s %d structural settings, %d (structure, post-filter) configurations' % (n, summaries[n]['n_struct'], summaries[n]['n_eval']) for n in names) + '.\n')
    for n in names:
        out.extend(sections[n])
    raw = {'README': np.array(RAW_README), 'columns': np.array(['x_nm', 'y_nm', 'photons', 'sx_nm'])}
    for n in names:
        raw.update(raws[n])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    C.save_npz_deterministic(args.raw_out, raw)
    digest = hashlib.sha256(open(args.raw_out, 'rb').read()).hexdigest()
    out.append('Detections of every row on the development and test frames: `results/smlm/picasso_ladder_raw.npz` (sets: %s; sha256 `%s`; format in its README entry).\n' % (', '.join(names), digest))
    with open(args.out, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    print('wrote', args.out, 'and', args.raw_out, '(%.0f s)' % (time.time() - t_all))


def _job(name, argd):
    ns = argparse.Namespace(**argd)
    return run_set(name, ns, log=lambda s: print(s, flush=True))


if __name__ == '__main__':
    main()
