"""
DECODE, a GPU-trained learned localizer, on the four 2D single-molecule challenge sets, against the sliding inverse scale space and Picasso
==========================================================================================================================================

The reviewer's question "DECODE (Speiser et al. 2021) is the standard learned localizer": `results/smlm/decode_run_spec.md` specified a run for the authors' GPU machine; this
script ingests a run that was made on a free cloud GPU (the protocol, versions and budget in `results/smlm/decode_run_plan.md`, written BEFORE the final run) and scores it with the
harness of `experiments/smlm_common.py` on the SAME test frames as the sliding ISS (`exp_smlm_iss.py`), Picasso (`exp_smlm_picasso.py`) and the CPU-scale CNN (`exp_smlm_learned.py`).

What was run (details and deviations: the plan and the "Protocol and deviations" section of the report)
    * DECODE at commit 4c0f38c3 (the version the spec was written against), one model per set, trained on DECODE's own simulator with the camera constants of the challenge (the same numbers
      the ISS and Picasso are given), a PSF that is a cubic-spline calibration built (no bead stack exists for the 2D modality) from the radial single-emitter profile measured on the
      DEVELOPMENT frames of MT0.N1.LD and their ground truth, the photon statistics and the background of the development frames of the set's noise class, the emitter density of the set.
      The test frames' ground truth enters no choice. The training budget was a fraction of the hours the authors train for: it is reported next to every number.
    * The checkpoint (one of three snapshots per set), the axis convention (DECODE's first coordinate is the image ROW axis) and, here, the probability threshold and a constant position offset
      are chosen on the development frames only: the threshold twice (by E and by J) as every other method, by `smlm_decode_ingest.evaluate`.
    * Candidates: every DECODE detection with probability >= 0.05 (DECODE's own post-processing with its detection threshold lowered), the development and the test frames of the set.

Usage
    python experiments/exp_smlm_decode.py --ingest RUN_DIR [--sets ...]    (RUN_DIR: the files returned by the GPU run: <set>_dev.csv, <set>_test.csv, <set>_select.json,
                                                                            <set>_train_summary.json, versions.txt, timing.json, gpu.txt and, optionally, run_extra.json = {"versions_text": [lines]}; needs the SMLM data and results/smlm/learned_raw.npz;
                                                                            writes decode_raw.npz and decode.md)
    python experiments/exp_smlm_decode.py --from-raw                       (regenerates decode.md from decode_raw.npz alone; numpy and scipy only; the report is byte-identical)
Environment: DECODE_OUT (directory of the raw and the md, default results/smlm), SMLM_DATA.

The raw holds, per set, the ground truth of the development and test frames, the candidate detections of DECODE on both (frame, x, y, photons, prob), the stored test detections of the
comparison rows, and the selection records of the GPU run; the thresholds, offsets and test detections are RECOMPUTED from the candidates by the report and checked against the stored
test detections (an inconsistent raw raises).
"""
import argparse
import hashlib
import json
import math
import os
import sys

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, HERE)
import smlm_common as C                       # noqa: E402
import smlm_decode_ingest as IN               # noqa: E402
from smlm_decode_export import protocol       # noqa: E402

OUT = os.environ.get('DECODE_OUT', os.path.join(REPO, 'results', 'smlm'))
RAW = os.path.join(OUT, 'decode_raw.npz')
MD = os.path.join(OUT, 'decode.md')
LEARNED_RAW = os.path.join(REPO, 'results', 'smlm', 'learned_raw.npz')
SETS = C.SETS
TAG = 'DECODE'
OPS = ('E', 'J')
BOOT, BLOCK = 1000, C.DEFAULT_BLOCK
CAND_COLS = ('frame', 'x_nm', 'y_nm', 'photons', 'prob')
CAND_DTYPE = np.dtype([(c, 'f8') for c in CAND_COLS])
COMPARATORS = (('ISS-0', 'sliding ISS, a priori default (no operating point chosen)'), ('ISS-E', 'sliding ISS, E-tuned on the development labels, offset removed'),
               ('R0', 'Picasso, defaults (R0)'), ('R2E+c', 'Picasso, every option tuned for E on the development labels + offset correction (R2E+c)'),
               ('CNN-a', 'CPU-scale CNN, supervised on the development labels (a)'), ('CNN-b', 'CPU-scale CNN, trained on simulation only (b)'))


# ----------------------------------------------------------------------------------------------------------------------------------
# conventions
# ----------------------------------------------------------------------------------------------------------------------------------
def decode_px_to_nm(x_px, y_px):
    """DECODE's pixel coordinates -> (x_nm, y_nm) in the corner convention of the ground truth. DECODE's pixel centres are at the integers and its FIRST coordinate runs along
    the first image axis, the ROWS (the PSF is drawn at frames[n, x_ix, y_ix], the post-processing adds the offsets to a meshgrid with 'ij' indexing); the ground truth's x is the COLUMN
    axis. So x_nm = (y_px + 0.5) * 100 and y_nm = (x_px + 0.5) * 100. (The GPU run chooses this convention or the plain one by the development frames; the choice is recorded.)"""
    return (np.asarray(y_px, dtype=np.float64) + 0.5) * C.PIXEL_NM, (np.asarray(x_px, dtype=np.float64) + 0.5) * C.PIXEL_NM


def to_struct(rows):
    """(n, 5) float array (frame, x_nm, y_nm, photons, prob) -> the structured array `smlm_decode_ingest` works on."""
    rows = np.asarray(rows, dtype=np.float64).reshape(-1, 5)
    out = np.zeros(len(rows), dtype=CAND_DTYPE)
    for i, c in enumerate(CAND_COLS):
        out[c] = rows[:, i]
    return out


def read_candidates(path):
    """A candidate CSV of the GPU run as an (n, 5) float array."""
    return np.array([[r[c] for c in CAND_COLS] for r in IN.read_csv(path)], dtype=np.float64).reshape(-1, 5)


def pack_cand(rows, n_frames):
    """Candidates -> (off, flat): off (n_frames + 1,) offsets, flat (m, 4) = x_nm, y_nm, photons, prob of the rows of frame i at flat[off[i]:off[i + 1]] (rows sorted by frame, stable)."""
    rows = np.asarray(rows, dtype=np.float64).reshape(-1, 5)
    fr = rows[:, 0].astype(np.int64)
    if len(fr) and (fr.min() < 0 or fr.max() >= n_frames):
        raise ValueError('frame index outside 0..%d' % (n_frames - 1))
    order = np.argsort(fr, kind='stable')
    rows, fr = rows[order], fr[order]
    off = np.zeros(n_frames + 1, dtype=np.int64)
    np.add.at(off, fr + 1, 1)
    return np.cumsum(off), rows[:, 1:5].copy()


def unpack_cand(off, flat):
    """Inverse of `pack_cand`: the (n, 5) candidate array."""
    off = np.asarray(off)
    fr = np.repeat(np.arange(len(off) - 1), np.diff(off))
    return np.c_[fr, np.asarray(flat, dtype=np.float64)]


def gt_array(g):
    """A structured ground-truth array (smlm_common.GT_DTYPE) from a dict of columns."""
    a = np.zeros(len(g['frame']), dtype=C.GT_DTYPE)
    a['id'] = np.arange(1, len(a) + 1)
    for k in ('frame', 'x_nm', 'y_nm', 'z_nm', 'photons'):
        a[k] = g[k]
    return a


# ----------------------------------------------------------------------------------------------------------------------------------
# the raw
# ----------------------------------------------------------------------------------------------------------------------------------
def load_npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def set_result(name, dev_cand, test_cand, gt):
    """The development-chosen operating points of DECODE on one set: `smlm_decode_ingest.evaluate` on the stored candidates. Returns dict(op -> dict(thr, off, dev, test, test_dets))."""
    return IN.evaluate(name, to_struct(dev_cand), to_struct(test_cand), gt)


def build_raw(sets_in, meta, comparators=None):
    """The arrays of the raw from per-set inputs: sets_in[name] = dict(gt=dict of columns for the dev+test frames, dev=(n,5) candidates, test=(n,5) candidates, select=json str, train=json str);
    comparators[name] = dict(tag -> list of per-frame (n, 4) detection arrays) and the strings/arrays under extra keys (`cmp_extra`). The thresholds and the test detections are computed here
    and stored too (the ingest key names), so that `smlm_decode_ingest` and the report can be checked against each other."""
    arrays = {}
    for name, s in sets_in.items():
        (da, db), (ta, tb) = protocol(name)['dev'], protocol(name)['test']
        arrays['%s|dev_frames' % name], arrays['%s|test_frames' % name] = np.array([da, db]), np.array([ta, tb])
        for k in ('frame', 'x_nm', 'y_nm', 'z_nm', 'photons'):
            arrays['%s|gt|%s' % (name, k)] = np.asarray(s['gt'][k])
        for split, (a, b) in (('dev', (da, db)), ('test', (ta, tb))):
            off, flat = pack_cand(s[split], b - a)
            arrays['%s|%s|%s|cand|off' % (name, TAG, split)], arrays['%s|%s|%s|cand|flat' % (name, TAG, split)] = off, flat
        res = set_result(name, s['dev'], s['test'], gt_array(s['gt']))
        for op, o in res.items():
            off, flat = IN.pack(o['test_dets'])
            arrays['%s|%s-%s|test|off' % (name, TAG, op)], arrays['%s|%s-%s|test|flat' % (name, TAG, op)] = off, flat
        arrays['%s|select' % name] = np.array(s.get('select', ''))
        arrays['%s|train' % name] = np.array(s.get('train', ''))
        arrays['%s|pmin' % name] = np.array(s['pmin'])
        for tag, dets in (comparators or {}).get(name, {}).items():
            if tag.endswith('|config') or tag.endswith('|secs'):
                arrays['%s|cmp|%s' % (name, tag)] = np.asarray(dets)
                continue
            off, flat = C.pack_detections(dets, columns=(0, 1, 2, 3))
            arrays['%s|cmp|%s|off' % (name, tag)], arrays['%s|cmp|%s|flat' % (name, tag)] = off, flat
    arrays['meta'] = np.array(json.dumps(meta, sort_keys=True))
    return arrays


def cand_of(raw, name, split):
    return unpack_cand(raw['%s|%s|%s|cand|off' % (name, TAG, split)], raw['%s|%s|%s|cand|flat' % (name, TAG, split)])


def gt_of(raw, name):
    return gt_array({k: raw['%s|gt|%s' % (name, k)] for k in ('frame', 'x_nm', 'y_nm', 'z_nm', 'photons')})


def check_raw(raw, sets):
    """The stored operating-point detections equal those recomputed from the stored candidates (thresholds chosen on the development frames); raises otherwise. Returns the results."""
    out = {}
    for name in sets:
        res = set_result(name, cand_of(raw, name, 'dev'), cand_of(raw, name, 'test'), gt_of(raw, name))
        for op in OPS:
            off, flat = IN.pack(res[op]['test_dets'])
            if not (np.array_equal(off, raw['%s|%s-%s|test|off' % (name, TAG, op)]) and np.array_equal(flat, raw['%s|%s-%s|test|flat' % (name, TAG, op)])):
                raise ValueError('%s %s: the stored test detections differ from those recomputed from the candidates' % (name, op))
        out[name] = res
    return out


# ----------------------------------------------------------------------------------------------------------------------------------
# report
# ----------------------------------------------------------------------------------------------------------------------------------
def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    return '\n'.join(out + ['| ' + ' | '.join(str(c) for c in r) + ' |' for r in rows])


def fmt_ci(v, ci, scale=100.0, nd=1):
    return '%.*f [%.*f, %.*f]' % (nd, scale * v, nd, scale * ci[0], nd, scale * ci[1])


def pct(v, ci, nd=1):
    return fmt_ci(v, ci, 100.0, nd)


def f1(v, ci):
    return fmt_ci(v, ci, 1.0, 1)


def verdict(d, ci):
    """'win' if the paired interval of (DECODE minus comparator) excludes 0 in DECODE's favour, 'loss' against it, else 'tie'."""
    return 'win' if ci[0] > 0 else ('loss' if ci[1] < 0 else 'tie')


def cmp_dets(raw, name, tag):
    k = '%s|cmp|%s' % (name, tag)
    if k + '|off' not in raw:
        return None
    return C.unpack_detections(raw[k + '|off'], raw[k + '|flat'])


def method_table(raw, name, res):
    """Ordered (key, label, detections) of the methods scored on one set."""
    ms = [('DECODE-E', 'DECODE, probability threshold and offset chosen for E on the development frames', res['E']['test_dets']),
          ('DECODE-J', 'DECODE, probability threshold and offset chosen for J on the development frames', res['J']['test_dets'])]
    for tag, lab in COMPARATORS:
        d = cmp_dets(raw, name, tag)
        if d is not None:
            ms.append((tag, lab, d))
    return ms


def report(raw, boot=BOOT):
    meta = json.loads(str(raw['meta']))
    sets = meta['sets']
    results = check_raw(raw, sets)
    L = ['# DECODE, a GPU-trained learned localizer, against the sliding inverse scale space and Picasso on the 2D challenge sets', '',
         'Generated by `experiments/exp_smlm_decode.py --from-raw` from `decode_raw.npz`; every number below is printed by the script from that file (the ground truth, the candidate detections of DECODE on the '
         'development and the test frames and the detections of the comparison rows are stored in it; the operating points are recomputed from the candidates and checked against the stored test detections). '
         'The protocol, versions and budget were written to `decode_run_plan.md` before the final run. Frame indices are 0-based stack indices, positions in nm (corner convention), metrics pooled over the frames '
         'with 95 %% intervals of a bootstrap in blocks of %d consecutive frames (%d resamples, seed 0); differences are paired (the same resampled blocks). '
         'DECODE commit %s; %s.' % (BLOCK, boot, meta.get('decode_commit', '?'), meta.get('software', '?')), '']
    R, dets_by = {}, {}
    for name in sets:
        gt = gt_of(raw, name)
        a, b = (int(v) for v in raw['%s|test_frames' % name])
        idx = np.arange(a, b)
        pmin = float(raw['%s|pmin' % name])
        R[name], dets_by[name] = {}, {}
        for key, lab, dets in method_table(raw, name, results[name]):
            assert len(dets) == len(idx), (name, key, len(dets), len(idx))
            dets_by[name][key] = dets
            r = C.score(dets, gt, idx)
            r2 = C.score(dets, gt, idx, border_nm=450.0, photon_min=pmin)
            R[name][key] = (lab, r, r2, C.bootstrap_ci(r, B=boot, seed=0, block=BLOCK), C.bootstrap_ci(r2, B=boot, seed=0, block=BLOCK))
    L += summary_lines(raw, R, sets, boot, meta)
    for name in sets:
        L += set_lines(raw, name, R[name], results[name], dets_by[name], boot, meta)
    L += protocol_lines(raw, meta)
    return '\n'.join(L) + '\n'


def summary_lines(raw, R, sets, boot, meta):
    L = ['## Summary: efficiency E and Jaccard index J on the test frames (point estimates; intervals in the per-set tables)', '',
         'E = 100 - sqrt((100 - J)^2 + RMSE^2) (J in %%, RMSE in nm), all ground-truth emitters counted. DECODE and the tuned rows (ISS E-tuned + offset, Picasso R2E+c, CNN (a)) all use the '
         'development labels: DECODE for its PSF profile, photon statistics, checkpoint, probability threshold and a constant offset; the ISS a priori default, Picasso R0 and CNN (b) are label-free on the frames. '
         'DECODE\'s training budget was %s per set on a cloud GPU (%s), a small fraction of what the authors use.' % (meta.get('budget_text', '?'), meta.get('gpu', '?')), '']
    keys = [('DECODE-E', 'DECODE (E-tuned threshold)'), ('DECODE-J', 'DECODE (J-tuned threshold)'), ('ISS-E', 'ISS E-tuned + offset'), ('R2E+c', 'Picasso R2E+c'), ('CNN-a', 'CNN (a) supervised'),
            ('ISS-0', 'ISS a priori default'), ('R0', 'Picasso defaults R0'), ('CNN-b', 'CNN (b) simulation-trained')]
    rows = []
    for k, lab in keys:
        row = [lab]
        for s in sets:
            row.append('%.1f (%.1f)' % (R[s][k][1]['efficiency'], 100 * R[s][k][1]['jaccard']) if k in R[s] else '-')
        rows.append(row)
    full = [s_ for s_ in sets if all(k in R[s_] for k, _ in keys[:4])]
    if len(full) > 1:
        rows.append(['mean E over the %d sets: DECODE E-tuned / ISS E-tuned / Picasso R2E+c' % len(full), ' / '.join('%.1f' % np.mean([R[s_][k][1]['efficiency'] for s_ in full]) for k in ('DECODE-E', 'ISS-E', 'R2E+c'))] + [''] * (len(sets) - 1))
    L += [table(['method'] + [s + ' E (J %)' for s in sets], rows), '']
    L += ['Paired difference of E (DECODE E-tuned minus comparator, all emitters) with its 95 % interval and the verdict (win: the interval excludes 0 in DECODE\'s favour; loss: against it; tie: it contains 0). '
          'The first two rows compare like with like (both tuned on the development labels).', '']
    rows = []
    for kb, note in (('ISS-E', 'like-for-like (both use the development labels)'), ('R2E+c', 'like-for-like (both use the development labels)'), ('CNN-a', 'like-for-like (both use the development labels)'),
                     ('ISS-0', 'DECODE has the labels'), ('R0', 'DECODE has the labels')):
        row = ['DECODE - %s' % kb, note]
        for s in sets:
            if kb in R[s]:
                ra, rb = R[s]['DECODE-E'][1], R[s][kb][1]
                d = C.bootstrap_diff_ci(ra['per_frame'], rb['per_frame'], B=boot, seed=0, block=BLOCK)['efficiency']
                row.append('%s %s' % (f1(ra['efficiency'] - rb['efficiency'], d), verdict(ra['efficiency'] - rb['efficiency'], d)))
            else:
                row.append('-')
        rows.append(row)
    L += [table(['comparison', 'fairness'] + list(sets), rows), '']
    tallies = []
    for kb in ('ISS-E', 'R2E+c'):
        v = [verdict(0, C.bootstrap_diff_ci(R[s]['DECODE-E'][1]['per_frame'], R[s][kb][1]['per_frame'], B=boot, seed=0, block=BLOCK)['efficiency']) for s in sets if kb in R[s]]
        tallies.append('against %s: %d win(s), %d tie(s), %d loss(es)' % (kb, v.count('win'), v.count('tie'), v.count('loss')))
    L += ['Tally over the sets of the like-for-like comparisons (E, all emitters, DECODE E-tuned minus comparator): %s.' % '; '.join(tallies), '']
    return L


def op_line(name, op, o):
    return '%s-tuned: probability >= %.2f, offset removed (%.1f, %.1f) nm, development J %.1f %%, RMSE %.1f nm, E %.1f' % (
        op, o['thr'], o['off'][0], o['off'][1], 100 * o['dev']['jaccard'], o['dev']['rmse_nm'], o['dev']['efficiency'])


def set_lines(raw, name, R, res, dets_by, boot, meta):
    gt = gt_of(raw, name)
    a, b = (int(v) for v in raw['%s|test_frames' % name])
    idx = np.arange(a, b)
    tg = gt[(gt['frame'] >= a) & (gt['frame'] < b)]
    sel = json.loads(str(raw['%s|select' % name])) if str(raw['%s|select' % name]) else {}
    tr = json.loads(str(raw['%s|train' % name])) if str(raw['%s|train' % name]) else {}
    L = ['## %s' % name, '',
         'Test frames %d-%d (%d frames, %.3f ground-truth emitters per frame). DECODE candidates: %d on the development frames, %d on the test frames (probability >= 0.05). Development-chosen operating points: %s; %s.\n' % (
             a, b - 1, b - a, len(tg) / float(b - a), len(cand_of(raw, name, 'dev')), len(cand_of(raw, name, 'test')), op_line(name, 'E', res['E']), op_line(name, 'J', res['J']))]
    L.append('### Training, selection and cost\n')
    if tr or sel:
        rows = [['training epochs run (stopped by the clock; budget %s s)' % tr.get('cfg', {}).get('budget_s', '?'), tr.get('epochs_done', '?')], ['training wall clock (s)', '%.0f' % tr['seconds'] if 'seconds' in tr else '?'],
                ['simulated training density (emitters per frame)', tr.get('cfg', {}).get('emitter_av', '?')], ['photon statistics of the simulation (mean, sd)', '(%.0f, %.0f)' % (tr['cfg']['photons_mean'], tr['cfg']['photons_sd']) if 'cfg' in tr else '?'],
                ['background of the simulation (photons, +-10 %)', '%.2f' % tr['cfg']['bg_photons'] if 'cfg' in tr else '?'], ['snapshot chosen on the development frames', sel.get('chosen_snapshot', '?')],
                ['axis convention chosen on the development frames', 'x_nm from DECODE y (rows are DECODE x)' if sel.get('swap_chosen') else 'plain', ],
                ['best development E of the two conventions (swap / plain)', '%.1f / %.1f' % (sel['best_dev_E_by_axis_convention']['x=DECODE_y (swap)'], sel['best_dev_E_by_axis_convention']['x=DECODE_x (no swap)']) if sel else '?'],
                ['development E of the chosen snapshot: the run\'s own scorer / the project\'s (E-tuned point)', '%.1f / %.1f' % (sel['chosen_dev']['E'], res['E']['dev']['efficiency']) if sel.get('chosen_dev') else '?']]
        L.append(table(['item', 'value'], rows))
        if sel.get('all'):
            L.append('\nDevelopment frames, best probability threshold of each snapshot (offset-corrected E, as the notebook\'s own scorer computed it; the convention chosen above):\n')
            rws = [[r['snap'], r['n_candidates'], '%.2f' % r['best']['thr'], '%.1f' % (100 * r['best']['J']), '%.1f' % r['best']['rmse'], '%.1f' % r['best']['E']] for r in sel['all'] if r['swap'] == sel['swap_chosen']]
            L.append(table(['snapshot', 'candidates', 'threshold', 'dev J %', 'dev RMSE nm', 'dev E'], rws))
        L.append('')
    L.append('### Test frames: all ground-truth emitters (primary setting)\n')
    L.append(table(['method', 'J %', 'precision %', 'recall %', 'RMSE nm', 'efficiency', 'bias (dx, dy) nm', 'detected / frame'],
                   [[m, pct(r['jaccard'], ci['jaccard']), pct(r['precision'], ci['precision']), pct(r['recall'], ci['recall']), f1(r['rmse_nm'], ci['rmse_nm']), f1(r['efficiency'], ci['efficiency']),
                     '(%.1f, %.1f)' % (r['bias_x_nm'], r['bias_y_nm']), '%.3f' % (r['n_det'] / len(idx))] for m, r, r2, ci, ci2 in R.values()]))
    L.append('\nChallenge assessment setting (emitters within 450 nm of the border and the dimmest quarter of the activations dropped after the pairing):\n')
    L.append(table(['method', 'J %', 'RMSE nm', 'efficiency'], [[m, pct(r2['jaccard'], ci2['jaccard']), f1(r2['rmse_nm'], ci2['rmse_nm']), f1(r2['efficiency'], ci2['efficiency'])] for m, r, r2, ci, ci2 in R.values()]))
    L.append('\n### Paired differences (DECODE minus comparator, same frames, same resampled blocks), 95 % intervals\n')
    rows = []
    for ka in ('DECODE-E', 'DECODE-J'):
        for kb, note in (('ISS-E', 'both use the development labels'), ('R2E+c', 'both use the development labels'), ('CNN-a', 'both use the development labels'), ('ISS-0', 'DECODE has the labels'),
                         ('R0', 'DECODE has the labels'), ('CNN-b', 'DECODE has the labels')):
            if kb not in R:
                continue
            ra, rb = R[ka], R[kb]
            d = C.bootstrap_diff_ci(ra[1]['per_frame'], rb[1]['per_frame'], B=boot, seed=0, block=BLOCK)
            d2 = C.bootstrap_diff_ci(ra[2]['per_frame'], rb[2]['per_frame'], B=boot, seed=0, block=BLOCK)
            rows.append(['%s - %s' % (ka, kb), note, pct(ra[1]['jaccard'] - rb[1]['jaccard'], d['jaccard']), f1(ra[1]['rmse_nm'] - rb[1]['rmse_nm'], d['rmse_nm']), f1(ra[1]['efficiency'] - rb[1]['efficiency'], d['efficiency']),
                         f1(ra[2]['efficiency'] - rb[2]['efficiency'], d2['efficiency'])])
    L.append(table(['comparison', 'fairness', 'dJ % (all)', 'dRMSE nm (all)', 'dE (all)', 'dE (challenge)'], rows))
    L.append('')
    L += recall_lines(R, dets_by, tg, idx)
    L += z_lines(R, dets_by, tg, idx)
    L += threshold_lines(raw, name, tg, idx, res)
    return L


def recall_lines(R, dets_by, gt, idx):
    """Recall by the photon count of the ground-truth emitter, in bins of the test frames' own photon quantiles (descriptive)."""
    ph = np.asarray(gt['photons'])
    edges = np.quantile(ph, [0.0, 0.125, 0.25, 0.5, 0.75, 1.0])
    edges[-1] = np.inf
    cols = ['%.0f-%s photons (n = %d)' % (lo, '%.0f' % hi if np.isfinite(hi) else 'max', int(((ph >= lo) & (ph < hi)).sum())) for lo, hi in zip(edges[:-1], edges[1:])]
    rows = []
    for key, (m, r, r2, ci, ci2) in R.items():
        mm = C.matched_errors(dets_by[key], gt, idx)
        p_ = ph[mm['gt_row']]
        rows.append([m] + ['%.1f' % (100.0 * mm['matched'][(p_ >= lo) & (p_ < hi)].mean()) for lo, hi in zip(edges[:-1], edges[1:])])
    return ['### Recall by the photon count of the emitter (test frames; bins are quantiles of the test frames\' ground-truth photons)\n', table(['method (recall %)'] + cols, rows), '']


def z_lines(R, dets_by, gt, idx):
    """Recall and RMSE of the matched pairs by |z| (the PSF broadens with |z|)."""
    if not np.any(gt['z_nm'] != 0):
        return []
    z = np.abs(np.asarray(gt['z_nm']))
    edges = [0.0, 100.0, 200.0, 350.0, np.inf]
    cols = ['|z| %s-%s nm (n = %d)' % ('%.0f' % lo, '%.0f' % hi if np.isfinite(hi) else 'max', int(((z >= lo) & (z < hi)).sum())) for lo, hi in zip(edges[:-1], edges[1:])]
    rows = []
    for key, (m, r, r2, ci, ci2) in R.items():
        mm = C.matched_errors(dets_by[key], gt, idx)
        zz = z[mm['gt_row']]
        cells = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s_ = (zz >= lo) & (zz < hi)
            ok = s_ & mm['matched']
            cells.append('%.1f / %.1f' % (100.0 * mm['matched'][s_].mean(), math.sqrt(float(np.mean(mm['dx'][ok] ** 2 + mm['dy'][ok] ** 2)))) if ok.any() else '- / -')
        rows.append([m] + cells)
    return ['### Recall % / RMSE nm of the matched pairs by the axial position of the emitter (test frames; descriptive)\n', table(['method'] + cols, rows), '']


def threshold_lines(raw, name, gt, idx, res):
    """Descriptive only: the test-frame efficiency of DECODE at every threshold of the grid (the E-chosen offset removed). The thresholds were chosen on the development frames; the last column uses the
    test labels and is not an operating point."""
    cand = cand_of(raw, name, 'test')
    off = res['E']['off']
    rows, es = [], []
    for t in IN.THRESHOLDS:
        d = IN.shift(IN.per_frame(to_struct(cand), len(idx), t), off)
        es.append(C.score(d, gt_of(raw, name), idx)['efficiency'])
    best = max(range(len(es)), key=lambda i: -np.inf if math.isnan(es[i]) else es[i])
    rows.append(['DECODE, offset of the E point'] + [('**%.1f**' if t == res['E']['thr'] else '%.1f') % e for t, e in zip(IN.THRESHOLDS, es)] + ['%.2f (E %.1f)' % (IN.THRESHOLDS[best], es[best])])
    return ['### Threshold sensitivity on the test frames (descriptive only; the E-chosen threshold in bold)\n', table(['network'] + ['thr %.2f' % t for t in IN.THRESHOLDS] + ['best in hindsight'], rows),
            '\nThe last column uses the test labels and is not an operating point.\n']


def protocol_lines(raw, meta):
    L = ['## Protocol and deviations', '']
    for s in meta.get('protocol_text', []):
        L.append('* ' + s)
    L += ['', '## Versions and cost', '']
    for s in meta.get('versions_text', []):
        L.append('* ' + s)
    return L + ['']


# ----------------------------------------------------------------------------------------------------------------------------------
# ingest of a GPU run
# ----------------------------------------------------------------------------------------------------------------------------------
def read_text(path):
    return open(path).read().strip() if os.path.exists(path) else ''


def comparators_from_learned(name, path=LEARNED_RAW):
    """The stored test detections of the comparison rows of `learned_raw.npz` (ISS a priori default, ISS E-tuned, Picasso R0 and R2E+c, the CPU CNN (a) and (b)), 4 columns each, plus the ISS config."""
    z = load_npz(path)
    out = {}
    for tag in ('ISS-0', 'ISS-E', 'R0', 'R2E+c'):
        out[tag] = C.unpack_detections(z['%s|cmp|%s|off' % (name, tag)], z['%s|cmp|%s|flat' % (name, tag)])
    for tag, k in (('CNN-a', 'sup'), ('CNN-b', 'sim')):
        key = '%s|%s' % (name, k)
        if key + '|sel' in z:
            sel = json.loads(str(z[key + '|sel']))
            d = [c[c[:, 3] >= sel['thr']] if len(c) else c for c in C.unpack_detections(z[key + '|off'], z[key + '|flat'])]
            if sel['tag'].startswith('sup'):
                d = IN.shift(d, np.array(sel['offset']))
            out[tag] = d
    return out, z


DECODE_COMMIT = '4c0f38c33681d1c3005bc18bd99a11f46d30f1c1'


def load_json(path):
    t = read_text(path)
    return json.loads(t) if t else {}


def compose_meta(run_dir, sets):
    """The text of the report that comes from the files returned by the GPU run (versions, GPU, per-set epochs and seconds) and the fixed statements of the plan."""
    ver = read_text(os.path.join(run_dir, 'versions.txt')).splitlines()
    gpu = read_text(os.path.join(run_dir, 'gpu.txt')).splitlines()
    gpu_name = gpu[1].split(',')[0] if len(gpu) > 1 else 'unknown GPU'
    timing = load_json(os.path.join(run_dir, 'timing.json'))
    trs = {s_: load_json(os.path.join(run_dir, '%s_train_summary.json' % s_)) for s_ in sets}
    secs = [t['seconds'] for t in trs.values() if t]
    eps = {s_: t.get('epochs_done') for s_, t in trs.items() if t}
    meta = dict(tag=TAG, sets=list(sets), thresholds=list(IN.THRESHOLDS), boot=BOOT, block=BLOCK, decode_commit=DECODE_COMMIT[:8] + ' (' + DECODE_COMMIT + ')',
                software=' '.join(ver[:1]) if ver else '?', gpu='%s; one GPU per set, two sets trained at a time' % gpu_name,
                budget_text=('%.0f-%.0f s (%s epochs of 10,000 simulated frames)' % (min(secs), max(secs), '/'.join(str(eps[s_]) for s_ in sets if s_ in eps))) if secs else '?')
    meta['protocol_text'] = [
        'The plan (`decode_run_plan.md`) was written before the run (about 09:22 UTC on 2026-11-09; the run was submitted at 09:23:25 UTC). Checkpoint, axis convention, probability threshold and offset were chosen on the development frames only; '
        'the test frames were run once, with the selected snapshot. One training run (one configuration) per set.',
        'Inputs of DECODE: the camera constants of the dataset pages (QE 0.90 set; DECODE\'s default is 1), the quiet-frame background, the photon statistics and the 97 % |z| quantile of the development ground truth of the set, the emitter density of '
        'the set, and a PSF calibration built from the radial single-emitter profile of the development frames of MT0.N1.LD (51 x 51 px, 41 z cells of 50 nm, 1.42 photons per ground-truth photon, symmetric). A symmetric PSF model measured on 360 '
        'single-emitter frames is a handicap that a bead calibration would not have; the HD sets reuse the PSF and the background of the LD sets of their noise class.',
        'Deviations from `decode_run_spec.md`: (1) the pinned commit needs torch 1.7.1 (its `environment.yaml`): used, with cuda 11.0, on Tesla T4 GPUs; (2) the spec\'s `PostProcessingParam.raw_th: 0.05` is not read by the commit\'s `infer.py` (it hard-codes 0.1): '
        'the run uses DECODE\'s own inference pipeline in a script that copies `infer.py` with the threshold 0.05; (3) the PSF is option 1 of the spec in a radially symmetric form; (4) axial positions uniform on [0, 600] nm (the PSF is even in z); (5) '
        '`max_number_targets` 40 / 60 instead of 250; (6) the learning rate decays by a factor 10 over the wall-clock budget instead of DECODE\'s step schedule, and training is stopped by the clock; (7) the x and y of DECODE are the row and column '
        'axes (x_nm = (y_px + 0.5) 100): the run chose the convention on the development frames and recorded the alternative.',
        'Not comparable with the leaderboard results of DECODE\'s authors (hidden test data of the challenge; the MT0 sets are simulated training sets that share one structure of three microtubules and are not independent datasets). A network trained '
        'on emitters scattered uniformly over the field has not seen the microtubule structure; the ISS and Picasso do not use it either.',
        'Cost: the ISS needs no training and ~1 s (LD) to ~5 s (HD) per frame on one CPU core; DECODE trained for %s per set here (hours on a GPU in the spec) and then fits a frame in milliseconds on a GPU.' % meta['budget_text']]
    meta['versions_text'] = ['DECODE commit %s (%s); packages: %s' % (DECODE_COMMIT, ' '.join(ver[:1]), '; '.join(' '.join(l.split()).replace('.json', '') for l in ver[1:] if l))] + [
        'GPU: %s' % ' / '.join(gpu[1:3]) if len(gpu) > 2 else 'GPU: ?']
    for s_ in sets:
        t = trs[s_]
        if t:
            cfg = t.get('cfg', {})
            meta['versions_text'].append('%s: %s epochs in %.0f s (budget %.0f s), seed %s, one training run; simulated density %s per frame, z in [0, %.0f] nm, max targets %s' % (
                s_, t.get('epochs_done'), t['seconds'], cfg.get('budget_s', float('nan')), cfg.get('seed'), cfg.get('emitter_av'), cfg.get('z_hi', float('nan')), cfg.get('max_targets')))
    if timing:
        meta['versions_text'].append('GPU session: %.0f s in total (environment installation, data, calibration, training of two waves, selection and inference); earlier short timing runs of the same pipeline: 860 s and 640 s' % timing.get('total_s', float('nan')))
    meta['versions_text'] += list(load_json(os.path.join(run_dir, 'run_extra.json')).get('versions_text', []))          # lines the person who ran it adds (total GPU seconds of the account, notebook bookkeeping)
    return meta


def ingest(run_dir, sets, out_path=RAW, log=print, gt_loader=None, learned_path=LEARNED_RAW):
    """Read the files returned by the GPU run, build the raw and save it. `gt_loader(name)` returns the structured ground truth of a set (default: the SMLM data of `smlm_common`)."""
    meta = compose_meta(run_dir, sets)
    gt_loader = gt_loader or (lambda name: C.load_gt(C.data_path(name, 'gt')))
    sets_in, comps = {}, {}
    cmp_extra = {}
    for name in sets:
        (da, db), (ta, tb) = protocol(name)['dev'], protocol(name)['test']
        gt = gt_loader(name)
        g = gt[(gt['frame'] >= da) & (gt['frame'] < tb)]
        dev = read_candidates(os.path.join(run_dir, '%s_dev.csv' % name))
        test = read_candidates(os.path.join(run_dir, '%s_test.csv' % name))
        comps[name], z = comparators_from_learned(name, learned_path)
        cmp_extra[name] = np.asarray(z['%s|cmp|ISS-E|config' % name])
        sets_in[name] = dict(gt={k: np.asarray(g[k]) for k in ('frame', 'x_nm', 'y_nm', 'z_nm', 'photons')}, dev=dev, test=test, pmin=float(z['%s|pmin' % name]),
                             select=read_text(os.path.join(run_dir, '%s_select.json' % name)), train=read_text(os.path.join(run_dir, '%s_train_summary.json' % name)))
        log('%s: %d dev and %d test candidates' % (name, len(dev), len(test)))
    arrays = build_raw(sets_in, meta, comps)
    for name in sets:
        arrays['%s|cmp|ISS-E|config' % name] = cmp_extra[name]
    C.save_npz_deterministic(out_path, arrays)
    log('wrote %s (%.1f MB)' % (out_path, os.path.getsize(out_path) / 1e6))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ingest', metavar='RUN_DIR')
    ap.add_argument('--from-raw', action='store_true')
    ap.add_argument('--sets', nargs='*', default=list(SETS))
    A = ap.parse_args()
    if A.ingest:
        ingest(A.ingest, A.sets)
    if A.from_raw or A.ingest:
        raw = load_npz(RAW)
        text = report(raw)
        with open(MD, 'w') as f:
            f.write(text)
        print('wrote', MD)
    if not (A.ingest or A.from_raw):
        ap.print_help()


if __name__ == '__main__':
    main()
