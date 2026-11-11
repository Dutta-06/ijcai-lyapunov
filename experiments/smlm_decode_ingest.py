"""
Score the detections that a GPU method (DECODE) returned on the exported frames, with the harness of the other methods
======================================================================================================================

Input, one CSV per set and split (dev, test): columns  frame,x_nm,y_nm,photons,prob  — `frame` the 0-based index WITHIN the exported TIFF (dev.tif / test.tif), positions in nm in the
CORNER convention of the ground truth (the centre of pixel j at (j + 0.5) * 100 nm; DECODE's pixel coordinates have the pixel centre at the integer, so x_nm = (x_px + 0.5) * 100),
`prob` the emitter's detection probability (all candidates down to prob 0.05 or lower), one row per candidate.

The detection threshold on `prob` is chosen on the DEVELOPMENT frames only, twice (by the efficiency E and by the Jaccard index J, as every other method in the study), a constant
position offset (the mean signed error of the matched development pairs) is estimated there and removed, and the TEST frames are only scored. Output: the scores, the paired
comparison is made by `exp_smlm_iss.py` from the saved raw (results/smlm/decode_raw.npz, keys '<set>|<op>|test|off' and '|flat' in the packing of the other raws).

    python experiments/smlm_decode_ingest.py --set MT0.N1.LD --dev dev.csv --test test.csv [--tag DECODE]

A convention check is printed: the mean signed error of matched development pairs; a value near (+-50, +-50) nm means that the pixel convention was not converted.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, HERE)
import smlm_common as C     # noqa: E402
from smlm_decode_export import protocol    # noqa: E402

THRESHOLDS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98)
RAW = os.path.join(HERE, '..', 'results', 'smlm', 'decode_raw.npz')


def read_csv(path):
    raw = np.genfromtxt(path, delimiter=',', names=True)
    need = {'frame', 'x_nm', 'y_nm', 'photons', 'prob'}
    if not need <= set(raw.dtype.names):
        raise ValueError('%s: columns %s, need %s' % (path, raw.dtype.names, sorted(need)))
    return np.atleast_1d(raw)


def per_frame(rows, n_frames, thr):
    """Detections (x_nm, y_nm, photons) per frame with prob >= thr."""
    sel = rows[rows['prob'] >= thr]
    fr = sel['frame'].astype(int)
    if len(fr) and (fr.min() < 0 or fr.max() >= n_frames):
        raise ValueError('frame index outside 0..%d' % (n_frames - 1))
    return [np.c_[sel['x_nm'][fr == i], sel['y_nm'][fr == i], sel['photons'][fr == i]] for i in range(n_frames)]


def shift(dets, off):
    return [np.c_[d[:, 0] - off[0], d[:, 1] - off[1], d[:, 2:]] if len(d) else d for d in dets]


def evaluate(name, dev_rows, test_rows, gt, root=None):
    """Choose the thresholds on the development frames, score the test frames. Returns dict(op -> dict(thr, off, dev, test, test_dets))."""
    (da, db), (ta, tb) = protocol(name)['dev'], protocol(name)['test']
    out = {}
    cand = []
    for thr in THRESHOLDS:
        dets = per_frame(dev_rows, db - da, thr)
        # frames of the development file are 0-based within the file; the ground truth is indexed by the stack frame
        r0 = C.score(dets, gt, np.arange(da, db))
        off = np.nan_to_num(np.array([r0['bias_x_nm'], r0['bias_y_nm']]))
        dc = shift(dets, off)
        r = C.score(dc, gt, np.arange(da, db))
        cand.append((thr, off, r))
    pick = {'E': max(cand, key=lambda c: (round(c[2]['efficiency'], 6), c[0])), 'J': max(cand, key=lambda c: (round(c[2]['jaccard'], 6), -c[2]['rmse_nm']))}
    for op, (thr, off, r) in pick.items():
        tdets = shift(per_frame(test_rows, tb - ta, thr), off)
        out[op] = dict(thr=thr, off=off, dev=r, test=C.score(tdets, gt, np.arange(ta, tb)), test_dets=tdets, uncorrected_bias=off)
    return out


def pack(dets):
    off = np.r_[0, np.cumsum([len(d) for d in dets])].astype(np.int64)
    flat = np.concatenate([d[:, :2] if len(d) else np.zeros((0, 2)) for d in dets])
    flat = np.c_[flat, np.concatenate([d[:, 2] if len(d) else np.zeros(0) for d in dets]), np.zeros(len(flat))]
    return off, flat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', required=True, choices=C.SETS)
    ap.add_argument('--dev', required=True)
    ap.add_argument('--test', required=True)
    ap.add_argument('--tag', default='DECODE')
    ap.add_argument('--save', action='store_true', help='merge the chosen operating points into results/smlm/decode_raw.npz')
    A = ap.parse_args()
    gt = C.load_gt(C.data_path(A.set, 'gt'))
    res = evaluate(A.set, read_csv(A.dev), read_csv(A.test), gt)
    for op, o in res.items():
        print('%s %s %s-tuned: prob >= %g, offset (%.1f, %.1f) nm | dev J %.3f RMSE %.1f E %.1f | test J %.3f P %.3f R %.3f RMSE %.1f E %.1f' % (
            A.tag, A.set, op, o['thr'], o['off'][0], o['off'][1], o['dev']['jaccard'], o['dev']['rmse_nm'], o['dev']['efficiency'], o['test']['jaccard'], o['test']['precision'],
            o['test']['recall'], o['test']['rmse_nm'], o['test']['efficiency']))
        if max(abs(o['uncorrected_bias'])) > 30:
            print('   WARNING: the development pairs have a mean signed error of (%.1f, %.1f) nm before the offset correction: check the pixel convention (+-50 nm = not converted)' % tuple(o['uncorrected_bias']))
    if A.save:
        old = dict(np.load(RAW)) if os.path.exists(RAW) else {}
        for op, o in res.items():
            off, flat = pack(o['test_dets'])
            old['%s|%s-%s|test|off' % (A.set, A.tag, op)] = off
            old['%s|%s-%s|test|flat' % (A.set, A.tag, op)] = flat
        np.savez_compressed(RAW, **old)
        print('saved', RAW)


if __name__ == '__main__':
    main()
