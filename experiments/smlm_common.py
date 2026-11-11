"""
Evaluation harness for 2D single-molecule localization benchmarks, with a Picasso baseline
==========================================================================================

Everything a localization method needs to be scored under identical inputs and metrics on the simulated 2D microtubule training sets
of the 2016 localization-microscopy challenge (64 x 64 px frames, 100 nm pixels, EMCCD noise; ground truth public): the loaders for
the image stack and the ground truth, the matching and the metrics (Jaccard, RMSE, precision, recall, efficiency), a frame bootstrap,
and a Picasso runner / tuner (identification and fit box, Gaussian pre-filter, median padding, fit method, post-filters, a persistent
worker for many runs). Run as a script (``python smlm_common.py --picasso-worker ...``) it is the Picasso worker process.

Conventions (the places where a comparison can silently go wrong)
-----------------------------------------------------------------
* Positions are in nm. x is the image COLUMN axis, y the ROW axis. (0, 0) is the OUTER CORNER of the top-left pixel: the centre of the
  pixel with 0-based indices (row i, column j) is at ((j + 0.5) * 100 nm, (i + 0.5) * 100 nm), and the field is [0, 6400] nm in both
  axes. The ground truth follows this corner convention; a pixel-index coordinate p (pixel centre = integer, the convention of Picasso
  and of most fitting code) is converted with x_nm = (p + 0.5) * 100. Getting this wrong shifts every localization by 50 nm.
* Frame indices are 0-BASED stack indices everywhere in this module (frame 0 = first image of the stack). The ground-truth file numbers its
  frames from 1; ``load_gt`` converts (frame = file frame - 1). The offset was checked against the images (see ``results/smlm/calibration.md``).
* Images are in DN (digital numbers). Camera: baseline 100 DN, EM gain 300, 45 e-/ADU, QE 0.90, so 1 photon = 6.0 DN above the baseline;
  ``to_photons`` returns (DN - 100) / 6 (photons, QE included).
* Matching is a per-frame Hungarian assignment (maximum number of pairs within the tolerance, then minimum total distance); a pair at exactly
  the tolerance distance is a match. Every metric is pooled over the frames given (sum of TP, FP, FN), not averaged per frame.
* Efficiency E = 100 - sqrt((100 - J)^2 + (alpha * RMSE)^2), J the Jaccard index in percent, RMSE the lateral root-mean-square error of the
  matched pairs in nm and alpha = 1 per nm (the lateral value of the challenge's ranking score; E = 100 is perfect).
* The challenge's own assessment additionally dropped ground-truth emitters closer than 450 nm to the field border and kept only the 75 %
  brightest ground-truth activations (applied after the pairing). ``score`` does neither by default; ``border_nm`` and ``photon_min`` switch them on.
* A detection with a non-finite coordinate never matches and is a false positive in EVERY setting (also when the border exclusion is on).
* The RMSE contains the constant offset of a symmetric fit against the simulated emitter (see ``results/smlm/calibration.md``: the offset
  depends on the estimator and grows with |z|; the median of a pixel-integrated Gaussian fit pooled over |z| is about (+11, -9) nm). ``score`` reports ``bias_x_nm``, ``bias_y_nm`` and
  ``rmse_centered_nm`` (RMSE about the mean error of the matched pairs); ``shift_detections`` applies a constant correction estimated elsewhere
  (e.g. on development frames).

Where the rules come from (what is verified and what is from memory)
-------------------------------------------------------------------
Source of the assessment: D. Sage et al., "Super-resolution fight club: assessment of 2D and 3D single-molecule localization microscopy
software", Nature Methods 16:387-395 (2019), Methods; and the pages of the 2016 challenge, now the "Metric" page of the Super-Resolution
Microscopy Hub, https://srm.epfl.ch/srm/assessment/index.html (figure ``metric-2d.png``; https://srm.epfl.ch/srm/metric.html redirects there).
* VERIFIED against that page (read 2026-11-07): lateral tolerance 250 nm for 2D (no axial tolerance); pairing of each localization with a
  ground-truth emitter; Recall = TP / (TP + FN), Precision = TP / (TP + FP), Jaccard = TP / (TP + FP + FN); RMSE_lat = sqrt(mean over the TP of
  (dx^2 + dy^2)), i.e. the 2D Euclidean error pooled over the matched pairs; mean signed errors Delta_x, Delta_y of the pairs.
* FROM MEMORY (moderate confidence, NOT on that page and not checked against the paper, which was not reachable): the efficiency
  E = 100 - sqrt((100 - J)^2 + alpha^2 RMSE^2) with alpha = 1 per nm for the lateral ranking; the exclusion of ground-truth emitters within
  450 nm of the field border; the use of only the 75 % brightest ground-truth activations (here: the quantile of the WHOLE set's activations,
  ``photon_threshold``); that both exclusions are applied after the pairing; the Hungarian algorithm as the pairing method; whether unpaired
  detections in the border zone are dropped (``score`` does so, and counts non-finite ones as false positives). An absolute efficiency
  therefore carries these uncertainties (e.g. a per-axis instead of a 2D RMSE would raise it); differences between methods scored here
  under the same rules do not.
* Consequence of "after the pairing" (documented behaviour, tested): an excluded dim emitter can take (steal) the detection of a kept emitter
  within the tolerance when the detection is closer to the dim one; the pair is then dropped with the dim emitter and the kept emitter counts
  as a false negative although the detection lies within the tolerance of it. This is rare: in the baseline ladder
  (``results/smlm/picasso_baseline.md``) it affects 0-2 of the 379 kept emitters of the LD test frames and 19-60 of the 3314 of the HD test
  frames (0.6-1.8 %); the referee counted about 1 % of the ground-truth rows (``review/REVIEW_roundSM1.md``).

Intervals: ``bootstrap_ci`` and ``bootstrap_diff_ci`` resample circular BLOCKS of ``DEFAULT_BLOCK`` = 10 consecutive frames by default. An
emitter stays on for several consecutive frames (consecutive-frame runs of one fluorophore average 2.27 frames in MT0.N1.LD, 2.26 in
MT0.N1.HD, see ``results/smlm/calibration.md``), so neighbouring frames are not independent; resampling single frames (``block=1``, the i.i.d.
bootstrap, available and to be labelled as such) gives intervals that are too narrow (an independent referee measured a coverage of 61-69 %
instead of 95 % for emitters persisting about 3 frames, 82 % on 39 disjoint windows of the real stack, against 90 % with blocks of 10; see
``review/REVIEW_roundSM1.md`` and ``tests/test_smlm_common.py``). Only the variability of the frame sample is captured, not that of other
datasets or of the simulation, and the development and test frames come from the same three simulated microtubules: ``positions.csv`` of a set (fetched
by ``data/smlm/fetch.py``) lists the exact fluorophore sites of those microtubules (no frame information), so a method that learns the structure or uses
these sites as a prior can exploit what Picasso cannot, and has to say so. The tuning effort on the development frames is part of a comparison: a method
tuned over a handful of parameters is compared against a baseline tuned over its documented options (see ``experiments/exp_smlm_picasso.py``, which
reports how many configurations were scored).

Protocol used by the experiments: frames 0-999 of a set are the DEVELOPMENT frames (tuning, calibration), frames 1000-1499 the TEST frames of
the spike. Nothing is tuned on the test frames.

Dropping in another method (data from ``data/smlm/fetch.py``, in $SMLM_DATA)::

    import numpy as np, smlm_common as C
    frames = C.load_frames(C.data_path('MT0.N1.LD', 'stack'), *C.TEST_FRAMES)       # (500, 64, 64) DN; C.to_photons(frames) for photons
    gt = C.load_gt(C.data_path('MT0.N1.LD', 'gt'))
    dets = my_method(frames)               # list of 500 arrays (n_i, 2+): x_nm, y_nm (corner convention; C.pixels_to_nm converts pixel indices)
    res = C.score(dets, gt, np.arange(*C.TEST_FRAMES))                              # jaccard, rmse_nm, precision, recall, efficiency, bias_*, per_frame
    ci = C.bootstrap_ci(res)                                                        # 95 % block-bootstrap intervals (blocks of 10 frames)
    diff = C.bootstrap_diff_ci(res, other_res)                                      # paired difference of two methods scored on the same frames
"""
import io
import json
import math
import os
import subprocess
import sys
import tempfile
import time
import zipfile
from collections import namedtuple

import numpy as np
from scipy.optimize import linear_sum_assignment

PIXEL_NM = 100.0                  # camera pixel (nm)
FIELD_NM = 6400.0                 # field of view (64 px)
BASELINE_DN = 100.0               # camera baseline (DN)
EM_GAIN = 300.0
E_PER_ADU = 45.0                  # electrons per analog-digital unit
QE = 0.90
DN_PER_PHOTON = QE * EM_GAIN / E_PER_ADU          # 6.0
TOL_NM = 250.0                    # the challenge's lateral matching tolerance
ALPHA_PER_NM = 1.0                # lateral efficiency weight of the challenge
GT_FRAME_BASE = 1                 # the ground-truth file numbers frames from 1
DEV_FRAMES = (0, 1000)            # [start, stop) development frames (tuning, calibration)
TEST_FRAMES = (1000, 1500)        # [start, stop) test frames of the spike
DEFAULT_DATA = '/home/user/smlm_data'
SETS = ('MT0.N1.LD', 'MT0.N1.HD', 'MT0.N2.LD', 'MT0.N2.HD')

GT_DTYPE = np.dtype([('id', 'i8'), ('frame', 'i8'), ('x_nm', 'f8'), ('y_nm', 'f8'), ('z_nm', 'f8'), ('photons', 'f8')])
PER_FRAME_DTYPE = np.dtype([('frame', 'i8'), ('tp', 'i8'), ('fp', 'i8'), ('fn', 'i8'), ('n_det', 'i8'), ('n_gt', 'i8'),
                            ('sse', 'f8'), ('sx', 'f8'), ('sy', 'f8')])
Match = namedtuple('Match', 'tp fp fn pairs')


# ----------------------------------------------------------------------------------------------------------------------------------
# data access
# ----------------------------------------------------------------------------------------------------------------------------------
def data_path(name='MT0.N1.LD', kind='stack', root=None):
    """Path of a downloaded file: kind 'stack' (zip of the TIFF stack), 'gt' (activations.csv, the ground truth), 'positions', 'summary'."""
    root = root or os.environ.get('SMLM_DATA', DEFAULT_DATA)
    files = {'stack': 'sequence-%s-2D-Exp-as-stack.zip' % name, 'gt': 'activations.csv', 'positions': 'positions.csv',
             'summary': 'stats_activation_sumary.csv'}
    return os.path.join(root, name, files[kind])


def _tiff_bytes(path, max_bytes=1 << 30):
    """Bytes of a TIFF file, or of the single TIFF member of a zip (decompressed in memory, nothing is extracted to disk)."""
    path = str(path)
    if path.lower().endswith('.zip'):
        with zipfile.ZipFile(path) as z:
            members = [i for i in z.infolist() if not i.is_dir() and i.filename.lower().endswith(('.tif', '.tiff'))]
            if len(members) != 1:
                raise ValueError('expected exactly one TIFF member in %s, found %d' % (path, len(members)))
            if members[0].file_size > max_bytes:
                raise ValueError('TIFF member of %s is %d bytes (limit %d)' % (path, members[0].file_size, max_bytes))
            return z.read(members[0])
    with open(path, 'rb') as f:
        return f.read()


def load_frames(path, start=0, stop=None):
    """Frames ``start`` .. ``stop - 1`` (0-based, stop exclusive; stop=None: to the end) of a TIFF stack as float64 (n, 64, 64) in DN.

    ``path`` is a multi-page TIFF or a zip holding exactly one TIFF stack (read in memory with Pillow, no extraction to disk; the stack of
    the challenge sets is an ImageJ-style file whose page table sits at the end, hence the whole member is decompressed, ~2 s for 19,996
    frames). Raises IndexError when the stack has fewer than ``stop`` frames (no silent truncation)."""
    from PIL import Image
    if start < 0 or (stop is not None and stop < start):
        raise ValueError('need 0 <= start <= stop, got start=%r stop=%r' % (start, stop))
    im = Image.open(io.BytesIO(_tiff_bytes(path)))
    w, h = im.size
    out, k = [], start
    while stop is None or k < stop:
        try:
            im.seek(k)
        except EOFError:
            if stop is None and out:
                break
            if stop is None and not out:
                raise IndexError('stack has no frame %d' % start)
            raise IndexError('stack has only %d frames, requested stop=%d' % (k, stop))
        a = np.asarray(im)
        if a.ndim != 2:
            raise ValueError('expected single-channel frames, got array of shape %r' % (a.shape,))
        out.append(a.astype(np.float64))
        k += 1
    if not out:
        return np.empty((0, h, w))
    return np.stack(out)


def to_photons(frames):
    """Photons per pixel from DN: (DN - 100) / 6 (baseline 100 DN; total gain QE * EM gain / (e-/ADU) = 6 DN per photon)."""
    return (np.asarray(frames, dtype=np.float64) - BASELINE_DN) / DN_PER_PHOTON


_GT_HEADER = ['ground-truth', 'frame', 'xnano', 'ynano', 'znano', 'intensity']


def load_gt(path, frame_base=GT_FRAME_BASE):
    """Ground truth (``activations.csv``: one row per emitter ON event) as a structured array sorted by (frame, id).

    Fields: id (activation number, 1-based as in the file), frame (0-BASED stack index = file frame - ``frame_base``; the file counts from 1),
    x_nm, y_nm (nm, corner convention: x column axis, y row axis, field [0, 6400]), z_nm (nm, ignored by 2D methods), photons (emitted photons
    of the activation in that frame, incl. partial-frame activations)."""
    with open(path) as f:
        head = [s.strip().lower() for s in f.readline().split(',')]
    if head != _GT_HEADER:
        raise ValueError('unexpected ground-truth header %r (expected %r)' % (head, _GT_HEADER))
    raw = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    gt = np.empty(len(raw), dtype=GT_DTYPE)
    gt['id'] = raw[:, 0].astype(np.int64)
    fr = raw[:, 1]
    if not np.all(fr == np.round(fr)):
        raise ValueError('non-integer frame numbers in %s' % path)
    gt['frame'] = fr.astype(np.int64) - frame_base
    gt['x_nm'], gt['y_nm'], gt['z_nm'], gt['photons'] = raw[:, 2], raw[:, 3], raw[:, 4], raw[:, 5]
    if len(gt) and gt['frame'].min() < 0:
        raise ValueError('negative frame index after the base shift (frame_base=%r)' % frame_base)
    return gt[np.lexsort((gt['id'], gt['frame']))]


def photon_threshold(gt, keep=0.75):
    """Photon count above which the brightest ``keep`` fraction of the dataset's ground-truth activations lie (the challenge kept 0.75).
    Pass the whole dataset's ``gt`` (not a subset of frames) so that the threshold does not depend on the frames evaluated."""
    return float(np.quantile(gt['photons'], 1.0 - keep))


def run_lengths(gt):
    """Persistence of the emitters: the ground truth has one row per activation per frame, and a fluorophore that stays on appears in consecutive frames with
    identical coordinates. Returns dict(runs (number of runs of one fluorophore in consecutive frames), mean (mean run length in frames), lengths (the
    length of every run), ge3 (share of the runs lasting 3 or more frames), dist (number of runs of length 1 .. 8))."""
    key = np.round(gt['x_nm'] * 1000).astype(np.int64) * 10 ** 7 + np.round(gt['y_nm'] * 1000).astype(np.int64)
    order = np.lexsort((gt['frame'], key))
    a, kk = gt[order], key[order]
    new = np.ones(len(a), bool)
    new[1:] = (kk[1:] != kk[:-1]) | (a['frame'][1:] != a['frame'][:-1] + 1)
    lengths = np.bincount(np.cumsum(new) - 1)
    return dict(runs=len(lengths), mean=float(lengths.mean()), lengths=lengths, ge3=float(np.mean(lengths >= 3)), dist=np.bincount(lengths, minlength=9)[1:9])


def pixels_to_nm(xy_px, pixel_nm=PIXEL_NM):
    """Pixel-index coordinates (pixel centre = integer, x along the columns, y along the rows; Picasso and most fitting code) -> nm in the
    ground-truth convention (origin at the outer corner of the first pixel): nm = (p + 0.5) * pixel_nm."""
    return (np.asarray(xy_px, dtype=np.float64) + 0.5) * pixel_nm


def nm_to_pixels(xy_nm, pixel_nm=PIXEL_NM):
    """Inverse of ``pixels_to_nm``: ground-truth nm -> pixel-index coordinates (pixel centre = integer): p = nm / pixel_nm - 0.5."""
    return np.asarray(xy_nm, dtype=np.float64) / pixel_nm - 0.5


# ----------------------------------------------------------------------------------------------------------------------------------
# matching and metrics
# ----------------------------------------------------------------------------------------------------------------------------------
def _xy(a):
    """(n, 2) float64 positions from an (n, 2+) array (extra columns ignored); an empty input is accepted in any shape."""
    a = np.asarray(a, dtype=np.float64)
    if a.size == 0:
        return np.empty((0, 2))
    if a.ndim != 2 or a.shape[1] < 2:
        raise ValueError('expected an (n, 2+) array of positions, got shape %r' % (a.shape,))
    return a[:, :2]


def match(det_xy_nm, gt_xy_nm, tol_nm=TOL_NM):
    """Hungarian pairing of detections with ground-truth emitters of one frame (positions in nm, columns x, y).

    Among all one-to-one pairings of pairs closer than ``tol_nm`` (inclusive) it takes one with the maximum number of pairs and, among those,
    the minimum total distance (scipy ``linear_sum_assignment`` on the distance matrix, out-of-tolerance entries priced above any feasible total).
    A detection with a non-finite coordinate never matches. Returns Match(tp, fp, fn, pairs); ``pairs`` is an (tp, 2) int array of
    (detection index, ground-truth index). A second detection of the same emitter is a false positive; an unpaired emitter a false negative."""
    det, gt = _xy(det_xy_nm), _xy(gt_xy_nm)
    n, m = len(det), len(gt)
    empty = np.empty((0, 2), dtype=np.int64)
    if n == 0 or m == 0:
        return Match(0, n, m, empty)
    with np.errstate(invalid='ignore'):
        d = np.sqrt(((det[:, None, :] - gt[None, :, :]) ** 2).sum(axis=2))
        ok = d <= tol_nm
    if not ok.any():
        return Match(0, n, m, empty)
    big = 2.0 * (min(n, m) + 1) * max(tol_nm, 1e-12) + 1.0
    r, c = linear_sum_assignment(np.where(ok, d, big))
    keep = ok[r, c]
    pairs = np.stack([r[keep], c[keep]], axis=1).astype(np.int64)
    tp = len(pairs)
    return Match(tp, n - tp, m - tp, pairs)


def efficiency(jaccard, rmse_nm, alpha_per_nm=ALPHA_PER_NM):
    """E = 100 - sqrt((100 - J)^2 + (alpha RMSE)^2), J = Jaccard index in percent (``jaccard`` is the fraction, 0..1), RMSE in nm."""
    return 100.0 - math.sqrt((100.0 - 100.0 * jaccard) ** 2 + (alpha_per_nm * rmse_nm) ** 2)


def _div(a, b):
    return a / b if b > 0 else float('nan')


def pooled_metrics(tp, fp, fn, sse, sx=0.0, sy=0.0, alpha_per_nm=ALPHA_PER_NM):
    """Pooled metrics from summed counts: jaccard, precision, recall (fractions), rmse_nm (sqrt(sse / tp)), efficiency, bias_x_nm, bias_y_nm
    (mean signed detection - truth error of the matched pairs), rmse_centered_nm (RMSE about that mean offset, i.e. sqrt(rmse^2 - bias_x^2 - bias_y^2): the
    spread of the matched errors without their constant shift). 0/0 is NaN (e.g. the RMSE and the efficiency when nothing was matched)."""
    tp, fp, fn = float(tp), float(fp), float(fn)
    jac = _div(tp, tp + fp + fn)
    rmse = math.sqrt(sse / tp) if tp > 0 else float('nan')
    eff = efficiency(jac, rmse, alpha_per_nm) if (tp > 0 and not math.isnan(jac)) else float('nan')
    bx, by = _div(sx, tp), _div(sy, tp)
    rmse_c = math.sqrt(max(sse / tp - bx * bx - by * by, 0.0)) if tp > 0 else float('nan')      # RMSE about the mean offset (spread without the bias)
    return dict(jaccard=jac, rmse_nm=rmse, precision=_div(tp, tp + fp), recall=_div(tp, tp + fn), efficiency=eff,
                bias_x_nm=bx, bias_y_nm=by, rmse_centered_nm=rmse_c)


def score(detections, gt, frame_indices, tol_nm=TOL_NM, border_nm=0.0, photon_min=None, alpha_per_nm=ALPHA_PER_NM):
    """Score per-frame detections against the ground truth.

    detections: list with one (n, 2+) array per evaluated frame, columns x_nm, y_nm (nm, corner convention, see the module docstring);
    gt: structured array from ``load_gt`` (any superset of the frames); frame_indices: the 0-based stack index of each detection list
    (unique, same length). Pairing is per frame (``match``), using all detections and all ground-truth emitters of the frame. Optional
    exclusions are applied AFTER the pairing, as in the challenge (from memory, see the module docstring): ``border_nm`` > 0 drops
    ground-truth emitters (and their pairs) closer than that to the field border and also the unpaired detections that lie in that zone;
    ``photon_min`` drops ground-truth activations dimmer than that (see ``photon_threshold``) and their pairs. With the defaults every
    emitter and every detection counts. A detection with a non-finite coordinate is never paired and is a false positive in every
    setting (it does not vanish with the border exclusion). Because the pairing precedes the exclusions, an excluded emitter can take the
    detection of a kept emitter (the "stealing" case: the detection is closer to the dim or border emitter); the kept emitter is then a
    false negative although the detection lies within the tolerance of it.

    Returns a dict: jaccard, precision, recall (fractions), rmse_nm, efficiency (see ``efficiency``), n_det and n_gt (counted in the metrics),
    tp, fp, fn, bias_x_nm, bias_y_nm, rmse_centered_nm (matched pairs, detection - truth), n_det_all, n_gt_all (before the exclusions), and per_frame, a
    structured array (PER_FRAME_DTYPE) for ``bootstrap_ci``. NaN where a ratio is 0/0."""
    frames = np.asarray(frame_indices, dtype=np.int64).ravel()
    if len(detections) != len(frames):
        raise ValueError('%d detection lists for %d frame indices' % (len(detections), len(frames)))
    if len(np.unique(frames)) != len(frames):
        raise ValueError('duplicate frame indices')
    order = np.argsort(gt['frame'], kind='stable')
    g = gt[order]
    lo = np.searchsorted(g['frame'], frames, side='left')
    hi = np.searchsorted(g['frame'], frames, side='right')
    use_border = border_nm is not None and border_nm > 0
    use_photon = photon_min is not None
    if use_photon and 'photons' not in (gt.dtype.names or ()):
        raise ValueError('photon_min needs a ground truth with a photons field')
    per = np.zeros(len(frames), dtype=PER_FRAME_DTYPE)
    per['frame'] = frames
    n_det_all = n_gt_all = 0

    def inside(p):
        return (p[:, 0] >= border_nm) & (p[:, 0] <= FIELD_NM - border_nm) & (p[:, 1] >= border_nm) & (p[:, 1] <= FIELD_NM - border_nm)

    for i in range(len(frames)):
        gf = g[lo[i]:hi[i]]
        gxy = np.stack([gf['x_nm'], gf['y_nm']], axis=1) if len(gf) else np.empty((0, 2))
        dxy = _xy(detections[i])
        n_det_all += len(dxy)
        n_gt_all += len(gxy)
        m = match(dxy, gxy, tol_nm)
        keep_g = np.ones(len(gxy), bool)
        if use_border:
            keep_g &= inside(gxy)
        if use_photon:
            keep_g &= gf['photons'] >= photon_min
        det_paired = np.zeros(len(dxy), bool)
        gt_paired = np.zeros(len(gxy), bool)
        if m.tp:
            det_paired[m.pairs[:, 0]] = True
            gt_paired[m.pairs[:, 1]] = True
        sel = keep_g[m.pairs[:, 1]] if m.tp else np.zeros(0, bool)
        pairs = m.pairs[sel]
        tp = len(pairs)
        fn = int((keep_g & ~gt_paired).sum())
        fp_mask = ~det_paired
        if use_border and len(dxy):
            fp_mask &= inside(dxy) | ~np.isfinite(dxy).all(axis=1)           # a non-finite detection is a false positive, never dropped
        fp = int(fp_mask.sum())
        if tp:
            diff = dxy[pairs[:, 0]] - gxy[pairs[:, 1]]
            per['sse'][i] = float((diff ** 2).sum())
            per['sx'][i], per['sy'][i] = float(diff[:, 0].sum()), float(diff[:, 1].sum())
        per['tp'][i], per['fp'][i], per['fn'][i] = tp, fp, fn
        per['n_det'][i], per['n_gt'][i] = tp + fp, tp + fn
    tot = {k: float(per[k].sum()) for k in ('tp', 'fp', 'fn', 'sse', 'sx', 'sy')}
    res = pooled_metrics(tot['tp'], tot['fp'], tot['fn'], tot['sse'], tot['sx'], tot['sy'], alpha_per_nm)
    res.update(n_det=int(per['n_det'].sum()), n_gt=int(per['n_gt'].sum()), tp=int(tot['tp']), fp=int(tot['fp']), fn=int(tot['fn']),
               n_det_all=int(n_det_all), n_gt_all=int(n_gt_all), per_frame=per)
    return res


_BOOT_KEYS = ('jaccard', 'rmse_nm', 'precision', 'recall', 'efficiency')
BOOT_KEYS_EXT = _BOOT_KEYS + ('rmse_centered_nm', 'bias_x_nm', 'bias_y_nm')          # ``keys`` for intervals of the bias and the centred RMSE too
DEFAULT_BLOCK = 10                  # frames per resampled block: emitters stay on for ~2.3 frames, so neighbouring frames are not independent


def _per_frame_array(x):
    return x['per_frame'] if isinstance(x, dict) else x


def _boot_indices(F, B, seed, block):
    """Yield (first replicate number, index array (n, F)) of bootstrap resamples of F frames: plain draws with replacement (block = 1, the
    i.i.d. bootstrap), or circular blocks of ``block`` consecutive frames (block > 1)."""
    rng = np.random.default_rng(seed)
    block = int(max(1, block))
    nb = -(-F // block)
    chunk = max(1, min(B, 2_000_000 // max(F, 1)))
    for b0 in range(0, B, chunk):
        n = min(chunk, B - b0)
        if block == 1:
            idx = rng.integers(0, F, size=(n, F))
        else:
            starts = rng.integers(0, F, size=(n, nb))
            idx = ((starts[:, :, None] + np.arange(block)[None, None, :]) % F).reshape(n, nb * block)[:, :F]
        yield b0, idx


_SUM_COLS = ('tp', 'fp', 'fn', 'sse', 'sx', 'sy')


def _sum_columns(p):
    return {k: np.asarray(p[k], dtype=np.float64) for k in _SUM_COLS}


def _replicate_metrics(cols, idx, alpha_per_nm):
    tp, fp, fn, sse, sx, sy = (cols[k][idx].sum(axis=1) for k in _SUM_COLS)
    return [pooled_metrics(tp[j], fp[j], fn[j], sse[j], sx[j], sy[j], alpha_per_nm) for j in range(len(idx))]


def _percentile_ci(reps, level, keys=_BOOT_KEYS):
    q = [100.0 * (1.0 - level) / 2.0, 100.0 * (1.0 + level) / 2.0]
    out = {}
    for k in keys:
        v = reps[k][~np.isnan(reps[k])]
        out[k] = (float(np.percentile(v, q[0])), float(np.percentile(v, q[1]))) if len(v) else (float('nan'), float('nan'))
    return out


def bootstrap_ci(per_frame_stats, B=1000, seed=0, block=DEFAULT_BLOCK, level=0.95, alpha_per_nm=ALPHA_PER_NM, keys=_BOOT_KEYS):
    """Percentile confidence intervals of the pooled metrics by resampling frames with replacement, in circular blocks of consecutive frames.

    per_frame_stats: the ``per_frame`` array of ``score`` (or the dict ``score`` returned). Each of the B replicates draws as many frames as
    there are, in circular blocks of ``block`` consecutive frames (default ``DEFAULT_BLOCK`` = 10: an emitter stays on for about 2.3 frames, so
    neighbouring frames are correlated and single-frame resampling, ``block=1`` = the i.i.d. bootstrap, undercovers; label such intervals as
    i.i.d.), and recomputes the pooled metrics. Returns {metric: (lo, hi)} for ``keys`` (default jaccard, rmse_nm, precision, recall,
    efficiency; ``BOOT_KEYS_EXT`` adds rmse_centered_nm, bias_x_nm, bias_y_nm), the central ``level`` interval (percentiles 2.5 / 97.5 for
    0.95); NaN replicates (0/0) are ignored. Only the variability of the frame sample is captured, not that of other datasets or of the simulation."""
    p = _per_frame_array(per_frame_stats)
    F = len(p)
    if F == 0:
        raise ValueError('no frames to resample')
    cols = _sum_columns(p)
    reps = {k: np.empty(B) for k in keys}
    for b0, idx in _boot_indices(F, B, seed, block):
        for j, m in enumerate(_replicate_metrics(cols, idx, alpha_per_nm)):
            for k in keys:
                reps[k][b0 + j] = m[k]
    return _percentile_ci(reps, level, keys)


def bootstrap_diff_ci(per_frame_a, per_frame_b, B=1000, seed=0, block=DEFAULT_BLOCK, level=0.95, alpha_per_nm=ALPHA_PER_NM, keys=_BOOT_KEYS):
    """Paired bootstrap interval of the metric differences A - B of two methods scored on the SAME frames (the ``per_frame`` arrays of two
    ``score`` calls with identical ``frame_indices``): both are evaluated on every resample of the frames (circular blocks of ``block``
    frames, default ``DEFAULT_BLOCK``), so the frame-to-frame variability that the methods share cancels and the interval is much tighter than
    the difference of the two separate intervals. Returns {metric: (lo, hi)}; an interval that excludes 0 is a difference larger than the
    frame-sampling noise (not covered: other datasets, the simulation, and the shared structure of development and test frames)."""
    a, b = _per_frame_array(per_frame_a), _per_frame_array(per_frame_b)
    if len(a) == 0 or len(a) != len(b) or not np.array_equal(a['frame'], b['frame']):
        raise ValueError('the two methods must be scored on the same (non-empty) list of frames')
    ca, cb = _sum_columns(a), _sum_columns(b)
    reps = {k: np.empty(B) for k in keys}
    for b0, idx in _boot_indices(len(a), B, seed, block):
        for j, (ma, mb) in enumerate(zip(_replicate_metrics(ca, idx, alpha_per_nm), _replicate_metrics(cb, idx, alpha_per_nm))):
            for k in keys:
                reps[k][b0 + j] = ma[k] - mb[k]
    return _percentile_ci(reps, level, keys)


def bootstrap_pooled_ci(columns, stat, B=1000, seed=0, block=DEFAULT_BLOCK, level=0.95):
    """Interval of any pooled statistic of per-frame columns, by the same resampling as ``bootstrap_ci`` (the same seed gives the same frame
    resamples, so intervals of different statistics, and of different methods, are paired).

    columns: {name: (F,) array of per-frame sums}; stat: function of {name: (n,) array of the resample sums} returning (n,) values (e.g.
    ``lambda s: np.sqrt(s['sse'] / s['tp'])``, the RMSE over a subset of the true positives). Returns (lo, hi) (NaN replicates ignored)."""
    cols = {k: np.asarray(v, dtype=np.float64) for k, v in columns.items()}
    F = len(next(iter(cols.values())))
    if F == 0:
        raise ValueError('no frames to resample')
    reps = np.empty(B)
    for b0, idx in _boot_indices(F, B, seed, block):
        with np.errstate(invalid='ignore', divide='ignore'):
            reps[b0:b0 + len(idx)] = stat({k: c[idx].sum(axis=1) for k, c in cols.items()})
    return _percentile_ci({'v': reps}, level, ('v',))['v']


def bootstrap_pooled_diff_ci(columns_a, columns_b, stat, B=1000, seed=0, block=DEFAULT_BLOCK, level=0.95):
    """Paired interval of stat(A) - stat(B) for two sets of per-frame columns over the same F frames (see ``bootstrap_pooled_ci``)."""
    ca = {k: np.asarray(v, dtype=np.float64) for k, v in columns_a.items()}
    cb = {k: np.asarray(v, dtype=np.float64) for k, v in columns_b.items()}
    F = len(next(iter(ca.values())))
    if F == 0 or any(len(v) != F for v in list(ca.values()) + list(cb.values())):
        raise ValueError('the two sets of columns must cover the same (non-empty) list of frames')
    reps = np.empty(B)
    for b0, idx in _boot_indices(F, B, seed, block):
        with np.errstate(invalid='ignore', divide='ignore'):
            reps[b0:b0 + len(idx)] = stat({k: c[idx].sum(axis=1) for k, c in ca.items()}) - stat({k: c[idx].sum(axis=1) for k, c in cb.items()})
    return _percentile_ci({'v': reps}, level, ('v',))['v']


# ----------------------------------------------------------------------------------------------------------------------------------
# Picasso
# ----------------------------------------------------------------------------------------------------------------------------------
PICASSO_DEFAULTS = dict(
    min_net_gradient=5000.0,    # Picasso's default (graphical and command-line)
    box=7,                      # Picasso's default box side length (px, odd): identification, and the fit box unless fit_box is given
    fit='gausslq',              # Picasso's default 2D fit: least-squares elliptical Gaussian (Levenberg-Marquardt)
    baseline=BASELINE_DN, gain=EM_GAIN, sensitivity=E_PER_ADU, qe=QE, pixel_nm=PIXEL_NM,
    eps=None, max_it=None,      # None: Picasso's own convergence schedule
    # ---- options of Picasso beyond its default pipeline (the rungs of a baseline ladder); every one is OFF at its default value, which gives the
    # ---- default pipeline bit for bit
    fit_box=None,               # side of the box fitted around each identified spot (px, odd); None: the identification box ``box``
    gauss_sigma=None,           # Picasso's "Gaussian Filter Sigma" (px): every frame is smoothed before identification (the fit uses the raw frame); None / 0: off
    pad=0,                      # median padding: every frame is padded by ``pad`` px on each side with its own median, Picasso identifies in the padded frame
                                # (it cannot identify in a border zone of box // 2 px at the top / left and box // 2 + 1 at the bottom / right), positions are
                                # shifted back and spots whose fitted centre lies outside the field are dropped
    # ---- post-filters on the returned localizations (rows of the table; None: off; strict inequalities; the width window needs BOTH sx and sy inside)
    min_photons=None,           # photons (QE-corrected, as returned) above this
    min_sigma_nm=None,          # |sx| and |sy| (fitted Gaussian widths, nm) above this
    max_sigma_nm=None,          # |sx| and |sy| below this
)
DET_COLUMNS = ('x_nm', 'y_nm', 'photons', 'sx_nm', 'sy_nm', 'net_gradient')    # columns of the arrays run_picasso returns
NG_COL = DET_COLUMNS.index('net_gradient')
PH_COL, SX_COL, SY_COL = DET_COLUMNS.index('photons'), DET_COLUMNS.index('sx_nm'), DET_COLUMNS.index('sy_nm')
POST_FILTER_KEYS = ('min_net_gradient', 'min_photons', 'min_sigma_nm', 'max_sigma_nm')      # the filters of ``filter_detections`` (None / absent: off)
_TABLE_COLUMNS = ('frame', 'x', 'y', 'photons', 'sx', 'sy', 'bg', 'lpx', 'lpy', 'net_gradient')


def detections_from_table(tab, n_frames, pixel_nm=PIXEL_NM, qe=QE):
    """Per-frame (n, 6) arrays (DET_COLUMNS) from a Picasso localization table (dict of arrays frame, x, y, photons, sx, sy, net_gradient).

    Picasso's x, y are pixel-index coordinates (pixel centre = integer; x along the columns): nm = (p + 0.5) * pixel_nm, the corner convention
    of the ground truth. Its 'photons' are photoelectrons (QE is not divided out since Picasso 0.6) and are divided by ``qe`` here to be photons;
    sx, sy (Gaussian widths) are converted from px to nm."""
    frame = np.asarray(tab['frame']).astype(np.int64)
    cols = np.stack([(np.asarray(tab['x'], np.float64) + 0.5) * pixel_nm, (np.asarray(tab['y'], np.float64) + 0.5) * pixel_nm,
                     np.asarray(tab['photons'], np.float64) / qe, np.asarray(tab['sx'], np.float64) * pixel_nm,
                     np.asarray(tab['sy'], np.float64) * pixel_nm, np.asarray(tab['net_gradient'], np.float64)], axis=1) if len(frame) \
        else np.empty((0, len(DET_COLUMNS)))
    order = np.argsort(frame, kind='stable')
    frame, cols = frame[order], cols[order]
    if len(frame) and (frame.min() < 0 or frame.max() >= n_frames):
        raise ValueError('localization frame index outside 0..%d' % (n_frames - 1))
    lo = np.searchsorted(frame, np.arange(n_frames), side='left')
    hi = np.searchsorted(frame, np.arange(n_frames), side='right')
    return [cols[a:b] for a, b in zip(lo, hi)]


def filter_detections(dets, filt):
    """Rows of per-frame (n, 6) detection arrays (DET_COLUMNS, e.g. from ``run_picasso``) that pass the post-filters in ``filt``.

    filt: dict (or the Picasso parameter dict itself) with any of ``POST_FILTER_KEYS``; None or a missing key switches that filter off.
    min_net_gradient: net gradient > value; min_photons: photons > value; min_sigma_nm / max_sigma_nm: |sx| and |sy| both above / below the
    value (a NaN width fails a window). All filters are strict inequalities and are combined with AND. Returns a list of arrays (an array
    that passes entirely is returned as it is, not copied)."""
    ng, ph, lo, hi = (filt.get(k) for k in POST_FILTER_KEYS)
    if ng is None and ph is None and lo is None and hi is None:
        return list(dets)
    out = []
    for d in dets:
        if len(d) == 0:
            out.append(d)
            continue
        keep = np.ones(len(d), bool)
        if ng is not None:
            keep &= d[:, NG_COL] > ng
        if ph is not None:
            keep &= d[:, PH_COL] > ph
        if lo is not None or hi is not None:
            sx, sy = np.abs(d[:, SX_COL]), np.abs(d[:, SY_COL])
            if lo is not None:
                keep &= (sx > lo) & (sy > lo)
            if hi is not None:
                keep &= (sx < hi) & (sy < hi)
        out.append(d if keep.all() else d[keep])
    return out


def median_pad(frames, pad):
    """(n, h, w) frames padded by ``pad`` px on every side with each frame's own median (computed in float, cast to the dtype of the frames, i.e.
    truncated for integer frames). A constant border keeps Picasso's local-maximum search and its Gaussian pre-filter free of edge effects."""
    frames = np.asarray(frames)
    pad = int(pad)
    if pad <= 0:
        return frames
    n, h, w = frames.shape
    out = np.empty((n, h + 2 * pad, w + 2 * pad), dtype=frames.dtype)
    out[:] = np.median(frames.reshape(n, -1), axis=1).astype(frames.dtype)[:, None, None]
    out[:, pad:-pad, pad:-pad] = frames
    return out


def _merge_params(params):
    p = dict(PICASSO_DEFAULTS)
    p.update(params or {})
    unknown = set(p) - set(PICASSO_DEFAULTS)
    if unknown:
        raise ValueError('unknown Picasso parameters %s' % sorted(unknown))
    for k in ('box', 'fit_box'):
        v = p[k]
        if v is None and k == 'fit_box':
            continue
        if v is None or int(v) != v or int(v) < 3 or int(v) % 2 == 0:
            raise ValueError('%s must be an odd integer >= 3 (px), got %r' % (k, v))
    if int(p['pad']) != p['pad'] or p['pad'] < 0:
        raise ValueError('pad must be a non-negative integer (px), got %r' % (p['pad'],))
    if p['gauss_sigma'] is not None and not p['gauss_sigma'] >= 0:
        raise ValueError('gauss_sigma must be None or >= 0 (px), got %r' % (p['gauss_sigma'],))
    return p


def _worker_frames(frames):
    """The stack as Picasso gets it: uint16 when the values are camera integers (as a TIFF is loaded by Picasso), else float32."""
    frames = np.asarray(frames)
    if np.all(frames == np.round(frames)) and frames.min() >= 0 and frames.max() <= 65535:
        return np.round(frames).astype(np.uint16)
    return frames.astype(np.float32)


def _picasso_localize(frames, p, PL):
    """Picasso's identification and fit on an in-memory stack (n, h, w) with the (merged) parameters ``p``; needs ``PL`` = picasso.localize.
    Returns the localization table as a dict of arrays (_TABLE_COLUMNS), positions in pixel-index coordinates of ``frames``.

    Identification (``PL.identify``: local maxima over a box of ``box`` px whose net gradient exceeds min_net_gradient, optionally on frames smoothed
    with Picasso's Gaussian filter) and fit (``PL.fit``: ``fit_box`` px around each spot, ``fit`` method) are called separately, which is what
    ``PL.localize_frames`` does with one box for both; with fit_box = box, no padding and no Gaussian filter the result is the same table."""
    box = int(p['box'])
    fit_box = int(p['fit_box'] or box)
    pad = int(p['pad'])
    h, w = frames.shape[1:]
    work = median_pad(frames, pad)
    movie = PL._InMemoryMovie(work, [])
    ids, _ = PL.identify(movie, p['min_net_gradient'], box, threaded=False, gaussian_filter_sigma=(p['gauss_sigma'] or None))
    hb = fit_box // 2
    ids = ids[(ids['x'] >= hb) & (ids['x'] < work.shape[2] - hb) & (ids['y'] >= hb) & (ids['y'] < work.shape[1] - hb)].reset_index(drop=True)
    if len(ids) == 0:
        return {k: np.zeros(0, dtype=np.uint32 if k == 'frame' else np.float32) for k in _TABLE_COLUMNS}
    cam = {'Baseline': p['baseline'], 'Sensitivity': p['sensitivity'], 'Gain': p['gain'], 'Qe': p['qe'], 'Pixelsize': p['pixel_nm']}
    locs, _ = PL.fit(movie=movie, camera_info=cam, identifications=ids, box=fit_box, fitting_method=p['fit'], eps=p['eps'], max_it=p['max_it'],
                     multiprocess=False)
    tab = {k: np.asarray(locs[k]) for k in _TABLE_COLUMNS}
    if pad:
        x, y = tab['x'].astype(np.float64) - pad, tab['y'].astype(np.float64) - pad
        keep = (x >= -0.5) & (x <= w - 0.5) & (y >= -0.5) & (y <= h - 0.5)                  # a spot centre outside the field is not an emitter of the field
        tab = {k: v[keep] for k, v in tab.items()}
        tab['x'], tab['y'] = x[keep], y[keep]
    return tab


def _picasso_run_and_save(frames, params, out_path, warmup_frames, picasso, PL):
    import warnings
    p = _merge_params(params)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t0 = time.perf_counter()
        _picasso_localize(frames[:max(1, min(warmup_frames, len(frames)))], p, PL)           # JIT compilation, not timed below
        warm = time.perf_counter() - t0
        t0 = time.perf_counter()
        tab = _picasso_localize(frames, p, PL)
        run = time.perf_counter() - t0
    flags = PL.parse_gauss_code(p['fit'])
    tol, maxit = PL.gauss_schedule(flags['mle'], flags['use_gpu'], p['eps'], p['max_it']) if flags else (float('nan'), -1)
    np.savez(out_path, n_frames=len(frames), run_s=run, warmup_s=warm, picasso_version=getattr(picasso, '__version__', 'unknown'),
             fit_tolerance=tol, fit_max_iterations=maxit, **tab)


def _picasso_worker(frames_path, out_path, params, warmup_frames=8):
    """Run Picasso's identification + fit on a frame stack saved as .npy and save the localization table (.npz). Needs picasso."""
    import picasso
    from picasso import localize as PL
    _picasso_run_and_save(_worker_frames(np.load(frames_path)), params, out_path, warmup_frames, picasso, PL)


_SERVER_TAG = '@@smlm-reply@@ '


def _picasso_server(frames_path):
    """A worker that stays alive: it loads the stack once, then answers one JSON request per stdin line ({"out": path, "params": {...}} or
    {"quit": true}) with one tagged JSON line on stdout; everything else Picasso or numba print is ignored by the client."""
    import picasso
    from picasso import localize as PL
    frames = _worker_frames(np.load(frames_path))
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        req = json.loads(line)
        if req.get('quit'):
            break
        try:
            _picasso_run_and_save(frames, req['params'], req['out'], int(req.get('warmup_frames', 8)), picasso, PL)
            reply = {'ok': True}
        except Exception:
            import traceback
            reply = {'ok': False, 'error': traceback.format_exc()[-3000:]}
        sys.stdout.write(_SERVER_TAG + json.dumps(reply) + '\n')
        sys.stdout.flush()


def _picasso_importable():
    try:
        import picasso  # noqa: F401
        return True
    except Exception:
        return False


def _json_default(o):
    return o.item() if hasattr(o, 'item') else str(o)


def _read_locs(path, n_frames, p, info=None):
    """Per-frame detection arrays (DET_COLUMNS), post-filtered with the parameters in ``p``, from a worker's .npz."""
    z = np.load(path, allow_pickle=False)
    tab = {k: z[k] for k in ('frame', 'x', 'y', 'photons', 'sx', 'sy', 'net_gradient')}
    if info is not None:
        info.update(run_s=float(z['run_s']), warmup_s=float(z['warmup_s']), picasso_version=str(z['picasso_version']),
                    n_loc=int(len(tab['frame'])), n_frames=int(z['n_frames']), fit_tolerance=float(z['fit_tolerance']),
                    fit_max_iterations=int(z['fit_max_iterations']))
    dets = detections_from_table(tab, n_frames, p['pixel_nm'], p['qe'])
    return filter_detections(dets, dict(p, min_net_gradient=None))                  # the net-gradient threshold was applied by the identification


def _worker_env():
    env = dict(os.environ)
    for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        env.setdefault(k, '1')
    return env


def run_picasso(frames_dn, params=None, python=None, info=None):
    """Run Picasso (identification by minimum net gradient and box size, then its least-squares Gaussian fit) on frames in DN.

    frames_dn: (n, 64, 64) array in DN (e.g. ``load_frames``); params: dict overriding ``PICASSO_DEFAULTS``: min_net_gradient, box (identification;
    also the fit box unless fit_box is given), fit_box, fit ('gausslq', 'gaussmle', ...), gauss_sigma (Picasso's Gaussian filter, identification
    only), pad (median padding, px), the post-filters min_photons / min_sigma_nm / max_sigma_nm, the camera baseline / gain / sensitivity (e-/ADU) /
    qe, pixel_nm, eps, max_it; the camera parameters are passed to Picasso as they are (baseline 100, EM gain 300, 45 e-/ADU, QE 0.9; Picasso uses
    them for photon counts and precision, the identification works on the raw DN). The defaults of every option reproduce Picasso's default pipeline.
    python: interpreter that can import picasso (default: the environment variable SMLM_PICASSO_PYTHON, else the current interpreter if it
    can import picasso); the work then runs in a subprocess of it (single-threaded) through this very file. info: optional dict filled with
    run_s (identification + fit seconds of the timed run, JIT compilation excluded), warmup_s, picasso_version, n_loc (rows before the
    post-filters), fit_tolerance and fit_max_iterations (Picasso's convergence schedule of the fit method).
    Returns a list of n arrays, one per frame (0-based position in ``frames_dn``), (m, 6) with columns DET_COLUMNS: x_nm, y_nm (nm, corner
    convention), photons (QE-corrected), sx_nm, sy_nm (fitted Gaussian widths), net_gradient. Nothing is filtered after Picasso except by the
    post-filter options. For many runs on the same frames use ``PicassoSession`` (one worker process, no repeated start-up)."""
    frames = np.asarray(frames_dn)
    if frames.ndim != 3:
        raise ValueError('frames must be (n, height, width)')
    p = _merge_params(params)
    python = python or os.environ.get('SMLM_PICASSO_PYTHON')
    with tempfile.TemporaryDirectory(prefix='smlm_picasso_') as tmp:
        fin, fout = os.path.join(tmp, 'frames.npy'), os.path.join(tmp, 'locs.npz')
        np.save(fin, frames)
        if python is None and _picasso_importable():
            _picasso_worker(fin, fout, p)
        elif python is None:
            raise RuntimeError('picasso cannot be imported here; set SMLM_PICASSO_PYTHON to an interpreter that has it')
        else:
            cmd = [python, '-I', os.path.abspath(__file__), '--picasso-worker', fin, fout, json.dumps(p, default=_json_default)]
            r = subprocess.run(cmd, env=_worker_env(), capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError('Picasso worker failed (%s):\n%s' % (r.returncode, (r.stderr or r.stdout)[-3000:]))
        return _read_locs(fout, len(frames), p, info)


class PicassoSession:
    """A Picasso worker that stays alive over many runs on the SAME frames (no start-up and no JIT compilation per run).

    ``with PicassoSession(frames_dn, python) as s: dets = s.run(params, info)``; ``run`` takes what ``run_picasso`` takes (without the frames)
    and returns the same list of per-frame arrays. Without ``python`` (and without $SMLM_PICASSO_PYTHON) the work is done in this process when
    picasso can be imported. The worker is a child process of this one, ended by ``close`` (it is the only process this class ever signals)."""

    def __init__(self, frames_dn, python=None):
        frames = np.asarray(frames_dn)
        if frames.ndim != 3:
            raise ValueError('frames must be (n, height, width)')
        self.n_frames = len(frames)
        self._tmp = tempfile.TemporaryDirectory(prefix='smlm_picasso_')
        self._k = 0
        self._proc = None
        self._frames = None
        python = python or os.environ.get('SMLM_PICASSO_PYTHON')
        try:
            if python is None and _picasso_importable():
                import picasso
                from picasso import localize as PL
                self._frames, self._picasso, self._PL = _worker_frames(frames), picasso, PL
            elif python is None:
                raise RuntimeError('picasso cannot be imported here; set SMLM_PICASSO_PYTHON to an interpreter that has it')
            else:
                fin = os.path.join(self._tmp.name, 'frames.npy')
                np.save(fin, frames)
                self._proc = subprocess.Popen([python, '-I', os.path.abspath(__file__), '--picasso-server', fin], env=_worker_env(), stdin=subprocess.PIPE,
                                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        except Exception:
            self._tmp.cleanup()
            raise

    def run(self, params=None, info=None):
        if self._proc is None and self._frames is None:
            raise RuntimeError('the Picasso session is closed')
        p = _merge_params(params)
        self._k += 1
        out = os.path.join(self._tmp.name, 'locs_%d.npz' % self._k)
        if self._proc is None:
            _picasso_run_and_save(self._frames, p, out, 8, self._picasso, self._PL)
        else:
            self._proc.stdin.write(json.dumps({'out': out, 'params': p}, default=_json_default) + '\n')
            self._proc.stdin.flush()
            tail = []
            while True:
                line = self._proc.stdout.readline()
                if not line:
                    raise RuntimeError('Picasso worker ended unexpectedly (exit status %s):\n%s' % (self._proc.poll(), ''.join(tail)[-3000:]))
                if line.startswith(_SERVER_TAG):
                    reply = json.loads(line[len(_SERVER_TAG):])
                    break
                tail.append(line)
                tail = tail[-40:]
            if not reply.get('ok'):
                raise RuntimeError('Picasso worker failed:\n%s' % reply.get('error', ''))
        try:
            return _read_locs(out, self.n_frames, p, info)
        finally:
            try:
                os.remove(out)
            except OSError:
                pass

    def close(self):
        proc, self._proc = self._proc, None
        if proc is not None:
            try:
                proc.stdin.write(json.dumps({'quit': True}) + '\n')
                proc.stdin.flush()
                proc.stdin.close()
            except (OSError, ValueError):
                pass
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.terminate()                                    # our own child, by its handle
                proc.wait(timeout=30)
            if proc.stdout is not None:
                proc.stdout.close()
        if getattr(self, '_tmp', None) is not None:
            self._tmp.cleanup()
            self._tmp = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


def tune_picasso(dev_frames, dev_gt, grid, params=None, frame_indices=None, objective='efficiency', tol_nm=TOL_NM, border_nm=0.0,
                 python=None, verify_frames=0):
    """Choose Picasso's minimum net gradient on development frames (the rung "net gradient tuned" of a baseline ladder).

    dev_frames: (n, 64, 64) DN; dev_gt: ground truth (``load_gt``) whose 0-based frame field refers to the stack; frame_indices: stack index of
    each dev frame (default 0..n-1); grid: at most 20 candidate thresholds; params: any other Picasso option (structure and post-filters of
    ``run_picasso``) held fixed while the threshold is tuned. Identification keeps the local maxima whose net gradient exceeds the
    threshold and every kept spot is fitted independently, so ONE run at the smallest grid value contains the output of every larger threshold
    (filter ``net_gradient > t``); ``verify_frames`` > 0 re-runs the first that many frames from scratch at the best and at the largest
    threshold and raises if the filtered output differs. ``objective`` is the metric maximised ('efficiency' (default), 'jaccard', ...).
    Returns (best_threshold, curve), curve a list of dicts (min_net_gradient, jaccard, rmse_nm, precision, recall, efficiency, n_det, tp, fp, fn)
    in grid order. Never pass test frames here. (For the other options see ``tune_postfilters`` and the ladder in exp_smlm_picasso.py.)"""
    grid = sorted({float(t) for t in grid})
    if not 1 <= len(grid) <= 20:
        raise ValueError('grid must hold between 1 and 20 distinct values (got %d)' % len(grid))
    n = len(dev_frames)
    frame_indices = np.arange(n) if frame_indices is None else np.asarray(frame_indices)
    p = _merge_params(params)
    p['min_net_gradient'] = grid[0]
    full = run_picasso(dev_frames, p, python)
    curve = []
    for t in grid:
        dets = [d[d[:, NG_COL] > t] for d in full]
        r = score(dets, dev_gt, frame_indices, tol_nm=tol_nm, border_nm=border_nm)
        row = {'min_net_gradient': t}
        row.update({k: r[k] for k in ('jaccard', 'rmse_nm', 'precision', 'recall', 'efficiency', 'n_det', 'tp', 'fp', 'fn')})
        curve.append(row)
    vals = np.array([c[objective] if not math.isnan(c[objective]) else -np.inf for c in curve])
    best_i = int(np.argmax(vals))
    best = curve[best_i]['min_net_gradient']
    if verify_frames:
        m = min(int(verify_frames), n)
        for t in sorted({best, grid[-1]}):
            q = dict(p, min_net_gradient=t)
            direct = run_picasso(dev_frames[:m], q, python)
            filt = [d[d[:, NG_COL] > t] for d in full[:m]]
            for a, b in zip(direct, filt):
                if a.shape != b.shape or not np.allclose(a, b, equal_nan=True):
                    raise RuntimeError('filtering one low-threshold run differs from a run at threshold %g' % t)
    return best, curve


_TUNE_KEYS = ('jaccard', 'precision', 'recall', 'rmse_nm', 'efficiency', 'rmse_centered_nm', 'bias_x_nm', 'bias_y_nm', 'tp', 'fp', 'fn', 'n_det_all')


def tune_postfilters(candidates, gt, frame_indices, grids, start=None, objective='efficiency', max_rounds=3, cache=None, tol_nm=TOL_NM,
                     border_nm=0.0, photon_min=None):
    """Coordinate ascent of post-filter thresholds on the development frames of one Picasso configuration.

    candidates: per-frame (n, 6) arrays (DET_COLUMNS) of ONE low-threshold ``run_picasso`` / ``PicassoSession.run`` (identification keeps the local
    maxima above its threshold and every spot is fitted independently, so a run at the smallest threshold contains the output of every larger
    one; ``filter_detections`` of this table equals a direct run with the filter, see ``tune_picasso``). grids: ordered dict {filter: candidate
    values, None = off}, keys from ``POST_FILTER_KEYS`` in the order of the line searches (e.g. min_net_gradient, min_photons, min_sigma_nm,
    max_sigma_nm). start: initial filter dict (default: every filter off). Each round searches every filter's grid in turn (the others held at their
    current values) and moves to a value only if the objective improves strictly (ties keep the current value and, within a grid, the earlier
    value); the ascent stops after ``max_rounds`` rounds or when a round changes nothing. objective: 'efficiency' (default), 'jaccard', ...;
    NaN counts as -inf. cache: optional dict shared between calls on the SAME candidates (e.g. the efficiency and the Jaccard ascents), which
    then never scores a filter combination twice. Nothing but ``frame_indices`` frames of ``gt`` is looked at: pass development frames only.
    Returns dict(filters (the final filter dict, every key of ``grids``), value, metrics (dict of the final configuration: jaccard, precision,
    recall, rmse_nm, efficiency, rmse_centered_nm, bias_x_nm, bias_y_nm, tp, fp, fn, n_det_all), rounds, converged (the last round changed
    nothing), first (dict(filters, value, metrics) after
    the first line search, i.e. with only the first filter tuned), n_eval (distinct filter combinations in ``cache``), cache)."""
    keys = list(grids)
    bad = set(keys) - set(POST_FILTER_KEYS)
    if bad:
        raise ValueError('unknown filters %s' % sorted(bad))
    cache = {} if cache is None else cache
    cur = {k: None for k in keys}
    cur.update(start or {})

    def ev(f):
        key = tuple(f.get(k) for k in keys)
        if key not in cache:
            r = score(filter_detections(candidates, f), gt, frame_indices, tol_nm=tol_nm, border_nm=border_nm, photon_min=photon_min)
            cache[key] = {k: r[k] for k in _TUNE_KEYS}
        return cache[key]

    def val(m):
        v = m[objective]
        return -np.inf if v is None or (isinstance(v, float) and math.isnan(v)) else float(v)

    best = val(ev(cur))
    first = None
    rounds = 0
    changed = True
    for rnd in range(int(max_rounds)):
        rounds = rnd + 1
        changed = False
        for ki, k in enumerate(keys):
            for v in grids[k]:
                if v == cur[k]:
                    continue
                q = dict(cur)
                q[k] = v
                x = val(ev(q))
                if x > best + 1e-12:
                    best, cur, changed = x, q, True
            if rnd == 0 and ki == 0:
                first = dict(filters=dict(cur), value=best, metrics=dict(ev(cur)))
        if not changed:
            break
    return dict(filters=dict(cur), value=best, metrics=dict(ev(cur)), rounds=rounds, converged=not changed, first=first, n_eval=len(cache), cache=cache)


# ----------------------------------------------------------------------------------------------------------------------------------
# position offset, matched errors, detections on disk
# ----------------------------------------------------------------------------------------------------------------------------------
def shift_detections(dets, dx_nm, dy_nm):
    """Copies of per-frame detection arrays with the constant offset (dx, dy) in nm SUBTRACTED from x_nm, y_nm (columns 0, 1): the correction of
    the mean signed error (detection - truth) of a symmetric fit, as estimated on development frames (``score`` -> bias_x_nm, bias_y_nm).
    Other columns are kept."""
    out = []
    for d in dets:
        d = np.array(d, dtype=np.float64, copy=True)
        if d.ndim == 2 and d.shape[0] and d.shape[1] >= 2:
            d[:, 0] -= dx_nm
            d[:, 1] -= dy_nm
        out.append(d)
    return out


MATCH_DTYPE = np.dtype([('gt_row', 'i8'), ('frame', 'i8'), ('pos', 'i8'), ('matched', '?'), ('det', 'i8'), ('dx', 'f8'), ('dy', 'f8')])


def matched_errors(detections, gt, frame_indices, tol_nm=TOL_NM):
    """One record per ground-truth row of the evaluated frames: whether it is matched (``match``, per frame) and the error (detection - truth,
    nm; NaN when unmatched). Fields: gt_row (index into ``gt``), frame (stack index), pos (position of the frame in ``frame_indices``),
    matched, det (row of the matched detection in ``detections[pos]``, -1 if unmatched), dx, dy. Rows are ordered by frame position, then by the order of ``gt`` within the frame (stable), the same for every call with the
    same ``gt`` and ``frame_indices``, so the true positives of several methods can be intersected row by row (e.g. the RMSE on the true positives
    ALL of them find, free of the survivor bias of a stricter method). Pairing without any exclusion, as ``score`` with its defaults."""
    frames = np.asarray(frame_indices, dtype=np.int64).ravel()
    if len(detections) != len(frames):
        raise ValueError('%d detection lists for %d frame indices' % (len(detections), len(frames)))
    if len(np.unique(frames)) != len(frames):
        raise ValueError('duplicate frame indices')
    order = np.argsort(gt['frame'], kind='stable')
    g = gt[order]
    lo = np.searchsorted(g['frame'], frames, side='left')
    hi = np.searchsorted(g['frame'], frames, side='right')
    parts = []
    for i in range(len(frames)):
        gf = g[lo[i]:hi[i]]
        n = len(gf)
        rec = np.zeros(n, dtype=MATCH_DTYPE)
        rec['gt_row'], rec['frame'], rec['pos'] = order[lo[i]:hi[i]], frames[i], i
        rec['dx'] = rec['dy'] = np.nan
        rec['det'] = -1
        if n:
            gxy = np.stack([gf['x_nm'], gf['y_nm']], axis=1)
            dxy = _xy(detections[i])
            m = match(dxy, gxy, tol_nm)
            if m.tp:
                diff = dxy[m.pairs[:, 0]] - gxy[m.pairs[:, 1]]
                rec['matched'][m.pairs[:, 1]] = True
                rec['det'][m.pairs[:, 1]] = m.pairs[:, 0]
                rec['dx'][m.pairs[:, 1]], rec['dy'][m.pairs[:, 1]] = diff[:, 0], diff[:, 1]
        parts.append(rec)
    return np.concatenate(parts) if parts else np.zeros(0, dtype=MATCH_DTYPE)


def pack_detections(dets, columns=(0, 1, 2, 3)):
    """Per-frame detection arrays -> (off, flat): off (F + 1,) int64 offsets and flat (M, len(columns)) float64 rows (default columns x_nm, y_nm,
    photons, sx_nm), the rows of frame i being flat[off[i]:off[i + 1]]. The compact form stored in results/smlm/picasso_ladder_raw.npz."""
    cols = list(columns)
    blocks = [np.asarray(d, dtype=np.float64)[:, cols] if np.size(d) else np.empty((0, len(cols))) for d in dets]
    off = np.concatenate([[0], np.cumsum([len(b) for b in blocks])]).astype(np.int64)
    flat = np.concatenate(blocks) if blocks else np.empty((0, len(cols)))
    return off, flat.reshape(-1, len(cols))


def unpack_detections(off, flat):
    """Inverse of ``pack_detections``: the list of F per-frame arrays (the columns kept, x_nm and y_nm first, so they can be given to ``score``)."""
    off = np.asarray(off)
    return [np.asarray(flat)[off[i]:off[i + 1]] for i in range(len(off) - 1)]


def load_ladder_detections(path, set_name, row, split):
    """Detections of one row of the baseline ladder from ``results/smlm/picasso_ladder_raw.npz``: ``row`` 'R0', 'R1E', 'R1J', 'R2E', 'R2J' (as detected) or the
    same with '+c' (the constant offset of '<set>|<rung>|offset_nm' subtracted), ``split`` 'dev' or 'test'. Returns (dets, frame_indices): per-frame arrays
    (x_nm, y_nm, photons, sx_nm) and the stack indices of the frames, ready for ``score(dets, gt, frame_indices)``."""
    z = np.load(path, allow_pickle=False)
    off, flat = z['%s|%s|%s|off' % (set_name, row, split)], z['%s|%s|%s|flat' % (set_name, row, split)]
    a, b = (int(v) for v in z['%s|%s_frames' % (set_name, split)])
    dets = unpack_detections(off, flat)
    if len(dets) != b - a:
        raise ValueError('%d detection lists for frames %d-%d' % (len(dets), a, b - 1))
    return dets, np.arange(a, b)


def save_npz_deterministic(path, arrays):
    """``np.savez`` without time stamps: sorted keys, stored (uncompressed), fixed zip date, so the same arrays always give the same bytes."""
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_STORED) as z:
        for k in sorted(arrays):
            zi = zipfile.ZipInfo(k + '.npy', date_time=(1980, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_STORED
            zi.external_attr = 0o644 << 16
            buf = io.BytesIO()
            np.lib.format.write_array(buf, np.asanyarray(arrays[k]), allow_pickle=False)
            z.writestr(zi, buf.getvalue())


if __name__ == '__main__':
    if len(sys.argv) == 5 and sys.argv[1] == '--picasso-worker':
        _picasso_worker(sys.argv[2], sys.argv[3], json.loads(sys.argv[4]))
    elif len(sys.argv) == 3 and sys.argv[1] == '--picasso-server':
        _picasso_server(sys.argv[2])
    else:
        sys.exit('usage: smlm_common.py --picasso-worker frames.npy out.npz params_json | --picasso-server frames.npy   (import the module for everything else)')
