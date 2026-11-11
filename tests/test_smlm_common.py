"""Fast tests (no downloads) of the localization-benchmark harness `experiments/smlm_common.py`, the baseline ladder helpers of
`experiments/exp_smlm_picasso.py` and the downloader `data/smlm/fetch.py`: the matching, the metrics, the bootstrap, the loaders, the
post-filters, the offset correction, the packing of detections and the downloader on synthetic data with known answers. Only the
end-to-end tests of Picasso's options need Picasso: they run when SMLM_PICASSO_PYTHON names an interpreter that has it, else they are skipped."""
import hashlib
import http.server
import importlib.util
import math
import os
import sys
import threading
import zipfile

import numpy as np
import pytest
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
import smlm_common as C


def _gt(rows):
    """Ground truth array from (frame0, x_nm, y_nm, photons) tuples."""
    g = np.zeros(len(rows), dtype=C.GT_DTYPE)
    for i, (f, x, y, ph) in enumerate(rows):
        g[i] = (i + 1, f, x, y, 0.0, ph)
    return g


# ------------------------------------------------------------------ match
def test_match_perfect_and_small_offsets():
    gt = np.array([[1000., 1000.], [3000., 3000.], [5000., 1000.]])
    det = gt[[2, 0, 1]] + np.array([[10., -20.], [0., 0.], [-30., 40.]])           # shuffled, each within 50 nm
    m = C.match(det, gt)
    assert (m.tp, m.fp, m.fn) == (3, 0, 0)
    assert sorted(map(tuple, m.pairs)) == [(0, 2), (1, 0), (2, 1)]                  # (detection, ground-truth) indices


def test_match_tolerance_is_inclusive_and_sharp():
    gt = np.array([[1000., 1000.]])
    assert C.match(np.array([[1250., 1000.]]), gt, 250.0).tp == 1                     # exactly at the tolerance: a match
    assert C.match(np.array([[1250. + 1e-6, 1000.]]), gt, 250.0).tp == 0              # a hair beyond: not
    m = C.match(np.array([[1000. + 150., 1000. + 200.]]), gt, 250.0)                  # distance 250 along a diagonal (3-4-5)
    assert (m.tp, m.fp, m.fn) == (1, 0, 0)
    m = C.match(np.array([[1000. + 150., 1000. + 200.1]]), gt, 250.0)
    assert (m.tp, m.fp, m.fn) == (0, 1, 1)
    assert C.match(np.array([[1100., 1000.]]), gt, tol_nm=50.0).tp == 0                # the tolerance argument is honoured


def test_match_empty_frames():
    e = np.empty((0, 2))
    assert tuple(C.match(e, e))[:3] == (0, 0, 0)
    assert tuple(C.match(e, np.array([[1., 2.], [3., 4.]])))[:3] == (0, 0, 2)          # nothing detected: all false negatives
    assert tuple(C.match(np.array([[1., 2.]]), e))[:3] == (0, 1, 0)                    # nothing true: all false positives
    assert tuple(C.match([], [[1., 2.]]))[:3] == (0, 0, 1)                             # plain empty lists are accepted
    assert C.match(e, e).pairs.shape == (0, 2)


def test_match_duplicate_detections_count_once():
    gt = np.array([[1000., 1000.]])
    det = np.array([[1010., 1000.], [1000., 1030.], [995., 1000.]])                    # three detections of one emitter
    m = C.match(det, gt)
    assert (m.tp, m.fp, m.fn) == (1, 2, 0)
    assert tuple(m.pairs[0]) == (2, 0)                                                 # the nearest one is the true positive


def test_match_uses_extra_columns_only_as_payload_and_ignores_nan():
    gt = np.array([[1000., 1000.]])
    det = np.array([[1000., 1000., 5e3, 9.], [np.nan, 1000., 1., 1.]])                # extra columns; one failed fit
    m = C.match(det, gt)
    assert (m.tp, m.fp, m.fn) == (1, 1, 0)
    with pytest.raises(ValueError):
        C.match(np.array([1., 2., 3.]), gt)                                            # a 1-D array is ambiguous: refused


def test_match_is_optimal_not_greedy():
    """Greedy nearest-first pairing would take the 90 nm pair and lose the second pair; the assignment keeps both."""
    gt = np.array([[0., 0.], [200., 0.]])
    det = np.array([[110., 0.], [-140., 0.]])                                           # d(det0, gt1) = 90 is the smallest overall
    m = C.match(det, gt, 250.0)
    assert (m.tp, m.fp, m.fn) == (2, 0, 0)
    assert sorted(map(tuple, m.pairs)) == [(0, 1), (1, 0)]


def test_match_maximises_pairs_before_minimising_distance():
    """One detection close to the first emitter and within tolerance of the second; a second detection only near the first emitter: the
    cardinality-2 pairing wins even though it has the larger total distance."""
    gt = np.array([[0., 0.], [240., 0.]])
    det = np.array([[10., 0.], [-200., 0.]])
    m = C.match(det, gt, 250.0)
    assert m.tp == 2 and sorted(map(tuple, m.pairs)) == [(0, 1), (1, 0)]


# ------------------------------------------------------------------ score
def test_score_known_answer():
    gt = _gt([(0, 1000., 1000., 4000.), (0, 3000., 3000., 4000.),         # frame 0: two emitters
              (1, 2000., 2000., 4000.),                                    # frame 1: one emitter
              (5, 500., 500., 4000.)])                                     # frame 5 is not evaluated
    dets = [np.array([[1030., 1040.], [3000., 3000.], [5000., 5000.]]),    # frame 0: 2 TP (errors 50 nm and 0), 1 FP
            np.array([[2400., 2000.]])]                                     # frame 1: 400 nm away: FP, and the emitter is a FN
    r = C.score(dets, gt, [0, 1])
    assert (r['tp'], r['fp'], r['fn']) == (2, 2, 1)
    assert r['jaccard'] == pytest.approx(2 / 5)
    assert r['precision'] == pytest.approx(2 / 4)
    assert r['recall'] == pytest.approx(2 / 3)
    assert r['rmse_nm'] == pytest.approx(math.sqrt((30 ** 2 + 40 ** 2 + 0) / 2))     # sqrt(mean squared distance of the pairs)
    assert r['bias_x_nm'] == pytest.approx(15.0) and r['bias_y_nm'] == pytest.approx(20.0)
    assert r['rmse_centered_nm'] == pytest.approx(math.sqrt(15 ** 2 + 20 ** 2))                 # spread about the mean offset: sqrt(var dx + var dy)
    assert (r['n_det'], r['n_gt'], r['n_det_all'], r['n_gt_all']) == (4, 3, 4, 3)
    assert r['efficiency'] == pytest.approx(100 - math.sqrt((100 - 40) ** 2 + r['rmse_nm'] ** 2))
    assert list(r['per_frame']['frame']) == [0, 1] and list(r['per_frame']['tp']) == [2, 0]


def test_efficiency_formula():
    assert C.efficiency(1.0, 0.0) == pytest.approx(100.0)
    assert C.efficiency(0.8, 20.0) == pytest.approx(100 - math.sqrt(20 ** 2 + 20 ** 2))
    assert C.efficiency(0.8, 20.0, alpha_per_nm=0.5) == pytest.approx(100 - math.sqrt(20 ** 2 + 10 ** 2))


def test_score_empty_frames_and_degenerate_cases():
    gt = _gt([(0, 1000., 1000., 4000.)])
    r = C.score([np.empty((0, 2)), np.empty((0, 2))], gt, [0, 1])                     # frame 1 has no ground truth and no detection
    assert (r['tp'], r['fp'], r['fn']) == (0, 0, 1)
    assert r['jaccard'] == 0.0 and r['recall'] == 0.0
    assert math.isnan(r['precision']) and math.isnan(r['rmse_nm']) and math.isnan(r['efficiency'])      # 0/0 is NaN, not 0
    r = C.score([np.array([[1., 1.]])], gt, [3])                                        # a frame without ground truth: pure false positive
    assert (r['tp'], r['fp'], r['fn']) == (0, 1, 0)
    r = C.score([np.empty((0, 2))], gt, [3])                                            # nothing at all
    assert math.isnan(r['jaccard'])


def test_score_duplicate_detections_and_input_checks():
    gt = _gt([(0, 1000., 1000., 4000.)])
    r = C.score([np.array([[1000., 1000.], [1001., 1000.], [999., 1002.]])], gt, [0])
    assert (r['tp'], r['fp'], r['fn']) == (1, 2, 0) and r['jaccard'] == pytest.approx(1 / 3)
    with pytest.raises(ValueError):
        C.score([np.empty((0, 2))], gt, [0, 1])                                          # lengths differ
    with pytest.raises(ValueError):
        C.score([np.empty((0, 2))] * 2, gt, [0, 0])                                      # duplicate frame index


def test_score_border_and_photon_exclusions_follow_the_challenge_rules():
    gt = _gt([(0, 300., 3000., 4000.),            # inside the 450 nm border zone
              (0, 3000., 3000., 100.),            # dim
              (0, 5000., 5000., 4000.)])          # kept
    dets = [np.array([[310., 3000.],              # pairs with the border emitter: dropped with it
                      [3000., 3000.],             # pairs with the dim emitter: dropped with it
                      [5010., 5000.],             # TP
                      [200., 200.],               # unpaired, in the border zone: dropped
                      [4000., 1000.]])]           # unpaired inside the field: FP
    r = C.score(dets, gt, [0], border_nm=450.0, photon_min=1000.0)
    assert (r['tp'], r['fp'], r['fn']) == (1, 1, 0)
    assert (r['n_det_all'], r['n_gt_all']) == (5, 3)
    r0 = C.score(dets, gt, [0])
    assert (r0['tp'], r0['fp'], r0['fn']) == (3, 2, 0)                                    # defaults: everything counts
    # a missed excluded emitter is not a false negative; a missed kept one is
    r = C.score([np.empty((0, 2))], gt, [0], border_nm=450.0, photon_min=1000.0)
    assert (r['tp'], r['fp'], r['fn']) == (0, 0, 1)
    # photon_threshold: the dimmest quarter is below it
    g = _gt([(0, 1., 1., float(v)) for v in range(1, 101)])
    assert C.photon_threshold(g, 0.75) == pytest.approx(np.quantile(np.arange(1, 101), 0.25))


def test_score_frame_selection_is_by_value_not_position():
    gt = _gt([(7, 1000., 1000., 4000.), (9, 2000., 2000., 4000.)])
    r = C.score([np.array([[2000., 2000.]])], gt, [9])
    assert (r['tp'], r['fp'], r['fn']) == (1, 0, 0)                                        # frame 7's emitter belongs to another frame


# ------------------------------------------------------------------ bootstrap
def _per_frame(tp, fp, fn, sse=None):
    p = np.zeros(len(tp), dtype=C.PER_FRAME_DTYPE)
    p['tp'], p['fp'], p['fn'] = tp, fp, fn
    p['sse'] = 0.0 if sse is None else sse
    return p


def test_bootstrap_degenerate_when_all_frames_are_equal():
    p = _per_frame([2] * 50, [1] * 50, [1] * 50, [200.0] * 50)
    ci = C.bootstrap_ci(p, B=200, seed=3)
    for k, v in dict(jaccard=0.5, precision=2 / 3, recall=2 / 3, rmse_nm=10.0).items():
        assert ci[k][0] == pytest.approx(v) and ci[k][1] == pytest.approx(v)


def test_bootstrap_is_reproducible_and_brackets_the_estimate():
    rng = np.random.default_rng(1)
    n = 400
    tp = rng.integers(0, 2, n)
    fn = 1 - tp
    fp = rng.integers(0, 2, n)
    p = _per_frame(tp, fp, fn, tp * rng.uniform(100, 900, n))
    a, b = C.bootstrap_ci(p, B=300, seed=5), C.bootstrap_ci(p, B=300, seed=5)
    assert a == b                                                                         # same seed, same intervals
    assert a != C.bootstrap_ci(p, B=300, seed=6)
    est = C.pooled_metrics(tp.sum(), fp.sum(), fn.sum(), p['sse'].sum())
    for k in ('jaccard', 'precision', 'recall', 'rmse_nm', 'efficiency'):
        assert a[k][0] <= est[k] <= a[k][1]
    # the analytic width of a proportion: recall = mean(tp) over frames with one emitter each (fn = 1 - tp), sd = sqrt(r(1-r)/n)
    r = tp.mean()
    half = 1.96 * math.sqrt(r * (1 - r) / n)
    iid = C.bootstrap_ci(p, B=300, seed=5, block=1)                                       # independent frames: the i.i.d. bootstrap (block=1) is the right one here
    assert (iid['recall'][1] - iid['recall'][0]) / 2 == pytest.approx(half, rel=0.25)
    assert (a['recall'][1] - a['recall'][0]) / 2 == pytest.approx(half, rel=0.4)           # the default block bootstrap is about as wide on independent frames


def test_bootstrap_accepts_score_output_and_blocks():
    gt = _gt([(f, 1000., 1000., 4000.) for f in range(30)])
    dets = [np.array([[1000., 1000.]]) if f % 3 else np.empty((0, 2)) for f in range(30)]
    r = C.score(dets, gt, range(30))
    ci1 = C.bootstrap_ci(r, B=200)
    ci5 = C.bootstrap_ci(r['per_frame'], B=200, block=5)
    assert set(ci1) == {'jaccard', 'rmse_nm', 'precision', 'recall', 'efficiency'}
    assert ci1['recall'][0] <= 20 / 30 <= ci1['recall'][1] and ci5['recall'][0] <= 20 / 30 <= ci5['recall'][1]
    assert ci5['precision'] == pytest.approx((1.0, 1.0))                                  # no false positives in any resample
    with pytest.raises(ValueError):
        C.bootstrap_ci(_per_frame([], [], []))


def test_paired_bootstrap_difference():
    rng = np.random.default_rng(2)
    n = 300
    base_tp = rng.integers(0, 2, n)
    a = _per_frame(base_tp, np.zeros(n, int), 1 - base_tp, base_tp * 100.0)
    same = C.bootstrap_diff_ci(a, a, B=200)
    assert all(v == (0.0, 0.0) for v in same.values())                                    # a method against itself: exactly zero
    better = base_tp.copy()
    better[rng.random(n) < 0.3] = 1                                                        # B finds every emitter A finds, and some more
    b = _per_frame(better, np.zeros(n, int), 1 - better, better * 100.0)
    d = C.bootstrap_diff_ci(b, a, B=400, seed=1)
    assert d['recall'][0] > 0 and d['jaccard'][0] > 0                                      # paired: a clear, one-sided difference
    est = C.pooled_metrics(better.sum(), 0, (1 - better).sum(), (better * 100.0).sum())['recall'] - base_tp.mean()
    assert d['recall'][0] <= est <= d['recall'][1]
    sep = C.bootstrap_ci(b, B=400, seed=1)['recall'], C.bootstrap_ci(a, B=400, seed=1)['recall']
    assert (d['recall'][1] - d['recall'][0]) < (sep[0][1] - sep[0][0]) + (sep[1][1] - sep[1][0])       # tighter than the two separate intervals
    with pytest.raises(ValueError):
        C.bootstrap_diff_ci(a, _per_frame(base_tp[:-1], np.zeros(n - 1, int), 1 - base_tp[:-1]))      # different frame lists are refused


def test_pixel_nm_conversions_are_inverse_and_use_the_corner_convention():
    assert np.allclose(C.pixels_to_nm([0.0, 31.5]), [50.0, 3200.0])                       # the centre of the first pixel is at 50 nm
    p = np.array([[10.25, 3.0], [0.0, 63.0]])
    assert np.allclose(C.nm_to_pixels(C.pixels_to_nm(p)), p)
    assert np.allclose(C.nm_to_pixels([6400.0, 0.0]), [63.5, -0.5])                       # the field edges are half a pixel outside the outer pixel centres


# ------------------------------------------------------------------ loaders
def _write_stack(path, frames):
    ims = [Image.fromarray(f.astype(np.uint16)) for f in frames]
    ims[0].save(path, save_all=True, append_images=ims[1:], format='TIFF')


def test_load_frames_from_tiff_and_zip(tmp_path):
    rng = np.random.default_rng(0)
    stack = rng.integers(100, 5000, size=(7, 64, 64))
    tif = tmp_path / 'stack.tif'
    _write_stack(tif, stack)
    path = str(tif)
    a = C.load_frames(path, 0, 7)
    assert a.dtype == np.float64 and a.shape == (7, 64, 64)
    assert np.array_equal(a, stack)
    assert np.array_equal(C.load_frames(path, 2, 5), stack[2:5])                           # 0-based, stop exclusive
    assert np.array_equal(C.load_frames(path, 3), stack[3:])                               # stop=None: to the end
    assert C.load_frames(path, 4, 4).shape == (0, 64, 64)                                  # empty slice
    with pytest.raises(IndexError):
        C.load_frames(path, 0, 8)                                                          # no silent truncation
    with pytest.raises(ValueError):
        C.load_frames(path, 5, 2)
    zpath = tmp_path / 'seq.zip'
    with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_DEFLATED) as z:
        z.write(tif, 'inside/sequence-as-stack.tif')
    assert np.array_equal(C.load_frames(str(zpath), 1, 4), stack[1:4])
    two = tmp_path / 'two.zip'
    with zipfile.ZipFile(two, 'w') as z:
        z.write(tif, 'a.tif')
        z.write(tif, 'b.tif')
    with pytest.raises(ValueError):
        C.load_frames(str(two), 0, 1)                                                      # ambiguous archive


def test_load_frames_reads_big_endian_stacks(tmp_path):
    """The challenge stacks are big-endian 16-bit TIFFs (mode I;16B): write one by hand and read it back."""
    import struct
    stack = np.arange(3 * 4 * 5, dtype=np.uint16).reshape(3, 4, 5) * 7 + 100
    h, w, n = 4, 5, 3
    data_len = h * w * 2
    ifd_len = 2 + 9 * 12 + 4
    blob = bytearray(b'MM\x00\x2a' + struct.pack('>I', 8))
    off = 8
    for k in range(n):
        data_off = off + ifd_len
        nxt = data_off + data_len if k < n - 1 else 0
        ent = [(256, 4, w), (257, 4, h), (258, 3, 16), (259, 3, 1), (262, 3, 1), (273, 4, data_off), (277, 3, 1), (278, 3, h), (279, 4, data_len)]
        blob += struct.pack('>H', len(ent))
        for tag, typ, val in ent:
            blob += struct.pack('>HHI', tag, typ, 1) + (struct.pack('>HH', val, 0) if typ == 3 else struct.pack('>I', val))
        blob += struct.pack('>I', nxt) + stack[k].astype('>u2').tobytes()
        off = nxt
    p = tmp_path / 'be.tif'
    p.write_bytes(bytes(blob))
    a = C.load_frames(str(p), 0, 3)
    assert a.shape == (3, 4, 5) and np.array_equal(a, stack)


def test_to_photons():
    dn = np.array([[100.0, 106.0], [700.0, 94.0]])
    assert np.allclose(C.to_photons(dn), [[0.0, 1.0], [100.0, -1.0]])
    assert C.DN_PER_PHOTON == pytest.approx(6.0)
    assert C.to_photons(np.array([700], dtype=np.uint16))[0] == pytest.approx(100.0)       # integer input is converted, not wrapped


def test_load_gt_converts_frames_to_zero_based_and_sorts(tmp_path):
    p = tmp_path / 'activations.csv'
    p.write_text('Ground-truth,frame,xnano,ynano,znano,intensity \n'
                 '    1,  9039, 43.325, 839.830, 340.704, 2473.663 \n'
                 '    2,     1, 10.000, 20.000, -5.000, 100.5 \n'
                 '    3,  9039, 6399.000, 6399.500, 0.000, 50.0 \n')
    g = C.load_gt(str(p))
    assert g.dtype.names[:4] == ('id', 'frame', 'x_nm', 'y_nm') and 'photons' in g.dtype.names
    assert list(g['frame']) == [0, 9038, 9038]                                             # file frame 1 is stack index 0
    assert list(g['id']) == [2, 1, 3]                                                      # sorted by (frame, id)
    assert g['x_nm'][1] == pytest.approx(43.325) and g['photons'][0] == pytest.approx(100.5) and g['z_nm'][0] == pytest.approx(-5.0)
    assert C.load_gt(str(p), frame_base=0)['frame'].min() == 1
    bad = tmp_path / 'positions.csv'
    bad.write_text('43.32457, 839.82958, 340.70413\n39.00505, 850.19918, 331.85586\n')
    with pytest.raises(ValueError):
        C.load_gt(str(bad))                                                                # the fluorophore list is not a ground truth


# ------------------------------------------------------------------ Picasso table conversion (no Picasso needed)
def test_picasso_table_to_detections_uses_the_pixel_centre_convention():
    tab = dict(frame=np.array([2, 0, 2]), x=np.array([10.0, 31.5, 0.0]), y=np.array([20.0, 0.0, 63.0]),
               photons=np.array([450.0, 900.0, 90.0]), sx=np.array([1.0, 1.1, 1.2]), sy=np.array([1.0, 1.0, 1.0]),
               net_gradient=np.array([7e3, 8e3, 9e3]))
    d = C.detections_from_table(tab, 4)
    assert len(d) == 4 and [len(a) for a in d] == [1, 0, 2, 0]                             # frames 1 and 3 have no detections
    assert d[0][0, 0] == pytest.approx(3200.0) and d[0][0, 1] == pytest.approx(50.0)       # pixel index 31.5 -> (31.5 + 0.5) * 100 nm
    assert d[0][0, 2] == pytest.approx(1000.0)                                              # photoelectrons / QE
    assert d[2][0, 0] == pytest.approx(1050.0) and d[2][0, 1] == pytest.approx(2050.0)
    assert d[2][0, 3] == pytest.approx(100.0) and d[2][0, C.NG_COL] == pytest.approx(7e3)
    assert C.detections_from_table({k: v[:0] for k, v in tab.items()}, 2)[1].shape == (0, 6)
    with pytest.raises(ValueError):
        C.detections_from_table(tab, 2)                                                    # a frame index beyond the stack
    with pytest.raises(ValueError):
        C.tune_picasso(np.zeros((1, 64, 64)), _gt([]), list(range(21)))                    # the budget is 20 values


# ------------------------------------------------------------------ the downloader (local server only, nothing from the internet)
def _load_fetch():
    spec = importlib.util.spec_from_file_location('smlm_fetch', os.path.join(os.path.dirname(__file__), '..', 'data', 'smlm', 'fetch.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _LocalServer:
    """HTTP/1.1 server for one payload with HEAD and Range support; ``cut_first`` closes the first full-file GET after a third of the body."""

    def __init__(self, payload, cut_first=False):
        outer = self
        self.payload, self.starts, self.cut = payload, [], cut_first

        class Handler(http.server.BaseHTTPRequestHandler):
            protocol_version = 'HTTP/1.1'

            def log_message(self, *a):
                pass

            def do_HEAD(self):
                self.send_response(200)
                self.send_header('Content-Length', str(len(outer.payload)))
                self.send_header('Accept-Ranges', 'bytes')
                self.end_headers()

            def do_GET(self):
                rng = self.headers.get('Range')
                start = int(rng.split('=')[1].split('-')[0]) if rng else 0
                outer.starts.append(start)
                body = outer.payload[start:]
                self.send_response(206 if rng else 200)
                self.send_header('Content-Length', str(len(body)))
                if rng:
                    self.send_header('Content-Range', 'bytes %d-%d/%d' % (start, len(outer.payload) - 1, len(outer.payload)))
                self.end_headers()
                if outer.cut and start == 0:
                    outer.cut = False
                    self.wfile.write(body[:len(body) // 3])
                    self.wfile.flush()
                    self.close_connection = True
                    return
                self.wfile.write(body)

        self.httpd = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        self.url = 'http://127.0.0.1:%d/file.bin' % self.httpd.server_address[1]
        threading.Thread(target=self.httpd.serve_forever, daemon=True).start()

    def close(self):
        self.httpd.shutdown()
        self.httpd.server_close()


@pytest.fixture
def local_payload(monkeypatch):
    for k in ('no_proxy', 'NO_PROXY'):
        monkeypatch.setenv(k, '127.0.0.1,localhost')
    payload = np.random.default_rng(0).integers(0, 256, 400_000, dtype=np.uint8).tobytes()
    servers = []

    def make(**kw):
        servers.append(_LocalServer(payload, **kw))
        return servers[-1]
    yield payload, make
    for sv in servers:
        sv.close()


def test_fetch_head_hash_and_resume_from_a_part_file(tmp_path, local_payload):
    F = _load_fetch()
    payload, make = local_payload
    sv = make()
    total, _, _ = F.head(sv.url)
    assert total == len(payload)
    dest = str(tmp_path / 'sub' / 'file.bin')
    os.makedirs(os.path.dirname(dest))
    with open(dest + '.part', 'wb') as f:
        f.write(payload[:123_456])                                                         # an interrupted earlier download
    F.download(sv.url, dest, total, log=lambda *a: None)
    assert open(dest, 'rb').read() == payload and not os.path.exists(dest + '.part')
    assert sv.starts == [123_456]                                                          # only the missing tail was requested
    assert F.sha256_of(dest) == hashlib.sha256(payload).hexdigest()


def test_fetch_survives_a_dropped_connection(tmp_path, local_payload):
    F = _load_fetch()
    payload, make = local_payload
    sv = make(cut_first=True)
    dest = str(tmp_path / 'file.bin')
    F.download(sv.url, dest, len(payload), log=lambda *a: None)
    assert open(dest, 'rb').read() == payload
    assert sv.starts[0] == 0 and len(sv.starts) >= 2 and sv.starts[1] > 0                  # the second request resumed where the first stopped


def test_fetch_refuses_incomplete_or_oversized_parts_and_repo_data_dir(tmp_path, local_payload):
    F = _load_fetch()
    payload, make = local_payload
    sv = make()
    dest = str(tmp_path / 'file.bin')
    with open(dest + '.part', 'wb') as f:
        f.write(payload + b'extra')                                                        # longer than the announced size: discarded, not trusted
    F.download(sv.url, dest, len(payload), log=lambda *a: None)
    assert open(dest, 'rb').read() == payload
    with pytest.raises(SystemExit):
        F.main(['--data-dir', os.path.join(F.REPO, 'raw_data_must_not_live_here'), '--head-only'])      # raw data never goes into the repository


def test_fetch_zip_members(tmp_path):
    F = _load_fetch()
    z = tmp_path / 's.zip'
    with zipfile.ZipFile(z, 'w', zipfile.ZIP_DEFLATED) as zz:
        zz.writestr('stack.tif', b'x' * 1000)
    assert F.zip_members(str(z)) == [('stack.tif', 1000)]


# ==================================================================================================================================
# the rules of the assessment: non-finite detections, the "stealing" case (SM1-7)
# ==================================================================================================================================
def test_non_finite_detections_are_false_positives_in_every_setting():
    gt = _gt([(0, 3000., 3000., 4000.)])
    dets = [np.array([[3010., 3000.], [np.nan, 3000.], [3000., np.inf], [-np.inf, np.nan]])]
    for kw in ({}, dict(border_nm=450.0), dict(photon_min=1000.0), dict(border_nm=450.0, photon_min=1000.0)):
        r = C.score(dets, gt, [0], **kw)
        assert (r['tp'], r['fp'], r['fn']) == (1, 3, 0), kw                       # the three non-finite ones count, also when the border exclusion is on
        assert r['n_det'] == 4 and r['n_det_all'] == 4
    # a finite unpaired detection inside the border zone still vanishes with the border exclusion (the challenge's rule), a non-finite one does not
    r = C.score([np.array([[100., 3000.], [np.nan, np.nan]])], gt, [0], border_nm=450.0)
    assert (r['tp'], r['fp'], r['fn']) == (0, 1, 1)
    r = C.score([np.array([[100., 3000.], [np.nan, np.nan]])], gt, [0])
    assert (r['tp'], r['fp'], r['fn']) == (0, 2, 1)


def test_an_excluded_dim_emitter_steals_the_detection_of_a_kept_one():
    """Documented behaviour (the pairing precedes the exclusions): the detection is closer to the dim emitter, so it pairs with that one; the pair is
    then dropped and the kept emitter is a false negative although the detection lies within the tolerance of it."""
    gt = _gt([(0, 3000., 3000., 100.),               # dim, excluded by photon_min
              (0, 3200., 3000., 4000.)])             # kept, 150 nm from the detection
    dets = [np.array([[3050., 3000.]])]              # 50 nm from the dim emitter, 150 nm from the kept one: both within 250 nm
    r0 = C.score(dets, gt, [0])
    assert (r0['tp'], r0['fp'], r0['fn']) == (1, 0, 1)                                   # no exclusion: one pair (with the nearer, dim emitter), one miss
    r = C.score(dets, gt, [0], photon_min=1000.0)
    assert (r['tp'], r['fp'], r['fn']) == (0, 0, 1)                                      # the kept emitter is lost to the excluded one
    # without the dim emitter the same detection is a true positive of the kept one
    r = C.score(dets, _gt([(0, 3200., 3000., 4000.)]), [0], photon_min=1000.0)
    assert (r['tp'], r['fp'], r['fn']) == (1, 0, 0)
    # the border exclusion steals in the same way
    gtb = _gt([(0, 300., 3000., 4000.), (0, 480., 3000., 4000.)])                        # 300 nm: inside the 450 nm border zone (excluded); 480 nm: kept
    rb = C.score([np.array([[330., 3000.]])], gtb, [0], border_nm=450.0)
    assert (rb['tp'], rb['fp'], rb['fn']) == (0, 0, 1)


# ==================================================================================================================================
# the block bootstrap (SM1-6)
# ==================================================================================================================================
def test_run_lengths_count_consecutive_frames_of_one_fluorophore():
    gt = _gt([(0, 100., 100., 1.), (1, 100., 100., 1.), (2, 100., 100., 1.),            # one fluorophore on for 3 frames
              (5, 100., 100., 1.),                                                         # the same site again after a gap: a new run of 1
              (1, 900., 900., 1.), (2, 900., 900., 1.),                                    # another fluorophore: a run of 2
              (7, 500., 500., 1.)])                                                        # a run of 1
    rl = C.run_lengths(gt)
    assert rl['runs'] == 4 and sorted(rl['lengths']) == [1, 1, 2, 3]
    assert rl['mean'] == pytest.approx(7 / 4) and rl['ge3'] == pytest.approx(1 / 4) and list(rl['dist'][:3]) == [2, 1, 1]


def test_default_bootstrap_block_is_ten_and_the_iid_bootstrap_stays_available():
    assert C.DEFAULT_BLOCK == 10
    rng = np.random.default_rng(3)
    n = 120
    tp = rng.integers(0, 3, n)
    p = _per_frame(tp, rng.integers(0, 2, n), rng.integers(0, 2, n), tp * rng.uniform(100, 900, n))
    assert C.bootstrap_ci(p, B=100, seed=1) == C.bootstrap_ci(p, B=100, seed=1, block=10)
    assert C.bootstrap_ci(p, B=100, seed=1) != C.bootstrap_ci(p, B=100, seed=1, block=1)
    a, b = p, _per_frame(tp[::-1], rng.integers(0, 2, n), rng.integers(0, 2, n), tp[::-1] * 300.0)
    b['frame'] = a['frame']
    assert C.bootstrap_diff_ci(a, b, B=100, seed=1) == C.bootstrap_diff_ci(a, b, B=100, seed=1, block=10)
    assert C.bootstrap_diff_ci(a, b, B=100, seed=1) != C.bootstrap_diff_ci(a, b, B=100, seed=1, block=1)


def _persistent_toy(rng, F=300, rate=0.35, q=0.6, mean_len=3.0):
    """Per-frame (tp, fn) of a detector with a per-EMITTER detectability q: emitters switch on at random frames (``rate`` new ones per frame), stay on for a
    geometric number of consecutive frames (mean ``mean_len``, like the 2.3 frames of the real sets) and are found in all of their frames or in none."""
    tp = np.zeros(F)
    fn = np.zeros(F)
    for f in range(F):
        for _ in range(rng.poisson(rate)):
            length = rng.geometric(1.0 / mean_len)
            hit = rng.random() < q
            sl = slice(f, min(F, f + length))
            (tp if hit else fn)[sl] += 1
    return tp, fn


def test_block_bootstrap_covers_persistent_emitters_where_the_iid_bootstrap_undercovers():
    """Emitters that stay on for several frames make neighbouring frames dependent: the interval of the i.i.d. frame bootstrap is too narrow (coverage far below
    95 %), blocks of 10 frames cover clearly better (the independent referee measured 64 % against 91 % for F = 500 on this kind of toy)."""
    q = 0.6                                                        # the true pooled recall (= Jaccard here): detectability is independent of the on-time
    rng = np.random.default_rng(11)
    cover = {1: 0, 10: 0}
    reps = 120
    for _ in range(reps):
        tp, fn = _persistent_toy(rng, q=q)
        p = _per_frame(tp, np.zeros_like(tp), fn)
        for block in cover:
            lo, hi = C.bootstrap_ci(p, B=150, seed=0, block=block)['recall']
            cover[block] += lo <= q <= hi
    iid, blk = cover[1] / reps, cover[10] / reps
    assert iid < 0.80                                              # undercoverage of the i.i.d. bootstrap (nominal 0.95)
    assert blk > iid + 0.12 and blk > 0.78                         # the block bootstrap is clearly better


def test_pooled_bootstrap_helpers_reproduce_the_metric_bootstrap_and_accept_other_statistics():
    rng = np.random.default_rng(5)
    n = 90
    tp = rng.integers(0, 3, n)
    a = _per_frame(tp, rng.integers(0, 2, n), rng.integers(0, 2, n), tp * rng.uniform(100, 900, n))
    b = _per_frame(rng.integers(0, 3, n), rng.integers(0, 2, n), rng.integers(0, 2, n))
    b['sse'] = b['tp'] * 400.0
    b['frame'] = a['frame']
    rmse = lambda s: np.sqrt(s['sse'] / s['tp'])                                      # noqa: E731
    ci = C.bootstrap_ci(a, B=120, seed=2)['rmse_nm']
    assert C.bootstrap_pooled_ci(dict(sse=a['sse'], tp=a['tp']), rmse, B=120, seed=2) == pytest.approx(ci)
    d = C.bootstrap_diff_ci(a, b, B=120, seed=2)['rmse_nm']
    assert C.bootstrap_pooled_diff_ci(dict(sse=a['sse'], tp=a['tp']), dict(sse=b['sse'], tp=b['tp']), rmse, B=120, seed=2) == pytest.approx(d)
    with pytest.raises(ValueError):
        C.bootstrap_pooled_diff_ci(dict(sse=a['sse'], tp=a['tp']), dict(sse=b['sse'][:-1], tp=b['tp'][:-1]), rmse)
    # extended keys: the bias and the centred RMSE are resampled too
    a['sx'] = a['tp'] * 10.0
    a['sy'] = a['tp'] * -5.0
    ext = C.bootstrap_ci(a, B=100, seed=2, keys=C.BOOT_KEYS_EXT)
    assert set(ext) == set(C.BOOT_KEYS_EXT)
    assert ext['bias_x_nm'][0] == pytest.approx(10.0) and ext['bias_y_nm'][1] == pytest.approx(-5.0)       # a constant error per true positive
    assert ext['rmse_centered_nm'][1] < ext['rmse_nm'][1]
    assert set(C.bootstrap_ci(a, B=50)) == set(C._BOOT_KEYS)                           # the default keys are unchanged


# ==================================================================================================================================
# post-filters, tuning of the post-filters, offsets, packing (the baseline ladder)
# ==================================================================================================================================
def _det(rows):
    """(n, 6) detection array in DET_COLUMNS order from (x, y, photons, sx, sy, net_gradient) tuples."""
    return np.array(rows, dtype=float).reshape(-1, len(C.DET_COLUMNS))


def test_filter_detections_strict_inequalities_abs_widths_and_nan():
    d = _det([(1, 1, 500., 130., 130., 6000.),
              (2, 2, 100., 130., 130., 6000.),          # few photons
              (3, 3, 500., -300., 120., 6000.),         # negative width: |sx| = 300
              (4, 4, 500., 130., np.nan, 6000.),        # NaN width
              (5, 5, 500., 130., 130., 4000.)])         # low net gradient
    assert len(C.filter_detections([d], {})[0]) == 5 and len(C.filter_detections([d], {'min_net_gradient': None})[0]) == 5
    assert list(C.filter_detections([d], {'min_net_gradient': 4000.0})[0][:, 0]) == [1, 2, 3, 4]               # strictly above: 4000 itself fails
    assert list(C.filter_detections([d], {'min_photons': 100.0})[0][:, 0]) == [1, 3, 4, 5]
    assert list(C.filter_detections([d], {'max_sigma_nm': 200.0})[0][:, 0]) == [1, 2, 5]                       # |-300| fails; NaN fails
    assert list(C.filter_detections([d], {'min_sigma_nm': 125.0})[0][:, 0]) == [1, 2, 5]                       # both axes above 125 nm: row 3 has sy = 120; NaN fails
    both = C.filter_detections([d], {'min_net_gradient': 5000.0, 'min_photons': 200.0, 'min_sigma_nm': 100.0, 'max_sigma_nm': 400.0})[0]
    assert list(both[:, 0]) == [1, 3]
    e = np.empty((0, 6))
    assert C.filter_detections([e, d], {'min_photons': 1e9})[0].shape == (0, 6)                                # empty frames pass through
    assert C.filter_detections([d], C._merge_params({}))[0] is not None                                        # the Picasso parameter dict is accepted as a filter dict
    assert len(C.filter_detections([d], dict(C.PICASSO_DEFAULTS, min_net_gradient=None))[0]) == 5


def test_median_pad_and_parameter_validation():
    fr = np.arange(2 * 4 * 6, dtype=np.uint16).reshape(2, 4, 6) + 10
    out = C.median_pad(fr, 2)
    assert out.shape == (2, 8, 10) and out.dtype == np.uint16
    assert np.array_equal(out[:, 2:-2, 2:-2], fr)
    med = np.median(fr.reshape(2, -1), axis=1)
    assert out[0, 0, 0] == int(med[0]) and out[1, 7, 9] == int(med[1])                         # the frame's own median, truncated to the integer type
    assert np.all(out[0, :2] == out[0, 0, 0]) and np.all(out[0, :, -2:] == out[0, 0, 0])
    assert C.median_pad(fr, 0) is fr or np.array_equal(C.median_pad(fr, 0), fr)
    assert C._merge_params(None) == C.PICASSO_DEFAULTS
    for bad in (dict(box=6), dict(fit_box=4), dict(box=1), dict(pad=-1), dict(pad=1.5), dict(gauss_sigma=-1.0), dict(nonsense=1)):
        with pytest.raises(ValueError):
            C._merge_params(bad)
    assert C._merge_params(dict(box=5, fit_box=9, pad=4, gauss_sigma=1.5, fit='gaussmle'))['fit_box'] == 9


def _toy_candidates(F=60):
    """One ground-truth emitter per frame at (1000, 1000); candidates: the true spot (ng 6000-7800, 130 nm wide, 20 nm off), in every sixth frame instead a wide spot
    (sx 500 nm, 200 nm off, ng 9000: the net gradient cannot separate it), and a false peak elsewhere (ng 500-2900, 100 photons)."""
    gt = _gt([(f, 1000., 1000., 4000.) for f in range(F)])
    dets = []
    for f in range(F):
        rows = [(1020., 1000., 3000., 130., 130., 6000. + 30 * f)]
        if f % 6 == 0:
            rows = [(1200., 1000., 3000., 500., 500., 9000.)]
        rows.append((4000., 4000., 100., 200., 200., 500. + 40 * f))                  # ng 500 .. 2860
        dets.append(_det(rows))
    return gt, dets


def test_tune_postfilters_finds_the_known_optimum_and_shares_its_cache():
    gt, dets = _toy_candidates()
    idx = np.arange(60)
    grids = {'min_net_gradient': [1000., 2000., 3000., 3500., 5000., 7000.], 'min_photons': [None, 50., 200.], 'min_sigma_nm': [None, 100.],
             'max_sigma_nm': [None, 200., 250., 600.]}
    cache = {}
    rJ = C.tune_postfilters(dets, gt, idx, grids, objective='jaccard', cache=cache)
    assert rJ['metrics']['jaccard'] == pytest.approx(1.0)                                  # every spot true, every false peak removed
    assert rJ['filters']['min_net_gradient'] == 3000.0 or rJ['filters']['min_photons'] is not None   # by net gradient (> 2860) or by photons: the first sufficient value wins
    assert rJ['filters']['max_sigma_nm'] is None                                          # the wide spots are true positives for the Jaccard index: not removed
    assert rJ['first']['filters']['min_photons'] is None and rJ['converged'] and rJ['n_eval'] == len(cache)
    n_after_J = len(cache)
    rE = C.tune_postfilters(dets, gt, idx, grids, objective='efficiency', cache=cache)
    assert rE['filters']['max_sigma_nm'] == 200.0                                          # the efficiency prefers dropping the badly localized wide spots (first grid value that does)
    assert rE['metrics']['jaccard'] == pytest.approx(50 / 60) and rE['metrics']['rmse_nm'] == pytest.approx(20.0)
    assert rE['value'] > rJ['metrics']['efficiency'] + 30
    assert len(cache) > n_after_J and rE['n_eval'] == len(cache)                          # one shared cache: combinations scored for the Jaccard ascent are not scored again
    # the value of the final configuration is what score() says for the filtered detections
    q = C.score(C.filter_detections(dets, rE['filters']), gt, idx)
    assert q['efficiency'] == pytest.approx(rE['metrics']['efficiency'])
    with pytest.raises(ValueError):
        C.tune_postfilters(dets, gt, idx, {'min_bananas': [1.0]})
    # max_rounds bounds the number of passes
    assert C.tune_postfilters(dets, gt, idx, grids, objective='jaccard', max_rounds=1)['rounds'] == 1


def test_shift_detections_subtracts_a_constant_and_keeps_other_columns():
    d = _det([(1000., 2000., 500., 130., 130., 6000.), (3000., 4000., 600., 140., 140., 7000.)])
    out = C.shift_detections([d, np.empty((0, 6)), np.empty(0)], 10.0, -5.0)
    assert np.allclose(out[0][:, 0], [990., 2990.]) and np.allclose(out[0][:, 1], [2005., 4005.]) and np.allclose(out[0][:, 2:], d[:, 2:])
    assert d[0, 0] == 1000.0                                                              # the input is not modified
    assert out[1].shape == (0, 6) and out[2].size == 0
    # an offset estimated on matched pairs removes the mean error: scoring the shifted detections gives bias 0 and the centred RMSE as RMSE
    gt = _gt([(0, 1000., 1000., 4000.), (1, 3000., 3000., 4000.)])
    dets = [np.array([[1012., 991.]]), np.array([[3008., 2995.]])]
    r = C.score(dets, gt, [0, 1])
    r2 = C.score(C.shift_detections(dets, r['bias_x_nm'], r['bias_y_nm']), gt, [0, 1])
    assert r2['bias_x_nm'] == pytest.approx(0.0, abs=1e-9) and r2['bias_y_nm'] == pytest.approx(0.0, abs=1e-9)
    assert r2['rmse_nm'] == pytest.approx(r['rmse_centered_nm'])


def test_matched_errors_agree_with_score_and_expose_the_matched_detection():
    rng = np.random.default_rng(8)
    rows = [(int(f), float(rng.uniform(0, 6400)), float(rng.uniform(0, 6400)), float(rng.uniform(100, 8000))) for f in rng.integers(0, 40, 70)]
    gt = _gt(rows)
    frames = np.arange(5, 35)
    dets = []
    for f in frames:
        g = gt[gt['frame'] == f]
        pts = [[x + rng.normal(0, 80), y + rng.normal(0, 80), 1.0 * k] for k, (x, y) in enumerate(zip(g['x_nm'], g['y_nm'])) if rng.random() < 0.7]
        pts.append([float(rng.uniform(0, 6400)), float(rng.uniform(0, 6400)), 99.0])
        dets.append(np.array(pts))
    r = C.score(dets, gt, frames)
    me = C.matched_errors(dets, gt, frames)
    assert me.dtype == C.MATCH_DTYPE and len(me) == r['n_gt_all']
    assert int(me['matched'].sum()) == r['tp']
    assert float((me['dx'][me['matched']] ** 2 + me['dy'][me['matched']] ** 2).sum()) == pytest.approx(r['rmse_nm'] ** 2 * r['tp'])
    assert float(me['dx'][me['matched']].mean()) == pytest.approx(r['bias_x_nm']) and float(me['dy'][me['matched']].mean()) == pytest.approx(r['bias_y_nm'])
    assert np.all(np.isnan(me['dx'][~me['matched']])) and np.all(me['det'][~me['matched']] == -1)
    for rec in me[me['matched']][:10]:                                                   # det indexes the detection of the frame, gt_row the ground-truth row
        d = dets[rec['pos']][rec['det']]
        assert d[0] - gt['x_nm'][rec['gt_row']] == pytest.approx(rec['dx']) and gt['frame'][rec['gt_row']] == rec['frame'] == frames[rec['pos']]
    # same ground truth and frames: the same rows in the same order for another method, so true positives can be intersected
    me2 = C.matched_errors([d[:1] for d in dets], gt, frames)
    assert np.array_equal(me['gt_row'], me2['gt_row']) and np.array_equal(me['pos'], me2['pos'])
    with pytest.raises(ValueError):
        C.matched_errors(dets[:-1], gt, frames)


def test_pack_unpack_roundtrip_and_deterministic_npz(tmp_path):
    dets = [_det([(1., 2., 3., 4., 5., 6.), (7., 8., 9., 10., 11., 12.)]), np.empty((0, 6)), _det([(13., 14., 15., 16., 17., 18.)]), np.empty(0)]
    off, flat = C.pack_detections(dets)
    assert off.dtype == np.int64 and list(off) == [0, 2, 2, 3, 3] and flat.shape == (3, 4)
    back = C.unpack_detections(off, flat)
    assert [len(b) for b in back] == [2, 0, 1, 0] and np.array_equal(back[0], dets[0][:, :4]) and np.array_equal(back[2], dets[2][:, :4])
    off0, flat0 = C.pack_detections([np.empty((0, 6))] * 3)
    assert list(off0) == [0, 0, 0, 0] and flat0.shape == (0, 4)
    gt = _gt([(0, 1., 2., 4000.)])
    r = C.score(back[:1], gt, [0])                                                       # the unpacked arrays are scoreable as they are (x_nm, y_nm first)
    assert r['tp'] + r['fn'] == 1
    arrays = {'b|flat': flat, 'a|off': off, 'text': np.array('hello {"json": 1}'), 'cols': np.array(['x_nm', 'y_nm'])}
    p1, p2 = tmp_path / 'one.npz', tmp_path / 'two.npz'
    C.save_npz_deterministic(str(p1), arrays)
    import time as _t
    _t.sleep(1.1)                                                                         # a different clock second must not change a byte
    C.save_npz_deterministic(str(p2), dict(reversed(list(arrays.items()))))              # nor does the insertion order
    assert p1.read_bytes() == p2.read_bytes()
    z = np.load(p1, allow_pickle=False)
    pl = tmp_path / 'ladder.npz'
    C.save_npz_deterministic(str(pl), {'S|R2E+c|test|off': off, 'S|R2E+c|test|flat': flat, 'S|test_frames': np.array([1000, 1004])})
    got, fidx = C.load_ladder_detections(str(pl), 'S', 'R2E+c', 'test')
    assert list(fidx) == [1000, 1001, 1002, 1003] and [len(g) for g in got] == [2, 0, 1, 0] and np.array_equal(got[0], dets[0][:, :4])
    C.save_npz_deterministic(str(pl), {'S|R2E+c|test|off': off, 'S|R2E+c|test|flat': flat, 'S|test_frames': np.array([1000, 1003])})
    with pytest.raises(ValueError):
        C.load_ladder_detections(str(pl), 'S', 'R2E+c', 'test')                          # a frame range that does not match the packed frames is refused
    assert sorted(z.files) == sorted(arrays) and np.array_equal(z['b|flat'], flat) and str(z['text']) == 'hello {"json": 1}' and list(z['cols']) == ['x_nm', 'y_nm']


# ==================================================================================================================================
# the ladder script: declared grids and helpers
# ==================================================================================================================================
def _load_ladder():
    spec = importlib.util.spec_from_file_location('exp_smlm_picasso', os.path.join(os.path.dirname(__file__), '..', 'experiments', 'exp_smlm_picasso.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_ladder_grids_are_declared_consistently():
    E = _load_ladder()
    assert len(E.struct_list()) == 3 * 3 * 3 * 2 * 2 == 108 and len({E.skey(s) for s in E.struct_list()}) == 108
    assert E.skey(E.DEFAULT_STRUCT) in {E.skey(s) for s in E.struct_list()}                    # the default structure is a point of the grid, so R2 contains R1
    for noise in ('N1', 'N2'):
        assert set(E.GRID_R1[noise]) <= set(E.FILTER_GRIDS[noise]['min_net_gradient'])        # the R1 grid lies inside the R2 net-gradient grid
        assert E.DEFAULT_NG in E.GRID_R1[noise] and len(E.GRID_R1[noise]) <= 20
        assert list(E.FILTER_GRIDS[noise]) == ['min_net_gradient', 'min_photons', 'min_sigma_nm', 'max_sigma_nm'] and set(E.FILTER_GRIDS[noise]) <= set(C.POST_FILTER_KEYS)
        assert all(g == sorted(g, key=lambda v: (v is not None, v)) for g in E.FILTER_GRIDS[noise].values())      # ascending (None first): the ties of the ascent go to the earlier value
    iss = E.iss_pic_grid()
    if iss is not None:                                                                       # the N2 grid of the earlier tuning is the one of the sliding-ISS experiment
        assert E.GRID_R1['N2'] == iss['N2'] and E.GRID_R1['N1'] == iss['N1']
    assert E.noise_of('MT0.N2.HD') == 'N2' and E.noise_of('MT0.N1.LD') == 'N1'
    assert set(E.PROTOCOL) == set(C.SETS) and E.PROTOCOL['MT0.N2.LD'] == dict(dev=(0, 1000), test=(1000, 1500)) and E.PROTOCOL['MT0.N1.HD'] == dict(dev=(0, 100), test=(100, 500))
    assert E.describe({}) == 'box 7 px, raw frames, no padding, gausslq, net gradient > 5000'
    assert 'ident 5 / fit 9 px' in E.describe(dict(box=5, fit_box=9, gauss_sigma=1.5, pad=4, fit='gaussmle', min_photons=200.0, max_sigma_nm=300.0))


def test_argmax_stability_of_a_threshold_curve():
    E = _load_ladder()
    F = 80
    rng = np.random.default_rng(4)
    clear = _per_frame(np.ones(F), np.zeros(F), np.zeros(F), np.full(F, 100.0))                 # threshold 1: every emitter found, error 10 nm
    worse = _per_frame(np.ones(F) * (rng.random(F) < 0.5), np.zeros(F), np.ones(F), np.full(F, 100.0))
    out = E.argmax_stability({1.0: clear, 2.0: worse, 3.0: worse}, B=50, seed=1)
    assert out['efficiency'][0] == (1.0, 1.0) and out['jaccard'][0] == (1.0, 1.0)                  # a clear winner wins every resample
    assert [t for t, _ in out['efficiency']] == [1.0, 2.0, 3.0] and sum(v for _, v in out['efficiency']) == pytest.approx(1.0)
    # two thresholds that differ only by noise share the resamples
    a = _per_frame(rng.integers(0, 2, F), np.zeros(F, int), rng.integers(0, 2, F), np.full(F, 100.0))
    b = _per_frame(rng.integers(0, 2, F), np.zeros(F, int), rng.integers(0, 2, F), np.full(F, 100.0))
    sh = dict(E.argmax_stability({1.0: a, 2.0: b}, B=200, seed=2)['jaccard'])
    assert 0.05 < sh[1.0] < 0.95 and sh[1.0] + sh[2.0] == pytest.approx(1.0)


def test_ladder_helpers_stealing_common_tp_and_neighbours():
    E = _load_ladder()
    gt = _gt([(0, 3000., 3000., 100.), (0, 3200., 3000., 4000.), (1, 3000., 3000., 4000.)])
    dets = [np.array([[3050., 3000.]]), np.array([[3000., 3000.]])]
    stolen, kept = E.count_stealing(dets, gt, np.array([0, 1]), pmin=1000.0)
    assert (stolen, kept) == (1, 2)                                                           # frame 0: the kept emitter lost its only detection to the dim one
    me = C.matched_errors(dets, gt, [0, 1])
    c = E.common_columns(me, me['matched'], 2)
    assert np.allclose(c['tp'], [1, 1]) and np.allclose(c['sse'], [50. ** 2, 0.0])
    frac, n = E.neighbour_fractions(gt, np.array([0, 1]))
    assert n == 3 and np.allclose(frac, [2 / 3, 2 / 3])                                       # the two emitters of frame 0 are 200 nm apart
    results = {E.skey(s): {'struct': s, 'E': {'metrics': {'jaccard': 0.5, 'efficiency': 10.0, 'rmse_nm': 30.0}}} for s in E.struct_list()}
    results[E.skey(E.DEFAULT_STRUCT)]['E']['metrics'] = {'jaccard': 0.4, 'efficiency': 8.0, 'rmse_nm': 31.0}
    chosen = collections_first(results)
    rows = E.ablation_rows(results, chosen, 'E', dict(photons=dict(jaccard=0.1, efficiency=1.0, rmse_nm=1.0), width=dict(jaccard=0.2, efficiency=2.0, rmse_nm=2.0)))
    labels = [r[0] for r in rows]
    assert labels[0] == 'chosen configuration' and any(l.startswith('photon filter -> off') for l in labels) and labels[-1].startswith('every structural option -> default')
    assert rows[-1][1]['efficiency'] == 8.0 and rows[1][2] == '' or rows[1][2] == 'already the default'


def collections_first(results):
    return next(iter(results))


# ==================================================================================================================================
# the downloader against a damaged local file (SM1-8)
# ==================================================================================================================================
def _flip_one_byte(path):
    b = bytearray(open(path, 'rb').read())
    b[len(b) // 2] ^= 0xFF                                                               # the size stays the same
    open(path, 'wb').write(bytes(b))


def test_fetch_repairs_a_same_size_corrupted_local_file(tmp_path, local_payload):
    """A flipped byte keeps the size. The recorded digest then no longer matches the file, but the server still delivers the recorded content: the local file is
    damaged, so it is replaced - not reported as 'the server content changed'."""
    F = _load_fetch()
    payload, make = local_payload
    sv = make()
    dest = str(tmp_path / 'file.bin')
    open(dest, 'wb').write(payload)
    rec = {'sha256': hashlib.sha256(payload).hexdigest()}
    assert F.ensure_file(sv.url, dest, len(payload), rec, log=lambda *a: None) == (rec['sha256'], 'verified') and sv.starts == []        # intact: nothing downloaded
    _flip_one_byte(dest)
    assert F.sha256_of(dest) != rec['sha256'] and os.path.getsize(dest) == len(payload)
    digest, action = F.ensure_file(sv.url, dest, len(payload), rec, log=lambda *a: None)
    assert (digest, action) == (rec['sha256'], 'repaired')
    assert open(dest, 'rb').read() == payload and not os.path.exists(dest + '.part') and sv.starts == [0]                                  # one complete fresh GET
    # a stale, possibly damaged '.part' of an earlier run is not trusted for the comparison
    _flip_one_byte(dest)
    open(dest + '.part', 'wb').write(payload[:1000] + b'garbage')
    assert F.ensure_file(sv.url, dest, len(payload), rec, log=lambda *a: None)[1] == 'repaired' and open(dest, 'rb').read() == payload


def test_fetch_stops_with_server_content_changed_only_when_a_fresh_download_also_differs(tmp_path, local_payload):
    F = _load_fetch()
    payload, make = local_payload
    other = bytes(reversed(payload))                                                      # the server now delivers other content of the same size
    sv = make()
    sv.payload = other
    dest = str(tmp_path / 'file.bin')
    rec = {'sha256': hashlib.sha256(payload).hexdigest()}
    bad = bytearray(payload)
    bad[7] ^= 1
    open(dest, 'wb').write(bytes(bad))                                                    # damaged local file AND changed server content
    with pytest.raises(F.ContentChanged, match='server content changed'):
        F.ensure_file(sv.url, dest, len(payload), rec, log=lambda *a: None)
    assert open(dest, 'rb').read() == bytes(bad) and not os.path.exists(dest + '.part')   # nothing was replaced, no leftover
    # the local file equals the server's new content, the record is old: also a real change
    open(dest, 'wb').write(other)
    with pytest.raises(F.ContentChanged):
        F.ensure_file(sv.url, dest, len(payload), rec, log=lambda *a: None)
    # --accept-new records the digest of the FRESH download
    digest, action = F.ensure_file(sv.url, dest, len(payload), rec, accept_new=True, log=lambda *a: None)
    assert (digest, action) == (hashlib.sha256(other).hexdigest(), 'accepted-new') and open(dest, 'rb').read() == other
    # no local file at all: the fresh download is compared with the record in the same way
    os.remove(dest)
    with pytest.raises(F.ContentChanged):
        F.ensure_file(sv.url, dest, len(payload), rec, log=lambda *a: None)
    assert not os.path.exists(dest)


def test_fetch_accept_new_never_records_the_digest_of_a_damaged_local_file(tmp_path, local_payload):
    F = _load_fetch()
    payload, make = local_payload
    sv = make()                                                                           # the server still delivers the recorded content
    dest = str(tmp_path / 'file.bin')
    open(dest, 'wb').write(payload)
    _flip_one_byte(dest)
    damaged = F.sha256_of(dest)
    rec = {'sha256': hashlib.sha256(payload).hexdigest()}
    digest, action = F.ensure_file(sv.url, dest, len(payload), rec, accept_new=True, log=lambda *a: None)
    assert digest == rec['sha256'] != damaged and action == 'repaired' and open(dest, 'rb').read() == payload
    # a local file of the wrong size is replaced by a fresh download; a file without a record is taken as it is
    open(dest, 'wb').write(payload[:10])
    assert F.ensure_file(sv.url, dest, len(payload), rec, log=lambda *a: None) == (rec['sha256'], 'downloaded')
    assert F.ensure_file(sv.url, dest, len(payload), None, log=lambda *a: None) == (rec['sha256'], 'recorded')


def _patch_fetch_to_local(F, monkeypatch, url, kind, rel, size):
    monkeypatch.setattr(F, 'files_for', lambda name: {kind: (rel, url)})
    monkeypatch.setattr(F, 'EXPECTED_BYTES', {'MT0.N1.LD': {kind: size}})


def test_fetch_main_repairs_and_keeps_the_recorded_digest(tmp_path, local_payload, monkeypatch, capsys):
    F = _load_fetch()
    payload, make = local_payload
    sv = make()
    _patch_fetch_to_local(F, monkeypatch, sv.url, 'activations', 'MT0.N1.LD/activations.csv', len(payload))
    data, ck = tmp_path / 'data', tmp_path / 'checksums.json'
    dest = data / 'MT0.N1.LD' / 'activations.csv'
    dest.parent.mkdir(parents=True)
    dest.write_bytes(payload)
    _flip_one_byte(str(dest))
    import json as _json
    ck.write_text(_json.dumps({'MT0.N1.LD/activations.csv': {'sha256': hashlib.sha256(payload).hexdigest(), 'bytes': len(payload)}}))
    F.main(['--sets', 'MT0.N1.LD', '--data-dir', str(data), '--checksums', str(ck), '--accept-new'])           # even with --accept-new
    assert dest.read_bytes() == payload
    rec = _json.loads(ck.read_text())['MT0.N1.LD/activations.csv']
    assert rec['sha256'] == hashlib.sha256(payload).hexdigest()
    assert 'repaired' in capsys.readouterr().out


def test_fetch_zip_crc_error_is_a_clean_exit_and_records_nothing(tmp_path, local_payload, monkeypatch):
    F = _load_fetch()
    zpath = tmp_path / 'src.zip'
    with zipfile.ZipFile(zpath, 'w', zipfile.ZIP_STORED) as z:
        z.writestr('stack.tif', bytes(range(256)) * 80)
    raw = bytearray(zpath.read_bytes())
    raw[200] ^= 0xFF                                                                      # inside the stored member data: the CRC no longer matches
    payload, make = local_payload
    sv = make()
    sv.payload = bytes(raw)
    _patch_fetch_to_local(F, monkeypatch, sv.url, 'stack', 'MT0.N1.LD/stack.zip', len(raw))
    (tmp_path / 'bad.zip').write_bytes(bytes(raw))
    with pytest.raises(RuntimeError, match='CRC error'):
        F.zip_members(str(tmp_path / 'bad.zip'))
    data, ck = tmp_path / 'data', tmp_path / 'checksums.json'
    with pytest.raises(SystemExit) as e:
        F.main(['--sets', 'MT0.N1.LD', '--data-dir', str(data), '--checksums', str(ck)])
    assert 'CRC error' in str(e.value) and 'no digest recorded' in str(e.value)
    assert not ck.exists() and not (data / 'MT0.N1.LD' / 'stack.zip').exists()          # nothing recorded, the damaged file removed so that the next run downloads it anew
    with open(tmp_path / 'junk.zip', 'wb') as f:
        f.write(b'this is not a zip file at all')
    with pytest.raises(RuntimeError, match='not a valid zip'):
        F.zip_members(str(tmp_path / 'junk.zip'))


# ==================================================================================================================================
# Picasso's options end to end (needs Picasso: set SMLM_PICASSO_PYTHON to an interpreter that has it, else skipped)
# ==================================================================================================================================
_PICASSO_PY = os.environ.get('SMLM_PICASSO_PYTHON')
needs_picasso = pytest.mark.skipif(not _PICASSO_PY and not C._picasso_importable(), reason='Picasso is not available (set SMLM_PICASSO_PYTHON)')


def _synthetic_stack(positions_nm, photons=9000.0, sigma_px=1.3, n_frames=3, seed=0):
    """Frames in DN with Gaussian spots (pixel-integrated, sigma 1.3 px) of ``photons`` at the given positions (nm, corner convention) in every frame, background
    100 photons (700 DN) and the EMCCD gain of the harness (6 DN per photon), with a little Gaussian noise."""
    rng = np.random.default_rng(seed)
    ax = np.arange(64) + 0.5
    img = np.zeros((64, 64))
    for x, y in positions_nm:
        px = 0.5 * (np.vectorize(math.erf)((ax + 0.5 - 0.5 - x / 100.0 + 0.5) / (math.sqrt(2) * sigma_px)) - np.vectorize(math.erf)((ax - 0.5 - 0.5 - x / 100.0 + 0.5) / (math.sqrt(2) * sigma_px)))
        py = 0.5 * (np.vectorize(math.erf)((ax + 0.5 - 0.5 - y / 100.0 + 0.5) / (math.sqrt(2) * sigma_px)) - np.vectorize(math.erf)((ax - 0.5 - 0.5 - y / 100.0 + 0.5) / (math.sqrt(2) * sigma_px)))
        img += photons * np.outer(py, px)
    frames = C.BASELINE_DN + C.DN_PER_PHOTON * (100.0 + img)[None] + rng.normal(0, 20.0, (n_frames, 64, 64))
    return np.round(np.clip(frames, 0, 65535))


@needs_picasso
def test_picasso_options_padding_filters_and_session_end_to_end():
    interior, corner = (3213.0, 3187.0), (120.0, 140.0)                 # a spot in the middle and one 1.2 px from the corner (inside Picasso's blind border zone)
    frames = _synthetic_stack([interior, corner])
    gt = _gt([(f, *interior, 9000.0) for f in range(3)] + [(f, *corner, 9000.0) for f in range(3)])
    base = dict(min_net_gradient=2000.0)
    with C.PicassoSession(frames, python=_PICASSO_PY) as s:              # one worker for all runs (Picasso's JIT compilation is paid once)
        d0 = s.run(base)
        r0 = C.score(d0, gt, [0, 1, 2])
        assert r0['tp'] == 3 and r0['fn'] == 3                             # the corner spot is never identified by the default pipeline
        d1 = s.run(dict(base, pad=4))
        r1 = C.score(d1, gt, [0, 1, 2])
        assert r1['tp'] == 6 and r1['fp'] == 0                            # median padding removes the blind zone; positions are shifted back to the corner convention
        assert r1['rmse_nm'] < 60.0
        assert all(d[:, 0].min() >= 0 and d[:, 1].min() >= 0 for d in d1 if len(d))            # centres outside the field are dropped
        # box options: an identification box of 5 with a fit box of 9; the Gaussian filter and the fit method are accepted
        for extra in (dict(box=5, fit_box=9, pad=4), dict(box=7, gauss_sigma=1.5, pad=4, min_net_gradient=500.0), dict(fit='gaussmle')):
            r = C.score(s.run(dict(base, **extra)), gt, [0, 1, 2])
            assert r['tp'] >= 3 and r['rmse_nm'] < 80.0, extra
        # post-filters are rows of the table: a photon minimum above the spots' photons removes everything, a width window around 130 nm keeps them
        assert sum(len(d) for d in s.run(dict(base, pad=4, min_photons=1e9))) == 0
        assert sum(len(d) for d in s.run(dict(base, pad=4, min_sigma_nm=80.0, max_sigma_nm=250.0))) == 6
        assert sum(len(d) for d in s.run(dict(base, pad=4, max_sigma_nm=50.0))) == 0
        # one low-threshold run filtered equals a direct run
        low = s.run(dict(base, pad=4, min_net_gradient=300.0))
        direct = s.run(dict(base, pad=4, min_net_gradient=3000.0, min_photons=100.0))
        filt = C.filter_detections(low, dict(min_net_gradient=3000.0, min_photons=100.0))
        for a, b in zip(direct, filt):
            assert a.shape == b.shape and np.allclose(a, b, equal_nan=True)
        info = {}
        again = s.run(dict(base, pad=4), info)
        assert info['n_frames'] == 3 and info['n_loc'] == sum(len(d) for d in again)
    assert all(a.shape == b.shape and np.allclose(a, b) for a, b in zip(again, d1))
    # a process per run (run_picasso) gives the same detections as the session
    one = C.run_picasso(frames, dict(base, pad=4), python=_PICASSO_PY)
    assert all(a.shape == b.shape and np.array_equal(a, b) for a, b in zip(one, d1))
    with pytest.raises(RuntimeError, match='failed'):
        with C.PicassoSession(frames, python=_PICASSO_PY) as s:
            s.run(dict(base, fit='no-such-fit'))


# ==================================================================================================================================
# the calibration script: quiet frames, halo profile, estimators, windows (SM1-3, SM1-4, SM1-5, SM1-10)
# ==================================================================================================================================
def _load_calibration():
    spec = importlib.util.spec_from_file_location('exp_smlm_calibration', os.path.join(os.path.dirname(__file__), '..', 'experiments', 'exp_smlm_calibration.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_calibration_quiet_frames_and_the_camera_constants_from_them():
    K = _load_calibration()
    gt = _gt([(5, 1., 1., 1.), (6, 1., 1., 1.), (12, 1., 1., 1.)])
    q = K.quiet_frames(gt, 16)
    assert list(np.flatnonzero(q)) == [0, 1, 2, 3, 8, 9, 10, 14, 15]                      # no row in the frame or in either neighbour; the ends have one neighbour only
    assert list(K.calibration_frames(2000)[[0, 999, 1000, -1]]) == [0, 999, 1500, 1999]    # the test frames 1000-1499 are left out
    rng = np.random.default_rng(1)
    m = 200
    dn = np.round(700.0 + 90.0 * rng.standard_normal((m, 64, 64)))                       # photons: mean 100, sd 15, var 225, var / mean 2.25
    st = dn.astype(np.uint16)
    qs = K.quiet_stats(st, np.arange(m))
    assert qs['mean'] == pytest.approx(100.0, abs=4 * qs['se_mean']) and qs['se_mean'] < 0.2
    assert qs['var'] == pytest.approx(225.0, abs=4 * qs['se_var']) and qs['var_over_mean'] == pytest.approx(2.25, rel=0.02)
    assert qs['frame_sd'] == pytest.approx(qs['expect_sd'], rel=0.15)                    # independent pixels: no common mode
    assert all(abs(r[2] - qs['mean']) < 5 * r[3] for r in qs['radial'])                    # flat in space
    common = np.round(700.0 + 90.0 * rng.standard_normal((m, 64, 64)) + 6.0 * rng.standard_normal(m)[:, None, None] * 2.0).astype(np.uint16)
    assert K.quiet_stats(common, np.arange(m))['frame_sd'] > 3 * qs['expect_sd']             # a frame-level common mode would show up in the sd of the frame means


def _halo_psf(r, core=0.7, sigma=1.3, cut=25.0):
    """Radial profile (fraction of the emitter's photons per pixel): a Gaussian core with ``core`` of the photons plus a halo ~ 1 / (1 + (r / 5)^2) cut off at ``cut`` px."""
    yy, xx = np.mgrid[-60:61, -60:61]
    rr = np.hypot(xx, yy)
    h = np.where(rr < cut, 1.0 / (1.0 + (rr / 5.0) ** 2), 0.0)
    norm = h.sum()
    return core * np.exp(-0.5 * r ** 2 / sigma ** 2) / (2 * np.pi * sigma ** 2) + np.where(r < cut, (1 - core) / norm / (1.0 + (r / 5.0) ** 2), 0.0)


def test_halo_profile_recovers_a_known_radial_profile_and_its_cutoff():
    K = _load_calibration()
    rng = np.random.default_rng(2)
    yy, xx = np.mgrid[0:64, 0:64]
    n, base = 60, 100.0
    rows, frames = [], []
    for f in range(n):
        x, y = rng.uniform(1700, 4700, 2)
        P_ = 5000.0
        r = np.hypot(xx + 0.5 - x / 100.0, yy + 0.5 - y / 100.0)
        frames.append(np.round(C.BASELINE_DN + C.DN_PER_PHOTON * (base + P_ * _halo_psf(r))))
        rows.append((f, x, y, P_))
    stack = np.array(frames).astype(np.uint16)
    gt = _gt(rows)
    gt['z_nm'] = 0.0
    hp = K.halo_profile(stack, gt, np.arange(n), base, 1000.0)['all']
    # the cumulative fraction inside a radius of the profile above (computed on a fine grid of an infinite plane)
    yy2, xx2 = np.mgrid[-60:61, -60:61]
    rr = np.hypot(xx2, yy2)
    truth = [float(_halo_psf(rr)[rr < R].sum()) for R in K.HALO_RADII]
    assert hp['n'] == n and np.allclose(hp['cum'], truth, atol=0.04)
    assert hp['cum'][3] > hp['cum'][1] + 0.15 and hp['cum'][3] < 1.02                      # a halo beyond 5 px, and the profile ends at 25 px
    k = list(K.HALO_EDGES).index(26.0)                                                      # annulus 26-28 px: nothing beyond the cut-off
    assert abs(hp['prof_per_1000'][k]) < 0.02 and hp['prof_per_1000'][list(K.HALO_EDGES).index(12.0)] == pytest.approx(0.15, abs=0.03)       # the halo of 0.3 of the photons, 0.15 per 1000 photons per pixel at 13 px
    assert hp['image_total'] == pytest.approx(1.0, abs=0.06)                                # the whole 64 x 64 image holds the photons inside the frame


def test_calibration_offset_estimators_and_window_fits_recover_known_values():
    K = _load_calibration()
    rng = np.random.default_rng(3)
    # a symmetric Gaussian centred (+0.10, -0.05) px from the "ground truth": the fit of the stacked samples finds the centre and the width
    dx, dy = np.meshgrid(np.arange(-7, 8) + 0.37, np.arange(-7, 8) - 0.22)
    w = 0.02 * np.exp(-0.5 * (((dx - 0.10) ** 2 + (dy + 0.05) ** 2) / 1.3 ** 2)) + rng.normal(0, 1e-4, dx.shape)
    x0, y0, s = K.gaussian_profile_fit(dx.ravel(), dy.ravel(), w.ravel())
    assert (x0, y0) == pytest.approx((10.0, -5.0), abs=2.0) and s == pytest.approx(1.3, abs=0.02)
    # the iterated centroid of a spot 30 nm right and 20 nm above the ground-truth position, and the stacked samples around that position
    ax = np.arange(64) + 0.5
    gx, gy, sx, sy = 3210.0, 3190.0, 3240.0, 3170.0                                          # truth, and the spot centre in nm
    img = C.BASELINE_DN + C.DN_PER_PHOTON * (100.0 + 4000.0 * np.outer(np.exp(-0.5 * (ax - sy / 100) ** 2 / 1.3 ** 2), np.exp(-0.5 * (ax - sx / 100) ** 2 / 1.3 ** 2)) / (2 * np.pi * 1.3 ** 2))
    frames = np.round(img)[None].astype(np.uint16)
    gt = _gt([(0, gx, gy, 4000.0)])
    cen = K.per_emitter_centroids(frames, gt, np.array([0]), 100.0)
    assert cen[0] == pytest.approx((30.0, -20.0), abs=3.0)
    E_, DX, DY, W = K.stacked_samples(frames, gt, np.array([0]), 100.0, half=7)
    r_ = np.hypot(DX, DY)
    w_ = W[r_ < 6]
    assert 100 * (w_ * DX[r_ < 6]).sum() / w_.sum() == pytest.approx(30.0, abs=3.0) and w_.sum() == pytest.approx(1.0, abs=0.05)      # the whole PSF holds the photons
    # integrated-Gaussian fits in windows of every size recover width and centre of a noise-free spot
    for half in K.WINDOW_HALVES:
        yy, xx = K._grid(half)
        patch = 700.0 + 6.0 * 5000.0 * np.outer(K._int_profile(yy[:, 0], -0.3, 1.4), K._int_profile(xx[0], 0.2, 1.4))
        fit = K.fit_patch(patch, 'int', half)
        assert fit is not None and fit[0][3] == pytest.approx(1.4, abs=0.01) and (fit[0][1], fit[0][2]) == pytest.approx((0.2, -0.3), abs=0.01), half
    # the illumination table from a summary file
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'stats_activation_sumary.csv')
        open(path, 'w').write('feature,value \nCenter of power illumination X (nm), 3000.0\nCenter of power illumination X (nm), 3000.0\nDistance at half maximum illumination (nm), 3500.0\n')
        gt2 = np.zeros(300, dtype=C.GT_DTYPE)
        gt2['x_nm'] = np.r_[np.full(150, 3000.0), np.full(150, 5500.0)]
        gt2['y_nm'] = 3000.0
        gt2['photons'] = np.r_[np.full(150, 5000.0), np.full(150, 3000.0)]
        cen2, half2, rows2 = K.illumination_profile(gt2, path)
        assert cen2 == (3000.0, 3000.0) and half2 == 3500.0 and [r[2] for r in rows2] == [150, 150] and rows2[0][3] > rows2[1][3]
        assert K.illumination_profile(gt2, os.path.join(d, 'missing.csv')) is None
