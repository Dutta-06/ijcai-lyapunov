"""Tests of the export and ingest scripts for the GPU methods (experiments/smlm_decode_export.py, smlm_decode_ingest.py): a TIFF round trip through the harness's loader,
the convention check, and the dev-chosen threshold / test scoring on detections made from a synthetic ground truth."""
import os
import sys

import numpy as np
import pytest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'experiments')
sys.path.insert(0, HERE)
import smlm_common as C            # noqa: E402
import smlm_decode_export as EX    # noqa: E402
import smlm_decode_ingest as IN    # noqa: E402


def test_tiff_roundtrip_through_the_harness_loader(tmp_path):
    rng = np.random.default_rng(0)
    fr = rng.integers(0, 3000, (5, 64, 64)).astype(float)
    EX.write_tiff(str(tmp_path / 'x.tif'), fr)
    back = C.load_frames(str(tmp_path / 'x.tif'))
    assert back.shape == (5, 64, 64) and np.array_equal(back, fr)
    with pytest.raises(ValueError):
        EX.write_tiff(str(tmp_path / 'y.tif'), fr + 0.5)


def fake_gt(n_frames, per_frame=2, seed=1):
    rng = np.random.default_rng(seed)
    rows = []
    for f in range(n_frames):
        for _ in range(per_frame):
            rows.append((len(rows) + 1, f, rng.uniform(300, 6100), rng.uniform(300, 6100), 0.0, 3000.0))
    return np.array(rows, dtype=C.GT_DTYPE)


def csv_from_gt(gt, a, b, rng, offset=(0.0, 0.0), miss=0.2, n_false=1, jitter=20.0):
    """Candidates of a 'detector': each true emitter with probability 1 - miss (prob high), jittered, and false candidates with low prob; frames relative to a."""
    rows = []
    for g in gt[(gt['frame'] >= a) & (gt['frame'] < b)]:
        if rng.uniform() > miss:
            rows.append((g['frame'] - a, g['x_nm'] + offset[0] + jitter * rng.standard_normal(), g['y_nm'] + offset[1] + jitter * rng.standard_normal(), 3000.0, rng.uniform(0.7, 1.0)))
    for f in range(b - a):
        for _ in range(n_false):
            rows.append((f, rng.uniform(0, 6400), rng.uniform(0, 6400), 500.0, rng.uniform(0.05, 0.4)))
    return np.array(rows, dtype=[('frame', float), ('x_nm', float), ('y_nm', float), ('photons', float), ('prob', float)])


def test_ingest_chooses_thresholds_on_dev_removes_the_offset_and_scores_test():
    gt = fake_gt(1600)
    rng = np.random.default_rng(3)
    dev = csv_from_gt(gt, 0, 1000, rng, offset=(12.0, -9.0))
    test = csv_from_gt(gt, 1000, 1500, rng, offset=(12.0, -9.0))
    res = IN.evaluate('MT0.N1.LD', dev, test, gt)
    for op in ('E', 'J'):
        o = res[op]
        assert 0.4 <= o['thr'] <= 0.7                                  # the false candidates have prob < 0.4, the true ones > 0.7
        assert abs(o['off'][0] - 12.0) < 3 and abs(o['off'][1] + 9.0) < 3
        assert o['test']['precision'] > 0.98 and 0.7 < o['test']['recall'] < 0.9 and o['test']['rmse_nm'] < 40
    assert res['E']['test']['rmse_nm'] < 35                            # the offset is removed: the jitter 20 nm per axis gives ~28 nm
    # a pixel-convention slip of +50 nm is flagged by the bias before the correction
    dev_slip = csv_from_gt(gt, 0, 1000, np.random.default_rng(4), offset=(50.0, 50.0))
    test_slip = csv_from_gt(gt, 1000, 1500, np.random.default_rng(5), offset=(50.0, 50.0))
    r = IN.evaluate('MT0.N1.LD', dev_slip, test_slip, gt)
    assert max(abs(r['E']['uncorrected_bias'])) > 30


def test_ingest_validates_columns_and_frame_range(tmp_path):
    p = tmp_path / 'bad.csv'
    p.write_text('frame,x_nm,y_nm\n0,1,2\n')
    with pytest.raises(ValueError):
        IN.read_csv(str(p))
    rows = np.array([(5000, 10.0, 10.0, 1.0, 0.9)], dtype=[('frame', float), ('x_nm', float), ('y_nm', float), ('photons', float), ('prob', float)])
    with pytest.raises(ValueError):
        IN.per_frame(rows, 100, 0.5)
    off, flat = IN.pack([np.array([[1.0, 2.0, 3.0]]), np.zeros((0, 3)), np.array([[4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])])
    assert list(off) == [0, 1, 1, 3] and flat.shape == (3, 4) and flat[2, 2] == 9.0
