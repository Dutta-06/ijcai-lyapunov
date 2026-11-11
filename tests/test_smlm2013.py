"""Tests of the loader and protocol of the 2013 challenge's training sets (`experiments/smlm2013.py`) and of the helpers of `experiments/exp_smlm2013.py`:
synthetic tiny files with known answers (per-frame TIFF zips in both layouts, both byte orders, the ground-truth text, the field of the set, the convention checks, tiles, the noise law,
the width fit, the raw-file format) and, when the downloaded data are present (`data/smlm/fetch_2013.py`, $SMLM_DATA_2013), facts of the real sets: counts, ranges, the frame base,
one-frame activations, and a Picasso run on 5 frames of each set matched to the ground truth within tolerance (this one needs SMLM_PICASSO_PYTHON, else it is skipped)."""
import json
import os
import struct
import sys
import zipfile

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'experiments'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import smlm_common as C
import smlm2013 as S


# ------------------------------------------------------------------ helpers: single-page 16-bit TIFFs written by hand (Pillow cannot write the 16-bit modes here)
def _tiff16(arr, big_endian=True):
    arr = np.asarray(arr)
    h, w = arr.shape
    e = '>' if big_endian else '<'
    data = arr.astype(e + 'u2').tobytes()
    ent = [(256, 4, w), (257, 4, h), (258, 3, 16), (259, 3, 1), (262, 3, 1), (273, 4, 0), (277, 3, 1), (278, 3, h), (279, 4, len(data))]
    ifd_len = 2 + 12 * len(ent) + 4
    ent[5] = (273, 4, 8 + ifd_len)
    blob = bytearray((b'MM\x00\x2a' if big_endian else b'II\x2a\x00') + struct.pack(e + 'I', 8) + struct.pack(e + 'H', len(ent)))
    for tag, typ, val in ent:
        blob += struct.pack(e + 'HHI', tag, typ, 1) + (struct.pack(e + 'HH', val, 0) if typ == 3 else struct.pack(e + 'I', val))
    blob += struct.pack(e + 'I', 0) + data
    return bytes(blob)


def _zip_of_frames(path, stack, prefix='', big_endian=True, first=1, extra=True):
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        if prefix:
            z.writestr(prefix, b'')                                                         # a directory member
        if extra:
            z.writestr(prefix + 'readme.txt', b'not a frame')
        for k in range(len(stack)):
            z.writestr('%s%05d.tif' % (prefix, k + first), _tiff16(stack[k], big_endian))
    return str(path)


GT_TEXT = ('Ground-truth,frame,xnano,ynano,znano,intensity\n'
           '    1,  3,  1000.50,  2000.25,  30.0,  1500.5\n'
           '    2,  1,   100.00,   200.00,  10.0,   900.0\n'
           '    3,  3,  5000.00,  3000.00,  50.0,  2500.0\n')


# ------------------------------------------------------------------ loader
@pytest.mark.parametrize('prefix', ['', 'sequence/'])
@pytest.mark.parametrize('big', [True, False])
def test_load_frames_per_frame_tiffs_are_zero_based_in_both_layouts_and_byte_orders(tmp_path, prefix, big):
    rng = np.random.default_rng(1)
    stack = rng.integers(100, 3000, size=(6, 5, 7))
    p = _zip_of_frames(tmp_path / 'seq.zip', stack, prefix, big)
    a = S.load_frames(p, 0, 6)
    assert a.dtype == np.float64 and a.shape == (6, 5, 7) and np.array_equal(a, stack)
    assert np.array_equal(S.load_frames(p, 2, 4), stack[2:4])                              # frame k is the TIFF numbered k + 1; stop exclusive
    assert np.array_equal(S.load_frames(p, 3), stack[3:])                                  # stop=None: to the end
    assert S.load_frames(p, 4, 4).shape[0] == 0
    assert S.n_frames_in(p) == 6
    with pytest.raises(IndexError):
        S.load_frames(p, 0, 7)                                                             # no silent truncation
    with pytest.raises(ValueError):
        S.load_frames(p, 5, 2)


def test_load_frames_refuses_gaps_and_a_wrong_first_number(tmp_path):
    stack = np.arange(4 * 3 * 3).reshape(4, 3, 3) + 100
    p = tmp_path / 'gap.zip'
    with zipfile.ZipFile(p, 'w') as z:
        for k in (1, 2, 4):
            z.writestr('%05d.tif' % k, _tiff16(stack[k - 1]))
    with pytest.raises(ValueError):
        S.load_frames(str(p), 0, 2)
    with pytest.raises(ValueError):
        S.n_frames_in(str(p))
    q = _zip_of_frames(tmp_path / 'from0.zip', stack, first=0)
    with pytest.raises(ValueError):
        S.load_frames(q, 0, 2)                                                             # numbered 0..3: the files of these sets start at 1


def test_load_frames_refuses_frames_of_different_shapes(tmp_path):
    p = tmp_path / 'mixed.zip'
    with zipfile.ZipFile(p, 'w') as z:
        z.writestr('00001.tif', _tiff16(np.ones((3, 3))))
        z.writestr('00002.tif', _tiff16(np.ones((3, 4))))
    with pytest.raises(ValueError):
        S.load_frames(str(p), 0, 2)


def test_to_photons_uses_the_pages_conversion():
    assert S.BASELINE_DN == 100.0 and S.DN_PER_PHOTON == 1.0
    assert np.allclose(S.to_photons(np.array([[100.0, 101.5], [600.0, 90.0]])), [[0.0, 1.5], [500.0, -10.0]])
    assert S.to_photons(np.array([300], dtype=np.uint16))[0] == 200.0                      # integers are converted, not wrapped


# ------------------------------------------------------------------ ground truth
def test_parse_gt_equals_the_harness_loader_and_is_zero_based_and_sorted(tmp_path):
    g = S.parse_gt(GT_TEXT)
    csv = tmp_path / 'activation.csv'
    csv.write_text(GT_TEXT)
    h = C.load_gt(str(csv))
    assert g.dtype == h.dtype and all(np.array_equal(g[k], h[k]) for k in g.dtype.names)    # the same parse as the 2016 harness
    assert list(g['frame']) == [0, 2, 2] and list(g['id']) == [2, 1, 3]                      # file frame 1 is stack index 0; sorted by (frame, id)
    assert g['x_nm'][1] == 1000.5 and g['y_nm'][1] == 2000.25 and g['z_nm'][1] == 30.0 and g['photons'][1] == 1500.5
    assert S.parse_gt(GT_TEXT, frame_base=0)['frame'].min() == 1
    with pytest.raises(ValueError):
        S.parse_gt('43.3, 839.8, 340.7\n39.0, 850.1, 331.8\n')                              # a position list is not a ground truth
    assert len(S.parse_gt('Ground-truth,frame,xnano,ynano,znano,intensity\n')) == 0


def test_load_gt_from_zip_finds_activation_csv_wherever_it_sits(tmp_path):
    for member in ('activation.csv', 'fluorophores/activation.csv'):                        # the layouts of the tubulin and of the bundle sets
        p = tmp_path / ('gt_%d.zip' % len(member))
        with zipfile.ZipFile(p, 'w') as z:
            z.writestr(member, GT_TEXT)
            z.writestr(member.replace('activation', 'frames/00003').replace('.csv', '.csv'), 'Ground-truth,frame,xnano,ynano,znano,intensity\n    1,  3,  1.0, 2.0, 3.0, 4.0\n')
        g = S.load_gt(str(p))
        assert len(g) == 3 and list(g['frame']) == [0, 2, 2]
    bad = tmp_path / 'two.zip'
    with zipfile.ZipFile(bad, 'w') as z:
        z.writestr('a/activation.csv', GT_TEXT)
        z.writestr('b/activation.csv', GT_TEXT)
    with pytest.raises(ValueError):
        S.load_gt(str(bad))                                                                # ambiguous
    none = tmp_path / 'none.zip'
    with zipfile.ZipFile(none, 'w') as z:
        z.writestr('oracle/projection.tif', b'x')
    with pytest.raises(ValueError):
        S.load_gt(str(none))                                                               # the oracle archive holds no localizations


def test_gt_frame_files_collect_the_per_frame_csvs(tmp_path):
    p = tmp_path / 'f.zip'
    with zipfile.ZipFile(p, 'w') as z:
        z.writestr('activation.csv', GT_TEXT)
        z.writestr('frames/00001.csv', 'Ground-truth,frame,xnano,ynano,znano,intensity\n 1, 1, 100.00, 200.00, 10.0, 900.0\n')
        z.writestr('frames/00003.csv', 'Ground-truth,frame,xnano,ynano,znano,intensity\n 1, 3, 1000.50, 2000.25, 30.0, 1500.5\n 2, 3, 5000.00, 3000.00, 50.0, 2500.0\n')
    a, b = S.load_gt(str(p)), S.gt_frame_files(str(p))
    assert len(a) == len(b) == 3
    assert all(np.array_equal(np.sort(a[k]), np.sort(b[k])) for k in ('frame', 'x_nm', 'y_nm', 'z_nm', 'photons'))
    assert len(S.gt_in_frames(a, 0, 1)) == 1 and len(S.gt_in_frames(a, 1, 3)) == 2


# ------------------------------------------------------------------ field and scoring
def test_field_and_spec_of_the_four_sets():
    assert S.field_nm('Tubulins_I') == 38400.0 and S.field_nm('Tubulins_II') == 38400.0
    assert S.field_nm('Bundled_Tubes_High_Density') == 6400.0 and S.field_nm('Bundled_Tubes_Long_Sequence') == 6400.0
    for name in S.SETS:
        sp = S.spec(name)
        assert sp['size'] * sp['pixel_nm'] == S.field_nm(name)
    with pytest.raises(KeyError):
        S.spec('MT0.N1.LD')


def test_score_uses_the_field_of_the_set_for_the_border_and_restores_the_harness_field():
    gt = np.zeros(3, dtype=C.GT_DTYPE)
    gt['frame'] = 0
    gt['x_nm'] = [100.0, 20000.0, 38300.0]                                                  # in the border zone, mid-field, in the border zone of the far edge
    gt['y_nm'] = [20000.0, 20000.0, 20000.0]
    gt['photons'] = 1000.0
    dets = [np.c_[gt['x_nm'], gt['y_nm']]]
    r = S.score(dets, gt, [0], 'Tubulins_I', border_nm=450.0)
    assert (r['tp'], r['fp'], r['fn']) == (1, 0, 0)                                          # only the mid-field emitter is counted on the 38,400 nm field
    r64 = C.score(dets, gt, [0], border_nm=450.0)                                           # the harness's 6,400 nm field drops the mid-field one as well
    assert r64['tp'] == 0
    assert C.FIELD_NM == 6400.0
    with pytest.raises(RuntimeError):
        with S.field('Tubulins_I'):
            assert C.FIELD_NM == 38400.0
            raise RuntimeError
    assert C.FIELD_NM == 6400.0                                                             # restored after an exception too
    gb = gt.copy()
    gb['x_nm'] = [100.0, 3000.0, 6350.0]
    gb['y_nm'] = 3000.0
    db = [np.c_[gb['x_nm'], gb['y_nm']]]
    a, b = S.score(db, gb, [0], 'Bundled_Tubes_Long_Sequence', border_nm=450.0), C.score(db, gb, [0], border_nm=450.0)
    assert (a['tp'], a['fp'], a['fn']) == (b['tp'], b['fp'], b['fn']) == (1, 0, 0)          # on a 6,400 nm field it is the harness's score


def test_sigma_page_is_a_reading_of_the_page():
    s0 = 258.21 / (2 * np.sqrt(2 * np.log(2)))
    assert S.sigma_page_nm(0.0, 'Tubulins_I') == pytest.approx(s0)
    assert S.sigma_page_nm(500.0, 'Tubulins_I') == pytest.approx(2 * s0)                    # 2 x FWHM at 500 nm (tubulins)
    assert S.sigma_page_nm(300.0, 'Bundled_Tubes_Long_Sequence') == pytest.approx(2 * s0)   # and at 300 nm (bundles)
    assert S.sigma_page_nm(-150.0, 'Tubulins_I') == pytest.approx(S.sigma_page_nm(150.0, 'Tubulins_I'))


# ------------------------------------------------------------------ the convention checks
def _bump_frames(n, ny, nx, rows, pixel_nm, rng, amp=400.0):
    """n frames of noise around 140 DN with a bump (3 x 3 plateau of height amp) at each ground-truth row (frame, x_nm, y_nm) of the corner convention."""
    fr = rng.normal(140.0, 5.0, size=(n, ny, nx))
    for f, x, y in rows:
        ix, iy = int(np.floor(x / pixel_nm)), int(np.floor(y / pixel_nm))
        fr[f, iy - 1:iy + 2, ix - 1:ix + 2] += amp
    return fr


def test_excess_at_ground_truth_peaks_at_the_true_frame_offset_and_pixel():
    rng = np.random.default_rng(3)
    rows = [(f, 1250.0 + 700.0 * f, 2250.0 + 500.0 * ((f * 3) % 4)) for f in range(8)]
    fr = _bump_frames(8, 40, 40, rows, 150.0, rng)
    gt = np.zeros(len(rows), dtype=C.GT_DTYPE)
    gt['frame'] = [r[0] for r in rows]
    gt['x_nm'] = [r[1] for r in rows]
    gt['y_nm'] = [r[2] for r in rows]
    gt['photons'] = 1000.0
    ex = S.excess_at_ground_truth(fr, gt, 0, pixel_nm=150.0)
    assert ex[0][0] > 350.0 and all(abs(ex[o][0]) < 25.0 for o in (-2, -1, 1, 2))           # the right frame stands out sharply
    g2 = gt.copy()
    g2['frame'] = g2['frame'] + 1                                                           # a ground truth with the wrong frame base
    ex2 = S.excess_at_ground_truth(fr, g2, 0, pixel_nm=150.0)
    assert ex2[-1][0] > 350.0 and abs(ex2[0][0]) < 25.0                                      # ... peaks at offset -1 instead
    bright = S.excess_at_ground_truth(fr, gt, 0, pixel_nm=150.0, min_photons=2000.0)
    assert bright[0][1] == 0                                                                # the photon cut selects rows


def test_check_conventions_separates_corner_centre_axes_and_frames():
    rng = np.random.default_rng(5)
    p = 150.0
    gt = np.zeros(60, dtype=C.GT_DTYPE)
    gt['frame'] = np.repeat(np.arange(12), 5)
    gt['x_nm'] = rng.uniform(1000.0, 30000.0, 60)
    gt['y_nm'] = rng.uniform(1000.0, 30000.0, 60) * 0.5 + 8000.0
    gt['photons'] = 1000.0
    dets = []
    for f in range(12):
        g = gt[gt['frame'] == f]
        dets.append(np.c_[g['x_nm'] + rng.normal(0, 8, len(g)) + 3.0, g['y_nm'] + rng.normal(0, 8, len(g)) - 2.0])       # a symmetric estimator: ~ (+3, -2) nm off
    cv = S.check_conventions('Tubulins_I', dets, gt, np.arange(12))
    assert cv['corner']['tp'] == 60 and abs(cv['corner']['bias_x_nm'] - 3.0) < 4.0 and abs(cv['corner']['bias_y_nm'] + 2.0) < 4.0
    assert cv['centre']['bias_x_nm'] == pytest.approx(cv['corner']['bias_x_nm'] - 0.5 * p) and cv['centre']['bias_y_nm'] == pytest.approx(cv['corner']['bias_y_nm'] - 0.5 * p)
    assert cv['swapped']['tp'] < 10 and cv['flip_x']['tp'] < 10 and cv['flip_y']['tp'] < 10
    assert cv['frames'][0] == 60 and all(cv['frames'][o] < 10 for o in (-2, -1, 1, 2))      # the true frame base stands out


# ------------------------------------------------------------------ the declared protocol
def test_protocol_is_declared_consistently():
    for name in S.SETS:
        (a, b), (c, d) = S.PROTOCOL[name]['dev'], S.PROTOCOL[name]['test']
        n = S.SPEC[name]['n_frames']
        assert 0 <= a < b <= c < d <= n                                                     # development before test, no overlap, inside the set
        assert (d - c) >= 60 and (d - c) % 10 == 0                                          # whole blocks of the block bootstrap
    assert S.PROTOCOL['Tubulins_I'] == S.PROTOCOL['Tubulins_II']                            # the same frames for the two noise levels
    assert S.CAL_FROM == {'Bundled_Tubes_High_Density': 'Bundled_Tubes_Long_Sequence'}


# ------------------------------------------------------------------ the helpers of the experiment
@pytest.fixture(scope='module')
def X():
    import exp_smlm2013 as mod
    return mod


def test_declarations_of_the_experiment_are_consistent(X):
    import exp_smlm_iss as ISS
    import exp_smlm_picasso as PIC
    assert X.DEFAULT == dict(cfg='hp5h', kappa=5.0) and X.KAPPAS == ISS.KAPPAS and X.CONFIGS is ISS.CONFIGS
    assert X.SIGMA_REF == 1.4 and X.WIDTH_REL[0] * 1.4 == pytest.approx(1.0) and X.WIDTH_REL[1] * 1.4 == pytest.approx(3.2) and X.REX_REL * 1.4 == pytest.approx(2.0)
    assert [v * 1.4 for v in X.SMAX_REL] == pytest.approx(list(ISS.SMAX)) and X.SMAX_REL[-1] == pytest.approx(X.WIDTH_REL[1])         # the last width bound = no filter
    assert X.GRID_R1 == PIC.GRID_R1['N2'] and len(X.GRID_R1) == 20
    for k, g in X.FILTER_GRID.items():
        vals = [v for v in g if v is not None]
        assert vals == sorted(vals) and len(set(vals)) == len(vals)
    assert X.PIC_NG_LOW == 100.0 and set(X.FILTER_GRID) == set(C.POST_FILTER_KEYS)
    assert {k: X.RECIPE_OF[k] for k in S.SETS} == {'Tubulins_I': 'sparse', 'Tubulins_II': 'sparse', 'Bundled_Tubes_Long_Sequence': 'sparse', 'Bundled_Tubes_High_Density': 'crowded'}
    assert X.RECIPES['sparse']['snr_min'] == min(X.KAPPAS) and X.RECIPES['crowded']['snr_min'] is None and X.RECIPES['crowded']['cap'] == 400
    cam = X.camera('Tubulins_I')
    assert cam == dict(baseline=100.0, gain=1.0, sensitivity=1.0, qe=1.0, pixel_nm=150.0)
    C._merge_params(cam)                                                                    # accepted by the Picasso wrapper


def test_tile_grid_and_ownership_partition_the_axis(X):
    assert X.tile_grid(256, 64, 8) == [(0, 72, 0), (56, 136, 1), (120, 200, 2), (184, 256, 3)]
    assert X.tile_grid(64, 64, 8) == [(0, 64, 0)]                                           # a 64 px frame is one tile
    assert X.tile_grid(100, 64, 8) == [(0, 72, 0), (56, 100, 1)]                            # a shorter last core
    for n in (64, 100, 256):
        nt = len(X.tile_grid(n, 64, 8))
        x = np.linspace(-0.5, n - 0.5 - 1e-9, 4001)
        own = X.owner_index(x, 64, nt)
        assert own.min() == 0 and own.max() == nt - 1 and np.all(np.diff(own) >= 0)         # every position has exactly one owner, in order
        for k, (lo, hi, _) in enumerate(X.tile_grid(n, 64, 8)):
            xs = x[own == k]
            assert xs.min() >= lo - 0.5 - 1e-9 and xs.max() <= hi - 0.5 + 1e-9               # the owner's tile (core + margin) contains what it owns, with room to spare
    assert X.owner_index(63.49, 64, 4) == 0 and X.owner_index(63.51, 64, 4) == 1            # the core boundary is at x + 1/2 = 64


def test_merge_tiles_merges_close_detections_of_different_tiles_only(X):
    a = np.array([[10.0, 10.0, 500.0, 1.0], [20.0, 10.0, 900.0, 1.0], [20.4, 30.0, 100.0, 1.0]])
    b = np.array([[20.3, 10.2, 300.0, 1.0], [40.0, 40.0, 700.0, 1.0], [20.5, 30.2, 800.0, 1.0]])
    m = X.merge_tiles([a, b], 1.0)
    pts = {(round(r[0], 1), round(r[1], 1)): r[2] for r in m}
    assert len(m) == 4 and (20.0, 10.0) in pts and (20.3, 10.2) not in pts                   # the weaker of the pair (900 vs 300 photons) is dropped
    assert (20.5, 30.2) in pts and (20.4, 30.0) not in pts                                   # the other pair: the stronger one (800 vs 100) is kept
    c = np.array([[5.0, 5.0, 100.0, 1.0], [5.3, 5.1, 200.0, 1.0]])
    assert len(X.merge_tiles([c], 1.0)) == 2 and len(X.merge_tiles([c, np.zeros((0, 4))], 1.0)) == 2     # two detections of ONE tile are never merged
    assert X.merge_tiles([], 1.0).shape == (0, 4)


def test_pack_unpack_post_and_the_float32_width_bound(X):
    dets = [np.array([[1.5, 2.5, 100.0, 0.7], [3.0, 4.0, 50.0, 2.302]]), np.zeros((0, 4)), np.array([[9.0, 9.5, 10.0, 1.0]])]
    off, flat = X.pack(dets)
    assert list(off) == [0, 2, 2, 3] and flat.dtype == np.float32 and flat.shape == (3, 4)
    back = X.unpack(off, flat)
    assert [len(d) for d in back] == [2, 0, 1] and back[0].dtype == np.float64 and np.allclose(back[0], dets[0], rtol=1e-6)
    bound = 2.302176543
    stored = X.unpack(*X.pack([np.array([[1.0, 1.0, 5.0, bound]])]))
    assert stored[0][0, 3] != bound                                                         # float32 rounding moves the width ...
    assert len(X.post(stored, bound)[0]) == 1                                                # ... but a detection AT the bound survives the "no filter" bound
    assert len(X.post(back, 1.0)[0]) == 1 and len(X.post(back, 0.5)[0]) == 0


def test_save_npz_is_deterministic_and_roundtrips(X, tmp_path):
    arrays = {'a': np.arange(10.0), 'b|c': np.array('text'), 'gt': np.zeros(3, dtype=C.GT_DTYPE), 'off': np.arange(4, dtype=np.int64)}
    arrays['gt']['x_nm'] = [1.0, 2.0, 3.0]
    p1, p2 = str(tmp_path / 'one.npz'), str(tmp_path / 'two.npz')
    X.save_npz(p1, arrays)
    X.save_npz(p2, dict(reversed(list(arrays.items()))))                                    # the order of the keys does not matter
    assert open(p1, 'rb').read() == open(p2, 'rb').read()
    back = X.load_npz(p1)
    assert set(back) == set(arrays) and np.array_equal(back['a'], arrays['a']) and str(back['b|c']) == 'text' and np.array_equal(back['gt'], arrays['gt'])


def test_parts_assemble_into_one_raw_file(X, tmp_path):
    d = str(tmp_path / 'parts')
    X.save_part('Tubulins_I', 'meta', {'n': np.array(3)}, d)
    X.save_part('Tubulins_I', 'isdev.hp5h', {'frames': np.array([0, 4])}, d)
    X.save_part('Bundled_Tubes_Long_Sequence', 'cal', {'a': np.array(1.41)}, d)
    assert X.part_names('Tubulins_I', d) == ['isdev.hp5h', 'meta'] and X.have_part('Tubulins_I', 'meta', d) and not X.have_part('Tubulins_II', 'meta', d)
    p = str(tmp_path / 'raw.npz')
    X.assemble(X.SETS, d, p, log=lambda s: None)
    R = X.load_raw(p)
    assert int(R['Tubulins_I']['meta|n']) == 3 and list(R['Tubulins_I']['isdev.hp5h|frames']) == [0, 4] and float(R['Bundled_Tubes_Long_Sequence']['cal|a']) == 1.41
    assert json.loads(str(R['sets'])) == ['Tubulins_I', 'Bundled_Tubes_Long_Sequence']
    q = str(tmp_path / 'raw2.npz')
    X.assemble(X.SETS, d, q, log=lambda s: None)
    assert open(p, 'rb').read() == open(q, 'rb').read()                                     # the same parts give the same bytes


def _synthetic_noise_frames(rng, n=24, ny=64, nx=64, a=1.41, b0=150.0, emitters=2):
    """Frames of a smooth background blob (10 to 200 photons, a scale of 18 px, like the autofluorescence) with Gaussian noise of variance a x + b0, and two bright blobs (emitters) per frame."""
    yy, xx = np.mgrid[0:ny, 0:nx]
    bg = 10.0 + 190.0 * np.exp(-0.5 * ((xx - nx / 2) ** 2 + (yy - ny / 2) ** 2) / 18.0 ** 2)
    fr = bg[None] + rng.normal(size=(n, ny, nx)) * np.sqrt(a * bg[None] + b0)
    for f in range(n):
        for _ in range(emitters):
            cx, cy = rng.uniform(5, nx - 5), rng.uniform(5, ny - 5)
            fr[f] += 800.0 * np.exp(-0.5 * ((xx - cx) ** 2 + (yy - cy) ** 2) / 1.0 ** 2)
    return fr


def test_noise_law_recovers_the_constant_of_a_known_noise_model(X):
    rng = np.random.default_rng(7)
    fr = _synthetic_noise_frames(rng)
    rows, b0 = X.noise_law(fr, a=1.41)
    assert 0.92 * 150.0 < b0 < 1.12 * 150.0                                                  # the robust high-frequency variance, emitters and all (the curvature of a smooth blob adds a few percent)
    clean = _synthetic_noise_frames(rng, emitters=0)
    rows_c, b0_c = X.noise_law(clean, a=1.41)
    af, bf, span = X.free_slope(rows_c)
    assert span > 100.0 and 1.2 < af < 1.65 and 0.9 * 150.0 < b0_c < 1.12 * 150.0             # the free slope of a wide-lever-arm set without emitters is the page's slope
    flat = rng.normal(size=(24, 64, 64)) * np.sqrt(60.0) + 40.0
    rows2, b02 = X.noise_law(flat, a=1.41)
    assert abs(b02 + 1.41 * 40.0 - 60.0) < 0.1 * 60.0                                        # variance 60 at level 40: b0 = 60 - a x


def test_background_surface_is_the_smoothed_temporal_median(X):
    rng = np.random.default_rng(9)
    base = np.linspace(0, 100, 32)[None, :] * np.ones((32, 1))
    fr = base[None] + rng.normal(0, 5.0, size=(31, 32, 32))
    fr[5, 10:14, 10:14] += 5000.0                                                            # an emitter in one frame does not move the median
    s = X.background_surface(fr)
    assert s.shape == (32, 32) and np.abs(s[4:-4, 4:-4] - base[4:-4, 4:-4]).max() < 3.0


def test_fit_width_and_sigma0_estimate_recover_known_widths(X):
    from scipy.special import erf
    rng = np.random.default_rng(11)
    half = 5
    ax = np.arange(-half, half + 1, dtype=float)
    U, V = np.meshgrid(ax, ax)

    def patch(sig, xc, yc, N=3000.0):
        pr = lambda c: 0.5 * (erf((ax + 0.5 - c) / (np.sqrt(2) * sig)) - erf((ax - 0.5 - c) / (np.sqrt(2) * sig)))
        return 60.0 + 2.0 * U - 1.0 * V + N * np.outer(pr(yc), pr(xc)) + rng.normal(0, 4.0, U.shape)

    for sig in (0.7, 1.0, 1.3):
        t, rms = X.fit_width(patch(sig, 0.2, -0.3), half, 0.0, 0.0)
        assert abs(t[3] - sig) < 0.05 * sig + 0.02 and abs(t[1] - 0.2) < 0.1 and abs(t[2] + 0.3) < 0.1 and rms < 8.0
    # sigma0_estimate: the selection rules (bright, isolated, interior)
    n, ny, nx, px = 3, 40, 40, 100.0
    frames = np.zeros((n, ny, nx))
    gt = np.zeros(6, dtype=C.GT_DTYPE)
    spots = [(0, 10.0, 10.0, 3000.0), (0, 30.0, 30.0, 3000.0), (1, 10.0, 10.0, 3000.0), (1, 14.0, 10.0, 3000.0), (2, 20.0, 20.0, 100.0), (2, 1.0, 20.0, 3000.0)]
    for i, (f, x, y, ph) in enumerate(spots):
        gt[i] = (i + 1, f, (x + 0.5) * px, (y + 0.5) * px, 0.0, ph)
        pr = lambda c, ax_: 0.5 * (erf((ax_ + 0.5 - c) / (np.sqrt(2) * 0.9)) - erf((ax_ - 0.5 - c) / (np.sqrt(2) * 0.9)))
        frames[f] += ph * np.outer(pr(y, np.arange(ny)), pr(x, np.arange(nx)))
    frames += rng.normal(0, 1.0, frames.shape)
    out = X.sigma0_estimate(frames, gt, 0, px, half=5, iso_px=8.0, bright=0.4)
    assert len(out) == 2                                                                      # frame 0: both isolated; frame 1: the pair 4 px apart excluded; dim and edge emitters excluded
    assert np.all(np.abs(out[:, 0] - 0.9) < 0.05) and np.all(np.abs(out[:, 1]) < 10.0)


def _init_for_synthetic(X, shape, bg_level, a, b0, sigma0):
    cal = dict(bg=np.full(shape, bg_level), a=a, b0=b0, sigma0=sigma0)
    X._init('Bundled_Tubes_Long_Sequence', 'hp5', cal, dict(kinds=('snr',)))


def test_tiled_iss_finds_each_emitter_once_across_tile_boundaries(X):
    """A frame of 40 x 96 px with emitters on both sides of two tile-core boundaries (cores of 40 px): every emitter is found once, within a pixel fraction, with no duplicate."""
    from src.core import imaging as IM
    rng = np.random.default_rng(13)
    ny, nx, sigma, bgl = 40, 96, 1.0, 30.0
    P = IM.PixelPSF(ny, nx, sigma)
    xs = [12.0, 39.2, 40.8, 60.0, 79.3, 80.7, 90.0]                                         # core boundaries at x = 39.5 and 79.5
    c = np.array([[x, 10.0 + 4.0 * k % 20 + 5.0] for k, x in enumerate(xs)])
    frame = IM.simulate_frame(P, c, np.full(len(xs), 3000.0), bgl, rng, read_sd=8.0)
    old = X.RECIPES['sparse']
    X.RECIPES['sparse'] = dict(old, tile=(40, 8))
    try:
        _init_for_synthetic(X, (ny, nx), bgl, 1.0, 64.0, sigma)
        res, n_ent, n_cap, secs, cpu = X.run_frame((0, frame))
    finally:
        X.RECIPES['sparse'] = old
    d = res[('snr', 5.0)]
    assert len(d) == len(xs)                                                                  # nothing lost at a boundary, nothing counted twice
    truth = (c + 0.5) * 100.0
    for t in truth:
        dist = np.hypot(d[:, 0] - t[0], d[:, 1] - t[1])
        assert dist.min() < 40.0 and np.sum(dist < 100.0) == 1                               # each emitter within 0.4 px, exactly one detection near it
    assert n_cap == 0 and n_ent >= len(xs)


# ------------------------------------------------------------------ the real data (skipped when absent)
def _have(name):
    return S.data_available(name)


needs_data = pytest.mark.skipif(not all(_have(s) for s in S.SETS), reason='the 2013 data are not downloaded (data/smlm/fetch_2013.py)')
_PICASSO_PY = os.environ.get('SMLM_PICASSO_PYTHON')
needs_picasso = pytest.mark.skipif(not _PICASSO_PY and not C._picasso_importable(), reason='Picasso is not available (set SMLM_PICASSO_PYTHON)')


@needs_data
@pytest.mark.parametrize('name', S.SETS)
def test_real_counts_ranges_and_one_frame_activations(name):
    sp = S.SPEC[name]
    gt = S.load_gt(S.data_path(name, 'gt'))
    assert S.n_frames_in(S.data_path(name, 'stack')) == sp['n_frames']
    assert len(gt) == sp['n_gt']
    assert gt['frame'].min() == 0 and gt['frame'].max() == sp['n_frames'] - 1             # file frames 1..N become stack indices 0..N-1
    fld = S.field_nm(name)
    assert gt['x_nm'].min() > 0 and gt['x_nm'].max() < fld and gt['y_nm'].min() > 0 and gt['y_nm'].max() < fld
    assert 10.0 < gt['z_nm'].min() and gt['z_nm'].max() < 300.0 and gt['photons'].min() > 100.0
    rl = C.run_lengths(gt)
    assert rl['runs'] == len(gt) and rl['mean'] == 1.0                                       # every activation is one frame long
    per = np.bincount(gt['frame'], minlength=sp['n_frames'])
    assert per.mean() == pytest.approx(sp['n_gt'] / sp['n_frames'])
    for split in ('dev', 'test'):
        a, b = S.PROTOCOL[name][split]
        assert per[a:b].mean() > 0.9 * per.mean() and per[a:b].mean() < 1.1 * per.mean()                # the ranges are typical stretches of the sequence (empty frames occur in the long sequence)


@needs_data
@pytest.mark.parametrize('name', S.SETS)
def test_real_per_frame_files_equal_the_activation_file(name):
    a, b = S.load_gt(S.data_path(name, 'gt')), S.gt_frame_files(S.data_path(name, 'gt'))
    assert len(a) == len(b)
    o1, o2 = np.lexsort((a['y_nm'], a['x_nm'], a['frame'])), np.lexsort((b['y_nm'], b['x_nm'], b['frame']))
    assert all(np.array_equal(a[k][o1], b[k][o2]) for k in ('frame', 'x_nm', 'y_nm', 'z_nm', 'photons'))


@needs_data
def test_real_tubulins_share_one_ground_truth_and_the_bundles_one_structure():
    t1, t2 = S.load_gt(S.data_path('Tubulins_I', 'gt')), S.load_gt(S.data_path('Tubulins_II', 'gt'))
    assert all(np.array_equal(t1[k], t2[k]) for k in t1.dtype.names)
    hd, ls = S.load_gt(S.data_path('Bundled_Tubes_High_Density', 'gt')), S.load_gt(S.data_path('Bundled_Tubes_Long_Sequence', 'gt'))
    assert np.array_equal(np.sort(hd['x_nm']), np.sort(ls['x_nm'])) and np.array_equal(np.sort(hd['y_nm']), np.sort(ls['y_nm']))        # the same 81,049 sites, other frames and photons


@needs_data
@pytest.mark.parametrize('name', S.SETS)
def test_real_frames_are_16_bit_camera_counts_on_the_pages_baseline(name):
    a, b = S.PROTOCOL[name]['dev']
    fr = S.load_frames(S.data_path(name, 'stack'), a, a + 5)
    n = S.SPEC[name]['size']
    assert fr.shape == (5, n, n) and fr.min() >= 0 and fr.max() <= 16383                     # 14-bit saturation of the page
    assert np.all(fr == np.round(fr))
    assert 120.0 < np.median(fr) < 180.0                                                     # baseline 100 DN plus the autofluorescence
    assert S.to_photons(fr).min() >= -100.0                                                 # DN is clipped at 0 (4.7e-5 of the pixels of Tubulins_II, none elsewhere)


@needs_data
@pytest.mark.parametrize('name', S.SETS)
def test_real_frame_base_by_excess_intensity(name):
    a, b = S.PROTOCOL[name]['dev']
    gt = S.load_gt(S.data_path(name, 'gt'))
    fr = S.load_frames(S.data_path(name, 'stack'), a, a + 10)
    ex = S.excess_at_ground_truth(fr, gt, a, pixel_nm=S.SPEC[name]['pixel_nm'], min_photons=float(np.median(gt['photons'])))
    wrong = max(ex[o][0] for o in (-2, -1, 1, 2))
    assert ex[0][0] > (1.3 if name.endswith('High_Density') else 2.0) * wrong                # offset 0 stands out (the dense set's static structure leaves excess at any offset)


@needs_data
@needs_picasso
@pytest.mark.parametrize('name', S.SETS)
def test_real_picasso_matches_the_ground_truth_within_tolerance_on_five_frames(name):
    import exp_smlm2013 as X
    a, b = S.PROTOCOL[name]['dev']
    gt = S.load_gt(S.data_path(name, 'gt'))
    fr = S.load_frames(S.data_path(name, 'stack'), a, a + 5)
    dets = C.run_picasso(fr, dict(X.camera(name), min_net_gradient=1000.0), python=_PICASSO_PY)
    cv = S.check_conventions(name, dets, gt, np.arange(a, a + 5))
    px = S.SPEC[name]['pixel_nm']
    tp = cv['corner']['tp']
    assert tp >= (60 if name.endswith('High_Density') else 70 if name.startswith('Tub') else 15)
    assert abs(cv['corner']['bias_x_nm']) < 10.0 and abs(cv['corner']['bias_y_nm']) < 10.0           # the corner convention without any shift
    if name.endswith('High_Density'):
        # one emitter per 18 pixels: any detection finds SOME ground-truth emitter within 250 nm after an exchange of the axes or a frame shift, so the counts do not discriminate; the mean error does
        assert cv['centre']['bias_x_nm'] < -0.25 * px and cv['centre']['bias_y_nm'] < -0.25 * px         # the pairing finds nearer neighbours, so less than -p/2, but clearly negative
        assert cv['frames'][0] == max(cv['frames'].values())
    else:
        assert cv['centre']['bias_x_nm'] < -0.4 * px and cv['centre']['bias_y_nm'] < -0.4 * px            # the other convention is about -p/2 off
        assert cv['swapped']['tp'] < 0.4 * tp and cv['flip_x']['tp'] < 0.4 * tp and cv['flip_y']['tp'] < 0.4 * tp
        assert all(cv['frames'][o] < 0.3 * tp for o in (-2, -1, 1, 2))                                   # and the frame base is 0


# ------------------------------------------------------------------ the downloader (no network)
def _fetch_module():
    import importlib.util
    path = os.path.join(os.path.dirname(__file__), '..', 'data', 'smlm', 'fetch_2013.py')
    spec = importlib.util.spec_from_file_location('fetch_2013', path)
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, os.path.dirname(path))
    spec.loader.exec_module(mod)
    return mod


def test_fetch_2013_urls_sizes_and_recorded_digests_are_consistent():
    F13 = _fetch_module()
    assert F13.SETS == S.SETS and set(F13.EXPECTED_BYTES) == set(S.SETS)
    for name in S.SETS:
        files = F13.files_for(name, F13.KINDS)
        assert set(files) == set(F13.KINDS) == set(F13.EXPECTED_BYTES[name])
        for kind, (rel, url) in files.items():
            assert rel == '%s/%s.zip' % (name, kind)
            assert url == 'https://srm.epfl.ch/srm/Data/challenge-2D-simulation/%s/Data/%s.zip' % (name, kind)       # the sub-page's link ../../../Data/... resolved
    assert F13.DEFAULT_KINDS == ('sequence', 'fluorophores')                                                       # the ground truth is in fluorophores.zip, not oracle.zip
    book = json.load(open(F13.CHECKSUMS))
    for name in S.SETS:
        for kind in F13.DEFAULT_KINDS:
            rec = book['%s/%s.zip' % (name, kind)]
            assert len(rec['sha256']) == 64 and rec['bytes'] == F13.EXPECTED_BYTES[name][kind] and rec['url'] == F13.files_for(name, F13.KINDS)[kind][1]
    for name in S.SETS:                                                                                            # the page facts name the frame size the loader uses
        assert '%d x %d px' % (S.SPEC[name]['size'], S.SPEC[name]['size']) in F13.PAGE_FACTS[name]['frame']
        assert ('pixel %g nm' % S.SPEC[name]['pixel_nm']) in F13.PAGE_FACTS[name]['frame']


def test_fetch_2013_reads_member_facts_from_synthetic_archives(tmp_path):
    F13 = _fetch_module()
    stack = np.arange(3 * 4 * 5).reshape(3, 4, 5) + 100
    p = _zip_of_frames(tmp_path / 'sequence.zip', stack, 'sequence/')
    sf = F13.sequence_facts(p)
    assert sf[0] == 3 and sf[1] == 1 and sf[2] == 3 and sf[3] is True and sf[4:6] == (4, 5)                         # 3 frames numbered 1..3, 4 rows x 5 columns
    q = tmp_path / 'fluorophores.zip'
    with zipfile.ZipFile(q, 'w') as z:
        z.writestr('fluorophores/activation.csv', GT_TEXT)
        z.writestr('fluorophores/frames/00001.csv', 'x\n')
        z.writestr('fluorophores/frames/00003.csv', 'x\n')
    gf = F13.gt_facts(str(q))
    assert gf[0] == 'fluorophores/activation.csv' and gf[2] == 3 and gf[3] == 1 and gf[4] == 3 and gf[5] == 2
    assert F13.tiff_members([('a/00010.tif', 5), ('a/readme.txt', 5), ('a/00002.tif', 5), ('a/00003.tif', 0)]) == [(2, 'a/00002.tif'), (10, 'a/00010.tif')]
